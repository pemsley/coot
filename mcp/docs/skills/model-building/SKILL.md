---
name: coot-model-building
description: "Best Practices for Model-Building Tools and Refinement"
---

# Key Lessons from Chain A Refinement Session

## Workflow Checklist - Follow This for Model-Building and Refinement

Best practices for fixing any issue:
1. ☐ Center on the interesting residue: `coot.set_go_to_atom_chain_residue_atom_name(chain, resno, "CA")`
     or interesting position `coot.set_rotation_centre(x,y,z)`
2. ☐ Check current metrics (Rama/correlation/overlaps)
3. ☐ Make checkpoint if trying something experimental
4. ☐ Apply fix
5. ☐ Re-check metrics to confirm improvement
6. ☐ If worse, restore checkpoint

## Jiggle Fitting a Whole Chain

### When to use jiggle fitting

Use `fit_chain_to_map_by_random_jiggle_and_blur()` when a whole chain needs to be
placed or re-placed in the map — before any residue-level refinement.
Note, however, that the fit is only local. If the whole chain is (say) 10 or more Angstroms
away from the real position, then `fit_chain_to_map_by_random_jiggle_and_blur()` will not work.
The per-residue density correlation graph is the key diagnostic: if the mean correlation across
the whole chain is near zero with high variance (many residues below 0.4, many negative),
the chain is in the wrong place and jiggle fitting is the right first step.

**Do not attempt residue-level fixes (rotamers, pepflips, refinement) until the chain
is correctly placed.** Those tools cannot recover a globally misplaced chain.

### Diagnosing misplacement from the correlation graph

Before jiggling, always get a fresh per-residue correlation to use as a true baseline:

```python
stats = coot.map_to_model_correlation_stats_per_residue_range_py(
    imol, chain_id, imol_map, 1, 0)
corrs = [entry[1][1] for entry in stats[0]]
mean_corr = sum(corrs) / len(corrs)
n_below = sum(1 for c in corrs if c < 0.4)
print(f"Mean correlation: {mean_corr:.3f}, residues below 0.4: {n_below}/{len(corrs)}")
```

**Interpretation:**
- Mean > 0.7, few residues below 0.4 → chain is well placed, proceed to residue-level work
- Mean 0.4–0.7, moderate number below 0.4 → chain placed but needs refinement
- Mean near zero, high variance, many negative → **chain is in the wrong place**, jiggle fit first

**Why negative correlations occur in cryo-EM maps:** cryo-EM maps are Fourier-based
(no F000 term), so the mean density is zero. Negative density exists between atomic
features as a mathematical consequence of the reconstruction — not a physical reality.
A misplaced chain within the molecular envelope samples this real structured signal
(both positive peaks and negative troughs) from the wrong part of the map, producing
scattered correlations around zero. This is distinct from sitting in solvent, which
gives low-magnitude, low-variance correlations.

### Running the jiggle fit

Set the refinement map first, then use the blur variant for cryo-EM:

```python
coot.set_imol_refinement_map(imol_map)

result = coot.fit_chain_to_map_by_random_jiggle_and_blur(
    imol,          # model molecule
    chain_id,      # e.g. "N"
    n_trials,      # 200-500; more trials = better chance of finding global optimum
    jiggle_scale,  # 1.0-2.0; scale of random perturbations
    blur_factor    # 100-200; Fourier filtering makes the map smoother, helping
                   # the fit escape local minima and find the correct placement
)
print(f"Jiggle fit score: {result:.3f}")
```

**Return values:** the function returns a score (higher = better fit) or a large
negative number (e.g. −100, −999) when it cannot improve on the current placement.
These large negative returns are not errors — they mean the current position is already
a local optimum, or the function had difficulty scoring. Do not treat them as failures
without checking the actual correlation afterwards.

**Blur factor guidance:** the blur smooths the map by a Fourier filter, effectively
low-pass filtering it. A higher blur factor produces a smoother map with broader
features, which helps the chain find the correct region before fine-grained fitting.
For cryo-EM maps, blur factors of 100–200 are typical starting points.

### Confirming the result

Always re-run the per-residue correlation after jiggling and compare to the baseline:

```python
stats_after = coot.map_to_model_correlation_stats_per_residue_range_py(
    imol, chain_id, imol_map, 1, 0)
corrs_after = [entry[1][1] for entry in stats_after[0]]
mean_after = sum(corrs_after) / len(corrs_after)
n_below_after = sum(1 for c in corrs_after if c < 0.4)
print(f"AFTER  — mean: {mean_after:.3f}, below 0.4: {n_below_after}")
```

A successful jiggle fit for a correctly placed chain should bring the mean correlation
above 0.6, with most residues above 0.7. Remaining poor-fitting residues in loop
regions are expected and addressed by subsequent residue-level refinement.

### Visualising the before/after result

Use the `coot-inline-graphs` skill to render a before/after bar chart with secondary
structure overlays. Fetch secondary structure from PDBe
(`/api/pdb/entry/secondary_structure/{pdb_id}`) rather than relying solely on Coot's
header secondary structure detection, which may miss strands — particularly single-
residue assignments and short strands in nanobody/VHH-type folds.

### Complete workflow summary

```python
# 1. Get fresh baseline correlation (after any manual moves, before jiggle)
stats_before = coot.map_to_model_correlation_stats_per_residue_range_py(
    imol, chain_id, imol_map, 1, 0)
corrs_before = [e[1][1] for e in stats_before[0]]
print(f"BEFORE — mean: {sum(corrs_before)/len(corrs_before):.3f}")

# 2. Set refinement map and jiggle fit
coot.set_imol_refinement_map(imol_map)
result = coot.fit_chain_to_map_by_random_jiggle_and_blur(imol, chain_id, 500, 1.0, 200.0)
print(f"Jiggle score: {result:.3f}")  # large negative = ignore, check correlation anyway

# 3. Get after correlation and compare
stats_after = coot.map_to_model_correlation_stats_per_residue_range_py(
    imol, chain_id, imol_map, 1, 0)
corrs_after = [e[1][1] for e in stats_after[0]]
print(f"AFTER  — mean: {sum(corrs_after)/len(corrs_after):.3f}")

# 4. Only proceed to residue-level work once mean > 0.6
```

## Refinement Best Practices

### 1. Make Checkpoints Before Model Changes
**CRITICAL: Always create a checkpoint before making significant model changes.**

Use `make_backup_checkpoint()` before any operation that might need to be reverted:
- Adding/deleting residues
- Adding ligands or waters
- Major refinement operations
- Any experimental model building

```python
# Create a named checkpoint before risky operation
checkpoint_idx = coot.make_backup_checkpoint(0, "before adding OXT")

# Try the operation
coot.add_OXT_to_residue(0, "A", 93, "")
result = coot.refine_residues_py(0, [["A", 93, ""]])

# Check if it worked - if not, restore
if result_is_bad:
    coot.restore_to_backup_checkpoint(0, checkpoint_idx)
```

**Why checkpoints are better than undo:**
- `apply_undo()` only steps back one operation at a time
- Checkpoints let you jump back to a specific point
- Named checkpoints are self-documenting
- Multiple checkpoints allow comparing different approaches

```python
# Compare two different approaches
checkpoint_before = coot.make_backup_checkpoint(0, "original state")

# Try approach 1
coot.auto_fit_best_rotamer(0, "A", 42, "", "", 1, 1, 0.01)
results = coot.refine_residues_py(0, [["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])
score_approach1 = check_correlation(0, "A", 42)
checkpoint_approach1 = coot.make_backup_checkpoint(0, "after approach 1")

# Restore and try approach 2
coot.restore_to_backup_checkpoint(0, checkpoint_before)
coot.pepflip(0, "A", 42, "", "")
results = coot.refine_residues_py(0, [["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])
score_approach2 = check_correlation(0, "A", 42)

# Keep the better result
if score_approach1 > score_approach2:
    coot.restore_to_backup_checkpoint(0, checkpoint_approach1)
```

## MANDATORY: Complete Validation Workflow

**CRITICAL: You MUST check ALL validation metrics before AND after EVERY fix.**

Fixing only one problem (e.g., Ramachandran) while leaving others (rotamer, density fit) is a FAILED fix. A residue is only "fixed" when ALL metrics are acceptable.

### Before fixing ANY residue:
1. **ALWAYS** center on it: `coot.set_go_to_atom_chain_residue_atom_name(chain, resno, "CA")`
2. **ALWAYS** check ALL of these:
   - Ramachandran probability (from `all_molecule_ramachandran_score_py`)
   - Rotamer score (from `rotamer_graphs_py`)
   - Density correlation - both all-atom and side-chain (from `map_to_model_correlation_stats_per_residue_range_py`)
   - Atom overlaps involving this residue (from `molecule_atom_overlaps_py`)

### After fixing ANY residue:
3. **ALWAYS** re-check ALL the same metrics, Sometime residues/issues are just
     not fixable (that's what makes refinement and model-building tricky).
4. **ONLY MOVE** on to the next residue/issue unless you have tried to make all of
     these are true:
   - Ramachandran probability > 0.02 (preferably > 0.1)
   - Rotamer score > 1.0% (preferably > 5%)
   - Density correlation > 0.7 (all-atom and side-chain, preferably > 0.8)
   - No severe clashes (< 2.0 Å cubed overlap volume)

### If ANY metric is still bad after your first fix:
5. **MUST** try additional fixes:
   - Bad rotamer → `auto_fit_best_rotamer()`, and try experiment with following that up
     with refine_residues_py() for that residue and its upstream and downstream
     neighbours (if any).
   - Poor density fit → try alternative rotamers, check for missing atoms
   - Persistent clashes → refine with the addition of spatial neighbors using
     `residues_near_residue()`
6. **NEVER** declare a residue "fixed" based on only one metric improving
7. **ALWAYS** re-validate after each additional fix

### Acceptable Reasons to Stop (without perfect metrics):
- You've tried multiple approaches and documented them
- The best achievable metrics are still recorded
- You've created checkpoints to compare approaches
- You explain why the residue remains problematic (e.g., poor density, crystal contact)

## Example of CORRECT Workflow

```python
# 1. ALWAYS center on problem residue first
coot.set_go_to_atom_chain_residue_atom_name("A", 41, "CA")

# 2. Get ALL metrics BEFORE
rama_data = [r for r in coot.all_molecule_ramachandran_score_py(0)[5:]
             if r[1] == ['A', 41, '']][0]
rama_prob_before = rama_data[2]

rotamer_data = [r for r in coot.rotamer_graphs_py(0)
                if r[0] == 'A' and r[1] == 41][0]
rotamer_score_before = rotamer_data[3]

corr_data = [s for s in coot.map_to_model_correlation_stats_per_residue_range_py(0, "A", 1, 1, 0)[0]
             if s[0][1] == 41][0]
correlation_before = corr_data[1][1]

overlaps_before = [o for o in coot.molecule_atom_overlaps_py(0, 30)
                   if (o['atom-1-spec'][1:3] == ['A', 41] or
                       o['atom-2-spec'][1:3] == ['A', 41])]

print(f"BEFORE: Rama={rama_prob_before:.4f}, Rotamer={rotamer_score_before:.2f}%, Corr={correlation_before:.3f}, Clashes={len(overlaps_before)}")

# 3. Apply first fix (e.g., pepflip for backbone)
coot.pepflip(0, "A", 41, "", "")
coot.refine_residues_py(0, [["A", 40, ""], ["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])

# 4. Check ALL metrics AFTER first fix
rama_prob_after = [r for r in coot.all_molecule_ramachandran_score_py(0)[5:] 
                   if r[1] == ['A', 41, '']][0][2]
rotamer_score_after = [r for r in coot.rotamer_graphs_py(0) 
                       if r[0] == 'A' and r[1] == 41][0][3]
correlation_after = [s for s in coot.map_to_model_correlation_stats_per_residue_range_py(0, "A", 1, 1, 0)[0] 
                     if s[0][1] == 41][0][1][1]

print(f"AFTER:  Rama={rama_prob_after:.4f}, Rotamer={rotamer_score_after:.2f}%, Corr={correlation_after:.3f}")

# 5. If rotamer or correlation still bad, DON'T STOP - fix them!
if rotamer_score_after < 1.0:
    print("Rotamer still bad - trying auto_fit_best_rotamer")
    coot.auto_fit_best_rotamer(0, "A", 41, "", "", 1, 1, 0.01)
    coot.refine_residues_py(0, [["A", 40, ""], ["A", 41, ""], ["A", 42, ""]])
    
    # 6. ALWAYS re-check after additional fixes
    rotamer_score_final = [r for r in coot.rotamer_graphs_py(0) 
                           if r[0] == 'A' and r[1] == 41][0][3]
    correlation_final = [s for s in coot.map_to_model_correlation_stats_per_residue_range_py(0, "A", 1, 1, 0)[0] 
                         if s[0][1] == 41][0][1][1]
    print(f"FINAL:  Rotamer={rotamer_score_final:.2f}%, Corr={correlation_final:.3f}")

# 7. Only NOW can you move to the next residue
```

## Example of WRONG Workflow (DO NOT DO THIS)

```python
# ❌ WRONG: Checking only Ramachandran
coot.pepflip(0, "A", 41, "", "")
coot.refine_residues_py(0, [["A", 40, ""], ["A", 41, ""], ["A", 42, ""]])
rama_after = coot.all_molecule_ramachandran_score_py(0)[5][39][2]
print(f"Ramachandran improved to {rama_after}")
# MOVES ON without checking rotamer or density fit - WRONG!

# ❌ WRONG: Not centering on residue
# Goes straight to fix without set_go_to_atom_chain_residue_atom_name()

# ❌ WRONG: Not checking metrics before the fix
# How do you know if it improved if you don't know what it was before?

# ❌ WRONG: Declaring success with bad rotamer
rama = 0.30  # Good!
rotamer = 0.0001  # TERRIBLE!
correlation = 0.59  # POOR!
print("Residue fixed!")  # NO IT ISN'T!
```

## Why This Matters

A residue with:
- ✅ Good Ramachandran (0.30)
- ❌ Terrible rotamer (0.01%)
- ❌ Poor density fit (0.59)

is NOT fixed. The side chain is clearly wrong. The backbone geometry might be OK, but the model is still incorrect.

**ALL metrics must be acceptable before moving on.**



### 2. Extend Selection Around Problem Residues
**Don't refine problem residues in isolation** - include neighboring residues for context.

- ❌ **Bad:** `refine_residues_py(0, [["A", 41, ""]])`  - Often fails to correct the model
- ✅ **Good:** `refine_residues_py(0, [["A", 40, ""], ["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])`

**Recommended approach:**
- For single problem residue: include ±1 or ±2 neighbors
- For consecutive problem residues: include ±1 neighbor on each end
- Larger regions (±3-4 residues) can sometimes help severe issues

**Example from session:**
- Residues 41-42 had severe Ramachandran outliers
- Refining just 41-42 failed
- Refining 40-43 succeeded: Residue 41 improved from p=0.00004 to p=0.308

**Neighboring Residues:**
- You can use coot.residues_near_residue() to find residues that are close in space
  but distant in sequence, so that they can be added to the residue selection for
  refinement.

### 3. Include Spatial Neighbours, Not Just Sequence Neighbours

**Critical insight:** Residues that are close in 3D space affect each other during refinement, even if they're far apart in sequence.

Coot's refinement includes spatially neighbouring atoms in the non-bonded contact interactions, but only the selected residues can move during minimization. If a nearby (but unselected) residue is in the wrong position, it will "push" your selected residues away via non-bonded contact penalties - potentially pushing them out of correct density to avoid the clash with the incorrectly-placed neighbour.

**Diagnostic workflow:**

1. **Check for clashes** after refinement:
```python
overlaps = coot.molecule_atom_overlaps_py(0, 50)
for o in overlaps:
    spec1, spec2 = o['atom-1-spec'], o['atom-2-spec']
    vol = o['overlap-volume']
    if vol > 0.5:  # Significant clash
        print(f"{spec1[1]}/{spec1[2]} {spec1[4]} - {spec2[1]}/{spec2[2]} {spec2[4]}: {vol:.2f}")
```

2. **If a problem residue clashes with a distant residue**, fix the distant residue first:
```python
# Example: A/2 has poor correlation (0.13) and clashes with A/89
# First fix A/89:
coot.auto_fit_best_rotamer(0, "A", 89, "", "", 1, 1, 0.01)
results = coot.refine_residues_py(0, [["A", 88, ""], ["A", 89, ""], ["A", 90, ""]])

# Then re-refine A/2 INCLUDING A/89 as a spatial neighbour:
results = coot.refine_residues_py(0, [["A", 1, ""], ["A", 2, ""], ["A", 3, ""], ["A", 89, ""]])
# A/2 correlation improved: 0.13 → 0.81
```

**Why this matters:**
- Coot's refinement "feels" spatial neighbours via non-bonded contact terms
- But only selected residues can move during minimization
- A/89 was in a wrong position (correlation 0.050) and pushing A/2 away
- A/2 moved out of its correct density to reduce the non-bonded penalty with A/89
- Fixing A/89 first put it in the right place, so it no longer pushed A/2 incorrectly

**Real example:**
```
Before: A/2 correlation = 0.131, A/89 correlation = 0.050
        A/2 CA ↔ A/89 CZ clash: 1.06 Ų

After fixing A/89 first, then refining together:
        A/2 correlation = 0.805, A/89 correlation = 0.928
        No clash
```

### 4. Iterative Refinement Strategy
Sometimes multiple rounds of refinement with different selections help:

1. **First pass:** Refine larger region to establish general geometry
2. **Second pass:** Refine smaller region to fine-tune specific problem
3. **Check validation** after each step
4. **Restore checkpoint** if results get worse

**Example workflow:**
```python
# Create checkpoint first!
checkpoint = coot.make_backup_checkpoint(0, "before iterative refinement")

# First: larger region
results_1 = coot.refine_residues_py(0, [["A", i, ""] for i in range(40, 44)])
check_validation()  # Did it help?

# Second: targeted refinement
results_2 = coot.refine_residues_py(0, [["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])
check_validation()  # Better or worse?

# If worse:
coot.restore_to_backup_checkpoint(0, checkpoint)
```

### 5. Measure Before and After
**Always validate changes objectively** using:
- Ramachandran probabilities
- Density correlation (all-atom and side-chain)
- Geometry statistics
```python
def check_residue_validation(imol, chain_id, resno):
    """Check both Ramachandran and density correlation"""
    # Get Ramachandran
    rama_data = coot.all_molecule_ramachandran_score_py(imol)
    residue_data = rama_data[5]
    rama_score = None
    for r in residue_data:
        if r[1][0] == chain_id and r[1][1] == resno:
            rama_score = r[2]
            break

    # Get density correlation
    corr_data = coot.map_to_model_correlation_stats_per_residue_range_py(
        imol, chain_id, 1, 1, 1
    )
    all_atom_corr = None
    sidechain_corr = None

    for r in corr_data[0]:
        if r[0][1] == resno:
            all_atom_corr = r[1][1]
            break

    for r in corr_data[1]:
        if r[0][1] == resno:
            sidechain_corr = r[1][1]
            break

    return {
        'residue': resno,
        'rama_prob': rama_score,
        'all_atom_corr': all_atom_corr,
        'sidechain_corr': sidechain_corr
    }

# Usage with checkpoint
checkpoint = coot.make_backup_checkpoint(0, "before refinement test")
before = check_residue_validation(0, "A", 41)

coot.refine_residues_py(0, [["A", 40, ""], ["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])
coot.accept_moving_atoms_py()
after = check_residue_validation(0, "A", 41)

# Compare and decide
if after['all_atom_corr'] > before['all_atom_corr']:
    # Keep it!
    pass
else:
    # Revert to checkpoint
    coot.restore_to_backup_checkpoint(0, checkpoint)
```

### 6. Use Checkpoints and Undo Liberally
**Don't be afraid to revert changes:**

- `make_backup_checkpoint()` / `restore_to_backup_checkpoint()` - for jumping back to a specific state
- `apply_undo()` - for stepping back one operation at a time

Use checkpoints when:
- Starting a new model-building task
- About to try something experimental
- Before a series of related operations

Use undo when:
- The last single operation made things worse
- Quick single-step revert needed

### 7. Auto-fit Rotamer for Side-chain Issues
**For poor side-chain density correlation**, try `auto_fit_best_rotamer()` first:
```python
# Create checkpoint first
checkpoint = coot.make_backup_checkpoint(0, "before rotamer fitting")

# Check if it's a side-chain problem
validation = check_residue_validation(0, "A", 89)
if validation['sidechain_corr'] < 0.5:
    # Try auto-fit rotamer
    score = coot.auto_fit_best_rotamer(0, "A", 89, "", "", 1, 1, 0.01)
    
    if score > 0:  # Positive score is good
        # Check improvement
        after = check_residue_validation(0, "A", 89)
        if after['sidechain_corr'] > validation['sidechain_corr']:
            # Success! (e.g., 0.034 → 0.900)
            pass
        else:
            coot.restore_to_backup_checkpoint(0, checkpoint)
    else:
        # Negative score means failure
        coot.restore_to_backup_checkpoint(0, checkpoint)
```

### 8. Set Refinement to Synchronous Mode
**Always call this at the start** to make refinement complete immediately:
```python
coot.set_refinement_immediate_replacement(1)
```

Without this, refinement may be asynchronous and difficult to control programmatically.

### 9. Navigate to Residue Before Working
**Bring residue to screen center** so you can watch the refinement:
```python
coot.set_go_to_atom_molecule(0)
coot.set_go_to_atom_chain_residue_atom_name("A", 41, "CA")
```

This helps with:
- Visual inspection of the problem
- Seeing the refinement in real-time
- Verifying the result makes geometric sense

### 10. Flipping peptides

If the Ramachandran Plot is poor, try using `coot.pepflip(imol, chain_id, res_no, ins_code, alt_conf)` followed by a refinement of the residues in the extended region.

### 11. Flipping side-chains terminal Chi-angle

If the Rotamer score is poor, try using `coot.do_180_degree_side_chain_flip()` to improve the Rotamer score. It is occasionally useful.

### 12. Delete Waters Displaced by Newly Built Residues

**CRITICAL: After adding any new residues, always check for and delete water molecules that are now too close to the new protein atoms.**

Waters are placed by earlier refinement cycles that had no knowledge of atoms you are about to build. Once new residues are placed, any water within ~3.5 Å is almost certainly spurious — it was filling density that now belongs to the protein model.

**Always do this immediately after building new residues, before refinement.**

```python
# After building new residues (e.g. A21-A25), check each for clashing waters
newly_built = [21, 22, 23, 24, 25]
waters_to_delete = set()

for resno in newly_built:
    neighbours = coot.residues_near_residue_py(0, ["A", resno, ""], 3.5)
    for n in neighbours:
        chain, nresno, inscode = n[0], n[1], n[2]
        name = coot.residue_name_py(0, chain, nresno, inscode)
        if name == "HOH":
            waters_to_delete.add((chain, nresno, inscode))
            print(f"Water {chain}{nresno} within 3.5A of new A{resno} - will delete")

for chain, resno, inscode in waters_to_delete:
    coot.delete_residue(0, chain, resno, inscode)
    print(f"Deleted HOH {chain}{resno}")
```

**Why 3.5 Å?** This is a generous hydrogen-bond distance. Any water oxygen closer than this to a protein heavy atom is physically incompatible with the new model.

**Don't just check for hard clashes** — even waters that don't show up in `molecule_atom_overlaps_py()` should be deleted if they are within 3.5 Å, because the overlap check uses van der Waals radii and may miss close-but-not-overlapping waters that are nonetheless chemically unreasonable.


### 13. Building N-terminal Extensions with Enforced Secondary Structure

When building backwards (toward the N-terminus) into a known secondary structure element, use `add_terminal_residue_using_phi_psi()` with the target φ/ψ angles and enforce the conformation further by refining with `set_secondary_structure_restraints_type()`.

**α-helix:** φ = −57°, ψ = −47°
**β-strand:** φ = −120°, ψ = +120°

**Complete workflow for N-terminal helix extension:**

```python
# Step 1: Checkpoint
cp = coot.make_backup_checkpoint(0, "before N-terminal extension")

# Step 2: Build residues backwards one at a time with helix phi/psi
phi, psi = -57.0, -47.0
current_nterm = 26  # existing N-terminal residue of fragment
for i in range(4):  # build 4 residues: 25, 24, 23, 22
    r = coot.add_terminal_residue_using_phi_psi(0, "A", current_nterm, "auto", phi, psi)
    current_nterm -= 1

# Step 3: IMMEDIATELY check for and delete displaced waters
newly_built = list(range(current_nterm + 1, 26))
waters_to_delete = set()
for resno in newly_built:
    for n in coot.residues_near_residue_py(0, ["A", resno, ""], 3.5):
        if coot.residue_name_py(0, n[0], n[1], n[2]) == "HOH":
            waters_to_delete.add((n[0], n[1], n[2]))
for chain, resno, inscode in waters_to_delete:
    coot.delete_residue(0, chain, resno, inscode)
    print(f"Deleted displaced water {chain}{resno}")

# Step 4: Mutate from ALA placeholders to correct sequence
coot.mutate_residue_range(0, "A", 22, 25, "NLLND")  # sequence from UniProt etc.

# Step 5: Autofit rotamers for mutated residues
for resno in newly_built:
    coot.auto_fit_best_rotamer(0, "A", resno, "", "", 1, 1, 0.01)

# Step 6: Refine with helix restraints - always turn off afterwards
coot.set_secondary_structure_restraints_type(1)
coot.refine_residues_py(0, [["A", r, ""] for r in range(22, 28)])
coot.accept_moving_atoms_py()
coot.set_secondary_structure_restraints_type(0)  # ALWAYS turn off when done
```

**Key rules for secondary structure restraints:**
- Always call `set_secondary_structure_restraints_type(0)` after refinement to clear them
- Extend the refinement selection 1-2 residues beyond the newly built region into the existing chain, to provide proper geometric context at the junction
- Helix restraints (type 1), strand restraints (type 2), no restraints (type 0)

```python
coot.set_secondary_structure_restraints_type(1)  # alpha helix
coot.set_secondary_structure_restraints_type(2)  # beta strand
coot.set_secondary_structure_restraints_type(0)  # off - ALWAYS restore to this
```


## NCS-Related Chains: Always Use Copy, Never Rebuild

**CRITICAL RULE: When the user asks you to fix a region in one chain using a corrected region
in another chain, ALWAYS use `copy_residue_range_from_ncs_master_to_others()` first.
Never attempt to rebuild the target chain manually.**

### Trigger Phrases That Must Invoke This Approach

Any of these mean: use the NCS copy tool immediately:
- *"Use the fixed B chain to fix chain A"*
- *"Apply the same fix to the other chain"*
- *"The equivalent region in chain X needs fixing"*
- *"Propagate this fix to the NCS-related chain"*
- *"Now do the same for chain A/B/C..."*
- Any mention of NCS, related chains, or symmetry-equivalent regions after a fix

### NCS Copy Workflow

```python
# 1. Check NCS is detected — do this first
ncs_chains = coot.ncs_chain_ids_py(imol)
print(ncs_chains)  # e.g. [['A', 'B']] — master chain listed first

# 2. Set the fixed chain as NCS master
coot.ncs_control_change_ncs_master_to_chain_id(imol, "B")  # B is the fixed chain

# 3. Copy the fixed region to all NCS-related chains using the NCS operator
coot.copy_residue_range_from_ncs_master_to_others(imol, "B", start_resno, end_resno)

# 4. Refine the target chain in its own density to let it settle
specs = [["A", r, ""] for r in range(start_resno - 1, end_resno + 2)
         if coot.residue_name_py(imol, "A", r, "") not in [None, "False", ""]]
coot.refine_residues_py(imol, specs)
```

**Why this is always better than manual rebuilding:**
- The NCS operator correctly handles the rotation and translation between chains —
  NCS-related chains are in different parts of the unit cell, so you cannot simply
  copy coordinates without applying the NCS transformation
- It propagates correct backbone geometry that is already known to be good
- It is fast — one function call instead of dozens of pepflips and refinement cycles
- Refinement against the target chain's own density then handles the small differences
  between NCS-related chains

**Only fall back to manual rebuilding if:**
- `ncs_chain_ids_py()` returns `False` (no NCS detected)
- The chains have genuinely different sequences in the region of interest
- The NCS copy produces clearly wrong geometry (check correlations immediately after)

---

## Per-Atom Density Probing: The Most Powerful Diagnostic

**Use `density_at_point()` per backbone atom to diagnose backbone problems before
attempting any fix. This is more informative than per-residue correlation alone.**

```python
def probe_backbone(imol, imol_map, imol_diff, chain, resnos):
    sigma      = coot.map_sigma_py(imol_map)
    sigma_diff = coot.map_sigma_py(imol_diff)
    backbone   = [" N  ", " CA ", " C  ", " O  "]

    for resno in resnos:
        atoms = coot.residue_info_py(imol, chain, resno, "")
        if not isinstance(atoms, list): continue
        name = coot.residue_name_py(imol, chain, resno, "")
        coords = {a[0][0]: a[2] for a in atoms}
        for atname in backbone:
            xyz = coords.get(atname)
            if xyz is None: continue
            x, y, z = xyz
            sig_main = coot.density_at_point(imol_map,  x, y, z) / sigma
            sig_diff = coot.density_at_point(imol_diff, x, y, z) / sigma_diff
            flag = " << LOW" if sig_main < 0.8 else ""
            if sig_main < 0.0: flag = " << NEGATIVE"
            print(f"  {chain}/{resno} {name} {atname.strip():4s}  "
                  f"{sig_main:.2f}σ  diff={sig_diff:+.2f}σ{flag}")
```

### Interpreting Per-Atom Density Values

| Pattern | Diagnosis | Fix |
|---|---|---|
| O near 0σ, N/CA/C good | Peptide bond flipped — carbonyl O in empty space | `pepflip()` at this residue |
| CA near 0σ or negative | Atom genuinely in wrong position | Register shift or delete/rebuild |
| Large negative diff at CA | Atom displaced — strong negative difference density | Delete and rebuild or NCS copy |
| Large positive diff nearby | Missing atoms or unmodelled density | Build missing residues |
| All atoms < 1σ | Entire residue misplaced | Delete and rebuild |

**The pepflip signature:** near-zero or negative carbonyl O with good N, CA, C.
This is the most common backbone error and the fastest to fix.

```python
# If O < 0.5σ but CA > 2σ — it's a pepflip
coot.pepflip(imol, chain_id, resno, "", "")
specs = [[chain_id, r, ""] for r in range(resno-2, resno+3)
         if coot.residue_name_py(imol, chain_id, r, "") not in [None, "False", ""]]
coot.refine_residues_py(imol, specs)
# Re-probe to confirm O has moved into density
```

---

## Mainchain vs Sidechain Correlation: Distinguish Before Acting

**Always check mainchain and sidechain correlation separately before deciding on a fix.**

```python
# atom_mask_mode: 0=all, 1=mainchain only, 2=sidechain only
residue_specs = [[chain_id, resno, ""]]

mc = coot.map_to_model_correlation_py(imol, residue_specs, [], 1, imol_map)
sc = coot.map_to_model_correlation_py(imol, residue_specs, [], 2, imol_map)
aa = coot.map_to_model_correlation_py(imol, residue_specs, [], 0, imol_map)

print(f"All-atom: {aa:.3f}  Mainchain: {mc:.3f}  Sidechain: {sc:.3f}")
```

| Result | Conclusion | Fix |
|---|---|---|
| MC poor, SC good | Backbone error (pepflip, register) | Per-atom density probe, then pepflip |
| SC poor, MC good | Sidechain rotamer wrong | `auto_fit_best_rotamer()` |
| Both poor | Whole residue misplaced | Delete/rebuild or NCS copy |

---

## Register Shift Diagnosis

When a stretch of residues has consistently poor correlations flanked by good residues,
**always test for a register shift before attempting pepflips**.

A sequence of consecutive Ramachandran outliers with poor density for the sidechains, but
acceptable or even good coorelations for the mainchain atoms, is the hallmark of a
register error, not individual residue problems.

### How to Test for a Register Shift

```python
# Test if B[n] matches A[n-1] (chain B one residue ahead of A)
for n in range(start, end):
    nb = coot.residue_name_py(imol, "B", n, "")
    na = coot.residue_name_py(imol, "A", n-1, "")
    match = "==" if nb == na else "!="
    print(f"B/{n} {nb}  vs  A/{n-1} {na}  {match}")

# Also compare phi/psi angles — a run of 3+ matching phi/psi confirms the shift
rama = coot.all_molecule_ramachandran_score_py(imol)
per_res = [r for r in rama[5] if r != -1]
rama_b = {r[1][1]: (r[0][0], r[0][1]) for r in per_res if r[1][0] == "B"}
rama_a = {r[1][1]: (r[0][0], r[0][1]) for r in per_res if r[1][0] == "A"}

for n in range(start, end):
    phi_b, psi_b = rama_b.get(n, (0, 0))
    phi_a, psi_a = rama_a.get(n-1, (0, 0))
    dphi = min(abs(phi_b-phi_a), 360-abs(phi_b-phi_a))
    dpsi = min(abs(psi_b-psi_a), 360-abs(psi_b-psi_a))
    match = "GOOD" if dphi < 30 and dpsi < 30 else "bad"
    print(f"B/{n} vs A/{n-1}: dphi={dphi:.0f} dpsi={dpsi:.0f}  {match}")
```

### Also Verify Chain Geometry

```python
import math
def dist(a, b): return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

gap = dist(c_prev, n_next)
# ~1.33 A = direct peptide bond (no gap)
# ~3.8  A = one residue missing
# ~9.6  A = two residues missing
```

---

## Spurious Inserted Residues: Delete and Renumber

A residue with near-zero or negative correlation, missing sidechain atoms, and a
zero Ramachandran score, and highly distorted bonds and angles is almost certainly
a phantom insertion placed to compensate for a genuine deletion elsewhere in the loop.
This is called an "out of register" error - at some point in the model there will
be a residue squeezed in to where it should not be and elsewhere there will be a residue
that is stretched across where 2 residues should be.

**Signs of a spurious residue:**
- Correlation < 0.2 or negative
- Ramachandran score = 0.00000
- Missing sidechain atoms
- Flanking residues fit well
- CA→CA distance to neighbour is too short (~3.2 Å instead of ~3.8 Å)

**The correct fix is delete + renumber, NOT pepflip + refine:**

```python
coot.make_backup_checkpoint(imol, "before delete spurious residue")
coot.delete_residue(imol, chain_id, resno, "")
coot.renumber_residue_range(imol, chain_id, resno+1, last_resno, -1)

# Fix residue types if sequence labels are now wrong
for n, target_restype in corrections.items():
    coot.mutate(imol, chain_id, n, "", target_restype)

for n in affected_resnos:
    coot.auto_fit_best_rotamer(imol, chain_id, n, "", "", imol_map, 1, 0.01)

specs = [[chain_id, r, ""] for r in range(resno-2, last_affected+2)
         if coot.residue_name_py(imol, chain_id, r, "") not in [None, "False", ""]]
coot.refine_residues_py(imol, specs)
```

**After renumbering, always check:**
1. Bond geometry: C→N ~1.33 Å, CA→CA ~3.8 Å
2. Sequence alignment vs the other NCS chain
3. Correlations throughout the corrected region

---

## Anisotropic B-Factors: Always Handle Safely

B-factor fields can be lists (anisotropic) not floats (isotropic).
Always use this helper — failing to do so causes a TypeError:

```python
def get_b(atom):
    b = atom[1][1]
    return b[0] if isinstance(b, list) else b
```

---

## Chain Break Geometry: Measure Before Building

```python
import math
def dist(a, b): return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

c_prev  = get_coords(imol, chain, resno_before)['C']
n_next  = get_coords(imol, chain, resno_after)['N']
cn_dist = dist(c_prev, n_next)
print(f"C->N: {cn_dist:.2f} A")
# ~1.33 A: direct bond — no residue missing, do not attempt to build
# ~3.8  A: one residue missing
# ~9.6  A: two residues missing
```

If C→N ≈ 1.33 Å, there is no room — the deletion is genuine, the model is correct.

---

## Know When to Stop Automated Methods

When pepflips and rotamer fitting cycle without convergence after 2-3 attempts, stop.
Use the per-atom density probe to diagnose the true problem.

Oscillation (O good → CA bad → O bad → CA good) means the residue is genuinely
misplaced and automated tools cannot find the correct basin. Options:
- Use `protein_db_loops_py()` for proper fragment-database loop fitting
- Use NCS copy if a related chain is already correct (always prefer this)

**Do not use `multi_residue_torsion_fit_py` as a substitute for proper loop fitting.**
It is useful for initial placement of manually built placeholder atoms only.

---

## Loop Fitting with protein_db_loops_py

Use only when NCS copy is not available and the loop is genuinely missing (not a
register error).

```python
# Anchors must exist; missing residues need placeholder atoms first
specs = [
    [chain_id, anchor_start,   ""],  # last good residue before gap
    [chain_id, anchor_start+1, ""],  # last good residue before gap
    [chain_id, missing_1,      ""],  # built with placeholder atoms
    [chain_id, missing_2,      ""],  # built with placeholder atoms
    [chain_id, anchor_end,     ""],  # first good residue after gap
]

result = coot.protein_db_loops_py(
    imol, specs, imol_map,
    20,    # nfrags — 20 is thorough
    True   # preserve_residue_names
)
```

---

## MANDATORY: Complete Validation Workflow

**Check ALL metrics before AND after EVERY fix.**

### Before fixing ANY residue:
1. Center on it
2. Check Ramachandran, rotamer, all-atom corr, sidechain corr, overlaps
3. For backbone problems: check per-atom density with `density_at_point()`
4. For NCS structures: check if the related chain has the correct conformation

### After fixing ANY residue:
- Ramachandran probability > 0.02 (preferably > 0.1)
- Rotamer score > 1.0% (preferably > 5%)
- Density correlation > 0.7 all-atom and sidechain (preferably > 0.8)
- No severe clashes (overlap volume < 2.0 Å³)
- Backbone atoms all > 1σ in 2mFo-DFc map

### Acceptable Reasons to Stop Without Perfect Metrics:
- Per-atom probe shows atoms are in correct density despite low per-residue correlation
- Remaining issues are at gap boundaries (Ramachandran nan is expected there)
- Residue is genuinely disordered (high B-factors consistent with low correlation)
- Multiple approaches documented and best achievable state recorded

---

## Refinement Best Practices

### Make Checkpoints Before Changes

```python
checkpoint_idx = coot.make_backup_checkpoint(imol, "before pepflip A/262")
# ... make changes ...
# If worse: coot.restore_to_backup_checkpoint(imol, checkpoint_idx)
# If checkpoint fails: coot.set_undo_molecule(imol); coot.apply_undo()
```

### Extend Selection Around Problem Residues

```python
specs = [[chain_id, r, ""] for r in range(resno-2, resno+3)
         if coot.residue_name_py(imol, chain_id, r, "") not in [None, "False", ""]]
coot.refine_residues_py(imol, specs)
```

### Sphere Refinement for Spatial Context

```python
central = [chain_id, resno, ""]
neighbours = coot.residues_near_residue_py(imol, central, 4.5)
coot.refine_residues_py(imol, neighbours + [central])
```

### Set Refinement to Synchronous Mode

```python
coot.set_refinement_immediate_replacement(1)  # once at session start
coot.set_imol_refinement_map(imol_map)
```

### Delete Waters Displaced by Newly Built Residues

```python
for resno in newly_built_resnos:
    for n in coot.residues_near_residue_py(imol, [chain, resno, ""], 3.5):
        if coot.residue_name_py(imol, n[0], n[1], n[2]) == "HOH":
            coot.delete_residue(imol, n[0], n[1], n[2])
```

### Secondary Structure Restraints

```python
coot.set_secondary_structure_restraints_type(1)  # helix (phi=-57, psi=-47)
coot.set_secondary_structure_restraints_type(2)  # strand (phi=-120, psi=+120)
coot.refine_residues_py(imol, specs)
coot.set_secondary_structure_restraints_type(0)  # ALWAYS turn off afterwards
```

---

## Mainchain Geometry Diagnostic: Bond, Angle and Omega Z-scores

**Use this tool to identify distorted backbone geometry. It is fast, independent of the
refinement environment, and often reveals register errors and misplaced residues that
correlation and Ramachandran scores alone miss.**

The **omega torsion angle** (the C-CA-C-N dihedral across the peptide bond, ideal = 180°
for trans) is the most sensitive indicator. A wildly wrong omega (e.g. 9°, 60°, 123°)
is definitive proof that the backbone is misplaced — it cannot be a real conformation.

### Complete Diagnostic Function

```python
import math

def dist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

def angle_deg(a, b, c):
    ba = tuple(a[i]-b[i] for i in range(3))
    bc = tuple(c[i]-b[i] for i in range(3))
    cos_a = sum(ba[i]*bc[i] for i in range(3)) / (
        math.sqrt(sum(x*x for x in ba)) * math.sqrt(sum(x*x for x in bc)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_a))))

def omega_deg(ca1, c1, n2, ca2):
    def cross(a, b):
        return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
    def sub(a, b): return tuple(a[i]-b[i] for i in range(3))
    def norm(v):
        l = math.sqrt(sum(x*x for x in v))
        return tuple(x/l for x in v)
    def dot(a, b): return sum(a[i]*b[i] for i in range(3))
    b1 = sub(c1, ca1); b2 = sub(n2, c1); b3 = sub(ca2, n2)
    n1 = norm(cross(b1, b2)); n2v = norm(cross(b2, b3))
    m1 = cross(n1, norm(b2))
    return math.degrees(math.atan2(dot(m1, n2v), dot(n1, n2v)))

# Engh & Huber ideal values (sigma = expected standard deviation)
IDEAL = {
    'N-CA':   (1.458, 0.019),
    'CA-C':   (1.525, 0.021),
    'C-N':    (1.329, 0.014),   # peptide bond
    'C-O':    (1.231, 0.020),
    'N-CA-C': (111.2, 2.8),
    'CA-C-N': (116.2, 2.0),
    'C-N-CA': (121.7, 1.8),
    'CA-C-O': (120.8, 1.7),
}

def get_coords(imol, chain, resno):
    atoms = coot.residue_info_py(imol, chain, resno, "")
    if not isinstance(atoms, list) or atoms is False:
        return {}
    return {a[0][0].strip(): a[2] for a in atoms}

def check_mainchain_geometry(imol, chain, res_start, res_end, z_threshold=3.0):
    """
    Check mainchain bond lengths, angles and omega torsions against Engh & Huber ideals.
    Returns list of (residue_tag, parameter, ideal, actual, z_score, kind) tuples,
    sorted by |Z| descending.
    """
    distortions = []
    prev = None

    for resno in range(res_start, res_end + 1):
        name = coot.residue_name_py(imol, chain, resno, "")
        if not name or name in ["False", "", False]:
            prev = None
            continue
        c = get_coords(imol, chain, resno)
        if not c or "CA" not in c:
            prev = None
            continue

        tag = f"{chain}/{resno} {name}"

        # Intra-residue bonds
        for bond, (a1, a2) in [("N-CA",("N","CA")), ("CA-C",("CA","C")), ("C-O",("C","O"))]:
            if a1 in c and a2 in c:
                d = dist(c[a1], c[a2])
                ideal, sigma = IDEAL[bond]
                z = (d - ideal) / sigma
                if abs(z) >= z_threshold:
                    distortions.append((tag, bond, ideal, d, z, "bond"))

        # Intra-residue angles
        for ang_name, (a1,a2,a3) in [("N-CA-C",("N","CA","C")), ("CA-C-O",("CA","C","O"))]:
            if all(a in c for a in (a1,a2,a3)):
                ang = angle_deg(c[a1], c[a2], c[a3])
                ideal, sigma = IDEAL[ang_name]
                z = (ang - ideal) / sigma
                if abs(z) >= z_threshold:
                    distortions.append((tag, ang_name, ideal, ang, z, "angle"))

        # Inter-residue geometry (peptide bond to previous residue)
        if prev:
            prev_name = coot.residue_name_py(imol, chain, resno-1, "")
            ptag = f"{chain}/{resno-1} {prev_name}→{tag}"

            if "C" in prev and "N" in c:
                # Peptide C-N bond length
                d = dist(prev["C"], c["N"])
                ideal, sigma = IDEAL["C-N"]
                z = (d - ideal) / sigma
                if abs(z) >= z_threshold:
                    distortions.append((ptag, "C-N(pept)", ideal, d, z, "bond"))

                # Angles across the peptide bond
                if "CA" in prev:
                    ang = angle_deg(prev["CA"], prev["C"], c["N"])
                    ideal, sigma = IDEAL["CA-C-N"]
                    z = (ang - ideal) / sigma
                    if abs(z) >= z_threshold:
                        distortions.append((ptag, "CA-C-N", ideal, ang, z, "angle"))

                if "CA" in c:
                    ang = angle_deg(prev["C"], c["N"], c["CA"])
                    ideal, sigma = IDEAL["C-N-CA"]
                    z = (ang - ideal) / sigma
                    if abs(z) >= z_threshold:
                        distortions.append((ptag, "C-N-CA", ideal, ang, z, "angle"))

                    # Omega torsion — flag deviation from 180° > 15°
                    if "CA" in prev:
                        om = omega_deg(prev["CA"], prev["C"], c["N"], c["CA"])
                        dev = abs(abs(om) - 180.0)
                        if dev > 15.0:
                            # Express as Z using sigma~5° for trans peptides
                            distortions.append((ptag, "omega", 180.0, om, dev/5.0, "torsion"))

        prev = c

    distortions.sort(key=lambda x: -abs(x[4]))
    return distortions
```

### Usage

```python
# Check a region — results sorted worst first
dist_b = check_mainchain_geometry(imol, "B", 256, 268, z_threshold=3.0)

for res, param, ideal, actual, z, kind in dist_b:
    unit = "Å" if kind == "bond" else "°"
    flag = " ***" if abs(z) > 5 else ""
    print(f"  {res:35s}  {param:12s}  ideal={ideal:.3f}{unit}  actual={actual:.3f}{unit}  z={z:+.2f}{flag}")
```

### Interpreting the Results

| Distortion | Likely cause | Fix |
|---|---|---|
| omega far from 180° (> 20°) | Misplaced backbone / register error | Delete + renumber or NCS copy |
| omega near 0° (cis peptide) | May be genuine (rare) — check density | Visual inspection; pepflip if in empty density |
| N-CA-C angle severely wrong | Atom in completely wrong position | Per-atom density probe, then delete/rebuild |
| N-CA or CA-C bond short | Compressed geometry from register error | Delete + renumber |
| CA-C-N or C-N-CA angle wrong | Distorted peptide junction | Pepflip or rebuild |

**Omega outliers are the fastest route to finding register errors.** A run of 3+ consecutive
omega outliers in the same region is essentially diagnostic of a misplaced loop.

### Always render the results as an interactive SVG widget

After running `check_mainchain_geometry()`, render the results using `visualize:show_widget`
so each distorted residue is a clickable block that navigates to it or triggers a fix.
Use `c-red` for |Z| > 5 or omega deviation > 20°, `c-amber` for |Z| 3–5. See the
`coot-essential-api` skill for the widget rendering guidelines.

---

## Key Takeaways

**For NCS structures: always check for a related chain that is correctly built before
attempting manual rebuilding. One NCS copy call beats dozens of pepflips.**

**Use per-atom density probing before pepflipping.** The pattern of which backbone
atoms are in low density tells you exactly what is wrong and what fix to apply.

**Distinguish mainchain from sidechain problems** using `atom_mask_mode=1` and `2`
before deciding on a fix strategy.

**Measure C→N and CA→CA distances before building.** If C→N ≈ 1.33 Å there is no
gap and no room to build — the deletion is genuine.

**Context matters in refinement.** Include neighbouring residues, both sequence and
spatial neighbours, for best results.

**Know when to stop.** If automated methods cycle without convergence after 2-3
attempts, use the density probe to decide whether to delete, rebuild, or accept.


## Key Takeaway

**Context matters in refinement.** Including neighboring residues provides the geometric and density context needed for refinement algorithms to find better solutions, especially for severe outliers.

**Always checkpoint before changes.** Use `make_backup_checkpoint()` before any significant model modification so you can easily revert if needed.

**Delete displaced waters immediately** after building new residues — don't wait until refinement reveals clashes.

**Always turn off secondary structure restraints** (`coot.set_secondary_structure_restraints_type(0)`) after the refinement that uses them, so they don't inadvertently affect subsequent work.
