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


## Key Takeaway

**Context matters in refinement.** Including neighboring residues provides the geometric and density context needed for refinement algorithms to find better solutions, especially for severe outliers.

**Always checkpoint before changes.** Use `make_backup_checkpoint()` before any significant model modification so you can easily revert if needed.

**Delete displaced waters immediately** after building new residues — don't wait until refinement reveals clashes.

**Always turn off secondary structure restraints** (`coot.set_secondary_structure_restraints_type(0)`) after the refinement that uses them, so they don't inadvertently affect subsequent work.
