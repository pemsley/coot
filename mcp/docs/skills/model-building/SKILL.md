---
name: coot-model-building
description: "Best Practices for Model-Building Tools and Refinement"
---

# Key Lessons from Chain A Refinement Session

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

### 2. Always Accept Regularizement After Refining
**Critical:** Call `coot.accept_moving_atoms_py()` after every `refine_residues_py()` call.
- Without this, subsequent refinements may fail silently
- The refinement creates intermediate atoms that must be accepted to update the model
```python
results = coot.refine_residues_py(0, [["A", 41, ""]])
```

### 3. Extend Selection Around Problem Residues
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

### 3a. Include Spatial Neighbours, Not Just Sequence Neighbours

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
coot.refine_residues_py(0, [["A", 88, ""], ["A", 89, ""], ["A", 90, ""]])
coot.accept_moving_atoms_py()

# Then re-refine A/2 INCLUDING A/89 as a spatial neighbour:
coot.refine_residues_py(0, [["A", 1, ""], ["A", 2, ""], ["A", 3, ""], ["A", 89, ""]])
coot.accept_moving_atoms_py()
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
coot.refine_residues_py(0, [["A", i, ""] for i in range(40, 44)])
coot.accept_moving_atoms_py()
check_validation()  # Did it help?

# Second: targeted refinement
coot.refine_residues_py(0, [["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])
coot.accept_moving_atoms_py()
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

## Key Takeaway

**Context matters in refinement.** Including neighboring residues provides the geometric and density context needed for refinement algorithms to find better solutions, especially for severe outliers.

**Always checkpoint before changes.** Use `make_backup_checkpoint()` before any significant model modification so you can easily revert if needed.
