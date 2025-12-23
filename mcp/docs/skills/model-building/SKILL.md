---
name: coot-model-building
description: "Best Practices for Model-Building Tools and Refinement"
---

# Key Lessons from Chain A Refinement Session

## Refinement Best Practices

### 1. Always Accept Regularizement After Refining
**Critical:** Call `coot.accept_regularizement()` after every `refine_residues_py()` call.
- Without this, subsequent refinements may fail silently
- The refinement creates intermediate atoms that must be accepted to update the model
```python
coot.refine_residues_py(0, [["A", 41, ""]])
coot.accept_regularizement()  # Essential!
```

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

### 3. Iterative Refinement Strategy
Sometimes multiple rounds of refinement with different selections help:

1. **First pass:** Refine larger region to establish general geometry
2. **Second pass:** Refine smaller region to fine-tune specific problem
3. **Check validation** after each step
4. **Undo** if results get worse

**Example workflow:**
```python
# First: larger region
coot.refine_residues_py(0, [["A", i, ""] for i in range(40, 44)])
coot.accept_regularizement()
check_validation()  # Did it help?

# Second: targeted refinement
coot.refine_residues_py(0, [["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])
coot.accept_regularizement()
check_validation()  # Better or worse?

# If worse:
coot.apply_undo()
```

### 4. Measure Before and After
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

# Usage
before = check_residue_validation(0, "A", 41)
coot.refine_residues_py(0, [["A", 40, ""], ["A", 41, ""], ["A", 42, ""], ["A", 43, ""]])
coot.accept_regularizement()
after = check_residue_validation(0, "A", 41)

# Compare and decide
if after['all_atom_corr'] > before['all_atom_corr']:
    # Keep it!
    pass
else:
    # Revert
    coot.apply_undo()
```

### 5. Use Undo Liberally
**Don't be afraid to revert changes** - `apply_undo()` is your friend.

- If validation metrics get worse, undo immediately
- If refinement creates new problems, undo
- You can try multiple approaches and keep the best result

### 6. Auto-fit Rotamer for Side-chain Issues
**For poor side-chain density correlation**, try `auto_fit_best_rotamer()` first:
```python
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
            coot.apply_undo()
    else:
        # Negative score means failure
        coot.apply_undo()
```

### 7. Set Refinement to Synchronous Mode
**Always call this at the start** to make refinement complete immediately:
```python
coot.set_refinement_immediate_replacement(1)
```

Without this, refinement may be asynchronous and difficult to control programmatically.

### 8. Navigate to Residue Before Working
**Bring residue to screen center** so you can watch the refinement:
```python
coot.set_go_to_atom_molecule(0)
coot.set_go_to_atom_chain_residue_atom_name("A", 41, "CA")
```

This helps with:
- Visual inspection of the problem
- Seeing the refinement in real-time
- Verifying the result makes geometric sense

### 9. Flipping peptides

If the Ramachandran Plot is poor, try using `coot.pepflip(imol, chain_id, res_no, ins_code, alt_conf)` followed by a refinement of the residues in the extended region.

### 10. Flipping side-chains terminal Chi-angle

If the Rotamer score is poor, try using `coot.do_180_degree_side_chain_flip()` to improve the Rotamer score. It is occassionally useful.

## Key Takeaway

**Context matters in refinement.** Including neighboring residues provides the geometric and density context needed for refinement algorithms to find better solutions, especially for severe outliers.
