---
name: coot-refinement
description: "Best practices for protein structure refinement and validation in Coot. Use when performing (1) Residue refinement operations, (2) Model building and fitting, (3) Rotamer fixing, (4) Scripted/automated refinement workflows, (5) Validation and correlation checking."
---

# Coot Refinement Best Practices

This skill provides guidance for effective and safe structure refinement in Coot, based on lessons learned from crashes and workflow optimization.

## Critical Configuration

### Set Immediate Replacement Mode for Scripting

**ALWAYS** call these two before any refinement operations in scripts:

```python
coot.set_refinement_immediate_replacement(1)
# Makes refinement synchronous — results are committed immediately.

coot.set_imol_refinement_map(imol_map)
# Tells Coot which map to refine against.
# CRITICAL: without this, refine_residues_py() silently fails (returns -2 with no atoms moved).
# Call once per session, or whenever the active map changes.
```

**Why:** Enables synchronous operation where refinement directly updates coordinates. Without this:
- Refinement happens asynchronously in background threads
- Need to call `coot.set_refinement_immediate_replacement()` before using refinement functions (just once is enough)
  which should remove the risk of threading conflicts and crashes, race conditions between refinement and rendering.

**When to use:** Any time you're scripting refinement operations (refine_residues_py, refine_zone, etc.)

## Safe Refinement Workflow

```python
# 1. Enable immediate replacement
coot.set_refinement_immediate_replacement(1)

# 2. Set the refinement map
coot.set_imol_refinement_map(imol_map)

# 3. Refine residues (preferred method)
# Pass a list of residue specs: [["A", 42, ""], ["A", 43, ""], ...]
residue_specs = [["A", resno, ""] for resno in range(start_resno, end_resno + 1)]
result = coot.refine_residues_py(imol, residue_specs)

# 4. so-called "sphere refine" is often useful because it allows movement/improvement
# of residues that are close in space but distant in sequence.
central_residue_spec = ['A', 12, '']
neighbours = coot.residues_near_residue_py(imol, central_residue_spec, 4.5)  # radius in Angstroms
residue_specs = neighbours + [central_residue_spec]
result = coot.refine_residues_py(imol, residue_specs)
```

### Hard-Mode-Soft-Mode-Hard-Mode

Some problems are of the kind where the geometry restraints "get in the way" of the
atoms moving to the correct minimum. In such cases, you can try letting the refinement
"go soft" - ie.

 o Mediocre Refinement
 o increase the map weight (weight_matrix) by a factor of 10 or so.
 o Refine the same residues again
 o Restore the previous map weight
 o Refine the same residues again

This can sometimes be helpful to remove bad geometry or bad fit - atom clashes in particular.

Additionally, or perhaps at the same time, you can try to reduce LJ epsilon
using coot.set_refinement_lennard_jones_epsilon() by a factor of 100 or so.
And then restore it and re-refine the same residues.

### Zone Refinement

```python
# Enable immediate replacement
coot.set_refinement_immediate_replacement(1)

# Refine zone
coot.refine_zone(
    imol,
    chain_id,
    start_resno,
    end_resno,
    ""  # alt conf
)
```

## Validation and Correlation

### Check Residue Quality

```python
# Per-residue density correlation for a chain (1 residue per range, all atoms)
all_atom, _sidechain = coot.map_to_model_correlation_stats_per_residue_range_py(
    imol, chain_id, imol_map, 1, 0)
# Each entry is [[chain_id, resno, ins_code], [n_points, correlation]]

# Worst-fitting residues: lowest all-atom correlation
# (GLY correlations are unreliable - see the coot-correlations skill)
def worst_residues(imol, imol_map, chain_id, n):
    """Return the n residues with the lowest all-atom density correlation."""
    all_atom, _sidechain = coot.map_to_model_correlation_stats_per_residue_range_py(
        imol, chain_id, imol_map, 1, 0)
    # drop empty entries (no grid points) and NaN
    scored = [(spec, corr) for spec, (n_points, corr) in all_atom
              if n_points > 0 and corr == corr]
    return sorted(scored, key=lambda sc: sc[1])[:n]

# Returns a list of ([chain_id, resno, ins_code], correlation), worst first
worst = worst_residues(imol, imol_map, chain_id, n_residues)
```

### Correlation Targets

Good fit:
- All-atom correlation: > 0.8
- Side-chain correlation: > 0.8

Poor fit (needs attention):
- All-atom correlation: < 0.5
- Side-chain correlation: < 0.5

## Common Crash Scenarios to Avoid

### Threading/Rendering Conflicts

**Problem:** Rapid-fire refinement operations while graphics are rendering can cause memory corruption in GTK rendering pipeline.

**Symptoms:**
```
nanov2_guard_corruption_detected
gdk_gl_texture_new_from_builder
gtk_gl_area_snapshot
```

**Solution:** Use `set_refinement_immediate_replacement(1)` for synchronous operation

### Asynchronous Refinement Without Accept

## Handling refine_residues_py() Return Values

`refine_residues_py()` returns `['', status, lights]` where status is a GSL minimiser code:
- **0** (`GSL_SUCCESS`): converged — done.
- **-2** (`GSL_CONTINUE`): not yet converged — **call refine_residues_py() again** (once or twice more as needed).
- **27** (`GSL_ENOPROG`): no progress — stop, refinement is stuck.

```python
# Correct pattern for robust refinement:
for _ in range(3):
    result = coot.refine_residues_py(imol, residue_specs)
    if result and result[1] != -2:
        break
lights = result[2] if result else []
```

## Multi-Line Code Limitation

Coot only returns values if code is a single line:

**Doesn't work:**
```python
x = 5
y = 10
x + y  # Won't return value
```

**Workaround:**
```python
# Call 1: Define function
def calculate():
    x = 5
    y = 10
    return x + y

# Call 2: Execute function
calculate()  # Returns 15
```

## Systematic Residue Fixing Workflow

This uses `worst_residues()` from the "Check Residue Quality" section above.

```python
def residue_correlation(imol, imol_map, residue_spec):
    """All-atom correlation of one residue (from the per-residue statistics)."""
    chain_id, resno, ins_code = residue_spec
    all_atom, _sidechain = coot.map_to_model_correlation_stats_per_residue_range_py(
        imol, chain_id, imol_map, 1, 0)
    for spec, (n_points, corr) in all_atom:
        if spec == [chain_id, resno, ins_code]:
            return corr
    return None

# 1. Enable immediate replacement and set the refinement map
coot.set_refinement_immediate_replacement(1)
coot.set_imol_refinement_map(imol_map)

# 2. Checkpoint (restore_to_backup_checkpoint(imol, checkpoint) if it gets worse),
#    then find the worst residues
checkpoint = coot.make_backup_checkpoint(imol, "before fixing worst residues")
worst = worst_residues(imol, imol_map, "A", 10)

# 3. For each poor residue:
for residue_spec, corr in worst:
    chain_id, resno, ins_code = residue_spec

    # Try a rotamer fix
    coot.auto_fit_best_rotamer(imol, chain_id, resno, ins_code, "", imol_map, 1, 0.01)

    # Refine in context (sphere refinement)
    neighbours = coot.residues_near_residue_py(imol, residue_spec, 4.5)
    coot.refine_residues_py(imol, neighbours + [residue_spec])

    # Check improvement (same measure as before)
    new_corr = residue_correlation(imol, imol_map, residue_spec)
    print(f"{chain_id}/{resno}: {corr:.3f} → {new_corr:.3f}")
```

## Refining with Secondary Structure Restraints

When adding residues in a known secondary structure conformation (e.g. a helix
or strand), use `set_secondary_structure_restraints_type()` to maintain that
geometry during real-space refinement. Without this, `refine_residues_py()` 
treats residues independently and the conformation can distort away from the
intended geometry if the density doesn't strongly support it.

### Secondary structure restraint types
- `0` — no restraints (default)
- `1` — alpha helix (restrains i→i+4 hydrogen bond geometry)
- `2` — beta strand

### Pattern: always bracket the refinement call
```python
# Set BEFORE refinement
coot.set_secondary_structure_restraints_type(1)  # 1 = alpha helix

coot.refine_residues_py(imol, residue_specs)

# Reset AFTER refinement - critical, or all subsequent refinements
# will use helix restraints unintentionally
coot.set_secondary_structure_restraints_type(0)
```

### Example: add and refine a 6-residue helix
```python
coot.set_refinement_immediate_replacement(1)
coot.set_secondary_structure_restraints_type(1)

residue_specs = [["A", resno + i, ""] for i in range(6)]
coot.refine_residues_py(imol, residue_specs)

coot.set_secondary_structure_restraints_type(0)
```

### Common mistake
Calling `refine_residues_py()` directly after building helical residues
without setting the restraint type — the helix geometry will not be
maintained during refinement.

## Common Functions

### Refinement
- `auto_fit_best_rotamer(imol, chain_id, resno, ins_code, alt_conf, imol_map, clash_flag, lowest_probability)` - Fix rotamer
- `refine_residues_py(imol, residue_specs)` - Refine a list of residue specs (add `residues_near_residue_py(imol, spec, radius)` for "sphere" refinement)
- `refine_zone(imol, chain_id, start_resno, end_resno, alt_conf)` - Refine residue range
- `set_refinement_immediate_replacement(istate)` - Enable synchronous refinement

### Validation
- `map_to_model_correlation_stats_per_residue_range_py(imol, chain_id, imol_map, n_per_range, exclude_NOC)` - Per-residue density correlation for a chain

### Other
- `pepflip(imol, chain_id, resno, ins_code, alt_conf)` - Flip peptide
