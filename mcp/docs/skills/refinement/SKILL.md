---
name: coot-refinement
description: "Best practices for protein structure refinement and validation in Coot. Use when performing (1) Residue refinement operations, (2) Model building and fitting, (3) Rotamer fixing, (4) Scripted/automated refinement workflows, (5) Validation and correlation checking."
---

# Coot Refinement Best Practices

This skill provides guidance for effective and safe structure refinement in Coot, based on lessons learned from crashes and workflow optimization.

## Critical Configuration

### Set Immediate Replacement Mode for Scripting

**ALWAYS** call this before any refinement operations in scripts:

```python
coot.set_refinement_immediate_replacement(1)
```

**Why:** Enables synchronous operation where refinement directly updates coordinates. Without this:
- Refinement happens asynchronously in background threads
- Need to call `coot.set_refinement_immediate_replacement()` before using refinement functions (just once is enough)
  which should remove the risk of threading conflicts and crashes, race conditions between refinement and rendering.

**When to use:** Any time you're scripting refinement operations (refine_residues_py, refine_zone, etc.)

## Safe Refinement Workflow

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
neigbs = coot.residues_near_residue(imol, central_residue_spec)
residue_spec = neigbs
residue_specs.append(central_residue_spec)
result = coot.refine_residues_py(imol, residue_specs)

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
# Get density correlation for specific residue
correlation = coot.density_correlation_analysis_scm(imol, chain_id, resno, ins_code)
# Returns dict with 'all-atom' and 'side-chain' correlations

# Get worst residues
worst = coot.get_n_residues_with_worst_density_fit(imol, n_residues)
# Returns list of [chain_id, resno, inscode, correlation]
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

```python
# 1. Enable immediate replacement
coot.set_refinement_immediate_replacement(1)

# 2. Find worst residues
worst = coot.get_n_residues_with_worst_density_fit(0, 10)

# 3. For each poor residue:
for residue in worst:
    chain_id, resno, inscode, corr = residue

    # Build CID
    cid = f"//{chain_id}/{resno}"

    # Try rotamer fix
    coot.auto_fit_best_rotamer(cid, "", 0, 1, 1, 0.1)

    # Refine in context
    coot.refine_residues_using_atom_cid(0, cid, "SPHERE", 4000)

    # Check improvement
    new_corr = coot.density_correlation_analysis_scm(0, chain_id, resno, inscode)
    print(f"{cid}: {corr:.3f} → {new_corr['all-atom']:.3f}")
```

## Common Functions

### Refinement
- `auto_fit_best_rotamer(cid, alt_conf, imol, imol_map, use_rama, rama_weight)` - Fix rotamer
- `refine_residues_using_atom_cid(imol, cid, mode, radius)` - Sphere/zone refinement
- `refine_zone(imol, chain, start, end, alt_conf)` - Refine residue range
- `set_refinement_immediate_replacement(istate)` - Enable synchronous refinement

### Validation
- `density_correlation_analysis_scm(imol, chain, resno, inscode)` - Get correlation
- `get_n_residues_with_worst_density_fit(imol, n)` - Find problem residues

### Other
- `pepflip(imol, atom_cid, alt_conf)` - Flip peptide
