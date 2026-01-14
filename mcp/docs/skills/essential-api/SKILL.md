---
name: coot-essential-api
description: "API documentation to be loaded at startup"
---
# Coot Essential API Functions

This document contains the core Coot API functions needed for typical validation and model-building workflows. Reading these function signatures at session start eliminates the need for searching.

## Setup & Configuration

```python
coot.set_refinement_immediate_replacement(1)
# CRITICAL: Call this before any refinement operations to make them synchronous
# Without this, refinement results may not be available immediately
```

## Molecule Management

```python
coot.is_valid_model_molecule(imol) -> int  # Returns 1 if valid model, 0 otherwise
coot.is_valid_map_molecule(imol) -> int    # Returns 1 if valid map, 0 otherwise
coot.molecule_name(imol) -> str            # Returns the molecule filename/description
coot.n_chains(imol) -> int                 # Returns number of chains
coot.coot_version() -> str                 # Returns Coot version string
coot.load_tutorial_model_and_data()        # Loads tutorial RNase structure + maps
```

## Chain and Residue Information

```python
# Requires: import coot_utils
coot_utils.chain_ids(imol) -> list         # Returns list of chain IDs, e.g., ['A', 'B']

# Direct C++ functions (preferred when possible)
coot.chain_id_py(imol, chain_index) -> str # Get chain ID by index
```

## Navigation

```python
coot.set_go_to_atom_chain_residue_atom_name(chain_id, resno, atom_name) -> int
# Centers view on specified atom. Returns 1 on success.
# Example: coot.set_go_to_atom_chain_residue_atom_name("A", 42, "CA")

coot.closest_atom_simple_py() -> list
# Returns: [imol, chain_id, resno, ins_code, atom_name, alt_conf]
# Gets the atom closest to screen center across all displayed molecules

coot.closest_atom_py(imol) -> list
# Same as above but for specific molecule

coot.active_atom_spec_py() -> list
# Returns the currently "active" atom specification
```

## Residue Inspection

```python
coot.residue_info_py(imol, chain_id, resno, ins_code) -> list
# Returns detailed atom information for a residue
#
# Parameters:
#   imol: Model molecule index
#   chain_id: Chain identifier (e.g., "A")
#   resno: Residue number
#   ins_code: Insertion code (use "" if none)
#
# Returns: List of atom entries, each containing:
#   [[atom_name, alt_conf], [occupancy, b_factor, element, ?], [x, y, z], atom_index]
#
# Example output for a complete CYS:
#   [[' N  ', ''], [1.0, 12.5, ' N', ''], [x, y, z], 100],
#   [[' CA ', ''], [1.0, 11.2, ' C', ''], [x, y, z], 101],
#   [[' CB ', ''], [1.0, 14.3, ' C', ''], [x, y, z], 102],
#   [[' SG ', ''], [1.0, 18.1, ' S', ''], [x, y, z], 103],  # Sulfur!
#   [[' C  ', ''], [1.0, 10.8, ' C', ''], [x, y, z], 104],
#   [[' O  ', ''], [1.0, 11.0, ' O', ''], [x, y, z], 105]

# Check for missing atoms in a residue
atoms = coot.residue_info_py(0, "A", 72, "")
atom_names = [a[0][0].strip() for a in atoms]
print(f"Atoms present: {atom_names}")

# Expected atoms for common residues
expected_atoms = {
    'CYS': ['N', 'CA', 'CB', 'SG', 'C', 'O'],
    'ILE': ['N', 'CA', 'CB', 'CG1', 'CG2', 'CD1', 'C', 'O'],
    'GLY': ['N', 'CA', 'C', 'O'],
    'PHE': ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'C', 'O'],
}

# Find missing atoms
res_type = coot.residue_name(0, "A", 72, "")
if res_type in expected_atoms:
    missing = [a for a in expected_atoms[res_type] if a not in atom_names]
    if missing:
        print(f"WARNING: Missing atoms in {res_type}: {missing}")
```

## Validation - Density Fit

```python
coot.map_to_model_correlation_stats_per_residue_range_py(
    imol,           # Model molecule number
    chain_id,       # Chain identifier (e.g., "A")
    imol_map,       # Map molecule number
    n_per_range,    # Residues per window (use 1 for per-residue)
    exclude_NOC     # 0=include all atoms, 1=exclude backbone N,O,C
) -> list
# Returns: [[all_atom_stats], [sidechain_stats]]
# Each stats list: [[residue_spec, [n_points, correlation]], ...]
# residue_spec = [chain_id, resno, ins_code]

# Example - find worst fitting residues:
stats = coot.map_to_model_correlation_stats_per_residue_range_py(0, "A", 1, 1, 0)
all_atom = stats[0]
worst = sorted(all_atom, key=lambda x: x[1][1])[:5]  # 5 worst by correlation
```

## Validation - Geometry

```python
coot.all_molecule_ramachandran_score_py(imol) -> list
# Returns: [score, n_residues, ..., per_residue_data]
# per_residue_data: [[[phi, psi], residue_spec, probability, [prev, curr, next]], ...]
# LOW probability = BAD (outlier)

coot.rotamer_graphs_py(imol) -> list
# Returns: [[chain_id, resno, ins_code, score_percentage, resname], ...]
# LOW score = BAD rotamer

coot.molecule_atom_overlaps_py(imol, n_pairs) -> list
# Returns worst n_pairs atom overlaps (use -1 for all)
# Each overlap: {
#   'atom-1-spec': [imol, chain, resno, ins, atom_name, alt],
#   'atom-2-spec': [imol, chain, resno, ins, atom_name, alt],
#   'overlap-volume': float  # in Å³, >5.0 is severe
# }
```

## Validation - Unmodeled Density

```python
coot.find_blobs_py(imol_model, imol_map, sigma_cutoff) -> list
# Finds unmodeled density blobs
# Returns: [[position, score], ...]
# position has .x(), .y(), .z() methods
# Use sigma_cutoff=3.0 for difference maps, 1.0 for 2mFo-DFc
# Higher score = larger/stronger blob
```

## Refinement

```python
coot.refine_residues_py(imol, residue_specs) -> list
# Real-space refinement of specified residues
# residue_specs = [["A", 42, ""], ["A", 43, ""], ...]  # [chain, resno, ins_code]
# Returns refinement status

coot.accept_moving_atoms_py() -> list
# Accept the current refinement/regularization
# Returns: ['', status, [[metric_name, description, value], ...]]
# Call this after refine_residues_py()
```

## Model Building - Rotamers

```python
coot.auto_fit_best_rotamer(
    imol,              # Model molecule
    chain_id,          # Chain (e.g., "A")
    resno,             # Residue number
    ins_code,          # Insertion code (usually "")
    altloc,            # Alt conf (usually "")
    imol_map,          # Map for density scoring
    clash_flag,        # 1=check clashes, 0=ignore
    lowest_probability # Minimum rotamer probability (e.g., 0.01)
) -> float
# Returns new rotamer score, or -99.9 if residue has no rotamers (GLY, ALA)
```

## Model Building - Backbone

```python
coot.pepflip(imol, chain_id, resno, ins_code, altloc)
# Flips the peptide bond at specified residue
# Use for fixing cis/trans peptide issues or Ramachandran outliers
# or other false minimum backbone conformations.
# Follow with refinement of surrounding residues
```

## Typical Validation & Fix Workflow

```python
# 1. Setup
coot.set_refinement_immediate_replacement(1)

# 2. Check what's loaded
for i in range(10):
    if coot.is_valid_model_molecule(i):
        print(f"Model {i}: {coot.molecule_name(i)}")
    if coot.is_valid_map_molecule(i):
        print(f"Map {i}: {coot.molecule_name(i)}")

# 3. Validate density fit
stats = coot.map_to_model_correlation_stats_per_residue_range_py(0, "A", 1, 1, 0)
worst = sorted(stats[0], key=lambda x: x[1][1])[:10]

# 4. Check for clashes
overlaps = coot.molecule_atom_overlaps_py(0, 30)
severe = [o for o in overlaps if o['overlap-volume'] > 5.0]

# 5. Fix bad rotamers
coot.auto_fit_best_rotamer(0, "A", 89, "", "", 1, 1, 0.01)
coot.refine_residues_py(0, [["A", 89, ""]])
coot.accept_moving_atoms_py()

# 6. Fix backbone issues
coot.pepflip(0, "A", 41, "", "")
coot.refine_residues_py(0, [["A", 40, ""], ["A", 41, ""], ["A", 42, ""]])
coot.accept_moving_atoms_py()

# 7. Re-validate
overlaps_after = coot.molecule_atom_overlaps_py(0, 10)
```

## Important Notes

1. **Always call `set_refinement_immediate_replacement(1)` first** - makes refinement synchronous
2. **Always call `accept_moving_atoms_py()` after refinement** - commits the changes
3. **Use `coot.*_py()` functions directly** - faster than `coot_utils` wrappers
4. **Import coot_utils only when needed** - for convenience functions like `chain_ids()`
5. **The `coot` module is auto-imported** - no import statement needed
