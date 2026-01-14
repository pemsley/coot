---
name: coot-correlations
description: "Using Density-Fit Correlations in Coot"
---

# Improved MCP Function Documentation for Coot

## Map to Model Correlation Functions

### Overview
These functions are used to assess how well a molecular model fits into electron density maps. They are essential for model validation and identifying poorly-fitted regions.

---

### `map_to_model_correlation_stats_per_residue_range_py()`

**Purpose**: Find residues with poor density fit by analyzing correlation statistics across the entire chain.

**Use Case**: "Which side chain is the worst fitting to density?" - Use this function to get comprehensive density fit statistics for all residues.

**Function Signature**:
```python
PyObject *map_to_model_correlation_stats_per_residue_range_py(
    int imol,                                  # Model molecule number
    const std::string &chain_id,               # Chain identifier
    int imol_map,                              # Map molecule number  
    unsigned int n_residue_per_residue_range,  # Number of residues per range (typically 1)
    short int exclude_NOC_flag                 # Exclude N, O, C atoms (1=yes, 0=no)
)
```

**Parameters**:
- `imol`: The model molecule index
- `chain_id`: The chain identifier (e.g., "A", "B")
- `imol_map`: The map molecule index to correlate against
- `n_residue_per_residue_range`: Window size for averaging (use 1 for per-residue stats)
- `exclude_NOC_flag`: Whether to exclude backbone N, O, C atoms (1 to exclude, 0 to include)

**IMPORTANT BUG**: When `exclude_NOC_flag=0`, the side-chain correlation array (stats[1]) returns all zeros/NaN values. To get valid side-chain correlations, you must use `exclude_NOC_flag=1`.

**Returns**: 
A list with two elements:
- `stats[0]`: All-atom correlations - `[[chain, resno, ins_code], [n_points, correlation]]`
- `stats[1]`: Side-chain correlations (only valid when `exclude_NOC_flag=1`)

**Example Usage**:
```python
# Get per-residue correlation stats for chain A
# Use exclude_NOC_flag=1 to get valid side-chain correlations
stats = map_to_model_correlation_stats_per_residue_range_py(
    imol=0,           # Model molecule
    chain_id="A",     # Chain A
    imol_map=1,       # Map molecule
    n_residue_per_residue_range=1,  # Per-residue statistics
    exclude_NOC_flag=1  # MUST be 1 for valid side-chain correlations
)

all_atom = stats[0]      # All-atom correlations
sidechain = stats[1]     # Side-chain correlations

# Find residues with poor all-atom correlation
for res in all_atom:
    resno = res[0][1]
    n_points = res[1][0]
    corr = res[1][1]
    if corr < 0.7:
        print(f"Poor fit: residue {resno}, correlation={corr:.3f}")
```

---

### `map_to_model_correlation_py()`

**Purpose**: Calculate the overall correlation for specific residues and their neighbors.

**Use Case**: Evaluate the fit of a specific region after refinement.

**Function Signature**:
```python
PyObject *map_to_model_correlation_py(
    int imol,
    PyObject *residue_specs,        # List of residue specs to evaluate
    PyObject *neighb_residue_specs, # Neighboring residues to exclude from grid
    unsigned short int atom_mask_mode,  # Which atoms to include (see below)
    int imol_map
)
```

**Atom Mask Modes**:
- `0`: All atoms
- `1`: Main-chain atoms if standard amino acid, else all atoms
- `2`: Side-chain atoms if standard amino acid, else all atoms  
- `3`: Side-chain atoms excluding CB if standard amino acid, else all atoms
- `4`: Main-chain atoms if standard amino acid, else nothing
- `5`: Side-chain atoms if standard amino acid, else nothing
- `10`: Atom radius dependent on B-factor

**Returns**: Float - correlation coefficient between model and map

**Example Usage**:
```python
# Evaluate side-chain fit for residues 40-44
residue_specs = [[chain_id, res_no, ins_code] for res_no in range(40, 45)]
correlation = map_to_model_correlation_py(
    imol=1,
    residue_specs=residue_specs,
    neighb_residue_specs=[],  # No neighbors to exclude
    atom_mask_mode=2,  # Side-chain atoms only
    imol_map=2
)
print(f"Side-chain correlation: {correlation}")
```

---

### `map_to_model_correlation_stats_py()`

**Purpose**: Get detailed statistics (mean, std dev, etc.) for map-model correlation.

**Function Signature**:
```python
PyObject *map_to_model_correlation_stats_py(
    int imol,
    PyObject *residue_specs,
    PyObject *neighb_residue_specs,
    unsigned short int atom_mask_mode,
    int imol_map
)
```

**Returns**: Statistics object with mean, std dev, min, max correlation values

---

### `map_to_model_correlation_per_residue_py()`

**Purpose**: Get correlation values individually for each specified residue.

**Function Signature**:
```python
PyObject *map_to_model_correlation_per_residue_py(
    int imol,
    PyObject *residue_specs,
    unsigned short int atom_mask_mode,
    int imol_map
)
```

**Returns**: List of (residue_spec, correlation) pairs

**Example Usage**:
```python
# Get per-residue correlations for a chain
residues = get_residues_in_chain_py(imol=1, chain_id="A")
correlations = map_to_model_correlation_per_residue_py(
    imol=1,
    residue_specs=residues,
    atom_mask_mode=0,  # All atoms
    imol_map=2
)

# Find worst 10 residues
worst_10 = sorted(correlations, key=lambda x: x[1])[:10]
for spec, corr in worst_10:
    print(f"Residue {spec}: correlation = {corr}")
```

---

## Validation Functions

### `all_molecule_ramachandran_score_py()`

**Purpose**: Comprehensive Ramachandran validation for an entire molecule.

**Use Case**: "Validate the backbone geometry" or "Find Ramachandran outliers"

**Function Signature**:
```python
PyObject *all_molecule_ramachandran_score_py(int imol)
```

**Returns**: List containing:
- Overall statistics
- Per-residue data: `[[phi, psi], residue_spec, score, [prev_res, this_res, next_res]]`

**Score Interpretation**:
- **High scores (>1.0)**: GOOD - highly favored geometry
- **Low scores (<0.01)**: BAD - outliers/unfavored regions
- Lower probability = worse geometry

**Example Usage**:
```python
rama_data = all_molecule_ramachandran_score_py(1)
residues = rama_data[5]  # Get per-residue data

# Find worst outlier (minimum score)
worst = min(residues, key=lambda x: x[2])
chain, resno = worst[1][1], worst[1][2]
score = worst[2]
print(f"Worst outlier: {chain} {resno}, score = {score}")
```

---

### `all_molecule_rotamer_score_py()`

**Purpose**: Comprehensive rotamer validation for side chains.

**Use Case**: "Check rotamer quality" or "Find unusual side-chain conformations"

**Function Signature**:
```python
PyObject *all_molecule_rotamer_score_py(int imol)
```

**Returns**: `[overall_score, n_residues]`

**Example Usage**:
```python
score, n_residues = all_molecule_rotamer_score_py(1)
print(f"Overall rotamer score: {score} for {n_residues} residues")
```

---

### `rotamer_graphs_py()`

**Purpose**: Get detailed rotamer information for each residue.

**Use Case**: "Find the worst rotamer outlier"

**Function Signature**:
```python
PyObject *rotamer_graphs_py(int imol)
```

**Returns**: List of `[chain_id, resno, ins_code, score_percentage, resname]`

**Score Interpretation**:
- **High scores (>50%)**: GOOD rotamers
- **Low scores (<5%)**: BAD rotamers - poor conformations
- 0.0 or very low: Severe outliers

**Example Usage**:
```python
rotamers = rotamer_graphs_py(1)

# Find worst rotamer (excluding missing atoms)
valid_rotamers = [r for r in rotamers if r[3] > 0]
worst = min(valid_rotamers, key=lambda x: x[3])

chain, resno, score = worst[0], worst[1], worst[3]
print(f"Worst rotamer: {chain} {resno}, score = {score}%")

# Go to worst rotamer
coot.set_go_to_atom_chain_residue_atom_name(chain, resno, 'CA')
```

---

### `deviant_geometry()`

**Purpose**: Check for unusual bond lengths, angles, and other geometric outliers.

**Function Signature**:
```python
void deviant_geometry(int imol)
```

**Returns**: None (displays results in GUI or console)

---

## Common Validation Workflow

```python
# 1. Load tutorial data
load_tutorial_model_and_data()

# 2. Run comprehensive validation
imol = 1  # Model molecule
imol_map = 2  # Map molecule

# Ramachandran validation
rama_data = all_molecule_ramachandran_score_py(imol)
worst_rama = min(rama_data[5], key=lambda x: x[2])
print(f"Worst Ramachandran: {worst_rama[1]}, score={worst_rama[2]}")

# Rotamer validation  
rotamers = rotamer_graphs_py(imol)
worst_rot = min([r for r in rotamers if r[3] > 0], key=lambda x: x[3])
print(f"Worst rotamer: {worst_rot[0]} {worst_rot[1]}, score={worst_rot[3]}%")

# Density fit validation
fit_stats = map_to_model_correlation_stats_per_residue_range_py(
    imol, "A", imol_map, 1, 0
)

# Geometry validation
deviant_geometry(imol)
```

---

## Refinement Functions

### `auto_fit_best_rotamer()`

**Purpose**: Automatically fit the best rotamer for a residue.

**Function Signature**:
```python
float auto_fit_best_rotamer(
    int imol_coords,
    const char *chain_id,
    int resno,
    const char *insertion_code,
    const char *altloc,
    int imol_map,
    int clash_flag,        # 1 to check clashes, 0 to ignore
    float lowest_probability  # Minimum acceptable probability (e.g., 0.01)
)
```

**Returns**: The new rotamer probability score

**Example Usage**:
```python
# Fix worst rotamer
new_score = auto_fit_best_rotamer(
    imol_coords=1,
    chain_id='A',
    resno=91,
    insertion_code='',
    altloc='',
    imol_map=2,
    clash_flag=1,  # Check for clashes
    lowest_probability=0.01
)
print(f"New rotamer score: {new_score}%")
```

---

### `refine_zone()`

**Purpose**: Real-space refinement of a residue range.

**Function Signature**:
```python
void refine_zone(
    int imol,
    const char *chain_id,
    int resno_start,
    int resno_end,
    const char *altconf
)
```

**Example Usage**:
```python
# Refine 5 residues around residue 42
refine_zone(1, 'A', 40, 44, '')
accept_regularizement()  # Accept the refinement
```

---

### `pepflip()`

**Purpose**: Flip a peptide bond (useful for fixing cis/trans peptides).

**Function Signature**:
```python
void pepflip(
    int imol,
    const char *chain_id,
    int resno,
    const char *ins_code,
    const char *altconf
)
```

**Example Usage**:
```python
# Flip peptide at residue 41
pepflip(1, 'A', 41, '', '')
refine_zone(1, 'A', 40, 44, '')
accept_regularizement()
```

---

## Tips for Better Documentation

### For Search Queries

When users ask questions like:
- **"Which side chain is worst fitting?"** → Use `map_to_model_correlation_stats_per_residue_range_py()`
- **"Find Ramachandran outliers"** → Use `all_molecule_ramachandran_score_py()` and look for minimum scores
- **"Check rotamer quality"** → Use `rotamer_graphs_py()` and look for minimum scores
- **"Validate the model"** → Combine Ramachandran, rotamer, density fit, and geometry checks

### Key Concepts

1. **Scores vs. Statistics**:
   - Ramachandran scores: Lower = worse (outliers have low probability)
   - Rotamer scores: Lower = worse (percentage probability)
   - Correlation: Higher = better fit to density

2. **GLY Correlation Warning**:
   - **GLY residue correlations are unreliable and should be treated with skepticism**
   - GLY has only 4 backbone atoms (N, CA, C, O) and no sidechain
   - When neighbouring residue atoms are masked out during correlation calculation, very few grid points remain
   - This leads to unreliable/meaningless correlation values for GLY
   - **Recommendation**: When identifying poorly-fitted residues, filter out GLY residues or verify GLY problems by visual inspection before attempting fixes
   
   ```python
   # Example: Filter out GLY when finding problem residues
   stats = coot.map_to_model_correlation_stats_per_residue_range_py(0, "A", 1, 1, 0)
   poor_residues = []
   for res in stats[0]:
       resno = res[0][1]
       corr = res[1][1]
       res_name = coot.residue_name(0, "A", resno, "")
       if corr < 0.7 and res_name != "GLY":  # Skip GLY
           poor_residues.append((resno, res_name, corr))
   ```

3. **Atom Mask Modes**:
   - Use mode 2 for side-chain-only analysis
   - Use mode 0 for all-atom analysis
   - Use mode 1 for main-chain analysis

4. **Workflow**:
   - Always validate BEFORE and AFTER refinement
   - Fix worst outliers first (excluding GLY correlation issues)
   - Re-validate after each fix

---

## Function Categories

### Load/Display
- `load_tutorial_model_and_data()` - Load example data
- `set_mol_displayed()` - Show/hide molecules
- `scale_zoom()` - Zoom in/out

### Navigation
- `set_go_to_atom_chain_residue_atom_name()` - Center on atom
- `rotate_x_scene()`, `rotate_y_scene()`, `rotate_z_scene()` - Rotate view

### Validation
- `all_molecule_ramachandran_score_py()` - Backbone validation
- `rotamer_graphs_py()` - Side-chain validation
- `map_to_model_correlation_stats_per_residue_range_py()` - Density fit
- `deviant_geometry()` - Geometry validation

### Refinement
- `refine_zone()` - Real-space refinement
- `auto_fit_best_rotamer()` - Fix rotamers
- `pepflip()` - Flip peptides
- `accept_regularizement()` - Accept refinement

### Model Building
- `mutate_residue_range()` - Change residue types
- `add_terminal_residue()` - Extend chains
- `delete_residue()` - Remove residues



