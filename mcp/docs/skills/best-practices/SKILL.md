---
name: coot-best-practices
description: "Best Practices for using Coot MCP"
---
# Coot Python API Best Practices

## Overview

This skill provides best practices for interacting with Coot's Python API through the MCP server. Following these guidelines ensures optimal performance, correct usage, and reliable results.

## Startup Procedure

**CRITICAL: On every "Coot Mode" start, before doing anything else:**

1. Read the `/mnt/skills/user/coot-essential-api/SKILL.md` file
2. Extract all function names mentioned in that file
3. Call `get_function_descriptions()` with the complete list of function names
4. This loads the essential API documentation into context, providing immediate access to the ~25 core functions needed for typical validation and model-building workflows

This startup procedure eliminates the need to search for basic functions during the session and ensures you have the foundational API ready to use.

Example:
```python
# After reading coot-essential-api/SKILL.md, call:
Coot:get_function_descriptions([
    "set_refinement_immediate_replacement",
    "is_valid_model_molecule",
    "is_valid_map_molecule",
    # ... all other functions from the essential API
])
```

Only after completing this startup should you proceed with the user's task.


## Critical Rule: Prefer C++ Functions Over Python Wrappers

**ALWAYS use `coot.*_py()` functions directly instead of `coot_utils.*` equivalents when they are simple passthroughs.**

### Why?

1. **Performance**: C++ functions are significantly faster (no Python overhead)
2. **Import requirements**: `coot` is auto-imported, `coot_utils` requires explicit import
3. **Simplicity**: Direct access to the core API without unnecessary abstraction layers
4. **Reliability**: Fewer layers means fewer potential points of failure

## Module Import Requirements

### Auto-imported
- **`coot`** - The core C++/SWIG binding is automatically available
- No import statement needed

### Requires explicit import
- **`coot_utils`** - Python utility library built on top of `coot`
- Must execute: `import coot_utils` before using any of its functions

**Example of the problem:**
```python
# This will fail with NameError if coot_utils not imported
coot_utils.chain_ids(0)

# Solution:
import coot_utils
coot_utils.chain_ids(0)  # Now works
```

## Function Naming Conventions

### C++ Functions (SWIG bindings)
- End with `_py()` suffix
- Examples: `closest_atom_simple_py()`, `is_valid_model_molecule()`, `chain_id_py()`
- These are the **core functions** - prefer these

### Python Wrapper Functions
- No `_py()` suffix
- Examples: `closest_atom()`, `closest_atom_simple()`, `chain_ids()`
- Only use when they provide **genuine convenience**

## Specific Function Guidance

### ✅ CORRECT: Getting Closest Atom

```python
# Get closest atom across all displayed molecules
atom_spec = coot.closest_atom_simple_py()
# Returns: [imol, chain-id, resno, ins-code, atom-name, alt-conf, [x, y, z]]

# Get closest atom in specific molecule
atom_spec = coot.closest_atom_py(0)
# Returns: [imol, chain-id, resno, ins-code, atom-name, alt-conf, [x, y, z]]

# Get raw closest atom (no CA substitution)
atom_spec = coot.closest_atom_raw_py()
```

### ❌ INCORRECT: Don't use coot_utils for simple passthroughs

```python
# DON'T DO THIS - unnecessary import and no added value
import coot_utils
atom_spec = coot_utils.closest_atom(0)  # Just calls coot.closest_atom_py()

# DON'T DO THIS EITHER
atom_spec = coot_utils.closest_atom_simple()  # Just calls coot.closest_atom_simple_py()
```

### ✅ CORRECT: When to use coot_utils

Use `coot_utils` functions when they provide **genuine convenience or abstraction**:

```python
import coot_utils

# chain_ids() is a convenience wrapper that constructs a list
# It calls coot.chain_id_py() in a loop and builds a list
chains = coot_utils.chain_ids(0)  # Returns: ['A', 'B']

# Without coot_utils, you'd have to do:
n = coot.n_chains(0)
chains = [coot.chain_id_py(0, i) for i in range(n)]
```

### Checking Molecule Validity

```python
# ✅ CORRECT: Direct C++ function
if coot.is_valid_model_molecule(0):
    print("Molecule 0 is a valid model")

if coot.is_valid_map_molecule(1):
    print("Molecule 1 is a valid map")

# ❌ INCORRECT: Don't use coot_utils for this
import coot_utils
if coot_utils.valid_model_molecule_qm(0):  # Unnecessary
    pass
```

### Getting Molecule Information

```python
# ✅ CORRECT: Direct access
name = coot.molecule_name(0)
n_chains = coot.n_chains(0)

# ✅ CORRECT: When coot_utils adds value
import coot_utils
chains = coot_utils.chain_ids(0)  # Convenience wrapper
```
## MMDB Atom Selection Syntax

When using functions like `new_molecule_by_atom_selection()` or `superpose_with_atom_selection()`, use MMDB atom selection strings to specify which atoms to include.

### Format
```
//chn/seq(res).ic/atm[elm]:aloc
```

### Components

- **`//`** - Model specifier (typically `//` for single-model structures, or `/1/` for model 1)
- **`chn`** - Chain ID (e.g., `A`, `B`, `X`)
  - Multiple chains can be specified with commas: `A,B,C`
- **`seq`** - Residue number or range:
  - Single: `50`
  - Range: `10-20`
- **`res`** - Residue name in parentheses (e.g., `(HIS)`, `(ALA)`, `(GLY)`)
- **`ic`** - Insertion code
- **`atm`** - Atom name (e.g., `CA`, `N`, `O`)
- **`elm`** - Element in square brackets (e.g., `[C]`, `[N]`)
- **`aloc`** - Alternate location indicator

All components are optional - you only need to specify what you want to filter.

### Examples
```python
# Select entire chain
"//A"                          # All atoms in chain A
"//A,B,C"                      # All atoms in chains A, B, and C

# Select residue range
"//A/12-130"                   # Residues 12-130 in chain A
"//A/12-130/CA"                # CA atoms from residues 12-130 in chain A

# Select specific residue type
"//B/10-20(GLY)"               # GLY residues 10-20 in chain B
"//A/*(HIS)"                   # All HIS residues in chain A

# Select specific atom
"//A/50/CA"                    # CA atom of residue 50 in chain A
"//A/50(HIS)/CA"               # CA atom of HIS 50 in chain A

# Select multiple chains in one selection
"//X,Y,1,2"                    # All atoms in chains X, Y, 1, and 2
```

### Usage Examples
```python
import coot

# Create a new molecule with chains A and B
imol_ab = coot.new_molecule_by_atom_selection(0, "//A,B")

# Create a new molecule with CA atoms from residues 10-50 in chain A
imol_ca = coot.new_molecule_by_atom_selection(0, "//A/10-50/CA")

# Create molecule with transcription factor (chains X, Y) and DNA (chains 1, 2)
imol_complex = coot.new_molecule_by_atom_selection(0, "//X,Y,1,2")

# Superpose using atom selection
coot.superpose_with_atom_selection(
    imol1=0,
    imol2=1,
    mmdb_atom_sel_str_1="//A/10-100",
    mmdb_atom_sel_str_2="//A/10-100",
    move_imol2_copy_flag=0
)
```

## Code Execution Patterns

### Single-line expressions
Single-line expressions return their evaluated value:

```python
coot.is_valid_model_molecule(0)  # Returns: 1 or 0
```

### Multi-line code blocks
Multi-line blocks require explicit return or final expression:

```python
# ❌ This returns None (print doesn't return a value)
mols = []
for i in range(3):
    mols.append(i)
print(mols)

# ✅ CORRECT: Return the value or use final expression
mols = []
for i in range(3):
    mols.append(i)
mols  # Final expression is returned

# ✅ ALSO CORRECT: List comprehension (single expression)
[i for i in range(3)]
```

## MMDB Atom Selection Syntax

When using functions like `new_molecule_by_atom_selection()` or `superpose_with_atom_selection()`, use MMDB atom selection strings to specify which atoms to include.

### Format
```
//chn/seq(res).ic/atm[elm]:aloc
```

### Components

- **`//`** - Model specifier (typically `//` for single-model structures, or `/1/` for model 1)
- **`chn`** - Chain ID (e.g., `A`, `B`, `X`)
- **`seq`** - Residue number or range:
  - Single: `50`
  - Range: `10-20`
- **`res`** - Residue name in parentheses (e.g., `(HIS)`, `(ALA)`, `(GLY)`)
- **`ic`** - Insertion code
- **`atm`** - Atom name (e.g., `CA`, `N`, `O`)
- **`elm`** - Element in square brackets (e.g., `[C]`, `[N]`)
- **`aloc`** - Alternate location indicator

All components are optional - you only need to specify what you want to filter.

### Examples
```python
# Select entire chain
"//A"                          # All atoms in chain A

# Select residue range
"//A/12-130"                   # Residues 12-130 in chain A
"//A/12-130/CA"                # CA atoms from residues 12-130 in chain A

# Select specific residue type
"//B/10-20(GLY)"               # GLY residues 10-20 in chain B
"//A/(HIS)"                    # All HIS residues in chain A

# Select specific atom
"//A/50/CA"                    # CA atom of residue 50 in chain A
"//A/50(HIS)/CA"               # CA atom of HIS 50 in chain A

# Multiple chains (create separate selections and merge)
imol_a = coot.new_molecule_by_atom_selection(imol, "//A")
imol_b = coot.new_molecule_by_atom_selection(imol, "//B")
```

### Usage Example
```python
import coot

# Create a new molecule with only chain A
imol_chain_a = coot.new_molecule_by_atom_selection(0, "//A")

# Create a new molecule with CA atoms from residues 10-50
imol_ca = coot.new_molecule_by_atom_selection(0, "//A/10-50/CA")

# Superpose using atom selection
coot.superpose_with_atom_selection(
    imol1=0,
    imol2=1,
    mmdb_atom_sel_str_1="//A/10-100/CA",
    mmdb_atom_sel_str_2="//A/10-100/CA",
    move_imol2_copy_flag=0
)
```


## Common Tasks Reference

### Loading Tutorial Data

```python
# ✅ CORRECT function name
coot.load_tutorial_model_and_data()

# ❌ INCORRECT function names that don't exist
# coot.tutorial_model_and_data()  # Wrong!
```

### Listing Molecules

```python
# Check molecules 0-5
[(i, coot.is_valid_model_molecule(i), coot.is_valid_map_molecule(i)) for i in range(6)]

# Get molecule names for valid molecules
for i in range(10):
    if coot.is_valid_model_molecule(i):
        print(f"Model {i}: {coot.molecule_name(i)}")
    elif coot.is_valid_map_molecule(i):
        print(f"Map {i}: {coot.molecule_name(i)}")
```

### Working with Chain IDs

```python
# ✅ CORRECT: Use coot_utils for convenience
import coot_utils
chains = coot_utils.chain_ids(0)  # Returns: ['A', 'B', 'C']

# Iterate over chains
for chain in chains:
    print(f"Chain {chain}")
```

### Getting Active/Closest Residue

```python
# ✅ Get closest atom across displayed molecules
atom = coot.closest_atom_simple_py()
if atom:
    imol, chain, resno, ins, atom_name, alt, coords = atom[0], atom[1], atom[2], atom[3], atom[4], atom[5], atom[6]
    
# ✅ Get active residue (with potential CA substitution)
import coot_utils
active = coot_utils.active_residue()
if active:
    imol, chain, resno, ins, atom_name, alt = active
```

## Density Fit Analysis

### Map Correlation Functions

```python
# Get correlation for specific residues
import coot_utils

# Single residue
residue_spec = ["A", 42, ""]
correlation = coot.density_score_residue_py(0, residue_spec, 1)

# Per-residue correlation for a range
residue_specs = [["A", i, ""] for i in range(40, 50)]
results = coot.map_to_model_correlation_per_residue_py(0, residue_specs, 0, 1)
# Returns: [(residue_spec, correlation), ...]

# Main function for "which residue fits worst?"
stats = coot.map_to_model_correlation_stats_per_residue_range_py(
    0, "A", 1, 100, 1  # imol, chain, start, end, imol_map
)
```

## Zoom and View Settings

When adjusting the view in Coot, remember that **higher zoom values mean the molecule appears larger on screen** (i.e., zoomed in), while **lower values show more of the scene** (zoomed out). Typical ranges:

- **150-300**: Whole-molecule overview (appropriate for ribbons, surfaces, overall architecture)
- **50-100**: Domain or region level
- **20-50**: Residue-level detail (inspecting side chains, density fit, rotamers)

Small proteins like RNase A (~124 residues) may appear compact even at zoom 200, while larger complexes will fill the screen at lower zoom values. When presenting a ribbon diagram or other overview representation, consider turning off the bond representation (`coot.set_mol_displayed(imol, 0)`) and hiding electron density maps (`coot.set_map_displayed(imol_map, 0)`) to reduce visual clutter while keeping the ribbon mesh visible.

Use `coot.zoom_factor()` to query the current zoom level and `coot.set_zoom(value)` to set it. For interactive exploration, users can also adjust zoom with the scroll wheel.


## Performance Considerations

### Function Call Overhead

```python
# ❌ SLOW: Multiple function calls through Python wrapper
import coot_utils
for i in range(1000):
    atom = coot_utils.closest_atom(0)  # Unnecessary indirection

# ✅ FAST: Direct C++ calls
for i in range(1000):
    atom = coot.closest_atom_py(0)
```

### When Python Wrappers Are Worth It

Python wrappers are valuable when they:
1. **Aggregate multiple C++ calls** (e.g., `chain_ids()` calls `chain_id_py()` in a loop)
2. **Transform data** into more convenient formats
3. **Provide meaningful abstractions** that simplify complex operations
4. **Add error handling** or validation logic

## Decision Tree

```
Need to call a Coot function?
│
├─ Does it require coot_utils for convenience features?
│  └─ YES → import coot_utils and use it
│     Examples: chain_ids(), active_residue()
│
└─ NO → Use coot.*_py() directly
   Examples: closest_atom_simple_py(), is_valid_model_molecule()
```

## Quick Reference Table

| Task | ❌ Avoid | ✅ Use Instead | Reason |
|------|---------|---------------|--------|
| Get closest atom (all molecules) | `coot_utils.closest_atom_simple()` | `coot.closest_atom_simple_py()` | Direct C++, no import needed |
| Get closest atom (specific mol) | `coot_utils.closest_atom(imol)` | `coot.closest_atom_py(imol)` | Direct C++, no import needed |
| Get chain IDs | Multiple C++ calls | `coot_utils.chain_ids(imol)` | Convenience wrapper adds value |
| Check if valid model | `coot_utils.valid_model_molecule_qm()` | `coot.is_valid_model_molecule(imol)` | Direct C++, clearer name |
| Get molecule name | N/A | `coot.molecule_name(imol)` | Direct C++ only |
| Load tutorial data | `coot.tutorial_model_and_data()` | `coot.load_tutorial_model_and_data()` | Correct function name |

## Common Mistakes to Avoid

### 1. Using coot_utils without import
```python
# ❌ Will fail with NameError
chains = coot_utils.chain_ids(0)

# ✅ Import first
import coot_utils
chains = coot_utils.chain_ids(0)
```

### 2. Using coot_utils when unnecessary
```python
# ❌ Unnecessary indirection
import coot_utils
atom = coot_utils.closest_atom_simple()

# ✅ Direct and faster
atom = coot.closest_atom_simple_py()
```

### 3. Wrong function names
```python
# ❌ Function doesn't exist
coot.tutorial_model_and_data()

# ✅ Correct name
coot.load_tutorial_model_and_data()
```


### 4. Getting output from multi-line code
```python
# ✅ Use print() to see output - it appears in stdout
result = []
for i in range(5):
    result.append(i)
print(result)  # Output: [0, 1, 2, 3, 4]

# ❌ A bare expression at the end of multi-line code does NOT return a value
result = []
for i in range(5):
    result.append(i)
result  # Returns None - this doesn't work!

## API Discovery Tools

### Using search_coot_functions

The `search_coot_functions` tool is your primary method for finding Coot functions. It supports powerful search patterns:

**Space-separated words = Logical AND**
```python
# Find functions containing ALL these words
search_coot_functions("map model correlation")
# Returns functions like: map_to_model_correlation_stats_per_residue_range_py

search_coot_functions("residue range chain")
# Returns functions dealing with residue ranges in chains

search_coot_functions("min max residue")
# Returns functions with all three words (not just any one)
```

**Single words = Simple search**
```python
search_coot_functions("correlation")  # All functions with "correlation"
search_coot_functions("validation")   # All functions with "validation"
search_coot_functions("rotamer")      # All functions with "rotamer"
```

**Common search patterns:**
- `"map correlation"` - density fit functions
- `"residue validation"` - geometry checking
- `"chain residue"` - chain/residue operations
- `"ligand environment"` - ligand analysis
- `"ramachandran"` - backbone validation
- `"density fit"` - map fitting functions

### ✅ CORRECT Search Strategy

```python
# Looking for functions to get residues in a chain
search_coot_functions("chain residue")  # Logical AND

# Looking for min/max residue number functions
search_coot_functions("min max residue")  # All three words required

# Looking for correlation analysis
search_coot_functions("correlation residue")  # Both words required
```

### ❌ INCORRECT Search Strategy

```python
# DON'T use grep with pipe (|) when you mean AND
# This searches for min OR max OR residue (logical OR)
# Use search_coot_functions with spaces instead
```

### When to use each discovery tool

1. **search_coot_functions(pattern)** - First choice
   - Use space-separated words for AND logic
   - Returns max 40 results with documentation
   - Best for targeted searches

2. **list_coot_categories()** - For browsing
   - Returns: ['load', 'read', 'display', 'refinement', 'validation', 'ligand', 'util']
   - Use when you want to explore a general area

3. **get_functions_in_category(category)** - For comprehensive lists
   - Returns all functions in a category (50-200 functions)
   - Use after identifying the right category

### Search Tips

- Start with **2-3 specific words** that describe what you need
- If too many results, add more words to narrow down
- If no results, try synonyms or broader terms
- Common terms: validation, correlation, residue, chain, map, model, ligand, fit, geometry

## Summary

1. **Always prefer `coot.*_py()` functions** when they're simple passthroughs
2. **Only use `coot_utils` functions** when they add genuine convenience
3. **Remember `coot` is auto-imported**, `coot_utils` is not
4. **Use single-line expressions** when possible for cleaner returns
5. **Check function names** - `load_tutorial_model_and_data()` not `tutorial_model_and_data()`
6. **Use `search_coot_functions` with space-separated words** for AND logic when searching the API

Following these practices ensures optimal performance and correct API usage when working with Coot through the MCP server.
