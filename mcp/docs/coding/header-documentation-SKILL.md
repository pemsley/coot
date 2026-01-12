---
name: coot-cpp-docs
description: "C++ documentation style conventions for Coot codebase"
---

# Coot C++ Documentation Style

## Doxygen Comment Style

When writing documentation comments for Coot C++ header files, use `//!` as the leading token, NOT `/** */` or `///`.

### ✅ CORRECT

```cpp
//! \brief Accept refined atoms into the main molecule
//!
//! This function waits for refinement to complete and commits
//! the atomic positions.
//!
//! \param imol The molecule index
//! \return PyObject* containing the results
PyObject *accept_moving_atoms_py();
```

### ❌ INCORRECT

```cpp
/**
 * @brief Accept refined atoms into the main molecule
 *
 * This function waits for refinement to complete...
 */
PyObject *accept_moving_atoms_py();
```

```cpp
/// \brief Accept refined atoms into the main molecule
///
/// This function waits for refinement to complete...
PyObject *accept_moving_atoms_py();
```

## Common Doxygen Commands

Use backslash style (`\brief`) rather than at-sign style (`@brief`):

- `\brief` - Short description
- `\param` - Parameter documentation
- `\return` - Return value documentation
- `\c` - Inline code formatting (e.g., `\c Py_False`)
- `\code{.py}` / `\endcode` - Python code blocks
- `\code{.cpp}` / `\endcode` - C++ code blocks

## Example: Full Function Documentation

```cpp
//! \brief Get the worst-fitting residues by density correlation
//!
//! Analyzes all residues in a chain and returns correlation statistics
//! for each residue against the electron density map.
//!
//! \param imol Model molecule index
//! \param chain_id Chain identifier (e.g., "A")
//! \param imol_map Map molecule index
//! \param n_per_range Residues per window (use 1 for per-residue stats)
//! \param exclude_NOC 0 to include all atoms, 1 to exclude backbone N,O,C
//!
//! \return PyObject* containing [[all_atom_stats], [sidechain_stats]]
//!
//! Example usage:
//! \code{.py}
//! stats = coot.map_to_model_correlation_stats_per_residue_range_py(0, "A", 1, 1, 0)
//! all_atom = stats[0]
//! worst = sorted(all_atom, key=lambda x: x[1][1])[:5]
//! \endcode
PyObject *map_to_model_correlation_stats_per_residue_range_py(
    int imol,
    const std::string &chain_id,
    int imol_map,
    unsigned int n_per_range,
    short int exclude_NOC);
```
