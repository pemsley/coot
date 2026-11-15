---
name: converting-coot-scripts-to-chapi
description: Use when converting or porting a desktop Coot Python script (one that does `import coot` and calls coot.read_pdb, coot.*_py functions, superpose, refine_zone, accept_moving_atoms, set_rotation_centre, etc.) into a script for the Coot headless API (chapi, `coot_headless_api`, `molecules_container_t`), or when a Coot GUI script must run without graphics.
---

# Converting Coot Scripts to chapi

## Overview

Desktop Coot scripts call module-level functions that act on global state and show
their results in the GUI. chapi scripts call methods on one `molecules_container_t`
object, have no graphics, and **return** their results (often as JSON strings).
The source of truth for chapi is `api/molecules-container.hh`; a function is callable
from Python only if it is bound in `api/molecules-container-nanobind.cc`.

## Quick Reference

| Desktop Coot | chapi |
|---|---
| `import coot` | `import coot_headless_api`<br>`chapi = coot_headless_api.molecules_container_t(False)` (arg is `verbose`) |
| `coot.f(...)` | `chapi.f(...)` - same name for many basics (`read_pdb`, `get_types_in_molecule`, `set_imol_refinement_map`) |
| `read_cif_dictionary(fn)` | `chapchapi.import_cif_dictionary(fn, chapi.get_imol_enc_any())` |
| `read_ccp4_map(fn, 0)` | `chapi.read_ccp4_map(fn, False)` (bool `is_a_difference_map`) |
| `write_pdb_file(imol, fn)` | `chapi.write_coordinates(imol, fn)` |
| `refine_zone(...)` + `accept_moving_atoms()` | `chapi.refine_residue_range(imol, chain, start, end, 500)` or `chapi.refine_residues_using_atom_cid(imol, cid, mode, 500)` - no accept step; use `n_cycles` 500 |
| `superpose(imol_ref, imol_mov, 0)` (whole molecule) | `chapi.SSM_superpose(imol_ref, chain_ref, imol_mov, chain_mov)` is per chain: take the chain from the script's context, or from `chapi.get_chains_in_model(imol)`, and say which you chose. For known ranges use `add_lsq_superpose_match(...)` + `lsq_superpose(ref, mov)` |
| `map_to_model_correlation_per_residue_py(imol, specs, 0, imol_map)` | `mc.density_correlation_analysis(imol, imol_map)`: same per-residue correlation (atom mask mode 0, all atoms), but for every residue. Pick out the wanted specs from `.cviv`/`.rviv` in Python. |
| view/graphics calls (`set_rotation_centre`, `set_go_to_atom_*`, `set_zoom`, colours, display toggles) | drop them |
| `something_py(python_list, ...)` that displays results | usually a chapi function taking a **JSON string** and returning a JSON string: build a dict, `json.dumps()`, `json.loads()` the result |

`molecules_container_t()` already calls `geometry_init_standard()` - don't call it again.

**Return values** are one of three kinds:
- plain values (`int` imol, `bool`, lists of strings)
- JSON strings - `json.loads()` them
- bound C++ objects. Example: `validation_information_t` → `.cviv` (chains) → `.rviv` (residues) → `.residue_spec.chain_id/.res_no/.ins_code` and `.function_value`. Look up the field names in the `nb::class_<...>` blocks of the nanobind file.

## Finding the Equivalent

1. `grep -n "<keyword>" api/molecules-container.hh` - names often match, but sometimes differ (e.g. Coot function names might have a trailing "_py"); if a match is not found then search by concept, not by the Coot name.
2. Confirm the binding: `grep -n '"<name>"' api/molecules-container-nanobind.cc`. Not bound → not callable from Python, even if it is in the header.
3. Read the doc comment for the input format (JSON keys, CID syntax, enum strings such as `"TRIPLE"`).
4. **Check whether the chapi function does more than the Coot one.** A chapi function may already do a step that the Coot script did by hand. Delete that step rather than translating it.
5. If there is no equivalent, say so in a comment in the script and tell the user. Don't substitute a different metric without flagging it.

## Example (desktop `test-cfc.py` → `test-chapi-cfc.py`)

```python
# Desktop: superpose by hand, then the GUI displays the clusters
coot.read_cif_dictionary(cif_file)
for mi in mol_infos:
    coot.superpose(imol_0, mi[0], 0)
coot.new_chemical_feature_clustering_py(mol_infos)

# chapi: the function superposes internally (if told the molecules are not
# superposed) and returns JSON; nothing is displayed or moved
chapi.import_cif_dictionary(cif_file, chapi.get_imol_enc_any())
mol_info = {"molecules": mol_infos,               # [[imol, "RES"], ...]
            "molecules_are_superposed": False,
            "superposition_radius": 10.0}
results = json.loads(chapi.chemical_feature_clustering(json.dumps(mol_info)))
```

## Common Mistakes

- **Translating a step that the chapi function already does.** In the example above, the superposition loop goes away. Check the doc comment first.
- **Guessing arguments that desktop Coot defaulted.** Examples are chain IDs for `SSM_superpose` and `n_cycles`. Find a value used in `api/test-molecules-container.cc` or `api/test_chapi.py`, or ask.
- **Leaving GUI-era file checks in place.** Example: requiring `BaseFeatures.fdef` in the current directory, which chapi now finds via `$RDKIT_BASEFEATURES` and `$prefix/share/RDKit/Data`.
- **Expecting results on screen.** chapi returns them; print or `json.loads` them.
- **Assuming in-place moves.** Some chapi functions move the model (`SSM_superpose`, `lsq_superpose`); others report positions in another frame without moving anything. Read the docs.

## Known Gaps (no chapi equivalent yet)

| Desktop Coot | Nearest chapi | Difference |
|---|---|---|
| `map_to_model_correlation_py(imol, specs, ...)` (one value for a set of residues) | `density_correlation_analysis(imol, imol_map)` | per residue, not one combined value |
| `map_to_model_correlation_per_residue_py(...)` with atom mask mode other than 0 (main-chain, side-chain, ...) | `density_correlation_analysis(imol, imol_map)` | mask mode 0 (all atoms) only |
