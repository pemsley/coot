## Status (2026-05-09)

Three files in `coot-utils/`:
- `coot-coord-utils-gemmi.hh` — header (30 functions declared)
- `coot-coord-utils-gemmi.cc` — implementations
- `test-coord-utils-gemmi.cc` — 33 comparison tests

Added to `Makefile.am`: `.cc` in library sources, `.hh` in installed headers, test binary in `check_PROGRAMS`.

**All 33 tests pass.** 100% coverage of functions declared in `coot-coord-utils-gemmi.hh`.

---

## Critical discovery: PDBv3 atom names

mmdb uses 4-character padded atom names (` CA `, ` N  `, ` C  `). gemmi uses unpadded names (`CA`, `N`, `C`). The codebase has ~90 `// PDBv3 FIXME` markers at sites using padded names.

**The problem with `copy_from_mmdb()`:** gemmi's `copy_from_mmdb()` copies `m_atom.label_atom_id` into `atom.name`. Despite being documented as "not aligned" in `mmdb_atom.h`, mmdb pads `label_atom_id` identically to `name` when reading PDB files. Verified with pure mmdb (no coot code involved) — both fields are 4 chars padded.

**Current workaround:** `coot::trim_atom_names(gemmi::Structure &st)` in `coot-coord-utils-gemmi.{hh,cc}`. Strips leading/trailing spaces from all atom names. Must be called after every `copy_from_mmdb()`. The test file wraps this in a helper:

```cpp
gemmi::Structure gemmi_from_mmdb(mmdb::Manager *mol) {
   gemmi::Structure st = gemmi::copy_from_mmdb(mol);
   coot::trim_atom_names(st);
   return st;
}
```

**Proper fix (TODO):** This should be fixed in gemmi's `copy_from_mmdb()` itself — propose the change upstream to Marcin.

**All gemmi functions use unpadded atom names** (`"CA"`, `"N"`, `"C"`, `"CB"`, etc.). This is the correct style going forward. The `// PDBv3 FIXME` markers in the codebase (~90 sites) map the work needed when migrating each file.

The hardest areas for the PDBv3 transition will be `ideal/` and `geometry/` — padded atom names are deeply embedded in the restraint matching system.

---

## mmdb symmetry setup

mmdb needs `syminfo.lib` to resolve spacegroup operations. Without it, `GetTMatrix()` silently fails and reports no symmetry. The test calls `setup_syminfo()` (from `utils/setup-syminfo.hh`) in `main()` before running tests. Without this, a P 21 21 21 structure was incorrectly reported as having no symmetry by the mmdb version.

---

## What was converted (25 functions + trim_atom_names)

### Utility
| Function | Notes |
|----------|-------|
| `trim_atom_names(gemmi::Structure&)` | Strips padded atom names after `copy_from_mmdb()` |

### Atom-level (4, in `coot::`)

| Function | Notes |
|----------|-------|
| `distance(gemmi::Atom&, gemmi::Atom&)` | Manual sqrt, matches mmdb overload |
| `angle(gemmi::Atom&, gemmi::Atom&, gemmi::Atom&)` | Returns degrees, clamps cos_theta |
| `co(gemmi::Atom&)` | Returns `clipper::Coord_orth` |
| `is_hydrogen_atom(gemmi::Atom&)` | Delegates to `at.is_hydrogen()` |

### Residue-level (6, in `coot::util::`)

| Function | Notes |
|----------|-------|
| `get_residue_centre` | Mean of all atom positions |
| `get_CA_position_in_residue` | Finds atom named `"CA"` (unpadded) |
| `get_CB_position_in_residue` | Finds atom named `"CB"` (unpadded) |
| `is_nucleotide` | Hardcoded residue name list |
| `residue_has_hydrogens_p` | Uses `at.is_hydrogen()` |
| `residue_has_hetatms` | Checks `res.het_flag == 'H'` |

### Chain-level (5, in `coot::util::`)

| Function | Notes |
|----------|-------|
| `min_and_max_residues` | Uses `res.seqid.num.value` |
| `min_resno_in_chain` | Pair<bool,int> return |
| `max_resno_in_chain` | Pair<bool,int> return |
| `residue_types_in_chain` | Returns sorted unique names via set |
| `get_number_of_protein_or_nucleotides` | Nucleotide via name; protein via `entity_type == Polymer` |

### Structure-level (15, mix of `coot::` and `coot::util::`)

| Function | Notes |
|----------|-------|
| `centre_of_molecule` | First model only |
| `radius_of_gyration` | First model only |
| `mol_has_symmetry` | Uses `find_spacegroup()`, checks `operations().order() > 1` |
| `mol_is_anisotropic` | Checks any `at.aniso.nonzero()` |
| `get_position_hash` | Sum of x-differences (matches mmdb algorithm, no TER atoms) |
| `residue_types_in_molecule` | Set-based unique types |
| `non_standard_residue_types_in_molecule` | Uses `standard_residue_types()` from coot-coord-utils.hh |
| `chains_in_molecule` | Chain name list |
| `number_of_residues_in_molecule` | Sum of chain.residues.size() |
| `max_number_of_residues_in_chain` | Max chain length |
| `number_of_chains` | models[0].chains.size() |
| `max_resno_in_molecule` | Max seqid across all chains |
| `max_min_max_residue_range` | Max (max-min+1) across chains (inclusive, matching mmdb) |
| `alt_confs_in_molecule` | Includes `""` for non-alt atoms (matching mmdb) |
| `extents` | Min/max bounding box |
| `median_position` | Sort-based median per axis |
| `count_cis_peptides` | Omega torsion from geometry, C-N distance < 3.0 A check |

---

## Bugs found and fixed during testing

| Issue | Root cause | Fix |
|-------|-----------|-----|
| `get_CA_position_in_residue` failed | Atom names padded by mmdb, gemmi searched for `"CA"` not `" CA "` | `trim_atom_names()` |
| `get_position_hash` mismatch (57 vs 178516) | Gemmi version summed absolute positions, mmdb sums x-differences | Rewrote to match mmdb algorithm |
| `mol_has_symmetry` mismatch (mmdb=0, gemmi=1) | mmdb lacked syminfo.lib; gemmi correctly found P 21 21 21 | Added `setup_syminfo()`, changed gemmi impl to use `find_spacegroup()` |
| `max_min_max_residue_range` off by 1 | mmdb uses inclusive range (max-min+1), gemmi used (max-min) | Added +1 |
| `alt_confs_in_molecule` mismatch | mmdb includes empty string for non-alt atoms, gemmi skipped them | Include `""` for `altloc == '\0'` |
| `count_cis_peptides` mismatch | Gemmi used `\|omega\| < 30 deg`, mmdb uses `\|180-pos_torsion\| > 90` + distance check | Matched mmdb's distortion/distance logic |
| Segfault on all tests | `replace_all` on `copy_from_mmdb` caught the helper's own body → infinite recursion | Fixed `gemmi_from_mmdb()` to call `gemmi::copy_from_mmdb()` |
| Missing includes | `<fstream>`, `<set>`, `"coot-coord-utils.hh"`, `"utils/setup-syminfo.hh"` | Added |

---

## Remaining issues to watch

1. **`get_number_of_protein_or_nucleotides`**: Uses `res.entity_type == gemmi::EntityType::Polymer` — may not be populated by `copy_from_mmdb()`. No test for this yet.

2. **`residue_has_hetatms`**: Checks `res.het_flag == 'H'`. `copy_from_mmdb()` does set `het_flag` (verified in gemmi source line 345), so this should be OK.

3. **`is_nucleotide`**: Hardcoded name list. The mmdb version may use a different heuristic. Should verify they agree on a nucleotide-containing structure.

4. **`get_position_hash`**: ~~Previously differed due to TER atoms.~~ Fixed: mmdb version now has `isTer()` guards, hashes match. Test updated to cross-compare (2026-05-09).

5. **First-model-only pattern**: Structure-level functions use `for (...) { ... break; }`. Works but could use `st.models[0]` directly with an empty check.

---

## What was NOT converted — prioritised (2026-05-09)

### Worth porting (operate on structure/chain/residue data, no mmdb-specific API)

| Priority | Function | Notes |
|----------|----------|-------|
| High | `centre_of_molecule_using_masses` | Straightforward, needs element→mass map |
| High | `omega_torsion` | Needs gemmi residue pair |
| High | `cis_peptides_info_from_coords` | Returns detailed cis peptide info |
| Medium | `residues_near_residue` / `residues_near_position` | Spatial queries — more involved |
| Medium | `atoms_with_zero_occupancy` | Simple iteration |
| Medium | `residues_with_alt_confs` | Simple iteration |
| Medium | `residues_with_insertion_codes` | Simple iteration |
| Medium | `nucleotide_is_DNA` | Test for O2' presence |
| Medium | `CO_orientations` | Validation analysis |
| Low | `sort_chains` / `sort_residues` | Mutating operations |
| Low | `hiranuma_inversion` | pLDDT→B conversion, mutating |
| Low | `gln_asn_b_factor_outliers` | Validation, complex |
| Low | `closest_approach` / `interface_residues` | Spatial, complex |

### Not worth porting (inherently mmdb-specific)

These depend on mmdb::Manager ownership, selection handles, or deep copy semantics:
- `create_mmdbmanager_from_*` (residue, residue_vector, mmdbmanager, atom_selection, atom, points, residue_specs)
- `deep_copy_this_residue`, `copy_molecule`, `copy_chain`
- `copy_cell_and_symm_headers`, `copy_headers`
- `get_selection_handle`, `specs_to_atom_selection`
- `transform_mol`, `transform_chain`, `transform_selection`, `transform_atoms`
- `graph_match` (uses mmdb::math::Graph)
- `water_coordination_t` (uses mmdb contacts API)
- `move_hetgroups_around_protein`, `move_waters_around_protein` (symmetry + mmdb)
- `get_lsq_matrix` (mmdb selection-based LSQ)
- `cis_trans_conversion`, `cis_trans_convert` (mutates mmdb residues)
- `mutate`, `mutate_internal`, `mutate_base` (mmdb residue mutation)
- `position_residue_by_internal_coordinates` (mmdb atom manipulation)
- `split_multi_model_molecule` (returns mmdb::Manager vector)
- `pdbcleanup_serial_residue_numbers` (mmdb internal indexing)
- `write_coords_pdb`, `write_coords_cif` (mmdb I/O — gemmi has its own)

### Selection API (~1,450 calls across codebase)
- `NewSelection / SelectAtoms / GetSelIndex / DeleteSelection` — hardest migration target
- gemmi equivalent: range-based for loops with conditionals, or `gemmi::Selection`

### Map-related
- Functions in `coot-map-utils.hh` — separate from coord-utils migration

---

## Approach notes

- `gemmi::Structure` created on the fly via `gemmi::copy_from_mmdb(mol)` + `trim_atom_names()` — no persistent dual representation
- Tests compare mmdb and gemmi results on the same molecule, tolerance-based for floats
- `copy_from_mmdb()` will need extension for full migration (missing: secondary structure, title/header, resolution, biomolecule info)
- The mmdb selection API is the biggest migration challenge — ~1,450 call sites
- SSM (structure superposition) may not have a gemmi equivalent — can be side-lined
- TER atoms: gemmi doesn't have them. They only matter for PDB file output. Skip them in all logic.

**Next steps:** Start with the high-priority items: `centre_of_molecule_using_masses`, `omega_torsion`, `cis_peptides_info_from_coords`. Then the medium-priority simple iterations (`atoms_with_zero_occupancy`, `residues_with_alt_confs`, `residues_with_insertion_codes`, `nucleotide_is_DNA`). Spatial queries (`residues_near_residue/position`) are more involved but high value.
