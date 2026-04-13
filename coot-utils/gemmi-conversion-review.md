## Status (2026-04-10)

Three files in `coot-utils/`:
- `coot-coord-utils-gemmi.hh` — header
- `coot-coord-utils-gemmi.cc` — implementations
- `test-coord-utils-gemmi.cc` — 26 comparison tests

Added to `Makefile.am`: `.cc` in library sources, `.hh` in installed headers, test binary in `check_PROGRAMS`.

**All 26 tests pass.** Build playground: `build-for-claude/coot-utils/`

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

4. **`get_position_hash`**: mmdb and gemmi hashes differ numerically (TER atoms affect mmdb's x-difference chain). Test verifies gemmi self-consistency rather than cross-comparing. This is fine — the hash just needs to detect when atoms move.

5. **First-model-only pattern**: Structure-level functions use `for (...) { ... break; }`. Works but could use `st.models[0]` directly with an empty check.

---

## What was NOT converted (~75+ functions remaining in coot-coord-utils.hh)

### Residue operations (high value for migration)
- `get_residue(chain_id, resno, ins_code, mol)` — residue lookup
- `get_atom(atom_spec, mol)` — atom lookup by spec
- `deep_copy_residue` — residue duplication
- `get_residue_name(chain_id, resno, ins_code, mol)` — name lookup
- `residues_near_residue` — spatial queries
- `residues_near_position` — spatial queries
- `get_waters` — water identification

### Atom selection / searching
- `atoms_with_zero_occ` — find zero-occupancy atoms
- `fill_residue_alt_confs` — enumerate alt confs per residue
- `contact_info` — atom contacts

### Whole-molecule operations
- `sort_chains`, `sort_residues` — reordering
- `shift_to_origin` — coordinate transformation
- `copy_and_delete_hydrogens` — hydrogen stripping
- `move_waters_to_last_chain`

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

**Next steps:** Prioritize converting residue-lookup functions (`get_residue`, `get_atom`) and spatial queries (`residues_near_residue/position`) — most commonly used across the codebase.
