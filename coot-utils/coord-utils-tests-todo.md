# coot-coord-utils.cc — functions needing tests

**18** `coot::util` functions defined in `coot-utils/coot-coord-utils.cc` that have
**no coverage** in either `test-coord-utils.cc` (C++ unit tests) or `python-tests/`
(integration tests via the `coot.*` scripting API). (Was 19; `get_residue_alt_confs`
removed — it's `#if 0`'d here and actually lives in geometry/mol-utils.cc. See row 16.)

Of the 18, **3 are BLOCKED** on gemmi gaining graph matching (#1–#3), so the
**testable-now set is 15** (#4–#15, #17–#19).

## CALLER-SCOPE PRIORITISATION (2026-06-20)
Only functions reached from the libcootapi/headless build (api/coot-utils/ideal/ligand/...) are in
the current migration scope. Functions called ONLY from `src/` (the GUI layer) are NOT compiled into
coot_headless_api/libcootapi, so defer them (like #15 chain_atoms_segid, now being removed). Functions
with NO external callers are internal/comparators or near-dead — investigate before any work.
(`current-state/` is a duplicate snapshot tree — ignored in this analysis.)

- **In scope — DO THESE:** #5 close_residues (api, ideal), #12 chiral_4th_atom (api),
  #17 shift_to_origin (coot-utils, ligand), #18 translate_close_to_origin (coot-utils),
  #19 position_residue_by_internal_coordinates (coot-utils). [#13 nucleotide_is_DNA already DONE.]
- **src-only — DEFER:** #7 mol_by_symmetry, #9 set_mol_cell, #10 rotate_atom_about,
  #11 standardize_peptide_C_N_distances.
- **No external callers (investigate / maybe dead):** #4 compare_residues (comparator),
  #6 residues_near_residues_for_residues, #8 symmetry_move_atoms, #14 chains_in_atom_selection.

LESSON (user, 2026-06-20): check caller scope BEFORE picking a function to migrate/test — don't
prioritise src-only functions. chain_atoms_segid (#15) was such a trap (1 caller, src NCS code).

Excluded from this list: functions that shouldn't be unit-tested directly
(comparators, ctors), molecule copy/link plumbing that's indirectly exercised, and
the four already covered transitively by python-tests (`cis_trans_conversion`,
`cis_trans_convert`, `move_waters_around_protein`, `nucleotide_to_nucleotide`).

These are all mmdb-typed and on the path to gemmi migration ([[project_mmdb_to_gemmi]]),
so prefer writing them as **mmdb-vs-gemmi differential tests** (run both backends on the
same structure, compare outputs) rather than hand-rolled expected values — see
`test-coord-utils-gemmi.cc` for the pattern.

Header = coot-coord-utils.hh, def = coot-coord-utils.cc.

| # | function | header decl | .cc def | notes |
|---|----------|-------------|---------|-------|
| 1 | `graph_match` | `graph_match_info_t graph_match(mmdb::Residue *res_moving, mmdb::Residue *res_reference, bool apply_rtop_flag, ...)` | 2023 | **BLOCKED — no gemmi equivalent. gemmi graph matching to be added separately/later ([[project_gemmi_graph_matching]]); decided 2026-06-15. Defer test until then.** |
| 2 | `match_names` | `void match_names(mmdb::Residue *res_moving_names)` | 2207 | **BLOCKED — graph-match family (see #1)** |
| 3 | `matching_atoms` | `std::vector<mmdb::Atom *> matching_atoms(mmdb::Residue *residue)` | 1686 | **BLOCKED — graph-match family (see #1)** |
| 4 | `compare_residues` | `bool compare_residues(const std::pair<mmdb::Residue *, int> &a, ...)` | 5432 | predicate; borderline (used as comparator) |
| 5 | `close_residues` | `... close_residues(mmdb::Manager *mol1, mmdb::Manager *mol2, float dist)` | 7319 | spatial query across two mols |
| 6 | `residues_near_residues_for_residues` | `std::map<mmdb::Residue*, std::set<mmdb::Residue*>> residues_near_residues_for_residues(const std::map<...> &all_molecule_map, ...)` | 601 | neighbour map |
| 7 | `mol_by_symmetry` | `mmdb::Manager *mol_by_symmetry(mmdb::Manager *mol, ...)` | 7119 | returns a new molecule |
| 8 | `symmetry_move_atoms` | `... symmetry_move_atoms(const std::vector<clipper::Coord_orth> &protein_coords, ...)` | 7654 | clipper-typed in/out |
| 9 | `set_mol_cell` | `bool set_mol_cell(mmdb::Manager *mol, clipper::Cell cell)` | 7791 | mutates cell |
| 10 | `rotate_atom_about` | `void rotate_atom_about(const clipper::Coord_orth &direction, ...)` | 6258 | geometry edit |
| 11 | `standardize_peptide_C_N_distances` | `void standardize_peptide_C_N_distances(const std::vector<std::pair<mmdb::Atom*, mmdb::Atom*>> &C_N_pairs)` | 6277 | geometry edit |
| 12 | `chiral_4th_atom` | `mmdb::Atom *chiral_4th_atom(mmdb::Residue *residue_p, mmdb::Atom *at_centre, ...)` | 8557 | returns atom |
| 13 | `nucleotide_is_DNA` | `bool nucleotide_is_DNA(mmdb::Residue *r)` | 5122 | **DONE** — gemmi twin in coot-coord-utils-gemmi.{hh,cc} (unpadded "O2'"/"O2*"); test `test_nucleotide_is_DNA` PASS |
| 14 | `chains_in_atom_selection` | `std::vector<mmdb::Chain *> chains_in_atom_selection(mmdb::Manager *mol, int model_number, const std::string &atom_selection)` | 1993 | selection query |
| ~~15~~ | ~~`chain_atoms_segid`~~ | — | 1105 | **REMOVED from scope. src-only (1 caller: src/molecule-class-info-ncs.cc:1423). User is deleting it and replacing with `chain_residues_segid()` (does not exist yet). My gemmi twin + test were reverted. The inconsistent-segID case never occurs in real data.** |
| ~~16~~ | ~~`get_residue_alt_confs`~~ | — | — | **NOT a coot-coord-utils.cc function. Parser false positive: the def at .cc:3191 is inside `#if 0` ("we already have one of these in geometry"). The LIVE function is `geometry/mol-utils.cc:31` (decl mol-utils.hh:36), used 12× in reduce.cc. If it needs a test, that belongs with mol-utils, not here. REMOVED from this list.** |
| 17 | `shift_to_origin` | `clipper::Coord_frac shift_to_origin(mmdb::Manager *mol)` | 7661/7917 | mmdb overload here; the clipper overload was moved to coot-coord-utils-mmdb-free.cc |
| 18 | `translate_close_to_origin` | `void translate_close_to_origin(mmdb::Manager *mol)` | (mmdb overload) | clipper overload moved to coot-coord-utils-mmdb-free.cc |
| 19 | `position_residue_by_internal_coordinates` | `class position_residue_by_internal_coordinates { ... }` (ctor + `move_moving_residue()`) | 8447/8481 | a class, not a free fn; test via construct + move_moving_residue() |

## Suggested ordering (easy/scalar first)
13 `nucleotide_is_DNA`, 15 `chain_atoms_segid`, 14 `chains_in_atom_selection`,
17 `shift_to_origin`, 18 `translate_close_to_origin` → then the geometry edits
(10, 11), spatial queries (5, 6), symmetry (7, 8, 9), match family (1, 2, 3, 4),
12 `chiral_4th_atom`, 19 `position_residue_by_internal_coordinates`.
Resolve 16 `get_residue_alt_confs` (dead) separately.
