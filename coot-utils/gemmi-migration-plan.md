# Plan: mmdb → gemmi conversion

A roadmap for migrating Coot's molecular structure backend from mmdb2 to gemmi.

Companion to `gemmi-conversion-review.md` (which documents the 26 functions
already converted in `coot-coord-utils-gemmi.{hh,cc}`).

---

## Scope by the numbers

| Pattern | Call sites | Difficulty |
|---|---|---|
| Selection API (`NewSelection` / `SelectAtoms` / `GetSelIndex` / `SeekContacts`) | ~1,613 | high |
| `PutUDData` / `GetUDData` | ~632 | high (no native gemmi equivalent) |
| `atom->residue` back-pointer | ~172 | medium |
| `// PDBv3 FIXME` (padded atom names) | ~90 | medium |
| `coot-coord-utils.hh` functions still to port | ~75 | low–medium |

Counts are grep-derived from the main coot tree (excluding `ccp4mg-utils/`).

---

## Three cross-cutting design problems

These need a chosen approach **before** bulk conversion, otherwise each file
gets re-debated.

### Problem A: `atom->residue` (no gemmi equivalent)

`mmdb::Atom` carries a back-pointer to its `Residue` (and via that, `Chain`).
`gemmi::Atom` does not. Four approaches, in preference order:

1. **Restructure loops to carry residue context** — change
   `for (atom : flat_atom_list)` to `for (chain) for (res) for (atom)`. The
   residue is then in the enclosing scope. Most gemmi-native, removes the
   back-pointer requirement entirely. Works for the majority of the 172 sites.

2. **`gemmi::CRA` triples** — when you must carry an atom out of its scope
   (return value, container), carry `{Chain*, Residue*, Atom*}` instead of
   `Atom*`. Define `using coot::cra_t = gemmi::CRA;` for grep-ability.

3. **`NeighborSearch::Mark`** — when atoms come from a spatial query, the
   `Mark` already stores `chain_idx` / `residue_idx` / `atom_idx` and has
   `to_cra(model)`. Use this directly; don't synthesize a back-pointer.

4. **`gemmi::Model::find_cra(atom*)` linear-scan helper** — last-resort
   retrofit. Only acceptable for cold paths.

**Decision rule:** if you find yourself reaching for #4, the surrounding code
probably needs restructuring (#1) instead.

#### Inventory of `Atom*`-returning functions whose callers use `->residue`

The 259 C++ uses of `atom->residue` (and `atom->residue->chain`) trace back
to four kinds of source. Categorising them is essential because each needs a
different conversion strategy.

##### Category A — single-atom lookup functions returning `mmdb::Atom *`

These take a spec/CID and return a single pointer. Callers routinely use
`->residue` on the result. **These are the ones whose signature must change**
(return `gemmi::CRA` or `std::optional<gemmi::CRA>` instead of a bare
`Atom*`):

| Function | Header |
|---|---|
| `coot::util::get_atom(atom_spec_t, mmdb::Manager*)` | `coot-utils/coot-coord-utils.hh:773` |
| `coot::util::get_atom_using_fuzzy_search(atom_spec_t, mmdb::Manager*)` | `coot-utils/coot-coord-utils.hh:777` |
| `coot::util::get_atom(res_1, res_2, atom_name_quad, idx)` | `coot-utils/coot-coord-utils.hh:580` |
| `chiral_4th_atom(residue_p, at_centre, ...)` | `coot-utils/coot-coord-utils.hh:596` |
| `intelligent_this_residue_mmdb_atom(res_p)` | `coot-utils/coot-coord-utils.hh:1247` |
| `get_first_atom_with_atom_name(name, ...)` | `coords/mmdb-extras.hh:127` |
| `previous_non_riding_atom(atom_vector, ...)` | `coot-utils/coot-shelx.hh:280` |
| `atom_spec_t::get_atom(mmdb::Manager*)` | `geometry/residue-and-atom-specs.hh:199` |
| `coot::molecule_t::cid_to_atom(cid)` | `api/coot-molecule.hh:621` |
| `coot::molecule_t::get_atom(atom_spec)` | `api/coot-molecule.hh:655` |
| `molecules_container_t::get_atom_using_cid(imol, cid)` | `api/molecules-container.hh:1366` |
| `molecules_container_t::get_atom(imol, atom_spec)` | `api/molecules-container.hh:1385` |
| `graphics_info_t::get_atom(imol, spec)` | `src/graphics-info.h:1347` |
| `graphics_info_t::find_atom_in_moving_atoms(spec)` | `src/graphics-info.h:1494` |
| `graphics_info_t::get_moving_atom(pick_info)` | `src/graphics-info.h:1499` |

Note: `coot::util::get_atom(spec, mmdb::Residue*)` at
`coot-coord-utils.hh:781` is **not** in this list — the residue is already
the input, so the back-pointer is trivially known. Its gemmi equivalent can
keep returning `gemmi::Atom*`.

##### Category B — `mmdb::Atom**` selection arrays

These don't return a single `Atom*` but a flat array. Callers iterate with
`selection[i]->residue`. Top offenders:

- `coords/Bond_lines.cc` (~23 sites)
- `coot-utils/merge-atom-selections.cc` (~12 sites)
- `ideal/torsion-bonds.cc` (~6 sites)
- `coot-utils/coot-coord-utils.cc` (~10 sites, mostly internal selection use)

The mmdb sources are:

- `mmdb::Manager::GetSelIndex(handle, atom_selection, n_atoms)`
- `mmdb::Manager::SeekContacts(...)` → `mmdb::Contact*` whose `id1`/`id2`
  index back into an `atom_selection`

These don't translate function-by-function. They restructure into nested
`for (chain) for (res) for (atom)` loops or `gemmi::NeighborSearch` —
whose `Mark::to_cra(model)` gives back chain + residue + atom explicitly.

##### Category C — refinement's internal `atom[i]` array

`restraints_container_t::atom` is an `mmdb::Atom**` member of the refinement
container. The entry point is

- `restraints_container_t::get_atom(int i) const`
  (`ideal/simple-restraint.hh:2396`)

used throughout `ideal/` as `atom[i]->residue`,
`atom[restraint.atom_index_1]->residue`, etc. ~80+ sites across
`simple-restraint.cc`, `ng.cc`, `mods.cc`, `make-restraints.cc`,
`trans-peptide.cc`, `crankshaft.cc`.

Strategy: store `std::vector<gemmi::CRA>` (or parallel
`std::vector<int> residue_index, chain_index`) instead of `mmdb::Atom**`.
`get_atom(i)` returns the `Atom*` slot of the CRA; add a new
`get_residue(i)` / `get_cra(i)` for the rest. Combined with the
`atom.serial`-as-stable-index decision in Problem B below, this becomes a
clean refactor.

##### Category D — `at->residue->chain` (one extra hop)

A handful of sites need the chain too, e.g.

- `coot-utils/coot-coord-utils.cc:5268` — `at->residue->chain == moving_chain`
- `analysis/daca.cc:1082` — `at->residue->chain == reference_residue_p->chain`

CRA covers these natively — `cra.chain == ...` — without a second lookup.

##### Heat map

Per-file `atom->residue` counts (top 10), useful for prioritising:

| File | Count |
|---|---|
| `coords/Bond_lines.cc` | 23 |
| `coot-utils/coot-coord-utils-glyco.cc` | 21 |
| `src/graphics-info-defines.cc` | 13 |
| `ideal/simple-restraint.cc` | 12 |
| `coot-utils/merge-atom-selections.cc` | 12 |
| `ideal/ng.cc` | 10 |
| `coot-utils/coot-coord-utils.cc` | 10 |
| `src/molecule-class-info-other.cc` | 9 |
| `coot-utils/stack-and-pair.cc` | 9 |
| `ideal/trans-peptide.cc` | 8 |

### Problem B: `PutUDData` / `GetUDData` replacement (~632 sites)

`gemmi::Atom` has these spare-ish slots (from `model.hpp` lines 129–144):

| Slot | Type | Default use | Hijackable? |
|---|---|---|---|
| `serial` | `int` | PDB serial number | Yes — already standard idiom for "atom index handle" |
| `flag` | `char` | "custom flag" | Yes — boolean / small-enum tags |
| `tls_group_id` | `short` | TLS group | Yes if no TLS in play |
| `fraction` | `float` | Refmac D fraction | Yes if no D refinement |

UDD usage in coot falls into three categories:

| Category | Example | Replacement |
|---|---|---|
| Atom-index handle (refinement) | `atom[i]->PutUDData(udd_atom_index_handle, i)` then `at->GetUDData(...,ai)` to map back into a flat coord array | Use `atom.serial = i`. Stable for the duration of the refinement. The site at `ideal/simple-restraint.cc:4367` is the canonical example. |
| Boolean / enum flags | `at->PutUDData(uddHnd, 0/1)` for "fixed" / "moving" / "in residue X" | Either `atom.flag` (one char) or external `std::vector<bool>` indexed by `atom.serial` |
| Multi-handle data | Separate handles for bond / angle / torsion membership | External `std::unordered_map<int, T>` keyed by `atom.serial`, or parallel `std::vector<T>` sized to atom count |

**Decision rule:** make `atom.serial` the canonical stable index; everything
else is an external container keyed on it. Never use raw `gemmi::Atom*` as a
map key — vector reallocation invalidates the pointer.

**Risk:** `atom.serial` is also written by PDB parsing, so any code that reads
"PDB serial" must do so before refinement repurposes the field. Document this
clearly. Alternatively: introduce `coot::set_atom_indices(structure)` and
treat `serial` as ours from that point. Note that `serial` also matters for
CONECT records, ligand parsing and PDB round-tripping — repurposing it is not
free.

#### Float-typed UDD (`RegisterUDReal`) — 6 sites

In addition to the integer handles, coot uses `RegisterUDReal` to attach a
float per atom in 6 places:

| Site | Purpose |
|---|---|
| `coot-utils/coot-shelx-ins.cc:119` | Riding-atom negative U value |
| `coot-utils/q-score.hh:264` | Q-score per atom |
| `coot-utils/hole.cc:58` | Atom radius (hole calculation) |
| `coords/Bond_lines.cc:5439` | Rainbow colour position |
| `coords/Bond_lines.cc:5504` | B-factor fraction point |
| `src/molecule-class-info.cc:10577` | B-factor scale (registered on hierarchy, not atom) |

`gemmi::Atom::fraction` (float, "custom value, one use is Refmac's
ccp4_deuterium_fraction") is the natural single-payload slot — coot does no D
refinement, so it is free. Good for the single-float uses (Q-score, radius,
rainbow position). For more than one float per atom at once, the same
external-vector pattern applies: `std::vector<float>` indexed by the stable
atom index.

#### Design option: extending `gemmi::Atom` upstream

Three shapes were considered for a richer per-atom user slot. Memory cost is
the deciding factor (cryo-EM structures reach ~1M atoms):

| Approach | Per-atom cost (unused) | Per-atom cost (used) | Verdict |
|---|---|---|---|
| `int user_defined_int` | 4 B always | 4 B | Narrow (one slot) but trivial; doesn't solve multi-handle |
| `void *user_data` | 8 B always | 8 B + caller payload | Most flexible; C-API-standard "hang anything here"; best upstream pitch |
| `std::unordered_map<int,int> user_int_map` per atom | ~56 B always (empty) | ~56 B + ~32 B/entry | Faithful to mmdb's model but worst layout: ~56 MB of empty maps on a 1M-atom structure, pointer indirection + allocator pressure per access. Rejected. |

The per-atom map is the *same data shape* as the recommended external
container but with the worst possible memory layout — it taxes every atom
that has no UDD at all. The external `std::unordered_map` keyed on the stable
atom index is pay-only-for-what-you-use and needs no gemmi change.

**Plan of record:**

1. Propose `void *user_data` (or `uint32_t user_tag`) upstream to gemmi,
   pitched generically ("one stable per-atom hook for downstream tools,
   equivalent to a single mmdb `UDR_ATOM` handle") — not as a coot-specific
   feature. Best chance of acceptance; serves any tool.
2. Design the migration **as if** it lands. If accepted, the heaviest single
   handle (`udd_atom_index_handle`) and the single-float cases hang off it
   cleanly without touching `serial`.
3. Fall back to `atom.serial` (int) / `atom.fraction` (float) only if the
   upstream proposal stalls.
4. Either way, the multi-handle pattern (refinement registering several int
   handles on the same atom simultaneously) still uses an external
   `std::vector<T>` / `std::unordered_map<int,T>` keyed by the stable atom
   index. No gemmi change displaces that.

### Problem C: Selection API (~1,613 sites)

mmdb's `NewSelection` / `SelectAtoms` / `GetSelIndex` returns `mmdb::Atom**`
flat arrays. Replacements by pattern:

| mmdb pattern | gemmi replacement |
|---|---|
| Whole-molecule atom scan | nested `for (chain) for (res) for (atom)` |
| Single chain | `for (auto &res : chain.residues)` |
| Spatial / contacts | `gemmi::NeighborSearch::find_neighbors` (returns `Mark*`) |
| Atom-spec lookup | `chain.find_residue_group(SeqId)` then `res.find_atom(name, altloc)` |
| Text-spec ("/A/12-30/CA") | `gemmi::Selection` from `select.hpp` |
| "All atoms in residue list" | direct iteration over `res.atoms` |
| `SeekContacts` | `NeighborSearch::populate()` once, then `find_neighbors` per query atom |

**Decision rule:** convert by pattern in batches, not file-by-file. After
Problems A and B are solved, most selection sites become trivial.

---

## Other problem corners

- **SSM (structure superposition)** — no gemmi equivalent. Either keep the
  mmdb sub-dependency just for SSM, or replace with clipper or a hand-rolled
  Kabsch. Park this until last.
- **TER atoms** — gemmi has no TER atoms in memory. Chain ends are implicit
  (end of `chain.residues`, or the polymer / non-polymer boundary determined
  by `entity_type`). On PDB **output**, the writer synthesises TER records
  from that structural information. From `~/gemmi/include/gemmi/to_pdb.hpp`:

  ```cpp
  bool ter_records = true;         // Write TER records (default on)
  bool numbered_ter = true;        // TER gets its own serial number
  bool ter_ignores_type = false;   // Put TER after last atom in chain
                                   // (even water) — reverts to the older
                                   // "TER after the last atom regardless"
                                   // behaviour
  ```

  Default placement: TER goes after the last polymer residue of a chain;
  non-polymer (waters, ligands) follow without a preceding TER. Set
  `ter_ignores_type = true` for the older behaviour.

  **Migration implication:** mmdb code routinely does `if (!at->isTer())`
  checks to skip dummy TER atoms in iteration — these are unnecessary and
  meaningless when iterating gemmi residues. They should just be **removed**
  during conversion, not translated. Also: anything that depends on
  `entity_type` being correctly populated for sensible TER placement should
  verify that `copy_from_mmdb` (or the upstream PDB read) sets it.
- **PDBv3 atom names** — already understood, ~90 FIXME markers. Resolved
  per-file as conversion proceeds. All gemmi code uses unpadded names
  natively.
- **`copy_from_mmdb` gaps** — secondary structure, headers, resolution,
  biomolecule info not currently copied. Extend or accept loss for now.
- **`label_atom_id` padding bug** — `trim_atom_names()` workaround in place;
  upstream fix to gemmi still TODO (propose to Marcin).

---

## Phase plan

### Phase 1 — Finish the utility layer in `coot-coord-utils-gemmi`

~25 more functions, in this order (each with a comparison test against mmdb):

1. `get_residue(chain_id, res_no, ins_code, st)` — wraps
   `chain.find_residue_group(SeqId)`
2. `get_atom(atom_spec_t, st)` — trims spec name, then `res.find_atom`
3. `residues_near_position(pt, st, radius)` — `NeighborSearch::find_neighbors`
4. `residues_near_residue(res, st, radius)` — same, looped over reference
   residue's atoms
5. `residues_near_residues(st, dist)` — full neighbour map
6. `get_waters`, `get_hetgroups`, `atoms_with_zero_occupancy`
7. `deep_copy_residue`, `copy_and_delete_hydrogens`
8. `sort_chains`, `sort_residues`, `shift_to_origin`,
   `move_waters_to_last_chain`

### Phase 2 — Lock down the cross-cutting decisions

- Write a short "gemmi conversion conventions" doc covering Problems A/B/C.
- Add `coot::set_atom_indices(gemmi::Structure&)` helper that walks the
  structure assigning `atom.serial = global_index`.
- Add convenience `coot::cra_t` alias.
- Pick one small file (e.g. `analysis/stats.hh` or a chunk of
  `coot-utils/peak-search.cc`) and convert it end-to-end as a worked example
  demonstrating the conventions. This becomes the reference for everyone.

### Phase 3 — File-by-file migration, easiest first

Order by reverse-dependency depth (leaves first):

1. `analysis/` — validation, mostly read-only
2. remaining `coot-utils/` files
3. `density-contour/` — already mostly clipper
4. `MoleculesToTriangles/` — coordinate-driven
5. `ligand/` — heavy selection use
6. `lidia-core/` — partly RDKit already
7. `geometry/` — restraints, padded names everywhere
8. `ideal/` — refinement, heaviest UDD use; hardest
9. `api/` — `molecules_container_t` method-by-method
10. `src/` — GUI layer, last

For each file: keep mmdb at the public interface initially (caller passes
`mmdb::Manager*`, the function does `copy_from_mmdb` internally). Once enough
callers are converted, flip the signature to take `gemmi::Structure&`.

### Phase 4 — Storage migration in `molecule_t` / `molecule_class_info_t`

Add `gemmi::Structure structure_;` alongside the existing
`mmdb::Manager *atom_sel.mol`.

Three strategies:

- **Dual-storage with sync** (start here): both representations live, gemmi
  rebuilt from mmdb on modification. Lets gemmi paths be validated against
  mmdb.
- **gemmi-primary, mmdb lazy** (when ≥80% of read paths converted): gemmi is
  authoritative; mmdb built on demand for unconverted code.
- **mmdb-only legacy mode** (escape hatch): for SSM and any unported code.

### Phase 5 — Cut mmdb

Remove the conversion layer. Drop mmdb dependency (or keep it isolated to
SSM).

---

## What to do this week

1. Write `get_residue` / `get_atom` / `residues_near_position` /
   `residues_near_residue` + tests (~1–2 days)
2. Write the conventions doc + `set_atom_indices` helper (~half day)
3. Convert one small file end-to-end as a worked example (~1 day)
