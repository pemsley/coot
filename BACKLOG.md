# Backlog

## Tier D: per-bond/per-angle SQL lookups for the AceDRG tables

The current SQLite backend (`coot-utils/acedrg-sqlite-tables.cc`) only
accelerates *loading*: `prefetch_for_molecule()` pulls whole hash-buckets
from `acedrg.sqlite` and feeds them into gemmi's in-memory caches through
the `insert_bond_row()`/`insert_angle_row()` seam (gemmi PR #438), after
which gemmi's `search_bond_multilevel()`/`search_angle_multilevel()` walk
those maps exactly as before. Bucket granularity is the ceiling: the
aromatic (sp²-C, sp²-C, in-ring) bucket alone is ~159k rows in the CCP4-9
COD data, so the measured 17× speedup on trivial molecules degrades to
~2-5× on real aromatic ligands (ATP 3.3×, A1DET 2.0×). The selectivity
that would shrink the fetch lives in the COD atom-type and neighbour
columns, which only exist mid-`fill_bond` — not at prefetch time.

Tier D replaces the walk itself: at each fallback level (12 for bonds, 6
for angles incl. the 5D wildcard), issue one indexed SELECT returning ≤1
row (~10 µs) instead of scanning the bucket. ~530 queries ≈ 5 ms per
molecule against the current 1.5-5 s. Detailed implementation plan (May
2026, written pre-split for the gemmi fork):
`/hard-disk-partition-2/files/claude-play/notes/2026-05-31-acedrg-per-bond-sql.md`.

What survives from that plan unchanged: the per-level SQL statements and
prepared-statement cache, the level semantics (incl. the 5D wildcard LIKE
trick and the ordered-map tie-break), the index pre-flight (Task 1), the
byte-identical 16-monomer harness discipline, and the open question about
where `min_observations_*` thresholds apply.

What the coot/gemmi split changes: the plan's Tasks 3-7 all edit gemmi's
`acedrg_tables.cpp`, which we no longer patch. The per-level lookup code
moves into `coot-utils/acedrg-sqlite-tables.cc`; the dispatch inside
`search_*_multilevel()` needs a NEW gemmi seam, since the fine-grained
COD atom types (`a1_type_f` etc.) are computed inside gemmi mid-search
and coot cannot reproduce them:

- Proposed seam (superseding the per-level idea — full spec in the next section): whole-function delegation. Gemmi defines plain-data `AcedrgBondSearchArgs`/`AcedrgAngleSearchArgs` structs and one optional `std::function` hook per search function; the hook is authoritative when set, byte-identical fallback when unset. SQLite-free — same upstreamability argument as the #438 insert-row seam.
- Coot implements the hook with the plan's prepared-statement SQL against `acedrg.sqlite` (built by `coot-make-acedrg-sqlite`; verify the indexes there cover the per-level WHERE clauses — plan Task 1).
- Once verified byte-identical, `prefetch_for_molecule()` becomes dead weight for bonds/angles and is dropped (plan Task 7); HRS-backed levels (bond 9-11, angle 6) stay in-memory — those tables are small and eagerly loaded.
- Fallback duplication of the ladder wholly inside coot (no seam) was considered and rejected: it would fork several hundred lines of intricate fallback logic and still couldn't get the mid-search atom types.

Sequencing: (1) design + prototype the search seam on a gemmi fork branch, (2) propose it upstream alongside nudging #437/#438, (3) implement the coot-side lookups against the agreed signature, (4) harness, timing, drop prefetch. Real-world payoff on this box: pyrogen per-dictionary time in the tautomer pipeline (44 dictionaries/run) is dominated by exactly this prefetch cost.

## Tier D hand-off spec: the AceDRG search seam (gemmi side)

Hand-off spec for the gemmi session (2026-09-14). Defines the seam gemmi
must grow so that coot can implement per-bond/per-angle SQL lookups
(Tier D above) without gemmi depending on SQLite. Written against
`/home/paule/Projects/gemmi/fork-5/gemmi`, branch `master`, commit
`2e4998e9` (upstream official-master merged; includes the
`insert_bond_row`/`insert_angle_row` seam of gemmi PR #438). All
file:line references below are to that tree.

The old (pre-split) implementation plan
`/hard-disk-partition-2/files/claude-play/notes/2026-05-31-acedrg-per-bond-sql.md`
has its Tasks 3-7 superseded by this spec; its SQL/index/harness material
is still useful to the coot side.

### Design decision: whole-function delegation, not per-level hooks

The May plan proposed per-level lookup hooks. That no longer fits: the
current `search_bond_multilevel()` (src/acedrg_tables.cpp:2509) entangles
control flow with in-memory map *presence* — the `start_level` gating
(:2633-2643), the in-ring Y/N fallback probe against
`bond_full_4prefix_keys_` (:2534-2547), and the `has_a1_class_only`
post-search override branch (:2837-2867). `search_angle_multilevel()`
(:3040) likewise has the same-hash property-swap retry loop (:3106-3122)
and the 5D relaxed match (:3192-3211). Per-level hooks would need probe
hooks for all of these.

Instead: gemmi grows ONE optional hook per search function. When the hook
is set, the search function computes its derived keys exactly as today,
calls the hook, and returns the hook's answer. All fallback semantics
(levels, gating, probes, overrides, thresholds) become the hook
implementer's responsibility, verified byte-identical by the harness.
When the hook is unset — or declines — behaviour is bit-identical to
today. This keeps the gemmi patch small, dependency-free, and
upstreamable on the same argument as PR #438.

### Gemmi-side deliverable

**1. Two plain-data args structs (public, in acedrg_tables.hpp).**
Field values are exactly the locals the search functions already compute.

```cpp
// Arguments handed to an external bond-search hook. Atoms are already
// ordered by order_two_atoms() (hash, then cod_main, then id — see
// src/acedrg_tables.cpp:1215); "1"/"2" below mean left/right in that
// order, matching the bond-table column convention.
struct AcedrgBondSearchArgs {
  int ha1 = 0, ha2 = 0;        // hashing values, ordered
  std::string hybr_comb;       // "SP2_SP2" etc., two names sorted, '_'-joined
  bool same_ring = false;      // are_in_same_ring(a1, a2) — RAW, pre-probe
  std::string a1_nb2, a2_nb2;  // nb2_symb
  std::string a1_nb,  a2_nb;   // nb1nb2_sp
  std::string a1_type, a2_type;    // cod_main
  std::string a1_class, a2_class;  // cod_class_no_charge
  int num_th = 3;              // min_observations_bond, or 1 for As/Ge pairs
};

// Arguments for an external angle-search hook. ha1 is the CENTER atom's
// hash; flanks are ordered by order_two_atoms() (min = "2", max = "3"),
// matching the angle-table column convention.
struct AcedrgAngleSearchArgs {
  int ha1 = 0, ha2 = 0, ha3 = 0;   // center, flank_min, flank_max
  int ring_val = 0;                // angle_ring_size(center, fmin, fmax)
  std::string center_hybr;         // hybridization_to_string(center) — for the 5D relaxed match
  std::string value_key;           // cat(ring_val, ':', hc_'_'-joined-sorted-flank-hybrs)
  std::string a1_root, a2_root, a3_root;   // cod_root
  std::string a1_nb2,  a2_nb2,  a3_nb2;    // nb2_symb
  std::string a1_nb,   a2_nb,   a3_nb;     // nb_symb
  std::string a1_type, a2_type, a3_type;   // cod_main
  int min_obs = 3;                 // min_obs_eff (see :3091)
};
```

**2. Two hook members + setters on AcedrgTables.**

```cpp
// External search hooks (e.g. a SQL backend living outside gemmi).
// Contract: return true when the hook is AUTHORITATIVE — `out` then
// carries the complete search result, including the no-match sentinel
// (bond: out.level == -1; angle: default CodStats with NaN value, and
// *out_level untouched i.e. left at the caller's initial value).
// Return false to decline (backend unavailable, error): gemmi then
// falls through to the normal in-memory walk.
std::function<bool(const AcedrgBondSearchArgs&, CodStats& out)>
    external_bond_search;
std::function<bool(const AcedrgAngleSearchArgs&, CodStats& out,
                   int* out_level)> external_angle_search;
```

Public members (like `verbose`) or private with setters — implementer's
choice; public members are simpler and match the class's existing style.

**3. Dispatch points.**

- `search_bond_multilevel()` — after the derived-key block ends at src/acedrg_tables.cpp:2565 (`num_th` including the As/Ge exception) and before the verbose print/key construction at :2567. Build `AcedrgBondSearchArgs` from the locals `ha1, ha2, hybr_comb`, the raw `are_in_same_ring(a1, a2)` result (NOT the possibly-swapped `in_ring` string — move the dispatch above the probe at :2534, or re-derive), `a1_nb2 … a2_class`, `num_th`. Then `if (external_bond_search) { CodStats out; if (external_bond_search(args, out)) return out; }`.
- `search_angle_multilevel()` — after the type strings at :3076, before the verbose print at :3078. `min_obs` field = `min_obs_eff` per :3091. `if (external_angle_search) { CodStats out; if (external_angle_search(args, out, out_level)) return out; }`.
- No other call sites change. `search_bond_hrs`, `search_bond_en`, metal bonds, torsions, prot-hydr distances all stay in-memory (their tables are small and eagerly loaded).

**4. Acceptance criteria (gemmi side).**

- Hooks unset ⇒ zero behaviour change (byte-identical `gemmi drg` output on the 16-monomer harness, ASCII backend).
- A trivial test hook that returns false ⇒ still byte-identical.
- Compiles without sqlite3 anywhere in gemmi.
- Commit in the style of `70087d54` (the insert-row seam), PR-able upstream alongside #438.

### Semantics inventory (what the coot-side hook must reproduce)

The hook must return exactly what the in-memory walk would. Read the two
functions in full before implementing; the inventory below is the map,
not the territory.

**Bond (src/acedrg_tables.cpp:2509-2878):**

- In-ring probe (:2533-2547): `in_ring` = same_ring ? "Y" : "N", but if no rows exist for `(ha1, ha2, hybr_comb, in_ring)` and rows exist for the flipped value, flip it. In SQL: two EXISTS probes.
- Exact-codClass pre-check (:2581-2599): row at (full 10-part key + a1_class + a2_class) with `count >= num_th` returns as level 0. If the a1_class node exists but the a2_class node doesn't, set `has_a1_class_only` — it changes the endgame (below).
- `start_level` gating (:2633-2643): computed from key presence (has_a2_type ⇒ 1, has_a1_type ⇒ 2, a2_nb ⇒ 4, a1_nb ⇒ 5, a2_nb2 ⇒ 7, a1_nb2 ⇒ 8, in_ring ⇒ 10, hybr ⇒ 11). Presence = EXISTS probes in SQL. Note gating is not merely an optimisation: starting at level N skips the more-specific levels even if their queries would coincidentally aggregate something.
- Levels (loop :2654-2835), with the containers they walk:
  - 0: `bond_idx_1d_[key_8][a1_type][a2_type]` — `.front()` of the vector (FIRST-LOADED row: SQL must preserve/order by insertion rowid). Threshold: `vs.count >= num_th` (observation count).
  - 1: fronts of every a2_type group under a1_type, PLUS fronts of (other-a1_type, ==a2_type) groups; aggregate. Threshold: number of contributing entries `>= num_th`.
  - 2: fronts of (a1_type != given, a2_type == given) groups; aggregate; entry-count threshold.
  - 3: all rows at `bond_idx_2d_[key_4][a1_nb2][a2_nb2][a1_nb][a2_nb]`; aggregate; entry-count = row count.
  - 4: all rows under a1_nb (any a2_nb) plus rows under (other a1_nb, ==a2_nb); aggregate.
  - 5: only rows under (other a1_nb, ==a2_nb); aggregate.
  - 6: every row under `[key_4][a1_nb2][a2_nb2]`; aggregate.
  - 7: every row under `[key_4][a1_nb2]` (any a2_nb2 …) plus rows under (other a1_nb2, ==a2_nb2); aggregate.
  - 8: only rows under (other a1_nb2, ==a2_nb2); aggregate.
  - 9/10/11: `bond_hasp_2d_[key_4]` / `bond_hasp_1d_[ha1|ha2|hybr]` / `bond_hasp_0d_[ha1|ha2]` vectors; aggregate; threshold is mere presence.
- Endgame (:2837-2870): normally return the first level that met threshold. But when `has_a1_class_only`: prefer exact (a1_type, a2_type) front if `count >= num_th` (as level 0), else exact (a1_nb, a2_nb) aggregate regardless of threshold (as level 3), else the matched level. No match at all ⇒ `level = -1` sentinel.
- Aggregation formula (`aggregate_stats`, :2359-2384): weighted mean `SUM(v*n)/SUM(n)`; pooled sigma `sqrt(|SUM((n-1)*s² + n*v²) − mean*SUM(v*n)| / (N−1))` (divide by N when N==1); count = `SUM(n)`. Single-entry input returns the entry unchanged. Expressible as SQL SUMs with the final arithmetic in C++.

**Angle (src/acedrg_tables.cpp:3040-3226):**

- All levels are flat compound-key lookups returning `.front()` (first-loaded row) — no aggregation anywhere. Threshold `front().count >= min_obs` at levels 1D-5D; 6D is presence-only.
- Same-hash retry (:3106-3122): when ha1==ha2 or ha1==ha3, levels 1D-4D are tried a second time with the center's property columns swapped with the equal-hash flank's. The 5D/6D levels are OUTSIDE the retry loop.
- Levels and keys: 1D = 16-part key (reported level 0), 2D = 13-part (1), 3D = 10-part (2), 4D = 7-part (3), 5D = `[ha1][ha2][ha3][value_key]` (4), 6D = `ha1|ha2|ha3` (6). Level 5 is never reported.
- 5D relaxed match (:3192-3211): ONLY when the exact value_key is absent (an under-threshold exact entry falls through without relaxing). Scan value_keys in LEXICOGRAPHIC order (std::map order — SQL: ORDER BY value_key) for the first with prefix `ring_val:` whose first '_'-delimited hybr token equals `center_hybr`; try only that one, then stop.
- Miss ⇒ default `CodStats()` (NaN value), out_level untouched.

**Cross-cutting:**

- front() everywhere means the SQLite DB must preserve per-key row order matching the ASCII load order (rowid ASC from the existing builder is exactly that — verify `coot-make-acedrg-sqlite` inserts in file order).
- Level numbers matter downstream: `fill_bond` (:2480) treats `level >= 9` like HRS (count-only acceptance) and returns the level as the approx-level; get them exactly right.
- Verification: the 16-monomer byte-identical harness (`/hard-disk-partition-2/files/claude-play/experiments/compare_sqlite_vs_ascii.py`) is the gate at every step, exactly as in the May plan. Timing targets also as there (A1DET from ~5.6 s to under ~1 s).

### Tier D open questions

- Public members vs setters for the hooks (cosmetic; decide in gemmi session).
- Whether coot precomputes per-level aggregate tables at DB build time (turns aggregation levels into single-row SELECTs) or aggregates at query time with SQL SUMs. Query-time is simpler and the aggregation levels are rarely reached; start there.
- The May plan's open question about threshold semantics is answered by the code read: bond levels 1-8 threshold on contributing-entry count, level 0 on observation count, 9-11 on presence; angles threshold on the front row's observation count. The inventory above is authoritative as of commit 2e4998e9 — re-verify if gemmi master moves.

## Replace quick_protein_donor_acceptors with a dictionary-aware donor/acceptor lookup

`coot::quick_protein_donor_acceptors` (geometry/protein-donor-acceptors.{hh,cc})
is a hard-coded (residue-name, atom-name) → hb_type map for the standard amino
acids only. It is used by the environment-distances display ("Constructor G" in
coords/Bond_lines.cc, ~line 3000) to decide whether a close contact is drawn as
a hydrogen bond (purple) or a plain contact.

Because it knows nothing about ligands, any non-protein atom pair
fails the lookup and gets the fallback colour — so ligand H-bonds are never
properly classified (e.g. the t3ND4 lactam O27 ↔ ASP main-chain O pair,
acceptor + acceptor, 2026-09-13). The replacement should get hb types from the
dictionaries: `protein_geometry::get_h_bond_type(atom_name, monomer_name,
imol_enc)` already maps `_chem_comp_atom.type_energy` → hb_t via the energy
lib, and `protein_geom_p` is in scope at the call site. The quick table could
remain as a fast path for protein atoms, or go entirely.

Quick-table defects all fixed 2026-09-13: THR " OG1", THR in the main-chain list, HOH " O  " = HB_BOTH (explicit entry after the loop, overriding the loop's blanket HB_ACCEPTOR), and `n_res_types` now comes from `std::size(l)`. Consider `#include <iterator>` for `std::size` portability (GCC's `<iostream>` provides it transitively; libc++ may not).
