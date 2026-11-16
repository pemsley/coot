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

## Speed up average-map generation by checking the grids

`on_make_an_average_map_ok_button_clicked()` (src/glade-callbacks-main-window.cc)
currently uses the general `coot::util::average_map()` route, which does not
require the contributing maps to share a grid sampling. It could be speeded up
by first checking that the grids of the contributing maps match, and if so using
the `regen_map` route (as `regen_map_py()` in src/c-interface-maps.cc does),
which averages in place and presumes matching grids.

## Add property histograms behind the Layla QED desirability curves

The Layla QED tab draws the 8 ADS desirability curves d(x) (layla/ui.cpp,
`draw_desirability_curve()`; parameters in layla/qed.cpp). At some later stage
we would like to draw, behind each curve, the histogram of that property over
the reference set of 771 oral drugs from Bickerton et al. 2012 (as in Fig. 1 of
the paper) — light-blue bars under the curve.

Where the histogram data is: the paper's supplementary Excel workbook
`NIHMS50746-supplement-QEDExampleExcel.xlsx`, sheet "Desirability Functions",
773 rows = header + 771 drugs. Columns: A=Drug name, then the raw per-drug
properties B=MW, C=ALOGP, D=HBA, E=HBD, F=PSA, G=ROTB, H=AROM, I=ALERTS
(J-Q are the individual desirabilities, R/S the QED scores). This column order
matches the `QED::QEDPropName` enum. Histogramming columns B-I reproduces the
Fig. 1 bars. (The SuppInfo PDF and the main-text PDF do NOT contain this raw
data — only the fitted ADS parameters in Suppl. Table 2.)

To ship the bars we'd need to bake the 771-drug values (or precomputed bin
counts per property) into the code/data, since we can't depend on the Excel
file at runtime.


## Use a constant view direction for specular under orthographic projection

The specular term in the lit shaders (`moleculestotriangles*.shader`,
`meshes*.shader`, `instanced-objects.shader`,
`instanced-meshes-with-shadows.shader`, `texture-meshes-with-shadows.shader`)
computes `view_dir = normalize(eye_position_in_molecule_coordinates_space - frag_pos)`
from a finite eye position, set by `Shader::setup_eye_position()` as
`graphics_info_t::eye_position · view_rotation + rotation_centre`.

That is right for a perspective projection, but Coot's default is orthographic
(`perspective_projection_flag == false`), where the true view direction is the
same for every fragment: the view-space +z axis rotated into molecule space.
Using a finite eye there makes the highlights shift a little with zoom (the
right-drag zoom scales `eye_position`), and a surface close to the eye gets a
wider spread of view directions than it should.

Fix: when the projection is orthographic, send a direction instead of a
position (e.g. a `view_direction_in_molecule_coordinates_space` uniform, the
(0,0,1) view vector times the view rotation) and use it directly as `view_dir`.
`mvp_is_perspective()` in `Shader.hh` / the `is_perspective_projection`
uniform (now set by the Mesh, TextureMesh and Instanced_Markup_Mesh draw
functions) can select between the two in the shader.

Deferred 2026-10-13 after the specular/eye-position fixes, as the remaining
error is small.

## "Open Recent..." menu

The usual desktop "Open Recent" assumes one-file-equals-one-document, which
Coot breaks: a working session is a *set* of files (coordinates, reflection
data or a map, often a ligand dictionary, sometimes several of each). So two
tiers are worth having.

### Tier 1 — Recent files (flat MRU)

A per-file most-recently-used list, newest first, grouped or icon-tagged by
type (coordinates / data / map / dictionary). Expected from the menu label,
cheap to implement, and the right tool when the user just wants to pull one
more file into an existing session. Things to get right:

- **Stat before offering.** Crystallography pipelines move and delete files
  constantly; grey out or annotate missing entries rather than failing on
  click.
- **Dedupe by canonical path** and show enough of the path to disambiguate —
  people have forty `refmac.mtz` in forty directories.
- **Remember the MTZ columns.** An `.mtz` isn't self-describing; reopening it
  usefully means recording *which columns* were used last time (FWT/PHWT,
  F/SIGF, …) so we don't re-prompt the column dialog every time. That column
  choice must be part of the recent entry, not just the path. This is the
  single biggest difference from a text editor's Open Recent.

### Tier 2 — Recent sessions (file-set + view state)

A "Recent Sessions" submenu where each entry restores the whole working set
the user had. This is a nicer front-end onto the existing state-script
machinery (`0-coot.state.*` / `state.scm`): each entry is a saved/auto-saved
state script, clicking it replays the script. Label entries meaningfully
(derived name — PDB id or coord filename — timestamp, maybe molecule count),
not the opaque script path. This is the "offer to restore, don't impose"
behaviour from the startup discussion: it restores *inputs* the user
knowingly reopens, not an in-memory snapshot that could silently diverge from
real modelling work.

Illustrative menu shape:

```
File
  Open Recent
     Coordinates ▶   3pdz.pdb, 6g7k.pdb, …
     Maps/Data   ▶   refmac.mtz (FWT/PHWT), …
     Dictionaries▶   LIG.cif, …
  Recent Sessions ▶  3pdz — 2 maps  (10:14 today)
                     6g7k            (yesterday)
```

### Cross-cutting concerns

- **Moved projects.** Crystallographers rsync whole directories between
  machines, so sessions with absolute paths break. Store paths relative to a
  project root (and remember that root) to survive the move.
- **Where the list lives.** Per-user (`~/.coot` / preferences), capped
  (~10–15), pruned of dead entries lazily, so it survives the blank-startup
  model.
- **Privacy / clutter.** Provide "Clear Recent" and ideally a way to not
  record a given open (people demo with confidential structures).
- **Don't double-load.** If a file/session is already open, offer to
  raise/switch rather than loading a second copy.

Minimum that delivers most of the value: Tier 1 flat MRU *with remembered MTZ
columns* — self-contained, and it kills the re-prompting annoyance. Tier 2 is
the bigger conceptual win but is low-risk once the Tier 1 metadata format can
carry a "this entry is a session" variant.

## Edit Backbone Torsions Cancel button does not clear the moving atoms

Low priority.

`on_edit_backbone_torsion_cancel_button_clicked()` (src/glade-callbacks.cc)
has a comment saying that `clear_moving_atoms_object()` is "done as part of
window destroy callback", but the button only hides the dialog, so the dialog
is never destroyed and the moving atoms are never cleared. The other call it
makes, `destroy_edit_backbone_rama_plot()`, ends up in
`graphics_info_t::destroy_edit_backbone_rama_plot()` (src/graphics-info.cc),
which is an empty function.

A window-manager close does clear them: since 2026-11-01 that goes through
`on_edit_backbone_torsions_dialog_close_request()`, which calls
`clear_moving_atoms_object()` and hides the dialog (the old `destroy` handler
became a `close-request` handler so the builder dialog can be reopened).

Fix: call `clear_moving_atoms_object()` in the Cancel callback, matching
`on_edit_chi_angles_cancel_button_clicked()`, which calls
`clear_up_moving_atoms()`.

## Helix placement scoring (ligand/helix-placement.cc)

Found 2026-11-13 while reading
`coot::helix_placement::place_alpha_helix_near_kc_version()`, the function
behind `place_helix_here()`.

### Changed 2026-11-13, not yet compiled or tested

- `score_helix_position()` and `decompose_helix_by_cbeta()` looped with
  `ires < max_residue_number()`, which skipped the last residue
  (`max_residue_number()` is inclusive, see mini-mol/mini-mol.hh). Both now use
  `<=`.
- `score_helix_position()` now samples with
  `coot::util::density_at_point_by_linear_interpolation()` instead of the cubic
  `density_at_point()`.

Both changes shift the scores a little, and the cut-offs that use them were
tuned with the old behaviour: the `final_score < 0` rejection, the CB
mean/stddev test (`descrimination_ratio` 0.5), and the "Reverse" helix
trimming (the 2.0 and 2.5 factors). `score_helix_position()` is also the
scorer in `optimize_rotation_fit()`, `place_alpha_helix_near()`, `fit_strand()`
and `find_strand_candidates_by_shift_sampling()`, so helix and strand placement
both need re-testing.

### Still open

- **The score is a sum, not a mean.** `score_helix_position()` adds up
  occupancy-weighted density over all atoms, so it scales with the number of
  residues. Scores for the 20-residue and 9-residue attempts in
  `place_helix_here()` are not comparable, and neither are scores before and
  after `trim_and_grow()`. Dividing by the summed occupancy would fix that, but
  the callers that compare against absolute values need checking first.
- **Mixed interpolation.** `score_atoms()` (the CB test) and `score_residue()`
  (the end trimming in `trim_end()`) still use cubic interpolation, while
  `score_helix_position()` is now linear. Decide whether they should follow.
- **Comment disagrees with the code.** In
  `place_alpha_helix_near_kc_version()` the comment says the centre movement is
  limited to 1.5A, but `max_allowed_move` is 1.0.

## Implement the post-intermediate-atoms-moved hook

`register_post_intermediate_atoms_moved_hook()` (`src/c-interface-refine.cc`,
declared in `src/cc-interface.hh`) stores a Python function that is meant to be
run when the intermediate (moving) atoms have moved. It has never worked, and
the header doc now says that it is not implemented. Nothing in `python/` or
`scheme/` registers a hook.

What is missing, all in `src/graphics-info-intermediate.cc`:

- **The call site.** The call in `drag_refine_refine_intermediate_atoms()` is commented out with "put this somewhere else (after the refinement has finished)". That dates from the "Threaded refinement works" commit (0d87edebe).
- **The runner.** `run_post_intermediate_atoms_moved_hook_maybe()` builds the intermediate-atoms bonds representation and stops at "do something with python". It never calls the stored function.
- **Reference counting.** `graphics_info_t::register_post_intermediate_atoms_moved_hook()` stores the `PyObject *` without `Py_INCREF` (and does not `Py_XDECREF` a previous hook), so a lambda passed in can be freed while Coot still holds the pointer. It should also check `PyCallable_Check()` and drop its debug print.

To decide before implementing:

- **What the function is called with:** no arguments, or the bonds representation that the runner already builds.
- **When it is called:** after each refinement step is drawn, or once when the refinement finishes.

The call must be made on the main thread (with the GIL), so the place for it is
where the intermediate atoms are redrawn after a refinement step, not inside
the refinement thread. Check the result of the call for a Python error and
print it, so that a broken hook does not fail silently.

## Draw molecular symmetry

`add_molecular_symmetry()`, `add_molecular_symmetry_from_mtrix_from_file()` and
`add_molecular_symmetry_from_mtrix_from_self_file()` (`src/cc-interface.hh`),
and the `REMARK 350 BIOMT` reader that runs when a model is read
(`molecule_class_info_t::add_molecular_symmetry_matrices()`,
`src/molecule-class-info-ncs.cc`), all fill
`molecule_class_info_t::molecular_symmetry_matrices`. Nothing reads it, so
nothing is drawn. The header docs now say so.

- **Drawing.** The only drawing code is the `display_molecular_symmetry` lambda in `molecule_class_info_t::display_symmetry_bonds()` (`src/molecule-class-info.cc`). That is immediate-mode OpenGL (`glBegin(GL_LINES)`), its call is commented out, and `display_symmetry_bonds()` itself is never called. The current renderer draws crystal symmetry with `molecule_class_info_t::draw_symmetry()` and `Mesh::draw_symmetry()` and has no molecular-symmetry path. It needs one: a mesh per operator, or instanced drawing of the molecule's bonds mesh with one transformation per operator. The old code gave each copy its own colour.
- **The MTRIX origin.** `add_molecular_symmetry_from_mtrix_from_file()` (`src/molecule-class-info-molecular-symmetry.cc`) stores half of the MTRIX translation as the origin. For an operator `x' = Rx + t` the fixed point `o` satisfies `(I - R) o = t`, and `o = t/2` holds only for a 2-fold. The general case needs a least-squares solve, because `I - R` is singular along the rotation axis. A screw component cannot be expressed as a rotation about a point at all, so storing the full `clipper::RTop_orth` for each operator (as MTRIX and BIOMT give it) would be simpler and more general than matrix-plus-origin.
- **Remove the old code** (`display_symmetry_bonds()`) when the new drawing is in place.

## Review the nomenclature-errors prompt on reading a model

`show_fix_nomenclature_errors_gui()` (`src/c-interface-build-gui.cc`) returns at
once ("this function crashes, I don't know why - fix later", from the GTK3 port
in 2019). The dialog that it would show, `fix_nomenclature_errors_dialog` with
`nomenclature_errors_vbox`, exists only in `coot-gtk2.glade`; no GTK4 `.ui`
file has it, so removing the `return` would not bring it back.

The consequence: `set_nomenclature_errors_on_read()` defaults to "prompt", and
the prompt is this function. So by default, when a model with nomenclature
errors is read, nothing is shown and nothing is corrected. Only "auto-correct"
does anything. The menu item "Fix Nomenclature Errors..." still works.

First decide whether the dialog is still useful. Then one of:

- **Rebuild it** in `ui/coot.ui`. The OK and Cancel handlers (`on_fix_nomenclature_errors_ok_button_clicked` and `on_fix_nomenclature_errors_cancel_button_clicked` in `src/glade-callbacks.cc`) are still there and would serve it as they are.
- **Replace it** with a message that gives the number of errors and points to the menu item.
- **Remove it:** delete the function, the two handlers and the "prompt" mode, and choose what the default mode should then be ("ignore" is what happens now).

## Tidy up the API for Coot 2

The API is not to change for the current release. Coot 2 can change it, and
these are the things to tidy up then. They were found in the API documentation
review (November 2026), and each is documented in its header as it behaves now.

### Parameters that are accepted and never used

Either make the parameter do something or remove it.

- `make_link`: `link_name`, `length`.
- `hole`: `n_runs`, `show_probe_radius_graph_flag`.
- `map_contours_as_triangles`: `contour_level`.
- `sharpen_with_gompertz_scaling`: `gompertz_factor`.
- `handle_read_draw_probe_dots_unformatted`: `show_clash_gui_flag`.
- `lsq_improve`: `n_res`, `dist_crit`. This one is cheap: `coot::lsq_improve` has `set_n_res_for_frag()` and `set_crit_close()`, and `molecule_class_info_t::lsq_improve()` does not call them. The built-in values are 6 residues and 6 A.
- `regularize_zone_with_score_py` and `regularize_zone_with_score_scm`: `altconf`.
- `add_terminal_residue_using_phi_psi`: `residue_type` (an ALA is always built).
- `make_dots`: `atom_radius_scale`.
- `monomer_lib_3_letter_codes_matching`: `allow_minimal_descriptions_flag`.
- `get_symmetry_bonds_colour`: `imol`.
- `display_pae_from_file_in_a_dialog`: `imol`.
- `servalcat_refine`: `mask_map` (the `--mask` argument is commented out of the command line in `src/cc-interface-map-utils.cc`).

### Stubs: implement or remove

- `alignment_results_scm` and `alignment_results_py`: return False and nothing else. `molecule_class_info_t::align_on_chain()` could supply the result.
- `qq_plot_map_and_model_scm/_py`: the body is inside `HAVE_GOOCANVAS`, which a GTK4 build never defines. The calculation (`coot::util::qq_plot_for_map_over_model()`) is still there; it needs somewhere to plot, or it could return the points.
- `global_phasing_screen`: does its work as a side effect and always returns False.
- `spherical_density_overlap` (Scheme only): returns `#f` and writes `map-fragment-at-origin.map` into the current directory.

### Functions that do nothing useful

- `deviant_geometry()` and `molecule_class_info_t::find_deviant_geometry()`: marked "Deletable for Coot 2" in the source.
- The multi-refine toolbar functions (`toolbar_multi_refine_continue()` and the six others in the same block of `src/c-interface.h`): documented as not functional in GTK4 Coot.
- Setters whose value is never read, each documented as having no effect: `set_cryo_em_refinement`, `set_save_coordinates_in_original_directory`, `set_add_terminal_residue_do_rigid_body_refine`, `set_button_label_for_external_refinement`, `set_add_ccp4i_projects_to_file_dialogs`, `set_file_chooser_overwrite`, `set_filter_fileselection_filenames`, `set_unpathed_backup_file_names`, `set_python_draw_function`, `set_limit_aniso`, `set_show_aniso`, `set_use_stroke_characters`, `set_model_display_radius`.

## Sequence mismatches: chains without a sequence, and repeated assignment

Seen on 2026-11-16 when testing `alignment_mismatches_py()` on
`tutorial-modern.pdb`. Both come from
`molecule_class_info_t::residue_mismatches()` or the sequence assignment, not
from the Python wrapper.

- **A chain with no assigned sequence is reported.** A sequence was assigned to chain A only, with `assign_sequence_from_string()`. The result also had entries for chain B, including "Delete B 94, 95, 96" when the sequence of chain A matched exactly. Find out which sequence chain B is aligned against, and whether that is intended.
- **Assigning a sequence twice doubles the results.** After a second `assign_sequence_from_string()` for chain A there were four alignments, not two, and the chain B deletions were listed twice. The second sequence seems to be added to the first, not to replace it.

## Implement topological_equivalence_chiral_centres()

`topological_equivalence_chiral_centres(residue_type)` (`src/cc-interface.hh`,
body in `src/c-interface-ligands.cc`) is a stub: it returns an empty vector, and
the header doc says so. Nothing calls it.

It is meant to return the atom ids of the chiral centres of a residue type,
found by topological equivalence analysis of the dictionary: the bonds and the
elements only, not the chiral restraints that the dictionary lists. An atom is a
centre if it has four neighbours (or three and an implicit hydrogen) that are
all topologically different.

- **Where it goes.** The comment on the stub says that a better design would put this in `coot::protein_geometry`, with the API function a thin wrapper.
- **What is there to build on.** `protein_geometry::filter_chiral_centres()` already removes dictionary chiral restraints that are not real centres (those with more than one hydrogen). With RDKit (enhanced-ligand builds) the dictionary can be made into an RDKit molecule and the centres found from the canonical atom ranks.
- **A use for it.** Checking the chiral restraints of a dictionary: a centre found here with no chiral restraint, or a restraint on an atom that is not a centre, points to a bad dictionary.

## API documentation bug review: what was not covered

The bugs found while writing the API documentation were fixed on branch
`fix-api-docs-bugs` (November 2026). `api-docs-bug-list.md` has the entry-by-entry
record. This is what that work left undone. The larger items have their own
sections in this file: the post-intermediate-atoms-moved hook, drawing
molecular symmetry, the nomenclature-errors prompt, the Coot 2 API tidy-up, the
sequence mismatches and `topological_equivalence_chiral_centres()`.

### Not decided

- `overlap_ligands`: it graph-matches every residue in `imol_ligand`, moves the best one, and deletes all the other residues of `imol_ligand`. Is the deletion intended?
- `execute_refmac_real`: the Guile branches are guarded by `USE_GUILE_XXX`, which nothing defines, so it always takes the Python path and a build without Python would silently do nothing.

### Fixed but not tested

All of it compiles and installs in the `build-arch-for-claude` sandbox. The tests were headless Python scripts, so these were not run:

- **Anything in the GUI.** In particular the difference-map-peaks dialog with "around the model only" (`only_around_protein_flag`), and the run-state-file dialog.
- **The Scheme functions.** `model_composition_statistics_scm`, `glyco_tree_scm`, `get_residue_alt_confs_scm`, the Guile branches of `run_state_file()` and `run_state_file_maybe()`, and the 6-item spec in `atom_spec_from_scm_expression`.
- **The view matrix of `raster3d()`, `povray()` and `renderman()`.** The numbers in the `.r3d` header look right for a 90 degree rotation, but nothing was rendered. It may need transposing.
- **The FLEV distances.** `fle_view_set_water_dist_max()` and `fle_view_set_h_bond_dist_max()` are now passed to `pli::fle_view_with_rdkit_internal()`, but `fle_view()` opens a dialog and was not run.
- **`probe_available_p()`** with a real `molprobity.probe` (a dummy executable on the `PATH` was used).
- **`servalcat_refine`** and anything else that needs an external program.
- **The chapi build.** `build-chapi-arch-for-claude` was not rebuilt. `pli::fle_view_with_rdkit_internal()` has two new arguments with defaults, so the chapi caller should compile unchanged.
- **Glycan trees on a real glycoprotein.** `print_glyco_tree()` and `glyco_tree_py()` were run only on a two-residue tree (FUC, NAG), so branches and ASN-rooted trees are untested.

### Known and left

- **`generate_local_self_restraints()` in chapi** (`coot::molecule_t`, `api/coot-molecule-refine.cc`) still clears the extra bond restraints on entry. The desktop version no longer does.
- **Self restraints accumulate.** A second call of `generate_local_self_restraints()` for the same chain duplicates the Geman-McClure restraints. Accepted for now: clear them first (`delete_all_extra_restraints()`, or `clear_extra_restraints()` in chapi). The fix would be to remove the restraints between atoms of the selection before adding the new ones.
- **`check_waters_baddies` in AND mode** (`molecule_class_info_t::find_water_baddies_AND()`) still ignores `part_occ_contact_flag`. OR mode uses it.
- **The run-state-file dialog** (`run_state_file_dialog`) previews the Python state file, but in a build with Guile and Python its OK button runs the Scheme one. It is only reached from `src/scm-boot-guile.cc`.
- **The desktop `clear_extra_restraints(int imol)`** (`src/c-interface-kk.cc`) is not declared in any header.
- **`model_composition_statistics_py()`:** the atom overlap score was 0.0 for two of the four test models (those with a residue type that has no dictionary). It is the value that `atom_overlap_score()` gives; find out why. The "good rotamer" cut-off (a probability of 0.3% or more) is a first choice, and C-beta deviation outliers are not counted.
- **`glyco_tree_t::print()`** indents the root the most and the children less, which reads backwards.

## Dialogs shown without a transient parent (macOS full-screen bug)

On macOS, when the main window is full-screen, a dialog that is mapped before
it has been made transient for the main window is itself made full-screen, and
its top row of widgets is off the top of the screen. Seen with the Restraints
Editor and with Go To Atom (both now fixed); Display Manager, which is made
transient before it is shown, was fine.

The fix is per-dialog: call `set_transient_for_main_window(dialog)` before the
`gtk_widget_set_visible(dialog, TRUE)` / `gtk_window_present()`. Where all the
show paths go through one `wrapped_create_*()` function, put the call there (as
was done in `wrapped_create_goto_atom_window()`). Also make sure that the
top-level in the `.ui` file does not have `visible` set to `True` (it is shown
at builder-load time otherwise) - the Restraints Editor was the only one.

The lists below are from a script pass over `ui/*.ui` and `src/*.cc` (14 Nov
2026), with the show call for each hit then read by hand. Line numbers are as
of that date. Only Go To Atom and the Restraints Editor have been confirmed by
running Coot.

### Shown with no transient parent

Navigation / display:

- `global_map_properties_window` - `create-menu-item-actions.cc:4502`
- `display_control_window_glade` - `c-interface.cc:6233` (the scripting path,
  `post_display_control_window()`, only; the button path is fine)

Model building:

- `add_alt_conf_dialog` - `c-interface-build-gui.cc:610`
- `align_and_mutate_dialog` - `create-menu-item-actions.cc:1368`
- `change_chain_id_dialog` - `create-menu-item-actions.cc:945`
- `merge_molecules_dialog` - `c-interface-build-gui.cc:800`
- `replace_fragment_dialog` - `create-menu-item-actions.cc:1100`
- `other_model_tools_dialog` - `create-menu-item-actions.cc:3286`
- `fixed_atom_dialog` - `glade-callbacks-main-window.cc:105`,
  `glade-callbacks.cc:1026,5805`
- `pointer_atom_type_dialog` - `c-interface-build-gui.cc:329`
- `edit_backbone_torsions_dialog` - `graphics-info-gui.cc:2820`
- `chi_angle_help_dialog` - `glade-callbacks.cc:3120`
- `residue_type_chooser_window` - `graphics-info-defines.cc:1557,1590`,
  `c-interface-build-side-chain.cc:661`
- `nucleic_acid_base_chooser_dialog` - `graphics-info-defines.cc:1554`
- `undo_molecule_chooser_dialog` - `graphics-info.cc:4519,4586`
- `recover_coordinates_dialog`, `nothing_to_recover_dialog` -
  `c-interface-build-gui.cc:645,653`

Baton / skeleton:

- `baton_dialog` - `graphics-info-gui.cc:765`
- `baton_mode_make_skeleton_dialog` - `graphics-info.cc:3912`
- `skeleton_dialog` - `graphics-info.cc:3901`

Ligands / waters:

- `find_ligand_many_atoms_dialog` - `c-interface-ligands-gui.cc:680`
- `new_ligands_info_dialog`, `no_new_ligands_info_dialog` -
  `c-interface-ligands.cc:1081,1085`
- `ligand_check_dialog` - `c-interface-ligands-gui.cc:1168`
- `monomer_search_dialog` - `create-menu-item-actions.cc:559`
- `find_waters_dialog` - `c-interface-waters-gui.cc:261`
- `interesting_waters_by_difference_map_check_dialog` -
  `graphics-info-gui.cc:3437`

Validation / refinement:

- `bad_chiral_volumes_dialog`, `no_bad_chiral_volumes_dialog` -
  `graphics-info-gui.cc:3619`, `graphics-info.cc:5919`
- `chiral_restraints_problem_dialog` - `graphics-info.cc:5922`
- `no_restraints_info_dialog` - `graphics-info-modelling.cc:1513`,
  `graphics-info-graphs.cc:1442`
- `download_monomers_dialog` - `graphics-info-modelling.cc:399`
- `no_cif_dictionary_bonds_dialog` - `graphics-info.cc:514`
- `rigid_body_refinement_failed_dialog` - `graphics-info-modelling.cc:3213`
- `add_terminal_residue_finds_none_dialog` - `graphics-info-modelling.cc:3495`
- `nothing_bad_dialog` - `glade-callbacks.cc:3978`
- `simple_refmac_dialog` - `simple-refmac.cc:67`
- `simple_refmac_filechooser_dialog` - `glade-callbacks.cc:6835`

Maps / other:

- `map_partition_by_chain_dialog` - `cc-interface-map-utils.cc:709`
- `hole_dialog` - `create-menu-item-actions.cc:2034`
- `phs_info_box` - `cmtz-interface.cc:678`
- `run_state_file_dialog` - `c-interface.cc:8155`

File choosers:

- `cif_dictionary_filechooser_dialog` - `c-interface-gui.cc:193`
- `export_map_file_chooser_dialog` - `glade-callbacks.cc:6367`
- `save_symmetry_coords_filechooser_dialog` - `graphics-info-defines.cc:1974`
- `save_restraint_chooserdialog` (Restraints Editor "Save") -
  `restraints-editor.cc:1379`

### Not yet traced

The script found no show call for these, so they need checking by hand:
`active_map_window`, `aniso_window`, `centre_atom_label_dialog`,
`check_chiral_volumes_dialog`, `check_waters_diff_map_dialog`,
`close_molecule_dialog`, `delete_item_dialog`, `doc_urls_dialog`, `fps_window`,
`hints_dialog`, `lsq_plane_dialog`, `no_symmetry_warning_dialog`,
`popup_info_window`, `ramachandran_plot_differences_dialog`,
`read_pdb_recentre_dialog`, `rotation_centre_cube_size_window`,
`run_refmac_dialog`, `run_refmac_help_dialog`, `run_script_filechooser_dialog`,
`save_coordinates_filechooser_dialog`, `skeletonization_box_radius_window`,
`skeletonization_level_window`, `smooth_scroll_window`, `stereo_dialog`,
`superpose_nonsense_dialog`, `symmetry_controller_dialog`, `zoom_dialog`,
`single_map_properties_dialog` (own `.ui` file).

### Also unchecked

- The ~42 dialogs that the script passed: it only tested that there is a
  transient call somewhere in the same function as the show call, not that the
  call is for that dialog or that it comes before the show.
- Windows made in code (`gtk_dialog_new()`, `gtk_window_new()`) and dialogs
  made from Python.
