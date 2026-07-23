# mmdb-shim architecture

A drop-in replacement for the MMDB2 C++ API (`mmdb::*`) backed by
[gemmi](https://gemmi.readthedocs.io). Coot is built against this
shim instead of libmmdb2: the public headers under `include/mmdb2/` present
MMDB's classes, but every object is a thin wrapper over a live
`gemmi::Structure`. The goal is to remove the MMDB2 dependency while reproducing
MMDB's observable behaviour exactly.


---

## 1. Why a wrapper tree, not a live proxy

A live proxy over `gemmi::Structure` (translate every call on demand) is not
possible, because Coot:

- **writes MMDB public data members directly** — `atom->x`, `->occupancy`,
  `residue->seqNum`, etc. (hundreds of sites);
- relies on **UDData** (`GetUDData`/`PutUDData`, hundreds of call sites) and the
  **handle-based selection engine** — neither of which gemmi has;
- keys ~1000 containers on raw `mmdb::T*` and assumes **pointer identity and
  stability** across edits (`AddAtom`, mid-list insert, `delete atom`).

So the shim **owns an MMDB-shaped wrapper hierarchy** over the gemmi structure
and **reimplements selection + UDData**. gemmi is used only where it does the
real work well (I/O, symmetry, neighbour search, chemistry) — see §6.

Field access is handled by making every MMDB public scalar an **accessor**:
`atom->x` was mechanically rewritten to `atom->x()` — this was a **one-time,
in-place edit of the Coot source tree** (§2), not a header trick. `x()` returns
`gemmi::Atom::pos.x` by reference so both reads and writes flow through to gemmi.
Narrower gemmi types (float `occ`/`b_iso`, `signed char charge`, single-char
`altloc`) use value-get + `set_*`.

---

## 2. The one-time Coot source rewrite (`->x` → `->x()`)

MMDB exposes atom/residue state as **public data members**; gemmi does not, and a
member cannot transparently become a method. So before Coot could build against
the shim, its source tree was rewritten once: every access to a rewritten MMDB
field was turned into a call, e.g.

```
atom->x            ->  atom->x()
atom->x = 1.2      ->  atom->x() = 1.2      // ref-returning accessor: write still works
strcpy(at->name,s) ->  at->SetAtomName(s)   // char[] fields go through a setter
```

**This is a real, committed change to the Coot `.cc`/`.hh` sources** — not
something the headers do at compile time. The shim then defines the matching
accessors (§1).

**The tool.** `mmdb-recon/ast/mmdb_tool.cpp` is a clang **libTooling /
AST-matcher** program (the same Phase-0 recon tool that emitted
`mmdb-recon/ast/mmdb_surface.json`; run with `--rewrite-fields` it becomes a
`RefactoringTool`). Build it with `mmdb-recon/ast/build.sh` (needs
`brew install llvm` — Apple clang lacks the libTooling dev headers). It is
**class-aware**: a config table (`kRules` / `kSetters` in `mmdb_tool.cpp`) lists
exactly which fields on which classes to rewrite —

- `mmdb::Atom`: `x y z occupancy tempFactor charge altLoc serNum`, the ESD
  fields `sigX…sigTemp`, the aniso tensor `u11…u23` (all ref-returning);
  `name`→`GetAtomName`, `element`→`GetElementName`, `residue`→`GetResidue`;
  `strcpy(name…)`→`SetAtomName`, `strcpy(element…)`→`SetElementName`.
- `mmdb::Residue`: `name`→`GetResName`, `seqNum`→`GetSeqNum`,
  `insCode`→`GetInsCode`, `index`→`GetIndex`; `strcpy(name…)`→`SetResName`.

Fields not in the table are left untouched — they surface as compile errors when
building Coot against the shim and are added to the table (or the shim) as they
appear.

**Driver + guardrails** (`rewrite.sh` at the repo root):

- **Rewrite against *real* MMDB, build against the shim.** The tool matches
  `atom->x` as a real-mmdb *field*; if the compile DB pointed at the shim, `->x`
  would already be a method and nothing would match. `rewrite.sh` aborts if the
  `compile_commands.json` contains `mmdb-shim/include`.
- **Dedup across TUs.** A header is parsed once per includer; point-insertions of
  `()` don't conflict in clang's `Replacements`, so without a guard they *stack*
  (`x()()()…`). `firstRewriteAt(file,offset)` records each site so only the first
  wins; `rewrite.sh` also greps the diff for `()()` as a safety net.
- **Do it on a dedicated branch** — it touches hundreds of files; the tool writes
  all edits at the end of the run (after parsing every TU), and it is reversible
  with `git checkout`.

The synthetic fixture in `mmdb-shim/test/rewrite/` exercises the Atom+Residue
rules in the standalone suite.

---

## 3. File layout

Public API headers (what Coot includes) are macro-guarded veneers:

```
include/mmdb2/mmdb_manager.h, mmdb_atom.h, …   (14 headers)
    #ifdef COOT_USE_MMDB_SHIM   -> #include "_shim_impl.hh"   (this shim)
    #else                       -> #include_next <mmdb2/...>  (real MMDB)
```

`#include_next` lets the same tree fall through to real MMDB with the macro off,
for A/B comparison. The implementation is header-only (so the accessor rewrite
and inline hot paths cross TU boundaries) and split into four dependency-ordered
layers, each opening its own `namespace mmdb`:

| Header | Contents |
|---|---|
| `_shim_types.hh` | typedefs, enums, `UDStore` base, forward decls, leaf record classes (LINK/CisPep/Cryst/Helix/Sheet/SymOps/Title), free helpers |
| `_shim_hierarchy.hh` | `Atom` / `Residue` / `Chain` / `Model` class definitions |
| `_shim_manager.hh` | `Manager` (MMDB Root + CoorManager + SelManager) and nested `Selection` / `UDReg` |
| `_shim_inline.hh` | `g()` resolvers, UDData helpers, all out-of-line method bodies, the `detail::` selection matchers, free functions |
| `_shim_impl.hh` | umbrella: includes the four layers in order, then the two subsystems below |

Two subsystems have their own headers, included last (they need a complete
`Atom`/`Residue`):

- `_mmcif_impl.hh` — `mmdb::mmcif::*` over `gemmi::cif`.
- `_graph_impl.hh` — `mmdb::math::{Vertex,Edge,Graph,GraphMatch}` graph matching.

Two heavyweight operations are compiled once in `src/*.cc` (built into
`lib/libmmdbshim.a`) rather than inlined, to keep gemmi's large read/write and
neighbour-search headers out of the ~230 Coot TUs that include
`<mmdb2/mmdb_manager.h>`:

- `src/io.cc` — `Manager::Read*/Write*` via gemmi.
- `src/contacts.cc` — `Manager::SeekContacts` / `SelectNeighbours` via gemmi
  `NeighborSearch`.

---

## 4. The wrapper node and `g()` resolution

Each hierarchy class (`Atom`/`Residue`/`Chain`/`Model`) is a stable node holding:

- `Manager* mgr` — owning manager;
- a **parent pointer** (`res` / `chain` / `model` / — ; null ⇒ *detached*);
- a **cached sibling index** (`ai` / `ri` / `ci` / `mi`);
- a canonical `std::vector<child*>` — this doubles as the **identity cache** and
  as the `PPAtom`/`PPResidue` table MMDB hands back;
- a `gemmi::T _local` used only while detached.

`g()` resolves a wrapper to its live gemmi object **by index through the parent
chain**:

```cpp
gemmi::Atom& Atom::g() { return res->chain->model->mgr->st
                                .models[mi].chains[ci].residues[ri].atoms[ai]; }
```

Index-based (not pointer-based) resolution is the key trick: when a
`std::vector` push reallocates, sibling wrappers stay valid — only the affected
container's indices are patched. This is what gives Coot the pointer stability it
assumes while gemmi's vectors move underneath.

**Detached construction.** Coot's ubiquitous idiom is `new mmdb::Atom; set
fields…; residue->AddAtom(at)`. A wrapper with no parent stores its data in
`_local`; `g()` returns `_local`. `Add*()` copies `_local` into the parent's
gemmi vector, sets the parent pointer + index, and cascades `mgr` down the
subtree that was built while detached.

---

## 5. Memory management

MMDB owns hierarchy nodes with `new`/`delete` and Coot exploits this directly
(`delete atom;`, `delete residue_p;`). The shim matches that contract:

**Ownership.**
- The `Manager` owns the one `gemmi::Structure st` — the single source of truth.
- **Atoms and residues are individually heap-allocated** (`Manager::newAtom` /
  `newRes`) and tracked in `std::set<Atom*> _atom_allocs` / `_res_allocs`, so a
  single `delete atom;` frees exactly one node.
- Chains and models are pooled in `std::deque<Chain>/<Model>` (stable addresses,
  never individually deleted by Coot).
- gemmi-derived metadata records (LINK/CISPEP/HELIX/SHEET/author) live in
  per-`Manager` `std::deque` pools; the `Model` containers hold bare pointers
  into them.

**Deletion of a live atom** (`~Atom`, i.e. Coot's `delete atom;`): detach from
the parent residue's wrapper+gemmi vectors (kept in lockstep), remove from the
manager/model flat atom lists and any selections, reindex trailing siblings, then
drop from `_atom_allocs`. No dangling pointer survives.

**Deferred residue deletion.** MMDB's `DeleteResidue` (and `delete residue_p;`)
is *deferred*: it frees the residue's atoms and **tombstones** the slot (nulls
the wrapper pointer, keeps a gemmi placeholder residue) without changing
`nResidues`, until `FinishStructEdit`/`TrimResidueTable`. Coot depends on this
(e.g. `change_chain_id` iterates to the original count while deleting). The
placeholder keeps surviving siblings' `ri` valid; `_compact_residues()` later
drops the null wrapper slot and its gemmi placeholder in lockstep and reindexes.

**Teardown.** `~Manager` sets `_bulk_free = true` and frees the alloc sets
wholesale; `~Atom`/`~Residue` see the flag and skip all the detach bookkeeping
(memory is going away anyway), avoiding O(n²) teardown.

**Returned arrays.** `SeekContacts` returns `Contact*` allocated with `new[]`;
the caller `delete[]`s it — the MMDB contract. Borrowed C-strings from accessors
(`GetAtomName`, UDData strings) are backed by per-object buffers / `deque<string>`
and are never freed by Coot.

---

## 6. How gemmi is used (the boundaries)

The rule is **prefer gemmi**: never reimplement what gemmi provides. Mapping:

| MMDB surface | gemmi backing |
|---|---|
| PDB / mmCIF read + write | `read_pdb_file` / `read_structure_file` / `write_pdb` / `make_mmcif_document`; `merge_chain_parts()` on read for MMDB one-chain-per-ID parity |
| coordinates, occ, B, aniso, altLoc, element, serial | `gemmi::Atom` fields via `g()` |
| cell / fractionalisation | `gemmi::UnitCell` (`fractionalize`/`orthogonalize`) |
| symmetry (`GetTMatrix`, sym ops) | `SpaceGroup::operations()` composed with the cell frac↔orth transforms |
| neighbour / contact search | `gemmi::NeighborSearch` (grid-accelerated) |
| element data (`getVdWaalsRadius`, `isMetal`) | `gemmi::Element` |
| residue classification (`isAminoacid`/`isSolvent`/`isSugar`/`Get1LetterCode`) | `gemmi::find_tabulated_residue` |
| LINK / CISPEP / HELIX / SHEET / authors | `Structure::{connections,cispeps,helices,sheets}`, `meta.authors` |
| `mmdb::mmcif::*` | `gemmi::cif` (Document/Block/Loop/Table) |

`Manager::_load_metadata()` reshapes gemmi's structure-level metadata into the
MMDB per-`Model` record containers on load.

---

## 7. What the shim implements itself (no gemmi equivalent)

These are hand-written because gemmi has nothing to map onto. They are kept
minimal.

- **UDData slots.** Each object derives from `UDStore` (contiguous
  `int`/`double`/`string` vectors + selection-membership bits). `Manager` holds a
  1-based handle registry (`ud_regs`); registering a handle assigns a per-(type,
  kind) slot; `PutUDData`/`GetUDData` index into the object's vectors.
  Handle `0` means "not registered" — Coot relies on `if (h==0) Register…`.

- **Selection engine.** `Manager` owns a vector of `Selection`s (1-based
  handles). `Select` walks the wrapper tree filtering on model/chain/(seqNum,
  insCode)-range/resname/atom-name/element/altLoc, supporting `SKEY_*` set
  combinators and MMDB's `"!X"` negation and `"*"` wildcard. Membership is
  mirrored in each object's `_inSel` bits so deletion can scrub stale pointers.
  CID strings (`"/1/A/10-20/CA"`) are parsed by a pragmatic tokeniser into that
  same filter call.

- **Pointer-stable wrapper identity.** The whole `vector<child*>` +
  index-locator scheme in §4 exists because gemmi's containers move and have no
  spare identity field; MMDB code cannot tolerate that.

- **Contacts (illustrative).** The common path *does* use gemmi
  `NeighborSearch` — but where gemmi doesn't fit, the shim fills the gap
  in-house. `SeekContacts` with an MMDB `TMatrix` needs contacts against an
  arbitrarily **symmetry-transformed** copy of the second atom set; gemmi's
  search indexes the model's own atoms, not a caller-supplied transformed set. So
  that path (`contacts_transformed` in `src/contacts.cc`) transforms the set and
  runs a small hand-rolled **uniform-grid** neighbour search. The
  `NeighborSearch::Mark`→wrapper mapping goes back through the parallel tree
  (`chains[chain_idx].residues[residue_idx].atoms[atom_idx]`).

- **Graph / subgraph matching** (`_graph_impl.hh`). `mmdb::math::Graph`/
  `GraphMatch` is a branch-and-bound maximum-common-**induced**-subgraph matcher
  (element + bond-type constrained). gemmi has no subgraph isomorphism, so this
  is genuinely shim-owned (design: `../mmdb-graph-matching-for-gemmi.md`).

- **Sequence alignment** (`math::Alignment`). A faithful re-creation of MMDB's
  Needleman-Wunsch with **identity scoring** (match=1/mismatch=0/linear gap).
  gemmi *has* an aligner, but only with substitution-matrix scoring, which would
  change which residues Coot calls mutations vs indels — so matching the MMDB
  baseline requires MMDB's specific scoring, which gemmi does not offer.

- **PDB-column name alignment.** MMDB returns 4-char aligned atom names
  (`" CA "`) and 2-char right-justified upper-case elements (`" C"`, `"NA"`);
  gemmi stores trimmed/mixed-case. `GetAtomName`/`GetElementName` re-pad, and
  matchers trim both sides, because Coot's specs and colour/element tests depend
  on the padded form.

---

## 8. External dependency walls

MMDB is also used by libraries Coot links that are **precompiled against real
MMDB 2.0.22** (different ABI). Passing shim objects across those boundaries is
UB. Handled by shadowing headers (the shim include dir precedes theirs on `-I`):

- **`include/gemmi/mmdb.hpp`** shadows gemmi's own MMDB bridge. The real one
  copies field-by-field through real-mmdb public fields (incompatible with
  accessors). The shim's `copy_to_mmdb`/`copy_from_mmdb` are trivial because
  `Manager` already owns a `gemmi::Structure`.

- **`include/ssm/ssm_align.h`** shadows libssm (structural superposition) with a
  self-contained **no-op `ssm::Align`** — `align()` returns `RC_NoHits` with an
  identity matrix, so callers take their no-superposition path. SSM is being
  removed; its real headers pull the full real-mmdb binary-serialization API.

- **clipper `MMDBAtom_list`**: clipper's prebuilt dylib reads `mmdb::Atom` x/y/z/
  occ/B as *data fields* → garbage on shim atoms. The fix lives in Coot
  (`coot-utils/mmdb-to-clipper-atom-list.hh`): build the `clipper::Atom_list`
  through the shim's *method* accessors. `clipper::Atom_list` is pure clipper with
  no mmdb ABI dependency.

---

## 9. Intentionally inert / mocked

Documented in-source; each is a deliberate decision, not an omission:

- **Secondary-structure assignment** (`CalcSecStructure`) returns the non-OK code
  — gemmi's DSSP has SS prediction disabled upstream, so there is no gemmi SS to
  forward and residues keep `SSE_None`.
- **`MakeBonds`** is a no-op — Coot's only caller recomputes bonds from geometry
  and ignores the MMDB bond table.
- **`Cryst::GetTMatrix`** on a bare `Cryst` is identity (no structure ref);
  `Manager::GetTMatrix` is the live symmetry path.
- `PutPDBString`, `InitMatType`, `SetFlag` are no-ops (behaviour is fixed by the
  gemmi reader/writer).

---

## 10. Build and test

- **Standalone suite** (seconds): `bash mmdb-shim/build.sh` — compiles the
  `test/*.cc` against the shim (core, hierarchy, UDData, selection, contacts,
  I/O, audit, leaf integration, rewriter) plus a macro-off real-MMDB coexistence
  check.
- **Static lib**: `bash mmdb-shim/build-shimlib.sh` → `lib/libmmdbshim.a`
  (rebuild only when `src/*.cc` change; header-only edits don't need it).
- **Full Coot** (the real integration test): `cd ~/lmb/build-coot-mmdb-shim &&
  cmake --build . -j6` (never bare `-j`). Reaches 100% building libcootapi, the
  CLI tools, `test-molecules-container`, and the `coot_headless_api` python
  module.
