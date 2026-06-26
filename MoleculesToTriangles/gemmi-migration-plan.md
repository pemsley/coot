# MoleculesToTriangles → gemmi migration plan

A focused conversion plan for the `MoleculesToTriangles/` rendering subsystem, as a
self-contained pilot of the larger mmdb2 → gemmi backend migration. This subsystem is a
good test bed precisely because it is *standalone*: nothing inside it reaches into other
Coot code, and it returns no mmdb objects — only `coot::simple_mesh_t` meshes.

Reference trees: gemmi at `~/gemmi/include/gemmi/`, mmdb at `~/compile/mmdb2-2.0.22/`.

---

## 1. Scope by the numbers

~70 files reference mmdb across two libraries:

| Library | Dir | Role |
|---|---|---|
| `libMoleculesToTrianglesCXXClasses` | `CXXClasses/` | molecule wrapper, selections, colour, representations, primitives |
| `libMoleculesToTrianglesCXXSurface` | `CXXSurface/` | molecular / electrostatic surface generation |

mmdb symbol frequency (top): `mmdb::Atom` ×264, `mmdb::Manager` ×132, `mmdb::Residue` ×34,
`mmdb::Chain` ×21, `SELECTION_KEY` ×16, `mmdb::Model` ×13, `SSE_*` ×27, `UDR_ATOM` ×6,
`getVdWaalsRadius` ×5.

mmdb-heaviest files: `CXXClasses/MolecularRepresentation.cpp` (107), `CXXClasses/MyMolecule.cpp`
(72), `CXXClasses/CompoundSelection.cpp` (61), `CXXSurface/CXXSurfaceMaker.cpp` (37),
`CXXClasses/tubes.cc` (37), `CXXSurface/CXXUtils.cpp` (21), `CXXSurface/CXXCreator.cpp` (21).

---

## 2. The external API contract (must keep working)

The rest of Coot calls into M2T from `api/coot-molecule-moltris.cc`,
`src/molecular-mesh-generator.cc`, `src/molecular-mesh-generator-mol-tris.cc`,
`src/molecule-class-info.h`, `docking/intermolecular-energy.cc`. The surface it depends on:

| Entry point | Currently takes | Returns |
|---|---|---|
| `MyMolecule(mmdb::Manager*, int ssFlag)` | mmdb manager | wrapped molecule |
| `MolecularRepresentationInstance::create(MyMolecule, ColorScheme, selStr, style)` | selection **string** | representation → primitives |
| `ColorScheme::colorByElementScheme()` / `…Secondary` / `…BFactor` / `…ChainsScheme…` | — | colour scheme |
| `CompoundSelection(std::string)` | mmdb-style CID string | selection object |
| `CXXCreator(mmdb::Manager*, int selHnd)` | manager **+ mmdb selHnd** | clipper NXmap (electrostatics) |
| `CXXUtils::assignCharge(mmdb::Manager*, int selHnd, CXXChargeTable*)` | manager **+ mmdb selHnd** | writes charges |
| output | — | `coot::simple_mesh_t` (pos/normal/colour + triangles) |

Two of these — `CXXCreator` and `CXXUtils::assignCharge` — take an **mmdb selection handle
from the caller** (`atom_sel.SelectionHandle`). They are the most entangled with the caller and
are sequenced last (the electrostatics path). Everything else takes a selection *string* or no
mmdb at all, so the boundary is clean.

Selection strings in use are mmdb CID, e.g. `"//A/*.*/*:*"`, `"/*/*/*.*/N,CA,C,O,H"`,
`"/*/*/(WAT,HOH)"`.

---

## 3. Strategy

**Gemmi-native internals, mmdb adapter at the boundary.** `MyMolecule` becomes the owner of a
`gemmi::Structure`. The existing `MyMolecule(mmdb::Manager*)` constructor stays as a thin adapter
that calls `gemmi::copy_from_mmdb()` once; a new `MyMolecule(gemmi::Structure)` constructor is
added for the future when callers feed gemmi directly. Every internal consumer (selections,
colour, representations, surface) is rewritten against gemmi. The external signatures above are
preserved; only the *insides* change.

This matches the established migration pattern (on-the-fly `copy_from_mmdb`, see the coot-utils
plan) and keeps the rest of Coot building unchanged throughout.

---

## 4. The three hard subsystems (and their gemmi replacements)

mmdb gives M2T three things "for free" that gemmi must be made to provide. All three have a
gemmi path — confirmed present in the headers — so the migration is feasible end-to-end.

### 4A. Selection engine — built on `gemmi::Selection`

mmdb provides: CID string parsing (`Select(h, STYPE_ATOM, "/*/*/...", SKEY_NEW)`), numbered
selection handles (`int selHnd`), per-atom membership (`atom->isInSelection(h)`), boolean
combination (`SKEY_NEW/AND/OR/CLR`), and inversion. The whole colour/representation machinery is
built on this: `ColorScheme::prepareForMMDB` turns each colour rule into a handle, then
`colorForAtom` asks every atom `isInSelection(ruleHandle)`; representation draw-loops gate each
atom/bond on `isInSelection(selHnd)`.

gemmi already supplies most of this — use it rather than rebuild it:

- **`gemmi::Selection` (`select.hpp`)** parses a CID string in its constructor
  (`Selection(const std::string&)`, CCP4 pdbcur syntax — the same family as mmdb's) and exposes
  `matches(const Atom&)` / `matches(const CRA&)` plus filtered iteration
  (`models()/chains()/residues()/atoms()`). It supports comma-separated **name lists** with `*`
  and `!` inversion for chain ids, **residue names** (the `/mdl/chn/*(res).ic/at[el]:aloc` form,
  i.e. `(WAT,HOH)` and `(!ALA,…)`), atom names, elements and altlocs, **sequence ranges**, and
  occupancy/B-factor inequalities (`q`/`b`).
- **Per-atom membership needs no handle.** `isInSelection(selHnd)` becomes `sel.matches(atom)`. A
  colour rule simply *holds a `gemmi::Selection`* and tests each atom — no `NewSelection`,
  `MakeSelIndex`, or `DeleteSelection` lifecycle to manage.
- **`pymol_select.hpp`** additionally provides a full boolean AST (`AndNode/OrNode/NotNode`) over
  `FlatAtom`, if a PyMOL-style path is ever wanted.

The only genuine residual is that a single `gemmi::Selection` is **conjunctive** (its fields AND
together) — it cannot OR two *different* selections. So `CompoundSelection`'s `&`/`|`/`!`/`{}`
combinator stays, but its job shrinks to combining per-atom booleans:

- `CompoundSelection`'s string tokeniser (`setSelectionString`, the `& | ! {}` parser) is already
  mmdb-independent — keep it.
- Each leaf becomes a `gemmi::Selection`:
  - `MMDBStringPrimitive` (CID) → `gemmi::Selection(selectionString)` **with no translation**.
    Empirically verified (probe over `reference-structures/1e9g.pdb`, 10448 atoms): gemmi parses
    the raw mmdb CID strings M2T uses and returns identical counts — `//A/*.*/*:*` → 4672 (chain A),
    `/*/*/(WAT,HOH,OH2,H2O)` → 1023 (= HOH atoms), `/*/*/*/*`/`""` → 10448. gemmi accepts the `*.*`
    residue wildcard and `(resname,…)` group lists; list inversion `!N,C,O,H` works too. (The only
    forms that fail are *gemmi-native* shapes with an empty residue field, e.g. `///N,CA,C,O,H`.)
  - `MMDBSubsetTypePrimitive` MAIN/SIDE/WATER/AMINOACIDS/NUCLEICACIDS/MONOMERS/ALL → the existing
    CID strings in `CompoundSelection.cpp:388` feed straight into `gemmi::Selection`, unchanged.
  - `MMDBSecondaryTypePrimitive` SSE_Helix/Strand/None → the one leaf *not* expressible in CID;
    test our gemmi per-residue SS (4B) in a small custom leaf predicate.
- The leaf `handleInMMDB` becomes `evaluate(structure) -> std::vector<bool>` (membership keyed by
  the stable `gemmi::Atom::serial` that `copy_from_mmdb` populates); NEW/AND/OR/CLR become bitset
  ops, invert flips the bitset. Where laziness is preferred, hold the `gemmi::Selection`(s) and
  evaluate `matches()` per atom on demand instead of materialising the bitset.

This layer is self-contained and **unit-testable in isolation** against a `gemmi::Structure`
(per-clause atom counts must match the current mmdb path), so it is built and verified first.

### 4B. Secondary structure

mmdb: `model->CalcSecStructure(true)` fills `residue->SSE` (`SSE_Helix/Strand/Bulge/None`); the
HELIX/SHEET header path is `secondary_structure_header_to_residue_sse()` in `MyMolecule.cpp`.
`drawRibbon` switches ribbon width/primitive on `residue->SSE`.

gemmi: `dssp.hpp` `DsspCalculator::calculate_secondary_structure(NeighborSearch&,
Topo::ChainInfo&)` and the free `calculate_dssp(...)`, producing `gemmi::SecondaryStructure`
(Helix/Strand/Loop/…). Heavier than mmdb (needs a `Topo` and a `NeighborSearch`, optionally
hydrogens) but gives better assignments.

**Plan:** add a `gemmi_sse.hh/.cc` helper storing per-residue SS in a side-map keyed by
`{chain, seqid}` (or in `gemmi::Residue::flag`, a spare char). Provide the same three buckets
`drawRibbon` needs (helix / strand / other). Keep the HELIX/SHEET header path by reading
`gemmi::Structure::helices` / `sheets` (gemmi parses these). `USE_HEADER_INFO` /
`CALC_SECONDARY_STRUCTURE` / `DONT_USE` flag semantics are preserved.

### 4C. Bonds and VdW radius

mmdb: `MakeBonds(true)` then `atom->GetBonds(AtomBond*, n)`; `MyMolecule::identifyBonds` also adds
explicit peptide bonds. `mmdb::getVdWaalsRadius(element)` for radii.

gemmi: `vdw_radius(El)` / `Element::vdw_r()` is a direct radius replacement. For connectivity,
`neighbor.hpp` `NeighborSearch` + a covalent-radius distance test reproduces mmdb's distance-based
`MakeBonds`; `bond_idx.hpp` `BondIndex` gives `are_linked`. (mmdb MakeBonds is distance+element
based, so a NeighborSearch + `covalent_radius` sum tolerance is the faithful equivalent and does
not need the monomer library.)

**Plan:** a `gemmi_bonds.hh/.cc` helper building a `std::vector<std::pair<int,int>>` (atom-serial
pairs) via NeighborSearch, plus a `bonded_atoms(serial)` adjacency for the stick/cylinder draw
loops. Peptide-bond augmentation is ported directly (backbone CA/N/C distance test).

---

## 5. Smaller corners

| mmdb usage | site | gemmi replacement |
|---|---|---|
| `AtomStat` centre of selection | `MyMolecule::centreOfSelectionHandle` | mean of selected atom `pos` |
| `SeekContacts` / `SelectNeighbours` | `MolecularRepresentation` H-bonds; `CXXCreator` 35 Å context | `NeighborSearch` / `ContactSearch` (`contact.hpp`) |
| UDD float `"PerAtomRadius"` | `CXXUtils.cpp:57/103/114` etc. | side-table `serial → float`, or `gemmi::Atom::fraction` spare |
| UDD int `"tmp_atom_int"` | `CXXSurface.cpp:509` | use `atom.serial` directly |
| `atom->charge` write (partial charge) | `CXXUtils::assignCharge` | side-table `serial → float` (gemmi `charge` is `signed char` formal charge — cannot hold it) |
| `atom->serNum` ordering | `CXXNewHood`, `CXXCircleNode`, `CXXQADSurface` | `atom.serial` |
| `segID` | `MyMolecule::identifySegments` | gemmi has no segID; use `subchain` or empty |
| `isAminoacid`/`isNucleotide`/`isNucleotideChain` | several | `gemmi::find_tabulated_residue` / `ResidueInfo::is_amino_acid()` / `is_nucleic_acid()` |
| `WritePDBASCII` / `ReadCoorFile` | `MyMolecule` | `gemmi::write_pdb` / `gemmi::read_structure` |

---

## 6. Phase plan (dependency-ordered)

**Phase 0 — leaf replacements (no architecture commitment).** `gemmi_atom_radius` (vdw +
united-atom table port of `CXXUtils::assignUnitedAtomRadius`), residue-type predicates, centre-of-
atoms. Unit-testable.

**Phase 1 — selection engine (4A).** New gemmi-native `CompoundSelection::evaluate(structure)`
returning membership bitsets; CID-dialect translator; MAIN/SIDE/WATER/… predicates. Unit test
against a known PDB (counts per clause must match the current mmdb path exactly).

**Phase 2 — `MyMolecule` core.** Hold `gemmi::Structure`; mmdb-adapter constructor via
`copy_from_mmdb`; port `processCoords`, `identifyBonds` (4C), `identifySegments`,
`identifyDishyBases`, centre/write. SS via 4B.

**Phase 3 — colour.** `ColorScheme` / `ColorRule` / `AtomPropertyRampColorRule` /
`SecondaryColorScheme` / `SolidColorRule`: `colorForAtom(const gemmi::Atom*)`, `prepareForMMDB`
→ `prepare(structure)` building rule bitsets.

**Phase 4 — `MolecularRepresentation`.** All `draw*` methods onto gemmi: spheres, sticks,
cylinders, new-sticks, ribbon (uses 4B SS), calphas, dishy/stick bases, H-bonds, surface chunking.

**Phase 5 — surfaces (`CXXSurface/`).** `CXXSurfaceMaker`, `CXXUtils`, `CXXSurface::assignAtom`,
`CXXQADSurface`, `CXXNewHood`/`CXXCircle*` (use `atom.serial`, `pos`, `element`). Replace the
"PerAtomRadius" UDD with the side-table.

**Phase 6 — electrostatics (most entangled, last).** `CXXCreator` + `CXXUtils::assignCharge`
take an mmdb `selHnd` from the caller (`coot-molecule-moltris.cc`). Either (a) keep these on mmdb
initially behind a converted façade, or (b) change the two call sites to pass a gemmi selection.
Partial charges go to the `serial → float` side-table.

Each phase keeps the build green; the mmdb constructor remains the live entry until Phase 6.

---

## 7. Decisions (settled)

1. **Boundary strategy — "feed gemmi now" (decided).** `MyMolecule` takes a `gemmi::Structure`
   directly; the external callers (`api/coot-molecule-moltris.cc`, `src/molecular-mesh-generator*.cc`)
   change to pass gemmi. No mmdb left at the seam; the mmdb→gemmi conversion (if any) lives in the
   caller, not inside M2T.
2. **Secondary structure — stub for now (decided).** Do **not** implement an SS calculator yet.
   Replace `model->CalcSecStructure(true)` / `residue->SSE` with a clearly-marked stub + TODO note;
   `drawRibbon` and the SSE selection leaf consume that placeholder until the real source (gemmi
   `dssp.hpp` or a CA-geometry calc) is filled in later. **See `// TODO GEMMI-SSE` markers.**
3. **Selection membership** — evaluate `gemmi::Selection::matches(CRA)` per atom during the draw
   iteration (no `selHnd`, no stored bit); materialise a `std::vector<bool>` keyed by `atom.serial`
   only where a whole-model membership set is convenient.
4. **Build** — `CMakeLists-gemmi.txt` only (not Makefile.am). New M2T sources are added to the
   `cootmoleculestotriangles` library; `gemmi::gemmi_cpp` must be added to that target's
   `target_link_libraries` (currently gemmi is linked only at the `cootapi` level). `-DUSE_GEMMI`
   is already global.
5. **Electrostatics path (Phase 6)** — convert the two mmdb-`selHnd` call sites (`CXXCreator`,
   `CXXUtils::assignCharge`) as part of "feed gemmi now"; partial charges → `serial → float`
   side-table.
