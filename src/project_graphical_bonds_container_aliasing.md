---
name: graphical_bonds_container shallow-copy aliasing
description: graphical_bonds_container uses raw owning pointers with no copy ctor — shallow copies cause use-after-free crashes
type: project
---

`graphical_bonds_container` (coords/graphical-bonds-container.hh) manages multiple raw owning pointer arrays (`bonds_`, `rotamer_markups`, `cis_peptide_markups`, `atom_centres_`, etc.) with manual `new[]`/`delete[]` in `clear_up()`, but has no copy constructor, move constructor, or assignment operator. The compiler-generated defaults do shallow copies.

**Why:** `regularize_object_bonds_box` and `moving_atoms_molecule.bonds_box` are both `graphical_bonds_container` instances. At graphics-info.cc:2459, a shallow copy aliases their raw pointers. When `regularize_object_bonds_box.clear_up()` runs (on any bond regeneration), `moving_atoms_molecule.bonds_box` is left with dangling pointers. The render path (`draw_hud_geometry_bars()`) then reads `moving_atoms_molecule.bonds_box.rotamer_markups` — corrupted `residue_spec_t` strings cause `std::bad_alloc` / SIGABRT.

**How to apply:** Whenever `regularize_object_bonds_box` is updated (clear_up + reassign), `moving_atoms_molecule.bonds_box` must be re-assigned immediately after. This was fixed in `rotate_intermediate_atoms_maybe()`, `rotate_intermediate_atoms_round_screen_z()`, and `rotate_intermediate_atoms_round_screen_x()` (2026-04-24). If adding new code that calls `regularize_object_bonds_box.clear_up()`, always follow with `moving_atoms_molecule.bonds_box = regularize_object_bonds_box;`. The proper long-term fix is to give `graphical_bonds_container` Rule-of-Five semantics or replace the raw arrays with `std::vector`.
