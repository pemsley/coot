# Backlog

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
