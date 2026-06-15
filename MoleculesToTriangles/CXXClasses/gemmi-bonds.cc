/*
 * MoleculesToTriangles/CXXClasses/gemmi-bonds.cc
 *
 * gemmi-native covalent bond determination.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "gemmi-bonds.hh"
#include <gemmi/neighbor.hpp>
#include <gemmi/elem.hpp>

std::vector<coot::m2t::bond_t>
coot::m2t::make_bonds(gemmi::Structure &st, double tolerance) {

   std::vector<bond_t> bonds;
   if (st.models.empty()) return bonds;
   gemmi::Model &model = st.models[0];

   // Assign a stable 0-based serial in iteration order (matches gemmi-selection).
   int idx = 0;
   for (gemmi::Chain &chain : model.chains)
      for (gemmi::Residue &res : chain.residues)
         for (gemmi::Atom &atom : res.atoms)
            atom.serial = idx++;

   // Search radius: largest plausible bond is ~2*max_covalent_r + tolerance.
   // 4.0 A comfortably covers it and keeps grid spacing sensible.
   const double max_search = 4.0;
   gemmi::NeighborSearch ns(model, st.cell, max_search);
   ns.populate(true); // include hydrogens

   for (gemmi::Chain &chain : model.chains) {
      for (gemmi::Residue &res : chain.residues) {
         for (gemmi::Atom &atom : res.atoms) {
            double r1 = atom.element.covalent_r();
            std::vector<gemmi::NeighborSearch::Mark*> marks =
               ns.find_neighbors(atom, 0.1, max_search);
            for (gemmi::NeighborSearch::Mark *m : marks) {
               gemmi::CRA cra = m->to_cra(model);
               if (!cra.atom) continue;
               if (cra.atom->serial <= atom.serial) continue; // each pair once
               double r2 = cra.atom->element.covalent_r();
               double max_bond = r1 + r2 + tolerance;
               double d = atom.pos.dist(cra.atom->pos);
               if (d > 0.4 && d <= max_bond)
                  bonds.push_back({atom.serial, cra.atom->serial});
            }
         }
      }
   }
   return bonds;
}

std::vector<std::vector<int> >
coot::m2t::bonds_adjacency(int n_atoms, const std::vector<bond_t> &bonds) {
   std::vector<std::vector<int> > adj(n_atoms);
   for (const bond_t &b : bonds) {
      if (b.serial_1 >= 0 && b.serial_1 < n_atoms &&
          b.serial_2 >= 0 && b.serial_2 < n_atoms) {
         adj[b.serial_1].push_back(b.serial_2);
         adj[b.serial_2].push_back(b.serial_1);
      }
   }
   return adj;
}
