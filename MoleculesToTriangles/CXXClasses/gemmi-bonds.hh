/*
 * MoleculesToTriangles/CXXClasses/gemmi-bonds.hh
 *
 * gemmi-native covalent bond determination (replaces mmdb MakeBonds / GetBonds
 * for the MoleculesToTriangles rendering subsystem).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef MOLECULESTOTRIANGLES_GEMMI_BONDS_HH
#define MOLECULESTOTRIANGLES_GEMMI_BONDS_HH

#include <vector>
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      // A bond between two atoms, identified by their stable serial (the
      // iteration-order index assigned over model 0 - the same convention used
      // by gemmi-selection's matching_atoms()).
      struct bond_t {
         int serial_1;
         int serial_2;
      };

      // Determine covalent bonds for model 0 of `st` by a distance criterion:
      // atoms i,j are bonded if covalent_r(i)+covalent_r(j)+tolerance >= |ri-rj|
      // (and they are not the same atom and not closer than 0.4 A). This
      // reproduces mmdb's distance-based MakeBonds, including backbone peptide
      // bonds. Uses gemmi::NeighborSearch for O(N) spatial lookup.
      //
      // Side effect: assigns atom.serial = iteration index (0-based) over model 0,
      // so the returned serials index a contiguous [0, n_atoms) range.
      std::vector<bond_t> make_bonds(gemmi::Structure &st, double tolerance = 0.4);

      // Adjacency list: adjacency[serial] = serials bonded to that atom.
      std::vector<std::vector<int> > bonds_adjacency(int n_atoms,
                                                     const std::vector<bond_t> &bonds);
   }
}

#endif // MOLECULESTOTRIANGLES_GEMMI_BONDS_HH
