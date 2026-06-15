/*
 * MoleculesToTriangles/CXXSurface/CXXUtils-gemmi.hh
 *
 * gemmi-native twin of the CXXUtils atom-radius logic (united-atom radii).
 * The original stored per-atom radius via mmdb UDD; here we look it up directly
 * by (atom name, residue name). Lives alongside the original.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXUtils_gemmi_hh
#define CXXUtils_gemmi_hh

#include <string>
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      // United-atom radius by (unpadded) atom name + residue name; 1.8 if no match.
      // Reuses the existing CXXUtils::unitedAtomRadii table (keyed by trimmed names).
      double get_atom_radius(const std::string &atom_name, const std::string &residue_name);
      double get_atom_radius(const gemmi::Atom &atom, const gemmi::Residue &residue);
   }
}

#endif // CXXUtils_gemmi_hh
