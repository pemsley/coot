/*
 * MoleculesToTriangles/CXXSurface/CXXUtils-gemmi.cc
 *
 * gemmi-native united-atom radius lookup. Reuses the original
 * CXXUtils::unitedAtomRadii table (a static POD array), re-keyed on trimmed
 * names so gemmi's unpadded atom/residue names match.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "CXXUtils-gemmi.hh"
#include "CXXUtils.h"   // reuse CXXUtils::unitedAtomRadii / nAtomRadii (static table)
#include <map>

namespace {
   std::string trim(const std::string &s) {
      size_t a = s.find_first_not_of(' ');
      size_t b = s.find_last_not_of(' ');
      if (a == std::string::npos || b == std::string::npos) return "";
      return s.substr(a, (b - a) + 1);
   }

   // residueName -> (atomName -> radius), keys trimmed. Built once.
   const std::map<std::string, std::map<std::string, float> > &radii_map() {
      static const std::map<std::string, std::map<std::string, float> > m = []() {
         std::map<std::string, std::map<std::string, float> > built;
         for (int i = 0; i < CXXUtils::nAtomRadii; i++) {
            std::string res = trim(CXXUtils::unitedAtomRadii[i].residueName);
            std::string atom = trim(CXXUtils::unitedAtomRadii[i].atomName);
            built[res][atom] = CXXUtils::unitedAtomRadii[i].radius;
         }
         return built;
      }();
      return m;
   }
}

double coot::m2t::get_atom_radius(const std::string &atom_name, const std::string &residue_name) {
   const auto &m = radii_map();
   auto rit = m.find(residue_name);
   if (rit != m.end()) {
      auto ait = rit->second.find(atom_name);
      if (ait != rit->second.end()) return ait->second;
   }
   auto git = m.find("*");
   if (git != m.end()) {
      auto ait = git->second.find(atom_name);
      if (ait != git->second.end()) return ait->second;
   }
   return 1.8;
}

double coot::m2t::get_atom_radius(const gemmi::Atom &atom, const gemmi::Residue &residue) {
   return get_atom_radius(std::string(atom.name), std::string(residue.name));
}
