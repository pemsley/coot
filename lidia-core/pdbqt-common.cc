/*
 * lidia-core/pdbqt-common.cc
 *
 * Shared, RDKit-free PDBQT helpers: line formatting and element-only typing.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option) any
 * later version.
 */

#include <cstdio>
#include "utils/coot-utils.hh"
#include "lidia-core/pdbqt-common.hh"

std::string
coot::pdbqt::atom_line(const writable_atom_t &wa, bool is_het) {
   mmdb::Atom *at = wa.at;
   mmdb::Residue *res = at->residue;
   char line[128];
   std::string alt(at->altLoc);
   std::string ins(res ? res->GetInsCode() : "");
   std::snprintf(line, sizeof(line),
                 "%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f    %6.3f %-2s",
                 is_het ? "HETATM" : "ATOM",
                 wa.serial,
                 at->GetAtomName(),
                 alt.c_str(),
                 res ? res->GetResName() : "UNK",
                 res ? res->GetChainID() : "",
                 res ? res->GetSeqNum() : 0,
                 ins.c_str(),
                 at->x, at->y, at->z,
                 at->occupancy, at->tempFactor,
                 wa.charge,
                 wa.ad_type.c_str());
   return std::string(line);
}

std::string
coot::pdbqt::ad_type_from_element(const std::string &element) {
   std::string e = coot::util::remove_leading_spaces(element);
   if (e == "C")  return "C";
   if (e == "N")  return "N";
   if (e == "O")  return "OA";
   if (e == "S")  return "SA";
   if (e == "H")  return "HD"; // conservative: keep as polar
   if (e == "P")  return "P";
   if (e == "F")  return "F";
   if (e == "CL" || e == "Cl") return "Cl";
   if (e == "BR" || e == "Br") return "Br";
   if (e == "I")  return "I";
   if (e == "ZN" || e == "Zn") return "Zn";
   if (e == "MG" || e == "Mg") return "Mg";
   if (e == "CA" || e == "Ca") return "Ca";
   if (e == "FE" || e == "Fe") return "Fe";
   if (e == "MN" || e == "Mn") return "Mn";
   return coot::util::capitalise(e);
}
