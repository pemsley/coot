/*
 * MoleculesToTriangles/CXXClasses/gemmi-selection.cc
 *
 * gemmi-native replacement for the mmdb-based CompoundSelection.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "gemmi-selection.hh"

#include <iostream>
#include <algorithm>

namespace coot {
   namespace m2t {

      // The named subset CID strings, copied verbatim from the old
      // CompoundSelection (MMDBSubsetTypePrimitive::subsetTypesArray). gemmi
      // parses these mmdb-style CID strings directly.
      const std::vector<std::pair<std::string, std::string> > &subset_types() {
         static const std::vector<std::pair<std::string, std::string> > t = {
            {"MAIN",         "/*/*/*.*/N,CA,C,O,H"},
            {"SIDE",         "/*/*/*.*/!N,C,O,H"},
            {"WATER",        "/*/*/(WAT,HOH,OH2,H2O)"},
            {"MONOMERS",     "/*/*/(!ALA,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR,WAT,HOH,THP,SEP,TPO,TYP,PTR,OH2,H2O)"},
            {"AMINOACIDS",   "/*/*/(ALA,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR,THP,SEP,TPO,TYP,PTR,MSE)"},
            {"NUCLEICACIDS", "/*/*/(DG,DA,DC,DT,DU,A,G,T,C,U)"},
            {"ALL",          "/*/*/*.*/*:*"}
         };
         return t;
      }

      // Mirrors MMDBSecondaryTypePrimitive::secondaryTypesArray. The integer
      // codes are placeholders (SSE is stubbed - see TODO GEMMI-SSE).
      const std::vector<std::pair<std::string, int> > &secondary_types() {
         static const std::vector<std::pair<std::string, int> > t = {
            {"SSE_None",   0},
            {"SSE_Helix",  1},
            {"SSE_Strand", 2}
         };
         return t;
      }

      static std::string trim(const std::string &s) {
         size_t a = s.find_first_not_of(' ');
         size_t b = s.find_last_not_of(' ');
         if (a == std::string::npos || b == std::string::npos) return "";
         return s.substr(a, (b - a) + 1);
      }
   }
}

// ----- cid_leaf_t -----

coot::m2t::cid_leaf_t::cid_leaf_t(const std::string &cid_in) : cid(cid_in) {
   try {
      sel = gemmi::Selection(cid_in);
   } catch (const std::exception &e) {
      std::cout << "WARNING:: gemmi::Selection failed to parse CID \"" << cid_in
                << "\": " << e.what() << std::endl;
      // leave sel default-constructed but flag it dead via an impossible chain id
      sel = gemmi::Selection();
      sel.chain_ids.all = false;
      sel.chain_ids.inverted = false;
      sel.chain_ids.list = "\x01"; // a chain id that cannot occur
   }
}

bool coot::m2t::cid_leaf_t::matches(const gemmi::CRA &cra) const {
   bool m = sel.matches(cra);
   return invert ? !m : m;
}

// ----- sse_leaf_t -----

coot::m2t::sse_leaf_t::sse_leaf_t(const std::string &name_in) : sse_type(0), name(name_in) {
   const auto &tab = secondary_types();
   for (const auto &p : tab)
      if (p.first == name_in) sse_type = p.second;
}

bool coot::m2t::sse_leaf_t::matches(const gemmi::CRA &cra) const {
   // Secondary structure is stamped on gemmi::Residue::flag ('H'/'E'/'L') by
   // coot::m2t::assign_secondary_structure (DSSP). sse_type follows
   // secondary_types(): SSE_None=0, SSE_Helix=1, SSE_Strand=2.
   bool m = false;
   if (cra.residue) {
      char f = cra.residue->flag;
      if (sse_type == 1)      m = (f == 'H');                 // SSE_Helix
      else if (sse_type == 2) m = (f == 'E');                 // SSE_Strand
      else                    m = (f != 'H' && f != 'E');     // SSE_None
   } else {
      m = (sse_type == 0);
   }
   return invert ? !m : m;
}

// ----- compound_selection_t -----

coot::m2t::compound_selection_t::compound_selection_t(const std::string &selection_text_in) {
   parse(selection_text_in);
}

// Ported from CompoundSelection::setSelectionString - the '& | ! { }' mini-language.
void coot::m2t::compound_selection_t::parse(const std::string &text) {
   selection_text = trim(text);
   invert = false;
   pairs.clear();

   int subClauseLevel = 0;
   std::string subClause;
   bool doInvert = false;
   bool inCID = false;
   combine_t combineRule = combine_t::NEW;

   auto finish_token = [&](const std::string &raw) {
      std::string sub = trim(raw);
      if (sub.empty()) return;
      std::shared_ptr<selection_leaf_t> leaf;
      bool is_secondary = false, is_subset = false;
      for (const auto &p : secondary_types()) if (p.first == sub) is_secondary = true;
      for (const auto &p : subset_types())    if (p.first == sub) is_subset = true;
      if (is_secondary) {
         leaf = std::make_shared<sse_leaf_t>(sub);
      } else if (is_subset) {
         std::string cid;
         for (const auto &p : subset_types()) if (p.first == sub) cid = p.second;
         leaf = std::make_shared<cid_leaf_t>(cid);
      } else {
         leaf = std::make_shared<cid_leaf_t>(sub);
      }
      leaf->invert = doInvert;
      pairs.push_back({combineRule, leaf});
      doInvert = false;
   };

   for (auto it = selection_text.begin(); it != selection_text.end(); ++it) {
      char c = *it;
      if (subClauseLevel > 0) {
         if (c == '}') {
            subClauseLevel--;
            if (subClauseLevel == 0) {
               auto leaf = std::make_shared<compound_selection_t>(subClause);
               leaf->invert = doInvert;
               pairs.push_back({combineRule, leaf});
               doInvert = false;
               subClause.clear();
            } else {
               subClause += c;
            }
         } else if (c == '{') {
            subClauseLevel++;
            subClause += c;
         } else {
            subClause += c;
         }
      } else if (inCID) {
         if (c == '&' || c == '|' || c == '\0' || c == ' ') {
            finish_token(subClause);
            inCID = false;
            if (c == '&') combineRule = combine_t::AND;
            else if (c == '|') combineRule = combine_t::OR;
         } else {
            subClause += c;
         }
      } else if (c == '{') {
         subClause.clear();
         subClauseLevel++;
      } else if (c == '!') {
         doInvert = true;
      } else if (c == '&') {
         combineRule = combine_t::AND;
      } else if (c == '|') {
         combineRule = combine_t::OR;
      } else if (c == ' ' || c == '\t' || c == '\0' || c == '\n' || c == '\r') {
         // skip
      } else {
         subClause = std::string(1, c);
         inCID = true;
      }
   }
   if (!trim(subClause).empty())
      finish_token(subClause);
}

bool coot::m2t::compound_selection_t::matches(const gemmi::CRA &cra) const {
   bool result = false;
   for (const auto &pr : pairs) {
      bool b = pr.second->matches(cra); // leaf applies its own invert
      switch (pr.first) {
         case combine_t::NEW: result = b;             break;
         case combine_t::AND: result = result && b;   break;
         case combine_t::OR:  result = result || b;   break;
      }
   }
   return invert ? !result : result;
}

std::vector<bool> coot::m2t::compound_selection_t::matching_atoms(gemmi::Model &model) const {
   std::vector<bool> result;
   int idx = 0;
   for (gemmi::Chain &chain : model.chains) {
      for (gemmi::Residue &res : chain.residues) {
         for (gemmi::Atom &atom : res.atoms) {
            gemmi::CRA cra{&chain, &res, &atom};
            atom.serial = idx; // stable per-atom index into this bitset
            result.push_back(matches(cra));
            idx++;
         }
      }
   }
   return result;
}

int coot::m2t::compound_selection_t::count_matching_atoms(gemmi::Model &model) const {
   int n = 0;
   for (gemmi::Chain &chain : model.chains)
      for (gemmi::Residue &res : chain.residues)
         for (gemmi::Atom &atom : res.atoms) {
            gemmi::CRA cra{&chain, &res, &atom};
            if (matches(cra)) n++;
         }
   return n;
}

std::string coot::m2t::compound_selection_t::format() const {
   std::string s;
   for (const auto &pr : pairs) {
      switch (pr.first) {
         case combine_t::NEW: s += " ";   break;
         case combine_t::AND: s += " & "; break;
         case combine_t::OR:  s += " | "; break;
      }
      if (pr.second->invert) s += "! ";
      s += pr.second->format();
   }
   return s;
}
