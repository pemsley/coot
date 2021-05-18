/* coot-utils/bonded-pairs.cc
 * 
 * Copyright 2010, 2011, 2012 by The University of Oxford
 * Copyright 2015 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include <algorithm>
#include <stdlib.h> // for abs()

#include "geometry/residue-and-atom-specs.hh"
#include "bonded-pairs.hh"

bool
coot::bonded_pair_container_t::linked_already_p(mmdb::Residue *r1, mmdb::Residue *r2) const {

   bool r = 0;
   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      if (((bonded_residues[i].res_1 == r1) &&
           (bonded_residues[i].res_2 == r2)) ||
          ((bonded_residues[i].res_1 == r2) &&
           (bonded_residues[i].res_2 == r1))) {
         r = 1;
         break;
      }
   }
   return r;
}

bool
coot::bonded_pair_container_t::try_add(const coot::bonded_pair_t &bp) {

   bool found = false;
   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      if ( (bonded_residues[i].res_1 == bp.res_1 &&
            bonded_residues[i].res_2 == bp.res_2) ||
           (bonded_residues[i].res_1 == bp.res_2 &&
            bonded_residues[i].res_2 == bp.res_1) ) {
         found = true;
         break;
      }
   }

   if (! found) {
      bonded_residues.push_back(bp);
   }
   return found;
}

std::ostream&
coot::operator<<(std::ostream &s, coot::bonded_pair_container_t bpc) {

   s << "Bonded Pair Container contains " << bpc.bonded_residues.size() << " bonded residues"
     << "\n";

   for (unsigned int i=0; i<bpc.bonded_residues.size(); i++)
      s << "   " << i << "  [\""
        << bpc[i].link_type << "\" "
        << bpc[i].res_1->GetChainID() << " "    << bpc[i].res_1->GetSeqNum() << " "
        << bpc[i].res_1->GetInsCode() << " to " << bpc[i].res_2->GetChainID() << " "
        << bpc[i].res_2->GetSeqNum() << " "     << bpc[i].res_2->GetInsCode() << "]"
        << "   " << bpc[i]
        << "\n";

   return s; 
}

std::ostream&
coot::operator<<(std::ostream &s, coot::bonded_pair_t bp) {
   s << "[:";
   s << bp.link_type << " ";
   if (bp.res_1)
      s << bp.res_1->GetChainID() <<  " " << bp.res_1->GetSeqNum() << " " << bp.res_1->GetInsCode();
   s << " ";
   if (bp.res_2)
      s << bp.res_2->GetChainID() <<  " " << bp.res_2->GetSeqNum() << " " << bp.res_2->GetInsCode();
   s << "]";
   s << " fixed-flags: " << bp.is_fixed_first << " " << bp.is_fixed_second;
   return s;
}


void
coot::bonded_pair_t::reorder_as_needed() {

   if (res_2->GetSeqNum() < res_1->GetSeqNum()) {
      std::string chain_id_1_i = res_1->GetChainID();
      std::string chain_id_2_i = res_2->GetChainID();
      if (chain_id_1_i == chain_id_2_i) {
         if (res_1->isAminoacid()) {
            if (res_2->isAminoacid()) {
               mmdb::Residue *r = res_1;
               res_1 = res_2;
               res_2 = r;
               std::swap(is_fixed_first, is_fixed_second);
            }
         }
         if (res_1->isNucleotide()) {
            if (res_2->isNucleotide()) {
               mmdb::Residue *r = res_1;
               res_1 = res_2;
               res_2 = r;
               std::swap(is_fixed_first, is_fixed_second);
            }
         }
      }
   }
   // check here for ins code ordering
}

void
coot::bonded_pair_t::apply_chem_mods(const coot::protein_geometry &geom) {

   int imol = protein_geometry::IMOL_ENC_ANY;
   
   if (res_2 && res_2) { 
      try { 
         // apply the mods given the link type

         // get the chem mods for each residue (can throw a runtime
         // error if there is one - (not an error).
         // 
         std::pair<chem_mod, chem_mod> mods = geom.get_chem_mods_for_link(link_type);
         std::string res_1_name = res_1->GetResName();
         std::string res_2_name = res_2->GetResName();
         for (unsigned int i=0; i<mods.first.atom_mods.size(); i++) {
            if (mods.first.atom_mods[i].function == CHEM_MOD_FUNCTION_DELETE) {
               std::string atom_name = mods.first.atom_mods[i].atom_id;
               std::string at_name = geom.atom_id_expand(atom_name, res_1_name, imol);
               delete_atom(res_1, at_name);
            }
         }
         for (unsigned int i=0; i<mods.second.atom_mods.size(); i++) {
            if (mods.second.atom_mods[i].function == CHEM_MOD_FUNCTION_DELETE) {
               std::string atom_name = mods.second.atom_mods[i].atom_id;
               std::string at_name = geom.atom_id_expand(atom_name, res_2_name, imol);
               delete_atom(res_2, at_name);
            }
         }
      }
      catch (const std::runtime_error &rte) {
         // it's OK if we don't find a chem mod for this link
      }
   }
}

void
coot::bonded_pair_container_t::apply_chem_mods(const protein_geometry &geom) {

   std::vector<coot::bonded_pair_t>::iterator it;
   for (it=bonded_residues.begin(); it != bonded_residues.end(); it++) {
      it->apply_chem_mods(geom);
   }
} 


void
coot::bonded_pair_t::delete_atom(mmdb::Residue *res, const std::string &atom_name) {
   
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   bool deleted = false;
   res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (at) {
         std::string at_name(at->name);
         if (at_name == atom_name) {
            delete at;
            at = NULL;
            deleted = true;
         }
      }
   }

   if (deleted)
      res->TrimAtomTable();
}

void
coot::bonded_pair_container_t::reorder() {

   for (unsigned int i=0; i<bonded_residues.size(); i++)
      bonded_residues[i].reorder_as_needed();
}


// remove residue 1-3 bonds if residue 1-2 or 2-3 bonds exist.
void
coot::bonded_pair_container_t::filter() {

   reorder(); // 201501123 - don't reorder. Why? because reorder changes the residue order and
                 // the residues are currently in the correct order (with say, the NAG as res_1
                 // and the ASN as res_2).  Why did I think I needed this?
                 // 20170422 use reorder. But reorder now only reorders protein and nucleotides.
                 // reorder is needed for too-distant-by-residue-numbering comparison filter below.

   std::vector<bonded_pair_t> new_bonded_residues;
   bool debug = false;

   if (debug) {
      std::cout << "DEBUG::: bonded_pair_container_t::filter(): we have these bonded pairs" << std::endl;
      for (unsigned int i=0; i<bonded_residues.size(); i++) {
         const bonded_pair_t &bp_i = bonded_residues[i];
         std::cout << i << "   "
                   << residue_spec_t(bp_i.res_1) << " - to - "
                   << residue_spec_t(bp_i.res_2) << std::endl;
      }
   }

   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      bool keep_this = true;
      const bonded_pair_t &bp_i = bonded_residues[i];
      if (bp_i.res_1 && bp_i.res_2) {
         int resno_delta_i = bp_i.res_2->GetSeqNum() - bp_i.res_1->GetSeqNum();
         if (abs(resno_delta_i) > 1) {
            std::string chain_id_1_i = bp_i.res_1->GetChainID();
            std::string chain_id_2_i = bp_i.res_2->GetChainID();
            if (chain_id_1_i != chain_id_2_i) {

               if (bp_i.link_type == "SS" || bp_i.link_type == "disulf") {
                  // this might be OK then
               } else {
                  keep_this = false; // 20170422 covalently linked residues and carbohydrate must be in
                                     // the same chain now
               }
            } else {
               for (unsigned int j=0; j<bonded_residues.size(); j++) {
                  if (j!=i) {
                     const bonded_pair_t &bp_j = bonded_residues[j];
                     int resno_delta_j = bp_j.res_2->GetSeqNum() - bp_j.res_1->GetSeqNum();

                     // We need to ask if the j'th linked-pair has a residue in
                     // common with the i'th linked pair. If so, which one? If it's
                     // the first one, then we are interested in comparing linked residues
                     // with residue numbers greater than the first, and likewise, if it's
                     // the second residue, then we are interested in comparing the
                     // residue numbers of the first.
                     //
                     // Is the j'th linked-pair more reasonable than the i'th linked-pair?
                     // as judged by the residue number difference being smaller
                     //
                     if (((resno_delta_i > 0) && (bp_i.res_1 == bp_j.res_1)) ||
                         ((resno_delta_i < 0) && (bp_i.res_2 == bp_j.res_2))) {
                     
                        if (abs(resno_delta_j) < abs(resno_delta_i)) {
                           std::string chain_id_1_j = bp_j.res_1->GetChainID();
                           std::string chain_id_2_j = bp_j.res_2->GetChainID();
                           if (chain_id_1_j == chain_id_2_j) {
                              if (chain_id_1_j == chain_id_1_i) {
                                 if (bp_i.link_type == "CIS" || bp_i.link_type == "TRANS" || bp_i.link_type == "PTRANS") {
                                    if (bp_j.link_type == "CIS" || bp_j.link_type == "TRANS" || bp_j.link_type == "PTRANS") {
                                       keep_this = false;
                                       if (debug)
                                          std::cout << ":::::::::::::::::::::: delete bonded pair "
                                                    << residue_spec_t(bp_i.res_1) << " - to = "
                                                    << residue_spec_t(bp_i.res_2) << " because "
                                                    << residue_spec_t(bp_j.res_1) << " - to - "
                                                    << residue_spec_t(bp_j.res_2) << " is closer " << std::endl;
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      if (keep_this)
         new_bonded_residues.push_back(bp_i);
   }

   bonded_residues = new_bonded_residues;

   // We can't do this because closer_exists_p() is not static
   //bonded_residues.erase(std::remove_if(bonded_residues.begin(), bonded_residues.end(), closer_exists_p),
   // bonded_residues.end());
   
}

bool
coot::bonded_pair_container_t::closer_exists_p(const coot::bonded_pair_t &bp_in) const {

   bool e = false;

   if (bp_in.res_1 && bp_in.res_2) { 
      int resno_delta_i = bp_in.res_2->GetSeqNum() - bp_in.res_1->GetSeqNum();
      if (abs(resno_delta_i) > 1) { // needs more sophisticated test for ins-code linked residues
         std::string chain_id_1_i = bp_in.res_1->GetChainID();
         std::string chain_id_2_i = bp_in.res_1->GetChainID();
         if (chain_id_1_i == chain_id_2_i) {
            for (unsigned int j=0; j<bonded_residues.size(); j++) {
               const bonded_pair_t &bp_j = bonded_residues[j];
               int resno_delta_j = bp_j.res_2->GetSeqNum() - bp_j.res_1->GetSeqNum();
               if (abs(resno_delta_j) < resno_delta_i) {
                  std::string chain_id_1_j = bp_j.res_1->GetChainID();
                  std::string chain_id_2_j = bp_j.res_1->GetChainID();
                  if (chain_id_1_j == chain_id_2_j) {
                     if (chain_id_1_i == chain_id_1_j) {
                        e = true;
                     }
                  }
               }
               if (e)
                  break;
            }
         }
      }
   }

   return e;

} 
