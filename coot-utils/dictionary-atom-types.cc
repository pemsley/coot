/* lidia-core/dictionary-atom-types.cc
 * 
 * Copyright 2016 by Medical Research Council
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

#ifndef MAKE_ENHANCED_LIGAND_TOOLS
int main(int argc, char **argv) {return 0;}
#else 
#include "lidia-core/cod-atom-types.hh"

#include <map>
#include <algorithm>

#include "geometry/protein-geometry.hh"
// #include "coot-utils/coot-coord-utils.hh" out of order now
#include "lidia-core/rdkit-interface.hh"
#include "coot-coord-utils.hh"

// ordered so that the table is written in decreasing order
bool string_int_pair_sorter(const std::pair<std::string, unsigned int> &p1,
                            const std::pair<std::string, unsigned int> &p2) {

   return (p2.second < p1.second);
} 

int find_string_in_vector(const std::vector<std::pair<std::string, unsigned int> > &v,
                          const std::string &t) {

   int r = -1;
   for (unsigned int i=0; i<v.size(); i++) {
      if (v[i].first == t) {
         return i;
      } 
   }
   return r;
} 


int main(int argc, char **argv) {

   int r = 0;
   bool debug = false;

   coot::protein_geometry geom;
   geom.set_verbose(false);

   int read_number = 0;
   for (int i=1; i<argc; i++) {
      std::string file_name = argv[i];
      coot::read_refmac_mon_lib_info_t status =
         geom.init_refmac_mon_lib(file_name, read_number);
      read_number++;
   }


   if (debug)
      std::cout << "Examining the atoms in the " << geom.size() << " entries."<< std::endl;

   // key: refmac enerygy type
   // data: vector of COD atom type pairs (type and n occurances)
   std::map<std::string, std::vector<std::pair<std::string, unsigned int> > > atom_map;
   std::map<std::string, std::vector<std::pair<std::string, unsigned int> > >::const_iterator it;

   for (unsigned int i=0; i<geom.size(); i++) {
      const coot::dictionary_residue_restraints_t &r = geom.get_monomer_restraints(i);
      // std::cout << "adding atoms of " << r.residue_info.comp_id << " to atoms map" << std::endl;

      const std::string &comp_id = r.residue_info.comp_id;
      bool idealised_flag = true;
      int imol = 0; // dummy
      mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, imol, idealised_flag);

      if (! mol) {
         std::cout << "Null mol from mol_from_dictionary() for " <<  comp_id << std::endl;
      } else {
         
         mmdb::Residue *residue_p = coot::util::get_first_residue(mol);

         if (! residue_p) {
            // pretty strange
            std::cout << "Null residue from mol from mol_from_dictionary() for "
                      << comp_id << std::endl;
         } else { 

            try { 
               RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, geom);

               cod::atom_types_t t;
               std::vector<cod::atom_type_t> v = t.get_cod_atom_types(rdkm);
               if (v.size() == r.atom_info.size()) {
                  // std::cout << comp_id << " was good" << std::endl;
                  for (unsigned int iat=0; iat<r.atom_info.size(); iat++) {
                     const std::string &key = r.atom_info[iat].type_energy;
                     if (0) 
                        std::cout << comp_id << "  "
                                  << std::setw(4) << r.atom_info[iat].atom_id << "  "
                                  << std::setw(6) << r.atom_info[iat].type_energy << "   " 
                                  << v[iat].level_4 << "\n";

                     it = atom_map.find(key);
                     if (it == atom_map.end()) {
                        std::pair<std::string, unsigned int> p(v[iat].level_4, 1);
                        atom_map[key].push_back(p);
                     } else { 
                        int idx = find_string_in_vector(atom_map[key], v[iat].level_4);
                        if (idx == -1) { // not found
                           std::pair<std::string, unsigned int> p(v[iat].level_4, 1);
                           atom_map[key].push_back(p);
                        } else {
                           atom_map[key][idx].second++;
                        }
                     }
                  }
               } else {
                  std::cout << comp_id << " was not good" << v.size() << " " << r.atom_info.size()
                            << std::endl;
               }
            }
            catch (const std::runtime_error &rte) {
               std::cout << "error:: rte " << rte.what() << " " << comp_id << std::endl;
            } 
            catch (const std::exception &e) {
               std::cout << "error:: exception " << e.what() << " " << comp_id << std::endl;
            } 
         }
      }
   }

   bool output_atom_map = false;
   if (output_atom_map) { 
      std::map<std::string, std::vector<std::pair<std::string, unsigned int> > >::iterator itv;
      for (itv=atom_map.begin(); itv!=atom_map.end(); itv++) {
         std::sort(itv->second.begin(), itv->second.end(), string_int_pair_sorter);
         std::cout << itv->first << "       ";
         for (unsigned int i=0; i<itv->second.size(); i++) { 
            std::cout << "   " << itv->second[i].first  << " " << itv->second[i].second;
         }
         std::cout << "\n";
      }
   }

   // lightweight - this only tells us which refmac atom types there
   // are for a given COD atom type.  That's interesting, but not
   // enough, we need to know the number of ocurrances of that type,
   // so that we can pick the most popular one! Which gives us
   // COD-type -> Refmac-Energy-Type
   //
   // reverse the mapping: the COD type is the key and the data is a
   // vector of refmac energy types (of which should only be one, of
   // course)

   if (0) { 
      std::map<std::string, std::vector<std::pair<std::string, unsigned int> > >::iterator itv;
      std::map<std::string, std::vector<std::string> > cod_map;
      std::map<std::string, std::vector<std::string> >::const_iterator it_cod;
      for (itv=atom_map.begin(); itv!=atom_map.end(); itv++) {
         for (unsigned int i=0; i<itv->second.size(); i++) {
            cod_map[itv->second[i].first].push_back(itv->first);
         }
      }

      for (it_cod=cod_map.begin(); it_cod!=cod_map.end(); it_cod++) {
         if (it_cod->second.size() != 1) {
            std::cout << "strange cod: " << std::left << std::setw(20) << it_cod->first << "     ";
            for (unsigned int i=0; i<it_cod->second.size(); i++) { 
               std::cout << std::setw(5) << it_cod->second[i] << " ";
            }
            std::cout << "\n";
         }
      }
   }

   // reverse atom map:
   // key: COD atom type
   // data: vector of refmac energy type (type and n occurances)
   std::map<std::string, std::vector<std::pair<std::string, unsigned int> > > reverse_atom_map;

   if (debug)
      std::cout << "Reverse atom map " << std::endl;

   for (unsigned int i=0; i<geom.size(); i++) {
      const coot::dictionary_residue_restraints_t &r = geom.get_monomer_restraints(i);

      const std::string &comp_id = r.residue_info.comp_id;
      bool idealised_flag = true;
      int imol = 0; // dummy
      mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, imol, idealised_flag);

      if (! mol) {
         std::cout << "Null mol from mol_from_dictionary() for " <<  comp_id << std::endl;
      } else {
         
         mmdb::Residue *residue_p = coot::util::get_first_residue(mol);

         if (! residue_p) {
            // pretty strange
            std::cout << "Null residue from mol from mol_from_dictionary() for "
                      << comp_id << std::endl;
         } else { 

            try { 
               RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, geom);
               cod::atom_types_t t;
               std::vector<cod::atom_type_t> v = t.get_cod_atom_types(rdkm);
               if (v.size() == r.atom_info.size()) {
                  for (unsigned int iat=0; iat<r.atom_info.size(); iat++) {
                     const std::string &te = r.atom_info[iat].type_energy;
                     const std::string &key = v[iat].level_4;

                     it = reverse_atom_map.find(key);

                     if (it == reverse_atom_map.end()) {
                        std::pair<std::string, unsigned int> p(te, 1);
                        reverse_atom_map[key].push_back(p);
                     } else {
                        // is the type_energy there already?
                        int idx = find_string_in_vector(reverse_atom_map[key], te);
                        if (idx == -1) { // not found
                           std::pair<std::string, unsigned int> p(te, 1);
                           reverse_atom_map[key].push_back(p);
                        } else {
                           reverse_atom_map[key][idx].second++;
                        }
                     }
                  }
               }
            }
            catch (const std::runtime_error &rte) {
               std::cout << "error:: rte " << rte.what() << " " << comp_id << std::endl;
            } 
            catch (const std::exception &e) {
               std::cout << "error:: exception " << e.what() << " " << comp_id << std::endl;
            } 
         }
      }
   }

   bool output_reverse_atom_map = true;
   if (output_reverse_atom_map) { 
      std::map<std::string, std::vector<std::pair<std::string, unsigned int> > >::iterator itv;
      for (itv=reverse_atom_map.begin(); itv!=reverse_atom_map.end(); itv++) {
         std::sort(itv->second.begin(), itv->second.end(), string_int_pair_sorter);
         std::cout << "cod type " << std::left << std::setw(20) << itv->first
                   << "       energy-lib-types: ";
         for (unsigned int i=0; i<itv->second.size(); i++) { 
            std::cout << "   " << itv->second[i].first  << " " << itv->second[i].second;
         }
         std::cout << "\n";
      }
   }
   
   return r;
}

#endif // ENHANCED_LIGAND_TOOLS
