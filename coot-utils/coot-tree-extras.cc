/* coot-utils/coot-tree-extras.cc
 * 
 * Copyright 2009 by The University of Oxford
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

#include <queue>
#include <sstream>
#include <algorithm>
#include <string.h>

#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"
#include "atom-tree.hh"

// constructor can throw an exception
// 
// the constructor should not throw an exception if there are no
// tree in the restraints.  It should instead try the bonds in
// the restraints.  If there are no bonds then it throws an
// exception.
// 
coot::atom_tree_t::atom_tree_t(const coot::dictionary_residue_restraints_t &rest,
                               const coot::minimol::residue &res_in,
                               const std::string &altconf) {

   made_from_minimol_residue_flag = 1;
   mmdb::Residue *residue_p = coot::GetResidue(res_in);
   n_selected_atoms = 0;
   atom_selection = NULL;
   construct_internal(rest, residue_p, altconf);
}


coot::atom_tree_t::atom_tree_t(const dictionary_residue_restraints_t &rest,
                               const std::vector<std::vector<int> > &contact_indices,
                               int base_atom_index,
                               const minimol::residue &res_in,
                               const std::string &altconf) {

   made_from_minimol_residue_flag = 1;
   residue = coot::GetResidue(res_in);
   n_selected_atoms = 0;
   atom_selection = NULL;
   fill_name_map(altconf);
   // don't forget that contact_indices can be generated from the
   // bonds in the dictionary (not only distance-based)
   fill_atom_vertex_vec_using_contacts(contact_indices, base_atom_index);
}

// constructor, given a list of bonds and a base atom index.
//
// Used for multi-residue tree generation.
//
coot::atom_tree_t::atom_tree_t(const std::vector<std::vector<int> > &contact_indices,
                               int base_atom_index,
                               mmdb::Manager *mol,
                               int selection_handle) {

   made_from_minimol_residue_flag = 0;
   residue = NULL; // we can't use (single) residue, we have an
                   // arbitrary atom selection;
   mol->GetSelIndex(selection_handle, atom_selection, n_selected_atoms);
   bool r = fill_atom_vertex_vec_using_contacts_by_atom_selection(contact_indices,
                                                                  atom_selection,
                                                                  n_selected_atoms,
                                                                  base_atom_index);
}


// can throw an exception.
coot::atom_tree_t::atom_tree_t(const coot::dictionary_residue_restraints_t &rest,
                               mmdb::Residue *res,
                               const std::string &altconf) {

   made_from_minimol_residue_flag = 0;
   n_selected_atoms = 0;
   atom_selection = NULL;
   construct_internal(rest, res, altconf);

}

// the constructor, given a list of bonds and a base atom index.
// Used perhaps as the fallback when the above raises an
// exception.
coot::atom_tree_t::atom_tree_t(const std::vector<std::vector<int> > &contact_indices,
                               int base_atom_index,
                               mmdb::Residue *res,
                               const std::string &altconf) {

   made_from_minimol_residue_flag = 0;
   n_selected_atoms = 0;
   atom_selection = NULL;
   if (! res) {
      std::string mess = "null residue in alternate atom_tree_t constructor";
      throw std::runtime_error(mess);
   } else {
      residue = res;
      fill_name_map(altconf);
      fill_atom_vertex_vec_using_contacts(contact_indices, base_atom_index);
   }

}



void coot::atom_tree_t::construct_internal(const coot::dictionary_residue_restraints_t &rest,
                                           mmdb::Residue *res,
                                           const std::string &altconf) {

   auto residue_has_deuterium_atoms = [] (mmdb::Residue *residue) {
                                         int nResidueAtoms = residue->GetNumberOfAtoms();
                                         bool has_deuterium_atoms = false;
                                         for (int iat=0; iat<nResidueAtoms; iat++) {
                                            mmdb::Atom *atom_p = residue->GetAtom(iat);
                                            if (! atom_p->isTer()) {
                                               std::string atom_ele(atom_p->element);
                                               if (atom_ele == " D") {
                                                  has_deuterium_atoms = true;
                                                  break;
                                               }
                                            }
                                         }
                                         return has_deuterium_atoms;
                                      };

   if (! res) {
      std::string mess = "Null residue in atom tree constructor";
      throw std::runtime_error(mess);
   }
   residue = res; // class data now, used in rotate_about() - so that
                  // we don't have to pass the mmdb::Residue * again.

   if (rest.tree.size() == 0) {
      std::string mess = "atom_tree_t()::construct_internal(): No tree in restraints for " + rest.comp_id();
      throw std::runtime_error(mess);
   }

   // Fill bonds
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   res->GetAtomTable(residue_atoms, n_residue_atoms);
   // fill the bonds
   bool has_deuterium_atoms = residue_has_deuterium_atoms(residue);
   for (unsigned int i=0; i<rest.bond_restraint.size(); i++) {
      int idx1 = -1;
      int idx2 = -1;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         std::string atom_name = residue_atoms[iat]->name;
         std::string atom_altl = residue_atoms[iat]->altLoc;
         //          std::cout << "comparing :" << atom_name << ": with :" << rest.bond_restraint[i].atom_id_1()
         //                    << ":" << std::endl;
         if (atom_name == rest.bond_restraint[i].atom_id_1_4c())
            if (atom_altl == "" || atom_altl == altconf)
               idx1 = iat;
         if (atom_name == rest.bond_restraint[i].atom_id_2_4c())
            if (atom_altl == "" || atom_altl == altconf)
               idx2 = iat;
         // OK, we have found them, no need to go on searching.
         if ((idx1 != -1) && (idx2 != -1))
            break;
      }

      if (has_deuterium_atoms) {

         // std::cout << ":::::: in construct_internal() has_deuterium_atoms" << std::endl;

         // same again with dictionary atom name changes
         for (int iat=0; iat<n_residue_atoms; iat++) {
            std::string atom_name = residue_atoms[iat]->name;
            std::string atom_altl = residue_atoms[iat]->altLoc;
            //          std::cout << "comparing :" << atom_name << ": with :" << rest.bond_restraint[i].atom_id_1()
            //                    << ":" << std::endl;
            std::string bond_restraint_atom_name_1 = rest.bond_restraint[i].atom_id_1_4c();
            std::string bond_restraint_atom_name_2 = rest.bond_restraint[i].atom_id_2_4c();
            if (bond_restraint_atom_name_1[0] == 'H') bond_restraint_atom_name_1[0] = 'D';
            if (bond_restraint_atom_name_1[1] == 'H') bond_restraint_atom_name_1[1] = 'D';
            if (bond_restraint_atom_name_2[0] == 'H') bond_restraint_atom_name_2[0] = 'D';
            if (bond_restraint_atom_name_2[1] == 'H') bond_restraint_atom_name_2[1] = 'D';
            if (atom_name == bond_restraint_atom_name_1)
               if (atom_altl == "" || atom_altl == altconf)
                  idx1 = iat;
            if (atom_name == bond_restraint_atom_name_2)
               if (atom_altl == "" || atom_altl == altconf)
                  idx2 = iat;
            // OK, we have found them, no need to go on searching.
            if ((idx1 != -1) && (idx2 != -1))
               break;
         }
      }

      // std::cout << "bond restraint " << i << " " << rest.bond_restraint[i] << " " << idx1 << " " << idx2 << std::endl;

      if ((idx1 != -1) && (idx2 != -1)) {
         bonds.push_back(std::pair<int,int>(idx1, idx2));
      }
   }

   fill_name_map(altconf);

   bool success_vertex = fill_atom_vertex_vec(rest, res, altconf);
   if (! success_vertex) {
      std::string mess = "Failed to fill atom vector from cif atom tree - bad tree?";
      throw std::runtime_error(mess);
   }

   bool success_torsion = fill_torsions(rest, res, altconf);

   // if you want print out the atom tree, do it in fill_atom_vertex_vec()
}


void
coot::atom_tree_t::fill_name_map(const std::string &altconf) {

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);

   // atom-name -> index, now class variable
   // std::map<std::string, int, std::less<std::string> > name_to_index;
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      std::string atom_altl = residue_atoms[iat]->altLoc;
      if (false)
         std::cout << "debug:: in fill_name_map(): comparing altconf of this atom :" << atom_altl
                   << ": to (passed arg) :" << altconf << ": or blank" << std::endl; 
      if (atom_altl == "" || atom_altl == altconf) {
         // std::cout << "assigning " << atom_name << " map_index_t of " << iat << std::endl;
         name_to_index[atom_name] = map_index_t(iat);
      }
   }
}



bool
coot::atom_tree_t::fill_torsions(const coot::dictionary_residue_restraints_t &rest, mmdb::Residue *res,
                                 const std::string &altconf) {


   bool r = 0;
   int n_torsions_inserted = 0;
   if (rest.torsion_restraint.size() > 0) {
      std::vector<coot::atom_index_quad> quads;
      for (unsigned int itr=0; itr<rest.torsion_restraint.size(); itr++) { 
         coot::dict_torsion_restraint_t tr = rest.torsion_restraint[itr];
         std::pair<bool,coot::atom_index_quad> quad_pair = get_atom_index_quad(tr, res, altconf);
         if (quad_pair.first)
            quads.push_back(quad_pair.second);
      }

      // Now we have a set of atom index quads, put them in the
      // atom_vertex_vec, at the position of the second atom in the
      // torsion (the position of the first atom in the rotation
      // vector).

      // std::cout << " debug:: " << quads.size() << " quads" << std::endl;
      for (unsigned int iquad=0; iquad<quads.size(); iquad++) {
         bool inserted = 0;
         int n_atom_vertex = atom_vertex_vec.size();
         for (int iv=0; iv<n_atom_vertex; iv++) {
            if (iv == quads[iquad].index2) {
               for (unsigned int ifo=0; ifo<atom_vertex_vec[iv].forward.size(); ifo++) { 
                  if (atom_vertex_vec[iv].forward[ifo] == quads[iquad].index3) {

                     // now check that the forward atom of this
                     // forward atom is index4
                     int this_forward = atom_vertex_vec[iv].forward[ifo];
                     for (unsigned int ifo2=0; ifo2<atom_vertex_vec[this_forward].forward.size(); ifo2++) {
                        if (atom_vertex_vec[this_forward].forward[ifo2] == quads[iquad].index4) {
                           atom_vertex_vec[iv].torsion_quad.first = 1;
                           atom_vertex_vec[iv].torsion_quad.second = quads[iquad];
                           if (0) {
                              std::cout << "            DEBUG:: inserting torsion " << iquad
                                        << " into vertex " << iv << std::endl;
                           }
                           r = 1;
                           inserted = 1;
                           n_torsions_inserted++;
                        }
                     } 
                  }
               }
            }
            if (inserted)
               break;
         }
      } 
   }
   // std::cout << "DEBUG:: inserted " << n_torsions_inserted << " torsions" << std::endl;
   return r;
}

// Don't forget that contact_indices can be generated from the bonds
// in the dictionary (not only distance-based).
// 
bool
coot::atom_tree_t::fill_atom_vertex_vec_using_contacts(const std::vector<std::vector<int> > &contact_indices,
                                                       int base_atom_index) {


   bool r=0; // return 1 on successful fill

   mmdb::Atom **residue_atoms = nullptr;
   int n_residue_atoms = 0;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   atom_vertex_vec.resize(n_residue_atoms);

   r = fill_atom_vertex_vec_using_contacts_by_atom_selection(contact_indices, residue_atoms,
                                                             n_residue_atoms, base_atom_index);
   return r;
}

bool
coot::atom_tree_t::fill_atom_vertex_vec_using_contacts_by_atom_selection(const std::vector<std::vector<int> > &contact_indices,
                                                                         mmdb::PPAtom residue_atoms,
                                                                         int n_atoms,
                                                                         int base_atom_index) {

   bool debug = false;
   if (debug)
      std::cout << ":::::::::::: fill_atom_vertex_vec_using_contacts_by_atom_selection() --- start -- with n_atoms "
                << n_atoms << std::endl;

   bool r = false;
   coot::atom_vertex av;
   atom_vertex_vec.resize(n_atoms);
   av.connection_type = coot::atom_vertex::START;
   atom_vertex_vec[base_atom_index] = av;

   // fail to set up
   if (contact_indices.size() == 0)
      return 0;

   if (debug) {
      std::cout << "debug:: starting fill_atom_vertex_vec_using_contacts_by_atom_selection(): --- start --- " << std::endl;
      std::cout << " debug:: =========== contact indices ======= " << std::endl;
      for (unsigned int ic1=0; ic1<contact_indices.size(); ic1++) {
         std::cout << " index " << ic1 << " : ";
         for (unsigned int ic2=0; ic2<contact_indices[ic1].size(); ic2++)
            std::cout << contact_indices[ic1][ic2] << " ";
         std::cout << std::endl;
      }
   }


   std::queue<int> q;
   q.push(base_atom_index);
   std::vector<int> done; // things that were in the queue that are
                          // not done (on poping). So that we don't
                          // loop endlessly doing atom indices that we
                          // have already considered.

   while (q.size()) {
      int this_base_atom = q.front();
      // now what are the forward atoms of av?
      std::vector<int> av_contacts = contact_indices[this_base_atom]; // size check above

      for (unsigned int iav=0; iav<av_contacts.size(); iav++) {

         int i_forward = av_contacts[iav];

         // Check that this forward atom is not already in the forward
         // atoms of this_base_atom
         bool ifound_forward = 0;
         for (unsigned int ifo=0; ifo<atom_vertex_vec[this_base_atom].forward.size(); ifo++) {
            if (atom_vertex_vec[this_base_atom].forward[ifo] == av_contacts[iav]) {
               ifound_forward = 1;
               break;
            }
         }

         if (! ifound_forward) {
            // Add the contact as a forward atom of this_base_atom but
            // only if the forward atom does not have a forward atom
            // which is this_base_atom (and this keeps the tree going in
            // up and not back again).
            bool ifound_forward_forward = 0;
            for (unsigned int ifo=0; ifo<atom_vertex_vec[i_forward].forward.size(); ifo++) {
               if (atom_vertex_vec[i_forward].forward[ifo] == this_base_atom) {
                  ifound_forward_forward = 1;
                  break;
               }
            }
            if (!ifound_forward_forward) {
               atom_vertex_vec[this_base_atom].forward.push_back(av_contacts[iav]);
            }
         }

         // add the forward atoms to the queue if they are not already
         // in the done list.
         bool in_done = 0;
         for (unsigned int idone=0; idone<done.size(); idone++) {
            if (done[idone] == av_contacts[iav]) {
               in_done = 1;
               break;
            }
         }
         if (!in_done)
            q.push(av_contacts[iav]);
      }

      // if the forward atoms of this_atom_index do not have a back
      // atom, add one (but not, of course, if forward atom is the
      // start point of the tree).
      for (unsigned int iav=0; iav<av_contacts.size(); iav++) {
         if (atom_vertex_vec[av_contacts[iav]].backward.size() == 0) {
            if (atom_vertex_vec[av_contacts[iav]].connection_type != coot::atom_vertex::START)
               atom_vertex_vec[av_contacts[iav]].backward.push_back(this_base_atom);
         }
      }
      q.pop();
      done.push_back(this_base_atom);
      r = 1; // return success
   }

   // print out the name_to_index map
   if (debug) {
      std::cout << "==== atom indexes ===\n";
      for(std::map<std::string, coot::map_index_t>::const_iterator it = name_to_index.begin();
          it != name_to_index.end(); ++it)
         std::cout << "Atom :" << it->first << ": Index: " << it->second.index() << '\n';
   }

   // print out the atom tree
   if (debug) {
      std::cout << "debug:: ==== atom_vertex_vec === "
                << "from fill_atom_vertex_vec_using_contacts_by_atom_selection "
                << std::endl;
      for (unsigned int iv=0; iv<atom_vertex_vec.size(); iv++) {
         std::cout << "   atom_vertex_vec[" << iv << "] forward atom ";
         for (unsigned int ifo=0; ifo<atom_vertex_vec[iv].forward.size(); ifo++)
            std::cout << atom_vertex_vec[iv].forward[ifo] << " ";
         std::cout << std::endl;
      }
   }

   return r;
}

// This used to throw an exception.  But now it does not because the
// catcher was empty (we don't want to know about unfilled hydrogen
// torsions).  And an empty catch is bad.  So now we return a pair,
// the first of which defines whether the torsion is good or not.
//
std::pair<bool,coot::atom_index_quad>
coot::atom_tree_t::get_atom_index_quad(const coot::dict_torsion_restraint_t &tr,
                                       mmdb::Residue *res,
                                       const std::string &altconf) const {
   bool success = 0;
   coot::atom_index_quad quad(-1,-1,-1,-1);
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(   residue_atoms[iat]->name);
      std::string atom_altconf(residue_atoms[iat]->altLoc);
      if (atom_name == tr.atom_id_1_4c())
         if (atom_altconf == "" || atom_altconf == altconf)
            quad.index1 = iat;
      if (atom_name == tr.atom_id_2_4c())
         if (atom_altconf == "" || atom_altconf == altconf)
            quad.index2 = iat;
      if (atom_name == tr.atom_id_3_4c())
         if (atom_altconf == "" || atom_altconf == altconf)
            quad.index3 = iat;
      if (atom_name == tr.atom_id_4_4c())
         if (atom_altconf == "" || atom_altconf == altconf)
            quad.index4 = iat;
   }
   if ((quad.index1 == -1) || (quad.index2 == -1) ||
       (quad.index3 == -1) || (quad.index4 == -1)) {
      success = 0;
   } else {
      success = 1;
   }
   return std::pair<bool,coot::atom_index_quad> (success,quad);
}


bool
coot::atom_tree_t::fill_atom_vertex_vec(const coot::dictionary_residue_restraints_t &rest,
                                        mmdb::Residue *res,
                                        const std::string &altconf) {

   bool debug = 0;

   // Note that we add an extra check on the forward and back atom
   // indices before adding them.  This is a sanity check - if we come
   // here with a big tree but a small residue, then the forward and
   // back atom indices could go beyond the limits of the size of the
   // atom_vertex_vec.  Which would be badness.  Hmm.. I am not sure
   // about this now. The indices are into the residue table, are they
   // not?  And the residue table corresponds here to residue_atoms
   // and n_residue_atoms. Hmm.. maybe that was not the problem.

   bool retval = 0; // fail initially.

   bool found_start = 0;
   int rest_tree_start_index = -1;
   for (unsigned int i=0; i<rest.tree.size(); i++) {
      if (rest.tree[i].connect_type == "START") { 
         found_start = 1;
         rest_tree_start_index = i;
         break;
      }
   }

   // print out the name_to_index map
   if (debug) {
      std::cout << "==== atom indexes ===\n";
      for(std::map<std::string, coot::map_index_t>::const_iterator it = name_to_index.begin();
          it != name_to_index.end(); ++it)
         std::cout << "Atom :" << it->first << ": Index: " << it->second.index() << '\n';
   }
 

   if (found_start) {

      // atom vertex is based on residue atoms (not dict).
      int n_residue_atoms = res->GetNumberOfAtoms();
      if (debug)
         std::cout << "in fill_atom_vertex_vec(), residue has " << n_residue_atoms
                   << " atoms " << std::endl;
      atom_vertex_vec.resize(n_residue_atoms);

      // fill atom_vertex_vec, converting dictionary atom names to
      // residue indices.
      // 
      for (unsigned int itree=0; itree<rest.tree.size(); itree++) {
         coot::map_index_t atom_id_index = name_to_index[rest.tree[itree].atom_id];
         if (atom_id_index.is_assigned()) {
            retval = 1;
            int idx = atom_id_index.index();
            coot::map_index_t atom_back_index =
               name_to_index[rest.tree[itree].atom_back];
            if (rest.tree[itree].atom_back != "") { 
               // if connect_type is START then the back atom is n/a -> UNASSIGNED
               if (atom_back_index.is_assigned()) { 
                  if (atom_back_index.index() < n_residue_atoms) { 
                     // there should be only one back atom for each vertex
                     atom_vertex_vec[idx].backward.push_back(atom_back_index.index());
                     // 
                     // make a synthetic forward atom
                     add_unique_forward_atom(atom_back_index.index(), idx);
                  }
               }
            }
               
            coot::map_index_t atom_forward_index =
               name_to_index[rest.tree[itree].atom_forward];

            if (atom_forward_index.is_assigned()) {
               if (atom_forward_index.index() < n_residue_atoms)
                  add_unique_forward_atom(idx, atom_forward_index.index());
            }

            if (debug) { 
               std::cout << "   itree " << itree << " :"
                         << rest.tree[itree].atom_id << ": :"
                      << rest.tree[itree].atom_back << ": :"
                         << rest.tree[itree].atom_forward << ": -> "
                         << atom_id_index.index() << " "
                         << atom_back_index.index() << " "
                         << atom_forward_index.index()
                         << std::endl;
            }
            
            atom_vertex_vec[idx].connection_type = coot::atom_vertex::STANDARD;
            if (rest.tree[itree].connect_type == "START")
               atom_vertex_vec[idx].connection_type = coot::atom_vertex::START;
            if (rest.tree[itree].connect_type == "END")
               atom_vertex_vec[idx].connection_type = coot::atom_vertex::END;
         } else {
            if (debug)
               std::cout << "DEBUG:: in fill_atom_vertex_vec() "
                         << " no index assignment for "
                         << rest.tree[itree].atom_id 
                         << std::endl;
         } 
      }
   } else {
      if (debug)
         std::cout << "DEBUG:: in fill_atom_vertex_vec() no found start"
                   << std::endl;
   } 

   // print out the atom tree
   if (debug) {
      std::cout << "debug:: ==== atom_vertex_vec (atom_tree) === from fill_atom_vertex_vec "
                << "======== size " <<  atom_vertex_vec.size() << std::endl;
      for (unsigned int iv=0; iv<atom_vertex_vec.size(); iv++) {
         std::cout << "   atom_vertex_vec[" << iv << "] forward atoms ("
                   << atom_vertex_vec[iv].forward.size() << ") ";
         for (unsigned int ifo=0; ifo<atom_vertex_vec[iv].forward.size(); ifo++) 
            std::cout << atom_vertex_vec[iv].forward[ifo] << " ";
         std::cout << std::endl;
      }
   }
   if (debug)
      std::cout << "fill_atom_vertex_vec() returns " << retval << std::endl;
   return retval;
}




// Add forward_atom_index as a forward atom of this_index - but only
// if forward_atom_index is not already a forward atom of this_index.
void
coot::atom_tree_t::add_unique_forward_atom(int this_index, int forward_atom_index) { 

   bool ifound = 0;
   for (unsigned int ifo=0; ifo<atom_vertex_vec[this_index].forward.size(); ifo++) {
      if (atom_vertex_vec[this_index].forward[ifo] == forward_atom_index) {
         ifound = 1;
         break;
      }
   }

   std::pair<int, std::vector<coot::map_index_t> >
      forward_atoms_of_forward_atom_index_pair = get_forward_atoms(map_index_t(forward_atom_index), map_index_t(forward_atom_index));

   if (0) 
      std::cout << " in add_unique_forward_atom get_forward_atoms called "
                << forward_atoms_of_forward_atom_index_pair.first 
                << " times - indices size: "
                << forward_atoms_of_forward_atom_index_pair.second.size() << std::endl;
   
   if (0) { 
      std::cout << "debug:: forward atoms of " << forward_atom_index << ":"; 
      for (unsigned int i=0; i<forward_atoms_of_forward_atom_index_pair.second.size(); i++)
         std::cout << " " << forward_atoms_of_forward_atom_index_pair.second[i].index();
      std::cout << std::endl;
   }

   for (unsigned int i=0; i<forward_atoms_of_forward_atom_index_pair.second.size(); i++)
      if (this_index == forward_atoms_of_forward_atom_index_pair.second[i].index()) {
         ifound = 1;
         if (0) { 
            std::cout << " reject this attempt to add synthetic forward atom because "
                      << this_index << " is a forward atom of " << forward_atom_index
                      << std::endl;
         }
      }

   if (ifound == 0) {
//       std::cout << "So, new forward: to this_index " << this_index << " added forward "
//                 << forward_atom_index << std::endl;
      atom_vertex_vec[this_index].forward.push_back(forward_atom_index);
   } 
}

std::pair<unsigned int, unsigned int>
coot::atom_tree_t::fragment_sizes(const std::string &atom1,
                                  const std::string &atom2,
                                  bool reversed_flag) {

   coot::map_index_t index2 = name_to_index[atom1];
   coot::map_index_t index3 = name_to_index[atom2];
   std::vector<map_index_t> m = get_unique_moving_atom_indices(atom1, atom2, reversed_flag);
   std::vector<map_index_t> c = complementary_indices(m, index2, index3);

   return std::pair<unsigned int, unsigned int>(m.size(), c.size());
} 

std::vector<coot::map_index_t>
coot::atom_tree_t::get_unique_moving_atom_indices(const std::string &atom1,
                                                  const std::string &atom2,
                                                  bool reversed_flag) {

   std::vector<coot::map_index_t> unique_moving_atoms;

   bool debug = false;
   // OK, so when the user clicks atom2 then atom1 (as the middle 2
   // atoms), then they implictly want the fragment revsersed,
   // relative to the atom order internally.  In that case, set
   // internal_reversed, and only if it is not set and the
   // reversed_flag flag *is* set, do the reversal (if both are set
   // they cancel each other out).  Now we do not pre-reverse the
   // indices in the calling function. The reverse is done here.
   bool internal_reversed = false;

   if (debug)
      std::cout << "rotate_about() " << atom1 << " " << atom2 << std::endl;

   coot::map_index_t index2 = name_to_index[atom1];
   coot::map_index_t index3 = name_to_index[atom2];

   if ((atom_vertex_vec[index2.index()].forward.size() == 0) && 
       (atom_vertex_vec[index3.index()].forward.size() == 0)) {
      std::string s = "Neither index2 ";
      s += coot::util::int_to_string(index2.index());
      s += " nor index3 ";
      s += coot::util::int_to_string(index3.index());
      s += " has forward atoms!";
      throw std::runtime_error(s);
   } 
       
   if (index2.is_assigned()) { 
      if (index3.is_assigned()) {

         // is index3 in the forward atoms of index2?
         // if not, then swap and test again.
         //
         bool index3_is_forward = 0;

         for (unsigned int ifo=0; ifo<atom_vertex_vec[index2.index()].forward.size(); ifo++) {
            if (atom_vertex_vec[index2.index()].forward[ifo] == index3.index()) {
               index3_is_forward = 1;
               break;
            }
         }

         if (! index3_is_forward) {
            // perhaps index2 is the forward atom of index3?
            bool index2_is_forward = 0;
            for (unsigned int ifo=0; ifo<atom_vertex_vec[index3.index()].forward.size(); ifo++) {
               if (atom_vertex_vec[index3.index()].forward[ifo] == index2.index()) {
                  index2_is_forward = 1;
                  break;
               }
            }

            // if index2 *was* the forward atom of index3, then swap
            // around index2 and 3 into "standard" order.
            // 
            if (index2_is_forward) {
               if (debug)
                  std::cout << "    index2 was the forward atom of index3 - swapping "
                            << "and setting internal_reversed" << std::endl;
               std::swap(index2, index3);
               index3_is_forward = 1;
               internal_reversed = 1;
            }
         }

         // OK, try again
         if (index3_is_forward) {

            std::pair<int, std::vector<coot::map_index_t> > moving_atom_indices_pair =
               get_forward_atoms(index3, index3);
            std::vector<coot::map_index_t> moving_atom_indices =
               moving_atom_indices_pair.second;
            if (debug) 
               std::cout << " in rotate_about() get_forward_atoms() called "
                         << moving_atom_indices_pair.first
                         << " times - indices size: "
                         << moving_atom_indices_pair.second.size() << std::endl;

            // Maybe a synthetic forward atom was made, and later on
            // in the dictionary, it was a real (normal) forward atom
            // of a different atom.  In that case there can be 2
            // copies of the atom in moving_atom_indices.  So let's
            // filter them one.  Only the first copy.
            // 
            std::vector<coot::map_index_t> unique_moving_atom_indices =
               uniquify_atom_indices(moving_atom_indices);

            if (debug) {
               std::cout << "  number of moving atoms based on atom index " << index3.index()
                         << " is " << moving_atom_indices.size() << std::endl;
               for (unsigned int imov=0; imov<unique_moving_atom_indices.size(); imov++) {
                  std::cout << "now pre-reverse  moving atom[" << imov << "] is "
                            << unique_moving_atom_indices[imov].index() << std::endl;
               }
            }

            // internal_reversed reversed_flag   action
            //      0                0             -
            //      0                1          complementary_indices
            //      1                0          complementary_indices
            //      1                1             -

            bool xor_reverse = reversed_flag ^ internal_reversed;
            if (xor_reverse) {
               unique_moving_atom_indices = complementary_indices(unique_moving_atom_indices,
                                                                  index2, index3);
               if (debug) {
                  for (unsigned int imov=0; imov<unique_moving_atom_indices.size(); imov++) {
                     std::cout << "now post-reverse  moving atom[" << imov << "] is "
                               << unique_moving_atom_indices[imov].index() << std::endl;
                  }
               }
            }

            unique_moving_atoms = unique_moving_atom_indices;
         }
      }
   }
   
   return unique_moving_atoms;

}

// give the user access, so that they know which atoms are moving.
std::vector<int>
coot::atom_tree_t::get_moving_atom_indices(const std::string &atom1,
                                           const std::string &atom2,
                                           bool reversed_flag) {
   std::vector<int> r;
   std::vector<coot::map_index_t> m = get_unique_moving_atom_indices(atom1, atom2, reversed_flag);
   for (unsigned int i=0; i<m.size(); i++) { 
      if (m[i].is_assigned())
         r.push_back(m[i].index());
   }
   return r;
}


// Throw exception on unable to rotate atoms.
//
// Return the new torsion angle (use the embedded torsion on index2 if you can)
// Otherwise return -1.0;
// 
// The angle is in degrees.
// 
double
coot::atom_tree_t::rotate_about(const std::string &atom1, const std::string &atom2,
                                double angle, bool reversed_flag) {

   bool debug = 0;
   double new_torsion = 0.0;
   // OK, so when the user clicks atom2 then atom1 (as the middle 2
   // atoms), then they implictly want the fragment revsersed,
   // relative to the atom order internally.  In that case, set
   // internal_reversed, and only if it is not set and the
   // reversed_flag flag *is* set, do the reversal (if both are set
   // they cancel each other out).  Now we do not pre-reverse the
   // indices in the calling function. The reverse is done here.
   bool internal_reversed = 0;

   if (debug)
      std::cout << "rotate_about() " << atom1 << " " << atom2 << std::endl;

   coot::map_index_t index2 = name_to_index[atom1];
   coot::map_index_t index3 = name_to_index[atom2];

   if ((atom_vertex_vec[index2.index()].forward.size() == 0) && 
       (atom_vertex_vec[index3.index()].forward.size() == 0)) {
      std::string s = "Neither index2 ";
      s += coot::util::int_to_string(index2.index());
      s += " nor index3 ";
      s += coot::util::int_to_string(index3.index());
      s += " has forward atoms!";
      throw std::runtime_error(s);
   } 
       
   if (index2.is_assigned()) { 
      if (index3.is_assigned()) {

         // is index3 in the forward atoms of index2?
         // if not, then swap and test again.
         //
         bool index3_is_forward = 0;

         for (unsigned int ifo=0; ifo<atom_vertex_vec[index2.index()].forward.size(); ifo++) {
            if (atom_vertex_vec[index2.index()].forward[ifo] == index3.index()) {
               index3_is_forward = 1;
               break;
            }
         }

         if (! index3_is_forward) {
            // perhaps index2 is the forward atom of index3?
            bool index2_is_forward = 0;
            for (unsigned int ifo=0; ifo<atom_vertex_vec[index3.index()].forward.size(); ifo++) {
               if (atom_vertex_vec[index3.index()].forward[ifo] == index2.index()) {
                  index2_is_forward = 1;
                  break;
               }
            }

            // if index2 *was* the forward atom of index3, then swap
            // around index2 and 3 into "standard" order.
            // 
            if (index2_is_forward) {
               if (debug)
                  std::cout << "    index2 was the forward atom of index3 - swapping "
                            << "and setting internal_reversed" << std::endl;
               std::swap(index2, index3);
               index3_is_forward = 1;
               internal_reversed = 1;
            }
         }

         // OK, try again
         if (index3_is_forward) {

            std::pair<int, std::vector<coot::map_index_t> > moving_atom_indices_pair =
               get_forward_atoms(index3, index3);
            std::vector<coot::map_index_t> moving_atom_indices =
               moving_atom_indices_pair.second;
            if (debug) 
               std::cout << " in rotate_about() get_forward_atoms() called "
                         << moving_atom_indices_pair.first
                         << " times - indices size: "
                         << moving_atom_indices_pair.second.size() << std::endl;

            // Maybe a synthetic forward atom was made, and later on
            // in the dictionary, it was a real (normal) forward atom
            // of a different atom.  In that case there can be 2
            // copies of the atom in moving_atom_indices.  So let's
            // filter them one.  Only the first copy.
            // 
            std::vector<coot::map_index_t> unique_moving_atom_indices =
               uniquify_atom_indices(moving_atom_indices);

            if (debug) {
               std::cout << "  number of moving atoms based on atom index " << index3.index()
                         << " is " << moving_atom_indices.size() << std::endl;
               for (unsigned int imov=0; imov<unique_moving_atom_indices.size(); imov++) {
                  std::cout << "now pre-reverse  moving atom[" << imov << "] is "
                            << unique_moving_atom_indices[imov].index() << std::endl;
               }
            }

            // internal_reversed reversed_flag   action
            //      0                0             -
            //      0                1          complementary_indices
            //      1                0          complementary_indices
            //      1                1             -

            bool xor_reverse = reversed_flag ^ internal_reversed;
            if (xor_reverse) {
               unique_moving_atom_indices = complementary_indices(unique_moving_atom_indices,
                                                                  index2, index3);
               if (debug)
                  for (unsigned int imov=0; imov<unique_moving_atom_indices.size(); imov++) {
                     std::cout << "now post-reverse  moving atom[" << imov << "] is "
                               << unique_moving_atom_indices[imov].index() << std::endl;
                  }
            }

            // so now we have a set of moving atoms:
            // set up the coordinates for the rotation, and rotate
            mmdb::PPAtom residue_atoms = 0;
            int n_residue_atoms;
            residue->GetAtomTable(residue_atoms, n_residue_atoms);
            mmdb::Atom *at2 = residue_atoms[index2.index()];
            mmdb::Atom *at3 = residue_atoms[index3.index()];
            clipper::Coord_orth base_atom_pos(at2->x, at2->y, at2->z);
            clipper::Coord_orth    third_atom(at3->x, at3->y, at3->z);;
            clipper::Coord_orth direction = third_atom - base_atom_pos;
            if (xor_reverse) {
               direction = base_atom_pos - third_atom;
               base_atom_pos = third_atom;
            }

            if (debug)
               std::cout << " calling rotate_internal() vector "
                         <<  direction.format() << " with base position "
                         << base_atom_pos.format() << " by " << angle
                         << std::endl;
            rotate_internal(unique_moving_atom_indices, direction, base_atom_pos,
                            clipper::Util::d2rad(angle));

            // set the new_torsion (return value) if possible.
            // 
            if (atom_vertex_vec[index2.index()].torsion_quad.first) {
               new_torsion = quad_to_torsion(index2);
            } 
         }
      } else {
         throw std::runtime_error("ERROR:: rotate_about(): index3 not assigned");
      } 
   } else {
      throw std::runtime_error("ERROR:: rotate_about(): index2 not assigned");
   } 
   return new_torsion;
}

// Throw exception on unable to rotate atoms.
//
// Return the new torsion angle (use the embedded torsion on index2 if you can)
// Otherwise return -1.0;
//
// We pass atoms 2 and 3 of the torsion
// 
double
coot::atom_tree_t::rotate_about(int index2, int index3, double angle, bool reversed_flag) {

   bool debug = false;
   double new_torsion = 0.0;

   // throw exception if not sane passed indices
   // 
   if (index2 == -1) {
      std::string s = "Ooops! rotate_about() Bad atom index: index-2";
      throw std::runtime_error(s);
   } 
   if (index3 == -1) {
      std::string s = "Ooops! rotate_about() Bad atom index: index-3";
      throw std::runtime_error(s);
   } 

   if ((atom_vertex_vec[index2].forward.size() == 0) && 
       (atom_vertex_vec[index3].forward.size() == 0)) {
      std::string s = "Neither index2 ";
      s += coot::util::int_to_string(index2);
      s += " nor index3 ";
      s += coot::util::int_to_string(index3);
      s += " has forward atoms!";
      throw std::runtime_error(s);
   }
   // OK, so when the user clicks atom2 then atom1 (as the middle 2
   // atoms), then they implictly want the fragment revsersed,
   // relative to the atom order internally.  In that case, set
   // internal_reversed, and only if it is not set and the
   // reversed_flag flag *is* set, do the reversal (if both are set
   // they cancel each other out).  Now we do not pre-reverse the
   // indices in the calling function. The reverse is done here.
   bool internal_reversed = false;

   // is index3 in the forward atoms of index2?
   // if not, then swap and test again.
   //
   bool index3_is_forward = false;

   for (unsigned int ifo=0; ifo<atom_vertex_vec[index2].forward.size(); ifo++) {
      if (atom_vertex_vec[index2].forward[ifo] == index3) {
         index3_is_forward = true;
         break;
      }
   }

   if (! index3_is_forward) {
      // perhaps index2 is the forward atom of index3?
      bool index2_is_forward = 0;
      for (unsigned int ifo=0; ifo<atom_vertex_vec[index3].forward.size(); ifo++) {
         if (atom_vertex_vec[index3].forward[ifo] == index2) {
            index2_is_forward = true;
            break;
         }
      }

      // if index2 *was* the forward atom of index3, then swap
      // around index2 and 3 into "standard" order.
      // 
      if (index2_is_forward) {
         if (debug)
            std::cout << "    index2 was the forward atom of index3 - swapping "
                      << "and setting internal_reversed" << std::endl;
         std::swap(index2, index3);
         index3_is_forward = true;
         internal_reversed = true;
      }
   }

   if (index3_is_forward) {
      std::pair<int, std::vector<coot::map_index_t> > moving_atom_indices_pair =
         get_forward_atoms(map_index_t(index3), map_index_t(index3));
      std::vector<coot::map_index_t> moving_atom_indices = moving_atom_indices_pair.second;

      if (debug) 
         std::cout << " in rotate_about(int, int) get_forward_atoms() called "
                   << moving_atom_indices_pair.first
                   << " times - indices size: "
                   << moving_atom_indices_pair.second.size() << std::endl;

      // Maybe a synthetic forward atom was made, and later on
      // in the dictionary, it was a real (normal) forward atom
      // of a different atom.  In that case there can be 2
      // copies of the atom in moving_atom_indices.  So let's
      // filter them one.  Only the first copy.
      // 
      std::vector<coot::map_index_t> unique_moving_atom_indices =
         uniquify_atom_indices(moving_atom_indices);

      if (debug) {
         std::cout << "  number of moving atoms based on atom index " << index3
                   << " is " << moving_atom_indices.size() << std::endl;
         for (unsigned int imov=0; imov<unique_moving_atom_indices.size(); imov++) {
            std::cout << "now pre-reverse  moving atom[" << imov << "] is "
                      << unique_moving_atom_indices[imov].index() << std::endl;
         }
      }

      // internal_reversed reversed_flag   action
      //      0                0             -
      //      0                1          complementary_indices
      //      1                0          complementary_indices
      //      1                1             -

      bool xor_reverse = reversed_flag ^ internal_reversed;
      if (xor_reverse) {
         unique_moving_atom_indices = complementary_indices(unique_moving_atom_indices,
                                                            map_index_t(index2), map_index_t(index3));
         if (debug)
            for (unsigned int imov=0; imov<unique_moving_atom_indices.size(); imov++) {
               std::cout << "now post-reverse  moving atom[" << imov << "] is "
                         << unique_moving_atom_indices[imov].index() << std::endl;
            }
      }

      mmdb::Atom *at2 = NULL; 
      mmdb::Atom *at3 = NULL;

      if (atom_selection) { 
         at2 = atom_selection[index2];
         at3 = atom_selection[index3];
      } 
      if (residue) {
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms;
         residue->GetAtomTable(residue_atoms, n_residue_atoms);
         at2 = residue_atoms[index2];
         at3 = residue_atoms[index3];
      } 

      if (at2 && at3) { 
         clipper::Coord_orth base_atom_pos(at2->x, at2->y, at2->z);
         clipper::Coord_orth    third_atom(at3->x, at3->y, at3->z);;
         clipper::Coord_orth direction = third_atom - base_atom_pos;
         if (xor_reverse) {
            direction = base_atom_pos - third_atom;
            base_atom_pos = third_atom;
         }

         if (debug)
            std::cout << " calling rotate_internal() vector "
                      <<  direction.format() << " with base position "
                      << base_atom_pos.format() << " by " << angle
                      << std::endl;
         rotate_internal(unique_moving_atom_indices, direction, base_atom_pos,
                         clipper::Util::d2rad(angle));

         // set the new_torsion (return value) if possible.
         // 
         if (atom_vertex_vec[index2].torsion_quad.first) {
            new_torsion = quad_to_torsion(map_index_t(index2));
         }
      } else {
         std::cout << "ERROR:: null atom rotate_about() - this should not happen"
                   << std::endl;
      } 
   }
   return new_torsion;
}



// return the torsion angle in degrees.
double
coot::atom_tree_t::quad_to_torsion(const coot::map_index_t &index2) const {

   // std::cout << " this bond has a chi angle assigned" << std::endl;
   coot::atom_index_quad quad = atom_vertex_vec[index2.index()].torsion_quad.second;
   clipper::Coord_orth co[4];
   mmdb::PPAtom residue_atoms;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   co[0] = clipper::Coord_orth(residue_atoms[quad.index1]->x,
                               residue_atoms[quad.index1]->y,
                               residue_atoms[quad.index1]->z);
   co[1] = clipper::Coord_orth(residue_atoms[quad.index2]->x,
                               residue_atoms[quad.index2]->y,
                               residue_atoms[quad.index2]->z);
   co[2] = clipper::Coord_orth(residue_atoms[quad.index3]->x,
                               residue_atoms[quad.index3]->y,
                               residue_atoms[quad.index3]->z);
   co[3] = clipper::Coord_orth(residue_atoms[quad.index4]->x,
                               residue_atoms[quad.index4]->y,
                               residue_atoms[quad.index4]->z);
   double ar = clipper::Coord_orth::torsion(co[0], co[1], co[2], co[3]);
   double new_torsion = clipper::Util::rad2d(ar);
   return new_torsion;
} 

// Back atoms, not very useful - it is only a single path.
// 
// some sort of recursion here
std::vector<coot::map_index_t>
coot::atom_tree_t::get_back_atoms(const coot::map_index_t &index) const {

   std::vector<coot::map_index_t> v;

   if (atom_vertex_vec[index.index()].connection_type == coot::atom_vertex::END)
      return v;
   
   // all of these
   for (unsigned int ib=0; ib<atom_vertex_vec[index.index()].backward.size(); ib++) {
      v.push_back(map_index_t(atom_vertex_vec[index.index()].backward[ib]));
   }

   // and the back atoms of back atoms.
   for (unsigned int ib=0; ib<atom_vertex_vec[index.index()].backward.size(); ib++) {
      int index_back = atom_vertex_vec[index.index()].backward[ib];
      coot::map_index_t back_index(index_back); // ho ho
      std::vector<coot::map_index_t> nv = 
         coot::atom_tree_t::get_back_atoms(map_index_t(index_back));
      for (unsigned int inv=0; inv<nv.size(); inv++)
         v.push_back(nv[inv]);
   }
   return v;
}


// with forward recursion
// 
std::pair<int, std::vector<coot::map_index_t> > 
coot::atom_tree_t::get_forward_atoms(const coot::map_index_t &base_index,
                                     const coot::map_index_t &index) const {

   bool debug = 0;
   
   std::vector<coot::map_index_t> v;
   int n_forward_count = 1; // this time at least

   // If atom_vertex_vec was not filled, then we should not index into
   // it with index.index():  Stops a crash, at least.
   //
   int n_atom_vertex = atom_vertex_vec.size(); // (fixed unsigned int vs int warning)
   if (index.index() >= n_atom_vertex)
      return std::pair<int, std::vector<coot::map_index_t> > (n_forward_count, v);


   // how can this happen?
   if (atom_vertex_vec[index.index()].connection_type == coot::atom_vertex::START) {
      // std::cout << " in get_forward_atoms index " <<  index.index() << " at start" << std::endl;
      return std::pair<int, std::vector<coot::map_index_t> > (n_forward_count, v);
   } 

   // Bizarre atom trees in the refmac dictionary - the END atom is
   // marks the C atom - not the O.  How can that be sane?
   // 
//    if (atom_vertex_vec[index.index()].connection_type == coot::atom_vertex::END) { 
//       // std::cout << " in get_forward_atoms index " <<  index.index() << " at end" << std::endl;
//       return v;
//    } 
   
   if (debug) {
      std::cout << "DEBUG::get_forward_atoms(): " << index.index() << " has "
                << atom_vertex_vec[index.index()].forward.size()
                << " forward atoms ";
      for (unsigned int ifo=0; ifo<atom_vertex_vec[index.index()].forward.size(); ifo++)
         std::cout << " " << atom_vertex_vec[index.index()].forward[ifo];
      std::cout << std::endl;
   }
   
   for (unsigned int ifo=0; ifo<atom_vertex_vec[index.index()].forward.size(); ifo++) { 
      v.push_back(map_index_t(atom_vertex_vec[index.index()].forward[ifo]));
      if (debug) 
         std::cout << " adding to forward atoms " << atom_vertex_vec[index.index()].forward[ifo]
                   << " which is atom_vertex_vec[" << index.index() << "].forward[" << ifo << "]"
                   << std::endl;
   }

   // and the forward atoms of the forward atoms
   for (unsigned int ifo=0; ifo<atom_vertex_vec[index.index()].forward.size(); ifo++) {
      int index_forward = atom_vertex_vec[index.index()].forward[ifo];
      coot::map_index_t forward_index(index_forward);
      std::pair<int, std::vector<coot::map_index_t> > nv_pair;
      if (base_index.index() != forward_index.index())
         nv_pair = coot::atom_tree_t::get_forward_atoms(base_index, forward_index);
      std::vector<coot::map_index_t> nv = nv_pair.second;
      n_forward_count += nv_pair.first;
      for (unsigned int inv=0; inv<nv.size(); inv++) {
         if (nv[inv].index() != base_index.index()) { 
            if (debug) 
               std::cout << " adding item  " << inv << " which has index " << nv[inv].index()
                         << " as a forward atom of " << index.index() << std::endl;
            v.push_back(nv[inv]);
         }
      }
   }

   return std::pair<int, std::vector<coot::map_index_t> > (n_forward_count, uniquify_atom_indices(v));
}



std::vector<coot::map_index_t>
coot::atom_tree_t::uniquify_atom_indices(const std::vector<coot::map_index_t> &vin) const {

   std::vector<coot::map_index_t> v;

   for (unsigned int iv=0; iv<vin.size(); iv++) {
      // vin[iv] is in v?
      bool found = 0;

      for (unsigned int ii=0; ii<v.size(); ii++) {
         if (vin[iv] == v[ii]) {
            found = 1;
            break;
         } 
      }

      if (found == 0)
         v.push_back(vin[iv]);
   } 
   return v; 
} 



// Return the complementary indices c.f. the moving atom set, but do
// not include index2 or index3 in the returned set (they do not move
// even with the reverse flag (of course)).
std::vector<coot::map_index_t> 
coot::atom_tree_t::complementary_indices(const std::vector<coot::map_index_t> &moving_atom_indices,
                                         const coot::map_index_t &index2,
                                         const coot::map_index_t &index3) const {

   std::vector<coot::map_index_t> v;
   for (int ivert=0; ivert<int(atom_vertex_vec.size()); ivert++) {
      bool ifound = 0;
      for (unsigned int im=0; im<moving_atom_indices.size(); im++) {
         if (moving_atom_indices[im].index() == ivert) {
            ifound = 1;
            break;
         }
      }
      if (ifound == 0)
         if (index2.index() != ivert)
            if (index3.index() != ivert)
               v.push_back(map_index_t(ivert));
   }
   
   return v;
} 




// so now we have a set of moving and non-moving atoms:
//
// the angle is in radians. 
void
coot::atom_tree_t::rotate_internal(std::vector<coot::map_index_t> moving_atom_indices,
                                   const clipper::Coord_orth &dir,
                                   const clipper::Coord_orth &base_atom_pos,
                                   double angle) {

   bool debug = 0;

   if (debug)
      std::cout << "in rotate_internal with " << moving_atom_indices.size()
                << " moving atoms " << std::endl;
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;

   if (residue) {
      residue->GetAtomTable(residue_atoms, n_residue_atoms);
   } else {
      // this atom_tree_t was generated from an (multi-residue) atom
      // selection then, just copy what was passed (well, generated
      // from the atom-selection-handle)
      residue_atoms = atom_selection;
      n_residue_atoms = n_selected_atoms;
   } 
   
   for (unsigned int im=0; im<moving_atom_indices.size(); im++) {
      int idx = moving_atom_indices[im].index();
      mmdb::Atom *at = residue_atoms[idx];
      clipper::Coord_orth po(at->x, at->y, at->z);
      clipper::Coord_orth pt = coot::util::rotate_around_vector(dir, po, base_atom_pos, angle);
      if (debug)
         std::cout << "  rotate_internal() moving atom number " << im << " " << at->name
                   << " from\n    " << at->x << "," << at->y << "," << at->z << " to "
                   << pt.format() << std::endl;
      at->x = pt.x();
      at->y = pt.y();
      at->z = pt.z();
   }
}

// can throw an exception
//
// pass an atom name quad.
double
coot::atom_tree_t::set_dihedral(const std::string &atom1, const std::string &atom2,
                                const std::string &atom3, const std::string &atom4,
                                double angle) {

    // debugging
    if (false) {
       std::map<std::string, map_index_t, std::less<std::string> >::const_iterator it;
       for (it=name_to_index.begin(); it!= name_to_index.end(); ++it) {
          std::cout << "set_dihedral() :" << it->first << ": -> " <<  it->second.index() << std::endl;
       }
   }

   coot::map_index_t i1 = name_to_index[atom1];
   coot::map_index_t i2 = name_to_index[atom2];
   coot::map_index_t i3 = name_to_index[atom3];
   coot::map_index_t i4 = name_to_index[atom4];

   if (i1.is_assigned() && i2.is_assigned() && i3.is_assigned() && i4.is_assigned()) {
      if (false)
         std::cout << "in atom_tree_t::set_dihedral() calling set_dihedral with indices  "
                   << i1.index() << " " << i2.index() << " " << i3.index() << " " << i4.index()
                   << " to angle " << angle << std::endl;
      return set_dihedral(i1, i2, i3, i4, angle); // can throw an exception
   } else {
      std::string mess = "Atom name(s) not found in residue. ";

      std::vector<std::string> unassigned;
      if (! i1.is_assigned()) unassigned.push_back(atom1);
      if (! i2.is_assigned()) unassigned.push_back(atom2);
      if (! i3.is_assigned()) unassigned.push_back(atom3);
      if (! i4.is_assigned()) unassigned.push_back(atom4);
      if (unassigned.size() > 0) {
         mess += "Unassigned atoms: ";
         for (unsigned int i=0; i<unassigned.size(); i++) {
            mess += "\"";
            mess += unassigned[i];
            mess += "\"  ";
         }
      }
      throw std::runtime_error(mess);
      // no return from this path.
   }
}


// can throw an exception.
//
// i1, i2, i3, i4 are not checked for having been assigned.  We
// presume that that has been done before we get here.
//
// pass an atom index quad.
//
double
coot::atom_tree_t::set_dihedral(const coot::map_index_t &i1,
                                const coot::map_index_t &i2,
                                const coot::map_index_t &i3,
                                const coot::map_index_t &i4,
                                double angle) {

   double dihedral_angle = 0.0;
   try {
      coot::atom_index_quad iq(i1.index(), i2.index(), i3.index(), i4.index());

      double current_dihedral_angle = -1000; // unset
      if (residue)
         current_dihedral_angle = iq.torsion(residue);
      if (atom_selection)
         current_dihedral_angle = iq.torsion(atom_selection, n_selected_atoms);
      // this should not happen
      if (current_dihedral_angle == -1000)
         throw std::runtime_error("bad current_dihedral_angle, no resiude or selection?");

      double diff = angle - current_dihedral_angle;
      if (diff > 360.0)
         diff -= 360.0;
      if (diff < -360.0)
         diff += 360.0;
      // take out this try/catch when done.
      rotate_about(i2.index(), i3.index(), diff, 0);
      dihedral_angle = iq.torsion(residue);
      if (0)
         std::cout << "   current, target, diff new "
                   << current_dihedral_angle << "  " << angle << "  " << diff << "  "
                   << dihedral_angle << std::endl;
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
      std::string mess = "Torsion failure for index ";
      mess += util::int_to_string(i1.index());
      mess += " to ";
      mess += util::int_to_string(i2.index());
      mess += " to ";
      mess += util::int_to_string(i3.index());
      mess += " to ";
      mess += util::int_to_string(i4.index());
      mess += ": rotate_about() fails.";
      throw std::runtime_error(mess);
   }
   return dihedral_angle;
}

// this can throw an exception
//
// return the angle of the torsion (should be the same as was set)
//
double
coot::atom_tree_t::set_dihedral(const coot::atom_quad &quad, double angle,
                                bool reverse_flag) {

   if (false) {
      std::cout << "debug:: in set_dihedral() A "
                << quad.atom_1 << " " << quad.atom_2 << " "
                << quad.atom_3 << " " << quad.atom_4 << std::endl;
      int idx1 = get_index(quad.atom_1).index();
      int idx2 = get_index(quad.atom_2).index();
      int idx3 = get_index(quad.atom_3).index();
      int idx4 = get_index(quad.atom_4).index();
      std::cout << "debug:: in set_dihedral() B " << idx1 << " " << idx2 << " " << idx3 << "  " << idx4 << std::endl;
   }

   double dihedral_angle = 0.0;
   double current_dihedral_angle = quad.torsion();
   double diff = angle - current_dihedral_angle;
   if (diff >  360.0) diff -= 360.0;
   if (diff < -360.0) diff += 360.0;
   int ind_2 = get_index(quad.atom_2).index();
   int ind_3 = get_index(quad.atom_3).index();
   if (ind_2 == -1) throw std::runtime_error("set_dihedral(quad) missing atom 2");
   if (ind_3 == -1) throw std::runtime_error("set_dihedral(quad) missing atom 3");
   // std::cout << "rotate_about " << ind_2 << " " << ind_3 << " " << diff << std::endl;
   rotate_about(ind_2, ind_3, diff, reverse_flag);

   dihedral_angle = quad.torsion();
   return dihedral_angle;
}


// this can throw an exception
// 
// return the set of angles - should be the same that they were
// set to (for validation).
std::vector<double>
coot::atom_tree_t::set_dihedral_multi(const std::vector<tree_dihedral_info_t> &di) {

   // tree_dihedral_info_t is an atom_name_quad and a angle in degrees.
   //
   std::vector<double> v(di.size());
   for (unsigned int id=0; id<di.size(); id++) {
      try {
         v[id] = set_dihedral(di[id].quad.atom_name(0), di[id].quad.atom_name(1),
                              di[id].quad.atom_name(2), di[id].quad.atom_name(3),
                              di[id].dihedral_angle);
      }
      catch (const std::exception &e) {
         std::cout << "WARNING:: " << e.what() << std::endl;
      }
   }
   return v;
}

std::vector<double>
coot::atom_tree_t::set_dihedral_multi(const std::vector<tree_dihedral_quad_info_t> &quads) {

   std::vector<double> v(quads.size());
   for (unsigned int iquad=0; iquad<quads.size(); iquad++) {
      coot::map_index_t index2 = get_index(quads[iquad].quad.atom_2);
      coot::map_index_t index3 = get_index(quads[iquad].quad.atom_3);
      if (false)
         std::cout << "testing if  " << quads[iquad].fixed_atom_index.index() << " "
                   << coot::atom_spec_t(atom_selection[quads[iquad].fixed_atom_index.index()]) << " "
                   << "is a forward atom of "
                   << index2.index() << " "
                   << coot::atom_spec_t(atom_selection[index2.index()]) << " "
                   << std::endl;
      bool iff = in_forward_atoms(index2, quads[iquad].fixed_atom_index);

      bool index3_is_forward_of_index2 = 0;
      for (unsigned int ifo=0; ifo<atom_vertex_vec[index2.index()].forward.size(); ifo++) {
         if (atom_vertex_vec[index2.index()].forward[ifo] == index3.index()) {
            index3_is_forward_of_index2 = true;
            break;
         }
      }
      if (0)
         std::cout << "   result: " << iff << "  index3_is_forward_of_index2: "
                   << index3_is_forward_of_index2 << std::endl;

      bool xor_reverse = ! (iff ^ index3_is_forward_of_index2);
      v[iquad] = set_dihedral(quads[iquad].quad, quads[iquad].dihedral_angle, xor_reverse);
   }
   return v;
}

bool
coot::atom_tree_t::in_forward_atoms(const coot::map_index_t &bond_atom_index,
                                    const coot::map_index_t &fixed) const {

   bool status = false;

   if (fixed.is_assigned()) {
      std::pair<int, std::vector<coot::map_index_t> >  p =
         get_forward_atoms(bond_atom_index, bond_atom_index);
      if (std::find(p.second.begin(), p.second.end(), fixed) != p.second.end())
         status = true;
   }
   return status;
}


coot::map_index_t
coot::atom_tree_t::get_index(mmdb::Atom *atom) const {

   coot::map_index_t idx;
   // std::cout << "debug:: in get_index() residue is " << residue << " and atom_selection is " << atom_selection << std::endl;
   if (residue) {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         if (residue_atoms[iat] == atom) {
            idx = map_index_t(iat);
            break;
         }
      }
   }
   if (atom_selection) {
      for (int iat=0; iat<n_selected_atoms; iat++) {
         if (false)
            std::cout << "::: in get_index() comparing atom " << atom << " " << coot::atom_spec_t(atom)
                      << " with atom-seletion-atom: " <<  iat << " " << atom_selection[iat] << " "
                      << coot::atom_spec_t(atom_selection[iat]) << std::endl;
         if (atom_selection[iat] == atom) {
            idx = map_index_t(iat);
            break;
         }
      }
   }

   return idx;
}



// For use with above function, where the class constructor is made
// with a minimol::residue.
// 
coot::minimol::residue
coot::atom_tree_t::GetResidue() const {

   return coot::minimol::residue(residue);
}

std::ostream&
coot::operator<<(std::ostream &o, coot::atom_tree_t::tree_dihedral_info_t t) { 

   o << "[dihedral-info: " << t.quad << " " << t.dihedral_angle << "]";
   return o;
}


