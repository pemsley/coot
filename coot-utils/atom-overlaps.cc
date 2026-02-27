
/* coot-utils/atom-overlaps.cc
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

#include <string.h>

#include "mmdb2/mmdb_uddata.h"
#include <stdlib.h> // needed from abs()

#include <iomanip>


#ifdef HAVE_CXX_THREAD
#include <thread>
#include <future>
#endif

#include "utils/coot-utils.hh"
#include "atom-overlaps.hh"
#include "coot-coord-utils.hh"
#include "geometry/main-chain.hh"

#include "utils/logging.hh"
extern logging logger;

coot::atom_overlaps_container_t::atom_overlaps_container_t(mmdb::Residue *res_central_in,
                                                           const std::vector<mmdb::Residue *> &neighbours_in,
                                                           mmdb::Manager *mol_in,
                                                           const protein_geometry *geom_p_in) {
   geom_p = geom_p_in;
   res_central = res_central_in;
   neighbours = neighbours_in;
   mol = mol_in;
   imol_enc = protein_geometry::IMOL_ENC_ANY;
   clash_spike_length = 0.5;
   init();

}

// is this used?
coot::atom_overlaps_container_t::atom_overlaps_container_t(mmdb::Residue *res_central_in,
                                                           mmdb::Residue *neighbour,
                                                           mmdb::Manager *mol_in,
                                                           const protein_geometry *geom_p_in) {
   geom_p = geom_p_in;
   res_central = res_central_in;
   neighbours.push_back(neighbour);
   mol = mol_in;
   imol_enc = protein_geometry::IMOL_ENC_ANY; // hack
   clash_spike_length = 0.5;
   init();

}

// this one for contact dots (around central ligand)
//
// this can throw a std::out_of_range (missing residue from dictionary)
//
// clash spike_length should be 0.5;
coot::atom_overlaps_container_t::atom_overlaps_container_t(mmdb::Residue *res_central_in,
                                                           const std::vector<mmdb::Residue *> &neighbours_in,
                                                           mmdb::Manager *mol_in,
                                                           int imol_enc_in,
                                                           const protein_geometry *geom_p_in,
                                                           double clash_spike_length_in,
                                                           double probe_radius_in) {
   probe_radius = probe_radius_in;
   geom_p = geom_p_in;
   res_central = res_central_in;
   neighbours = neighbours_in;
   mol = mol_in;
   imol_enc = imol_enc_in,
   clash_spike_length = clash_spike_length_in;
   init();

}

// this can throw a std::out_of_range (missing residue from dictionary)
// /default args for clash_spike_length_in and probe_radius_in.
coot::atom_overlaps_container_t::atom_overlaps_container_t(mmdb::Manager *mol_in,
                                                           const protein_geometry *geom_p_in,
                                                           bool ignore_water_contacts_flag_in,
                                                           double clash_spike_length_in,
                                                           double probe_radius_in) {

   geom_p = geom_p_in;
   res_central = 0;
   mol = mol_in;
   imol_enc = protein_geometry::IMOL_ENC_ANY; // hack
   clash_spike_length = 0.5;
   probe_radius = probe_radius_in;
   ignore_water_contacts_flag = ignore_water_contacts_flag_in;
   init_for_all_atom();
}



// this can throw a std::out_of_range (missing residue from dictionary)
//
void
coot::atom_overlaps_container_t::init() {

   overlap_mode = CENTRAL_RESIDUE;

   udd_residue_index_handle = -1; // unset

   have_dictionary = false; // initially.
   molecule_has_hydrogens = false; // initially

   if (res_central) {

      std::string cres_name = res_central->GetResName();
      std::pair<bool, dictionary_residue_restraints_t> d = geom_p->get_monomer_restraints(cres_name, imol_enc);
      if (! d.first) {
         std::cout << "WARNING:: (or ERROR::) in atom_overlaps_container_t::init() Failed to get dictionary for "
                   << cres_name << std::endl;
      } else {
         // Happy path
         central_residue_dictionary = d.second;

         if (false)
            std::cout << "central_residue_dictionary has " << central_residue_dictionary.atom_info.size()
                      << " atoms and " << central_residue_dictionary.bond_restraint.size()
                      << " bond restraints " << std::endl;

         neighb_dictionaries.resize(neighbours.size());
         have_dictionary = true;
         for (unsigned int i=0; i<neighbours.size(); i++) {
            std::string residue_name = neighbours[i]->GetResName();
            // std::cout << "testing for " << residue_name << std::endl;
            d = geom_p->get_monomer_restraints(residue_name, imol_enc);
            if (! d.first) {
               std::cout << "WARNING:: Overlap fail. Failed to get dictionary for name "
                         << residue_name << std::endl;
               have_dictionary = false;
               break;
            } else {
               // happy path
               neighb_dictionaries[i] = d.second;
            }
         }
      }

      if (have_dictionary) {
         fill_ligand_atom_neighbour_map(); // and add radius

         mark_donors_and_acceptors();
      }
   }

}

// this can throw a std::out_of_range (missing residue from dictionary)
//
void
coot::atom_overlaps_container_t::init_for_all_atom() {

   // 20250225-PE the molecule could be re-used in a new atom_overlaps_container_t
   // so let's check if we have a udd for hb_type and if so, we can skip this
   // init() - because the molecule has already been setup.
   // int udd_test_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "hb_typeadsfadf");
   // std::cout << "::::::::::::::::::::::::: udd_test_handle " << udd_test_handle << std::endl;

   overlap_mode = ALL_ATOM;

   int udd_handle_t2 = mol->GetUDDHandle(mmdb::UDR_HIERARCHY, "setup-hydrogen-types");
   // std::cout << "::::: udd_handle_t2 " << udd_handle_t2 << std::endl;
   if (udd_handle_t2 == 0) {
      udd_handle_t2 = mol->RegisterUDInteger(mmdb::UDR_HIERARCHY, "setup-hydrogen-types");
      mol->PutUDData(udd_handle_t2, 11);
   } else {
      // std::cout << ":::: udd_handle_t2 found - returning early" << std::endl;
      // return;
   }

   // only do this once
   udd_residue_index_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "neighb-residue-index");

   // neighbours is a misnomer in this case - it is merely a list of all residue

   have_dictionary = true; // initially.
   molecule_has_hydrogens = false; // initially

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int neigh_idx = neighbours.size();
                  neighbours.push_back(residue_p);

                  // for lookups in setup_env_residue_atoms_radii()
                  mmdb::Atom **residue_atoms = 0;
                  int n_residue_atoms;
                  residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                  for (int iat=0; iat<n_residue_atoms; iat++) {
                     mmdb::Atom *at = residue_atoms[iat];
                     at->PutUDData(udd_residue_index_handle, neigh_idx);
                  }

                  std::string residue_name(residue_p->GetResName());
                  std::map<std::string, dictionary_residue_restraints_t>::const_iterator it;
                  it = dictionary_map.find(residue_name);
                  if (it == dictionary_map.end()) {
                     std::pair<bool, dictionary_residue_restraints_t> d =
                        geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);
                     if (! d.first) {
                        std::cout << "WARNING::Failed to get dictionary for " << residue_name << std::endl;
                        std::cout << "WARNING:: turning off have_dictionary" << std::endl;
                        have_dictionary = false;
                     } else {
                        dictionary_map[residue_name] = d.second;
                     }
                  }
               }
            }
         }
      }
   }

   // now mark donors and acceptors.
   //
   udd_h_bond_type_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "hb_type");
   // std::cout << "::::: registered ud integer handle: " << udd_h_bond_type_handle
   //           << " for mol " << mol << std::endl;

   mark_donors_and_acceptors_for_neighbours(udd_h_bond_type_handle);

}

void
coot::atom_overlaps_container_t::mark_donors_and_acceptors() {

   // now mark donors and acceptors.
   //
   udd_h_bond_type_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "hb_type");

   mark_donors_and_acceptors_central_residue(udd_h_bond_type_handle);

   mark_donors_and_acceptors_for_neighbours(udd_h_bond_type_handle);
}

void
coot::atom_overlaps_container_t::mark_donors_and_acceptors_central_residue(int udd_h_bond_type_handle) {

   if (res_central) {
      mmdb::PAtom *central_residue_atoms = 0;
      int n_central_residue_atoms;
      res_central->GetAtomTable(central_residue_atoms, n_central_residue_atoms);
      for (int iat=0; iat<n_central_residue_atoms; iat++) {
         mmdb::Atom *at = central_residue_atoms[iat];
         std::string atom_name(at->name);
         std::string ele = at->element;
         if (ele == " H") {
            molecule_has_hydrogens = true;
            // Hydrogens have energy type "H" from Refmac and acedrg, that doesn't
            // tell us if this atom is a donor hydrogen.
            // So, find the atom to which the H is attached and if that is a donor then this
            // is a hydrogen bond hydrogen.
            std::string heavy_neighb_of_H_atom =
               central_residue_dictionary.get_bonded_atom(atom_name);
            if (! heavy_neighb_of_H_atom.empty()) {
               std::string neigh_energy_type = central_residue_dictionary.type_energy(heavy_neighb_of_H_atom);
               energy_lib_atom neighb_ela = geom_p->get_energy_lib_atom(neigh_energy_type);
               hb_t neighb_hb_type = neighb_ela.hb_type;
               if (neighb_hb_type == coot::HB_DONOR) {
                  // std::cout << "----- adding ligand HB_HYDROGEN udd " << atom_spec_t(at) << std::endl;
                  at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN); // hb_t -> int
               }
               if (neighb_hb_type == coot::HB_BOTH) {
                  // std::cout << "----- adding ligand HB_HYDROGEN udd " << atom_spec_t(at) << std::endl;
                  at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN); // hb_t -> int
               }
               // additionally NR5 can be a donor (because the real type should be NR15)
               if (neigh_energy_type == "NR5") {
                  at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN);
               }
            }
         } else {
            std::string energy_type = central_residue_dictionary.type_energy(atom_name);
            energy_lib_atom ela = geom_p->get_energy_lib_atom(energy_type);
            hb_t hb_type = ela.hb_type;
            at->PutUDData(udd_h_bond_type_handle, hb_type); // hb_t -> int
            if (false)
               std::cout << "marked ligand atom " << coot::atom_spec_t(at) << " as hb_type " << hb_type
                         << " note HB_DONOR " << HB_DONOR << " HB_ACCEPTOR " << HB_ACCEPTOR << std::endl;
            if (energy_type == "NR5") {
              // test here that there is a hydrogen atom bonded to this one
              if (true) {
                at->PutUDData(udd_h_bond_type_handle, coot::HB_DONOR); // hack, for above reasons
              }
            }
         }
      }
   }
}

// this can throw a std::out_of_range
//
const coot::dictionary_residue_restraints_t &
coot::atom_overlaps_container_t::get_dictionary(mmdb::Residue *r, unsigned int idx) const {

   if (false)
      std::cout << ":::::::::::::::::::::::: in get_dictionary() the dictionary_map is of size "
                << dictionary_map.size() << " and idx is " << idx
                << " overlap_mode " << overlap_mode  << " vs ALL_ATOM " << ALL_ATOM << std::endl;

   if (overlap_mode == ALL_ATOM) {

      std::string res_name = r->GetResName();
      std::map<std::string, dictionary_residue_restraints_t>::const_iterator it =
         dictionary_map.find(res_name);

      if (it == dictionary_map.end()) {
         std::cout << "========= hideous failure in get_dictionary for type " << res_name
                   << " using " << dictionary_map.size() << " dictionary entries" << std::endl;
         std::string mess = "dictionary index out of range for ";
         mess += res_name;
         std::out_of_range oor(mess);
         throw oor;
      }
      return it->second;
   } else {
      return neighb_dictionaries[idx];
   }
}

// this can throw a std::out_of_range
//
void
coot::atom_overlaps_container_t::mark_donors_and_acceptors_for_neighbours(int udd_h_bond_type_handle) {

   // std::cout << "mark_donors_and_acceptors_for_neighbours() " << udd_h_bond_type_handle << std::endl;

   if (false) {
      std::cout << "mark_donors_and_acceptors_for_neighbours() here are the neigbhors"
                << std::endl;
      for (unsigned int i=0; i<neighbours.size(); i++) {
         std::cout << "    " << coot::residue_spec_t(neighbours[i]) << std::endl;
      }
   }

   for (unsigned int i=0; i<neighbours.size(); i++) {
      // const dictionary_residue_restraints_t &dict = neighb_dictionaries[i];
      try {
         const dictionary_residue_restraints_t &dict = get_dictionary(neighbours[i], i);
         mmdb::PAtom *residue_atoms = 0;
         int n_residue_atoms;
         neighbours[i]->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *n_at = residue_atoms[iat];
            std::string atom_name(n_at->name);
            std::string ele = n_at->element;
            if (ele == " H") {
               molecule_has_hydrogens = true;
               // as above
               std::string heavy_neighb_of_H_atom = dict.get_bonded_atom(atom_name);
               if (! heavy_neighb_of_H_atom.empty()) {
                  std::string neigh_energy_type = dict.type_energy(heavy_neighb_of_H_atom);
                  energy_lib_atom neighb_ela = geom_p->get_energy_lib_atom(neigh_energy_type);
                  hb_t neighb_hb_type = neighb_ela.hb_type;
                  if (false)
                     std::cout << "222 mark_donors_and_acceptors_for_neighbours() "
                               << "heavy_neighb_of_H_atom " << heavy_neighb_of_H_atom << " neigh_hb_type for H-atom "
                               << atom_name << " is " << neighb_hb_type << " with neighb_ela type " << neighb_ela.type << std::endl;
                  if (neighb_hb_type == coot::HB_DONOR) {
                     // std::cout << "----- adding env HB_HYDROGEN udd " << atom_spec_t(n_at) << std::endl;
                     n_at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN); // hb_t -> int
                  }
                  if (neighb_hb_type == coot::HB_BOTH) {
                     // std::cout << "----- adding env HB_HYDROGEN udd " << atom_spec_t(n_at) << std::endl;
                     n_at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN); // hb_t -> int
                  }
                  // additionally NR5 can be a donor (because the real type should be NR15)
                  if (neigh_energy_type == "NR5") {
                    n_at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN);
                  }
               } else {
                  std::cout << "ERROR:: mark_donors_and_acceptors_for_neighbours() oops empty heavy for H-atom "
                            << atom_name << std::endl;
               }
            } else {
               std::string energy_type = dict.type_energy(atom_name);
               energy_lib_atom ela = geom_p->get_energy_lib_atom(energy_type);
               hb_t hb_type = ela.hb_type;
               if (false)
                  std::cout << "----- adding env hb_type " << hb_type << " udd " << atom_spec_t(n_at) << " "
                            << " using handle " << udd_h_bond_type_handle << std::endl;
               if (n_at->PutUDData(udd_h_bond_type_handle, hb_type) == mmdb::UDDATA_Ok) { // hb_t -> int
                  // std::cout << "putting uddata for hb_type " << hb_type << " OK" << std::endl;
               } else {
                  std::cout << "ERROR:: mark_donors_and_acceptors_for_neighbours() putting uddata for hb_type "
                            << hb_type << " fail" << std::endl;
               }
               if (energy_type == "NR5") {
                 // test here that there is a hydrogen atom bonded to this one
                 if (true) {
                   n_at->PutUDData(udd_h_bond_type_handle, coot::HB_DONOR); // hack, for above reasons
                 }
               }
               if (false)
                  std::cout << "........... type_energy for atom name " << atom_name
                            << " in dictionary for comp-id \"" << dict.residue_info.comp_id
                            << "\" put type " << hb_type << std::endl;
            }
         }
      }
      catch (const std::out_of_range &ex) {
         std::cout << "Opps " << ex.what() << std::endl;
      }
   }
}

// and radius
void
coot::atom_overlaps_container_t::fill_ligand_atom_neighbour_map() {

   mmdb::realtype max_dist = 2.3;
   if (mol) {
      mmdb::Contact *pscontact = NULL;
      int n_contacts;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
         for (int j=0; j<4; j++)
            my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      res_central->GetAtomTable(residue_atoms, n_residue_atoms);

      mol->SeekContacts(residue_atoms, n_residue_atoms,
                        residue_atoms, n_residue_atoms,
                        0, max_dist,
                        0, // in same residue
                        pscontact, n_contacts,
                        0, &my_matt, i_contact_group); // makes reverses also
      if (n_contacts > 0) {
         if (pscontact) {
            for (int i=0; i<n_contacts; i++) {
               mmdb::Atom *neighb_at = residue_atoms[pscontact[i].id2];
               double radius = get_vdw_radius_ligand_atom(neighb_at);
               std::pair<mmdb::Atom *, double> p(neighb_at, radius);
               ligand_atom_neighbour_map[pscontact[i].id1].push_back(p);
            }
         }
      }
   }
}




void
coot::atom_overlaps_container_t::make_overlaps() {

   if (false)
      std::cout << ":::::::::::::::::: make_overlaps()" << std::endl;

   if (! have_dictionary) {
      std::cout << "WARNING:: make_overlaps(): No dictionary!" << std::endl;
      return;
   }

   // This is completely non-clever about finding which atoms are close
   // to each other - it checks all neighbour atoms against all
   // central-residue atoms.

   double dist_crit = 2.3; // Everything that is bonded is less than this
   double dist_crit_sqrd = (2*dist_crit) * (2*dist_crit);

   // so that we can have a sorted list of baddies

   class baddie_attribs_t {
   public:
      mmdb::Atom *cr_at;
      mmdb::Atom *n_at;
      float r_1;
      float r_2;
      float d;
      float o;
      bool is_hydrogen_bond;
      bool hydrogen_atom_is_first_atom;
      baddie_attribs_t(mmdb::Atom *cr_at, mmdb::Atom *n_at, float r_1, float r_2, float d, float o,
                       bool is_hydrogen_bond, bool hydrogen_atom_is_first_atom) :
         cr_at(cr_at), n_at(n_at), r_1(r_1), r_2(r_2), d(d), o(o),
         is_hydrogen_bond(is_hydrogen_bond), hydrogen_atom_is_first_atom(hydrogen_atom_is_first_atom) {}
      static bool sorter(const baddie_attribs_t &b1, const baddie_attribs_t &b2) {
         return b2.o < b1.o;
      }
   };

   mmdb::PAtom *central_residue_atoms = 0;
   int n_central_residue_atoms;
   if (! res_central) return; // it isn't if we have all atom
   res_central->GetAtomTable(central_residue_atoms, n_central_residue_atoms);
   std::string res_central_name(res_central->GetResName());
   std::vector<std::pair<std::string, std::string> > bonds_for_cr_at =
      geom_p->get_bonded_and_1_3_angles(res_central_name, protein_geometry::IMOL_ENC_ANY);

   std::vector<baddie_attribs_t> baddies; // and not so baddies

   for (int j=0; j<n_central_residue_atoms; j++) {

      mmdb::Atom *cr_at = central_residue_atoms[j];
      // std::cout << "Surface points for ligand atom " << atom_spec_t(cr_at) << std::endl;
      clipper::Coord_orth co_cr_at = co(cr_at);
      double r_1 = get_vdw_radius_ligand_atom(cr_at);

      for (unsigned int i=0; i<neighbours.size(); i++) {

         std::string res_name_nat(neighbours[i]->GetResName());
         std::vector<std::pair<std::string, std::string> > bonds_for_n_at =
            geom_p->get_bonded_and_1_3_angles(res_name_nat, protein_geometry::IMOL_ENC_ANY);

         mmdb::PAtom *residue_atoms = 0;
         int n_residue_atoms = 0;
         neighbours[i]->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *n_at = residue_atoms[iat];
            clipper::Coord_orth co_n_at = co(n_at);

            double ds = (co_cr_at - co_n_at).lengthsq();
            if (ds < dist_crit_sqrd) {

               // duck out if this is linked

               if (is_linked(cr_at, n_at))
                  continue;
               if (is_angle_related_via_link(cr_at, n_at, bonds_for_cr_at, bonds_for_n_at))
                  continue;

               double r_2 = get_vdw_radius_neighb_atom(n_at, i);
               double d = sqrt(ds);

               // first is yes/no, second is if the H is on the ligand
               //
               // std::pair<bool, bool> might_be_h_bond_flag =
               // is_h_bond_H_and_acceptor(cr_at, n_at, udd_h_bond_type_handle);
               h_bond_info_t hbi(cr_at, n_at, udd_h_bond_type_handle);

               if (d < (r_1 + r_2 + probe_radius)) {
                  double o = get_overlap_volume(d, r_2, r_1);
                  bool h_bond_flag = false;
                  if (hbi.is_h_bond_H_and_acceptor) {
                     h_bond_flag = true;
                     if (hbi.H_is_first_atom_flag) {
                        baddies.push_back(baddie_attribs_t(cr_at, n_at, r_1, r_2, d, o, true, hbi.H_is_first_atom_flag));
                        // std::cout << "INFO:: " << atom_spec_t(cr_at) << "" << " and " << atom_spec_t(n_at)
                        // << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d
                        // << " overlap " << o << " IS H-Bond (ligand donor)" << std::endl;
                     } else {
                        baddies.push_back(baddie_attribs_t(cr_at, n_at, r_1, r_2, d, o, true, hbi.H_is_first_atom_flag));
                        // std::cout << "INFO:: " << atom_spec_t(cr_at) << "   " << " and " << atom_spec_t(n_at)
                        // << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d
                        // << " overlap " << o << " IS H-Bond (ligand acceptor)" << std::endl;
                     }
                  } else {
                     baddies.push_back(baddie_attribs_t(cr_at, n_at, r_1, r_2, d, o, false, false));
                     // std::cout << "INFO:: " << atom_spec_t(cr_at) << "" << " and " << atom_spec_t(n_at)
                     // << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d
                     // << " overlap " << o << std::endl;
                  }
                  atom_overlap_t ao(j, cr_at, n_at, r_1, r_2, o);
                  // atom_overlap_t ao2(-1, n_at, cr_at, r_2, r_1, o);
                  ao.is_h_bond = h_bond_flag;
                  overlaps.push_back(ao);
                  // overlaps.push_back(ao2);
               } else {
                  if (hbi.is_h_bond_H_and_acceptor) {
                     if (d < (dist_crit + 0.5)) {
                        // std::cout << "INFO:: " << atom_spec_t(cr_at) << "" << " and " << atom_spec_t(n_at)
                        //           << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d
                        //           << " but might be h-bond anyway (is this strange?)"
                        //           << std::endl;
                        logger.log(log_t::INFO, atom_spec_t(cr_at).format(), "and", atom_spec_t(n_at).format(),
                                   "r_1", r_1, "and r_2", r_2, "and d", d,
                                   "but might be h-bond anyway (is this strange?)");
                        // double o = 0;
                        // atom_overlap_t ao(cr_at, n_at, r_1, r_2, o);
                        // ao.is_h_bond = true;
                        // overlaps.push_back(ao);
                     }
                  }
               }
            }
         }
      }
   }

   std::sort(baddies.begin(), baddies.end(), baddie_attribs_t::sorter);
   for(auto const &b : baddies) {
      // std::cout << "INFO:: make_overlaps(): baddie: "
      //           << atom_spec_t(b.cr_at) << "" << " and " << atom_spec_t(b.n_at)
      //           << " r_1 " << std::left << std::setw(4) << b.r_1 << " and r_2 "
      //           << std::left << std::setw(4) << b.r_2
      //           <<  " and d " << std::setw(4) << b.d << " overlap " << b.o;
      // // << " IS H-Bond (ligand donor)" << std::endl;
      // if (b.is_hydrogen_bond) {
      //    if (b.hydrogen_atom_is_first_atom)
      //       std::cout << " is H-bond (ligand donor)";
      //    else
      //       std::cout << " is H-bond (ligand acceptor)";
      // }
      // std::cout << std::endl;
      std::string hb_str;
      if (b.is_hydrogen_bond) {
         if (b.hydrogen_atom_is_first_atom)
            hb_str = "is H-bond (ligand donor)";
         else
            hb_str = "is H-bond (ligand acceptor)";
      }
      logger.log(log_t::INFO, "make_overlaps(): baddie: " + atom_spec_t(b.cr_at).format() +
                 " and " + atom_spec_t(b.n_at).format() +
                 " r_1 " + std::to_string(b.r_1) + " and r_2 " + std::to_string(b.r_2) +
                 " and d " + std::to_string(b.d) + " overlap " + std::to_string(b.o) + " " + hb_str);
   }
}

bool
coot::atom_overlaps_container_t::kludge_filter(mmdb::Atom *at_1, mmdb::Atom *at_2) const {

   // C2 (NAG) and NE2 (ASN) are angle related - but we don't know that unless we
   // generate full restraints - which we don't do.

   bool reject = false;
   if (at_1->residue->chain == at_2->residue->chain) {
      std::string res_name_1(at_1->residue->GetResName());
      if (res_name_1 == "ASN") {
         std::string res_name_2(at_2->residue->GetResName());
         if (res_name_2 == "NAG") {
            std::string atom_name_1(at_1->name);
            if (atom_name_1 == " NE2") {
               std::string atom_name_2(at_2->name);
               if (atom_name_2 == " C2 ") {
                  reject = true;
               }
            }
         }
      }
      if (res_name_1 == "NAG") {
         std::string res_name_2(at_2->residue->GetResName());
         if (res_name_2 == "ASN") {
            std::string atom_name_1(at_1->name);
            if (atom_name_1 == " C2 ") {
               std::string atom_name_2(at_2->name);
               if (atom_name_2 == " NE2") {
                  reject = true;
               }
            }
         }
      }
   }

   return reject;

}


// modifies overlaps
void
coot::atom_overlaps_container_t::make_all_atom_overlaps() {

   if (! have_dictionary) {
      std::cout << "No dictionary" << std::endl;
      return;
   }

   if (mol) {
      mmdb::realtype max_dist = 1.75 + 1.75 + 2.0 * probe_radius; // max distance for an interaction
      mmdb::realtype min_dist = 0.01;

      for (int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {

         // std::cout << "make_all_atom_overlaps() imod " << imod << std::endl;

         mmdb::Model *model_p = mol->GetModel(imod);
         if (! model_p) continue;
         std::string selection_string = "/" + std::to_string(imod);
         int i_sel_hnd = mol->NewSelection(); // d
         mol->SelectAtoms (i_sel_hnd, imod, "*",
                           mmdb::ANY_RES, // starting resno, an int
                           "*", // any insertion code
                           mmdb::ANY_RES, // ending resno
                           "*", // ending insertion code
                           "*", // any residue name
                           "*", // atom name
                           "*", // elements
                           "*"  // alt loc.
                           );
         bool exclude_mc_flag = true;

         long i_contact_group = 1;
         mmdb::mat44 my_matt;
         mmdb::SymOps symm;
         for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
               my_matt[i][j] = 0.0;
         for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

         mmdb::Atom **atom_selection = 0;
         int n_selected_atoms;
         mol->GetSelIndex(i_sel_hnd, atom_selection, n_selected_atoms);
         setup_env_residue_atoms_radii(i_sel_hnd); // fills neighb_atom_radius
         mmdb::Contact *pscontact = NULL;
         int n_contacts;
         mol->SeekContacts(atom_selection, n_selected_atoms,
                           atom_selection, n_selected_atoms,
                           0.01, max_dist,
                           0, // 0: in same residue also?
                           pscontact, n_contacts,
                           0, &my_matt, i_contact_group);
         if (false) {
            std::cout << "found " << n_selected_atoms << " selected atoms" << std::endl;
            std::cout << "found " << n_contacts << " all-atom contacts" << std::endl;
         }
         if (n_contacts > 0) {
            if (pscontact) {
               overlaps.reserve(1000);
               std::map<std::string, std::vector<std::pair<std::string, std::string> > > bonded_neighbours;
               std::map<std::string, std::vector<std::vector<std::string> > > ring_list_map;
               for (int i=0; i<n_contacts; i++) {
                  const int &ii = pscontact[i].id1;
                  const int &jj = pscontact[i].id2;
                  if (ii < jj) {
                     mmdb::Atom *at_1 = atom_selection[ii];
                     mmdb::Atom *at_2 = atom_selection[jj];
                     if (clashable_alt_confs(at_1, at_2)) {

                        // also check links
                        atom_interaction_type ait =
                           bonded_angle_or_ring_related(mol, at_1, at_2, exclude_mc_flag,
                                                        &bonded_neighbours,   // updated by fn.
                                                        &ring_list_map        // updated by fn.
                                                        );
                        if (ait == CLASHABLE) {
                           if (false)
                              std::cout << "debug:: clashable: " << atom_spec_t(at_1) << " " << atom_spec_t(at_2)
                                        << std::endl;
                           double r_1 = get_vdw_radius_neighb_atom(ii);
                           double r_2 = get_vdw_radius_neighb_atom(jj);
                           clipper::Coord_orth co_at_1 = co(at_1);
                           clipper::Coord_orth co_at_2 = co(at_2);
                           double ds = (co_at_1 - co_at_2).lengthsq();
                           double d = std::sqrt(ds);
                           if (false) { // debug
                              std::string resname_1(at_1->GetResName());
                              std::string resname_2(at_2->GetResName());
                              if (resname_1 == "HOH") {
                                 if (resname_2 == "HOH") {
                                    std::cout << "HOH-HOH " << d << " " << r_1 << " " << r_2 << " " << r_1 + r_2 << std::endl;
                                 }
                              }
                           }
                           if (d < (r_1 + r_2)) {
                              if (! kludge_filter(at_1, at_2)) {
                                 double ovl = get_overlap_volume(d, r_2, r_1);
                                 atom_overlap_t ao(pscontact[i].id2, at_1, at_2, r_1 ,r_2, ovl);
                                 overlaps.push_back(ao);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
         mol->DeleteSelection(i_sel_hnd);
      }
   }

   // biggest first
   sort_overlaps();
}

void
coot::atom_overlaps_container_t::sort_overlaps() {

   std::sort(overlaps.begin(), overlaps.end(), overlap_sorter);

}


// static
bool
coot::atom_overlaps_container_t::overlap_sorter(const atom_overlap_t &ao1, const atom_overlap_t &ao2) {

   return (ao2.overlap_volume < ao1.overlap_volume);

}



// first is yes/no, second is if the H is on the ligand
//
// also allow water to be return true values
//
// (ligand_atom and env_atom variable names are not appropriate for all-atom case)
//
// static
std::pair<bool, bool>
coot::atom_overlaps_container_t::is_h_bond_H_and_acceptor(mmdb::Atom *ligand_atom,
                                                          mmdb::Atom *env_atom,
                                                          int udd_h_bond_type_handle
                                                          ) {

   bool status = false;
   bool H_on_ligand = false;
   unsigned int n_h = 0;
   bool h_on_ligand_atom = false;
   bool h_on_env_atom    = false;

   int hb_1 = -1;
   int hb_2 = -1;
   if (ligand_atom->GetUDData(udd_h_bond_type_handle, hb_1) == mmdb::UDDATA_Ok) {
      if (env_atom->GetUDData(udd_h_bond_type_handle, hb_2) == mmdb::UDDATA_Ok) {

         if (false) // testing
            std::cout << "    hb_1 " << hb_1 << " for " << coot::atom_spec_t(ligand_atom)
                      << " and hb_2 " << hb_2 << " for neighb-atom " << coot::atom_spec_t(env_atom)
                      << std::endl;

         if (hb_1 == HB_HYDROGEN) {
            if (hb_2 == HB_ACCEPTOR || hb_2 == HB_BOTH) {
               status = true;
               H_on_ligand = true;
            }
         }

         if (hb_1 == HB_ACCEPTOR || hb_1 == HB_BOTH) {
            if (hb_2 == HB_HYDROGEN) {
               status = true;
               H_on_ligand = false;
            }
         }
      } else {
         // it's not bad.  Some H atoms don't have hydrogen UDD.
         if (false)
            std::cout << "   bad get of uud_h_bond info for env atom " << coot::atom_spec_t(env_atom)
                      << std::endl;
      }

      if (status == false) {
         // allow HOH to H-bond
         std::string resname_1 = ligand_atom->GetResName();
         std::string resname_2 = env_atom->GetResName();
         if (resname_1 == "HOH")
            if (hb_2 == HB_ACCEPTOR || hb_2 == HB_DONOR || hb_2 == HB_BOTH || hb_2 == HB_HYDROGEN)
               status = true;
         if (resname_2 == "HOH")
            if (hb_1 == HB_ACCEPTOR || hb_1 == HB_DONOR || hb_1 == HB_BOTH || hb_1 == HB_HYDROGEN)
               status = true;
      }
   } else {
      if (false)
         std::cout << "   bad get of uud_h_bond info for ligand atom " << coot::atom_spec_t(ligand_atom)
                   << std::endl;
   }
   return std::pair<bool, bool> (status, H_on_ligand);
}

std::string
coot::atom_overlaps_container_t::h_bond_info_t::format() const {

   std::string s;
   s += "is_H-bond-H_and_acceptor: ";
   s += is_h_bond_H_and_acceptor ? "true" : "false";
   s += " is_h_bond_donor_and_acceptor: ";
   s += is_h_bond_donor_and_acceptor ? "true" : "false";
   s += " H_is_first_atom_flag ";
   s += H_is_first_atom_flag ?  "true" : "false";
   s += " H_is_second_atom_flag ";
   s += H_is_second_atom_flag ?  "true" : "false";
   s += " donor_is_second_atom_flag: ";
   s += donor_is_second_atom_flag ?  "true" : "false";

   return s;
}


coot::atom_overlaps_container_t::h_bond_info_t::h_bond_info_t(mmdb::Atom *ligand_atom,
                                                              mmdb::Atom *env_atom,
                                                              int udd_h_bond_type_handle) {

   is_h_bond_H_and_acceptor = false;
   is_h_bond_donor_and_acceptor = false;
   H_is_first_atom_flag = false;
   H_is_second_atom_flag = false;
   donor_is_second_atom_flag = false;

   int hb_1 = -1;
   int hb_2 = -1;
   if (ligand_atom->GetUDData(udd_h_bond_type_handle, hb_1) == mmdb::UDDATA_Ok) {
      if (env_atom->GetUDData(udd_h_bond_type_handle, hb_2) == mmdb::UDDATA_Ok) {

         if (false) // 20240325-PE used to help debug the error message below (not an error, I now think)
            // std::cout << "INFO::  h_bond_info_t() env_atom->GetUDData(udd_h_bond_type_handle, hb_2) worked with handle "
            //           << udd_h_bond_type_handle << " ligand-atom: " << coot::atom_spec_t(ligand_atom) << " env_atom: "
            //           << env_atom << " " << coot::atom_spec_t(env_atom) << std::endl;
            logger.log(log_t::INFO, "h_bond_info_t() env_atom->GetUDData(udd_h_bond_type_handle, hb_2) worked with handle",
                       udd_h_bond_type_handle, "ligand-atom:", coot::atom_spec_t(ligand_atom).format(),
                       "env_atom:", coot::atom_spec_t(env_atom).format());
#if 0
        if (ligand_atom->GetSeqNum() == 901 || env_atom->GetSeqNum() == 901) {
          std::cout << "h_bond_info_t constructor " << atom_spec_t(ligand_atom)   << " "
                    << atom_spec_t(env_atom) << " " << hb_1 << " " << hb_2 << "\n";
        }
#endif

         if (hb_1 == HB_HYDROGEN) {
            if (hb_2 == HB_ACCEPTOR || hb_2 == HB_BOTH) {
               is_h_bond_H_and_acceptor = true;
               H_is_first_atom_flag = true;
            }
         }

         if (hb_1 == HB_ACCEPTOR || hb_1 == HB_BOTH) {
            if (hb_2 == HB_HYDROGEN) {
               if (false) {
                  std::cout << "in h_bond_info_t() constructor " << atom_spec_t(ligand_atom) << " " << atom_spec_t(env_atom)
                            << " --E-- " << hb_1 << " " << hb_2 << std::endl;
                  std::cout << "in h_bond_info_t() constructor found ACCEPTOR and H-atom " << coot::atom_spec_t(ligand_atom)
                            << " " << coot::atom_spec_t(env_atom) << std::endl;
               }
               is_h_bond_H_and_acceptor = true;
               H_is_second_atom_flag = true;
            }
         }

         if (hb_1 == HB_DONOR || hb_1 == HB_BOTH) {
            if (hb_2 == HB_ACCEPTOR || hb_2 == HB_BOTH) {
               is_h_bond_donor_and_acceptor = true;
            }
         }

         if (hb_1 == HB_ACCEPTOR || hb_1 == HB_BOTH) {
            if (hb_2 == HB_DONOR || hb_2 == HB_BOTH) {
               is_h_bond_donor_and_acceptor = true;
            }
         }
      } else {
         // something to do with failure to add udd_h_bond_type_handle for env_atoms from multiple MODELs?

         // 20240325-PE I no longer think this is an error - it's just that the UDData was not set for this atom.
         if (false)
            std::cout << "ERROR:: h_bond_info_t() env_atom->GetUDData(udd_h_bond_type_handle, hb_2) failed with handle "
                      << udd_h_bond_type_handle << " ligand-atom: " << coot::atom_spec_t(ligand_atom) << " env_atom: "
                      << env_atom << " " << coot::atom_spec_t(env_atom) << std::endl;
      }

      if (is_h_bond_donor_and_acceptor == false) {
         // as before, allow HOH to H-bond
         std::string resname_1 = ligand_atom->GetResName();
         std::string resname_2 = env_atom->GetResName();
         if (resname_1 == "HOH") {
            if (hb_2 == HB_ACCEPTOR || hb_2 == HB_DONOR || hb_2 == HB_BOTH || hb_2 == HB_HYDROGEN) {
               is_h_bond_H_and_acceptor = true;
            }
         }
         if (resname_2 == "HOH") {
            if (hb_1 == HB_ACCEPTOR || hb_1 == HB_DONOR || hb_1 == HB_BOTH || hb_1 == HB_HYDROGEN){
               is_h_bond_H_and_acceptor = true;
            }
         }
      }
   }
}


// in A^3
double
coot::atom_overlaps_container_t::get_overlap_volume(const double &d, const double &r_1, const double &r_2) const {

   // V = π /(12d) * (r1+r2-d) ^2 * (d^2+2d(r1+r2)-3*(r1-r2)^2)

   double V = (M_PI/(12*d)) * (r_1+r_2-d) * (r_1+r_2-d) * (d*d + 2.0*d*(r_1+r_2) - 3*(r_1-r_2)*(r_1-r_2));
   return V;
}


// this only works for ligand contacts (not all-atom)
//
double
coot::atom_overlaps_container_t::get_vdw_radius_ligand_atom(mmdb::Atom *at) {

   double r = 2.5;

   std::map<mmdb::Atom *, double>::const_iterator it = central_residue_atoms_vdw_radius_map.find(at);
   if (it == central_residue_atoms_vdw_radius_map.end()) {
      // we need to add it then

      // What's the energy type of Atom at?
      //
      // PDBv3 FIXME - change from 4-char
      std::string te = central_residue_dictionary.type_energy(at->GetAtomName());
      if (! te.empty()) {
         std::map<std::string, double>::const_iterator it_type =
            type_to_vdw_radius_map.find(te);
         if (it_type == type_to_vdw_radius_map.end()) {
            // didn't find it. so look it up and add it.
            if (geom_p)
               r = type_energy_to_radius(te);
            if (false)
               std::cout << "setting ligand atom map: type_to_vdw_radius_h_bond_type_map["
                         << std::setw(4) << te << "] to " << r << std::endl;
            type_to_vdw_radius_map[te] = r;
         } else {
            r = it_type->second;
         }
         central_residue_atoms_vdw_radius_map[at] = r;
      } else {
         std::cout << "failed to find type-energy for atom " << atom_spec_t(at) << std::endl;
      }
   } else {
      if (false)
         std::cout << "radius for atom " << atom_spec_t(at) << " was found in map: value: "
                   << it->second << std::endl;
      r = it->second;
   }

   return r;
}


double
coot::atom_overlaps_container_t::get_vdw_radius_neighb_atom(int idx_neigh_atom) const {

   // no index checking (a bit cowboy?)
   //
   double r = neighb_atom_radius[idx_neigh_atom];
   return r;

}

double
coot::atom_overlaps_container_t::get_vdw_radius_neighb_atom(mmdb::Atom *at, unsigned int idx_res) {

   double r = 1.5;

   std::map<mmdb::Atom *, double>::const_iterator it = neighbour_atoms_vdw_radius_map.find(at);
   if (it == neighbour_atoms_vdw_radius_map.end()) {
      // we need to add it then

      // What's the energy type of Atom at?
      //
      std::string te = neighb_dictionaries[idx_res].type_energy(at->GetAtomName());
      std::map<std::string, double>::const_iterator it_type = type_to_vdw_radius_map.find(te);
      if (it_type == type_to_vdw_radius_map.end()) {
         // didn't find te in types map. so look it up from the dictionary and add to the types map
         r = geom_p->get_energy_lib_atom(te).vdw_radius;
         type_to_vdw_radius_map[te] = r;
      } else {
         r = it_type->second;
      }
      neighbour_atoms_vdw_radius_map[at] = r;
   } else {
      r = it->second;
   }

   return r;
}



coot::hb_t
coot::atom_overlaps_container_t::get_h_bond_type(mmdb::Atom *at) {

   hb_t type = HB_UNASSIGNED;
   std::string atom_name = at->name;
   std::string res_name = at->GetResName();
   type = geom_p->get_h_bond_type(atom_name, res_name, protein_geometry::IMOL_ENC_ANY); // heavyweight

   return type;
}


#include "fib-sphere.hh"

// probably not the function that you're looking for.
// Maybe delete this?
void
coot::atom_overlaps_container_t::contact_dots_for_overlaps() const {

   double spike_length = 0.5;
   double clash_dist = 0.4;

   std::vector<clipper::Coord_orth> sphere_points = fibonacci_sphere(450); // dot density
   std::vector<clipper::Coord_orth> H_sphere_points = fibonacci_sphere(270); // less than above

   for (unsigned int i=0; i<overlaps.size(); i++) {

      // std::cout << "considering overlap idx: " << i << std::endl;

      clipper::Coord_orth pt_at_1(overlaps[i].atom_1->x,
                                  overlaps[i].atom_1->y,
                                  overlaps[i].atom_1->z);
      clipper::Coord_orth pt_at_2(overlaps[i].atom_2->x,
                                  overlaps[i].atom_2->y,
                                  overlaps[i].atom_2->z);
      const double &r_1 = overlaps[i].r_1;
      const double &r_2 = overlaps[i].r_2;
      const int &idx =    overlaps[i].ligand_atom_index;
      double r_2_sqrd = r_2 * r_2;
      double r_2_plus_prb_squard  = r_2_sqrd + 2 * r_2 * probe_radius + probe_radius * probe_radius;
      double r_2_minux_prb_squard = r_2_sqrd - 2 * r_2 * probe_radius + probe_radius * probe_radius;
      //
      bool done_atom_name = false;

      for (unsigned int j=0; j<sphere_points.size(); j++) {

         const clipper::Coord_orth &pt(sphere_points[j]);
         clipper::Coord_orth pt_at_surface = r_1 * pt + pt_at_1;
         double d_sqrd = (pt_at_2 - pt_at_surface).lengthsq();

         // std::cout << "comparing " << sqrt(d_sqrd) << " vs " << sqrt(r_2_plus_prb_squard)
         // << " with r_2 " << r_2 << std::endl;

         if (d_sqrd > r_2_plus_prb_squard) {

            bool draw_it = ! is_inside_another_ligand_atom(idx, pt_at_surface);

            if (false) // debugging
               if (std::string(overlaps[i].atom_1->name) != " HO3")
                  draw_it = false;

            if (draw_it) {

               std::string type = "wide-contact";
               bool only_pt = true; // not a spike

               if (d_sqrd < r_2_sqrd)
                  type = "close-contact";

               if (d_sqrd < (r_2_sqrd - 2 * r_2 * clash_dist + clash_dist * clash_dist)) {
                  type = "clash";
                  only_pt = false;
               }

               if (overlaps[i].is_h_bond) {
                  type = "H-bond";
               }

               if (false) // debugging
                  if (! done_atom_name) {
                     std::cout << "   spike for atom " << coot::atom_spec_t(overlaps[i].atom_1)
                               << std::endl;
                     done_atom_name = true;
                  }

               clipper::Coord_orth pt_spike_inner = pt_at_surface;
               if (! only_pt) {
                  std::cout << "considering overlap idx: " << i << " "
                            << atom_spec_t(overlaps[i].atom_1) << " to "
                            << atom_spec_t(overlaps[i].atom_2) << std::endl;

                  clipper::Coord_orth vect_to_pt_1 = pt_at_1 - pt_at_surface;
                  clipper::Coord_orth vect_to_pt_1_unit(vect_to_pt_1.unit());

                  // these days, spikes project away from the atom, not inwards
                  //
                  // clipper::Coord_orth pt_spike_inner = pt_at_surface + spike_length * vect_to_pt_1_unit;
                  pt_spike_inner = pt_at_surface - spike_length * vect_to_pt_1_unit;
               }

               // on the surface of atom_1 inside the sphere of atom_2
               if (false)
                  std::cout << "spike "
                            << type << " "
                            << pt_at_surface.x() << " "
                            << pt_at_surface.y() << " "
                            << pt_at_surface.z() << " to "
                            << pt_spike_inner.x() << " "
                            << pt_spike_inner.y() << " "
                            << pt_spike_inner.z()
                            << std::endl;
            }
         }
      }
   }
}

bool
coot::atom_overlaps_container_t::is_inside_another_ligand_atom(int idx,
                                                               const clipper::Coord_orth &dot_pt) const {
   bool r = false;

   if (idx >= 0) {

      std::map<int, std::vector<std::pair<mmdb::Atom *, double> > >::const_iterator it;
      it = ligand_atom_neighbour_map.find(idx);

      if (it != ligand_atom_neighbour_map.end()) {

         const std::vector<std::pair<mmdb::Atom *, double> > &v = it->second;
         for (unsigned int i=0; i<v.size(); i++) {
            clipper::Coord_orth pt = co(v[i].first);
            double dist_sqrd = (dot_pt - pt).lengthsq();

            const double &radius_other = v[i].second;
            if (dist_sqrd < radius_other * radius_other) {
               r = true;
               break;
            }
         }
      } else {
         std::cout << "Opps! Missing in ligand_atom_neighbour_map: idx " << idx << std::endl;
      }
   }
   return r;
}
bool
coot::atom_overlaps_container_t::is_inside_another_ligand_atom(int idx,
                                                               const clipper::Coord_orth &probe_pos,
                                                               const clipper::Coord_orth &dot_pt) const {

   // for better speed, change the atom -> radius map to atom index -> radius vector
   // (and directly look up the radius).

   bool consider_probe_also = false; // if this is true, then don't make surface points for inner cusps
   bool r = false;

   if (idx >= 0) {

      const std::vector<std::pair<mmdb::Atom *, double> > &v = ligand_atom_neighbour_map.find(idx)->second;
      for (unsigned int i=0; i<v.size(); i++) {
         clipper::Coord_orth pt = co(v[i].first);
         double dist_sqrd = (dot_pt - pt).lengthsq(); // should be r_1 squared, right?
         double radius_other = v[i].second;
         radius_other += probe_radius;
         if (dist_sqrd < radius_other * radius_other) {
            r = true;
            break;
         }
      }
   }
   return r;
}


coot::atom_overlaps_dots_container_t
coot::atom_overlaps_container_t::contact_dots_for_ligand(double dot_density_in) { // or residue

   atom_overlaps_dots_container_t ao;
   if (!have_dictionary) {
      std::cout << "WARNING:: contact_dots_for_ligand() no dictionary " << std::endl;
      return ao;
   }

   bool add_vdw_dots = true; // pass this

   mmdb::realtype max_dist = 4.0; // max distance for an interaction

   bool excl_mc_flag = true; // exclude main-chain to main-chain interactions also

   if (mol) {
      mmdb::Contact *pscontact = NULL;
      int n_contacts;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
         for (int j=0; j<4; j++)
            my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mmdb::Atom **ligand_residue_atoms = 0;
      int n_ligand_residue_atoms = 0;
      res_central->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);

      std::vector<residue_spec_t> env_residue_specs(neighbours.size());
      for (unsigned int i=0; i<neighbours.size(); i++)
         env_residue_specs[i] = residue_spec_t(neighbours[i]);
      add_residue_neighbour_index_to_neighbour_atoms();

      int mask_mode = 0; // all atoms
      int i_sel_hnd_env_atoms = specs_to_atom_selection(env_residue_specs, mol, mask_mode);

      mmdb::Atom **env_residue_atoms = 0;
      int n_env_residue_atoms = 0;
      mol->GetSelIndex(i_sel_hnd_env_atoms, env_residue_atoms, n_env_residue_atoms);
      setup_env_residue_atoms_radii(i_sel_hnd_env_atoms);

      mol->SeekContacts(ligand_residue_atoms, n_ligand_residue_atoms,
                        env_residue_atoms, n_env_residue_atoms,
                        0, max_dist,
                        1, // 0: in same residue also?
                        pscontact, n_contacts,
                        0, &my_matt, i_contact_group);

      if (n_contacts > 0) {
         if (pscontact) {

            int   atom_n_sphere_dots = static_cast<int>(450 * dot_density_in);
            int H_atom_n_sphere_dots = static_cast<int>(200 * dot_density_in);  // less than above
            std::vector<clipper::Coord_orth> sphere_points = fibonacci_sphere(atom_n_sphere_dots);
            std::vector<clipper::Coord_orth> H_sphere_points = fibonacci_sphere(H_atom_n_sphere_dots);

            // which atoms are close to which other atoms?
            std::map<int, std::vector<int> > contact_map; // these atoms can have nbc interactions
            std::map<int, std::vector<int> > bonded_map;  // these atoms are bonded and can mask

            // the dots of an atom
            // which atom names of which residues are bonded or 1-3 related? (update the map
            // as you find and add new residue types in bonded_angle_or_ring_related().
            std::map<std::string, std::vector<std::pair<std::string, std::string> > > bonded_neighbours;
            // similar thinking: update the ring list map
            std::map<std::string, std::vector<std::vector<std::string> > > ring_list_map;

            for (int i=0; i<n_contacts; i++) {
               atom_interaction_type ait =
                  bonded_angle_or_ring_related(mol, // also check links
                                               ligand_residue_atoms[pscontact[i].id1],
                                               env_residue_atoms[pscontact[i].id2], excl_mc_flag,
                                               &bonded_neighbours,   // updatedby fn.
                                               &ring_list_map        // updatedby fn.
                                               );
               if (ait == CLASHABLE) {
                  contact_map[pscontact[i].id1].push_back(pscontact[i].id2);
               } else {
                  if (ait == BONDED) {
                     bonded_map[pscontact[i].id1].push_back(pscontact[i].id2);
                  }
               }
            }
            // std::cout << "done contact map" << std::endl;


            for (int iat=0; iat<n_ligand_residue_atoms; iat++) {

               mmdb::Atom *cr_at = ligand_residue_atoms[iat];
               clipper::Coord_orth pt_at_1 = co(cr_at);

               double r_1 = get_vdw_radius_ligand_atom(cr_at);

               std::vector<clipper::Coord_orth> &sphere_points_for_atom = sphere_points;
               if (std::string(cr_at->element) == " H")
                  sphere_points_for_atom = H_sphere_points;

               for (unsigned int j=0; j<sphere_points_for_atom.size(); j++) {

                  const clipper::Coord_orth &pt(sphere_points[j]);
                  clipper::Coord_orth pt_at_surface = r_1 * pt + pt_at_1;
                  bool draw_it = ! is_inside_another_ligand_atom(iat, pt_at_surface);

                  bool draw_it_2 = ! is_inside_an_env_atom_to_which_its_bonded(iat,
                                                                               bonded_map[iat],
                                                                               env_residue_atoms,
                                                                               pt_at_surface);

                  if (draw_it && draw_it_2) {

                     // is the point on the surface of cr_at inside the sphere
                     // of an environment atom?
                     // If so, we want to know which one to which it was closest

                     double biggest_overlap = -1; // should be positive if we get a hit
                     mmdb::Atom *atom_with_biggest_overlap = 0;
                     double r_2_for_biggest_overlap = 0;

                     // now check which atom this is clashing with (if any) and pick the
                     // one with the biggest overlap
                     //
                     const std::vector<int> &v = contact_map[iat];
                     for (unsigned int jj=0; jj<v.size(); jj++) {

                        mmdb::Atom *neighb_atom = env_residue_atoms[v[jj]];

                        // turn off main-chain to main-chain interactions (to match the dots from probe)
                        if (is_main_chain_p(cr_at)) // this might be true for amino-acid residues
                           if (is_main_chain_p(neighb_atom))
                              continue;

                        double r_2 = get_vdw_radius_neighb_atom(v[jj]);
                        double r_2_sqrd = r_2 * r_2;

                        // note that we want to check against the outside of the probe
                        // sphere - not the middle of the probe sphere - and
                        // thus we will pick up more wide contacts (like molprobity probe)

                        double rmp = 1.6;  // radius multiplier
                        double r_2_plus_prb_squard =
                           r_2_sqrd + 2 * r_2 * rmp * probe_radius + rmp * rmp * probe_radius * probe_radius;

                        clipper::Coord_orth pt_na = co(neighb_atom);
                        double d_sqrd = (pt_na - pt_at_surface).lengthsq();

                        if (false)
                           std::cout << " for atom "
                                     << atom_spec_t(env_residue_atoms[v[jj]])
                                     << " comparing " << d_sqrd << " vs "
                                     << r_2_plus_prb_squard << " with r_2 " << r_2
                                     << " probe_radius " << probe_radius << std::endl;

                        if (d_sqrd < r_2_plus_prb_squard) {

                           // OK it was close to something.
                           double delta_d_sqrd = r_2_plus_prb_squard - d_sqrd;
                           if (delta_d_sqrd > biggest_overlap) {
                              biggest_overlap = delta_d_sqrd;
                              atom_with_biggest_overlap = neighb_atom;
                              r_2_for_biggest_overlap = r_2;
                           }
                        }
                     }

                     if (atom_with_biggest_overlap) {
                        double d_surface_pt_to_atom_sqrd =
                           (co(atom_with_biggest_overlap)-pt_at_surface).lengthsq();
                        double d_surface_pt_to_atom = sqrt(d_surface_pt_to_atom_sqrd);
                        double overlap_delta = r_2_for_biggest_overlap - d_surface_pt_to_atom;

                        // first is yes/no, second is H-is-on-ligand?
                        // std::pair<bool, bool> might_be_h_bond_flag =
                        // is_h_bond_H_and_acceptor(cr_at, atom_with_biggest_overlap, udd_h_bond_type_handle);
                        h_bond_info_t hbi(cr_at, atom_with_biggest_overlap, udd_h_bond_type_handle);
                        bool is_h_bond = false;
                        if (hbi.is_h_bond_H_and_acceptor) is_h_bond = true;
                        if (hbi.is_h_bond_donor_and_acceptor)  is_h_bond = true;

                        std::pair<std::string, std::string> c_type_col =
                           overlap_delta_to_contact_type(overlap_delta, is_h_bond);
                        const std::string &c_type = c_type_col.first;
                        const std::string &col    = c_type_col.second;

                        if (false)
                           std::cout << "............ here in contact_dots_for_ligand() "
                                     << coot::atom_spec_t(cr_at) << " "
                                     << coot::atom_spec_t(atom_with_biggest_overlap) << " "
                                     << "with is_h_bond " << is_h_bond
                                     << " " <<  c_type << " " << col << std::endl;

                        if (false)
                           std::cout << "spike check "
                                     << c_type << " "
                                     << pt_at_surface.x() << " "
                                     << pt_at_surface.y() << " "
                                     << pt_at_surface.z() << " to "
                                     << pt_at_surface.x() << " "
                                     << pt_at_surface.y() << " "
                                     << pt_at_surface.z()
                                     << std::endl;

                        if (c_type != "clash") {
                           atom_overlaps_dots_container_t::dot_t dot(overlap_delta, col, pt_at_surface);
                           ao.dots[c_type].push_back(dot);
                        } else {
                           clipper::Coord_orth vect_to_pt_1 = pt_at_1 - pt_at_surface;
                           clipper::Coord_orth vect_to_pt_1_unit(vect_to_pt_1.unit());
                           // these days, spikes project away from the atom, not inwards
                           // clipper::Coord_orth pt_spike_inner =
                           // pt_at_surface + clash_spike_length * vect_to_pt_1_unit;
                           clipper::Coord_orth pt_spike_inner =
                              pt_at_surface -
                              0.34* sqrt(biggest_overlap) * clash_spike_length * vect_to_pt_1_unit;
                           std::pair<clipper::Coord_orth, clipper::Coord_orth> p(pt_at_surface,
                                                                                 pt_spike_inner);
                           ao.clashes.positions.push_back(p);
                        }

                     } else {

                        // no environment atom was close to this ligand atom, so just add
                        // a surface point

                        if (false)
                           std::cout << "spike-surface "
                                     << pt_at_surface.x() << " "
                                     << pt_at_surface.y() << " "
                                     << pt_at_surface.z() << std::endl;

                        if (add_vdw_dots) {
                           atom_overlaps_dots_container_t::dot_t dot(0, "grey", pt_at_surface);
                           ao.dots["vdw-surface"].push_back(dot);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return ao;
}

// dot_density_in is default 0.5.
// bool make_vdw_surface is default false.
//
coot::atom_overlaps_dots_container_t
coot::atom_overlaps_container_t::all_atom_contact_dots(double dot_density_in,
                                                       bool make_vdw_surface) {

   atom_overlaps_dots_container_t ao;

   if (mol) {
      mmdb::realtype max_dist = 1.75 + 1.75 + 2.0 * probe_radius; // max distance for an interaction
      mmdb::realtype min_dist = 0.01;
      for (int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {

         int i_sel_hnd = mol->NewSelection(); // d
         mol->SelectAtoms (i_sel_hnd, imod, "*",
                           mmdb::ANY_RES, // starting resno, an int
                           "*", // any insertion code
                           mmdb::ANY_RES, // ending resno
                           "*", // ending insertion code
                           "*", // any residue name
                           "*", // atom name
                           "*", // elements
                           "*"  // alt loc.
                           );

         unsigned int n_threads = get_max_number_of_threads();

         n_threads = 1;

         if (n_threads == 0) {
            ao = all_atom_contact_dots_internal_single_thread(dot_density_in, mol, i_sel_hnd, i_sel_hnd,
                                                              min_dist, max_dist, make_vdw_surface);
         } else {
            ao = all_atom_contact_dots_internal_multi_thread(dot_density_in, mol, i_sel_hnd, i_sel_hnd,
                                                             min_dist, max_dist, make_vdw_surface);
         }

#if 0 // single thread
      // set ao using non-threaded version
      //
         ao = all_atom_contact_dots_internal_single_thread(dot_density_in, mol, i_sel_hnd, i_sel_hnd,
                                                           min_dist, max_dist, make_vdw_surface);

#endif //

         mol->DeleteSelection(i_sel_hnd);
      }
   }
   return ao;
}



coot::atom_overlaps_dots_container_t
coot::atom_overlaps_container_t::all_atom_contact_dots_internal_single_thread(double dot_density_in,
                                                                              mmdb::Manager *mol,
                                                                              int i_sel_hnd_1,
                                                                              int i_sel_hnd_2,
                                                                              mmdb::realtype min_dist,
                                                                              mmdb::realtype max_dist,
                                                                              bool make_vdw_surface) {

   bool exclude_mc_flag = true; // pass this?

   coot::atom_overlaps_dots_container_t ao;

   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   mmdb::Atom **atom_selection = 0;
   int n_selected_atoms;
   mol->GetSelIndex(i_sel_hnd_1, atom_selection, n_selected_atoms);
   setup_env_residue_atoms_radii(i_sel_hnd_1);
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   mol->SeekContacts(atom_selection, n_selected_atoms,
                     atom_selection, n_selected_atoms,
                     0.01, max_dist,
                     0, // 0: in same residue also?
                     pscontact, n_contacts,
                     0, &my_matt, i_contact_group);
   if (false) {
      std::cout << "found " << n_selected_atoms << " selected atoms" << std::endl;
      std::cout << "found " << n_contacts << " all-atom contacts" << std::endl;
   }

   int n_sphere_points        = static_cast<int>(dot_density_in * 900);
   int n_sphere_points_for_H  = static_cast<int>(dot_density_in * 440);
   std::vector<clipper::Coord_orth>   sphere_points = fibonacci_sphere(n_sphere_points);
   std::vector<clipper::Coord_orth> H_sphere_points = fibonacci_sphere(n_sphere_points_for_H);

   if (n_contacts > 0) {
      if (pscontact) {
         // which atoms are close to which other atoms?
         std::map<int, std::vector<int> > contact_map; // these atoms can have nbc interactions
         std::map<int, std::vector<int> > bonded_map;  // these atoms are bonded and can mask
         // the dots of an atom
         // which atom names of which residues are bonded or 1-3 related? (update the map
         // as you find and add new residue types in bonded_angle_or_ring_related().
         std::map<std::string, std::vector<std::pair<std::string, std::string> > > bonded_neighbours;
         // similar thinking: update the ring list map
         std::map<std::string, std::vector<std::vector<std::string> > > ring_list_map;

         // initialize contact_map and bonded_map - not sure if this speeds things up - but it
         // doesn't seem to slow things down.
         //
         std::vector<int> dum;
         for (int iat=0; iat<n_selected_atoms; iat++)
            contact_map[iat] = dum;
         for (int iat=0; iat<n_selected_atoms; iat++)
            bonded_map[iat] = dum;

         for (int i=0; i<n_contacts; i++) {
            if (pscontact[i].id1 < pscontact[i].id2) {
               if (clashable_alt_confs(atom_selection[pscontact[i].id1], atom_selection[pscontact[i].id2])) {
                  atom_interaction_type ait =
                     bonded_angle_or_ring_related(mol, // also check links
                                                  atom_selection[pscontact[i].id1],
                                                  atom_selection[pscontact[i].id2], exclude_mc_flag,
                                                  &bonded_neighbours,   // updatedby fn.
                                                  &ring_list_map        // updatedby fn.
                                                  );
                  if (ait == CLASHABLE) {
                     contact_map[pscontact[i].id1].push_back(pscontact[i].id2);
                     contact_map[pscontact[i].id2].push_back(pscontact[i].id1);
                  } else {
                     if (ait == BONDED) {
                        bonded_map[pscontact[i].id1].push_back(pscontact[i].id2);
                        bonded_map[pscontact[i].id2].push_back(pscontact[i].id1);
                     }
                  }
               }
            }
         }

         for (int iat=0; iat<n_selected_atoms; iat++) {

            mmdb::Atom *at = atom_selection[iat];

            // if (!(atom_spec_t(at) == debug_spec))
            // continue;
            // if ((atom_spec_t(at) == debug_spec) || (atom_spec_t(at) == debug_spec_2)) {
            //                   // don't stop this atom
            // } else {
            // continue;
            // }

            clipper::Coord_orth pt_at_1 = co(at);
            double dot_density = dot_density_in;
            double r_1 = get_vdw_radius_neighb_atom(iat);

            std::vector<clipper::Coord_orth> &sphere_points_for_atom = sphere_points;
            if (std::string(at->element) == " H")
               sphere_points_for_atom = H_sphere_points;

            // if (std::string(at->element) == " H")
               // dot_density *=0.66; // so that surface dots on H atoms don't appear (weirdly) more fine
            if (std::string(at->element) == " H")
               sphere_points_for_atom = H_sphere_points;

            for (unsigned int j=0; j<sphere_points_for_atom.size(); j++) {
               const clipper::Coord_orth &pt(sphere_points_for_atom[j]);
               clipper::Coord_orth pt_at_surface = r_1 * pt + pt_at_1;

               if (false)
                  std::cout << "pt_at_surface "
                            << pt_at_surface.x() << " "
                            << pt_at_surface.y() << " "
                            << pt_at_surface.z() << "\n";
               bool draw_it = ! is_inside_another_atom_to_which_its_bonded(iat, at,
                                                                           pt_at_surface,
                                                                           bonded_map[iat],
                                                                           atom_selection,
                                                                           neighb_atom_radius);

               if (draw_it) {

                  double biggest_overlap = -1; // should be positive if we get a hit
                  mmdb::Atom *atom_with_biggest_overlap = 0;
                  double r_2_for_biggest_overlap = 0;

                  // now check which atom this is clashing with (if any) and pick the
                  // one with the biggest overlap
                  //
                  const std::vector<int> &v = contact_map[iat];
                  for (unsigned int jj=0; jj<v.size(); jj++) {
                     mmdb::Atom *neighb_atom = atom_selection[v[jj]];
                     double r_2 = get_vdw_radius_neighb_atom(v[jj]);
                     double r_2_sqrd = r_2 * r_2;
                     double r_2_plus_prb_squard = r_2_sqrd + 2 * r_2 * probe_radius +
                        4 * probe_radius * probe_radius;
                     clipper::Coord_orth pt_na = co(neighb_atom);
                     double d_sqrd = (pt_na - pt_at_surface).lengthsq();

                     if (false)
                        std::cout << " for " << atom_spec_t(at) << " " << atom_spec_t(neighb_atom)
                                  << " comparing " << d_sqrd << " vs " << r_2_plus_prb_squard
                                  << std::endl;
                     if (d_sqrd < r_2_plus_prb_squard) {
                        // a contact dot on something
                        double delta_d_sqrd = r_2_plus_prb_squard - d_sqrd;
                        if (delta_d_sqrd > biggest_overlap) {
                           biggest_overlap = delta_d_sqrd;
                           atom_with_biggest_overlap = neighb_atom;
                           r_2_for_biggest_overlap = r_2;
                        }
                     }
                  }

                  if (atom_with_biggest_overlap) {
                     double d_surface_pt_to_atom_sqrd =
                        (co(atom_with_biggest_overlap) - pt_at_surface).lengthsq();
                     double d_surface_pt_to_atom = sqrt(d_surface_pt_to_atom_sqrd);
                     double overlap_delta = r_2_for_biggest_overlap - d_surface_pt_to_atom;
                     // first is yes/no, second is H-is-on-ligand?
                     // allow waters to H-bond (without being an H)
                     // std::pair<bool, bool> might_be_h_bond_flag =
                     // is_h_bond_H_and_acceptor(at, atom_with_biggest_overlap, udd_h_bond_type_handle);
                     h_bond_info_t hbi(at, atom_with_biggest_overlap, udd_h_bond_type_handle);
                     bool is_h_bond = false;
                     if (hbi.is_h_bond_H_and_acceptor)
                        is_h_bond = true;

                     std::pair<std::string, std::string> c_type_col =
                        overlap_delta_to_contact_type(overlap_delta, is_h_bond);
                     const std::string &c_type = c_type_col.first;
                     const std::string &col    = c_type_col.second;

                     clipper::Coord_orth pt_spike_inner = pt_at_surface;
                     if (c_type == "clash") {
                        clipper::Coord_orth vect_to_pt_1 = pt_at_1 - pt_at_surface;
                        clipper::Coord_orth vect_to_pt_1_unit(vect_to_pt_1.unit());
                        pt_spike_inner = pt_at_surface -
                           0.3 * sqrt(biggest_overlap) * clash_spike_length * vect_to_pt_1_unit;
                     }

                     if (c_type != "clash") {

                        // draw dot if these are atoms from different residues or this is not
                        // a wide contact
                        if (at->residue != atom_with_biggest_overlap->residue) {
                           atom_overlaps_dots_container_t::dot_t dot(overlap_delta, col, pt_at_surface);
                           ao.dots[c_type].push_back(dot);
                        }
                     } else {
                        // clash
                        clipper::Coord_orth vect_to_pt_1 = pt_at_1 - pt_at_surface;
                        clipper::Coord_orth vect_to_pt_1_unit(vect_to_pt_1.unit());
                        // these days, spikes project away from the atom, not inwards
                        // clipper::Coord_orth pt_spike_inner =
                        // pt_at_surface + clash_spike_length * vect_to_pt_1_unit;
                        clipper::Coord_orth pt_spike_inner =
                           pt_at_surface -
                           0.34* sqrt(biggest_overlap) * clash_spike_length * vect_to_pt_1_unit;
                        std::pair<clipper::Coord_orth, clipper::Coord_orth> p(pt_at_surface,
                                                                              pt_spike_inner);
                        ao.clashes.positions.push_back(p);
                     }

                  } else {

                     // no environment atom was close to this ligand atom, so just add
                     // a surface point

                     if (make_vdw_surface) {
                        atom_overlaps_dots_container_t::dot_t dot(0, "grey", pt_at_surface);
                        ao.dots["vdw-surface"].push_back(dot);
                     }
                  }
               }
            }
         }
      }
   }
   return ao;
}




coot::atom_overlaps_dots_container_t
coot::atom_overlaps_container_t::all_atom_contact_dots_internal_multi_thread(double dot_density_in,
                                                                             mmdb::Manager *mol,
                                                                             int i_sel_hnd_1,
                                                                             int i_sel_hnd_2,
                                                                             mmdb::realtype min_dist,
                                                                             mmdb::realtype max_dist,
                                                                             bool make_vdw_surface) {

   coot::atom_overlaps_dots_container_t ao;

   bool exclude_mc_flag = true;

   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   mmdb::Atom **atom_selection = 0;
   int n_selected_atoms = 0;
   mol->GetSelIndex(i_sel_hnd_1, atom_selection, n_selected_atoms);
   setup_env_residue_atoms_radii(i_sel_hnd_1);
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   mol->SeekContacts(atom_selection, n_selected_atoms,
                     atom_selection, n_selected_atoms,
                     0.01, max_dist,
                     0, // 0: in same residue also?
                     pscontact, n_contacts,
                     0, &my_matt, i_contact_group);
   if (false) {
      std::cout << "found " << n_selected_atoms << " selected atoms" << std::endl;
      std::cout << "found " << n_contacts << " all-atom contacts" << std::endl;
   }
   if (n_contacts > 0) {
      if (pscontact) {
         // which atoms are close to which other atoms?
         std::map<int, std::vector<int> > contact_map; // these atoms can have nbc interactions
         std::map<int, std::vector<int> > bonded_map;  // these atoms are bonded and can mask
         // the dots of an atom
         // which atom names of which residues are bonded or 1-3 related? (update the map
         // as you find and add new residue types in bonded_angle_or_ring_related().
         std::map<std::string, std::vector<std::pair<std::string, std::string> > > bonded_neighbours;
         // similar thinking: update the ring list map
         std::map<std::string, std::vector<std::vector<std::string> > > ring_list_map;

         // initialize contact_map and bonded_map - not sure if this speeds things up - but it
         // doesn't seem to slow things down.
         //
         for (int iat=0; iat<n_selected_atoms; iat++)
            contact_map[iat].reserve(12);
         for (int iat=0; iat<n_selected_atoms; iat++)
            bonded_map[iat].reserve(4);

         for (int i=0; i<n_contacts; i++) {
            if (pscontact[i].id1 < pscontact[i].id2) {
               if (clashable_alt_confs(atom_selection[pscontact[i].id1], atom_selection[pscontact[i].id2])) {
                  atom_interaction_type ait =
                     bonded_angle_or_ring_related(mol, // also check links
                                                  atom_selection[pscontact[i].id1],
                                                  atom_selection[pscontact[i].id2], exclude_mc_flag,
                                                  &bonded_neighbours,   // updatedby fn.
                                                  &ring_list_map        // updatedby fn.
                                                  );
                  if (ait == CLASHABLE) {
                     contact_map[pscontact[i].id1].push_back(pscontact[i].id2);
                     contact_map[pscontact[i].id2].push_back(pscontact[i].id1);
                  } else {
                     if (ait == BONDED) {
                        bonded_map[pscontact[i].id1].push_back(pscontact[i].id2);
                        bonded_map[pscontact[i].id2].push_back(pscontact[i].id1);
                     }
                  }
               }
            }
         }


         if (false) { // show contact map:
            std::cout << "done contact_map and bonded_map " << std::endl;
            std::cout << "     contact map: " << contact_map.size() << std::endl;
            std::map<int, std::vector<int> >::const_iterator it;
            for (it=contact_map.begin(); it!=contact_map.end(); ++it) {
               for (unsigned int jj=0; jj<it->second.size(); jj++) {
                  std::cout << "   " << atom_spec_t(atom_selection[it->first]) << " contact to "
                            << atom_spec_t(atom_selection[it->second[jj]]) << std::endl;
               }
            }
            std::cout << "     bonded map: " << bonded_map.size() << std::endl;
            for (it=bonded_map.begin(); it!=bonded_map.end(); ++it) {
               for (unsigned int jj=0; jj<it->second.size(); jj++) {
                  std::cout << "   " << atom_spec_t(atom_selection[it->first]) << " bonded to "
                            << atom_spec_t(atom_selection[it->second[jj]]) << std::endl;
               }
            }
         }

         unsigned int n_threads = get_max_number_of_threads();
         n_threads = 1;
         std::vector<std::thread> threads;
         unsigned int n_per_thread = n_selected_atoms/n_threads;
         // std::cout << "n per thread " << n_per_thread << std::endl;
         std::vector<atom_overlaps_dots_container_t> results_container_vec(n_threads);

         for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {
            int iat_start = i_thread * n_per_thread;
            int iat_end = iat_start + n_per_thread;
            // for the last thread, set the end atom index
            if (i_thread == (n_threads - 1))
               iat_end = n_selected_atoms; // for loop uses iat_start and tests for < iat_end

            if (false) // useful for debugging
               std::cout << "thread: " << i_thread << " from atom index " << iat_start << " to "
                         << iat_end << std::endl;

            results_container_vec[i_thread] = atom_overlaps_dots_container_t(n_per_thread);
            // 20230814-PE I saw a crash here today. It seems to be crashing on the creation
            // of the thread - I don't know what that means.
            //
            threads.push_back(std::thread(contacts_for_atoms, iat_start, iat_end,
                                          atom_selection, contact_map, bonded_map,
                                          neighb_atom_radius, udd_h_bond_type_handle,
                                          molecule_has_hydrogens, probe_radius,
                                          dot_density_in, clash_spike_length, make_vdw_surface,
                                          &results_container_vec[i_thread]));

         }

         for (unsigned int i_thread=0; i_thread<n_threads; i_thread++)
            threads.at(i_thread).join();

         if (false) {
            for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {
               std::unordered_map<std::string, std::vector<atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
               for (it =results_container_vec[i_thread].dots.begin();
                    it!=results_container_vec[i_thread].dots.end();
                    ++it) {
                  std::cout << "thread  " << i_thread << " results size "
                            << it->first << " " << it->second.size() << std::endl;
               }
            }
         }

         for (unsigned int i_thread=0; i_thread<n_threads; i_thread++)
            ao.add(results_container_vec[i_thread]);

         if (false) { // debugging
            std::cout << "consolidated" << std::endl;
            std::unordered_map<std::string, std::vector<atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
            for (it=ao.dots.begin(); it!=ao.dots.end(); ++it)
               std::cout << " consolidated size "
                         << it->first << " " << it->second.size() << std::endl;
         }
      }
   }

   return ao;
}

// put results in ao
//
void
coot::atom_overlaps_container_t::contacts_for_atoms(int iat_start, int iat_end,
                                                    mmdb::Atom **atom_selection,
                                                    const std::map<int, std::vector<int> > &contact_map,
                                                    const std::map<int, std::vector<int> > &bonded_map,
                                                    const std::vector<double> &neighb_atom_radius,
                                                    int udd_h_bond_type_handle,
                                                    bool molecule_has_hydrogens,
                                                    double probe_radius,
                                                    double dot_density_in,
                                                    double clash_spike_length,
                                                    bool make_vdw_surface,
                                                    coot::atom_overlaps_dots_container_t *ao) {

   for (int iat=iat_start; iat<iat_end; iat++) {

      ao->add(contacts_for_atom(iat, atom_selection, contact_map, bonded_map, neighb_atom_radius,
                                udd_h_bond_type_handle, molecule_has_hydrogens,
                                probe_radius, dot_density_in, clash_spike_length, make_vdw_surface));
   }
}


float
coot::atom_overlaps_container_t::score() {
   float s = 0.0;
   unsigned int nos = overlaps.size();
   if (nos > 0) {
      for (unsigned int i=0; i<overlaps.size(); i++) {
         const atom_overlap_t &o(overlaps[i]);
         s += o.overlap_volume;
      }
      s /= static_cast<float>(nos);
      s *= 1000.0; // score per 1000 atoms
   }
   return s;
}


// static
coot::atom_overlaps_dots_container_t
coot::atom_overlaps_container_t::contacts_for_atom(int iat,
                                                   mmdb::Atom **atom_selection,
                                                   const std::map<int, std::vector<int> > &contact_map,
                                                   const std::map<int, std::vector<int> > &bonded_map,
                                                   const std::vector<double> &neighb_atom_radius,
                                                   int udd_h_bond_type_handle,
                                                   bool molecule_has_hydrogens,
                                                   double probe_radius,
                                                   double dot_density_in,
                                                   double clash_spike_length,
                                                   bool make_vdw_surface) {
   atom_overlaps_dots_container_t ao;

   mmdb::Atom *at = atom_selection[iat];

   clipper::Coord_orth pt_at_1 = co(at);

   int n_sphere_points       = static_cast<int>(dot_density_in * 900);
   int n_sphere_points_for_H = static_cast<int>(dot_density_in * 440);

   // Yikes! We don't need to calculate these for every atom! They can be pre-calculated
   // at the time of construction of this object - and made into a const variable ref that is passed to
   // this function.
   //
   std::vector<clipper::Coord_orth>   sphere_points = fibonacci_sphere(n_sphere_points);
   std::vector<clipper::Coord_orth> H_sphere_points = fibonacci_sphere(n_sphere_points_for_H); // less than above

   std::vector<clipper::Coord_orth> sphere_points_for_atom;
   if (std::string(at->element) == " H")
      sphere_points_for_atom = H_sphere_points;
   else
      sphere_points_for_atom = sphere_points;

   double r_1 = neighb_atom_radius[iat];

   for (unsigned int j=0; j<sphere_points_for_atom.size(); j++) {
      const clipper::Coord_orth &pt(sphere_points_for_atom[j]);
      clipper::Coord_orth pt_at_surface = r_1 * pt + pt_at_1;

      if (false) // yes these are points on a full sphere
         std::cout << "pt_at_surface "
                   << pt_at_surface.x() << " "
                   << pt_at_surface.y() << " "
                   << pt_at_surface.z() << "\n";

      std::map<int, std::vector<int> >::const_iterator it = bonded_map.find(iat);
      bool draw_it = ! is_inside_another_atom_to_which_its_bonded(iat, at,
                                                                  pt_at_surface,
                                                                  it->second,
                                                                  atom_selection,
                                                                  neighb_atom_radius);

      if (draw_it) { // it's on the vdw surface at least...

         double biggest_overlap = -1; // should be positive if we get a hit
         mmdb::Atom *atom_with_biggest_overlap = 0;
         double r_2_for_biggest_overlap = 0;

         // now check which atom this is clashing with (if any) and pick the
         // one with the biggest overlap
         //
         std::map<int, std::vector<int> >::const_iterator it_contact_map = contact_map.find(iat);
         const std::vector<int> &v = it_contact_map->second;

         for (unsigned int jj=0; jj<v.size(); jj++) {
            mmdb::Atom *neighb_atom = atom_selection[v[jj]];
            // no contacts for atoms in the same residue
            if (neighb_atom->residue == at->residue) continue;
            double r_2 = neighb_atom_radius[v[jj]];
            double r_2_sqrd = r_2 * r_2;
            double r_2_plus_prb_squard = r_2_sqrd + 2 * r_2 * probe_radius +
               4 * probe_radius * probe_radius;
            clipper::Coord_orth pt_na = co(neighb_atom);
            double d_sqrd = (pt_na - pt_at_surface).lengthsq();

            if (d_sqrd < r_2_plus_prb_squard) {
               // a contact dot on something
               double delta_d_sqrd = r_2_plus_prb_squard - d_sqrd;
               if (delta_d_sqrd > biggest_overlap) {
                  biggest_overlap = delta_d_sqrd;
                  atom_with_biggest_overlap = neighb_atom;
                  r_2_for_biggest_overlap = r_2;
               }
            }
         }

         if (atom_with_biggest_overlap) {
            double d_surface_pt_to_atom_sqrd = (co(atom_with_biggest_overlap) - pt_at_surface).lengthsq();
            double d_surface_pt_to_atom = sqrt(d_surface_pt_to_atom_sqrd);
            double overlap_delta = r_2_for_biggest_overlap - d_surface_pt_to_atom;
            // first is yes/no, second is H-is-on-ligand?
            // allow waters to H-bond (without being an H)
            // std::pair<bool, bool> might_be_h_bond_flag =
            // is_h_bond_H_and_acceptor(at, atom_with_biggest_overlap, udd_h_bond_type_handle);
            h_bond_info_t hbi(at, atom_with_biggest_overlap, udd_h_bond_type_handle);

            // std::pair<std::string, std::string> c_type_col =
            // overlap_delta_to_contact_type(overlap_delta, is_h_bond);

            // std::string c_type;
            // std::string col;

            // test_get_type(overlap_delta, is_h_bond, &c_type, &col);

            // test_get_type(overlap_delta, is_h_bond, &c_type, &col);

            // std::pair<std::string, std::string> ("wide-contact", "blue");
            // was
            // overlap_delta_to_contact_type(overlap_delta, is_h_bond);

            // const std::string &c_type = c_type_col.first;
            // const std::string &col    = c_type_col.second;

            std::pair<std::string, std::string> c_type_col = overlap_delta_to_contact_type(overlap_delta, hbi, molecule_has_hydrogens);
            const std::string &c_type = c_type_col.first;
            const std::string &col    = c_type_col.second;

            if (false) { //ddebug
               int h_bond_type = -999;
               if (at->GetUDData(udd_h_bond_type_handle, h_bond_type) == mmdb::UDDATA_Ok) {
                  std::cout << "at: " << atom_spec_t(at) << " " << atom_spec_t(atom_with_biggest_overlap)
                            << " delta " << std::setw(10) << overlap_delta << " hbi: " << hbi.format() << " "
                            << c_type << " " << col << " h_bond_type " << h_bond_type << std::endl;
               }
            }


            if (c_type != "clash") {

               // draw dot if these are atoms from different residues or this is not
               // a wide contact
               if (at->residue != atom_with_biggest_overlap->residue) {
                  atom_overlaps_dots_container_t::dot_t dot(overlap_delta, col, pt_at_surface);
                  ao.dots[c_type].push_back(dot);
               } else {
                  // this should not happen now that we have the test/continue above
               }
            } else {
               // clash
               clipper::Coord_orth vect_to_pt_1 = pt_at_1 - pt_at_surface;
               clipper::Coord_orth vect_to_pt_1_unit(vect_to_pt_1.unit());
               // these days, spikes project away from the atom, not inwards
               // clipper::Coord_orth pt_spike_inner =
               // pt_at_surface + clash_spike_length * vect_to_pt_1_unit;
               clipper::Coord_orth pt_spike_inner =
                  pt_at_surface -
                  0.34* sqrt(biggest_overlap) * clash_spike_length * vect_to_pt_1_unit;
               std::pair<clipper::Coord_orth, clipper::Coord_orth> p(pt_at_surface,
                                                                     pt_spike_inner);
               ao.clashes.positions.push_back(p);
            }

         } else {

            // no environment atom was close to this ligand atom, so just add
            // a surface point

            if (make_vdw_surface) {
               // we are outputing this point then...
               // std::cout << "vdw " << pt_at_surface.x()
               //           << " " << pt_at_surface.y()
               //           << " " << pt_at_surface.z() << "\n";
               atom_overlaps_dots_container_t::dot_t dot(0, "grey", pt_at_surface);
               ao.dots["vdw-surface"].push_back(dot);
            }
         }
      }
   }
   return ao;
}

bool
coot::atom_overlaps_container_t::clashable_alt_confs(mmdb::Atom *at_1, mmdb::Atom *at_2) const {

   bool r = true;

   std::string alt_conf_1 = at_1->altLoc;
   std::string alt_conf_2 = at_2->altLoc;

   if (alt_conf_1.empty()) {
      return true;
   } else {
      if (alt_conf_2.empty()) {
         return true;
      } else {
         return (alt_conf_1 == alt_conf_2);
      }
   }

   return r;
}



// for all-atom contacts
// static (because parallel)
bool
coot::atom_overlaps_container_t::is_inside_another_atom_to_which_its_bonded(int atom_idx,
                                                                            mmdb::Atom *at,
                                                                            const clipper::Coord_orth &pt_at_surface,
                                                                            const std::vector<int> &bonded_neighb_indices,
                                                                            mmdb::Atom **atom_selection,
                                                                            const std::vector<double> &neighb_atom_radius) {

   bool r = false;
   // double r_1 = get_vdw_radius_neighb_atom(atom_idx);
   double r_1 = neighb_atom_radius[atom_idx];

   for (unsigned int i=0; i<bonded_neighb_indices.size(); i++) {
      mmdb::Atom *clash_neighb = atom_selection[bonded_neighb_indices[i]];
      clipper::Coord_orth pt_clash_neigh = co(clash_neighb);
      // double r_2 = get_vdw_radius_neighb_atom(bonded_neighb_indices[i]);
      double r_2 = neighb_atom_radius[bonded_neighb_indices[i]];
      double r_2_sqrd = r_2 * r_2;
      double d_sqrd = (pt_at_surface - pt_clash_neigh).lengthsq();
//       std::cout << "debug this_atom " << atom_spec_t(at) << " has r_1 " << r_1
//                 << " clash atom " << atom_spec_t(clash_neighb) << " has r_2 " << r_2 << std::endl;
      if (d_sqrd < r_2_sqrd) {
         r = true;
         break;
      }
   }

   return r;
}

bool
coot::atom_overlaps_container_t::is_inside_an_env_atom_to_which_its_bonded(int idx,
                                                                           const std::vector<int> &bonded_neighb_indices,
                                                                           mmdb::Atom **env_residue_atoms,
                                                                           const clipper::Coord_orth &pt_at_surface) {

   bool r = false;
   double r_1 = get_vdw_radius_neighb_atom(idx);
   for (unsigned int i=0; i<bonded_neighb_indices.size(); i++) {
      mmdb::Atom *env_atom = env_residue_atoms[bonded_neighb_indices[i]];
      // std::cout << "testing env_atom: " << atom_spec_t(env_atom) << std::endl;
      clipper::Coord_orth pt_env_atom = co(env_atom);
      double r_2 = 1.6;
      std::string ele(env_atom->element);
      if (ele == " H")
         r_2 = 0.97;
      double r_2_sqrd = r_2 * r_2;
      double d_sqrd = (pt_at_surface - pt_env_atom).lengthsq();
      if (d_sqrd < r_2_sqrd) {
         // std::cout << "   " << pt_at_surface.format() << " is inside atom " << atom_spec_t(env_atom) << std::endl;
         r = true;
         break;
      }
   }

   return r;
}


// static
void
coot::atom_overlaps_container_t::test_get_type(double delta, bool is_h_bond,
                                               std::string *c_type_p, std::string *col_p) {

   *c_type_p = std::string("wide-contact");
   *col_p    = std::string("blue");
}



// return H-bond, or wide-contact or close-contact or small-overlap or big-overlap or clash
// We add the mapping to colour here so that it's clearer the relationship betwen the overlap_delta
// and the colour.  One contact type can have 2 or more colours.
//
// static
std::pair<std::string, std::string>
coot::atom_overlaps_container_t::overlap_delta_to_contact_type(double delta, bool is_h_bond) {

// from the Word 1999 paper:
// pale green for H-bonds
// green (narrow gaps) or yellow (slight overlaps, < 0.2) for good contacts
// blue for wider gaps > 0.25
// orange and red for unfavourable (0.25 to 0.4)
// hot pink for >= 0.4

   std::string type = "wide-contact";
   std::string colour = "sky";

   if (is_h_bond) {
      if (delta >= -0.15) { // not 0, so that we turn small overlaps to H-bond dots
         delta -= 0.8;
         if (delta > 0.4) {
            type = "clash";
            colour = "hotpink";
         } else {
            type = "H-bond";
            colour = "darkpurple";
         }
      }
   } else {

      if (delta > -0.3) {
         type = "close-contact";
         colour = "royalblue";
      }

      if (delta > -0.2) {
         type = "close-contact";
         colour = "sea";
      }

      if (delta > -0.1) {
         type = "small-overlap";
         colour = "green";
      }

      if (delta > 0.15) {
         type = "small-overlap";
         colour = "yellow";
      }

      if (delta > 0.25) {
         type = "small-overlap";
         colour = "orange";
      }

      if (delta > 0.35) {
         type = "small-overlap";
         colour = "orangered";
      }

      if (delta > 0.40) {        // was 0.30, not enough red
         type = "big-overlap";
         colour = "red";
      }

      if (delta > 0.5) {        // Word: 0.4
         type = "clash";
         colour = "hotpink";
      }
   }
   return std::pair<std::string, std::string> (type, colour);

}

// a h_bond_info_t contains information about donor-acceptor as well as H-acceptor
//
std::pair<std::string, std::string>
coot::atom_overlaps_container_t::overlap_delta_to_contact_type(double delta,
                                                               const h_bond_info_t &hbi,
                                                               bool molecule_has_hydrogens_flag) {

// from the Word 1999 paper:
// pale green for H-bonds
// green (narrow gaps) or yellow (slight overlaps, < 0.2) for good contacts
// blue for wider gaps > 0.25
// orange and red for unfavourable (0.25 to 0.4)
// hot pink for >= 0.4

   // If the contact is a hydrogen-bond, then
   //   it's either clash or green-pillow
   // else
   //   rainbow colours
   //
   // Note: I am not sure that I like this for the case of H-bond donor and acceptor
   // with no Hydrogen atoms in the model.

   std::string type = "wide-contact";
   std::string colour = "sky";

   bool done = false;

   // return std::make_pair(std::string("H-bond"), std::string("sea"));


   if (hbi.is_h_bond_H_and_acceptor) {
      done = true;
      delta -= 2.0; // Hydrogen atom can get to within 1.5A of the acceptor atom
      if (delta >= -0.15) { // not 0, so that we turn small overlaps to H-bond dots
         if (delta > 0.4) {
            type = "clash";
            colour = "hotpink";
         } else {
            type = "H-bond";
            colour = "darkpurple";
         }
      } else {
        type = "H-bond";
        colour = "darkpurple";
      }
   } else {

      if (! molecule_has_hydrogens_flag) {
         if (hbi.is_h_bond_donor_and_acceptor) {
            done = true;
            delta -= 0.4; // where does this number come from?
            if (delta > 0.4) {
               type = "clash";
               colour = "hotpink";
            } else {
               type = "H-bond";
               colour = "greentint";
            }
         }
      }

   }

   if (! done) {

     // "not hydrogen bond" colouring

      if (delta > -0.3) {
         type = "close-contact";
         colour = "royalblue";
      }

      if (delta > -0.2) {
         type = "close-contact";
         colour = "sea";
      }

      if (delta > -0.1) {
         type = "small-overlap";
         colour = "green";
      }

      if (delta > 0.10) {
         type = "small-overlap";
         colour = "yellow";
      }

      if (delta > 0.18) {
         type = "small-overlap";
         colour = "orange";
      }

      if (delta > 0.25) {
         type = "small-overlap";
         colour = "orangered";
      }

      if (delta > 0.3) {
         type = "big-overlap";
         colour = "red";
      }

      if (delta > 0.4) {         // Word: 0.4
         type = "clash";
         colour = "hotpink";
      }
   }

   if (false)
      std::cout << "overlap_delta_to_contact_type() " << delta
                << " H-bond-H-and-acceptor-status " << hbi.is_h_bond_H_and_acceptor
                << " returns " << type << " " << colour << std::endl;
   return std::pair<std::string, std::string> (type, colour);

}


void
coot::atom_overlaps_container_t::add_residue_neighbour_index_to_neighbour_atoms() {

   udd_residue_index_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "neighb-residue-index");
   for (unsigned int i=0; i<neighbours.size(); i++) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      neighbours[i]->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         at->PutUDData(udd_residue_index_handle, int(i));
      }
   }
}



// fill std::vector<double> env_residue_radii;
//
// in the case of all-atom env_residues is all atoms
void
coot::atom_overlaps_container_t::setup_env_residue_atoms_radii(int i_sel_hnd_env_atoms) {

   if (! neighb_atom_radius.empty()) return;

   if (!have_dictionary) {
      std::cout << "setup_env_residue_atoms_radii() no dictionary " << std::endl;
   }
   double r = 1.5;
   mmdb::Atom **env_residue_atoms = 0;
   int n_env_residue_atoms;
   mol->GetSelIndex(i_sel_hnd_env_atoms, env_residue_atoms, n_env_residue_atoms);
   neighb_atom_radius.resize(n_env_residue_atoms);

   if (false) {
      for (int i=0; i<n_env_residue_atoms; i++) {
         mmdb::Atom *at = env_residue_atoms[i];
         std::cout << "atom at " << i << " " << coot::atom_spec_t(at) << std::endl;
      }
   }

   for (int i=0; i<n_env_residue_atoms; i++) {
      mmdb::Atom *at = env_residue_atoms[i];
      // std::cout << "loop that breaks " << i << " " << at << " " << at->GetAtomName() << std::endl;
      int residue_index = -1;
      int ierr = at->GetUDData(udd_residue_index_handle, residue_index);

      if (ierr == mmdb::UDDATA_Ok) {

         // const dictionary_residue_restraints_t &rest = neighb_dictionaries[residue_index];
         try {

            // in ALL_ATOM mode the lookup uses the residue name, not the index.
            // there should be two functions because this is confusing.
            if (false)
               std::cout << ":::::::::::::::::: in setup_env_residue_atoms_radii() the dictionary_map is of size "
                         << dictionary_map.size() << " and residue_index is " << residue_index << std::endl;

            mmdb::Residue *res = at->residue;
            const dictionary_residue_restraints_t &rest = get_dictionary(res, residue_index);
            // is rest sane?
            // std::cout << "debug:: rest has " << rest.atom_info.size() << " atoms "<< std::endl;
            // std::cout << "debug:: rest has " << rest.bond_restraint.size() << " bond restraints "<< std::endl;
            std::string atom_name = at->GetAtomName();
            std::string te = rest.type_energy(atom_name);
            if (! te.empty()) {
               std::map<std::string, double>::const_iterator it_type = type_to_vdw_radius_map.find(te);
               if (it_type == type_to_vdw_radius_map.end()) {
                  if (geom_p)
                     r = type_energy_to_radius(te);
                  type_to_vdw_radius_map[te] = r;
                  // std::cout << "debug type_to_vdw_radius_map " << std::setw(4) << te << "   " << r << std::endl;
               } else {
                  r = it_type->second;
               }
               neighb_atom_radius[i] = r;
            }
         }
         catch (const std::out_of_range &ex) {
            std::cout << "OOpps " << ex.what() << std::endl;
         }
      } else {
         std::cout << "ERROR:: failed to get UDData for residue index" << std::endl;
      }
   }
}

double
coot::atom_overlaps_container_t::type_energy_to_radius(const std::string &te) const {

   double r = 1.2; // 1.17 probe values (from the paper):
   double r_arom  = 1.05; // 1.0 (paper)
   double r_polar = 1.05; // 1.0 (paper) 1.05 (measured)
   // C-H H   1.17
   // arom H  1.0
   // polar H 1.0
   // didn't find te in types map. so look it up from the dictionary and add to the types map
   if (te[0] == 'H') {

      // HCH1 HCH2 HCH3 HCR6 HCR5 are C-H (default)
      //
      if (te == "HNH1") r = r_polar;
      if (te == "HNT3") r = r_polar;
      if (te == "HOH1") r = r_polar; // HOH H
      if (te == "HNC1") r = r_polar; // e.g. ARG
      if (te == "HNC2") r = r_polar; // ARG
      if (te == "HNH1") r = r_polar; // mainchain H
      if (te == "H")    r = r_polar; // unset
      if (te == "HNH2") r = r_polar; // ND2 of ASN
      if (te == "HNR5") r = r_arom;  // arom
   } else {
      r = geom_p->get_energy_lib_atom(te).vdw_radius;
   }
   return r;
}

// We need here a flag for main-chain -> mainchain interactions also.
// Currently only consider side-chain -> side-chain and side-chain -> main-chain.
//
// Needs clear thinking.
//
// If atoms are bonded/1-3/ring related (or (note to self) possibly also 1-4?) then the atoms
// cannot clash and only these atoms can be neighbouring atoms that are tested to see if
// surface points are within the sphere of other atoms.
//
// so, either clashable or bonded or ignored
//
coot::atom_overlaps_container_t::atom_interaction_type
coot::atom_overlaps_container_t::bonded_angle_or_ring_related(mmdb::Manager *mol,
                                                              mmdb::Atom *at_1,
                                                              mmdb::Atom *at_2,
                                                              bool exclude_mainchain_also,
                                                              std::map<std::string, std::vector<std::pair<std::string, std::string> > > *bonded_neighbours,
                                                              std::map<std::string, std::vector<std::vector<std::string> > > *ring_list_map) {

   // At the N terminus, the N has attached Hydrogen atoms H1, H2, H3.  These should not clash
   // with the N Nitrogen atom.

   // atom_spec_t debug_spec("A", 383, "", " OG ", "");
   atom_interaction_type ait = CLASHABLE;
   mmdb::Residue *res_1 = at_1->GetResidue();
   mmdb::Residue *res_2 = at_2->GetResidue();

   if (res_1 != res_2) {
      if (are_bonded_residues(res_1, res_2)) { // just checks seqNums, should check residue index
         if (is_main_chain_p(at_1)) {
            if (is_main_chain_p(at_2)) {
               ait = BONDED; // :-) everything mainchain is bonded to each other (for dot masking)
            } else {
               // this is a bit hacky - don't allow clashes to CDs in proline from mainchain atoms
               // (in other residues)
               std::string res_name_2 = res_2->GetResName();
               if (res_name_2 == "PRO") {
                  std::string at_name_1 = at_1->GetAtomName();
                  std::string at_name_2 = at_2->GetAtomName();
                  if (at_name_2 == " CD ") {
                     ait = BONDED;
                  } else {
                     ait = CLASHABLE;
                  }
               } else {
                  ait = CLASHABLE;
               }
            }
         } else {

            // ---- at_1 is not main-chain
            if (is_main_chain_p(at_2)) {
               // perhaps the PRO CD - C atoms were the other way around (to above)
               std::string at_name_1 = at_1->GetAtomName();
               if (at_name_1 == " CD ") {
                  std::string res_name_2 = res_2->GetResName();
                  if (res_name_2 == "PRO") {
                     ait = BONDED;
                  } else {
                     ait = CLASHABLE;
                  }
               } else {
                  ait = CLASHABLE;
               }
            } else {
               ait = BONDED;
            }
         }
      } else {
         std::string res_name_1 = res_1->GetResName();
         std::string res_name_2 = res_2->GetResName();

         if (res_name_1 == res_name_2) {
            if (res_name_1 == "HOH") {
               if (ignore_water_contacts_flag) {
                  ait = IGNORED;
               } else {
                  ait = CLASHABLE;
               }
            }
         } else {
            ait = CLASHABLE;
         }

         if (false) { // debug
            std::cout << "------------------ in bonded_angle_or_ring_related() with res names " << res_name_1 << " " << res_name_2
                      << " and ignore_water_contacts_flag " << ignore_water_contacts_flag << " ait: " << ait << std::endl;
         }

      }
   } else {

      // ------------------- same residue ------------------------
      //
      std::string res_name = res_1->GetResName();
      std::vector<std::pair<std::string, std::string> > bps;
      std::string atom_name_1 = at_1->name;
      std::string atom_name_2 = at_2->name;
      std::map<std::string, std::vector<std::pair<std::string, std::string> > >::const_iterator it;
      it = bonded_neighbours->find(res_name);
      if (it == bonded_neighbours->end()) {
         bps = geom_p->get_bonded_and_1_3_angles(res_name, protein_geometry::IMOL_ENC_ANY);
         // bonded_neighbours[res_name] = bps; // when arg is reference
         // bonded_neighbours->insert(res_name, bps);
         (*bonded_neighbours)[res_name] = bps;
      } else {
         bps = it->second;
      }
      for (unsigned int ipr=0; ipr<bps.size(); ipr++) {

         if (atom_name_1 == bps[ipr].first) {
            if (atom_name_2 == bps[ipr].second) {
               ait = BONDED;
               break;
            }
         }
         if (atom_name_2 == bps[ipr].first) {
            if (atom_name_1 == bps[ipr].second) {
               ait = BONDED;
               break;
            }
         }
      }
      if (ait == CLASHABLE) { // i.e. so far, unset by this block
         // std::cout << "checking for in-same-ring for " << coot::atom_spec_t(at_1) << " " << coot::atom_spec_t(at_2) << std::endl;
         bool ringed = in_same_ring(at_1, at_2, *ring_list_map); // update ring_list_map
         if (ringed) ait = BONDED;
      }
      if (ait == CLASHABLE) {
         // N-terminus Hs?

         // Don't clash on H-H interaction or H-CA either.
         //
         if (res_1->isNTerminus()) {

            std::string at_name_1 = at_1->GetAtomName();
            std::string at_name_2 = at_2->GetAtomName();
            if (at_name_1 == " H1 " || at_name_1 == " H2 " || at_name_1 == " H3 ") {
               if (at_name_2 == " H1 " || at_name_2 == " H2 " || at_name_2 == " H3 " ||
                   at_name_2 == " CA " || at_name_2 == " N  ") {
                  ait = BONDED;
               }
            } else {
               if (at_name_1 == " CA " || at_name_1 == " N  ") {
                  if (at_name_2 == " H1 " || at_name_2 == " H2 " || at_name_2 == " H3 ")
                     ait = BONDED;
               }
            }
         }
      }
   }

   if (ait == CLASHABLE) {
      // maybe it was a link
      if (is_linked(at_1, at_2) || is_ss_bonded_or_CYS_CYS_SGs(at_1, at_2)) {
         ait = BONDED;
      } else {
         std::vector<std::pair<std::string, std::string> > bonds_for_at_1;
         std::vector<std::pair<std::string, std::string> > bonds_for_at_2;
         std::string res_name_1(at_1->GetResName());
         std::string res_name_2(at_2->GetResName());
         std::map<std::string, std::vector<std::pair<std::string, std::string> > >::const_iterator it_1;
         std::map<std::string, std::vector<std::pair<std::string, std::string> > >::const_iterator it_2;
         it_1 = bonded_neighbours->find(res_name_1);
         it_2 = bonded_neighbours->find(res_name_2);
         if (it_1 != bonded_neighbours->end()) bonds_for_at_1 = it_1->second;
         if (it_2 != bonded_neighbours->end()) bonds_for_at_2 = it_2->second;
         if (is_angle_related_via_link(at_1, at_2, bonds_for_at_1, bonds_for_at_2))
            ait = BONDED;
      }
   }
   return ait;
}

bool
coot::atom_overlaps_container_t::is_angle_related_via_link(mmdb::Atom *at_1,
                                                           mmdb::Atom *at_2,
                                     const std::vector<std::pair<std::string, std::string> > &bonds_for_at_1,
                                     const std::vector<std::pair<std::string, std::string> > &bonds_for_at_2) const {
   bool status = false;
   if (! at_1) return false;
   if (! at_2) return false;

   mmdb::Model *model_p_1 = at_1->GetModel();
   mmdb::Model *model_p_2 = at_2->GetModel();

   if (model_p_2 != model_p_1) return false;

   if (false) { // are you sure that bonds_for_at_1 and bonds_for_at_2 were created with the correct arguments?
      for (auto const &b : bonds_for_at_1)
         std::cout << std::string(at_1->GetAtomName()) << " "
                   << std::string(at_2->GetAtomName()) << " "
                   << " bonds for at_1 " << b.first << " " << b.second << std::endl;
      for (auto const &b : bonds_for_at_2)
         std::cout << std::string(at_1->GetAtomName()) << " "
                   << std::string(at_2->GetAtomName()) << " "
                   << " bonds for at_2 " << b.first << " " << b.second << std::endl;
   }

   if (model_p_1) {
      int n_links = model_p_1->GetNumberOfLinks();
      if (n_links > 0) {
         for (int i_link=1; i_link<=n_links; i_link++) {
            mmdb::Link *link = model_p_1->GetLink(i_link);
            if (link) {
               std::pair<atom_spec_t, atom_spec_t> linked_atoms = link_atoms(link, model_p_1);
               atom_spec_t spec_1(at_1);
               atom_spec_t spec_2(at_2);

               if (spec_1 == linked_atoms.first) {
                  // is spec_1/at_1 linked to an atom that is bonded to at_2?
                  std::string linked_atom_2_name = linked_atoms.second.atom_name;
                  for (unsigned int i=0; i<bonds_for_at_2.size(); i++) {
                     const std::string &bond_atom_1_name = bonds_for_at_2[i].first;
                     const std::string &bond_atom_2_name = bonds_for_at_2[i].second;
                     if (bond_atom_1_name == linked_atom_2_name)
                        if (bond_atom_2_name == spec_2.atom_name)
                           status = true;
                     if (bond_atom_2_name == linked_atom_2_name)
                        if (bond_atom_1_name == spec_2.atom_name)
                           status = true;
                     if (status) break;
                  }
               }

               if (spec_1 == linked_atoms.second) {
                  // is spec_1/at_1 linked to an atom that is bonded to at_2?
                  std::string linked_atom_1_name = linked_atoms.first.atom_name;
                  for (unsigned int i=0; i<bonds_for_at_2.size(); i++) {
                     const std::string &bond_atom_1_name = bonds_for_at_2[i].first;
                     const std::string &bond_atom_2_name = bonds_for_at_2[i].second;
                     if (bond_atom_1_name == linked_atom_1_name)
                        if (bond_atom_2_name == spec_2.atom_name)
                           status = true;
                     if (bond_atom_2_name == linked_atom_1_name)
                        if (bond_atom_1_name == spec_2.atom_name)
                           status = true;
                     if (status) break;
                  }
               }

               // and now the other way round.
               if (spec_2 == linked_atoms.first) {
                  // is spec_2/at_2 linked to an atom that is bonded to at_1?
                  std::string linked_atom_1_name = linked_atoms.first.atom_name;
                  for (unsigned int i=0; i<bonds_for_at_1.size(); i++) {
                     const std::string &bond_atom_1_name = bonds_for_at_1[i].first;
                     const std::string &bond_atom_2_name = bonds_for_at_1[i].second;
                     if (bond_atom_1_name == linked_atom_1_name) {
                        if (bond_atom_2_name == spec_1.atom_name)
                           status = true;
                     }
                     if (bond_atom_2_name == linked_atom_1_name) {
                        if (bond_atom_1_name == spec_1.atom_name)
                           status = true;
                     }
                     if (status) break;
                  }
               }

               if (spec_2 == linked_atoms.second) {
                  // is spec_2/at_2 linked to an atom that is bonded to at_1?
                  std::string linked_atom_1_name = linked_atoms.first.atom_name;
                  for (unsigned int i=0; i<bonds_for_at_1.size(); i++) {
                     const std::string &bond_atom_1_name = bonds_for_at_1[i].first;
                     const std::string &bond_atom_2_name = bonds_for_at_1[i].second;
                     if (bond_atom_1_name == linked_atom_1_name) {
                        if (bond_atom_2_name == spec_1.atom_name)
                           status = true;
                     }
                     if (bond_atom_2_name == linked_atom_1_name) {
                        if (bond_atom_1_name == spec_1.atom_name)
                           status = true;
                     }
                     if (status) break;
                  }
               }
            }
            if (status) break;
         }
      }
   }
   return status;
}

bool
coot::atom_overlaps_container_t::is_linked(mmdb::Atom *at_1,
                                           mmdb::Atom *at_2) const {

   bool status = false;
   if (! at_1) return false;
   if (! at_2) return false;

   mmdb::Model *model_p_1 = at_1->GetModel();
   mmdb::Model *model_p_2 = at_2->GetModel();

   if (model_p_2 != model_p_1) return false;

   if (model_p_1) {
      int n_links = model_p_1->GetNumberOfLinks();
      if (n_links > 0) {
         for (int i_link=1; i_link<=n_links; i_link++) {
            mmdb::Link *link = model_p_1->GetLink(i_link);
            if (link) {
               std::pair<atom_spec_t, atom_spec_t> atoms = link_atoms(link, model_p_1);
               atom_spec_t spec_1(at_1);
               atom_spec_t spec_2(at_2);
               if (spec_1 == atoms.first) {
                  if (spec_2 == atoms.second) {
                     status = true;
                     break;
                  }
               }
               if (spec_2 == atoms.first) {
                  if (spec_1 == atoms.second) {
                     status = true;
                     break;
                  }
               }
            }
         }
      }
   }
   return status;
}

bool
coot::atom_overlaps_container_t::is_ss_bonded_or_CYS_CYS_SGs(mmdb::Atom *at_1,
                                                             mmdb::Atom *at_2) const {

   // There is no mmdb class for SSBOND! - must fix.

   bool status = false;
   std::string res_name_1 = at_1->residue->GetResName();
   if (res_name_1 == "CYS") {
      std::string res_name_2 = at_2->residue->GetResName();
      if (res_name_2 == "CYS") {
         std::string atom_name_1 = at_1->GetAtomName();
         if (atom_name_1 == " SG ") {
            std::string atom_name_2 = at_2->GetAtomName();
            if (atom_name_2 == " SG ") {
               status = true;
            }
         }
      }
   }

   return status;
}

bool
coot::atom_overlaps_container_t::is_ss_bonded(mmdb::Residue *residue_p) const {

   // There is no mmdb class for SSBOND! - must fix.

   bool status = false;
   if (residue_p) {
      std::string res_name = residue_p->GetResName();
      if (res_name == "CYS") {
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         // check SS bonds here
         status = true;
      }
   }

   return status;
}



// this modifies ring_list_map, so it is not const.
bool
coot::atom_overlaps_container_t::in_same_ring(mmdb::Atom *at_1, mmdb::Atom *at_2,
                                              std::map<std::string, std::vector<std::vector<std::string> > > &ring_list_map) const {

   bool same = false;
   mmdb::Residue *res_1 = at_1->GetResidue();
   mmdb::Residue *res_2 = at_2->GetResidue();
   if (res_1 == res_2) {
      std::string at_name_1 = at_1->GetAtomName();
      std::string at_name_2 = at_2->GetAtomName();
      try {
         const dictionary_residue_restraints_t &dict = get_dictionary(res_1, 0);
         std::string res_name = at_1->GetResName();

         // old "calculate ring list each time" way
         // same = dict.in_same_ring(at_name_1, at_name_2);

         std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator it;
         it = ring_list_map.find(res_name);
         if (it != ring_list_map.end()) {
            // same = in_same_ring(at_name_1, at_name_2, it->second);
            same = dict.in_same_ring(at_name_1, at_name_2, it->second);
         } else {
            // on reflection, 5-membered rings are probably not needed because in those cases all atoms are
            // bonded or 1-3 angles
            std::vector<std::vector<std::string> > ring_list;
            if (res_name == "HIS") {
               ring_list = his_ring_list();
            } else {
               if (res_name == "PHE" || res_name == "TYR") {
                  ring_list = phe_ring_list();
               } else {
                  if (res_name == "TRP") {
                     ring_list = trp_ring_list();
                  } else {
                     if (res_name == "PRO") {
                        ring_list = pro_ring_list();
                     } else {
                        ring_list = dict.get_ligand_ring_list();
                     }
                  }
               }
            }
            ring_list_map[res_name] = ring_list;
            // same = in_same_ring(at_name_1, at_name_2, ring_list);
            same = dict.in_same_ring(at_name_1, at_name_2, ring_list);
         }
      }
      catch (const std::out_of_range &ex) {
         std::cout << "OOpps " << ex.what() << std::endl;
         logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__),
                    "Ooops", ex.what());
      }
   }
   return same;
}

// moved to dictionary_residue_restraints_t
//
// bool
// coot::atom_overlaps_container_t::in_same_ring(const std::string &atom_name_1,
//                                               const std::string &atom_name_2,
//                                               const std::vector<std::vector<std::string> > &ring_list) const {

//    bool match = false;
//    for (unsigned int i=0; i<ring_list.size(); i++) {
//       unsigned int n_match = 0;
//       for (unsigned int j=0; j<ring_list[i].size(); j++) {
//          if (ring_list[i][j] == atom_name_1)
//             n_match++;
//          if (ring_list[i][j] == atom_name_2)
//             n_match++;
//       }
//       if (n_match == 2) {
//          match = true;
//          break;
//       }
//    }
//    return match;
// }

std::vector<std::vector<std::string> >
coot::atom_overlaps_container_t::phe_ring_list() const {
   std::vector<std::vector<std::string> > v;
   std::vector<std::string> vi(6);
   vi[0] = " CG ";
   vi[1] = " CD1";
   vi[2] = " CD2";
   vi[3] = " CE1";
   vi[4] = " CE2";
   vi[5] = " CZ ";
   v.push_back(vi);
   return v;
}


std::vector<std::vector<std::string> >
coot::atom_overlaps_container_t::his_ring_list() const {
   std::vector<std::vector<std::string> > v;
   std::vector<std::string> vi(5);
   vi[0] = " CG ";
   vi[1] = " ND1";
   vi[2] = " CD2";
   vi[3] = " NE2";
   vi[4] = " CE1";
   v.push_back(vi);
   return v;
}

std::vector<std::vector<std::string> >
coot::atom_overlaps_container_t::trp_ring_list() const {

   std::vector<std::vector<std::string> > v;
   std::vector<std::string> vi(5);
   std::vector<std::string> vi2(6);
   vi[0] = " CG ";
   vi[1] = " CD1";
   vi[2] = " NE1";
   vi[3] = " CE2";
   vi[4] = " CD2";
   vi2[0] = " CE2";
   vi2[1] = " CD2";
   vi2[2] = " CE3";
   vi2[3] = " CZ3";
   vi2[4] = " CH2";
   vi2[5] = " CZ2";
   v.push_back(vi);
   v.push_back(vi2);
   return v;
}

std::vector<std::vector<std::string> >
coot::atom_overlaps_container_t::pro_ring_list() const {

   std::vector<std::vector<std::string> > v;
   std::vector<std::string> vi(5);
   vi[0] = " CA ";
   vi[1] = " CB";
   vi[2] = " CG";
   vi[3] = " CD";
   vi[4] = " N";
   v.push_back(vi);
   return v;
}


bool
coot::atom_overlaps_container_t::are_bonded_residues(mmdb::Residue *res_1, mmdb::Residue *res_2) const {

   bool r = false;

   if (res_1) {
      if (res_2) {
         if (res_1->GetChain() == res_2->GetChain()) {

            // simple

            if (abs(res_1->GetSeqNum() - res_2->GetSeqNum()) < 2) {
               std::string res_name_1 = res_1->GetResName();
               std::string res_name_2 = res_2->GetResName();
               if (res_name_1 != "HOH")
                  if (res_name_2 != "HOH")
                     r = true;
            }
         }

         if (! r) {
            // perhaps they are next to each other by serial number but not sequence number
            int ser_1 = res_1->index;
            int ser_2 = res_2->index;
            if (ser_2 > ser_1) {
               if ((ser_2 - ser_1) == 1) {
                  if (! res_1->isCTerminus()) {
                     // check that the linking atoms are close enough
                     // residue type-specific check. not for here
                     r = true;
                  }
               }
            }
            if (ser_1 > ser_2) { // backwards indexing
               if ((ser_1 - ser_2) == 1) {
                  if (! res_2->isCTerminus()) {
                     r = true;
                  }
               }
            }
         }
      }
   }
   return r;
}

// this should be a vector of derived symmetry_atom class really.
// the second atom needs symmetry info (if we are going to use it).
//
std::vector<coot::atom_overlap_t>
coot::atom_overlaps_container_t::symmetry_contacts(float d) {

   std::vector<atom_overlap_t> v;

   int n_symm = mol->GetNumberOfSymOps();
   int shift_lim = 2;
   mmdb::realtype min_contact_dist = d;
   mmdb::mat44 my_matt;

   residue_spec_t spec(res_central);

   int selHnd_1 = mol->NewSelection(); // d
   int selHnd_2 = mol->NewSelection(); // d
   mmdb::Atom **atom_selection_res;
   mmdb::Atom **atom_selection_all;
   int n_selection_atoms;
   int n_all_atoms;
   spec.select_atoms(mol, selHnd_1, mmdb::SKEY_NEW);
   mol->SelectAtoms(selHnd_2, 1, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*"," * ","*","*" );
   mol->GetSelIndex(selHnd_1, atom_selection_res, n_selection_atoms);
   mol->GetSelIndex(selHnd_2, atom_selection_all, n_all_atoms);

   for (int x_shift= -shift_lim; x_shift <= shift_lim; x_shift++) {
      for (int y_shift= -shift_lim; y_shift <= shift_lim; y_shift++) {
         for (int z_shift= -shift_lim; z_shift <= shift_lim; z_shift++) {
            for (int i_symm=0; i_symm < n_symm; i_symm++) {
               if (! (x_shift == 0 && y_shift == 0 && z_shift == 0 && i_symm == 0)) {
                  int i_status = mol->GetTMatrix(my_matt, i_symm, x_shift, y_shift, z_shift);

                  if (i_status == 0) { // Happy

                     mmdb::Contact *contact = NULL;
                     int ncontacts = 0;
                     long i_contact_group = 1;

                     mol->SeekContacts(atom_selection_res, n_selection_atoms,
                                       atom_selection_all, n_all_atoms,
                                       0.001, min_contact_dist,
                                       0,
                                       contact, ncontacts,
                                       0, &my_matt, i_contact_group);
                     if (ncontacts) {

                        // symm_trans_t st(i_symm, x_shift, y_shift, z_shift);
                        for (int i=0; i< ncontacts; i++) {
                           mmdb::Atom *at_1 = atom_selection_res[contact[i].id1];
                           mmdb::Atom *at_2 = atom_selection_all[contact[i].id2];
                           atom_overlap_t aop(at_1, at_2);
                           v.push_back(aop);
                        }
                     }
                  }
               }
            }
         }
      }
   }

   mol->DeleteSelection(selHnd_1);
   mol->DeleteSelection(selHnd_2);

   return v;
}
