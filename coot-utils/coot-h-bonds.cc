/* coot-utils/coot-h-bonds.cc
 * 
 * Copyright 2010 by The University of Oxford
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

#include <chrono>

#include <sstream>
#include "coot-h-bonds.hh"


std::ostream &
coot::operator<<(std::ostream &s, h_bond hb) {

   bool H_is_water = false;
   std::string res_name = hb.hb_hydrogen->GetResName();
   if (res_name == "HOH")
      H_is_water = true;

   if (hb.hb_hydrogen) {
      s << "H: "<< coot::atom_spec_t(hb.hb_hydrogen) << " ";
      if (H_is_water)
         s << " (HOH) ";
   } else {
      s << "H: NULL ";
   }
   if (! H_is_water)
      s << "donor: " << coot::atom_spec_t(hb.donor);
   s << " acceptor: " << coot::atom_spec_t(hb.acceptor);
   if (hb.donor_neigh)
      s << " donor_neigh: " << coot::atom_spec_t(hb.donor_neigh);
   else
      s << " donor_neigh: NULL ";

   if (hb.acceptor_neigh)
      s << " acceptor_neigh: " << coot::atom_spec_t(hb.acceptor_neigh);
   else
      s << " acceptor_neigh: NULL [problem!?]";

   s << " dist: " << hb.dist
     << " ligand-atom-is-donor?: " << hb.ligand_atom_is_donor;
   return s;
}

// typically the atom selection selHnd_1 is for the ligand and
// selHnd_2 is for everything (else).
// 
std::vector<coot::h_bond>
coot::h_bonds::get(int selHnd_1, int selHnd_2, mmdb::Manager *mol, const coot::protein_geometry &geom, int imol) {

   bool debug = 0;
   
   mmdb::realtype min_dist = 2.4; // H-bonds are longer than this
   mmdb::realtype max_dist = 3.9; // H-bonds are shorter than this
   
   std::vector<coot::h_bond> v;

   int hb_type_udd_handle = mark_donors_and_acceptors(selHnd_1, selHnd_2, mol, geom, imol); // using UDD data

   // What is the nearest neighbour of the atoms in mol?
   // 
   std::map<mmdb::Atom *, std::vector<std::pair<mmdb::Atom *, float> > > neighbour_map =
      make_neighbour_map(selHnd_1, selHnd_2, mol);

   mmdb::PPAtom sel_1_atoms = 0;
   mmdb::PPAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mmdb::mat44 my_matt;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
         my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   
   mol->GetSelIndex   (selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex   (selHnd_2, sel_2_atoms, n_sel_2_atoms);

   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;

   auto tp_0 = std::chrono::high_resolution_clock::now();

   mol->SeekContacts(sel_1_atoms, n_sel_1_atoms,
                     sel_2_atoms, n_sel_2_atoms,
                     min_dist, max_dist,
                     0, // seqDist 0 -> also in same res.
                     pscontact, n_contacts,
                     0, &my_matt, i_contact_group);

   auto tp_1 = std::chrono::high_resolution_clock::now();

   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();

   // small
   // std::cout << "------------------- timing: " << d10 <<  " milliseconds for SeekContacts "
   //            << std::endl;

   if (debug)
      std::cout << "h_bonds:: get(): found n_contacts between atom selections: " << n_contacts << std::endl;

   if (n_contacts > 0) {
      if (pscontact) {
         // we want to stop water-water Hbonds being added twice.
         // so store the hashed key here.
         std::set<std::string> existing_h_bond_keys;
         for (int i_contact=0; i_contact<n_contacts; i_contact++) {
            mmdb::Atom *at_1 = sel_1_atoms[pscontact[i_contact].id1];
            mmdb::Atom *at_2 = sel_2_atoms[pscontact[i_contact].id2];

            // are they donor and acceptor?
            //
            int hb_type_1 = coot::HB_UNASSIGNED;
            int hb_type_2 = coot::HB_UNASSIGNED;

            at_1->GetUDData(hb_type_udd_handle, hb_type_1);
            at_2->GetUDData(hb_type_udd_handle, hb_type_2);

            bool match = false;
            bool swap = false;
            bool ligand_atom_is_donor = true;

            if ((hb_type_1 == coot::HB_DONOR    || hb_type_1 == coot::HB_BOTH) &&
                (hb_type_2 == coot::HB_ACCEPTOR || hb_type_2 == coot::HB_BOTH)) {
               match = true;
            }

            if ((hb_type_1 == coot::HB_ACCEPTOR || hb_type_1 == coot::HB_BOTH) &&
                (hb_type_2 == coot::HB_DONOR    || hb_type_2 == coot::HB_BOTH)) {
               match = true;
               swap = true;
               ligand_atom_is_donor = false;
            }

            if (false) {
               std::cout << "atom 1st selection: " << coot::atom_spec_t(at_1) << " has hb_type "
                         << hb_type_1 << std::endl;
               std::cout << "atom 2nd selection: " << coot::atom_spec_t(at_2) << " has hb_type "
                         << hb_type_2 << std::endl;
               std::cout << ".... makes match: " << match << std::endl;
            }

            if (match) {

               if (false)
                  std::cout << "MATCH (pre-swap) H-bond donor acceptor match "
                            << coot::atom_spec_t(at_1) << " " << at_1->GetResName()
                            << " to " << coot::atom_spec_t(at_2) << " " << at_2->GetResName() << "\n";

               if (swap)
                  std::swap(at_1, at_2);

               // donor first: xxx_1
               // 
               std::vector<std::pair<mmdb::Atom *, float> > nm_1 = neighbour_map[at_1];
               std::vector<std::pair<mmdb::Atom *, float> > nm_2 = neighbour_map[at_2];
               
               std::string res_type_1 = at_1->GetResName();
               std::string res_type_2 = at_2->GetResName();

               double angle_1_0 = -1; // assigned later
               double angle_2_0 = -1; // assigned later
               
               // std::cout << "   #donor-neighbs: " << nm_1.size()
               // << " #acceptor-neighbs: " << nm_2.size()
               // << " res_type_1 " << res_type_1
               // << " res_type_2 " << res_type_2
               // << std::endl;

               if (((nm_1.size() > 0) || (res_type_1 == "HOH")) &&
                   ((nm_2.size() > 0) || (res_type_2 == "HOH"))) {

                  // only check the angles if there is not already a
                  // bond in the list between these atoms.
                  // 
                  coot::h_bond hb(at_1, at_2);
                  std::vector<coot::h_bond>::const_iterator it = std::find(v.begin(), v.end(), hb);
                  if (it == v.end()) {

                     bool neighbour_angles_are_good = 1;  // initially
                     double dist = coot::distance(at_1, at_2);
                     for (unsigned int ii=0; ii<nm_1.size(); ii++) {

                        double angle_1 = coot::angle(nm_1[ii].first, at_1, at_2);

                        if (angle_1_0 < 0) {  // as yet unset
                           coot::residue_spec_t donor_res_spec(at_1->GetResidue());
                           coot::residue_spec_t neigh_res_spec(nm_1[ii].first->GetResidue());
                           if (donor_res_spec == neigh_res_spec)
                              angle_1_0 = angle_1;
                        }

                        if (angle_1 < 90) {
                           neighbour_angles_are_good = 0;
                           if (debug)
                              std::cout << "Rejecting (donor) "
                                        << coot::atom_spec_t(at_1) << "..."
                                        << coot::atom_spec_t(at_2) << "  because angle "
                                        << coot::atom_spec_t(nm_1[ii].first) << "-"
                                        << coot::atom_spec_t(at_1) << "-"
                                        << coot::atom_spec_t(at_2) << " is " << angle_1
                                        << "\n";

                           break;
                        }
                     }

                     // We don't loop the second angle.  Is that a bad thing?

                     double angle_2 = -1;
                     if (nm_2.size()) // not a HOH
                        angle_2 = coot::angle(at_1, at_2, nm_2[0].first);

                     if (angle_2_0 < 0) {  // as yet unset
                        angle_2_0 = angle_2;
                     }

                     // Reject if bad angles (but don't reject if the
                     // bonded residue is a HOH (in that case there is
                     // no angle_2)).
                     //
                     if ((angle_2 < 90) && ((res_type_1 != "HOH") && (res_type_2 != "HOH"))) {
                        neighbour_angles_are_good = false;
                        if (debug)
                           std::cout << "Rejecting (acceptor) "
                                     << coot::atom_spec_t(at_1) << "..."
                                     << coot::atom_spec_t(at_2) << "  because angle "
                                     << coot::atom_spec_t(at_1) << "-"
                                     << coot::atom_spec_t(at_2)
                                     << coot::atom_spec_t(nm_2[0].first)
                                     << " is " << angle_2
                                     << "\n";
                     }
                  
                     if (neighbour_angles_are_good) {

                        if (nm_1.size())
                           hb.donor_neigh = nm_1[0].first;
                        if (nm_2.size())
                           hb.acceptor_neigh = nm_2[0].first;
                        hb.dist = dist;
                        hb.angle_1 = angle_1_0;
                        hb.angle_2 = angle_2_0;
                        hb.ligand_atom_is_donor = ligand_atom_is_donor; // selHnd_1 is presumed "ligand"

                        bool add_it = true;
                        std::string at_name_1(at_1->GetAtomName());
                        std::string at_name_2(at_1->GetAtomName());
                        std::string ch_id_1 = std::string(at_1->GetChainID());
                        std::string ch_id_2 = std::string(at_2->GetChainID());
                        std::string rn_1 = std::to_string(at_1->GetSeqNum());
                        std::string rn_2 = std::to_string(at_2->GetSeqNum());
                        std::string key_1 = at_name_1 + ch_id_1 + rn_1 + at_name_2 + ch_id_2 + rn_2;
                        std::string key_2 = at_name_2 + ch_id_2 + rn_2 + at_name_1 + ch_id_1 + rn_1;
                        if (existing_h_bond_keys.find(key_1) != existing_h_bond_keys.end()) add_it = false;
                        if (existing_h_bond_keys.find(key_2) != existing_h_bond_keys.end()) add_it = false;
                        if (add_it) {
                           existing_h_bond_keys.insert(key_1);
                           existing_h_bond_keys.insert(key_2);
                        }

                        // Don't add h-bonds that somehow have not had
                        // donor_neigh and acceptor_neigh assigned.
                        //
                        // It's OK for donor_neigh to be null/unset if
                        // the donor is an HOH.
                        //
                        // Otherwise, do add it.
                        // 
                        if (false)
                           std::cout << "Adding an H-bond " << coot::atom_spec_t(hb.donor)
                                     << " to " << coot::atom_spec_t(hb.acceptor) << std::endl;

                        if (add_it) v.push_back(hb);

                     }
                  }
               }
            } else {
               // std::cout << "No H-bond donor acceptor match " << coot::atom_spec_t(at_1)
               // << " to " << coot::atom_spec_t(at_2) << "\n";
            }
         }
      }
   }

   std::sort(v.begin(), v.end());

   if (false) {
      std::cout << "returning these h bonds: " << std::endl;
      for (unsigned int i=0; i<v.size(); i++)
         std::cout << "   " << i << "  " << v[i] << std::endl;
   }

   return v;
}

// Use hydrogen->acceptor distances, not donor and acceptor distances,
// i.e. model (both selections) has full hydrogens.
//
// ligand is the selHnd_1 and selHnd_2 is everything (including ligand (usually) AFAICS).
//
//
//          H
//         /
//        /              A
//   DD--D                \                        .
//        \                \                       .
//         \               AA
//         DD
//
// D-H-A  > 90 degrees
// H-A-AA > 90 degrees
// D-A-AA > 90 degrees
// H-A < 2.5A
// D-A < 3.9A
//
// Do not consider internal hydrogen bonds (check that the residues
// are different in H-bond analysis).
//
std::vector<coot::h_bond>
coot::h_bonds::get_mcdonald_and_thornton(int selHnd_1, int selHnd_2, mmdb::Manager *mol,
                                         const protein_geometry &geom, int imol,
                                         mmdb::realtype max_dist) {
   bool debug = false;
   std::vector<coot::h_bond> v;
   // (and mark HB hydrogens too)

   int hb_type_udd_handle = mark_donors_and_acceptors(selHnd_1, selHnd_2, mol, geom, imol); // using UDD data

   // These distance are from the acceptor to the H - not the donor
   mmdb::realtype min_dist = 0.1; // H-bonds are longer than this

   mmdb::PPAtom sel_1_atoms = 0;
   mmdb::PPAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mmdb::mat44 my_matt;
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   mol->GetSelIndex   (selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex   (selHnd_2, sel_2_atoms, n_sel_2_atoms);

   mol->SeekContacts(sel_1_atoms, n_sel_1_atoms,
                     sel_2_atoms, n_sel_2_atoms,
                     min_dist, max_dist,
                     0, // seqDist 0 -> also in same res.
                     pscontact, n_contacts,
                     0, &my_matt, i_contact_group);

   if (n_contacts > 0) {
      if (pscontact) {

         // What is the nearest neighbour of the atoms in mol?
         //
         std::map<mmdb::Atom *, std::vector<std::pair<mmdb::Atom *, float> > > neighbour_map =
            make_neighbour_map(selHnd_1, selHnd_2, mol);

         for (int i_contact=0; i_contact<n_contacts; i_contact++) {
            mmdb::Atom *at_1 = sel_1_atoms[pscontact[i_contact].id1];
            mmdb::Atom *at_2 = sel_2_atoms[pscontact[i_contact].id2];

            // move on if these are interacting atoms
            std::string alt_conf_1 = at_1->altLoc;
            std::string alt_conf_2 = at_2->altLoc;
            if (!alt_conf_1.empty() && ! alt_conf_2.empty())
               if (alt_conf_1 != alt_conf_2)
                  continue;

            if (at_1->residue != at_2->residue) {

               // are they HB_HYDROGEN and HB_ACCEPTOR?
               //
               int hb_type_1 = -999; // was coot::HB_UNASSIGNED;
               int hb_type_2 = -999; // was coot::HB_UNASSIGNED;

               at_1->GetUDData(hb_type_udd_handle, hb_type_1);
               at_2->GetUDData(hb_type_udd_handle, hb_type_2);

               if (debug) // checking this? Are the types HB_UNASSIGNED?
                  std::cout << "DEBUG:: in get_mcdonald_and_thornton() "
                            << coot::atom_spec_t(at_1) << " "
                            << coot::atom_spec_t(at_2) << "   "
                            << hb_type_1 << " " << hb_type_2 << " vs HB_UNASSIGNED: "
                            << coot::HB_UNASSIGNED << std::endl;

               // hydrogen on ligand
               //
               if (hb_type_1 == coot::HB_HYDROGEN) {
                  if (hb_type_2 == coot::HB_ACCEPTOR ||
                      hb_type_2 == coot::HB_BOTH) {

                     std::vector<std::pair<mmdb::Atom *, float> > nb_1 = neighbour_map[at_1];
                     std::vector<std::pair<mmdb::Atom *, float> > nb_2 = neighbour_map[at_2];
                     std::pair<bool, coot::h_bond> b_hbond =
                        make_h_bond_from_ligand_hydrogen(at_1, at_2, nb_1, nb_2);
                     if (b_hbond.first)
                        v.push_back(b_hbond.second);
                  }
               }

               // hydrogen on environment (protein) residue
               //
               // Allow a special alternative case where the acceptor
               // is on the ligand and the donor is a water (because
               // waters may not (probably do not) have hydrogens.
               //
               if (hb_type_1 == coot::HB_ACCEPTOR ||
                   hb_type_1 == coot::HB_BOTH) {
                  if (hb_type_2 == coot::HB_HYDROGEN || std::string(at_2->GetResName()) == "HOH") {

                     std::vector<std::pair<mmdb::Atom *, float> > nb_1 = neighbour_map[at_1];
                     std::vector<std::pair<mmdb::Atom *, float> > nb_2 = neighbour_map[at_2];
                     std::pair<bool, coot::h_bond> b_hbond =
                        make_h_bond_from_environment_residue_hydrogen(at_1, at_2, nb_1, nb_2);

                     if (b_hbond.first) {
                        if (debug)
                           std::cout << "DEBUG:: ===> in get_m&d: pushing back b_hbond "
                                     << b_hbond.second << std::endl;
                        v.push_back(b_hbond.second);
                     } else {
                        if (debug)
                           std::cout << "DEBUG:: reject" << std::endl;
                     }
                  }
               }
            }
         }
      }
   }
   return v;
}

// Shared McDonald & Thornton assessment for a putative H-bond H...A, given the
// heavy-atom neighbours of H (the donor D candidates) and of A (the acceptor
// antecedents AA). Following the "augment, don't filter" policy:
//
//   gate (first of the pair):  some D with D...A < 3.9 A, and D-H...A >= 90 deg
//   recorded, not gated:       H...A distance, H...A-AA and D...A-AA (the
//                              smallest over the heavy antecedents), with a
//                              geometry_warnings entry for each of M&T's
//                              limits that is not met (H...A > 2.5, angle < 90)
//
// The acceptor-side angles used to reject the bond. They now only warn, because
// they depend on hydrogen positions that Coot placed by rule, and a phosphate
// oxygen accepting a backbone N-H at 81 deg (FAD in 5S9S) is still worth
// reporting.
//
bool
coot::h_bonds::assess_mcdonald_and_thornton(h_bond &bond, mmdb::Atom *H, mmdb::Atom *A,
                                            const std::vector<mmdb::Atom *> &donors,
                                            const std::vector<mmdb::Atom *> &antecedents) const {

   const double D_A_max   = 3.9;  // McDonald & Thornton
   const double H_A_ideal = 2.5;  // M&T's H...A limit - a warning, not a gate
   const double angle_min = 90.0;

   bond.dist = coot::distance(H, A);

   // the donor: the heavy neighbour of H with the best D-H...A angle
   mmdb::Atom *D = nullptr;
   double best_angle_1 = -1.0;
   for (unsigned int i=0; i<donors.size(); i++) {
      double a = coot::angle(donors[i], H, A);
      if (a > best_angle_1) {
         best_angle_1 = a;
         D = donors[i];
      }
   }
   if (! D) return false;
   bond.donor = D;
   bond.angle_1 = best_angle_1;
   bond.donor_acceptor_dist = coot::distance(D, A);

   bool ok = (bond.donor_acceptor_dist < D_A_max) && (bond.angle_1 >= angle_min);

   // the acceptor antecedents: report the least favourable angles
   bond.angle_2 = -1.0;
   bond.angle_3 = -1.0;
   for (unsigned int i=0; i<antecedents.size(); i++) {
      double a2 = coot::angle(H, A, antecedents[i]);
      double a3 = coot::angle(D, A, antecedents[i]);
      if (bond.angle_2 < 0.0 || a2 < bond.angle_2) {
         bond.angle_2 = a2;
         bond.acceptor_neigh = antecedents[i];
      }
      if (bond.angle_3 < 0.0 || a3 < bond.angle_3)
         bond.angle_3 = a3;
   }

   auto warn = [&bond] (const std::string &what, double value, double limit, bool too_small) {
      std::ostringstream ss;
      ss.precision(too_small ? 0 : 2);
      ss << std::fixed << what << " " << value << (too_small ? " < " : " > ") << limit;
      bond.geometry_warnings.push_back(ss.str());
   };
   if (bond.dist > H_A_ideal)                       warn("H...A distance", bond.dist, H_A_ideal, false);
   if (bond.angle_2 >= 0.0 && bond.angle_2 < angle_min) warn("H...A-AA angle", bond.angle_2, angle_min, true);
   if (bond.angle_3 >= 0.0 && bond.angle_3 < angle_min) warn("D...A-AA angle", bond.angle_3, angle_min, true);
   return ok;
}

// heavy atoms only: the 1.8 A neighbour search also picks up sibling hydrogens
// (H-N-H, H-C-H) and an acceptor hydroxyl's own hydrogen, which are not M&T's
// D or AA.
std::vector<mmdb::Atom *>
coot::h_bonds::heavy_neighbours(const std::vector<std::pair<mmdb::Atom *, float> > &nb) {
   std::vector<mmdb::Atom *> v;
   for (unsigned int i=0; i<nb.size(); i++)
      if (! is_hydrogen_atom(nb[i].first))
         v.push_back(nb[i].first);
   return v;
}

// return an h_bond and whether it passes the gate (D...A < 3.9 A, D-H...A >= 90).
//
std::pair<bool, coot::h_bond>
coot::h_bonds::make_h_bond_from_ligand_hydrogen(mmdb::Atom *at_1, // H on ligand
                                                mmdb::Atom *at_2, // acceptor on residue
                                                const std::vector<std::pair<mmdb::Atom *, float> > &nb_1,
                                                const std::vector<std::pair<mmdb::Atom *, float> > &nb_2) const {

   coot::h_bond bond(at_1, at_2, 1); // ligand atom is Hydrogen
   bool ok = assess_mcdonald_and_thornton(bond, at_1, at_2, heavy_neighbours(nb_1), heavy_neighbours(nb_2));
   return std::pair<bool, coot::h_bond> (ok, bond);
}

// return an h_bond and whether it passes the gate.
//
// either at_1 is acceptor on the ligand
//        at_2 is H atom on the residue
// or     at_1 is acceptor on the ligand
//        at_2 is O of water, standing in for its (unmodelled) hydrogen:
//        then the gate is O...A < water_dist_max and the acceptor-side
//        angles are measured from the water O.
//
std::pair<bool, coot::h_bond>
coot::h_bonds::make_h_bond_from_environment_residue_hydrogen(mmdb::Atom *at_1, // acceptor on ligand
                                                             mmdb::Atom *at_2, // H on residue (or water O)
                                                             const std::vector<std::pair<mmdb::Atom *, float> > &nb_1,
                                                             const std::vector<std::pair<mmdb::Atom *, float> > &nb_2) const {

   double water_dist_max = 3.25; // pass this

   bool ligand_atom_is_H_flag = false;
   h_bond bond(at_2, at_1, ligand_atom_is_H_flag); // H atom goes first for this constructor

   if (std::string(at_2->GetResName()) == "HOH") {
      bond.dist = coot::distance(at_1, at_2);
      bond.donor = at_2;
      bond.acceptor = at_1;
      bond.donor_acceptor_dist = bond.dist;
      bond.angle_1 = -1.0; // no hydrogen
      bond.angle_2 = -1.0;
      bond.angle_3 = -1.0;
      std::vector<mmdb::Atom *> antecedents = heavy_neighbours(nb_1);
      for (unsigned int i=0; i<antecedents.size(); i++) {
         double a = coot::angle(at_2, at_1, antecedents[i]);
         if (bond.angle_2 < 0.0 || a < bond.angle_2) {
            bond.angle_2 = a;
            bond.angle_3 = a;
            bond.acceptor_neigh = antecedents[i];
         }
      }
      if (bond.angle_2 >= 0.0 && bond.angle_2 < 90.0) {
         std::ostringstream ss;
         ss.precision(0);
         ss << std::fixed << "O(water)...A-AA angle " << bond.angle_2 << " < 90";
         bond.geometry_warnings.push_back(ss.str());
      }
      bool ok = (bond.dist < water_dist_max);
      return std::pair<bool, h_bond> (ok, bond);
   }

   bool ok = assess_mcdonald_and_thornton(bond, at_2, at_1, heavy_neighbours(nb_2), heavy_neighbours(nb_1));
   bond.acceptor = at_1;
   return std::pair<bool, h_bond> (ok, bond);
}





// using UDD data
//
// return the UDD handle
int
coot::h_bonds::mark_donors_and_acceptors(int selHnd_1, int selHnd_2, mmdb::Manager *mol,
                                         const coot::protein_geometry &geom, int imol) {

   bool debug = false;
   mmdb::PPAtom sel_1_atoms = 0;
   mmdb::PPAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mol->GetSelIndex(selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex(selHnd_2, sel_2_atoms, n_sel_2_atoms);
   int udd_h_bond_type_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "hb_type");

   for (int i=0; i<n_sel_1_atoms; i++) {
      std::string name = sel_1_atoms[i]->name;
      std::string res_name = sel_1_atoms[i]->GetResName();

      int h_bond_type = geom.get_h_bond_type(name, res_name, imol);

      sel_1_atoms[i]->PutUDData(udd_h_bond_type_handle, h_bond_type);
      if (debug)
         std::cout << "   h_bonds:: "
                   << sel_1_atoms[i]->GetChainID() << " "
                   << sel_1_atoms[i]->GetSeqNum()  << " "
                   << sel_1_atoms[i]->GetResName() << " "
                   << "name: " << name << " marked as " << h_bond_type << "\n";
   }

   if (selHnd_1 != selHnd_2) {
      for (int i=0; i<n_sel_2_atoms; i++) {
         std::string name = sel_2_atoms[i]->name;
         std::string res_name = sel_2_atoms[i]->GetResName();
         int h_bond_type = geom.get_h_bond_type(name, res_name, imol);
         sel_2_atoms[i]->PutUDData(udd_h_bond_type_handle, h_bond_type);
      if (debug)
         std::cout << "h_bonds:: "
                   << sel_2_atoms[i]->GetChainID() << " "
                   << sel_2_atoms[i]->GetSeqNum() << " "
                   << sel_2_atoms[i]->GetResName() << " "
                   << "name: " << name << " marked as " << h_bond_type << "\n";
      }
   }

   return udd_h_bond_type_handle;
}


// What is the nearest neighbour of the atoms in mol?
//
std::map<mmdb::Atom *, std::vector<std::pair<mmdb::Atom *, float> > >
coot::h_bonds::make_neighbour_map(int selHnd_1, int selHnd_2, mmdb::Manager *mol) {

   std::map<mmdb::Atom *, std::vector<std::pair<mmdb::Atom *, float> > > atom_map;
   mmdb::PPAtom sel_1_atoms = 0;
   mmdb::PPAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mol->GetSelIndex   (selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex   (selHnd_2, sel_2_atoms, n_sel_2_atoms);

   mmdb::mat44 my_matt;
   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;

   mol->SeekContacts(sel_1_atoms, n_sel_1_atoms,
                     sel_1_atoms, n_sel_1_atoms,
                     0.1, 1.8,
                     0, // seqDist 0 -> also in same res.
                     pscontact, n_contacts,
                     0, &my_matt, i_contact_group);

   if (n_contacts) {
      if (pscontact) {
         for (int i_contact=0; i_contact<n_contacts; i_contact++) {
            clipper::Coord_orth pt_1(sel_1_atoms[pscontact[i_contact].id1]->x,
                                     sel_1_atoms[pscontact[i_contact].id1]->y,
                                     sel_1_atoms[pscontact[i_contact].id1]->z);
            clipper::Coord_orth pt_2(sel_1_atoms[pscontact[i_contact].id2]->x,
                                     sel_1_atoms[pscontact[i_contact].id2]->y,
                                     sel_1_atoms[pscontact[i_contact].id2]->z);
            float d = clipper::Coord_orth::length(pt_1, pt_2);
            coot::residue_spec_t res_1(sel_1_atoms[pscontact[i_contact].id1]->GetResidue());
            coot::residue_spec_t res_2(sel_1_atoms[pscontact[i_contact].id2]->GetResidue());

            // neighbours of donor (or acceptors) have to be in the
            // same residue as the donor (or acceptor).
            //
            if (res_1 == res_2) {
               std::pair<mmdb::Atom *, float> p(sel_1_atoms[pscontact[i_contact].id2], d);
               atom_map[sel_1_atoms[pscontact[i_contact].id1]].push_back(p);
            }
         }
      }
   }
   delete [] pscontact;
   pscontact = NULL;


   mol->SeekContacts(sel_2_atoms, n_sel_2_atoms,
                     sel_2_atoms, n_sel_2_atoms,
                     0.1, 1.7,
                     0, // seqDist 0 -> also in same res.
                     pscontact, n_contacts,
                     0, &my_matt, i_contact_group);

   if (n_contacts) {
      if (pscontact) {
         for (int i_contact=0; i_contact<n_contacts; i_contact++) {
            clipper::Coord_orth pt_1(sel_2_atoms[pscontact[i_contact].id1]->x,
                                     sel_2_atoms[pscontact[i_contact].id1]->y,
                                     sel_2_atoms[pscontact[i_contact].id1]->z);
            clipper::Coord_orth pt_2(sel_2_atoms[pscontact[i_contact].id2]->x,
                                     sel_2_atoms[pscontact[i_contact].id2]->y,
                                     sel_2_atoms[pscontact[i_contact].id2]->z);
            coot::residue_spec_t res_1(sel_2_atoms[pscontact[i_contact].id1]->GetResidue());
            coot::residue_spec_t res_2(sel_2_atoms[pscontact[i_contact].id2]->GetResidue());

            // neighbours of donor (or acceptors) have to be in the
            // same residue as the donor (or acceptor).
            //
            if (res_1 == res_2) {

               float d = clipper::Coord_orth::length(pt_1, pt_2);
               std::pair<mmdb::Atom *, float> p(sel_2_atoms[pscontact[i_contact].id2], d);

               // only add p if is not already in the atom map vector for this atom:
               // (this relies on the doubles matching :) but it seems to work...
               //
               std::vector<std::pair<mmdb::Atom *, float> >::const_iterator it =
                  std::find(atom_map[sel_2_atoms[pscontact[i_contact].id1]].begin(),
                            atom_map[sel_2_atoms[pscontact[i_contact].id1]].end(), p);
               if (it == atom_map[sel_2_atoms[pscontact[i_contact].id1]].end())
                  atom_map[sel_2_atoms[pscontact[i_contact].id1]].push_back(p);
            }
         }
      }
   }
   delete [] pscontact;
   pscontact = NULL;
   

   // sort the neighbour by distance 
   //
   std::map<mmdb::Atom *, std::vector<std::pair<mmdb::Atom *, float> > >::iterator it;
   for (it=atom_map.begin(); it != atom_map.end(); it++) {

      //  You can't use a const_iterator to sort this vector's innards.
      // 
      coot::h_bonds::atom_sorter as(it->first);
      std::sort(it->second.begin(), it->second.end(), as);
   }

   bool debug = false;

   // were they sorted correctly?  Debug
   //
   if (debug) {
      for (it=atom_map.begin(); it != atom_map.end(); it++) {
         std::cout << coot::atom_spec_t(it->first) << "\n";
         for (unsigned int i=0; i<it->second.size(); i++) {
            std::cout << "    DEUBG:: " << coot::atom_spec_t(it->second[i].first) << " "
                      << it->second[i].second << "\n";
         }
      }
   }

   return atom_map;
}


// check that some (formally, at least one) of the atoms have a
// defined HB status (energy_lib_atom hb_t).
// 
// return the hb_type_udd_handle as second.
// 
std::pair<bool, int>
coot::h_bonds::check_hb_status(int selhnd, mmdb::Manager *mol, const protein_geometry &geom, int imol) {

   bool status = false;
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;

   int hb_type = HB_UNASSIGNED;

   int hb_type_udd_handle = mark_donors_and_acceptors(selhnd, -1, mol, geom, imol); // using UDD data

   mol->GetSelIndex(selhnd, residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      at->GetUDData(hb_type_udd_handle, hb_type);
      if (0)
         std::cout << "   " << atom_spec_t(at) << " " << hb_type << std::endl;
      if (hb_type != HB_UNASSIGNED)
         status = true;
   }
   return std::pair<bool, int> (status, hb_type_udd_handle);
}
