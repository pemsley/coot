/*
 * ideal/ng.cc
 * 
 * Copyright 2019 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#include <iomanip>

#include "compat/coot-sysdep.h"

#include "coot-utils/bonded-pairs.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "simple-restraint.hh"
#include "coot-utils/contacts-by-bricks.hh"
#include "coot-utils/stack-and-pair.hh"

#include "utils/logging.hh"
extern logging logger;


int
coot::restraints_container_t::make_restraints_ng(int imol,
                                                 const coot::protein_geometry &geom,
                                                 coot::restraint_usage_Flags flags_in,
                                                 bool do_residue_internal_torsions,
                                                 bool do_trans_peptide_restraints,
                                                 float rama_plot_target_weight,
                                                 bool do_rama_plot_restraints,
                                                 bool do_auto_helix_restraints,
                                                 bool do_auto_strand_restraints,
                                                 bool do_auto_h_bond_restraints,
                                                 coot::pseudo_restraint_bond_type sec_struct_pseudo_bonds,
                                                 bool do_link_restraints,
                                                 bool do_flank_restraints) {


   // std::cout << "..................     make_restraints_ng() --- start --- n_atoms: " << n_atoms << std::endl;

   // debug_sets();

   bool console_output_for_restraints_generation_timings = false;

   if (! thread_pool_p) {
         std::cout << "ERROR:: " << __FUNCTION__ << " --- thread pool was not set! ---------"
                   << std::endl;
         std::cout << "ERROR:: make_restraints_ng() stops before it starts" << std::endl;
         return -1;
   }

   auto tp_0 = std::chrono::high_resolution_clock::now();
   restraints_usage_flag = flags_in;
   bool just_chirals = false;
   if (restraints_usage_flag == coot::CHIRAL_VOLUMES) just_chirals = true;
   rama_plot_weight = rama_plot_target_weight;
   if (n_atoms > 0) {
      mark_OXT(geom);
      make_monomer_restraints(imol, geom, do_residue_internal_torsions);
      auto tp_1 = std::chrono::high_resolution_clock::now();

      // This should be a set, not a vector, i.e.
      // std::map<mmdb::Residue *, std::set<mmdb::Residue *> > residue_link_set_map;

      // Fix that in the lab.
      std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > residue_link_vector_map;

      // This should be a trivial class (that also contains the link type)
      // not a pair (you can keep the variable name though)

#if 0
      // Heres's a clue:  We should add them both ways when they are made,
      // so it's cleaner to check if they are there later on in the code (say,
      // rama code)
      //
      class linked_residue_pair {
      public:
         mmdb::Residue *r1;
         mmdb::Residue *r2;
         std::string link_type;
         linked_residue_pair(mmdb::Residue *r1_in, mmdb::Residue *r2_in) : r1(r1_in), r2(r2_in) {}
         linked_residue_pair(mmdb::Residue *r1_in, mmdb::Residue *r2_in, const std::string &s) : r1(r1_in), r2(r2_in), link_type(s) {}
         bool match_p(mmdb::Residue *test_pair_1, mmdb::Residue *test_pair_2) const {
            if (r1 == test_pair_1)
               if (r2 == test_pair_2)
                  return true;
            return false;
         }
         // this class will be used in a std::map, so it needs operator==() and
         // operator<()
         bool operator==(const linked_residue_pair &lrp) {
            return match_p(lrp.r1, lrp.r2);
         }
         bool operator<(const linked_residue_pair &lrp) {
            return (lrp.r1 < r1);
         }
      };
#endif

      std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > residue_pair_link_set;

      // make_link_restraints_ng makes flanking residues restraints also
      if (residues_vec.size() > 0) // always true?
         make_link_restraints_ng(geom,
                                 do_rama_plot_restraints, do_trans_peptide_restraints,
                                 &residue_link_vector_map, // fill this
                                 &residue_pair_link_set);  // fill this

      auto tp_2 = std::chrono::high_resolution_clock::now();

      auto tp_3 = std::chrono::high_resolution_clock::now();

      if (! just_chirals)
         raic.init(restraints_vec);

      auto tp_4 = std::chrono::high_resolution_clock::now();

      // the non-bonded contact atoms at the end of the atoms array
      // don't have (forward) neighbours, so we don't want to
      // calculate NBC restraints for them (and don't use
      // them when splitting non-bonded contacts into range sets).
      //

      if (! just_chirals)
         non_bonded_contacts_atom_indices.resize(n_atoms_limit_for_nbc);

      // the non-threaded version has a different limit on the
      // non_bonded_contacts_atom_indices (so, out of range if you use it?)

      if (! thread_pool_p) {
         std::cout << "ERROR:: ------- " << __FUNCTION__ << " - thread pool was not set! ---------"
                   << std::endl;
         // and yet we continue... that's bad news.
         std::cout << "ERROR:: Bad things will now happen" << std::endl;
      }

      if (! just_chirals)
         make_non_bonded_contact_restraints_using_threads_ng(imol, geom);

      auto tp_5 = std::chrono::high_resolution_clock::now();

      if (do_rama_plot_restraints)
         make_rama_plot_restraints_ng(residue_link_vector_map, residue_pair_link_set, geom);

      if (console_output_for_restraints_generation_timings) {
         auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
         auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
         auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
         auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
         auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
         std::cout << "------------------ timings: for make_restraints_ng(): monomers: "
                   << d10 << " links: " << d21 << " flank: " << d32 << " raic: " << d43 << " nbc: " << d54
                   << " milliseconds " << std::endl;
      }

      // probably not needed because plane restraints don't come in a bunch (I don't know why)
      // disperse_plane_restraints();

      if (! just_chirals) {
         if (residues_vec.size() > 1)
            if (sec_struct_pseudo_bonds == coot::HELIX_PSEUDO_BONDS)
               make_helix_pseudo_bond_restraints();

         if (residues_vec.size() > 1)
            if (sec_struct_pseudo_bonds == coot::STRAND_PSEUDO_BONDS)
               make_strand_pseudo_bond_restraints();

         if (residues_vec.size() > 1)
            if (do_auto_helix_restraints)
               make_helix_pseudo_bond_restraints_from_res_vec_auto();

         if (residues_vec.size() > 1)
            if (do_auto_h_bond_restraints)
               make_h_bond_restraints_from_res_vec_auto(geom, imol);

         if (residues_vec.size() > 1)
            make_base_pairing_and_stacking_restraints_ng(imol, geom);
      }

      make_df_restraints_indices();
      make_distortion_electron_density_ranges();

      // Now we don't do this here - we do it after all the restraints have been added -
      // called from graphics-info-modelling.cc
      //
      // std::cout << "DEBUG:: make-restraints(): analysis of bad geometry in input model" << std::endl;
      // analyze_for_bad_restraints(); // bonds and non-bonded.

      // info();  - are the NBCs correct?

   }

   for (unsigned int i=0; i<restraints_vec.size(); i++)
      restraints_vec[i].restraint_index = i;
   return size();
}


void
coot::restraints_container_t::make_rama_plot_restraints_ng(const std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > &residue_link_vector_map,
                                                           const std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > &residue_pair_link_set,
                                                           const protein_geometry &geom) {

   bool debug = false;
   // std::cout << "############################ debug #########################  " << std::endl;
   // std::cout << "################# make_rama_plot_restraints_ng() ###########  " << std::endl;
   std::map<mmdb::Residue *, std::vector<mmdb::Residue *> >::const_iterator it_map;

   if (debug) {
      for (it_map=residue_link_vector_map.begin(); it_map != residue_link_vector_map.end(); ++it_map) {
         mmdb::Residue *key = it_map->first;
         std::cout << "key " << residue_spec_t(key) << " link vector map: " << std::endl;
         const std::vector<mmdb::Residue *> &v = it_map->second;
         for (unsigned int i=0; i<v.size(); i++) {
            std::cout << "    " << i << " " << residue_spec_t(v[i]) << std::endl;
         }
      }
      std::set<std::pair<mmdb::Residue *, mmdb::Residue *> >::const_iterator it_set;
      for (it_set=residue_pair_link_set.begin(); it_set!=residue_pair_link_set.end(); ++it_set) {
         mmdb::Residue *r_1 = it_set->first;
         mmdb::Residue *r_2 = it_set->second;
         std::cout << " paired: " << residue_spec_t(r_1) << " " << residue_spec_t(r_2) << std::endl;
      }
   }

   // of course residues_vec is sorted by now

   int n_vec_rama_residues = residues_vec.size() -1; // the residue at the end will not have the atoms
                                                       // for the next residue

   if (residues_vec.size() > 2) {
      for (int i=1; i<n_vec_rama_residues; i++) {
         mmdb::Residue *residue_prev_p = residues_vec[i-1].second;
         mmdb::Residue *residue_this_p = residues_vec[i  ].second;
         mmdb::Residue *residue_next_p = residues_vec[i+1].second;
         bool f1 = residues_vec[i-1].first; // fixed flags
         bool f2 = residues_vec[i  ].first;
         bool f3 = residues_vec[i+1].first;
         if (f1 && f2 && f3)
            continue;
         // So there are some moving atoms if we are here
         int index_p = residue_prev_p->index;
         int index_t = residue_this_p->index;
         int index_n = residue_next_p->index;
         if (false)
            std::cout << "residues "
                      << residue_spec_t(residue_prev_p) << " "
                      << residue_spec_t(residue_this_p) << " "
                      << residue_spec_t(residue_next_p) << " "
                      << "residue indices: " << index_p << " " << index_t << " " << index_n << std::endl;
         if ((index_t-index_p) == 1) {
            if ((index_n-index_t) == 1) {

               // OK so the residues are next to each other in the chain.
               // Are they protein?
               std::string rn_1(residue_prev_p->GetResName());
               std::string rn_2(residue_this_p->GetResName());
               std::string rn_3(residue_next_p->GetResName());
               if (util::is_standard_amino_acid_name(rn_1)) {
                  if (util::is_standard_amino_acid_name(rn_2)) {
                     if (util::is_standard_amino_acid_name(rn_3)) {

                        // rama restraints should not be added to CIS, PCIS. Only TRANS and PTRANS

                        // were they linked by peptide bonds?
                        mmdb::Residue *key = residue_this_p;
                        it_map = residue_link_vector_map.find(key);
                        if (it_map != residue_link_vector_map.end()) {
                           const std::vector<mmdb::Residue *> &v = it_map->second;
                           std::vector<mmdb::Residue *>::const_iterator it_vec_prev;
                           std::vector<mmdb::Residue *>::const_iterator it_vec_next;
                           it_vec_prev = std::find(v.begin(), v.end(), residue_prev_p);
                           it_vec_next = std::find(v.begin(), v.end(), residue_next_p);
                           if (it_vec_prev != v.end()) {
                              if (it_vec_next != v.end()) {
                                 // OK, they were linked by peptide bonds
                                 std::string link_type = find_peptide_link_type_ng(residue_prev_p, residue_this_p, geom);

                                 if (link_type == "TRANS" || link_type == "PTRANS") {
                                    rama_triple_t triple(residue_prev_p, residue_this_p, residue_next_p, link_type, f1, f2, f3);
                                    add_rama(triple, geom);
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

}

void
coot::restraints_container_t::make_base_pairing_and_stacking_restraints_ng(int imol, const coot::protein_geometry &geom) {

   // maybe useful for to change to DEBUG
   bool console_output_for_restraints_generation_timings = false;

   auto tp_6 = std::chrono::high_resolution_clock::now();
   stack_and_pair sp(mol, residues_vec);
   std::vector<stack_and_pair::stacked_planes_info_t> stacked_residues = sp.stacked_residues(mol);

   // using spec indirection - this is not fast
   extra_restraints_t extra_restraints;
   for (std::size_t i=0; i<stacked_residues.size(); i++) {
      parallel_planes_t ppr(residue_spec_t(stacked_residues[i].res_1),
                            residue_spec_t(stacked_residues[i].res_2),
                            stacked_residues[i].atom_names_1,
                            stacked_residues[i].atom_names_2,
                            "", "");
      extra_restraints.parallel_plane_restraints.push_back(ppr);
   }

   auto tp_7 = std::chrono::high_resolution_clock::now();

   // base pairs:
   bool all_atoms_are_moving_flag = false; // can be true in the future.
   std::vector<stack_and_pair::paired_residues_info_t> pr =
      sp.paired_residues(mol, residues_vec, all_atoms_are_moving_flag, geom, imol);

   auto tp_8 = std::chrono::high_resolution_clock::now();

   unsigned int n_base_pairing_bonds = 0;
   for (std::size_t i=0; i<pr.size(); i++) {
      for (std::size_t j=0; j<pr[i].atom_pair_vec.size(); j++) {
         mmdb::Atom *at_1 = pr[i].atom_pair_vec[j].first;
         mmdb::Atom *at_2 = pr[i].atom_pair_vec[j].second;
         int index_1 = -1;
         int index_2 = -1;
         at_1->GetUDData(udd_atom_index_handle, index_1);
         at_2->GetUDData(udd_atom_index_handle, index_2);
         std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);
         if (!fixed_flags[0] || !fixed_flags[1]) {
            // This should be hydrogen bond type (whatever that means)
            // the target distance depends on the pair and should be calculated in
            // paired_residues and made part of the atom "pair" (-> tuple/class)
            //
            // typical values: 2.92, 2.83, 2.87
            add(BOND_RESTRAINT, index_1, index_2, fixed_flags, 2.88, 0.08, 1.2);
            n_base_pairing_bonds++;
         }
      }
   }

   // std::cout << "INFO:: Made " << n_base_pairing_bonds << " base pairing Hydrogen bonds"
   //           << std::endl;
   logger.log(log_t::DEBUG, "Made", n_base_pairing_bonds, "base pairing hydrogen bonds");

   auto tp_9 = std::chrono::high_resolution_clock::now();

   if (extra_restraints.has_restraints())
      add_extra_restraints(imol, "from make_base_pairing_and_stacking_restraints_ng()", extra_restraints, geom);

   if (console_output_for_restraints_generation_timings) {
      auto tp_10 = std::chrono::high_resolution_clock::now();
      auto d76  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_7  - tp_6).count();
      auto d87  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_8  - tp_7).count();
      auto d98  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_9  - tp_8).count();
      auto d109 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_10 - tp_9).count();
      std::cout << "------------------ timings: for make_restraints_ng(): stacking and pairing: "
                << d76 << " " << d87 << " " << d98 << " " << d109 << " ms" << std::endl;
   }

}


// residue_link_count_map should be residue_link_info_map - the data of the map should be
// a vector of linked residues (we can get the index of those residues if needed)
// but more importantly perhaps, filter them out of the residues near the A and C residues
// of our (say) A-B-C residue triplet.
//
void
coot::restraints_container_t::make_flanking_atoms_restraints_ng(const coot::protein_geometry &geom,
                                                                std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > *residue_link_vector_map_p,
                                                                std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > *residue_pair_link_set_p,
                                                                bool do_rama_plot_restraints,
                                                                bool do_trans_peptide_restraints) {

   bool debug = false;

   link_restraints_counts flank_restraints("flank");

   // residue 20 is being refined, we need to make flanking link restraints to residues 19 and 21 using
   // fixed_neighbours_set
   //
   // residues_vec contains residue 20
   //
   // fixed_neighbours_set are the (fixed) neighbours of residue 20.
   //
   // mol contains residue 20 and the residues around it (made by the function that
   // creates the restraints object)

   if (debug) { // debugging
      std::cout << "########## make_flanking_atoms_restraints_ng() residue_link_count_map size "
                << residue_link_vector_map_p->size() << std::endl;
      std::map<mmdb::Residue *, std::vector<mmdb::Residue *> >::const_iterator it;
      for (it=residue_link_vector_map_p->begin();
                it!=residue_link_vector_map_p->end(); it++) {
              std::cout << "        residue_link_vector_map key: "
                             << residue_spec_t(it->first) << " " << std::endl;
      }
   }

   if (debug) {
      std::cout << "####### make_flanking_atoms_restraints_ng() debugging fixed_neighbours_set() " << std::endl;
      std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
      for (it=fixed_neighbours_set.begin(); it!=fixed_neighbours_set.end(); it++) {
              mmdb::Residue *residue_p = it->first;
              std::cout << "\n\n...fixed-neighbour for residue " << residue_spec_t(residue_p) << " index " << residue_p->index << std::endl;
              const std::set<mmdb::Residue *> &s = it->second;
              std::set<mmdb::Residue *>::const_iterator its;
              for (its=s.begin(); its!=s.end(); its++) {
                 mmdb::Residue *neighb = *its;
                 std::cout << "      neighb: " << residue_spec_t(neighb) << " " << neighb << " index " << neighb->index << std::endl;
              }
      }
      std::cout << "####### done make_flanking_atoms_restraints_ng() debugging fixed_neighbours_set() " << std::endl;
   }

   int idx_reference_index_handle = mol->GetUDDHandle(mmdb::UDR_RESIDUE, "index from reference residue");
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
   for (it=fixed_neighbours_set.begin(); it!=fixed_neighbours_set.end(); it++) {
      mmdb::Residue *residue_p = it->first;
      unsigned int n_link_for_residue = 0;
      std::map<mmdb::Residue *, std::vector<mmdb::Residue *> >::const_iterator itm = residue_link_vector_map_p->find(residue_p);
      if (itm != residue_link_vector_map_p->end()) {
         n_link_for_residue = itm->second.size();
      }

      // was it polymer-linked at both ends? (we've only made _polymer_ links so far)
      //
      if (n_link_for_residue < 2) {

         int index_for_residue;
         residue_p->GetUDData(idx_reference_index_handle, index_for_residue);
         const std::set<mmdb::Residue *> &s = it->second;
         std::set<mmdb::Residue *>::const_iterator its;
         for (its=s.begin(); its!=s.end(); its++) {
            mmdb::Residue *neighb = *its;

            // Let's say that we refine residue 21,22,23.
            // Residue 21 needs to be flanking-linked to residue 20.
            // fixed_neighbours_set can contain residue
            // 22 as a neighbour of 21. We don't want to link 22 to 21 here.
            // because we have done it in the polymer function. Let's
            // see if we can find 21-22 here in residue_link_vector_map.
            //
            // If we are refining just one residue, there will be nothing in the
            // residue_link_vector_map, so only skip this one (continue)
            // if we can find the key and the neighb is in the vector (set) of
            // residues to whicih this residue is already linked.

            if (neighb->chain != residue_p->chain) continue;

            // this kills links for insertion codes, but also kills links
            // across gaps (this is the feature that we want).
            //
            int rn_1 = neighb->GetSeqNum();
            int rn_2 = residue_p->GetSeqNum();
            int rn_delta = abs(rn_2 - rn_1);
            if (rn_delta != 1) {
               std::string ins_1(residue_p->GetInsCode());
               std::string ins_2(neighb->GetInsCode());
               if (ins_1.empty())
                  if (ins_2.empty())
                     continue;
            }

            std::map<mmdb::Residue *, std::vector<mmdb::Residue *> >::const_iterator itm;
            itm = residue_link_vector_map_p->find(residue_p);
            if (itm != residue_link_vector_map_p->end())
               if (std::find(itm->second.begin(), itm->second.end(), neighb) != itm->second.end())
                  continue;

            int index_delta = neighb->index - residue_p->index;

            if (index_delta == -1 || index_delta == 1) {
               // std::cout << "        fixed neigb: " << residue_spec_t(*its) << std::endl;

               std::pair<bool, mmdb::Residue *> pair_1(true,  neighb);
               std::pair<bool, mmdb::Residue *> pair_2(false, residue_p);

               // if this is a link onto the N of this residue, then we need
               // to call try_make_peptide_link_ng with the arguments the other way
               // around
               //
               if (index_delta == 1) std::swap(pair_1, pair_2);

               if (debug)
                  std::cout << "         making flanking-peptide link: "
                            << residue_spec_t(pair_1.second) << " fixed: " << pair_1.first << "  "
                            << residue_spec_t(pair_2.second) << " fixed: " << pair_2.first << "  "
                            << std::endl;

               std::pair<bool, link_restraints_counts> link_result =
                  try_make_peptide_link_ng(geom,
                                           pair_1,
                                           pair_2,
                                           do_trans_peptide_restraints);

               if (! link_result.first)
                  link_result = try_make_phosphodiester_link_ng(geom, pair_1, pair_2);

               if (link_result.first) {
                  flank_restraints.add(link_result.second);
                  (*residue_link_vector_map_p)[residue_p].push_back(neighb);
                  (*residue_link_vector_map_p)[neighb   ].push_back(residue_p);

                  std::pair<mmdb::Residue *, mmdb::Residue *> p1(residue_p, neighb);
                  std::pair<mmdb::Residue *, mmdb::Residue *> p2(neighb, residue_p);
                  residue_pair_link_set_p->insert(p1);
                  residue_pair_link_set_p->insert(p2);
               }
            }
         }
      }
   }

   // flank_restraints.report();
}

// static
void
coot::restraints_container_t::make_non_bonded_contact_restraints_workpackage_ng(int ithread,
                                                                                int imol,
                                                                                const coot::protein_geometry &geom,
                                                                                const std::vector<std::set<int> > &bonded_atom_indices,
                                                                                const reduced_angle_info_container_t &raic,
                                                                                const std::vector<std::set<unsigned int> > &vcontacts,
                                                                                std::pair<unsigned int, unsigned int> atom_index_range_pair,
                                                                                const std::set<int> &fixed_atom_indices,
                                                                                const std::vector<std::string> &energy_type_for_atom,
                                                                                bool extended_atom_mode,
                                                                                mmdb::PPAtom atom,
                                                                                const std::vector<bool> &atom_is_metal,
                                                                                const std::vector<bool> &atom_is_hydrogen,
                                                                                const std::vector<bool> &H_atom_parent_atom_is_donor_vec,
                                                                                const std::vector<bool> &atom_is_acceptor_vec,
                                                                                std::vector<std::set<int> > *non_bonded_contacts_atom_indices_p,
                                                                                std::vector<simple_restraint> *nbc_restraints_fragment_p,
                                                                                std::atomic<unsigned int> &done_count) {


   if (false) {
      get_print_lock();
      std::cout << "H_atom_parent_atom_is_donor_vec vector size " << H_atom_parent_atom_is_donor_vec.size() << std::endl;
      for (unsigned int i=0; i<H_atom_parent_atom_is_donor_vec.size(); i++)
         std::cout << "H_atom_parent_atom_is_donor_vec " << i << " " << H_atom_parent_atom_is_donor_vec[i] << std::endl;
      release_print_lock();
   }

   // make_fixed_flags() will have to be done in place using fixed_atom_indices

   // think about is_in_same ring without a cache. Hmm.

   // pointer to reference
   std::vector<std::set<int> > &non_bonded_contacts_atom_indices = *non_bonded_contacts_atom_indices_p;
   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map_cache;

   // bool extended_atom_mode = false; // turn this on if there are no Hydrogen atoms in the model
        // extended_atom_mode = ! model_has_hydrogen_atoms;

   // there is no need to modify a pointer. It can return a value
   //
   auto tweak_neighbours_for_rn_diff_is_1 = [] (const std::vector<bool> &fixed_atom_flags,
                                                const std::string &atom_name_1,
                                                const std::string &atom_name_2,
                                                const int &res_no_1,
                                                const int &res_no_2,
                                                bool second_is_pro,
                                                double *dist_min_p) {

                                               double &dist_min(*dist_min_p);

                                               bool strange_exception = false;
                                               if (fixed_atom_flags.size()) {
                                                  if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
                                                     if (atom_name_1 == " O  ")
                                                        if (atom_name_2 == " CA ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " CA ")
                                                        if (atom_name_2 == " O  ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " N  ")
                                                        if (atom_name_2 == " CB ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " CB ")
                                                        if (atom_name_2 == " N  ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " C  ")
                                                        if (atom_name_2 == " CB ")
                                                           strange_exception = true;
                                                  }
                                               }
                                               if (strange_exception)
                                                  dist_min = 2.7;

                                               // Strange that these are not marked as 1-4 related.  Fix here...
                                               // HA-CA-N-C can be down to ~2.4A.
                                               // HA-CA-C-N can be down to ~2.41A.
                                               if (res_no_2 > res_no_1) {
                                                  if (atom_name_1 == " C  ") {
                                                     if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
                                                        strange_exception = true;
                                                        dist_min = 2.4;
                                                     }
                                                  }
                                                  if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
                                                     if (atom_name_2 == " N  ") {
                                                        strange_exception = true;
                                                        dist_min = 2.41;
                                                     }
                                                  }

                                                  // Another 1-4:
                                                  // N-terminal hydrogens should not have NBCs to C of next residue
                                                  if (atom_name_2 == " C  ") {
                                                     if (atom_name_1 == " H1 " || atom_name_1 == " H2 " || atom_name_1 == " H3 ") {
                                                        dist_min = 2.15;
                                                     }
                                                  }

                                                  if (atom_name_1 == " N  ") {
                                                     if (atom_name_2 == " H  ") {
                                                        strange_exception = true;
                                                        dist_min = 2.4;
                                                     }
                                                  }

                                                  // CA(n) - CD(n+1) for XXX-PRO needs to be reduced by 0.8 or so
                                                  if (second_is_pro) {
                                                     if (atom_name_1 == " CA ") {
                                                        if (atom_name_2 == " CD ") {
                                                           dist_min = 3.1;
                                                        }
                                                     }
                                                  }

                                               } else {
                                                  if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
                                                     if (atom_name_2 == " C  ") {
                                                        strange_exception = true;
                                                        dist_min = 2.4;
                                                     }
                                                  }
                                                  if (atom_name_1 == " N  ") {
                                                     if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
                                                        strange_exception = true;
                                                        dist_min = 2.41;
                                                     }
                                                  }
                                                  if (atom_name_2 == " N  ") {
                                                     if (atom_name_1 == " H  ") {
                                                        strange_exception = true;
                                                        dist_min = 2.4;
                                                     }
                                                  }
                                               }
                                            };

   auto tweak_neighbours_for_rn_diff_is_2 = [] (const std::vector<bool> &fixed_atom_flags,
                                                const std::string &atom_name_1,
                                                const std::string &atom_name_2,
                                                double *dist_min_p) {

                                               double &dist_min(*dist_min_p);
                                               bool strange_exception = false;
                                               if (fixed_atom_flags.size()) {
                                                  if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
                                                     if (atom_name_1 == " C  ")
                                                        if (atom_name_2 == " N  ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " N  ")
                                                        if (atom_name_2 == " C  ")
                                                           strange_exception = true; // 3.1 would be enough

                                                     if (strange_exception)
                                                        dist_min = 2.7;
                                                  }
                                               }
                                               return dist_min;
                                            };

   // note to self: this is the multi-threaded version

   auto debug_print = [energy_type_for_atom] (const std::string &remark,
                                              int i, int j, mmdb::Atom *at_1, mmdb::Atom *at_2,
                                              const std::vector<bool> &fixed_atom_flags, double dist_min) {

                         get_print_lock();
                         std::cout << "Adding NBC: " << remark << " " << std::setw(4) << i << " " << std::setw(4) << j << " "
                                   << atom_spec_t(at_1) << " " << atom_spec_t(at_2) << " types: "
                                   << std::setw(4) << energy_type_for_atom[i] << " "
                                   << std::setw(4) << energy_type_for_atom[j]
                                   << " fixed-flags: " << fixed_atom_flags[0] << " " << fixed_atom_flags[1]
                                   << " dist-min: " << dist_min <<  "\n";
                         release_print_lock();
                      };

   for (unsigned int i=atom_index_range_pair.first; i<atom_index_range_pair.second; i++) {

      mmdb::Atom *at_1 = atom[i];
      if (! at_1) {
         std::cout << "ERROR:: make_non_bonded_contact_restraints_workpackage_ng()"
                   << " null atom at index " << i << " in range " << atom_index_range_pair.first
                   << " " << atom_index_range_pair.second << " with n_atoms (bonded_atom_indices size()) "
                   << bonded_atom_indices.size() << std::endl;
         continue;
      }
      if (at_1->isTer()) continue;

      const std::set<unsigned int> &n_set = vcontacts[i];
      std::string alt_conf_1(at_1->altLoc);
      // std::cout << "base atom: " << atom_spec_t(at_1) << std::endl;

      // std::cout << "Here with i " << i << " which has " << n_set.size() << " neighbours " << std::endl;
      std::set<unsigned int>::const_iterator it;
      for (it=n_set.begin(); it!=n_set.end(); ++it) {

         const unsigned int &j = *it;
         if (j < i) /* only add NBC one way round */
            continue;

         if (bonded_atom_indices[i].find(j) != bonded_atom_indices[i].end())
         continue;

         // for updating non-bonded contacts
         if (non_bonded_contacts_atom_indices[i].find(j) != non_bonded_contacts_atom_indices[i].end())
            continue;

         mmdb::Atom *at_2 = atom[j];
         std::string alt_conf_2(at_2->altLoc);

         {
            bool at_1_is_fixed_flag = false;
            bool at_2_is_fixed_flag = false;
            if (fixed_atom_indices.find(i) != fixed_atom_indices.end()) at_1_is_fixed_flag = true;
            if (fixed_atom_indices.find(j) != fixed_atom_indices.end()) at_2_is_fixed_flag = true;

            if (false)
               std::cout << "make_non_bonded_contact_restraints_workpackage_ng considering "
                         << atom_spec_t(at_1) << " " << atom_spec_t(at_2) << " : "
                         << fixed_atom_indices.size() << " fixed flags: " << at_1_is_fixed_flag << " " << at_2_is_fixed_flag
                         << std::endl;
         }

         if (!alt_conf_1.empty())        // alt confs don't see each other
            if (!alt_conf_2.empty())
               if (alt_conf_1 != alt_conf_2)
                  continue;

         if (fixed_atom_indices.find(i) != fixed_atom_indices.end())
            if (fixed_atom_indices.find(j) != fixed_atom_indices.end())
               continue;

         std::string res_name_1 = at_1->GetResName();
         std::string res_name_2 = at_2->GetResName();
         int res_no_1 = at_1->GetSeqNum();
         int res_no_2 = at_2->GetSeqNum();
         bool second_is_pro = false;
         if (res_name_2 == "PRO") second_is_pro = true; // residues are sorted and j > i
         if (res_name_2 == "HYP") second_is_pro = true;

         std::string element_1 = at_1->element;
         std::string element_2 = at_2->element;
         const std::string &type_1 = energy_type_for_atom[i];
         const std::string &type_2 = energy_type_for_atom[j];

         // std::vector<bool> fixed_atom_flags = make_fixed_flags(i, j);
         std::vector<bool> fixed_atom_flags(2, false);
         if (fixed_atom_indices.find(i) != fixed_atom_indices.end()) fixed_atom_flags[0] = true;
         if (fixed_atom_indices.find(j) != fixed_atom_indices.end()) fixed_atom_flags[1] = true;

         double dist_min = 3.4;

         bool in_same_residue_flag = (at_1->residue == at_2->residue);
         bool in_same_ring_flag = true;

         // part of this test is not needed.
         if (at_2->residue != at_1->residue) {
            in_same_ring_flag    = false;
            in_same_residue_flag = false;
         }

         std::string atom_name_1(at_1->GetAtomName());
         std::string atom_name_2(at_2->GetAtomName());

         if (in_same_ring_flag)
            in_same_ring_flag = is_in_same_ring(imol, at_2->residue, residue_ring_map_cache,
                                                atom_name_1, atom_name_2, geom);

         // this doesn't check 1-4 over a moving->non-moving peptide link (see comment above function)
         // because the non-moving atom doesn't have angle restraints.
         // 
         // 20191122-PE It should do now because we have bond restraints in a raic (and other fixes).
         //
         bool is_1_4_related = raic.is_1_4(i, j, fixed_atom_flags);

         if (false) { // debug 1-4s
            get_print_lock();
            std::cout << "make_non_bonded_contact_restraints_workpackage_ng() here with at_1 " << atom_spec_t(at_1)
                      << " at_2 " << atom_spec_t(at_2) << " is_1_4_related: " << is_1_4_related << std::endl;
            release_print_lock();
         }

         bool mc_atoms_tandem = false;
         bool mc_CC_atoms_tandem = false;
         if (! is_1_4_related) {
            // hacky case for C in a helix. Also N.
            // (because fixed atoms don't have angle restraints - so raic.is_1_4_related()
            // will not work)

            if (atom_name_1 == " C  ")
               if (atom_name_2 == " C  ")
                  if (at_2->residue->index - at_1->residue->index == 1) {
                     mc_atoms_tandem = true;
                     mc_CC_atoms_tandem = true;
                  }
            if (atom_name_1 == " N  ")
               if (atom_name_2 == " N  ")
                  if (at_2->residue->index - at_1->residue->index == -1)
                     mc_atoms_tandem = true;
            // down-weight CA-CA in cis-peptide bonds
            if (atom_name_1 == " CA ")
               if (atom_name_2 == " CA ")
                  if (true) // cis-peptide test
                     mc_atoms_tandem = true;
         }

         if (mc_atoms_tandem)
            is_1_4_related = true;

         if (is_1_4_related) {

            if (mc_atoms_tandem) dist_min = 2.99;

            if (in_same_ring_flag) {
               dist_min = 2.64; // was 2.7 but c.f. guanine ring distances
            } else {
               if (mc_CC_atoms_tandem) {
                  dist_min = 3.05; // in a helix, more elsewhere.
               } else {
                  dist_min = 2.82; // N-N, but captures all others, including H-H
               }
            }

            if (atom_is_hydrogen[i] && atom_is_hydrogen[j]) {
               dist_min = 2.42; // shold depend on energy type
            } else {
               if (atom_is_hydrogen[i])
                  dist_min -= 0.6;  // was 0.4
               if (atom_is_hydrogen[j])
                  dist_min -= 0.6;  // was 0.4
            }

         } else {
            // not 1-4 related
            std::pair<bool, double> nbc_dist = geom.get_nbc_dist_v2(type_1, type_2,
                                                                    element_1, element_2,
                                                                    atom_is_metal[i],
                                                                    atom_is_metal[j],
                                                                    extended_atom_mode,
                                                                    in_same_residue_flag,
                                                                    in_same_ring_flag);

            if (nbc_dist.first) {

               // In a helix O(n) is close to C(n+1), we should allow it.
               //
               bool is_O_C_1_5_related = check_for_O_C_1_5_relation(at_1, at_2);

               if (is_O_C_1_5_related) {
                  dist_min = 2.84;
               } else {

                  dist_min = nbc_dist.second;

                  if (false) debug_print("0", i, j, at_1, at_2, fixed_atom_flags, dist_min);

                  // Perhaps we don't have angle restraints to both atoms because one
                  // of the atoms is fixed (and thus miss that these have a 1-4 relationship).
                  // e.g. O(n) [moving] -> CA(n+1) [fixed]
                  //
                  // (this test will fail on insertion codes)
                  //

                  int rn_diff = abs(res_no_2 - res_no_1);

                  if (rn_diff == 1) {
                     // tweak dist_min
                     tweak_neighbours_for_rn_diff_is_1(fixed_atom_flags, atom_name_1, atom_name_2, res_no_1, res_no_2,
                                                       second_is_pro, &dist_min);
                  }
                  // this should be signed - not both ways, I think.
                  if (rn_diff == 2) {
                     // tweak dist_min
                     tweak_neighbours_for_rn_diff_is_2(fixed_atom_flags, atom_name_1, atom_name_2, &dist_min);
                  }

               }
            } else {
               // missing NBC distance so use this fallback. short/standard value
               dist_min = 2.8;
            }

         }

         if (false) debug_print("A", i, j, at_1, at_2, fixed_atom_flags, dist_min);

         bool is_H_non_bonded_contact = false;

         if (atom_is_hydrogen[i]) {
            is_H_non_bonded_contact = true;
            if (H_atom_parent_atom_is_donor_vec[i])
               if (atom_is_acceptor_vec[j])
                  dist_min -= 0.7;
         }

         if (false) debug_print("B", i, j, at_1, at_2, fixed_atom_flags, dist_min);

         if (atom_is_hydrogen[j]) {
            is_H_non_bonded_contact = true;
            if (H_atom_parent_atom_is_donor_vec[j])
               if (atom_is_acceptor_vec[i])
                  dist_min -= 0.7;
         }

         if (false) debug_print("C", i, j, at_1, at_2, fixed_atom_flags, dist_min);

         // if they are both hydrogens 1.2 + 1.2, the min dist for molprobity should be 2.1 <smiley face>
         if (atom_is_hydrogen[i])
            if (atom_is_hydrogen[j])
               dist_min -= 0.2; // aim for 2.2 so that 2.1 is a bit energetically unfavourable

         if (false) debug_print("D", i, j, at_1, at_2, fixed_atom_flags, dist_min);

         non_bonded_contacts_atom_indices[i].insert(j);
         simple_restraint::nbc_function_t nbcf = simple_restraint::LENNARD_JONES;
         // simple_restraint::nbc_function_t nbcf = simple_restraint::HARMONIC;
         simple_restraint r(NON_BONDED_CONTACT_RESTRAINT,
                            nbcf, i, j,
                            is_H_non_bonded_contact,
                            fixed_atom_flags, dist_min);
         nbc_restraints_fragment_p->push_back(r);

         if (false) debug_print("-end-", i, j, at_1, at_2, fixed_atom_flags, dist_min);
      }
   }

   done_count += 1; // atomic
}

void
coot::restraints_container_t::make_non_bonded_contact_restraints_using_threads_ng(int imol,
                                                                                  const coot::protein_geometry &geom) {

   bool console_output_for_restraints_generation_timings = false;

   auto tp_0 = std::chrono::high_resolution_clock::now();
   std::vector<std::string> energy_type_for_atom(n_atoms);
   std::vector<bool> H_atom_parent_atom_is_donor_vec(n_atoms, false);
   std::vector<bool> atom_is_acceptor_vec(n_atoms, false);

#if 0
   // bonded_atom_indices is constructed as bonds and angles are added (before calling add()).
   for (int i=0; i<n_atoms; i++) {
      const std::set<int> &s = bonded_atom_indices[i];
      std::cout << i << " : ";
      for (std::set<int>::const_iterator it=s.begin(); it!=s.end(); it++)
         std::cout << *it << " ";
      std::cout << "\n";
   }
#endif

   // needs timing test - might be slow (it isn't)
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      if (! at->isTer()) {
         std::string et = get_type_energy(imol, at, geom);
         energy_type_for_atom[i] = et;
         if (H_parent_atom_is_donor(at))
            H_atom_parent_atom_is_donor_vec[i] = true;
         if (is_acceptor(et, geom))
            atom_is_acceptor_vec[i] = true;
      }
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   if (console_output_for_restraints_generation_timings)
      std::cout << "------------------ timings: make_non_bonded_contact_restraints_ng(): energy types "
                << d10 << " milliseconds\n";

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map_cache;

   std::set<unsigned int> fixed_atom_flags_set; // signed to unsigned conversion - bleugh.
   std::set<int>::const_iterator it;
   for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); ++it)
      fixed_atom_flags_set.insert(*it);

   auto tp_2 = std::chrono::high_resolution_clock::now();
   float dist_max = 8.0;
   // I think that contacts_by_bricks should take a set of ints
   contacts_by_bricks cb(atom, n_atoms, fixed_atom_flags_set);
   cb.set_dist_max(dist_max);
   std::vector<std::set<unsigned int> > vcontacts;
   cb.find_the_contacts(&vcontacts);

   auto tp_2a = std::chrono::high_resolution_clock::now();

   if (n_threads == 0) n_threads = 1;

   // n_threads is a class data member
   std::vector<std::pair<unsigned int, unsigned int> > start_stop_pairs_vec;
   start_stop_pairs_vec.reserve(n_threads);

   // n_atoms is int, n_threads is unsigned int, what a mess
   //
   int n_per_thread = n_atoms_limit_for_nbc/static_cast<int>(n_threads);
   if ((n_per_thread * static_cast<int>(n_threads)) < n_atoms_limit_for_nbc)
      n_per_thread += 1;

   auto tp_2b = std::chrono::high_resolution_clock::now();

   for (unsigned int i=0; i<n_threads; i++) {
      unsigned int start = n_per_thread * i;
      unsigned int stop  = n_per_thread * (i+1); // test uses <
      if (stop > static_cast<unsigned int>(n_atoms_limit_for_nbc))
        stop = static_cast<unsigned int>(n_atoms_limit_for_nbc);
      std::pair<unsigned int, unsigned int> p(start, stop);
      start_stop_pairs_vec.push_back(p);
   }

   if (false) { // debugging
      std::cout << "n_per_thread " << n_per_thread << std::endl;
      std::cout << "n_atoms_limit_for_nbc " << n_atoms_limit_for_nbc << std::endl;
      for (std::size_t ii=0; ii<start_stop_pairs_vec.size(); ii++) {
         std::cout << "start_stop_pairs_vec: " << ii << ": "
                   << start_stop_pairs_vec[ii].first << " "
                   << start_stop_pairs_vec[ii].second << std::endl;
      }
   }

   auto tp_2c = std::chrono::high_resolution_clock::now();

   std::vector<std::vector<simple_restraint> > nbc_restraints;
   nbc_restraints.resize(n_threads);
   for (std::size_t i=0; i<n_threads; i++)
      nbc_restraints[i].reserve(n_per_thread*80);

   auto d2a = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2a - tp_2 ).count();
   auto d2b = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2b - tp_2a).count();
   auto d2c = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2c - tp_2b).count();
   if (console_output_for_restraints_generation_timings)
      std::cout << "------------------ timings: make_non_bonded_contact_restraints_ng(): find_the_contacts(): "
                << d2a << " start-stop-reserve: " << d2b << " start-stop-push: " << d2c << " milliseconds\n";

   std::atomic<unsigned int> done_count(0);
   bool use_extended_atom_mode = ! model_has_hydrogen_atoms;

   auto tp_3 = std::chrono::high_resolution_clock::now();
   for (std::size_t i=0; i<n_threads; i++) {
      thread_pool_p->push(make_non_bonded_contact_restraints_workpackage_ng,
                          imol,
                          std::cref(geom),
                          std::cref(bonded_atom_indices),
                          std::cref(raic),
                          std::cref(vcontacts),
                          start_stop_pairs_vec[i],
                          std::cref(fixed_atom_indices),
                          std::cref(energy_type_for_atom),
                          use_extended_atom_mode,
                          atom,
                          std::cref(atom_is_metal),
                          std::cref(atom_is_hydrogen),
                          std::cref(H_atom_parent_atom_is_donor_vec),
                          std::cref(atom_is_acceptor_vec),
                          &non_bonded_contacts_atom_indices,
                          &(nbc_restraints[i]),
                          std::ref(done_count));
   }

   auto tp_4 = std::chrono::high_resolution_clock::now();

   // wait for the thread pool to empty - not the right way.
   while (done_count != n_threads) {
      std::this_thread::sleep_for(std::chrono::microseconds(10000));
   }
   auto tp_5 = std::chrono::high_resolution_clock::now();

   // OK now we are back to this thread, add the nbc_restraints to restraints_vec
   //
   int n_extra_restraints = 0;
   for (std::size_t i=0; i<n_threads; i++)
      n_extra_restraints += nbc_restraints[i].size();
   restraints_vec.reserve(restraints_vec.size() + n_extra_restraints);

   // do I want to use std::move here?
   // for (std::size_t i=0; i<n_threads; i++)
   // // restraints_vec.insert(restraints_vec.end(), nbc_restraints[i].begin(), nbc_restraints[i].end());
   // std::move(nbc_restraints[i].begin(), nbc_restraints[i].end(), std::back_inserter(restraints_vec));

   for (std::size_t i=0; i<n_threads; i++) {
      restraints_vec.insert(restraints_vec.end(),
                            std::make_move_iterator(nbc_restraints[i].begin()),
                            std::make_move_iterator(nbc_restraints[i].end()));
      nbc_restraints[i].erase(nbc_restraints[i].begin(), nbc_restraints[i].end());
   }

   auto tp_6 = std::chrono::high_resolution_clock::now();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
   auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
   auto d65 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_6 - tp_5).count();

   if (console_output_for_restraints_generation_timings)
      std::cout << "------------------ timings: non-bonded contacts "  << d32
                << " dispatching threads: " << d43 << " waiting: " << d54
                << " adding NBCs to restraints " << d65 << std::endl;

}


unsigned int
coot::restraints_container_t::make_non_bonded_contact_restraints_ng(int imol,
                                                                    const coot::protein_geometry &geom) {

   unsigned int n_nbc_restraints = 0;

   // std::cout << "make_non_bonded_contact_restraints_ng() " << size() << " "  << std::endl;

   // potentially multithreadable.

   // Use raic to know what is 1-4 related

   // relies on non_bonded_contacts_atom_indices begin allocated to the correct size (n_atoms).

   // raic can be a member of the class
   // to make it: reduced_angle_info_container_t raic(restraints_vec);

   bool extended_atom_mode = false; // turn this on if the model has no Hydrogen atoms;
        extended_atom_mode = ! model_has_hydrogen_atoms;

   float dist_max = 8.0;

   // cache the energy types:
   // maybe this could be a class variable? How long does it take
   // to fill it? ~1ms for rnase.
   //
   auto tp_0 = std::chrono::high_resolution_clock::now();
   std::vector<std::string> energy_type_for_atom(n_atoms);
   std::vector<bool> H_atom_parent_atom_is_donor_vec(n_atoms, false);
   std::vector<bool> atom_is_acceptor_vec(n_atoms, false);

   // needs timing test - might be slow (it isn't)
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      std::string et = get_type_energy(imol, at, geom);
      energy_type_for_atom[i] = et;
      if (H_parent_atom_is_donor(at))
         H_atom_parent_atom_is_donor_vec[i] = true;
      if (is_acceptor(et, geom))
         atom_is_acceptor_vec[i] = true;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
   // std::cout << "   info:: make_non_bonded_contact_restraints_ng(): energy types " << d10 << " microseconds\n";

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map_cache;

   std::set<unsigned int> fixed_atom_flags_set; // signed to unsigned conversion - bleugh.
   std::set<int>::const_iterator it;
   for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); ++it)
      fixed_atom_flags_set.insert(*it);

   // I think that contacts_by_bricks should take a set of ints
   contacts_by_bricks cb(atom, n_atoms, fixed_atom_flags_set);
   cb.set_dist_max(dist_max);
   std::vector<std::set<unsigned int> > vcontacts;
   cb.find_the_contacts(&vcontacts);

   for (std::size_t i=0; i<vcontacts.size(); i++) {
      const std::set<unsigned int> &n_set = vcontacts[i];
      mmdb::Atom *at_1 = atom[i];
      std::set<unsigned int>::const_iterator it_inner;
      for (it_inner=n_set.begin(); it_inner!=n_set.end(); ++it_inner) {

         const unsigned int &j = *it_inner;

         if (bonded_atom_indices[i].find(j) != bonded_atom_indices[i].end())
            continue;

         // for updating non-bonded contacts
         if (non_bonded_contacts_atom_indices[i].find(j) != non_bonded_contacts_atom_indices[i].end())
            continue;

         mmdb::Atom *at_2 = atom[j];

         if (fixed_atom_indices.find(i) != fixed_atom_indices.end())
            if (fixed_atom_indices.find(j) != fixed_atom_indices.end())
               continue;

         if (j < i) /* only add NBC one way round */
            continue;

         std::string res_name_1 = at_1->GetResName();
         std::string res_name_2 = at_2->GetResName();
         int res_no_1 = at_1->GetSeqNum();
         int res_no_2 = at_2->GetSeqNum();

         const std::string &type_1 = energy_type_for_atom[i];
         const std::string &type_2 = energy_type_for_atom[j];

         std::vector<bool> fixed_atom_flags = make_fixed_flags(i, j);

         double dist_min = 3.4;

         bool in_same_residue_flag = (at_1->residue == at_2->residue);
         bool in_same_ring_flag = true;

         // part of this test is not needed.
         if (at_2->residue != at_1->residue) {
            in_same_ring_flag    = false;
            in_same_residue_flag = false;
         }

         std::string atom_name_1 = at_1->GetAtomName();
         std::string atom_name_2 = at_2->GetAtomName();
         std::string element_1   = at_1->element;
         std::string element_2   = at_2->element;

         if (in_same_ring_flag) {

            // in_same_ring_flag = restraints_map[at_2->residue].second.in_same_ring(atom_name_1,
            //                                                                       atom_name_2);

            in_same_ring_flag = is_in_same_ring(imol, at_2->residue,
                                                residue_ring_map_cache,
                                                atom_name_1, atom_name_2, geom);
         }

         // this doesn't check 1-4 over a moving->non-moving peptide link (see comment above function)
         // because the non-moving atom doesn't have angle restraints.
         //
         bool is_1_4_related = raic.is_1_4(i, j, fixed_atom_flags);

         if (false)
            std::cout << "here with at_1 " << atom_spec_t(at_1) << " at_2 " << atom_spec_t(at_2)
                      << " is_1_4_related " << is_1_4_related << std::endl;

         bool mc_atoms_tandem = false;
         if (! is_1_4_related) {
            // hacky case for C in a helix.
            // Also N.

            if (atom_name_1 == " C  ")
               if (atom_name_2 == " C  ")
                  if (at_2->residue->index - at_1->residue->index == 1) {
                     if (false)
                        std::cout << "DEBUG:: in make_non_bonded_contact_restraints_ng() C to C neighbs "
                                  << at_1->residue->index << " " << at_2->residue->index
                                  << std::endl;
                     mc_atoms_tandem = true;
                  }
            if (atom_name_1 == " N  ")
               if (atom_name_2 == " N  ")
                  if (at_2->residue->index - at_1->residue->index == -1) {
                     mc_atoms_tandem = true;
                     std::cout << "-------- Here 2 " << at_1->residue->index << " " << at_2->residue->index
                               << std::endl;
                  }
         }

         if (mc_atoms_tandem)
            is_1_4_related = true;

         if (is_1_4_related) {
            dist_min = 2.64; // was 2.7 but c.f. guanine ring distances

            if (mc_atoms_tandem) dist_min = 2.99;

            if (is_hydrogen(at_1))
               dist_min -= 0.7;
            if (is_hydrogen(at_2))
               dist_min -= 0.7;
         } else {

            std::pair<bool, double> nbc_dist = geom.get_nbc_dist_v2(type_1, type_2,
                                                                    element_1, element_2,
                                                                    atom_is_metal[i],
                                                                    atom_is_metal[j],
                                                                    extended_atom_mode,
                                                                    in_same_residue_flag,
                                                                    in_same_ring_flag);

            if (nbc_dist.first) {

               // In a helix O(n) is close to C(n+1), we should allow it.
               // 
               bool is_O_C_1_5_related = check_for_O_C_1_5_relation(at_1, at_2);

               if (is_O_C_1_5_related) {
                  dist_min = 2.84;
               } else {

                  // Perhaps we don't have angle restraints to both atoms because one
                  // of the atoms is fixed (and thus miss that these have a 1-4 relationship).
                  // e.g. O(n) [moving] -> CA(n+1) [fixed]
                  // 
                  // (this test will fail on insertion codes)
                  //

                  bool strange_exception = false;
                  int rn_diff = abs(res_no_2 - res_no_1);
                  if (rn_diff == 1) {
                     if (fixed_atom_flags.size()) {
                        if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
                           if (atom_name_1 == " O  ")
                              if (atom_name_2 == " CA ") 
                                 strange_exception = true;
                           if (atom_name_1 == " CA ")
                              if (atom_name_2 == " O  ")
                                 strange_exception = true;
                           if (atom_name_1 == " N  ")
                              if (atom_name_2 == " CB ")
                                 strange_exception = true;
                           if (atom_name_1 == " CB ")
                              if (atom_name_2 == " N  ")
                                 strange_exception = true;
                           if (atom_name_1 == " C  ")
                              if (atom_name_2 == " CB ")
                                 strange_exception = true;
                        }
                     }
                     if (strange_exception)
                        dist_min = 2.7;

                     // Strange that these are not marked as 1-4 related.  Fix here...
                     // HA-CA-N-C can be down to ~2.4A.
                     // HA-CA-C-N can be down to ~2.41A.
                     if (res_no_2 > res_no_1) {
                        if (atom_name_1 == " C  ") {
                           if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
                              strange_exception = true;
                              dist_min = 2.4;
                           }
                        }
                        if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
                           if (atom_name_2 == " N  ") {
                              strange_exception = true;
                              dist_min = 2.41;
                           }
                        }
                        if (atom_name_1 == " N  ") {
                           if (atom_name_2 == " H  ") {
                              strange_exception = true;
                              dist_min = 2.4;
                           }
                        }
                     } else {
                        if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
                           if (atom_name_2 == " C  ") {
                              strange_exception = true;
                              dist_min = 2.4;
                           }
                        }
                        if (atom_name_1 == " N  ") {
                           if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
                              strange_exception = true;
                              dist_min = 2.41;
                           }
                        }
                        if (atom_name_2 == " N  ") {
                           if (atom_name_1 == " H  ") {
                              strange_exception = true;
                              dist_min = 2.4;
                           }
                        }
                     }
                  }
                  if (rn_diff == 2) {
                     if (fixed_atom_flags.size()) {
                        if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
                           if (atom_name_1 == " C  ")
                              if (atom_name_2 == " N  ")
                                 strange_exception = true;
                           if (atom_name_1 == " N  ")
                              if (atom_name_2 == " C  ")
                                 strange_exception = true; // 3.1 would be enough

                           if (strange_exception)
                              dist_min = 2.7;
                        }
                     }
                  }

                  if (! strange_exception)
                     dist_min = nbc_dist.second;
               }
            } else {
               // short/standard value for no nbc value calculated
               // std::cout << ".......... using short/standard 2.8" << std::endl;
               dist_min = 2.8;
            }
         }

         bool is_H_non_bonded_contact = false;

         if (is_hydrogen(at_1)) {
            is_H_non_bonded_contact = true;
            if (H_parent_atom_is_donor(at_1))
               if (is_acceptor(type_2, geom))
                  dist_min -= 0.7;
         }
         if (is_hydrogen(at_2)) {
            is_H_non_bonded_contact = true;
           if (H_parent_atom_is_donor(at_2))
               if (is_acceptor(type_1, geom))
                  dist_min -= 0.7;
         }

         non_bonded_contacts_atom_indices[i].insert(j);
         simple_restraint::nbc_function_t nbcf = simple_restraint::LENNARD_JONES;
         simple_restraint r(NON_BONDED_CONTACT_RESTRAINT,
                            nbcf, i, j,
                            is_H_non_bonded_contact,
                            fixed_atom_flags, dist_min);
         if (false) // debug
            std::cout << "Adding NBC " << i << " " << *it << " " << energy_type_for_atom[i] << " "
                              << energy_type_for_atom[j] << " "
                              << is_H_non_bonded_contact << " "
                              << fixed_atom_flags[0] << " " << fixed_atom_flags[1] << " "
                              << dist_min <<  "\n";

         n_nbc_restraints++;
         r.n_atoms_from_all_restraints = n_atoms; // for debugging crash in non-bonded contact
                                                  // restraints
         r.restraints_index = size(); // likewise
         restraints_vec.push_back(r); // use push_back_restraint

      }
   }

   make_df_restraints_indices();
   // is the following needed? I doubt.
   make_distortion_electron_density_ranges();

   return n_nbc_restraints;

}

coot::restraints_container_t::link_restraints_counts
coot::restraints_container_t::make_link_restraints_for_link_ng(const coot::new_linked_residue_t &nlr,
                                                               const coot::protein_geometry &geom) {

   bool do_trans_peptide_restraints = false; // we are not linking peptides with this function

   // order switch is done (if needed) in the constructor of new_linked_residue_t

   return make_link_restraints_for_link_ng(nlr.link_type,
                                           nlr.res_1,
                                           nlr.res_2,
                                           nlr.is_fixed_first,
                                           nlr.is_fixed_second,
                                           do_trans_peptide_restraints,
                                           geom);

}


coot::restraints_container_t::link_restraints_counts
coot::restraints_container_t::make_link_restraints_for_link_ng(const std::string &link_type,
                                                               mmdb::Residue *res_1,
                                                               mmdb::Residue *res_2,
                                                               bool is_fixed_first_residue,
                                                               bool is_fixed_second_residue,
                                                               bool do_trans_peptide_restraints,
                                                               const coot::protein_geometry &geom) {

   link_restraints_counts lrc;

   // restraints_usage_flag = BONDS_ANGLES_AND_CHIRALS;

   if (false)
      std::cout << "make_link_restraints_for_link_ng(): "
                << " type: " << link_type << " "
                << residue_spec_t(res_1) << " " << residue_spec_t(res_2) << " "
                << is_fixed_first_residue << " " << is_fixed_second_residue << " "
                << "restraints_usage_flag: " << restraints_usage_flag << std::endl;

   if (restraints_usage_flag & BONDS_MASK)
      lrc.n_link_bond_restr += add_link_bond(link_type,
                                             res_1, res_2,
                                             is_fixed_first_residue,
                                             is_fixed_second_residue,
                                             geom);

   if (restraints_usage_flag & ANGLES_MASK)
      lrc.n_link_angle_restr += add_link_angle(link_type,
                                               res_1, res_2,
                                               is_fixed_first_residue,
                                               is_fixed_second_residue,
                                               geom);

   if (restraints_usage_flag & TORSIONS_MASK)
      lrc.n_link_torsion_restr += add_link_torsion(link_type, res_1, res_2,
                                                   is_fixed_first_residue, is_fixed_second_residue, geom);


   if (restraints_usage_flag & TRANS_PEPTIDE_MASK) {
      if (do_trans_peptide_restraints) {
         lrc.n_link_trans_peptide += add_link_trans_peptide(res_1, res_2,
                                                            is_fixed_first_residue,
                                                            is_fixed_second_residue, false); // don't add if cis
      } else {
         if (false) // debug (we don't want to try to add trans-pep restraints for SS bonds (for example))
            std::cout << "make_link_restraints_for_link_ng(): trans-pep flag off "
                      << residue_spec_t(res_1) << " " << residue_spec_t(res_2) << " "
                      << restraints_usage_flag << std::endl;
      }
   }

   if (restraints_usage_flag & PLANES_MASK)
      lrc.n_link_plane_restr += add_link_plane(link_type,
                                               res_1, res_2,
                                               is_fixed_first_residue,
                                               is_fixed_second_residue,
                                               geom);

   return lrc;
}

std::string
coot::restraints_container_t::find_peptide_link_type_ng(mmdb::Residue *res_1,
                                                        mmdb::Residue *res_2,
                                                        const coot::protein_geometry &geom) const {

   std::string link_type;

   std::string t1;
   std::string t2;

   std::string residue_type_1 = res_1->name;
   std::string residue_type_2 = res_2->name;

   for (unsigned int idr=0; idr<geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_1)) {
         t1 = geom[idr].second.residue_info.group;
         break;
      }
   }
   for (unsigned int idr=0; idr<geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_2)) {
         t2 = geom[idr].second.residue_info.group;
         break;
      }
   }

   // This should be true because we do a test on amino acid type at the start of try_make_peptide_link_ng
   //
   if (t1 == "L-peptide" || t1 == "D-peptide" || t1 == "M-peptide" || t1 == "P-peptide" || t1 == "peptide") {
      if (t2 == "L-peptide" || t2 == "D-peptide" || t2 == "M-peptide" || t2 == "P-peptide" || t2 == "peptide") {
         if (residue_type_2 == "PRO" || residue_type_2 == "HYP") {
            link_type = "PTRANS";
         } else {
            link_type = "TRANS";
         }
      }
   }

   // CIS and PCIS will need a torsion check around
   // r1::CA - r1::C - r2::N - r2::CA

   mmdb::Atom **residue_1_atoms = 0;
   mmdb::Atom **residue_2_atoms = 0;
   int n_residue_1_atoms;
   int n_residue_2_atoms;
   mmdb::Atom *ca_r_1 = 0;
   mmdb::Atom *c_r_1  = 0;
   mmdb::Atom *n_r_2  = 0;
   mmdb::Atom *ca_r_2 = 0;
   res_1->GetAtomTable(residue_1_atoms, n_residue_1_atoms);
   res_2->GetAtomTable(residue_2_atoms, n_residue_2_atoms);
   for (int iat=0; iat<n_residue_1_atoms; iat++) {
      mmdb::Atom *at_1 = residue_1_atoms[iat];
      std::string at_name_1(at_1->GetAtomName());
      if (at_name_1 == " C  ")// PDBv3 FIXE
         c_r_1 = at_1;
      if (at_name_1 == " CA ") // PDBv3 FIXE
         ca_r_1 = at_1;
   }
   for (int iat=0; iat<n_residue_2_atoms; iat++) {
      mmdb::Atom *at_2 = residue_2_atoms[iat];
      std::string at_name_2(at_2->GetAtomName());
      if (at_name_2 == " N  ")// PDBv3 FIXE
         n_r_2 = at_2;
      if (at_name_2 == " CA ")// PDBv3 FIXE
         ca_r_2 = at_2;
   }
   if (ca_r_1 && c_r_1 && n_r_2 && ca_r_2) {
      clipper::Coord_orth pt_1 = co(ca_r_1);
      clipper::Coord_orth pt_2 = co(c_r_1);
      clipper::Coord_orth pt_3 = co(n_r_2);
      clipper::Coord_orth pt_4 = co(ca_r_2);
      double torsion = clipper::Coord_orth::torsion(pt_1, pt_2, pt_3, pt_4);
      if (torsion > M_PI)
         torsion -= 2.0 * M_PI;
      // std::cout << "torsion " << clipper::Util::rad2d(torsion) << std::endl;
      double ninety_degrees = 0.5 * M_PI;
      if (torsion < ninety_degrees && torsion > -ninety_degrees) {
         if (link_type ==  "TRANS") link_type =  "CIS";
         if (link_type == "PTRANS") link_type = "PCIS";
      }
   }
   return link_type;
}



std::pair<bool, coot::restraints_container_t::link_restraints_counts>
coot::restraints_container_t::try_make_peptide_link_ng(const coot::protein_geometry &geom,
                                                       std::pair<bool, mmdb::Residue *> res_1_pair,
                                                       std::pair<bool, mmdb::Residue *> res_2_pair,
                                                       bool do_trans_peptide_restraints) {

   if (false)
      std::cout << "try_make_peptide_link_ng():         "
                << residue_spec_t(res_1_pair.second) << " " << residue_spec_t(res_2_pair.second)
                << " fixed-1: " << res_1_pair.first << " fixed-2: " << res_2_pair.first
                << " trans " << do_trans_peptide_restraints
                << std::endl;

   mmdb::Residue *res_1 = res_1_pair.second;
   mmdb::Residue *res_2 = res_2_pair.second;
   std::string res_name_1(res_1->GetResName());
   std::string res_name_2(res_2->GetResName());
   bool status = false;
   link_restraints_counts lrc;

   // find_peptide_link_type_ng doesn't test the geometry -
   // just the residues types, return empty on the groups being peptide linked
   //
   std::string link_type = find_peptide_link_type_ng(res_1, res_2, geom);
   if (! link_type.empty()) {
      {
         mmdb::Atom **residue_1_atoms = 0;
         mmdb::Atom **residue_2_atoms = 0;
         int n_residue_1_atoms;
         int n_residue_2_atoms;
         res_1->GetAtomTable(residue_1_atoms, n_residue_1_atoms);
         res_2->GetAtomTable(residue_2_atoms, n_residue_2_atoms);
         for (int iat_1=0; iat_1<n_residue_1_atoms; iat_1++) {
            mmdb::Atom *at_1 = residue_1_atoms[iat_1];
            std::string at_name_1(at_1->GetAtomName());
            if (at_name_1 == " C  ") { // PDBv3 FIXE
               std::string alt_conf_1(at_1->altLoc);
               for (int iat_2=0; iat_2<n_residue_2_atoms; iat_2++) {
                  mmdb::Atom *at_2 = residue_2_atoms[iat_2];
                  std::string at_name_2(at_2->GetAtomName());
                  if (at_name_2 == " N  ") { // PDBv3 FIXE
                     std::string alt_conf_2(at_2->altLoc);
                     if (alt_conf_1 == alt_conf_2 || alt_conf_1.empty() || alt_conf_2.empty()) {

                        bool is_fixed_first_residue  = res_1_pair.first;
                        bool is_fixed_second_residue = res_2_pair.first;

                        if (false)
                           std::cout << "###################### adding link bond: "
                                     << atom_spec_t(at_1) << " " << atom_spec_t(at_2) << std::endl;

                        lrc = make_link_restraints_for_link_ng(link_type, res_1, res_2,
                                                               is_fixed_first_residue,
                                                               is_fixed_second_residue,
                                                               do_trans_peptide_restraints,
                                                               geom);
                        status = true;
                     }
                  }
               }
            }
         }
      }
   }

   return std::pair<bool, coot::restraints_container_t::link_restraints_counts>(status, lrc);

}

#include "coot-utils/coot-coord-extras.hh"

std::pair<bool, coot::restraints_container_t::link_restraints_counts>
coot::restraints_container_t::try_make_phosphodiester_link_ng(const coot::protein_geometry &geom,
                                                              std::pair<bool, mmdb::Residue *> res_1_pair,
                                                              std::pair<bool, mmdb::Residue *> res_2_pair) {

   // if this ever happens in real life, this can be enabled:
   //
   bool use_distance_cut_off = false; // distance test for residues that are not tandem in sequence
                                      // but are tandem in index.
   const float distance_cut_off = 5.0; // residues with O3' and P more thant this that are
                                 // not tandem in sequence will not be linked
   const float distance_cut_off_srd = distance_cut_off * distance_cut_off;

   mmdb::Residue *res_1 = res_1_pair.second;
   mmdb::Residue *res_2 = res_2_pair.second;
   std::string res_name_1(res_1->GetResName());
   std::string res_name_2(res_2->GetResName());
   bool status = false;
   link_restraints_counts lrc;

   if (util::is_nucleotide_by_dict(res_1, geom)) {
      if (util::is_nucleotide_by_dict(res_2, geom)) {

         mmdb::Atom **residue_1_atoms = 0;
         mmdb::Atom **residue_2_atoms = 0;
         int n_residue_1_atoms;
         int n_residue_2_atoms;
         res_1->GetAtomTable(residue_1_atoms, n_residue_1_atoms);
         res_2->GetAtomTable(residue_2_atoms, n_residue_2_atoms);
         for (int iat_1=0; iat_1<n_residue_1_atoms; iat_1++) {
            mmdb::Atom *at_1 = residue_1_atoms[iat_1];
            std::string at_name_1(at_1->GetAtomName());
            if (at_name_1 == " O3'") { // PDBv3 FIXME
               std::string alt_conf_1(at_1->altLoc);
               for (int iat_2=0; iat_2<n_residue_2_atoms; iat_2++) {
                  mmdb::Atom *at_2 = residue_2_atoms[iat_2];
                  std::string at_name_2(at_2->GetAtomName());
                  if (at_name_2 == " P  ") { // PDBv3 FIXE
                     std::string alt_conf_2(at_2->altLoc);
                     if (alt_conf_1 == alt_conf_2 || alt_conf_1.empty() || alt_conf_2.empty()) {

                        if (use_distance_cut_off) {
                           int res_no_1 = res_1->GetSeqNum();
                           int res_no_2 = res_2->GetSeqNum();
                           if ((res_no_2-res_no_1) > 1) {
                              clipper::Coord_orth pt_1 = coot::co(at_1);
                              clipper::Coord_orth pt_2 = coot::co(at_2);
                              if ((pt_2-pt_1).lengthsq() > distance_cut_off_srd)
                                 continue;
                           }
                        }

                        // find_peptide_link_type_ng doesn't test the geometry -
                        // just the residues types
                        std::string link_type = "p";
                        bool is_fixed_first_residue  = res_1_pair.first;
                        bool is_fixed_second_residue = res_2_pair.first;
                        lrc = make_link_restraints_for_link_ng(link_type, res_1, res_2,
                                                               is_fixed_first_residue,
                                                               is_fixed_second_residue,
                                                               false, geom);
                        status = true;
                     }
                  }
               }
            }
         }
      }
   }

   return std::pair<bool, coot::restraints_container_t::link_restraints_counts>(status, lrc);

}


void
coot::restraints_container_t::disperse_plane_restraints() {

#if 0 // I don't think this function is necessary, because the plane restraints get
      // dispersed as they are allocated to threads

   if (n_threads > 1) {
      std::pair<unsigned int, unsigned int> restraints_limit_for_planes(999999999, 0);
      bool found = false;
      for (unsigned int i=0; i<restraints_vec.size(); i++) {
         const simple_restraint &restraint = restraints_vec[i];
         if (restraint.restraint_type == PLANE_RESTRAINT) {
            found = true;
            if (i < restraints_limit_for_planes.first)
               restraints_limit_for_planes.first = i;
            if (i > restraints_limit_for_planes.second)
               restraints_limit_for_planes.second = i;
         }
      }

      if (found) {
         for (unsigned int idx=restraints_limit_for_planes.first; idx<=restraints_limit_for_planes.second; idx++) {
            const simple_restraint &restraint = restraints_vec[idx];
            if (restraint.restraint_type == PLANE_RESTRAINT) {
               unsigned int idx_base = idx - restraints_limit_for_planes.first;
               unsigned int idx_from = idx_base * n_threads;
               if (idx != idx_from) {
                  if (idx_from < restraints_vec.size()) {
                     std::cout << "swapping " << idx_from << " " << idx << std::endl;
                     std::swap(restraints_vec[idx_from], restraints_vec[idx]);
                  }
               }
            }
         }
      }
   }
#endif
}

bool
coot::restraints_container_t::N_and_C_are_close_ng(mmdb::Residue *res_1,
                                                   mmdb::Residue *res_2,
                                                   float d_crit) const {

   bool status = false;

   mmdb::Atom *at_1 = res_1->GetAtom(" C  ");
   mmdb::Atom *at_2 = res_2->GetAtom(" N  ");
   if (at_1) {
      if (at_2) {
         float dd =
            (at_1->x - at_2->x) * (at_1->x - at_2->x) +
            (at_1->y - at_2->y) * (at_1->y - at_2->y) +
            (at_1->z - at_2->z) * (at_1->z - at_2->z);
         if (dd < d_crit * d_crit)
            status = true;
      }
   }
   /// No linking if there was a TER in res_1
   if (status) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      res_1->GetAtomTable(residue_atoms, n_residue_atoms);
      for(int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (at->isTer()) {
            status = false;
            break;
         }
      }
   }
   return status;
}

bool
coot::restraints_container_t::O3prime_and_P_are_close_ng(mmdb::Residue *res_1,
                                                         mmdb::Residue *res_2,
                                                         float d_crit) const {

   bool status = false;

   mmdb::Atom *at_1 = res_1->GetAtom(" O3'"); // PDBv3 FIXME
   mmdb::Atom *at_2 = res_2->GetAtom(" P  ");
   if (at_1) {
      if (at_2) {
         float dd =
            (at_1->x - at_2->x) * (at_1->x - at_2->x) +
            (at_1->y - at_2->y) * (at_1->y - at_2->y) +
            (at_1->z - at_2->z) * (at_1->z - at_2->z);
         if (dd < d_crit * d_crit)
            status = true;
      }
   }
   /// No linking if there was a TER in res_1
   if (status) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      res_1->GetAtomTable(residue_atoms, n_residue_atoms);
      for(int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (at->isTer()) {
            status = false;
            break;
         }
      }
   }
   return status;
}


void
coot::restraints_container_t::make_polymer_links_ng(const coot::protein_geometry &geom,
                                                    bool do_rama_plot_restraints,
                                                    bool do_trans_peptide_restraints,
                                                    std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > *residue_link_vector_map_p,
                                                    std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > *residue_pair_link_set_p) {

   // make_polymer_links_ng uses the residues_vec (which is sorted in initialization)
   // so, only the passed active atoms are considered for polymer links

   // std::cout << "INFO:: making link polymer links " << std::endl;

   std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > &residue_pair_link_set = *residue_pair_link_set_p;

   int idx_reference_index_handle = mol->GetUDDHandle(mmdb::UDR_RESIDUE, "index from reference residue");

   link_restraints_counts accum_links;

   if (residues_vec.size() == 1)
      return; // nothing to link

   // the residues are sorted by now, we can't have a forward link on the last residue of residues_vec
   //
   int n=residues_vec.size()-1;
   for (int i1=0; i1<n; i1++) {
      int i2=i1+1;

      if (residues_vec[i1].first && residues_vec[i2].first)
         continue;

      mmdb::Residue *res_1 = residues_vec[i1].second;
      mmdb::Residue *res_2 = residues_vec[i2].second;

      std::string res_name_1(res_1->GetResName());
      std::string res_name_2(res_2->GetResName());
      if (res_name_1 == "HOH") continue; // break?
      if (res_name_2 == "HOH") continue;

      if (res_1->chain == res_2->chain) {
         int serial_delta = res_2->index - res_1->index;

         // *this* is the serial_delta we should be checking
         int ref_index_1 = residues_vec[i1].second->index;
         int ref_index_2 = residues_vec[i2].second->index;

         serial_delta = ref_index_2 - ref_index_1;

         // but I don't think that this correctly deals with
         // antibodies that have missing residue numbers - but
         // are linked with a peptide: 19-20-22-23

         if (serial_delta == 1) {
            // possibly polymer

            int res_no_delta = res_2->GetSeqNum() - res_1->GetSeqNum();

            std::pair<bool, link_restraints_counts> results;
            // first test for normal peptide residue
            //
            if (res_no_delta == 1) {

               // this can fail if there are missing link (C,N) atoms
               results = try_make_peptide_link_ng(geom, residues_vec[i1], residues_vec[i2],
                                                  do_trans_peptide_restraints);
            } else {
               if (res_no_delta == 0) { // insertion code?
                  std::string ins_code_1(res_1->GetInsCode());
                  std::string ins_code_2(res_2->GetInsCode());
                  if (ins_code_1 != ins_code_2) {
                     results = try_make_peptide_link_ng(geom, residues_vec[i1], residues_vec[i2],
                                                        do_trans_peptide_restraints);
                  }
               } else {

                  // antibody case: H212 is peptide linked to H216 5d6c (no LINK record)
                  if (N_and_C_are_close_ng(res_1, res_2, 2.0))
                     results = try_make_peptide_link_ng(geom, residues_vec[i1], residues_vec[i2],
                                                        do_trans_peptide_restraints);

               }
            }

            if (results.first) {
               accum_links.add(results.second);
            } else {

               if (res_no_delta < 2) {
                  results = try_make_phosphodiester_link_ng(geom, residues_vec[i1], residues_vec[i2]);
                  if (results.first) {
                     accum_links.add(results.second);
                  }
               } else {
                  if (O3prime_and_P_are_close_ng(res_1, res_2, 2.0)) {
                     results = try_make_phosphodiester_link_ng(geom, residues_vec[i1], residues_vec[i2]);
                     if (results.first)
                        accum_links.add(results.second);
                  }
               }
            }

            if (results.first) {

               // make a note that those residues were linked

               (*residue_link_vector_map_p)[res_1].push_back(res_2);
               (*residue_link_vector_map_p)[res_2].push_back(res_1);

               std::pair<mmdb::Residue *, mmdb::Residue *> residue_pair_f(res_1, res_2);
               std::pair<mmdb::Residue *, mmdb::Residue *> residue_pair_b(res_2, res_1);
               residue_pair_link_set.insert(residue_pair_f);
               residue_pair_link_set.insert(residue_pair_b);
            }
         }
      }
   }

   if (verbose_geometry_reporting != QUIET)
      accum_links.report();

}


void
coot::restraints_container_t::add_header_metal_link_bond_ng(const coot::atom_spec_t &atom_spec_1,
                                                            const coot::atom_spec_t &atom_spec_2,
                                                            double  dist) {

   // we know the specs from the link line - we need to find the atom indices
   //
   // If speed is needed here, we could group up all the atom metal bonds and call this
   // function as a vector - then we could cache the atom indices, rather than forget
   // them every time this function is called.
   // At the moment, if it's not slow, don't speed it up.

   int index_1 = -1;
   int index_2 = -1;

   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      int rn_at = at->GetSeqNum();
      if (rn_at == atom_spec_1.res_no) {
         atom_spec_t s(at);
         if (s == atom_spec_1) {
            index_1 = i;
            continue;
         }
      }
      if (rn_at == atom_spec_2.res_no) {
         atom_spec_t s(at);
         if (s == atom_spec_2) {
            index_2 = i;
         }
      }
      if ((index_1 != -1) && (index_2 != -1))
         break;
   }

   if (false)
      std::cout << "debug:: in add_header_metal_link_bond_ng() indices: " << index_1 << " " << index_2
                << std::endl;

   if (index_1 != -1) {
      if (index_2 != -1) {
         bonded_atom_indices[index_1].insert(index_2);
         bonded_atom_indices[index_2].insert(index_1);
         std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);
         add(BOND_RESTRAINT, index_1, index_2, fixed_flags, dist, 0.1, 1.2);
      }
   }

}

void
coot::restraints_container_t::make_header_metal_links_ng(const coot::protein_geometry &geom) {
   // we need to use the links passed to the constructor - currently they are not saved.


   int imol = protein_geometry::IMOL_ENC_ANY;
   for (std::size_t i=0; i<links.size(); i++) {
     const mmdb::Link &link = links[i];
      atom_spec_t a1(link.chainID1, link.seqNum1, link.insCode1, link.atName1, link.aloc1);
      atom_spec_t a2(link.chainID2, link.seqNum2, link.insCode2, link.atName2, link.aloc2);
      if ((a1.alt_conf == a2.alt_conf) || a1.alt_conf.empty() || a2.alt_conf.empty()) {
         std::string rn_1 = link.resName1;
         std::string rn_2 = link.resName2;

         bool is_oxygen_a1   = false;
         bool is_oxygen_a2   = false;
         bool is_nitrogen_a1 = false;
         bool is_nitrogen_a2 = false;
         bool is_sulfur_a1   = false;
         bool is_sulfur_a2   = false;

         std::pair<bool, dict_atom> da_1 = geom.get_monomer_atom_info(rn_1, a1.atom_name, imol);
         std::pair<bool, dict_atom> da_2 = geom.get_monomer_atom_info(rn_2, a2.atom_name, imol);

         if (da_1.first) {
            if (da_1.second.type_symbol == "O") is_oxygen_a1   = true;
            if (da_1.second.type_symbol == "S") is_sulfur_a1   = true;
            if (da_1.second.type_symbol == "N") is_nitrogen_a1 = true;
         }
         if (da_2.first) {
            if (da_2.second.type_symbol == "O") is_oxygen_a2   = true;
            if (da_2.second.type_symbol == "S") is_sulfur_a2   = true;
            if (da_2.second.type_symbol == "N") is_nitrogen_a2 = true;
         }

         if (false)
            std::cout << "here with link "
                      << a1 << " is_oxygen_a1 " << is_oxygen_a1 << " "
                      << a2 << " is_oxygen_a2 " << is_oxygen_a2 << std::endl;

         if (is_oxygen_a1) {
            std::map<std::string, double>::const_iterator it = geom.metal_O_map.find(rn_2);
            if (it != geom.metal_O_map.end()) {
               add_header_metal_link_bond_ng(a1, a2, it->second);
            }
         }
         if (is_nitrogen_a1) {
            std::map<std::string, double>::const_iterator it = geom.metal_N_map.find(rn_2);
            if (it != geom.metal_N_map.end()) {
               // std::cout << "Yay! Make metal N bond restraint " << a1 << " to metal " << a2 << std::endl;
               add_header_metal_link_bond_ng(a1, a2, it->second);
            }
         }
         if (is_sulfur_a1) {
            std::map<std::string, double>::const_iterator it = geom.metal_S_map.find(rn_2);
            if (it != geom.metal_S_map.end()) {
               // std::cout << "Yay! Make metal S bond restraint " << a1 << " to metal " << a2 << std::endl;
               add_header_metal_link_bond_ng(a1, a2, it->second);
            }
         }
         if (is_oxygen_a2) {
            std::map<std::string, double>::const_iterator it = geom.metal_O_map.find(rn_1);
            if (it != geom.metal_O_map.end()) {
               add_header_metal_link_bond_ng(a1, a2, it->second);
            }
         }
         if (is_nitrogen_a2) {
            std::map<std::string, double>::const_iterator it = geom.metal_N_map.find(rn_1);
            if (it != geom.metal_N_map.end()) {
               // std::cout << "Yay! Make metal N bond restraint " << a2 << " to metal " << a1 << std::endl;
               add_header_metal_link_bond_ng(a1, a2, it->second);
            }
         }
         if (is_sulfur_a2) {
            std::map<std::string, double>::const_iterator it = geom.metal_S_map.find(rn_1);
            if (it != geom.metal_S_map.end()) {
               // std::cout << "Yay! Make metal S bond restraint " << a2 << " to metal " << a1 << std::endl;
               add_header_metal_link_bond_ng(a1, a2, it->second);
            }
         }
      }
   }
   // std::cout << "------------------------- done make_header_metal_links_ng() " << std::endl;
}


// non polymer links and carbohydrate links.
// And metal links
coot::restraints_container_t::link_restraints_counts
coot::restraints_container_t::make_other_types_of_link(const coot::protein_geometry &geom,
                                                       const std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > &residue_link_vector_map,
                                                       const std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > &residue_pair_link_set) {

   // ------------- Generic/non-polymer -------------------

   // carbohydrate links are generated using LINKs in the header, not here, make_link_restraints_from_links()

   auto is_member = [] (const std::vector<std::string> &names, const std::string &n) {
      return std::find(names.begin(), names.end(), n) != names.end();
   };


   // 21/10/2025-PE I no longer want carbohydrate to be linked by distance (and this is a link-by-distance function)

   link_restraints_counts lrc("CHO-SS-and-other");

   if (false) {
      std::cout << "debug----------- fixed_atom_indices " << std::endl;
      std::set<int>::const_iterator it;
      for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); ++it)
         std::cout << "   fixed atoms: " << *it << std::endl;
   }

   make_header_metal_links_ng(geom);

   new_linked_residue_list_t nlrs;
   std::vector<std::string> cho_residue_types = {"MAN", "BMA", "NAG", "FUC", "GLC"};

   std::set<unsigned int> fixed_atom_flags_set;
   {
      std::set<int>::const_iterator it;
      for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); ++it)
         fixed_atom_flags_set.insert(*it);
   }
   contacts_by_bricks cb(atom, n_atoms, fixed_atom_flags_set);
   cb.set_dist_max(3.0);
   std::vector<std::set<unsigned int> > vcontacts;
   bool only_between_different_residues_flag = true;
   cb.find_the_contacts(&vcontacts, only_between_different_residues_flag);

   std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > tested_but_nothing;

   bool debug_vcontacts = false;
   if (debug_vcontacts) {
      for (std::size_t i=0; i<vcontacts.size(); i++) {
         const std::set<unsigned int> &n_set = vcontacts[i];
         if (! n_set.empty()) {
           mmdb::Atom *at_1 = atom[i];
            std::cout << "vcontact for " << i << " " << atom_spec_t(at_1) << ": ";
            std::set<unsigned int>::const_iterator it_inner;
            for (it_inner=n_set.begin(); it_inner!=n_set.end(); ++it_inner) {
               mmdb::Atom *at_2 = atom[*it_inner];
               std::cout << " " << atom_spec_t(at_2);
            }
            std::cout << std::endl;
         }
      }
   }

   for (std::size_t i=0; i<vcontacts.size(); i++) {
      const std::set<unsigned int> &n_set = vcontacts[i];
      mmdb::Atom *at_1 = atom[i];
      if (strcmp(at_1->element, " H") == 0) continue;
      // if (! is_fully_linked_ng(at_1->residue, residue_link_count_map)) {
      if (true) {
         std::set<unsigned int>::const_iterator it;
         for (it=n_set.begin(); it!=n_set.end(); ++it) {
            mmdb::Atom *at_2 = atom[*it];

            if (strcmp(at_2->element, " H") == 0) continue;

            mmdb::Residue *res_1 = at_1->residue;
            mmdb::Residue *res_2 = at_2->residue;

            std::string res_name_1(res_1->GetResName());
            std::string res_name_2(res_2->GetResName());

            if (res_name_1 == "HOH") continue; // waters don't make links that are not
            if (res_name_2 == "HOH") continue; // handled in make_header_metal_links_ng()

            if (res_1 == res_2) {
               // these should not be here - they should be filtered out in contacts_by_bricks.
               // At some state test if this is still needed.
               continue;
            }

            // now delete carbohydrate links - they will be done by header
            if (is_member(cho_residue_types, res_name_1)) continue;
            if (is_member(cho_residue_types, res_name_2)) continue;

            // Keitaro wanted this test deleted. OK. 20220917-PE
            // if (at_2->GetChain() == at_1->GetChain()) {
            //
            if (true) {

               // no links for non-moving bumping residues
               //
               if (fixed_atom_indices.find(i) != fixed_atom_indices.end())
                  if (fixed_atom_indices.find(*it) != fixed_atom_indices.end())
                     continue;

               // at_2 might have a lower residue number than at_1.

               // At this stage most of the contacts will be peptide contacts (generally speaking)
               // How can I know if the residues of these atoms have already been linked
               // (quite likely they will have been)

               std::pair<mmdb::Residue *, mmdb::Residue *> test_residue_pair(res_1, res_2);

               if (residue_pair_link_set.find(test_residue_pair) == residue_pair_link_set.end())  {
                  // it was a new pair.

                  if (false) {
                     std::cout << "                 Here in make_other_types_of_link() A with at_1 at_2:    "
                               << coot::atom_spec_t(at_1) << " " << coot::atom_spec_t(at_2) << std::endl;
                     std::cout << "                 failed to find these residues in the polyer-linked set: "
                               << residue_spec_t(res_1) << " "
                               << residue_spec_t(res_2) << std::endl;
                  }

                  if (false) {
                     std::cout << "debug:: Here are the pairs in the linked set " << std::endl;
                     std::set<std::pair<mmdb::Residue *, mmdb::Residue *> >::const_iterator it_inner;
                     for (it_inner=residue_pair_link_set.begin(); it_inner!=residue_pair_link_set.end(); ++it_inner)
                        std::cout << "   " << residue_spec_t(it_inner->first) << " " << residue_spec_t(it_inner->second)
                                  << std::endl;
                  }

                  if (nlrs.already_added_p(res_1, res_2)) {
                     // std::cout << "      already added " << std::endl;
                     continue;
                  }

                  // add a check here for res_1 and res_2 already considered, but not added.
                  // continue, in that case.
                  //
                  std::pair<mmdb::Residue *, mmdb::Residue *> pair_for_nothing_test(res_1, res_2);
                  if (std::find(tested_but_nothing.begin(),
                                tested_but_nothing.end(),
                                pair_for_nothing_test) != tested_but_nothing.end()) {
                     // std::cout << "     this has already been tested and nothing was found" << std::endl;
                     continue;
                  }

                  if (false)
                     std::cout << "debug:: ---- for atom pair " << atom_spec_t(at_1) << " " << atom_spec_t(at_2)
                               << " considering residues " << residue_spec_t(res_1) << " to " << residue_spec_t(res_2)
                               << std::endl;

                  // Returns first (link_type) as "" if not found, second is order switch flag
                  std::pair<std::string, bool> lt = find_link_type_2022(res_1, res_2, geom);

                  if (false)
                     std::cout << "debug:: found " << lt.first.size() << " link types" << std::endl;

                  if (! lt.first.empty()) {
                     // std::cout << "DEBUG:: make_other_types_of_link() \"" << lt.first << "\""
                     //           << " for " << atom_spec_t(at_1) << " " << atom_spec_t(at_2)
                     //           << " detected" << std::endl;
                     logger.log(log_t::DEBUG, "make_other_types_of_link()", lt.first, "for",
                                atom_spec_t(at_1).format(), atom_spec_t(at_2).format(), "detected");
                  }

                  if (! lt.first.empty()) {
                     // this is not the place to make peptide links, event though find_link_type_complicado()
                     // will return a peptide link (for links that were not made before (perhaps because
                     // missing atoms)).
                     if ((lt.first != "TRANS") && (lt.first != "PTRANS") && (lt.first != "CIS") && (lt.first != "PCIS") && (lt.first != "p")) {
                        if (false)
                           std::cout << "DEBUG:: make_other_types_of_link(): now making a link restraint "
                                     << residue_spec_t(res_1) << " " << residue_spec_t(res_2)
                                     << " with type " << lt.first << " and order switch " << lt.second
                                     << std::endl;
                        bool fixed_1 = (fixed_atom_flags_set.find(i)   != fixed_atom_flags_set.end());
                        bool fixed_2 = (fixed_atom_flags_set.find(*it) != fixed_atom_flags_set.end());

                        nlrs.insert(at_1, at_2, res_1, res_2, fixed_1, fixed_2, lt.first, lt.second);
                     }
                  } else {
                     std::pair<mmdb::Residue *, mmdb::Residue *> pair_for_nothing_test_1(res_1, res_2);
                     std::pair<mmdb::Residue *, mmdb::Residue *> pair_for_nothing_test_2(res_2, res_1);
                     tested_but_nothing.push_back(pair_for_nothing_test_1);
                     tested_but_nothing.push_back(pair_for_nothing_test_2);
                  }
               }
            }
         }
      }
   }

   for (std::size_t i=0; i<nlrs.size(); i++)
      lrc.add(make_link_restraints_for_link_ng(nlrs[i], geom));

   return lrc;
}


// fill residue_link_vector_map_p and residue_pair_link_set_p
void
coot::restraints_container_t::make_link_restraints_ng(const coot::protein_geometry &geom,
                                                      bool do_rama_plot_restraints,
                                                      bool do_trans_peptide_restraints,
                                                      std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > *residue_link_vector_map_p,
                                                      std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > *residue_pair_link_set_p) {

   bool show_timings = false; // maybe pass this - or put it in the class
   auto tp_0 = std::chrono::high_resolution_clock::now();

   make_polymer_links_ng(geom, do_rama_plot_restraints, do_trans_peptide_restraints,
                         residue_link_vector_map_p, residue_pair_link_set_p);

   auto tp_1 = std::chrono::high_resolution_clock::now();

   // bonded_atom_indices is a vector of what is bonded or angle-bonded to each atom.

   make_flanking_atoms_restraints_ng(geom,
                                     residue_link_vector_map_p,
                                     residue_pair_link_set_p,
                                     do_rama_plot_restraints,  // not done here.
                                     do_trans_peptide_restraints);

   auto tp_2 = std::chrono::high_resolution_clock::now();

   // std::cout << "debug:: calling make_other_types_of_link() with links size " << links.size() << std::endl;
   link_restraints_counts others = make_other_types_of_link(geom, *residue_link_vector_map_p, *residue_pair_link_set_p);

   bonded_pair_container_t header_based_links = make_link_restraints_from_links(geom);

   if (verbose_geometry_reporting != QUIET)
      others.report();

   if (show_timings) {
      auto tp_3 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
      auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
      std::cout << "------------------ timings: for make_link_restraints_ng(): polymers: "
                << d10 << " flanking: " << d21 << " others: " << d32 << " milliseconds " << std::endl;
   }
}

void
coot::restraints_container_t::analyze_for_bad_restraints(analyze_bad_restraints_mode mode) {

   double interesting_distortion_limit = 5.0;
   analyze_for_bad_restraints(     CHIRAL_VOLUME_RESTRAINT,     interesting_distortion_limit);
   analyze_for_bad_restraints(              BOND_RESTRAINT,     interesting_distortion_limit);
   if (mode == BAD_RESTRAINT_ANALYSIS_INCLUDE_ANGLES)
      analyze_for_bad_restraints(          ANGLE_RESTRAINT,     interesting_distortion_limit);
   analyze_for_bad_restraints(NON_BONDED_CONTACT_RESTRAINT,     interesting_distortion_limit);
   analyze_for_bad_restraints(GEMAN_MCCLURE_DISTANCE_RESTRAINT, interesting_distortion_limit);

}

void
coot::restraints_container_t::analyze_for_bad_restraints(restraint_type_t r_type, double interesting_distortion_limit) {

   unsigned int n_baddies = 20;

   std::vector<std::tuple<unsigned int, double, double, double> > distortions; // index n_z bl_delta, target-value distortion
   for (unsigned int i=0; i<restraints_vec.size(); i++) {
      const simple_restraint &rest = restraints_vec[i];

      if (rest.restraint_type == r_type) {

         // distortion and bond length delta
         std::pair<double, double> distortion_pair = rest.distortion(atom, lennard_jones_epsilon); // 2nd arg is not used for bonds
         double distortion_score = distortion_pair.first;

         if (distortion_score >= interesting_distortion_limit) {
            double n_z = sqrt(distortion_score);
            double delta = distortion_pair.second;
            std::tuple<unsigned int, double, double, double> p(i, n_z, delta, distortion_score);
            distortions.push_back(p);
         }
      }
   }

   if (false)
      std::cout << "DEBUG:: in analyze_for_bad_restraints() n-restraints " << size() << " type " << r_type
                << " count " << distortions.size() << std::endl;

   auto distortion_sorter_lambda = // big distortions at the top
     [] (const std::tuple<unsigned int, double, double, double> &p1,
         const std::tuple<unsigned int, double, double, double> &p2) {
     return (std::get<3>(p1) > std::get<3>(p2)); };

   std::sort(distortions.begin(), distortions.end(), distortion_sorter_lambda);
   if (distortions.size() < n_baddies)
      n_baddies = distortions.size();
   std::cout << std::setw(3);
   for (unsigned int i=0; i<n_baddies; i++) {
      const std::tuple<unsigned int, double, double, double> &d = distortions[i];
      const simple_restraint &rest = restraints_vec[std::get<0>(d)];
      mmdb::Atom *at_1 = atom[rest.atom_index_1];
      mmdb::Atom *at_2 = atom[rest.atom_index_2];

      std::cout << std::setprecision(3);

      if (r_type == CHIRAL_VOLUME_RESTRAINT) {
         mmdb::Atom *at_c = atom[rest.atom_index_centre];
         // std::cout << "INFO:: Model: Bad Chiral Volume: "
         //           << atom_spec_t(at_c)
         //           << " delta "      << std::setw(6) << std::get<2>(d)
         //           << " target "     << std::setw(5) <<rest.target_chiral_volume
         //           << " sigma "      << std::setw(4) << rest.sigma
         //           << " distortion " << std::setw(4) << std::right << std::get<3>(d) << "\n";
         logger.log(log_t::INFO, {"Model: Bad Chiral Volume",
                                  atom_spec_t(at_c).format(), "delta", std::get<2>(d), "target", rest.target_chiral_volume,
                                  "sigma", rest.sigma, "distortion", std::get<3>(d)});
      }

      if (r_type == BOND_RESTRAINT) {
         // std::cout << "INFO:: Model: Bad Bond: " << std::setw(5)
         //           << atom_spec_t(at_1) << " to " << atom_spec_t(at_2)
         //           << " nZ "         << std::right << std::setprecision(2) << std::fixed << std::get<1>(d)
         //           << " delta "      << std::right << std::setprecision(3) << std::fixed << std::setw(6) << std::get<2>(d)
         //           << " target "     << std::right << std::setprecision(3) << std::fixed << rest.target_value
         //           << " sigma "      << std::right << std::setprecision(3) << std::fixed << rest.sigma
         //           << " distortion " << std::setw(6) << std::get<3>(d);
         if (rest.is_hydrogen_bond) {
            logger.log(log_t::INFO, {"Model: Bad Bond:",
                                     atom_spec_t(at_1).format(), "to", atom_spec_t(at_2).format(),
                                     "nZ", std::get<1>(d), "delta", std::get<2>(d), "target", rest.target_value,
                                     "sigma", rest.sigma, "distortion", std::get<3>(d)});
         } else {
            logger.log(log_t::INFO, {"Model: Bad Bond:",
                                     atom_spec_t(at_1).format(), "to", atom_spec_t(at_2).format(),
                                     "nZ", std::get<1>(d), "delta", std::get<2>(d), "target", rest.target_value,
                                     "sigma", rest.sigma, "distortion", std::get<3>(d), "H"});
         }
      }

      if (r_type == ANGLE_RESTRAINT) {
         // std::cout << "INFO:: Model: Bad Angle: " << std::setw(5)
         //           << atom_spec_t(at_1) << " to " << atom_spec_t(at_2)
         //           << " nZ "         << std::right << std::setprecision(2) << std::fixed << std::get<1>(d)
         //           << " delta "      << std::right << std::setprecision(3) << std::fixed << std::setw(6) << std::get<2>(d)
         //           << " target "     << std::right << std::setprecision(3) << std::fixed << rest.target_value
         //           << " sigma "      << std::right << std::setprecision(3) << std::fixed << rest.sigma
         //           << " distortion " << std::setw(6) << std::get<3>(d) << "\n";
         logger.log(log_t::INFO, {"Model: Bad Angle:",
                                  atom_spec_t(at_1).format(), "to", atom_spec_t(at_2).format(),
                                  "nZ", std::get<1>(d), "delta", std::get<2>(d), "target", rest.target_value,
                                  "sigma", rest.sigma, "distortion", std::get<3>(d)});
      }

      if (r_type == NON_BONDED_CONTACT_RESTRAINT) {
         // std::cout << "INFO:: Model: Bad Non-Bonded Contact: " << std::setw(5)
         //           << atom_spec_t(at_1) << " " << std::setw(3) << at_1->GetResName()
         //           << " to " << atom_spec_t(at_2) << " " << std::setw(3) << at_2->GetResName()
         //           << " delta "      << std::get<2>(d)
         //           << " target "     << rest.target_value
         //           << " distortion " << std::get<3>(d) << "\n";
         logger.log(log_t::INFO, {"Model: Bad Non-Bonded Contact:",
                                  atom_spec_t(at_1).format(), at_1->GetResName(),
                                  "to", atom_spec_t(at_2).format(), at_2->GetResName(),
                                  "delta", std::get<2>(d), "target", rest.target_value,
                                  "distortion", std::get<3>(d)});
      }

      if (r_type == GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
         // std::cout << "INFO:: Model: Bad GM Distance: " << std::setw(5)
         //           << atom_spec_t(at_1) << " to " << atom_spec_t(at_2) << " "
         //           << " delta "      << std::get<2>(d)
         //           << " target "     << rest.target_value
         //           << " distortion " << std::get<3>(d) << "\n";
         logger.log(log_t::INFO, {"Model: Bad GM Distance:",
                                  atom_spec_t(at_1).format(), "to", atom_spec_t(at_2).format(),
                                  "delta", std::get<2>(d), "target", rest.target_value,
                                  "distortion", std::get<3>(d)});
      }
   }
}
