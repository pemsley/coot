//
#include <iostream>
#include <iomanip>
#include <string>
#include <deque>
#include <thread>

#include "utils/split-indices.hh"
#include "utils/coot-fasta.hh"
#include "cootaneer/buccaneer-prot.h"  // for clipper's globuarise()
#include "mini-mol/mini-mol.hh"
#include "coot-utils/coot-map-utils.hh"
#include "scored-node.hh"
#include "ligand.hh"
#include "utils/coot-utils.hh"
#include "side-chain-densities.hh"

coot::minimol::molecule
get_flood_molecule(const clipper::Xmap<float> &xmap, float rmsd_cut_off) {

   bool debug = true;
   coot::ligand lig;
   float flood_atom_mask_radius = 1.1;

   lig.set_cluster_size_check_off();
   lig.set_chemically_sensible_check_off();
   lig.set_sphericity_test_off();

   lig.set_map_atom_mask_radius(flood_atom_mask_radius);
   lig.set_water_to_protein_distance_limits(10.0, 1.5);

   lig.import_map_from(xmap);

   lig.flood2(rmsd_cut_off);
   coot::minimol::molecule water_mol = lig.water_mol();

   if (debug) {
      std::string output_pdb = "flood-mol.pdb";
      water_mol.write_file(output_pdb, 30.0);
      lig.output_map("find-waters-masked-flooded.map");
   }

   return water_mol;
}

std::vector<std::pair<unsigned int, unsigned int> >
atom_pairs_within_distance(mmdb::Manager *mol_in,
                           mmdb::Atom **atom_selection, int n_selected_atoms,
                           double trans_peptide_CA_CA_dist,
                           double trans_peptide_CA_CA_dist_variation) {

   // set class members mol, atom_selection and n_sel_atoms.
   mmdb::Manager *mol = mol_in;  // the peaks in the map - some of which are CAs hopefully.

   std::vector<std::pair<unsigned int, unsigned int> > v;
   if (mol) {
      int uddHnd = mol->RegisterUDInteger(mmdb::UDR_ATOM, "index");
      if (uddHnd<0)  {
         std::cout << " atom bonding registration failed.\n";
      } else {
         for (int i=0; i< n_selected_atoms; i++) {
            mmdb::Atom *at = atom_selection[i];
            at->PutUDData(uddHnd, i); // is this needed any more?
         }

         mmdb::Contact *pscontact = NULL;
         int n_contacts;
         long i_contact_group = 1;
         mmdb::mat44 my_matt;
         for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
               my_matt[i][j] = 0.0;
         for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
         //
         mmdb::realtype local_dist_max = trans_peptide_CA_CA_dist + trans_peptide_CA_CA_dist_variation;
         mmdb::realtype local_dist_min = trans_peptide_CA_CA_dist - trans_peptide_CA_CA_dist_variation;

         std::cout << "debug:: SeekContacts with distance limits "
                   << local_dist_min << " " << local_dist_max << std::endl;

         mol->SeekContacts(atom_selection, n_selected_atoms,
                           atom_selection, n_selected_atoms,
                           local_dist_min, local_dist_max,
                           0,        // seqDist 0 -> in same res also
                           pscontact, n_contacts,
                           0, &my_matt, i_contact_group);

         if (n_contacts > 0) {
            if (pscontact) {
               for (int i=0; i<n_contacts; i++) {
                  if (pscontact[i].id1 < pscontact[i].id2) {
                     mmdb::Atom *at_1 = atom_selection[pscontact[i].id1];
                     mmdb::Atom *at_2 = atom_selection[pscontact[i].id2];
                     int idx_1, idx_2;
                     at_1->GetUDData(uddHnd, idx_1);
                     at_2->GetUDData(uddHnd, idx_2);
                     std::pair<unsigned int, unsigned int> p(idx_1, idx_2);
                      v.push_back(p);
                  }
               }
            }
         }
         std::cout << "INFO:: found " << n_contacts << " potential distance pairs " << std::endl;
         std::cout << "INFO:: made  " << v.size() << " uniqued distance pairs " << std::endl;
      }
   }
   return v;

}

std::pair<unsigned int, coot::scored_node_t>
spin_score(unsigned int idx_1, unsigned int idx_2, mmdb::Atom **atom_selection,
           const clipper::Xmap<float> &xmap, float map_rmsd) {

   auto index_to_pos = [atom_selection] (unsigned int idx) {
                          mmdb::Atom *at = atom_selection[idx];
                          return clipper::Coord_orth(at->x, at->y, at->z);
                       };

   float scale_CO       =  1.95;
   float scale_CO_low   = -0.6;
   float scale_CO_anti  = -0.15;
   float scale_N        =  1.2;
   float scale_perp     = -0.8;
   float scale_mid      =  1.6;
   float scale_non_line =  0.4;

   mmdb::Atom *at_1 = atom_selection[idx_1];
   mmdb::Atom *at_2 = atom_selection[idx_2];

   const clipper::Coord_orth pos_1 = index_to_pos(idx_1);
   const clipper::Coord_orth pos_2 = index_to_pos(idx_2);

   // draw a line between pos_1 and pos_2

   // find a point A that is 1.56 down the stick and 1.57 away from the line
   // find a point B that is 1.9  down the stick and 1.91 away from the line
   // find a point C that is 1.9  down the stick and 1.91 away from the line
   //      in the opposite direction to B.

   // spin this points A, B and C around the pos_1 - pos_2 line and the score is
   // rho(A) - rho(B) - rho(C)

   // Note to self: also "above" and "below" the peptide plane we expect to have
   // little to no density - so that should be added to the scoring system too.

   clipper::Coord_orth arb(0,0,1);
   clipper::Coord_orth diff_p(pos_2 - pos_1);
   clipper::Coord_orth diff_p_unit(diff_p.unit());

   clipper::Coord_orth perp(clipper::Coord_orth::cross(arb, diff_p));
   clipper::Coord_orth perp_unit(perp.unit());

   clipper::Coord_orth double_perp(clipper::Coord_orth::cross(diff_p, perp));
   clipper::Coord_orth double_perp_unit(double_perp.unit());

   double along_CA_CA_pt_O = 1.53; // the C is lower down than the O.
   double along_CA_CA_pt_for_perp = 2.33;

   double along_CA_CA_pt_N = 2.5;
   double ideal_peptide_length = 3.81;

   // we don't want the peptide to be scrunched up on one side of a
   // "long" peptide... let the atom positions expand along a long peptide. 
   double diff_p_len = sqrt(diff_p.lengthsq());
   double f_ca_ca_o = along_CA_CA_pt_O * diff_p_len/ideal_peptide_length;
   // double f_ca_ca_c = along_CA_CA_pt_C * diff_p_len/ideal_peptide_length;
   double f_ca_ca_n = along_CA_CA_pt_N * diff_p_len/ideal_peptide_length;
   double f_ca_ca_pt_for_perp = along_CA_CA_pt_for_perp * diff_p_len/ideal_peptide_length;

   // clipper::Coord_orth rel_line_pt_C(diff_p_unit * f_ca_ca_c + perp_unit * 0.7);
   // clipper::Coord_orth rel_line_pt_N(diff_p_unit * f_ca_ca_n - perp_unit * 0.5);

   // there is good density 1.9A away from the mid-line in the direction of the CO
   // (at the O). 
   // there is little density 3.7A away from the mid-line in the direction of the CO
   clipper::Coord_orth rel_line_pt_O(      diff_p_unit * f_ca_ca_o + perp_unit * 1.89);
   clipper::Coord_orth rel_line_pt_O_low(  diff_p_unit * f_ca_ca_o + perp_unit * 3.9);
   clipper::Coord_orth rel_line_pt_CO_anti(diff_p_unit * f_ca_ca_o - perp_unit * 0.5);
   clipper::Coord_orth rel_line_pt_N(      diff_p_unit * f_ca_ca_n - perp_unit * 0.3);
   clipper::Coord_orth rel_line_pt_perp1(diff_p_unit * f_ca_ca_pt_for_perp  + double_perp_unit * 1.85);
   clipper::Coord_orth rel_line_pt_perp2(diff_p_unit * f_ca_ca_pt_for_perp  - double_perp_unit * 1.72);


   // Idea from Andrea T (originally from George, I understand)
   // 
   // There should be density for a hydrogen bond acceptor in the
   // extension of the direction to the N from the CA-CA line.
   // However, we need to spin-search this a bit because it's often somewhat off axis.
   // 
   // Currently doesn't work.
   clipper::Coord_orth rel_line_pt_N_accpt(diff_p_unit * f_ca_ca_n - perp_unit * 3.0);
   clipper::Coord_orth rel_line_pt_N_accpt_off(diff_p_unit * (f_ca_ca_n + 0.5) - perp_unit * 3.0);
   
   float rho_at_1 = coot::util::density_at_point(xmap, pos_1);
   float rho_at_2 = coot::util::density_at_point(xmap, pos_2);

   // the mid-point between CAs should have density too.
   clipper::Coord_orth pt_mid(pos_1 * 0.50 + pos_2 * 0.5);
   float rho_mid = coot::util::density_at_point(xmap, pt_mid);

   int n_steps = 36;
   float best_score = -999;

   float rho_CO_best      = -999; // for testing scoring
   float rho_CO_low_best  = -999;
   float rho_CO_anti_best = -999;
   float rho_N_best       = -999;
   float rho_perp_1_best  = -999;
   float rho_perp_2_best  = -999;
   double alpha_best = -1; // negative indicates this was not set - which is a bad thing.

   // // these can be optimized with machine learning?
   // float scale_CO       =  0.5;
   // float scale_CO_low   = -0.6;
   // float scale_CO_anti  = -0.1;
   // float scale_perp     = -0.7;
   // float scale_mid      =  1.6;
   // float scale_non_line =  1.0;
   // float scale_N        = -0.0;

   // unsigned int idx_test_1 = 1228;
   // unsigned int idx_test_2 = 1238;
   // unsigned int idx_test_1 = 101;
   // unsigned int idx_test_2 = 109;
   
   unsigned int idx_test_1 = 131;
   unsigned int idx_test_2 = 139;
   
   for (int i=0; i< int(n_steps); i++) { 
      double alpha = 2 * M_PI * double(i)/double(n_steps);

      // direction position orig-shift angle
      // 
      clipper::Coord_orth p_CO = coot::util::rotate_around_vector(diff_p_unit,
                                                            pos_1 + rel_line_pt_O,
                                                            pos_1, alpha);
      
      clipper::Coord_orth p_CO_low = coot::util::rotate_around_vector(diff_p_unit,
                                                                pos_1 + rel_line_pt_O_low,
                                                                pos_1, alpha);
      
      clipper::Coord_orth p_CO_anti = coot::util::rotate_around_vector(diff_p_unit,
                                                                 pos_1 + rel_line_pt_CO_anti,
                                                                 pos_1, alpha);
      
      clipper::Coord_orth p_N = coot::util::rotate_around_vector(diff_p_unit,
                                                           pos_1 + rel_line_pt_N,
                                                           pos_1, alpha);
      
      clipper::Coord_orth p_2 = coot::util::rotate_around_vector(diff_p_unit,
                                                           pos_1 + rel_line_pt_perp1,
                                                           pos_1, alpha);
      clipper::Coord_orth p_3 = coot::util::rotate_around_vector(diff_p_unit,
                                                           pos_1 + rel_line_pt_perp2,
                                                           pos_1, alpha);

      // clipper::Coord_orth p_N_acceptor_1 = util::rotate_round_vector(diff_p_unit,
      //                                                                      pos_1 + rel_line_pt_N_accpt,
      //                                                                      pos_1, alpha);
      // clipper::Coord_orth p_N_acceptor_2 = util::rotate_round_vector(diff_p_unit,
      //                                                                      pos_1 + rel_line_pt_N_accpt,
      //                                                                      pos_1, (alpha + 0.26));
      // clipper::Coord_orth p_N_acceptor_3 = util::rotate_round_vector(diff_p_unit,
      //                                                                      pos_1 + rel_line_pt_N_accpt,
      //                                                                      pos_1, (alpha - 0.26));
      // clipper::Coord_orth p_N_acceptor_4 = util::rotate_round_vector(diff_p_unit,
      //                                                                      pos_1 + rel_line_pt_N_accpt_off,
      //                                                                      pos_1, alpha);
      // clipper::Coord_orth p_N_acceptor_5 = util::rotate_round_vector(diff_p_unit,
      //                                                                      pos_1 + rel_line_pt_N_accpt_off,
      //                                                                      pos_1, (alpha + 0.26));
      // clipper::Coord_orth p_N_acceptor_6 = util::rotate_round_vector(diff_p_unit,
      //                                                                      pos_1 + rel_line_pt_N_accpt_off,
      //                                                                      pos_1, (alpha - 0.26));
      
      
      float rho_CO      = coot::util::density_at_point(xmap, p_CO);
      float rho_CO_low  = coot::util::density_at_point(xmap, p_CO_low);
      float rho_CO_anti = coot::util::density_at_point(xmap, p_CO_anti);
      float rho_N       = coot::util::density_at_point(xmap, p_N);
      float rho_perp_1  = coot::util::density_at_point(xmap, p_2);
      float rho_perp_2  = coot::util::density_at_point(xmap, p_3);

      // float rho_acceptor_1 = util::density_at_point(xmap, p_N_acceptor_1);
      // float rho_acceptor_2 = util::density_at_point(xmap, p_N_acceptor_2);
      // float rho_acceptor_3 = util::density_at_point(xmap, p_N_acceptor_3);
      // float rho_acceptor_4 = util::density_at_point(xmap, p_N_acceptor_4);
      // float rho_acceptor_5 = util::density_at_point(xmap, p_N_acceptor_5);
      // float rho_acceptor_6 = util::density_at_point(xmap, p_N_acceptor_6);

      // float rho_acceptor_best = rho_acceptor_1;
      // if (rho_acceptor_2 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_2;
      // if (rho_acceptor_3 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_3;
      // if (rho_acceptor_4 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_4;
      // if (rho_acceptor_5 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_5;
      // if (rho_acceptor_6 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_6;

      float this_score =
         scale_CO      * rho_CO      +
         scale_CO_low  * rho_CO_low  + 
         scale_CO_anti * rho_CO_anti +
         scale_N       * rho_N       +
         scale_perp    * rho_perp_1  +
         scale_perp    * rho_perp_2;       
      // scale_N_accpt * rho_acceptor_best


      if (idx_1 == idx_test_1 && idx_2 == idx_test_2) {
         std::cout << "debug_pos:: CO     " << p_CO.x() << " " << p_CO.y() << " " << p_CO.z()
                   << " " << rho_CO << std::endl;
         std::cout << "debug_pos:: CO_low " << p_CO_low.x() << " " << p_CO_low.y()
                   << " " << p_CO_low.z() << " " << rho_CO_low << std::endl;
         std::cout << "debug_pos:: CO_anti " << p_CO_anti.x() << " " << p_CO_anti.y()
                   << " " << p_CO_anti.z() << " " << rho_CO_anti << std::endl;
         std::cout << "debug_pos:: perp-1 " << p_2.x() << " " << p_2.y() << " " << p_2.z()
                   << " " << rho_perp_1 << std::endl;
         std::cout << "debug_pos:: perp-2 " << p_3.x() << " " << p_3.y() << " " << p_3.z()
                   << " " << rho_perp_2 << std::endl;
         std::cout << "debug_pos:: N " << p_N.x() << " " << p_N.y() << " " << p_N.z()
                   << " " << rho_N << std::endl;
      } 

      if (this_score > best_score) { 
         best_score       = this_score;
         rho_CO_best      = rho_CO;
         rho_CO_low_best  = rho_CO_low;
         rho_CO_anti_best = rho_CO_anti;
         rho_N_best       = rho_N;
         rho_perp_1_best  = rho_perp_1;
         rho_perp_2_best  = rho_perp_2;
         alpha_best       = alpha;
      }
   }

   bool output_density_values = false;
   if (output_density_values) {
      std::cout << "debug-rho:: CO     "  << rho_CO_best      << " "
                << "debug-rho:: CO_low "  << rho_CO_low_best  << " "
                << "debug-rho:: CO_anti " << rho_CO_anti_best << " "
                << "debug-rho:: perp-1 "  << rho_perp_1_best  << " "
                << "debug-rho:: perp-2 "  << rho_perp_2_best  << " "
                << "debug-rho:: N "       << rho_N_best << std::endl;
   }


   float non_line_equal_density_penalty_1 = rho_at_1 + rho_at_2 - 2 * rho_mid;
   float non_line_equal_density_penalty = 
      - non_line_equal_density_penalty_1 * non_line_equal_density_penalty_1/map_rmsd;

   best_score += scale_mid * rho_mid;
   best_score += scale_non_line * non_line_equal_density_penalty;

   coot::scored_node_t best_node(idx_2, best_score, alpha_best);

   bool using_test_model = false;
   // for testing
   if (using_test_model) {
      std::string atom_name_1(at_1->name);
      std::string atom_name_2(at_2->name);
      if (atom_name_1 == " CA ") {
         if (atom_name_2 == " CA ") {
            if ((at_1->GetSeqNum() + 1) == at_2->GetSeqNum())
               best_node.udd_flag = true;
         }
      }
   }
   return std::pair<unsigned int, coot::scored_node_t> (idx_1, best_node);
}


// return sorted scores
// 
std::vector<std::pair<unsigned int, coot::scored_node_t> > 
spin_score_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &atom_pairs_within_distance,
                 const clipper::Xmap<float> &xmap,
                 mmdb::Manager *mol, mmdb::Atom **atom_selection, int n_selected_atoms) {

   unsigned int n_top = 400; // top-scoring spin-score pairs for tracing. Pass this, or compute it

   // apwd : atom (index) pairs within distance

   std::vector<std::pair<unsigned int, coot::scored_node_t> > scores;
   if (! mol) return scores;

   bool debug = false;

   scores.resize(atom_pairs_within_distance.size()*2); // results go here

   std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
   float map_rmsd = sqrt(mv.second);

   unsigned int n_atom_pairs = atom_pairs_within_distance.size();
   unsigned int n_threads = coot::get_max_number_of_threads();
   std::vector<std::pair<unsigned int, unsigned int> > air = coot::atom_index_ranges(n_atom_pairs, n_threads);

   auto spin_score_workpackage = [] (std::pair<unsigned int, unsigned int> range,
                                     const std::vector<std::pair<unsigned int, unsigned int> > &atom_pairs_within_distance,
                                     mmdb::Atom **atom_selection,
                                     const clipper::Xmap<float> &xmap,
                                     float map_rmsd,
                                     std::vector<std::pair<unsigned int, coot::scored_node_t> > &scores) { // fill the scores

                                    for (unsigned int i=range.first; i<range.second; i++) {
                                       const unsigned &at_idx_1 = atom_pairs_within_distance[i].first;
                                       const unsigned &at_idx_2 = atom_pairs_within_distance[i].second;
                                       scores[2*i  ] = spin_score(at_idx_1, at_idx_2, atom_selection, xmap, map_rmsd);
                                       scores[2*i+1] = spin_score(at_idx_2, at_idx_1, atom_selection, xmap, map_rmsd);
                                       scores[2*i  ].second.reverse_spin_score = std::make_pair(true, scores[2*i+1].second.spin_score);
                                       scores[2*i+1].second.reverse_spin_score = std::make_pair(true, scores[2*i  ].second.spin_score);
                                    }
                                 };

#if 0
   for (unsigned int i=0; i<n_atom_pairs; i++) {
      const unsigned &at_idx_1 = atom_pairs_within_distance[i].first;
      const unsigned &at_idx_2 = atom_pairs_within_distance[i].second;
      scores[2*i ]  = spin_score(at_idx_1, at_idx_2, atom_selection, xmap, map_rmsd);
      scores[2*i+1] = spin_score(at_idx_2, at_idx_1, atom_selection, xmap, map_rmsd);
   }
#endif

   std::vector<std::thread> threads;
   for (unsigned int i=0; i<air.size(); i++) {
      const auto &range = air[i];
      threads.push_back(std::thread(spin_score_workpackage, range, std::cref(atom_pairs_within_distance), atom_selection,
                                    std::cref(xmap), map_rmsd, std::ref(scores)));
   }
   for (unsigned int i=0; i<air.size(); i++)
      threads[i].join();
   

   std::sort(scores.begin(), scores.end(), coot::scored_node_t::sort_pair_scores);
   if (scores.size() > n_top)
      scores.resize(n_top);


   if (debug) {

      // I want to see a histogram of the scores
      for (const auto &score : scores) {
         std::cout << "score " << score.first << " " << score.second.atom_idx << " "
                   << score.second.spin_score << " " << score.second.alpha << "\n";
      }

      auto index_to_pos = [atom_selection] (unsigned int idx) {
                             mmdb::Atom *at = atom_selection[idx];
                             return clipper::Coord_orth(at->x, at->y, at->z);
                          };
      auto index_to_name = [atom_selection] (unsigned int idx) {
                              mmdb::Atom *at = atom_selection[idx];
                              return std::string(at->GetAtomName());
                           };
      std::cout << "---- sorted scores ----- " << n_top <<  std::endl;
      for (unsigned int i=0; i<n_top; i++) {
         const std::string &at_name_1 = index_to_name(scores[i].first);
         const std::string &at_name_2 = index_to_name(scores[i].second.atom_idx);
         int res_no_1 = atom_selection[scores[i].first          ]->GetSeqNum();
         int res_no_2 = atom_selection[scores[i].second.atom_idx]->GetSeqNum();
         std::string chain_id_1 = atom_selection[scores[i].first          ]->GetChainID();
         std::string chain_id_2 = atom_selection[scores[i].second.atom_idx]->GetChainID();

         std::cout << "sorted spin scores " << " " << std::setw(4) << scores[i].first << " "
                   << " to " << std::setw(4) << scores[i].second.atom_idx << " "
                   << at_name_1 << " " << std::setw(3) << res_no_1 << " " << chain_id_1 << " "
                   << at_name_2 << " " << std::setw(3) << res_no_2 << " " << chain_id_2 << " "
                   << index_to_pos(scores[i].first).format() << " "
                   << index_to_pos(scores[i].second.atom_idx).format() << " "
                   << scores[i].second.spin_score << std::endl;
      }
   }

   return scores;
}

#include <fstream>

//  edit the passed tree
void
remove_tree_front(std::deque<std::pair<unsigned int, coot::scored_node_t> > &tree, unsigned int until_idx) {

   // std::cout << "::::::: remove_tree_front start " << tree.size() << " with until_idx " << until_idx << std::endl;

   bool made_a_deletion = true; // to start the loop

   while (made_a_deletion) {
      if (tree.empty())
         break;
      made_a_deletion = false;
      for (const auto &item : tree) {
         if (item.first == until_idx)
            break;
         if (item.first != until_idx) {
            tree.pop_front();
            made_a_deletion = true;
            break;
         }
      }
   }
   // std::cout << ":::::::::::::::::::: remove_tree_front end   " << tree.size() << std::endl;

}

//  edit the passed tree
void
remove_tree_back(std::deque<std::pair<unsigned int, coot::scored_node_t> > &tree, unsigned int until_idx) {

   // std::cout << "::::::: remove_tree_back start, now size " << tree.size() << " with until_idx " << until_idx << std::endl;

   // for (const auto &item : tree)
   // std::cout << "   Pre-Tree " << item.first << " " << item.second.atom_idx << " " << item.second.name << std::endl;
   
   bool made_a_deletion = true; // to start the loop

   while (made_a_deletion) {
      if (tree.empty())
         break;
      made_a_deletion = false;
      const auto &item = tree.back();
      if (item.second.atom_idx == until_idx) {
         break;
      } else {
         tree.pop_back();
         made_a_deletion = true;
      }
   }
   std::cout << "---" << std::endl;
   
   // for (const auto &item : tree)
   // std::cout << "  Post-Tree " << item.first << " " << item.second.atom_idx << " " << item.second.name << std::endl;
   
   // std::cout << "::::::: remove_tree_back done, now size " << tree.size() << std::endl;
}

mmdb::Manager *
make_fragments(std::vector<std::pair<unsigned int, coot::scored_node_t> > &scored_pairs,
               mmdb::Atom **atom_selection,
               const clipper::Xmap<float> &xmap) { // pass the atom_selection for debugging

   // maybe it will be faster to calculate which residues are connected to which other residues:
   //
   // std::map<unsigned int, std::vector<std::pair<unsigned int, coot::scored_node_t> > > neighbour_map
   // that can be calculated in parallel (with locking)
   //

   // can the end of one tree be added to another tree?
   // (keep going until no more merging is possible)
   //
   // edit trees

   typedef std::deque<std::pair<unsigned int, coot::scored_node_t> > tree_t;

   auto distance_check_min = [atom_selection] (int atom_index_1, int atom_index_2, double dist_min) {
                                clipper::Coord_orth pt_1 = coot::co(atom_selection[atom_index_1]);
                                clipper::Coord_orth pt_2 = coot::co(atom_selection[atom_index_2]);
                                double dist_min_sqrd = dist_min * dist_min;
                                double dist_sqrd = (pt_1-pt_2).lengthsq();
                                // std::cout << "distance check " << std::sqrt(dist_sqrd) << std::endl;
                                return dist_sqrd > dist_min_sqrd;
                         };


   auto grow_trees = [distance_check_min] (std::vector<std::pair<unsigned int, coot::scored_node_t> > &scored_pairs) {

                        std::vector<std::deque<std::pair<unsigned int, coot::scored_node_t> > > trees; // return this

                        // and keep a note of which trees are are edited copies of which other trees.
                        // Note that (say) 3 is a copy of 1, and 1 is a copy of 3.
                        // return this also.
                        std::map<unsigned int, std::set<unsigned int> >  tree_copies;

                        for (unsigned int i=0; i<scored_pairs.size(); i++) {
                           const auto &scored_pair = scored_pairs[i];

                           // std::cout << "------ considering pair " << scored_pair.first << " " << scored_pair.second.atom_idx
                           //           << " " << scored_pair.second.name << std::endl;

                           // can I add this to a tree?
                           bool added_to_a_tree = false;
                           for (unsigned int j=0; j<trees.size(); j++) {

                              auto &tree = trees[j];
                              unsigned int new_tree_index = trees.size(); // this will be the index of a new tree

                              // first, try to add to each of the ends of the tree
                              // (ideal case)

                              const auto &tree_node_b = tree.back();

                              if (tree_node_b.second.atom_idx == scored_pair.first) {
                                 if (tree_node_b.first == scored_pair.second.atom_idx) {
                                    // skip the trivial reverse
                                 } else {
                                    if (false)
                                       std::cout << "Push back! current-back: " << tree.back().first << " " << tree.back().second.atom_idx
                                                 << " new-back: " << scored_pair.first << " " << scored_pair.second.atom_idx
                                                 << std::endl;
                                    // CA(n-1) to CA(n+1) is shortest for helicies (5.3)
                                    // if (distance_check_min(tree_node.first, scored_pair.second.atom_idx, 5.0)) { // 5.3 with some room for noise
                                    if (distance_check_min(tree_node_b.first, scored_pair.second.atom_idx, 5.0)) { // 5.3 with some room for noise
                                       tree.push_back(scored_pair);
                                       added_to_a_tree = true;
                                    }
                                 }
                              }

                              if (! added_to_a_tree) {
                                 const auto &tree_node_f = tree.front();
                                 if (tree_node_f.first == scored_pair.second.atom_idx) {
                                    if (tree_node_f.second.atom_idx == scored_pair.first) {
                                       // skip the reverse peptide of the current front
                                    } else {
                                       if (distance_check_min(tree_node_f.second.atom_idx, scored_pair.first, 5.0)) {
                                          if (false)
                                             std::cout << "Push front! current-front: " << tree.front().first << " " << tree.front().second.atom_idx
                                                       << " new-front: " << scored_pair.first << " " << scored_pair.second.atom_idx
                                                       << std::endl;
                                          tree.push_front(scored_pair);
                                          added_to_a_tree = true;
                                       }
                                    }
                                 }
                              }

                              // OK, can I connect it part-way through a tree?
                              // If that is the case, then make a new tree
                              // by copying the current one, and chopping off
                              // the part of the tree that is connecting
                              // to the residue that has the same connection
                              // as the residue that we are adding.
                              //
                              if (! added_to_a_tree) {
                                 // std::deque<std::pair<unsigned int, coot::scored_node_t> > it;
                                 for (auto it=tree.begin(); it!=tree.end(); ++it) {
                                    if (it->first == scored_pair.second.atom_idx) {
                                       if (distance_check_min(it->second.atom_idx, scored_pair.first, 5.0)) {
                                          auto new_tree = tree;
                                          remove_tree_front(new_tree, it->first);
                                          added_to_a_tree = true;
                                          new_tree.push_front(scored_pair); // maybe push_front()
                                          trees.push_back(new_tree);
                                          tree_copies[j].insert(new_tree_index);
                                          tree_copies[new_tree_index].insert(j);
                                          std::cout << "Middle-tree - path A " << std::endl;
                                          break;
                                       }
                                    }
                                 }
                                 if (! added_to_a_tree) {
                                    for (auto it=tree.begin(); it!=tree.end(); ++it) {
                                       if (it->second.atom_idx == scored_pair.first) {
                                          if (distance_check_min(it->first, scored_pair.second.atom_idx, 5.0)) {
                                             auto new_tree = tree;
                                             remove_tree_back(new_tree, it->first); // check this
                                             new_tree.push_back(scored_pair); // maybe push_back()
                                             added_to_a_tree = true;
                                             trees.push_back(new_tree);
                                             tree_copies[j].insert(new_tree_index);
                                             tree_copies[new_tree_index].insert(j);
                                             std::cout << "Middle-tree - path B " << std::endl;
                                             break;
                                          }
                                       }
                                    }
                                 }
                              }
                              if (added_to_a_tree) { // trees iterator is broken let's jump out of that loop
                                 break;
                              }

                           }
                           if (! added_to_a_tree) {
                              // OK, so start a new tree
                              // std::cout << ".......... new tree " << scored_pair.first << " " << scored_pair.second.atom_idx << std::endl;
                              std::deque<std::pair<unsigned int, coot::scored_node_t> > new_path;
                              new_path.push_back(scored_pair);
                              trees.push_back(new_path);
                           }
                        }
                        return std::make_pair(trees, tree_copies);
                     };


   class scored_tree_t {
   public:
      unsigned int index; // index in the vector of trees (same as vector of scored_trees)
      tree_t tree;
      double score; // best score - backwards score was checked.
      bool marked_for_deletion;
      std::string chain_id;
      scored_tree_t(unsigned int idx, const tree_t &tree_in, double s) : index(idx), tree(tree_in), score(s) { marked_for_deletion = false; }
      scored_tree_t() { index = 0; marked_for_deletion = false; score = 0.0; }
   };

   // are these peptides trees in consitent directions?
   auto merge_tree_to_back_of_tree = [] (scored_tree_t &changeable_tree, const scored_tree_t &tree_to_be_added) {
                                        for (const auto &item : tree_to_be_added.tree)
                                           changeable_tree.tree.push_back(item);
                                        changeable_tree.score += tree_to_be_added.score;
                                     };

   auto merge_tree_to_front_of_tree = [] (scored_tree_t &changeable_tree, const scored_tree_t &tree_to_be_added) {
                                         tree_t::const_reverse_iterator it;
                                         for(it=tree_to_be_added.tree.rbegin(); it!=tree_to_be_added.tree.rend(); ++it)
                                            changeable_tree.tree.push_front(*it);
                                         changeable_tree.score += tree_to_be_added.score;
                                      };

   auto join_the_trees = [merge_tree_to_back_of_tree, merge_tree_to_front_of_tree] (std::vector<scored_tree_t> &scored_trees) {

                            bool made_a_merge = true; // start the loop
                            while (made_a_merge) {
                               made_a_merge = false;
                               for (unsigned int i=0; i<scored_trees.size(); i++) {
                                  auto &tree_i = scored_trees[i];
                                  for (unsigned int j=0; j<scored_trees.size(); j++) {
                                     if (i != j) {
                                        const auto &tree_j = scored_trees[j];

                                        if (tree_i.tree.front().first == tree_j.tree.back().second.atom_idx) {
                                           if (tree_i.tree.front().second.atom_idx == tree_j.tree.back().first) {
                                              // trivial backwards - skip
                                           } else {
                                              if (true)
                                                 std::cout << "Merge these - front! " << std::setw(3) << i << " " << std::setw(3) << j
                                                           << " with tree_i size " << std::setw(3) << tree_i.tree.size()
                                                           << " and tree_j size " << std::setw(3) << tree_j.tree.size() << std::endl;
                                              made_a_merge = true;
                                              merge_tree_to_front_of_tree(tree_i, tree_j);
                                              std::vector<scored_tree_t>::iterator it_erase = scored_trees.begin() + j;
                                              scored_trees.erase(it_erase);
                                           }
                                        }
                                        if (!made_a_merge) {
                                           if (tree_i.tree.back().second.atom_idx == tree_j.tree.front().first) {
                                              if (tree_i.tree.back().first == tree_j.tree.front().second.atom_idx) {
                                                 // trivial backwards - skip
                                              } else {
                                                 if (true)
                                                    std::cout << "Merge these - back!  " << std::setw(3) << i << " " << std::setw(3) << j
                                                              << " with tree_i size " << std::setw(3) << tree_i.tree.size() 
                                                              << " and tree_j size " << std::setw(3) << tree_j.tree.size() << std::endl;
                                                 made_a_merge = true;
                                                 merge_tree_to_back_of_tree(tree_i, tree_j);
                                                 std::vector<scored_tree_t>::iterator it_erase = scored_trees.begin() + j;
                                                 scored_trees.erase(it_erase);
                                              }
                                           }
                                        }
                                     }
                                     if (made_a_merge) break; // messed up indexing
                                  }
                                  if (made_a_merge) break;
                               }
                            }
                         };

   auto delete_singleton_Ns = [] (coot::minimol::fragment &frag) {
                                 for (unsigned int jj=0; jj<frag.residues.size(); jj++) {
                                    auto &residue = frag.residues[jj];
                                    unsigned int n_atoms = residue.n_atoms();
                                    if (n_atoms == 1) {
                                       auto &atom = residue[0];
                                       if (atom.name == " N  ") {
                                          // delete all the atoms (just one atom)
                                          residue.atoms.clear();
                                       }
                                    }
                                 }
                              };

   auto spin_search_N_and_CB = [] (const clipper::Coord_orth &N_pos,
                                   const clipper::Coord_orth &CB_pos,
                                   const clipper::Coord_orth &CA_pos,
                                   const clipper::Coord_orth &C_pos,
                                   const clipper::Xmap<float> &xmap) {

                                  // rotate about C-CA to find best position of the N and CB

                                  clipper::Coord_orth dir = CA_pos - C_pos;
                                  const clipper::Coord_orth &origin_shift = CA_pos;
                                  float sum_best = -999.9;
                                  double ar_best = 0.0;
                                  for (double a=0.0; a<360.0; a+=3.0) {
                                     double ar = a * M_PI/180.0;
                                     clipper::Coord_orth N_pos_new  = coot::util::rotate_around_vector(dir,  N_pos, origin_shift, ar);
                                     clipper::Coord_orth CB_pos_new = coot::util::rotate_around_vector(dir, CB_pos, origin_shift, ar);
                                     float d1 = coot::util::density_at_point(xmap,  N_pos_new);
                                     float d2 = coot::util::density_at_point(xmap, CB_pos_new);
                                     float sum = d1 + d2;
                                     if (sum > sum_best) {
                                        sum_best = sum;
                                        ar_best = ar;
                                     }
                                  }
                                  clipper::Coord_orth N_pos_new  = coot::util::rotate_around_vector(dir,  N_pos, origin_shift, ar_best);
                                  clipper::Coord_orth CB_pos_new = coot::util::rotate_around_vector(dir, CB_pos, origin_shift, ar_best);
                                  return std::make_pair(N_pos_new, CB_pos_new);
                               };

   auto invent_Ns_and_CBs_by_spin_search = [spin_search_N_and_CB] (coot::minimol::fragment &frag, const clipper::Xmap<float> &xmap) {
                                              double torsion_relative_to_N = 111.5 * M_PI/180.0;
                                              for (unsigned int jj=0; jj<frag.residues.size(); jj++) {
                                                 auto &residue = frag.residues[jj];
                                                 bool found_N  = false;
                                                 bool found_CB = false;
                                                 for (unsigned int iat=0; iat<residue.atoms.size(); iat++) {
                                                    const auto &atom = residue.atoms[iat];
                                                    if (atom.name == " N  ") {
                                                       found_N = true;
                                                       break;
                                                    }
                                                 }
                                                 if (!found_N && !found_CB) {
                                                    bool found_C  = false;
                                                    bool found_CA = false;
                                                    clipper::Coord_orth CA_pos;
                                                    clipper::Coord_orth  C_pos;
                                                    for (unsigned int iat=0; iat<residue.atoms.size(); iat++) {
                                                       const auto &atom = residue.atoms[iat];
                                                       if (atom.name == " C  ") {
                                                          C_pos = atom.pos;
                                                          found_C = true;
                                                       }
                                                       if (atom.name == " CA ") {
                                                          CA_pos = atom.pos;
                                                          found_CA = true;
                                                       }
                                                    }
                                                    if (found_C && found_CA) {
                                                       double torsion = 0.6;
                                                       double bond_length = 1.482;
                                                       double angle = 109.63 * M_PI/180.0;
                                                       clipper::Coord_orth p0(0,0,0);
                                                       // update p0 if possible:
                                                       if ((jj+1) < frag.residues.size()) {
                                                          const auto &next_residue = frag.residues[jj+1];
                                                          for (unsigned int iat=0; iat<next_residue.atoms.size(); iat++) {
                                                             const auto &atom = next_residue.atoms[iat];
                                                             if (atom.name == " N  ") {
                                                                p0 = atom.pos;
                                                             }
                                                          }
                                                       }
                                                       clipper::Coord_orth  N_pos(p0, C_pos, CA_pos, bond_length, angle, torsion);
                                                       clipper::Coord_orth CB_pos(p0, C_pos, CA_pos, bond_length, angle, torsion + torsion_relative_to_N);
                                                       auto new_positions = spin_search_N_and_CB(N_pos, CB_pos, CA_pos, C_pos, xmap);
                                                       coot::minimol::atom  N_at(" N  ", " N", new_positions.first,  "", 30.0);
                                                       coot::minimol::atom CB_at(" CB ", " C", new_positions.second, "", 30.0);
                                                       residue.addatom( N_at);
                                                       residue.addatom(CB_at);
                                                    }
                                                 }
                                              }
                                           };

   auto invent_extra_Ns_if_needed = [] (coot::minimol::fragment &frag) {
                                       for (unsigned int jj=0; jj<frag.residues.size(); jj++) {
                                          auto &residue = frag.residues[jj];
                                          bool found_N  = false;
                                          for (unsigned int iat=0; iat<residue.atoms.size(); iat++) {
                                             const auto &atom = residue.atoms[iat];
                                             if (atom.name == " N  ") {
                                                found_N = true;
                                                break;
                                             }
                                          }
                                          if (! found_N) {
                                             bool found_C  = false;
                                             bool found_CA = false;
                                             clipper::Coord_orth ca_pos;
                                             clipper::Coord_orth c_pos;
                                             for (unsigned int iat=0; iat<residue.atoms.size(); iat++) {
                                                const auto &atom = residue.atoms[iat];
                                                if (atom.name == " C  ") {
                                                   c_pos = atom.pos;
                                                   found_C = true;
                                                }
                                                if (atom.name == " CA ") {
                                                   ca_pos = atom.pos;
                                                   found_CA = true;
                                                }
                                             }
                                             if (found_C && found_CA) {
                                                double torsion = 0.6;
                                                double bond_length = 1.482;
                                                double angle = 109.63 * M_PI/180.0;
                                                clipper::Coord_orth p0(0,0,0);
                                                // update p0 if possible:
                                                if ((jj+1) < frag.residues.size()) {
                                                   const auto &next_residue = frag.residues[jj+1];
                                                   for (unsigned int iat=0; iat<next_residue.atoms.size(); iat++) {
                                                      const auto &atom = next_residue.atoms[iat];
                                                      if (atom.name == " N  ") {
                                                         p0 = atom.pos;
                                                      }
                                                   }
                                                }
                                                clipper::Coord_orth N_pos(p0, c_pos, ca_pos, bond_length, angle, torsion);
                                                coot::minimol::atom N_at(" N  ", " N", N_pos, "", 30.0);
                                                residue.addatom(N_at);
                                             }
                                          }
                                       }
                                    };

   // big trees at the top
   auto tree_length_sorter = [] (const std::deque<std::pair<unsigned int, coot::scored_node_t> > &tree_1,
                                 const std::deque<std::pair<unsigned int, coot::scored_node_t> > &tree_2) {
                                return tree_2.size() < tree_1.size();
                             };

   // sort scored_trees by score - best at the top of the list
   auto scored_tree_sorter = [] (const scored_tree_t &t_1, const scored_tree_t &t_2) {
                                return t_1.score > t_2.score;
                             };


   auto print_trees = [] (const std::vector<tree_t> &trees) {
                         for (unsigned int i=0; i<trees.size(); i++) {
                            const auto &tree = trees[i];
                            if (tree.size() > 1) {
                               std::cout << "--------- Tree pre-merge " << i << " size: " << tree.size()
                                         << " -------" << std::endl;
                               for (const auto &item : tree) {
                                  std::cout << "   Tree pre-merge " << item.first << " " << item.second.atom_idx << " "
                                            << item.second.name << std::endl;
                               }
                            }
                         }
                      };

   // not the same function as merge_fragments() uses (that needs overlapping fragments)
   //
   auto merge_fragments_by_joining_close_C_and_N = [] (coot::minimol::molecule &mol) {
                                                  };

   auto compare_tree_scores = [] (const tree_t &tree_trace, unsigned int tree_index) {
                                 double sum_forward  = 0.0;
                                 double sum_backward = 0.0;
                                 unsigned int n_forward  = 0;
                                 unsigned int n_backward = 0;
                                 for (unsigned int i=0; i<tree_trace.size(); i++) {
                                    const auto &pep = tree_trace[i];
                                    sum_forward += pep.second.spin_score;
                                    n_forward++;
                                    if (pep.second.reverse_spin_score.first) {
                                       n_backward++;
                                       sum_backward += pep.second.reverse_spin_score.second;
                                    }
                                 }
                                 if (false)
                                    std::cout << "INFO:: tree score for tree with index " << tree_index << " "
                                              << sum_forward << " from " << n_forward << " points and "
                                              << sum_backward << " from " << n_backward << " points, peptide-delta: "
                                              << (sum_forward-sum_backward)/static_cast<double>(n_forward) << std::endl;
                                 if (false) { // debug. Why do many 2 point trees give a score that is the same
                                             // for forward and backward? If the indexing of the second is the
                                             // reverse of the first, then we are failing to prevent trivial
                                             // node addition somewhere. OK, it's in join_the_trees().
                                    if (n_forward == 2) {
                                       const auto &pep_1 = tree_trace[0];
                                       const auto &pep_2 = tree_trace[1];
                                       std::cout << "node-1: " << std::setw(3) << pep_1.first << " node-2: " << std::setw(3) << pep_1.second.atom_idx << " scores: "
                                                 << pep_1.second.spin_score << " " << pep_1.second.reverse_spin_score.second << std::endl;
                                       std::cout << "node-1: " << std::setw(3) << pep_2.first << " node-2: " << std::setw(3) << pep_2.second.atom_idx << " scores: "
                                                 << pep_2.second.spin_score << " " << pep_2.second.reverse_spin_score.second << std::endl;
                                    }
                                 }
                              };

   auto index_to_chain_id = [] (unsigned int idx) {
                               char ch_id = 'A' + idx;
                               std::string chain_id(1, ch_id);
                               
                               // multiple chains
                               if (idx >= 26) {
                                  unsigned int i = idx;
                                  unsigned int j = 0;
                                  while (i > 25) {
                                     i -= 26;
                                     chain_id = std::string(1, 'A' + j) + std::string(1, 'A' + i);
                                     j += 1;
                                  }
                               }
                               return chain_id;
                            };

   auto add_chain_ids_to_scored_trees = [index_to_chain_id] (std::vector<scored_tree_t> &scored_trees) {
                                           for (unsigned int i=0; i<scored_trees.size(); i++) {
                                              scored_trees[i].chain_id = index_to_chain_id(i);
                                           }
                                        };

   auto tree_to_mc = [atom_selection] (const scored_tree_t &scored_tree, unsigned int ith_trace) {

                        std::string chain_id = scored_tree.chain_id;

                        coot::minimol::fragment fragment(chain_id);
                        for (unsigned int i=0; i<scored_tree.tree.size(); i++) {
                           const auto &pep = scored_tree.tree[i];
                           clipper::Coord_orth pt_1 = coot::co(atom_selection[pep.first]);
                           clipper::Coord_orth pt_2 = coot::co(atom_selection[pep.second.atom_idx]);
                           // std::swap(pt_1, pt_2);
                           clipper::Coord_orth diff_p = pt_2 - pt_1;
                           clipper::Coord_orth diff_p_unit(diff_p.unit());
                           clipper::Coord_orth arb(0,0,1);
                           clipper::Coord_orth perp(clipper::Coord_orth::cross(arb, diff_p));
                           clipper::Coord_orth perp_unit(perp.unit());
                           clipper::Coord_orth double_perp(clipper::Coord_orth::cross(diff_p, perp));
                           clipper::Coord_orth double_perp_unit(double_perp.unit());
                           double diff_p_len = std::sqrt(diff_p.lengthsq());
                           double ideal_peptide_length = 3.81;
                           double along_CA_CA_pt_O = 1.53; // the C is lower down than the O.
                           double along_CA_CA_pt_C = 1.46;
                           double along_CA_CA_pt_N = 2.44;
                           double f_ca_ca_o = along_CA_CA_pt_O * diff_p_len/ideal_peptide_length;
                           double f_ca_ca_c = along_CA_CA_pt_C * diff_p_len/ideal_peptide_length;
                           double f_ca_ca_n = along_CA_CA_pt_N * diff_p_len/ideal_peptide_length;
                           clipper::Coord_orth rel_line_pt_C(diff_p_unit * f_ca_ca_c + perp_unit * 0.48);
                           clipper::Coord_orth rel_line_pt_N(diff_p_unit * f_ca_ca_n - perp_unit * 0.47);
                           clipper::Coord_orth rel_line_pt_O(diff_p_unit * f_ca_ca_o + perp_unit * 1.89);
                           clipper::Coord_orth O_position_raw = pt_1 + rel_line_pt_O;
                           clipper::Coord_orth N_position_raw = pt_1 + rel_line_pt_N;
                           clipper::Coord_orth C_position_raw = pt_1 + rel_line_pt_C;
                           // now rotate O_position about diff_p
                           clipper::Coord_orth rotated_position_for_O = coot::util::rotate_around_vector(diff_p, O_position_raw, pt_1, pep.second.alpha);
                           clipper::Coord_orth rotated_position_for_N = coot::util::rotate_around_vector(diff_p, N_position_raw, pt_1, pep.second.alpha);
                           clipper::Coord_orth rotated_position_for_C = coot::util::rotate_around_vector(diff_p, C_position_raw, pt_1, pep.second.alpha);

                           coot::minimol::atom  O(" O  ", " O", rotated_position_for_O, "", 1.0, 30.0f);
                           coot::minimol::atom  N(" N  ", " N", rotated_position_for_N, "", 1.0, 30.0f);
                           coot::minimol::atom  C(" C  ", " C", rotated_position_for_C, "", 1.0, 30.0f);
                           coot::minimol::atom CA(" CA ", " C", pt_1, "", 30.0f);
                           coot::minimol::residue &r_s = fragment[i+2];  // large number first, so we don't get a resize between the lines.
                           coot::minimol::residue &r_f = fragment[i+1];
                           if (r_f.is_empty()) {
                              coot::minimol::residue r_1(i+1, "ALA");
                              r_1.addatom(C);
                              r_1.addatom(O);
                              r_1.addatom(CA);
                              fragment[i+1] = r_1;
                           } else {
                              r_f.addatom(CA);
                              r_f.addatom(C);
                              r_f.addatom(O);
                           }
                           if (r_s.is_empty()) {
                              coot::minimol::residue r_2(i+2, "ALA");
                              r_2.addatom(N);
                              fragment[i+2] = r_2;
                           } else {
                              r_s.addatom(N);
                           }
                           
                        }
                        return fragment;
                     };

   auto debug_trees = [atom_selection] (const std::vector<tree_t> &trees) {
                         auto index_to_colour = [] (unsigned int idx) {
                                                   coot::colour_holder c(0.9, 0.2, 0.2);
                                                   float f = 0.38f * idx;
                                                   c.rotate_by(f);
                                                   return c;
                                                };

                         auto index_to_pos = [atom_selection] (unsigned int idx) {
                                                mmdb::Atom *at = atom_selection[idx];
                                                return clipper::Coord_orth(at->x, at->y, at->z);
                                             };
                         std::ofstream f("debug-trace.points");
                         if (f) {
                            for (unsigned int i=0; i<trees.size(); i++) {
                               if (i > 4) continue;
                               coot::colour_holder col = index_to_colour(i);
                               const auto &tree = trees[i];
                               std::cout << "------- Tree " << i << " size: " << tree.size() << " ------" << std::endl;
                               if (tree.size() < 4) {
                               } else {
                                  for (const auto &item : tree) {
                                     clipper::Coord_orth pt_1 = index_to_pos(item.first);
                                     clipper::Coord_orth pt_2 = index_to_pos(item.second.atom_idx);
                                     std::cout << "   Tree " << item.first << " " << item.second.atom_idx << " "
                                               << item.second.name << std::endl;
                                     f << "trace-point " << col.red << " " << col.green << " " << col.blue << " "
                                       << pt_1.x() << " " << pt_1.y() << " " << pt_1.z() << " "
                                       << pt_2.x() << " " << pt_2.y() << " " << pt_2.z() << "\n";
                                  }
                               }
                            }
                         }
                         f.close();
                      };

   auto debug_oxygen_positions = [] (const std::vector<clipper::Coord_orth> &oxygen_positions) {
                                    coot::minimol::molecule mol(oxygen_positions, "DUM", " DUM", "A");
                                    mol.write_file("oxygen-positions.pdb", 10);
                                    std::ofstream f("oxygen-positions.table");
                                    for (const auto &item : oxygen_positions) {
                                       f << "oxygen-position " << item.x() << " " << item.y() << " " << item.z() << "\n";
                                    }
                                    f.close();
                                 };

   // pass the chain_id for debugging
   auto make_CB_ideal_pos = [] (const coot::minimol::residue &res, const std::string &chain_id) {
                               bool found = false;
                               clipper::Coord_orth pos;
                               // don't add one if there is one already there.
                               std::pair<bool, coot::minimol::atom> CB = res.get_atom(" CB ");
                               if (! CB.first) {
                                  auto CA = res.get_atom(" CA ");
                                  auto C  = res.get_atom(" C  ");
                                  auto N  = res.get_atom(" N  ");
                                  if (CA.first) {
                                     if (C.first) {
                                        if (N.first) {
                                           clipper::Coord_orth C_to_N = N.second.pos - C.second.pos;
                                           clipper::Coord_orth C_to_N_mid_point(0.5 * (N.second.pos + C.second.pos));
                                           clipper::Coord_orth CA_to_CN_mid_point = C_to_N_mid_point - CA.second.pos;
                                           clipper::Coord_orth CA_to_CN_mid_point_uv(CA_to_CN_mid_point.unit());
                                           clipper::Coord_orth perp(clipper::Coord_orth::cross(C_to_N, CA_to_CN_mid_point));
                                           clipper::Coord_orth perp_uv(perp.unit());
                                           // guess and fiddle these - good enough
                                           clipper::Coord_orth CB_pos(CA.second.pos + 1.21 * perp_uv - 0.95 * CA_to_CN_mid_point_uv);
                                           pos = CB_pos;
                                           found = true;
                                        } else {
                                           std::cout << "INFO:: make_CB_ideal_pos(): sad residue " << res << " in chain " << chain_id << ", has no N " << std::endl;
                                        }
                                     } else {
                                        std::cout << "INFO:: make_CB_ideal_pos(): sad residue " << res << " in chain " << chain_id << ", has no C " << std::endl;
                                     }
                                  } else {
                                     std::cout << "INFO:: make_CB_ideal_pos(): sad residue " << res << " in chain " << chain_id << ", has no CA " << std::endl;
                                  }
                               } else {
                                  std::cout << "INFO:: make_CB_ideal_pos(): residue " << res << " in chain " << chain_id << " " << res.seqnum
                                            << ", already has a CB " << std::endl;
                               }
                               return std::make_pair(found, pos);
                            };

   auto add_CBs_to_residues = [make_CB_ideal_pos] (coot::minimol::fragment &frag) {
                                 const std::string &frag_id = frag.fragment_id;
                                 for (int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
                                    if (frag[ires].n_atoms() > 0) {
                                       std::pair<bool, coot::minimol::atom> CB_from_res = frag[ires].get_atom(" CB ");
                                       if (! CB_from_res.first) {
                                          std::pair<bool, clipper::Coord_orth> cb_pos = make_CB_ideal_pos(frag[ires], frag_id);
                                          if (cb_pos.first) {
                                             coot::minimol::atom CB(" CB ", " C", cb_pos.second, "", 1.0, 30.0);
                                             // std::cout << "debug:: calling addatom for residue " << ires << " and self index " << frag[ires].seqnum
                                             // << " in chain " << frag_id << std::endl;
                                             frag[ires].addatom(CB);
                                             // std::cout << "debug:: add_CBs_to_residues(): added CB for " << frag_id << " " << ires << std::endl;
                                          } else {
                                             std::cout << "sadge: add_CBs_to_residues(): CB bool is false: " << frag_id << " " << ires
                                                       << " with n-atoms " << frag[ires].n_atoms() << std::endl;
                                          }
                                       }
                                    }
                                 }
                              };

   auto print_the_tree_copies = [] (std::vector<std::deque<std::pair<unsigned int, coot::scored_node_t> > > &trees,
                                    const std::map<unsigned int, std::set<unsigned int> > &tree_copies) {
                                   std::cout << "------- Tree Copies ------" << std::endl;
                                   for (const auto &entry : tree_copies) {
                                      const unsigned int &key = entry.first;
                                      const std::set<unsigned int> &set = entry.second;
                                      std::cout << "Tree " << key << " (" << trees[key].size() << ") : ";
                                      std::set<unsigned int>::const_iterator it;
                                      for (it=set.begin(); it!=set.end(); ++it) {
                                         std::cout << *it << " (" << trees[*it].size() << ") ";
                                      }
                                      std::cout << "\n";
                                   }
                                };

   // check the revers score also, and use that if it is better
   //
   auto get_tree_score = [] (const tree_t &tree) {
                            double sum_forward = 0.0;
                            double sum_backward = 0.0;
                            for (unsigned int i=0; i<tree.size(); i++) {
                               const auto &pep = tree[i];
                               sum_forward += pep.second.spin_score;
                               if (pep.second.reverse_spin_score.first)
                                  sum_backward += pep.second.reverse_spin_score.second;
                            }
                            if (sum_backward > sum_forward)
                               sum_forward = sum_backward;
                            return sum_forward;
                          };

   auto make_scored_trees = [get_tree_score] (const std::vector<tree_t> &trees) {
                               std::vector<scored_tree_t> scored_trees(trees.size());
                               for (unsigned int i=0; i<trees.size(); i++) {
                                  const tree_t &tree = trees[i];
                                  double score = get_tree_score(tree);
                                  scored_tree_t scored_tree(i, tree, score);
                                  scored_trees[i] = scored_tree;
                               }
                               return scored_trees;
                            };

   // we have made copies of trees - and allowed each to grow.
   // Now we need to remove the trees that are copies of better scoring trees. Score by sum of peptide spin score - allowing
   // backwards scoring also
   // 
   auto filter_trees = [] (std::vector<scored_tree_t> &scored_trees,
                           const std::map<unsigned int, std::set<unsigned int> > &tree_copies) {

                          for (auto &st : scored_trees) {
                             std::map<unsigned int, std::set<unsigned int> >::const_iterator it = tree_copies.find(st.index);
                             if (it != tree_copies.end()) {
                                const std::set<unsigned int> &set = it->second;
                                std::set<unsigned int>::const_iterator it_set;
                                for (it_set=set.begin(); it_set!=set.end(); ++it_set) {
                                   const unsigned int &idx_other_tree = *it_set;
                                   if (scored_trees[idx_other_tree].score > st.score) {
                                      st.marked_for_deletion = true;
                                      std::cout << "marking tree " << st.index << " for deleteion" << std::endl;
                                   }
                                }
                             } else {
                                // that tree had no copies
                             }
                          }
                          auto tree_eraser = [] (const scored_tree_t &st) { return st.marked_for_deletion; };
                          scored_trees.erase(std::remove_if(scored_trees.begin(), scored_trees.end(), tree_eraser), scored_trees.end());
                       };

   auto find_chains_that_overlap_other_chains = [] (mmdb::Manager *mol, float big_overlap_fraction_limit) {
                                                   mmdb::realtype local_dist_max = 3.0;
                                                   std::map<std::string, std::set<std::string> > delete_worse_chain_map;
                                                   std::map<std::string, std::set<std::string> > mergeable_chains_map;
                                                   std::map<std::string, std::map<std::string, unsigned int> > contact_counts_for_chain_pair;

                                                   std::map<std::string, unsigned int> ca_count_per_chain;

                                                   int SelHnd = mol->NewSelection();
                                                   mol->SelectAtoms(SelHnd, 0, "*",
                                                                    mmdb::ANY_RES, "*",
                                                                    mmdb::ANY_RES, "*",
                                                                    "*", " CA ", " C", "");

                                                   mmdb::Atom **atom_selection = 0;
                                                   int n_selected_atoms = 0;
                                                   mol->GetSelIndex(SelHnd, atom_selection, n_selected_atoms);

                                                   if (true) {
                                                      for (int iat=0; iat<n_selected_atoms; iat++) {
                                                         mmdb::Atom *at = atom_selection[iat];
                                                         std::cout << "selected-atoms " << iat << " " << coot::atom_spec_t(at) << " at " << coot::co(at).format() << std::endl;
                                                      }
                                                   }

                                                   mmdb::Model *model_p = mol->GetModel(1);
                                                   if (model_p) {
                                                      int n_chains = model_p->GetNumberOfChains();
                                                      for (int ichain=0; ichain<n_chains; ichain++) {
                                                         mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                                         std::string chain_id(chain_p->GetChainID());
                                                         int n_res = chain_p->GetNumberOfResidues();
                                                         ca_count_per_chain[chain_id] = n_res;
                                                      }
                                                   }
                                                   
                                                   mmdb::Contact *pscontact = NULL;
                                                   int n_contacts;
                                                   long i_contact_group = 1;
                                                   mmdb::mat44 my_matt;
                                                   for (int i=0; i<4; i++)
                                                      for (int j=0; j<4; j++)
                                                         my_matt[i][j] = 0.0;
                                                   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
                                                   //
                                                   mol->SeekContacts(atom_selection, n_selected_atoms,
                                                                     atom_selection, n_selected_atoms,
                                                                     0.0, local_dist_max, // min, max distances
                                                                     1,        // seqDist 0 -> in same res also
                                                                     pscontact, n_contacts,
                                                                     0, &my_matt, i_contact_group);

                                                   std::cout << "debug:: Selection found " << n_selected_atoms << " CAs" << std::endl;
                                                   if (n_contacts > 0) {
                                                      std::cout << "debug:: Selection had " << n_contacts << " CA-CA contacts" << std::endl;
                                                      if (pscontact) {
                                                         for (int i=0; i<n_contacts; i++) {
                                                            mmdb::Atom *at_1 = atom_selection[pscontact[i].id1];
                                                            mmdb::Atom *at_2 = atom_selection[pscontact[i].id2];
                                                            if (at_1->GetChain() == at_2->GetChain()) {
                                                            } else {
                                                               std::string chain_id_1(at_1->GetChainID());
                                                               std::string chain_id_2(at_2->GetChainID());
                                                               contact_counts_for_chain_pair[chain_id_1][chain_id_2]++;
                                                               std::cout << "debug:: contact_counts_for_chain_pair " << chain_id_1 << " " << chain_id_2
                                                                         << " : " << contact_counts_for_chain_pair[chain_id_1][chain_id_2] << std::endl;
                                                            }
                                                         }
                                                      }
                                                   }

                                                   if (false) { // debug
                                                      for(const auto &chain_id_1 : contact_counts_for_chain_pair) {
                                                         const std::string &key = chain_id_1.first;
                                                         const std::map<std::string, unsigned int> &c_map = chain_id_1.second;
                                                         for(const auto &chain_id_2 : c_map) {
                                                            std::cout << "contact-count: " << key << " "
                                                                      << chain_id_2.first << " " << chain_id_2.second << std::endl;
                                                         }
                                                      }
                                                   }

                                                   for(const auto &item_1 : contact_counts_for_chain_pair) {
                                                      const std::string &chain_id_1 = item_1.first;
                                                      const std::map<std::string, unsigned int> &c_map = item_1.second;
                                                      for(const auto &item_2 : c_map) {
                                                         const std::string &chain_id_2 = item_2.first;
                                                         unsigned int count = item_2.second;
                                                         float fraction_overlap_1 = static_cast<float>(count)/static_cast<float>(ca_count_per_chain[chain_id_1]);
                                                         float fraction_overlap_2 = static_cast<float>(count)/static_cast<float>(ca_count_per_chain[chain_id_2]);
                                                         std::cout << "contact-count: " << chain_id_1 << " " << chain_id_2 << " " << count << " "
                                                                   << fraction_overlap_1 << " from " << ca_count_per_chain[chain_id_1] << " CAs " 
                                                                   << fraction_overlap_2 << " from " << ca_count_per_chain[chain_id_2] << " CAs " << std::endl;
                                                         if (fraction_overlap_1 > big_overlap_fraction_limit)
                                                            delete_worse_chain_map[chain_id_1].insert(chain_id_2);
                                                      }
                                                   }
                                                   return delete_worse_chain_map;
                                                };

   auto make_scores_for_chain_ids = [] (const std::vector<scored_tree_t> &scored_trees) {
                                       std::map<std::string, double> scores_for_chain_ids;
                                       std::cout << "debug:: in make_scores_for_chain_ids() scored_trees size: " << scored_trees.size() << std::endl;
                                       for (unsigned int i=0; i<scored_trees.size(); i++) {
                                          const auto &scored_tree = scored_trees[i];
                                          std::cout << "tree: " << i << "  chain_id " << scored_tree.chain_id << " score " << scored_tree.score << std::endl;
                                          scores_for_chain_ids[scored_tree.chain_id] = scored_tree.score;
                                       }
                                       return scores_for_chain_ids;
                                    };

   // get rid of lower scoring chains that overlay better scoring chains
   //
   // Hmmm...   but what if 70% of the chains are overlapped? Then, we'd want to merge, not delete. Hmm.
   //
   auto filter_chains = [] (mmdb::Manager *mol,
                            const std::map<std::string, std::set<std::string> > &delete_worse_chain_map,
                            const std::map<std::string, double> &scores_for_chain_ids) {
                           std::set<std::string> delete_these_chains;
                           std::map<std::string, std::set<std::string> >::const_iterator it;
                           for (it=delete_worse_chain_map.begin(); it!=delete_worse_chain_map.end(); ++it) {
                              const std::string &chain_id_1 = it->first;
                              double score_1 = scores_for_chain_ids.at(chain_id_1); // a bit naughty?
                              const std::set<std::string> &chain_id_set = it->second;
                              std::set<std::string>::const_iterator it_set;
                              for (it_set=chain_id_set.begin(); it_set!=chain_id_set.end(); ++it_set) {
                                 const std::string &chain_id_2 = *it_set;
                                 double score_2 = scores_for_chain_ids.at(chain_id_2); // a bit naughty again?
                                 if (score_1 > score_2) delete_these_chains.insert(chain_id_2);
                              }
                           }

                           // now delete the chains in delete_these_chains
                           while (! delete_these_chains.empty()) {
                              const std::string &current_chain_id = *delete_these_chains.begin();
                              delete_these_chains.erase(delete_these_chains.begin());
                              mmdb::Model *model_p = mol->GetModel(1);
                              if (model_p) {
                                 std::cout << "DeleteChain for " << current_chain_id << std::endl;
                                 model_p->DeleteChain(current_chain_id.c_str());
                              }
                           }
                           mol->FinishStructEdit();
                        };

   // -------------------------------------------------------------
   // -------------------------------------------------------------
   //  main line
   // -------------------------------------------------------------
   // -------------------------------------------------------------

   unsigned int top_n_trees = 104;
   std::cout << "finding connections" << std::endl;

   std::pair<std::vector<tree_t>, std::map<unsigned int, std::set<unsigned int> >  > trees_pair = grow_trees(scored_pairs); // make and grow

   std::vector<std::deque<std::pair<unsigned int, coot::scored_node_t> > > &trees = trees_pair.first;
   std::map<unsigned int, std::set<unsigned int> > &tree_copies = trees_pair.second;

   if (true)
      print_the_tree_copies(trees, tree_copies);

   std::sort(trees.begin(), trees.end(), tree_length_sorter);
   unsigned int n_trees = trees.size();
   std::cout << "pre-tree merge: n_trees: " << n_trees << std::endl;
   if (true)
      print_trees(trees);

   std::vector<scored_tree_t> scored_trees = make_scored_trees(trees);

   std::cout << ":::::::::::::::::::::::::::: before filter_trees scored_tree.size() " << scored_trees.size()  << std::endl;
   filter_trees(scored_trees, tree_copies);
   std::cout << "::::::::::::::::::::::::::::   post filter_trees scored_tree.size() " << scored_trees.size()  << std::endl;

   // can the end of one tree be added to another tree?
   // (keep going until no more merging is possible)
   //
   // join_the_trees(scored_trees); // edit trees

   n_trees = scored_trees.size();
   std::cout << "post-tree merge: n_trees: " << n_trees << std::endl;

   std::sort(scored_trees.begin(), scored_trees.end(), scored_tree_sorter);

   add_chain_ids_to_scored_trees(scored_trees);

   n_trees = scored_trees.size();
   std::cout << "post-tree sort: n_trees: " << n_trees << std::endl;

   for (const auto &item : scored_trees) {
      std::cout << "Sorted-Scored-Tree " << item.index << " " << item.score << std::endl;
   }

   coot::minimol::molecule m;
   for (unsigned int i=0; i<scored_trees.size(); i++) {
      if (i > top_n_trees) continue;
      try {
         const auto &scored_tree = scored_trees[i];
         coot::minimol::fragment mc_fragment = tree_to_mc(scored_tree, i);
         delete_singleton_Ns(mc_fragment);
         invent_Ns_and_CBs_by_spin_search(mc_fragment, xmap);
         invent_extra_Ns_if_needed(mc_fragment);
         add_CBs_to_residues(mc_fragment);

         coot::minimol::molecule fragmol(mc_fragment);
         std::string fn = "mc-fragment-" + std::to_string(i) + ".pdb";
         fragmol.write_file(fn, 30.0);

         unsigned int frag_idx = m.fragment_for_chain(mc_fragment.fragment_id);
         m[frag_idx] = mc_fragment;
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }

   std::string fn = "all-poly-ALA-fragments.pdb";
   m.write_file(fn, 30.0);

   if (true) // debug the trees, write them out so that I can see them (c.f. the real model)
      debug_trees(trees);

   if (true) {
      for (unsigned int i=0; i<trees.size(); i++) {
         const auto &tree = trees[i];
         compare_tree_scores(tree, i);
      }
   }

   mmdb::Manager *mol = m.pcmmdbmanager();
   float big_overlap_fraction_limit = 0.3;
   std::map<std::string, std::set<std::string> > overlapping_chains_map =
      find_chains_that_overlap_other_chains(mol, big_overlap_fraction_limit);
   std::map<std::string, double> scores_for_chain_ids = make_scores_for_chain_ids(scored_trees);

   if (true) {
      std::cout << "debug::scores for chain_ids" << std::endl;
      std::map<std::string, double>::const_iterator it;
      for(it=scores_for_chain_ids.begin(); it!=scores_for_chain_ids.end(); ++it)
         std::cout << "scores-for-chain_ids " << it->first << " " << it->second << std::endl;
   }

   filter_chains(mol, overlapping_chains_map, scores_for_chain_ids);

   if (false) {
      for(const auto &item : overlapping_chains_map) {
         std::cout << "overlapping chains " << item.first << " : ";
         const std::set<std::string> &set = item.second;
         std::set<std::string>::const_iterator it;
         for (it=set.begin(); it!=set.end(); ++it) {
            std::cout << *it << " ";
         }
         std::cout << "\n";
      }
   }

   return mol;

}

void
globularize(mmdb::Manager *mol, const clipper::Xmap<float> &xmap) {

   // find the centre
   clipper::Coord_orth sum(0,0,0);

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_atoms = 0;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            int n_atoms_in_res = residue_p->GetNumberOfAtoms();
            for (int iat=0; iat<n_atoms_in_res; iat++) {
               mmdb::Atom *at = residue_p->GetAtom(iat);
               if (! at->isTer()) {
                  sum += coot::co(at);
                  n_atoms++;
               }
            }
         }
      }

      if (n_atoms > 0) {
         double f = 1.0/static_cast<double>(n_atoms);
         clipper::Coord_orth mean(f * sum); // fine for cryo, not good for crystallography

         // hack for testing tutorial modern
         mean = clipper::Coord_orth(60, 5, 12);

         const clipper::Spacegroup spgr = xmap.spacegroup();
         const clipper::Cell       cell = xmap.cell();

         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               int n_atoms_in_res = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms_in_res; iat++) {
                  mmdb::Atom *at = residue_p->GetAtom(iat);
                  if (! at->isTer()) {
                     clipper::Coord_orth current_position_xyz = coot::co(at);
                     clipper::Coord_orth current_position_frac(current_position_xyz.coord_frac(cell));
                     double d_best_sq = (current_position_xyz-mean).lengthsq();
                     clipper::Coord_orth pos_best = current_position_xyz;
                     bool updated = false;
                     for(int i_tr_x=-3; i_tr_x<=3; i_tr_x++) {
                        for(int i_tr_y=-3; i_tr_y<=3; i_tr_y++) {
                           for(int i_tr_z=-3; i_tr_z<=3; i_tr_z++) {
                              for (int i_symm=0; i_symm<spgr.num_symops(); i_symm++) {
                                 clipper::Coord_frac t(i_tr_x, i_tr_y, i_tr_z);
                                 clipper::Coord_frac cpt(current_position_frac + t);
                                 clipper::Coord_frac cfs(spgr.symop(i_symm) * cpt);
                                 clipper::Coord_orth p(cfs.coord_orth(cell));
                                 double d_this_sq = (p-mean).lengthsq();
                                 if (d_this_sq < d_best_sq) {
                                    d_best_sq = d_this_sq;
                                    pos_best = p;
                                    updated = true;
                                 }
                              }
                           }
                        }
                     }
                     if (updated) {
                        at->x = pos_best.x();
                        at->y = pos_best.y();
                        at->z = pos_best.z();
                     }
                  }
               }
            }
         }
      }
   }
}


std::vector<std::pair<unsigned int, coot::scored_node_t> >
make_some_fake_scored_pairs() {

   std::vector<std::pair<unsigned int, coot::scored_node_t> > v;

   coot::scored_node_t node_1(1, 1.0, 0.0, false, "node_1");
   coot::scored_node_t node_2(2, 1.0, 0.0, false, "node_2");
   coot::scored_node_t node_3(3, 1.0, 0.0, false, "node_3");
   coot::scored_node_t node_4(4, 1.0, 0.0, false, "node_4");
   coot::scored_node_t node_5(5, 1.0, 0.0, false, "node_5");
   coot::scored_node_t node_6(6, 1.0, 0.0, false, "node_6");
   coot::scored_node_t node_7(7, 1.0, 0.0, false, "node_7");
   coot::scored_node_t node_8(8, 1.0, 0.0, false, "node_8");

   // coot::scored_node_t node_9(3, 1.0, 0.0, false, "branch-node");
   // coot::scored_node_t node_10(8, 1.0, 0.0, false, "connected-to-branch-node");

   v.push_back(std::make_pair(0, node_1));
   v.push_back(std::make_pair(1, node_2));
   v.push_back(std::make_pair(2, node_3));
   v.push_back(std::make_pair(3, node_4));
   v.push_back(std::make_pair(4, node_5));
   v.push_back(std::make_pair(5, node_6));
   v.push_back(std::make_pair(6, node_7));
   v.push_back(std::make_pair(7, node_8));

   v.push_back(std::make_pair(9, coot::scored_node_t(10, 1.0, 0.0, false, "node_branch")));
   v.push_back(std::make_pair(10, coot::scored_node_t(15, 1.0, 0.0, false, "node_branch")));


   if (false) { // should make new trees (without mid-add)
      v.clear();
      v.push_back(std::make_pair(7, node_8));
      v.push_back(std::make_pair(6, node_7));
      v.push_back(std::make_pair(5, node_6));
      v.push_back(std::make_pair(4, node_5));
      v.push_back(std::make_pair(3, node_4));
      v.push_back(std::make_pair(2, node_3));
      v.push_back(std::make_pair(1, node_2));
      v.push_back(std::make_pair(0, node_1));
   }

   if (false) { // should add to front
      v.clear();

      // to add to the front:
      // the front_node.first should == this_scored_node.second.atom_idx

      v.push_back(std::make_pair(1, coot::scored_node_t(0, 1.0, 0.0, false, "node_1")));
      v.push_back(std::make_pair(2, coot::scored_node_t(1, 1.0, 0.0, false, "node_2")));
      v.push_back(std::make_pair(3, coot::scored_node_t(2, 1.0, 0.0, false, "node_3")));
      v.push_back(std::make_pair(4, coot::scored_node_t(3, 1.0, 0.0, false, "node_3")));
      v.push_back(std::make_pair(5, coot::scored_node_t(4, 1.0, 0.0, false, "node_3")));
      v.push_back(std::make_pair(6, coot::scored_node_t(5, 1.0, 0.0, false, "node_3")));

      v.push_back(std::make_pair(7, coot::scored_node_t(3, 1.0, 0.0, false, "node_branch")));
      
   }

   return v;
}

#include "backrub-rotamer.hh"

mmdb::Manager *
find_connected_fragments(const coot::minimol::molecule &flood_mol,
                         const clipper::Xmap<float> &xmap) {

   double variation = 0.5; // make bigger at lower resolutions (maybe up to 0.5?) pass this

   // somewhere here - not sure before or after action_mol is created, I want to globularize the molecule
   mmdb::Manager *action_mol = flood_mol.pcmmdbmanager();
   globularize(action_mol, xmap); // move around the atoms so they they are arranged in space in a sphere rather
                                  // than in a strip of the map (the asymmetric unit).

   action_mol->WritePDBASCII("flood-mol-globularized.pdb");
   
   mmdb::Atom **atom_selection = 0; // member data - cleared on destruction
   int n_selected_atoms = 0;
   int selhnd = action_mol->NewSelection();
   action_mol->SelectAtoms(selhnd, 0, "*",
                           mmdb::ANY_RES, // starting resno, an int
                           "*", // any insertion code
                           mmdb::ANY_RES, // ending resno
                           "*", // ending insertion code
                           "*", // any residue name
                           "*", // atom name
                           "*", // elements
                           "");
   action_mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
   std::cout << "INFO:: selected " << n_selected_atoms << " for distance pair check" << std::endl;
   std::vector<std::pair<unsigned int, unsigned int> > apwd =
      atom_pairs_within_distance(action_mol, atom_selection, n_selected_atoms, 3.81, variation);
   std::cout << "spin_score_pairs..." << std::endl;
   // the first of the scores is the index of the first atom
   std::vector<std::pair<unsigned int, coot::scored_node_t> > scores =
      spin_score_pairs(apwd, xmap, action_mol, atom_selection, n_selected_atoms);
   std::cout << "spin_score_pairs done" << std::endl;

   // scores = make_some_fake_scored_pairs();
   mmdb::Manager *mol = make_fragments(scores, atom_selection, xmap);
   action_mol->DeleteSelection(selhnd);
   return mol;
}


// mutate mol
void
apply_sequence_to_fragments(mmdb::Manager *mol, const clipper::Xmap<float> &xmap, const coot::fasta_multi &fam, const coot::protein_geometry &pg) {

   unsigned int n_sequences = fam.size();
   std::cout << "debug:: apply_sequence_to_fragments(): with n_sequences " << n_sequences << std::endl;
   if (n_sequences > 0) {
      for (unsigned int idx=0; idx<n_sequences; idx++) {
         std::string sequence = fam[idx].sequence;
         std::cout << "debug sequence: " << sequence << std::endl;
         const std::string &name = fam[idx].name;

         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               if (n_res > 5) {
                  std::string chain_id(chain_p->GetChainID());
                  int resno_start = chain_p->GetResidue(0)->GetSeqNum();
                  int resno_end   = chain_p->GetResidue(n_res-1)->GetSeqNum();

                  coot::side_chain_densities scd;
                  std::pair<std::string, std::vector<mmdb::Residue *> > a_run_of_residues =
                     scd.setup_test_sequence(mol, chain_id, resno_start, resno_end, xmap);
                  std::cout << "debug:: a run of residues: " << a_run_of_residues.first << " with " << a_run_of_residues.second.size() << " residues" << std::endl;
                  if (a_run_of_residues.first.empty()) {
                     scd.test_sequence(a_run_of_residues.second, xmap, name, sequence);
                  } else {
                     std::cout << "ERROR:: when generating a run of residues " << std::endl;
                     std::cout << a_run_of_residues.first << std::endl;
                  }
                  coot::side_chain_densities::results_t new_sequence_result = scd.get_result();
                  std::string new_sequence = new_sequence_result.sequence;
                  std::cout << "debug:: new_sequence " << new_sequence << std::endl;
                  if (! new_sequence.empty()) {
                     int sl = new_sequence.length();
                     if (sl == n_res) {
                        for (int ires=0; ires<n_res; ires++) {
                           mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                           if (residue_p) {
                              char letter = new_sequence[ires];
                              std::string new_residue_type = coot::util::single_letter_to_3_letter_code(letter);
                              coot::util::mutate(residue_p, new_residue_type);
                           }
                        }
                     }
                  }
               }
            }
            coot::backrub_molecule(mol, &xmap, pg);
         }
      }
   }
}

void proc(const clipper::Xmap<float> &xmap, const coot::fasta_multi &fam, float rmsd_cut_off) {

   coot::minimol::molecule flood_molecule = get_flood_molecule(xmap, rmsd_cut_off);
   mmdb::Manager *mol = find_connected_fragments(flood_molecule, xmap);
   coot::protein_geometry pg;
   pg.init_standard();
   apply_sequence_to_fragments(mol, xmap, fam, pg); // mutate mol
   mol->WritePDBASCII("sequenced.pdb");

}

#include <clipper/ccp4/ccp4_map_io.h>

int main(int argc, char **argv) {
   int status = 0;

   std::string hklin_file_name = "rnasa-1.8-all_refmac1.mtz";
   std::string f_col_label = "FWT";
   std::string phi_col_label = "PHWT";

   // hklin_file_name = "../src/tm-A-1.0.mtz";
   // f_col_label   = "FC";
   // phi_col_label = "PHIC";

   std::string pir_file_name = "../src/rnase.pir";
   coot::fasta_multi fam;
   fam.read(pir_file_name);

   bool test_from_mtz = true;

   if (test_from_mtz) {
      if (! f_col_label.empty()) {
         if (! phi_col_label.empty()) {
            try {
               std::cout << "Read mtz file " << hklin_file_name << " " << f_col_label << " " << phi_col_label << std::endl;
               clipper::Xmap<float> xmap;
               bool use_weights = false;
               bool is_diff_map = false;
               bool stat = coot::util::map_fill_from_mtz(&xmap, hklin_file_name, f_col_label, phi_col_label, "",
                                                         use_weights, is_diff_map);
               if (stat) {
                  proc(xmap, fam, 1.5);
               }
            }
            catch (clipper::Message_fatal &rte) {
               std::cout << "ERROR::" << rte.text() << std::endl;
            }
         }
      }
   } else {
      try {
         std::string map_file_name = "emd_22898.map";
         if (coot::file_exists(map_file_name)) {
            clipper::CCP4MAPfile file;
            file.open_read(map_file_name);
            clipper::Xmap<float> xmap;
            file.import_xmap(xmap);
            proc(xmap, fam, 5.5);
         } else {
            std::cout << "No such file " << map_file_name << std::endl;
         }
      }
      catch (clipper::Message_fatal &rte) {
         std::cout << "ERROR::" << rte.text() << std::endl;
      }
   }

   return status;
}
//
