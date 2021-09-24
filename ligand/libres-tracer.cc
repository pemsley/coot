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
#include "coot-utils/merge-atom-selections.hh"
#include "coot-utils/merge-C-and-N-terminii.hh"
#include "scored-node.hh"
#include "ligand.hh"
#include "utils/coot-utils.hh"
#include "side-chain-densities.hh"
#include "analysis/stats.hh"
#include "ideal/simple-restraint.hh"
#include "ideal/regularize-minimol.hh"
#include "rama-rsr-extend-fragments.hh"
#include "libres-tracer.hh"

typedef std::deque<std::pair<unsigned int, coot::scored_node_t> > tree_t;

class scored_tree_t {
public:
   unsigned int index; // index in the vector of trees (same as vector of scored_trees)
   std::string chain_id;
   tree_t tree;
   double forward_score;
   double backward_score;
   bool marked_for_deletion;
   bool checked;
   std::set<unsigned int> live_progenitor_index_set;
   scored_tree_t(unsigned int idx, const std::string &chain_id_in, const tree_t &tree_in, const std::set<unsigned int> &lpis, double sf, double sb) :
      index(idx), chain_id(chain_id_in), tree(tree_in), forward_score(sf), backward_score(sb), live_progenitor_index_set(lpis)  {
      marked_for_deletion = false;
      checked = false; }
   scored_tree_t() { index = 0; marked_for_deletion = false; checked = false; forward_score = 0.0; backward_score = 0.0; }
};

coot::minimol::molecule
get_flood_molecule(const clipper::Xmap<float> &xmap, float rmsd_cut_off, float flood_atom_mask_radius) {

   bool debug = true;
   coot::ligand lig;

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


float spin_score_best_so_far = 0.0;

std::pair<unsigned int, coot::scored_node_t>
spin_score(unsigned int idx_1, unsigned int idx_2, mmdb::Atom **atom_selection,
           const clipper::Xmap<float> &xmap, float map_rmsd) {

   auto index_to_pos = [atom_selection] (unsigned int idx) {
                          mmdb::Atom *at = atom_selection[idx];
                          return clipper::Coord_orth(at->x, at->y, at->z);
                       };

   float inv_rmsd = 1.0f/map_rmsd;

   // so that density values beyond ~3 are scored less well and beyond 4 are ~0.1
   // This should stop disulfides being marked as the best peptide due to their
   // strong density. This punishes strong density.
   //
   // the passed value is direct rho from the map
   //
   auto roll_off = [inv_rmsd] (const float rho) {
                      float f = inv_rmsd * rho;
                      // return f + f * f * f * f * f * -0.003f; // was 0.002
                      // return f + f * f * f * f * f * -0.0003f; // was 0.002
                      // std::cout << "debug f: " << f << std::endl;
                      return f;
                   };

#if 0 // original scales
   float scale_CO       =  1.95;
   float scale_CO_low   = -0.6;
   float scale_CO_anti  = -0.15;
   float scale_N        =  1.2;
   float scale_perp     = -0.8;
   float scale_mid      =  1.6;  // not used
   float scale_non_line =  0.4;  // not used
#endif

   float scale_CO       =  1.4;
   float scale_CO_low   = -0.8;
   float scale_CO_anti  = -0.3;
   float scale_N        =  1.0;
   float scale_N_low    = -1.0;
   float scale_perp     = -0.9;
   float scale_mid      =  1.6;  // not used
   float scale_non_line =  0.4;  // not used

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
   clipper::Coord_orth rel_line_pt_O_low(  diff_p_unit * f_ca_ca_o + perp_unit * 3.2); // was 3.7
   clipper::Coord_orth rel_line_pt_CO_anti(diff_p_unit * f_ca_ca_o * 0.9 - perp_unit * 0.6);
   clipper::Coord_orth rel_line_pt_N(      diff_p_unit * f_ca_ca_n - perp_unit * 0.3);
   clipper::Coord_orth rel_line_pt_N_low(  diff_p_unit * f_ca_ca_n - perp_unit * 1.45);
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

      clipper::Coord_orth p_N_low = coot::util::rotate_around_vector(diff_p_unit,
                                                                     pos_1 + rel_line_pt_N_low,
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
      float rho_N_low   = coot::util::density_at_point(xmap, p_N_low);
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
         scale_CO      * roll_off(rho_CO)      +
         scale_CO_low  * roll_off(rho_CO_low)  +
         scale_CO_anti * roll_off(rho_CO_anti) +
         scale_N       * roll_off(rho_N)       +
         scale_N_low   * roll_off(rho_N_low)   +
         scale_perp    * roll_off(rho_perp_1)  +
         scale_perp    * roll_off(rho_perp_2);
      // scale_N_accpt * rho_acceptor_best

      if (false) { // write out finger-print points
         std::string fn = "fingerprint-point-" + std::to_string(idx_1) + "-" + std::to_string(idx_2) + ".table";
         std::ofstream f(fn.c_str(), std::ios_base::app);
         // std::ofstream f("fingerprint-point");
         if (f) {
            f << "idx_1 " << idx_1 << " " << idx_2 << "\n";
            f << "CO      " << p_CO.x()      << " " << p_CO.y()      << " " << p_CO.z()      << "\n";
            f << "CO_low  " << p_CO_low.x()  << " " << p_CO_low.y()  << " " << p_CO_low.z()  << "\n";
            f << "CO_anti " << p_CO_anti.x() << " " << p_CO_anti.y() << " " << p_CO_anti.z() << "\n";
            f << "N       " << p_N.x()       << " " << p_N.y()       << " " << p_N.z()       << "\n";
            f << "N_low   " << p_N_low.x()   << " " << p_N_low.y()   << " " << p_N_low.z()   << "\n";
            f << "perp-1  " << p_2.x()       << " " << p_2.y()       << " " << p_2.z()       << "\n";
            f << "perp-2  " << p_3.x()       << " " << p_3.y()       << " " << p_3.z()       << "\n";
            f << "pos-1   " << pos_1.x()     << " " << pos_1.y()     << " " << pos_1.z()     << "\n";
            f << "pos-2   " << pos_2.x()     << " " << pos_2.y()     << " " << pos_2.z()     << "\n";
            f.close();
         }
      }

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

   // why am I changing the density value after the scoring?

   float non_line_equal_density_penalty_1 = roll_off(rho_at_1 + rho_at_2 - 2 * rho_mid);
   float non_line_equal_density_penalty =
      - non_line_equal_density_penalty_1 * non_line_equal_density_penalty_1/map_rmsd;

   best_score += scale_mid * roll_off(rho_mid);
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
make_spin_scored_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &atom_pairs_within_distance,
                       unsigned int n_top,
                       const clipper::Xmap<float> &xmap,
                       mmdb::Manager *mol, mmdb::Atom **atom_selection, int n_selected_atoms) {

   // apwd : atom (index) pairs within distance

   std::vector<std::pair<unsigned int, coot::scored_node_t> > scores;
   if (! mol) return scores;

   bool debug = false;

   unsigned int n_atom_pairs = atom_pairs_within_distance.size();
   scores.resize(n_atom_pairs * 2); // results go here

   std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
   float map_rmsd = sqrt(mv.second);

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
                                       auto spin_score_1 = spin_score(at_idx_1, at_idx_2, atom_selection, xmap, map_rmsd);
                                       auto spin_score_2 = spin_score(at_idx_2, at_idx_1, atom_selection, xmap, map_rmsd);
                                       scores[2*i  ] = spin_score_1;
                                       scores[2*i+1] = spin_score_2;
                                       scores[2*i  ].second.reverse_spin_score = std::make_pair(true, scores[2*i+1].second.spin_score);
                                       scores[2*i+1].second.reverse_spin_score = std::make_pair(true, scores[2*i  ].second.spin_score);
                                    }
                                 };

#if 1
   for (unsigned int i=0; i<n_atom_pairs; i++) {
      const unsigned &at_idx_1 = atom_pairs_within_distance[i].first;
      const unsigned &at_idx_2 = atom_pairs_within_distance[i].second;
      scores[2*i ]  = spin_score(at_idx_1, at_idx_2, atom_selection, xmap, map_rmsd);
      scores[2*i+1] = spin_score(at_idx_2, at_idx_1, atom_selection, xmap, map_rmsd);
   }
#endif

#if 0
   std::cout << "debug:: n_atom_pairs: " << n_atom_pairs << std::endl;
   std::vector<std::thread> threads;
   for (unsigned int i=0; i<air.size(); i++) {
      const auto &range = air[i];
      std::cout << "thread " << i << " has range " << range.first << " " << range.second << std::endl;
      threads.push_back(std::thread(spin_score_workpackage, range, std::cref(atom_pairs_within_distance), atom_selection,
                                    std::cref(xmap), map_rmsd, std::ref(scores)));
   }
   for (unsigned int i=0; i<air.size(); i++)
      threads[i].join();
#endif

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
   // std::cout << "---" << std::endl;

   // for (const auto &item : tree)
   // std::cout << "  Post-Tree " << item.first << " " << item.second.atom_idx << " " << item.second.name << std::endl;

   // std::cout << "::::::: remove_tree_back done, now size " << tree.size() << std::endl;
}



// get rid of less chains
//
// Hmmm...   but what if 70% of the chains are overlapped? Then, we'd want to merge, not delete. Hmm.
// Well, that is not a concern for this function.
//
void
filter_similar_chains(mmdb::Manager *mol,
                      const std::map<std::string, std::set<std::string> > &delete_worse_chain_map) {

   std::set<std::string> delete_these_chains;
   std::map<std::string, std::set<std::string> >::const_iterator it;
   for (it=delete_worse_chain_map.begin(); it!=delete_worse_chain_map.end(); ++it) {
      const std::set<std::string> &chain_id_set = it->second;
      std::set<std::string>::const_iterator it_set;
      for (it_set=chain_id_set.begin(); it_set!=chain_id_set.end(); ++it_set) {
         const std::string &chain_id_2 = *it_set;
         delete_these_chains.insert(chain_id_2);
      }
   }

   // now delete the chains in delete_these_chains
   mmdb::Model *model_p = mol->GetModel(1);
   if (model_p) {
      unsigned int count = 0; // for output formatting
      while (! delete_these_chains.empty()) {
         const std::string &current_chain_id = *delete_these_chains.begin();
         delete_these_chains.erase(delete_these_chains.begin());
         if (count == 0)
            std::cout << "filter_similar_chains(): DeleteChain";
         std::cout << " " << current_chain_id;
         model_p->DeleteChain(current_chain_id.c_str());
         count++;
         if (count == 30) {
            std::cout << "\n";
            count = 0;
         }
      }
      if (count != 0) std::cout << "\n";
   }
   mol->FinishStructEdit();
}


#include <cmath>

// return a map of sets of chain_ids to be deleted.
std::map<std::string, std::set<std::string> >
find_chains_that_overlap_other_chains(mmdb::Manager *mol, float big_overlap_fraction_limit,
                                      const std::map<std::string, double> &scores_for_chains) {

   auto choose_deletable_chain_from_pair = [] (unsigned int overlap_count,
                                               unsigned int ca_count_for_chain_1,
                                               unsigned int ca_count_for_chain_2,
                                               const std::string &chain_id_1, const std::string &chain_id_2,
                                               const std::map<std::string, double> &scores_for_chains) {

                                              // chain_id_2 is a similar to chain_id_1 but lower scoring, so it should be deleted.
                                              // return std::make_pair(chain_id_1,  chain_id_2); i.e. chain_id_2 is a deletable version of ~chain_id_1

                                              if (overlap_count == (ca_count_for_chain_2 - 1))
                                                 if (overlap_count < (ca_count_for_chain_1 - 1))
                                                    return std::make_pair(chain_id_1,  chain_id_2);

                                              if (overlap_count == (ca_count_for_chain_1 - 1))
                                                 if (overlap_count < (ca_count_for_chain_2 - 1))
                                                    return std::make_pair(chain_id_2,  chain_id_1);

                                              const double &score_1 = scores_for_chains.at(chain_id_1);
                                              const double &score_2 = scores_for_chains.at(chain_id_2);
                                              if (score_2 < score_1) {
                                                 return std::make_pair(chain_id_1,  chain_id_2);
                                              } else {
                                                 return std::make_pair(chain_id_2,  chain_id_1);
                                              }
                                           };

   // keep a record of the residue number differences
   class res_no_delta_stats_t {
   public:
      int sum;
      int sum_sq;
      int count;
      res_no_delta_stats_t() : sum(0), sum_sq(0), count(0) {}
      void add(mmdb::Atom *at_1, mmdb::Atom *at_2) {
         // It is usually the case that the residue numbers are the same
         int res_no_1 = at_1->residue->GetSeqNum();
         int res_no_2 = at_2->residue->GetSeqNum();
         int d = abs(res_no_2-res_no_1);
         // std::cout << "debug:: add() " << res_no_1 << " " << res_no_2 << " d " << d << std::endl;
         sum += d;
         sum_sq += d * d;
         count++;
      }
      float average() const { if (count == 0) return -999.9f; else return static_cast<float>(sum)/static_cast<float>(count); }
      float variance() const { float m = average(); float ms = static_cast<float>(sum_sq)/static_cast<float>(count); return ms - m*m; }
      bool chains_go_in_the_same_direction() const {
         if (count == 0) return true;
         float v = variance();
         float s = std::sqrt(v);
         // float a = average();
         // std::cout << "debug:: res_no_delta_stats_t: count " << count << " mean " << a << " v " << v << " s " << s << std::endl;
         return (s < 3.0f);
      }
   };

   mmdb::realtype local_dist_max = 1.0;
   std::map<std::string, std::set<std::string> > delete_worse_chain_map;
   std::map<std::string, std::set<std::string> > mergeable_chains_map;
   // std::map<std::string, std::map<std::string, unsigned int> > contact_counts_for_chain_pair;
   std::map<std::string, std::map<std::string, res_no_delta_stats_t> > contact_counts_for_chain_pair;

   std::map<std::string, unsigned int> ca_count_per_chain;

   int SelHnd = mol->NewSelection(); // d
   mol->SelectAtoms(SelHnd, 0, "*",
                    mmdb::ANY_RES, "*",
                    mmdb::ANY_RES, "*",
                    "*", " CA ", " C", "");

   mmdb::Atom **atom_selection = 0;
   int n_selected_atoms = 0;
   mol->GetSelIndex(SelHnd, atom_selection, n_selected_atoms);

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
   auto tp_0 = std::chrono::high_resolution_clock::now();
   mol->SeekContacts(atom_selection, n_selected_atoms,
                     atom_selection, n_selected_atoms,
                     0.0, local_dist_max, // min, max distances
                     1,        // seqDist 0 -> in same res also
                     pscontact, n_contacts,
                     0, &my_matt, i_contact_group);

   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   std::cout << "Timings: for SeekContacts() in find_chains_that_overlap_other_chains(): " << d10 << " milliseconds" << std::endl;

   std::cout << "debug:: Selection found " << n_selected_atoms << " CAs" << std::endl;
   if (n_contacts > 0) {
      // std::cout << "debug:: Selection had " << n_contacts << " CA-CA contacts" << std::endl;
      if (pscontact) {
         std::cout << "debug:: Selection found " << n_contacts << " contacts" << std::endl;
         for (int i=0; i<n_contacts; i++) {
            mmdb::Atom *at_1 = atom_selection[pscontact[i].id1];
            mmdb::Atom *at_2 = atom_selection[pscontact[i].id2];
            if (at_1->GetChain() == at_2->GetChain()) {
            } else {
               std::string chain_id_1(at_1->GetChainID());
               std::string chain_id_2(at_2->GetChainID());
               contact_counts_for_chain_pair[chain_id_1][chain_id_2].add(at_1, at_2);

               if (n_contacts%1000==0)
                  std::cout << "debug:: contact_counts_for_chain_pair " << chain_id_1 << " " << chain_id_2
                            << " : " << contact_counts_for_chain_pair[chain_id_1][chain_id_2].sum
                            << contact_counts_for_chain_pair[chain_id_1][chain_id_2].count << " "
                            << contact_counts_for_chain_pair[chain_id_1][chain_id_2].average()
                            << std::endl;
            }
         }
      }
   }

   if (false) { // debug
      for (const auto &chain_id_1 : contact_counts_for_chain_pair) {
         const std::string &key = chain_id_1.first;
         const std::map<std::string, res_no_delta_stats_t> &s_map = chain_id_1.second;
         for(const auto &stats : s_map) {

            if (false)
               std::cout << "contact-count: " << key << " "
                         << stats.first << " " << stats.second.sum << " " << stats.second.count << std::endl;
         }
      }
   }

   auto tp_2 = std::chrono::high_resolution_clock::now();

   for (const auto &item_1 : contact_counts_for_chain_pair) {
      const std::string &chain_id_1 = item_1.first;
      const std::map<std::string, res_no_delta_stats_t> &s_map = item_1.second;
      for (const auto &item_2 : s_map) {
         const std::string &chain_id_2 = item_2.first;
         const res_no_delta_stats_t &rnds = item_2.second;
         float fraction_overlap_1 = static_cast<float>(rnds.count)/static_cast<float>(ca_count_per_chain[chain_id_1]);
         float fraction_overlap_2 = static_cast<float>(rnds.count)/static_cast<float>(ca_count_per_chain[chain_id_2]);
         bool something_is_deletable = false; // when true, one of these is deletable

         if (fraction_overlap_1 > 0.5 || fraction_overlap_2 > 0.5) {

            bool same_direction = rnds.chains_go_in_the_same_direction();
            if (same_direction) {

               // chains go in same direction

               if (fraction_overlap_1 > big_overlap_fraction_limit || fraction_overlap_2 > big_overlap_fraction_limit)
                  something_is_deletable = true;

               if (rnds.count == (static_cast<int>(ca_count_per_chain[chain_id_1]) - 1))
                  something_is_deletable = true;

               if (rnds.count == (static_cast<int>(ca_count_per_chain[chain_id_2]) - 1))
                  something_is_deletable = true;

               if (false)
                  std::cout << "contact-count: " << chain_id_1 << " " << chain_id_2 << " " << rnds.count << " "
                            << fraction_overlap_1 << " from " << ca_count_per_chain[chain_id_1] << " CAs " 
                            << fraction_overlap_2 << " from " << ca_count_per_chain[chain_id_2] << " CAs "
                            << "something_is_deletable " << something_is_deletable << std::endl;

               if (something_is_deletable) {
                  // the second is a deletable chain copy of the key/first
                  auto dc = choose_deletable_chain_from_pair(rnds.count,
                                                             ca_count_per_chain[chain_id_1],
                                                             ca_count_per_chain[chain_id_2],
                                                             chain_id_1, chain_id_2, scores_for_chains);
                  delete_worse_chain_map[dc.first].insert(dc.second);
                  if (false) {
                     if (dc.first.length() == 1) {
                        std::cout << "contact-count: " << chain_id_1 << " " << chain_id_2 << " " << rnds.count << " "
                                  << fraction_overlap_1 << " from " << ca_count_per_chain[chain_id_1] << " CAs " 
                                  << fraction_overlap_2 << " from " << ca_count_per_chain[chain_id_2] << " CAs "
                                  << " makring " << dc.second << " for deletion,"
                                  << " scores " << scores_for_chains.at(chain_id_1) << " " << scores_for_chains.at(chain_id_2)
                                  << std::endl;
                     }
                  }
               }
            } else {

               // delete the shorter one.
               // if (ca_count_per_chain[chain_id_1] > ca_count_per_chain[chain_id_2])
               // delete_worse_chain_map[chain_id_1].insert(chain_id_2);
               // else
               // delete_worse_chain_map[chain_id_2].insert(chain_id_1);

            }
         }
      }
   }
   auto tp_3 = std::chrono::high_resolution_clock::now();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   std::cout << "Timings: find_chains_that_overlap_other_chains() post-processing: part 1: " << d21 << " milliseconds" << std::endl;
   std::cout << "Timings: find_chains_that_overlap_other_chains() post-processing: part 2: " << d32 << " milliseconds" << std::endl;
   mol->DeleteSelection(SelHnd);
   return delete_worse_chain_map;
};


void
find_chains_that_overlap_other_chains_inner(mmdb::Manager *mol, float big_overlap_fraction_limit,
                                            const std::vector<std::string> &chain_ids,
                                            std::pair<unsigned int, unsigned int> start_start,
                                            std::map<std::string, std::set<std::string> > &results) { // fill results, natch

   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   //

   mmdb::Atom **atom_selection = 0; // member data - cleared on destruction
   int n_selected_atoms = 0;
   int selhnd = mol->NewSelection(); // d
   mol->SelectAtoms(selhnd, 0, "*",
                           mmdb::ANY_RES, // starting resno, an int
                           "*", // any insertion code
                           mmdb::ANY_RES, // ending resno
                           "*", // ending insertion code
                           "*", // any residue name
                           "*", // atom name
                           "*", // elements
                           "");
   mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

   mol->SeekContacts(atom_selection, n_selected_atoms,
                     atom_selection, n_selected_atoms,
                     0.0f, 1.0f,
                     1,        // seqDist 0 -> in same res also
                     pscontact, n_contacts,
                     0, &my_matt, i_contact_group);

   if (n_contacts > 0) {
      if (pscontact) {
         for (int i=0; i<n_contacts; i++) {
         }
      }
   }
   mol->DeleteSelection(selhnd);
}



std::map<std::string, std::set<std::string> >
find_chains_that_overlap_other_chains_v2(mmdb::Manager *mol, float big_overlap_fraction_limit,
                                         const std::map<std::string, double> &scores_for_chains) {

   std::map<std::string, std::set<std::string> > delete_worse_chain_map;

   // first what are the chain ids, how many are there - so that we can split them up.
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      std::vector<std::string> chain_ids;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         chain_ids.push_back(std::string(chain_p->GetChainID()));
      }

      std::cout << "There were " << chain_ids.size() << " chain ids" << std::endl;

      std::vector<std::pair<unsigned int, unsigned int> > ranges = coot::atom_index_ranges(chain_ids.size(), 32);
      std::vector<std::map<std::string, std::set<std::string> > > results(ranges.size());

   }
   return delete_worse_chain_map;
}



mmdb::Manager *
make_fragments(std::vector<std::pair<unsigned int, coot::scored_node_t> > &scored_pairs,
               mmdb::Atom **atom_selection,
               const clipper::Xmap<float> &xmap,
               unsigned int top_n_fragments,
               bool cryo_em_peptide_seeds_mode) { // pass the atom_selection for debugging

   // maybe it will be faster to calculate which residues are connected to which other residues:
   //
   // std::map<unsigned int, std::vector<std::pair<unsigned int, coot::scored_node_t> > > neighbour_map
   // that can be calculated in parallel (with locking)
   //

   // can the end of one tree be added to another tree?
   // (keep going until no more merging is possible)
   //
   // edit trees

   auto distance_check_min = [atom_selection] (int atom_index_1, int atom_index_2, double dist_min) {
                                mmdb::Atom *at_1 = atom_selection[atom_index_1];
                                mmdb::Atom *at_2 = atom_selection[atom_index_2];
                                clipper::Coord_orth pt_1 = coot::co(at_1);
                                clipper::Coord_orth pt_2 = coot::co(at_2);
                                double dist_min_sqrd = dist_min * dist_min;
                                double dist_sqrd = (pt_1-pt_2).lengthsq();
                                if (false)
                                   std::cout << "debug distance_check_min: " << atom_index_1 << " " << atom_index_2 << " dist: " << sqrt(dist_sqrd)
                                             << std::endl;
                                return dist_sqrd > dist_min_sqrd;
                         };



   auto node_is_already_in_tree = [] (const std::pair<unsigned int, coot::scored_node_t> &scored_pair, const tree_t &tree) {
                                     bool found = false;
                                     for(const auto &item : tree) {
                                        if (item.first == scored_pair.first) {
                                           if (item.second.atom_idx == scored_pair.second.atom_idx) {
                                              found = true;
                                              break;
                                           }
                                        }
                                     }
                                     return found;
                                  };

   auto node_is_similar_to_another_node = [atom_selection] (const std::pair<unsigned int, coot::scored_node_t> &scored_pair, const tree_t &tree) {
                                             bool similar = false;
                                             double dist_close = 3.0;
                                             const double dist_close_sqrd = dist_close * dist_close;
                                             clipper::Coord_orth p_1_sc = coot::co(atom_selection[scored_pair.first]);
                                             clipper::Coord_orth p_2_sc = coot::co(atom_selection[scored_pair.second.atom_idx]);
                                             for (const auto &item : tree) {
                                                clipper::Coord_orth p_1_item = coot::co(atom_selection[item.first]);
                                                clipper::Coord_orth p_2_item = coot::co(atom_selection[item.second.atom_idx]);

                                                if (false) {
                                                   double d1 = std::sqrt((p_1_item-p_1_sc).lengthsq());
                                                   double d2 = std::sqrt((p_2_item-p_2_sc).lengthsq());
                                                   std::cout << "debug d1 d2 " << d1 << " " << d2 << std::endl;
                                                }

                                                if ((p_1_item-p_1_sc).lengthsq() < dist_close_sqrd) {
                                                   if ((p_2_item-p_2_sc).lengthsq() < dist_close_sqrd) {
                                                      similar = true;
                                                      break;
                                                   }
                                                }
                                                if ((p_1_item-p_2_sc).lengthsq() < dist_close_sqrd) {
                                                   if ((p_2_item-p_1_sc).lengthsq() < dist_close_sqrd) {
                                                      similar = true;
                                                      break;
                                                   }
                                                }
                                             }
                                             return similar;
                                          };

   std::atomic<bool> tree_addition_lock(false);

   auto get_tree_addition_lock = [&tree_addition_lock] () {
                            bool unlocked = false;
                            while (! tree_addition_lock.compare_exchange_weak(unlocked, true)) {
                               std::this_thread::sleep_for(std::chrono::nanoseconds(100));
                               unlocked = false;
                            }
                         };

   auto release_tree_addition_lock = [&tree_addition_lock] () {
                                        tree_addition_lock = false;
                                     };

   // same again for the "dont_add" protection
   std::atomic<bool> dont_add_lock(false);

   auto get_dont_add_lock = [&dont_add_lock] () {
                            bool unlocked = false;
                            while (! dont_add_lock.compare_exchange_weak(unlocked, true)) {
                               std::this_thread::sleep_for(std::chrono::nanoseconds(10)); // or faster?  (was 100ns)
                               unlocked = false;
                            }
                         };

   auto release_dont_add_lock = [&dont_add_lock] () {
                                   dont_add_lock = false;
                                };
   class addition_t {
   public:
      enum end_t { FRONT, BACK };
      std::size_t tree_index;
      std::size_t scored_pair_index;
      end_t end;
      addition_t(const std::size_t &tree_index, const std::size_t &scored_pair_index, end_t end_in) :
         tree_index(tree_index), scored_pair_index(scored_pair_index), end(end_in) {}
   };

   auto get_additions_for_this_range_of_trees = [distance_check_min, node_is_already_in_tree, node_is_similar_to_another_node, get_dont_add_lock, release_dont_add_lock]
      (const std::vector<scored_tree_t> &trees,
       std::pair<unsigned int, unsigned int> start_stop,
       const std::vector<std::pair<unsigned int, coot::scored_node_t> > &scored_pairs,
       std::map<std::size_t, std::set<std::size_t> > &dont_adds, // by ref
       std::vector<addition_t> &additions) { // fill the results

                                                    for (unsigned int j=start_stop.first; j<start_stop.second; j++) {
                                                       const auto &scored_tree = trees[j];
                                                       if (scored_tree.marked_for_deletion) continue;

                                                       const auto &tree_node_b = scored_tree.tree.back();
                                                       const auto &tree_node_f = scored_tree.tree.front();

                                                       for (unsigned int i=0; i<scored_pairs.size(); i++) {

                                                          // std::cout << "considering scored-pair " << i << " of " << scored_pairs.size() << std::endl;
                                                          const auto &scored_pair = scored_pairs[i];

                                                          if (false) {
                                                             std::cout << "raw-compare: "
                                                                       << std::setw(4) << scored_pair.first << " " << std::setw(4) << scored_pair.second.atom_idx << " to tree-node-f "
                                                                       << std::setw(4) << tree_node_f.first << " " << std::setw(4) << tree_node_f.second.atom_idx << std::endl;
                                                             std::cout << "raw-compare: "
                                                                       << std::setw(4) << scored_pair.first << " " << std::setw(4) << scored_pair.second.atom_idx << " to tree-node-b "
                                                                       << std::setw(4) << tree_node_b.first << " " << std::setw(4) << tree_node_b.second.atom_idx << std::endl;
                                                          }

                                                          if (tree_node_b.second.atom_idx == scored_pair.first) {
                                                             if (tree_node_b.first == scored_pair.second.atom_idx) {
                                                                // skip the trivial reverse
                                                             } else {

#if 0
                                                                // I had written: // maybe this access needs locking?
                                                                // And after a week of testing, I discovered that it hung here in one case. So Add the locking.
                                                                //
                                                                get_dont_add_lock();
                                                                bool dont_add = (dont_adds[j].find(i) == dont_adds[j].end());
                                                                release_dont_add_lock();

#endif
                                                                // but maybe I can get away with not needing the lock here?
                                                                auto it_end_for_dont_add = dont_adds[j].end();
                                                                bool dont_add_flag = (dont_adds[j].find(i) == it_end_for_dont_add);

                                                                if (dont_add_flag) { // OK to add then!
                                                                   // CA(n-1) to CA(n+1) is shortest for helicies (5.3)
                                                                   if (distance_check_min(tree_node_b.first, scored_pair.second.atom_idx, 5.0)) { // 5.3 with some room for noise
                                                                      if (! node_is_already_in_tree(scored_pair, scored_tree.tree)) {
                                                                         if (! node_is_similar_to_another_node(scored_pair, scored_tree.tree)) {
                                                                            addition_t addition(j, i, addition_t::BACK);
                                                                            if (false)
                                                                               std::cout << "adding to tree back with tree-index: " << j << " this-tree size " << scored_tree.tree.size()
                                                                                         << " tree_node_f indices " << tree_node_f.first << " " << tree_node_f.second.atom_idx
                                                                                         << " scored-pair indices : " << scored_pair.first << " " << scored_pair.second.atom_idx
                                                                                         << " spin-score: " << scored_pair.second.spin_score << std::endl;

                                                                            additions.push_back(addition);
                                                                            get_dont_add_lock();
                                                                            dont_adds[j].insert(i);
                                                                            release_dont_add_lock();
                                                                         }
                                                                      }
                                                                   }
                                                                }
                                                             }
                                                          }
                                                       }
                                                    }
                                                 };

   auto filter_out_old_traces_from_progenitors = [] (std::vector<scored_tree_t> &traces) {

                                                    // I need only check the traces that I haven't checked before. Hmm

                                                    for (unsigned int i=0; i<traces.size(); i++) {
                                                       auto &trace = traces[i];
                                                       if (trace.marked_for_deletion) continue;
                                                       // if (trace.checked) continue; // hmm - I think I should be able to do this, but if I enable this line, I get different results!
                                                       trace.checked = true; // it is now
                                                       unsigned int this_trace_length = trace.tree.size();
                                                       // std::cout << "debug:: i " << i << " this_trace_length: " << this_trace_length
                                                       // << " trace.live_progenitor_index_set() " << trace.live_progenitor_index_set.size() << std::endl;
                                                       std::set<unsigned int>::const_iterator it_set;
                                                       for (it_set=trace.live_progenitor_index_set.begin(); it_set!=trace.live_progenitor_index_set.end(); ++it_set) {
                                                          const unsigned int &progenitor_index = *it_set;
                                                          unsigned int progenitor_trace_length = traces[progenitor_index].tree.size();
                                                          if (false)
                                                             std::cout << "debug:: trace " << i << " comparing this_trace_length " << this_trace_length
                                                                       << " with progenitor_trace_length " << progenitor_trace_length
                                                                       << " delta " << this_trace_length-progenitor_trace_length << std::endl;
                                                          if (! traces[progenitor_index].marked_for_deletion) {
                                                             if (this_trace_length > (progenitor_trace_length + 4)) {
                                                                traces[progenitor_index].marked_for_deletion = true;
                                                                if (false)
                                                                   std::cout << "DEBUG:: trace " << i << " filter_out_old_traces_from_progenitors() marking trace "
                                                                          << progenitor_index << " for deletion" << std::endl;
                                                             }
                                                          }
                                                       }
                                                    }
                                                 };


   // this can be abstracted (with a bit of edting)
   //
   auto trace_length_histogram = [] (const std::vector<scored_tree_t> &sorted_traces) {

                                    auto stars = [] (unsigned int n) {
                                                    unsigned int n_stars = static_cast<int>(static_cast<float>(n)/30.0);
                                                    for (unsigned int i=0; i<n_stars; i++)
                                                       std::cout << "*";
                                                    };
                                    unsigned int trace_length_max = 0;
                                    for (const auto &trace : sorted_traces) {
                                       unsigned int l = trace.tree.size();
                                       if (l > trace_length_max) {
                                          trace_length_max = l;
                                       }
                                    }
                                    // so now we have set trace_length_max
                                    const unsigned int n_bins = 20;
                                    std::vector<unsigned int> length_bins(n_bins, 0);
                                    for (const auto &trace : sorted_traces) {
                                       unsigned int l = trace.tree.size();
                                       float f  = static_cast<float>(l)/static_cast<float>(trace_length_max);
                                       float fb = static_cast<unsigned int>(f * static_cast<float>(n_bins));
                                       unsigned int bin_index = static_cast<unsigned int>(fb);
                                       if (bin_index == 20) bin_index = 19;
                                       length_bins[bin_index]++;
                                    }
                                    std::cout << ":::: Trace Length Histogram:\n";
                                    for (unsigned int i=0; i<n_bins; i++) {
                                       float abl = (static_cast<float>(i) + 0.5) * static_cast<float>(trace_length_max) / static_cast<float>(n_bins);
                                       int abl_i = static_cast<int>(abl);
                                       std::cout << std::setw(2) << abl_i << " : " << std::setw(5) << length_bins[i] << " ";
                                       stars(length_bins[i]);
                                       std::cout << "\n";
                                    }
                                 };

   auto index_to_chain_id = [] (unsigned int  i) {
                               std::string s;

                               if (i < 26) {
                                  s = std::string(1, 'A' + i);
                               } else {
                                  if (i < 27 * 26) {
                                     unsigned int idx_1 = i/26;
                                     unsigned int idx_2 = i - 26 * idx_1;
                                     std::string s1 = std::string(1, 'A' + idx_1 - 1);
                                     std::string s2 = std::string(1, 'A' + idx_2);
                                     s = s1 + s2;
                                  } else {
                                     // 702 or more
                                     unsigned int idx_1 = i/(27*26);
                                     unsigned int r = i - idx_1 * (27*26);
                                     unsigned int idx_2 = r/26;
                                     std::string s1 = std::string(1, 'A' + idx_1 - 1);
                                     std::string s2 = std::string(1, 'A' + idx_2 - 1);
                                     if (idx_2 == 0) {
                                        s2 = "_";
                                        unsigned int idx_3 = r;
                                        std::string s3 = std::string(1, 'A' + idx_3);
                                        s = s1 + s2 + s3;
                                     } else {
                                        unsigned int idx_3 = r - idx_2 * 26;
                                        std::string s3 = std::string(1, 'A' + idx_3);
                                        s = s1 + s2 + s3;
                                     }
                                  }
                               }
                               return s;
                            };

   auto tree_report = [trace_length_histogram] (const std::vector<scored_tree_t> &trees) {

                         // put this in a detached thread?
                         // unsigned int n_trees = trees.size();
                         // (trees can be changing size as we do this) - so be a bit careful
                         // std::vector<scored_tree_t> trees_copy(n_trees);
                         // for (unsigned int i=0; i<n_trees; i++)
                         //     trees_copy[i] = trees[i];
                         // trace_length_histogram(trees_copy);
                         // std::thread t(trace_length_histogram, std::cref(trees_copy));
                         // t.detach();
                         // no. crashy.

                         trace_length_histogram(trees);
                   };

   // make *sorted* peptide seeds
   auto make_peptide_seeds = [atom_selection, index_to_chain_id] (const std::vector<std::pair<unsigned int, coot::scored_node_t> > &scored_pairs_in,
                                                                  unsigned int max_number_of_seeds) {

                                auto seed_sorter = [] (const scored_tree_t &st1, const scored_tree_t &st2) {
                                                      return st2.forward_score < st1.forward_score;
                                                };

                                auto index_to_pos = [atom_selection] (unsigned int idx) {
                                                       mmdb::Atom *at = atom_selection[idx];
                                                       return clipper::Coord_orth(at->x, at->y, at->z);
                                                    };

                                auto distance_check_min = [atom_selection] (int atom_index_1, int atom_index_2, double dist_min) {
                                                             mmdb::Atom *at_1 = atom_selection[atom_index_1];
                                                             mmdb::Atom *at_2 = atom_selection[atom_index_2];
                                                             clipper::Coord_orth pt_1 = coot::co(at_1);
                                                             clipper::Coord_orth pt_2 = coot::co(at_2);
                                                             double dist_min_sqrd = dist_min * dist_min;
                                                             double dist_sqrd = (pt_1-pt_2).lengthsq();
                                                             if (false)
                                                                std::cout << "debug distance_check_min: " << atom_index_1 << " " << atom_index_2 << " dist: " << sqrt(dist_sqrd)
                                                                          << std::endl;
                                                             return dist_sqrd > dist_min_sqrd;
                                                          };

                                std::cout << "PROGRESS:: making peptide seeds" << std::endl;

                                // a seed is a pair of scored_pairs in a scored_tree_t;
                                // sorted by goodness and the return vector is truncated so that it is not insanely long
                                //
                                std::vector<scored_tree_t> seeds; // return this

                                for (unsigned int i=0; i<scored_pairs_in.size(); i++) {
                                   const std::pair<unsigned int, coot::scored_node_t> &base = scored_pairs_in[i];
                                   for (unsigned int j=0; j<scored_pairs_in.size(); j++) {
                                      if (i != j) {
                                         const std::pair<unsigned int, coot::scored_node_t> &jth_spin_pair = scored_pairs_in[j];
                                         // does the jth spin pair connect to the base?
                                         if (jth_spin_pair.first == base.second.atom_idx) {
                                            // avoid the trivial reverse peptide
                                            if (jth_spin_pair.second.atom_idx != base.first) {
                                               if (distance_check_min(jth_spin_pair.second.atom_idx, base.first, 5.0)) {
                                                  double score = base.second.spin_score + jth_spin_pair.second.spin_score;
                                                  scored_tree_t st;
                                                  st.tree.push_back(base);
                                                  st.tree.push_back(jth_spin_pair);
                                                  st.forward_score = score;
                                                  seeds.push_back(st);
                                               }
                                            }
                                         }
                                      }
                                   }
                                }

                                std::sort(seeds.begin(), seeds.end(), seed_sorter);

                                if (false) {
                                   if (seeds.size() > max_number_of_seeds) {
                                      int every_nth = 1 + seeds.size()/max_number_of_seeds;
                                      std::vector<scored_tree_t> sampled_seeds;
                                      for (unsigned int i=0; i<seeds.size(); i += every_nth) {
                                         sampled_seeds.push_back(seeds[i]);
                                      }
                                      seeds = sampled_seeds;
                                   }
                                }

                                if (false)
                                   for (unsigned int i=0; i<seeds.size(); i++)
                                      std::cout << "sorted seed " << i <<  " "
                                                << seeds[i].tree[0].first << " "
                                             << seeds[i].tree[0].second.atom_idx << " "
                                                << seeds[i].tree[1].first << " "
                                                << seeds[i].tree[1].second.atom_idx << " "
                                                << seeds[i].forward_score << "\n";

                                if (seeds.size() > max_number_of_seeds)
                                   seeds.resize(max_number_of_seeds);

                                for (unsigned int i=0; i<seeds.size(); i++)
                                   seeds[i].chain_id = index_to_chain_id(i);

                                if (true) {
                                   std::ofstream f("debug-seed-points");
                                   if (f) {
                                      for (unsigned int i=0; i<seeds.size(); i++) {
                                         clipper::Coord_orth pt_1 = index_to_pos(seeds[i].tree[0].first);
                                         clipper::Coord_orth pt_2 = index_to_pos(seeds[i].tree[0].second.atom_idx);
                                         clipper::Coord_orth pt_3 = index_to_pos(seeds[i].tree[1].first);
                                         clipper::Coord_orth pt_4 = index_to_pos(seeds[i].tree[1].second.atom_idx);
                                         f << "seed " << i << " "
                                           << pt_1.x() << " " << pt_1.y() << " " << pt_1.z() << " "
                                           << pt_2.x() << " " << pt_2.y() << " " << pt_2.z() << " "
                                           << pt_3.x() << " " << pt_3.y() << " " << pt_3.z() << " "
                                           << pt_4.x() << " " << pt_4.y() << " " << pt_4.z() << " "
                                           << seeds[i].forward_score << "\n";
                                      }
                                   } else {
                                      std::cout << "failed to open debug-seed-points" << std::endl;
                                   }
                                }

                                std::cout << "PROGRESS:: done peptide seeds" << std::endl;

                                return seeds;
                             };

   auto grow_trees_v4 = [get_additions_for_this_range_of_trees, tree_report,
                         filter_out_old_traces_from_progenitors] (const std::vector<std::pair<unsigned int, coot::scored_node_t> > &scored_pairs_in) {

                           auto there_are_additions_to_be_made = [] (const std::vector<std::vector<addition_t> > &additions_vec) {
                                                                    bool status = false;
                                                                    for (auto it = additions_vec.begin(); it != additions_vec.end(); ++it) {
                                                                       if (! it->empty()) {
                                                                          status = true;
                                                                          break;
                                                                       }
                                                                    }
                                                                    return status;
                                                                 };

                           unsigned int max_number_of_trees = 14000; // so that we don't keep going for millions of trees

                           std::cout << "grow_trees_v4" << std::endl;

                           auto scored_pairs = scored_pairs_in;

                           std::vector<scored_tree_t> trees; // return this

                           // return this also.
                           // The vast majority of trees will be sub-trees of the final trees. We could
                           // detect those and delete them before returning. Should be cheap to do.
                           std::map<unsigned int, std::set<unsigned int> > tree_copies;

                           // don't add peptides to the original/unextented tree the next round after a peptide was added
                           std::map<std::size_t, std::set<std::size_t> > dont_adds;

                           // for scored_pair, can I add this scored_pair to each/any of the current set of trees?
                           // if I can, make a note.
                           //
                           // If I can't then make a tree and add it to the tree list

                           std::size_t loop_begin_tree_size = trees.size();
                           std::size_t scored_pairs_front_offset = 0;
                           unsigned int n_threads = 16;

                           do {
                              loop_begin_tree_size = trees.size();

                              std::vector<std::vector<addition_t> > additions_vec;

                              auto tp_0 = std::chrono::high_resolution_clock::now();

                              // make a non-thread version of this for testing/timing
                              //
                              std::vector<std::pair<unsigned int, unsigned int> > ranges = coot::atom_index_ranges(loop_begin_tree_size, n_threads);
                              additions_vec.resize(ranges.size());
                              std::vector<std::thread> threads;
                              for (unsigned int ir=0; ir<ranges.size(); ir++)
                                 threads.push_back(std::thread(get_additions_for_this_range_of_trees, std::cref(trees), ranges[ir], std::cref(scored_pairs),
                                                               std::ref(dont_adds), std::ref(additions_vec[ir])));
                              for (unsigned int ir=0; ir<ranges.size(); ir++)
                                 threads[ir].join();

#if 1
                              auto tp_1 = std::chrono::high_resolution_clock::now();
                              auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
                              // std::cout << "Timings: for these " << trees.size() << " trees: " << d10 << " microseconds" << std::endl;
#endif
                              if (there_are_additions_to_be_made(additions_vec)) {

                                 for (auto it = additions_vec.begin(); it != additions_vec.end(); ++it) {
                                    auto &additions(*it);
                                    for (const auto &addition : additions) {
                                       scored_tree_t new_tree = trees[addition.tree_index];
                                       const auto &scored_pair = scored_pairs[addition.scored_pair_index];
                                       if (addition.end == addition_t::FRONT)
                                          new_tree.tree.push_front(scored_pair);
                                       else
                                          new_tree.tree.push_back(scored_pair);
                                       unsigned int new_index = trees.size();
                                       std::set<unsigned int> &progenitor_list = trees[addition.tree_index].live_progenitor_index_set; // and then we add to it.
                                       progenitor_list.insert(addition.tree_index);
                                       new_tree.index = new_index;
                                       trees.push_back(new_tree);
                                       tree_copies[addition.tree_index].insert(new_index);
                                    }
                                 }

                                 // I need to add here - I think - something to purge the early trees if later trees are n (4?) or more
                                 // scored pairs longer - beacuse to be 4 or more pairs longer means that we clearly went the right path
                                 // and we don't need to keep the stubby root tree.

                                 filter_out_old_traces_from_progenitors(trees); // edit trees

                              } else {

                                 if (! scored_pairs.empty()) {
                                    // make a new tree with this peptide.
                                    if (scored_pairs_front_offset < scored_pairs.size()) {
                                       std::vector<std::pair<unsigned int, coot::scored_node_t> >::iterator it = scored_pairs.begin() + scored_pairs_front_offset;
                                       // scored_pairs.erase(scored_pairs.begin());
                                       scored_tree_t new_tree;
                                       unsigned int new_index = trees.size();
                                       new_tree.index = new_index;
                                       new_tree.tree.push_back(*it);
                                       scored_pairs_front_offset++;
                                       trees.push_back(new_tree);
                                    }
                                 }
                              }

                              if (true)
                                 if ((trees.size() %100) == 0) {
                                    std::cout << "Trees size " << trees.size() << std::endl;
                                    tree_report(trees);
                                 }

                           } while ((trees.size() > loop_begin_tree_size) && (trees.size() < max_number_of_trees) );

                           // now erase the earlier copies (which we since extended), starting from the end.
                           // all the keys need to be deleted. Let's make a vector of the tree indices to be deleted.
                           //
                           // Using a set would have beeen more concise.
                           //

                           // debugging counts;
                           unsigned int n_delete_from_copies = 0;
                           unsigned int n_delete_from_short_peptides = 0;
                           unsigned int n_delete_from_marked_for_deletion = 0;

                           std::vector<unsigned int> delete_these_indices;

                           // add to the list-of-deletions those trees that are 3 peptides or less
                           // (some of these are useful, but there is a lot of "noise")
                           for (unsigned int i=0; i<trees.size(); i++)
                              if (trees[i].tree.size() < 4)
                                 if (std::find(delete_these_indices.begin(), delete_these_indices.end(), i) == delete_these_indices.end()) {
                                    delete_these_indices.push_back(i);
                                    n_delete_from_short_peptides++;
                                 }

                           // it's a bad sign if all the deletions come from short peptides - so allow them, so that we have
                           // something to build on at least! We are not in a good state if this happens.
                           if (n_delete_from_short_peptides == trees.size())
                              delete_these_indices.clear();

                           std::map<unsigned int, std::set<unsigned int> >::const_iterator it_map;
                           for (it_map=tree_copies.begin(); it_map!=tree_copies.end(); ++it_map)
                              if (std::find(delete_these_indices.begin(), delete_these_indices.end(), it_map->first) == delete_these_indices.end()) {
                                 delete_these_indices.push_back(it_map->first);
                                 n_delete_from_copies++;
                              }

                           // add to the list-of-delete_these_indices those that are marked for deletion
                           for (unsigned int i=0; i<trees.size(); i++)
                              if (trees[i].marked_for_deletion) {
                                 if (std::find(delete_these_indices.begin(), delete_these_indices.end(), i) == delete_these_indices.end())
                                    delete_these_indices.push_back(i);
                                 n_delete_from_marked_for_deletion ++;
                              }

                           // sort delete_these_indices so that high values are first

                           if (true) {
                              std::cout << "DEBUG:: n_delete_from_copies: " << n_delete_from_copies
                                        << " n_delete_from_short_peptides " << n_delete_from_short_peptides
                                        << " n_delete_from_marked_for_deletion " << n_delete_from_marked_for_deletion
                                        << std::endl;
                           }

                           if (delete_these_indices.size() == trees.size()) {
                              trees.clear();
                              std::cout << "No traces remain" << std::endl;
                           } else {

                              auto index_sorter = [] (unsigned int i1, unsigned int i2) { return i2<i1; };
                              std::sort(delete_these_indices.begin(), delete_these_indices.end(), index_sorter);

                              // high indices are first in the list
                              std::cout << "INFO:: grow_trees_v4() Now to delete " << delete_these_indices.size() << " of the " << trees.size()
                                        << " trees" << std::endl;
                              for (unsigned int i=0; i<delete_these_indices.size(); i++) {
                                 if (false) // debugged now - I should have used a set!
                                    std::cout << "debug:: delete_these_indices: " << i << " -> index " << delete_these_indices[i]
                                              << " trees.size() " << trees.size() << std::endl;
                                 auto it = trees.begin() + delete_these_indices[i];
                                 trees.erase(it);
                              }
                              std::cout << "INFO:: " << trees.size() << " trees remain." << std::endl;
                           }
                           return trees;
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
                                                       coot::minimol::atom  N_at(" N  ", " N", new_positions.first,  "", 10.0);
                                                       coot::minimol::atom CB_at(" CB ", " C", new_positions.second, "", 10.0);
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
                                                coot::minimol::atom N_at(" N  ", " N", N_pos, "", 10.0);
                                                residue.addatom(N_at);
                                             }
                                          }
                                       }
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


   auto trace_to_main_chain = [atom_selection] (const scored_tree_t &scored_tree, unsigned int ith_trace) {

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
                           clipper::Coord_orth rel_line_pt_C(diff_p_unit * f_ca_ca_c + perp_unit * 0.48);  // was 0.48
                           clipper::Coord_orth rel_line_pt_N(diff_p_unit * f_ca_ca_n - perp_unit * 0.47);  // was 0.47
                           clipper::Coord_orth rel_line_pt_O(diff_p_unit * f_ca_ca_o + perp_unit * 1.89);
                           clipper::Coord_orth O_position_raw = pt_1 + rel_line_pt_O;
                           clipper::Coord_orth N_position_raw = pt_1 + rel_line_pt_N;
                           clipper::Coord_orth C_position_raw = pt_1 + rel_line_pt_C;
                           // now rotate O_position about diff_p
                           clipper::Coord_orth rotated_position_for_O = coot::util::rotate_around_vector(diff_p, O_position_raw, pt_1, pep.second.alpha);
                           clipper::Coord_orth rotated_position_for_N = coot::util::rotate_around_vector(diff_p, N_position_raw, pt_1, pep.second.alpha);
                           clipper::Coord_orth rotated_position_for_C = coot::util::rotate_around_vector(diff_p, C_position_raw, pt_1, pep.second.alpha);

                           coot::minimol::atom  O(" O  ", " O", rotated_position_for_O, "", 1.0, 10.0f);
                           coot::minimol::atom  N(" N  ", " N", rotated_position_for_N, "", 1.0, 10.0f);
                           coot::minimol::atom  C(" C  ", " C", rotated_position_for_C, "", 1.0, 10.0f);
                           coot::minimol::atom CA(" CA ", " C", pt_1, "", 10.0f);
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
                              r_f.addatom(O);
                              r_f.addatom(C);
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
                                                   coot::colour_holder c(0.2, 0.9, 0.2);
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
                               if (i > 2224) continue;
                               coot::colour_holder col = index_to_colour(i);
                               const auto &tree = trees[i];
                               std::cout << "------- Tree " << i << " size: " << tree.size() << " ------" << std::endl;
                               if (tree.size() < 1) {
                               } else {
                                  for (const auto &item : tree) {
                                     clipper::Coord_orth pt_1 = index_to_pos(item.first);
                                     clipper::Coord_orth pt_2 = index_to_pos(item.second.atom_idx);
                                     std::cout << "   Tree-b " << item.first << " " << item.second.atom_idx << " "
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
                                             coot::minimol::atom CB(" CB ", " C", cb_pos.second, "", 1.0, 10.0f);
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

   // return the reverse score also, and use that if it is better
   //
   auto get_tree_scores = [] (const scored_tree_t &tree) {
                             double sum_forward = 0.0;
                             double sum_backward = 0.0;
                             for (unsigned int i=0; i<tree.tree.size(); i++) {
                                const auto &pep = tree.tree[i];
                                sum_forward += pep.second.spin_score;
                                if (pep.second.reverse_spin_score.first)
                                   sum_backward += pep.second.reverse_spin_score.second;
                             }
                             return std::make_pair(sum_forward, sum_backward);
                          };

   // add a chain-id also
   auto make_scored_trees = [get_tree_scores, index_to_chain_id] (const std::vector<scored_tree_t> &unscored_trees) {
                               std::vector<scored_tree_t> scored_trees(unscored_trees.size());
                               for (unsigned int i=0; i<unscored_trees.size(); i++) {
                                  const scored_tree_t &tree = unscored_trees[i]; //  teehee!
                                  std::pair<double, double> score = get_tree_scores(tree);
                                  std::string chain_id = index_to_chain_id(i);
                                  if (false)
                                     std::cout << "tree with index " << std::setw(3) << i << " and chain-id " << chain_id
                                               << " has length " << tree.tree.size()
                                               << " and scores " << score.first << " " << score.second << std::endl;
                                  scored_tree_t scored_tree(i, chain_id, tree.tree, tree.live_progenitor_index_set, score.first, score.second);
                                  scored_trees[i] = scored_tree;
                               }
                               return scored_trees;
                            };

   auto make_scores_for_chain_ids = [] (const std::vector<scored_tree_t> &scored_trees) {
                                       std::map<std::string, double> scores_for_chain_ids;
                                       std::cout << "debug:: in make_scores_for_chain_ids() scored_trees size: " << scored_trees.size() << std::endl;
                                       for (unsigned int i=0; i<scored_trees.size(); i++) {
                                          const auto &scored_tree = scored_trees[i];
                                          if (scored_tree.marked_for_deletion) continue;
                                          if (false)
                                             std::cout << "tree: " << i << "  chain_id " << scored_tree.chain_id << " score " << scored_tree.forward_score
                                                       << " from " << scored_tree.tree.size() << " peptides " << std::endl;
                                          scores_for_chain_ids[scored_tree.chain_id] = scored_tree.forward_score;
                                       }
                                       return scores_for_chain_ids;
                                    };

   // return a minimol of mainchain+CB fragments
   auto scored_trees_to_multi_fragment = [trace_to_main_chain,
                                          delete_singleton_Ns,
                                          invent_Ns_and_CBs_by_spin_search,
                                          invent_extra_Ns_if_needed,
                                          add_CBs_to_residues] (const std::vector<scored_tree_t> &scored_trees,
                                             unsigned int top_n_traces,
                                             const clipper::Xmap<float> &xmap) {

                                            std::cout << "debug:: scored_trees_to_multi_fragment() was given " << scored_trees.size() << " scored traces"
                                                      << std::endl;

                                            coot::minimol::molecule m; // return this
                                            m.set_cell(xmap.cell());
                                            m.set_spacegroup(xmap.spacegroup().symbol_hm());
                                            unsigned int n_trees_added = 0;
                                            for (unsigned int i=0; i<scored_trees.size(); i++) {
                                               if (i > top_n_traces) continue;
                                               try {
                                                  const auto &scored_tree = scored_trees[i];
                                                  if (scored_tree.marked_for_deletion) {
                                                     std::cout << "debug:: scored_trees_to_multi_fragment() tree " << i << " marked for deletion" << std::endl;
                                                     continue;
                                                  }

                                                  if (false)
                                                     std::cout << "debug:: scored_trees_to_multi_fragment() trace " << i << " has length "
                                                               << scored_tree.tree.size() << std::endl;

                                                  if (scored_tree.tree.size() > 0) {
                                                     coot::minimol::fragment mc_fragment = trace_to_main_chain(scored_tree, i);

                                                     delete_singleton_Ns(mc_fragment);
                                                     invent_Ns_and_CBs_by_spin_search(mc_fragment, xmap);
                                                     invent_extra_Ns_if_needed(mc_fragment);
                                                     add_CBs_to_residues(mc_fragment);

                                                     if (false) { // debugging fragments
                                                        coot::minimol::molecule fragmol(mc_fragment);
                                                        std::string chain_id = mc_fragment.fragment_id;
                                                        std::string fn = "mc-initial-build-fragment-" + chain_id + ".pdb";
                                                        fragmol.write_file(fn, 30.0);
                                                        std::cout << "for chain " << i << " " << chain_id << " with forward_score " << scored_tree.forward_score
                                                                  << " made " << fn << std::endl;
                                                     }

                                                     unsigned int frag_idx = m.fragment_for_chain(mc_fragment.fragment_id);
                                                     m[frag_idx] = mc_fragment;
                                                     n_trees_added++;
                                                  }
                                               }
                                               catch (const std::runtime_error &rte) {
                                                  std::cout << "ERROR:: " << rte.what() << std::endl;
                                               }
                                            }
                                            std::cout << "debug:: scored_trees_to_multi_fragment() molecule m has " << n_trees_added << " fragments" << std::endl;
                                            return m;
                                         };

   auto debug_longest_tree = [atom_selection] (const std::vector<std::deque<std::pair<unsigned int, coot::scored_node_t> > > &trees) {

                                auto index_to_pos = [atom_selection] (int idx) {
                                                       mmdb::Atom *at = atom_selection[idx];
                                                       return clipper::Coord_orth(at->x, at->y, at->z);
                                                    };

                                unsigned int n_top_traces = 800;

                                // index and length
                                std::vector<std::pair<std::size_t, std::size_t> > tree_length_vec;
                                for (unsigned int i=0; i<trees.size(); i++) {
                                   const auto &tree = trees[i];
                                   tree_length_vec.push_back(std::make_pair(i, tree.size()));
                                }

                                auto tls = [] (const std::pair<std::size_t, std::size_t> &t1,
                                               const std::pair<std::size_t, std::size_t> &t2) { return t2.second < t1.second; };

                                std::sort(tree_length_vec.begin(), tree_length_vec.end(), tls);

                                // for (const auto &item : tree_length_vec)
                                // std::cout << "debug_longest_tree(): longest index: " << item.first << " length: " << item.second << std::endl;

                                std::ofstream f("longest-trace.points");
                                for (unsigned int j=0; j<tree_length_vec.size(); j++) {
                                   if (j >= n_top_traces) continue;
                                   const auto &item_il = tree_length_vec[j];
                                   std::size_t tree_index = item_il.first;
                                   const tree_t &tree = trees[tree_index];
                                   coot::colour_holder col(0.3 + static_cast<float>(j) * 0.1f, 0.4, 0.5 - static_cast<float>(j) * 0.1);
                                   for (unsigned int i=0; i<tree.size(); i++) {
                                      const auto &item = tree[i];
                                      clipper::Coord_orth pt_1 = index_to_pos(item.first);
                                      clipper::Coord_orth pt_2 = index_to_pos(item.second.atom_idx);
                                      // std::cout << "   debug_longest_tree(): Tree-a " << j << " " << item.first << " " << item.second.atom_idx << " "
                                      //           << item.second.name << std::endl;
                                      f << "longest-trace-point " << col.red << " " << col.green << " " << col.blue << " "
                                        << pt_1.x() << " " << pt_1.y() << " " << pt_1.z() << " "
                                        << pt_2.x() << " " << pt_2.y() << " " << pt_2.z() << "\n";
                                   }
                                }
                                f.close();
                             };

   auto write_pdbs_for_each_of_the_chains = [] (mmdb::Manager *mol) {
                                               coot::minimol::molecule m(mol);
                                               for (unsigned int i=0; i<m.fragments.size(); i++) {
                                                  coot::minimol::molecule mf(m.fragments[i]);
                                                  std::string fn_pdb = "mc-fragment-" + m.fragments[i].fragment_id + ".pdb";
                                                  std::string fn_cif = "mc-fragment-" + m.fragments[i].fragment_id + ".cif";
                                                  mf.write_file(fn_pdb, 30.0f);
                                                  mf.write_cif_file(fn_cif);
                                               }
                                            };

   // Lots of the traces still look alike. There are hundreds or even thousands of very similar traces
   //
   // I want to thin that down to about 10 or so for each trace cluster
   auto eigenfilter_traces = [&scored_pairs, atom_selection] (std::vector<scored_tree_t> &traces) {

                                class position_and_eigenvectors_and_value_t {
                                public:
                                   clipper::Coord_orth pos;
                                   std::vector<clipper::Coord_orth> eigenvectors;
                                   std::vector<double> eigenvalues;
                                };

                                std::vector<position_and_eigenvectors_and_value_t> position_and_eigenvectors_and_value(traces.size());
                                for (unsigned int i=0; i<traces.size(); i++) {
                                   std::vector<clipper::Coord_orth> points;
                                   const auto &trace(traces[i]);
                                   for (auto it_pep = trace.tree.begin(); it_pep != trace.tree.end(); ++it_pep) {
                                      const auto &peptide = *it_pep;
                                      const int &atom_idx = peptide.second.atom_idx;
                                      mmdb::Atom *at = atom_selection[atom_idx];
                                      clipper::Coord_orth pt = coot::co(at);
                                      points.push_back(pt);
                                   }

                                   // calculate mean and eigen vectors of that.
                                }
                             };

   // -------------------------------------------------------------
   // -------------------------------------------------------------
   //  main line of make_fragments()
   // -------------------------------------------------------------
   // -------------------------------------------------------------

   if (cryo_em_peptide_seeds_mode) {

      top_n_fragments = 400;
      unsigned int max_number_of_seeds = top_n_fragments;
      auto seeds = make_peptide_seeds(scored_pairs, max_number_of_seeds);

      std::cout << "debug:: cryo_em_peptide_seeds_mode: top_n_fragments is " << top_n_fragments << std::endl;
      coot::minimol::molecule m = scored_trees_to_multi_fragment(seeds, top_n_fragments, xmap);
      std::cout << "DEBUG:: scored_trees_to_multi_fragment() was given " << seeds.size() << " seeds and returned "
                << m.fragments.size() << " fragments and " << m.count_atoms() << " atoms" << std::endl;
      mmdb::Manager *mol = m.pcmmdbmanager();

      return mol;

   } else {

      auto tp_0 = std::chrono::high_resolution_clock::now();

      std::cout << "finding connections" << std::endl;
      std::cout << "PROGRESS calling grow_trees_v4()" << std::endl;
      std::cout << "Grow trees..." << std::endl;
      std::vector<scored_tree_t> trees = grow_trees_v4(scored_pairs); // make and grow, not scored yet!
      auto tp_1 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "Timings: from make_fragments() grow_trees_v4()                     " << d10 << " milliseconds" << std::endl;

      // debug_longest_tree(trees);

      std::cout << "INFO:: :::::: after grow_trees, there are " << trees.size() << " trees" << std::endl;

      if (trees.empty()) {
         return 0; // null pointer - nothing to build on!
      } else {

         // sort traces by length
         auto trace_sorter = [] (const scored_tree_t &tree_1, const scored_tree_t &tree_2) { return tree_2.tree.size() < tree_1.tree.size(); };
         std::sort(trees.begin(), trees.end(), trace_sorter); // sorted by length
         std::vector<scored_tree_t> &sorted_traces = trees; //so that it's clear they have been sorted!

         eigenfilter_traces(sorted_traces); // function not yet implemented. let's call them traces. 

         // debug_trees(sorted_trees);

         trace_length_histogram(sorted_traces);

         std::vector<scored_tree_t> scored_traces = make_scored_trees(sorted_traces); // add a chain-id also.

         if (true) {
            for (const auto &st : scored_traces) {
               double delta = st.forward_score - st.backward_score;
               std::cout << "INFO:: scored_trace " << st.chain_id << " has forward_score " << std::setw(5) << st.forward_score
                         << " backward_score " << st.backward_score << " delta " << delta << std::endl;
            }
         }

         auto tp_2 = std::chrono::high_resolution_clock::now();
         coot::minimol::molecule m = scored_trees_to_multi_fragment(scored_traces, top_n_fragments, xmap);
         auto tp_3 = std::chrono::high_resolution_clock::now();
         std::string fn = "stage-1-post-make_fragments-all-fragments.pdb";
         m.write_file(fn, 30.0);

         float big_overlap_fraction_limit = 0.89;
         mmdb::Manager *mol = m.pcmmdbmanager();
         auto scores_for_chain_ids = make_scores_for_chain_ids(scored_traces);
         auto tp_4 = std::chrono::high_resolution_clock::now();
         std::cout << "PROGRESS calling find_chains_that_overlap_other_chains()" << std::endl;
         auto deletable_chains_map = find_chains_that_overlap_other_chains(mol, big_overlap_fraction_limit, scores_for_chain_ids);
         auto tp_5 = std::chrono::high_resolution_clock::now();
         auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
         auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
         std::cout << "Timings: from make_fragments() for scored_trees_to_multi_fragment()        " << d32 << " milliseconds" << std::endl;
         std::cout << "Timings: from make_fragments() for find_chains_that_overlap_other_chains() " << d54 << " milliseconds" << std::endl;

         if (false) {
            std::map<std::string, std::set<std::string> >::const_iterator it;
            for (it = deletable_chains_map.begin(); it != deletable_chains_map.end(); ++it) {
               std::cout << "deletable_chains_map: " << it->first << " : ";
               for (const auto &item : it->second)
                  std::cout << " " << item;
               std::cout << std::endl;
            }
         }
         auto tp_6 = std::chrono::high_resolution_clock::now();
         filter_similar_chains(mol, deletable_chains_map); // don't think - just do.
         auto tp_7 = std::chrono::high_resolution_clock::now();
         auto d76 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_7 - tp_6).count();
         std::cout << "Timings: from make_fragments() for scored_trees_to_multi_fragment()        " << d76 << " milliseconds" << std::endl;

         mol->WritePDBASCII("stage-2-post-overlapping-chain-filter.pdb");

         if (false)
            write_pdbs_for_each_of_the_chains(mol);
         return mol;
      }
   }

}


void
globularize(mmdb::Manager *mol, const clipper::Xmap<float> &xmap, const clipper::Coord_orth &hack_centre, bool use_hack_centre) {

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

         if (use_hack_centre) mean = hack_centre;

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

// name this?
bool
this_and_next_are_close(mmdb::Residue *residue_this_p, mmdb::Residue *residue_next_p) {

   bool status = false;

   if (residue_this_p && residue_next_p) {
      int n_atoms_in_this_res = residue_this_p->GetNumberOfAtoms();
      for (int iat=0; iat<n_atoms_in_this_res; iat++) {
         mmdb::Atom *at_this = residue_this_p->GetAtom(iat);
         if (! at_this->isTer()) {
            std::string atom_name(at_this->GetAtomName());
            if (atom_name == " C  ") {
               clipper::Coord_orth c_pos_this = coot::co(at_this);
               int n_atoms_in_next_res = residue_next_p->GetNumberOfAtoms();
               for (int jat=0; jat<n_atoms_in_next_res; jat++) {
                  mmdb::Atom *at_next = residue_next_p->GetAtom(jat);
                  if (! at_next->isTer()) {
                     std::string atom_name_next(at_next->GetAtomName());
                     if (atom_name_next == " N  ") {
                        clipper::Coord_orth n_pos_next = coot::co(at_next);
                        double dd = (n_pos_next-c_pos_this).lengthsq();
                        double d = std::sqrt(dd);
                        std::cout << "debug this_and_next_are_close() "
                                  << coot::atom_spec_t(at_this) << " " << coot::atom_spec_t(at_next) << " " << d << std::endl;
                        if (d < 5.0) { // !! what number here?
                           status = true;
                        }
                        break;
                     }
                  }
               }
               if (status) break;
            }
         }
      }
   }

   std::cout << "debug:: this_and_next_are_close(): "  << coot::residue_spec_t(residue_this_p) << " " << coot::residue_spec_t(residue_next_p)
             << " close-status: " << status << std::endl;

   return status;
}


void
bring_together_consecutive_C_and_N_by_symmetry_transformation(mmdb::Manager *mol, const clipper::Xmap<float> &xmap) {

   class simple_symm_info_t {
   public:
      int isym;
      int i_tr_x;
      int i_tr_y;
      int i_tr_z;
      simple_symm_info_t(int isym, int i_tr_x, int i_tr_y, int i_tr_z) : isym(isym), i_tr_x(i_tr_x), i_tr_y(i_tr_y), i_tr_z(i_tr_z) {}
      simple_symm_info_t() : isym(0), i_tr_x(0), i_tr_y(0), i_tr_z(0) {}
   };

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      const clipper::Spacegroup &spgr = xmap.spacegroup();
      const clipper::Cell       &cell = xmap.cell();

      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<(nres-1); ires++) {
            mmdb::Residue *residue_this_p = chain_p->GetResidue(ires);
            mmdb::Residue *residue_next_p = chain_p->GetResidue(ires+1);
            if (residue_this_p && residue_next_p) {

               int res_no_this = residue_this_p->GetSeqNum();
               int res_no_next = residue_next_p->GetSeqNum();
               if (res_no_next != (res_no_this + 1))
                  continue;

               int n_atoms_in_this_res = residue_this_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms_in_this_res; iat++) {
                  mmdb::Atom *at_this = residue_this_p->GetAtom(iat);
                  if (! at_this->isTer()) {
                     std::string atom_name(at_this->GetAtomName());
                     if (atom_name == " C  ") {
                        clipper::Coord_orth c_pos_this = coot::co(at_this);
                        int n_atoms_in_next_res = residue_next_p->GetNumberOfAtoms();
                        for (int jat=0; jat<n_atoms_in_next_res; jat++) {
                           mmdb::Atom *at_next = residue_next_p->GetAtom(jat);
                           if (! at_next->isTer()) {
                              std::string atom_name_next(at_next->GetAtomName());
                              if (atom_name_next == " N  ") {
                                 clipper::Coord_orth n_pos_next = coot::co(at_next);
                                 clipper::Coord_orth n_pos_frac(n_pos_next.coord_frac(cell));
                                 double dd = (n_pos_next-c_pos_this).lengthsq();
                                 double d = std::sqrt(dd);

                                 if (d > 5.0) { // !! what number here?
                                    std::cout << "bring_together_consecutive_C_and_N_by_symmetry_transformation() found a discrepancy that needs fixing "
                                              << coot::residue_spec_t(residue_this_p) << " " << coot::residue_spec_t(residue_next_p)
                                              << " " << d << std::endl;
                                    clipper::Coord_orth pos_best = n_pos_next;
                                    double d_best_sqrd = dd;
                                    simple_symm_info_t sym_best;
                                    for(int i_tr_x=-3; i_tr_x<=3; i_tr_x++) {
                                       for(int i_tr_y=-3; i_tr_y<=3; i_tr_y++) {
                                          for(int i_tr_z=-3; i_tr_z<=3; i_tr_z++) {
                                             for (int i_symm=0; i_symm<spgr.num_symops(); i_symm++) {
                                                clipper::Coord_frac t(i_tr_x, i_tr_y, i_tr_z);
                                                // clipper::Coord_frac cpt(n_pos_frac + t);
                                                // clipper::Coord_frac cfs(spgr.symop(i_symm) * cpt);
                                                clipper::Coord_frac cfs(spgr.symop(i_symm) * n_pos_frac + t);
                                                clipper::Coord_orth p(cfs.coord_orth(cell));
                                                double d_this_sqrd = (c_pos_this -p).lengthsq();
                                                if (false)
                                                   std::cout << "sym for test: " <<  i_symm << " " << i_tr_x << " " << i_tr_y << " " << i_tr_z
                                                             << " dist " << std::sqrt(d_this_sqrd) << std::endl;
                                                if (d_this_sqrd < d_best_sqrd) {
                                                   d_best_sqrd = d_this_sqrd;
                                                   pos_best = p;
                                                   sym_best = simple_symm_info_t(i_symm, i_tr_x, i_tr_y, i_tr_z);
                                                   if (false)
                                                      std::cout << "updating sym_best " <<  i_symm << " " << i_tr_x << " " << i_tr_y << " " << i_tr_z
                                                                << " " << std::sqrt(d_best_sqrd) << std::endl;
                                                }
                                             }
                                          }
                                       }
                                    }
                                    if (d_best_sqrd < dd) {

                                       if (sqrt(d_best_sqrd) > 10.0) { // don't bother moving stuff if it doesn't anneal

                                          std::cout << "bring_together_consecutive_C_and_N_by_symmetry_transformation() no move stuff "
                                                    << " because best " << std::sqrt(d_best_sqrd) << " orig " << d << " doesn't anneal "
                                                    << std::endl;
                                          
                                       } else {

                                          // happy path

                                          // now lets make a vector (now a set) of residues in the same fragment as this residue_next_p
                                          std::cout << "bring_together_consecutive_C_and_N_by_symmetry_transformation() Let's move some stuff "
                                                    << " because " << std::sqrt(d_best_sqrd) << " beats " << d << std::endl;
                                          std::set<mmdb::Residue *> fragment_residues;

                                          if (true) {
                                             for (int i_res_frag=(ires+1); i_res_frag<(nres-1); i_res_frag++) {
                                                mmdb::Residue *r_this_p = chain_p->GetResidue(i_res_frag);
                                                std::cout << "chain_p residue sanity check " << coot::residue_spec_t(r_this_p) << std::endl;
                                             }
                                          }
                                       
                                          for (int i_res_frag=(ires+1); i_res_frag<(nres-1); i_res_frag++) {
                                             mmdb::Residue *r_this_p = chain_p->GetResidue(i_res_frag);
                                             mmdb::Residue *r_next_p = chain_p->GetResidue(i_res_frag+1);

                                             int rres_no_this = r_this_p->GetSeqNum();
                                             int rres_no_next = r_next_p->GetSeqNum();
                                             if (rres_no_next != (rres_no_this + 1))
                                                break;
                                             
                                             if (this_and_next_are_close(r_this_p, r_next_p)) {
                                                std::cout << "inserting 1 into fragment_residues: " << coot::residue_spec_t(r_this_p) << std::endl;
                                                std::cout << "inserting 2 into fragment_residues: " << coot::residue_spec_t(r_next_p) << std::endl;
                                                fragment_residues.insert(r_this_p);
                                                fragment_residues.insert(r_next_p);
                                             } else {
                                                break;
                                             }
                                          }

                                          std::cout << "bring_together_consecutive_C_and_N_by_symmetry_transformation() fragment_residues size "
                                                    << fragment_residues.size() << std::endl;

                                          if (true) {
                                             std::set<mmdb::Residue *>::iterator it;
                                             for (it=fragment_residues.begin(); it!=fragment_residues.end(); ++it) {
                                                mmdb::Residue *residue_p(*it);
                                                std::cout << "deubg:: fragment residue " << coot::residue_spec_t(residue_p) << std::endl;
                                             }
                                          }

                                          if (! fragment_residues.empty()) {
                                             std::set<mmdb::Residue *>::iterator it;
                                             for (it=fragment_residues.begin(); it!=fragment_residues.end(); ++it) {
                                                mmdb::Residue *residue_p(*it);
                                                int n_atoms = residue_p->GetNumberOfAtoms();
                                                for (int iat_res=0; iat_res<n_atoms; iat_res++) {
                                                   mmdb::Atom *at_res = residue_p->GetAtom(iat_res);
                                                   clipper::Coord_orth at_pos = coot::co(at_res);
                                                   clipper::Coord_orth at_pos_frac(at_pos.coord_frac(cell));
                                                   clipper::Coord_frac t(sym_best.i_tr_x, sym_best.i_tr_y, sym_best.i_tr_z);
                                                   // clipper::Coord_frac cpt(at_pos_frac + t);
                                                   // clipper::Coord_frac cfs(spgr.symop(sym_best.isym) * cpt);
                                                   clipper::Coord_frac cfs(spgr.symop(sym_best.isym) * at_pos_frac + t);
                                                   clipper::Coord_orth p(cfs.coord_orth(cell));
                                                   std::cout << "Moving " << coot::atom_spec_t(at_res) << " from "
                                                             << at_res->x << " " << at_res->y << " " << at_res->z << " to "
                                                             << p.x() << " " << p.y() << " " << p.z() << std::endl;
                                                   at_res->x = p.x(); at_res->y = p.y(); at_res->z = p.z(); 
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
                         const clipper::Xmap<float> &xmap,
                         double variation,
                         unsigned int n_top_spin_pairs,
                         unsigned int n_top_fragments,
                         bool cryo_em_mode,
                         std::pair<bool, clipper::Coord_orth> hack_centre) {

   auto debug_scored_spin_pairs = [] (const std::vector<std::pair<unsigned int, coot::scored_node_t> > &scored_pairs,
                                      mmdb::Atom **atom_selection, int n_selected_atoms) {

                                     double max_score = 26.0; // for 1gwd
                                     // xmax_score = 250.0; // for emd-22898

                                     auto score_to_colour = [max_score] (double score) {
                                                               float f = score/max_score;
                                                               if (f < 0.0) f = 0.0;
                                                               if (f > 1.0) f = 1.0;
                                                               if (f < 0.7) f = 0.0;
                                                               float ff = -1.6 * f;
                                                               coot::colour_holder ch(0.2, 0.7, 0.3);
                                                               ch.rotate_by(ff);
                                                               return ch;
                                                            };

                                     std::ofstream f("debug-scored-peptides.table");
                                     for (unsigned int i=0; i<scored_pairs.size(); i++) {
                                        // if (i> 10) continue;
                                        const auto &scored_pair = scored_pairs[i];
                                        clipper::Coord_orth pt_1 = coot::co(atom_selection[scored_pair.first]);
                                        clipper::Coord_orth pt_2 = coot::co(atom_selection[scored_pair.second.atom_idx]);
                                        double score = scored_pair.second.spin_score;
                                        coot::colour_holder ch = score_to_colour(score);
                                        if (score > 2.0) {
                                           f << "scored-peptide idx_1 " << scored_pair.first << " idx_2 " << scored_pair.second.atom_idx << " "
                                             << std::setw(9) << pt_1.x() << " " << std::setw(9) << pt_1.y() << " " << std::setw(9) << pt_1.z() << " "
                                             << std::setw(9) << pt_2.x() << " " << std::setw(9) << pt_2.y() << " " << std::setw(9) << pt_2.z()
                                             << "  score: " << score <<  " col " << ch.red << " " << ch.green << " " << ch.blue
                                             << "  for score-ratio " << score/max_score << "\n";
                                        }
                                     }
                                  };

   // somewhere here - not sure before or after action_mol is created, I want to globularize the molecule
   mmdb::Manager *action_mol = flood_mol.pcmmdbmanager();

   globularize(action_mol, xmap, hack_centre.second, hack_centre.first); // move around the atoms so they they are arranged in space in a sphere rather
                                                                // than in a strip of the map (the asymmetric unit).

   action_mol->WritePDBASCII("flood-mol-globularized.pdb");

   mmdb::Atom **atom_selection = 0; // member data - cleared on destruction
   int n_selected_atoms = 0;
   int selhnd = action_mol->NewSelection(); // d
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

   std::cout << "PROGRESS:: calling make_spin_scored_pairs() using " << apwd.size() << " atom pairs within distance" << std::endl;
   // the first of the scores is the index of the first atom

   std::vector<std::pair<unsigned int, coot::scored_node_t> > scored_pairs =
      make_spin_scored_pairs(apwd, n_top_spin_pairs, xmap, action_mol, atom_selection, n_selected_atoms); // sorted
   std::cout << "spin_score_pairs done" << std::endl;
   if (true)
      debug_scored_spin_pairs(scored_pairs, atom_selection, n_selected_atoms);

   // scores = make_some_fake_scored_pairs();
   mmdb::Manager *mol = make_fragments(scored_pairs, atom_selection, xmap, n_top_fragments, cryo_em_mode);

   float acell[6];
   acell[0] = xmap.cell().descr().a();
   acell[1] = xmap.cell().descr().b();
   acell[2] = xmap.cell().descr().c();
   acell[3] = clipper::Util::rad2d(xmap.cell().descr().alpha());
   acell[4] = clipper::Util::rad2d(xmap.cell().descr().beta());
   acell[5] = clipper::Util::rad2d(xmap.cell().descr().gamma());
   std::string spacegroup_str_hm = xmap.spacegroup().symbol_hm();
   mol->SetCell(acell[0], acell[1], acell[2], acell[3], acell[4], acell[5]);
   mol->SetSpaceGroup(spacegroup_str_hm.c_str());

   action_mol->DeleteSelection(selhnd); // finished with this now

   return mol;
}

// twisted trans or cis
bool
peptide_is_twisted(mmdb::Residue *residue_with_CO, mmdb::Residue *residue_with_N, double deformation_limit_deg=30.0) {

   bool status = false;
   mmdb::Atom *CA_1 = residue_with_CO->GetAtom(" CA ");
   mmdb::Atom *C_1  = residue_with_CO->GetAtom(" C  ");
   mmdb::Atom *N_2  = residue_with_N->GetAtom(" N  ");
   mmdb::Atom *CA_2 = residue_with_N->GetAtom(" CA ");
   if (CA_1 && C_1 && N_2 && CA_2) {
      clipper::Coord_orth ca_1_pt = coot::co(CA_1);
      clipper::Coord_orth c_1_pt  = coot::co(C_1);
      clipper::Coord_orth n_2_pt  = coot::co(N_2);
      clipper::Coord_orth ca_2_pt = coot::co(CA_2);
      double torsion = clipper::Coord_orth::torsion(ca_1_pt, c_1_pt, n_2_pt, ca_2_pt);
      double torsion_deg = clipper::Util::rad2d(torsion);
      if (torsion_deg > (-180.0 + deformation_limit_deg))
         if (torsion_deg < (180.0 - deformation_limit_deg))
            status = true;
      if (status) {
         std::cout << "Torsion check  "
                   << coot::residue_spec_t(residue_with_CO) << " "
                   << coot::residue_spec_t(residue_with_N) << " torsion "
                   << std::setw(7) << torsion << " in degs: " << std::setw(7) << torsion_deg << " Baddie\n";
      }
   } else {
      std::cout << "ERROR:: peptide_is_twisted(): missing atoms torsion " << std::endl;
   }
   return status;
};


void
delete_chains_that_have_twisted_trans_peptides(mmdb::Manager *mol, unsigned int n_per_chain_max, double deformation_limit_deg=30.0) {

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);

   if (model_p) {
      std::set<mmdb::Chain *> this_chain_is_not_twisted_set; // don't keep testing unmodified chains.
      bool keep_looping = false; // unless we delete something. Only one chain delete per loop is allowed
      do {
         keep_looping = false;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (this_chain_is_not_twisted_set.find(chain_p) != this_chain_is_not_twisted_set.end()) continue;
            bool twisted_chain = false;

            unsigned int n_twisted_peptide = 0;
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<(n_res-1); ires++) {
               mmdb::Residue *residue_this_p = chain_p->GetResidue(ires);
               if (residue_this_p) {
                  mmdb::Residue *residue_next_p = chain_p->GetResidue(ires+1);
                  if (residue_next_p) {
                     bool twisted = peptide_is_twisted(residue_this_p, residue_next_p, deformation_limit_deg);
                     if (twisted) {
                        n_twisted_peptide++;
                        twisted_chain = true;
                     }
                  }
               }
            }
            if (! twisted_chain)
               this_chain_is_not_twisted_set.insert(chain_p);

            if (n_twisted_peptide > n_per_chain_max) {
               std::string chain_id(chain_p->GetChainID());
               std::cout << "INFO:: delete_chains_that_have_twisted_trans_peptides(): deleting chain " << chain_id << std::endl;
               model_p->DeleteChain(ichain);
               mol->FinishStructEdit();
               keep_looping = true;
               break;
            }
         }
      } while (keep_looping);
   }
}

void
delete_chains_that_are_too_short(mmdb::Manager *mol, int n_res_min) {

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);

   if (model_p) {

      // what is the maxium chain length?
      int max_chain_length = 0;
      int n_chains_start = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains_start; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         if (n_res > max_chain_length)
            max_chain_length = n_res;
      }

      // don't kill off everything - even if told to.
      if (n_res_min > max_chain_length)
         n_res_min = max_chain_length;

      bool keep_looping = false; // unless we delete something. Only one chain delete per loop is allowed
      do {
         keep_looping = false;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (chain_p) {  // probably not needed now that I've fixed the bug
               int n_res = chain_p->GetNumberOfResidues();
               if (n_res < n_res_min) {
                  std::string chain_id(chain_p->GetChainID());
                  std::cout << "INFO:: delete_chains_that_are_too_short(): deleting chain " << chain_id << std::endl;
                  model_p->DeleteChain(ichain);
                  mol->FinishStructEdit();
                  keep_looping = true;
                  break;
               }
            }
         }
      } while (keep_looping);
   }
}


// mutate mol, returning the map of the chain id to the index in the fam of the sequnece that
// it was assigned
//
std::map<std::string, unsigned int>
apply_sequence_to_fragments(mmdb::Manager *mol_in, const clipper::Xmap<float> &xmap, const coot::fasta_multi &fam,
                            const coot::protein_geometry &pg) {

   // This needs some clear thinking about what to do when there are 2 sequnces in the fam.
   // which sequnce should be applied?
   // Perhaps the sequence loop should be inside the chain loop, and then we can gather the scores
   // from each of the chains?  Also, try to consider if the chain is sufficiently ambiguous that
   // it is not assignable - get_result() should provide that information.

   auto print_atoms_in_chain = [] (const std::string &tag, mmdb::Chain *chain_p) {
                                  int n_res = chain_p->GetNumberOfResidues();
                                  for (int ires=0; ires<n_res; ires++) {
                                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                     std::string res_name(residue_p->GetResName());
                                     if (residue_p) {
                                        int n_atoms = residue_p->GetNumberOfAtoms();
                                        for (int iat=0; iat<n_atoms; iat++) {
                                           mmdb::Atom *at = residue_p->GetAtom(iat);
                                           if (! at->isTer()) {
                                              std::cout << tag << " " << "atom: " << coot::atom_spec_t(at) << " " << res_name << std::endl;
                                           }
                                        }
                                     }
                                  }
                               };

   std::map<std::string, unsigned int> chain_id_to_fam_index;

   unsigned int n_sequences = fam.size();
   std::cout << "debug:: apply_sequence_to_fragments(): with n_sequences " << n_sequences << std::endl;
   if (n_sequences > 0) {
      for (unsigned int idx=0; idx<n_sequences; idx++) {
         std::string sequence = fam[idx].sequence;
         std::cout << "debug sequence: " << sequence << std::endl;
         const std::string &name = fam[idx].name;

         int imod = 1;
         mmdb::Model *model_p = mol_in->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               std::string chain_id(chain_p->GetChainID());
               int n_res = chain_p->GetNumberOfResidues();
               if (n_res > 5) {
                  int resno_start = chain_p->GetResidue(0)->GetSeqNum();
                  int resno_end   = chain_p->GetResidue(n_res-1)->GetSeqNum();

                  coot::side_chain_densities scd;

                  std::pair<mmdb::Chain *, mmdb::Manager *> chain_mol_pair = coot::util::copy_chain(chain_p); // add it into a molecule hierarchy
                  
                  std::pair<std::string, std::vector<mmdb::Residue *> > a_run_of_residues =
                     scd.setup_test_sequence(chain_mol_pair.second, chain_id, resno_start, resno_end, xmap);

                  scd.fill_residue_blocks(a_run_of_residues.second, xmap); // return fast if already filled, uses atomic locking.

                  std::cout << "debug:: apply_sequence_to_fragments(): a run of residues: " << a_run_of_residues.first << " with "
                            << a_run_of_residues.second.size() << " residues" << std::endl;
                  if (a_run_of_residues.first.empty()) {
                     scd.test_sequence(a_run_of_residues.second, xmap, name, sequence);
                  } else {
                     std::cout << "ERROR:: apply_sequence_to_fragments(): when generating a run of residues:" << std::endl;
                     std::cout << "        " << a_run_of_residues.first << std::endl;
                  }

                  // I want get_result() to give me the the rotamer name of each of the sidechains too!
                  bool probable_only = true;
                  bool print_sequencing_solutions_flag = true;
                  coot::side_chain_densities::results_t new_sequence_result = scd.get_result(probable_only, print_sequencing_solutions_flag);
                  std::string new_sequence = new_sequence_result.sequence;
                  std::cout << "debug:: apply_sequence_to_fragments(): new_sequence " << new_sequence << std::endl;
                  if (! new_sequence.empty()) {
                     int sl = new_sequence.length();
                     if (sl == n_res) {

                        // OK, so attach this chain-ID to ths fam index
                        chain_id_to_fam_index[chain_id] = idx;

                        // Find the offset of the sequence - so that I can get the residue numbers correct (even if the chains
                        // are fragmented and overlaying here)
                        std::string::size_type seq_pos = sequence.find(new_sequence);
                        
                        for (int ires=0; ires<n_res; ires++) {
                           mmdb::Residue *residue_p = chain_mol_pair.first->GetResidue(ires);
                           if (residue_p) {
                              char letter = new_sequence[ires];
                              std::string new_residue_type = coot::util::single_letter_to_3_letter_code(letter);
                              coot::util::mutate(residue_p, new_residue_type);
                           } else {
                              std::cout << "WARNING:: null residue in apply_sequence_to_fragments() " << ires << std::endl;
                           }
                        }
                        // now copy the atoms of chain_copy_p into chain_p
                        bool do_finishstructedit = true;
                        coot::util::replace_chain_contents_with_atoms_from_chain(chain_p, mol_in, chain_mol_pair.first, do_finishstructedit); // (to_chain, from_chain)
                        std::cout << "debug:: apply_sequence_to_fragments(): seq_pos " << seq_pos << std::endl;
                        if (seq_pos != std::string::npos) {
                           n_res = chain_p->GetNumberOfResidues();
                           for (int ires=0; ires<n_res; ires++) {
                              mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                              residue_p->seqNum = seq_pos + ires + 1;
                              // std::cout << "new seqnum for residue with index " << ires << " " << residue_p->seqNum << std::endl;
                           }
                           mol_in->FinishStructEdit();
                        }

                        // std::string fn = "post-mutate-and-renumber-" + std::to_string(ichain) + ".pdb";
                        // mol_in->WritePDBASCII(fn.c_str());
                        
                     } else {
                        std::cout << "WARNING:: Ooops sl != n_res " << sl << " " << n_res << std::endl;
                     }
                  } else {
                     std::cout << "INFO:: no clear solution found for chain " << chain_id << std::endl;
                  }
               } else {
                  std::cout << "WARNING:: Skip trace n_res > 5 test failed " << n_res << " for chain " << chain_id << std::endl;
               }
            }

            // we can't do backrub here - we have not filtered out the overlapping chains.
            // Hmm... maybe there should be "ignore clashes" mode?
            // coot::backrub_molecule(mol_in, &xmap, pg); 
         }
      }
   }
   return chain_id_to_fam_index;
}

#include "residue_by_phi_psi.hh"

// working_mol is so that the graphics can show were we are in the building. Change count_p when working_mol gets updated.
// (but don't do that if count_p is null).
//
void res_tracer_proc(const clipper::Xmap<float> &xmap, const coot::fasta_multi &fam, double variation,
                     unsigned int n_top_spin_pairs, unsigned int n_top_fragments,
                     float rmsd_cut_off_for_flood, float flood_atom_mask_radius, float weight, unsigned int n_phi_psi_trials,
                     bool with_ncs, watch_res_tracer_data_t *watch_res_tracer_data_p) {

   auto get_mol_edit_lock = [watch_res_tracer_data_p] (const std::string &locker) mutable {
                               std::cout << "debug:: locker: " << locker
                                         << " trying to get the lock with mol_edit_lock: " << watch_res_tracer_data_p->mol_edit_lock << std::endl;
                               bool unlocked = false;
                               while (! watch_res_tracer_data_p->mol_edit_lock.compare_exchange_weak(unlocked, true)) {
                                  std::this_thread::sleep_for(std::chrono::microseconds(100));
                                  unlocked = false;
                               }
                               std::cout << "debug:: locker: " << locker << " got the lock" << std::endl;
                            };

   auto release_mol_edit_lock = [watch_res_tracer_data_p] (const std::string &locker) mutable {
                                   watch_res_tracer_data_p->mol_edit_lock = false;
                                   std::cout << "debug:: locker: " << locker << " released the lock" << std::endl;
                                };

   auto density_and_omega_based_trim_chain_terminii = [] (mmdb::Manager *mol, const clipper::Xmap<float> &xmap) {

                                               int imod = 1;
                                               mmdb::Model *model_p = mol->GetModel(imod);

                                               if (model_p) {

                                                  // find the average density for each residue and remove the terminal residues of each chain
                                                  // if they are low outliers
                                                  std::map<coot::residue_spec_t, coot::stats::single> residue_density_stats;
                                                  int n_chains = model_p->GetNumberOfChains();
                                                  for (int ichain=0; ichain<n_chains; ichain++) {
                                                     mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                                     int n_res = chain_p->GetNumberOfResidues();
                                                     for (int ires=0; ires<n_res; ires++) {
                                                        mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                                        if (residue_p) {
                                                           coot::residue_spec_t res_spec(residue_p);
                                                           coot::stats::single stats;
                                                           int n_atoms = residue_p->GetNumberOfAtoms();
                                                           for (int iat=0; iat<n_atoms; iat++) {
                                                              mmdb::Atom *at = residue_p->GetAtom(iat);
                                                              if (! at->isTer()) {
                                                                 clipper::Coord_orth pt = coot::co(at);
                                                                 double dv = coot::util::density_at_point(xmap, pt);
                                                                 stats.add(dv);
                                                              }
                                                           }
                                                           residue_density_stats[res_spec] = stats;
                                                           // std::cout << "storing density stats " << res_spec << " " << stats.mean() << std::endl;
                                                        }
                                                     }
                                                  }

                                                  double sum = 0.0;
                                                  unsigned int n = 0;
                                                  std::map<coot::residue_spec_t, coot::stats::single>::const_iterator it;
                                                  for (it=residue_density_stats.begin(); it!=residue_density_stats.end(); ++it) {
                                                     n += it->second.size();
                                                     sum += it->second.mean() * static_cast<float>(it->second.size());
                                                  }
                                                  double mean_density = sum/static_cast<float>(n);
                                                  double density_level_crit = 0.9 * mean_density;

                                                  auto test_density_and_delete = [density_level_crit, residue_density_stats] (mmdb::Residue *residue_p,
                                                                                                                              mmdb::Chain *chain_p,
                                                                                                                              int index) {
                                                                            bool deleted = false;
                                                                            if (residue_p) {
                                                                               std::map<coot::residue_spec_t, coot::stats::single>::const_iterator it;
                                                                               coot::residue_spec_t res_spec(residue_p);
                                                                               it = residue_density_stats.find(res_spec);
                                                                               if (it != residue_density_stats.end()) {
                                                                                  const auto &this_residue_stats = it->second;
                                                                                  if (this_residue_stats.mean() < density_level_crit) {
                                                                                     chain_p->DeleteResidue(index);
                                                                                     std::cout << "Density-based Residue Trim: " << res_spec << std::endl;
                                                                                     deleted = true;
                                                                                  }
                                                                               } else {
                                                                                  std::cout << "failed to find " << res_spec << " in residue density stats" << std::endl;
                                                                               }
                                                                            }
                                                                            return deleted;
                                                                         };

                                                  // residue_with_CO is residue 1 and residue_with_N is residue 2 (say)
                                                  //
                                                  auto test_omega_and_delete = [] (mmdb::Residue *residue_with_CO, mmdb::Residue *residue_with_N,
                                                                                   int index_CO_residue,
                                                                                   int index_N_residue,
                                                                                   mmdb::Chain *chain_p) {
                                                                                  bool status = false; // no deletion
                                                                                  mmdb::Atom *CA_1 = residue_with_CO->GetAtom(" CA ");
                                                                                  mmdb::Atom *C_1  = residue_with_CO->GetAtom(" C  ");
                                                                                  mmdb::Atom *N_2  = residue_with_N->GetAtom(" N  ");
                                                                                  mmdb::Atom *CA_2 = residue_with_N->GetAtom(" CA ");
                                                                                  if (CA_1 && C_1 && N_2 && CA_2) {
                                                                                     clipper::Coord_orth ca_1_pt = coot::co(CA_1);
                                                                                     clipper::Coord_orth c_1_pt  = coot::co(C_1);
                                                                                     clipper::Coord_orth n_2_pt  = coot::co(N_2);
                                                                                     clipper::Coord_orth ca_2_pt = coot::co(CA_2);
                                                                                     double torsion = clipper::Coord_orth::torsion(ca_1_pt, c_1_pt, n_2_pt, ca_2_pt);
                                                                                     double torsion_deg = clipper::Util::rad2d(torsion);
                                                                                     if (torsion_deg > -160.0)
                                                                                        if (torsion_deg < 160.0)
                                                                                           status = true;
                                                                                     if (true) { // debug
                                                                                        std::cout << "Torsion check  "
                                                                                                  << coot::residue_spec_t(residue_with_CO) << " "
                                                                                                  << coot::residue_spec_t(residue_with_N) << " torsion "
                                                                                                  << torsion << " in degs: " << torsion_deg;
                                                                                        if (status)
                                                                                           std::cout << " Baddie\n";
                                                                                        else
                                                                                           std::cout << "\n";
                                                                                     }
                                                                                     if (status) {
                                                                                        chain_p->DeleteResidue(index_N_residue);
                                                                                        chain_p->DeleteResidue(index_CO_residue);
                                                                                     }
                                                                                  } else {
                                                                                     std::cout << "ERROR:: test_omega_and_delete() missing atoms torsion " << std::endl;
                                                                                  }
                                                                                  return status;
                                                                               };


                                                  for (int ichain=0; ichain<n_chains; ichain++) {
                                                     mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                                     if (chain_p) {
                                                        bool continue_looping = true;
                                                        while (continue_looping) {
                                                           continue_looping = false;
                                                           int n_res = chain_p->GetNumberOfResidues();

                                                           if (n_res > 10000) {
                                                              std::cout << "ERROR:: Trapped a bug! n_res " << n_res << " ichain " << ichain << " chain_p" << chain_p << std::endl;
                                                              continue_looping = false;
                                                              break;
                                                           }

                                                           if (n_res > 0) {
                                                              mmdb::Residue *residue_p = chain_p->GetResidue(n_res-1);
                                                              // std::cout << "debug:: residue_p (end of chain) " << coot::residue_spec_t(residue_p) << std::endl;
                                                              bool deleted = test_density_and_delete(residue_p, chain_p, n_res-1);
                                                              if (deleted) {
                                                                 residue_p = 0;
                                                                 continue_looping = true;
                                                              }

                                                              if (! deleted) {
                                                                 // try the front end
                                                                 mmdb::Residue *residue_front_p = chain_p->GetResidue(0);
                                                                 if (residue_p) {
                                                                    bool start_was_deleted = test_density_and_delete(residue_front_p, chain_p, 0);
                                                                    if (start_was_deleted) {
                                                                       deleted = true;
                                                                       continue_looping = true;
                                                                       residue_front_p = 0;
                                                                    }
                                                                 }
                                                              }

                                                              if (! deleted) {
                                                                 mmdb::Residue *prev_residue_with_CO = chain_p->GetResidue(n_res-2); // "(which has CO)" I mean
                                                                 if (residue_p && prev_residue_with_CO) {
                                                                    if (false)
                                                                       std::cout << "calling A torsion check for "
                                                                                 << coot::residue_spec_t(prev_residue_with_CO) << " " << coot::residue_spec_t(residue_p)
                                                                                 << std::endl;
                                                                    deleted = test_omega_and_delete(prev_residue_with_CO, residue_p, n_res-2, n_res-1, chain_p);
                                                                    if (deleted) {
                                                                       continue_looping = true;
                                                                       residue_p = 0;
                                                                    }
                                                                 }
                                                              }

                                                              if (! deleted) {
                                                                 if (n_res > 1) {
                                                                    mmdb::Residue *residue_front_p = chain_p->GetResidue(0);
                                                                    mmdb::Residue *next_residue_with_N = chain_p->GetResidue(1);
                                                                    if (residue_front_p && next_residue_with_N) {
                                                                       if (false)
                                                                          std::cout << "calling B torsion check for "
                                                                                    << coot::residue_spec_t(residue_front_p) << " "
                                                                                    << coot::residue_spec_t(next_residue_with_N) << std::endl;
                                                                       deleted = test_omega_and_delete(residue_front_p, next_residue_with_N, 0, 1, chain_p);
                                                                       if (deleted) {
                                                                          continue_looping = true;
                                                                       }
                                                                    }
                                                                 }
                                                              }

                                                              if (deleted)
                                                                 mol->FinishStructEdit();
                                                           }
                                                        }
                                                     }
                                                  }
                                                  mol->FinishStructEdit();
                                                  coot::util::pdbcleanup_serial_residue_numbers(mol);
                                               }
                                            };

   auto refine_isolated_chain = [&xmap] (mmdb::Chain *chain_p, mmdb::Manager *mol_for_this_chain, const coot::protein_geometry &geom,
                                         ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight) {
                                   std::vector<std::pair<bool, mmdb::Residue *> > residues;
                                   int n_res = chain_p->GetNumberOfResidues();
                                   unsigned int n_atoms_for_refinement = 0; // debug crash
                                   for (int ires=0; ires<n_res; ires++) {
                                      mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                      if (residue_p) {
                                         residues.push_back(std::make_pair(false, residue_p));

                                         if (true) {
                                            int n_atoms = residue_p->GetNumberOfAtoms();
                                            for (int iat=0; iat<n_atoms; iat++) {
                                               mmdb::Atom *at = residue_p->GetAtom(iat);
                                               if (! at->isTer()) {
                                                  n_atoms_for_refinement++;
                                               }
                                            }
                                         }
                                      }
                                   }

                                   if (false) {
                                      std::string chain_id(chain_p->GetChainID());
                                      std::cout << "in refine_isolated_chain(): chain " << chain_id
                                                << " n_atoms_for_refinement: " << n_atoms_for_refinement << std::endl;
                                   }

                                   std::vector<mmdb::Link> links;
                                   std::vector<coot::atom_spec_t> fixed_atom_specs;
                                   coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;
                                   coot::restraints_container_t restraints(residues, links, geom, mol_for_this_chain, fixed_atom_specs, &xmap);
                                   restraints.thread_pool(thread_pool_p, n_threads);
                                   restraints.set_quiet_reporting();
                                   coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
                                   bool do_internal_torsions = false;
                                   restraints.add_map(weight);
                                   int imol = 0;
                                   restraints.make_restraints(imol, geom, flags, do_internal_torsions, false, 0, 0, true, true, false, pseudos);
                                   restraints.minimize(flags);
                                };


   // these scores are density at the residue atom position scores, not peptide scores
   //
   auto make_density_fit_scores_for_chains =  [] (mmdb::Manager *mol, const clipper::Xmap<float> &xmap) {
                                                 std::map<std::string, double> scores_for_chain_ids;
                                                 int imod = 1;
                                                 mmdb::Model *model_p = mol->GetModel(imod);
                                                 if (model_p) {
                                                    int n_chains = model_p->GetNumberOfChains();
                                                    for (int ichain=0; ichain<n_chains; ichain++) {
                                                       mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                                       std::string chain_id(chain_p->GetChainID());
                                                       int n_res = chain_p->GetNumberOfResidues();
                                                       for (int ires=0; ires<n_res; ires++) {
                                                          mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                                          if (residue_p) {
                                                             int n_atoms = residue_p->GetNumberOfAtoms();
                                                             for (int iat=0; iat<n_atoms; iat++) {
                                                                mmdb::Atom *at = residue_p->GetAtom(iat);
                                                                if (! at->isTer()) {
                                                                   clipper::Coord_orth pt = coot::co(at);
                                                                   float d = coot::util::density_at_point(xmap, pt);
                                                                   scores_for_chain_ids[chain_id] += d;
                                                                }
                                                             }
                                                          }
                                                       }
                                                    }
                                                 }
                                                 return scores_for_chain_ids;
                                              };

   auto write_chain = [](mmdb::Chain *chain_p, mmdb::Manager *mol, const std::string &file_name, bool as_cif=false) {
                         std::pair<mmdb::Chain *, mmdb::Manager *> chain_mol_pair = coot::util::copy_chain(chain_p); // add it into a molecule hierarchy
                          mmdb::Manager *mol_for_chain_copy = chain_mol_pair.second;
                          if (as_cif)
                             mol_for_chain_copy->WriteCIFASCII(file_name.c_str());
                          else
                             mol_for_chain_copy->WritePDBASCII(file_name.c_str());
                          delete mol_for_chain_copy;
                      };

   // refinement now uses mol to determine the residues to be used for the environment.
   // So refining a chain will cause explosion when there are overlapping chains.
   //
   // So each chain should be refined in isolation and the coordinates copied over - meh.
   //
   auto refine_chain = [refine_isolated_chain]
      (mmdb::Chain *chain_p, mmdb::Manager *mol, const coot::protein_geometry &geom,
       ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight) {

                          std::pair<mmdb::Chain *, mmdb::Manager *> chain_mol_pair = coot::util::copy_chain(chain_p); // add it into a molecule hierarchy
                          mmdb::Chain *chain_copy_p = chain_mol_pair.first;
                          mmdb::Manager *mol_for_chain_copy = chain_mol_pair.second;

                          refine_isolated_chain(chain_copy_p, mol_for_chain_copy, geom, thread_pool_p, n_threads, weight);

                          // now copy the atoms of chain_copy_p into mol_for_chain_copy
                          coot::util::copy_atoms_from_chain_to_chain(chain_copy_p, chain_p);

                          delete mol_for_chain_copy;
                       };

   auto refine_chains = [refine_chain] (const std::vector<unsigned int> &chain_indices, mmdb::Manager *mol,
                            const coot::protein_geometry &geom, ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight) {

                           mmdb::Model *model_p = mol->GetModel(1);
                           if (model_p) {
                              for (auto chain_index : chain_indices) {
                                 mmdb::Chain *chain_p = model_p->GetChain(chain_index);
                                 refine_chain(chain_p, mol, geom, thread_pool_p, n_threads, weight);
                              }
                           }
                        };

   auto multi_refine_individual_chains = [refine_chains] (mmdb::Manager *mol, const coot::protein_geometry &geom,
                                                          ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight) {

                                            mmdb::Model *model_p = mol->GetModel(1);
                                            if (model_p) {
                                               unsigned int n_chains = model_p->GetNumberOfChains();
                                               std::vector<std::vector<unsigned int> > chain_indices;
                                               coot::split_indices(&chain_indices, n_chains, 40);

                                               std::vector<std::thread> threads;
                                               for (unsigned int i=0; i<chain_indices.size(); i++) {
                                                  threads.push_back(std::thread(refine_chains, std::cref(chain_indices[i]), mol, std::cref(geom), thread_pool_p, n_threads, weight));
                                               }
                                               for (unsigned int i=0; i<chain_indices.size(); i++)
                                                  threads[i].join();
                                            }
                                         };

   auto is_simply_insertable_fragment = [] (mmdb::Chain *master_chain, mmdb::Chain *insertable_chain) {

                                           // insertable_chain is a contiguous chain, but master_chain may not be.
                                           // (hmmm... how am I going to make sure that that is true?)

                                           bool status = true;
                                           int res_no_min_insert =  9999;
                                           int res_no_max_insert = -9999;
                                           int n_res_insert = insertable_chain->GetNumberOfResidues();
                                           for (int ires=0; ires<n_res_insert; ires++) {
                                              mmdb::Residue *residue_p = insertable_chain->GetResidue(ires);
                                              if (residue_p) {
                                                 int res_no = residue_p->GetSeqNum();
                                                 if (res_no < res_no_min_insert) res_no_min_insert = res_no;
                                                 if (res_no > res_no_max_insert) res_no_max_insert = res_no;
                                              }
                                           }

                                           int n_res_master = master_chain->GetNumberOfResidues();
                                           for (int ires=0; ires<n_res_master; ires++) {
                                              mmdb::Residue *residue_p = master_chain->GetResidue(ires);
                                              if (residue_p) {
                                                 int res_no = residue_p->GetSeqNum();
                                                 if (res_no >= res_no_min_insert) {
                                                    if (res_no <= res_no_max_insert) {
                                                       status = false;
                                                    }
                                                 }
                                              }
                                           }
                                           return status;
                                        };

   auto get_chains = [] (const std::string &chain_id_1, const std::string &chain_id_2, mmdb::Manager *mol) {
                           mmdb::Chain *chain_1 = 0;
                           mmdb::Chain *chain_2 = 0;
                           mmdb::Model *model_p = mol->GetModel(1);
                           if (model_p) {
                              int n_chains = model_p->GetNumberOfChains();
                              for (int ichain=0; ichain<n_chains; ichain++) {
                                 mmdb::Chain *this_chain_p = model_p->GetChain(ichain);
                                 if (this_chain_p) {
                                    std::string this_chain_id(this_chain_p->GetChainID());
                                    if (this_chain_id == chain_id_1) chain_1 = this_chain_p;
                                    if (this_chain_id == chain_id_2) chain_2 = this_chain_p;
                                 }
                              }
                           }
                           return std::make_pair(chain_1, chain_2);
                   };

   auto add_or_insert_residue = [] (mmdb::Residue *residue_p, mmdb::Chain *chain_p) {
                                   int res_no_insert = residue_p->GetSeqNum();
                                   int n_res_master = chain_p->GetNumberOfResidues();
                                   mmdb::Residue *best_contender = 0;
                                   int idx_best_contender = -1;
                                   for (int ires=0; ires<n_res_master; ires++) {
                                      mmdb::Residue *contender_residue_p = chain_p->GetResidue(ires);
                                      if (contender_residue_p) {
                                         int res_no_contender = contender_residue_p->GetSeqNum();
                                         if (res_no_contender > res_no_insert) {
                                            best_contender = residue_p;
                                            idx_best_contender = ires;
                                            break;
                                         }
                                      }
                                   }
                                   if (best_contender) {
                                      chain_p->InsResidue(residue_p, idx_best_contender);
                                   } else {
                                      // OK, we were not able to insert - let's add
                                      chain_p->AddResidue(residue_p);
                                   }
                                };

   auto simply_insert_chain = [add_or_insert_residue] (mmdb::Chain *chain_master, mmdb::Chain *insertable_chain, mmdb::Manager *mol) {

                                 int n_res_insert = insertable_chain->GetNumberOfResidues();
                                 for (int ires=0; ires<n_res_insert; ires++) {
                                    mmdb::Residue *residue_p = insertable_chain->GetResidue(ires);
                                    if (residue_p) {
                                       mmdb::Residue *rc = coot::util::deep_copy_this_residue(residue_p);
                                       if (rc) {
                                          add_or_insert_residue(rc, chain_master);
                                       }
                                    }
                                 }
                                 mol->FinishStructEdit();
                              };
   
   // I need something more complex because I need to use the model_p to get the chains from the chain-ids.
   // auto chain_length_sorter = [] (const std::string &chain_id_1, const std::string &chain_id_2) {
   // bool status = true;
   // return status;
   // };

   class chain_length_sorter {
   public:
      mmdb::Model *model_p;
      explicit chain_length_sorter(mmdb::Model *model_p_in) : model_p(model_p_in) {}
      bool operator()(const std::string &chain_id_1, const std::string &chain_id_2) const {
         bool status = false;
         mmdb::Chain *chain_1 = model_p->GetChain(chain_id_1.c_str());
         mmdb::Chain *chain_2 = model_p->GetChain(chain_id_2.c_str());
         if (chain_1 && chain_2) {
            int n_res_1 = chain_1->GetNumberOfResidues();
            int n_res_2 = chain_2->GetNumberOfResidues();
            if (n_res_1 > n_res_2)
               return true;
         }
         return status;
      }
   };

   // insert residues from other chain if they are not in master chain
   //
   auto insert_residues_from_other_chains_using_chains = [add_or_insert_residue] (mmdb::Manager *mol,
                                                                                  mmdb::Chain *master_chain, mmdb::Chain *other_chain) {
                                                            if (other_chain) {
                                                               int n_res = other_chain->GetNumberOfResidues();
                                                               if (n_res > 0) {
                                                                  int min_res_no =  9999;
                                                                  int max_res_no = -9999;
                                                                  for (int ires=0; ires<n_res; ires++) {
                                                                     mmdb::Residue *residue_p = other_chain->GetResidue(ires);
                                                                     if (residue_p) {
                                                                        int res_no = residue_p->GetSeqNum();
                                                                        if (res_no > max_res_no) max_res_no = res_no;
                                                                        if (res_no < min_res_no) min_res_no = res_no;
                                                                     }
                                                                  }

                                                                  for (int res_no = min_res_no; res_no <= max_res_no; res_no++) {
                                                                     mmdb::Residue *residue_p = master_chain->GetResidue(res_no, "");
                                                                     if (! residue_p) {
                                                                        // OK can we find this residue other_chain?
                                                                        mmdb::Residue *rem_residue_p = other_chain->GetResidue(res_no, "");
                                                                        if (rem_residue_p) {
                                                                           mmdb::Residue *residue_copy_p = coot::util::deep_copy_this_residue(rem_residue_p);
                                                                           add_or_insert_residue(residue_copy_p, master_chain);
                                                                           mol->FinishStructEdit();
                                                                           break;
                                                                        }
                                                                     }
                                                                  }
                                                               }
                                                            }
                                                         };

   auto insert_residues_from_other_chains = [add_or_insert_residue] (mmdb::Manager *mol, const std::vector<std::string> &chain_ids) {

                                               int imol = 1;
                                               mmdb::Model *model_p = mol->GetModel(imol);
                                               if (!model_p) return;

                                               // Which chains are left?
                                               int n_chains = model_p->GetNumberOfChains();
                                               int min_res_no =  99999;
                                               int max_res_no = -99999;
                                               std::vector<mmdb::Chain *> remaining_chains;
                                               for (int ichain=0; ichain<n_chains; ichain++) {
                                                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                                  std::string chain_id(chain_p->GetChainID()); // do I need to check that this is in chain_ids?
                                                  if (chain_p) {
                                                     if (std::find(chain_ids.begin(), chain_ids.end(), chain_id) != chain_ids.end()) {
                                                        remaining_chains.push_back(chain_p);
                                                        std::cout << "Remaining chain (for this sequence) " << chain_id << std::endl;
                                                        int n_res = chain_p->GetNumberOfResidues();
                                                        for (int ires=0; ires<n_res; ires++) {
                                                           mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                                           if (residue_p) {
                                                              int res_no = residue_p->GetSeqNum();
                                                              if (res_no > max_res_no) max_res_no = res_no;
                                                              if (res_no < min_res_no) min_res_no = res_no;
                                                           }
                                                        }
                                                     }
                                                  }
                                               }
                                               if (remaining_chains.size() > 1) {
                                                  int n_remaining_chains = remaining_chains.size();
                                                  mmdb::Chain *master_chain = remaining_chains[0];
                                                  for (int res_no = min_res_no; res_no <= max_res_no; res_no++) {
                                                     mmdb::Residue *residue_p = master_chain->GetResidue(res_no, "");
                                                     if (! residue_p) {
                                                        // OK can we find this residue in another chain?
                                                        for (int ich=1; ich<n_remaining_chains; ich++) {
                                                           mmdb::Chain *rem_chain_p = remaining_chains[ich];
                                                           if (rem_chain_p) {
                                                              mmdb::Residue *rem_residue_p = rem_chain_p->GetResidue(res_no, "");
                                                              if (rem_residue_p) {
                                                                 mmdb::Residue *residue_copy_p = coot::util::deep_copy_this_residue(rem_residue_p);
                                                                 add_or_insert_residue(residue_copy_p, master_chain);
                                                                 mol->FinishStructEdit();
                                                                 break;
                                                              }
                                                           }
                                                        }
                                                     }
                                                  }
                                                  // now delete all other chains other than the master chain
                                                  for (int ich=1; ich<n_remaining_chains; ich++) {
                                                     mmdb::Chain *rem_chain_p = remaining_chains[ich];
                                                     std::string chain_id(rem_chain_p->GetChainID());
                                                     model_p->DeleteChain(chain_id.c_str());
                                                  }
                                                  mol->FinishStructEdit();
                                               }
                                            };

   auto merge_and_delete_sequenced_fragments_into_single_chain = [get_chains, is_simply_insertable_fragment, simply_insert_chain, insert_residues_from_other_chains,
                                               &xmap, &fam] (mmdb::Manager *mol,
                                                             const std::map<std::string, unsigned int> &chain_id_to_fam_index) {
                                                 // non-trivial function

                                                 mmdb::Model *model_p = mol->GetModel(1);
                                                 if (! model_p) return;

                                                 // first convert the fam indexing
                                                 std::map<std::string, unsigned int>::const_iterator it;
                                                 std::map<unsigned int, std::vector<std::string> > fam_index_to_chain_ids;
                                                 for (it=chain_id_to_fam_index.begin(); it!=chain_id_to_fam_index.end(); ++it) {
                                                    const std::string &chain_id = it->first;
                                                    unsigned int fam_index = it->second;
                                                    fam_index_to_chain_ids[fam_index].push_back(chain_id);
                                                 }

                                                 // now sort the chain ids in those vectors by length (number of residues)
                                                 std::map<unsigned int, std::vector<std::string> >::iterator it_seq;
                                                 for (it_seq=fam_index_to_chain_ids.begin(); it_seq!=fam_index_to_chain_ids.end(); ++it_seq) {
                                                    std::vector<std::string> &chain_ids = it_seq->second;
                                                    std::sort(chain_ids.begin(), chain_ids.end(), chain_length_sorter(model_p));
                                                 }

                                                 for (it_seq=fam_index_to_chain_ids.begin(); it_seq!=fam_index_to_chain_ids.end(); ++it_seq) {
                                                    const std::vector<std::string> &chain_ids(it_seq->second);
                                                    std::set<std::string> done_chains;
                                                    for (const auto &chain_id_1 : chain_ids) {
                                                       if (done_chains.find(chain_id_1) == done_chains.end()) {
                                                          for (const auto &chain_id_2 : chain_ids) {
                                                             if (done_chains.find(chain_id_2) == done_chains.end()) {
                                                                if (chain_id_1 != chain_id_2) {
                                                                   std::pair<mmdb::Chain *, mmdb::Chain *> chain_pair = get_chains(chain_id_1, chain_id_2, mol);
                                                                   if (chain_pair.first && chain_pair.second) {
                                                                      std::cout << "    testing chains " << chain_id_1 << " " << chain_id_2 << std::endl;
                                                                      if (is_simply_insertable_fragment(chain_pair.first, chain_pair.second)) {
                                                                         std::cout << "     insert chain " << chain_id_2 << " into " << chain_id_1 << std::endl;
                                                                         simply_insert_chain(chain_pair.first, chain_pair.second, mol);
                                                                         model_p->DeleteChain(chain_id_2.c_str());
                                                                         mol->FinishStructEdit();
                                                                         done_chains.insert(chain_id_2);
                                                                      }
                                                                   }
                                                                }
                                                             }
                                                          }
                                                       }
                                                    }

                                                    insert_residues_from_other_chains(mol, chain_ids);
                                                 }
                                              };


   class chain_overap_info_t {
   public:
      bool status;
      int n_residues_first_chain;
      int n_residues_second_chain;
      int n_match_res_no_and_res_type;
      chain_overap_info_t() { status = false; n_residues_first_chain = 0; n_residues_second_chain = 0; n_match_res_no_and_res_type = 0; }
      explicit chain_overap_info_t(bool status_in) : status(status_in) { n_residues_first_chain = 0; n_residues_second_chain = 0; n_match_res_no_and_res_type = 0; }
      void set_match_info(int n_residues_1st_chain, int n_residues_2nd_chain, int n_matches) {
         n_residues_first_chain  = n_residues_1st_chain;
         n_residues_second_chain = n_residues_2nd_chain;
         n_match_res_no_and_res_type = n_matches;
      }
      std::tuple<bool, float, float> get_match_fraction() const {
         bool c_status = false;
         bool f_1 = 0.0f;
         bool f_2 = 0.0f;
         if (n_residues_first_chain > 0) {
            if (n_residues_second_chain > 0) {
               f_1 = static_cast<float>(n_match_res_no_and_res_type)/static_cast<float>(n_residues_first_chain);
               f_2 = static_cast<float>(n_match_res_no_and_res_type)/static_cast<float>(n_residues_second_chain);
               c_status = true;
            }
         }
         return std::tuple<bool, float, float>(c_status, f_1, f_2);
      }
   };

   auto get_chains_overlap_info = [] (mmdb::Chain *chain_first, mmdb::Chain *chain_second, mmdb::Manager *mol, float big_overlap_fraction_limit) {

                                     bool status = false;
                                     mmdb::Contact *pscontact = NULL;
                                     int n_contacts;
                                     long i_contact_group = 1;
                                     mmdb::mat44 my_matt;
                                     for (int i=0; i<4; i++)
                                        for (int j=0; j<4; j++)
                                           my_matt[i][j] = 0.0;
                                     for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
                                     //

                                     std::string chain_first_chain_id(chain_first->GetChainID());
                                     std::string chain_second_chain_id(chain_second->GetChainID());
                                     mmdb::Atom **atom_selection_1 = 0; // member data - cleared on destruction
                                     mmdb::Atom **atom_selection_2 = 0; // member data - cleared on destruction
                                     int n_selected_atoms_1 = 0;
                                     int n_selected_atoms_2 = 0;
                                     int selhnd_1 = mol->NewSelection();
                                     int selhnd_2 = mol->NewSelection();
                                     mol->SelectAtoms(selhnd_1, 0, chain_first_chain_id.c_str(),
                                                      mmdb::ANY_RES, // starting resno, an int
                                                      "*", // any insertion code
                                                      mmdb::ANY_RES, // ending resno
                                                      "*", // ending insertion code
                                                      "*", // any residue name
                                                      "*", // atom name
                                                      "*", // elements
                                                      "");
                                     mol->SelectAtoms(selhnd_2, 0, chain_second_chain_id.c_str(),
                                                      mmdb::ANY_RES, // starting resno, an int
                                                      "*", // any insertion code
                                                      mmdb::ANY_RES, // ending resno
                                                      "*", // ending insertion code
                                                      "*", // any residue name
                                                      "*", // atom name
                                                      "*", // elements
                                                      "");
                                     mol->GetSelIndex(selhnd_1, atom_selection_1, n_selected_atoms_1);
                                     mol->GetSelIndex(selhnd_2, atom_selection_2, n_selected_atoms_2);

                                     if (n_selected_atoms_1 == 0) return chain_overap_info_t(false);
                                     if (n_selected_atoms_2 == 0) return chain_overap_info_t(false);

                                     mol->SeekContacts(atom_selection_1, n_selected_atoms_1,
                                                       atom_selection_2, n_selected_atoms_2,
                                                       0.0f, 1.0f,
                                                       1,        // seqDist 0 -> in same res also
                                                       pscontact, n_contacts,
                                                       0, &my_matt, i_contact_group);

                                     unsigned int n_match_res_no_and_res_type = 0;

                                     if (n_contacts > 0) {
                                        float contact_fract_1 = static_cast<float>(n_contacts)/static_cast<float>(n_selected_atoms_1);
                                        float contact_fract_2 = static_cast<float>(n_contacts)/static_cast<float>(n_selected_atoms_2);
                                        if (contact_fract_1 > big_overlap_fraction_limit) status = true;
                                        if (contact_fract_2 > big_overlap_fraction_limit) status = true;
                                        if (false)
                                           std::cout << "debug:: in chains_substantially_overlay() "
                                                     << chain_first_chain_id << " " << chain_second_chain_id << " "
                                                     << n_contacts << " " << n_selected_atoms_1 << " " << n_selected_atoms_2 << " "
                                                     << contact_fract_1 << " " <<  contact_fract_2 << std::endl;
                                        if (pscontact) {
                                           for (int i=0; i<n_contacts; i++) {
                                              mmdb::Atom *at_1 = atom_selection_1[pscontact[i].id1];
                                              mmdb::Atom *at_2 = atom_selection_2[pscontact[i].id2];
                                              int res_no_1 = at_1->GetSeqNum();
                                              int res_no_2 = at_2->GetSeqNum();
                                              if (res_no_2 == res_no_1) {
                                                 std::string res_type_1(at_1->GetResName());
                                                 std::string res_type_2(at_2->GetResName());
                                                 if (res_type_2 == res_type_1) {
                                                    std::string atom_name_1(at_1->GetAtomName());
                                                    std::string atom_name_2(at_2->GetAtomName());
                                                    if (atom_name_1 == " CA ")
                                                       if (atom_name_2 == " CA ")
                                                          n_match_res_no_and_res_type += 1;
                                                 }
                                              }
                                           }
                                        }
                                     }
                                     mol->DeleteSelection(selhnd_1);
                                     mol->DeleteSelection(selhnd_2);
                                     if (false)
                                        std::cout << "debug:: chains_substantially_overlay(): "
                                                  << chain_first << " " << chain_second << " "
                                                  << chain_first_chain_id << " " << chain_second_chain_id
                                                  << " " << status << std::endl;
                                     chain_overap_info_t coi(status);
                                     coi.set_match_info(chain_first->GetNumberOfResidues(), chain_second->GetNumberOfResidues(), n_match_res_no_and_res_type);
                                     return coi;
                                  };

   auto chains_substantially_overlap = [get_chains_overlap_info] (mmdb::Chain *chain_first, mmdb::Chain *chain_second, mmdb::Manager *mol,
                                                                  float big_overlap_fraction_limit) {
                                          chain_overap_info_t coi = get_chains_overlap_info(chain_first, chain_second, mol, big_overlap_fraction_limit);
                                          return coi.status;
                                       };

   // this version we deal with NCS. A fragment doesn't always want to try to be placed in the longest chain - it could be the NCS related
   // copy of this chain
   //
   auto merge_and_delete_sequenced_fragments = [get_chains, get_chains_overlap_info, insert_residues_from_other_chains_using_chains]
      (mmdb::Manager *mol, const std::map<std::string, unsigned int> &chain_id_to_fam_index) {

                                                  mmdb::Model *model_p = mol->GetModel(1);
                                                  if (! model_p) return;

                                                  // first convert the fam indexing
                                                  std::map<std::string, unsigned int>::const_iterator it;
                                                  std::map<unsigned int, std::vector<std::string> > fam_index_to_chain_ids;
                                                  for (it=chain_id_to_fam_index.begin(); it!=chain_id_to_fam_index.end(); ++it) {
                                                     const std::string &chain_id = it->first;
                                                     unsigned int fam_index = it->second;
                                                     fam_index_to_chain_ids[fam_index].push_back(chain_id);
                                                  }

                                                  // now sort the chain ids in those vectors by length (number of residues)
                                                  std::map<unsigned int, std::vector<std::string> >::iterator it_seq;
                                                  for (it_seq=fam_index_to_chain_ids.begin(); it_seq!=fam_index_to_chain_ids.end(); ++it_seq) {
                                                     std::vector<std::string> &chain_ids = it_seq->second;
                                                     std::sort(chain_ids.begin(), chain_ids.end(), chain_length_sorter(model_p));
                                                  }

                                                  for (it_seq=fam_index_to_chain_ids.begin(); it_seq!=fam_index_to_chain_ids.end(); ++it_seq) {
                                                     const std::vector<std::string> &chain_ids(it_seq->second);

                                                     // first we want to delete fragments - say there are 2 or more fragments that overlay to
                                                     // to some large extent - there are 2 cases:
                                                     // 1: the sequences match
                                                     // 2: the sequences don't match
                                                     //
                                                     //
                                                     // This is only considering chains that were sequnced from the same sequence.
                                                     // Models that have multiple sequences will need a different merge/delete algorithm

                                                     if (chain_ids.size() > 1) {
                                                        for (const auto &chain_id_1 : chain_ids) {
                                                           for (const auto &chain_id_2 : chain_ids) {
                                                              if (chain_id_1 != chain_id_2) {
                                                                 std::cout << "........................ debug " << chain_id_1 << " " << chain_id_2 << std::endl;
                                                                 std::pair<mmdb::Chain *, mmdb::Chain *> chain_pair = get_chains(chain_id_1, chain_id_2, mol);
                                                                 if (chain_pair.first && chain_pair.second) {
                                                                    float bof = 0.5; // big overlap fraction
                                                                    chain_overap_info_t coi = get_chains_overlap_info(chain_pair.first, chain_pair.second, mol, bof);
                                                                    std::cout << "........................ debug " << chain_id_1 << " " << chain_id_2
                                                                              << " coi status " << coi.status << std::endl;
                                                                    if (coi.status) {
                                                                       auto results = coi.get_match_fraction();
                                                                       if (std::get<0>(results)) {
                                                                          float frac_1 = std::get<1>(results);
                                                                          float frac_2 = std::get<2>(results);
                                                                          if (frac_1 > 0.0) {
                                                                             if (frac_2 > 0.0) {
                                                                                if (coi.n_match_res_no_and_res_type > 1) {
                                                                                   std::cout << "........................ debug " << chain_id_1 << " " << chain_id_2
                                                                                             << " insert_residues_from_other_chains_using_chains() " << std::endl;
                                                                                   insert_residues_from_other_chains_using_chains(mol, chain_pair.first, chain_pair.second);
                                                                                   std::cout << "Delete Chain " << chain_id_2 << std::endl;
                                                                                   model_p->DeleteChain(chain_id_2.c_str());
                                                                                   mol->FinishStructEdit();
                                                                                }
                                                                             }
                                                                          }

                                                                          if (coi.n_match_res_no_and_res_type == 0) {
                                                                             // OK, they overlapped but not overlapping residue numbers.
                                                                             // i.e. we have a disagreement about the trace
                                                                             std::cout << "........................ debug " << chain_id_1 << " " << chain_id_2
                                                                                       << " overlapping but disagree about directionn - delete chain "
                                                                                       << chain_id_2 << std::endl;
                                                                             model_p->DeleteChain(chain_id_2.c_str());
                                                                          }
                                                                       }
                                                                    }
                                                                 }
                                                              }
                                                           }
                                                        }
                                                     }
                                                  }
                                               };

   auto is_poly_ala_chain = [] (mmdb::Chain *chain_p) {
                               bool status = true;
                               int n_res = chain_p->GetNumberOfResidues();
                               for (int ires=0; ires<n_res; ires++) {
                                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                  if (residue_p) {
                                     std::string res_name(residue_p->GetResName());
                                     if (res_name != "ALA") {
                                        status = false;
                                        break;
                                     }
                                  }
                               }
                               return status;
                            };

   auto delete_poly_ala_chains_that_overlap_on_sequenced_chains = [is_poly_ala_chain, chains_substantially_overlap] (mmdb::Manager *mol,
                                                                                                                     float big_overlap_fraction_limit) {

                                                                     int imodel = 1;
                                                                     mmdb::Model *model_p = mol->GetModel(imodel);
                                                                     if (model_p) {

                                                                        bool continue_looping = false;
                                                                        do {
                                                                           continue_looping = false;
                                                                           int n_chains = model_p->GetNumberOfChains();
                                                                           for (int ichain=0; ichain<n_chains; ichain++) {
                                                                              mmdb::Chain *i_chain_p = model_p->GetChain(ichain);
                                                                              std::string chain_id(i_chain_p->GetChainID());
                                                                              bool poly_ala = is_poly_ala_chain(i_chain_p);
                                                                              if (poly_ala) {
                                                                                 for (int jchain=0; jchain<n_chains; jchain++) {
                                                                                    if (ichain != jchain) {
                                                                                       mmdb::Chain *j_chain_p = model_p->GetChain(jchain);
                                                                                       if (chains_substantially_overlap(i_chain_p, j_chain_p, mol,
                                                                                                                        big_overlap_fraction_limit)) {
                                                                                          std::cout << "INFO:: delete poly-ala chain " << chain_id << std::endl;
                                                                                          model_p->DeleteChain(ichain);
                                                                                          mol->FinishStructEdit();
                                                                                          continue_looping = true;
                                                                                          break;
                                                                                       }
                                                                                    }
                                                                                 }
                                                                              }
                                                                              if (continue_looping) break;
                                                                           }
                                                                        } while (continue_looping);
                                                                     }
                                                                  };


   auto score_atom_positions = [] (const std::vector<clipper::Coord_orth> &atom_positions,
                                   const clipper::Coord_orth &CA_this,
                                   const clipper::Coord_orth &CA_next,
                                   const clipper::Xmap<float> &xmap) {

                                  // Add peptide fingerprinting!

                                  float sum = 0.0f;
                                  for (const auto &pos : atom_positions)
                                     sum += coot::util::density_at_point(xmap, pos);
                                  return sum;
                               };

   auto apply_pepflips = [score_atom_positions] (mmdb::Manager *mol, const clipper::Xmap<float> &xmap) {
                            int imod = 1;
                            mmdb::Model *model_p = mol->GetModel(imod);
                            if (model_p) {
                               int n_chains = model_p->GetNumberOfChains();
                               for (int ichain=0; ichain<n_chains; ichain++) {
                                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                  int n_res = chain_p->GetNumberOfResidues();
                                  for (int ires=0; ires<(n_res-1); ires++) {
                                     mmdb::Residue *residue_this_p = chain_p->GetResidue(ires);
                                     mmdb::Residue *residue_next_p = chain_p->GetResidue(ires+1);
                                     if (residue_this_p) {
                                        if (residue_next_p) {
                                           mmdb::Atom *at_CA_this = residue_this_p->GetAtom(" CA ");
                                           mmdb::Atom *at_C_this  = residue_this_p->GetAtom(" C  ");
                                           mmdb::Atom *at_O_this  = residue_this_p->GetAtom(" O  ");
                                           mmdb::Atom *at_N_next  = residue_next_p->GetAtom(" N  ");
                                           mmdb::Atom *at_CA_next = residue_next_p->GetAtom(" CA ");
                                           if (at_CA_this && at_C_this && at_O_this && at_N_next && at_CA_next) {
                                              clipper::Coord_orth pos_CA_this = coot::co(at_CA_this);
                                              clipper::Coord_orth pos_C_this  = coot::co(at_C_this);
                                              clipper::Coord_orth pos_O_this  = coot::co(at_O_this);
                                              clipper::Coord_orth pos_N_next  = coot::co(at_N_next);
                                              clipper::Coord_orth pos_CA_next = coot::co(at_CA_next);

                                              std::vector<clipper::Coord_orth> atom_positions;
                                              atom_positions.push_back(pos_C_this);
                                              atom_positions.push_back(pos_O_this);
                                              atom_positions.push_back(pos_N_next);
                                              float d_current = score_atom_positions(atom_positions, pos_CA_this, pos_CA_next, xmap);

                                              // spin C,O,N position 180 around CA-CA
                                              //
                                              clipper::Coord_orth axis = pos_CA_next - pos_CA_this;
                                              pos_C_this = coot::util::rotate_around_vector(axis, pos_C_this, pos_CA_this, M_PI);
                                              pos_O_this = coot::util::rotate_around_vector(axis, pos_O_this, pos_CA_this, M_PI);
                                              pos_N_next = coot::util::rotate_around_vector(axis, pos_N_next, pos_CA_this, M_PI);

                                              atom_positions.clear();
                                              atom_positions.push_back(pos_C_this);
                                              atom_positions.push_back(pos_O_this);
                                              atom_positions.push_back(pos_N_next);
                                              float d_spun = score_atom_positions(atom_positions, pos_CA_this, pos_CA_next, xmap);

                                              if (d_spun > d_current) {
                                                 at_C_this->x = pos_C_this.x(); at_C_this->y = pos_C_this.y(); at_C_this->z = pos_C_this.z();
                                                 at_O_this->x = pos_O_this.x(); at_O_this->y = pos_O_this.y(); at_O_this->z = pos_O_this.z();
                                                 at_N_next->x = pos_N_next.x(); at_N_next->y = pos_N_next.y(); at_N_next->z = pos_N_next.z();
                                              }
                                           }
                                        }
                                     }
                                  }
                               }
                            }
                         };

   // this could/should be a general-purpose function
   //
   auto rsr_molecule = [refine_isolated_chain] (mmdb::Manager *mol, const clipper::Xmap<float> &xmap, const coot::protein_geometry &geom,
                           ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight) {
                          int imod = 1;
                          mmdb::Model *model_p = mol->GetModel(imod);
                          if (model_p) {
                             int n_chains = model_p->GetNumberOfChains();
                             for (int ichain=0; ichain<n_chains; ichain++) {
                                mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                refine_isolated_chain(chain_p, mol, geom, thread_pool_p, n_threads, weight);
                             }
                          }
                       };

   // ---------------------------------------------------------------------------------------

   coot::protein_geometry geom;
   geom.init_standard();
   bool cryo_em_peptide_seeds_mode = false;

   coot::minimol::molecule flood_molecule = get_flood_molecule(xmap, rmsd_cut_off_for_flood, flood_atom_mask_radius);


   // hack for testing tutorial modern
   // mean = clipper::Coord_orth(60, 5, 12);
   // mean = clipper::Coord_orth(0, 20, 19); // 1gwd
   // mean = clipper::Coord_orth(109, 109, 110); // emd-22898

   // clipper::Coord_orth hack_centre = clipper::Coord_orth(109, 109, 110); // emd-22898
   // clipper::Coord_orth hack_centre = clipper::Coord_orth(40, 9, 13); // tutorial modern
   // clipper::Coord_orth hack_centre = clipper::Coord_orth(0, 20, 19); // 1gwd
   // clipper::Coord_orth hack_centre = clipper::Coord_orth(109, 109, 110); // emd-22898
   //bool use_hack_centre = true;

   // this mol is mostly filtered, but can have some overlapping fragments
   // it has mainchain + CBs. It is not refined.
   std::pair<bool, clipper::Coord_orth> hack_centre(true, clipper::Coord_orth(0, 20, 19));
   mmdb::Manager *mol = find_connected_fragments(flood_molecule, xmap, variation, n_top_spin_pairs, n_top_fragments, cryo_em_peptide_seeds_mode, hack_centre);

   auto update_working_mol = [get_mol_edit_lock, release_mol_edit_lock] (watch_res_tracer_data_t *watch_res_tracer_data_p, mmdb::Manager *mol_to_copy,
                                                                         const std::string &locker) mutable {
                                if (watch_res_tracer_data_p) {
                                   get_mol_edit_lock(locker);
                                   watch_res_tracer_data_p->working_mol->Copy(mol_to_copy, mmdb::MMDBFCM_All);
                                   watch_res_tracer_data_p->update_flag = true;
                                   watch_res_tracer_data_p->update_count++;
                                   release_mol_edit_lock(locker);
                                }
                             };

   if (mol) {

      update_working_mol(watch_res_tracer_data_p, mol, "update-1");
      mol->WritePDBASCII("stage-3-post-find-connected-fragments.pdb");
      unsigned int update_count = 0; // functions that modify mol, update this count

      unsigned int n_threads = coot::get_max_number_of_threads();
      ctpl::thread_pool thread_pool(n_threads);
      auto tp_0 = std::chrono::high_resolution_clock::now();

      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         int n_chains_post_find_connected_fragments = model_p->GetNumberOfChains();
#if 0
         for (int ichain=0; ichain<n_chains_post_find_connected_fragments; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::cout << "Refine chain " << chain_p->GetChainID() << std::endl;
            refine_chain(chain_p, mol, geom, &thread_pool, n_threads, weight);
         }
#endif
#if 1
         std::cout << "multi-refine chains: " << n_chains_post_find_connected_fragments << " chains..." << std::endl;
         multi_refine_individual_chains(mol, geom, &thread_pool, n_threads, weight);
#endif
         auto tp_1 = std::chrono::high_resolution_clock::now();
         auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
         std::cout << "Timings: proc(): for the initial RSR of the " << n_chains_post_find_connected_fragments << " chains:          "
                   << d10 << " milliseconds" << std::endl;
         update_working_mol(watch_res_tracer_data_p, mol, "update-2");

         mol->WritePDBASCII("stage-4-post-refine-of-connected-fragments.pdb");

         density_and_omega_based_trim_chain_terminii(mol, xmap);

         update_working_mol(watch_res_tracer_data_p, mol, "update-3");

         std::cout << "PROGRESS:: density_and_omega_based_trim_chain_terminii() done " << std::endl;

         mol->WritePDBASCII("stage-5-post-omega-and-density-check-trim.pdb");

         auto write_chains = [write_chain] (const std::string &file_name_stub, mmdb::Manager *mol) {
                                mmdb::Model *model_p = mol->GetModel(1);
                                int n_chains = model_p->GetNumberOfChains();
                                for (int ichain=0; ichain<n_chains; ichain++) {
                                   mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                   std::string chain_id(chain_p->GetChainID());
                                   std::cout << "Write " << file_name_stub << " "<< chain_id << std::endl;
                                   std::string file_name = file_name_stub + "-" + chain_id + ".pdb";
                                   write_chain(chain_p, mol, file_name);
                                }
                             };

         if (cryo_em_peptide_seeds_mode) {

            // I can filter the chains of mol here I think
            // cluster_filter_peptide_seeds(mol, map); cluster on the sum of the distances between the CA positions?
            // I want to see that histogram!

            rama_rsr_extend_fragments(mol, xmap, &thread_pool, n_threads, weight, n_phi_psi_trials, geom, &update_count);
            write_chains("post-extensions", mol);
            mol->WritePDBASCII("c-stage-5-post-extensions.pdb");
            int n_res_min = 3;
            delete_chains_that_are_too_short(mol, n_res_min);
            mol->WritePDBASCII("c-stage-5b-post-delete-short-fragments.pdb");
            coot::merge_C_and_N_terminii_0_gap(mol);
            mol->WritePDBASCII("c-stage-6-post-merge-C_and_N_terminii.pdb");
            std::map<std::string, unsigned int> chain_id_to_fam_index_map = apply_sequence_to_fragments(mol, xmap, fam, geom);
            mol->WritePDBASCII("c-stage-7-post-sequence-fragments.pdb");
            float big_overlap_fraction_limit_for_poly_ala = 0.2;
            delete_poly_ala_chains_that_overlap_on_sequenced_chains(mol, big_overlap_fraction_limit_for_poly_ala);
            merge_and_delete_sequenced_fragments(mol, chain_id_to_fam_index_map);
            n_chains_post_find_connected_fragments = model_p->GetNumberOfChains();
            std::cout << "multi-refine chains: " << n_chains_post_find_connected_fragments << " chains..." << std::endl;
            multi_refine_individual_chains(mol, geom, &thread_pool, n_threads, weight);
            mol->WritePDBASCII("c-stage-8-post-merge-sequenced-fragments.pdb");
            bring_together_consecutive_C_and_N_by_symmetry_transformation(mol, xmap);
            mol->WritePDBASCII("c-stage-9-post-sequencing-bring-together-consecutive-Cs-and-Ns.pdb");

         } else {

            int n_res_min = 2;
            delete_chains_that_are_too_short(mol, n_res_min);

            update_working_mol(watch_res_tracer_data_p, mol, "update-4");

            mol->WritePDBASCII("stage-6-post-delete-chains-that-are-too-short.pdb");

            if (false) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  std::string chain_id(chain_p->GetChainID());
                  std::cout << "Write post-trim chain " << chain_id << std::endl;
                  std::string file_name = "post-trim-chain-" + chain_id + ".pdb";
                  write_chain(chain_p, mol, file_name);
               }
            }

            // now filter out close chains again, after refinement chains will have/may have converged.
            float big_overlap_fraction_limit = 0.89;
            std::map<std::string, double> scores_for_chain_ids = make_density_fit_scores_for_chains(mol, xmap);
            auto deletable_chains_map = find_chains_that_overlap_other_chains(mol, big_overlap_fraction_limit, scores_for_chain_ids);
            filter_similar_chains(mol, deletable_chains_map); // don't think - just do.

            update_working_mol(watch_res_tracer_data_p, mol, "update-5");

            if (false) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  std::string chain_id(chain_p->GetChainID());
                  std::cout << "Write post-filter-2 chain " << chain_id << std::endl;
                  std::string file_name = "post-filter-2-chain-" + chain_id + ".pdb";
                  write_chain(chain_p, mol, file_name);
               }
            }

            mol->WritePDBASCII("stage-7-post-filter-2-scoring-real-atoms.pdb");

            unsigned int n_twisted_peptides_per_chain_max = 1; // delete chains with more than one
            double omega_twist_max = 33.3; // max degrees from trans
            delete_chains_that_have_twisted_trans_peptides(mol, n_twisted_peptides_per_chain_max, omega_twist_max);

            update_working_mol(watch_res_tracer_data_p, mol, "update-6");

            mol->WritePDBASCII("stage-8-post-twisted-peptide-chain-filter.pdb");

            if (false) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  std::string chain_id(chain_p->GetChainID());
                  std::cout << "Write post-twisted-peptide-chain-filter chain " << chain_id << std::endl;
                  std::string file_name = "post-twisted-peptide-chain-filter-chain-" + chain_id + ".cif";
                  write_chain(chain_p, mol, file_name, true);
               }
            }

            auto tp_3 = std::chrono::high_resolution_clock::now();
            rama_rsr_extend_fragments(mol, xmap, &thread_pool, n_threads, weight, n_phi_psi_trials, geom, &update_count);
            auto tp_4 = std::chrono::high_resolution_clock::now();

            update_working_mol(watch_res_tracer_data_p, mol, "update-7");

            // this is the right thing to do generally of course, but not when I am trying to debug tracing.
            if (false)
               coot::renumber_chains_start_at_least_at_1(mol); // make this verbose, so I can track the residue numbers?

            auto tp_5 = std::chrono::high_resolution_clock::now();
            mol->WritePDBASCII("stage-9-post-extensions.pdb");

            coot::merge_C_and_N_terminii_0_gap(mol);

            mol->WritePDBASCII("stage-10-post-merge-C_and_N_terminii.pdb");

            auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
            std::cout << "Timings: proc(): Rama-RSR chain extensions: " << d43 << " milliseconds" << std::endl;

            // now I want to merge fragments (non-trivial)
            //
            // first I want a sorted list of mergeable chains (sort by distance of terminii)

            std::map<std::string, unsigned int> chain_id_to_fam_index_map = apply_sequence_to_fragments(mol, xmap, fam, geom);

            update_working_mol(watch_res_tracer_data_p, mol, "update-8");

            mol->WritePDBASCII("stage-11-post-sequence-fragments.pdb");
            mol->WriteCIFASCII("stage-11-post-sequence-fragments.cif");

            if (true) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  std::string chain_id(chain_p->GetChainID());
                  std::cout << "Write post-sequence-fragment chain " << chain_id << std::endl;
                  std::string file_name = "post-sequence-fragment-chain-" + chain_id + ".cif";
                  write_chain(chain_p, mol, file_name, true);
               }
            }

            float big_overlap_fraction_limit_for_poly_ala = 0.2;
            delete_poly_ala_chains_that_overlap_on_sequenced_chains(mol, big_overlap_fraction_limit_for_poly_ala);

            update_working_mol(watch_res_tracer_data_p, mol, "update-9");

            mol->WritePDBASCII("stage-12-post-delete-overlapping-poly-ala.pdb");
            mol->WriteCIFASCII("stage-12-post-delete-overlapping-poly-ala.cif");

            if (true) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  std::string chain_id(chain_p->GetChainID());
                  std::cout << "Write post-delete-overlapping-poly-ala " << chain_id << std::endl;
                  std::string file_name = "post-delete-overlapping-poly-ala-" + chain_id + ".cif";
                  write_chain(chain_p, mol, file_name, true);
               }
            }

            std::cout << "::::::::::::::::: merge and delete sequenced fragments into single chain " << std::endl;
            if (! with_ncs)
               merge_and_delete_sequenced_fragments_into_single_chain(mol, chain_id_to_fam_index_map);
            else
               merge_and_delete_sequenced_fragments(mol, chain_id_to_fam_index_map);

            update_working_mol(watch_res_tracer_data_p, mol, "update-10");

            mol->WritePDBASCII("stage-13-post-merge-sequence-fragments.pdb");

            if (false) {
               // print out the residue numbers so that I can see that they are sane
               // (they are)
               std::cout << "Sanity check" << std::endl;
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int n_res = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<n_res; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        std::cout << "   " << coot::residue_spec_t(residue_p) << std::endl;
                     }
                  }
               }
            }

            // Look for big gaps between N and C of consecutive residues - in those cases there is a fragment
            // drawn in a different asymmetric unit.

            bring_together_consecutive_C_and_N_by_symmetry_transformation(mol, xmap);

            update_working_mol(watch_res_tracer_data_p, mol, "update-11");

            mol->WritePDBASCII("stage-14-post-sequencing-bring-together-consecutive-Cs-and-Ns.pdb");

            // delete_poly_ala_chains_that_overlap_on_sequenced_chains(mol);

            // mol->WritePDBASCII("stage-14-post-delete-overlapping-poly-als.pdb");

            std::cout << "backrub-1" << std::endl;
            coot::backrub_molecule(mol, &xmap, geom);
            std::cout << "backrub-2" << std::endl;
            coot::backrub_molecule(mol, &xmap, geom);
            //
            mol->WritePDBASCII("stage-15-post-backrub-chain.pdb");

            update_working_mol(watch_res_tracer_data_p, mol, "update-12");

            // resolve clashing side-chains by deletion?

            // add calculation of difference map and pepflips here. - maybe a difference map is not needed
            // Just a fingerprinted flip search - which should be fast.

            std::cout << "pepflips" << std::endl;
            apply_pepflips(mol, xmap);

            update_working_mol(watch_res_tracer_data_p, mol, "update-13");

            mol->WritePDBASCII("stage-16-post-pepflips.pdb");

            std::cout << "RSR" << std::endl;
            coot::util::pdbcleanup_serial_residue_numbers(mol);
            rsr_molecule(mol, xmap, geom, &thread_pool, n_threads, weight); // chain by chain actually

            std::cout << "debug:: in res_tracer_proc() update_working_mol post-rsr" << std::endl;
            update_working_mol(watch_res_tracer_data_p, mol, "update-14");

            mol->WritePDBASCII("stage-17-post-rsr.pdb");
         }
      }

      std::cout << "--- done proc() ---" << std::endl;
   }
   watch_res_tracer_data_p->finished = true;

}
