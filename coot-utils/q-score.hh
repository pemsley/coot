#ifndef SIMPLE_ATOM_NEIGHBOURS_HH
#define SIMPLE_ATOM_NEIGHBOURS_HH

#include <vector>
#include <map>
#include <chrono>
#include <thread>
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>

#include "analysis/stats.hh"
#include "utils/split-indices.hh"
#include "coot-coord-utils.hh"
#include "coot-map-utils.hh"
#include "coot-density-stats.hh"

namespace coot {

   class q_score_t {

      mmdb::Manager *mol;
      std::map<mmdb::Atom *, std::vector<clipper::Coord_orth> > atom_neighbours;
      int sel_hnd;

      bool is_closer_to_neighbour(const clipper::Coord_orth &pos, float radius,
                                  const std::vector<clipper::Coord_orth> &neighb_pos) const {

         bool status = false;
         const double d_crit = 2.0;
         const double d2_crit_sqd = d_crit * d_crit;
         std::vector<clipper::Coord_orth>::const_iterator it;
         for (it=neighb_pos.begin(); it!=neighb_pos.end(); ++it) {
            double ll = (pos - *it).lengthsq();
            if (ll < d2_crit_sqd) {
               if (ll < radius * radius) {
                  status = true;
                  break;
               }
            }
         }
         return status;
      }

      float calc_model_value(const float &radius, const float &map_mean, const float &map_rmsd) const {
         float A = map_mean + 10.0 * map_rmsd;
         float B = map_mean -  1.0 * map_rmsd;
         float sigma = 0.6; // Angstroms
         float z = radius/sigma;
         float model = A * std::exp(-0.5 * z*z) + B;
         return model;
      }

      class atom_q_score_info_t {
      public:
         bool is_valid;
         double q_score;
         unsigned int n_points;
         atom_q_score_info_t(bool status, double q, unsigned int n) : is_valid(status), q_score(q), n_points(n) {}
      };

      double squared(const double &f) const { return f * f; }

      // first is if the return value is valid.
      atom_q_score_info_t calc_q_score(const std::vector<std::pair<float, float> > &u_v_pairs) const {
         stats::single u_stats;
         stats::single v_stats;
         std::vector<std::pair<float, float> >::const_iterator it;
         for (it=u_v_pairs.begin(); it!=u_v_pairs.end(); ++it) {
            u_stats.add(it->first);
            v_stats.add(it->second);
         }
         double u_mean = u_stats.mean();
         double v_mean = v_stats.mean();
         double running_top = 0.0;
         double running_delta_u = 0.0;
         double running_delta_v = 0.0;
         unsigned int n_points = u_v_pairs.size();
         for (it=u_v_pairs.begin(); it!=u_v_pairs.end(); ++it) {
            running_top += (it->first - u_mean) * (it->second - v_mean);
            running_delta_u += squared(it->first  - u_mean);
            running_delta_v += squared(it->second - v_mean);
         }
         // double av_delta_u = running_delta_u/static_cast<double>(u_v_pairs.size());
         // double av_delta_v = running_delta_v/static_cast<double>(u_v_pairs.size());
         double av_top = running_top/static_cast<double>(u_v_pairs.size());
         double b1 = running_delta_u/static_cast<double>(u_v_pairs.size());
         double b2 = running_delta_v/static_cast<double>(u_v_pairs.size());
         double q_score = av_top/(std::sqrt(b1) * std::sqrt(b2));
         bool status = true;
         return atom_q_score_info_t(status, q_score, n_points);
      }

   public:
      explicit q_score_t(mmdb::Manager *mol) : mol(mol) {
         if (mol) {
            float max_dist = 4.0;

            int imod = 1;
            sel_hnd = mol->NewSelection(); // d
            mol->SelectAtoms (sel_hnd, imod, "*",
                              mmdb::ANY_RES, // starting resno, an int
                              "*", // any insertion code
                              mmdb::ANY_RES, // ending resno
                              "*", // ending insertion code
                              "*", // any residue name
                              "*", // atom name
                              "*", // elements
                              "*"  // alt loc.
                              );
            mmdb::PPAtom atom_selection = 0;
            int n_selected_atoms = 0;
            mol->GetSelIndex(sel_hnd, atom_selection, n_selected_atoms);

            mmdb::Contact *pscontact = NULL;
            int n_contacts;
            long i_contact_group = 1;
            mmdb::mat44 my_matt;
            mmdb::SymOps symm;
            for (int i=0; i<4; i++)
               for (int j=0; j<4; j++)
                  my_matt[i][j] = 0.0;
            for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

            mol->SeekContacts(atom_selection, n_selected_atoms,
                              atom_selection, n_selected_atoms,
                              0, max_dist,
                              0, // in same residue
                              pscontact, n_contacts,
                              0, &my_matt, i_contact_group); // makes reverses also
            if (n_contacts > 0) {
               if (pscontact) {
                  for (int i=0; i<n_contacts; i++) {
                     mmdb::Atom *at        = atom_selection[pscontact[i].id1];
                     mmdb::Atom *neighb_at = atom_selection[pscontact[i].id2];
                     clipper::Coord_orth neighb_pos = coot::co(neighb_at);
                     atom_neighbours[at].push_back(neighb_pos);
                  }
               }
            }
         }
      }

      void close() {
         mol->DeleteSelection(sel_hnd);
      }

      void score_atom_range(int atom_idx_begin, int atom_idx_end,
                            mmdb::Atom **atom_selection,
                            const clipper::Xmap<float> &xmap,
                            float map_mean, float map_rmsd,
                            const std::vector<clipper::Coord_orth>  &directions,
                            int n_radius, float radius_max,
                            std::vector<float> *q_score_results_p) {

         std::vector<clipper::Coord_orth> empty_v; // microwave ovens.
         for (int iat=atom_idx_begin; iat!=atom_idx_end; iat++) {
            mmdb::Atom *at = atom_selection[iat];
            if (!at->isTer()) {
               std::vector<clipper::Coord_orth> &v = empty_v;
               std::map<mmdb::Atom *, std::vector<clipper::Coord_orth> >::const_iterator it =
                  atom_neighbours.find(at);
               if (it != atom_neighbours.end()) v = it->second;
               clipper::Coord_orth atom_pos = coot::co(at);
               std::vector<std::pair<float, float> > u_v_pairs;

               for (int i_rad=0; i_rad<n_radius; i_rad++) {
                  float radius = static_cast<float>(i_rad) * radius_max / static_cast<float>(n_radius);
                  float model_value = calc_model_value(radius, map_mean, map_rmsd);
                  for (const auto &dir : directions) {
                     clipper::Coord_orth pt = atom_pos + radius * dir;
                     if (! is_closer_to_neighbour(pt, radius, v)) {
                        float d = coot::util::density_at_point_by_linear_interpolation(xmap, pt);
                        std::pair<float, float> u_v(d, model_value);
                        u_v_pairs.push_back(u_v);
                     }
                  }
               }
               if (! u_v_pairs.empty()) {
                  atom_q_score_info_t atom_q_score = calc_q_score(u_v_pairs);
                  if (atom_q_score.is_valid && atom_q_score.n_points >= 4) {
                     double q_score = atom_q_score.q_score;
                     q_score_results_p->at(iat) = q_score;

                     if (false) {
                        // debug using density stats:
                        util::density_correlation_stats_info_t dcsi;
                        for (const auto &uv : u_v_pairs) dcsi.add(uv.first, uv.second);
                        double cc = dcsi.correlation();

                        // 100 ms for 7vvl
                        if (false)
                           std::cout << "at: "
                                     << at->GetAtomName() << " " << at->GetResName() << " " << at->GetChainID() << " " << at->GetSeqNum()
                                     << " q-score " << q_score << " cc " << cc << std::endl;

                     }
                  }
               }
            }
         }
      }

      void calc(const clipper::Xmap<float> &xmap, float map_mean, float map_variance) {

         int n_radius = 21;
         double radius_max = 2.1;

         std::vector<clipper::Coord_orth> directions =
            { clipper::Coord_orth(-1, -1, -1), clipper::Coord_orth(-1, -1,  1),
              clipper::Coord_orth(-1,  1, -1), clipper::Coord_orth(-1,  1,  1),
              clipper::Coord_orth( 1, -1, -1), clipper::Coord_orth( 1, -1,  1),
              clipper::Coord_orth( 1,  1, -1), clipper::Coord_orth( 1,  1,  1)};
         for (auto &dir : directions)
            dir = (1.0/sqrt(3.0)) * dir;

         float map_rmsd = std::sqrt(map_variance);
         mmdb::PPAtom atom_selection = 0;
         int n_selected_atoms = 0;
         mol->GetSelIndex(sel_hnd, atom_selection, n_selected_atoms);

         std::vector<float> q_score_results(n_selected_atoms, -1111.1);
         unsigned int n_threads = 6;
         std::vector<std::pair<unsigned int, unsigned int> > airs =
            atom_index_ranges(n_selected_atoms, n_threads);

         std::vector<std::thread> threads; // non-trivial to use threading atm
         // fill q_score_results
         for (const auto &air : airs)
            score_atom_range(air.first, air.second, atom_selection, xmap, map_mean, map_rmsd,
                             directions, n_radius, radius_max, &q_score_results);

         // for (unsigned int i_thread=0; i_thread<threads.size(); i_thread++) threads[i_thread].join();

         // Now batch the atom-based results into residues

         // auto tp_0 = std::chrono::high_resolution_clock::now();

         // add a UDD for every atom:
         int udd_q_score = mol->RegisterUDReal(mmdb::UDR_ATOM, "Q Score");
         std::map<residue_spec_t, stats::single> residue_q_scores;
         for (int iat=0; iat<n_selected_atoms; iat++) {
            mmdb:: Atom *at = atom_selection[iat];
            residue_spec_t res_spec(at->GetResidue());
            float q_score = q_score_results[iat];
            at->PutUDData(udd_q_score, q_score);
            if (false)
               std::cout << "results per atom " << atom_spec_t(at) << " " << q_score << " B-factor " << at->tempFactor
                         << std::endl;
            if (q_score > -1000.0)
               residue_q_scores[res_spec].add(q_score);
         }
         std::map<residue_spec_t, stats::single>::const_iterator it;
         if (false)
            for (it=residue_q_scores.begin(); it!=residue_q_scores.end(); ++it)
               std::cout << "   " << it->first.chain_id << " " << it->first.res_no
                         << " n-atoms " << it->second.size() << " mean " << it->second.mean() << std::endl;

         // auto tp_1 = std::chrono::high_resolution_clock::now();
         // auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
         // std::cout << "# Timing-for-consolidation " << d10 << " milliseconds" << std::endl;

      }
   };

}

#endif // SIMPLE_ATOM_NEIGHBOURS_HH
