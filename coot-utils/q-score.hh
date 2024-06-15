#ifndef SIMPLE_ATOM_NEIGHBOURS_HH
#define SIMPLE_ATOM_NEIGHBOURS_HH

#include <vector>
#include <map>
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>

#include "analysis/stats.hh"
#include "coot-coord-utils.hh"
#include "coot-map-utils.hh"
#include "coot-density-stats.hh"

namespace coot {

   class q_score_t {

      mmdb::Manager *mol;
      std::map<mmdb::Atom *, std::vector<clipper::Coord_orth> > atom_neighbours;
      int sel_hnd;

      bool is_closer_to_neighbour(const clipper::Coord_orth &pos, float radius, const std::vector<clipper::Coord_orth> &neighb_pos) const {

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

      double calc_q_score(const std::vector<std::pair<float, float> > &u_v_pairs) const {
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
         for (it=u_v_pairs.begin(); it!=u_v_pairs.end(); ++it) {
            running_top += (it->first - u_mean) * (it->second - v_mean);
            running_delta_u += fabs(it->first  - u_mean);
            running_delta_v += fabs(it->second - v_mean);
         }
         double av_delta_u = running_delta_u/static_cast<double>(u_v_pairs.size());
         double av_delta_v = running_delta_v/static_cast<double>(u_v_pairs.size());
         double av_top = running_top/static_cast<double>(u_v_pairs.size());
         double q_score = av_top/(av_delta_u * av_delta_v);
         return q_score;
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
                     const int &id2 = pscontact[i].id2;
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

      void calc(clipper::Xmap<float> &xmap, float map_mean, float map_variance) {

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
         // std::map<mmdb::Atom *, std::vector<clipper::Coord_orth> >::const_iterator it;
         // for (it=atom_neighbours.begin(); it!=atom_neighbours.end(); ++it) {
         //   mmdb::Atom *at = it->first;
         //   const std::vector<clipper::Coord_orth> &v = it->second;
         mmdb::PPAtom atom_selection = 0;
         int n_selected_atoms = 0;
         std::vector<clipper::Coord_orth> empty_v; // microwave ovens.
         mol->GetSelIndex(sel_hnd, atom_selection, n_selected_atoms);
         for (int iat=0; iat<n_selected_atoms; iat++) {
            mmdb::Atom *at = atom_selection[iat];
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
            double q_score = calc_q_score(u_v_pairs);

            // debug using density stats:
            util::density_correlation_stats_info_t dcsi;
            for (const auto &uv : u_v_pairs) dcsi.add(uv.first, uv.second);
            double cc = dcsi.correlation();

            // 100 ms for 7vvl
            std::cout << "at: " << at->GetAtomName() << " " << at->GetResName() << " " << at->GetChainID() << " " << at->GetSeqNum()
                      << " q-score " << q_score << " cc " << cc << std::endl;

         }
      }
   };

}

#endif // SIMPLE_ATOM_NEIGHBOURS_HH
