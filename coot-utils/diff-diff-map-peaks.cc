
#include <clipper/core/map_interp.h>
#include <clipper/contrib/skeleton.h>

#include "diff-diff-map-peaks.hh"

std::vector<std::pair<clipper::Coord_orth, float>> coot::diff_diff_map_peaks(const clipper::Xmap<float> &m1,
                                                                             const clipper::Xmap<float> &m2,
                                                                             float base_level) {
   std::vector<std::pair<clipper::Coord_orth, float> > v_plus;
   std::vector<std::pair<clipper::Coord_orth, float> > v_neg;

   clipper::Skeleton_basic::Neighbours neighb(m1, 0.25, 1.75); // 3x3x3 cube, not centre

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = m1.first(); !ix.last(); ix.next())  { // iterator index.
      float v1 = m1[ix];
      float v2 = m2[ix];
      float d = v2 - v1;
      // std::cout << ix.coord().format() << " d " << d << std::endl;
      if (d > base_level) {
         if (fabsf(v2) < fabsf(v1)) {
            bool is_peak = true; // unless we find a higher neighboour
            for (int i=0; i<neighb.size(); i++) {
               clipper::Coord_grid c_g = ix.coord() + neighb[i];
               float n1 = m1.get_data(c_g);
               float n2 = m2.get_data(c_g);
               float nd = n2 - n1;
               if (d < nd) {
                  is_peak = false;
                  break;
               }
            }
            if (is_peak) {
               // std::cout << ix.coord().format() << " d pos " << d << std::endl;
               clipper::Coord_frac cf = ix.coord().coord_frac(m1.grid_sampling());
               clipper::Coord_orth co = cf.coord_orth(m1.cell());
               std::pair p = std::make_pair(co, d);
               v_plus.push_back(p);
            }
         }
      }

      if (d < -base_level) {
         if (fabsf(v2) < fabsf(v1)) {
            bool is_peak = true; // negative peak, that is.
            for (int i=0; i<neighb.size(); i++) {
               clipper::Coord_grid c_g = ix.coord() + neighb[i];
               float n1 = m1.get_data(c_g);
               float n2 = m2.get_data(c_g);
               float nd = n2 - n1;
               if (d > nd) {
                  is_peak = false;
                  break;
               }
            }
            if (is_peak) {
               // std::cout << ix.coord().format() << " d neg " << d << std::endl;
               clipper::Coord_frac cf = ix.coord().coord_frac(m1.grid_sampling());
               clipper::Coord_orth co = cf.coord_orth(m1.cell());
               std::pair p = std::make_pair(co, d);
               v_plus.push_back(p);
            }
         }
      }
   }
   return v_plus;
}

#include "coot-coord-utils.hh"

std::vector<std::pair<clipper::Coord_orth, float> >
coot::move_peaks_to_around_position(const clipper::Coord_orth &screen_centre,
                                const clipper::Spacegroup &spacegroup,
                                const clipper::Cell &cell,
                                const std::vector<std::pair<clipper::Coord_orth, float> > &v_in) {

   std::vector<std::pair<clipper::Coord_orth, float> > v = v_in;
   int n = spacegroup.num_symops();
   for (unsigned int i = 0; i < v_in.size(); i++) {
      const auto &pos = v_in[i].first;
      double d_best = std::sqrt((pos-screen_centre).lengthsq()); // current
      bool improved = false;
      clipper::RTop_orth rtop_best(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), clipper::Coord_orth(0,0,0));
      for (int isym=0; isym<n; isym++) {
         for (int x_shift = -1; x_shift<2; x_shift++) {
            for (int y_shift = -1; y_shift<2; y_shift++) {
               for (int z_shift = -1; z_shift<2; z_shift++) {
                  clipper::Coord_frac cell_shift = clipper::Coord_frac(x_shift, y_shift, z_shift);
                  clipper::RTop_orth orthop =
                     clipper::RTop_frac(spacegroup.symop(isym).rot(),
                                        spacegroup.symop(isym).trn() + cell_shift).rtop_orth(cell);
                  clipper::Coord_orth t_point = pos.transform(orthop);
                  double dd = (t_point - screen_centre).lengthsq();
                  if (dd < d_best * d_best) {
                     d_best = std::sqrt(dd);
                     rtop_best = orthop;
                     improved = true;
                  }
               }
            }
         }
      }
      if (improved) {
         clipper::Coord_orth t_point = pos.transform(rtop_best);
         v[i].first = t_point;
      }
   }
   return v;

}
