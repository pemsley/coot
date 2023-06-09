#ifndef RAIL_POINTS_HH
#define RAIL_POINTS_HH

// 20230429-PE on the merge with the gtk3 branch, I suppose that this can be replaced by the api version.
// Of which this is a copy.

#include <vector>

namespace api {
   class rail_points_t {
   public:
      int model_rail_points_delta; // for the latest change, I mean
      int   map_rail_points_delta;
      float rmsd_of_difference_map;
      explicit rail_points_t(float rmsd) {
         model_rail_points_delta = 0;
         map_rail_points_delta = 0;
         rmsd_of_difference_map = rmsd;
      }
      rail_points_t(float rmsd_diff_map_current, const rail_points_t &rail_points_prev) {
         model_rail_points_delta = 0;
         rmsd_of_difference_map = rmsd_diff_map_current;
         map_rail_points_delta = rail_points_delta(rail_points_prev);
      }
      int rail_points_delta(const rail_points_t &prev) {
         float fudge = 2.4; // 20230117-PE makes 1000 rail points equal ~1% in R-factor for the tutorial data
         return int(100000.0 * fudge * (prev.rmsd_of_difference_map - rmsd_of_difference_map));
      }
      static int total(const std::vector<rail_points_t> &rail_point_history) {
         int sum = 0;
         for (const auto &item : rail_point_history) {
            sum += item.map_rail_points_delta;
         }
         return sum;
      }
   };
}

#endif // RAIL_POINTS_HH
