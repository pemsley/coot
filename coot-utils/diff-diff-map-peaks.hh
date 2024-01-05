
#ifndef COOT_DIFF_DIFF_MAP_PEAKS_HH
#define COOT_DIFF_DIFF_MAP_PEAKS_HH

#include <clipper/core/xmap.h>
#include <mmdb2/mmdb_manager.h>

namespace coot {
   std::vector<std::pair<clipper::Coord_orth, float> > diff_diff_map_peaks(const clipper::Xmap<float> &m1,
                                                                           const clipper::Xmap<float> &m2,
                                                                           float base_level);

   std::vector<std::pair<clipper::Coord_orth, float> >
   move_peaks_to_around_position(const clipper::Coord_orth &screen_centre,
                                const clipper::Spacegroup &sg,
                                const clipper::Cell &cell,
                                const std::vector<std::pair<clipper::Coord_orth, float> > &v_in);

   // perhaps I should have a wrapper function?
}

#endif
