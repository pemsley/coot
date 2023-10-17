#ifndef FIND_WATER_BADDIES_HH
#define FIND_WATER_BADDIES_HH

#include <vector>
#include "geometry/residue-and-atom-specs.hh"
#include "atom-selection-container.hh"
#include <clipper/core/xmap.h>

namespace coot {
   std::vector <atom_spec_t>
   find_water_baddies_OR(atom_selection_container_t asc,
                         float b_factor_lim, const clipper::Xmap<float> &xmap_in,
                         float map_in_sigma,
                         float outlier_sigma_level,
                         float min_dist, float max_dist,
                         short int ignore_part_occ_contact_flag,
                         short int ignore_zero_occ_flag);

}


#endif // FIND_WATER_BADDIES_HH
