

#ifndef TRANSFER_OCCLUSIONS_HH
#define TRANSFER_OCCLUSIONS_HH

#include "density-contour-triangles.hh"

namespace coot {

   // I guess that this should live somewhere else? occlusion.cc perhaps.
   // But density_contour_triangles_container_t is not used in occlusion.hh

   void
   transfer_occlusions(const std::vector<coot::augmented_position> &positions,
                       coot::density_contour_triangles_container_t *tri_con_p);

}

#endif // TRANSFER_OCCLUSIONS_HH
