
#ifndef t3D_TEXTURE_HH
#define t3D_TEXTURE_HH

#include "density-contour/CIsoSurface.h"

class three_d_texture_t {

   const unsigned int n_width_raw = 100;
   unsigned int point_count[101][101][101];
   unsigned int n_points; // that were added to the grid
   float n_points_f; // above as float
   void init_point_count();
   float min_x;
   float min_y;
   float min_z;
   float inv_range;

public:
   three_d_texture_t() { n_points = 0; init_point_count(); }
   three_d_texture_t(const std::vector<coot::density_contour_triangles_container_t> &draw_vectors,
                     const clipper::Coord_orth &centre, float box_radius);

   float get_density(const clipper::Coord_orth &point) const; // maybe a glm::vec3?
   void fill_occlusions(coot::density_contour_triangles_container_t &contours);
};



#endif // t3D_TEXTURE_HH
