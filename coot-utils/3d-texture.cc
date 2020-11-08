
#include <map>

#include "3d-texture.hh"

void
three_d_texture_t::init_point_count() {
   
   // initialize
   unsigned int n_width = 101;
   for (unsigned int i=0; i<n_width; i++) {
      for (unsigned int j=0; j<n_width; j++) {
         for (unsigned int k=0; k<n_width; k++) {
            point_count[i][j][k] = 0;
         }
      }
   }

}

three_d_texture_t::three_d_texture_t(const std::vector<coot::density_contour_triangles_container_t> &draw_vectors,
                                     const clipper::Coord_orth &centre,
                                     float box_radius) {

   init_point_count();
   min_x = centre.x() - box_radius;
   min_y = centre.y() - box_radius;
   min_z = centre.z() - box_radius;
   float range = 2.0 * box_radius;
   inv_range = 1.0/range;
   for (unsigned int i=0; i<draw_vectors.size(); i++) {
      const std::vector<clipper::Coord_orth> v = draw_vectors[i].points;
      for (auto p : v) { // oh how I chuckled when I accidently wrote this in Python syntax
         float frac_x = (p.x()-min_x) * inv_range;
         float frac_y = (p.y()-min_y) * inv_range;
         float frac_z = (p.z()-min_z) * inv_range;
         int if_x = static_cast<int>(100.0f * frac_x);
         int if_y = static_cast<int>(100.0f * frac_y);
         int if_z = static_cast<int>(100.0f * frac_z);
         std::cout << "frac_x " << frac_x << " frac_y " << frac_y << " frac_z " << frac_z << std::endl;
         std::cout << "if_x " << if_x << " if_y " << if_y << " if_z " << if_z << std::endl;
         if (frac_x >= 0.0)
            if (frac_y >= 0.0)
               if (frac_z >= 0.0)
                  if (frac_x < 1.0)
                     if (frac_y < 1.0)
                        if (frac_z < 1.0)
                           point_count[if_x][if_y][if_z] += 1;
      }
      n_points += draw_vectors.size();
   }

   n_points_f = static_cast<float>(n_points);
}

float
three_d_texture_t::get_density(const clipper::Coord_orth &p) const {

   float r = 0.0f;
   float frac_x = (static_cast<float>(p.x())-min_x) * inv_range;
   float frac_y = (static_cast<float>(p.y())-min_y) * inv_range;
   float frac_z = (static_cast<float>(p.z())-min_z) * inv_range;
   int if_x = static_cast<int>(100.0f * frac_x);
   int if_y = static_cast<int>(100.0f * frac_y);
   int if_z = static_cast<int>(100.0f * frac_z);

   std::cout << "de-indexing " << if_x << " " << if_y << " " << if_z << std::endl;
 
   if (frac_x >= 0.0)
      if (frac_y >= 0.0)
         if (frac_z >= 0.0)
            if (frac_x < 1.0)
               if (frac_y < 1.0)
                  if (frac_z < 1.0) {
                     float density = 1000.0 * static_cast<float>(point_count[if_x][if_y][if_z]);
                     r = density / n_points_f; // is this the right normalization?
                     std::cout << "density: " << density << " r " << r
                               << " n_points_f " << n_points_f << std::endl;
                  }

   return r;
}

void
three_d_texture_t::fill_occlusions(coot::density_contour_triangles_container_t &contours) {

   contours.occlusion_factor.resize(contours.points.size(), 0.0f);
   for (unsigned int i=0; i<contours.points.size(); i++) {
      const clipper::Coord_orth &pt(contours.points[i]);
      float d = get_density(pt);
      contours.occlusion_factor[i] = d;
      std::cout << "occlusion_factor " << i << " " << d << "\n";
   }
   
}
