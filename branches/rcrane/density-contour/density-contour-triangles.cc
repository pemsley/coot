
#include <algorithm>

#include "CIsoSurface.h"

void 
coot::density_contour_triangles_container_t::depth_sort(const clipper::Coord_orth &back_plane_point,
							const clipper::Coord_orth &front_plane_point) { 

   clipper::Coord_orth back_front = front_plane_point - back_plane_point;
   double bf_squared = back_front.lengthsq();
   if (bf_squared < 0.000001)
      bf_squared = 0.000001;
   for (unsigned int i=0; i<point_indices.size(); i++) {
      clipper::Coord_orth h_back = point_indices[i].mid_point - back_plane_point;
      double dot = clipper::Coord_orth::dot(back_front, h_back);
      point_indices[i].back_front_projection_distance = dot * dot / bf_squared;
   }
   std::sort(point_indices.begin(), point_indices.end());
	     
}

void
coot::density_contour_triangles_container_t::calculate_normals() {

   // for each triangle
   //    for each of the 3 vertices
   //       add the flat shading normal to the sum of normals of that vertex
   // then normalise the normals

   std::vector<clipper::Coord_orth> sum_normals(normals.size());
   clipper::Coord_orth zero(0,0,0);
   for (unsigned int i=0; i<sum_normals.size(); i++)
      sum_normals[i] = zero;

   for (unsigned int i=0; i<point_indices.size(); i++) { 
      for (int j=0; j<3; j++) {
	 sum_normals[point_indices[i].pointID[j]] += point_indices[i].normal_for_flat_shading;
      }
   }
   for (unsigned int i=0; i<points.size(); i++)
      normals[i] = clipper::Coord_orth(sum_normals[i].unit());

}
