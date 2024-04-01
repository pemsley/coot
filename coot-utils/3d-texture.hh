/*
 * coot-utils/3d-texture.hh
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

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
