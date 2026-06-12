/*
 * coot-utils/map-cap.hh
 *
 * Copyright 2026 by Medical Research Council
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

#ifndef MAP_CAP_HH
#define MAP_CAP_HH

#include <clipper/core/xmap.h>
#include <clipper/core/coords.h>
#include "simple-mesh.hh"

namespace coot {

   //! Make a map cap mesh: a 2D cross-section (marching squares) combined with the
   //! 3D isosurface trimmed behind the clipping plane, stitched together at the boundary.
   //!
   //! @param xmap the electron density map
   //! @param contour_level the contour level
   //! @param base_point the origin of the cap plane
   //! @param x_axis_uv the x-axis unit vector of the cap plane
   //! @param y_axis_uv the y-axis unit vector of the cap plane
   //! @param x_axis_step_size grid spacing along x
   //! @param y_axis_step_size grid spacing along y
   //! @param n_x_axis_points number of grid points along x
   //! @param n_y_axis_points number of grid points along y
   //! @param isosurface_vertices the pre-computed 3D isosurface vertices
   //! @param isosurface_triangles the pre-computed 3D isosurface triangles
   //! @param other_map_for_colouring_p optional map for colouring (nullptr to use xmap)
   //! @param other_map_for_colouring_min_value min value for colour ramp
   //! @param other_map_for_colouring_max_value max value for colour ramp
   //! @param radial_map_colour_saturation saturation for colour ramp
   //! @return the combined cap mesh
   simple_mesh_t make_map_cap_mesh(const clipper::Xmap<float> &xmap,
                                   float contour_level,
                                   const clipper::Coord_orth &base_point,
                                   const clipper::Coord_orth &x_axis_uv,
                                   const clipper::Coord_orth &y_axis_uv,
                                   double x_axis_step_size,
                                   double y_axis_step_size,
                                   unsigned int n_x_axis_points,
                                   unsigned int n_y_axis_points,
                                   const std::vector<api::vnc_vertex> &isosurface_vertices,
                                   const std::vector<g_triangle> &isosurface_triangles,
                                   const clipper::Xmap<float> *other_map_for_colouring_p = nullptr,
                                   float other_map_for_colouring_min_value = 0.0f,
                                   float other_map_for_colouring_max_value = 1.0f,
                                   float radial_map_colour_saturation = 0.5f);

}

#endif // MAP_CAP_HH
