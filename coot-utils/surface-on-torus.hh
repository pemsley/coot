/*
 * coot-utils/surface-on-torus.hh
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
#ifndef SURFACE_ON_TORUS_HH
#define SURFACE_ON_TORUS_HH

#include <vector>
#include <glm/glm.hpp>
#include "simple-mesh.hh"

namespace coot {

   // Convert a height value to a vertex colour.
   // Replace the body of this function to implement a real colour ramp.
   // this is not a great function name and may need to change for disambiguation
   // later. h is the raw (unscaled) height value, height_scale is the same
   // scale factor passed to make_surface_on_torus().
   glm::vec4 height_to_colour(float h, float height_scale);

   // Generate a simple_mesh_t by mapping a 2D height grid onto a torus surface.
   //
   // The height_data grid is square (N x N) and maps uniformly over 360 degrees
   // on each axis. height_data[row][col] where row indexes the psi axis (around
   // the tube cross-section) and col indexes the phi axis (around the ring).
   //
   // R is the major radius (centre of tube to centre of torus).
   // r is the minor radius (tube radius). The height values displace vertices
   // radially outward from the base torus surface.
   //
   // Vertex normals are computed analytically from the displaced surface using
   // central differences with periodic wrapping.

   simple_mesh_t make_surface_on_torus(const std::vector<std::vector<float> > &height_data,
                                       float R, float r, float height_scale);

}

#endif // SURFACE_ON_TORUS_HH
