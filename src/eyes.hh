/*
 * src/eyes.hh
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

#include <vector>
#include "generic-vertex.hh"
#include "coot-utils/g_triangle.hh"


// construct these around the origin - the instanced matrices will put them
// in the right place.
std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
make_spherical_surface_circular_patch(float atom_radius,
                                      float solid_theta,
                                      float h_scale,
                                      float v_scale,
                                      unsigned int n_slices);


std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
make_curved_bar_patch(float atom_radius, float theta_start, float theta_end, float height,
                      unsigned int n_slices, float curve);
