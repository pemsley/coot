/*
 * src/map_triangle.hh
 *
 * Copyright 2022 by Medical Research Council
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
#ifndef MAP_TRIANGLE_HH
#define MAP_TRIANGLE_HH

#include "glm/glm.hpp"
#include "coot-utils/g_triangle.hh"

class map_triangle_t : public g_triangle {
public:
   glm::vec3 mid_point;
   float back_front_projection_distance;
   map_triangle_t(const unsigned int &a0,
                  const unsigned int &a1,
                  const unsigned int &a2,
                  const glm::vec3 &mp) : g_triangle(a0, a1, a2), mid_point(mp), back_front_projection_distance(0.0f) {}
   map_triangle_t(const g_triangle &gt, const glm::vec3 &mp) : g_triangle(gt), mid_point(mp), back_front_projection_distance(0.0f) {}
};


#endif // MAP_TRIANGLE_HH
