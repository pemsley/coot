/*
 * src/extra-distance-restraint-markup.hh
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
#ifndef EXTRA_DISTANCE_RESTRAINT_MARKUP_HH
#define EXTRA_DISTANCE_RESTRAINT_MARKUP_HH

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

// this is converted from the bonds vector in extra-restraints-representation.hh

class extra_distance_restraint_markup_instancing_data_t {
public:
   float width;
   float length;
   glm::vec3 position; // position of the base (which will be at an atom)
   glm::mat3 orientation;
   glm::vec4 colour;
};

#endif // EXTRA_DISTANCE_RESTRAINT_MARKUP_HH
