/*
 * coot-utils/prideout-octasphere.hh
 *
 * Copyright 2023 by Medical Research Council
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

#ifndef PRIDEOUT_OCTASPHERE_HH
#define PRIDEOUT_OCTASPHERE_HH

#include <vector>

#include <glm/glm.hpp>

#include "coot-utils/g_triangle.hh"

std::pair<std::vector<glm::vec3>, std::vector<g_triangle> >
tessellate_octasphere_patch(unsigned int num_subdivisions);



#endif // PRIDEOUT_OCTASPHERE_HH
