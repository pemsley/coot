/*
 * src/atom-label-info.hh
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

#ifndef ATOM_LABEL_INFO_HH
#define ATOM_LABEL_INFO_HH

#include <string>
#include <glm/glm.hpp>

// not really an atom label.

class atom_label_info_t {
public:
   std::string label;
   glm::vec3 position;
   glm::vec4 colour;
   atom_label_info_t(const std::string &l, const glm::vec3 &p, const glm::vec4 c) : label(l), position(p), colour(c) {}
};


#endif // ATOM_LABEL_INFO_HH

