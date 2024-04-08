/*
 * src/matrix-utils.hh
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

#ifndef COOT_MATRIX_UTILS_HH
#define COOT_MATRIX_UTILS_HH

#include <clipper/core/coords.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "coot-utils/coot-coord-utils.hh"

clipper::Mat33<double> glm_to_mat33(const glm::mat4 &quat_mat);

glm::quat coot_quaternion_to_glm(const coot::util::quaternion &q);

coot::util::quaternion glm_to_coot_quaternion(const glm::quat &q);

#endif // COOT_MATRIX_UTILS_HH
