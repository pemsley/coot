/*
 * density-contour/transfer-occlusions.hh
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */



#ifndef TRANSFER_OCCLUSIONS_HH
#define TRANSFER_OCCLUSIONS_HH

#include "density-contour-triangles.hh"

namespace coot {

   // I guess that this should live somewhere else? occlusion.cc perhaps.
   // But density_contour_triangles_container_t is not used in occlusion.hh

   void
   transfer_occlusions(const std::vector<coot::augmented_position> &positions,
                       coot::density_contour_triangles_container_t *tri_con_p);

}

#endif // TRANSFER_OCCLUSIONS_HH
