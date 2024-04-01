/*
 * coot-utils/arc-info.hh
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

#ifndef ARC_INFO_TYPE_HH
#define ARC_INFO_TYPE_HH

#include <clipper/core/coords.h>
#include <mmdb2/mmdb_manager.h>

namespace coot {
// can throw an exception (e.g. null pointers, overlapping atoms)
   //
   class arc_info_type {
   public:
      float delta; // the difference between the starting angle delta
                   // and the end - using the orientation_matrix, means
                   // that the start angle is 0.
      clipper::Coord_orth start_point;
      clipper::Coord_orth start_dir;
      clipper::Coord_orth normal;
      clipper::Mat33<double> orientation_matrix;
      arc_info_type(mmdb::Atom *at_1, mmdb::Atom *at_2, mmdb::Atom *at_3);
   };

}

#endif // ARC_INFO_TYPE_HH
