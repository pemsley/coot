/*
 * coot-utils/segmap.hh
 *
 * Copyright 2019 by Medical Research Council
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

#include "clipper/core/coords.h"
#include "clipper/core/xmap.h"
#include "clipper/ccp4/ccp4_map_io.h"

// Another go at segmenting maps with symmetry considerations

namespace coot {

   // I want to give the symmetry as a n-fold a vector and centre
   //
   // or: as n-fold and presume that the axis is the z axis and goes through the middle
   // 
   // or: as n-fold and presume that the axis is the z axis and goes through the origin

   class segmap {
      std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float > > find_peaks(float c) const;
      clipper::Xmap<float> flood_from_peaks(const std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float > > &peaks,
					    float cut_off_for_flooding);
   public:
      segmap(const clipper::Xmap<float> &xmap_in) : xmap(xmap_in) {}
      const clipper::Xmap<float> &xmap;
      void proc(bool do_write_flag, const std::string &file_name);
      // remove "dust" that is smaller than 5.0A across (by default)
      void dedust(float vol=5.0);

   };

}
