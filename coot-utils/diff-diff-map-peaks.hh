/*
 * coot-utils/diff-diff-map-peaks.hh
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

#ifndef COOT_DIFF_DIFF_MAP_PEAKS_HH
#define COOT_DIFF_DIFF_MAP_PEAKS_HH

#include <clipper/core/xmap.h>
#include <mmdb2/mmdb_manager.h>

namespace coot {
   std::vector<std::pair<clipper::Coord_orth, float> > diff_diff_map_peaks(const clipper::Xmap<float> &m1,
                                                                           const clipper::Xmap<float> &m2,
                                                                           float base_level);

   std::vector<std::pair<clipper::Coord_orth, float> >
   move_peaks_to_around_position(const clipper::Coord_orth &screen_centre,
                                const clipper::Spacegroup &sg,
                                const clipper::Cell &cell,
                                const std::vector<std::pair<clipper::Coord_orth, float> > &v_in);

   // perhaps I should have a wrapper function?
}

#endif
