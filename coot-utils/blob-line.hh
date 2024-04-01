/*
 * coot-utils/blob-line.hh
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
#ifndef BLOB_LINE_HH
#define BLOB_LINE_HH

#include <utility>
#include <clipper/core/xmap.h>

namespace coot {

std::pair<bool, clipper::Coord_orth>
find_peak_along_line_favour_front(const clipper::Coord_orth &p1,
                                  const clipper::Coord_orth &p2,
                                  float contour_level,
                                  const clipper::Xmap<float> &xmap);

}


#endif // BLOB_LINE_HH
