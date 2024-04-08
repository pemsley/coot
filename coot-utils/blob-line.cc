/*
 * coot-utils/blob-line.cc
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

#include "coot-map-utils.hh"
#include "blob-line.hh"

std::pair<bool, clipper::Coord_orth>
coot::find_peak_along_line_favour_front(const clipper::Coord_orth &p1,
                                        const clipper::Coord_orth &p2,
                                        float contour_level,
                                        const clipper::Xmap<float> &xmap) {

   float best_score = -9999999.9;
   clipper::Coord_orth pbest;
   int istep_max = 500;
   bool point_set = false;

   for (int istep=0; istep<=istep_max; istep++) {
      float fr = float(istep)/float(istep_max);
      clipper::Coord_orth pc = p1 + fr*(p2-p1);
      float d = util::density_at_point(xmap, pc);
      if (d > contour_level) {
         // OK, so the point we want is somewhere in this peak
         for (int jstep=istep; jstep<=istep_max; jstep++) {
            fr = float(jstep)/float(istep_max);
            pc = p1 + fr*(p2-p1);
            d = util::density_at_point(xmap, pc);
            if (d > contour_level) {
               if (d> best_score) {
                  best_score = d;
                  pbest = pc;
                  point_set = true;
               }
            } else {
               // the front peak is over (now below the contour level), we have the pbest.
               break;
            }
         }
         break;
      }
   }
   return std::make_pair(point_set, pbest);
}
