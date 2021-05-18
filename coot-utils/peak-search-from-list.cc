/* coot-utils/peak-search.cc
 * 
 * Copyright 2009 by Kevin Keating
 * Copyright 2009 by the University of Oxford
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "peak-search.hh"
 
std::vector<clipper::Coord_orth>
coot::peak_search::get_peaks_from_list(const clipper::Xmap<float> &xmap,
                                       const clipper::Coord_orth &center,
                                       float radius,
                                       const std::vector<clipper::Coord_orth> &peaks) {

   std::vector<clipper::Coord_orth> r;
     
   //move the search center near the origin
   std::vector<clipper::Coord_orth> sampled_protein_coords;
   sampled_protein_coords.push_back(center);
   std::vector<int> iprotein_trans = find_protein_to_origin_translations(sampled_protein_coords, xmap);
   clipper::Mat33<double> m_identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
   clipper::Coord_frac cell_shift(iprotein_trans[0], iprotein_trans[1], iprotein_trans[2]);
   clipper::RTop_frac rtf(m_identity, cell_shift);
   clipper::Coord_frac cell_shift_rev(-iprotein_trans[0], -iprotein_trans[1], -iprotein_trans[2]);
   clipper::RTop_frac rtf_rev(m_identity, cell_shift_rev);

   clipper::RTop_orth orthop = rtf.rtop_orth(xmap.cell());
   clipper::Coord_orth translated_center = center.transform(orthop);
   
     
   //prepare the transformation operation to move points back to the original center location
   clipper::RTop_orth orthop_rev = rtf_rev.rtop_orth(xmap.cell());
    
   //go through each peak and find any translations or rotations of it that are near the translated center
   for (unsigned int i=0; i<peaks.size(); i++) {
      clipper::Coord_orth pt = peaks[i];
      std::pair<bool, clipper::Coord_orth> pair = sym_shift_test(pt, xmap.spacegroup(), xmap.cell(),
                                                                 center, radius, iprotein_trans,
                                                                 translated_center, orthop_rev);
      if (pair.first)
         r.push_back(pair.second);
   }
   return r;
}


std::pair<bool, clipper::Coord_orth>
coot::peak_search::sym_shift_test(const clipper::Coord_orth &pt,
                                  clipper::Spacegroup spacegroup,
                                  clipper::Cell cell,
                                  const clipper::Coord_orth &center,
                                  float radius,
                                  const std::vector<int> &iprotein_trans,
                                  const clipper::Coord_orth &translated_center,
                                  const clipper::RTop_orth &orthop_rev) const {

   clipper::Coord_orth r;
   bool done = 0;
   
   int shift_range = 2;
   int nsyms = spacegroup.num_symops();

   clipper::Coord_frac cell_shift_pt;
   for (int isym=0; isym<nsyms; isym++) {
      for (int x_shift = -shift_range; x_shift<=shift_range; x_shift++) { 
         for (int y_shift = -shift_range; y_shift<=shift_range; y_shift++) { 
            for (int z_shift = -shift_range; z_shift<=shift_range; z_shift++) {
               cell_shift_pt = clipper::Coord_frac(x_shift, y_shift, z_shift);
               clipper::RTop_orth orthop =
                  clipper::RTop_frac(spacegroup.symop(isym).rot(),
                                     spacegroup.symop(isym).trn() +
                                     cell_shift_pt).rtop_orth(cell);
               clipper::Coord_orth t_point = pt.transform(orthop);
               double t_dist = clipper::Coord_orth::length(t_point, translated_center);
               if (t_dist < radius) {
                  //if this translation and rotation of the point is within the radius, then add it to r
                  r = t_point.transform(orthop_rev);
                  done = 1;
               }
            }
         }
      }
   }
   return std::pair<bool, clipper::Coord_orth> (done, r);
}

   
 
