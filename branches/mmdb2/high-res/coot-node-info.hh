/* high-res/coot-node-info.hh
 * 
 * Copyright 2004  The University of York
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

#include "clipper/core/symop.h"

namespace coot { 

   // Store the index of a neighbour of this node and how to get to it.
   // 
   class node_info { 
   public:
      short int symm_trans_needed_flag;
      int index;
      clipper::RTop_orth rtop;
      node_info(CMMDBManager *mol, int index_in, int isym, int ix, int iy, int iz) {

	 // construct the RTop_orth rtop here:
	 // 
	 index = index_in;
	 symm_trans_needed_flag = 1;

	 mat44 my_matt;
	 mol->GetTMatrix(my_matt, isym, ix, iy, iz);
	 for (int i=0; i<4; i++) 
	    for (int j=0; j<4; j++) 
	       if (fabs(my_matt[i][j]) < 0.0001)
		  my_matt[i][j] = 0.0;
	 clipper::Mat33<double> clipper_mat(my_matt[0][0], my_matt[0][1], my_matt[0][2],
					    my_matt[1][0], my_matt[1][1], my_matt[1][2],
					    my_matt[2][0], my_matt[2][1], my_matt[2][2]);
	 clipper::Coord_orth  cco(my_matt[0][3], my_matt[1][3], my_matt[2][3]);
	 rtop = clipper::RTop_orth(clipper_mat, cco);
      }
      node_info(int index_in) { 
	 symm_trans_needed_flag = 0;
	 index = index_in;
      }
   };
}
