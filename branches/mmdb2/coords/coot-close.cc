// -*-c++-*-
/* coords/coot-close.cc
 * 
 * Copyright 2006 by The University of York
 * Author: Paul Emsley
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


#include <vector>
#include "mmdb-extras.h"
#include "coot-close.hh"

// clipper::Coord_orth
// closest_approach(const clipper::Coord_orth &moving_point, 
// 		 const clipper::Coord_orth &reference_point,
// 		 const atom_selection_container_t &asc) { 

//    return closest_approach(moving_point, reference_point, asc.mol);
// }


clipper::Coord_orth
closest_approach(const clipper::Coord_orth &moving_point, 
		 const clipper::Coord_orth &reference_point,
		 CMMDBManager *mol) { 

   clipper::Coord_orth pos(-1.0, -1.0, -1.0);
   clipper::Coord_orth trans_pos;
   double closest_dist_sq = 99999999999.9;
   double d_sq;
   mat44 my_matt;
   int err;
   int nsymm = mol->GetNumberOfSymOps();
   // std::cout << "nsymm: " << nsymm << std::endl;

   for (int x_shift= -2; x_shift<= 2; x_shift++) { 
      for (int y_shift= -2; y_shift<= 2; y_shift++) { 
	 for (int z_shift= -2; z_shift<= 2; z_shift++) {
	    for (int ii=0; ii<nsymm; ii++) { 
	       err = mol->GetTMatrix(my_matt, ii, x_shift, y_shift, z_shift);
	       
	       if (err != 0) {
		  std::cout << "ERROR:: something BAD with closest_approach's GetTMatrix()\n";
	       } else {
		  
		  clipper::Mat33<double> clipper_mat(my_matt[0][0], my_matt[0][1], my_matt[0][2],
						     my_matt[1][0], my_matt[1][1], my_matt[1][2],
						     my_matt[2][0], my_matt[2][1], my_matt[2][2]);
		  clipper::Coord_orth  cco(my_matt[0][3], my_matt[1][3], my_matt[2][3]);
		  clipper::RTop_orth rtop(clipper_mat, cco);
		  
		  trans_pos = moving_point.transform(rtop);
		  
		  d_sq = (trans_pos - reference_point).lengthsq();
		  if (d_sq < closest_dist_sq) { 
		     closest_dist_sq = d_sq;
		     pos = trans_pos;
		  }
	       } 
	    }
	 }
      }
   }
   return pos;
} 


clipper::RTop_orth
closest_approach_transformation(const clipper::Coord_orth &moving_point, 
				const clipper::Coord_orth &reference_point,
				CMMDBManager *mol) {

   clipper::RTop_orth r;
   clipper::Coord_orth trans_pos;
   double closest_dist_sq = 99999999999.9;
   double d_sq;
   mat44 my_matt;
   int err;
   int nsymm = mol->GetNumberOfSymOps();
   // std::cout << "nsymm: " << nsymm << std::endl;

   for (int x_shift= -2; x_shift<= 2; x_shift++) { 
      for (int y_shift= -2; y_shift<= 2; y_shift++) { 
	 for (int z_shift= -2; z_shift<= 2; z_shift++) {
	    for (int ii=0; ii<nsymm; ii++) { 
	       err = mol->GetTMatrix(my_matt, ii, x_shift, y_shift, z_shift);
	       
	       if (err != 0) {
		  std::cout << "ERROR:: something BAD with closest_approach's GetTMatrix()\n";
	       } else {
		  
		  clipper::Mat33<double> clipper_mat(my_matt[0][0], my_matt[0][1], my_matt[0][2],
						     my_matt[1][0], my_matt[1][1], my_matt[1][2],
						     my_matt[2][0], my_matt[2][1], my_matt[2][2]);
		  clipper::Coord_orth  cco(my_matt[0][3], my_matt[1][3], my_matt[2][3]);
		  clipper::RTop_orth rtop(clipper_mat, cco);
		  
		  trans_pos = moving_point.transform(rtop);
		  
		  d_sq = (trans_pos - reference_point).lengthsq();
		  if (d_sq < closest_dist_sq) { 
		     closest_dist_sq = d_sq;
		     r = rtop;
		  }
	       } 
	    }
	 }
      }
   }
   return r;
}
