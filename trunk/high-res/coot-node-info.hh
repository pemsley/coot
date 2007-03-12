
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
