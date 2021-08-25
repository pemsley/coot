/* src/coot-colour.cc
 * 
 * Copyright 2016 by Medical Research Council
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
#include <iostream>
#include <math.h>

#include "compat/coot-sysdep.h"

#include "coot-colour.hh"

namespace coot { 
   std::ostream& operator<<(std::ostream &s, colour_t col) {
      s << col.col[0] << " " << col.col[1] << " " << col.col[2];
      return s;
   }

   std::vector<float> colour_t::convert_to_hsv() const {

      std::vector<float> hsv(3);

      float maxc = -1.0;
      float minc = 9.0;

      for (unsigned int i=0; i<3; i++) {
	 if (maxc < col[i]) maxc = col[i];
	 if (minc > col[i]) minc = col[i];
      }
      hsv[2] = maxc;

      if (minc == maxc) {
	 hsv[0] = 0.0;
	 hsv[1] = 0.0;
	 hsv[2] = maxc;
      } else { 

	 float range = maxc - minc;
	 hsv[1] = range/maxc;
	 float rc = (maxc - col[0]) / range;
	 float gc = (maxc - col[1]) / range;
	 float bc = (maxc - col[2]) / range;
	 if (col[0] == maxc) {
	    hsv[0] = bc-gc;
	 } else {
	    if (col[1]==maxc) {
	       hsv[0] = 2.0+rc-bc;
	    } else {
	       hsv[0] = 4.0 + gc-rc;
	    }
	 }
	 hsv[0] = hsv[0]/6.0- floorf(hsv[0]/6.0);
      }
      return hsv;
   }

   void colour_t::convert_from_hsv(const std::vector<float> &hsv) {

      if (hsv[1] == 0.0) {
	 col[0] = hsv[2];
	 col[1] = hsv[2];
	 col[2] = hsv[2];
      } else {
	 float fi = floorf(hsv[0]*6.0);
	 float f  = (hsv[0]*6.0) - fi;
	 float p = hsv[2]*(1.0 - hsv[1]);
	 float q = hsv[2]*(1.0 - hsv[1]*f);
	 float t = hsv[2]*(1.0 - hsv[1]*(1.0-f));

	 int i = int(fi);
	 switch (i) {

	 case 0:
	 case 6:
	    col[0] = hsv[2]; 
	    col[1] = t; 
	    col[2] = p;
	    break;

	 case 1:
	    col[0] = q;
	    col[1] = hsv[2]; 
	    col[2] = p;
	    break;

	 case 2:
	    col[0] = p;
	    col[1] = hsv[2]; 
	    col[2] = t;
	    break;

	 case 3:
	    col[0] = p;
	    col[1] = q; 
	    col[2] = hsv[2];
	    break;

	 case 4:
	    col[0] = t;
	    col[1] = p; 
	    col[2] = hsv[2];
	    break;

	 case 5:
	    col[0] = hsv[2];
	    col[1] = p; 
	    col[2] = q;
	    break;
	 }
      }
   }

   // this should be rotate_to()
   void colour_t::rotate(float amount) {
      std::vector<float> hsv = convert_to_hsv();
      hsv[0] += amount;
      if (hsv[0] > 1.0) hsv[0] -= 1.0;
      convert_from_hsv(hsv);
   }
}
   
