/* coords/cos-sin.cc
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

/* given the sine of an angle, return the cosine, 
   e.g sqrt(3)/2 returns 1/2.

   We are only interested in values of angles between 0-90 degrees (pi/2)
   so that (for now) we don't consider negative sine/cosines.

   We create a table before this function is used which contains a
   number of reference points, which will be used for simple
   interpolation here.

*/

#include <stdlib.h>
#include <iostream>

#include <math.h>

#include "cos-sin.h"

// initialize the statics in class cos-sin.  They don't matter because
// the object is reinitiallised in init() [globjects.cc].

int    cos_sin::cos_to_sine_table_steps = 1000;
int    cos_sin::is_table_filled         = 0; 
float *cos_sin::cos_to_sine_table       = NULL;


float 
cos_sin::operator()(float v) const { 
   
   int whole_part;

   float frac_part; 
   float tmp; 
   float a1, a2, a_interp; 

   /* Recall that the sin is symmetric about pi/2 and that
      we only have cos values between 0 and pi. */
   if (v < 0.0) {
      v = -v;

   } 

   if (v > 1.0) { 
      std::cout << "Impossible cosine: " << v << std::endl;
      exit(1);
   }

   /* Check that there are values in the table */
   if (is_table_filled == FALSE) { 
     std::cout << "Need to call construct_cos_to_sin_table() first"
	   << std::endl;
      exit(1);
   } 


   tmp = (v*cos_to_sine_table_steps);

   whole_part = (int) tmp;

   frac_part = tmp - (float) whole_part; 

   if (frac_part != 0.0) { 
      /* The usual case - interpolation needed */
      a1 = cos_to_sine_table[whole_part]; 
      a2 = cos_to_sine_table[whole_part+1];
      
      a_interp = a1 + frac_part*(a2-a1); 
      return a_interp; 

   } else { 

      /* simply lookup and return the whole_part */
      return cos_to_sine_table[whole_part];
   }
 
}

/* cos_steps typically 10000 */
void 
cos_sin::fillTable (int cos_steps) { 

   int i; 
   float val, table_val;
   
   /* assign the class's holder */
   cos_to_sine_table_steps = cos_steps; 

   cos_to_sine_table = new float[cos_steps + 1];

   for(i=0; i<= cos_steps; i++) { 
      val = (float) (i/(float)cos_steps); 
      table_val = (float) sin(acos(val)); 
      cos_to_sine_table[i] = table_val;
   }

   is_table_filled = TRUE;

}
      
cos_sin::cos_sin(int isteps) {

   fillTable(isteps);
   
}

cos_sin::cos_sin() {

   // do nothing because this has already been inited in init (globjects). 
   //
   
}
   

cos_sin::~cos_sin() {

   //std::cout << "deleting cos_sin table (disabled)" << std::endl;

      // delete [] cos_to_sine_table;
} 

   
