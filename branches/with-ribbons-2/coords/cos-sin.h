/* coords/cos-sin.h
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


#ifndef COS_SIN_H
#define COS_SIN_H

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

class cos_sin { 

   static int cos_to_sine_table_steps;
   static int is_table_filled; 

   //
   static float *cos_to_sine_table; 

 public:
   
   void fillTable(int nsteps); 
   float operator()(float) const; 

   cos_sin();
   cos_sin(int); 
   ~cos_sin();

   void check_table(void) const; 

};

#endif // COS_SIN_H   
 

