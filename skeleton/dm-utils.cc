/*
 * skeleton/dm-utils.cc
 *
 * Copyright 2007 by University of York
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

#include <iostream>

// fortran indexing.  The (key and ptr) arrays need to be of size
// n+1, since we use key[n].
//
// The sorted postions are 1..n.
// 
void 
shsorti(float *key, int *ptr, int n) {

   int l=1;
   int i;
   int i1, i2;
   int p1, p2; 

   while (l<n) l *=2;
   std::cout << "l set to " << l << std::endl;
      
   while (l > 1) {
      l /= 2;
      for (i=1; i<=n-l; i++) {
	 for (i1=i; i1>=1; i1--) {
	    i2 = i1 + 1;
	    p1 = ptr[i1];
	    p2 = ptr[i2];

	    if (key[p2] >= key[p1]) break;
	    
	    std::cout << "assigning ptr [" << i1 << "] as " << p2 << std::endl;
	    std::cout << "assigning ptr [" << i2 << "] as " << p1 << std::endl;
	    ptr[i1] = p2;
	    ptr[i2] = p1;
	 }
      }
   }
}



