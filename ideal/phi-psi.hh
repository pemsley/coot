/*
 * ideal/phi-psi.hh
 * 
 * Copyright 2017 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#ifndef PHI_PSI_T
#define PHI_PSI_T

#include <iosfwd>

namespace coot {

   // consolidate with phi_psi_pair (i.e. make that disappear)
   class phi_psi_t {
      // and now with tau!
   public:
      phi_psi_t(float a, float b) {
	 phi = a;
	 psi = b;
	 tau = -20; // unset value
      }
      phi_psi_t(float a, float b, float c) {
	 phi = a;
	 psi = b;
	 tau = c;
      }
      float phi;
      float psi;
      float tau;
      friend std::ostream &operator<<(std::ostream &s, const phi_psi_t &pp);
   };
   std::ostream &operator<<(std::ostream &s, const phi_psi_t &pp);

}

#endif // PHI_PSI_T
