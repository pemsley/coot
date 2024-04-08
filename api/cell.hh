/*
 * api/cell.hh
 * 
 * Copyright 2020 by Medical Research Council
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef COOT_API_CELL_HH
#define COOT_API_CELL_HH

namespace coot {
   //! Simple Cell class for transfering symmetry info.
   class Cell {
   public:
      //! the unit cell lengths and angles in radians.
      float a, b, c, alpha, beta, gamma;
      Cell(float a, float b, float c, float alpha, float beta, float gamma) :
         a(a), b(b), c(c), alpha(alpha), beta(beta), gamma(gamma) {}
      Cell() { a = -1; b = -1; c = -1; alpha = -1; beta = -1; gamma = -1; }
   };

}


#endif // COOT_API_CELL_HH
