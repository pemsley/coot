/*
 * api/api-cell.hh
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
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#ifndef API_CELL_HH
#define API_CELL_HH

namespace api {

   //! POD cell type for moorhen.
   //!
   //! the angles are in radians.
   class cell_t {
   public:
      //! a
      float a;
      //! b
      float b;
      //! c
      float c;
      //! alpha
      float alpha;
      //! beta
      float beta;
      //! gamma
      float gamma;
      //! is_set
      bool is_set;
      cell_t() { is_set = false; }
      cell_t(float a, float b, float c, float alpha, float beta, float gamma) :
         a(a), b(b), c(c), alpha(alpha), beta(beta), gamma(gamma), is_set(true) {}
   };

}

#endif // API_CELL_HH
