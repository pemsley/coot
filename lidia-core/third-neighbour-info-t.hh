/*
 * lidia-core/third-neighbour-info-t.hh
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


#ifndef THIRD_NEIGHBOUR_INFO_T 
#define THIRD_NEIGHBOUR_INFO_T

#include "use-rdkit.hh"

namespace cod { 
   class third_neighbour_info_t {
   public:
      const RDKit::Atom *atom_p;
      std::string ele;
      unsigned int degree;
      third_neighbour_info_t() {
	 degree = 0;
      }
      third_neighbour_info_t(const RDKit::Atom *atom_p_in,
			     const std::string &e,
			     unsigned int d) {
	 atom_p = atom_p_in;
	 ele = e;
	 degree = d;
      }
      bool operator==(const third_neighbour_info_t &t) const {
	 return (t.atom_p == atom_p);
      }
      bool operator<(const third_neighbour_info_t &t) const {
	 return (t.atom_p < atom_p);
      }
   };

}

#endif // THIRD_NEIGHBOUR_INFO_T

