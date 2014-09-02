/* lidia-core/lbg-shared.hh
 * 
 * Copyright 2011 by The University of Oxford
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

#ifndef LBG_SHARED_HH
#define LBG_SHARED_HH

namespace coot {


   class bash_distance_t {
   public:
      double dist;
      bool limited;
      bash_distance_t() {
	 limited = 0;
	 dist = -1;
      }
      bash_distance_t(double d) {
	 limited = 1;
	 dist = d;
      }
      bool unlimited() const {
	 return !limited;
      } 
      friend std::ostream& operator<< (std::ostream& s, const bash_distance_t &bd);
   };
   std::ostream& operator<< (std::ostream& s, const bash_distance_t &bd);
}

#endif // LBG_SHARED_HH
