/* lbg/some-coot-utils.hh
 * 
 * Author: Paul Emsley
 * Copyright 2010 by The University of Oxford
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

#include <string>
#include <vector>

namespace coot {

   bool is_directory_p(const std::string &dir);

   namespace util { 
      std::string int_to_string(int i);
      std::string float_to_string(float f);
      std::vector<std::string> split_string_no_blanks(const std::string &string_in,
						      const std::string &splitter);
      std::string downcase(const std::string &s);
   }
}
