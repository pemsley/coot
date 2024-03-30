/*
 * utils/win-compat.hh
 *
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
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


#include <string>

namespace coot {

   std::string get_fixed_font();
   bool is_regular_file(const std::string &file_name);
   bool is_dir_or_link( const std::string &file_name);
   std::string uri_to_file_name(const std::string &file_name);
   int rename_win(const char *old_filename, const char *new_filename);

}
