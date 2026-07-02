/* layla/url_fetch.hpp
 *
 * Copyright 2026 by Medical Research Council
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

#ifndef LAYLA_URL_FETCH_HPP
#define LAYLA_URL_FETCH_HPP

#include <string>

// Layla's own minimal libcurl interface. Layla must not depend on src/, so this
// does not use src/'s coot_get_url(); it is self-contained (Layla already links
// libcurl). Modelled on the curl usage in layla/python_utils.cpp.

namespace coot {
   namespace layla {

      // Fetch url to a string. Returns an empty string on failure.
      std::string get_url_as_string(const std::string &url);

      // Fetch url to output_file_name (written via a .tmp file then renamed).
      // Returns 0 on success, non-zero on failure (the file is not created on
      // failure). A response shorter than min_bytes is treated as failure
      // (guards against error pages saved as content).
      int get_url(const std::string &url, const std::string &output_file_name,
                  long min_bytes = 1);
   }
}

#endif // LAYLA_URL_FETCH_HPP
