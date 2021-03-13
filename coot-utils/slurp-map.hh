/* coot-utils/slurp-map.hh
 *
 * Copyright 2019 by Medical Research Council
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


#include <sys/stat.h>

#include <string>

#include <clipper/ccp4/ccp4_map_io.h>
#include "utils/ctpl.h"

namespace coot {
   namespace util {

      bool is_basic_em_map_file(const std::string &file_name);

      // inf check_only is true, then just read the header, check that it is sane
      // and return that status (don't touch the xmap). Otherwise, fill the xmap.
      //
      bool slurp_fill_xmap_from_map_file(const std::string &file_name,
                                         clipper::Xmap<float> *xmap_p,
                                         bool check_only=false);

      bool slurp_parse_xmap_data(char *data, clipper::Xmap<float> *xmap_p,
                                 bool check_only=false);

   }

}
