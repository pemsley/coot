/*
 * geometry/mol-utils-2.hh
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
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

// # put this in mol-utils.cc when that arrives in this branch

#ifndef MOL_UTILS_2_HH
#define MOL_UTILS_2_HH

#include <utility>
#include <map>
#include <string>

#include "mmdb2/mmdb_manager.h"

namespace coot {

   std::map<std::string, std::pair<int, int> > get_residue_number_limits(mmdb::Manager *mol);

}


#endif // MOL_UTILS_2_HH
