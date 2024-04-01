/*
 * geometry/srs-interface.hh
 *
 * Copyright 2013 by Medical Research Council
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

#ifndef SRS_INTERFACE_HH
#define SRS_INTERFACE_HH

#include <string>

#define MONOMER_DIR_STR "COOT_CCP4SRS_DIR"

namespace coot { 
   int get_min_match(const int &n1, const float similarity);
   std::string get_srs_dir(); // use environment variables or ccp4/prefix-dir fall-back
}

#endif // SRS_INTERFACE_HH
