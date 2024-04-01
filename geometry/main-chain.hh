/*
 * geometry/main-chain.hh
 *
 * Copyright 2016 by Medical Research Council
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

#ifndef MAIN_CHAIN_HH
#define MAIN_CHAIN_HH

#include <string>
#include <mmdb2/mmdb_manager.h>

namespace coot {
   
   // return 0 or 1
   bool is_main_chain_p(mmdb::Atom *at);

   // return 0 or 1
   bool is_hydrogen_p(mmdb::Atom *at);

   // return 0 or 1
   bool is_main_chain_or_cb_p(mmdb::Atom *at);

   // return 0 or 1
   bool is_main_chain_p(const std::string &atom_name);

   // return 0 or 1
   bool is_main_chain_or_cb_p(const std::string &atom_name);

}


#endif // MAIN_CHAIN_HH

