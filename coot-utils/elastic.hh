/*
 * coot-utils/elastic.hh
 *
 * Copyright 2012 by University of York
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


#ifndef ELASTIC_HH
#define ELASTIC_HH

#include <mmdb2/mmdb_manager.h>
#include <string>
#include <vector>

namespace coot { 
   class elastic_network_item_t {
   public:
      mmdb::Atom *at_1;
      mmdb::Atom *at_2;
      double spring_constant;
      elastic_network_item_t(mmdb::Atom *at1, mmdb::Atom *at2, double c) {
	 at_1 = at1;
	 at_2 = at2;
	 spring_constant = c;
      }
      elastic_network_item_t() {
	 at_1 = NULL;
	 at_2 = NULL;
      }
   };
   
   class elastic_network_model_t {
      std::vector<elastic_network_item_t> d;
   public:
      elastic_network_model_t() {}
      elastic_network_model_t(mmdb::Manager *mol, int atom_selection_handle,
			      mmdb::realtype min_dist,
			      mmdb::realtype max_dist,
			      unsigned int max_n_distances);
   };

   void test_elastic();
}


#endif // ELASTIC_HH
