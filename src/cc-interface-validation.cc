/*
 * src/cc-interface-validation.cc
 *
 * Copyright 2025 by Medical Research Council
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

#include "cc-interface.hh"
#include "coot-utils/json.hpp"
using json = nlohmann::json;

void read_interesting_places_json_file(const std::string &file_name) {

   if (coot::file_exists(file_name)) {
      
      std::string s;
      std::fstream f(file_name);
      f.seekg(0, std::ios::end);
      s.reserve(f.tellg());
      f.seekg(0, std::ios::beg);
      s.assign((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
      json j = json::parse(s);
      unsigned int n_outer = j.size();
      std::cout << "Found " << n_outer << " parts in the json" << std::endl;
      json j_items = j["items"];
      unsigned int n_items = j_items.size();
      std::cout << "items has " << n_items << " items" << std::endl;
      for (std::size_t i=0; i<n_items; i++) {
	 const json &j_item = j_items[i];
	 json::const_iterator it = j_item.find(std::string("type"));
	 if (it != j_item.end()) {
	    std::string type = it.value();
	    if (type == "Chiral") {
	    }
	    if (type == "Bad Density Fit") {
	    }
	    if (type == "Blob") {
	       it = j_item.find(std::string("position"));
	       if (it != j_item.end()) {
		  const json &j_pos = *it;
		  unsigned int l = j_pos.size();
		  if (l == 3) {
		     std::cout << "Found a position" << std::endl;
		     const json &x_item = j_pos[0];
		     const json &y_item = j_pos[1];
		     const json &z_item = j_pos[2];
		     float x = x_item.get<float>();
		     float y = y_item.get<float>();
		     float z = z_item.get<float>();
		  }
	       }
	    }
	 }
      }
   } else {
      std::cout << "File does not exist " << file_name << std::endl;
   }
}
