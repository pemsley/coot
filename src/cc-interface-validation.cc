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
#include "geometry/residue-and-atom-specs.hh"
using json = nlohmann::json;

void show_interesting_positions_dialog(const std::vector<coot::atom_spec_t> &atom_specs,
                                       const std::vector<coot::residue_spec_t> &residue_specs,
                                       const std::vector<coot::Cartesian> &positions) {

}

void read_interesting_places_json_file(const std::string &file_name) {

   if (coot::file_exists(file_name)) {

      // the button labels go inthe the user_data of the specs
      std::vector<coot::atom_spec_t> atom_specs;
      std::vector<coot::residue_spec_t> residue_specs;
      std::vector<coot::Cartesian> positions;

      auto get_atom_spec = [] (const json &j) {
         coot::atom_spec_t spec;
         if (j.size() == 5) {
            const json &chain_id_item  = j[0];
            const json &res_no_item    = j[1];
            const json &ins_code_item  = j[2];
            const json &atom_name_item = j[3];
            const json &alt_conf_item  = j[4];
            std::string chain_id  =  chain_id_item.get<std::string>();
            int res_no            =    res_no_item.get<int>();
            std::string ins_code  =  ins_code_item.get<std::string>();
            std::string atom_name = atom_name_item.get<std::string>();
            std::string alt_conf  =  alt_conf_item.get<std::string>();
            spec = coot::atom_spec_t(chain_id, res_no, ins_code, atom_name, alt_conf);
         }
         return spec;
      };

      auto get_residue_spec = [] (const json &j) {
         coot::residue_spec_t spec;
         if (j.size() == 5) {
            const json &chain_id_item  = j[0];
            const json &res_no_item    = j[1];
            const json &ins_code_item  = j[2];
            const json &atom_name_item = j[3];
            const json &alt_conf_item  = j[4];
            std::string chain_id  =  chain_id_item.get<std::string>();
            int res_no            =    res_no_item.get<int>();
            std::string ins_code  =  ins_code_item.get<std::string>();
            spec = coot::residue_spec_t(chain_id, res_no, ins_code);
         }
         return spec;
      };

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
         json::const_iterator it_1 = j_item.find(std::string("position-type"));
         if (it_1 != j_item.end()) {
            std::string type = it_1.value();
            std::string label;
            bool position_is_set = false;
            coot::Cartesian position(0,0,0);
            json::const_iterator it_2 = j_item.find(std::string("label"));
            if (it_2 != j_item.end()) {
               label = it_2.value();
            }
            if (type == "by-atom-spec") {
               json::const_iterator it_3 = j_item.find(std::string("atom-spec"));
               if (it_3 != j_item.end()) {
                  const json &j_atom_spec = *it_3;
                  coot::atom_spec_t atom_spec = get_atom_spec(j_atom_spec);
                  if (atom_spec.chain_id != "unset") {
                     atom_spec.string_user_data = label;
                     atom_specs.push_back(atom_spec);
                  }
               }
            }
            if (type == "by-residue-spec") {
               json::const_iterator it_3 = j_item.find(std::string("residue-spec"));
               if (it_3 != j_item.end()) {
                  const json &j_residue_spec = *it_3;
                  coot::residue_spec_t residue_spec = get_residue_spec(j_residue_spec);
                  if (residue_spec.chain_id != "unset") {
                     residue_spec.string_user_data = label;
                     residue_specs.push_back(residue_spec);
                  }
               }
            }
            if (type == "by-coordinates") {
               json::const_iterator it_3 = j_item.find(std::string("position"));
               if (it_3 != j_item.end()) {
                  const json &j_pos = *it_3;
                  unsigned int l = j_pos.size();
                  if (l == 3) {
                     std::cout << "Found a position" << std::endl;
                     const json &x_item = j_pos[0];
                     const json &y_item = j_pos[1];
                     const json &z_item = j_pos[2];
                     float x = x_item.get<float>();
                     float y = y_item.get<float>();
                     float z = z_item.get<float>();
                     coot::Cartesian c(x,y,z);
                     positions.push_back(c);
                  }
               }
            }
         }
      }
      if (! atom_specs.empty() || ! residue_specs.empty() || positions.empty()) {
         show_interesting_positions_dialog(atom_specs, residue_specs, positions);
      }
   } else {
      std::cout << "File does not exist " << file_name << std::endl;
   }
}
