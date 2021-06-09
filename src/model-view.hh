/* src/model-view.hh
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2015 by Medical Research Council
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
#ifndef MODEL_VIEW_HH
#define MODEL_VIEW_HH

#include <string>
#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot { 
   class model_view_residue_button_info_t {
   public:
      model_view_residue_button_info_t(){} // for new allocator

      model_view_residue_button_info_t(const std::string &lab,
				       mmdb::Residue *res) : button_label(lab) {
	 residue_spec = residue_spec_t(res);
      } 
      std::string button_label;
      // mmdb::Residue *residue; No. This can go out of date.
      residue_spec_t residue_spec;
   };

   class model_view_atom_tree_item_info_t {
   public:
      // model_view_atom_tree_item_info_t() {}  // not needed?
      model_view_atom_tree_item_info_t(const std::string &label,
				       mmdb::Residue *res) : button_label(label) {
	 residue_spec = residue_spec_t(res);
      }
      std::string button_label;
      // mmdb::Residue *residue;
      residue_spec_t residue_spec;
   };

   class model_view_atom_tree_chain_t {
   public:
      model_view_atom_tree_chain_t() {} // for new allocator
      model_view_atom_tree_chain_t(const std::string &chain_id_in) : chain_id(chain_id_in) {
      }
      void add_residue(const model_view_atom_tree_item_info_t &res) {
	 tree_residue.push_back(res);
      } 
      std::vector<model_view_atom_tree_item_info_t> tree_residue;
      std::string chain_id;
   };

   // old 
   class model_view_atom_button_info_t {
   public:
      model_view_atom_button_info_t() {} // for new allocator
      model_view_atom_button_info_t(const std::string &label,
				    mmdb::Atom *atom_in) : button_label(label) {
	 atom = atom_in;
      }
      std::string button_label;
      mmdb::Atom *atom;
   };

}


#endif // MODEL_VIEW_HH

