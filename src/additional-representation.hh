/* src/additional-representation.hh
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
#ifndef ADDITIONAL_REPRESENTATION_HH
#define ADDITIONAL_REPRESENTATION_HH

namespace coot {

   // representation_types
   enum { SIMPLE_LINES, STICKS, BALL_AND_STICK, LIQUORICE, SURFACE };


   class additional_representations_t { 
   public:
      bool show_it;
      int bonds_box_type;
      int representation_type;
      float bond_width;
      float sphere_radius;
      bool draw_atom_spheres_flag;
      bool draw_hydrogens_flag;
      graphical_bonds_container bonds_box;
      atom_selection_info_t atom_sel_info;
      mmdb::Manager *mol;
      int display_list_handle;
      void update_self() {
	 if (representation_type != BALL_AND_STICK || representation_type != LIQUORICE) {
	    fill_bonds_box();
	 }
      }
      void update_self_display_list_entity(int handle_in) {
	 display_list_handle = handle_in;
      }
      void fill_bonds_box();
      void core (mmdb::Manager *mol_in,
		 int representation_type_in,
		 int bonds_box_type_in,
		 float bond_width_in,
		 bool draw_hydrogens_flag_in,
		 const atom_selection_info_t &atom_sel_info_in) {
	 show_it = true;
	 mol = mol_in;
	 bond_width = bond_width_in;
	 representation_type = representation_type_in;
	 bonds_box_type = bonds_box_type_in;
	 draw_hydrogens_flag = draw_hydrogens_flag_in;
	 // draw_atom_spheres_flag = draw_atom_spheres_flag_in;
	 atom_sel_info = atom_sel_info_in;
	 fill_bonds_box();
      }
      additional_representations_t(mmdb::Manager *mol_in,
				   int representation_type_in,
				   int bonds_box_type_in,
				   float bond_width_in,
				   float sphere_radius_in,
				   bool draw_spheres_flag_in,
				   bool draw_hydrogens_flag_in,
				   const atom_selection_info_t &atom_sel_info_in) {
	 core(mol_in, representation_type_in, bonds_box_type_in, bond_width_in,
	      draw_hydrogens_flag_in, atom_sel_info_in);
	 sphere_radius = sphere_radius_in;
	 draw_atom_spheres_flag = draw_spheres_flag_in;
      }
      
      // on changind the outside (molecule_class_info_t's mol) we need
     // to change that of the additional_representations too.
     void change_mol(mmdb::Manager *mol_in) { 
       mol = mol_in;
     } 
     void clear() { 
       show_it = 0;
     } 
     std::string info_string() const;
     void add_display_list_handle(int handle) { 
       display_list_handle = handle;
     } 
   };


}

#endif // ADDITIONAL_REPRESENTATION_HH

