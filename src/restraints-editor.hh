
/* src/restraints-editor.hh
 * 
 * Copyright 2008 by The University of Oxford
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

#ifndef RESTRAINTS_EDITOR_HH
#define RESTRAINTS_EDITOR_HH

#include <gtk/gtk.h>
#include "geometry/protein-geometry.hh"

namespace coot { 
   class restraints_editor {
      class view_and_store {
      public:
	 GtkTreeView *view;
	 GtkTreeStore *store;
	 view_and_store(GtkTreeView *view_in,
			GtkTreeStore *store_in) {
	    view = view_in;
	    store = store_in;
	 }
	 view_and_store() {
	    view = NULL;
	    store = NULL;
	 }
      };
      GtkBuilder *builder; // set in setup_builder();
      void setup_builder();
      GtkWidget *widget_from_builder(const std::string &widget_name);
   public:
      enum {TREE_TYPE_INFO, TREE_TYPE_CHIRALS, TREE_TYPE_TORSIONS, TREE_TYPE_ANGLES,
	    TREE_TYPE_BONDS, TREE_TYPE_PLANES, TREE_TYPE_ATOMS};
      enum { UNKNOWN_TYPE= -456723 };
      restraints_editor() {
         builder = nullptr;
	 dialog = NULL;
	 is_valid_flag = false;
	 max_number_of_atoms_in_plane = -1;
         setup_builder(); // assigns builder
      }
      void fill_dialog(const dictionary_residue_restraints_t &restraints); // set is_valid_flag
      dictionary_residue_restraints_t make_restraint() const;
      bool is_valid() const { return is_valid_flag; }
      bool matches_dialog(GtkWidget *w) const {
	 // std::cout << " comparing " << dialog << " vs " << w << std::endl;
	 return (w == dialog);}
      void delete_restraint(GtkWidget *w);
      void add_restraint(GtkWidget *w);
      GtkWidget *get_dialog() { return dialog; }

   private:
      enum { TORSION_COL_PERIODICTY = 7,
      CHIRAL_COL_SIGN = 5 };
      GtkWidget *dialog;
      bool is_valid_flag;
      view_and_store view_and_store_atoms;
      view_and_store view_and_store_info;
      view_and_store view_and_store_bonds;
      view_and_store view_and_store_angles;
      view_and_store view_and_store_torsions;
      view_and_store view_and_store_chirals;
      view_and_store view_and_store_planes;

//       GtkTreeStore *tree_store_atoms;
//       GtkTreeStore *tree_store_angles;
//       GtkTreeStore *tree_store_info;
//       GtkTreeStore *tree_store_bonds;
//       GtkTreeStore *tree_store_torsions;
//       GtkTreeStore *tree_store_chirals;
//       GtkTreeStore *tree_store_planes;
      
      int max_number_of_atoms_in_plane;
      static int get_column_type(int tree_type, int column_number, int max_n_plane_atoms);
      void fill_info_tree_data(GtkWidget *restraints_editor_dialog,
			       const coot::dictionary_residue_restraints_t &restaints);
      void fill_atom_tree_data(GtkWidget *restraints_editor_dialog,
			       const coot::dictionary_residue_restraints_t &restaints);
      void fill_bond_tree_data(GtkWidget *restraints_editor_dialog,
			       const coot::dictionary_residue_restraints_t &restaints);
      void fill_angle_tree_data(GtkWidget *restraints_editor_dialog,
				const coot::dictionary_residue_restraints_t &restaints);
      void fill_torsion_tree_data(GtkWidget *restraints_editor_dialog,
				  const coot::dictionary_residue_restraints_t &restaints);
      void fill_chiral_tree_data(GtkWidget *restraints_editor_dialog,
				 const coot::dictionary_residue_restraints_t &restaints);
      void fill_plane_tree_data(GtkWidget *restraints_editor_dialog,
				const coot::dictionary_residue_restraints_t &restaints);
      void add_plane_cell_renderer(GtkTreeView *tree_view,
				   GtkTreeStore *store, 
				   const std::string &column_title, int pos,
				   int tree_type,
				   int max_n_plane_atoms);
      GtkCellRenderer *add_cell_renderer(GtkTreeView *tree_view,
					 GtkTreeStore *store, 
					 const std::string &column_title, int pos,
					 int tree_type);
      GtkTreeStore *make_tree_store_for_planes(int natoms);
      static void cell_edited_callback (GtkCellRendererText *cell,
					gchar               *path_string,
					gchar               *new_text,
					gpointer             user_data);

      std::vector<coot::dict_bond_restraint_t> 
      get_bond_restraints() const;
      std::vector<coot::dict_angle_restraint_t> 
      get_angle_restraints() const;
      std::vector<coot::dict_torsion_restraint_t> 
      get_torsion_restraints() const;
      std::vector<coot::dict_chiral_restraint_t> 
      get_chiral_restraints() const;
      std::vector<coot::dict_plane_restraint_t> 
      get_plane_restraints() const;
      // also return if a flag we have partial charges
      std::pair<bool, std::vector<coot::dict_atom> > get_atom_info() const;
      // also return a flag if the dict_chem_comp_t was valid
      std::pair<bool, coot::dict_chem_comp_t> get_residue_info() const;
      GtkTreeStore *get_tree_store_by_notebook_page(gint current_page_index) const;
      GtkTreeView *get_tree_view_by_notebook_page(gint current_page_index) const;
      std::string make_chiral_volume_string(int chiral_sign) const;
      int chiral_volume_string_to_chiral_sign(const std::string &chiral_vol_string) const;
      
   };
}

#endif // RESTRAINTS_EDITOR_HH
