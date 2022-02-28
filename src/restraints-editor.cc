/* src/restraints-editor.cc
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

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#include <vector>
#include <string>
#include <iostream>

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>  // strcpy
#include <gtk/gtk.h>

#include "interface.h"
#include "support.h"

#include "geometry/protein-geometry.hh"
#include "utils/coot-utils.hh"

#include "coot-fileselections.h"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"

#include "restraints-editor.hh"
#include "restraints-editor-c.h"

#include "graphics-info.h"

std::vector<coot::restraints_editor> graphics_info_t::restraints_editors;

void
coot::restraints_editor::fill_dialog(const coot::dictionary_residue_restraints_t &restraints) { 
   dialog = create_restraints_editor_dialog(); // defined in interface.h
//    std::cout << "restraints editor saving "
// 	     << dialog << std::endl;
   fill_info_tree_data   (dialog, restraints);
   fill_atom_tree_data   (dialog, restraints);
   fill_bond_tree_data   (dialog, restraints);
   fill_angle_tree_data  (dialog, restraints);
   fill_torsion_tree_data(dialog, restraints);
   fill_chiral_tree_data (dialog, restraints);
   fill_plane_tree_data  (dialog, restraints);
   gtk_widget_show (dialog);
   gtk_window_present(GTK_WINDOW(dialog));
   is_valid_flag = 1;
}

void
coot::restraints_editor::fill_atom_tree_data(GtkWidget *restraints_editor_dialog,
					     const coot::dictionary_residue_restraints_t &restraints) { 

   GtkWidget *atoms_treeview = lookup_widget(restraints_editor_dialog, "atoms_treeview");

   GtkTreeView *tv_atoms = GTK_TREE_VIEW(atoms_treeview);
   GtkTreeStore *tree_store_atoms = gtk_tree_store_new (4, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
							// G_TYPE_FLOAT  // for partial change
                                                        G_TYPE_INT // for formal change
                                                        );

   view_and_store_atoms.store = tree_store_atoms;
   view_and_store_atoms.view  = tv_atoms;
   
   GtkTreeIter toplevel;
   gtk_tree_view_set_model(tv_atoms, GTK_TREE_MODEL(tree_store_atoms));

   for (unsigned int i=0; i<restraints.atom_info.size(); i++) {
      gtk_tree_store_append(tree_store_atoms, &toplevel, NULL);

      // std::pair<bool, float> pcp = restraints.atom_info[i].partial_charge;
      // float pc = pcp.second;
      // if (! pcp.first) pc = -999.99;

      std::pair<bool, int> fcp = restraints.atom_info[i].formal_charge;
      int fc = 0;
      if (fcp.first) fc = fcp.second;

      gtk_tree_store_set(tree_store_atoms, &toplevel,
			 0, restraints.atom_info[i].atom_id_4c.c_str(),
			 1, restraints.atom_info[i].type_symbol.c_str(),
			 2, restraints.atom_info[i].type_energy.c_str(),
			 3, fc,
			 -1);
   }

   int tree_type = coot::restraints_editor::TREE_TYPE_ATOMS;
   add_cell_renderer(tv_atoms, tree_store_atoms, "Atom Name",      0, tree_type);
   add_cell_renderer(tv_atoms, tree_store_atoms, "Element",        1, tree_type);
   add_cell_renderer(tv_atoms, tree_store_atoms, "Energy Type",    2, tree_type);
   // add_cell_renderer(tv_atoms, tree_store_atoms, "Partial Charge", 3, tree_type);  not today.
   add_cell_renderer(tv_atoms, tree_store_atoms, "Formal Charge", 3, tree_type);
}

GtkCellRenderer *
coot::restraints_editor::add_cell_renderer(GtkTreeView *tree_view,
					   GtkTreeStore *tree_store,
					   const std::string &column_title,
					   int pos,
					   int tree_type) {

   
  GtkCellRenderer *cell_renderer = gtk_cell_renderer_text_new();
  GtkTreeViewColumn *column = gtk_tree_view_column_new();
  gtk_tree_view_column_set_title(column, column_title.c_str());
  gtk_tree_view_append_column(tree_view, column);
  gtk_tree_view_column_pack_start(column, cell_renderer, TRUE);
  gtk_tree_view_column_add_attribute(column, cell_renderer, "text", pos);
  g_object_set(cell_renderer, "editable", TRUE, NULL);

  gtk_tree_view_column_set_sort_column_id(column, pos);

  g_object_set_data (G_OBJECT (cell_renderer), "column", GINT_TO_POINTER (pos));
  g_object_set_data (G_OBJECT (cell_renderer), "tree_type", GINT_TO_POINTER (tree_type));
  g_signal_connect(cell_renderer, "edited", (GCallback) cell_edited_callback,
		   (gpointer) tree_store);
  return cell_renderer;
}

void
coot::restraints_editor::add_plane_cell_renderer(GtkTreeView *tree_view,
						 GtkTreeStore *tree_store,
						 const std::string &column_title,
						 int pos,
						 int tree_type,
						 int max_n_plane_atoms) {

   GtkCellRenderer *cell_renderer = add_cell_renderer(tree_view, tree_store, column_title,
						      pos, tree_type);
   g_object_set_data (G_OBJECT (cell_renderer), "max_n_plane_atoms",
		      GINT_TO_POINTER (max_n_plane_atoms));
}

// static
void
coot::restraints_editor::cell_edited_callback (GtkCellRendererText *cell,
					       gchar               *path_string,
					       gchar               *new_text,
					       gpointer             user_data) {

   // Path is telling us the row.  What is the column?

   GtkTreeModel *model = (GtkTreeModel *) user_data;

   GtkTreePath *path = gtk_tree_path_new_from_string(path_string);
   gint column = GPOINTER_TO_INT (g_object_get_data (G_OBJECT (cell), "column"));
   gint tree_type = GPOINTER_TO_INT (g_object_get_data (G_OBJECT (cell), "tree_type"));
   gint max_n_plane_atoms = GPOINTER_TO_INT (g_object_get_data (G_OBJECT (cell), "max_n_plane_atoms"));
   GtkTreeIter iter;
   int data_type = get_column_type(tree_type, column, max_n_plane_atoms);
   
   gtk_tree_model_get_iter (model, &iter, path);

   if (data_type == G_TYPE_FLOAT) {
   
      float f = atof (new_text);
//       std::cout << "shoving in a float " << f << " to the model " << std::endl;
//       std::cout << " to column " << column << " path_string "<< path_string << std::endl;
      gtk_tree_store_set (GTK_TREE_STORE (model), &iter,
			  column, f,
			  -1);
   } 
   if (data_type == G_TYPE_INT) {
      int ii = atoi (new_text);
//       std::cout << "shoving in an int " << ii << " to the model " << std::endl;
//       std::cout << " to column " << column << " path_string "<< path_string << std::endl;
      gtk_tree_store_set (GTK_TREE_STORE (model), &iter,
			  column, ii,
			  -1);
   } 
   if (data_type == G_TYPE_STRING) {
//       std::cout << "shoving in an string " << new_text << " to the model " << std::endl;
//       std::cout << " to column " << column << " path_string "<< path_string << std::endl;
      gtk_tree_store_set (GTK_TREE_STORE (model), &iter,
			  column, new_text,
			  -1);
   } 
   
}

// static
int
coot::restraints_editor::get_column_type(int tree_type, int column_number, int max_n_plane_atoms) {

   int r = coot::restraints_editor::UNKNOWN_TYPE;
   int max_plane_col_no = max_n_plane_atoms;
   if (tree_type == coot::restraints_editor::TREE_TYPE_INFO) {
      switch (column_number) {
      case(4):
	 r = G_TYPE_INT;
	 break;
      case(5):
	 r = G_TYPE_INT;
	 break;
      default:
	 r = G_TYPE_STRING;
      }
   }
   if (tree_type == coot::restraints_editor::TREE_TYPE_ATOMS) {
      switch(column_number) {
      case(3):
	 r = G_TYPE_INT; // now formal chage (was partial charge)
	 break;
      default:
	 r = G_TYPE_STRING;
      }
   }
   if (tree_type == coot::restraints_editor::TREE_TYPE_CHIRALS) {
      switch(column_number) {
      case(5):
	 // r = G_TYPE_INT; // 20091019 chiral entity is now a string,
	 // from Andrew Leslie comment.
	 r = G_TYPE_STRING;
	 break;
      default:
	 r = G_TYPE_STRING;
      } 
   } 
   if (tree_type == coot::restraints_editor::TREE_TYPE_TORSIONS) {
      switch(column_number) {
      case(5):
	 r = G_TYPE_FLOAT;
	 break;
      case(6):
	 r = G_TYPE_FLOAT;
	 break;
      case(7):
	 r = G_TYPE_INT;
	 break;
      default:
	 r = G_TYPE_STRING;
      } 
   } 
   if (tree_type == coot::restraints_editor::TREE_TYPE_ANGLES) {
      switch(column_number) {
      case(3):
	 r = G_TYPE_FLOAT;
	 break;
      case(4):
	 r = G_TYPE_FLOAT;
	 break;
      default:
	 r = G_TYPE_STRING;
      } 
   } 
   if (tree_type == coot::restraints_editor::TREE_TYPE_BONDS) {
      switch(column_number) {
      case(3):
	 r = G_TYPE_FLOAT;
	 break;
      case(4):
	 r = G_TYPE_FLOAT;
	 break;
      default:
	 r = G_TYPE_STRING;
      } 
   } 
   if (tree_type == coot::restraints_editor::TREE_TYPE_PLANES) {
      if (column_number > max_plane_col_no) { 
	 r = G_TYPE_FLOAT;
      } else {
	 r = G_TYPE_STRING;
      }
   } 
   return r;
} 


void coot::restraints_editor::fill_angle_tree_data(GtkWidget *restraints_editor_dialog,
						   const coot::dictionary_residue_restraints_t &restraints) {
   
   GtkWidget *angles_treeview = lookup_widget(restraints_editor_dialog, 
					      "angles_treeview");
   
   GtkTreeView *tv_angles = GTK_TREE_VIEW(angles_treeview);
   GtkTreeStore *tree_store_angles =
      gtk_tree_store_new (5, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			  G_TYPE_FLOAT, G_TYPE_FLOAT);
   
   view_and_store_angles.store = tree_store_angles;
   view_and_store_angles.view  = tv_angles;
   GtkTreeIter   toplevel;
   
   gtk_tree_view_set_model(tv_angles, GTK_TREE_MODEL(tree_store_angles));
   
   for (unsigned int i=0; i<restraints.angle_restraint.size(); i++) { 
      gtk_tree_store_append(tree_store_angles, &toplevel, NULL);
      gtk_tree_store_set(tree_store_angles, &toplevel,
			 0, restraints.angle_restraint[i].atom_id_1_4c().c_str(),
			 1, restraints.angle_restraint[i].atom_id_2_4c().c_str(),
			 2, restraints.angle_restraint[i].atom_id_3_4c().c_str(),
			 3, restraints.angle_restraint[i].angle(),
			 4, restraints.angle_restraint[i].esd(),
			 -1);
   }

   int tree_type = coot::restraints_editor::TREE_TYPE_ANGLES;
   add_cell_renderer(tv_angles, tree_store_angles, "Atom Name 1", 0, tree_type);
   add_cell_renderer(tv_angles, tree_store_angles, "Atom Name 2", 1, tree_type);
   add_cell_renderer(tv_angles, tree_store_angles, "Atom Name 3", 2, tree_type);
   add_cell_renderer(tv_angles, tree_store_angles, "Angle",       3, tree_type);
   add_cell_renderer(tv_angles, tree_store_angles, "ESD",         4, tree_type);
}

void
coot::restraints_editor::fill_info_tree_data(GtkWidget *restraints_editor_dialog,
			 const coot::dictionary_residue_restraints_t &restraints) {

   GtkWidget *info_treeview = lookup_widget(restraints_editor_dialog, 
					    "info_treeview");
   
   GtkTreeView *tv_info = GTK_TREE_VIEW(info_treeview);
   GtkTreeStore *tree_store_info =
      gtk_tree_store_new (7, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			  G_TYPE_INT, G_TYPE_INT, G_TYPE_STRING);

   view_and_store_info.view = tv_info;
   view_and_store_info.store = tree_store_info;
   
   GtkTreeIter   toplevel;
   
   gtk_tree_view_set_model(tv_info, GTK_TREE_MODEL(tree_store_info));
   std::string tlc = restraints.residue_info.three_letter_code;
   if (tlc.length() == 0) {
      std::cout << "WARNING:: three_letter_code blank/unset." << std::endl;
      std::string subcomp = restraints.residue_info.comp_id;
      tlc = restraints.residue_info.comp_id;
      if (tlc.length() > 3)
	 tlc = tlc.substr(0,3);
      std::cout << "WARNING:: resetting three_letter_code to " << tlc << std::endl;
   } 
   
   gtk_tree_store_append(tree_store_info, &toplevel, NULL);
   gtk_tree_store_set(tree_store_info, &toplevel,
		      0, restraints.residue_info.comp_id.c_str(),
		      1, tlc.c_str(),
		      2, restraints.residue_info.name.c_str(),
		      3, restraints.residue_info.group.c_str(),
		      4, restraints.residue_info.number_atoms_all,
		      5, restraints.residue_info.number_atoms_nh,
		      6, restraints.residue_info.description_level.c_str(),
		      -1);

   int tree_type = coot::restraints_editor::TREE_TYPE_INFO;
   add_cell_renderer(tv_info, tree_store_info, "Comp ID",       0, tree_type);
   add_cell_renderer(tv_info, tree_store_info, "3LetCode",      1, tree_type);
   add_cell_renderer(tv_info, tree_store_info, "Name                  ", 2, tree_type);
   add_cell_renderer(tv_info, tree_store_info, "Group",         3, tree_type);
   add_cell_renderer(tv_info, tree_store_info, "# Non-H Atoms", 4, tree_type);
   add_cell_renderer(tv_info, tree_store_info, "# H Atoms",     5, tree_type);
   add_cell_renderer(tv_info, tree_store_info, "Desc Lev",      6, tree_type);
}

void
coot::restraints_editor::fill_bond_tree_data(GtkWidget *restraints_editor_dialog,
			 const coot::dictionary_residue_restraints_t &restraints) {

   GtkWidget *bonds_treeview = lookup_widget(restraints_editor_dialog, "bonds_treeview");
   GtkTreeView *tv_bonds = GTK_TREE_VIEW(bonds_treeview);
   GtkTreeStore *tree_store_bonds =
      gtk_tree_store_new (5, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			  G_TYPE_FLOAT, G_TYPE_FLOAT);
   view_and_store_bonds.view = tv_bonds;
   view_and_store_bonds.store = tree_store_bonds;

   GtkTreeIter   toplevel;
   gtk_tree_view_set_model(tv_bonds, GTK_TREE_MODEL(tree_store_bonds));
   
   for (unsigned int i=0; i<restraints.bond_restraint.size(); i++) { 
      gtk_tree_store_append(tree_store_bonds, &toplevel, NULL);
      try {

	 // std::cout << "here in fill_bond_tree_data() " << i << " " << restraints.bond_restraint[i] << std::endl;
	 // 
	 // we still want to see the atom names, even if there is no value_esd() or value_dist() (which can throw runtime
	 // errors if they are not present (reasonably enough)
	 //
	 double value_esd = 0.0;
	 double value_dist = 0.0;

	 try {
	    value_dist = restraints.bond_restraint[i].value_dist();
	    value_esd  = restraints.bond_restraint[i].value_esd();
	 } 
	 catch (const std::runtime_error &rte) { } // it's OK if these are not updated.
	 
	 gtk_tree_store_set(tree_store_bonds, &toplevel,
			    0, restraints.bond_restraint[i].atom_id_1_4c().c_str(),
			    1, restraints.bond_restraint[i].atom_id_2_4c().c_str(),
			    2, restraints.bond_restraint[i].type().c_str(),
			    3, value_dist,
			    4, value_esd,
			    -1);
      }
      catch (const std::runtime_error &rte) {
	 
	 // do nothing, it's not really an error if the dictionary
	 // doesn't have target geometry (the bonding description came
	 // from a Chemical Component Dictionary entry for example).
	 //
	 // but we don't want to set the store in that case.

	 std::cout << "caught rte: " << rte.what() << std::endl;
      } 
   }

   int tree_type = coot::restraints_editor::TREE_TYPE_BONDS;
   add_cell_renderer(tv_bonds, tree_store_bonds, "Atom Name 1", 0, tree_type);
   add_cell_renderer(tv_bonds, tree_store_bonds, "Atom Name 2", 1, tree_type);
   add_cell_renderer(tv_bonds, tree_store_bonds, "Type",        2, tree_type);
   add_cell_renderer(tv_bonds, tree_store_bonds, "Bond Length", 3, tree_type);
   add_cell_renderer(tv_bonds, tree_store_bonds, "ESD",         4, tree_type);
}

void
coot::restraints_editor::fill_torsion_tree_data(GtkWidget *restraints_editor_dialog,
			 const coot::dictionary_residue_restraints_t &restraints) {

   GtkWidget *torsions_treeview = lookup_widget(restraints_editor_dialog, 
						"torsions_treeview");

   GtkTreeView *tv_torsions = GTK_TREE_VIEW(torsions_treeview);
   GtkTreeStore *tree_store_torsions =
      gtk_tree_store_new (8, G_TYPE_STRING,
			  G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			  G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_INT);
   view_and_store_torsions.view = tv_torsions;
   view_and_store_torsions.store = tree_store_torsions;
   
   GtkTreeIter   toplevel;
   gtk_tree_view_set_model(tv_torsions, GTK_TREE_MODEL(tree_store_torsions));

   for (unsigned int i=0; i<restraints.torsion_restraint.size(); i++) { 
       gtk_tree_store_append(tree_store_torsions, &toplevel, NULL);
       gtk_tree_store_set(tree_store_torsions, &toplevel,
			  0, restraints.torsion_restraint[i].id().c_str(),
			  1, restraints.torsion_restraint[i].atom_id_1_4c().c_str(),
			  2, restraints.torsion_restraint[i].atom_id_2_4c().c_str(),
			  3, restraints.torsion_restraint[i].atom_id_3_4c().c_str(),
			  4, restraints.torsion_restraint[i].atom_id_4_4c().c_str(),
			  5, restraints.torsion_restraint[i].angle(),
			  6, restraints.torsion_restraint[i].esd(),
			  coot::restraints_editor::TORSION_COL_PERIODICTY, restraints.torsion_restraint[i].periodicity(),
			  -1);
   }

   int tree_type = coot::restraints_editor::TREE_TYPE_TORSIONS;
    add_cell_renderer(tv_torsions, tree_store_torsions, "Tors ID",     0, tree_type);
    add_cell_renderer(tv_torsions, tree_store_torsions, "Atom Name 1", 1, tree_type);
    add_cell_renderer(tv_torsions, tree_store_torsions, "Atom Name 2", 2, tree_type);
    add_cell_renderer(tv_torsions, tree_store_torsions, "Atom Name 3", 3, tree_type);
    add_cell_renderer(tv_torsions, tree_store_torsions, "Atom Name 4", 4, tree_type);
    add_cell_renderer(tv_torsions, tree_store_torsions, "Torsion",     5, tree_type);
    add_cell_renderer(tv_torsions, tree_store_torsions, "ESD",         6, tree_type);
    add_cell_renderer(tv_torsions, tree_store_torsions, "Period",      7, tree_type);
}

void
coot::restraints_editor::fill_chiral_tree_data(GtkWidget *restraints_editor_dialog,
			 const coot::dictionary_residue_restraints_t &restraints) {

   GtkWidget *chirals_treeview = lookup_widget(restraints_editor_dialog, 
					     "chirals_treeview");

   GtkTreeView *tv_chirals = GTK_TREE_VIEW(chirals_treeview);
   GtkTreeStore *tree_store_chirals =
      gtk_tree_store_new (6, G_TYPE_STRING,
			  G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			  G_TYPE_STRING);
   view_and_store_chirals.view = tv_chirals;
   view_and_store_chirals.store = tree_store_chirals;
   
   GtkTreeIter   toplevel;
   gtk_tree_view_set_model(tv_chirals, GTK_TREE_MODEL(tree_store_chirals));

   for (unsigned int i=0; i<restraints.chiral_restraint.size(); i++) { 
      std::string chiral_volume_string = 
	  make_chiral_volume_string(restraints.chiral_restraint[i].volume_sign);
      gtk_tree_store_append(tree_store_chirals, &toplevel, NULL);
      gtk_tree_store_set(tree_store_chirals, &toplevel,
			 0, restraints.chiral_restraint[i].Chiral_Id().c_str(),
			 1, restraints.chiral_restraint[i].atom_id_c_4c().c_str(),
			 2, restraints.chiral_restraint[i].atom_id_1_4c().c_str(),
			 3, restraints.chiral_restraint[i].atom_id_2_4c().c_str(),
			 4, restraints.chiral_restraint[i].atom_id_3_4c().c_str(),
			 coot::restraints_editor::CHIRAL_COL_SIGN, chiral_volume_string.c_str(),
			 -1);
   }

   int tree_type = coot::restraints_editor::TREE_TYPE_CHIRALS;
   add_cell_renderer(tv_chirals, tree_store_chirals, "Chrial ID",    0, tree_type);
   add_cell_renderer(tv_chirals, tree_store_chirals, "Centre Atom ", 1, tree_type);
   add_cell_renderer(tv_chirals, tree_store_chirals, "Atom 1",       2, tree_type);
   add_cell_renderer(tv_chirals, tree_store_chirals, "Atom 2",       3, tree_type);
   add_cell_renderer(tv_chirals, tree_store_chirals, "Atom 3",       4, tree_type);
   add_cell_renderer(tv_chirals, tree_store_chirals, "Sign",         coot::restraints_editor::CHIRAL_COL_SIGN, tree_type);
}

void
coot::restraints_editor::fill_plane_tree_data(GtkWidget *restraints_editor_dialog,
			  const coot::dictionary_residue_restraints_t &restraints) {
   
   GtkWidget *planes_treeview = lookup_widget(restraints_editor_dialog, 
					     "planes_treeview");

   GtkTreeView *tv_planes = GTK_TREE_VIEW(planes_treeview);
   
    // parse lanes to get this
    max_number_of_atoms_in_plane = -1;
    for (unsigned int iplane=0; iplane<restraints.plane_restraint.size(); iplane++) {
       if (restraints.plane_restraint[iplane].n_atoms() > max_number_of_atoms_in_plane)
	  max_number_of_atoms_in_plane = restraints.plane_restraint[iplane].n_atoms();
    } 

    if (max_number_of_atoms_in_plane > 0) { 
       GtkTreeStore *tree_store_planes = make_tree_store_for_planes(max_number_of_atoms_in_plane);
       if (tree_store_planes) { 
       
	  view_and_store_planes.view = tv_planes;
	  view_and_store_planes.store = tree_store_planes;
	  
	  GtkTreeIter   toplevel;
	  gtk_tree_view_set_model(tv_planes, GTK_TREE_MODEL(tree_store_planes));
	  int esd_col_no = max_number_of_atoms_in_plane + 1;

	  for (unsigned int iplane=0; iplane<restraints.plane_restraint.size(); iplane++) {
	     gtk_tree_store_append(tree_store_planes, &toplevel, NULL);

	     // HACK HACK! // FIXME
  	     gtk_tree_store_set(tree_store_planes, &toplevel,
  				esd_col_no, restraints.plane_restraint[iplane].dist_esd(0),
  				-1);
 	     gtk_tree_store_set(tree_store_planes, &toplevel,
 				0, restraints.plane_restraint[iplane].plane_id.c_str(),
 				-1);
	     for (int iat=0; iat<restraints.plane_restraint[iplane].n_atoms(); iat++) { 
 		gtk_tree_store_set(tree_store_planes, &toplevel,
 				   iat+1, restraints.plane_restraint[iplane][iat].first.c_str(),
 				   -1);
	     }
	  }

	  add_plane_cell_renderer(tv_planes, tree_store_planes, "Plane ID", 0,
				  coot::restraints_editor::TREE_TYPE_PLANES,
				  max_number_of_atoms_in_plane);
	  int col_no = 1;
	  for (int i=1; i<=max_number_of_atoms_in_plane; i++) {
	     std::string atom_n_str = "Atom ";
	     atom_n_str += coot::util::int_to_string(i);
	     add_plane_cell_renderer(tv_planes, tree_store_planes, atom_n_str.c_str(), col_no,
				     coot::restraints_editor::TREE_TYPE_PLANES,
				     max_number_of_atoms_in_plane);
	     col_no++;
	  }
	  add_plane_cell_renderer(tv_planes, tree_store_planes, "ESD", esd_col_no,
				  coot::restraints_editor::TREE_TYPE_PLANES,
				  max_number_of_atoms_in_plane);
       }
    }
}

GtkTreeStore *
coot::restraints_editor::make_tree_store_for_planes(int natoms) {

   GtkTreeStore *tree_store_planes_local = NULL;
   if (natoms == 3)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, 
			     G_TYPE_FLOAT);
   if (natoms == 4)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 5)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 6)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 7)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 8)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 9)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 10)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 11)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 12)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 13)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 14)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 15)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 16)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 17)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 18)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, 
			     G_TYPE_FLOAT);
   if (natoms == 19)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, 
			     G_TYPE_FLOAT);
   if (natoms == 20)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 21)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 22)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 23)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 24)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT);
   if (natoms == 25)
      tree_store_planes_local =
	 gtk_tree_store_new (natoms+2, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_STRING,
			     G_TYPE_FLOAT);
   return tree_store_planes_local;
} 

coot::dictionary_residue_restraints_t
coot::restraints_editor::make_restraint() const {

   coot::dictionary_residue_restraints_t rest("s", 0);
   // read the tree stores: 

   std::vector<coot::dict_bond_restraint_t> bonds = get_bond_restraints();
   std::vector<coot::dict_angle_restraint_t> angles = get_angle_restraints();
   std::vector<coot::dict_torsion_restraint_t> torsions = get_torsion_restraints();
   std::vector<coot::dict_chiral_restraint_t> chirals = get_chiral_restraints();
   std::vector<coot::dict_plane_restraint_t> planes = get_plane_restraints();
   std::pair<bool, std::vector <coot::dict_atom> > atom_info = get_atom_info();
   std::pair<bool, coot::dict_chem_comp_t> residue_info = get_residue_info();

   if (residue_info.first) { 
      rest.residue_info      = residue_info.second;
      rest.atom_info         = atom_info.second;
      rest.bond_restraint    = bonds;
      rest.angle_restraint   = angles;
      rest.torsion_restraint = torsions;
      rest.chiral_restraint  = chirals;
      rest.assign_chiral_volume_targets();
      rest.plane_restraint   = planes;

      bool p_c_state = atom_info.first;
      rest.set_has_partial_charges(p_c_state);
   }

   return rest;
}

std::vector<coot::dict_bond_restraint_t> 
coot::restraints_editor::get_bond_restraints() const {

   std::vector<coot::dict_bond_restraint_t> r;
   // e.g. tree_store_bonds:
   GtkTreeIter  iter;

   bool v = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(view_and_store_bonds.store), &iter);
   while (v) {
      std::string atom1("");
      std::string atom2("");
      std::string type("");
      float dist = -1.0;
      float esd  = -1.0;
      for (int col_no=0; col_no<5; col_no++) {
	 int col_type = get_column_type(coot::restraints_editor::TREE_TYPE_BONDS, col_no, -1);
	 if (col_type == G_TYPE_STRING) { 
	    gchar *place_string_here = NULL;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_bonds.store), &iter, col_no, &place_string_here, -1);
	    // std::cout << "got string :" << place_string_here << ": from col_no " << col_no << std::endl;
	    if (col_no == 0)
	       if (place_string_here)
		  atom1 = place_string_here;
	    if (col_no == 1)
	       if (place_string_here)
		  atom2 = place_string_here;
	    if (col_no == 2)
	       if (place_string_here)
		  type = place_string_here;
	 }
	 if (col_type == G_TYPE_FLOAT) {
	    float f;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_bonds.store), &iter, col_no, &f, -1);
	    // std::cout << "got float :" << f << ": from col_no " << col_no << std::endl;
	    if (col_no == 3)
	       dist = f;
	    if (col_no == 4)
	       esd = f;
	 }
      }
      v = gtk_tree_model_iter_next (GTK_TREE_MODEL(view_and_store_bonds.store), &iter);
      if ((atom1.length() > 0) &&
	  (atom2.length() > 0) &&
 	  (type.length() > 0)  &&
 	  (dist > 0.0)         &&
 	  (esd > 0.0)) { 
 	 coot::dict_bond_restraint_t rest(atom1, atom2, type, dist, esd, 0.0, 0.0, false);
// 	 std::cout << "added a bond restraint ";
// 	 std::cout << ":" << atom1 << ": :" << atom2 << ": :" << type << ": " << dist << " " << esd << std::endl;

 	 r.push_back(rest);
      }
   }
   return r;
} 

std::vector<coot::dict_angle_restraint_t> 
coot::restraints_editor::get_angle_restraints() const {
   std::vector<coot::dict_angle_restraint_t> r;
   // e.g. view_and_store_angles.store:
   GtkTreeIter  iter;

   bool v = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(view_and_store_angles.store), &iter);
   while (v) {
      std::string atom1("");
      std::string atom2("");
      std::string atom3("");
      float angle = -1.0;
      float esd   = -1.0;
      for (int col_no=0; col_no<5; col_no++) {
	 int col_type = get_column_type(coot::restraints_editor::TREE_TYPE_ANGLES, col_no, -1);
	 if (col_type == G_TYPE_STRING) { 
	    gchar *place_string_here;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_angles.store), &iter, col_no, &place_string_here, -1);
	    // std::cout << "got string :" << place_string_here << ": from col_no " << col_no << std::endl;
	    if (col_no == 0)
	       atom1 = place_string_here;
	    if (col_no == 1)
	       atom2 = place_string_here;
	    if (col_no == 2)
	       atom3 = place_string_here;
	 }
	 if (col_type == G_TYPE_FLOAT) {
	    float f;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_angles.store), &iter, col_no, &f, -1);
	    // std::cout << "got float :" << f << ": from col_no " << col_no << std::endl;
	    if (col_no == 3)
	       angle = f;
	    if (col_no == 4)
	       esd = f;
	 }
      }
      v = gtk_tree_model_iter_next (GTK_TREE_MODEL(view_and_store_angles.store), &iter);
      if ((atom1.length() > 0) &&
	  (atom2.length() > 0) &&
	  (atom3.length() > 0) &&
 	  (angle > 0.0)  &&
 	  (esd   > 0.0)) { 
 	 coot::dict_angle_restraint_t rest(atom1, atom2, atom3, angle, esd);
 	 // std::cout << "added a angle restraint ";
	 // std::cout << ":" << atom1 << ": :" << atom2 << ": :" << atom3 << ": " << angle << " " << esd << std::endl;

 	 r.push_back(rest);
      }
   }

   return r;
}

std::vector<coot::dict_torsion_restraint_t> 
coot::restraints_editor::get_torsion_restraints() const {
   std::vector<coot::dict_torsion_restraint_t> r;
   // e.g. view_and_store_torsions.store:
   GtkTreeIter  iter;

   // Perhaps there were no plane restraints?
   //
   if (! view_and_store_torsions.store)
      return r;

   bool v = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(view_and_store_torsions.store), &iter);
   while (v) {
      std::string torsion_id("");
      std::string atom1("");
      std::string atom2("");
      std::string atom3("");
      std::string atom4("");
      float torsion = -999.0;
      float esd   = -1.0;
      int period = -1;
      for (int col_no=0; col_no<8; col_no++) {
	 int col_type = get_column_type(coot::restraints_editor::TREE_TYPE_TORSIONS, col_no, -1);
	 if (col_type == G_TYPE_STRING) { 
	    gchar *place_string_here;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_torsions.store), &iter, col_no, &place_string_here, -1);
	    // std::cout << "got string :" << place_string_here << ": from col_no " << col_no << std::endl;
	    if (col_no == 0)
	       torsion_id = place_string_here;
	    if (col_no == 1)
	       atom1 = place_string_here;
	    if (col_no == 2)
	       atom2 = place_string_here;
	    if (col_no == 3)
	       atom3 = place_string_here;
	    if (col_no == 4)
	       atom4 = place_string_here;
	 }
	 if (col_type == G_TYPE_FLOAT) {
	    float f;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_torsions.store), &iter, col_no, &f, -1);
	    // std::cout << "got float :" << f << ": from col_no " << col_no << std::endl;
	    if (col_no == 5)
	       torsion = f;
	    if (col_no == 6)
	       esd = f;
	 }
	 if (col_type == G_TYPE_INT) {
	    int ii;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_torsions.store), &iter, col_no, &ii, -1);
	    if (col_no == coot::restraints_editor::TORSION_COL_PERIODICTY)
	       period = ii;
	 }
      }
      v = gtk_tree_model_iter_next (GTK_TREE_MODEL(view_and_store_torsions.store), &iter);

      if ((atom1.length() > 0) &&
	  (atom2.length() > 0) &&
	  (atom3.length() > 0) &&
	  (atom4.length() > 0) &&
	  (torsion_id.length() > 0) &&
 	  (torsion > -990)     &&
 	  (esd   > -1.0)        && 
	  (period > -1)) { 
 	 coot::dict_torsion_restraint_t rest(torsion_id, atom1, atom2, atom3, atom4, torsion, esd, period);
 	 // std::cout << "added a torsion restraint ";
	 // std::cout << ":" << torsion_id << ": :" << atom1 << ": :" << atom2 << ": :" << atom3 << ": " << atom4 << ": "
	 // << torsion << " " << esd << " " << period << std::endl;

 	 r.push_back(rest);
      }
   }

   return r;
}

// This should be replaced when plane restraints are a tree.
// 
std::vector<coot::dict_plane_restraint_t>
coot::restraints_editor::get_plane_restraints() const {

   std::vector<coot::dict_plane_restraint_t> v;
   int mpa = max_number_of_atoms_in_plane; 
   int n_cols = mpa + 2 ; // plane_id and esd

   // Perhaps there were no plane restraints?
   //
   if (! view_and_store_planes.store)
      return v;

   GtkTreeIter  iter;
   bool b = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(view_and_store_planes.store), &iter);
   while (b) {
      std::string plane_id("");
      std::vector<std::string> atoms;
      float esd = - 1.0; // test for non negative at end
      for (int col_no=0; col_no<n_cols; col_no++) {
	 int col_type = get_column_type(coot::restraints_editor::TREE_TYPE_PLANES, col_no, mpa);
//  	 std::cout << "col_no " << col_no << " of " << n_cols << " is of type "
//  		   << col_type << std::endl;
	 if (col_type == G_TYPE_STRING) {
	    gchar *place_string_here;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_planes.store), &iter, col_no, &place_string_here, -1);
	    if (place_string_here) { 
	       if (col_no == 0) {
		  plane_id = place_string_here;
	       } else {
		  atoms.push_back(place_string_here);
	       }
	    }
	 }
	 if (col_type == G_TYPE_FLOAT) {
	    float val;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_planes.store), &iter, col_no, &val, -1);
	    esd = val;
	 }
      }

      bool made_restraint = 0; 
      if (plane_id != "") {
	 if (atoms.size() > 0) {
	    if (esd > 0.0) {
	       made_restraint = 1;
	       coot::dict_plane_restraint_t r(plane_id, atoms, esd); 
	       v.push_back(r);
	    } 
	 } 
      }
      if (made_restraint) {
	 // std::cout << "Made restraint from :" << plane_id << ": " << atoms.size() << " atoms "
	 // << "with esd " << esd << std::endl;
      } else {
	 // std::cout << "No restraint from :" << plane_id << ": " << atoms.size() << " atoms "
	 // << "with esd " << esd << std::endl;
      } 
      b = gtk_tree_model_iter_next (GTK_TREE_MODEL(view_and_store_planes.store), &iter);
   }
   return v;
}

std::vector<coot::dict_chiral_restraint_t> 
coot::restraints_editor::get_chiral_restraints() const {
   std::vector<coot::dict_chiral_restraint_t> r;
   // e.g. view_and_store_chirals.store:
   GtkTreeIter  iter;

   bool v = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(view_and_store_chirals.store), &iter);
   while (v) {
      std::string chiral_id("");
      std::string atomc("");
      std::string atom1("");
      std::string atom2("");
      std::string atom3("");
      int sign = -1111;
      for (int col_no=0; col_no<6; col_no++) {
	 int col_type = get_column_type(coot::restraints_editor::TREE_TYPE_CHIRALS, col_no, -1);
	 if (col_type == G_TYPE_STRING) { 
	    gchar *place_string_here;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_chirals.store), &iter, col_no, &place_string_here, -1);
	    // std::cout << "got string :" << place_string_here << ": from col_no " << col_no << std::endl;
	    if (col_no == 0)
	       chiral_id = place_string_here;
	    if (col_no == 1)
	       atomc = place_string_here;
	    if (col_no == 2)
	       atom1 = place_string_here;
	    if (col_no == 3)
	       atom2 = place_string_here;
	    if (col_no == 4)
	       atom3 = place_string_here;
	    // new style string (not int) 20091019
	    if (col_no == coot::restraints_editor::CHIRAL_COL_SIGN)
	       sign = chiral_volume_string_to_chiral_sign(place_string_here);
	 }
      }
      v = gtk_tree_model_iter_next (GTK_TREE_MODEL(view_and_store_chirals.store), &iter);

      if ((atom1.length() > 0) &&
	  (atom2.length() > 0) &&
	  (atom3.length() > 0) &&
	  (atomc.length() > 0) &&
	  (chiral_id.length() > 0) &&
	  (sign > -999)) { 
 	 coot::dict_chiral_restraint_t rest(chiral_id, atomc, atom1, atom2, atom3, sign);
//  	 std::cout << "added a chiral restraint ";
// 	 std::cout << ":" << chiral_id << ": :" << atomc << ": :" << atom1 << ": :"
// 		   << atom2 << ": " << atom3 << ": " << sign << std::endl;

 	 r.push_back(rest);
      }
   }

   return r;
}

std::string
coot::restraints_editor::make_chiral_volume_string(int chiral_sign) const {
  
  std::string s = coot::protein_geometry::make_chiral_volume_string(chiral_sign);
  return s;
}
 
int 
coot::restraints_editor::chiral_volume_string_to_chiral_sign(const std::string &chiral_vol_string) const {

  int i = coot::protein_geometry::chiral_volume_string_to_chiral_sign(chiral_vol_string);
  return i;
}


std::pair<bool, std::vector <coot::dict_atom> >
coot::restraints_editor::get_atom_info() const {

   std::vector <coot::dict_atom> r;
   bool have_partial_charges_flag = 0; 
   
   // e.g. view_and_store_atoms.store:
   GtkTreeIter  iter;

   bool v = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(view_and_store_atoms.store), &iter);
   while (v) {
      std::string atom_name("");
      std::string atom_element("");
      std::string energy_type("");
      float partial_charge = -100.0; // not used now
      int formal_charge = 0;
      for (int col_no=0; col_no<4; col_no++) {
	 int col_type = get_column_type(coot::restraints_editor::TREE_TYPE_ATOMS, col_no, -1);
	 if (col_type == G_TYPE_STRING) { 
	    gchar *place_string_here;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_atoms.store), &iter, col_no, &place_string_here, -1);
	    // std::cout << "got string :" << place_string_here << ": from col_no " << col_no << std::endl;
	    if (col_no == 0)
	       atom_name = place_string_here;
	    if (col_no == 1)
	       atom_element = place_string_here;
	    if (col_no == 2)
	       energy_type = place_string_here;
	 }
	 if (col_type == G_TYPE_FLOAT) {
	    // float f;
	    // gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_atoms.store), &iter, col_no, &f, -1);
            int i;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_atoms.store), &iter, col_no, &i, -1);
	    // if (col_no == 3)
            // partial_charge = f;
            if (col_no == 3)
               formal_charge = i;
	 }
      }
      v = gtk_tree_model_iter_next (GTK_TREE_MODEL(view_and_store_atoms.store), &iter);

      if ((atom_name.length()    > 0) &&
	  (atom_element.length() > 0) &&
 	  (energy_type.length()  > 0)) {
	 std::pair<bool, float> part_charge_pair(false, partial_charge); // invalid
	 if (partial_charge > -99.9) {
	    part_charge_pair.first = true; // make valid
	    have_partial_charges_flag = 1;
	 }

         std::pair<bool, int> formal_charge_pair(true, formal_charge);
 	 coot::dict_atom info(atom_name, atom_name, atom_element, energy_type, part_charge_pair);
         info.formal_charge = formal_charge_pair; // not part of a constructor. Maybe change that in future.
//  	 std::cout << "added a atom info ";
// 	 std::cout << ":" << atom_name << ": :" << atom_element << ": :" << energy_type << ": " << partial_charge
// 		   << " " << part_charge_pair.first << std::endl;

 	 r.push_back(info);
      }
   }
   return std::pair<bool, std::vector<coot::dict_atom> > (have_partial_charges_flag, r);
}


std::pair<bool, coot::dict_chem_comp_t>
coot::restraints_editor::get_residue_info() const {

   coot::dict_chem_comp_t info;
   bool proper = 0;

   // e.g. view_and_store_atoms.store:
   GtkTreeIter  iter;

   // Currently with an input mol file, Acedrg writes out dictionaries with no proper name
   // (uses . for name):
   //
   // when the dialog is filled then, there is no string added to the name field.
   // when we read it we get a blank and hence the "Incomprehensible comp-id' message.
   // We don't want that.
   // 

   bool v = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(view_and_store_info.store), &iter);
   if (v) { 
      std::string comp_id("");
      std::string tlc("");
      std::string name("");
      std::string group("");
      std::string description_level("not-set");
      int n_atoms = -1;
      int n_H_atoms = -1;
   
      for (int col_no=0; col_no<7; col_no++) {
	 int col_type = get_column_type(coot::restraints_editor::TREE_TYPE_INFO, col_no, -1);
	 if (col_type == G_TYPE_STRING) { 
	    gchar *place_string_here;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_info.store), &iter, col_no,
			       &place_string_here, -1);
	    if (col_no == 0)
	       comp_id = place_string_here;
	    if (col_no == 1)
	       tlc = place_string_here;
	    if (col_no == 2)
	       name = place_string_here;
	    if (col_no == 3)
	       group = place_string_here;
	    if (col_no == 6)
	       description_level = place_string_here;
	 }
	 if (col_type == G_TYPE_INT) {
	    int ii;
	    gtk_tree_model_get(GTK_TREE_MODEL(view_and_store_info.store), &iter, col_no, &ii, -1);
	    if (col_no == 4)
	       n_atoms = ii;
	    if (col_no == 5)
	       n_H_atoms = ii;
	 }
      }

      if (tlc.length() == 0) {
	 std::cout << "WARNING:: get_residue_info() three_letter_code blank/unset." << std::endl;
	 std::cout << "WARNING:: get_residue_info() resetting three_letter_code to "
		   << comp_id << std::endl;
	 tlc = comp_id;
      }
      
      if ((comp_id.length() > 0) &&
	  (tlc.length()     > 0) &&
	  // (name.length()    > 0) && // see comments above.
	  // 	  (group.length()   > 0) && // libcheck from SMILES has group of "." in cif file, '
	                                    // which after the cif reader is "".
	  (description_level != "not-set") &&
	  (n_atoms > 0)          &&
	  (n_H_atoms > 0)) {

// 	 std::cout << "DEBUG:: makeing a coot::dict_chem_comp_t with comp_id "
// 		   << comp_id << std::endl;
	 
	 coot::dict_chem_comp_t res_info(comp_id, tlc, name, group, n_atoms,
					 n_H_atoms, description_level);
	 proper = true;
	 info = res_info;
      } else {
	 std::cout << "WARNING:: Incomprehensible chem_comp!\n";
	 std::cout << "comp_id :" << comp_id << ":, tlc :" << tlc << ":, name :" << name
		   << ":, group :" << group << ": n_atoms: "
		   << n_atoms << " n_H_atoms: " << n_H_atoms << " desc-level :"
		   << description_level << ":" << std::endl;
      }
   }

   // std::cout << "debug:: in get_residue_info() comp_id is " << info.comp_id << std::endl;

   return std::pair<bool, coot::dict_chem_comp_t> (proper, info);
}


void apply_restraint_by_widget(GtkWidget *w) {

   // std::cout << "Make dictionary reisdue restraints container here" << std::endl;

   int imol_restraints = coot::protein_geometry::IMOL_ENC_ANY; // FIXME

   // which restraints editor has dialog w?
   graphics_info_t g;
   coot::restraints_editor re = g.get_restraints_editor(w);
   if (re.is_valid()) {
      coot::dictionary_residue_restraints_t r = re.make_restraint();
      // do something with r.
      std::string filename = "coot-tmp-restraints.cif";
      r.write_cif(filename);
      coot::protein_geometry *pg = g.Geom_p();
      std::string type = r.residue_info.comp_id;
      bool v = pg->replace_monomer_restraints(type, imol_restraints, r);
      g.redraw_molecules_with_residue(type);
      if (v)
	 std::cout << "INFO:: restraints for \"" << type << "\" were replaced" << std::endl;
      else
	 std::cout << "INFO:: restraints for \"" << type << "\" were added " << std::endl;
   }
}

void restraints_editor_save_restraint_by_widget(GtkWidget *w) {

   graphics_info_t g;
   coot::restraints_editor re = g.get_restraints_editor(w);
   if (re.is_valid()) {
      GtkWidget *ww = create_save_restraint_chooserdialog();
      coot::dictionary_residue_restraints_t r = re.make_restraint();
      std::string filename = "monomer-";

      filename += r.residue_info.comp_id;
      filename += ".cif";
#if (GTK_MAJOR_VERSION == 1) || ((GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION < 10))
#else
      gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(w), TRUE);
#endif      
      gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(ww), filename.c_str());
      add_ccp4i_project_optionmenu(ww, COOT_CIF_DICTIONARY_FILE_SELECTION);
      add_filename_filter_button(ww, COOT_CIF_DICTIONARY_FILE_SELECTION);
      coot::dictionary_residue_restraints_t *ptr = new coot::dictionary_residue_restraints_t("", 0);
      *ptr = r;
      g_object_set_data(G_OBJECT(ww), "restraints", (gpointer) ptr);
      gtk_widget_show(ww);
   }
}

void save_monomer_restraints_by_widget(GtkDialog *chooser) {

   // recall the restraints come from reading the entries in the dialog
   // 
   const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(chooser));
   coot::dictionary_residue_restraints_t *t =
      (coot::dictionary_residue_restraints_t *) g_object_get_data (G_OBJECT (chooser), "restraints");
   t->write_cif(filename);
}


void restraints_editor_delete_restraint_by_widget(GtkWidget *w) {

   graphics_info_t g;
   coot::restraints_editor re = g.get_restraints_editor(w);
   if (re.is_valid()) {
      re.delete_restraint(w);
   }
}

void
coot::restraints_editor::delete_restraint(GtkWidget *w) {
   
   //first find the active tab in the notebook.  That will give us the
   //model and view
      
   GtkTreeIter   iter;
   
   GtkWidget *nb = lookup_widget(w, "restraints_editor_notebook");
   GtkNotebook *notebook = GTK_NOTEBOOK(nb);
   gint current_page_index = gtk_notebook_get_current_page(notebook);
   if (current_page_index != -1) { 
      GtkTreeStore *tree_store = get_tree_store_by_notebook_page(current_page_index);
      GtkTreeView *tree_view = get_tree_view_by_notebook_page(current_page_index);
      GtkTreeSelection *tree_selection = gtk_tree_view_get_selection(tree_view);
      if (tree_store) {
	 if (tree_selection) {
	    GtkTreeModel *model = GTK_TREE_MODEL(tree_store);
	    gtk_tree_selection_get_selected(tree_selection, &model, &iter);
	    gtk_tree_store_remove(tree_store, &iter);
	 }
      } 
   }
}

void restraints_editor_add_restraint_by_widget(GtkWidget *w) {
   graphics_info_t g;
   coot::restraints_editor re = g.get_restraints_editor(w);
   if (re.is_valid()) {
      re.add_restraint(w);
   }
} 

void
coot::restraints_editor::add_restraint(GtkWidget *w) {

   //first find the active tab in the notebook.  That will give us the
   //model and view
   GtkTreeIter   iter;
   
   GtkWidget *nb = lookup_widget(w, "restraints_editor_notebook");
   GtkNotebook *notebook = GTK_NOTEBOOK(nb);
   gint current_page_index = gtk_notebook_get_current_page(notebook);
   if (current_page_index != -1) { 
      GtkTreeStore *tree_store = get_tree_store_by_notebook_page(current_page_index);
      GtkTreeView *tree_view = get_tree_view_by_notebook_page(current_page_index);
      GtkTreeSelection *tree_selection = gtk_tree_view_get_selection(tree_view);
      if (tree_store) {
	 if (tree_selection) {
	    GtkTreeModel *model = GTK_TREE_MODEL(tree_store);
	    gtk_tree_selection_get_selected(tree_selection, &model, &iter);

	    GtkTreeIter *parent = NULL;
	    // now add a restraint under the current line:
	    gtk_tree_store_append(tree_store, &iter, parent);
	 }
      }
   }
}

GtkTreeStore *
coot::restraints_editor::get_tree_store_by_notebook_page(gint current_page_index) const {

   GtkTreeStore *tree_store = NULL;
   switch(current_page_index) {
   case(0):
      tree_store = view_and_store_info.store;
      break;
   case(1):
      tree_store = view_and_store_atoms.store;
      break;
   case(2):
      tree_store = view_and_store_bonds.store;
      break;
   case(3):
      tree_store = view_and_store_angles.store;
      break;
   case(4):
      tree_store = view_and_store_torsions.store;
      break;
   case(5):
      tree_store = view_and_store_chirals.store;
      break;
   case(6):
      tree_store = view_and_store_planes.store;
      break;
   }
   return tree_store;
}

GtkTreeView *
coot::restraints_editor::get_tree_view_by_notebook_page(gint current_page_index) const {

   GtkTreeView *tree_view = NULL;

   switch(current_page_index) {
   case(0):
      tree_view = view_and_store_info.view;
      break;
   case(1):
      tree_view = view_and_store_atoms.view;
      break;
   case(2):
      tree_view = view_and_store_bonds.view;
      break;
   case(3):
      tree_view = view_and_store_angles.view;
      break;
   case(4):
      tree_view = view_and_store_torsions.view;
      break;
   case(5):
      tree_view = view_and_store_chirals.view;
      break;
   case(6):
      tree_view = view_and_store_planes.view;
      break;
   }
   return tree_view;
}

