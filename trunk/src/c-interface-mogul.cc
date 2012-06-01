
#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "mogul-interface.hh"

#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface-mogul.hh"

#include "interface.h"

void
undefined_function_markup(int imol, const char *chain_id, int res_no, const char *ins_code, const char *mogul_out_file_name) {

   coot::mogul m;
   m.parse(mogul_out_file_name);
   m.set_max_z_badness(5.0);
   graphics_info_t g;

   if (is_valid_model_molecule(imol)) { 
      CResidue *residue_p = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p == NULL) {
	 std::cout << "WARNING:: no such residue" << std::endl;
      } else { 
	 if (m.n_items() > 0) {
	    show_mogul_geometry_dialog(m);
	    int new_obj = new_generic_object_number("Mogul Validation");
	    PPCAtom residue_atoms = 0;
	    int n_residue_atoms;
	    residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (unsigned int i=0; i<m.n_items(); i++) { 
   	       if (m[i].type == coot::mogul_item::BOND) {
		  // mogul indexes (from sdf-indices) are 1-based, atoms in residue are 0-based.
		  int idx_1 = m[i].idx_1 - 1;
		  int idx_2 = m[i].idx_2 - 1;
		  if (idx_1 >= 0 && idx_1 < n_residue_atoms) { 
		     if (idx_2 >= 0 && idx_2 < n_residue_atoms) {
			CAtom *at_1 = residue_atoms[idx_1];
			CAtom *at_2 = residue_atoms[idx_2];
			std::string hex_colour = m[i].colour();
			to_generic_object_add_line(new_obj, hex_colour.c_str(), 5,
						   at_1->x, at_1->y, at_1->z,
						   at_2->x, at_2->y, at_2->z);
		     }
		  }
	       }
	       
	       if (m[i].type == coot::mogul_item::ANGLE) {
		  // mogul indexes (from sdf-indices) are 1-based, atoms in residue are 0-based.
		  int idx_1 = m[i].idx_1 - 1;
		  int idx_2 = m[i].idx_2 - 1;
		  int idx_3 = m[i].idx_3 - 1;
		  if (idx_1 >= 0 && idx_1 < n_residue_atoms) { 
		     if (idx_2 >= 0 && idx_2 < n_residue_atoms) {
			if (idx_3 >= 0 && idx_3 < n_residue_atoms) {
			   CAtom *at_1 = residue_atoms[idx_1];
			   CAtom *at_2 = residue_atoms[idx_2];
			   CAtom *at_3 = residue_atoms[idx_3];
			   std::string hex_colour = m[i].colour();
			   try { 
			      coot::arc_info_type angle_info(at_1, at_2, at_3);
			      to_generic_object_add_arc(new_obj, hex_colour.c_str(), 5,
							angle_info.start,
							angle_info.end,
							angle_info.start_point.x(),
							angle_info.start_point.y(),
							angle_info.start_point.z(),
							angle_info.start_dir.x(),
							angle_info.start_dir.y(),
							angle_info.start_dir.z(),
							angle_info.normal.x(),
							angle_info.normal.y(),
							angle_info.normal.z());
			   }
			   catch (std::runtime_error rte) {
			      std::cout << "WARNING:: " << rte.what() << std::endl;
			   }
			}
		     }
		  }
	       }
	    }
	    set_display_generic_object(new_obj, 1);
	    graphics_draw();
	 }
      }
   }
} 


// The mogul output file refers to a specific residue (we need the
// residue to get the atom order).  The restraints can then apply to
// all residues of that type.
// 
int update_restraints_using_undefined_function(int imol, const char *chain_id, int res_no, const char *ins_code,
					       const char *monomer_type, const char *mogul_out_file_name) {

   int s = 0;
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) { 
      CResidue *residue_p = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) { 
	 coot::mogul m(mogul_out_file_name);
	 coot::dictionary_residue_restraints_t new_restraints =
	    m.make_restraints(residue_p, monomer_type, *g.Geom_p());
	 s = g.Geom_p()->replace_monomer_restraints_conservatively(monomer_type, new_restraints);
      }
   }
   return s;
}

void show_mogul_geometry_dialog(const coot::mogul &m) {

   GtkWidget *w = wrapped_create_mogul_geometry_dialog(m);
   if (w)
      gtk_widget_show(w);
}

GtkWidget *wrapped_create_mogul_geometry_dialog(const coot::mogul &m) {

   GtkWidget *w = create_mogul_geometry_dialog();
   // fill w here.

   GtkTreeView *mogul_bonds_treeview    = GTK_TREE_VIEW(lookup_widget(w, "mogul_bonds_treeview"));
   GtkTreeView *mogul_angles_treeview   = GTK_TREE_VIEW(lookup_widget(w, "mogul_angles_treeview"));
   GtkTreeView *mogul_torsions_treeview = GTK_TREE_VIEW(lookup_widget(w, "mogul_torsions_treeview"));

//    GtkTreeStore *tree_store_atoms = gtk_tree_store_new (,
// 							G_TYPE_STRING, G_TYPE_STRING, // atom names
// 							G_TYPE_FLOAT, // value from model
// 							G_TYPE_FLOAT, // median
// 							G_TYPE_FLOAT, // std_dev

// 							);

   return w;
}

