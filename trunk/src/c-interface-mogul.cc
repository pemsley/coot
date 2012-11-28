
#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "coot-utils.hh"
#include "mogul-interface.hh"

#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface-mogul.hh"

#include "interface.h"

void
mogul_markup(int imol, const char *chain_id, int res_no, const char *ins_code, const char *mogul_out_file_name) {

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
	    show_mogul_geometry_dialog(m, residue_p);
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
int update_restraints_using_mogul(int imol, const char *chain_id, int res_no, const char *ins_code,
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

void show_mogul_geometry_dialog(const coot::mogul &m, CResidue *residue) {

   if (graphics_info_t::use_graphics_interface_flag) { 
#ifdef HAVE_GOOCANVAS
      GtkWidget *w = wrapped_create_mogul_geometry_dialog(m, residue);
      if (w) { 
	 gtk_widget_show(w);
      } 
#endif // HAVE_GOOCANVAS
   }
}

GtkCellRenderer *
coot::mogul_results_add_cell_renderer(GtkTreeView *tree_view,
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

   gtk_tree_view_column_set_sort_column_id(column, pos);

   char *s = new char[column_title.length() + 1];
   strcpy(s, column_title.c_str());
   g_object_set_data (G_OBJECT (cell_renderer), "column", GINT_TO_POINTER (pos));
   g_object_set_data (G_OBJECT (cell_renderer), "tree_type", GINT_TO_POINTER (tree_type));
   return cell_renderer;
}


GtkWidget
*wrapped_create_mogul_geometry_dialog(const coot::mogul &m, CResidue *residue) {

   GtkWidget *w = create_mogul_geometry_dialog();

   if (residue) { 

      // fill w here.

      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      residue->GetAtomTable(residue_atoms, n_residue_atoms);
      
      GtkTreeView *mogul_bonds_treeview    = GTK_TREE_VIEW(lookup_widget(w, "mogul_bonds_treeview"));
      GtkTreeView *mogul_angles_treeview   = GTK_TREE_VIEW(lookup_widget(w, "mogul_angles_treeview"));
      GtkTreeView *mogul_torsions_treeview = GTK_TREE_VIEW(lookup_widget(w, "mogul_torsions_treeview"));

      coot::fill_mogul_bonds_tab(      mogul_bonds_treeview, w, m, residue);
      coot::fill_mogul_angles_tab(    mogul_angles_treeview, w, m, residue);
      coot::fill_mogul_torsions_tab(mogul_torsions_treeview, w, m, residue);
   }
   coot::goograph *goograph = new coot::goograph;
   g_object_set_data(G_OBJECT(w), "goograph", goograph);
   return w;
}

void
coot::fill_mogul_bonds_tab(GtkTreeView *mogul_bonds_treeview, GtkWidget *dialog,
			   const coot::mogul &m, CResidue *r) {

   // We want to see: atom-name-1 atom-name-2 value mean median std-dev z
   // 
   GtkTreeStore *tree_store_bonds = gtk_tree_store_new(7, G_TYPE_STRING, G_TYPE_STRING,
						       G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT,
						       G_TYPE_FLOAT, G_TYPE_FLOAT);

   GtkTreeView *tv_bonds = GTK_TREE_VIEW(mogul_bonds_treeview);
   gtk_tree_view_set_model(tv_bonds, GTK_TREE_MODEL(tree_store_bonds));
   GtkTreeIter   toplevel;

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   r->GetAtomTable(residue_atoms, n_residue_atoms);

   for (unsigned int i=0; i<m.n_items(); i++) {
      const coot::mogul_item &item = m[i];
      if (item.type == coot::mogul_item::BOND) {
	 int idx_1 = m[i].idx_1-1;
	 int idx_2 = m[i].idx_2-1;
	 CAtom *at_1 = residue_atoms[idx_1];
	 CAtom *at_2 = residue_atoms[idx_2];
	 std::string atom_name_1 = at_1->name;
	 std::string atom_name_2 = at_2->name;
	 
	 gtk_tree_store_append(tree_store_bonds, &toplevel, NULL);
	 gtk_tree_store_set(tree_store_bonds, &toplevel,
			    0, atom_name_1.c_str(),
			    1, atom_name_2.c_str(),
			    2, m[i].value,
			    3, m[i].mean,
			    4, m[i].median,
			    5, m[i].std_dev,
			    6, m[i].z,
			    -1);
      }
   }

   int tree_type = 0; // coot::mogul::TREE_TYPE_BONDS;
   coot::mogul_results_add_cell_renderer(tv_bonds, tree_store_bonds, "Atom 1", 0, tree_type);
   coot::mogul_results_add_cell_renderer(tv_bonds, tree_store_bonds, "Atom 2", 1, tree_type);
   coot::mogul_results_add_cell_renderer(tv_bonds, tree_store_bonds, "Value",  2, tree_type);
   coot::mogul_results_add_cell_renderer(tv_bonds, tree_store_bonds, "Mean",   3, tree_type);
   coot::mogul_results_add_cell_renderer(tv_bonds, tree_store_bonds, "Median", 4, tree_type);
   coot::mogul_results_add_cell_renderer(tv_bonds, tree_store_bonds, "ESD",    5, tree_type);
   coot::mogul_results_add_cell_renderer(tv_bonds, tree_store_bonds, "z",      6, tree_type);

   GtkTreeSelection *sel_bonds = gtk_tree_view_get_selection(GTK_TREE_VIEW(mogul_bonds_treeview));
   g_signal_connect(sel_bonds,  "changed", (GCallback) coot::on_mogul_bonds_selection_changed,  dialog);
   coot::mogul *mcp = new coot::mogul(m);
   coot::minimol::residue *mmres_p = new coot::minimol::residue(r);
   g_object_set_data(G_OBJECT(sel_bonds), "mogul", mcp); 
   g_object_set_data(G_OBJECT(sel_bonds), "residue", mmres_p); 
} 

void
coot::on_mogul_bonds_selection_changed(GtkTreeSelection *treeselection,
				       gpointer          user_data) {

   GtkTreeIter  iter;
   GtkTreeModel *model; 
   std::string altconf = "";
   gboolean r = gtk_tree_selection_get_selected(treeselection, &model, &iter);
   if (r) {
      gchar *atom_id_1, *atom_id_2;
      gchar *value, *mean, *median, *esd, *z;
      coot::mogul *mogul_p = static_cast<mogul *>(g_object_get_data(G_OBJECT(treeselection), "mogul"));
      coot::minimol::residue *mmres_p = static_cast<minimol::residue *>(g_object_get_data(G_OBJECT(treeselection), "residue"));
      gtk_tree_model_get(model, &iter,
			 0, &atom_id_1,
			 1, &atom_id_2,
			 2, &value,
			 3, &mean,
			 4, &median,
			 5, &esd,
			 6, &z,
			 -1);
      bool ifound_1 = false;
      bool ifound_2 = false;
      clipper::Coord_orth pos_1(1,0,1);
      clipper::Coord_orth pos_2(0,0,1);
      double x, y, zz;
      for (unsigned int i=0; i<mmres_p->n_atoms(); i++) {
	 const minimol::atom &at = mmres_p->atoms[i];
	 if (at.name == atom_id_1) {
	    pos_1 == at.pos;
	    x  = at.pos.x();
	    y  = at.pos.y();
	    zz = at.pos.z();
	    pos_1 = clipper::Coord_orth(x, y, zz);
	    ifound_1 = true;
	 } 
	 if (at.name == atom_id_2) {
	    x  = at.pos.x();
	    y  = at.pos.y();
	    zz = at.pos.z();
	    pos_2 = clipper::Coord_orth(x, y, zz);
	    ifound_2 = true;
	 } 
      }
      if (ifound_1 && ifound_2) {
	 clipper::Coord_orth p(0.5*(pos_1.x() + pos_2.x()),
			       0.5*(pos_1.y() + pos_2.y()),
			       0.5*(pos_1.z() + pos_2.z()));
	 set_rotation_centre(p.x(), p.y(), p.z());

	 std::vector<std::string> atom_ids;
	 atom_ids.push_back(atom_id_1);
	 atom_ids.push_back(atom_id_2);

	 GtkWidget *mogul_geometry_dialog = static_cast<GtkWidget *> (user_data);
	 if (mogul_geometry_dialog) {
	    update_mogul_histogram_dialog(mogul_geometry_dialog, *mogul_p, atom_ids, mmres_p, altconf);
	 } else {
	    std::cout << "null mogul_geometry_dialog" << std::endl;
	 }
      }
   }
}

void
coot::fill_mogul_angles_tab(GtkTreeView *mogul_angles_treeview, GtkWidget *dialog,
			    const coot::mogul &m, CResidue *r) {

   // We want to see: atom-name-1 atom-name-2 atom-name-3 value mean median std-dev z
   // 
   GtkTreeStore *tree_store_angles = gtk_tree_store_new(8, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
							G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT,
							G_TYPE_FLOAT, G_TYPE_FLOAT);

   GtkTreeView *tv_angles = GTK_TREE_VIEW(mogul_angles_treeview);
   gtk_tree_view_set_model(tv_angles, GTK_TREE_MODEL(tree_store_angles));
   GtkTreeIter   toplevel;

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   r->GetAtomTable(residue_atoms, n_residue_atoms);

   for (unsigned int i=0; i<m.n_items(); i++) {
      const coot::mogul_item &item = m[i];
      if (item.type == coot::mogul_item::ANGLE) {
	 int idx_1 = m[i].idx_1-1;
	 int idx_2 = m[i].idx_2-1;
	 int idx_3 = m[i].idx_3-1;
	 CAtom *at_1 = residue_atoms[idx_1];
	 CAtom *at_2 = residue_atoms[idx_2];
	 CAtom *at_3 = residue_atoms[idx_3];
	 std::string atom_name_1 = at_1->name;
	 std::string atom_name_2 = at_2->name;
	 std::string atom_name_3 = at_3->name;
	 
	 gtk_tree_store_append(tree_store_angles, &toplevel, NULL);
	 gtk_tree_store_set(tree_store_angles, &toplevel,
			    0, atom_name_1.c_str(),
			    1, atom_name_2.c_str(),
			    2, atom_name_3.c_str(),
			    3, m[i].value,
			    4, m[i].mean,
			    5, m[i].median,
			    6, m[i].std_dev,
			    7, m[i].z,
			    -1);
      }
   }

   int tree_type = 1; // coot::mogul::TREE_TYPE_ANGLES
   coot::mogul_results_add_cell_renderer(tv_angles, tree_store_angles, "Atom 1", 0, tree_type);
   coot::mogul_results_add_cell_renderer(tv_angles, tree_store_angles, "Atom 2", 1, tree_type);
   coot::mogul_results_add_cell_renderer(tv_angles, tree_store_angles, "Atom 3", 2, tree_type);
   coot::mogul_results_add_cell_renderer(tv_angles, tree_store_angles, "Value",       3, tree_type);
   coot::mogul_results_add_cell_renderer(tv_angles, tree_store_angles, "Mean",        4, tree_type);
   coot::mogul_results_add_cell_renderer(tv_angles, tree_store_angles, "Median",      5, tree_type);
   coot::mogul_results_add_cell_renderer(tv_angles, tree_store_angles, "ESD",         6, tree_type);
   coot::mogul_results_add_cell_renderer(tv_angles, tree_store_angles, "z",           7, tree_type);

   GtkTreeSelection *sel_angles = gtk_tree_view_get_selection(GTK_TREE_VIEW(mogul_angles_treeview));
   g_signal_connect(sel_angles, "changed", (GCallback) coot::on_mogul_angles_selection_changed, dialog);
   coot::mogul *mcp = new coot::mogul(m);
   coot::minimol::residue *mmres_p = new coot::minimol::residue(r);
   g_object_set_data(G_OBJECT(sel_angles), "mogul", mcp); 
   g_object_set_data(G_OBJECT(sel_angles), "residue", mmres_p); 

}

void
coot::on_mogul_angles_selection_changed(GtkTreeSelection *treeselection,
					gpointer          user_data) {

   GtkTreeIter  iter;
   GtkTreeModel *model; 
   gboolean r = gtk_tree_selection_get_selected(treeselection, &model, &iter);
   if (r) {
      gchar *atom_id_1, *atom_id_2, *atom_id_3;
      gchar *value, *mean, *median, *esd, *z;
      coot::mogul *mogul_p = static_cast<mogul *>(g_object_get_data(G_OBJECT(treeselection), "mogul"));
      coot::minimol::residue *mmres_p =
	 static_cast<minimol::residue *>(g_object_get_data(G_OBJECT(treeselection), "residue"));
      std::string altconf = "";
      gtk_tree_model_get(model, &iter,
			 0, &atom_id_1,
			 1, &atom_id_2,
			 2, &atom_id_3,
			 3, &value,
			 4, &mean,
			 5, &median,
			 6, &esd,
			 7, &z,
			 -1);
      bool ifound_1 = false;
      bool ifound_2 = false;
      bool ifound_3 = false;
      clipper::Coord_orth pos_1(1,0,1);
      clipper::Coord_orth pos_2(0,0,1);
      clipper::Coord_orth pos_3(0,0,1);
      double x, y, zz;
      for (unsigned int i=0; i<mmres_p->n_atoms(); i++) {
	 const minimol::atom &at = mmres_p->atoms[i];
	 if (at.name == atom_id_1) {
	    pos_1 == at.pos;
	    x  = at.pos.x();
	    y  = at.pos.y();
	    zz = at.pos.z();
	    pos_1 = clipper::Coord_orth(x, y, zz);
	    ifound_1 = true;
	 } 
	 if (at.name == atom_id_2) {
	    x  = at.pos.x();
	    y  = at.pos.y();
	    zz = at.pos.z();
	    pos_2 = clipper::Coord_orth(x, y, zz);
	    ifound_2 = true;
	 } 
	 if (at.name == atom_id_3) {
	    x  = at.pos.x();
	    y  = at.pos.y();
	    zz = at.pos.z();
	    pos_3 = clipper::Coord_orth(x, y, zz);
	    ifound_3 = true;
	 } 
      }
      if (ifound_1 && ifound_2 && ifound_3) {
	 double scale = 0.3333;
	 double x_ = scale * (pos_3.x() + pos_2.x() + pos_1.x());
	 double y_ = scale * (pos_3.y() + pos_2.y() + pos_1.y());
	 double z_ = scale * (pos_3.z() + pos_2.z() + pos_1.z());
	 if (0) 
	    std::cout << "atoms: "
		      << pos_1.format() << "\n       "
		      << pos_2.format() << "\n       "
		      << pos_3.format() << " ----> " << x_ << " " << y_ << " " << z_
		      << "\n";
	 set_rotation_centre(x_, y_, z_);
	 std::vector<std::string> atom_ids;
	 atom_ids.push_back(atom_id_1);
	 atom_ids.push_back(atom_id_2);
	 atom_ids.push_back(atom_id_3);
	 GtkWidget *mogul_geometry_dialog = static_cast<GtkWidget *> (user_data);
						
	 if (mogul_geometry_dialog) {
	    update_mogul_histogram_dialog(mogul_geometry_dialog, *mogul_p, atom_ids, mmres_p, altconf);
	 } else {
	    std::cout << "null mogul_geometry_dialog" << std::endl;
	 } 
      }
   }
}


void
coot::fill_mogul_torsions_tab(GtkTreeView *mogul_torsions_treeview, GtkWidget *dialog,
			      const coot::mogul &m, CResidue *r) {
   
   // We want to see: atom-name-1 atom-name-2 atom-name-3 atom-name-4 value mean median std-dev z
   // 
   GtkTreeStore *tree_store_torsions = gtk_tree_store_new(9,
							G_TYPE_STRING, G_TYPE_STRING,
							G_TYPE_STRING, G_TYPE_STRING,
							G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT,
							G_TYPE_FLOAT, G_TYPE_FLOAT);

   GtkTreeView *tv_torsions = GTK_TREE_VIEW(mogul_torsions_treeview);
   gtk_tree_view_set_model(tv_torsions, GTK_TREE_MODEL(tree_store_torsions));
   GtkTreeIter   toplevel;

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   r->GetAtomTable(residue_atoms, n_residue_atoms);

   for (unsigned int i=0; i<m.n_items(); i++) {
      const coot::mogul_item &item = m[i];
      if (item.type == coot::mogul_item::TORSION) {
	 int idx_1 = m[i].idx_1-1;
	 int idx_2 = m[i].idx_2-1;
	 int idx_3 = m[i].idx_3-1;
	 int idx_4 = m[i].idx_4-1;
	 CAtom *at_1 = residue_atoms[idx_1];
	 CAtom *at_2 = residue_atoms[idx_2];
	 CAtom *at_3 = residue_atoms[idx_3];
	 CAtom *at_4 = residue_atoms[idx_4];
	 std::string atom_name_1 = at_1->name;
	 std::string atom_name_2 = at_2->name;
	 std::string atom_name_3 = at_3->name;
	 std::string atom_name_4 = at_4->name;
	 
	 gtk_tree_store_append(tree_store_torsions, &toplevel, NULL);
	 gtk_tree_store_set(tree_store_torsions, &toplevel,
			    0, atom_name_1.c_str(),
			    1, atom_name_2.c_str(),
			    2, atom_name_3.c_str(),
			    3, atom_name_4.c_str(),
			    4, m[i].value,
			    5, m[i].mean,
			    6, m[i].median,
			    7, m[i].std_dev,
			    8, m[i].z,
			    -1);
      }
   }

   int tree_type = 2; // coot::mogul::TREE_TYPE_TORSION;
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "Atom Name 1", 0, tree_type);
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "Atom Name 2", 1, tree_type);
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "Atom Name 3", 2, tree_type);
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "Atom Name 4", 3, tree_type);
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "Value",       4, tree_type);
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "Mean",        5, tree_type);
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "Median",      6, tree_type);
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "ESD",         7, tree_type);
   coot::mogul_results_add_cell_renderer(tv_torsions, tree_store_torsions, "z",           8, tree_type);

}


void
coot::update_mogul_histogram_dialog(GtkWidget *geometry_dialog,
				    const mogul &m, const std::vector<std::string> &atom_ids,
				    coot::minimol::residue *r, const std::string &altconf) { 

#ifdef HAVE_GOOCANVAS   

   int ifound = 0;
      
   // bonds
   if (atom_ids.size() == 2) {
      std::vector<int> indices(2,0);
      for (unsigned int i=0; i<r->n_atoms(); i++) {
	 coot::minimol::atom &at = (*r)[i];
	 std::string atom_name = at.name;
	 std::string atom_alt_conf = at.altLoc;
	 for (unsigned int j=0; j<atom_ids.size(); j++) { 
	    if (atom_name == atom_ids[j]) {
	       if (altconf == atom_alt_conf) { 
		  indices[j] = i+1; // the mogul atom index
		  ifound++;
	       }
	    }
	 }
      }

      if (ifound == 2) {
	 try {
	    coot::mogul_item item = m.get_bond_item(indices);
	    coot::goograph *gg = static_cast<coot::goograph *> (g_object_get_data(G_OBJECT(geometry_dialog),
										  "goograph"));
	    std::string title = "Bond distribution ";
	    title += atom_ids[0];
	    title += " - ";
	    title += atom_ids[1];
	    mogul_histogram_for_item(gg, item, "Bond length", title);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: " << rte.what() << std::endl;
	 } 
      }
   }

   // angles
   if (atom_ids.size() == 3) {

      std::vector<int> indices(3,0);
      for (unsigned int i=0; i<r->n_atoms(); i++) {
	 coot::minimol::atom &at = (*r)[i];
	 std::string atom_name = at.name;
	 std::string atom_alt_conf = at.altLoc;
	 for (unsigned int j=0; j<atom_ids.size(); j++) { 
	    if (atom_name == atom_ids[j]) {
	       if (altconf == atom_alt_conf) { 
		  indices[j] = i+1; // the mogul atom index
		  ifound++;
	       }
	    }
	 }
      }

      if (ifound == 3) {
	 try {
	    coot::mogul_item item = m.get_angle_item(indices);
	    coot::goograph *gg = static_cast<coot::goograph *> (g_object_get_data(G_OBJECT(geometry_dialog),
										  "goograph"));

	    std::string title = "Angle distribution -";
	    title += atom_ids[0];
	    title += " - ";
	    title += atom_ids[1];
	    title += " - ";
	    title += atom_ids[2];
	    title += " - ";
	    mogul_histogram_for_item(gg, item, "Angle", title);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: " << rte.what() << std::endl;
	 } 
      } 
   }
#endif // HAVE_GOOCANVAS   
}

#ifdef HAVE_GOOCANVAS
void
coot::mogul_histogram_for_item(coot::goograph *gg, const coot::mogul_item &item,
			       const std::string &x_axis_label, const std::string &title) {
   
	    
   if (gg) { 
      std::vector<std::pair<double, double> > data;
      for (unsigned int i=0; i<item.distribution.counts.size(); i++) {
	 std::pair<double, double> p(item.distribution.bin_start + i*item.distribution.bin_width,
				     item.distribution.counts[i]);
	 data.push_back(p);
      }

      if (data.size()) { 
	 double min_x = +1e20;
	 double max_x = -1e20;
	 double min_y = +1e20;
	 double max_y = -1e20;
	 for (unsigned int i=0; i<data.size(); i++) { 
	    if (data[i].first < min_x)
	       min_x = data[i].first;
	    if (data[i].first > max_x)
	       max_x = data[i].first;
	    if (data[i].second < min_y)
	       min_y = data[i].second;
	    if (data[i].second > max_y)
	       max_y = data[i].second;
	 }

	 gg->clear_traces_and_annotations();
	 int trace_id = gg->trace_new();

	 gg->set_data(trace_id, data);
	 gg->set_trace_type(trace_id, coot::graph_trace_info_t::PLOT_TYPE_BAR);
	 gg->set_plot_title(title);
	 gg->set_axis_label(coot::goograph::X_AXIS, x_axis_label);
	 gg->set_axis_label(coot::goograph::Y_AXIS, "Counts");
	 gg->set_extents(coot::goograph::Y_AXIS, 0, max_y);

	 if (item.value < min_x)
	    gg->set_extents(coot::goograph::X_AXIS, 0.99*item.value, max_x);
	 if (item.value > max_x)
	    gg->set_extents(coot::goograph::X_AXIS, min_x, 1.01*item.value);

	 bool dashed = true;
	 double p3x = item.value - 0.01 * (max_x - min_x);
	 lig_build::pos_t p1(item.value, 0);
	 lig_build::pos_t p2(item.value, max_y);
	 lig_build::pos_t p3(p3x, min_y+0.95*(max_y-min_y));
	 gg->add_annotation_line(p1, p2, "#880000", 2, dashed, true, false);
	 std::string s = "Value from \nmodel: ";
	 s += coot::util::float_to_string_using_dec_pl(item.value, 2);
	 gg->add_annotation_text(s, p3, "#880000", "");

	 gg->draw_graph();
	 gg->show_dialog();
      }
   } 
} 
#endif // HAVE_GOOCANVAS
