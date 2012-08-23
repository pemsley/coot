/* lbg/lbg.cc
 * 
 * Copyright 2010, 2011, 2012 by The University of Oxford
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

#ifdef HAVE_GOOCANVAS

#ifdef USE_PYTHON
#include <Python.h> // this is here get round header warnings
#endif

#include <sys/types.h>  // for stating
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <algorithm>

#ifdef MAKE_ENTERPRISE_TOOLS
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <boost/python.hpp>
#endif

#include <cairo.h>
#if CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#include "lbg.hh"
#include "lbg-drag-and-drop.hh"
#include "qed-interface.hh" // interface to silicos-it biscu-it python function

#if ( ( (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 11) ) || GTK_MAJOR_VERSION > 2)
// 
lbg_info_t *
lbg(lig_build::molfile_molecule_t mm,
    std::pair<bool, coot::residue_spec_t> ligand_spec_pair,
    CMMDBManager *mol,
    const std::string &view_name,
    const std::string &molecule_file_name,
    int imol,
    bool use_graphics_interface_flag,
    bool stand_alone_flag_in,
    int (*get_url_func_pointer_in) (const char *s1, const char *s2),
    void (*prodrg_import_function_pointer) (std::string file_name),
    void (*sbase_import_function_pointer) (std::string comp_id),
    std::string (*get_drug_mdl_file_function_pointer_in) (std::string drug_name)) {
   
   lbg_info_t *lbg = NULL; // failure return value.
   std::string glade_file = "lbg.glade";

   std::string glade_file_full = coot::package_data_dir();
   glade_file_full += "/";
   glade_file_full += glade_file;

   bool glade_file_exists = 0;
   struct stat buf;
   int err = stat(glade_file_full.c_str(), &buf);
   if (! err)
      glade_file_exists = 1;

   if (! glade_file_exists) { 
      std::cout << "ERROR:: glade file " << glade_file_full << " not found" << std::endl;
   } else {

      // If we are using the graphics interface then we want non-null from the builder.
      // add_from_file_status should be good or we are in trouble.
      // 
      // If not, we need not call gtk_builder_add_from_file().
      // 

      if (use_graphics_interface_flag) {
	 
	 GtkBuilder *builder = gtk_builder_new ();
	 guint add_from_file_status =
	    gtk_builder_add_from_file (builder, glade_file_full.c_str(), NULL);

	 if (! add_from_file_status) {

	    // Handle error...
	    
	    std::cout << "ERROR:: gtk_builder_add_from_file() \"" << glade_file_full
		      << "\" failed." << std::endl;
	    if (builder) {
	       std::cout << "ERROR:: where builder was non-null" << std::endl;
	    } else {
	       std::cout << "ERROR:: where builder was NULL" << std::endl;
	    }
	    
	 } else {

	    // Happy Path
	 
	    lbg = new lbg_info_t(imol);
	    if (stand_alone_flag_in)
	       lbg->set_stand_alone();
	    int init_status = lbg->init(builder);
	    
	    // normal built-in and using-graphics path...
	    
	    if (! init_status) {
	       std::cout << "ERROR:: lbg init failed." << std::endl;
	    } else {
	       // happy path
	       if (use_graphics_interface_flag)
		  gtk_builder_connect_signals (builder, lbg->canvas);
	       g_object_unref (G_OBJECT (builder));
	       if (ligand_spec_pair.first)
		  lbg->set_ligand_spec(ligand_spec_pair.second);
	       
	       widgeted_molecule_t wmol = lbg->import_mol_file(mm, molecule_file_name, mol);
	       lbg->render_from_molecule(wmol);
	       lbg->update_descriptor_attributes();
	       lbg->save_molecule();
	    }
	 }
      }

      if (! use_graphics_interface_flag) {
	 lbg = new lbg_info_t(imol);
	 lbg->no_graphics_mode(); // sets flag for "don't display a window"
	 int init_status = lbg->init(NULL); // setup canvas, but not gui.
	 if (ligand_spec_pair.first)
	    lbg->set_ligand_spec(ligand_spec_pair.second);
	 
	 widgeted_molecule_t wmol = lbg->import_mol_file(mm, molecule_file_name, mol);
	 lbg->render_from_molecule(wmol);
	 lbg->update_descriptor_attributes();
      }
      
      if (get_url_func_pointer_in != NULL) {
	 lbg->set_curl_function(get_url_func_pointer_in);
      }
      if (prodrg_import_function_pointer) {
	 lbg->set_prodrg_import_function(prodrg_import_function_pointer);
      }
      if (sbase_import_function_pointer) {
	 lbg->set_sbase_import_function(sbase_import_function_pointer);
      }
      if (get_drug_mdl_file_function_pointer_in) {
	 lbg->set_get_drug_mdl_file_function(get_drug_mdl_file_function_pointer_in);
      } 
   }
   return lbg;
}
#endif // GTK_VERSION

	     

GtkWidget *get_canvas_from_scrolled_win(GtkWidget *canvas) {

   return canvas;
}

#ifdef MAKE_ENTERPRISE_TOOLS
// 
RDKit::RWMol
lbg_info_t::rdkit_mol(const widgeted_molecule_t &mol) const {

   RDKit::RWMol m;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

   // we need to account for atoms that are closed.  So we make
   // atom_index[idx] to return the index of the atom in the rdkit
   // molecule.  i.e. converting from (possibly closed) atoms in a
   // lig_build::molecule
   // 
   std::vector<int> atom_index(mol.atoms.size());
   int at_no = 0;

   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      if (! mol.atoms[iat].is_closed()) { 
	 RDKit::Atom *at = new RDKit::Atom;
	 at->setAtomicNum(tbl->getAtomicNumber(mol.atoms[iat].element));

	 // add the name to "at" too here (if you can).
	 // 
	 try {
	    at->setProp("name", std::string(mol.atoms[iat].get_atom_name()));
	    std::string ele_capped =
	       coot::util::capitalise(coot::util::remove_leading_spaces(mol.atoms[iat].element));
	    int atomic_number = tbl->getAtomicNumber(ele_capped);
	    at->setAtomicNum(atomic_number);
	    at->setMass(tbl->getAtomicWeight(atomic_number));
	    m.addAtom(at, false, true); // take ownership
	    std::string dum;
	    at->setProp("lbg_atom_index", coot::util::int_to_string(iat));
	    if (at->hasProp("lbg_atom_index")) {
	       // std::cout << "atom has lbg_atom_index property " << std::endl;
	       at->getProp("lbg_atom_index", dum);
	       // std::cout << "   lbg_atom_index value is " << dum << std::endl;
	    } else {
	       // std::cout << "atom does not have lbg_atom_index property " << std::endl;
	    } 
	    atom_index[iat] = at_no++;
	 }
	 catch (boost::bad_lexical_cast blc) {
	    std::cout << "on making atom bad_lexical_cast " << iat << " "
		      << blc.what() << std::endl;
	 }
	 catch (std::exception rte) {
	    std::cout << "on making atom " << iat << " " << rte.what() << std::endl;
	 }
      } else {
	 // std::cout << "debug:: atom " << iat << " is closed " << std::endl;
      } 
   }

   for (unsigned int ib=0; ib<mol.bonds.size(); ib++) {
      if (! mol.bonds[ib].is_closed()) {
	 RDKit::Bond::BondType type = convert_bond_type(mol.bonds[ib].get_bond_type());
	 RDKit::Bond *bond = new RDKit::Bond(type);
	 int idx_1 = mol.bonds[ib].get_atom_1_index();
	 int idx_2 = mol.bonds[ib].get_atom_2_index();
	 if (!mol.atoms[idx_1].is_closed() && !mol.atoms[idx_2].is_closed()) {
	    int idx_1_rdkit = atom_index[idx_1];
	    int idx_2_rdkit = atom_index[idx_2];
	    bond->setBeginAtomIdx(idx_1_rdkit);
	    bond->setEndAtomIdx(  idx_2_rdkit);
	    if (type == RDKit::Bond::AROMATIC) { 
	       bond->setIsAromatic(true);
	       m[idx_1_rdkit]->setIsAromatic(true);
	       m[idx_2_rdkit]->setIsAromatic(true);
	    } 
	    m.addBond(bond, true); // take ownership
	 }
      } else {
	 // std::cout << "debug:: bond " << ib << " is closed " << std::endl;
      } 
   }
   return m;
}
#endif

#ifdef MAKE_ENTERPRISE_TOOLS
RDKit::Bond::BondType
lbg_info_t::convert_bond_type(const lig_build::bond_t::bond_type_t &t) const {

   // There are lots more in RDKit::Bond::BondType!
   // 
   RDKit::Bond::BondType bt = RDKit::Bond::UNSPECIFIED;
   if (t == lig_build::bond_t::SINGLE_BOND)
      bt = RDKit::Bond::SINGLE;
   if (t == lig_build::bond_t::DOUBLE_BOND)
      bt = RDKit::Bond::DOUBLE;
   if (t == lig_build::bond_t::TRIPLE_BOND)
      bt = RDKit::Bond::TRIPLE;
   return bt;
}
#endif


#ifdef MAKE_ENTERPRISE_TOOLS
// this can throw a std::exception
// 
std::string
lbg_info_t::get_smiles_string_from_mol_rdkit() const {

   RDKit::RWMol rdkm = rdkit_mol(mol);
   coot::rdkit_mol_sanitize(rdkm);
   RDKit::ROMol *rdk_mol_with_no_Hs = RDKit::MolOps::removeHs(rdkm);
   std::string s = RDKit::MolToSmiles(*rdk_mol_with_no_Hs);
   delete rdk_mol_with_no_Hs;
   return s;
}
#endif

void
lbg_info_t::write_mdl_molfile_using_default_file_name() const {

   mol.write_mdl_molfile(mdl_file_name);
}


void
lbg_info_t::untoggle_others_except(GtkToggleToolButton *button_toggled_on) {
   
   std::map<std::string, GtkToggleToolButton *>::const_iterator it;
   for (it=widget_names.begin(); it!=widget_names.end(); it++) {
      if (it->second != button_toggled_on) {
	 if (gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(it->second))) {
	    gtk_toggle_tool_button_set_active(GTK_TOGGLE_TOOL_BUTTON(it->second), FALSE);
	 }
      }
   }
}

/* This handles button presses on the canvas */
static
gboolean
on_canvas_button_press_old (GtkWidget *widget, GdkEventButton *event)
{
   int x_as_int, y_as_int;
   GdkModifierType state;
   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   GtkObject *obj = GTK_OBJECT(widget);
   if (obj) {
      gpointer gp = gtk_object_get_user_data(obj);
      if (gp) {
	 lbg_info_t *l = static_cast<lbg_info_t *> (gp);
	 l->set_mouse_pos_at_click(x_as_int, y_as_int); // save for dragging
	 // std::cout << "mouse_at_click: " << x_as_int << " " << y_as_int << std::endl;
	 
	 if (0)
	    std::cout << "   on click: scale_correction  " << l->mol.scale_correction.first
		      << " " << l->mol.scale_correction.second
		      << " centre_correction " << l->mol.centre_correction << std::endl;

	    
	 if (l->in_delete_mode_p()) { 
	    l->handle_item_delete(event);
	 } else {
	    l->handle_item_add(event);
	 }
      }
   }
   return TRUE;
}

static bool
on_canvas_button_press(GooCanvasItem  *item,
		       GooCanvasItem  *target_item,
		       GdkEventButton *event,
		       gpointer        user_data) {

   int x_as_int = -1, y_as_int = -1;
   GdkModifierType state;
   if (! event) {
      std::cout << "on_canvas_button_press_new() error NULL event!" << std::endl;
   } else { 
      x_as_int = int(event->x);
      y_as_int = int(event->y);
   }

   if (! target_item) { 
      lbg_info_t *l =
	 static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
      if (!l) {
	 std::cout << "on button-click: NULL lbg-info from item" << std::endl;
      } else {
	 // std::cout << "debug:: handle item add " << event << std::endl;
	 if (l->in_delete_mode_p()) {
	    // std::cout << "calling handle_item_delete() " << event << std::endl;
	    l->handle_item_delete(event);
	 } else { 
	    // std::cout << "calling handle_item_add() " << event << std::endl;
	    l->handle_item_add(event);
	 }
      }
      
   } else {

      // std::cout << "on_canvas_button_press_new() non-NULL target_item " << target_item
      // << std::endl;
      
      coot::residue_spec_t *spec_p =
	 (coot::residue_spec_t *) g_object_get_data (G_OBJECT (target_item), "spec");
      lbg_info_t *l = NULL;
      if (item)
	 l = static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
   
#ifdef MAKE_ENTERPRISE_TOOLS      
      if (spec_p) { 
	 std::cout << "clicked on " << *spec_p << std::endl;
	 if (event->type==GDK_2BUTTON_PRESS) { // double click
	    int imol = spec_p->int_user_data;
	    std::cout << "imol: " << imol << std::endl;
	    if (! l) {
	       std::cout << "oops failed to get lbg_info_t pointer" << std::endl;
	    } else {
	       std::pair<bool, coot::residue_spec_t> ligand_spec_pair = l->get_ligand_spec();
	       std::cout << "ligand pair first: " << ligand_spec_pair.first << std::endl;
	       if (ligand_spec_pair.first) {
		  std::cout << " ligand_spec: " << ligand_spec_pair.second << std::endl;
		  std::cout << "residue_spec: " << *spec_p << std::endl;
		  orient_view(imol,  ligand_spec_pair.second, *spec_p);
	       }
	    }
	 } 
      } else {
	 // std::cout << "null spec" << std::endl;
      }
#endif       

      if (!l) {
	 // std::cout << "null lbg-info from target_item" << std::endl;
      } else {

	 // note: it is not an error if l is null, that can happen
	 // when clicking on the atom names for example - they don't
	 // have a lbg-info attached (at the moment).
	 
	 l->set_mouse_pos_at_click(x_as_int, y_as_int); // save for dragging
	 
	 if (0) 
	    std::cout << "   on click: scale_correction  " << l->mol.scale_correction.first
		      << " " << l->mol.scale_correction.second
		      << " centre_correction " << l->mol.centre_correction << std::endl;
      
      
	 if (l->in_delete_mode_p()) { 
	    l->handle_item_delete(event);
	 } else {
	    l->handle_item_add(event);
	 }
      }
      
   } 
   
   return TRUE;
}
   
   

static bool
on_canvas_button_release(GooCanvasItem  *item,
			     GooCanvasItem  *target_item,
			     GdkEventButton *event,
			     gpointer        user_data) {

   // target_item is null (usually?)
   
   if (item) { 
      lbg_info_t *l =
	 static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
      if (l) {
	 l->clear_button_down_bond_addition();
      } else {
	 // this should not happen
	 std::cout << "======= null lbg_info_t pointer " << std::endl;
      } 
   } else {
      // std::cout << "======= null target_item " << std::endl;
   } 
   return TRUE;
}


void
lbg_info_t::clear_button_down_bond_addition() {

   button_down_bond_addition = 0;
   
   // 20120122 why were we saving here? It was giving a double-save on
   // every modification, thus leading to the need to double tap
   // delete to undo a modification...
   // 
   // std::cout << "save_molecule() from clear_button_down_bond_addition() " << std::endl;
   // save_molecule();
} 
   

static gboolean
on_canvas_motion (GtkWidget *widget, GdkEventMotion *event) {
   int x_as_int=0, y_as_int=0;
   //   if (event->is_hint) {
   if (1) {
      GdkModifierType state;
      gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
      GtkObject *obj = GTK_OBJECT(widget);
      if (obj) {
	 gpointer gp = gtk_object_get_user_data(obj);
	 // std::cout << "got gp: " << gp << std::endl;
	 if (gp) { 
	    lbg_info_t *l = static_cast<lbg_info_t *> (gp);
	    l->handle_drag(state, x_as_int, y_as_int);
	 } else {
	    std::cout << "Failed to get gpointer from object on canvas motion"
		      << std::endl;
	 } 
      } else {
	 std::cout << "Failed to get GtkObject from widget on canvas motion"
		   << std::endl;
      }
   }
   return TRUE;
}



static bool
on_canvas_motion_new(GooCanvasItem  *item,
		     GooCanvasItem  *target_item,
		     GdkEventMotion *event,
		     gpointer        user_data) {
   
   int x_as_int = -1, y_as_int = -1;
   GdkModifierType state;
   if (! event) {
      std::cout << "on_canvas_button_press_new() error NULL event!" << std::endl;
   } else { 
      x_as_int = int(event->x);
      y_as_int = int(event->y);
   }

   if (! target_item) { 
      // std::cout << "canvas_motion: NULL target item" << std::endl;
      lbg_info_t *l =
	 static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
      if (!l) {
	 // std::cout << "canvas_motion: null lbg-info from target_item" << std::endl;
      }
   } else {

#ifdef MAKE_ENTERPRISE_TOOLS
      coot::residue_spec_t *spec_p =
	 static_cast<coot::residue_spec_t *> (g_object_get_data (G_OBJECT (target_item), "spec"));
      int *add_rep_handle_p =
	 static_cast<int *> (g_object_get_data (G_OBJECT (target_item), "add_rep_handle"));

      if (spec_p && add_rep_handle_p) {
	 int add_rep_handle = *add_rep_handle_p;
	 // std::cout << "moused over " << *spec_p << std::endl;
	 if (add_rep_handle >= 0) {
	    int imol = spec_p->int_user_data;
	    // std::cout << "show add_rep_handle: " << add_rep_handle << " for "
	    // << *spec_p << std::endl;
	    set_show_additional_representation(imol, add_rep_handle, 1);
	    all_additional_representations_off_except(imol, add_rep_handle, 0);
	 } 
      }
#endif      
      
   }

   if (! item) {
      // std::cout << "canvas_motion: NULL item!" << std::endl;
   } else {
      // std::cout << "canvas_motion:  item:" << item << std::endl;
      lbg_info_t *l =
	 static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
      if (!l) {
	 // std::cout << "canvas_motion: null lbg-info from item!" << std::endl;
      } else {
	 // std::cout << "Yay - got an lbg_info_t pointer " << l << std::endl;
	 guint state = event->state;
	 GdkModifierType g_state = static_cast<GdkModifierType> (state);
	 int x_as_int = int(event->x);
	 int y_as_int = int(event->y);

//  	 std::cout << ":::::::::::: state: " << state 
//  		   << " " << (state & GDK_BUTTON1_MASK)
// 		   << " button_down_bond_addition " << l->button_down_bond_addition_state()
//  		   << std::endl;
	 
	 l->handle_drag(g_state, x_as_int, y_as_int);
      }
   }
   return TRUE;
}


void
lbg_info_t::handle_drag(GdkModifierType state, int x_as_int, int y_as_int) {

   bool highlight_status = item_highlight_maybe(x_as_int, y_as_int);
   if (! highlight_status) {
      if (state & GDK_BUTTON1_MASK) {
	 if (button_down_bond_addition) {
	    rotate_latest_bond(x_as_int, y_as_int);
	 } else { 
	    drag_canvas(x_as_int, y_as_int);
	 }
      }
   } else {
      if (button_down_bond_addition) {
	 if (highlight_data.single_atom()) 
	    extend_latest_bond(); // checks for sensible atom to snap to.
      }
   }
}

void
lbg_info_t::drag_canvas(int mouse_x, int mouse_y) {

   double delta_x = (double(mouse_x) - mouse_at_click.x) * 1;
   double delta_y = (double(mouse_y) - mouse_at_click.y) * 1;

   mouse_at_click.x = mouse_x;
   mouse_at_click.y = mouse_y;

   // GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   // double top_left_left = GOO_CANVAS(canvas)->hadjustment->lower - delta_x;
   // double top_left_top  = GOO_CANVAS(canvas)->vadjustment->lower - delta_y;
   // std::cout << "scrolling to " << top_left_left << " " << top_left_top << std::endl;
   // goo_canvas_scroll_to(GOO_CANVAS(canvas), top_left_left, top_left_top);

   lig_build::pos_t delta(delta_x, delta_y);
   canvas_drag_offset += delta;

   // should santize delta *here* (and check for not 0,0) before
   // passing to these functions.

   if (is_sane_drag(delta)) {
      clear_and_redraw(delta);
   }
}

void
lbg_info_t::clear_and_redraw(const lig_build::pos_t &delta) {

   clear_canvas();
   if (delta.non_zero()) {
      widgeted_molecule_t new_mol = translate_molecule(delta); // and do a canvas update
      translate_residue_circles(delta);
      render_from_molecule(new_mol);
      update_descriptor_attributes();
   } else {
      std::cout << "==== delta is zero path ==== " << std::endl;
      widgeted_molecule_t saved_mol = mol;
      render_from_molecule(saved_mol);
      update_descriptor_attributes();
   }
   draw_all_flev_annotations();
}


widgeted_molecule_t
lbg_info_t::translate_molecule(const lig_build::pos_t &delta) {  // and do a canvas update

   // we can't translate mol, because that gets wiped in render_from_molecule.
   widgeted_molecule_t new_mol = mol;
   new_mol.translate(delta);
   return new_mol;
}

void
lbg_info_t::translate_residue_circles(const lig_build::pos_t &delta) {
   
   for (unsigned int i=0; i<residue_circles.size(); i++) 
      residue_circles[i].pos += delta;

}

bool
lbg_info_t::is_sane_drag(const lig_build::pos_t &delta) const {

   bool sane = 0;
   if (delta.length() > 0.5) { // don't move with a delta of 0,0.
      if (delta.length() < 50 ) { // don't move with an absurd delta.
	 sane = 1;
      }
   }
   return sane;
} 


void
lbg_info_t::rotate_latest_bond(int x_mouse, int y_mouse) {
   
  // If no highlighted data, simply rotate the bond.  (Note here that
   // the canvas item and the atom coordinates are separately rotated
   // (but result in the same position).
   //

   bool debug = 0;

   if (mol.bonds.size() > 0) { 
      if (mol.atoms.size() > 0) {
	 GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));
	 
	 widgeted_bond_t &bond = mol.bonds.back();
	 widgeted_atom_t &atom = mol.atoms.back();

	 double x_cen = penultimate_atom_pos.x;
	 double y_cen = penultimate_atom_pos.y;
	 double x_at = atom.atom_position.x;
	 double y_at = atom.atom_position.y;
	 double theta_rad = atan2(-(y_mouse-y_cen), (x_mouse-x_cen));
	 double theta_target = theta_rad/DEG_TO_RAD;
	 double theta_current_rad = atan2(-(y_at-y_cen), (x_at-x_cen));
	 double theta_current = theta_current_rad/DEG_TO_RAD;
	 double degrees = theta_target - theta_current;
	 
	 lig_build::pos_t new_atom_pos =
	    atom.atom_position.rotate_about(x_cen, y_cen, -degrees);

	 bool is_close_to_atom = mol.is_close_to_non_last_atom(new_atom_pos);

	 if (is_close_to_atom) {

	    std::cout << " rotate prevented -- too close to atom " << std::endl;

	 }  else {

	    // move ultimate_atom_index atom to new_atom_pos
	    atom.atom_position = new_atom_pos;


	    if (penultimate_atom_index == UNASSIGNED_INDEX) {

	       std::cout << "ERROR:: trapped UNASSIGNED_INDEX for penultimate_atom_index in "
			 << "rotate_latest_bond()" << std::endl;
	       
	    } else {

	       bond.rotate_canvas_item(x_cen, y_cen, -degrees);
	       
	    }
	 }
      }
   }
}

void
lbg_info_t::extend_latest_bond() {

   // we need to check that the highlighted atom is not the atom that
   // we just added the bond to, because (1) that will always be the
   // case just after the button down and the mouse moves with the
   // button down and (2) we don't want to draw a bond back to the
   // atom we just added a bond to.

   // use highlight_data

   // std::cout << "in extend_latest_bond() " << std::endl;

   if (mol.bonds.size() > 0) {
      if (mol.atoms.size() > 0) {

	 widgeted_bond_t &bond = mol.bonds.back();
	 widgeted_atom_t &atom = mol.atoms.back();
	 if (highlight_data.has_contents()) { 
	    if (highlight_data.single_atom()) {

	       int atom_index = highlight_data.get_atom_index();
	       if (atom_index != penultimate_atom_index) { 
		  if (atom_index != ultimate_atom_index) {
		     
		     // std::cout << "extend_latest_bond() " << std::endl;
   
		     GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

		     mol.close_atom(ultimate_atom_index, root);
	       
		     // and add a new bond.
		     // 
		     lig_build::bond_t::bond_type_t bt =
			addition_mode_to_bond_type(canvas_addition_mode);
		     widgeted_bond_t b(penultimate_atom_index, atom_index,
				       mol.atoms[penultimate_atom_index],
				       mol.atoms[atom_index],
				       bt, root);
		     // std::cout << "adding bond! " << b << std::endl;
		     mol.add_bond(b); // can reject addition if bond
				      // between the given atom
				      // indices already exists.
		     // std::cout << "Now we have " << mol.n_open_bonds() << " bonds" << std::endl;
		  }
	       }
	    }
	 }
      }
   }
} 


bool
lbg_info_t::item_highlight_maybe(int x, int y) {

   bool highlight_status = 0;
   remove_bond_and_atom_highlighting();
   std::pair<bool, widgeted_bond_t> bond_info = mol.highlighted_bond_p(x, y);
   if (bond_info.first) {
      highlight_status = 1;
      highlight_bond(bond_info.second, in_delete_mode_p());
   } else {
      std::pair<int, widgeted_atom_t> atom_info = mol.highlighted_atom_p(x, y);
      if (atom_info.first != UNASSIGNED_INDEX) {
	 highlight_status = 1;
	 highlight_atom(atom_info.second, atom_info.first, in_delete_mode_p());
      } 
   }
   return highlight_status;
}

void
lbg_info_t::remove_bond_and_atom_highlighting() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   if (highlight_data.has_contents())
      highlight_data.undisplay();
}

// set highlight_data
void
lbg_info_t::highlight_bond(const lig_build::bond_t &bond, bool delete_mode) {

   int i1 = bond.get_atom_1_index();
   int i2 = bond.get_atom_2_index();
   std::pair<int, int> bond_indices(i1, i2);
   lig_build::pos_t A = mol.atoms[i1].atom_position;
   lig_build::pos_t B = mol.atoms[i2].atom_position;

   std::string col = "#20cc20";
   if (delete_mode)
      col = "#D03030";
   
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
   GooCanvasItem *h_line =
      goo_canvas_polyline_new_line(root,
				   A.x, A.y,
				   B.x, B.y,
				   "line-width", 7.0,
				   "stroke-color", col.c_str(),
				   NULL);
   highlight_data = highlight_data_t(h_line, bond_indices, A, B);
}

// set highlight_data
void
lbg_info_t::highlight_atom(const lig_build::atom_t &atom, int atom_index, bool delete_mode) {

   lig_build::pos_t A = atom.atom_position;
   std::string col = "#20cc20";
   if (delete_mode)
      col = "#D03030";
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));
   double width = 10; // atom.name.length() * 3;
   double height = 14;
   double x1 = A.x - width/2;
   double y1 = A.y - height/2;

   // std::cout << x1 << " " << x2 << " "<< width << " " << height << std::endl;

   GooCanvasItem *rect_item =
      goo_canvas_rect_new (root, x1, y1, width, height,
			   "line-width", 1.0,
			   "stroke-color", col.c_str(),
			   NULL);
   highlight_data = highlight_data_t(rect_item, A, atom_index);
}


void
lbg_info_t::handle_item_add(GdkEventButton *event) {

   bool changed_status = 0;

   int x_as_int, y_as_int;
   GdkModifierType state;
   bool shift_is_pressed = 0;
   if (event->state & GDK_SHIFT_MASK)
      shift_is_pressed = 1;

   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   if (canvas_addition_mode == lbg_info_t::PENTAGON)
      changed_status = try_stamp_polygon(5, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::HEXAGON)
      changed_status = try_stamp_polygon(6, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::HEXAGON_AROMATIC)
      changed_status = try_stamp_hexagon_aromatic(x_as_int, y_as_int, shift_is_pressed);
   if (canvas_addition_mode == lbg_info_t::TRIANGLE)
      changed_status = try_stamp_polygon(3, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::SQUARE)
      changed_status = try_stamp_polygon(4, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::HEPTAGON)
      changed_status = try_stamp_polygon(7, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::OCTAGON)
      changed_status = try_stamp_polygon(8, x_as_int, y_as_int, shift_is_pressed, 0);

   if (is_atom_element(canvas_addition_mode)) {
      changed_status = try_change_to_element(canvas_addition_mode);
   }

   if (is_bond(canvas_addition_mode)) {
      changed_status = try_add_or_modify_bond(canvas_addition_mode, x_as_int, y_as_int);
   }

   if (changed_status) {
      save_molecule();
      update_descriptor_attributes();
   }
}

void
lbg_info_t::update_descriptor_attributes() {
#ifdef MAKE_ENTERPRISE_TOOLS
   try {
      RDKit::RWMol rdkm = rdkit_mol(mol);
      coot::rdkit_mol_sanitize(rdkm);
      update_statusbar_smiles_string(rdkm);
      update_qed(rdkm);
      update_alerts(rdkm);
   }
   catch (std::exception e) {
      std::cout << "WARNING::" << e.what() << std::endl;

      // SMILES string
      if (lbg_statusbar) { 
	 std::string status_string;
	 guint statusbar_context_id =
	    gtk_statusbar_get_context_id(GTK_STATUSBAR(lbg_statusbar), status_string.c_str());
	 gtk_statusbar_push(GTK_STATUSBAR(lbg_statusbar),
			    statusbar_context_id,
			    status_string.c_str());
#ifdef MAKE_ENTERPRISE_TOOLS
	 // QED progress bar
	 gtk_label_set_text(GTK_LABEL(lbg_qed_text_label), "");
	 gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(lbg_qed_progressbar), 0);
	 // alerts
	 clear_canvas_alerts();
#endif	 
      }
   }
#endif      
} 

void
lbg_info_t::handle_item_delete(GdkEventButton *event) {

   int x_as_int, y_as_int;
   GdkModifierType state;
   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

   if (highlight_data.has_contents()) {
      if (highlight_data.single_atom()) {
	 mol.close_atom(highlight_data.get_atom_index(), root);
      } else {
	 int ind_1 = highlight_data.get_bond_indices().first;
	 int ind_2 = highlight_data.get_bond_indices().second;
	 int bond_index = mol.get_bond_index(ind_1, ind_2);
	 mol.close_bond(bond_index, root, 1);
      }
      save_molecule();
      std::cout << "in handle_item_delete() calling update_descriptor_attributes() "
		<< std::endl;
      update_descriptor_attributes();
   }
}


// that's canvas_addition_mode (from the button press).
bool
lbg_info_t::is_bond(int addition_mode) const {
   bool r = 0; 
   if (addition_mode == lbg_info_t::ADD_SINGLE_BOND)
      r = 1;
   if (addition_mode == lbg_info_t::ADD_DOUBLE_BOND)
      r = 1;
   if (addition_mode == lbg_info_t::ADD_TRIPLE_BOND)
      r = 1;
   if (addition_mode == lbg_info_t::ADD_STEREO_OUT_BOND)
      r = 1;
   return r;
}

bool
lbg_info_t::is_atom_element(int addition_mode) const {

   bool r = 0;
   if (addition_mode == lbg_info_t::ATOM_C)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_N)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_O)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_S)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_P)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_F)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_CL)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_I)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_BR)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_X)
      r = 1;
   return r;
}

std::string
lbg_info_t::to_element(int addition_mode) const {

   std::string r = "";
   if (addition_mode == lbg_info_t::ATOM_C)
      r = "C";
   if (addition_mode == lbg_info_t::ATOM_N)
      r = "N";
   if (addition_mode == lbg_info_t::ATOM_O)
      r = "O";
   if (addition_mode == lbg_info_t::ATOM_S)
      r = "S";
   if (addition_mode == lbg_info_t::ATOM_P)
      r = "P";
   if (addition_mode == lbg_info_t::ATOM_F)
      r = "F";
   if (addition_mode == lbg_info_t::ATOM_CL)
      r = "Cl";
   if (addition_mode == lbg_info_t::ATOM_I)
      r = "I";
   if (addition_mode == lbg_info_t::ATOM_BR)
      r = "Br";
   if (addition_mode == lbg_info_t::ATOM_X)
      r = atom_X;  // wrong.
   return r;
}


bool
lbg_info_t::try_change_to_element(int addition_element_mode) {

   bool changed_status = 0;
   if (highlight_data.has_contents()) {
      if (highlight_data.single_atom()) {
	 int atom_index = highlight_data.get_atom_index();
	 if (atom_index != UNASSIGNED_INDEX) { 
	    std::string new_ele = to_element(addition_element_mode);
	    changed_status = mol.atoms[atom_index].change_element(new_ele);
	    std::string fc = font_colour(addition_element_mode);
	    if (changed_status) {
	       change_atom_element(atom_index, new_ele, fc);
	    }
	 }
      }
   }
   return changed_status;
}

void
lbg_info_t::change_atom_id_maybe(int atom_index) {

   std::string ele = mol.atoms[atom_index].element;
   std::string fc = font_colour(ele);
   change_atom_element(atom_index, ele, fc);

} 

bool
lbg_info_t::change_atom_element(int atom_index, std::string new_ele, std::string fc) {

   bool changed_status = 0;
   std::vector<int> local_bonds = mol.bonds_having_atom_with_atom_index(atom_index);
   lig_build::pos_t pos = mol.atoms[atom_index].atom_position;

   // 20110410: Old/simple
   // std::string atom_id = mol.make_atom_id_by_using_bonds(new_ele, local_bonds);
   
   bool gl_flag = false; // not a GL render engine
   lig_build::atom_id_info_t atom_id_info = 
      mol.make_atom_id_by_using_bonds(atom_index, new_ele, local_bonds, gl_flag);
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

   changed_status = mol.atoms[atom_index].update_atom_id_maybe(atom_id_info, fc, root);

   if (changed_status) { 
      for (unsigned int ib=0; ib<local_bonds.size(); ib++) {
	 // make a new line for the bond 
	 //
	 int index_1 = mol.bonds[local_bonds[ib]].get_atom_1_index();
	 int index_2 = mol.bonds[local_bonds[ib]].get_atom_2_index();
	 lig_build::atom_t atom_1 = mol.atoms[index_1];
	 lig_build::atom_t atom_2 = mol.atoms[index_2];

	 mol.bonds[local_bonds[ib]].make_new_canvas_item(atom_1, atom_2, root);
      }
   }
   return changed_status;
}


bool
lbg_info_t::try_add_or_modify_bond(int canvas_addition_mode, int x_mouse, int y_mouse) {

   // button_down_bond_addition is check later so that we can
   // distinguish between canvas dragging and new bond rotation.

   bool changed_status = 0;
   if (! highlight_data.has_contents()) {
      if (mol.is_empty()) {
	 // calling this sets penultimate_atom_pos and latest_bond_canvas_item
	 try_stamp_bond_anywhere(canvas_addition_mode, x_mouse, y_mouse);
	 changed_status = 1; // try_stamp_bond_anywhere always modifies
	 button_down_bond_addition = 1;
      }
   } else {
      if (highlight_data.single_atom()) {
	 int atom_index = highlight_data.get_atom_index();
	 if (atom_index != UNASSIGNED_INDEX) {
	    // add_bond_to_atom() sets latest_bond_canvas_item
	    changed_status = add_bond_to_atom(atom_index, canvas_addition_mode);
	    button_down_bond_addition = changed_status;
	 }
      } else {
	 // highlighted item was a bond then.
	 int ind_1 = highlight_data.get_bond_indices().first;
	 int ind_2 = highlight_data.get_bond_indices().second;
	 int bond_index = mol.get_bond_index(ind_1, ind_2);
	 lig_build::bond_t::bond_type_t bt =
	    addition_mode_to_bond_type(canvas_addition_mode);

	 if (0)
	    std::cout << "....... canvas_addition_mode: " << canvas_addition_mode << " " 
		      << "bond type: " << bt << " c.f. " << lbg_info_t::ADD_STEREO_OUT_BOND
		      << " and " << lig_build::bond_t::OUT_BOND
		      << std::endl;
	 
	 if (bond_index != UNASSIGNED_INDEX) {
	    // we need to pass the atoms so that we know if and how to
	    // shorten the bonds (canvas items) to account for atom names.
	    lig_build::atom_t at_1 = mol.atoms[ind_1];
	    lig_build::atom_t at_2 = mol.atoms[ind_2];
	    GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
	    if (canvas_addition_mode == lbg_info_t::ADD_TRIPLE_BOND)
	       mol.bonds[bond_index].change_bond_order(at_1, at_2, 1, root);
	    else {

	       // (bt is OUT_BOND if we want to convert an OUT_BOND to an IN_BOND).
	       // 
	       if (bt == lig_build::bond_t::OUT_BOND) {

		  if (mol.bonds[bond_index].get_bond_type() != lig_build::bond_t::OUT_BOND &&
		      mol.bonds[bond_index].get_bond_type() != lig_build::bond_t::IN_BOND) {
		     

		     // if it is not currently a stereo out bond, convert
		     // it to a stereo out bond. If it is a stereo_out
		     // bond, convert it with change_bond_order().
		     //
		     mol.bonds[bond_index].set_bond_type(bt);
		     mol.bonds[bond_index].make_new_canvas_item(at_1, at_2, root);
		     
		  } else {

		     // was already a stereo bond
// 		     std::cout << "!!!!! was already a stereo bond, with current: "
// 			       << mol.bonds[bond_index].get_bond_type() << std::endl;
		     if (mol.bonds[bond_index].get_bond_type() == lig_build::bond_t::IN_BOND) { 
			highlight_data.swap_bond_indices();
		     }
			 
		     mol.bonds[bond_index].change_bond_order(at_1, at_2, root); // out to in and vice
             		                                                      // versa (with direction change)
		  } 
	       } else {

		  // Not to a wedge bond

		  // conventional/usual change
		  mol.bonds[bond_index].change_bond_order(at_1, at_2, root); // single to double
				  		                             // or visa versa
	       }
	    } 

	    // Now that we have modified this bond, the atoms
	    // compromising the bond may need to have their atom name
	    // changed (e.g.) N -> NH,
	    // 
	    // we do this for both of the atoms of the bond.
	    //
	    change_atom_id_maybe(ind_1);
	    change_atom_id_maybe(ind_2);
	    changed_status = 1;
	 }
      }
   }
   return changed_status;
}

lig_build::bond_t::bond_type_t
lbg_info_t::addition_mode_to_bond_type(int canvas_addition_mode) const {

   lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
   if (canvas_addition_mode == lbg_info_t::ADD_DOUBLE_BOND)
      bt = lig_build::bond_t::DOUBLE_BOND;
   if (canvas_addition_mode == lbg_info_t::ADD_TRIPLE_BOND)
      bt = lig_build::bond_t::TRIPLE_BOND;
   if (canvas_addition_mode == lbg_info_t::ADD_STEREO_OUT_BOND)
      bt = lig_build::bond_t::OUT_BOND;
   return bt;
}

bool
lbg_info_t::add_bond_to_atom(int atom_index, int canvas_addition_mode) {

   bool changed_status = 0;
   std::vector<int> bonds = mol.bonds_having_atom_with_atom_index(atom_index);

   switch (bonds.size()) {

   case 0:
      add_bond_to_atom_with_0_neighbours(atom_index, canvas_addition_mode);
      changed_status = 1;
      break;

   case 1:
      add_bond_to_atom_with_1_neighbour(atom_index, canvas_addition_mode, bonds[0]);
      changed_status = 1;
      break;

   case 2:
      add_bond_to_atom_with_2_neighbours(atom_index, canvas_addition_mode, bonds);
      changed_status = 1;
      break;

   case 3:
      add_bond_to_atom_with_3_neighbours(atom_index, canvas_addition_mode, bonds);
      changed_status = 1;
      break;

   default:
      std::cout << "not handled yet: " << bonds.size() << " neighbouring bonds" << std::endl;
   }
   return changed_status;
}


std::string
lbg_info_t::font_colour(int addition_element_mode) const {
   std::string font_colour = "black";
   if (addition_element_mode == lbg_info_t::ATOM_O)
      font_colour = "red";
   if (addition_element_mode == lbg_info_t::ATOM_N)
      font_colour = "blue";
   if (addition_element_mode == lbg_info_t::ATOM_S)
      font_colour = "#999900";
   if (addition_element_mode == lbg_info_t::ATOM_F)
      font_colour = "#00aa00";
   if (addition_element_mode == lbg_info_t::ATOM_CL)
      font_colour = "#229900";
   if (addition_element_mode == lbg_info_t::ATOM_I)
      font_colour = "#220077";
   return font_colour;
}


std::string
lbg_info_t::font_colour(const std::string &ele) const {
   std::string font_colour = "#111111";

   if (ele == "N")
      font_colour = "blue";
   if (ele == "O") 
      font_colour = "red";
   if (ele == "S") 
      font_colour = "#888800";
   if (ele == "F") 
      font_colour = "#006600";
   if (ele == "CL") 
      font_colour = "#116600";
   if (ele == "I") 
      font_colour = "#220066";
   
   return font_colour;
}

void
lbg_info_t::add_bond_to_atom_with_0_neighbours(int atom_index, int canvas_addition_mode) {

   // certain change
   
   widgeted_atom_t atom = mol.atoms[atom_index];
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   lig_build::pos_t current_atom_pos = atom.atom_position;
   lig_build::pos_t a_bond(SINGLE_BOND_CANVAS_LENGTH, 0.0);
   a_bond.rotate(60);
   lig_build::pos_t new_atom_pos = current_atom_pos + a_bond;
      
   widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
   int new_index = mol.add_atom(new_atom).second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   widgeted_bond_t b(atom_index, new_index, atom, new_atom, bt, root);
   mol.add_bond(b);

   change_atom_id_maybe(atom_index);
   change_atom_id_maybe(new_index);


   penultimate_atom_pos = atom.atom_position;
   penultimate_atom_index = atom_index;
   ultimate_atom_index    = new_index;
   
}


void
lbg_info_t::add_bond_to_atom_with_1_neighbour(int atom_index, int canvas_addition_mode,
					      int bond_index) {
   
   int index_1 = mol.bonds[bond_index].get_atom_1_index();
   int index_2 = mol.bonds[bond_index].get_atom_2_index();

   int other_atom_index = index_1;
   if (index_1 == atom_index)
      other_atom_index = index_2;
   
   widgeted_atom_t atom = mol.atoms[atom_index];

   lig_build::pos_t pos_1 = mol.atoms[atom_index].atom_position;
   lig_build::pos_t pos_2 = mol.atoms[other_atom_index].atom_position;

   lig_build::pos_t diff = pos_1 - pos_2;
   lig_build::pos_t du = diff.unit_vector();
   lig_build::pos_t candidate_new_vec_1 = du.rotate( 60);
   lig_build::pos_t candidate_new_vec_2 = du.rotate(-60);

   lig_build::pos_t candidate_pos_1 = pos_1 + candidate_new_vec_1 * SINGLE_BOND_CANVAS_LENGTH;
   lig_build::pos_t candidate_pos_2 = pos_1 + candidate_new_vec_2 * SINGLE_BOND_CANVAS_LENGTH;

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   std::vector<int> avoid_atoms(2);
   avoid_atoms[0] = atom_index;
   avoid_atoms[1] = other_atom_index;

   std::pair<bool,double> d1 = mol.dist_to_other_atoms_except(avoid_atoms, candidate_pos_1);
   std::pair<bool,double> d2 = mol.dist_to_other_atoms_except(avoid_atoms, candidate_pos_2);

   lig_build::pos_t new_atom_pos = candidate_pos_1;
   if (d1.first && d2.first) {
      if (d2.second > d1.second)
	 new_atom_pos = candidate_pos_2;
   }

   // Make an atom, add it, make a bond (using that atom and atom_index) and add it.


   // GooCanvasItem *ci = canvas_line_bond(pos_1, new_atom_pos, root, canvas_addition_mode);
   
   widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
   int new_index = mol.add_atom(new_atom).second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   widgeted_bond_t b(atom_index, new_index, atom, new_atom, bt, root);
   mol.add_bond(b);

   change_atom_id_maybe(atom_index);
   change_atom_id_maybe(new_index);

   penultimate_atom_pos = atom.atom_position;
   penultimate_atom_index = atom_index;
   ultimate_atom_index    = new_index;
   
}

void
lbg_info_t::add_bond_to_atom_with_2_neighbours(int atom_index, int canvas_addition_mode,
					       const std::vector<int> &bond_indices) {

   widgeted_atom_t atom = mol.atoms[atom_index];
   int atom_index_1 = mol.bonds[bond_indices[0]].get_atom_1_index();
   int atom_index_2 = mol.bonds[bond_indices[0]].get_atom_2_index();
   int atom_index_3 = mol.bonds[bond_indices[1]].get_atom_1_index();
   int atom_index_4 = mol.bonds[bond_indices[1]].get_atom_2_index();

   int A_index = atom_index_1;
   if (atom_index_1 == atom_index)
      A_index = atom_index_2;
   int B_index = atom_index_3;
   if (atom_index_3 == atom_index)
      B_index = atom_index_4;

   lig_build::pos_t A_neighb = mol.atoms[A_index].atom_position;
   lig_build::pos_t B_neighb = mol.atoms[B_index].atom_position;
   lig_build::pos_t atom_pos = mol.atoms[atom_index].atom_position;

   lig_build::pos_t mp =
      lig_build::pos_t::mid_point(A_neighb, B_neighb);
   lig_build::pos_t mpa = atom_pos - mp;
   lig_build::pos_t mpa_unit = mpa.unit_vector();

   lig_build::pos_t new_atom_pos = atom_pos + mpa_unit * SINGLE_BOND_CANVAS_LENGTH;
   
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
   widgeted_atom_t at(new_atom_pos, "C", 0, NULL);
   int new_index = mol.add_atom(at).second;
   
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   widgeted_bond_t b(atom_index, new_index, atom, at, bt, root);
   mol.add_bond(b);


   // Now, what about the atom that has been added to?  Now that we
   // have a extra bond, that atoms name could go from NH -> N (or C
   // -> C+ for 5 bonds).
   //
   change_atom_id_maybe(atom_index); 
   change_atom_id_maybe(new_index);
   penultimate_atom_pos = atom.atom_position;
   penultimate_atom_index = atom_index;
   ultimate_atom_index    = new_index;
}


void
lbg_info_t::add_bond_to_atom_with_3_neighbours(int atom_index, int canvas_addition_mode,
					       const std::vector<int> &bond_indices) {

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   // are (at least) 2 of the bonds attached to atom_index terminated
   // at their end?
   //
   // If so, we shall remove and replace those bonds before we add
   // this new one.
   //
   std::pair<bool, std::vector<int> > pr = have_2_stubs_attached_to_atom(atom_index, bond_indices);
   std::vector<int> attached_bonds = pr.second;
   
   if (pr.first) {
      int l = attached_bonds.size();
      
      widgeted_bond_t bond_to_core = orthogonalise_2_bonds(atom_index, attached_bonds, bond_indices);
      // now add a third

      int p1_idx = bond_to_core.get_other_index(atom_index);
      lig_build::pos_t p1 = mol.atoms[p1_idx].atom_position;
      
      lig_build::pos_t central_atom_pos = mol.atoms[atom_index].atom_position;
      lig_build::pos_t existing_bond_dir = central_atom_pos - p1;
      lig_build::pos_t ebd_uv = existing_bond_dir.unit_vector();

      lig_build::pos_t new_atom_pos = central_atom_pos + ebd_uv * SINGLE_BOND_CANVAS_LENGTH;
      widgeted_atom_t atom(new_atom_pos, "C", 0, NULL);
      int new_atom_index = mol.add_atom(atom).second;
      lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
      widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], atom, bt, root);
      mol.add_bond(b);
      change_atom_id_maybe(atom_index);
      change_atom_id_maybe(new_atom_index);

   } else {
      // bond_indices are bonds have an atom that is atom_index.
      squeeze_in_a_4th_bond(atom_index, canvas_addition_mode, bond_indices);
   } 

   
}

void
lbg_info_t::squeeze_in_a_4th_bond(int atom_index, int canvas_addition_mode,
				  const std::vector<int> &bond_indices) {

   
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
   std::vector<double> angles = get_angles(atom_index, bond_indices);
   std::vector<double> sorted_angles = angles;
   std::sort(sorted_angles.begin(), sorted_angles.end());
   if (sorted_angles[0] > 115.0) { // smallest angle, i.e. each perfect 120.
      if (all_closed_rings(atom_index, bond_indices)) { 
	 lig_build::pos_t centre = mol.atoms[atom_index].atom_position;
	 bool found_from_pos = 0;
	 lig_build::pos_t from_pos(0,0); // unfilled
	 for (unsigned int i=0; i<bond_indices.size(); i++) { 
	    int idx = mol.bonds[bond_indices[i]].get_other_index(atom_index);
	    lig_build::pos_t pos = mol.atoms[idx].atom_position;
	    lig_build::pos_t diff = pos - centre;
	    double ori = diff.axis_orientation();
	    if (ori > -20) { 
	       if (ori < 90) { 
		  from_pos = pos;
		  found_from_pos = 1;
	       }
	    }
	 }
	 if (! found_from_pos) {
	    // use the first one.
	    int idx = mol.bonds[bond_indices[0]].get_other_index(atom_index);
	    from_pos = mol.atoms[idx].atom_position;
	 }

	 // now rotate from_pos by 60 degrees
	 lig_build::pos_t diff = from_pos - centre;
	 lig_build::pos_t d_uv = diff.unit_vector();
	 lig_build::pos_t d_uv_60 = d_uv.rotate(-60); // upside down canvas
	 lig_build::pos_t new_atom_pos = centre + d_uv_60 * SINGLE_BOND_CANVAS_LENGTH * 0.7;
	 widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
	 int new_atom_index = mol.add_atom(new_atom).second;
	 lig_build::bond_t::bond_type_t bt =
	    addition_mode_to_bond_type(canvas_addition_mode);
	 widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index],
			   new_atom, bt, root);
	 mol.add_bond(b);
	 change_atom_id_maybe(atom_index);
	 change_atom_id_maybe(new_atom_index);
	 ultimate_atom_index    = new_atom_index;

      } else {
	 
	 // OK, there were not centres on all side of all the bonds.
	 // Go place the new atom not towards a ring centre

	 std::cout << "Go place a bond not towards a ring centre" << std::endl;
	 lig_build::pos_t new_atom_pos = get_new_pos_not_towards_ring_centres(atom_index, bond_indices);
	 widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
	 int new_atom_index = mol.add_atom(new_atom).second;
	 lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
	 widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], new_atom, bt, root);
	 mol.add_bond(b);
	 change_atom_id_maybe(atom_index);
	 change_atom_id_maybe(new_atom_index);
	 ultimate_atom_index    = new_atom_index;
	 
      }
      

   } else {
      lig_build::pos_t new_atom_pos = new_pos_by_bisection(atom_index, bond_indices, angles, root);
      widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
      int new_atom_index = mol.add_atom(new_atom).second;
      lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
      widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], new_atom, bt, root);
      mol.add_bond(b);
      change_atom_id_maybe(atom_index);
      change_atom_id_maybe(new_atom_index);
      ultimate_atom_index    = new_atom_index;
   } 

   widgeted_atom_t atom = mol.atoms[atom_index];
   penultimate_atom_pos = atom.atom_position;
   penultimate_atom_index = atom_index;

}


// bond_indices is a vector of bond indices that have an atom with
// index atom_index.
//
bool
lbg_info_t::all_closed_rings(int atom_index, const std::vector<int> &bond_indices) const {

   bool status = 0;

   if (bond_indices.size() > 3) {
      std::vector<lig_build::pos_t> centres = get_centres_from_bond_indices(bond_indices);
      if (centres.size() > 2)
	 status = 1;
   }
   return status;
}

std::vector<lig_build::pos_t>
lbg_info_t::get_centres_from_bond_indices(const std::vector<int> &bond_indices) const {

   std::vector<lig_build::pos_t> centres;
   for (unsigned int ib=0; ib<bond_indices.size(); ib++) {
      if (mol.bonds[bond_indices[ib]].have_centre_pos()) {
	 lig_build::pos_t test_centre = mol.bonds[bond_indices[ib]].centre_pos();

	 // add test_centre to centres only if it wasnt there already.
	 bool found_centre = 0;
	 for (unsigned int j=0; j<centres.size(); j++) {
	    if (test_centre.close_point(centres[j])) {
	       found_centre = 1; // it was already there
	       break;
	    }
	 }
	 if (! found_centre)
	    centres.push_back(test_centre);
      }
   }
   return centres;
} 

lig_build::pos_t
lbg_info_t::get_new_pos_not_towards_ring_centres(int atom_index,
						 const std::vector<int> &bond_indices) const { 

   lig_build::pos_t centre = mol.atoms[atom_index].atom_position;
   lig_build::pos_t p;
   std::vector<lig_build::pos_t> centres = get_centres_from_bond_indices(bond_indices);
   if (centres.size() == 2) { 
      lig_build::pos_t diff = centres[1] - centres[0];
      lig_build::pos_t diff_uv_90 = diff.unit_vector().rotate(90);
      lig_build::pos_t extra = diff_uv_90 * SINGLE_BOND_CANVAS_LENGTH;
      lig_build::pos_t candidate_1 = centre + extra;
      lig_build::pos_t candidate_2 = centre - extra;
      lig_build::pos_t best_candidate;
      double d_1 = 10000;
      double d_2 = 10000;
      for (unsigned int i=0; i<mol.atoms.size(); i++) { 
	 double dt_1 = lig_build::pos_t::length(candidate_1, mol.atoms[i].atom_position);
	 double dt_2 = lig_build::pos_t::length(candidate_2, mol.atoms[i].atom_position);
	 if (dt_1 < d_1)
	    d_1 = dt_1;
	 if (dt_2 < d_2)
	    d_2 = dt_2;
      }
      if (d_1 < 9999) { 
	 p = candidate_1;
	 if (d_2 > d_1)
	    p = candidate_2;
      }
   } else {
      // build away from ring_centre
      for (unsigned int i=0; i<bond_indices.size(); i++) { 
	 if (mol.bonds[bond_indices[i]].have_centre_pos()) {
	    lig_build::pos_t ring_centre = mol.bonds[bond_indices[i]].centre_pos();
	    int other_index = mol.bonds[bond_indices[i]].get_other_index(atom_index);
	    lig_build::pos_t p1 = mol.atoms[other_index].atom_position;
	    lig_build::pos_t p2 = mol.atoms[atom_index ].atom_position;
	    lig_build::pos_t bond_dir = p2 - p1;
	    lig_build::pos_t bond_dir_uv = bond_dir.unit_vector();
	    p = centre + bond_dir_uv * SINGLE_BOND_CANVAS_LENGTH * 0.8;
	    break;
	 }
      }
   }

   return p;
}


lig_build::pos_t
lbg_info_t::new_pos_by_bisection(int atom_index,
				 const std::vector<int> &bond_indices,
				 const std::vector<double> &angles,
				 GooCanvasItem *root) const {

   lig_build::pos_t p(0,0);
   std::vector<double> sorted_angles = angles;
   std::sort(sorted_angles.begin(), sorted_angles.end());
   lig_build::pos_t centre = mol.atoms[atom_index].atom_position;
   for (unsigned int i=0; i<bond_indices.size(); i++) { 
      int j = i+1;
      if (j == bond_indices.size())
	 j = 0;
      int idx_1 = mol.bonds[bond_indices[i]].get_other_index(atom_index);
      int idx_2 = mol.bonds[bond_indices[j]].get_other_index(atom_index);
      lig_build::pos_t pos_1 = mol.atoms[idx_1].atom_position;
      lig_build::pos_t pos_2 = mol.atoms[idx_2].atom_position;
      lig_build::pos_t d_1 = pos_1 - centre;
      lig_build::pos_t d_2 = pos_2 - centre;
      double dot = lig_build::pos_t::dot(d_1, d_2);
      double cos_theta = dot/(d_1.length() * d_2.length());
      double angle = acos(cos_theta)/DEG_TO_RAD;
      if (angle > (sorted_angles.back()-0.01)) {
	 lig_build::pos_t mp_diff = lig_build::pos_t::mid_point(pos_1, pos_2) - centre;
	 lig_build::pos_t mp_diff_uv = mp_diff.unit_vector();
	 p = centre + mp_diff_uv * SINGLE_BOND_CANVAS_LENGTH * 0.9;
	 break;
      }
   }

   return p;
}

std::vector<double>
lbg_info_t::get_angles(int atom_index, const std::vector<int> &bond_indices) const {

   std::vector<double> v(bond_indices.size());

   lig_build::pos_t centre = mol.atoms[atom_index].atom_position;

   for (unsigned int i=0; i<bond_indices.size(); i++) {
      int j = i+1;
      if (j == bond_indices.size())
	 j = 0;
      int idx_1 = mol.bonds[bond_indices[i]].get_other_index(atom_index);
      int idx_2 = mol.bonds[bond_indices[j]].get_other_index(atom_index);
      lig_build::pos_t pos_1 = mol.atoms[idx_1].atom_position;
      lig_build::pos_t pos_2 = mol.atoms[idx_2].atom_position;
      lig_build::pos_t d_1 = pos_1 - centre;
      lig_build::pos_t d_2 = pos_2 - centre;
      double dot = lig_build::pos_t::dot(d_1, d_2);
      double cos_theta = dot/(d_1.length() * d_2.length());
      double angle = acos(cos_theta)/DEG_TO_RAD;
      v[i] = angle;
   }
   return v;
}


// Delete the last 2 atoms of stub_attached_atoms, regenerate new
// atoms.
//
// bond_indices are the indices of the bonds attached to atom_index
// (the "central" atom).
// 
widgeted_bond_t
lbg_info_t::orthogonalise_2_bonds(int atom_index,
				  const std::vector<int> &stub_attached_atoms,
				  const std::vector<int> &bond_indices) { 

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   // this is the set of bonds that are just stubs to a C or so
   // (branch is no more than 1 atom).
   //
   int l = stub_attached_atoms.size();

   if (l < 2) {
      std::cout << "ERROR:: trapped l = " << l << " in orthogonalise_2_bonds()"
		<< std::endl;
      widgeted_bond_t dummy;
      return dummy;
   }
   
   int ind_1 = stub_attached_atoms[l-1];
   int ind_2 = stub_attached_atoms[l-2];
   std::string ele_1 = mol.atoms[ind_1].element;
   std::string ele_2 = mol.atoms[ind_2].element;
   mol.close_atom(ind_1, root);
   mol.close_atom(ind_2, root);

   lig_build::pos_t central_atom_pos = mol.atoms[atom_index].atom_position;

   // Find the bond that does not have an atom in stub_attached_atoms.
   bool found_bond = 0;
   widgeted_bond_t bond_to_core;
   for (unsigned int i=0; i<bond_indices.size(); i++) { 
      int idx_1 = mol.bonds[bond_indices[i]].get_atom_1_index();
      int idx_2 = mol.bonds[bond_indices[i]].get_atom_2_index();
      if (std::find(stub_attached_atoms.begin(), stub_attached_atoms.end(), idx_1) != stub_attached_atoms.end()) {
	 if (std::find(stub_attached_atoms.begin(), stub_attached_atoms.end(), idx_2) != stub_attached_atoms.end()) {
	    found_bond = 1;
	    bond_to_core = mol.bonds[bond_indices[i]];
	 }
      }
   }

   if (! found_bond) {

      // OK, pathological case of all of the atoms attached to the
      // central atom being stub atoms?
      //
      bond_to_core = mol.bonds[bond_indices[0]];
      found_bond = 1;
   }

   // p1 is the position of the atom that is bonded to the central
   // atom, i.e. is towards the core of the ligand.  The vector from
   // p1 to the central atom is the vector from which all the others
   // are based.
   //
   int p1_idx = bond_to_core.get_other_index(atom_index);
   lig_build::pos_t p1 = mol.atoms[p1_idx].atom_position;
   lig_build::pos_t existing_bond_dir = central_atom_pos - p1;
   
   lig_build::pos_t ebd_uv = existing_bond_dir.unit_vector();
   lig_build::pos_t ebd_uv_90 = ebd_uv.rotate(90);

   lig_build::pos_t new_atom_pos_1 = central_atom_pos + ebd_uv_90 * SINGLE_BOND_CANVAS_LENGTH;
   lig_build::pos_t new_atom_pos_2 = central_atom_pos - ebd_uv_90 * SINGLE_BOND_CANVAS_LENGTH;
   widgeted_atom_t atom_1(new_atom_pos_1, ele_1, 0, NULL);
   widgeted_atom_t atom_2(new_atom_pos_2, ele_2, 0, NULL);
   int new_atom_index_1 = mol.add_atom(atom_1).second;
   int new_atom_index_2 = mol.add_atom(atom_2).second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   widgeted_bond_t b_1(atom_index, new_atom_index_1, mol.atoms[atom_index], atom_1, bt, root);
   widgeted_bond_t b_2(atom_index, new_atom_index_2, mol.atoms[atom_index], atom_2, bt, root);
   mol.add_bond(b_1);
   mol.add_bond(b_2);

   return bond_to_core; // so that we don't have to recalculate it in caller.

}



// return a status and a vector of atoms (bonded to atom_index) having
// only one bond.
// 
std::pair<bool, std::vector<int> > 
lbg_info_t::have_2_stubs_attached_to_atom(int atom_index, const std::vector<int> &bond_indices) const {

   std::vector<int> v;
   for (unsigned int i=0; i<bond_indices.size(); i++) { 
      int other_index = mol.bonds[bond_indices[i]].get_other_index(atom_index);
      // now, does other_index have only one bond?
      std::vector<int> local_bonds = mol.bonds_having_atom_with_atom_index(other_index);
      if (local_bonds.size() == 1) {
	 v.push_back(other_index);
      }
   }
   bool status = 0;
   if (v.size() > 1)
      status = 1;
   return std::pair<bool, std::vector<int> > (status, v);
}



// Set the addition mode and turn off other toggle buttons (if
// needed).
void
lbg_handle_toggle_button(GtkToggleToolButton *tb, GtkWidget *canvas, int mode) {

   gpointer gp = gtk_object_get_user_data(GTK_OBJECT(canvas));
   lbg_info_t *l = static_cast<lbg_info_t *> (gp);

   if (l) { 
      if (gtk_toggle_tool_button_get_active(tb)) {

	 l->untoggle_others_except(tb);
	 l->canvas_addition_mode = mode;
      } else {
	 l->canvas_addition_mode = lbg_info_t::NONE;
      }
   } else {
      std::cout << "ERROR:: in lbg_handle_toggle_button() failed to get lbg pointer!"
		<< std::endl;
   } 
}


bool
lbg_info_t::try_stamp_hexagon_aromatic(int x_cen, int y_cen, bool shift_is_pressed) {

   bool is_aromatic = 1;
   return try_stamp_polygon(6, x_cen, y_cen, shift_is_pressed, is_aromatic);
}

// shift_is_pressed implies spiro.  Return changed status.
bool
lbg_info_t::try_stamp_polygon(int n_edges, int x_cen, int y_cen, bool is_spiro, bool is_aromatic) {

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));
   bool changed_status = 0;

   if (mol.is_empty()) {
      stamp_polygon_anywhere(n_edges, x_cen, y_cen, is_aromatic, root);
      changed_status = 1; // stamp_polygon_anywhere() always modifies
   } else { 
      if (highlight_data.has_contents()) {
	 std::vector<int> new_atoms =
	    try_stamp_polygon_using_highlighted_data(n_edges, is_spiro, is_aromatic, root);
	 if (new_atoms.size() > 0)
	    changed_status = 1;
	 if (highlight_data.single_atom()) {
	    if (new_atoms.size() > 0) {
	       // we need to add the spur bond, if this is not a spiro compound

	       if (! is_spiro) {

		  // what are the atom indices of the spur bond?

		  int highlighted_atom_index  = highlight_data.get_atom_index();
		  if (highlighted_atom_index != new_atoms[0]) { 
		     // std::cout << "adding spur bond...!" << std::endl;
		     
		     if (highlighted_atom_index != UNASSIGNED_INDEX) {
			lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
			widgeted_bond_t bond(highlighted_atom_index, new_atoms[0],
					     mol.atoms[highlighted_atom_index],
					     mol.atoms[new_atoms[0]],
					     bt, root);
			mol.add_bond(bond);
		     }
		  }
	       }
	    }
	 }
      }
   }
   return changed_status;
}

// always modifies.  Always sets penultimate_atom_pos.
// 
void
lbg_info_t::try_stamp_bond_anywhere(int canvas_addition_mode, int x_mouse, int y_mouse) {

   bool changed_status = 0;

   double dx = double(x_mouse);
   double dy = double(y_mouse);
   lig_build::pos_t new_atom_pos(dx, dy);

    widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
    int new_atom_index = mol.add_atom(new_atom).second;
    add_bond_to_atom(new_atom_index, canvas_addition_mode);
    penultimate_atom_pos = new_atom_pos;
    penultimate_atom_index = new_atom_index;
}


double
lbg_info_t::radius(int n_edges) const {

   double angle_step = 360.0/double(n_edges);
   double r = SINGLE_BOND_CANVAS_LENGTH/(2*sin(angle_step*DEG_TO_RAD*0.5));
   return r;
}

// always modifies
void
lbg_info_t::stamp_polygon_anywhere(int n_edges, int x_cen, int y_cen,
				   bool is_aromatic,
				   GooCanvasItem *root) {

   double angle_off = 0;
   lig_build::polygon_position_info_t ppi(x_cen, y_cen, angle_off);
   stamp_polygon(n_edges, ppi, is_aromatic, root);
}


// is_aromatic means (for now) that double and single bonds alternate.
// 
// always modifies
// 
std::vector<int>
lbg_info_t::stamp_polygon(int n_edges, lig_build::polygon_position_info_t ppi,
			  bool is_aromatic, GooCanvasItem *root) {

   // use the next two to construct the returned value.
   std::vector<int> atom_indices;
      
   double angle_step = 360.0/double(n_edges);
   double dx_cen = double(ppi.pos.x);
   double dy_cen = double(ppi.pos.y);
   lig_build::pos_t centre(dx_cen, dy_cen);
   std::vector<std::pair<lig_build::pos_t, std::pair<bool,int> > > atom_pos_index;

   // This makes the bond length the same, no matter how many edges.
   //
   double r = radius(n_edges);

   // offset so that the first bond lies along the X axis:
   // 
   // double internal_angle_offset = -angle_step/2.0;
   // double orientation_correction = 180;
   
   double internal_angle_offset = 0;
   double orientation_correction = 180;
   if (! ppi.apply_internal_angle_offset_flag) { 
      internal_angle_offset = 0;
      orientation_correction = 0;
   }
   
   for (int i=0; i<n_edges; i++) {
      double theta_deg =
	 (i * angle_step + orientation_correction + internal_angle_offset + ppi.angle_offset);
      double theta = theta_deg * DEG_TO_RAD;
      double pt_1_x = r * sin(theta) + dx_cen;
      double pt_1_y = r * cos(theta) + dy_cen;
      lig_build::pos_t pt(pt_1_x, pt_1_y);
      GooCanvasItem *aci = 0; // initially carbon, no canvas item
      widgeted_atom_t at(pt, "C", 0, aci);
      std::pair<bool,int> atom_index = mol.add_atom(at);
      std::pair<lig_build::pos_t, std::pair<bool,int> > pr(pt, atom_index);
      atom_pos_index.push_back(pr);
      if (0) { // debug polygon vertex placement
	 std::string stroke_colour = get_stroke_colour(i,atom_pos_index.size());
	 GooCanvasItem *rect_item  =
	    goo_canvas_rect_new (root,
				 pt.x - 5.0, 
				 pt.y - 5.0,
				 10.0, 10.0,
				 "line-width", 7.0,
				 "stroke-color", stroke_colour.c_str(),
				 NULL);
      }
   }
   
   for (unsigned int i=0; i<atom_pos_index.size(); i++) {
      int j = i + 1;
      if (j == atom_pos_index.size())
	 j = 0;
      int i_index = atom_pos_index[i].second.second;
      int j_index = atom_pos_index[j].second.second;


      // if either the first or the second atoms were new atoms then
      // draw the bond (and of course, of the first and second atom
      // indexes were both from pre-existing atoms then *don't* draw
      // the bond).
      // 
      if (atom_pos_index[i].second.first || atom_pos_index[j].second.first) { 
	 lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
	 if (is_aromatic) {
	    if (((i/2)*2) == i)
	       bt = lig_build::bond_t::DOUBLE_BOND;
	 }
	 widgeted_bond_t bond(i_index, j_index,
			      mol.atoms[atom_pos_index[i].second.second],
			      mol.atoms[atom_pos_index[j].second.second], centre, bt, root);
	 mol.add_bond(bond);

      }
   }

   // transfer from atom_pos_index to atom_indices
   atom_indices.resize(atom_pos_index.size());
   for (unsigned int i=0; i<atom_pos_index.size(); i++)
      atom_indices[i] = atom_pos_index[i].second.second;
   
   return atom_indices;
}


std::vector<int>
lbg_info_t::try_stamp_polygon_using_highlighted_data(int n_edges,
						     bool spiro_flag,
						     bool is_aromatic,
						     GooCanvasItem *root) {

   std::vector<int> new_atoms;
   
   double alpha = 2 * M_PI/double(n_edges);
   double corr_r = radius(n_edges) * cos(alpha/2.0);
   double rad_std = radius(n_edges);
   
   lig_build::polygon_position_info_t ppi =
      highlight_data.get_new_polygon_centre(n_edges, spiro_flag, rad_std, corr_r, mol);
   if (ppi.can_stamp) 
      new_atoms = stamp_polygon(n_edges, ppi, is_aromatic, root);
   return new_atoms;
}

lig_build::polygon_position_info_t
lbg_info_t::highlight_data_t::get_new_polygon_centre(int n_edges, bool spiro_flag,
						     const double &rad_std,
						     const double &cor_radius,
						     const widgeted_molecule_t &mol) const {

   lig_build::polygon_position_info_t ppi;

   if (n_atoms_ == 2)
      ppi = get_new_polygon_centre_using_2_atoms(n_edges, cor_radius);

   if (n_atoms_ == 1)
      ppi = get_new_polygon_centre_using_1_atom(n_edges, spiro_flag, rad_std, cor_radius, mol);

   return ppi;
}


lig_build::polygon_position_info_t
lbg_info_t::highlight_data_t::get_new_polygon_centre_using_2_atoms(int n_edges,
								   const double &cor_radius) const {

   // we pass the cor_radius, how much to move the centre of the new
   // polygon from the mid point of "this" bond (which is not the
   // radius because the bond mid-point is "inside" the circle for the
   // new polygon points.

   double x_cen = 0;
   double y_cen = 0;

   // Make a bisector (vector) of the bond.  The new polygon centre
   // should lie on that bisector.
   //
   // The distance along the bisector depends on the number of edges.
   //
   // Bond is A->B, M is the midpoint.  M->P is the bisector.  C, the
   // centre of the new polygon lies along M->P, at r from M. Let's
   // say that M->P is a unit vector, d is |A->M|.
   //
   //
   //              A
   //              |\                                   . 
   //              | \                                  .
   //              |  \                                 .
   //              |   \                                .
   //            d |    \                               .
   //              |     \                              .
   //              |      \                             .
   //            M |-------\ P ................... C
   //              |   1 
   //              |
   //              |
   //              |
   //              |
   //              B
   //
   // So the question here is: where is P?
   // 
   //    Take a unit vector along A->B.  Rotate it 90 degress about M.
   //
   //

   lig_build::pos_t a = pos_1_;
   lig_build::pos_t b = pos_2_;
   lig_build::pos_t m = lig_build::pos_t::mid_point(a,b);
   lig_build::pos_t ab_unit = (b - a).unit_vector();
   lig_build::pos_t mp_unit = ab_unit.rotate(90);
   lig_build::pos_t c_1 = m + mp_unit * cor_radius;
   lig_build::pos_t c_2 = m + mp_unit * cor_radius;

   lig_build::pos_t c = c_1;
   if (has_ring_centre_flag) {
      double l1 = lig_build::pos_t::length(c_1, ring_centre);
      double l2 = lig_build::pos_t::length(c_2, ring_centre);
      if (l1 < l2)
	 c = c_2;
   }

   //  0    for a triangle : 120
   //  45   for a square   : 90
   //  0    for pentagon   : 72
   // 30    for hexagon    : 60
   // 51.43 for pentagon   : 51.43
   // 22.5 for octagon : 45

   // orientation correction
   // 
   // double oc = 0;

   lig_build::pos_t ab = b - a;

   // canvas is upside down:
   lig_build::pos_t ab_cor(ab.x, -ab.y);
   
   double ao =  ab_cor.axis_orientation();

   // std::cout << "axis orientation of " << ab_cor << " is " << ao << std::endl;

   double angle_step = 360.0/double(n_edges);
   double internal_angle_offset = -angle_step/2.0;

   double angle_off = ao + internal_angle_offset;
   
   lig_build::polygon_position_info_t ppi(c, angle_off);
   return ppi;
   
}

// This creates a polygon on a spur (unless spiro_flag is set)
// 
lig_build::polygon_position_info_t
lbg_info_t::highlight_data_t::get_new_polygon_centre_using_1_atom(int n_edges,
								  bool spiro_flag,
								  const double &radius,
								  const double &cor_radius,
								  const widgeted_molecule_t &mol) const {

   lig_build::polygon_position_info_t ppi(pos_1_, 0);
   lig_build::pos_t A = pos_1_;

   int atom_index = get_atom_index();
   std::vector<int> bv = mol.bonds_having_atom_with_atom_index(atom_index);

   // std::cout << "in get_new_polygon_centre_using_1_atom() bv size is " << bv.size()
   // << std::endl;

   if (bv.size() == 2) {
      std::vector<lig_build::pos_t> neighbours;
      lig_build::pos_t test_pos = mol.atoms[mol.bonds[bv[0]].get_atom_1_index()].atom_position;
      if (! test_pos.near_point(A, 2))
	 neighbours.push_back(test_pos);
      test_pos = mol.atoms[mol.bonds[bv[0]].get_atom_2_index()].atom_position;
      if (! test_pos.near_point(A, 2))
	 neighbours.push_back(test_pos);
      test_pos = mol.atoms[mol.bonds[bv[1]].get_atom_1_index()].atom_position;
      if (! test_pos.near_point(A, 2))
	 neighbours.push_back(test_pos);
      test_pos = mol.atoms[mol.bonds[bv[1]].get_atom_2_index()].atom_position;
      if (! test_pos.near_point(A, 2))
	 neighbours.push_back(test_pos);

      if (neighbours.size() == 2) {
	 lig_build::pos_t mp =
	    lig_build::pos_t::mid_point(neighbours[0], neighbours[1]);
	 lig_build::pos_t mpa = A - mp;
	 lig_build::pos_t mpa_unit = mpa.unit_vector();

	 lig_build::pos_t new_centre = A + mpa_unit * radius;
	 if (! spiro_flag) {
	    new_centre += mpa_unit * SINGLE_BOND_CANVAS_LENGTH;
	 } 

	 ppi.pos = new_centre;
	 lig_build::pos_t bond_vector = neighbours[1] - neighbours[0];

	 // upside-down canvas
	 double angle = -mpa.axis_orientation();
	 // std::cout << "spur orientation " << angle << std::endl;
	 
	 ppi.angle_offset = angle - 90;
	 ppi.apply_internal_angle_offset_flag = 0; // don't internally correct orientation
      }
   }

   // now, say we are adding a ring to the end of a chain of single bonds
   // 
   if (bv.size() == 1) {
      // We need the coordinates of the bond atom, We have one of
      // those (A), what is the other?
      lig_build::pos_t pos_other = mol.atoms[mol.bonds[bv[0]].get_atom_1_index()].atom_position;
      if (A.close_point(pos_other))
	 pos_other = mol.atoms[mol.bonds[bv[0]].get_atom_2_index()].atom_position;
      lig_build::pos_t other_to_A = (A - pos_other);
      lig_build::pos_t other_to_A_unit = other_to_A.unit_vector();
      lig_build::pos_t new_centre = A + other_to_A_unit * radius;
      double angle = -other_to_A_unit.axis_orientation();
      ppi.pos = new_centre;
      ppi.angle_offset = angle - 90;
      ppi.apply_internal_angle_offset_flag = 0; // don't internally correct orientation
   }
   return ppi;
}


#if ( ( (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 11) ) || GTK_MAJOR_VERSION > 2)
bool
lbg_info_t::save_togglebutton_widgets(GtkBuilder *builder) {

   std::vector<std::string> w_names;
   // w_names.push_back("charge_toggle_toolbutton");

   w_names.push_back("single_toggle_toolbutton");
   w_names.push_back("double_toggle_toolbutton");
   w_names.push_back("triple_toggle_toolbutton");
   w_names.push_back("stereo_out_toggle_toolbutton");
   w_names.push_back("c3_toggle_toolbutton");
   w_names.push_back("c4_toggle_toolbutton");
   w_names.push_back("c5_toggle_toolbutton");
   w_names.push_back("c6_toggle_toolbutton");
   w_names.push_back("c6_arom_toggle_toolbutton");
   w_names.push_back("c7_toggle_toolbutton");
   w_names.push_back("c8_toggle_toolbutton");
   w_names.push_back("carbon_toggle_toolbutton");
   w_names.push_back("nitrogen_toggle_toolbutton");
   w_names.push_back("oxygen_toggle_toolbutton");
   w_names.push_back("sulfur_toggle_toolbutton");
   w_names.push_back("phos_toggle_toolbutton");
   w_names.push_back("fluorine_toggle_toolbutton");
   w_names.push_back("chlorine_toggle_toolbutton");
   w_names.push_back("bromine_toggle_toolbutton");
   w_names.push_back("iodine_toggle_toolbutton");
   w_names.push_back("other_element_toggle_toolbutton");
   w_names.push_back("delete_item_toggle_toolbutton");

   for (unsigned int i=0; i<w_names.size(); i++) {
      GtkToggleToolButton *tb =
	 GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object (builder, w_names[i].c_str()));
      widget_names[w_names[i]] = tb;
   }
   return TRUE;
} 
#endif // GTK_VERSION


// perhaps we shouldn't use this function for both the user-click
// clear and the canvas-drag clear...  Hmm... (because in the latter
// case, we don't need to regnerate the SMILES string and the QED
// score.
// 
void
lbg_info_t::clear() {

   clear_canvas();
   // and that clears the alerts group
   alert_group = NULL;

   key_group = NULL;
   
   // clear the molecule
   mol.clear();

   update_descriptor_attributes();
}

void
lbg_info_t::clear_canvas() {

   // clear the canvas
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

   gint n_children = goo_canvas_item_get_n_children (root);
   for (int i=0; i<n_children; i++)
      goo_canvas_item_remove_child(root, 0);
}




#if ( ( (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 11) ) || GTK_MAJOR_VERSION > 2)
bool
lbg_info_t::init(GtkBuilder *builder) {

   GtkWidget *pe_test_function_button = NULL;

   if (use_graphics_interface_flag) {
      
      lbg_window = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_window"));
      if (! lbg_window) {

	 lbg_window = NULL;
	 about_dialog = NULL; 
	 lbg_search_combobox = NULL;
	 open_dialog = NULL;
	 save_as_dialog = NULL;
	 lbg_export_as_pdf_dialog = NULL;
	 lbg_export_as_png_dialog = NULL;
	 lbg_sbase_search_results_dialog = NULL;
	 lbg_sbase_search_results_vbox = NULL;
	 lbg_smiles_dialog = NULL;
	 lbg_smiles_entry = NULL;
	 lbg_statusbar = NULL;
	 lbg_toolbar_layout_info_label = NULL;
	 lbg_atom_x_dialog = NULL;
	 lbg_atom_x_entry = NULL;
	 canvas = NULL;
	 return false; // boo.

      } else {

	 gtk_widget_show (lbg_window);
	 about_dialog =    GTK_WIDGET (gtk_builder_get_object (builder, "lbg_aboutdialog"));
	 open_dialog     = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_open_filechooserdialog"));
	 save_as_dialog  = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_save_as_filechooserdialog"));
	 lbg_sbase_search_results_dialog = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_sbase_search_results_dialog"));
	 lbg_sbase_search_results_vbox = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_sbase_search_results_vbox"));
	 lbg_export_as_pdf_dialog =      GTK_WIDGET (gtk_builder_get_object (builder, "lbg_export_as_pdf_filechooserdialog"));
	 lbg_export_as_png_dialog =      GTK_WIDGET (gtk_builder_get_object (builder, "lbg_export_as_png_filechooserdialog"));
	 lbg_smiles_dialog =             GTK_WIDGET(gtk_builder_get_object(builder, "lbg_smiles_dialog"));
	 lbg_smiles_entry =              GTK_WIDGET(gtk_builder_get_object(builder, "lbg_smiles_entry"));
	 lbg_search_combobox =           GTK_WIDGET(gtk_builder_get_object(builder, "lbg_search_combobox"));
	 lbg_statusbar =                 GTK_WIDGET(gtk_builder_get_object(builder, "lbg_statusbar"));
	 lbg_toolbar_layout_info_label = GTK_WIDGET(gtk_builder_get_object(builder, "lbg_toolbar_layout_info_label"));
	 lbg_atom_x_dialog =             GTK_WIDGET(gtk_builder_get_object(builder, "lbg_atom_x_dialog"));
	 lbg_atom_x_entry =              GTK_WIDGET(gtk_builder_get_object(builder, "lbg_atom_x_entry"));
	 lbg_qed_hbox =                  GTK_WIDGET(gtk_builder_get_object(builder, "lbg_qed_hbox"));
	 lbg_qed_progressbar =           GTK_WIDGET(gtk_builder_get_object(builder, "lbg_qed_progressbar"));
	 lbg_qed_text_label =            GTK_WIDGET(gtk_builder_get_object(builder, "lbg_qed_text_label"));
	 lbg_alert_hbox =                GTK_WIDGET(gtk_builder_get_object(builder, "lbg_alert_hbox"));
	 lbg_alert_hbox_outer =          GTK_WIDGET(gtk_builder_get_object(builder, "lbg_alert_hbox_outer"));
	 lbg_show_alerts_checkbutton =   GTK_WIDGET(gtk_builder_get_object(builder, "lbg_show_alerts_checkbutton"));
	 lbg_alert_name_label =          GTK_WIDGET(gtk_builder_get_object(builder, "lbg_alert_name_label"));
	 lbg_get_drug_dialog =           GTK_WIDGET(gtk_builder_get_object(builder, "lbg_get_drug_dialog"));
	 lbg_get_drug_entry =            GTK_WIDGET(gtk_builder_get_object(builder, "lbg_get_drug_entry"));
	 pe_test_function_button =       GTK_WIDGET(gtk_builder_get_object(builder, "pe_test_function_button"));
	 lbg_flip_rotate_hbox =          GTK_WIDGET(gtk_builder_get_object(builder, "lbg_flip_rotate_hbox"));

	 gtk_label_set_text(GTK_LABEL(lbg_toolbar_layout_info_label), "---");
      }
   } 

   canvas = goo_canvas_new();

   if (0) {  // hide rdkit stuff
      GtkWidget *ww = GTK_WIDGET(gtk_builder_get_object(builder, "lbg_qed_and_alert_hbox"));
      gtk_widget_hide(lbg_flip_rotate_hbox);
      gtk_widget_hide(lbg_alert_hbox_outer);
      gtk_widget_hide(ww);
   }

   if (use_graphics_interface_flag) { 
      gtk_widget_set(GTK_WIDGET(canvas), "bounds-padding", 50.0, NULL);
      gtk_object_set_user_data(GTK_OBJECT(canvas), (gpointer) this);
      if (0) 
	 std::cout << ":::::: attached this lbg_info_t pointer to canvas: " << this
		   << " to " << canvas << std::endl;
   
      save_togglebutton_widgets(builder);
      GtkWidget *lbg_scrolled_win =
	 GTK_WIDGET(gtk_builder_get_object (builder, "lbg_scrolledwindow"));
      gtk_container_add(GTK_CONTAINER(lbg_scrolled_win), canvas);
      goo_canvas_set_bounds (GOO_CANVAS (canvas), 0, 0, 1200, 700);
   }

   GooCanvas *gc = GOO_CANVAS(canvas);

   if (0) 
      std::cout << "   after setting bounds: canvas adjustments: h-lower " 
		<< gc->hadjustment->lower << " h-upper-page "
		<< gc->hadjustment->upper - gc->hadjustment->page_size  << " v-lower "
		<< gc->vadjustment->lower  << " v-upper-page"
		<< gc->vadjustment->upper - gc->vadjustment->page_size  << " h-page "
		<< gc->hadjustment->page_size  << " "
		<< gc->vadjustment->page_size << std::endl;

   if (use_graphics_interface_flag) { 
   
      g_object_set_data_full(G_OBJECT(goo_canvas_get_root_item(GOO_CANVAS(canvas))),
			     "lbg-info", this, NULL);

      g_signal_connect(G_OBJECT(goo_canvas_get_root_item(GOO_CANVAS(canvas))),
		       "button_press_event",
		       G_CALLBACK(on_canvas_button_press), NULL);

      g_signal_connect(G_OBJECT(goo_canvas_get_root_item(GOO_CANVAS(canvas))),
		       "button_release_event",
		       G_CALLBACK(on_canvas_button_release), NULL);

      g_signal_connect(G_OBJECT(goo_canvas_get_root_item(GOO_CANVAS(canvas))),
		       "motion_notify_event",
		       G_CALLBACK(on_canvas_motion_new), NULL);
      
      //    g_signal_connect(G_OBJECT(goo_canvas_get_root_item(GOO_CANVAS(canvas))),
      // 		    "key_press_event",
      // 		    G_CALLBACK(on_key_press_event), NULL);

      // setup_lbg_drag_and_drop(lbg_window); // this lbg_info_t is not attached to this object.
      setup_lbg_drag_and_drop(canvas);
      gtk_widget_show(canvas);
      add_search_combobox_text();
   }

   // ------------ sbase ---------------------------
   // 
   init_sbase(".");

   
   // ------------ (old-style) watch for new files from coot ---------------------
   //
   if (is_stand_alone()) 
      int timeout_handle = gtk_timeout_add(500, watch_for_mdl_from_coot, this);

   // Hack in a button for PE to test stuff
   if (getenv("COOT_LBG_TEST_FUNCTION") != NULL)
      gtk_widget_show(pe_test_function_button);

   // if we don't have rdkit or python then we don't want to see qed progress bar
   // or the "show alerts" (because we can't match to the alert patterns).
   
#ifdef MAKE_ENTERPRISE_TOOLS
#ifdef USE_PYTHON   
   // all, with QED
#else   
   gtk_widget_hide(lbg_qed_hbox);
#endif    
#else
   gtk_widget_hide(lbg_qed_hbox);
   gtk_widget_hide(lbg_alert_hbox_outer);
   gtk_widget_hide(lbg_show_alerts_checkbutton); // perhaps this should be in the lbg_alert_hbox_outer?
#endif    

   return true;

}
#endif // GTK_VERSION

void
lbg_info_t::setup_lbg_drag_and_drop(GtkWidget *lbg_window) {

   // setup drag and drop
   int n_dnd_targets = 2;
   GtkTargetEntry target_list[] = {
      { (gchar *) "STRING",     0, TARGET_STRING },
      { (gchar *) "text/plain", 0, TARGET_STRING }}; // deprecated conversions from string constant
   // to char *.  GTK+ problem AFAICS.
   GtkDestDefaults dest_defaults = GTK_DEST_DEFAULT_ALL; // we don't need ALL, but this
                                                         // removes int->GtkDestDefaults
                                                         // conversion problems with |.
   

   gtk_drag_dest_set(GTK_WIDGET(lbg_window), /* widget that will accept a drop */
		     dest_defaults,
		     target_list,            /* lists of target to support */
		     n_dnd_targets,          
		     GDK_ACTION_COPY);       /* what to do with data after dropped */
   gtk_drag_dest_add_uri_targets(GTK_WIDGET(lbg_window));
   // if something was dropped
   g_signal_connect (GTK_WIDGET(lbg_window), "drag-drop",
		     G_CALLBACK (on_lbg_drag_drop), NULL);
   // what to we do with dropped data...
   g_signal_connect (GTK_WIDGET(lbg_window), "drag-data-received",
		     G_CALLBACK(on_lbg_drag_data_received), NULL);
   
} 
      

std::string
lbg_info_t::get_smiles_string_from_mol() const {

#ifdef MAKE_ENTERPRISE_TOOLS
   return get_smiles_string_from_mol_rdkit();
#else
   return get_smiles_string_from_mol_openbabel();
#endif
} 

std::string
lbg_info_t::get_smiles_string_from_mol_openbabel() const {

   std::string s;
   mol.write_mdl_molfile(".lbg-tmp-smiles-mol");
   std::string system_string = "babel -i mol .lbg-tmp-smiles-mol -o smi .lbg-tmp-smiles.smi";
   int status = system(system_string.c_str());
   if (status == 0) {
      std::string file_text = file_to_string(".lbg-tmp-smiles.smi");
      std::string::size_type ihash = file_text.find("#");
      std::string smiles_text = file_text;
      if (ihash != std::string::npos) {
	 s = file_text.substr(0, ihash);
      }
   }
   return s;
}

#ifdef MAKE_ENTERPRISE_TOOLS

void
lbg_info_t::update_qed(const RDKit::RWMol &rdkm) {

#ifdef USE_PYTHON   
   if (rdkm.getNumAtoms() == 0) {
      // non-interesting case first
      gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(lbg_qed_progressbar), 0);
      gtk_label_set_text(GTK_LABEL(lbg_qed_text_label), "");
   } else {
      bool all_set = false;
      double qed = 0.0;
      if (silicos_it_qed_default_func) { 
	 qed = get_qed(silicos_it_qed_default_func, rdkm);
	 if (qed > 0)
	    all_set = true;
      }

      if (all_set) { 
	 std::string s = coot::util::float_to_string_using_dec_pl(qed, 3);
	 gtk_label_set_text(GTK_LABEL(lbg_qed_text_label), s.c_str());
	 gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(lbg_qed_progressbar), qed);
      } else {
	 gtk_label_set_text(GTK_LABEL(lbg_qed_text_label), "");
	 gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(lbg_qed_progressbar), 0);
      }
   }
#endif
}
#endif

#ifdef MAKE_ENTERPRISE_TOOLS
void
lbg_info_t::update_alerts(const RDKit::RWMol &rdkm) {

   if (rdkm.getNumAtoms() == 0) {
      gtk_widget_hide(lbg_alert_hbox);
   } else { 

      if (show_alerts_user_control) {
	 
	 std::vector<alert_info_t> v = alerts(rdkm);

	 if (v.size()) {
	    GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));
	    gtk_label_set_text(GTK_LABEL(lbg_alert_name_label), v[0].smarts_name.c_str());

	    gtk_widget_show(lbg_alert_hbox);

	    if (false) // debugging
	       for (unsigned int i=0; i<v.size(); i++) 
		  std::cout << "ALERT: " << v[i].smarts << " " << v[i].smarts_name << std::endl;

	    // if root gets renewed then this needs to be renewed too (so
	    // clear the alert group (and set it to null) when you clear the
	    // root.
	    //
	    if (alert_group == NULL) {
	       alert_group = goo_canvas_group_new (root, NULL);
	    } else {
	       clear_canvas_alerts();
	    } 
      
      
	    for (unsigned int imatch=0; imatch<v.size(); imatch++) {
	       for (unsigned int imatch_atom=0; imatch_atom<v[imatch].matches.size(); imatch_atom++) {
		  unsigned int rdkmol_idx = v[imatch].matches[imatch_atom].second;
		  // std::cout << "   " << rdkmol_idx;
		  // idx is some function of rdkmol_idx - look it up
		  std::string lbg_atom_index_str;
		  int lbg_atom_index = rdkmol_idx;
		  try { 
		     RDKit::ATOM_SPTR at_p = rdkm[rdkmol_idx];
		     at_p->getProp("lbg_atom_index", lbg_atom_index_str);
		     int lbg_atom_index = coot::util::string_to_int(lbg_atom_index_str);
		     lig_build::pos_t pos = mol.atoms[lbg_atom_index].atom_position;
		     double radius = 14;
		     double line_width = 0;
		     std::string col = "";
		     GooCanvasItem *circle = 
			goo_canvas_ellipse_new(alert_group,
					       pos.x, pos.y, radius, radius,
					       "line_width", line_width,
					       "fill-color-rgba", 0xffaa0050,
					       NULL);
		     goo_canvas_item_lower(circle, NULL); // to the bottom
		  }
		  catch (KeyErrorException kee) {
		     // actually, perhaps we don't get this we should fail, else we risk a crash?
		     // lbg_atom_index = rdkmol_idx;
		     std::cout << "ERROR:: in update_alerts() failed to get atom name " 
			       << kee.what() << std::endl;
		  }
	       }
	    }
	 } else {
	    gtk_widget_hide(lbg_alert_hbox);
	    clear_canvas_alerts();
	 }
      }
   }
}
#endif

void
lbg_info_t::clear_canvas_alerts() {
   if (alert_group) {
      gint n_children = goo_canvas_item_get_n_children (alert_group);
      // std::cout << "removing " << n_children << " children " << std::endl;
      for (int i=0; i<n_children; i++)
	 goo_canvas_item_remove_child(alert_group, 0);
   }
} 



// static
gboolean
lbg_info_t::watch_for_mdl_from_coot(gpointer user_data) {

   lbg_info_t *l = static_cast<lbg_info_t *> (user_data);
   
   // the pdb file has the atom names.  When we read the mol file, we
   // will stuff the atom names into the widgeted_molecule (by
   // matching the coordinates).
   //

   std::string coot_dir = l->get_flev_analysis_files_dir();
   std::string coot_ccp4_dir = coot_dir + "/coot-ccp4";
   
   // in use (non-testing), coot_dir will typically be ".";
   
   std::string ready_file    = coot_ccp4_dir + "/.coot-to-lbg-mol-ready";

   struct stat buf;

   // std::cout << "watching reading file: " << ready_file << std::endl;
   
   int err = stat(ready_file.c_str(), &buf);
   if (! err) {
      time_t m = buf.st_mtime;
      if (m > l->coot_mdl_ready_time) {
	 if (l->coot_mdl_ready_time != 0) {
	    l->read_files_from_coot();
	 }
	 l->coot_mdl_ready_time = m;
      }
   }
   return 1; // keep running
}


void
lbg_info_t::handle_read_draw_coords_mol_and_solv_acc(const std::string &coot_pdb_file,
						     const std::string &coot_mdl_file,
						     const std::string &sa_file) {
   
   CMMDBManager *flat_pdb_mol = get_cmmdbmanager(coot_pdb_file);
   lig_build::molfile_molecule_t mm;
   mm.read(coot_mdl_file);
   widgeted_molecule_t wmol = import_mol_file(mm, coot_mdl_file, flat_pdb_mol);
   std::vector<solvent_accessible_atom_t> solvent_accessible_atoms =
      read_solvent_accessibilities(sa_file);
   wmol.map_solvent_accessibilities_to_atoms(solvent_accessible_atoms);
   render_from_molecule(wmol);
}


CMMDBManager *
lbg_info_t::get_cmmdbmanager(const std::string &file_name) const {

   CMMDBManager *mol = new CMMDBManager;

   int err = mol->ReadCoorFile(file_name.c_str());
   if (err) {
      std::cout << "WARNING:: Problem reading coordinates file " << file_name << std::endl;
      delete mol;
      mol = NULL;
   } 

   return mol;

}



std::string
lbg_info_t::get_flev_analysis_files_dir() const {

   // in use (non-testing), coot_dir will typically be ".";
   
   std::string coot_dir = "../../build-coot-ubuntu-64bit/src"; 
   coot_dir = "../../build-lucid/src";

   const char *t_dir = getenv("COOT_LBG_OUT_DIR");
   if (t_dir)
      coot_dir = t_dir;

   // std::cout << "DEBUG:: get_flev_analysis_files_dir: " << coot_dir << std::endl;
//    if (t_dir)
//       std::cout << "DEBUG:: t_dir: " << t_dir << std::endl;
//    else 
//       std::cout << "DEBUG:: t_dir: NULL" << std::endl;
   return coot_dir;
}

void
lbg_info_t::read_files_from_coot() {

   std::string coot_dir = get_flev_analysis_files_dir();
	    
   std::string coot_ccp4_dir = coot_dir + "/coot-ccp4";
   // in use (non-testing), coot_dir will typically be ".";
	    
   std::string coot_mdl_file = coot_ccp4_dir + "/.coot-to-lbg-mol";
   std::string coot_pdb_file = coot_ccp4_dir + "/.coot-to-lbg-pdb";
   std::string ready_file    = coot_ccp4_dir + "/.coot-to-lbg-mol-ready";
   std::string sa_file       = coot_dir      + "/coot-tmp-fle-view-solvent-accessibilites.txt";

   std::cout << "trying to read " << coot_pdb_file << " "
	     << coot_mdl_file << std::endl;

   handle_read_draw_coords_mol_and_solv_acc(coot_pdb_file,
					    coot_mdl_file,
					    sa_file);
}


void
lbg_info_t::add_search_combobox_text() const {

   bool done_set_active = 0;
   GtkTreeIter   iter;
   GtkListStore *list_store_similarities =
      gtk_list_store_new (1, G_TYPE_STRING);
   gtk_combo_box_set_model(GTK_COMBO_BOX(lbg_search_combobox), GTK_TREE_MODEL(list_store_similarities));
   for (unsigned int i=0; i<6; i++) {
      double f = 0.75 + i*0.05;
      std::string s = coot::util::float_to_string(f);
      gtk_list_store_append(GTK_LIST_STORE(list_store_similarities), &iter);
      gtk_list_store_set(GTK_LIST_STORE(list_store_similarities), &iter,
			 0, s.c_str(), -1);
      // the default search similarity is 0.95 initially, so if we are
      // close to that, then make this the active widget.
      if (f < (search_similarity + 0.005)) { 
	 if (f > (search_similarity - 0.005)) {
	    done_set_active = 1;
	    gtk_combo_box_set_active(GTK_COMBO_BOX(lbg_search_combobox), i);
	 }
      }
   }
   if (! done_set_active)
      gtk_combo_box_set_active(GTK_COMBO_BOX(lbg_search_combobox), 0);

   GtkCellRenderer *renderer = gtk_cell_renderer_text_new ();
   gtk_cell_layout_pack_start (GTK_CELL_LAYOUT (lbg_search_combobox), renderer, TRUE);
   gtk_cell_layout_add_attribute (GTK_CELL_LAYOUT (lbg_search_combobox), renderer, "text", 0);
   

}

std::string
lbg_info_t::get_stroke_colour(int i, int n) const {

   std::string s("fedcba9876543210");

   double f =  1 + 12 * double(i)/double(n);
   int f_int = int(f);

   char c = s[f_int];
   std::string r = "#";
   for (i=0; i<6; i++)
      r += c;
   return r;
}

void
lbg_info_t::render_from_molecule(const widgeted_molecule_t &mol_in) {

   // std::cout << "render_from_molecule... " << std::endl;

   make_saves_mutex = 0; // stop saving changes (restored at end)
   clear();
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   
   int re_index[mol_in.atoms.size()]; // map from mol_in atom indexing
				      // the this molecule atom
				      // indexing.
   for (unsigned int i=0; i<mol_in.atoms.size(); i++)
      re_index[i] = UNASSIGNED_INDEX;

   // add in atoms
   for (unsigned int iat=0; iat<mol_in.atoms.size(); iat++) {
      if (!mol_in.atoms[iat].is_closed()) { 
	 GooCanvasItem *ci = NULL;
	 lig_build::pos_t pos = mol_in.atoms[iat].atom_position;
	 widgeted_atom_t new_atom = widgeted_atom_t(pos,
						    mol_in.atoms[iat].element,
						    mol_in.atoms[iat].charge,
						    ci);
	 double sa = mol_in.atoms[iat].get_solvent_accessibility();
	 new_atom.set_atom_name(mol_in.atoms[iat].get_atom_name());
	 if (0)
	    std::cout << "in render_from_molecule() atom with name :"
		      << mol_in.atoms[iat].get_atom_name()
		      << ": has solvent_accessibility " << sa << std::endl;
	 new_atom.add_solvent_accessibility(sa);
	 new_atom.bash_distances = mol_in.atoms[iat].bash_distances;
	 re_index[iat] = mol.add_atom(new_atom).second;
	 
	 if (0)
	    // clever c++
	    std::cout << " in render_from_molecule: old sa: "
		      << mol_in.atoms[iat].get_solvent_accessibility() << " new: "
		      << mol.atoms[re_index[iat]].get_solvent_accessibility() << std::endl;
// 	 if (sa > 0) {
// 	    std::cout << "draw solvent accessibility " << sa << " at " << pos << std::endl;
// 	    draw_solvent_accessibility_of_atom(pos, sa, root);
// 	 }
      }
   }


   // add in bonds
   for (unsigned int ib=0; ib<mol_in.bonds.size(); ib++) {
      int idx_1 = re_index[mol_in.bonds[ib].get_atom_1_index()];
      int idx_2 = re_index[mol_in.bonds[ib].get_atom_2_index()];
      if ((idx_1 != UNASSIGNED_INDEX) && (idx_2 != UNASSIGNED_INDEX)) { 
	 lig_build::bond_t::bond_type_t bt = mol_in.bonds[ib].get_bond_type();
	 if (mol_in.bonds[ib].have_centre_pos()) {
	    lig_build::pos_t centre_pos = mol_in.bonds[ib].centre_pos();
	    widgeted_bond_t bond(idx_1, idx_2, mol.atoms[idx_1], mol.atoms[idx_2], centre_pos, bt, root);
	    mol.add_bond(bond);
	 } else {
	    widgeted_bond_t bond(idx_1, idx_2, mol.atoms[idx_1], mol.atoms[idx_2], bt, root);
	    mol.add_bond(bond);
	 }
      }
   }

   // redo the atoms, this time with widgets.
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {

      std::vector<int> local_bonds = mol.bonds_having_atom_with_atom_index(iat);
      std::string ele = mol.atoms[iat].element;
      bool gl_flag = false; // not a GL render engine
      lig_build::atom_id_info_t atom_id_info =
	 mol.make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);

      std::string fc = font_colour(ele);
      if (ele != "C") 
	 mol.atoms[iat].update_atom_id_forced(atom_id_info, fc, root);
   }

   // for input_coords_to_canvas_coords() to work:
   //

   draw_substitution_contour();
   
   mol.centre_correction = mol_in.centre_correction;
   mol.scale_correction  = mol_in.scale_correction;
   mol.mol_in_min_y = mol_in.mol_in_min_y;
   mol.mol_in_max_y = mol_in.mol_in_max_y;
   
   // 
   make_saves_mutex = 1; // allow saves again.

}

void
lbg_info_t::undo() {

   save_molecule_index--;
   if (save_molecule_index >= 0) { 
      widgeted_molecule_t saved_mol = previous_molecules[save_molecule_index];
      // std::cout << "undo... reverting to save molecule number " << save_molecule_index << std::endl;
      render_from_molecule(saved_mol);
      update_descriptor_attributes();
   } else {
      clear();
   } 
} 

void
lbg_info_t::delete_hydrogens() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   // I don't understand what is going on here. Why can't I do this?
   // mol.delete_hydrogens();  render(mol); ?  -> Blank canvas
   widgeted_molecule_t copy_mol = mol;
   copy_mol.delete_hydrogens(root);
   save_molecule();
   render_from_molecule(copy_mol);
} 



void
lbg_info_t::write_pdf(const std::string &file_name) const { 

#if CAIRO_HAS_PDF_SURFACE

   cairo_surface_t *surface;
   cairo_t *cr;

   std::pair<lig_build::pos_t, lig_build::pos_t> extents = mol.ligand_extents();
   double pos_x = (extents.second.x + 220.0);
   double pos_y = (extents.second.y + 220.0);
   surface = cairo_pdf_surface_create(file_name.c_str(), pos_x, pos_y);
   cr = cairo_create (surface);

   /* Place it in the middle of our 9x10 page. */
   // cairo_translate (cr, 20, 130);
   cairo_translate (cr, 2, 13);

   goo_canvas_render (GOO_CANVAS(canvas), cr, NULL, 1.0);
   cairo_show_page (cr);
   cairo_surface_destroy (surface);
   cairo_destroy (cr);

#else
   std::cout << "No PDF (no PDF Surface in Cairo)" << std::endl;
#endif    

}

void
lbg_info_t::write_png(const std::string &file_name) const {

   std::pair<lig_build::pos_t, lig_build::pos_t> extents = mol.ligand_extents();

   gdouble scale =  goo_canvas_get_scale(GOO_CANVAS(canvas));
   // doesn't do anything
   // goo_canvas_set_scale(GOO_CANVAS(canvas), 4*scale);
   // goo_canvas_update(GOO_CANVAS(canvas));
   int size_x = (int)extents.second.x + 220; // or so... (ideally should reside circle-based).
   int size_y = (int)extents.second.y + 220;
   
   cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, size_x, size_y);
   cairo_t *cr = cairo_create (surface);

   // 1.0 is item's visibility threshold to see if they should be rendered
   // (everything should be rendered).
   // 
   goo_canvas_render (GOO_CANVAS(canvas), cr, NULL, 1.0); 
   cairo_surface_write_to_png(surface, file_name.c_str());
   cairo_surface_destroy (surface);
   cairo_destroy (cr);
   goo_canvas_set_scale(GOO_CANVAS(canvas), scale);
   // goo_canvas_update(GOO_CANVAS(canvas));
}



void
lbg_info_t::save_molecule() {

   if (make_saves_mutex) { 
      previous_molecules.push_back(mol);
      save_molecule_index = previous_molecules.size() - 1;
      if (0)
	 std::cout << "saved molecule to index: " << save_molecule_index << std::endl;
      
   } else {
      std::cout << "debug:: save_molecule() excluded" << std::endl;
   }
}

void
lbg_info_t::import_mol_from_file(const std::string &file_name) {

   bool try_as_mdl_mol = false;
#ifndef MAKE_ENTERPRISE_TOOLS
   // fallback
   try_as_mdl_mol = true; 
#else    
   try { 
      RDKit::RWMol *m = RDKit::Mol2FileToMol(file_name);
      if (m) {
	 // molfile molecules don't know about aromatic bonds, we need
	 // to kekulize now.
	 RDKit::MolOps::Kekulize(*m); // non-const reference?
	 double weight_for_3d_distances = 0.4; 
	 int iconf = 0; // current conformer number
	 int iconf_new = coot::add_2d_conformer(m, weight_for_3d_distances);
	 // could we add another conformer (which is 2d), if so, use
	 // that for layout instead.
	 if (iconf_new > iconf) {
	    iconf = iconf_new;
	 } else {
	    std::cout << "WARNING:: import_mol_from_file() failed to make 2d conformer "
		      << std::endl;
	 }
	 lig_build::molfile_molecule_t mm = coot::make_molfile_molecule(*m, iconf);
	 CMMDBManager *mol = NULL; // no atom names to transfer
	 widgeted_molecule_t wmol = import_mol_file(mm, file_name, mol);
	 render_from_molecule(wmol);
	 std::cout << "in import_mol_from_file() calling update_descriptor_attributes() "
		   << std::endl;
	 update_descriptor_attributes();

      } else {
	 // should throw an exception before getting here, I think.
	 std::cout << "Null m in import_mol_from_file() " << std::endl;
	 try_as_mdl_mol = true;
      } 
   }
   catch (RDKit::FileParseException rte) {
      std::cout << "WARNING:: " << rte.message() << std::endl;
      try_as_mdl_mol = true;
   }
   catch (RDKit::BadFileException &e) {
      std::cout << "WARNING:: Bad file " << file_name << " " << e.message() << std::endl;
      try_as_mdl_mol = true;
   }
   catch (std::runtime_error rte) {
      std::cout << "WARNING runtime_error in mol_to_asc_rdkit() " << rte.what() << std::endl;
      try_as_mdl_mol = true;
   } 
   catch (std::exception e) {
      std::cout << "WARNING:: import_mol_from_file: exception: " << e.what() << std::endl;
      try_as_mdl_mol = true;
   }
#endif // MAKE_ENTERPRISE_TOOLS

   if (try_as_mdl_mol) { 
      // read as an MDL mol file
      // 
      CMMDBManager *mol = NULL; // no atom names to transfer
      lig_build::molfile_molecule_t mm;
      mm.read(file_name);
      widgeted_molecule_t wmol = import_mol_file(mm, file_name, mol);
      render_from_molecule(wmol);
      std::cout << "in import_mol_from_file() 2 calling update_descriptor_attributes() "
		<< std::endl;
      update_descriptor_attributes();
   }
   
}

// pdb_mol is the pdb representation of the (flat) ligand - and it has
// the atom names.  We will add the atom names into mol by matching
// coordinates.
// 
widgeted_molecule_t
lbg_info_t::import_mol_file(const lig_build::molfile_molecule_t &mol_in,
			    const std::string &file_name,
			    CMMDBManager *pdb_mol) {

   widgeted_molecule_t new_mol(mol_in, pdb_mol);
   mdl_file_name = file_name;
   save_molecule();
   return new_mol;
}


void
lbg_info_t::read_draw_residues(const std::string &file_name) {

   bool show_dynamics = 0;

   // read residues need to happen after the ligand has been placed on the canvas
   // i.e. mol.input_coords_to_canvas_coords() is called in read_residues();
   // 
   std::vector<lbg_info_t::residue_circle_t> file_residue_circles = read_residues(file_name);

   double max_dist_water_to_ligand_atom  = 3.2; // don't draw waters that are far from ligand
   double max_dist_water_to_protein_atom = 3.2; // don't draw waters that are not somehow 
                                                // attached to the protein.

   residue_circles = filter_residue_waters(file_residue_circles,
					   max_dist_water_to_ligand_atom,
					   max_dist_water_to_protein_atom);

   std::vector<int> primary_indices = get_primary_indices();

   initial_residues_circles_layout(); // twiddle residue_circles
   
   // cycles of optimisation
   // 
   std::vector<residue_circle_t> current_circles = residue_circles;
   std::vector<int> empty_vector; // for add_rep_handles, can't do it from a file.
   if (show_dynamics)
      draw_residue_circles(current_circles, empty_vector);

   // For debugging the minimizer (vs original positions).
   // 
   current_circles = offset_residues_from_orig_positions();
   
   for (int iround=0; iround<120; iround++) {
      std::cout << ":::::::::::::::::::: minimization round " << iround
		<< " :::::::::::::::::::::::::::::::::::::::::::\n";
      std::pair<int, std::vector<lbg_info_t::residue_circle_t> > new_c =
	 optimise_residue_circle_positions(residue_circles, current_circles, primary_indices);
      current_circles = new_c.second;
      if (show_dynamics)
	 draw_residue_circles(current_circles, empty_vector);
      if (new_c.first == GSL_ENOPROG)
	 break;
      if (new_c.first == GSL_SUCCESS) { 
	 break;
      }
   }

   // save the class member data
   residue_circles = current_circles;

   // has the current solution problems due to residues too close to the ligand?
   std::pair<bool, std::vector<int> > problem_status = solution_has_problems_p();
   std::cout << "::::::::: problem status: " << problem_status.first << std::endl;
   for (unsigned int ip=0; ip<problem_status.second.size(); ip++) {
      std::cout << ":::::::::::: "
		<< residue_circles[problem_status.second[ip]].residue_label << " "
		<< residue_circles[problem_status.second[ip]].residue_type << " "
		<< residue_circles[problem_status.second[ip]].pos
		<< std::endl;
   }
   
   if (problem_status.second.size()) {
      
      // fiddle with residue_circles and reoptimise.
      // 
      reposition_problematics_and_reoptimise(problem_status.second,
					     primary_indices);
   }
   
   draw_all_flev_annotations(); // contour debugging.
}

std::vector<lbg_info_t::residue_circle_t>
lbg_info_t::filter_residue_waters(const std::vector<lbg_info_t::residue_circle_t> &r_in,
				  double max_dist_water_to_ligand_atom,
				  double max_dist_water_to_protein_atom) const {

   std::vector<lbg_info_t::residue_circle_t> v;

   for (unsigned int i=0; i<r_in.size(); i++) { 
      if (r_in[i].residue_type != "HOH") {
	 v.push_back(r_in[i]);
      } else { 
	 if (r_in[i].water_dist_to_protein > max_dist_water_to_protein_atom) {
	    // pass
	 } else {
	    bool all_too_long_distance = 1;
	    std::vector<lbg_info_t::bond_to_ligand_t> sasisfactory_bonds;
	    for (unsigned int j=0; j<r_in[i].bonds_to_ligand.size(); j++) {
	       if (r_in[i].bonds_to_ligand[j].bond_length < max_dist_water_to_ligand_atom) { 
		  all_too_long_distance = 0;
		  sasisfactory_bonds.push_back(r_in[i].bonds_to_ligand[j]);
	       }
	    }
	    if (all_too_long_distance) {
	       // pass
	    } else {
	       // replace input waters with sasisfactory_bonds.
	       // 
	       lbg_info_t::residue_circle_t n = r_in[i];
	       n.bonds_to_ligand = sasisfactory_bonds;
	       v.push_back(n);
	    } 
	 } 
      }
   }

   return v;
} 


// fiddle with the position of some residues in residue_circles
//
void
lbg_info_t::reposition_problematics_and_reoptimise(const std::vector<int> &problematics,
						   const std::vector<int> &primary_indices) {

   std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair =
      mol.ligand_extents();
   lbg_info_t::ligand_grid grid(l_e_pair.first, l_e_pair.second);
   grid.fill(mol);
   for (unsigned int ip=0; ip<problematics.size(); ip++) {

      // a trapped residue is now treated as a primary (bonding)
      // residue with a distance to the "bad" solution of 3.5A (this
      // distance does not matter) - for initial positioning.  It is
      // not added to the primaries, so this bond length is not
      // optimised.
      
      std::vector<std::pair<lig_build::pos_t, double> > attachment_points;
      lig_build::pos_t attach_pos = residue_circles[problematics[ip]].pos;
      attach_pos+= lig_build::pos_t (5* double(rand())/double(RAND_MAX),
				     5* double(rand())/double(RAND_MAX));
      
      std::pair<lig_build::pos_t, double> p(attach_pos, 3.5);

      attachment_points.push_back(p);
      initial_primary_residue_circles_layout(grid, problematics[ip],
					     attachment_points);
   }

   
   // set the initial positions for optimisation, now that we have
   // fiddled with residue_circles
   std::vector<residue_circle_t> current_circles = residue_circles;

   for (int iround=0; iround<120; iround++) {
      std::pair<int, std::vector<lbg_info_t::residue_circle_t> > new_c =
	 optimise_residue_circle_positions(residue_circles, current_circles, primary_indices);
      current_circles = new_c.second;
      if (new_c.first == GSL_ENOPROG)
	 break;
      if (new_c.first == GSL_SUCCESS)
	 break;
   }
   // now transfer results from the updating/tmp current_circles to the class data item:
   residue_circles = current_circles;
}

void
lbg_info_t::recentre_considering_residue_centres() {

   std::pair<bool,lig_build::pos_t> tl = get_residue_circles_top_left();

   std::cout << "----------------------- top left correction " << tl.first << " " << tl.second << std::endl;
  
    if (tl.first) { 
       top_left_correction = lig_build::pos_t(-tl.second.x+40, -tl.second.y+30);
       clear_and_redraw(top_left_correction);
    }
} 


// for debugging the minimization vs original positions
std::vector<lbg_info_t::residue_circle_t>
lbg_info_t::offset_residues_from_orig_positions() {

   std::vector<residue_circle_t> current_circles = residue_circles;
   for (unsigned int i=0; i<current_circles.size(); i++) {
      double delta_x = 5 * double(rand())/double(RAND_MAX);
      double delta_y = 5 * double(rand())/double(RAND_MAX);
      current_circles[i].pos.x += delta_x;
      current_circles[i].pos.y += delta_y;
   }
   return current_circles;
}

void
lbg_info_t::draw_all_flev_annotations() {

   draw_all_flev_residue_attribs();
   draw_all_flev_ligand_annotations();
} 

void
lbg_info_t::draw_all_flev_residue_attribs() {

   if (draw_flev_annotations_flag) {
      draw_residue_circles(residue_circles, additional_representation_handles);
      draw_bonds_to_ligand();
      draw_stacking_interactions(residue_circles);
   }
}

void
lbg_info_t::draw_all_flev_ligand_annotations() {

   if (draw_flev_annotations_flag) {
      draw_substitution_contour();
      draw_solvent_accessibility_of_atoms();
   } 
} 

std::vector<int>
lbg_info_t::get_primary_indices() const {

   std::vector<int> primary_indices;  // this primary_indices needs to
                                      // get passed to the
                                      // primary_indices used in
                                      // residue cirlce optimization.
   
   for(unsigned int ic=0; ic<residue_circles.size(); ic++) {
      if (residue_circles[ic].is_a_primary_residue()) {
         primary_indices.push_back(ic);
      }
   }
   return primary_indices;
} 


// twiddle residue_circles, taking in to account the residues that
// bond to the ligand (including stacking residues).
// 
void
lbg_info_t::initial_residues_circles_layout() {

   // when we move a primary, we want to know it's index in
   // residue_circles, because that's what we really want to move.
   // 
   // std::vector<std::pair<int, lbg_info_t::residue_circle_t> > primaries;

   // now a class data member because it is used in the layout penalty
   // function (we want to have nice bond lengths for residues bonded
   // to the ligand).
   // 
   std::vector<int> primary_indices;  // this primary_indices needs to
				      // get passed to the
				      // primary_indices used in
				      // residue circle optimization.
   
   for(unsigned int ic=0; ic<residue_circles.size(); ic++) {
      if (residue_circles[ic].is_a_primary_residue()) {
	 primary_indices.push_back(ic);
      }
   }

   // primaries get placed first and are checked for non-crossing
   // ligand interaction bonds (they have penalty scored added if they
   // cross).
   //
   try { 
      std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair =
	 mol.ligand_extents();
      lbg_info_t::ligand_grid grid(l_e_pair.first, l_e_pair.second);
      grid.fill(mol);

      for (unsigned int iprimary=0; iprimary<primary_indices.size(); iprimary++) {
	 int idx = primary_indices[iprimary];
	 std::vector<std::pair<lig_build::pos_t, double> > attachment_points =
	    residue_circles[idx].get_attachment_points(mol);
	 initial_primary_residue_circles_layout(grid, idx, attachment_points);
      }
      // show_mol_ring_centres();

      position_non_primaries(grid, primary_indices); // untrap residues as needed.
   }
   catch (std::runtime_error rte) {
      std::cout << rte.what() << std::endl;
   }
}

// fiddle with the position of the residue_circles[primary_index].
//
void
lbg_info_t::initial_primary_residue_circles_layout(const lbg_info_t::ligand_grid &grid,
						   int primary_index,
						   const std::vector<std::pair<lig_build::pos_t, double> > &attachment_points) {

   if (0) 
      std::cout << "DEBUG:: starting initial_primary_residue_circles_layout() primary_index " << primary_index
		<< " " << residue_circles[primary_index].residue_label << " "
		<< residue_circles[primary_index].residue_type
		<< " has position " << residue_circles[primary_index].pos
		<< std::endl;
   
   if (0)
      std::cout << " =========== adding quadratic for residue " 
		<< residue_circles[primary_index].residue_label
		<< " ============================"
		<< std::endl;
   lbg_info_t::ligand_grid primary_grid = grid;
	 
   // attachment points are points on the ligand, in ligand
   // coordinates to which this primary residue circle is
   // attached (often only one attachment, but can be 2
   // sometimes).
   primary_grid.add_quadratic(attachment_points);

   if (0)
      if (primary_index == 0) 
	 show_grid(grid);
   
   lig_build::pos_t best_pos = primary_grid.find_minimum_position();

   // OK, consider the case where there are 2 residues bonding to the
   // same atom in the ligand.  They wil both be given exactly the
   // same best_pos and they will subsequently refine together,
   // resulting in one residue sitting on top of another - bad news.
   //
   // So let's shift the residue a bit in the direction that it came
   // from so that the residues don't refine together.
   //
   lig_build::pos_t a_to_b_uv = (residue_circles[primary_index].pos - best_pos).unit_vector();
   
   residue_circles[primary_index].pos = best_pos + a_to_b_uv * 4;

   if (0)
      std::cout << "DEBUG::  ending initial_primary_residue_circles_layout() primary_index "
		<< primary_index
		<< " " << residue_circles[primary_index].residue_label << " "
		<< residue_circles[primary_index].residue_type
		<< " has position " << residue_circles[primary_index].pos
		<< std::endl;
   
}

// untrap residues as needed.
void
lbg_info_t::position_non_primaries(const lbg_info_t::ligand_grid &grid,
				   const std::vector<int> &primary_indices) {

   std::vector<lbg_info_t::grid_index_t> already_positioned;
   
   for (unsigned int irc=0; irc<residue_circles.size(); irc++) {
      // if not a primary...
      if (std::find(primary_indices.begin(), primary_indices.end(), irc) ==
	  primary_indices.end()) {
	 
	 // if this point has rejection underneath it, find the
	 // position in grid that doesn't have any rejection.
	 //

	 std::pair<lbg_info_t::grid_index_t, lig_build::pos_t> pos =
	    grid.find_nearest_zero(residue_circles[irc].pos, already_positioned);
	 residue_circles[irc].pos = pos.second;
	 if (pos.first.is_valid_p())
	    already_positioned.push_back(pos.first);
			 
      }
   }
}


// attachment points are points on the ligand, in ligand coordinates
// to which this primary residue circle is attached (often only one
// attachment, but can be 2 sometimes).
//
void
lbg_info_t::ligand_grid::add_quadratic(const std::vector<std::pair<lig_build::pos_t, double> > &attachment_points) {

   if (attachment_points.size()) {
      double scale_by_n_attach = 1.0/double(attachment_points.size());
      
      for (unsigned int iattach=0; iattach<attachment_points.size(); iattach++) {
	 for (unsigned int ix=0; ix<x_size(); ix++) {
	    for (unsigned int iy=0; iy<y_size(); iy++) {
	       lig_build::pos_t pos = to_canvas_pos(ix, iy);
	       double d2 = (pos-attachment_points[iattach].first).lengthsq();
	       double val = 0.00002 * d2 * scale_by_n_attach;
	       grid_[ix][iy] += val;
	    }
	 }
      }
   }
}

// can throw an exception
lig_build::pos_t
lbg_info_t::ligand_grid::find_minimum_position() const {

   double best_pos_score = 1000000;
   lig_build::pos_t best_pos;
   for (unsigned int ix=0; ix<x_size(); ix++) {
      for (unsigned int iy=0; iy<y_size(); iy++) {
	 if (grid_[ix][iy] < best_pos_score) {
	    best_pos_score = grid_[ix][iy];
	    best_pos = to_canvas_pos(ix,iy);
	 }
      }
   }
   if (best_pos_score > (1000000-1))
      throw std::runtime_error("failed to get minimum position from ligand grid");
   return best_pos;
}

// actually, not exactly zero but something small.
//
// Don't return a grid-point/position that matches anything in
// already_positioned.
// 
std::pair<lbg_info_t::grid_index_t, lig_build::pos_t>
lbg_info_t::ligand_grid::find_nearest_zero(const lig_build::pos_t &pos,
					   const std::vector<lbg_info_t::grid_index_t> &already_positioned) const {

   lig_build::pos_t p;
   lbg_info_t::grid_index_t rgi; // initially invalid
   double shortest_dist = 43e23;
   double crit = 0.05; // less than this is effectively zero.

   try { 
      lbg_info_t::grid_index_t gi=grid_pos_nearest(pos);
      if (grid_[gi.i()][gi.j()] < crit) {
	 p = pos; // fine, no change
      } else {
	 // search for someplace else
	 for (unsigned int ix=0; ix<x_size(); ix++) {
	    for (unsigned int iy=0; iy<y_size(); iy++) {
	       if (grid_[ix][iy] < crit) {
		  lig_build::pos_t gp = to_canvas_pos(ix, iy);
		  double d = (gp - pos).lengthsq();
		  if (d < shortest_dist) {
		     lbg_info_t::grid_index_t candidate_grid_index(ix, iy);
		     // This is OK if there is no other previous
		     // solution at the same position (if there is,
		     // keep trying, of course).
		     if (std::find(already_positioned.begin(), already_positioned.end(),
				   candidate_grid_index) == already_positioned.end()) { 
			shortest_dist = d;
			p = gp;
			rgi = candidate_grid_index;
		     }
		  } 
	       }
	    }
	 }
      }
   }
   catch (std::runtime_error rte) {
      // the pos was off the grid.  It won't be trapped inside the
      // ligand, so just return what we were given.
      p = pos;
   } 
   return std::pair<lbg_info_t::grid_index_t, lig_build::pos_t> (rgi, p);
} 


lig_build::pos_t
lbg_info_t::ligand_grid::to_canvas_pos(const double &ii, const double &jj) const {

   lig_build::pos_t p(ii*scale_fac, jj*scale_fac);
   p += top_left;
   return p;
}

// can throw a std::runtime_error.
// 
lbg_info_t::grid_index_t
lbg_info_t::ligand_grid::grid_pos_nearest(const lig_build::pos_t &pos) const {

   lig_build::pos_t p = pos - top_left;
   int idx_x = int(p.x/scale_fac+0.5);
   int idx_y = int(p.y/scale_fac+0.5);

   if ((idx_x < 0) || (idx_x >= x_size()) || (idx_y < 0) || (idx_y >= y_size()))
       throw std::runtime_error("out of grid index");
   
   return lbg_info_t::grid_index_t(idx_x, idx_y);
} 

 


std::pair<int, int>
lbg_info_t::ligand_grid::canvas_pos_to_grid_pos(const lig_build::pos_t &pos) const {

   lig_build::pos_t p = pos - top_left;
   int ix = (int)(p.x/scale_fac);
   int iy = (int)(p.y/scale_fac);
   return std::pair<int, int> (ix, iy);
}

// for marching squares, ii and jj are the indices of the bottom left-hand side.
int 
lbg_info_t::ligand_grid::square_type(int ii, int jj, float contour_level) const {

   int square_type = lbg_info_t::ligand_grid::MS_NO_SQUARE;
   if ((ii+1) >= x_size_) {
      return lbg_info_t::ligand_grid::MS_NO_SQUARE;
   } else { 
      if ((jj+1) >= y_size_) {
	 return lbg_info_t::ligand_grid::MS_NO_SQUARE;
      } else {
	 float v00 = get(ii, jj);
	 float v01 = get(ii, jj+1);
	 float v10 = get(ii+1, jj);
	 float v11 = get(ii+1, jj+1);
	 if (v00 > contour_level) { 
	    if (v01 > contour_level) { 
	       if (v10 > contour_level) { 
		  if (v11 > contour_level) {
		     return lbg_info_t::ligand_grid::MS_NO_CROSSING;
		  }
	       }
	    }
	 }
	 if (v00 < contour_level) { 
	    if (v01 < contour_level) { 
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_NO_CROSSING;
		  }
	       }
	    }
	 }

	 // OK, so it is not either of the trivial cases (no
	 // crossing), there are 14 other variants.
	 // 
	 if (v00 < contour_level) { 
	    if (v01 < contour_level) { 
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_NO_CROSSING;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_1_1;
		  }
	       } else {
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_1_0;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_1_0_and_1_1;
		  }
	       }
	    } else {

	       // 0,1 is up
	       
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_1;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_1_and_1_1;
		  }
	       } else {
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_1_and_1_0;      // hideous valley
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_1_and_1_0_and_1_1; // (only 0,0 down)
		  }
	       }
	    }
	 } else {

	    // 0,0 is up
	    
	    if (v01 < contour_level) { 
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_0;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_1_1; // another hideous valley
		  }
	       } else {
		  // 1,0 is up
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_1_0;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_1_0_and_1_1; // 0,1 is down
		  }
	       }
	    } else {

	       // 0,1 is up
	       
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1_and_1_1; // 1,0 is down
		  }
	       } else {
		  // if we get here, this test must pass.
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1_and_1_0; // only 1,1 is down
		  }
	       }
	    }
	 } 
      }
   }
   return square_type;
} 

lbg_info_t::contour_fragment::contour_fragment(int ms_type,
					       const float &contour_level, 
					       const lbg_info_t::grid_index_t &grid_index_prev,
					       const lbg_info_t::grid_index_t &grid_index,
					       const lbg_info_t::ligand_grid &grid) {

   int ii_next = lbg_info_t::grid_index_t::INVALID_INDEX;
   int jj_next = lbg_info_t::grid_index_t::INVALID_INDEX;

   float v00 = grid.get(grid_index.i(),   grid_index.j());
   float v01 = grid.get(grid_index.i(),   grid_index.j()+1);
   float v10 = grid.get(grid_index.i()+1, grid_index.j());
   float v11 = grid.get(grid_index.i()+1, grid_index.j()+1);

   float frac_x1 = -1; 
   float frac_y1 = -1;
   float frac_x2 = -1;  // for hideous valley
   float frac_y2 = -1;

   lbg_info_t::contour_fragment::coordinates c1(0.0, X_AXIS_LOW);
   lbg_info_t::contour_fragment::coordinates c2(Y_AXIS_LOW, 0.0);
   cp_t p(c1,c2);
   
   switch (ms_type) {

   case lbg_info_t::ligand_grid::MS_UP_0_0:
   case lbg_info_t::ligand_grid::MS_UP_0_1_and_1_0_and_1_1:

      frac_x1 = (v00-contour_level)/(v00-v10);
      frac_y1 = (v00-contour_level)/(v00-v01);
      c1 = lbg_info_t::contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_LOW, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case lbg_info_t::ligand_grid::MS_UP_0_1:
   case lbg_info_t::ligand_grid::MS_UP_0_0_and_1_0_and_1_1:
      
      // these look fine
      // std::cout << " ----- case MS_UP_0,1 " << std::endl;
      frac_y1 = (v00 - contour_level)/(v00-v01);
      frac_x2 = (v01 - contour_level)/(v01-v11);
      c1 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_LOW, frac_y1);
      c2 = lbg_info_t::contour_fragment::coordinates(frac_x2, X_AXIS_HIGH);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

      
   case lbg_info_t::ligand_grid::MS_UP_1_0:
   case lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1_and_1_1:
      
      // std::cout << " ----- case MS_UP_1,0 " << std::endl;
      frac_x1 = (contour_level - v00)/(v10-v00);
      frac_y1 = (contour_level - v10)/(v11-v10);
      c1 = lbg_info_t::contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_HIGH, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

      
      
   case lbg_info_t::ligand_grid::MS_UP_1_1:
   case lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1_and_1_0:

      // std::cout << " ----- case MS_UP_1,1 " << std::endl;
      frac_x1 = (v01-contour_level)/(v01-v11);
      frac_y1 = (v10-contour_level)/(v10-v11);
      c1 = lbg_info_t::contour_fragment::coordinates(frac_x1, X_AXIS_HIGH);
      c2 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_HIGH, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1:
   case lbg_info_t::ligand_grid::MS_UP_1_0_and_1_1:
      
      // std::cout << " ----- case MS_UP_0,0 and 0,1 " << std::endl;
      frac_x1 = (v00-contour_level)/(v00-v10);
      frac_x2 = (v01-contour_level)/(v01-v11);
      c1 = lbg_info_t::contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = lbg_info_t::contour_fragment::coordinates(frac_x2, X_AXIS_HIGH);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case lbg_info_t::ligand_grid::MS_UP_0_0_and_1_0:
   case lbg_info_t::ligand_grid::MS_UP_0_1_and_1_1:
      
      // std::cout << " ----- case MS_UP_0,0 and 1,0 " << std::endl;
      frac_y1 = (v00-contour_level)/(v00-v01);
      frac_y2 = (v10-contour_level)/(v10-v11);
      c1 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_LOW,  frac_y1);
      c2 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_HIGH, frac_y2);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;
      

   default:
      std::cout << "ERROR:: unhandled square type: " << ms_type << std::endl;

   } 

} 

std::vector<std::pair<lig_build::pos_t, double> >
lbg_info_t::residue_circle_t::get_attachment_points(const widgeted_molecule_t &mol) const {

   // debug
//    std::cout << "........ in get_attachment_points() mol has " << mol.atoms.size()
// 	     << " atoms " << std::endl;
//    for (unsigned int i=0; i<mol.atoms.size(); i++) { 
//       std::cout << "    " << i << "   " << mol.atoms[i] << " "
// 		<< mol.atoms[i].get_atom_name() <<  std::endl;
//    }
   
   
   std::vector<std::pair<lig_build::pos_t, double> > v;

   for (unsigned int i=0; i<bonds_to_ligand.size(); i++) {
      if (bonds_to_ligand[i].is_set()) {
	 try {
	    lig_build::pos_t pos =
	       mol.get_atom_canvas_position(bonds_to_ligand[i].ligand_atom_name);
	    if (bonds_to_ligand[i].is_set()) { 
	       std::pair<lig_build::pos_t, double> p(pos, bonds_to_ligand[i].bond_length);
	       v.push_back(p);
	    }
	 }
	 catch (std::runtime_error rte) {
	    std::cout << "WARNING:: " << rte.what() << std::endl;
	 }
      }
   }

   // with a ring system on the ligand, that is.
   // 
   if (has_ring_stacking_interaction()) {
      try {
	 lig_build::pos_t pos = mol.get_ring_centre(ligand_ring_atom_names);
	 double stacking_dist = 4.5; // A
	 std::pair<lig_build::pos_t, double> p(pos, stacking_dist);
	 v.push_back(p);
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }

   if (get_stacking_type() == lbg_info_t::residue_circle_t::CATION_PI_STACKING) {
      try {
	 std::string at_name = get_ligand_cation_atom_name();
	 double stacking_dist = 4.2; // a bit shorter, because we
				     // don't have to go from the
				     // middle of a ring system, the
				     // target point (an atom) is more
				     // accessible.  Still 3.5 was too
				     // short (->crowded) when
				     // showing 3 cation-pi interactions
				     // on a alkylated N.
	 
	 lig_build::pos_t pos = mol.get_atom_canvas_position(at_name);
	 std::pair<lig_build::pos_t, double> p(pos, stacking_dist);
	 v.push_back(p);
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }

   return v;
}


// return 1 on solution having problems, 0 for no problems, also
// return a list of the residue circles with problems.
// 
std::pair<bool, std::vector<int> >
lbg_info_t::solution_has_problems_p() const {

   std::vector<int> v;
   bool status = 0;

   // scale factor of 1.2 is too small, 2.0 seems good for H102 in 2aei 
   double crit_dist_2 = 2.0 * SINGLE_BOND_CANVAS_LENGTH * SINGLE_BOND_CANVAS_LENGTH;
   double crit_dist_wide_2 = 4 * crit_dist_2;

   // a problem can be 2 atoms < crit_dist_2 or
   //                  5 atoms < crit_dist_wide_2
   // 
   // we need the wider check to find atoms that are enclosed by a
   // system of 10 or so atoms.
   // 

   for (unsigned int i=0; i<residue_circles.size(); i++) {
      int n_close = 0;
      int n_close_wide = 0;
      for (unsigned int j=0; j<mol.atoms.size(); j++) {
	 double d2 = (residue_circles[i].pos-mol.atoms[j].atom_position).lengthsq();
//  	 std::cout << "comparing " << d2 << " " << crit_dist_2 << " "
//  		   << residue_circles[i].residue_label
//  		   << std::endl;
	 if (d2 < crit_dist_2) {
	    n_close++;
	    if (n_close > 1) {
	       v.push_back(i);
	       status = 1;
	       break;
	    }
	 }

	 if (d2 < crit_dist_wide_2) {
	    n_close_wide++;
	    if (n_close_wide > 5) {
	       v.push_back(i);
	       status = 1;
	       break;
	    } 
	 } 
      }
   }
   return std::pair<bool, std::vector<int> > (status, v);
} 
 



void
lbg_info_t::show_grid(const lbg_info_t::ligand_grid &grid) {

   int n_objs = 0;
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   for (int ix=0; ix<grid.x_size(); ix+=1) {
      for (int iy=0; iy<grid.y_size(); iy+=1) {
	 lig_build::pos_t pos = grid.to_canvas_pos(ix, iy);
	 double intensity = grid.get(ix, iy);
	 int int_as_int = int(255 * intensity);
	 std::string colour = grid_intensity_to_colour(int_as_int);
 	 GooCanvasItem *rect =
 	    goo_canvas_rect_new(root, pos.x, pos.y, 3.5, 3.5,
 				"line-width", 1.0,
 				"fill_color", colour.c_str(),
				"stroke-color", colour.c_str(),
 				NULL);
	 n_objs++;
      }
   }

   if (1) {  // debugging
      grid.show_contour(root, 0.1);
      grid.show_contour(root, 0.2);
      grid.show_contour(root, 0.3);
      grid.show_contour(root, 0.4);
      grid.show_contour(root, 0.5);
      grid.show_contour(root, 0.6);
      grid.show_contour(root, 0.7);
      grid.show_contour(root, 0.8);
      grid.show_contour(root, 0.9);
      grid.show_contour(root, 1.0);
   }
   
}

void
lbg_info_t::ligand_grid::avoid_ring_centres(std::vector<std::vector<std::string> > ring_atoms_list,
					    const widgeted_molecule_t &mol) {

   // For the substitution contour we blob in a circle of radius
   // 1/(2*sin(180/n_ring_atoms)) bond lengths.
   
   for (unsigned int iring=0; iring<ring_atoms_list.size(); iring++) {
      try { 
	 lig_build::pos_t centre = mol.get_ring_centre(ring_atoms_list[iring]);
	 int n_atoms = ring_atoms_list.size();
	 double radius = 1/(2*sin(M_PI/double(n_atoms))) * 1.5; // in "A" or close
	 add_for_accessibility(radius, centre);
      }
      catch (std::runtime_error rte) {
	 std::cout << "Opps - failed to find ring centre for ring atom name "
		   << iring << std::endl;
      }
   }
}

void
lbg_info_t::ligand_grid::show_contour(GooCanvasItem *root, float contour_level) const {

   std::vector<widgeted_atom_ring_centre_info_t> dummy_unlimited_atoms;
   std::vector<std::vector<std::string> > dummy_ring_atom_names;   
   show_contour(root, contour_level, dummy_unlimited_atoms, dummy_ring_atom_names);
}

void
lbg_info_t::ligand_grid::show_contour(GooCanvasItem *root, float contour_level,
				      const std::vector<widgeted_atom_ring_centre_info_t> &unlimited_atoms,
				      const std::vector<std::vector<std::string> > &ring_atoms_list) const {


   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", "#880000",
						NULL);
   bool debug = 0;
   int ii=0;
   int jj=0;

   // fill the ring centre vector, if the unlimited atom have ring centres
   std::vector<std::pair<bool, lig_build::pos_t> > ring_centres(unlimited_atoms.size());
   for (unsigned int i=0; i<unlimited_atoms.size(); i++) { 
      try {
	 lig_build::pos_t p;
      }
      catch (std::runtime_error rte) {
	 std::cout << rte.what() << std::endl;
      } 
   }

   std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > line_fragments;

   lbg_info_t::grid_index_t grid_index_prev(0,0);

   for (int ix=0; ix<x_size(); ix+=1) {
      for (int iy=0; iy<y_size(); iy+=1) {
	 int ms_type = square_type(ix, iy, contour_level);

	 lbg_info_t::grid_index_t grid_index(ix,iy);

	 if ((ms_type != MS_NO_CROSSING) && (ms_type != MS_NO_SQUARE)) { 
	    lbg_info_t::contour_fragment cf(ms_type, contour_level,
					    grid_index_prev,
					    grid_index,
					    *this); // sign of bad architecture
	    if (cf.coords.size() == 1) {

	       if (debug)
		  std::cout << "plot contour ("
			    << cf.get_coords(ix, iy, 0).first << " "
			    << cf.get_coords(ix, iy, 0).second << ") to ("
			    << cf.get_coords(ix, iy, 1).first << " "
			    << cf.get_coords(ix, iy, 1).second << ")" << std::endl;
			
	       std::pair<double, double> xy_1 = cf.get_coords(ix, iy, 0);
	       std::pair<double, double> xy_2 = cf.get_coords(ix, iy, 1);
	       lig_build::pos_t pos_1 = to_canvas_pos(xy_1.first, xy_1.second);
	       lig_build::pos_t pos_2 = to_canvas_pos(xy_2.first, xy_2.second);

	       lig_build::pos_t p1 = to_canvas_pos(cf.get_coords(ix, iy, 0).first,
						   cf.get_coords(ix, iy, 0).second);
	       lig_build::pos_t p2 = to_canvas_pos(cf.get_coords(ix, iy, 1).first,
						   cf.get_coords(ix, iy, 1).second);
	       std::pair<lig_build::pos_t, lig_build::pos_t> fragment_pair(p1, p2);

	       // Now filter out this fragment pair if it is too close
	       // to an unlimited_atom_positions
	       bool plot_it = 1;
	       double dist_crit = 4.0 * LIGAND_TO_CANVAS_SCALE_FACTOR;

	       
	       for (unsigned int i=0; i<unlimited_atoms.size(); i++) { 
// 		  lig_build::pos_t p = to_canvas_pos(unlimited_atom_positions[i].x,
// 						     unlimited_atom_positions[i].y);
		  
		  lig_build::pos_t p = unlimited_atoms[i].atom.atom_position;

		  // if this atom has a ring centre, use the ring
		  // centre to atom vector to unplot vectors only in a
		  // particular direction.
		  // 
		  if ((p - p1).lengthsq() < (dist_crit * dist_crit)) {
		     if (1) { // for debugging
			if (unlimited_atoms[i].has_ring_centre_flag) {
			   // std::cout << " atom " << i << " has ring_centre ";
			   lig_build::pos_t d_1 =
			      unlimited_atoms[i].ring_centre - unlimited_atoms[i].atom.atom_position;
			   lig_build::pos_t d_2 = unlimited_atoms[i].atom.atom_position - p1;
			   double cos_theta =
			      lig_build::pos_t::dot(d_1, d_2)/(d_1.length()*d_2.length());
			   // std::cout << " cos_theta " << cos_theta << " for unlimited atom " << i << std::endl;
			   if (cos_theta > 0.3) { // only cut in the "forwards" direction

// 			      std::cout << " cutting by ring-centred unlimited atom " << i << " "
// 				 << unlimited_atoms[i].atom.get_atom_name()
// 					<< std::endl;
			      
			      plot_it = 0;
			      break;
			   }
			   // std::cout << std::endl;
			
			} else {
			   plot_it = 0;
//  			   std::cout << " cutting by unlimited atom " << i << " "
//  				     << unlimited_atoms[i].atom.get_atom_name()
//  				     << std::endl;
			   break;
			}
		     }
		  }
		  
	       } // end unlimited atoms loop

	       if (plot_it)
		  line_fragments.push_back(fragment_pair);
	    } 
	 }
      }
   }

   std::vector<std::vector<lig_build::pos_t> > contour_lines = make_contour_lines(line_fragments);

   plot_contour_lines(contour_lines, root);

   // check the orientation of the canvas
   if (0) { 
      lig_build::pos_t grid_ori = to_canvas_pos(0.0, 0.0);
      goo_canvas_rect_new (group,
			   grid_ori.x, grid_ori.y, 5.0, 5.0,
			   "line-width", 1.0,
			   "stroke-color", "green",
			   "fill_color", "blue",
			   NULL);
   }
} 

std::vector<std::vector<lig_build::pos_t> >
lbg_info_t::ligand_grid::make_contour_lines(const std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > &line_fragments) const {

   // Look for neighboring cells so that a continuous line fragment is
   // generated. That is fine for closed contours, but this handled
   // badly contour line fragments that go over the edge, they will be
   // split into several line fragment - depending on which fragment
   // is picked first.
   //
   // Perhaps the last items in the line_frag_queue should be (any)
   // points that are on the edge of the grid - because they are the
   // best starting place for dealing with contour line fragments that
   // go over the edge.  If you do this, then the input arg
   // line_fragments should carry with it info about this line
   // fragment being on an edge (generated in show_contour())..
   //
   
   std::vector<std::vector<lig_build::pos_t> > v; // returned item.

   std::deque<std::pair<lig_build::pos_t, lig_build::pos_t> > line_frag_queue;
   for (unsigned int i=0; i<line_fragments.size(); i++)
      line_frag_queue.push_back(line_fragments[i]);

   while (!line_frag_queue.empty()) {

      // start a new line fragment
      // 
      std::vector<lig_build::pos_t> working_points;
      std::pair<lig_build::pos_t, lig_build::pos_t> start_frag = line_frag_queue.back();
      working_points.push_back(start_frag.first);
      working_points.push_back(start_frag.second);
      line_frag_queue.pop_back();
      bool line_fragment_terminated = 0;

      while (! line_fragment_terminated) { 
	 // Is there a fragment in line_frag_queue that has a start that
	 // is working_points.back()?
	 //
	 bool found = 0;
	 std::deque<std::pair<lig_build::pos_t, lig_build::pos_t> >::iterator it;
	 for (it=line_frag_queue.begin(); it!=line_frag_queue.end(); it++) { 
	    if (it->first.near_point(working_points.back(), 0.1)) {
	       working_points.push_back(it->second);
	       line_frag_queue.erase(it);
	       found = 1;
	       break;
	    }
	 }
	 if (! found)
	    line_fragment_terminated = 1;
      }

      v.push_back(working_points);
   }
   
   return v;
} 

void
lbg_info_t::ligand_grid::plot_contour_lines(const std::vector<std::vector<lig_build::pos_t> > &contour_lines, GooCanvasItem *root) const {

   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color",
						// "#111111", not so dark
						"#777777",
						NULL);
   GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 1.5, 2.5);

   for (unsigned int i=0; i<contour_lines.size(); i++) { 
      for (int j=0; j<(contour_lines[i].size()-1); j++) {
	 
	 goo_canvas_polyline_new_line(group,
				      contour_lines[i][j].x,
				      contour_lines[i][j].y,
				      contour_lines[i][j+1].x,
				      contour_lines[i][j+1].y,
				      "line_width", 1.0,
				      "line-dash", dash,
				      NULL);
	 
      }
   }
   goo_canvas_line_dash_unref(dash);
   
} 




void
lbg_info_t::show_mol_ring_centres() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   std::vector<lig_build::pos_t> c = mol.get_ring_centres();
   for (unsigned int i=0; i<c.size(); i++) {
      GooCanvasItem *rect_item = goo_canvas_rect_new (root,
						      c[i].x, c[i].y, 4.0, 4.0,
						      "line-width", 1.0,
						      "stroke-color", "blue",
						      "fill_color", "blue",
						      NULL);
   }
}

void
lbg_info_t::show_unlimited_atoms(const std::vector<widgeted_atom_ring_centre_info_t> &ua) {
   
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   for (unsigned int i=0; i<ua.size(); i++) {
      GooCanvasItem *rect_item =
	 goo_canvas_rect_new(root,
			     ua[i].atom.atom_position.x -6.0,
			     ua[i].atom.atom_position.y -6.0,
			     12.0, 12.0,
			     "line-width", 1.0,
			     "stroke-color", "lightblue",
			     "fill_color", "lightblue",
			     NULL);
   }
}

void
lbg_info_t::show_ring_centres(std::vector<std::vector<std::string> > ring_atoms_list,
			      const widgeted_molecule_t &mol) {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   for (unsigned int iring=0; iring<ring_atoms_list.size(); iring++) {
      try { 
	 lig_build::pos_t ring_centre = mol.get_ring_centre(ring_atoms_list[iring]);
	 GooCanvasItem *rect_item =
	    goo_canvas_rect_new(root,
				ring_centre.x -6.0,
				ring_centre.y -6.0,
				12.0, 12.0,
				"line-width", 1.0,
				"stroke-color", "purple",
				"fill_color", "purple",
				NULL);

      }
      catch (std::runtime_error rte) {
	 std::cout << "Opps - failed to find ring centre for ring atom name "
		   << iring << std::endl;
      }
   }

} 



// convert val [0->255] to a hex colour string in red
// 
std::string
lbg_info_t::grid_intensity_to_colour(int val) const {

   std::string colour = "#000000";
   int step = 8;
      
   // higest value -> ff0000
   // lowest value -> ffffff

   int inv_val = 255 - val;
   if (inv_val < 0)
      inv_val = 0;
   
   int val_high = (inv_val & 240)/16;
   int val_low  = inv_val &  15;

   std::string s1 = sixteen_to_hex_let(val_high);
   std::string s2 = sixteen_to_hex_let(val_low);

   colour = "#ff";
   colour += s1;
   colour += s2;
   colour += s1;
   colour += s2;

   // std::cout << "val " << val << " returns colour " << colour << "\n";
   return colour;
}

std::string
lbg_info_t::sixteen_to_hex_let(int v) const {

   std::string l = "0";
   if (v == 15)
      l = "f";
   if (v == 14)
      l = "e";
   if (v == 13)
      l = "d";
   if (v == 12)
      l = "c";
   if (v == 11)
      l = "b";
   if (v == 10)
      l = "a";
   if (v == 9)
      l = "9";
   if (v == 8)
      l = "8";
   if (v == 7)
      l = "7";
   if (v == 6)
      l = "6";
   if (v == 5)
      l = "5";
   if (v == 4)
      l = "4";
   if (v == 3)
      l = "3";
   if (v == 2)
      l = "2";
   if (v == 1)
      l = "1";

   return l;

}

lbg_info_t::ligand_grid::ligand_grid(const lig_build::pos_t &low_x_and_y,
				     const lig_build::pos_t &high_x_and_y) {

   double extra_extents = 100;
   top_left     = low_x_and_y  - lig_build::pos_t(extra_extents, extra_extents);
   bottom_right = high_x_and_y + lig_build::pos_t(extra_extents, extra_extents);
   scale_fac = 5; // seems good
   double delta_x = bottom_right.x - top_left.x;
   double delta_y = bottom_right.y - top_left.y;
   if (0)
      std::cout << " in making grid, got delta_x and delta_y "
		<< delta_x << " " << delta_y << std::endl;
   x_size_ = int(delta_x/scale_fac+1);
   y_size_ = int(delta_y/scale_fac+1);

   std::vector<double> tmp_y (y_size_, 0.0);
   grid_.resize(x_size_);
   for (unsigned int i=0; i<x_size_; i++)
      grid_[i] = tmp_y;
   if (0) 
      std::cout << "--- ligand_grid constructor, grid has extents "
		<< x_size_ << " " << y_size_ << " real " << grid_.size()
		<< " " << grid_[0].size()
		<< std::endl;
   
}

// arg is not const reference because get_ring_centres() caches the return value.
void
lbg_info_t::ligand_grid::fill(widgeted_molecule_t mol) {

   double exp_scale = 0.0011;
   double rk = 3000.0;

   // int grid_extent = 15; // 10, 12 is not enough
   int grid_extent = 50 ; // untraps 2wot residues?
   
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      for (int ipos_x= -grid_extent; ipos_x<=grid_extent; ipos_x++) {
	 for (int ipos_y= -grid_extent; ipos_y<=grid_extent; ipos_y++) {
	    std::pair<int, int> p = canvas_pos_to_grid_pos(mol.atoms[iat].atom_position);
	    int ix_grid = ipos_x + p.first;
	    int iy_grid = ipos_y + p.second;
	    if ((ix_grid >= 0) && (ix_grid < x_size())) {
	       if ((iy_grid >= 0) && (iy_grid < y_size())) {
		  double d2 = (to_canvas_pos(ix_grid, iy_grid) - mol.atoms[iat].atom_position).lengthsq();
		  double val =  rk * exp(-0.5*exp_scale*d2);
		  grid_[ix_grid][iy_grid] += val;
	       } else {
// 		  std::cout << "ERROR:: out of range in y: " << ix_grid << "," << iy_grid << " "
// 			    << "and grid size: " << x_size() << "," << y_size() << std::endl;
	       }
	    } else {
// 	       std::cout << "ERROR:: out of range in x: " << ix_grid << "," << iy_grid << " "
// 			 << "and grid size: " << x_size() << "," << y_size() << std::endl;
	    }
	 }
      }
   }

   std::vector<lig_build::pos_t> mol_ring_centres = mol.get_ring_centres();

   // std::cout << "DEBUG:: found " << mol_ring_centres.size() << " ring centres " << std::endl;
   
   for (unsigned int ir=0; ir<mol_ring_centres.size(); ir++) { 
      for (int ipos_x= -10; ipos_x<=10; ipos_x++) {
	 for (int ipos_y= -10; ipos_y<=10; ipos_y++) {
	    std::pair<int, int> p = canvas_pos_to_grid_pos(mol_ring_centres[ir]);
	    int ix_grid = ipos_x + p.first;
	    int iy_grid = ipos_y + p.second;
	    if ((ix_grid >= 0) && (ix_grid < x_size())) {
	       if ((iy_grid >= 0) && (iy_grid < y_size())) {
		  double d2 = (to_canvas_pos(ix_grid, iy_grid) - mol_ring_centres[ir]).lengthsq();
		  double val = rk * exp(-0.5* exp_scale * d2);
		  grid_[ix_grid][iy_grid] += val;
	       }
	    }
	 }
      }
   }
   normalize(); // scaled peak value to 1.
}


void
lbg_info_t::ligand_grid::add_for_accessibility(double bash_dist, const lig_build::pos_t &atom_pos) {

   bool debug = 0;
   int grid_extent = 45;

   double inv_scale_factor = 1.0/double(LIGAND_TO_CANVAS_SCALE_FACTOR);
   
   for (int ipos_x= -grid_extent; ipos_x<=grid_extent; ipos_x++) {
      for (int ipos_y= -grid_extent; ipos_y<=grid_extent; ipos_y++) {
	 std::pair<int, int> p = canvas_pos_to_grid_pos(atom_pos);
	 int ix_grid = ipos_x + p.first;
	 int iy_grid = ipos_y + p.second;
	 if ((ix_grid >= 0) && (ix_grid < x_size())) {
	    if ((iy_grid >= 0) && (iy_grid < y_size())) {
	       
	       double d2 = (to_canvas_pos(ix_grid, iy_grid) - atom_pos).lengthsq();
	       d2 *= (inv_scale_factor * inv_scale_factor);
	       double val = substitution_value(d2, bash_dist);
	       if (debug)
		  if (val > 0)
		     std::cout << "adding " << val << " to grid " << ix_grid << " " << iy_grid
			       << " from " << sqrt(d2) << " vs " << bash_dist << std::endl;
	       grid_[ix_grid][iy_grid] += val;
	    }
	 }
      }
   }
}

void
lbg_info_t::ligand_grid::add_for_accessibility_no_bash_dist_atom(double scale,
								 const lig_build::pos_t &atom_pos) { 

   bool debug = 1;
   int grid_extent = 40;

   double inv_scale_factor = 1.0/double(LIGAND_TO_CANVAS_SCALE_FACTOR);
   
   for (int ipos_x= -grid_extent; ipos_x<=grid_extent; ipos_x++) {
      for (int ipos_y= -grid_extent; ipos_y<=grid_extent; ipos_y++) {
	 std::pair<int, int> p = canvas_pos_to_grid_pos(atom_pos);
	 int ix_grid = ipos_x + p.first;
	 int iy_grid = ipos_y + p.second;
	 if ((ix_grid >= 0) && (ix_grid < x_size())) {
	    if ((iy_grid >= 0) && (iy_grid < y_size())) {
	       
	       double d2 = (to_canvas_pos(ix_grid, iy_grid) - atom_pos).lengthsq();
	       d2 *= (inv_scale_factor * inv_scale_factor);
	       // a triangle function, 1 at the atom centre, 0 at 1.5A and beyond
	       //
	       double dist_at_zero = 2.6;
	       
	       double val = 0.0;
	       if (d2< dist_at_zero * dist_at_zero)
		  val = - 1.0/dist_at_zero * sqrt(d2) + 1.0;
	       grid_[ix_grid][iy_grid] += val;
	    }
	 }
      }
   }
}



// scale peak value to 1.0
void
lbg_info_t::ligand_grid::normalize() {

   double max_int = 0.0;

   // std::cout << "normalizing grid " << x_size() << " by " << y_size() << std::endl;
   for (int ix=0; ix<x_size(); ix++) {
      for (int iy=0; iy<y_size(); iy++) {
	 double intensity = grid_[ix][iy];
	 if (intensity > max_int)
	    max_int = intensity;
      }
   }
   if (max_int > 0.0) {
      double sc_fac = 1.0/max_int;
      for (int ix=0; ix<x_size(); ix++) {
	 for (int iy=0; iy<y_size(); iy++) {
	    grid_[ix][iy] *= sc_fac;
	 }
      }
   }
}



// minimise layout energy
std::pair<int, std::vector<lbg_info_t::residue_circle_t> >
lbg_info_t::optimise_residue_circle_positions(const std::vector<lbg_info_t::residue_circle_t> &r,
					      const std::vector<lbg_info_t::residue_circle_t> &c,
					      const std::vector<int> &primary_indices) const { 

   if (r.size() > 0) {
      if (c.size() == r.size()) { 
	 optimise_residue_circles orc(r, c, mol, primary_indices);
	 int status = orc.get_gsl_min_status();
	 return orc.solution();
      }
   }
   std::vector<lbg_info_t::residue_circle_t> dv; // dummy 
   return std::pair<int, std::vector<lbg_info_t::residue_circle_t> > (0, dv);
}


std::vector<lbg_info_t::residue_circle_t>
lbg_info_t::read_residues(const std::string &file_name) const {
      
   std::vector<residue_circle_t> v;
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {

      std::cout << "opened residue_circle file: " << file_name << std::endl;
      
      std::vector<std::string> lines;
      std::string line;
      while (std::getline(f, line)) { 
	 lines.push_back(line);
      }

      for (unsigned int i=0; i<lines.size(); i++) {
	 std::vector<std::string> words = coot::util::split_string_no_blanks(lines[i], " ");
	 
	 // debug input
	 if (0) { 
	    std::cout << i << " " << lines[i] << "\n";
	    for (unsigned int j=0; j<words.size(); j++) { 
	       std::cout << "  " << j << " " << words[j] << " ";
	    }
	    std::cout << "\n";
	 }

	 if (words.size() > 5) {
	    if (words[0] == "RES") {
	       try { 
		  double pos_x = lig_build::string_to_float(words[1]);
		  double pos_y = lig_build::string_to_float(words[2]);
		  double pos_z = lig_build::string_to_float(words[3]);
		  std::string res_type = words[4];
		  std::string label = words[5];
		  coot::residue_spec_t dum_spec;
		  residue_circle_t rc(pos_x, pos_y, pos_z, dum_spec, res_type, label);
		  clipper::Coord_orth cp(pos_x, pos_y, pos_z);
		  lig_build::pos_t pos = mol.input_coords_to_canvas_coords(cp);
		  if (res_type == "HOH") {
		     if (words.size() > 7) {
			try {
			   double dist_to_protein = lig_build::string_to_float(words[7]);
			   rc.set_water_dist_to_protein(dist_to_protein);
			}
			catch (std::runtime_error rte) {
			}
		     } 
		  } 
		  rc.set_canvas_pos(pos);
		  v.push_back(residue_circle_t(rc)); // why is there a constructor here?
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }
	 if (words.size() == 6) {
	    if (words[0] == "BOND") { // written with space
	       try {
		  double bond_l = lig_build::string_to_float(words[3]);
		  std::string atom_name = lines[i].substr(5,4);
		  int bond_type = lig_build::string_to_int(words[5]);
		  lbg_info_t::bond_to_ligand_t btl(atom_name, bond_l);
		  btl.bond_type = bond_type;
// 		  std::cout << "adding bond " << v.size() << " to :"
// 			    << atom_name << ": " << bond_l << std::endl;
		  if (v.size())
		     v.back().add_bond_to_ligand(btl);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }

	 if (words.size() == 3) {
	    if (words[0] == "SOLVENT_ACCESSIBILITY_DIFFERENCE") {
	       try {
		  double se_holo = lig_build::string_to_float(words[1]);
		  double se_apo  = lig_build::string_to_float(words[2]);
		  if (v.size()) { 
		     v.back().set_solvent_exposure_diff(se_holo, se_apo);
		  }
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }

	 if (words.size() > 2) {
	    if (words[0] == "STACKING") {
	       std::vector<std::string> ligand_ring_atom_names;
	       std::string type = words[1];
	       std::string::size_type l = lines[i].length();
	       int n_atoms_max = 20;
	       for (unsigned int iat=0; iat<n_atoms_max; iat++) {
		  int name_index = 19+iat*6;
		  if (l > (name_index+4)) {
		     std::string ss = lines[i].substr(name_index, 4);
		     if (0)
			std::cout << "   Found ligand ring stacking atom :"
				  << ss << ":" << std::endl;
		     ligand_ring_atom_names.push_back(ss);
		  }
	       }
	       if (type == "pi-pi")
		  if (v.size())
		     v.back().set_stacking(type, ligand_ring_atom_names, "");
	       if (type == "pi-cation")
		  if (v.size())
		     v.back().set_stacking(type, ligand_ring_atom_names, "");
	       if (type == "cation-pi") {
		  std::string ligand_cation_atom_name = lines[i].substr(19,4);
		  // std::cout << "DEBUG:: on reading residue info file found ligand cation "
		  // << "name :" << ligand_cation_atom_name << ":" << std::endl;
		  v.back().set_stacking(type, ligand_ring_atom_names, ligand_cation_atom_name);
	       }
	    }
	 }
      }
   }
   std::cout << "found " << v.size() << " residue centres" << std::endl;
   return v;
} 

// "must take exactly one argument" problem
// 
// std::ostream &
// lbg_info_t::operator<<(std::ostream &s, residue_circle_t rc) {

//    s << "res-circ{" << rc.pos << " " << rc.label << " with bond_to_ligand length "
//      << bond_to_ligand.bond_length << "}";

// }


// if you don't have add_rep_handles, then pass a vector or size 0.
// 
void
lbg_info_t::draw_residue_circles(const std::vector<residue_circle_t> &l_residue_circles,
				 const std::vector<int> &add_rep_handles) {

   double max_dist_water_to_ligand_atom  = 3.3; // don't draw waters that are far from ligand
   double max_dist_water_to_protein_atom = 3.3; // don't draw waters that are not somehow 
                                                // attached to the protein.

   if (draw_flev_annotations_flag) { 
      bool draw_solvent_exposures = 1;
      try { 
	 lig_build::pos_t ligand_centre = mol.get_ligand_centre();
	 GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

	 if (draw_solvent_exposures)
	    for (unsigned int i=0; i<l_residue_circles.size(); i++)
	       draw_solvent_exposure_circle(l_residue_circles[i], ligand_centre, root);

	 for (unsigned int i=0; i<l_residue_circles.size(); i++) {
	    lig_build::pos_t pos = l_residue_circles[i].pos;
	    int add_rep_handle = -1; // default, no handle
	    if (add_rep_handles.size() == l_residue_circles.size())
	       add_rep_handle = add_rep_handles[i];

	    draw_residue_circle_top_layer(l_residue_circles[i], ligand_centre, add_rep_handle);
	 }
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING:: draw_residue_circles: " << rte.what() << std::endl;
      }
   }
}

// static
/* This handles button presses on the canvas */
// extern "C" G_MODULE_EXPORT 
static bool
on_residue_circle_clicked(GooCanvasItem  *item,
			  GooCanvasItem  *target_item,
			  GdkEventButton *event,
			  gpointer        user_data) { 

   gboolean r = true;

   std::cout << "on_residue_circle_clicked()!!!!! " << std::endl;

   if (target_item) { 
      coot::residue_spec_t *spec_p = (coot::residue_spec_t *) g_object_get_data (G_OBJECT (target_item), "spec");
   
      if (spec_p) { 
      
	 std::cout << "clicked on " << *spec_p << std::endl;
      
      } else {
	 std::cout << "null spec" << std::endl;
      }
   } else {
      std::cout << "NULL target item" << std::endl;
   } 

   return r;
}


// if you don't have an add_rep_handle, then pass -1 (something negative)
// 
void
lbg_info_t::draw_residue_circle_top_layer(const residue_circle_t &residue_circle,
					  const lig_build::pos_t &ligand_centre,
					  int add_rep_handle) {


   if (0)
      std::cout << "   adding cirles " << residue_circle.residue_type
		<< " at init pos " << residue_circle.pos << " and canvas_drag_offset "
		<< canvas_drag_offset << std::endl;

   lig_build::pos_t circle_pos = residue_circle.pos;
      
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", "#111111",
						NULL);

   // don't draw waters that don't have bonds to the ligand (or the
   // bonds to the ligand atoms are too long, or the water is too far
   // from any protein atom).
   //
   if (residue_circle.residue_type == "HOH") {
      if (residue_circle.bonds_to_ligand.size() == 0) {
	 return;
      }
   }

   // Capitalise the residue type (takes less space than upper case).
   std::string rt = residue_circle.residue_type.substr(0,1);
   rt += coot::util::downcase(residue_circle.residue_type.substr(1));

   // fill colour and stroke colour of the residue circle
   std::pair<std::string, std::string> col = get_residue_circle_colour(residue_circle.residue_type);
   double line_width = 1.0;
   if (col.second != "#111111") // needs checking, FIXME
      line_width = 3.0;
   GooCanvasItem *circle = NULL;
   GooCanvasItem *text_1 = NULL;
   GooCanvasItem *text_2 = NULL;
   

   if (col.first != "") {
      circle = goo_canvas_ellipse_new(group,
				      circle_pos.x, circle_pos.y,
				      standard_residue_circle_radius,
				      standard_residue_circle_radius,
				      "line_width", line_width,
				      "fill-color",   col.first.c_str(),
				      "stroke-color", col.second.c_str(),
				      NULL);
   } else {
      circle = goo_canvas_ellipse_new(group,
				      circle_pos.x, circle_pos.y,
				      standard_residue_circle_radius,
				      standard_residue_circle_radius,
				      "line_width", line_width,
				      "stroke-color", col.second.c_str(),
				      NULL);
   }
   
   text_1 = goo_canvas_text_new(group, rt.c_str(),
				circle_pos.x, circle_pos.y-6, -1,
				GTK_ANCHOR_CENTER,
				"font", "Sans 9",
				"fill_color", "#111111",
				NULL);

   text_2 = goo_canvas_text_new(group, residue_circle.residue_label.c_str(),
				circle_pos.x, circle_pos.y+6.5, -1,
				GTK_ANCHOR_CENTER,
				"font", "Sans 7",
				"fill_color", "#111111",
				NULL);

   if (circle && text_1 && text_2) {
      // so that if a residue is clicked on (or moused over or
      // whatever, then the target_item has some interesting attached
      // data for use in the callback)
      coot::residue_spec_t *sp_p_1 = new coot::residue_spec_t(residue_circle.spec);
      coot::residue_spec_t *sp_p_2 = new coot::residue_spec_t(residue_circle.spec);
      coot::residue_spec_t *sp_p_3 = new coot::residue_spec_t(residue_circle.spec);
      sp_p_1->int_user_data = imol;
      sp_p_2->int_user_data = imol;
      sp_p_3->int_user_data = imol;
      g_object_set_data_full(G_OBJECT(circle), "spec", sp_p_1, g_free);
      g_object_set_data_full(G_OBJECT(text_1), "spec", sp_p_2, g_free);
      g_object_set_data_full(G_OBJECT(text_2), "spec", sp_p_3, g_free);
      int *add_rep_handle_p_1 = new int(add_rep_handle);
      int *add_rep_handle_p_2 = new int(add_rep_handle);
      int *add_rep_handle_p_3 = new int(add_rep_handle);
      g_object_set_data_full(G_OBJECT(circle), "add_rep_handle", add_rep_handle_p_1, g_free);
      g_object_set_data_full(G_OBJECT(text_1), "add_rep_handle", add_rep_handle_p_2, g_free);
      g_object_set_data_full(G_OBJECT(text_2), "add_rep_handle", add_rep_handle_p_3, g_free);
   }
}


// Return the fill colour and the stroke colour.
// 
std::pair<std::string, std::string>
lbg_info_t::get_residue_circle_colour(const std::string &residue_type) const {

   std::string fill_colour = "";
   std::string stroke_colour = "#111111";

   std::string grease = "#ccffbb";
   std::string purple = "#eeccee";
   std::string red    = "#cc0000";
   std::string blue   = "#0000cc";
   std::string metalic_grey = "#d9d9d9";

   if (residue_type == "ALA")
      fill_colour = grease;
   if (residue_type == "TRP")
      fill_colour = grease;
   if (residue_type == "PHE")
      fill_colour = grease;
   if (residue_type == "LEU")
      fill_colour = grease;
   if (residue_type == "PRO")
      fill_colour = grease;
   if (residue_type == "ILE")
      fill_colour = grease;
   if (residue_type == "VAL")
      fill_colour = grease;
   if (residue_type == "MET")
      fill_colour = grease;
   if (residue_type == "MSE")
      fill_colour = grease;

   if (residue_type == "GLY")
      fill_colour = purple;
   if (residue_type == "ASP")
      fill_colour = purple;
   if (residue_type == "ASN")
      fill_colour = purple;
   if (residue_type == "CYS")
      fill_colour = purple;
   if (residue_type == "GLN")
      fill_colour = purple;
   if (residue_type == "GLU")
      fill_colour = purple;
   if (residue_type == "HIS")
      fill_colour = purple;
   if (residue_type == "LYS")
      fill_colour = purple;
   if (residue_type == "LYS")
      fill_colour = purple;
   if (residue_type == "ARG")
      fill_colour = purple;
   if (residue_type == "SER")
      fill_colour = purple;
   if (residue_type == "THR")
      fill_colour = purple;
   if (residue_type == "TYR")
      fill_colour = purple;
   if (residue_type == "HOH")
      fill_colour = "white";

   if (residue_type == "ASP")
      stroke_colour = red;
   if (residue_type == "GLU")
      stroke_colour = red;
   if (residue_type == "LYS")
      stroke_colour = blue;
   if (residue_type == "ARG")
      stroke_colour = blue;
   if (residue_type == "HIS")
      stroke_colour = blue;

   // metals
   if (residue_type == "ZN") 
      fill_colour = metalic_grey;
   if (residue_type == "MG")
      fill_colour = metalic_grey;
   if (residue_type == "NA")
      fill_colour = metalic_grey;
   if (residue_type == "CA")
      fill_colour = metalic_grey;
   if (residue_type == "K")
      fill_colour = metalic_grey;
	 
   return std::pair<std::string, std::string> (fill_colour, stroke_colour);
}




// solvent exposure difference of the residue due to ligand binding
void 
lbg_info_t::draw_solvent_exposure_circle(const residue_circle_t &residue_circle,
					 const lig_build::pos_t &ligand_centre,
					 GooCanvasItem *group) {

   if (residue_circle.residue_type != "HOH") { 
      if (residue_circle.se_diff_set()) {
	 std::pair<double, double> se_pair = residue_circle.solvent_exposures();
	 double radius_extra = (se_pair.second - se_pair.first) * 19;  // was 18, was 14, was 22.
	 if (radius_extra > 0) {
	    lig_build::pos_t to_lig_centre_uv = (ligand_centre - residue_circle.pos).unit_vector();
	    lig_build::pos_t se_circle_centre = residue_circle.pos - to_lig_centre_uv * radius_extra;

	    std::string fill_colour = get_residue_solvent_exposure_fill_colour(radius_extra);
	    double r = standard_residue_circle_radius + radius_extra;

	    GooCanvasItem *circle = goo_canvas_ellipse_new(group,
							   se_circle_centre.x, se_circle_centre.y,
							   r, r,
							   "line_width", 0.0,
							   "fill-color", fill_colour.c_str(),
							   NULL);
	 }
      }
   } 
}

std::string
lbg_info_t::get_residue_solvent_exposure_fill_colour(double r) const {

   std::string colour = "#8080ff";
   double step = 0.7;
   if (r > step)
      colour = "#e0e0ff";
   if (r > step * 2)
      colour = "#d8d8ff";
   if (r > step * 3)
      colour = "#d0d0ff";
   if (r > step * 4)
      colour = "#c0c8ff";
   if (r > step * 5)
      colour = "#b0c0ff";
   if (r > step * 6)
      colour = "#a0b8ff";
   if (r > step * 7)
      colour = "#90b0ff";
   if (r > step * 8)
      colour = "#80a8ff";
   if (r > step * 9)
      colour = "#70a0ff";

   return colour;
} 


std::vector<solvent_accessible_atom_t>
lbg_info_t::read_solvent_accessibilities(const std::string &file_name) const {

   // return this
   std::vector<solvent_accessible_atom_t> solvent_accessible_atoms;
   
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {

      std::cout << "reading solvent accessibilites file: " << file_name << std::endl;
      
      std::vector<std::string> lines;
      std::string line;
      while (std::getline(f, line)) { 
	 lines.push_back(line);
      }

      for (unsigned int i=0; i<lines.size(); i++) {
	 std::vector<std::string> words = coot::util::split_string_no_blanks(lines[i], " ");
	 if (words.size() > 5) {
	    
	    if (words[0] == "ATOM:") {
	       std::string atom_name = lines[i].substr(5,4);
	       try {
		  double pos_x = lig_build::string_to_float(words[2]);
		  double pos_y = lig_build::string_to_float(words[3]);
		  double pos_z = lig_build::string_to_float(words[4]);
		  double sa    = lig_build::string_to_float(words[5]);
		  clipper::Coord_orth pt(pos_x, pos_y, pos_z);
		  
		  if (0) 
		     std::cout << "got atom name :" << atom_name << ": and pos "
			       << pt.format() << " and accessibility: " << sa
			       << std::endl;
		  solvent_accessible_atom_t saa(atom_name, pt, sa);
		  solvent_accessible_atoms.push_back(saa);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }

	 if (words[0] == "BASH:") {
	    if (words.size() == 2) {
	       try {
		  if (words[1] == "unlimited") {
		     if (solvent_accessible_atoms.size())
			solvent_accessible_atoms.back().add_unlimited();
		  } else {
		     if (solvent_accessible_atoms.size()) { 
			double bash_dist = lig_build::string_to_float(words[1]);
			solvent_accessible_atoms.back().add_bash_dist(bash_dist);
		     }
		  }
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    } 
	 }
      }
   }
   return solvent_accessible_atoms;
}

void
lbg_info_t::draw_solvent_accessibility_of_atoms() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      lig_build::pos_t pos = mol.atoms[iat].atom_position;
      double sa = mol.atoms[iat].get_solvent_accessibility();
      if (sa  > 0) 
	 draw_solvent_accessibility_of_atom(pos, sa, root);
   }
}

void
lbg_info_t::draw_solvent_accessibility_of_atom(const lig_build::pos_t &pos, double sa,
					       GooCanvasItem *root) {

   if (draw_flev_annotations_flag) { 
      int n_circles = int(sa*40) + 1;    // needs fiddling?
      if (n_circles> 10) n_circles = 10; // needs fiddling?
      
      for (unsigned int i=0; i<n_circles; i++) { 
	 double rad =  LIGAND_TO_CANVAS_SCALE_FACTOR/23.0 * 3.0 * double(i+1); // needs fiddling?
	 GooCanvasItem *circle = goo_canvas_ellipse_new(root,
							pos.x, pos.y,
							rad, rad,
							"line_width", 0.0,
							"fill-color-rgba", 0x5555cc30,
							NULL);
      }
   }
}

void
lbg_info_t::draw_substitution_contour() {

   bool debug = 0;

   if (draw_flev_annotations_flag) { 
      if (mol.atoms.size() > 0) {


	 // first of all, do we have any bash distances for the atoms of this molecule?
	 bool have_bash_distances = 0;
	 for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
	    if (mol.atoms[iat].bash_distances.size()) {
	       have_bash_distances = 1;
	       break;
	    }
	 }

	 // If we don't have bash distances, then don't grid and contour anything.  If
	 // we do, we do....
	 // 
	 if (have_bash_distances) {
      
	    GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
	    try { 

	       std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair =
		  mol.ligand_extents();
	       ligand_grid grid(l_e_pair.first, l_e_pair.second);
   
	       if (0) { // debug
		  for (unsigned int i=0; i<mol.atoms.size(); i++) { 
		     std::cout << "in draw_substitution_contour() atom " << i << " "
			       << mol.atoms[i].get_atom_name()
			       << " has "
			       << mol.atoms[i].bash_distances.size() << " bash distances"
			       << std::endl;
		     for (unsigned int j=0; j<mol.atoms[i].bash_distances.size(); j++) {
			std::cout << "  " << mol.atoms[i].bash_distances[j];
		     }
		     if (mol.atoms[i].bash_distances.size())
			std::cout << std::endl;
		  }
	       }

	       std::vector<widgeted_atom_ring_centre_info_t> unlimited_atoms;
	       for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
		  int n_bash_distances = 0;
		  double sum_bash = 0.0;
		  bool unlimited = 0;
		  int n_unlimited = 0;
	       
		  if (mol.atoms[iat].bash_distances.size()) { 
		     for (unsigned int j=0; j<mol.atoms[iat].bash_distances.size(); j++) {
			if (! mol.atoms[iat].bash_distances[j].unlimited()) {
			   sum_bash += mol.atoms[iat].bash_distances[j].dist;
			   n_bash_distances++;
			} else {
			   unlimited = 1;
			   n_unlimited++;
			} 
		     }

		     if (! unlimited) { 
			if (n_bash_distances > 0) {
			   double bash_av = sum_bash/double(n_bash_distances);
			   if (debug)
			      std::cout << "   none unlimited, using bash_av " << bash_av
					<< " for atom " << mol.atoms[iat].get_atom_name()
					<< std::endl;
			   grid.add_for_accessibility(bash_av, mol.atoms[iat].atom_position);
			}
		     } else {

			// Now, were more than half of the bash distances
			// unlimited?  If yes, then this is unlimited.

			if (double(n_unlimited)/double(mol.atoms[iat].bash_distances.size()) > 0.5) { 
	    
			   // just shove some value in to make the grid values vaguely
			   // correct - the unlimited atoms properly assert themselves
			   // in the drawing of the contour (that is, the selection of
			   // the contour fragments).
			   //
			   grid.add_for_accessibility(1.2, mol.atoms[iat].atom_position);
			   // 	       std::cout << "adding unlimited_atom_position "
			   // 			 << iat << " "
			   // 			 << mol.atoms[iat].atom_position

			   // not elegant because no constructor for
			   // widgeted_atom_ring_centre_info_t (because no simple
			   // constructor for lig_build::atom_t).
			   // 
			   widgeted_atom_ring_centre_info_t ua(mol.atoms[iat]);
			   unlimited_atoms.push_back(ua);
			} else {

			   // treat as a limited:
			   double bash_av =
			      (sum_bash + 4.0 * n_unlimited)/double(n_bash_distances+n_unlimited);
			   if (debug)
			      std::cout << "   few unlimited, as limited using bash_av "
					<< bash_av << " for atom "
					<< mol.atoms[iat].get_atom_name()
					<< std::endl;
			
			   grid.add_for_accessibility(bash_av, mol.atoms[iat].atom_position);
			} 
		     }
		  
		  } else {

		     // we don't get here currently, now that there is an
		     // outer test for having bash distances. Current way
		     // is OK, I think.
		  
		     // there were no bash distances - what do we do?  Leaving out
		     // atoms means gaps over the ligand - bleugh.  Shove some
		     // value in?  1.0?  if they are not hydrogens, of course.
		     // 
		     if (mol.atoms[iat].element != "H") // checked.
			grid.add_for_accessibility_no_bash_dist_atom(1.0, mol.atoms[iat].atom_position);
		  } 
	       }

	       // Put some values around the ring centres too.
	       // 
	       grid.avoid_ring_centres(ring_atoms_list, mol);

	       // show_grid(grid);

	       // std::vector<widgeted_atom_ring_centre_info_t> dummy_unlimited_atoms;
	       grid.show_contour(root, 0.5, unlimited_atoms, ring_atoms_list);
	       // debug
	       // show_unlimited_atoms(unlimited_atoms);
	       // show_ring_centres(ring_atoms_list, mol);

	       if (debug) { 
		  std::cout << "Here are the "<< unlimited_atoms.size()
			    << " unlimited atoms: " << std::endl;
		  for (unsigned int iat=0; iat<unlimited_atoms.size(); iat++)
		     std::cout << "   " << unlimited_atoms[iat] << std::endl;
	       }

	    }
	    catch (std::runtime_error rte) {
	       std::cout << rte.what() << std::endl;
	    }
	 }
      }
   }
}


// r_squared the squared distance between the atom and the grid point.
//
double
lbg_info_t::ligand_grid::substitution_value(double r_squared, double bash_dist) const {

   double D = bash_dist;
   double r = sqrt(r_squared);

   if (bash_dist < 1) {
      // this is not in C&L, so this is an enhancement of their
      // algorithm.  Without this there is non-smoothness as the
      // contour follows the interface between 0 and 1.

      // So we add a nice smooth slope (similar to that of
      // conventional bash distances (below).  However we add a slope
      // between +/- 0.2 of the bash distance.
      //
      double small = 0.2;
      if (r > (D + small)) { 
	 return 0;
      } else {
	 if (r < (D - small)) {
	    return 1;
	 } else {
	    double m = (r-(D-small))/(2*small);
	    return (0.5 * (1 + cos(M_PI * m)));
	 } 
      }

   } else { 

      if (r<1)
	 return 1; 
      if (r<(D-1)) {
	 return 1;
      } else {
	 if (r > D) { 
	    return 0;
	 } else {
	    // double v = 0.5*(1 + cos(0.5 *M_PI * (D-r))); // C&L - eh? typo/bug
	    double v = 0.5 * (1.0 + cos(M_PI * (D-1-r)));
	    return v;
	 }
      }
   }
}




void
lbg_info_t::draw_stacking_interactions(const std::vector<residue_circle_t> &rc) {


   for (unsigned int ires=0; ires<rc.size(); ires++) {
      int st = rc[ires].get_stacking_type();
      if (rc[ires].has_ring_stacking_interaction()) {
	 std::vector<std::string> ligand_ring_atom_names =
	    rc[ires].get_ligand_ring_atom_names();
	 if ((st == lbg_info_t::residue_circle_t::PI_PI_STACKING) ||
	     (st == lbg_info_t::residue_circle_t::PI_CATION_STACKING)) {
	     
	    try {
	       lig_build::pos_t lc = mol.get_ring_centre(ligand_ring_atom_names);
	       draw_annotated_stacking_line(lc, rc[ires].pos, st);
	    }
	    catch (std::runtime_error rte) {
	       std::cout << rte.what() << std::endl;
	    }
	 }
      }
      if (st == lbg_info_t::residue_circle_t::CATION_PI_STACKING) {
	 std::string at_name = rc[ires].get_ligand_cation_atom_name();
	 lig_build::pos_t atom_pos = mol.get_atom_canvas_position(at_name);
	 draw_annotated_stacking_line(atom_pos, rc[ires].pos, st);
      }
   }
}

void
lbg_info_t::draw_annotated_stacking_line(const lig_build::pos_t &ligand_ring_centre,
					 const lig_build::pos_t &residue_pos,
					 int stacking_type) {	

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   double ligand_target_shortening_factor = 3;
   if (stacking_type == lbg_info_t::residue_circle_t::CATION_PI_STACKING)
      ligand_target_shortening_factor = 6;
   lig_build::pos_t a_to_b_uv = (ligand_ring_centre - residue_pos).unit_vector();
   lig_build::pos_t A = residue_pos + a_to_b_uv * 20;
   lig_build::pos_t B = ligand_ring_centre;
   // just short of the middle of the ring
   lig_build::pos_t C = B - a_to_b_uv * ligand_target_shortening_factor; 
   lig_build::pos_t mid_pt = (A + B) * 0.5;
   lig_build::pos_t close_mid_pt_1 = mid_pt - a_to_b_uv * 17;
   lig_build::pos_t close_mid_pt_2 = mid_pt + a_to_b_uv * 17;

   gboolean start_arrow = 0;
   gboolean end_arrow = 1;
   std::string stroke_colour = "#008000";
   GooCanvasItem *group = goo_canvas_group_new (root,
						"stroke-color", stroke_colour.c_str(),
						NULL);

   std::vector<lig_build::pos_t> hex_and_ring_centre(2);
   hex_and_ring_centre[0] = mid_pt - a_to_b_uv * 7;
   hex_and_ring_centre[1] = mid_pt + a_to_b_uv * 7;

   // for cations, we draw a 6-member ring and a "+" text.
   //
   // for pi-pi, we draw 2 6-member rings
   // 
   for(int ir=0; ir<2; ir++) {

      bool do_ring = 0;
      bool do_anion = 0;
      
      if (stacking_type == lbg_info_t::residue_circle_t::PI_PI_STACKING) {
	 do_ring = 1; 
      } 
      if (stacking_type == lbg_info_t::residue_circle_t::PI_CATION_STACKING) {
	 if (ir == 0) {
	    do_ring = 0;
	    do_anion = 1;
	 } else {
	    do_ring = 1;
	    do_anion = 0;
	 }
      }
      if (stacking_type == lbg_info_t::residue_circle_t::CATION_PI_STACKING) {
	 if (ir == 0) {
	    do_ring = 1;
	    do_anion = 0;
	 } else {
	    do_ring = 0;
	    do_anion = 1;
	 }
      }
      
      if (do_ring) { 
	 double angle_step = 60;
	 double r = 8; // radius
	 for (int ipt=1; ipt<=6; ipt++) {
	    int ipt_0 = ipt - 1;
	    double theta_deg_1 = 30 + angle_step * ipt;
	    double theta_1 = theta_deg_1 * DEG_TO_RAD;
	    lig_build::pos_t pt_1 =
	       hex_and_ring_centre[ir] + lig_build::pos_t(sin(theta_1), cos(theta_1)) * r;
      
	    double theta_deg_0 = 30 + angle_step * ipt_0;
	    double theta_0 = theta_deg_0 * DEG_TO_RAD;
	    lig_build::pos_t pt_0 =
	       hex_and_ring_centre[ir] + lig_build::pos_t(sin(theta_0), cos(theta_0)) * r;
	    goo_canvas_polyline_new_line(group,
					 pt_1.x, pt_1.y,
					 pt_0.x, pt_0.y,
					 "line_width", 1.8,
					 NULL);
	 }
	 // Now the circle in the annotation aromatic ring:
	 goo_canvas_ellipse_new(group,
				hex_and_ring_centre[ir].x,
				hex_and_ring_centre[ir].y,
				4.0, 4.0,
				"line_width", 1.0,
				NULL);
      }
      
      if (do_anion) {
	 // the "+" symbol for the anion
	 // 
	 GooCanvasItem *text_1 = goo_canvas_text_new(group,
						     "+",
						     hex_and_ring_centre[ir].x,
						     hex_and_ring_centre[ir].y,
						     -1,
						     GTK_ANCHOR_CENTER,
						     "font", "Sans 12",
						     "fill_color", stroke_colour.c_str(),
						     NULL);
      }
   }


   GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 2.5, 2.5);
   GooCanvasItem *item_1 =
      goo_canvas_polyline_new_line(group,
				   A.x, A.y,
				   close_mid_pt_1.x, close_mid_pt_1.y,
				   "line-width", 2.5,
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   NULL);

   GooCanvasItem *item_2 =
      goo_canvas_polyline_new_line(group,
				   close_mid_pt_2.x, close_mid_pt_2.y,
				   C.x, C.y,
				   "line-width", 2.5,
				   "line-dash", dash,
				   // "end_arrow",   end_arrow,
				   "stroke-color", stroke_colour.c_str(),
				   NULL);

   // Now the circle blob at the centre of the aromatic ligand ring:
   if (stacking_type != lbg_info_t::residue_circle_t::CATION_PI_STACKING)
      goo_canvas_ellipse_new(group,
			     B.x, B.y,
			     3.0, 3.0,
			     "line_width", 1.0,
			     "fill_color", stroke_colour.c_str(),
			     NULL);

}

void
lbg_info_t::draw_bonds_to_ligand() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   GooCanvasLineDash *dash_dash     = goo_canvas_line_dash_new (2, 2.1, 2.1);
   GooCanvasLineDash *dash_unbroken = goo_canvas_line_dash_new (1, 3);

   GooCanvasLineDash *dash = dash_dash; // re-assigned later (hopefully)
   for (unsigned int ic=0; ic<residue_circles.size(); ic++) { 
      if (residue_circles[ic].bonds_to_ligand.size()) {
	 for (unsigned int ib=0; ib<residue_circles[ic].bonds_to_ligand.size(); ib++) { 

	    try {

	       lig_build::pos_t pos = residue_circles[ic].pos;
	       std::string at_name = residue_circles[ic].bonds_to_ligand[ib].ligand_atom_name;
	       lig_build::pos_t lig_atom_pos = mol.get_atom_canvas_position(at_name);
	       lig_build::pos_t rc_to_lig_atom = lig_atom_pos - pos;
	       lig_build::pos_t rc_to_lig_atom_uv = rc_to_lig_atom.unit_vector();
	       lig_build::pos_t B = lig_atom_pos;

	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type != coot::fle_ligand_bond_t::BOND_COVALENT)
		  B -= rc_to_lig_atom_uv * 8; // put the end point (say) of a hydrogen bond
		                              // some pixels away from the atom centre - for better
              		                      // aesthetics.
		  
	       lig_build::pos_t A = pos + rc_to_lig_atom_uv * 20;

	       // some colours
	       std::string blue = "blue";
	       std::string green = "darkgreen";
	       std::string lime = "#888820";
	       std::string stroke_colour = "#111111"; // unset

	       // arrows (acceptor/donor) and stroke colour (depending
	       // on mainchain or sidechain interaction)
	       // 
	       gboolean start_arrow = 0;
	       gboolean   end_arrow = 0;
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_DONOR_SIDECHAIN) {
		  end_arrow = 1;
		  stroke_colour = green;
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_DONOR_MAINCHAIN) {
		  end_arrow = 1;
		  stroke_colour = blue;
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_SIDECHAIN) {
		  start_arrow = 1;
		  stroke_colour = green;
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_MAINCHAIN) {
		  start_arrow = 1;
		  stroke_colour = blue;
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::METAL_CONTACT_BOND) {
		  stroke_colour = "#990099";
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::BOND_COVALENT) {
		  dash = goo_canvas_line_dash_new (2, 2.7, 0.1);
		  stroke_colour = "#bb00bb";
	       }
	       
	       if (residue_circles[ic].residue_type == "HOH") { 
		  stroke_colour = lime;
		  start_arrow = 0;
		  end_arrow = 0;
		  dash = dash_dash;
	       }
	       GooCanvasItem *item = goo_canvas_polyline_new_line(root,
								  A.x, A.y,
								  B.x, B.y,
								  "line-width", 2.5,
								  "line-dash", dash,
 								  "start_arrow", start_arrow,
 								  "end_arrow",   end_arrow,
								  "stroke-color", stroke_colour.c_str(),
								  NULL);
	    }
	    catch (std::runtime_error rte) {
	       std::cout << "WARNING:: " << rte.what() << std::endl;
	    }
	 }
	 
      } else {
	 if (0) 
	    std::cout << "... no bond to ligand from residue circle "
		      << residue_circles[ic].residue_label << " "
		      << residue_circles[ic].residue_type << std::endl;
      }
   }
}

std::string
lbg_info_t::file_to_string(const std::string &file_name) const {

   std::string s;
   std::string line;
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {
      while (std::getline(f, line)) { 
	 s += line;
      }
   }
   return s;
}

void
lbg_info_t::show_key() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   // class variable.
   if (! key_group)
      key_group = goo_canvas_group_new (root, NULL);

   GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 2.5, 2.5);

   double dy  =  30.0;
   double dxc = 250.0; // column offset
   lig_build::pos_t co(dxc, 0);
   lig_build::pos_t ro(0, dy);
   lig_build::pos_t tl(100, 500); // top left of key

   // Lines
   lig_build::pos_t A       = tl;
   lig_build::pos_t B       = tl + lig_build::pos_t(40,  0);
   lig_build::pos_t AB_txt  = tl + lig_build::pos_t(50,  0);
   lig_build::pos_t Ad      = A +  ro;
   lig_build::pos_t Bd      = B +  ro;
   lig_build::pos_t ABd_txt = AB_txt + ro;
   lig_build::pos_t C       = tl + lig_build::pos_t(0,  2*dy);
   lig_build::pos_t D       = tl + lig_build::pos_t(40, 2*dy);
   lig_build::pos_t CD_txt  = tl + lig_build::pos_t(50, 2*dy);
   lig_build::pos_t Cd      = C +  lig_build::pos_t(0,  dy);
   lig_build::pos_t Dd      = D +  lig_build::pos_t(0,  dy);
   lig_build::pos_t CDd_txt = CD_txt + lig_build::pos_t(0,  dy);

   lig_build::pos_t WatA    = tl + ro * 4;
   lig_build::pos_t WatB    = WatA + lig_build::pos_t(40,  0);
   lig_build::pos_t Wat_txt = WatA + lig_build::pos_t(50,  0);
   
   lig_build::pos_t MetalA    = tl + ro * 5;
   lig_build::pos_t MetalB    = MetalA + lig_build::pos_t(40,  0);
   lig_build::pos_t Metal_txt = MetalA + lig_build::pos_t(50,  0);
   // Circles
   lig_build::pos_t E       = tl + co;
   lig_build::pos_t F       = E + ro;
   lig_build::pos_t G       = E + ro * 2;
   lig_build::pos_t H       = E + ro * 3;
   lig_build::pos_t I       = E + ro * 4;
   lig_build::pos_t J       = E + ro * 5;
   lig_build::pos_t Et      = E + lig_build::pos_t(25, 0);
   lig_build::pos_t Ft      = F + lig_build::pos_t(25, 0);
   lig_build::pos_t Gt      = G + lig_build::pos_t(25, 0);
   lig_build::pos_t Ht      = H + lig_build::pos_t(25, 0);
   lig_build::pos_t It      = I + lig_build::pos_t(25, 0);
   lig_build::pos_t Jt      = J + lig_build::pos_t(25, 0);
   // Rest
   lig_build::pos_t K       = tl + co + lig_build::pos_t(150, 0);
   lig_build::pos_t L       = K + ro * 1.5;
   lig_build::pos_t M       = K + ro * 3;
   lig_build::pos_t Kt      = K + lig_build::pos_t(25, 0);
   lig_build::pos_t Lt      = L + lig_build::pos_t(25, 0);
   lig_build::pos_t Mt      = M + lig_build::pos_t(25, 0);
   
   std::string stroke_colour = "#008000";
   stroke_colour = "#0000cc"; // blue
   gboolean end_arrow = 1;
   GooCanvasItem *item_1 =
      goo_canvas_polyline_new_line(key_group,
				   A.x, A.y,
				   B.x, B.y,
				   "line-width", 2.5,
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   "end_arrow",   end_arrow,
				   NULL);
   GooCanvasItem *text_1 = goo_canvas_text_new(key_group,
					       "Main-chain acceptor",
					       AB_txt.x, AB_txt.y,
					       -1,
					       GTK_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);

   GooCanvasItem *item_2 =
      goo_canvas_polyline_new_line(key_group,
				   Ad.x, Ad.y,
				   Bd.x, Bd.y,
				   "line-width", 2.5,
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   "start_arrow",   end_arrow,
				   NULL);
   
   GooCanvasItem *text_2 = goo_canvas_text_new(key_group,
					       "Main-chain donor",
					       ABd_txt.x, ABd_txt.y,
					       -1,
					       GTK_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);
   
   stroke_colour = "#008000"; // green
   GooCanvasItem *item_3 =
      goo_canvas_polyline_new_line(key_group,
				   C.x, C.y,
				   D.x, D.y,
				   "line-width", 2.5,
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   "end_arrow",   end_arrow,
				   NULL);

   GooCanvasItem *text_3 = goo_canvas_text_new(key_group,
					       "Side-chain acceptor",
					       CD_txt.x, CD_txt.y,
					       -1,
					       GTK_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);
   
   GooCanvasItem *item_4 =
      goo_canvas_polyline_new_line(key_group,
				   Cd.x, Cd.y,
				   Dd.x, Dd.y,
				   "line-width", 2.5,
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   "start_arrow",   end_arrow,
				   NULL);

   GooCanvasItem *text_4 = goo_canvas_text_new(key_group,
					       "Side-chain donor",
					       CDd_txt.x, CDd_txt.y,
					       -1,
					       GTK_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);

   std::string lime = "#888820";
   GooCanvasItem *item_5 =
      goo_canvas_polyline_new_line(key_group,
				   WatA.x, WatA.y,
				   WatB.x, WatB.y,
				   "line-width", 2.5,
				   "line-dash", dash,
				   "stroke-color", lime.c_str(),
				   NULL);

   GooCanvasItem *text_5 = goo_canvas_text_new(key_group,
					       "H-bond to Water",
					       Wat_txt.x, Wat_txt.y,
					       -1,
					       GTK_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);
   GooCanvasItem *item_6 =
      goo_canvas_polyline_new_line(key_group,
				   MetalA.x, MetalA.y,
				   MetalB.x, MetalB.y,
				   "line-width", 2.5,
				   "line-dash", dash,
				   "stroke-color",  "#990099",
				   NULL);

   GooCanvasItem *text_6 = goo_canvas_text_new(key_group,
					       "Metal Bond",
					       Metal_txt.x, Metal_txt.y,
					       -1,
					       GTK_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);

   
   // circles
   double key_residue_radius = 13;
   double line_width = 1.0;
   std::string grease = "#ccffbb";
   std::string purple = "#eeccee";
   std::string red    = "#cc0000";
   std::string blue   = "#0000cc";
   std::string metalic_grey = "#d9d9d9";

   GooCanvasItem *circle_grease = 
      goo_canvas_ellipse_new(key_group,
			     E.x, E.y, key_residue_radius, key_residue_radius,
			     "line_width", line_width,
			     "fill-color", grease.c_str(),
			     NULL);
   GooCanvasItem *circle_polar = 
      goo_canvas_ellipse_new(key_group,
			     F.x, F.y, key_residue_radius, key_residue_radius,
			     "line_width", line_width,
			     "fill-color", purple.c_str(),
			     NULL);
   GooCanvasItem *circle_acidic = 
      goo_canvas_ellipse_new(key_group,
			     G.x, G.y, key_residue_radius, key_residue_radius,
			     "stroke-color", red.c_str(),
			     "fill-color", purple.c_str(),
			     "line_width", 2.0,
			     NULL);
   GooCanvasItem *circle_basic = 
      goo_canvas_ellipse_new(key_group,
			     H.x, H.y, key_residue_radius, key_residue_radius,
			     "fill-color", purple.c_str(),
			     "stroke-color", blue.c_str(),
			     "line_width", 2.0,
			     NULL);
   GooCanvasItem *circle_water = 
      goo_canvas_ellipse_new(key_group,
			     I.x, I.y, key_residue_radius, key_residue_radius,
			     "line_width", line_width,
			     "fill-color", "white",
			     NULL);
   GooCanvasItem *circle_metal = 
      goo_canvas_ellipse_new(key_group,
			     J.x, J.y, key_residue_radius, key_residue_radius,
			     "line_width", line_width,
			     "fill-color", metalic_grey.c_str(),
			     NULL);

   GooCanvasItem *text_grease = goo_canvas_text_new(key_group,
						    "Grease",
						    Et.x, Et.y,
						    -1,
						    GTK_ANCHOR_WEST,
						    "font", "Sans 10",
						    "fill_color", "#111111",
						    NULL);
   GooCanvasItem *text_polar = goo_canvas_text_new(key_group,
						   "Polar",
						   Ft.x, Ft.y,
						   -1,
						   GTK_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);
   GooCanvasItem *text_acidic = goo_canvas_text_new(key_group,
						   "Acidic",
						   Gt.x, Gt.y,
						   -1,
						   GTK_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);
   GooCanvasItem *text_basic = goo_canvas_text_new(key_group,
						   "Basic",
						   Ht.x, Ht.y,
						   -1,
						   GTK_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);
   GooCanvasItem *text_water = goo_canvas_text_new(key_group,
						   "Water",
						   It.x, It.y,
						   -1,
						   GTK_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);
   
   GooCanvasItem *text_metal = goo_canvas_text_new(key_group,
						   "Metal",
						   Jt.x, Jt.y,
						   -1,
						   GTK_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);

   // Solvent accessibility of atoms
   bool save_state = draw_flev_annotations_flag;
   draw_flev_annotations_flag = true;
   lig_build::pos_t K_dx = K + lig_build::pos_t(-4, 0); // K + delta x to make objects line up...
   draw_solvent_accessibility_of_atom(K_dx, 0.15, key_group);
   draw_flev_annotations_flag = save_state;
   GooCanvasItem *text_solvent_atom = goo_canvas_text_new(key_group,
							  "Solvent exposure",
							  Kt.x, Kt.y,
							  -1,
							  GTK_ANCHOR_WEST,
							  "font", "Sans 10",
							  "fill_color", "#111111",
							  NULL);
   
   // solvent exposue ligand difference on environment/receptor residue
   GooCanvasItem *circle_residue_exposure_1 = 
      goo_canvas_ellipse_new(key_group,
			     L.x-5.0, L.y-5.0, key_residue_radius * 1.3 , key_residue_radius * 1.3,
			     "line_width", 0.0,
			     "fill-color", "#bbbbff",
			     NULL);
   GooCanvasItem *circle_residue_exposure_2 = 
      goo_canvas_ellipse_new(key_group,
			     L.x, L.y, key_residue_radius * 0.9, key_residue_radius * 0.9,
			     "line_width", line_width,
			     "fill-color", "white",
			     NULL);
   GooCanvasItem *text_solvent_residue = goo_canvas_text_new(key_group,
							     "Residue protection",
							     Lt.x, Lt.y,
							     -1,
							     GTK_ANCHOR_WEST,
							     "font", "Sans 10",
							     "fill_color", "#111111",
							     NULL);

   std::vector<std::vector<lig_build::pos_t> > contour_lines;
   std::vector<lig_build::pos_t> contour_set;

   for (unsigned int i=0; i<11; i++) {
      double delta_x = 15 * sin(i * 2 * M_PI * 0.1);
      double delta_y = 12 * cos(i * 2 * M_PI * 0.1);
      lig_build::pos_t p = M + lig_build::pos_t(delta_x-3, delta_y);
      contour_set.push_back(p);
   }
   contour_lines.push_back(contour_set);
   lig_build::pos_t dummy_pos;
   ligand_grid grid(dummy_pos, dummy_pos);
   grid.plot_contour_lines(contour_lines, key_group);

   GooCanvasItem *text_substitution_contour = goo_canvas_text_new(key_group,
								  "Substitution Contour",
								  Mt.x, Mt.y,
								  -1,
								  GTK_ANCHOR_WEST,
								  "font", "Sans 10",
								  "fill_color", "#111111",
								  NULL);
}


void
lbg_info_t::hide_key() {

   if (key_group) { 
      gint n_children = goo_canvas_item_get_n_children (key_group);
      for (int i=0; i<n_children; i++)
	 goo_canvas_item_remove_child(key_group, 0);
      key_group = NULL;
   }
}

bool
lbg_info_t::annotate(const std::vector<std::pair<coot::atom_spec_t, float> > &s_a_v,
		     const std::vector<coot::fle_residues_helper_t> &centres,
		     const std::vector<int> &additional_representation_handles_in,
		     const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
		     const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
		     const coot::flev_attached_hydrogens_t &ah,
		     const coot::pi_stacking_container_t &pi_stack_info,
		     const coot::dictionary_residue_restraints_t &restraints) {


   if (0) {
      for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	 std::cout << "  =============== lbg::annotate() bond to ligand " << ib << " "
		   << bonds_to_ligand[ib].ligand_atom_spec.atom_name
		   << " by residue " << bonds_to_ligand[ib].res_spec.chain << " "
		   << bonds_to_ligand[ib].res_spec.resno << " type: "
		   << bonds_to_ligand[ib].bond_type
		   << std::endl;
      }
   }

   if (0) { 
      std::cout << "--------------------------------------------------------------" << std::endl;
      std::cout << "======== lbg_info_t::annotate() here are bash distances for atoms:" << std::endl;
      std::map<std::string, std::vector<coot::bash_distance_t> >::const_iterator it;
      for (it=ah.atom_bashes.begin(); it!=ah.atom_bashes.end(); it++) {
	 std::cout << it->first << " " << it->second.size() << " bashes " << std::endl;
	 for (unsigned int i=0; i<it->second.size(); i++) { 
	    std::cout << "   " << it->second[i] << std::endl;
	 }
      }
      std::cout << "--------------------------------------------------------------" << std::endl;
   }

   bool r = false;
   std::vector<solvent_accessible_atom_t> solvent_accessible_atoms =
      convert(s_a_v, ah);
   widgeted_molecule_t new_mol = mol;
   new_mol.map_solvent_accessibilities_to_atoms(solvent_accessible_atoms);

   // fill class data item std::vector<residue_circle_t> residue_circles;
   // 
   residue_circles.clear();
   for (unsigned int i=0; i<centres.size(); i++) {

      if (0)
	 std::cout << "debugg:: in lbg_info_t::annotate() handling circle " << i << " of "
		   << centres.size() << std::endl;
      
      std::string label = centres[i].spec.chain;
      label += coot::util::int_to_string(centres[i].spec.resno);
      label += centres[i].spec.insertion_code;
      
      residue_circle_t circle(centres[i].centre.x(),
			      centres[i].centre.y(),
			      centres[i].centre.z(),
			      centres[i].spec,
			      centres[i].residue_name,
			      label);
      
      clipper::Coord_orth cp(centres[i].centre.x(),
			     centres[i].centre.y(),
			     centres[i].centre.z());
      lig_build::pos_t pos = mol.input_coords_to_canvas_coords(cp);
      circle.set_canvas_pos(pos);
      if (centres[i].residue_name == "HOH") {
	 for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	    if (bonds_to_ligand[ib].res_spec == centres[i].spec) {
	       circle.set_water_dist_to_protein(bonds_to_ligand[ib].water_protein_length);
	    }
	 }
      }

      // now add any H-bonds to the ligand:
      // 
      for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	 if (bonds_to_ligand[ib].res_spec == centres[i].spec) {
	    double bond_l = bonds_to_ligand[ib].bond_length;
	    std::string ligand_atom_name = bonds_to_ligand[ib].ligand_atom_spec.atom_name;
	    lbg_info_t::bond_to_ligand_t btl(ligand_atom_name, bond_l);
	    btl.bond_type = bonds_to_ligand[ib].bond_type;
	    circle.add_bond_to_ligand(btl);
	 }
      }

      // Now add any pi-stacking interactions to/from the ligand:
      //
      for (unsigned int istack=0; istack<pi_stack_info.stackings.size(); istack++) {
	 coot::residue_spec_t spec(pi_stack_info.stackings[istack].res);
	 if (spec == centres[i].spec) {
	    if (pi_stack_info.stackings[istack].type == coot::pi_stacking_instance_t::PI_PI_STACKING) {
	       std::vector<std::string> lra = pi_stack_info.stackings[istack].ligand_ring_atom_names;
	       circle.set_stacking("pi-pi", lra, "");
	    }
	    if (pi_stack_info.stackings[istack].type == coot::pi_stacking_instance_t::PI_CATION_STACKING) {
	       std::vector<std::string> lra = pi_stack_info.stackings[istack].ligand_ring_atom_names;
	       circle.set_stacking("pi-cation", lra, "");
	    }
	    // cation in ligand, PHE (say)
	    if (pi_stack_info.stackings[istack].type == coot::pi_stacking_instance_t::CATION_PI_STACKING) {
	       std::vector<std::string> lra_null;
	       circle.set_stacking("cation-pi", lra_null, pi_stack_info.stackings[istack].ligand_cationic_atom_name);
	    }
	 }
      }

      // solvent exposure difference annotation.
      //
      // std::vector<coot::solvent_exposure_difference_helper_t> &sed,
      //
      for (unsigned int ised=0; ised<sed.size(); ised++) { 
	 if (sed[ised].res_spec == centres[i].spec) {
	    circle.set_solvent_exposure_diff(sed[ised].exposure_fraction_holo,
					     sed[ised].exposure_fraction_apo);
	 }
      }
      residue_circles.push_back(circle);
   }

   // additional_representation_handles transfer
   additional_representation_handles = additional_representation_handles_in;


   // ligand ring centres (some of which may be aromatic and are in pi_stack_info too).
   //
   // However, for the substitution contour and the initial layout
   // ring avoidance, we blob in a circle of radius 1/(2*sin(180/n_ring_atoms)) bond lengths.

   // sets the object variable
   ring_atoms_list = restraints.get_ligand_ring_list();
   
// just checking that it was passed correctly - it is, for the
// moment, but do check 1x8b ligand.
//    
//    for (unsigned int i=0; i<ring_atoms_list.size(); i++) {
//       std::cout << "ring list " << i << "   ";
//       for (unsigned int j=0; j<ring_atoms_list[i].size(); j++) { 
// 	 std::cout << ring_atoms_list[i][j] << "  ";
//       }
//       std::cout << std::endl;
//    }

   render_from_molecule(new_mol);
   refine_residue_circle_positions();

   // has the current solution problems due to residues too close to the ligand?
   std::pair<bool, std::vector<int> > problem_status = solution_has_problems_p();
   std::cout << "::::::::: problem status: " << problem_status.first << std::endl;
   for (unsigned int ip=0; ip<problem_status.second.size(); ip++) {
      std::cout << ":::::::::::: "
		<< residue_circles[problem_status.second[ip]].residue_label << " "
		<< residue_circles[problem_status.second[ip]].residue_type << " "
		<< residue_circles[problem_status.second[ip]].pos
		<< std::endl;
   }
   
   if (problem_status.second.size()) {
      // fiddle with residue_circles and reoptimise.
      std::vector<int> primary_indices = get_primary_indices();
      reposition_problematics_and_reoptimise(problem_status.second, primary_indices);
   }

   recentre_considering_residue_centres();
   return r;
}


std::vector<solvent_accessible_atom_t>
lbg_info_t::convert(const std::vector<std::pair<coot::atom_spec_t, float> > &s_a_v,
		    const coot::flev_attached_hydrogens_t &ah) const {

   std::vector<solvent_accessible_atom_t> r(s_a_v.size());

   for (unsigned int i=0; i<s_a_v.size(); i++) {
      std::string name = s_a_v[i].first.atom_name;
      r[i].atom_name = name;
      r[i].solvent_accessibility = s_a_v[i].second;

      std::map<std::string, std::vector<coot::bash_distance_t> >::const_iterator it =
	 ah.atom_bashes.find(name);
      if (it != ah.atom_bashes.end()) 
	 r[i].bash_distances = it->second;
   }

   return r;
} 


void
lbg_info_t::update_statusbar_smiles_string(const std::string &smiles_string) const {

   // "memory leak" here because we don't dispose of
   // (gtk_statusbar_pop()) the previous label.
   //
   std::string status_string;
   if (lbg_statusbar) { 
      status_string = " SMILES:  ";
      status_string += smiles_string;
      guint statusbar_context_id =
	 gtk_statusbar_get_context_id(GTK_STATUSBAR(lbg_statusbar), status_string.c_str());
      gtk_statusbar_push(GTK_STATUSBAR(lbg_statusbar),
			 statusbar_context_id,
			 status_string.c_str());
   }
}

void
lbg_info_t::update_statusbar_smiles_string() const {

   std::string status_string;
   try {
      std::string s = get_smiles_string_from_mol();
      update_statusbar_smiles_string(s);

   }
   catch (std::exception rte) {
      std::cout << rte.what() << std::endl;
   }
}

#ifdef MAKE_ENTERPRISE_TOOLS
void
lbg_info_t::update_statusbar_smiles_string(const RDKit::RWMol &mol) const {

   std::string s;
   if (mol.getNumAtoms() > 0) {
      // new
      RDKit::ROMol mol_copy(mol);
      s = RDKit::MolToSmiles(mol_copy);
   } 
   update_statusbar_smiles_string(s);
}
#endif

void
lbg_info_t::show_atom_x_dialog() {

   // show dialog here.
   if (lbg_atom_x_dialog && lbg_atom_x_entry) {
      gtk_entry_set_text(GTK_ENTRY(lbg_atom_x_entry), atom_X.c_str());
      gtk_widget_show(lbg_atom_x_dialog);
   } else {
      std::cout << "lbg_atom_x_dialog null " << std::endl;
   } 
} 

void
lbg_info_t::set_atom_x_string(const std::string &s) {
   atom_X = s;
}


// uses lbg_get_drug_entry
void
lbg_info_t::get_drug_using_entry_text() {

   const char *txt = gtk_entry_get_text(GTK_ENTRY(lbg_get_drug_entry));
   get_drug(txt);
} 

// get mol file and load it
void
lbg_info_t::get_drug(const std::string &drug_name) {

   if (get_drug_mdl_file_function_pointer) {
      try {
	 std::string file_name = get_drug_mdl_file_function_pointer(drug_name);
	 import_mol_from_file(file_name);
	 save_molecule();
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      } 
   } else {
      std::cout << "Null get-drug function" << std::endl;
   } 
}


#ifdef MAKE_ENTERPRISE_TOOLS   
#ifdef USE_PYTHON   
PyObject *
lbg_info_t::get_callable_python_func(const std::string &module_name,
				     const std::string &function_name) const {

   PyObject *extracted_func = NULL;

   //
   // Build the name object
   PyObject *pName = PyString_FromString(module_name.c_str());
   // Load the module object
   PyObject *pModule = PyImport_Import(pName);
   if (pModule == NULL) {
      // This happens, for example, when we can't pick up the rdkit installation properly
      // 
      std::cout << "Null pModule" << std::endl;
   } else { 
      // pDict is a borrowed reference 
      PyObject *pDict = PyModule_GetDict(pModule);
      if (! PyDict_Check(pDict)) {
	 std::cout << "pDict is not a dict"<< std::endl;
      } else {
	 // pFunc is also a borrowed reference 
	 PyObject *pFunc = PyDict_GetItemString(pDict, function_name.c_str());
	 if (! pFunc) {
	    std::cout << "pFunc is NULL" << std::endl;
	 } else { 
	    if (PyCallable_Check(pFunc)) {
	       extracted_func = pFunc;
	       if (1)
		  std::cout << "found " << module_name << " " << function_name << " function"
			    << std::endl;
	    }
	 }
      }
   }
   return extracted_func;
}
#endif
#endif


void
lbg_info_t::pe_test_function() {

#ifdef MAKE_ENTERPRISE_TOOLS   
   std::cout << "PE test function" << std::endl;

   // RECAP mol

   // first get the function:
   // 
   // PyObject *recap_func = get_callable_python_func("rdkit.Chem.Recap", "RecapDecompose");
   PyObject *recap_func = get_callable_python_func("recap_wrapper", "wrap_recap");

   if (PyCallable_Check(recap_func)) {
      PyObject *arg_list = PyTuple_New(1);
      RDKit::RWMol rdkm = rdkit_mol(mol);
      coot::rdkit_mol_sanitize(rdkm);
      
      RDKit::ROMol *mol_copy_p = new RDKit::ROMol(rdkm);
      boost::shared_ptr<RDKit::ROMol> xx(mol_copy_p);
      boost::python::object obj(xx);
      PyObject *rdkit_mol_py = obj.ptr();
      PyTuple_SetItem(arg_list, 0, rdkit_mol_py);
      PyObject *result = PyEval_CallObject(recap_func, arg_list);
      if (! result) {
	 std::cout << "Null result in pe test function" << std::endl;
      } else {
	 std::cout << "Result: " << result << std::endl;
      }
   } 
#endif 
} 



#endif // HAVE_GOOCANVAS
