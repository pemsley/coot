/* lbg/lbg.cc
 *
 * Copyright 2010, 2011, 2012 by The University of Oxford
 * Copyright 2012, 2013, 2014, 2015, 2016 by Medical Research Council
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

#include "compat/coot-sysdep.h"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include <RDGeneral/versions.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#ifdef USE_PYTHON
#include <boost/python.hpp>
#endif
#endif

#include <cairo.h>
#if CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#if CAIRO_HAS_PS_SURFACE
#include <cairo-ps.h>
#endif
#include <cairo-svg.h>
#include "lbg.hh"
#include "lbg-drag-and-drop.hh"
#include "qed-interface.hh" // interface to silicos-it biscu-it python function

#include <gdk/gdkkeysyms.h> // for keyboarding.

#if ( ( (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 11) ) || GTK_MAJOR_VERSION > 2)
//
lbg_info_t *
lbg(lig_build::molfile_molecule_t mm,
    std::pair<bool, coot::residue_spec_t> ligand_spec_pair,
    mmdb::Manager *mol,
    const std::string &view_name,
    const std::string &molecule_file_name,
    int imol,
    coot::protein_geometry *geom_p_in,
    bool use_graphics_interface_flag,
    bool stand_alone_flag_in,
    int (*get_url_func_pointer_in) (const char *s1, const char *s2),
    void (*prodrg_import_function_pointer) (std::string file_name, std::string comp_id),
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

	 GtkBuilder *builder = gtk_builder_new();
	 guint add_from_file_status =
	    gtk_builder_add_from_file(builder, glade_file_full.c_str(), NULL);

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

	    lbg = new lbg_info_t(imol, geom_p_in);
	    if (stand_alone_flag_in) {
	       lbg->set_stand_alone();
	    }
	    int init_status = lbg->init(builder);

	    // normal built-in and using-graphics path...

	    if (! init_status) {
	       std::cout << "ERROR:: lbg init failed." << std::endl;
	    } else {

	       // Happy Path

	       gtk_builder_connect_signals (builder, lbg->canvas);
	       g_object_unref (G_OBJECT (builder));
	       if (ligand_spec_pair.first)
		  lbg->set_ligand_spec(ligand_spec_pair.second);

	       mmdb::Residue *residue_p = coot::get_first_residue_helper_fn(mol);
	       if (residue_p) {
		  std::string res_name = residue_p->GetResName();
		  gtk_label_set_text(GTK_LABEL(lbg->lbg_toolbar_layout_info_label), res_name.c_str());
	       }

	       widgeted_molecule_t wmol = lbg->import_mol_file(mm, molecule_file_name, mol);

	       lbg->import_from_widgeted_molecule(wmol);
	       lbg->render();
	       lbg->update_descriptor_attributes();
               lbg->update_apply_button_sensitivity();
	       lbg->save_molecule();
	    }
	 }
      }

      if (! use_graphics_interface_flag) {
	 lbg = new lbg_info_t(imol, geom_p_in);
	 lbg->no_graphics_mode(); // sets flag for "don't display a window"
	 int init_status = lbg->init(NULL); // setup canvas, but not gui.
	 if (ligand_spec_pair.first)
	    lbg->set_ligand_spec(ligand_spec_pair.second);

	 widgeted_molecule_t wmol = lbg->import_mol_file(mm, molecule_file_name, mol);
	 lbg->import_from_widgeted_molecule(wmol);
	 lbg->render();
	 lbg->update_descriptor_attributes(); // why?
      }

      if (lbg) {
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
	 } else {
	    // if (lbg->lbg_get_drug_menuitem)
	    // gtk_widget_set_sensitive(GTK_WIDGET(lbg->lbg_get_drug_menuitem), FALSE);
	 }
      }
   }
   return lbg;
}
#endif // GTK_VERSION

void
lbg_info_t::new_lbg_window() {

   lig_build::molfile_molecule_t blank_mol;
   std::pair<bool, coot::residue_spec_t> ligand_spec_pair(0, 0);
   mmdb::Manager *mol = NULL;
   std::string view_name;
   std::string molecule_file_name;
   int imol = -1;
   bool use_graphics_interface_flag = true;

   lbg(blank_mol, ligand_spec_pair, mol, view_name, molecule_file_name,
       imol,
       geom_p,
       use_graphics_interface_flag,
       is_stand_alone(),
       get_url_func_ptr,
       prodrg_import_func_ptr,
       sbase_import_func_ptr,
       get_drug_mdl_file_function_pointer);
}


// when feeling righteous, remove this and put it into geometry and
// remove the one in coot-utils
//
mmdb::Residue *
coot::get_first_residue_helper_fn(mmdb::Manager *mol) {

   mmdb::Residue *res = NULL;
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
	 mmdb::Chain *chain_p;

	 int n_chains = model_p->GetNumberOfChains();
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    chain_p = model_p->GetChain(i_chain);
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::Residue *residue_p;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       if (residue_p) {
		  res = residue_p;
		  break;
	       }
	    }
	    if (res)
	       break;
	 }
      }
   }
   return res;
}



GtkWidget *get_canvas_from_scrolled_win(GtkWidget *canvas) {

   return canvas;
}

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
//
// widgeted_molecule_t mol needs to (store and?) transfer chiral info to the RDKit molecule.
// Sounds tricky.
//
RDKit::RWMol
lbg_info_t::rdkit_mol(const widgeted_molecule_t &mol) const {

   RDKit::RWMol m;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   std::pair<double, lig_build::pos_t> scale_and_correction = mol.current_scale_and_centre();

   if (0) {
      std::cout << "------ in rdkit_mol(widgeted_molecule_t) atoms.size() " << mol.atoms.size() << std::endl;
      std::cout << "------ in rdkit_mol(widgeted_molecule_t) bonds.size() " << mol.bonds.size() << std::endl;
      std::cout << ":::::::; in lbg_info_t::rdkit_mol() mol's scale correction is "
		<< mol.scale_correction.first << " " << mol.scale_correction.second
		<< std::endl;
      std::cout << ":::::::; in lbg_info_t::rdkit_mol() mol's centre correction is "
		<< mol.centre_correction << std::endl;
      std::cout << ":::::::; in lbg_info_t::rdkit_mol() mol's current scale and correction is "
		<< scale_and_correction.first << " " << scale_and_correction.second << std::endl;
   }

   double local_scale = 1.0;
   lig_build::pos_t local_offset(0,0);
   if (scale_and_correction.first > 0) {
      local_scale  = 1.0/scale_and_correction.first;
      local_offset = scale_and_correction.second;
   }

   // we need to account for atoms that are closed.  So we make
   // atom_index[idx] to return the index of the atom in the rdkit
   // molecule.  i.e. converting from (possibly closed) atoms in a
   // lig_build::molecule
   //
   std::vector<int> atom_index(mol.atoms.size());
   int at_no = 0;

   RDKit::Conformer *conf = new RDKit::Conformer(mol.num_atoms());

   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      if (! mol.atoms[iat].is_closed()) {
	 RDKit::Atom *at = new RDKit::Atom;
	 at->setAtomicNum(tbl->getAtomicNumber(mol.atoms[iat].element));
	 at->setFormalCharge(mol.atoms[iat].charge);

	 // add the name to "at" too here (if you can).
	 //
	 try {
	    at->setProp("name", std::string(mol.atoms[iat].get_atom_name()));
	    std::string ele_capped =
	       coot::util::capitalise(coot::util::remove_leading_spaces(mol.atoms[iat].element));
	    int atomic_number = tbl->getAtomicNumber(ele_capped);
	    at->setAtomicNum(atomic_number);
	    // at->setMass(tbl->getAtomicWeight(atomic_number));
	    at->setIsotope(0);

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

	    // canvas positions are upsidedown in y
	    RDGeom::Point3D pos(+local_scale * (mol.atoms[iat].atom_position.x - local_offset.x),
				-local_scale * (mol.atoms[iat].atom_position.y - local_offset.y), 0);
	    conf->setAtomPos(at_no, pos);
	    atom_index[iat] = at_no++;
	 }
	 catch (const boost::bad_lexical_cast &blc) {
	    std::cout << "on making atom bad_lexical_cast " << iat << " "
		      << blc.what() << std::endl;
	 }
	 catch (const std::exception &rte) {
	    std::cout << "on making atom " << iat << " " << rte.what() << std::endl;
	 }
      } else {
	 // std::cout << "debug:: atom " << iat << " is closed " << std::endl;
      }
   }

   for (unsigned int ib=0; ib<mol.bonds.size(); ib++) {
      if (! mol.bonds[ib].is_closed()) {
	 RDKit::Bond::BondType type = convert_bond_type(mol.bonds[ib].get_bond_type());
	 RDKit::Bond::BondDir   dir = convert_bond_dir (mol.bonds[ib].get_bond_type());
	 RDKit::Bond *bond = new RDKit::Bond(type);
	 if (dir != RDKit::Bond::NONE)
	    bond->setBondDir(dir);
	 int idx_1 = mol.bonds[ib].get_atom_1_index();
	 int idx_2 = mol.bonds[ib].get_atom_2_index();
	 if (0)
	    std::cout << "   rdkit_mol(wmol) convert bond for atoms "
		      << idx_1 << " " << idx_2
		      << " from type " << mol.bonds[ib].get_bond_type()
		      << " to type " << type << std::endl;
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
	 } else {
	    delete bond;
	 }
      } else {
	 // std::cout << "debug:: bond " << ib << " is closed " << std::endl;
      }
   }
   m.addConformer(conf, true);
   return m;
}
#endif

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
RDKit::Bond::BondType
lbg_info_t::convert_bond_type(const lig_build::bond_t::bond_type_t &t) const {

   // There are lots more in RDKit::Bond::BondType!
   //
   RDKit::Bond::BondType bt = RDKit::Bond::UNSPECIFIED; // lig_build::bond_t::BOND_UNDEFINED
   if (t == lig_build::bond_t::SINGLE_BOND)
      bt = RDKit::Bond::SINGLE;
   if (t == lig_build::bond_t::IN_BOND)
      bt = RDKit::Bond::SINGLE;
   if (t == lig_build::bond_t::OUT_BOND)
      bt = RDKit::Bond::SINGLE;
   if (t == lig_build::bond_t::DOUBLE_BOND)
      bt = RDKit::Bond::DOUBLE;
   if (t == lig_build::bond_t::TRIPLE_BOND)
      bt = RDKit::Bond::TRIPLE;

   // PUTS_IT_BACK?
   if (t == lig_build::bond_t::AROMATIC_BOND)
       bt = RDKit::Bond::AROMATIC;

   if (t == lig_build::bond_t::DELOC_ONE_AND_A_HALF)
      bt = RDKit::Bond::ONEANDAHALF;

   return bt;
}
#endif

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
RDKit::Bond::BondDir
lbg_info_t::convert_bond_dir(const lig_build::bond_t::bond_type_t &t) const {

   RDKit::Bond::BondDir bd = RDKit::Bond::NONE;
   if (t == lig_build::bond_t::IN_BOND)
      bd = RDKit::Bond::BEGINDASH;
   if (t == lig_build::bond_t::OUT_BOND)
      bd = RDKit::Bond::BEGINWEDGE;

   return bd;
}

#endif


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
// this can throw a std::exception
//
std::string
lbg_info_t::get_smiles_string_from_mol_rdkit() const {

   // here mol is a widgeted_molecule_t (no idea about aromaticity)
   //
   // What about chiral atoms, eh?
   //
   RDKit::RWMol rdkm = rdkit_mol(mol);
   coot::rdkit_mol_sanitize(rdkm);
   return get_smiles_string(rdkm);

//    RDKit::ROMol *rdk_mol_with_no_Hs = RDKit::MolOps::removeHs(rdkm);
//    bool doIsomericSmiles = true;
//    std::string s = RDKit::MolToSmiles(*rdk_mol_with_no_Hs, doIsomericSmiles);
//    delete rdk_mol_with_no_Hs;
//    return s;

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


static bool
on_canvas_button_press(GooCanvasItem  *item,
		       GooCanvasItem  *target_item,
		       GdkEventButton *event,
		       gpointer        user_data) {

   // currently the scroll wheel events are consumed before we get here - by the scrolled window perhaps

   int x_as_int = -1, y_as_int = -1;
   GdkModifierType state;
   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   bool button_1_is_pressed = false;
   if (state & GDK_BUTTON1_MASK)
      button_1_is_pressed = true;

   if (! event) {
      std::cout << "on_canvas_button_press() error NULL event!" << std::endl;
   } else {
      x_as_int = int(event->x);
      y_as_int = int(event->y);
   }

   if (! target_item) {
      lbg_info_t *l =
	 static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
      if (!l) {
	 // std::cout << "on button-click: NULL lbg-info from item" << std::endl;
      } else {
	 if (l->in_delete_mode_p()) {
	    l->handle_item_delete(event);
	 } else {
	    l->handle_item_add(event);
	 }

	 if (event) {
	    // std::cout << "here with event button: " << event->button << std::endl;

	    if (event->state & GDK_CONTROL_MASK) {
	       if (event->button == 4) { // scroll up
		  gdouble scale =  goo_canvas_get_scale(GOO_CANVAS(l->canvas));
		  scale *= 1.2;
		  goo_canvas_set_scale(GOO_CANVAS(l->canvas), 1.2);
	       }
	       if (event->button == 5) { // scroll down
		  gdouble scale =  goo_canvas_get_scale(GOO_CANVAS(l->canvas));
		  scale /= 1.2;
		  goo_canvas_set_scale(GOO_CANVAS(l->canvas), 1.2);
	       }
	    }
	 }

      }

   } else {

      coot::residue_spec_t *spec_p =
	 static_cast<coot::residue_spec_t *>(g_object_get_data (G_OBJECT (target_item), "spec"));
      clipper::Coord_orth *pos_p =
	 static_cast<clipper::Coord_orth *> (g_object_get_data (G_OBJECT (target_item), "position"));

      lbg_info_t *l = NULL;
      if (item)
	 l = static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));


//#ifdef MAKE_ENHANCED_LIGAND_TOOLS
      if (spec_p) {
	 // std::cout << "on_canvas_button_press(): clicked on " << *spec_p << std::endl;
	 if (event->type==GDK_2BUTTON_PRESS) { // double click
	    int imol = spec_p->int_user_data;
	    if (! l) {
	       // std::cout << "oops failed to get lbg_info_t pointer" << std::endl;
	    } else {
	       std::pair<bool, coot::residue_spec_t> ligand_spec_pair = l->get_ligand_spec();
	       // std::cout << "ligand pair first: " << ligand_spec_pair.first << std::endl;
	       if (ligand_spec_pair.first) {
		  // std::cout << " ligand_spec: " << ligand_spec_pair.second << std::endl;
		  // std::cout << "residue_spec: " << *spec_p << std::endl;
		  l->orient_view(imol,  ligand_spec_pair.second, *spec_p);
	       }
	    }
	 }
      }

      if (l) {
	 if (pos_p) {
	    l->set_rotation_centre(*pos_p);  // see graphics-c-interface-functions.hh
	 }
      }

//#endif

      if (!l) {
	 std::cout << "null lbg-info from target_item" << std::endl;
      } else {

	 bool handled = false;
#ifdef MAKE_ENHANCED_LIGAND_TOOLS

	 // did we pick on a molecule bond (does not consider annotation bonds).
	 handled = l->handle_bond_picking_maybe();
#endif

	 if (! handled) {

	    // note: it is not an error if l is null, that can happen
	    // when clicking on the atom names for example - they don't
	    // have a lbg-info attached (at the moment).

	    l->set_mouse_pos_at_click(x_as_int, y_as_int); // save for dragging

	    if (0)
	       std::cout << "   on click: scale_correction  "
			 << l->mol.scale_correction.first  << " "
			 << l->mol.scale_correction.second << " centre_correction "
			 << l->mol.centre_correction << std::endl;

	    if (l->in_delete_mode_p()) {
	       l->handle_item_delete(event);
	    } else {
	       l->handle_item_add(event);
	    }
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


	 // missing button up from a recently deleted canvas item?
	 if (0)
	    std::cout << "      here in on_canvas_button_release() clear button down "
		      << std::endl;

	 l->clear_button_down_bond_addition();
      } else {
	 // this should not happen
	 std::cout << "======= null lbg_info_t pointer " << std::endl;
      }
   } else {
      // std::cout << "======= non null target_item " << std::endl;
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

   most_recent_bond_made_new_atom_flag = false;
   atom_index_of_atom_to_which_the_latest_bond_was_added = -1; // no more rotations are allowed
}


static bool
on_canvas_motion(GooCanvasItem  *item,
		 GooCanvasItem  *target_item,
		 GdkEventMotion *event,
		 gpointer        user_data) {

   int x_as_int = -1, y_as_int = -1;
   GdkModifierType state;
   if (! event) {
      // std::cout << "on_canvas_button_press_new() error NULL event!" << std::endl;
   } else {
      x_as_int = int(event->x);
      y_as_int = int(event->y);
   }

   if (! target_item) {
      // std::cout << "on_canvas_motion(): NULL target item" << std::endl;
      lbg_info_t *l =
	 static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
      if (!l) {
	 // std::cout << "on_canvas_motion(): null lbg-info from target_item" << std::endl;
      }
   } else {

      lbg_info_t *l =
	 static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
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
	    l->set_show_additional_representation(imol, add_rep_handle, 1);
	    l->all_additional_representations_off_except(imol, add_rep_handle, 0);
	 }
      }
   }

   if (! item) {
      // std::cout << "on_canvas_motion_new(): NULL item!" << std::endl;
   } else {
      // std::cout << "on_canvas_motion_new():  item:" << item << std::endl;
      lbg_info_t *l =
	 static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (item), "lbg-info"));
      if (!l) {
	 // std::cout << "canvas_motion: null lbg-info from item!" << std::endl;
      } else {
	 if (event) {
	    guint state = event->state;
	    GdkModifierType g_state = static_cast<GdkModifierType> (state);
	    int x_as_int = int(event->x);
	    int y_as_int = int(event->y);

	    l->handle_drag(g_state, x_as_int, y_as_int);
	 }
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
	    // std::cout << "handle_drag(): calling rotate_latest_bond()" << std::endl;
	    rotate_latest_bond(x_as_int, y_as_int);
	 } else {
	    drag_canvas(x_as_int, y_as_int);
	 }
      }
   } else {
      if (button_down_bond_addition) {
	 if (highlight_data.single_atom()) {

	    most_recent_drag_atom = highlight_data;
	    // std::cout << "handle_drag(): extend_latest_bond_maybe()" << std::endl;
	    bool added_status = extend_latest_bond_maybe(); // checks for sensible atom to snap to.

	    if (!added_status) {
	       rotate_latest_bond(x_as_int, y_as_int);
	    }
	 }
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
      import_from_widgeted_molecule(new_mol);
      render();
      // std::cout << "calling update_descriptor_attributes() from clear_and_redraw() " << std::endl;
      // update_descriptor_attributes();
   } else {
      // this path gets called when the "Env Residues" button is pressed.
      // std::cout << "==== delta is zero path ==== " << std::endl;
      widgeted_molecule_t saved_mol = mol;
      // render_from_molecule(saved_mol);

      import_from_widgeted_molecule(saved_mol);
      render();
   }
   draw_all_flev_annotations();
}

void
lbg_info_t::clear_and_redraw() {

   clear_canvas();
   widgeted_molecule_t saved_mol = mol;
   import_from_widgeted_molecule(saved_mol);
   render();
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
   // (but should result in the same position).
   //

   // We need to know if the last atom was newly created as a result
   // of adding the bond we are now rotating. If it was, then we rotate
   // this new atom as we rotate the bond. If it was not, then we need to
   // create an atom (at the end of this bond) as a result of rotating this
   // bond.
   //

   bool debug = 0;

   // std::cout << "======================== start rotate_latest_bond() ============"
   // << std::endl;

   if (mol.bonds.size() > 0) {
      if (mol.atoms.size() > 0) {
	 GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

	 // most_recent_drag_atom is the rotating atom (the atom to which
	 // the newly created bond was attached).
	 // atom_index_of_atom_to_which_the_latest_bond_was_added is the
	 // other (rotating-bond-centre) atom.
	 //
	 if (false)
	    std::cout << "rotate_latest_bond(): with most_recent_drag_atom_index() "
		      << most_recent_drag_atom.get_atom_index() << " made-new-flag: "
		      << most_recent_bond_made_new_atom_flag << " "
		      << "atom_index_of_atom_to_which_the_latest_bond_was_added "
		      << atom_index_of_atom_to_which_the_latest_bond_was_added
		      << std::endl;

	 // widgeted_bond_t &bond = mol.bonds.back();

	 //
	 // widgeted_atom_t &atom = mol.atoms.back();

	 // idx is the index of the atom that is the other index of
	 // a newly created bond (often is a new (and last) atom
	 // but maybe not in "pathological" case)
	 //
	 unsigned int idx = ultimate_atom_index;

	 if (false)
	    std::cout << "rotate_latest_bond(): added-to-atom-idx: "
		      << atom_index_of_atom_to_which_the_latest_bond_was_added
		      << " and idx (\"to-atom\" index) " << idx << std::endl;

	 if (atom_index_of_atom_to_which_the_latest_bond_was_added != UNASSIGNED_INDEX) {

	    if (most_recent_bond_made_new_atom_flag) {

	       if (true) { // deubugging

		  if (false)
		     std::cout << "--------------------- rotate_bond block ------------"
			       << std::endl;

		  // we want to rotate atom idx

		  widgeted_atom_t &atom = mol.atoms[idx];

		  double x_cen = mol.atoms[atom_index_of_atom_to_which_the_latest_bond_was_added].atom_position.x;
		  double y_cen = mol.atoms[atom_index_of_atom_to_which_the_latest_bond_was_added].atom_position.y;

		  double x_at = atom.atom_position.x;
		  double y_at = atom.atom_position.y;
		  double theta_rad = atan2(-(y_mouse-y_cen), (x_mouse-x_cen));
		  double theta_target = theta_rad/DEG_TO_RAD;
		  double theta_current_rad = atan2(-(y_at-y_cen), (x_at-x_cen));
		  double theta_current = theta_current_rad/DEG_TO_RAD;
		  double degrees = theta_target - theta_current;

		  lig_build::pos_t new_atom_pos =
		     atom.atom_position.rotate_about(x_cen, y_cen, -degrees);

		  if (true) {

		     // move ultimate_atom_index atom to new_atom_pos
		     atom.atom_position = new_atom_pos;

		     if (penultimate_atom_index == UNASSIGNED_INDEX) {

			std::cout << "ERROR:: trapped UNASSIGNED_INDEX for penultimate_atom_index in "
				  << "rotate_latestbond()" << std::endl;

		     } else {

			unsigned int idx_2 = atom_index_of_atom_to_which_the_latest_bond_was_added;
			unsigned int bond_index = mol.get_bond_index(idx, idx_2);

			if (bond_index == lig_build::UNASSIGNED_BOND_INDEX) {
			   std::cout << "failed to find bond between " << idx << " " << idx_2 << std::endl;
			} else {
			   if (false) {
			      std::cout << "rotate_latest_bond(): rotate bond idx " << bond_index << std::endl;
			      std::cout << "rotate_latest_bond(): rotate_canvas_item for bond " << bond_index
					<< " atom indices "
					<< mol.bonds[bond_index].get_atom_1_index() << " "
					<< mol.bonds[bond_index].get_atom_2_index() << " "
					<< " about " << x_cen << " " << y_cen
					<< " -degrees: " << -degrees << std::endl;
			   }
			   mol.bonds[bond_index].rotate_canvas_item(x_cen, y_cen, -degrees);
			}
		     }
		  }
	       } // debug

	    } else { // of most_recent_bond_made_new_atom_flag test

	       if (false)
		  std::cout << "--------------------- cut the bond and rotate block ------------"
			    << std::endl;

	       // we don't want to rotate atom idx, we want to
	       // delete the latest bond (because it connect the "wrong" atom)
	       // create a new atom,
	       // rotate the bond to the new atom.

	       // atom idx is the moving atom

	       if (false)
		  std::cout << "rotate_latest_bond(): close bond between atom index "
			    << atom_index_of_atom_to_which_the_latest_bond_was_added
			    << " and " << idx << std::endl;

	       // create a new position, a new atom, add the atom, make a bond, add the bond.
	       //
	       widgeted_atom_t &atom = mol.atoms[idx];

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

		  // std::cout << "rotate_latest_bond(): rotate prevented -- too close to atom "
		  // << std::endl;

	       }  else {

		  unsigned int bond_index =
		     mol.get_bond_index(idx, atom_index_of_atom_to_which_the_latest_bond_was_added);

		  bool close_status = mol.close_bond(bond_index, root, true);
		  if (false)
		     std::cout << "rotate_latest_bond(): close bond between atom index "
			       << atom_index_of_atom_to_which_the_latest_bond_was_added
			       << " and " << idx << " close-status: " << close_status
			       << std::endl;

		  // make a new atom and bond
		  widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
		  std::pair<bool, int> checked_add = mol.add_atom(new_atom);
		  if (checked_add.first) {

		     int new_atom_index = checked_add.second;
		     if (false)
			std::cout << "rotate_latest_bond(): added atom! with new_atom_index "
				  << new_atom_index << std::endl;
		     lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
		     widgeted_atom_t atom = mol.atoms[atom_index_of_atom_to_which_the_latest_bond_was_added];
		     // add a bond to the atom at the rotation centre
		     std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_first_atom = mol.make_other_connections_to_first_atom_info(bond_index);
		     std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_second_atom = mol.make_other_connections_to_second_atom_info(bond_index);
		     widgeted_bond_t b(atom_index_of_atom_to_which_the_latest_bond_was_added,
				       new_atom_index, atom, new_atom,
				       false, true, bt,
				       other_connections_to_first_atom,
				       other_connections_to_second_atom, root);
		     int new_bond_idx = mol.add_bond(b);
		     if (new_bond_idx != lig_build::UNASSIGNED_BOND_INDEX) {
			penultimate_atom_index = atom_index_of_atom_to_which_the_latest_bond_was_added;
			ultimate_atom_index = new_atom_index;
			most_recent_bond_made_new_atom_flag = true;
		     }
		  }
	       }
	    }
	 }
      }
   }
}

// Return "we-created-a-new-bond" status.
//
bool
lbg_info_t::extend_latest_bond_maybe() {

   // we need to check that the highlighted atom is not the atom that
   // we just added the bond to, because (1) that will always be the
   // case just after the button down and the mouse moves with the
   // button down and (2) we don't want to draw a bond back to the
   // atom we just added a bond to.

   bool status = false;

   if (mol.bonds.size() > 0) {
      if (mol.atoms.size() > 0) {

	 if (highlight_data.has_contents()) {
	    if (highlight_data.single_atom()) {

	       widgeted_bond_t &bond = mol.bonds.back();
	       widgeted_atom_t &atom = mol.atoms[atom_index_of_atom_to_which_the_latest_bond_was_added];

	       int atom_index_snap_to = highlight_data.get_atom_index();
	       if (atom_index_snap_to != atom_index_of_atom_to_which_the_latest_bond_was_added) {
		  if (atom_index_snap_to != ultimate_atom_index) {

		     GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

		     mol.close_atom(ultimate_atom_index, root);

		     // and add a new bond.
		     //
		     lig_build::bond_t::bond_type_t bt =
			addition_mode_to_bond_type(canvas_addition_mode);

		     std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;

		     widgeted_bond_t b(atom_index_snap_to,
				       atom_index_of_atom_to_which_the_latest_bond_was_added,
				       mol.atoms[atom_index_snap_to],
				       mol.atoms[atom_index_of_atom_to_which_the_latest_bond_was_added],
				       false, false,
				       bt, empty, empty,
				       root);
		     int new_bond_idx = mol.add_bond(b); // can reject addition if bond
				                         // between the given atom
				                         // indices already exists.
		     std::cout << "extend_latest_bond_maybe() adding bond! " << b
			       << " new_bond_index: " << new_bond_idx << std::endl;

		     // std::cout << "Now we have " << mol.n_open_bonds() << " bonds" << std::endl;
		     status = true;
		  }
	       }
	    }
	 }
      }
   }
   return status;
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

gboolean
lbg_info_t::on_highlight_key_press_event (GooCanvasItem *item,
					  GooCanvasItem *target,
					  GdkEventKey *event,
					  gpointer data) {

//    std::cout << "on_highlight_key_press_event: item   " << item   << std::endl;
//    std::cout << "on_highlight_key_press_event: target " << target << std::endl;
//    std::cout << "on_highlight_key_press_event: event  " << event  << std::endl;
//    std::cout << "on_highlight_key_press_event: data   " << data   << std::endl;

   if (item) {
      GooCanvas *canvas_l = goo_canvas_item_get_canvas(item);
      GooCanvasItem *root = goo_canvas_get_root_item(canvas_l);
      lbg_info_t *l = static_cast<lbg_info_t *> (g_object_get_data (G_OBJECT (root), "lbg-info"));
      if (l) {
	 if (event) {
	    bool ctrl_is_pressed = event->state & GDK_CONTROL_MASK;
	    l->handle_key_press_button_toggle(event->keyval, ctrl_is_pressed);
	 }
      }
   }
   return FALSE;
}

void
lbg_info_t::handle_key_press_button_toggle(int keyval, bool ctrl_is_pressed) {

   switch (keyval) {

   case GDK_KEY_n:
      lbg_toggle_button_my_toggle(widget_names["nitrogen_toggle_toolbutton"]);
      break;

   case GDK_KEY_c:
      lbg_toggle_button_my_toggle(widget_names["carbon_toggle_toolbutton"]);
      break;

   case GDK_KEY_o:
      lbg_toggle_button_my_toggle(widget_names["oxygen_toggle_toolbutton"]);
      break;

   case GDK_KEY_p:
      lbg_toggle_button_my_toggle(widget_names["phos_toggle_toolbutton"]);
      break;

   case GDK_KEY_i:
      lbg_toggle_button_my_toggle(widget_names["iodine_toggle_toolbutton"]);
      break;

   case GDK_KEY_s:
      if  (ctrl_is_pressed)
	 lbg_toggle_button_my_toggle(widget_names["sulfur_toggle_toolbutton"]);
      else
	 lbg_toggle_button_my_toggle(widget_names["single_toggle_toolbutton"]);
      break;

   case GDK_KEY_d:
      lbg_toggle_button_my_toggle(widget_names["double_toggle_toolbutton"]);
      break;

   case GDK_KEY_1:
      lbg_toggle_button_my_toggle(widget_names["single_toggle_toolbutton"]);
      break;

   case GDK_KEY_2:
      lbg_toggle_button_my_toggle(widget_names["double_toggle_toolbutton"]);
      break;

   case GDK_KEY_3:
      lbg_toggle_button_my_toggle(widget_names["c3_toggle_toolbutton"]);
      break;

   case GDK_KEY_4:
      lbg_toggle_button_my_toggle(widget_names["c4_toggle_toolbutton"]);
      break;

   case GDK_KEY_5:
      lbg_toggle_button_my_toggle(widget_names["c5_toggle_toolbutton"]);
      break;

   case GDK_KEY_6:
      lbg_toggle_button_my_toggle(widget_names["c6_toggle_toolbutton"]);
      break;

   case GDK_KEY_7:
      lbg_toggle_button_my_toggle(widget_names["c7_toggle_toolbutton"]);
      break;

   case GDK_KEY_8:
      lbg_toggle_button_my_toggle(widget_names["c8_toggle_toolbutton"]);
      break;

   case GDK_KEY_b:
      // lbg_toggle_button_my_toggle(lbg_bromine_toggle_toolbutton);
      lbg_toggle_button_my_toggle(widget_names["c6_arom_toggle_toolbutton"]);
      break;

   default:
      std::cout << "lbg::handle_key_press_button_toggle():: nothing we know" << std::endl;
   }
}

void
lbg_info_t::lbg_toggle_button_my_toggle(GtkToggleToolButton *tb) {
   if (gtk_toggle_tool_button_get_active(tb))
      gtk_toggle_tool_button_set_active(tb, FALSE);
   else
      gtk_toggle_tool_button_set_active(tb, TRUE);
}


// set highlight_data
void
lbg_info_t::highlight_bond(const lig_build::bond_t &bond, bool delete_mode) {

   int i1 = bond.get_atom_1_index();
   int i2 = bond.get_atom_2_index();
   std::pair<int, int> bond_indices(i1, i2);
   lig_build::pos_t A = mol.atoms[i1].atom_position;
   lig_build::pos_t B = mol.atoms[i2].atom_position;

   // std::string col = "#20cc20";
   guint col = 0x20cc2080;
   if (delete_mode)
      col = 0xdd101080;

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
   GooCanvasItem *group = goo_canvas_group_new (root, NULL);
   GooCanvasItem *h_line =
      goo_canvas_polyline_new_line(group,
				   A.x, A.y,
				   B.x, B.y,
				   "line-width", 7.0, // in highlight_bond()
				   "stroke-color-rgba", col,
				   "can-focus", 1,
				   NULL);
   double rad = 8;
   GooCanvasItem *circle_1 = goo_canvas_ellipse_new(group, A.x, A.y, rad, rad,
						    "line_width", 0.0,
						    "fill-color-rgba", col,
						    NULL);
   GooCanvasItem *circle_2 = goo_canvas_ellipse_new(group, B.x, B.y, rad, rad,
						    "line_width", 0.0,
						    "fill-color-rgba", col, NULL);

   goo_canvas_grab_focus(GOO_CANVAS(canvas), h_line);

   g_signal_connect (h_line, "key_press_event",
		     G_CALLBACK (on_highlight_key_press_event), NULL);

   highlight_data = highlight_data_t(group, bond_indices, A, B);
   if (bond.have_centre_pos())
      highlight_data.set_ring_centre(bond.centre_pos());
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
   double x1 = A.x - width/2.0;
   double y1 = A.y - height/2.0;

   // std::cout << x1 << " " << x2 << " "<< width << " " << height << std::endl;

   GooCanvasItem *rect_item =
      goo_canvas_rect_new (root, x1, y1, width, height,
			   "line-width", 2.0, // in highlight_atom()
			   "stroke-color", col.c_str(),
			   "can-focus", 1,
			   NULL);

   goo_canvas_grab_focus(GOO_CANVAS(canvas), rect_item);

   g_signal_connect(rect_item, "key_press_event",
		     G_CALLBACK (on_highlight_key_press_event), NULL);

   highlight_data = highlight_data_t(rect_item, A, atom_index);
}


void
lbg_info_t::handle_item_add(GdkEventButton *event) {

   bool changed_status = 0;

   int x_as_int, y_as_int;
   GdkModifierType state;
   bool shift_is_pressed = false;
   if (event->state & GDK_SHIFT_MASK)
      shift_is_pressed = true;

   // set button_1_is_pressed:
   bool button_1_is_pressed = false;
   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   if (state & GDK_BUTTON1_MASK)
      button_1_is_pressed = true;

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
      changed_status = try_add_or_modify_bond(canvas_addition_mode, x_as_int, y_as_int,
					      button_1_is_pressed);
   }

   if (canvas_addition_mode == lbg_info_t::CHARGE) {
      changed_status = handle_charge_change();
   }

   if (changed_status) {
      save_molecule();
      update_descriptor_attributes();
      update_apply_button_sensitivity();
   }
}

void
lbg_info_t::update_apply_button_sensitivity() {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   if (use_graphics_interface_flag) {
      try {
	 RDKit::RWMol rdkm = rdkit_mol(mol);
	 coot::rdkit_mol_sanitize(rdkm);
         // restore sensitive (if needed)
         gtk_widget_set_sensitive(GTK_WIDGET(lbg_apply_button), TRUE);
      }
      catch (const std::exception &e) {
	 std::cout << "WARNING:: Oops non-sane molecule " << e.what() << std::endl;
         gtk_widget_set_sensitive(GTK_WIDGET(lbg_apply_button), FALSE);
      }
   }
#endif
}

void
lbg_info_t::update_descriptor_attributes() {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   if (use_graphics_interface_flag) {
      try {
	 RDKit::RWMol rdkm = rdkit_mol(mol);
	 coot::rdkit_mol_sanitize(rdkm);
	 update_statusbar_smiles_string(rdkm);
	 update_qed(rdkm);
	 update_alerts(rdkm);
      }
      catch (const std::exception &e) {
	 std::cout << "WARNING:: from update_descriptor_attributes() " << e.what() << std::endl;

	 if (lbg_statusbar) {
	    std::string status_string;
	    guint statusbar_context_id =
	       gtk_statusbar_get_context_id(GTK_STATUSBAR(lbg_statusbar), status_string.c_str());
	    gtk_statusbar_push(GTK_STATUSBAR(lbg_statusbar),
			       statusbar_context_id,
			       status_string.c_str());
	    // QED progress bar
	    gtk_label_set_text(GTK_LABEL(lbg_qed_text_label), "");
	    gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(lbg_qed_progressbar), 0);
	    reset_qed_properties_progress_bars();
	    // alerts
	    gtk_widget_hide(lbg_alert_hbox);
	    clear_canvas_alerts();
	 }
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
	 unsigned int ind_1 = highlight_data.get_bond_indices().first;
	 unsigned int ind_2 = highlight_data.get_bond_indices().second;
	 unsigned int bond_index = mol.get_bond_index(ind_1, ind_2);

	 mol.close_bond(bond_index, root, 1);
      }
      save_molecule();
      update_descriptor_attributes();
      update_apply_button_sensitivity();
   }
}

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
bool
lbg_info_t::handle_bond_picking_maybe() {

   bool handled = false;
   if (bond_pick_pending) {
      if (highlight_data.has_contents()) {
	 if (!highlight_data.single_atom()) {
	    int ind_1 = highlight_data.get_bond_indices().first;
	    int ind_2 = highlight_data.get_bond_indices().second;
	    int bond_index = mol.get_bond_index(ind_1, ind_2);
	    std::cout << "picked bond index: " << bond_index << std::endl;
	    pending_action_data_picked_bond_index = bond_index;
	    handled = true;
	    bond_pick_pending = false;
	    atom_pick_pending = true;
	 }
      }
   }

   if (atom_pick_pending) {
      if (highlight_data.has_contents()) {
	 if (highlight_data.single_atom()) {
	    int idx = highlight_data.get_atom_index();

	    RDKit::RWMol rdkm = rdkit_mol(mol);

	    RDKit::ROMol *main_sans_R
	       = coot::split_molecule(rdkm, pending_action_data_picked_bond_index, idx);

	    if (main_sans_R) {

	    }
	    delete main_sans_R;

	    handled = true;
	 }
      }
   }
   return handled;
}
#endif



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
   if (addition_mode == lbg_info_t::ATOM_H)
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
   if (addition_mode == lbg_info_t::ATOM_H)
      r = "H";
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
lbg_info_t::change_atom_id_maybe(unsigned int atom_index) {

   std::string ele = mol.atoms[atom_index].element;
   std::string fc = font_colour(ele);
   change_atom_element(atom_index, ele, fc);

}


bool
lbg_info_t::change_atom_element(unsigned int atom_index,
				std::string new_ele,
				std::string font_colour) {

   bool changed_status = 0;
   std::vector<unsigned int> local_bonds = mol.bonds_having_atom_with_atom_index(atom_index);
   lig_build::pos_t pos = mol.atoms[atom_index].atom_position;

   // 20110410: Old/simple
   // std::string atom_id = mol.make_atom_id_by_using_bonds(new_ele, local_bonds);

   bool gl_flag = false; // not a GL render engine
   lig_build::atom_id_info_t atom_id_info =
      mol.make_atom_id_by_using_bonds(atom_index, new_ele, local_bonds, gl_flag);
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

   changed_status = mol.atoms[atom_index].update_atom_id_maybe(atom_id_info,
							       font_colour, root);

   if (changed_status) {
      for (unsigned int ib=0; ib<local_bonds.size(); ib++) {
	 // make a new line for the bond
	 //
	 int index_1 = mol.bonds[local_bonds[ib]].get_atom_1_index();
	 int index_2 = mol.bonds[local_bonds[ib]].get_atom_2_index();
	 const lig_build::atom_t &atom_1 = mol.atoms[index_1];
	 const lig_build::atom_t &atom_2 = mol.atoms[index_2];

	 bool shorten_first  = false;
	 bool shorten_second = false;

	 std::pair<bool, bool> shorten = mol.shorten_flags(local_bonds[ib]);
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;

	 mol.bonds[local_bonds[ib]].make_new_canvas_item(atom_1, atom_2,
							 shorten.first, shorten.second,
							 empty, empty, root);
      }
   }
   return changed_status;
}


bool
lbg_info_t::try_add_or_modify_bond(int canvas_addition_mode,
				   int x_mouse, int y_mouse,
				   bool button_1_is_pressed) {


   // button_down_bond_addition is check later so that we can
   // distinguish between canvas dragging and new bond rotation.

   bool changed_status = 0;
   if (! highlight_data.has_contents()) {
      if (mol.is_empty()) {
	 // calling this sets penultimate_atom_pos and latest_bond_canvas_item
	 try_stamp_bond_anywhere(canvas_addition_mode, x_mouse, y_mouse);
	 changed_status = 1; // try_stamp_bond_anywhere always modifies
	 if (button_1_is_pressed)
	    button_down_bond_addition = 1;
      }
   } else {
      if (highlight_data.single_atom()) {
	 int atom_index = highlight_data.get_atom_index();
	 if (atom_index != UNASSIGNED_INDEX) {

	    // add_bond_to_atom() sets latest_bond_canvas_item
	    //
	    std::pair<bool, bool> changed_status_pair =
	       add_bond_to_atom(atom_index, canvas_addition_mode);
	    changed_status = changed_status_pair.first;

	    // we set button_down_bond_addition so that we can rotate
	    // the bond.  But we only want to (be able to) rotate the
	    // bond if button_1 is being pressed.
	    if (button_1_is_pressed)
	       button_down_bond_addition = changed_status;

	    if (changed_status_pair.second)
	       most_recent_bond_made_new_atom_flag = true;

	    if (false)
	       std::cout << "in try_add_or_modify_bond() 2 most_recent_bond_made_new_atom_flag "
			 << most_recent_bond_made_new_atom_flag << std::endl;

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

	    std::pair<bool, bool> shorten = mol.shorten_flags(bond_index);

	    // 20160909 the shortening is calculated outside of the bond
	    //          so I guess we don't need to pass the atoms.  Hmmm...
	    // // we need to pass the atoms so that we know if and how to
	    // // shorten the bonds (canvas items) to account for atom names.
	    //
	    lig_build::atom_t at_1 = mol.atoms[ind_1];
	    lig_build::atom_t at_2 = mol.atoms[ind_2];
	    GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

	    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
	       other_connections_to_second_atom =
	       mol.make_other_connections_to_second_atom_info(bond_index);
	    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
	       other_connections_to_first_atom =
	       mol.make_other_connections_to_first_atom_info(bond_index);

	    if (canvas_addition_mode == lbg_info_t::ADD_TRIPLE_BOND) {
	       mol.bonds[bond_index].change_bond_order(at_1, at_2, shorten.first, shorten.second, true,
						       other_connections_to_first_atom,
						       other_connections_to_second_atom,
						       root);
	    } else {

	       lig_build::bond_t::bond_type_t current_bt = mol.bonds[bond_index].get_bond_type();

	       // bt is the bond type we want to change to (or add)
	       //
	       // (bt is OUT_BOND if we want to convert an OUT_BOND to an IN_BOND (because it
	       // depends only on the "Stereo" button being depressed)).
	       //

	       if (bt == lig_build::bond_t::OUT_BOND) {

		  if (current_bt != lig_build::bond_t::OUT_BOND &&
		      current_bt != lig_build::bond_t::IN_BOND) {

		     // if it is not currently a stereo out bond, convert
		     // it to a stereo out bond. If it is a stereo_out
		     // bond, convert it with change_bond_order().
		     //
		     std::pair<bool, bool> shorten = mol.shorten_flags(bond_index);
		     mol.bonds[bond_index].set_bond_type(bt);
		     mol.bonds[bond_index].make_new_canvas_item(at_1, at_2,
								shorten.first, shorten.second,
								other_connections_to_first_atom,
								other_connections_to_second_atom,
								root);

		  } else {

		     if (mol.bonds[bond_index].get_bond_type() == lig_build::bond_t::IN_BOND) {
			highlight_data.swap_bond_indices(); // what does this do?
			std::swap(shorten.first, shorten.second);
		     }

		     mol.bonds[bond_index].change_bond_order(at_1, at_2,
							     shorten.first, shorten.second,
							     other_connections_to_first_atom,
							     other_connections_to_second_atom,
							     root); // out to in and vice
             		                                            // versa (with direction change)
		  }

	       } else {

		  // Not to a wedge bond

		  if (current_bt == lig_build::bond_t::OUT_BOND ||
		      current_bt == lig_build::bond_t::IN_BOND) {

		     // We want to change from wedge to single
		     mol.bonds[bond_index].set_bond_type(lig_build::bond_t::SINGLE_BOND);
		     std::pair<bool, bool> shorten = mol.shorten_flags(bond_index);
		     mol.bonds[bond_index].make_new_canvas_item(at_1, at_2,
								shorten.first, shorten.second,
								other_connections_to_first_atom,
								other_connections_to_second_atom,
								root);

		  } else {

		     // conventional/usual change
		     mol.bonds[bond_index].change_bond_order(at_1, at_2,
							     shorten.first, shorten.second,
							     other_connections_to_first_atom,
							     other_connections_to_second_atom,
							     root); // single to double
				  		                    // or visa versa
		  }
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

// return "was changed" and new-atom-created status pair
std::pair<bool, bool>
lbg_info_t::add_bond_to_atom(unsigned int atom_index, int canvas_addition_mode) {

   bool changed_status = false;
   bool created_new_atom_status = false;

   std::vector<unsigned int> bonds = mol.bonds_having_atom_with_atom_index(atom_index);

   switch (bonds.size()) {

   case 0:
      created_new_atom_status = add_bond_to_atom_with_0_neighbours(atom_index, canvas_addition_mode);
      changed_status = true;
      break;

   case 1:
      created_new_atom_status = add_bond_to_atom_with_1_neighbour(atom_index, canvas_addition_mode, bonds[0]);
      changed_status = true;
      break;

   case 2:
      created_new_atom_status = add_bond_to_atom_with_2_neighbours(atom_index, canvas_addition_mode, bonds);
      changed_status = true;
      break;

   case 3:
      created_new_atom_status = add_bond_to_atom_with_3_neighbours(atom_index, canvas_addition_mode, bonds);
      changed_status = true;
      break;

   default:
      std::cout << "not handled yet: " << bonds.size() << " neighbouring bonds" << std::endl;
   }
   // we might want to drag-rotate about atom_index shortly.
   if (changed_status)
      atom_index_of_atom_to_which_the_latest_bond_was_added = atom_index;

   return std::pair<bool, bool> (changed_status, created_new_atom_status);
}


std::string
lbg_info_t::font_colour(int addition_element_mode) const {
   std::string font_colour = "black";
   if (addition_element_mode == lbg_info_t::ATOM_O)
      font_colour = "red";
   if (addition_element_mode == lbg_info_t::ATOM_N)
      font_colour = "#0000aa";
   if (addition_element_mode == lbg_info_t::ATOM_S)
      font_colour = "#669900";
   if (addition_element_mode == lbg_info_t::ATOM_F)
      font_colour = "#00aa00";
   if (addition_element_mode == lbg_info_t::ATOM_CL)
      font_colour = "#229900";
   if (addition_element_mode == lbg_info_t::ATOM_BR)
      font_colour = "#A52A2A";
   if (addition_element_mode == lbg_info_t::ATOM_I)
      font_colour = "#220088";
   return font_colour;
}


std::string
lbg_info_t::font_colour(const std::string &ele) const {
   std::string font_colour = "#111111";

   if (ele == "N")
      font_colour = "#0000aa";
   if (ele == "O")
      font_colour = "red";
   if (ele == "S")
      font_colour = "#888800";
   if (ele == "P")
      font_colour = "#DD9500"; // orange
   if (ele == "F")
      font_colour = "#006600";
   if (ele == "CL")
      font_colour = "#116600";
   if (ele == "Cl") // mol files should have the second character in lower case, I think.
      font_colour = "#116600";
   if (ele == "Br")
      font_colour = "#A52A2A";
   if (ele == "I")
      font_colour = "#220066";

   return font_colour;
}

bool
lbg_info_t::add_bond_to_atom_with_0_neighbours(unsigned int atom_index, int canvas_addition_mode) {

   // certain change

   widgeted_atom_t atom = mol.atoms[atom_index];
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   lig_build::pos_t current_atom_pos = atom.atom_position;
   lig_build::pos_t a_bond(SINGLE_BOND_CANVAS_LENGTH, 0.0);
   a_bond.rotate(60);
   lig_build::pos_t new_atom_pos = current_atom_pos + a_bond;

   widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
   std::pair<bool, int> add_atom_status = mol.add_atom(new_atom);
   int new_index = add_atom_status.second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
   widgeted_bond_t b(atom_index, new_index, atom, new_atom, true, true, bt, empty, empty, root);
   mol.add_bond(b);

   change_atom_id_maybe(atom_index);
   change_atom_id_maybe(new_index);


   penultimate_atom_pos = atom.atom_position;
   penultimate_atom_index = atom_index;
   ultimate_atom_index    = new_index;

   return add_atom_status.first;
}


bool
lbg_info_t::add_bond_to_atom_with_1_neighbour(unsigned int atom_index, int canvas_addition_mode,
					      unsigned int bond_index) {

   unsigned int index_1 = mol.bonds[bond_index].get_atom_1_index();
   unsigned int index_2 = mol.bonds[bond_index].get_atom_2_index();

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

   std::vector<unsigned int> avoid_atoms(2);
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
   std::pair<bool, int> add_atom_status = mol.add_atom(new_atom);
   int new_index = add_atom_status.second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   bool shorten_first  = false;
   bool shorten_second = true; // CH3 superatom
   std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
   widgeted_bond_t b(atom_index, new_index, atom, new_atom, shorten_first, shorten_second, bt, empty, empty, root);
   mol.add_bond(b);

   change_atom_id_maybe(atom_index);
   change_atom_id_maybe(new_index);

   penultimate_atom_pos = atom.atom_position;
   penultimate_atom_index = atom_index;
   ultimate_atom_index    = new_index;

   return add_atom_status.second;
}

bool
lbg_info_t::add_bond_to_atom_with_2_neighbours(unsigned int atom_index,
					       int canvas_addition_mode,
					       const std::vector<unsigned int> &bond_indices) {

   widgeted_atom_t atom = mol.atoms[atom_index];
   unsigned int atom_index_1 = mol.bonds[bond_indices[0]].get_atom_1_index();
   unsigned int atom_index_2 = mol.bonds[bond_indices[0]].get_atom_2_index();
   unsigned int atom_index_3 = mol.bonds[bond_indices[1]].get_atom_1_index();
   unsigned int atom_index_4 = mol.bonds[bond_indices[1]].get_atom_2_index();

   unsigned int A_index = atom_index_1;
   if (atom_index_1 == atom_index)
      A_index = atom_index_2;
   unsigned int B_index = atom_index_3;
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
   std::pair<bool, int> checked_add = mol.add_atom(at);
   if (false) {
      if (checked_add.first)
	 std::cout << "debug:: created was a new atom: " << checked_add.second << std::endl;
      else
	 std::cout << "debug:: That connected to an already-existing atom: "
		   << checked_add.second << std::endl;
   }

   int new_index = checked_add.second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
   widgeted_bond_t b(atom_index, new_index, atom, at, false, true, bt, empty, empty, root);
   int new_bond_index = mol.add_bond(b);

   // Now, what about the atom that has been added to?  Now that we
   // have a extra bond, that atoms name could go from NH -> N (or C
   // -> C+ for 5 bonds).
   //
   change_atom_id_maybe(atom_index);
   change_atom_id_maybe(new_index);
   penultimate_atom_pos = atom.atom_position;
   penultimate_atom_index = atom_index;
   ultimate_atom_index    = new_index;
   return checked_add.first;
}


bool
lbg_info_t::add_bond_to_atom_with_3_neighbours(unsigned int atom_index,
					       int canvas_addition_mode,
					       const std::vector<unsigned int> &bond_indices) {

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   // are (at least) 2 of the bonds attached to atom_index terminated
   // at their end?
   //
   // If so, we shall remove and replace those bonds before we add
   // this new one.
   //
   std::pair<bool, std::vector<unsigned int> > pr = have_2_stubs_attached_to_atom(atom_index, bond_indices);
   std::vector<unsigned int> attached_bonds = pr.second;

   if (pr.first) {
      unsigned int l = attached_bonds.size();

      widgeted_bond_t bond_to_core = orthogonalise_2_bonds(atom_index, attached_bonds, bond_indices);
      // now add a third

      int p1_idx = bond_to_core.get_other_index(atom_index);
      lig_build::pos_t p1 = mol.atoms[p1_idx].atom_position;

      lig_build::pos_t central_atom_pos = mol.atoms[atom_index].atom_position;
      lig_build::pos_t existing_bond_dir = central_atom_pos - p1;
      lig_build::pos_t ebd_uv = existing_bond_dir.unit_vector();

      lig_build::pos_t new_atom_pos = central_atom_pos + ebd_uv * SINGLE_BOND_CANVAS_LENGTH;
      widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
      std::pair<bool, int> add_atom_status = mol.add_atom(new_atom);
      int new_atom_index = add_atom_status.second;
      lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
      bool shorten_first  = false;
      bool shorten_second = false;
      if (mol.atoms[atom_index].element != "C") shorten_first = true;
      std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
      widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], new_atom,
			shorten_first, shorten_second, bt, empty, empty, root);
      mol.add_bond(b);
      change_atom_id_maybe(atom_index);
      change_atom_id_maybe(new_atom_index);
      return add_atom_status.first;

   } else {
      // bond_indices are bonds have an atom that is atom_index.
      squeeze_in_a_4th_bond(atom_index, canvas_addition_mode, bond_indices);
      return false; // this virtually never happens to be true, but fix one day
   }
}

void
lbg_info_t::squeeze_in_a_4th_bond(unsigned int atom_index, int canvas_addition_mode,
				  const std::vector<unsigned int> &bond_indices) {

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
	 bool shorten_first  = false;
	 bool shorten_second = false;
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
	 widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index],
			   new_atom, shorten_first, shorten_second, bt, empty, empty, root);
	 mol.add_bond(b);
	 change_atom_id_maybe(atom_index);
	 change_atom_id_maybe(new_atom_index);
	 ultimate_atom_index    = new_atom_index;

      } else {

	 // OK, there were not centres on all side of all the bonds.
	 // Go place the new atom not towards a ring centre

	 lig_build::pos_t new_atom_pos = get_new_pos_not_towards_ring_centres(atom_index,
									      bond_indices);
	 widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
	 int new_atom_index = mol.add_atom(new_atom).second;
	 lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
	 bool shorten_first  = false;
	 bool shorten_second = false;
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
	 widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], new_atom,
			   shorten_first, shorten_second, bt, empty, empty, root);
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
      bool shorten_first  = false;
      bool shorten_second = false;
      std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
      widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], new_atom,
			shorten_first, shorten_second, bt, empty, empty, root);
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
lbg_info_t::all_closed_rings(unsigned int atom_index,
			     const std::vector<unsigned int> &bond_indices) const {

   bool status = false;

   if (bond_indices.size() > 3) {
      std::vector<lig_build::pos_t> centres = get_centres_from_bond_indices(bond_indices);
      if (centres.size() > 2)
	 status = true;
   }
   return status;
}

std::vector<lig_build::pos_t>
lbg_info_t::get_centres_from_bond_indices(const std::vector<unsigned int> &bond_indices) const {

   std::vector<lig_build::pos_t> centres;
   for (unsigned int ib=0; ib<bond_indices.size(); ib++) {
      if (mol.bonds[bond_indices[ib]].have_centre_pos()) {
	 lig_build::pos_t test_centre = mol.bonds[bond_indices[ib]].centre_pos();

	 // add test_centre to centres only if it wasnt there already.
	 bool found_centre = false;
	 for (unsigned int j=0; j<centres.size(); j++) {
	    if (test_centre.close_point(centres[j])) {
	       found_centre = true; // it was already there
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
lbg_info_t::get_new_pos_not_towards_ring_centres(unsigned int atom_index,
						 const std::vector<unsigned int> &bond_indices) const {

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

      if (centres.size() > 0) {
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
      } else {

	 std::vector<lig_build::pos_t> neighbours;
	 p = lig_build::pos_t(100,100);
	 for (unsigned int i=0; i<bond_indices.size(); i++) {
	    int other_index = mol.bonds[bond_indices[i]].get_other_index(atom_index);
	    lig_build::pos_t p1 = mol.atoms[other_index].atom_position;
	    neighbours.push_back(p1);
	 }

	 std::pair<unsigned int, unsigned int> best_pair(0,0);
	 double delta_length = 0;
	 if (neighbours.size() > 0) {
	    for (unsigned int i=0; i<(neighbours.size()-1); i++) {
	       lig_build::pos_t delta = neighbours[i+1] - neighbours[i];
	       double len = delta.length();
	       if (len > delta_length) {
		  delta_length = len;
		  best_pair.first  = i+1;
		  best_pair.second = i;
	       }
	    }

	    // now additionally test last to first
	    lig_build::pos_t delta = neighbours[0] - neighbours.back();
	    double len = delta.length();
	    if (len > delta_length) {
	       best_pair.first  = 0;
	       best_pair.second = neighbours.size()-1; // the last one
	    }

	    if (delta_length > 0) {
	       lig_build::pos_t b1 = neighbours[best_pair.second] + neighbours[best_pair.first];
	       lig_build::pos_t b2(0.5*b1.x, 0.5*b1.y);
	       lig_build::pos_t bond_dir = b2 - centre;
	       lig_build::pos_t bond_dir_uv = bond_dir.unit_vector();
	       p = centre + bond_dir_uv * SINGLE_BOND_CANVAS_LENGTH * 0.8;
	    }
	 }
      }
   }
   return p;
}


lig_build::pos_t
lbg_info_t::new_pos_by_bisection(unsigned int atom_index,
				 const std::vector<unsigned int> &bond_indices,
				 const std::vector<double> &angles,
				 GooCanvasItem *root) const {

   lig_build::pos_t p(0,0);
   std::vector<double> sorted_angles = angles;
   std::sort(sorted_angles.begin(), sorted_angles.end());
   lig_build::pos_t centre = mol.atoms[atom_index].atom_position;
   for (unsigned int i=0; i<bond_indices.size(); i++) {
      unsigned int j = i+1;
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
lbg_info_t::get_angles(unsigned int atom_index, const std::vector<unsigned int> &bond_indices) const {

   std::vector<double> v(bond_indices.size());

   lig_build::pos_t centre = mol.atoms[atom_index].atom_position;

   for (unsigned int i=0; i<bond_indices.size(); i++) {
      unsigned int j = i+1;
      if (j == bond_indices.size())
	 j = 0;
      unsigned int idx_1 = mol.bonds[bond_indices[i]].get_other_index(atom_index);
      unsigned int idx_2 = mol.bonds[bond_indices[j]].get_other_index(atom_index);
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
lbg_info_t::orthogonalise_2_bonds(unsigned int atom_index,
				  const std::vector<unsigned int> &stub_attached_atoms,
				  const std::vector<unsigned int> &bond_indices) {

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
   bool shorten_first  = false;
   bool shorten_second = false;
   if (mol.atoms[atom_index].element != "C") shorten_first = true;
   if (atom_2.element != "C") shorten_second = true;
   std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
   widgeted_bond_t b_1(atom_index, new_atom_index_1, mol.atoms[atom_index], atom_1, shorten_first, shorten_second, bt, empty, empty, root);
   widgeted_bond_t b_2(atom_index, new_atom_index_2, mol.atoms[atom_index], atom_2, shorten_first, shorten_second, bt, empty, empty, root);
   mol.add_bond(b_1);
   mol.add_bond(b_2);

   return bond_to_core; // so that we don't have to recalculate it in caller.

}



// return a status and a vector of atoms (bonded to atom_index) having
// only one bond.
//
std::pair<bool, std::vector<unsigned int> >
lbg_info_t::have_2_stubs_attached_to_atom(unsigned int atom_index,
					  const std::vector<unsigned int> &bond_indices) const {

   std::vector<unsigned int> v;
   for (unsigned int i=0; i<bond_indices.size(); i++) {
      int other_index = mol.bonds[bond_indices[i]].get_other_index(atom_index);
      // now, does other_index have only one bond?
      std::vector<unsigned int> local_bonds = mol.bonds_having_atom_with_atom_index(other_index);
      if (local_bonds.size() == 1) {
	 v.push_back(other_index);
      }
   }
   bool status = 0;
   if (v.size() > 1)
      status = 1;
   return std::pair<bool, std::vector<unsigned int> > (status, v);
}



// Set the addition mode and turn off other toggle buttons (if
// needed).
void
lbg_handle_toggle_button(GtkToggleToolButton *tb, GtkWidget *canvas, int mode) {

   gpointer gp = g_object_get_data(G_OBJECT(canvas), "lbg");
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
			bool shorten_first  = false;
			bool shorten_second = false;
			std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
			widgeted_bond_t bond(highlighted_atom_index, new_atoms[0],
					     mol.atoms[highlighted_atom_index],
					     mol.atoms[new_atoms[0]],
					     shorten_first, shorten_second,
					     bt, empty, empty, root);
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

bool
lbg_info_t::handle_charge_change() {

   bool changed_status = false;
   if (highlight_data.has_contents()) {
      if (highlight_data.single_atom()) {
	 int atom_index = highlight_data.get_atom_index();
	 if (atom_index != UNASSIGNED_INDEX) {
	    int charge = mol.atoms[atom_index].charge;
	    int pre_charge = charge;
	    if (charge >= 2) {
	       charge = -2;
	    } else {
	       charge += 1;
	    }
	    mol.atoms[atom_index].charge = charge;
	    if (false)
	       std::cout << "change charge on " << mol.atoms[atom_index]
			 << " from " << pre_charge << " to " << charge
			 << std::endl;

	    // 20160701:
	    // prep for calling update_atom_id_forced():
	    // no longer change_atom_id_maybe(atom_index);
	    //
	    GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
	    std::string ele = mol.atoms[atom_index].element;
	    std::string fc = font_colour(ele);
	    bool gl_flag = false; // not a GL render engine
	    std::vector<unsigned int> local_bonds = mol.bonds_having_atom_with_atom_index(atom_index);
	    lig_build::atom_id_info_t atom_id_info =
	       mol.make_atom_id_by_using_bonds(atom_index, ele, local_bonds, gl_flag);

	    mol.atoms[atom_index].update_atom_id_forced(atom_id_info, fc, root);
	    changed_status = true;
	 }
      }
   }

  // update_descriptor_attributes() is done by calling function (e.g. handle_item_add()).

  return changed_status;
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
				 "line-width", 7.0, // in stamp_polygon_anywhere() (debug)
				 "stroke-color", stroke_colour.c_str(),
				 NULL);
      }
   }

   for (unsigned int i=0; i<atom_pos_index.size(); i++) {
      unsigned int j = i + 1;
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
	 bool shorten_first  = false;
	 bool shorten_second = false;
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
	 widgeted_bond_t bond(i_index, j_index,
			      mol.atoms[atom_pos_index[i].second.second],
			      mol.atoms[atom_pos_index[j].second.second],
			      shorten_first, shorten_second,
			      centre, bt, empty, empty, root);
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
   lig_build::pos_t c_2 = m - mp_unit * cor_radius;

   lig_build::pos_t c = c_1;
   bool swapped_centre_direction = false;

   if (has_ring_centre_flag) {
      double l1 = lig_build::pos_t::length(c_1, ring_centre);
      double l2 = lig_build::pos_t::length(c_2, ring_centre);
      if (l1 < l2) {
	 c = c_2;
	 swapped_centre_direction = true;
      }
   }

   //  0    for a triangle : 120
   //  45   for a square   : 90
   //  0    for pentagon   : 72
   // 30    for hexagon    : 60
   // 51.43 for pentagon   : 51.43
   // 22.5  for octagon    : 45

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
   if (swapped_centre_direction)
      angle_off -= 180;

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

   unsigned int atom_index = get_atom_index();
   std::vector<unsigned int> bv = mol.bonds_having_atom_with_atom_index(atom_index);

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
   w_names.push_back("lbg_charge_toggle_toolbutton");

   // undo and clear, cut, smiles

   for (unsigned int i=0; i<w_names.size(); i++) {
      GtkToggleToolButton *tb =
	 GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object(builder, w_names[i].c_str()));
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
lbg_info_t::clear(bool do_descriptor_updates) {

   clear_canvas();
   // and that clears the alerts group
   alert_group = NULL;

   key_group = NULL;

   // clear the molecule
   mol.clear();

   if (do_descriptor_updates) {
      update_descriptor_attributes();
      update_apply_button_sensitivity();
   }
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
	 lbg_export_as_svg_dialog = NULL;
	 lbg_sbase_search_results_dialog = NULL;
	 lbg_sbase_search_results_vbox = NULL;
	 lbg_smiles_dialog = NULL;
	 lbg_smiles_entry = NULL;
	 lbg_statusbar = NULL;
	 lbg_toolbar_layout_info_label = NULL;
	 lbg_atom_x_dialog = NULL;
	 lbg_atom_x_entry = NULL;
	 lbg_clean_up_2d_toolbutton = NULL;
	 lbg_search_database_frame = NULL;
	 lbg_import_from_smiles_dialog = NULL;
	 lbg_import_from_smiles_entry = NULL;
	 lbg_view_rotate_entry = NULL;
	 lbg_get_drug_menuitem = NULL;
	 lbg_scale_spinbutton = NULL;
	 for (unsigned int i=0; i<8; i++)
	    lbg_qed_properties_progressbars[i] = NULL;
	 lbg_srs_search_results_scrolledwindow = NULL;
	 lbg_srs_search_results_vbox = NULL;

	 canvas = NULL;
	 return false; // boo.

      } else {

	 gtk_widget_show (lbg_window);
	 about_dialog =    GTK_WIDGET (gtk_builder_get_object (builder, "lbg_aboutdialog"));
	 open_dialog     = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_open_filechooserdialog"));
	 save_as_dialog  = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_save_as_filechooserdialog"));
         lbg_apply_button = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_apply_button"));
	 lbg_sbase_search_results_dialog = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_sbase_search_results_dialog"));
	 lbg_sbase_search_results_vbox = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_sbase_search_results_vbox"));
	 lbg_export_as_pdf_dialog =      GTK_WIDGET (gtk_builder_get_object (builder, "lbg_export_as_pdf_filechooserdialog"));
	 lbg_export_as_png_dialog =      GTK_WIDGET (gtk_builder_get_object (builder, "lbg_export_as_png_filechooserdialog"));
	 lbg_export_as_svg_dialog =      GTK_WIDGET (gtk_builder_get_object (builder, "lbg_export_as_svg_filechooserdialog"));
	 lbg_smiles_dialog =             GTK_WIDGET(gtk_builder_get_object(builder, "lbg_smiles_dialog"));
	 lbg_smiles_entry =              GTK_WIDGET(gtk_builder_get_object(builder, "lbg_smiles_entry"));
	 lbg_import_from_smiles_dialog = GTK_WIDGET(gtk_builder_get_object(builder, "lbg_import_from_smiles_dialog"));
	 lbg_import_from_smiles_entry  = GTK_WIDGET(gtk_builder_get_object(builder, "lbg_import_from_smiles_entry"));
	 lbg_import_from_comp_id_dialog = GTK_WIDGET(gtk_builder_get_object(builder, "lbg_import_from_comp_id_dialog"));
	 lbg_import_from_comp_id_entry  = GTK_WIDGET(gtk_builder_get_object(builder, "lbg_import_from_comp_id_entry"));
	 lbg_import_from_comp_id_hydrogens_checkbutton  = GTK_WIDGET(gtk_builder_get_object(builder, "lbg_import_from_comp_id_hydrogens_checkbutton"));
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
	 lbg_qed_properties_vbox =       GTK_WIDGET(gtk_builder_get_object(builder, "lbg_qed_properties_vbox"));
	 lbg_alert_name_label =          GTK_WIDGET(gtk_builder_get_object(builder, "lbg_alert_name_label"));
	 lbg_get_drug_dialog =           GTK_WIDGET(gtk_builder_get_object(builder, "lbg_get_drug_dialog"));
	 lbg_get_drug_entry =            GTK_WIDGET(gtk_builder_get_object(builder, "lbg_get_drug_entry"));
	 lbg_get_drug_menuitem =         GTK_WIDGET(gtk_builder_get_object(builder, "lbg_get_drug_menuitem"));
	 pe_test_function_button =       GTK_WIDGET(gtk_builder_get_object(builder, "pe_test_function_button"));
	 lbg_flip_rotate_hbox =          GTK_WIDGET(gtk_builder_get_object(builder, "lbg_flip_rotate_hbox"));
	 lbg_clean_up_2d_toolbutton =    GTK_WIDGET(gtk_builder_get_object(builder, "lbg_clean_up_2d_toolbutton"));
	 lbg_search_database_frame =     GTK_WIDGET(gtk_builder_get_object(builder, "lbg_search_database_frame"));
	 lbg_view_rotate_entry     =     GTK_WIDGET(gtk_builder_get_object(builder, "lbg_view_rotate_entry"));
	 lbg_scale_spinbutton      =     GTK_WIDGET(gtk_builder_get_object(builder, "lbg_scale_spinbutton"));
	 lbg_srs_search_results_vbox =   GTK_WIDGET(gtk_builder_get_object(builder, "lbg_srs_search_results_vbox"));
	 lbg_srs_search_results_scrolledwindow = GTK_WIDGET(gtk_builder_get_object(builder, "lbg_srs_search_results_scrolledwindow"));

	 for (unsigned int i=0; i<8; i++) {
	    std::string name = "qed_properties_" + coot::util::int_to_string(i) + "_progressbar";
	    lbg_qed_properties_progressbars[i] = GTK_WIDGET(gtk_builder_get_object(builder, name.c_str()));
	 }

	 gtk_label_set_text(GTK_LABEL(lbg_toolbar_layout_info_label), "---");

      }
   }

   canvas = goo_canvas_new();

   // this works, but need gui access
   //
   goo_canvas_set_scale(GOO_CANVAS(canvas), 1.0);

   // gtk_widget_set_flags(canvas, GTK_CAN_FOCUS);
   std::cout << "FIXME gtk_widget_set_flags GTK_CAN_FOCUS!"  << std::endl;

   GooCanvasItem *root_item = goo_canvas_get_root_item(GOO_CANVAS(canvas));
   g_object_set(G_OBJECT(root_item), "line_width", 1.5, NULL); // thank you Damon Chaplin
   g_object_set(G_OBJECT(root_item), "line_cap", CAIRO_LINE_CAP_ROUND, NULL);

   if (lbg_scale_spinbutton) {
                                     // value lower upper step_inc page_inc page_size (which should be 0)
      GtkAdjustment *adj = gtk_adjustment_new(1.0, 0.2, 12.0, 0.2, 0.2, 0.0);
      gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(lbg_scale_spinbutton), GTK_ADJUSTMENT(adj));
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(lbg_scale_spinbutton), 1.0);

      g_object_set_data(G_OBJECT(lbg_scale_spinbutton), "user_data", canvas);
      g_signal_connect(adj, "value-changed", G_CALLBACK(lbg_scale_adj_changed), lbg_scale_spinbutton);

   }


#ifdef HAVE_CCP4SRS
   gtk_widget_show(lbg_search_database_frame);
   gtk_widget_hide(lbg_srs_search_results_scrolledwindow); // show this when SRS results
#else
   // ... we don't have ccp4 srs
   gtk_widget_hide(lbg_srs_search_results_scrolledwindow);
#endif // HAVE_CCP4SRS

   if (use_graphics_interface_flag) {
      // gtk_widget_set(GTK_WIDGET(canvas), "bounds-padding", 50.0, NULL);
      std::cout << "FIXME! bounds padding!" << std::endl;
      g_object_set_data(G_OBJECT(canvas), "lbg", (gpointer) this);

      save_togglebutton_widgets(builder);
      GtkWidget *lbg_scrolled_win =
	 GTK_WIDGET(gtk_builder_get_object (builder, "lbg_scrolledwindow"));
      gtk_container_add(GTK_CONTAINER(lbg_scrolled_win), canvas);
      goo_canvas_set_bounds (GOO_CANVAS (canvas), 0, 0, 1200, 1200);
   }

   GooCanvas *gc = GOO_CANVAS(canvas);

#if 0
   if (0)
      std::cout << "   after setting bounds: canvas adjustments: h-lower "
		<< gc->hadjustment->lower << " h-upper-page "
		<< gc->hadjustment->upper - gc->hadjustment->page_size  << " v-lower "
		<< gc->vadjustment->lower  << " v-upper-page"
		<< gc->vadjustment->upper - gc->vadjustment->page_size  << " h-page "
		<< gc->hadjustment->page_size  << " "
		<< gc->vadjustment->page_size << std::endl;
#endif

   if (use_graphics_interface_flag) {

      gtk_widget_add_events(GTK_WIDGET(canvas),
 			    GDK_KEY_PRESS_MASK     |
 			    GDK_KEY_RELEASE_MASK);

      GooCanvasItem *root_item = goo_canvas_get_root_item(GOO_CANVAS(canvas));

      g_object_set_data_full(G_OBJECT(root_item), "lbg-info", this, NULL);
      g_signal_connect(G_OBJECT(root_item), "button_press_event",
		       G_CALLBACK(on_canvas_button_press), NULL);
      g_signal_connect(G_OBJECT(root_item), "button_release_event",
		       G_CALLBACK(on_canvas_button_release), NULL);
      g_signal_connect(G_OBJECT(root_item), "motion_notify_event",
		       G_CALLBACK(on_canvas_motion), NULL);

      // setup_lbg_drag_and_drop(lbg_window); // this lbg_info_t is not attached to this object.
      setup_lbg_drag_and_drop(canvas);
      gtk_widget_show(canvas);
      add_search_combobox_text();
   }


   // ------------ (old-style) watch for new files from coot ---------------------
   //
   // Gtk3 porting: comment this out for now.
   // if (is_stand_alone())
   // int timeout_handle = gtk_timeout_add(500, watch_for_mdl_from_coot, this);

   // Hack in a button (or hide the hbox) for PE to test stuff
   //
   if (use_graphics_interface_flag) {
      if (getenv("COOT_LBG_TEST_FUNCTION") != NULL) {
         gtk_widget_show(pe_test_function_button);
      }
   }

   // if we don't have rdkit or python then we don't want to see qed progress bar
   // or the "show alerts" (because we can't match to the alert patterns).

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#ifdef USE_PYTHON
   // all, with QED

   if (! silicos_it_qed_default_func) { // set in init
      gtk_widget_hide(lbg_qed_hbox);
      gtk_widget_hide(lbg_qed_properties_vbox);
   }

#else

   gtk_widget_hide(lbg_qed_hbox);
   gtk_widget_hide(lbg_qed_properties_vbox);

#endif
#else
   gtk_widget_hide(lbg_qed_hbox);
   gtk_widget_hide(lbg_alert_hbox_outer);
   gtk_widget_hide(lbg_show_alerts_checkbutton); // perhaps this should be in the
                                                 // lbg_alert_hbox_outer?
   gtk_widget_hide(lbg_qed_properties_vbox);
   gtk_widget_hide(lbg_clean_up_2d_toolbutton);
#endif

   return true;

}
#endif // GTK_VERSION


// void lbg_scale_adj_changed(GtkWidget *widget, GtkSpinButton *spinbutton)

// extern "C" G_MODULE_EXPORT
void
lbg_scale_adj_changed(GtkWidget *widget, GtkSpinButton *spinbutton) {

   float f = 1.0; // gtk_spin_button_get_value_as_float(spinbutton);
   std::cout << "FIXME why oh why is gtk_spin_button_get_value_as_float() missing now!?"
             << std::endl;
   gpointer user_data = g_object_get_data(G_OBJECT(spinbutton), "user_data");
   if (user_data) {
      GtkWidget *canvas = GTK_WIDGET(user_data);
      lbg_info_t *l = static_cast<lbg_info_t *> (g_object_get_data(G_OBJECT(canvas), "lbg"));
      if (l) {
	 l->scale_canvas(f);
      }
   }
}



void
lbg_info_t::scale_canvas(double sf) {

   // std::cout << "setting scale " << sf << std::endl;
   if (sf < 0.1) sf = 0.1;
   goo_canvas_set_scale(GOO_CANVAS(canvas), sf);

}


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

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   std::string s;
   try {
      s = get_smiles_string_from_mol_rdkit();
   }
   catch (const std::exception &e) {
      std::cout << "get_smiles_string_from_mol() " << e.what() << std::endl;
   }
   return s;
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

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

void
lbg_info_t::update_qed(const RDKit::RWMol &rdkm) {

   // see setup_silicos_it_qed_default_func()

#ifdef USE_PYTHON

   if (rdkm.getNumAtoms() == 0) {
      // non-interesting case first
      gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(lbg_qed_progressbar), 0);
      gtk_label_set_text(GTK_LABEL(lbg_qed_text_label), "");
      std::vector<std::pair<double, double> > dummy; // resets progressbars to 0
      update_qed_properties(dummy);
   } else {
      bool all_set = false;
      double qed = 0.0;
      if (silicos_it_qed_default_func) {
	 qed = get_qed(silicos_it_qed_default_func, rdkm);
	 if (qed > 0)
	    all_set = true;

	 // get the values and their desirabilites (0->1)
	 std::vector<std::pair<double, double> > properties =
	    get_qed_properties(silicos_it_qed_properties_func,
			       silicos_it_qed_pads, rdkm);
	 update_qed_properties(properties);

      } else {

	 // If you are reading this: are you sure that Biscu-it has
	 // been installed - and works properly in just python?

	 // std::cout << "null silicos_it_qed_default_func " << std::endl;
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

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
void
lbg_info_t::update_qed_properties(const std::vector<std::pair<double, double> > &properties) {

   if (properties.size() == 8) {
      std::vector<std::string> ideal_values(8);
      ideal_values[0] = "300";
      ideal_values[1] = "3";
      ideal_values[2] = "3"; // acceptors
      ideal_values[3] = "1"; // donors
      ideal_values[4] = "55"; // PSA
      ideal_values[5] = "4"; // rotatable bonds
      ideal_values[6] = "2"; // arom
      ideal_values[7] = "0"; // alerts
      for (unsigned int i=0; i<8; i++) {
	 if (false)
	    std::cout << "desirability " << i << " "
		      << properties[i].first << " "
		      << properties[i].second << " "
		      << lbg_qed_properties_progressbars[i] << std::endl;
	 if (properties[i].second >= 0) {
	    if (properties[i].second <= 1) {
	       std::string s;
	       if (i == 2 || i == 3 || i == 5 || i == 6 || i == 7)
		  s = coot::util::int_to_string(int(properties[i].first));
	       else
		  s = coot::util::float_to_string(properties[i].first);

	       gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(lbg_qed_properties_progressbars[i]),
					     properties[i].second);
	       gtk_progress_bar_set_text(GTK_PROGRESS_BAR(lbg_qed_properties_progressbars[i]),
					 s.c_str());

	       // tooltips here?
	       GtkTooltips *tooltips = gtk_tooltips_new();
	       std::string m = "Ideal: " + ideal_values[i];
	       gtk_tooltips_set_tip(tooltips, GTK_WIDGET(lbg_qed_properties_progressbars[i]), m.c_str(), NULL);
	    }
	 }
      }
   } else {
      reset_qed_properties_progress_bars();
   }
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
void
lbg_info_t::reset_qed_properties_progress_bars() {

   for (unsigned int i=0; i<8; i++) {
      gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(lbg_qed_properties_progressbars[i]), 0);
      gtk_progress_bar_set_text(GTK_PROGRESS_BAR(lbg_qed_properties_progressbars[i]), "");
   }
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
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
					       "fill-color-rgba", 0xffbb3350,
					       NULL);

		  }
		  catch (const KeyErrorException &kee) {
		     // actually, perhaps we don't get this we should fail, else we risk a crash?
		     // lbg_atom_index = rdkmol_idx;
		     std::cout << "ERROR:: in update_alerts() failed to get atom name "
			       << kee.what() << std::endl;
		  }
	       }
	    }
	    goo_canvas_item_lower(alert_group, NULL); // to the bottom
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

      // g_object_set(alert_group, "visibility", GOO_CANVAS_ITEM_INVISIBLE, NULL);

      gint n_children = goo_canvas_item_get_n_children (alert_group);
      // std::cout << "removing " << n_children << " children " << std::endl;
      for (int i=0; i<n_children; i++) {
	 // no, don't delete it because if we delete it then we don't
	 // get a button up signal when button-1 is released
	 //
         // goo_canvas_item_remove_child(alert_group, 0);
	 GooCanvasItem *child = goo_canvas_item_get_child(alert_group, i);
	 g_object_set(child, "visibility", GOO_CANVAS_ITEM_INVISIBLE, NULL);
      }
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

   mmdb::Manager *flat_pdb_mol = get_cmmdbmanager(coot_pdb_file);
   lig_build::molfile_molecule_t mm;
   mm.read(coot_mdl_file);
   widgeted_molecule_t wmol = import_mol_file(mm, coot_mdl_file, flat_pdb_mol);
   std::vector<solvent_accessible_atom_t> solvent_accessible_atoms =
      read_solvent_accessibilities(sa_file);
   wmol.map_solvent_accessibilities_to_atoms(solvent_accessible_atoms);

   import_from_widgeted_molecule(wmol);
   render();

}


mmdb::Manager *
lbg_info_t::get_cmmdbmanager(const std::string &file_name) const {

   mmdb::Manager *mol = new mmdb::Manager;
   mmdb::ERROR_CODE err = mol->ReadCoorFile(file_name.c_str());
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

// update the internal class variable widgeted_molecule_t mol from mol_in
//
void
lbg_info_t::import_from_widgeted_molecule(const widgeted_molecule_t &mol_in) {

   // updates mol

   make_saves_mutex = 0; // stop saving changes (restored at end)
   clear(false); // clear and no descriptor updates
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   int re_index[mol_in.atoms.size()]; // map from mol_in atom indexing
				      // the this molecule atom
				      // indexing.

   if (false)
      std::cout << "import_from_widgeted_molecule with " << mol_in.atoms.size()
		<< " atoms" << std::endl;

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
	 if (false)
	    std::cout << "render from molecule " << iat
		      << " id: " << mol_in.atoms[iat].atom_id
		      << " name: " << mol_in.atoms[iat].atom_name
		      << " charge:" << mol_in.atoms[iat].charge
		      << " pos: " << mol_in.atoms[iat].atom_position
		      << std::endl;

	 double sa = mol_in.atoms[iat].get_solvent_accessibility();
	 new_atom.set_atom_name(mol_in.atoms[iat].get_atom_name());

	 if (0)
	    std::cout << "in render_from_molecule() atom with name :"
		      << mol_in.atoms[iat].get_atom_name()
		      << ": has solvent_accessibility " << sa << std::endl;

	 new_atom.add_solvent_accessibility(sa);
	 new_atom.bash_distances = mol_in.atoms[iat].bash_distances;
	 std::pair<bool, int> s = mol.add_atom(new_atom);
	 re_index[iat] = s.second;

	 if (false)
	    if (! s.first)
	       std::cout << "render_from_molecule() atom " << iat
			 << " was close to atom " << s.second << " "
			 << mol_in.atoms[iat].atom_position << "  vs.  "
			 << mol.atoms[s.second].atom_position
			 << std::endl;

	 if (0)
	    // clever c++
	    std::cout << " in render_from_molecule: old sa: "
		      << mol_in.atoms[iat].get_solvent_accessibility() << " new: "
		      << mol.atoms[re_index[iat]].get_solvent_accessibility() << std::endl;
      }
   }

   // add in bonds
   for (unsigned int ib=0; ib<mol_in.bonds.size(); ib++) {
      if (! mol_in.bonds[ib].is_closed()) {
	 int idx_1 = re_index[mol_in.bonds[ib].get_atom_1_index()];
	 int idx_2 = re_index[mol_in.bonds[ib].get_atom_2_index()];
	 if ((idx_1 != UNASSIGNED_INDEX) && (idx_2 != UNASSIGNED_INDEX)) {
	    if (idx_1 != idx_2) { // 20130111 just in case there are problems from re-indexing.
	                          // better to not have the bond than to have and atom bonded to itself.
	       lig_build::bond_t::bond_type_t bt = mol_in.bonds[ib].get_bond_type();
	       if (mol_in.bonds[ib].have_centre_pos()) {

		  std::pair<bool, bool> shorten = mol.shorten_flags(ib);
		  if (display_atom_names || display_atom_numbers) {
		     shorten.first = true;
		     shorten.second = true;
		  }
		  lig_build::pos_t centre_pos = mol_in.bonds[ib].centre_pos();
		  std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
		  widgeted_bond_t bond(idx_1, idx_2, mol.atoms[idx_1], mol.atoms[idx_2],
				       shorten.first, shorten.second,
				       centre_pos, bt, empty, empty, root);
		  mol.add_bond(bond);

		  if (false) { // debug ring centres
		     if (bt == lig_build::bond_t::DOUBLE_BOND) {
			lig_build::pos_t pos = centre_pos;
			pos += mol_in.atoms[mol_in.bonds[ib].get_atom_1_index()].atom_position;
			pos += mol_in.atoms[mol_in.bonds[ib].get_atom_2_index()].atom_position;
			pos = pos * 0.333333;
			GooCanvasItem *item =
			   goo_canvas_ellipse_new(root,
						  pos.x, pos.y,
						  6.0, 6.0,
						  "line-width", 1.0,
						  "stroke-color-rgba", 0xffffffaa,
						  "fill_color_rgba", 0x992299aa,
						  NULL);
		     }
		  }


	       } else {
		  // bond with no ring centre
		  bool shorten_first  = false;
		  bool shorten_second = false;
		  if (mol.atoms[idx_1].element != "C")
		     shorten_first = true;
		  if (mol.atoms[idx_2].element != "C")
		     shorten_second = true;

		  if (display_atom_names || display_atom_numbers) {
		     shorten_first  = true;
		     shorten_second = true;
		  }

		  std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
		  widgeted_bond_t bond(idx_1, idx_2, mol.atoms[idx_1], mol.atoms[idx_2],
				       shorten_first, shorten_second, bt, empty, empty, root);
		  mol.add_bond(bond);
	       }
	    }
	 }
      }
   }

   // for input_coords_to_canvas_coords() to work:
   //
   mol.centre_correction = mol_in.centre_correction;
   mol.scale_correction  = mol_in.scale_correction;
   mol.mol_in_min_y = mol_in.mol_in_min_y;
   mol.mol_in_max_y = mol_in.mol_in_max_y;

   //
   make_saves_mutex = 1; // allow saves again.

}

void
lbg_info_t::render() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   if (false)
      std::cout << "------------ render_from_molecule() with display_atom_names " << display_atom_names
		<< " display_atom_numbers "  << display_atom_numbers << " "
		<< mol.atoms.size() << " atoms " << mol.bonds.size() << " bonds" << std::endl;

   if (! display_atom_names && ! display_atom_numbers) {

      // shorten C superatom bonds (we have to do this after all the bonds have been added)
      //
      for (unsigned int ibond=0; ibond<mol.bonds.size(); ibond++) {
	 int idx_1 = mol.bonds[ibond].get_atom_1_index();
	 int idx_2 = mol.bonds[ibond].get_atom_2_index();
	 std::pair<bool, bool> shorten = mol.shorten_flags(ibond);
	 // why do we pass the atom_t if we lose it?
	 const lig_build::atom_t &at_1 = mol.atoms[idx_1];
	 const lig_build::atom_t &at_2 = mol.atoms[idx_2];
	 if (shorten.first || shorten.second) {

	    if (mol.bonds[ibond].get_bond_type() == lig_build::bond_t::OUT_BOND) {

	       if (shorten.second) {
		  std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
		  mol.bonds[ibond].update(at_1, at_2, shorten.first, shorten.second,
					  empty, empty, root);
	       } else {
		  std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_second_atom =
		     mol.make_other_connections_to_second_atom_info(ibond);
		  std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
		  mol.bonds[ibond].update(at_1, at_2, shorten.first, shorten.second,
					  empty, other_connections_to_second_atom,
					  root);
	       }
	    } else {
	       std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
	       mol.bonds[ibond].update(at_1, at_2, shorten.first, shorten.second,
				       empty, empty, root);
	    }
	 } else {

	    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_first_atom =
	       mol.make_other_connections_to_first_atom_info(ibond);
	    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_second_atom =
	       mol.make_other_connections_to_second_atom_info(ibond);
	    mol.bonds[ibond].update(at_1, at_2, shorten.first, shorten.second,
				    other_connections_to_first_atom,
				    other_connections_to_second_atom, root);
	 }
      }
   }

   // redo the atoms, this time with widgets.
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {

      if (display_atom_names || display_atom_numbers) {

	 lig_build::atom_id_info_t atom_id_info = mol.atoms[iat].atom_name;
	 if (display_atom_numbers)
	    atom_id_info = mol.atoms[iat].element + ":" + coot::util::int_to_string(iat+1); // mol file numbering, 1-indexed
	 atom_id_info.size_hint = -1;
	 const std::string &ele = mol.atoms[iat].element;
	 std::string fc = font_colour(ele);
	 mol.atoms[iat].update_atom_id_forced(atom_id_info, fc, root);

      } else {
	 std::vector<unsigned int> local_bonds = mol.bonds_having_atom_with_atom_index(iat);
	 const std::string &ele = mol.atoms[iat].element;
	 bool gl_flag = false; // not a GL render engine
	 lig_build::atom_id_info_t atom_id_info =
	    mol.make_atom_id_by_using_bonds(iat, ele, local_bonds, gl_flag);

	 std::string fc = font_colour(ele);
	 if (ele != "C") {
	    mol.atoms[iat].update_atom_id_forced(atom_id_info, fc, root);
	 } else {
	    std::vector<unsigned int> bonds = mol.bonds_having_atom_with_atom_index(iat);
	    if (bonds.size() == 1) {
	       // CH3 superatom
	       mol.atoms[iat].update_atom_id_forced(atom_id_info, fc, root);
	    }
	 }
      }
   }

   draw_substitution_contour();

}

void
lbg_info_t::undo() {

   save_molecule_index--;
   if (save_molecule_index >= 0) {
      widgeted_molecule_t saved_mol = previous_molecules[save_molecule_index];
      // std::cout << "undo... reverting to save molecule number " << save_molecule_index << std::endl;
      import_from_widgeted_molecule(saved_mol);
      render();
      update_descriptor_attributes();
      update_apply_button_sensitivity();
   } else {
      clear(true);
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
   import_from_widgeted_molecule(copy_mol);
   render();

}



void
lbg_info_t::write_pdf(const std::string &file_name) const {

#if CAIRO_HAS_PDF_SURFACE

   cairo_surface_t *surface;
   cairo_t *cr;

   std::pair<lig_build::pos_t, lig_build::pos_t> extents = mol.ligand_extents();
   double pos_x = (extents.second.x + 220.0);
   double pos_y = (extents.second.y + 220.0);

   if (key_group) {
      pos_y += 240;
      // pos_x += 50;
      pos_x += 150;
   }
   surface = cairo_pdf_surface_create(file_name.c_str(), pos_x, pos_y);
   cr = cairo_create (surface);
   cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);

   /* Place it in the middle of our 9x10 page. */
   // cairo_translate (cr, 20, 130);
   cairo_translate(cr, 2, 13);

   goo_canvas_render(GOO_CANVAS(canvas), cr, NULL, 1.0);
   cairo_show_page(cr);
   cairo_surface_destroy(surface);
   cairo_destroy(cr);

#else
   std::cout << "No PDF (no PDF Surface in Cairo)" << std::endl;
#endif

}

void
lbg_info_t::write_ps(const std::string &file_name) const {

#if CAIRO_HAS_PS_SURFACE
   cairo_surface_t *surface;
   cairo_t *cr;

   std::pair<lig_build::pos_t, lig_build::pos_t> extents = mol.ligand_extents();
   double pos_x = (extents.second.x + 220.0);
   double pos_y = (extents.second.y + 220.0);

   if (key_group) {
      pos_y += 240;
      // pos_x += 50;
      pos_x += 150;
   }
   surface = cairo_ps_surface_create(file_name.c_str(), pos_x, pos_y);
   cr = cairo_create (surface);
   cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);

   /* Place it in the middle of our 9x10 page. */
   // cairo_translate (cr, 20, 130);
   cairo_translate(cr, 2, 13);

   goo_canvas_render(GOO_CANVAS(canvas), cr, NULL, 1.0);
   cairo_show_page(cr);
   cairo_surface_destroy(surface);
   cairo_destroy(cr);

#else
   std::cout << "No PS (no PS Surface in Cairo)" << std::endl;
#endif
}


void
lbg_info_t::write_svg(const std::string &file_name) const {

#if CAIRO_HAS_SVG_SURFACE

   std::pair<lig_build::pos_t, lig_build::pos_t> extents = mol.ligand_extents();
   double pos_x = (extents.second.x + 220.0);
   double pos_y = (extents.second.y + 220.0);

   if (key_group) {
      pos_y += 240;
      pos_x += 50;
   }

   cairo_surface_t *surface = cairo_svg_surface_create(file_name.c_str(), pos_x, pos_y);
   cairo_t *cr = cairo_create(surface);
   cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);

   /* Place it in the middle of our 9x10 page. */
   cairo_translate(cr, 2, 13);

   goo_canvas_render(GOO_CANVAS(canvas), cr, NULL, 1.0);
   cairo_show_page(cr);
   cairo_surface_destroy(surface);
   cairo_destroy(cr);
#endif

}

void
lbg_info_t::write_png(const std::string &file_name) {

   std::pair<lig_build::pos_t, lig_build::pos_t> extents = mol.ligand_extents();
   std::pair<lig_build::pos_t, lig_build::pos_t> extents_flev = flev_residues_extents();

   gdouble scale =  goo_canvas_get_scale(GOO_CANVAS(canvas));

   // std::cout << "goo_canvas_get_scale() " << scale << std::endl;
   // std::cout << "extents_flev: " << extents_flev.first << " " << extents_flev.second
   //           << std::endl;

   int size_x = int(extents.second.x) + 30;
   int size_y = int(extents.second.y) + 30;

   if (extents_flev.second.x > size_x) size_x = extents_flev.second.x;
   if (extents_flev.second.y > size_y) size_y = extents_flev.second.y;

   // the above extents give us the centre of residue circle - we want to see the whole circle,
   // so expand things a bit.

   size_x += 40;
   size_y += 40;

   if (key_group) {
      // std::cout << "adding space for (just) the key" << std::endl;
      size_y += 260; // *scale?
      size_x += 150;
   }
   // std::cout << "new size_x: " << size_x << std::endl;
   // std::cout << "new size_y: " << size_y << std::endl;

   cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, size_x, size_y);
   cairo_t *cr = cairo_create (surface);

   {

      int new_width  = size_x * scale;
      int new_height = size_y * scale;
      cairo_content_t content = CAIRO_CONTENT_COLOR_ALPHA;
      cairo_surface_t *new_surface = cairo_surface_create_similar(surface, content, new_width, new_height);
      cairo_t *cr = cairo_create (new_surface);
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);

      /* Scale *before* setting the source surface (1) */
      cairo_scale(cr, scale, scale);
      cairo_set_source_surface(cr, surface, 0, 0);
      cairo_pattern_set_extend(cairo_get_source(cr), CAIRO_EXTEND_REFLECT);
      /* Replace the destination with the source instead of overlaying */
      cairo_set_operator (cr, CAIRO_OPERATOR_SOURCE);
      /* Do the actual drawing */
      cairo_paint (cr);
      goo_canvas_render(GOO_CANVAS(canvas), cr, NULL, 1.0);
      cairo_surface_write_to_png(new_surface, file_name.c_str());
      cairo_surface_destroy(new_surface);
      cairo_destroy(cr);
   }
   cairo_surface_destroy(surface);
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
lbg_info_t::import_molecule_from_file(const std::string &file_name) { // mol or cif

   std::string ext = coot::util::file_name_extension(file_name);
   if (ext == ".cif") {
      import_molecule_from_cif_file(file_name);
   } else {
      if (ext == ".smi" || ext == ".smiles") {
	 import_mol_from_smiles_file(file_name);
      } else {
	 import_mol_from_file(file_name);
      }
   }
}

// Let's have a wrapper around sbase_import_func_ptr so that lbg-search doesn't need
// to ask if sbase_import_func_ptr is valid or not.
//
void
lbg_info_t::import_srs_monomer(const std::string &comp_id) {

   // this is passed from coot as it calls lbg()
   //
   if (! stand_alone_flag) {
      if (sbase_import_func_ptr) {
	 sbase_import_func_ptr(comp_id);
      } else {
	 std::cout << "ERROR:: null sbase_import_func_ptr" << std::endl;
      }
   } else {
      // here we are in lidia
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#ifdef HAVE_CCP4SRS
      coot::protein_geometry pg;

      const char *d1 = getenv("COOT_CCP4SRS_DIR"); // we should provide an "installed" dir
                                                   // if COOT_CCP4SRS_DIR is not set
      std::string srs_dir;
      if (d1)
	 srs_dir = d1;
      pg.init_ccp4srs(srs_dir);
      bool s = pg.fill_using_ccp4srs(comp_id);
      std::cout << "DEBUG:: pg.fill_using_ccp4srs() returned " << s << std::endl;
      if (s) {
	 int imol_local = 0; // dummy
	 std::pair<bool, coot::dictionary_residue_restraints_t> p =
	    pg.get_monomer_restraints(comp_id, imol_local);
	 if (p.first) {
	    bool show_hydrogens_flag = false;
	    import_via_rdkit_from_restraints_dictionary(p.second, show_hydrogens_flag);
	 } else {
	    std::cout << "ERROR:: bad extraction of restriants from srs " << comp_id << std::endl;
	 }
      }
#endif // HAVE_CCP4SRS
#endif // MAKE_ENHANCED_LIGAND_TOOLS
   }
}



void
lbg_info_t::import_molecule_from_cif_file(const std::string &file_name) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   coot::protein_geometry pg;
   pg.init_refmac_mon_lib(file_name, 43);
   std::vector<std::string> types = pg.monomer_types();
   if (types.size() > 0) {
      int imol = 0; // dummy
      std::string comp_id = types.back();
      std::pair<bool, coot::dictionary_residue_restraints_t> p =
	 pg.get_monomer_restraints(comp_id, imol);
      if (p.first) {
	 bool show_hydrogens_flag = false;
	 import_via_rdkit_from_restraints_dictionary(p.second, show_hydrogens_flag);
      }
   } else {
      std::cout << "import_molecule_from_cif_file() no types from "
		<< file_name << std::endl;
   }

#endif
}


// Mol, sdf or mol2, that is.
//
// What happens when this fails?  File name is missing?
// File is null?
// File is not a molecule format?
//
// 20140322: Let's throw a runtime_error.
//
void
lbg_info_t::import_mol_from_file(const std::string &file_name) {

   // if we have rdkit, try to read as an Mol2 file, if that fails try
   // to read as a Mol file.  If that fails try to read with my parser.
   //
   // if we don't have rdkit of course, just use my parser.

   bool try_as_mdl_mol = false;
#ifndef MAKE_ENHANCED_LIGAND_TOOLS
   // fallback
   try_as_mdl_mol = true;
#else
   try {
      RDKit::RWMol *m = RDKit::Mol2FileToMol(file_name);
      coot::set_3d_conformer_state(m);
      if (m) {
	 rdkit_mol_post_read_handling(m, file_name);
      } else {
	 // should throw an exception before getting here, I think.
	 std::cout << "Null m in import_mol_from_file() " << std::endl;
	 try_as_mdl_mol = true;
      }
   }
   catch (const RDKit::FileParseException &rte) {

      try {
	 bool sanitize = true;
	 bool removeHs = false;
	 bool strict_parsing = true;

	 // strict_parsing is not in the MolFileToMol() interface for old RDKits.
	 RDKit::RWMol *m = RDKit::MolFileToMol(file_name, sanitize, removeHs);
	 unsigned int iconf = 0;

	 if (!m) {
	    std::string s = "Null molecule from MolFile file \"";
	    s += file_name;
	    s += "\"";
	    throw std::runtime_error(s);
	 } else {
	    coot::set_3d_conformer_state(m);
	    if (coot::has_zero_coords(m, 0)) { // test for (,0,0,0) - not 2d test
	       iconf = RDDepict::compute2DCoords(*m, NULL, true);

	       // bond wedge fiddle here removed
	    }
	 }
	 rdkit_mol_post_read_handling(m, file_name, iconf);
      }
      catch (const RDKit::FileParseException &rte) {
	 try_as_mdl_mol = true;
      }
      catch (const RDKit::MolSanitizeException &rte) {
	 try_as_mdl_mol = true; // e.g. charges wrong.
      }
      catch (const Invar::Invariant &rte) {
	 try_as_mdl_mol = false; // e.g. bad input file
      }
   }
   catch (const RDKit::BadFileException &e) {
      std::cout << "WARNING:: Bad file " << file_name << " " << e.message() << std::endl;
      try_as_mdl_mol = true;
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING runtime_error in mol_to_asc_rdkit() " << rte.what() << std::endl;
      try_as_mdl_mol = true;
   }
   catch (const std::exception &e) {
      std::cout << "WARNING:: import_mol_from_file: exception: " << e.what() << std::endl;
      try_as_mdl_mol = true;
   }
#endif // MAKE_ENHANCED_LIGAND_TOOLS

   if (try_as_mdl_mol) {
      std::cout << "..................... using my mdl parser.... " << std::endl;
      // read as an MDL mol file
      //
      mmdb::Manager *mol = NULL; // no atom names to transfer
      lig_build::molfile_molecule_t mm;
      mm.read(file_name);
      widgeted_molecule_t wmol = import_mol_file(mm, file_name, mol);
      import_from_widgeted_molecule(wmol);
      render();
      update_descriptor_attributes();
      update_apply_button_sensitivity();
   }
}

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
// iconf is default arg value 0.
void
lbg_info_t::rdkit_mol_post_read_handling(RDKit::RWMol *m, const std::string &file_name, unsigned int iconf) {

   // 20160909-PE Paolo Tosco explains how atom names are stored in an RDKit molecule from sdf.
   if (false) {
      if (m) {
	 unsigned int n_atoms = m->getNumAtoms();
	 for (unsigned int iat=0; iat<n_atoms; iat++) {
	    RDKit::ATOM_SPTR at_p = (*m)[iat];
	    std::string name;
	    try {
	       at_p->getProp("molFileAlias", name);
	       std::cout << "name: " << name << std::endl;
	    }
	    catch (const KeyErrorException &kee) {
	    }
	 }
      }
   }

   try {
      // molfile molecules don't know about aromatic bonds, we need
      // to kekulize now.

      RDKit::MolOps::Kekulize(*m); // non-const reference?
      double weight_for_3d_distances = 0.4;

      int n_confs = m->getNumConformers();
      if (n_confs > 0) {
	 if (m->getConformer(iconf).is3D()) {
	    int iconf_local = coot::add_2d_conformer(m, weight_for_3d_distances); // 3d -> 2d
	    // (if not 3d, do nothing)
	    if (false)
	       std::cout << "rdkit_mol_post_read_handling() add_2d_conformer returned "
			 << iconf_local << std::endl;
	    if (iconf_local == -1)
	       std::cout << "WARNING:: import_mol_from_file() failed to make 2d conformer "
			 << std::endl;
	    iconf = iconf_local;
	 }

	 //    Old way (Pre-Feb 2013) goving via a molfile_molecule_t
	 //
	 //    lig_build::molfile_molecule_t mm = coot::make_molfile_molecule(*m, iconf);
	 //    mmdb::Manager *mol = NULL; // no atom names to transfer
	 //    widgeted_molecule_t wmol = import_mol_file(mm, file_name, mol);

	 const RDKit::Conformer conformer = m->getConformer(iconf);
	 RDKit::WedgeMolBonds(*m, &conformer);

	 widgeted_molecule_t wmol = import_rdkit_mol(m, iconf);
	 mdl_file_name = file_name;
	 import_from_widgeted_molecule(wmol);
	 render();
	 update_descriptor_attributes();
         update_apply_button_sensitivity();

      } else {
	 std::cout << "WARNING:: molecule from " << file_name << " had 0 conformers " << std::endl;
      }
   }
   catch (const std::exception &e) {
      // Vitamin B12 DB00115 throws a std::exception because it can't kekulize mol
      std::cout << "WARNING:: " << e.what() << " on reading " << file_name << std::endl;

      std::string status_string = "  Can't kekulize mol from " + file_name;
      guint statusbar_context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR(lbg_statusbar),
								status_string.c_str());
      gtk_statusbar_push(GTK_STATUSBAR(lbg_statusbar),
			 statusbar_context_id,
			 status_string.c_str());
   }
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

// pdb_mol is the pdb representation of the (flat) ligand - and it has
// the atom names.  We will add the atom names into mol by matching
// coordinates.
//
widgeted_molecule_t
lbg_info_t::import_mol_file(const lig_build::molfile_molecule_t &mol_in,
			    const std::string &file_name,
			    mmdb::Manager *pdb_mol) {

   widgeted_molecule_t new_mol(mol_in, pdb_mol);
   mdl_file_name = file_name;
   save_molecule();
   return new_mol;
}

void
lbg_info_t::import_mol_from_smiles_file(const std::string &file_name) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   if (coot::file_exists(file_name)) {

      std::string file_contents = file_to_string(file_name);
      std::string::size_type isplit = file_contents.find_first_of(' ');
      if (isplit != std::string::npos)
	 file_contents = file_contents.substr(0, isplit);
      isplit = file_contents.find_first_of('\n');
      if (isplit != std::string::npos)
	 file_contents = file_contents.substr(0, isplit);

      try {
	 RDKit::RWMol *m = RDKit::SmilesToMol(file_contents);
	 if (m) {
	    RDDepict::compute2DCoords(*m, NULL, true);
	    coot::set_3d_conformer_state(m); // checks for null
	    rdkit_mol_post_read_handling(m, file_name);
	 } else {
	    // should throw an exception before getting here, I think.
	    std::cout << "Null m in import_mol_from_file() " << std::endl;
	 }
      }
      catch (const RDKit::FileParseException &rte) {
	 std::cout << "WARNING:: something bad in import_mol_from_file() " << rte.what() << std::endl;
      }
   } else {
      std::cout << "WARNING:: No file: " << file_name << std::endl;
   }

#endif // MAKE_ENHANCED_LIGAND_TOOLS
}

void
lbg_info_t::import_mol_from_smiles_string(const std::string &smiles) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   try {
      RDKit::RWMol *m = RDKit::SmilesToMol(smiles);
      if (m) {
	 RDDepict::compute2DCoords(*m, NULL, true);
	 coot::set_3d_conformer_state(m); // checks for null
	 rdkit_mol_post_read_handling(m, "from-SMILES-string");
      } else {
	 // should throw an exception before getting here, I think.
	 std::cout << "Null m in import_mol_from_smiles_string() " << std::endl;
      }
   }
   catch (const RDKit::FileParseException &rte) {
      std::cout << "WARNING:: something bad in import_mol_from_smiles_string() " << rte.what() << std::endl;
   }


#else
   std::cout << "You need enhanced ligand tools version" << std::endl;

#endif // MAKE_ENHANCED_LIGAND_TOOLS
}

void
lbg_info_t::import_mol_from_comp_id(const std::string &comp_id,
				    bool show_hydrogens_status) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   int imol = 0; // dummy

   coot::protein_geometry pg;
   int dynamic_add_status = pg.try_dynamic_add(comp_id, 1);
   coot::dictionary_residue_restraints_t dict;
   std::pair<bool, coot::dictionary_residue_restraints_t> p =
      pg.get_monomer_restraints(comp_id, imol);
   bool have_dict = true;
   if (p.first) {
     have_dict = true;
     dict = p.second;
   } else {
#ifdef HAVE_CCP4SRS
     std::cout << "Trying to load " << comp_id << " from CCP4 SRS" << std::endl;
     pg.init_ccp4srs(".");
     bool status = pg.fill_using_ccp4srs(comp_id);
     if (status) {
	p = pg.get_monomer_restraints(comp_id, imol);
	if (p.first) {
	   dict = p.second;
	}
     } else {
	std::cout << "Failed to load " << comp_id << " using SRS " << std::endl;
     }
#endif // HAVE_CCP4SRS
   }

   if (have_dict) {
      import_via_rdkit_from_restraints_dictionary(dict, show_hydrogens_status);
   }
#else
   std::cout << "You need enhanced ligand tools version" << std::endl;
#endif  // MAKE_ENHANCED_LIGAND_TOOLS
}

// #include "GraphMol/MolPickler.h" VERSION problems.

void
lbg_info_t::import_via_rdkit_from_restraints_dictionary(const coot::dictionary_residue_restraints_t &dict, bool show_hydrogens_status) {

#ifdef  MAKE_ENHANCED_LIGAND_TOOLS

   // double try/catch because we want to separate rdkit molecule construction
   // problems from undelocalise/sanitize/compute2D problems.
   //
   try {
      RDKit::RWMol m = coot::rdkit_mol(dict);

      try {
	 coot::undelocalise(&m);
	 if (! show_hydrogens_status)
	    coot::remove_non_polar_Hs(&m);
	 unsigned int n_mol_atoms = m.getNumAtoms();
	 for (unsigned int iat=0; iat<n_mol_atoms; iat++)
	    m[iat]->calcImplicitValence(true);

	 // 20160702 use add_2d_conformer() instead now?

	 coot::rdkit_mol_sanitize(m);
	 RDKit::MolOps::Kekulize(m); // non-const reference?
	 bool canonOrient=true;
	 bool clearConfs=true;
	 unsigned int nFlipsPerSample=3;
	 unsigned int nSamples=200;
	 int sampleSeed=10;
	 bool permuteDeg4Nodes=true;

	 unsigned int conf_id = RDDepict::compute2DCoords(m, NULL,
							  canonOrient,
							  clearConfs,
							  nFlipsPerSample,
							  nSamples,
							  sampleSeed,
							  permuteDeg4Nodes);
	 RDKit::Conformer conf = m.getConformer(conf_id);

	 RDKit::WedgeMolBonds(m, &conf);

	 if (false)
	    std::cout << "..... n_confs B " << m.getNumConformers()
		      << " with new 2D conf_id " << conf_id
		      << " 3d-flag: " << m.getConformer(conf_id).is3D() << std::endl;

	 gtk_label_set_text(GTK_LABEL(lbg_toolbar_layout_info_label), dict.residue_info.comp_id.c_str());
	 rdkit_mol_post_read_handling(&m, "from-comp-id");
      }
      catch (const RDKit::MolSanitizeException &e) {
	 // calcImplicitValence() can make this happend
	 std::cout << "ERROR:: on Sanitize (inner) " << e.what() << std::endl;
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: (inner) import_via_rdkit_from_restraints_dictionary "
		   << rte.what() << std::endl;
	 coot::debug_rdkit_molecule(&m);

	 // boost::python::api::object obj = RDKit::MolToBinary(m);

// 	 {
// 	    std::string res;
// 	    RDKit::MolPickler::pickleMol(m,res);
// 	 }
      }
   }
   // we don't have an rdkit molecule for these catches.
   catch (const RDKit::MolSanitizeException &e) {
      // calcImplicitValence() can make this happend
      std::cout << "ERROR:: on Sanitize" << e.what() << std::endl;
   }
   catch (const std::runtime_error &rte) {
      std::cout << "ERROR:: " << rte.what() << std::endl;
   }

#endif // MAKE_ENHANCED_LIGAND_TOOLS
}




#ifdef MAKE_ENHANCED_LIGAND_TOOLS
widgeted_molecule_t
lbg_info_t::import_rdkit_mol(RDKit::ROMol *rdkm, int iconf) const {

   // see that this returns a widgeted_molecule_t, it doesn't fill the class's mol data item.

   // transfer atom names if you can.

   widgeted_molecule_t m;

   int n_conf  = rdkm->getNumConformers();
   if (iconf < n_conf) {
      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

      RDKit::Conformer &conf = rdkm->getConformer(iconf);
      unsigned int n_mol_atoms = rdkm->getNumAtoms();

      // determine the centre correction
      double sum_x = 0;
      double sum_y = 0;
      double min_y = 9e9;
      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	 RDKit::ATOM_SPTR at_p = (*rdkm)[iat];
	 RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
	 sum_x += r_pos.x;
	 sum_y += r_pos.y;
	 if (r_pos.y < min_y)
	    min_y = r_pos.y;
      }

      // set the scale correction
      //
      std::vector<double> bond_lengths;
      for (unsigned int i=0; i<rdkm->getNumBonds(); i++) {
	 const RDKit::Bond *bond_p = rdkm->getBondWithIdx(i);
	 int idx_1 = bond_p->getBeginAtomIdx();
	 int idx_2 = bond_p->getEndAtomIdx();
	 if ( (*rdkm)[idx_1]->getAtomicNum() != 1) {
	    if ( (*rdkm)[idx_2]->getAtomicNum() != 1) {
	       RDGeom::Point3D &r_pos_1 = conf.getAtomPos(idx_1);
	       RDGeom::Point3D &r_pos_2 = conf.getAtomPos(idx_2);
	       clipper::Coord_orth p1(r_pos_1.x, r_pos_1.y, r_pos_1.z);
	       clipper::Coord_orth p2(r_pos_2.x, r_pos_2.y, r_pos_2.z);
	       double l = clipper::Coord_orth::length(p1, p2);
	       bond_lengths.push_back(l);
	    }
	 }
      }
      if (bond_lengths.size() > 0) {
	 std::sort(bond_lengths.begin(), bond_lengths.end());
	 int index = bond_lengths.size()/2;
	 double bll = bond_lengths[index];
	 double scale = 1.0/bll;
	 m.scale_correction.first = 1;
	 m.scale_correction.second = scale;
      }


      if (n_mol_atoms > 0) {
	 double centre_x = sum_x/double(n_mol_atoms);
	 double centre_y = sum_y/double(n_mol_atoms);
	 m.centre_correction = lig_build::pos_t(centre_x, centre_y);
	 m.mol_in_min_y = min_y;
      }

      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	 RDKit::ATOM_SPTR at_p = (*rdkm)[iat];
	 RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
	 std::string name = "";
	 try {
	    at_p->getProp("name", name);
	 }
	 catch (const KeyErrorException &kee) {
	    // we don't need to see these.  We get them when reading an mdl file
	    // (for example).
	    // std::cout << "caught no-name for atom exception in import_rdkit_mol(): "
	    // <<  kee.what() << std::endl;
	 }
	 try {
	    at_p->getProp("molFileAlias", name);
	 }
	 catch (const KeyErrorException &kee) { }
	 clipper::Coord_orth cp(r_pos.x , r_pos.y, r_pos.z);
	 lig_build::pos_t pos = m.input_coords_to_canvas_coords(cp);
	 int n = at_p->getAtomicNum();
	 std::string element = tbl->getElementSymbol(n);
	 int charge = at_p->getFormalCharge();
	 widgeted_atom_t mol_at(pos, element, charge, NULL);

	 mol_at.charge = charge;
	 if (! name.empty())
	    mol_at.atom_name = name;
	 std::pair<bool,int> added = m.add_atom(mol_at);
      }

      unsigned int n_bonds = rdkm->getNumBonds();
      for (unsigned int ib=0; ib<n_bonds; ib++) {
	 const RDKit::Bond *bond_p = rdkm->getBondWithIdx(ib);
	 int idx_1 = bond_p->getBeginAtomIdx();
	 int idx_2 = bond_p->getEndAtomIdx();
	 lig_build::bond_t::bond_type_t bt = coot::convert_bond_type(bond_p->getBondType());

	 try {
	    const widgeted_atom_t &wat1 = m[idx_1];
	    const widgeted_atom_t &wat2 = m[idx_2];
	    bool shorten_first  = false;
	    bool shorten_second = false;
	    if ( (*rdkm)[idx_1]->getAtomicNum() != 6) {
	       shorten_first = true;
	    } else {
	       // how many bonds does this atom have?
	       // if (mol[idx_1].getNumBonds() == 1)
	       //  shorten_first = true;
	    }
	    if ( (*rdkm)[idx_2]->getAtomicNum() != 6) {
	       shorten_second = true;
	    } else {
	    }

	    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;

	    widgeted_bond_t bond(idx_1, idx_2, wat1, wat2, shorten_first, shorten_second, bt, empty, empty, NULL);
	    RDKit::Bond::BondDir bond_dir = bond_p->getBondDir();
	    if (false)
	       std::cout << "bond " << ib << ":  type " << bt
			 << " between " << idx_1 << " at "
			 << conf.getAtomPos(idx_1)
			 << " and " << idx_2 << " at "
			 << conf.getAtomPos(idx_2)
			 << " dir " << bond_dir << std::endl;
	    if (bond_dir != RDKit::Bond::NONE) {
	       if (bond_dir == RDKit::Bond::BEGINWEDGE) {
		  bond.set_bond_type(lig_build::bond_t::OUT_BOND);
	       }
	       if (bond_dir == RDKit::Bond::BEGINDASH)
		  bond.set_bond_type(lig_build::bond_t::IN_BOND);
	    }
	    m.add_bond(bond);
	 }
	 catch (...) {
	    std::cout << "WARNING:: problem. scrambled input molecule? numbers of atoms: ";
	    std::cout << rdkm->getNumAtoms() << " vs "
		      << m.get_number_of_atoms_including_hydrogens() << std::endl;
	 }
      }

      m.assign_ring_centres();
   }
   return m;
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

void
lbg_info_t::clean_up_2d_representation() {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   if (use_graphics_interface_flag) {
      try {
	 RDKit::RWMol rdkm = rdkit_mol(mol);
	 coot::rdkit_mol_sanitize(rdkm);
	 RDKit::MolOps::Kekulize(rdkm); // non-const reference?
	 bool canonOrient=true;
	 bool clearConfs=true;
	 unsigned int nFlipsPerSample=3;
	 unsigned int nSamples=200;
	 int sampleSeed=10;
	 bool permuteDeg4Nodes=true;
	 RDDepict::compute2DCoords(rdkm, NULL,
				   canonOrient,
				   clearConfs,
				   nFlipsPerSample,
				   nSamples,
				   sampleSeed,
				   permuteDeg4Nodes);
	 int iconf = 0;
	 RDKit::Conformer conf = rdkm.getConformer(iconf);
	 RDKit::WedgeMolBonds(rdkm, &conf);
	 lig_build::molfile_molecule_t mm =
	    coot::make_molfile_molecule(rdkm, iconf);
	 mmdb::Manager *mol = NULL; // no atom names to transfer
	 widgeted_molecule_t wmol(mm, mol);
	 import_from_widgeted_molecule(wmol);
	 render();
	 update_descriptor_attributes();
	 save_molecule();
      }
      catch (const std::exception &e) {
	 std::cout << "WARNING:: clean_up_2d_representation() " << e.what() << std::endl;
      }
   }
#endif // MAKE_ENHANCED_LIGAND_TOOLS
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
   catch (const std::exception &rte) {
      std::cout << rte.what() << std::endl;
   }
}

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
std::string
lbg_info_t::get_smiles_string(const RDKit::ROMol &mol) const {

   std::string s;
   if (mol.getNumAtoms() > 0) {
      // new
      RDKit::ROMol mol_copy(mol);
      s = RDKit::MolToSmiles(mol_copy);
   }
   return s;
}
#endif

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
void
lbg_info_t::update_statusbar_smiles_string(const RDKit::ROMol &mol) const {
   std::string s = get_smiles_string(mol);
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
   if (txt)
      get_drug(std::string(txt));
   else
      std::cout << "ERROR:: null text in get_drug_using_entry_text() " << std::endl;
}

// get mol file and load it
void
lbg_info_t::get_drug(const std::string &drug_name) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   bool status = false;
   std::string status_string;

   if (get_drug_mdl_file_function_pointer) {

      try {
         if (get_drug_mdl_file_function_pointer) {
            std::cout << "DEBUG:: Using get_drug_mdl_file_function_pointer" << std::endl;
            // this could fail for SSL reasons. Try to dig out the libcurl error
            std::string file_name = get_drug_mdl_file_function_pointer(drug_name);
            if (file_name.empty()) {
               std::cout << "WARNING:: in get_drug(): empty mol file name." << std::endl;
            } else {
               status = true;
               import_mol_from_file(file_name);
               save_molecule();
            }
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "WARNING:: " << rte.what() << std::endl;
      }

   }

   if (!status) {

      // If we have failed so far, try again with lidia.fetch.
      // But import rdkit first so that import.fetch doesn't crash - !
      std::cout << "DEBUG:: --- start import rdkit/lidia.fetch --- " << std::endl;
      PyObject *pName_rdkit = PyString_FromString("rdkit");
      std::cout << "DEBUG:: --- we have pName_rdkit --- " << pName_rdkit << std::endl;
      PyObject *pModule = PyImport_Import(pName_rdkit);
      if (! pModule) return;
      std::cout << "DEBUG:: --- we have imported rdkit --- " << pName_rdkit << std::endl;

      PyObject *pName = PyString_FromString("lidia.fetch");
      std::cout << "DEBUG:: with lidia.fetch pName " << pName << std::endl;
      pModule = PyImport_Import(pName);
      if (pModule == NULL) {
         std::cout << "NULL pModule" << std::endl;
      } else {
         // std::cout << "Hooray - found pModule" << std::endl;

         PyObject *pDict = PyModule_GetDict(pModule);
         if (! PyDict_Check(pDict)) {
            std::cout << "DEBUG:: pDict is not a dict" << std::endl;
         } else {

            PyObject *pFunc = PyDict_GetItemString(pDict, "fetch_molecule");
            if (PyCallable_Check(pFunc)) {
               std::cout << "DEBUG:: fetch_molecule is a callable function" << std::endl;
               PyObject *arg_list = PyTuple_New(1);
               PyObject *drug_name_py = PyString_FromString(drug_name.c_str());
               PyTuple_SetItem(arg_list, 0, drug_name_py);
               std::cout << "DEBUG:: fetch_molecule called with arg " << drug_name << std::endl;
               PyObject *result_py = PyEval_CallObject(pFunc, arg_list);
               std::cout << "DEBUG:: fetch_molecule got result " << result_py << std::endl;

               if (result_py) {

                  if (PyString_Check(result_py)) {
                     std::cout << "fetch_molecule result was a string " << std::endl;
                     std::string file_name = PyString_AsString(result_py);
                     std::cout << "fetch_molecule result was a file_name " << file_name << std::endl;

                     try {
                        bool sanitize = true;
                        bool removeHs = false;
                        RDKit::RWMol *m = RDKit::MolFileToMol(file_name, sanitize, removeHs);
                        unsigned int iconf = 0;

                        // we only want to compute coords if there are no coords
                        //
                        if (coot::has_zero_coords(m, 0))
                        iconf = RDDepict::compute2DCoords(*m, NULL, true);
                        if (m->getNumConformers() == 0)
                        iconf = RDDepict::compute2DCoords(*m, NULL, true);

                        if (m) {
                           rdkit_mol_post_read_handling(m, file_name, iconf);
                           status = true;
                        }
                     }

                     catch (const RDKit::FileParseException &fpe) {
                        std::cout << "WARNING::" << fpe.what() << std::endl;
                        status_string = "File parsing error on reading " + file_name;
                     }
                  } else {
                     std::cout << "Null result" << std::endl;
                     status_string = "Null result on fetching " + drug_name;
                  }
               } else {
                  std::cout << "results_py was not a string" << std::endl;
                  status_string = "Failed to download molecule " + drug_name;
               }

            } else {
               std::cout << "Non-callable pFunc" << std::endl;
               status_string = "Non-callable pFunc";
            }
         }
      }
   }

   if (! status) {
      // problem
      clear(true);
      if (status_string.empty())
      status_string = "  Problem downloading molecule " + drug_name;
      guint statusbar_context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR(lbg_statusbar),
      status_string.c_str());
      gtk_statusbar_push(GTK_STATUSBAR(lbg_statusbar),
      statusbar_context_id,
      status_string.c_str());

   }
#else
   std::cout << "WARNING:: Not compiled with enhanced-ligand-tools - Nothing doing Hermione"
	     << std::endl;
#endif // MAKE_ENHANCED_LIGAND_TOOLS
}


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
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
	       if (true)
		  std::cout << "DEBUG:: found " << module_name << " " << function_name << " function"
			    << std::endl;
	    }
	 }
      }
   }
   return extracted_func;
}
#endif
#endif

// std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
// lbg_info_t::make_other_connections_to_second_atom(unsigned int bond_index) const {

//    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > v;
//    return v;
// }



// flipping
void
lbg_info_t::flip_molecule(int axis) {

   widgeted_molecule_t new_mol = mol;
   new_mol.flip(axis);
   // render_from_molecule(new_mol); // wipes mol
   // update_descriptor_attributes();

   import_from_widgeted_molecule(new_mol);
   render();
}

void
lbg_info_t::rotate_z_molecule(double degrees) {

   widgeted_molecule_t new_mol = mol;
   new_mol.rotate_z(degrees);
   // update_descriptor_attributes();
   import_from_widgeted_molecule(new_mol);
   render();
}

// in degrees (used in on_lbg_view_rotate_apply_button_clicked
// callback).
void
lbg_info_t::rotate_z_molecule(const std::string &angle_str) {

   try {
      double angle = coot::util::string_to_double(angle_str);
      rotate_z_molecule(angle);
   }
   catch (const std::exception &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
}



void
lbg_info_t::pe_test_function() {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   std::cout << "PE test function" << std::endl;

   std::cout << "identify bond..." << std::endl;
   bond_pick_pending = true;


#ifdef USE_PYTHON
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

#endif
}


#endif // HAVE_GOOCANVAS
