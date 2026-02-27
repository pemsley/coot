/* src/c-interface-build.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by Bernhard Lohkamp
 * Copyright 2008 by Kevin Cowtan
 * Copyright 2007, 2008, 2009, 2010, 2011 The University of Oxford
 * Copyright 2016 by Medical Research Council
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#define HAVE_CIF  // will become unnessary at some stage.

#include <sys/types.h> // for stating
#include <sys/stat.h>
#include <string.h> // strncmp
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#include <windows.h>
#endif


#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

#include "globjects.h" //includes gtk/gtk.h



#include "graphics-info.h"

#include "coot-utils/coot-coord-utils.hh"
#include "utils/coot-fasta.hh"

#include "skeleton/BuildCas.h"
#include "ligand/helix-placement.hh"
#include "ligand/fast-ss-search.hh"

#include "utils/coot-utils.hh"  // for is_member_p
#include "coot-utils/coot-map-heavy.hh"  // for fffear

#include "guile-fixups.h"


#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"

#include "ligand/ligand.hh" // for rigid body fit by atom selection.

#include "cmtz-interface.hh" // for valid columns mtz_column_types_info_t
#include "c-interface-mmdb.hh"
#include "c-interface-scm.hh"
#include "c-interface-python.hh"
#include "new-molecule-by-symmetry-matrix.hh"


#ifdef USE_DUNBRACK_ROTAMERS
#include "ligand/dunbrack.hh"
#else 
#include "ligand/richardson-rotamer.hh"
#endif

#include "ligand/backrub-rotamer.hh"
#include "rotamer-search-modes.hh"

#include "protein_db/protein_db_utils.h"
#include "protein_db-interface.hh"

#include "cootilus/cootilus-build.h"

#include "c-interface-refine.hh"

#include "utils/logging.hh"
extern logging logger;




// --------------------------------------------------------------
//                 symmetry
// --------------------------------------------------------------

/* for shelx FA pdb files, there is no space group.  So allow the user
   to set it.  This can be initted with a HM symbol or a symm list for
   clipper */
short int set_space_group(int imol, const char *spg) {

   short int r = 0;

   // we should test that this is a clipper string here...
   // does it contain X Y and Z chars?
   // 
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].set_mmdb_symm(spg);
   }

   return r; 
}

//! \brief set the unit cell for a given molecule
//!
//! @return  the success status of the setting (1 good, 0 fail).
int set_unit_cell_and_space_group(int imol, float a, float b, float c,
				  float alpha, float beta, float gamma,
				  const char *sp_in) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::pair<std::vector<float>, std::string> cs_pair;
      cs_pair.second = sp_in;
      cs_pair.first.resize(6);
      cs_pair.first[0] = a;     cs_pair.first[1] = b;    cs_pair.first[2] = c;
      cs_pair.first[3] = alpha; cs_pair.first[4] = beta; cs_pair.first[5] = gamma;
      g.molecules[imol].set_mmdb_cell_and_symm(cs_pair); // no return value
      status = 1; // hmm.
   }
   return status;
}

//! \brief set the unit cell for a given molecule using the cell of moecule imol_from
//!
//! @return  the success status of the setting (1 good, 0 fail).
int set_unit_cell_and_space_group_using_molecule(int imol, int imol_from) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_model_molecule(imol_from)) {
	 graphics_info_t g;
	 std::pair<std::vector<float>, std::string> cs_pair = g.molecules[imol_from].get_cell_and_symm();
	 g.molecules[imol].set_mmdb_cell_and_symm(cs_pair); // void
	 status =1 ;
      }
   }
   return status;
}



// 20250612-PE old-style
void setup_save_symmetry_coords() {

   graphics_info_t::in_save_symmetry_define = 1;
   std::string s = "Now click on a symmetry atom";
   graphics_info_t g;
   g.add_status_bar_text(s);
   pick_cursor_maybe();

}

void
on_save_symm_coords_filechooser_dialog_response(GtkDialog *dialog, int response) {

   if (response == GTK_RESPONSE_YES) {
      // zooop off to c-interface-gui.cc, which  then calls save_symmetry_coords() below
      save_symmetry_coords_from_filechooser(GTK_WIDGET(dialog)); // not ideal casting
   }
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

#include "coot-fileselections.h"

void save_symmetry_coords_based_on_position() {

   graphics_info_t g;
   coot::Symm_Atom_Pick_Info_t sap = g.symmetry_atom_close_to_screen_centre();
   if (sap.success == GL_TRUE) {
      if (is_valid_model_molecule(sap.imol)) {
         // 20250612-PE We don't use the save_symmetry_coords_filechooser_dialog
         // because it doesn't have buttons (I don't known why).
         // so, do instead like on_save_coords_save_button_clicked() (overlay button callback)

         GtkWindow *parent_window = GTK_WINDOW(graphics_info_t::get_main_window());
         GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
         GtkWidget *file_chooser_dialog =
            gtk_file_chooser_dialog_new("Save Symm Coordinates",
                                        parent_window,
                                        action,
                                        ("_Cancel"),
                                        GTK_RESPONSE_CANCEL,
                                        ("_Save"),
                                        GTK_RESPONSE_YES,
                                        NULL);

            // used like this coot::Symm_Atom_Pick_Info_t *symm_info =
            // (coot::Symm_Atom_Pick_Info_t *) g_object_get_data(G_OBJECT(filechooser), "symm_info");

            coot::Symm_Atom_Pick_Info_t *sap_p = new coot::Symm_Atom_Pick_Info_t(sap);
            g_object_set_data(G_OBJECT(file_chooser_dialog), "symm_info", sap_p);
            g_signal_connect(file_chooser_dialog, "response",
                             G_CALLBACK(on_save_symm_coords_filechooser_dialog_response), NULL);
            gtk_widget_set_visible(file_chooser_dialog, TRUE);
            set_file_for_save_filechooser(file_chooser_dialog);
            add_filename_filter_button(file_chooser_dialog, COOT_SAVE_COORDS_FILE_SELECTION);

      }
   } else {
      std::string mess = "Not centred on symmetry atom";
      add_status_bar_text(mess);
      logger.log(log_t::WARNING, logging::function_name_t("save_symmetry_coords_based_on_position"), mess);
   }
}


void save_symmetry_coords(int imol,
			  const char *filename,
			  int symop_no,
			  int shift_a,
			  int shift_b,
			  int shift_c,
			  int pre_shift_to_origin_na,
			  int pre_shift_to_origin_nb,
			  int pre_shift_to_origin_nc) {

   // Copy the coordinates molecule manager
   // Transform them
   // write them out

   if (imol >= 0) {
      if (imol < graphics_info_t::n_molecules()) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    mmdb::Manager *mol2 = new mmdb::Manager;
	    mol2->Copy(graphics_info_t::molecules[imol].atom_sel.mol, mmdb::MMDBFCM_All);

	    atom_selection_container_t asc = make_asc(mol2);
	    mmdb::mat44 mat;
	    mmdb::mat44 mat_origin_shift;

	    mol2->GetTMatrix(mat_origin_shift, 0,
			     -pre_shift_to_origin_na,
			     -pre_shift_to_origin_nb,
			     -pre_shift_to_origin_nc);
			     
	    mol2->GetTMatrix(mat, symop_no, shift_a, shift_b, shift_c);

	    clipper::RTop_orth to_origin_rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
					      clipper::Coord_orth(mat_origin_shift[0][3],
								  mat_origin_shift[1][3],
								  mat_origin_shift[2][3]));
	    
	    for (int i=0; i<asc.n_selected_atoms; i++) {
	       clipper::Coord_orth co;
	       clipper::Coord_orth trans_pos; 

	       clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
						  mat[1][0], mat[1][1], mat[1][2],
						  mat[2][0], mat[2][1], mat[2][2]);
	       clipper::Coord_orth  cco(mat[0][3], mat[1][3], mat[2][3]);
	       clipper::RTop_orth rtop(clipper_mat, cco);
	       co = clipper::Coord_orth(asc.atom_selection[i]->x, 
					asc.atom_selection[i]->y, 
					asc.atom_selection[i]->z);
	       trans_pos = co.transform(to_origin_rtop);
	       trans_pos = trans_pos.transform(rtop);
	       asc.atom_selection[i]->x = trans_pos.x();
	       asc.atom_selection[i]->y = trans_pos.y();
	       asc.atom_selection[i]->z = trans_pos.z();
	    } 
	    asc.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
	    asc.mol->FinishStructEdit();
	    
	    mmdb_manager_delete_conect(mol2);
	    int ierr = -1;
	    if (coot::is_mmcif_filename(filename))
	       ierr = mol2->WriteCIFASCII(filename);
	    else
	       ierr = mol2->WritePDBASCII(filename);
	    if (ierr) {
	       std::cout << "WARNING:: WritePDBASCII to " << filename << " failed." << std::endl;
	       std::string s = "WARNING:: WritePDBASCII to file ";
	       s += filename;
	       s += " failed.";
	       graphics_info_t g;
	       g.add_status_bar_text(s);
	    } else {
	       // std::cout << "INFO:: Wrote symmetry atoms to " << filename << "." << std::endl;
	       logger.log(log_t::INFO, logging::function_name_t("save_symmetry_coords"), "Wrote symmetry atoms to", filename);
	       std::string s = "INFO:: Wrote symmetry atoms to file ";
	       s += filename;
	       s += ".";
	       graphics_info_t g;
	       g.add_status_bar_text(s);
	    }
	    
	    std::vector<std::string> command_strings;
	    command_strings.push_back("save-symmetry-coords");
	    command_strings.push_back(coot::util::int_to_string(imol));
	    command_strings.push_back(single_quote(filename));
	    command_strings.push_back(coot::util::int_to_string(symop_no));
	    command_strings.push_back(coot::util::int_to_string(shift_a));
	    command_strings.push_back(coot::util::int_to_string(shift_b));
	    command_strings.push_back(coot::util::int_to_string(shift_c));
	    command_strings.push_back(coot::util::int_to_string(pre_shift_to_origin_na));
	    command_strings.push_back(coot::util::int_to_string(pre_shift_to_origin_nb));
	    command_strings.push_back(coot::util::int_to_string(pre_shift_to_origin_nc));
	    add_to_history(command_strings);
	 }
      }
   }
}

/*! \brief create a new molecule (molecule number is the return value)
  from imol. 

The rotation/translation matrix components are given in *orthogonal*
coordinates.

Allow a shift of the coordinates to the origin before symmetry
expansion is apllied.

Return -1 on failure. */ 
int new_molecule_by_symmetry(int imol,
			     const char *name_in,
			     double m11, double m12, double m13, 
			     double m21, double m22, double m23, 
			     double m31, double m32, double m33, 
			     double tx, double ty, double tz,
			     int pre_shift_to_origin_na,
			     int pre_shift_to_origin_nb,
			     int pre_shift_to_origin_nc) {
   int istate = -1;
   if (is_valid_model_molecule(imol)) { 
      std::pair<bool, clipper::Cell> cell_info = graphics_info_t::molecules[imol].cell();

      mmdb::Manager *mol_orig = graphics_info_t::molecules[imol].atom_sel.mol;
      // test if returned molecule is non-null
      std::string name = "Symmetry copy of ";
      name += coot::util::int_to_string(imol);
      if (std::string(name_in) != "")
	 name = name_in;
      mmdb::Manager *mol_symm = new_molecule_by_symmetry_matrix_from_molecule(mol_orig,
									     m11, m12, m13,
									     m21, m22, m23,
									     m31, m32, m33,
									     tx, ty, tz,
									     pre_shift_to_origin_na,
									     pre_shift_to_origin_nb,
									     pre_shift_to_origin_nc);

      if (mol_symm) { 
	 int imol_new = graphics_info_t::create_molecule(); 
	 atom_selection_container_t asc = make_asc(mol_symm);
	 graphics_info_t g;
	 g.molecules[imol_new].install_model(imol_new, asc, g.Geom_p(), name, 1);
	 g.molecules[imol].set_have_unsaved_changes_from_outside();
	 update_go_to_atom_window_on_new_mol();
	 graphics_draw();
	 istate = imol_new;
      } else { 
	 std::cout << "WARNING:: molecule " << imol << " does not have a proper cell " 
		   << std::endl;
      }
   } else { 
	 std::cout << "WARNING:: molecule " << imol << " is not a valid model molecule " 
		   << std::endl;
   }
   return istate; 
} 



int new_molecule_by_symop(int imol, const char *symop_string, 
			  int pre_shift_to_origin_na,
			  int pre_shift_to_origin_nb,
			  int pre_shift_to_origin_nc) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) { 
      std::pair<bool, clipper::Cell> cell_info = graphics_info_t::molecules[imol].cell();
      if (cell_info.first) { 
	 coot::symm_card_composition_t sc(symop_string);
	 std::cout << symop_string << " ->\n" 
		   << sc.x_element[0] << " " << sc.y_element[0] << " " << sc.z_element[0] << "\n"
		   << sc.x_element[1] << " " << sc.y_element[1] << " " << sc.z_element[1] << "\n"
		   << sc.x_element[2] << " " << sc.y_element[2] << " " << sc.z_element[2] << "\n"
		   << "translations: "
		   << sc.trans_frac(0) << " "
		   << sc.trans_frac(1) << " "
		   << sc.trans_frac(2) << std::endl;
	 std::cout << "pre-trans: "
		   << pre_shift_to_origin_na << " "
		   << pre_shift_to_origin_nb << " " 
		   << pre_shift_to_origin_nc << std::endl;

	 // those matrix elements are in fractional coordinates.  We want to
	 // pass a components for an rtop_orth.

	 clipper::Mat33<double> mat(sc.x_element[0], sc.y_element[0], sc.z_element[0], 
				    sc.x_element[1], sc.y_element[1], sc.z_element[1], 
				    sc.x_element[2], sc.y_element[2], sc.z_element[2]);

	 clipper::Vec3<double> vec(sc.trans_frac(0),
				   sc.trans_frac(1),
				   sc.trans_frac(2));

	 clipper::RTop_frac rtop_frac(mat, vec);
	 clipper::RTop_orth rtop_orth = rtop_frac.rtop_orth(cell_info.second);
	 clipper::Mat33<double> orth_mat = rtop_orth.rot();
	 clipper::Coord_orth    orth_trn(rtop_orth.trn());

	 std::string new_mol_name = "SymOp_";
	 new_mol_name += symop_string;
	 new_mol_name += "_Copy_of_";
	 new_mol_name += coot::util::int_to_string(imol);

	 imol_new =  new_molecule_by_symmetry(imol,
					      new_mol_name.c_str(),
					      orth_mat(0,0), orth_mat(0,1), orth_mat(0,2), 
					      orth_mat(1,0), orth_mat(1,1), orth_mat(1,2), 
					      orth_mat(2,0), orth_mat(2,1), orth_mat(2,2), 
					      orth_trn.x(), orth_trn.y(), orth_trn.z(),
					      pre_shift_to_origin_na,
					      pre_shift_to_origin_nb,
					      pre_shift_to_origin_nc);
      }
   }
   return imol_new;
}


int
new_molecule_by_symmetry_with_atom_selection(int imol, 
					     const char *name,
					     const char *mmdb_atom_selection_string,
					     double m11, double m12, double m13, 
					     double m21, double m22, double m23, 
					     double m31, double m32, double m33, 
					     double tx, double ty, double tz,
					     int pre_shift_to_origin_na,
					     int pre_shift_to_origin_nb,
					     int pre_shift_to_origin_nc) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol_orig = graphics_info_t::molecules[imol].atom_sel.mol;
      int SelectionHandle = mol_orig->NewSelection();
      mol_orig->Select(SelectionHandle, mmdb::STYPE_ATOM,
		       mmdb_atom_selection_string,
		       mmdb::SKEY_OR);

      mmdb::Manager *mol =
	 coot::util::create_mmdbmanager_from_atom_selection(mol_orig, SelectionHandle);

      mmdb::Manager *new_mol = new_molecule_by_symmetry_matrix_from_molecule(mol,
									    m11, m12, m13,
									    m21, m22, m23,
									    m31, m32, m33,
									    tx, ty, tz,
									    pre_shift_to_origin_na,
									    pre_shift_to_origin_nb,
									    pre_shift_to_origin_nc);

      delete mol; // done with it
      if (new_mol) {
	 imol_new = graphics_info_t::create_molecule();
	 atom_selection_container_t asc = make_asc(new_mol);
	 graphics_info_t g;
	 g.molecules[imol_new].install_model(imol_new, asc, g.Geom_p(), name, 1);
	 g.molecules[imol].set_have_unsaved_changes_from_outside();
	 update_go_to_atom_window_on_new_mol();
	 graphics_draw();
      }
      mol_orig->DeleteSelection(SelectionHandle);
   }
   return imol_new;
}

mmdb::Manager *new_molecule_by_symmetry_matrix_from_molecule(mmdb::Manager *mol,
							    double m11, double m12, double m13, 
							    double m21, double m22, double m23, 
							    double m31, double m32, double m33, 
							    double tx, double ty, double tz,
							    int pre_shift_to_origin_na,
							    int pre_shift_to_origin_nb,
							    int pre_shift_to_origin_nc) {
   mmdb::Manager *new_mol = 0;

   try {
      std::pair<clipper::Cell, clipper::Spacegroup> cell_info = coot::util::get_cell_symm(mol);
      std::vector<int> pre_shift(3);
      pre_shift[0] = pre_shift_to_origin_na;
      pre_shift[1] = pre_shift_to_origin_nb;
      pre_shift[2] = pre_shift_to_origin_nc;
      clipper::Mat33<double> mat(m11, m12, m13, m21, m22, m23, m31, m32, m33);
      clipper::Vec3<double> vec(tx, ty, tz);
      clipper::RTop_orth rtop_orth(mat, vec);
      clipper::RTop_frac rtop_frac = rtop_orth.rtop_frac(cell_info.first);
      if (0) { 
	 std::cout << "DEBUG:: tx,ty,tz:  " << tx << " " << ty << " " << tz << std::endl;
	 std::cout << "DEBUG:: mol_by_symmetry() passed args:\n" << cell_info.first.format()
		   << std::endl << rtop_frac.format() << "     "
		   << pre_shift_to_origin_na << " "
		   << pre_shift_to_origin_nb << " "
		   << pre_shift_to_origin_nc << " "
		   << std::endl;
      }
      new_mol = coot::mol_by_symmetry(mol, cell_info.first, rtop_frac, pre_shift);
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
   } 
   return new_mol;
}

/*! \brief return the number of symmetry operators for the given molecule

return -1 on no-symmetry for molecule or inappropriate imol molecule number */
int n_symops(int imol) {

   int r = -1;

   if (is_valid_model_molecule(imol)) {
      std::pair<std::vector<float>, std::string> cs =
	 graphics_info_t::molecules[imol].get_cell_and_symm();
      if (cs.second.length() >0) { 
	 r = graphics_info_t::molecules[imol].atom_sel.mol->GetNumberOfSymOps();
      }
   };

   if (is_valid_map_molecule(imol)) {
      r = graphics_info_t::molecules[imol].xmap.spacegroup().num_symops();
   } 
   return r;
}

/* This function works by active symm atom. */
int move_reference_chain_to_symm_chain_position() {

   graphics_info_t g;
   return g.move_reference_chain_to_symm_chain_position();
} 



#ifdef USE_GUILE
/*! \brief return the pre-shift as a list of fraction or scheme false
  on failure  */
SCM origin_pre_shift_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      try { 
	 clipper::Coord_frac cf = coot::util::shift_to_origin(mol);
	 r = SCM_EOL;
	 r = scm_cons(scm_from_int(int(round(cf.w()))), r);
	 r = scm_cons(scm_from_int(int(round(cf.v()))), r);
	 r = scm_cons(scm_from_int(int(round(cf.u()))), r);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << rte.what() << std::endl;
      } 
   } 
   return r;
} 
#endif  /* USE_GUILE */

#ifdef USE_PYTHON
/*! \brief return the pre-shift as a list of fraction or python false
  on failure  */
PyObject *origin_pre_shift_py(int imol) {
   
   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      try { 
	 clipper::Coord_frac cf = coot::util::shift_to_origin(mol);
	 r = PyList_New(0);
	 PyList_Append(r, PyLong_FromLong(int(round(cf.u()))));
	 PyList_Append(r, PyLong_FromLong(int(round(cf.v()))));
	 PyList_Append(r, PyLong_FromLong(int(round(cf.w()))));
      }
      catch (const std::runtime_error &rte) {
	 std::cout << rte.what() << std::endl;
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif  /* USE_PYTHON */


/*  ----------------------------------------------------------------------- */
/*                  cis <-> trans conversion                                */
/*  ----------------------------------------------------------------------- */
void do_cis_trans_conversion_setup(int istate) {

   if (istate == 1) {

      graphics_info_t::in_cis_trans_convert_define = 1;
      pick_cursor_maybe(); // depends on ctrl key for rotate

   } else {
      graphics_info_t::in_cis_trans_convert_define = 0;
      normal_cursor(); // depends on ctrl key for rotate
   }
}

// scriptable interface:
//
void
cis_trans_convert(int imol, const char *chain_id, int resno, const char *inscode) {

   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *standard_residues_mol = g.standard_residues_asc.mol;
      g.molecules[imol].cis_trans_conversion(chain_id, resno, inscode, standard_residues_mol);
      graphics_draw();
   }
}


/*  ----------------------------------------------------------------------- */
/*                  reverse direction                                       */
/*  ----------------------------------------------------------------------- */
void
setup_reverse_direction(short int istate) {

   graphics_info_t::in_reverse_direction_define = istate;
   graphics_info_t g;
   if (istate == 1) {
      g.pick_cursor_maybe();
      g.add_status_bar_text("Click on an atom in the fragment that you want to reverse");
      g.pick_pending_flag = 1;
   } else {
      g.normal_cursor();
   }

}

