/* src/c-interface.h
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009 by The University of Oxford
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

/* svn $Id: c-interface.h 1458 2007-01-26 20:20:18Z emsley $ */

/*! \file 
  \brief Coot Scripting Interface

  Here is a list of all the scripting interface functions. They are
  described/formatted in c/python format.

  Usually coot is compiled with the guile interpreter, and in this
  case these function names and usage are changed a little, e.g.:

  c-format:
  chain_n_residues("A", 1)

  scheme format:
  (chain-n-residues "A" 1)

  Note the prefix usage of the parenthesis and the lack of comma to
  separate the arguments.

*/

#ifndef C_INTERFACE_H
#define C_INTERFACE_H

/*
  The following extern stuff here because we want to return the
  filename from the file entry box.  That code (e.g.) 
  on_ok_button_coordinates_clicked (callback.c), is written and
  compiled in c.
 
  But, we need that function to set the filename in mol_info, which 
  is a c++ class.
 
  So we need to have this function external for c++ linking.
 
*/

/* Francois says move this up here so that things don't get wrapped
   twice in C-declarations inside gmp library. Hmm! */
#ifdef __cplusplus
#ifdef USE_GUILE
#include <cstdio> // for std::FILE in gmp.h for libguile.h
#include <libguile.h>		/* for SCM type (returned by safe_scheme_command) */
#endif /*  USE_GUILE */
// BL says:: ok then we put python here too
#ifdef USE_PYTHON
#include <Python.h>
#endif /*  PYTHON */
#endif /* c++ */

#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS extern
#define END_C_DECLS     
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS

/* Fix this on a rainy day. */
/* #ifdef __cplusplus */
/* #include "mini-mol.hh" */
/* #endif //  __cplusplus */

/* For similar reason to above, we can't have this here (mmdb poisoning) */
/* #ifdef __cplusplus */
/* #include "sequence-view.hh" */
/* #endif */

#define COOT_SCHEME_DIR "COOT_SCHEME_DIR"

/*  ------------------------------------------------------------------------ */
/*                         Startup Functions:                                */
/*  ------------------------------------------------------------------------ */
#ifdef USE_GUILE
void try_load_scheme_extras_dir();
#endif // USE_GUILE
#ifdef USE_PYTHON
void try_load_python_extras_dir();
#endif // USE_PYTHON

/*!  \brief tell coot that you prefer to run python scripts if/when
  there is an option to do so. */
void set_prefer_python();

/*  ------------------------------------------------------------------------ */
/*                         File system Functions:                            */
/*  ------------------------------------------------------------------------ */
/*  File system Utility function: maybe there is a better place for it... */

 
/*  Return like mkdir: mkdir returns zero on success, or -1 if an error */
/*  occurred */

/*  if it already exists as a dir, return 0 of course.  Perhaps should */
/*  be called "is_directory?_if_not_make_it". */

/* section File System Functions */
/*!  \name File System Functions */
/* \{ */

/*! \brief make a directory dir (if it doesn't exist) and return error code

   If it can be created, create the directory dir, return the success status
   like mkdir: mkdir 

   @return zero on success, or -1 if an  error  occurred.
   If dir already exists as a directory, return 0 of course.
 */
int make_directory_maybe(const char *dir);

/*! \brief Show Paths in Display Manager?

    Some people don't like to see the full path names in the display manager
    here is the way to turn them off, with an argument of 1.
*/
void set_show_paths_in_display_manager(int i);

/*! \brief return the internal state

   What is the internal flag? 

   @return 1 for "yes, display paths" , 0 for not
 */
int show_paths_in_display_manager_state();

/*! \brief add an extension to be treated as coordinate files 
*/
void add_coordinates_glob_extension(const char *ext);

/*! \brief add an extension to be treated as data (reflection) files 
*/
void add_data_glob_extension(const char *ext);

/*! \brief add an extension to be treated as geometry dictionary files 
*/
void add_dictionary_glob_extension(const char *ext);

/*! \brief add an extension to be treated as geometry map files 
*/
void add_map_glob_extension(const char *ext);

/*! \brief remove an extension to be treated as coordinate files 
*/
void remove_coordinates_glob_extension(const char *ext);

/*! \brief remove an extension to be treated as data (reflection) files 
*/
void remove_data_glob_extension(const char *ext);

/*! \brief remove an extension to be treated as geometry dictionary files 
*/
void remove_dictionary_glob_extension(const char *ext);

/*! \brief remove an extension to be treated as geometry map files 
*/
void remove_map_glob_extension(const char *ext);

/*! \brief sort files in the file selection by date?

  some people like to have their files sorted by date by default */
void set_sticky_sort_by_date(); 

/*! \brief do not sort files in the file selection by date?

  removes the sorting of files by date */
void unset_sticky_sort_by_date(); 

/*! \brief on opening a file selection dialog, pre-filter the files.

set to 1 to pre-filter, [0 (off, non-pre-filtering) is the default */
void set_filter_fileselection_filenames(int istate);

void set_file_selection_dialog_size(GtkWidget *w);

/*! \brief, return the state of the above variable */
int filter_fileselection_filenames_state();


/*! \brief display the open coordinates dialog */
void open_coords_dialog();

#ifdef COOT_USE_GTK2_INTERFACE
void on_filename_filter_toggle_button_toggled (GtkButton       *button,
					      gpointer         user_data);
#else
void on_filename_filter_toggle_button_toggled_gtk1 (GtkButton       *button,
						    gpointer         user_data);
#endif 
void add_filename_filter(GtkWidget *fileselection);


// where data type:
// 0 coords
// 1 mtz etc
// 2 maps
// (return the button)
GtkWidget *add_filename_filter_button(GtkWidget *fileselection, 
				      short int type);

#ifdef COOT_USE_GTK2_INTERFACE
void add_filechooser_filter_button(GtkWidget *fileselection, 
				      short int data_type);

void add_filechooser_extra_filter_button(GtkWidget *fileselection, 
				      const gchar *name,
                                      const gchar *name2);
#endif // GTK2

gboolean on_filename_filter_key_press_event (GtkWidget       *widget,
					     GdkEventKey     *event,
					     gpointer         user_data);

/* a c callable wrapper to the graphics_info_t function */
void fill_option_menu_with_coordinates_options(GtkWidget *option_menu, 
					       GtkSignalFunc signal_func,
					       int imol_active_position);
void fill_option_menu_with_coordinates_options_unsaved_first(GtkWidget *option_menu, 
							     GtkSignalFunc signal_func,
							     int imol_active_position);
GtkWidget *coot_file_chooser();

GtkWidget *coot_dataset_chooser();

GtkWidget *coot_map_name_chooser();

GtkWidget *coot_save_coords_chooser();

GtkWidget *coot_cif_dictionary_chooser();

GtkWidget *coot_run_script_chooser();

GtkWidget *coot_save_state_chooser();

GtkWidget *coot_save_symmetry_chooser();

GtkWidget *coot_screendump_chooser();

void set_directory_for_coot_file_chooser(GtkWidget *w);

const char *coot_file_chooser_file_name(GtkWidget *widget);

void set_filename_for_filechooserselection(GtkWidget *widget, const gchar *name);

/* some BL functions for gtk2 */

/*! \brief this flag set chooser as default for windows, otherwise use
  selector 0 is selector 1 is chooser */
#ifdef COOT_USE_GTK2_INTERFACE
void set_file_chooser_selector(int istate);
int file_chooser_selector_state();

void set_file_chooser_overwrite(int istate);
int file_chooser_overwrite_state();
#endif
/* \} */


/*  ---------------------------------------------------------------------- */
/*                     widget utilities                                    */
/*  ---------------------------------------------------------------------- */
/* section Widget Utilities */
/*! \name Widget Utilities */
/* \{ */
/* return negative if fail */
float get_positive_float_from_entry(GtkEntry *w); 

#ifdef COOT_USE_GTK2_INTERFACE
void handle_filename_filter_gtk2(GtkWidget *widget);
#else 
void handle_filename_filter_gtk1(GtkWidget *widget);
#endif 

void set_transient_and_position(int window_type, GtkWidget *window);

/*! \brief create a dialog with information

  create a dialog with information string txt.  User has to click to
  dismiss it, but it is not modal (nothing in coot is modal). */
void info_dialog(const char *txt); 

GtkWidget *main_menubar();
GtkWidget *main_statusbar();
GtkWidget *main_toolbar();

/* \} */

/*  -------------------------------------------------------------------- */
/*                   mtz and data handling utilities                     */
/*  -------------------------------------------------------------------- */
/* section MTZ and data handling utilities */
/*! \name  MTZ and data handling utilities */
/* \{ */
/* We try as .phs and .cif files first */

/*! \brief given a filename, try to read it as a data file 

   We try as .phs and .cif files first */
void manage_column_selector(const char *filename);
void manage_refmac_column_selection(GtkWidget *w);
void fill_f_optionmenu_with_expert_options(GtkWidget *f_optionmenu);
void handle_column_label_make_fourier(GtkWidget *column_label_window);
void wrapped_create_run_refmac_dialog();

/* \} */

/*  -------------------------------------------------------------------- */
/*                     Molecule Functions       :                        */
/*  -------------------------------------------------------------------- */
/* section Molecule Info Functions */
/*! \name Molecule Info Functions */
/* \{ */

/*! \brief  the number of residues in chain chain_id and molecule number imol
  @return the number of residues 
*/
int chain_n_residues(const char *chain_id, int imol); 
/* @return status, less than -9999 is for failure (eg. bad imol); */
float molecule_centre_internal(int imol, int iaxis);
/*! \brief return the rename from a residue serial number

   @return NULL (scheme False) on failure. */
char *resname_from_serial_number(int imol, const char *chain_id, 
				 int serial_num);

/*! \brief a residue seqnum (normal residue number) from a residue
  serial number

   @return < -9999 on failure */
int  seqnum_from_serial_number(int imol, const char *chain_id, 
			       int serial_num);

/*! \brief the insertion code of the residue.

   @return NULL (scheme False) on failure. */
char *insertion_code_from_serial_number(int imol, const char *chain_id, int serial_num);

/*! \brief the chain_id (string) of the ichain-th chain
  molecule number imol  
   @return the chain-id */
/* char *chain_id(int imol, int ichain); */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM
chain_id_scm(int imol, int ichain); 
#endif
#ifdef USE_PYTHON
PyObject *
chain_id_py(int imol, int ichain); 
#endif
#endif


/*! \brief  number of chains in molecule number imol 

   @return the number of chains*/
int n_chains(int imol);

/*! \brief is this a solvent chain? [Raw function]

   This is a raw interface function, you should generally not use
   this, but instead use (is-solvent-chain? imol chain-id)

   @return -1 on error, 0 for no, 1 for is "a solvent chain".  We
   wouldn't want to be doing rotamer searches and the like on such a
   chain.
 */
int is_solvent_chain_p(int imol, const char *chain_id);

/*! \brief sort the chain ids of the imol-th molecule in lexographical order */
void sort_chains(int imol);	


/*! \brief simply print secondardy structure info to the
  terminal/console.  In future, this could/should return the info.  */
void print_header_secondary_structure_info(int imol);


/*  Placeholder only.

    not documented, it doesn't work yet, because CalcSecStructure()
    creates SS type on the residues, it does not build and store
    CHelix, CStrand, CSheet records. */
void write_header_secondary_structure_info(int imol, const char *file_name);


/*! \brief copy molecule imol

@return the new molecule number.
Return -1 on failure to copy molecule (out of range, or molecule is
closed) */
int copy_molecule(int imol);

/*! \brief Experimental interface for Ribosome People. 

Ribosome People have many chains in their pdb file, they prefer segids
to chainids (chainids are only 1 character).  But coot uses the
concept of chain ids and not seg-ids.  mmdb allow us to use more than
one char in the chainid, so after we read in a pdb, let's replace the
chain ids with the segids. Will that help? */
int exchange_chain_ids_for_seg_ids(int imol);

/* \} */

/*  -------------------------------------------------------------------- */
/*                     Library/Utility Functions:                        */
/*  -------------------------------------------------------------------- */

/* section Library and Utility Functions */
/*! \name Library and Utility Functions */
/* \{ */

/*! \brief the coot version string 

   @return something like "coot-0.1.3".  New versions of coot will
   always be lexographically greater than previous versions. */
char *coot_version();

#ifdef __cplusplus
#ifdef USE_GUILE
SCM coot_sys_build_type_scm();
#endif
#ifdef USE_PYTHON
PyObject *coot_sys_build_type_py();
#endif /* USE_PYTHON */
#endif /* c++ */

/*! \brief return the subversion revision number of this build.
 
Used in finding updates.  */
int svn_revision(); 


/*! \brief return the name of molecule number imol

 @return 0 if not a valid name ( -> False in scheme) 
 e.g. "/a/b/c.pdb" for "d/e/f.mtz FWT PHWT" */ 
const char *molecule_name(int imol);
/*! \brief set the molecule name of the imol-th molecule */
void set_molecule_name(int imol, const char *new_name);
GtkWidget *main_window(); 
gboolean coot_checked_exit(int retval); 
/*! \brief exit from coot, give return value retval back to invoking
  process. */
void coot_real_exit(int retval); 
void coot_clear_backup_or_real_exit(int retval);

#ifdef __cplusplus
#ifdef USE_GUILE
void run_clear_backups(int retval);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
void run_clear_backups_py(int retval);
#endif /* USE_PYTHON */
#endif  /* c++ */

void fill_about_window(GtkWidget *widget);
void add_coot_references_button(GtkWidget *widget);
GtkWidget *wrapped_create_coot_references_dialog();
#ifdef COOT_USE_GTK2_INTERFACE
void fill_references_notebook(GtkToolButton *toolbutton, int reference_id);
#endif
 
/*! \brief What is the molecule number of first coordinates molecule?

   return -1 when there is none. */
int first_coords_imol(); 	

/*! \brief What is the molecule number of first unsaved coordinates molecule?

   return -1 when there is none. */
int first_unsaved_coords_imol(); 	

/* \} */

/*  -------------------------------------------------------------------- */
/*                    More Library/Utility Functions:                    */
/*  -------------------------------------------------------------------- */
/* section Graphics Utility Functions */
/*! \name Graphics Utility Functions */
/* \{ */

/*! \brief set the bond lines to be antialiased */
void set_do_anti_aliasing(int state);
/*! \brief return the flag for antialiasing the bond lines */
int do_anti_aliasing_state();

/*! \brief turn the GL lighting on (state = 1) or off (state = 0) 

   slows down the display of simple lines
*/
void set_do_GL_lighting(int state);
/*! \brief return the flag for GL lighting */
int do_GL_lighting_state();

/*! \brief shall we start up the Gtk and the graphics window? 

   if passed the command line argument --no-graphics, coot will not start up gtk
   itself.

   An interface function for Ralf.
*/
short int use_graphics_interface_state();
short int python_at_prompt_at_startup_state();

/*! \brief start Gtk (and graphics) 

   This function is useful if it was not started already (which can be
   achieved by using the command line argument --no-graphics).

   An interface for Ralf */
void start_graphics_interface(); 

/*! \brief "Reset" the view
 

  return 1 if we moved, else return 0.

   centre on last-read molecule with zoom 100. If we are there, then
   go to the previous molecule, if we are there, then go to the origin. */
int reset_view();

/*! \brief return the number of molecules (coordinates molecules and
  map molecules combined) that are currently in coot

  @return the number of molecules (closed molecules are not counted) */
int graphics_n_molecules(); 

/* return either 1 (yes, there is at least one hydrogen) or 0 (no
   hydrogens, or no such molecule).  The scripting interface to this
   does not have the _raw suffix and returns a scheme or python
   boolean True or False.  */
int molecule_has_hydrogens_raw(int imol); 

/* a testing/debugging function.  Used in a test to make sure that the
   outside number of a molecule (the vector index) is the same as that
   embedded in the molecule description object.  Return -1 on
   non-valid passed imol. */
int own_molecule_number(int imol);

int next_map_for_molecule(int imol); /* return a map number */

/*! \brief Spin spin spin (or not) */
void toggle_idle_spin_function(); 

/*! \brief how far should we rotate when (auto) spinning? Fast
  computer? set this to 0.1  */
void set_idle_function_rotate_angle(float f);  // degrees

float idle_function_rotate_angle();

/* pass back the newly created molecule number */
/*! \brief a synonym for read-pdb.  Read the coordinates from
  filename (can be pdb, cif or shelx format)  */
int handle_read_draw_molecule(const char *filename);

/*! \brief read coordinates from filename with option to not recentre.

   set recentre_on_read_pdb_flag to 0 if you don't want the view to
   recentre on the new coordinates. */
int handle_read_draw_molecule_with_recentre(const char *filename, 
					    int recentre_on_read_pdb_flag);

/*! \brief read coordinates from filename and recentre the new
  molecule at the screen rotation centre. */
int handle_read_draw_molecule_and_move_molecule_here(const char *filename);

/*! \brief read coordinates from filename */
int read_pdb(const char *filename); // cc4mg function name

/*! \brief some programs produce PDB files with ATOMs where there
  should be HETATMs.  This is a function to assign HETATMs as per the
  PDB definition. */
int assign_hetatms(int imol);

/*! \brief replace the parts of molecule number imol that are
  duplicated in molecule number imol_frag */
int replace_fragment(int imol_target, int imol_fragment, const char *atom_selection);

/*! \brief replace pdb.  Fail if molecule_number is not a valid model molecule.
  Return -1 on failure.  Else return molecule_number  */
int clear_and_update_model_molecule_from_file(int molecule_number, 
					      const char *file_name);

/* Used in execute_rigid_body_refine */
/* Fix this on a rainy day. */
/* atom_selection_container_t  */
/* make_atom_selection(int imol, const coot::minimol::molecule &mol);  */

/*! \brief dump the current screen image to a file.  Format ppm 

You can use this, in conjunction with spinning and view moving functions to 
make movies */
void screendump_image(const char *filename);

void add_is_difference_map_checkbutton(GtkWidget *fileselection); 
/* the callback for the above: */
void
on_read_map_difference_map_toggle_button_toggled (GtkButton       *button,
						  gpointer         user_data);

void add_recentre_on_read_pdb_checkbutton(GtkWidget *fileselection); 
/* the callback for the above: */
void
on_recentre_on_read_pdb_toggle_button_toggled (GtkButton       *button,
					       gpointer         user_data);

/* \} */

/*  -------------------------------------------------------------------- */
/*                     Testing Interface:                                */
/*  -------------------------------------------------------------------- */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM test_internal_scm(); 
SCM test_internal_single_scm(); 
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
PyObject *test_internal_py(); 
PyObject *test_internal_single_py(); 
#endif	/* USE_PYTHON */
#endif	/* __cplusplus */


/*  --------------------------------------------------------------------- */
/*                      Interface Preferences                             */
/*  --------------------------------------------------------------------- */
/* section Interface Preferences */
/*! \name   Interface Preferences */
/*! \{ */

/*! \brief Some people (like Phil Evans) don't want to scroll their
  map with the mouse-wheel.

  To turn off mouse wheel recontouring call this with istate value of 0  */
void set_scroll_by_wheel_mouse(int istate);
/*! \brief return the internal state of the scroll-wheel map contouring */
int scroll_by_wheel_mouse_state();

/*! \brief set the default inital contour for 2FoFc-style map

in sigma */
void set_default_initial_contour_level_for_map(float n_sigma);

/*! \brief set the default inital contour for FoFc-style map

in sigma */
void set_default_initial_contour_level_for_difference_map(float n_sigma);

/*! \brief print the view matrix to the console, useful for molscript,
  perhaps */
void print_view_matrix();		/* print the view matrix */

float get_view_matrix_element(int row, int col); /* used in (view-matrix) command */

/*! \brief internal function to get an element of the view quaternion.
  The whole quaternion is returned by the scheme function
  view-quaternion  */
float get_view_quaternion_internal(int element);

/*! \brief Set the view quaternion */
void set_view_quaternion(float i, float j, float k, float l);

/*! \brief Given that we are in chain current_chain, apply the NCS
  operator that maps current_chain on to next_ncs_chain, so that the
  relative view is preserved.  For NCS skipping. */
void apply_ncs_to_view_orientation(int imol, const char *current_chain, const char *next_ncs_chain);
/*! \brief as above, but shift the screen centre also.  */
void apply_ncs_to_view_orientation_and_screen_centre(int imol, 
						     const char *current_chain, 
						     const char *next_ncs_chain);

void set_fps_flag(int t);
int  get_fps_flag();

/*! \brief set a flag: is the origin marker to be shown? 1 for yes, 0
  for no. */
void set_show_origin_marker(int istate);
/*! \brief return the origin marker shown? state */
int  show_origin_marker_state();

/*! \brief hide the vertical modelling toolbar in the GTK2 version */
void hide_modelling_toolbar();
/*! \brief show the vertical modelling toolbar in the GTK2 version
  (the toolbar is shown by default) */
void show_modelling_toolbar();

/*! \brief show all available icons in the modelling toolbar (same as MFR dialog) */
void show_model_toolbar_all_icons();
/*! \brief show only a selection of icons in the modelling toolbar */
void show_model_toolbar_main_icons();

void toolbar_popup_menu(GtkToolbar *toolbar, 
		    GdkEventButton *event_button,
		    gpointer user_data);

void set_model_toolbar_docked_position_callback(GtkWidget *w, gpointer user_data);

/*! \brief reattach the modelling toolbar to the last attached position */
void reattach_modelling_toolbar();

/*! \brief to swap sides of the Model/Fit/Refine toolbar
  0 (default) is right, 1 is left, 2 is top, 3 is bottom */
void set_model_toolbar_docked_position(int state);

/*! \brief reparent the Model/Fit/Refine dialog so that it becomes
  part of the main window, next to the GL graphics context */
int suck_model_fit_dialog();
int suck_model_fit_dialog_bl();

/* return the dialog if it exists, else null */
GtkWidget *close_model_fit_dialog(GtkWidget *dialog_hbox);
/* use this from the scripting layer to say something to the user (popup). */
GtkWidget *popup_window(const char *s);

/*! \brief Put text s into the status bar.

  use this to put info for the user in the statusbar (less intrusive
  than popup). */
void add_status_bar_text(const char *s);

void set_model_fit_refine_dialog_stays_on_top(int istate);
int model_fit_refine_dialog_stays_on_top_state();

void save_accept_reject_dialog_window_position(GtkWidget *acc_reg_dialog);
void set_accept_reject_dialog(GtkWidget *w); /* used by callbacks to unset the widget */

/* functions to dock the accept/reject dialog to the toolbar */
void set_accept_reject_dialog_docked(int state);
int accept_reject_dialog_docked_state();

/* functions to show/hide i.e. make sensitive the docked accept/reject toolbar */
void set_accept_reject_dialog_docked_show(int state);
int accept_reject_dialog_docked_show_state();

/* functions for the refinement toolbar style */
void set_model_toolbar_style(int state);
int model_toolbar_style_state();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                           mouse buttons                                  */
/*  ----------------------------------------------------------------------- */
/* section Mouse Buttons */
/*! \name   Mouse Buttons */
/* \{ */

/* Note, when you have set these, there is no way to turn them of
   again (other than restarting). */
void quanta_buttons(); 
void quanta_like_zoom();


/* -------------------------------------------------------------------- */
/*    Ctrl for rotate or pick: */
/* -------------------------------------------------------------------- */
/*! \brief Alternate mode for rotation

Prefered by some, including Dirk Kostrewa.  I don't think this mode
works properly yet */
void set_control_key_for_rotate(int state);
/*! \brief return the control key rotate state */
int control_key_for_rotate_state();

/*! \brief Put the blob under the cursor to the screen centre.  Check only
positive blobs.  Useful function if bound to a key.

The refinement map must be set.  (We can't check all maps because they
are not (or may not be) on the same scale).

   @return 1 if successfully found a blob and moved there.
   return 0 if no move.
*/
int blob_under_pointer_to_screen_centre();

/* \} */

/*  --------------------------------------------------------------------- */
/*                      Cursor Functions:                                 */
/*  --------------------------------------------------------------------- */
/* section Cursor Function */
/*! \name Cursor Function */
/*! \{ */
void normal_cursor();
void fleur_cursor();
void pick_cursor_maybe();
void rotate_cursor();

/*! \brief let the user have a different pick cursor

sometimes (the default) GDK_CROSSHAIR is hard to see, let the user set
their own */
void set_pick_cursor_index(int icursor_index);

/*! \} */


/*  --------------------------------------------------------------------- */
/*                      Model/Fit/Refine Functions:                       */
/*  --------------------------------------------------------------------- */
/* section Model/Fit/Refine Functions  */
/*! \name Model/Fit/Refine Functions  */
/* \{ */
/*! \brief display the Model/Fit/Refine dialog */
void post_model_fit_refine_dialog();
GtkWidget *wrapped_create_model_fit_refine_dialog(); 
void update_model_fit_refine_dialog_menu(GtkWidget *widget);
void update_model_fit_refine_dialog_buttons(GtkWidget *widget);
void unset_model_fit_refine_dialog();
void unset_refine_params_dialog();
/*! \brief display the Display Manager dialog */
void show_select_map_dialog();
/*! \brief Allow the changing of Model/Fit/Refine button label from
  "Rotate/Translate Zone" */
void set_model_fit_refine_rotate_translate_zone_label(const char *txt);
/*! \brief Allow the changing of Model/Fit/Refine button label from
  "Place Atom at Pointer" */
void set_model_fit_refine_place_atom_at_pointer_label(const char *txt);

/* other tools */
GtkWidget *wrapped_create_other_model_tools_dialog();
void unset_other_modelling_tools_dialog();

/*! \brief display the Other Modelling Tools dialog */
void post_other_modelling_tools_dialog();

/*! \brief shall atoms with zero occupancy be moved when refining? (default 1, yes) */
void set_refinement_move_atoms_with_zero_occupancy(int state);
/*! \brief return the state of "shall atoms with zero occupancy be moved when refining?" */
int refinement_move_atoms_with_zero_occupancy_state();

GtkWidget *wrapped_create_fast_ss_search_dialog();

/*! \} */

/*  --------------------------------------------------------------------- */
/*                      backup/undo functions:                            */
/*  --------------------------------------------------------------------- */
/* section Backup Functions */
/*! \name Backup Functions */
/* \{ */
/* c-interface-build functions */

/*! \brief make backup for molecule number imol */
void make_backup(int imol);

/*! \brief turn off backups for molecule number imol */
void turn_off_backup(int imol); 
/*! \brief turn on backups for molecule number imol */
void turn_on_backup(int imol);
/*! \brief return the backup state for molecule number imol 

 return 0 for backups off, 1 for backups on, -1 for unknown */
int  backup_state(int imol);
void apply_undo();		/* "Undo" button callback */
void apply_redo();

/*! \brief set the molecule number imol to be marked as having unsaved changes */
void set_have_unsaved_changes(int imol);

/*! \brief does molecule number imol have unsaved changes?
 @return -1 on bad imol, 0 on no unsaved changes, 1 on has unsaved changes */
int have_unsaved_changes_p(int imol);


/*! \brief set the molecule to which undo operations are done to
  molecule number imol */
void set_undo_molecule(int imol);

/*! \brief show the Undo Molecule chooser - i.e. choose the molecule
  to which the "Undo" button applies. */
void show_set_undo_molecule_chooser();

/* called by above */
GtkWidget *wrapped_create_undo_molecule_chooser_dialog();

/*! \brief set the state for adding paths to backup file names

  by default directories names are added into the filename for backup
  (with / to _ mapping).  call this with state=1 to turn off directory
  names  */
void set_unpathed_backup_file_names(int state);
/*! \brief return the state for adding paths to backup file names*/
int  unpathed_backup_file_names_state();
/* \} */

/*  --------------------------------------------------------------------- */
/*                         recover session:                               */
/*  --------------------------------------------------------------------- */
/* section Recover Session Function */
/*! \name  Recover Session Function */
/* \{ */
/*! \brief recover session

   After a crash, we provide this convenient interface to restore the
   session.  It runs through all the molecules with models and looks
   at the coot backup directory looking for related backup files that
   are more recent that the read file. (Not very good, because you
   need to remember which files you read in before the crash - should
   be improved.) */
void recover_session();
void execute_recover_session(GtkWidget *w);
/* \} */
   
/*  ---------------------------------------------------------------------- */
/*                       map functions:                                    */
/*  ---------------------------------------------------------------------- */
/* section Map Functions */
/*! \name  Map Functions */
/* \{ */

/*! \brief fire up a GUI, which asks us which model molecule we want
  to calc phases from.  On "OK" button there, we call
  map_from_mtz_by_refmac_calc_phases() */
void calc_phases_generic(const char *mtz_file_name);

/*! \brief Calculate SFs (using refmac optionally) from an MTZ file
  and generate a map. Get F and SIGF automatically (first of their
  type) from the mtz file.

@return the new molecule number, -1 on a problem. */
int map_from_mtz_by_refmac_calc_phases(const char *mtz_file_name, 
				       const char *f_col, 
				       const char *sigf_col, 
				       int imol_coords);


/*! \brief Calculate SFs from an MTZ file and generate a map. 
 @return the new molecule number. */
int map_from_mtz_by_calc_phases(const char *mtz_file_name, 
				const char *f_col, 
				const char *sigf_col, 
				int imol_coords);

gdouble* get_map_colour(int imol);

void add_on_map_colour_choices(GtkWidget *w);

/* the callback set on the submenu items in the above function */
void map_colour_mol_selector_activate (GtkMenuItem     *menuitem,
				       gpointer         user_data);
void my_delete_menu_items(GtkWidget *widget, void *data);

/* similarly for the scrollwheel */
void add_on_map_scroll_whell_choices(GtkWidget *menu);
void map_scroll_wheel_mol_selector_activate (GtkMenuItem     *menuitem,
					     gpointer         user_data);

/*! \brief set the map that is moved by changing the scroll wheel and
  change_contour_level(). */
void set_scroll_wheel_map(int imap);
/*! \brief return the molecule number to which the mouse scroll wheel
  is attached */
int scroll_wheel_map();
/*! \brief save previous colour map for molecule number imol */
void save_previous_map_colour(int imol);
/*! \brief restore previous colour map for molecule number imol */
void restore_previous_map_colour(int imol);

/*! \brief set the state of immediate map upate on map drag.

By default, it is on (t=1).  On slower computers it might be better to
set t=0. */
void set_active_map_drag_flag(int t);
/*! \brief return the state of the dragged map flag  */
short int get_active_map_drag_flag(); 

/*! \brief set the colour of the last (highest molecule number) map */
void set_last_map_colour(double f1, double f2, double f3);
/*! \brief set the colour of the imolth map */
void set_map_colour(int imol, float red, float green, float blue);

void handle_map_colour_change     (int map_no, gdouble[4]);
void handle_symmetry_colour_change(int mol,    gdouble[4]);
void fill_single_map_properties_dialog(GtkWidget *window, int imol);

void set_contour_level_absolute(int imol_map, float level);
void set_contour_level_in_sigma(int imol_map, float level);

/*! \brief set the sigma step of the last map to f sigma */
void set_last_map_sigma_step(float f);
void set_contour_sigma_button_and_entry(GtkWidget *window, int imol);
void set_contour_by_sigma_step_maybe(GtkWidget *window, int imol);
/*! \brief set the contour level step

   set the contour level step of molecule number imol to f and
   variable state (setting state to 0 turns off contouring by sigma
   level)  */
void set_contour_by_sigma_step_by_mol(float f, short int state, int imol); 

/*! \brief return the resolution of the data for molecule number imol.
   Return negative number on error, otherwise resolution in A (eg. 2.0) */
float data_resolution(int imol);

void solid_surface(int imap, short int on_off_flag);

/*! \brief export (write to disk) the map of molecule number imol to
  filename.  
  
  Return 0 on failure, 1 on success. */
int export_map(int imol, const char *filename);

/* return the new molecule number. The cell is given in Angstroms and
   the angles in degrees.  The ref_space_group can be a H-M symbol or
   a colon-separated string of symmetry operators.
*/
int transform_map_raw(int imol, 
		      double r00, double r01, double r02, 
		      double r10, double r11, double r12, 
		      double r20, double r21, double r22, 
		      double t0, double t1, double t2, 
		      double pt0, double pt1, double pt2, 
		      double box_half_size, 		     
		      const char *ref_space_group,
		      double cell_a, double cell_b, double cell_c,
		      double alpha, double beta, double gamma);

/* return the new molecule number (or should I do it in place?) */
int rotate_map_round_screen_axis_x(float r_degrees); 
int rotate_map_round_screen_axis_y(float r_degrees); 
int rotate_map_round_screen_axis_z(float r_degrees); 

/*! \brief make a difference map, taking map_scale * imap2 from imap1,
  on the grid of imap1.  Return the new molecule number.  
  Return -1 on failure. */
int difference_map(int imol1, int imol2, float map_scale);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief make an average map from the map_number_and_scales (which
  is a list of pairs (list map-number scale-factor)) (the scale
  factors are typically 1.0 of course). The output map is in the same
  grid as the first (valid) map.  Return -1 on failure to make an
  averaged map, otherwise return the new map molecule number. */
int average_map_scm(SCM map_number_and_scales);
#endif 
#ifdef USE_PYTHON
/*! \brief make an average map from the map_number_and_scales (which
  is a list of pairs [map_number, scale_factor] (the scale factors
  are typically 1.0 of course). The output map is in the same
  grid as the first (valid) map.  Return -1 on failure to make an
  averaged map, otherwise return the new map molecule number. */
int average_map_py(PyObject *map_number_and_scales);
#endif /* USE_PYTHON */
#endif /* c++ */

/* \} */

/*  ----------------------------------------------------------------------- */
/*                         (density) iso level increment entry */
/*  ----------------------------------------------------------------------- */
/* section Density Increment */
/*! \name  Density Increment */
/* \{ */

char* get_text_for_iso_level_increment_entry(int imol); /* const gchar *text */
char* get_text_for_diff_map_iso_level_increment_entry(int imol); /* const gchar *text */

/* void set_iso_level_increment(float val); */
/*! \brief set the contour scroll step (in absolute e/A3) for
  2Fo-Fc-style maps to val

The is only activated when scrolling by sigma is turned off */
void set_iso_level_increment(float val);
float get_iso_level_increment();
void set_iso_level_increment_from_text(const char *text, int imol);

/*! \brief set the contour scroll step for difference map (in absolute
  e/A3) to val

The is only activated when scrolling by sigma is turned off */
void set_diff_map_iso_level_increment(float val);
float get_diff_map_iso_level_increment();
void set_diff_map_iso_level_increment_from_text(const char *text, int imol);

/*  find the molecule that the single map dialog applies to and set
    the contour level and redraw */
void single_map_properties_apply_contour_level_to_map(GtkWidget *w);

void set_map_sampling_rate_text(const char *text);

/*! \brief set the map sampling rate (default 1.5)

Set to something like 2.0 or 2.5 for more finely sampled maps.  Useful
for baton-building low resolution maps. */
void set_map_sampling_rate(float r);
char* get_text_for_map_sampling_rate_text();

/*! \brief return the map sampling rate */
float get_map_sampling_rate();

/*! \brief set the map that has its contour level changed by the
  scrolling the mouse wheel to molecule number imol */
void set_scrollable_map(int imol); 
/*! \brief change the contour level of the current map by a step

if is_increment=1 the contour level is increased.  If is_increment=0
the map contour level is decreased.
 */
void change_contour_level(short int is_increment); /* else is decrement.  */

/*! \brief set the contour level of the map with the highest molecule
    number to level */
void set_last_map_contour_level(float level);
/*! \brief set the contour level of the map with the highest molecule
    number to n_sigma sigma */
void set_last_map_contour_level_by_sigma(float n_sigma);

/*! \brief create a lower limit to the "Fo-Fc-style" map contour level changing 

  (default 1 on) */
void set_stop_scroll_diff_map(int i);
/*! \brief create a lower limit to the "2Fo-Fc-style" map contour level changing 

  (default 1 on) */
void set_stop_scroll_iso_map(int i);

/*! \brief set the actual map level changing limit 

   (default 0.0) */
void set_stop_scroll_iso_map_level(float f); 

/*! \brief set the actual difference map level changing limit 

   (default 0.0) */
void set_stop_scroll_diff_map_level(float f);

/*! \brief set the scale factor for the Residue Density fit analysis */
void set_residue_density_fit_scale_factor(float f);


/*! \} */

/*  ------------------------------------------------------------------------ */
/*                         density stuff                                     */
/*  ------------------------------------------------------------------------ */
/* section Density Functions */
/*! \name  Density Functions */
/*! \{ */
/*! \brief draw the lines of the chickenwire density in width w */
void set_map_line_width(int w);
/*! \brief return the width in which density contours are drawn */
int map_line_width_state();

/*! \brief make a map from an mtz file (simple interface)

 given mtz file mtz_file_name and F column f_col and phases column
 phi_col and optional weight column weight_col (pass use_weights=0 if
 weights are not to be used).  Also mark the map as a difference map
 (is_diff_map=1) or not (is_diff_map=0) because they are handled
 differently inside coot.

 @return -1 on error, else return the new molecule number */
int make_and_draw_map(const char *mtz_file_name, 
		      const char *f_col, const char *phi_col, 
		      const char *weight,
		      int use_weights, int is_diff_map); 

/*! \brief as the above function, execpt set refmac parameters too

 pass along the refmac column labels for storage (not used in the
 creation of the map)

 @return -1 on error, else return imol */
int  make_and_draw_map_with_refmac_params(const char *mtz_file_name, 
		       const char *a, const char *b, const char *weight,
					  int use_weights, int is_diff_map,
					  short int have_refmac_params,
					  const char *fobs_col,
					  const char *sigfobs_col,
					  const char *r_free_col,
					  short int sensible_f_free_col);

/*! \brief as the above function, except set expert options too.

*/
/* Note to self, we need to save the reso limits in the state file  */
int make_and_draw_map_with_reso_with_refmac_params(const char *mtz_file_name, 
						   const char *a, const char *b, 
						   const char *weight,
						   int use_weights, int is_diff_map,
						   short int have_refmac_params,
						   const char *fobs_col,
						   const char *sigfobs_col,
						   const char *r_free_col,
						   short int sensible_f_free_col,
						   short int is_anomalous,
						   short int use_reso_limits,
						   float low_reso_limit,
						   float high_reso_lim);

#ifdef __cplusplus
#ifdef USE_GUILE
SCM refmac_parameters_scm(int imol);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *refmac_parameters_py(int imol);
#endif	/* USE_PYTHON */

#endif	/* __cplusplus */


/*! \brief does the mtz file have the columms that we want it to have? */
int valid_labels(const char *mtz_file_name, const char *f_col, 
		 const char *phi_col, 
		 const char *weight_col, 
		 int use_weights);

/* We need to know if an mtz file has phases.  If it doesn't then we */
/*  go down a (new 20060920) different path. */
int mtz_file_has_phases_p(const char *mtz_file_name); 

int is_mtz_file_p(const char *filename);

int cns_file_has_phases_p(const char *cns_file_name); 

/*! \brief read MTZ file filename and from it try to make maps

Useful for reading the output of refmac.  The default labels (FWT/PHWT
and DELFWT/PHDELFWT) can be changed using ...[something] 

 @return the molecule number for the new map 
*/
int auto_read_make_and_draw_maps(const char *filename); 
/*! \brief set the flag to do a difference map (too) on auto-read MTZ */
void set_auto_read_do_difference_map_too(int i);
/*! \brief return the flag to do a difference map (too) on auto-read MTZ 

   @return 0 means no, 1 means yes. */
int auto_read_do_difference_map_too_state();
/*! \brief set the exected MTZ columns for Auto-reading MTZ file. 

  Not every program uses the default refmac labels (FWT/PHWT) for its
  MTZ file.  Here we can tell coot to expect other labels, 

  e.g. (set-auto-read-column-labels "2FOFCWT" "PH2FOFCWT" 0) */
void set_auto_read_column_labels(const char *fwt, const char *phwt, 
				 int is_for_diff_map_flag);

char* get_text_for_density_size_widget(); /* const gchar *text */
void set_density_size_from_widget(const char *text);

/*! \brief set the extent of the box/radius of electron density contours */
void set_map_radius(float f);

/*! \brief another (old) way of setting the radius of the map */
void set_density_size(float f);

void set_map_radius_slider_max(float f); 

/*! \brief Give me this nice message str when I start coot */
void set_display_intro_string(const char *str);

/*! \brief return the extent of the box/radius of electron density contours */
float get_map_radius();

/*! \brief not everone likes coot's esoteric depth cueing system

  Pass an argument istate=1 to turn it off

 (this function is currently disabled). */
void set_esoteric_depth_cue(int istate);

/*! \brief native depth cueing system

  return the state of the esoteric depth cueing flag */
int  esoteric_depth_cue_state();

/*! \brief not everone lies coot's default difference map colouring.  

   Pass an argument i=1 to swap the difference map colouring so that
   red is positve and green is negative. */
void set_swap_difference_map_colours(int i);
int swap_difference_map_colours_state();
/*! \brief post-hoc set the map of molecule number imol to be a
  difference map
  @return success status, 0 -> failure (imol does not have a map) */
int set_map_is_difference_map(int imol);

int map_is_difference_map(int imol);

/*! \brief Add another contour level for the last added map.  

  Currently, the map must have been generated from an MTZ file.
  @return the molecule number of the new molecule or -1 on failure */
int another_level();

/*! \brief Add another contour level for the given map.

  Currently, the map must have been generated from an MTZ file.
  @return the molecule number of the new molecule or -1 on failure */
int another_level_from_map_molecule_number(int imap);

/*! \brief return the scale factor for the Residue Density fit analysis */
float residue_density_fit_scale_factor();

/*! \brief return the density at the given point for the given
  map. Return 0 for bad imol */
float density_at_point(int imol, float x, float y, float z);


#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return sigma for the given map.  Return scheme False if not
  a valid map molecule number. */
SCM map_sigma_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief return sigma for the given map.  Return Python False if not
  a valid map molecule number. */
PyObject *map_sigma_py(int imol);
#endif /*USE_PYTHON */
#endif  /* c++ */

/* \} */
 

/*  ------------------------------------------------------------------------ */
/*                         Parameters from map:                              */
/*  ------------------------------------------------------------------------ */
/* section Parameters from map */
/*! \name  Parameters from map */
/* \{ */

/*! \brief return the mtz file that was use to generate the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
const char *mtz_hklin_for_map(int imol_map);

/*! \brief return the FP column in the file that was use to generate
  the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
const char *mtz_fp_for_map(int imol_map);

/*! \brief return the phases column in mtz file that was use to generate
  the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
const char *mtz_phi_for_map(int imol_map);

/*! \brief return the weight column in the mtz file that was use to
  generate the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say) or no weights were used. */
const char *mtz_weight_for_map(int imol_map);

/*! \brief return flag for whether weights were used that was use to
  generate the map

  return 0 when no weights were used or there is no mtz file
  associated with that map. */
short int mtz_use_weight_for_map(int imol_map);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the parameter that made the map, #f or something
  like ("xxx.mtz" "FPH" "PHWT" "" #f) */
SCM map_parameters_scm(int imol);
/*! \brief return the parameter that made the map, #f or something
  like (45 46 47 90 90 120), angles in degress */
SCM cell_scm(int imol);
/*! \brief return the parameter of the molecule, something 
  like (45 46 47 90 90 120), angles in degress */
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the parameter that made the map, False or something
  like ["xxx.mtz", "FPH", "PHWT", "", False] */
PyObject *map_parameters_py(int imol);
/*! \brief return the parameter that made the map, False or something
  like [45, 46, 47, 90, 90, 120], angles in degress */
PyObject *cell_py(int imol);
#endif /* USE_PYTHON */
#endif /* c++ */
/*! \} */


/*  ------------------------------------------------------------------------ */
/*                         Write PDB file:                                   */
/*  ------------------------------------------------------------------------ */
/* section PDB Functions */
/*! \name  PDB Functions */
/* \{ */

/*! \brief write molecule number imol as a PDB to file file_name */
/*  return 0 on success, 1 on error. */
int write_pdb_file(int imol, const char *file_name); 

/*! \brief write molecule number imol's residue range as a PDB to file
  file_name */
/*  return 0 on success, 1 on error. */
int write_residue_range_to_pdb_file(int imol, const char *chainid, 
				    int resno_start, int resno_end, 
				    const char *filename);

/*! \brief save all modified coordinates molecules to the default
  names and save the state too. */
int quick_save(); 
/*! \} */


/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */
/* section Refmac Functions */
/*! \name  Refmac Functions */
/* \{ */
void execute_refmac(GtkWidget *window); /* lookup stuff here. */
/*  this is the option menu callback - does nothing. */
void refmac_molecule_button_select(GtkWidget *item, GtkPositionType pos); 
int set_refmac_molecule(int imol); /* used by callback.c */
void fill_option_menu_with_refmac_options(GtkWidget *optionmenu);
void fill_option_menu_with_refmac_methods_options(GtkWidget *optionmenu);
void fill_option_menu_with_refmac_phase_input_options(GtkWidget *optionmenu);
void fill_option_menu_with_refmac_labels_options(GtkWidget *optionmenu);
void fill_option_menu_with_refmac_file_labels_options(GtkWidget *optionmenu);
void fill_option_menu_with_refmac_ncycle_options(GtkWidget *optionmenu);

void update_refmac_column_labels_frame(GtkWidget *optionmenu, 
				       GtkWidget *fobs_menu, GtkWidget *fiobs_menu, GtkWidget *fpm_menu,
				       GtkWidget *f_free_menu,
				       GtkWidget *phases_menu, GtkWidget *fom_menu, GtkWidget *hl_menu);


void free_memory_run_refmac(GtkWidget *window); 

/*! \brief set counter for runs of refmac so that this can be used to
  construct a unique filename for new output */
void set_refmac_counter(int imol, int refmac_count);
/*! \brief the name for refmac 

 @return a stub name used in the construction of filename for refmac output */
const char *refmac_name(int imol);

/* some methods to get refmac run parameters */
int get_refmac_refinement_method(void);
void set_refmac_refinement_method(int method);
int get_refmac_phase_input(void);
void set_refmac_phase_input(int phase_flag);
void set_refmac_use_tls(int state);
int refmac_use_tls_state(void);
void set_refmac_use_twin(int state);
int refmac_use_twin_state(void);
void set_refmac_use_sad(int state);
int refmac_use_sad_state(void);
int get_refmac_ncycles(void);
void set_refmac_ncycles(int no_cycles);
void add_refmac_ncycle_no(int cycle);
void set_refmac_use_ncs(int state);
int refmac_use_ncs_state(void);
void set_refmac_use_intensities(int state);
int refmac_use_intensities_state(void);
int refmac_imol_coords(void);
void add_refmac_sad_atom(const char *atom_name, float fp, float fpp, float lambda);
void add_refmac_sad_atom_fp(const char *atom_name, float fp, float fpp);
void add_refmac_sad_atom_lambda(const char *atom_name, float lambda);
void clear_refmac_sad_atoms();
void store_refmac_mtz_file_label(GtkWidget *label);
GtkWidget *get_refmac_mtz_file_label(void);
void fill_refmac_sad_atom_entry(GtkWidget *widget);
short int get_refmac_used_mtz_file_state();
void set_refmac_used_mtz_file(int state);
const gchar *get_saved_refmac_file_filename(void);
void set_stored_refmac_file_mtz_filename(int imol, const char *mtz_filename);
void save_refmac_params_to_map(int imol_map,
			       const char *mtz_filename,
			       const char *fobs_col,
			       const char *sigfobs_col,
			       const char *r_free_col,
			       int r_free_flag_sensible);
void save_refmac_phase_params_to_map(int imol_map,
			     	     const char *phi,
				     const char *fom,
				     const char *hla,
				     const char *hlb,
				     const char *hlc,
				     const char *hld);
#ifdef __cplusplus
#ifdef USE_GUILE
SCM get_refmac_sad_atom_info_scm();
#endif /* GUILE */
#ifdef USE_PYTHON
PyObject *get_refmac_sad_atom_info_py();
#endif /* PYTHON */
#endif /* c++ */

/*! \brief swap the colours of maps 

  swap the colour of maps imol1 and imol2.  Useful to some after
  running refmac, so that the map to be build into is always the same
  colour*/
void swap_map_colours(int imol1, int imol2);
/*! \brief flag to enable above

call this with istate=1 */
void set_keep_map_colour_after_refmac(int istate);

/*! \brief the keep-map-colour-after-refmac internal state

  @return 1 for "yes", 0 for "no"  */
int keep_map_colour_after_refmac_state();

/* refmac vresion testing, returns 1 for new refmac (>5.3) otherwise 0 */
int refmac_runs_with_nolabels(void);

/* \} */

/*  --------------------------------------------------------------------- */
/*                      symmetry                                          */
/*  --------------------------------------------------------------------- */
/* section Symmetry Functions */
/*! \name Symmetry Functions */
/* \{ */
char* get_text_for_symmetry_size_widget(); /* const gchar *text */
void set_symmetry_size_from_widget(const char *text); 
/*! \brief set the size of the displayed symmetry */
void set_symmetry_size(float f); 
double* get_symmetry_bonds_colour(int imol);
/*! \brief is symmetry master display control on? */
short int get_show_symmetry(); /* master */
/*! \brief set display symmetry, master controller */
void set_show_symmetry_master(short int state);
/*! \brief set display symmetry for molecule number mol_no

   pass with state=0 for off, state=1 for on */
void set_show_symmetry_molecule(int mol_no, short int state);
/*! \brief display symmetry as CAs?


   pass with state=0 for off, state=1 for on */
void symmetry_as_calphas(int mol_no, short int state); 
/*! \brief what is state of display CAs for molecule number mol_no?

   return state=0 for off, state=1 for on
*/
short int get_symmetry_as_calphas_state(int imol);

/*! \brief set the colour map rotation (i.e. the hue) for the symmetry
    atoms of molecule number imol  */
void set_symmetry_molecule_rotate_colour_map(int imol, int state);

/*! \brief should there be colour map rotation (i.e. the hue) change
    for the symmetry atoms of molecule number imol?

   return state=0 for off, state=1 for on
*/
int symmetry_molecule_rotate_colour_map_state(int imol);

void set_symmetry_colour_by_symop(int imol, int state);
void set_symmetry_whole_chain(int imol, int state);
void set_symmetry_atom_labels_expanded(int state);
GtkWidget *wrapped_create_show_symmetry_window();
void symmetry_colour_adjustment_changed (GtkAdjustment *adj, 
					 GtkWidget *window); 
GtkWidget *symmetry_molecule_controller_dialog();

/* used by destroy callback, needed because there should only be one of these. */
void set_symmetry_controller_dialog_widget(GtkWidget *w); 

/*! \brief molecule number imol has a unit cell? 

   @return 1 on "yes, it has a cell", 0 for "no" */
int has_unit_cell_state(int imol); 


/*! \brief Undo symmetry view. Translate back to main molecule from
  this symmetry position.  */
int undo_symmetry_view();

/*! \brief return the molecule number.

@return -1 if there is no molecule with symmetry displayed.  */
int first_molecule_with_symmetry_displayed();

/*! \brief save the symmetry coordinates of molecule number imol to
  filename 

Allow a shift of the coordinates to the origin before symmetry
expansion is apllied (this is how symmetry works in Coot
internals). */
void save_symmetry_coords(int imol,
			  const char *filename, 
			  int symop_no, 
			  int shift_a, 
			  int shift_b, 
			  int shift_c,
			  int pre_shift_to_origin_na,
			  int pre_shift_to_origin_nb,
			  int pre_shift_to_origin_nc);

/*! \brief create a new molecule (molecule number is the return value)
  from imol.  

The rotation/translation matrix components are given in *orthogonal*
coordinates.

Allow a shift of the coordinates to the origin before symmetry
expansion is aplied.  

Pass "" as the name-in and a name will be constructed for you.

Return -1 on failure. */ 
int new_molecule_by_symmetry(int imol,
			     const char *name,
			     double m11, double m12, double m13, 
			     double m21, double m22, double m23, 
			     double m31, double m32, double m33, 
			     double tx, double ty, double tz,
			     int pre_shift_to_origin_na,
			     int pre_shift_to_origin_nb,
			     int pre_shift_to_origin_nc);



/*! \brief create a new molecule (molecule number is the return value)
  from imol, but only for atom that match the
  mmdb_atom_selection_string.

The rotation/translation matrix components are given in *orthogonal*
coordinates.

Allow a shift of the coordinates to the origin before symmetry
expansion is aplied.  

Pass "" as the name-in and a name will be constructed for you.

Return -1 on failure. */ 
int new_molecule_by_symmetry_with_atom_selection(int imol, 
						 const char *name,
						 const char *mmdb_atom_selection_string,
						 double m11, double m12, double m13, 
						 double m21, double m22, double m23, 
						 double m31, double m32, double m33, 
						 double tx, double ty, double tz,
						 int pre_shift_to_origin_na,
						 int pre_shift_to_origin_nb,
						 int pre_shift_to_origin_nc);
			     

/*! \brief create a new molecule (molecule number is the return value)
  from imol.  
*/
int new_molecule_by_symop(int imol, const char *symop_string,
			  int pre_shift_to_origin_na,
			  int pre_shift_to_origin_nb,
			  int pre_shift_to_origin_nc);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the pre-shift (the shift that translates the centre
  of the molecule as close as possible to the origin) as a list of
  ints or scheme false on failure  */
SCM origin_pre_shift_scm(int imol);
#endif  /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the pre-shift (the shift that translates the centre
  of the molecule as close as possible to the origin) as a list of
  ints or Python false on failure  */
PyObject *origin_pre_shift_py(int imol);
#endif  /* USE_PYTHON */
#endif 

void setup_save_symmetry_coords();

void save_symmetry_coords_from_fileselection(GtkWidget *fileselection);

/*! \brief set the space group for a coordinates molecule

 for shelx FA pdb files, there is no space group.  So allow the user
   to set it.  This can be initted with a HM symbol or a symm list for
   clipper.  Return the succes status of the setting. */
short int set_space_group(int imol, const char *spg);

/*! \brief set the cell shift search size for symmetry searching.

When the coordinates for one (or some) symmetry operator are missing
(which happens sometimes, but rarely), try changing setting this to 2
(default is 1).  It slows symmetry searching, which is why it is not
set to 2 by default.  */
void set_symmetry_shift_search_size(int shift);

/* \} */ /* end of symmetry functions */

/*  ------------------------------------------------------------------- */
/*                    file selection                                    */
/*  ------------------------------------------------------------------- */
/* section File Selection Functions */
/*! \name File Selection Functions */
/* \{ */ /* start of file selection functions */

/* so that we can save/set the directory for future fileselections
   (i.e. the new fileselection will open in the directory that the
   last one ended in) */
void set_directory_for_fileselection(GtkWidget *coords_fileselection1); 
void save_directory_from_fileselection(const GtkWidget *fileselection);
void save_directory_for_saving_from_fileselection(const GtkWidget *fileselection);
void set_file_for_save_fileselection(GtkWidget *fileselection);

/* we include thes functions for the chooser here */
#ifdef COOT_USE_GTK2_INTERFACE
void set_directory_for_filechooser(GtkWidget *coords_fileselection1); 
void save_directory_from_filechooser(const GtkWidget *fileselection);
void save_directory_for_saving_from_filechooser(const GtkWidget *fileselection);
#endif

#ifdef __cplusplus
#ifdef USE_GUILE
/* Return the default file name suggestion (that would come up in the
   save coordinates dialog) or scheme false if imol is not a valid
   model molecule. */
SCM save_coords_name_suggestion_scm(int imol); 
#endif /*  USE_GUILE */
#ifdef USE_PYTHON
/* Return the default file name suggestion (that would come up in the
   save coordinates dialog) or Python false if imol is not a valid
   model molecule. */
PyObject *save_coords_name_suggestion_py(int imol); 
#endif /*  USE_PYTHON */
#endif /*  __cplusplus */

/* Eleanor likes to sort her files by date when selecting a file */

/* return the button. */
GtkWidget *add_sort_button_fileselection(GtkWidget *fileselection); 

void add_ccp4i_project_optionmenu(GtkWidget *fileselection, int file_selector_type);
void add_ccp4i_projects_to_optionmenu(GtkWidget *optionmenu, int file_selector_type, GtkSignalFunc func);
void add_ccp4i_project_shortcut(GtkWidget *fileselection);
void option_menu_refmac_ccp4i_project_signal_func(GtkWidget *item, GtkPositionType pos);
void run_refmac_ccp4i_option_menu_signal_func(GtkWidget *item, GtkPositionType pos);
void clear_refmac_ccp4i_project();
GtkWidget *lookup_file_selection_widgets(GtkWidget *item, int file_selector_type);

/* We wrote this button/callback by hand, most of the rest are in
   callbacks.c  */

#ifdef COOT_USE_GTK2_INTERFACE
void fileselection_sort_button_clicked( GtkWidget *sort_button,
					GtkWidget *file_list); 
#else
void fileselection_sort_button_clicked_gtk1( GtkWidget *sort_button,
					     GtkCList  *file_list); 
#endif

void push_the_buttons_on_fileselection(GtkWidget *filter_button, 
				       GtkWidget *sort_button,
				       GtkWidget *fileselection);

/* \} */ /* end of file selection functions */

/*  -------------------------------------------------------------------- */
/*                     history                                           */
/*  -------------------------------------------------------------------- */
/* section History Functions */
/*! \name  History Functions */

/* \{ */ /* end of file selection functions */
/* We don't want this exported to the scripting level interface,
   really... (that way lies madness, hehe). oh well... */

/*! \brief print the history in scheme format */
void print_all_history_in_scheme();
/*! \brief print the history in python format */
void print_all_history_in_python();

/*! \brief set a flag to show the text command equivalent of gui
  commands in the console as they happen. 

  1 for on, 0 for off. */
void set_console_display_commands_state(short int istate);
/*! \brief set a flag to show the text command equivalent of gui
  commands in the console as they happen in bold and colours.

  colour_flag: pass  1 for on, 0 for off.

  colour_index 0 to 7 inclusive for various different colourings.
 */
void set_console_display_commands_hilights(short int bold_flag, short int colour_flag, int colour_index);

/* \} */

/*  --------------------------------------------------------------------- */
/*                  state (a graphics_info thing)                         */
/*  --------------------------------------------------------------------- */
/* info */
/*! \name State Functions */
/* \{ */

/*! \brief save the current state to the default filename */
void save_state(); 

/*! \brief save the current state to file filename */
void save_state_file(const char *filename);

/*! \brief set the default state file name (default 0-coot.state.scm) */
/* set the filename */
void set_save_state_file_name(const char *filename);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief the save state file name 

  @return the save state file name*/
SCM save_state_file_name_scm();
#endif
#ifdef USE_PYTHON
/*! \brief the save state file name 

  @return the save state file name*/
PyObject *save_state_file_name_py();
#endif /* USE_PYTHON */
#endif	/* c++ */

/* only to be used in callbacks.c, don't export */
const char *save_state_file_name_raw();

/*! \brief set run state file status

0: never run it
1: ask to run it
2: run it, no questions */
void set_run_state_file_status(short int istat);
/*! \brief run the state file (reading from default filenname) */
void run_state_file();		/* just do it */
#ifdef USE_PYTHON
void run_state_file_py();		/* just do it */
#endif // USE_PYTHON
/*! \brief run the state file depending on the state variables */
void run_state_file_maybe();	/* depending on the above state variables */

GtkWidget *wrapped_create_run_state_file_dialog();
#ifdef USE_PYTHON
GtkWidget *wrapped_create_run_state_file_dialog_py();
#endif /* USE_PYTHON */

/* \} */

/*  -------------------------------------------------------------------- */
/*                     virtual trackball                                 */
/*  -------------------------------------------------------------------- */
/* subsection Virtual Trackball */
/*! \name The Virtual Trackball */
/* \{ */

#define VT_FLAT 1
#define VT_SPHERICAL 2
/*! \brief How should the mouse move the view?

mode=1 for "Flat", mode=2 for "Spherical Surface"  */
void vt_surface(int mode); 
/*! \brief return the mouse view status mode

mode=1 for "Flat", mode=2 for "Spherical Surface"  */
int  vt_surface_status();

/* \} */


/*  --------------------------------------------------------------------- */
/*                      clipping                                          */
/*  --------------------------------------------------------------------- */
/* section Clipping Functions */
/*! \name  Clipping Functions */
/* \{ */
void do_clipping1_activate();
void clipping_adjustment_changed (GtkAdjustment *adj, GtkWidget *window);

void set_clipping_back( float v);
void set_clipping_front(float v);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                         Unit Cell                                        */
/*  ----------------------------------------------------------------------- */
/* section Unit Cell interface */
/*! \name  Unit Cell interface */
/* \{ */

/*! \brief return the stage of show unit cell for molecule number imol */
short int get_show_unit_cell(int imol);

/*! \brief set the state of show unit cell for all molecules

1 for displayed   
0 for undisplayed */
void set_show_unit_cells_all(short int istate);

/*! \brief set the state of show unit cell for the particular molecule number imol 

1 for displayed   
0 for undisplayed */
void set_show_unit_cell(int imol, short int istate);

void set_unit_cell_colour(float red, float green, float blue);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                         Colour                                           */
/*  ----------------------------------------------------------------------- */
/* section Colour */
/*! \name  Colour */
/* \{ */

/* set the colour merge ratio (a fraction 0.0 to 1.0) */
void set_symmetry_colour_merge(float v);

/*! \brief set the hue change step on reading a new molecule */
void set_colour_map_rotation_on_read_pdb(float f); 

/*! \brief shall the hue change step be used?

 @param i 0 for no, 1 for yes */
void set_colour_map_rotation_on_read_pdb_flag(short int i); 

/*! \brief shall the colour map rotation apply only to C atoms?

 @param i 0 for no, 1 for yes */
void set_colour_map_rotation_on_read_pdb_c_only_flag(short int i);

/*! \brief colour molecule number imol by chain type */
void set_colour_by_chain(int imol); 

/*! \brief colour molecule number imol by molecule */
void set_colour_by_molecule(int imol); 

/* get the value of graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag */
int get_colour_map_rotation_on_read_pdb_c_only_flag();

/*! \brief set the symmetry colour base */
void set_symmetry_colour(float r, float g, float b); 

/* \} */

/*  Section Map colour*/
/*! \name   Map colour*/
/* \{ */
/*! \brief set the colour map rotation (hue change) for maps

   default: for maps is 14 degrees. */
void set_colour_map_rotation_for_map(float f); /* "global"/default */

/* widget work */
GtkWidget *wrapped_create_coords_colour_control_dialog();

/*! \brief set the colour map rotation for molecule number imol

theta is in degrees */
void  set_molecule_bonds_colour_map_rotation(int imol, float theta);

/*! \brief Get the colour map rotation for molecule number imol */
float get_molecule_bonds_colour_map_rotation(int imol); 
/* \} */

/*  ----------------------------------------------------------------------- */
/*                         Anisotropic Atoms */
/*  ----------------------------------------------------------------------- */
/* section Anisotropic Atoms Interface */
/*! \name  Anisotropic Atoms Interface */
/* \{ */
/*  we use the text interface to this in callback.c rather */
/*  than getting the float directly. */
 
float get_limit_aniso();           /* not a function of the molecule */

short int get_show_limit_aniso();  /* not a function of the molecule */

short int get_show_aniso();       /*  not a function of the molecule */

void set_limit_aniso(short int state);

void set_aniso_limit_size_from_widget(const char *text);

void set_show_aniso(int state); 

char *get_text_for_aniso_limit_radius_entry();

void set_aniso_probability(float f);

float get_aniso_probability();

/* \} */

/*  ---------------------------------------------------------------------- */
/*                         Display Functions                               */
/*  ---------------------------------------------------------------------- */
/* section Display Functions */
/*! \name  Display Functions */
/* \{ */
/*  currently doesn't get seen when the window starts due to */
/*  out-of-order issues. */

/*! \brief set the window size */
void   set_graphics_window_size(int x_size, int y_size);
/*! \brief set the window position */
void   set_graphics_window_position(int x_pos, int y_pos);
void store_graphics_window_position(int x_pos, int y_pos); /*  "configure_event" callback */

/* a general purpose version of the above, where we pass a widget flag */
void store_window_position(int window_type, GtkWidget *w);
void store_window_size(int window_type, GtkWidget *w);

/*! \brief draw a frame */
void graphics_draw(); 	/* and wrapper interface to gtk_widget_draw(glarea)  */

/*! \brief try to turn on Zalman stereo mode  */
void zalman_stereo_mode();
/*! \brief try to turn on stereo mode  */
void hardware_stereo_mode();
/*! \brief what is the stero state?

  @return 1 for in hardware stereo, 2 for side by side stereo, else return 0. */
int  stereo_mode_state();
/*! \brief try to turn on mono mode  */
void mono_mode();
/*! \brief turn on side bye side stereo mode */
void side_by_side_stereo_mode(short int use_wall_eye_mode);

/* DTI stereo mode - undocumented, secret interface for testing, currently. 
state should be 0 or 1. */
/* when it works, call it dti_side_by_side_stereo_mode() */
void set_dti_stereo_mode(short int state);

/*! \brief how much should the eyes be separated in stereo mode? 

   @param f the angular difference (in multiples of 4.5 degrees) */
void set_hardware_stereo_angle_factor(float f);
/*! \brief return the hardware stereo angle factor */
float hardware_stereo_angle_factor_state();

/*! \brief set position of Model/Fit/Refine dialog */
void set_model_fit_refine_dialog_position(int x_pos, int y_pos);
/*! \brief set position of Display Control dialog */
void set_display_control_dialog_position(int x_pos, int y_pos);
/*! \brief set position of Go To Atom dialog */
void set_go_to_atom_window_position(int x_pos, int y_pos);
/*! \brief set position of Delete dialog */
void set_delete_dialog_position(int x_pos, int y_pos);
/*! \brief set position of the Rotate/Translate Residue Range dialog */
void set_rotate_translate_dialog_position(int x_pos, int y_pos);
/*! \brief set position of the Accept/Reject dialog */
void set_accept_reject_dialog_position(int x_pos, int y_pos);
/*! \brief set position of the Ramachadran Plot dialog */
void set_ramachandran_plot_dialog_position(int x_pos, int y_pos);
/*! \brief set edit chi angles dialog position */
void set_edit_chi_angles_dialog_position(int x_pos, int y_pos);

/* \} */

/*  ---------------------------------------------------------------------- */
/*                         Smooth "Scrolling" */
/*  ---------------------------------------------------------------------- */
/* section Smooth Scrolling */
/*! \name  Smooth Scrolling */
/* \{ */

/*! \brief set smooth scrolling 

  @param v use v=1 to turn on smooth scrolling, v=0 for off (default on). */
void set_smooth_scroll_flag(int v);

/*! \brief return the smooth scrolling state */
int  get_smooth_scroll();

void set_smooth_scroll_steps_str(const char * t);

/*  useful exported interface */
/*! \brief set the number of steps in the smooth scroll

   Set more steps (e.g. 50) for more smoothness (default 10).*/
void set_smooth_scroll_steps(int i); 

char  *get_text_for_smooth_scroll_steps();

void  set_smooth_scroll_limit_str(const char *t);

/*  useful exported interface */
/*! \brief do not scroll for distances greater this limit */
void  set_smooth_scroll_limit(float lim); 

char *get_text_for_smooth_scroll_limit();

/* \} */


/*  ---------------------------------------------------------------------- */
/*                         Font Size */
/*  ---------------------------------------------------------------------- */
/* section Font Size */
/*! \name  Font Size */
/* \{ */

/*! \brief set the font size

  @param i 1 (small) 2 (medium, default) 3 (large) */
void set_font_size(int i);

/*! \brief return the font size

  @return 1 (small) 2 (medium, default) 3 (large) */
int get_font_size();

/*! \brief set the colour of the atom label font - the arguments are
  in the range 0->1 */
void set_font_colour(float red, float green, float blue);
/* \} */

/*  ---------------------------------------------------------------------- */
/*                         Rotation Centre                                 */
/*  ---------------------------------------------------------------------- */
/* section Rotation Centre */
/*! \name  Rotation Centre */
/* \{ */
void set_rotation_centre_size_from_widget(const gchar *text); /* and redraw */
void set_rotation_centre_size(float f); /* and redraw (maybe) */
gchar *get_text_for_rotation_centre_cube_size(); 
short int recentre_on_read_pdb(); 
void set_recentre_on_read_pdb(short int);
void set_rotation_centre(float x, float y, float z);
// The redraw happens somewhere else...
void set_rotation_centre_internal(float x, float y, float z); 
float rotation_centre_position(int axis); /* only return one value: x=0, y=1, z=2 */
/* \} */

/*  ---------------------------------------------------------------------- */
/*                         orthogonal axes                                 */
/*  ---------------------------------------------------------------------- */

/* section Orthogonal Axes */
/*! \name Orthogonal Axes */
/* \{ */
/* Draw the axes in the top left? 
  
0 off, 1 on */
void set_draw_axes(int i);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  utility function                                        */
/*  ----------------------------------------------------------------------- */
/* section Atom Selection Utilities */
/*! \name  Atom Selection Utilities */
/* \{ */

/* does not account for alternative conformations properly */
/* return -1 if atom not found. */
int atom_index(int imol, const char *chain_id, int iresno, const char *atom_id); 
/* Refine zone needs to be passed atom indexes (which it then converts */
/* to residue numbers - sigh).  So we need a function to get an
   atom. Return -1 on failure */
/* index from a given residue to use with refine_zone().  Return -1 on failure */
int atom_index_first_atom_in_residue(int imol, const char *chain_id, 
				     int iresno, const char *ins_code);
/* For rotamers, we are given a residue spec (and altconf), we need
  the index of the first atom of this type, no atom name is given,
  hence we cannot use full_atom_spec_to_atom_index(). */
int atom_index_first_atom_in_residue_with_altconf(int imol, 
						  const char *chain_id, 
						  int iresno, 
						  const char *ins_code, 
						  const char *alt_conf);
float median_temperature_factor(int imol);
float average_temperature_factor(int imol);
float standard_deviation_temperature_factor(int imol);

void clear_pending_picks(); 
char *centre_of_mass_string(int imol);
#ifdef USE_PYTHON
char *centre_of_mass_string_py(int imol);
#endif // PYTHON
/*! \brief set the default temperature factor for newly created atoms
  (initial default 20) */
void set_default_temperature_factor_for_new_atoms(float new_b);
/*! \brief return the default temperature factor for newly created atoms */
float default_new_atoms_b_factor();

/*! \brief reset temperature factor for all moved atoms to the default
  for new atoms (usually 30) */
void set_reset_b_factor_moved_atoms(int state);
/*! \brief return the state if temperature factors shoudl be reset for
  moved atoms */
int get_reset_b_factor_moved_atoms_state();

/*! \brief set a numberical attibute to the atom with the given specifier.

Attributes can be "x", "y","z", "B", "occ" and the attribute val is a floating point number*/
int set_atom_attribute(int imol, const char *chain_id, int resno, const char *ins_code, const char *atom_name, const char*alt_conf, const char *attribute_name, float val);

/*! \brief set a string attibute to the atom with the given specifier.

Attributes can be "atom-name", "alt-conf", "element" or "segid". */
int set_atom_string_attribute(int imol, const char *chain_id, int resno, const char *ins_code, const char *atom_name, const char*alt_conf, const char *attribute_name, const char *attribute_value);

/*! \brief set lots of atom attributes at once by-passing the rebonding and redrawing of the above 2 functions  */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
int set_atom_attributes(SCM attribute_expression_list);
#endif 

#ifdef USE_PYTHON
int set_atom_attributes_py(PyObject *attribute_expression_list);
#endif 
#endif // __cplusplus

void set_residue_name(int imol, const char *chain_id, int res_no, const char *ins_code, const char *new_residue_name);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                            skeletonization                               */
/*  ----------------------------------------------------------------------- */
/* section Skeletonization Interface */
/*! \name  Skeletonization Interface */
/* \{ */
void skel_greer_on(); 
void skel_greer_off(); 

void skel_foadi_on(); 
void skel_foadi_off(); 

void skeletonize_map_by_optionmenu(GtkWidget *optionmenu);
void skeletonize_map_single_map_maybe(GtkWidget *window, int imol); 

GtkWidget *wrapped_create_skeleton_dialog();

/*! \brief skeletonize molecule number imol

   the prune_flag should almost  always be 0.  */
int skeletonize_map(int prune_flag, int imol);

/*! \brief undisplay the skeleton on molecule number imol */
int unskeletonize_map(int imol); 

void set_initial_map_for_skeletonize(); /* set graphics_info variable
					   for use in callbacks.c */

/*! \brief set the skeleton search depth, used in baton building

  For high resolution maps, you need to search deeper down the skeleton tree.  This 
  limit needs to be increased to 20 or so for high res maps (it is 10 by default)  */
void set_max_skeleton_search_depth(int v); /* for high resolution
					      maps, change to 20 or
					      something (default 10). */

/* set the radio buttons in the frame to the be on or off for the map
   that is displayed in the optionmenu (those menu items "active"
   callbacks (graphics_info::skeleton_map_select change
   g.map_for_skeletonize).  */
void set_on_off_skeleton_radio_buttons(GtkWidget *skeleton_frame);
void set_on_off_single_map_skeleton_radio_buttons(GtkWidget *skeleton_frame, int imol);


/*  ----------------------------------------------------------------------- */
/*                  skeletonization level widgets                           */
/*  ----------------------------------------------------------------------- */

gchar *get_text_for_skeletonization_level_entry(); 

void set_skeletonization_level_from_widget(const char *txt); 

gchar *get_text_for_skeleton_box_size_entry(); 

void set_skeleton_box_size_from_widget(const char *txt); 

/*! \brief the box size (in Angstroms) for which the skeleton is displayed */
void set_skeleton_box_size(float f);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                        Skeleton                                          */
/*  ----------------------------------------------------------------------- */
/* section Skeleton Colour */
/*! \name  Skeleton Colour */
/* \{ */
void handle_skeleton_colour_change(int mol, gdouble* map_col);
void set_skeleton_colour(int imol, float r, float g, float b);

gdouble* get_skeleton_colour(); 

/* \} */

/*  ----------------------------------------------------------------------- */
/*                         read a ccp4 map                                  */
/*  ----------------------------------------------------------------------- */
/* section Read Maps */
/*! \name  Read Maps */
/* \{ */

/*! \brief read a CCP4 map or a CNS map (despite the name). */
int handle_read_ccp4_map(const char* filename, int is_diff_map_flag); 
/* \} */

/*  ----------------------------------------------------------------------- */
/*                        save coordintes                                   */
/*  ----------------------------------------------------------------------- */
/* section Save Coordinates */
/*! \name  Save Coordinates */
/* \{ */


void save_coordinates_using_widget(GtkWidget *widget); /* do a get_user_data for
					     the molecule and a lookup
					     of the entry? to find the
					     filename in c-interface,
					     not in the callback.c  */
/*! \brief save coordinates of molecule number imol in filename

  @return status 1 is good (success), 0 is fail. */
int save_coordinates(int imol, const char *filename);

void set_save_coordinates_in_original_directory(int i);

/* not really a button select, its a menu item select */
/* not productive */
void save_molecule_coords_button_select(GtkWidget *item, GtkPositionType pos); 

/* access to graphics_info_t::save_imol for use in callback.c */
int save_molecule_number_from_option_menu();
/* access from callback.c, not to be used in scripting, I suggest.
   Sets the *save* molecule number */
void set_save_molecule_number(int imol);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                        .phs file reading                                 */
/*  ----------------------------------------------------------------------- */
/* section Read Phases File Functions */
/*! \name  Read Phases File Functions */
/* \{ */

/*! \brief read phs file use coords to get cell and symm to make map 

uses pending data to make the map.

*/
void
read_phs_and_coords_and_make_map(const char *pdb_filename);

/*! \brief read a phs file, the cell and symm information is from
  previously read (most recently read) coordinates file

 For use with phs data filename provided on the command line */
int 
read_phs_and_make_map_using_cell_symm_from_previous_mol(const char *phs_filename);


/*! \brief read phs file and use a previously read molecule to provide
  the cell and symmetry information

@return the new molecule number, return -1 if problem creating the map
(e.g. not phs data, file not found etc).  */
int
read_phs_and_make_map_using_cell_symm_from_mol(const char *phs_filename, int imol);

int
read_phs_and_make_map_using_cell_symm_from_mol_using_implicit_phs_filename(int imol);

/*! \brief read phs file use coords to use cell and symm to make map */
int
read_phs_and_make_map_using_cell_symm(const char *phs_file_name,
				      const char *hm_spacegroup, float a, float b, float c,
				      float alpha, float beta, float gamma); /*!< in degrees */

/*! \brief read a phs file and use the cell and symm in molecule
  number imol and use the resolution limits reso_lim_high (in Angstroems).

@param imol is the molecule number of the reference (coordinates)
molecule from which the cell and symmetry can be obtained.

@param phs_file_name is the name of the phs data file.

@param reso_lim_high is the high resolution limit in Angstroems.

@param reso_lim_low the low resoluion limit (currently ignored).  */
int 
read_phs_and_make_map_with_reso_limits(int imol, const char* phs_file_name,
				       float reso_lim_low, float reso_lim_high); 

/* work out the spacegroup from the given symm operators, e.g. return "P 1 21 1" 
given "x,y,z ; -x,y+1/2,-z" */
/* char * */
/* spacegroup_from_operators(const char *symm_operators_in_clipper_format);  */

void
graphics_store_phs_filename(const gchar *phs_filename); 

const char* graphics_get_phs_filename();

short int possible_cell_symm_for_phs_file(); 

gchar *get_text_for_phs_cell_chooser(int imol, char *field); 

/* \} */

/*  ----------------------------------------------------------------------- */
/*                                  Movement                                */
/*  ----------------------------------------------------------------------- */
/* section Graphics Move */
/*! \name Graphics Move */
/* \{ */
/*! \brief undo last move  */
void undo_last_move(); // suggested by Frank von Delft

/*! \brief translate molecule number imol by (x,y,z) in Angstroms  */
void translate_molecule_by(int imol, float x, float y, float z);

/*! \brief transform molecule number imol by the given rotation
  matrix, then translate by (x,y,z) in Angstroms  */
void transform_molecule_by(int imol, 
			   float m11, float m12, float m13,
			   float m21, float m22, float m23,
			   float m31, float m32, float m33,
			   float x, float y, float z);

void transform_zone(int imol, const char *chain_id, int resno_start, int resno_end, const char *ins_code,
		    float m11, float m12, float m13,
		    float m21, float m22, float m23,
		    float m31, float m32, float m33,
		    float x, float y, float z);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                        go to atom widget                                 */
/*  ----------------------------------------------------------------------- */
/* section Go To Atom Widget Functions */
/*! \name Go To Atom Widget Functions */
/* \{ */
GtkWidget *wrapped_create_goto_atom_window();

/*! \brief Post the Go To Atom Window */
void post_go_to_atom_window();
void fill_go_to_atom_window(GtkWidget *widget);

/* gchar *get_text_for_go_to_atom_chain_entry();  */
/* gchar *get_text_for_go_to_atom_residue_entry(); */
/* gchar *get_text_for_go_to_atom_atom_name_entry(); */

int go_to_atom_molecule_number();
char *go_to_atom_chain_id();
char *go_to_atom_atom_name();
int go_to_atom_residue_number();
char *go_to_atom_ins_code();
char *go_to_atom_alt_conf();


/*! \brief set the go to atom specification

   It seems important for swig that the char * arguments are const
   char *, not const gchar * (or else we get wrong type of argument
   error on (say) "A"

@return the success status of the go to.  0 for fail, 1 for success.
*/
int set_go_to_atom_chain_residue_atom_name(const char *t1_chain_id, int iresno, 
					   const char *t3_atom_name);
int set_go_to_atom_chain_residue_atom_name_no_redraw(const char *t1, int iresno, const char *t3, 
						     short int make_the_move_flag);

/* FIXME one day */
/* #ifdef __cplusplus */
/* int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec); */
/* #endif // __cplusplus */

int set_go_to_atom_chain_residue_atom_name_strings(const gchar *t1, 
						   const gchar *t2, 
						   const gchar *txt); 


int goto_next_atom_maybe_new(GtkWidget *window);
int goto_previous_atom_maybe_new(GtkWidget *window);

/*! \brief update the Go To Atom widget entries to atom closest to
  screen centre. */
void update_go_to_atom_from_current_position(); 

/* used by keypress (return) callbacks */

/*  read the widget values and apply them to the graphics */
 
int apply_go_to_atom_values(GtkWidget * window);


/* moving gtk function out of build functions, delete_atom() updates
   the go to atom atom list on deleting an atom  */
void update_go_to_atom_residue_list(int imol);

/*  return an atom index */
/*! \brief what is the atom index of the given atom? */
int atom_spec_to_atom_index(int mol, char *chain, int resno, char *atom_name); 

/*! \brief what is the atom index of the given atom? */
int full_atom_spec_to_atom_index(int imol, const char *chain, int resno,
				 const char *inscode, const char *atom_name,
				 const char *altloc); 

/*! \brief update the Go To Atom window */
void update_go_to_atom_window_on_changed_mol(int imol);

/*! \brief update the Go To Atom window.  This updates the option menu
  for the molecules. */
void update_go_to_atom_window_on_new_mol();

void update_go_to_atom_window_on_other_molecule_chosen(int imol);

/*! \brief set the molecule for the Go To Atom

   For dynarama callback sake. The widget/class knows which
   molecule that it was generated from, so in order to go to the
   molecule from dynarama, we first need to the the molecule - because
   set_go_to_atom_chain_residue_atom_name() does not mention the
   molecule (see "Next/Previous Residue" for reasons for that).  This
   function simply calls the graphics_info_t function of the same
   name.

   Also used in scripting, where go-to-atom-chain-residue-atom-name
   does not mention the molecule number. 

   20090914-PE set-go-to-atom-molecule can be used in a script and it
   should change the go-to-atom-molecule in the Go To Atom dialog (if
   it is being displayed).  This does mean, of course that using the
   ramachandran plot to centre on atoms will change the Go To Atom
   dialog.  Maybe that is surprising (maybe not).

*/
void set_go_to_atom_molecule(int imol); 

int go_to_atom_molecule_optionmenu_active_molecule(GtkWidget *widget); 


void save_go_to_atom_widget(GtkWidget *widget); /* store in a static */
void unset_go_to_atom_widget(); /* unstore the static go_to_atom_window */


void clear_atom_list(GtkWidget *atom_gtklist);

void apply_go_to_atom_from_widget(GtkWidget *widget);

void
on_go_to_atom_residue_list_select_child (GtkList         *list,
					 GtkWidget       *widget,
					 gpointer         user_data);

#ifdef COOT_USE_GTK2_INTERFACE
/* stuff moved into graphics_info_t */
#else
void on_go_to_atom_residue_tree_selection_changed_gtk1 (GtkList         *gtktree,
							gpointer         user_data);
#endif

/* Nothing calls this? */
/* void */
/* on_go_to_atom_residue_list_selection_changed (GtkList         *list, */
/* 						   gpointer         user_data); */


#ifdef COOT_USE_GTK2_INTERFACE

#else
void on_go_to_atom_atom_list_selection_changed_gtk1 (GtkList         *list,
						gpointer         user_data);
void
on_go_to_atom_residue_list_unselect_child (GtkList         *list,
					   GtkWidget       *widget,
					   gpointer         user_data);
#endif
/* \} */


/*  ----------------------------------------------------------------------- */
/*                  autobuilding control                                    */
/*  ----------------------------------------------------------------------- */
/* section AutoBuilding functions (Defunct) */
/* void autobuild_ca_on();  - moved to junk */

void autobuild_ca_off(); 

void test_fragment(); 

void do_skeleton_prune();

int test_function(int i, int j);


/*  ----------------------------------------------------------------------- */
/*                  map and molecule control                                */
/*  ----------------------------------------------------------------------- */
/* section Map and Molecule Control */
/*! \name Map and Molecule Control */
/* \{ */
void save_display_control_widget_in_graphics(GtkWidget *widget); 

GtkWidget *wrapped_create_display_control_window();

/*! \brief display the Display Constrol window  */
void post_display_control_window(); 

void add_map_display_control_widgets(); 
void add_mol_display_control_widgets(); 
void add_map_and_mol_display_control_widgets(); 

void reset_graphics_display_control_window(); 
void close_graphics_display_control_window(); /* destroy widget */

/*! \brief make the map displayed/undisplayed, 0 for off, 1 for on */
void set_map_displayed(int imol, int state);
/*! \brief make the coordinates molecule displayed/undisplayed, 0 for off, 1 for on */
void set_mol_displayed(int imol, int state);
/*! \brief make the coordinates molecule active/inactve (clickable), 0
  for off, 1 for on */
void set_mol_active(int imol, int state);


/*! \brief return the display state of molecule number imol 

 @return 1 for on, 0 for off
*/
int mol_is_displayed(int imol); 
/*! \brief return the active state of molecule number imol
 @return 1 for on, 0 for off */
int mol_is_active(int imol); 
/*! \brief return the display state of molecule number imol
 @return 1 for on, 0 for off */
int map_is_displayed(int imol); 

/*! \brief if on_or_off is 0 turn off all maps displayed, for other
  values of on_or_off turn on all maps */
void set_all_maps_displayed(int on_or_off);

/*! \brief if on_or_off is 0 turn off all models displayed and active,
  for other values of on_or_off turn on all models. */
void set_all_models_displayed_and_active(int on_or_off);

/*! \brief return the spacegroup of molecule number imol 

@return "No Spacegroup" when the spacegroup of a molecule has not been
set.*/
char *show_spacegroup(int imol);


#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief return a list of symmetry operators as strings - or scheme false if
  that is not possible. */
SCM symmetry_operators_scm(int imol);
/* take the return value from above and return a xHM symbol (for
   testing currently) */
SCM symmetry_operators_to_xHM_scm(SCM symmetry_operators);
#endif 

#ifdef USE_PYTHON
/*! \brief return a list of symmetry operators as strings - or Python False if
  that is not possible. */
PyObject *symmetry_operators_py(int imol);
/* take the return value from above and return a xHM symbol (for
   testing currently) */
PyObject *symmetry_operators_to_xHM_py(PyObject *symmetry_operators);
#endif /* USE_PYTHON */
#endif /* c++ */


/* \} */

/*  ----------------------------------------------------------------------- */
/*                         Merge Molecules                                  */
/*  ----------------------------------------------------------------------- */
/* section Merge Molecules */
GtkWidget *wrapped_create_merge_molecules_dialog();
void do_merge_molecules_gui();
void do_merge_molecules(GtkWidget *dialog);
void fill_vbox_with_coordinates_options(GtkWidget *vbox,
					GtkSignalFunc checkbox_callback_func);
void merge_molecules_menu_item_activate(GtkWidget *item, 
					GtkPositionType pos);
void on_merge_molecules_check_button_toggled (GtkToggleButton *togglebutton,
					      gpointer         user_data);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM merge_molecules(SCM add_molecules, int imol);
#endif

#ifdef USE_PYTHON
PyObject *merge_molecules_py(PyObject *add_molecules, int imol);
#endif /* PYTHON */
#endif	/* c++ */


/*  ----------------------------------------------------------------------- */
/*                         Mutate Sequence and Loops GUI                    */
/*  ----------------------------------------------------------------------- */
/* section Mutate Sequence and Loops GUI */
GtkWidget *wrapped_create_mutate_sequence_dialog();
void do_mutate_sequence(GtkWidget *dialog); 
void mutate_sequence_molecule_menu_item_activate(GtkWidget *item, 
						 GtkPositionType pos);
/* void fill_chain_option_menu(GtkWidget *chain_option_menu, int imol); */
/* the generic form of the above - also used by superpose chain optionmenu */
/* void fill_chain_option_menu_with_callback(GtkWidget *chain_option_menu, 
					  int imol,
					  GtkSignalFunc callback); */
void mutate_sequence_chain_option_menu_item_activate (GtkWidget *item,
						      GtkPositionType pos);

GtkWidget *wrapped_fit_loop_dialog();
void fit_loop_from_widget(GtkWidget *w);

/*  ----------------------------------------------------------------------- */
/*                         Align and Mutate GUI                             */
/*  ----------------------------------------------------------------------- */
/* section Align and Mutate */
/*! \name  Align and Mutate */
/* \{ */
GtkWidget *wrapped_create_align_and_mutate_dialog();
/* return the handled_state, so that we know if we should kill the dialog or not */
int do_align_mutate_sequence(GtkWidget *w);
void align_and_mutate_molecule_menu_item_activate(GtkWidget *item, 
						  GtkPositionType pos);
void align_and_mutate_chain_option_menu_item_activate (GtkWidget *item,
						       GtkPositionType pos);

/*! \brief aligand and mutate the given chain to the given sequence  */
void align_and_mutate(int imol, const char *chain_id, const char *fasta_maybe, short int renumber_residues_flag);
/*! \brief set the penalty for affine gap and space when aligning */
void set_alignment_gap_and_space_penalty(float wgap, float wspace);


/* What are these functions?  consider deleting them - we have alignment_mismatches_* . */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM alignment_results_scm(int imol, const char* chain_id, const char *seq); 
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *alignment_results_py(int imol, const char* chain_id, const char *seq); 
#endif /* USE_PYTHON */
#endif  /* c++ */

/* \} */

/*  ----------------------------------------------------------------------- */
/*                         Renumber residue range                           */
/*  ----------------------------------------------------------------------- */
/* section Renumber Residue Range */
/*! \name Renumber Residue Range */

/* \{ */
/*! \brief renumber the given residue range by offset residues */
int renumber_residue_range(int imol, const char *chain_id, 
			   int start_res, int last_res, int offset);

GtkWidget *wrapped_create_renumber_residue_range_dialog();
void renumber_residues_from_widget(GtkWidget *window);

/*! \brief change chain id, residue number or insertion code for given
  residue  */
int change_residue_number(int imol, const char *chain_id, int current_resno, const char *current_inscode, int new_resno, const char *new_inscode);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                         Change chain id                                  */
/*  ----------------------------------------------------------------------- */
/* section Change Chain ID */
GtkWidget *wrapped_create_change_chain_id_dialog();
void change_chain_id_by_widget(GtkWidget *w);
void change_chain_ids_mol_option_menu_item_activate(GtkWidget *item,
						    GtkPositionType pos);
void change_chain_ids_chain_menu_item_activate(GtkWidget *item,
					       GtkPositionType pos);

void  change_chain_id(int imol, const char *from_chain_id, const char *to_chain_id, 
		      short int use_res_range_flag, int from_resno, int to_resno);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/* Paul fill me in please */
SCM change_chain_id_with_result_scm(int imol, const char *from_chain_id, const char *to_chain_id,
                                         short int use_res_range_flag, int from_resno, int to_resno);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *change_chain_id_with_result_py(int imol, const char *from_chain_id, const char *to_chain_id, short int use_res_range_flag, int from_resno, int to_resno);
#endif /* USE_PYTHON */
#endif  /* c++ */

/*  ----------------------------------------------------------------------- */
/*                  scripting                                               */
/*  ----------------------------------------------------------------------- */
/* section Scripting Interface */

/*! \name Scripting Interface */
/* \{ */

/*! \brief Can we run probe (was the executable variable set
  properly?) (predicate).

@return 1 for yes, 2 for no  */
int probe_available_p();
#ifdef USE_PYTHON
int probe_available_p_py();
#endif

/*! \brief do nothing - compatibility function */
void post_scripting_window(); 

/*! \brief pop-up a scripting window for scheming */
void post_scheme_scripting_window(); 

/*! \brief pop-up a scripting window for pythoning */
void post_python_scripting_window(); 

/* called from c-inner-main */
void run_command_line_scripts();

void setup_guile_window_entry(GtkWidget *entry); 
void setup_python_window_entry(GtkWidget *entry); 

/*  Check if this is needed still, I think not. */
#ifdef USE_GUILE
void guile_window_enter_callback( GtkWidget *widget,
				  GtkWidget *entry ); 
#endif /* USE_GUILE */

#ifdef USE_PYTHON
void python_window_enter_callback( GtkWidget *widget,
				   GtkWidget *entry ); 
#endif /* USE_PYTHON */

void set_guile_gui_loaded_flag(); 
void set_python_gui_loaded_flag(); 
void set_found_coot_gui(); 
void set_found_coot_python_gui(); 

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Monomer                                                 */
/*  ----------------------------------------------------------------------- */
/* section monomers */

/*! \name Monomer */
/* \{ */
/* Accession code */
void handle_get_accession_code(GtkWidget *widget); 

/* in here we check if libcheck is available (if scripting is available) */
GtkWidget *wrapped_create_libcheck_monomer_dialog();

/* Libcheck monomer code */
void handle_get_libcheck_monomer_code(GtkWidget *widget); 

/*! \brief import libcheck monomer give the 3-letter code. 

@return the new molecule number, if not -1 (error). */
int get_monomer(const char *three_letter_code);

/* Use the protein geometry dictionary to retrieve a set of
   coordinates quickly from cif data read in from the RCSB's Chemical
   Component Library.  There are no restraints from this method
   though. */
int get_monomer_from_dictionary(const char *three_letter_code, int idealised_flag);


int
handle_make_monomer_search(const char *text, GtkWidget *viewport);



/*  Don't let this be seen by standard c, since I am using a std::string */
/*  and now we make it return a value, which we can decode in the calling
    function. I make a dummy version for when GUILE is not being used in case
    there are functions in the rest of the code that call safe_scheme_command
    without checking if there is USE_GUILE first.
    importing ability for python modules, use_namespace to maintain the 
    namespace of the module, returns 1 if not running, 0 on success, -1 
    when error importing (no further information) */

/*! \brief run script file */
void run_script       (const char *filename);
void run_guile_script (const char *filename);
void run_python_script(const char *filename);
int import_python_module(const char *module_name, int use_namespace);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  regularize/refine                                       */
/*  ----------------------------------------------------------------------- */
/* section Regularization and Refinement */
/*! \name  Regularization and Refinement */
/* \{ */

void do_regularize(short int state); /* pass 0 for off (unclick togglebutton) */
void do_refine(short int state);

/* and a gui function used by those functions */
void do_regularize_kill_delete_dialog();

/*! \brief add a restraint on peptides to make them planar 

  This adds a 5 atom restraint that includes both CA atoms of the
  peptide.  Use this rather than editting the mon_lib_list.cif file. */
void add_planar_peptide_restraints(); 

/*! \brief remove restraints on peptides to make them planar. */
void remove_planar_peptide_restraints();

/* return 1 if planar peptide restraints are on, 0 if off */
int planar_peptide_restraints_state();

/*! \brief add restraints on the omega angle of the peptides

  (that is the torsion round the peptide bond).  Omega angles that are
  closer to 0 than to 180 will be refined as cis peptides (and of
  course if omega is greater than 90 then the peptide will be refined
  as a trans peptide (this is the normal case). */
void add_omega_torsion_restriants(); 

/*! \brief remove omega restraints on CIS and TRANS linked residues. */
void remove_omega_torsion_restriants(); 

/*! \brief set immediate replacement mode for refinement and
  regularization.  You need this (call with istate=1) if you are
  scripting refinement/regularization  */
void set_refinement_immediate_replacement(int istate);

/*! \brief query the state of the immediate replacement mode */
int  refinement_immediate_replacement_state();

/*! \brief set the number of frames for which the selected residue
  range flashes 

 On fast computers, this can be set to higher than the default for
 more aesthetic appeal. */
void set_residue_selection_flash_frames_number(int i);

/*! \brief accept the new positions of the regularized or refined residues 

    If you are scripting refinement and/or regularization, this is the
    function that you need to call after refine-zone or regularize-zone.  */
void accept_regularizement();
void clear_up_moving_atoms();	/* remove the molecule and bonds */
void clear_moving_atoms_object(); /* just get rid of just the bonds (redraw done here). */
/* now we use */
void fill_option_menu_with_refine_options(GtkWidget *option_menu);

GtkWidget *wrapped_create_refine_params_dialog(); 

void do_torsions_toggle(GtkWidget *button);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM refine_residues_scm(int imol, SCM r); /* presumes the alt_conf is "". */
SCM refine_residues_with_alt_conf_scm(int imol, SCM r, const char *alt_conf); /* to be renamed later. */
SCM regularize_residues_scm(int imol, SCM r); /* presumes the alt_conf is "". */
SCM regularize_residues_with_alt_conf_scm(int imol, SCM r, const char *alt_conf); 
#endif
#ifdef USE_PYTHON
PyObject *refine_residues_py(int imol, PyObject *r);  /* to be renamed later. */
PyObject *refine_residues_with_alt_conf_py(int imol, PyObject *r, const char *alt_conf);  /* to be renamed later. */
#endif /* PYTHON */
#endif /* c++ */

/*! \brief turn on (or off) torsion restraints 

   Pass with istate=1 for on, istate=0 for off.
*/
void set_refine_with_torsion_restraints(int istate); 
/*! \brief return the state of above */
int refine_with_torsion_restraints_state(); 
void set_refine_params_toggle_buttons(GtkWidget *button);


/*! \brief set the relative weight of the geometric terms to the map terms 

 The default is 60.

 The higher the number the more weight that is given to the map terms
 but the resulting chi squared values are higher).  This will be
 needed for maps generated from data not on (or close to) the absolute
 scale or maps that have been scaled (for example so that the sigma
 level has been scaled to 1.0).

*/
void set_matrix(float f);

/*! \brief return the relative weight of the geometric terms to the map terms. */
float matrix_state();

/* Now the refinement weight can be set from an entry in the refine_params_dialog. */
void set_refinemenent_weight_from_entry(GtkWidget *entry);


/*! \brief change the +/- step for autoranging (default is 1)

Auto-ranging alow you to select a range from one button press, this
allows you to set the number of residues either side of the clicked
residue that becomes the selected zone */
void set_refine_auto_range_step(int i);

/*! \brief set the heuristic fencepost for the maximum number of
  residues in the refinement/regularization residue range  

  Default is 20

*/
void set_refine_max_residues(int n);

/*! \brief refine a zone based on atom indexing */
void refine_zone_atom_index_define(int imol, int ind1, int ind2); 

/*! \brief refine a zone

 presumes that imol_Refinement_Map has been set */
void refine_zone(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);

/*! \brief refine a zone using auto-range

 presumes that imol_Refinement_Map has been set */
void refine_auto_range(int imol, const char *chain_id, int resno1, const char *altconf);

/*! \brief regularize a zone

@return a status, whether the regularisation was done or not.  0 for no, 1 for yes.
  */
int regularize_zone(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);

/*! \brief set the number of refinement steps applied to the
  intermediate atoms each frame of graphics.

  smaller numbers make the movement of the intermediate atoms slower,
  smoother, more elegant.

  Default: 50. */
void set_dragged_refinement_steps_per_frame(int v);

/*! \brief return the number of steps per frame in dragged refinement */
int dragged_refinement_steps_per_frame();

/*! \brief allow refinement of intermediate atoms after dragging,
  before displaying (default: 0, off).

   An attempt to do something like xfit does, at the request of Frank
   von Delft.

   Pass with istate=1 to enable this option. */
void set_refinement_refine_per_frame(int istate);

/*! \brief query the state of the above option */
int refinement_refine_per_frame_state(); 

/*! \brief turn on Ramachandran angles refinement in refinement and regularization */
void set_refine_ramachandran_angles(int state);

int refine_ramachandran_angles_state();

void set_numerical_gradients(int istate);


/*! \brief correct the sign of chiral volumes before commencing refinement?
   
   Do we want to fix chiral volumes (by moving the chiral atom to the
   other side of the chiral plane if necessary).  Default yes
   (1). Note: doesn't work currently. */
void set_fix_chiral_volumes_before_refinement(int istate);

/*! \brief query the state of the above option */
void check_chiral_volumes(int imol);

void check_chiral_volumes_from_widget(GtkWidget *window); 
void fill_chiral_volume_molecule_option_menu(GtkWidget *w);
void chiral_volume_molecule_option_menu_item_select(GtkWidget *item, GtkPositionType pos);

/*! \brief For experienced Cooters who don't like Coot nannying about
  chiral volumes during refinement. */
void set_show_chiral_volume_errors_dialog(short int istate);

/*! \brief set the type of secondary structure restraints

0 no sec str restraints

1 alpha helix restraints

2 beta strand restraints */
void set_secondary_structure_restraints_type(int itype);

/*! \brief return the secondary structure restraints type */
int secondary_structure_restraints_type();

/*! \brief the molecule number of the map used for refinement 

   @return the map number, if it has been set or there is only one
   map, return -1 on no map set (ambiguous) or no maps.
*/
   
int imol_refinement_map();	/* return -1 on no map */

/*! \brief set the molecule number of the map to be used for
  refinement/fitting.

   @return imol on success, -1 on failure*/
int set_imol_refinement_map(int imol);	/* returns imol on success, otherwise -1 */

/*! \brief Does the residue exist? (Raw function) 

   @return 0 on not-exist, 1 on does exist.
*/
int does_residue_exist_p(int imol, char *chain_id, int resno, char *inscode); 

/*  ----------------------------------------------------------------------- */
/*                  Restraints editor                                       */
/*  ----------------------------------------------------------------------- */
#if (GTK_MAJOR_VERSION > 1) 
GtkWidget *wrapped_create_residue_editor_select_monomer_type_dialog();
void show_restraints_editor(const char *monomer_type);
void show_restraints_editor_by_index(int menu_item_index);
#endif

void write_restraints_cif_dictionary(const char *monomer_type, const char *file_name);

/* \} */

/*  ----------------------------------------------------------------------- */
/*               Simplex Refinement                                         */
/*  ----------------------------------------------------------------------- */
/*! \name Simplex Refinement Interface */
/*! \{ */

/*! \brief refine residue range using simplex optimization */
void
fit_residue_range_to_map_by_simplex(int res1, int res2, char *altloc, char *chain_id, int imol, int imol_for_map); 

/*! \brief simply score the residue range fit to map */
float
score_residue_range_fit_to_map(int res1, int res2, char *altloc, char *chain_id, int imol, int imol_for_map); 
/*! \} */

/*  ----------------------------------------------------------------------- */
/*               Nomenclature Errors                                        */
/*  ----------------------------------------------------------------------- */
/*! \name Nomenclature Errors */
/* \{ */
/*! \brief fix nomenclature errors in molecule number imol

   @return the number of resides altered. */
int fix_nomenclature_errors(int imol);

/* \} */

/*  ----------------------------------------------------------------------- */
/*               Move Molecule Here                                        */
/*  ----------------------------------------------------------------------- */
/*! \name move molecule here (wrapper to scheme function) */
/* { */
GtkWidget *wrapped_create_move_molecule_here_dialog();
void move_molecule_here_by_widget(GtkWidget *w);
int move_molecule_to_screen_centre_internal(int imol);
void fill_move_molecule_here_dialog(GtkWidget *w);

/* } */

/*  ----------------------------------------------------------------------- */
/*               Atom info                                                  */
/*  ----------------------------------------------------------------------- */
/* section Atom Info Interface */
/*! \name Atom Info  Interface */
/* \{ */

/*! \brief output to the terminal the Atom Info for the give atom specs
 */
void
output_atom_info_as_text(int imol, const char *chain_id, int resno,
			 const char *ins_code, const char *atname,
			 const char *altconf);

/* \} */

/*  ----------------------------------------------------------------------- */
/*               (Eleanor's) Residue info                                   */
/*  ----------------------------------------------------------------------- */
/* section Residue Info */
/*! \name Residue Info */
/* \{ */
/* Similar to above, we need only one click though. */
void do_residue_info_dialog();
void output_residue_info_dialog    (int atom_index, int imol); /* widget version */
/* scripting version */
void residue_info_dialog(int imol, const char *chain_id, int resno, const char *ins_code); 
int residue_info_dialog_is_displayed();
void output_residue_info_as_text(int atom_index, int imol); /* text version */
/* functions that uses mmdb_manager functions/data types moved to graphics_info_t */

void apply_residue_info_changes(GtkWidget *widget);
void do_distance_define();
void do_angle_define();
void do_torsion_define();
void residue_info_apply_all_checkbutton_toggled();
GtkWidget *wrapped_create_residue_info_dialog();
void clear_residue_info_edit_list();

/* a graphics_info_t function wrapper: */
void residue_info_release_memory(GtkWidget *widget); 
void unset_residue_info_widget(); 
void clear_simple_distances();
void clear_last_simple_distance();
GtkWidget *wrapped_create_geometry_dialog();
void store_geometry_dialog(GtkWidget *w);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  residue enviroment                                      */
/*  ----------------------------------------------------------------------- */
/* section Residue Environment Functions */
/*! \name Residue Environment Functions */
/* \{ */
void fill_environment_widget(GtkWidget *widget);
void execute_environment_settings(GtkWidget *widget);
void toggle_environment_show_distances(GtkToggleButton *button); 
/*! \brief show environment distances.  If state is 0, distances are
  turned off, otherwise distances are turned on. */
void set_show_environment_distances(int state);
void set_show_environment_distances_bumps(int state);
void set_show_environment_distances_h_bonds(int state);
/*! \brief show the state of display of the  environment distances  */
int show_environment_distances_state();
/*! \brief min and max distances for the environment distances */
void set_environment_distances_distance_limits(float min_dist, float max_dist);
/* \} */



/*  ----------------------------------------------------------------------- */
/*                  pointer distances                                      */
/*  ----------------------------------------------------------------------- */
/* section Pointer Functions */
/*! \name Pointer Functions */
/* \{ */
void fill_pointer_distances_widget(GtkWidget *widget);
void execute_pointer_distances_settings(GtkWidget *widget);
void toggle_pointer_distances_show_distances(GtkToggleButton *button); 
/*! \brief turn on (or off) the pointer distance by passing 1 (or 0). */
void set_show_pointer_distances(int istate);
/*! \brief show the state of display of the  pointer distances  */
int  show_pointer_distances_state();
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  zoom                                                    */
/*  ----------------------------------------------------------------------- */
/* section Zoom Functions */
/*! \name Zoom Functions */
/* \{ */
/*! \brief scale the view by f 

   external (scripting) interface (with redraw) 
    @param f the smaller f, the bigger the zoom, typical value 1.3*/
void scale_zoom(float f);  
/* internal interface */
void scale_zoom_internal(float f);  
/*! \brief return the current zoom factor */
float zoom_factor(); 

/*! \brief set smooth scroll with zoom 
   @param i 0 means no, 1 means yes: (default 0) */
void set_smooth_scroll_do_zoom(int i);
/* default 1 (on) */
/*! \brief return the state of the above system */
int      smooth_scroll_do_zoom();   
float    smooth_scroll_zoom_limit(); 
void set_smooth_scroll_zoom_limit(float f);
void set_zoom_adjustment(GtkWidget *w);
void set_zoom(float f);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  CNS data stuff                                          */
/*  ----------------------------------------------------------------------- */
/*! \name CNS Data Functions */
/* \{ */
/*! \brief read CNS data (currently only a placeholder)  */
int handle_cns_data_file(const char *filename, int imol);

/*! \brief read CNS data (currently only a placeholder)

a, b,c are in Angstroems.  alpha, beta, gamma are in degrees.  spg is
the space group info, either ;-delimited symmetry operators or the
space group name*/
int handle_cns_data_file_with_cell(const char *filename, int imol, float a, float b, float c, float alpha, float beta, float gamma, const char *spg_info);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  cif stuff                                               */
/*  ----------------------------------------------------------------------- */
/* section mmCIF Functions */
/*! \name mmCIF Functions */
/* dataset stuff */
/* \{ */
int auto_read_cif_data_with_phases(const char *filename); 
int read_cif_data_with_phases_sigmaa(const char *filename); 
int read_cif_data_with_phases_diff_sigmaa(const char *filename); 
int read_cif_data(const char *filename, int imol_coords);
int read_cif_data_2fofc_map(const char *filename, int imol_coords);  
int read_cif_data_fofc_map(const char *filename, int imol_coords);  
int read_cif_data_with_phases_fo_fc(const char *filename); 
int read_cif_data_with_phases_2fo_fc(const char *filename);
int read_cif_data_with_phases_nfo_fc(const char *filename,
				     int map_type);
int read_cif_data_with_phases_fo_alpha_calc(const char *filename);


/*                  cif (geometry) dictionary                            */
void handle_cif_dictionary(const char *filename);
void read_cif_dictionary(const char *filename);
int write_connectivity(const char* monomer_name, const char *filename);
/*! \brief open the cif dictionary file selector dialog */
void open_cif_dictionary_file_selector_dialog(); 

#ifdef __cplusplus
#ifdef USE_GUILE
SCM non_standard_residue_names_scm(int imol);
#endif
#ifdef USE_PYTHON
PyObject *non_standard_residue_names_py(int imol);
#endif /* USE_PYTHON */
#endif /* c++ */

/* Use the environment variable COOT_REFMAC_LIB_DIR to find cif files
   in subdirectories and import them all. */
void import_all_refmac_cifs(); 
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  SHELX stuff                                             */
/*  ----------------------------------------------------------------------- */
/* section SHELXL Functions */
/*! \name SHELXL Functions */
/* \{ */
/*! \brief read a SHELXL .ins file */
int read_shelx_ins_file(const char *filename, short int recentre_flag);
/*! \brief write a SHELXL .ins file for molecule number imol */
int write_shelx_ins_file(int imol, const char *filename);
/* for shelx fcf file that needs to be filtered: */
int handle_shelx_fcf_file_internal(const char *filename);
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief @return the chain id for the given residue.  Return #f if
  can't do it/fail. */
SCM chain_id_for_shelxl_residue_number(int imol, int resno);
#endif /* USE_GUILE */
/* return 1 for yes, 0 for invalid imol or no. */
int is_shelx_molecule(int imol);

#ifdef USE_PYTHON
/*! \brief @return the chain id for the given residue.  Return Py_False if
  can't do it/fail. */
PyObject *chain_id_for_shelxl_residue_number_py(int imol, int resno);
#endif /* USE_PYTHON */

void add_shelx_string_to_molecule(int imol, const char *string);
#endif /* c++ */

/* \} */
/*  ------------------------------------------------------------------------ */
/*                         Validation:                                       */
/*  ------------------------------------------------------------------------ */
/* section Validation Functions */
/*! \name Validation Functions */
/* \{ */
void deviant_geometry(int imol);
short int is_valid_model_molecule(int imol);
short int is_valid_map_molecule(int imol);

void free_geometry_graph(GtkWidget *dialog); /* free the lines in the widget  */
void unset_geometry_graph(GtkWidget *dialog); /* set the graphics info
						 static to NULL, so
						 that we on longer try
						 to update the
						 widget*/

void add_on_validation_graph_mol_options(GtkWidget *menu, const char *type_in);
void my_delete_validaton_graph_mol_option(GtkWidget *widget, void *);
void validation_graph_b_factor_mol_selector_activate (GtkMenuItem     *menuitem,
						      gpointer         user_data);
void validation_graph_geometry_mol_selector_activate (GtkMenuItem     *menuitem,
						      gpointer         user_data);
void validation_graph_omega_mol_selector_activate (GtkMenuItem     *menuitem,
						   gpointer         user_data);
void validation_graph_rotamer_mol_selector_activate (GtkMenuItem     *menuitem,
						   gpointer         user_data);
void validation_graph_density_fit_mol_selector_activate (GtkMenuItem     *menuitem,
						   gpointer         user_data);
void gln_and_asn_b_factor_outlier_mol_selector_activate (GtkMenuItem     *menuitem,
							 gpointer         user_data);
void validation_graph_ncs_diffs_mol_selector_activate (GtkMenuItem     *menuitem,
						       gpointer         user_data);

void probe_mol_selector_activate (GtkMenuItem     *menuitem,
				  gpointer         user_data);

/* These are called right at the beginning (main) */
/* old style not-generic menu initialization */
/* void create_initial_validation_graph_b_factor_submenu(GtkWidget *window1); */
/* void create_initial_validation_graph_geometry_submenu(GtkWidget *window1); */
/* void create_initial_validation_graph_omega_submenu(GtkWidget *window1); */

/*! \brief generate a list of difference map peaks */
void difference_map_peaks(int imol, int imol_coords, float level, int do_positive_level_flag, int do_negative_level_flag); 

void difference_map_peaks_by_widget(GtkWidget *dialog);
void set_difference_map_peaks_widget(GtkWidget *w);

void clear_diff_map_peaks();
GtkWidget *wrapped_create_generate_diff_map_peaks_dialog();

/*! \brief Make a gui for GLN adn ASN B-factor outiers, compairing the
  O and N temperatur factors difference to the distribution of
  temperature factors from the other atoms.  */
void gln_asn_b_factor_outliers(int imol);
#ifdef USE_PYTHON
void gln_asn_b_factor_outliers_py(int imol);
#endif /*  USE_PYTHON */ 

#ifdef __cplusplus 
#ifdef USE_PYTHON 
/*! \brief return a list of map peaks of molecule number imol_map
  above n_sigma.  There will be cluster filtering of the map peaks.
  Return a list of 3d cartestian coordinates or Python False if
  imol_map is not suitable for peak picking. */
PyObject *map_peaks_py(int imol_map, float n_sigma);
PyObject *map_peaks_near_point_py(int imol_map, float n_sigma, float x, float y, float z, float radius);
PyObject *map_peaks_near_point_from_list_py(int imol_map, PyObject *peak_list, float x, float y, float z, float radius);
#endif /*  USE_PYTHON */

#ifdef USE_GUILE
/*! \brief return a list of map peaks of molecule number imol_map
  above n_sigma.  There will be cluster filtering of the map peaks.
  Return a list of 3d cartestian coordinates or scheme false if
  imol_map is not suitable for peak picking. */
SCM map_peaks_scm(int imol_map, float n_sigma);
SCM map_peaks_near_point_scm(int imol_map, float n_sigma, float x, float y, float z, float radius);
#endif  /* USE_GUILE */


/* does this live here really? */
#ifdef USE_GUILE
SCM get_torsion_scm(int imol, SCM atom_spec_1, SCM atom_spec_2, SCM atom_spec_3, SCM atom_spec_4);
#endif  /* USE_GUILE */

#ifdef USE_PYTHON
PyObject *get_torsion_py(int imol, PyObject *atom_spec_1, PyObject *atom_spec_2, PyObject *atom_spec_3, PyObject *atom_spec_4);
#endif  /* USE_PYTHON */


#endif /* __cplusplus  */

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  ramachandran plot                                       */
/*  ----------------------------------------------------------------------- */
/* section Ramachandran Plot Functions */
/*! \name Ramachandran Plot Functions */
/* \{ */
/* Note for optionmenu from the main window menubar, we should use
   code like this, rather than the map_colour/attach_scroll_wheel code
   (actually, they are mostly the same, differing only the container
   delete code). */

/*! \brief Ramachandran plot for molecule number imol */
void do_ramachandran_plot(int imol);

/*! \brief set the number of biggest difference arrows on the Kleywegt
  plot.  */
void set_kleywegt_plot_n_diffs(int n_diffs);


/*  the the menu */
void add_on_rama_choices();

/*! \brief set the contour levels for theremachandran plot, default
  values are 0.02 (prefered) 0.002 (allowed) */
void set_ramachandran_plot_contour_levels(float level_prefered, float level_allowed);
/*! \brief set the ramachandran plot background block size. 

  Smaller is smoother but slower.  Should be divisible exactly into
  360.  Default value is 10. */
void set_ramachandran_plot_background_block_size(float blocksize) ;

void my_delete_ramachandran_mol_option(GtkWidget *widget, void *);

/* call with value non-zero for on, 0 for off/not. */

/* This should not be used for scripting. */
 
/*  If called with 0, it checks to see if it was previously non-zero, */
/*  if so, then it does a get_user_data to find the pointer to the */
/*  object and deletes it. */
void set_dynarama_is_displayed(GtkWidget *dynarama_widget, int imol);
GtkWidget *dynarama_is_displayed_state(int imol);

/*  return -1 on error. */
int get_mol_from_dynarama(GtkWidget *window);

void set_moving_atoms(double phi, double psi);

void accept_phi_psi_moving_atoms();
void setup_edit_phi_psi(short int state);	/* a button callback */


/* no need to export this to scripting interface */
void setup_dynamic_distances(short int state); 

void destroy_edit_backbone_rama_plot();

/*! \brief  2 molecule ramachandran plot (NCS differences) a.k.a. A Kleywegt Plot. */
void ramachandran_plot_differences(int imol1, int imol2);

/*! \brief  A chain-specific Kleywegt Plot. */
void ramachandran_plot_differences_by_chain(int imol1, int imol2, 
					    const char *a_chain, const char *b_chain);

GtkWidget *wrapped_ramachandran_plot_differences_dialog();
int  do_ramachandran_plot_differences_by_widget(GtkWidget *w); /* return status */
void fill_ramachandran_plot_differences_option_menu_with_chain_options(GtkWidget *chain_optionmenu, 
								       int is_first_mol_flag);
void ramachandran_plot_differences_mol_option_menu_activate_first(GtkWidget *item, GtkPositionType pos);
void ramachandran_plot_differences_mol_option_menu_activate_second(GtkWidget *item, GtkPositionType pos);
void ramachandran_plot_differences_chain_option_menu_activate_first(GtkWidget *item, GtkPositionType pos);
void ramachandran_plot_differences_chain_option_menu_activate_second(GtkWidget *item, GtkPositionType pos);
/* \} */

/*  ----------------------------------------------------------------------- */
/*           sequence_view                                                  */
/*  ----------------------------------------------------------------------- */
/*! \name Sequence View Interface  */
/* \{ */
/*! \brief display the sequence view dialog for molecule number imol */
void do_sequence_view(int imol);
void add_on_sequence_view_choices();
void set_sequence_view_is_displayed(GtkWidget *widget, int imol); 
/* \} */

/*  ----------------------------------------------------------------------- */
/*           rotate moving atoms peptide                                    */
/*  ----------------------------------------------------------------------- */
void change_peptide_carbonyl_by(double angle);/*  in degrees. */
void change_peptide_peptide_by(double angle);  /* in degress */
void execute_setup_backbone_torsion_edit(int imol, int atom_index);
void setup_backbone_torsion_edit(short int state);

void set_backbone_torsion_peptide_button_start_pos(int ix, int iy);
void change_peptide_peptide_by_current_button_pos(int ix, int iy);
void set_backbone_torsion_carbonyl_button_start_pos(int ix, int iy);
void change_peptide_carbonyl_by_current_button_pos(int ix, int iy);

/*  ----------------------------------------------------------------------- */
/*                  atom labelling                                          */
/*  ----------------------------------------------------------------------- */
/* The guts happens in molecule_class_info_t, here is just the
   exported interface */
/* section Atom Labelling */
/*! \name Atom Labelling */
/* \{ */
/*  Note we have to search for " CA " etc */
int    add_atom_label(int imol, char *chain_id, int iresno, char *atom_id); 
int remove_atom_label(int imol, char *chain_id, int iresno, char *atom_id); 
void remove_all_atom_labels();

void set_label_on_recentre_flag(int i); /* 0 for off, 1 or on */

int centre_atom_label_status();

/*! \brief use brief atom names for on-screen labels

 call with istat=1 to use brief labels, istat=0 for normal labels */
void set_brief_atom_labels(int istat);

/*! \brief the brief atom label state */
int brief_atom_labels_state();
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  scene rotation                                          */
/*  ----------------------------------------------------------------------- */
/* section Screen Rotation */
/*! \name Screen Rotation */
/* stepsize in degrees */
/* \{ */
/*! \brief rotate view round y axis stepsize degrees for nstep such steps */
void rotate_y_scene(int nsteps, float stepsize); 
/*! \brief rotate view round x axis stepsize degrees for nstep such steps */
void rotate_x_scene(int nsteps, float stepsize);
/*! \brief rotate view round z axis stepsize degrees for nstep such steps */
void rotate_z_scene(int nsteps, float stepsize);

/*! \brief Bells and whistles rotation 

    spin, zoom and translate.

    where axis is either x,y or z,
    stepsize is in degrees, 
    zoom_by and x_rel etc are how much zoom, x,y,z should 
            have changed by after nstep steps.
*/
void spin_zoom_trans(int axis, int nstep, float stepsize, float zoom_by, 
		     float x_rel, float y_rel, float z_rel);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Views                                                   */
/*  ----------------------------------------------------------------------- */
/*! \brief return the view number */
int add_view_here(const char *view_name); 
/*! \brief return the view number */
int add_view_raw(float rcx, float rcy, float rcz, float quat1, float quat2, 
		 float quat3, float quat4, float zoom, const char *view_name);
void play_views();
void remove_this_view();
/*! \brief the view with the given name */
int remove_named_view(const char *view_name);
/*! \brief the given view number */
void remove_view(int view_number);
/* go to the first view.  if snap_to_view_flag, go directly, else
   smooth twisty path.*/
int go_to_first_view(int snap_to_view_flag);
int go_to_view_number(int view_number, int snap_to_view_flag);
int add_spin_view(const char *view_name, int n_steps, float degrees_total);
/*! \brief Add a view description/annotation to the give view number */
void add_view_description(int view_number, const char *description);
/*! \brief add a view (not add to an existing view) that *does*
  something (e.g. displays or undisplays a molecule) rather than move
  the graphics.

  @return the view number for this (new) view.
 */
int add_action_view(const char *view_name, const char *action_function);
/*! \brief add an action view after the view of the given view number

  @return the view number for this (new) view.
 */
int insert_action_view_after_view(int view_number, const char *view_name, const char *action_function);
int n_views(); 

/*! \brief save views to view_file_name */
void save_views(const char *view_file_name);

float views_play_speed(); 
void set_views_play_speed(float f);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief return the name of the given view, if view_number does not
  specify a view return scheme value False */

SCM view_name(int view_number);
SCM view_description(int view_number);
void go_to_view(SCM view);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
/*! \brief return the name of the given view, if view_number does not
  specify a view return Python value False */
PyObject *view_name_py(int view_number);
PyObject *view_description_py(int view_number);
void go_to_view_py(PyObject *view);
#endif /* USE_PYTHON */
#endif	/* __cplusplus */

/*! \brief Clear the view list */
void clear_all_views();

/* movies */
void set_movie_file_name_prefix(const char *file_name);
void set_movie_frame_number(int frame_number);
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM movie_file_name_prefix();
#endif
#ifdef USE_PYTHON
PyObject *movie_file_name_prefix_py();
#endif /* USE_PYTHON */
#endif /* c++ */
int movie_frame_number();
void set_make_movie_mode(int make_movies_flag);

/*  ----------------------------------------------------------------------- */
/*                  graphics background colour                              */
/*  ----------------------------------------------------------------------- */
/* section Background Colour */
/*! \name Background Colour */
/* \{ */

/*! \brief set the background colour

 red, green and blue are numbers between 0.0 and 1.0 */
void set_background_colour(double red, double green, double blue); 

/*! \brief re draw the background colour when switching between mono and stereo */
void redraw_background();

/*! \brief is the background black (or nearly black)?

@return 1 if the background is black (or nearly black), 
else return 0. */
int  background_is_black_p();
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  ligand fitting stuff                                    */
/*  ----------------------------------------------------------------------- */
/* section Ligand Fitting Functions */
/*! \name   Ligand Fitting Functions */
/* \{ */

/*! \brief set the fraction of atoms which must be in positive density
  after a ligand fit */
void set_ligand_acceptable_fit_fraction(float f);

/*! \brief set the default sigma level that the map is searched to
  find potential ligand sites */
void set_ligand_cluster_sigma_level(float f); /* default 2.2 */

/*! \brief set the number of conformation samples 

    big ligands require more samples.  Default 10.*/
void set_ligand_flexible_ligand_n_samples(int i); /* default 50: Really? */
void set_ligand_verbose_reporting(int i); /* 0 off (default), 1 on */

/*! \brief search the top n sites for ligands.

   Default 10. */
void set_find_ligand_n_top_ligands(int n); /* fit the top n ligands,
					      not all of them, default
					      10. */
/*! \brief how shall we treat the waters during ligand fitting? 

   pass with istate=1 for waters to mask the map in the same way that
   protein atoms do.
   */
void set_find_ligand_mask_waters(int istate);

/* get which map to search, protein mask and ligands from button and
   then do it*/

/*  extract the sigma level and stick it in */
/*  graphics_info_t::ligand_cluster_sigma_level */
void set_ligand_cluster_sigma_level_from_widget(GtkWidget *button);

/*! \brief set the protein molecule for ligand searching */
void set_ligand_search_protein_molecule(int imol);
/*! \brief set the map molecule for ligand searching */
void set_ligand_search_map_molecule(int imol_map);
/*! \brief add a rigid ligand molecule to the list of ligands to search for
  in ligand searching */
void add_ligand_search_ligand_molecule(int imol_ligand);
/*! \brief add a flexible ligand molecule to the list of ligands to search for
  in ligand searching */
void add_ligand_search_wiggly_ligand_molecule(int imol_ligand);

#ifdef __cplusplus
#ifdef USE_GUILE
SCM execute_ligand_search();  
#endif
#ifdef USE_PYTHON
PyObject *execute_ligand_search_py();  
#endif /* USE_PYTHON */
#endif /* __cplusplus */
void free_ligand_search_user_data(GtkWidget *button); 
void add_ligand_clear_ligands(); 

/*! \brief this sets the flag to have expert option ligand entries in
  the Ligand Searching dialog */
void ligand_expert(); 

/*! \brief display the find ligands dialog, if maps, coords and ligands are available */
void do_find_ligands_dialog();

/* Widget functions */

int fill_ligands_dialog(GtkWidget *dialog);
int fill_ligands_dialog_map_bits(GtkWidget *dialog, short int diff_maps_only_flag);	
int fill_ligands_dialog_protein_bits(GtkWidget *dialog);	
int fill_ligands_dialog_ligands_bits(GtkWidget *dialog);	
/*  we need to delete the find_ligand_dialog when we are done, so  */
/*  add this pointer as user data. */
void do_find_ligand_many_atoms_in_ligands(GtkWidget *find_ligand_dialog); 
/* these I factored out, they can be used for the waters dialog too */
int fill_ligands_dialog_map_bits_by_dialog_name(GtkWidget *find_ligand_dialog,
						const char *dialog_name, 
						short int diff_maps_only_flag); 
int fill_ligands_dialog_protein_bits_by_dialog_name(GtkWidget *find_ligand_dialog,
						    const char *dialog_name); 
int fill_vbox_with_coords_options_by_dialog_name(GtkWidget *find_ligand_dialog,
						 const char *dialog_name,
						 short int have_ncs_flag);

void fill_ligands_sigma_level_entry(GtkWidget *dialog);
void fill_ligands_expert_options(GtkWidget *find_ligand_dialog);
void set_ligand_expert_options_from_widget(GtkWidget *button);

/*! \brief "Template"-based matching.  Overlap the first residue in
  imol_ligand onto the residue specified by the reference parameters.
  Use graph matching, not atom names.  

@return success status, False = failed to find residue in either
imol_ligand or imo_ref, 

otherwise return the RT operator */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM overlap_ligands(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
SCM analyse_ligand_differences(int imol_ligand, int imol_ref, const char *chain_id_ref,
			       int resno_ref);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *overlap_ligands_py(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
PyObject *analyse_ligand_differences_py(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
#endif /* PYTHON*/
#endif	/* __cplusplus */

void execute_get_mols_ligand_search(GtkWidget *button); 
/*  info is stored in graphics_info_t beforehand */

/* This has pointers to Coord_orths poked into it, let's clear them
   up. */
void  free_blob_dialog_memory(GtkWidget *w);

/*! \brief flip the ligand (usually active residue) around its eigen vectors
   to the next flip number.  Immediate replacement (like flip
   peptide). */
void flip_ligand(int imol, const char *chain_id, int resno);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  water fitting                                           */
/*  ----------------------------------------------------------------------- */
/* section Water Fitting Functions */
/*! \name Water Fitting Functions */
/* \{ */

/*! Renumber the waters of molecule number imol with consecutive numbering */
void renumber_waters(int imol); 
void wrapped_create_find_waters_dialog();
void fill_find_waters_dialog(GtkWidget *find_ligand_dialog);
/* interface fluff */
void execute_find_waters(GtkWidget *ok_button);  

/*! \brief find waters */
void execute_find_waters_real(int imol_for_map,
			      int imol_for_protein,
			      short int new_waters_mol_flag, 
			      float sigma_cut_off);

void find_waters(int imol_for_map,
		 int imol_for_protein,
		 short int new_waters_mol_flag, 
		 float sigma_cut_off,
		 short int show_blobs_dialog);


/*! \brief move waters of molecule number imol so that they are around the protein.

@return the number of moved waters. */
int move_waters_to_around_protein(int imol);


/*! \brief return the maximum minimum distance of any water atom to
  any protein atom - used in validation of
  move_waters_to_around_protein() funtion.*/
float max_water_distance(int imol);

char *get_text_for_find_waters_sigma_cut_off();
void set_value_for_find_waters_sigma_cut_off(float f); 
/* #ifdef __cplusplus */
/* Just too painful... */
/* void wrapped_create_big_blobs_dialog(const std::vector<Cartesian> &blobs); */
/* #endif */
void on_big_blob_button_clicked(GtkButton *button,
				gpointer user_data);

/*! \brief set the limit of interesting variance, above which waters
  are listed (otherwise ignored)

default 0.12. */
void set_water_check_spherical_variance_limit(float f);

/*! \brief set ligand to protein distance limits

 f1 is the minimum distance, f2 is the maximum distance */
void set_ligand_water_to_protein_distance_limits(float f1, float f2);

/*! \brief set the number of cycles of water searching */
void set_ligand_water_n_cycles(int i); 
void set_write_peaksearched_waters(); 

/*! \brief find blobs  */
void execute_find_blobs(int imol_model, int imol_for_map, float cut_off, short int interactive_flag);

void execute_find_blobs_from_widget(GtkWidget *dialog);

GtkWidget *wrapped_create_unmodelled_blobs_dialog();



/* \} */

/*  ----------------------------------------------------------------------- */
/*                  bond representation                                     */
/*  ----------------------------------------------------------------------- */
/* section Bond Representation */
/*! \name Bond Representation */
/* \{ */

/*! \brief set the default thickness for bonds (e.g. in ~/.coot)  */
void set_default_bond_thickness(int t);

/*! \brief set the thickness of the bonds in molecule number imol to t pixels  */
void set_bond_thickness(int imol, float t);
/*! \brief set the thickness of the bonds of the intermediate atoms to t pixels  */
void set_bond_thickness_intermediate_atoms(float t);
void set_unbonded_atom_star_size(float f);

/*! \brief get the default thickness for bonds*/
int get_default_bond_thickness();

/*! \brief set status of drawing zero occupancy markers.

  default status is 1. */
void set_draw_zero_occ_markers(int status);

/*! \brief set the hydrogen drawing state. istat = 0 is hydrogens off,
  istat = 1: show hydrogens */
void set_draw_hydrogens(int imol, int istat);

/*! \brief the state of draw hydrogens for molecule number imol.  

return -1 on bad imol.  */
int draw_hydrogens_state(int imol);

/*! \brief draw molecule number imol as CAs */
void graphics_to_ca_representation   (int imol);
/*! \brief draw molecule number imol as CA + ligands */
void graphics_to_ca_plus_ligands_representation   (int imol);
/*! \brief draw molecule number imol with no waters */
void graphics_to_bonds_no_waters_representation(int imol);
/*! \brief draw molecule number imol with normal bonds */
void graphics_to_bonds_representation(int mol);
/*! \brief draw molecule number imol with CA bonds in secondary
  structure representation and ligands */
void graphics_to_ca_plus_ligands_sec_struct_representation(int imol);
/*! \brief draw molecule number imol with bonds in secondary structure
  representation */
void graphics_to_sec_struct_bonds_representation(int imol); 
/*! \brief draw molecule number imol in Jones' Rainbow */
void graphics_to_rainbow_representation(int imol);
/*! \brief draw molecule number imol coloured by B-factor */
void graphics_to_b_factor_representation(int imol);
/*! \brief draw molecule number imol coloured by B-factor, CA + ligands */
void graphics_to_b_factor_cas_representation(int imol);
/*! \brief draw molecule number imol coloured by occupancy */
void graphics_to_occupancy_representation(int imol);
/*! \brief what is the bond drawing state of molecule number imol  */
int graphics_molecule_bond_type(int imol); 
/*! \brief scale the colours for colour by b factor representation */
int set_b_factor_bonds_scale_factor(int imol, float f);


GtkWidget *wrapped_create_bond_parameters_dialog();
void apply_bond_parameters(GtkWidget *w);

/*! \brief make a ball and stick representation of imol given atom selection

e.g. (make-ball-and-stick 0 "/1" 0.15 0.25 1) */
int make_ball_and_stick(int imol,
			const char *atom_selection_str,
			float bond_thickness, float sphere_size,
			int do_spheres_flag);
/*! \brief clear ball and stick representation of molecule number imol */
int clear_ball_and_stick(int imol);

/* \brief display/undisplay the given additional representation  */
void set_show_additional_representation(int imol, int representation_number, int on_off_flag);

/* \brief display/undisplay all the additional representations for the given molecule  */
void set_show_all_additional_representations(int imol, int on_off_flag);

/* delete a given additional representation */
void delete_additional_representation(int imol, int representation_number);

/* return the index of the additional representation.  Return -1 on error */
int additional_representation_by_string(int imol,  const char *atom_selection, 
					int representation_type, 
					int bonds_box_type,
					float bond_width,
					int draw_hydrogens_flag);

/* return the index of the additional representation.  Return -1 on error */
int additional_representation_by_attributes(int imol,  const char *chain_id, 
					    int resno_start, int resno_end, 
					    const char *ins_code,
					    int representation_type, 
					    int bonds_box_type,
					    float bond_width,
					    int draw_hydrogens_flag);

GtkWidget *wrapped_create_add_additional_representation_gui();
void add_additional_representation_by_widget(GtkWidget *w);
void add_reps_molecule_option_menu_item_select(GtkWidget *item, GtkPositionType pos);


#ifdef __cplusplus

#ifdef USE_GUILE
SCM additional_representation_info_scm(int imol); 
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *additional_representation_info_py(int imol); 
#endif	/* USE_PYTHON */

#endif	/* __cplusplus */





/* \} */

/*  ----------------------------------------------------------------------- */
/*                  dots display                                            */
/*  ----------------------------------------------------------------------- */

/*! \name Dots Representation */
/* \{ */
/*! \brief display dotted surface

return a generic objects handle (which can be used to remove later) */
int dots(int imol,
	  const char *atom_selection_str,
	  float dot_density, float sphere_size_scale);

/*! \brief clear dots in imol with dots_handle */
void clear_dots(int imol, int dots_handle);

/*! \brief return the number of dots sets for molecule number imol */
int n_dots_sets(int imol); 
/* \} */


/*  ----------------------------------------------------------------------- */
/*                  pepflip                                                 */
/*  ----------------------------------------------------------------------- */
/* section Pep-flip Interface */
/*! \name Pep-flip Interface */
/* \{ */
void do_pepflip(short int state); /* sets up pepflip, ready for atom pick. */
/*! \brief pepflip the given residue */
/* the residue with CO, for scripting interface. */
void pepflip(int imol, const char *chain_id, int resno, const char *inscode, 
	     const char *altconf); 
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  rigid body refinement                                   */
/*  ----------------------------------------------------------------------- */
/* section Rigid Body Refinement Interface */
/*! \name Rigid Body Refinement Interface */
/* \{ */
/* a gui-based interface: setup the rigid body refinement.*/
void do_rigid_body_refine(short int state);	/* set up for atom picking */

/*! \brief setup rigid body refine zone

   where we set the atom selection
   holders according to the arguments and then call
   execute_rigid_body_refine() */
void rigid_body_refine_zone(int reso_start, int resno_end, 
			    const char *chain_id, int imol);

void
rigid_body_refine_by_atom_selection(int imol, 
				    const char *atom_selection_string);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief rigid body refine using residue ranges.  residue_ranges is
    a list of residue ranges.  A residue range is (list chain-id
    resno-start resno-end). */
SCM rigid_body_refine_by_residue_ranges_scm(int imol, SCM residue_ranges); 
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *
rigid_body_refine_by_residue_ranges_py(int imol, PyObject *residue_ranges); 
#endif /* USE_PYTHON */
#endif /* __cplusplus */

void execute_rigid_body_refine(short int auto_range_flag); /* atom picking has happened.
				     Actually do it */


/*! \brief set rigid body fraction of atoms in positive density

 @param f in the range 0.0 -> 1.0 (default 0.75) */
void set_rigid_body_fit_acceptable_fit_fraction(float f);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  dynamic map                                             */
/*  ----------------------------------------------------------------------- */
/* section Dynamic Map */
/*! \name Dynamic Map */
/* \{ */
void   toggle_dynamic_map_display_size();
void   toggle_dynamic_map_sampling();
void   set_map_dynamic_map_sampling_checkbutton(GtkWidget *checkbutton);
void   set_map_dynamic_map_display_size_checkbutton(GtkWidget *checkbutton);
/* scripting interface: */
void set_dynamic_map_size_display_on();
void set_dynamic_map_size_display_off();
int get_dynamic_map_size_display();
void set_dynamic_map_sampling_on();
void set_dynamic_map_sampling_off();
int get_dynamic_map_sampling();
void set_dynamic_map_zoom_offset(int i);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  build one residue by phi/psi search                     */
/*  ----------------------------------------------------------------------- */
/* section Add Terminal Residue Functions */
/*! \name Add Terminal Residue Functions */
/* \{ */
void do_add_terminal_residue(short int state);
/*  execution of this is in graphics_info_t because it uses a CResidue */
/*  in the interface and we can't have that in c-interface.h */
/*  (compilation of coot_wrap_guile goes mad on inclusion of */
/*   mmdb_manager.h) */
void set_add_terminal_residue_n_phi_psi_trials(int n); 
/* Add Terminal Residues actually build 2 residues, this allows us to
   see both residues - default is 0 (off). */
void set_add_terminal_residue_add_other_residue_flag(int i);
void set_terminal_residue_do_rigid_body_refine(short int v); 
int add_terminal_residue_immediate_addition_state(); 

/*! \brief set immediate addition of terminal residue

call with i=1 for immediate addtion */
void set_add_terminal_residue_immediate_addition(int i);

/*! \brief Add a terminal residue

return 0 on failure, 1 on success */
int add_terminal_residue(int imol, char *chain_id, int residue_number,
			 char *residue_type, int immediate_add); 

/*! \brief set the residue type of an added terminal residue.   */
void set_add_terminal_residue_default_residue_type(const char *type);
/*! \brief set a flag to run refine zone on terminal residues after an
  addition.  */
void set_add_terminal_residue_do_post_refine(short int istat); 
/*! \brief what is the value of the previous flag? */
int add_terminal_residue_do_post_refine_state();

#ifdef __cplusplus 
#ifdef USE_GUILE
SCM find_terminal_residue_type(int imol, const char *chain_id, int resno);
#endif 
#ifdef USE_PYTHON
PyObject *find_terminal_residue_type_py(int imol, const char *chain_id, int resno);
#endif /* PYTHON */
#endif /* c++ */

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  delete residue                                          */
/*  ----------------------------------------------------------------------- */
/* section Delete Residues */
/*! \name Delete Residues */
/* in build */
/* by graphics */
/* \{ */
void delete_atom_by_atom_index(int imol, int index, short int do_delete_dialog);
void delete_residue_by_atom_index(int imol, int index, short int do_delete_dialog);
void delete_residue_hydrogens_by_atom_index(int imol, int index, short int do_delete_dialog);

/*! \brief delete residue range */
void delete_residue_range(int imol, const char *chain_id, int resno_start, int end_resno);

/*! \brief delete residue  */
void delete_residue(int imol, const char *chain_id, int resno, const char *inscode); 
/*! \brief delete residue with altconf  */
void delete_residue_with_altconf(int imol, const char *chain_id, int resno, const char *inscode, const char *altloc); 
/*! \brief delete hydrogen atoms in residue  */
void delete_residue_hydrogens(int imol, const char *chain_id, int resno, const char *inscode, const char *altloc); 
/*! \brief delete atom in residue */
void delete_atom(int imol, const char *chain_id, int resno, const char *ins_code, const char *at_name, const char *altloc);
/*! \brief delete all atoms in residue that are not main chain or CB */
void delete_residue_sidechain(int imol, const char *chain_id, int resno, const char*ins_code, 
			      short int do_delete_dialog);
/* toggle callbacks */
void set_delete_atom_mode();
void set_delete_residue_mode();
void set_delete_residue_zone_mode();
void set_delete_residue_hydrogens_mode();
void set_delete_water_mode();
void set_delete_sidechain_mode();
short int delete_item_mode_is_atom_p(); /* (predicate) a boolean */
short int delete_item_mode_is_residue_p(); /* predicate again */
short int delete_item_mode_is_water_p();
short int delete_item_mode_is_sidechain_p();
void store_delete_item_widget(GtkWidget *widget);
void clear_pending_delete_item(); /* for when we cancel with picking an atom */
void clear_delete_item_widget();
void store_delete_item_widget_position();
short int delete_item_widget_is_being_shown();
short int delete_item_widget_keep_active_on();

/* We need to set the pending delete flag and that can't be done in
   callback, so this wrapper does it */
GtkWidget *wrapped_create_delete_item_dialog();

/* utility function, moving widget work out of c-interface-build.cc */
void delete_object_handle_delete_dialog(short int do_delete_dialog);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  rotate/translate buttons                                */
/*  ----------------------------------------------------------------------- */
/* section Rotate/Translate Buttons */
/*  sets flag for atom selection clicks */
void do_rot_trans_setup(short int state); 
void do_rot_trans_adjustments(GtkWidget *dialog);
void rot_trans_reset_previous();
void set_rotate_translate_zone_rotates_about_zone_centre(int istate);
void set_rot_trans_object_type(short int rt_type); /* zone, chain, mol */
int get_rot_trans_object_type();

/*  ----------------------------------------------------------------------- */
/*                  cis <-> trans conversion                                */
/*  ----------------------------------------------------------------------- */
void do_cis_trans_conversion_setup(int istate);
void cis_trans_convert(int imol, const char *chain_id, int resno, const char *altconf);

#ifdef __cplusplus	/* need this wrapper, else gmp.h problems in callback.c */
#ifdef USE_GUILE 
/*! \brief return cis_peptide info for imol.

Return a SCM list object of (residue1 residue2 omega) */
SCM cis_peptides(int imol);
#endif // GUILE
#ifdef USE_PYTHON
/*! \brief return cis_peptide info for imol.

Return a Python list object of [residue1, residue2, omega] */
PyObject *cis_peptides_py(int imol);
#endif // PYTHON
#endif 


/*  ----------------------------------------------------------------------- */
/*                  db-main                                                 */
/*  ----------------------------------------------------------------------- */
/* section Mainchain Building Functions */

/*! \name Mainchain Building Functions */
/* \{ */
void do_db_main(short int state); 
/*! \brief CA -> mainchain conversion */
void db_mainchain(int imol,
		  const char *chain_id,
		  int iresno_start,
		  int iresno_end,
		  const char *direction_string);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  close molecule                                          */
/*  ----------------------------------------------------------------------- */
/* section Close Molecule FUnctions */
/*! \name Close Molecule FUnctions */
/* \{ */

void close_molecule(int imol);

/* get the molecule to delete from the optionmenu */
void close_molecule_by_widget(GtkWidget *optionmenu);
void fill_close_option_menu_with_all_molecule_options(GtkWidget *optionmenu);
/* The callback for the above menuitems */
void close_molecule_item_select(GtkWidget *item, GtkPositionType pos); 

/* New version of close molecule */
void new_close_molecules(GtkWidget *window);
GtkWidget *wrapped_create_new_close_molecules_dialog();

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  rotamers                                                */
/*  ----------------------------------------------------------------------- */
/* section Rotatmer Functions */
/*! \name Rotatmer Functions */
/* \{ */

/* functions defined in c-interface-build */

/*! \brief set the mode of rotamer search, options are (ROTAMERSEARCHAUTOMATIC),  
  (ROTAMERSEARCHLOWRES) (aka. "backrub rotamers), 
  (ROTAMERSEARCHHIGHRES) (with rigid body fitting) */
void set_rotamer_search_mode(int mode);

void setup_rotamers(short int state);

/*  display the rotamer option and display the most likely in the graphics as a */
/*  moving_atoms_asc */
void do_rotamers(int atom_index, int imol);

void show_rotamers_dialog(int imol, const char *chain_id, int resno, const char *ins_code, const char *altconf);

/*! \brief For Dunbrack rotamers, set the lowest probability to be
   considered.  Set as a percentage i.e. 1.00 is quite low.  For
   Richardson Rotamers, this has no effect. */
void set_rotamer_lowest_probability(float f);

/*! \brief set a flag: 0 is off, 1 is on */
void set_rotamer_check_clashes(int i); 

/*! \brief auto fit by rotamer search.

   return the score, for some not very good reason.  clash_flag
   determines if we use clashes with other residues in the score for
   this rotamer (or not).  It would be cool to call this from a script
   that went residue by residue along a (newly-built) chain (now available). */
float auto_fit_best_rotamer(int resno, 
			    const char *altloc, 
			    const char *insertion_code, 
			    const char *chain_id, int imol_coords, int imol_map, 
			    int clash_flag, float lowest_probability);

/*! \brief set the clash flag for rotamer search

   And this functions for [pre-setting] the variables for
   auto_fit_best_rotamer called interactively (using a graphics_info_t
   function). 0 off, 1 on.*/
void set_auto_fit_best_rotamer_clash_flag(int i); /*  */
/* currently stub function only */
float rotamer_score(int imol, const char *chain_id, int res_no, const char *insertion_code, 
		    const char *alt_conf);
void setup_auto_fit_rotamer(short int state);	/* called by the Auto Fit button call
				   back, set's in_auto_fit_define. */

/*! \brief return the number of rotamers for this residue - return -1
  on no residue found.*/
int n_rotamers(int imol, const char *chain_id, int resno, const char *ins_code);
/*! \brief set the residue specified to the rotamer number specifed. */
int set_residue_to_rotamer_number(int imol, const char *chain_id, int resno, const char *ins_code, int rotamer_number);

#ifdef __cplusplus
#ifdef USE_GUILE
SCM get_rotamer_name_scm(int imol, const char *chain_id, int resno, const char *ins_code);
#endif 
#ifdef USE_PYTHON
PyObject *get_rotamer_name_py(int imol, const char *chain_id, int resno, const char *ins_code);
#endif /* USE_GUILE */
#endif /* c++ */


/*! \brief fill all the residues of molecule number imol that have
   missing atoms.

To be used to remove the effects of chainsaw.  */
void fill_partial_residues(int imol);

void fill_partial_residue(int imol, const char *chain_id, int resno, const char* inscode);

#ifdef __cplusplus
#ifdef USE_GUILE
SCM missing_atom_info_scm(int imol);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *missing_atom_info_py(int imol);
#endif /* USE_PYTHON */
#endif /* __cplusplus */


/* Used for unsetting the rotamer dialog when it gets destroyed. */
void
set_graphics_rotamer_dialog(GtkWidget *w);


#ifdef __cplusplus	/* need this wrapper, else gmp.h problems in callback.c */
#ifdef USE_GUILE 
/*! \brief Activate rotamer graph analysis for molecule number imol.  

Return rotamer info - function used in testing.  */
SCM rotamer_graphs(int imol);
#endif // USE_GUILE
#ifdef USE_PYTHON
/*! \brief Activate rotamer graph analysis for molecule number imol.  

Return rotamer info - function used in testing.  */
PyObject *rotamer_graphs_py(int imol);
#endif /* USE_PYTHON */
#endif /* c++ */

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  180 degree flip                                         */
/*  ----------------------------------------------------------------------- */
/*! \name 180 Flip Side chain */
/* \{ */

/*! \brief rotate 180 degrees round the last chi angle */
void do_180_degree_side_chain_flip(int imol, const char* chain_id, int resno, 
				   const char *inscode, const char *altconf);

void setup_180_degree_flip(short int state); 
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  mutate                                                  */
/*  ----------------------------------------------------------------------- */
/* section Mutate Functions */
/*! \name Mutate Functions */
/* \{ */

/* c-interface-build */
void setup_mutate(short int state);
/*! \brief Mutate then fit to map

 that we have a map define is checked first */
void setup_mutate_auto_fit(short int state);

void do_mutation(const char *type, short int is_stub_flag);

/* auto-mutate stuff */
short int progressive_residues_in_chain_check(const char *chain_id, int imol);

/*! \brief mutate a given residue 

target_res_type is a three-letter-code.

Return 1 on a good mutate. */
int mutate(int imol, const char *chain_id, int ires, const char *inscode,  const char *target_res_type);

/*! \brief mutate a base. return success status, 1 for a good mutate. */
int mutate_base(int imol, const char *chain_id, int res_no, const char *ins_code, const char *res_type);


/*! \brief Do you want Coot to automatically run a refinement after
  every mutate and autofit?

 1 for yes, 0 for no. */
void set_mutate_auto_fit_do_post_refine(short int istate);

/*! \brief what is the value of the previous flag? */
int mutate_auto_fit_do_post_refine_state();

/*! \brief Do you want Coot to automatically run a refinement after
  every rotamer autofit?

 1 for yes, 0 for no. */
void set_rotamer_auto_fit_do_post_refine(short int istate);

/*! \brief what is the value of the previous flag? */
int rotamer_auto_fit_do_post_refine_state();


/*! \brief an alternate interface to mutation of a singe residue.

 @return 1 on success, 0 on failure 

  ires-ser is the serial number of the residue, not the seqnum 
  There 2 functions don't make backups, but mutate() does - CHECKME
   Hence mutate() is for use as a "one-by-one" type and the following
   2 by wrappers that muate either a residue range or a whole chain

   Note that the target_res_type is a char, not a string (or a char *).  
   So from the scheme interface you'd use (for example) hash
   backslash A for ALA.  */


int mutate_single_residue_by_serial_number(int ires_ser, 
					   const char *chain_id,
					   int imol, char target_res_type);
/* ires is the seqnum of the residue (conventional) */
int mutate_single_residue_by_seqno(int ires, const char *inscode,
				   const char *chain_id,
				   int imol, char target_res_type);

/* an internal function - not useful for scripting: */

void do_base_mutation(const char *type);

/*! \brief set a flag saying that the residue chosen by mutate or
  auto-fit mutate should only be added as a stub (mainchain + CB) */
void set_residue_type_chooser_stub_state(short int istat);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  alternate conformation                                  */
/*  ----------------------------------------------------------------------- */
/* section Alternative Conformation */
/*! \name Alternative Conformation */
/* c-interface-build function */
/*! \{ */

short int alt_conf_split_type_number();
void set_add_alt_conf_split_type_number(short int i);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief add an alternative conformer to a residue.  Add it in
  conformation rotamer number rotamer_number.  

Return the new alt_conf chain_id on sucess, scheme false on fail */
SCM add_alt_conf_scm(int imol, const char*chain_id, int res_no, const char *ins_code, 
		     const char *alt_conf, int rotamer_number);
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief add an alternative conformer to a residue.  Add it in
  conformation rotamer number rotamer_number.  

Return the new alt_conf chain_id on sucess, python False on fail */
PyObject *add_alt_conf_py(int imol, const char*chain_id, int res_no, const char *ins_code, 
		     const char *alt_conf, int rotamer_number);
#endif	/* USE_PYTHON */
#endif /* __cplusplus */

void setup_alt_conf_with_dialog(GtkWidget *dialog); 
void unset_add_alt_conf_dialog(); /* set the static dialog holder in
				     graphics info to NULL */
void unset_add_alt_conf_define(); /* turn off pending atom pick */
void altconf();			/* temporary debugging interface. */
void set_add_alt_conf_new_atoms_occupancy(float f); /* default 0.5 */
void set_show_alt_conf_intermediate_atoms(int i);
int  show_alt_conf_intermediate_atoms_state();
void zero_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2);
void fill_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2);
void set_b_factor_residue_range(int imol, const char *chain_id, int ires1, int ires2, float bval);
void reset_b_factor_residue_range(int imol, const char *chain_id, int ires1, int ires2);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  pointer atoms                                           */
/*  ----------------------------------------------------------------------- */
/* section Pointer Atom Functions */
/*! \name Pointer Atom Functions */
/* c-interface-build */
/*! \{ */

void place_atom_at_pointer();
/* which calls the following gui function (if using non dummies) */
void place_atom_at_pointer_by_window(); 
void place_typed_atom_at_pointer(const char *type);

/* ! \brief set pointer atom is a water (HOH) */
void set_pointer_atom_is_dummy(int i);
void fill_place_atom_molecule_option_menu(GtkWidget *optionmenu);
void display_where_is_pointer(); /* print the coordinates of the
				    pointer to the console */
/*! \brief Return the current pointer atom molecule, create a pointer
  atom molecule if necessary (i.e. when the user has not set it).  */
int create_pointer_atom_molecule_maybe(); 
/*! \brief Return the current pointer atom molecule */
int pointer_atom_molecule(); 
void set_pointer_atom_molecule(int imol);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  baton mode                                              */
/*  ----------------------------------------------------------------------- */
/*! \name Baton Build Interface Functions */
/* section Baton Build Functions */
/* c-interface-build */
/*! \{ */
/*! \brief toggle so that mouse movement moves the baton not rotates the view. */
void set_baton_mode(short int i); /* Mouse movement moves the baton not the view? */
/*! \brief draw the baton or not */
int try_set_draw_baton(short int i); /* draw the baton or not */
/*! \brief accept the baton tip position - a prime candidate for a key binding */
void accept_baton_position();	/* put an atom at the tip */
/*! \brief move the baton tip position - another prime candidate for a key binding */
void baton_try_another();
/*! \brief shorten the baton length */
void shorten_baton();
/*! \brief lengthen the baton */
void lengthen_baton();
/*! \brief delete the most recently build CA position */
void baton_build_delete_last_residue();
/*! \brief set the parameters for the start of a new baton-built fragment. direction can either 
     be "forwards" or "backwards" */
void set_baton_build_params(int istart_resno, const char *chain_id, const char *direction); 
void set_baton_build_params_from_widget(GtkWidget *params_dialog);
void baton_mode_calculate_skeleton(GtkWidget *window);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  post baton mode                                         */
/*  ----------------------------------------------------------------------- */
/* section Post-Baton Functions */
/* c-interface-build */
/* Reverse the direction of a the fragment of the clicked on
   atom/residue.  A fragment is a consequitive range of residues -
   where there is a gap in the numbering, that marks breaks between
   fragments in a chain.  There also needs to be a distance break - if
   the CA of the next/previous residue is more than 5A away, that also
   marks a break. Thow away all atoms in fragment other than CAs.*/
void reverse_direction_of_fragment(int imol, const char *chain_id, int resno);
void setup_reverse_direction(short int i);


/*  ----------------------------------------------------------------------- */
/*                  terminal OXT atom                                       */
/*  ----------------------------------------------------------------------- */
/* section Terminal OXT Atom */
/*! \name Terminal OXT Atom */
/* c-interface-build */
/*! \{ */
short int add_OXT_to_residue(int imol, int reso, const char *insertion_code, const char *chain_id);
GtkWidget *wrapped_create_add_OXT_dialog();
void apply_add_OXT_from_widget(GtkWidget *w);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  crosshairs                                              */
/*  ----------------------------------------------------------------------- */
/* section Crosshairs Interface */
/*! \name Crosshairs  Interface */
/*! \{ */
/*! \brief draw the distance crosshairs, 0 for off, 1 for on. */
void set_draw_crosshairs(short int i);
/* so that we display the crosshairs with the radiobuttons in the
   right state, return draw_crosshairs_flag */
short int draw_crosshairs_state();
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Edit Chi Angles                                         */
/*  ----------------------------------------------------------------------- */
/* section Edit Chi Angles */
/*! \name  Edit Chi Angles */
/* \{ */
/* c-interface-build functions */
void setup_edit_chi_angles(short int state); 

void rotate_chi(float am);

/*! \brief show torsions that rotate hydrogens in the torsion angle
  manipulation dialog.  Note that this may be needed if, in the
  dictionary cif file torsion which have as a 4th atom both a hydrogen
  and a heavier atom bonding to the 3rd atom, but list the 4th atom as
  a hydrogen (not a heavier atom). */
void set_find_hydrogen_torsions(short int state);
void set_graphics_edit_current_chi(int ichi); /* button callback */
void unset_moving_atom_move_chis();

/*! \brief display the edit chi angles gui for the given residue 

 return a status of 0 if it failed to fined the residue, 
 return a value of 1 if it worked. */
int edit_chi_angles(int imol, const char *chain_id, int resno, 
		     const char *ins_code, const char *altconf);

int set_show_chi_angle_bond(int imode);

/* not for user consumption, this finds (from itself) the residue type
   and calls the graphics_info_t function. */
void fill_chi_angles_vbox(GtkWidget *vbox);

/* a callback from the callbacks.c, setting the state of
   graphics_info_t::edit_chi_angles_reverse_fragment flag */
void set_edit_chi_angles_reverse_fragment_state(short int istate);

/* No need for this to be exported to scripting */
/*! \brief beloved torsion general at last makes an entrance onto the
  Coot scene...  */
void setup_torsion_general(short int state);
/* No need for this to be exported to scripting */
void toggle_torsion_general_reverse();

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Backrub                                                 */
/*  ----------------------------------------------------------------------- */
/*! \name Backrubbing function */
/*! \{ */
/* \brief Do a back-rub rotamer search (with autoaccept). 

@return the success status, 0 for fail, 1 for successful fit.  */
int backrub_rotamer(int imol, const char *chain_id, int res_no, 
		    const char *ins_code, const char *alt_conf);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Mask                                                    */
/*  ----------------------------------------------------------------------- */
/* section Masks */
/*! \name Masks */
/*! \{ */
/* The idea is to generate a new map that has been masked by some
   coordinates. */
/*! \brief  generate a new map that has been masked by some coordinates

        (mask-map-by-molecule map-no mol-no invert?)  creates and
        displays a masked map, cuts down density where the coordinates
        are (invert is 0).  If invert? is 1, cut the density down
        where there are no atoms atoms.  */
int mask_map_by_molecule(int map_mol_no, int coord_mol_no, short int invert_flag);

int mask_map_by_atom_selection(int map_mol_no, int coords_mol_no, const char *mmdb_atom_selection, short int invert_flag); 

/*! \brief set the atom radius for map masking */
void set_map_mask_atom_radius(float rad);
/*! \brief get the atom radius for map masking */
float map_mask_atom_radius();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  check waters interface                                  */
/*  ----------------------------------------------------------------------- */
/* section Check Waters Interface */
/*! \name check Waters Interface */
/* interactive check by b-factor, density level etc. */
/*! \{ */
GtkWidget *wrapped_create_check_waters_dialog();
void set_check_waters_b_factor_limit(float f);
void set_check_waters_map_sigma_limit(float f);
void set_check_waters_min_dist_limit(float f);
void set_check_waters_max_dist_limit(float f);
void check_waters_molecule_menu_item_activate(GtkWidget *item, 
					      GtkPositionType pos);
void check_water_by_difference_maps_option_menu_item_select(GtkWidget *item, 
							    GtkPositionType pos);
void do_check_waters_by_widget(GtkWidget *dialog);
void store_checked_waters_baddies_dialog(GtkWidget *dialog);

GtkWidget *wrapped_checked_waters_baddies_dialog(int imol, float b_factor_lim, 
						 float map_sigma_lim, 
						 float min_dist, float max_dist,
						 short int part_occ_contact_flag,
						 short int zero_occ_flag,
						 short int logical_operator_and_or_flag);

/*! \brief Delete waters that are fail to meet the given criteria. */
void delete_checked_waters_baddies(int imol, float b_factor_lim, 
				   float map_sigma_lim, 
				   float min_dist, float max_dist,
				   short int part_occ_contact_flag,
				   short int zero_occ_flag,
				   short int logical_operator_and_or_flag);

/* difference map variance check  */
void check_waters_by_difference_map(int imol_waters, int imol_diff_map, 
				    int interactive_flag); 
/* results widget are in graphics-info.cc  */
/* Let's give access to the sigma level (default 4) */
float check_waters_by_difference_map_sigma_level_state();
void set_check_waters_by_difference_map_sigma_level(float f);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return a list of waters that are coordinated with at least
  coordination_number of other atoms at distances less than or equal
  to dist_max.  Return scheme false on not able to make a list,
  otherwise a list of atoms and neighbours */
SCM highly_coordinated_waters_scm(int imol, int coordination_number, float dist_max);
#endif
#ifdef USE_PYTHON
/*! \brief return a list of waters that are coordinated with at least
  coordination_number of other atoms at distances less than or equal
  to dist_max */
PyObject *highly_coordinated_waters_py(int imol, int coordination_number, float dist_max);
#endif
#endif


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Least squares                                           */
/*  ----------------------------------------------------------------------- */
/* section Least-Squares matching */
/*! \name Least-Squares matching */
/*! \{ */
void clear_lsq_matches();
void add_lsq_match(int reference_resno_start, 
		   int reference_resno_end,
		   const char *chain_id_reference,
		   int moving_resno_start, 
		   int moving_resno_end,
		   const char *chain_id_moving,
		   int match_type); /* 0: all
				       1: main
				       2: CA 
				    */
#ifdef __cplusplus
#ifdef USE_GUILE
void add_lsq_atom_pair_scm(SCM atom_spec_ref, SCM atom_spec_moving);
#endif
#ifdef USE_PYTHON
void add_lsq_atom_pair_py(PyObject *atom_spec_ref, PyObject *atom_spec_moving);
#endif
#endif /* __cplusplus */

#ifdef __cplusplus	/* need this wrapper, else gmp.h problems in callback.c */
#ifdef USE_GUILE
/* Return an rtop pair (proper list) on good match, else #f */
SCM apply_lsq_matches(int imol_reference, int imol_moving);
#endif
#ifdef USE_PYTHON
/* Return an rtop pair (proper list) on good match, else False */
PyObject *apply_lsq_matches_py(int imol_reference, int imol_moving);
#endif // PYTHON
#endif /* __cplusplus */

/* poor old python programmers... */
int apply_lsq_matches_simple(int imol_reference, int imol_moving);
		    
/* section Least-Squares plane interface */
void setup_lsq_deviation(int state);
void setup_lsq_plane_define(int state); 
GtkWidget *wrapped_create_lsq_plane_dialog();
void unset_lsq_plane_dialog(); /* callback from destroy of widget */
void remove_last_lsq_plane_atom();

GtkWidget *wrapped_create_least_squares_dialog();
int apply_lsq_matches_by_widget(GtkWidget *lsq_dialog); /* return 1 for good fit */
void lsq_ref_mol_option_menu_changed(GtkWidget *item, GtkPositionType pos);
void lsq_mov_mol_option_menu_changed(GtkWidget *item, GtkPositionType pos);
void lsq_reference_chain_option_menu_item_activate(GtkWidget *item,
						   GtkPositionType pos);
void lsq_moving_chain_option_menu_item_activate(GtkWidget *item,
						GtkPositionType pos);
void fill_lsq_option_menu_with_chain_options(GtkWidget *chain_optionmenu, 
					     int is_reference_structure_flag);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  trim                                                    */
/*  ----------------------------------------------------------------------- */
/* section Molecule Trimming Interface */
/*! \name Trim */
/*! \{ */

/* a c-interface-build function */
/*! \brief cut off (delete or give zero occupancy) atoms in the given
  molecule if they are below the given map (absolute) level. */
void trim_molecule_by_map(int imol_coords, int imol_map, 
			  float map_level, int delete_or_zero_occ_flag);

/*! \} */


/*  ------------------------------------------------------------------------ */
/*                       povray/raster3d interface                           */
/*  ------------------------------------------------------------------------ */
/* make the text input to external programs */
/*! \name External Ray-Tracing */
/* \{ */

/*! \brief create a r3d file for the current view */
void raster3d(const char *rd3_filename);
void povray(const char *filename);
void renderman(const char *rib_filename);
/* a wrapper for the (scheme) function that makes the image, callable
   from callbacks.c  */
void make_image_raster3d(const char *filename);
void make_image_povray(const char *filename);
#ifdef USE_PYTHON
void make_image_raster3d_py(const char *filename);
void make_image_povray_py(const char *filename);
#endif /* USE_PYTHON */

/*! \brief set the bond thickness for the Raster3D representation  */
void set_raster3d_bond_thickness(float f);
/*! \brief set the atom radius for the Raster3D representation  */
void set_raster3d_atom_radius(float f);
/*! \brief set the density line thickness for the Raster3D representation  */
void set_raster3d_density_thickness(float f);
/*! \brief set the flag to show atoms for the Raster3D representation  */
void set_renderer_show_atoms(int istate);
/*! \brief set the bone (skeleton) thickness for the Raster3D representation  */
void set_raster3d_bone_thickness(float f);
/*! \brief turn off shadows for raster3d output - give argument 0 to turn off  */
void set_raster3d_shadows_enabled(int state);
/*! \brief run raster3d and display the resulting image.  */
void raster_screen_shot(); /* run raster3d or povray and guile */
                           /* script to render and display image */
#ifdef USE_PYTHON
void raster_screen_shot_py(); /* run raster3d or povray and python */
#endif
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  citation notice                                         */
/*  ----------------------------------------------------------------------- */
void citation_notice_off();

/*  ----------------------------------------------------------------------- */
/*                  Superpose                                               */
/*  ----------------------------------------------------------------------- */
/* section Superposition (SSM) */
/*! \name Superposition (SSM) */
/* \{ */

/*! \brief simple interface to superposition. 

Superpose all residues of imol2 onto imol1.  imol1 is reference, we
can either move imol2 or copy it to generate a new molecule depending
on the vaule of move_imol2_flag (1 for move 0 for copy). */
void superpose(int imol1, int imol2, short int move_imol2_flag); 


/*! \brief chain-based interface to superposition. 

Superpose the given chains of imol2 onto imol1.  imol1 is reference,
we can either move imol2 or copy it to generate a new molecule
depending on the vaule of move_imol2_flag (1 for move 0 for copy). */
void superpose_with_chain_selection(int imol1, int imol2, 
				    const char *chain_imol1,
				    const char *chain_imol2,
				    int chain_used_flag_imol1,
				    int chain_used_flag_imol2,
				    short int move_imol2_copy_flag);

/*! \brief detailed interface to superposition. 

Superpose the given atom selection (specified by the mmdb atom
selection strings) of imol2 onto imol1.  imol1 is reference, we can
either move imol2 or copy it to generate a new molecule depending on
the vaule of move_imol2_flag (1 for move 0 for copy). */
void superpose_with_atom_selection(int imol1, int imol2,
				   const char *mmdb_atom_sel_str_1, 
				   const char *mmdb_atom_sel_str_2,
				   short int move_imol2_copy_flag);

void execute_superpose(GtkWidget *w);
GtkWidget *wrapped_create_superpose_dialog(); /* used by callback */
void fill_superpose_option_menu_with_chain_options(GtkWidget *chain_optionmenu, 
 						   int is_reference_structure_flag);
 
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  NCS                                                     */
/*  ----------------------------------------------------------------------- */
/* section NCS */
/*! \name NCS */

/* \{ */
/*! \brief set drawing state of NCS ghosts for molecule number imol   */
void set_draw_ncs_ghosts(int imol, int istate);
/*! \brief return the drawing state of NCS ghosts for molecule number
  imol.  Return -1 on imol is a bad molecule or no ghosts.  */
int draw_ncs_ghosts_state(int imol);

/*! \brief set bond thickness of NCS ghosts for molecule number imol   */
void set_ncs_ghost_bond_thickness(int imol, float f);
/*! \brief update ghosts for molecule number imol */
void ncs_update_ghosts(int imol); // update ghosts
/*! \brief make NCS map */
int make_dynamically_transformed_ncs_maps(int imol_model, int imol_map, 
					  int overwrite_maps_of_same_name_flag);
int make_dynamically_transformed_ncs_maps_by_widget(GtkWidget *dialog);
GtkWidget *wrapped_create_ncs_maps_dialog();
void make_ncs_ghosts_maybe(int imol);
/*! \brief Add NCS matrix */
void add_ncs_matrix(int imol, const char *this_chain_id, const char *target_chain_id,
		    float m11, float m12, float m13, 
		    float m21, float m22, float m23, 
		    float m31, float m32, float m33, 
		    float t1,  float t2,  float t3);

void clear_ncs_ghost_matrices(int imol);

/*! \brief add an NCS matrix for strict NCS molecule representation

for CNS strict NCS usage: expand like normal symmetry does  */
int add_strict_ncs_matrix(int imol,
			  const char *this_chain_id,
			  const char *target_chain_id,
			  float m11, float m12, float m13, 
			  float m21, float m22, float m23, 
			  float m31, float m32, float m33, 
			  float t1,  float t2,  float t3);

/*! \brief return the state of NCS ghost molecules for molecule number imol   */
int show_strict_ncs_state(int imol);
/*! \brief set display state of NCS ghost molecules for molecule number imol   */
void set_show_strict_ncs(int imol, int state);
/*! \brief At what level of homology should we say that we can't see homology
   for NCS calculation? (default 0.8) */
void set_ncs_homology_level(float flev);
/* for a single copy */
/*! \brief Copy single NCS chain */
void copy_chain(int imol, const char *from_chain, const char *to_chain);
/* do multiple copies */
/*! \brief Copy chain from master to all related NCS chains */
void copy_from_ncs_master_to_others(int imol, const char *chain_id);
/* void copy_residue_range_from_ncs_master_to_others(int imol, const char *master_chain_id, int residue_range_start, int residue_range_end); */
/*! \brief Copy residue range to all related NCS chains.  If the
  target residues do not exist in the peer chains, then create
  them. */
void copy_residue_range_from_ncs_master_to_others(int imol, const char *master_chain_id, 
						  int residue_range_start, int residue_range_end);
GtkWidget *wrapped_create_ncs_control_dialog();	
/*! \brief change the NCS master chain  (by number)*/
void ncs_control_change_ncs_master_to_chain(int imol, int ichain); 
/*! \brief change the NCS master chain  (by chain_id)*/
void ncs_control_change_ncs_master_to_chain_id(int imol, const char *chain_id); 
void ncs_control_change_ncs_master_to_chain_update_widget(GtkWidget *w, int imol, int ichain); 
/*! \brief display the NCS master chain  */
void ncs_control_display_chain(int imol, int ichain, int state);

void set_ncs_matrix_type(int flag);
int get_ncs_matrix_state();

#ifdef __cplusplus
#ifdef USE_GUILE
/* Return e.g. ("B" "A" '(((1 "") (1 "") 0.4) ((2 "") (2 "") 0.3))
   i.e. ncs-related-chain its-master-chain-id and a list of residue
   info: (residue number matches: (this-resno this-inscode
   matching-mater-resno matching-master-inscode
   rms-atom-position-differences))) */
SCM ncs_chain_differences_scm(int imol, const char *master_chain_id);

/*! \brief return something like: '(("A" "B")) or '(("A" "C" "E") ("B"
  "D" "F")). The master chain goes in first. 

   If imol does not have NCS ghosts, return #f */
SCM ncs_chain_ids_scm(int imol);
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
/* Return e.g. ["B", "A", [[[1, ""], [1, ""], 0.4], [[2, ""], [2, ""], 0.3]]]
   i.e. ncs-related-chain its-master-chain-id and a list of residue
   info: [residue number matches: [this-resno, this-inscode
   matching-master-resno, matching-master-inscode, 
   rms-atom-position-differences]] */
PyObject *ncs_chain_differences_py(int imol, const char *master_chain_id);

/*! \brief return something like: [["A", "B"]] or [["A", "C", "E"], ["B",
  "D", "F"]]. The master chain goes in first.

   If imol does not have NCS ghosts, return #f */
PyObject *ncs_chain_ids_py(int imol);
#endif  /* USE_PYTHON */

#ifdef USE_GUILE
/* return #f on bad imol or a list of ghosts on good imol.  Can
   include NCS rtops if they are available, else the rtops are #f */
SCM ncs_ghosts_scm(int imol);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
/* return False on bad imol or a list of ghosts on good imol.  Can
   include NCS rtops if they are available, else the rtops are False */
PyObject *ncs_ghosts_py(int imol);
#endif	/* USE_PYTHON */


#endif	/* __cplusplus */

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Autobuild helices and strands                           */
/*  ----------------------------------------------------------------------- */

#define FIND_SECSTRUC_NORMAL 0
#define FIND_SECSTRUC_STRICT 1
#define FIND_SECSTRUC_HI_RES 2
#define FIND_SECSTRUC_LO_RES 3

/*! \name Helices and Strands*/
/*! \{ */
/*! \brief add a helix 

   Add a helix somewhere close to this point in the map, try to fit
   the orientation. Add to a molecule called "Helix", create it if
   needed.  Create another moecule called "Reverse Helix" if the helix
   orientation isn't completely unequivocal.

   @return the index of the new molecule.*/
int place_helix_here();

/*! \brief add a strands

   Add a strand close to this point in the map, try to fit
   the orientation. Add to a molecule called "Strand", create it if needed.
   n_residues is the estimated number of residues in the strand.

   n_sample_strands is the number of strands from the database tested
   to fit into this strand density.  8 is a suggested number.  20 for
   a more rigourous search, but it will be slower.

   @return the index of the new molecule.*/
int place_strand_here(int n_residues, int n_sample_strands);


/*! \brief show the strand placement gui.

  Choose the python version in there, if needed.  Call scripting
  function, display it in place, don't return a widget. */
void   place_strand_here_dialog(); 


/*! \brief autobuild helices

   Find secondary structure in the current map.
   Add to a molecule called "Helices", create it if
   needed.

   @return the index of the new molecule.*/
int find_helices();

/*! \brief autobuild strands

   Find secondary structure in the current map.
   Add to a molecule called "Strands", create it if
   needed.

   @return the index of the new molecule.*/
int find_strands();

/*! \brief autobuild secondary structure

   Find secondary structure in the current map.
   Add to a molecule called "SecStruc", create it if needed.

   @return the index of the new molecule.*/
int find_secondary_structure(
    short int use_helix,  int helix_length,  int helix_target,
    short int use_strand, int strand_length, int strand_target );

/*! \brief autobuild secondary structure

   Find secondary structure local to current view in the current map.
   Add to a molecule called "SecStruc", create it if needed.

   @return the index of the new molecule.*/
int find_secondary_structure_local(
    short int use_helix,  int helix_length,  int helix_target,
    short int use_strand, int strand_length, int strand_target,
    float radius );

/* \} */



/*  ----------------------------------------------------------------------- */
/*             New Molecule by Various Selection                            */
/*  ----------------------------------------------------------------------- */
/*! \brief create a new molecule that consists of only the residue of
  type residue_type in molecule number imol

@return the new molecule number, -1 means an error. */
int new_molecule_by_residue_type_selection(int imol, const char *residue_type);

/*! \brief create a new molecule that consists of only the atoms specified
  by the mmdb atoms selection string in molecule number imol

@return the new molecule number, -1 means an error. */
int new_molecule_by_atom_selection(int imol, const char* atom_selection);

/*! \brief create a new molecule that consists of only the atoms 
  within the given radius (r) of the given position.

@return the new molecule number, -1 means an error. */
int new_molecule_by_sphere_selection(int imol, float x, float y, float z, 
				     float r, short int allow_partial_residues);


/*  ----------------------------------------------------------------------- */
/*                  Miguel's orientation axes matrix                         */
/*  ----------------------------------------------------------------------- */
/* section Miguel's orientation axes matrix */

void
set_axis_orientation_matrix(float m11, float m12, float m13,
			    float m21, float m22, float m23,
			    float m31, float m32, float m33);

void
set_axis_orientation_matrix_usage(int state);



/*  ----------------------------------------------------------------------- */
/*                  RNA/DNA                                                 */
/*  ----------------------------------------------------------------------- */
/* section RNA/DNA */
/*! \name RNA/DNA */

/* \{ */
/*!  \brief create a molecule of idea nucleotides 

use the given sequence (single letter code)

RNA_or_DNA is either "RNA" or "DNA"

form is either "A" or "B"

@return the new molecule number or -1 if a problem */
int ideal_nucleic_acid(const char *RNA_or_DNA, const char *form,
		       short int single_stranged_flag,
		       const char *sequence);

GtkWidget *wrapped_nucleotide_builder_dialog();
void ideal_nucleic_acid_by_widget(GtkWidget *builder_dialog);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE

/*! \brief return #f if residue not found,  otherwise
 (list (list phosphate-distance puckered-atom out-of-plane-distance plane-distortion) chain-id resno ins-code)

 (where plane-distortion is for the other 4 atoms in the plane (I think)).
 
 and if there is no following residue, then the phosphate distance
 cannot be calculated, so the (inner) list is null (not filled).
*/
SCM pucker_info_scm(int imol, SCM residue_spec, int do_pukka_pucker_check);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return False if residue not found,  otherwise
 [[phosphate_distance, puckered_atom, out_of_plane_distance, plane_distortion], chain_id, resno, ins_code]

 (where plane_distortion is for the other 4 atoms in the plane (I think)).
 
 and if there is no following residue, then the phosphate distance
 cannot be calculated, so the (inner) list is null (not filled).
*/
PyObject *pucker_info_py(int imol, PyObject *residue_spec, int do_pukka_pucker_check);
#endif /* USE_PYTHON */
#endif /*  __cplusplus */

/*! \brief Return a molecule that contains a residue that is the WC pair
   partner of the clicked/picked/selected residue */
int watson_crick_pair(int imol, const char * chain_id, int resno);

/* not for user level */
void setup_base_pairing(int state);


/* \} */

/*  ----------------------------------------------------------------------- */
/*                  sequence (assignment)                                   */
/*  ----------------------------------------------------------------------- */
/* section Sequence (Assignment) */
/*! \name Sequence (Assignment) */
/* \{ */

/*! \brief Print the sequence to the console of the given molecule */
void print_sequence_chain(int imol, const char *chain_id);
/*! \brief Assign a FASTA sequence to a given chain in the  molecule */
void assign_fasta_sequence(int imol, const char *chain_id_in, const char *seq);
/*! \brief Assign a PIR sequence to a given chain in the molecule.  If
  the chain of the molecule already had a chain assigned to it, then
  this will overwrite that old assignment with the new one. */
void assign_pir_sequence(int imol, const char *chain_id_in, const char *seq);
/* I don't know what this does. */
void assign_sequence(int imol_model, int imol_map, const char *chain_id);
/*! \brief Assign a sequence to a given molecule from (whatever) sequence
  file. */
void assign_sequence_from_file(int imol, const char *file);
/*! \brief Assign a sequence to a given molecule from a simple string */
void assign_sequence_from_string(int imol, const char *chain_id_in, const char *seq);
/*! \brief Delete all the sequences from a given molecule */
void delete_all_sequences_from_molecule(int imol);
/*! \brief Delete the sequence for a given chain_id from a given molecule */
void delete_sequence_by_chain_id(int imol, const char *chain_id_in);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief return the sequence info that has been assigned to molecule
  number imol. return as a list of dotted pairs (list (cons chain-id
  seq)).  To be used in constructing the cootaneer gui.  Return Scheme
  False when no sequence has been assigned to imol. */
SCM sequence_info(int imol);

/*! \brief do a internal alignment of all the assigned sequences,
  return a list of mismatches that need to be made to model number
  imol to match the input sequence.

Return a list of mutations deletions insetions.
Return scheme false on failure to align (e.g. not assigned sequence)
and the empty list on no alignment mismatches.*/
SCM alignment_mismatches_scm(int imol);
#endif /* USE_GUILE */

#ifdef USE_PYTHON
/*! \brief return the sequence info that has been assigned to molecule
  number imol. return as a list of dotted pairs [[chain-id, seq]].  To
  be used in constructing the cootaneer gui.  Return False when no
  sequence has been assigned. */
PyObject *sequence_info_py(int imol);
/*! \brief 

  do a internal alignment of all the assigned sequences,
  return a list of mismatches that need to be made to model number
  imol to match the input sequence.

Return a list of mutations deletions insetions.
Return  False on failure to align (e.g. not assigned sequence)
and the empty list on no alignment mismatches.*/
PyObject *alignment_mismatches_py(int imol);
#endif /* USE_PYTHON */
#endif /* C++ */
/* \} */
 
/*  ----------------------------------------------------------------------- */
/*                  Surfaces                                                */
/*  ----------------------------------------------------------------------- */
/* section Surfaces */
/*! \name Surfaces */
/* \{ */
/*! \brief draw surface of molecule number imol 

if state = 1 draw the surface (normal representation goes away)

if state = 0 don't draw surface */
void do_surface(int imol, int istate);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  FFfearing                                               */
/*  ----------------------------------------------------------------------- */
/*! \name FFFearing */
/* \{ */
/*! \brief fffear search model in molecule number imol_model in map
   number imol_map */
int fffear_search(int imol_model, int imol_map);
/*! \brief set and return the fffear angular resolution in degrees */
void set_fffear_angular_resolution(float f); 
/*! \brief return the fffear angular resolution in degrees */
float fffear_angular_resolution();
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  remote control                                          */
/*  ----------------------------------------------------------------------- */
/* section Remote Control */
/*! \name Remote Control */
/* \{ */
/*! \brief try to make socket listener */
void make_socket_listener_maybe(); 
int coot_socket_listener_idle_func(GtkWidget *w);
void set_coot_listener_socket_state_internal(int sock_state);

void set_socket_string_waiting(const char *s);

void set_remote_control_port(int port_number);
int get_remote_control_port_number();


/* tooltip */
void set_tip_of_the_day_flag(int state);

/* \} */
/*  ----------------------------------------------------------------------- */
/*                  Display lists                                           */
/*  ----------------------------------------------------------------------- */
/* section Display Lists for Maps */
/*! \brief Should display lists be used for maps? It may speed things
  up if these are turned on (or off) - depends on graphics card and
  drivers.  Pass 1 for on, 0 for off. */
void set_display_lists_for_maps(int i);

/*! \brief return the state of display_lists_for_maps.   */
int display_lists_for_maps_state();

/* update the maps to the current position - rarely needed */
void update_maps();

/*  ----------------------------------------------------------------------- */
/*                  Preferences Notebook                                    */
/*  ----------------------------------------------------------------------- */
/* section Preferences */
void preferences();
void show_preferences();
void clear_preferences();
void set_mark_cis_peptides_as_bad(int istate);
int show_mark_cis_peptides_as_bad_state();
#if (GTK_MAJOR_VERSION > 1)
void show_hide_preferences_tabs(GtkToggleToolButton *toggletoolbutton, int preference_type);
#endif // GTK_MAJOR_VERSION
void update_preference_gui();
void make_preferences_internal();
void make_preferences_internal_default();
void reset_preferences();
void save_preferences();
void preferences_internal_change_value_int(int preference_type, int ivalue);
void preferences_internal_change_value_int2(int preference_type, int ivalue1, int ivalue2);
void preferences_internal_change_value_float(int preference_type, float fvalue);
void preferences_internal_change_value_float3(int preference_type, 
					float fvalue1, float fvalue2, float fvalue3);
//void preferences_internal_change_value_vector_add_remove(int preference_type, int ivalue, int add_remove_flag);
void show_model_toolbar_icon(int pos);
void hide_model_toolbar_icon(int pos);
void fill_preferences_model_toolbar_icons(GtkWidget *preferences,
				     	  GtkWidget *scrolled_window);
void update_model_toolbar_icons_menu();

/*  ----------------------------------------------------------------------- */
/*                  Browser Help                                            */
/*  ----------------------------------------------------------------------- */
/*! \name Browser Interface */
/* \{ */
/*! \brief try to open given url in Web browser */
void browser_url(const char *url);
/*! \brief set command to open the web browser, 

examples are "open" or "mozilla" */
void set_browser_interface(const char *browser);

/*! \brief the search interface

find words, construct a url and open it. */
void handle_online_coot_search_request(const char *entry_text);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Generic Objects                                         */
/*  ----------------------------------------------------------------------- */
/*! \name Generic Objects */
/* \{ */

/*! \brief create a new generic object with name objname and return the index 
   of the object */
int new_generic_object_number(const char *objname);

/*! \brief add line to generic object object_number */
void to_generic_object_add_line(int object_number, 
				const char *colour,
				int line_width,
				float from_x1, 
				float from_y1, 
				float from_z1, 
				float to_x2, 
				float to_y2, 
				float to_z2);

/*! \brief add a dashed line to generic object object_number 

dash_density is number of dashes per Angstrom.*/
void to_generic_object_add_dashed_line(int object_number, 
				       const char *colour,
				       int line_width,
				       float dash_density,
				       float from_x1, 
				       float from_y1, 
				       float from_z1, 
				       float to_x2, 
				       float to_y2, 
				       float to_z2); 

/*! \brief add point to generic object object_number */
void to_generic_object_add_point(int object_number, 
				 const char *colour,
				 int point_width,
				 float from_x1, 
				 float from_y1, 
				 float from_z1);

/*! \brief add a display list handle generic object */
void to_generic_object_add_display_list_handle(int object_number, int display_list_id); 

/*! \brief set the display status of object number object_number, 

  when they are created, by default objects are not displayed, so we
  generally need this function.  */
void set_display_generic_object(int object_number, short int istate);

/*! \brief is generic display object displayed?

  @return 1 for yes, otherwise 0  */
int generic_object_is_displayed_p(int object_number);

/*! \brief return the index of the object with name name, if not, return -1; */
int generic_object_index(const char *name);

/*! \brief what is the name of generic object number obj_number? 

 @return 0 (NULL) (scheme False)  on obj_number not available */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM generic_object_name_scm(int obj_number);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *generic_object_name_py(int obj_number);
#endif /* USE_PYTHON */
#endif /*  __cplusplus */

/*! \brief return the number of generic display objects */
int number_of_generic_objects();

/*! \brief print to the console the name and display status of the
  generic display objects */
void generic_object_info(); 

/*! \brief does generic display object number obj_no have things to
  display? (predicate name)

@return 0 for no things, 1 for things. */
short int generic_object_has_objects_p(int obj_no); 

/*! \brief close generic object, clear the lines/points etc, not
  available for buttons/displaying etc */
void close_generic_object(int object_number);

/*! \brief has the generic object been closed? 

   @return 1 for yes, 0 othersize
*/
short int is_closed_generic_object_p(int object_number);

/*! \brief clear out the lines and points from object_number, but keep
  it displayable (not closed). */
void generic_object_clear(int object_number);

/*! \brief a kludgey thing, so that the generic objects gui can be
  called from a callback.  */
void generic_objects_gui_wrapper();



/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Molprobity interface                                    */
/*  ----------------------------------------------------------------------- */
/*! \name Molprobity Interface */
/* \{ */
/*! \brief pass a filename that contains molprobity's probe output in XtalView 
format */
void handle_read_draw_probe_dots(const char *dots_file);

/*! \brief pass a filename that contains molprobity's probe output in unformatted
format */
void handle_read_draw_probe_dots_unformatted(const char *dots_file, int imol, int show_clash_gui_flag);


/*! \brief shall we run molprobity for on edit chi angles intermediate atoms? */
void set_do_probe_dots_on_rotamers_and_chis(short int state);
/*! \brief return the state of if run molprobity for on edit chi
  angles intermediate atoms? */
short int do_probe_dots_on_rotamers_and_chis_state();
/*! \brief shall we run molprobity after a refinement has happened? */
void set_do_probe_dots_post_refine(short int state);
/*! \brief show the state of shall we run molprobity after a
  refinement has happened? */
short int do_probe_dots_post_refine_state();

/*! \brief make an attempt to convert pdb hydrogen name to the name
  used in Coot (and the refmac dictionary, perhaps). */
char *unmangle_hydrogen_name(const char *pdb_hydrogen_name);

/*! \brief set the radius over which we can run interactive probe,
  bigger is better but slower. 

  default is 6.0 */
void set_interactive_probe_dots_molprobity_radius(float r);

/*! \brief return the radius over which we can run interactive probe.
*/
float interactive_probe_dots_molprobity_radius();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the parsed user mod fields from the PDB file
  file_name (output by reduce most likely) */
SCM user_mods_scm(const char *file_name);
#endif // USE_GUILE
#ifdef USE_PYTHON
/*! \brief return the parsed user mod fields from the PDB file
  file_name (output by reduce most likely) */
PyObject *user_mods_py(const char *file_name);
#endif // USE_PYTHON
#endif	/* c++ */

/* \} */


/*  ----------------------------------------------------------------------- */
/*           Sharpen                                                        */
/*  ----------------------------------------------------------------------- */
/*! \name Map Sharpening Interface */
/* \{ */
/*! \brief Sharpen map imol by b_factor (note (of course) that positive numbers 
    blur the map).  */
void sharpen(int imol, float b_factor);

GtkWidget *wrapped_create_map_shapening_dialog();
void map_sharpening_map_select(GtkWidget *item, GtkPositionType pos);
void map_sharpening_value_changed (GtkAdjustment *adj, GtkWidget *window);
int fill_option_menu_with_map_options(GtkWidget *option_menu, GtkSignalFunc signalfunc);
/*! \brief set the limit of the b-factor map sharpening slider (default 30) */
void set_map_sharpening_scale_limit(float f);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*           Intermediate Atom Manipulation                                 */
/*  ----------------------------------------------------------------------- */

/*! \name Intermediate Atom Manipulation Interface */
/* \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM drag_intermediate_atom_scm(SCM atom_spec, SCM position);
#endif 
#ifdef USE_PYTHON
PyObject *drag_intermediate_atom_py(PyObject *atom_spec, PyObject *position);
#endif 
#endif /* c++ */
/* \} */

/*  ----------------------------------------------------------------------- */
/*           Fixed Atom Manipulation                                        */
/*  ----------------------------------------------------------------------- */

/*! \name Marking Fixed Atom Interface */
/* \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM mark_atom_as_fixed_scm(int imol, SCM atom_spec, int state);
#endif 
#ifdef USE_PYTHON
PyObject *mark_atom_as_fixed_py(int imol, PyObject *atom_spec, int state);
#endif 
#endif /* c++ */

void setup_fixed_atom_pick(short int ipick, short int is_unpick);

/*! \brief clear all fixed atoms */
void clear_all_fixed_atoms(int imol);
void store_fixed_atom_dialog(GtkWidget *w);
GtkWidget *wrapped_create_fixed_atom_dialog();
void clear_fixed_atoms_all();

/* produce debugging output from problematic atom picking  */
void set_debug_atom_picking(int istate);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Partial Charge                                          */
/*  ----------------------------------------------------------------------- */
/*! \name Partial Charges */
/* \{ */
/*! \brief show the partial charges for the residue of the given specs
   (charges are read from the dictionary) */
void show_partial_charge_info(int imol, const char *chain_id, int resno, const char *ins_code); 
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  EM Interface                                            */
/*  ----------------------------------------------------------------------- */
/*! \name EM interface */
/* \{ */
/*! \brief Scale the cell, for use with EM maps, where the cell needs
   to be adjusted.  Use like:  (scale-cell 2 1.012 1.012 1.012). Return error 
   status, 1 means it worked, 0 means it did not work. */
int scale_cell(int imol_map, float fac_u, float fac_v, float fac_w); 
/* \} */


/*  ----------------------------------------------------------------------- */
/*                  CCP4i Interface                                         */
/*  ----------------------------------------------------------------------- */
/*! \name CCP4mg Interface */
/* \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return a list of pairs of strings, the project names and
  the directory.  Include aliases. */
SCM ccp4i_projects_scm();
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return a list of pairs of strings, the project names and
  the directory.  Include aliases. */
PyObject *ccp4i_projects_py();
#endif /* USE_PYTHON */
#endif /* c++ */

/*! \brief write a ccp4mg picture description file */
void write_ccp4mg_picture_description(const char *filename);
/*! \brief get element colour for imol as Python formatted list char*/
char *get_atom_colour_from_mol_no(int imol, const char *element);

/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Dipoles                                                 */
/*  ----------------------------------------------------------------------- */
/*! \name Dipoles */
/* \{ */
void delete_dipole(int imol, int dipole_number);
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief generate a dipole from all atoms in the given
  residues. Return the dipole description */
SCM add_dipole_for_residues_scm(int imol, SCM residue_specs);
/*! \brief return the dipole number */
SCM add_dipole_scm(int imol, const char* chain_id, int res_no, const char *ins_code); 
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief generate a dipole from all atoms in the given residues. */
PyObject *add_dipole_py(int imol, const char* chain_id, int res_no, 
			const char *ins_code); 
/*! \brief add a dipole given a set of residues.  Return a dipole
  description. */
PyObject *add_dipole_for_residues_py(int imol, PyObject *residue_specs);
#endif /* USE_PYTHON */
#endif /* c++ */
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Laplacian                                               */
/*  ----------------------------------------------------------------------- */
/*! \name Aux functions */
/* \{ */
/*! \brief Create the "Laplacian" (-ve second derivative) of the given map. */
int laplacian (int imol);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  Tips                                                    */
/*  ----------------------------------------------------------------------- */
/*! \name Tips Interface */
/* \{ */
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  PKGDATADIR                                              */
/*  ----------------------------------------------------------------------- */
/*! \name PKGDATADIR */
/* \{ */
#ifdef __cplusplus
#ifdef USE_PYTHON
PyObject *get_pkgdatadir_py();
#endif /* USE_PYTHON */
#ifdef USE_GUILE
SCM get_pkgdatadir_scm();
#endif
#endif /*  __cplusplus */
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  SMILES                                                  */
/*  ----------------------------------------------------------------------- */
/*! \name SMILES */
/* \{ */
/*! \brief display the SMILES string dialog */
void do_smiles_gui();
/* \} */
/*  ----------------------------------------------------------------------- */
/*                  Fun                                                     */
/*  ----------------------------------------------------------------------- */
/* section Fun */
void do_tw();

/*  ----------------------------------------------------------------------- */
/*                  Phenix Support                                          */
/*  ----------------------------------------------------------------------- */
/*! \name PHENIX Support */
/* \{ */
/*! \brief set the button label of the external Refinement program */
void set_button_label_for_external_refinement(const char *button_label);
/* \} */


/*  ----------------------------------------------------------------------- */
/*                  Text                                                    */
/*  ----------------------------------------------------------------------- */
/*! \name Graphics Text */
/* \{ */
/*! \brief Put text at x,y,z  

@return a text handle

size variable is currently ignored.*/
int place_text(const char*text, float x, float y, float z, int size);

/*! \brief Remove text */
void remove_text(int text_handle);
/* \} */

/*  ----------------------------------------------------------------------- */
/*                  PISA Interface                                      */
/*  ----------------------------------------------------------------------- */
/*! \name PISA Interaction */
/* \{ */
/*! \brief return the molecule number of the interacting
  residues. Return -1 if no new model was created.  */
int pisa_interaction(int imol_1, int imol_2);
/* \} */


/*  ----------------------------------------------------------------------- */
/*                  Jiggle fit                                              */
/*  ----------------------------------------------------------------------- */
/*! \name Jiggle Fit */
/* \{ */
/*!  \brief jiggle fit to the current refinment map.  return < -100 if
  not possible, else return the new best fit for this residue.  */
float fit_to_map_by_random_jiggle(int imol, const char *chain_id, int resno, const char *ins_code,
				  int n_trials,
				  float jiggle_scale_factor);
/* \} */


/*  ----------------------------------------------------------------------- */
/*                  SBase interface                                         */
/*  ----------------------------------------------------------------------- */
/*! \name SBase interface */
/* \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return a list of compoundIDs of in SBase of which the
  given string is a substring of the compound name */
SCM matching_compound_names_from_sbase_scm(const char *compound_name_fragment);
#endif
#endif

/*! \brief return the new molecule number of the monomer.

The monomer will have chainid "A" and residue number 1.

Return -1 on failure to get monomer. */
int get_sbase_monomer(const char *comp_id);
/* \} */


/*  ----------------------------------------------------------------------- */
/* Multirefine interface (because in guile-gtk there is no way to
		    insert toolbuttons into the toolbar) so this
		    rather kludgy interface.  It should go when we
		    move to guile-gnome, I think. */
/*  ----------------------------------------------------------------------- */
void toolbar_multi_refine_stop();
void toolbar_multi_refine_continue();
void toolbar_multi_refine_cancel();
void set_visible_toolbar_multi_refine_stop_button(short int state);
void set_visible_toolbar_multi_refine_continue_button(short int state);
void set_visible_toolbar_multi_refine_cancel_button(short int state);
/* button_type is one of "stop", "continue", "cancel"
   state is 1 for on, 0 for off. */
void toolbar_multi_refine_button_set_sensitive(const char *button_type, short int state);


/*  ----------------------------------------------------------------------- */
/*                  update self                                             */
/*  ----------------------------------------------------------------------- */
/* this function is here because it is called by c_inner_main() (ie. need a c interface). */
void run_update_self_maybe(); // called when --update-self given at command line

/*  ----------------------------------------------------------------------- */
/*                  experimental                                            */
/*  ----------------------------------------------------------------------- */
void nsv(int imol);
void sequence_view_old_style(int imol);
#ifdef __cplusplus
#ifdef USE_GUILE
void user_defined_click_scm(int n_clicks, SCM func);
#endif
#ifdef USE_PYTHON
void user_defined_click_py(int n_clicks, PyObject *func);
#endif /* PYTHON */
#endif /* c++ */

#endif /* C_INTERFACE_H */
END_C_DECLS

