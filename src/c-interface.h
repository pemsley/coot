/* src/c-interface.h
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 by The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
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

/* svn $Id: c-interface.h 1458 2007-01-26 20:20:18Z emsley $ */

/*! \file
  \brief Coot Scripting Interface - General

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

// Python conditionally compiled test is needed for WebAssembly build
#ifdef USE_PYTHON
#include "Python.h"
#endif

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
#include <cstdio> /* for std::FILE in gmp.h for libguile.h */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#include <libguile.h>
#pragma GCC diagnostic pop
#else
#include <string> /* for std::string; included (sic!) in above for guile */
#endif /*  USE_GUILE */
#endif /* c++ */

#ifdef USE_PYTHON
#include "Python.h"
#endif

#include <gtk/gtk.h>

#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS


#define COOT_SCHEME_DIR "COOT_SCHEME_DIR"
#define COOT_PYTHON_DIR "COOT_PYTHON_DIR"

/*  ------------------------------------------------------------------------ */
/*                         Startup Functions:                                */
/*  ------------------------------------------------------------------------ */
#ifdef USE_GUILE
void try_load_scheme_extras_dir();
#endif /* USE_GUILE */
#ifdef USE_PYTHON

void try_load_python_extras_dir();
#endif /* USE_PYTHON */

/* section Startup Functions */
/*!  \name Startup Functions */
/*! \{ */
/*!  \brief tell coot that you prefer to run python scripts if/when
  there is an option to do so. */
void set_prefer_python();

/*! \brief the python-prefered mode.

This is available so that the scripting functions know whether on not
to put themselves onto in as menu items.

If you consider using this, consider in preference use_gui_qm == 2,
which is used elsewhere to stop python functions adding to the gui,
when guile-gtk functions have alread done so.  We should clean up this
(rather obscure) interface at some stage.

@return 1 for python is prefered, 0 for not. */
int prefer_python();

/*! \} */

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
/*! \{ */

/*! \brief Show Paths in Display Manager?

    Some people don't like to see the full path names in the display manager
    here is the way to turn them off, with an argument of 1.
*/
void set_show_paths_in_display_manager(int i);

/*! \brief return the internal state

   What is the internal flag?

   @return 1 for "yes, display paths", 0 for not
 */
int show_paths_in_display_manager_state();

/*! \brief add an extension to be treated as coordinate files
   @param ext the extension to be added
*/
void add_coordinates_glob_extension(const char *ext);

/*! \brief add an extension to be treated as data (reflection) files
   @param ext the extension to be added
*/
void add_data_glob_extension(const char *ext);

/*! \brief add an extension to be treated as geometry dictionary files
   @param ext the extension to be added
*/
void add_dictionary_glob_extension(const char *ext);

/*! \brief add an extension to be treated as geometry map files
*/
void add_map_glob_extension(const char *ext);

/*! \brief remove an extension to be treated as coordinate files
   @param ext the extension to be added
*/
void remove_coordinates_glob_extension(const char *ext);

/*! \brief remove an extension to be treated as data (reflection) files
   @param ext the extension to be removed
*/
void remove_data_glob_extension(const char *ext);

/*! \brief remove an extension to be treated as geometry dictionary files
   @param ext the extension to be removed
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

/*! \brief, return the state of the above variable */
int filter_fileselection_filenames_state();

/*! \brief is the given file name suitable to be read as coordinates? */
short int file_type_coords(const char *file_name);

/*! \brief display the open coordinates dialog */
void open_coords_dialog();

/*! \brief this flag set chooser as default for windows, otherwise use
  selector 0 is selector 1 is chooser */


/* --- CHECKME - do these need to be here? --- */
void set_file_chooser_selector(int istate);
int file_chooser_selector_state();
void set_file_chooser_overwrite(int istate);
int file_chooser_overwrite_state();

/*! \brief show the export map GUI */
void export_map_gui(short int export_map_fragment);

/*! \} */


/*! \name Widget Utilities */
/*! \{ */

/*! \brief set the main window title.

function added for Lothar Esser */
void set_main_window_title(const char *s);

/*! \} */

/*  -------------------------------------------------------------------- */
/*                   mtz and data handling utilities                     */
/*  -------------------------------------------------------------------- */
/* section MTZ and data handling utilities */
/*! \name  MTZ and data handling utilities */
/*! \{ */
/* We try as .phs and .cif files first */

/*! \brief given a filename, try to read it as a data file

   We try as .phs and .cif files first */
void manage_column_selector(const char *filename);

/*! \} */

/*  -------------------------------------------------------------------- */
/*                     Molecule Functions       :                        */
/*  -------------------------------------------------------------------- */
/* section Molecule Info Functions */
/*! \name Molecule Info Functions */
/*! \{ */

/*! \brief the number of residues in chain chain_id and molecule number imol
  @return the number of residues
*/
int chain_n_residues(const char *chain_id, int imol);
/*! \brief internal function for molecule centre

@return status, less than -9999 is for failure (eg. bad imol); */
float molecule_centre_internal(int imol, int iaxis);
/*! \brief a residue seqnum (normal residue number) from a residue
  serial number

   @return < -9999 on failure */
int  seqnum_from_serial_number(int imol, const char *chain_id,
			       int serial_num);

/*! \brief the insertion code of the residue.

   @return NULL (scheme False) on failure. */
char *insertion_code_from_serial_number(int imol, const char *chain_id, int serial_num);

#ifdef __cplusplus
#ifdef USE_PYTHON
PyObject *python_representation_kk(int imol);
#endif
#endif

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

/*! \brief return the number of models in molecule number imol

useful for NMR or other such multi-model molecules.

return the number of models or -1 if there was a problem with the
given molecule.
*/
int n_models(int imol);

/*! \brief get the number of chains in molecule number imol

  @param imol is the molecule index
  @return the number of chains
*/
int n_chains(int imol);

#ifdef USE_PYTHON
/*! \brief get the chain ids of molecule number imol

  @param imol is the molecule index
  @return a list of the the chain ids or False on failure
*/
PyObject *get_chain_ids_py(int imol);
#endif

/*! \brief is this a solvent chain? [Raw function]

   This is a raw interface function, you should generally not use
   this, but instead use (is-solvent-chain? imol chain-id)

   This wraps the mmdb function isSolventChain().

   @param imol is the molecule index
   @param chain_id is the chain id (e.g. "A" or "B")
   @return -1 on error, 0 for no, 1 for is "a solvent chain".  We
   wouldn't want to be doing rotamer searches and the like on such a
   chain.

 */
int is_solvent_chain_p(int imol, const char *chain_id);

/*! \brief is this a protein chain? [Raw function]

   This is a raw interface function, you should generally not use
   this, but instead use (is-protein-chain? imol chain-id)

   @return -1 on error, 0 for no, 1 for is "a protein chain".  We
   wouldn't want to be doing rotamer searches and the like on such a
   chain.

   This wraps the mmdb function isAminoacidChain().
 */
int is_protein_chain_p(int imol, const char *chain_id);

/*! \brief is this a nucleic acid chain? [Raw function]

   This is a raw interface function, you should generally not use
   this, but instead use (is-nucleicacid-chain? imol chain-id)

   @return -1 on error, 0 for no, 1 for is "a nucleicacid chain".  We
   wouldn't want to be doing rotamer searches and the like on such a
   chain.

   This wraps the mmdb function isNucleotideChain().
   For completeness.
 */
int is_nucleotide_chain_p(int imol, const char *chain_id);


/*! \brief return the number of residues in the molecule,

return -1 if this is a map or closed.
 */
int n_residues(int imol);

/*! \brief return the ATOMs of residues in the molecule,

return -1 if this is a map or closed. HETATMs are not counted.
 */
int n_atoms(int imol);


/* Does this work? */
/*! \brief return a list of the remarks of hte molecule number imol
  */
/* list remarks(int imol); */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM remarks_scm(int imol);
/* return a list or scheme false */
SCM residue_centre_scm(int imol, const char *chain_id, int resno, const char *ins_code);
#endif
#ifdef USE_PYTHON
/* return a list or python false */
PyObject *remarks_py(int imol);
PyObject *residue_centre_py(int imol, const char *chain_id, int resno, const char *ins_code);
#endif
#endif


#ifdef __cplusplus
#ifdef USE_GUILE
SCM model_composition_statistics_scm(int imol);
#endif

PyObject *model_composition_statistics_py(int imol);
#endif


/*! \brief sort the chain ids of the imol-th molecule in lexographical order */
void sort_chains(int imol);

/*! \brief sort the residues of the imol-th molecule */
void sort_residues(int imol);

/*! \brief a gui dialog showing remarks header info (for a model molecule). */
void remarks_dialog(int imol);

/*! \brief simply print secondary structure info to the
  terminal/console.  In future, this could/should return the info.  */
void print_header_secondary_structure_info(int imol);

/*! \brief get the secondary structure from the header
 *
 * @param imol the molecule index
 * @return a dictionary of header info
 * Returns: {'helices': [...], 'strands': [...]}
 *
 * Each helix dict contains:
 *  serNum, helixID, initChainID, initSeqNum, endChainID, endSeqNum, length, comment
 *
 * Each strand dict contains:
 *  SheetID, strandNo, initChainID, initSeqNum, endChainID, endSeqNum
 */
PyObject *get_header_secondary_structure_info(int imol);


/*! \brief add secondary structure info to the
  internal representation of the model */
void add_header_secondary_structure_info(int imol);


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

/*! \brief Copy a molecule with addition of a ligand and a deletion of
  current ligand.

  This function is used when adding a new (modified) ligand to a
  structure.  It creates a new molecule that is a copy of the current
  molecule except that the new ligand is added and the current
  ligand/residue is deleted.

 */
int add_ligand_delete_residue_copy_molecule(int imol_ligand_new,
					    const char *chain_id_ligand_new,
					    int resno_ligand_new,
					    int imol_current,
					    const char *chain_id_ligand_current,
					    int resno_ligand_current);

/*! \brief Experimental interface for Ribosome People.

Ribosome People have many chains in their pdb file, they prefer segids
to chainids (chainids are only 1 character).  But coot uses the
concept of chain ids and not seg-ids.  mmdb allow us to use more than
one char in the chainid, so after we read in a pdb, let's replace the
chain ids with the segids. Will that help? */
int exchange_chain_ids_for_seg_ids(int imol);

/*! \brief show the remarks browser */
void show_remarks_browswer();


/*! \} */

/*  -------------------------------------------------------------------- */
/*                     Library/Utility Functions:                        */
/*  -------------------------------------------------------------------- */

/* section Library and Utility Functions */
/*! \name Library and Utility Functions */
/*! \{ */

#ifdef __cplusplus

#ifdef USE_GUILE
SCM coot_sys_build_type_scm();
#endif
#ifdef USE_PYTHON
PyObject *coot_sys_build_type_py();
#endif /* USE_PYTHON */
#endif /* c++ */

/*! \brief return the git revision count for for this build.
  */
int git_revision_count();
/*! \brief an alias to git_revision_count() for backwards compatibility  */
int svn_revision();


/*! \brief return the name of molecule number imol

 @return 0 if not a valid name ( -> False in scheme)
 e.g. "/a/b/c.pdb" for "d/e/f.mtz FWT PHWT" */
const char *molecule_name(int imol);
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the molecule name without file extension */
SCM molecule_name_stub_scm(int imol, int include_path_flag);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the molecule name without file extension */
PyObject *molecule_name_stub_py(int imol, int include_path_flag);
#endif /* USE_PYTHON */
#endif	/* __cplusplus */
/*! \brief set the molecule name of the imol-th molecule */
void set_molecule_name(int imol, const char *new_name);
int coot_checked_exit(int retval);
/*! \brief exit from coot, give return value retval back to invoking
  process. */
void coot_real_exit(int retval);
/*! \brief exit without writing a state file */
void coot_no_state_real_exit(int retval);

/*! \brief exit coot doing clear-backup maybe */
void coot_clear_backup_or_real_exit(int retval);
/*! \brief exit coot, write a state file */
void coot_save_state_and_exit(int retval, int save_state_flag);


#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief run clear-backups */
void run_clear_backups(int retval);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
void run_clear_backups_py(int retval);
#endif /* USE_PYTHON */
#endif  /* c++ */


/*! \brief What is the molecule number of first coordinates molecule?

   return -1 when there is none. */
int first_coords_imol();

/*! \brief molecule number of first small (<400 atoms) molecule.

return -1 on no such molecule
  */
int first_small_coords_imol();

/*! \brief What is the molecule number of first unsaved coordinates molecule?

   return -1 when there is none. */
int first_unsaved_coords_imol();

/*! \brief convert the structure factors in cif_file_name to an mtz
  file.

  Return 1 on success. Return 0 on a file without Rfree, return
  -1 on complete failure to write a file. */
int mmcif_sfs_to_mtz(const char *cif_file_name, const char *mtz_file_name);

/*! \} */

/*  -------------------------------------------------------------------- */
/*                    More Library/Utility Functions:                    */
/*  -------------------------------------------------------------------- */
/* section Graphics Utility Functions */
/*! \name Graphics Utility Functions */
/*! \{ */

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

/*! \brief set the GUI dark mode state
 */
void set_use_dark_mode(short int state);

/*! \brief is the python interpreter at the prompt?

@return 1 for yes, 0 for no.*/
short int python_at_prompt_at_startup_state();

/*! \brief "Reset" the view


  return 1 if we moved, else return 0.

   centre on last-read molecule with zoom 100. If we are there, then
   go to the previous molecule, if we are there, then go to the origin. */
int reset_view();

/*! \brief set the view rotation scale factor

 Useful/necessary for high resolution displayed, where, without this factor
 the view doesn't rotate enough */
void set_view_rotation_scale_factor(float f);

/*! \brief return the number of molecules (coordinates molecules and
  map molecules combined) that are currently in coot

  @return the number of molecules (closed molecules are not counted) */
int get_number_of_molecules();

/*! \brief As above, return the number of molecules (coordinates molecules and
  map molecules combined) that are currently in coot.

  This is the old name for the function.

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

/*! \brief Spin spin spin (or not) */
void toggle_idle_spin_function();

/*! \brief Rock (not roll) (self-timed) */
void toggle_idle_rock_function();
/* used by above to set the angle to rotate to (time dependent) */
double get_idle_function_rock_target_angle();


/*! \brief Settings for the inevitable discontents who dislike the
   default rocking rates (defaults 1 and 1)  */
void set_rocking_factors(float width_scale, float frequency_scale);

/*! \brief how far should we rotate when (auto) spinning? Fast
  computer? set this to 0.1  */
void set_idle_function_rotate_angle(float f);  /* degrees */

/*! \brief what is the idle function rotation angle? */
float idle_function_rotate_angle();

/*! \brief make a model molecule from the give file name.

* If the file updates, then the model will be updated. */
int make_updating_model_molecule(const char *filename);

/* or better still, use the json file from refmac ... but not yet. */
/* void updating_refmac_refinement_files(const char *updating_refmac_refinement_files_json_file_name); */

/* used by above, no API for this - also, not yet */
/* int updating_refmac_refinement_json_timeout_function(gpointer data); */


/*! \brief show the updating maps gui

this function is called from callbacks.c and calls a python gui function
*/
void show_calculate_updating_maps_pythonic_gui();

/*! \brief enable reading PDB/pdbx files with duplicate sequence numbers */
void allow_duplicate_sequence_numbers();

/*! \brief shall we convert nucleotides to match the old dictionary
  names?

Usually (after 2006 or so) we do not want to do this (given current
Coot architecture).  Coot should handle the residue synonyms
transparently.

default off (0).

 */
void set_convert_to_v2_atom_names(short int state);


#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief - get the name state of the input model. Return false if there was an erro with the molecule index */
SCM get_input_model_was_cif_state_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief - get the name state of the input model */
PyObject *get_input_molecule_was_in_mmcif_state_py(int imol);
#endif
#endif


/*! \brief some programs produce PDB files with ATOMs where there
  should be HETATMs.  This is a function to assign HETATMs as per the
  PDB definition. */
int assign_hetatms(int imol);

/*! \brief if this is not a standard group, then turn the atoms to HETATMs.

Return 1 on atoms changes, 0 on not. Return -1 if residue not found.
*/
int hetify_residue(int imol, const char * chain_id, int resno, const char *ins_code);

/*! \brief residue has HETATMs?

return 1 if all atoms of the specified residue are HETATMs, else,
return 0.  If residue not found, return -1. */
int residue_has_hetatms(int imol, const char * chain_id, int resno, const char *ins_code);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief - get the specs for hetgroups - waters are not counted as het-groups. */
SCM het_group_residues_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief - get the specs for hetgroups - waters are not counted as het-groups. */
PyObject *het_group_residues_py(int imol);
#endif
#endif

/*! \brief return the number of non-hydrogen atoms in the given
  het-group (comp-id).

Return -1 on comp-id not found in dictionary.  */
int het_group_n_atoms(const char *comp_id);

/*! \brief replace the parts of molecule number imol that are
  duplicated in molecule number imol_frag */
int replace_fragment(int imol_target, int imol_fragment, const char *atom_selection);

/*! \brief copy the given residue range from the reference chain to the target chain

resno_range_start and resno_range_end are inclusive. */
int copy_residue_range(int imol_target,    const char *chain_id_target,
		       int imol_reference, const char *chain_id_reference,
		       int resno_range_start, int resno_range_end);

/*! \brief replace the given residues from the reference molecule to the target molecule
*/
#ifdef __cplusplus
#ifdef USE_GUILE
int replace_residues_from_mol_scm(int imol_target,
				 int imol_ref,
				 SCM residue_specs_list_ref_scm);
#endif /* USE_GUILE */

#ifdef USE_PYTHON
int replace_residues_from_mol_py(int imol_target,
				 int imol_ref,
				 PyObject *residue_specs_list_ref_py);
#endif /* USE_PYTHON */
#endif	/* __cplusplus */


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

/*! \brief give a warning dialog if density it too dark (blue) */
void check_for_dark_blue_density();

/* is this a good place for this function? */

/*! \brief sets the density map of the given molecule to be drawn as a
  (transparent) solid surface. */
void set_draw_solid_density_surface(int imol, short int state);

/*! \brief toggle for standard lines representation of map.

  This turns off/on standard lines representation of map.  transparent
  surface is another representation type.

  If you want to just turn off a map, don't use this, use
  set_map_displayed().

  */
void set_draw_map_standard_lines(int imol, short int state);

/*! \brief set the opacity of density surface representation of the
  given map.

0.0 is totally transparent, 1.0 is completely opaque and (because the
objects are no longer depth sorted) considerably faster to render.
0.3 is a reasonable number.

 */
void set_solid_density_surface_opacity(int imol, float opacity);

float get_solid_density_surface_opacity(int imol);

/*! \brief set the flag to do flat shading rather than smooth shading
  for solid density surface.

Default is 1 (on. */
void set_flat_shading_for_solid_density_surface(short int state);

/*! \} */

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

/*! \brief turn off (0) or on (1) auto recontouring (on screen centre change) (default it on) */
void  set_auto_recontour_map(int state);

/*! \brief return the auto-recontour state */
int get_auto_recontour_map();

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
						     const char *next_ncs_chain,
						     short int forward_flag);

/*! \brief set show frame-per-second flag */
void set_show_fps(int t);

/*! \brief the old name for set_show_fps() */
void set_fps_flag(int t);
/*! \brief set the state of show frames-per-second flag */
int  get_fps_flag();

/*! \brief set show frame-per-second flag */
void set_show_fps(int t);

/*! \brief set a flag: is the origin marker to be shown? 1 for yes, 0
  for no. */
void set_show_origin_marker(int istate);
/*! \brief return the origin marker shown? state */
int  show_origin_marker_state();

/*! \brief hide the horizontal main toolbar in the GTK2 version */
void hide_main_toolbar();
/*! \brief show the horizontal main toolbar in the GTK2 version
  (the toolbar is shown by default) */
void show_main_toolbar();

/*! \brief reparent the Model/Fit/Refine dialog so that it becomes
  part of the main window, next to the GL graphics context */
int suck_model_fit_dialog();
int suck_model_fit_dialog_bl();

/*! \brief model-fit-refine dialog stays on top */
void set_model_fit_refine_dialog_stays_on_top(int istate);
/*! \brief return the state model-fit-refine dialog stays on top */
int model_fit_refine_dialog_stays_on_top_state();



/*! \} */

/*  ----------------------------------------------------------------------- */
/*                           mouse buttons                                  */
/*  ----------------------------------------------------------------------- */
/* section Mouse Buttons */
/*! \name   Mouse Buttons */
/*! \{ */

/*! \brief quanta-like buttons

Note, when you have set these, there is no way to turn them of
   again (other than restarting). */
void quanta_buttons();
/*! \brief quanta-like zoom buttons

Note, when you have set these, there is no way to turn them of
   again (other than restarting). */
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

Not useful for MCP. For interactive use only.

   @return 1 if successfully found a blob and moved there.
   return 0 if no move.
*/
int blob_under_pointer_to_screen_centre();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return scheme false or a list of molecule number and an atom spec  */
SCM select_atom_under_pointer_scm();
#endif

#ifdef USE_PYTHON
/*! \brief return Python false or a list of molecule number and an atom spec
 *
 * Not useful for MCP. For interactive use only.
*/
PyObject *select_atom_under_pointer_py();
#endif
#endif /* __cplusplus */

/*! \} */

/*  --------------------------------------------------------------------- */
/*                      Cursor Functions:                                 */
/*  --------------------------------------------------------------------- */
/* section Cursor Function */
/*! \name Cursor Function */
/*! \{ */
/*! \brief normal cursor */
void normal_cursor();
/*! \brief fleur cursor */
void fleur_cursor();
/*! \brief pick cursor maybe */
void pick_cursor_maybe();
/*! \brief rotate cursor */
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
/*! \{ */

/*! \brief display the Display Manager dialog */
void show_select_map_frame();
/*! \brief Allow the changing of Model/Fit/Refine button label from
  "Rotate/Translate Zone" */
void set_model_fit_refine_rotate_translate_zone_label(const char *txt);
/*! \brief Allow the changing of Model/Fit/Refine button label from
  "Place Atom at Pointer" */
void set_model_fit_refine_place_atom_at_pointer_label(const char *txt);


/*! \brief shall atoms with zero occupancy be moved when refining? (default 1, yes) */
void set_refinement_move_atoms_with_zero_occupancy(int state);
/*! \brief return the state of "shall atoms with zero occupancy be moved
  when refining?" */
int refinement_move_atoms_with_zero_occupancy_state();

/*! \} */

/*  --------------------------------------------------------------------- */
/*                      backup/undo functions:                            */
/*  --------------------------------------------------------------------- */
/* section Backup Functions */
/*! \name Backup Functions */
/*! \{ */
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

/*! \brief apply undo - the "Undo" button callback
 *
 * undo the most recent modification on the model
 * set in set_undo_molecule().
 *
 * @return 1 on succesful undo, 0 on failed to undo.
 */
int apply_undo();		/* "Undo" button callback */

/*! \brief apply redo - the "Redo" button callback */
int apply_redo();

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

/*! \brief set the state for adding paths to backup file names

  by default directories names are added into the filename for backup
  (with / to _ mapping).  call this with state=1 to turn off directory
  names  */
void set_unpathed_backup_file_names(int state);
/*! \brief return the state for adding paths to backup file names*/
int  unpathed_backup_file_names_state();

/*! \brief set the state for adding paths to backup file names

  by default directories names are added into the filename for backup
  (with / to _ mapping).  call this with state=1 to turn off directory
  names  */
void set_decoloned_backup_file_names(int state);
/*! \brief return the state for adding paths to backup file names*/
int  decoloned_backup_file_names_state();


/*! \brief return the state for compression of backup files*/
int  backup_compress_files_state();

/*! \brief set if backup files will be compressed or not using gzip */
void  set_backup_compress_files(int state);

/*! \brief Make a backup for a model molecule
 *
 * @param imol the model molecule index
 * @description a description that goes along with this back point
 * @return the index of the backup, or -1 on failure
 */
int make_backup_checkpoint(int imol, const char *description);

/*! \brief Restore molecule from backup
 * 
 * restore model @p imol to checkpoint backup @p backup_index
 *
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 * @return the index of the backup, or -1 on failure
 */
int restore_to_backup_checkpoint(int imol, int backup_index);

#ifdef USE_PYTHON
/*! \brief Compare current model to backup
 * 
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 * @return a Python dict, with 2 items, a "status" which is either "ok" 
 *         or "error" or "bad-index". The other key is "moved-residues-list",
 *         the value for which is a list of residue specs for residues
 *         that have at least one atom in a different place (which might be empty).
 */
PyObject *compare_current_model_to_backup(int imol, int backup_index);
#endif

/*! \brief Print the history info
 * 
 */
void print_backup_history_info(int imol);

#ifdef USE_PYTHON
/*! \brief Get backup info
 * 
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 * @return a Python list of the given description (str)
 *         and a timestamp (str).
 */
PyObject *get_backup_info(int imol, int backup_index);
#endif

/*! \} */

/*  --------------------------------------------------------------------- */
/*                         recover session:                               */
/*  --------------------------------------------------------------------- */
/* section Recover Session Function */
/*! \name  Recover Session Function */
/*! \{ */
/*! \brief recover session

   After a crash, we provide this convenient interface to restore the
   session.  It runs through all the molecules with models and looks
   at the coot backup directory looking for related backup files that
   are more recent that the read file. (Not very good, because you
   need to remember which files you read in before the crash - should
   be improved.) */
void recover_session();
/*! \} */

/*  ---------------------------------------------------------------------- */
/*                       map functions:                                    */
/*  ---------------------------------------------------------------------- */
/* section Map Functions */
/*! \name  Map Functions */
/*! \{ */

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

#ifdef USE_PYTHON
/*! \brief Calculate structure factors and make a 2FoFC map and a Fo-Fc map updating the given
   molecule numbers for those maps - if thase molecule ids are not valid maps, them generate
   new maps (return the model number information in the returned object) */
PyObject *calculate_maps_and_stats_py(int imol_model,
                                      int imol_map_with_data_attached,
                                      int imol_map_2fofc,
                                      int imol_map_fofc);
#endif

/*! \brief Calculate structure factors from the model and update the given difference
           map accordingly */
void sfcalc_genmap(int imol_model, int imol_map_with_data_attached, int imol_updating_difference_map);

/*! \brief As above, calculate structure factors from the model and update the given difference
           map accordingly - but difference map gets updated automatically on modification of
           the imol_model molecule */
void set_auto_updating_sfcalc_genmap(int imol_model, int imol_map_with_data_attached, int imol_updating_difference_map);

/*! \brief As above, calculate structure factors from the model and update the given difference
           map accordingly - but the 2fofc and difference map get updated automatically on modification of
           the imol_model molecule */
void set_auto_updating_sfcalc_genmaps(int imol_model, int imol_map_with_data_attached, int imol_updating_2fofc_map, int imol_updating_difference_map);


/* gdouble* get_map_colour(int imol); delete on merge 20220228-PE */

#ifdef __cplusplus
#ifdef USE_GUILE
SCM get_map_colour_scm(int imol);
#endif
#ifdef USE_PYTHON
PyObject *get_map_colour_py(int imol);
#endif
#endif



/*! \brief set the map that is moved by changing the scroll wheel and
  change_contour_level(). */
void set_scroll_wheel_map(int imap);
/*! \brief return the molecule number to which the mouse scroll wheel
  is attached */
/*! \brief set the map that has its contour level changed by the
  scrolling the mouse wheel to molecule number imol (same as set_scroll_wheel_map()). */
void set_scrollable_map(int imol);
/*! \brief the contouring of which map is altered when the scroll wheel changes? */
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

/*! \brief set the colour of the imolth map using a (7-character) hex colour */
void set_map_hexcolour(int imol, const char *hex_colour);

/*! \brief set the contour level, direct control */
void set_contour_level_absolute(int imol_map, float level);
/*! \brief set the contour level, direct control in r.m.s.d. (if you like that sort of thing) */
void set_contour_level_in_sigma(int imol_map, float level);

/*! \brief get the contour level */
float get_contour_level_absolute(int imol);

/*! \brief get the contour level in rmd above 0. */
float get_contour_level_in_sigma(int imol);

/*! \brief set the sigma step of the last map to f sigma */
void set_last_map_sigma_step(float f);
/*! \brief set the contour level step

   set the contour level step of molecule number imol to f and
   variable state (setting state to 0 turns off contouring by sigma
   level)  */
void set_contour_by_sigma_step_by_mol(int imol, float f, short int state);

/*! \brief return the resolution of the data for molecule number imol.
   Return negative number on error, otherwise resolution in A (eg. 2.0) */
float data_resolution(int imol);

/*! \brief return the resolution set in the header of the
  model/coordinates file.  If this number is not available, return a
  number less than 0.  */
float model_resolution(int imol);

/*! \brief export (write to disk) the map of molecule number imol to
  filename.

  Return 0 on failure, 1 on success. */
int export_map(int imol, const char *filename);
/*! \brief export a fragment of the map about (x,y,z)  */
int export_map_fragment(int imol, float x, float y, float z, float radius, const char *filename);

/*! convenience function, called from callbacks.c */
void export_map_fragment_with_text_radius(int imol, const char *radius_text, const char *filename);

/*! \brief export a fragment of the map about (x,y,z)  */
int export_map_fragment_with_origin_shift(int imol, float x, float y, float z, float radius, const char *filename);

/*! \brief tmp interface for Hamish */
int export_map_fragment_to_plain_file(int imol, float x, float y, float z, float radius, const char *filename);

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


/*! \brief make a difference map, taking map_scale * imap2 from imap1,
  on the grid of imap1.  Return the new molecule number.
  Return -1 on failure. */
int difference_map(int imol1, int imol2, float map_scale);

/*! \brief by default, maps that are P1 and have 90 degree angles
           are considered as maps without symmetry (i.e. EM maps).
           In some cases though P1 maps do/should have symmetry -
           and this is the means by you can tell Coot that.
    @param imol is the moleculle number to be acted on
    @param state the desired state, a value of 1 turns on map symmetry
*/
void set_map_has_symmetry(int imol, int state);

/*! \brief make a new map (a copy of map_no) that is in the cell,
  spacegroup and gridding of the map in reference_map_no.

Return the new map molecule number - return -1 on failure */
int reinterp_map(int map_no, int reference_map_no);

/*! \brief make a new map (a copy of map_no) that is in the cell,
  spacegroup and a multiple of the sampling of the input map (a
  sampling factor of more than 1 makes the output maps smoother) */
int smooth_map(int map_no, float sampling_multiplier);

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

/*! \brief Somewhat similar to the above function, except in this
case we overwrite the imol_map and we also presume that the
grid sampling of the contributing maps match. This makes it
much faster to generate than an average map.
*/
void regen_map_py(int imol_map, PyObject *map_number_and_scales);
#endif /* USE_PYTHON */
#endif /* c++ */

/* \} */

/*  ----------------------------------------------------------------------- */
/*                         (density) iso level increment entry */
/*  ----------------------------------------------------------------------- */
/* section Density Increment */
/*! \name  Density Increment */
/*! \{ */

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

/*! \brief return difference maps iso-map level increment  */
float get_diff_map_iso_level_increment();
/*! \brief set the difference maps iso-map level increment  */
void set_diff_map_iso_level_increment_from_text(const char *text, int imol);

/*! \brief sampling rate

find the molecule for which the single map dialog applies and set
    the contour level and redraw */
void set_map_sampling_rate_text(const char *text);

/*! \brief set the map sampling rate (default 1.5)

Set to something like 2.0 or 2.5 for more finely sampled maps.  Useful
for baton-building low resolution maps. */
void set_map_sampling_rate(float r);

/* MOVE-ME to c-interface-gtk-widgets.h */
char* get_text_for_map_sampling_rate_text();

/*! \brief return the map sampling rate */
float get_map_sampling_rate();

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

/*! \brief the function is a synonym of the above function - which now has an archaic-style
            name
*/
int read_mtz(const char *mtz_file_name,
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

/*! \brief make a map molecule from the give file name.

 If the file updates, then the map will be updated. */
int make_updating_map(const char *mtz_file_name,
		      const char *f_col, const char *phi_col,
		      const char *weight,
		      int use_weights, int is_diff_map);


void stop_updating_molecule(int imol);

#ifdef __cplusplus
#ifdef USE_GUILE
SCM refmac_parameters_scm(int imol);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *refmac_parameters_py(int imol);
#endif	/* USE_PYTHON */

#endif	/* __cplusplus */


/*! \brief does the mtz file have phases? */
/* We need to know if an mtz file has phases.  If it doesn't then we */
/*  go down a (new 20060920) different path. */
int mtz_file_has_phases_p(const char *mtz_file_name);

/*! \brief is the given filename an mtz file? */
int is_mtz_file_p(const char *filename);

/*! \brief does the given file have cns phases? */
int cns_file_has_phases_p(const char *cns_file_name);

void wrapped_auto_read_make_and_draw_maps(const char *filename);

void set_auto_read_do_difference_map_too(int i);
/*! \brief return the flag to do a difference map (too) on auto-read MTZ

   @return 0 means no, 1 means yes. */

int auto_read_do_difference_map_too_state();
/*! \brief set the expected MTZ columns for Auto-reading MTZ file.

  Not every program uses the default refmac labels ("FWT"/"PHWT") for
  its MTZ file.  Here we can tell coot to expect other labels so that
  coot can "Auto-open" such MTZ files.

  e.g. (set-auto-read-column-labels "2FOFCWT" "PH2FOFCWT" 0) */
 void set_auto_read_column_labels(const char *fwt, const char *phwt,
				 int is_for_diff_map_flag);


/* MOVE-ME to c-interface-gtk-widgets.h */
char* get_text_for_density_size_widget(); /* const gchar *text */
void set_density_size_from_widget(const char *text);

/* MOVE-ME to c-interface-gtk-widgets.h */
char *get_text_for_density_size_em_widget();
void set_density_size_em_from_widget(const char *text);

/*! \brief set the extent of the box/radius of electron density contours for x-ray maps */
void set_map_radius(float f);

/*! \brief set the extent of the box/radius of electron density contours for EM map*/
void set_map_radius_em(float radius);

/*! \brief another (old) way of setting the radius of the map */
void set_density_size(float f);

void set_map_radius_slider_max(float f);

/*! \brief Give me this nice message str when I start coot */
void set_display_intro_string(const char *str);

/*! \brief return the extent of the box/radius of electron density contours */
float get_map_radius();

/*! \brief not everyone likes coot's esoteric depth cueing system

  Pass an argument istate=1 to turn it off

 (this function is currently disabled). */
void set_esoteric_depth_cue(int istate);

/*! \brief native depth cueing system

  return the state of the esoteric depth cueing flag */
int  esoteric_depth_cue_state();

/*! \brief not everone likes coot's default difference map colouring.

   Pass an argument i=1 to swap the difference map colouring so that
   red is positive and green is negative. */
void set_swap_difference_map_colours(int i);
int swap_difference_map_colours_state();

/*! \brief post-hoc set the map of molecule number imol to be a
  difference map
  @return success status, 0 -> failure (imol does not have a map) */
int set_map_is_difference_map(int imol, short int bool_flag);

/*! \brief map is difference map? */
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
float density_at_point(int imol_map, float x, float y, float z);

/*! \} */


/*  ------------------------------------------------------------------------ */
/*                         Parameters from map:                              */
/*  ------------------------------------------------------------------------ */
/* section Parameters from map */
/*! \name  Parameters from map */
/*! \{ */

/*! \brief return the mtz file that was use to generate the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
const char *mtz_hklin_for_map(int imol_map);

/*! \brief return the FP column in the file that was use to generate
  the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say).

 Caller should dispose of returned pointer.
*/
const char *mtz_fp_for_map(int imol_map);

/*! \brief return the phases column in mtz file that was use to generate
  the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say).
 Caller should dispose of returned pointer.
*/
const char *mtz_phi_for_map(int imol_map);

/*! \brief return the weight column in the mtz file that was use to
  generate the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say) or no weights were used.
 Caller should dispose of returned pointer.
*/
const char *mtz_weight_for_map(int imol_map);

/*! \brief return flag for whether weights were used that was use to
  generate the map

  return 0 when no weights were used or there is no mtz file
  associated with that map. */
short int mtz_use_weight_for_map(int imol_map);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the parameter that made the map,

@return false or a string
  like ("xxx.mtz" "FPH" "PHWT" "" False) */
SCM map_parameters_scm(int imol);
/*! \brief return the parameter that made the map,

@return false or a list
  like (45 46 47 90 90 120), angles in degress */
SCM cell_scm(int imol);
/*! \brief return the parameter of the molecule, something
  like (45 46 47 90 90 120), angles in degress */
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the parameter that made the map,

@return False or something
  like ["xxx.mtz", "FPH", "PHWT", "", False] */
PyObject *map_parameters_py(int imol);
/*! \brief return the parameter that made the map,

@return False or something
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
/*! \{ */

/*! \brief write molecule number imol as a PDB to file file_name */
/*  return 0 on success, 1 on error. */
int write_pdb_file(int imol, const char *file_name);

/*! \brief write molecule number imol as a mmCIF to file file_name */
/*  return 0 on success, 1 on error. */
int write_cif_file(int imol, const char *file_name);

/*! \brief write molecule number imol's residue range as a PDB to file
  file_name */
/*  return 0 on success, 1 on error. */
int write_residue_range_to_pdb_file(int imol, const char *chainid,
				    int resno_start, int resno_end,
				    const char *filename);

/*  return 0 on success, -1 on error. */
int write_chain_to_pdb_file(int imol, const char *chainid, const char *filename);


/*! \brief save all modified coordinates molecules to the default
  names and save the state too. */
int quick_save();

/*! \brief return the state of the write_conect_records_flag.
  */
int get_write_conect_record_state();

/*! \brief set the flag to write (or not) conect records to the PDB file.
*/
void set_write_conect_record_state(int state);

/*! \} */



/*  ------------------------------------------------------------------------ */
/*                         Info Dialog                                       */
/*  ------------------------------------------------------------------------ */
/* section Info Dialog */
/*! \name  Info Dialog */
/*! \{ */

/*! \brief create a dialog with information

  create a dialog with information string txt.  User has to click to
  dismiss it, but it is not modal (nothing in coot is modal). */
void info_dialog(const char *txt);

/*! \brief create a dialog with information and print to console

  as info_dialog but print to console as well.  */
void info_dialog_and_text(const char *txt);

/*! \brief as above, create a dialog with information

This dialog is left-justified and can use markup such as angled bracketted tt or i
*/
void info_dialog_with_markup(const char *txt);


/*! \} */


/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */
/* section Refmac Functions */
/*! \name  Refmac Functions */
/*! \{ */
/*! \brief set counter for runs of refmac so that this can be used to
  construct a unique filename for new output */
void set_refmac_counter(int imol, int refmac_count);

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

/*! \} */

/*  --------------------------------------------------------------------- */
/*                      symmetry                                          */
/*  --------------------------------------------------------------------- */
/* section Symmetry Functions */
/*! \name Symmetry Functions */
/*! \{ */
char* get_text_for_symmetry_size_widget(); /* const gchar *text */

/* MOVE-ME to c-interface-gtk-widgets.h */
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

/*! \brief set symmetry colour by symop mode */
void set_symmetry_colour_by_symop(int imol, int state);
/*! \brief set symmetry colour for the chain */
void set_symmetry_whole_chain(int imol, int state);
/*! \brief set use expanded symmetry atom labels */
void set_symmetry_atom_labels_expanded(int state);

/*! \brief molecule number imol has a unit cell?

   @return 1 on "yes, it has a cell", 0 for "no" */
int has_unit_cell_state(int imol);

/* a gui function really */
void add_symmetry_on_to_preferences_and_apply();

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

/*! \brief return the number of symmetry operators for the given molecule

return -1 on no-symmetry for molecule or inappropriate imol number */
int n_symops(int imol);

/* This function works by active symm atom. */
int move_reference_chain_to_symm_chain_position();

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

/*! \brief set the space group for a coordinates molecule

 for shelx FA pdb files, there is no space group.  So allow the user
   to set it.  This can be initted with a HM symbol or a symm list for
   clipper.

This will only work on model molecules.

@return the success status of the setting  (1 good, 0 fail). */
short int set_space_group(int imol, const char *spg);

//! \brief set the unit cell for a given model molecule
//!
//! Angles in degress, cell lengths in Angstroms.
//!
//! @return  the success status of the setting (1 good, 0 fail).
int set_unit_cell_and_space_group(int imol, float a, float b, float c, float alpha, float beta, float gamma, const char *space_group);

//! \brief set the unit cell for a given model molecule using the cell of moecule imol_from
//!
//! This will only work on model molecules.
//! @return  the success status of the setting (1 good, 0 fail).
int set_unit_cell_and_space_group_using_molecule(int imol, int imol_from);

/*! \brief set the cell shift search size for symmetry searching.

When the coordinates for one (or some) symmetry operator are missing
(which happens sometimes, but rarely), try changing setting this to 2
(default is 1).  It slows symmetry searching, which is why it is not
set to 2 by default.  */
void set_symmetry_shift_search_size(int shift);

/*! \} */ /* end of symmetry functions */

/*  ------------------------------------------------------------------- */
/*                    file selection                                    */
/*  ------------------------------------------------------------------- */
/* section File Selection Functions */
/*! \name File Selection Functions */
/*! \{ */ /* start of file selection functions */

/* so that we can save/set the directory for future fileselections
   (i.e. the new fileselection will open in the directory that the
   last one ended in) */


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

/*! \} */ /* end of file selection functions */

/*  -------------------------------------------------------------------- */
/*                     history                                           */
/*  -------------------------------------------------------------------- */
/* section History Functions */
/*! \name  History Functions */

/*! \{ */ /* end of file selection functions */
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

/*! \} */

/*  --------------------------------------------------------------------- */
/*                  state (a graphics_info thing)                         */
/*  --------------------------------------------------------------------- */
/* info */
/*! \name State Functions */
/*! \{ */

/*! \brief save the current state to the default filename */
void save_state();

/*! \brief save the current state to file filename */
void save_state_file(const char *filename);

/*! \brief save the current state to file filename */
void save_state_file_py(const char *filename);

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
#endif /* USE_PYTHON */
/*! \brief run the state file depending on the state variables */
void run_state_file_maybe();	/* depending on the above state variables */


/*! \} */

/*  -------------------------------------------------------------------- */
/*                     virtual trackball                                 */
/*  -------------------------------------------------------------------- */
/* subsection Virtual Trackball */
/*! \name The Virtual Trackball */
/*! \{ */

#define VT_FLAT 1
#define VT_SPHERICAL 2

//! @param mode 1 for "Flat", 2 for "Spherical Surface" 
//!
void vt_surface(int mode);

//! @return the status, mode=1 for "Flat", mode=2 for "Spherical Surface"
int  vt_surface_status();

/*! \} */


/*  --------------------------------------------------------------------- */
/*                      clipping                                          */
/*  --------------------------------------------------------------------- */
/* section Clipping Functions */
/*! \name  Clipping Functions */
/*! \{ */

//! increase the amount of clipping, that is (independent of projection matrix)
void increase_clipping_front();

//! increase the amount of clipping, that is (independent of projection matrix)
void increase_clipping_back();

//! decrease the amount of clipping, that is (independent of projection matrix)
void decrease_clipping_front();

//! decrease the amount of clipping, that is (independent of projection matrix)
void decrease_clipping_back();

//! set clipping plane back  - this goes in differnent directions for orthographics vs perspective 
void set_clipping_back(float v);

//! set clipping plane front - this goes in differnent directions for orthographics vs perspective  
void set_clipping_front(float v);

//! get clipping plane front 
float get_clipping_plane_front();

//! get clipping plane back 
float get_clipping_plane_back();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Unit Cell                                        */
/*  ----------------------------------------------------------------------- */
/* section Unit Cell interface */
/*! \name  Unit Cell interface */
/*! \{ */

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
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Colour                                           */
/*  ----------------------------------------------------------------------- */
/* section Colour */
/*! \name  Colour */
/*! \{ */

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

/*! \brief colour molecule number imol by chain type */
void set_colour_by_ncs_chain(int imol, short int goodsell_mode);

/*! \brief colour molecule number imol by chain type, goodsell-like colour scheme */
void set_colour_by_chain_goodsell_mode(int imol);

/*! \brief set the goodsell chain colour colour wheel step (default 0.22) */
void set_goodsell_chain_colour_wheel_step(float s);

/*! \brief colour molecule number imol by molecule */
void set_colour_by_molecule(int imol);

/* get the value of graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag */
int get_colour_map_rotation_on_read_pdb_c_only_flag();

/*! \brief set the symmetry colour base */
void set_symmetry_colour(float r, float g, float b);

/*! \} */

/*  Section Map colour*/
/*! \name   Map colour*/
/*! \{ */
/*! \brief set the colour map rotation (hue change) for maps

   default: for maps is 14 degrees. */
void set_colour_map_rotation_for_map(float f); /* "global"/default */

/*! \brief set the colour map rotation for molecule number imol

theta is in degrees */
void set_molecule_bonds_colour_map_rotation(int imol, float theta);

/*! \brief Get the colour map rotation for molecule number imol */
float get_molecule_bonds_colour_map_rotation(int imol);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Anisotropic Atoms */
/*  ----------------------------------------------------------------------- */
/* section Anisotropic Atoms Interface */
/*! \name  Anisotropic Atoms Interface */
/*! \{ */
/*  we use the text interface to this in callback.c rather */
/*  than getting the float directly. */

/*! \brief get the aniso radius limit */
float get_limit_aniso();           /* not a function of the molecule */

/*! \brief get show the aniso limit */
short int get_show_limit_aniso();  /* not a function of the molecule */

/*! \brief return show-aniso-atoms state  - FIXME- per molecule */
short int get_show_aniso();       /*  not a function of the molecule */

/*! \brief set the aniso atom limit */
void set_limit_aniso(short int state);

/*! \brief does nothing */
void set_show_aniso(int state);

/*! \brief set show aniso atoms */
void set_show_aniso_atoms(int imol, int state);

/*! \brief set show aniso atoms as ortep */
void set_show_aniso_atoms_as_ortep(int imol, int state);

/* DELETE-ME */
void set_aniso_limit_size_from_widget(const char *text);

/* DELETE-ME .h */
char *get_text_for_aniso_limit_radius_entry();

/*! DELETE-ME */
void set_aniso_probability(float f);

/*! DELETE-ME  */
float get_aniso_probability();

/*! \} */

/*  ---------------------------------------------------------------------- */
/*                         Display Functions                               */
/*  ---------------------------------------------------------------------- */
/* section Display Functions */
/*! \name  Display Functions */
/*! \{ */
/*  currently doesn't get seen when the window starts due to */
/*  out-of-order issues. */

/*! \brief set the window size */
void   set_graphics_window_size(int x_size, int y_size);
/*! \brief set the window size as gtk_widget (flag=1) or gtk_window (flag=0) */
void   set_graphics_window_size_internal(int x_size, int y_size, int as_widget_flag);
/*! \brief set the graphics window position */
void   set_graphics_window_position(int x_pos, int y_pos);
/*! \brief store the graphics window position */
void store_graphics_window_position(int x_pos, int y_pos); /*  "configure_event" callback */

/*! \brief store the graphics window position and size to zenops-graphics-window-size-and-postion.scm in
 *         the preferences directory. */
void graphics_window_size_and_position_to_preferences();

/*! \brief draw a frame */
void graphics_draw(); 	/* and wrapper interface to gtk_widget_draw(glarea)  */

/*! \brief try to turn on Zalman stereo mode  */
void zalman_stereo_mode();
/*! \brief try to turn on stereo mode  */
void hardware_stereo_mode();


/*! \brief set the stereo mode (the relative view of the eyes)

0 is 2010-mode
1 is modern mode
*/
void set_stereo_style(int mode);

/*! \brief what is the stero state?

  @return 1 for in hardware stereo, 2 for side by side stereo, else return 0. */
int  stereo_mode_state();
/*! \brief try to turn on mono mode  */
void mono_mode();

/*! \brief turn on side bye side stereo mode
 *
 * @param use_wall_eye_mode 1 mean wall-eyed, 0 means cross-eyed
 * */
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

void set_model_display_radius(int state, float radius);

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
/*! \brief set rotamer selection dialog position */
void set_rotamer_selection_dialog_position(int x_pos, int y_pos);

/*! \} */

/*  ---------------------------------------------------------------------- */
/*                         Smooth "Scrolling" */
/*  ---------------------------------------------------------------------- */
/* section Smooth Scrolling */
/*! \name  Smooth Scrolling */
/*! \{ */

/*! \brief set smooth scrolling

  @param v use v=1 to turn on smooth scrolling, v=0 for off (default on). */
void set_smooth_scroll_flag(int v);

/*! \brief return the smooth scrolling state */
int  get_smooth_scroll();

/* MOVE-ME to c-interface-gtk-widgets.h */
void set_smooth_scroll_steps_str(const char * t);

/*  useful exported interface */
/*! \brief set the number of steps in the smooth scroll

   Set more steps (e.g. 50) for more smoothness (default 10).*/
void set_smooth_scroll_steps(int i);

/* MOVE-ME to c-interface-gtk-widgets.h */
char  *get_text_for_smooth_scroll_steps();

/* MOVE-ME to c-interface-gtk-widgets.h */
void  set_smooth_scroll_limit_str(const char *t);

/*  useful exported interface */
/*! \brief do not scroll for distances greater this limit */
void  set_smooth_scroll_limit(float lim);

char *get_text_for_smooth_scroll_limit();

/*! \} */


/*  ---------------------------------------------------------------------- */
/*                         Font Size */
/*  ---------------------------------------------------------------------- */
/* section Font Parameters */
/*! \name  Font Parameters */
/*! \{ */

/*! \brief set the font size

  @param i 1 (small) 2 (medium, default) 3 (large) */
void set_font_size(int i);

/*! \brief return the font size

  @return 1 (small) 2 (medium, default) 3 (large) */
int get_font_size();

/*! \brief set the colour of the atom label font - the arguments are
  in the range 0->1 */
void set_font_colour(float red, float green, float blue);

/*! \brief set use stroke characters */
void set_use_stroke_characters(int state);

/*! \} */

/*  ---------------------------------------------------------------------- */
/*                         Rotation Centre                                 */
/*  ---------------------------------------------------------------------- */
/* section Rotation Centre */
/*! \name  Rotation Centre */
/*! \{ */

/* 20220723-PE I agree with my comments from earlier - these should not be here */
/* MOVE-ME to c-interface-gtk-widgets.h */
void set_rotation_centre_size_from_widget(const gchar *text); /* and redraw */
/* MOVE-ME to c-interface-gtk-widgets.h */
gchar *get_text_for_rotation_centre_cube_size();

/*! \brief set the rotation centre marker size */
void set_rotation_centre_size(float f); /* and redraw (maybe) */

/*! \brief set the rotation centre marker size */
void set_user_defined_rotation_centre_crosshairs_size_scale_factor(float f);

/*! \brief set rotation centre colour

This is the colour for a dark background - if the background colour is not dark,
then the cross-hair colour becomes the inverse colour */
void set_rotation_centre_cross_hairs_colour(float r, float g, float b, float alpha);

/*! \brief return the recentre-on-pdb state */
short int recentre_on_read_pdb();
/*! \brief set the recentre-on-pdb state */
void set_recentre_on_read_pdb(short int);

/*! \brief set the rotation centre */
void set_rotation_centre(float x, float y, float z);
/* The redraw happens somewhere else... */
void set_rotation_centre_internal(float x, float y, float z);
float rotation_centre_position(int axis); /* only return one value: x=0, y=1, z=2 */
/*! \brief centre on the ligand of the "active molecule", if we are
  already there, centre on the next hetgroup (etc) */
void go_to_ligand();

#ifdef USE_PYTHON
#ifdef __cplusplus
PyObject *go_to_ligand_py();
#endif
#endif

/*! \brief go to the ligand that has more than n_atom_min atoms */
void set_go_to_ligand_n_atoms_limit(int n_atom_min);

/*! \brief rotate the view so that the next main-chain atoms are oriented
 in the same direction as the previous - hence side-chain always seems to be
"up" - set this mode to 1 for reorientation-mode - and 0 for off (standard translation)
*/
void set_reorienting_next_residue_mode(int state);

/*! \} */

/*  ---------------------------------------------------------------------- */
/*                         orthogonal axes                                 */
/*  ---------------------------------------------------------------------- */

/* section Orthogonal Axes */
/*! \name Orthogonal Axes */
/*! \{ */
/* Draw the axes in the top left?

0 off, 1 on */
void set_draw_axes(int i);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  utility function                                        */
/*  ----------------------------------------------------------------------- */
/* section Atom Selection Utilities */
/*! \name  Atom Selection Utilities */
/*! \{ */

#ifdef __cplusplus /* protection from use in callbacks.c, else compilation probs */
#ifdef USE_PYTHON
/* Get model molecule list */
PyObject *get_model_molecule_list_py();
#endif
#endif

#ifdef __cplusplus
#ifdef USE_GUILE
/* Get model molecule list */
SCM get_model_molecule_list_scm();
#endif
#endif

/* does not account for alternative conformations properly */
/* return -1 if atom not found. */
int atom_index(int imol, const char *chain_id, int iresno, const char *atom_id);
/* using alternative conformations properly ?! */
/* return -1 if atom not found. */
int atom_index_full(int imol, const char *chain_id, int iresno, const char *inscode, const char *atom_id, const char *altconf);
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
/*! \brief return the minimum residue number for imol chain chain_id */
int min_resno_in_chain(int imol, const char *chain_id);
/*! \brief return the maximum residue number for imol chain chain_id */
int max_resno_in_chain(int imol, const char *chain_id);
/*! \brief return the median temperature factor for imol */
float median_temperature_factor(int imol);
/*! \brief return the average temperature factor for the atoms in imol */
float average_temperature_factor(int imol);
/*! \brief return the standard deviation of the atom temperature factors for imol */
float standard_deviation_temperature_factor(int imol);

/*! \brief clear pending picks (stop coot thinking that the user is about to pick an atom).  */
void clear_pending_picks();
char *centre_of_mass_string(int imol);
#ifdef USE_PYTHON
char *centre_of_mass_string_py(int imol);
#endif
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

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
void set_temperature_factors_for_atoms_in_residue_scm(int imol, SCM residue_spec_scm, float bf);
#endif
#endif

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM get_residue_alt_confs_scm(int imol, const char *chain_id, int res_no, const char *ins_code);
#endif
#endif

#ifdef __cplusplus /* protection from use in callbacks.c, else compilation probs */
#ifdef USE_PYTHON
/*! \brief Return either False (on failure) or a list of alt-conf strings (might be [""]) */
PyObject *get_residue_alt_confs_py(int imol, const char *chain_id, int res_no, const char *ins_code);
#endif
#endif



/*! \brief swap atom alt-confs */
int swap_atom_alt_conf(int imol, const char *chain_id, int res_no, const char *ins_code,
                       const char *atom_name, const char*alt_conf);

/*! \brief swap atom alt-confs */
int swap_residue_alt_confs(int imol, const char *chain_id, int res_no, const char *ins_code);

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
#endif /* __cplusplus */

/*! \brief set the residue name of the specified residue */
void set_residue_name(int imol, const char *chain_id, int res_no, const char *ins_code, const char *new_residue_name);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                            skeletonization                               */
/*  ----------------------------------------------------------------------- */
/* section Skeletonization Interface */
/*! \name  Skeletonization Interface */
/*! \{ */
void skel_greer_on();
void skel_greer_off();

/*! \brief skeletonize molecule number imol

   the prune_flag should almost  always be 0.

   NOTE:: The arguments to have been reversed for coot 0.8.3 and later
   (now the molecule number comes first).

    */
int skeletonize_map(int imol, short int prune_flag);

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



/*  ----------------------------------------------------------------------- */
/*                  skeletonization level widgets                           */
/*  ----------------------------------------------------------------------- */

/* MOVE-ME to c-interface-gtk-widgets.h */
gchar *get_text_for_skeletonization_level_entry();

/* MOVE-ME to c-interface-gtk-widgets.h */
void set_skeletonization_level_from_widget(const char *txt);

/* MOVE-ME to c-interface-gtk-widgets.h */
gchar *get_text_for_skeleton_box_size_entry();

/* MOVE-ME to c-interface-gtk-widgets.h */
void set_skeleton_box_size_from_widget(const char *txt);


/*! \brief the box size (in Angstroms) for which the skeleton is displayed */
void set_skeleton_box_size(float f);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                        save coordinates                                  */
/*  ----------------------------------------------------------------------- */
/* section Save Coordinates */
/*! \name  Save Coordinates */
/*! \{ */


/*! \brief save coordinates of molecule number imol in filename

  @return status 1 is good (success), 0 is fail. */
int save_coordinates(int imol, const char *filename);

/*! \brief set save coordinates in the starting directory */
void set_save_coordinates_in_original_directory(int i);

/* access to graphics_info_t::save_imol for use in callback.c */
int save_molecule_number_from_option_menu();
/* access from callback.c, not to be used in scripting, I suggest.
   Sets the *save* molecule number */
void set_save_molecule_number(int imol);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                        .phs file reading                                 */
/*  ----------------------------------------------------------------------- */
/* section Read Phases File Functions */
/*! \name  Read Phases File Functions */
/*! \{ */

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

// 20220723-PE MOVE-ME!
void
graphics_store_phs_filename(const gchar *phs_filename);

short int possible_cell_symm_for_phs_file();

/* MOVE-ME to c-interface-gtk-widgets.h */
gchar *get_text_for_phs_cell_chooser(int imol, const char *field);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                                  Movement                                */
/*  ----------------------------------------------------------------------- */
/* section Graphics Move */
/*! \name Graphics Move */
/*! \{ */
/*! \brief undo last move  */
void undo_last_move(); /* suggested by Frank von Delft */

/*! \brief translate molecule number imol by (x,y,z) in Angstroms  */
void translate_molecule_by(int imol, float x, float y, float z);

/*! \brief transform molecule number imol by the given rotation
  matrix, then translate by (x,y,z) in Angstroms  */
void transform_molecule_by(int imol,
			   float m11, float m12, float m13,
			   float m21, float m22, float m23,
			   float m31, float m32, float m33,
			   float x, float y, float z);

/*! \brief transform fragment of molecule number imol by the given rotation
  matrix, then translate by (x,y,z) in Angstroms  */
void transform_zone(int imol, const char *chain_id, int resno_start, int resno_end, const char *ins_code,
		    float m11, float m12, float m13,
		    float m21, float m22, float m23,
		    float m31, float m32, float m33,
		    float x, float y, float z);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                        go to atom widget                                 */
/*  ----------------------------------------------------------------------- */
/* section Go To Atom Widget Functions */
/*! \name Go To Atom Widget Functions */
/*! \{ */

/*! \brief Post the Go To Atom Window */
void post_go_to_atom_window();

/*! \brief the go-to-atom molecule number */
int go_to_atom_molecule_number();
/*! \brief the go-to-atom chain-id */
char *go_to_atom_chain_id();
/*! \brief the go-to-atom atom name */
char *go_to_atom_atom_name();
/*! \brief the go-to-atom residue number */
int go_to_atom_residue_number();
/*! \brief the go-to-atom insertion code */
char *go_to_atom_ins_code();
/*! \brief the go-to-atom alt conf */
char *go_to_atom_alt_conf();


/*! \brief set the go to atom specification

   It seems important for swig that the `char *` arguments are `const
   char *`, not `const gchar *` (or else we get wrong type of argument
   error on (say) "A"

@return the success status of the go to.  0 for fail, 1 for success.
*/
int set_go_to_atom_chain_residue_atom_name(const char *t1_chain_id, int iresno,
					   const char *t3_atom_name);

/*! \brief set the go to (full) atom specification

   It seems important for swig that the `char *` arguments are `const
   char *`, not `const gchar *` (or else we get wrong type of argument
   error on (say) "A"

@return the success status of the go to.  0 for fail, 1 for success.
*/
int set_go_to_atom_chain_residue_atom_name_full(const char *chain_id,
						int resno,
						const char *ins_code,
						const char *atom_name,
						const char *alt_conf);
/*! \brief set go to atom but don't redraw */
int set_go_to_atom_chain_residue_atom_name_no_redraw(const char *t1, int iresno, const char *t3,
						     short int make_the_move_flag);

// MOVE-ME!
int set_go_to_atom_chain_residue_atom_name_strings(const gchar *t1,
						   const gchar *t2,
						   const gchar *txt);


/*! \brief update the Go To Atom widget entries to atom closest to
  screen centre. */
void update_go_to_atom_from_current_position();


/* moving gtk function out of build functions, delete_atom() updates
   the go to atom atom list on deleting an atom  */
void update_go_to_atom_residue_list(int imol);

/*  return an atom index */
/*! \brief what is the atom index of the given atom? */
int atom_spec_to_atom_index(int mol, const char *chain, int resno, const char *atom_name);

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

/* MOVE-ME to c-interface-gtk-widgets.h */
void unset_go_to_atom_widget(); /* unstore the static go_to_atom_window */



/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  autobuilding control                                    */
/*  ----------------------------------------------------------------------- */
/* section AutoBuilding functions (Defunct) */
/* void autobuild_ca_on();  - moved to junk */

void autobuild_ca_off();

void test_fragment();

void do_skeleton_prune();

int test_function(int i, int j);

#ifdef __cplusplus
#ifdef USE_GUILE
SCM test_function_scm(SCM i, SCM j);
#endif
#ifdef USE_PYTHON
PyObject *test_function_py(PyObject *i, PyObject *j);
#endif /* PYTHON */
#endif


/*                    glyco tools test  */
void glyco_tree_test();

#ifdef __cplusplus
#ifdef USE_GUILE
SCM glyco_tree_scm(int imol, SCM active_residue_scm);
SCM glyco_tree_residues_scm(int imol, SCM active_residue_scm);
SCM glyco_tree_internal_distances_fn_scm(int imol, SCM residue_spec, const std::string &file_name); // testing function
SCM glyco_tree_residue_id_scm(int imol, SCM residue_spec_scm);
SCM glyco_tree_compare_trees_scm(int imol_1, SCM res_spec_1, int imol_2, SCM res_spec_2);
SCM glyco_tree_matched_residue_pairs_scm(int imol_1, SCM res_spec_1, int imol_2, SCM res_spec_2);
#endif
#ifdef USE_PYTHON
PyObject *glyco_tree_py(int imol, PyObject *active_residue_py);
PyObject *glyco_tree_residues_py(int imol, PyObject *active_residue_py);
PyObject *glyco_tree_internal_distances_fn_py(int imol, PyObject *residue_spec, const std::string &file_name); // testing function
PyObject *glyco_tree_residue_id_py(int imol, PyObject *residue_spec_py);
PyObject *glyco_tree_compare_trees_py(int imol_1, PyObject *res_spec_1, int imol_2, PyObject *res_spec_2);
PyObject *glyco_tree_matched_residue_pairs_py(int imol_1, PyObject *res_spec_1, int imol_2, PyObject *res_spec_2);
#endif /* PYTHON */
#endif



/*  ----------------------------------------------------------------------- */
/*                  map and molecule control                                */
/*  ----------------------------------------------------------------------- */
/* section Map and Molecule Control */
/*! \name Map and Molecule Control */
/*! \{ */

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

/*! \brief from all the model molecules, display only imol

This stops flashing/delayed animations with many molecules */
void set_display_only_model_mol(int imol);

/*! \brief make the coordinates molecule active/inactve (clickable), 0
  for off, 1 for on */
void set_mol_active(int imol, int state);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief maps_list is a list of integers (map molecule numbers).

This interface is uses so that we don't get flashing when a map is turned off (using set_mol_displayed). */
void display_maps_scm(SCM maps_list);
#endif
#ifdef USE_PYTHON
void display_maps_py(PyObject *pyo);
#endif
#endif



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

/*! \brief only display the last model molecule
*/
void set_only_last_model_molecule_displayed();

/*! \brief display only the active mol and the refinement map */
void display_only_active();


#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the spacegroup as a string, return scheme false if unable to do so. */
SCM space_group_scm(int imol);
#endif
#ifdef USE_PYTHON
PyObject *space_group_py(int imol);
#endif
#endif

/*! \brief return the spacegroup of molecule number imol . Deprecated.

@return "No Spacegroup" when the spacegroup of a molecule has not been
set. */
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


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Merge Molecules                                  */
/*  ----------------------------------------------------------------------- */
/* section Merge Molecules */

/*! \brief merge molecules

@return a pair, the first item of which is a status (1 is good) the second is
a list of merge-infos (one for each of the items in add_molecules). If the
molecule of an add_molecule item is just one residue, return a spec for the
new residue, if it is many residues return a chain id.

the first argument is a list of molecule numbers and the second is the target
   molecule into which the others should be merged  */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM merge_molecules(SCM add_molecules, int imol);
void set_merge_molecules_ligand_spec_scm(SCM ligand_spec_scm);
#endif

#ifdef USE_PYTHON


/*! \brief merge molecules

@return a pair, the first item of which is a status (1 is good) the second is
a list of merge-infos (one for each of the items in add_molecules). If the
molecule of an add_molecule item is just one residue, return a spec for the
new residue, if it is many residues return a chain id.

the first argument is a list of molecule numbers and the second is the target
   molecule into which the others should be merged  */
PyObject *merge_molecules_py(PyObject *add_molecules, int imol);
void set_merge_molecules_ligand_spec_py(PyObject *ligand_spec_py);
#endif /* PYTHON */
#endif	/* c++ */


/*  ----------------------------------------------------------------------- */
/*                         Align and Mutate GUI                             */
/*  ----------------------------------------------------------------------- */
/* section Align and Mutate */
/*! \name  Align and Mutate */
/*! \{ */

/*! \brief align and mutate the given chain to the given sequence  */
void align_and_mutate(int imol, const char *chain_id, const char *fasta_maybe, short int renumber_residues_flag);
/*! \brief set the penalty for affine gap and space when aligning, defaults -3.0 and -0.4 */
void set_alignment_gap_and_space_penalty(float wgap, float wspace);


/* What are these functions?  consider deleting them - we have alignment_mismatches_* . */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM alignment_results_scm(int imol, const char* chain_id, const char *seq);
/*! \brief return the residue spec of the nearest residue by sequence
  numbering.

@return  scheme false if not possible */
SCM nearest_residue_by_sequence_scm(int imol, const char* chain_id, int resno, const char *ins_code);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *alignment_results_py(int imol, const char* chain_id, const char *seq);
/*! \brief return the residue spec of the nearest residue by sequence
  numbering.  Return Python False if not possible */
PyObject *nearest_residue_by_sequence_py(int imol, const char* chain_id, int resno, const char *ins_code);
#endif /* USE_PYTHON */
#endif  /* c++ */

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Renumber residue range                           */
/*  ----------------------------------------------------------------------- */
/* section Renumber Residue Range */
/*! \name Renumber Residue Range */

/*! \{ */
/*! \brief renumber the given residue range by offset residues */
int renumber_residue_range(int imol, const char *chain_id,
			   int start_res, int last_res, int offset);


/*! \brief change chain id, residue number or insertion code for given
  residue  */
int change_residue_number(int imol, const char *chain_id, int current_resno, const char *current_inscode, int new_resno, const char *new_inscode);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Change chain id                                  */
/*  ----------------------------------------------------------------------- */
/* section Change Chain ID */
/*! \name Change Chain ID */
/*! \{ */

/*! \brief change the chain id of the specified residue */
void  change_chain_id(int imol, const char *from_chain_id, const char *to_chain_id,
		      short int use_res_range_flag, int from_resno, int to_resno);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM change_chain_id_with_result_scm(int imol, const char *from_chain_id, const char *to_chain_id,
                                         short int use_res_range_flag, int from_resno, int to_resno);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *change_chain_id_with_result_py(int imol, const char *from_chain_id, const char *to_chain_id, short int use_res_range_flag, int from_resno, int to_resno);
#endif /* USE_PYTHON */
#endif  /* c++ */
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  scripting                                               */
/*  ----------------------------------------------------------------------- */
/* section Scripting Interface */

/*! \name Scripting Interface */
/*! \{ */

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


/* called from c-inner-main */
void run_command_line_scripts();

void set_guile_gui_loaded_flag();
void set_python_gui_loaded_flag();
void set_found_coot_gui();
void set_found_coot_python_gui();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Monomer                                                 */
/*  ----------------------------------------------------------------------- */
/* section monomers */

/*! \name Monomer */
/*! \{ */

int get_monomer_for_molecule_by_index(int dict_idx, int imol_enc);


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
/*! \brief guile run script file */
void run_guile_script (const char *filename);
/*! \brief run python script file */
void run_python_script(const char *filename);
/*! \brief import python module */
int import_python_module(const char *module_name, int use_namespace);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return a list of compoundIDs in the dictionary of which the
  given string is a substring of the compound name */
SCM matching_compound_names_from_dictionary_scm(const char *compound_name_fragment,
						short int allow_minimal_descriptions_flag);

/*! \brief return the monomer name

  return scheme false if not found */
SCM comp_id_to_name_scm(const char *comp_id);
#endif /* USE_GUILE */

/*! \brief try to auto-load the dictionary for comp_id from the refmac monomer library.

   return 0 on failure.
   return 1 on successful auto-load.
   return 2 on already-read.
   */
int auto_load_dictionary(const char *comp_id);
/*! \brief as above, but dictionary is cleared and re-read if it already exists */
int reload_dictionary(const char *comp_id);

/*! \brief add residue name to the list of residue names that don't
  get auto-loaded from the Refmac dictionary. */
void add_non_auto_load_residue_name(const char *s);
/*! \brief remove residue name from the list of residue names that don't
  get auto-loaded from the Refmac dictionary. */
void remove_non_auto_load_residue_name(const char *s);

#ifdef USE_PYTHON
/*! \brief return a list of compoundIDs in the dictionary which the
  given string is a substring of the compound name */
PyObject *matching_compound_names_from_dictionary_py(const char *compound_name_fragment,
						     short int allow_minimal_descriptions_flag);
/*! \brief return the monomer name

  return python false if not found */
PyObject *comp_id_to_name_py(const char *comp_id);
#endif /* USE_PYTHON */
#endif /*__cplusplus */


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  regularize/refine                                       */
/*  ----------------------------------------------------------------------- */
/* section Regularization and Refinement */
/*! \name  Regularization and Refinement */
/*! \{ */

void do_regularize(short int state); /* pass 0 for off (unclick togglebutton) */
void do_refine(short int state);

/*! \brief add a restraint on peptides to make them planar

  This adds a 5 atom restraint that includes both CA atoms of the
  peptide.  Use this rather than editting the mon_lib_list.cif file. */
void add_planar_peptide_restraints();

/*! \brief remove restraints on peptides to make them planar. */
void remove_planar_peptide_restraints();

/*! \brief make the planar peptide restraints tight

Useful when refining models with cryo-EM maps */
void make_tight_planar_peptide_restraints();


/* return 1 if planar peptide restraints are on, 0 if off */
int planar_peptide_restraints_state();


/*! \brief add a restraint on peptides to keep trans peptides trans

i.e. omega in trans-peptides is restraints to 180 degrees.
 */
void set_use_trans_peptide_restraints(short int on_off_state);

/*! \brief add restraints on the omega angle of the peptides

  (that is the torsion round the peptide bond).  Omega angles that are
  closer to 0 than to 180 will be refined as cis peptides (and of
  course if omega is greater than 90 then the peptide will be refined
  as a trans peptide (this is the normal case). */
void add_omega_torsion_restriants();

/*! \brief remove omega restraints on CIS and TRANS linked residues. */
void remove_omega_torsion_restriants();

/*! \brief add or remove auto H-bond restraints */
void set_refine_hydrogen_bonds(int state);


/*! \brief set immediate replacement mode for refinement and regularization
 *
 * This can enable synchronous refinement (with istate = 1).
 * You need this (call with istate=1) if you are
 * scripting refinement/regularization
 *
 * @param istate set the state of immediate-refinemnt 
 * */
void set_refinement_immediate_replacement(int istate);

/*! \brief query the state of the immediate replacement mode */
int  refinement_immediate_replacement_state();

void set_refine_use_noughties_physics(short int state);

int get_refine_use_noughties_physics_state();

/*! \brief set the number of frames for which the selected residue
  range flashes

 On fast computers, this can be set to higher than the default for
 more aesthetic appeal. */
void set_residue_selection_flash_frames_number(int i);

/*! \brief accept the new positions of the regularized or refined residues

    If you are scripting refinement and/or regularization, this is not the
    function that you need to call after refine-zone or regularize-zone.
    If you are using Python, use accept_moving_atoms_py() and that will
    provide a return value that may be of some use.
*/
void c_accept_moving_atoms();

/*! \brief a hideously-named alias for `c_accept_moving_atoms()`  */
void accept_regularizement();

/*! \brief clear up moving atoms 
 */
void clear_up_moving_atoms();	/* remove the molecule and bonds */

/*! \brief remove just the bonds

   A redraw is done.
 */
void clear_moving_atoms_object(); /* just get rid of just the bonds (redraw done here). */

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */


/*! \brief If there is a refinement on-going already, we don't want to start a new one

The is the means to ask if that is the case. This needs a scheme wrapper to provide refinement-already-ongoing?
The question is translated to "are the intermediate atoms being displayed?" so that might be a more
accurate function name than the current one.

@return 1 for yes, 0 for no.
*/
short int refinement_already_ongoing_p();

#ifdef USE_GUILE
/*! \brief refine residues, r is a list of residue specs.

 @return refinement results, which consists of

   1 - an information string (in case of error)

   2 - the progress variable (from GSL)

   3 - refinement results for each particular geometry type (bonds, angles etc.)

 */
SCM refine_residues_scm(int imol, SCM r); /* presumes the alt_conf is "". */
SCM refine_residues_with_alt_conf_scm(int imol, SCM r, const char *alt_conf);
SCM refine_residues_with_modes_with_alt_conf_scm(int imol, SCM residues_spec_list_scm,
						 const char *alt_conf,
						 SCM mode_1,
						 SCM mode_2,
						 SCM mode_3);
SCM regularize_residues_scm(int imol, SCM r); /* presumes the alt_conf is "". */
SCM regularize_residues_with_alt_conf_scm(int imol, SCM r, const char *alt_conf);
#endif
#ifdef USE_PYTHON
/*! \brief refine the residues in the given residue spec list

@return the refinement summary statistics  */
PyObject *refine_residues_py(int imol, PyObject *r);  /* presumes the alt_conf is "". */
PyObject *refine_residues_with_modes_with_alt_conf_py(int imol, PyObject *r, const char *alt_conf,
						      PyObject *mode_1,
						      PyObject *mode_2,
						      PyObject *mode_3);
PyObject *refine_residues_with_alt_conf_py(int imol, PyObject *r, const char *alt_conf);
PyObject *regularize_residues_py(int imol, PyObject *r);  /* presumes the alt_conf is "". */
PyObject *regularize_residues_with_alt_conf_py(int imol, PyObject *r, const char *alt_conf);
#endif /* PYTHON */
#endif /* c++ */

/* Used by on_accept_reject_refinement_reject_button_clicked() */
void stop_refinement_internal();

void set_refinement_use_soft_mode_nbc_restraints(short int flag);

/*! \brief shiftfield B-factor refinement */
void shiftfield_b_factor_refinement(int imol);

/*! \brief shiftfield xyz refinement */
void shiftfield_xyz_factor_refinement(int imol);

/*! \brief turn on (or off) torsion restraints

   Pass with istate=1 for on, istate=0 for off.
*/
void set_refine_with_torsion_restraints(int istate);
/*! \brief return the state of above */
int refine_with_torsion_restraints_state();

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

/*! \brief return the relative weight of the geometric terms to the map terms.

A more sensible name for the matrix_state() function) */
float get_map_weight();

float estimate_map_weight(int imol_map);


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
/*! \brief repeat the previous (user-selected) refine zone */
void repeat_refine_zone(); /* use stored atom indices to re-run the refinement using the same atoms as previous */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM refine_zone_with_score_scm(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
SCM regularize_zone_with_score_scm(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
#endif /* guile */
#ifdef USE_PYTHON
PyObject *refine_zone_with_score_py(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
PyObject *regularize_zone_with_score_py(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
#endif /* PYTHON */
#endif /* c++ */

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

  Default: 80. */
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

/*! \brief - the elasticity of the dragged atom in refinement mode.

Default 0.33

 Bigger numbers mean bigger movement of the other atoms.*/
void set_refinement_drag_elasticity(float e);

/*! \brief turn on Ramachandran angles refinement in refinement and regularization */
/*! name consistent with set_refine_with_torsion_restraints() !?  */
void set_refine_ramachandran_angles(int state);
void set_refine_ramachandran_torsion_angles(int state);

/*! \brief change the target function type  */
void set_refine_ramachandran_restraints_type(int type);
/*! \brief change the target function weight

a big number means bad things  */
void set_refine_ramachandran_restraints_weight(float w);

/*! ramachandran restraints weight

@return weight as a float */
float refine_ramachandran_restraints_weight();

/* not ready yet \brief set the weight for torsion restraints (default 1.0)*/
void set_torsion_restraints_weight(double w);

/*! \brief set the state for using rotamer restraints "drive" mode

1 in on, 0 is off (off by default) */
void set_refine_rotamers(int state);

void set_refinement_geman_mcclure_alpha_from_text(int combobox_item_idx, const char *t);
void set_refinement_lennard_jones_epsilon_from_text(int combobox_item_idx, const char *t);
void set_refinement_ramachandran_restraints_weight_from_text(int combobox_item_idx, const char *t);
void set_refinement_overall_weight_from_text(const char *t);
void set_refinement_torsion_weight_from_text(int combobox_item_index, const char *t);
void set_refine_params_dialog_more_control_frame_is_active(int state);


int refine_ramachandran_angles_state();

void set_numerical_gradients(int istate);

void set_debug_refinement(int state);


/*! \brief correct the sign of chiral volumes before commencing refinement?

   Do we want to fix chiral volumes (by moving the chiral atom to the
   other side of the chiral plane if necessary).  Default yes
   (1). Note: doesn't work currently. */
void set_fix_chiral_volumes_before_refinement(int istate);

/*! \brief query the state of the above option */
void check_chiral_volumes(int imol);

#ifdef __cplusplus
#ifdef USE_GUILE
SCM chiral_volume_errors_scm(int imol);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *chiral_volume_errors_py(int imol);
#endif	/* USE_PYTHON */
#endif	/* __cplusplus */


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
int does_residue_exist_p(int imol, const char *chain_id, int resno, const char *inscode);

/*! \brief delete the restraints for the given comp_id (i.e. residue name)

@return success status (0 is failed, 1 is success)
*/
int delete_restraints(const char *comp_id);

/*! \brief add a user-define bond restraint

   this extra restraint is used when the given atoms are selected in
   refinement or regularization.

   @return the index of the new restraint.

   @return -1 when the atoms were not found and no extra bond
   restraint was stored.  */

int add_extra_bond_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2, double bond_dist, double esd);

/*! \brief add a user-define GM distance restraint

   this extra restraint is used when the given atoms are selected in
   refinement or regularization.

   @return the index of the new restraint.

   @return -1 when the atoms were not found and no extra bond
   restraint was stored.  */

int add_extra_geman_mcclure_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2, double bond_dist, double esd);
#ifdef __cplusplus
#ifdef USE_GUILE
int add_extra_bond_restraints_scm(int imol, SCM extra_bond_restraints_scm);
#endif // USE_GUILE
#ifdef USE_PYTHON
int add_extra_bond_restraints_py(int imol, PyObject *extra_bond_restraints_py);
#endif // USE_GUILE
#endif

void set_show_extra_distance_restraints(short int state);

int add_extra_angle_restraint(int imol,
				const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1,
				const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2,
				const char *chain_id_3, int res_no_3, const char *ins_code_3, const char *atom_name_3, const char *alt_conf_3,
				double torsion_angle, double esd);
int add_extra_torsion_restraint(int imol,
				const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1,
				const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2,
				const char *chain_id_3, int res_no_3, const char *ins_code_3, const char *atom_name_3, const char *alt_conf_3,
				const char *chain_id_4, int res_no_4, const char *ins_code_4, const char *atom_name_4, const char *alt_conf_4,
				double torsion_angle, double esd, int period);
int add_extra_start_pos_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, double esd);

int add_extra_target_position_restraint(int imol,
					const char *chain_id,
					int res_no,
					const char *ins_code,
					const char *atom_name,
 					const char *alt_conf, float x, float y, float z, float weight);

/*! \brief clear out all the extra/user-defined restraints for molecule number imol  */
void delete_all_extra_restraints(int imol);

/*! \brief clear out all the extra/user-defined restraints for this residue in molecule number imol  */
void delete_extra_restraints_for_residue(int imol, const char *chain_id, int res_no, const char *ins_code);

#ifdef __cplusplus
#ifdef USE_GUILE
void delete_extra_restraints_for_residue_spec_scm(int imol, SCM residue_spec_in);
#endif // USE_GUILE
#ifdef USE_PYTHON
void delete_extra_restraints_for_residue_spec_py(int imol, PyObject *residue_spec_in_py);
#endif // USE_PYTHON
#endif // __cplusplus

void delete_extra_restraints_worse_than(int imol, float n_sigma);

/*! read in prosmart (typically) extra restraints */
void add_refmac_extra_restraints(int imol, const char *file_name);

void set_show_extra_restraints(int imol, int state);
int extra_restraints_are_shown(int imol);

/*! \brief often we don't want to see all prosmart restraints, just the (big) violations */
void set_extra_restraints_prosmart_sigma_limits(int imol, double limit_high, double limit_low);

/*! \brief generate external distance local self restraints */
void generate_local_self_restraints(int imol, const char *chain_id, float local_dist_max);

/*! \brief generate external distance all-molecule self restraints */
void generate_self_restraints(int imol, float local_dist_max);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief generate external distance self restraints for selected residues */
void generate_local_self_restraints_by_residues_scm(int imol, SCM residue_specs, float local_dist_max);
#endif // USE_GUILE
#ifdef USE_PYTHON
void generate_local_self_restraints_by_residues_py(int imol, PyObject *residue_specs, float local_dist_max);
#endif // USE_PYTHON
#endif // __cplusplus


/*! \brief proSMART interpolated restraints for model morphing  */
void write_interpolated_extra_restraints(int imol_1, int imol_2, int n_steps, const char *file_name_stub);

/*! \brief proSMART interpolated restraints for model morphing and write interpolated model

interpolation_mode is currently dummy - in due course I will addd torion angle interpolation.
*/
void write_interpolated_models_and_extra_restraints(int imol_1, int imol_2, int n_steps, const char *file_name_stub,
						    int interpolation_mode);

void set_show_parallel_plane_restraints(int imol, int state);
int parallel_plane_restraints_are_shown(int imol);
void add_parallel_plane_restraint(int imol,
				  const char *chain_id_1, int re_no_1, const char *ins_code_1,
				  const char *chain_id_2, int re_no_2, const char *ins_code_2);
void set_extra_restraints_representation_for_bonds_go_to_CA(int imol, short int state);


#ifdef __cplusplus
#ifdef USE_GUILE
/* restraint_spec is something like (list 'bond spec-1 spec-2)

  spec-1 and spec-2 do not have to be in the order that the bond was created.  */
void delete_extra_restraint_scm(int imol, SCM restraint_spec);
SCM list_extra_restraints_scm(int imol);
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
/* restraint_spec is something like ['bond', spec_1, spec_2]

  spec_1 and spec_2 do not have to be in the order that the bond was created.  */
void delete_extra_restraint_py(int imol, PyObject *restraint_spec);
PyObject *list_extra_restraints_py(int imol);
#endif /* USE_PYTHON */
#endif /*  __cplusplus */

/*! \brief set use only extra torsion restraints for torsions */
void set_use_only_extra_torsion_restraints_for_torsions(short int state);
/*! \brief return only-use-extra-torsion-restraints-for-torsions state */
int use_only_extra_torsion_restraints_for_torsions_state();

void clear_all_atom_pull_restraints();

/*! \brief set auto-clear atom pull restraint */
void set_auto_clear_atom_pull_restraint(int state);

/*! \brief get auto-clear atom pull restraint state */
int  get_auto_clear_atom_pull_restraint_state();

/*! \brief iscrease the proportional editing radius*/
void increase_proportional_editing_radius();

/*! \brief descrease the proportional editing radius*/
void decrease_proportional_editing_radius();


/*  ----------------------------------------------------------------------- */
/*                  Restraints editor                                       */
/*  ----------------------------------------------------------------------- */

/*! \} */

/*  ----------------------------------------------------------------------- */
/*               Simplex Refinement                                         */
/*  ----------------------------------------------------------------------- */
/*! \name Simplex Refinement Interface */
/*! \{ */

/*! \brief refine residue range using simplex optimization */
void
fit_residue_range_to_map_by_simplex(int res1, int res2, const char *altloc, const char *chain_id, int imol, int imol_for_map);

/*! \brief simply score the residue range fit to map */
float
score_residue_range_fit_to_map(int res1, int res2, const char *altloc, const char *chain_id, int imol, int imol_for_map);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*               Nomenclature Errors                                        */
/*  ----------------------------------------------------------------------- */
/*! \name Nomenclature Errors */
/*! \{ */
/*! \brief fix nomenclature errors in molecule number imol

   @return the number of resides altered. */
int fix_nomenclature_errors(int imol);

/*! \brief set way nomenclature errors should be handled on reading
  coordinates.

  mode should be "auto-correct", "ignore", "prompt".  The
  default is "prompt" */
void set_nomenclature_errors_on_read(const char *mode);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*               Atom info                                                  */
/*  ----------------------------------------------------------------------- */
/* section Atom Info Interface */
/*! \name Atom Info  Interface */
/*! \{ */

/*! \brief output to the terminal the Atom Info for the give atom specs
 */
void
output_atom_info_as_text(int imol, const char *chain_id, int resno,
			 const char *ins_code, const char *atname,
			 const char *altconf);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*               (Eleanor's) Residue info                                   */
/*  ----------------------------------------------------------------------- */
/* section Residue Info */
/*! \name Residue Info */
/*! \{ */
/* Similar to above, we need only one click though. */
void do_residue_info_dialog();

/* MOVE-ME to c-interface-gtk-widgets.h */
 void output_residue_info_dialog    (int imol, int atom_index); /* widget version */
/* scripting version */
/*! \brief show residue info dialog for given residue */
void residue_info_dialog(int imol, const char *chain_id, int resno, const char *ins_code);
int residue_info_dialog_is_displayed();
void output_residue_info_as_text(int atom_index, int imol); /* text version */
/* functions that uses mmdb_manager functions/data types moved to graphics_info_t */

void do_distance_define();
void do_angle_define();
void do_torsion_define();
void residue_info_apply_all_checkbutton_toggled();
void clear_residue_info_edit_list();

/* a graphics_info_t function wrapper: */
void unset_residue_info_widget();
void clear_measure_distances();
void clear_last_measure_distance();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*               GUI edit functions                                   */
/*  ----------------------------------------------------------------------- */
/* section Edit Fuctions */
/*! \name Edit Fuctions */
/*! \{ */

void  do_edit_copy_molecule();
void  do_edit_copy_fragment();
void  do_edit_replace_fragment();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  residue environment                                     */
/*  ----------------------------------------------------------------------- */
/* section Residue Environment Functions */
/*! \name Residue Environment Functions */
/*! \{ */
/*! \brief show environment distances.  If state is 0, distances are
  turned off, otherwise distances are turned on. */
void set_show_environment_distances(int state);
/*! \brief show bumps environment distances.  If state is 0, bump distances are
  turned off, otherwise bump distances are turned on. */
void set_show_environment_distances_bumps(int state);
/*! \brief show H-bond environment distances.  If state is 0, bump distances are
  turned off, otherwise H-bond distances are turned on. */
void set_show_environment_distances_h_bonds(int state);
/*! \brief show the state of display of the  environment distances  */
int show_environment_distances_state();
/*! \brief min and max distances for the environment distances */
void set_environment_distances_distance_limits(float min_dist, float max_dist);

/*! \brief show the environment distances with solid modelling */
void set_show_environment_distances_as_solid(int state);

/*! \brief Label the atom on Environment Distances start/change */
void set_environment_distances_label_atom(int state);

/*! \brief Label the atoms in the residues around the central residue */
void label_neighbours();

/*! \brief Label the atoms in the central residue */
void label_atoms_in_residue();

/*! \brief Label the atoms with their B-factors */
void set_show_local_b_factors(short int state);

/*! \brief Add a geometry distance between points in a given molecule

@return the distance between the points

*/
double add_geometry_distance(int imol_1, float x_1, float y_1, float z_1, int imol_2, float x_2, float y_2, float z_2);
#ifdef __cplusplus
#ifdef USE_GUILE
double add_atom_geometry_distance_scm(int imol_1, SCM atom_spec_1, int imol_2, SCM atom_spec_2);
#endif
#ifdef USE_PYTHON
double add_atom_geometry_distance_py(int imol_1, PyObject *atom_spec_1, int imol_2, PyObject *atom_spec_2);
#endif
#endif /* __cplusplus */

/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  pointer position                                        */
/*  ----------------------------------------------------------------------- */
/* section Pointer Position Function */
/*! \name Pointer Position Function */
/*! \{ */
/*! \brief return the [x,y] position of the pointer in fractional coordinates.

the origin is top-left.
may return false if pointer is not available */
#ifdef __cplusplus
#ifdef USE_PYTHON
PyObject *get_pointer_position_frac_py();
#endif // USE_PYTHON
#endif	/* c++ */
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  pointer distances                                      */
/*  ----------------------------------------------------------------------- */
/* section Pointer Functions */
/*! \name Pointer Functions */
/*! \{ */
/*! \brief turn on (or off) the pointer distance by passing 1 (or 0). */
void set_show_pointer_distances(int istate);
/*! \brief show the state of display of the  pointer distances  */
int  show_pointer_distances_state();
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  zoom                                                    */
/*  ----------------------------------------------------------------------- */
/* section Zoom Functions */
/*! \name Zoom Functions */
/*! \{ */
/*! \brief scale the view by f

   Values outside the range 0.5 to 1.8 have no effect.
   external (scripting) interface (with redraw)
    @param f the smaller f, the bigger the zoom, typical value 1.3.
    */
void scale_zoom(float f);
/* internal interface */
void scale_zoom_internal(float f);
/*! \brief return the current zoom factor i.e. get_zoom_factor() */
float zoom_factor();

/*! \brief set smooth scroll with zoom
   @param i 0 means no, 1 means yes: (default 0) */
void set_smooth_scroll_do_zoom(int i);
/* default 1 (on) */
/*! \brief return the state of the above system */
int      smooth_scroll_do_zoom();
float    smooth_scroll_zoom_limit();
void set_smooth_scroll_zoom_limit(float f);

/*! \brief set the zoom factor (absolute value) - maybe should be called set_zoom_factor() */
void set_zoom(float f);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  CNS data stuff                                          */
/*  ----------------------------------------------------------------------- */
/*! \name CNS Data Functions */
/*! \{ */
/*! \brief read CNS data (currently only a placeholder)  */
int handle_cns_data_file(const char *filename, int imol);

/*! \brief read CNS data (currently only a placeholder)

a, b,c are in Angstroems.  alpha, beta, gamma are in degrees.  spg is
the space group info, either ;-delimited symmetry operators or the
space group name*/
int handle_cns_data_file_with_cell(const char *filename, int imol, float a, float b, float c, float alpha, float beta, float gamma, const char *spg_info);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  cif stuff                                               */
/*  ----------------------------------------------------------------------- */
/* section mmCIF Functions */
/*! \name mmCIF Functions */
/* dataset stuff */
/*! \{ */
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

int write_connectivity(const char* monomer_name, const char *filename);
/*! \brief open the cif dictionary file selector dialog */
void open_cif_dictionary_file_selector_dialog();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief non-standard residue/monomer names (note HOH is not non-standard). */
SCM non_standard_residue_names_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief non-standard residue/monomer names (note HOH is not non-standard). */
PyObject *non_standard_residue_names_py(int imol);
#endif /* USE_PYTHON */
#endif /* c++ */

/* Use the environment variable COOT_REFMAC_LIB_DIR to find cif files
   in subdirectories and import them all. */
void import_all_refmac_cifs();

int read_small_molecule_cif(const char *file_name);

int read_small_molecule_data_cif(const char *file_name);

int read_small_molecule_data_cif_and_make_map_using_coords(const char *file_name,
							   int imol_coords);

/*! \} */
/*  ------------------------------------------------------------------------ */
/*                         Validation:                                       */
/*  ------------------------------------------------------------------------ */
/* section Validation Functions */
/*! \name Validation Functions */
/*! \{ */
void deviant_geometry(int imol);
short int is_valid_model_molecule(int imol);
short int is_valid_map_molecule(int imol);


/*! \brief generate a list of difference map peaks

peaks within max_closeness (2.0 A typically) of a larger peak are not
listed.

the flag around_model_only_flag limits the peak list to those only within 4A
of the selected model (useful for maps with molecular symmetry).

*/
void difference_map_peaks(int imol, int imol_coords, float level, float max_closeness, int do_positive_level_flag, int do_negative_level_flag, int around_model_only_flag);

/*! \brief set the max closeness (i.e. no smaller peaks can be within
   max_closeness of a larger peak)

In the GUI for difference map peaks, there is not a means to set the
max_closeness, so here is a means to set it and query it. */
void set_difference_map_peaks_max_closeness(float m);
float difference_map_peaks_max_closeness();

void clear_diff_map_peaks();

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
PyObject *map_peaks_around_molecule_py(int imol_map, float sigma, int negative_also_flag, int imol_coords);

/* BL says:: this probably shouldnt be here but cluster with KK code */
PyObject *screen_vectors_py();
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

/*! \brief set the given torsion the given residue. tors is in
  degrees.  Return the resulting torsion (also in degrees). */
SCM set_torsion_scm(int imol, const char *chain_id, int res_no, const char *insertion_code,
		    const char *alt_conf,
		    const char *atom_name_1,
		    const char *atom_name_2,
		    const char *atom_name_3,
		    const char *atom_name_4, double tors);

/*! \brief create a multi-residue torsion dialog (user manipulation of torsions) */
void multi_residue_torsion_scm(int imol, SCM residues_specs_scm);


#endif  /* USE_GUILE */


#ifdef USE_PYTHON
PyObject *get_torsion_py(int imol, PyObject *atom_spec_1, PyObject *atom_spec_2, PyObject *atom_spec_3, PyObject *atom_spec_4);

/*! \brief set the given torsion the given residue. tors is in
  degrees.  Return the resulting torsion (also in degrees). */
PyObject *set_torsion_py(int imol, const char *chain_id, int res_no, const char *insertion_code,
		         const char *alt_conf,
		         const char *atom_name_1,
		         const char *atom_name_2,
		         const char *atom_name_3,
		         const char *atom_name_4, double tors);

/*! \brief create a multi-residue torsion dialog (user manipulation of torsions) */
void multi_residue_torsion_py(int imol, PyObject *residues_specs_py);

#endif  /* USE_PYTHON */


#endif /* __cplusplus  */

/* These functions are called from callbacks.c */
void clear_multi_residue_torsion_mode();
void set_multi_residue_torsion_reverse_mode(short int mode);
void show_multi_residue_torsion_dialog(); /* show the rotatable bonds dialog */
void setup_multi_residue_torsion();  /* show the pick dialog */

/*! \brief return the atom overlap score */
float atom_overlap_score(int imol);


/*! \brief set the state of showing chiral volume outlier markers - of a model molecule that is,
   not the intermediate atoms (derived from restraints) */
void set_show_chiral_volume_outliers(int imol, int state);

/*! \brief set the state of showing non-bonded contact markers - of a model molecule that is,
   not the intermediate atoms (derived from restraints) */
void set_show_chiral_volume_outliers(int imol, int state);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  ramachandran plot                                       */
/*  ----------------------------------------------------------------------- */
/* section Ramachandran Plot Functions */
/*! \name Ramachandran Plot Functions */
/*! \{ */
/* Note for optionmenu from the main window menubar, we should use
   code like this, rather than the map_colour/attach_scroll_wheel code
   (actually, they are mostly the same, differing only the container
   delete code). */

/*! \brief Ramachandran plot for molecule number imol */
void do_ramachandran_plot(int imol);

/*! \brief set the number of biggest difference arrows on the Kleywegt
  plot.  */
void set_kleywegt_plot_n_diffs(int n_diffs);

/*! \brief set the contour levels for the ramachandran plot, default
  values are 0.02 (prefered) 0.002 (allowed) */
void set_ramachandran_plot_contour_levels(float level_prefered, float level_allowed);
/*! \brief set the ramachandran plot background block size.

  Smaller is smoother but slower.  Should be divisible exactly into
  360.  Default value is 10. */
void set_ramachandran_plot_background_block_size(float blocksize) ;

/*! \brief set the psi axis for the ramachandran plot. Default (0) from -180
 to 180. Alternative (1) from -120 to 240.  */
void set_ramachandran_psi_axis_mode(int mode);
int ramachandran_psi_axis_mode();

void set_moving_atoms(double phi, double psi);

/*! \brief this does the same as `accept_moving_atoms()`
*/
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
/*! \} */

/*  ----------------------------------------------------------------------- */
/*           sequence_view                                                  */
/*  ----------------------------------------------------------------------- */
/*! \name Sequence View Interface  */
/*! \{ */
/*! \brief display the sequence view dialog for molecule number imol */
void sequence_view(int imol);

/*! \brief old name for the above function */
void do_sequence_view(int imol);

/*!  \brief update the sequnce view current position highlight based on active atom */
void update_sequence_view_current_position_highlight_from_active_atom();

/*! \} */

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
/*! \{ */
/*  Note we have to search for " CA " etc */
int    add_atom_label(int imol, const char *chain_id, int iresno, const char *atom_id);
int remove_atom_label(int imol, const char *chain_id, int iresno, const char *atom_id);
void remove_all_atom_labels();

void set_label_on_recentre_flag(int i); /* 0 for off, 1 or on */

int centre_atom_label_status();

/*! \brief use brief atom names for on-screen labels

 call with istat=1 to use brief labels, istat=0 for normal labels */
void set_brief_atom_labels(int istat);

/*! \brief the brief atom label state */
int brief_atom_labels_state();

/*! \brief set if brief atom labels should have seg-ids also */
void set_seg_ids_in_atom_labels(int istat);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  scene rotation                                          */
/*  ----------------------------------------------------------------------- */
/* section Screen Rotation */
/*! \name Screen Rotation */
/* stepsize in degrees */
/*! \{ */
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

/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  scene rotation                                          */
/*  ----------------------------------------------------------------------- */
/* section Screen Translation */
/*! \name  Screen Translation */
/*! \{ */
/*! \brief translate rotation centre relative to screen axes for nsteps */
void translate_scene_x(int nsteps);
/*! \brief translate rotation centre relative to screen axes for nsteps */
void translate_scene_y(int nsteps);
/*! \brief translate rotation centre relative to screen axes for nsteps */
void translate_scene_z(int nsteps);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Views                                                   */
/*  ----------------------------------------------------------------------- */
/* section Views Interface */
/*! \name Views Interface */
/*! \{ */
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
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  movies                                                  */
/*  ----------------------------------------------------------------------- */
/* movies */
/* section Movies Interface */
/*! \name Movies Interface */
/*! \{ */
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
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  graphics background colour                              */
/*  ----------------------------------------------------------------------- */
/* section Background Colour */
/*! \name Background Colour */
/*! \{ */

/*! \brief set the background colour

 red, green and blue are numbers between 0.0 and 1.0 */
void set_background_colour(double red, double green, double blue);

/*! \brief re draw the background colour when switching between mono and stereo */
void redraw_background();

/*! \brief is the background black (or nearly black)?

@return 1 if the background is black (or nearly black),
else return 0. */
int  background_is_black_p();
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  ligand fitting stuff                                    */
/*  ----------------------------------------------------------------------- */
/* section Ligand Fitting Functions */
/*! \name   Ligand Fitting Functions */
/*! \{ */

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

void set_find_ligand_do_real_space_refinement(short int state);

/*! \brief allow multiple ligand solutions per cluster.

The first limit is the fraction of the top scored positions that go on
to correlation scoring (closer to 1 means less and faster - default
0.7).

The second limit is the fraction of the top correlation score that is
considered interesting.  Limits the number of solutions displayed to
user. Default 0.9.

There is currently no chi-angle set redundancy filtering - I suspect
that there should be.

Nino-mode.

*/
void set_find_ligand_multi_solutions_per_cluster(float lim_1, float lim_2);

/*! \brief how shall we treat the waters during ligand fitting?

   pass with istate=1 for waters to mask the map in the same way that
   protein atoms do.
   */
void set_find_ligand_mask_waters(int istate);

/* get which map to search, protein mask and ligands from button and
   then do it*/

/*  extract the sigma level and stick it in */
/*  graphics_info_t::ligand_cluster_sigma_level */

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

/*! \brief  Allow the user a scripting means to find ligand at the rotation centre */
void set_find_ligand_here_cluster(int state);

void execute_ligand_search();

#ifdef __cplusplus
#ifdef USE_GUILE
SCM execute_ligand_search_scm();
#endif
#ifdef USE_PYTHON
PyObject *execute_ligand_search_py();
#endif /* USE_PYTHON */
#endif /* __cplusplus */
void add_ligand_clear_ligands();

/* conformers added to cc-interface because it uses a std::vector internally.  */


/*! \brief this sets the flag to have expert option ligand entries in
  the Ligand Searching dialog */
void ligand_expert();

/*! \brief display the find ligands dialog

   if maps, coords and ligands are available, that is.
*/
void do_find_ligands_dialog();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief Overlap residue with "template"-based matching.

  Overlap the first residue in
  imol_ligand onto the residue specified by the reference parameters.
  Use graph matching, not atom names.

@return success status, False = failed to find residue in either
imol_ligand or imo_ref.  If success, return the RT operator.
*/
SCM overlap_ligands(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
SCM analyse_ligand_differences(int imol_ligand, int imol_ref, const char *chain_id_ref,
			       int resno_ref);
SCM compare_ligand_atom_types_scm(int imol_ligand, int imol_ref, const char *chain_id_ref,
				  int resno_ref);
#endif /* USE_GUILE */
void match_ligand_torsions(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
#ifdef USE_PYTHON
PyObject *overlap_ligands_py(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
PyObject *analyse_ligand_differences_py(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
PyObject *compare_ligand_atom_types_py(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
#endif /* PYTHON*/
#endif	/* __cplusplus */

/*! \brief Match ligand atom names

  By using graph matching, make the names of the atoms of the
  given ligand/residue match those of the reference residue/ligand as
  closely as possible - where there would be an atom name clash, invent
  a new atom name.
 */
void match_ligand_atom_names(int imol_ligand,
			     const char *chain_id_ligand, int resno_ligand, const char *ins_code_ligand,
			     int imol_reference, const char *chain_id_reference,
			     int resno_reference, const char *ins_code_reference);

/*! \brief Match ligand atom names to a reference ligand type (comp_id)

  By using graph matching, make the names of the atoms of the
  given ligand/residue match those of the reference ligand from the
  geometry store as closely as possible. Where there would be an
  atom name clash, invent a new atom name.

  This doesn't create a new dictionary for the selected ligand -
  and that's a big problem (see match_residue_and_dictionary).
 */
void match_ligand_atom_names_to_comp_id(int imol_ligand, const char *chain_id_ligand, int resno_ligand, const char *ins_code_ligand, const char *comp_id_ref);

/* Transfer as many atom names as possible from the reference ligand
   to the given ligand.  The atom names are determined from graph
   matching the reference ligand onto the given ligand.

   Function needs to be written.  Non-trivial (the atom graph matching
   is OK, and the atom pairs straightforwardly determined, but what
   should be done with the atom names that are matched from the
   reference ligand, but also a different atom of the same name is not
   matched?).
 */
/* void tranfer_atom_names(int imol_ligand, const char *chain_id_ligand, int res_no_ligand, const char *ins_code_ligand, */
/* 			int imol_reference, const char *chain_id_reference, int res_no_reference, const char *ins_code_reference); */



/* Just pondering - just a stub currently.

   For use with exporting ligands from the 2D sketcher to the main
   coot window. It should find the residue that this residue is
   sitting on top of that is in a molecule that has lot of atoms
   (i.e. is a protein) and create a new molecule that is a copy of the
   molecule without the residue/ligand that this (given) ligand
   overlays - and a copy of this given ligand.  */
int exchange_ligand(int imol_lig, const char *chain_id_lig, int resno_lig, const char *ins_code_lig);



/*! \brief flip the ligand (usually active residue) around its eigen vectors
   to the next flip number.  Immediate replacement (like flip
   peptide). */
void flip_ligand(int imol, const char *chain_id, int resno);

void jed_flip(int imol, const char *chain_id, int res_no, const char *ins_code, const char *atom_name, const char *alt_conf, short int invert_selection);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  water fitting                                           */
/*  ----------------------------------------------------------------------- */
/* section Water Fitting Functions */
/*! \name Water Fitting Functions */
/*! \{ */

/*! \brief create a dialog for water fitting */
void show_create_find_waters_dialog();

/*! \brief Renumber the waters of molecule number imol with consecutive numbering */
void renumber_waters(int imol);

/*! \brief find waters */
void execute_find_waters_real(int imol_for_map,
			      int imol_for_protein,
			      short int new_waters_mol_flag,
			      float rmsd_cut_off);

void find_waters(int imol_for_map,
		 int imol_for_protein,
		 short int new_waters_mol_flag,
		 float rmsd_cut_off,
		 short int show_blobs_dialog);


/*! \brief move waters of molecule number imol so that they are around the protein.

@return the number of moved waters. */
int move_waters_to_around_protein(int imol);

/*! \brief move all hetgroups (including waters) of molecule number
  imol so that they are around the protein.
 */
void move_hetgroups_to_around_protein(int imol);

/*! \brief return the maximum minimum distance of any water atom to
  any protein atom - used in validation of
  move_waters_to_around_protein() funtion.*/
float max_water_distance(int imol);

char *get_text_for_find_waters_sigma_cut_off();
void set_value_for_find_waters_sigma_cut_off(float f);

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

/*! \brief find blobs
 *
 * Not useful for MCP. For interactive use only.
 *
 * */
void execute_find_blobs(int imol_model, int imol_for_map, float cut_off, short int interactive_flag);

/* there is also a c++ interface to find blobs, which returns a vector
   of pairs (currently) */

/*! \brief split the given water and fit to map.

If refinement map is not defined, don't do anything.

If there is more than one atom in the specified resiue, don't do
anything.

If the given atom does not have an alt conf of "", don't do anything.

 @param imol the index of the molecule
 @param chain_id the chain id
 @param res_no the residue number
 @param ins_code the insertion code of the residue

 */
void split_water(int imol, const char *chain_id, int res_no, const char *ins_code);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  bond representation                                     */
/*  ----------------------------------------------------------------------- */
/* section Bond Representation */
/*! \name Bond Representation */
/*! \{ */

/*! \brief set the default thickness for bonds (e.g. in ~/.coot)  */
void set_default_bond_thickness(int t);

/*! \brief set the thickness of the bonds in molecule number imol to t pixels  */
void set_bond_thickness(int imol, float t);
/*! \brief set the thickness of the bonds of the intermediate atoms to t pixels  */
void set_bond_thickness_intermediate_atoms(float t);

/*! \brief allow lines that are further away to be thinner */
void set_use_variable_bond_thickness(short int state);

/*! \brief set bond colour for molecule */
void set_bond_colour_rotation_for_molecule(int imol, float f);

/*! \brief set default for the drawing of atoms in stick mode (default is on (1)) */
void set_draw_stick_mode_atoms_default(short int state);


/*! \brief get the bond colour for molecule.

Return -1 on err (bad molecule number) */
float get_bond_colour_rotation_for_molecule(int imol);

void set_unbonded_atom_star_size(float f);

/*! \brief set the default representation type (default 1).*/
void set_default_representation_type(int type);

/*! \brief get the default thickness for bonds*/
int get_default_bond_thickness();

/*! \brief set status of drawing zero occupancy markers.

  default status is 1. */
void set_draw_zero_occ_markers(int status);


/*! \brief set status of drawing cis-peptide markups

  default status is 1. */
void set_draw_cis_peptide_markups(int status);



/*! \brief set the hydrogen drawing state. istat = 0 is hydrogens off,
  istat = 1: show hydrogens */
void set_draw_hydrogens(int imol, int istat);

/*! \brief the state of draw hydrogens for molecule number imol.

return -1 on bad imol.  */
int draw_hydrogens_state(int imol);

/*! \brief draw little coloured balls on atoms

turn off with state = 0

turn on with state = 1 */
void set_draw_stick_mode_atoms(int imol, short int state);

/*! \brief set the state for drawing missing resiude loops

For taking screenshots, we often don't want to see them.
*/
void set_draw_missing_residues_loops(short int state);

/*! \brief draw molecule number imol as CAs */
void graphics_to_ca_representation   (int imol);
/*! \brief draw molecule number imol coloured by chain */
void graphics_to_colour_by_chain(int imol);
/*! \brief draw molecule number imol as CA + ligands */
void graphics_to_ca_plus_ligands_representation   (int imol);
/*! \brief draw molecule number imol as CA + ligands + sidechains*/
void graphics_to_ca_plus_ligands_and_sidechains_representation   (int imol);
/*! \brief draw molecule number imol with no waters */
void graphics_to_bonds_no_waters_representation(int imol);
/*! \brief draw molecule number imol with normal bonds */
void graphics_to_bonds_representation(int mol);
/*! \brief draw molecule with colour-by-molecule colours */
void graphics_to_colour_by_molecule(int imol);
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
/*! \brief draw molecule number imol in CA+Ligands mode coloured by user-defined atom colours */
void graphics_to_user_defined_atom_colours_representation(int imol);
/*! \brief draw molecule number imol all atoms coloured by user-defined atom colours */
void graphics_to_user_defined_atom_colours_all_atoms_representation(int imol);
/*! \brief what is the bond drawing state of molecule number imol  */
int get_graphics_molecule_bond_type(int imol);
/*! \brief scale the colours for colour by b factor representation */
int set_b_factor_bonds_scale_factor(int imol, float f);
/*! \brief change the representation of the model molecule closest to
  the centre of the screen */
void change_model_molecule_representation_mode(int up_or_down);

/* not today void set_ca_bonds_loop_params(float p1, float p2, float p3); */

/*! \brief make the carbon atoms for molecule imol be grey
 */
void set_use_grey_carbons_for_molecule(int imol, short int state);
/*! \brief set the colour for the carbon atoms

can be not grey if you desire, r, g, b in the range 0 to 1.
 */
void set_grey_carbon_colour(int imol, float r, float g, float b);

/* undocumented feature for development. */
void set_draw_moving_atoms_restraints(int state);

/* undocumented feature for development. */
short int get_draw_moving_atoms_restraints();

/*! \brief make a ball and stick representation of imol given atom selection

e.g. (make-ball-and-stick 0 "/1" 0.15 0.25 1) */
int make_ball_and_stick(int imol,
			const char *atom_selection_str,
			float bond_thickness, float sphere_size,
			int do_spheres_flag);
/*! \brief clear ball and stick representation of molecule number imol */
int clear_ball_and_stick(int imol);

/*! \brief set the model molecule representation stye 0 for ball-and-stick/licorice (default) and 1 for ball */
void set_model_molecule_representation_style(int imol, unsigned int mode);

/*! \brief set show a ribbon/mesh for a given molecule */
void set_show_molecular_representation(int imol, int mesh_index, short int state);

/* removed from API brief display/undisplay the given additional representation  */
void set_show_additional_representation(int imol, int representation_number, int on_off_flag);

/*! \brief display/undisplay all the additional representations for the given molecule  */
void set_show_all_additional_representations(int imol, int on_off_flag);

/*! removed from API brief undisplay all the additional representations for the given
   molecule, except the given representation number (if it is off, leave it off)  */
void all_additional_representations_off_except(int imol, int representation_number,
					       short int ball_and_sticks_off_too_flag);

/*! removed from API brief delete a given additional representation */
void delete_additional_representation(int imol, int representation_number);

/*! removed from API brief return the index of the additional representation.  Return -1 on error */
int additional_representation_by_string(int imol,  const char *atom_selection,
					int representation_type,
					int bonds_box_type,
					float bond_width,
					int draw_hydrogens_flag);

/*   representation_types: */
/*   enum { coot::SIMPLE_LINES, coot::STICKS, coot::BALL_AND_STICK, coot::SURFACE };

  bonds_box_type:
  enum {  UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2, COLOUR_BY_CHAIN_BONDS=3,
	  CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
	  BONDS_NO_HYDROGENS=15,
	  CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
	  CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
	  COLOUR_BY_MOLECULE_BONDS=8,
	  COLOUR_BY_RAINBOW_BONDS=9, COLOUR_BY_B_FACTOR_BONDS=10,
	  COLOUR_BY_OCCUPANCY_BONDS=11};

*/

/*! \brief return the index of the additional representation.
  @return -1 on error.
 */
int additional_representation_by_attributes(int imol,  const char *chain_id,
					    int resno_start, int resno_end,
					    const char *ins_code,
					    int representation_type,
					    int bonds_box_type,
					    float bond_width,
					    int draw_hydrogens_flag);


#ifdef __cplusplus

#ifdef USE_GUILE
SCM additional_representation_info_scm(int imol);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *additional_representation_info_py(int imol);
#endif	/* USE_PYTHON */

#endif	/* __cplusplus */

/* Turn on nice animated ligand interaction display.

turn on with arg 1.

turn off with arg 0. */
void set_flev_idle_ligand_interactions(int state);

/* Toggle for animated ligand interaction display above */
void toggle_flev_idle_ligand_interactions();

void calculate_hydrogen_bonds(int imol);

void set_draw_hydrogen_bonds(int state);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  dots display                                            */
/*  ----------------------------------------------------------------------- */

/*! \name Dots Representation */
/*! \{ */
/*! \brief display dotted surface

return a generic objects handle (which can be used to remove later) */
int dots(int imol,
	 const char *atom_selection_str,
	 const char *dots_object_name,
	 float dot_density, float sphere_size_scale);


/*! \brief set the colour of the surface dots of the imol-th molecule
  to be the given single colour

  r,g,b are values between 0.0 and 1.0 */
void set_dots_colour(int imol, float r, float g, float b);

/*! \brief no longer set the dots of molecule imol to a single colour

i.e. go back to element-based colours. */
void unset_dots_colour(int imol);

/*! \brief clear dots in imol with dots_handle */
void clear_dots(int imol, int dots_handle);

/*! \brief clear the first dots object for imol with given name */
void clear_dots_by_name(int imol, const char *dots_object_name);

/*! \brief return the number of dots sets for molecule number imol */
int n_dots_sets(int imol);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  pepflip                                                 */
/*  ----------------------------------------------------------------------- */
/* section Pep-flip Interface */
/*! \name Pep-flip Interface */
/*! \{ */
void do_pepflip(short int state); /* sets up pepflip, ready for atom pick. */

/*! \brief pepflip (flip the peptide) of the given residue
 *
 *  Rotate the the carbonyl C and O atom of this residue and the N of the
 *  next residue around a vector between the two CA atoms by 180 degrees.
 *  This is often a useful modelling operation to create a different hypothesis
 *  about the orientation of the main-chain atoms - that can then be used
 *  for refinement. This can sometimes allow the model to be removed from
 *  local minima of backbone conformations.
 *
 *  @param imol is the index of the model molecule
 *  @param chain_id is the chain-id
 *  @param res_no is the residue number (the residue that has the C and O atoms)
 *  @param inscode the insertion code (typically "")
 *  @param altconf the altconf (typically "")
 *
 */
void pepflip(int imol, const char *chain_id, int resno, const char *inscode,
	     const char *altconf);

int pepflip_intermediate_atoms();
int pepflip_intermediate_atoms_other_peptide();

#ifdef __cplusplus
#ifdef USE_GUILE
SCM pepflip_using_difference_map_scm(int imol_coords, int imol_difference_map, float n_sigma);
#endif
#ifdef USE_PYTHON
PyObject *pepflip_using_difference_map_py(int imol_coords, int imol_difference_map, float n_sigma);
#endif
#endif
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  rigid body refinement                                   */
/*  ----------------------------------------------------------------------- */
/* section Rigid Body Refinement Interface */
/*! \name Rigid Body Refinement Interface */
/*! \{ */
/* a gui-based interface: setup the rigid body refinement.*/
void do_rigid_body_refine(short int state);	/* set up for atom picking */

/*! \brief setup rigid body refine zone

   where we set the atom selection
   holders according to the arguments and then call
   execute_rigid_body_refine() */
void rigid_body_refine_zone(int imol, const char *chain_id, int reso_start, int resno_end);

void
rigid_body_refine_by_atom_selection(int imol, const char *atom_selection_string);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief rigid body refine using residue ranges.  residue_ranges is
    a list of residue ranges.  A residue range is (list chain-id
    resno-start resno-end). */
SCM rigid_body_refine_by_residue_ranges_scm(int imol, SCM residue_ranges);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief rigid body refine using residue ranges.  residue_ranges is
    a list of residue ranges.  A residue range is [chain_id,
    resno_start, resno_end]. */
PyObject *
rigid_body_refine_by_residue_ranges_py(int imol, PyObject *residue_ranges);
#endif /* USE_PYTHON */
#endif /* __cplusplus */

void execute_rigid_body_refine(short int auto_range_flag); /* atom picking has happened.
				     Actually do it */


/*! \brief set rigid body fraction of atoms in positive density

 @param f in the range 0.0 -> 1.0 (default 0.75) */
void set_rigid_body_fit_acceptable_fit_fraction(float f);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  dynamic map                                             */
/*  ----------------------------------------------------------------------- */
/* section Dynamic Map */
/*! \name Dynamic Map */
/*! \{ */
void   toggle_dynamic_map_display_size();
void   toggle_dynamic_map_sampling();
/* scripting interface: */
void set_dynamic_map_size_display_on();
void set_dynamic_map_size_display_off();
int get_dynamic_map_size_display();
void set_dynamic_map_sampling_on();
void set_dynamic_map_sampling_off();
int get_dynamic_map_sampling();
void set_dynamic_map_zoom_offset(int i);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  build one residue by phi/psi search                     */
/*  ----------------------------------------------------------------------- */
/* section Add Terminal Residue Functions */
/*! \name Add Terminal Residue Functions */
/*! \{ */
void do_add_terminal_residue(short int state);
/*  execution of this is in graphics_info_t because it uses a mmdb::Residue */
/*  in the interface and we can't have that in c-interface.h */
/*  (compilation of coot_wrap_guile goes mad on inclusion of */
/*   mmdb_manager.h) */
void set_add_terminal_residue_n_phi_psi_trials(int n);
/* Add Terminal Residues actually build 2 residues, this allows us to
   see both residues - default is 0 (off). */
void set_add_terminal_residue_add_other_residue_flag(int i);
void set_add_terminal_residue_do_rigid_body_refine(short int v);
void set_terminal_residue_do_rigid_body_refine(short int v); /* remove this for 0.9, wraps above */
void set_add_terminal_residue_debug_trials(short int debug_state);
int add_terminal_residue_immediate_addition_state();

/*! \brief set immediate addition of terminal residue

call with i=1 for immediate addtion */
void set_add_terminal_residue_immediate_addition(int i);

/*! \brief Add a terminal residue
 
  Some text here that should be a detailed-description

   @param residue_type can be "auto" 
   @param immediate_add is recommended to be 1.
   @return 0 on failure, 1 on success
*/
int add_terminal_residue(int imol, const char *chain_id, int residue_number,
                          const char *residue_type, int immediate_add);

/*! \brief Add a residue to a chain or at the end of a fragment

  This can be used to fill a gap of one residue or to fill a gap
  of multiple residues by being called several times. Probably
  RSR refinement would be useful after each call to this function in
  such a case.

  @param imol the molecule index
  @param chain_id the chain ID
  @param residue_number the residue number (of the existing residue to attach to)
  @param residue_type the type for new residue, can be "auto"
  @param immediate_add is recommended to be 1

   @return 0 on failure, 1 on success
*/
int add_residue_by_map_fit(int imol, const char *chain_id, int residue_number,
                           const char *residue_type, int immediate_add);

/*! \brief Add a terminal nucleotide

No fitting is done
*/
int add_nucleotide(int imol, const char *chain_id, int res_no);


/*! \brief Add a terminal residue using given phi and psi angles

  @param imol the molecule index
  @param chain_id the chain ID
  @param residue_number the residue number (of the existing residue to attach to
  @param residue_type can be "auto"
  @param phi is phi in degrees
  @param psi is psi in degrees
  @return the success status, 0 on failure, 1 on success
 */
int add_terminal_residue_using_phi_psi(int imol, const char *chain_id, int res_no,
				       const char *residue_type, float phi, float psi);

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

/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  scripting a residue with atoms                          */
/*  ----------------------------------------------------------------------- */
/* section Add A Residue Functions */
/*! \name  Add A Residue Functions */
/*! \{ */
#ifdef __cplusplus

/*! \brief add a residue with atoms in scripting

  @return the number of atoms added
*/
#ifdef USE_PYTHON
int add_residue_with_atoms_py(int imol, PyObject *residue_spec, const std::string &res_name, PyObject *list_of_atoms);
#endif
#endif /* c++ */

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  delete residue                                          */
/*  ----------------------------------------------------------------------- */
/* section Delete Residues */
/*! \name Delete Residues */
/* in build */
/* by graphics */
/*! \{ */
void delete_atom_by_atom_index(int imol, int index, short int do_delete_dialog);
void delete_residue_by_atom_index(int imol, int index, short int do_delete_dialog);
void delete_residue_hydrogens_by_atom_index(int imol, int index, short int do_delete_dialog);

/*! \brief delete residue range */
void delete_residue_range(int imol, const char *chain_id, int resno_start, int end_resno);

/*! \brief delete residue
 *
 * @param imol the molecule index
 * @param chain_id the chain id
 * @param res_no the residue number
 * @param inscode the insertion code
 *
 * @return 0 on failure to delete, return 1 on residue successfully deleted
 *
 * */
int delete_residue(int imol, const char *chain_id, int res_no, const char *inscode);

/*! \brief delete residue with altconf  */
void delete_residue_with_full_spec(int imol, int imodel, const char *chain_id, int resno, const char *inscode, const char *altloc);
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief delete residues in the residue spec list */
void delete_residues_scm(int imol, SCM residue_specs_scm);
#endif
#ifdef USE_PYTHON
/*! \brief delete residues in the residue spec list */
void delete_residues_py(int imol, PyObject *residue_specs_py);
#endif
#endif	/* c++ */
/*! \brief delete hydrogen atoms in residue  */
void delete_residue_hydrogens(int imol, const char *chain_id, int resno, const char *inscode, const char *altloc);
/*! \brief delete atom in residue */
void delete_atom(int imol, const char *chain_id, int resno, const char *ins_code, const char *at_name, const char *altloc);
/*! \brief delete all atoms in residue that are not main chain or CB */
void delete_residue_sidechain(int imol, const char *chain_id, int resno, const char*ins_code,
			      short int do_delete_dialog);
/*! \brief delete all hydrogens in molecule,

   @return number of hydrogens deleted. */
int delete_hydrogen_atoms(int imol);

/*! \brief delete all hydrogens in molecule,

   @return number of hydrogens deleted. */
int delete_hydrogens(int imol);

/*! \brief delete all waters in molecule,

   @return number of waters deleted. */
int delete_waters(int imol);

void post_delete_item_dialog();



/* toggle callbacks */
void set_delete_atom_mode();
void set_delete_residue_mode();
void set_delete_residue_zone_mode();
void set_delete_residue_hydrogens_mode();
void set_delete_water_mode();
void set_delete_sidechain_mode();
void set_delete_sidechain_range_mode();
void set_delete_chain_mode();
short int delete_item_mode_is_atom_p(); /* (predicate) a boolean */
short int delete_item_mode_is_residue_p(); /* predicate again */
short int delete_item_mode_is_water_p();
short int delete_item_mode_is_sidechain_p();
short int delete_item_mode_is_sidechain_range_p();
short int delete_item_mode_is_chain_p();
void clear_pending_delete_item(); /* for when we cancel with picking an atom */


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  rotate/translate buttons                                */
/*  ----------------------------------------------------------------------- */
/* section Rotate/Translate Buttons */
/*  sets flag for atom selection clicks */
void do_rot_trans_setup(short int state);
void rot_trans_reset_previous();
void set_rotate_translate_zone_rotates_about_zone_centre(int istate);
void set_rot_trans_object_type(short int rt_type); /* zone, chain, mol */
int get_rot_trans_object_type();

/*  ----------------------------------------------------------------------- */
/*                  cis and trans info and conversion                       */
/*  ----------------------------------------------------------------------- */
void do_cis_trans_conversion_setup(int istate);
void cis_trans_convert(int imol, const char *chain_id, int resno, const char *altconf);

#ifdef __cplusplus	/* need this wrapper, else gmp.h problems in callback.c */
#ifdef USE_GUILE
/*! \brief return cis_peptide info for imol.

Return a SCM list object of (residue1 residue2 omega) */
SCM cis_peptides(int imol);
SCM twisted_trans_peptides(int imol);
#endif /* GUILE */
#ifdef USE_PYTHON
/*! \brief return cis_peptide info for imol.

Return a Python list object of [residue1, residue2, omega] */
PyObject *cis_peptides_py(int imol);
PyObject *twisted_trans_peptides_py(int imol);
#endif /* PYTHON */
#endif

/*! \brief cis-trans convert the active residue of the active atom in the
    inermediate atoms, and continue with the refinement  */
int cis_trans_convert_intermediate_atoms();


/*  ----------------------------------------------------------------------- */
/*                  db-main                                                 */
/*  ----------------------------------------------------------------------- */
/* section Mainchain Building Functions */

/*! \name Mainchain Building Functions */
/*! \{ */
void do_db_main(short int state);
/*! \brief CA -> mainchain conversion

direction is either "forwards" or "backwards"

See also the function below.

return the new molecule number */
int db_mainchain(int imol,
		 const char *chain_id,
		 int iresno_start,
		 int iresno_end,
		 const char *direction);

/*! \brief CA-Zone to Mainchain for a fragment based on the given residue.

Both directions are built. This is the modern interface.
 */
int db_mainchains_fragment(int imol, const char *chain_id, int res_no);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  close molecule                                          */
/*  ----------------------------------------------------------------------- */
/* section Close Molecule Functions */
/*! \name Close Molecule Functions */
/*! \{ */

/*! \brief close the molecule */
void close_molecule(int imol);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  rotamers                                                */
/*  ----------------------------------------------------------------------- */
/* section Rotamer Functions */
/*! \name Rotamer Functions */
/*! \{ */

/* functions defined in c-interface-build */

/*! \brief set the mode of rotamer search, options are (ROTAMERSEARCHAUTOMATIC),
  (ROTAMERSEARCHLOWRES) (aka. "backrub rotamers"),
  (ROTAMERSEARCHHIGHRES) (with rigid body fitting) */
void set_rotamer_search_mode(int mode);

/*! \ brief get the mode of rotamer search */
int rotamer_search_mode_state();

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
float auto_fit_best_rotamer(int imol_coords,
                            const char *chain_id,
                            int resno,
			    const char *insertion_code,
			    const char *altloc,
			    int imol_map, int clash_flag, float lowest_probability);

/*! auto-fit the rotamer for the active residue */
float auto_fit_rotamer_active_residue();

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
int set_residue_to_rotamer_number(int imol, const char *chain_id, int resno, const char *ins_code,
				  const char *alt_conf, int rotamer_number);

/*! \brief set the residue specified to the rotamer name specified.

(rotamer names are the Richardson rotamer names.)

return value is 0 if atoms were not moved (e.g. because rotamer-name was not know)
*/
int set_residue_to_rotamer_name(int imol, const char *chain_id, int resno, const char *ins_code,
				const char *alt_conf, const char *rotamer_name);

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

/*! \brief Fill amino acid residues

do backrub rotamer search for residues, but don't do refinement
*/
void simple_fill_partial_residues(int imol);

#ifdef __cplusplus
#ifdef USE_GUILE
SCM missing_atom_info_scm(int imol);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *missing_atom_info_py(int imol);
#endif /* USE_PYTHON */
#endif /* __cplusplus */



#ifdef __cplusplus	/* need this wrapper, else gmp.h problems in callback.c */
#ifdef USE_GUILE
/*! \brief Activate rotamer graph analysis for molecule number imol.

Return rotamer info - function used in testing.  */
SCM rotamer_graphs(int imol);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief Activate rotamer graph analysis for molecule number imol.

Return rotamer info - function used in testing.  */
PyObject *rotamer_graphs_py(int imol);
#endif /* USE_PYTHON */
#endif /* c++ */

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  180 degree flip                                         */
/*  ----------------------------------------------------------------------- */
/*! \name 180 Flip Side chain */
/*! \{ */

/*! \brief rotate 180 degrees around the last chi angle */
void do_180_degree_side_chain_flip(int imol, const char* chain_id, int resno,
				   const char *inscode, const char *altconf);

void setup_180_degree_flip(short int state);

/* ! \brief side-chain 180 flip the terminal chi angle of the residue of the active atom */
int side_chain_flip_180_intermediate_atoms();


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  mutate                                                  */
/*  ----------------------------------------------------------------------- */
/* section Mutate Functions */
/*! \name Mutate Functions */
/*! \{ */

/* c-interface-build */
void setup_mutate(short int state);
/*! \brief Mutate then fit to map

 that we have a map define is checked first */
void setup_mutate_auto_fit(short int state);

void do_mutation(const char *type, short int is_stub_flag);

/*! \brief display a dialog that allows the choice of residue type to which to mutate
 */
void mutate_active_residue();

/* auto-mutate stuff */
short int progressive_residues_in_chain_check(const char *chain_id, int imol);

/*! \brief mutate a given residue

target_res_type is a three-letter-code.

Return 1 on a good mutate. */
int mutate(int imol, const char *chain_id, int ires, const char *inscode,  const char *target_res_type);

/*! \brief mutate a base. return success status, 1 for a good mutate. */
int mutate_base(int imol, const char *chain_id, int res_no, const char *ins_code, const char *res_type);

/* push the residues along a bit

e.g. if nudge_by is 1, then the sidechain of residue 20 is moved up
onto what is currently residue 21.  The mainchain numbering and atoms is not changed.

@return 0 for failure to nudge (becauese not all the residues were in the range)
        and 1 for success.
*/
int nudge_residue_sequence(int imol, const char *chain_id, int res_no_range_start, int res_no_range_end, int nudge_by, short int nudge_residue_numbers_also);

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

   This  function doesnt make backups, but `mutate()` does - CHECKME
   Hence `mutate()` is for use as a "one-by-one" type and the following
   2 by wrappers that muate either a residue range or a whole chain.

   @param ires_ser is the serial number of the residue, not the seqnum
   @param chain_id is the chain-id
   @param imol is the index of the model molecule
   @param target_res_type is the single-letter-code for the target residue

   Note that the target_res_type is a char, not a string (or a char *).
   So from the scheme interface you'd use (for example) hash
   backslash A for ALA.

   @return 1 on success, 0 on failure

*/
int mutate_single_residue_by_serial_number(int ires_ser,
					   const char *chain_id,
					   int imol, char target_res_type);

/*!  \brief ires is the seqnum of the residue (conventional) */
int mutate_single_residue_by_seqno(int imol, const char *chain_id, int ires, const char *inscode,
				   char target_res_type);

/*! \brief mutate and auto-fit

(Move this and the above function into cc-interface.hh one day)
 */
int mutate_and_autofit_residue_range(int imol, const char *chain_id, int start_res_no, int stop_res_no,
                                     const char *sequence);

/* an internal function - not useful for scripting: */

void do_base_mutation(const char *type);

/*! \brief set a flag saying that the residue chosen by mutate or
  auto-fit mutate should only be added as a stub (mainchain + CB) */
void set_residue_type_chooser_stub_state(short int istat);

void handle_residue_type_chooser_entry_chose_type(const char *entry_text, short int stub_mode);


/*! \} */

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
SCM add_alt_conf_scm(int imol, const char *chain_id, int res_no, const char *ins_code,
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

void unset_add_alt_conf_dialog(); /* set the static dialog holder in
				     graphics info to NULL */
void unset_add_alt_conf_define(); /* turn off pending atom pick */
void altconf();			/* temporary debugging interface. */
void set_add_alt_conf_new_atoms_occupancy(float f); /* default 0.5 */
float get_add_alt_conf_new_atoms_occupancy();
void set_show_alt_conf_intermediate_atoms(int i);
int  show_alt_conf_intermediate_atoms_state();
void zero_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2);
void fill_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2);
void set_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2, float occ);
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
/*! \brief accept the baton tip position

  a prime candidate for a key binding */
void accept_baton_position();	/* put an atom at the tip */
/*! \brief move the baton tip position - another prime candidate for a key binding */
void baton_tip_try_another();
/*! \brief move the baton tip to the previous position*/
void baton_tip_previous();
/*! \brief shorten the baton length */
void shorten_baton();
/*! \brief lengthen the baton */
void lengthen_baton();
/*! \brief delete the most recently build CA position */
void baton_build_delete_last_residue();
/*! \brief set the parameters for the start of a new baton-built fragment. direction can either
     be "forwards" or "backwards" */
void set_baton_build_params(int istart_resno, const char *chain_id, const char *direction);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  post baton mode                                         */
/*  ----------------------------------------------------------------------- */
/* section Post-Baton Functions */
/* c-interface-build */
/*! \brief Reverse the direction of a the fragment of the clicked on
   atom/residue.

    A fragment is a consecutive range of residues -
   where there is a gap in the numbering, that marks breaks between
   fragments in a chain.  There also needs to be a distance break - if
   the CA of the next/previous residue is more than 5A away, that also
   marks a break. Throw away all atoms in fragment other than CAs.*/
void reverse_direction_of_fragment(int imol, const char *chain_id, int resno);
void setup_reverse_direction(short int i);


/*  ----------------------------------------------------------------------- */
/*                  terminal OXT atom                                       */
/*  ----------------------------------------------------------------------- */
/* section Terminal OXT Atom */
/*! \name Terminal OXT Atom */
/* c-interface-build */
/*! \{ */
short int add_OXT_to_residue(int imol, const char *chain_id, int reso, const char *insertion_code);

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
/*! \{ */
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
void set_moving_atom_move_chis();

/*! \brief display the edit chi angles gui for the given residue

 return a status of 0 if it failed to fined the residue,
 return a value of 1 if it worked. */
int edit_chi_angles(int imol, const char *chain_id, int resno,
		     const char *ins_code, const char *altconf);

int set_show_chi_angle_bond(int imode);

/* a callback from the callbacks.c, setting the state of
   graphics_info_t::edit_chi_angles_reverse_fragment flag */
void set_edit_chi_angles_reverse_fragment_state(short int istate);

/* No need for this to be exported to scripting */
/*! \brief beloved torsion general at last makes an entrance onto the
  Coot scene...  */
void setup_torsion_general(short int state);
/* No need for this to be exported to scripting */
void toggle_torsion_general_reverse();

void setup_residue_partial_alt_locs(short int state);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Backrub                                                 */
/*  ----------------------------------------------------------------------- */
/*! \name Backrubbing function */
/*! \{ */
/*! \brief Do a back-rub rotamer search (with autoaccept).

@return the success status, 0 for fail, 1 for successful fit.  */
int backrub_rotamer(int imol, const char *chain_id, int res_no,
		    const char *ins_code, const char *alt_conf);

/*! \brief apply rotamer backrub to the active atom of the intermediate atoms */
int backrub_rotamer_intermediate_atoms();

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

/*! \brief mask map by atom selection */
int mask_map_by_atom_selection(int map_mol_no, int coords_mol_no, const char *mmdb_atom_selection, short int invert_flag);

/*! \brief make chain masked maps

   needs to return a list of values
 */
int make_masked_maps_split_by_chain(int imol, int imol_map);

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
void set_check_waters_b_factor_limit(float f);
void set_check_waters_map_sigma_limit(float f);
void set_check_waters_min_dist_limit(float f);
void set_check_waters_max_dist_limit(float f);


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
/*! \brief return an improper list first is list of metals, second is
  list of waters that are coordinated with at least
  coordination_number of other atoms at distances less than or equal
  to dist_max.  Return scheme false on not able to make a list,
  otherwise a list of atoms and neighbours.  Can return scheme false
  if imol is not a valid molecule. */
SCM highly_coordinated_waters_scm(int imol, int coordination_number, float dist_max);

SCM metal_coordination_scm(int imol, float dist_max);
#endif
#ifdef USE_PYTHON
/*! \brief return a list first of waters, second metals that are
  coordinated with at least coordination_number of other atoms at
  distances less than or equal to dist_max. Return Python false on
  not able to make a list, otherwise a list of atoms and neighours.
  can return Python False if imol is not a valid molecule.  */
PyObject *highly_coordinated_waters_py(int imol, int coordination_number, float dist_max);
PyObject *metal_coordination_py(int imol, float dist_max);
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
/*! \brief apply the LSQ matches

@return an rtop pair (proper list) on good match, else false */
SCM apply_lsq_matches(int imol_reference, int imol_moving);
SCM get_lsq_matrix_scm(int imol_reference, int imol_moving);
#endif
#ifdef USE_PYTHON
/* Return an rtop pair (proper list) on good match, else False */
PyObject *apply_lsq_matches_py(int imol_reference, int imol_moving);
PyObject *get_lsq_matrix_py(int imol_reference, int imol_moving);
#endif /* PYTHON */
#endif /* __cplusplus */

/* poor old python programmers... */
int apply_lsq_matches_simple(int imol_reference, int imol_moving);

/* section Least-Squares plane interface */
void setup_lsq_deviation(int state);
void setup_lsq_plane_define(int state);
void unset_lsq_plane_dialog(); /* callback from destroy of widget */
void remove_last_lsq_plane_atom();

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

/*! \brief trim the molecule by the value in the B-factor column.

If an atom in a residue has a "B-factor" above (or below, if keep_higher is true) limit, then the whole residue is deleted */
void trim_molecule_by_b_factor(int imol, float limit, short int keep_higher);

/*! \brief convert the value in the B-factor column (typically pLDDT for AlphaFold models) to a temperature factor */
void pLDDT_to_b_factor(int imol);

/*! \} */


/*  ------------------------------------------------------------------------ */
/*                       povray/raster3d interface                           */
/*  ------------------------------------------------------------------------ */
/* make the text input to external programs */
/*! \name External Ray-Tracing */
/*! \{ */

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
/*! \brief set the flag to show waters as spheres for the Raster3D
representation. 1 show as spheres, 0 the usual stars. */
void set_raster3d_water_sphere(int istate);
/*! \brief set the font size (as a string) for raster3d*/
void set_raster3d_font_size(const char *size_in);
/*! \brief run raster3d and display the resulting image.  */
void raster_screen_shot(); /* run raster3d or povray and guile */
                           /* script to render and display image */
#ifdef USE_PYTHON
void raster_screen_shot_py(); /* run raster3d or povray and python */
#endif
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  citation notice                                         */
/*  ----------------------------------------------------------------------- */
void citation_notice_off();

/*  ----------------------------------------------------------------------- */
/*                  Superpose                                               */
/*  ----------------------------------------------------------------------- */
/* section Superposition (SSM) */
/*! \name Superposition (SSM) */
/*! \{ */

/*! \brief simple interface to superposition.
 *
 * @param imol1 the reference model index
 * @param imol2 the index of the superposed molecule

   Superpose all residues of imol2 onto imol1.  imol1 is reference, we
   can either move imol2 or copy it to generate a new molecule depending
   on the vaule of move_imol2_flag (1 for copy 0 for move). */
void superpose(int imol1, int imol2, short int move_imol2_flag);


/*! \brief chain-based interface to superposition.
 *
 * @param imol1 the reference model index
 * @param imol2 the index of the superposed molecule
 * @param chain_imol1 the chain_id of imol1
 * @param chain_imol2 the chain_id of imol2
 * @param chain_used_flag_imol1 should the chain-id be used for imol1 (1 for yes, 0 for no)
 * @param chain_used_flag_imol2 should the chain-id be used for imol2 (1 for yes, 0 for no)
 * @param move_imol2_coppy_flag flag to control if imol2 is
 *        copied (1) or moved (0)
 
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
   the vaule of move_imol2_flag (1 for move 0 for copy).

   @return the index of the superposed molecule - which could either be a
   new molecule (if move_imol2_flag was 1) or the imol2 or -1 (signifying
   failure to do the SSM superposition).
 *
 * @param imol1 the reference model index
 * @param imol2 the index of the superposed molecule
 * @param mmdb_atom_sel_str_1 the mmdb-format atom selection for imol1
 * @param mmdb_atom_sel_str_2 the mmdb-format atom selection for imol2
 * @param move_imol2_coppy_flag flag to control if imol2 is
 *        copied (1) or moved (0)
*/
int superpose_with_atom_selection(int imol1, int imol2,
				  const char *mmdb_atom_sel_str_1,
				  const char *mmdb_atom_sel_str_2,
				  short int move_imol2_copy_flag);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  NCS                                                     */
/*  ----------------------------------------------------------------------- */
/* section NCS */
/*! \name NCS */

/*! \{ */
/*! \brief set drawing state of NCS ghosts for molecule number imol   */
void set_draw_ncs_ghosts(int imol, int istate);
/*! \brief return the drawing state of NCS ghosts for molecule number
  imol.  Return -1 on imol is a bad molecule or no ghosts.  */
int draw_ncs_ghosts_state(int imol);

/*! \brief set bond thickness of NCS ghosts for molecule number imol   */
void set_ncs_ghost_bond_thickness(int imol, float f);
/*! \brief update ghosts for molecule number imol */
void ncs_update_ghosts(int imol); /* update ghosts */
/*! \brief make NCS map */
int make_dynamically_transformed_ncs_maps(int imol_model, int imol_map,
					  int overwrite_maps_of_same_name_flag);
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
int add_strict_ncs_from_mtrix_from_self_file(int imol);

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
/*! \brief Copy residue range to all related NCS chains.

  If the
  target residues do not exist in the peer chains, then create
  them. */
void copy_residue_range_from_ncs_master_to_others(int imol, const char *master_chain_id,
						  int residue_range_start, int residue_range_end);
#ifdef __cplusplus
#ifdef USE_GUILE

/*! \brief Copy chain from master to specified related NCS chains */
void copy_from_ncs_master_to_specific_other_chains_scm(int imol, const char *chain_id, SCM other_chain_id_list_scm);

/*! \brief Copy residue range to selected NCS chains

   If the target residues do not exist in the peer chains, then create
   them.
*/
/*! \brief return a list of NCS masters or scheme false */
SCM ncs_master_chains_scm(int imol);
void copy_residue_range_from_ncs_master_to_chains_scm(int imol, const char *master_chain_id,
						      int residue_range_start, int residue_range_end,
						      SCM chain_id_list);
/*! \brief Copy chain from master to a list of NCS chains */
void copy_from_ncs_master_to_chains_scm(int imol, const char *master_chain_id,
					SCM chain_id_list);
#endif
#ifdef USE_PYTHON
/*! \brief Copy chain from master to specified other NCS chains */
void copy_from_ncs_master_to_specific_other_chains_py(int imol, const char *chain_id, PyObject *other_chain_id_list_py);

PyObject *ncs_master_chains_py(int imol);
void copy_residue_range_from_ncs_master_to_chains_py(int imol, const char *master_chain_id,
						     int residue_range_start, int residue_range_end,
						     PyObject *chain_id_list);
void copy_from_ncs_master_to_chains_py(int imol, const char *master_chain_id,
				       PyObject *chain_id_list);
#endif
#endif

/*! \brief change the NCS master chain  (by number)*/
void ncs_control_change_ncs_master_to_chain(int imol, int ichain);
/*! \brief change the NCS master chain  (by chain_id)*/
void ncs_control_change_ncs_master_to_chain_id(int imol, const char *chain_id);
/*! \brief display the NCS master chain  */
void ncs_control_display_chain(int imol, int ichain, int state);

void set_ncs_matrix_type(int flag);
int get_ncs_matrix_state();

#ifdef __cplusplus
#ifdef USE_GUILE
/* Return the NCS differences as a list.

   e.g. ("B" "A" '(((1 "") (1 "") 0.4) ((2 "") (2 "") 0.3))
   i.e. ncs-related-chain its-master-chain-id and a list of residue
   info: (residue number matches: (this-resno this-inscode
   matching-mater-resno matching-master-inscode
   rms-atom-position-differences))) */
SCM ncs_chain_differences_scm(int imol, const char *master_chain_id);

/*! \brief Return the ncs chains id for the given molecule.

  return something like: '(("A" "B")) or '(("A" "C" "E") ("B"
  "D" "F")). The master chain goes in first.

   If imol does not have NCS ghosts, return scheme false.
*/
SCM ncs_chain_ids_scm(int imol);
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
/* Return the NCS differences as a list.

   e.g. ["B", "A", [[[1, ""], [1, ""], 0.4], [[2, ""], [2, ""], 0.3]]]
   i.e. ncs_related_chain its_master_chain_id and a list of residue
   info: [residue number matches: [this_resno, this_inscode,
   matching_master_resno, matching_master_inscode,
   rms_atom_position_differences]] */
PyObject *ncs_chain_differences_py(int imol, const char *master_chain_id);

/*! \brief Return the ncs chains id for the given molecule.

  return something like: [["A", "B"]] or [["A", "C", "E"], ["B",
  "D", "F"]]. The master chain goes in first.

   If imol does not have NCS ghosts, return python False.
*/
PyObject *ncs_chain_ids_py(int imol);
#endif  /* USE_PYTHON */

#ifdef USE_GUILE
/*! \brief get the NCS ghost description

@return false on bad imol or a list of ghosts on good imol.  Can
   include NCS rtops if they are available, else the rtops are False */
SCM ncs_ghosts_scm(int imol);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
/*! \brief Get the NCS ghosts description

@return False on bad imol or a list of ghosts on good imol.  Can
   include NCS rtops if they are available, else the rtops are False */
PyObject *ncs_ghosts_py(int imol);
#endif	/* USE_PYTHON */


#endif	/* __cplusplus */

/*! \} */

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


void set_place_helix_here_fudge_factor(float ff);


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

/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Autobuild nucleotides                                   */
/*  ----------------------------------------------------------------------- */

/*! \name Nucleotides*/
/*! \{ */

/*! \brief autobuild nucleic acid chains

   Find secondary structure local to current view in the current map.
   Add to a molecule called "NuclAcid", create it if needed.

   @return the index of the new molecule.*/

int find_nucleic_acids_local( float radius );

/*! \} */


/*  ----------------------------------------------------------------------- */
/*             New Molecule by Various Selection                            */
/*  ----------------------------------------------------------------------- */
/*! \name New Molecule by Section Interface */
/*! \{ */
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


#ifdef __cplusplus
#ifdef USE_PYTHON
/*! \brief create a new molecule that consists of only the atoms
  of the specified list of residues
@return the new molecule number, -1 means an error. */
int new_molecule_by_residue_specs_py(int imol, PyObject *residue_spec_list_py);
#endif /* USE_PYTHON */

#ifdef USE_GUILE
/*! \brief create a new molecule that consists of only the atoms
  of the specified list of residues
@return the new molecule number, -1 means an error. */
int new_molecule_by_residue_specs_scm(int imol, SCM residue_spec_list_scm);
#endif /* USE_GUILE */
#endif /* __cplusplus */

/*! \} */


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

/*! \{ */
/*!  \brief create a molecule of idea nucleotides

use the given sequence (single letter code)

RNA_or_DNA is either "RNA" or "DNA"

form is either "A" or "B"

@return the new molecule number or -1 if a problem */
int ideal_nucleic_acid(const char *RNA_or_DNA, const char *form,
		       short int single_stranged_flag,
		       const char *sequence);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE

/*! \brief get the pucker info for the specified residue

 @return scheme false if residue not found,  otherwise
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
/*! \brief add base pairs for the given residue range, modify molecule imol by creating a new chain */
int watson_crick_pair_for_residue_range(int imol, const char * chain_id, int resno_start, int resno_end);

/* not for user level */
void setup_base_pairing(int state);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  sequence file (assignment)                              */
/*  ----------------------------------------------------------------------- */
/* section Sequence File (Assignment/Association) */
/*! \name Sequence File (Assignment/Association) */
/*! \{ */

/*! \brief Print the sequence to the console of the given molecule */
void print_sequence_chain(int imol, const char *chain_id);

/*! \brief optionally write the sequence to the file for the given molecule,
    optionally in PIR format */
void print_sequence_chain_general(int imol, const char *chain_id,
                                  short int pir_format,
                                  short int file_output,
                                  const char *file_name);

/*! \brief Assign a FASTA sequence to a given chain in the  molecule */
void assign_fasta_sequence(int imol, const char *chain_id_in, const char *seq);
/*! \brief Assign a PIR sequence to a given chain in the molecule.  If
  the chain of the molecule already had a chain assigned to it, then
  this will overwrite that old assignment with the new one. */
void assign_pir_sequence(int imol, const char *chain_id_in, const char *seq);
/* I don't know what this does. */
void assign_sequence(int imol_model, int imol_map, const char *chain_id);
/*! \brief Assign a sequence to a given molecule from (whatever) sequence
  file by alignment. */
void assign_sequence_from_file(int imol, const char *file);
/*! \brief Assign a sequence to a given molecule from a simple string */
void assign_sequence_from_string(int imol, const char *chain_id_in, const char *seq);
/*! \brief Delete all the sequences from a given molecule */
void delete_all_sequences_from_molecule(int imol);
/*! \brief Delete the sequence for a given chain_id from a given molecule */
void delete_sequence_by_chain_id(int imol, const char *chain_id_in);

/*! \brief Associate the sequence to the molecule - to be used later for sequence assignment (.c.f assign_pir_sequence)   */
void associate_sequence_from_file(int imol, const char *file_name);

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
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Surfaces                                                */
/*  ----------------------------------------------------------------------- */
/* section Surface Interface */
/*! \name Surface Interface */
/*! \{ */
/*! \brief draw surface of molecule number imol

if state = 1 draw the surface (normal representation goes away)

if state = 0 don't draw surface */
void do_surface(int imol, int istate);
int molecule_is_drawn_as_surface_int(int imol); /* predicate */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief draw the surface of the imolth molecule clipped to the
  residues given by residue_specs.

  residue_specs must not contain spec for waters (you wouldn't want to
  surface over waters anyway).
 */
void do_clipped_surface_scm(int imol, SCM residue_specs);
#endif /*  USE_GUILE */
#ifdef USE_PYTHON
/*! \brief draw the surface of the imolth molecule clipped to the
  residues given by residue_specs */
void do_clipped_surface_py(int imol, PyObject *residue_specs);
#endif /*  USE_PYTHON */
#endif	/* __cplusplus */

/*! make molecular surface for given atom selection

    per-chain functions can be added later */
void make_molecular_surface(int imol, const char *selection_string);

/*! make electrostatics surface for given atom selection

  per-chain functions can be added later */
void make_electrostatic_surface(int imol, const char *selection_string);

void set_electrostatic_surface_charge_range(float v);
float get_electrostatic_surface_charge_range();

/*! \brief simple on/off screendoor transparency at the moment, an
  opacity > 0.0 will turn on screendoor transparency (stippling). */
void set_transparent_electrostatic_surface(int imol, float opacity);

/*! \brief return 1.0 for non transparent and 0.5 if screendoor
  transparency has been turned on. */
float get_electrostatic_surface_opacity(int imol);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  FFfearing                                               */
/*  ----------------------------------------------------------------------- */
/*! \name FFFearing */
/*! \{ */
/*! \brief fffear search model in molecule number imol_model in map
   number imol_map */
int fffear_search(int imol_model, int imol_map);
/*! \brief set and return the fffear angular resolution in degrees */
void set_fffear_angular_resolution(float f);
/*! \brief return the fffear angular resolution in degrees */
float fffear_angular_resolution();
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  remote control                                          */
/*  ----------------------------------------------------------------------- */
/* section Remote Control */
/*! \name Remote Control */
/*! \{ */
/*! \brief try to make socket listener */
void make_socket_listener_maybe();
void set_coot_listener_socket_state_internal(int sock_state);

/*! \brief feed the main thread a scheme script to evaluate */
void set_socket_string_waiting(const char *s);
/*! \brief feed the main thread a python script to evaluate */
void set_socket_python_string_waiting(const char *s);

void set_remote_control_port(int port_number);
int get_remote_control_port_number();


/* tooltip */
void set_tip_of_the_day_flag(int state);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Display lists                                           */
/*  ----------------------------------------------------------------------- */
/* section Display Lists for Maps */
/*! \name Display Lists for Maps */
/*! \brief Should display lists be used for maps? It may speed things
  up if these are turned on (or off) - depends on graphics card and
  drivers.  Pass 1 for on, 0 for off. */
void set_display_lists_for_maps(int i);

/*! \brief return the state of display_lists_for_maps.   */
int display_lists_for_maps_state();

/* update the maps to the current position - rarely needed */
void update_maps();
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Browser Help                                            */
/*  ----------------------------------------------------------------------- */
/*! \name Browser Interface */
/*! \{ */
/*! \brief try to open given url in Web browser */
void browser_url(const char *url);
/*! \brief set command to open the web browser,

examples are "open" or "mozilla" */
void set_browser_interface(const char *browser);

/*! \brief the search interface

find words, construct a url and open it. */
void handle_online_coot_search_request(const char *entry_text);
/*! \} */

// #include "c-interface-generic-objects.h"


/*  ----------------------------------------------------------------------- */
/*                  Molprobity interface                                    */
/*  ----------------------------------------------------------------------- */
/*! \name Molprobity Interface */
/*! \{ */
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

/* state is 1 for on and 0 for off */
void set_do_coot_probe_dots_during_refine(short int state);

/* get state: 1 for on and 0 for off */
short int get_do_coot_probe_dots_during_refine();


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
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the parsed user mod fields from the PDB file
  file_name (output by reduce most likely) */
PyObject *user_mods_py(const char *file_name);
#endif /* USE_PYTHON */
#endif	/* c++ */

/*! \} */


/*  ----------------------------------------------------------------------- */
/*           Sharpen                                                        */
/*  ----------------------------------------------------------------------- */
/*! \name Map Sharpening Interface */
/*! \{ */
/*! \brief Sharpen map imol by b_factor (note (of course) that positive numbers
    blur the map).  */
void sharpen(int imol, float b_factor);
void sharpen_with_gompertz_scaling(int imol, float b_factor, short int try_gompertz, float gompertz_factor);

/*! \brief set the limit of the b-factor map sharpening slider (default 30) */
void set_map_sharpening_scale_limit(float f);
/*! \} */
/* ---------------------------------------------------------------------------- */
/*	Density Map Kurtosis							*/
/* ----------------------------------------------------------------------------	*/
float optimal_B_kurtosis(int imol);


/*  ----------------------------------------------------------------------- */
/*           Intermediate Atom Manipulation                                 */
/*  ----------------------------------------------------------------------- */

/*! \name Intermediate Atom Manipulation Interface */
/*! \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM drag_intermediate_atom_scm(SCM atom_spec, SCM position);
#endif
#ifdef USE_PYTHON
PyObject *drag_intermediate_atom_py(PyObject *atom_spec, PyObject *position);

//! \brief add a target position for an intermediate atom and refine
//
// A function requested by Hamish.
// This aplies to intermediate atoms (add_extra_target_position_restraint)
// does not. This activates refinement after the restraint is added (add_extra_target_position_restraint
// does not).
PyObject *add_target_position_restraint_for_intermediate_atom_py(PyObject *atom_spec, PyObject *position);

// and the multiple-atom version of that (so that they can be applied at the same time)
PyObject *add_target_position_restraints_for_intermediate_atoms_py(PyObject *atom_spec_position_list);
#endif
#endif /* c++ */
/*! \} */

/*  ----------------------------------------------------------------------- */
/*           Fixed Atom Manipulation                                        */
/*  ----------------------------------------------------------------------- */

/*! \name Marking Fixed Atom Interface */
/*! \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM mark_atom_as_fixed_scm(int imol, SCM atom_spec, int state);
int mark_multiple_atoms_as_fixed_scm(int imol, SCM atom_spec_list, int state);
#endif
#ifdef USE_PYTHON
PyObject *mark_atom_as_fixed_py(int imol, PyObject *atom_spec, int state);
int mark_multiple_atoms_as_fixed_py(int imol, PyObject *atom_spec_list, int state);
#endif
#endif /* c++ */

void setup_fixed_atom_pick(short int ipick, short int is_unpick);

/*! \brief clear all fixed atoms */
void clear_all_fixed_atoms(int imol);
void clear_fixed_atoms_all();

/* produce debugging output from problematic atom picking  */
void set_debug_atom_picking(int istate);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Partial Charge                                          */
/*  ----------------------------------------------------------------------- */
/*! \name Partial Charges */
/*! \{ */
/*! \brief show the partial charges for the residue of the given specs
   (charges are read from the dictionary) */
void show_partial_charge_info(int imol, const char *chain_id, int resno, const char *ins_code);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  EM Interface                                            */
/*  ----------------------------------------------------------------------- */
/*! \name EM interface */
/*! \{ */
/*! \brief Scale the cell, for use with EM maps, where the cell needs
   to be adjusted.  Use like:  (scale-cell 2 1.012 1.012 1.012). Return error
   status, 1 means it worked, 0 means it did not work. */
int scale_cell(int imol_map, float fac_u, float fac_v, float fac_w);

/* create a number of maps by segmenting the given map, above the
   (absolute) low_level.  New maps are on the same grid as the input
   map.  */
void segment_map(int imol_map, float low_level);

void segment_map_multi_scale(int imol_map, float low_level, float b_factor_inc, int n_rounds);

/*! \brief make a map histogram */
void map_histogram(int imol_map);

/*! \brief ignore pseudo-zeros when calculationg maps stats (default 1 = true) */
void set_ignore_pseudo_zeros_for_map_stats(short int state);

/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  CCP4i Interface                                         */
/*  ----------------------------------------------------------------------- */
/*! \name CCP4mg Interface */
/*! \{ */
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

/*! \brief allow the user to not add ccp4i directories to the file choosers

use state=0 to turn it off */
void set_add_ccp4i_projects_to_file_dialogs(short int state);

/*! \brief write a ccp4mg picture description file */
void write_ccp4mg_picture_description(const char *filename);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Dipoles                                                 */
/*  ----------------------------------------------------------------------- */
/*! \name Dipoles */
/*! \{ */
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
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Patterson                                               */
/*  ----------------------------------------------------------------------- */
/*! \brief Make a patterson molecule

\return a new molecule number or -1 on failure */
int make_and_draw_patterson(const char *mtz_file_name,
			    const char *f_col,
			    const char *sigf_col);
/*! \brief Make a patterson molecule

\return a new molecule number or -1 on failure */
int make_and_draw_patterson_using_intensities(const char *mtz_file_name,
					      const char *i_col,
					      const char *sigi_col);

/*  ----------------------------------------------------------------------- */
/*                  Laplacian                                               */
/*  ----------------------------------------------------------------------- */
/*! \name Aux functions */
/*! \{ */
/*! \brief Create the "Laplacian" (-ve second derivative) of the given map. */
int laplacian (int imol);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  PKGDATADIR                                              */
/*  ----------------------------------------------------------------------- */
/*! \name PKGDATADIR */
/*! \{ */
#ifdef __cplusplus
#ifdef USE_PYTHON
PyObject *get_pkgdatadir_py();
#endif /* USE_PYTHON */
#ifdef USE_GUILE
// note: built-ins: (%package-data-dir) and %guile-build-info
SCM get_pkgdatadir_scm();
#endif
#endif /*  __cplusplus */
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  SMILES                                                  */
/*  ----------------------------------------------------------------------- */
/*! \name SMILES */
/*! \{ */
/*! \brief display the SMILES string dialog */
void do_smiles_gui();
/*! \} */
/*  ----------------------------------------------------------------------- */
/*                  Fun                                                     */
/*  ----------------------------------------------------------------------- */
/* section Fun */
void do_tw();

/*  ----------------------------------------------------------------------- */
/*                  Phenix Support                                          */
/*  ----------------------------------------------------------------------- */
/*! \name PHENIX Support */
/*! \{ */
/*! \brief set the button label of the external Refinement program */
void set_button_label_for_external_refinement(const char *button_label);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Text                                                    */
/*  ----------------------------------------------------------------------- */
/*! \name Graphics Text */
/*! \{ */
/*! \brief Put text at x,y,z

@return a text handle

size variable is currently ignored.*/
int place_text(const char*text, float x, float y, float z, int size);

/*! \brief Remove "3d" text item */
void remove_text(int text_handle);

void edit_text(int text_handle, const char *new_text);

/*! \brief return the closest text that is with r A of the given
  position.  If no text item is close, then return -1 */
int text_index_near_position(float x, float y, float z, float r);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  PISA Interface                                      */
/*  ----------------------------------------------------------------------- */
/*! \name PISA Interaction */
/*! \{ */
/*! \brief return the molecule number of the interacting
  residues. Return -1 if no new model was created. Old, not very useful. */
int pisa_interaction(int imol_1, int imol_2);
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief the scripting interface, called from parsing the PISA XML
   interface description

   An interface_description_scm is a record detailing the interface.
   A record contains the bsa, asa, and 2 molecule records.  Molecule
   records contain list of residue records.  The interface (dots) is
   be made from these lists of residue records. Note of course that
   imol_2 (or 1) can be a symmetry copy of (part of) mol_1 (or 2).

   Return the dot indexes (currently -1)

*/
SCM handle_pisa_interfaces_scm(SCM interfaces_description_scm);

/* internal function */
SCM pisa_molecule_record_residues(SCM molecule_record_1);
SCM pisa_molecule_record_chain_id(SCM molecule_record_1);
void add_pisa_interface_bond_scm(int imol_1, int imol_2, SCM pisa_bond_scm,
				 int interface_number);


/* clear out and undisplay all pisa interface descriptions. */
void pisa_clear_interfaces();
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief the scripting interface, called from parsing the PISA XML
   interface description

   An interface_description_scm is a record detailing the interface.
   A record contains the bsa, asa, and 2 molecule records.  Molecule
   records contain list of residue records.  The interface (dots) is
   be made from these lists of residue records. Note of course that
   imol_2 (or 1) can be a symmetry copy of (part of) mol_1 (or 2).

   Return the dot indexes (currently -1)

*/
PyObject *handle_pisa_interfaces_py(PyObject *interfaces_description_py);

/* internal function */
/* PyObject *pisa_molecule_record_residues_py(PyObject *molecule_record_1); */
/* PyObject *pisa_molecule_record_chain_id_py(PyObject *molecule_record_1); */
void add_pisa_interface_bond_py(int imol_1, int imol_2, PyObject *pisa_bond_py,
                                 int interface_number);

/* clear out and undisplay all pisa interface descriptions. */
void pisa_clear_interfaces();
#endif /* USE_PYTHON */
#endif /* c++ */


/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Jiggle fit                                              */
/*  ----------------------------------------------------------------------- */
/*! \name Jiggle Fit */
/*! \{ */
/*!  \brief jiggle fit to the current refinment map.  return < -100 if
  not possible, else return the new best fit for this residue.  */
float fit_to_map_by_random_jiggle(int imol, const char *chain_id, int resno, const char *ins_code,
                                  int n_trials, float jiggle_scale_factor);

/*!  \brief jiggle fit the molecule to the current refinment map.  return < -100 if
  not possible, else return the new best fit for this molecule.  */
float fit_molecule_to_map_by_random_jiggle(int imol, int n_trials, float jiggle_scale_factor);
/*!  \brief jiggle fit the molecule to the current refinment map.  return < -100 if
  not possible, else return the new best fit for this molecule - create a map that is blurred
  by the given factor for fitting  */
float fit_molecule_to_map_by_random_jiggle_and_blur(int imol, int n_trials, float jiggle_scale_factor, float map_blur_factor);

/*!  \brief jiggle fit the chain to the current refinment map.  return < -100 if
  not possible, else return the new best fit for this chain.  */
float fit_chain_to_map_by_random_jiggle(int imol, const char *chain_id, int n_trials, float jiggle_scale_factor);

/*!  \brief jiggle fit the chain to the current refinment map
 *
 * Use a map that is blurred by the give factor for fitting.
 * @return < -100 if not possible, else return the new best fit for this chain.  */
float fit_chain_to_map_by_random_jiggle_and_blur(int imol, const char *chain_id, int n_trials, float jiggle_scale_factor, float map_blur_factor);

/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  SBase interface                                         */
/*  ----------------------------------------------------------------------- */
/*! \name SBase interface */
/*! \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return a list of compoundIDs of in SBase of which the
  given string is a substring of the compound name */
SCM matching_compound_names_from_sbase_scm(const char *compound_name_fragment);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return a list of compoundIDs of in SBase of which the
  given string is a substring of the compound name */
PyObject *matching_compound_names_from_sbase_py(const char *compound_name_fragment);
#endif /* USE_PYTHON */
#endif /*__cplusplus */

/*! \brief return the new molecule number of the monomer.

The monomer will have chainid "A" and residue number 1.

Return -1 on failure to get monomer. */
int get_ccp4srs_monomer_and_dictionary(const char *comp_id);

/*! \brief same as above but using old name for back-compatibility */
int get_sbase_monomer(const char *comp_id);

/*! \} */


/* Needs a/the correct section */
/* add a linked residue based purely on dictionary template.
   For addition of NAG to ASNs typically.

   This doesn't work with residues with alt confs.

   Link type is the refmac dictionary link type (e.g. "ASN-NAG").

   return success status (0 = fail).
*/
int add_linked_residue(int imol, const char *chain_id, int resno, const char *ins_code,
		       const char *new_residue_comp_id, const char *link_type, int n_trials);
#ifdef __cplusplus
#ifdef USE_GUILE
// mode is either 1: add  2: add and fit  3: add, fit and refine
SCM add_linked_residue_scm(int imol, const char *chain_id, int resno, const char *ins_code,
			   const char *new_residue_comp_id, const char *link_type, int mode);
#endif
#ifdef USE_PYTHON
PyObject *add_linked_residue_py(int imol, const char *chain_id, int resno, const char *ins_code,
				const char *new_residue_comp_id, const char *link_type, int mode);
#endif
#endif
void set_add_linked_residue_do_fit_and_refine(int state);

/*  ----------------------------------------------------------------------- */
/*               Flattened Ligand Environment View  Interface               */
/*  ----------------------------------------------------------------------- */
/*! \name FLE-View */
/*! \{ */

void fle_view(int imol, const char *chain_id, int res_no, const char *ins_code, float dist_max);

/* delete these other functions  */

void fle_view_internal(int imol, const char *chain_id, int res_no,
		       const char *ins_code,
		       int imol_ligand_fragment,
		       const char *prodrg_output_flat_mol_file_name,
		       const char *prodrg_output_flat_pdb_file_name,
		       const char *prodrg_output_3d_pdb_file_name,
		       const char *prodrg_output_dict_cif_file_name);
/* for command-line operation */
void fle_view_internal_to_png(int imol, const char *chain_id, int res_no,
			      const char *ins_code,
			      int imol_ligand_fragment,
			      const char *prodrg_output_flat_mol_file_name,
			      const char *prodrg_output_flat_pdb_file_name,
			      const char *prodrg_output_3d_pdb_file_name,
			      const char *prodrg_output_dict_cif_file_name,
			      int output_to_png_file_flag,
			      const char *png_file_name);

void fle_view_with_rdkit(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius);
void fle_view_with_rdkit_to_png(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *png_file_name);
void fle_view_with_rdkit_to_svg(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *svg_file_name);

void fle_view_with_rdkit_internal(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *file_format, const char *file_name);

/*! \brief set the maximum considered distance to water

default 3.25 A.  */
void fle_view_set_water_dist_max(float dist_max);
/*! \brief set the maximum considered hydrogen bond distance

default 3.9 A.  */
void fle_view_set_h_bond_dist_max(float h_bond_dist_max);

/*! \brief Add hydrogens to specificied residue

@return success status.

use RDKit for enterprise version
 */
int sprout_hydrogens(int imol, const char *chain_id, int res_no, const char *ins_code);

/*! \} */


/*  ----------------------------------------------------------------------- */
/*               LSQ-improve                                                */
/*  ----------------------------------------------------------------------- */
/*! \name LSQ-improve */
/*! \{ */
/*! \brief an slightly-modified implementation of the "lsq_improve"
  algorithm of Kleywegt and Jones (1997).

  Note that if a residue selection is specified in the residue
  selection(s), then the first residue of the given range must exist
  in the molecule (if not, then mmdb will not select any atoms from
  that molecule).

  Kleywegt and Jones set n_res to 4 and dist_crit to 6.0.

 */
void lsq_improve(int imol_ref, const char *ref_selection,
		 int imol_moving, const char *moving_selection,
		 int n_res, float dist_crit);
/*! \} */



/*  ----------------------------------------------------------------------- */
/* Multirefine interface (because in guile-gtk there is no way to
   insert toolbuttons into the toolbar) so this
   rather kludgy interface.  It should go when we
   move to guile-gnome, I think.                                            */
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

/*! \brief load tutorial model and data
 *
 * Loads an example dataset - the sample is an RNase structure (model and maps) and is
 * used for learning and testing.
 *
 * This is the standard Coot tutorial dataset for practicing model building
 * and validation.
 *
 * */
void load_tutorial_model_and_data();


/*  ----------------------------------------------------------------------- */
/*                         single-model view                                */
/*  ----------------------------------------------------------------------- */
/*! \name single-model view */
/*! \{ */
/*! \brief put molecule number imol to display only model number imodel */
void single_model_view_model_number(int imol, int imodel);
/*! \brief the current model number being displayed

return 0 on non-multimodel-molecule. */
int single_model_view_this_model_number(int imol);
/*! \brief change the representation to the next model number to be displayed

return 0 on non-multimodel-molecule.
*/
int single_model_view_next_model_number(int imol);
/*! \brief change the representation to the previous model number to be displayed

return 0 on non-multimodel-molecule. */
int single_model_view_prev_model_number(int imol);
/*! \} */



/*  ----------------------------------------------------------------------- */
/*                  update self                                             */
/*  ----------------------------------------------------------------------- */
/* this function is here because it is called by c_inner_main() (ie. need a c interface). */
void run_update_self_maybe(); /* called when --update-self given at command line */

/*  ----------------------------------------------------------------------- */
/*                    keyboarding mode                                      */
/*  ----------------------------------------------------------------------- */
void show_go_to_residue_keyboarding_mode_window();
void handle_go_to_residue_keyboarding_mode(const char *text); /* should this be here? */

/*  ----------------------------------------------------------------------- */
/*                    graphics ligand view                                  */
/*  ----------------------------------------------------------------------- */
/*! \name graphics 2D ligand view */
/*! \{ */
/*! \brief set the graphics ligand view state

 (default is 1 (on)). */
void set_show_graphics_ligand_view(int state);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  experimental                                            */
/*  ----------------------------------------------------------------------- */
/*! \name Experimental */
/*! \{ */

void fetch_and_superpose_alphafold_models_using_active_molecule();

// void add_ligand_builder_menu_item_maybe(); // what does this do?

/*!  \brief display the ligand builder dialog */
void start_ligand_builder_gui();

/*  ----------------------------------------------------------------------- */
/*                  end                                                     */
/*  ----------------------------------------------------------------------- */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM all_molecule_rotamer_score(int imol);
SCM all_molecule_ramachandran_score(int imol); /* a stub currently */
#endif /* USE_GUILE */

#ifdef USE_PYTHON
/**
 * @brief Compute rotamer score for an entire molecule and return result as a Python object.
 *
 * This wrapper computes the rotamer score information for molecule number @p imol
 * and returns a Python list containing the numeric score and the number of
 * rotamer-bearing residues.
 *
 * Parameters
 * ----------
 * @param imol
 *     Model (molecule) index to analyze.
 *
 * Return value
 * ------------
 * Returns a NEW reference to a Python object. Two possible outcomes:
 *
 * - Success: a Python list of length 2 (PyList), with elements:
 *     0 : float — overall rotamer score (PyFloat)
 *     1 : int   — number of rotamer residues considered (PyLong)
 *
 * - Failure / invalid model index: Py_False (Python False).
 *   The implementation INCREFs Py_False before returning, so the caller receives
 *   a new reference in the failure case as well.
 *
 * Reference counting
 * -----------------
 * The returned PyObject* is a new reference. The caller is responsible for
 * DECREFing it when finished.
 *
 * Notes
 * -----
 * - Callers should detect the failure case by testing with PyBool_Check (Py_False).
 * - Ensure the Python GIL is held when calling this function from non-Python threads.
 */
PyObject *all_molecule_rotamer_score_py(int imol);
#endif /* USE_PYTHON */

#ifdef USE_PYTHON
/**
 * @brief Compute overall and per-residue Ramachandran statistics for a molecule and
 *        return the results as a Python object.
 *
 * This wrapper gathers the Ramachandran score information computed for molecule
 * number @p imol and returns a Python list with six elements describing the
 * overall scores and per-residue details.
 *
 * Parameters
 * ----------
 * @param imol
 *     Model (molecule) index to analyze.
 *
 * Return value
 * ------------
 * @return a NEW reference to a Python object. Two possible outcomes:
 *
 * - Success: a Python list of length 6 (PyList), with elements:
 *     0 : float   — overall Ramachandran score (PyFloat)
 *     1 : int     — number of residues considered (PyLong)
 *     2 : float   — Ramachandran score restricted to non-secondary-structure residues (PyFloat)
 *     3 : int     — number of residues used for the non-secondary-structure score (PyLong)
 *     4 : int     — number of zero-score residues (PyLong)
 *     5 : list    — info_by_residue: a list with one entry per residue (length == number of residues).
 *                     Each entry is either:
 *                       - a list of four items:
 *                           [ phi_psi_list, residue_spec_py, residue_score, res_names_list ]
 *                             * phi_psi_list: list of two floats [phi, psi]
 *                             * residue_spec_py: Python representation of the residue spec (see residue_spec_to_py)
 *                             * residue_score: float (PyFloat)
 *                             * res_names_list: list of three strings [prev_res_name, this_res_name, next_res_name]
 *                       - the integer -1 as a placeholder if per-residue info could not be computed for that index.
 *
 * - Failure / invalid model index: Py_False (Python False). The implementation INCREFs Py_False before returning,
 *   so the caller receives a new reference in this case as well.
 *
 * Reference counting
 * -----------------
 * The returned PyObject* is a new reference. The caller is responsible for DECREFing it when finished.
 *
 * Notes
 * -----
 * - Callers should check the return with PyBool_Check to detect the failure case (Py_False).
 * - Some per-residue entries may be -1 (an integer) if residue data was unavailable.
 * - The function must be called with appropriate GIL handling if invoked from non-Python threads.
 */

PyObject *all_molecule_ramachandran_score_py(int imol); /* a stub currently */
#endif /* USE_PYTHON */

#ifdef USE_PYTHON
/**
 * @brief Return the Ramachandran region annotation for a molecule as a Python list.
 *
 * This wrapper returns per-residue region information computed for molecule @p imol.
 * It queries the internal Ramachandran scoring machinery and returns a Python list
 * of (residue_spec, region_int) pairs for residues that lie in the computed region.
 *
 * Parameters
 * ----------
 * @param imol
 *     Model (molecule) index to query.
 *
 * Return value
 * ------------
 * Returns a NEW reference to a Python object. There are two possible outcomes:
 *
 * - Success: a Python list (PyList) of length N > 0, where each element is a 2-tuple:
 *     ( residue_spec_py, region_code )
 *     * residue_spec_py: Python representation of the residue spec (as produced by residue_spec_to_py).
 *     * region_code: integer (PyLong) — the integer label associated with that residue's Ramachandran region
 *       (as provided by the underlying rama score/region computation).
 *
 * - Failure / no region entries / invalid model index: Py_False (Python False).
 *   The implementation INCREFs Py_False before returning, so the caller receives a new reference
 *   in this case as well.
 *
 * Reference counting
 * -----------------
 * The returned PyObject* is a new reference. The caller is responsible for DECREFing it when finished.
 *
 * Notes
 * -----
 * - Callers should detect the failure/empty case by testing with PyBool_Check (Py_False).
 * - The exact meaning of the integer region_code is defined by the internal Ramachandran scoring code;
 *   consult the implementation or documentation for interpretation of region codes.
 * - Ensure the Python GIL is held when calling this function from non-Python threads.
 */
PyObject *all_molecule_ramachandran_region_py(int imol);
#endif /* USE_PYTHON */
#endif /* __cplusplus */

/*! \brief globularize the molecule.

This is not guaranteed to generate the correct biological entity, but will bring together
molecules (chains/domains) that are dispersed throughout the unit cell.

@param imol the molecule index.
*/
void globularize(int imol);

#ifdef __cplusplus
#ifdef USE_GUILE
/*!

    20100616 This doesn't get into the doxygen documentation for some
    reason I can't figure out.

    \fn user_defined_click_py(int n_clicks, SCM func);

    \brief run a user defined function

      Define a function func which runs after the user has made
      n_clicked atom picks.  func is called with a list of atom
      specifiers - with leading molecule number.
*/
void user_defined_click_scm(int n_clicks, SCM func);
#endif
#ifdef USE_PYTHON
void user_defined_click_py(int n_clicks, PyObject *func);
#endif /* PYTHON */
#endif /* __cplusplus */

/*! \} */

#endif /* C_INTERFACE_H */
END_C_DECLS
