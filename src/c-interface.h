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
#include "pytypedefs.h"
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



#define COOT_SCHEME_DIR "COOT_SCHEME_DIR"
#define COOT_PYTHON_DIR "COOT_PYTHON_DIR"

/*  ------------------------------------------------------------------------ */
/*                         Startup Functions:                                */
/*  ------------------------------------------------------------------------ */
#ifdef USE_GUILE
/*! \brief load the .scm files in the directories listed in the
  environment variable \c COOT_SCHEME_EXTRAS_DIR

  The variable is a list of directories separated by ":" (";" on
  Windows). Nothing happens if the variable is not set. */
void try_load_scheme_extras_dir();
#endif /* USE_GUILE */
#ifdef USE_PYTHON

/*! \brief load the .py files in the directories listed in the
  environment variable \c COOT_PYTHON_EXTRAS_DIR

  The variable is a list of directories separated by ":" (";" on
  Windows). Nothing happens if the variable is not set. */
void try_load_python_extras_dir();
#endif /* USE_PYTHON */

/* section Startup Functions */
/*!  \name Startup Functions */
/*! \{ */
/*!  \brief tell coot that you prefer to run python scripts if/when
  there is an option to do so. */
void set_prefer_python();

/*! \brief the python-prefered mode.

This is available so that the scripting functions know whether or not
to put themselves into menus as menu items.

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

    Some people don't like to see the full path names in the display
    manager. By default (0) only the file name is shown for model
    molecules; an argument of 1 turns on display of the full path.

    @param i 1 to show full paths, 0 to show just the file name
*/
void set_show_paths_in_display_manager(int i);

/*! \brief return the state of showing full paths in the Display Manager

   @return 1 for "yes, display paths", 0 for not
 */
int show_paths_in_display_manager_state();

/*! \brief add an extension to be treated as coordinate files

   The extension should include the leading dot, e.g. ".pdb" (the
   default list includes ".pdb", ".pdb.gz", ".ent", ".cif" and ".mmcif").

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

/*! \brief add an extension to be treated as map files
   @param ext the extension to be added (including the leading dot)
*/
void add_map_glob_extension(const char *ext);

/*! \brief remove an extension to be treated as coordinate files
   @param ext the extension to be removed
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

/*! \brief remove an extension to be treated as map files
   @param ext the extension to be removed
*/
void remove_map_glob_extension(const char *ext);

/*! \brief sort files in the file selection by date?

  some people like to have their files sorted by date by default */
void set_sticky_sort_by_date();

/*! \brief do not sort files in the file selection by date?

  removes the sorting of files by date */
void unset_sticky_sort_by_date();

/*! \brief on opening a file selection dialog, pre-filter the files.

set to 1 to pre-filter, 0 for no pre-filtering. The default is 1 (on).

@param istate 1 for on, 0 for off */
void set_filter_fileselection_filenames(int istate);

/*! \brief return the state of pre-filtering in the file selection dialog

@return 1 for pre-filtering on, 0 for off */
int filter_fileselection_filenames_state();

/*! \brief is the given file name suitable to be read as coordinates?

The test is on the file name extension only, which is compared
against the list of coordinates extensions (see
add_coordinates_glob_extension()).

@return 1 if the extension is a coordinates extension, 0 if not */
short int file_type_coords(const char *file_name);

/*! \brief display the open coordinates dialog */
void open_coords_dialog();

/* this flag set chooser as default for windows, otherwise use
  selector 0 is selector 1 is chooser */


/* --- CHECKME - do these need to be here? --- */
/*! \brief set the file chooser style

@param istate 1 (the default) to use the file chooser, 0 for the old-style file selector */
void set_file_chooser_selector(int istate);
/*! \brief return the file chooser style

@return 1 for the file chooser, 0 for the old-style file selector */
int file_chooser_selector_state();
/*! \brief set the file chooser overwrite flag

@param istate 1 (the default) for overwrite-protect, 0 for overwrite */
void set_file_chooser_overwrite(int istate);
/*! \brief return the file chooser overwrite flag

@return 1 for overwrite-protect, 0 for overwrite */
int file_chooser_overwrite_state();

/*! \brief show the export map GUI

@param export_map_fragment if non-zero, the dialog also shows the
radius entry, so that a map fragment (rather than the whole map) is
exported */
void export_map_gui(short int export_map_fragment);

/*! \} */


/*! \name Widget Utilities */
/*! \{ */

/*! \brief set the main window title.

function added for Lothar Esser. Does nothing if there is no
graphics interface or if \c s is empty.

@param s the new title */
void set_main_window_title(const char *s);

/*! \brief set the state of the validation graphs box
 *
 * By "docked" I mean, in the main window. The alternative
 * is a floating dialog.
 *
 * @param state 0 is not docked, 1 is docked
 */
void set_validation_graphs_is_docked(short int state);

/*! \} */

/*  -------------------------------------------------------------------- */
/*                   mtz and data handling utilities                     */
/*  -------------------------------------------------------------------- */
/* section MTZ and data handling utilities */
/*! \name  MTZ and data handling utilities */
/*! \{ */
/* We try as .phs and .cif files first */

/*! \brief given a filename, try to read it as a data file

   We try as .phs and .cif files first. If the file needs column
   label selection (e.g. an MTZ file), the column label selector
   dialog is displayed. Does nothing if there is no graphics
   interface.

   @param filename the reflection data file name */
void manage_column_selector(const char *filename);

/*! \} */

/*  -------------------------------------------------------------------- */
/*                     Molecule Functions       :                        */
/*  -------------------------------------------------------------------- */
/* section Molecule Info Functions */
/*! \name Molecule Info Functions */
/*! \{ */

/*! \brief the number of residues in chain chain_id and molecule number imol

  Only the first model is considered.

  @param chain_id the chain id
  @param imol the model molecule index
  @return the number of residues, or -1 if imol is not a valid model
  molecule or the chain was not found
*/
int chain_n_residues(const char *chain_id, int imol);
/*! \brief internal function for molecule centre

The centre is the mean position of the atoms of the molecule.

@param imol the model molecule index
@param iaxis the axis: 0 for x, 1 for y, 2 for z
@return the coordinate (in Å) of the centre along the given axis, or
a value less than -9999 for failure (e.g. bad imol or iaxis) */
float molecule_centre_internal(int imol, int iaxis);
/*! \brief a residue seqnum (normal residue number) from a residue
  serial number

   The serial number is the (0-based) index of the residue in the
   chain (in the first model).

   @param imol the model molecule index
   @param chain_id the chain id
   @param serial_num the residue serial number
   @return the residue number, < -9999 on failure */
int  seqnum_from_serial_number(int imol, const char *chain_id,
			       int serial_num);

/*! \brief the insertion code of the residue.

   @param imol the model molecule index
   @param chain_id the chain id
   @param serial_num the (0-based) index of the residue in the chain
   @return the insertion code (an empty string if the residue has
   none), NULL (scheme False) on failure. */
char *insertion_code_from_serial_number(int imol, const char *chain_id, int serial_num);

#ifdef __cplusplus
#ifdef USE_PYTHON
/*! \brief return a Python nested-list representation of (the first
  model of) molecule number imol

  The result is a list containing one model, which is a list of chains,
  each of the form [chain_id, residues]. Each residue is
  [res_no, ins_code, res_name, atoms] and each atom is
  [[atom_name, alt_conf], [occupancy, b_factor, element, seg_id], [x, y, z]],
  where b_factor is a list [B_iso, u11, u22, u33, u12, u13, u23] for
  anisotropic atoms.

  @return the list, or False if imol is not a valid model molecule */
PyObject *python_representation_kk(int imol);
#endif
#endif

/*! \brief the chain_id (string) of the ichain-th chain
  molecule number imol

   @param imol the model molecule index
   @param ichain the (0-based) chain index in the first model
   @return the chain-id, or False if the molecule or chain is not found */
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

@param imol the model molecule index
@return the number of models or -1 if there was a problem with the
given molecule.
*/
int n_models(int imol);


/*! \brief split an NMR model or other such multi-model molecule
 into multiple (separate) molecules - all in MODEL 1.

 @param imol the model molecule index
 @return the list of new molecule indices (an empty list on failure). */
#ifdef USE_PYTHON
PyObject *split_multi_model_molecule_py(int imol);
#endif


/*! \brief get the number of chains in molecule number imol

  Only the first model is considered.

  @param imol is the molecule index
  @return the number of chains, or -1 if imol is not a valid model molecule
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

   @return 0 for no (or on error: chain not found or imol not a valid
   model molecule), 1 for is "a nucleicacid chain".  Note that, unlike
   is_protein_chain_p(), this does not return -1 on error.

   This wraps the mmdb function isNucleotideChain().
   For completeness.
 */
int is_nucleotide_chain_p(int imol, const char *chain_id);


/*! \brief return the number of residues in the molecule,

All residues (including waters and ligands) of all models are counted.

@return the number of residues, -1 if this is a map or closed.
 */
int n_residues(int imol);

/*! \brief return the number of ATOMs in the molecule,

All models are counted. HETATMs (and TER records) are not counted.

@return the number of atoms, -1 if this is a map or closed.
 */
int n_atoms(int imol);


/* Does this work? */
/*! \brief return a list of the remarks of the molecule number imol

  Each item in the list is a [remark_number, remark_text] pair.

  @return the list of remarks (an empty list on failure)
  */
/* list remarks(int imol); */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM remarks_scm(int imol);
/* return a list or scheme false */
/*! \brief return the centre of the given residue

  The centre is the mean position of the atoms of the residue.

  @return a list [x, y, z] (in Å) or False if the residue was not found */
SCM residue_centre_scm(int imol, const char *chain_id, int resno, const char *ins_code);
#endif
#ifdef USE_PYTHON
/* return a list or python false */
/*! \brief return a list of the remarks of the molecule number imol

  Each item in the list is a [remark_number, remark_text] pair.

  @return the list of remarks, or False if imol is not a valid model molecule */
PyObject *remarks_py(int imol);
/*! \brief return the centre of the given residue

  The centre is the mean position of the atoms of the residue.

  @return a list [x, y, z] (in Å) or False if the residue was not found */
PyObject *residue_centre_py(int imol, const char *chain_id, int resno, const char *ins_code);
#endif
#endif


#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief model composition statistics - not yet implemented

  The statistics are calculated but not returned: this currently
  always returns an empty list. */
SCM model_composition_statistics_scm(int imol);
#endif

/*! \brief model composition statistics - not yet implemented

  The statistics are calculated but not returned: this currently
  always returns False. */
PyObject *model_composition_statistics_py(int imol);
#endif


/*! \brief sort the chain ids of the imol-th molecule in lexographical order */
void sort_chains(int imol);

/*! \brief sort the residues of the imol-th molecule

  The residues of each chain (in every model) are sorted into
  residue-number order. */
void sort_residues(int imol);

/*! \brief a gui dialog showing remarks header info (for a model molecule). */
void remarks_dialog(int imol);

/*! \brief simply print secondary structure info to the
  terminal/console.

  See get_header_secondary_structure_info() for a version that returns
  the info. */
void print_header_secondary_structure_info(int imol);

/*! \brief get the secondary structure from the header
 *
 * The HELIX and SHEET records of the first model are used.
 *
 * @param imol the molecule index
 * @return a dictionary of header info, or False if imol is not a valid
 * model molecule.
 * Returns: {'helices': [...], 'strands': [...]}
 * (a key is absent if there are no helices or no strands)
 *
 * Each helix dict contains:
 *  serNum, helixID, initChainID, initSeqNum, endChainID, endSeqNum, length, comment
 *  (comment only if set)
 *
 * Each strand dict contains:
 *  SheetID, strandNo, initChainID, initSeqNum, endChainID, endSeqNum
 */
PyObject *get_header_secondary_structure_info(int imol);


/*! \brief add secondary structure info to the
  internal representation of the model

  The secondary structure is calculated (by mmdb) and HELIX and SHEET
  records are generated from it. Nothing is done if the model already
  has helix or sheet records. */
void add_header_secondary_structure_info(int imol);


/*  Placeholder only.

    not documented, it doesn't work yet, because CalcSecStructure()
    creates SS type on the residues, it does not build and store
    CHelix, CStrand, CSheet records. */
void write_header_secondary_structure_info(int imol, const char *file_name);


/*! \brief copy molecule imol

Both model and map molecules can be copied.

@return the new molecule number.
Return -1 on failure to copy molecule (out of range, or molecule is
closed) */
int copy_molecule(int imol);

/*! \brief Copy a molecule with addition of a ligand and a deletion of
  current ligand.

  This function is used when adding a new (modified) ligand to a
  structure.  It creates a new molecule that is a copy of the current
  molecule except that the atoms of the current ligand/residue are
  replaced by (copies of) the atoms of the new ligand (and the residue
  takes the new ligand's residue name). Residues are specified by
  chain id and residue number only (no insertion code).

  @param imol_ligand_new the molecule index of the new ligand
  @param chain_id_ligand_new the chain id of the new ligand
  @param resno_ligand_new the residue number of the new ligand
  @param imol_current the molecule index of the molecule to be copied
  @param chain_id_ligand_current the chain id of the residue to be replaced
  @param resno_ligand_current the residue number of the residue to be replaced
  @return the index of the new molecule, or -1 on failure.
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
chain ids with the segids. Will that help?

The atoms of each model are regrouped into new chains, one for each
run of consecutive atoms with the same segid, and the new chain ids
are made from the segids. The original chains are deleted.

@param imol the model molecule index
@return 0 (the current implementation does not report whether
anything was changed). */
int exchange_chain_ids_for_seg_ids(int imol);

/*! \brief show the remarks browser

The remarks of the molecule chosen in the remarks browser molecule
chooser are shown (see remarks_dialog()). */
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
/*! \brief return the system build type string of this Coot build

  This is the build-system value \c COOT_SYS_BUILD_TYPE (made from the
  OS, system type, python and GTK versions), used for the updating of
  Coot. */
SCM coot_sys_build_type_scm();
#endif
#ifdef USE_PYTHON
/*! \brief return the system build type string of this Coot build

  This is the build-system value \c COOT_SYS_BUILD_TYPE (made from the
  OS, system type, python and GTK versions), used for the updating of
  Coot. */
PyObject *coot_sys_build_type_py();
#endif /* USE_PYTHON */
#endif /* c++ */

/*! \brief return the git revision count for for this build.
  */
int git_revision_count();
/*! \brief an alias to git_revision_count() for backwards compatibility  */
int svn_revision();


/*! \brief return the name of molecule number imol

 For a map molecule made from an MTZ file, the name includes the
 column labels, e.g. "d/e/f.mtz FWT PHWT"; for a model it is
 e.g. "/a/b/c.pdb".

 @return 0 if not a valid molecule ( -> False in scheme) */
const char *molecule_name(int imol);
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the molecule name without file extension

 Everything from the last ".pdb" in the name onwards is removed.

 @param imol the molecule index
 @param include_path_flag if 1, keep the directory part of the name,
 otherwise strip it
 @return the name stub, an empty string if imol is not a valid molecule */
SCM molecule_name_stub_scm(int imol, int include_path_flag);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the molecule name without file extension

 Everything from the last ".pdb" in the name onwards is removed.

 @param imol the molecule index
 @param include_path_flag if 1, keep the directory part of the name,
 otherwise strip it
 @return the name stub, an empty string if imol is not a valid molecule */
PyObject *molecule_name_stub_py(int imol, int include_path_flag);
#endif /* USE_PYTHON */
#endif	/* __cplusplus */
/*! \brief set the molecule name of the imol-th molecule */
void set_molecule_name(int imol, const char *new_name);
/*! \brief exit from coot, checking first for unsaved changes

 The command history is written. If there are no unsaved changes,
 coot exits (via coot_real_exit()); otherwise the unsaved-changes dialog is
 shown and coot does not exit.

 @param retval the exit status for the invoking process
 @return 1 (only reached when there were unsaved changes) */
int coot_checked_exit(int retval);
/*! \brief exit from coot, give return value retval back to invoking
  process.

  The state file and history are written first (see
  coot_save_state_and_exit()). */
void coot_real_exit(int retval);
/*! \brief exit without writing a state file

  The command history is written. */
void coot_no_state_real_exit(int retval);

/*! \brief exit coot doing clear-backup maybe

  In a Guile build this runs run_clear_backups(); otherwise it simply
  calls coot_real_exit(). */
void coot_clear_backup_or_real_exit(int retval);
/*! \brief exit coot, optionally writing a state file

  Waits for any running refinement to finish, then (if
  save_state_flag is set) saves the state file and the history, closes
  all molecules and exits.

  @param retval the exit status for the invoking process
  @param save_state_flag 1 to write the state file and history, 0 not to */
void coot_save_state_and_exit(int retval, int save_state_flag);


#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief run clear-backups

  Runs the scheme function (clear-backups-maybe), and exits coot
  (coot_real_exit()) if that did not need to show the backup dialog.

  @param retval the exit status for the invoking process */
void run_clear_backups(int retval);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief run clear-backups

  Runs the python function clear_backups_maybe(), and exits coot
  (coot_real_exit()) if that did not need to show the backup dialog.

  @param retval the exit status for the invoking process */
void run_clear_backups_py(int retval);
#endif /* USE_PYTHON */
#endif  /* c++ */


/*! \brief What is the molecule number of first coordinates molecule?

   @return -1 when there is none. */
int first_coords_imol();

/*! \brief molecule number of first small (<400 atoms) molecule.

@return -1 on no such molecule
  */
int first_small_coords_imol();

/*! \brief What is the molecule number of first unsaved coordinates molecule?

   @return -1 when there is none. */
int first_unsaved_coords_imol();

/*! \brief convert the structure factors in cif_file_name to an mtz
  file.

  Whichever of I, I_anom, F, F_anom, phi/fom, D, HL ABCD and the
  free-R flags are found in the mmCIF file are written to the MTZ file.

  @param cif_file_name the input mmCIF structure factors file name
  @param mtz_file_name the output MTZ file name
  @return 1 on success. Return 0 on a file without Rfree, return
  -1 on complete failure to write a file. */
int mmcif_sfs_to_mtz(const char *cif_file_name, const char *mtz_file_name);

/*! \} */

/*  -------------------------------------------------------------------- */
/*                    More Library/Utility Functions:                    */
/*  -------------------------------------------------------------------- */
/* section Graphics Utility Functions */
/*! \name Graphics Utility Functions */
/*! \{ */

/*! \brief set the bond lines to be antialiased

   @param state 1 for on, 0 for off (the default) */
void set_do_anti_aliasing(int state);
/*! \brief return the flag for antialiasing the bond lines */
int do_anti_aliasing_state();

/*! \brief turn the GL lighting on (state = 1) or off (state = 0)

   This function is no longer meaningful and does nothing.
*/
void set_do_GL_lighting(int state);
/*! \brief return the flag for GL lighting (0, since set_do_GL_lighting() no
  longer changes it) */
int do_GL_lighting_state();

/*! \brief shall we start up the Gtk and the graphics window?

   if passed the command line argument --no-graphics, coot will not start up gtk
   itself.

   An interface function for Ralf.

   @return 1 if the graphics interface is in use, 0 if not
*/
short int use_graphics_interface_state();

/*! \brief set the GUI dark mode state

 Sets the GTK "prefer dark theme" setting.

 @param state 1 for dark mode, 0 for not
 */
void set_use_dark_mode(short int state);

/*! \brief is the python interpreter at the prompt?

This is set by the --python command line option.

@return 1 for yes, 0 for no.*/
short int python_at_prompt_at_startup_state();

/*! \brief "Reset" the view

  Centre on the first displayed model molecule. If we are already
  centred there (within 0.1 Å), then centre on the next displayed
  model molecule (wrapping around to the first). The zoom is not
  changed and nothing happens if there are no displayed model molecules.

  @return 0 (the current implementation always returns 0). */
int reset_view();

/*! \brief set the view rotation scale factor

 Useful/necessary for high resolution displayed, where, without this factor
 the view doesn't rotate enough

 @param f the scale factor applied to mouse-drag view rotation (default 1.0) */
void set_view_rotation_scale_factor(float f);

/*! \brief return the number of molecules (coordinates molecules and
  map molecules combined) that are currently in coot

  @return the number of molecule slots. Closed molecules still occupy
  their slot and so are counted: this is one more than the highest
  molecule index, suitable as a loop limit (test each index with
  is_valid_model_molecule() or is_valid_map_molecule()). */
int get_number_of_molecules();

/*! \brief As above, return the number of molecules (coordinates molecules and
  map molecules combined) that are currently in coot.

  This is the old name for the function.

  @return the number of molecule slots (closed molecules are counted) */
int graphics_n_molecules();

/*! \brief does molecule number imol have hydrogen (or deuterium) atoms?

   The scripting interface to this does not have the _raw suffix and
   returns a scheme or python boolean True or False.

   @return either 1 (yes, there is at least one hydrogen) or 0 (no
   hydrogens, or no such molecule). */
int molecule_has_hydrogens_raw(int imol);

/* a testing/debugging function.  Used in a test to make sure that the
   outside number of a molecule (the vector index) is the same as that
   embedded in the molecule description object.  Return -1 on
   non-valid passed imol. */
int own_molecule_number(int imol);

/*! \brief Spin spin spin (or not)

  Toggles the continuous spinning of the view.
  See set_idle_function_rotate_angle(). */
void toggle_idle_spin_function();

/*! \brief Rock (not roll) (self-timed)

  Toggles the rocking of the view. See set_rocking_factors(). */
void toggle_idle_rock_function();
/* used by above to set the angle to rotate to (time dependent) */
/* (no longer used: currently always returns 0.0) */
double get_idle_function_rock_target_angle();


/*! \brief Settings for the inevitable discontents who dislike the
   default rocking rates (defaults 1 and 1)

   @param width_scale the scale factor for the rocking amplitude
   @param frequency_scale the scale factor for the rocking frequency */
void set_rocking_factors(float width_scale, float frequency_scale);

/*! \brief how far should we rotate when (auto) spinning? Fast
  computer? set this to 0.1

  @param f the spin speed factor (default 1.0). Although nominally in
  degrees, in the current implementation the view is rotated by
  0.004 * f radians per frame. */
void set_idle_function_rotate_angle(float f);  /* degrees */

/*! \brief what is the idle function rotation angle?

  @return the spin speed factor set by set_idle_function_rotate_angle() */
float idle_function_rotate_angle();

/*! \brief make a model molecule from the give file name.

* If the file updates, then the model will be updated (the file is
* checked every 500 ms).
*
* @param filename the coordinates file name
* @return 1 on success, 0 if the file could not be read as a model */
int make_updating_model_molecule(const char *filename);

/* or better still, use the json file from refmac ... but not yet. */
/* void updating_refmac_refinement_files(const char *updating_refmac_refinement_files_json_file_name); */

/* used by above, no API for this - also, not yet */
/* int updating_refmac_refinement_json_timeout_function(gpointer data); */


/*! \brief show the updating maps gui

this function is called from callbacks.c and calls a python gui function
*/
void show_calculate_updating_maps_pythonic_gui();

/*! \brief enable reading PDB/pdbx files with duplicate sequence numbers

(This is already on by default.) */
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
/*! \brief - was the input model read from an mmCIF file?

  @return true if the model was read from an mmCIF file, false if not
  or if there was an error with the molecule index */
SCM get_input_model_was_cif_state_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief - was the input model read from an mmCIF file?

  @return True if the model was read from an mmCIF file, False if not
  or if imol is not a valid model molecule */
PyObject *get_input_molecule_was_in_mmcif_state_py(int imol);
#endif
#endif


/*! \brief some programs produce PDB files with ATOMs where there
  should be HETATMs.  This is a function to assign HETATMs as per the
  PDB definition.

  The atoms of every residue whose name is not a PDB standard residue
  type are marked as HETATMs.

  @param imol the model molecule index
  @return the number of atoms in the non-standard residues (i.e. the
  atoms set as HETATMs), 0 if imol is not a valid model molecule */
int assign_hetatms(int imol);

/*! \brief if this is not a standard group, then turn the atoms to HETATMs.

@return 1 on atoms changes, 0 on not. Return -1 if residue not found
(or imol is not a valid model molecule).
*/
int hetify_residue(int imol, const char * chain_id, int resno, const char *ins_code);

/*! \brief residue has HETATMs?

@return 1 if any atom of the specified residue is a HETATM, else,
return 0.  If residue not found (or it has no atoms), return -1. */
int residue_has_hetatms(int imol, const char * chain_id, int resno, const char *ins_code);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief - get the specs for hetgroups - waters are not counted as het-groups.

  A het-group is a residue (in the first model) with at least one
  HETATM and a residue name other than HOH.

  @return a list of residue specs (an empty list on failure) */
SCM het_group_residues_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief - get the specs for hetgroups - waters are not counted as het-groups.

  A het-group is a residue (in the first model) with at least one
  HETATM and a residue name other than HOH.

  @return a list of residue specs, or False if imol is not a valid
  model molecule */
PyObject *het_group_residues_py(int imol);
#endif
#endif

/*! \brief return the number of non-hydrogen atoms in the given
  het-group (comp-id).

Return -1 on comp-id not found in dictionary.  */
int het_group_n_atoms(const char *comp_id);

/*! \brief replace the parts of molecule number imol that are
  duplicated in molecule number imol_frag

  The atoms of imol_fragment chosen by atom_selection are copied into
  imol_target: matching atoms in imol_target have their positions
  replaced, and atoms not found in imol_target are added (creating
  residues and chains as needed).

  @param imol_target the model molecule to be modified
  @param imol_fragment the model molecule from which atoms are taken
  @param atom_selection an mmdb atom selection string (for imol_fragment);
  several selections can be combined with "||"
  @return 1 on success, 0 on failure (e.g. invalid molecule indices) */
int replace_fragment(int imol_target, int imol_fragment, const char *atom_selection);

/*! \brief copy the given residue range from the reference chain to the target chain

resno_range_start and resno_range_end are inclusive. Residues that
already exist in the target chain have their atoms replaced; missing
residues are added.

@return 0 on failure (invalid molecule or chain). Note that the
current implementation also returns 0 on success. */
int copy_residue_range(int imol_target,    const char *chain_id_target,
		       int imol_reference, const char *chain_id_reference,
		       int resno_range_start, int resno_range_end);

/*! \brief replace the given residues from the reference molecule to the target molecule

The atoms of the given residues of imol_ref are copied into
imol_target (as for replace_fragment()).

@param imol_target the model molecule to be modified
@param imol_ref the model molecule from which the residues are taken
@param residue_specs_list_ref_scm a list of residue specs (for imol_ref)
@return 1 on success, 0 on failure
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

  The coordinates of model molecule molecule_number are replaced by
  those read from file_name.

  @param molecule_number the model molecule index
  @param file_name the coordinates file name
  @return -1 on failure.  Else return molecule_number  */
int clear_and_update_model_molecule_from_file(int molecule_number,
					      const char *file_name);

/* Used in execute_rigid_body_refine */
/* Fix this on a rainy day. */
/* atom_selection_container_t  */
/* make_atom_selection(int imol, const coot::minimol::molecule &mol);  */

/*! \brief dump the current screen image to a file.  Format tga
*
* make a copy of the screen image and write it to the file system
*
* You can use this, in conjunction with spinning and view moving functions to
* make movies
*
* @param tga_filename the output file name; ".tga" is appended if the
* name does not already end in ".tga" */
void screendump_image(const char *tga_filename);

/*! \brief give a warning dialog if density it too dark (blue)

  The dialog is shown (once) if the background is black and a
  displayed map is too blue. */
void check_for_dark_blue_density();

/* is this a good place for this function? */

/*! \brief sets the density map of the given molecule to be drawn as a
  (transparent) solid surface.

  @param imol the map molecule index
  @param state 1 for on, 0 for off */
void set_draw_solid_density_surface(int imol, short int state);

/*! \brief toggle for standard lines representation of map.

  This turns off/on standard lines representation of map.  transparent
  surface is another representation type.

  If you want to just turn off a map, don't use this, use
  set_map_displayed().

  @param imol the map molecule index
  @param state 1 for on, 0 for off
  */
void set_draw_map_standard_lines(int imol, short int state);

/*! \brief set the opacity of density surface representation of the
  given map.

0.0 is totally transparent, 1.0 (the default) is completely opaque and
(because the objects are no longer depth sorted) considerably faster
to render. 0.3 is a reasonable number.

@param imol the map molecule index
@param opacity the opacity, in the range 0.0 to 1.0
 */
void set_solid_density_surface_opacity(int imol, float opacity);

/*! \brief get the opacity of density surface representation of the
  given map.

@param imol the map molecule index
@return the opacity, or -1 if imol is not a valid map molecule */
float get_solid_density_surface_opacity(int imol);

/*! \brief set the flag to do flat shading rather than smooth shading
  for solid density surface.

Default is 1 (on).

@param state 1 for flat shading, 0 for smooth shading */
void set_flat_shading_for_solid_density_surface(short int state);

/*! \} */

/*  -------------------------------------------------------------------- */
/*                     Testing Interface:                                */
/*  -------------------------------------------------------------------- */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief run the built-in internal test suite (scheme interface)

  @return scheme true on success (or if Coot was built without
  BUILT_IN_TESTING), scheme false if a test failed */
SCM test_internal_scm();
/*! \brief run the single built-in internal test (scheme interface)

  @return scheme true on success (or if Coot was built without
  BUILT_IN_TESTING), scheme false if the test failed */
SCM test_internal_single_scm();
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief run the built-in internal test suite

  @return Python True on success (or if Coot was built without
  BUILT_IN_TESTING), False if a test failed */
PyObject *test_internal_py();
/*! \brief run the single built-in internal test (the one under development)

  @return Python True on success (or if Coot was built without
  BUILT_IN_TESTING), False if the test failed */
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

  To turn off mouse wheel recontouring call this with istate value of 0

  @param istate 1 for on (the default), 0 for off */
void set_scroll_by_wheel_mouse(int istate);
/*! \brief return the internal state of the scroll-wheel map contouring

  @return 1 for on, 0 for off */
int scroll_by_wheel_mouse_state();

/*! \brief turn off (0) or on (1) auto recontouring (on screen centre change) (default it on) */
void  set_auto_recontour_map(int state);

/*! \brief return the auto-recontour state

  @return 1 for on, 0 for off */
int get_auto_recontour_map();

/*! \brief set the default initial contour level for 2Fo-Fc-style maps

  This is applied to maps created after this call.

  @param n_sigma the contour level in multiples of the map r.m.s.d. (default 1.5) */
void set_default_initial_contour_level_for_map(float n_sigma);

/*! \brief set the default initial contour level for Fo-Fc-style (difference) maps

  This is applied to maps created after this call.

  @param n_sigma the contour level in multiples of the map r.m.s.d. (default 3.0) */
void set_default_initial_contour_level_for_difference_map(float n_sigma);

/*! \brief print the view matrix to the console

  Currently not functional: it no longer reflects the current view
  (the view is now held as a glm quaternion) and prints an identity
  matrix. */
void print_view_matrix();		/* print the view matrix */

/*! \brief get an element of the view matrix (used by the scripting view-matrix function)

  Currently not functional: it no longer reflects the current view
  and returns the element of an identity matrix.

  @param row the row index (0-3)
  @param col the column index (0-3)
  @return the matrix element */
float get_view_matrix_element(int row, int col); /* used in (view-matrix) command */

/*! \brief internal function to get an element of the view quaternion.
  The whole quaternion is returned by the scheme function
  view-quaternion

  Currently not functional: this always returns 0. */
float get_view_quaternion_internal(int element);

/*! \brief Set the view quaternion

  i, j and k are the vector (imaginary) part and l is the scalar
  (real) part. The quaternion is rejected (and the view not
  changed) if its magnitude is less than 0.5. The view is redrawn. */
void set_view_quaternion(float i, float j, float k, float l);

/*! \brief Given that we are in chain current_chain, apply the NCS
  operator that maps current_chain on to next_ncs_chain, so that the
  relative view is preserved.  For NCS skipping.

  Currently not functional: the implementation is disabled (pending
  a rewrite for the quaternion-based view) and this does nothing. */
void apply_ncs_to_view_orientation(int imol, const char *current_chain, const char *next_ncs_chain);
/*! \brief as above, but shift the screen centre also.

  Currently not functional: the implementation is disabled (pending
  a rewrite for the quaternion-based view) and this does nothing. */
void apply_ncs_to_view_orientation_and_screen_centre(int imol,
						     const char *current_chain,
						     const char *next_ncs_chain,
						     short int forward_flag);

/*! \brief set show frame-per-second flag

  @param t 1 for on, 0 for off (the default) */
void set_show_fps(int t);

/*! \brief the old name for set_show_fps() */
void set_fps_flag(int t);
/*! \brief return the state of the show frames-per-second flag

  @return 1 for on, 0 for off */
int  get_fps_flag();

/*! \brief set show frame-per-second flag */
void set_show_fps(int t);

/*! \brief set a flag: is the origin marker to be shown? 1 for yes, 0
  for no. (default 1) */
void set_show_origin_marker(int istate);
/*! \brief return the origin marker shown? state

  @return 1 for shown, 0 for not shown */
int  show_origin_marker_state();

/*! \brief hide the horizontal main toolbar

  Does nothing (other than printing a message) if there is no
  "main_toolbar" widget in the user interface. */
void hide_main_toolbar();
/*! \brief show the horizontal main toolbar

  Does nothing (other than printing a message) if there is no
  "main_toolbar" widget in the user interface. */
void show_main_toolbar();

/*! \brief reparent the Model/Fit/Refine dialog so that it becomes
  part of the main window, next to the GL graphics context

  No longer functional: this does nothing and returns 0. */
int suck_model_fit_dialog();
/*! \brief an alternative (handlebox) version of suck_model_fit_dialog()

  No longer functional: this does nothing and returns 0. */
int suck_model_fit_dialog_bl();

/*! \brief set the flag for "model-fit-refine dialog stays on top"

  @param istate 1 for on (the default), 0 for off */
void set_model_fit_refine_dialog_stays_on_top(int istate);
/*! \brief return the state model-fit-refine dialog stays on top */
int model_fit_refine_dialog_stays_on_top_state();

/* Legacy functions for the accept/reject dialog docking - no longer functional but
   retained for backwards compatibility with user startup scripts */
/*! \brief set the accept/reject dialog docked state - no longer functional */
void set_accept_reject_dialog_docked(int state);
/*! \brief set the accept/reject dialog docked show state - no longer functional */
void set_accept_reject_dialog_docked_show(int state);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                           mouse buttons                                  */
/*  ----------------------------------------------------------------------- */
/* section Mouse Buttons */
/*! \name   Mouse Buttons */
/*! \{ */

/*! \brief quanta-like buttons

  Swaps the roles of mouse buttons 1 and 2.

Note, when you have set these, there is no way to turn them off
   again (other than restarting). */
void quanta_buttons();
/*! \brief quanta-like zoom buttons

Note, when you have set these, there is no way to turn them off
   again (other than restarting). */
void quanta_like_zoom();


/* -------------------------------------------------------------------- */
/*    Ctrl for rotate or pick: */
/* -------------------------------------------------------------------- */
/*! \brief Alternate mode for rotation

Prefered by some, including Dirk Kostrewa.  I don't think this mode
works properly yet

  @param state 1 means Ctrl + left-mouse rotates the view (the
  default), 0 means Ctrl + left-mouse picks */
void set_control_key_for_rotate(int state);
/*! \brief return the control key rotate state */
int control_key_for_rotate_state();

/*! \brief Put the blob under the cursor to the screen centre.  Check only
positive blobs.  Useful function if bound to a key.

The refinement map must be set (and displayed).  (We can't check all maps because they
are not (or may not be) on the same scale).

Not useful for MCP. For interactive use only.

   @return currently always 0 (the result of the search is not passed back).
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
/*! \brief set the normal (pointer) cursor in the graphics window

  Cursor changing is not currently implemented in the GTK4 build. */
void normal_cursor();
/*! \brief fleur cursor

  Currently does nothing other than redraw. */
void fleur_cursor();
/*! \brief set the pick cursor, if Ctrl-for-rotate mode is on

  Cursor changing is not currently implemented in the GTK4 build. */
void pick_cursor_maybe();
/*! \brief rotate cursor (the same as normal_cursor()) */
void rotate_cursor();

/*! \brief let the user have a different pick cursor

sometimes (the default) GDK_CROSSHAIR is hard to see, let the user set
their own

  @param icursor_index a GdkCursorType value. This has no effect in the GTK4 build. */
void set_pick_cursor_index(int icursor_index);

/*! \} */


/*  --------------------------------------------------------------------- */
/*                      Model/Fit/Refine Functions:                       */
/*  --------------------------------------------------------------------- */
/* section Model/Fit/Refine Functions  */
/*! \name Model/Fit/Refine Functions  */
/*! \{ */

/*! \brief display the "Select Map for Fitting" (refinement map) chooser

  If no refinement map has been set, the first map molecule is made
  the refinement map. */
void show_select_map_frame();
/*! \brief Allow the changing of Model/Fit/Refine button label from
  "Rotate/Translate Zone" */
void set_model_fit_refine_rotate_translate_zone_label(const char *txt);
/*! \brief Allow the changing of Model/Fit/Refine button label from
  "Place Atom at Pointer" */
void set_model_fit_refine_place_atom_at_pointer_label(const char *txt);


/*! \brief shall atoms with zero occupancy be moved when refining? (default 1, yes)

  @param state 1 for yes, 0 for no */
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

 return 0 for backups off, 1 for backups on, -1 for unknown (e.g. imol is not a valid model molecule) */
int  backup_state(int imol);

/*! \brief apply undo - the "Undo" button callback
 *
 * undo the most recent modification on the model
 * set in set_undo_molecule() (or the only modified molecule, if
 * there is just one). If the choice of molecule is ambiguous, the
 * Undo Molecule chooser dialog is shown. Undo is not applied to an
 * undisplayed molecule.
 *
 * @return 1 on succesful undo, 0 on failed to undo.
 */
int apply_undo();		/* "Undo" button callback */

/*! \brief apply redo - the "Redo" button callback

  @return 1 on successful redo, 0 on failure to redo. */
int apply_redo();

/*! \brief set the molecule number imol to be marked as having unsaved changes */
void set_have_unsaved_changes(int imol);

/*! \brief does molecule number imol have unsaved changes?
 @return -1 on bad imol, 0 on no unsaved changes, 1 on has unsaved changes */
int have_unsaved_changes_p(int imol);

/*! \brief set the molecule to which undo operations are done to
  molecule number imol

  Ignored if imol is out of range or the molecule has no model. */
void set_undo_molecule(int imol);

/*! \brief show the Undo Molecule chooser - i.e. choose the molecule
  to which the "Undo" button applies. */
void show_set_undo_molecule_chooser();

/*! \brief set the state for adding paths to backup file names

  by default directories names are added into the filename for backup
  (with / to _ mapping).  call this with state=1 to turn off directory
  names

  Note: this flag is not currently used when making backup file names
  in Coot (it has no effect). See set_decoloned_backup_file_names(). */
void set_unpathed_backup_file_names(int state);
/*! \brief return the state for adding paths to backup file names*/
int  unpathed_backup_file_names_state();

/*! \brief set the state for "decoloned" backup file names

  When on, the backup file name is made from the molecule name as
  shown in the Display Manager (i.e. without the directory, unless
  paths are shown there) and colons in the time-stamp part of the
  name are replaced by underscores. This is on by default on
  Windows, off otherwise.

  @param state 1 for on, 0 for off */
void set_decoloned_backup_file_names(int state);
/*! \brief return the state for "decoloned" backup file names */
int  decoloned_backup_file_names_state();


/*! \brief return the state for compression of backup files*/
int  backup_compress_files_state();

/*! \brief set if backup files will be compressed or not using gzip

  @param state 1 for compressed (the default), 0 for uncompressed */
void  set_backup_compress_files(int state);

/*! \brief Make a backup for a model molecule
 *
 * @param imol the model molecule index
 * @param description a description that goes along with this backup point
 * @return the index of the backup (to be used in restore_to_backup_checkpoint()), or -1 on failure
 */
int make_backup_checkpoint(int imol, const char *description);

/*! \brief Restore molecule from backup
 *
 * restore model @p imol to checkpoint backup @p backup_index
 *
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to (as returned by make_backup_checkpoint())
 * @return the index of the backup, or -1 on failure
 */
int restore_to_backup_checkpoint(int imol, int backup_index);

#ifdef USE_PYTHON
/*! \brief Compare current model to backup
 *
 * @param imol the model molecule index
 * @param backup_index the backup index to compare with
 * @return a Python dict. The "status" key is either "ok",
 *         "fail" (imol is not a valid model molecule) or "bad-index".
 *         When the status is "ok", the other key is "moved-residues-list",
 *         the value for which is a list of residue specs for residues
 *         that have at least one atom in a different place (which might be empty).
 */
PyObject *compare_current_model_to_backup(int imol, int backup_index);
#endif

/*! \brief Print the backup history info for molecule imol to the console
 *
 * For each backup: the index, the molecule number, the backup file name,
 * the molecule name, the description and the time.
 */
void print_backup_history_info(int imol);

#ifdef USE_PYTHON
/*! \brief Get backup info
 *
 * @param imol the model molecule index
 * @param backup_index the backup index to query
 * @return a Python list of the given description (str)
 *         and a timestamp (str), or an empty list if imol is not
 *         a valid model molecule.
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

/*! \brief calculate phases using refmac and make a map

  Uses the first F and SIGF columns found in mtz_file_name, then
  fires up a GUI (via the scripting function
  refmac-for-phases-and-make-map) which asks which model molecule to
  calculate phases from; on selection, refmac is run and a map made.
  Does nothing if the file does not exist or has no F or SIGF columns. */
void calc_phases_generic(const char *mtz_file_name);

/*! \brief Calculate SFs (using refmac optionally) from an MTZ file
  and generate a map.

  Not implemented: this currently does nothing and always returns -1.

@return the new molecule number, -1 on a problem. */
int map_from_mtz_by_refmac_calc_phases(const char *mtz_file_name,
				       const char *f_col,
				       const char *sigf_col,
				       int imol_coords);


/*! \brief Calculate SFs from an MTZ file and generate a map.

 Structure factors are calculated from the model imol_coords and
 combined with the observed amplitudes to make a 2Fo-Fc-style map.

 @param mtz_file_name the MTZ file containing the observed data
 @param f_col the F column label
 @param sigf_col the SIGF column label
 @param imol_coords the model molecule from which phases are calculated
 @return the new molecule number, -1 on failure. */
int map_from_mtz_by_calc_phases(const char *mtz_file_name,
				const char *f_col,
				const char *sigf_col,
				int imol_coords);

#ifdef USE_PYTHON
/*! \brief Calculate structure factors (with bulk solvent correction) and update
   the given 2mFo-DFc map and Fo-Fc map

   Both imol_map_2fofc and imol_map_fofc must be existing valid maps (new maps
   are not created). The observed data are taken from imol_map_2fofc (which needs to
   have been read from an MTZ with Fobs/SigFobs/R-free data attached);
   imol_map_with_data_attached is currently not used for the calculation.
   The maps are recontoured at their current sigma levels and the R-factors are
   shown in the status bar.

   @return False on failure, otherwise a list of
   [r_factor, free_r_factor, bulk_solvent_volume, bulk_correction, table]
   where table is a list of [invresolsq, scale, lack_of_closure] items. */
PyObject *calculate_maps_and_stats_py(int imol_model,
                                      int imol_map_with_data_attached,
                                      int imol_map_2fofc,
                                      int imol_map_fofc);
#endif

/*! \brief Calculate structure factors from the model and update the given difference
           map accordingly

  @param imol_model the model molecule
  @param imol_map_with_data_attached a map molecule that has Fobs/SigFobs and R-free flags
         attached (i.e. read from an MTZ file with refmac parameters)
  @param imol_updating_difference_map the map to be overwritten; it must be a difference map */
void sfcalc_genmap(int imol_model, int imol_map_with_data_attached, int imol_updating_difference_map);

/*! \brief As above, calculate structure factors from the model and update the given difference
           map accordingly - but difference map gets updated automatically on modification of
           the imol_model molecule

  Only one such updating difference map can be active at a time. */
void set_auto_updating_sfcalc_genmap(int imol_model, int imol_map_with_data_attached, int imol_updating_difference_map);

/*! \brief As above, calculate structure factors from the model and update the given difference
           map accordingly - but the 2fofc and difference map get updated automatically on modification of
           the imol_model molecule

  imol_updating_difference_map must be a difference map. */
void set_auto_updating_sfcalc_genmaps(int imol_model, int imol_map_with_data_attached, int imol_updating_2fofc_map, int imol_updating_difference_map);


/* gdouble* get_map_colour(int imol); delete on merge 20220228-PE */

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the colours of map molecule imol

  Currently not functional (needs fixing): it always returns scheme false. */
SCM get_map_colour_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief return the colours of map molecule imol

  @return False if imol is not a valid map molecule, otherwise a list
  of two [red, green, blue] lists (values 0.0 to 1.0): the main map colour
  and the secondary (e.g. negative difference map) colour */
PyObject *get_map_colour_py(int imol);
#endif
#endif



/*! \brief set the map that is moved by changing the scroll wheel and
  change_contour_level().

  Ignored if imap is not a valid map molecule. */
void set_scroll_wheel_map(int imap);
/*! \brief set the map that has its contour level changed by the
  scrolling the mouse wheel to molecule number imol (same as set_scroll_wheel_map()). */
void set_scrollable_map(int imol);
/*! \brief the contouring of which map is altered when the scroll wheel changes?

  @return the molecule number of the scroll-wheel map */
int scroll_wheel_map();
/*! \brief save previous colour map for molecule number imol

  (used so that the map colour can be restored on Cancel) */
void save_previous_map_colour(int imol);
/*! \brief restore previous colour map for molecule number imol */
void restore_previous_map_colour(int imol);

/*! \brief set the immediate-map-update-on-map-drag state

By default, it is on (t=1).  On slower computers it might be better to
set t=0. */
void set_active_map_drag_flag(int t);
/*! \brief return the state of the dragged map flag  */
short int get_active_map_drag_flag();

/*! \brief set the colour of the last (highest molecule number) map

  The colour components are in the range 0.0 to 1.0 (values outside
  this range are clamped). */
void set_last_map_colour(double f1, double f2, double f3);

/*! \brief set the colour of the imolth map

  The colour components are in the range 0.0 to 1.0. */
void set_map_colour(int imol, float red, float green, float blue);

/*! \brief set the colour of the imolth map using a (7-character) hex colour

  @param hex_colour a colour string of the form "#rrggbb" (an invalid string gives grey)
  @param imol the map molecule index
*/
void set_map_hexcolour(int imol, const char *hex_colour);

/*! \brief  Make the maps 25% brighter */
void brighten_maps();

/*! \brief set the contour level, direct control

  @param imol_map the map molecule
  @param level the contour level in absolute map units (e.g. e/Å^3) */
void set_contour_level_absolute(int imol_map, float level);
/*! \brief set the contour level, direct control in r.m.s.d. (if you like that sort of thing)

  @param imol_map the map molecule
  @param level the contour level in multiples of the map r.m.s.d. */
void set_contour_level_in_sigma(int imol_map, float level);

/*! \brief get the contour level

  @return the contour level in absolute map units, or 0 if imol is not a valid map molecule */
float get_contour_level_absolute(int imol);

/*! \brief get the contour level in r.m.s.d. above 0.

  @return the contour level divided by the map r.m.s.d., or 0 if imol
  is not a valid map molecule */
float get_contour_level_in_sigma(int imol);

/*! \brief set the sigma step of the last map to f sigma

  This also turns on contouring by sigma steps for that map. */
void set_last_map_sigma_step(float f);
/*! \brief set the contour level step

   set the contour level step of molecule number imol to f and
   variable state (setting state to 0 turns off contouring by sigma
   level)

  @param imol the map molecule
  @param f the step in multiples of the map r.m.s.d. (ignored if state is 0)
  @param state 1 to contour by sigma steps, 0 to use the absolute step
  set by set_iso_level_increment() (or set_diff_map_iso_level_increment()) */
void set_contour_by_sigma_step_by_mol(int imol, float f, short int state);

/*! \brief return the resolution of the data for molecule number imol.
   Return negative number on error, otherwise resolution in A (eg. 2.0) */
float data_resolution(int imol);

/*! \brief return the resolution set in the header of the
  model/coordinates file.  If this number is not available, return a
  number less than 0.  */
float model_resolution(int imol);

/*! \brief export (write to disk) the map of molecule number imol to
  filename (in CCP4 map format).

  Return 0 on failure, 1 on success. */
int export_map(int imol, const char *filename);
/*! \brief export a fragment of the map about (x,y,z)

  Writes a CCP4 map file of the box of grid points within radius Å
  (along each axis) of the point (x,y,z), keeping the original cell.

  @return 1 if imol is a valid map molecule (even if the file
  writing failed), 0 otherwise */
int export_map_fragment(int imol, float x, float y, float z, float radius, const char *filename);

/*! \brief convenience function, called from callbacks.c

  Export a fragment of the map about the current rotation centre;
  radius_text is converted to an integer number of Å. */
void export_map_fragment_with_text_radius(int imol, const char *radius_text, const char *filename);

/*! \brief export a fragment of the map about (x,y,z), shifted to the origin

  Writes a CCP4 map file of a new box (of side 2*radius Å) with its
  origin at (0,0,0) and the point (x,y,z) at its centre. The density
  is masked to a sphere (of about 0.92 * radius) with a soft edge;
  outside that the density is set to 0.

  @return 1 if imol is a valid map molecule, 0 otherwise */
int export_map_fragment_with_origin_shift(int imol, float x, float y, float z, float radius, const char *filename);

/*! \brief export a fragment of the map about (x,y,z) as a plain-text file

  Writes the density values of the box of grid points within radius Å
  of (x,y,z) as text, 6 values per line, without a header.

  @return 1 if imol is a valid map molecule, 0 otherwise */
int export_map_fragment_to_plain_file(int imol, float x, float y, float z, float radius, const char *filename);

/*! \brief transform a map, making a new map molecule

  The density of map imol is transformed by the operator (rotation
  r00..r22 and translation t0,t1,t2 in Å) into a new map with the
  given space group and cell. The new density is generated in a box
  about the target point (pt0,pt1,pt2), extending box_half_size Å in
  each direction.

  @param imol the map molecule index
  @param r00 rotation matrix element (row 0, column 0)
  @param r01 rotation matrix element (row 0, column 1)
  @param r02 rotation matrix element (row 0, column 2)
  @param r10 rotation matrix element (row 1, column 0)
  @param r11 rotation matrix element (row 1, column 1)
  @param r12 rotation matrix element (row 1, column 2)
  @param r20 rotation matrix element (row 2, column 0)
  @param r21 rotation matrix element (row 2, column 1)
  @param r22 rotation matrix element (row 2, column 2)
  @param t0 the x component of the translation (in Å)
  @param t1 the y component of the translation (in Å)
  @param t2 the z component of the translation (in Å)
  @param pt0 the x coordinate of the target point (in Å)
  @param pt1 the y coordinate of the target point (in Å)
  @param pt2 the z coordinate of the target point (in Å)
  @param box_half_size the half-width of the box of new density (in Å)
  @param ref_space_group a H-M symbol or a colon-separated string of symmetry operators
  @param cell_a the cell length a in Å
  @param cell_b the cell length b in Å
  @param cell_c the cell length c in Å
  @param alpha the cell angle alpha in degrees
  @param beta the cell angle beta in degrees
  @param gamma the cell angle gamma in degrees
  @return the new molecule number, or -1 if imol is not a valid map molecule
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
  Return -1 on failure.

  The new map is marked as a difference map. */
int difference_map(int imol1, int imol2, float map_scale);

/*! \brief by default, maps that are P1 and have 90 degree angles
           are considered as maps without symmetry (i.e. EM maps).
           In some cases though P1 maps do/should have symmetry -
           and this is the means by you can tell Coot that.
    @param imol is the moleculle number to be acted on
    @param state the desired state, a value of 1 turns on map symmetry,
           0 turns it off (the map is treated as an EM map)
*/
void set_map_has_symmetry(int imol, int state);

/*! \brief make a new map (a copy of map_no) that is in the cell,
  spacegroup and gridding of the map in reference_map_no.

Return the new map molecule number - return -1 on failure */
int reinterp_map(int map_no, int reference_map_no);

/*! \brief make a new map (a copy of map_no) that is in the cell,
  spacegroup and a multiple of the sampling of the input map (a
  sampling factor of more than 1 makes the output maps smoother)

  @return the new map molecule number, or -1 on failure */
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
  averaged map, otherwise return the new map molecule number.

  The result is the weighted sum of the maps divided by the sum of the scale factors. */
int average_map_py(PyObject *map_number_and_scales);

/*! \brief Somewhat similar to the above function, except in this
case we overwrite the imol_map and we also presume that the
grid sampling of the contributing maps match. This makes it
much faster to generate than an average map.

The map values are replaced by the weighted sum of the
contributing maps (not divided by the sum of the scale factors).
This function does not itself recontour or redraw.
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

/*! \brief (GUI helper) return the iso-level increment as text

  imol is ignored. The returned string is malloc()ed and should be freed by the caller. */
char* get_text_for_iso_level_increment_entry(int imol); /* const gchar *text */
/*! \brief (GUI helper) return the difference-map iso-level increment as text

  imol is ignored. The returned string is malloc()ed and should be freed by the caller. */
char* get_text_for_diff_map_iso_level_increment_entry(int imol); /* const gchar *text */

/* void set_iso_level_increment(float val); */
/*! \brief set the contour scroll step (in absolute e/A3) for
  2Fo-Fc-style maps to val

The is only activated when scrolling by sigma is turned off

  (default 0.05) */
void set_iso_level_increment(float val);
/*! \brief return the contour scroll step (in absolute map units) for 2Fo-Fc-style maps */
float get_iso_level_increment();
/*! \brief (GUI helper) set the contour scroll step for 2Fo-Fc-style maps from text

  imol is ignored. */
void set_iso_level_increment_from_text(const char *text, int imol);

/*! \brief set the contour scroll step for difference map (in absolute
  e/A3) to val

The is only activated when scrolling by sigma is turned off

  (default 0.005) */
void set_diff_map_iso_level_increment(float val);

/*! \brief return difference maps iso-map level increment  */
float get_diff_map_iso_level_increment();
/*! \brief (GUI helper) set the difference maps iso-map level increment from text

  imol is ignored. */
void set_diff_map_iso_level_increment_from_text(const char *text, int imol);

/*! \brief set the map sampling rate from text

    This is for use by a GUI callback - not for user use. Values
    outside the range 1 to 100 are replaced by 1.5. */
void set_map_sampling_rate_text(const char *text);

/*! \brief set the map sampling rate (default 2.5)
 *
 * set_map_sampling_rate(float r)
 *
 * Set the Shannon Limit multiplier for map sampling (default 2.5).
 * This applies to maps made (e.g. from MTZ files) after this call.
 *
 * The Shannon (Nyquist) sampling theorem requires that a map be sampled
 * at least twice per resolution cycle (d_min/2) to faithfully represent
 * all frequencies present. This parameter sets the multiplier above that
 * theoretical minimum — so a value of 1.5 means maps are sampled at
 * 1.5 × the Shannon limit (i.e., grid spacing = d_min / 3.0).
 *
 * Higher values (e.g. 2.0–2.5) produce more finely sampled maps, which
 * can be useful for low-resolution baton-building or visual inspection,
 * since they are more visually attractive-looking,
 * at the cost of larger map files and slower computation.
 *
 * Set to something like 2.0 or 2.5 for more finely sampled maps.  Useful
 * for baton-building low resolution maps. */
void set_map_sampling_rate(float r);

/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief (GUI helper) return the map sampling rate as text

  The returned string is malloc()ed and should be freed by the caller. */
char* get_text_for_map_sampling_rate_text();

/*! \brief return the map sampling rate */
float get_map_sampling_rate();

/*! \brief change the contour level of the current (scroll-wheel) map by a step

if is_increment=1 the contour level is increased.  If is_increment=0
the map contour level is decreased.

The step is the absolute iso-level increment (see
set_iso_level_increment()). Note that for difference maps the level is
currently always increased (by the difference-map increment),
whatever the value of is_increment.
 */
void change_contour_level(short int is_increment); /* else is decrement.  */

/*! \brief set the contour level of the map with the highest molecule
    number to level (in absolute map units) */
void set_last_map_contour_level(float level);
/*! \brief set the contour level of the map with the highest molecule
    number to n_sigma sigma */
void set_last_map_contour_level_by_sigma(float n_sigma);

/*! \brief create a lower limit to the "Fo-Fc-style" map contour level changing

  When on, scrolling down will not take the difference map contour
  level below the level set by set_stop_scroll_diff_map_level().

  @param i 1 for on (the default), 0 for off */
void set_stop_scroll_diff_map(int i);
/*! \brief create a lower limit to the "2Fo-Fc-style" map contour level changing

  When on, scrolling down will not take the map contour level below
  the level set by set_stop_scroll_iso_map_level() (does not apply to
  Patterson maps).

  @param i 1 for on (the default), 0 for off */
void set_stop_scroll_iso_map(int i);

/*! \brief set the actual map level changing limit

   (in absolute map units, default 0.0) */
void set_stop_scroll_iso_map_level(float f);

/*! \brief set the actual difference map level changing limit

   (in absolute map units, default 0.0) */
void set_stop_scroll_diff_map_level(float f);

/*! \brief set the scale factor for the Residue Density fit analysis (default 1.0) */
void set_residue_density_fit_scale_factor(float f);


/*! \} */

/*  ------------------------------------------------------------------------ */
/*                         density stuff                                     */
/*  ------------------------------------------------------------------------ */
/* section Density Functions */
/*! \name  Density Functions */
/*! \{ */
/*! \brief draw the lines of the chickenwire density in width w (default 1) */
void set_map_line_width(int w);
/*! \brief return the width in which density contours are drawn */
int map_line_width_state();

/*! \brief make a map from an mtz file (simple interface)

 given mtz file mtz_file_name and F column f_col and phases column
 phi_col and optional weight column weight_col (pass use_weights=0 if
 weights are not to be used).  Also mark the map as a difference map
 (is_diff_map=1) or not (is_diff_map=0) because they are handled
 differently inside coot.

 A non-difference map becomes the scroll-wheel map.

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

 @param mtz_file_name the MTZ file name
 @param F_phi the F column label
 @param phi_col the phases column label
 @param weight_col the weights (e.g. FOM) column label
 @param use_weights 1 to use the weights column, 0 otherwise
 @param is_diff_map 1 if this is a difference map, 0 otherwise
 @param have_refmac_params currently ignored (the refmac parameters are always stored)
 @param fobs_col the Fobs column label
 @param sigfobs_col the SigFobs column label
 @param r_free_col the R-free flag column label
 @param sensible_f_free_col 1 if the R-free column is to be used, 0 otherwise

 @return -1 on error, else return imol */
int  make_and_draw_map_with_refmac_params(const char *mtz_file_name,
		                          const char *F_phi, const char *phi_col, const char *weight_col,
					  int use_weights, int is_diff_map,
					  short int have_refmac_params,
					  const char *fobs_col,
					  const char *sigfobs_col,
					  const char *r_free_col,
					  short int sensible_f_free_col);

/*! \brief as the above function, except set expert options too.

 @param mtz_file_name the MTZ file name
 @param a the F column label
 @param b the phases column label
 @param weight the weights (e.g. FOM) column label
 @param use_weights 1 to use the weights column, 0 otherwise
 @param is_diff_map 1 if this is a difference map, 0 otherwise
 @param have_refmac_params if 1, the refmac column labels (fobs_col, sigfobs_col,
        r_free_col, sensible_f_free_col) are stored with the map
 @param fobs_col the Fobs column label
 @param sigfobs_col the SigFobs column label
 @param r_free_col the R-free flag column label
 @param sensible_f_free_col 1 if the R-free column is to be used, 0 otherwise
 @param is_anomalous 1 if the map is an anomalous map (phases are shifted by 90 degrees)
 @param use_reso_limits 1 to apply the resolution limits, 0 to use all the data
 @param low_reso_limit the low resolution limit in Å
 @param high_reso_lim the high resolution limit in Å
 @return the new molecule number, or -1 on error
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

 If the file updates, then the map will be updated (the file is checked every 0.5 s).
 Use stop_updating_molecule() to stop.

 @return currently always 1 (not the molecule number) */
int make_updating_map(const char *mtz_file_name,
		      const char *f_col, const char *phi_col,
		      const char *weight,
		      int use_weights, int is_diff_map);


/*! \brief stop the updating of molecule imol

  For a map molecule, stop watching its MTZ file (see make_updating_map());
  for a model molecule, stop watching its coordinates file. */
void stop_updating_molecule(int imol);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the refmac parameters of map molecule imol

  @return a list (empty if imol is not a valid map or has no refmac
  parameters). See refmac_parameters_py(). */
SCM refmac_parameters_scm(int imol);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
/*! \brief return the refmac parameters of map molecule imol

  Used to test that the state of a refmac map is saved correctly.

  @return a list (empty if imol is not a valid map or has no refmac
  parameters) of: mtz file name, F column, phase column, weight
  column, use-weights flag, is-difference-map flag, 1, Fobs column,
  SigFobs column, R-free column, R-free-is-sensible flag */
PyObject *refmac_parameters_py(int imol);
#endif	/* USE_PYTHON */

#endif	/* __cplusplus */


/*! \brief does the mtz file have phases?

  @return 1 if the MTZ file has at least one phase column, 0 otherwise */
/* We need to know if an mtz file has phases.  If it doesn't then we */
/*  go down a (new 20060920) different path. */
int mtz_file_has_phases_p(const char *mtz_file_name);

/*! \brief is the given filename an mtz file?

  @return 1 if the file exists and has at least one F column, 0 otherwise */
int is_mtz_file_p(const char *filename);

/*! \brief does the given file have cns phases?

  @return 1 if the start of the file looks like CNS reflection data with phases, 0 otherwise */
int cns_file_has_phases_p(const char *cns_file_name);

/*! \brief auto-read an MTZ (or CNS) file, making maps from the
  recognised column labels (as for "Auto Open MTZ") */
void wrapped_auto_read_make_and_draw_maps(const char *filename);

/*! \brief set the flag to make a difference map (too) on auto-read MTZ

  @param i 1 for yes (the default), 0 for no */
void set_auto_read_do_difference_map_too(int i);
/*! \brief return the flag to do a difference map (too) on auto-read MTZ

   @return 0 means no, 1 means yes. */

int auto_read_do_difference_map_too_state();
/*! \brief set the expected MTZ columns for Auto-reading MTZ file.

  Not every program uses the default refmac labels ("FWT"/"PHWT") for
  its MTZ file.  Here we can tell coot to expect other labels so that
  coot can "Auto-open" such MTZ files. Each call adds another pair of
  labels to the list of those that are tried.

  e.g. (set-auto-read-column-labels "2FOFCWT" "PH2FOFCWT" 0) */
 void set_auto_read_column_labels(const char *fwt, const char *phwt,
				 int is_for_diff_map_flag);


/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief (GUI helper) return the x-ray map radius as text */
char* get_text_for_density_size_widget(); /* const gchar *text */
/*! \brief (GUI helper) set the x-ray map radius from text

  Values outside the range 0 to 1999.9 Å are ignored. */
void set_density_size_from_widget(const char *text);

/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief (GUI helper) return the EM map radius as text */
char *get_text_for_density_size_em_widget();
/*! \brief (GUI helper) set the EM map radius from text

  Values outside the range 0 to 19999.9 Å are ignored. */
void set_density_size_em_from_widget(const char *text);

/*! \brief set the extent of the box/radius of electron density contours for x-ray maps

  @param f the radius in Å (default 20) */
void set_map_radius(float f);

/*! \brief set the extent of the box/radius of electron density contours for EM map

  @param radius the radius in Å (default 80) */
void set_map_radius_em(float radius);

/*! \brief another (old) way of setting the radius of the map (the same as set_map_radius()) */
void set_density_size(float f);

/*! \brief set the maximum value of the map radius slider (default 50)

  This value is not currently used by the GUI. */
void set_map_radius_slider_max(float f);

/*! \brief Give me this nice message str when I start coot

  The message is shown in the status bar. */
void set_display_intro_string(const char *str);

/*! \brief return the extent of the box/radius of electron density contours for x-ray maps (in Å) */
float get_map_radius();

/*! \brief return the extent of the box/radius of electron density contours for EM maps (in Å) */
float get_map_radius_em();

/*! \brief not everyone likes coot's esoteric depth cueing system

  Pass an argument istate=0 to turn it off (it is on by default)

 (this function is currently disabled). */
void set_esoteric_depth_cue(int istate);

/*! \brief native depth cueing system

  return the state of the esoteric depth cueing flag */
int  esoteric_depth_cue_state();

/*! \brief not everone likes coot's default difference map colouring.

   Pass an argument i=1 to swap the difference map colouring so that
   red is positive and green is negative.

   This is used when map colours are next set; existing maps are not
   re-coloured by this call. */
void set_swap_difference_map_colours(int i);
/*! \brief return the state of the swap-difference-map-colours flag

  @return 1 if swapped (red is positive), 0 otherwise (the default) */
int swap_difference_map_colours_state();

/*! \brief post-hoc set the map of molecule number imol to be a
  difference map

  @param bool_flag 1 to mark the map as a difference map, 0 to mark it as a normal map
  @return success status, 0 -> failure (imol does not have a map)
  @param imol the map molecule index
*/
int set_map_is_difference_map(int imol, short int bool_flag);

/*! \brief map is difference map?

  @return 1 if imol is a valid map molecule that is a difference map, 0 otherwise */
int map_is_difference_map(int imol);

/*! \brief Add another contour level for the last added map.

  The map used is the refinement map, or, if that is not set, the
  non-difference map with the highest molecule number.

  Currently, the map must have been generated from an MTZ file.
  @return the molecule number of the new molecule or -1 on failure */
int another_level();

/*! \brief Add another contour level for the given map.

  The map is regenerated from its MTZ file as a new molecule, which is
  contoured 1 r.m.s.d. higher.

  Currently, the map must have been generated from an MTZ file.
  @return the molecule number of the new molecule or -1 on failure */
int another_level_from_map_molecule_number(int imap);

/*! \brief return the scale factor for the Residue Density fit analysis */
float residue_density_fit_scale_factor();

/*! \brief return the density at the given point for the given
  map. Return -999.9 for bad imol

  @param imol_map the map molecule
  @param x the x coordinate in Å
  @param y the y coordinate in Å
  @param z the z coordinate in Å
  @return the (linearly interpolated) density value */
float density_at_point(int imol_map, float x, float y, float z);

/*! \} */


/*  ------------------------------------------------------------------------ */
/*                         Parameters from map:                              */
/*  ------------------------------------------------------------------------ */
/* section Parameters from map */
/*! \name  Parameters from map */
/*! \{ */

/*! \brief return the mtz file that was used to generate the map

  @param imol_map the map molecule index
  @return a newly-allocated copy of the mtz file name. This is an
  empty string (not a null pointer) when \c imol_map is not a valid
  map molecule or the map was not generated from an mtz file (it was
  read from a CCP4 map file, say). Caller should dispose of the
  returned pointer. */
const char *mtz_hklin_for_map(int imol_map);

/*! \brief return the FP column in the file that was used to generate
  the map

  @param imol_map the map molecule index
  @return a newly-allocated copy of the amplitude column label. This
  is an empty string (not a null pointer) when \c imol_map is not a
  valid map molecule or there is no mtz file associated with that map
  (it was generated from a CCP4 map file, say). Caller should dispose
  of the returned pointer.
*/
const char *mtz_fp_for_map(int imol_map);

/*! \brief return the phases column in mtz file that was used to generate
  the map

  @param imol_map the map molecule index
  @return a newly-allocated copy of the phase column label. This is
  an empty string (not a null pointer) when \c imol_map is not a
  valid map molecule or there is no mtz file associated with that map
  (it was generated from a CCP4 map file, say). Caller should dispose
  of the returned pointer.
*/
const char *mtz_phi_for_map(int imol_map);

/*! \brief return the weight column in the mtz file that was used to
  generate the map

  @param imol_map the map molecule index
  @return a newly-allocated copy of the weight column label. This is
  an empty string (not a null pointer) when \c imol_map is not a
  valid map molecule, there is no mtz file associated with that map
  (it was generated from a CCP4 map file, say) or no weights were
  used. Caller should dispose of the returned pointer.
*/
const char *mtz_weight_for_map(int imol_map);

/*! \brief return flag for whether weights were used that was use to
  generate the map

  @param imol_map the map molecule index
  @return 1 if weights were used, 0 when no weights were used, there
  is no mtz file associated with that map or \c imol_map is not a
  valid map molecule. */
short int mtz_use_weight_for_map(int imol_map);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the parameters that made the map

@param imol the map molecule index
@return false if imol is not a valid map molecule, else a list of
  mtz file name, F column, phase column, weight column and use-weights
  flag, like ("xxx.mtz" "FPH" "PHWT" "" \#f) */
SCM map_parameters_scm(int imol);
/*! \brief return the unit cell of the molecule

@param imol a model or map molecule index
@return false (if imol is not a valid molecule or has no cell) or a
  list like (45 46 47 90 90 120), lengths in Angstroms, angles in
  degrees */
SCM cell_scm(int imol);
/* return the parameter of the molecule, something
  like (45 46 47 90 90 120), angles in degress */
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the parameters that made the map

@param imol the map molecule index
@return False if imol is not a valid map molecule, else a list of
  mtz file name, F column, phase column, weight column and use-weights
  flag, like ["xxx.mtz", "FPH", "PHWT", "", False] */
PyObject *map_parameters_py(int imol);
/*! \brief return the unit cell of the molecule

@param imol a model or map molecule index
@return False (if imol is not a valid molecule or has no cell) or a
  list like [45, 46, 47, 90, 90, 120], lengths in Angstroms, angles
  in degrees */
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

/*! \brief write molecule number imol as a PDB to file file_name

  If the file name has a SHELX extension (.ins, .res or .hat) a SHELX
  ins file is written instead (and the return value is then 1).

  @param imol the model molecule index
  @param file_name the output file name
  @return 0 on success, non-zero on error. Note that 0 is also
  returned if \c imol is not a valid model molecule (nothing is
  written). */
/*  return 0 on success, 1 on error. */
int write_pdb_file(int imol, const char *file_name);

/*! \brief write molecule number imol as a mmCIF to file file_name

  @param imol the model molecule index
  @param file_name the output file name
  @return 0 on success, non-zero on error. Note that 0 is also
  returned if \c imol is not a valid model molecule (nothing is
  written). */
/*  return 0 on success, 1 on error. */
int write_cif_file(int imol, const char *file_name);

/*! \brief write molecule number imol's residue range as a PDB to file
  file_name

  If \c resno_end is less than \c resno_start the two are swapped.
  CONECT records are removed unless the write-CONECT-records flag is
  set (see \c set_write_conect_record_state()).

  @param imol the model molecule index
  @param chainid the chain id
  @param resno_start the first residue number of the range
  @param resno_end the last residue number of the range
  @param filename the output file name
  @return 0 on success, non-zero on error (including when \c imol is
  not a valid model molecule). */
/*  return 0 on success, 1 on error. */
int write_residue_range_to_pdb_file(int imol, const char *chainid,
				    int resno_start, int resno_end,
				    const char *filename);

/*! \brief write the given chain of molecule number imol as a PDB to
  file filename

  @param imol the model molecule index
  @param chainid the chain id
  @param filename the output file name
  @return 0 on success, non-zero (an mmdb error code) on write
  failure. Note that 0 is also returned if \c imol is not a valid
  model molecule (nothing is written). */
/*  return 0 on success, -1 on error. */
int write_chain_to_pdb_file(int imol, const char *chainid, const char *filename);


/*! \brief save all modified coordinates molecules to the default
  names and save the state too.

  Each model molecule that has unsaved changes is written using its
  default save name. The state script is written to the XDG state
  directory ($XDG_STATE_HOME if set, otherwise
  ~/.local/state/Coot).

  @return 0 (always) */
int quick_save();

/*! \brief return the state of the write_conect_records_flag.

  @return 1 if CONECT records are written, 0 if not (the default) */
int get_write_conect_record_state();

/*! \brief set the flag to write (or not) conect records to the PDB file.

  This is used when saving coordinates from the Save Coordinates
  dialog and by \c write_residue_range_to_pdb_file().

  @param state 1 to write CONECT records, any other value for no (the
  default is no) */
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


/*! \brief created an ephemeral label in the graphics window
 *
 * the text stays on screen for about 2 seconds.
 *
 * @param txt the text
*/
void ephemeral_overlay_label(const char *txt);


/*! \} */


/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */
/* section Refmac Functions */
/*! \name  Refmac Functions */
/*! \{ */
/*! \brief set counter for runs of refmac so that this can be used to
  construct a unique filename for new output

  @param imol the molecule index
  @param refmac_count the new value of the counter */
void set_refmac_counter(int imol, int refmac_count);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the stored refmac SAD atom info

  Not implemented for Scheme: this always returns an empty list. */
SCM get_refmac_sad_atom_info_scm();
#endif /* GUILE */
#ifdef USE_PYTHON
/*! \brief return the stored refmac SAD atoms (to be used by refmac
  with the SAD option)

  @return a list of [atom_name, fp, fpp, lambda] entries, e.g.
  [["SE", -8.0, -4.0, None]], where unset values are None */
PyObject *get_refmac_sad_atom_info_py();
#endif /* PYTHON */
#endif /* c++ */

/*! \brief swap the colours of maps

  swap the colour of maps imol1 and imol2.  Useful to some after
  running refmac, so that the map to be build into is always the same
  colour

  Nothing happens unless both are valid map molecules. */
void swap_map_colours(int imol1, int imol2);
/*! \brief flag to enable above

  @param istate 1 to swap the pre- and post-refmac map colours so
  that the new map keeps the colour of the old one, 0 for off (the
  default) */
void set_keep_map_colour_after_refmac(int istate);

/*! \brief the keep-map-colour-after-refmac internal state

  @return 1 for "yes", 0 for "no"  */
int keep_map_colour_after_refmac_state();

/*! \brief test the refmac version

  Calls the scripting function \c get_refmac_version().

  @return 0 if refmac is older than 5.4 (or the version could not be
  determined), 1 if it is 5.4, 2 if it is 5.5 or newer (supports SAD
  and twin refinement) */
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
/*! \brief set the size of the displayed symmetry

  Symmetry-related atoms are shown within this radius of the screen
  centre. The symmetry of all molecules is updated and the graphics
  redrawn.

  @param f the radius in Angstroms (default 13) */
void set_symmetry_size(float f);
/*! \brief return the symmetry bonds colour

  @param imol ignored (the symmetry colour is not per-molecule)
  @return a newly-allocated array of 4 doubles of which the first 3
  (red, green, blue) are set. Caller should free it. */
double* get_symmetry_bonds_colour(int imol);
/*! \brief is symmetry master display control on?

  @return 1 for on, 0 for off */
short int get_show_symmetry(); /* master */
/*! \brief set display symmetry, master controller

  @param state 1 for on, 0 for off */
void set_show_symmetry_master(short int state);
/*! \brief set display symmetry for molecule number mol_no

   pass with state=0 for off, state=1 for on */
void set_show_symmetry_molecule(int mol_no, short int state);
/*! \brief display symmetry as CAs?


   pass with state=0 for off, state=1 for on */
void symmetry_as_calphas(int mol_no, short int state);
/*! \brief what is state of display CAs for molecule number mol_no?

   return state=0 for off, state=1 for on, -1 if imol is not a valid
   model molecule
*/
short int get_symmetry_as_calphas_state(int imol);

/*! \brief set the colour map rotation (i.e. the hue) for the symmetry
    atoms of molecule number imol

    @param imol the model molecule index
    @param state 1 for on, 0 for off */
void set_symmetry_molecule_rotate_colour_map(int imol, int state);

/*! \brief should there be colour map rotation (i.e. the hue) change
    for the symmetry atoms of molecule number imol?

   return state=0 for off, state=1 for on, -1 if imol is not a valid
   model molecule
*/
int symmetry_molecule_rotate_colour_map_state(int imol);

/*! \brief set symmetry colour by symop mode

  Colour the symmetry-related atoms according to their symmetry
  operator. This only has an effect when the graphics interface is in
  use.

  @param imol the model molecule index
  @param state 1 for on, 0 for off */
void set_symmetry_colour_by_symop(int imol, int state);
/*! \brief set symmetry display to show whole chains

  When on, symmetry-related chains that come near the screen centre
  are displayed in full (the "Display Near Chains" option), rather
  than only the atoms within the symmetry radius. This only has an
  effect when the graphics interface is in use.

  @param imol the model molecule index
  @param state 1 for on, 0 for off (the default) */
void set_symmetry_whole_chain(int imol, int state);
/*! \brief set use expanded symmetry atom labels

  @param state 1 for on, 0 for off (the default). In the current code
  this flag is only reflected in the Symmetry dialog and the saved
  state. */
void set_symmetry_atom_labels_expanded(int state);

/*! \brief molecule number imol has a unit cell?

   @return 1 on "yes, it has a cell", 0 for "no" */
int has_unit_cell_state(int imol);

/* a gui function really */
void add_symmetry_on_to_preferences_and_apply();

/*! \brief Undo symmetry view. Translate back to main molecule from
  this symmetry position.

  Uses the first molecule with symmetry displayed (see
  \c first_molecule_with_symmetry_displayed()).

  @return 0 (always) */
int undo_symmetry_view();

/*! \brief return the molecule number of the first model molecule
  that has a cell and symmetry and is displaying symmetry

@return -1 if there is no molecule with symmetry displayed.  */
int first_molecule_with_symmetry_displayed();

/*! \brief save the symmetry coordinates of molecule number imol to
  filename

Allow a shift of the coordinates to the origin before symmetry
expansion is applied (this is how symmetry works in Coot
internals): the coordinates are first translated by
-(pre_shift_to_origin_na, pre_shift_to_origin_nb,
pre_shift_to_origin_nc) unit cells, then symmetry operator
\c symop_no with the unit cell translation (shift_a, shift_b,
shift_c) is applied.

The file is written as mmCIF if the file name has an mmCIF
extension, otherwise as PDB.

@param imol the model molecule index
@param filename the output file name
@param symop_no the symmetry operator index (0 is the identity)
@param shift_a unit cell translation along a
@param shift_b unit cell translation along b
@param shift_c unit cell translation along c
@param pre_shift_to_origin_na pre-shift in unit cells along a
@param pre_shift_to_origin_nb pre-shift in unit cells along b
@param pre_shift_to_origin_nc pre-shift in unit cells along c */
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
coordinates (the translation in Angstroms).

Allow a shift of the coordinates to the origin before symmetry
expansion is applied: the operator is applied about the pre-shift
point, i.e. the coordinates are translated by -(na, nb, nc) unit
cells, transformed, then translated back by +(na, nb, nc) unit cells.

Pass "" as the name-in and a name will be constructed for you.

Return -1 on failure (e.g. imol is not a valid model molecule or has
no cell). */
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
coordinates (the translation in Angstroms).

Allow a shift of the coordinates to the origin before symmetry
expansion is applied (as for \c new_molecule_by_symmetry()).

Unlike \c new_molecule_by_symmetry(), the name is used as given (no
name is constructed if "" is passed).

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
  from imol by applying a symmetry operator

  The new molecule is named "SymOp_<symop>_Copy_of_<imol>".

  @param imol the model molecule index (it must have a cell)
  @param symop_string the symmetry operator in fractional coordinates,
  e.g. "-x+1/2,-y,z+1/2"
  @param pre_shift_to_origin_na pre-shift in unit cells along a (as for
  \c new_molecule_by_symmetry())
  @param pre_shift_to_origin_nb pre-shift in unit cells along b
  @param pre_shift_to_origin_nc pre-shift in unit cells along c
  @return the new molecule index, or -1 on failure
*/
int new_molecule_by_symop(int imol, const char *symop_string,
			  int pre_shift_to_origin_na,
			  int pre_shift_to_origin_nb,
			  int pre_shift_to_origin_nc);

/*! \brief return the number of symmetry operators for the given molecule

  @param imol a model or map molecule index

return -1 on no-symmetry for molecule or inappropriate imol number */
int n_symops(int imol);

/*! \brief move the chain of the reference molecule to the position
  of the symmetry-related copy

  Works on the symmetry atom closest to the screen centre: the chain
  containing that atom in the original molecule is moved to the
  symmetry-related position. Only works when the graphics interface
  is in use.

  @return 0 (always) */
/* This function works by active symm atom. */
int move_reference_chain_to_symm_chain_position();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the pre-shift (the shift that translates the centre
  of the molecule as close as possible to the origin) as a list of
  ints (unit cell translations along a, b and c) or scheme false on
  failure  */
SCM origin_pre_shift_scm(int imol);
#endif  /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the pre-shift (the shift that translates the centre
  of the molecule as close as possible to the origin) as a list of
  ints (unit cell translations along a, b and c) or Python False on
  failure  */
PyObject *origin_pre_shift_py(int imol);
#endif  /* USE_PYTHON */
#endif

/*! \brief start the "save symmetry coordinates" mode

  The user is asked to click on a symmetry atom. */
void setup_save_symmetry_coords();

/*! \brief set the space group for a coordinates molecule

 for shelx FA pdb files, there is no space group.  So allow the user
   to set it.  This can be initted with a HM symbol or a symm list for
   clipper.

This will only work on model molecules.

@return the success status of the setting  (1 good, 0 fail). Success
means that the space group of the molecule, read back after setting,
is identical to \c spg. */
short int set_space_group(int imol, const char *spg);

//! \brief set the unit cell for a given model molecule
//!
//! Angles in degress, cell lengths in Angstroms.
//!
//! @return  the success status of the setting (1 good, 0 fail).
//! Currently 1 is returned for any valid model molecule.
int set_unit_cell_and_space_group(int imol, float a, float b, float c, float alpha, float beta, float gamma, const char *space_group);

//! \brief set the unit cell and space group for a given model molecule using those of molecule imol_from
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
/*! \brief return the default file name suggestion (that would come
   up in the save coordinates dialog) or scheme false if imol is not a
   valid model molecule. */
SCM save_coords_name_suggestion_scm(int imol);
#endif /*  USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the default file name suggestion (that would come
   up in the save coordinates dialog) or Python False if imol is not a
   valid model molecule. */
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

  1 for on (the default), 0 for off. */
void set_console_display_commands_state(short int istate);
/*! \brief set a flag to show the text command equivalent of gui
  commands in the console as they happen in bold and colours.

  bold_flag: pass 1 for on (the default), 0 for off.

  colour_flag: pass  1 for on, 0 for off (the default).

  colour_index 0 to 7 inclusive for various different colourings
  (the ANSI terminal colours: 0 black, 1 red, 2 green, 3 yellow,
  4 blue (the default), 5 magenta, 6 cyan, 7 white).
 */
void set_console_display_commands_hilights(short int bold_flag, short int colour_flag, int colour_index);

/*! \} */

/*  --------------------------------------------------------------------- */
/*                  state (a graphics_info thing)                         */
/*  --------------------------------------------------------------------- */
/* info */
/*! \name State Functions */
/*! \{ */

/*! \brief scale up graphics - now available in scripting */
void scale_up_graphics();

/*! \brief scale down graphics - now available in scripting */
void scale_down_graphics();

/*! \brief save the current state to the default filename

  The state is written to the XDG state directory
  ($XDG_STATE_HOME if set, otherwise ~/.local/state/Coot) as
  0-coot.state.py and, in Guile-enabled builds, also as a Scheme
  script using the save-state file name. */
void save_state();

/*! \brief save the current state to file filename

  The state is written as a Scheme script in Guile-enabled builds,
  otherwise as a Python script. */
void save_state_file(const char *filename);

/*! \brief save the current state to file filename as a Python script */
void save_state_file_py(const char *filename);

/*! \brief set the default state file name (default 0-coot.state.scm
  in Guile-enabled builds, otherwise 0-coot.state.py) */
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
1: ask to run it (the default)
2: run it, no questions */
void set_run_state_file_status(short int istat);
/*! \brief run the state file (reading from default filenname)

  Runs 0-coot.state.scm (Guile-enabled builds) or 0-coot.state.py
  from the current directory, if it exists. */
void run_state_file();		/* just do it */
#ifdef USE_PYTHON
/*! \brief run the Python state file 0-coot.state.py from the current
  directory, if it exists */
void run_state_file_py();		/* just do it */
#endif /* USE_PYTHON */
/*! \brief run the state file depending on the state variables

  See \c set_run_state_file_status(). */
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

//! \brief set the virtual trackball mode
//!
//! @param mode 1 for "Flat", 2 for "Spherical Surface"
//!
void vt_surface(int mode);

//! \brief get the virtual trackball mode
//!
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

//! @param v in perspective mode, the distance from the camera to the
//! back clipping plane (ignored unless it is greater than 1.01 times
//! the eye z-position and less than 1000); in orthographic mode a
//! dimensionless scale factor (default 1.0)
void set_clipping_back(float v);

//! set clipping plane front - this goes in differnent directions for orthographics vs perspective
//! @param v in perspective mode, the distance from the camera to the
//! front clipping plane (ignored unless it is greater than 2 and
//! less than 0.99 times the eye z-position); in orthographic mode a
//! dimensionless scale factor (default 1.0)

void set_clipping_front(float v);

//! get clipping plane front
//! @return the perspective near-plane distance in perspective mode,
//! otherwise the orthographic front clipping factor
float get_clipping_plane_front();

//! get clipping plane back
//! @return the perspective far-plane distance in perspective mode,
//! otherwise the orthographic back clipping factor
float get_clipping_plane_back();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Unit Cell                                        */
/*  ----------------------------------------------------------------------- */
/* section Unit Cell interface */
/*! \name  Unit Cell interface */
/*! \{ */

/*! \brief return the state of show unit cell for molecule number imol */
//! @param imol the molecule index (it is not range-checked)
//! @return 1 for displayed, 0 for undisplayed
short int get_show_unit_cell(int imol);

//! set the state of show unit cell for all molecules
//! @param istate 1 for displayed, 0 for undisplayed
void set_show_unit_cells_all(short int istate);

//! set the state of show unit cell for the particular molecule number imol
//! @param imol is the molecule index
//! @param istate 1 for displayed, 0 for undisplayed
void set_show_unit_cell(int imol, short int istate);

//! set unit cell colour
//! @param red the red component (0 to 1)
//! @param green the green component (0 to 1)
//! @param blue the blue component (0 to 1)
void set_unit_cell_colour(float red, float green, float blue);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Colour                                           */
/*  ----------------------------------------------------------------------- */
/* section Colour */
/*! \name  Colour */
/*! \{ */

//! Set the symmetry colour merge
//!
//! The symmetry atom colour is molecule_colour * (1 - v) + symmetry_colour * v.
//! @param v the colour merge ratio (a fraction 0.0 to 1.0, default 0.5)
void set_symmetry_colour_merge(float v);

/*! \brief Set the hue change step on reading a new molecule */
//! Molecule n gets a hue rotation of (n+1) times this step.
//! @param f the hue change step in degrees (default 21)
void set_colour_map_rotation_on_read_pdb(float f);



/*! \brief should the hue change be applied to newly-read molecules?

  @param i 0 for no, 1 for yes (the default) */
void set_colour_map_rotation_on_read_pdb_flag(short int i);

/*! \brief shall the colour map rotation apply only to C atoms?

  Molecules currently coloured by chain are re-coloured.

  @param i 0 for no, 1 for yes (the default) */
void set_colour_map_rotation_on_read_pdb_c_only_flag(short int i);

/*! \brief Colour by chain */
//! @param imol the molecule index
void set_colour_by_chain(int imol);

/*! \brief colour molecule number imol by ncs chain type */
//! @param imol the molecule index
//! @param goodsell_mode 0 for no, 1 for yes
void set_colour_by_ncs_chain(int imol, short int goodsell_mode);

/*! \brief colour molecule number imol by chain type, goodsell-like colour scheme */
//! @param imol the molecule index
void set_colour_by_chain_goodsell_mode(int imol);

/*! \brief Set the goodsell chain colour colour wheel step  */
//! Molecules already drawn in goodsell mode are not re-coloured.
//! @param s the step size, default 0.221
void set_goodsell_chain_colour_wheel_step(float s);

/*! \brief Colour by molecule */
//! @param imol the molecule index
void set_colour_by_molecule(int imol);

/*! \brief get the colour-map-rotation-on-read-pdb C-only flag

  @return 1 if the hue rotation applies only to carbon atoms, else 0 */
/* get the value of graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag */
int get_colour_map_rotation_on_read_pdb_c_only_flag();

/*! \brief set the symmetry colour base */
//! The symmetry bonds are regenerated and redrawn.
//! @param r the red component (0 to 1)
//! @param g the green component (0 to 1)
//! @param b the blue component (0 to 1)
void set_symmetry_colour(float r, float g, float b);

/*! \} */

/*  Section Map colour*/
/*! \name   Map colour*/
/* \{ */
/*! \brief Set the colour map rotation (hue change) for maps

  Map molecule n is given a hue rotation of n times this step.

  @param f the hue change step in degrees, default for maps is 31 degrees. */
void set_colour_map_rotation_for_map(float f);

/*! \brief Set the colour map rotation

  The bonds of the molecule are regenerated and redrawn.

  @param imol the model molecule index
  @param theta is in degrees */
void set_molecule_bonds_colour_map_rotation(int imol, float theta);

/*! \brief Get the colour map rotation */
//! @param imol the molecule index
//! @return the rotation in degrees, or -1 if imol is not a valid model molecule
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

/*! \brief Get the aniso radius limit

  @return the radius in Angstroms (default 12). Note that the radius
  limit is not used by the current renderer. */
float get_limit_aniso();           /* not a function of the molecule */

/*! \brief Get the state of the aniso radius limit flag

  @return 1 for on, 0 for off (the default) */
short int get_show_limit_aniso();  /* not a function of the molecule */

/*! \brief return show-aniso-atoms state  - FIXME- per molecule

  This global flag can no longer be set (see \c set_show_aniso()), so
  it is always 0; use \c set_show_aniso_atoms() per molecule. */
short int get_show_aniso();       /*  not a function of the molecule */

/*! \brief Set the aniso atom radius limit flag

  Note that the radius limit is not used by the current renderer. */
//! @param state 0 for no, 1 for yes
void set_limit_aniso(short int state);

/*! \brief does nothing

  Use \c set_show_aniso_atoms() instead. */
void set_show_aniso(int state);

/*! \brief set show aniso atoms */
//! @param imol the molecule index
//! @param state 0 for no, 1 for yes
void set_show_aniso_atoms(int imol, int state);

/*! \brief set show aniso atoms as ortep */
//! Turning this on also turns on aniso atom display for the molecule.
//! @param imol the molecule index
//! @param state 0 for no, 1 for yes
void set_show_aniso_atoms_as_ortep(int imol, int state);

/*! \brief set show aniso atoms as empty ellipsoids

  The ellipsoids are drawn as (thicker) principal-axis bands around a
  small atom sphere, rather than as solid ellipsoids.

  @param imol the model molecule index
  @param state 0 for no, 1 for yes */
void set_show_aniso_atoms_as_empty(int imol, int state);

/* DELETE-ME */
void set_aniso_limit_size_from_widget(const char *text);

/* DELETE-ME .h */
char *get_text_for_aniso_limit_radius_entry();

/* DELETE-ME */
/*! \brief set the probability level of the displayed aniso ellipsoids

  @param f the probability as a fraction (default 0.5) */
void set_aniso_probability(float f);

/* DELETE-ME */
/*! \brief get the probability level of the displayed aniso ellipsoids

  @return the probability as a fraction (default 0.5) */
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

/*! \brief set the window size

  @param x_size the width in pixels
  @param y_size the height in pixels */
void   set_graphics_window_size(int x_size, int y_size);
/*! \brief set the window size as gtk_widget (flag=1) or gtk_window (flag=0) */
void   set_graphics_window_size_internal(int x_size, int y_size, int as_widget_flag);
/*! \brief set the graphics window position

  Not currently implemented for GTK4 (nothing is moved). */
void   set_graphics_window_position(int x_pos, int y_pos);
/*! \brief store the graphics window position */
void store_graphics_window_position(int x_pos, int y_pos); /*  "configure_event" callback */

/*! \brief store the graphics window position and size to
 *         xenops-graphics.scm and xenops-graphics.py in the
 *         preferences directory ($HOME/.coot). */
void graphics_window_size_and_position_to_preferences();

/*! \brief draw a frame */
void graphics_draw(); 	/* and wrapper interface to gtk_widget_draw(glarea)  */

/*! \brief try to turn on Zalman stereo mode  */
void zalman_stereo_mode();

void hardware_stereo_mode();

/*! \brief try to turn on stereo mode  */
void hardware_stereo_mode();

/*! \brief set the stereo mode (the relative view of the eyes)

0 is 2010-mode
1 is modern mode (the default)
*/
void set_stereo_style(int mode);

/*! \brief what is the stero state?

  @return the display mode: 0 for mono, 1 for hardware stereo, 2 for
  (cross-eyed) side by side stereo, 3 for DTI side by side stereo,
  4 for wall-eyed side by side stereo, 5 for Zalman stereo. */
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

/*! \brief set the stereo angle
 *
 * @param angle: stereo angle in degrees - default is 6 degrees
 * */
void set_stereo_angle(float angle);

/*! \brief return the hardware stereo angle factor

  Obsolete: this always returns 0. */
float hardware_stereo_angle_factor_state();

/*! \brief set the model display radius limit

  @param state 1 to limit the display of the model to within
  \c radius of the screen centre, 0 for no limit (the default)
  @param radius the radius in Angstroms (default 15) */
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

/*! \brief return the smooth scrolling state

  @return 1 for on, 0 for off */
int  get_smooth_scroll();

/* MOVE-ME to c-interface-gtk-widgets.h */
void set_smooth_scroll_steps_str(const char * t);

/*  useful exported interface */
/*! \brief set the number of steps in the smooth scroll

   Set more steps (e.g. 50) for more smoothness (default 20).*/
void set_smooth_scroll_steps(int i);

/* MOVE-ME to c-interface-gtk-widgets.h */
char  *get_text_for_smooth_scroll_steps();

/* MOVE-ME to c-interface-gtk-widgets.h */
void  set_smooth_scroll_limit_str(const char *t);

/*  useful exported interface */
/*! \brief do not scroll for distances greater this limit

  @param lim the limit in Angstroms (default 20) */
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

/*! \brief set use stroke characters

  @param state 1 for on, 0 for off (the default). The flag is stored
  but not used by the current renderer. */
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
/*! \brief set the rotation centre cross-hairs size scale factor from
  the text of a GUI entry, and redraw

  Values outside the range 0 to 1000 are rejected and 1.0 is used instead.

  @param text the entry text, interpreted as a number */
void set_rotation_centre_size_from_widget(const gchar *text); /* and redraw */
/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief return the current rotation centre cross-hairs size scale
  factor formatted as text (for a GUI entry)

  @return a newly allocated (malloc) string */
gchar *get_text_for_rotation_centre_cube_size();

/*! \brief set the rotation centre marker size

  This sets the same scale factor as
  \c set_user_defined_rotation_centre_crosshairs_size_scale_factor()
  and redraws.

  @param f the cross-hairs size scale factor (default 0.05) */
void set_rotation_centre_size(float f); /* and redraw (maybe) */

/*! \brief set the rotation centre marker (cross-hairs) size scale factor, and redraw

  @param f the scale factor (default 0.05) */
void set_user_defined_rotation_centre_crosshairs_size_scale_factor(float f);

/*! \brief set rotation centre colour

This is the colour for a dark background - if the background colour is not dark,
then the cross-hair colour becomes the inverse colour */
void set_rotation_centre_cross_hairs_colour(float r, float g, float b, float alpha);

/*! \brief return the recentre-on-pdb state

  @return 1 if the view is recentred on newly read coordinates, 0 if not */
short int recentre_on_read_pdb();
/*! \brief set the recentre-on-pdb state

  Should the view be centred on a newly read coordinates molecule?
  1 for yes, 0 for no (default 1). */
void set_recentre_on_read_pdb(short int);

/*! \brief set the rotation centre

  The view is moved to the new centre and redrawn.

  @param x the x coordinate in Angstroms
  @param y the y coordinate in Angstroms
  @param z the z coordinate in Angstroms */
void set_rotation_centre(float x, float y, float z);
/* The redraw happens somewhere else... */
/*! \brief set the rotation centre without redrawing (internal use) */
void set_rotation_centre_internal(float x, float y, float z);
/*! \brief return one component of the rotation centre

  @param axis 0 for x, 1 for y, 2 for z
  @return the coordinate in Angstroms (0.0 for any other value of axis) */
float rotation_centre_position(int axis); /* only return one value: x=0, y=1, z=2 */
/*! \brief centre on the ligand of the "active molecule", if we are
  already there, centre on the next hetgroup (etc) */
void go_to_ligand();

#ifdef USE_PYTHON
#ifdef __cplusplus
/*! \brief Python version of go_to_ligand()

  @return the new rotation centre as a list [x, y, z]. The view is
  moved by an animation that has only started when this returns.
  If no ligand was found the returned position is not meaningful. */
PyObject *go_to_ligand_py();
#endif
#endif

/*! \brief set the minimum number of atoms a residue must have to be
  considered a ligand by go_to_ligand()

  @param n_atom_min ligands must have at least this many atoms (default 6) */
void set_go_to_ligand_n_atoms_limit(int n_atom_min);

/*! \brief rotate the view so that the next main-chain atoms are oriented
 in the same direction as the previous - hence side-chain always seems to be
"up" - set this mode to 1 for reorientation-mode - and 0 for off (standard translation)

 The default is 0 (off).
*/
void set_reorienting_next_residue_mode(int state);

/*! \} */

/*  ---------------------------------------------------------------------- */
/*                         orthogonal axes                                 */
/*  ---------------------------------------------------------------------- */

/* section Orthogonal Axes */
/*! \name Orthogonal Axes */
/*! \{ */
/*! \brief draw the orthogonal axes in the top left of the graphics?

@param i 0 for off, 1 for on (default 1) */
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
/*! \brief get the model molecule list

  @return a list of the molecule indices of all the valid model (coordinates) molecules */
PyObject *get_model_molecule_list_py();
#endif
#endif

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief get the model molecule list

  @return a list of the molecule indices of all the valid model (coordinates) molecules */
SCM get_model_molecule_list_scm();
#endif
#endif

/* does not account for alternative conformations properly */
/* return -1 if atom not found. */
/*! \brief return the atom index of the specified atom

  The insertion code and alt conf are taken to be "".

  @return the index of the atom in the molecule's atom selection, or -1
  if the atom was not found */
int atom_index(int imol, const char *chain_id, int iresno, const char *atom_id);
/* using alternative conformations properly ?! */
/* return -1 if atom not found. */
/*! \brief return the atom index of the fully-specified atom

  @param imol the model molecule index
  @param chain_id the chain id
  @param iresno the residue number
  @param inscode the insertion code
  @param atom_id the atom name
  @param altconf the alt conf
  @return the index of the atom in the molecule's atom selection, or -1
  if the atom was not found */
int atom_index_full(int imol, const char *chain_id, int iresno, const char *inscode, const char *atom_id, const char *altconf);
/* Refine zone needs to be passed atom indexes (which it then converts */
/* to residue numbers - sigh).  So we need a function to get an
   atom. Return -1 on failure */
/* index from a given residue to use with refine_zone().  Return -1 on failure */
/*! \brief return the atom index of the first atom in the given residue

  @return the atom index, or -1 on failure (e.g. invalid molecule or
  residue not found) */
int atom_index_first_atom_in_residue(int imol, const char *chain_id,
				     int iresno, const char *ins_code);
/* For rotamers, we are given a residue spec (and altconf), we need
  the index of the first atom of this type, no atom name is given,
  hence we cannot use full_atom_spec_to_atom_index(). */
/*! \brief return the atom index of the first atom in the given residue
  that has the given alt conf

  @return the atom index, or -1 on failure */
int atom_index_first_atom_in_residue_with_altconf(int imol,
						  const char *chain_id,
						  int iresno,
						  const char *ins_code,
						  const char *alt_conf);
/*! \brief return the minimum residue number for imol chain chain_id

  @return the minimum residue number, or 999997 on failure (invalid
  molecule or chain not found) */
int min_resno_in_chain(int imol, const char *chain_id);
/*! \brief return the maximum residue number for imol chain chain_id

  @return the maximum residue number, or -99999 on failure (invalid
  molecule or chain not found) */
int max_resno_in_chain(int imol, const char *chain_id);
/*! \brief return the median temperature factor for imol

  @return the median atomic B-factor, or -1 if imol is not a valid model molecule */
float median_temperature_factor(int imol);
/*! \brief return the average temperature factor for the atoms in imol

  @return the mean atomic B-factor, or -1 if imol is not a valid model molecule */
float average_temperature_factor(int imol);
/*! \brief return the standard deviation of the atom temperature factors for imol

  @return the standard deviation, or -1 if imol is not a valid model molecule */
float standard_deviation_temperature_factor(int imol);

/*! \brief clear pending picks (stop coot thinking that the user is about to pick an atom).  */
void clear_pending_picks();
/*! \brief return the centre of mass of molecule imol as a string

  The format is "(x y z)" when Coot is built with Guile, otherwise
  "[x,y,z]".

  @return a newly allocated string, or NULL if imol is not a valid model molecule */
char *centre_of_mass_string(int imol);
#ifdef USE_PYTHON
/*! \brief return the centre of mass of molecule imol as a string in
  Python list format "[x,y,z]"

  @return a newly allocated string, or NULL if imol is not a valid model molecule */
char *centre_of_mass_string_py(int imol);
#endif
/*! \brief set the default temperature factor for newly created atoms
  (initial default 30) */
void set_default_temperature_factor_for_new_atoms(float new_b);
/*! \brief return the default temperature factor for newly created atoms */
float default_new_atoms_b_factor();

/*! \brief reset temperature factor for all moved atoms to the default
  for new atoms (usually 30)

  @param state 1 for on, 0 for off (default 0) */
void set_reset_b_factor_moved_atoms(int state);
/*! \brief return the state if temperature factors should be reset for
  moved atoms

  @return 1 for on, 0 for off */
int get_reset_b_factor_moved_atoms_state();

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief set the temperature factor of all the atoms in the specified residue

  No redraw is done.

  @param imol the model molecule index
  @param residue_spec_scm the residue spec
  @param bf the new B-factor
*/
void set_temperature_factors_for_atoms_in_residue_scm(int imol, SCM residue_spec_scm, float bf);
#endif
#endif

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief not yet implemented in the Scheme API - currently always returns an empty list */
SCM get_residue_alt_confs_scm(int imol, const char *chain_id, int res_no, const char *ins_code);
#endif
#endif

#ifdef __cplusplus /* protection from use in callbacks.c, else compilation probs */
#ifdef USE_PYTHON
/*! \brief Return either False (on failure) or a list of alt-conf strings (might be [""]) */
PyObject *get_residue_alt_confs_py(int imol, const char *chain_id, int res_no, const char *ins_code);
#endif
#endif



/*! \brief swap atom alt-confs

  Swap the alt conf labels of the first two atoms with the name atom_name
  in the specified residue. The molecule is redrawn.

  Note: the alt_conf argument is currently not used.

  @return currently always 0 */
int swap_atom_alt_conf(int imol, const char *chain_id, int res_no, const char *ins_code,
                       const char *atom_name, const char*alt_conf);

/*! \brief swap the alt-confs of all the atoms in the specified residue

  The residue must have exactly two different alt confs (an atom with no
  alt conf counts as one of them); the two labels are exchanged on every
  atom of the residue. A backup is made and the molecule is redrawn.

  @return currently always 0 */
int swap_residue_alt_confs(int imol, const char *chain_id, int res_no, const char *ins_code);

/*! \brief set a numerical attribute of the atom with the given specifier.

Attributes can be "x", "y", "z" (in Angstroms), "B" (or "b") and "occ"
and the attribute val is a floating point number. Only the first matching atom
is changed. The molecule is redrawn.

@return currently always 0 */
int set_atom_attribute(int imol, const char *chain_id, int resno, const char *ins_code, const char *atom_name, const char*alt_conf, const char *attribute_name, float val);

/*! \brief set a string attribute of the atom with the given specifier.

Attributes can be "atom-name", "alt-conf", "element" or "segid".
Only the first matching atom is changed. The molecule is redrawn.

@return currently always 0 */
int set_atom_string_attribute(int imol, const char *chain_id, int resno, const char *ins_code, const char *atom_name, const char*alt_conf, const char *attribute_name, const char *attribute_value);

/*! \brief set lots of atom attributes at once by-passing the rebonding and redrawing of the above 2 functions

  attribute_expression_list is a list of 8-element lists:
  (imol chain-id res-no ins-code atom-name alt-conf attribute-name attribute-value).
  String attribute values set "atom-name", "alt-conf", "element" or "segid";
  numerical values set "x", "y", "z", "B" (or "b") or "occ". Each molecule is
  backed up and rebonded once and there is a single redraw at the end.

  @return currently always 0 */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
int set_atom_attributes(SCM attribute_expression_list);
#endif

#ifdef USE_PYTHON
/*! \brief set lots of atom attributes at once (Python version)

  attribute_expression_list is a list of 8-element lists:
  [imol, chain_id, res_no, ins_code, atom_name, alt_conf, attribute_name, attribute_value].
  See set_atom_attributes().

  @return currently always 0 */
int set_atom_attributes_py(PyObject *attribute_expression_list);
#endif
#endif /* __cplusplus */

/*! \brief set the residue name of the specified residue

  The residue is renamed in all models. A backup is made and the molecule is redrawn. */
void set_residue_name(int imol, const char *chain_id, int res_no, const char *ins_code, const char *new_residue_name);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                            skeletonization                               */
/*  ----------------------------------------------------------------------- */
/* section Skeletonization Interface */
/*! \name  Skeletonization Interface */
/*! \{ */
/*! \brief turn on the (old) Greer skeleton drawing flag for the first
  non-difference map (legacy) */
void skel_greer_on();
/*! \brief turn off the (old) Greer skeleton drawing flag for all
  non-difference maps (legacy) */
void skel_greer_off();

/*! \brief skeletonize molecule number imol

   the prune_flag should almost  always be 0.

   NOTE:: The arguments to have been reversed for coot 0.8.3 and later
   (now the molecule number comes first).

   The skeleton level is set to the map mean + 1.5 rmsd. Nothing is done if
   the map already has a skeleton.

   @param imol the map molecule index
   @param prune_flag if non-zero, the skeleton is segmented by connectivity
   and coloured by segment; if 0, it is coloured by level
   @return currently always 0
    */
int skeletonize_map(int imol, short int prune_flag);

/*! \brief undisplay the skeleton on molecule number imol

   @return imol */
int unskeletonize_map(int imol);

/*! \brief if no map has yet been chosen for skeletonization, choose the
  first map molecule (for use by the GUI) */
void set_initial_map_for_skeletonize(); /* set graphics_info variable
					   for use in callbacks.c */

/*! \brief set the skeleton search depth, used in baton building

  For high resolution maps, you need to search deeper down the skeleton tree.  This
  limit needs to be increased to 20 or so for high res maps (it is 10 by default)

  @param v the search depth */
void set_max_skeleton_search_depth(int v); /* for high resolution
					      maps, change to 20 or
					      something (default 10). */



/*  ----------------------------------------------------------------------- */
/*                  skeletonization level widgets                           */
/*  ----------------------------------------------------------------------- */

/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief return the current skeletonization level as text (for a GUI entry)

  @return a newly allocated (malloc) string */
gchar *get_text_for_skeletonization_level_entry();

/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief set the skeletonization level from the text of a GUI entry

  Values that are not between 0 and 9999.9 are replaced by 0.2. The
  skeletons of all non-difference maps are updated and the graphics redrawn. */
void set_skeletonization_level_from_widget(const char *txt);

/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief return the current skeleton box size as text (for a GUI entry)

  @return a newly allocated (malloc) string */
gchar *get_text_for_skeleton_box_size_entry();

/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief set the skeleton box size from the text of a GUI entry

  Values that are not between 0 and 9999.9 are replaced by 0.2. */
void set_skeleton_box_size_from_widget(const char *txt);


/*! \brief the box size (in Angstroms) for which the skeleton is displayed

  The skeletons of all non-difference maps are updated.

  @param f the box size (radius) in Angstroms (default 40) */
void set_skeleton_box_size(float f);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                        save coordinates                                  */
/*  ----------------------------------------------------------------------- */
/* section Save Coordinates */
/*! \name  Save Coordinates */
/*! \{ */


/*! \brief save coordinates of molecule number imol in filename

  The format is chosen from the file name extension: mmCIF for
  mmCIF file names, a SHELX .ins/.res file for SHELX extensions,
  otherwise PDB format.

  @param imol the model molecule index
  @param filename the output file name
  @return 0 on success, non-zero on failure. Note that 0 is also returned
  if imol is not a valid model molecule (nothing is written). */
int save_coordinates(int imol, const char *filename);

/*! \brief set save coordinates in the starting directory

  Note: this flag is currently not used. */
void set_save_coordinates_in_original_directory(int i);

/* access to graphics_info_t::save_imol for use in callback.c */
/*! \brief return the molecule number selected for saving (GUI internal) */
int save_molecule_number_from_option_menu();
/* access from callback.c, not to be used in scripting, I suggest.
   Sets the *save* molecule number */
/*! \brief set the molecule number to be saved (GUI internal) */
void set_save_molecule_number(int imol);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                        .phs file reading                                 */
/*  ----------------------------------------------------------------------- */
/* section Read Phases File Functions */
/*! \name  Read Phases File Functions */
/*! \{ */

/*! \brief read phs file use coords to get cell and symm to make map

uses pending data to make the map: the phs file name is the one previously
stored with graphics_store_phs_filename(); the cell and symmetry are taken
from the coordinates file pdb_filename. A new map molecule is created (and
removed again if the cell or space group cannot be read).

*/
void
read_phs_and_coords_and_make_map(const char *pdb_filename);

/*! \brief read a phs file, the cell and symm information is from
  previously read (most recently read) coordinates file

 For use with phs data filename provided on the command line

 @return the new map molecule number, or -1 if there is no model molecule
 or the cell and symmetry could not be found */
int
read_phs_and_make_map_using_cell_symm_from_previous_mol(const char *phs_filename);


/*! \brief read phs file and use a previously read molecule to provide
  the cell and symmetry information

@return the new molecule number, return -1 if the cell and symmetry could not
be obtained from imol. Note that problems reading the phs file itself are not
detected (a new molecule index is still returned). */
int
read_phs_and_make_map_using_cell_symm_from_mol(const char *phs_filename, int imol);

/*! \brief as read_phs_and_make_map_using_cell_symm_from_mol() but the phs
  file name is the one previously stored with graphics_store_phs_filename()

@return the new molecule number, or -1 if the cell and symmetry could not be
obtained from imol */
int
read_phs_and_make_map_using_cell_symm_from_mol_using_implicit_phs_filename(int imol);

/*! \brief read phs file use coords to use cell and symm to make map

  The cell and space group are given explicitly.

  @param phs_file_name the phs file name
  @param hm_spacegroup the Hermann-Mauguin space group symbol, e.g. "P 21 21 21"
  @param a cell length in Angstroms
  @param b cell length in Angstroms
  @param c cell length in Angstroms
  @param alpha cell angle in degrees
  @param beta cell angle in degrees
  @param gamma cell angle in degrees
  @return the new map molecule number */
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

@param reso_lim_low the low resoluion limit (currently ignored).

@return the new map molecule number, or -1 on failure (the cell and symmetry
could not be obtained from imol, or the map could not be made). */
int
read_phs_and_make_map_with_reso_limits(int imol, const char* phs_file_name,
				       float reso_lim_low, float reso_lim_high);

/* work out the spacegroup from the given symm operators, e.g. return "P 1 21 1"
given "x,y,z ; -x,y+1/2,-z" */
/* char * */
/* spacegroup_from_operators(const char *symm_operators_in_clipper_format);  */

// 20220723-PE MOVE-ME!
/*! \brief store the phs file name for later use by the phs-reading
  functions that use an implicit phs file name */
void
graphics_store_phs_filename(const gchar *phs_filename);

/*! \brief could there be a cell and symmetry available for a phs file?

  @return 1 if any molecule has been read (or created), 0 if there are no molecules */
short int possible_cell_symm_for_phs_file();

/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief return a cell or symmetry field of molecule imol as text, for
  the phs cell chooser GUI

  @param imol the molecule index
  @param field one of "symm", "a", "b", "c", "alpha", "beta", "gamma"
  @return a newly allocated (malloc) string */
gchar *get_text_for_phs_cell_chooser(int imol, const char *field);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                                  Movement                                */
/*  ----------------------------------------------------------------------- */
/* section Graphics Move */
/*! \name Graphics Move */
/*! \{ */
/*! \brief undo last move

  Move the view (rotation centre) back to the previous rotation centre.  */
void undo_last_move(); /* suggested by Frank von Delft */

/*! \brief translate molecule number imol by (x,y,z) in Angstroms  */
void translate_molecule_by(int imol, float x, float y, float z);

/*! \brief transform molecule number imol by the given rotation
  matrix, then translate by (x,y,z) in Angstroms

  The matrix elements are given in row order (m11, m12, m13 is the first row). */
void transform_molecule_by(int imol,
			   float m11, float m12, float m13,
			   float m21, float m22, float m23,
			   float m31, float m32, float m33,
			   float x, float y, float z);

/*! \brief transform fragment of molecule number imol by the given rotation
  matrix, then translate by (x,y,z) in Angstroms

  The fragment is the residues of chain_id (in the first model) with residue
  numbers from resno_start to resno_end (inclusive) whose insertion code matches
  ins_code. The matrix elements are given in row order. A backup is made first. */
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
/*! \brief the go-to-atom chain-id

  @return a newly allocated (malloc) string */
char *go_to_atom_chain_id();
/*! \brief the go-to-atom atom name

  @return a newly allocated (malloc) string */
char *go_to_atom_atom_name();
/*! \brief the go-to-atom residue number */
int go_to_atom_residue_number();
/*! \brief the go-to-atom insertion code

  @return a newly allocated (malloc) string */
char *go_to_atom_ins_code();
/*! \brief the go-to-atom alt conf

  @return a newly allocated (malloc) string */
char *go_to_atom_alt_conf();


/*! \brief set the go to atom specification

   It seems important for swig that the `char *` arguments are `const
   char *`, not `const gchar *` (or else we get wrong type of argument
   error on (say) "A"

   The molecule is the current go-to-atom molecule (see set_go_to_atom_molecule()).
   The atom name may have an alt conf appended after a comma, e.g. " CA ,B".
   On success the view is centred on the atom.

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
/*! \brief set go to atom but don't redraw

   @param t1 the chain id
   @param iresno the residue number
   @param t3 the atom name, optionally with an alt conf after a comma
   @param make_the_move_flag if non-zero, centre on the atom; if 0, only set the
   go-to-atom specification
   @return 1 on success (always 1 if make_the_move_flag is 0), 0 for fail */
int set_go_to_atom_chain_residue_atom_name_no_redraw(const char *t1, int iresno, const char *t3,
						     short int make_the_move_flag);

// MOVE-ME!
/*! \brief as set_go_to_atom_chain_residue_atom_name() but the residue
  number t2 is given as a string */
int set_go_to_atom_chain_residue_atom_name_strings(const gchar *t1,
						   const gchar *t2,
						   const gchar *txt);


/*! \brief update the Go To Atom widget entries to atom closest to
  screen centre.

  The go-to-atom molecule and atom are set from the active atom
  (the atom closest to the screen centre) and the view is centred on it. */
void update_go_to_atom_from_current_position();


/* moving gtk function out of build functions, delete_atom() updates
   the go to atom atom list on deleting an atom  */
/*! \brief update the residue list of the Go To Atom window (currently does nothing) */
void update_go_to_atom_residue_list(int imol);

/*  return an atom index */
/*! \brief what is the atom index of the given atom?

  Any insertion code and alt conf match; the first matching atom is used.

  @return the atom index, or -1 if the atom was not found */
int atom_spec_to_atom_index(int mol, const char *chain, int resno, const char *atom_name);

/*! \brief what is the atom index of the given atom?

  @return the atom index, or -1 if the atom was not found */
int full_atom_spec_to_atom_index(int imol, const char *chain, int resno,
				 const char *inscode, const char *atom_name,
				 const char *altloc);

/*! \brief update the Go To Atom window

  Call this when molecule imol has changed (e.g. atoms added, deleted or renamed). */
void update_go_to_atom_window_on_changed_mol(int imol);

/*! \brief update the Go To Atom window.  This updates the option menu
  for the molecules. */
void update_go_to_atom_window_on_new_mol();

/*! \brief update the Go To Atom window when a different molecule (imol)
  has been chosen */
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
/*! \brief forget the stored Go To Atom window (GUI internal) */
void unset_go_to_atom_widget(); /* unstore the static go_to_atom_window */



/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  autobuilding control                                    */
/*  ----------------------------------------------------------------------- */
/* section AutoBuilding functions (Defunct) */
/* void autobuild_ca_on();  - moved to junk */

/*! \brief turn off the (defunct) CA autobuild flag */
void autobuild_ca_off();

/*! \brief developer test function */
void test_fragment();

/*! \brief prune the skeletons of all skeletonized non-difference maps,
  segmenting them by connectivity at the current skeleton level */
void do_skeleton_prune();

/*! \brief developer test function - what it does changes from time to time */
int test_function(int i, int j);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief developer test function */
SCM test_function_scm(SCM i, SCM j);
#endif
#ifdef USE_PYTHON
/*! \brief developer test function */
PyObject *test_function_py(PyObject *i, PyObject *j);
#endif /* PYTHON */
#endif


/*                    glyco tools test  */
/*! \brief developer test function for glyco trees (uses the active residue) */
void glyco_tree_test();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief Scheme version of glyco_tree_py() */
SCM glyco_tree_scm(int imol, SCM active_residue_scm);
/*! \brief Scheme version of glyco_tree_residues_py() */
SCM glyco_tree_residues_scm(int imol, SCM active_residue_scm);
/*! \brief Scheme version of glyco_tree_internal_distances_fn_py() (testing function) */
SCM glyco_tree_internal_distances_fn_scm(int imol, SCM residue_spec, const std::string &file_name); // testing function
/*! \brief Scheme version of glyco_tree_residue_id_py() */
SCM glyco_tree_residue_id_scm(int imol, SCM residue_spec_scm);
/*! \brief Scheme version of glyco_tree_compare_trees_py() */
SCM glyco_tree_compare_trees_scm(int imol_1, SCM res_spec_1, int imol_2, SCM res_spec_2);
/*! \brief Scheme version of glyco_tree_matched_residue_pairs_py() */
SCM glyco_tree_matched_residue_pairs_scm(int imol_1, SCM res_spec_1, int imol_2, SCM res_spec_2);
#endif
#ifdef USE_PYTHON
/*! \brief build the glyco tree from the given residue (incomplete)

  @return currently always False */
PyObject *glyco_tree_py(int imol, PyObject *active_residue_py);
/*! \brief return the residue specs of the residues in the glyco tree
  that contains the given residue

  @return a list of residue specs, or False if imol is not a valid model molecule */
PyObject *glyco_tree_residues_py(int imol, PyObject *active_residue_py);
/*! \brief write the internal distances of the glyco tree based on
  residue_spec to file_name (testing function)

  @return False */
PyObject *glyco_tree_internal_distances_fn_py(int imol, PyObject *residue_spec, const std::string &file_name); // testing function
/*! \brief return the glyco-tree identification of the given residue

  @return a list [level, prime-flag, res-type, link-type, parent-res-type,
  parent-residue-spec] where prime-flag is "prime", "non-prime" or "unset";
  or False on failure */
PyObject *glyco_tree_residue_id_py(int imol, PyObject *residue_spec_py);
/*! \brief compare the glyco trees based on the two given residues

  @return True if the trees match, False otherwise */
PyObject *glyco_tree_compare_trees_py(int imol_1, PyObject *res_spec_1, int imol_2, PyObject *res_spec_2);
/*! \brief return the matched residue pairs of the glyco trees based on
  the two given residues

  @return a list of [residue-spec-1, residue-spec-2] pairs, or False if
  there were no matches (or on failure) */
PyObject *glyco_tree_matched_residue_pairs_py(int imol_1, PyObject *res_spec_1, int imol_2, PyObject *res_spec_2);
#endif /* PYTHON */
#endif



/*  ----------------------------------------------------------------------- */
/*                  map and molecule control                                */
/*  ----------------------------------------------------------------------- */
/* section Map and Molecule Control */
/*! \name Map and Molecule Control */
/*! \{ */

/*! \brief display the Display Control window  */
void post_display_control_window();

/*! \brief regenerate the map entries of the Display Control window (GUI internal) */
void add_map_display_control_widgets();
/*! \brief regenerate the model entries of the Display Control window (GUI internal) */
void add_mol_display_control_widgets();
/*! \brief regenerate the model and map entries of the Display Control window (GUI internal) */
void add_map_and_mol_display_control_widgets();

/*! \brief forget the stored Display Control window (GUI internal) */
void reset_graphics_display_control_window();
/*! \brief hide the Display Control window (it is not destroyed) */
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
/*! \brief display only the maps in the given list

  pyo is a list of integers (map molecule numbers). Those maps are displayed
  and all other maps are undisplayed, with a single redraw. */
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
/*! \brief return the display state of map molecule number imol
 @return 1 for on, 0 for off (also 0 if imol is not a valid map molecule) */
int map_is_displayed(int imol);

/*! \brief if on_or_off is 0 turn off all maps displayed, for other
  values of on_or_off turn on all maps

  Note: this does nothing when Coot is running without graphics. */
void set_all_maps_displayed(int on_or_off);

/*! \brief if on_or_off is 0 turn off all models displayed and active,
  for other values of on_or_off turn on all models. */
void set_all_models_displayed_and_active(int on_or_off);

/*! \brief only display the last model molecule

  All other model molecules are undisplayed and made inactive.
*/
void set_only_last_model_molecule_displayed();

/*! \brief display only the active mol

  The model molecule of the active atom is displayed and made active; all other
  model molecules are undisplayed and made inactive. Maps are not changed. */
void display_only_active();


#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the spacegroup as a string, return scheme false if unable to do so. */
SCM space_group_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief return the spacegroup as a string, return Python False if unable to do so. */
PyObject *space_group_py(int imol);
#endif
#endif

/*! \brief return the spacegroup of molecule number imol . Deprecated.

@return "No spacegroup" when the spacegroup of a molecule has not been
set, an empty string if imol is not a valid model or map molecule. */
char *show_spacegroup(int imol);




#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief return a list of symmetry operators as strings - or scheme false if
  that is not possible. */
SCM symmetry_operators_scm(int imol);
/* take the return value from above and return a xHM symbol (for
   testing currently) */
/*! \brief convert a list of symmetry operator strings (as returned by
  symmetry_operators_scm()) to a Hermann-Mauguin space group symbol

  @return the symbol, or scheme false if no space group could be made */
SCM symmetry_operators_to_xHM_scm(SCM symmetry_operators);
#endif

#ifdef USE_PYTHON
/*! \brief return a list of symmetry operators as strings - or Python False if
  that is not possible. */
PyObject *symmetry_operators_py(int imol);
/* take the return value from above and return a xHM symbol (for
   testing currently) */
/*! \brief convert a list of symmetry operator strings (as returned by
  symmetry_operators_py()) to a Hermann-Mauguin space group symbol

  @return the symbol, or Python False if no space group could be made */
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
   molecule into which the others should be merged

   Molecules in add_molecules that are not valid model molecules, or are imol
   itself, are ignored. The merged-in molecules are undisplayed and made
   inactive (they are not closed). */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM merge_molecules(SCM add_molecules, int imol);
/*! \brief store a ligand residue spec for the merge-molecules code (currently stored but not used) */
void set_merge_molecules_ligand_spec_scm(SCM ligand_spec_scm);
#endif

#ifdef USE_PYTHON


/*! \brief merge molecules

@return a list, the first item of which is a status (1 is good); the
following items are the merge-infos (one for each merged item), e.g.
[1, "C", "D"]. If the molecule of an add_molecule item is just one residue,
the merge-info is a spec for the new residue, if it is many residues it is a
chain id.

the first argument is a list of molecule numbers and the second is the target
   molecule into which the others should be merged

   Molecules in add_molecules that are not valid model molecules, or are imol
   itself, are ignored. The merged-in molecules are undisplayed and made
   inactive (they are not closed). */
PyObject *merge_molecules_py(PyObject *add_molecules, int imol);
/*! \brief store a ligand residue spec for the merge-molecules code (currently stored but not used) */
void set_merge_molecules_ligand_spec_py(PyObject *ligand_spec_py);

/*! \brief split a multi-model ligand molecule (e.g. docking poses) and merge each
   conformer into its own copy of a protein molecule.

   imol_ligand is a multi-model ligand molecule; imol_protein is the protein.
   For each model (pose) in imol_ligand a fresh copy of imol_protein is made and
   the pose is merged into it, giving one protein+ligand complex molecule per pose.

   @return a list with one entry per pose: [imol_complex, merge_info], where
   merge_info is as returned by merge_molecules_py() ([status, spec-or-chain, ...])
   and describes where the ligand ended up in the complex. An empty list
   is returned if either molecule is not a valid model molecule. */
PyObject *split_multi_model_molecule_and_merge_py(int imol_ligand, int imol_protein);
#endif /* PYTHON */
#endif	/* c++ */


/*  ----------------------------------------------------------------------- */
/*                         Align and Mutate GUI                             */
/*  ----------------------------------------------------------------------- */
/* section Align and Mutate */
/*! \name  Align and Mutate */
/*! \{ */

/*! \brief align and mutate the given chain to the given sequence

  @param imol the model molecule index
  @param chain_id the chain to be mutated
  @param fasta_maybe the target sequence, either in FASTA format (starting
  with a "> name" line) or as plain sequence text
  @param renumber_residues_flag if non-zero, renumber the residues of the
  chain according to the alignment with the sequence */
void align_and_mutate(int imol, const char *chain_id, const char *fasta_maybe, short int renumber_residues_flag);
/*! \brief set the penalty for affine gap and space when aligning, defaults -3.0 and -0.4 */
void set_alignment_gap_and_space_penalty(float wgap, float wspace);


/* What are these functions?  consider deleting them - we have alignment_mismatches_* . */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief not implemented - always returns scheme false */
SCM alignment_results_scm(int imol, const char* chain_id, const char *seq);
/*! \brief return the residue spec of the nearest residue by sequence
  numbering.

@return  scheme false if not possible */
SCM nearest_residue_by_sequence_scm(int imol, const char* chain_id, int resno, const char *ins_code);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief not implemented - always returns Python False */
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
/*! \brief renumber the given residue range by offset residues

  The range is inclusive (start_res and last_res may be given in either order).
  Nothing is done if the renumbering would clash with residues outside the range.

  @return 1 if at least one residue was renumbered, 0 otherwise */
int renumber_residue_range(int imol, const char *chain_id,
			   int start_res, int last_res, int offset);


/*! \brief change residue number and insertion code for given
  residue

  @return 1 if imol is a valid model molecule (the change was attempted),
  -1 otherwise */
int change_residue_number(int imol, const char *chain_id, int current_resno, const char *current_inscode, int new_resno, const char *new_inscode);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                         Change chain id                                  */
/*  ----------------------------------------------------------------------- */
/* section Change Chain ID */
/*! \name Change Chain ID */
/*! \{ */

/*! \brief change the chain id of the specified residues

  @param imol the model molecule index
  @param from_chain_id the current chain id
  @param to_chain_id the new chain id
  @param use_res_range_flag if 0, change the chain id of the whole chain
  (which fails if to_chain_id already exists in the molecule); if 1, only
  the residues from from_resno to to_resno are moved to chain to_chain_id
  @param from_resno the first residue of the range
  @param to_resno the last residue of the range */
void  change_chain_id(int imol, const char *from_chain_id, const char *to_chain_id,
		      short int use_res_range_flag, int from_resno, int to_resno);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief as change_chain_id() but return the result

  @return a list (status message) where status is 1 on success and 0 on
  failure, or scheme false if imol is not a valid model molecule */
SCM change_chain_id_with_result_scm(int imol, const char *from_chain_id, const char *to_chain_id,
                                         short int use_res_range_flag, int from_resno, int to_resno);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief as change_chain_id() but return the result

  @return a list [status, message] where status is 1 on success and 0 on
  failure, or False if imol is not a valid model molecule */
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

@return 1 for yes, 0 for no, -1 for not yet known. Note that nothing
currently sets this state, so -1 is returned. */
int probe_available_p();
#ifdef USE_PYTHON
/*! \brief Python version of probe_available_p() */
int probe_available_p_py();
#endif

/*! \brief do nothing - compatibility function

  (If Coot is built with Guile, this calls post_scheme_scripting_window().) */
void post_scripting_window();

/*! \brief pop-up a scripting window for scheming */
void post_scheme_scripting_window();


/* called from c-inner-main */
/*! \brief run the scripts, commands and accession-code fetches given on
  the command line (internal), then clear them so they are not run again */
void run_command_line_scripts();

/*! \brief note that the Guile GUI scripts have been loaded (internal) */
void set_guile_gui_loaded_flag();
/*! \brief note that the Python GUI scripts have been loaded (internal) */
void set_python_gui_loaded_flag();
/*! \brief note that the Scheme scripting GUI code was found and loaded (internal) */
void set_found_coot_gui();
/*! \brief note that the Python scripting GUI code was found and loaded (internal) */
void set_found_coot_python_gui();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Monomer                                                 */
/*  ----------------------------------------------------------------------- */
/* section monomers */

/*! \name Monomer */
/*! \{ */

/*! \brief create a new molecule from (idealised coordinates of) a
  dictionary entry

  @param dict_idx the index of the entry in the monomer dictionary
  @param imol_enc the molecule number for which the dictionary entry is
  defined (or the "any molecule" value)
  @return the new molecule index, or -1 on failure */
int get_monomer_for_molecule_by_index(int dict_idx, int imol_enc);


/*  Don't let this be seen by standard c, since I am using a std::string */
/*  and now we make it return a value, which we can decode in the calling
    function. I make a dummy version for when GUILE is not being used in case
    there are functions in the rest of the code that call safe_scheme_command
    without checking if there is USE_GUILE first.
    importing ability for python modules, use_namespace to maintain the
    namespace of the module, returns 1 if not running, 0 on success, -1
    when error importing (no further information) */

/*! \brief run script file

  If the file name ends in ".py" it is run as a Python script, otherwise
  as a Scheme script. */
void run_script       (const char *filename);
/*! \brief guile run script file */
void run_guile_script (const char *filename);
/*! \brief run python script file */
void run_python_script(const char *filename);
/*! \brief import python module

  @param module_name the module name
  @param use_namespace if non-zero use "import module_name", otherwise
  "from module_name import *"
  @return 0 on success, -1 on error importing, 1 if Python is not available */
int import_python_module(const char *module_name, int use_namespace);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return a list of compoundIDs in the dictionary of which the
  given string is a substring of the compound name */
SCM matching_compound_names_from_dictionary_scm(const char *compound_name_fragment,
						short int allow_minimal_descriptions_flag);

/*! \brief return the monomer name

  If needed, an attempt is made to load the dictionary for comp_id first.

  return scheme false if not found */
SCM comp_id_to_name_scm(const char *comp_id);
#endif /* USE_GUILE */

/*! \brief try to auto-load the dictionary for comp_id from the refmac monomer library.

   return 0 on failure.
   return 1 on successful auto-load.
   return 2 on already-read.
   */
int auto_load_dictionary(const char *comp_id);
/*! \brief as above, but dictionary is re-read even if it already exists

   @return 1 on success, 0 on failure */
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

  Unlike the Scheme version, no attempt is made to load a dictionary.

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

/*! \brief start (or cancel) interactive regularize-zone atom picking

  With state 1, Coot switches to the pick cursor and waits for the
  user to click 2 atoms (in the same molecule) that define the zone
  to be regularized.

  @param state 1 to start picking, 0 to cancel (restore the normal cursor) */
void do_regularize(short int state); /* pass 0 for off (unclick togglebutton) */

/*! \brief start (or cancel) interactive refine-zone atom picking

  With state 1, Coot waits for the user to click 2 atoms (or 1 atom
  and then press the A key for an auto-zone). If no refinement map
  has been set, the map-selection frame is shown first; if there is
  still no refinement map, a warning dialog is shown and picking
  is not started.

  @param state 1 to start picking, 0 to cancel */
void do_refine(short int state);

/*! \brief add a restraint on peptides to make them planar

  This adds a 5 atom plane restraint (CA, C, O of the first residue
  and N, CA of the second; esd 0.08 Å) to the TRANS and PTRANS
  links, so it includes both CA atoms of the peptide.  Use this
  rather than editing the mon_lib_list.cif file. */
void add_planar_peptide_restraints();

/*! \brief remove the 5-atom planar peptide restraints (added by
  add_planar_peptide_restraints()) from the TRANS and PTRANS links. */
void remove_planar_peptide_restraints();

/*! \brief make the planar peptide restraints tight

  Sets the esd of the 5-atom planar peptide restraint of the TRANS
  link to 0.03 Å. The restraint must already exist (see
  add_planar_peptide_restraints()).

Useful when refining models with cryo-EM maps */
void make_tight_planar_peptide_restraints();


/*! \brief query the state of the planar peptide restraints

  @return 1 if the TRANS link contains the 5-atom planar peptide
  restraint, 0 if not */
int planar_peptide_restraints_state();


/*! \brief add a restraint on peptides to keep trans peptides trans

i.e. omega in trans-peptides is restrained to 180 degrees.

@param on_off_state 1 for on, 0 for off (default on)
 */
void set_use_trans_peptide_restraints(short int on_off_state);

/*! \brief add restraints on the omega angle of the peptides

  (that is the torsion round the peptide bond).  An omega torsion
  restraint (esd 5 degrees) is added to the dictionary links:
  180 degrees for TRANS and PTRANS links, 0 degrees for CIS and PCIS
  links, so cis-linked peptides are refined as cis and trans-linked
  peptides (the normal case) as trans. */
void add_omega_torsion_restraints();

/*! \brief remove omega restraints on CIS and TRANS linked residues. */
void remove_omega_torsion_restriants();

/*! \brief add or remove auto H-bond restraints

  When on, hydrogen bond restraints are generated automatically
  during refinement.

  @param state 1 for on, 0 for off (default off) */
void set_refine_hydrogen_bonds(int state);


/*! \brief set immediate replacement mode for refinement and regularization
 *
 * This can enable synchronous refinement (with istate = 1).
 * You need this (call with istate=1) if you are
 * scripting refinement/regularization: the refine/regularize
 * functions then wait for the refinement to finish before
 * returning, and some (e.g. refine_residues_py()) also accept the
 * new coordinates immediately.
 *
 * @param istate 1 for on, 0 for off (default 0)
 * */
void set_refinement_immediate_replacement(int istate);

/*! \brief query the state of the immediate replacement mode

  @return 1 for on, 0 for off */
int  refinement_immediate_replacement_state();

/*! \brief set "noughties physics" for dragged-atom refinement

  In this mode, dragging an atom of the intermediate atoms moves the
  other atoms by a shear-style displacement (or, with Ctrl pressed,
  just the dragged atom) rather than running the threaded refinement
  as the atom is dragged.

  @param state 1 for on, 0 for off (default off) */
void set_refine_use_noughties_physics(short int state);

/*! \brief query the "noughties physics" state

  @return 1 for on, 0 for off */
int get_refine_use_noughties_physics_state();

/*! \brief set the number of frames for which the selected residue
  range flashes

 On fast computers, this can be set to higher than the default (2) for
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

   Discard the intermediate (moving) atoms and their graphics object
   without accepting them, i.e. the model is not changed.
 */
void clear_up_moving_atoms();	/* remove the molecule and bonds */

/*! \brief remove just the bonds (the graphics object) of the moving atoms

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

   or \#f if restraints could not be set up (e.g. no refinement map
   has been set or no residues were found).
 */
SCM refine_residues_scm(int imol, SCM r); /* presumes the alt_conf is "". */
/*! \brief refine residues (a list of residue specs) using the given alt conf

  @return as refine_residues_scm() */
SCM refine_residues_with_alt_conf_scm(int imol, SCM r, const char *alt_conf);
/*! \brief refine residues using the given alt conf

  The mode arguments are currently ignored; this is the same as
  refine_residues_with_alt_conf_scm(). */
SCM refine_residues_with_modes_with_alt_conf_scm(int imol, SCM residues_spec_list_scm,
						 const char *alt_conf,
						 SCM mode_1,
						 SCM mode_2,
						 SCM mode_3);
/*! \brief regularize residues, r is a list of residue specs (the alt conf is "")

  Regularization uses geometric restraints only (no map).

  @return refinement results as for refine_residues_scm(), or \#f */
SCM regularize_residues_scm(int imol, SCM r); /* presumes the alt_conf is "". */
/*! \brief regularize residues (a list of residue specs) using the given alt conf

  @return refinement results as for refine_residues_scm(), or \#f */
SCM regularize_residues_with_alt_conf_scm(int imol, SCM r, const char *alt_conf);
#endif
#ifdef USE_PYTHON
/*! \brief refine the residues in the given residue spec list
 *
 * @param imol is the molecule index
 * @param rl is a Python list of residue specs, where a residue spec is a list of
 *  [`chain_id`, `res_no`, `ins_code`]
 *
 *  When using this function from scripting, make sure that
 *  set_refinement_immediate_replacement(1) is called first.
 *
 *  @return False if restraints could not be set up (or no refinement map has been
 *  set, see set_imol_refinement_map()), or a list of 3 elements on success:
 *  - [0] info_text (string): refinement information text
 *  - [1] progress (int): refinement progress indicator.
 *          0: GSL_SUCCESS: means refinement was successfully completed
 *         -2: GSL_CONTINUE: means refinement was successful, but didn't terminate (so more cycles needed)
 *         27: GSL_ENOPROG: iteration is not making progress towards solution
 *  - [2] lights (list or False): False if empty, otherwise a list of
 *    [`name`, `label`, `value`] triples where `name` and `label` are strings
 *    and `value` is a float  */
PyObject *refine_residues_py(int imol, PyObject *rl);  /* presumes the alt_conf is "". */

/*! \brief refine the residues in the given residue spec list, with modes and alt conf

  If mode_1 is the string "soft-mode/hard-mode", no refinement is
  done (that mode is not currently implemented) and False is
  returned; mode_2 and mode_3 are ignored. Otherwise this is the
  same as refine_residues_with_alt_conf_py().

  @return as refine_residues_py() */
PyObject *refine_residues_with_modes_with_alt_conf_py(int imol, PyObject *r, const char *alt_conf,
						      PyObject *mode_1,
						      PyObject *mode_2,
						      PyObject *mode_3);
/*! \brief refine the residues in the given residue spec list using the given alt conf

  @return as refine_residues_py() */
PyObject *refine_residues_with_alt_conf_py(int imol, PyObject *r, const char *alt_conf);
/*! \brief regularize the residues in the given residue spec list (the alt conf is "")

  Regularization uses geometric restraints only (no map).

  @return False if restraints could not be set up, otherwise a list as
  described for refine_residues_py() */
PyObject *regularize_residues_py(int imol, PyObject *r);  /* presumes the alt_conf is "". */
/*! \brief regularize the residues in the given residue spec list using the given alt conf

  @return as regularize_residues_py() */
PyObject *regularize_residues_with_alt_conf_py(int imol, PyObject *r, const char *alt_conf);
#endif /* PYTHON */
#endif /* c++ */

/*! \brief stop a running (threaded) refinement

  Used by on_accept_reject_refinement_reject_button_clicked(). Waits
  until the refinement has stopped; it does not clear up the moving
  atoms. */
void stop_refinement_internal();

/*! \brief use soft (harmonic approximation) non-bonded contact restraints

  @param flag 1 for on, 0 for off (default off) */
void set_refinement_use_soft_mode_nbc_restraints(short int flag);

/*! \brief shiftfield B-factor refinement

  Refine the B-factors of model molecule imol using the observed
  data (Fobs, SigFobs and R-free flags) associated with the current
  refinement map (see set_imol_refinement_map()), so that map must
  have been made from reflection data.

  @param imol the model molecule index */
void shiftfield_b_factor_refinement(int imol);

/*! \brief shiftfield xyz refinement

  Not implemented yet - this function does nothing. */
void shiftfield_xyz_factor_refinement(int imol);

/*! \brief turn on (or off) torsion restraints

   Pass with istate=1 for on, istate=0 for off (default off).
*/
void set_refine_with_torsion_restraints(int istate);
/*! \brief return the state of torsion restraints in refinement (1 for on, 0 for off) */
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

/*! \brief estimate a suitable map weight for refinement against the given map

  The estimate is 15/rmsd of the map (multiplied by 0.35 for EM maps).
  The current map weight is not changed - use set_matrix() for that.

  @param imol_map the map molecule index
  @return the estimated weight, or -1 if imol_map is not a valid map */
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

/*! \brief refine a zone based on atom indexing

  @param imol the model molecule index
  @param ind1 the index (in the molecule's atom selection) of an atom in the first residue
  @param ind2 the index of an atom in the last residue */
void refine_zone_atom_index_define(int imol, int ind1, int ind2);

/*! \brief refine a zone

 Refine the residues from resno1 to resno2 (inclusive) in the given
 chain (insertion codes are assumed to be ""). Both end residues must
 exist.

 presumes that imol_Refinement_Map has been set

 @param altconf the alt conf of the atoms to refine ("" for none)
 @param imol the model molecule index
 @param chain_id the chain id
 @param resno1 the first residue number of the zone
 @param resno2 the last residue number of the zone
*/
void refine_zone(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
/*! \brief repeat the previous (user-selected) refine zone */
void repeat_refine_zone(); /* use stored atom indices to re-run the refinement using the same atoms as previous */
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief refine a zone and return the refinement results

  @return refinement results as for refine_residues_scm(), or \#f if
  either end residue was not found or restraints could not be set up */
SCM refine_zone_with_score_scm(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
/*! \brief regularize a zone and return the refinement results

  altconf is currently ignored.

  @return refinement results as for refine_residues_scm(), or \#f on failure */
SCM regularize_zone_with_score_scm(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
#endif /* guile */
#ifdef USE_PYTHON
/*! \brief refine a zone and return the refinement results

  Refine the residues from resno1 to resno2 (inclusive) in the given
  chain (insertion codes are assumed to be "") against the
  refinement map.

  @return refinement results as for refine_residues_py(), or False if
  either end residue was not found or restraints could not be set up */
PyObject *refine_zone_with_score_py(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
/*! \brief regularize a zone and return the refinement results

  altconf is currently ignored.

  @return refinement results as for refine_residues_py(), or False on failure */
PyObject *regularize_zone_with_score_py(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);
#endif /* PYTHON */
#endif /* c++ */

/*! \brief refine a zone using auto-range

 Refine the residues around residue resno1 (whose CA atom is used),
 using the auto-range step (see set_refine_auto_range_step()).

 presumes that imol_Refinement_Map has been set */
void refine_auto_range(int imol, const char *chain_id, int resno1, const char *altconf);

/*! \brief regularize a zone

 Regularize (geometry only, no map) the residues from resno1 to
 resno2 in the given chain (insertion codes are assumed to be "").
 altconf is currently ignored.

@return a status, whether the regularisation was done or not.  0 for no, 1 for yes.
  */
int regularize_zone(int imol, const char *chain_id, int resno1, int resno2, const char *altconf);

/*! \brief set the number of refinement steps applied to the
  intermediate atoms each frame of graphics.

  smaller numbers make the movement of the intermediate atoms slower,
  smoother, more elegant.

  Default: 20. */
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

Default 0.25

 Bigger numbers mean bigger movement of the other atoms.*/
void set_refinement_drag_elasticity(float e);

/*! \brief turn on Ramachandran angles refinement in refinement and regularization

  Also shows (or hides) the "Rama" indicator in the main toolbar.

  @param state 1 for on, 0 for off (default off) */
/* name consistent with set_refine_with_torsion_restraints() !?  */
void set_refine_ramachandran_angles(int state);
/*! \brief turn on Ramachandran angles refinement - an alias for
  set_refine_ramachandran_angles() with a better name

  @param state 1 for on, 0 for off */
void set_refine_ramachandran_torsion_angles(int state);

/*! \brief change the target function type

  @param type 0 for the "ZO" Ramachandran restraints (this also resets
  the Ramachandran restraints weight to 1.0), 1 for the log(Ramachandran
  probability) restraints (the default) */
void set_refine_ramachandran_restraints_type(int type);
/*! \brief change the target function weight

The default is 1.0. A big number means bad things (the refinement
may fail to converge). */
void set_refine_ramachandran_restraints_weight(float w);

/*! \brief ramachandran restraints weight

@return weight as a float */
float refine_ramachandran_restraints_weight();

/*! \brief set the weight for torsion restraints (default 1.0)*/
void set_torsion_restraints_weight(double w);

/*! \brief set the state for using rotamer restraints "drive" mode

1 in on, 0 is off (off by default) */
void set_refine_rotamers(int state);

/*! \brief set the Geman-McClure alpha from text (used by the refinement parameters dialog)

  @param combobox_item_idx the index of the item in the dialog's combobox
  @param t the value as text */
void set_refinement_geman_mcclure_alpha_from_text(int combobox_item_idx, const char *t);
/*! \brief set the Lennard-Jones epsilon from text (used by the refinement parameters dialog)

  A running refinement is restarted with the new value.

  @param combobox_item_idx the index of the item in the dialog's combobox
  @param t the value as text */
void set_refinement_lennard_jones_epsilon_from_text(int combobox_item_idx, const char *t);
/*! \brief set the Ramachandran restraints weight from text (used by the refinement parameters dialog)

  A running refinement is restarted with the new value.

  @param combobox_item_idx the index of the item in the dialog's combobox
  @param t the value as text */
void set_refinement_ramachandran_restraints_weight_from_text(int combobox_item_idx, const char *t);
/*! \brief set the overall map weight (c.f. set_matrix()) from text

  A running refinement is restarted with the new value.

  @param t the value as text */
void set_refinement_overall_weight_from_text(const char *t);
/*! \brief set the torsion restraints weight from text (used by the refinement parameters dialog)

  A running refinement is restarted with the new value.

  @param combobox_item_index the index of the item in the dialog's combobox
  @param t the value as text */
void set_refinement_torsion_weight_from_text(int combobox_item_index, const char *t);
/*! \brief record whether the "more control" frame of the refinement
  parameters dialog is visible (1) or not (0) */
void set_refine_params_dialog_more_control_frame_is_active(int state);


/*! \brief return the state of Ramachandran restraints in refinement (1 for on, 0 for off) */
int refine_ramachandran_angles_state();

/*! \brief use numerical gradients in refinement (for debugging)

  @param istate 1 for on, 0 for off (default off) */
void set_numerical_gradients(int istate);

/*! \brief turn on (or off) debugging output for refinement

  @param state 1 for on, 0 for off */
void set_debug_refinement(int state);


/*! \brief correct the sign of chiral volumes before commencing refinement?

   Do we want to fix chiral volumes (by moving the chiral atom to the
   other side of the chiral plane if necessary).  Default yes
   (1). Note: doesn't work currently - the setting is stored but not used. */
void set_fix_chiral_volumes_before_refinement(int istate);

/*! \brief check the chiral volumes of molecule imol

  Shows a dialog listing the atoms with chiral volume errors (and,
  if there were residues with no chiral restraints in the dictionary,
  a dialog listing those residue types).

  @param imol the model molecule index */
void check_chiral_volumes(int imol);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return a list of atom specs of the atoms with chiral volume errors

  @return a list of atom specs (possibly empty), or \#f if imol is not a
  valid model molecule */
SCM chiral_volume_errors_scm(int imol);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return a list of atom specs of the atoms with chiral volume errors

  @return a list of atom specs (possibly empty), or False if imol is not a
  valid model molecule */
PyObject *chiral_volume_errors_py(int imol);
#endif	/* USE_PYTHON */
#endif	/* __cplusplus */


/*! \brief For experienced Cooters who don't like Coot nannying about
  chiral volumes during refinement.

  @param istate 0 to stop the chiral volume errors dialog being shown,
  1 to show it (default 1) */
void set_show_chiral_volume_errors_dialog(short int istate);

/*! \brief set the type of secondary structure restraints

0 no sec str restraints

1 alpha helix restraints

2 beta strand restraints.

Call this before refine_residues_py() to maintain secondary structure
geometry during real-space refinement. Reset to 0 after refinement to
avoid affecting subsequent refinement operations.

*/
void set_secondary_structure_restraints_type(int itype);

/*! \brief return the secondary structure restraints type

  @return 0 (none), 1 (alpha helix) or 2 (beta strand), as for
  set_secondary_structure_restraints_type() */
int secondary_structure_restraints_type();

/*! \brief the molecule number of the map used for refinement

   @return the map number, if it has been set or there is only one
   (non-difference) map, otherwise return -1 on no map set (ambiguous) or no maps.
*/
int imol_refinement_map();	/* return -1 on no map */

/*! \brief set the molecule number of the map to be used for
  refinement/fitting.

  @param imol the map molecule index
  @return imol on success, -1 on failure (e.g. imol is not a valid map)
*/
int set_imol_refinement_map(int imol);	/* returns imol on success, otherwise -1 */

/*! \brief Does the residue exist? (Raw function)

   @return 0 on not-exist, 1 on does exist.
*/
int does_residue_exist_p(int imol, const char *chain_id, int resno, const char *inscode);

/*! \brief delete the restraints for the given comp_id (i.e. residue name)

Only a dictionary that applies to all molecules is deleted (not one
that was read for a specific molecule).

@return success status (0 is failed, 1 is success)
*/
int delete_restraints(const char *comp_id);

/*! \brief add a user-define bond restraint

   this extra restraint is used when the given atoms are selected in
   refinement or regularization.

   @param imol the model molecule index
   @param chain_id_1 the chain id of the first atom
   @param res_no_1 the residue number of the first atom
   @param ins_code_1 the insertion code of the first atom
   @param atom_name_1 the atom name of the first atom
   @param alt_conf_1 the alt conf of the first atom ("" for none)
   @param chain_id_2 the chain id of the second atom
   @param res_no_2 the residue number of the second atom
   @param ins_code_2 the insertion code of the second atom
   @param atom_name_2 the atom name of the second atom
   @param alt_conf_2 the alt conf of the second atom ("" for none)
   @param bond_dist the target distance (in Å)
   @param esd the esd of the target distance (in Å)

   @return the index of the new restraint.

   @return -1 when the atoms were not found and no extra bond
   restraint was stored.  */

int add_extra_bond_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2, double bond_dist, double esd);

/*! \brief add a user-define GM distance restraint

   A Geman-McClure distance restraint is like a bond restraint but
   its penalty is robust (it flattens out) for large deviations.

   this extra restraint is used when the given atoms are selected in
   refinement or regularization.

   @param imol the model molecule index
   @param chain_id_1 the chain id of the first atom
   @param res_no_1 the residue number of the first atom
   @param ins_code_1 the insertion code of the first atom
   @param atom_name_1 the atom name of the first atom
   @param alt_conf_1 the alt conf of the first atom ("" for none)
   @param chain_id_2 the chain id of the second atom
   @param res_no_2 the residue number of the second atom
   @param ins_code_2 the insertion code of the second atom
   @param atom_name_2 the atom name of the second atom
   @param alt_conf_2 the alt conf of the second atom ("" for none)
   @param bond_dist the target distance (in Å)
   @param esd the esd of the target distance (in Å)

   @return the index of the new restraint.

   @return -1 when the atoms were not found and no extra bond
   restraint was stored.  */

int add_extra_geman_mcclure_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2, double bond_dist, double esd);
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief add several extra bond restraints in one go

  @param extra_bond_restraints_scm a list of restraint descriptions,
  each of which is a 4-element list: (atom-spec-1 atom-spec-2
  target-distance esd)
  @return the number of well-formed restraint descriptions (restraints
  whose atoms are not found are not stored but are still counted)
  @param imol the model molecule index
*/
int add_extra_bond_restraints_scm(int imol, SCM extra_bond_restraints_scm);
#endif // USE_GUILE
#ifdef USE_PYTHON
/*! \brief add several extra bond restraints in one go

  @param extra_bond_restraints_py a list of restraint descriptions,
  each of which is a 4-element list: [atom_spec_1, atom_spec_2,
  target_distance, esd]
  @return the number of well-formed restraint descriptions (restraints
  whose atoms are not found are not stored but are still counted)
  @param imol the model molecule index
*/
int add_extra_bond_restraints_py(int imol, PyObject *extra_bond_restraints_py);
#endif // USE_GUILE
#endif

/*! \brief show (or hide) extra distance restraints (in all molecules)

  @param state 1 for on, 0 for off (default on) */
void set_show_extra_distance_restraints(short int state);

/*! \brief add a user-defined angle restraint

  The atom specs of the 3 atoms are given by chain id, residue number,
  insertion code, atom name and alt conf. Note that the restraint is
  stored even if the atoms are not found.

  @param imol the model molecule index
  @param torsion_angle the target angle (in degrees) - despite the parameter name,
  this is an angle, not a torsion
  @param esd the esd of the target angle (in degrees)
  @return the index of the new restraint, or -1 if imol is not a valid
  model molecule
  @param chain_id_1 the chain id of the first atom
  @param res_no_1 the residue number of the first atom
  @param ins_code_1 the insertion code of the first atom
  @param atom_name_1 the atom name of the first atom
  @param alt_conf_1 the alt conf ("" for none) of the first atom
  @param chain_id_2 the chain id of the second atom
  @param res_no_2 the residue number of the second atom
  @param ins_code_2 the insertion code of the second atom
  @param atom_name_2 the atom name of the second atom
  @param alt_conf_2 the alt conf ("" for none) of the second atom
  @param chain_id_3 the chain id of the third atom
  @param res_no_3 the residue number of the third atom
  @param ins_code_3 the insertion code of the third atom
  @param atom_name_3 the atom name of the third atom
  @param alt_conf_3 the alt conf ("" for none) of the third atom
*/
int add_extra_angle_restraint(int imol,
				const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1,
				const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2,
				const char *chain_id_3, int res_no_3, const char *ins_code_3, const char *atom_name_3, const char *alt_conf_3,
				double torsion_angle, double esd);
/*! \brief add a user-defined torsion restraint

  The atom specs of the 4 atoms are given by chain id, residue number,
  insertion code, atom name and alt conf. Note that the restraint is
  stored even if the atoms are not found.

  @param imol the model molecule index
  @param torsion_angle the target torsion angle (in degrees)
  @param esd the esd of the target torsion (in degrees)
  @param period the period of the torsion restraint
  @return the index of the new restraint, or -1 if imol is not a valid
  model molecule
  @param chain_id_1 the chain id of the first atom
  @param res_no_1 the residue number of the first atom
  @param ins_code_1 the insertion code of the first atom
  @param atom_name_1 the atom name of the first atom
  @param alt_conf_1 the alt conf ("" for none) of the first atom
  @param chain_id_2 the chain id of the second atom
  @param res_no_2 the residue number of the second atom
  @param ins_code_2 the insertion code of the second atom
  @param atom_name_2 the atom name of the second atom
  @param alt_conf_2 the alt conf ("" for none) of the second atom
  @param chain_id_3 the chain id of the third atom
  @param res_no_3 the residue number of the third atom
  @param ins_code_3 the insertion code of the third atom
  @param atom_name_3 the atom name of the third atom
  @param alt_conf_3 the alt conf ("" for none) of the third atom
  @param chain_id_4 the chain id of the fourth atom
  @param res_no_4 the residue number of the fourth atom
  @param ins_code_4 the insertion code of the fourth atom
  @param atom_name_4 the atom name of the fourth atom
  @param alt_conf_4 the alt conf ("" for none) of the fourth atom
*/
int add_extra_torsion_restraint(int imol,
				const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1,
				const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2,
				const char *chain_id_3, int res_no_3, const char *ins_code_3, const char *atom_name_3, const char *alt_conf_3,
				const char *chain_id_4, int res_no_4, const char *ins_code_4, const char *atom_name_4, const char *alt_conf_4,
				double torsion_angle, double esd, int period);
/*! \brief add a restraint that keeps the given atom close to its starting position

  If the atom already has a start-position restraint, that restraint is updated.

  @param imol the model molecule index
  @param esd the esd of the restraint (in Å)
  @return the index of the restraint, or -1 if the atom was not found
  @param chain_id_1 the chain id of the first atom
  @param res_no_1 the residue number of the first atom
  @param ins_code_1 the insertion code of the first atom
  @param atom_name_1 the atom name of the first atom
  @param alt_conf_1 the alt conf ("" for none) of the first atom
*/
int add_extra_start_pos_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, double esd);

/*! \brief add a restraint that pulls the given atom to the target position (x, y, z)

  Target position restraints act like atom pull restraints in
  refinement but are stored with the other extra restraints. No
  redraw is done.

  @param imol the model molecule index
  @param chain_id the chain id
  @param res_no the residue number
  @param ins_code the insertion code
  @param atom_name the atom name
  @param alt_conf the alt conf ("" for none)
  @param x the x coordinate of the target position (in Å)
  @param y the y coordinate of the target position (in Å)
  @param z the z coordinate of the target position (in Å)
  @param weight the weight of the restraint
  @return 1 on success, -1 if the atom was not found (or imol is not a
  valid model molecule) */
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
/*! \brief clear out all the extra/user-defined restraints for the given residue spec in molecule number imol */
void delete_extra_restraints_for_residue_spec_scm(int imol, SCM residue_spec_in);
#endif // USE_GUILE
#ifdef USE_PYTHON
/*! \brief clear out all the extra/user-defined restraints for the given residue spec in molecule number imol */
void delete_extra_restraints_for_residue_spec_py(int imol, PyObject *residue_spec_in_py);
#endif // USE_PYTHON
#endif // __cplusplus

/*! \brief delete extra bond restraints that are badly violated

  Extra bond restraints for which |target distance - current distance|/esd
  is n_sigma or more are deleted (other types of extra restraints are not
  affected).

  @param imol the model molecule index
  @param n_sigma the deviation (in units of the restraint esd) at or beyond
  which restraints are deleted */
void delete_extra_restraints_worse_than(int imol, float n_sigma);

/*! \brief read in prosmart (typically) extra restraints

  The restraints (in Refmac external-restraints format) are added to the
  extra restraints of molecule imol. */
void add_refmac_extra_restraints(int imol, const char *file_name);

/*! \brief show (or hide) the extra restraints of molecule imol

  @param state 1 for on, 0 for off
  @param imol the model molecule index
*/
void set_show_extra_restraints(int imol, int state);
/*! \brief are the extra restraints of molecule imol being shown?

  @return 1 for yes, 0 for no (or if imol is not a valid model molecule) */
int extra_restraints_are_shown(int imol);

/*! \brief often we don't want to see all prosmart restraints, just the (big) violations

  An extra bond restraint is drawn only if its signed deviation
  (current distance - target distance)/esd is less than or equal to
  the lower limit or greater than or equal to the upper limit, so for
  example (-2.0, 2.0) shows only restraints deviating by 2 sigma or more.

  Note: despite the parameter names in this declaration, the second
  argument (\c limit_high here) is used as the lower limit and the
  third argument (\c limit_low here) as the upper limit, i.e. pass
  the lower limit first. */
void set_extra_restraints_prosmart_sigma_limits(int imol, double limit_high, double limit_low);

/*! \brief generate external distance local self restraints

  Geman-McClure distance restraints (esd 0.05 Å) are generated between
  all non-hydrogen atom pairs in the given chain that are within
  local_dist_max of each other (excluding bonded and angle-related
  pairs), with the current distances as targets. Existing extra bond
  restraints are cleared.

  @param imol the model molecule index
  @param chain_id the chain id
  @param local_dist_max the maximum distance (in Å) */
void generate_local_self_restraints(int imol, const char *chain_id, float local_dist_max);

/*! \brief generate external distance all-molecule self restraints

  As generate_local_self_restraints(), but for all atoms of molecule imol.

  @param imol the model molecule index
  @param local_dist_max the maximum distance (in Å) */
void generate_self_restraints(int imol, float local_dist_max);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief generate external distance self restraints for selected residues

  As generate_local_self_restraints(), but for the atoms of the residues
  in residue_specs (a list of residue specs). */
void generate_local_self_restraints_by_residues_scm(int imol, SCM residue_specs, float local_dist_max);
#endif // USE_GUILE
#ifdef USE_PYTHON
/*! \brief generate external distance self restraints for selected residues

  As generate_local_self_restraints(), but for the atoms of the residues
  in residue_specs (a list of residue specs).

  @param imol the model molecule index
  @param residue_specs a list of residue specs
  @param local_dist_max the maximum distance (in Å) */
void generate_local_self_restraints_by_residues_py(int imol, PyObject *residue_specs, float local_dist_max);
#endif // USE_PYTHON
#endif // __cplusplus


/*! \brief proSMART interpolated restraints for model morphing

  Write restraint files interpolated between the extra restraints of
  molecules imol_1 and imol_2. n_steps must be greater than 2 and less
  than 5000.

  @param imol_1 the model molecule with the starting extra restraints
  @param imol_2 the model molecule with the final extra restraints
  @param n_steps the number of interpolation steps (and restraint files)
  @param file_name_stub the stub for the output file names */
void write_interpolated_extra_restraints(int imol_1, int imol_2, int n_steps, const char *file_name_stub);

/*! \brief proSMART interpolated restraints for model morphing and write interpolated model

n_steps must be greater than 2 and less than 5000.

interpolation_mode is currently dummy - in due course I will addd torion angle interpolation.
*/
void write_interpolated_models_and_extra_restraints(int imol_1, int imol_2, int n_steps, const char *file_name_stub,
						    int interpolation_mode);

/*! \brief show (or hide) the parallel plane restraints of molecule imol

  @param state 1 for on, 0 for off
  @param imol the model molecule index
*/
void set_show_parallel_plane_restraints(int imol, int state);
/*! \brief are the parallel plane restraints of molecule imol being shown?

  @return 1 for yes, 0 for no (or if imol is not a valid model molecule) */
int parallel_plane_restraints_are_shown(int imol);
/*! \brief add a parallel plane restraint between 2 residues

  For nucleotides the base atoms are used; for PHE, TYR, TRP and ARG
  the planar side-chain atoms are used.

  @param imol the model molecule index
  @param chain_id_1 the chain id of the first residue
  @param re_no_1 the residue number of the first residue
  @param ins_code_1 the insertion code of the first residue
  @param chain_id_2 the chain id of the second residue
  @param re_no_2 the residue number of the second residue
  @param ins_code_2 the insertion code of the second residue
*/
void add_parallel_plane_restraint(int imol,
				  const char *chain_id_1, int re_no_1, const char *ins_code_1,
				  const char *chain_id_2, int re_no_2, const char *ins_code_2);
/*! \brief draw the extra bond restraints of molecule imol between CA atoms

  @param state 1 to draw restraints between different residues as dashed
  lines between the residues' CA atoms, 0 to draw them between the
  restrained atoms
  @param imol the model molecule index
*/
void set_extra_restraints_representation_for_bonds_go_to_CA(int imol, short int state);


#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief delete an extra restraint

  restraint_spec is something like (list 'bond spec-1 spec-2); 'angle
  (3 specs), 'torsion (4 specs) and 'start-pos (1 spec) are also
  understood.

  spec-1 and spec-2 do not have to be in the order that the bond was created.  */
void delete_extra_restraint_scm(int imol, SCM restraint_spec);
/*! \brief list the extra restraints of molecule imol

  @return a list of restraints, e.g. ('bond spec-1 spec-2 dist esd),
  ('angle spec-1 spec-2 spec-3 angle esd), ('torsion spec-1 spec-2
  spec-3 spec-4 torsion esd period) and ('start-pos spec esd), or \#f
  if there are none */
SCM list_extra_restraints_scm(int imol);
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief delete an extra restraint

  restraint_spec is something like ['bond', spec_1, spec_2]; "angle"
  (3 specs), "torsion" (4 specs) and "start pos" (1 spec) are also
  understood.

  spec_1 and spec_2 do not have to be in the order that the bond was created.  */
void delete_extra_restraint_py(int imol, PyObject *restraint_spec);
/*! \brief list the extra restraints of molecule imol

  @return a list of restraints, e.g. ["bond", spec_1, spec_2, dist, esd],
  ["angle", spec_1, spec_2, spec_3, angle, esd], ["torsion", spec_1,
  spec_2, spec_3, spec_4, torsion, esd, period] and ["start pos",
  spec, esd], or False if there are none */
PyObject *list_extra_restraints_py(int imol);
#endif /* USE_PYTHON */
#endif /*  __cplusplus */

/*! \brief set use only extra torsion restraints for torsions

  @param state 1 for on, 0 for off (default off) */
void set_use_only_extra_torsion_restraints_for_torsions(short int state);
/*! \brief return only-use-extra-torsion-restraints-for-torsions state */
int use_only_extra_torsion_restraints_for_torsions_state();

/*! \brief clear all atom pull restraints (and redraw) */
void clear_all_atom_pull_restraints();

/*! \brief set auto-clear atom pull restraint

  @param state 1 for on, 0 for off (default on) */
void set_auto_clear_atom_pull_restraint(int state);

/*! \brief get auto-clear atom pull restraint state */
int  get_auto_clear_atom_pull_restraint_state();

/*! \brief increase the proportional editing radius*/
void increase_proportional_editing_radius();

/*! \brief decrease the proportional editing radius*/
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

/*! \brief refine residue range using simplex optimization

  The atoms of residues res1 to res2 in chain chain_id (with the given
  alt loc) of molecule imol are fitted to the map imol_for_map as a
  rigid body by simplex optimization. A backup is made first.

  @param res1 the first residue number
  @param res2 the last residue number
  @param altloc the alt loc of the atoms to fit
  @param chain_id the chain id
  @param imol the model molecule index
  @param imol_for_map the map molecule index */
void
fit_residue_range_to_map_by_simplex(int res1, int res2, const char *altloc, const char *chain_id, int imol, int imol_for_map);

/*! \brief simply score the residue range fit to map

  The score is the sum over the selected atoms of the map density at
  the atom position multiplied by the atom occupancy.

  @param res1 the first residue number
  @param res2 the last residue number
  @param altloc the alt loc of the atoms to score
  @param chain_id the chain id
  @param imol the model molecule index
  @param imol_for_map the map molecule index
  @return the score, or 0 if no atoms were selected or the molecules
  are not valid */
float
score_residue_range_fit_to_map(int res1, int res2, const char *altloc, const char *chain_id, int imol, int imol_for_map);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*               Nomenclature Errors                                        */
/*  ----------------------------------------------------------------------- */
/*! \name Nomenclature Errors */
/*! \{ */
/*! \brief fix nomenclature errors in molecule number imol

   Atom names in the side chains of PHE, TYR, ASP, GLU, LEU and VAL
   residues are swapped where they do not follow the naming
   convention. A backup is made first.

   @return the number of residues altered. */
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

/*! \brief output to the terminal the Atom Info for the given atom specs

  The atom name, model, chain, residue number and name, occupancy,
  B-factor, element and position are printed, followed by the side-chain
  chi angles of the residue.
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
/*! \brief start residue info picking: the user is asked to click on an atom

  If there are pending (un-applied) residue info edits, a warning dialog
  is shown instead. */
void do_residue_info_dialog();

/* MOVE-ME to c-interface-gtk-widgets.h */
/*! \brief show the residue info dialog for the residue of the given atom

  @param imol the model molecule index
  @param atom_index the index of an atom in the molecule's atom selection */
 void output_residue_info_dialog    (int imol, int atom_index); /* widget version */
/* scripting version */
/*! \brief show residue info dialog for given residue */
void residue_info_dialog(int imol, const char *chain_id, int resno, const char *ins_code);
/*! \brief is the residue info dialog displayed?

  @return 1 for yes, 0 for no */
int residue_info_dialog_is_displayed();
/*! \brief output the residue info for the residue of the given atom as text

  Note the argument order: atom_index first, then imol. */
void output_residue_info_as_text(int atom_index, int imol); /* text version */
/* functions that uses mmdb_manager functions/data types moved to graphics_info_t */

/*! \brief show the distance labels

  @param state 1 for on, 0 for off (default on)
 * */
void set_show_distance_labels(short int state);

/*! \brief start a distance measurement: the user is asked to click on 2 atoms */
void do_distance_define();
/*! \brief start an angle measurement: the user is asked to click on 3 atoms */
void do_angle_define();
/*! \brief start a torsion measurement: the user is asked to click on 4 atoms */
void do_torsion_define();
/*! \brief callback for the residue info "apply to all" checkbutton - currently does nothing */
void residue_info_apply_all_checkbutton_toggled();
/*! \brief clear the list of pending residue info edits */
void clear_residue_info_edit_list();

/* a graphics_info_t function wrapper: */
/*! \brief forget the residue info dialog (i.e. note that it is no longer displayed) */
void unset_residue_info_widget();
/*! \brief clear all the distance, angle and torsion measurements */
void clear_measure_distances();
/*! \brief clear the last distance measurement */
void clear_last_measure_distance();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*               GUI edit functions                                   */
/*  ----------------------------------------------------------------------- */
/* section Edit Fuctions */
/*! \name Edit Fuctions */
/*! \{ */

/*! \brief show the GUI for copying a molecule */
void  do_edit_copy_molecule();
/*! \brief show the GUI for copying a fragment (an atom selection) of a molecule to a new molecule */
void  do_edit_copy_fragment();
/*! \brief show the GUI for replacing a fragment of a molecule with the
  corresponding atoms from a reference molecule */
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
/*! \brief show H-bond environment distances.  If state is 0, H-bond distances are
  turned off, otherwise H-bond distances are turned on. */
void set_show_environment_distances_h_bonds(int state);
/*! \brief show the state of display of the  environment distances

  @return 1 if environment distances are shown, 0 if not */
int show_environment_distances_state();
/*! \brief min and max distances for the environment distances

  @param min_dist the minimum distance (in Å, default 0.0)
  @param max_dist the maximum distance (in Å, default 3.2) */
void set_environment_distances_distance_limits(float min_dist, float max_dist);

/*! \brief show the environment distances with solid modelling

  @param state 1 for on, 0 for off */
void set_show_environment_distances_as_solid(int state);

/*! \brief Label the atom on Environment Distances start/change

  @param state 1 for on, 0 for off (default off) */
void set_environment_distances_label_atom(int state);

/*! \brief Label the atoms in the residues around the central residue

  The central residue is the residue of the active atom; for each
  neighbouring residue within 4 Å, the closest atoms are labelled. */
void label_neighbours();

/*! \brief Label the atoms in the central residue

  The central residue is the residue of the active atom. */
void label_atoms_in_residue();

/*! \brief Label the atoms with their B-factors

  The non-hydrogen atoms within 8 Å of the screen centre in the
  molecule of the active atom are labelled.

  @param state 1 for on, 0 for off */
void set_show_local_b_factors(short int state);

/*! \brief Add a geometry distance between points in a given molecule

The distance is added to the distance measurements in the graphics.
imol_1 and imol_2 are not used.

@return the distance between the points

*/
double add_geometry_distance(int imol_1, float x_1, float y_1, float z_1, int imol_2, float x_2, float y_2, float z_2);
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief Add a geometry distance between 2 atoms (given by atom specs)

@return the distance between the atoms (in Å), or -1 if either atom was
not found or either molecule is not a valid model molecule */
double add_atom_geometry_distance_scm(int imol_1, SCM atom_spec_1, int imol_2, SCM atom_spec_2);
#endif
#ifdef USE_PYTHON
/*! \brief Add a geometry distance between 2 atoms (given by atom specs)

@return the distance between the atoms (in Å), or -1 if either atom was
not found or either molecule is not a valid model molecule */
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

the origin is top-left: the values are the last known pointer position
in the graphics window divided by the window width and height.
may return false if pointer is not available (e.g. there is no graphics interface) */
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
/*! \brief show the state of display of the  pointer distances

  @return 1 for on, 0 for off */
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
/*! \brief scale the view by f without a redraw (internal interface)

  Values outside the range 0.5 to 1.8 have no effect. */
void scale_zoom_internal(float f);
/*! \brief return the current zoom factor i.e. get_zoom_factor()

  The default is 100. */
float zoom_factor();

/*! \brief set smooth scroll with zoom
   @param i 0 means no, 1 means yes: (default 0) */
void set_smooth_scroll_do_zoom(int i);
/* default 0 (off) */
/*! \brief return the state of the above system */
int      smooth_scroll_do_zoom();
/*! \brief return the smooth scroll zoom limit (default 30.0) */
float    smooth_scroll_zoom_limit();
/*! \brief set the smooth scroll zoom limit (default 30.0) */
void set_smooth_scroll_zoom_limit(float f);

/*! \brief set the zoom factor (absolute value) - maybe should be called set_zoom_factor()

  A redraw is done. The default zoom factor is 100. */
void set_zoom(float f);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  CNS data stuff                                          */
/*  ----------------------------------------------------------------------- */
/*! \name CNS Data Functions */
/*! \{ */
/*! \brief read CNS data and make a map

  The CNS reflection file (F and phi) is read and a map is calculated
  in a new molecule, using the cell and space group of the model
  molecule imol.

  @param filename the CNS reflection file name
  @param imol a model molecule that provides the cell and space group
  @return the new map molecule index, or -1 on failure */
int handle_cns_data_file(const char *filename, int imol);

/*! \brief read CNS data and make a map, using the given cell and space group

a, b,c are in Angstroems.  alpha, beta, gamma are in degrees.  spg is
the space group info, either ;-delimited symmetry operators or the
space group name. imol is not used.

Note: the gamma argument is currently ignored (alpha is used in its place).

@return the new map molecule index */
int handle_cns_data_file_with_cell(const char *filename, int imol, float a, float b, float c, float alpha, float beta, float gamma, const char *spg_info);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  cif stuff                                               */
/*  ----------------------------------------------------------------------- */
/* section mmCIF Functions */
/*! \name mmCIF Functions */
/* dataset stuff */
/*! \{ */
/*! \brief read a reflection mmCIF file with phases and make maps

  Makes a sigmaA-weighted map and a difference sigmaA map (both in new
  molecules) from the Fobs and Fcalc/phases in the file.

  @return the molecule index of the sigmaA map, or -1 on failure */
int auto_read_cif_data_with_phases(const char *filename);
/*! \brief read a reflection mmCIF file with phases and make a sigmaA-weighted map

  @return the new map molecule index, or -1 on failure */
int read_cif_data_with_phases_sigmaa(const char *filename);
/*! \brief read a reflection mmCIF file with phases and make a sigmaA-weighted difference map

  @return the new map molecule index, or -1 on failure */
int read_cif_data_with_phases_diff_sigmaa(const char *filename);
/*! \brief read a reflection mmCIF file (Fobs) and make a sigmaA map

  Structure factors are calculated from the model molecule imol_coords.

  @return the new map molecule index, or -1 on failure */
int read_cif_data(const char *filename, int imol_coords);
/*! \brief read a reflection mmCIF file (Fobs) and make a 2Fo-Fc map

  Structure factors are calculated from the model molecule imol_coords.
  Note: currently this makes the same (sigmaA) map as read_cif_data().

  @return the new map molecule index, or -1 on failure */
int read_cif_data_2fofc_map(const char *filename, int imol_coords);
/*! \brief read a reflection mmCIF file (Fobs) and make an Fo-Fc map

  Structure factors are calculated from the model molecule imol_coords.

  @return the new map molecule index, or -1 on failure */
int read_cif_data_fofc_map(const char *filename, int imol_coords);
/*! \brief read a reflection mmCIF file with phases (Fobs and Fcalc, phicalc)
  and make an Fo-Fc map

  @return the new map molecule index, or -1 on failure */
int read_cif_data_with_phases_fo_fc(const char *filename);
/*! \brief read a reflection mmCIF file with phases (Fobs and Fcalc, phicalc)
  and make a 2Fo-Fc map

  @return the new map molecule index, or -1 on failure */
int read_cif_data_with_phases_2fo_fc(const char *filename);
/*! \brief read a reflection mmCIF file with phases (Fobs and Fcalc, phicalc)
  and make a map of the given type

  @param filename the reflection mmCIF file name
  @param map_type 1 for 2Fo-Fc, 2 for Fo-Fc, 3 for Fo with calculated phases
  @return the new map molecule index, or -1 on failure */
int read_cif_data_with_phases_nfo_fc(const char *filename,
				     int map_type);
/*! \brief read a reflection mmCIF file with phases and make a map using Fo
  and the calculated phases

  @return the new map molecule index, or -1 on failure */
int read_cif_data_with_phases_fo_alpha_calc(const char *filename);

/*! \brief write a connectivity file for the given monomer

  The file (a "RESIDUE"/"CONECT" style listing of the bonded atoms of
  each atom) is made from the bonds in the monomer's dictionary.

  @param monomer_name the residue type (comp_id)
  @param filename the output file name
  @return 1 on success, 0 on failure (e.g. no dictionary for the monomer) */
int write_connectivity(const char* monomer_name, const char *filename);
/*! \brief open the cif dictionary file selector dialog (if the graphics interface is in use) */
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

/*! \brief import all the monomer dictionaries in the Refmac monomer library

   Use the environment variable COOT_REFMAC_LIB_DIR to find cif files
   in subdirectories (of data/monomers) and import them all. */
void import_all_refmac_cifs();

/*! \brief read a small-molecule cif file and make a model molecule

   Symmetry display is turned on.

   @return the new molecule index, or -1 on failure */
int read_small_molecule_cif(const char *file_name);

/*! \brief read the reflection data from a small-molecule cif file and make maps

   If maps can be calculated from the data, a sigmaA map and a
   difference sigmaA map are made in new molecules.

   @return the index of the new (first) molecule, or -1 if the data
   could not be read */
int read_small_molecule_data_cif(const char *file_name);

/*! \brief read a SHELX-style small-molecule cif (e.g. from the COD) and make maps

   The model is read from the embedded SHELX res file and, from the embedded
   SHELX hkl reflection list, a (model-phased) 2Fo-Fc map and an Fo-Fc
   difference map are made.  Return the model molecule index (-1 on failure). */
int read_small_molecule_cif_and_make_map(const char *file_name);

/*! \brief read the reflection data from a small-molecule cif and make maps
  phased by the model imol_coords

  Structure factors are calculated from the model, and a 2Fo-Fc map and
  an Fo-Fc difference map are made in new molecules.

  @return the molecule index of the 2Fo-Fc map, or -1 on failure */
int read_small_molecule_data_cif_and_make_map_using_coords(const char *file_name,
							   int imol_coords);

/*! \} */
/*  ------------------------------------------------------------------------ */
/*                         Validation:                                       */
/*  ------------------------------------------------------------------------ */
/* section Validation Functions */
/*! \name Validation Functions */
/*! \{ */
/*! \brief (legacy) deviant geometry analysis for molecule number imol

  Currently this has no visible effect: it builds restraints for each
  chain of the molecule but reports nothing.

  @param imol the model molecule index */
void deviant_geometry(int imol);
/*! \brief is imol a valid model molecule?

  @param imol the molecule index
  @return 1 if imol is a molecule index that has a model (coordinates), 0 otherwise */
short int is_valid_model_molecule(int imol);
/*! \brief is imol a valid map molecule?

  @param imol the molecule index
  @return 1 if imol is a molecule index that has a map, 0 otherwise */
short int is_valid_map_molecule(int imol);


/*! \brief generate a list of difference map peaks

Peaks are found in the difference map imol (nothing happens if imol is
not a difference map). Peaks within max_closeness (2.0 A typically) of
a larger peak are not listed. In graphical mode the peaks are shown in
the difference map peaks panel (or a "No difference map peaks" dialog
is shown); the peaks are also written to the log.

If imol_coords is a valid model molecule, peaks are moved (by symmetry)
to be close to the model.

Note: around_model_only_flag is currently ignored - the underlying peak
search forces it off.

@param imol the difference map molecule index
@param imol_coords the model molecule index (may be invalid, e.g. -1)
@param level the peak search level in map rmsd (sigma) units
@param max_closeness the minimum separation (in A) between listed peaks
@param do_positive_level_flag 1 to find positive peaks
@param do_negative_level_flag 1 to find negative peaks
@param around_model_only_flag (currently ignored) intended to limit peaks
       to those within 4A of the model
*/
void difference_map_peaks(int imol, int imol_coords, float level, float max_closeness, int do_positive_level_flag, int do_negative_level_flag, int around_model_only_flag);

/*! \brief set the max closeness (i.e. no smaller peaks can be within
   max_closeness of a larger peak)

In the GUI for difference map peaks, there is not a means to set the
max_closeness, so here is a means to set it and query it.

@param m the max closeness in A (default 2.0) */
void set_difference_map_peaks_max_closeness(float m);
/*! \brief return the max closeness (in A) used for difference map peaks (default 2.0) */
float difference_map_peaks_max_closeness();

/*! \brief clear the difference map peaks (e.g. the peak markers in the graphics) */
void clear_diff_map_peaks();

/*! \brief find GLN and ASN B-factor outliers

  Compares the B-factor difference of the side-chain O and N atoms of
  GLN and ASN residues to the distribution of B-factors of the other
  atoms (a Z score). Only works in graphical mode: the outliers are
  printed and listed in an "interesting things" dialog with buttons to
  flip the side chain (or a "no outliers" dialog is shown).

  @param imol the model molecule index */
void gln_asn_b_factor_outliers(int imol);
#ifdef USE_PYTHON
/*! \brief old Python variant of gln_asn_b_factor_outliers()

  The outliers are printed; the dialog is only made when Coot is
  compiled with PyGTK, otherwise use gln_asn_b_factor_outliers(). */
void gln_asn_b_factor_outliers_py(int imol);
#endif /*  USE_PYTHON */

#ifdef __cplusplus
#ifdef USE_PYTHON
/*! \brief return a list of map peaks of molecule number imol_map
  above n_sigma.

  Only positive peaks are found. Clusters of grid points above the
  level are reduced to one peak, but there is no further distance
  filtering of the peaks.

  @param imol_map the map molecule index
  @param n_sigma the level in map rmsd (sigma) units
  @return a list of [x, y, z] cartesian coordinates or Python False if
  imol_map is not a valid map molecule. */
PyObject *map_peaks_py(int imol_map, float n_sigma);
/*! \brief return a list of positive map peaks near a point

  Peaks above n_sigma are moved by symmetry to be close to the given
  point and those within radius of the point are returned.

  @param imol_map the map molecule index
  @param n_sigma the level in map rmsd (sigma) units
  @param x the x coordinate of the point
  @param y the y coordinate of the point
  @param z the z coordinate of the point
  @param radius the search radius (A)
  @return a list of [x, y, z, density_value] items (the density value is
  in map units, not sigma) or Python False if imol_map is not a valid map. */
PyObject *map_peaks_near_point_py(int imol_map, float n_sigma, float x, float y, float z, float radius);
/*! \brief filter a given list of peaks to those near a point

  Each peak (an [x, y, z] list) is tested over the symmetry operators
  and nearby cell translations of the map; peaks with a copy within
  radius of the point are returned (as that copy). No map peak search
  is done.

  @param imol_map the map molecule index (provides the cell and symmetry)
  @param peak_list a list of [x, y, z] peak positions
  @param x the x coordinate of the point
  @param y the y coordinate of the point
  @param z the z coordinate of the point
  @param radius the search radius (A)
  @return a list of [x, y, z] positions or Python False if imol_map is not
  a valid map. */
PyObject *map_peaks_near_point_from_list_py(int imol_map, PyObject *peak_list, float x, float y, float z, float radius);
/*! \brief return a list of map peaks around a model molecule

  Peaks are moved by symmetry to be close to the model and are filtered
  by difference_map_peaks_max_closeness().

  @param imol_map the map molecule index
  @param sigma the level in map rmsd (sigma) units
  @param negative_also_flag 1 to also find negative peaks
  @param imol_coords the model molecule index
  @return a list of [density_value, [x, y, z]] items or Python False if
  either molecule index is not valid. */
PyObject *map_peaks_around_molecule_py(int imol_map, float sigma, int negative_also_flag, int imol_coords);

/* BL says:: this probably shouldnt be here but cluster with KK code */
/*! \brief return the screen axes in world coordinates

  @return a list of 3 vectors [screen_x, screen_y, screen_z], each a list
  of 3 numbers, or Python None if there is no graphics. */
PyObject *screen_vectors_py();
#endif /*  USE_PYTHON */

#ifdef USE_GUILE
/*! \brief return a list of map peaks of molecule number imol_map
  above n_sigma.  Only positive peaks; clusters of grid points are
  reduced to one peak.
  Return a list of 3d cartestian coordinates or scheme false if
  imol_map is not a valid map molecule. */
SCM map_peaks_scm(int imol_map, float n_sigma);
/*! \brief return a list of positive map peaks above n_sigma within
  radius (A) of the point (x, y, z): a list of (x y z density-value)
  items, or scheme false if imol_map is not a valid map molecule. */
SCM map_peaks_near_point_scm(int imol_map, float n_sigma, float x, float y, float z, float radius);
#endif  /* USE_GUILE */


/* does this live here really? */
#ifdef USE_GUILE
/*! \brief return the torsion angle (in degrees) defined by the 4 given
  atom specs, or scheme false if the molecule or atoms are not found */
SCM get_torsion_scm(int imol, SCM atom_spec_1, SCM atom_spec_2, SCM atom_spec_3, SCM atom_spec_4);

/*! \brief set the given torsion the given residue. tors is in
  degrees.  Return the resulting torsion (also in degrees).

  The torsion is set using the dictionary of the residue type (a
  backup is made). Returns scheme false for an invalid molecule and
  -999.9 if the residue or its dictionary is not found. */
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
/*! \brief return the torsion angle defined by 4 atoms

  @param imol the model molecule index
  @param atom_spec_1 the first atom spec
  @param atom_spec_2 the second atom spec
  @param atom_spec_3 the third atom spec
  @param atom_spec_4 the fourth atom spec
  @return the torsion in degrees, or Python False if imol is not a
  valid model molecule or (some of) the atoms are not found. */
PyObject *get_torsion_py(int imol, PyObject *atom_spec_1, PyObject *atom_spec_2, PyObject *atom_spec_3, PyObject *atom_spec_4);

/*! \brief set the given torsion the given residue. tors is in
  degrees.  Return the resulting torsion (also in degrees).

  The torsion is set using the dictionary of the residue type (so the
  dictionary must be available); a backup is made.

  @return the resulting torsion, -999.9 if the residue or its
  dictionary is not found, or Python False if imol is not a valid
  model molecule. */
PyObject *set_torsion_py(int imol, const char *chain_id, int res_no, const char *insertion_code,
		         const char *alt_conf,
		         const char *atom_name_1,
		         const char *atom_name_2,
		         const char *atom_name_3,
		         const char *atom_name_4, double tors);

/*! \brief create a multi-residue torsion dialog (user manipulation of torsions)

  @param imol the model molecule index
  @param residues_specs_py a list of residue specs */
void multi_residue_torsion_py(int imol, PyObject *residues_specs_py);

#endif  /* USE_PYTHON */


#endif /* __cplusplus  */

/* These functions are called from callbacks.c */
/*! \brief leave multi-residue torsion mode */
void clear_multi_residue_torsion_mode();
/*! \brief set the multi-residue torsion reverse-fragment mode (1 for on, 0 for off) */
void set_multi_residue_torsion_reverse_mode(short int mode);
/*! \brief show the rotatable bonds dialog for the picked residues */
void show_multi_residue_torsion_dialog(); /* show the rotatable bonds dialog */
/*! \brief start picking residues for multi-residue torsion manipulation (shows the pick dialog) */
void setup_multi_residue_torsion();  /* show the pick dialog */

/*! \brief return the atom overlap score

  The score is the mean overlap volume of the atom overlaps of the
  molecule, multiplied by 1000 (0 if there are no overlaps).

  @param imol the model molecule index
  @return the score, or -1 if imol is not a valid model molecule */
float atom_overlap_score(int imol);


/*! \brief set the state of showing chiral volume outlier markers - of a model molecule that is,
   not the intermediate atoms (derived from restraints)

  @param imol the model molecule index
  @param state 1 for on, 0 for off */
void set_show_chiral_volume_outliers(int imol, int state);

/* Note: this is a duplicate declaration of set_show_chiral_volume_outliers().
   (The function for non-bonded contact markers is
   set_show_non_bonded_contact_baddies_markers(imol, state)). */
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

/*! \brief Ramachandran plot for molecule number imol

  Note: this is the old (goocanvas) plot, which is not compiled in
  current builds, so this function does nothing. */
void do_ramachandran_plot(int imol);

/*! \brief set the number of biggest difference arrows on the Kleywegt
  plot (default 50).  */
void set_kleywegt_plot_n_diffs(int n_diffs);

/*! \brief set the contour levels for the ramachandran plot, default
  values are 0.02 (prefered) 0.002 (allowed)

  These levels are also used as the preferred/allowed thresholds when
  scoring the Ramachandran plot of a molecule. */
void set_ramachandran_plot_contour_levels(float level_prefered, float level_allowed);
/*! \brief set the ramachandran plot background block size.

  Smaller is smoother but slower.  Should be divisible exactly into
  360.  Default value is 2. */
void set_ramachandran_plot_background_block_size(float blocksize) ;

/*! \brief set the psi axis for the ramachandran plot. Default (0) from -180
 to 180. Alternative (1) from -120 to 240.

 Note: this currently does nothing. */
void set_ramachandran_psi_axis_mode(int mode);
/*! \brief return the psi axis mode of the ramachandran plot (currently always 0) */
int ramachandran_psi_axis_mode();

/*! \brief set the phi/psi of the residue being edited (moving atoms)

  Note: this is part of the old phi/psi editing tool and currently does nothing. */
void set_moving_atoms(double phi, double psi);

/*! \brief this does the same as `accept_moving_atoms()`

  (and then clears the moving atoms object)
*/
void accept_phi_psi_moving_atoms();

/*! \brief set the pick mode for phi/psi editing

  @param state 1 to start picking (click on an atom in the residue), 0 to stop */
void setup_edit_phi_psi(short int state);	/* a button callback */

/* no need to export this to scripting interface */
void setup_dynamic_distances(short int state);

/*! \brief destroy the edit-backbone Ramachandran plot */
void destroy_edit_backbone_rama_plot();

/*! \brief  2 molecule ramachandran plot (NCS differences) a.k.a. A Kleywegt Plot.

  Note: not yet converted to the current GUI - this function currently only
  prints an error. */
void ramachandran_plot_differences(int imol1, int imol2);

/*! \brief  A chain-specific Kleywegt Plot.

  Note: not yet converted to the current GUI - this function currently does nothing. */
void ramachandran_plot_differences_by_chain(int imol1, int imol2,
					    const char *a_chain, const char *b_chain);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*           sequence_view                                                  */
/*  ----------------------------------------------------------------------- */
/*! \name Sequence View Interface  */
/*! \{ */
/*! \brief display the sequence view dialog for molecule number imol

  The sequence view is added to the sequence view pane of the main
  window. Clicking on a residue centres the view on that residue. */
void sequence_view(int imol);

/*! \brief old name for the above function */
void do_sequence_view(int imol);

/*!  \brief update the sequnce view current position highlight based on active atom */
void update_sequence_view_current_position_highlight_from_active_atom();

/*! \brief remove the sequence view for molecule number imol from the
  sequence view pane (the pane is hidden if it becomes empty) */
void remove_sequence_view_from_sequence_view_box(int imol);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*           rotate moving atoms peptide                                    */
/*  ----------------------------------------------------------------------- */
/*! \brief (old backbone edit) rotate the carbonyl of the moving atoms by angle degrees

  Note: this currently does nothing (the old canvas-based tool is not compiled). */
void change_peptide_carbonyl_by(double angle);/*  in degrees. */
/*! \brief (old backbone edit) rotate the peptide of the moving atoms by angle degrees

  Note: this currently does nothing (the old canvas-based tool is not compiled). */
void change_peptide_peptide_by(double angle);  /* in degress */
/*! \brief make the moving atoms for backbone torsion editing around the
  peptide containing the given atom (atom_index in imol) */
void execute_setup_backbone_torsion_edit(int imol, int atom_index);
/*! \brief set the pick mode for backbone torsion editing

  @param state 1 to start picking (click on an atom in the peptide), 0 to stop.
  Not available while there are moving atoms. */
void setup_backbone_torsion_edit(short int state);

/*! \brief (old backbone edit) record the mouse start position for the peptide drag */
void set_backbone_torsion_peptide_button_start_pos(int ix, int iy);
/*! \brief (old backbone edit) rotate the peptide by 0.05 degrees per pixel of
  x mouse movement since the start position */
void change_peptide_peptide_by_current_button_pos(int ix, int iy);
/*! \brief (old backbone edit) record the mouse start position for the carbonyl drag */
void set_backbone_torsion_carbonyl_button_start_pos(int ix, int iy);
/*! \brief (old backbone edit) rotate the carbonyl by 0.05 degrees per pixel of
  x mouse movement since the start position */
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
/*! \brief add a label to the given atom

  @param imol the model molecule index
  @param chain_id the chain id
  @param iresno the residue number
  @param atom_id the atom name, e.g. " CA "
  @return the atom index of the labelled atom, -1 if the atom was not
  found (or 0 if imol is not a valid model molecule) */
int    add_atom_label(int imol, const char *chain_id, int iresno, const char *atom_id);
/*! \brief remove the label from the given atom

  @param imol the model molecule index (must be valid - it is not checked)
  @param chain_id the chain id
  @param iresno the residue number
  @param atom_id the atom name, padded as in the PDB file, e.g. " CA "
  @return the atom index, or -1 if the atom was not found */
int remove_atom_label(int imol, const char *chain_id, int iresno, const char *atom_id);
/*! \brief remove all atom labels in all molecules */
void remove_all_atom_labels();

/*! \brief label the central atom when recentring (e.g. on go-to-atom)?

  @param i 1 for on (default), 0 for off */
void set_label_on_recentre_flag(int i); /* 0 for off, 1 or on */

/*! \brief return the label-on-recentre state (1 for on, 0 for off) */
int centre_atom_label_status();

/*! \brief use brief atom names for on-screen labels

 call with istat=1 to use brief labels, istat=0 for normal labels (default) */
void set_brief_atom_labels(int istat);

/*! \brief the brief atom label state */
int brief_atom_labels_state();

/*! \brief set if brief atom labels should have seg-ids also

 @param istat 1 for on, 0 for off (default) */
void set_seg_ids_in_atom_labels(int istat);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  scene rotation                                          */
/*  ----------------------------------------------------------------------- */
/* section Screen Rotation */
/*! \name Screen Rotation */
/* stepsize in degrees */
/*! \{ */
/*! \brief rotate view round y axis stepsize degrees for nstep such steps

  Note: currently this does nothing (disabled since the move to the new graphics). */
void rotate_y_scene(int nsteps, float stepsize);
/*! \brief rotate view round x axis stepsize degrees for nstep such steps

  Note: currently this does nothing (disabled since the move to the new graphics). */
void rotate_x_scene(int nsteps, float stepsize);
/*! \brief rotate view round z axis stepsize degrees for nstep such steps

  Note: currently this does nothing (disabled since the move to the new graphics). */
void rotate_z_scene(int nsteps, float stepsize);

/*! \brief Bells and whistles rotation

    spin, zoom and translate.

    where axis is either x,y or z (1, 2 or 3),
    stepsize is in degrees,
    zoom_by and x_rel etc are how much zoom, x,y,z should
            have changed by after nstep steps.

    Note: currently this does nothing (disabled since the move to the new graphics).
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
/*! \brief translate rotation centre relative to screen axes for nsteps

  Note: currently a placeholder that does nothing. */
void translate_scene_x(int nsteps);
/*! \brief translate rotation centre relative to screen axes for nsteps

  Note: currently a placeholder that does nothing. */
void translate_scene_y(int nsteps);
/*! \brief translate rotation centre relative to screen axes for nsteps

  Note: currently a placeholder that does nothing. */
void translate_scene_z(int nsteps);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Views                                                   */
/*  ----------------------------------------------------------------------- */
/* section Views Interface */
/*! \name Views Interface */
/*! \{ */
/*! \brief add the current view (rotation centre, orientation and zoom)
  to the list of views

  @param view_name the name of the view
  @return the view number (the index of the new view) */
int add_view_here(const char *view_name);
/*! \brief add a view from explicit values

  @param rcx the x coordinate of the rotation centre
  @param rcy the y coordinate of the rotation centre
  @param rcz the z coordinate of the rotation centre
  @param quat1 the first quaternion component, passed in glm::quat (w, x, y, z) order
  @param quat2 the second quaternion component
  @param quat3 the third quaternion component
  @param quat4 the fourth quaternion component
  @param zoom the zoom
  @param view_name the name of the view
  @return the view number (the index of the new view) */
int add_view_raw(float rcx, float rcy, float rcz, float quat1, float quat2,
		 float quat3, float quat4, float zoom, const char *view_name);
/*! \brief play the views - animate (interpolate) from each view to the next

  The speed is controlled by set_views_play_speed(). */
void play_views();
/*! \brief remove the view that matches the current view (if any) */
void remove_this_view();
/*! \brief remove the (first) view with the given name

  @return 0 (always) */
int remove_named_view(const char *view_name);
/*! \brief remove the given view number (nothing happens if there is no such view) */
void remove_view(int view_number);
/*! \brief go to the first view.

  @param snap_to_view_flag if 1 go directly, else move along a smooth path
  @return 0 (always) */
int go_to_first_view(int snap_to_view_flag);
/*! \brief go to the given view number.

  @param view_number the view number
  @param snap_to_view_flag if 1 go directly, else move along a smooth path
  @return 0 (always) */
int go_to_view_number(int view_number, int snap_to_view_flag);
/*! \brief add a spin view (a rotation about the screen y axis) to the list of views

  @param view_name the name of the view
  @param n_steps the number of steps
  @param degrees_total the total rotation in degrees
  @return the view number of the new view */
int add_spin_view(const char *view_name, int n_steps, float degrees_total);
/*! \brief Add a view description/annotation to the give view number */
void add_view_description(int view_number, const char *description);
/*! \brief add a view (not add to an existing view) that *does*
  something (e.g. displays or undisplays a molecule) rather than move
  the graphics.

  The action function (a string) is stored with the view;
  play_views() and go_to_view_number() do not move the graphics for an
  action view.

  @return the view number for this (new) view.
 */
int add_action_view(const char *view_name, const char *action_function);
/*! \brief add an action view after the view of the given view number

  If view_number is beyond the end of the list, the view is added at the end.

  @return the view number for this (new) view.
 */
int insert_action_view_after_view(int view_number, const char *view_name, const char *action_function);
/*! \brief return the number of views */
int n_views();

/*! \brief save views to view_file_name

  The views are written as script commands (nothing is written if there are no views). */
void save_views(const char *view_file_name);

/*! \brief return the views play speed (default 10.0) */
float views_play_speed();
/*! \brief set the views play speed - larger is faster (default 10.0) */
void set_views_play_speed(float f);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief return the name of the given view, if view_number does not
  specify a view return scheme value False */

SCM view_name(int view_number);
/*! \brief return the description of the given view, or scheme false if there
  is no such view or it has no description */
SCM view_description(int view_number);
/*! \brief go to the given view (not implemented for scheme) */
void go_to_view(SCM view);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
/*! \brief return the name of the given view, if view_number does not
  specify a view return Python value False */
PyObject *view_name_py(int view_number);
/*! \brief return the description of the given view, or Python False if there
  is no such view or it has no description */
PyObject *view_description_py(int view_number);
/*! \brief animate to the given view

  @param view a list [quaternion, rotation_centre, zoom, name] where
  quaternion is a list of 4 numbers, rotation_centre a list of 3
  numbers and name a string */
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
/*! \brief set the movie frame file name prefix (default "movie_")

  Frame images are written as prefix + 5-digit frame number + ".ppm". */
void set_movie_file_name_prefix(const char *file_name);
/*! \brief set the number of the next movie frame (default 0) */
void set_movie_frame_number(int frame_number);
#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
/*! \brief return the movie file name prefix */
SCM movie_file_name_prefix();
#endif
#ifdef USE_PYTHON
/*! \brief return the movie file name prefix */
PyObject *movie_file_name_prefix_py();
#endif /* USE_PYTHON */
#endif /* c++ */
/*! \brief return the number of the next movie frame */
int movie_frame_number();
/*! \brief set movie mode: when on, an image is dumped for every redraw

  @param make_movies_flag 1 for on, 0 for off (default) */
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
  after a ligand fit

  @param f the fraction (default 0.75); values outside 0 to 1 are ignored */
void set_ligand_acceptable_fit_fraction(float f);

/*! \brief set the default sigma level that the map is searched to
  find potential ligand sites

  @param f the level in map rmsd (sigma) units (default 1.0) */
void set_ligand_cluster_sigma_level(float f); /* default 1.0 */

/*! \brief set the number of conformation samples

    big ligands require more samples.  Default 50.*/
void set_ligand_flexible_ligand_n_samples(int i); /* default 50: Really? */
/*! \brief set verbose reporting for ligand fitting

  @param i 1 for on, 0 for off (default) */
void set_ligand_verbose_reporting(int i); /* 0 off (default), 1 on */

/*! \brief search the top n sites for ligands.

   Default 10. */
void set_find_ligand_n_top_ligands(int n); /* fit the top n ligands,
					      not all of them, default
					      10. */

/*! \brief real-space refine the ligand solutions after fitting?

  @param state 1 for on (default), 0 for off */
void set_find_ligand_do_real_space_refinement(short int state);

/*! \brief allow multiple ligand solutions per cluster.

The first limit is the fraction of the top scored positions that go on
to correlation scoring (closer to 1 means less and faster - default
0.7).

The second limit is the fraction of the top correlation score that is
considered interesting.  Limits the number of solutions displayed to
user. Default 0.9. (Note: this value is currently stored but not used:
the search uses a fixed value of 0.9.)

There is currently no chi-angle set redundancy filtering - I suspect
that there should be.

Nino-mode.

*/
void set_find_ligand_multi_solutions_per_cluster(float lim_1, float lim_2);

/*! \brief how shall we treat the waters during ligand fitting?

   pass with istate=1 for waters to mask the map in the same way that
   protein atoms do (default 0).
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

/*! \brief  Allow the user a scripting means to find ligand at the rotation centre

  @param state 1: fit only to the cluster at the rotation centre (rather than
  searching the whole map), 0: search the map (default) */
void set_find_ligand_here_cluster(int state);

/*! \brief run the ligand search

  Uses the protein, map and ligand molecules set by
  set_ligand_search_protein_molecule(), set_ligand_search_map_molecule()
  and add_ligand_search_ligand_molecule() (or
  add_ligand_search_wiggly_ligand_molecule()). A masked map molecule is
  created and each solution is added as a new molecule. */
void execute_ligand_search();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief run the ligand search and return a list of the new
  (solution) molecule indices */
SCM execute_ligand_search_scm();
#endif
#ifdef USE_PYTHON
/*! \brief run the ligand search (see execute_ligand_search())

  @return a list of the molecule indices of the solutions (possibly empty) */
PyObject *execute_ligand_search_py();
#endif /* USE_PYTHON */
#endif /* __cplusplus */
/*! \brief clear the list of ligands to search for in ligand searching */
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

  Each residue in imol_ligand is graph-matched (not by atom names) onto
  the residue specified by the reference parameters; the best match is
  moved onto the reference and all other residues in imol_ligand are
  deleted.

@return success status, False = failed to find residue in either
imol_ligand or imo_ref.  If success, return the RT operator and
the match info (distance score and number of matched atoms).
*/
SCM overlap_ligands(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
/*! \brief like overlap_ligands() but without moving or deleting anything */
SCM analyse_ligand_differences(int imol_ligand, int imol_ref, const char *chain_id_ref,
			       int resno_ref);
/*! \brief compare the dictionary atom (energy) types of graph-matched
  atoms of the first residue of imol_ligand and the reference residue
  (for testing). Return the number of mismatches, or scheme false. */
SCM compare_ligand_atom_types_scm(int imol_ligand, int imol_ref, const char *chain_id_ref,
				  int resno_ref);
#endif /* USE_GUILE */
/*! \brief rotate the torsions of the first residue of imol_ligand to
  match those of the reference residue

  The torsions come from the dictionary of the reference residue type,
  which must be available.

  @param imol_ligand the molecule index of the ligand to be changed
  @param imol_ref the molecule index of the reference residue
  @param chain_id_ref the chain id of the reference residue
  @param resno_ref the residue number of the reference residue */
void match_ligand_torsions(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
#ifdef USE_PYTHON
/*! \brief Overlap residue with "template"-based matching.

  Each residue in imol_ligand is graph-matched (not by atom names) onto
  the residue specified by the reference parameters. The best match is
  moved onto the reference and all other residues in imol_ligand are
  deleted.

  @param imol_ligand the molecule index of the ligand(s) to be moved
  @param imol_ref the molecule index of the reference residue
  @param chain_id_ref the chain id of the reference residue
  @param resno_ref the residue number of the reference residue
  @return False on failure, else [rtop, [dist_score, n_matched_atoms]]
  where rtop is [rotation_matrix_as_9_numbers, translation_as_3_numbers] */
PyObject *overlap_ligands_py(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
/*! \brief like overlap_ligands_py() but without moving or deleting anything

  @return False on failure, else [rtop, [dist_score, n_matched_atoms]] */
PyObject *analyse_ligand_differences_py(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref);
/*! \brief compare the dictionary atom (energy) types of graph-matched
  atoms of the first residue of imol_ligand and the reference residue

  For testing that pyrogen generates consistent atom types.

  @return the number of mismatched atom types, or False on failure (or if
  no atoms were matched) */
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



/* Just pondering - just a stub currently (it returns -1).

   For use with exporting ligands from the 2D sketcher to the main
   coot window. It should find the residue that this residue is
   sitting on top of that is in a molecule that has lot of atoms
   (i.e. is a protein) and create a new molecule that is a copy of the
   molecule without the residue/ligand that this (given) ligand
   overlays - and a copy of this given ligand.  */
int exchange_ligand(int imol_lig, const char *chain_id_lig, int resno_lig, const char *ins_code_lig);



/*! \brief flip the ligand (usually active residue) around its eigen vectors
   to the next flip number.  Immediate replacement (like flip
   peptide).

   Repeated calls cycle through the 4 flip states. */
void flip_ligand(int imol, const char *chain_id, int resno);

/*! \brief JED-flip: flip a fragment of a residue around a rotatable bond

  Uses the non-CONST, non-ring dictionary torsions that involve the given
  atom (as the second or third torsion atom), choosing the one with the
  smallest fragment. That fragment is rotated by 360/period degrees.
  Any problem is reported in the status bar.

  @param imol the model molecule index
  @param chain_id the chain id
  @param res_no the residue number
  @param ins_code the insertion code
  @param atom_name the atom name (as in the PDB file, e.g. " C5 ")
  @param alt_conf the alt conf of the atom
  @param invert_selection if 1, move the other fragment */
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

/*! \brief find waters

  This is find_waters() with show_blobs_dialog = 1. */
void execute_find_waters_real(int imol_for_map,
			      int imol_for_protein,
			      short int new_waters_mol_flag,
			      float rmsd_cut_off);

/*! \brief find waters

  The map is masked by imol_for_protein (including its waters) and
  waters are placed in the peaks above the cut-off (using the water
  distance limits, variance limit and number of cycles set by the
  functions below).

  @param imol_for_map the map molecule index
  @param imol_for_protein the model molecule index
  @param new_waters_mol_flag 1: put the waters in a new molecule, 0: add them to imol_for_protein
  @param rmsd_cut_off the cut-off level in map rmsd (sigma) units
  @param show_blobs_dialog 1 to show the "big blobs" (unmodelled density) results in graphical mode */
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
  move_waters_to_around_protein() funtion.

  @return the distance (A), or -1 on failure (e.g. bad imol, no waters) */
float max_water_distance(int imol);

/*! \brief return the find-waters sigma cut-off as text (newly allocated - the caller should free it) */
char *get_text_for_find_waters_sigma_cut_off();
/*! \brief set the find-waters sigma cut-off (default 1.4) */
void set_value_for_find_waters_sigma_cut_off(float f);

/*! \brief set the limit of interesting variance, above which waters
  are listed (otherwise ignored)

default 0.12. */
void set_water_check_spherical_variance_limit(float f);

/*! \brief set water to protein distance limits for water finding

 f1 is the minimum distance (default 2.4 A), f2 is the maximum distance (default 3.2 A) */
void set_ligand_water_to_protein_distance_limits(float f1, float f2);

/*! \brief set the number of cycles of water searching (default 3) */
void set_ligand_water_n_cycles(int i);
/*! \brief write the raw peak-searched water positions in subsequent water finding (for debugging) */
void set_write_peaksearched_waters();

/*! \brief find blobs
 *
 * Not useful for MCP. For interactive use only.
 *
 * @param imol_model the model molecule index (used to mask the map)
 * @param imol_for_map the map molecule index
 * @param cut_off the cut-off level in map rmsd (sigma) units
 * @param interactive_flag 1 to show the blobs dialog
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

/*! \brief set the default thickness for bonds (e.g. in ~/.coot)

  @param t the thickness (default 5) */
void set_default_bond_thickness(int t);

/*! \brief set the thickness of the bonds in molecule number imol to t pixels  */
void set_bond_thickness(int imol, float t);
/*! \brief set the thickness of the bonds of the intermediate atoms to t pixels  */
void set_bond_thickness_intermediate_atoms(float t);

/*! \brief allow lines that are further away to be thinner

  @param state 1 for on, 0 for off (default) */
void set_use_variable_bond_thickness(short int state);

/*! \brief set bond colour for molecule

  @param imol the model molecule index
  @param f the colour rotation (a hue rotation in degrees, 0 to 360) */
void set_bond_colour_rotation_for_molecule(int imol, float f);

/*! \brief set default for the drawing of atoms in stick mode (default is on (1)) */
void set_draw_stick_mode_atoms_default(short int state);


/*! \brief get the bond colour for molecule.

Return -1 on err (bad molecule number) */
float get_bond_colour_rotation_for_molecule(int imol);

/*! \brief set the size of the stars drawn for unbonded atoms (default 0.5) */
void set_unbonded_atom_star_size(float f);

/*! \brief set the default representation type (default 1).

  The type is a bond-colour mode (bonds box type): e.g. 1 for normal bonds,
  2 for CA bonds, 3 for colour by chain. */
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

/*! \brief get hydrogen bonds
 *
 * \details
 * For the returned value, an "atom" here looks like:
 *
 *    PyDict_SetItemString(at_py, "x", PyFloat_FromDouble(at->x));
 *    PyDict_SetItemString(at_py, "y", PyFloat_FromDouble(at->y));
 *    PyDict_SetItemString(at_py, "z", PyFloat_FromDouble(at->z));
 *    PyDict_SetItemString(at_py, "charge",       PyFloat_FromDouble(at->charge));
 *    PyDict_SetItemString(at_py, "occ",          PyFloat_FromDouble(at->occupancy));
 *    PyDict_SetItemString(at_py, "b_iso",        PyFloat_FromDouble(at->tempFactor));
 *    PyDict_SetItemString(at_py, "element",      myPyString_FromString(at->element));
 *    PyDict_SetItemString(at_py, "name",         myPyString_FromString(at->name));
 *    PyDict_SetItemString(at_py, "model",        PyFloat_FromDouble(at->GetModelNum()));
 *    PyDict_SetItemString(at_py, "chain",        myPyString_FromString(at->GetChainID()));
 *    PyDict_SetItemString(at_py, "altLoc",       myPyString_FromString(at->altLoc));
 *    PyDict_SetItemString(at_py, "residue_name", myPyString_FromString(at->GetResidue()->GetResName()));
 *
 *    For the returned value, a hydrogen bond looks like this:
 *
 *    PyList_SetItem(l, 0, hb_hydrogen_py);       // an atom
 *    PyList_SetItem(l, 1, donor_py);             // an atom
 *    PyList_SetItem(l, 2, acceptor_py);          // an atom
 *    PyList_SetItem(l, 3, donor_neigh_py);       // an atom, possibly None
 *    PyList_SetItem(l, 4, acceptor_neigh_py);    // an atom, possibly None
 *
 *    PyList_SetItem(l, 5, PyFloat_FromDouble(h_bond.angle_1));
 *    PyList_SetItem(l, 6, PyFloat_FromDouble(h_bond.angle_2));
 *    PyList_SetItem(l, 7, PyFloat_FromDouble(h_bond.angle_3));
 *    PyList_SetItem(l, 8, PyFloat_FromDouble(h_bond.dist));
 *
 *    PyList_SetItem(l,  9, PyBool_FromLong(h_bond.ligand_atom_is_donor));
 *    PyList_SetItem(l, 10, PyBool_FromLong(h_bond.hydrogen_is_ligand_atom));
 *    PyList_SetItem(l, 11, PyBool_FromLong(h_bond.bond_has_hydrogen_flag));
 *
 * @param imol the molecule index
 * @param selection_1  the atom selection of the "from" atoms
 * @param selection_2  the atom selection of the "to" atoms (treated as the "ligand").
 *              Note that often atom_selection_1 and atom_selection_2 are the same,
 *              e.g. "//A"
 * @param mcdonald_and_thornton_algoritnm use 0 if the model does not have hydrogen atoms
                                          use 1 if the model has hydrogen atoms.
 * @return the hydrogen bonds as a python list object, or False if
 *         imol is not a valid model molecule
 *
 */
PyObject *get_hydrogen_bonds_py(int imol, const char *selection_1, const char *selection_2, short int mcdonald_and_thornton_algoritnm);

/*! \brief draw little coloured balls on atoms

turn off with state = 0

turn on with state = 1 */
void set_draw_stick_mode_atoms(int imol, short int state);

/*! \brief set the state for drawing missing resiude loops
 *
 * Used, for example, when taking screenshots, we often
 * don't want to see them in such cases.
 * Or maybe there's just too many of them to be useful
 *
 * @param state the draw state (0 for "off", 1 for "on" (default))
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
/*! \brief draw molecule number imol all atoms coloured by user-defined atom colours
 *
 * Use this function after using set_user_defined_atom_colour_by_selection_py()
 * and/or set_user_defined_atom_colour_py().
 * When atom selection colouring has been created or updated, then calling this function
 * actually forces the regeneration and drawing of the molecule with the new colour scheme.
 *
 * @param imol the molecule index
 * */
void graphics_to_user_defined_atom_colours_all_atoms_representation(int imol);
/*! \brief what is the bond drawing state of molecule number imol

  @return the bond-colour mode (bonds box type, e.g. 1 for normal bonds,
  2 for CA bonds), or -1 if imol is not a valid model molecule */
int get_graphics_molecule_bond_type(int imol);
/*! \brief scale the colours for colour by b factor representation

  @return 1 on success, 0 if imol is not a valid model molecule */
int set_b_factor_bonds_scale_factor(int imol, float f);
/*! \brief change the representation of the model molecule closest to
  the centre of the screen

  The molecule is that of the active atom.

  @param up_or_down 1 to step up, -1 to step down through the representations */
void change_model_molecule_representation_mode(int up_or_down);

/* not today void set_ca_bonds_loop_params(float p1, float p2, float p3); */

/*! \brief make the carbon atoms for molecule imol be grey

  (or the colour set by set_grey_carbon_colour())

  @param imol the model molecule index
  @param state 1 for on, 0 for off
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

e.g. (make-ball-and-stick 0 "/1" 0.15 0.25 1)

Note: this uses the old display-list graphics and currently has no visible
effect - use set_model_molecule_representation_style() or additional
representations instead.

@return imol */
int make_ball_and_stick(int imol,
			const char *atom_selection_str,
			float bond_thickness, float sphere_size,
			int do_spheres_flag);
/*! \brief clear ball and stick representation of molecule number imol

  (the old display-list ball and stick objects)

  @return 0 */
int clear_ball_and_stick(int imol);

/*! \brief set the model molecule representation stye 0 for ball-and-stick/licorice (default) and 1 for ball

  (2 for van der Waals balls) */
void set_model_molecule_representation_style(int imol, unsigned int mode);

/*! \brief set show a ribbon/mesh for a given molecule

  @param imol the model molecule index
  @param mesh_index the index of the molecular representation mesh of that molecule
  @param state 1 to show, 0 to hide */
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

/*! \brief add an additional representation for a residue range

  @param imol the model molecule index
  @param chain_id the chain id
  @param resno_start the first residue number
  @param resno_end the last residue number
  @param ins_code the insertion code
  @param representation_type 0: simple lines, 1: sticks, 2: ball and stick,
         3: liquorice, 4: surface (see coot::additional_representations_t)
  @param bonds_box_type the bond colour mode (see the list above, e.g. 1 for normal bonds)
  @param bond_width the bond width
  @param draw_hydrogens_flag 1 to draw hydrogens, 0 not
  @return the index of the additional representation, -1 on error.
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
/*! \brief return information about the additional representations of
  molecule number imol, or scheme false if imol is not a valid model molecule */
SCM additional_representation_info_scm(int imol);
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
/*! \brief return information about the additional representations of a molecule

  @return a list of [index, info_string, is_shown, bond_width] items, or
  Python False if imol is not a valid model molecule */
PyObject *additional_representation_info_py(int imol);
#endif	/* USE_PYTHON */

#endif	/* __cplusplus */

/* Turn on nice animated ligand interaction display.

turn on with arg 1.

turn off with arg 0.

(Note: in GTK4 Coot turning them on currently does nothing.) */
void set_flev_idle_ligand_interactions(int state);

/* Toggle for animated ligand interaction display above */
void toggle_flev_idle_ligand_interactions();

/*! \brief calculate the hydrogen bonds of molecule number imol (all atoms)
  and make a mesh of them for display (see set_draw_hydrogen_bonds()) */
void calculate_hydrogen_bonds(int imol);

/*! \brief set the drawing state of the hydrogen bonds made by calculate_hydrogen_bonds()

  @param state 1 for on, 0 for off */
void set_draw_hydrogen_bonds(int state);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  dots display                                            */
/*  ----------------------------------------------------------------------- */

/*! \name Dots Representation */
/*! \{ */
/*! \brief display a dotted (van der Waals) surface for an atom selection

  Dots are placed on a sphere around each selected atom (radius by
  element) and dots that lie inside another selected atom are omitted.

  @param imol the model molecule index
  @param atom_selection_str an mmdb atom selection string, e.g. "//A/10-20"
  @param dots_object_name the name for the dots object (can be used with
         \c clear_dots_by_name())
  @param dot_density the dot density: 1.0 gives dots every 5 degrees,
         larger values give denser dots
  @param sphere_size_scale currently not used (the radius is not scaled)
  @return the dots handle (for use with \c clear_dots()), or -1 on failure
          (e.g. invalid molecule) */
int dots(int imol,
	 const char *atom_selection_str,
	 const char *dots_object_name,
	 float dot_density, float sphere_size_scale);


/*! \brief set the colour of the surface dots of the imol-th molecule
  to be the given single colour

  The colour is applied when dots are created, so this affects dots
  made after this call, not existing ones.

  r,g,b are values between 0.0 and 1.0 */
void set_dots_colour(int imol, float r, float g, float b);

/*! \brief no longer set the dots of molecule imol to a single colour

i.e. go back to element-based colours (for dots created after this call). */
void unset_dots_colour(int imol);

/*! \brief clear dots in imol with dots_handle

  @param imol the model molecule index
  @param dots_handle the handle returned by \c dots() */
void clear_dots(int imol, int dots_handle);

/*! \brief clear the first dots object for imol with given name */
void clear_dots_by_name(int imol, const char *dots_object_name);

/*! \brief return the number of dots sets for molecule number imol

  Cleared dots sets are still counted.

  @return the number of dots sets, or -1 if imol is out of range */
int n_dots_sets(int imol);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  pepflip                                                 */
/*  ----------------------------------------------------------------------- */
/* section Pep-flip Interface */
/*! \name Pep-flip Interface */
/*! \{ */
/*! \brief set up (state 1) or cancel (state 0) a pepflip, ready for an atom pick (GUI use) */
void do_pepflip(short int state); /* sets up pepflip, ready for atom pick. */

/*! \brief pepflip (flip the peptide) of the given residue
 *
 *  Rotate the carbonyl C and O atom of this residue and the N of the
 *  next residue around a vector between the two CA atoms by 180 degrees.
 *  This is often a useful modelling operation to create a different hypothesis
 *  about the orientation of the main-chain atoms - that can then be used
 *  for refinement. This can sometimes allow the model to be removed from
 *  local minima of backbone conformations.
 *
 *  @param imol is the index of the model molecule
 *  @param chain_id is the chain-id
 *  @param resno is the residue number (the residue that has the C and O atoms)
 *  @param inscode the insertion code (typically "")
 *  @param altconf the altconf (typically "")
 *
 */
void pepflip(int imol, const char *chain_id, int resno, const char *inscode,
	     const char *altconf);

/*! \brief pepflip the intermediate (refining) atoms near the screen centre

  Finds the intermediate atom closest to the rotation centre (within 2 Å).
  If that atom is an N, the peptide between the previous residue and this
  one is flipped, otherwise the peptide between this residue and the next.
  The running refinement is then restarted.

  @return 1 if a flip was made, 0 otherwise */
int pepflip_intermediate_atoms();

/*! \brief pepflip the "other" peptide of the intermediate-atoms residue near the screen centre

  As \c pepflip_intermediate_atoms(), but flips the peptide on the other
  side of the residue of the closest intermediate atom.

  @return 1 if a flip was made, 0 otherwise */
int pepflip_intermediate_atoms_other_peptide();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief suggest peptide flips using a difference map

  @param imol_coords the model molecule index
  @param imol_difference_map the index of a difference map
  For each peptide, the difference-map value at the flipped O position
  minus that at the current O position is compared with the distribution
  of such differences for random pairs of points; a flip is suggested
  when it exceeds mean + n_sigma * sd of that distribution.

  @param n_sigma the cut-off, in standard deviations of that distribution
  @return a list of residue specs for suggested flips (empty if the inputs
          are invalid or imol_difference_map is not a difference map) */
SCM pepflip_using_difference_map_scm(int imol_coords, int imol_difference_map, float n_sigma);
#endif
#ifdef USE_PYTHON
/*! \brief suggest peptide flips using a difference map

  The model is not changed.

  @param imol_coords the model molecule index
  @param imol_difference_map the index of a difference map
  For each peptide, the difference-map value at the flipped O position
  minus that at the current O position is compared with the distribution
  of such differences for random pairs of points; a flip is suggested
  when it exceeds mean + n_sigma * sd of that distribution.

  @param n_sigma the cut-off, in standard deviations of that distribution
  @return a list of residue specs for suggested flips (empty if the inputs
          are invalid or imol_difference_map is not a difference map) */
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
/*! \brief set up (state 1) or cancel (state 0) rigid body refinement, ready for 2 atom picks (GUI use) */
void do_rigid_body_refine(short int state);	/* set up for atom picking */

/*! \brief rigid body refine a residue range

   Sets the residue-range atom selection from the arguments and then
   calls \c execute_rigid_body_refine(). The fit is made against the
   refinement map (see \c set_imol_refinement_map()); if that is not set
   nothing is fitted. The fitted atoms become intermediate atoms that
   need to be accepted, unless immediate replacement is on.

   @param imol the model molecule index
   @param chain_id the chain id
   @param reso_start the first residue number of the range
   @param resno_end the last residue number of the range */
void rigid_body_refine_zone(int imol, const char *chain_id, int reso_start, int resno_end);

/*! \brief rigid body refine the atoms of an atom selection

   The selected atoms are fitted as a rigid body into the refinement map
   (see \c set_imol_refinement_map()), with the map masked by the
   non-selected atoms. The result is presented as intermediate atoms to
   be accepted, unless immediate replacement is on.

   @param imol the model molecule index
   @param atom_selection_string an mmdb atom selection string, e.g. "//A/10-20" */
void
rigid_body_refine_by_atom_selection(int imol, const char *atom_selection_string);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief rigid body refine using residue ranges.  residue_ranges is
    a list of residue ranges.  A residue range is (list chain-id
    resno-start resno-end).

    The ranges are fitted together as one rigid body against the
    refinement map.

    @return \#t on success, \#f on failure */
SCM rigid_body_refine_by_residue_ranges_scm(int imol, SCM residue_ranges);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief rigid body refine using residue ranges.  residue_ranges is
    a list of residue ranges.  A residue range is [chain_id,
    resno_start, resno_end].

    The ranges are fitted together as one rigid body against the
    refinement map (see \c set_imol_refinement_map()).

    @return True on success, False on failure (e.g. no refinement map,
            bad input or no fit found) */
PyObject *
rigid_body_refine_by_residue_ranges_py(int imol, PyObject *residue_ranges);
#endif /* USE_PYTHON */
#endif /* __cplusplus */

/*! \brief run the rigid body refinement after the atoms have been picked (GUI use)

   @param auto_range_flag if 1, the range is determined automatically from
          the first picked atom; if 0, the range is between the 2 picked atoms */
void execute_rigid_body_refine(short int auto_range_flag); /* atom picking has happened.
				     Actually do it */


/*! \brief set rigid body fraction of atoms in positive density

 The minimum fraction of atoms that must be in positive density for a
 rigid body fit to be accepted.

 @param f in the range 0.0 -> 1.0 (default 0.75); values outside this
        range are ignored */
void set_rigid_body_fit_acceptable_fit_fraction(float f);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  dynamic map                                             */
/*  ----------------------------------------------------------------------- */
/* section Dynamic Map */
/*! \name Dynamic Map */
/*! \{ */
/*! \brief toggle dynamic map size display (see \c set_dynamic_map_size_display_on()) */
void   toggle_dynamic_map_display_size();
/*! \brief toggle dynamic map sampling (see \c set_dynamic_map_sampling_on()) */
void   toggle_dynamic_map_sampling();
/* scripting interface: */
/*! \brief turn on dynamic map size display

  When dynamic map sampling is also on, the map contouring radius is
  scaled up by the sampling step, so that a coarser map covers a larger
  region when zoomed out. Default off. */
void set_dynamic_map_size_display_on();
/*! \brief turn off dynamic map size display */
void set_dynamic_map_size_display_off();
/*! \brief return the dynamic map size display state (1 on, 0 off) */
int get_dynamic_map_size_display();
/*! \brief turn on dynamic map sampling

  The map is contoured with a coarser grid sampling step as the view is
  zoomed out. Default off. */
void set_dynamic_map_sampling_on();
/*! \brief turn off dynamic map sampling */
void set_dynamic_map_sampling_off();
/*! \brief return the dynamic map sampling state (1 on, 0 off) */
int get_dynamic_map_sampling();
/*! \brief set the dynamic map zoom offset

  The offset is added to the zoom when calculating the dynamic sampling
  step (step = 1 + int(0.009 * (zoom + offset)), max 15). Default 0.

  @param i the zoom offset */
void set_dynamic_map_zoom_offset(int i);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  build one residue by phi/psi search                     */
/*  ----------------------------------------------------------------------- */
/* section Add Terminal Residue Functions */
/*! \name Add Terminal Residue Functions */
/*! \{ */
/*! \brief set up (state 1) or cancel (state 0) terminal residue addition, ready for an atom pick (GUI use)

  If the refinement map is not set, the map-selection dialog is shown instead. */
void do_add_terminal_residue(short int state);
/*  execution of this is in graphics_info_t because it uses a mmdb::Residue */
/*  in the interface and we can't have that in c-interface.h */
/*  (compilation of coot_wrap_guile goes mad on inclusion of */
/*   mmdb_manager.h) */
/*! \brief set the number of random phi/psi trials for terminal residue addition

  @param n the number of trials (default 5000) */
void set_add_terminal_residue_n_phi_psi_trials(int n);
/* Add Terminal Residues actually build 2 residues, this allows us to
   see both residues - default is 0 (off). */
/*! \brief keep both residues built by terminal residue addition

  The phi/psi search actually builds 2 residues; with this flag set both
  are added to the model (the second one only if a residue with that
  number does not already exist).

  @param i 1 for on, 0 for off (default 0) */
void set_add_terminal_residue_add_other_residue_flag(int i);
/*! \brief set the add-terminal-residue rigid body refine flag

  Note that this flag is currently not used by the terminal residue
  addition code.

  @param v 1 for on, 0 for off (default 0) */
void set_add_terminal_residue_do_rigid_body_refine(short int v);
/*! \brief deprecated - use \c set_add_terminal_residue_do_rigid_body_refine() */
void set_terminal_residue_do_rigid_body_refine(short int v); /* remove this for 0.9, wraps above */
/*! \brief write out the trial solutions of terminal residue addition as PDB files (for debugging)

  @param debug_state 1 for on, 0 for off (default 0) */
void set_add_terminal_residue_debug_trials(short int debug_state);
/*! \brief return the state of the immediate addition flag for terminal residues (1 on, 0 off) */
int add_terminal_residue_immediate_addition_state();

/*! \brief set immediate addition of terminal residue

call with i=1 for immediate addition (the default), 0 means the new
residue is shown as intermediate atoms to be accepted */
void set_add_terminal_residue_immediate_addition(int i);

/*! \brief Add a terminal residue

  The new residue is fitted to the refinement map by a random phi/psi
  search, so the refinement map must be set (see
  \c set_imol_refinement_map()). The residue given must be at a terminus
  and must have a CA atom.

   @param imol the model molecule index
   @param chain_id the chain id
   @param residue_number the residue number of the existing terminal residue
   @param residue_type the 3-letter code of the new residue, or "auto",
          which uses the sequence assigned to the chain (if any), otherwise ALA
   @param immediate_add is recommended to be 1. If 0, the new residue
          is shown as intermediate atoms to be accepted.
   @return 0 on failure, 1 on success
*/
int add_terminal_residue(int imol, const char *chain_id, int residue_number,
                          const char *residue_type, int immediate_add);

/*! \brief Add a residue to a chain or at the end of a fragment

  This can be used to fill a gap of one residue or to fill a gap
  of multiple residues by being called several times. Probably
  RSR refinement would be useful after each call to this function in
  such a case.

  This is a synonym for \c add_terminal_residue() (and so needs the
  refinement map to be set).

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

No fitting is done: an ideal A-form (RNA) or B-form (DNA) single-stranded
nucleotide is added to the terminus. The base type is taken from the
sequence assigned to the chain if available.

  @param imol the model molecule index
  @param chain_id the chain id
  @param res_no the residue number of the existing terminal nucleotide
  @return 1 if imol is a valid model molecule (even if no nucleotide
          could be added), 0 otherwise
*/
int add_nucleotide(int imol, const char *chain_id, int res_no);


/*! \brief Add a terminal residue using given phi and psi angles

  No map fitting is done. The new residue is built as ALA (mainchain
  and CB) at the N- or C-terminus as appropriate.

  @param imol the molecule index
  @param chain_id the chain ID
  @param res_no the residue number (of the existing residue to attach to)
  @param residue_type is currently ignored (an ALA is built)
  @param phi is phi in degrees
  @param psi is psi in degrees
  @return the success status, 0 on failure, 1 on success
 */
int add_terminal_residue_using_phi_psi(int imol, const char *chain_id, int res_no,
				       const char *residue_type, float phi, float psi);

/*! \brief set the residue type of an added terminal residue.

  @param type a 3-letter code, or "auto" (the default) to use the
         assigned sequence (ALA if none) */
void set_add_terminal_residue_default_residue_type(const char *type);
/*! \brief set a flag to run refine zone on terminal residues after an
  addition.

  This applies to interactive (clicked) terminal residue additions.

  @param istat 1 for on, 0 for off (default 0) */
void set_add_terminal_residue_do_post_refine(short int istat);
/*! \brief what is the value of the previous flag? */
int add_terminal_residue_do_post_refine_state();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief find the residue type of a residue to be added at a terminus

  Uses the sequence alignment of the chain to the sequence assigned to it.

  @param imol the model molecule index
  @param chain_id the chain id
  @param resno the residue number of the residue to be added
  @return the 3-letter residue type, or \#f if it cannot be determined */
SCM find_terminal_residue_type(int imol, const char *chain_id, int resno);
#endif
#ifdef USE_PYTHON
/*! \brief find the residue type of a residue to be added at a terminus

  Uses the sequence alignment of the chain to the sequence assigned to it.

  @param imol the model molecule index
  @param chain_id the chain id
  @param resno the residue number of the residue to be added
  @return the 3-letter residue type, or False if it cannot be determined */
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

  If the residue given by residue_spec does not exist it is created
  (and its chain too, if needed); the atoms are then added to it.

  @param imol the model molecule index
  @param residue_spec the residue spec of the residue to add to or create
  @param res_name the residue name (used if the residue is created)
  @param list_of_atoms a list of atoms, each of the form
         [[atom_name, alt_conf], [occupancy, b_factor, element, segid], [x, y, z]]
         (the same form as an item of \c residue_info_py(), with an
         isotropic B-factor; an optional 4th item is ignored). Atoms not
         in this form are skipped.
  @return intended to be the number of atoms added, but currently always 0
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
/*! \brief delete the atom with the given atom index (internal/GUI use)

  @param imol the model molecule index
  @param index the atom index in the molecule's atom selection
  @param do_delete_dialog currently ignored */
void delete_atom_by_atom_index(int imol, int index, short int do_delete_dialog);
/*! \brief delete the residue of the atom with the given atom index (internal/GUI use)

  If the atom has an alt conf (or the molecule has more than one model)
  only the atoms of that alt conf in that model are deleted.

  @param imol the model molecule index
  @param index the atom index in the molecule's atom selection
  @param do_delete_dialog currently ignored */
void delete_residue_by_atom_index(int imol, int index, short int do_delete_dialog);
/*! \brief delete the hydrogen atoms of the residue of the atom with the given atom index (internal/GUI use)

  @param imol the model molecule index
  @param index the atom index in the molecule's atom selection
  @param do_delete_dialog currently ignored */
void delete_residue_hydrogens_by_atom_index(int imol, int index, short int do_delete_dialog);

/*! \brief delete residue range

  All atoms (all alt confs) of the residues in the range are deleted.
  The start and end may be given in either order.

  @param imol the model molecule index
  @param chain_id the chain id
  @param resno_start the first residue number of the range
  @param end_resno the last residue number of the range */
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

/*! \brief delete the atoms of a residue that have the given alt conf

  Only atoms whose alt conf matches altloc exactly are deleted ("" matches
  atoms with no alt conf). The occupancy of the deleted atoms is added to
  the remaining atoms of the same name.

  @param imol the model molecule index
  @param imodel the model number (or mmdb::MinInt4 for all models)
  @param chain_id the chain id
  @param resno the residue number
  @param inscode the insertion code
  @param altloc the alt conf of the atoms to delete */
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
/*! \brief delete hydrogen atoms in residue

  Hydrogen and deuterium atoms of the residue (in the first model) are
  deleted. The altloc argument is currently not used. */
void delete_residue_hydrogens(int imol, const char *chain_id, int resno, const char *inscode, const char *altloc);
/*! \brief delete atom in residue

  If the atom is the last atom of the residue, the whole residue is deleted.

  @param imol the model molecule index
  @param chain_id the chain id
  @param resno the residue number
  @param ins_code the insertion code
  @param at_name the atom name (PDB 4-character form, e.g. " CA ")
  @param altloc the alt conf (typically "") */
void delete_atom(int imol, const char *chain_id, int resno, const char *ins_code, const char *at_name, const char *altloc);
/*! \brief delete all atoms in residue that are not main chain or CB

  @param do_delete_dialog currently ignored
  @param imol the model molecule index
  @param chain_id the chain id
  @param resno the residue number
  @param ins_code the insertion code
*/
void delete_residue_sidechain(int imol, const char *chain_id, int resno, const char*ins_code,
			      short int do_delete_dialog);
/*! \brief delete all hydrogens in molecule,

   Synonym for \c delete_hydrogens().

   @return number of hydrogens deleted. */
int delete_hydrogen_atoms(int imol);

/*! \brief delete all hydrogens in molecule,

   Hydrogen and deuterium atoms are deleted.

   @return number of hydrogens deleted. */
int delete_hydrogens(int imol);

/*! \brief delete all waters in molecule,

   All residues named HOH are deleted. Note that no backup (undo point)
   is made.

   @return number of water atoms deleted (the same as the number of waters
           when they have no hydrogen atoms). */
int delete_waters(int imol);

/*! \brief (no longer does anything) */
void post_delete_item_dialog();



/* toggle callbacks */
/*! \brief set the delete-item mode to atom (GUI use) */
void set_delete_atom_mode();
/*! \brief set the delete-item mode to residue (GUI use) */
void set_delete_residue_mode();
/*! \brief set the delete-item mode to residue zone (GUI use) */
void set_delete_residue_zone_mode();
/*! \brief set the delete-item mode to residue hydrogens (GUI use) */
void set_delete_residue_hydrogens_mode();
/*! \brief set the delete-item mode to water (GUI use) */
void set_delete_water_mode();
/*! \brief set the delete-item mode to side chain (GUI use) */
void set_delete_sidechain_mode();
/*! \brief set the delete-item mode to side chain range (GUI use) */
void set_delete_sidechain_range_mode();
/*! \brief set the delete-item mode to chain (GUI use) */
void set_delete_chain_mode();
/*! \brief is the delete-item mode atom? (1 for yes, 0 for no) */
short int delete_item_mode_is_atom_p(); /* (predicate) a boolean */
/*! \brief is the delete-item mode residue? (1 for yes, 0 for no) */
short int delete_item_mode_is_residue_p(); /* predicate again */
/*! \brief is the delete-item mode water? (1 for yes, 0 for no) */
short int delete_item_mode_is_water_p();
/*! \brief is the delete-item mode side chain? (1 for yes, 0 for no) */
short int delete_item_mode_is_sidechain_p();
/*! \brief is the delete-item mode side chain range? (1 for yes, 0 for no) */
short int delete_item_mode_is_sidechain_range_p();
/*! \brief is the delete-item mode chain? (1 for yes, 0 for no) */
short int delete_item_mode_is_chain_p();
/*! \brief clear the pending delete-item atom, residue, residue zone and residue hydrogens modes (GUI use) */
void clear_pending_delete_item(); /* for when we cancel with picking an atom */


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  rotate/translate buttons                                */
/*  ----------------------------------------------------------------------- */
/* section Rotate/Translate Buttons */
/*  sets flag for atom selection clicks */
/*! \brief set up (state 1) or cancel (state 0) rotate/translate zone, ready for 2 atom picks (GUI use) */
void do_rot_trans_setup(short int state);
/*! \brief reset the stored previous values of the rotate/translate adjustments (GUI use) */
void rot_trans_reset_previous();
/*! \brief set whether rotate/translate zone rotates about the zone centre

  @param istate 1: rotate about the centre of the moving atoms; 0 (the
         default): rotate about the rotation-origin (clicked) atom */
void set_rotate_translate_zone_rotates_about_zone_centre(int istate);
/*! \brief set the rotate/translate object type

  @param rt_type one of ROT_TRANS_TYPE_RESIDUE (11), ROT_TRANS_TYPE_ZONE (12,
         the default), ROT_TRANS_TYPE_CHAIN (13), ROT_TRANS_TYPE_MOLECULE (14)
         or ROT_TRANS_TYPE_MULTI_RANGE (15) */
void set_rot_trans_object_type(short int rt_type); /* zone, chain, mol */
/*! \brief return the rotate/translate object type (see \c set_rot_trans_object_type()) */
int get_rot_trans_object_type();

/*  ----------------------------------------------------------------------- */
/*                  cis and trans info and conversion                       */
/*  ----------------------------------------------------------------------- */
/*! \brief set up (istate 1) or cancel (istate 0) cis/trans conversion, ready for an atom pick (GUI use) */
void do_cis_trans_conversion_setup(int istate);
/*! \brief cis/trans convert the peptide of the given residue

  Note that despite its name, the 4th argument is used as the insertion
  code (typically "").

  @param imol the model molecule index
  @param chain_id the chain id
  @param resno the residue number
  @param altconf the insertion code */
void cis_trans_convert(int imol, const char *chain_id, int resno, const char *altconf);

#ifdef __cplusplus	/* need this wrapper, else gmp.h problems in callback.c */
#ifdef USE_GUILE
/*! \brief return cis_peptide info for imol.

Return a SCM list object of (residue1 residue2 omega), where residue1 and
residue2 are residue specs and omega is in degrees */
SCM cis_peptides(int imol);
/*! \brief return twisted trans peptide info for imol

Twisted trans peptides have 90 < |omega| < 150 degrees.

Return a SCM list object of (residue1 residue2 omega), omega in degrees */
SCM twisted_trans_peptides(int imol);
#endif /* GUILE */
#ifdef USE_PYTHON
/*! \brief return cis_peptide info for imol.

Return a Python list object of [residue1, residue2, omega], where residue1
and residue2 are residue specs and omega is in degrees (an empty list
for an invalid molecule) */
PyObject *cis_peptides_py(int imol);
/*! \brief return twisted trans peptide info for imol

Twisted trans peptides have 90 < |omega| < 150 degrees.

Return a Python list object of [residue1, residue2, omega], where residue1
and residue2 are residue specs and omega is in degrees */
PyObject *twisted_trans_peptides_py(int imol);
#endif /* PYTHON */
#endif

/*! \brief cis-trans convert the active residue of the active atom in the
    intermediate atoms, and continue with the refinement

    The intermediate atom closest to the screen centre (within 2 Å) is used;
    a neighbouring peptide that is cis is preferred. The trans-peptide
    restraint is added or removed as appropriate.

    @return currently always 0 */
int cis_trans_convert_intermediate_atoms();


/*  ----------------------------------------------------------------------- */
/*                  db-main                                                 */
/*  ----------------------------------------------------------------------- */
/* section Mainchain Building Functions */

/*! \name Mainchain Building Functions */
/*! \{ */
/*! \brief set up (state 1) or cancel (state 0) CA-zone to mainchain conversion, ready for 2 atom picks (GUI use) */
void do_db_main(short int state);
/*! \brief CA -> mainchain conversion

Builds mainchain for the CA atoms of the given residue range by fitting
fragments from a database of reference structures. The result is put in
a new molecule (named "mainchain-" + direction).

direction is either "forwards" or "backwards" (anything other than
"backwards" is treated as forwards)

See also the function below.

  @param imol the model molecule index (containing the CA atoms)
  @param chain_id the chain id
  @param iresno_start the first residue number of the range
  @param iresno_end the last residue number of the range
  @param direction "forwards" or "backwards"

return the new molecule number, or -1 on failure */
int db_mainchain(int imol,
		 const char *chain_id,
		 int iresno_start,
		 int iresno_end,
		 const char *direction);

/*! \brief CA-Zone to Mainchain for a fragment based on the given residue.

Both directions are built. This is the modern interface.

The fragment is the set of residues connected to the given residue by
CA-CA distances of 4.5 Å or less. Two new molecules are created (one for
each direction).

  @param imol the model molecule index
  @param chain_id the chain id
  @param res_no the residue number of a residue in the fragment
  @return the index of the new "forwards" molecule, or -1 on failure
 */
int db_mainchains_fragment(int imol, const char *chain_id, int res_no);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  close molecule                                          */
/*  ----------------------------------------------------------------------- */
/* section Close Molecule Functions */
/*! \name Close Molecule Functions */
/*! \{ */

/*! \brief close the molecule

  Works for both model and map molecules.

  @param imol the molecule index */
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
  (ROTAMERSEARCHHIGHRES) (with rigid body fitting)

  The values are ROTAMERSEARCHAUTOMATIC = 0 (the default),
  ROTAMERSEARCHHIGHRES = 1 and ROTAMERSEARCHLOWRES = 2. In automatic mode,
  backrub rotamers are used when the map resolution is worse than 2.9 Å.
  Other values are ignored. */
void set_rotamer_search_mode(int mode);

/*! \brief get the mode of rotamer search (see \c set_rotamer_search_mode()) */
int rotamer_search_mode_state();

/*! \brief set up (state 1) or cancel (state 0) the rotamers dialog, ready for an atom pick (GUI use) */
void setup_rotamers(short int state);

/*  display the rotamer option and display the most likely in the graphics as a */
/*  moving_atoms_asc */
/*! \brief show the rotamer selection dialog for the residue of the given atom (GUI use)

  @param atom_index the atom index in the molecule's atom selection
  @param imol the model molecule index */
void do_rotamers(int atom_index, int imol);

/*! \brief show the rotamer selection dialog for the given residue (GUI use)

  @param imol the model molecule index
  @param chain_id the chain id
  @param resno the residue number
  @param ins_code the insertion code
  @param altconf the alt conf */
void show_rotamers_dialog(int imol, const char *chain_id, int resno, const char *ins_code, const char *altconf);

/*! \brief For Dunbrack rotamers, set the lowest probability to be
   considered.  Set as a percentage i.e. 1.00 is quite low.  For
   Richardson Rotamers, this has no effect. */
void set_rotamer_lowest_probability(float f);

/*! \brief set the flag for checking clashes in interactive rotamer fitting

  @param i 1 for on (the default), 0 for off */
void set_rotamer_check_clashes(int i);

/*! \brief auto fit by rotamer search.

   return the score, for some not very good reason.  clash_flag
   determines if we use clashes with other residues in the score for
   this rotamer (or not).  It would be cool to call this from a script
   that went residue by residue along a (newly-built) chain (now available).

   The search mode is set by \c set_rotamer_search_mode(). If imol_map
   is not a valid map, the rotamer is chosen by clash score only.

   @param imol_coords the model molecule index
   @param chain_id the chain id
   @param resno the residue number
   @param insertion_code the insertion code
   @param altloc the alt conf
   @param imol_map the map molecule index
   @param clash_flag 1 to include clashes in the score, 0 not to
   @param lowest_probability the lowest rotamer probability to consider
   @return the fit score, or -999.9 if imol_coords is not a valid model
           molecule */
float auto_fit_best_rotamer(int imol_coords,
                            const char *chain_id,
                            int resno,
			    const char *insertion_code,
			    const char *altloc,
			    int imol_map, int clash_flag, float lowest_probability);

/*! \brief auto-fit the rotamer for the active residue

  Uses the refinement map, clash checking and a lowest probability of 2.0.

  @return the fit score, or -1 if there is no active atom */
float auto_fit_rotamer_active_residue();

/*! \brief set the clash flag for rotamer search

   And this functions for [pre-setting] the variables for
   auto_fit_best_rotamer called interactively (using a graphics_info_t
   function). 0 off, 1 on.*/
void set_auto_fit_best_rotamer_clash_flag(int i); /*  */
/* currently stub function only */
/*! \brief return the rotamer probability of the given residue

  @return the rotamer probability, or 0 if it could not be determined
          (e.g. missing atoms, GLY or ALA, or residue not found) */
float rotamer_score(int imol, const char *chain_id, int res_no, const char *insertion_code,
		    const char *alt_conf);
/*! \brief set up (state 1) or cancel (state 0) rotamer auto-fit, ready for an atom pick (GUI use) */
void setup_auto_fit_rotamer(short int state);	/* called by the Auto Fit button call
				   back, set's in_auto_fit_define. */

/*! \brief return the number of rotamers for this residue - return -1
  on no residue found.*/
int n_rotamers(int imol, const char *chain_id, int resno, const char *ins_code);
/*! \brief set the residue specified to the rotamer number specifed.

   @param rotamer_number the rotamer index, from 0 to \c n_rotamers() - 1
   @return 1 if the atoms were moved, 0 otherwise
   @param imol the model molecule index
   @param chain_id the chain id
   @param resno the residue number
   @param ins_code the insertion code
   @param alt_conf the alt conf ("" for none)
*/
int set_residue_to_rotamer_number(int imol, const char *chain_id, int resno, const char *ins_code,
				  const char *alt_conf, int rotamer_number);

/*! \brief set the residue specified to the rotamer name specified

Note that the rotamer names are the Richardson rotamer names.

   @param imol the molecule index
   @param chain_id the chain-id
   @param resno the residue number
   @param ins_code the insertion code
   @param alt_conf the alt-conf
   @param rotamer_name the name of the rotamer
   @return value is 0 if atoms were not moved (e.g. because rotamer-name was not know)
*/
int set_residue_to_rotamer_name(int imol, const char *chain_id, int resno, const char *ins_code,
				const char *alt_conf, const char *rotamer_name);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the (Richardson) rotamer name of the given residue

  @return the rotamer name, or \#f if the residue is not found */
SCM get_rotamer_name_scm(int imol, const char *chain_id, int resno, const char *ins_code);
#endif
#ifdef USE_PYTHON
/*! \brief return the (Richardson) rotamer name of the given residue

  The residue's atoms with no alt conf are used.

  @return the rotamer name, or False if the residue is not found */
PyObject *get_rotamer_name_py(int imol, const char *chain_id, int resno, const char *ins_code);
#endif /* USE_GUILE */
#endif /* c++ */


/*! \brief fill all the residues of molecule number imol that have
   missing atoms.

Note that not all the atoms can be filled by this method. It is for filling
side-chain atom. The main chain atoms cannot be filled by this method.

To be used to remove the effects of chainsaw.

If the refinement map is set, the side chains are fitted by backrub
rotamer search and then the filled residues are refined (and the result
accepted). Otherwise the atoms are added without fitting and the
map-selection dialog is shown.

@param imol the molecule index

*/
void fill_partial_residues(int imol);

/*! \brief fill the missing side-chain atoms of the specified residue

Note that not all the atoms can be filled by this method. It is for filling
side-chain atom. The main chain atoms cannot be filled by this method.

If a mainchain, CA, C or N is missing, the best course is, if a map is available
for fitting, to delete the residue and add a residue back, then mutate it back
from ALA to whatever it used to be (if needed) and then auto-fit rotamer.

The refinement map must be set: the side chain is fitted by backrub
rotamer search and the residue is then refined (and the result accepted).
If no refinement map is set, nothing is done and the map-selection
dialog is shown.

@param imol the molecule index
@param chain_id the chain-id
@param resno the residue number
@param inscode the residue's inscode (typically "")

*/
void fill_partial_residue(int imol, const char *chain_id, int resno, const char* inscode);

/*! \brief Fill amino acid residues

do backrub rotamer search for residues, but don't do refinement

Needs the refinement map to be set; if it is not, nothing is done.
*/
void simple_fill_partial_residues(int imol);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the residues of molecule imol that have missing (non-hydrogen) atoms

  @return a list of (chain_id resno inscode) items, or \#f for an invalid molecule */
SCM missing_atom_info_scm(int imol);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return the residues of molecule imol that have missing (non-hydrogen) atoms

  Residues are checked against their dictionary entries.

  @return a list of [chain_id, resno, inscode] items, or False for an
          invalid molecule */
PyObject *missing_atom_info_py(int imol);
#endif /* USE_PYTHON */
#endif /* __cplusplus */



#ifdef __cplusplus	/* need this wrapper, else gmp.h problems in callback.c */
#ifdef USE_GUILE
/*! \brief Activate rotamer graph analysis for molecule number imol.

Return rotamer info - function used in testing.

@return a list of (chain_id resno inscode probability rotamer_name) items,
        or \#f if there are no results */
SCM rotamer_graphs(int imol);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief Activate rotamer graph analysis for molecule number imol.

Return rotamer info - function used in testing.

@return a list of [chain_id, resno, inscode, probability, rotamer_name]
        items, or False if there are no results */
PyObject *rotamer_graphs_py(int imol);
#endif /* USE_PYTHON */
#endif /* c++ */

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  180 degree flip                                         */
/*  ----------------------------------------------------------------------- */
/*! \name 180 Flip Side chain */
/*! \{ */

/*! \brief rotate 180 degrees around the last chi angle

  e.g. to flip the side chain amide of ASN or GLN, or the ring of HIS. */
void do_180_degree_side_chain_flip(int imol, const char* chain_id, int resno,
				   const char *inscode, const char *altconf);

/*! \brief set up (state 1) or cancel (state 0) a 180 degree side chain flip, ready for an atom pick (GUI use) */
void setup_180_degree_flip(short int state);

/*! \brief side-chain 180 flip the terminal chi angle of the residue of the active atom
  in the intermediate atoms, and continue with the refinement

  @return 1 if there are intermediate atoms (even if the active atom was
          not found in them), 0 otherwise */
int side_chain_flip_180_intermediate_atoms();


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  mutate                                                  */
/*  ----------------------------------------------------------------------- */
/* section Mutate Functions */
/*! \name Mutate Functions */
/*! \{ */

/* c-interface-build */
/*! \brief set up (state 1) or cancel (state 0) mutation, ready for an atom pick (GUI use) */
void setup_mutate(short int state);
/*! \brief Mutate then fit to map

 that we have a map define is checked first: if the refinement map is not
 set, the map-selection dialog is shown instead (GUI use) */
void setup_mutate_auto_fit(short int state);

/*! \brief mutate the residue of the active atom (GUI use)

  @param type the 3-letter code of the target residue type
  @param is_stub_flag 1 to add only a stub (mainchain + CB), 0 for the full side chain */
void do_mutation(const char *type, short int is_stub_flag);

/*! \brief display a dialog that allows the choice of residue type to which to mutate
 */
void mutate_active_residue();

/* auto-mutate stuff */
/*! \brief check that residue numbers increase along the chain

  @param chain_id the chain id
  @param imol the model molecule index
  @return 1 if the residue numbers in the chain (in the first model) are
          strictly increasing, 0 otherwise (or if the chain is not found) */
short int progressive_residues_in_chain_check(const char *chain_id, int imol);

/*! \brief mutate a given residue

target_res_type is a three-letter-code.

A backup is made.

Return 1 on a good mutate, 0 on failure (or invalid molecule), -1 if the
residue was not found. */
int mutate(int imol, const char *chain_id, int ires, const char *inscode,  const char *target_res_type);

/*! \brief mutate a base. return success status, 1 for a good mutate.

  @param res_type the new base type: "A", "C", "G", "T" or "U" (converted
         to the DNA names if the residue is DNA), or "DA", "DC", "DG", "DT"
  @param imol the model molecule index
  @param chain_id the chain id
  @param res_no the residue number
  @param ins_code the insertion code
*/
int mutate_base(int imol, const char *chain_id, int res_no, const char *ins_code, const char *res_type);

/*! \brief push the residues along a bit

e.g. if nudge_by is 1, then the sidechain of residue 20 is moved up
onto what is currently residue 21.  The mainchain numbering and atoms is not changed.

The residue types are shifted within the range only: the first nudge_by
residues of the range keep their types. res_no_range_start must be less
than res_no_range_end.

@param nudge_residue_numbers_also if non-zero, the residue numbers in the
       range are also decreased by nudge_by

@return 0 for failure to nudge (because not all the residues were in the range)
        and 1 for success.
@param imol the model molecule index
@param chain_id the chain id
@param res_no_range_start the first residue number of the range
@param res_no_range_end the last residue number of the range
@param nudge_by the number of residues by which the residue types are shifted
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

   This function doesn't make backups, but `mutate()` does.
   Hence `mutate()` is for use as a "one-by-one" type and the following
   2 by wrappers that mutate either a residue range or a whole chain.

   If the residue is already of the target type, nothing is changed (and
   1 is returned).

   @param ires_ser is the serial number of the residue (the 0-based index
          of the residue in the chain, in the first model), not the seqnum
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

/*!  \brief mutate a single residue (given by seqnum) to the type given by a single-letter code

  ires is the seqnum of the residue (conventional)

  @param imol the model molecule index
  @param chain_id the chain id
  @param ires the residue number
  @param inscode the insertion code
  @param target_res_type the single-letter code of the target residue type
  @return 1 on success, 0 on failure to mutate, -1 if imol is invalid or
          the residue was not found */
int mutate_single_residue_by_seqno(int imol, const char *chain_id, int ires, const char *inscode,
				   char target_res_type);

/*! \brief mutate and auto-fit

(Move this and the above function into cc-interface.hh one day)

Mutates the residues start_res_no to stop_res_no (with blank insertion
codes) to the given sequence and then auto-fits the rotamer of each
against the refinement map (by clash score only if that map is not set).
Nothing is done if the length of the sequence does not match the
residue range.

  @param imol the model molecule index
  @param chain_id the chain id
  @param start_res_no the first residue number of the range
  @param stop_res_no the last residue number of the range
  @param sequence the target sequence as single-letter codes, one per residue
  @return currently always 0
 */
int mutate_and_autofit_residue_range(int imol, const char *chain_id, int start_res_no, int stop_res_no,
                                     const char *sequence);

/* an internal function - not useful for scripting: */

/*! \brief mutate the base of the previously picked nucleotide (internal GUI use)

  @param type the base type, e.g. "A", "C", "G", "T" or "U" */
void do_base_mutation(const char *type);

/*! \brief set a flag saying that the residue chosen by mutate or
  auto-fit mutate should only be added as a stub (mainchain + CB) */
void set_residue_type_chooser_stub_state(short int istat);

/*! \brief mutate (and auto-fit) the residue of the active atom using the residue
  type typed into the residue-type chooser entry (GUI use)

  @param entry_text text whose first character is used as a single-letter
         amino acid code
  @param stub_mode 1 to add only a stub (mainchain + CB), 0 for the full side chain */
void handle_residue_type_chooser_entry_chose_type(const char *entry_text, short int stub_mode);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  alternate conformation                                  */
/*  ----------------------------------------------------------------------- */
/* section Alternative Conformation */
/*! \name Alternative Conformation */
/* c-interface-build function */
/*! \{ */

/*! \brief return the alt conf split type (see \c set_add_alt_conf_split_type_number()) */
short int alt_conf_split_type_number();
/*! \brief set the alt conf split type

  @param i 0: split at CA (partial split, mainchain not split); 1: split
         the whole residue (the default); 2: split a residue range */
void set_add_alt_conf_split_type_number(short int i);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief add an alternative conformer to a residue.  Add it in
  conformation rotamer number rotamer_number.

  Note that rotamer_number is currently ignored: the residue is simply
  split (according to the split type, see
  \c set_add_alt_conf_split_type_number()) and the new atoms are given
  the occupancy set by \c set_add_alt_conf_new_atoms_occupancy().

  @param imol the model molecule index
  @param chain_id the chain id
  @param res_no the residue number
  @param ins_code the insertion code
  @param alt_conf the alt conf of the atoms to be split (typically "")
  @param rotamer_number currently ignored
  @return the new alt conf on success, scheme false on fail
*/
SCM add_alt_conf_scm(int imol, const char *chain_id, int res_no, const char *ins_code,
		     const char *alt_conf, int rotamer_number);
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief add an alternative conformer to a residue.  Add it in
  conformation rotamer number rotamer_number.

  Note that rotamer_number is currently ignored: the residue is simply
  split (according to the split type, see
  \c set_add_alt_conf_split_type_number()) and the new atoms are given
  the occupancy set by \c set_add_alt_conf_new_atoms_occupancy().

  @param imol the model molecule index
  @param chain_id the chain id
  @param res_no the residue number
  @param ins_code the insertion code
  @param alt_conf the alt conf of the atoms to be split (typically "")
  @param rotamer_number currently ignored
  @return the new alt conf (e.g. "B") on success, python False on fail */
PyObject *add_alt_conf_py(int imol, const char*chain_id, int res_no, const char *ins_code,
		     const char *alt_conf, int rotamer_number);
#endif	/* USE_PYTHON */
#endif /* __cplusplus */

/*! \brief forget the add-alt-conf dialog (internal GUI use) */
void unset_add_alt_conf_dialog(); /* set the static dialog holder in
				     graphics info to NULL */
/*! \brief cancel a pending add-alt-conf atom pick (GUI use) */
void unset_add_alt_conf_define(); /* turn off pending atom pick */
/*! \brief show the add-alt-conf dialog (GUI use) */
void altconf();			/* temporary debugging interface. */
/*! \brief set the occupancy of the new atoms made by adding an alt conf

  @param f the occupancy (default 0.5) */
void set_add_alt_conf_new_atoms_occupancy(float f); /* default 0.5 */
/*! \brief return the occupancy of the new atoms made by adding an alt conf */
float get_add_alt_conf_new_atoms_occupancy();
/*! \brief set whether the new alt conf atoms are shown as intermediate atoms

  When off (and the residue has the atoms needed), the new alt conf atoms
  are added directly to the model; when on, they are shown as intermediate
  atoms.

  @param i 1 for on, 0 for off (default 0) */
void set_show_alt_conf_intermediate_atoms(int i);
/*! \brief return the state of the show-alt-conf-intermediate-atoms flag (1 on, 0 off) */
int  show_alt_conf_intermediate_atoms_state();
/*! \brief set the occupancy of the atoms in the residue range to 0.0 */
void zero_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2);
/*! \brief set the occupancy of the atoms in the residue range to 1.0 */
void fill_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2);
/*! \brief set the occupancy of the atoms in the residue range

  @param imol the model molecule index
  @param chain_id the chain id
  @param ires1 the first residue number of the range
  @param ires2 the last residue number of the range
  @param occ the new occupancy */
void set_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2, float occ);
/*! \brief set the B-factor of the atoms in the residue range

  @param imol the model molecule index
  @param chain_id the chain id
  @param ires1 the first residue number of the range
  @param ires2 the last residue number of the range
  @param bval the new B-factor (in Å^2) */
void set_b_factor_residue_range(int imol, const char *chain_id, int ires1, int ires2, float bval);
/*! \brief reset the B-factor of the atoms in the residue range to the default B-factor for new atoms

  (see \c set_default_temperature_factor_for_new_atoms()) */
void reset_b_factor_residue_range(int imol, const char *chain_id, int ires1, int ires2);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  pointer atoms                                           */
/*  ----------------------------------------------------------------------- */
/* section Pointer Atom Functions */
/*! \name Pointer Atom Functions */
/* c-interface-build */
/*! \{ */

/*! \brief place an atom at the pointer (the screen centre)

  If the pointer atom is a dummy (see \c set_pointer_atom_is_dummy()), a
  dummy atom is added to the pointer-atom molecule; otherwise the atom-type
  dialog is shown. */
void place_atom_at_pointer();
/* which calls the following gui function (if using non dummies) */
/*! \brief show the pointer-atom type dialog (GUI use) */
void place_atom_at_pointer_by_window();
/*! \brief place an atom of the given type at the pointer (the screen centre)

  The atom is added to the molecule of the active atom (which must be
  displayed).

  @param type the atom type, e.g. "Water", an element/ion name, "SO4" or "PO4" */
void place_typed_atom_at_pointer(const char *type);

/*! \brief set whether pointer atoms are dummy atoms

  @param i 1 for dummy atoms (placed directly with no type dialog), 0 (the
         default) to choose the atom type with a dialog */
void set_pointer_atom_is_dummy(int i);
/*! \brief print the coordinates of the pointer (the screen centre) to the console */
void display_where_is_pointer(); /* print the coordinates of the
				    pointer to the console */
/*! \brief Return the current pointer atom molecule, create a pointer
  atom molecule if necessary (i.e. when the user has not set it).

  An existing molecule named "Pointer Atoms" is used if there is one,
  otherwise a new (empty) molecule with that name is created. */
int create_pointer_atom_molecule_maybe();
/*! \brief Return the current pointer atom molecule

  @return the molecule set by \c set_pointer_atom_molecule(), or -1 if
          it has not been set */
int pointer_atom_molecule();
/*! \brief set the molecule to which pointer atoms are added

  @param imol a valid model molecule index (otherwise ignored) */
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
/*! \brief draw the baton or not

  Turning the baton on needs a skeletonized map.

  @param i 1 to try to turn the baton on, 0 to turn it off
  @return the resulting baton-drawing state (1 drawn, 0 not) */
int try_set_draw_baton(short int i); /* draw the baton or not */
/*! \brief accept the baton tip position

  a prime candidate for a key binding */
void accept_baton_position();	/* put an atom at the tip */
/*! \brief move the baton tip position - another prime candidate for a key binding

  Cycles to the next candidate CA position (wrapping round to the first). */
void baton_tip_try_another();
/*! \brief move the baton tip to the previous position*/
void baton_tip_previous();
/*! \brief shorten the baton length (by a factor of 0.952) */
void shorten_baton();
/*! \brief lengthen the baton (by a factor of 1.05) */
void lengthen_baton();
/*! \brief delete the most recently build CA position */
void baton_build_delete_last_residue();
/*! \brief set the parameters for the start of a new baton-built fragment. direction can either
     be "forwards" or "backwards"

  @param istart_resno the residue number of the first baton-built residue
  @param chain_id the chain id for the new residues
  @param direction "forwards" or "backwards" (anything else is treated
         as an unknown direction) */
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
   fragments in a chain.  (Note that, currently, only numbering gaps are
   used - no CA-CA distance check is made.) Throw away all atoms in
   fragment other than CAs.

   @param imol the model molecule index
   @param chain_id the chain id
   @param resno the residue number of a residue in the fragment */
void reverse_direction_of_fragment(int imol, const char *chain_id, int resno);
/*! \brief set up (1) or cancel (0) reverse-direction-of-fragment, ready for an atom pick (GUI use) */
void setup_reverse_direction(short int i);


/*  ----------------------------------------------------------------------- */
/*                  terminal OXT atom                                       */
/*  ----------------------------------------------------------------------- */
/* section Terminal OXT Atom */
/*! \name Terminal OXT Atom */
/* c-interface-build */
/*! \{ */
/*! \brief add an OXT atom to the given residue

  The residue needs N, CA, C and O atoms and must not already have an OXT.

  @param imol the model molecule index
  @param chain_id the chain id
  @param reso the residue number
  @param insertion_code the insertion code (typically "")
  @return 1 on success, 0 on failure, -1 if imol is not a valid model molecule */
short int add_OXT_to_residue(int imol, const char *chain_id, int reso, const char *insertion_code);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  crosshairs                                              */
/*  ----------------------------------------------------------------------- */
/* section Crosshairs Interface */
/*! \name Crosshairs  Interface */
/*! \{ */
/*! \brief draw the distance crosshairs, 0 for off, 1 for on.

  The crosshair ticks are at 1.54 Å (C-C bond), 2.7 Å (H-bond) and
  3.8 Å (CA-CA). Turning them on prints that key to the console and
  redraws. Default: off. */
void set_draw_crosshairs(short int i);
/*! \brief return the crosshairs drawing state

  so that we display the crosshairs with the radiobuttons in the
  right state.

  @return 1 if the crosshairs are drawn, 0 if not */
short int draw_crosshairs_state();
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Edit Chi Angles                                         */
/*  ----------------------------------------------------------------------- */
/* section Edit Chi Angles */
/*! \name  Edit Chi Angles */
/*! \{ */
/* c-interface-build functions */
/*! \brief start (or cancel) the "edit chi angles" atom-pick mode

  With state 1, the next atom click selects the residue whose chi
  angles are to be edited (the clicked atom determines which side of
  the torsion moves). With state 0, the pick mode is cancelled. */
void setup_edit_chi_angles(short int state);

/*! \brief rotate the currently selected chi angle

  Only acts while in chi-angle (or general torsion) editing mode. The
  chi chosen with \c set_graphics_edit_current_chi() is changed by
  20 x \c am degrees. If no chi has been chosen, the residue is
  flashed to indicate that one needs to be selected. */
void rotate_chi(float am);

/*! \brief show torsions that rotate hydrogens in the torsion angle
  manipulation dialog.  Note that this may be needed if, in the
  dictionary cif file torsion which have as a 4th atom both a hydrogen
  and a heavier atom bonding to the 3rd atom, but list the 4th atom as
  a hydrogen (not a heavier atom).

  @param state 1 for on, 0 for off (default off) */
void set_find_hydrogen_torsions(short int state);
/*! \brief select which chi angle is edited by mouse motion (a button callback)

  @param ichi the chi angle number (starting at 1); 0 turns off chi editing mode */
void set_graphics_edit_current_chi(int ichi); /* button callback */
/*! \brief stop the keyboard (1, 2, 3 ...) from putting the intermediate atoms
  into rotate-chi mode */
void unset_moving_atom_move_chis();
/*! \brief allow the keyboard (1, 2, 3 ...) to put the intermediate atoms
  into rotate-chi mode */
void set_moving_atom_move_chis();

/*! \brief display the edit chi angles gui for the given residue

 The edit uses the first atom of the residue; altconf is used as the
 alt conf of the atoms to be edited.

 @return a status of 0 if it failed to find the residue (or imol is
 not a valid model molecule), 1 if it worked. */
int edit_chi_angles(int imol, const char *chain_id, int resno,
		     const char *ins_code, const char *altconf);

/*! \brief set the flag to show the "flash" bond (the rotatable bond) of
  the chi angle being edited

  @param imode 1 for on, 0 for off (default off)
  @return 0 always */
int set_show_chi_angle_bond(int imode);

/* a callback from the callbacks.c, setting the state of
   graphics_info_t::edit_chi_angles_reverse_fragment flag */
/*! \brief set the flag to reverse the moving fragment (i.e. which side of the
  bond moves) when editing chi angles

  @param istate 1 for reversed, 0 for normal (default 0) */
void set_edit_chi_angles_reverse_fragment_state(short int istate);

/* No need for this to be exported to scripting */
/*! \brief beloved torsion general at last makes an entrance onto the
  Coot scene...

  With state 1, start the general-torsion atom-pick mode: the user then
  clicks 4 atoms that define the torsion to be edited (the order of the
  clicked atoms affects which part moves). With state 0, cancel the
  pick mode. */
void setup_torsion_general(short int state);
/* No need for this to be exported to scripting */
/*! \brief toggle which end of the fragment moves in general torsion editing */
void toggle_torsion_general_reverse();

/*! \brief start (state 1) or cancel (state 0) the atom-pick mode to choose a
  residue for adding partial alt confs

  After the pick, the edit-chi-angles dialog is shown for the residue in
  "partial alt locs" mode. */
void setup_residue_partial_alt_locs(short int state);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Backrub                                                 */
/*  ----------------------------------------------------------------------- */
/*! \name Backrubbing function */
/*! \{ */
/*! \brief Do a back-rub rotamer search (with autoaccept).

The model is changed in place. Needs the refinement map to be set
(\c set_imol_refinement_map()).

@param imol the model molecule index
@param chain_id the chain id
@param res_no the residue number
@param ins_code the insertion code
@param alt_conf the alt conf of the residue
@return the success status, 0 for fail, 1 for successful fit.  */
int backrub_rotamer(int imol, const char *chain_id, int res_no,
		    const char *ins_code, const char *alt_conf);

/*! \brief apply rotamer backrub to the active atom of the intermediate atoms

The residue used is that of the intermediate (moving) atom closest to
the screen centre (within 2 Å); it needs both a preceding and a
following residue in the intermediate atoms. Needs the refinement map
to be set. Refinement of the intermediate atoms is then restarted.

@return 1 on success, 0 on failure */
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
        where there are no atoms atoms.

        Waters are not used for masking unless the find-ligand "mask
        waters" flag is set. The mask radius around each atom is 2.0 Å
        unless set with \c set_map_mask_atom_radius().

        @param map_mol_no the map molecule index
        @param coord_mol_no the model molecule index
        @param invert_flag 0 to remove density at the atoms, 1 to remove density
               away from the atoms
        @return the index of the new map molecule, or -1 on failure */
int mask_map_by_molecule(int map_mol_no, int coord_mol_no, short int invert_flag);

/*! \brief mask map by atom selection

  Create a new map from map_mol_no, masked by the atoms of
  coords_mol_no that match the given mmdb atom selection string. The
  new map is contoured at 0.99 x the contour level of the input map.

  @param map_mol_no the map molecule index
  @param coords_mol_no the model molecule index
  @param mmdb_atom_selection the mmdb-format atom selection string, e.g. "//A/1-10"
  @param invert_flag 0 to remove density at the selected atoms, 1 to keep
         only the density at the selected atoms
  @return the index of the new map molecule, or -1 on failure */
int mask_map_by_atom_selection(int map_mol_no, int coords_mol_no, const char *mmdb_atom_selection, short int invert_flag);

/*! \brief make chain masked maps

   Create one new map per chain in (the first model of) imol, each
   containing only the density of imol_map within 3.3 Å of the atoms of
   that chain. The new maps are named "Masked Map for Chain X".

   @param imol the model molecule index
   @param imol_map the map molecule index
   @return 0 (always - the new map molecule indices are not returned)
 */
int make_masked_maps_split_by_chain(int imol, int imol_map);

/*! \brief set the atom radius (in Å) for map masking

  Used by \c mask_map_by_molecule() and \c mask_map_by_atom_selection()
  when the value is positive. */
void set_map_mask_atom_radius(float rad);
/*! \brief get the atom radius for map masking

  @return the radius in Å, or -99 if it has not been set (in which case
  the default of 2.0 Å is used) */
float map_mask_atom_radius();

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  check waters interface                                  */
/*  ----------------------------------------------------------------------- */
/* section Check Waters Interface */
/*! \name check Waters Interface */
/* interactive check by b-factor, density level etc. */
/*! \{ */
/*! \brief set the B-factor limit for the check-waters dialog (default 80.0) */
void set_check_waters_b_factor_limit(float f);
/*! \brief set the map level limit (in map rmsd) for the check-waters dialog (default 1.0) */
void set_check_waters_map_sigma_limit(float f);
/*! \brief set the minimum distance limit (in Å) for the check-waters dialog (default 2.3) */
void set_check_waters_min_dist_limit(float f);
/*! \brief set the maximum distance limit (in Å) for the check-waters dialog (default 3.5) */
void set_check_waters_max_dist_limit(float f);


/*! \brief Delete waters that are fail to meet the given criteria.

  Uses the refinement map (\c set_imol_refinement_map()) for the
  density test.

  @param imol the model molecule index
  @param b_factor_lim waters with a B-factor above this are flagged
  @param map_sigma_lim waters in density below this level (in map rmsd) are flagged
  @param min_dist waters closer than this (in Å) to their nearest (non-hydrogen) atom are flagged
  @param max_dist waters further than this (in Å) from their nearest (non-hydrogen) atom are flagged
  @param part_occ_contact_flag if non-zero, the distance tests are skipped
  @param zero_occ_flag if non-zero, the distance tests are skipped for waters with zero occupancy
  @param logical_operator_and_or_flag 0 for AND (a water must fail all the
         criteria), 1 for OR (a water failing any criterion is deleted).
         In OR mode, a negative b_factor_lim, min_dist or max_dist disables
         that test, as does a map_sigma_lim below -50. */
void delete_checked_waters_baddies(int imol, float b_factor_lim,
				   float map_sigma_lim,
				   float min_dist, float max_dist,
				   short int part_occ_contact_flag,
				   short int zero_occ_flag,
				   short int logical_operator_and_or_flag);

/* difference map variance check  */
/*! \brief check waters using the difference map

  For each water oxygen, the sum of squares of the difference map
  density within 1.5 Å is calculated; waters whose value is more than
  the check-waters difference-map sigma level (see
  \c set_check_waters_by_difference_map_sigma_level()) standard
  deviations above the mean are flagged.

  @param imol_waters the model molecule index
  @param imol_diff_map the map molecule index, which must be a difference map
  @param interactive_flag if non-zero, show a dialog listing the suspicious waters */
void check_waters_by_difference_map(int imol_waters, int imol_diff_map,
				    int interactive_flag);
/* results widget are in graphics-info.cc  */
/* Let's give access to the sigma level (default 4) */
/*! \brief return the sigma level used by \c check_waters_by_difference_map() (default 3.5) */
float check_waters_by_difference_map_sigma_level_state();
/*! \brief set the sigma level used by \c check_waters_by_difference_map() (default 3.5) */
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

/*! \brief print the metal coordination distances of molecule imol to the console

  Only does something if the molecule has symmetry (a cell and space
  group).

  @return scheme false (always) */
SCM metal_coordination_scm(int imol, float dist_max);
#endif
#ifdef USE_PYTHON
/*! \brief return a list first of metals, second of waters that are
  coordinated with at least coordination_number of other atoms at
  distances less than or equal to dist_max.

  The metals list is of [atom_spec, element] pairs (metals are
  identified amongst the waters using a 4.0 Å contact distance), the
  waters list is of [central_atom_spec, [neighbour_atom_specs]] pairs.

  @param imol the model molecule index
  @param coordination_number the minimum number of contacts
  @param dist_max the maximum contact distance in Å
  @return [metals, waters], or False if imol is not a valid molecule.
  Note that False is also returned when no highly coordinated waters
  are found. */
PyObject *highly_coordinated_waters_py(int imol, int coordination_number, float dist_max);
/*! \brief print the metal coordination distances of molecule imol to the console

  Only does something if the molecule has symmetry (a cell and space
  group).

  @param imol the model molecule index
  @param dist_max the maximum contact distance in Å
  @return False (always) */
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
/*! \brief clear the list of LSQ matches */
void clear_lsq_matches();
/*! \brief add a residue-range match to the list of LSQ matches

  match_type: 0 all atoms (main-chain only if the residue types
  differ), 1 main-chain, 2 CA (P for nucleotides), 3 N, CA and C,
  4 N, CA, CB and C. */
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
/*! \brief add a single atom pair to the list of LSQ matches */
void add_lsq_atom_pair_scm(SCM atom_spec_ref, SCM atom_spec_moving);
#endif
#ifdef USE_PYTHON
/*! \brief add a single atom pair to the list of LSQ matches

  @param atom_spec_ref the atom spec of the atom in the reference molecule
  @param atom_spec_moving the atom spec of the atom in the moving molecule */
void add_lsq_atom_pair_py(PyObject *atom_spec_ref, PyObject *atom_spec_moving);
#endif
#endif /* __cplusplus */

#ifdef __cplusplus	/* need this wrapper, else gmp.h problems in callback.c */
#ifdef USE_GUILE
/*! \brief apply the LSQ matches

@return an rtop pair (proper list) on good match, else false */
SCM apply_lsq_matches(int imol_reference, int imol_moving);
/*! \brief get the LSQ matrix from the current matches, without moving imol_moving

@return an rtop pair (proper list) on good match, else false */
SCM get_lsq_matrix_scm(int imol_reference, int imol_moving);
#endif
#ifdef USE_PYTHON
/*! \brief apply the LSQ matches

  Calculate the LSQ transformation from the current list of matches and
  apply it to imol_moving (which also gets the cell and symmetry of
  imol_reference).

  @return an rtop pair [[m11, m12, ... m33], [t1, t2, t3]] on good match, else False */
PyObject *apply_lsq_matches_py(int imol_reference, int imol_moving);
/*! \brief get the LSQ matrix from the current matches, without moving imol_moving

  @return an rtop pair [[m11, m12, ... m33], [t1, t2, t3]] on good match, else False */
PyObject *get_lsq_matrix_py(int imol_reference, int imol_moving);
#endif /* PYTHON */
#endif /* __cplusplus */

/* poor old python programmers... */
/*! \brief apply the LSQ matches (simple interface)

  As \c apply_lsq_matches_py() but just returns the status.

  @return 1 on a good match (imol_moving is moved), 0 on failure */
int apply_lsq_matches_simple(int imol_reference, int imol_moving);

/* section Least-Squares plane interface */
/*! \brief start (state 1) or stop (state 0) the atom-pick mode that reports
  the deviation of picked atoms from the LSQ plane */
void setup_lsq_deviation(int state);
/*! \brief start (state 1) or stop (state 0) the atom-pick mode that adds
  picked atoms to the LSQ plane definition */
void setup_lsq_plane_define(int state);
/*! \brief clear the LSQ plane dialog and its atoms (a callback from the
  destroy of the widget) */
void unset_lsq_plane_dialog(); /* callback from destroy of widget */
/*! \brief remove the most recently added atom of the LSQ plane

  Does nothing if there is only one atom (or none). */
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
  molecule if they are below the given map (absolute) level.

  @param imol_coords the model molecule index
  @param imol_map the map molecule index
  @param map_level the absolute map level (not rmsd)
  @param delete_or_zero_occ_flag 0 to delete the atoms, 1 to set their occupancy to zero */
void trim_molecule_by_map(int imol_coords, int imol_map,
			  float map_level, int delete_or_zero_occ_flag);

/*! \brief trim the molecule by the value in the B-factor column.

If an atom in a residue has a "B-factor" above (or below, if keep_higher is true) limit, then the whole residue is deleted */
void trim_molecule_by_b_factor(int imol, float limit, short int keep_higher);

/*! \brief convert the value in the B-factor column (typically pLDDT for AlphaFold models) to a temperature factor

  The new B-factor is 2 x (100 - pLDDT), with a minimum of 2.0. */
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
/*! \brief create a povray (.pov) file for the current view */
void povray(const char *filename);
/*! \brief create a RenderMan (RIB) file for the current view */
void renderman(const char *rib_filename);
/* a wrapper for the (scheme) function that makes the image, callable
   from callbacks.c  */
/*! \brief write filename.r3d for the current view, then run Raster3D's
  render on it (via the scripting layer) to make the image filename
  and display it */
void make_image_raster3d(const char *filename);
/*! \brief write filename.pov for the current view, then run povray on it
  (via the scripting layer) to make an image and display it */
void make_image_povray(const char *filename);
#ifdef USE_PYTHON
/*! \brief Python-specific version of \c make_image_raster3d() */
void make_image_raster3d_py(const char *filename);
/*! \brief Python-specific version of \c make_image_povray() */
void make_image_povray_py(const char *filename);
#endif /* USE_PYTHON */

/*! \brief set the bond thickness for the Raster3D representation (default 0.18) */
void set_raster3d_bond_thickness(float f);
/*! \brief set the atom radius for the Raster3D representation (default 0.25) */
void set_raster3d_atom_radius(float f);
/*! \brief set the density line thickness for the Raster3D representation (default 0.015) */
void set_raster3d_density_thickness(float f);
/*! \brief set the flag to show atoms for the Raster3D representation

  @param istate 1 for on, 0 for off (default 1) */
void set_renderer_show_atoms(int istate);
/*! \brief set the bone (skeleton) thickness for the Raster3D representation (default 0.05) */
void set_raster3d_bone_thickness(float f);
/*! \brief turn off shadows for raster3d output - give argument 0 to turn off

  @param state 1 for shadows on (the default), 0 for off */
void set_raster3d_shadows_enabled(int state);
/*! \brief set the flag to show waters as spheres for the Raster3D
representation. 1 show as spheres, 0 the usual stars. */
void set_raster3d_water_sphere(int istate);
/*! \brief set the font size (as a string) for raster3d (default "4") */
void set_raster3d_font_size(const char *size_in);
/*! \brief run raster3d and display the resulting image.

  Calls the scripting function render-image (Scheme, if available) or
  render_image() (Python), which writes coot.r3d, renders it and opens
  the resulting image. */
void raster_screen_shot(); /* run raster3d or povray and guile */
                           /* script to render and display image */
#ifdef USE_PYTHON
/*! \brief run raster3d and display the resulting image (using the Python render_image() function) */
void raster_screen_shot_py(); /* run raster3d or povray and python */
#endif
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  citation notice                                         */
/*  ----------------------------------------------------------------------- */
/*! \brief turn off the citation notice

  (The notice is off by default and this flag is no longer used
  elsewhere.) */
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
 * @param move_imol2_flag 1 to make a superposed copy of imol2 as a new
 *        molecule, 0 to move imol2 itself

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
 * @param move_imol2_copy_flag flag to control if imol2 is
 *        copied (1) or moved (0)

Superpose the given chains of imol2 onto imol1.  imol1 is reference,
we can either move imol2 or copy it to generate a new molecule
depending on the vaule of move_imol2_copy_flag (1 for copy 0 for move). */
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
   the vaule of move_imol2_copy_flag (1 for copy 0 for move).

   @return the index of the superposed molecule - which could either be a
   new molecule (if move_imol2_copy_flag was 1) or the imol2 or -1 (signifying
   failure to do the SSM superposition).
 *
 * @param imol1 the reference model index
 * @param imol2 the index of the superposed molecule
 * @param mmdb_atom_sel_str_1 the mmdb-format atom selection for imol1
 * @param mmdb_atom_sel_str_2 the mmdb-format atom selection for imol2
 * @param move_imol2_copy_flag flag to control if imol2 is
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
/*! \brief set drawing state of NCS ghosts for molecule number imol

  @param imol the model molecule index
  @param istate 1 for on, 0 for off */
void set_draw_ncs_ghosts(int imol, int istate);
/*! \brief return the drawing state of NCS ghosts for molecule number
  imol.

  @return 1 if the ghosts are drawn, 0 if not, -1 if imol is not a
  valid model molecule.  */
int draw_ncs_ghosts_state(int imol);

/*! \brief set bond thickness of NCS ghosts for molecule number imol   */
void set_ncs_ghost_bond_thickness(int imol, float f);
/*! \brief update ghosts for molecule number imol */
void ncs_update_ghosts(int imol); /* update ghosts */
/*! \brief make NCS map

  Make a map for each NCS ghost of imol_model by transforming
  imol_map, named "Map <imol_map> <ghost name>". NCS-averaged maps are
  also made (as they are by default). The NCS operators are calculated
  first if needed.

  @param imol_model the model molecule index (with NCS)
  @param imol_map the map molecule index
  @param overwrite_maps_of_same_name_flag if 1, reuse an existing map
         molecule with the same name rather than creating a new one
  @return the number of maps made */
int make_dynamically_transformed_ncs_maps(int imol_model, int imol_map,
					  int overwrite_maps_of_same_name_flag);
/*! \brief calculate the NCS ghost operators for imol if the molecule has
  NCS and they have not yet been calculated */
void make_ncs_ghosts_maybe(int imol);
/*! \brief Add NCS matrix

  Add an NCS ghost for the atoms of chain this_chain_id, using the
  given operator (rotation matrix m and translation t in Å) that maps
  this_chain_id onto target_chain_id. */
void add_ncs_matrix(int imol, const char *this_chain_id, const char *target_chain_id,
		    float m11, float m12, float m13,
		    float m21, float m22, float m23,
		    float m31, float m32, float m33,
		    float t1,  float t2,  float t3);

/*! \brief remove all the NCS ghosts (and their matrices) from molecule imol */
void clear_ncs_ghost_matrices(int imol);

/*! \brief add an NCS matrix for strict NCS molecule representation

for CNS strict NCS usage: expand like normal symmetry does

@return 1 on success, 0 if imol is not a valid model molecule */
int add_strict_ncs_matrix(int imol,
			  const char *this_chain_id,
			  const char *target_chain_id,
			  float m11, float m12, float m13,
			  float m21, float m22, float m23,
			  float m31, float m32, float m33,
			  float t1,  float t2,  float t3);
/*! \brief add strict NCS matrices from the MTRIX records of molecule imol's
  own coordinates file

  @return 0 (always) */
int add_strict_ncs_from_mtrix_from_self_file(int imol);

/*! \brief return the state of NCS ghost molecules for molecule number imol

  @return the show-strict-NCS flag, or 0 if imol is not a valid model molecule */
int show_strict_ncs_state(int imol);
/*! \brief set display state of NCS ghost molecules for molecule number imol   */
void set_show_strict_ncs(int imol, int state);
/*! \brief At what level of homology should we say that we can't see homology
   for NCS calculation? (default 0.7) */
void set_ncs_homology_level(float flev);

/* for a single copy */
/*! \brief Copy single NCS chain

  Replace chain to_chain with an NCS-transformed copy of from_chain.
  This needs NCS ghosts in which from_chain is the master (target) of
  to_chain. */
void copy_chain(int imol, const char *from_chain, const char *to_chain);
/* do multiple copies */
/*! \brief Copy chain from master to all related NCS chains

  @param imol the model molecule index
  @param chain_id the master chain id */
void copy_from_ncs_master_to_others(int imol, const char *chain_id);
/*! \brief Copy residue range to all related NCS chains.

  If the
  target residues do not exist in the peer chains, then create
  them.

  This also makes master_chain_id the NCS master chain. */
void copy_residue_range_from_ncs_master_to_others(int imol, const char *master_chain_id,
						  int residue_range_start, int residue_range_end);
#ifdef __cplusplus
#ifdef USE_GUILE

/*! \brief Copy chain from master to specified related NCS chains */
void copy_from_ncs_master_to_specific_other_chains_scm(int imol, const char *chain_id, SCM other_chain_id_list_scm);

/*! \brief return a list of NCS masters or scheme false */
SCM ncs_master_chains_scm(int imol);
/*! \brief Copy residue range to selected NCS chains

   If the target residues do not exist in the peer chains, then create
   them.
*/
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

/*! \brief return a list of the NCS master chain ids, or False if there are none */
PyObject *ncs_master_chains_py(int imol);
/*! \brief Copy residue range to selected NCS chains

  @param imol the model molecule index
  @param master_chain_id the NCS master chain id
  @param residue_range_start the first residue number of the range
  @param residue_range_end the last residue number of the range
  @param chain_id_list a list of the chain ids to copy to */
void copy_residue_range_from_ncs_master_to_chains_py(int imol, const char *master_chain_id,
						     int residue_range_start, int residue_range_end,
						     PyObject *chain_id_list);
/*! \brief Copy chain from master to a list of NCS chains

  @param imol the model molecule index
  @param master_chain_id the NCS master chain id
  @param chain_id_list a list of the chain ids to copy to */
void copy_from_ncs_master_to_chains_py(int imol, const char *master_chain_id,
				       PyObject *chain_id_list);
#endif
#endif

/*! \brief change the NCS master chain  (by number)

  @param ichain the index of the chain in the molecule (starting at 0)
  @param imol the model molecule index
*/
void ncs_control_change_ncs_master_to_chain(int imol, int ichain);
/*! \brief change the NCS master chain  (by chain_id)*/
void ncs_control_change_ncs_master_to_chain_id(int imol, const char *chain_id);
/*! \brief set the display state of the NCS ghost for a chain

  Note that this also sets the display state of all the NCS ghosts of
  the molecule (as \c set_draw_ncs_ghosts() does).

  @param imol the model molecule index
  @param ichain the index of the chain in the molecule (starting at 0)
  @param state 1 for on, 0 for off */
void ncs_control_display_chain(int imol, int ichain, int state);

/*! \brief set the method used to determine NCS matrices

  @param flag 0 for SSM (falls back to LSQ if Coot was compiled without
  SSM), 1 for LSQ (the default), 2 for LSQ2 */
void set_ncs_matrix_type(int flag);
/*! \brief return the method used to determine NCS matrices (0: SSM, 1: LSQ, 2: LSQ2) */
int get_ncs_matrix_state();

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief Return the NCS differences as a list.

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
/*! \brief Return the NCS differences as a list.

   e.g. ["B", "A", [[[1, ""], [1, ""], 0.4], [[2, ""], [2, ""], 0.3]]]
   i.e. ncs_related_chain its_master_chain_id and a list of residue
   info: [residue number matches: [this_resno, this_inscode,
   matching_master_resno, matching_master_inscode,
   rms_atom_position_differences]]

   When there is more than one NCS-related chain, the chain id, master
   chain id and residue info list for each are appended to the same
   (flat) list.

   @return the list, or False if imol is not a valid model molecule or
   there are no differences */
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

  Each ghost is described as [name, chain_id, target_chain_id, rtop,
  display_flag].

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

   Add a helix somewhere close to the screen centre in the refinement
   map, try to fit the orientation. A new molecule called "Helix-N" is
   created. Another new molecule called "Reverse-Helix-N" is created if
   the helix orientation isn't completely unequivocal.

   @return the index of the new molecule (the reverse helix, if one
   was made), -1 if no refinement map has been set or the helix
   could not be placed.*/
int place_helix_here();

/*! \brief add a strands

   Add a strand close to the screen centre in the refinement map, try
   to fit the orientation. A new molecule called "Strand-N" is created
   and the strand is then refined.

   n_residues is the estimated number of residues in the strand.

   n_sample_strands is the number of strands from the database tested
   to fit into this strand density.  8 is a suggested number.  20 for
   a more rigourous search, but it will be slower.

   @return the index of the new molecule, -1 on failure.*/
int place_strand_here(int n_residues, int n_sample_strands);


/*! \brief set the fudge factor for helix (and strand) placement map-level
  limits (multiplicative, default 1.0) */
void set_place_helix_here_fudge_factor(float ff);


/*! \brief show the strand placement gui.

  Choose the python version in there, if needed.  Call scripting
  function, display it in place, don't return a widget. */
void   place_strand_here_dialog();


/*! \brief autobuild helices

   Find helices in the refinement map (the whole map is searched).
   A new molecule called "SecStruc" is created (drawn as a CA trace).

   @return the index of the new molecule, -1 on failure.*/
int find_helices();

/*! \brief autobuild strands

   Find strands in the refinement map (the whole map is searched).
   A new molecule called "SecStruc" is created (drawn as a CA trace).

   @return the index of the new molecule, -1 on failure.*/
int find_strands();

/*! \brief autobuild secondary structure

   Find secondary structure in the refinement map (the whole map is
   searched). A new molecule called "SecStruc" is created (drawn as a
   CA trace).

   @param use_helix 1 to search for helices, 0 not to
   @param helix_length the helix length (in residues)
   @param helix_target one of FIND_SECSTRUC_NORMAL, FIND_SECSTRUC_STRICT,
          FIND_SECSTRUC_HI_RES or FIND_SECSTRUC_LO_RES
   @param use_strand 1 to search for strands, 0 not to
   @param strand_length the strand length (in residues)
   @param strand_target one of the FIND_SECSTRUC_ values (as for helix_target)
   @return the index of the new molecule, -1 on failure.*/
int find_secondary_structure(
    short int use_helix,  int helix_length,  int helix_target,
    short int use_strand, int strand_length, int strand_target );

/*! \brief autobuild secondary structure

   Find secondary structure in the refinement map. Parameters as for
   \c find_secondary_structure(). Note that the radius is currently
   ignored: the whole map is searched.
   A new molecule called "SecStruc" is created (drawn as a CA trace).

   @return the index of the new molecule, -1 on failure.*/
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

   Find nucleic acid chains within radius (Å) of the screen centre in
   the refinement map. Add to a molecule called "NuclAcid", create it
   if needed.

   @return the index of the "NuclAcid" molecule (even if nothing was
   found), -1 if there is no refinement map or the nautilus library
   file was not found.*/

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

  Several selections can be combined (as a logical OR) by separating
  them with "||". The view is moved to the centre of the new molecule.

@return the new molecule number, -1 means an error. */
int new_molecule_by_atom_selection(int imol, const char* atom_selection);

/*! \brief create a new molecule that consists of only the atoms
  within the given radius (r) of the given position.

  @param imol the model molecule index
  @param x the x coordinate of the centre (Å)
  @param y the y coordinate of the centre (Å)
  @param z the z coordinate of the centre (Å)
  @param r the radius in Å
  @param allow_partial_residues 1 to select only the atoms within the
         sphere, 0 to select whole residues that have atoms within the
         sphere

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

/*! \brief set the orientation matrix for the screen orientation axes

  Note: this matrix is only used by the legacy OpenGL axes drawing and
  currently has no effect. */
void
set_axis_orientation_matrix(float m11, float m12, float m13,
			    float m21, float m22, float m23,
			    float m31, float m32, float m33);

/*! \brief set the flag to use the axis orientation matrix (1 for on, 0 for off, default off)

  Note: this currently has no effect (see \c set_axis_orientation_matrix()). */
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

single_stranged_flag is 1 for a single strand, 0 for a double strand

The new molecule is centred at the screen centre.

@return the new molecule number or -1 if a problem */
int ideal_nucleic_acid(const char *RNA_or_DNA, const char *form,
		       short int single_stranged_flag,
		       const char *sequence);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE

/*! \brief get the pucker info for the specified residue

 @return scheme false if residue not found (or is a protein residue),
 otherwise, if do_pukka_pucker_check is 1,
 (list phosphate-distance puckered-atom out-of-plane-distance plane-distortion)

 (where plane-distortion is for the other 4 atoms in the plane (I think)).

 and if there is no following residue, then the phosphate distance
 cannot be calculated, so the list is null (not filled).

 If do_pukka_pucker_check is 0, the list is
 (puckered-atom out-of-plane-distance plane-distortion), with the
 phosphate distance prepended if there is a following residue.
*/
SCM pucker_info_scm(int imol, SCM residue_spec, int do_pukka_pucker_check);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief get the pucker info for the specified residue

 @return False if residue not found (or is a protein residue), otherwise,
 if do_pukka_pucker_check is 1,
 [phosphate_distance, puckered_atom, out_of_plane_distance, plane_distortion]

 (where plane_distortion is for the other 4 atoms in the plane (I think)).

 and if there is no following residue, then the phosphate distance
 cannot be calculated, so the list is empty (not filled).

 If do_pukka_pucker_check is 0, the list is
 [puckered_atom, out_of_plane_distance, plane_distortion], with the
 phosphate distance prepended if there is a following residue.
*/
PyObject *pucker_info_py(int imol, PyObject *residue_spec, int do_pukka_pucker_check);
#endif /* USE_PYTHON */
#endif /*  __cplusplus */

/*! \brief Return a molecule that contains a residue that is the WC pair
   partner of the clicked/picked/selected residue

   A new molecule called "WC partner" is created.

   @return -1 (currently always - the new molecule index is not returned) */
int watson_crick_pair(int imol, const char * chain_id, int resno);
/*! \brief add base pairs for the given residue range, modify molecule imol by creating a new chain

  @return the status from the molecule's base-pair addition (0 if imol
  is not a valid model molecule) */
int watson_crick_pair_for_residue_range(int imol, const char * chain_id, int resno_start, int resno_end);

/* not for user level */
/*! \brief start (state 1) or stop (state 0) the atom-pick mode for adding a
  base pair */
void setup_base_pairing(int state);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  sequence file (assignment)                              */
/*  ----------------------------------------------------------------------- */
/* section Sequence File (Assignment/Association) */
/*! \name Sequence File (Assignment/Association) */
/*! \{ */

/*! \brief Print the sequence to the console of the given molecule

  The sequence of the given chain is printed in FASTA format. */
void print_sequence_chain(int imol, const char *chain_id);

/*! \brief optionally write the sequence to the file for the given molecule,
    optionally in PIR format

    @param imol the model molecule index
    @param chain_id the chain id
    @param pir_format 1 for PIR format, 0 for FASTA format
    @param file_output 1 to write to file_name, 0 to print to the console
    @param file_name the output file name (used if file_output is 1) */
void print_sequence_chain_general(int imol, const char *chain_id,
                                  short int pir_format,
                                  short int file_output,
                                  const char *file_name);

/*! \brief Assign a FASTA sequence to a given chain in the  molecule

  seq is in FASTA format: a "> name" line followed by the sequence. */
void assign_fasta_sequence(int imol, const char *chain_id_in, const char *seq);
/*! \brief Assign a PIR sequence to a given chain in the molecule.  If
  the chain of the molecule already had a chain assigned to it, then
  this will overwrite that old assignment with the new one. */
void assign_pir_sequence(int imol, const char *chain_id_in, const char *seq);
/* I don't know what this does. */
/*! \brief (incomplete) assign the sequence of a chain using the map

  Currently this only passes the sequence(s) already assigned to
  chain_id to a side-chain scorer and does not change the model. */
void assign_sequence(int imol_model, int imol_map, const char *chain_id);
/*! \brief Assign a sequence to a given molecule from (whatever) sequence
  file by alignment.

  Each chain is aligned against each of the sequences in the file; a
  chain is assigned the best-aligned sequence if the alignment is good
  enough. Previously assigned sequences are cleared. */
void assign_sequence_from_file(int imol, const char *file);
/*! \brief Assign a sequence to a given molecule from a simple string

  The sequence is also assigned to the NCS-related chains of chain_id_in. */
void assign_sequence_from_string(int imol, const char *chain_id_in, const char *seq);
/*! \brief Delete all the sequences from a given molecule */
void delete_all_sequences_from_molecule(int imol);
/*! \brief Delete the sequence for a given chain_id from a given molecule */
void delete_sequence_by_chain_id(int imol, const char *chain_id_in);

/*! \brief Associate the sequence to the molecule - to be used later for sequence assignment (.c.f assign_pir_sequence)

  The file is read as PIR if its extension is ".pir", otherwise as
  FASTA. The sequence is stored with a blank chain id. */
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
  number imol. return as a list of pairs [[chain_id, seq]].  To
  be used in constructing the cootaneer gui.  Return False when no
  sequence has been assigned. */
PyObject *sequence_info_py(int imol);
/*! \brief

  do a internal alignment of all the assigned sequences,
  return a list of mismatches that need to be made to model number
  imol to match the input sequence.

Return  False on failure to align (e.g. not assigned sequence).
Note: the construction of the list of mutations, deletions and
insertions is currently disabled, so on a successful alignment an
empty list is returned.*/
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

Note: this function currently does nothing; use
\c make_molecular_surface() or \c make_electrostatic_surface() instead. */
void do_surface(int imol, int istate);
/*! \brief obsolete predicate: returns 1 (always) */
int molecule_is_drawn_as_surface_int(int imol); /* predicate */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief draw the surface of the imolth molecule clipped to the
  residues given by residue_specs.

  residue_specs must not contain spec for waters (you wouldn't want to
  surface over waters anyway).

  Note: this currently has no effect.
 */
void do_clipped_surface_scm(int imol, SCM residue_specs);
#endif /*  USE_GUILE */
#ifdef USE_PYTHON
/*! \brief draw the surface of the imolth molecule clipped to the
  residues given by residue_specs

  Note: this currently has no effect. */
void do_clipped_surface_py(int imol, PyObject *residue_specs);
#endif /*  USE_PYTHON */
#endif	/* __cplusplus */

/*! \brief make molecular surface for given atom selection

    The surface (coloured by chain) is added as a generic display object.

    per-chain functions can be added later

    @param imol the model molecule index
    @param selection_string the mmdb-format atom selection, e.g. "//A" */
void make_molecular_surface(int imol, const char *selection_string);

/*! \brief make electrostatics surface for given atom selection

  The surface (coloured by electrostatic potential) is added as a
  generic display object.

  per-chain functions can be added later

  @param imol the model molecule index
  @param selection_string the mmdb-format atom selection, e.g. "//A" */
void make_electrostatic_surface(int imol, const char *selection_string);

/*! \brief set the electrostatic surface charge range (colour scale, default 0.5)

  Note: this is not used by \c make_electrostatic_surface(). */
void set_electrostatic_surface_charge_range(float v);
/*! \brief get the electrostatic surface charge range (default 0.5) */
float get_electrostatic_surface_charge_range();

/*! \brief simple on/off screendoor transparency at the moment, an
  opacity > 0.0 (and < 1.0) will turn on screendoor transparency (stippling).

  Note: the flag is not currently used when drawing surfaces. */
void set_transparent_electrostatic_surface(int imol, float opacity);

/*! \brief return 1.0 for non transparent and 0.5 if screendoor
  transparency has been turned on (-1 if imol is not a valid model molecule). */
float get_electrostatic_surface_opacity(int imol);


/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  FFfearing                                               */
/*  ----------------------------------------------------------------------- */
/*! \name FFFearing */
/*! \{ */
/*! \brief fffear search model in molecule number imol_model in map
   number imol_map

   An FFFEAR-style rotation/translation search of the model through the map,
   using the angular step set by \c set_fffear_angular_resolution(). The
   result is installed as a new map molecule ("FFFear search results") holding
   the search score; the model itself is not moved.

   @param imol_model the index of the search model molecule
   @param imol_map the index of the map molecule to be searched
   @return the molecule index of the new results map, or -1 if either
   molecule index is not valid. */
int fffear_search(int imol_model, int imol_map);
/*! \brief set the fffear angular resolution in degrees

   @param f the angular step of the rotation search in degrees (default 15) */
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
/*! \brief try to make socket listener

   If remote control has been requested (e.g. with the \c --port command
   line option) and no listener is running yet, open a TCP socket on
   localhost (127.0.0.1 only) on the port given by
   \c get_remote_control_port_number() and poll it once a second for
   JSON-RPC requests. Otherwise do nothing. */
void make_socket_listener_maybe();
/*! \brief internal use: record the state of the listener socket

   @param sock_state the socket state flag (currently not used by the listener) */
void set_coot_listener_socket_state_internal(int sock_state);

/*! \brief feed the main thread a scheme script to evaluate

   The string is evaluated from an idle callback on the main thread. */
void set_socket_string_waiting(const char *s);
/*! \brief feed the main thread a python script to evaluate

   The string is evaluated from an idle callback on the main thread. */
void set_socket_python_string_waiting(const char *s);

/*! \brief set the port number used by the remote-control socket listener

   This does not open the socket; see \c make_socket_listener_maybe().

   @param port_number the TCP port number on localhost */
void set_remote_control_port(int port_number);
/*! \brief return the port number of the remote-control socket listener */
int get_remote_control_port_number();


/* tooltip */
/*! \brief set the "tip of the day" preference

   @param state 0 for off, any other value for on (default 1) */
void set_tip_of_the_day_flag(int state);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Display lists                                           */
/*  ----------------------------------------------------------------------- */
/* section Display Lists for Maps */
/*! \name Display Lists for Maps */
/*! \{ */
/*! \brief Should display lists be used for maps?

  Obsolete: display lists are no longer used for map drawing and this
  function now does nothing.

  @param i ignored */
void set_display_lists_for_maps(int i);

/*! \brief return the state of display_lists_for_maps.

  @return the display-lists flag (0, as it can no longer be changed) */
int display_lists_for_maps_state();

/*! \brief update the maps to the current position - rarely needed

  Recontour all map molecules around the current rotation centre. */
void update_maps();
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Browser Help                                            */
/*  ----------------------------------------------------------------------- */
/*! \name Browser Interface */
/*! \{ */
/*! \brief try to open given url in Web browser

In builds with Guile the URL is passed to the browser command set by
\c set_browser_interface() via a system call; otherwise the Python
function \c open_url() is used (and the browser command is not used). */
void browser_url(const char *url);
/*! \brief set command to open the web browser,

examples are "open" or "mozilla" (default "firefox -remote") */
void set_browser_interface(const char *browser);

/*! \brief the search interface

Split \c entry_text into words, construct a Google search URL for those
words restricted to the Coot web site, and open it with \c browser_url(). */
void handle_online_coot_search_request(const char *entry_text);
/*! \} */

// #include "c-interface-generic-objects.h"


/*  ----------------------------------------------------------------------- */
/*                  Molprobity interface                                    */
/*  ----------------------------------------------------------------------- */
/*! \name Molprobity Interface */
/*! \{ */
/*! \brief pass a filename that contains molprobity's probe output in XtalView
format

The contacts are drawn as generic display objects (lines and points),
grouped by contact type (e.g. "wide contact", "small overlap", "H-bonds"). */
void handle_read_draw_probe_dots(const char *dots_file);

/*! \brief pass a filename that contains molprobity's probe output in unformatted
format

The dots are drawn as generic display objects, one per contact type;
previously-drawn probe contact objects are cleared first.

@param dots_file the file name of the probe output (\c -unformated, i.e. colon-separated)
@param imol the model molecule to which the dots refer
@param show_clash_gui_flag currently has no effect */
void handle_read_draw_probe_dots_unformatted(const char *dots_file, int imol, int show_clash_gui_flag);


/*! \brief shall we run molprobity for on edit chi angles intermediate atoms?

@param state 1 for on, 0 for off (default 0) */
void set_do_probe_dots_on_rotamers_and_chis(short int state);
/*! \brief return the state of if run molprobity for on edit chi
  angles intermediate atoms? */
short int do_probe_dots_on_rotamers_and_chis_state();
/*! \brief shall we run molprobity after a refinement has happened?

When on, probe dots are calculated when refinement results are accepted.

@param state 1 for on, 0 for off (default 0) */
void set_do_probe_dots_post_refine(short int state);
/*! \brief show the state of shall we run molprobity after a
  refinement has happened? */
short int do_probe_dots_post_refine_state();

/*! \brief set whether Coot's own (internal) contact dots are calculated
  and drawn during interactive refinement

  @param state 1 for on and 0 for off (default 0) */
void set_do_coot_probe_dots_during_refine(short int state);

/*! \brief return whether Coot's own contact dots are calculated during
  refinement

  @return 1 for on and 0 for off */
short int get_do_coot_probe_dots_during_refine();


/*! \brief make an attempt to convert pdb hydrogen name to the name
  used in Coot (and the refmac dictionary, perhaps).

  For a 4-character name whose first character is a digit 1-4 or '*', that
  character is moved to the end (e.g. "1HB " becomes " HB1" and "1HG1"
  becomes "HG11").

  @return a newly-allocated string (allocated with \c new[]) */
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
    blur the map).

    The map's structure factors (from the original data, or calculated from
    the map the first time) are scaled and the map is regenerated in place.
    Sharpening is always applied to the original structure factors, so
    successive calls do not accumulate.

    @param imol the map molecule index
    @param b_factor the B-factor in Å^2: negative sharpens, positive blurs, 0 restores the original map */
void sharpen(int imol, float b_factor);
/*! \brief Sharpen map imol by b_factor, optionally also down-weighting
    weak reflections using a Gompertz function of F/sigF

    The Gompertz scaling is only possible when the map has F/sigF data
    associated with it (e.g. a map made with refmac parameters).

    @param imol the map molecule index
    @param b_factor the B-factor in Å^2 (negative sharpens, positive blurs)
    @param try_gompertz 1 to apply the Gompertz F/sigF scaling, 0 not to
    @param gompertz_factor currently ignored (the Gompertz parameters are fixed) */
void sharpen_with_gompertz_scaling(int imol, float b_factor, short int try_gompertz, float gompertz_factor);

/*! \brief set the limit of the b-factor map sharpening slider (default 200)

    This also sets the search range (+/- the limit) for \c optimal_B_kurtosis(). */
void set_map_sharpening_scale_limit(float f);
/*! \} */
/* ---------------------------------------------------------------------------- */
/*	Density Map Kurtosis							*/
/* ----------------------------------------------------------------------------	*/
/*! \brief find the sharpening B-factor that maximises the kurtosis of the map

    A golden-section search over B in the range +/- the map sharpening scale
    limit (see \c set_map_sharpening_scale_limit()). The map is re-sharpened
    at each trial and is left sharpened by the last trial B-factor. The
    result is cached for the molecule; later calls return the cached value
    without searching again.

    @param imol the map molecule index
    @return the optimal B-factor in Å^2 (0 if imol is not a valid map) */
float optimal_B_kurtosis(int imol);


/*  ----------------------------------------------------------------------- */
/*           Intermediate Atom Manipulation                                 */
/*  ----------------------------------------------------------------------- */

/*! \name Intermediate Atom Manipulation Interface */
/*! \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief move an intermediate (refining) atom to the given position and
  continue the refinement

  @param atom_spec an atom spec, e.g. '("A" 81 "" " CA " "")
  @param position a list of 3 numbers (x y z) in Å
  @return \#f */
SCM drag_intermediate_atom_scm(SCM atom_spec, SCM position);
#endif
#ifdef USE_PYTHON
/*! \brief move an intermediate (refining) atom to the given position and
  continue the refinement

  Does nothing (other than print a warning) if there are no intermediate
  atoms, i.e. no refinement is in progress.

  @param atom_spec an atom spec: [chain_id, res_no, ins_code, atom_name, alt_conf],
         e.g. ["A", 81, "", " CA ", ""]
  @param position a list of 3 floats [x, y, z] in Å
  @return True if the atom spec and position were well-formed, False otherwise */
PyObject *drag_intermediate_atom_py(PyObject *atom_spec, PyObject *position);

/*! \brief add a target position for an intermediate atom and refine

  A function requested by Hamish. Unlike \c add_extra_target_position_restraint(),
  this applies to the intermediate (refining) atoms, and it (re)starts the
  refinement after the target-position (atom pull) restraint is added.

  @param atom_spec an atom spec: [chain_id, res_no, ins_code, atom_name, alt_conf]
  @param position a list of 3 floats [x, y, z] in Å
  @return True if the atom spec and position were well-formed, False otherwise */
PyObject *add_target_position_restraint_for_intermediate_atom_py(PyObject *atom_spec, PyObject *position);

/*! \brief the multiple-atom version of
  \c add_target_position_restraint_for_intermediate_atom_py()

  All the restraints are added before the refinement is restarted (so that
  the refinement is not stopped and started for each one). Nothing is done if
  no refinement is in progress.

  @param atom_spec_position_list a list of [atom_spec, [x, y, z]] pairs
  @return False (always) */
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
/*! \brief mark (or unmark) an atom as fixed in refinement

  @param imol the model molecule index
  @param atom_spec an atom spec, e.g. '("A" 81 "" " CA " "")
  @param state 1 to fix the atom, 0 to unfix it
  @return \#f */
SCM mark_atom_as_fixed_scm(int imol, SCM atom_spec, int state);
/*! \brief mark (or unmark) a list of atoms as fixed in refinement

  @param imol the model molecule index
  @param atom_spec_list a list of atom specs
  @param state 1 to fix the atoms, 0 to unfix them
  @return the length of atom_spec_list */
int mark_multiple_atoms_as_fixed_scm(int imol, SCM atom_spec_list, int state);
#endif
#ifdef USE_PYTHON
/*! \brief mark (or unmark) an atom as fixed in refinement

  @param imol the model molecule index
  @param atom_spec an atom spec: [chain_id, res_no, ins_code, atom_name, alt_conf]
  @param state 1 to fix the atom, 0 to unfix it
  @return True if the atom spec was well-formed, False otherwise */
PyObject *mark_atom_as_fixed_py(int imol, PyObject *atom_spec, int state);
/*! \brief mark (or unmark) a list of atoms as fixed in refinement

  @param imol the model molecule index
  @param atom_spec_list a list of atom specs [chain_id, res_no, ins_code, atom_name, alt_conf]
  @param state 1 to fix the atoms, 0 to unfix them
  @return the number of well-formed atom specs that were processed */
int mark_multiple_atoms_as_fixed_py(int imol, PyObject *atom_spec_list, int state);
#endif
#endif /* c++ */

/*! \brief enter (or leave) the mode where the next atom picks fix or
  unfix atoms

  @param ipick 1 to start picking, 0 to stop
  @param is_unpick 1 if picked atoms are to be unfixed, 0 if they are to be fixed */
void setup_fixed_atom_pick(short int ipick, short int is_unpick);

/*! \brief clear all fixed atoms

  @param imol the model molecule index */
void clear_all_fixed_atoms(int imol);
/*! \brief clear the fixed atoms of all model molecules */
void clear_fixed_atoms_all();

/*! \brief produce debugging output from problematic atom picking

  @param istate 1 for on, 0 for off (default 0) */
void set_debug_atom_picking(int istate);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Partial Charge                                          */
/*  ----------------------------------------------------------------------- */
/*! \name Partial Charges */
/*! \{ */
/*! \brief show the partial charges for the residue of the given specs
   (charges are read from the dictionary)

   Currently not implemented: this checks that the dictionary for the
   residue type is available but displays nothing. */
void show_partial_charge_info(int imol, const char *chain_id, int resno, const char *ins_code);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  EM Interface                                            */
/*  ----------------------------------------------------------------------- */
/*! \name EM interface */
/*! \{ */
/*! \brief Scale the cell, for use with EM maps, where the cell needs
   to be adjusted.  Use like:  (scale-cell 2 1.012 1.012 1.012).

   The cell lengths a, b and c are multiplied by fac_u, fac_v and fac_w
   (the angles and the grid sampling are unchanged) and the map is recontoured.

   @param imol_map the map molecule index
   @param fac_u the scale factor for a
   @param fac_v the scale factor for b
   @param fac_w the scale factor for c
   @return 0 (currently the success status is not set, so 0 is returned
   whether or not the cell was scaled) */
int scale_cell(int imol_map, float fac_u, float fac_v, float fac_w);

/*! \brief create a number of maps by segmenting the given map

   Segment the density above the (absolute) low_level into connected regions.
   Each segment (up to 300) is installed as a new map molecule ("Map N Segment M")
   on the same grid as the input map, with the contour level of the input map.

   @param imol_map the map molecule index
   @param low_level the absolute density threshold (not in units of rmsd) */
void segment_map(int imol_map, float low_level);

/*! \brief create maps by "scale-space" segmentation of the given map

   As \c segment_map(), but segments are merged by progressively blurring
   the map over n_rounds rounds. At most 8 new segment maps are created.

   @param imol_map the map molecule index
   @param low_level the absolute density threshold (not in units of rmsd)
   @param b_factor_inc the Gaussian blurring increment applied in each round
   @param n_rounds the number of rounds of blurring */
void segment_map_multi_scale(int imol_map, float low_level, float b_factor_inc, int n_rounds);

/*! \brief make a map histogram

   Display a histogram of the density values of the map in a dialog
   (only in builds with goocanvas and when the graphics interface is in use). */
void map_histogram(int imol_map);

/*! \brief ignore pseudo-zeros when calculationg maps stats (default 1 = true)

   @param state 1 to ignore pseudo-zero grid points in the map mean and
   rmsd calculations, 0 to include them */
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

use state=0 to turn it off (default 1). Note that the current file choosers
do not consult this setting. */
void set_add_ccp4i_projects_to_file_dialogs(short int state);

/*! \brief write a ccp4mg picture description file

The file describes the current view (centre and zoom), the background colour,
the bond width and the model and map molecules. */
void write_ccp4mg_picture_description(const char *filename);

/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Dipoles                                                 */
/*  ----------------------------------------------------------------------- */
/*! \name Dipoles */
/*! \{ */
/*! \brief delete a dipole

  @param imol the model molecule index
  @param dipole_number the dipole number (as returned by the add_dipole functions) */
void delete_dipole(int imol, int dipole_number);
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief generate a dipole from all atoms in the given
  residues. Return the dipole description

  The partial charges are taken from the dictionary.

  @return a list (dipole_number (x y z)) where (x y z) is the dipole vector
  (dipole_number is -1 on failure), or \#f if imol is not a valid model molecule */
SCM add_dipole_for_residues_scm(int imol, SCM residue_specs);
/*! \brief generate a dipole from all atoms in the given residue

  @return a list (dipole_number (x y z)) where (x y z) is the dipole vector
  (dipole_number is -1 on failure), or \#f if imol is not a valid model molecule */
SCM add_dipole_scm(int imol, const char* chain_id, int res_no, const char *ins_code);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief generate a dipole from all atoms in the given residue

  The partial charges are taken from the dictionary.

  @return a list [dipole_number, [x, y, z]] where [x, y, z] is the dipole vector
  (dipole_number is -1 on failure), or False if imol is not a valid model molecule */
PyObject *add_dipole_py(int imol, const char* chain_id, int res_no,
			const char *ins_code);
/*! \brief add a dipole given a set of residues.  Return a dipole
  description.

  @param imol the model molecule index
  @param residue_specs a list of residue specs
  @return a list [dipole_number, [x, y, z]] where [x, y, z] is the dipole vector
  (dipole_number is -1 on failure), or False if imol is not a valid model molecule */
PyObject *add_dipole_for_residues_py(int imol, PyObject *residue_specs);
#endif /* USE_PYTHON */
#endif /* c++ */
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  Patterson                                               */
/*  ----------------------------------------------------------------------- */
/*! \brief Make a patterson molecule

Calculate a Patterson map from the amplitudes in an MTZ file and install it
as a new map molecule.

@param mtz_file_name the MTZ file name
@param f_col the amplitude (F) column label
@param sigf_col the sigma(F) column label
@return a new molecule number or -1 on failure */
int make_and_draw_patterson(const char *mtz_file_name,
			    const char *f_col,
			    const char *sigf_col);
/*! \brief Make a patterson molecule using intensities

Calculate a Patterson map from the intensities in an MTZ file and install it
as a new map molecule.

@param mtz_file_name the MTZ file name
@param i_col the intensity (I) column label
@param sigi_col the sigma(I) column label
@return a new molecule number or -1 on failure */
int make_and_draw_patterson_using_intensities(const char *mtz_file_name,
					      const char *i_col,
					      const char *sigi_col);

/*  ----------------------------------------------------------------------- */
/*                  Laplacian                                               */
/*  ----------------------------------------------------------------------- */
/*! \name Aux functions */
/*! \{ */
/*! \brief Create the "Laplacian" (-ve second derivative) of the given map.

@param imol the map molecule index
@return the molecule index of the new map, or -1 if imol is not a valid map */
int laplacian (int imol);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  PKGDATADIR                                              */
/*  ----------------------------------------------------------------------- */
/*! \name PKGDATADIR */
/*! \{ */
#ifdef __cplusplus
#ifdef USE_PYTHON
/*! \brief return the Coot package data directory (e.g. \c $prefix/share/coot) as a string */
PyObject *get_pkgdatadir_py();
#endif /* USE_PYTHON */
#ifdef USE_GUILE
// note: built-ins: (%package-data-dir) and %guile-build-info
/*! \brief return the Coot package data directory (e.g. \c $prefix/share/coot) as a string */
SCM get_pkgdatadir_scm();
#endif
#endif /*  __cplusplus */
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  SMILES                                                  */
/*  ----------------------------------------------------------------------- */
/*! \name SMILES */
/*! \{ */
/*! \brief display the SMILES string dialog

This runs the old PyGTK/guile-gtk scripted dialog, so in current builds
(which have neither) it does nothing. */
void do_smiles_gui();
/*! \} */
/*  ----------------------------------------------------------------------- */
/*                  Fun                                                     */
/*  ----------------------------------------------------------------------- */
/* section Fun */
/*! \brief does nothing */
void do_tw();

/*  ----------------------------------------------------------------------- */
/*                  Phenix Support                                          */
/*  ----------------------------------------------------------------------- */
/*! \name PHENIX Support */
/*! \{ */
/*! \brief set the button label of the external Refinement program

(The label is stored, but no current GUI element uses it.) */
void set_button_label_for_external_refinement(const char *button_label);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Text                                                    */
/*  ----------------------------------------------------------------------- */
/*! \name Graphics Text */
/*! \{ */
/*! \brief Put text at x,y,z

@param text the text to display
@param x the x coordinate in Å
@param y the y coordinate in Å
@param z the z coordinate in Å
@param size currently ignored
@return a text handle */
int place_text(const char*text, float x, float y, float z, int size);

/*! \brief Remove "3d" text item

@param text_handle the handle returned by \c place_text() */
void remove_text(int text_handle);

/*! \brief change the string of a "3d" text item

@param text_handle the handle returned by \c place_text()
@param new_text the replacement text */
void edit_text(int text_handle, const char *new_text);

/*! \brief return the closest text that is with r A of the given
  position.  If no text item is close, then return -1

  @return the index of the closest text item within r Å of (x, y, z), or -1 */
int text_index_near_position(float x, float y, float z, float r);
/*! \} */

/*  ----------------------------------------------------------------------- */
/*                  PISA Interface                                      */
/*  ----------------------------------------------------------------------- */
/*! \name PISA Interaction */
/*! \{ */
/*! \brief return the molecule number of the interacting
  residues. Return -1 if no new model was created. Old, not very useful.

  The residues of each molecule that are within 4 Å of the other molecule
  are copied into new molecules ("interacting residues from N"): one for
  imol_1 and (if there are any) one for imol_2.

  @return the molecule index of the new molecule of residues from imol_1 */
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
/*! \brief internal function: return the residues field of a PISA molecule record */
SCM pisa_molecule_record_residues(SCM molecule_record_1);
/*! \brief internal function: return the chain-id field of a PISA molecule record */
SCM pisa_molecule_record_chain_id(SCM molecule_record_1);
/*! \brief internal function: draw a PISA interface bond

  The bond is added to a generic display object for that interface and bond
  type ("H-bonds-interface-N", "salt-bridges-interface-N", "SS-bonds-interface-N"
  or "Covalent-interface-N"), which are created if needed.

  @param imol_1 the molecule index of the first atom
  @param imol_2 the molecule index of the second atom
  @param pisa_bond_scm a list of the bond type ('h-bonds, 'salt-bridges,
  'cov-bonds or 'ss-bonds), an atom spec in imol_1 and an atom spec in imol_2
  @param interface_number the interface number, used to name the generic objects */
void add_pisa_interface_bond_scm(int imol_1, int imol_2, SCM pisa_bond_scm,
				 int interface_number);


/*! \brief clear out and undisplay all pisa interface descriptions.

  Currently not implemented: this does nothing. */
void pisa_clear_interfaces();
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief the scripting interface, called from parsing the PISA XML
   interface description

   interfaces_description_py is a list of interface descriptions; each is a
   list of 6 items: [molecules, bonds, area, solv_en, pvalue, stab_en], where
   molecules is a list of 2 molecule records [imol, symop, molecule_dictionary].
   The interface (dots) is made from the residues of the molecule records.
   Note of course that imol_2 (or 1) can be a symmetry copy of (part of)
   mol_1 (or 2). If any interfaces were found, an interfaces dialog is shown.

   @return the dot indexes (currently -1)

*/
PyObject *handle_pisa_interfaces_py(PyObject *interfaces_description_py);

/* internal function */
/* PyObject *pisa_molecule_record_residues_py(PyObject *molecule_record_1); */
/* PyObject *pisa_molecule_record_chain_id_py(PyObject *molecule_record_1); */
/*! \brief internal function: draw a PISA interface bond

  The bond is added to a generic display object for that interface and bond
  type ("H-bonds-interface-N", "salt-bridges-interface-N", "SS-bonds-interface-N"
  or "Covalent-interface-N"), which are created if needed.

  @param imol_1 the molecule index of the first atom
  @param imol_2 the molecule index of the second atom
  @param pisa_bond_py a list of the bond type ("h-bonds", "salt-bridges",
  "cov-bonds" or "ss-bonds"), an atom spec in imol_1 and an atom spec in imol_2
  @param interface_number the interface number, used to name the generic objects */
void add_pisa_interface_bond_py(int imol_1, int imol_2, PyObject *pisa_bond_py,
                                 int interface_number);

/*! \brief clear out and undisplay all pisa interface descriptions.

  Currently not implemented: this does nothing. */
void pisa_clear_interfaces();
#endif /* USE_PYTHON */
#endif /* c++ */


/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  Jiggle fit                                              */
/*  ----------------------------------------------------------------------- */
/*! \name Jiggle Fit */
/*! \{ */

/*!  \brief jiggle fit to the current refinment map
 *
 * Rigid-body fit the residue to the refinement map (see
 * \c set_imol_refinement_map()) by scoring randomly rotated and translated
 * trial positions; the best-scoring position is applied to the residue.
 *
 * @param imol the model molecule index
 * @param chain_id the chain id
 * @param resno the residue number
 * @param ins_code the insertion code
 * @param n_trials the number of random trials
 * @param jiggle_scale_factor scales the size of the random rotations and translations
 * @return -999 if not possible (e.g. no refinement map or residue not found),
 * else return the new best fit (density score) for this residue.  */
float fit_to_map_by_random_jiggle(int imol, const char *chain_id, int resno, const char *ins_code,
                                  int n_trials, float jiggle_scale_factor);

/*!  \brief jiggle fit the molecule to the current refinment map.

  The trials are scored using the main-chain (and nucleotide base/backbone)
  atoms; the best transformation is applied to all the chains of the molecule.

  @param imol the model molecule index
  @param n_trials the number of random trials
  @param jiggle_scale_factor scales the size of the random rotations and translations
  @return -999 if not possible (e.g. the refinement map is not set), else
  return the new best fit for this molecule.  */
float fit_molecule_to_map_by_random_jiggle(int imol, int n_trials, float jiggle_scale_factor);
/*!  \brief jiggle fit the molecule to the current refinment map
 *
 * As \c fit_molecule_to_map_by_random_jiggle(), but the fit is first done
 * against a copy of the refinement map that is blurred by map_blur_factor,
 * followed by a short (12-trial) fit against the unblurred map.
 *
 * @param imol the model molecule index
 * @param n_trials the number of random trials
 * @param jiggle_scale_factor scales the size of the random rotations and translations
 * @param map_blur_factor the B-factor (Å^2) used to blur the map
 * @return -100 if not possible, else return the new best fit for this molecule
 * (from the final fit against the unblurred map) */
float fit_molecule_to_map_by_random_jiggle_and_blur(int imol, int n_trials, float jiggle_scale_factor, float map_blur_factor);

/*!  \brief jiggle fit the chain to the current refinment map.

  @param imol the model molecule index
  @param chain_id the chain id
  @param n_trials the number of random trials
  @param jiggle_scale_factor scales the size of the random rotations and translations
  @return currently always -999 (the fit score is not passed back) */
float fit_chain_to_map_by_random_jiggle(int imol, const char *chain_id, int n_trials, float jiggle_scale_factor);

/*!  \brief jiggle fit the chain to the current refinment map
 *
 * Use a map that is blurred by the give factor for fitting.
 *
 * @param imol the model molecule index
 * @param chain_id the chain id
 * @param n_trials the number of random trials
 * @param jiggle_scale_factor scales the size of the random rotations and translations
 * @param map_blur_factor the B-factor (Å^2) used to blur the map
 * @return currently always -100 (the fit score is not passed back) */
float fit_chain_to_map_by_random_jiggle_and_blur(int imol, const char *chain_id, int n_trials, float jiggle_scale_factor, float map_blur_factor);

/*! \brief Patterson overlap plus phased translationn function MR-like local fitting
 *
 * Use the imol_refinement map. The search is centred on the current
 * rotation (screen) centre; the coordinates of molecule imol are replaced
 * by those of the best solution.
 *
 * @param imol the molecule index
 * @param n_top_rotations use only the top n_top_rotations rotation solutions
 * @param n_top_translations use only the top n_top_translation translation solutions
 * */
void molecular_replacement_fit_about_screen_centre(int imol, int n_top_rotations, int n_top_translations);

/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  SBase interface                                         */
/*  ----------------------------------------------------------------------- */
/*! \name SBase interface */
/*! \{ */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return a list of compoundIDs in SBase (CCP4 SRS) of which the
  given string is a substring of the compound name */
SCM matching_compound_names_from_sbase_scm(const char *compound_name_fragment);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
/*! \brief return a list of compoundIDs in SBase (CCP4 SRS) of which the
  given string is a substring of the compound name */
PyObject *matching_compound_names_from_sbase_py(const char *compound_name_fragment);
#endif /* USE_PYTHON */
#endif /*__cplusplus */

/*! \brief return the new molecule number of the monomer.

Get the monomer from the CCP4 SRS (SBase), install it as a new molecule
at the screen centre, and fill the dictionary for it from the SRS too.

The monomer will have chainid "A" and residue number 1.

Return -1 on failure to get monomer. */
int get_ccp4srs_monomer_and_dictionary(const char *comp_id);

/*! \brief same as above but using old name for back-compatibility */
int get_sbase_monomer(const char *comp_id);

/*! \} */


/* Needs a/the correct section */
/*! \brief add a linked residue based purely on dictionary template.

   For addition of NAG to ASNs typically.

   This doesn't work with residues with alt confs.

   The dictionary for new_residue_comp_id is read if needed. If the
   refinement map is set, the new residue (and the residue it is linked to)
   are then torsion-fitted to the map.

   @param imol the model molecule index
   @param chain_id the chain id of the residue to which the new residue is linked
   @param resno the residue number of that residue
   @param ins_code the insertion code of that residue
   @param new_residue_comp_id the residue type of the new residue, e.g. "NAG"
   @param link_type the dictionary link type (chem_link id), e.g. "pyr-ASN"
   @param n_trials the number of trials in the torsion fit
   @return currently always 0 (the success status is not set)
*/
int add_linked_residue(int imol, const char *chain_id, int resno, const char *ins_code,
		       const char *new_residue_comp_id, const char *link_type, int n_trials);
#ifdef __cplusplus
#ifdef USE_GUILE
// mode is either 1: add  2: add and fit  3: add, fit and refine
/*! \brief add a linked residue based on the dictionary template, and
  optionally fit and refine it

  @param imol the model molecule index
  @param chain_id the chain id of the residue to which the new residue is linked
  @param resno the residue number of that residue
  @param ins_code the insertion code of that residue
  @param new_residue_comp_id the residue type of the new residue, e.g. "NAG"
  @param link_type the dictionary link type (chem_link id), e.g. "pyr-ASN"
  @param mode 1: add, 2: add and fit (to the refinement map), 3: add, fit and refine
  @return the residue spec of the new residue (in modes 2 and 3, if it was
  added), otherwise \#f */
SCM add_linked_residue_scm(int imol, const char *chain_id, int resno, const char *ins_code,
			   const char *new_residue_comp_id, const char *link_type, int mode);
#endif
#ifdef USE_PYTHON
/*! \brief add a linked residue based on the dictionary template, and
  optionally fit and refine it

  Whether the new residue is fitted and refined (two rounds, if the
  refinement map is set) is controlled by
  \c set_add_linked_residue_do_fit_and_refine(), not by mode.

  @param imol the model molecule index
  @param chain_id the chain id of the residue to which the new residue is linked
  @param resno the residue number of that residue
  @param ins_code the insertion code of that residue
  @param new_residue_comp_id the residue type of the new residue, e.g. "NAG"
  @param link_type the dictionary link type (chem_link id), e.g. "pyr-ASN"
  @param mode currently ignored
  @return the residue spec of the new residue if fit-and-refine is on and the
  residue was added, otherwise False */
PyObject *add_linked_residue_py(int imol, const char *chain_id, int resno, const char *ins_code,
				const char *new_residue_comp_id, const char *link_type, int mode);
#endif
#endif
/*! \brief set whether \c add_linked_residue_py() fits and refines the new residue

  @param state 1 for on, 0 for off (default 1) */
void set_add_linked_residue_do_fit_and_refine(int state);

/*  ----------------------------------------------------------------------- */
/*               Flattened Ligand Environment View  Interface               */
/*  ----------------------------------------------------------------------- */
/*! \name FLE-View */
/*! \{ */

/*! \brief show the 2D Flattened Ligand Environment View of the given ligand

The depiction (made with RDKit) of the ligand and its interactions with the
surrounding residues is displayed as an SVG in a dialog.

@param imol the model molecule index
@param chain_id the chain id of the ligand
@param res_no the residue number of the ligand
@param ins_code the insertion code of the ligand
@param dist_max the radius (in Å) within which residues are considered to be in the environment */
void fle_view(int imol, const char *chain_id, int res_no, const char *ins_code, float dist_max);

/*! \brief obsolete: this function has been removed and now does nothing (use \c fle_view()) */
void fle_view_with_rdkit(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius);
/*! \brief obsolete: this function has been removed and now does nothing (no file is written) */
void fle_view_with_rdkit_to_png(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *png_file_name);
/*! \brief obsolete: this function has been removed and now does nothing (no file is written) */
void fle_view_with_rdkit_to_svg(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *svg_file_name);

/*! \brief obsolete: this function has been removed and now does nothing */
void fle_view_with_rdkit_internal(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *file_format, const char *file_name);

/*! \brief set the maximum considered distance to water

default 3.25 A. (Note that the current \c fle_view() uses its own
internal value and does not consult this setting.) */
void fle_view_set_water_dist_max(float dist_max);
/*! \brief set the maximum considered hydrogen bond distance

default 3.9 A. (Note that the current \c fle_view() uses its own
internal value and does not consult this setting.) */
void fle_view_set_h_bond_dist_max(float h_bond_dist_max);

/*! \brief Add hydrogens to specificied residue

The hydrogens are added using the dictionary for the residue type (a
backup is made first). This needs RDKit (an "enhanced ligand tools" build);
otherwise it fails. On failure, the reason is shown in an info dialog.

@return success status: 1 on success, 0 on failure.
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

  The moving molecule is transformed in place (a backup is made first).

  @param imol_ref the reference model molecule index
  @param ref_selection an mmdb atom selection string for the reference molecule
  @param imol_moving the moving model molecule index
  @param moving_selection an mmdb atom selection string for the moving molecule
  @param n_res currently ignored (the implementation uses fragments of 6 residues)
  @param dist_crit currently ignored

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
/*! \brief stop the scripted multi-residue refinement

  Sets the scripting-layer continue flag to false and makes the "continue"
  and "cancel" buttons available. */
void toolbar_multi_refine_stop();
/*! \brief continue the scripted multi-residue refinement

  Sets the scripting-layer continue flag to true and re-adds the
  multi-refine idle function. */
void toolbar_multi_refine_continue();
/*! \brief cancel the scripted multi-residue refinement

  Sets the scripting-layer continue flag to false and hides the
  multi-refine toolbar buttons. */
void toolbar_multi_refine_cancel();
/*! \brief show or hide the multi-refine "stop" toolbar button

  @param state 1 to show, 0 to hide */
void set_visible_toolbar_multi_refine_stop_button(short int state);
/*! \brief show or hide the multi-refine "continue" toolbar button

  This also makes the "cancel" button insensitive.

  @param state 1 to show, 0 to hide */
void set_visible_toolbar_multi_refine_continue_button(short int state);
/*! \brief show or hide the multi-refine "cancel" toolbar button

  @param state 1 to show, 0 to hide */
void set_visible_toolbar_multi_refine_cancel_button(short int state);
/*! \brief set the sensitivity of a multi-refine toolbar button

  @param button_type one of "stop", "continue", "cancel"
  @param state 1 for sensitive, 0 for insensitive */
void toolbar_multi_refine_button_set_sensitive(const char *button_type, short int state);

/*! \brief load tutorial model and data
 *
 * Loads an example dataset - the sample is an RNase structure (model and maps) and is
 * used for learning and testing.
 *
 * The model (\c tutorial-modern.pdb) and the 2mFo-DFc (FWT/PHWT) and
 * difference (DELFWT/PHDELWT) maps (from \c rnasa-1.8-all_refmac1.mtz) are
 * read from the \c data directory of the Coot package data directory.
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
/*! \brief put molecule number imol to display only model number imodel

@param imol the model molecule index
@param imodel the (1-based) model number, or 0 to display all models */
void single_model_view_model_number(int imol, int imodel);
/*! \brief the current model number being displayed

@return the model number, or 0 if all models are displayed or imol is not
a valid model molecule. */
int single_model_view_this_model_number(int imol);
/*! \brief change the representation to the next model number to be displayed

The cycle includes the "all models" state: ... model n, all models, model 1, ...

@return the new model number, 0 for "all models" (also returned for a
non-multimodel-molecule).
*/
int single_model_view_next_model_number(int imol);
/*! \brief change the representation to the previous model number to be displayed

The cycle includes the "all models" state: ... model 1, all models, model n, ...

@return the new model number, 0 for "all models" (also returned for a
non-multimodel-molecule). */
int single_model_view_prev_model_number(int imol);
/*! \} */



/*  ----------------------------------------------------------------------- */
/*                  update self                                             */
/*  ----------------------------------------------------------------------- */
/* this function is here because it is called by c_inner_main() (ie. need a c interface). */
/*! \brief if \c --update-self was given on the command line, run the
  scripting-layer \c update_self() function (only in builds with libcurl) */
void run_update_self_maybe(); /* called when --update-self given at command line */

/*  ----------------------------------------------------------------------- */
/*                    keyboarding mode                                      */
/*  ----------------------------------------------------------------------- */
/*! \brief show the keyboard "go to residue" entry window */
void show_go_to_residue_keyboarding_mode_window();
/*! \brief go to the residue described by text

  text is interpreted relative to the molecule (and chain) of the active
  atom: either a residue/atom specifier or a sequence triplet (3 residue
  one-letter codes). The rotation centre is moved to the matching atom.
  Nothing happens if there is no active atom. */
void handle_go_to_residue_keyboarding_mode(const char *text); /* should this be here? */

/*  ----------------------------------------------------------------------- */
/*                    graphics ligand view                                  */
/*  ----------------------------------------------------------------------- */
/*! \name graphics 2D ligand view */
/*! \{ */
/*! \brief set the graphics ligand view state

 @param state 1 for on, 0 for off (default is 1 (on)). */
void set_show_graphics_ligand_view(int state);
/*! \} */


/*  ----------------------------------------------------------------------- */
/*                  experimental                                            */
/*  ----------------------------------------------------------------------- */
/*! \name Experimental */
/*! \{ */

/*! \brief fetch and superpose AlphaFold models for the molecule of the active atom

  For each chain with a UniProt (UNP) DBREF record in the header, the
  AlphaFold model for that UniProt accession is fetched, superposed onto
  the chain and shown in CA + ligands representation. An info dialog is
  shown if no DBREF is found. Nothing happens if there is no active atom. */
void fetch_and_superpose_alphafold_models_using_active_molecule();

// void add_ligand_builder_menu_item_maybe(); // what does this do?

/*!  \brief display the ligand builder dialog

  This launches Layla (needs GTK 4.10 or later). */
void start_ligand_builder_gui();

/*  ----------------------------------------------------------------------- */
/*                  end                                                     */
/*  ----------------------------------------------------------------------- */
#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief return the rotamer score for the whole molecule

  @param imol the model molecule index
  @return a list (score n_rotamer_residues), where score is the sum of the
  log rotamer probabilities, or \#f if imol is not a valid model molecule */
SCM all_molecule_rotamer_score(int imol);
/*! \brief return the Ramachandran score for the whole molecule

  The Scheme version of \c all_molecule_ramachandran_score_py(): a 6-element
  list, but per-residue entries without residue names have 3 items
  (phi-psi residue-spec score), and \#f if imol is not a valid model molecule. */
SCM all_molecule_ramachandran_score(int imol);
#endif /* USE_GUILE */

#ifdef USE_PYTHON
/*! \brief return the rotamer score for the whole molecule

  The score uses the rotamer probability tables; residues that are always
  "pass" (e.g. GLY, PRO, ALA) are not counted.

  @param imol the model molecule index
  @return a list [score, n_rotamer_residues], where score is the sum of the
  log rotamer probabilities of the scored residues, or False if imol is not a
  valid model molecule. (If the rotamer tables are not available the list is
  [0.0, 0].) */
PyObject *all_molecule_rotamer_score_py(int imol);
#endif /* USE_PYTHON */

#ifdef USE_PYTHON
/*! \brief return the Ramachandran scores for the whole molecule

  @param imol the model molecule index
  @return False if imol is not a valid model molecule, otherwise a list of 6 items:
    - 0: the overall Ramachandran score (float)
    - 1: the number of residues scored (int)
    - 2: the score of the non-secondary-structure residues (float)
    - 3: the number of non-secondary-structure residues scored (int)
    - 4: the number of residues with (near) zero probability (int)
    - 5: a list with one entry per scored residue:
      [[phi, psi], residue_spec, probability, [prev_res_name, this_res_name, next_res_name]],
      or -1 if the neighbouring residues are not available.

  Note that currently only items 1 and 5 are calculated: items 0, 2, 3 and 4 are
  always 0. The per-residue probability is the (clipper) Ramachandran probability
  for the residue type (GLY, PRO, ILE/VAL, pre-PRO or general). */
PyObject *all_molecule_ramachandran_score_py(int imol);
#endif /* USE_PYTHON */

#ifdef USE_PYTHON
/*! \brief return the Ramachandran region of each residue of the molecule

  @param imol the model molecule index
  @return a list of (residue_spec, region) tuples, where region is a
  \c coot::rama_plot region code (preferred, allowed or outlier), or False if
  imol is not a valid model molecule or if the list is empty. Note that the
  current Ramachandran scoring does not fill the region list, so at present
  this always returns False. */
PyObject *all_molecule_ramachandran_region_py(int imol);
#endif /* USE_PYTHON */
#endif /* __cplusplus */

/*! \brief globularize the molecule.

This is not guaranteed to generate the correct biological entity, but will bring together
molecules (chains/domains) that are dispersed throughout the unit cell.
The molecule is modified in place (a backup is made first).

@param imol the molecule index.
*/
void globularize(int imol);

#ifdef __cplusplus
#ifdef USE_GUILE
/*! \brief run a user defined function

      Define a function func which runs after the user has made
      n_clicks atom picks.  func is called with one argument: a list
      of the picked atom specifiers, each with a leading model number,
      i.e. (model_number imol chain_id res_no ins_code atom_name alt_conf).

      @param n_clicks the number of atom picks (must be greater than 0)
      @param func the function to run
*/
void user_defined_click_scm(int n_clicks, SCM func);
#endif
#ifdef USE_PYTHON
/*! \brief run a user defined function

      Define a function func which runs after the user has made
      n_clicks atom picks. The picked atoms are labelled. func is called
      with n_clicks arguments, one per picked atom, each an atom specifier
      with a leading model number:
      [model_number, imol, chain_id, res_no, ins_code, atom_name, alt_conf].

      @param n_clicks the number of atom picks (must be greater than 0)
      @param func a callable
*/
void user_defined_click_py(int n_clicks, PyObject *func);
#endif /* PYTHON */
#endif /* __cplusplus */

/*! \} */

#endif /* C_INTERFACE_H */
