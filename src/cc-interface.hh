/* src/cc-interface.hh
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2008, 2009, 2010, 2011, 2012 by The University of Oxford
 * Copyright 2015 by Medical Research Council
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifndef CC_INTERFACE_HH
#define CC_INTERFACE_HH

#ifdef USE_PYTHON
#include "Python.h"
#endif

#include <gtk/gtk.h>
#include <optional>

#ifdef USE_GUILE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#include <libguile.h>
#pragma GCC diagnostic pop
#endif // USE_GUILE

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-density-stats.hh"

#include "ligand/dipole.hh"
#include "high-res/sequence-assignment.hh" // for residue_range_t

#include "coords/mmdb-extras.hh"
#include "coords/mmdb-crystal.hh"

#include "pli/flev-annotations.hh" // animated ligand interactions
#include "named-rotamer-score.hh"

#include "coords/phenix-geo.hh"
#include "gtk-utils.hh"

/*! \file
  \brief Coot Scripting Interface - General (C++ functions)
*/

namespace coot {

   //! alias path
   class alias_path_t {
   public:
      int index;
      std::string s;
      bool flag;
      alias_path_t(int index_in, const std::string &s_in, bool flag_in) : s(s_in) {
         index = index_in;
         flag = flag_in;
      }
   };

   //! pisa internal function
   //
   class pisa_interface_bond_info_t {
   public:
      pisa_interface_bond_info_t() {
         n_h_bonds = 0;
         n_salt_bridges = 0;
         n_cov_bonds = 0;
         n_ss_bonds = 0;
      }
      int n_h_bonds;
      int n_salt_bridges;
      int n_cov_bonds;
      int n_ss_bonds;
   };

   #ifdef USE_GUILE
   pisa_interface_bond_info_t get_pisa_interface_bond_info_scm(SCM bonds_info_scm);
   #endif
   #ifdef USE_PYTHON
   pisa_interface_bond_info_t get_pisa_interface_bond_info_py(PyObject *bonds_info_py);
   #endif

}

//! \brief Get the Git commit hash of the Coot build
//!
//! Returns the Git commit identifier for the version of Coot currently running.
//! Useful for bug reports and version tracking.
//!
//! \return String containing the Git commit hash (e.g., "a3f2c1d")
//!
//! Example:
//! \code{.cpp}
//! std::string commit = git_commit();
//! std::cout << "Coot version: commit " << commit << std::endl;
//! \endcode
std::string git_commit();

//! \brief Filter files by glob pattern
//!
//! Returns a list of files in the specified directory that match the given
//! glob pattern, useful for file selection dialogs.
//!
//! @param pre_directory Directory path to search
//! @param data_type Type of data to filter (specific values TBD)
//!
//! @return Vector of filename strings matching the criteria
std::vector<std::string> filtered_by_glob(const std::string &pre_directory, int data_type);


//! \brief Check if a string exists in a vector
//!
//! Searches for an exact match of the search string in the provided list.
//!
//! @param search String to search for
//! @param list Vector of strings to search within
//!
//! @return 1 if found, 0 if not found
//!
//! Example:
//! \code{.cpp}
//! std::vector<std::string> chains = {"A", "B", "C"};
//! if (string_member("B", chains)) {
//!     std::cout << "Chain B exists" << std::endl;
//! }
//! \endcode
short int string_member(const std::string &search, const std::vector<std::string> &list);

//! \brief Compare two strings
//!
//! Performs string comparison for sorting purposes.
//!
//! @param a First string
//! @param b Second string
//!
//! @return true if a < b in lexicographic order
//!
//! \note Useful as a comparator function for std::sort
bool compare_strings(const std::string &a, const std::string &b);

/*
std::string pre_directory_file_selection(GtkWidget *sort_button);
void filelist_into_fileselection_clist(GtkWidget *fileselection, const std::vector<std::string> &v);
*/

// These widget declarations don't belong in this file.
// void
// add_validation_mol_menu_item(int imol, const std::string &name, GtkWidget *menu, GtkSignalFunc callback);
// void create_initial_validation_graph_submenu_generic(GtkWidget *window1,
// 						     const std::string &menu_name,
// 						     const std::string &sub_menu_name);


// CaBLAM
std::vector<std::pair<coot::residue_spec_t, double> >
add_cablam_markup(int imol, const std::string &cablam_file_name);
#ifdef USE_GUILE
SCM add_cablam_markup_scm(int imol, const std::string &cablam_log_file_name);
#endif
#ifdef USE_PYTHON
//! \brief Add CaBLAM validation markup (Python interface)
//!
//! Reads CaBLAM output and adds colored markup to the model showing
//! backbone validation results.
//!
//! @param imol Model molecule number  
//! @param cablam_log_file_name Path to CaBLAM output file
//!
//! @return PyObject* - List of tuples [(residue_spec, score), ...]
//!         where residue_spec is [chain_id, resno, ins_code]
//!         and score is the CaBLAM validation score
//!
//! Example usage:
//! \code{.py}
//! # Read CaBLAM results and add markup
//! results = add_cablam_markup_py(1, "cablam_results.txt")
//!
//! # Find worst CaBLAM outliers
//! outliers = [r for r in results if r[1] < 0.5]
//! print(f"Found {len(outliers)} CaBLAM outliers")
//!
//! for (chain, resno, ins), score in outliers:
//!     print(f"  {chain} {resno}: score = {score:.3f}")
//! \endcode
PyObject *add_cablam_markup_py(int imol, const std::string &cablam_log_file_name);
#endif

/*  ---------------------------------------------------------------------- */
/*                       key bindings :                                    */
/*  ---------------------------------------------------------------------- */
//! \name Key Bindings
//! \{
void print_key_bindings();
//! \}


/*  ---------------------------------------------------------------------- */
/*                       go to atom   :                                    */
/*  ---------------------------------------------------------------------- */

// not for export
void set_rotation_centre(const clipper::Coord_orth &pos);

#ifdef USE_GUILE
//
//! \brief Navigate to the next atom in the sequence (Guile interface)
//!
//! Given the current atom position, return the specification for the next atom
//! in the sequence. This traverses atoms in order: within a residue, then to
//! the next residue, then to the next chain.
//!
//! @param chain_id Current chain identifier
//! @param resno Current residue number
//! @param ins_code Current insertion code
//! @param atom_name Current atom name
//!
//! @return SCM - List containing the next atom specification [chain_id, resno, ins_code, atom_name]
//!
//! \note Returns the current position if already at the last atom
//!
//! Example usage:
//! \code{.scm}
//! (define next-atom (goto-next-atom-maybe "A" 42 "" "CA"))
//! ;; next-atom will be something like ("A" 42 "" "C") or ("A" 43 "" "N")
//! \endcode
SCM goto_next_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code, const char *atom_name);

//! \brief Navigate to the previous atom in the sequence (Guile interface)
//!
//! Given the current atom position, return the specification for the previous atom
//! in the sequence. This traverses atoms in reverse order: within a residue, then to
//! the previous residue, then to the previous chain.
//!
//! @param chain_id Current chain identifier
//! @param resno Current residue number
//! @param ins_code Current insertion code
//! @param atom_name Current atom name
//!
//! @return SCM - List containing the previous atom specification [chain_id, resno, ins_code, atom_name]
//!
//! \note Returns the current position if already at the first atom
SCM goto_prev_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif

#ifdef USE_PYTHON
//! \brief Navigate to the next atom in the sequence (Python interface)
//!
//! Given the current atom position, return the specification for the next atom
//! in the sequence. Useful for iterating through atoms programmatically.
//!
//! @param chain_id Current chain identifier
//! @param resno Current residue number
//! @param ins_code Current insertion code (use "" if none)
//! @param atom_name Current atom name
//!
//! @return PyObject* - List [chain_id, resno, ins_code, atom_name] for next atom
//!
//! Example usage:
//! \code{.py}
//! # Navigate forward through atoms
//! current = ["A", 42, "", "CA"]
//! next_atom = goto_next_atom_maybe_py(*current)
//! print(f"Next atom: {next_atom}")
//! \endcode
PyObject *goto_next_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);

//! \brief Navigate to the previous atom in the sequence (Python interface)
//!
//! Given the current atom position, return the specification for the previous atom.
//! Useful for iterating backwards through atoms programmatically.
//!
//! @param chain_id Current chain identifier
//! @param resno Current residue number
//! @param ins_code Current insertion code (use "" if none)
//! @param atom_name Current atom name
//!
//! @return PyObject* - List [chain_id, resno, ins_code, atom_name] for previous atom
PyObject *goto_prev_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif

int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec);
int set_go_to_atom_from_res_spec(const coot::residue_spec_t &spec);
#ifdef USE_GUILE
int set_go_to_atom_from_res_spec_scm(SCM residue_spec);
int set_go_to_atom_from_atom_spec_scm(SCM residue_spec);
#endif
#ifdef USE_PYTHON
int set_go_to_atom_from_res_spec_py(PyObject *residue_spec);
int set_go_to_atom_from_atom_spec_py(PyObject *residue_spec);
#endif


// This is to make porting the active atom more easy for Bernhard.
// Return a class rather than a list, and rewrite the active-residue
// function use this atom-spec.  The first value of the pair indicates
// if an atom spec was found.
//
std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec();
#ifdef USE_PYTHON
//! \brief Get the currently active atom (Python interface)
//!
//! Returns the atom specification for the atom at the current rotation center.
//!
//! @return PyObject* - Tuple (found, (molecule_number, atom_spec)) where:
//!         - found: Boolean indicating if an atom exists at center
//!         - molecule_number: Integer molecule ID
//!         - atom_spec: List [chain_id, resno, ins_code, atom_name, alt_conf]
//!
//! Example usage:
//! \code{.py}
//! found, (imol, atom_spec) = active_atom_spec_py()
//! if found:
//!     chain, resno, ins, atom, alt = atom_spec
//!     print(f"Active: {chain} {resno} {atom} in molecule {imol}")
//! else:
//!     print("No active atom")
//! \endcode
PyObject *active_atom_spec_py();
#endif // USE_PYTHON


/*  ---------------------------------------------------------------------- */
/*                       symmetry functions:                               */
/*  ---------------------------------------------------------------------- */
// get the symmetry operators strings for the given molecule
//
#ifdef USE_GUILE

//! \name More Scheme Symmetry Functions
//! \{

//! \brief return the symmetry of the imolth molecule
//!
//!   Return as a list of strings the symmetry operators of the
//!   given molecule. If imol is a not a valid molecule, return an empty
//!   list.*/
SCM get_symmetry(int imol);
//! \}
#endif // USE_GUILE

#ifdef USE_PYTHON
//! \name More Python Symmetry Functions
//! \{

//! \brief return the symmetry of the imolth molecule
//!
//!   Return as a list of strings the symmetry operators of the
//!   given molecule. If imol is a not a valid molecule, return an empty
//!   list.*/
//! @return a python object as a list (or some other python container)
PyObject *get_symmetry_py(int imol);
//! \}

#endif // USE_PYTHON

//! \name More Symmetry Functions
//! \{

//! \brief return 1 if this residue clashes with the symmetry-related
//!  atoms of the same molecule.
//!
//! @param imol the molecule index
//! @param chain_id the chain id
//! @param res_no the residue number
//! @param ins_code the insertion code
//! @param clash_dist the clash distance cut-off - typically 3.6 Angstroms
//!
//! @return 0 means that it did not clash,
//! -1 means that the residue or molecule could not be found or that there
//!    was no cell and symmetry.
int clashes_with_symmetry(int imol, const char *chain_id, int res_no, const char *ins_code,
                          float clash_dist);

//! Add molecular symmetry
//!
//! You will need to know how to expand your point group molecular symmetry
//! to a set of 3x3 matrices. Call this function for every matrix.
//! @param imol the molecule index
void add_molecular_symmetry(int imol,
                            double r_00, double r_01, double r_02,
                            double r_10, double r_11, double r_12,
                            double r_20, double r_21, double r_22,
                            double about_origin_x,
                            double about_origin_y,
                            double about_origin_z);

//! Add molecular symmetry from MTRIX records from file
//!
//! Often molecular symmetry is descibed using MTRIX card in a PDB file header.
//! Use this function to extract and apply such molecular symmmetry
int add_molecular_symmetry_from_mtrix_from_file(int imol, const std::string &file_name);

//! Add molecular symmetry
//!
//! This is a convenience function for the above - where you don't need to
//! specify the PDB file name.
//!
//! @param imol the molecule index
//! @return 1 on success and 0 on failure
int add_molecular_symmetry_from_mtrix_from_self_file(int imol);

//! \}

/*  ---------------------------------------------------------------------- */
/*                       map functions:                                    */
/*  ---------------------------------------------------------------------- */
//! \name Extra Map Functions
//! \{

/*! \brief read MTZ file filename and from it try to make maps

Useful for reading the output of refmac.  The default labels (FWT/PHWT
and DELFWT/PHDELFWT) can be changed using ...[something]

 @return a list of molecule numbers for the new maps
*/
std::vector<int> auto_read_make_and_draw_maps(const char *filename);
/*! \brief set the flag to do a difference map (too) on auto-read MTZ */
std::vector<int> auto_read_make_and_draw_maps_from_mtz(const std::string &file_name);
std::vector<int> auto_read_make_and_draw_maps_from_cns(const std::string &file_name);


//! \brief does the mtz file have the columms that we want it to have?
//!
//! @param mtz_file_name the mtz file name
//! @param f_col desired f_col
//! @param phi_col desired phi_col
//! @param weight_col desired weight col
//! @param mtz_file_name use_weights_flag specifies if the weight_col
//!        flay should be used/is interesting
int valid_labels(const std::string &mtz_file_name, const std::string &f_col,
		 const std::string &phi_col,
		 const std::string &weight_col,
		 bool use_weights_flag);

/* ----- remove wiget functions from this header GTK-FIXME
void add_map_colour_mol_menu_item(int imol, const std::string &name,
				  GtkWidget *sub_menu, GtkSignalFunc callback);
void add_map_scroll_wheel_mol_menu_item(int imol,
					const std::string &name,
					GtkWidget *menu,
					GtkSignalFunc callback);
*/

//! \brief make a sharpened or blurred map
//!
//! blurred maps are generated by using a positive value of b_factor.
//!
//! @param imol the molecule index
//! @param b_factor is the B-factor to blur by (positive numbers blur)
//! @return the index of the map created by applying a b-factor
//!        to the given map. Return -1 on failure.
int sharpen_blur_map(int imol_map, float b_factor);

//! \brief make a sharpened or blurred map with resampling
//!
//! resampling factor might typically be 1.3
//!
//! blurred maps are generated by using a positive value of b_factor.
//!
//! @param imol the molecule index
//! @param b_factor is the B-factor to blur by (positive numbers blur)
//! @return the index of the map created by applying a b-factor
//!        to the given map. Return -1 on failure.
int sharpen_blur_map_with_resampling(int imol_map, float b_factor, float resample_factor);

// This (gui function) allows a progress bar, and should not be part of the documented API
void sharpen_blur_map_with_resampling_threaded_version(int imol_map, float b_factor, float resample_factor);

#ifdef USE_GUILE
//! \brief make many sharpened or blurred maps
//!
//! blurred maps are generated by using a positive value of b_factor.
//!
void multi_sharpen_blur_map_scm(int imol_map, SCM b_factors_list);
#endif

#ifdef USE_PYTHON
//! \brief make many sharpened or blurred maps
//!
//! blurred maps are generated by using a positive value of b_factor.
//!
//! @param imol the molecule index
//! @param b_factor_lisst is a list of B-factor to blur by (positive numbers blur)
void multi_sharpen_blur_map_py(int imol_map, PyObject *b_factors_list);
#endif

#ifdef USE_PYTHON
//! amplitude vs resolution data for graph
//!
//! @param imol the molecule index
//! @return a list of list of lists, for each item in the list the resolution is the 0th element (in reciprocal Angstroms squared)
//!  and the count is the element index 1 and the amplitude is in element index 2.
PyObject *amplitude_vs_resolution_py(int mol_map);
#endif

#ifdef USE_GUILE
//! amplitude vs resolution data for graph
//! @param imol the molecule index
//! @return a list of list pairs, resolution is the first element (in reciprocal Angstroms squared) and the amplitude is in the second element.
SCM amplitude_vs_resolution_scm(int mol_map);
#endif

//! \brief Flip the hand of the map
//!
//! in case it was accidentally generated on the wrong one.
//! @return the molecule number of the flipped map.
int flip_hand(int imol_map);

#ifndef SWIG
//! test function for analysis of multiple map
int analyse_map_point_density_change(const std::vector<int> &map_number_list, int imol_map_mask);
#endif

#ifdef USE_PYTHON
//! I don't know what this does
//!
//! Looks interesting though
//! @param map_number_list a list of map numbers
//! @param imol_map_mask the molecule index for the mask
//! @return 0 on failure, 1 on success
int analyse_map_point_density_change_py(PyObject *map_number_list, int imol_map_mask);
#endif

//! \brief Go to the centre of the molecule - for Cryo-EM Molecules
//!
//! and recontour at a sensible value.
//!
//! @param imol the molecule index
void go_to_map_molecule_centre(int imol_map);

//! \brief b-factor from map
//!
//! calculate structure factors and use the amplitudes to estimate
//! the B-factor of the data using a wilson plot using a low resolution
//! limit of 4.5A.
//! @return -1 when given a bad map or there were no data beyond 4.5A
//!
float b_factor_from_map(int imol_map);


#ifdef USE_GUILE
//! \brief return the colour triple of the imolth map
//!
//! (e.g.: (list 0.4 0.6 0.8). If invalid imol return scheme false.
//!
SCM map_colour_components(int imol);
#endif // GUILE

#ifdef USE_PYTHON
//! \brief return the colour triple of the imolth map
//
//! e.g.: [0.4, 0.6, 0.8]. If invalid imol return Py_False.
//!
//! @param imol the molecule index
//! @return the colour triple of the imolth map
PyObject *map_colour_components_py(int imol);
#endif // PYTHON

//! \brief read a CCP4 map or a CNS map (despite the name)
//!
//! @param filename is the file name
//! @param is_diff_map_flag is either 0 or 1 denoting if this is a
//!        difference map
//! @return the molecule index of the new map. Return -1 on failure
int read_ccp4_map(const std::string &filename, int is_diff_map_flag);

//! \brief same function as above - old name for the function. Deleted from the API at some stage
//! @param filename is the file name
//! @param is_diff_map_flag is either 0 or 1 denoting if this is a
//!        difference map
//! @return the molecule index of the new map. Return -1 on failure
int handle_read_ccp4_map(const std::string &filename, int is_diff_map_flag);

//! \brief this reads a EMDB bundle - I don't think they exist any more
int handle_read_emdb_data(const std::string &dir_name);

void show_map_partition_by_chain_dialog();

//! Use the function for scriptng
std::vector<int> map_partition_by_chain(int imol_map, int imol_model);

//! Use the function for use in the GUI (non-blocking, no results returned)
void map_partition_by_chain_threaded(int imol_map, int imol_model);

//! use (or not) vertex gradients for the specified map
//!
//! vertex gradients make the map look smoother but are slower
//! to calculate
//!
//! @param imol the molecule index
//! @param state 0 for no, 1 for yes
void set_use_vertex_gradients_for_map_normals(int imol, int state);

//! the map should be displayed and not a difference map
void use_vertex_gradients_for_map_normals_for_latest_map();

//! alias for the above (more canonical naming)
void set_use_vertex_gradients_for_map_normals_for_latest_map();


//! \}

#ifdef SWIG
#else

// non-SWIGable functions:

/*! \brief
We overwrite the imol_map and we also presume that the
grid sampling of the contributing maps match. This makes it
much faster to generate than an average map.
*/
void regen_map_internal(int imol_map, const std::vector<std::pair<int, float> > &weighted_map_indices);

// As above, but also create a new map
int make_weighted_map_simple_internal(const std::vector<std::pair<int, float> > &weighted_map_indices);
#endif

//! \brief colour map by other map
//!
//! maybe we need to spefify_other things like the colour table.
//! @param imol_map the molecule index
//! @param imol_map_used_for_colouring is the other map index
void
colour_map_by_other_map(int imol_map, int imol_map_used_for_colouring);

#ifdef USE_PYTHON
//! \brief  the colour_table should be a list of colours
//!
//! So, if the map has 4 entries covering the range from  0 to 1, then the table_bin_size would be 0.25
//! and the colour_table list would have 4 entries covering the range 0->0.25, 0.25->0.5, 0.5->0.75, 0.75->1.0
//!
void
colour_map_by_other_map_py(int imol_map, int imol_map_used_for_colouring, float table_bin_start, float table_bin_size,
                           PyObject *colour_table_list);

PyObject *export_molecule_as_x3d(int imol);

#endif

bool export_molecule_as_obj(int imol, const std::string &file_name);

// later this will need a handle as extra arg.
bool export_molecule_as_gltf(int imol, const std::string &file_name);

//! \brief  turn of colour map by other map
void colour_map_by_other_map_turn_off(int imol_map);

//! \brief Add map caps
void add_density_map_cap();

//! \brief colour meshes (e.g. Ribbon diagrams) by map
//!
//! scale might be 2 and offset 1 (for example)
void recolour_mesh_by_map(int imol_model, int imol_map, float scale, float offset);


//! \name Multi-Residue Torsion
//! \{
#ifdef USE_GUILE
/*! \brief fit residues

(note: fit to the current-refinement map)
*/
SCM multi_residue_torsion_fit_scm(int imol, SCM residues_specs_scm, int n_trials);
#endif // GUILE
void multi_residue_torsion_fit(int imol, const std::vector<coot::residue_spec_t> &specs, int n_trials);

#ifdef USE_PYTHON
/*! \brief fit residues

(note: fit to the current-refinement map)
*/
PyObject *multi_residue_torsion_fit_py(int imol, PyObject *residues_specs_py, int n_trials);
#endif // PYTHON
//! \}


// Where should this go?
void import_bild(const std::string &file_name);

// resolution in A.
void servalcat_fofc(int imol_model,
                    int imol_fofc_map, const std::string &half_map_1, const std::string &half_map_2,
                    float resolution);

//! resolution in A.
//!
void servalcat_refine(int imol_model,
                      const std::string &half_map_1, const std::string &half_map_2,
                      const std::string &mask_map, float resolution);

//! run acedrg link
void
run_acedrg_link_generation(const std::string &acedrg_link_command);

//! run generic process - doesn't work at the moment - on_completion_args
//! is wrongly interpretted.
void add_toolbar_subprocess_button(const std::string &button_label,
                                   const std::string &subprocess_command,
                                   PyObject *arg_list,
                                   PyObject *on_completion_function,
                                   PyObject *on_completion_args);


/*  ------------------------------------------------------------------------ */
/*                             Add an Atom                                   */
/*  ------------------------------------------------------------------------ */
//! \name Add an Atom
//! \{
//! \brief add an atom. element can be "Water", "Na", "K", "I" etc.
void add_an_atom(const std::string &element);
//! \}

/*  ------------------------------------------------------------------------ */
/*                             Nudge the B-factors                           */
/*  ------------------------------------------------------------------------ */
//! \name Nudge the B-factors
//! \{
//! \brief change the b-factors of the specified by a (small) amount
#ifdef USE_PYTHON
void nudge_the_temperature_factors_py(int imol, PyObject *residue_spec_py, float amount);
#endif
//! \}


/*  ------------------------------------------------------------------------ */
/*                         merge fragments                                   */
/*  ------------------------------------------------------------------------ */
//! \name Merge Fragments
//! \{
//! \brief merge fragments
//
//! each fragment is presumed to be in its own chain.
//
int merge_fragments(int imol);
//! \}

/*  ------------------------------------------------------------------------ */
/*                         delete items                                      */
/*  ------------------------------------------------------------------------ */
//! \name Delete Items
//! \{
//! \brief Delete Items

/*! \brief delete the chain  */
void delete_chain(int imol, const std::string  &chain_id);

/*! \brief delete the side chains in the chain  */
void delete_sidechains_for_chain(int imol, const std::string &chain_id);

//! \}


/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */
//! \name Execute Refmac
//! \{
//! \brief execute refmac
//
//! if swap_map_colours_post_refmac_flag is not 1 thenn imol_refmac_map is ignored.
//
void
execute_refmac_real(std::string pdb_in_filename,
                    std::string pdb_out_filename,
                    std::string mtz_in_filename,
                    std::string mtz_out_filename,
                    std::string cif_lib_filename, /* use "" for none */
                    std::string fobs_col_name,
                    std::string sigfobs_col_name,
                    std::string r_free_col_name,
                    short int have_sensible_free_r_flag,
                    short int make_molecules_flag,
                    std::string refmac_count_string,
                    int swap_map_colours_post_refmac_flag,
                    int imol_refmac_map,
                    int diff_map_flag,
                    int phase_combine_flag,
                    std::string phib_string,
                    std::string fom_string,
                    std::string ccp4i_project_dir);

/*! \brief the name for refmac

 @return a stub name used in the construction of filename for refmac output */
std::string refmac_name(int imol);

//! \}


/*  ------------------------------------------------------------------- */
/*                    file selection                                    */
/*  ------------------------------------------------------------------- */

namespace coot {
   //! str mtime for for attributes
   class str_mtime {
   public:
      str_mtime(std::string file_in, time_t mtime_in) {
         mtime = mtime_in;
         file = file_in;
      }
      str_mtime() {}
      time_t mtime;
      std::string file;
   };

   //! trivial helper function for file attributes
   class file_attribs_info_t {
   public:
      std::string directory_prefix;
      std::vector<str_mtime> file_mtimes;
   };
}


bool compare_mtimes(coot::str_mtime a, coot::str_mtime b);

std::vector<std::pair<std::string, std::string> > parse_ccp4i_defs(const std::string &filename);

std::string ccp4_project_directory(const std::string &ccp4_project_name);

/*  -------------------------------------------------------------------- */
/*                     history                                           */
/*  -------------------------------------------------------------------- */
#include "command-arg.hh"

void add_to_history(const std::vector<std::string> &ls);
void add_to_history_simple(const std::string &cmd);
void add_to_history_typed(const std::string &command,
                          const std::vector<coot::command_arg_t> &args);
std::string single_quote(const std::string &s);
std::string pythonize_command_name(const std::string &s);
std::string schemize_command_name(const std::string &s);
std::string languagize_command(const std::vector<std::string> &command_parts);

void add_to_database(const std::vector<std::string> &command_strings);


/*  ----------------------------------------------------------------------- */
/*                         Merge Molecules                                  */
/*  ----------------------------------------------------------------------- */
#include "api/merge-molecule-results-info-t.hh"
// return the status and vector of chain-ids of the new chain ids.
//
std::pair<int, std::vector<merge_molecule_results_info_t> > merge_molecules_by_vector(const std::vector<int> &add_molecules, int imol);

/*  ----------------------------------------------------------------------- */
/*                         Dictionaries                                     */
/*  ----------------------------------------------------------------------- */
//! \name Dictionary Functions
//! \{

/*                  cif (geometry) dictionary                            */
/* \brief return the number of bonds read (> 0 can be treated as success) */
int handle_cif_dictionary(const std::string &filename);
/* \brief synonym for above.

return the number of bonds read (> 0 can be treated as success) */
int read_cif_dictionary(const std::string &filename);

/* \brief return the number of bonds read (> 0 can be treated as success).
 Apply to the given molecule.

 imol_enc can be the model molecule number or
 IMOL_ENC_ANY = -999999, IMOL_ENC_AUTO = -999998, IMOL_ENC_UNSET = -999997

 */
int handle_cif_dictionary_for_molecule(const std::string &filename, int imol_enc, short int new_molecule_from_dictionary_cif_checkbutton_state);

//! dictionary entries
std::vector<std::string> dictionary_entries();

//! debug dictionary information
void debug_dictionary();

//! get types in molecule
//!
//! @param imol the molecule index
//! @return a vector of residue types
std::vector<std::string> get_types_in_molecule(int imol);

//! Get the SMILES for the given residue type
//!
//! @param comp_id is the residue type
//! @return the SMILES string
std::string SMILES_for_comp_id(const std::string &comp_id);

/*! \brief return a list of all the dictionaries read */
#ifdef USE_GUILE
SCM dictionaries_read();
SCM cif_file_for_comp_id_scm(const std::string &comp_id);
SCM dictionary_entries_scm();
SCM SMILES_for_comp_id_scm(const std::string &comp_id);
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *dictionaries_read_py();
PyObject *cif_file_for_comp_id_py(const std::string &comp_id);
PyObject *dictionary_entries_py();

//! Get the SMILES for the given residue type
//
//! @param comp_id is the residue type
//! @return the SMILES string or False on failure to find the
//!         residue type or SMILES string
PyObject *SMILES_for_comp_id_py(const std::string &comp_id);

#endif // PYTHON
//! \}


/*  ----------------------------------------------------------------------- */
/*                         Restraints                                       */
/*  ----------------------------------------------------------------------- */
//! \name  Restraints Interface
/// \{

#ifdef USE_GUILE
//! \brief return the monomer restraints for the given monomer_type,
//!       return scheme false on "restraints for monomer not found"
SCM monomer_restraints(const char *monomer_type);

//! \brief set the monomer restraints of the given monomer_type
//!
//!  @return scheme false or true for success or failure to set the
//!  restrains for monomer_type */
SCM set_monomer_restraints(const char *monomer_type, SCM restraints);
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *monomer_restraints_py(std::string monomer_type);
PyObject *monomer_restraints_for_molecule_py(std::string monomer_type, int imol);
PyObject *set_monomer_restraints_py(const char *monomer_type, PyObject *restraints);
#endif // USE_PYTHON

/*! \brief show restraints editor */
void show_restraints_editor(std::string monomer_type);

void show_restraints_editor_by_index(int menu_item_index);

/*! \brief write cif restraints for monomer */
void write_restraints_cif_dictionary(std::string monomer_type, std::string file_name);

//! \}

/*  ----------------------------------------------------------------------- */
/*                      list nomenclature errors                            */
/*  ----------------------------------------------------------------------- */
std::vector<std::pair<std::string, coot::residue_spec_t> >
list_nomenclature_errors(int imol);

#ifdef USE_GUILE
SCM list_nomenclature_errors_scm(int imol);
#endif // USE_GUILE
#ifdef USE_PYTHON
PyObject *list_nomenclature_errors_py(int imol);
#endif // USE_PYTHON

void
show_fix_nomenclature_errors_gui(int imol,
                                 const std::vector<std::pair<std::string, coot::residue_spec_t> > &nomenclature_errors);

/*  ----------------------------------------------------------------------- */
/*                  dipole                                                  */
/*  ----------------------------------------------------------------------- */
// This is here because it uses a C++ class, coot::dipole
//
#ifdef USE_GUILE
SCM dipole_to_scm(std::pair<coot::dipole, int> dp);
#endif // USE_GUILE
#ifdef USE_PYTHON
PyObject *dipole_to_py(std::pair<coot::dipole, int> dp);
#endif // USE_PYTHON


#ifdef USE_PYTHON
PyObject *coot_has_guile();
#endif

bool coot_can_do_lidia_p();


/* commands to run python commands from guile and vice versa */
/* we ignore return values for now */
#ifdef USE_PYTHON
PyObject *run_scheme_command(const char *scheme_command);
#endif // USE_PYTHON
#ifdef USE_GUILE
SCM run_python_command(const char *python_command);
#endif // USE_GUILE

// This is not inside a #ifdef USE_PYTHON because we want to use it
// from the guile level and USE_PYTHON is not passed as an argument to
// swig when generating coot_wrap_guile.cc.
//
// [Consider removing safe_python_command_by_char_star() which is
// conditionally compiled].
int pyrun_simple_string(const char *python_command);

#ifdef USE_GUILE
// Return a list describing a residue like that returned by
// residues-matching-criteria (list return-val chain-id resno ins-code)
// This is a library function really.  There should be somewhere else to put it.
// It doesn't need expression at the scripting level.
// return a null list on problem
SCM residue_spec_to_scm(const coot::residue_spec_t &res);
#endif

#ifdef USE_PYTHON
// Return a list describing a residue like that returned by
// residues-matching-criteria [return_val, chain_id, resno, ins_code]
// This is a library function really.  There should be somewhere else to put it.
// It doesn't need expression at the scripting level.
// return a null list on problem
PyObject *residue_spec_to_py(const coot::residue_spec_t &res);
#endif

#ifdef USE_PYTHON
PyObject *residue_spec_make_triple_py(PyObject *residue_spec_py);
#endif // USE_PYTHON

#ifdef USE_GUILE
coot::residue_spec_t residue_spec_from_scm(SCM residue_in);
#endif

#ifdef USE_PYTHON
coot::residue_spec_t residue_spec_from_py(PyObject *residue_in);
#endif

// return a spec for the first residue with the given type.
// test the returned spec for unset_p().
//
coot::residue_spec_t get_residue_by_type(int imol, const std::string &residue_type);

std::vector<coot::residue_spec_t> get_residue_specs_in_mol(int imol, const std::string &residue_type);

#ifdef USE_PYTHON
// Always returns a list
PyObject *get_residue_specs_in_mol_py(int imol, const std::string &residue_type);
#endif

#ifdef USE_GUILE
// return a residue spec or scheme false
SCM get_residue_by_type_scm(int, const std::string &residue_type);
#endif
#ifdef USE_PYTHON
// return a residue spec or Python False.
PyObject *get_residue_by_type_py(int, const std::string &residue_type);
#endif


/*  ----------------------------------------------------------------------- */
/*               Atom info                                                  */
/*  ----------------------------------------------------------------------- */

//! \name Atom Information functions
//! \{

#ifdef USE_GUILE
//! \brief output atom info in a scheme list for use in scripting
//!
//! in this format (list occ temp-factor element x y z).  Return empty
//! list if atom not found. */
SCM atom_info_string_scm(int imol, const char *chain_id, int resno,
                         const char *ins_code, const char *atname,
                         const char *altconf);
SCM molecule_to_pdb_string_scm(int imol);
#endif // USE_GUILE

/*! \brief return the rename from a residue serial number

   @return blank ("") on failure. */
std::string resname_from_serial_number(int imol, const char *chain_id, int serial_num);

//! \brief return the residue name of the specified residue
std::string residue_name(int imol, const std::string &chain_id, int resno, const std::string &ins_code);

//! \brief return the serial number of the specified residue
//!
//! @return -1 on failure to find the residue
//
int serial_number_from_residue_specs(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);


#ifdef USE_GUILE
//! \brief Return a list of atom info for each atom in the specified residue.
//!
//! output is like this:
//! (list
//!    (list (list atom-name alt-conf)
//!          (list occ temp-fact element)
//!          (list x y z)))
//!
//! occ can be a single number or a list of seven numbers of which the first is
//! the isotropic B.
//!
SCM residue_info(int imol, const char* chain_id, int resno, const char *ins_code);
SCM residue_name_scm(int imol, const char* chain_id, int resno, const char *ins_code);

//! chain fragments
SCM chain_fragments_scm(int imol, short int screen_output_also);

//! \brief generate a molecule from an s-expression
//!
//! return a molecule number, -1 on error
int add_molecule(SCM molecule_expression, const char *name);

//! \brief update a molecule from a s-expression
//!
//! And going the other way, given an s-expression, update
//! molecule_number by the given molecule.  Clear what's currently
//! there first though.
//!
int clear_and_update_molecule(int molecule_number, SCM molecule_expression);

//! \brief return specs of the atom close to screen centre
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf) for atom that is closest to the screen centre in any
//! displayed molecule.  If there are multiple models with the same
//! coordinates at the screen centre, return the attributes of the atom
//! in the highest number molecule number.
//!
//! return scheme false if no active residue
//!
SCM active_residue();

//! \brief return the specs of the closest displayed atom
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf (list x y z)) for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, potential CA
//! substition is not performed).  If there is no atom, or if imol is
//! not a valid model molecule, return scheme false.
//!
SCM closest_atom_simple_scm();

//! \brief return the specs of the closest atom in imolth molecule
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf (list x y z)) for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, no account is
//! taken of the displayed state of the molecule).  If there is no
//! atom, or if imol is not a valid model molecule, return scheme false.
//!
SCM closest_atom(int imol);

//! \brief return the specs of the closest atom to the centre of the screen
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf (list x y z)) for atom that is closest to the screen
//! for displayed molecules. If there is no atom, return scheme false.
//! Don't choose the CA of the residue if there is a CA in the residue
//! of the closest atom.
//! 201602015-PE: I add this now, but I have a feeling that I've done this
//! before.
SCM closest_atom_raw_scm();

//! \brief return residues near residue
//!
//! Return residue specs for residues that have atoms that are
//! closer than radius Angstroems to any atom in the residue
//! specified by res_in.
//!
SCM residues_near_residue(int imol, SCM residue_in_scm, float radius);

//! \brief return residues near the given residues
//!
//! Return residue specs for residues that have atoms that are
//! closer than radius Angstroems to any atom in the residue
//! specified by res_in.
//!
SCM residues_near_residues_scm(int imol, SCM residues_in, float radius);

//! \brief residues near residue
//!
//! @return residues within radius of pos (x,y,z) position
//!
//! Return a list of pairs of (imol, residue_spec).
//! pos is a list of 3 numbers.  (get imol from active-atom)
//!
SCM residues_near_position_scm(int imol, SCM pos, float radius);

//! \brief label the closest atoms in the residues that neighbour residue_spec
//!
void label_closest_atoms_in_neighbour_residues_scm(int imol, SCM residue_spec_scm, float radius);

#endif        /* USE_GUILE */

//! \brief find the active residue, find the near residues (within radius)
//! create a new molecule, run reduce on that, import hydrogens from
//! the result and apply them to the molecule of the active residue.
void hydrogenate_region(float radius);

//! \brief Add hydrogens to imol from the given pdb file
void add_hydrogens_from_file(int imol, std::string pdb_with_Hs_file_name);

//! \brief add hydrogen atoms to the specified residue
void add_hydrogen_atoms_to_residue(int imol, std::string chain_id, int res_no, std::string ins_code);

#ifdef USE_PYTHON
//! \brief add hydrogen atoms to the specified residue
void add_hydrogen_atoms_to_residue_py(int imol, PyObject *residue_spec_py);
#endif

/* Here the Python code for ATOM INFO */

//! \brief output atom info in a python list for use in scripting:
//!
//! in this format [occ, temp_factor, element, x, y, z].  Return empty
//! list if atom not found. */
#ifdef USE_PYTHON
PyObject *atom_info_string_py(int imol, const char *chain_id, int resno,
                              const char *ins_code, const char *atname,
                              const char *altconf);

//! \brief
//!
//! Return the molecule as a PDB string
PyObject *molecule_to_pdb_string_py(int imol);

//! \brief Get detailed atom information for a residue (Python interface)
//!
//! Returns per-atom information including coordinates, occupancy, B-factor,
//! and element for all atoms in the specified residue. Useful for inspecting
//! residue completeness and identifying missing atoms.
//!
//! @param imol Model molecule index
//! @param chain_id Chain identifier (e.g., "A")
//! @param resno Residue number
//! @param ins_code Insertion code (use "" if none)
//!
//! @return PyObject* - A list of atom information, one entry per atom:
//!   \code
//!   [
//!     [[atom_name, alt_conf], [occupancy, b_factor, element, seg_id], [x, y, z], atom_index],
//!     ...
//!   ]
//!   \endcode
//!   - \c atom_name (str): Atom name (e.g., " CA ", " SG ")
//!   - \c alt_conf (str): Alternate conformation identifier ("" if none)
//!   - \c occupancy (float): Atom occupancy (0.0-1.0)
//!   - \c b_factor (float or list of [b_iso, B11, B22, B33, B12, B13, B23]): Temperature factor
//!   - \c element (str): Element symbol (e.g., " C", " N", " S")
//!   - \c x, \c y, \c z (float): Cartesian coordinates in Ångstroms
//!   - \c atom_index (int): Internal atom index
//!
//! Example usage:
//! \code{.py}
//! # Check if a CYS residue has all expected atoms
//! atoms = coot.residue_info_py(0, "A", 72, "")
//! atom_names = [a[0][0].strip() for a in atoms]
//! print(f"Atoms present: {atom_names}")
//!
//! expected_cys = ['N', 'CA', 'CB', 'SG', 'C', 'O']
//! missing = [a for a in expected_cys if a not in atom_names]
//! if missing:
//!     print(f"Missing atoms: {missing}")
//!
//! # Get B-factors for all atoms
//! for atom in atoms:
//!     name = atom[0][0].strip()
//!     b_factor = atom[1][1]
//!     print(f"{name}: B={b_factor:.2f}")
//! \endcode
PyObject *residue_info_py(int imol, const char* chain_id, int resno, const char *ins_code);

//! \brief get the residue name
//!
//! @param imol Model molecule index
//! @param chain_id Chain identifier (e.g., "A")
//! @param resno Residue number
//! @param ins_code Insertion code (use "" if none)
//! @return residue name string or blank string on failure
//!
PyObject *residue_name_py(int imol, const char* chain_id, int resno, const char *ins_code);

// the expanded form of this is in c-interface.h
PyObject *residue_centre_from_spec_py(int imol,
                                      PyObject *spec_py);

PyObject *chain_fragments_py(int imol, short int screen_output_also);

#ifdef USE_PYTHON
void set_b_factor_residues_py(int imol, PyObject *residue_specs_b_value_tuple_list_py);
#endif

#ifdef USE_GUILE
void set_b_factor_residues_scm(int imol, SCM residue_specs_b_value_tuple_list_scm);
#endif

//! \}

//! \name Using S-expression molecules
//! \{

// And going the other way, given an python-expression, update
// molecule_number by the given molecule.  Clear what's currently
// there first though.
//
int clear_and_update_molecule_py(int molecule_number, PyObject *molecule_expression);
// return a molecule number, -1 on error
int add_molecule_py(PyObject *molecule_expression, const char *name);

//! \brief
//! Return a list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf] for atom that is closest to the screen centre.  If there
//! are multiple models with the same coordinates at the screen centre,
//! return the attributes of the atom in the highest number molecule
//! number.
//
//! return False if no active residue
//
PyObject *active_residue_py();

//! \brief return the spec of the closest displayed atom
//!
//! @return a list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf, [x, y, z]] for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, potential CA
//! substition is not performed).  If there is no atom, or if imol is
//! not a valid model molecule, return False.
//!
PyObject *closest_atom_simple_py();

//! \brief return closest atom in imolth molecule
//!
//! @param imol is the molecule index
//! @return a list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf, [x, y, z]] for the atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, no account is
//! taken of the displayed state of the imol molecule).  If there is no
//! atom, or if imol is not a valid model molecule, return False.
//
PyObject *closest_atom_py(int imol);

//! \brief return the specs of the closest atom to the centre of the screen
//!
//! @return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf (list x y z)) for atom that is closest to the screen
//! for displayed molecules. If there is no atom, return scheme false.
//! Don't choose the CA of the residue if there is a CA in the residue
//! of the closest atom
PyObject *closest_atom_raw_py();


//! \brief get the residues near a specified residue
//!
//! This is useful to select the residues for "Sphere" refinement.
//!
//! @param imol is the molecule index
//! @param residue_in is a residue spec [chain_id, res_no, insertion_code]
//! @param radius is the cut-off distance atoms of the surrounding residues
//!        if they are to be included in the residue selection.
//! @return a list of residue specs for residues that have atoms that are
//! closer than radius Angstroems to any atom in the residue
//! specified by residue_in.
//
PyObject *residues_near_residue_py(int imol, PyObject *residue_spec_in, float radius);

//! \brief get the residues near a specified list of residues
//!
//! @param imol is the molecule index
//! @param residues_in is a list of residue specs each of which is
//!        [chain_id, res_no, insertion_code]
//! @param radius is the cut-off distance atoms of the surrounding residues
//!        if they are to be included in the residue selection.
//! @return a list of residue specs for residues that have atoms that are
//! closer than radius Angstroems to any atom in the residue
//! specified by the input residue spec list.
//
PyObject *residues_near_residues_py(int imol, PyObject *residues_specs_in, float radius);

//! \brief
//! Return residue specs for residues that have atoms that are
//! closer than radius Angstroems to the given position.
//!
PyObject *residues_near_position_py(int imol, PyObject *pos_in, float radius);

//! \brief label the closest atoms in the residues that neighbour residue_spec
//!
void label_closest_atoms_in_neighbour_residues_py(int imol, PyObject *residue_spec_py, float radius);

//! \brief return a Python object for the bonds
//
PyObject *get_bonds_representation(int imol);

//! \brief replace any current non-drawn bonds with these - and regen bonds
void set_new_non_drawn_bonds(int imol, const std::string &cid);

//! \brief add to non-drawn bonds - and regen bonds
void add_to_non_drawn_bonds(int imol, const std::string &cid);

//! \brief clear the non-drawn bonds - force regen bonds to restore all
void clear_non_drawn_bonds(int imol);

//! \brief return a Python object for the radii of the atoms in the dictionary
//
PyObject *get_dictionary_radii();

//! \brief return a Python object for the representation of bump and hydrogen bonds oft
//          the specified residue
PyObject *get_environment_distances_representation_py(int imol, PyObject *residue_spec_py);

//! \brief return a Python object for the intermediate atoms bonds
//
PyObject *get_intermediate_atoms_bonds_representation();

#endif // USE_PYTHON

//! \brief return the continue-updating-refinement-atoms state
//
// 0 means off, 1 means on.
// Given the current wiring of the refinement, this is always 0, i.e. refine_residues()
// will return only after the atoms have finished moving.
int get_continue_updating_refinement_atoms_state();


//! \}

//! \name status bar string functions
//! \{
// status bar atom info text here?!
std::string atom_info_as_text_for_statusbar(int atom_index, int imol);
std::string atom_info_as_text_for_statusbar(int atom_index, int imol,
                                            const std::pair<symm_trans_t, Cell_Translation> &sts);
//! \}


/*  ----------------------------------------------------------------------- */
/*                  Refinement                                              */
/*  ----------------------------------------------------------------------- */

//! \name Refinement with specs
//! \{

//! \brief
//! a utility to return the specs of all the residues, each spec prefixed by the serial number
#ifdef USE_GUILE
SCM all_residues_with_serial_numbers_scm(int imol);
#endif
#ifdef USE_PYTHON
PyObject *all_residues_with_serial_numbers_py(int imol);
#endif


#ifdef SWIG
#else

//! \brief
//! regularize the given residues
//!
void regularize_residues(int imol, const std::vector<coot::residue_spec_t> &residues);
#endif

//! presumes that imol_Refinement_Map has been set
std::string mtz_file_name(int imol);

#ifdef USE_GUILE

//! \brief
//! Refine the given residue range
//!
SCM refine_zone_with_full_residue_spec_scm(int imol, const char *chain_id,
                                           int resno1,
                                           const char*inscode_1,
                                           int resno2,
                                           const char*inscode_2,
                                           const char *altconf);
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *refine_zone_with_full_residue_spec_py(int imol, const char *chain_id,
                                           int resno1,
                                           const char*inscode_1,
                                           int resno2,
                                           const char*inscode_2,
                                           const char *altconf);
#endif // USE_PYTHON

//! set display of rotamer markup during interactive real space refinement
void set_draw_moving_atoms_rota_markup(short int state);
//! set display of ramachandran markup during interactive real space refinement
void set_draw_moving_atoms_rama_markup(short int state);

//! the old names for the above functions:
void set_show_intermediate_atoms_rota_markup(short int state);
//! the old names for the above functions:
void set_show_intermediate_atoms_rama_markup(short int state);

//! the geters for the rota markup
int get_draw_moving_atoms_rota_markup_state();

//! the geters for the rama markup
int get_draw_moving_atoms_rama_markup_state();

//! the old names for the above functions:
int get_show_intermediate_atoms_rota_markup();

//! the old names for the above functions:
int get_show_intermediate_atoms_rama_markup();

void set_cryo_em_refinement(bool mode);
bool get_cryo_em_refinement();

#ifdef USE_GUILE
SCM accept_moving_atoms_scm();
#endif
#ifdef USE_PYTHON
//! Accept moving atoms
//!
//! \brief Accept refined/regularized atoms into the main molecule (Python interface)
//!
//! When scripting refinement with set_refinement_immediate_replacement(1),
//! call this function after refinement operations to ensure atoms are
//! committed. While immediate replacement mode should handle this
//! automatically, calling accept_moving_atoms_py() ensures reliable
//! synchronization.
//!
//! \return PyObject* with one of:
//!   - \c Py_False if no restraints were found (nothing to accept)
//!   - A Python list \c [info_text, progress, lights] on success:
//!     - \c info_text (str): Usually empty string
//!     - \c progress (int): GSL minimization status
//!       - 0 = GSL_SUCCESS (converged)
//!       - -2 = GSL_CONTINUE
//!       - 27 = GSL_ENOPROG (no progress)
//!     - \c lights (list): Refinement statistics as [[name, label, value], ...]
//!       - \c name (str): Restraint type (e.g., "Bonds", "Angles", 
//!         "Trans_peptide", "Planes", "Non-bonded", "Chirals")
//!       - \c label (str): Formatted string (e.g., "Bonds: 0.625")
//!       - \c value (float): Distortion value (lower is better)
//!
//! Example usage:
//! \code{.py}
//! coot.set_refinement_immediate_replacement(1)
//! coot.refine_residues_py(0, [["A", 42, ""]])
//! result = coot.accept_moving_atoms_py()
//!
//! if result:
//!     info, progress, lights = result
//!     for name, label, value in lights:
//!         print(f"{name}: {value:.3f}")
//! \endcode
PyObject *accept_moving_atoms_py();
#endif


#ifdef USE_PYTHON
void register_post_intermediate_atoms_moved_hook(PyObject *function_name);
#endif

void set_regenerate_bonds_needs_make_bonds_type_checked(bool state);
bool get_regenerate_bonds_needs_make_bonds_type_checked_state();

//! \}


/*  ----------------------------------------------------------------------- */
/*                  rigid body fitting (multiple residue ranges)            */
/*  ----------------------------------------------------------------------- */

//! \brief
//! return 0 on fail to refine (no sensible place to put atoms) and 1
//! on fitting happened.
int rigid_body_fit_with_residue_ranges(int imol, const std::vector<coot::high_res_residue_range_t> &ranges);

// Model morphing (average the atom shift by using shifts of the
// atoms within shift_average_radius A of the central residue).
//
// return 0 on fail to move atoms and 1 on fitting happened.
//
int morph_fit_all(int imol, float transformation_averaging_radius);

//! \brief Morph the given chain
int morph_fit_chain(int imol, std::string chain_id, float transformation_averaging_radius);
#ifdef USE_GUILE
int morph_fit_residues_scm(int imol, SCM residue_specs,       float transformation_averaging_radius);
#endif
#ifdef USE_PYTHON
int morph_fit_residues_py( int imol, PyObject *residue_specs, float transformation_averaging_radius);
#endif
//! \brief morph the given residues.
int morph_fit_residues(int imol, const std::vector<coot::residue_spec_t> &residue_specs,
                       float transformation_averaging_radius);

//! \brief morph transformation are based primarily on rigid body refinement
//! of the secondary structure elements.
int morph_fit_by_secondary_structure_elements(int imol, const std::string &chain_id);


/*  ----------------------------------------------------------------------- */
/*                  check water baddies                                     */
/*  ----------------------------------------------------------------------- */

std::vector<coot::atom_spec_t>
check_waters_baddies(int imol, float b_factor_lim, float map_sigma_lim, float min_dist, float max_dist, short int part_occ_contact_flag, short int zero_occ_flag, short int logical_operator_and_or_flag);

// blobs, returning position and volume
//
std::vector<std::pair<clipper::Coord_orth, double> >
find_blobs(int imol_model, int imol_map, float cut_off_density_level);

#ifdef USE_GUILE
//! \brief find blobs
SCM find_blobs_scm(int imol_model, int imol_map, float cut_off_density_level);
#endif

#ifdef USE_PYTHON
//! @brief Find regions of unmodeled electron density ("blobs") in a map
//! 
//! Identifies regions of significant electron density that are not explained by the current
//! atomic model. This is essential for discovering missing features such as waters, ligands,
//! alternative conformations, metal ions, or missing residues during structure validation
//! and refinement.
//! 
//! The function masks out density already explained by the model atoms, then searches for
//! contiguous regions of density above the specified sigma threshold. Each blob is 
//! characterized by its center position and an integrated volume/score representing the
//! strength of the feature.
//! 
//! @param imol_model  The model molecule index. Density explained by atoms in this model
//!                    will be masked out (excluded) from the search. Must be a valid model
//!                    molecule.
//! 
//! @param imol_map    The map molecule index to search for unmodeled density. This is 
//!                    typically a difference map (mFo-DFc) for most sensitive detection,
//!                    but can also be a regular map (2mFo-DFc). Must be a valid map molecule.
//! 
//! @param cut_off_sigma  The sigma threshold for blob detection (in units of map sigma).
//!                       Typical values:
//!                       - 3.5σ: Standard threshold for significant features in difference maps
//!                       - 2.5σ: More sensitive, finds weaker features (more false positives)
//!                       - 4.5σ: Conservative, only strong features (fewer false positives)
//!                       - 1.0σ: For regular maps (2mFo-DFc), lower threshold appropriate
//! 
//! @return PyObject* - Returns a Python list of blobs, or Py_False on error.
//! 
//! Return format (on success):
//! @code{.py}
//! [
//!   [[x1, y1, z1], volume1],  # First blob: [position_list, score]
//!   [[x2, y2, z2], volume2],  # Second blob
//!   ...
//! ]
//! @endcode
//! 
//! Each blob is represented as a 2-element list:
//! - Element 0: Position as [x, y, z] list (coordinates in Ångströms, orthogonal space)
//! - Element 1: Volume/score as float (integrated density strength)
//! 
//! Return value (on error):
//! - Py_False if imol_model is not a valid model molecule
//! - Py_False if imol_map is not a valid map molecule
//! 
//! @section interpretation Interpreting Results
//! 
//! Blob Score Interpretation (for difference maps at 3σ):
//! - **>50**: Large feature - likely missing ligand, cofactor, metal with waters, or several residues
//! - **10-50**: Medium feature - likely 1-3 water molecules or alternative conformation
//! - **3-10**: Small feature - likely single water or weak alternative conformation
//! - **<3**: Very weak feature - may be noise or very weak partial occupancy
//! 
//! Map Type Considerations:
//! - **Difference maps (mFo-DFc)**: Most sensitive for finding missing features
//!   - Positive blobs (>3σ): Missing atoms/features that should be added
//!   - Negative blobs (<-3σ): Incorrectly modeled atoms that should be removed/moved
//! - **Regular maps (2mFo-DFc)**: Less sensitive to model bias
//!   - Good for finding larger missing features
//!   - Use lower sigma threshold (0.5-1.5σ)
//! 
//! @section algorithm Algorithm Details
//! 
//! The blob detection algorithm:
//! 1. Imports the map specified by imol_map
//! 2. Masks the map using atoms from imol_model with a 1.9Å radius
//! 3. Waters can be optionally excluded from masking (controlled by find_ligand_mask_waters_flag)
//! 4. Searches for contiguous density regions above cut_off_sigma threshold
//! 5. Calculates center position and integrated volume for each blob
//! 6. Returns blobs sorted by significance
//! 
//! @section example Example Usage
//! 
//! @code{.py}
//! # Find blobs in difference map at 3 sigma
//! blobs = coot.find_blobs_py(
//!     imol_model=0,              # protein model
//!     imol_map=2,                # difference map (mFo-DFc)
//!     cut_off_sigma=3.0          # standard threshold
//! )
//! 
//! if blobs:
//!     print(f"Found {len(blobs)} blobs")
//!     
//!     # Process by size
//!     for position, score in blobs:
//!         x, y, z = position
//!         if score > 50:
//!             print(f"Large blob at ({x:.2f}, {y:.2f}, {z:.2f}) - "
//!                   f"score: {score:.2f} - likely ligand/metal")
//!         elif score > 10:
//!             print(f"Medium blob at ({x:.2f}, {y:.2f}, {z:.2f}) - "
//!                   f"score: {score:.2f} - likely waters")
//!         else:
//!             print(f"Small blob at ({x:.2f}, {y:.2f}, {z:.2f}) - "
//!                   f"score: {score:.2f} - check carefully")
//! 
//! # Sort blobs by score to prioritize investigation
//! sorted_blobs = sorted(blobs, key=lambda x: x[1], reverse=True)
//! 
//! # Investigate top 5 blobs
//! for position, score in sorted_blobs[:5]:
//!     x, y, z = position
//!     print(f"Priority blob: ({x:.2f}, {y:.2f}, {z:.2f}) score: {score:.2f}")
//!     # Navigate to this position in Coot to inspect visually
//!     coot.set_rotation_centre(x, y, z)
//! @endcode
//! 
//! @section advanced_usage Advanced Usage Patterns
//! 
//! Finding waters to place:
//! @code{.py}
//! blobs = coot.find_blobs_py(0, 2, 3.0)
//! for position, score in blobs:
//!     if 5 < score < 30:  # Typical water blob size
//!         x, y, z = position
//!         # Check hydrogen bonding potential before placing water
//!         # Could use coot.place_typed_atom_at_pointer("HOH")
//! @endcode
//! 
//! Finding large missing features:
//! @code{.py}
//! blobs = coot.find_blobs_py(0, 2, 3.0)
//! large_features = [(pos, score) for pos, score in blobs if score > 100]
//! for position, score in large_features:
//!     print(f"Large unmodeled density at {position}")
//!     print(f"  -> Check for missing residues, ligands, or nucleotides")
//! @endcode
//! 
//! Comparing difference map vs regular map:
//! @code{.py}
//! # Difference map - sensitive to model errors
//! diff_blobs = coot.find_blobs_py(0, 2, 3.0)  # mFo-DFc at 3σ
//! 
//! # Regular map - less model-biased
//! regular_blobs = coot.find_blobs_py(0, 1, 1.0)  # 2mFo-DFc at 1σ
//! 
//! # Features in both maps are high confidence
//! @endcode
//! 
//! @section validation Integration with Validation Workflows
//! 
//! Blob detection should be performed as part of comprehensive structure validation
//! alongside:
//! - Ramachandran analysis (local backbone geometry)
//! - Rotamer analysis (side-chain conformations)
//! - Atom overlap detection (packing problems)
//! - Density correlation analysis (model-to-map fit)
//! 
//! Complete validation example:
//! @code{.py}
//! # 1. Check model fit to map
//! rama = coot.all_molecule_ramachandran_score_py(0)
//! rotamers = coot.rotamer_graphs_py(0)
//! 
//! # 2. Check for packing problems
//! overlaps = coot.molecule_atom_overlaps_py(0, 30)
//! 
//! # 3. Find missing features (THIS FUNCTION)
//! blobs = coot.find_blobs_py(0, 2, 3.0)
//! 
//! # All three perspectives needed for complete validation
//! @endcode
//! 
//! @section notes Important Notes
//! 
//! - The returned blob positions are suitable for navigation with coot.set_rotation_centre()
//! - Blob scores are not directly comparable between different maps or molecules
//! - Very large blobs may indicate systematic problems (e.g., missing domains)
//! - Blob detection is computationally inexpensive and should be run routinely
//! - Results are affected by the find_ligand_mask_waters_flag global setting
//! - The masking radius is fixed at 1.9Å (reasonable for most cases)
//! 
//! @section see_also Related Functions
//! 
//! @see set_find_waters_sigma_cut_off() - Set default sigma for water finding
//! @see execute_find_waters_real() - Automated water placement
//! @see place_typed_atom_at_pointer() - Place atoms at blob positions
//! @see set_rotation_centre() - Navigate to blob positions
//! @see molecule_atom_overlaps_py() - Check for steric clashes
//! @see map_to_model_correlation_stats_per_residue_range_py() - Check model fit
//! 
//! @section references References
//! 
//! The blob detection algorithm is part of the ligand finding machinery in Coot.
//! For methodology details, see:
//! - Emsley & Cowtan (2004) Acta Cryst. D60, 2126-2132
//! - "Coot: model-building tools for molecular graphics"
//! 
//! @warning This function requires valid model and map molecules. Always check
//!          validity with is_valid_model_molecule() and is_valid_map_molecule()
//!          if molecule indices are not guaranteed.
//! 
//! @since This function has been available in Coot for many versions. The return
//!        format (list of [position, score] pairs) has been stable.
//! 
//! @note Reference counting: The returned PyObject is a new reference. The caller
//!       is responsible for DECREFing when finished. Py_False is properly INCREF'd
//!       before return on failure paths.
PyObject *find_blobs_py(int imol_model, int imol_map, float cut_off_sigma_density_level);
#endif

//! B-factor distribution histogram
void b_factor_distribution_graph(int imol);

/*  ----------------------------------------------------------------------- */
/*                  water chain                                             */
/*  ----------------------------------------------------------------------- */

//! \name Water Chain Functions
//! \{

#ifdef USE_GUILE
//! \brief return the chain id of the water chain from a shelx molecule.  Raw interface
//!
//!  @return scheme false if no chain or bad imol
SCM water_chain_from_shelx_ins_scm(int imol);
/*! \brief return the chain id of the water chain. Raw interface */
SCM water_chain_scm(int imol);
#endif

#ifdef USE_PYTHON
/* return the chain id of the water chain from a shelx molecule.  Raw interface.

Return False if no chain or bad imol*/
PyObject *water_chain_from_shelx_ins_py(int imol);
/*! \brief return the chain id of the water chain. Raw interface */
PyObject *water_chain_py(int imol);
#endif

//! \}


/*  ----------------------------------------------------------------------- */
/*                  interface utils                                          */
/*  ----------------------------------------------------------------------- */
//! \name Interface Utils
//! \{

/*! \brief Put text s into the status bar.

  use this to put info for the user in the statusbar (less intrusive
  than popup). */
void add_status_bar_text(const std::string &s);

//! set the logging level
//!
//! @param level is either "LOW" or "HIGH" or "DEBUGGING"
void set_logging_level(const std::string &level);//!


//! \}


/*  ----------------------------------------------------------------------- */
/*                  glyco tools                                             */
/*  ----------------------------------------------------------------------- */
//! \name Glyco Tools
//! \{

//! \brief print the glycosylation tree that contains the specified residue
void
print_glyco_tree(int imol, const std::string &chain_id, int resno, const std::string &ins_code);

//! \}


/*  ----------------------------------------------------------------------- */
/*                  variance map                                            */
/*  ----------------------------------------------------------------------- */
//! \name Variance Map
//! \{
//! \brief Make a variance map, based on the grid of the first map.
//!
//!  @return the molecule number of the new map.  Return -1 if unable to
//!   make a variance map.
int make_variance_map(const std::vector<int> &map_molecule_number_vec);
#ifdef USE_GUILE
int make_variance_map_scm(SCM map_molecule_number_list);
#endif
#ifdef USE_PYTHON
int make_variance_map_py(PyObject *map_molecule_number_list);
#endif
//! \}

/*  ----------------------------------------------------------------------- */
/*                  spin search                                             */
/*  ----------------------------------------------------------------------- */
//! \name Spin Search Functions
//! \{

void spin_search_by_atom_vectors(int imol_map, int imol, const std::string &chain_id, int resno, const std::string &ins_code, const std::pair<std::string, std::string> &direction_atoms_list, const std::vector<std::string> &moving_atoms_list);
#ifdef USE_GUILE
//! \brief for the given residue, spin the atoms in moving_atom_list
//!   around the bond defined by direction_atoms_list looking for the best
//!   fit to density of imol_map map of the first atom in
//!   moving_atom_list.  Works (only) with atoms in altconf ""
void spin_search(int imol_map, int imol, const char *chain_id, int resno, const char *ins_code, SCM direction_atoms_list, SCM moving_atoms_list);
//! \brief Spin N and CB (and the rest of the side chain if extant)
//!
//!  Sometime on N-terminal addition, then N ends up pointing the wrong way.
//!  The allows us to (more or less) interchange the positions of the CB and the N.
//!  angle is in degrees.
//!
void spin_N_scm(int imol, SCM residue_spec_scm, float angle);

//! \brief Spin search the density based on possible positions of CG of a side-chain.
//!
//! c.f. EM-Ringer
SCM CG_spin_search_scm(int imol_model, int imol_map);
#endif

#ifdef USE_PYTHON
//! for the given residue, spin the atoms in moving_atom_list...
//!
//!   around the bond defined by direction_atoms_list looking for the best
//!   fit to density of imom_map map of the first atom in
//!   moving_atom_list.  Works (only) with atoms in altconf ""
void spin_search_py(int imol_map, int imol, const char *chain_id, int resno, const char *ins_code, PyObject *direction_atoms_list, PyObject *moving_atoms_list);

//! Spin N and CB (and the rest of the side chain if extant)
//!
//!  Sometime on N-terminal addition, then N ends up pointing the wrong way.
//!  The allows us to (more or less) interchange the positions of the CB and the N.
//!  angle is in degrees.
//!
//! @param imol is the index of the model molecule
//! @param residue_spec is the specifier for the residue
//! @param angle is the rotation angle, in degrees, typically 120.
void spin_N_py(int imol, PyObject *residue_spec, float angle);

//! Spin search the density based on possible positions of CG of a side-chain
//!
//! @param imol_model is the index of the model molecule
//! @param imol_map is the index of the map molecule
//! @return either False (in the case of a failure) or a list of pairs of
//!         residue specifers and score - for each spinnable residue
//!         in the model.
PyObject *CG_spin_search_py(int imol_model, int imol_map);

#endif

//! \}


/*  ----------------------------------------------------------------------- */
/*                  monomer lib                                             */
/*  ----------------------------------------------------------------------- */
std::vector<std::pair<std::string, std::string> > monomer_lib_3_letter_codes_matching(const std::string &search_string, short int allow_minimal_descriptions_flag);


/*  ----------------------------------------------------------------------- */
/*                  mutate                                                  */
/*  ----------------------------------------------------------------------- */

int mutate_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end, const std::string &target_sequence);

int mutate_internal(int ires, const char *chain_id,
                    int imol, const std::string &target_res_type);
/* a function for multimutate to make a backup and set
   have_unsaved_changes_flag themselves */

//! \brief mutate active residue to single letter code slc
void mutate_active_residue_to_single_letter_code(const std::string &slc);

//! \brief show keyboard mutate frame
void show_keyboard_mutate_frame();

//! mutate by overlap
int mutate_by_overlap(int imol, const std::string &chain_id, int res_no, const std::string &new_type);


/*  ----------------------------------------------------------------------- */
/*                  ligands                                                 */
/*  ----------------------------------------------------------------------- */
coot::graph_match_info_t
overlap_ligands_internal(int imol_ligand, int imol_ref, const char *chain_id_ref,
                         int resno_ref, bool apply_rtop_flag);

//! \brief display the SMILES entry. This is the simple version - no dictionary
//! is generated.
void do_smiles_to_simple_3d_overlay_frame();

//! \brief get residues in the specified chain
//!
//! @param imol the molecule index
//! @param chain_id the specified chain-id
//!
//! @return a python list of residue specs for the residues in the given chain
PyObject *get_residues_in_chain_py(int imol, const std::string &chain_id);

//! does the specfied residue exist?
//!
//! @param imol the molecule index
//! @param spec is the residue spec to test for existance
//! @return 0 for no, 1 for yes, -1 for error
int residue_exists_py(int imol, PyObject *residue_spec_py);

/*  ----------------------------------------------------------------------- */
/*                  conformers (part of ligand search)                      */
/*  ----------------------------------------------------------------------- */

//! \brief make conformers of the ligand search molecules, each in its
//!  own molecule.

//! \name Extra Ligand Functions
//! \{

#ifdef USE_GUILE

//! make conformations
//!
//! as if for a ligand search
//!
//! Don't search the density.
//!
//! @return a list of new molecule numbers
SCM ligand_search_make_conformers_scm();
#endif

#ifdef USE_PYTHON
PyObject *ligand_search_make_conformers_py();

//! \brief get an rdkit molecule as a pickled string
//!
//! @param imol the index of the molecule
//! @param residue spec the residue specifier, e..g ['A', 11, ""]
//! @return pickled string. Return empty string on failure.
std::string get_rdkit_mol_base64_from_molecule(int imol, PyObject *residue_spec);

//! \brief and back the other way - import an RDKit mol in base64-encoded binary format
//!
//! @return the index of the new molecule - or -1 on failure
int molecule_from_rdkit_mol_base64(const std::string &rdkit_mol, PyObject *atom_name_list, const std::string &comp_id);

// make minimal restraints from mol (bonds and atoms)
int restraints_from_rdkit_mol_base64(const std::string &rdkit_mol_binary_base64, PyObject *atom_name_list_py, const std::string &comp_id);

#endif

std::vector<int> ligand_search_make_conformers_internal();

//! \}

/*  ----------------------------------------------------------------------- */
//                  animated ligand interactions
/*  ----------------------------------------------------------------------- */
void add_animated_ligand_interaction(int imol, const pli::fle_ligand_bond_t &lb);


/*  ----------------------------------------------------------------------- */
/*                  Cootaneer                                               */
/*  ----------------------------------------------------------------------- */
int cootaneer_internal(int imol_map, int imol_model, const coot::atom_spec_t &atom_spec);

//! \name Dock Sidechains
//! \{

#ifdef USE_GUILE
//! \brief cootaneer (i.e. assign sidechains onto mainchain model)
//!
//! atom_in_fragment_atom_spec is any atom spec in the fragment that should be
//! assigned with sidechains.
//!
//! @return the success status (0 is fail).
int cootaneer(int imol_map, int imol_model, SCM atom_in_fragment_atom_spec);
#endif

#ifdef USE_PYTHON
//! \brief cootaneer (i.e. assign sidechains onto mainchain model)
//!
//! atom_in_fragment_atom_spec is any atom spec in the fragment that should be
//! assigned with sidechains.
//!
//! @return the success status (0 is fail).
int cootaneer_py(int imol_map, int imol_model, PyObject *atom_in_fragment_atom_spec);
#endif

//! \}

/*  ----------------------------------------------------------------------- */
/*                  Sequence from Map                                       */
/*  ----------------------------------------------------------------------- */

//! \name Sequence from Map
//! \{
//! \brief 
//!
//! Use the map to estimate the sequence - you will need a decent map
//! 
//! @return the guessed sequence (empty is fail).
std::string sequence_from_map(int imol, const std::string &chain_id,
                              int resno_start, int resno_end, int imol_map);

void apply_sequence_to_fragment(int imol, const std::string &chain_id, int resno_start, int resno_end,
                                int imol_map, const std::string &file_name_for_sequences);

void assign_sequence_to_active_fragment();

//! \}

/*  ----------------------------------------------------------------------- */
/*                  Generic Objects                                         */
/*  ----------------------------------------------------------------------- */

// return a clean name and a flag to say that this was something that
// we were interested to make a graphics object from (rather than just
// header info)
std::pair<short int, std::string> is_interesting_dots_object_next_p(const std::vector<std::string> &vs);

/*  ----------------------------------------------------------------------- */
/*                  Generic Functions                                       */
/*  ----------------------------------------------------------------------- */
#ifdef USE_GUILE
SCM generic_string_vector_to_list_internal(const std::vector<std::string> &v);
SCM generic_int_vector_to_list_internal(const std::vector<int> &v);
std::vector<std::string> generic_list_to_string_vector_internal(SCM l);
SCM rtop_to_scm(const clipper::RTop_orth &rtop);
SCM inverse_rtop_scm(SCM rtop_scm);
// expects an expr of length 5, ie: (list chain-id res-no ins-cod atom-name alt-conf)
coot::atom_spec_t atom_spec_from_scm_expression(SCM expr);
SCM atom_spec_to_scm(const coot::atom_spec_t &spec);
#endif        /* USE_GUILE */

#ifdef USE_PYTHON
PyObject *generic_string_vector_to_list_internal_py(const std::vector<std::string>&v);
PyObject *generic_int_vector_to_list_internal_py(const std::vector<int> &v);
std::vector<std::string> generic_list_to_string_vector_internal_py(PyObject *l);
PyObject *rtop_to_python(const clipper::RTop_orth &rtop);
PyObject *inverse_rtop_py(PyObject *rtop_py);
coot::atom_spec_t atom_spec_from_python_expression(PyObject *expr);
PyObject *atom_spec_to_py(const coot::atom_spec_t &spec);
#endif // PYTHON

void set_display_control_button_state(int imol, const std::string &button_type, int state);

void fullscreen();
void unfullscreen();

//! Use left-mouse for view rotation
void set_use_trackpad(short int state);

//! this is an alias for the above
void set_use_primary_mouse_button_for_view_rotation(short int state);



/*  ----------------------------------------------------------------------- */
/*                  LIBCURL/Download                                        */
/*  ----------------------------------------------------------------------- */

/*! \brief if possible, read in the new coords getting coords via web.

(no return value because get-url-str does not return one).
 */
void get_coords_for_accession_code(const std::string &code);

// internal use (strings, not binaries).
std::string coot_get_url_as_string_internal(const char *url);

void stop_curl_download(const char *file_name); // stop curling the to file_name;

std::string get_drug_mdl_via_wikipedia_and_drugbank(std::string drugname);

//! \brief fetch and superpose AlphaFold models corresponding to model
//!
//! model must have Uniprot DBREF info in the header.
void fetch_and_superpose_alphafold_models(int imol);

//! \brief return the model number
int fetch_alphafold_model_for_uniprot_id(const std::string &uniprot_id);

//! \brief Loads up map frmo emdb
void fetch_emdb_map(const std::string &emd_accession_code);

//! \brief return the COD entry, return a molecule index
int fetch_cod_entry(const std::string &cod_entry_id);


/*  ----------------------------------------------------------------------- */
/*                  Functions for FLEV layout callbacks                     */
/*  ----------------------------------------------------------------------- */
// orient the graphics somehow so that the interaction between
// central_residue and neighbour_residue is perpendicular to screen z.
void orient_view(int imol,
                 const coot::residue_spec_t &central_residue_spec, // ligand typically
                 const coot::residue_spec_t &neighbour_residue_spec);

/*  \brief return a list of chiral centre ids as determined from topological
    equivalence analysis based on the bond info (and element names). */
std::vector<std::string>
topological_equivalence_chiral_centres(const std::string &residue_type);


/*  ----------------------------------------------------------------------- */
/*                  New Screendump                                          */
/*  ----------------------------------------------------------------------- */
void screendump_tga(const std::string &file_name);
void set_framebuffer_scale_factor(unsigned int sf);

/*  ----------------------------------------------------------------------- */
/*                  New Graphics Control                                    */
/*  ----------------------------------------------------------------------- */

//! \brief set use perspective mode
void set_use_perspective_projection(short int state);

//! \brief query if perspective mode is being used
int use_perspective_projection_state();

//! \brief set the perspective fov. Default 20 degrees.
void set_perspective_fov(float degrees);

//! \brief set use ambient occlusion
void set_use_ambient_occlusion(short int state);
//! \brief query use ambient occlusion
int use_ambient_occlusion_state();

//! \brief set use depth blur
void set_use_depth_blur(short int state);
//! \brief query use depth blur
int use_depth_blur_state();

//! \brief set use fog
void set_use_fog(short int state);

//! \brief query use fog
int use_fog_state();

//! \brief set use ourline
void set_use_outline(short int state);

//! \brief query use outline
int use_outline_state();

//! \brief set the map shininess
void set_map_shininess(int imol, float shininess);

//! \brief set the map specular strength
void set_map_specular_strength(int imol, float specular_strength);

//! \brief
void set_draw_normals(short int state);

//! \brief
int draw_normals_state();

//! \brief
void set_draw_mesh(int imol, int mesh_index, short int state);

//! \brief return -1 on unable to lookup mesh
int draw_mesh_state(int imol, int mesh_index);

//! \brief
void set_map_material_specular(int imol, float specular_strength, float shininess);

//! \brief
void set_model_material_specular(int imol, float specular_strength, float shininess);

//! \brief set the ambient
void set_model_material_ambient(int imol, float r, float g, float b, float alpha);

//! \brief
void set_model_material_diffuse(int imol, float r, float g, float b, float alpha);

//! \brief set the goodselliness (pastelization_factor) 0.3 is about right, but "the right value"
//!        depends on the renderer, may be some personal choice.
void set_model_goodselliness(float pastelization_factor);

//! \brief
void set_map_fresnel_settings(int imol, short int state, float bias, float scale, float power);

//! \brief
void reload_map_shader();

//! \brief
void reload_model_shader();

//! \brief
void set_atom_radius_scale_factor(int imol, float scale_factor);

//! \brief set use fancy rendering lighting
//!
//! Turn on framebuffer effects
//!
//! @param state where 1 mean turn on and 0 means turn off.
void set_use_fancy_lighting(short int state);

//! \brief set use simple lines for model molecule
//!
//! @param state where 1 mean turn on and 0 means turn off.
void set_use_simple_lines_for_model_molecules(short int state);

//! \brief
//!
//! @param state where 1 mean turn on and 0 means turn off.
void set_fresnel_colour(int imol, float red, float green, float blue, float opacity);

//! \brief
void set_focus_blur_z_depth(float z);

//! \brief set use depth blur
void set_use_depth_blur(short int state);

//! \brief set focus blur strength
void set_focus_blur_strength(float st);

//! \brief set shadow stren
//!
//! @param s is the shadow strength between 0 and 1.
void set_shadow_strength(float s);

//! \brief set the shadow resolution (1,2,3,4)
void set_shadow_resolution(int reso_multiplier);

//! \brief set shadow box size - default 66;
void set_shadow_box_size(float size);

//! \brief set SSAO kernel n samples
void set_ssao_kernel_n_samples(unsigned int n_samples);

//! \brief set SSAO strength
//!
//! screen-space ambient occlusionn
//!
//! @param strength is the SSAO strength between 0 and 1.
void set_ssao_strength(float strength);

//! \brief set SSAO radius
//!
//! screen-space ambient occlusionn
//! Doesn't do much. Not worth adjusting
//!
//! @param radius is the SSAO radius.
void set_ssao_radius(float radius);

//! \brief set SSAO bias
//!
//! screen-space ambient occlusionn
//! Doesn't do much. Not worth adjusting
//!
void set_ssao_bias(float bias);

//! \brief set SSAO blur size (0, 1, or 2)
void set_ssao_blur_size(unsigned int blur_size);

//! \brief set the shadow softness (1, 2 or 3)
void set_shadow_softness(unsigned int softness);

//! \brief set the shadow softness (1, 2 or 3)
void set_shadow_texture_resolution_multiplier(unsigned int m);

//! \brief adjust the effects shader output type (for debugging effects)
void set_effects_shader_output_type(unsigned int type);

//! \brief adjust the effects shader brightness
void set_effects_shader_brightness(float f);

//! \brief adjust the effects shader gamma
void set_effects_shader_gamma(float f);

//! \brief set bond smoothness (default 1 (not smooth))
void set_bond_smoothness_factor(unsigned int fac);

//! \brief set the draw state of the Ramachandran plot display during Real Space Refinement
void set_draw_gl_ramachandran_plot_during_refinement(short int state);

//! \brief set the FPS timing scale factor - default 0.0025
void set_fps_timing_scale_factor(float f);

//! \brief draw background image
void set_draw_background_image(bool state);

// testing function
void read_test_gltf_models();

//! \brief load a gltf model
//!
//! If the gltf
//! files does not exist, an empty model will be created
//!
//! @param gltf_file_name is the name of the gltf file to load
//! @return the model index of the loaded model.
int load_gltf_model(const std::string &gltf_file_name);

//! \brief set the model animation parameters
void set_model_animation_parameters(unsigned int model_index, float amplitude, float wave_numer, float freq);

//! \brief enable/disable the model animation (on or off)
void set_model_animation_state(unsigned int model_index, bool state);

//! \brief load a gltf model
void scale_model(unsigned int model_index, float scale_factor);

//! \brief reset the frame buffers
void reset_framebuffers();


/*  ----------------------------------------------------------------------- */
/*               Return Rotamer score (don't touch the model)               */
/*  ----------------------------------------------------------------------- */

//! \name Rotamer Scoring
//! \{

std::vector<coot::named_rotamer_score> score_rotamers(int imol,
                                                      const char *chain_id,
                                                      int res_no,
                                                      const char *ins_code,
                                                      const char *alt_conf,
                                                      int imol_map,
                                                      int clash_flag,
                                                      float lowest_probability);

#ifdef USE_GUILE
//! \brief Score rotamers for a residue (Guile interface)
//!
//! Returns a list of possible rotamer conformations with their scores.
//! Each rotamer is scored based on rotamer library probability and
//! (optionally) density fit.
//!
//! @param imol Model molecule number
//! @param chain_id Chain identifier
//! @param res_no Residue number
//! @param ins_code Insertion code
//! @param alt_conf Alternate conformation
//! @param imol_map Map for density scoring (-1 to skip)
//! @param clash_flag 1 to check clashes, 0 to skip
//! @param lowest_probability Minimum probability threshold (0.0-1.0)
//!
//! @return SCM - List of rotamer descriptions, each containing rotamer name,
//!         probability, and density score. Empty list if residue not found
//!         or no rotamers above threshold.
//!
//! \note The density score is only meaningful if imol_map is a valid map
//!
//! Example usage:
//! \code{.scm}
//! ;; Score rotamers for LEU 42 in chain A
//! (define rotamers (score-rotamers-scm 1 "A" 42 "" "" 2 1 0.01))
//! (for-each
//!   (lambda (rot)
//!     (format #t "Rotamer: ~a, Probability: ~a, Fit: ~a~%"
//!             (list-ref rot 0)  ; name
//!             (list-ref rot 1)  ; probability
//!             (list-ref rot 2))) ; density fit
//!   rotamers)
//! \endcode
SCM score_rotamers_scm(int imol,
                       const char *chain_id,
                       int res_no,
                       const char *ins_code,
                       const char *alt_conf,
                       int imol_map,
                       int clash_flag,
                       float lowest_probability);
#endif

#ifdef USE_PYTHON
//! \brief Score all rotamers for a residue (Python interface)
//!
//! **USEFUL FOR FIXING BAD ROTAMERS**
//!
//! Evaluates all possible rotamer conformations and returns them with scores.
//! This is the function to call before using auto_fit_best_rotamer.
//!
//! @param imol Model molecule number
//! @param chain_id Chain identifier
//! @param res_no Residue number
//! @param ins_code Insertion code (use "" if none)
//! @param alt_conf Alternate conformation (use "" for default)
//! @param imol_map Map molecule for density scoring (use -1 to ignore density)
//! @param clash_flag 1 to check for clashes with other atoms, 0 to skip
//! @param lowest_probability Filter: only return rotamers above this probability
//!
//! @return PyObject* - List of rotamer dictionaries, each containing:
//!         - 'name': Rotamer name (e.g., "mt-85")
//!         - 'probability': Rotamer library probability (0.0-1.0)
//!         - 'density_score': Fit to density (if map provided)
//!         - 'richardson_name': Rotamer name in Richardson notation
//!         Empty list if no suitable rotamers found.
//!
//! \note Rotamers are ranked by combined probability and density fit
//! \note Use clash_flag=1 to avoid rotamers that clash with nearby atoms
//!
//! Example usage:
//! \code{.py}
//! # Score rotamers for LEU 42, considering density and clashes
//! rotamers = score_rotamers_py(
//!     imol=1,
//!     chain_id="A",
//!     res_no=42,
//!     ins_code="",
//!     alt_conf="",
//!     imol_map=2,           # Use map 2 for density scoring
//!     clash_flag=1,         # Check for clashes
//!     lowest_probability=0.01  # Only show rotamers >1% probability
//! )
//!
//! print(f"Found {len(rotamers)} possible rotamers")
//! for rot in rotamers:
//!     print(f"{rot['name']}: "
//!           f"prob={rot['probability']:.1%}, "
//!           f"fit={rot['density_score']:.3f}")
//!
//! # The best rotamer is typically first in the list
//! if rotamers:
//!     best = rotamers[0]
//!     print(f"Best rotamer: {best['name']}")
//! \endcode
PyObject *score_rotamers_py(int imol,
                            const char *chain_id,
                            int res_no,
                            const char *ins_code,
                            const char *alt_conf,
                            int imol_map,
                            int clash_flag,
                            float lowest_probability);
#endif

//! \}

/*  ----------------------------------------------------------------------- */
/*               Use Cowtan's protein_db to discover loops                  */
/*  ----------------------------------------------------------------------- */
/*! \name protein-db */
/* \{ */
/*! \brief Cowtan's protein_db loops */
// return in the first pair, the imol of the new molecule generated
// from an atom selection of the imol_coords for the residue selection
// of the loop and the molecule number of the consolidated solutions
// (displayed in purple).  and the second of the outer pair, there is
// vector of molecule indices for each of the candidate loops.
//
// return -1 in the first of the pair on failure
//
std::pair<std::pair<int, int> , std::vector<int> >
protein_db_loops(int imol_coords,
                 const std::vector<coot::residue_spec_t> &residue_specs,
                 int imol_map, int nfrags, bool preserve_residue_names);
// so that we can create a "original loop" molecule from the atom
// specs picked (i.e. the atom selection string should extend over the
// range from the smallest residue number to the largest (in the same
// chain)).
std::string
protein_db_loop_specs_to_atom_selection_string(const std::vector<coot::residue_spec_t> &specs);
#ifdef USE_GUILE
SCM protein_db_loops_scm(int imol_coords, SCM residues_specs, int imol_map, int nfrags, bool preserve_residue_names);
#endif
#ifdef USE_PYTHON
PyObject *protein_db_loops_py(int imol_coords, PyObject *residues_specs, int imol_map, int nfrags, bool preserve_residue_names);
#endif
/* \} */


/* ------------------------------------------------------------------------- */
/*                      HOLE                                                 */
/* ------------------------------------------------------------------------- */
/*! \name Coot's Hole implementation */

/*! \brief starting point and end point, colour map multiplier and
  shall the probe radius graph be shown (dummy value currently).

 if export_dots_file_name string length is zero, then don't try to
 export the surface dots.

*/
void hole(int imol,
          float start_x, float start_y, float start_z,
          float   end_x, float   end_y, float   end_z,
          float colour_map_multiplier, float colour_map_offset,
          int n_runs, bool show_probe_radius_graph_flag,
          std::string export_surface_dots_file_name);


/* ------------------------------------------------------------------------- */
/*                      Gaussian Surface                                     */
/* ------------------------------------------------------------------------- */
/*! \name Coot's Gaussian Surface */

//! \brief The surface is separated into chains to make Generid Display Objects
//!        There is no colour, contour, grid or sigma or material control yet for Gaussian surfaces.
//!
//!        there should be one day...
int gaussian_surface(int imol);

//! \brief set the sigma for gaussian surface (default 4.0)
void set_gaussian_surface_sigma(float s);

//! \brief set the contour_level for gaussian surface (default 4.4)
void set_gaussian_surface_contour_level(float s);

//! \brief set the box_radius for gaussian surface (defautl 5)
void set_gaussian_surface_box_radius(float s);

//! \brief set the grid_scale for gaussian surface (default 0.7)
void set_gaussian_surface_grid_scale(float s);
//! \brief set the fft B-factor for gaussian surface. Use 0 for no B-factor (default 100)
void set_gaussian_surface_fft_b_factor(float f);

//! \brief set the chain colour mode for Gaussian surfaces
//! mode = 1 means each chain has its own colour
//! mode = 2 means the chain colour is determined from NCS/molecular symmetry (so
//!         that, in this mode, chains with the same sequence have the same colour
void set_gaussian_surface_chain_colour_mode(short int mode);

void show_gaussian_surface_overlay();


/* ------------------------------------------------------------------------- */
/*                      Acedrg for dictionary                                */
/* ------------------------------------------------------------------------- */
void make_acedrg_dictionary_via_CCD_dictionary(int imol, const coot::residue_spec_t &spec);

/* ------------------------------------------------------------------------- */
/*                      LINKs                                                */
/* ------------------------------------------------------------------------- */

//! \brief make a link between the specified atoms
void
make_link(int imol, const coot::atom_spec_t &spec_1, const coot::atom_spec_t &spec_2,
          const std::string &link_name, float length);
#ifdef USE_GUILE
void make_link_scm(int imol, SCM spec_1, SCM spec_2, const std::string&link_name, float length);
// return a list of the links in the given molecule.
// will return an empty list for non-valid (i.e. non-model) molecules
//
SCM link_info_scm(int imol);
#endif
#ifdef USE_PYTHON
void make_link_py(int imol, PyObject *spec_1, PyObject *spec_2, const std::string&link_name, float length);
// return a list of the links in the given molecule.
// will return an empty list for non-valid (i.e. non-model) molecules
//
PyObject *link_info_py(int imol);
#endif

void show_acedrg_link_interface_overlay();

/* ------------------------------------------------------------------------- */
/*                      Drag and drop                                        */
/* ------------------------------------------------------------------------- */

/*! \name  Drag and Drop Functions */
// \{
//! \brief handle the string that get when a file or URL is dropped.
int handle_drag_and_drop_string(const std::string &uri);
// \}


/* ------------------------------------------------------------------------- */
/*                      Map Display Control                                  */
/* ------------------------------------------------------------------------- */

/*! \name  Map Display Control */
// \{
//! \brief undisplay all maps except the given one
void undisplay_all_maps_except(int imol_map);
// \}


/* ------------------------------------------------------------------------- */
/*                      Map Contours                                         */
/* ------------------------------------------------------------------------- */

/*! \name Map Contouring Functions */

#ifdef USE_PYTHON
// \{
//! \brief return a list of pairs of vertices for the lines
PyObject *map_contours(int imol, float contour_level);

//! \brief return two lists: a list of vertices and a list of index-triples for connection
PyObject *map_contours_as_triangles(int imol, float contour_level);

// \}
#endif // USE_PYTHON

//! \brief enable radial map colouring
void set_radial_map_colouring_enabled(int imol, int state);

//! \brief radial map colouring centre
void set_radial_map_colouring_centre(int imol, float x, float y, float z);

//! \brief radial map colouring min
void set_radial_map_colouring_min_radius(int imol, float r);

//! \brief radial map colouring max
void set_radial_map_colouring_max_radius(int imol, float r);

//! \brief radial map colouring inverted colour map
void set_radial_map_colouring_invert(int imol, int invert_state);

//! \brief radial map colouring saturation
//!
//! saturation is a number between 0 and 1, typically 0.5
void set_radial_map_colouring_saturation(int imol, float saturation);



/* ------------------------------------------------------------------------- */
/*                      correlation maps                                     */
/* ------------------------------------------------------------------------- */

//! \name Map to Model Correlation
//! \{


//! \brief The atom radius is not passed as a parameter to correlation
//         functions, let's set it here (default is 1.5A)
//
void set_map_correlation_atom_radius(float r);

// Don't count the grid points of residues_specs that are in grid
// points of (potentially overlapping) neighbour_residue_spec.
//
#ifdef USE_GUILE
//! \brief atom-mask-mode is as follows:
// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-excluding CB if is standard amino-acid, else all atoms
// 4: main-chain atoms if is standard amino-acid, else nothing
// 5: side-chain atoms if is standard amino-acid, else nothing
// 10: atom radius is dependent atom atom B-factor
SCM map_to_model_correlation_scm(int imol,
                                 SCM residue_specs,
                                 SCM neighb_residue_specs,
                                 unsigned short int atom_mask_mode,
                                 int imol_map);

//! \brief Map-to-model correlation statistics (Guile interface)
//!
//! @param imol Model molecule number
//! @param residue_specs Scheme list of residue specs
//! @param neighb_residue_specs Scheme list of neighboring residue specs
//! @param atom_mask_mode Controls which atoms to include
//! @param imol_map Map molecule number
//!
//! @return SCM - Association list with statistics
SCM map_to_model_correlation_stats_scm(int imol,
                                       SCM residue_specs,
                                       SCM neighb_residue_specs,
                                       unsigned short int atom_mask_mode,
                                       int imol_map);
#endif

#ifdef USE_PYTHON
//! \brief Calculate map-to-model correlation (Python interface)
//!
//! Python wrapper for map_to_model_correlation. Evaluates the fit of specific
//! residues to the electron density map.
//!
//! @param imol Model molecule number
//! @param residue_specs Python list of residue specs [[chain_id, resno, ins_code], ...]
//! @param neighb_residue_specs Python list of neighboring residue specs to exclude
//! @param atom_mask_mode Controls which atoms to include (see atom_mask_mode values above)
//! @param imol_map Map molecule number
//!
//! @return PyObject* - correlation coefficient as a Python float
//!
//! \note Use atom_mask_mode=2 to evaluate side-chain fit specifically
//!
//! Example usage:
//! \code{.py}
//! # Evaluate side-chain fit for residues 40-44
//! residue_specs = [['A', res_no, ''] for res_no in range(40, 45)]
//! correlation = map_to_model_correlation_py(
//!     imol=1,
//!     residue_specs=residue_specs,
//!     neighb_residue_specs=[],
//!     atom_mask_mode=2,  # Side-chain atoms only
//!     imol_map=2
//! )
//! print(f"Side-chain correlation: {correlation}")
//! \endcode
PyObject *map_to_model_correlation_py(int imol,
                                      PyObject *residue_specs,
                                      PyObject *neighb_residue_specs,
                                      unsigned short int atom_mask_mode,
                                      int imol_map);

//! \brief Get map-to-model correlation statistics (Python interface)
//!
//! Returns detailed statistics about the correlation between model and map,
//! including mean, standard deviation, and range.
//!
//! @param imol Model molecule number
//! @param residue_specs Python list of residue specs
//! @param neighb_residue_specs Python list of neighboring residue specs
//! @param atom_mask_mode Controls which atoms to include
//! @param imol_map Map molecule number
//!
//! @return PyObject* - Dictionary with statistics:
//!         - 'mean': Mean correlation
//!         - 'std_dev': Standard deviation
//!         - 'min': Minimum correlation
//!         - 'max': Maximum correlation
//!
//! Example usage:
//! \code{.py}
//! stats = map_to_model_correlation_stats_py(1, residues, [], 0, 2)
//! print(f"Mean: {stats['mean']:.3f}")
//! print(f"Std Dev: {stats['std_dev']:.3f}")
//! print(f"Range: {stats['min']:.3f} to {stats['max']:.3f}")
//! \endcode
PyObject *map_to_model_correlation_stats_py(int imol,
                                      PyObject *residue_specs,
                                      PyObject *neighb_residue_specs,
                                      unsigned short int atom_mask_mode,
                                      int imol_map);

//! \brief Get density statistics per residue range (Python interface)
//!
//! **PRIMARY FUNCTION FOR FINDING POORLY-FITTED RESIDUES**
//!
//! This is the main function to use when asked "Which side chain is worst fitting to density?"
//! It analyzes correlation statistics for all residues in a chain and returns comprehensive
//! data for both all-atom and side-chain-only analysis.
//!
//! @param imol Model molecule number
//! @param chain_id Chain identifier (e.g., "A", "B")
//! @param imol_map Map molecule number
//! @param n_residue_per_residue_range Number of residues per analysis window:
//!        - Use 1 for per-residue statistics (most common)
//!        - Use 3 for smoothed statistics over 3-residue windows
//! @param exclude_backbone_NOC_flag Whether to exclude backbone N, O, C atoms:
//!        - 0: Include all atoms (for overall fit assessment)
//!        - 1: Exclude N, O, C (for side-chain-focused analysis)
//!
//! @return PyObject* - List of two lists (first is all_atom_stats, second is sidechain_stats):
//!         Each list contains a residue_spec [chain_id, res_no, ins_code]  and correlation_stats containing:
//!         - n_points: Number of grid points
//!         - correlation: Correlation coefficient
//!
//! If the residue does not have a side-chain then the number of grid points is 0 and the
//! correlation is nan.
//!
//! \note This function analyzes the ENTIRE chain at once, making it very efficient
//! \note Returns both all-atom and side-chain statistics in one call
//! \note Use the side-chain statistics to identify problem side chains specifically
//!
//! Example usage - Find worst-fitting side chain:
//! \code{.py}
//! # Get correlation statistics for all residues in chain A
//! all_atom_stats, sidechain_stats = map_to_model_correlation_stats_per_residue_range_py(
//!     imol=1,              # Model molecule
//!     chain_id="A",        # Chain A
//!     imol_map=2,          # Map molecule
//!     n_residue_per_residue_range=1,  # Per-residue (not averaged)
//!     exclude_NOC_flag=0   # Include all atoms
//! )
//!
//! # Find worst-fitting side chain
//! worst_residue = min(sidechain_stats.items(), key=lambda x: x[1][1])
//!
//! chain_id, resno, ins_code = worst_residue[0]
//! correlation = worst_residue[1][1]
//!
//! print(f"Worst side chain: {chain_id} {resno}, correlation = {correlation:.3f}")
//!
//! # Center on worst residue
//! set_go_to_atom_chain_residue_atom_name(chain_id, resno, 'CA')
//! \endcode
//!
//! Example usage - Compare main-chain vs side-chain fit:
//! \code{.py}
//! all_atom, sidechain = map_to_model_correlation_stats_per_residue_range_py(
//!     1, "A", 2, 1, 0
//! )
//!
//! for residue_spec in all_atom.keys():
//!     all_corr = all_atom[residue_spec]['correlation']
//!     side_corr = sidechain[residue_spec]['correlation']
//!
//!     if all_corr > 0.7 and side_corr < 0.5:
//!         print(f"Residue {residue_spec}: Good backbone, poor sidechain")
//! \endcode
PyObject *
map_to_model_correlation_stats_per_residue_range_py(int imol,
                                                    const std::string &chain_id,
                                                    int imol_map,
                                                    unsigned int n_residue_per_residue_range,
                                                    short int exclude_mainchain_NOC_flag);

#endif

// Map to Model Correlation Functions - Enhanced Doxygen Documentation
// 
// These functions assess how well a molecular model fits into electron density maps.
// Essential for model validation and identifying poorly-fitted regions.

//! \name Map to Model Correlation
//! \{

//! \brief Calculate the correlation between a map and model for specific residues
//!
//! This function calculates the map-to-model correlation for specified residues,
//! excluding grid points that overlap with neighboring residues.
//!
//! @param imol Model molecule number
//! @param residue_specs Vector of residue specifications to evaluate
//! @param neigh_residue_specs Vector of neighboring residues whose grid points should be excluded
//! @param atom_mask_mode Controls which atoms are included in the calculation:
//!        - 0: All atoms
//!        - 1: Main-chain atoms if standard amino acid, else all atoms
//!        - 2: Side-chain atoms if standard amino acid, else all atoms
//!        - 3: Side-chain atoms excluding CB if standard amino acid, else all atoms
//!        - 4: Main-chain atoms if standard amino acid, else nothing
//!        - 5: Side-chain atoms if standard amino acid, else nothing
//!        - 10: Atom radius dependent on B-factor
//! @param imol_map Map molecule number to correlate against
//! 
//! @return Correlation coefficient (float) between model and map
//!
//! \note Use this after refinement to evaluate if the fit improved
//!
//! Example:
//! \code{.cpp}
//! std::vector<coot::residue_spec_t> residues = {{chain_id, 40, ""}, {chain_id, 41, ""}};
//! std::vector<coot::residue_spec_t> neighbors;
//! float corr = map_to_model_correlation(1, residues, neighbors, 2, 2);
//! \endcode
//! \brief atom-mask-mode is as follows:
float
map_to_model_correlation(int imol,
                         const std::vector<coot::residue_spec_t> &residue_specs,
                         const std::vector<coot::residue_spec_t> &neigh_residue_specs,
                         unsigned short int atom_mask_mode,
                         int imol_map);

//! \brief Get detailed statistics for map-to-model correlation
//!
//! Returns comprehensive statistics including mean, standard deviation, minimum,
//! and maximum correlation values for the specified residues.
//!
//! @param imol Model molecule number
//! @param residue_specs Vector of residue specifications to evaluate
//! @param neigh_residue_specs Vector of neighboring residues to exclude from grid
//! @param atom_mask_mode Controls which atoms are included (see atom_mask_mode values)
//! @param imol_map Map molecule number
//!
//! @return coot::util::density_correlation_stats_info_t containing:
//!         - mean correlation
//!         - standard deviation
//!         - minimum correlation
//!         - maximum correlation
//!         - number of residues evaluated
//!
//! \note Useful for getting overall quality metrics for a region
//!
//! Example:
//! \code{.cpp}
//! auto stats = map_to_model_correlation_stats(1, residues, neighbors, 0, 2);
//! std::cout << "Mean correlation: " << stats.mean << std::endl;
//! std::cout << "Std dev: " << stats.std_dev << std::endl;
//! \endcode
coot::util::density_correlation_stats_info_t
map_to_model_correlation_stats(int imol,
                               const std::vector<coot::residue_spec_t> &residue_specs,
                               const std::vector<coot::residue_spec_t> &neigh_residue_specs,
                               unsigned short int atom_mask_mode,
                               int imol_map);
#ifndef SWIG

//! \brief Get map-to-model correlation per residue
//!
//! Returns individual correlation values for each specified residue.
//! Useful for identifying which specific residues fit poorly.
//!
//! @param imol Model molecule number
//! @param specs Vector of residue specifications
//! @param atom_mask_mode Controls which atoms to include
//! @param imol_map Map molecule number
//!
//! @return Vector of pairs: (residue_spec, correlation_value)
//!
//! \note The correlation values can be sorted to find worst-fitting residues
//!
//! Example:
//! \code{.cpp}
//! auto correlations = map_to_model_correlation_per_residue(1, specs, 0, 2);
//! // Sort by correlation (lowest first)
//! std::sort(correlations.begin(), correlations.end(),
//!           [](const auto &a, const auto &b) { return a.second < b.second; });
//! // First element is now the worst-fitting residue
//! std::cout << "Worst residue: " << correlations[0].first 
//!           << " correlation: " << correlations[0].second << std::endl;
//! \endcode
std::vector<std::pair<coot::residue_spec_t,float> >
map_to_model_correlation_per_residue(int imol, const std::vector<coot::residue_spec_t> &specs,
                                     unsigned short int atom_mask_mode,
                                     int imol_map);

//! \brief Get detailed density statistics per residue
//!
//! Returns comprehensive density statistics for each residue including
//! correlation, mean, standard deviation, and other metrics.
//!
//! @param imol Model molecule number
//! @param residue_specs Vector of residue specifications
//! @param atom_mask_mode Controls which atoms to include
//! @param atom_radius_for_masking Radius around atoms for masking (typically 1.5 Å)
//! @param imol_map Map molecule number
//!
//! @return Map of residue_spec to density_stats_info_t containing detailed statistics
//!
//! \note Provides more detailed information than correlation alone
//!
//! Example:
//! \code{.cpp}
//! auto stats_map = map_to_model_correlation_stats_per_residue(
//!     1, specs, 0, 1.5, 2
//! );
//! for (const auto &pair : stats_map) {
//!     std::cout << pair.first << ": "
//!               << "corr=" << pair.second.correlation << ", "
//!               << "mean=" << pair.second.mean_density << std::endl;
//! }
//! \endcode
std::map<coot::residue_spec_t, coot::util::density_stats_info_t>
map_to_model_correlation_stats_per_residue(int imol,
                                           const std::vector<coot::residue_spec_t> &residue_specs,
                                           unsigned short int atom_mask_mode,
                                           float atom_radius_for_masking,
                                           int imol_map);

//! \brief Get density statistics per residue range
//!
//! Analyzes correlation statistics for residue ranges (windows) along a chain.
//! The middle residue of each range represents the statistics for that window.
//!
//! @param imol Model molecule number
//! @param chain_id Chain identifier (e.g., "A", "B")
//! @param imol_map Map molecule number
//! @param n_residue_per_residue_range Number of residues per analysis window (typically 1 for per-residue)
//! @param exclude_NOC_flag Whether to exclude backbone N, O, C atoms (1=yes, 0=no)
//!
//! @return Pair of maps:
//!         - first: All-atom statistics per residue
//!         - second: Side-chain-only statistics per residue
//!
//! \note This is the primary function for comprehensive chain-wide validation
//! \note Use n_residue_per_residue_range=1 for per-residue statistics
//! \note Use n_residue_per_residue_range=3 for smoothed statistics
//!
//! Example:
//! \code{.cpp}
//! auto [all_atom_stats, sidechain_stats] =
//!     map_to_model_correlation_stats_per_residue_range(1, "A", 2, 1, 0);
//!
//! // Find worst-fitting residue (all atoms)
//! auto worst = std::min_element(
//!     all_atom_stats.begin(), all_atom_stats.end(),
//!     [](const auto &a, const auto &b) {
//!         return a.second.correlation < b.second.correlation;
//!     }
//! );
//! \endcode
std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
          std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
map_to_model_correlation_stats_per_residue_range(int imol, const std::string &chain_id, int imol_map,
                                                 unsigned int n_residue_per_residue_range,
                                                 short int exclude_NOC_flag);

#endif // not for swigging.

#ifdef USE_GUILE
//! \brief Map-to-model correlation per residue (Guile interface)
//!
//! @param imol Model molecule number
//! @param residue_specs Scheme list of residue specs
//! @param atom_mask_mode Controls which atoms to include
//! @param imol_map Map molecule number
//!
//! @return SCM - List of pairs: ((residue-spec correlation) ...)
SCM
map_to_model_correlation_per_residue_scm(int imol, SCM residue_specs,
                                         unsigned short int atom_mask_mode,
                                         int imol_map);

//! \brief Map-to-model correlation stats per residue (Guile interface)
//!
//! @param imol Model molecule number
//! @param residue_specs_scm Scheme list of residue specs
//! @param atom_mask_mode Controls which atoms to include
//! @param atom_radius_for_masking Radius for atom masking
//! @param imol_map Map molecule number
//!
//! @return SCM - Association list mapping residue specs to stats
SCM
map_to_model_correlation_stats_per_residue_scm(int imol,
                                               SCM residue_specs_scm,
                                               unsigned short int atom_mask_mode,
                                               float atom_radius_for_masking,
                                               int imol_map);

SCM map_to_model_correlation_stats_per_residue_range_scm(int imol, const std::string &chain_id, int imol_map,
                                                         unsigned int n_residue_per_residue_range,
                                                         short int exclude_NOC_flag);


//! \brief QQ plot of the model density correlation, reported per residue
//!
//! \brief atom-mask-mode is as follows:
// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-excluding CB if is standard amino-acid, else all atoms
// 4: main-chain atoms if is standard amino-acid, else nothing
// 5: side-chain atoms if is standard amino-acid, else nothing
//
SCM qq_plot_map_and_model_scm(int imol,
                              SCM residue_specs_scm,
                              SCM neigh_residue_specs_scm,
                              unsigned short int atom_mask_mode,
                              int imol_map);
#endif

#ifdef USE_PYTHON
//! \brief Get map-to-model correlation per residue (Python interface)
//!
//! Returns correlation values individually for each specified residue.
//! Essential for identifying which specific residues have poor density fit.
//!
//! @param imol Model molecule number
//! @param residue_specs Python list of residue specs
//! @param atom_mask_mode Controls which atoms to include
//! @param imol_map Map molecule number
//!
//! @return PyObject* - List of tuples: [(residue_spec, correlation), ...]
//!
//! \note Sort the results to find worst-fitting residues
//!
//! Example usage:
//! \code{.py}
//! # Get all residues in chain A
//! residues = get_residues_in_chain_py(1, "A")
//!
//! # Get per-residue correlations
//! correlations = map_to_model_correlation_per_residue_py(
//!     imol=1,
//!     residue_specs=residues,
//!     atom_mask_mode=0,  # All atoms
//!     imol_map=2
//! )
//!
//! # Find worst 10 residues
//! worst_10 = sorted(correlations, key=lambda x: x[1])[:10]
//! for spec, corr in worst_10:
//!     print(f"Residue {spec}: correlation = {corr:.3f}")
//! \endcode
PyObject *map_to_model_correlation_per_residue_py(int imol, PyObject *residue_specs,
                                                  unsigned short int atom_mask_mode,
                                                  int imol_map);

PyObject *qq_plot_map_and_model_py(int imol,
                              PyObject *residue_specs_py,
                              PyObject *neigh_residue_specs_py,
                              unsigned short int atom_mask_mode,
                              int imol_map);
#endif

#ifdef __cplusplus
#ifdef USE_GUILE
//! \brief Calculate density score for a residue (Guile interface)
//!
//! Computes a numerical score indicating how well a residue fits into
//! the electron density map. Higher scores indicate better fit.
//!
//! @param imol Model molecule number
//! @param residue_spec Scheme list [chain_id, resno, ins_code]
//! @param imol_map Map molecule number
//!
//! @return Float score (typically 0.0 to 1.0+, higher is better)
//!
//! \note This uses a simplified scoring compared to map_to_model_correlation
//! \note Useful for quick assessment but map_to_model_correlation_* functions
//!       provide more detailed analysis
//!
//! Example usage:
//! \code{.scm}
//! (define score (density-score-residue-scm 1 '("A" 42 "") 2))
//! (format #t "Density score: ~a~%" score)
//! \endcode
float density_score_residue_scm(int imol, SCM residue_spec, int imol_map);
#endif
#ifdef USE_PYTHON

//! \brief Calculate density score for a residue (Python interface)
//!
//! Computes how well a residue fits into the electron density map.
//! This is a simpler alternative to the more comprehensive
//! map_to_model_correlation functions.
//!
//! @param imol Model molecule number
//! @param residue_spec Python list [chain_id, resno, ins_code]
//! @param imol_map Map molecule number
//!
//! @return Float score - higher values indicate better fit to density
//!
//! \note For comprehensive density validation, use
//!       map_to_model_correlation_stats_per_residue_range_py() instead
//!
//! Example usage:
//! \code{.py}
//! # Score a single residue
//! score = density_score_residue_py(1, ["A", 42, ""], 2)
//! print(f"Density fit score: {score:.3f}")
//!
//! # Find residues with poor density fit
//! for resno in range(1, 100):
//!     score = density_score_residue_py(1, ["A", resno, ""], 2)
//!     if score < 0.5:
//!         print(f"Poor fit: A {resno}, score = {score:.3f}")
//! \endcode
float density_score_residue_py(int imol, PyObject *residue_spec, int imol_map);
#endif
#endif

//! \brief Simple density score for given residue (C++ interface)
//!
//! Calculates a basic density fit score for the specified residue.
//! This function provides a quick assessment of how well atoms fit density.
//!
//! @param imol Model molecule number
//! @param chain_id Chain identifier
//! @param res_no Residue number
//! @param ins_code Insertion code (use "" if none)
//! @param imol_map Map molecule number to score against
//!
//! @return Float score indicating density fit quality (higher is better)
//!
//! \note This is a simplified scoring function. For detailed validation
//!       including correlation statistics, use map_to_model_correlation
//!       functions instead.
//!
//! Example:
//! \code{.cpp}
//! float score = density_score_residue(1, "A", 42, "", 2);
//! if (score < 0.5) {
//!     std::cout << "Residue A 42 has poor density fit: " << score << std::endl;
//! }
//! \endcode
float density_score_residue(int imol, const char *chain_id, int res_no, const char *ins_code, int imol_map);


#ifdef USE_GUILE
//! \brief Get the mean value of a map (Guile interface)
//!
//! Returns the mean (average) density value across all grid points in the map.
//!
//! @param imol Map molecule number
//!
//! @return SCM - Mean value as a number, or #f if imol is not a valid map
//!
//! \note Useful for understanding map scale and detecting data problems
//!
//! Example usage:
//! \code{.scm}
//! (define mean (map-mean-scm 2))
//! (format #t "Map mean: ~a~%" mean)
//! \endcode
SCM map_mean_scm(int imol);
//! \brief Get the standard deviation (sigma) of a map (Guile interface)
//!
//! Returns the standard deviation of density values in the map.
//! This is the "sigma" used for contouring at "N sigma" levels.
//!
//! @param imol Map molecule number
//!
//! @return SCM - Standard deviation as a number, or #f if invalid map
//!
//! \note The contouring level "1.5 sigma" means 1.5 times this value above the mean
//!
//! Example usage:
//! \code{.scm}
//! (define sigma (map-sigma-scm 2))
//! (format #t "Contour at 1.5 sigma = ~a~%" (* 1.5 sigma))
//! \endcode
SCM map_sigma_scm(int imol);

//! \brief Get comprehensive map statistics (Guile interface)
//!
//! Returns detailed statistical measures for the map including mean,
//! standard deviation, skew, and kurtosis.
//!
//! @param imol Map molecule number
//!
//! @return SCM - List (mean std-dev skew kurtosis) or #f if invalid map
//!
//! \note Skew and kurtosis help identify if the map has unusual distributions
//!       that might indicate problems with the data
//!
//! Example usage:
//! \code{.scm}
//! (define stats (map-statistics-scm 2))
//! (if stats
//!     (let ((mean (list-ref stats 0))
//!           (sigma (list-ref stats 1))
//!           (skew (list-ref stats 2))
//!           (kurtosis (list-ref stats 3)))
//!       (format #t "Mean: ~a, Sigma: ~a, Skew: ~a, Kurtosis: ~a~%"
//!               mean sigma skew kurtosis)))
//! \endcode
SCM map_statistics_scm(int imol);
#endif

#ifdef USE_PYTHON
//! \brief Get the mean value of a map (Python interface)
//!
//! Returns the mean (average) density value for all grid points in the map.
//!
//! @param imol Map molecule number
//!
//! @return PyObject* - Float mean value, or False if imol is not a valid map
//!
//! Example usage:
//! \code{.py}
//! mean = map_mean_py(2)
//! if mean is not False:
//!     print(f"Map mean: {mean}")
//! \endcode
PyObject *map_mean_py(int imol);
//! \brief Get the standard deviation (sigma) of a map (Python interface)
//!
//! Returns the standard deviation of density values. This is the "sigma"
//! value used when you set contouring to "1.5 sigma".
//!
//! @param imol Map molecule number
//!
//! @return PyObject* - Float sigma value, or False if invalid map
//!
//! Example usage:
//! \code{.py}
//! sigma = map_sigma_py(2)
//! if sigma is not False:
//!     print(f"1.5 sigma contour level: {1.5 * sigma}")
//! \endcode
PyObject *map_sigma_py(int imol);
//! \brief Get comprehensive map statistics (Python interface)
//!
//! Returns detailed statistical information about the map distribution.
//!
//! @param imol Map molecule number
//!
//! @return PyObject* - List [mean, std_dev, skew, kurtosis] or False if invalid
//!         - mean: Average density value
//!         - std_dev: Standard deviation (sigma)
//!         - skew: Asymmetry of the distribution
//!         - kurtosis: "Tailedness" of the distribution
//!
//! \note Normal distributions have skew≈0 and kurtosis≈3
//! \note Large deviations may indicate data problems
//!
//! Example usage:
//! \code{.py}
//! stats = map_statistics_py(2)
//! if stats is not False:
//!     mean, sigma, skew, kurtosis = stats
//!     print(f"Map statistics:")
//!     print(f"  Mean: {mean:.3f}")
//!     print(f"  Sigma: {sigma:.3f}")
//!     print(f"  Skew: {skew:.3f}")
//!     print(f"  Kurtosis: {kurtosis:.3f}")
//!
//!     if abs(skew) > 1.0:
//!         print("  Warning: Unusual skew detected")
//! \endcode
PyObject *map_statistics_py(int imol);
#endif /*USE_PYTHON */


//! \}

/*  ----------------------------------------------------------------------- */
/*                  sequence (assignment)                                   */
/*  ----------------------------------------------------------------------- */
/* section Get Sequence  */
/*! \name Get Sequence */
/* \{ */
//! \brief get the sequence for chain_id in imol
std::string get_sequence_as_fasta_for_chain(int imol, const std::string &chain_id);

//! \brief write the sequence for imol as fasta
void write_sequence(int imol, const std::string &file_name);

//! \brief trace the given map and try to apply the sequence in
//! the given pir file
void res_tracer(int imol_map, const std::string &pir_file_name);


/* \} */


/* ------------------------------------------------------------------------- */
/*                      interesting positions list                           */
/* ------------------------------------------------------------------------- */
#ifdef USE_GUILE
void register_interesting_positions_list_scm(SCM pos_list);
#endif // USE_GUILE
#ifdef USE_PYTHON
void register_interesting_positions_list_py(PyObject *pos_list);
#endif // USE_PYTHON

/* ------------------------------------------------------------------------- */
/*                      all-molecule atom overlaps                           */
/* ------------------------------------------------------------------------- */
#ifdef USE_PYTHON
//! \brief get the atom overlaps for the molecule
//!
//! @param imol the molecule index
//! @param n_max_pairs the maximum number of atom pairs to return. Typically this
//!        should be 20 or 30. Use -1 (with caution!) to get all of the
//!        (poteentially thousands) of atom overlaps.
//! @return a list of dictionaries with contact information.
//!        The list is sorted by largest overlap first.
//!        Return False on failure.
//!
PyObject *molecule_atom_overlaps_py(int imol, int n_max_pairs);
#endif // USE_PYTHON
#ifdef USE_GUILE
SCM molecule_atom_overlaps_scm(int imol);
#endif // USE_GUILE

/* ------------------------------------------------------------------------- */
/*                      prodrg import function                               */
/* ------------------------------------------------------------------------- */

// 20230916-PE Kill this
//
//! \brief import given mdl file into prodrg or other 3d generation program
//!
//! the function passed to lbg, so that it calls it when a new
//! prodrg-in.mdl file has been made.  We no longer have a timeout
//! function waiting for prodrg-in.mdl to be updated/written.
//
void prodrg_import_function(std::string file_name, std::string comp_id);



/* ------------------------------------------------------------------------- */
/*                       SBase import function                               */
/* ------------------------------------------------------------------------- */

// 20230916-PE Kill this
//
//! \brief import molecule from CCP4 SRS (or SBase, as it used to be called).
//!
//! the function passed to lbg, so that it calls it when a new
//! SBase comp_id is required.  We no longer have a timeout
//! function waiting for prodrg-in.mdl to be updated/written.
//
void sbase_import_function(std::string comp_id);

/* ------------------------------------------------------------------------- */
/*                       Alignment functions (now C++)                       */
/* ------------------------------------------------------------------------- */

//! \brief align sequence to closest chain (compare across all chains
//!   in all molecules).
//!
//! Typically match_fraction is 0.95 or so.
//!
//! Return the molecule number and chain id if successful, return -1 as the
//! molecule number if not.
//!
std::pair<int, std::string>
align_to_closest_chain(std::string target_seq, float match_fraction);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_PYTHON
PyObject *align_to_closest_chain_py(std::string target_seq, float match_fraction);
#endif /* USE_PYTHON */
#ifdef USE_GUILE
SCM align_to_closest_chain_scm(std::string target_seq, float match_fraction);
#endif /* USE_GUILE */
#endif /* c++ */


#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
SCM spherical_density_overlap(SCM i_scm, SCM j_scm);
#endif // USE_GUILE
#endif // __cplusplus

void resolve_clashing_sidechains_by_deletion(int imol);

void resolve_clashing_sidechains_by_rebuilding(int imol);

/*  ----------------------------------------------------------------------- */
/*                  GUIL Utility Functions                                  */
/*  ----------------------------------------------------------------------- */
//! \brief make a simple text dialog.
void simple_text_dialog(const std::string &dialog_title, const std::string &text,
                        int geom_x, int geom_y);


/*  ----------------------------------------------------------------------- */
/*                  Phenix Functions                                        */
/*  ----------------------------------------------------------------------- */
//! \brief phenix GEO bonds representation
void graphics_to_phenix_geo_representation(int imol, int mode,
                                           const coot::phenix_geo::phenix_geo_bonds &g);
//! \brief phenix GEO bonds representation, read from file
void graphics_to_phenix_geo_representation(int imol, int mode,
                                           const std::string &geo_file_name);

/*  ----------------------------------------------------------------------- */
/*                  Client/Server                                        */
/*  ----------------------------------------------------------------------- */
//! \brief client/server functions
#ifdef USE_PYTHON
void set_python_draw_function(const std::string &command_string);
#endif // USE_PYTHON


/*  ----------------------------------------------------------------------- */
/*                  Pathology Plots                                         */
/*  ----------------------------------------------------------------------- */
#ifdef USE_PYTHON
PyObject *pathology_data(const std::string &mtz_file_name,
                         const std::string &fp_col,
                         const std::string &sigfp_col);
#endif // USE_PYTHON

/*  ----------------------------------------------------------------------- */
/*                  Utility Functions                                       */
/*  ----------------------------------------------------------------------- */

//! \brief encoding of ints
//!
// These functions are for storing the molecule number and (some other
// number) as an int and used with GPOINTER_TO_INT and GINT_TO_POINTER.
int encode_ints(int i1, int i2);
std::pair<int, int> decode_ints(int i);

//! \brief store username and password for the database.
void store_keyed_user_name(std::string key, std::string user_name, std::string passwd);

#ifdef USE_GUILE
std::vector<coot::residue_spec_t> scm_to_residue_specs(SCM s);
// and that backwards is scm_residue() above.
int key_sym_code_scm(SCM s_scm);
#endif // USE_GUILE
#ifdef USE_PYTHON
std::vector<coot::residue_spec_t> py_to_residue_specs(PyObject *s);
int key_sym_code_py(PyObject *po);
#endif // USE_PYTHON
#ifdef USE_GUILE
#ifdef USE_PYTHON
PyObject *scm_to_py(SCM s);
SCM py_to_scm(PyObject *o);
#endif // USE_GUILE
#endif // USE_PYTHON

#ifdef USE_GUILE
clipper::Spacegroup scm_symop_strings_to_space_group(SCM symop_string_list);
#endif

#ifdef USE_PYTHON
clipper::Spacegroup py_symop_strings_to_space_group(PyObject *symop_string_list);
#endif

//! \brief enable or diable sounds (coot needs to have been compiled with sounds of course)
void set_use_sounds(bool state);

//! \brief turn off sounds and particles and textures
void curmudgeon_mode();

//! easter egg 2023
void halloween();

void display_svg_from_file_in_a_dialog(const std::string &file_name);

void display_svg_from_string_in_a_dialog(const std::string &string, const std::string &title);

void display_pae_from_file_in_a_dialog(int imol, const std::string &file_name);

void read_interesting_places_json_file(const std::string &file_name);

//! return the section index (the middle section currently)
int setup_tomo_slider(int imol);
void tomo_section_view(int imol, int axis_id);
void set_tomo_section_view_section(int imol, int section_index);

//! set tomo picker is active
void set_tomo_picker_mode_is_active(short int state);

#ifdef USE_PYTHON
void tomo_map_analysis(int imol_map, PyObject *spot_positions);
void tomo_map_analysis_2(int imol_map, PyObject *spot_positions);
#endif

//! negative becomes positive and positive becomes negative.
//! Apply an offset so that most of the map is above zero.
//!
void reverse_map(int imol_map);

void read_positron_metadata(const std::string &z_data, const std::string &table);

#ifdef USE_PYTHON
PyObject *positron_pathway(PyObject *map_molecule_list_py, PyObject *pathway_points_py);
#endif

#ifdef USE_PYTHON
void positron_plot_py(const std::string &fn_z_csv, const std::string &fn_s_csv,
                      PyObject *base_map_index_list);
#endif

#ifdef USE_PYTHON
PyObject *global_phasing_screen(int imol, PyObject *screen_dict);
#endif


#ifdef SWIG
#else
void positron_plot_internal(const std::string &fn_z_csv, const std::string &fn_s_csv,
                            const std::vector<int> &base_map_index_list);
#endif

#endif // CC_INTERFACE_HH
