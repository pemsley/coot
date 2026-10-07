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

#include "geometry/residue-and-atom-specs.hh"
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

   //! \brief a CCP4i project path or alias entry (internal)
   //!
   //! Used when parsing the PROJECT_PATH and PROJECT_ALIAS records of a CCP4i
   //! directories file to make the list of project directories.
   //! \c index is the project index, \c s is the path (or alias) string
   //! and \c flag is the alias flag.
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

   //! \brief counts of the bond types across a PISA interface (internal)
   //!
   //! The numbers of hydrogen bonds, salt bridges, covalent bonds and
   //! disulfide bonds found in a PISA interface description.
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
   //! \brief count the bond types in a PISA interface bond list (internal, Scheme)
   //!
   //! @param bonds_info_scm a list of 3-element bond records whose first element is
   //!        one of the symbols \c h-bonds, \c salt-bridges, \c ss-bonds or \c cov-bonds.
   //!        Records that are not of length 3 are ignored.
   //! @return the counts of each bond type
   pisa_interface_bond_info_t get_pisa_interface_bond_info_scm(SCM bonds_info_scm);
   #endif
   #ifdef USE_PYTHON
   //! \brief count the bond types in a PISA interface bond list (internal, Python)
   //!
   //! @param bonds_info_py a list of 3-element bond records (lists) whose first element is
   //!        one of the strings "h-bonds", "salt-bridges", "ss-bonds" or "cov-bonds".
   //!        Records that are not of length 3 are ignored.
   //! @return the counts of each bond type
   pisa_interface_bond_info_t get_pisa_interface_bond_info_py(PyObject *bonds_info_py);
   #endif

}

//! \brief Get the Git commit hash of the Coot build
//!
//! Returns the Git commit identifier for the version of Coot currently running.
//! Useful for bug reports and version tracking. The value is recorded at
//! configure time (the output of \c git \c rev-parse \c HEAD).
//!
//! @return String containing the full Git commit hash, or
//!         "git-commit-unavailable" if the source tree was not a git repository
//!         when Coot was configured
//!
//! Example:
//! \code{.cpp}
//! std::string commit = git_commit();
//! std::cout << "Coot version: commit " << commit << std::endl;
//! \endcode
std::string git_commit();

//! \brief Filter files by glob pattern
//!
//! Returns the files in the specified directory whose names end in one of
//! the glob extensions registered for the given file type (e.g. the
//! coordinates extensions for coordinates files). Used by the file selection
//! dialogs. Not functional on MSVC builds (returns an empty vector).
//!
//! @param pre_directory Directory path to search
//! @param data_type the type of file selection: one of
//!        \c COOT_COORDS_FILE_SELECTION, \c COOT_DATASET_FILE_SELECTION,
//!        \c COOT_MAP_FILE_SELECTION, \c COOT_CIF_DICTIONARY_FILE_SELECTION,
//!        \c COOT_SAVE_COORDS_FILE_SELECTION or \c COOT_PHS_COORDS_FILE_SELECTION
//!        (other values give an empty result)
//!
//! @return sorted, de-duplicated vector of matching file names (each prefixed by
//!         \c pre_directory and "/")
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


//! \brief Add CaBLAM validation markup
//!
//! Parses a (MolProbity/Phenix) CaBLAM log file and, for each residue
//! flagged as a "CaBLAM Outlier", draws markup (pink points on the
//! carbonyl O atoms and hotpink lines to the projected CA positions) in a
//! generic display object called "xxCaBLAM" (which is cleared and reused if
//! it already exists). The log file is parsed by column position (lines of
//! length 90).
//!
//! @param imol the model molecule index
//! @param cablam_file_name the CaBLAM log file name
//! @return a vector of (residue spec, outlier level) pairs for the outlier residues.
//!         Empty if \c imol is not a valid model molecule or the file cannot be read.
std::vector<std::pair<coot::residue_spec_t, double> >
add_cablam_markup(int imol, const std::string &cablam_file_name);
#ifdef USE_GUILE
//! \brief Add CaBLAM validation markup (Scheme interface)
//!
//! As add_cablam_markup().
//!
//! @param imol the model molecule index
//! @param cablam_log_file_name the CaBLAM log file name
//! @return a list of (list residue-spec level) items, where residue-spec is
//!         (list \#t chain-id res-no ins-code). Empty list on failure.
SCM add_cablam_markup_scm(int imol, const std::string &cablam_log_file_name);
#endif
#ifdef USE_PYTHON
//! \brief Add CaBLAM validation markup (Python interface)
//!
//! As add_cablam_markup(): reads a CaBLAM log file and adds markup to
//! the graphics for the residues flagged as CaBLAM outliers.
//!
//! @param imol the model molecule index
//! @param cablam_log_file_name Path to the CaBLAM log file
//!
//! @return a list of [residue_spec, level] items, one for each CaBLAM outlier,
//!         where residue_spec is [chain_id, res_no, ins_code] and level is the
//!         number read from the outlier line of the log file.
//!         An empty list on failure.
//!
//! Example usage:
//! \code{.py}
//! results = coot.add_cablam_markup_py(0, "cablam.log")
//! print(f"Found {len(results)} CaBLAM outliers")
//! for (chain_id, res_no, ins_code), level in results:
//!     print(f"  {chain_id} {res_no}: {level}")
//! \endcode
PyObject *add_cablam_markup_py(int imol, const std::string &cablam_log_file_name);
#endif

/*  ---------------------------------------------------------------------- */
/*                       key bindings :                                    */
/*  ---------------------------------------------------------------------- */
//! \name Key Bindings
//! \{
//
//! \brief print the current key bindings to the terminal
//!
//! For each binding, writes a line to standard output with the modifier
//! (Ctrl or none), the key, the description of the action and the
//! binding type.
void print_key_bindings();
//! \}


/*  ---------------------------------------------------------------------- */
/*                       go to atom   :                                    */
/*  ---------------------------------------------------------------------- */

//! \brief set the rotation centre (internal, not for export)
//!
//! @param pos the new rotation centre (in Å, orthogonal coordinates)
void set_rotation_centre(const clipper::Coord_orth &pos);

#ifdef USE_GUILE
//! \brief Find the "next" atom for go-to-atom navigation (Scheme interface)
//!
//! This is the "next residue" step used by the Go To Atom navigation
//! (e.g. the space bar). In the current go-to-atom molecule
//! (go_to_atom_molecule_number()): if the rotation centre is close to the
//! given residue, the result is the CA (or C1' for nucleotides, or else the
//! first atom) of the following residue - moving on to the next chain (or the
//! residue after a gap) when needed and wrapping to the first atom of the
//! molecule at the end. If the rotation centre is not near the given residue,
//! the atom of the given residue itself is returned. If the given residue
//! does not exist (e.g. it was deleted), the next residue with a higher
//! residue number in that chain is used.
//!
//! This function only returns the specification - it does not move the view.
//!
//! @param chain_id Current chain identifier
//! @param resno Current residue number
//! @param ins_code Current insertion code
//! @param atom_name Current atom name (not used to choose the result)
//!
//! @return a list (chain-id res-no ins-code atom-name) for the next atom,
//!         or \#f if the go-to-atom molecule is not valid or no atom was found
SCM goto_next_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code, const char *atom_name);

//! \brief Find the "previous" atom for go-to-atom navigation (Scheme interface)
//!
//! The reverse of goto_next_atom_maybe_scm(): returns the CA (or C1', or
//! first atom) of the previous residue in the go-to-atom molecule, going back
//! to the last residue of the previous chain when needed. If the rotation
//! centre is not near the given residue, the atom of the given residue itself
//! is returned. This function does not move the view.
//!
//! @param chain_id Current chain identifier
//! @param resno Current residue number
//! @param ins_code Current insertion code
//! @param atom_name Current atom name (not used to choose the result)
//!
//! @return a list (chain-id res-no ins-code atom-name) for the previous atom,
//!         or \#f if the go-to-atom molecule is not valid or no atom was found
SCM goto_prev_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif

#ifdef USE_PYTHON
//! \brief Find the "next" atom for go-to-atom navigation (Python interface)
//!
//! This is the "next residue" step used by the Go To Atom navigation.
//! In the current go-to-atom molecule: if the rotation centre is close to the
//! given residue, the result is the CA (or C1' for nucleotides, or else the
//! first atom) of the following residue - moving on to the next chain (or the
//! residue after a gap) when needed and wrapping to the first atom of the
//! molecule at the end. If the rotation centre is not near the given residue,
//! the atom of the given residue itself is returned.
//!
//! This function only returns the specification - it does not move the view.
//!
//! @param chain_id Current chain identifier
//! @param resno Current residue number
//! @param ins_code Current insertion code (use "" if none)
//! @param atom_name Current atom name (not used to choose the result)
//!
//! @return a list [chain_id, res_no, ins_code, atom_name] for the next atom,
//!         or False if the go-to-atom molecule is not valid or no atom was found
//!
//! Example usage:
//! \code{.py}
//! next_atom = coot.goto_next_atom_maybe_py("A", 42, "", " CA ")
//! if next_atom:
//!     chain_id, res_no, ins_code, atom_name = next_atom
//! \endcode
PyObject *goto_next_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);

//! \brief Find the "previous" atom for go-to-atom navigation (Python interface)
//!
//! The reverse of goto_next_atom_maybe_py(): returns the CA (or C1', or
//! first atom) of the previous residue in the go-to-atom molecule, going back
//! to the last residue of the previous chain when needed. If the rotation
//! centre is not near the given residue, the atom of the given residue itself
//! is returned. This function does not move the view.
//!
//! @param chain_id Current chain identifier
//! @param resno Current residue number
//! @param ins_code Current insertion code (use "" if none)
//! @param atom_name Current atom name (not used to choose the result)
//!
//! @return a list [chain_id, res_no, ins_code, atom_name] for the previous atom,
//!         or False if the go-to-atom molecule is not valid or no atom was found
PyObject *goto_prev_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif

//! \brief centre the view on the given atom
//!
//! Sets the go-to-atom chain, residue, atom name and alt conf from
//! \c atom_spec and recentres on that atom in the go-to-atom molecule,
//! updating the graphics.
//!
//! @param atom_spec the atom to go to
//! @return 1 on success, 0 on failure (empty spec or atom not found)
int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec);
//! \brief centre the view on the given residue
//!
//! In the go-to-atom molecule, centres on the CA (or C1', or else the
//! first atom) of the given residue.
//!
//! @param spec the residue to go to
//! @return 1 on success, 0 on failure (invalid go-to-atom molecule or residue not found)
int set_go_to_atom_from_res_spec(const coot::residue_spec_t &spec);
#ifdef USE_GUILE
//! \brief centre the view on the given residue (Scheme interface)
//!
//! As set_go_to_atom_from_res_spec().
//!
//! @param residue_spec a residue spec, e.g. (list "A" 42 "")
//! @return 1 on success, 0 on failure
int set_go_to_atom_from_res_spec_scm(SCM residue_spec);
//! \brief centre the view on the given atom (Scheme interface)
//!
//! As set_go_to_atom_from_spec().
//!
//! @param residue_spec (despite the name) an atom spec, e.g. (list "A" 42 "" " CA " "")
//! @return 1 on success, 0 on failure
int set_go_to_atom_from_atom_spec_scm(SCM residue_spec);
#endif
#ifdef USE_PYTHON
//! \brief centre the view on the given residue (Python interface)
//!
//! As set_go_to_atom_from_res_spec().
//!
//! @param residue_spec a residue spec [chain_id, res_no, ins_code] (a 4-element
//!        spec with a leading item is also accepted)
//! @return 1 on success, 0 if the residue was not found, -1 if the spec could not be parsed
int set_go_to_atom_from_res_spec_py(PyObject *residue_spec);
//! \brief centre the view on the given atom (Python interface)
//!
//! As set_go_to_atom_from_spec().
//!
//! @param residue_spec (despite the name) an atom spec
//!        [chain_id, res_no, ins_code, atom_name, alt_conf] (a 6-element spec
//!        prefixed by the molecule number is also accepted)
//! @return 1 on success, 0 on failure
int set_go_to_atom_from_atom_spec_py(PyObject *residue_spec);
#endif


//! \brief get the active atom
//!
//! The active atom is the atom closest to the rotation centre, considering
//! all displayed, pickable model molecules.
//!
//! This is to make porting the active atom more easy for Bernhard.
//! Return a class rather than a list, and rewrite the active-residue
//! function use this atom-spec.
//!
//! @return a pair: the first value is true if an atom was found, the second is
//!         the (molecule index, atom spec) pair for that atom
std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec();
#ifdef USE_PYTHON
//! \brief Get the currently active atom (Python interface)
//!
//! The active atom is the atom closest to the rotation centre in the
//! displayed, pickable model molecules. Note that (despite the name) the
//! returned spec is a residue spec, not an atom spec.
//!
//! @return a tuple (found, (imol, residue_spec)) where:
//!         - found: True if an active atom was found, otherwise False
//!         - imol: the molecule index (-1 if not found)
//!         - residue_spec: [chain_id, res_no, ins_code] of the active atom's residue
//!
//! Example usage:
//! \code{.py}
//! found, (imol, res_spec) = coot.active_atom_spec_py()
//! if found:
//!     chain_id, res_no, ins_code = res_spec
//!     print(f"Active residue: {chain_id} {res_no} in molecule {imol}")
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
//! Return as a list of strings the symmetry operators of the
//! given molecule (model or map). If imol is a not a valid molecule,
//! return an empty list.
//!
//! @param imol the molecule index
//! @return a list of symmetry operator strings
SCM get_symmetry(int imol);
//! \}
#endif // USE_GUILE

#ifdef USE_PYTHON
//! \name More Python Symmetry Functions
//! \{

//! \brief return the symmetry of the imolth molecule
//!
//! Return as a list of strings the symmetry operators of the
//! given molecule (model or map). If imol is a not a valid molecule,
//! return an empty list.
//!
//! @param imol the molecule index
//! @return a Python list of symmetry operator strings
PyObject *get_symmetry_py(int imol);
//! \}

#endif // USE_PYTHON

//! \name More Symmetry Functions
//! \{

//! \brief return 1 if this residue clashes with the symmetry-related
//!  atoms of the same molecule.
//!
//! Symmetry-related atoms (from all the space group operators and
//! unit cell shifts of up to +/-2 cells) within \c clash_dist of any atom of
//! the residue count as a clash.
//!
//! @param imol the molecule index
//! @param chain_id the chain id
//! @param res_no the residue number
//! @param ins_code the insertion code
//! @param clash_dist the clash distance cut-off in Å - typically 3.6
//!
//! @return 1 if there was a clash, 0 means that it did not clash
//!   (also returned if the molecule has no symmetry operators),
//!   -1 means that the residue or molecule could not be found.
int clashes_with_symmetry(int imol, const char *chain_id, int res_no, const char *ins_code,
                          float clash_dist);

//! \brief Add molecular symmetry
//!
//! You will need to know how to expand your point group molecular symmetry
//! to a set of 3x3 matrices. Call this function for every matrix.
//! The matrix and origin are appended to the molecule's list of molecular
//! symmetry matrices and the graphics are redrawn. The coordinates are not
//! changed. (Note: the drawing of the molecular symmetry copies is currently
//! disabled in the molecule's drawing code.)
//!
//! @param imol the model molecule index
//! @param r_00 rotation matrix element (row 0, column 0)
//! @param r_01 rotation matrix element (row 0, column 1)
//! @param r_02 rotation matrix element (row 0, column 2)
//! @param r_10 rotation matrix element (row 1, column 0)
//! @param r_11 rotation matrix element (row 1, column 1)
//! @param r_12 rotation matrix element (row 1, column 2)
//! @param r_20 rotation matrix element (row 2, column 0)
//! @param r_21 rotation matrix element (row 2, column 1)
//! @param r_22 rotation matrix element (row 2, column 2)
//! @param about_origin_x x coordinate (Å) of the point about which the matrix operates
//! @param about_origin_y y coordinate (Å) of the point about which the matrix operates
//! @param about_origin_z z coordinate (Å) of the point about which the matrix operates
void add_molecular_symmetry(int imol,
                            double r_00, double r_01, double r_02,
                            double r_10, double r_11, double r_12,
                            double r_20, double r_21, double r_22,
                            double about_origin_x,
                            double about_origin_y,
                            double about_origin_z);

//! \brief Add molecular symmetry from MTRIX records from file
//!
//! Often molecular symmetry is described using MTRIX cards in a PDB file header.
//! Use this function to extract and apply such molecular symmetry: each
//! MTRIX operator is added as for add_molecular_symmetry(), with the
//! origin set to half of the MTRIX translation.
//!
//! @param imol the model molecule index
//! @param file_name the PDB file containing the MTRIX records
//! @return currently always 0
int add_molecular_symmetry_from_mtrix_from_file(int imol, const std::string &file_name);

//! \brief Add molecular symmetry from the MTRIX records of the molecule's own file
//!
//! This is a convenience function for the above - where you don't need to
//! specify the PDB file name: the file from which molecule \c imol was read
//! is used (if it still exists).
//!
//! @param imol the model molecule index
//! @return currently always 0
int add_molecular_symmetry_from_mtrix_from_self_file(int imol);

//! \}

/*  ---------------------------------------------------------------------- */
/*                       map functions:                                    */
/*  ---------------------------------------------------------------------- */
//! \name Extra Map Functions
//! \{

//! \brief read MTZ file filename and from it try to make maps
//!
//! Useful for reading the output of refmac.
//!
//! If the file is an MTZ file, auto_read_make_and_draw_maps_from_mtz() is
//! used, otherwise the file is tried as an extended CNS reflection file with
//! auto_read_make_and_draw_maps_from_cns(). Extra F/phi label pairs can be
//! added with set_auto_read_column_labels().
//!
//! @param filename the reflection file name
//! @return a vector of molecule indices for the new maps (empty if the file
//!         does not exist or no maps could be made)
std::vector<int> auto_read_make_and_draw_maps(const char *filename);
//! \brief make maps from an MTZ file, trying the standard column labels
//!
//! Makes maps from an MTZ file, trying these F/phi column-label pairs in
//! turn: FWT/PHWT, 2FOFCWT/PH2FOFCWT, DELFWT/PHDELWT (difference map),
//! FOFCWT/PHFOFCWT (difference map), FDM/PHIDM, FAN/PHAN (difference map),
//! F_ano/PHI_ano (difference map), F_early-late/PHI_early-late (difference map),
//! then any user-defined pairs (from set_auto_read_column_labels()). A map is
//! also made if the file has exactly one F and one phase column, and for each
//! "<prefix>.F_phi.F"/"<prefix>.F_phi.phi" pair. If the file contains F/SIGF and
//! R-free columns these are attached to the maps as the refinement data.
//!
//! @param file_name the MTZ file name
//! @return a vector of molecule indices for the new maps
std::vector<int> auto_read_make_and_draw_maps_from_mtz(const std::string &file_name);
//! \brief make maps from an extended CNS reflection file
//!
//! Makes two maps: from the "F2" columns (a normal map) and then from the
//! "F1" columns (a difference map). Files with an ".mtz" extension are
//! rejected.
//!
//! @param file_name the CNS reflection file name
//! @return a vector of molecule indices for the new maps (empty on failure)
std::vector<int> auto_read_make_and_draw_maps_from_cns(const std::string &file_name);


//! \brief does the mtz file have the columns that we want it to have?
//!
//! Column labels are matched either in full (e.g. "/crystal/dataset/FWT") or
//! by the part after the last slash (e.g. "FWT"). Anomalous difference (D)
//! columns are also accepted for \c f_col.
//!
//! @param mtz_file_name the mtz file name
//! @param f_col desired f_col
//! @param phi_col desired phi_col
//! @param weight_col desired weight col
//! @param use_weights_flag specifies if the weight_col
//!        should be checked too
//! @return 1 if the columns were found, 0 otherwise
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
//! The new map is named after the original with " Sharpen " or " Blur " and
//! the B-factor appended, and is contoured at 5 rmsd.
//!
//! @param imol_map the map molecule index
//! @param b_factor is the B-factor (in Å²) to blur by (positive numbers blur,
//!        negative numbers sharpen)
//! @return the index of the map created by applying a b-factor
//!        to the given map. Return -1 on failure.
int sharpen_blur_map(int imol_map, float b_factor);

//! \brief make a sharpened or blurred map with resampling
//!
//! resampling factor might typically be 1.3
//!
//! blurred maps are generated by using a positive value of b_factor.
//! The new map has the contour level of the original map.
//!
//! @param imol_map the map molecule index
//! @param b_factor is the B-factor (in Å²) to blur by (positive numbers blur,
//!        negative numbers sharpen)
//! @param resample_factor the factor by which the grid sampling is made finer
//!        (values of 1.0 or more; e.g. 1.5 for a "normal" X-ray map sampling)
//! @return the index of the map created by applying a b-factor
//!        to the given map. Return -1 on failure.
int sharpen_blur_map_with_resampling(int imol_map, float b_factor, float resample_factor);

//! \brief as sharpen_blur_map_with_resampling() but run in a thread (internal GUI function)
//!
//! This (gui function) allows a progress bar, and should not be part of the documented API.
//! The new map is created when the calculation has finished and it is then
//! also set as the refinement map. No molecule index is returned.
void sharpen_blur_map_with_resampling_threaded_version(int imol_map, float b_factor, float resample_factor);

#ifdef USE_GUILE
//! \brief make many sharpened or blurred maps
//!
//! blurred maps are generated by using a positive value of b_factor.
//! One new map (named "Map Blur <b>" or "Map Sharpen <b>") is created for each
//! B-factor, contoured at the original contour level times exp(-0.02 b).
//!
//! @param imol_map the map molecule index
//! @param b_factors_list a list of B-factors (in Å²)
void multi_sharpen_blur_map_scm(int imol_map, SCM b_factors_list);
#endif

#ifdef USE_PYTHON
//! \brief make many sharpened or blurred maps
//!
//! blurred maps are generated by using a positive value of b_factor.
//! One new map (named "Map Blur <b>" or "Map Sharpen <b>") is created for each
//! B-factor, contoured at the original contour level times exp(-0.02 b).
//!
//! @param imol_map the map molecule index
//! @param b_factors_list is a list of B-factors (floats, in Å²) to blur by
//!        (positive numbers blur)
void multi_sharpen_blur_map_py(int imol_map, PyObject *b_factors_list);
#endif

#ifdef USE_PYTHON
//! \brief amplitude vs resolution data for graph
//!
//! The map is converted to structure factors and binned by resolution.
//!
//! @param mol_map the map molecule index
//! @return a list of lists, one per resolution bin: element 0 is the resolution
//!  (in reciprocal Angstroms squared, 1/d²), element 1 is the count of reflections
//!  and element 2 is the average F² in that bin. Return False if mol_map is not
//!  a valid map.
PyObject *amplitude_vs_resolution_py(int mol_map);
#endif

#ifdef USE_GUILE
//! \brief amplitude vs resolution data for graph
//!
//! The map is converted to structure factors and binned by resolution.
//!
//! @param mol_map the map molecule index
//! @return a list of (list average-F² count resolution) items, one per resolution bin,
//!  where resolution is in reciprocal Angstroms squared (1/d²). Note that the element
//!  order differs from that of amplitude_vs_resolution_py(). Return an empty list
//!  if mol_map is not a valid map.
SCM amplitude_vs_resolution_scm(int mol_map);
#endif

//! \brief Flip the hand of the map
//!
//! in case it was accidentally generated on the wrong one.
//! A new map molecule ("Map <imol_map> Flipped Hand") is created, with the
//! contour level of the original map.
//!
//! @param imol_map the map molecule index
//! @return the molecule number of the flipped map, or -1 on failure.
int flip_hand(int imol_map);

#ifndef SWIG
//! \brief zero-dose extrapolation of a series of maps (test function)
//!
//! For a series of maps (e.g. from increasing doses), each map is multiplied
//! by the mask map, converted to structure factors and, for each reflection,
//! an exponential is fitted to the amplitudes as a function of the map's
//! position in the list; the amplitude is then extrapolated to the start of the
//! series. The resulting map (labelled "zde") is created as a new (EM) map molecule.
//! The maps are presumed to have the same grid.
//!
//! @param map_number_list the map molecule indices, in dose order
//! @param imol_map_mask the molecule index of the mask map
//! @return the molecule index of the new map, or -1 if there were no valid maps
int analyse_map_point_density_change(const std::vector<int> &map_number_list, int imol_map_mask);
#endif

#ifdef USE_PYTHON
//! \brief zero-dose extrapolation of a series of maps (test function, Python interface)
//!
//! As analyse_map_point_density_change().
//!
//! @param map_number_list a list of map numbers (ints), in dose order
//! @param imol_map_mask the molecule index for the mask
//! @return the molecule index of the new map, or -1 on failure
int analyse_map_point_density_change_py(PyObject *map_number_list, int imol_map_mask);
#endif

//! \brief Go to the centre of the molecule - for Cryo-EM Molecules
//!
//! and recontour at a sensible value.
//! The centre and the suggested contour level are estimated from the map
//! density; nothing happens if that fails.
//!
//! @param imol_map the map molecule index
void go_to_map_molecule_centre(int imol_map);

//! \brief b-factor from map
//!
//! calculate structure factors and use the amplitudes to estimate
//! the B-factor of the data using a wilson plot using a low resolution
//! limit of 4.5A.
//!
//! @param imol_map the map molecule index
//! @return the estimated B-factor, -1 when given a bad map, or 0 if there
//!         were too few data in the resolution range to fit
//!
float b_factor_from_map(int imol_map);


#ifdef USE_GUILE
//! \brief return the colour triple of the imolth map
//!
//! (e.g.: (list 0.4 0.6 0.8). If invalid imol return scheme false.
//!
//! @param imol the map molecule index
//! @return a list of the red, green and blue components (0 to 1)
SCM map_colour_components(int imol);
#endif // GUILE

#ifdef USE_PYTHON
//! \brief return the colour triple of the imolth map
//!
//! e.g.: [0.4, 0.6, 0.8]. If invalid imol return Py_False.
//!
//! @param imol the map molecule index
//! @return a list of the red, green and blue components (0 to 1) of the map colour
PyObject *map_colour_components_py(int imol);
#endif // PYTHON

//! \brief read a CCP4 map or a CNS map (despite the name)
//!
//! The new map becomes the scroll-wheel map.
//!
//! @param filename is the file name
//! @param is_diff_map_flag is either 0 or 1 denoting if this is a
//!        difference map
//! @return the molecule index of the new map. Return -1 on failure
int read_ccp4_map(const std::string &filename, int is_diff_map_flag);

//! \brief same function as above - old name for the function. Deleted from the API at some stage
//!
//! @param filename is the file name
//! @param is_diff_map_flag is either 0 or 1 denoting if this is a
//!        difference map
//! @return the molecule index of the new map. Return -1 on failure
int handle_read_ccp4_map(const std::string &filename, int is_diff_map_flag);

//! \brief this reads a EMDB bundle - I don't think they exist any more
//!
//! Reads all the "*.map" files in \c dir_name/map (as CCP4 maps) and all
//! the "*.ent" files in \c dir_name/fittedModels/PDB (as coordinates).
//!
//! @param dir_name the top directory of the bundle
//! @return currently always 0
int handle_read_emdb_data(const std::string &dir_name);

//! \brief show the "Map Partition by Chain" dialog (internal GUI function)
void show_map_partition_by_chain_dialog();

//! \brief partition a map by the chains of a model - use this function for scripting
//!
//! The map is split into one new map for each chain of the model
//! (labelled "Partioned map Chain <chain-id>"). The original map is then undisplayed.
//! This blocks until the partitioning is complete.
//!
//! @param imol_map the map molecule index
//! @param imol_model the model molecule index
//! @return a vector of the molecule indices of the new maps (empty on failure)
std::vector<int> map_partition_by_chain(int imol_map, int imol_model);

//! \brief partition a map by chain - use this function for use in the GUI
//!
//! As map_partition_by_chain(), but the calculation is run in a thread
//! (non-blocking, no results returned): the new maps are created (and the
//! original map undisplayed) when it has finished.
//!
//! @param imol_map the map molecule index
//! @param imol_model the model molecule index
void map_partition_by_chain_threaded(int imol_map, int imol_model);

//! \brief use (or not) vertex gradients for the specified map
//!
//! vertex gradients make the map look smoother but are slower
//! to calculate. The map contours are regenerated.
//!
//! @param imol the map molecule index
//! @param state 0 for no, 1 for yes
void set_use_vertex_gradients_for_map_normals(int imol, int state);

//! \brief turn on vertex gradients for the most recent map
//!
//! The map chosen is the most recent (highest-numbered) map molecule that
//! is displayed and is not a difference map.
void use_vertex_gradients_for_map_normals_for_latest_map();

//! \brief alias for the above (more canonical naming)
void set_use_vertex_gradients_for_map_normals_for_latest_map();


//! \}

#ifdef SWIG
#else

// non-SWIGable functions:

//! \brief overwrite a map with a weighted sum of other maps
//!
//! We overwrite the imol_map and we also presume that the
//! grid sampling of the contributing maps match. This makes it
//! much faster to generate than an average map.
//!
//! The map in \c imol_map is replaced by the weighted sum of the given
//! maps (sum of weight * map value at each grid point). Internal.
//!
//! @param imol_map the map molecule index to be overwritten
//! @param weighted_map_indices pairs of (map molecule index, weight)
void regen_map_internal(int imol_map, const std::vector<std::pair<int, float> > &weighted_map_indices);

//! \brief make a new weighted-sum map (internal)
//!
//! As regen_map_internal(), but the result is put in a new molecule (a copy of
//! the first map in the list).
//!
//! @param weighted_map_indices pairs of (map molecule index, weight)
//! @return the molecule index of the new map, or -1 if the list was empty
int make_weighted_map_simple_internal(const std::vector<std::pair<int, float> > &weighted_map_indices);
#endif

//! \brief colour map by other map
//!
//! maybe we need to specify other things like the colour table.
//! The contours of \c imol_map are coloured according to the density of
//! \c imol_map_used_for_colouring at each vertex.
//!
//! @param imol_map the molecule index
//! @param imol_map_used_for_colouring is the other map index
void
colour_map_by_other_map(int imol_map, int imol_map_used_for_colouring);

#ifdef USE_PYTHON
//! \brief colour a map by another map, using a colour table
//!
//! The colour_table should be a list of colours, each a list [r, g, b] (0 to 1).
//! So, if the colour table has 4 entries covering the range from 0 to 1, then
//! table_bin_start would be 0, the table_bin_size would be 0.25
//! and the colour_table list would have 4 entries covering the range 0->0.25, 0.25->0.5, 0.5->0.75, 0.75->1.0
//! Values in \c imol_map_used_for_colouring outside the table range take
//! the colour at the ends of the range.
//!
//! @param imol_map the map molecule index
//! @param imol_map_used_for_colouring the index of the map whose values determine the colours
//! @param table_bin_start the map value at the start of the colour table
//! @param table_bin_size the range of map values covered by each colour table entry
//! @param colour_table_list a list of [r, g, b] colours
void
colour_map_by_other_map_py(int imol_map, int imol_map_used_for_colouring, float table_bin_start, float table_bin_size,
                           PyObject *colour_table_list);

//! \brief export the contour triangles of a map as data for X3D
//!
//! @param imol the map molecule index
//! @return a list of 3 lists: [triangle_indices, vertices, normals], each a
//!         flat list (3 items per triangle/vertex/normal). The lists are empty
//!         if \c imol is not a valid map or there are no contour triangles.
PyObject *export_molecule_as_x3d(int imol);

#endif

//! \brief export a molecule as a Wavefront OBJ file
//!
//! For a map molecule the map contour triangles are exported. For a model
//! molecule the bonds representation is regenerated and its vertices and
//! triangles are written to the file.
//!
//! @param imol the map or model molecule index
//! @param file_name the output file name
//! @return true on success for a map molecule; false on failure (note: the
//!         model-molecule export currently always returns false)
bool export_molecule_as_obj(int imol, const std::string &file_name);

//! \brief export a molecule as a glTF file
//!
//! For a map molecule the map contour mesh is exported, for a model molecule
//! the first of its meshes (e.g. a ribbon representation) is exported (nothing
//! is written if the model has no meshes). The binary format (.glb) is used
//! unless \c file_name has the extension ".gltf".
//! Later this will need a handle as extra arg.
//!
//! @param imol the map or model molecule index
//! @param file_name the output file name
//! @return false if \c imol is not a valid molecule, otherwise the export status
bool export_molecule_as_gltf(int imol, const std::string &file_name);

//! \brief turn off colour map by other map
//!
//! @param imol_map the map molecule index
void colour_map_by_other_map_turn_off(int imol_map);

//! \brief Add map caps
//!
//! Adds a "cap" (a density slice) for the refinement map at the front clipping
//! plane of the current view. Requires a valid refinement map (see
//! set_imol_refinement_map()); nothing happens otherwise.
void add_density_map_cap();

//! \brief colour meshes (e.g. Ribbon diagrams) by map
//!
//! scale might be 2 and offset 1 (for example)
//!
//! At each mesh vertex the map value v is converted to f = v * scale - offset,
//! clamped to the range 0 to 1, and the vertex colour is set from f (from
//! red for f = 0 to green for f = 1).
//!
//! @param imol_model the model molecule index (whose meshes are recoloured)
//! @param imol_map the map molecule index
//! @param scale the scale factor applied to the map value
//! @param offset the offset subtracted from the scaled map value
void recolour_mesh_by_map(int imol_model, int imol_map, float scale, float offset);


//! \name Multi-Residue Torsion
//! \{
#ifdef USE_GUILE
//! \brief fit residues
//!
//! (note: fit to the current-refinement map)
//!
//! The torsions of the given residues (treated as a single fragment) are varied
//! over \c n_trials random trials to fit the refinement map, avoiding
//! neighbouring residues. The best-fitting coordinates replace those in the model.
//!
//! @param imol the model molecule index
//! @param residues_specs_scm a list of residue specs
//! @param n_trials the number of trials
//! @return \#t if the fit was run, \#f if imol was not a valid model or there
//!         was no valid refinement map
SCM multi_residue_torsion_fit_scm(int imol, SCM residues_specs_scm, int n_trials);
#endif // GUILE
//! \brief fit residues by varying their torsions
//!
//! (note: fit to the current-refinement map)
//!
//! The torsions of the given residues (treated as a single fragment) are varied
//! over \c n_trials random trials to fit the refinement map, avoiding
//! neighbouring residues (within 8 Å). The best-fitting coordinates replace
//! those in the model. Nothing happens if there is no valid refinement map.
//!
//! @param imol the model molecule index
//! @param specs the residues to fit
//! @param n_trials the number of trials
void multi_residue_torsion_fit(int imol, const std::vector<coot::residue_spec_t> &specs, int n_trials);

#ifdef USE_PYTHON
//! \brief fit residues
//!
//! (note: fit to the current-refinement map)
//!
//! As multi_residue_torsion_fit().
//!
//! @param imol the model molecule index
//! @param residues_specs_py a list of residue specs, e.g. [["A", 42, ""], ["A", 43, ""]]
//! @param n_trials the number of trials
//! @return True if the fit was run, False if imol was not a valid model or there
//!         was no valid refinement map
PyObject *multi_residue_torsion_fit_py(int imol, PyObject *residues_specs_py, int n_trials);
#endif // PYTHON
//! \}


//! \brief import a BILD file (cylinders only)
//!
//! Reads the \c .color and \c .cylinder commands of a (Chimera) BILD file and
//! adds the cylinders as a generic display object. Other BILD commands are
//! ignored. Where should this go?
//!
//! @param file_name the BILD file name
void import_bild(const std::string &file_name);

//! \brief Use servalcat for generation of Fo-Fc maps for cryo-EM data
//!
//! The model is written to a file, then "servalcat fofc" is run (in a
//! thread, so this does not block) with the given half maps. When it has
//! finished, the map in \c imol_fofc_map is updated (overwritten) from the
//! DELFWT/PHDELWT columns of the servalcat output. Requires servalcat to be
//! in the PATH.
//!
//! @param imol_model the model molecule index
//! @param imol_fofc_map the map molecule index for the difference map; if this
//!        is not a valid map, a new (empty) map molecule is created for it
//! @param half_map_1 the file name for half-map 1
//! @param half_map_2 the file name for half-map 2
//! @param resolution in A.
void servalcat_fofc(int imol_model,
                    int imol_fofc_map, const std::string &half_map_1, const std::string &half_map_2,
                    float resolution);

//! \brief Use servalcat for refinement for cryo-EM data
//!
//! Runs "servalcat refine_spa_norefmac" in a thread (so this does not block),
//! with the input and output files in a "servalcat-refine-<name>" prefix in the
//! XDG data directory. When it has finished, the refined model is read in as
//! a new molecule.
//!
//! @param imol_model is the model molecule index
//! @param half_map_1 is the file name for the half-map-1
//! @param half_map_2 is the file name for the half-map-2
//! @param mask_map is the file name for the mask (currently not passed to servalcat)
//! @param resolution in A.
//!
void servalcat_refine(int imol_model,
                      const std::string &half_map_1, const std::string &half_map_2,
                      const std::string &mask_map, float resolution);

//! \brief Use servalcat for refinement for x-ray data.
//!
//! This blocks until refinement has completed! This function has been designed
//! with scripting in mind - not interactivity!
//!
//! This presumes that the mtz for the data has already been associated with the map
//! (i.e. the map has Fobs, SigFobs and R-free columns set).
//!
//! This presumes that CCP4 has been setup correctly before invoking Coot
//! (\c CLIBD_MON must be set to an existing directory).
//!
//! Runs "servalcat refine_xtal_norefmac -s xray"; the input and output files
//! are written to the directory "coot-servalcat" with the given prefix.
//!
//! @param imol is the model molecule index
//! @param imol_map is the map molecule index
//! @param output_prefix is the prefix for the output
//! @param keyword_pairs_json a JSON string of keyword pairs to control the refinement,
//!        either an object (e.g. {"weight": "0.5"}) or a list of pairs
//!        (e.g. [["weight", "0.5"]]). Currently only "weight" is used.
//! @return the model index of the refined molecule - or -1 on failure
int servalcat_refine_xray_with_keywords(int imol, int imol_map, const std::string &output_prefix,
                                        const std::string &keyword_pairs_json);

//! \brief Use servalcat for refinement for x-ray data (asynchronous/non-blocking version).
//!
//! As servalcat_refine_xray_with_keywords(), but this does not block - it returns
//! immediately and the refined model is read in by an idle function when the
//! refinement has completed. A refinement-progress display is shown while
//! servalcat runs. This is the function to use for interactive use.
//!
//! This presumes that the mtz for the data has already been associated with the map.
//!
//! This presumes that CCP4 has been setup correctly before invoking Coot.
//!
//! @param imol is the model molecule index
//! @param imol_map is the map molecule index
//! @param output_prefix is the prefix for the output
//! @param keyword_pairs_json a JSON string of keyword pairs to control the refinement
//!        (as for servalcat_refine_xray_with_keywords())
void servalcat_refine_xray_with_keywords_async(int imol, int imol_map, const std::string &output_prefix,
                                               const std::string &keyword_pairs_json);

//! \brief run acedrg link generation
//!
//! Writes \c acedrg_link_command to the file "acedrg-link-in.txt" and runs
//! "acedrg -L acedrg-link-in.txt -o acedrg-link-from-coot" in a thread. On
//! success the resulting link dictionary (acedrg-link-from-coot_link.cif) is
//! read; on failure a dialog is shown pointing to the log
//! (acedrg-link-generation-output.log). Requires acedrg to be in the PATH.
//!
//! @param acedrg_link_command the acedrg link instructions (the contents of an acedrg -L input file)
void
run_acedrg_link_generation(const std::string &acedrg_link_command);

//! \brief add a button to the main toolbar that runs a subprocess
//!
//! When clicked, the button runs \c subprocess_command with the arguments
//! \c arg_list in a thread and redraws the graphics when it has finished.
//!
//! run generic process - doesn't work at the moment - on_completion_args
//! is wrongly interpretted. (Currently the call of \c on_completion_function
//! is disabled, so it is not run.)
//!
//! @param button_label the label for the toolbar button
//! @param subprocess_command the command to run
//! @param arg_list a list of string arguments for the command (non-string items are ignored)
//! @param on_completion_function the Python function to be called on completion
//! @param on_completion_args the arguments for on_completion_function
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
//! \brief add an atom at the rotation centre to the active molecule
//!
//! The atom is added to the molecule of the current active atom (as for
//! "Place Atom at Pointer"), at the current rotation centre. The target
//! molecule must be displayed. If there is an existing atom too close to
//! that position, the addition is disallowed and an info dialog is shown.
//! A "Water" is added to the molecule's water chain if there is one (else
//! to a new chain); other types are added as single-atom HETATM residues.
//!
//! @param element the atom type: "Water", or an element/ion name such as
//!        "Na", "K", "I", "Cl"; "SO4" and "PO4" are also handled
//!        (as multi-atom residues)
void add_an_atom(const std::string &element);
//! \}

/*  ------------------------------------------------------------------------ */
/*                             Nudge the B-factors                           */
/*  ------------------------------------------------------------------------ */
//! \name Nudge the B-factors
//! \{
//! \brief change the B-factors of the atoms of the specified residue by a (small) amount
//!
//! \c amount is added to the B-factor of every atom in the residue (it can
//! be negative); the resulting B-factors are not allowed to fall below 2.0.
//! A backup of the molecule is made first. Nothing happens if \c imol is
//! not a valid model molecule.
//!
//! @param imol the model molecule index
//! @param residue_spec_py the residue spec, e.g. \c ["A", 42, ""]
//! @param amount the B-factor change (in Å²)
#ifdef USE_PYTHON
void nudge_the_temperature_factors_py(int imol, PyObject *residue_spec_py, float amount);
#endif

//! \brief convert AlphaFold pLDDT to crystallographic B-factors
//!
//! AlphaFold models store pLDDT confidence scores (0-100) in the
//! B-factor column. This function converts them to crystallographic
//! B-factors using the Hiranuma et al. (2021) formula:
//!
//!   RMSD = 1.5 * exp(4 * (0.7 - pLDDT/100))
//!
//!   B = (8 * pi^2 / 3) * RMSD^2
//!
//! After conversion, high-confidence regions (pLDDT ~90) get
//! B-factors of ~8 A^2, while low-confidence regions (pLDDT ~50)
//! get B-factors of ~440 A^2.
//!
//! The B-factors of the molecule are modified in place and the graphics
//! are redrawn. Nothing happens if \c imol is not a valid model molecule.
//!
//! @param imol is the molecule index of the AlphaFold model
void hiranuma_inversion(int imol);

//! \}


/*  ------------------------------------------------------------------------ */
/*                         merge fragments                                   */
/*  ------------------------------------------------------------------------ */
//! \name Merge Fragments
//! \{
//! \brief merge fragments
//!
//! Each fragment is presumed to be in its own chain. The atom selections
//! (chains) of the molecule are merged together (a backup is made first),
//! the bonds are regenerated and the validation graphs are updated.
//!
//! @param imol the model molecule index
//! @return 1 if the merge was done, 0 if \c imol is not a valid model molecule
int merge_fragments(int imol);
//! \}

/*  ------------------------------------------------------------------------ */
/*                         delete items                                      */
/*  ------------------------------------------------------------------------ */
//! \name Delete Items
//! \{
//! \brief Delete Items

//! \brief delete the chain
//!
//! The chain is also removed from the geometry graphs and the Go To Atom
//! window is updated. If deleting the chain leaves the molecule empty it is
//! removed from the Display Manager.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id of the chain to be deleted
void delete_chain(int imol, const std::string  &chain_id);

//! \brief delete the side chains in the chain
//!
//! All atoms other than the main-chain atoms and CB are deleted from every
//! residue in the chain (in all models). A backup is made first.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
void delete_sidechains_for_chain(int imol, const std::string &chain_id);

//! \}


/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */
//! \name Execute Refmac
//! \{
//! \brief execute refmac
//!
//! This builds and runs (via the Python layer) a call to
//! \c refmac.run_refmac_by_filename() with these arguments, using the
//! current number of refmac cycles (graphics_info_t::refmac_ncycles).
//!
//! @param pdb_in_filename the input coordinates file name
//! @param pdb_out_filename the output coordinates file name
//! @param mtz_in_filename the input reflection data file name
//! @param mtz_out_filename the output reflection data file name
//! @param cif_lib_filename the restraints dictionary file name, use "" for none
//! @param fobs_col_name the Fobs column label
//! @param sigfobs_col_name the sigma(Fobs) column label
//! @param r_free_col_name the R-free flag column label (used only if
//!        \c have_sensible_free_r_flag is non-zero)
//! @param have_sensible_free_r_flag non-zero if the R-free column should be passed
//! @param make_molecules_flag 1 if the output model and maps should be read in
//!        as new molecules, 0 if not (use 0 when run in a sub-thread, as
//!        creating molecules there would try to update the graphics)
//! @param refmac_count_string the refmac run counter (as a string)
//! @param swap_map_colours_post_refmac_flag if this is not 1 then
//!        \c imol_refmac_map is ignored
//! @param imol_refmac_map the map molecule whose colour should be swapped
//!        with that of the new map
//! @param diff_map_flag non-zero to also make a difference map
//! @param phase_combine_flag phase combination mode: if 1 or 2, \c phib_string
//!        and \c fom_string are passed as the phase/FOM column labels
//! @param phib_string the phase column label
//! @param fom_string the figure-of-merit column label
//! @param ccp4i_project_dir the CCP4i project directory, where the log file is
//!        written; "" for the current directory
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

//! \brief the name for refmac
//!
//! The returned name is the refmac input coordinates file name for this
//! molecule: the stripped file name of the molecule with a "_refmac<n>"
//! counter suffix and "-pre.pdb" appended, e.g. "demo_refmac1-pre.pdb".
//!
//! @param imol the model molecule index (it is not checked for validity)
std::string refmac_name(int imol);

//! \}


/*  ------------------------------------------------------------------- */
/*                    file selection                                    */
/*  ------------------------------------------------------------------- */

namespace coot {
   //! \brief a file name and its modification time - for sorting files by date
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

   //! \brief trivial helper class for file attributes: a directory and its files' modification times
   class file_attribs_info_t {
   public:
      std::string directory_prefix;
      std::vector<str_mtime> file_mtimes;
   };
}


//! \brief comparison function for sorting by modification time (internal)
//!
//! @return true if \c a is more recent than \c b (so a sort puts the newest first)
bool compare_mtimes(coot::str_mtime a, coot::str_mtime b);

//! \brief parse the CCP4i project definitions file (internal)
//!
//! @param filename the CCP4i directories definitions file name
//! @return a vector of (project-name, directory) pairs. The first entry is
//!         (" - Current Dir - ", the current directory). The directory names
//!         have a trailing "/".
std::vector<std::pair<std::string, std::string> > parse_ccp4i_defs(const std::string &filename);

//! \brief return the directory of the given CCP4i project
//!
//! @param ccp4_project_name the CCP4i project name
//! @return the project directory (with a trailing "/"), or "" if the project
//!         is not found in the CCP4i definitions file
std::string ccp4_project_directory(const std::string &ccp4_project_name);

/*  -------------------------------------------------------------------- */
/*                     history                                           */
/*  -------------------------------------------------------------------- */
#include "command-arg.hh"

//! \brief add a command to the command history (internal)
//!
//! If console command display is enabled, the command is also printed to
//! the console (in scheme or python form).
//!
//! @param ls the command name followed by its arguments, as strings
void add_to_history(const std::vector<std::string> &ls);
//! \brief add a command without arguments to the command history (internal)
//!
//! @param cmd the command name (in scheme style, e.g. "display-where-is-pointer")
void add_to_history_simple(const std::string &cmd);
//! \brief add a command with typed arguments to the command history (internal)
//!
//! @param command the command name (in scheme style, e.g. "delete-chain")
//! @param args the arguments of the command
void add_to_history_typed(const std::string &command,
                          const std::vector<coot::command_arg_t> &args);
//! \brief return the string wrapped in double quotes (despite the name)
std::string single_quote(const std::string &s);
//! \brief convert a scheme-style command name to a python one (internal)
//!
//! "-" is replaced by "_"; "run-refmac-by-filename" becomes
//! "refmac.run_refmac_by_filename".
std::string pythonize_command_name(const std::string &s);
//! \brief convert a python-style command name to a scheme one (internal)
//!
//! A leading ".coot" is removed and "_" is replaced by "-".
std::string schemize_command_name(const std::string &s);
//! \brief convert the command parts to a command string in the scripting language (internal)
//!
//! Scheme form is used if Coot was built with Guile, otherwise Python.
//!
//! @param command_parts the command name followed by its arguments
//! @return the command string, or "" if there is no scripting language
std::string languagize_command(const std::vector<std::string> &command_parts);

//! \brief add the command to the MySQL session database (internal)
//!
//! This does nothing unless Coot was compiled with USE_MYSQL_DATABASE.
void add_to_database(const std::vector<std::string> &command_strings);


/*  ----------------------------------------------------------------------- */
/*                         Merge Molecules                                  */
/*  ----------------------------------------------------------------------- */
#include "api/merge-molecule-results-info-t.hh"
// return the status and vector of chain-ids of the new chain ids.
//
//! \brief merge the molecules \c add_molecules into molecule \c imol
//!
//! The added molecules (those that are valid model molecules and are not
//! \c imol itself) are undisplayed and made inactive. A single-residue
//! molecule (e.g. a ligand) may be added to an existing chain; otherwise new
//! chains (with new chain ids) are added.
//!
//! @param add_molecules the indices of the model molecules to be added
//! @param imol the model molecule index into which they are merged
//! @return a pair: the first is the status (1 for success, 0 for failure
//!         or nothing merged), the second is a vector describing the
//!         resulting merges (each item has the chain id, the residue spec
//!         and whether a whole chain was added)
std::pair<int, std::vector<merge_molecule_results_info_t> > merge_molecules_by_vector(const std::vector<int> &add_molecules, int imol);

/*  ----------------------------------------------------------------------- */
/*                         Dictionaries                                     */
/*  ----------------------------------------------------------------------- */
//! \name Dictionary Functions
//! \{

/*                  cif (geometry) dictionary                            */
//! \brief read a cif (restraints) dictionary file
//!
//! The dictionary is applied to all molecules (IMOL_ENC_ANY) and no new
//! molecule is created.
//!
//! @param filename the dictionary file name
//! @return the index of the (first) monomer read into the dictionary store,
//!         -1 on failure (so >= 0 can be treated as success).
//!         (Note: despite older documentation, this is not the number of bonds.)
int handle_cif_dictionary(const std::string &filename);
//! \brief read a cif (restraints) dictionary file - synonym for handle_cif_dictionary()
//!
//! @param filename the dictionary file name
//! @return the monomer index in the dictionary store, -1 on failure
int read_cif_dictionary(const std::string &filename);

//! \brief read a cif dictionary file - the callback used by the cif dictionary file import dialog
//!
//! If the dictionary was read successfully and
//! \c new_molecule_from_dictionary_cif_checkbutton_state is set, a new
//! model molecule of the monomer is also created - except when \c imol_enc is
//! a specific molecule, or is IMOL_ENC_AUTO and the residue type is a
//! non-auto-load ligand.
//!
//! @param filename the dictionary file name
//! @param imol_enc the model molecule number to which the dictionary applies, or
//!        IMOL_ENC_ANY = -999999, IMOL_ENC_AUTO = -999998, IMOL_ENC_UNSET = -999997
//! @param new_molecule_from_dictionary_cif_checkbutton_state 1 to generate a molecule
//!        of the monomer, 0 for not
//! @return the monomer index in the dictionary store, -1 on failure
int handle_cif_dictionary_for_molecule(const std::string &filename, int imol_enc, short int new_molecule_from_dictionary_cif_checkbutton_state);

//! \brief import a cif dictionary for the given molecule - the scripting interface
//!
//! This is handle_cif_dictionary_for_molecule() with the "generate a molecule"
//! state set to 1, so a new model molecule of the monomer may be created
//! (unless \c imol_enc is a specific molecule, or is IMOL_ENC_AUTO and the
//! residue type is a non-auto-load ligand).
//!
//! @param filename the dictionary file name
//! @param imol_enc the model molecule number, or
//!        IMOL_ENC_ANY = -999999, IMOL_ENC_AUTO = -999998, IMOL_ENC_UNSET = -999997
//! @return the monomer index in the dictionary store, -1 on failure
int read_cif_dictionary_for_molecule(const std::string &filename, int imol_enc);

//! \brief dictionary entries
//!
//! @return the comp-ids of all the monomer restraints currently in the
//!         dictionary store (there may be duplicates if a type has been
//!         read for more than one molecule)
std::vector<std::string> dictionary_entries();

//! \brief print debugging information about the dictionary store to the console
void debug_dictionary();

//! \brief get types in molecule
//!
//! @param imol the molecule index
//! @return a vector of the (unique) residue types in the molecule (empty if
//!         \c imol is not a valid model molecule)
std::vector<std::string> get_types_in_molecule(int imol);

//! \brief Get the SMILES for the given residue type
//!
//! A SMILES_CANONICAL descriptor in the dictionary is preferred, then SMILES.
//!
//! @param comp_id is the residue type
//! @return the SMILES string, or "" if the type or its SMILES descriptor is not found
std::string SMILES_for_comp_id(const std::string &comp_id);

/*! \brief return a list of all the dictionaries read */
#ifdef USE_GUILE
SCM dictionaries_read();
//! \brief return the file name of the dictionary cif file for the given residue type
//!
//! @return the file name as a string, or "" if not found
SCM cif_file_for_comp_id_scm(const std::string &comp_id);
//! \brief return a list of the comp-ids of the monomer restraints in the dictionary store
SCM dictionary_entries_scm();
//! \brief return the SMILES string for the given residue type
//!
//! @return the SMILES string ("" if not found)
SCM SMILES_for_comp_id_scm(const std::string &comp_id);
#endif // USE_GUILE


#ifdef USE_PYTHON
//! \brief return a list of the file names of all the dictionaries read
PyObject *dictionaries_read_py();
//! \brief return the file name of the dictionary cif file for the given residue type
//!
//! @param comp_id the residue type
//! @return the file name as a string, or "" if not found
PyObject *cif_file_for_comp_id_py(const std::string &comp_id);
//! \brief return a list of the comp-ids of the monomer restraints in the dictionary store
PyObject *dictionary_entries_py();

//! \brief Get the SMILES for the given residue type
//!
//! @param comp_id is the residue type
//! @return the SMILES string. If the residue type or SMILES string is not
//!         found, the empty string "" is returned (in practice False is not
//!         returned, because SMILES_for_comp_id() catches the lookup error).
PyObject *SMILES_for_comp_id_py(const std::string &comp_id);

#endif // PYTHON
//! \}


/*  ----------------------------------------------------------------------- */
/*                         Restraints                                       */
/*  ----------------------------------------------------------------------- */
//! \name  Restraints Interface
//! \{

#ifdef USE_GUILE
//! \brief return the monomer restraints for the given monomer_type,
//!       return scheme false on "restraints for monomer not found"
SCM monomer_restraints(const char *monomer_type);

//! \brief set the monomer restraints of the given monomer_type
//!
//! The restraints (in the format returned by monomer_restraints()) replace
//! those for monomer_type that apply to any molecule (IMOL_ENC_ANY).
//!
//! @return scheme true on success, scheme false on failure to set the
//!  restraints for monomer_type
SCM set_monomer_restraints(const char *monomer_type, SCM restraints);
#endif // USE_GUILE

#ifdef USE_PYTHON
//! \brief return the monomer restraints for the given monomer_type
//!
//! This is monomer_restraints_for_molecule_py() with imol = IMOL_ENC_ANY.
//!
//! @param monomer_type the residue type, e.g. "TYR"
//! @return a dictionary of the restraints or False if the restraints for the
//!         monomer were not found
PyObject *monomer_restraints_py(std::string monomer_type);
//! \brief return the monomer restraints for the given monomer_type and molecule
//!
//! The dictionary for the type is loaded (from the monomer library) if needed.
//! The returned dictionary has these keys:
//! - "_chem_comp": [comp_id, three_letter_code, name, group,
//!   number_atoms_all, number_atoms_nh, description_level]
//! - "_chem_comp_atom": list of [atom_id, type_symbol, type_energy,
//!   partial_charge, partial_charge_is_set]
//! - "_chem_comp_bond": list of [atom_id_1, atom_id_2, type, value_dist, value_esd]
//!   (the distance and esd are False if not set)
//! - "_chem_comp_angle": list of [atom_id_1, atom_id_2, atom_id_3, value_angle, value_esd]
//! - "_chem_comp_tor": list of [id, atom_id_1, atom_id_2, atom_id_3, atom_id_4,
//!   value_angle, value_esd, period]
//! - "_chem_comp_plane_atom": list of [plane_id, [atom_ids], esd]
//! - "_chem_comp_chir": list of [id, atom_id_centre, atom_id_1, atom_id_2,
//!   atom_id_3, volume_sign, esd]
//!
//! @param monomer_type the residue type
//! @param imol the model molecule index (or IMOL_ENC_ANY = -999999)
//! @return a dictionary of the restraints or False if the restraints for the
//!         monomer were not found
PyObject *monomer_restraints_for_molecule_py(std::string monomer_type, int imol);
//! \brief set the monomer restraints of the given monomer_type
//!
//! The restraints (a dictionary in the format returned by
//! monomer_restraints_py()) replace those for monomer_type that apply to any
//! molecule (IMOL_ENC_ANY).
//!
//! @param monomer_type the residue type
//! @param restraints the restraints dictionary
//! @return True on success, False on failure (e.g. \c restraints is not a dictionary)
PyObject *set_monomer_restraints_py(const char *monomer_type, PyObject *restraints);
#endif // USE_PYTHON

//! \brief show restraints editor
//!
//! Opens the restraints editor dialog for the given residue type, if its
//! restraints are in the dictionary store (and the graphics interface is in use).
//!
//! @param monomer_type the residue type
void show_restraints_editor(std::string monomer_type);

//! \brief show the restraints editor for the monomer type at the given menu index (internal, GUI callback)
//!
//! @param menu_item_index the index into the dictionary's list of monomer types
void show_restraints_editor_by_index(int menu_item_index);

//! \brief write cif restraints for monomer
//!
//! If the residue type is not found in the dictionary, a message is shown in
//! the status bar and nothing is written.
//!
//! @param monomer_type the residue type
//! @param file_name the output mmCIF file name
void write_restraints_cif_dictionary(std::string monomer_type, std::string file_name);

//! \}

/*  ----------------------------------------------------------------------- */
/*                      list nomenclature errors                            */
/*  ----------------------------------------------------------------------- */
//! \brief list the residues with nomenclature errors
//!
//! Nomenclature errors (e.g. swapped atom names in symmetric side chains)
//! are found by looking at the residue atom names and geometry (no changes
//! are made to the molecule).
//!
//! @param imol the model molecule index
//! @return a vector of (residue-type, residue-spec) pairs (empty if there are
//!         no errors or \c imol is not a valid model molecule)
std::vector<std::pair<std::string, coot::residue_spec_t> >
list_nomenclature_errors(int imol);

#ifdef USE_GUILE
//! \brief list the residues with nomenclature errors
//!
//! @return a list of residue specs (each of the form (\#t chain-id resno ins-code))
SCM list_nomenclature_errors_scm(int imol);
#endif // USE_GUILE
#ifdef USE_PYTHON
//! \brief list the residues with nomenclature errors
//!
//! @param imol the model molecule index
//! @return a list of residue specs [chain_id, res_no, ins_code] (an empty
//!         list if there are none)
PyObject *list_nomenclature_errors_py(int imol);
#endif // USE_PYTHON

//! \brief show a dialog to fix the given nomenclature errors (internal)
//!
//! Note: this is currently disabled - the function returns immediately.
void
show_fix_nomenclature_errors_gui(int imol,
                                 const std::vector<std::pair<std::string, coot::residue_spec_t> > &nomenclature_errors);

/*  ----------------------------------------------------------------------- */
/*                  dipole                                                  */
/*  ----------------------------------------------------------------------- */
// This is here because it uses a C++ class, coot::dipole
//
#ifdef USE_GUILE
//! \brief convert a dipole to a scheme list (internal)
//!
//! @return (list dipole-number (list x y z)) where (x y z) is the dipole vector
SCM dipole_to_scm(std::pair<coot::dipole, int> dp);
#endif // USE_GUILE
#ifdef USE_PYTHON
//! \brief convert a dipole to a python list (internal)
//!
//! @return [dipole_number, [x, y, z]] where [x, y, z] is the dipole vector
PyObject *dipole_to_py(std::pair<coot::dipole, int> dp);
#endif // USE_PYTHON


#ifdef USE_PYTHON
//! \brief was this Coot built with Guile?
//!
//! @return True if Coot was compiled with Guile (scheme) support, False otherwise
PyObject *coot_has_guile();
#endif

//! \brief can Coot run Lidia?
//!
//! @return true only if Coot was built with GooCanvas (HAVE_GOOCANVAS)
bool coot_can_do_lidia_p();


/* commands to run python commands from guile and vice versa */
/* we ignore return values for now */
#ifdef USE_PYTHON
//! \brief run a scheme command from python
//!
//! @param scheme_command the scheme expression
//! @return the result of the scheme command converted to a python object,
//!         or None if Coot was built without Guile
PyObject *run_scheme_command(const char *scheme_command);
#endif // USE_PYTHON
#ifdef USE_GUILE
//! \brief run a python command from scheme
//!
//! @param python_command the python expression
//! @return the result converted to a scheme object (unspecified if the
//!         result is None or Coot was built without Python)
SCM run_python_command(const char *python_command);
#endif // USE_GUILE

// This is not inside a #ifdef USE_PYTHON because we want to use it
// from the guile level and USE_PYTHON is not passed as an argument to
// swig when generating coot_wrap_guile.cc.
//
// [Consider removing safe_python_command_by_char_star() which is
// conditionally compiled].
//! \brief run the python command string with PyRun_SimpleString()
//!
//! @param python_command the python code
//! @return 0 on success, -1 if an exception was raised or Coot was built without Python
int pyrun_simple_string(const char *python_command);

#ifdef USE_GUILE
// Return a list describing a residue like that returned by
// residues-matching-criteria (list return-val chain-id resno ins-code)
// This is a library function really.  There should be somewhere else to put it.
// It doesn't need expression at the scripting level.
// return a null list on problem
//! \brief convert a residue spec to a scheme list (internal)
//!
//! @return (list \#t chain-id resno ins-code)
SCM residue_spec_to_scm(const coot::residue_spec_t &res);
#endif

#ifdef USE_PYTHON
// Return a list describing a residue like that returned by
// residues-matching-criteria [return_val, chain_id, resno, ins_code]
// This is a library function really.  There should be somewhere else to put it.
// It doesn't need expression at the scripting level.
// return a null list on problem
//! \brief convert a residue spec to a python list (internal)
//!
//! @return [chain_id, res_no, ins_code] (a 3-item list; it no longer has a
//!         leading True return value, unlike the scheme version)
PyObject *residue_spec_to_py(const coot::residue_spec_t &res);
#endif

#ifdef USE_PYTHON
//! \brief convert a residue spec to a 3-item list
//!
//! A 4-item spec (with a leading item, e.g. [True, chain_id, res_no, ins_code])
//! is reduced to [chain_id, res_no, ins_code].
//!
//! @param residue_spec_py the 3- or 4-item residue spec
//! @return [chain_id, res_no, ins_code]; if the input is not a list, the
//!         values of an unset residue spec are used
PyObject *residue_spec_make_triple_py(PyObject *residue_spec_py);
#endif // USE_PYTHON

#ifdef USE_GUILE
//! \brief convert a scheme residue spec (3 or 4 items, a 4-item spec has a leading item) to a residue spec (internal)
//!
//! @return the residue spec - unset if \c residue_in is not a list
coot::residue_spec_t residue_spec_from_scm(SCM residue_in);
#endif

//! \brief convert a python residue spec to a residue spec (internal)
//!
//! The input is [chain_id, res_no, ins_code] or [imol, chain_id, res_no, ins_code]
//! (for the 4-item form, an integer first item is stored in the spec's
//! \c int_user_data).
//!
//! @return the residue spec - unset (test with \c unset_p()) if the input is
//!         not a list of the right types
coot::residue_spec_t residue_spec_from_py(PyObject *residue_in);

// return a spec for the first residue with the given type.
// test the returned spec for unset_p().
//
//! \brief return the spec of the first residue of the given type (in the first model)
//!
//! @param imol the model molecule index
//! @param residue_type the residue type, e.g. "HEM"
//! @return the residue spec; test it with \c unset_p() - it is unset if
//!         not found or \c imol is not a valid model molecule
coot::residue_spec_t get_residue_by_type(int imol, const std::string &residue_type);

//! \brief return the specs of all the residues of the given type
//!
//! @param imol the model molecule index
//! @param residue_type the residue type
//! @return a vector of residue specs (empty if none found or \c imol is not
//!         a valid model molecule)
std::vector<coot::residue_spec_t> get_residue_specs_in_mol(int imol, const std::string &residue_type);

// Always returns a list
//! \brief return the specs of all the residues of the given type
//!
//! @param imol the model molecule index
//! @param residue_type the residue type
//! @return a list of residue specs [chain_id, res_no, ins_code] - always a
//!         list, possibly empty
PyObject *get_residue_specs_in_mol_py(int imol, const std::string &residue_type);

#ifdef USE_GUILE
// return a residue spec or scheme false
//! \brief return the spec of the first residue of the given type
//!
//! @return a residue spec or scheme false if not found
SCM get_residue_by_type_scm(int, const std::string &residue_type);
#endif

//! \brief get residue by type
//!
//! Find the first residue of the given type in the molecule (first model)
//!
//! @param imol the molecule index
//! @param residue_type the residue type requested
//! @return a residue spec [chain_id, res_no, ins_code] or Python False.
PyObject *get_residue_by_type_py(int imol, const std::string &residue_type);

//! \brief get the residue name of the specified residue
//!
//! @param imol the molecule index
//! @param residue_spec_py the residue spec
//! @return the residue name or blank on failure
std::string get_residue_name_py(int imol, PyObject *residue_spec_py);

//! \brief as above, but for use by callback
//!
//! @return the residue name or blank on failure
std::string get_residue_name(int imol, coot::residue_spec_t &res_spec);

//! \brief is the given residue the first residue in its chain? (for use by callback)
//!
//! Only the first model is checked.
//!
//! @return true if the residue is the first residue of its chain, false otherwise
//!         (or if \c imol is not a valid model molecule)
bool is_N_terminus(int imol, coot::residue_spec_t &res_spec);

//! \brief is the given residue the last residue in its chain? (for use by callback)
//!
//! @return true if the residue is the C-terminal residue of its chain, false
//!         otherwise (or if \c imol is not a valid model molecule)
bool is_C_terminus(int imol, coot::residue_spec_t &res_spec);


/*  ----------------------------------------------------------------------- */
/*               Atom info                                                  */
/*  ----------------------------------------------------------------------- */

//! \name Atom Information functions
//! \{

#ifdef USE_GUILE
//! \brief output atom info in a scheme list for use in scripting
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno the residue number
//! @param ins_code the insertion code ("" for none)
//! @param atname the atom name (PDB 4-character style, e.g. " CA ")
//! @param altconf the alt conf ("" for none)
//! @return a list in this format (list occ temp-factor element x y z).
//!         Return scheme false (\#f) if the atom is not found or imol is
//!         not a valid model molecule.
SCM atom_info_string_scm(int imol, const char *chain_id, int resno,
                         const char *ins_code, const char *atname,
                         const char *altconf);
//! \brief return the molecule as a PDB-format string
//!
//! @param imol the model molecule index
//! @return a string, or an empty list if imol is not a valid model molecule
SCM molecule_to_pdb_string_scm(int imol);
#endif // USE_GUILE

//! \brief return the residue name from a residue serial number
//!
//! The serial number is the index of the residue in the chain (starting
//! from 0) in the first model.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param serial_num the residue serial number in the chain
//! @return the residue name, or blank ("") on failure.
std::string resname_from_serial_number(int imol, const char *chain_id, int serial_num);

//! \brief return the residue name of the specified residue
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno the residue number
//! @param ins_code the insertion code ("" for none)
//! @return the residue name, or "" if the residue (or molecule) is not found
std::string residue_name(int imol, const std::string &chain_id, int resno, const std::string &ins_code);

//! \brief return the serial number of the specified residue
//!
//! The serial number is the (0-based) index of the residue in its chain.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param res_no the residue number
//! @param ins_code the insertion code ("" for none)
//! @return -1 on failure to find the residue
//
int serial_number_from_residue_specs(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);


#ifdef USE_GUILE
//! \brief Return a list of atom info for each atom in the specified residue.
//!
//! output is like this:
//! (list
//!    (list (list atom-name alt-conf)
//!          (list occ temp-fact element seg-id)
//!          (list x y z)))
//!
//! temp-fact can be a single number or a list of seven numbers (for
//! anisotropic atoms) of which the first is the isotropic B and the
//! rest are U11 U22 U33 U12 U13 U23.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno the residue number
//! @param ins_code the insertion code ("" for none)
//!
SCM residue_info(int imol, const char* chain_id, int resno, const char *ins_code);
//! \brief return the residue name of the specified residue
//!
//! @return the residue name, or scheme false (\#f) if not found
SCM residue_name_scm(int imol, const char* chain_id, int resno, const char *ins_code);

//! \brief chain fragments
//!
//! Calculates the fragment info of the chains of imol (and, if
//! screen_output_also is non-zero, writes it to the terminal).
//!
//! @return currently this always returns scheme false (\#f) - the
//!         fragment info is not converted to a scheme value.
SCM chain_fragments_scm(int imol, short int screen_output_also);

//! \brief generate a molecule from an s-expression
//!
//! A new molecule is created and the graphics are redrawn.
//!
//! @param molecule_expression the molecule as a scheme expression
//! @param name the name for the new molecule
//! @return a molecule number, -1 on error
int add_molecule(SCM molecule_expression, const char *name);

//! \brief update a molecule from a s-expression
//!
//! And going the other way, given an s-expression, update
//! molecule_number by the given molecule.  Clear what's currently
//! there first though.
//!
//! @param molecule_number the model molecule index
//! @param molecule_expression the molecule as a scheme expression
//! @return 1 on success, 0 on failure (invalid molecule or bad expression)
//!
int clear_and_update_molecule(int molecule_number, SCM molecule_expression);

//! \brief return specs of the atom close to screen centre
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf) for atom that is closest to the screen centre in any
//! displayed (and pickable) molecule. If the closest atom's residue
//! has a CA, the CA is returned. If there are multiple models with the same
//! coordinates at the screen centre, return the attributes of the atom
//! in the highest number molecule number.
//!
//! return scheme false if no active residue
//!
SCM active_residue();

//! \brief return the specs of the closest displayed atom
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf) for atom that is closest to the screen
//! centre in the displayed (and pickable) molecules (unlike
//! active-residue, potential CA substitution is not performed).
//! If there is no atom, return scheme false.
//!
SCM closest_atom_simple_scm();

//! \brief return the specs of the closest atom in imolth molecule
//!
//! Return a flat list of (list imol chain-id resno ins-code atom-name
//! alt-conf x y z) for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, no account is
//! taken of the displayed state of the molecule). If the residue of the
//! closest atom has a CA, the CA is returned.  If there is no
//! atom, or if imol is not a valid model molecule, return scheme false.
//!
SCM closest_atom(int imol);

//! \brief return the specs of the closest atom to the centre of the screen
//!
//! Return a flat list of (list imol chain-id resno ins-code atom-name
//! alt-conf x y z) for atom that is closest to the screen
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
//! @param imol the model molecule index
//! @param residue_in_scm a residue spec (list chain-id resno ins-code)
//! @param radius the distance cut-off in Å
//! @return a list of residue specs, each (list chain-id resno ins-code)
//!
SCM residues_near_residue(int imol, SCM residue_in_scm, float radius);

//! \brief return residues near the given residues
//!
//! For each of the given residues, find the residues that have atoms that are
//! closer than radius Angstroems to any atom in that residue.
//!
//! @param imol the model molecule index
//! @param residues_in a list of residue specs
//! @param radius the distance cut-off in Å
//! @return a list of (list key-residue-spec neighbour-residue-specs) items,
//!         one for each input residue that was found
//!
SCM residues_near_residues_scm(int imol, SCM residues_in, float radius);

//! \brief residues near position
//!
//! Return residue specs for residues (in the first model) that have an atom
//! closer than radius Å to the given position.
//!
//! @param imol the model molecule index (get imol from active-atom)
//! @param pos a list of 3 numbers (x y z)
//! @param radius the distance cut-off in Å
//! @return a list of residue specs (empty list on failure)
//!
SCM residues_near_position_scm(int imol, SCM pos, float radius);

//! \brief label the closest atoms in the residues that neighbour residue_spec
//!
//! For each residue that has an atom within radius Å of the specified
//! residue, label the atom of that residue that is closest to the
//! specified residue.
//!
//! @param imol the model molecule index
//! @param residue_spec_scm the central residue spec
//! @param radius the neighbour distance cut-off in Å
void label_closest_atoms_in_neighbour_residues_scm(int imol, SCM residue_spec_scm, float radius);

#endif        /* USE_GUILE */

//! \brief add hydrogens to the region around the active residue using reduce
//!
//! Find the active residue, find the near residues (within radius),
//! create a new molecule, run reduce on that, import hydrogens from
//! the result and apply them to the molecule of the active residue.
//!
//! The intermediate PDB files are written to the "coot-molprobity"
//! directory. Requires the reduce program.
//!
//! @param radius the radius (in Å) around the active residue
void hydrogenate_region(float radius);

//! \brief Add hydrogens to imol from the given pdb file
//!
//! The hydrogen atoms in the (first model of the) file are added to the
//! matching residues of imol; if a hydrogen atom of that name already
//! exists, its position is updated.
//!
//! @param imol the model molecule index
//! @param pdb_with_Hs_file_name the PDB file containing the hydrogen atoms
void add_hydrogens_from_file(int imol, std::string pdb_with_Hs_file_name);

//! \brief add hydrogen atoms to the specified residue
//!
//! Uses Coot's internal "reduce" implementation and the residue's dictionary.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param res_no the residue number
//! @param ins_code the insertion code ("" for none)
void add_hydrogen_atoms_to_residue(int imol, std::string chain_id, int res_no, std::string ins_code);

#ifdef USE_PYTHON
//! \brief add hydrogen atoms to the specified residue
//!
//! @param imol the model molecule index
//! @param residue_spec_py the residue spec, e.g. ["A", 42, ""]
void add_hydrogen_atoms_to_residue_py(int imol, PyObject *residue_spec_py);
#endif

/* Here the Python code for ATOM INFO */

//! \brief output atom info in a python list for use in scripting
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno the residue number
//! @param ins_code the insertion code ("" for none)
//! @param atname the atom name (PDB 4-character style, e.g. " CA ")
//! @param altconf the alt conf ("" for none)
//! @return a list in this format [occ, temp_factor, element, x, y, z].
//!         Return False if the atom is not found.
#ifdef USE_PYTHON
PyObject *atom_info_string_py(int imol, const char *chain_id, int resno,
                              const char *ins_code, const char *atname,
                              const char *altconf);

//! \brief return the molecule as a PDB-format string
//!
//! @param imol the model molecule index
//! @return a string, or False if imol is not a valid model molecule
PyObject *molecule_to_pdb_string_py(int imol);

//! \brief Get detailed atom information for a residue (Python interface)
//!
//! Returns per-atom information including coordinates, occupancy, B-factor,
//! and element for all atoms in the specified residue (in the first model).
//! Useful for inspecting residue completeness and identifying missing atoms.
//!
//! @param imol Model molecule index
//! @param chain_id Chain identifier (e.g., "A")
//! @param resno Residue number
//! @param ins_code Insertion code (use "" if none)
//!
//! @return PyObject* - A list of atom information, one entry per atom
//!   (False if the residue is not found):
//!   \code
//!   [
//!     [[atom_name, alt_conf], [occupancy, b_factor, element, seg_id], [x, y, z], atom_index],
//!     ...
//!   ]
//!   \endcode
//!   - \c atom_name (str): Atom name (e.g., " CA ", " SG ")
//!   - \c alt_conf (str): Alternate conformation identifier ("" if none)
//!   - \c occupancy (float): Atom occupancy (0.0-1.0)
//!   - \c b_factor (float or list of [b_iso, U11, U22, U33, U12, U13, U23]): Temperature factor
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
//! @return residue name string or False on failure
//!
PyObject *residue_name_py(int imol, const char* chain_id, int resno, const char *ins_code);

//! \brief return the centre of the specified residue
//!
//! The expanded form of this is in c-interface.h
//!
//! @param imol the model molecule index
//! @param spec_py the residue spec, e.g. ["A", 42, ""]
//! @return a list [x, y, z], or False if the residue is not found
PyObject *residue_centre_from_spec_py(int imol,
                                      PyObject *spec_py);

//! \brief chain fragments
//!
//! Calculates the fragment info of the chains of imol (and, if
//! screen_output_also is non-zero, writes it to the terminal).
//!
//! @return currently this always returns False - the fragment info is
//!         not converted to a Python value.
PyObject *chain_fragments_py(int imol, short int screen_output_also);

#ifdef USE_PYTHON
//! \brief set the B-factors of the atoms of the given residues
//!
//! All atoms of each residue are given the specified B-factor.
//! A backup is made and the bonds are regenerated.
//!
//! @param imol the model molecule index
//! @param residue_specs_b_value_tuple_list_py a list of tuples (residue_spec, b_factor),
//!        e.g. [(["A", 42, ""], 30.0), (["A", 43, ""], 35.0)]
void set_b_factor_residues_py(int imol, PyObject *residue_specs_b_value_tuple_list_py);
#endif

#ifdef USE_GUILE
//! \brief set the B-factors of the atoms of the given residues
//!
//! @param imol the model molecule index
//! @param residue_specs_b_value_tuple_list_scm a list of (list residue-spec b-factor) items
void set_b_factor_residues_scm(int imol, SCM residue_specs_b_value_tuple_list_scm);
#endif

//! \}

//! \name Using S-expression molecules
//! \{

//! \brief update a molecule from a Python expression
//!
//! Given a python-expression, update
//! molecule_number by the given molecule.  Clear what's currently
//! there first though.
//!
//! @param molecule_number the model molecule index
//! @param molecule_expression the molecule as a Python expression
//! @return 1 on success, 0 on failure
int clear_and_update_molecule_py(int molecule_number, PyObject *molecule_expression);
//! \brief generate a new molecule from a Python expression
//!
//! @param molecule_expression the molecule as a Python expression
//! @param name the name for the new molecule
//! @return a molecule number, -1 on error
int add_molecule_py(PyObject *molecule_expression, const char *name);

//! \brief return specs of the atom close to screen centre
//!
//! Return a list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf] for atom that is closest to the screen centre in the displayed
//! (and pickable) molecules. If the closest atom's residue has a CA, the
//! CA is returned. If there
//! are multiple models with the same coordinates at the screen centre,
//! return the attributes of the atom in the highest number molecule
//! number.
//!
//! @return False if no active residue
//!
//! \code{.py}
//! aa = coot.active_residue_py()
//! if aa:
//!     imol, chain_id, res_no, ins_code, atom_name, alt_conf = aa
//! \endcode
//
PyObject *active_residue_py();

//! \brief return the spec of the closest displayed atom
//!
//! @return a list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf] for atom that is closest to the screen
//! centre in the displayed (and pickable) molecules (unlike active-residue,
//! potential CA substitution is not performed).  If there is no atom,
//! return False.
//!
PyObject *closest_atom_simple_py();

//! \brief return closest atom in imolth molecule
//!
//! @param imol is the molecule index
//! @return a flat list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf, x, y, z] for the atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, no account is
//! taken of the displayed state of the imol molecule). If the residue of
//! the closest atom has a CA, the CA is returned.  If there is no
//! atom, or if imol is not a valid model molecule, return False.
//
PyObject *closest_atom_py(int imol);

//! \brief return the specs of the closest atom to the centre of the screen
//!
//! @return a flat list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf, x, y, z] for atom that is closest to the screen
//! for displayed molecules. If there is no atom, return False.
//! Don't choose the CA of the residue if there is a CA in the residue
//! of the closest atom
PyObject *closest_atom_raw_py();


//! \brief get the residues near a specified residue
//!
//! This is useful to select the residues for "Sphere" refinement.
//!
//! @param imol is the molecule index
//! @param residue_spec_in is a residue spec [chain_id, res_no, insertion_code]
//! @param radius is the cut-off distance (in Å) for atoms of the surrounding residues
//!        if they are to be included in the residue selection.
//! @return a list of residue specs for residues that have atoms that are
//! closer than radius Angstroems to any atom in the residue
//! specified by residue_in. The central residue is not included.
//
PyObject *residues_near_residue_py(int imol, PyObject *residue_spec_in, float radius);

//! \brief get the residues near a specified list of residues
//!
//! @param imol is the molecule index
//! @param residues_specs_in is a list of residue specs each of which is
//!        [chain_id, res_no, insertion_code]
//! @param radius is the cut-off distance (in Å) for atoms of the surrounding residues
//!        if they are to be included in the residue selection.
//! @return a list of [residue_spec, neighbour_residue_specs] pairs, one for each
//! input residue, where neighbour_residue_specs is the list of specs of residues
//! that have atoms closer than radius Å to any atom of residue_spec.
//! False if imol is not a valid model molecule or residues_specs_in is not a list.
//
PyObject *residues_near_residues_py(int imol, PyObject *residues_specs_in, float radius);

//! \brief get the residues near a position
//!
//! Return residue specs for residues (in the first model) that have atoms that are
//! closer than radius Angstroems to the given position.
//!
//! @param imol the model molecule index
//! @param pos_in the position as a list [x, y, z] (must be a list, not a tuple)
//! @param radius the cut-off distance in Å
//! @return a list of residue specs (empty list on failure)
//!
PyObject *residues_near_position_py(int imol, PyObject *pos_in, float radius);

//! \brief label the closest atoms in the residues that neighbour residue_spec
//!
//! For each residue that has an atom within radius Å of the specified
//! residue, label the atom of that residue that is closest to the
//! specified residue.
//!
//! @param imol the model molecule index
//! @param residue_spec_py the central residue spec, e.g. ["A", 42, ""]
//! @param radius the neighbour distance cut-off in Å
void label_closest_atoms_in_neighbour_residues_py(int imol, PyObject *residue_spec_py, float radius);

//! \brief return a Python object for the bonds
//!
//! @param imol the model molecule index
//! @return an 8-element tuple: ("atom-positions", atom_positions,
//! "bonds", bonds, "rama-goodness", rama_info, "cis-peptides", cis_peptides),
//! or False if imol is not a valid model molecule.
//
PyObject *get_bonds_representation(int imol);

//! \brief replace any current non-drawn bonds with these - and regen bonds
//!
//! Bonds to the atoms in the selection are not drawn.
//!
//! @param imol the model molecule index
//! @param cid an mmdb atom selection string; multiple selections can be
//!        combined with "||"
void set_new_non_drawn_bonds(int imol, const std::string &cid);

//! \brief add to non-drawn bonds - and regen bonds
//!
//! @param imol the model molecule index
//! @param cid an mmdb atom selection string; multiple selections can be
//!        combined with "||"
void add_to_non_drawn_bonds(int imol, const std::string &cid);

//! \brief clear the non-drawn bonds - force regen bonds to restore all
//!
//! @param imol the model molecule index
void clear_non_drawn_bonds(int imol);

//! \brief return a Python object for the radii of the atoms in the dictionary
//!
//! @return a dictionary keyed by residue type (comp-id) of the currently
//! loaded dictionary entries; each value is a dictionary of atom name to
//! van der Waals radius (in Å).
//
PyObject *get_dictionary_radii();

//! \brief return a Python object for the representation of bump and hydrogen bonds of
//!          the specified residue
//!
//! @param imol the model molecule index
//! @param residue_spec_py the residue spec, e.g. ["A", 42, ""]
//! @return the environment distances in the same tuple format as
//! get_bonds_representation(), or False if imol is not a valid model molecule.
PyObject *get_environment_distances_representation_py(int imol, PyObject *residue_spec_py);

//! \brief return a Python object for the intermediate atoms bonds
//!
//! @return the bonds of the intermediate (moving) atoms in the same tuple
//! format as get_bonds_representation(), or False if there are no
//! intermediate atoms.
//
PyObject *get_intermediate_atoms_bonds_representation();

#endif // USE_PYTHON

//! \brief return the continue-updating-refinement-atoms state
//!
//! @return 0 means off, 1 means on.
//
// Given the current wiring of the refinement, this is always 0, i.e. refine_residues()
// will return only after the atoms have finished moving.
int get_continue_updating_refinement_atoms_state();


//! \}

//! \name Status bar string functions
//! \{
//! \brief return a description of the given atom for the status bar (internal)
//!
//! @param atom_index the index of the atom in the atom selection of imol
//! @param imol the model molecule index
//! @return a string with molecule number and name, atom spec, occupancy,
//!         B-factor, element and position, or "" on failure
std::string atom_info_as_text_for_statusbar(int atom_index, int imol);
//! \brief return a description of the given symmetry atom for the status bar (internal)
//!
//! As above, but also includes the symmetry operator and cell translation.
//! (The position given is that of the original, not the symmetry-related, atom.)
std::string atom_info_as_text_for_statusbar(int atom_index, int imol,
                                            const std::pair<symm_trans_t, Cell_Translation> &sts);
//! \}


/*  ----------------------------------------------------------------------- */
/*                  Refinement                                              */
/*  ----------------------------------------------------------------------- */

//! \name Refinement with specs
//! \{

//! \brief return the specs of all the residues, each spec prefixed by the serial number
//!
//! The serial number is the (0-based) index of the residue in its chain.
//!
//! Python returns a list of [serial_number, chain_id, res_no, ins_code];
//! scheme returns a list of (serial-number \#t chain-id res-no ins-code).
//! False (\#f) is returned if imol is not a valid model molecule.
#ifdef USE_GUILE
SCM all_residues_with_serial_numbers_scm(int imol);
#endif
#ifdef USE_PYTHON
PyObject *all_residues_with_serial_numbers_py(int imol);
#endif


#ifdef SWIG
#else

//! \brief regularize the given residues
//!
//! Residues that are not found in imol are ignored.
//!
//! @param imol the model molecule index
//! @param residues the specs of the residues to regularize
//!
void regularize_residues(int imol, const std::vector<coot::residue_spec_t> &residues);
#endif

//! \brief return the MTZ file name from which the given map was made
//!
//! @param imol the map molecule index
//! @return the MTZ file name, or "" if imol is not a valid map molecule
//!         (or the map was not made from an MTZ file)
std::string mtz_file_name(int imol);

#ifdef USE_GUILE

//! \brief Refine the given residue range
//!
//! Refinement is against the map set by set_imol_refinement_map().
//! Note that the insertion codes are used only to check that the
//! residues exist; the range itself is refined with blank insertion codes.
//!
//! @return the refinement results (list info-text progress lights),
//!         or scheme false (\#f) on failure
//!
SCM refine_zone_with_full_residue_spec_scm(int imol, const char *chain_id,
                                           int resno1,
                                           const char*inscode_1,
                                           int resno2,
                                           const char*inscode_2,
                                           const char *altconf);
#endif // USE_GUILE

#ifdef USE_PYTHON
//! \brief Refine the given residue range
//!
//! Refinement is against the map set by set_imol_refinement_map().
//! Note that the insertion codes are used only to check that the
//! residues exist; the range itself is refined with blank insertion codes.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno1 the first residue number
//! @param inscode_1 the insertion code of the first residue
//! @param resno2 the last residue number
//! @param inscode_2 the insertion code of the last residue
//! @param altconf the alt conf ("" for none)
//! @return the refinement results [info_text, progress, lights] (as for
//!         accept_moving_atoms_py()), or False on failure
PyObject *refine_zone_with_full_residue_spec_py(int imol, const char *chain_id,
                                           int resno1,
                                           const char*inscode_1,
                                           int resno2,
                                           const char*inscode_2,
                                           const char *altconf);
#endif // USE_PYTHON

//! \brief set display of rotamer markup during interactive real space refinement
//!
//! @param state 1 for on, 0 for off
void set_draw_moving_atoms_rota_markup(short int state);
//! \brief set display of ramachandran markup during interactive real space refinement
//!
//! @param state 1 for on, 0 for off
void set_draw_moving_atoms_rama_markup(short int state);

//! \brief the old name for set_draw_moving_atoms_rota_markup()
void set_show_intermediate_atoms_rota_markup(short int state);
//! \brief the old name for set_draw_moving_atoms_rama_markup()
void set_show_intermediate_atoms_rama_markup(short int state);

//! \brief the getter for the rota markup state
int get_draw_moving_atoms_rota_markup_state();

//! \brief the getter for the rama markup state
int get_draw_moving_atoms_rama_markup_state();

//! \brief the old name for get_draw_moving_atoms_rota_markup_state()
int get_show_intermediate_atoms_rota_markup();

//! \brief the old name for get_draw_moving_atoms_rama_markup_state()
int get_show_intermediate_atoms_rama_markup();

//! \brief set the cryo-EM refinement flag
//!
//! Note: currently this only stores the flag; nothing else reads it.
void set_cryo_em_refinement(bool mode);
//! \brief get the cryo-EM refinement flag
bool get_cryo_em_refinement();

#ifdef USE_GUILE
//! \brief Accept refined/regularized atoms into the main molecule (scheme interface)
//!
//! Waits for any running refinement to finish first.
//!
//! @return (list info-text progress lights) or scheme false (\#f) if no
//!         restraints were found
SCM accept_moving_atoms_scm();
#endif
#ifdef USE_PYTHON
//! \brief Accept refined/regularized atoms into the main molecule (Python interface)
//!
//! Waits for any running refinement to finish first.
//!
//! When scripting refinement with set_refinement_immediate_replacement(1),
//! call this function after refinement operations to ensure atoms are
//! committed. While immediate replacement mode should handle this
//! automatically, calling accept_moving_atoms_py() ensures reliable
//! synchronization.
//!
//! @return PyObject* with one of:
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
//! \brief register a Python function to be run after the intermediate atoms have moved
//!
//! Note: currently the function is stored but not called.
void register_post_intermediate_atoms_moved_hook(PyObject *function_name);
#endif

//! \brief set whether making the moving atoms also regenerates the bonds of the
//! main molecule (without bonds to the moving atoms) - internal
//!
//! The default is true.
void set_regenerate_bonds_needs_make_bonds_type_checked(bool state);
//! \brief get the state set by set_regenerate_bonds_needs_make_bonds_type_checked() - internal
bool get_regenerate_bonds_needs_make_bonds_type_checked_state();

//! \}


/*  ----------------------------------------------------------------------- */
/*                  rigid body fitting (multiple residue ranges)            */
/*  ----------------------------------------------------------------------- */

//! \brief rigid body fit the given residue ranges as a single body
//!
//! The atoms of the ranges are fitted into the map set by
//! set_imol_refinement_map() (with the rest of the molecule masked out);
//! the result is presented as moving (intermediate) atoms.
//!
//! @param imol the model molecule index
//! @param ranges the residue ranges (chain id, start and end residue numbers)
//! @return 0 on fail to refine (no sensible place to put atoms) and 1
//! on fitting happened.
int rigid_body_fit_with_residue_ranges(int imol, const std::vector<coot::high_res_residue_range_t> &ranges);

//! \brief morph-fit all the residues of the molecule
//!
//! Model morphing (average the atom shift by using shifts of the
//! atoms within shift_average_radius A of the central residue).
//! Uses the map set by set_imol_refinement_map().
//!
//! @param imol the model molecule index
//! @param transformation_averaging_radius the averaging radius in Å
//! @return 0 on fail to move atoms and 1 on fitting happened.
//!
int morph_fit_all(int imol, float transformation_averaging_radius);

//! \brief Morph the given chain
//!
//! As morph_fit_all() but for just the given chain.
//!
//! @return 0 on failure, 1 on success
int morph_fit_chain(int imol, std::string chain_id, float transformation_averaging_radius);
#ifdef USE_GUILE
//! \brief morph the given residues (scheme interface)
int morph_fit_residues_scm(int imol, SCM residue_specs,       float transformation_averaging_radius);
#endif
#ifdef USE_PYTHON
//! \brief morph the given residues (Python interface)
//!
//! @param imol the model molecule index
//! @param residue_specs a list of residue specs, e.g. [["A", 42, ""], ["A", 43, ""]]
//! @param transformation_averaging_radius the averaging radius in Å
//! @return 0 on failure, 1 on success
int morph_fit_residues_py( int imol, PyObject *residue_specs, float transformation_averaging_radius);
#endif
//! \brief morph the given residues.
//!
//! Uses the map set by set_imol_refinement_map().
//!
//! @param imol the model molecule index
//! @param residue_specs the residues to be morphed
//! @param transformation_averaging_radius the averaging radius in Å
//! @return 0 on failure, 1 on success
int morph_fit_residues(int imol, const std::vector<coot::residue_spec_t> &residue_specs,
                       float transformation_averaging_radius);

//! \brief morph transformation are based primarily on rigid body refinement
//! of the secondary structure elements.
//!
//! Secondary structure header records are (re)calculated first.
//! Uses the map set by set_imol_refinement_map().
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @return 0 on failure, 1 on success
int morph_fit_by_secondary_structure_elements(int imol, const std::string &chain_id);


/*  ----------------------------------------------------------------------- */
/*                  check water baddies                                     */
/*  ----------------------------------------------------------------------- */

//! \brief find "bad" waters in a model molecule
//!
//! Tests the HOH/WAT residues of molecule imol against the density of the
//! refinement map (set with \c set_imol_refinement_map()), their B-factors
//! and their distances to other atoms. If the graphics interface is in use,
//! the "checked waters baddies" dialog is also shown.
//!
//! Note that this does not itself check that the refinement map is valid.
//!
//! @param imol the model molecule index
//! @param b_factor_lim waters with a B-factor above this are flagged.
//!        A negative value disables this test.
//! @param map_sigma_lim waters with density (in units of the map rmsd)
//!        below this are flagged. Pass a large negative value (e.g. -100) to
//!        disable this test.
//! @param min_dist waters closer than this (in Å) to their nearest
//!        (non-hydrogen) atom are flagged. A negative value disables this test.
//! @param max_dist waters further than this (in Å) from their nearest
//!        (non-hydrogen) atom are flagged. A negative value disables this test.
//! @param part_occ_contact_flag if non-zero, the distance tests are skipped
//! @param zero_occ_flag if non-zero, waters with zero occupancy are excluded
//!        from the distance tests
//! @param logical_operator_and_or_flag 0 means combine the criteria with
//!        logical AND, otherwise logical OR (a water is flagged if it fails
//!        any test)
//! @return a vector of the atom specs of the flagged water atoms (empty if
//!         imol is not a valid model molecule)
std::vector<coot::atom_spec_t>
check_waters_baddies(int imol, float b_factor_lim, float map_sigma_lim, float min_dist, float max_dist, short int part_occ_contact_flag, short int zero_occ_flag, short int logical_operator_and_or_flag);

//! \brief find blobs of unmodelled density
//!
//! The map is masked by the atoms of the model (mask radius 1.9 Å; waters
//! are masked or not according to \c set_find_ligand_mask_waters()) and then
//! clusters of grid points above the cut-off are found. Only clusters that are
//! too big to be a water are returned.
//!
//! @param imol_model the model molecule index
//! @param imol_map the map molecule index
//! @param cut_off_density_level the cut-off in units of the map rmsd (σ)
//! @return a vector of (blob centre, blob volume in Å^3) pairs, in order of
//!         decreasing summed density. Empty on failure.
// blobs, returning position and volume
//
std::vector<std::pair<clipper::Coord_orth, double> >
find_blobs(int imol_model, int imol_map, float cut_off_density_level);

#ifdef USE_GUILE
//! \brief find blobs of unmodelled density - scheme interface
//!
//! See \c find_blobs()
//!
//! @param imol_model the model molecule index
//! @param imol_map the map molecule index
//! @param cut_off_density_level the cut-off in units of the map rmsd (σ)
//! @return a list of blobs, each of the form (volume x y z), or \#f if
//!         either molecule index is invalid
SCM find_blobs_scm(int imol_model, int imol_map, float cut_off_density_level);
#endif

#ifdef USE_PYTHON
//! \brief Find regions of unmodelled density ("blobs") in a map
//!
//! Identifies regions of significant density that are not explained by the
//! current atomic model - e.g. missing ligands, cofactors or residues.
//!
//! The map is masked by the atoms of imol_model (mask radius 1.9 Å; whether
//! or not waters are used for masking is controlled by
//! \c set_find_ligand_mask_waters()). Clusters of grid points above the
//! cut-off are then found. Clusters that are small enough to be a water
//! (volume < 11 Å^3) are not returned - only the "big blobs" are.
//!
//! Only positive density is considered.
//!
//! @param imol_model the model molecule index. Density explained by the atoms
//!        of this model is masked out.
//! @param imol_map the map molecule index (typically a 2mFo-DFc or mFo-DFc map)
//! @param cut_off_sigma_density_level the cut-off in units of the map rmsd
//!        (σ), e.g. 1.0 for a 2mFo-DFc map or 3.0 for a difference map
//! @return a list of blobs, or False if imol_model is not a valid model
//!         molecule or imol_map is not a valid map molecule.
//!         Each blob is a 2-element list: [[x, y, z], volume], where the
//!         position is the blob centre (orthogonal coordinates in Å) and
//!         volume is the volume of the blob in Å^3. The blobs are ordered by
//!         decreasing summed density.
//!
//! \code{.py}
//! blobs = coot.find_blobs_py(0, 1, 1.0)
//! if blobs:
//!     for position, volume in blobs:
//!         print(position, volume)
//!     # go to the first (strongest) blob
//!     x, y, z = blobs[0][0]
//!     coot.set_rotation_centre(x, y, z)
//! \endcode
PyObject *find_blobs_py(int imol_model, int imol_map, float cut_off_sigma_density_level);
#endif

//! \brief show a B-factor distribution histogram of the given model molecule
//!
//! Only available if Coot was compiled with goocanvas; otherwise it does nothing.
//!
//! @param imol the model molecule index
void b_factor_distribution_graph(int imol);

/*  ----------------------------------------------------------------------- */
/*                  water chain                                             */
/*  ----------------------------------------------------------------------- */

//! \name Water Chain Functions
//! \{

#ifdef USE_GUILE
//! \brief return the chain id of the water chain from a shelx molecule.  Raw interface
//!
//! For a SHELX molecule, the water chain is taken to be the last chain of the
//! first model.
//!
//! @param imol the model molecule index
//! @return the chain id, or scheme false if no chain or bad imol
SCM water_chain_from_shelx_ins_scm(int imol);
//! \brief return the chain id of the water chain. Raw interface
//!
//! The water chain is the first chain (of the first model) that consists
//! only of HOH/WAT residues (for a SHELX molecule, see
//! \c water_chain_from_shelx_ins_scm()).
//!
//! @param imol the model molecule index
//! @return the chain id, or scheme false if there is no water chain or bad imol
SCM water_chain_scm(int imol);
#endif

#ifdef USE_PYTHON
//! \brief return the chain id of the water chain from a shelx molecule.  Raw interface
//!
//! For a SHELX molecule, the water chain is taken to be the last chain of the
//! first model.
//!
//! @param imol the model molecule index
//! @return the chain id, or False if no chain or bad imol
PyObject *water_chain_from_shelx_ins_py(int imol);
//! \brief return the chain id of the water chain. Raw interface
//!
//! The water chain is the first chain (of the first model) that consists
//! only of HOH/WAT residues (for a SHELX molecule, see
//! \c water_chain_from_shelx_ins_py()).
//!
//! @param imol the model molecule index
//! @return the chain id, or False if there is no water chain or bad imol
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

//! \brief set the logging level
//!
//! @param level is one of "LOW" (log messages are kept internally),
//!        "HIGH" (log messages are also written to the terminal) or
//!        "DEBUGGING" (as "HIGH", but debugging messages are also written).
//!        Any other value is ignored, with a warning.
void set_logging_level(const std::string &level);//!


//! \}


/*  ----------------------------------------------------------------------- */
/*                  glyco tools                                             */
/*  ----------------------------------------------------------------------- */
//! \name Glyco Tools
//! \{

//! \brief print the glycosylation tree that contains the specified residue
//!
//! The tree (rooted at an ASN, if there is one) is constructed from the
//! pyranose and ASN residues linked to the given residue.
//! Note: currently the tree is constructed, but not actually printed.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno the residue number
//! @param ins_code the insertion code
void
print_glyco_tree(int imol, const std::string &chain_id, int resno, const std::string &ins_code);

//! \brief add a named glyco tree (carbohydrate building)
//!
//! Add N-linked glycosylation starting at the given ASN residue, fitting into
//! the map in molecule imol_map.
//!
//! glycosylation_name is the type of glycosylation, one of:
//! "NAG-NAG-BMA", "high-mannose", "hybrid", "mammalian-biantennary" or
//! "plant-biantennary".
//!
//! @param imol the model molecule index
//! @param imol_map the map molecule index
//! @param glycosylation_name the type of glycosylation (see above)
//! @param chain_id the chain id of the ASN residue
//! @param res_no the residue number of the ASN residue
//! @param ins_code the insertion code of the ASN residue
void
add_named_glyco_tree(int imol, int imol_map, const std::string &glycosylation_name,
                     const std::string &chain_id, int res_no, const std::string &ins_code);

//! \}


/*  ----------------------------------------------------------------------- */
/*                  variance map                                            */
/*  ----------------------------------------------------------------------- */
//! \name Variance Map
//! \{
//! \brief Make a variance map, based on the grid of the first map.
//!
//! A new map molecule called "variance-map" is created. Invalid map
//! molecule indices in the list are ignored.
//!
//! @param map_molecule_number_vec the map molecule indices
//! @return the molecule number of the new map.  Return -1 if unable to
//!   make a variance map.
int make_variance_map(const std::vector<int> &map_molecule_number_vec);
#ifdef USE_GUILE
//! \brief Make a variance map - scheme interface
//!
//! @param map_molecule_number_list a list of map molecule indices
//! @return the molecule number of the new map, or -1 on failure
int make_variance_map_scm(SCM map_molecule_number_list);
#endif
#ifdef USE_PYTHON
//! \brief Make a variance map, based on the grid of the first map
//!
//! @param map_molecule_number_list a list of map molecule indices
//!        (non-integer items and invalid map molecules are ignored)
//! @return the molecule number of the new map, or -1 on failure
int make_variance_map_py(PyObject *map_molecule_number_list);
#endif
//! \}

/*  ----------------------------------------------------------------------- */
/*                  spin search                                             */
/*  ----------------------------------------------------------------------- */
//! \name Spin Search Functions
//! \{

//! \brief spin search - internal C++ version of \c spin_search()
//!
//! Atom names are 4-character PDB names, e.g. " CA ".
//!
//! @param imol_map the map molecule index
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno the residue number
//! @param ins_code the insertion code
//! @param direction_atoms_list the names of the two atoms that define the
//!        rotation axis
//! @param moving_atoms_list the names of the atoms to be rotated. The first
//!        of these is the one whose density is tested. Must not be empty.
void spin_search_by_atom_vectors(int imol_map, int imol, const std::string &chain_id, int resno, const std::string &ins_code, const std::pair<std::string, std::string> &direction_atoms_list, const std::vector<std::string> &moving_atoms_list);
#ifdef USE_GUILE
//! \brief for the given residue, spin the atoms in moving_atom_list
//!   around the bond defined by direction_atoms_list looking for the best
//!   fit to density of imol_map map of the first atom in
//!   moving_atom_list.  Works (only) with atoms in altconf ""
//!
//! The atoms are moved to the best-fitting position (a backup is made first).
//! Atom names are 4-character PDB names, e.g. " CA ".
//!
//! @param imol_map the map molecule index
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno the residue number
//! @param ins_code the insertion code
//! @param direction_atoms_list a list of 2 atom names
//! @param moving_atoms_list a list of atom names
void spin_search(int imol_map, int imol, const char *chain_id, int resno, const char *ins_code, SCM direction_atoms_list, SCM moving_atoms_list);
//! \brief Spin N and CB (and the rest of the side chain if extant)
//!
//!  Sometime on N-terminal addition, then N ends up pointing the wrong way.
//!  The allows us to (more or less) interchange the positions of the CB and the N.
//!  All atoms of the residue other than CA, C and O are rotated about the
//!  C-CA bond. A backup is made first.
//!
//! @param imol the model molecule index
//! @param residue_spec_scm the residue specifier
//! @param angle the rotation angle in degrees, typically 120.
void spin_N_scm(int imol, SCM residue_spec_scm, float angle);

//! \brief Spin search the density based on possible positions of CG of a side-chain.
//!
//! c.f. EM-Ringer
//!
//! See \c CG_spin_search_py(). The model is not changed.
//!
//! @param imol_model the model molecule index
//! @param imol_map the map molecule index
//! @return \#f on failure, or a list of (residue-spec delta-angle) items
SCM CG_spin_search_scm(int imol_model, int imol_map);
#endif

#ifdef USE_PYTHON
//! \brief for the given residue, spin the atoms in moving_atom_list...
//!
//!   around the bond defined by direction_atoms_list looking for the best
//!   fit to density of imol_map map of the first atom in
//!   moving_atom_list.  Works (only) with atoms in altconf ""
//!
//! The rotation is sampled in 3 degree steps and the atoms are then moved to
//! the best-fitting position (a backup is made first).
//! Atom names are 4-character PDB names, e.g. " CA ".
//!
//! @param imol_map the map molecule index
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno the residue number
//! @param ins_code the insertion code
//! @param direction_atoms_list a list of 2 atom names that define the
//!        rotation axis, e.g. [" CA ", " CB "]
//! @param moving_atoms_list a list of the atom names to be rotated, e.g.
//!        [" CG ", " CD1", " CD2"]
void spin_search_py(int imol_map, int imol, const char *chain_id, int resno, const char *ins_code, PyObject *direction_atoms_list, PyObject *moving_atoms_list);

//! \brief Spin N and CB (and the rest of the side chain if extant)
//!
//!  Sometime on N-terminal addition, then N ends up pointing the wrong way.
//!  The allows us to (more or less) interchange the positions of the CB and the N.
//!  All atoms of the residue other than CA, C and O are rotated about the
//!  C-CA bond. A backup is made first.
//!
//! @param imol is the index of the model molecule
//! @param residue_spec is the specifier for the residue, e.g. ["A", 1, ""]
//! @param angle is the rotation angle, in degrees, typically 120.
void spin_N_py(int imol, PyObject *residue_spec, float angle);

//! \brief Spin search the density based on possible positions of CG of a side-chain
//!
//! c.f. EM-Ringer. For every residue (other than PRO) in the first model
//! that has N, CA, CB and CG atoms, the N-CA-CB-CG torsion with the best
//! density fit is found and compared to the current torsion.
//! The model is not changed.
//!
//! @param imol_model is the index of the model molecule
//! @param imol_map is the index of the map molecule
//! @return either False (in the case of a failure) or a list of
//!         [delta_angle, residue_spec] pairs, where delta_angle is the
//!         (best-density torsion minus model torsion) in degrees, in the
//!         range -180 to 180.
PyObject *CG_spin_search_py(int imol_model, int imol_map);

#endif

//! \}


/*  ----------------------------------------------------------------------- */
/*                  monomer lib                                             */
/*  ----------------------------------------------------------------------- */
//! \brief search the monomer library descriptions
//!
//! The search string is split on spaces; a monomer matches if every
//! fragment is found (case-insensitively) in its description (name).
//!
//! @param search_string the text to search for
//! @param allow_minimal_descriptions_flag currently unused
//! @return a vector of (comp-id, name) pairs for the matching monomers
std::vector<std::pair<std::string, std::string> > monomer_lib_3_letter_codes_matching(const std::string &search_string, short int allow_minimal_descriptions_flag);


/*  ----------------------------------------------------------------------- */
/*                  mutate                                                  */
/*  ----------------------------------------------------------------------- */

//! \brief mutate a range of residues to the given sequence
//!
//! Each residue from res_no_start to res_no_end (no insertion codes) is
//! mutated to the corresponding residue type in target_sequence.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param res_no_start the first residue number of the range
//! @param res_no_end the last residue number of the range
//! @param target_sequence a string of single-letter residue codes; its length
//!        must match the number of residues in the range
//! @return 1 on success, 0 on failure (e.g. bad imol or length mismatch)
int mutate_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end, const std::string &target_sequence);

//! \brief mutate a residue specified by serial number - internal
//!
//! @param ires the serial number (index in the chain) of the residue
//! @param chain_id the chain id
//! @param imol the model molecule index
//! @param target_res_type the 3-letter code of the new residue type
//! @return 1 on success, 0 on failure
int mutate_internal(int ires, const char *chain_id,
                    int imol, const std::string &target_res_type);
/* a function for multimutate to make a backup and set
   have_unsaved_changes_flag themselves */

//! \brief mutate active residue to single letter code slc
//!
//! @param slc the single-letter code (case-insensitive) of the new residue type
void mutate_active_residue_to_single_letter_code(const std::string &slc);

//! \brief show keyboard mutate frame
//!
//! Internal GUI function: show the keyboard-mutate frame and give its
//! entry the focus.
void show_keyboard_mutate_frame();

//! \brief mutate by overlap
//!
//! Mutate the residue to the given type by overlapping the dictionary
//! idealized coordinates of the new type onto the current residue.
//! This works for non-standard residue types too, as long as dictionaries
//! for both the current and new types are available.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param res_no the residue number (no insertion code)
//! @param new_type the residue type (comp-id) of the new residue
//! @return 1 on success, 0 on failure
int mutate_by_overlap(int imol, const std::string &chain_id, int res_no, const std::string &new_type);


/*  ----------------------------------------------------------------------- */
/*                  ligands                                                 */
/*  ----------------------------------------------------------------------- */
//! \brief overlap ligands by graph matching - internal
//!
//! Graph-match the first residue (with atoms) of imol_ligand onto the
//! reference residue. If apply_rtop_flag is true, the best-matching ligand is
//! moved onto the reference and all other residues of imol_ligand are deleted.
//!
//! @param imol_ligand the molecule index of the moving ligand
//! @param imol_ref the molecule index of the reference
//! @param chain_id_ref the chain id of the reference residue
//! @param resno_ref the residue number of the reference residue
//! @param apply_rtop_flag whether or not to move the ligand
//! @return the graph-match info (its success member is false on failure)
coot::graph_match_info_t
overlap_ligands_internal(int imol_ligand, int imol_ref, const char *chain_id_ref,
                         int resno_ref, bool apply_rtop_flag);

//! \brief display the SMILES entry. This is the simple version - no dictionary
//! is generated.
void do_smiles_to_simple_3d_overlay_frame();

//! \brief get residues in the specified chain
//!
//! Only the first model is considered.
//!
//! @param imol the molecule index
//! @param chain_id the specified chain-id
//!
//! @return a python list of residue specs (e.g. ["A", 11, ""]) for the
//!         residues in the given chain. The list is empty if there are no
//!         such residues or imol is invalid.
PyObject *get_residues_in_chain_py(int imol, const std::string &chain_id);

//! \brief does the specified residue exist?
//!
//! @param imol the molecule index
//! @param residue_spec_py is the residue spec to test for existence
//! @return 0 for no, 1 for yes, -1 for error (bad imol)
int residue_exists_py(int imol, PyObject *residue_spec_py);

/*  ----------------------------------------------------------------------- */
/*                  conformers (part of ligand search)                      */
/*  ----------------------------------------------------------------------- */


//! \name Extra Ligand Functions
//! \{

#ifdef USE_GUILE

//! \brief make conformers of the ligand search molecules, each in its
//!  own molecule.
//!
//! As if for a ligand search: conformers are generated for the flexible
//! ligands (added with \c add_ligand_search_wiggly_ligand_molecule());
//! the number of samples is set by \c set_ligand_flexible_ligand_n_samples().
//!
//! Don't search the density.
//!
//! @return a list of new molecule numbers
SCM ligand_search_make_conformers_scm();
#endif

#ifdef USE_PYTHON
//! \brief make conformers of the ligand search molecules, each in its
//!  own molecule.
//!
//! As if for a ligand search: conformers are generated for the flexible
//! ligands (added with \c add_ligand_search_wiggly_ligand_molecule());
//! the number of samples is set by \c set_ligand_flexible_ligand_n_samples().
//! The density is not searched.
//!
//! @return a list of new molecule numbers
PyObject *ligand_search_make_conformers_py();

//! \brief get an RDKit molecule as a base64-encoded string
//!
//! The RDKit molecule is made from the residue and its dictionary, and is
//! pickled (RDKit binary format, with atom properties, so that the atom names
//! are kept) then base64-encoded. This is not a Python pickle.
//!
//! @param imol the index of the molecule
//! @param residue_spec the residue specifier, e.g. ['A', 11, ""]
//! @return the base64-encoded RDKit binary molecule. Return empty string on failure.
std::string get_rdkit_mol_base64_from_molecule(int imol, PyObject *residue_spec);

//! \brief and back the other way - import an RDKit mol in base64-encoded binary format
//!
//! A new model molecule is created from the first conformer of the RDKit
//! molecule.
//!
//! @param rdkit_mol the base64-encoded RDKit binary molecule
//! @param atom_name_list a list of atom names, one for each atom of the
//!        RDKit molecule (the length must match the number of atoms)
//! @param comp_id the residue type for the new residue
//! @return the index of the new molecule - or -1 on failure
int molecule_from_rdkit_mol_base64(const std::string &rdkit_mol, PyObject *atom_name_list, const std::string &comp_id);

//! \brief make minimal restraints (bonds and atoms) from an RDKit mol in
//! base64-encoded binary format
//!
//! The new dictionary is added to the restraints for all molecules
//! (and is also written to the file "test.cif").
//!
//! @param rdkit_mol_binary_base64 the base64-encoded RDKit binary molecule
//! @param atom_name_list_py a list of atom names, one for each atom of the
//!        RDKit molecule (ignored if the length does not match)
//! @param comp_id the residue type for the new dictionary
//! @return 1 on success, 0 on failure
// make minimal restraints from mol (bonds and atoms)
int restraints_from_rdkit_mol_base64(const std::string &rdkit_mol_binary_base64, PyObject *atom_name_list_py, const std::string &comp_id);

#endif

//! \brief make conformers of the ligand search molecules - internal
//!
//! @return a vector of the new molecule indices
std::vector<int> ligand_search_make_conformers_internal();

//! \}

/*  ----------------------------------------------------------------------- */
//                  animated ligand interactions
/*  ----------------------------------------------------------------------- */
//! \brief add an animated ligand interaction to molecule imol - internal
//!
//! @param imol the model molecule index
//! @param lb the ligand bond (interaction) to be drawn
void add_animated_ligand_interaction(int imol, const pli::fle_ligand_bond_t &lb);


/*  ----------------------------------------------------------------------- */
/*                  Cootaneer                                               */
/*  ----------------------------------------------------------------------- */
//! \brief cootaneer - internal
//!
//! @param imol_map the map molecule index
//! @param imol_model the model molecule index
//! @param atom_spec any atom in the fragment to be docked
//! @return the success status (0 is fail)
int cootaneer_internal(int imol_map, int imol_model, const coot::atom_spec_t &atom_spec);

//! \name Dock Sidechains
//! \{

#ifdef USE_GUILE
//! \brief cootaneer (i.e. assign sidechains onto mainchain model)
//!
//! atom_in_fragment_atom_spec is any atom spec in the fragment that should be
//! assigned with sidechains.
//!
//! The molecule must have a sequence assigned (e.g. using
//! \c assign_pir_sequence()); the sequence is applied only if the
//! confidence of the assignment is greater than 0.9.
//!
//! @param imol_map the map molecule index
//! @param imol_model the model molecule index
//! @param atom_in_fragment_atom_spec an atom spec of the form
//!        (chain-id res-no ins-code atom-name alt-conf)
//! @return the success status (0 is fail, -1 for a bad atom spec).
int cootaneer(int imol_map, int imol_model, SCM atom_in_fragment_atom_spec);
#endif

#ifdef USE_PYTHON
//! \brief cootaneer (i.e. assign sidechains onto mainchain model)
//!
//! atom_in_fragment_atom_spec is any atom spec in the fragment that should be
//! assigned with sidechains.
//!
//! The molecule must have a sequence assigned (e.g. using
//! \c assign_pir_sequence()); the sequence is applied only if the
//! confidence of the assignment is greater than 0.9.
//!
//! @param imol_map the map molecule index
//! @param imol_model the model molecule index
//! @param atom_in_fragment_atom_spec an atom spec, e.g. ["A", 10, "", " CA ", ""]
//!        (optionally prefixed by the molecule index)
//! @return the success status (0 is fail, -1 for a bad atom spec).
int cootaneer_py(int imol_map, int imol_model, PyObject *atom_in_fragment_atom_spec);
#endif

//! \}

/*  ----------------------------------------------------------------------- */
/*                  Sequence from Map                                       */
/*  ----------------------------------------------------------------------- */

//! \name Sequence from Map
//! \{
//! \brief guess the sequence of a fragment from the side-chain density
//!
//! Use the map to estimate the sequence - you will need a decent map
//!
//! If a sequence is found, the residues in the range are mutated to that
//! sequence and then the rotamers of the residues in the chain are
//! auto-fitted (using the refinement map).
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno_start the first residue number of the fragment
//! @param resno_end the last residue number of the fragment
//! @param imol_map the map molecule index
//! @return the guessed sequence (empty is fail).
std::string sequence_from_map(int imol, const std::string &chain_id,
                              int resno_start, int resno_end, int imol_map);

//! \brief find the best-matching sequence for a fragment and apply it
//!
//! The side-chain density of the fragment is tested against each sequence in
//! the file. If a sequence is found, the fragment is mutated, the side-chains
//! are fitted and the fragment is renumbered to match the sequence.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param resno_start the first residue number of the fragment
//! @param resno_end the last residue number of the fragment
//! @param imol_map the map molecule index
//! @param file_name_for_sequences a (multi-sequence) FASTA file
void apply_sequence_to_fragment(int imol, const std::string &chain_id, int resno_start, int resno_end,
                                int imol_map, const std::string &file_name_for_sequences);

//! \brief assign the sequence to the fragment containing the active atom
//!
//! As \c apply_sequence_to_fragment(), using the sequences associated with the
//! molecule and the refinement map.
void assign_sequence_to_active_fragment();

//! \}

/*  ----------------------------------------------------------------------- */
/*                  Generic Objects                                         */
/*  ----------------------------------------------------------------------- */

//! \brief is this (probe dots) header line something to make a graphics object from? - internal
//!
//! @param vs the 3 space-separated parts of the line
//! @return a pair: 1 and a clean name (e.g. "wide contact", "H-bonds") if this
//!         is an interesting dots object, otherwise 0 and an empty string
// return a clean name and a flag to say that this was something that
// we were interested to make a graphics object from (rather than just
// header info)
std::pair<short int, std::string> is_interesting_dots_object_next_p(const std::vector<std::string> &vs);

/*  ----------------------------------------------------------------------- */
/*                  Generic Functions                                       */
/*  ----------------------------------------------------------------------- */
#ifdef USE_GUILE
//! \brief convert a vector of strings to a scheme list - internal
SCM generic_string_vector_to_list_internal(const std::vector<std::string> &v);
//! \brief convert a vector of ints to a scheme list - internal
SCM generic_int_vector_to_list_internal(const std::vector<int> &v);
//! \brief convert a scheme list of strings to a vector of strings - internal
std::vector<std::string> generic_list_to_string_vector_internal(SCM l);
//! \brief convert an RTop to scheme - internal
//!
//! @return a list of the form ((m00 m01 m02 m10 m11 m12 m20 m21 m22) (tx ty tz))
SCM rtop_to_scm(const clipper::RTop_orth &rtop);
//! \brief invert an RTop
//!
//! @param rtop_scm an RTop of the form returned by \c rtop_to_scm()
//! @return the inverse, in the same form
SCM inverse_rtop_scm(SCM rtop_scm);
//! \brief convert a scheme atom spec to an atom spec - internal
//!
//! On failure, the string_user_data of the returned spec is "Bad Spec".
// expects an expr of length 5, ie: (list chain-id res-no ins-cod atom-name alt-conf)
coot::atom_spec_t atom_spec_from_scm_expression(SCM expr);
//! \brief convert an atom spec to scheme - internal
//!
//! @return a list of the form (int-user-data chain-id res-no ins-code atom-name alt-conf)
SCM atom_spec_to_scm(const coot::atom_spec_t &spec);
#endif        /* USE_GUILE */

#ifdef USE_PYTHON
//! \brief convert a vector of strings to a python list - internal
PyObject *generic_string_vector_to_list_internal_py(const std::vector<std::string>&v);
//! \brief convert a vector of ints to a python list - internal
PyObject *generic_int_vector_to_list_internal_py(const std::vector<int> &v);
//! \brief convert a python list of strings to a vector of strings - internal
std::vector<std::string> generic_list_to_string_vector_internal_py(PyObject *l);
//! \brief convert an RTop to python - internal
//!
//! @return a list of the form [[m00, m01, m02, m10, m11, m12, m20, m21, m22], [tx, ty, tz]]
PyObject *rtop_to_python(const clipper::RTop_orth &rtop);
//! \brief invert an RTop
//!
//! @param rtop_py an RTop of the form returned by \c rtop_to_python()
//! @return the inverse, in the same form
PyObject *inverse_rtop_py(PyObject *rtop_py);
//! \brief convert a python atom spec to an atom spec - internal
//!
//! The expression is a list [chain-id, res-no, ins-code, atom-name, alt-conf],
//! optionally prefixed by the molecule index (stored in int_user_data).
//! On failure, the string_user_data of the returned spec is "Bad Spec".
coot::atom_spec_t atom_spec_from_python_expression(PyObject *expr);
//! \brief convert an atom spec to python - internal
//!
//! @return a list of the form [int_user_data, chain_id, res_no, ins_code, atom_name, alt_conf]
PyObject *atom_spec_to_py(const coot::atom_spec_t &spec);
#endif // PYTHON

//! \brief set the state of a check button in the Display Manager - internal
//!
//! Does nothing if the graphics interface is not in use.
//!
//! @param imol the molecule index
//! @param button_type "Displayed" or "Active" ("Active" applies only to
//!        model molecules)
//! @param state the new state (1 for on, 0 for off)
void set_display_control_button_state(int imol, const std::string &button_type, int state);

//! \brief make the main window fullscreen
void fullscreen();
//! \brief undo fullscreen of the main window
void unfullscreen();

//! \brief Use left-mouse for view rotation
//!
//! @param state 1 to use the primary (left) mouse button for view rotation,
//!        0 to turn it off
void set_use_trackpad(short int state);

//! \brief this is an alias for the above
//!
//! @param state 1 to use the primary (left) mouse button for view rotation,
//!        0 to turn it off
void set_use_primary_mouse_button_for_view_rotation(short int state);



/*  ----------------------------------------------------------------------- */
/*                  LIBCURL/Download                                        */
/*  ----------------------------------------------------------------------- */

//! \brief if possible, read in the new coords getting coords via web.
//!
//! (no return value because get-url-str does not return one).
//!
//! This calls the scripting function get-ebi-pdb.
//!
//! @param code the PDB accession code
void get_coords_for_accession_code(const std::string &code);

//! \brief get the contents of a URL as a string - internal
//!
//! @param url the URL
//! @return the contents (empty, or possibly incomplete, on failure)
// internal use (strings, not binaries).
std::string coot_get_url_as_string_internal(const char *url);

//! \brief stop the download to the given file - internal
//!
//! @param file_name the name of the file being downloaded to
void stop_curl_download(const char *file_name); // stop curling the to file_name;

//! \brief get a drug molecule via Wikipedia and DrugBank
//!
//! This calls the scripting (python or scheme) get-drug-via-wikipedia function.
//!
//! @param drugname the name of the drug
//! @return the file name of the downloaded MDL mol file, or an empty string
//!         on failure
std::string get_drug_mdl_via_wikipedia_and_drugbank(std::string drugname);

//! \brief fetch and superpose AlphaFold models corresponding to model
//!
//! model must have Uniprot DBREF info in the header.
//!
//! For each chain with a UNP DBREF, the AlphaFold model is fetched, superposed
//! onto that chain and displayed as CA + ligands.
//!
//! @param imol the model molecule index
void fetch_and_superpose_alphafold_models(int imol);

//! \brief fetch the AlphaFold model for the given UniProt ID
//!
//! The model (and its PAE file, which is shown in a dialog) are downloaded
//! from the EBI AlphaFold database into the coot-download directory
//! (or are read from there, if they have already been downloaded).
//!
//! @param uniprot_id the UniProt accession code, e.g. "Q7N8I7"
//! @return the model number (-1 on failure)
int fetch_alphafold_model_for_uniprot_id(const std::string &uniprot_id);

//! \brief Loads up map from emdb
//!
//! This is an asynchronous function and will trigger a download subthread
//! and return immediately. If the map has already been downloaded (to the
//! download directory) it is read from there instead.
//!
//! @param emd_accession_code the EMDB accession code, without the "EMD-"
//!        prefix, e.g. "1234"
//!
void fetch_emdb_map(const std::string &emd_accession_code);

//! \brief fetch a COD entry
//!
//! The cif file is downloaded from the Crystallography Open Database (or read
//! from the coot-download cache directory) and read as a small-molecule cif,
//! also making maps from the reflection data it contains.
//!
//! @param cod_entry_id the COD entry id
//! @return the molecule index of the new model, or -1 on failure
int fetch_cod_entry(const std::string &cod_entry_id);


/*  ----------------------------------------------------------------------- */
/*                  Functions for FLEV layout callbacks                     */
/*  ----------------------------------------------------------------------- */
//! \brief orient the view for an interaction - internal
//!
//! Rotate the view about the screen Y axis so that the vector between the
//! two residues is more nearly perpendicular to the screen Z axis.
//! Does nothing if the graphics interface is not in use.
//!
//! @param imol the model molecule index
//! @param central_residue_spec the central residue (typically the ligand)
//! @param neighbour_residue_spec the neighbouring residue
// orient the graphics somehow so that the interaction between
// central_residue and neighbour_residue is perpendicular to screen z.
void orient_view(int imol,
                 const coot::residue_spec_t &central_residue_spec, // ligand typically
                 const coot::residue_spec_t &neighbour_residue_spec);

//! \brief return the chiral centres from topological equivalence analysis
//!
//! Note: not implemented - this currently always returns an empty vector.
//!
//! @param residue_type the residue type
/*  \brief return a list of chiral centre ids as determined from topological
    equivalence analysis based on the bond info (and element names). */
std::vector<std::string>
topological_equivalence_chiral_centres(const std::string &residue_type);


/*  ----------------------------------------------------------------------- */
/*                  New Screendump                                          */
/*  ----------------------------------------------------------------------- */

//! \brief - "save image" / "export image" / "screenshot" / "take a picture" → screendump_image()
//!
//! The scene is rendered offscreen at the framebuffer scale factor (see
//! \c set_framebuffer_scale_factor()) and written as a TGA image.
//! Unlike \c screendump_image(), ".tga" is not appended to the file name.
//!
//! @param file_name the output file name
void screendump_tga(const std::string &file_name);

//! \brief set the framebuffer scale factor
//!
//! The framebuffers are rebuilt and the graphics redrawn. This scale factor is
//! also used to make higher-resolution screendump images.
//!
//! @param sf the scale factor, e.g. 1 (normal) or 2 (double resolution)
void set_framebuffer_scale_factor(unsigned int sf);

/*  ----------------------------------------------------------------------- */
/*                  New Graphics Control                                    */
/*  ----------------------------------------------------------------------- */

//! \brief set use perspective mode
//!
//! @param state 1 for perspective projection, 0 for orthographic (the default)
void set_use_perspective_projection(short int state);

//! \brief query if perspective mode is being used
//!
//! @return 1 if perspective projection is in use, 0 if orthographic
int use_perspective_projection_state();

//! \brief set the perspective field of view
//!
//! Only has an effect when perspective projection is in use
//! (see \c set_use_perspective_projection()).
//!
//! @param degrees the field of view in degrees (default 26)
void set_perspective_fov(float degrees);

//! \brief set use ambient occlusion
//!
//! Screen-space ambient occlusion is applied only in the "fancy" rendering
//! mode (see \c set_use_fancy_lighting()).
//!
//! @param state 1 to turn on (the default), 0 to turn off
void set_use_ambient_occlusion(short int state);
//! \brief query use ambient occlusion
//!
//! @return 1 if ambient occlusion is on, 0 if off
int use_ambient_occlusion_state();

//! \brief set use depth blur
//!
//! Turn on or off depth-of-field (focus) blur.
//!
//! @param state 1 to turn on, 0 to turn off (the default)
void set_use_depth_blur(short int state);
//! \brief query use depth blur
//!
//! @return 1 if depth-of-field blur is on, 0 if off
int use_depth_blur_state();

//! \brief set use fog
//!
//! @param state 1 to turn on depth fog (the default), 0 to turn off
void set_use_fog(short int state);

//! \brief query use fog
//!
//! @return 1 if depth fog is on, 0 if off
int use_fog_state();

//! \brief set use outline
//!
//! @param state 1 to turn on the outline shader effect, 0 to turn off (the default)
void set_use_outline(short int state);

//! \brief query use outline
//!
//! @return 1 if the outline effect is on, 0 if off
int use_outline_state();

//! \brief set the map shininess
//!
//! Sets the shader shininess (specular exponent) for the map. Default 6.0.
//! Does nothing if \p imol is not a valid map molecule.
//!
//! @param imol the map molecule index
//! @param shininess the specular exponent
void set_map_shininess(int imol, float shininess);

//! \brief set the map specular strength
//!
//! Sets the shader specular strength for the map. Default 0.5.
//! Does nothing if \p imol is not a valid map molecule.
//!
//! @param imol the map molecule index
//! @param specular_strength the specular strength
void set_map_specular_strength(int imol, float specular_strength);

//! \brief set the draw state of mesh normals (for debugging)
//!
//! @param state 1 to draw the normals, 0 to not draw them
void set_draw_normals(short int state);

//! \brief query the draw state of mesh normals
//!
//! @return 1 if normals are being drawn, 0 if not
int draw_normals_state();

//! \brief set the draw state of a mesh attached to a molecule
//!
//! Does nothing if \p imol is not a valid model or map molecule or if
//! \p mesh_index is out of range.
//!
//! @param imol the molecule index (model or map)
//! @param mesh_index the index of the mesh in the molecule's list of meshes
//! @param state 1 to draw the mesh, 0 to not draw it
void set_draw_mesh(int imol, int mesh_index, short int state);

//! \brief query the draw state of a mesh attached to a molecule
//!
//! @param imol the molecule index (model or map)
//! @param mesh_index the index of the mesh in the molecule's list of meshes
//! @return 1 if the mesh is drawn, 0 if not, -1 if the mesh could not be
//!         looked up (invalid molecule or mesh index)
int draw_mesh_state(int imol, int mesh_index);

//! \brief set the default map material ambient
//!
//! The default material is applied to maps created after this call; existing
//! maps are not changed (use \c set_map_material_ambient() for those).
//! The colour components are in the range 0 to 1.
void set_default_map_material_ambient(float r, float g, float b, float alpha);

//! \brief set the default map material diffuse
//!
//! The default material is applied to maps created after this call; existing
//! maps are not changed (use \c set_map_material_diffuse() for those).
//! The colour components are in the range 0 to 1.
void set_default_map_material_diffuse(float r, float g, float b, float alpha);

//! \brief set the default model material ambient
//!
//! The default material is applied to model molecules created after this
//! call; existing models are not changed (use \c set_model_material_ambient()).
//! The colour components are in the range 0 to 1.
void set_default_model_material_ambient(float r, float g, float b, float alpha);

//! \brief set the default model material diffuse
//!
//! The default material is applied to model molecules created after this
//! call; existing models are not changed (use \c set_model_material_diffuse()).
//! The colour components are in the range 0 to 1.
void set_default_model_material_diffuse(float r, float g, float b, float alpha);

//! \brief set the default model material specularity
//!
//! Applies to model molecules created after this call.
//!
//! @param do_specularity 1 to turn specular highlights on, 0 to turn them off
//! @param specular_strength the specular strength
//! @param shininess the specular exponent
void set_default_model_material_specular(int do_specularity, float specular_strength, float shininess);

//! \brief set the map material ambient
//!
//! Does nothing if \p imol is not a valid map molecule.
//!
//! @param imol the map molecule index
//! @param r red component (0 to 1)
//! @param g green component (0 to 1)
//! @param b blue component (0 to 1)
//! @param alpha the alpha component (0 to 1)
void set_map_material_ambient(int imol, float r, float g, float b, float alpha);

//! \brief set the map material diffuse
//!
//! Does nothing if \p imol is not a valid map molecule.
//!
//! @param imol the map molecule index
//! @param r red component (0 to 1)
//! @param g green component (0 to 1)
//! @param b blue component (0 to 1)
//! @param alpha the alpha component (0 to 1)
void set_map_material_diffuse(int imol, float r, float g, float b, float alpha);

//! \brief set the specularity for a map
//!
//! This also turns specularity on for the map's material.
//! Does nothing if \p imol is not a valid map molecule.
//!
//! @param imol the map molecule index
//! @param specular_strength the specular strength
//! @param shininess the specular exponent
void set_map_material_specular(int imol, float specular_strength, float shininess);

//! \brief set the specularity for a model
//!
//! Does nothing if \p imol is not a valid model molecule.
//!
//! @param imol the model molecule index
//! @param specular_strength the specular strength
//! @param shininess the specular exponent
void set_model_material_specular(int imol, float specular_strength, float shininess);

//! \brief set the ambient material for a model
//!
//! Applied to the model's bond/atom meshes and its other meshes.
//! Does nothing if \p imol is not a valid model molecule.
//!
//! @param imol the model molecule index
//! @param r red component (0 to 1)
//! @param g green component (0 to 1)
//! @param b blue component (0 to 1)
//! @param alpha the alpha component (0 to 1)
void set_model_material_ambient(int imol, float r, float g, float b, float alpha);

//! \brief set the diffuse material for a model
//!
//! Applied to the model's bond/atom meshes and its other meshes.
//! Does nothing if \p imol is not a valid model molecule.
//!
//! @param imol the model molecule index
//! @param r red component (0 to 1)
//! @param g green component (0 to 1)
//! @param b blue component (0 to 1)
//! @param alpha the alpha component (0 to 1)
void set_model_material_diffuse(int imol, float r, float g, float b, float alpha);

//! \brief set the goodselliness (pastelization_factor) 0.3 is about right, but "the right value"
//!        depends on the renderer, may be some personal choice.
//!
//! Note that this also switches every model molecule to Goodsell-style
//! colour-by-chain mode (and redraws them).
//!
//! @param pastelization_factor the pastelization factor (default 0.3)
void set_model_goodselliness(float pastelization_factor);

//! \brief set the worm tube radius
//!
//! Sets the tube radius used for the worm representation of the given
//! model molecule. Does nothing if \c imol is not a valid model molecule.
//!
//! This needs to be set \e before the worms are calculated - it doesn't
//! do post-processing.
//!
//! @param imol the model molecule index
//! @param wtr the worm tube radius in Angstroms (the molecule default is 1.0)
void set_worm_tube_radius(int imol, float wtr);

//! \brief set the fresnel colour for a map
//!
//! Does nothing if \p imol is not a valid map molecule. The fresnel effect
//! itself is turned on with \c set_map_fresnel_settings().
//!
//! @param imol the map molecule index
//! @param red red component (0 to 1)
//! @param green green component (0 to 1)
//! @param blue blue component (0 to 1)
//! @param opacity the alpha component (0 to 1)
void set_fresnel_colour(int imol, float red, float green, float blue, float opacity);

//! \brief set the map fresnel lighting params
//!
//! The fresnel term in the map shader is
//! bias + scale * (1 + I.N)^power.
//! Does nothing if \p imol is not a valid map molecule.
//!
//! @param imol the map molecule index
//! @param state 1 to turn the fresnel effect on, 0 to turn it off (the default)
//! @param bias the fresnel bias (default 0.0)
//! @param scale the fresnel scale (default 0.3)
//! @param power the fresnel power (default 4.0)
void set_map_fresnel_settings(int imol, short int state, float bias, float scale, float power);

//! \brief hot reload map shaders
//!
//! Re-read and recompile the map and mesh shaders (map.shader and meshes.shader).
//! For shader development.
void reload_map_shader();

//! \brief hot reload model shaders
//!
//! Re-read and recompile the model shader (model.shader). For shader development.
void reload_model_shader();

//! \brief set the atom radius scale factor
//!
//! A multiplier on the radius of the atom balls (default 1.0); the bonds
//! are regenerated. Does nothing if \p imol is not a valid model molecule.
//!
//! @param imol the model molecule index
//! @param scale_factor the atom radius scale factor
void set_atom_radius_scale_factor(int imol, float scale_factor);

//! \brief set use fancy rendering lighting
//!
//! Turn on framebuffer effects (ambient occlusion, shadows etc.). When off,
//! the basic scene is rendered (the default).
//!
//! @param state where 1 mean turn on and 0 means turn off.
void set_use_fancy_lighting(short int state);

//! \brief set use simple lines for model molecule
//!
//! Applies to the model molecules currently loaded (molecules read later
//! are not affected); the bonds are regenerated where the state changes.
//!
//! @param state where 1 mean turn on and 0 means turn off.
void set_use_simple_lines_for_model_molecules(short int state);

//! \brief set the depth of the in-focus plane for depth-of-field blur
//!
//! @param z the in-focus depth in depth-buffer units (0 to 1, default 0.15)
void set_focus_blur_z_depth(float z);

//! \brief set use depth blur
//!
//! @param state 1 to turn on depth-of-field blur, 0 to turn off
void set_use_depth_blur(short int state);

//! \brief set focus blur strength
//!
//! @param st the depth-of-field blur strength (default 1.0)
void set_focus_blur_strength(float st);

//! \brief set shadow strength
//!
//! @param s is the shadow strength between 0 and 1 (default 0).
void set_shadow_strength(float s);

//! \brief set the shadow resolution
//!
//! The shadow texture is 1024 * \p reso_multiplier pixels square.
//! Values outside 1 to 7 are ignored. Default 2.
//! Equivalent to \c set_shadow_texture_resolution_multiplier().
//!
//! @param reso_multiplier the shadow texture resolution multiplier
void set_shadow_resolution(int reso_multiplier);

//! \brief set shadow box size
//!
//! @param size the shadow box size (default 120)
void set_shadow_box_size(float size);

//! \brief set SSAO kernel n samples
//!
//! The SSAO kernel samples are regenerated.
//!
//! @param n_samples the number of SSAO kernel samples (default 64)
void set_ssao_kernel_n_samples(unsigned int n_samples);

//! \brief set SSAO strength
//!
//! screen-space ambient occlusion
//!
//! @param strength is the SSAO strength between 0 and 1 (default 0.4).
void set_ssao_strength(float strength);

//! \brief set SSAO radius
//!
//! screen-space ambient occlusion
//! Doesn't do much. Not worth adjusting
//!
//! @param radius is the SSAO radius (default 30).
void set_ssao_radius(float radius);

//! \brief set SSAO bias
//!
//! screen-space ambient occlusion
//! Doesn't do much. Not worth adjusting
//!
//! @param bias the SSAO bias (default 0.02)
void set_ssao_bias(float bias);

//! \brief set SSAO blur size (0, 1, or 2)
//!
//! @param blur_size the SSAO blur size (default 1)
void set_ssao_blur_size(unsigned int blur_size);

//! \brief set the shadow softness (1, 2 or 3)
//!
//! @param softness the shadow softness (default 2)
void set_shadow_softness(unsigned int softness);

//! \brief set the shadow texture resolution multiplier
//!
//! The shadow texture is 1024 * \p m pixels square.
//! Values outside 1 to 7 are ignored. Default 2.
//!
//! @param m the shadow texture resolution multiplier
void set_shadow_texture_resolution_multiplier(unsigned int m);

//! \brief adjust the effects shader output type (for debugging effects)
//!
//! @param type 0: standard (the default), 1: input texture, 2: SSAO, 3: depth
void set_effects_shader_output_type(unsigned int type);

//! \brief adjust the effects shader brightness
//!
//! @param f the brightness (default 1.0)
void set_effects_shader_brightness(float f);

//! \brief adjust the effects shader gamma
//!
//! @param f the gamma (default 1.0)
void set_effects_shader_gamma(float f);

//! \brief set bond smoothness (default 1 (not smooth))
//!
//! Use `fac` 3 for screenshots. The bonds of all model molecules are regenerated.
//!
//! @param fac (1: coarse, 2: smooth, 3: fine)
void set_bond_smoothness_factor(unsigned int fac);

//! \brief increase bond smoothness
//!
//! and if it's currently at 3, reset back to 1. The bonds of all model
//! molecules are regenerated.
void toggle_bond_smoothness_factor();

//! \brief set the draw state of the Ramachandran plot display during Real Space Refinement
//!
//! @param state 1 to draw the in-graphics Ramachandran plot (the default), 0 to not draw it
void set_draw_gl_ramachandran_plot_during_refinement(short int state);

//! \brief set the FPS timing scale factor - default 0.002
//!
//! The scale factor converting milliseconds to screen height in the
//! frame-timing graph.
//!
//! @param f the scale factor
void set_fps_timing_scale_factor(float f);

//! \brief draw background image
//!
//! @param state true to draw the background image, false to not (the default)
void set_draw_background_image(bool state);

//! \brief internal testing function: read some test models
void read_test_gltf_models();

//! \brief load a gltf model
//!
//! If the gltf
//! files does not exist, an empty model will be created.
//! This also turns on continuous redrawing (so that models can be animated).
//!
//! @param gltf_file_name is the name of the gltf file to load
//! @return the model index of the loaded model (this is an index into
//!         the list of graphics models, not a molecule index).
int load_gltf_model(const std::string &gltf_file_name);

//! \brief set the model animation parameters
//!
//! Does nothing if \p model_index is out of range.
//!
//! @param model_index the model index (as returned by \c load_gltf_model())
//! @param amplitude the animation amplitude
//! @param wave_numer the animation wave number
//! @param freq the animation frequency
void set_model_animation_parameters(unsigned int model_index, float amplitude, float wave_numer, float freq);

//! \brief enable/disable the model animation (on or off)
//!
//! Does nothing if \p model_index is out of range.
//!
//! @param model_index the model index (as returned by \c load_gltf_model())
//! @param state true to animate, false to stop
void set_model_animation_state(unsigned int model_index, bool state);

//! \brief scale a (gltf) model
//!
//! Does nothing if \p model_index is out of range.
//!
//! @param model_index the model index (as returned by \c load_gltf_model())
//! @param scale_factor the scale factor
void scale_model(unsigned int model_index, float scale_factor);

//! \brief reset the frame buffers
//!
//! Recreate the frame buffers at the current size of the graphics window.
void reset_framebuffers();


/*  ----------------------------------------------------------------------- */
/*               Return Rotamer score (don't touch the model)               */
/*  ----------------------------------------------------------------------- */

//! \name Rotamer Scoring
//! \{

//! \brief Score rotamers for a residue (C++ interface)
//!
//! For each library rotamer of the residue type with probability above
//! \p lowest_probability, generate the rotamer and score it. The model is
//! not changed. The rotamers are in order of decreasing library probability.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param res_no the residue number
//! @param ins_code the insertion code
//! @param alt_conf the alternate conformation
//! @param imol_map the map molecule index (must be a valid map)
//! @param clash_flag 1 to calculate clash scores, 0 to skip (clash score is then -1)
//! @param lowest_probability the probability cut-off (in percent)
//! @return a vector of rotamer scores (each with name, rotamer probability
//!         (percent), clash score, per-atom densities and the summed side-chain
//!         density fit score). Empty if \p imol is not a valid model,
//!         \p imol_map is not a valid map or the residue is not found.
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
//! Each rotamer is scored for rotamer library probability, fit to density
//! and (optionally) clashes. The model is not changed.
//! The rotamers are in order of decreasing library probability.
//!
//! @param imol Model molecule number
//! @param chain_id Chain identifier
//! @param res_no Residue number
//! @param ins_code Insertion code
//! @param alt_conf Alternate conformation
//! @param imol_map Map for density scoring (must be a valid map, otherwise
//!        the empty list is returned)
//! @param clash_flag 1 to calculate clash scores, 0 to skip (clash score is then -1)
//! @param lowest_probability Minimum rotamer probability threshold (in percent)
//!
//! @return SCM - List of rotamer descriptions, each of the form
//!         (name probability atom-density-list density-fit clash-score),
//!         where probability is in percent, atom-density-list is a list of
//!         (atom-name density) pairs and density-fit is the sum of the
//!         densities at the side-chain atoms (excluding main-chain and CB).
//!         Empty list if the residue is not found, the map is not valid
//!         or there are no rotamers above the threshold.
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
//!             (list-ref rot 3))) ; density fit
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
//! Evaluates the library rotamer conformations (above the probability
//! threshold) and returns them with scores. The model is not changed.
//! This is the function to call before using auto_fit_best_rotamer.
//!
//! @param imol Model molecule number
//! @param chain_id Chain identifier
//! @param res_no Residue number
//! @param ins_code Insertion code (use "" if none)
//! @param alt_conf Alternate conformation (use "" for default)
//! @param imol_map Map molecule for density scoring (must be a valid map,
//!        otherwise an empty list is returned)
//! @param clash_flag 1 to calculate clash scores, 0 to skip (clash score is then -1)
//! @param lowest_probability Filter: only return rotamers above this probability (in percent)
//!
//! @return PyObject* - List of rotamers, each a list of the form
//!         [name, probability, density_fit, atom_densities, clash_score]
//!         where:
//!         - name: rotamer name (e.g. "mt")
//!         - probability: rotamer library probability (in percent)
//!         - density_fit: sum of the map densities at the side-chain atoms
//!           (excluding main-chain and CB atoms)
//!         - atom_densities: list of [atom_name, density] pairs
//!         - clash_score: clash score (-1 if clash_flag was 0)
//!
//!         Empty list if the residue is not found, the map is not valid
//!         or there are no rotamers above the threshold.
//!
//! \note Rotamers are ordered by decreasing library probability (not by fit)
//!
//! Example usage:
//! \code{.py}
//! # Score rotamers for LEU 42, considering density and clashes
//! rotamers = coot.score_rotamers_py(1, "A", 42, "", "", 2, 1, 0.01)
//!
//! print(f"Found {len(rotamers)} possible rotamers")
//! for rot in rotamers:
//!     name, prob, fit, atom_densities, clash = rot
//!     print(f"{name}: prob={prob:.1f}%, fit={fit:.3f}, clash={clash:.2f}")
//!
//! # pick the rotamer with the best fit to density
//! if rotamers:
//!     best = max(rotamers, key=lambda r: r[2])
//!     print(f"Best-fitting rotamer: {best[0]}")
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
//! \brief Cowtan's protein_db loops
//!
//! Search the ProteinDB fragment database (\c protein_db/protein.db in the
//! Coot package data directory) for main-chain fragments that fit the given
//! framework residues and the map. The residues given in \p residue_specs
//! (which must all be in the same chain) are the known flanking residues:
//! any residue number between the lowest and highest of the specs that is not
//! in the list (or not in the model) is treated as a gap to be filled.
//!
//! Several new molecules are created:
//!   - one molecule for each candidate loop ("Loop candidate #n"), undisplayed and inactive;
//!   - a consolidated molecule of all candidates ("All Loop candidates ", chain "Z",
//!     drawn in purple with thicker bonds);
//!   - a copy of the original residue range from \p imol_coords (made with
//!     protein_db_loop_specs_to_atom_selection_string()), undisplayed and inactive.
//!
//! @param imol_coords the model molecule index
//! @param residue_specs the framework residues (all in the same chain)
//! @param imol_map the map molecule index used to score the fragments
//! @param nfrags the number of candidate fragments to search for
//! @param preserve_residue_names if true, the candidate residues are named from the
//!        residue types of the matched database fragment; if false they are named "UNK"
//!
//! @return a pair: the first is (imol of the copy of the original loop region,
//!         imol of the consolidated candidates molecule), the second is the vector
//!         of molecule indices of the individual candidate loops.
//!         On failure (invalid model or map, or no candidates found) the first pair
//!         is (-1, -1) and the vector is empty.
std::pair<std::pair<int, int> , std::vector<int> >
protein_db_loops(int imol_coords,
                 const std::vector<coot::residue_spec_t> &residue_specs,
                 int imol_map, int nfrags, bool preserve_residue_names);
//! \brief convert loop residue specs to an mmdb atom selection string
//!
//! Used so that we can create an "original loop" molecule from the residue
//! specs picked: the selection extends over the range from the smallest residue
//! number to the largest (in the same chain), e.g. \c "//A/40-48".
//!
//! @param specs the residue specs (all must be in the same chain)
//! @return the atom selection string, or \c "////" if the specs are not all
//!         in a single chain (or empty)
std::string
protein_db_loop_specs_to_atom_selection_string(const std::vector<coot::residue_spec_t> &specs);
#ifdef USE_GUILE
//! \brief Cowtan's protein_db loops (Scheme interface)
//!
//! See protein_db_loops().
//!
//! @return a list: \c ((imol-original-loop imol-consolidated) (imol-candidate-0 imol-candidate-1 ...)),
//!         or \c \#f if \p residues_specs could not be converted to any residue specs
SCM protein_db_loops_scm(int imol_coords, SCM residues_specs, int imol_map, int nfrags, bool preserve_residue_names);
#endif

#ifdef USE_PYTHON
//! \brief Cowtan's protein_db loops
//!
//! return in the first pair, the imol of the new molecule generated
//! from an atom selection of the imol_coords for the residue selection
//! of the loop and the molecule number of the consolidated solutions
//! (displayed in purple).  and the second of the outer pair, there is
//! vector of molecule indices for each of the candidate loops.
//!
//! Use this to create hypotheses about where the atoms of the missing
//! residues could be. Often the top/first solution is the best one.
//! This fragment will then need to be patched back into molecule
//! imol_coords using copy_fragment().
//!
//! See protein_db_loops() for details.
//!
//! @param imol_coords the model molecule index
//! @param residues_specs a list of residue specs \c [chain_id, res_no, ins_code]
//!        for the framework residues (all in the same chain)
//! @param imol_map the map molecule index
//! @param nfrags the number of candidate fragments to search for
//! @param preserve_residue_names if True, use the residue types of the database
//!        fragments, otherwise name the residues "UNK"
//!
//! @return \c [[imol_original_loop, imol_consolidated], [imol_candidate_0, imol_candidate_1, ...]].
//!         On failure the first pair contains -1s (and the candidate list is empty);
//!         \c False is returned if \p residues_specs contains no valid residue specs.
PyObject *protein_db_loops_py(int imol_coords, PyObject *residues_specs, int imol_map, int nfrags, bool preserve_residue_names);
#endif

/* \} */


/* ------------------------------------------------------------------------- */
/*                      HOLE                                                 */
/* ------------------------------------------------------------------------- */
/*! \name Coot's Hole implementation */

//! \brief find the pore (channel) through a molecule between two points
//!
//! Run Coot's HOLE-like channel/pore analysis: a probe is passed along the
//! line from the start point to the end point and the maximum radius that fits
//! at each step is found. The pore surface is displayed as a generic display
//! object ("Probe surface"), a dialog with the probe radius as a function of
//! distance along the path is shown, and the radius profile is also written to
//! the file \c probe-radius.tab in the current directory. If a valid refinement
//! map is set (imol_refinement_map()) the map around the pore is also written
//! to \c hole.map.
//!
//! @param imol the model molecule index
//! @param start_x x coordinate of the start point (Å)
//! @param start_y y coordinate of the start point (Å)
//! @param start_z z coordinate of the start point (Å)
//! @param end_x x coordinate of the end point (Å)
//! @param end_y y coordinate of the end point (Å)
//! @param end_z z coordinate of the end point (Å)
//! @param colour_map_multiplier scale applied to the radius-to-colour mapping of the surface (the GUI uses 1.0)
//! @param colour_map_offset offset applied to the radius-to-colour mapping of the surface (the GUI uses 0.0)
//! @param n_runs currently unused
//! @param show_probe_radius_graph_flag currently unused (the probe radius text dialog is always shown)
//! @param export_surface_dots_file_name if not empty, write the surface dots to this file
//!        (one line per dot: x y z red green blue hex-colour)
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

//! \brief make a Gaussian surface for each chain of a model molecule
//!
//! The surface is separated into chains to make Generic Display Objects
//! (named "Gaussian Surface (Chain X)"). The surface parameters are taken from
//! the values set by set_gaussian_surface_sigma(), set_gaussian_surface_contour_level(),
//! set_gaussian_surface_box_radius(), set_gaussian_surface_grid_scale() and
//! set_gaussian_surface_fft_b_factor(); the colouring scheme is set by
//! set_gaussian_surface_chain_colour_mode(). Only the first model is used.
//!
//! @param imol the model molecule index
//! @return 0 (currently always)
int gaussian_surface(int imol);

//! \brief set the sigma for gaussian surface (default 4.4)
//!
//! @param s the sigma
void set_gaussian_surface_sigma(float s);

//! \brief set the contour_level for gaussian surface (default 4.0)
//!
//! @param s the contour level
void set_gaussian_surface_contour_level(float s);

//! \brief set the box_radius for gaussian surface (default 5.0)
//!
//! @param s the box radius
void set_gaussian_surface_box_radius(float s);

//! \brief set the grid_scale for gaussian surface (default 0.7)
//!
//! @param s the grid scale
void set_gaussian_surface_grid_scale(float s);
//! \brief set the fft B-factor for gaussian surface. Use 0 for no B-factor (default 100)
//!
//! @param f the B-factor
void set_gaussian_surface_fft_b_factor(float f);

//! \brief set the chain colour mode for Gaussian surfaces
//!
//! mode = 1 means each chain has its own colour (the default).
//! mode = 2 means the chain colour is determined from NCS/molecular symmetry (so
//!         that, in this mode, chains with the same sequence have the same colour).
//! (Any value other than 1 is currently treated as mode 2.)
//!
//! This affects surfaces made by subsequent calls to gaussian_surface().
//!
//! @param mode the colour mode
void set_gaussian_surface_chain_colour_mode(short int mode);

//! \brief set the opacity for a given molecule's gaussian_surface
//!
//! @param imol the model molecule index
//! @param opacity between 0.0 and 1.0 (default 1.0)
void set_gaussian_surface_opacity(int imol, float opacity);

//! \brief show the Gaussian surface overlay (GUI)
//!
//! Internal: shows the Gaussian surface frame, filled with the current surface
//! parameters and a model molecule chooser.
void show_gaussian_surface_overlay();

/* ------------------------------------------------------------------------- */
/*                      Cavities                                             */
/* ------------------------------------------------------------------------- */
/*! \name Coot's Cavities */
//! \brief find and display the cavities (internal pockets) of a model molecule
//!
//! The cavities are found using a grid of probe balls (probe radius 1.4 Å);
//! each non-trivial cavity is displayed as a semi-transparent Gaussian-surface
//! generic display object ("Cavity Surface #imol n"). As a side effect, the
//! files \c cavity-points.table and \c subpocket-points.table are written to the
//! current directory.
//!
//! @param imol the model molecule index
void show_cavities(int imol);


/* ------------------------------------------------------------------------- */
/*                      Acedrg for dictionary                                */
/* ------------------------------------------------------------------------- */
//! \brief make a dictionary for a residue using acedrg and the PDBe CCD entry
//!
//! The CCD mmCIF file for the residue's name is downloaded from the PDBe
//! (into the user's XDG data directory) and then acedrg is run on it in a
//! background thread (\c acedrg -r RES --noGeoOpt --coords -c file.cif). When
//! acedrg has finished, the output \c AcedrgOut.cif (in the current
//! directory) is read as a dictionary. Requires acedrg to be on the PATH.
//!
//! @param imol the model molecule index
//! @param spec the residue spec of the residue (its residue name is used as the comp-id)
void make_acedrg_dictionary_via_CCD_dictionary(int imol, const coot::residue_spec_t &spec);

/* ------------------------------------------------------------------------- */
/*                      LINKs                                                */
/* ------------------------------------------------------------------------- */

//! \brief make a link between the specified atoms
//!
//! A LINK record is added to the model containing the atoms (the atoms must be
//! in the same model). Any chem-mods defined in the dictionary for the resulting
//! link are applied (e.g. atom deletions). A backup is made first.
//!
//! @param imol the model molecule index
//! @param spec_1 the first atom
//! @param spec_2 the second atom
//! @param link_name currently unused
//! @param length currently unused
void
make_link(int imol, const coot::atom_spec_t &spec_1, const coot::atom_spec_t &spec_2,
          const std::string &link_name, float length);
#ifdef USE_GUILE
//! \brief make a link between the specified atoms (Scheme interface)
//!
//! See make_link(). \p spec_1 and \p spec_2 are Scheme atom specs;
//! \p link_name and \p length are currently unused.
void make_link_scm(int imol, SCM spec_1, SCM spec_2, const std::string&link_name, float length);
//! \brief return a list of the links in the given molecule (Scheme interface)
//!
//! @param imol the model molecule index
//! @return a list of \c (model-number atom-spec-1 atom-spec-2) items, one for each LINK.
//!         An empty list is returned for non-valid (i.e. non-model) molecules.
SCM link_info_scm(int imol);
#endif
#ifdef USE_PYTHON

//! \brief make a link between the specified atoms
//!
//! See make_link().
//!
//! @param imol the model molecule index
//! @param spec_1 the first atom spec, e.g. \c ["A", 42, "", " SG ", ""]
//! @param spec_2 the second atom spec
//! @param link_name currently unused
//! @param length currently unused
void make_link_py(int imol, PyObject *spec_1, PyObject *spec_2, const std::string&link_name, float length);

//! \brief return a list of the links in the given molecule.
//!
//! @param imol the model molecule index
//! @return a list of \c [model_number, atom_spec_1, atom_spec_2] items, one for each LINK.
//!         An empty list is returned for non-valid (i.e. non-model) molecules.
//!
PyObject *link_info_py(int imol);

//! \brief delete the links that contain the given residue
//!
//! Every LINK in which either of the linked atoms is in the given residue is deleted.
//!
//! @param imol the model molecule index
//! @param residue_spec_py the residue spec as 3 member list \c [chain_id, res_no, ins_code]
void delete_links_containing_residue_py(int imol, PyObject *residue_spec_py);

#endif

//! \brief show the acedrg link interface overlay (GUI)
//!
//! Internal: shows the acedrg link interface frame. A warning dialog is shown
//! if acedrg is not found on the PATH.
void show_acedrg_link_interface_overlay();

/* ------------------------------------------------------------------------- */
/*                      Drag and drop                                        */
/* ------------------------------------------------------------------------- */

/*! \name  Drag and Drop Functions */
// \{
//! \brief handle the string that get when a file or URL is dropped.
//!
//! An http:// or https:// URL is downloaded (into the "coot-download"
//! directory) and then read (as a dictionary, coordinates or an MTZ file), except
//! that a PDB-image .png URL is converted to an accession code and the entry is
//! fetched. A 4-character string is treated as a PDB accession code. Otherwise,
//! if the string is the name of an existing file it is read, or if it is a
//! \c file:/// URI with extension .cif, .pdb or .mtz, that file is read.
//!
//! @param uri the dropped string
//! @return 1 if the string was handled, 0 otherwise
int handle_drag_and_drop_string(const std::string &uri);
// \}


/* ------------------------------------------------------------------------- */
/*                      Map Display Control                                  */
/* ------------------------------------------------------------------------- */

/*! \name  Map Display Control */
// \{
//! \brief undisplay all maps except the given one
//!
//! The given map is displayed (if it was not already).
//!
//! @param imol_map the map molecule index
void undisplay_all_maps_except(int imol_map);
// \}


/* ------------------------------------------------------------------------- */
/*                      Map Contours                                         */
/* ------------------------------------------------------------------------- */

/*! \name Map Contouring Functions */

#ifdef USE_PYTHON
// \{
//! \brief return a list of pairs of vertices for the lines
//!
//! The map is contoured at \p contour_level in a box of the current map
//! radius around the current rotation centre.
//!
//! @param imol the map molecule index
//! @param contour_level the contour level (in map units, not rmsd)
//! @return a list of line segments, each \c [[x1,y1,z1],[x2,y2,z2]],
//!         or \c False if \p imol is not a valid map
PyObject *map_contours(int imol, float contour_level);

//! \brief return two lists: a list of vertices and a list of index-triples for connection
//!
//! The current triangle mesh of the map (as currently displayed) is returned.
//! Note that \p contour_level is currently not used - the map's current
//! contour level is used.
//!
//! @param imol the map molecule index
//! @param contour_level currently unused
//! @return \c [vertices, triangles] where vertices is a list of \c [x,y,z] and
//!         triangles a list of \c [i,j,k] vertex indices, or \c False if \p imol is not a valid map
PyObject *map_contours_as_triangles(int imol, float contour_level);

// \}
#endif // USE_PYTHON

//! \brief enable radial map colouring
//!
//! Colour the map by distance from a centre (see set_radial_map_colouring_centre()).
//! The map is recontoured if the state changes, so set the other radial
//! colouring parameters first.
//!
//! @param imol the map molecule index
//! @param state 1 to enable, 0 to disable
void set_radial_map_colouring_enabled(int imol, int state);

//! \brief radial map colouring centre
//!
//! The default is the centre of the unit cell.
//!
//! @param imol the map molecule index
//! @param x x coordinate of the centre (Å)
//! @param y y coordinate of the centre (Å)
//! @param z z coordinate of the centre (Å)
void set_radial_map_colouring_centre(int imol, float x, float y, float z);

//! \brief radial map colouring min
//!
//! Points closer than this to the centre get the colour at the start of the colour ramp.
//!
//! @param imol the map molecule index
//! @param r the minimum radius (Å)
void set_radial_map_colouring_min_radius(int imol, float r);

//! \brief radial map colouring max
//!
//! Points further than this from the centre get the colour at the end of the colour ramp.
//!
//! @param imol the map molecule index
//! @param r the maximum radius (Å)
void set_radial_map_colouring_max_radius(int imol, float r);

//! \brief radial map colouring inverted colour map
//!
//! @param imol the map molecule index
//! @param invert_state 1 to invert the colour ramp, 0 for normal
void set_radial_map_colouring_invert(int imol, int invert_state);

//! \brief radial map colouring saturation
//!
//! saturation is a number between 0 and 1, typically 0.5 (the default)
//!
//! @param imol the map molecule index
//! @param saturation the saturation
void set_radial_map_colouring_saturation(int imol, float saturation);



/* ------------------------------------------------------------------------- */
/*                      correlation maps                                     */
/* ------------------------------------------------------------------------- */

//! \name Map to Model Correlation
//! \{


//! \brief set the atom radius for map-to-model correlation functions
//!
//! The atom radius is not passed as a parameter to correlation
//! functions, so set it here (default is 1.5 Å).
//!
//! @param r the atom radius (Å)
void set_map_correlation_atom_radius(float r);

// Don't count the grid points of residues_specs that are in grid
// points of (potentially overlapping) neighbour_residue_spec.
//
#ifdef USE_GUILE
//! \brief map-to-model correlation (Scheme interface)
//!
//! See map_to_model_correlation() for the meaning of \p atom_mask_mode.
//! Grid points of \p residue_specs that are also covered by atoms of
//! \p neighb_residue_specs are not counted.
//!
//! @param imol the model molecule index
//! @param residue_specs Scheme list of residue specs
//! @param neighb_residue_specs Scheme list of neighbouring residue specs
//! @param atom_mask_mode controls which atoms are included
//! @param imol_map the map molecule index
//! @return the correlation coefficient as a real (nan on failure)
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
//! @param atom_mask_mode Controls which atoms to include (see map_to_model_correlation())
//! @param imol_map Map molecule number
//!
//! @return a list of 12 numbers: correlation, variance of the calculated (model)
//!         map, variance of the observed map, number of grid points, sum of
//!         calculated-map values, sum of observed-map values, Kolmogorov-Smirnov D
//!         of the sampled map values vs a normal with the whole-map mean and sd,
//!         Kolmogorov-Smirnov D vs a normal with the local mean and sd, map mean,
//!         local mean, map sd, local sd. (The map mean is -999 if \p imol_map is not a valid map.)
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
//! @param atom_mask_mode Controls which atoms to include (see map_to_model_correlation())
//! @param imol_map Map molecule number
//!
//! @return correlation coefficient as a Python float (nan on failure, e.g. invalid
//!         molecules or no grid points)
//!
//! \note Use atom_mask_mode=2 to evaluate side-chain fit specifically
//!
//! Example usage:
//! \code{.py}
//! # Evaluate side-chain fit for residues 40-44
//! residue_specs = [['A', res_no, ''] for res_no in range(40, 45)]
//! correlation = coot.map_to_model_correlation_py(1, residue_specs, [], 2, 2)
//! print(f"Side-chain correlation: {correlation}")
//! \endcode
PyObject *map_to_model_correlation_py(int imol,
                                      PyObject *residue_specs,
                                      PyObject *neighb_residue_specs,
                                      unsigned short int atom_mask_mode,
                                      int imol_map);

//! \brief Get map-to-model correlation statistics (Python interface)
//!
//! Returns the correlation and the underlying sums for the correlation
//! between the calculated (model) map and the observed map.
//!
//! @param imol Model molecule number
//! @param residue_specs Python list of residue specs
//! @param neighb_residue_specs Python list of neighboring residue specs
//! @param atom_mask_mode Controls which atoms to include (see map_to_model_correlation())
//! @param imol_map Map molecule number
//!
//! @return a list of 6 floats:
//!         \c [correlation, var_calc, var_map, n_points, sum_calc, sum_map]
//!         where "calc" is the map calculated from the model atoms and "map" the
//!         observed map. On failure n_points is 0 and the other values are 0 or nan.
//!
//! Example usage:
//! \code{.py}
//! stats = coot.map_to_model_correlation_stats_py(1, residues, [], 0, 2)
//! correlation, var_calc, var_map, n_points, sum_calc, sum_map = stats
//! print(f"Correlation: {correlation:.3f} from {int(n_points)} grid points")
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
//! It analyzes correlation statistics for all residues in a chain at once. For each
//! run of \p n_residue_per_residue_range consecutive residues, the correlation
//! between the observed map and a map calculated from the model is computed
//! (excluding grid points that are also covered by residues outside the run) and
//! reported against the middle residue of the run.
//!
//! @param imol Model molecule number
//! @param chain_id Chain identifier (e.g., "A", "B")
//! @param imol_map Map molecule number
//! @param n_residue_per_residue_range Number of residues per analysis window:
//!        - Use 1 for per-residue statistics (most common)
//!        - Use 3 for smoothed statistics over 3-residue windows
//! @param exclude_mainchain_NOC_flag Whether to split off the backbone:
//!        - 0: use all grid points; the first list has the all-atom statistics and
//!          the second list has empty statistics (n_points = 0)
//!        - 1: grid points within 1.8 Å of a backbone N, C, O (or H) atom go into the
//!          first list (i.e. it is then main-chain statistics) and the rest go into the
//!          second (side-chain) list
//!
//! @return a list of two lists \c [first_stats, sidechain_stats]. Each is a list of items
//!         \c [residue_spec, [n_points, correlation]] where residue_spec is
//!         \c [chain_id, res_no, ins_code]. Both lists are empty if \p imol or \p imol_map is not valid.
//!
//! If the residue does not have a side-chain then the number of grid points may be 0 and the
//! correlation is then nan.
//!
//! \note This function analyzes the ENTIRE chain at once, making it very efficient
//! \note Use \p exclude_mainchain_NOC_flag = 1 and the second list to identify problem side chains specifically
//!
//! Example usage - Find worst-fitting side chain:
//! \code{.py}
//! # Get correlation statistics for all residues in chain A
//! # (args: imol, chain_id, imol_map, n_residue_per_residue_range, exclude_mainchain_NOC_flag)
//! mainchain_stats, sidechain_stats = coot.map_to_model_correlation_stats_per_residue_range_py(1, "A", 2, 1, 1)
//!
//! # Find worst-fitting side chain (ignoring residues with no side-chain grid points)
//! scored = [item for item in sidechain_stats if item[1][0] > 0]
//! worst_residue = min(scored, key=lambda item: item[1][1])
//!
//! chain_id, resno, ins_code = worst_residue[0]
//! correlation = worst_residue[1][1]
//!
//! print(f"Worst side chain: {chain_id} {resno}, correlation = {correlation:.3f}")
//!
//! # Center on worst residue
//! coot.set_go_to_atom_chain_residue_atom_name(chain_id, resno, 'CA')
//! \endcode
//!
//! Example usage - Compare main-chain vs side-chain fit:
//! \code{.py}
//! mainchain, sidechain = coot.map_to_model_correlation_stats_per_residue_range_py(1, "A", 2, 1, 1)
//!
//! side_corr = {tuple(spec): stats[1] for spec, stats in sidechain}
//! for spec, stats in mainchain:
//!     mc_corr = stats[1]
//!     sc_corr = side_corr.get(tuple(spec))
//!     if sc_corr is not None and mc_corr > 0.7 and sc_corr < 0.5:
//!         print(f"Residue {spec}: Good backbone, poor sidechain")
//! \endcode
PyObject *
map_to_model_correlation_stats_per_residue_range_py(int imol,
                                                    const std::string &chain_id,
                                                    int imol_map,
                                                    unsigned int n_residue_per_residue_range,
                                                    short int exclude_mainchain_NOC_flag);

#endif

//! \}

// Map to Model Correlation Functions - Enhanced Doxygen Documentation
//
// These functions assess how well a molecular model fits into electron density maps.
// Essential for model validation and identifying poorly-fitted regions.

//! \name Map to Model Correlation
//! \{

//! \brief Calculate the correlation between a map and model for specific residues
//!
//! This function calculates the correlation between the observed map and a map
//! calculated from the selected atoms of the specified residues, over the grid
//! points within the atom radius (see set_map_correlation_atom_radius()),
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
//!        - 10: All atoms, with the atom radius dependent on the atom's B-factor
//! @param imol_map Map molecule number to correlate against
//!
//! @return Correlation coefficient (float) between model and map
//!         (nan on failure, e.g. invalid molecules or no grid points)
//!
//! \note In the current implementation, modes 4 and 5 select all atoms (rather
//!       than main-chain or side-chain atoms) for standard amino acids.
//! \note Use this after refinement to evaluate if the fit improved
//!
//! Example:
//! \code{.cpp}
//! std::vector<coot::residue_spec_t> residues = {{chain_id, 40, ""}, {chain_id, 41, ""}};
//! std::vector<coot::residue_spec_t> neighbors;
//! float corr = map_to_model_correlation(1, residues, neighbors, 2, 2);
//! \endcode
float
map_to_model_correlation(int imol,
                         const std::vector<coot::residue_spec_t> &residue_specs,
                         const std::vector<coot::residue_spec_t> &neigh_residue_specs,
                         unsigned short int atom_mask_mode,
                         int imol_map);

//! \brief Get the map-to-model correlation statistics for a set of residues
//!
//! A model map is calculated from the selected atoms of \p residue_specs and
//! compared with the map \p imol_map over the grid points that lie within the
//! map-correlation atom radius of those atoms (1.5 Å by default, see
//! set_map_correlation_atom_radius()). Grid points within that radius of the atoms
//! of \p neigh_residue_specs are excluded, so that density belonging to neighbouring
//! residues does not contribute.
//!
//! @param imol the model molecule index
//! @param residue_specs the residues to evaluate
//! @param neigh_residue_specs neighbouring residues whose atoms mask out grid points
//! @param atom_mask_mode which atoms of each residue are used:
//!        - 0: all atoms
//!        - 1: main-chain atoms if a standard amino acid, else all atoms
//!        - 2: side-chain atoms if a standard amino acid, else all atoms
//!        - 3: side-chain atoms excluding CB if a standard amino acid, else all atoms
//!        - 4: main-chain atoms if a standard amino acid, else nothing
//!        - 5: side-chain atoms if a standard amino acid, else nothing
//!        - 10: all atoms with a B-factor-dependent radius
//! @param imol_map the map molecule index
//!
//! @return a coot::util::density_correlation_stats_info_t: the accumulated sums
//!         (\c n, \c sum_x, \c sum_y, \c sum_xy, \c sum_sqrd_x, \c sum_sqrd_y,
//!         where x is the model-calculated map and y is \p imol_map) and the sampled
//!         map values (\c density_values). Use its \c correlation(), \c var_x() and
//!         \c var_y() member functions. If \p imol is not a valid model or
//!         \p imol_map is not a valid map (or no atoms were selected) the returned
//!         object is empty (\c n is 0).
//!
//! Example:
//! \code{.cpp}
//! auto stats = map_to_model_correlation_stats(1, residues, neighbors, 0, 2);
//! std::cout << "Correlation: " << stats.correlation() << " from " << stats.n << " points" << std::endl;
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
//! Returns individual correlation values (model-calculated map vs \p imol_map) for
//! each residue, using grid points within 1.5 Å of the selected atoms. Grid points
//! that are within that radius of atoms of more than one residue are excluded.
//! Useful for identifying which specific residues fit poorly.
//!
//! @param imol the model molecule index
//! @param specs the residues to evaluate. Note: with \p atom_mask_mode 0 the
//!        atom selection is the whole molecule, so a correlation is returned for
//!        every residue in the molecule, not only for those in \p specs.
//! @param atom_mask_mode which atoms to include (0: all atoms, 1: main-chain,
//!        2: side-chain, 3: side-chain excluding CB - for non-standard residues these
//!        use all atoms; 4: main-chain, 5: side-chain - for non-standard residues these
//!        use no atoms)
//! @param imol_map the map molecule index
//!
//! @return a vector of (residue_spec, correlation) pairs. Residues with fewer than 2
//!         contributing grid points are omitted. Empty if \p imol or \p imol_map is not
//!         valid.
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

//! \brief Get map density statistics per residue
//!
//! Despite the name, this does not calculate a correlation: for each residue it
//! accumulates the values of \p imol_map at the grid points within
//! \p atom_radius_for_masking of the residue's selected atoms. Grid points that are
//! within that radius of atoms of more than one residue are excluded.
//!
//! @param imol the model molecule index
//! @param residue_specs the residues to evaluate
//! @param atom_mask_mode which atoms to include (as for map_to_model_correlation_stats())
//! @param atom_radius_for_masking the masking radius around each atom (in Å, typically 1.5)
//! @param imol_map the map molecule index
//!
//! @return a map of residue_spec to coot::util::density_stats_info_t (members \c n,
//!         \c sum, \c sum_sq, \c sum_weight; use \c mean_and_variance() to get the
//!         mean and variance of the density). Empty if \p imol or \p imol_map is not
//!         valid.
//!
//! Example:
//! \code{.cpp}
//! auto stats_map = map_to_model_correlation_stats_per_residue(1, specs, 0, 1.5, 2);
//! for (const auto &pair : stats_map) {
//!     std::pair<double, double> mv = pair.second.mean_and_variance();
//!     std::cout << pair.first << ": mean=" << mv.first << " var=" << mv.second << std::endl;
//! }
//! \endcode
std::map<coot::residue_spec_t, coot::util::density_stats_info_t>
map_to_model_correlation_stats_per_residue(int imol,
                                           const std::vector<coot::residue_spec_t> &residue_specs,
                                           unsigned short int atom_mask_mode,
                                           float atom_radius_for_masking,
                                           int imol_map);

//! \brief Get map-to-model correlation statistics for runs of residues along a chain
//!
//! The chain is split into overlapping windows of \p n_residue_per_residue_range
//! residues; the correlation of the model-calculated map with \p imol_map over each
//! window is reported against the middle residue of that window. Windows that
//! include waters, nucleotides or HETATM residues are skipped. The atom mask radius
//! is 2.8 Å (and 1.8 Å around the main-chain N, O and C atoms when they are excluded).
//!
//! @param imol the model molecule index
//! @param chain_id the chain id (e.g. "A")
//! @param imol_map the map molecule index
//! @param n_residue_per_residue_range the number of residues per window (e.g. 1, 3 or 5)
//! @param exclude_NOC_flag if non-zero, also calculate side-chain statistics (excluding
//!        grid points near the main-chain N, O and C atoms)
//!
//! @return a pair of maps from (middle) residue spec to
//!         coot::util::density_correlation_stats_info_t:
//!         - first: all-atom statistics
//!         - second: side-chain statistics (only filled with data when
//!           \p exclude_NOC_flag is set, otherwise the entries have \c n = 0)
//!
//!         Both maps are empty if \p imol or \p imol_map is not valid.
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
//!         return a.second.correlation() < b.second.correlation();
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
//! See map_to_model_correlation_per_residue().
//!
//! @param imol the model molecule index
//! @param residue_specs a list of residue specs
//! @param atom_mask_mode which atoms to include (see map_to_model_correlation_per_residue())
//! @param imol_map the map molecule index
//!
//! @return a list of (residue-spec correlation) lists; empty on failure
SCM
map_to_model_correlation_per_residue_scm(int imol, SCM residue_specs,
                                         unsigned short int atom_mask_mode,
                                         int imol_map);

//! \brief Map density statistics per residue (Guile interface)
//!
//! See map_to_model_correlation_stats_per_residue().
//!
//! @param imol the model molecule index
//! @param residue_specs_scm a list of residue specs
//! @param atom_mask_mode which atoms to include
//! @param atom_radius_for_masking the masking radius around each atom (in Å)
//! @param imol_map the map molecule index
//!
//! @return a list of (residue-spec (mean variance)) items, where mean and variance are
//!         of the map density around the residue; empty on failure
SCM
map_to_model_correlation_stats_per_residue_scm(int imol,
                                               SCM residue_specs_scm,
                                               unsigned short int atom_mask_mode,
                                               float atom_radius_for_masking,
                                               int imol_map);

//! \brief Map-to-model correlation stats for runs of residues along a chain (Guile interface)
//!
//! See map_to_model_correlation_stats_per_residue_range() for the parameters.
//!
//! @return a list of two lists (all-atom, side-chain) whose items are
//!         (residue-spec (n correlation)), in residue-spec order.
SCM map_to_model_correlation_stats_per_residue_range_scm(int imol, const std::string &chain_id, int imol_map,
                                                         unsigned int n_residue_per_residue_range,
                                                         short int exclude_NOC_flag);


//! \brief QQ plot of the model density correlation, reported per residue
//!
//! Shows a "Difference Map QQ Plot" dialog (map quantiles vs reference normal
//! quantiles) for the density around the given residues. This only does anything if
//! Coot was compiled with goocanvas.
//!
//! @param imol the model molecule index
//! @param residue_specs_scm a list of residue specs
//! @param neigh_residue_specs_scm a list of neighbouring residue specs used for masking
//! @param atom_mask_mode which atoms to include:
//!        - 0: all-atoms
//!        - 1: main-chain atoms if is standard amino-acid, else all atoms
//!        - 2: side-chain atoms if is standard amino-acid, else all atoms
//!        - 3: side-chain atoms-excluding CB if is standard amino-acid, else all atoms
//!        - 4: main-chain atoms if is standard amino-acid, else nothing
//!        - 5: side-chain atoms if is standard amino-acid, else nothing
//! @param imol_map the map molecule index
//!
//! @return always \#f
SCM qq_plot_map_and_model_scm(int imol,
                              SCM residue_specs_scm,
                              SCM neigh_residue_specs_scm,
                              unsigned short int atom_mask_mode,
                              int imol_map);
#endif

#ifdef USE_PYTHON
//! \brief Get map-to-model correlation per residue (Python interface)
//!
//! See map_to_model_correlation_per_residue().
//!
//! @param imol the model molecule index
//! @param residue_specs a list of residue specs, e.g. [["A", 42, ""], ...]
//! @param atom_mask_mode which atoms to include (see map_to_model_correlation_per_residue())
//! @param imol_map the map molecule index
//!
//! @return a list of [residue_spec, correlation] items; an empty list on failure
//!
//! Example usage:
//! \code{.py}
//! residues = coot.get_residues_in_chain_py(1, "A")
//! correlations = coot.map_to_model_correlation_per_residue_py(1, residues, 0, 2)
//! # Find worst 10 residues
//! worst_10 = sorted(correlations, key=lambda x: x[1])[:10]
//! for spec, corr in worst_10:
//!     print(f"Residue {spec}: correlation = {corr:.3f}")
//! \endcode
PyObject *map_to_model_correlation_per_residue_py(int imol, PyObject *residue_specs,
                                                  unsigned short int atom_mask_mode,
                                                  int imol_map);

//! \brief QQ plot of the model density (Python interface)
//!
//! See qq_plot_map_and_model_scm(). This only does anything if Coot was compiled
//! with goocanvas.
//!
//! @return always False
PyObject *qq_plot_map_and_model_py(int imol,
                              PyObject *residue_specs_py,
                              PyObject *neigh_residue_specs_py,
                              unsigned short int atom_mask_mode,
                              int imol_map);
#endif

#ifdef __cplusplus
#ifdef USE_GUILE
//! \brief Simple density score for a residue (Guile interface)
//!
//! See density_score_residue().
//!
//! @param imol the model molecule index
//! @param residue_spec a residue spec, e.g. '("A" 42 "")
//! @param imol_map the map molecule index
//!
//! @return the occupancy-weighted sum of the map values at the atom positions
//!         (0.0 on failure)
//!
//! Example usage:
//! \code{.scm}
//! (define score (density-score-residue-scm 1 '("A" 42 "") 2))
//! (format #t "Density score: ~a~%" score)
//! \endcode
float density_score_residue_scm(int imol, SCM residue_spec, int imol_map);
#endif
#ifdef USE_PYTHON

//! \brief Simple density score for a residue (Python interface)
//!
//! See density_score_residue().
//!
//! @param imol the model molecule index
//! @param residue_spec a residue spec, e.g. ["A", 42, ""]
//! @param imol_map the map molecule index
//!
//! @return the occupancy-weighted sum of the map values at the atom positions
//!         (0.0 on failure)
//!
//! \note For density-fit validation, map_to_model_correlation_stats_per_residue_range_py()
//!       is usually more informative.
//!
//! Example usage:
//! \code{.py}
//! score = coot.density_score_residue_py(1, ["A", 42, ""], 2)
//! print(f"Density score: {score:.3f}")
//! \endcode
float density_score_residue_py(int imol, PyObject *residue_spec, int imol_map);
#endif
#endif

//! \brief Simple density score for given residue (C++ interface)
//!
//! The score is the sum over the residue's atoms of (map value at the atom
//! position × atom occupancy). It is in map units and is not normalised, so it
//! scales with the number of atoms and with the map scale - it is not a correlation.
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//! @param res_no the residue number
//! @param ins_code the insertion code (use "" if none)
//! @param imol_map the map molecule index
//!
//! @return the score, or 0.0 if the molecules are not valid or the residue is not found
//!
//! Example:
//! \code{.cpp}
//! float score = density_score_residue(1, "A", 42, "", 2);
//! \endcode
float density_score_residue(int imol, const char *chain_id, int res_no, const char *ins_code, int imol_map);


#ifdef USE_GUILE
//! \brief Get the mean value of a map (Guile interface)
//!
//! Returns the (cached) mean density value of the map.
//!
//! @param imol the map molecule index
//!
//! @return the mean as a number, or \#f if \p imol is not a valid map
//!
//! Example usage:
//! \code{.scm}
//! (define mean (map-mean-scm 2))
//! (format #t "Map mean: ~a~%" mean)
//! \endcode
SCM map_mean_scm(int imol);
//! \brief Get the standard deviation (sigma, rmsd) of a map (Guile interface)
//!
//! Returns the (cached) standard deviation of the density values in the map.
//! This is the "sigma" used for contouring at "N sigma" levels.
//!
//! @param imol the map molecule index
//!
//! @return the standard deviation as a number, or \#f if \p imol is not a valid map
//!
//! \note Contouring at "1.5 sigma" sets the absolute contour level to 1.5 times this value.
//!
//! Example usage:
//! \code{.scm}
//! (define sigma (map-sigma-scm 2))
//! (format #t "Contour at 1.5 sigma = ~a~%" (* 1.5 sigma))
//! \endcode
SCM map_sigma_scm(int imol);

//! \brief Get map statistics (Guile interface)
//!
//! Calculates statistics of all the (non-NaN) grid values of the map.
//!
//! @param imol the map molecule index
//!
//! @return a list (mean standard-deviation skew kurtosis) or \#f if \p imol is not a
//!         valid map. The skew is the (unnormalised) third central moment; the
//!         kurtosis is the fourth central moment divided by the variance squared (not
//!         the excess kurtosis, so it is about 3 for a normal distribution).
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
//! Returns the (cached) mean density value of the map.
//!
//! @param imol the map molecule index
//!
//! @return the mean as a float, or False if \p imol is not a valid map
//!
//! Example usage:
//! \code{.py}
//! mean = coot.map_mean_py(2)
//! if mean is not False:
//!     print(f"Map mean: {mean}")
//! \endcode
PyObject *map_mean_py(int imol);
//! \brief Get the standard deviation (sigma, rmsd) of a map (Python interface)
//!
//! Returns the (cached) standard deviation of the density values. This is the
//! "sigma" used when contouring at "N sigma": the absolute contour level is N times
//! this value.
//!
//! @param imol the map molecule index
//!
//! @return the sigma as a float, or False if \p imol is not a valid map
//!
//! Example usage:
//! \code{.py}
//! sigma = coot.map_sigma_py(2)
//! if sigma is not False:
//!     print(f"1.5 sigma contour level: {1.5 * sigma}")
//! \endcode
PyObject *map_sigma_py(int imol);
//! \brief Get map statistics (Python interface)
//!
//! Calculates statistics of all the (non-NaN) grid values of the map.
//!
//! @param imol the map molecule index
//!
//! @return a list [mean, std_dev, skew, kurtosis] or False if \p imol is not a valid map
//!         - mean: average density value
//!         - std_dev: standard deviation (sigma)
//!         - skew: the third central moment (not normalised by sigma cubed, so it
//!           depends on the map scale)
//!         - kurtosis: the fourth central moment divided by the variance squared
//!           (not excess kurtosis: about 3 for a normal distribution)
//!
//! Example usage:
//! \code{.py}
//! stats = coot.map_statistics_py(2)
//! if stats is not False:
//!     mean, sigma, skew, kurtosis = stats
//!     print(f"Mean: {mean:.3f} Sigma: {sigma:.3f} Skew: {skew:.3f} Kurtosis: {kurtosis:.3f}")
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
//! \brief get the sequence for chain_id in imol as FASTA
//!
//! @param imol the model molecule index
//! @param chain_id the chain id
//!
//! @return a FASTA-format string: a "> <molecule-name> <chain-id>" header line,
//!         a blank line, then the sequence
std::string get_sequence_as_fasta_for_chain(int imol, const std::string &chain_id);

//! \brief write the sequence for imol as fasta
//!
//! Writes the FASTA sequences of all the chains of the molecule to a file.
//!
//! @param imol the model molecule index
//! @param file_name the output file name
void write_sequence(int imol, const std::string &file_name);

//! \brief trace the given map and try to apply the sequence in
//! the given pir file
//!
//! Experimental. Creates a new (initially empty) model molecule and a copy of the map,
//! then runs the tracer in a background thread; the new model is updated in the
//! graphics as the trace progresses.
//!
//! @param imol_map the map molecule index
//! @param pir_file_name the file name of the sequence file (read as multi-FASTA/PIR)
void res_tracer(int imol_map, const std::string &pir_file_name);


/* \} */


/* ------------------------------------------------------------------------- */
/*                      interesting positions list                           */
/* ------------------------------------------------------------------------- */
#ifdef USE_GUILE
//! \brief register a list of user-defined interesting positions (Guile interface)
//!
//! @param pos_list a list of ((x y z) label) items. Malformed items are ignored.
//!        This replaces any previously registered list.
void register_interesting_positions_list_scm(SCM pos_list);
#endif // USE_GUILE
#ifdef USE_PYTHON
//! \brief register a list of user-defined interesting positions (Python interface)
//!
//! @param pos_list a list of [[x, y, z], label] items, where x, y and z must be
//!        floats and label is a string. Malformed items are ignored.
//!        This replaces any previously registered list.
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
//!        (potentially thousands) of atom overlaps.
//! @return a list of dictionaries with contact information (keys "atom-1-spec",
//!        "atom-2-spec", "overlap-volume", "radius-1", "radius-2").
//!        The list is sorted by largest overlap first.
//!        Return False on failure (invalid molecule or \p n_max_pairs < -1).
//!
PyObject *molecule_atom_overlaps_py(int imol, int n_max_pairs);
#endif // USE_PYTHON
#ifdef USE_GUILE
//! \brief get all the atom overlaps for the molecule (Guile interface)
//!
//! @param imol the molecule index
//! @return a list of (atom-spec-1 atom-spec-2 radius-1 radius-2 overlap-volume) items,
//!         sorted by largest overlap first; \#f if \p imol is not a valid model.
//!         If there are no overlaps because a dictionary was missing, a warning
//!         string is returned instead.
SCM molecule_atom_overlaps_scm(int imol);
#endif // USE_GUILE

/* ------------------------------------------------------------------------- */
/*                       Alignment functions (now C++)                       */
/* ------------------------------------------------------------------------- */

//! \brief align sequence to closest chain (compare across all chains
//!   in all molecules).
//!
//! Typically match_fraction is 0.95 or so. If a matching chain is found, the
//! sequence is assigned to that chain (as by assign_sequence_from_string()).
//!
//! Return the molecule number and chain id if successful, return -1 as the
//! molecule number if not.
//!
std::pair<int, std::string>
align_to_closest_chain(std::string target_seq, float match_fraction);

#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_PYTHON
//! \brief align sequence to closest chain (Python interface)
//!
//! See align_to_closest_chain().
//!
//! @return a list [imol, chain_id] on success, False on failure
PyObject *align_to_closest_chain_py(std::string target_seq, float match_fraction);
#endif /* USE_PYTHON */
#ifdef USE_GUILE
//! \brief align sequence to closest chain (Guile interface)
//!
//! See align_to_closest_chain().
//!
//! @return a list (imol chain-id) on success, \#f on failure
SCM align_to_closest_chain_scm(std::string target_seq, float match_fraction);
#endif /* USE_GUILE */
#endif /* c++ */


#ifdef __cplusplus/* protection from use in callbacks.c, else compilation probs */
#ifdef USE_GUILE
//! \brief experimental spherical density overlap test - internal
//!
//! @param i_scm the model molecule index
//! @param j_scm the map molecule index
//!
//! Writes a map fragment around the origin to "map-fragment-at-origin.map" and runs
//! a simple EMMA-style overlap test, printing the results.
//!
//! @return always \#f
SCM spherical_density_overlap(SCM i_scm, SCM j_scm);
#endif // USE_GUILE
#endif // __cplusplus

//! \brief resolve clashing side chains by deleting side chains
//!
//! For each pair of residues with an atom overlap volume greater than 2.0 Å³,
//! the side chain of the larger residue type is deleted. Waters are ignored.
//!
//! @param imol the model molecule index
void resolve_clashing_sidechains_by_deletion(int imol);

//! \brief resolve clashing side chains by rebuilding them
//!
//! For each pair of residues with an atom overlap volume greater than 2.0 Å³,
//! both side chains are deleted and rebuilt (filled as partial residues, using the
//! refinement map if one is set). Waters are ignored.
//!
//! @param imol the model molecule index
void resolve_clashing_sidechains_by_rebuilding(int imol);

/*  ----------------------------------------------------------------------- */
/*                  GUIL Utility Functions                                  */
/*  ----------------------------------------------------------------------- */
//! \brief make a simple text dialog.
//!
//! Shows a dialog with a read-only, word-wrapped text view and a Close button
//! (does nothing if the graphics interface is not in use).
//!
//! @param dialog_title the window title
//! @param text the text to display
//! @param geom_x the default window width (pixels)
//! @param geom_y the default window height (pixels)
void simple_text_dialog(const std::string &dialog_title, const std::string &text,
                        int geom_x, int geom_y);


/*  ----------------------------------------------------------------------- */
/*                  Phenix Functions                                        */
/*  ----------------------------------------------------------------------- */

//! \brief phenix GEO bonds representation
//!
//! Redraws the bonds of the molecule using the bonds of the given phenix geometry.
//! This function is not for scripting
void graphics_to_phenix_geo_representation(int imol, int mode,
                                           const coot::phenix_geo::phenix_geometry &g);

//! \brief phenix GEO bonds representation, read GEO info from file
//!
//! @param imol the molecule index
//! @param mode currently unused, so use 0
//! @param geo_file_name is the file name of the phenix_geo file
void graphics_to_phenix_geo_representation(int imol, int mode,
                                           const std::string &geo_file_name);

//! \brief validate using phenix geo bonds
//!
//! Typically this would be called shortly after
//! graphics_to_phenix_geo_representation(). Shows buttons in the validation pane
//! for the geometry restraints whose residual is greater than 4.4; clicking a button
//! centres on that restraint.
//!
//! @param imol the molecule index
//! @param geo_file_name is the file name of the phenix_geo file
void validate_using_phenix_geo_bonds(int imol, const std::string &geo_file_name);

/*  ----------------------------------------------------------------------- */
/*                  Client/Server                                        */
/*  ----------------------------------------------------------------------- */
#ifdef USE_PYTHON
//! \brief set a python command string to be run on each graphics draw
//!
//! Note: the draw code that ran this string is currently disabled, so setting it
//! has no effect.
void set_python_draw_function(const std::string &command_string);
#endif // USE_PYTHON


/*  ----------------------------------------------------------------------- */
/*                  Pathology Plots                                         */
/*  ----------------------------------------------------------------------- */
#ifdef USE_PYTHON
//! \brief get data for the data-pathology plots from an MTZ file
//!
//! @param mtz_file_name the MTZ file name
//! @param fp_col the F column label
//! @param sigfp_col the sigF column label
//!
//! @return a list [invresolsq_max, fp_vs_reso, fosf_vs_reso, sigf_vs_f, fosf_vs_f]
//!         where the last four are lists of (x, y) tuples: (1/d², F), (1/d², F/sigF),
//!         (F, sigF) and (F, F/sigF). If there are more than 20000 reflections the data
//!         are subsampled. Returns False on failure.
PyObject *pathology_data(const std::string &mtz_file_name,
                         const std::string &fp_col,
                         const std::string &sigfp_col);
#endif // USE_PYTHON

/*  ----------------------------------------------------------------------- */
/*                  Utility Functions                                       */
/*  ----------------------------------------------------------------------- */

//! \brief encoding of ints
//!
//! These functions are for storing the molecule number and (some other
//! number) as an int and used with GPOINTER_TO_INT and GINT_TO_POINTER.
//! The encoding is 1000 * i1 + i2, so i2 must be in the range 0-999.
int encode_ints(int i1, int i2);
//! \brief decode an int made by encode_ints()
//!
//! @return the pair (i1, i2)
std::pair<int, int> decode_ints(int i);

//! \brief store username and password for the database.
//!
//! @param key the key under which the user name and password are stored
//! @param user_name the user name
//! @param passwd the password
void store_keyed_user_name(std::string key, std::string user_name, std::string passwd);

#ifdef USE_GUILE
//! \brief convert a list of residue specs to a vector of residue specs - internal
//!
//! Items that are not valid residue specs are skipped.
std::vector<coot::residue_spec_t> scm_to_residue_specs(SCM s);
// and that backwards is scm_residue() above.
//! \brief convert a key name (e.g. "Return") to a GDK key-sym code - internal
//!
//! @return the key-sym code, or -1 if the name is not found or \p s_scm is not a string
int key_sym_code_scm(SCM s_scm);
#endif // USE_GUILE
#ifdef USE_PYTHON
//! \brief convert a list of residue specs to a vector of residue specs - internal
//!
//! Each spec may be [chain_id, res_no, ins_code] or a 4-item list with a leading
//! item that is skipped. Items that are not 3 or 4 long are skipped.
std::vector<coot::residue_spec_t> py_to_residue_specs(PyObject *s);
//! \brief convert a key name (e.g. "Return") to a GDK key-sym code - internal
//!
//! @return the key-sym code, or -1 if the name is not found or \p po is not a string
int key_sym_code_py(PyObject *po);
#endif // USE_PYTHON
#ifdef USE_GUILE
#ifdef USE_PYTHON
//! \brief convert a scheme object to a python object - internal
//!
//! Handles lists (recursively), booleans, integers, reals and strings; anything
//! else becomes None.
PyObject *scm_to_py(SCM s);
//! \brief convert a python object to a scheme object - internal
//!
//! Handles lists (recursively), booleans, integers, floats, strings and None
//! (which becomes unspecified); anything else becomes \#f.
SCM py_to_scm(PyObject *o);
#endif // USE_GUILE
#endif // USE_PYTHON

#ifdef USE_GUILE
//! \brief make a space group from a list of symmetry operator strings - internal
//!
//! @return the space group; it is null (\c is_null()) on failure
clipper::Spacegroup scm_symop_strings_to_space_group(SCM symop_string_list);
#endif

#ifdef USE_PYTHON
//! \brief make a space group from a list of symmetry operator strings - internal
//!
//! @return the space group; it is null (\c is_null()) on failure
clipper::Spacegroup py_symop_strings_to_space_group(PyObject *symop_string_list);
#endif

//! \brief enable or disable sounds (coot needs to have been compiled with sounds of course)
void set_use_sounds(bool state);

//! \brief turn off sounds and particles and textures
//!
//! Also turns off the happy-face residue markers and the "unhappy atom"
//! (bad non-bonded contact and chiral volume outlier) markers.
void curmudgeon_mode();

//! \brief easter egg 2023
//!
//! Adds a pumpkin graphics object.
void halloween();

//! \brief display an SVG file in a dialog
//!
//! @param file_name the SVG file name (nothing happens if the file does not exist)
void display_svg_from_file_in_a_dialog(const std::string &file_name);

//! \brief display an SVG string in a dialog
//!
//! Requires Coot to have been compiled with librsvg (otherwise does nothing).
//!
//! @param string the SVG document as a string
//! @param title the dialog title (prefixed by "Coot: ")
void display_svg_from_string_in_a_dialog(const std::string &string, const std::string &title);

//! \brief display a PAE (predicted aligned error) plot from a JSON file in a dialog
//!
//! Note: \p imol is currently not used - the molecule of the active atom is used
//! instead, and nothing is shown if there is no active atom.
//!
//! @param imol the model molecule index
//! @param file_name the PAE JSON file name (e.g. as downloaded from the AlphaFold DB)
void display_pae_from_file_in_a_dialog(int imol, const std::string &file_name);

//! \brief read an "interesting places" JSON file and show it as a dialog of buttons
//!
//! See read_interesting_places_json() for the format.
//!
//! @param file_name the JSON file name
void read_interesting_places_json_file(const std::string &file_name);

//! \brief read "interesting places" JSON and show it as a dialog of buttons
//!
//! The JSON is an object with a "title" and a list of "sections", each of which has
//! a "title" and a list of "items". Each item has a "position-type" (one of
//! "by-atom-spec", "by-atom-spec-pair", "by-residue-spec" or "by-coordinates"), a
//! "label", and respectively an "atom-spec", "atom-1-spec" and "atom-2-spec",
//! a "residue-spec" or a "position"; "badness" is optional. The buttons are attached
//! to the molecule of the active atom, so a model must be displayed.
//!
//! @param json_as_string the JSON as a string
void read_interesting_places_json(const std::string &json_as_string);

//! \brief set up the tomogram section slider for the given map
//!
//! Hides the main toolbar box and shows the slider, with a range of 0 to the number of
//! sections (along the third grid axis) minus 1, set to the middle section.
//!
//! @param imol the map molecule index
//! @return the section index (the middle section currently), or -1 if \p imol is not
//!         a valid map
int setup_tomo_slider(int imol);
//! \brief show a tomogram map as a section view
//!
//! Sets the zoom, clipping and rotation centre (the centre of the cell) and displays
//! the given section of the map.
//!
//! @param imol the map molecule index
//! @param axis_id the section index (despite the name, the section axis is currently
//!        always the third axis)
void tomo_section_view(int imol, int axis_id);
//! \brief set the tomogram section by moving the section slider
//!
//! @param imol the map molecule index (currently unused - the slider's molecule is used)
//! @param section_index the section index
void set_tomo_section_view_section(int imol, int section_index);

//! \brief set tomo picker is active
//!
//! @param state 1 to turn on picking of points in a tomogram section view, 0 to turn off
void set_tomo_picker_mode_is_active(short int state);

#ifdef USE_PYTHON
//! \brief experimental tomogram spot analysis - internal
//!
//! @param imol_map the map molecule index
//! @param spot_positions a list of dictionaries, each with a "position" key of [x, y, z]
void tomo_map_analysis(int imol_map, PyObject *spot_positions);
//! \brief experimental tomogram spot analysis (version 2) - internal
//!
//! Prints per-section scores for the spot columns and adds the spot positions as
//! a "Picked Points" generic display object.
//!
//! @param imol_map the map molecule index
//! @param spot_positions a list of dictionaries, each with a "position" key of [x, y, z]
void tomo_map_analysis_2(int imol_map, PyObject *spot_positions);
#endif

//! \brief reverse the sign of the map
//!
//! Negative becomes positive and positive becomes negative.
//! Apply an offset so that most of the map is above zero: each value f becomes
//! -f - (mean - 2.5 * variance).
//!
//! @param imol_map the map molecule index
void reverse_map(int imol_map);

//! \brief read positron metadata from two CSV files
//!
//! The metadata are appended to the existing positron metadata.
//!
//! @param z_data the file name of the CSV file of (x, y) latent-space coordinates
//! @param table the file name of the CSV file of the 6 parameters for each point
void read_positron_metadata(const std::string &z_data, const std::string &table);

#ifdef USE_PYTHON
//! \brief make positron maps along a pathway of latent-space points
//!
//! For each point the closest positron metadata point is found and a new map is made
//! (a copy of the first map, regenerated from the base maps weighted by that point's
//! parameters) and contoured at 0.02.
//!
//! @param map_molecule_list_py a list of 6 base map molecule indices
//! @param pathway_points_py a list of [x, y] points
//!
//! @return a list of the new map molecule indices
PyObject *positron_pathway(PyObject *map_molecule_list_py, PyObject *pathway_points_py);
#endif

#ifdef USE_PYTHON
//! \brief show the interactive positron plot dialog
//!
//! @param fn_z_csv the file name of the CSV file of latent-space coordinates
//! @param fn_s_csv the file name of the CSV file of the parameters for each point
//! @param base_map_index_list a list of the base map molecule indices
void positron_plot_py(const std::string &fn_z_csv, const std::string &fn_s_csv,
                      PyObject *base_map_index_list);
#endif

#ifdef USE_PYTHON
//! \brief show the Global Phasing (buster) geometry screen results in a dialog
//!
//! @param imol the model molecule index
//! @param screen_dict a dictionary of screen results with (optional) keys
//!        "Front page", "Bond length", "Bond angle", "Torsion", "Plane",
//!        "Ideal Contact", "Aniso SPH", "Aniso SPH nonb" and "Unhappy Atom"
//!
//! @return always False
PyObject *global_phasing_screen(int imol, PyObject *screen_dict);
#endif


/*! \brief Display a Ramachandran probability surface on a torus as a generic display object.
 *
 * R is the major radius (centre of tube to centre of torus),
 * r is the minor radius (tube radius),
 * height_scale controls the amplitude of the probability displacement.
 * Returns the generic object index, or -1 on failure. */
int show_ramachandran_surface_on_torus(float R, float r, float height_scale);

#ifdef SWIG
#else
//! \brief show the interactive positron plot dialog - internal
//!
//! See positron_plot_py().
void positron_plot_internal(const std::string &fn_z_csv, const std::string &fn_s_csv,
                            const std::vector<int> &base_map_index_list);
#endif

#endif // CC_INTERFACE_HH
