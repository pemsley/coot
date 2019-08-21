/* src/cc-interface.hh
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2008, 2009, 2010, 2011, 2012 by The University of Oxford
 * Copyright 2015 by Medical Research Council
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifndef CC_INTERFACE_HH
#define CC_INTERFACE_HH

#include <gtk/gtk.h>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-density-stats.hh"

#include "ligand/dipole.hh"
#include "high-res/sequence-assignment.hh" // for residue_range_t

#include "coords/mmdb-extras.h"
#include "coords/mmdb-crystal.h"

#include "pli/flev-annotations.hh" // animated ligand interactions
#include "named-rotamer-score.hh"

#include "coords/phenix-geo.hh"

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
      alias_path_t(int index_in, const std::string &s_in, bool flag_in) {
	 index = index_in;
	 s = s_in;
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


std::vector<std::string> filtered_by_glob(const std::string &pre_directory, 
					  int data_type);
/*  Return 1 if search appears in list, 0 if not) */
short int 
string_member(const std::string &search, const std::vector<std::string> &list); 
bool compare_strings(const std::string &a, const std::string &b); 

/*
std::string pre_directory_file_selection(GtkWidget *sort_button);
void filelist_into_fileselection_clist(GtkWidget *fileselection, const std::vector<std::string> &v);

GtkWidget *wrapped_nothing_bad_dialog(const std::string &label);

std::pair<short int, float> float_from_entry(GtkWidget *entry);
std::pair<short int, int>   int_from_entry(GtkWidget *entry);

void
add_validation_mol_menu_item(int imol, const std::string &name, GtkWidget *menu, GtkSignalFunc callback);
void create_initial_validation_graph_submenu_generic(GtkWidget *window1,
						     const std::string &menu_name,
						     const std::string &sub_menu_name);

// To be used to (typically) get the menu item text label from chain
// option menus (rather than the ugly/broken casting of
// GtkPositionType data.  A wrapper to a static graphics_info_t
// function.
std::string menu_item_label(GtkWidget *menu_item);

*/

std::string probe_dots_short_contact_name_to_expanded_name(const std::string &short_name);


/*  ---------------------------------------------------------------------- */
/*                       go to atom   :                                    */
/*  ---------------------------------------------------------------------- */

void set_rotation_centre(const clipper::Coord_orth &pos);

#ifdef USE_GUILE
//
// Pass the current values, return new values
SCM goto_next_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
SCM goto_prev_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif 

#ifdef USE_PYTHON

PyObject *goto_next_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
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
// return a tuple of (Py_Bool (number, atom_spec))
PyObject *active_atom_spec_py();
#endif // USE_PYTHON


/*  ---------------------------------------------------------------------- */
/*                       symmetry functions:                               */
/*  ---------------------------------------------------------------------- */
// get the symmetry operators strings for the given molecule
//
#ifdef USE_GUILE

//! \name  More Symmetry Functions
//! \{

//! \brief return the symmetry of the imolth molecule
//! 
//!   Return as a list of strings the symmetry operators of the
//!   given molecule. If imol is a not a valid molecule, return an empty
//!   list.*/
SCM get_symmetry(int imol);
#endif // USE_GUILE

#ifdef USE_PYTHON
// return a python object as a list (or some other python container)
PyObject *get_symmetry_py(int imol);
#endif // USE_PYTHON

//! \brief return 1 if this residue clashes with the symmetry-related
//!  atoms of the same molecule.
//! 
//! 0 means that it did not clash,
//! -1 means that the residue or molecule could not be found or that there
//!    was no cell and symmetry.
int clashes_with_symmetry(int imol, const char *chain_id, int res_no, const char *ins_code,
			  float clash_dist);
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
std::vector<int> auto_read_make_and_draw_maps_from_mtz(const char *filename);
std::vector<int> auto_read_make_and_draw_maps_from_cns(const char *filename);

/* ----- remove wiget functions from this header GTK-FIXME
void add_map_colour_mol_menu_item(int imol, const std::string &name,
				  GtkWidget *sub_menu, GtkSignalFunc callback);
 Actually this function is generic and could be renamed so.
void add_map_scroll_wheel_mol_menu_item(int imol, 
					const std::string &name,
					GtkWidget *menu, 
					GtkSignalFunc callback);
*/

//! \brief make a sharpened or blurred map
//!
//! blurred maps are generated by using a positive value of b_factor.
//!
//! @return the index of the map created by applying a b-factor
//!        to the given map. Return -1 on failure.
int sharpen_blur_map(int imol_map, float b_factor);

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
void multi_sharpen_blur_map_py(int imol_map, PyObject *b_factors_list);
#endif

#ifdef USE_PYTHON
PyObject *amplitude_vs_resolution_py(int mol_map);
#endif

#ifdef USE_GUILE
SCM amplitude_vs_resolution_scm(int mol_map);
#endif

//! \brief b-factor from map
//!
//! calculate structure factors and use the amplitudes to estimate
//! the B-factor of the data using a wilson plot using a low resolution
//! limit of 4.5A.
//! @return -1 when given a bad map or there were no data beyond 4.5A
//!
float b_factor_from_map(int imol_map);

//! \brief return the colour triple of the imolth map
//! 
//! (e.g.: (list 0.4 0.6 0.8). If invalid imol return scheme false.
//! 
#ifdef USE_GUILE
SCM map_colour_components(int imol);
#endif // GUILE

#ifdef USE_PYTHON
//! \brief return the colour triple of the imolth map
//
//! e.g.: [0.4, 0.6, 0.8]. If invalid imol return Py_False.
//
PyObject *map_colour_components_py(int imol);
#endif // PYTHON
//! \}


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
#include "merge-molecule-results-info-t.hh"
// return the status and vector of chain-ids of the new chain ids.
//
std::pair<int, std::vector<merge_molecule_results_info_t> > merge_molecules_by_vector(const std::vector<int> &add_molecules, int imol);

/*  ----------------------------------------------------------------------- */
/*                         Dictionaries                                     */
/*  ----------------------------------------------------------------------- */
//! \name Dictionary Functions
//! \{
std::vector<std::string> dictionary_entries();
void debug_dictionary();
// this can throw an exception
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
PyObject *SMILES_for_comp_id_py(const std::string &comp_id);
#endif // PYTHON
//! \}


/*  ----------------------------------------------------------------------- */
/*                         Restraints                                       */
/*  ----------------------------------------------------------------------- */
#ifdef USE_GUILE
//! \name  Restraints Interface
/// \{
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

#ifdef USE_GUILE
//! \name Atom Information functions
//! \{
//! \brief output atom info in a scheme list for use in scripting
//!
//! in this format (list occ temp-factor element x y z).  Return empty
//! list if atom not found. */
SCM atom_info_string_scm(int imol, const char *chain_id, int resno,
			 const char *ins_code, const char *atname,
			 const char *altconf);
#endif // USE_GUILE

/*! \brief return the rename from a residue serial number

   @return blank ("") on failure. */
std::string resname_from_serial_number(int imol, const char *chain_id, int serial_num);

//! \brief return the residue name of the specified residue
std::string residue_name(int imol, const std::string &chain_id, int resno, const std::string &ins_code);

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

//! \brief resdiues near residue
//!
//! @return residues within radius of pos (x,y,z) position
//! 
//! Return a list of pairs of (imol, residue_spec).
//! pos is a list of 3 numbers.  (get imol from active-atom)
//! 
SCM residues_near_position_scm(int imol, SCM pos, float radius);

#endif	/* USE_GUILE */

//! \brief find the active residue, find the near residues (within radius) 
//! create a new molecule, run reduce on that, import hydrogens from
//! the result and apply them to the molecule of the active residue.
void hydrogenate_region(float radius);

//! \brief Add hydrogens to imol from the given pdb file
void add_hydrogens_from_file(int imol, std::string pdb_with_Hs_file_name);

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
//! Return a list of atom info for each atom in the specified residue:
//
//! output is like this:
//! [
//!     [[atom-name,alt-conf]
//!      [occ,temp_fact,element]
//!      [x,y,z]]]
//!
PyObject *residue_info_py(int imol, const char* chain_id, int resno, const char *ins_code);
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
//! Return a list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf, [x, y, z]] for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, potential CA 
//! substition is not performed).  If there is no atom, or if imol is 
//! not a valid model molecule, return False.
//! 
PyObject *closest_atom_simple_py();

//! \brief return closest atom in imolth molecule
// 
//! Return a list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf, [x, y, z]] for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, no account is
//! taken of the displayed state of the molecule).  If there is no
//! atom, or if imol is not a valid model molecule, return False.
// 
PyObject *closest_atom_py(int imol);

//! \brief return the specs of the closest atom to the centre of the screen
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf (list x y z)) for atom that is closest to the screen
//! for displayed molecules. If there is no atom, return scheme false.
//! Don't choose the CA of the residue if there is a CA in the residue
//! of the closest atom
PyObject *closest_atom_raw_py();


//! \brief
// Return residue specs for residues that have atoms that are
// closer than radius Angstroems to any atom in the residue
// specified by residue_in.
// 
PyObject *residues_near_residue_py(int imol, PyObject *residue_in, float radius);

//! \brief
// Return residue specs for residues that have atoms that are
// closer than radius Angstroems to any atom in the residues
// specified by residues_in.
// 
PyObject *residues_near_residues_py(int imol, PyObject *residues_in, float radius);

//! \brief
//! Return residue specs for residues that have atoms that are
//! closer than radius Angstroems to the given position.
//!
PyObject *residues_near_position_py(int imol, PyObject *pos_in, float radius);

//! \brief return a Python object for the bonds
//
PyObject *get_bonds_representation(int imol);

//! \brief return a Python object for the radii of the atoms in the dictionary
//
PyObject *get_dictionary_radii();

//! \brief return a Python object for the representation of bump and hydrogen bonds of
//          the specified residue
PyObject *get_environment_distances_representation_py(int imol, PyObject *residue_spec_py);

//! \brief return a Python object for the intermediate atoms bonds
//
PyObject *get_intermediate_atoms_bonds_representation();

//! \brief return the continue-updating-refinement-atoms state
//
// 0 means off, 1 means on.
// Given the current wiring of the refinement, this is always 0, i.e. refine_residues()
// will return only after the atoms have finished moving.
int get_continue_updating_refinement_atoms_state();

#endif // USE_PYTHON

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


//! \brief
//! regularize the given residues
//!
void regularize_residues(int imol, const std::vector<coot::residue_spec_t> &residues);

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

void set_show_intermediate_atoms_rota_markup(short int state);
void set_show_intermediate_atoms_rama_markup(short int state);

void set_cryo_em_refinement(bool mode);
bool get_cryo_em_refinement();

#ifdef USE_PYTHON
void register_post_intermediate_atoms_moved_hook(PyObject *function_name);
#endif

//! \}


/*  ----------------------------------------------------------------------- */
/*                  rigid body fitting (multiple residue ranges)            */
/*  ----------------------------------------------------------------------- */

//! \brief
//! return 0 on fail to refine (no sensible place to put atoms) and 1
//! on fitting happened.
int rigid_body_fit_with_residue_ranges(int imol, const std::vector<coot::residue_range_t> &ranges);

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
PyObject *find_blobs_py(int imol_model, int imol_map, float cut_off_density_level);
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
int make_variance_map(std::vector<int> map_molecule_number_vec);
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
//! \brief for the given residue, spin the atoms in moving_atom_list...
//!
//!   around the bond defined by direction_atoms_list looking for the best
//!   fit to density of imom_map map of the first atom in
//!   moving_atom_list.  Works (only) with atoms in altconf ""
void spin_search_py(int imol_map, int imol, const char *chain_id, int resno, const char *ins_code, PyObject *direction_atoms_list, PyObject *moving_atoms_list);
//! \brief Spin N and CB (and the rest of the side chain if extant)
//!
//!  Sometime on N-terminal addition, then N ends up pointing the wrong way.
//!  The allows us to (more or less) interchange the positions of the CB and the N.
//!  angle is in degrees.
//!
void spin_N_py(int imol, PyObject *residue_spec, float angle);

//! \brief Spin search the density based on possible positions of CG of a side-chain
PyObject *CG_spin_search_py(int imol_model, int imol_map);

#endif

//! \}


/*  ----------------------------------------------------------------------- */
/*                  monomer lib                                             */
/*  ----------------------------------------------------------------------- */
std::vector<std::pair<std::string, std::string> > monomer_lib_3_letter_codes_matching(const std::string &search_string, short int allow_minimal_descriptions_flag);

void on_monomer_lib_search_results_button_press (GtkButton *button, gpointer user_data);
void on_monomer_lib_sbase_molecule_button_press (GtkButton *button, gpointer user_data);

/*  ----------------------------------------------------------------------- */
/*                  mutate                                                  */
/*  ----------------------------------------------------------------------- */
int mutate_internal(int ires, const char *chain_id,
		    int imol, std::string &target_res_type);
/* a function for multimutate to make a backup and set
   have_unsaved_changes_flag themselves */

/*  ----------------------------------------------------------------------- */
/*                  ligands                                                 */
/*  ----------------------------------------------------------------------- */
coot::graph_match_info_t
overlap_ligands_internal(int imol_ligand, int imol_ref, const char *chain_id_ref,
			 int resno_ref, bool apply_rtop_flag);


/*  ----------------------------------------------------------------------- */
/*                  conformers (part of ligand search)                      */
/*  ----------------------------------------------------------------------- */

#ifdef USE_GUILE
/*! \brief make conformers of the ligand search molecules, each in its
  own molecule.  

Don't search the density.

//! \name Extra Ligand Functions
//! \{

Return a list of new molecule numbers */
SCM ligand_search_make_conformers_scm();
#endif 

#ifdef USE_PYTHON
PyObject *ligand_search_make_conformers_py();
#endif

std::vector<int> ligand_search_make_conformers_internal();

//! \}

/*  ----------------------------------------------------------------------- */
//                  animated ligand interactions
/*  ----------------------------------------------------------------------- */
void add_animated_ligand_interaction(int imol, const coot::fle_ligand_bond_t &lb);


/*  ----------------------------------------------------------------------- */
/*                  Cootaneer                                               */
/*  ----------------------------------------------------------------------- */
int cootaneer_internal(int imol_map, int imol_model, const coot::atom_spec_t &atom_spec);

#ifdef USE_GUILE
//! \name Dock Sidechains
//! \{
//! \brief cootaneer (i.e. dock sidechains onto mainchain model)
//! 
//! atom_in_fragment_atom_spec is any atom spec in the fragment that should be
//! docked with sidechains.
//!
//! @return the success status (0 is fail).
int cootaneer(int imol_map, int imol_model, SCM atom_in_fragment_atom_spec);
#endif 

#ifdef USE_PYTHON
int cootaneer_py(int imol_map, int imol_model, PyObject *atom_in_fragment_atom_spec);
#endif

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
#endif	/* USE_GUILE */

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

/*  ----------------------------------------------------------------------- */
/*                  Abstraction of New molecule by symmetry functions       */
/*  ----------------------------------------------------------------------- */

mmdb::Manager *new_molecule_by_symmetry_matrix_from_molecule(mmdb::Manager *mol,
							    double m11, double m12, double m13, 
							    double m21, double m22, double m23, 
							    double m31, double m32, double m33, 
							    double tx, double ty, double tz,
							    int pre_shift_to_origin_na,
							    int pre_shift_to_origin_nb,
							    int pre_shift_to_origin_nc);


/*  ----------------------------------------------------------------------- */
/*                  LIBCURL/Download                                        */
/*  ----------------------------------------------------------------------- */

/*! \brief if possible, read in the new coords getting coords via web.

(no return value because get-url-str does not return one).
 */
   void get_coords_for_accession_code(const std::string &code);

#ifdef USE_LIBCURL
// return 0 on success.
int coot_get_url(const char *url, const char *file_name);
int coot_get_url_and_activate_curl_hook(const char *url, const char *file_name, short int do_hook_flag);
#ifdef USE_GUILE
// this handles URLs that are strings, not binaries. 
SCM coot_get_url_as_string(const char *url);
// for the callback of the update binary progress bar.  How much done
// is the file that I am downloading?
SCM curl_progress_info(const char *file_name);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
// this handles URLs that are strings, not binaries. 
PyObject *coot_get_url_as_string_py(const char *url);
// for the callback of the update binary progress bar.  How much done
// is the file that I am downloading? Not absolutely required for python
PyObject *curl_progress_info_py(const char *file_name);
#endif /* USE_PYTHON */
// internal use
size_t write_coot_curl_data(void *buffer, size_t size, size_t nmemb, void *userp);
// internal use
size_t write_coot_curl_data_to_file(void *buffer, size_t size, size_t nmemb, void *userp);
// internal use (strings, not binaries).
std::string coot_get_url_as_string_internal(const char *url);
void *wrapped_curl_easy_perform(void *data);
void stop_curl_download(const char *file_name); // stop curling the to file_name;

std::string get_drug_mdl_via_wikipedia_and_drugbank(std::string drugname);

#endif /* USE_LIBCURL */


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
/*                  Pisa internal                                           */
/*  ----------------------------------------------------------------------- */
clipper::Coord_orth 
make_complementary_dotted_surfaces(int imol_1, int imol_2, 
				   std::vector<coot::residue_spec_t> &r1, 
				   std::vector<coot::residue_spec_t> &r2);
#ifdef USE_GUILE
std::vector<coot::residue_spec_t>
residue_records_list_scm_to_residue_specs(SCM mol_1_residue_records,
					  const std::string &chain_id);
SCM symbol_value_from_record(SCM record_1, const std::string &symbol);
#endif 
#ifdef USE_PYTHON
std::vector<coot::residue_spec_t>
residue_records_list_py_to_residue_specs(PyObject *mol_1_residue_records,
					 const std::string &chain_id);
//PyObject *symbol_value_from_record(PyObject *record_1, const std::string &symbol);
#endif 

void
add_generic_object_bond(int imol1, int imol2,
			const coot::atom_spec_t &atom_spec_1,
			const coot::atom_spec_t &atom_spec_2,
			int generic_object_number,
			const std::string &colour);

void
pisa_interfaces_display_only(int imol_1, int imol_2, clipper::Coord_orth centre_pt);
std::string untangle_mmdb_chain_id_string(const std::string &mmdb_chain_id_in);


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
//! \brief return the scores of the rotamers for this residue.
//
// The density fit score is for side-chain atoms.
// return a list (possibly empty).
// 
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
// return a list (possibly empty)
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

/*! \brief starting piont and end point, colour map multiplier and
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


// GUI stuff
void probe_radius_graph_close_callback( GtkWidget *button,
 					GtkWidget *dialog);
void show_hole_probe_radius_graph(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length);
void show_hole_probe_radius_graph_basic(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length);
void show_hole_probe_radius_graph_goocanvas(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length);


/* ------------------------------------------------------------------------- */
/*                      LINKs                                                */
/* ------------------------------------------------------------------------- */

//! \brief make a link between the specified atoms
void
make_link(int imol, coot::atom_spec_t &spec_1, coot::atom_spec_t &spec_2,
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


/* ------------------------------------------------------------------------- */
/*                      Drag and drop                                        */
/* ------------------------------------------------------------------------- */

/*! \name  Drag and Drop Functions */
// \{
//! \brief handle the string that get when a file or URL is dropped.
int handle_drag_and_drop_string(const std::string &uri);
// \}


/* ------------------------------------------------------------------------- */
/*                      Map Contours                                         */
/* ------------------------------------------------------------------------- */

#ifdef USE_PYTHON
/*! \name Map Contouring Functions */
// \{
//! \brief return two lists: a list of vertices and a list of indices for connection
PyObject *map_contours(int imol, float contour_level);
// \}
#endif // USE_PYTHON



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
SCM map_to_model_correlation_stats_scm(int imol,
				       SCM residue_specs,
				       SCM neighb_residue_specs,
				       unsigned short int atom_mask_mode,
				       int imol_map);
#endif

#ifdef USE_PYTHON
PyObject *map_to_model_correlation_py(int imol,
				      PyObject *residue_specs,
				      PyObject *neighb_residue_specs,
				      unsigned short int atom_mask_mode,
				      int imol_map);
PyObject *map_to_model_correlation_stats_py(int imol,
				      PyObject *residue_specs,
				      PyObject *neighb_residue_specs,
				      unsigned short int atom_mask_mode,
				      int imol_map);
#endif

//! \brief atom-mask-mode is as follows:
// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-excluding CB if is standard amino-acid, else all atoms
// 4: main-chain atoms if is standard amino-acid, else nothing
// 5: side-chain atoms if is standard amino-acid, else nothing
// 10: atom radius is dependent atom atom B-factor
float
map_to_model_correlation(int imol, 
			 const std::vector<coot::residue_spec_t> &residue_specs,
			 const std::vector<coot::residue_spec_t> &neigh_residue_specs,
			 unsigned short int atom_mask_mode,
			 int imol_map);

//! \brief map to model density correlation stats
//! 
//! \brief atom-mask-mode is as follows:
// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-excluding CB if is standard amino-acid, else all atoms
// 4: main-chain atoms if is standard amino-acid, else nothing
// 5: side-chain atoms if is standard amino-acid, else nothing
//
coot::util::density_correlation_stats_info_t
map_to_model_correlation_stats(int imol,
			       const std::vector<coot::residue_spec_t> &residue_specs,
			       const std::vector<coot::residue_spec_t> &neigh_residue_specs,
			       unsigned short int atom_mask_mode,
			       int imol_map);
#ifndef SWIG

//! \brief map to model density correlation, reported per residue
//! 
//! \brief atom-mask-mode is as follows:
// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-excluding CB if is standard amino-acid, else all atoms
// 4: main-chain atoms if is standard amino-acid, else nothing
// 5: side-chain atoms if is standard amino-acid, else nothing
//
std::vector<std::pair<coot::residue_spec_t,float> >
map_to_model_correlation_per_residue(int imol, const std::vector<coot::residue_spec_t> &specs,
				     unsigned short int atom_mask_mode,
				     int imol_map);

//! \brief map to model density statistics, reported per residue
std::map<coot::residue_spec_t, coot::util::density_stats_info_t>
map_to_model_correlation_stats_per_residue(int imol,
					   const std::vector<coot::residue_spec_t> &residue_specs,
					   unsigned short int atom_mask_mode,
					   float atom_radius_for_masking,
					   int imol_map);
#endif // not for swigging.

#ifdef USE_GUILE
//! \brief map to model correlation
SCM
map_to_model_correlation_per_residue_scm(int imol, SCM residue_specs,
					 unsigned short int atom_mask_mode,
					 int imol_map);

//! \brief map to model stats
SCM
map_to_model_correlation_stats_per_residue_scm(int imol,
					       SCM residue_specs_scm,
					       unsigned short int atom_mask_mode,
					       float atom_radius_for_masking,
					       int imol_map);

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
float density_score_residue_scm(int imol, SCM residue_spec, int imol_map);
#endif 
#ifdef USE_PYTHON
float density_score_residue_py(int imol, PyObject *residue_spec, int imol_map);
#endif 
#endif 

/*! \brief simple density score for given residue (over-ridden by scripting function) */
float density_score_residue(int imol, const char *chain_id, int res_no, const char *ins_code, int imol_map);


#ifdef USE_GUILE
/*! \brief return sigma for the given map.  Return scheme False if not
  a valid map molecule number. */
SCM map_mean_scm(int imol);
SCM map_sigma_scm(int imol);
/*! \brief return either scheme false on non-a-map or list (mean, standard-deviation, skew, kurtosis) */
SCM map_statistics_scm(int imol);
#endif
#ifdef USE_PYTHON
/*! \brief return sigma for the given map.  Return Python False if not
  a valid map molecule number. */
PyObject *map_mean_py(int imol);
PyObject *map_sigma_py(int imol);
PyObject *map_statistics_py(int imol);
#endif /*USE_PYTHON */


//! \}

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
PyObject *molecule_atom_overlaps_py(int imol);
#endif // USE_PYTHON
#ifdef USE_GUILE
SCM molecule_atom_overlaps_scm(int imol);
#endif // USE_GUILE

/* ------------------------------------------------------------------------- */
/*                      prodrg import function                               */
/* ------------------------------------------------------------------------- */
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
					   const coot::phenix_geo_bonds &g);
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

#endif // CC_INTERFACE_HH
