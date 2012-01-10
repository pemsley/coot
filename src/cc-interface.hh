/* src/c-interface.h
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2008, 2009, 2010, 2011, 2012 by The University of Oxford
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

#include "coot-utils.hh"
#include "coot-coord-utils.hh"

#include "dipole.hh"
#include "sequence-assignment.hh" // for residue_range_t

#include "flev-annotations.hh" // animated ligand interactions


namespace coot {

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

   // pisa internal function
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

std::string probe_dots_short_contact_name_to_expanded_name(const std::string &short_name);

// To be used to (typically) get the menu item text label from chain
// option menus (rather than the ugly/broken casting of
// GtkPositionType data.  A wrapper to a static graphics_info_t
// function.
std::string menu_item_label(GtkWidget *menu_item);


/*  ---------------------------------------------------------------------- */
/*                       go to atom   :                                    */
/*  ---------------------------------------------------------------------- */

#ifdef USE_GUILE
// Bernie, no need to pythonize this, it's just to test the return
// values on pressing "next residue" and "previous residue" (you can
// if you wish of course).
//
// Pass the current values, return new values
SCM goto_next_atom_maybe(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
SCM goto_prev_atom_maybe(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif 

#ifdef USE_PYTHON
// but I 'want' to! Needed for python unittest!
PyObject *goto_next_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
PyObject *goto_prev_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif

int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec);
int set_go_to_atom_from_res_spec(const coot::residue_spec_t &spec);
#ifdef USE_GUILE
int set_go_to_atom_from_res_spec_scm(SCM residue_spec);
#endif 
#ifdef USE_PYTHON
int set_go_to_atom_from_res_spec_py(PyObject *residue_spec);
#endif 


// This is to make porting the active atom more easy for Bernhard.
// Return a class rather than a list, and rewrite the active-residue
// function use this atom-spec.  The first value of the pair indicates
// if an atom spec was found.
//
std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec();


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
//! \}

/*  ---------------------------------------------------------------------- */
/*                       map functions:                                    */
/*  ---------------------------------------------------------------------- */
//! \name Extra Map Functions
//! \{

void add_map_colour_mol_menu_item(int imol, const std::string &name,
				  GtkWidget *sub_menu, GtkSignalFunc callback);
/* Actually this function is generic and could be renamed so. */
void add_map_scroll_wheel_mol_menu_item(int imol, 
					const std::string &name,
					GtkWidget *menu, 
					GtkSignalFunc callback);

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


#ifdef USE_GUILE
/*! \brief fit residues

(note: fit to the current-refinement map)
*/
SCM multi_residue_torsion_fit_scm(int imol, SCM residues_specs_scm);
#endif // GUILE
void multi_residue_torsion_fit(int imol, const std::vector<coot::residue_spec_t> &specs);

#ifdef USE_PYTHON
/*! \brief fit residues

(note: fit to the current-refinement map)
*/
PyObject *multi_residue_torsion_fit_py(int imol, PyObject *residues_specs_py);
#endif // PYTHON


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
//! \}


/*  ------------------------------------------------------------------- */
/*                    file selection                                    */
/*  ------------------------------------------------------------------- */

namespace coot {
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
   
   // trivial helper function
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
// return the status and vector of chain-ids of the new chain ids.
// 
std::pair<int, std::vector<std::string> > merge_molecules_by_vector(const std::vector<int> &add_molecules, int imol);

/*  ----------------------------------------------------------------------- */
/*                         Dictionaries                                     */
/*  ----------------------------------------------------------------------- */
//! \name Dictionary Functions
//! \{
/*! \brief return a list of all the dictionaries read */
#ifdef USE_GUILE
SCM dictionaries_read();
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *dictionaries_read_py();
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
PyObject *monomer_restraints_py(const char *monomer_type);
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

/*  ----------------------------------------------------------------------- */
/*                  scripting                                               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_GUILE
SCM safe_scheme_command_test(const char *cmd);
SCM safe_scheme_command(const std::string &scheme_command);
#else 
void safe_scheme_command(const std::string &scheme_command); /* do nothing */
#endif // USE_GUILE

#ifdef USE_PYTHON
void safe_python_command(const std::string &python_command); 
void safe_python_command_by_char_star(const char *python_command);
PyObject *py_clean_internal(PyObject *obj);
PyObject *safe_python_command_with_return(const std::string &python_cmd);
PyObject *safe_python_command_test(const char *cmd);
void safe_python_command_with_unsafe_thread(const char *cmd);
#endif // PYTHON
/*  Is this a repeat of something?  I don't know. */
void run_generic_script(const std::vector<std::string> &cmd_strings);

#ifdef USE_GUILE
SCM coot_has_python_p();
#endif

#ifdef USE_PYTHON
PyObject *coot_has_guile();
#endif

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
SCM scm_residue(const coot::residue_spec_t &res);
#endif

#ifdef USE_PYTHON
// Return a list describing a residue like that returned by
// residues-matching-criteria [return_val, chain_id, resno, ins_code]
// This is a library function really.  There should be somewhere else to put it.
// It doesn't need expression at the scripting level.
// return a null list on problem
PyObject *py_residue(const coot::residue_spec_t &res);
#endif

#ifdef USE_GUILE
coot::residue_spec_t residue_spec_from_scm(SCM residue_in);
#endif

#ifdef USE_PYTHON
coot::residue_spec_t residue_spec_from_py(PyObject *residue_in);
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


#ifdef USE_GUILE
//! \brief Return a list of atom info for each atom in the specified residue
//! 
//! output is like this:
//! (list
//!    (list (list atom-name alt-conf)
//!          (list occ temp-fact element)
//!          (list x y z)))
//! 
SCM residue_info(int imol, const char* chain_id, int resno, const char *ins_code);
SCM residue_name(int imol, const char* chain_id, int resno, const char *ins_code);

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

//! \brief return the specs of the closest atom in imolth molecule
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf (list x y z)) for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, no account is
//! taken of the displayed state of the molecule).  If there is no
//! atom, or if imol is not a valid model molecule, return scheme false.
//! 
SCM closest_atom(int imol);

//! \brief return residues near residue
//! 
//! Return residue specs for residues that have atoms that are
//! closer than radius Angstroems to any atom in the residue
//! specified by res_in.
//! 
SCM residues_near_residue(int imol, SCM residue_in, float radius);

//! return residues within radius of pos (x,y,z) position
//! 
//! Return a list of pairs of (imol, residue_spec).
//! pos is a list of 3 numbers.  (get imol from active-atom)
//! 
SCM residues_near_position_scm(int imol, SCM pos, float radius);

#endif	/* USE_GUILE */

//! find the active residue, find the near residues (within radius) 
//! create a new molecule, run reduce on that, import hydrogens from
//! the result and apply them to the molecule of the active residue.
void hydrogenate_region(float radius);

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

//! \brief return closest atom in imolth molecule
// 
//! Return a list of [imol, chain-id, resno, ins-code, atom-name,
//! alt-conf, [x, y, z]] for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, no account is
//! taken of the displayed state of the molecule).  If there is no
//! atom, or if imol is not a valid model molecule, return #f.
// 
PyObject *closest_atom_py(int imol);

// Return residue specs for residues that have atoms that are
// closer than radius Angstroems to any atom in the residue
// specified by res_in.
// 
PyObject *residues_near_residue_py(int imol, PyObject *residue_in, float radius);
PyObject *residues_near_position_py(int imol, PyObject *pos_in, float radius);

#endif // USE_PYTHON

//! \}

//! \name status bar string functions
//! \{
// status bar atom info text here?!
std::string atom_info_as_text_for_statusbar(int atom_index, int imol, const char *symmetry_text = "");
//! \}


/*  ----------------------------------------------------------------------- */
/*                  Refinement                                              */
/*  ----------------------------------------------------------------------- */

//! \name Refinement with specs
//! \{


//! \brief refine a zone, allowing the specification of insertion
//!  codes for the residues too.
//! 
//! presumes that imol_Refinement_Map has been set
#ifdef USE_GUILE
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

//! \}


/*  ----------------------------------------------------------------------- */
/*                  rigid body fitting (multiple residue ranges)            */
/*  ----------------------------------------------------------------------- */
// return 0 on fail to refine (no sensible place to put atoms) and 1
// on fitting happened.
int rigid_body_fit_with_residue_ranges(int imol, const std::vector<coot::residue_range_t> &ranges);

/*  ----------------------------------------------------------------------- */
/*                  check water baddies                                     */
/*  ----------------------------------------------------------------------- */

std::vector<coot::atom_spec_t>
check_waters_baddies(int imol, float b_factor_lim, float map_sigma_lim, float min_dist, float max_dist, short int part_occ_contact_flag, short int zero_occ_flag, short int logical_operator_and_or_flag);

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
/* return the chain id of the water chain. Raw interface */
PyObject *water_chain_py(int imol);
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
#endif

#ifdef USE_PYTHON
//! \brief for the given residue, spin the atoms in moving_atom_list
//!   around the bond defined by direction_atoms_list looking for the best
//!   fit to density of imom_map map of the first atom in
//!   moving_atom_list.  Works (only) with atoms in altconf ""
void spin_search_py(int imol_map, int imol, const char *chain_id, int resno, const char *ins_code, PyObject *direction_atoms_list, PyObject *moving_atoms_list);
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
std::vector<int> execute_ligand_search_internal();
coot::graph_match_info_t
overlap_ligands_internal(int imol_ligand, int imol_ref, const char *chain_id_ref,
			 int resno_ref, bool apply_rtop_flag);

/*  ----------------------------------------------------------------------- */
//                  animated ligand interactions
/*  ----------------------------------------------------------------------- */
void add_animated_ligand_interaction(int imol, const coot::fle_ligand_bond_t &lb);


/*  ----------------------------------------------------------------------- */
/*                  Cootaneer                                               */
/*  ----------------------------------------------------------------------- */
int cootaneer_internal(int imol_map, int imol_model, coot::atom_spec_t &atom_spec);

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

CMMDBManager *new_molecule_by_symmetry_matrix_from_molecule(CMMDBManager *mol,
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
#ifdef USE_LIBCURL
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
#endif /* USE_LIBCURL */


/*  ----------------------------------------------------------------------- */
/*                  Functions for FLEV layout callbacks                     */
/*  ----------------------------------------------------------------------- */
// orient the graphics somehow so that the interaction between
// central_residue and neighbour_residue is perpendicular to screen z.
void orient_view(int imol,
		 const coot::residue_spec_t &central_residue_spec, // ligand typically
		 const coot::residue_spec_t &neighbour_residue_spec);


/*  return a list of chiral centre ids as determined from topological
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
		 int imol_map, int nfrags);
// so that we can create a "original loop" molecule from the atom
// specs picked (i.e. the atom selection string should extend over the
// range from the smallest residue number to the largest (in the same
// chain)).
std::string
protein_db_loop_specs_to_atom_selection_string(const std::vector<coot::residue_spec_t> &specs);
#ifdef USE_GUILE
SCM protein_db_loops_scm(int imol_coords, SCM residues_specs, int imol_map, int nfrags);
#endif 
#ifdef USE_PYTHON 
PyObject *protein_db_loops_py(int imol_coords, PyObject *residues_specs, int imol_map, int nfrags);
#endif 
/* \} */


/* ------------------------------------------------------------------------- */
/*                      HOLE                                                 */
/* ------------------------------------------------------------------------- */
/*! \name Coot's Hole implementation */

/*! \brief starting piont and end point, colour map multiplier and
  shall the probe radius graph be shown (dummy value currently). */
void hole(int imol,
	  float start_x, float start_y, float start_z,
	  float   end_x, float   end_y, float   end_z, 
	  float colour_map_multiplier, float colour_map_offset,
	  int n_runs, bool show_probe_radius_graph_flag);


// GUI stuff
void probe_radius_graph_close_callback( GtkWidget *button,
 					GtkWidget *dialog);
void show_hole_probe_radius_graph(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length);
void show_hole_probe_radius_graph_basic(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length);
void show_hole_probe_radius_graph_goocanvas(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length);


/* ------------------------------------------------------------------------- */
/*                      LINKs                                                */
/* ------------------------------------------------------------------------- */
void
make_link(int imol, coot::atom_spec_t &spec_1, coot::atom_spec_t &spec_2,
	  const std::string &link_name, float length);
#ifdef USE_GUILE
void make_link_scm(int imol, SCM spec_1, SCM spec_2, const std::string&link_name, float length);
#endif
#ifdef USE_PYTHON
void make_link_py(int imol, PyObject *spec_1, PyObject *spec_2, const std::string&link_name, float length);
#endif


/*! \name  Drag and Drop */
// \{
//! \brief handle the string that get when a file or URL is dropped.
int handle_drag_and_drop_string(const std::string &uri);
// \}



/* ------------------------------------------------------------------------- */
/*                      prodrg import function                               */
/* ------------------------------------------------------------------------- */
// the function passed to lbg, so that it calls it when a new
// prodrg-in.mdl file has been made.  We no longer have a timeout
// function waiting for prodrg-in.mdl to be updated/written.
// 
void prodrg_import_function(std::string file_name);



/* ------------------------------------------------------------------------- */
/*                       SBase import function                               */
/* ------------------------------------------------------------------------- */
// the function passed to lbg, so that it calls it when a new
// SBase comp_id is required.  We no longer have a timeout
// function waiting for prodrg-in.mdl to be updated/written.
// 
void sbase_import_function(std::string comp_id);



/*  ----------------------------------------------------------------------- */
/*                  GUIL Utility Functions                                  */
/*  ----------------------------------------------------------------------- */
void simple_text_dialog(const std::string &dialog_title, const std::string &text,
			std::pair<int, int> geom);

// gui nuts and bolts
void on_simple_text_dialog_close_button_pressed( GtkWidget *button,
						 GtkWidget *dialog);


/*  ----------------------------------------------------------------------- */
/*                  Utility Functions                                       */
/*  ----------------------------------------------------------------------- */
// These functions are for storing the molecule number and (some other
// number) as an int and used with GPOINTER_TO_INT and GINT_TO_POINTER.
int encode_ints(int i1, int i2);
std::pair<int, int> decode_ints(int i);

#ifdef USE_GUILE
std::vector<coot::residue_spec_t> scm_to_residue_specs(SCM s);
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
