/* src/c-interface.h
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#include "coot-coord-utils.hh"


std::vector<std::string> filtered_by_glob(const std::string &pre_directory, 
					  int data_type);
/*  Return 1 if search appears in list, 0 if not) */
short int 
string_member(const std::string &search, const std::vector<std::string> &list); 
bool compare_strings(const std::string &a, const std::string &b); 


std::string pre_directory_file_selection(GtkWidget *sort_button);
GtkWidget *wrapped_nothing_bad_dialog(const std::string &label);

std::pair<short int, float> float_from_entry(GtkWidget *entry);
std::pair<short int, int>   int_from_entry(GtkWidget *entry);

void
add_validation_mol_menu_item(int imol, const std::string &name, GtkWidget *menu, GtkSignalFunc callback);
void create_initial_validation_graph_submenu_generic(GtkWidget *window1,
						     const std::string &menu_name,
						     const std::string &sub_menu_name);

std::string probe_dots_short_contact_name_to_expanded_name(const std::string &short_name);


/*  ---------------------------------------------------------------------- */
/*                       go to atom   :                                    */
/*  ---------------------------------------------------------------------- */

int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec);

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

/*! \brief Return as a list of strings the symmetry operators of the
  given molecule. If imol is a not a valid molecule, return an empty
  list.*/
SCM get_symmetry(int imol);
#else 
#ifdef USE_PYTHON
// return a python object as a list (or some other python container)
#endif 
#endif

/*  ---------------------------------------------------------------------- */
/*                       map functions:                                    */
/*  ---------------------------------------------------------------------- */
void add_map_colour_mol_menu_item(int imol, const std::string &name,
				  GtkWidget *sub_menu, GtkSignalFunc callback);
/* Actually this function is generic and could be renamed so. */
void add_map_scroll_wheel_mol_menu_item(int imol, 
					const std::string &name,
					GtkWidget *menu, 
					GtkSignalFunc callback);

//! \brief return the colour of the imolth map (e.g.: (list 0.4 0.6
//0.8). If invalid imol return #f.
// 
#ifdef USE_GUILE
SCM map_colour_components(int imol);
#endif

/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */
// if swap_map_colours_post_refmac_flag is not 1 thenn imol_refmac_map is ignored.
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
		    std::string refmac_count_string, 
		    int swap_map_colours_post_refmac_flag,
		    int imol_refmac_map,
		    int diff_map_flag,
		    int phase_combine_flag,
		    std::string phib_string,
		    std::string fom_string,
		    std::string ccp4i_project_dir);


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

/*  -------------------------------------------------------------------- */
/*                     history                                           */
/*  -------------------------------------------------------------------- */

namespace coot { 
   class command_arg_t {
   public:
      enum coot_script_arg_type{INT, FLOAT, STRING};
      command_arg_t(int iin) {
	 i = iin;
	 type = INT;
      }
      command_arg_t(float fin) {
	 f = fin;
	 type = FLOAT;
      }
      command_arg_t(const clipper::String &sin) {
	 s = sin;
	 type = STRING;
      }
      command_arg_t(const char *sin) {
	 s = sin;
	 type = STRING;
      }
      coot_script_arg_type type;
      float f;
      int i;
      clipper::String s;
      std::string as_string() const {
	 std::string os("unknown-arg-type");
	 if (type == INT)
	    os = clipper::String(i);
	 if (type == FLOAT)
	    os = clipper::String(f);
	 if (type == STRING)
	    os = s;
	 return os;
      }
   };
}
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
int merge_molecules(const std::vector<int> &add_molecules, int imol);


/*  ----------------------------------------------------------------------- */
/*                         Dictionaries                                     */
/*  ----------------------------------------------------------------------- */
/*! \brief return a list of all the dictionaries read */
#ifdef USE_GUILE
SCM dictionaries_read();
#endif // USE_GUILE


/*  ----------------------------------------------------------------------- */
/*                  scripting                                               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_GUILE
SCM safe_scheme_command_test(const char *cmd);
SCM safe_scheme_command(const std::string &scheme_command);
#else 
void safe_scheme_command(const std::string &scheme_command); /* do nothing */
#endif // USE_GUILE

void safe_python_command(const std::string &python_command); 

/*  Is this a repeat of something?  I don't know. */
void run_generic_script(const std::vector<std::string> &cmd_strings);


/*  ----------------------------------------------------------------------- */
/*               Atom info                                                  */
/*  ----------------------------------------------------------------------- */

/*! \brief output atom info in a scheme list for use in scripting:

in this format (list occ temp-factor element x y z).  Return empty
list if atom not found. */
#ifdef USE_GUILE
const char *atom_info_string(int imol, const char *chain_id, int resno,
			     const char *ins_code, const char *atname,
			     const char *altconf);

//! \brief
// Return a list of atom info for each atom in the specified residue:
// 
// output is like this:
// (list
//    (list (list atom-name alt-conf)
//          (list occ temp-fact element)
//          (list x y z)))
// 
SCM residue_info(int imol, const char* chain_id, int resno, const char *ins_code);

// And going the other way, given an s-expression, update
// molecule_number by the given molecule.  Clear what's currently
// there first though.
//
int clear_and_update_molecule(int molecule_number, SCM molecule_expression);
// return a molecule number, -1 on error
int add_molecule(SCM molecule_expression, const char *name);

//! \brief 
// Return a list of (list imol chain-id resno ins-code atom-name
// alt-conf) for atom that is closest to the screen centre.  If there
// are multiple models with the same coordinates at the screen centre,
// return the attributes of the atom in the highest number molecule
// number.
//
// return #f if no active residue
// 
SCM active_residue();

#endif	/* USE_GUILE */


/*  ----------------------------------------------------------------------- */
/*                  monomer lib                                             */
/*  ----------------------------------------------------------------------- */
std::vector<std::pair<std::string, std::string> > monomer_lib_3_letter_codes_matching(const std::string &search_string, short int allow_minimal_descriptions_flag);

void on_monomer_lib_search_results_button_press (GtkButton *button, gpointer user_data);

/*  ----------------------------------------------------------------------- */
/*                  mutate                                                  */
/*  ----------------------------------------------------------------------- */
int mutate_internal(int ires, const char *chain_id,
		    int imol, std::string &target_res_type);
/* a function for multimutate to make a backup and set
   have_unsaved_changes_flag themselves */



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
SCM rtop_to_scm(const clipper::RTop_orth &rtop);
#endif	/* USE_GUILE */


#endif // CC_INTERFACE_HH
