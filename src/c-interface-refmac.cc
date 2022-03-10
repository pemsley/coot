/* src/c-interface-build-gui.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007,2008, 2009 by Bernhard Lohkamp
 * Copyright 2008 by Kevin Cowtan
 * Copyright 2007, 2008, 2009 The University of Oxford
 * Copyright 2015, 2016 by Medical Research Council
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


// 20220309-PE  there is nothing here now - this gui never worked properly and has been replaced.
// I don't want to fix up the widget references for gtk3 for something that is not used.

#include <iostream>
#include "cc-interface.hh"
#include "graphics-info.h"
#include "cc-interface-scripting.hh"
#include "c-interface.h" // for is_valid_map_molecule()
#include "python-3-interface.hh" // for myPyString_FromString()

void set_refmac_counter(int imol, int refmac_count) {

   graphics_info_t g;
   if (imol < g.n_molecules()) {
      g.molecules[imol].set_refmac_counter(refmac_count);
      std::cout << "INFO:: refmac counter of molecule number " << imol
               << " incremented to " << refmac_count << std::endl;
   } else {
      std::cout << "WARNING:: refmac counter of molecule number " << imol
               << " not incremented to " << refmac_count << std::endl;
   }
   std::string cmd = "set-refmac-counter";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(refmac_count);
   add_to_history_typed(cmd, args);
}

std::string refmac_name(int imol) {

   graphics_info_t g;
   std::string cmd = "refmac-name";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return g.molecules[imol].Refmac_in_name();
}

#ifdef USE_GUILE
// return a list of refmac parameters.  Used so that we can test that
// the save state of a refmac map works correctly.
SCM refmac_parameters_scm(int imol) {

   SCM r = SCM_EOL;
   if (is_valid_map_molecule(imol)) {
      std::vector<coot::atom_attribute_setting_help_t>
	 refmac_params = graphics_info_t::molecules[imol].get_refmac_params();
      if (refmac_params.size() > 0) {
	 // values have to go in in reverse order, as usual.
	 for (int i=(int(refmac_params.size())-1); i>=0; i--) {
	    if (refmac_params[i].type == coot::atom_attribute_setting_help_t::IS_STRING)
	       r = scm_cons(scm_from_locale_string(refmac_params[i].s.c_str()) ,r);
	    if (refmac_params[i].type == coot::atom_attribute_setting_help_t::IS_FLOAT)
	       r = scm_cons(scm_from_double(refmac_params[i].val) ,r);
	    if (refmac_params[i].type == coot::atom_attribute_setting_help_t::IS_INT)
	       r = scm_cons(scm_from_int(refmac_params[i].i) ,r);
	 }
      }
   }
   return r;
}

#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *refmac_parameters_py(int imol) {

   PyObject *r = PyList_New(0);
   if (is_valid_map_molecule(imol)) {
      std::vector<coot::atom_attribute_setting_help_t>
	 refmac_params = graphics_info_t::molecules[imol].get_refmac_params();
      if (refmac_params.size() > 0) {
	 // values have dont have to go in in reverse order.
	for (unsigned int i=0; i<refmac_params.size(); i++) {
	    if (refmac_params[i].type == coot::atom_attribute_setting_help_t::IS_INT)
	      PyList_Append(r, PyLong_FromLong(refmac_params[i].i));
	    if (refmac_params[i].type == coot::atom_attribute_setting_help_t::IS_FLOAT)
	      PyList_Append(r, PyFloat_FromDouble(refmac_params[i].val));
	    if (refmac_params[i].type == coot::atom_attribute_setting_help_t::IS_STRING)
	      PyList_Append(r, myPyString_FromString(refmac_params[i].s.c_str()));
	 }
      }
   }
   return r;
}
#endif	/* USE_PYTHON */


//       int slen = mtz_in_filename.length(); c
//       if (slen > 4) {
// 	 mtz_out_filename = mtz_in_filename.substr(0,slen - 4) + "-refmac-";
// 	 mtz_out_filename += g.int_to_string(g.molecules[imol_coords].Refmac_count());
// 	 mtz_out_filename += ".mtz";
//       } else {
// 	 mtz_out_filename = "post-refmac";
// 	 mtz_out_filename += g.int_to_string(g.molecules[imol_coords].Refmac_count());
// 	 mtz_out_filename += ".mtz";
//       }

// If ccp4i_project_dir is "", then carry on and put the log file in
// this directory.  If not, put it in the appropriate project dir. The
// pdb_in etc filename are manipulated in the calling routine.
//
// if swap_map_colours_post_refmac_flag is not 1 then imol_refmac_map is ignored.
//
// make_molecules_flag is 0 or 1: defining if run-refmac-by-filename
// function should create molecules (it should *not* create molecules
// if this is called in a sub-thread (because that will try to update
// the graphics from the subthread and a crash will result).
//
void
execute_refmac_real(std::string pdb_in_filename,
		    std::string pdb_out_filename,
		    std::string mtz_in_filename,
		    std::string mtz_out_filename,
		    std::string cif_lib_filename,
		    std::string fobs_col_name,
		    std::string sigfobs_col_name,
		    std::string r_free_col_name,
		    short int have_sensible_free_r_flag,
		    short int make_molecules_flag,
		    std::string refmac_count_str,
		    int swap_map_colours_post_refmac_flag,
		    int imol_refmac_map,
		    int diff_map_flag,
		    int phase_combine_flag,
		    std::string phib_string,
		    std::string fom_string,
		    std::string ccp4i_project_dir) {


   std::vector<std::string> cmds;

   cmds.push_back(std::string("run-refmac-by-filename"));
// BL says:: again debackslashing
   cmds.push_back(single_quote(coot::util::intelligent_debackslash(pdb_in_filename)));
   cmds.push_back(single_quote(coot::util::intelligent_debackslash(pdb_out_filename)));
   cmds.push_back(single_quote(coot::util::intelligent_debackslash(mtz_in_filename)));
   cmds.push_back(single_quote(coot::util::intelligent_debackslash(mtz_out_filename)));
   cmds.push_back(single_quote(coot::util::intelligent_debackslash(cif_lib_filename)));
   cmds.push_back(refmac_count_str);
   cmds.push_back(graphics_info_t::int_to_string(swap_map_colours_post_refmac_flag));
   cmds.push_back(graphics_info_t::int_to_string(imol_refmac_map));
   cmds.push_back(graphics_info_t::int_to_string(diff_map_flag));
   cmds.push_back(graphics_info_t::int_to_string(phase_combine_flag));

   std::string phase_combine_cmd;
   if (phase_combine_flag > 0 && phase_combine_flag < 3) {
#ifdef USE_GUILE
      phase_combine_cmd += "(cons ";
      phase_combine_cmd += single_quote(phib_string);
      phase_combine_cmd += " ";
      phase_combine_cmd += single_quote(fom_string);
      phase_combine_cmd += ")";
#else
#ifdef USE_PYTHON
      phase_combine_cmd += "[\'";
      phase_combine_cmd += phib_string;
      phase_combine_cmd += "\', ";
      phase_combine_cmd += single_quote(fom_string);
      phase_combine_cmd += "]";
#endif // USE_PYTHON
#endif // USE_GUILE
   } else {
      phase_combine_cmd += single_quote("dummy");
   }
   cmds.push_back(phase_combine_cmd);

   cmds.push_back(graphics_info_t::int_to_string(graphics_info_t::refmac_ncycles));
   cmds.push_back(graphics_info_t::int_to_string(make_molecules_flag));
   cmds.push_back(single_quote(coot::util::intelligent_debackslash(ccp4i_project_dir)));
   if (phase_combine_flag == 3 && fobs_col_name != "") {
     cmds.push_back(fobs_col_name);
     cmds.push_back(sigfobs_col_name);
   } else {
     cmds.push_back(single_quote(fobs_col_name));
     cmds.push_back(single_quote(sigfobs_col_name));
   }
   std::cout << "DEBUG:: in execute_refmac_real() ccp4i_project_dir :"
	     << single_quote(coot::util::intelligent_debackslash(ccp4i_project_dir))
	     << ":" << std::endl;

   if (have_sensible_free_r_flag) {
      cmds.push_back(single_quote(r_free_col_name));
   }

   graphics_info_t g;
   short int ilang = coot::STATE_SCM;
   std::string cmd;

#ifdef USE_PYTHON
#ifndef USE_GUILE
   ilang = coot::STATE_PYTHON;
#endif
#endif
   if (ilang == coot::STATE_PYTHON) {
      cmd = g.state_command(cmds, ilang);
#ifdef USE_PYTHON
      safe_python_command(cmd);
#endif
   } else {
      cmd = g.state_command(cmds, ilang);
      safe_scheme_command(cmd);
   }
}


#ifdef USE_GUILE
SCM get_refmac_sad_atom_info_scm() {

   return SCM_EOL;
}
#endif
