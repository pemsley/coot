/* src/c-interface-nucleotides.cc
 * 
 * Copyright 2008 The University of Oxford
 * Author: Paul Emsley
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
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"


#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"

#include "graphics-info.h"

#include "c-interface.h"
#include "guile-fixups.h"

#ifdef USE_GUILE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#endif // USE_GUILE


// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
// 20100813: Python.h needs to come before to stop"_POSIX_C_SOURCE" redefined problems
//
#ifdef USE_PYTHON
#if (PY_MINOR_VERSION > 4)
// no fixup needed
#else
#define Py_ssize_t int
#endif
#endif // USE_PYTHON

#include "cc-interface.hh"
#include "ligand/ideal-rna.hh"
#include "ligand/base-pairing.hh"

#ifdef USE_GUILE
SCM pucker_info_scm(int imol, SCM residue_spec_scm, int do_pukka_pucker_check) {

   // non-nucleic acid really
   std::vector<std::string> protein_residue_names = { "GLY", "ALA", "CYS", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU",
                                                      "MET", "MSE", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP",
                                                      "TYR", "MG", "ZN"};
   auto is_a_protein_residue = [&protein_residue_names] (mmdb::Residue *residue_p) {
                                  std::string res_name(residue_p->GetResName());
                                  return (std::find(protein_residue_names.begin(),
                                                    protein_residue_names.end(),
                                                    res_name) != protein_residue_names.end());
                               };
   std::string altconf = "";
   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_spec_from_scm(residue_spec_scm);
      mmdb::Residue *res_p = graphics_info_t::molecules[imol].get_residue(residue_spec);
      if (is_a_protein_residue(res_p)) {
      } else {
         if (res_p) {
            try {
               coot::pucker_analysis_info_t pi(res_p, altconf);
               mmdb::Residue *following_res =
                  graphics_info_t::molecules[imol].get_following_residue(residue_spec);
               if (do_pukka_pucker_check) {
                  if (following_res) {
                     // 		  std::cout << "   DEBUG:: " << coot::residue_spec_t(following_res)
                     // 			    << " follows " << residue_spec << std::endl;
                     try {
                        double phosphate_d = pi.phosphate_distance(following_res);
                        r = SCM_EOL;
                        r = scm_cons(scm_from_double(pi.plane_distortion), r);
                        r = scm_cons(scm_from_double(pi.out_of_plane_distance), r);
                        r = scm_cons(scm_from_locale_string(pi.puckered_atom().c_str()), r);
                        r = scm_cons(scm_from_double(phosphate_d), r);

                        // double dist_crit = xxx
                        // If C2', phosphate oop dist should be > dist_crit
                        // If C3', phosphate oop dist should be < dist_crit

                     }
                     catch (const std::runtime_error &phos_mess) {
                        std::cout << " Fail in Pucker analysis for "
                                  << coot::residue_spec_t(res_p) << " "
                                  << phos_mess.what() << std::endl;
                     }

                  } else {
                     r = SCM_EOL;
                  }
               } else {
                  // no pucker check
                  r = SCM_EOL;
                  r = scm_cons(scm_from_double(pi.plane_distortion), r);
                  r = scm_cons(scm_from_double(pi.out_of_plane_distance), r);
                  r = scm_cons(scm_from_locale_string(pi.puckered_atom().c_str()), r);
                  if (following_res) {
                     try {
                        double phosphate_d = pi.phosphate_distance(following_res);
                        r = scm_cons(scm_from_double(phosphate_d), r);
                     }
                     catch (const std::runtime_error &phos_mess) { }
                  }
               }
            }
            catch (const std::runtime_error &mess) {
               std::cout << " failed to find pucker for " << residue_spec << " "
                         << mess.what() << std::endl;
            }
         }
      }
   }
   return r;
}
#endif /* USE_GUILE */

#ifdef USE_PYTHON
PyObject *pucker_info_py(int imol, PyObject *residue_spec_py, int do_pukka_pucker_check) {

   // non-nucleic acid really
   std::vector<std::string> protein_residue_names = { "GLY", "ALA", "CYS", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU",
                                                      "MET", "MSE", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP",
                                                      "TYR", "MG", "ZN"};

   auto is_a_protein_residue = [&protein_residue_names] (mmdb::Residue *residue_p) {
                                  std::string res_name(residue_p->GetResName());
                                  return (std::find(protein_residue_names.begin(),
                                                    protein_residue_names.end(),
                                                    res_name) != protein_residue_names.end());
                               };

   std::string altconf = "";
   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_spec_from_py(residue_spec_py);
      mmdb::Residue *res_p = graphics_info_t::molecules[imol].get_residue(residue_spec);
      if (is_a_protein_residue(res_p)) {
      } else {
         if (res_p) {
            try {
               coot::pucker_analysis_info_t pi(res_p, altconf);
               mmdb::Residue *following_res =
                  graphics_info_t::molecules[imol].get_following_residue(residue_spec);
               if (do_pukka_pucker_check) {
                  if (following_res) {
                     // 		  std::cout << "   DEBUG:: " << coot::residue_spec_t(following_res)
                     // 			    << " follows " << residue_spec << std::endl;
                     try {
                        double phosphate_d = pi.phosphate_distance(following_res);
                        r = PyList_New(4);
                        PyList_SetItem(r, 0, PyFloat_FromDouble(phosphate_d));
                        PyList_SetItem(r, 1, myPyString_FromString(pi.puckered_atom().c_str()));
                        PyList_SetItem(r, 2, PyFloat_FromDouble(pi.out_of_plane_distance));
                        PyList_SetItem(r, 3, PyFloat_FromDouble(pi.plane_distortion));

                        // double dist_crit = 1.2;
                        // If C2', phosphate oop dist should be > dist_crit
                        // If C3', phosphate oop dist should be < dist_crit

                     }
                     catch (const std::runtime_error &phos_mess) {
                        std::cout << " Failed in pucker analysis for "
                                  << coot::residue_spec_t(following_res) << " "
                                  << phos_mess.what() << std::endl;
                     }

                  } else {
                     r = PyList_New(0);
                  }
               } else {
                  // no pucker check
                  r = PyList_New(3);
                  PyList_SetItem(r, 0, myPyString_FromString(pi.puckered_atom().c_str()));
                  PyList_SetItem(r, 1, PyFloat_FromDouble(pi.out_of_plane_distance));
                  PyList_SetItem(r, 2, PyFloat_FromDouble(pi.plane_distortion));
                  if (following_res) {
                     try {
                        double phosphate_d = pi.phosphate_distance(following_res);
                        PyList_Insert(r, 0, PyFloat_FromDouble(phosphate_d));
                     }
                     catch (const std::runtime_error &phos_mess) { }
                  }
               }
            }
            catch (const std::runtime_error &mess) {
               std::cout << " failed to find pucker for " << residue_spec << " "
                         << mess.what() << std::endl;
            }
         }
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif /* USE_PYTHON */



/*  \brief create a molecule of idea nucleotides 

use the given sequence (single letter code)

RNA_or_DNA is either "RNA" or "DNA"

form is either "A" or "B"

@return the new molecule number or -1 if a problem */
int ideal_nucleic_acid(const char *RNA_or_DNA, const char *form,
		       short int single_stranded_flag,
		       const char *sequence) {

   int istat = -1; 
   short int do_rna_flag = -1;
   short int form_flag = -1;

   float here_x = graphics_info_t::RotationCentre_x();
   float here_y = graphics_info_t::RotationCentre_y();
   float here_z = graphics_info_t::RotationCentre_z();

   std::string RNA_or_DNA_str(RNA_or_DNA);
   std::string form_str(form);

   if (RNA_or_DNA_str == "RNA")
      do_rna_flag = 1;
   if (RNA_or_DNA_str == "DNA")
      do_rna_flag = 0;

   if (form_str == "A")
      form_flag = 1;
   
   if (form_str == "B")
      form_flag = 1;

   if (! (form_flag > 0)) {
      std::cout << "Problem in nucleic acid form, use only either \"A\" or \"B\"."
		<< std::endl;
   } else {
      if (! (do_rna_flag >= 0)) {
	 std::cout << "Problem in nucleic acid type, use only either \"RNA\" or \"DNA\"."
		   << "You said: \"" << RNA_or_DNA << "\"" << std::endl;
      } else {
	 // proceed, input is good

	 std::string down_sequence(coot::util::downcase(sequence));
	 if (graphics_info_t::standard_residues_asc.read_success) {
	    coot::ideal_rna ir(RNA_or_DNA_str, form_str, single_stranded_flag,
			       down_sequence,
			       graphics_info_t::standard_residues_asc.mol);
	    ir.use_v3_names();
	    mmdb::Manager *mol = ir.make_molecule();

	    if (mol) {
	       std::pair<bool, clipper::Coord_orth> cm = coot::centre_of_molecule(mol);
	       graphics_info_t g;
	       if (cm.first) { 
		  clipper::Coord_orth mc = cm.second; // just for alias
		  int imol = graphics_info_t::create_molecule();
		  istat = imol;
		  std::string label = form_str;
		  label += "-form-";
		  label += RNA_or_DNA_str;
		  atom_selection_container_t asc = make_asc(mol);
		  graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), label, 1);
		  graphics_info_t::molecules[imol].translate_by(here_x-mc.x(), here_y-mc.y(), here_z-mc.z());
		  graphics_draw();
		  if (graphics_info_t::go_to_atom_window) {
		     g.update_go_to_atom_window_on_new_mol();
		     g.update_go_to_atom_window_on_changed_mol(imol);
		  }
	       } else {
		  std::cout << "WARNING:: ideal_nucleic_acid()::something bad happened "
			    << std::endl;
	       } 
	    }
	 } else {
	    std::string s("WARNING:: Can't proceed with Idea RNA - no standard residues!");
	    std::cout << s << std::endl;
	    graphics_info_t g;
	    g.add_status_bar_text(s);
	 } 
      }
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("ideal-nucleic-acid");
   command_strings.push_back(single_quote(RNA_or_DNA_str));
   command_strings.push_back(single_quote(form_str));
   command_strings.push_back(coot::util::int_to_string(single_stranded_flag));
   command_strings.push_back(single_quote(sequence));
   add_to_history(command_strings);

   return istat;
}


/* Return a molecule that contains a residue that is the WC pair
   partner of the clicked/picked/selected residue */
int watson_crick_pair(int imol, const char *chain_id, int resno) {

   int imol_return = -1;

   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *res = graphics_info_t::molecules[imol].get_residue(chain_id, resno, "");
      if (!res) {
	 std::cout << "Residue not found in " << imol << " " << chain_id << " " << resno
		   << std::endl;
      } else { 
	 mmdb::Residue *res_wc =
	    coot::watson_crick_partner(res, graphics_info_t::standard_residues_asc.mol);
	 if (res_wc) {
	    mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(res_wc);
	    if (mol) {
	       graphics_info_t g;
	       int imol_new = g.create_molecule();
	       atom_selection_container_t asc_wc = make_asc(mol);
	       g.molecules[imol_new].install_model(imol_new, asc_wc, g.Geom_p(), "WC partner", 1);
	       graphics_draw();
	    }
	 }
      } 
   } 

   return imol_return; 
} 


/* not for user level */
void setup_base_pairing(int state) {

   graphics_info_t g;
   if (state) { 
      g.in_base_paring_define = 1;
      pick_cursor_maybe();
   } else { 
      g.in_base_paring_define = 0;
      normal_cursor();
   } 

} 

/*! \brief add base pairs for the given residue range, modify molecule imol by creating a new chain */
int
watson_crick_pair_for_residue_range(int imol, const char * chain_id, int resno_start, int resno_end) {
   int status = 0;
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      status = g.molecules[imol].watson_crick_pair_for_residue_range(chain_id,
								     resno_start, resno_end,
								     g.standard_residues_asc.mol);
      graphics_draw();
   } 
   return status;
}
