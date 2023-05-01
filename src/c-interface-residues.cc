/* src/c-interface-residues.cc
 * 
 * Copyright 2012 by The University of Oxford
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

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"

#ifdef USE_GUILE
#include <cstddef> // define std::ptrdiff_t for clang
#include <libguile.h>
#endif

#include <vector>
#include "named-rotamer-score.hh"

#include "cc-interface.hh"
#include "graphics-info.h"

#include "c-interface.h" // for is_valid_model_molecule()

#include "c-interface-python.hh"

std::vector<coot::named_rotamer_score>
score_rotamers(int imol, 
	       const char *chain_id, 
	       int res_no, 
	       const char *ins_code, 
	       const char *alt_conf, 
	       int imol_map, 
	       int clash_flag,
	       float lowest_probability) {

   std::vector<coot::named_rotamer_score> v;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
	 graphics_info_t g;
	 const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
	 v = graphics_info_t::molecules[imol].score_rotamers(chain_id, res_no, ins_code,
							     alt_conf,
							     clash_flag, lowest_probability,
							     xmap, *g.Geom_p());
      }
   }
   return v; 
} 
 

#ifdef USE_GUILE
SCM score_rotamers_scm(int imol, 
		       const char *chain_id, 
		       int res_no, 
		       const char *ins_code, 
		       const char *alt_conf, 
		       int imol_map, 
		       int clash_flag, 
		       float lowest_probability) {

   SCM r = SCM_EOL;
   std::vector<coot::named_rotamer_score> v =
      score_rotamers(imol, chain_id, res_no, ins_code, alt_conf,
		     imol_map, clash_flag, lowest_probability);
   for (unsigned int i=0; i<v.size(); i++) {
      SCM name_scm  = scm_from_locale_string(v[i].name.c_str());
      SCM prob_scm  = scm_from_double(v[i].rotamer_probability_score);
      SCM fit_scm   = scm_from_double(v[i].density_fit_score);
      SCM clash_scm = scm_from_double(v[i].clash_score);
      SCM atom_list_scm = SCM_EOL;
      for (unsigned int iat=0; iat<v[i].density_score_for_atoms.size(); iat++) {
	 SCM p1 = scm_from_locale_string(v[i].density_score_for_atoms[iat].first.c_str());
	 SCM p2 = scm_from_double(v[i].density_score_for_atoms[iat].second);
	 SCM atom_item = scm_list_2(p1,p2);
	 atom_list_scm = scm_cons(atom_item, atom_list_scm);
      }
      atom_list_scm = scm_reverse(atom_list_scm);
      SCM item = scm_list_5(name_scm, prob_scm, atom_list_scm, fit_scm, clash_scm);
      
      r = scm_cons(item, r);
   }
   r = scm_reverse(r);
   return r;
} 
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
			    float lowest_probability) {
   
   std::vector<coot::named_rotamer_score> v =
      score_rotamers(imol, chain_id, res_no, ins_code, alt_conf,
		     imol_map, clash_flag, lowest_probability);
   PyObject *r = PyList_New(v.size());
   for (unsigned int i=0; i<v.size(); i++) { 
      PyObject *item = PyList_New(5);
      PyObject *name_py  = myPyString_FromString(v[i].name.c_str());
      PyObject *prob_py  = PyFloat_FromDouble(v[i].rotamer_probability_score);;
      PyObject *fit_py   = PyFloat_FromDouble(v[i].density_fit_score);;
      PyObject *clash_py = PyFloat_FromDouble(v[i].clash_score);;
      PyObject *atom_list_py = PyList_New(v[i].density_score_for_atoms.size());
      for (unsigned int iat=0; iat<v[i].density_score_for_atoms.size(); iat++) {
	 PyObject *atom_item = PyList_New(2);
	 PyObject *p0 = myPyString_FromString(v[i].density_score_for_atoms[iat].first.c_str());
	 PyObject *p1 = PyFloat_FromDouble(v[i].density_score_for_atoms[iat].second);
	 PyList_SetItem(atom_item, 0, p0);
	 PyList_SetItem(atom_item, 1, p1);
	 PyList_SetItem(atom_list_py, iat, atom_item);
      }
      PyList_SetItem(item, 0, name_py);
      PyList_SetItem(item, 1, prob_py);
      PyList_SetItem(item, 2, fit_py);
      PyList_SetItem(item, 3, atom_list_py);
      PyList_SetItem(item, 4, clash_py);
      PyList_SetItem(r, i, item);
   }
   return r;
} 
#endif


void
print_glyco_tree(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Residue *r = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	 
      if (r) {

	 std::vector<std::string> types_with_no_dictionary =
	    g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
	 for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	    g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

	 coot::glyco_tree_t t(r, mol, g.Geom_p());

      } 
   } 
} 


/* ------------------------------------------------------------------------- */
/*                      interesting positions list                           */
/* ------------------------------------------------------------------------- */

#ifdef USE_GUILE
// pass a list of (position,string) pairs
void register_interesting_positions_list_scm(SCM pos_list) {

   std::vector<std::pair<clipper::Coord_orth, std::string> > v;
   graphics_info_t g;

   if (scm_is_true(scm_list_p(pos_list))) {
      unsigned int pos_length = scm_to_int(scm_length(pos_list));
      for (unsigned int i=0; i<pos_length; i++) { 
	 SCM item = scm_list_ref(pos_list, scm_from_int(i));
	 if (scm_is_true(scm_list_p(item))) {
	    // pos, label pair
	    unsigned int item_length = scm_to_int(scm_length(item));
	    if (item_length == 2) {
	       SCM item_item_0 = scm_list_ref(item, 0);
	       SCM item_item_1 = scm_list_ref(item, 0);

	       if (scm_is_true(scm_list_p(item_item_0))) {
		  unsigned int l_p = scm_to_int(scm_length(item_item_0));
		  if (l_p == 3) {
		     SCM x = scm_list_ref(item_item_0, scm_from_int(0));
		     SCM y = scm_list_ref(item_item_0, scm_from_int(1));
		     SCM z = scm_list_ref(item_item_0, scm_from_int(2));
		     if (scm_number_p(x)) { 
			if (scm_number_p(y)) { 
			   if (scm_number_p(z)) {
			      clipper::Coord_orth pos(scm_to_double(x),
						      scm_to_double(y),
						      scm_to_double(z));

			      if (scm_is_true(scm_string_p(item_item_1))) { 
				 std::string s = scm_to_locale_string(item_item_1);
				 std::pair<clipper::Coord_orth, std::string> p(pos,s);
				 v.push_back(p);
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   } 
   
   g.register_user_defined_interesting_positions(v);
}
#endif // USE_GUILE

#ifdef USE_PYTHON
void register_interesting_positions_list_py(PyObject *pos_list) {

   std::vector<std::pair<clipper::Coord_orth, std::string> > v;

   if (PyList_Check(pos_list)) {
      unsigned int pos_length = PyObject_Length(pos_list);
      for (unsigned int i=0; i<pos_length; i++) {
	 PyObject *item = PyList_GetItem(pos_list, i);
	 if (PyList_Check(item)) {
	    unsigned int l_item = PyObject_Length(item);
	    if (l_item == 2) {
	       PyObject *item_item_0 = PyList_GetItem(item, 0);
	       PyObject *item_item_1 = PyList_GetItem(item, 1);

	       if (PyUnicode_Check(item_item_1)) {
		  if (PyList_Check(item_item_0)) {

		     unsigned int l_item_item = PyObject_Length(item_item_0);
		     if (l_item_item == 3) {
			PyObject *x = PyList_GetItem(item_item_0, 0);
			PyObject *y = PyList_GetItem(item_item_0, 1);
			PyObject *z = PyList_GetItem(item_item_0, 2);

			if (PyFloat_Check(x)) { 
			   if (PyFloat_Check(y)) { 
			      if (PyFloat_Check(z)) {

				 clipper::Coord_orth pos(PyFloat_AsDouble(x),
							 PyFloat_AsDouble(y),
							 PyFloat_AsDouble(z));
				 std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(item_item_1));
				 std::pair<clipper::Coord_orth, std::string> p(pos,s);
				 v.push_back(p);
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   graphics_info_t g;
   g.register_user_defined_interesting_positions(v);
}
#endif // USE_PYTHON 


#ifdef USE_GUILE
SCM glyco_tree_scm(int imol, SCM active_residue_scm) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {

      coot::residue_spec_t residue_spec = residue_spec_from_scm(active_residue_scm);
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(residue_spec);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      std::vector<std::string> types_with_no_dictionary =
	 g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
      for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	 g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);
      coot::glyco_tree_t t(residue_p, mol, g.Geom_p());
   }
   return r;
}
#endif

#ifdef USE_GUILE
SCM glyco_tree_residues_scm(int imol, SCM active_residue_scm) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {

      r = SCM_EOL;
      coot::residue_spec_t residue_spec = residue_spec_from_scm(active_residue_scm);
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(residue_spec);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      std::vector<std::string> types_with_no_dictionary =
	 g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
      for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	 g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);
      coot::glyco_tree_t t(residue_p, mol, g.Geom_p());
      std::vector<mmdb::Residue *> v_residues = t.residues(residue_spec);
      for (unsigned int i=0; i<v_residues.size(); i++) {
	 if (v_residues[i]) {
	    coot::residue_spec_t spec(v_residues[i]);
	    SCM spec_scm = residue_spec_to_scm(spec);
	    r = scm_cons(spec_scm, r);
	 }
      }
      r = scm_reverse(r);
   }

   return r;
}
#endif


#include "cc-interface.hh" // needed?

// need python version of this and glyco_tree_compare_trees_scm and
// glyco_tree_matched_residue_pairs_scm (although those are only
// useful for the results for the paper).
//
// BL says:: there is a python version below now...
#ifdef USE_GUILE
SCM glyco_tree_residue_id_scm(int imol, SCM residue_spec_scm) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_spec_from_scm(residue_spec_scm);
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(residue_spec);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      std::vector<std::string> types_with_no_dictionary =
	 g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
      for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	 g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);
      coot::glyco_tree_t t(residue_p, mol, g.Geom_p());
      coot::glyco_tree_t::residue_id_t id = t.get_id(residue_p);
      if (false)
	 std::cout << "got id " << id.level << " " << id.prime_arm_flag << " "
		   << id.res_type << std::endl;
      if (! id.res_type.empty()) {
	 SCM parent_spec_scm = residue_spec_to_scm(id.parent_res_spec);
	 SCM prime_flag_sym = scm_string_to_symbol(scm_from_locale_string("unset"));
	 if (id.prime_arm_flag == coot::glyco_tree_t::residue_id_t::PRIME)
	    prime_flag_sym = scm_string_to_symbol(scm_from_locale_string("prime"));
	 if (id.prime_arm_flag == coot::glyco_tree_t::residue_id_t::NON_PRIME)
	    prime_flag_sym = scm_string_to_symbol(scm_from_locale_string("non-prime"));

	 r = SCM_LIST6(scm_from_int(id.level),
		       prime_flag_sym,
		       scm_from_locale_string(id.res_type.c_str()),
		       scm_from_locale_string(id.link_type.c_str()),
		       scm_from_locale_string(id.parent_res_type.c_str()),
		       parent_spec_scm);
      }
   }
   return r;
}
#endif

#ifdef USE_GUILE
SCM glyco_tree_compare_trees_scm(int imol_1, SCM res_spec_1_scm, int imol_2, SCM res_spec_2_scm) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {
	 graphics_info_t g;

	 coot::residue_spec_t residue_spec_1 = residue_spec_from_scm(res_spec_1_scm);
	 mmdb::Residue *residue_1_p = g.molecules[imol_1].get_residue(residue_spec_1);
	 mmdb::Manager *mol_1 = g.molecules[imol_1].atom_sel.mol;
	 std::vector<std::string> types_with_no_dictionary =
	    g.molecules[imol_1].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
	 for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	    g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

	 coot::residue_spec_t residue_spec_2 = residue_spec_from_scm(res_spec_2_scm);
	 mmdb::Residue *residue_2_p = g.molecules[imol_2].get_residue(residue_spec_2);
	 mmdb::Manager *mol_2 = g.molecules[imol_2].atom_sel.mol;
	 types_with_no_dictionary = g.molecules[imol_2].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
	 for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	    g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

	 coot::glyco_tree_t t_1(residue_1_p, mol_1, g.Geom_p());
	 coot::glyco_tree_t t_2(residue_2_p, mol_2, g.Geom_p());

	 bool match = t_1.compare_trees(t_2.get_glyco_tree());

	 if (match)
	    r = SCM_BOOL_T;
      }
   }
   return r;
}
#endif

#ifdef USE_GUILE
SCM glyco_tree_matched_residue_pairs_scm(int imol_1, SCM res_spec_1_scm, int imol_2, SCM res_spec_2_scm) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {
	 graphics_info_t g;

	 coot::residue_spec_t residue_spec_1 = residue_spec_from_scm(res_spec_1_scm);
	 mmdb::Residue *residue_1_p = g.molecules[imol_1].get_residue(residue_spec_1);
	 mmdb::Manager *mol_1 = g.molecules[imol_1].atom_sel.mol;
	 std::vector<std::string> types_with_no_dictionary =
	    g.molecules[imol_1].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
	 for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	    g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

	 coot::residue_spec_t residue_spec_2 = residue_spec_from_scm(res_spec_2_scm);
	 mmdb::Residue *residue_2_p = g.molecules[imol_2].get_residue(residue_spec_2);
	 mmdb::Manager *mol_2 = g.molecules[imol_2].atom_sel.mol;
	 types_with_no_dictionary = g.molecules[imol_2].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
	 for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	    g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

	 coot::glyco_tree_t t_1(residue_1_p, mol_1, g.Geom_p());
	 coot::glyco_tree_t t_2(residue_2_p, mol_2, g.Geom_p());

	 std::vector<std::pair<coot::residue_spec_t, coot::residue_spec_t> > pv = t_1.matched_pairs(t_2.get_glyco_tree());
	 if (! pv.empty()) {
	    SCM rt = SCM_EOL;
	    for (unsigned int i=0; i<pv.size(); i++) {
	       SCM r_1_scm = residue_spec_to_scm(pv[i].first);
	       SCM r_2_scm = residue_spec_to_scm(pv[i].second);
	       SCM list_scm = scm_list_2(r_1_scm, r_2_scm);
	       rt = scm_cons(list_scm, rt);
	    }
	    r = rt;
	 }
      }
   }
   return r;
}
#endif



#ifdef USE_PYTHON
PyObject *glyco_tree_py(int imol, PyObject *active_residue_py) {

   // incomplete, doesn't work.
   // BL says:: should now

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {

      coot::residue_spec_t residue_spec = residue_spec_from_py(active_residue_py);
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(residue_spec);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      std::vector<std::string> types_with_no_dictionary =
            g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
      for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
         g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);
      coot::glyco_tree_t t(residue_p, mol, g.Geom_p());

   }

   if (PyBool_Check(r))
      Py_INCREF(r);

   return r;

}
#endif /* PYTHON */

#ifdef USE_PYTHON
PyObject *glyco_tree_residues_py(int imol, PyObject *active_residue_py) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {

      coot::residue_spec_t residue_spec = residue_spec_from_py(active_residue_py);
//      std::cout << "..... active residue spec: " << residue_spec << std::endl;
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(residue_spec);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      std::vector<std::string> types_with_no_dictionary =
            g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
      for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
         g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);
      coot::glyco_tree_t t(residue_p, mol, g.Geom_p());
      std::vector<mmdb::Residue *> v_residues = t.residues(residue_spec);
      // std::vector<coot::residue_spec_t> v(v_residues.size()); delete me
      r = PyList_New(v_residues.size());
      for (unsigned int i=0; i<v_residues.size(); i++) {
         // std::cout << "     " << i << " " << coot::residue_spec_t(v_residues[i]) << std::endl;
         coot::residue_spec_t spec(v_residues[i]);
         PyList_SetItem(r, i, residue_spec_to_py(spec));
      }
   }

   if (PyBool_Check(r))
      Py_INCREF(r);

   return r;

}
#endif /* PYTHON */

#ifdef USE_PYTHON
PyObject *glyco_tree_residue_id_py(int imol, PyObject *residue_spec_py) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_spec_from_py(residue_spec_py);
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(residue_spec);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      std::vector<std::string> types_with_no_dictionary =
            g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
      for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
         g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);
      coot::glyco_tree_t t(residue_p, mol, g.Geom_p());
      coot::glyco_tree_t::residue_id_t id = t.get_id(residue_p);
      std::cout << "got id " << id.level << " " << id.prime_arm_flag << " "
      << id.res_type << std::endl;
      if (! id.res_type.empty()) {
         PyObject *parent_spec_py = residue_spec_to_py(id.parent_res_spec);
         PyObject *prime_flag_sym = myPyString_FromString("unset");
         if (id.prime_arm_flag == coot::glyco_tree_t::residue_id_t::PRIME)
            prime_flag_sym = myPyString_FromString("prime");
         if (id.prime_arm_flag == coot::glyco_tree_t::residue_id_t::NON_PRIME)
            prime_flag_sym = myPyString_FromString("non-prime");

         PyObject *level = PyLong_FromLong(id.level);
         PyObject *res_type = myPyString_FromString(id.res_type.c_str());
         PyObject *link_type = myPyString_FromString(id.link_type.c_str());
         PyObject *parent = myPyString_FromString(id.parent_res_type.c_str());
         r = PyList_New(6);
         PyList_SetItem(r, 0, level);
         PyList_SetItem(r, 1, prime_flag_sym);
         PyList_SetItem(r, 2, res_type);
         PyList_SetItem(r, 3, link_type);
         PyList_SetItem(r, 4, parent);
         PyList_SetItem(r, 5, parent_spec_py);
      }
   }
   if (PyBool_Check(r))
        Py_INCREF(r);
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *glyco_tree_compare_trees_py(int imol_1, PyObject *res_spec_1_py,
                                      int imol_2, PyObject *res_spec_2_py) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {
         graphics_info_t g;

         coot::residue_spec_t residue_spec_1 = residue_spec_from_py(res_spec_1_py);
         mmdb::Residue *residue_1_p = g.molecules[imol_1].get_residue(residue_spec_1);
         mmdb::Manager *mol_1 = g.molecules[imol_1].atom_sel.mol;
         std::vector<std::string> types_with_no_dictionary =
               g.molecules[imol_1].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
         for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
            g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

         coot::residue_spec_t residue_spec_2 = residue_spec_from_py(res_spec_2_py);
         mmdb::Residue *residue_2_p = g.molecules[imol_2].get_residue(residue_spec_2);
         mmdb::Manager *mol_2 = g.molecules[imol_2].atom_sel.mol;
         types_with_no_dictionary = g.molecules[imol_2].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
         for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
            g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

         coot::glyco_tree_t t_1(residue_1_p, mol_1, g.Geom_p());
         coot::glyco_tree_t t_2(residue_2_p, mol_2, g.Geom_p());

         bool match = t_1.compare_trees(t_2.get_glyco_tree());

         if (match)
            r = Py_True;
      }
   }

   if (PyBool_Check(r))
      Py_XINCREF(r);

   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *glyco_tree_matched_residue_pairs_py(int imol_1, PyObject *res_spec_1_py,
                                              int imol_2, PyObject *res_spec_2_py) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {
         graphics_info_t g;

         coot::residue_spec_t residue_spec_1 = residue_spec_from_py(res_spec_1_py);
         mmdb::Residue *residue_1_p = g.molecules[imol_1].get_residue(residue_spec_1);
         mmdb::Manager *mol_1 = g.molecules[imol_1].atom_sel.mol;
         std::vector<std::string> types_with_no_dictionary =
               g.molecules[imol_1].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
         for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
            g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

         coot::residue_spec_t residue_spec_2 = residue_spec_from_py(res_spec_2_py);
         mmdb::Residue *residue_2_p = g.molecules[imol_2].get_residue(residue_spec_2);
         mmdb::Manager *mol_2 = g.molecules[imol_2].atom_sel.mol;
         types_with_no_dictionary = g.molecules[imol_2].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
         for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
            g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);

         coot::glyco_tree_t t_1(residue_1_p, mol_1, g.Geom_p());
         coot::glyco_tree_t t_2(residue_2_p, mol_2, g.Geom_p());

         std::vector<std::pair<coot::residue_spec_t, coot::residue_spec_t> > pv = t_1.matched_pairs(t_2.get_glyco_tree());
         if (! pv.empty()) {
            PyObject *rt = PyList_New(0);
            for (unsigned int i=0; i<pv.size(); i++) {
               PyObject *r_1_py = residue_spec_to_py(pv[i].first);
               PyObject *r_2_py = residue_spec_to_py(pv[i].second);
               PyObject *list_py = PyList_New(2);
               PyList_SetItem(list_py, 0, r_1_py);
               PyList_SetItem(list_py, 1, r_2_py);
               PyList_Append(rt, list_py);
            }
            r = rt;
         }
      }
   }
   if (PyBool_Check(r))
        Py_INCREF(r);
   return r;
}
#endif

#include "c-interface-scm.hh"

#ifdef USE_GUILE
SCM glyco_tree_internal_distances_fn_scm(int imol, SCM base_residue_spec_scm, const std::string &file_name) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      if (scm_is_true(scm_list_p(base_residue_spec_scm))) {
	 graphics_info_t g;
	 std::pair<bool, coot::residue_spec_t> spec = make_residue_spec(base_residue_spec_scm);
	 if (spec.first)
	    graphics_info_t::molecules[imol].glyco_tree_internal_distances_fn(spec.second, g.Geom_p(), file_name);
	 else
	    std::cout << "WARNING:: Failed to make residue spec "  << std::endl;
      }
   }
   return r;
}

#endif

#ifdef USE_PYTHON
PyObject *glyco_tree_internal_distances_fn_py(int imol, PyObject *base_residue_spec_py,
                                              const std::string &file_name) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      if (PyList_Check(base_residue_spec_py)) {
         graphics_info_t g;
         std::pair<bool, coot::residue_spec_t> spec = make_residue_spec_py(base_residue_spec_py);
         if (spec.first)
            graphics_info_t::molecules[imol].glyco_tree_internal_distances_fn(spec.second, g.Geom_p(), file_name);
         else
            std::cout << "WARNING:: Failed to make residue spec "  << std::endl;
      }
   }
   return r;
}

#endif

/*! \brief rotate the view so that the next main-chain atoms are oriented
 in the same direction as the previous - hence side-chain always seems to be
"up" - set this mode to 1 for reorientation-mode - and 0 for off (standard translation)
*/
void set_reorienting_next_residue_mode(int state) {

   graphics_info_t::reorienting_next_residue_mode = state;

}
