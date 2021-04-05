
#include "validation.hh"
#include "graphics-info.h"

#include "cc-interface.hh"  // for residue_spec_py
#include "coot-utils/c-beta-deviations.hh"


#ifdef USE_PYTHON
PyObject *c_beta_deviations_py(int imol) {

   PyObject *o = Py_False;

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {
	 std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >
	    residue_c_beta_map = coot::get_c_beta_deviations(mol);
	 std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >::const_iterator it;
	 // PyObject *dict_py = PyDict_New();
	 PyObject *outer_list_py = PyList_New(residue_c_beta_map.size());
	 unsigned int counter = 0;
	 for (it=residue_c_beta_map.begin(); it!=residue_c_beta_map.end(); it++) {
	    // multiple alt confs
	    mmdb::Residue *res_key = it->first;
	    const std::map<std::string, coot::c_beta_deviation_t> &value_map = it->second;
	    std::map<std::string, coot::c_beta_deviation_t>::const_iterator it_inner;
	    PyObject *map_dict_py = PyDict_New();
	    for (it_inner=value_map.begin(); it_inner!=value_map.end(); it_inner++) {
	       const std::string alt_conf_key = it_inner->first;
	       const coot::c_beta_deviation_t &cbd = it_inner->second;
	       PyDict_SetItemString(map_dict_py,
				    alt_conf_key.c_str(),
				    PyFloat_FromDouble(cbd.dist));
	    }
	    PyObject *residue_spec_py = residue_spec_to_py(coot::residue_spec_t(res_key));
	    PyObject *residue_pair_py = PyList_New(2);
	    PyList_SetItem(residue_pair_py, 0, residue_spec_py);
	    PyList_SetItem(residue_pair_py, 1, map_dict_py);
	    PyList_SetItem(outer_list_py, counter, residue_pair_py);
	    counter++;
	 }
	 o = outer_list_py;
      }
   }

   if (PyBool_Check(o))
      Py_INCREF(o);
   return o;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
SCM c_beta_deviations_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (graphics_info_t::is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {
	 r = SCM_EOL;
	 std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >
	    residue_c_beta_map = coot::get_c_beta_deviations(mol);
	 std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >::const_iterator it;
	 SCM outer_list_scm = SCM_EOL;
	 unsigned int counter = 0;
	 for (it=residue_c_beta_map.begin(); it!=residue_c_beta_map.end(); it++) {
	    // multiple alt confs
	    mmdb::Residue *res_key = it->first;
	    const std::map<std::string, coot::c_beta_deviation_t> &value_map = it->second;
	    std::map<std::string, coot::c_beta_deviation_t>::const_iterator it_inner;
	    SCM map_dict_scm = SCM_EOL;
	    SCM residue_spec_scm = residue_spec_to_scm(coot::residue_spec_t(res_key));
	    for (it_inner=value_map.begin(); it_inner!=value_map.end(); it_inner++) {
	       const std::string alt_conf_key = it_inner->first;
	       const coot::c_beta_deviation_t &cbd = it_inner->second;
	       SCM item_scm = scm_list_2(scm_from_locale_string(alt_conf_key.c_str()),
					 scm_from_double(cbd.dist));
	       map_dict_scm = scm_cons(item_scm, map_dict_scm);
	    }
	    SCM l_scm = scm_list_2(residue_spec_scm, scm_reverse(map_dict_scm));
	    r = scm_cons(l_scm, r);
	 }
      }
   }
   return r;

}
#endif // USE_GUILE
