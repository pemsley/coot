
#ifdef USE_PYTHON
#include <Python.h>
#endif
#include "rama_plot.hh"
#include "graphics-info.h"
#include "cc-interface.hh"

#include "coot-utils/coot-map-utils.hh"



/*! \brief read a CCP4 map or a CNS map (despite the name). */
int read_ccp4_map(const std::string &filename, int is_diff_map_flag) {
   return handle_read_ccp4_map(filename, is_diff_map_flag);
}

int
handle_read_ccp4_map(const std::string &filename, int is_diff_map_flag) {

   int istate = -1;
   if (true) {
      graphics_info_t g;
      int imol_new = graphics_info_t::create_molecule();

      istate = g.molecules[imol_new].read_ccp4_map(filename, is_diff_map_flag,
						   *graphics_info_t::map_glob_extensions);

      if (istate > -1) { // not a failure
	 g.scroll_wheel_map = imol_new;  // change the current scrollable map.
	 g.activate_scroll_radio_button_in_display_manager(imol_new);
      } else {
	 g.erase_last_molecule();
	 std::cout << "Read map " << filename << " failed" << std::endl;
	 std::string s = "Read map ";
	 s += filename;
	 s += " failed.";
	 g.add_status_bar_text(s);
      }
      g.graphics_draw();
   }
   return istate;
}


//! \brief map to model density statistics, reported per residue, the middle residue
//!        of a range of residues
std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> 
map_to_model_correlation_stats_per_residue_range(int imol, const std::string &chain_id, int imol_map,
                                                 unsigned int n_residue_per_residue_range) {

   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> m;

   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      if (g.is_valid_map_molecule(imol_map)) {
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
         m = coot::util::map_to_model_correlation_stats_per_residue_run(mol, chain_id, xmap, n_residue_per_residue_range);
      }
   }

   std::cout << "-------------------- debug map_to_model_correlation_stats_per_residue_run size " << m.size()
             << std::endl;
   return m;
}


#ifdef USE_PYTHON
PyObject *
map_to_model_correlation_stats_per_residue_range_py(int imol, const std::string &chain_id, int imol_map,
                                                    unsigned int n_residue_per_residue_range) {

   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> m =
      map_to_model_correlation_stats_per_residue_range(imol, chain_id, imol_map, n_residue_per_residue_range);

   unsigned int  n_items = m.size();
   PyObject *o = PyList_New(n_items);

   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;

   unsigned int count = 0;
   for (it=m.begin(); it!=m.end(); ++it) {
      const coot::residue_spec_t &spec(it->first);
      const coot::util::density_correlation_stats_info_t &stats(it->second);
      PyObject *spec_py = residue_spec_to_py(spec);
      PyObject *stats_py = PyList_New(2);
      PyList_SetItem(stats_py, 0, PyInt_FromLong(stats.n));
      PyList_SetItem(stats_py, 1, PyFloat_FromDouble(stats.correlation()));
      PyObject *list_py = PyList_New(2);
      PyList_SetItem(list_py, 0, spec_py);
      PyList_SetItem(list_py, 1, stats_py);
      PyList_SetItem(o, count, list_py);
      count++;
   }
   return o;
}
#endif

#ifdef USE_GUILE
SCM map_to_model_correlation_stats_per_residue_range_scm(int imol, const std::string &chain_id, int imol_map,
                                                         unsigned int n_residue_per_residue_range) {

   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> m =
      map_to_model_correlation_stats_per_residue_range(imol, chain_id, imol_map, n_residue_per_residue_range);

   SCM r = SCM_EOL;
   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;
   for (it=m.begin(); it!=m.end(); ++it) {
      const coot::residue_spec_t &spec(it->first);
      const coot::util::density_correlation_stats_info_t &stats(it->second);
      SCM spec_scm = residue_spec_to_scm(spec);
      SCM stats_scm = scm_list_2(SCM_MAKINUM(stats.n), scm_double2num(stats.correlation()));
      SCM item = scm_list_2(spec_scm, stats_scm);
      r = scm_cons(r, item);
   }
   return r;
}
#endif
