/*
 * src/cc-interface-map-utils.cc
 *
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#ifdef USE_PYTHON
#include <Python.h>
#endif
#include "graphics-info.h"
#include "cc-interface.hh"

#include "coot-utils/coot-map-utils.hh"

#include "utils/logging.hh"
extern logging logger;


/* ------------------------------------------------------------------------- */
/*                      Map Display Control                                  */
/* ------------------------------------------------------------------------- */

/*! \name  Map Display Control */
void undisplay_all_maps_except(int imol_map) {

   graphics_info_t g;
   int n = g.n_molecules();
   for (int i=0; i<n; i++) {
      if (g.is_valid_map_molecule(i)) {
         if (i == imol_map) {
            g.molecules[i].set_map_is_displayed(true); // just a state change
            set_display_control_button_state(i, "Displayed", true);
         } else {
            if (g.molecules[i].is_displayed_p()) {
               g.molecules[i].set_map_is_displayed(false);
               set_display_control_button_state(i, "Displayed", false);
            }
         }
      }
   }
   g.graphics_draw();
}


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

      istate = g.molecules[imol_new].read_ccp4_map(filename, is_diff_map_flag, *g.map_glob_extensions);

      if (istate > -1) { // not a failure
	 g.scroll_wheel_map = imol_new;  // change the current scrollable map.
	 g.activate_scroll_radio_button_in_display_manager(imol_new);
      } else {
	 g.erase_last_molecule();
	 // std::cout << "Read map " << filename << " failed" << std::endl;
         logger.log(log_t::INFO, "Read of map", filename, "failed");
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
std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
          std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
map_to_model_correlation_stats_per_residue_range(int imol, const std::string &chain_id, int imol_map,
                                                 unsigned int n_residue_per_residue_range,
                                                 short int exclude_NOC_flag) {

   bool exclude_NOC = false;
   if (exclude_NOC_flag)
      exclude_NOC = true;
   std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
             std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> > m;

   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      if (g.is_valid_map_molecule(imol_map)) {
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap(g.molecules[imol_map].xmap);
         m = coot::util::map_to_model_correlation_stats_per_residue_run(mol, chain_id, xmap, n_residue_per_residue_range, exclude_NOC);
      }
   }

   return m;
}


#ifdef USE_PYTHON
PyObject *
map_to_model_correlation_stats_per_residue_range_py(int imol, const std::string &chain_id, int imol_map,
                                                    unsigned int n_residue_per_residue_range,
                                                    short int exclude_NOC_flag) {

   std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
             std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> > m =
      map_to_model_correlation_stats_per_residue_range(imol, chain_id, imol_map, n_residue_per_residue_range, exclude_NOC_flag);

   unsigned int  n_aa_items = m.first.size();  // all atom and side chain
   unsigned int  n_sc_items = m.first.size();
   PyObject *o = PyList_New(2);
   PyObject *o_0 = PyList_New(n_aa_items);
   PyObject *o_1 = PyList_New(n_sc_items);

   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;

   unsigned int count = 0;
   for (it=m.first.begin(); it!=m.first.end(); ++it) {
      const coot::residue_spec_t &spec(it->first);
      const coot::util::density_correlation_stats_info_t &stats(it->second);
      PyObject *spec_py = residue_spec_to_py(spec);
      PyObject *stats_py = PyList_New(2);
      PyList_SetItem(stats_py, 0, PyLong_FromLong(stats.n));
      PyList_SetItem(stats_py, 1, PyFloat_FromDouble(stats.correlation()));
      PyObject *list_py = PyList_New(2);
      PyList_SetItem(list_py, 0, spec_py);
      PyList_SetItem(list_py, 1, stats_py);
      PyList_SetItem(o_0, count, list_py);
      count++;
   }
   count = 0;
   for (it=m.second.begin(); it!=m.second.end(); ++it) {
      const coot::residue_spec_t &spec(it->first);
      const coot::util::density_correlation_stats_info_t &stats(it->second);
      PyObject *spec_py = residue_spec_to_py(spec);
      PyObject *stats_py = PyList_New(2);
      PyList_SetItem(stats_py, 0, PyLong_FromLong(stats.n));
      PyList_SetItem(stats_py, 1, PyFloat_FromDouble(stats.correlation()));
      PyObject *list_py = PyList_New(2);
      PyList_SetItem(list_py, 0, spec_py);
      PyList_SetItem(list_py, 1, stats_py);
      PyList_SetItem(o_1, count, list_py);
      count++;
   }
   PyList_SetItem(o, 0, o_0);
   PyList_SetItem(o, 1, o_1);
   return o;
}
#endif

#ifdef USE_GUILE
SCM map_to_model_correlation_stats_per_residue_range_scm(int imol, const std::string &chain_id, int imol_map,
                                                         unsigned int n_residue_per_residue_range,
                                                         short int exclude_NOC_flag) {

   std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
             std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> > m =
      map_to_model_correlation_stats_per_residue_range(imol, chain_id, imol_map, n_residue_per_residue_range, exclude_NOC_flag);

   SCM r_0 = SCM_EOL;
   SCM r_1 = SCM_EOL;
   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;
   for (it=m.first.begin(); it!=m.first.end(); ++it) {
      const coot::residue_spec_t &spec(it->first);
      const coot::util::density_correlation_stats_info_t &stats(it->second);
      SCM spec_scm = residue_spec_to_scm(spec);
      SCM stats_scm = scm_list_2(scm_from_int(stats.n), scm_from_double(stats.correlation()));
      SCM item = scm_list_2(spec_scm, stats_scm);
      r_0 = scm_cons(r_0, item);
   }
   for (it=m.second.begin(); it!=m.second.end(); ++it) {
      const coot::residue_spec_t &spec(it->first);
      const coot::util::density_correlation_stats_info_t &stats(it->second);
      SCM spec_scm = residue_spec_to_scm(spec);
      SCM stats_scm = scm_list_2(scm_from_int(stats.n), scm_from_double(stats.correlation()));
      SCM item = scm_list_2(spec_scm, stats_scm);
      r_1 = scm_cons(r_1, item);
   }
   return scm_list_2(r_0, r_1);
}
#endif


#include "utils/subprocess.hpp"

// resolution in A.
// The map in imol_fofc_map gets overwritten.
void servalcat_fofc(int imol_model,
                    int imol_fofc_map, const std::string &half_map_1, const std::string &half_map_2,
                    float resolution) {

   // c.f. void sharpen_blur_map_with_resampling_threaded_version(int imol_map, float b_factor, float resample_factor);

   auto servalcat_fofc_func = [] (const std::string &half_map_1, const std::string &half_map_2,
                                  const std::string &pdb_file_name,
                                  const std::string &prefix,
                                  float resolution) {

      try {
         graphics_info_t g;
         std::string output_fn = prefix + std::string(".mtz");
         std::vector<std::string> cmd_list = {"servalcat", "fofc",
            "--halfmaps", half_map_1, half_map_2,
            "--trim", "--trim_mtz", "--resolution", std::to_string(resolution),
            "--model", pdb_file_name, "-o", prefix};
         subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
         if (true) {
            std::cout << "Data : " << obuf.buf.data() << std::endl;
            std::cout << "Data len: " << obuf.length << std::endl;
         }
         g.servalcat_fofc.second = output_fn;
         g.servalcat_fofc.first = true;
      }
      catch (const std::runtime_error &e) {
         // this happens when servalcat fails to run correctly (e.g. input error)
         std::cout << "WARNING:: runtime_error " << e.what() << std::endl;
      }
      catch (const std::exception &e) {
         std::cout << "WARNING:: exception " << e.what() << std::endl;
      }
      catch (...) {
         std::cout << "WARNING:: caught some other error" << std::endl;
      }
   };

   auto check_it = +[] (gpointer data) {

      graphics_info_t g;
      std::cout << "............... running check_it() " << g.servalcat_fofc.first << std::endl;

      if (g.servalcat_fofc.first) {
         const std::string &mtz_file_name = g.servalcat_fofc.second;
         g.servalcat_fofc.first = false; // turn it off
         int imol_map = GPOINTER_TO_INT(data);
         std::cout << "debug:: in check_it() with imol_map " << imol_map << std::endl;
         coot::mtz_to_map_info_t mmi;
         mmi.is_difference_map = true;
         mmi.mtz_file_name = mtz_file_name;
         mmi.f_col   = "DELFWT";
         mmi.phi_col = "PHDELWT";
         mmi.id = "Something or other - what should I put here?";
         std::cout << "............... calling update_self() mtz " << mtz_file_name << std::endl;
         g.molecules[imol_map].update_self(mmi);
         g.graphics_draw();
         return gboolean(FALSE);
      }
      return gboolean(TRUE);
   };

   graphics_info_t g;
   if (g.is_valid_model_molecule(imol_model)) {
      if (g.is_valid_map_molecule(imol_fofc_map)) {
         // fine
      } else {
         clipper::Xmap<float> xmap;
         bool is_em_map = true;
         std::string label = "diff map";
         imol_fofc_map = g.create_molecule();
         g.molecules[imol_fofc_map].install_new_map(xmap, label, is_em_map);
      }
      std::string prefix = g.molecules[imol_fofc_map].get_name();
      std::string model_file_name = std::string("servalcat-fofc-") + g.molecules[imol_model].get_name();
      g.molecules[imol_model].write_pdb_file(model_file_name);
      std::thread thread(servalcat_fofc_func, half_map_1, half_map_2, model_file_name, prefix, resolution);
      thread.detach();

      g.servalcat_fofc.first = false;
      GSourceFunc f = GSourceFunc(check_it);
      std::cout << "debug:: in servalcat_fofc() with imol_fofc_map " << imol_fofc_map << " as user data" << std::endl;
      g_timeout_add(400, f, GINT_TO_POINTER(imol_fofc_map));
   }
}


#include "read-molecule.hh"

// resolution in A.
void servalcat_refine(int imol_model,
                      const std::string &half_map_1, const std::string &half_map_2,
                      const std::string &mask_map, float resolution) {


   // c.f. void sharpen_blur_map_with_resampling_threaded_version(int imol_map, float b_factor, float resample_factor);

   auto servalcat_refine_func = [] (const std::string &half_map_1, const std::string &half_map_2,
                                    const std::string &mask_file_name,
                                    const std::string &input_pdb_file_name,
                                    const std::string &prefix,
                                    float resolution) {

      graphics_info_t g;
      try {
         xdg_t xdg;
         std::string output_pdb_file_name = prefix + std::string(".pdb");
         std::vector<std::string> cmd_list = {"servalcat", "refine_spa_norefmac",
            "--halfmaps", half_map_1, half_map_2,
            // "--mask", mask_file_name,
            "--resolution", std::to_string(resolution),
            "--model", input_pdb_file_name, "-o", prefix};
         subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
         if (true) {
            std::cout << "Data : " << obuf.buf.data() << std::endl;
            std::cout << "Data len: " << obuf.length << std::endl;
         }
         g.servalcat_refine.second = output_pdb_file_name;
      }
      catch (const std::runtime_error &e) {
         // this happens when servalcat fails to run correctly (e.g. input error)
         std::cout << "WARNING:: runtime_error " << e.what() << std::endl;
      }
      catch (const std::exception &e) {
         std::cout << "WARNING:: exception " << e.what() << std::endl;
      }
      catch (...) {
         std::cout << "WARNING:: caught some other error" << std::endl;
      }
      g.servalcat_refine.first = true; // mark refinement process has finished
   };

   auto check_it = +[] (gpointer data) {

      graphics_info_t g;
      std::cout << "debug:: running servalcat_refine() check_it() " << g.servalcat_fofc.first << std::endl;

      if (g.servalcat_refine.first) {
         const std::string &pdb_file_name = g.servalcat_refine.second;
         g.servalcat_refine.first = false; // turn it off

         std::cout << "............... servalcat_refine() check_it() read_pdb() " << pdb_file_name << std::endl;
         if (pdb_file_name.empty()) {
            std::cout << "servalcat_refine() check_it() something-went-wrong pdb_file_name is empty " << std::endl;
         } else {
            read_pdb(pdb_file_name);

            if (false) { // currently we are not reading maps
               coot::mtz_to_map_info_t mmi;
               // nmmi.is_difference_map = true;
               // mmi.mtz_file_name = mtz_file_name;
               mmi.f_col   = "DELFWT";
               mmi.phi_col = "PHDELWT";
               mmi.id = "Refine Something";
            }
            g.graphics_draw();
            return gboolean(FALSE);
         }
      }
      return gboolean(TRUE);
   };

   graphics_info_t g;
   if (g.is_valid_model_molecule(imol_model)) {
      xdg_t xdg;
      std::string n = g.molecules[imol_model].Refmac_name_stub();
      std::string prefix_dir = std::string("servalcat-refine-") + n;
      std::string prefix = (xdg.get_data_home() / prefix_dir).string();
      std::string pdb_file_name = prefix + "_input.pdb";
      g.molecules[imol_model].write_pdb_file(pdb_file_name);
      std::thread thread(servalcat_refine_func, half_map_1, half_map_2, mask_map, pdb_file_name, prefix, resolution);
      thread.detach();

      g.servalcat_refine.first = false;
      GSourceFunc f = GSourceFunc(check_it);
      g_timeout_add(400, f, nullptr);
   }

}

#include "c-interface-python.hh" // because we use display_python().
#include "python-3-interface.hh"

void add_toolbar_subprocess_button(const std::string &button_label,
                                   const std::string &subprocess_command,
                                   PyObject *arg_list,
                                   PyObject *on_completion_function, // should these be strings?
                                   PyObject *on_completion_args) {

   auto type_check = [] (PyObject *obj) {

      if (obj == nullptr) {
         return "NULL";
      }

      PyTypeObject* type = obj->ob_type;
      if (type == &PyLong_Type) {
         return "int"; // Python 3 integers are longs
      } else if (type == &PyFloat_Type) {
         return "float";
      } else if (type == &PyUnicode_Type) {
         return "str"; // Python 3 strings are unicode
      } else if (type == &PyBool_Type) {
         return "bool";
      } else if (type == &PyList_Type) {
         return "list";
      } else if (type == &PyTuple_Type) {
         return "tuple";
      } else if (type == &PyDict_Type) {
         return "dict";
      // I don't know what to do about Py_None
      // } else if (type == &PyNone_Type) {
      //    return "NoneType";
      } else if (type == &PyBytes_Type) {
         return "bytes";
      } else if (type == &PyByteArray_Type) {
         return "bytearray";
      } else {
         return type->tp_name; // Fallback to the type's name
      }
   };

   struct py_transfer_t {
      PyObject *on_completion_function;
      PyObject *on_completion_args;
      std::vector<std::string> cmd_list;
      bool termination_condtion;
   };

   if (PyList_Check(arg_list)) {
      unsigned int l1 = PyObject_Length(arg_list);
      std::vector<std::string> args;
      for (unsigned int i=0; i<l1; i++) {
         PyObject *o = PyList_GetItem(arg_list, i);
         if (PyUnicode_Check(o)) {
            std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(o));
            args.push_back(s);
         }
      }

      auto on_button_clicked = +[] (GtkButton *button, gpointer data) {

         auto my_subproc = [] (const std::vector<std::string> &cmd_list, py_transfer_t *pt) {
            try {
               subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
               if (true) {
                  std::cout << "Data : " << obuf.buf.data() << std::endl;
                  std::cout << "Data len: " << obuf.length << std::endl;
                  pt->termination_condtion = true;
               }
            }
            catch (...) {
               std::cout << "WARNING:: add_toolbar_subprocess_button() my_subproc() caught some error" << std::endl;
            }
         };

         auto check_it = +[] (gpointer data) {

            py_transfer_t *pt = reinterpret_cast<py_transfer_t *>(data);
            // std::cout << "check_it: " << pt->termination_condtion << std::endl;
            if (pt->termination_condtion) {
               PyObject *return_val = nullptr;
#if 0 // check version of python here - less than 3.13 is my guess
               return_val = PyEval_CallObject(pt->on_completion_function, pt->on_completion_args);
#endif
               std::cout << "DEBUG:: return_val " << return_val << std::endl;
               PyObject *error_thing = PyErr_Occurred();
               if (! error_thing) {
                  // std::cout << "INFO:: check_it() No Python error on callable check" << std::endl;
                  logger.log(log_t::INFO, "check_it() No Python error on callable check");
               } else {
                  std::cout << "ERROR:: while executing check_it() a python error occured " << std::endl;
                  PyObject *type, *value, *traceback;
                  PyErr_Fetch(&type, &value, &traceback);
                  PyErr_NormalizeException(&type, &value, &traceback);
                  PyObject *exception_string = PyObject_Repr(value);
                  const char *em = myPyString_AsString(exception_string);
                  std::cout << "ERROR:: " << em << std::endl;
                  Py_XDECREF(value);
                  Py_XDECREF(traceback);
                  Py_XDECREF(type);
               }

               if (return_val) {
                  PyObject *dp = display_python(return_val);
                  std::cout << "DEBUG:: return val as string: " << PyBytes_AS_STRING(PyUnicode_AsUTF8String(dp)) << std::endl;
               }
               graphics_info_t g;
               g.graphics_draw();
               return gboolean(FALSE);
            } else {
               return gboolean(TRUE);
            }
         };

         py_transfer_t *pt = reinterpret_cast<py_transfer_t *>(data);
         std::vector<std::string> command_list = pt->cmd_list;
         std::thread thread(my_subproc, command_list, pt);
         thread.detach();
         GSourceFunc f = GSourceFunc(check_it);
         g_timeout_add(400, f, pt);

      };

      int tuple_state = PyTuple_Check(on_completion_args);
      int unicode_state = PyUnicode_Check(on_completion_args);
      std::cout << "debug:: on_completion_args tuple-state: " << tuple_state << std::endl;
      std::cout << "debug:: on_completion_args unicode-state: " << unicode_state << std::endl;

      if (on_completion_args) {
         PyObject *dp = display_python(on_completion_args);
         if (dp)
            std::cout << "DEBUG:: on_completion_args: " << PyUnicode_AsUTF8String(dp) << std::endl;
         else
            std::cout << "DEBUG:: on_completion_args display_python null " << std::endl;
         PyObject *error_thing = PyErr_Occurred();
         if (! error_thing) {
            // std::cout << "INFO:: check_it() No Python error on printing on_completion_args" << std::endl;
            logger.log(log_t::INFO, "check_it() No Python error on printing on_completion_args");
         } else {
            std::cout << "ERROR:: while pringing on_completion_args a python error occured " << std::endl;
            PyObject *type, *value, *traceback;
            PyErr_Fetch(&type, &value, &traceback);
            PyErr_NormalizeException(&type, &value, &traceback);
            PyObject *exception_string = PyObject_Repr(value);
            const char *em = myPyString_AsString(exception_string);
            std::cout << "ERROR:: " << em << std::endl;
            Py_XDECREF(value);
            Py_XDECREF(traceback);
            Py_XDECREF(type);
         }
      }

      std::string oca_type = type_check(on_completion_args);
      std::cout << "oca_type " << oca_type << std::endl;

      py_transfer_t *py_transfer_p = new py_transfer_t;
      py_transfer_p->on_completion_function = on_completion_function;
      py_transfer_p->on_completion_args     = on_completion_args;
      py_transfer_p->cmd_list = args;
      py_transfer_p->cmd_list.insert(py_transfer_p->cmd_list.begin(), subprocess_command);
      py_transfer_p->termination_condtion = false;

      GtkWidget *button = gtk_button_new_with_label(button_label.c_str());
      g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(on_button_clicked), py_transfer_p);
      GtkWidget *toolbar = widget_from_builder("main_toolbar");
      GtkWidget *toolbar_hbox = widget_from_builder("main_window_toolbar_hbox");
      gtk_box_append(GTK_BOX(toolbar_hbox), button);

   }

}

#include "c-interface.h" // for is_valid_model_molecule()))

// resolution in A.
// The map in imol_fofc_map gets overwritten.
void emplacement_local(int imol_model,
                       const std::string &half_map_1, const std::string &half_map_2,
                       int search_centre_x,
                       int search_centre_y,
                       int search_centre_z,
                       float resolution) {

   // this doesn't actually work yet. It is a placeholder to be fixed up later.

   auto emplace_local_func = [] (const std::vector<std::string> &cmd_list ) {

      try {
         subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
         if (true) {
            std::cout << "Data : " << obuf.buf.data() << std::endl;
            std::cout << "Data len: " << obuf.length << std::endl;
         }
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING:: runtime_error " << e.what() << std::endl;
      }
      catch (const std::exception &e) {
         std::cout << "WARNING:: exception " << e.what() << std::endl;
      }
      catch (...) {
         std::cout << "WARNING:: caught some other error" << std::endl;
      }
   };

   auto check_it = +[] (gpointer data) {
      std::cout << "checking..." << std::endl;
   };

   if (is_valid_model_molecule(imol_model)) {
      char *ccp4_ev = getenv("CCP4");
      std::string model_file_name = "emplace_local_in.pdb";
      write_pdb_file(imol_model, model_file_name.c_str());
      if (ccp4_ev) {
         std::string ccp4_prfx(ccp4_ev);
         std::string em_placement_py_file = ccp4_prfx + "/ccp4-9/lib/python3.9/site-packages" +
            "/phaser_voyager/src/New_Voyager/scripts/emplace_local.py";
         std::vector<std::string> command_list = {
            "ccp4-python",
            em_placement_py_file,
            "--model_file",
            model_file_name,
            "--map1", half_map_1,
            "--map2", half_map_2,
            "--d_min", std::to_string(resolution),
            "--sphere_center",
            std::to_string(search_centre_x),
            std::to_string(search_centre_y),
            std::to_string(search_centre_z)
         };

         std::thread thread(emplace_local_func, command_list);
         thread.detach();

         GSourceFunc f = GSourceFunc(check_it);
         g_timeout_add(400, f, nullptr);

      }
   }

}


// in c-interface-gui.cc
// fill_combobox_with_map_options(combobox_1, callback);
// but not in a header
int fill_combobox_with_map_options(GtkWidget *combobox, GCallback signalfunc);

void show_map_partition_by_chain_dialog() {

   GtkWidget *dialog = widget_from_builder("map_partition_by_chain_dialog");

   // fill map partion dialog
   GtkWidget *combobox_1 = widget_from_builder("map_partition_by_chain_map_combobox");
   GtkWidget *combobox_2 = widget_from_builder("map_partition_by_chain_model_combobox");

   int imol_active = 0;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first)
      imol_active = active_atom.second.first;

   GCallback callback = G_CALLBACK(nullptr);
   graphics_info_t g;
   g.new_fill_combobox_with_coordinates_options(combobox_2, callback, imol_active);

   // in c-interface-gui.cc
   fill_combobox_with_map_options(combobox_1, callback);

   gtk_widget_set_visible(dialog, TRUE);

}

// use this for interactive
void map_partition_by_chain_threaded(int imol_map, int imol_model) {

   auto map_partition_func = [] (const clipper::Xmap<float> &xmap, mmdb::Manager *mol) {

      std::string &state_string = graphics_info_t::map_partition_results_state_string;
      graphics_info_t::map_partition_results_state = 1;
      graphics_info_t::map_partition_results = coot::util::partition_map_by_chain(xmap, mol, &state_string);
      graphics_info_t::map_partition_results_state = 0;
   };

   auto check_map_partition_results = +[] (gpointer data) {

      bool keep_going = TRUE;
      int imol_map = GPOINTER_TO_INT(data);
      // std::cout << "checking...  " << graphics_info_t::map_partition_results_state << std::endl;
      GtkWidget *label = widget_from_builder("partition_map_by_chain_status_label");
      if (label) {
         gtk_widget_set_visible(label, TRUE);
         gtk_label_set_text(GTK_LABEL(label), graphics_info_t::map_partition_results_state_string.c_str());
      }
      if (graphics_info_t::map_partition_results_state == 1) {
         // keep going
      } else {
         if (! graphics_info_t::map_partition_results.empty()) {
            bool is_em_map = graphics_info_t::molecules[imol_map].is_EM_map();
            for (const auto &mi : graphics_info_t::map_partition_results) {
               std::string chain_id = mi.first;
               int imol_for_map = graphics_info_t::create_molecule();
               std::string label = "Partitioned map Chain " + chain_id;
               graphics_info_t::molecules[imol_for_map].install_new_map(mi.second, label, is_em_map);
            }
            // Let's just check that the user didn't delete the original map in the meantime...
            if (graphics_info_t::is_valid_map_molecule(imol_map))
               graphics_info_t::molecules[imol_map].set_map_is_displayed(false);

            if (label) {
               gtk_label_set_text(GTK_LABEL(label), "");
               gtk_widget_set_visible(label, FALSE);
            }

            keep_going = FALSE; // turn off the timeout

            graphics_info_t::graphics_draw();
         }
      }
      return keep_going;
   };

   std::vector<int> v;
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol_model)) {
      if (g.is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
         mmdb::Manager *mol = graphics_info_t::molecules[imol_model].atom_sel.mol;
         std::thread t(map_partition_func, xmap, mol);
         t.detach();
         GSourceFunc cb = GSourceFunc(check_map_partition_results);
         g_timeout_add(1000, cb, GINT_TO_POINTER(imol_map));
      }
   }
}

// use this version for scripting
std::vector<int> map_partition_by_chain(int imol_map, int imol_model) {

   std::vector<int> v;
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol_model)) {
      if (g.is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
         mmdb::Manager *mol = graphics_info_t::molecules[imol_model].atom_sel.mol;
         std::string info_string;
         std::vector<std::pair<std::string, clipper::Xmap<float> > > maps_info =
            coot::util::partition_map_by_chain(xmap, mol, &info_string);
         if (! maps_info.empty()) {
            for (const auto &mi : maps_info) {
               std::string chain_id = mi.first;
               int imol_for_map = g.create_molecule();
               std::string label = "Partioned map Chain " + chain_id;
               bool is_em_map = g.molecules[imol_map].is_EM_map();
               g.molecules[imol_for_map].install_new_map(mi.second, label, is_em_map);
               v.push_back(imol_for_map);
            }
            g.molecules[imol_map].set_map_is_displayed(false);
         }
      }
   }
   g.graphics_draw();
   return v;
}
