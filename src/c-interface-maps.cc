/* src/c-interface-maps.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2008, 2009, 2011, 2012 by The University of Oxford
 * Copyright 2014, 2015 by Medical Research Council
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

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp> // for to_string()

#include "compat/coot-sysdep.h"

// for stat()
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <future>
#include <chrono>

#include "coot-utils/coot-map-utils.hh" // for variance map
#include "coot-utils/xmap-stats.hh"
#include "skeleton/BuildCas.h"
#include "skeleton/graphical_skel.h"
#include "coot-utils/xmap-stats.hh"

#ifdef HAVE_GOOCANVAS
#include "goograph/goograph.hh"
#endif

#include "c-interface-mmdb.hh"
#include "c-interface-python.hh"


#include "graphics-info.h"
#include "cc-interface.hh"
#include "c-interface.h"
#include "c-interface-gui.hh"
#include "c-interface-gtk-widgets.h"
#include "widget-headers.hh"

#include "analysis/kolmogorov.hh"
#include "analysis/stats.hh"

#include "read-molecule.hh" // 20230621-PE now with std::string args

#include "utils/logging.hh"
extern logging logger;


/*  ------------------------------------------------------------------------ */
/*                   Maps -                                                  */
/*  ------------------------------------------------------------------------ */
/*! \brief Calculate SFs from an MTZ file and generate a map.
 @return the new molecule number. */
int map_from_mtz_by_calc_phases(const char *mtz_file_name,
                                const char *f_col,
                                const char *sigf_col,
                                int imol_coords) {

   int ir = -1; // return value
   graphics_info_t g;
   if (is_valid_model_molecule(imol_coords)) {
      int imol_map = g.create_molecule();
      std::string m(mtz_file_name);
      std::string f(f_col);
      std::string s(sigf_col);
      atom_selection_container_t a = g.molecules[imol_coords].atom_sel;
      short int t = molecule_map_type::TYPE_2FO_FC;
      int istat = g.molecules[imol_map].make_map_from_mtz_by_calc_phases(imol_map,m,f,s,a,t);
      if (istat != -1) {
         graphics_draw();
         ir = imol_map;
      } else {
         ir = -1; // error
         graphics_info_t::erase_last_molecule();
      }
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("map-from-mtz-by-calc-phases");
   command_strings.push_back(mtz_file_name);
   command_strings.push_back(f_col);
   command_strings.push_back(sigf_col);
   command_strings.push_back(graphics_info_t::int_to_string(imol_coords));
   add_to_history(command_strings);
   return ir;
}

#include "cc-interface-scripting.hh"

/*! \brief fire up a GUI, which asks us which model molecule we want
  to calc phases from.  On "OK" button there, we call
  map_from_mtz_by_refmac_calc_phases() */
void calc_phases_generic(const char *mtz_file_name) {

   if (coot::file_exists(mtz_file_name)) {
      graphics_info_t g;
      coot::mtz_column_types_info_t r = coot::get_mtz_columns(mtz_file_name);
      if (r.f_cols.size() == 0) {
	 std::cout << "No Fobs found in " << mtz_file_name << std::endl;
	 std::string s =  "No Fobs found in ";
	 s += mtz_file_name;
	 g.add_status_bar_text(s);
      } else {
	 if (r.sigf_cols.size() == 0) {
	    std::cout << "No SigFobs found in " << mtz_file_name << std::endl;
	    std::string s =  "No SigFobs found in ";
	    s += mtz_file_name;
	    g.add_status_bar_text(s);
	 } else {
	    // normal path:
	    std::string f_obs_col = r.f_cols[0].column_label;
	    std::string sigfobs_col = r.sigf_cols[0].column_label;
	    std::vector<std::string> v;
	    v.push_back("refmac-for-phases-and-make-map");
	    // BL says:: dunno if we need the backslashing here, but just do it in case
	    v.push_back(coot::util::single_quote(coot::util::intelligent_debackslash(mtz_file_name)));
	    v.push_back(coot::util::single_quote(f_obs_col));
	    v.push_back(coot::util::single_quote(sigfobs_col));
	    std::string c = languagize_command(v);
	    std::cout << "command: " << c << std::endl;
#ifdef USE_GUILE
	    safe_scheme_command(c);
#else
#ifdef USE_PYTHON
	    safe_python_command(c);
#endif
#endif
	 }
      }
      std::vector<std::string> command_strings;
      command_strings.push_back("calc-phases-generic");
      command_strings.push_back(mtz_file_name);
      add_to_history(command_strings);
   }
}

/*! \brief Calculate SFs (using refmac optionally) from an MTZ file
  and generate a map. Get F and SIGF automatically (first of their
  type) from the mtz file.

@return the new molecule number, -1 on a problem. */
int map_from_mtz_by_refmac_calc_phases(const char *mtz_file_name,
				       const char *f_col,
				       const char *sigf_col,
				       int imol_coords) {

   int istat = -1;

   // ha! I was supposed to fill this in!

   std::vector<std::string> command_strings;
   command_strings.push_back("map-from-mtz-by-refmac-calc-phases");
   command_strings.push_back(mtz_file_name);
   command_strings.push_back(f_col);
   command_strings.push_back(sigf_col);
   command_strings.push_back(graphics_info_t::int_to_string(imol_coords));
   add_to_history(command_strings);
   return istat;
}

#ifdef USE_PYTHON
/*! \brief Calculate structure factors and make a 2FoFC map and a Fo-Fc map updating the given
   molecule numbers for those maps - if thase molecule ids are not valid maps, them generate
   new maps (return the model number information in the returned object) */
PyObject *calculate_maps_and_stats_py(int imol_model,
                                      int imol_map_with_data_attached,
                                      int imol_map_2fofc,
                                      int imol_map_fofc) {

   auto pythonize_stats = [] (const coot::util::sfcalc_genmap_stats_t &stats) {
                             PyObject *c = PyList_New(5);
                             PyList_SetItem(c, 0, PyFloat_FromDouble(stats.r_factor));
                             PyList_SetItem(c, 1, PyFloat_FromDouble(stats.free_r_factor));
                             PyList_SetItem(c, 2, PyFloat_FromDouble(stats.bulk_solvent_volume));
                             PyList_SetItem(c, 3, PyFloat_FromDouble(stats.bulk_correction));
                             unsigned int n_items = stats.loc_table.size();
                             PyObject *table_py = PyList_New(n_items);
                             for (unsigned int i=0; i<n_items; i++) {
                                const auto &item = stats.loc_table.items[i];
                                PyObject *item_py = PyList_New(3);
                                PyList_SetItem(item_py, 0, PyFloat_FromDouble(item.invresolsq));
                                PyList_SetItem(item_py, 1, PyFloat_FromDouble(item.scale));
                                PyList_SetItem(item_py, 2, PyFloat_FromDouble(item.lack_of_closure));
                                PyList_SetItem(table_py, i, item_py);
                             }
                             PyList_SetItem(c, 4, table_py);
                             return c;
   };

   auto make_status_bar_text = [] (const coot::util::sfcalc_genmap_stats_t &stats) {
      std::string s;
      s += "  R-factor: ";
      s += coot::util::float_to_string_using_dec_pl(100.0 * stats.r_factor, 2);
      s += " Free-R-factor: ";
      s += coot::util::float_to_string_using_dec_pl(100.0 * stats.free_r_factor, 2);
      return s;
   };

   coot::util::sfcalc_genmap_stats_t stats;

   // construct return value
   //
   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol_model)) {
      graphics_info_t g;
      if (is_valid_map_molecule(imol_map_2fofc)) {
         if (is_valid_map_molecule(imol_map_fofc)) {
            clipper::Xmap<float> &xmap_2fofc = g.molecules[imol_map_2fofc].xmap;
            clipper::Xmap<float> &xmap_fofc  = g.molecules[imol_map_fofc].xmap;
            stats = g.sfcalc_genmaps_using_bulk_solvent(imol_model, imol_map_2fofc, &xmap_2fofc, &xmap_fofc);
            g.molecules[imol_map_2fofc].set_mean_and_sigma(false, g.ignore_pseudo_zeros_for_map_stats);
            g.molecules[imol_map_fofc ].set_mean_and_sigma(false, g.ignore_pseudo_zeros_for_map_stats);
            float cls_2fofc = g.molecules[imol_map_2fofc].get_contour_level_by_sigma();
            float cls_fofc  = g.molecules[imol_map_fofc].get_contour_level_by_sigma();
            g.molecules[imol_map_2fofc].set_contour_level_by_sigma(cls_2fofc); // does an update
            g.molecules[imol_map_fofc].set_contour_level_by_sigma(cls_fofc);   // does an update
            std::string sbt = make_status_bar_text(stats);
            add_status_bar_text(sbt.c_str());
            r = pythonize_stats(stats);
         }
      }
   }

   graphics_info_t g;
   updating_model_molecule_parameters_t ummp(imol_model, imol_map_2fofc, imol_map_2fofc, imol_map_fofc);
   g.calculate_new_rail_points(ummp);

   graphics_draw();
   std::vector<coot::command_arg_t> commands;
   std::string cmd = "calculate-maps-and-stats";
   commands.push_back(imol_model);
   commands.push_back(imol_map_with_data_attached);
   commands.push_back(imol_map_2fofc);
   commands.push_back(imol_map_fofc);
   add_to_history_typed(cmd, commands);

   if (PyBool_Check(r))
      Py_XINCREF(r);
   return r;
}
#endif



/*  ----------------------------------------------------------------------- */
/*                  Display lists                                           */
/*  ----------------------------------------------------------------------- */
// 20230501-PE remove this
void set_display_lists_for_maps(int istat) {

}

int display_lists_for_maps_state() {

   return graphics_info_t::display_lists_for_maps_flag;
}

/* update the maps to the current position - rarely needed */
void update_maps() {
   for(int ii=0; ii<graphics_info_t::n_molecules(); ii++) {
      if (is_valid_map_molecule(ii)) {
         // std::cout << "DEBUG:: updating " << ii << std::endl;
         graphics_info_t::molecules[ii].update_map(graphics_info_t::auto_recontour_map_flag);
      }
   }
}

void swap_map_colours(int imol1, int imol2) {

   if (is_valid_map_molecule(imol1)) {
      if (is_valid_map_molecule(imol2)) {
         graphics_info_t g;
         std::pair<GdkRGBA, GdkRGBA> map_1_colours = g.molecules[imol1].get_map_colours();
         std::pair<GdkRGBA, GdkRGBA> map_2_colours = g.molecules[imol2].get_map_colours();
         short int main_or_secondary = 0; // main
         g.molecules[imol1].handle_map_colour_change(map_2_colours.first,
                                                     g.swap_difference_map_colours,
                                                     main_or_secondary,
                                                     g.get_rotation_centre_co(),
                                                     g.box_radius_xray);
         g.molecules[imol2].handle_map_colour_change(map_1_colours.first,
                                                     g.swap_difference_map_colours,
                                                     main_or_secondary,
                                                     g.get_rotation_centre_co(),
                                                     g.box_radius_xray);
         if (graphics_info_t::display_mode_use_secondary_p()) {
            g.make_gl_context_current(graphics_info_t::GL_CONTEXT_SECONDARY);
            main_or_secondary = 1; // secondary
            g.molecules[imol1].handle_map_colour_change(map_2_colours.second,
                                                        g.swap_difference_map_colours,
                                                        main_or_secondary,
                                                        g.get_rotation_centre_co(),
                                                        g.box_radius_xray);
            g.molecules[imol2].handle_map_colour_change(map_1_colours.second,
                                                        g.swap_difference_map_colours,
                                                        main_or_secondary,
                                                        g.get_rotation_centre_co(),
                                                        g.box_radius_xray);
            g.make_gl_context_current(graphics_info_t::GL_CONTEXT_MAIN);
         }
      }
   }
   std::string cmd = "swap-map-colours";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol1);
   args.push_back(imol2);
   add_to_history_typed(cmd, args);
}

void set_keep_map_colour_after_refmac(int istate) {
   std::string cmd = "set-keep-map-colour-after-refmac";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
   graphics_info_t::swap_pre_post_refmac_map_colours_flag = istate;
}

int keep_map_colour_after_refmac_state() {
   add_to_history_simple("keep_map_colour_after_refmac_state");
   return graphics_info_t::swap_pre_post_refmac_map_colours_flag;
}


void show_select_map_frame() {
   graphics_info_t g;
   g.show_select_map_frame();
   add_to_history_simple("show-select-map-frame");
}

int read_mtz(const char *mtz_file_name,
             const char *f_col, const char *phi_col,
             const char *weight,
             int use_weights, int is_diff_map) {

   return make_and_draw_map(mtz_file_name, f_col, phi_col, weight, use_weights, is_diff_map);

}


// return the new molecule number
int make_and_draw_map(const char* mtz_file_name,
                      const char *f_col, const char *phi_col,
                      const char *weight_col, int use_weights,
                      int is_diff_map) {

   graphics_info_t g;
   int imol = -1; // failure initially.
   struct stat buf;

   std::string f_col_str(f_col);
   std::string phi_col_str(phi_col);
   std::string weight_col_str("");
   if (use_weights)
      weight_col_str = std::string(weight_col);

   int status = stat(mtz_file_name, &buf);
   //
   if (status != 0) {
      std::cout << "WARNING:: Can't find file " << mtz_file_name << std::endl;
      if (S_ISDIR(buf.st_mode)) {
         std::cout << mtz_file_name << " is a directory! " << std::endl;
      }
   } else {

      if (false)
         std::cout << "valid_labels(" << mtz_file_name << ","
                   << f_col << ","
                   << phi_col << ","
                   << weight_col << ","
                   << use_weights << ") returns "
                   << valid_labels(mtz_file_name, f_col, phi_col, weight_col, use_weights)
                   << std::endl;

      if (valid_labels(mtz_file_name, f_col, phi_col, weight_col, use_weights)) {

         std::vector<std::string> command_strings;
         command_strings.push_back("make-and-draw-map");
         command_strings.push_back(single_quote(mtz_file_name));
         command_strings.push_back(single_quote(f_col));
         command_strings.push_back(single_quote(phi_col));
         command_strings.push_back(single_quote(weight_col));
         command_strings.push_back(graphics_info_t::int_to_string(use_weights));
         command_strings.push_back(graphics_info_t::int_to_string(is_diff_map));
         add_to_history(command_strings);

         // std::cout << "INFO:: making map from mtz filename " << mtz_file_name << std::endl;
         logger.log(log_t::INFO, std::string("Making map from mtz filename"), mtz_file_name);
         imol = g.create_molecule();
         std::string cwd = coot::util::current_working_dir();
         g.molecules[imol].map_fill_from_mtz(std::string(mtz_file_name),
                                             cwd,
                                             f_col_str,
                                             phi_col_str,
                                             weight_col_str,
                                             use_weights, is_diff_map,
                                             graphics_info_t::map_sampling_rate);
         // save the mtz file from where the map comes
         g.molecules[imol].store_refmac_mtz_filename(std::string(mtz_file_name));
         if (! is_diff_map)
            g.scroll_wheel_map = imol;
         graphics_draw();
         g.activate_scroll_radio_button_in_display_manager(imol);

      } else {
         std::cout << "WARNING:: label(s) not found in mtz file "
                   << mtz_file_name << " " << f_col_str << " "
                   <<  phi_col_str << " ";
         if (use_weights)
            std::cout << weight_col_str << std::endl;
         else
            std::cout << std::endl;
      }
   }
   return imol; // possibly -1
}

int make_and_draw_patterson(const char *mtz_file_name,
       const char *f_col,
       const char *sigf_col) {
   graphics_info_t g;
   int imol = g.create_molecule();
   int status = g.molecules[imol].make_patterson(mtz_file_name,
    f_col, sigf_col,
    g.map_sampling_rate);

   if (! status) {
      g.erase_last_molecule();
      imol = -1;
   } else {
      graphics_draw();
   }
   return imol;
}


/* return a new molecule number or -1 on failure */
int make_and_draw_patterson_using_intensities(const char *mtz_file_name,
         const char *i_col,
         const char *sigi_col) {

   graphics_info_t g;
   int imol = g.create_molecule();
   int status = g.molecules[imol].make_patterson_using_intensities(mtz_file_name,
      i_col, sigi_col,
      g.map_sampling_rate);

   if (! status) {
      g.erase_last_molecule();
      imol = -1;
   } else {
      graphics_draw();
   }
   return imol;
}



int  make_and_draw_map_with_refmac_params(const char *mtz_file_name,
     const char *a, const char *b,
     const char *weight,
     int use_weights, int is_diff_map,
     short int have_refmac_params,
     const char *fobs_col,
     const char *sigfobs_col,
     const char *r_free_col,
     short int sensible_f_free_col) {

   graphics_info_t g;
   int imol = -1;

   // this is order dependent.  the restore-state comand that is
   // constructed in make_and_draw_map checks to see if we
   // have_refmac_params, so we need to set them before making the map.
   //
   imol = make_and_draw_map(mtz_file_name, a, b, weight, use_weights, is_diff_map);
   if (is_valid_map_molecule(imol)) {
      g.molecules[imol].store_refmac_params(std::string(mtz_file_name),
       std::string(fobs_col),
       std::string(sigfobs_col),
       std::string(r_free_col),
       sensible_f_free_col);
      g.molecules[imol].set_refmac_save_state_commands(mtz_file_name, a, b, weight, use_weights, is_diff_map, fobs_col, sigfobs_col, r_free_col, sensible_f_free_col);
   }
   return imol;
}

// return imol, possibly -1;
int make_and_draw_map_with_reso_with_refmac_params(const char *mtz_file_name,
                                                   const char *f_col,
                                                   const char *phi_col,
                                                   const char *weight_col,
                                                   int use_weights, int is_diff_map,
                                                   short int have_refmac_params,
                                                   const char *fobs_col,
                                                   const char *sigfobs_col,
                                                   const char *r_free_col,
                                                   short int sensible_r_free_col,
                                                   short int is_anomalous_flag,
                                                   short int use_reso_limits,
                                                   float low_reso_limit,
                                                   float high_reso_limit) {

   graphics_info_t g;
   int imol = -1;

   struct stat buf;
   int status = stat(mtz_file_name, &buf);

   if (status != 0) {
      std::cout << "Error finding MTZ file " << mtz_file_name << std::endl;
      if (S_ISDIR(buf.st_mode)) {
         std::cout << mtz_file_name << " is a directory! " << std::endl;
      }
   } else {
      std::string map_type;
      if (is_diff_map)
         map_type = "difference";
      else
         map_type = "conventional";

      std::string mtz_file_name_str = mtz_file_name;

      std::cout << "INFO:: making " << map_type << " map from MTZ filename "
                << mtz_file_name_str << " using " << f_col << " "
                << phi_col << std::endl;

      if (valid_labels(mtz_file_name, f_col, phi_col, weight_col, use_weights)) {
         std::string weight_col_str("");
         if (use_weights)
            weight_col_str = std::string(weight_col);
         imol = g.create_molecule();
         float msr = graphics_info_t::map_sampling_rate;
         std::string cwd = coot::util::current_working_dir();
         g.molecules[imol].map_fill_from_mtz_with_reso_limits(mtz_file_name_str,
                                                              cwd,
                                                              std::string(f_col),
                                                              std::string(phi_col),
                                                              weight_col_str,
                                                              use_weights,
                                                              is_anomalous_flag,
                                                              is_diff_map,
                                                              use_reso_limits,
                                                              low_reso_limit,
                                                              high_reso_limit,
                                                              msr);
         if (have_refmac_params) {
            g.molecules[imol].store_refmac_params(std::string(mtz_file_name),
                                                  std::string(fobs_col),
                                                  std::string(sigfobs_col),
                                                  std::string(r_free_col),
                                                  sensible_r_free_col);
            g.molecules[imol].set_refmac_save_state_commands(mtz_file_name, f_col, phi_col,
                                                             weight_col, use_weights, is_diff_map,
                                                             fobs_col, sigfobs_col,
                                                             r_free_col, sensible_r_free_col);
         } else {
            // save at least the mtz file from where the map comes
            g.molecules[imol].store_refmac_mtz_filename(std::string(mtz_file_name));
         }
         if (! is_diff_map) {
            g.scroll_wheel_map = imol;
            g.activate_scroll_radio_button_in_display_manager(imol);
         }
         graphics_draw();
      } else {
         std::cout << "WARNING:: label(s) not found in MTZ file \""
                   << mtz_file_name << "\" \"" << f_col << "\" \""
                   <<  phi_col << "\" ";
         if (use_weights)
            std::cout << "\"" << weight_col << "\"";
         std::cout << std::endl;
      }
   }
   if (imol != -1) {

      // We reset some strings if we weren't given refmac params -
      // otherwise we quote garbage or unallocated memory.
      std::string weight_col_str;
      std::string fobs_col_str;
      std::string sigfobs_col_str;
      std::string r_free_col_str;
      if (weight_col)
         weight_col_str = single_quote(weight_col);
      else
         weight_col_str = single_quote("Weight:None-specified");

      if (! have_refmac_params) {
         fobs_col_str    = single_quote("Fobs:None-specified");
         sigfobs_col_str = single_quote("SigF:None-specified");
         r_free_col_str  = single_quote("RFree:None-specified");
         sensible_r_free_col = 0;
      } else {
         fobs_col_str    = single_quote(fobs_col);
         sigfobs_col_str = single_quote(sigfobs_col_str);
         r_free_col_str  = single_quote(r_free_col);
      }

      std::vector<std::string> command_strings;
      command_strings.push_back("make-and-draw-map-with-reso-with-refmac-params");
      command_strings.push_back(single_quote(coot::util::intelligent_debackslash(mtz_file_name)));
      command_strings.push_back(single_quote(f_col));
      command_strings.push_back(single_quote(phi_col));
      command_strings.push_back(weight_col_str);
      command_strings.push_back(graphics_info_t::int_to_string(use_weights));
      command_strings.push_back(graphics_info_t::int_to_string(is_diff_map));
      command_strings.push_back(graphics_info_t::int_to_string(have_refmac_params));
      command_strings.push_back(fobs_col_str);
      command_strings.push_back(sigfobs_col_str);
      command_strings.push_back(r_free_col_str);
      command_strings.push_back(graphics_info_t::int_to_string(sensible_r_free_col));
      command_strings.push_back(graphics_info_t::int_to_string(is_anomalous_flag));
      command_strings.push_back(graphics_info_t::int_to_string(use_reso_limits));
      command_strings.push_back(graphics_info_t::float_to_string(low_reso_limit));
      command_strings.push_back(graphics_info_t::float_to_string(high_reso_limit));
      add_to_history(command_strings);
   }
   return imol;
}

int make_updating_map(const char *mtz_file_name,
		      const char *f_col, const char *phi_col,
		      const char *weight_col,
		      int use_weights, int is_diff_map) {

   int status = 1;
   int imol = make_and_draw_map(mtz_file_name, f_col, phi_col, weight_col, use_weights, is_diff_map);

   if (is_valid_map_molecule(imol)) {
      // use a better constructor?
      updating_map_params_t *ump = new updating_map_params_t(imol, mtz_file_name,
							     f_col, phi_col,
							     weight_col,
							     use_weights, is_diff_map);
      graphics_info_t::molecules[imol].continue_watching_mtz = true;
      GSourceFunc f = GSourceFunc(graphics_info_t::molecules[imol].watch_mtz);
      guint updating_map_timeout_idx = g_timeout_add(500, f, ump);
   }

   return status;
}


void stop_updating_molecule(int imol) {

   if (is_valid_map_molecule(imol)) {
      // stop watch
      graphics_info_t::molecules[imol].continue_watching_mtz = false;
   }

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].continue_watching_coordinates_file = false;
   }

}


#include "cmtz-interface.hh"

std::vector<int> auto_read_make_and_draw_maps(const char *mtz_file_name) {

   std::vector<int> imol_vec;

   if (! coot::file_exists(mtz_file_name)) {
      std::cout << "WARNING:: file " << mtz_file_name << " does not exist" << std::endl;
   } else {
      if (is_mtz_file_p(mtz_file_name)) {
         imol_vec = auto_read_make_and_draw_maps_from_mtz(mtz_file_name);
      } else {
         imol_vec = auto_read_make_and_draw_maps_from_cns(mtz_file_name);
      }
   }
   return imol_vec;
}

std::vector<int> auto_read_make_and_draw_maps_from_cns(const std::string &mtz_file_name) {

   std::vector<int> imol_vec;
   int imol1 = -1;
   int imol2 = -1;
   if (coot::util::file_name_extension(mtz_file_name) != ".mtz") {

      // Otherwise try extended CNS format
      graphics_info_t g;
      float msr = graphics_info_t::map_sampling_rate;
      imol1 = g.create_molecule();
      bool success;
      success = g.molecules[imol1].map_fill_from_cns_hkl( mtz_file_name, "F2", 0, msr );
      if (success) {
         imol_vec.push_back(imol1);
         imol2 = g.create_molecule();
         success = g.molecules[imol2].map_fill_from_cns_hkl( mtz_file_name, "F1", 1, msr );
         if (success) {
            imol_vec.push_back(imol2);
            g.scroll_wheel_map = imol1;
            g.activate_scroll_radio_button_in_display_manager(imol1);
         } else {
            g.erase_last_molecule();
         }
      } else {
         g.erase_last_molecule();
      }
   }
   return imol_vec;
}

std::vector<int> auto_read_make_and_draw_maps_from_mtz(const std::string &mtz_file_name) {

   auto read_mtz_local = [mtz_file_name] (const std::string &f, const std::string &phi,
                                          const coot::mtz_column_types_info_t &mtz_col_types,
                                          bool is_diff) {

      if (false)
         std::cout << "::::::::::::::::::::::::: read_mtz_local() " << f << " " << phi
                   << " " << is_diff << std::endl;

      std::string fobs_col;
      std::string sig_fobs_col;
      std::string r_free_col;
      short int hrp = false; // have refmac params

      // is there Fobs data in that mtz file?
      // c.f. block in auto_read_mtz() in api/molecules_container.cc

      // std::cout << "debug:: in read_mtz_local() f: \"" << f << "\" phi: \"" << phi << "\"" << std::endl;

      for (unsigned int i=0; i<mtz_col_types.f_cols.size(); i++) {
         const std::string &f = mtz_col_types.f_cols[i].column_label;
         // example f: "/2vtq/1/FP"
         std::string  nd_f = coot::util::file_name_non_directory(f);
         std::string dir_f = coot::util::file_name_directory(f);
         for (unsigned int j=0; j<mtz_col_types.sigf_cols.size(); j++) {
            const std::string &sf = mtz_col_types.sigf_cols[j].column_label;
            std::string test_string = std::string(dir_f + std::string("SIG") + nd_f);
            if (sf == test_string) {
               // OK, can we find a R-free?
               if (!mtz_col_types.r_free_cols.empty()) {
                  fobs_col = f;
                  sig_fobs_col = sf;
                  r_free_col = mtz_col_types.r_free_cols[0].column_label;
                  hrp = 1;
               }
            }
         }
      }

      short int is_anom = 0;

      // if (f == "FAN") is_anom = 1;

      int imol = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name.c_str(),
                                                                f.c_str(),
                                                                phi.c_str(),
                                                                "",
                                                                0,       //    use_weights
                                                                is_diff, // is_diff_map,
                                                                hrp,     //   short int have_refmac_params,
                                                                fobs_col.c_str(),    //   const char *fobs_col,
                                                                sig_fobs_col.c_str(),    //   const char *sigfobs_col,
                                                                r_free_col.c_str(),    //   const char *r_free_col,
                                                                hrp,   //   short int sensible_f_free_col,
                                                                is_anom, //   short int is_anomalous_flag,
                                                                0,     //   short int use_reso_limits,
                                                                0,     //   float low_reso_limit,
                                                                0);    //   float high_reso_limit

      return imol;
   };

   graphics_info_t g;

   std::vector<coot::mtz_column_trials_info_t> auto_mtz_pairs;

   // Built-ins
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("FWT",          "PHWT",      false));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("2FOFCWT",      "PH2FOFCWT", false));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("DELFWT",       "PHDELWT",   true ));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("FOFCWT",       "PHFOFCWT",  true ));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("FDM",          "PHIDM",     false));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("FAN",          "PHAN",      true));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("F_ano",        "PHI_ano",   true));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("F_early-late","PHI_early-late", true));

   for (unsigned int i=0; i<g.user_defined_auto_mtz_pairs.size(); i++)
      auto_mtz_pairs.push_back(g.user_defined_auto_mtz_pairs[i]);

   std::vector<int> imols;

   coot::mtz_column_types_info_t r = coot::get_mtz_columns(mtz_file_name);

   for (unsigned int i=0; i<auto_mtz_pairs.size(); i++) {
      const coot::mtz_column_trials_info_t &b = auto_mtz_pairs[i];
      if (valid_labels(mtz_file_name, b.f_col, b.phi_col, "", 0)) {
         int imol = read_mtz_local(b.f_col, b.phi_col, r, b.is_diff_map);
	 if (is_valid_map_molecule(imol))
	    imols.push_back(imol);
      }
   }

   // 20221001-PE if there is one F and one PHI col, read that also (and it is not a difference map)
   if (r.f_cols.size() == 1) {
      if (r.phi_cols.size() == 1) {
         int imol = read_mtz_local(r.f_cols[0].column_label, r.phi_cols[0].column_label, r, false);
         imols.push_back(imol);
      }
   }

   for (unsigned int i=0; i<r.f_cols.size(); i++) {
      std::string s = r.f_cols[i].column_label;
      std::string::size_type idx = s.find(".F_phi.F");
      if (idx != std::string::npos) {
	 std::string prefix = s.substr(0, idx);
	 std::string trial_phi_col = prefix + ".F_phi.phi";
	 for (unsigned int j=0; j<r.phi_cols.size(); j++) {
	    if (r.phi_cols[j].column_label == trial_phi_col) {
	       std::string f_col   = r.f_cols[i].column_label;
	       std::string phi_col = r.phi_cols[j].column_label;
               int imol = read_mtz_local(f_col, phi_col, r, false);
               if (is_valid_map_molecule(imol))
                  imols.push_back(imol);
	    }
	 }
      }
   }

   return imols;
}

void wrapped_auto_read_make_and_draw_maps(const char *filename) {

   auto_read_make_and_draw_maps(filename);

}

int auto_read_make_and_draw_maps_old(const char *mtz_file_name) {


   int imol1 = -1;
   int imol2 = -1;
   graphics_info_t g;

   if (! coot::file_exists(mtz_file_name)) {
      std::cout << "WARNING:: file " << mtz_file_name << " does not exist" << std::endl;
      return -1;
   }

   if ( is_mtz_file_p(mtz_file_name) ) {

      // try MTZ file
      // list of standard column names
      const char coldefs[][40] = { "+FWT,PHWT,",
                                   "-DELFWT,PHDELWT,",
                                   "+2FOFCWT,PH2FOFCWT,",
                                   "-FOFCWT,PHFOFCWT,",
                                   "+FDM,PHIDM,",
                                   "+parrot.F_phi.F,parrot.F_phi.phi,",
                                   "+pirate.F_phi.F,pirate.F_phi.phi,",
                                   "-FAN,PHAN,"};
      /* "-DELFAN,PHDELAN," not sure about this last one,
         and maybe last one could be optional!? */
      const int nc = sizeof(coldefs)/sizeof(coldefs[0]);

      // make a list of column names to try: F, phase, and difference map flag
      std::vector<std::string> cols_f(2), cols_p(2), cols_w(2), cols_d(2);

      /* no longer compile this bit
      cols_f[0] = graphics_info_t::auto_read_MTZ_FWT_col;
      cols_p[0] = graphics_info_t::auto_read_MTZ_PHWT_col;
      cols_d[0] = "+";
      cols_f[1] = graphics_info_t::auto_read_MTZ_DELFWT_col;
      cols_p[1] = graphics_info_t::auto_read_MTZ_PHDELWT_col;
      cols_d[1] = "-";
      */

      for ( int ic = 0; ic < nc; ic++ ) {
    std::string s( coldefs[ic] );
    int c1 = s.find( "," );
    int c2 = s.find( ",", c1+1 );
    std::string f = s.substr(1,c1-1);
    std::string p = s.substr(c1+1,c2-c1-1);
    std::string w = s.substr(c2+1);
    std::string d = s.substr(0,1);
    if ( f != cols_f[0] && f != cols_f[1] ) {
       cols_f.push_back(f);
       cols_p.push_back(p);
       cols_w.push_back(w);
       cols_d.push_back(d);
    }
      }

      // try each set of columns in turn
      std::vector<int> imols;
      for (unsigned int ic = 0; ic < cols_f.size(); ic++ ) {
    int imol = -1;
    int w = (cols_w[ic] != "" ) ? 1 : 0;
    int d = (cols_d[ic] != "+") ? 1 : 0;

    if ( valid_labels( mtz_file_name, cols_f[ic].c_str(),
       cols_p[ic].c_str(), cols_w[ic].c_str(), w ) )
       imol = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name,
     cols_f[ic].c_str(),
     cols_p[ic].c_str(),
     cols_w[ic].c_str(),
     w, d,  //    use_weights,  is_diff_map,
     0,     //   short int have_refmac_params,
     "",    //   const char *fobs_col,
     "",    //   const char *sigfobs_col,
     "",    //   const char *r_free_col,
     0,     //   short int sensible_f_free_col,
     0,     //   short int is_anomalous_flag,
     0,     //   short int use_reso_limits,
     0,     //   float low_reso_limit,
     0);    //   float high_reso_limit
    if ( imol >= 0 ) imols.push_back( imol );
      }

      if ( imols.size() > 0 ) {
    imol1 = imols.front();
    imol2 = imols.back();
    g.scroll_wheel_map = imol1;
    g.activate_scroll_radio_button_in_display_manager(imol1);
      } else {
    GtkWidget *w = wrapped_nothing_bad_dialog("Failed to find any suitable F/phi columns in the MTZ file");
    gtk_widget_set_visible(w, TRUE);
      }

   } else {

      // don't try to read as a CNS file if it is called an .mtz file
      // (map_fill_from_cns_hkl can dump core if file is not properly
      // formatted).
      //
      if (coot::util::file_name_extension(mtz_file_name) != ".mtz") {

    // Otherwise try extended CNS format
    float msr = graphics_info_t::map_sampling_rate;
    imol1 = g.create_molecule();
    bool success;
    success = g.molecules[imol1].map_fill_from_cns_hkl( mtz_file_name, "F2", 0, msr );
    if (success) {
       imol2 = g.create_molecule();
       success = g.molecules[imol2].map_fill_from_cns_hkl( mtz_file_name, "F1", 1, msr );
       if (success) {
          g.scroll_wheel_map = imol1;
          g.activate_scroll_radio_button_in_display_manager(imol1);
       } else {
          g.erase_last_molecule();
       }
    } else {
       g.erase_last_molecule();
    }
      }
   }
   return imol2;
}

int auto_read_do_difference_map_too_state() {

   add_to_history_simple("auto-read-do-difference-map-too-state");
   int i = graphics_info_t::auto_read_do_difference_map_too_flag;
   return i;

}
void set_auto_read_do_difference_map_too(int i) {

   graphics_info_t::auto_read_do_difference_map_too_flag = i;
   std::string cmd = "set-auto-read-do-dfference-map-too";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
}

/*! \brief set the default inital contour for 2FoFc-style map

in sigma */
void set_default_initial_contour_level_for_map(float n_sigma) {

   graphics_info_t::default_sigma_level_for_map = n_sigma;
   std::string cmd = "set-default-initial-contour-level-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(n_sigma);
   add_to_history_typed(cmd, args);
}

/*! \brief set the default inital contour for FoFc-style map

in sigma */
void set_default_initial_contour_level_for_difference_map(float n_sigma) {

   graphics_info_t::default_sigma_level_for_fofc_map = n_sigma;
   std::string cmd = "set-default-initial-contour-level-for-difference-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(n_sigma);
   add_to_history_typed(cmd, args);
}


/*! \brief by default, maps that are P1 and have 90 degree angles
           are considered as maps without symmetry (i.e. EM maps).
           In some cases though P1 maps do/should have symmetry -
           and this is the means by you can tell Coot that.
    @param state 1 turns on map symmetry
*/
void set_map_has_symmetry(int imol, int state) {

   if (is_valid_map_molecule(imol)) {
      bool is_em_map = true;
      if (state) is_em_map = false;
      // graphics_info_t::molecules[imol].is_em_map_cached_flag = is_em_map;
      graphics_info_t::molecules[imol].set_map_has_symmetry(is_em_map);
   }

}





void set_map_line_width(int w) {
   graphics_info_t::map_line_width = w;
   // update the maps because they may be being draw as graphical
   // objects.
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++)
      graphics_info_t::molecules[imol].update_map(false);
   graphics_draw();
   std::string cmd = "set-map-line-width";
   std::vector<coot::command_arg_t> args;
   args.push_back(w);
   add_to_history_typed(cmd, args);

}

int map_line_width_state() {
   add_to_history_simple("map-line-width-state");
   return graphics_info_t::map_line_width;
}

int swap_difference_map_colours_state() {
  int ret = graphics_info_t::swap_difference_map_colours;
  return ret;
}

/* return success status 0 = failure (imol does not have a map) */
int set_map_is_difference_map(int imol, short int  bool_flag) {

   int istatus = 0;
   if (imol< graphics_n_molecules()) {
      if (graphics_info_t::molecules[imol].has_xmap()) {
	 graphics_info_t::molecules[imol].set_map_is_difference_map(bool_flag);
	 istatus = 1;
	 graphics_draw();
      } else {
         std::cout << "WARNING:: molecule " << imol << " does not have a map." <<  std::endl;
      }

   } else {
      std::cout << "WARNING:: No such molecule as " << imol << std::endl;
   }
   std::string cmd = "set-map-is-difference-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);

   return istatus;
}

int map_is_difference_map(int imol) {

   int istat = 0;
   if (is_valid_map_molecule(imol)) {
      istat = graphics_info_t::molecules[imol].is_difference_map_p();
   }
   return istat;
}


/* return the index of the new molecule or -1 on failure */
int another_level() {

   int istat = -1;
   int imap = -1;

   imap = imol_refinement_map();
   if (imap == -1) {
      for (int i=0; i<graphics_info_t::n_molecules(); i++) {
    if (is_valid_map_molecule(i)) {
       if (! graphics_info_t::molecules[i].is_difference_map_p()) {
          imap = i;
       }
    }
      }
   }

   if (imap > -1) {
      istat = another_level_from_map_molecule_number(imap);
   }
   // history elsewhere
   return istat;
}

int another_level_from_map_molecule_number(int imap) {
   int istat = -1;
   if (is_valid_map_molecule(imap)) {
      // create another map with the same parameters as imap and then
      // push up the contour level a sigma.
//       std::cout << "DEBUG:: calling make_and_draw_map_with_reso_with_refmac_params"
// 		<< std::endl;
      istat = make_and_draw_map_with_reso_with_refmac_params(
     graphics_info_t::molecules[imap].save_mtz_file_name.c_str(),
     graphics_info_t::molecules[imap].save_f_col.c_str(),
     graphics_info_t::molecules[imap].save_phi_col.c_str(),
          graphics_info_t::molecules[imap].save_weight_col.c_str(),
          graphics_info_t::molecules[imap].save_use_weights,
          graphics_info_t::molecules[imap].save_is_diff_map_flag,
     0, "None", "None", "None", 0, // refmac params
          graphics_info_t::molecules[imap].save_is_anomalous_map_flag,
          graphics_info_t::molecules[imap].save_use_reso_limits,
          graphics_info_t::molecules[imap].save_low_reso_limit,
          graphics_info_t::molecules[imap].save_high_reso_limit);

//       std::cout << "DEBUG:: istat in another_level_from_map_molecule_number "
// 		<< istat << std::endl;

      if (istat != -1) {

	 float map_sigma = graphics_info_t::molecules[istat].map_sigma();
	 float current_contour_level = graphics_info_t::molecules[istat].contour_level;
	 graphics_info_t::molecules[istat].set_contour_level(current_contour_level +
							     map_sigma*1.0);
	 graphics_info_t::molecules[istat].update_map(true);
	 graphics_draw();
      }
   }
   std::string cmd = "another-level-from-map-molecule-number";
   std::vector<coot::command_arg_t> args;
   args.push_back(imap);
   add_to_history_typed(cmd, args);
   return istat;
}

void set_map_radius_slider_max(float f) {
   graphics_info_t::map_radius_slider_max = f;
   std::string cmd = "set-map-radius-slider-max";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
}

void set_contour_level_absolute(int imol_map, float level) {

   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t::molecules[imol_map].set_contour_level(level);
   }
   graphics_draw();

   std::string cmd = "set-contour-level-absolute";
   std::vector<coot::command_arg_t> args;
   args.push_back(level);
   add_to_history_typed(cmd, args);
}

void set_contour_level_in_sigma(int imol_map, float level) {

   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t::molecules[imol_map].set_contour_level_by_sigma(level);
   }
   graphics_draw();
   std::string cmd = "set-contour-level-in-sigma";
   std::vector<coot::command_arg_t> args;
   args.push_back(level);
   add_to_history_typed(cmd, args);
}

/* \brief get the contour level */
float get_contour_level_absolute(int imol) {

   float r = 0;

   if (is_valid_map_molecule(imol)) {
      r = graphics_info_t::molecules[imol].contour_level;
   }
   return r;
}

/* \brief get the contour level in rmd above 0. */
float get_contour_level_in_sigma(int imol) {

   float r = 0;

   if (is_valid_map_molecule(imol)) {
      double s = graphics_info_t::molecules[imol].map_sigma();
      r = graphics_info_t::molecules[imol].contour_level/s;
   }
   return r;
}



void set_last_map_contour_level(float level) {

   graphics_info_t g;
   g.set_last_map_contour_level(level);
   std::string cmd = "set-last-map-contour-level";
   std::vector<coot::command_arg_t> args;
   args.push_back(level);
   add_to_history_typed(cmd, args);
}

void set_last_map_contour_level_by_sigma(float level) {

   graphics_info_t g;
   g.set_last_map_contour_level_by_sigma(level);
   std::string cmd = "set-last-map-contour-level-by-sigma";
   std::vector<coot::command_arg_t> args;
   args.push_back(level);
   add_to_history_typed(cmd, args);
}

void set_last_map_sigma_step(float f) {

   graphics_info_t g;
   g.set_last_map_sigma_step(f);
   std::string cmd = "set-last-map-sigma-step";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);

}

// -------------------------------------------------------------------------
//                        (density) iso level increment entry
// -------------------------------------------------------------------------
//

// imol is ignored.
//
char* get_text_for_iso_level_increment_entry(int imol) {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-6.4f", g.iso_level_increment);

   return text;

}

void set_iso_level_increment_from_text(const char *text, int imol) {

   float val;

   graphics_info_t g;

   val = atof(text);

   if ((val > 10000) || (val < -10000)) {
      std::cout << "Cannot interpret: " << text
                << ".  Assuming 0.05 for increment" << std::endl;
      val  = 0.05;

   }

   std::cout << "setting iso_level_increment to " << val << std::endl;
   g.iso_level_increment = val;

   graphics_draw();
}

void set_iso_level_increment(float val) {
   graphics_info_t g;
   g.iso_level_increment = val;
}

float get_iso_level_increment() {
  float ret = graphics_info_t:: iso_level_increment;
  return ret;
}

// imol is ignored.
//
char* get_text_for_diff_map_iso_level_increment_entry(int imol) {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-6.4f", g.diff_map_iso_level_increment);
   return text;

}

void set_diff_map_iso_level_increment_from_text(const char *text, int imol) {

   float val;
   graphics_info_t g;

   val = atof(text);

   if ((val > 10000) || (val < -10000)) {
      std::cout << "Cannot interpret: " << text
                << ".  Assuming 0.005 for increment" << std::endl;
      val  = 0.005;
   }
   g.diff_map_iso_level_increment = val;
   graphics_draw();
}

void set_diff_map_iso_level_increment(float val) {
   graphics_info_t::diff_map_iso_level_increment = val;
}

float get_diff_map_iso_level_increment() {
   float ret =graphics_info_t::diff_map_iso_level_increment;
   return ret;
}

void set_map_sampling_rate_text(const char *text) {

   float val;
   val = atof(text);

   if ((val > 100) || (val < 1)) {
      std::cout << "Nonsense value: " << text
                << ".  Assuming 1.5 for increment" << std::endl;
      val  = 1.5;
   }
   set_map_sampling_rate(val);

}

void set_map_sampling_rate(float r) {

   graphics_info_t g;
   g.map_sampling_rate = r;

}

char* get_text_for_map_sampling_rate_text() {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-5.4f", g.map_sampling_rate);
   return text;


}

float get_map_sampling_rate() {
   graphics_info_t g;
   return g.map_sampling_rate;
}


/* applies to the current map */
void change_contour_level(short int is_increment) { // else is decrement.

   graphics_info_t g;
   int s = g.scroll_wheel_map;

   if (is_valid_map_molecule(s)) {

      if (g.molecules[s].is_difference_map_p()) {
    g.molecules[s].contour_level += g.diff_map_iso_level_increment;
      } else {
    // normal case
    if (is_increment) {
       g.molecules[s].contour_level += g.iso_level_increment;
    } else {
       g.molecules[s].contour_level -= g.iso_level_increment;
    }
      }
      g.molecules[s].update_map(true);
      graphics_draw();
   }
}

void set_initial_map_for_skeletonize() {

   graphics_info_t::set_initial_map_for_skeletonize();

}

void set_max_skeleton_search_depth(int v) {
   graphics_info_t g;
   g.set_max_skeleton_search_depth(v);
}

/* Set the radio buttons in the frame to the be on or off for the map
   that is displayed in the optionmenu (those menu items "activate"
   callbacks (graphics_info::skeleton_map_select change
   g.map_for_skeletonize).  */
void set_on_off_skeleton_radio_buttons(GtkWidget *skeleton_frame) {

   graphics_info_t g;
   g.set_on_off_skeleton_radio_buttons(skeleton_frame);
}

// imol is used imap is ignored.
// You can fix this anachronism one day if you like.  FIXME.
//
void set_scrollable_map(int imol) {

   graphics_info_t g;
   if (is_valid_map_molecule(imol)) {
      g.set_scrollable_map(imol); // in graphics-info.h
   } else {
      std::cout << "WARNING:: " << imol << " is not a valid molecule" << " in set_scrollable_map\n";
   }
}

int
skeletonize_map(int imol, short int prune_flag) {

   graphics_info_t::skeletonize_map(imol, prune_flag);
   return 0;
}

int unskeletonize_map(int imol) {
   graphics_info_t::unskeletonize_map(imol);
   return imol;
}


void
do_skeleton_prune() {

   graphics_info_t g;
   float map_cutoff  = g.skeleton_level;

   for (int imol=0; imol<g.n_molecules(); imol++) {
      if (g.molecules[imol].has_xmap() &&
     !g.molecules[imol].xmap_is_diff_map) {

    if (g.molecules[imol].xskel_is_filled == 1) {

       BuildCas bc(g.molecules[imol].xmap, map_cutoff);

       // mark segments by connectivity
       //
       int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan,
          g.molecules[imol].xmap,
          map_cutoff);

       bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
       g.molecules[imol].update_clipper_skeleton();
    }
      }
   }
}
#ifdef USE_GUILE
SCM get_map_colour_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_map_molecule(imol)) {
      // std::pair<GdkRGBA, GdkRGBA> colours map_colours();
      // std::vector<float> colour_v = graphics_info_t::molecules[imol].map_colours();
      std::pair<GdkRGBA, GdkRGBA> colours = graphics_info_t::molecules[imol].get_map_colours();

      std::cout << "get_map_colour_scm() needs fixing " << std::endl; // FIXME
#if 0
      if (colour_v.size() > 2) {
         r = scm_list_3(scm_from_double(colour_v[0]),
                        scm_from_double(colour_v[1]),
                        scm_from_double(colour_v[2]));
      }
#endif
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *get_map_colour_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_map_molecule(imol)) {
      std::pair<GdkRGBA, GdkRGBA> colours = graphics_info_t::molecules[imol].get_map_colours();
      r = PyList_New(2);
      PyObject *col_1 = PyList_New(3);
      PyObject *col_2 = PyList_New(3);
      PyList_SetItem(col_1, 0, PyFloat_FromDouble(colours.first.red));
      PyList_SetItem(col_1, 1, PyFloat_FromDouble(colours.first.green));
      PyList_SetItem(col_1, 2, PyFloat_FromDouble(colours.first.blue));
      PyList_SetItem(col_2, 0, PyFloat_FromDouble(colours.second.red));
      PyList_SetItem(col_2, 1, PyFloat_FromDouble(colours.second.green));
      PyList_SetItem(col_2, 2, PyFloat_FromDouble(colours.second.blue));
      PyList_SetItem(r, 0, col_1);
      PyList_SetItem(r, 1, col_2);
   }
   if (PyBool_Check(r)) {
      Py_XINCREF(r);
   }
   return r;
}
#endif



/* give a warning dialog if density it too dark (blue) */
void check_for_dark_blue_density() {

   if (graphics_info_t::use_graphics_interface_flag) {
      for (int i=0; i<graphics_info_t::n_molecules(); i++) {
         if (graphics_info_t::molecules[i].has_xmap()) {
            if (graphics_info_t::molecules[i].is_displayed_p()) {
               if (background_is_black_p()) {
                  if (graphics_info_t::molecules[i].map_is_too_blue_p()) {
                     std::string s = "I suggest that you increase the brightness of the map\n";
                     s += " if this is for a presentation (blue projects badly).";
                     info_dialog(s.c_str());
                     break; // only make the dialog once
                  }
               }
            }
         }
      }
   }
}


#include "utils/coot-utils.hh"

int handle_read_emdb_data(const std::string &dir_name) {

   int status = 0;
   std::string map_dir = coot::util::append_dir_dir(dir_name, "map");
   std::string pdb_dir = coot::util::append_dir_dir(coot::util::append_dir_dir(dir_name, "fittedModels"), "PDB");
   std::vector<std::string> map_files = coot::util::glob_files(map_dir, "*.map");
   std::vector<std::string> pdb_files = coot::util::glob_files(pdb_dir, "*.ent");
   for (auto map_file : map_files)
      handle_read_ccp4_map(map_file, 0);
   for (auto pdb_file : pdb_files)
      read_pdb(pdb_file);

   return status;
}

void set_contour_by_sigma_step_by_mol(int imol, float f, short int state) {

   if (imol < graphics_info_t::n_molecules()) {
      if (imol >= 0) {
         if (graphics_info_t::molecules[imol].has_xmap()) {
            // NXMAP-FIXME
            graphics_info_t::molecules[imol].set_contour_by_sigma_step(f, state);
         }
      }
   }
}

int export_map(int imol, const char *filename) {

   int rv = 0; // fail
   if (is_valid_map_molecule(imol)) {
      try {
         clipper::CCP4MAPfile mapout;
         mapout.open_write(std::string(filename));
         mapout.export_xmap(graphics_info_t::molecules[imol].xmap);
         mapout.close_write();
         rv = 1;
      }
      catch (...) {
         std::cout << "WARNING:: CCP4 map writing error for " << filename << std::endl;
      }

   } else {
      graphics_info_t g;
      g.add_status_bar_text("WARNING:: Invalid map molecule number");
   }
   return rv;
}

int export_map_fragment(int imol, float x, float y, float z, float radius, const char *filename) {

   int rv = 0;
   if (is_valid_map_molecule(imol)) {
      graphics_info_t g;
      clipper::Coord_orth pos(x,y,z);
      g.molecules[imol].export_map_fragment(radius, pos, filename);
      rv = 1;
   }
   return rv;
}

int export_map_fragment_to_plain_file(int imol, float x, float y, float z, float radius, const char *file_name) {

   int rv = 0;
   if (is_valid_map_molecule(imol)) {
      graphics_info_t g;
      clipper::Coord_orth pos(x,y,z);
      g.molecules[imol].export_map_fragment_to_plain_file(radius, pos, file_name);
      rv = 1;
   }
   return rv;
}

/*! convenience function, called from callbacks.c */
void export_map_fragment_with_text_radius(int imol, const char *radius_text, const char *filename) {

   graphics_info_t g;
   coot::Cartesian rc = g.RotationCentre();
   float radius = coot::util::string_to_int(radius_text);
   export_map_fragment(imol, rc.x(), rc.y(), rc.z(), radius, filename);
}


/*! \brief export a fragment of the map about (x,y,z)  */
int export_map_fragment_with_origin_shift(int imol, float x, float y, float z, float radius, const char *filename) {

   int rv = 0;
   if (is_valid_map_molecule(imol)) {
      graphics_info_t g;
      clipper::Coord_orth pos(x,y,z);
      g.molecules[imol].export_map_fragment_with_origin_shift(radius, pos, filename);
      rv = 1;
   }
   return rv;
}

void set_ignore_pseudo_zeros_for_map_stats(short int state) {

   graphics_info_t::ignore_pseudo_zeros_for_map_stats = state;

}


void map_histogram(int imol_map) {

#if HAVE_GOOCANVAS
   if (graphics_info_t::use_graphics_interface_flag) {
      if (is_valid_map_molecule(imol_map)) {

         bool map_histograms_are_revealers = false;

         bool ignore_pseudo_zeros = false;

	 const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
	 unsigned int n_bins = 400;
         bool write_output_flag = false;
	 mean_and_variance<float> mv = map_density_distribution(xmap, n_bins, write_output_flag, ignore_pseudo_zeros);

         if (false) { // debug bin data
            std::cout << "mv data " << imol_map << std::endl;
            for (unsigned int i=0; i<mv.bins.size(); i++) {
               float constr_x = (static_cast<float>(i) + 0.5) * mv.bin_width + mv.min_density;
               std::cout << i << " constr-x " << constr_x << " " << mv.bins[i] << std::endl;
            }
         }

	 if (mv.bins.size() > 0) {
	    std::vector<std::pair<double, double> > data(mv.bins.size());
	    for (unsigned int ibin=0; ibin<mv.bins.size(); ibin++) {
	       double x = (ibin+0.5)*mv.bin_width + mv.min_density;
	       double y = mv.bins[ibin];
	       data[ibin] = std::pair<double, double> (x, y);
	    }

            const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
            unsigned int n_bins = 100;
            if (ignore_pseudo_zeros) {
               n_bins = 400;
            }
            mean_and_variance<float> mv = map_density_distribution(xmap, n_bins, false, ignore_pseudo_zeros);

            if (mv.bins.size() > 0) {
               std::vector<std::pair<double, double> > data(mv.bins.size());
               for (unsigned int ibin=0; ibin<mv.bins.size(); ibin++) {
                  double x = (ibin+0.5)*mv.bin_width + mv.min_density;
                  double y = mv.bins[ibin];
                  data[ibin] = std::pair<double, double> (x, y);
               }

               coot::goograph* g = new coot::goograph;
               int trace = g->trace_new();

               std::string title = "Density Histogram for map " + std::to_string(imol_map);
               g->set_plot_title(title);
               g->set_data(trace, data);
               g->set_axis_label(coot::goograph::X_AXIS, "Density Value");
               g->set_axis_label(coot::goograph::Y_AXIS, "Counts");
               g->set_trace_type(trace, coot::graph_trace_info_t::PLOT_TYPE_BAR);
               if (ignore_pseudo_zeros) {
                  float x_range_min = mv.mean-3.0*sqrt(mv.variance);
                  float x_range_max = mv.mean+3.0*sqrt(mv.variance);
                  g->set_extents(coot::goograph::X_AXIS, x_range_min, x_range_max);
                  g->set_extents(coot::goograph::Y_AXIS, 0, mv.histogram_max);
                  std::cout << "::::: set_extents() X: " << x_range_min << " " << x_range_max << "\n";
                  std::cout << "::::: set_extents() Y: " << mv.histogram_max << "\n";
               }
               g->show_dialog();
            }
         }
      }
   }
# endif // HAVE_GOOCANVAS
}



/* create a number of maps by segmenting the given map, above the
   (absolute) low_level.  New maps are on the same grid as the input
   map.  */
void segment_map(int imol_map, float low_level) {

   int max_segments = 300;

   if (is_valid_map_molecule(imol_map)) {
      const clipper::Xmap<float> &xmap_in = graphics_info_t::molecules[imol_map].xmap;
      coot::util::segment_map s;
      std::pair<int, clipper::Xmap<int> > segmented_map = s.segment(xmap_in, low_level);
      float contour_level = graphics_info_t::molecules[imol_map].get_contour_level();

      for (int iseg=0; (iseg<segmented_map.first) && iseg<max_segments; iseg++) {
	 // std::cout << "iseg: " << iseg << std::endl;
	 clipper::Xmap<float> xmap;
	 xmap.init(segmented_map.second.spacegroup(),
		   segmented_map.second.cell(),
		   segmented_map.second.grid_sampling());
	 clipper::Xmap_base::Map_reference_index ix;
	 for (ix = segmented_map.second.first(); !ix.last(); ix.next()) {
	    if (segmented_map.second[ix] == iseg)
	       xmap[ix] = xmap_in[ix];
	 }
	 int imol_new = graphics_info_t::create_molecule();
	 std::string name = "Map ";
	 name += coot::util::int_to_string(imol_map);
	 name += " Segment ";
	 name += coot::util::int_to_string(iseg);
	 bool is_em_map_flag = graphics_info_t::molecules[imol_map].is_EM_map();
	 graphics_info_t::molecules[imol_new].install_new_map(xmap, name, is_em_map_flag);
	 graphics_info_t::molecules[imol_new].set_contour_level(contour_level);
      }
   }
   graphics_draw();
}

void segment_map_multi_scale(int imol_map, float low_level, float b_factor_inc, int n_rounds) {

   int max_segments = 8;
   if (is_valid_map_molecule(imol_map)) {
      clipper::Xmap<float> &xmap_in = graphics_info_t::molecules[imol_map].xmap;
      coot::util::segment_map s;
      std::pair<int, clipper::Xmap<int> > segmented_map = s.segment(xmap_in, low_level, b_factor_inc, n_rounds);
      float contour_level = graphics_info_t::molecules[imol_map].get_contour_level();
      for (int iseg=0; (iseg<segmented_map.first) && iseg<max_segments; iseg++) {
	 clipper::Xmap<float> xmap;
	 xmap.init(segmented_map.second.spacegroup(),
		   segmented_map.second.cell(),
		   segmented_map.second.grid_sampling());
	 clipper::Xmap_base::Map_reference_index ix;
	 int n_points_in_map = 0;
	 for (ix = segmented_map.second.first(); !ix.last(); ix.next()) {
	    if (segmented_map.second[ix] == iseg) {
	       float f = xmap_in[ix];
	       xmap[ix] = f;
	       n_points_in_map++;
	    }
	 }
	 if (n_points_in_map) {
	    int imol_new = graphics_info_t::create_molecule();
	    std::string name = "Map ";
	    name += coot::util::int_to_string(imol_map);
	    name += " Segment ";
	    name += coot::util::int_to_string(iseg);
	    bool is_em_flag = graphics_info_t::molecules[imol_map].is_EM_map();
	    graphics_info_t::molecules[imol_new].install_new_map(xmap, name, is_em_flag);
	    graphics_info_t::molecules[imol_new].set_contour_level(contour_level);
	 }
      }
   }
   graphics_draw();
}



// the cell is given in Angstroms and the angles in degrees.
int transform_map_raw(int imol,
                      double r00, double r01, double r02,
                      double r10, double r11, double r12,
                      double r20, double r21, double r22,
                      double t0, double t1, double t2,
                      double pt1, double pt2, double pt3, double box_size,
                      const char *ref_space_group,
                      double cell_a, double cell_b, double cell_c,
                      double alpha, double beta, double gamma) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol)) {
      clipper::Mat33<double> m(r00, r01, r02, r10, r11, r12, r20, r21, r22);
      clipper::Coord_orth c(t0, t1, t2);
      clipper::RTop_orth rtop(m,c);
      // clipper::RTop_orth rtop_inv = rtop.inverse();
      clipper::Coord_orth pt(pt1, pt2, pt3);

      std::cout << "INFO:: in transforming map around target point "
                << pt.format() << std::endl;

      clipper::Spgr_descr sg_descr(ref_space_group);
      clipper::Spacegroup new_space_group(sg_descr);

      clipper::Cell_descr cell_d(cell_a, cell_b, cell_c,
                                 clipper::Util::d2rad(alpha),
                                 clipper::Util::d2rad(beta),
                                 clipper::Util::d2rad(gamma));
      clipper::Cell new_cell(cell_d);


      clipper::Xmap<float> new_map =
         coot::util::transform_map(graphics_info_t::molecules[imol].xmap,
                                   new_space_group, new_cell,
                                   rtop, pt, box_size);

      const coot::ghost_molecule_display_t ghost_info;
      // int is_diff_map_flag = graphics_info_t::molecules[imol].is_difference_map_p();
      // int swap_colours_flag = graphics_info_t::swap_difference_map_colours;
      bool ipz = graphics_info_t::ignore_pseudo_zeros_for_map_stats;
      mean_and_variance<float> mv = map_density_distribution(new_map, 40, false, ipz);
      std::string name = "Transformed map";
      imol_new = graphics_info_t::create_molecule();
      bool is_em_flag = graphics_info_t::molecules[imol].is_EM_map();
      graphics_info_t::molecules[imol_new].install_new_map(new_map, name, is_em_flag);
      graphics_draw();

   } else {
      std::cout << "WARNING:: molecule " << imol << " is not a valid map" << std::endl;
   }
   return imol_new;
}

/*! \brief make a difference map, taking map_scale * imap2 from imap1,
  on the grid of imap1.  Return the new molecule number.
  Return -1 on failure. */
int difference_map(int imol1, int imol2, float map_scale) {

   int r = -1;

   if (is_valid_map_molecule(imol1)) {
      if (is_valid_map_molecule(imol2)) {
 	 std::pair<clipper::Xmap<float>, float> dm =
 	    coot::util::difference_map(graphics_info_t::molecules[imol1].xmap,
 				       graphics_info_t::molecules[imol2].xmap,
 				       map_scale);
	 int imol = graphics_info_t::create_molecule();
	 std::string name = "difference-map";
	 // int swpcolf = graphics_info_t::swap_difference_map_colours;
	 bool is_em_flag = graphics_info_t::molecules[imol1].is_EM_map();
	 graphics_info_t::molecules[imol].install_new_map(dm.first, name, is_em_flag);
	 graphics_info_t::molecules[imol].set_map_is_difference_map(true);

	 r = imol;
	 graphics_draw();
      }
   }
   return r;
}


int reinterp_map(int map_no, int reference_map_no) {

   int r = -1;
   if (is_valid_map_molecule(map_no)) {
      if (is_valid_map_molecule(reference_map_no)) {
	 graphics_info_t g;
	 clipper::Xmap<float> new_map =
	    coot::util::reinterp_map(g.molecules[map_no].xmap,
				     g.molecules[reference_map_no].xmap);
	 int imol = graphics_info_t::create_molecule();
	 std::string name = "map ";
	 name += coot::util::int_to_string(map_no);
	 name += " re-interprolated to match ";
	 name += coot::util::int_to_string(reference_map_no);
	 bool is_em_flag = graphics_info_t::molecules[map_no].is_EM_map();
	 graphics_info_t::molecules[imol].install_new_map(new_map, name, is_em_flag);
	 r = imol;
	 graphics_draw();
      }
   }
   return r;
}

/*! \brief make a new map (a copy of map_no) that is in the cell,
  spacegroup and a multiple of the sampling of the input map (a
  sampling factor of more than 1 makes the output maps smoother) */
int smooth_map(int map_no, float sampling_multiplier) {

   int r = -1;

   if (is_valid_map_molecule(map_no)) {
      graphics_info_t g;
      clipper::Xmap<float> new_map =
         coot::util::reinterp_map(g.molecules[map_no].xmap, sampling_multiplier);
      int imol = graphics_info_t::create_molecule();
      std::string name = "map ";
      name += coot::util::int_to_string(map_no);
      name += " re-interprolated by factor ";
      name += coot::util::float_to_string(sampling_multiplier);
      bool is_em_flag = graphics_info_t::molecules[map_no].is_EM_map();
      graphics_info_t::molecules[imol].install_new_map(new_map, name, is_em_flag);
      r = imol;
      graphics_draw();
   }
   return r;
}


#ifdef USE_GUILE
/*! \brief make an average map from the map_number_and_scales (which
  is a list of pairs (list map-number scale-factor)) (the scale factor
  are typically 1.0 of course). */
int average_map_scm(SCM map_number_and_scales) {

   int r = -1;
   SCM n_scm = scm_length(map_number_and_scales);
   int n = scm_to_int(n_scm);
   bool is_em_flag = false;
   std::vector<std::pair<clipper::Xmap<float>, float> > maps_and_scales_vec;
   for (int i=0; i<n; i++) {
      SCM number_and_scale = scm_list_ref(map_number_and_scales, scm_from_int(i));
      SCM ns_scm = scm_length(number_and_scale);
      int ns = scm_to_int(ns_scm);
      if (ns == 2) {
	 SCM map_number_scm = scm_list_ref(number_and_scale, scm_from_int(0));
	 SCM map_scale_scm  = scm_list_ref(number_and_scale, scm_from_int(1));
	 if (scm_is_true(scm_integer_p(map_number_scm))) {
	    if (scm_is_true(scm_number_p(map_scale_scm))) {
	       int map_number = scm_to_int(map_number_scm);
	       if (is_valid_map_molecule(map_number)) {
	          float scale = scm_to_double(map_scale_scm);
	          std::pair<clipper::Xmap<float>, float> p(graphics_info_t::molecules[map_number].xmap, scale);
	          maps_and_scales_vec.push_back(p);
	          is_em_flag = graphics_info_t::molecules[map_number].is_EM_map();
	       } else {
	          std::cout << "Invalid map number " << map_number << std::endl;
	       }
	    } else {
	       std::cout << "Bad scale "
			 << scm_to_locale_string(display_scm(map_scale_scm))
			 << " ignoring map "
			 << scm_to_locale_string(display_scm(map_number_scm))
			 << std::endl;
	    }
	 } else {
	    std::cout << "Bad map number " << scm_to_locale_string(display_scm(map_number_scm))
	              << std::endl;
	 }

      }
   }
   if (maps_and_scales_vec.size() > 0) {
      clipper::Xmap<float> average_map = coot::util::average_map(maps_and_scales_vec);
      int imol = graphics_info_t::create_molecule();
      std::string name = "averaged-map";
      graphics_info_t::molecules[imol].install_new_map(average_map, name, is_em_flag);
      r = imol;
      graphics_draw();
   }
   return r;
}
#endif


#ifdef USE_PYTHON
/*! \brief make an average map from the map_number_and_scales (which
  is a list of pairs [map-number, scale-factor]) (the scale factors
  are typically 1.0 of course). */
int average_map_py(PyObject *map_number_and_scales) {

   int r = -1;
   int n = PyObject_Length(map_number_and_scales);
   std::vector<std::pair<clipper::Xmap<float>, float> > maps_and_scales_vec;
   bool is_em_flag = false;
   for (int i=0; i<n; i++) {
      PyObject *number_and_scale = PyList_GetItem(map_number_and_scales, i);
      int ns = PyObject_Length(number_and_scale);
      if (ns == 2) {
         PyObject *map_number_py = PyList_GetItem(number_and_scale, 0);
         PyObject *map_scale_py  = PyList_GetItem(number_and_scale, 1);
         if (PyLong_Check(map_number_py)) {
           if (PyFloat_Check(map_scale_py) || PyLong_Check(map_scale_py)) {
               int map_number = PyLong_AsLong(map_number_py);
               if (is_valid_map_molecule(map_number)) {
                  float scale = PyFloat_AsDouble(map_scale_py);
                  std::pair<clipper::Xmap<float>, float> p(graphics_info_t::molecules[map_number].xmap, scale);
                  maps_and_scales_vec.push_back(p);
                  is_em_flag = graphics_info_t::molecules[map_number].is_EM_map();
               } else {
                  std::cout << "Invalid map number " << map_number << std::endl;
               }
            } else {
             std::cout << "Bad scale " << PyUnicode_AsUTF8String(display_python(map_scale_py))   // FIXME
                              << std::endl;
            }
         } else {
           std::cout << "Bad map number " << PyUnicode_AsUTF8String(display_python(map_number_py))  // FIXME
                 << std::endl;
         }
      }
   }
   if (maps_and_scales_vec.size() > 0) {
      clipper::Xmap<float> average_map = coot::util::average_map(maps_and_scales_vec);
      int imol = graphics_info_t::create_molecule();
      std::string name = "averaged-map";
      graphics_info_t::molecules[imol].install_new_map(average_map, name, is_em_flag);
      r = imol;
      graphics_draw();
   }
   return r;
}
#endif



#ifdef USE_PYTHON
/*! \brief Somewhat similar to the above function, except in this
case we overwrite the imol_map and we also presume that the
grid sampling of the contributing maps match. This makes it
much faster to generate than an average map.
*/
void regen_map_py(int imol_map, PyObject *map_number_and_scales) {

   auto pyobject_to_map_index_and_scale_vec = [] (PyObject *map_number_and_scales) {
      std::vector<std::pair<int, float> > map_indices_and_scales_vec;
      int n = PyObject_Length(map_number_and_scales);
      for (int i=0; i<n; i++) {
         PyObject *number_and_scale = PyList_GetItem(map_number_and_scales, i);
         int ns = PyObject_Length(number_and_scale);
         if (ns == 2) {
            PyObject *map_number_py = PyList_GetItem(number_and_scale, 0);
            PyObject *map_scale_py  = PyList_GetItem(number_and_scale, 1);
            if (PyLong_Check(map_number_py)) {
               if (PyFloat_Check(map_scale_py) || PyLong_Check(map_scale_py)) {
                  int map_number = PyLong_AsLong(map_number_py);
                  if (is_valid_map_molecule(map_number)) {
                     float scale = PyFloat_AsDouble(map_scale_py);
                     std::pair<int, float> p(map_number, scale);
                     map_indices_and_scales_vec.push_back(p);
                  } else {
                     std::cout << "Invalid map number " << map_number << std::endl;
                  }
               } else {
                  std::cout << "Bad scale " << PyUnicode_AsUTF8String(display_python(map_scale_py))   // FIXME
                            << std::endl;
               }
            } else {
               std::cout << "Bad map number " << PyUnicode_AsUTF8String(display_python(map_number_py))  // FIXME
                         << std::endl;
            }
         }
      }
      return map_indices_and_scales_vec;
   };

   bool status = false;
   if (is_valid_map_molecule(imol_map)) {

      std::vector<std::pair<int, float> > map_indices_and_scales_vec =
         pyobject_to_map_index_and_scale_vec(map_number_and_scales);

      std::vector<std::pair<clipper::Xmap<float> *, float> > maps_and_scales_vec;

      graphics_info_t g;
      for (unsigned int i=0; i<map_indices_and_scales_vec.size(); i++) {
         int idx = map_indices_and_scales_vec[i].first;
         float w = map_indices_and_scales_vec[i].second;
         if (is_valid_map_molecule(idx)) {
            maps_and_scales_vec.push_back(std::make_pair(&g.molecules[idx].xmap, w));
         }
      }
      if (! maps_and_scales_vec.empty()) {
         coot::util::regen_weighted_map(&g.molecules[imol_map].xmap, maps_and_scales_vec);
         status = true;
      }
   }
}
#endif

/*! \brief
We overwrite the imol_map and we also presume that the
grid sampling of the contributing maps match. This makes it
much faster to generate than an average map.
*/
void regen_map_internal(int imol_map, const std::vector<std::pair<int, float> > &weighted_map_indices) {

   if (!weighted_map_indices.empty()) {
      graphics_info_t g;
      std::vector<std::pair<clipper::Xmap<float> *, float> > maps_and_scales_vec;
      for (unsigned int i=0; i<weighted_map_indices.size(); i++) {
         int idx = weighted_map_indices[i].first;
         float w = weighted_map_indices[i].second;
         std::pair<clipper::Xmap<float> *, float> p(&g.molecules[idx].xmap, w);
         maps_and_scales_vec.push_back(p);
      }
      coot::util::regen_weighted_map(&g.molecules[imol_map].xmap, maps_and_scales_vec);
   }
}

int make_weighted_map_simple_internal(const std::vector<std::pair<int, float> > &weighted_map_indices) {

   int imol = -1;

   if (!weighted_map_indices.empty()) {
      int imol_first = weighted_map_indices[0].first;
      imol = copy_molecule(imol_first);
      regen_map_internal(imol, weighted_map_indices);
   }
   return imol;
}

#ifdef USE_PYTHON
PyObject *positron_pathway(PyObject *map_molecule_list_py, PyObject *pathway_points_py) {

   // e.g.
   //  2,-4
   //  1,-3
   //  0,-2
   // -1,-1.4
   // -2,-0.3
   // -3, 0.8
   // -4, 2
   // -5, 3

   auto make_map = [] (const coot::positron_metadata_t &md,
                       const std::vector<int> &map_index_vec) {

      int imol = -1;
      if (! md.params.empty()) {
         if (md.params.size() == map_index_vec.size()) {
            PyObject *o = PyList_New(md.params.size());
            for (unsigned int i=0; i<md.params.size(); i++) {
               PyObject *item_py = PyList_New(2);
               PyList_SetItem(item_py, 0, PyLong_FromLong(map_index_vec[i]));
               PyList_SetItem(item_py, 1, PyFloat_FromDouble(md.params[i]));
               PyList_SetItem(o, i, item_py);
            }
            // imol = average_map_py(o);  // cubic interpolation
            int imol_first = map_index_vec[0];
            imol = copy_molecule(imol_first);
            // Use the new regen_map_internal()
            regen_map_py(imol, o);
         }
      }
      return imol;
   };

   float default_contour_level = 0.02;

   std::vector<int> new_map_index_list;
   if (PyList_Check(map_molecule_list_py)) {
      if (PyList_Check(pathway_points_py)) {
         int lmml = PyObject_Length(map_molecule_list_py);
         std::vector<int> map_index_list;
         for (int i=0; i<lmml; i++) {
            int ii = PyLong_AsLong(PyList_GetItem(map_molecule_list_py, i));
            map_index_list.push_back(ii);
         }
         if (map_index_list.size() == 6) { // the size of the params
            int l = PyObject_Length(pathway_points_py);
            for (int i=0; i<l; i++) {
               PyObject *x_y_point_py = PyList_GetItem(pathway_points_py, i);
               PyObject *x_py = PyList_GetItem(x_y_point_py, 0);
               PyObject *y_py = PyList_GetItem(x_y_point_py, 1);
               double x = PyFloat_AsDouble(x_py);
               double y = PyFloat_AsDouble(y_py);
               int idx_close = coot::get_closest_positron_metadata_point(graphics_info_t::positron_metadata, x, y);
               std::cout << "----------- i " << i << " idx_close " << idx_close << std::endl;
               if (idx_close != -1) {
                  coot::positron_metadata_t pmdi = graphics_info_t::positron_metadata[idx_close];
                  int imol_map_new = make_map(pmdi, map_index_list);
                  if (imol_map_new != -1) {
                     set_contour_level_absolute(imol_map_new, default_contour_level);
                     new_map_index_list.push_back(imol_map_new);
                  }
               }
            }
         }
      }
   }
   // convert new_map_index_list to a python list
   PyObject *new_map_index_list_py = PyList_New(new_map_index_list.size());
   for (unsigned int i=0; i<new_map_index_list.size(); i++) {
      PyObject *o = PyLong_FromLong(new_map_index_list[i]);
      PyList_SetItem(new_map_index_list_py, i, o);
   }
   return new_map_index_list_py;
}
#endif



/*  ----------------------------------------------------------------------- */
/*                  variance map                                            */
/*  ----------------------------------------------------------------------- */
//! \name Variance Map
//! \{
//! \brief Make a variance map, based on the grid of the first map.
//!
int make_variance_map(const std::vector<int> &map_molecule_number_vec) {

   int imol_map = -1;

   bool is_em_flag = false;
   std::vector<std::pair<clipper::Xmap<float>, float> > xmaps;
   for (unsigned int i=0; i<map_molecule_number_vec.size(); i++) {
      int imol = map_molecule_number_vec[i];
      is_em_flag = graphics_info_t::molecules[imol].is_EM_map();
      if (is_valid_map_molecule(imol)) {
         float scale = 1.0;
         xmaps.push_back(std::pair<clipper::Xmap<float>, float> (graphics_info_t::molecules[imol].xmap, scale));
      }
   }
   std::cout << "debug:: map_molecule_number_vec size " << map_molecule_number_vec.size() << std::endl;
   std::cout << "debug:: xmaps size " << xmaps.size() << std::endl;
   if (xmaps.size()) {
      clipper::Xmap<float> variance_map = coot::util::variance_map(xmaps);
      int imol = graphics_info_t::create_molecule();
      std::string name = "variance-map";
      graphics_info_t::molecules[imol].install_new_map(variance_map, name, is_em_flag);
      imol_map = imol;
      graphics_draw();
   }
   return imol_map;
}
//! \}


#ifdef USE_GUILE
int make_variance_map_scm(SCM map_molecule_number_list) {

   std::vector<int> v;
   SCM n_scm = scm_length(map_molecule_number_list);
   int n = scm_to_int(n_scm);
   for (int i=0; i<n; i++) {
      SCM mol_number_scm = scm_list_ref(map_molecule_number_list, scm_from_int(i));
      if (scm_is_true(scm_integer_p(mol_number_scm))) {
         int map_number = scm_to_int(mol_number_scm);
         if (is_valid_map_molecule(map_number))
            v.push_back(map_number);
      }
   }
   return make_variance_map(v);
}
#endif // USE_GUILE

#ifdef USE_PYTHON
int make_variance_map_py(PyObject *map_molecule_number_list) {

   std::vector<int> v;
   if (PyList_Check(map_molecule_number_list)) {
      int n = PyObject_Length(map_molecule_number_list);
      for (int i=0; i<n; i++) {
         PyObject *mol_number_py = PyList_GetItem(map_molecule_number_list, i);
         if (PyLong_Check(mol_number_py)) {
            int map_number = PyLong_AsLong(mol_number_py);
            if (is_valid_map_molecule(map_number)) {
               v.push_back(map_number);
            }
         }
      }
   }
   return make_variance_map(v);
}
#endif // USE_PYTHON


/* ------------------------------------------------------------------------- */
/*                      correllation maps                                    */
/* ------------------------------------------------------------------------- */

// The atom radius is not passed as a parameter to correlation
// functions, let's set it here (default is 1.5A)
void set_map_correlation_atom_radius(float r) {
   graphics_info_t::map_to_model_correlation_atom_radius = r;
}

// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-exclusing CB if is standard amino-acid, else all atoms
//
float
map_to_model_correlation(int imol,
    const std::vector<coot::residue_spec_t> &specs,
    const std::vector<coot::residue_spec_t> &neighb_specs,
    unsigned short int atom_mask_mode,
    int imol_map) {

   coot::util::density_correlation_stats_info_t dcs =
      map_to_model_correlation_stats(imol, specs, neighb_specs, atom_mask_mode, imol_map);
   return dcs.correlation();
}

// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-exclusing CB if is standard amino-acid, else all atoms
//
coot::util::density_correlation_stats_info_t
map_to_model_correlation_stats(int imol,
                               const std::vector<coot::residue_spec_t> &specs,
                               const std::vector<coot::residue_spec_t> &neighb_specs,
                               unsigned short int atom_mask_mode,
                               int imol_map) {

   coot::util::density_correlation_stats_info_t dcs;
   float atom_radius = graphics_info_t::map_to_model_correlation_atom_radius;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         clipper::Xmap<float> xmap_reference = graphics_info_t::molecules[imol_map].xmap;
         coot::map_stats_t map_stat_flag = coot::WITH_KOLMOGOROV_SMIRNOV_DIFFERENCE_MAP_TEST;
         dcs = coot::util::map_to_model_correlation_stats(mol, specs, neighb_specs,
                                                          atom_mask_mode,
                                                          atom_radius, xmap_reference,
                                                          map_stat_flag);
      }
   }
   return dcs;
}


#ifdef USE_GUILE
SCM map_to_model_correlation_scm(int imol,
    SCM residue_specs,
    SCM neighb_residue_specs,
    unsigned short int atom_mask_mode,
    int imol_map) {

   std::vector<coot::residue_spec_t> residues    = scm_to_residue_specs(residue_specs);
   std::vector<coot::residue_spec_t> nb_residues = scm_to_residue_specs(neighb_residue_specs);
   float c = map_to_model_correlation(imol, residues, nb_residues, atom_mask_mode, imol_map);
   SCM ret_val = scm_from_double(c);
   return ret_val;
}
#endif

#ifdef USE_GUILE
SCM map_to_model_correlation_stats_scm(int imol,
          SCM residue_specs,
          SCM neighb_residue_specs,
          unsigned short int atom_mask_mode,
          int imol_map) {

   std::vector<coot::residue_spec_t> residues    = scm_to_residue_specs(residue_specs);
   std::vector<coot::residue_spec_t> nb_residues = scm_to_residue_specs(neighb_residue_specs);
   coot::util::density_correlation_stats_info_t dcs
      = map_to_model_correlation_stats(imol, residues, nb_residues, atom_mask_mode, imol_map);

   double map_mean = -999.0;
   double map_sd   =    1.0;
   if (is_valid_map_molecule(imol_map)) {
      map_mean = graphics_info_t::molecules[imol_map].map_mean();
      map_sd   = graphics_info_t::molecules[imol_map].map_sigma();
   }
   double D = coot::stats::get_kolmogorov_smirnov_vs_normal(dcs.density_values,
       map_mean,
       map_sd);

   // Compare with the difference map at the ligand
   //
   coot::stats::single stats(dcs.density_values);
   double map_mean_at_ligand = stats.mean();
   double map_sd_at_ligand   = sqrt(stats.variance());
   double D2 = coot::stats::get_kolmogorov_smirnov_vs_normal(dcs.density_values,
        map_mean_at_ligand, map_sd_at_ligand);

   SCM ret_val = scm_list_n(scm_from_double(dcs.correlation()),
       scm_from_double(dcs.var_x()),
       scm_from_double(dcs.var_y()),
       scm_from_int(dcs.n),
       scm_from_double(dcs.sum_x),
       scm_from_double(dcs.sum_y),
       scm_from_double(D),
       scm_from_double(D2),
       scm_from_double(map_mean),
       scm_from_double(map_mean_at_ligand),
       scm_from_double(map_sd),
       scm_from_double(map_sd_at_ligand),
       SCM_UNDEFINED);

   return ret_val;
}
#endif


#ifdef USE_PYTHON
PyObject *map_to_model_correlation_py(int imol,
         PyObject *residue_specs,
         PyObject *neighb_residue_specs,
         unsigned short int atom_mask_mode,
         int imol_map) {

   std::vector<coot::residue_spec_t> residues    = py_to_residue_specs(residue_specs);
   std::vector<coot::residue_spec_t> nb_residues = py_to_residue_specs(neighb_residue_specs);
   float c = map_to_model_correlation(imol, residues, nb_residues,
         atom_mask_mode, imol_map);
   return PyFloat_FromDouble(c);
}
#endif

#ifdef USE_PYTHON
PyObject *map_to_model_correlation_stats_py(int imol,
       PyObject *residue_specs,
       PyObject *neighb_residue_specs,
       unsigned short int atom_mask_mode,
       int imol_map) {

   std::vector<coot::residue_spec_t> residues    = py_to_residue_specs(residue_specs);
   std::vector<coot::residue_spec_t> nb_residues = py_to_residue_specs(neighb_residue_specs);
   coot::util::density_correlation_stats_info_t dcs =
      map_to_model_correlation_stats(imol, residues, nb_residues, atom_mask_mode, imol_map);
   PyObject *r = PyList_New(6);
   PyList_SetItem(r, 0, PyFloat_FromDouble(dcs.correlation()));
   PyList_SetItem(r, 1, PyFloat_FromDouble(dcs.var_x()));
   PyList_SetItem(r, 2, PyFloat_FromDouble(dcs.var_y()));
   PyList_SetItem(r, 3, PyFloat_FromDouble(dcs.n));
   PyList_SetItem(r, 4, PyFloat_FromDouble(dcs.sum_x));
   PyList_SetItem(r, 5, PyFloat_FromDouble(dcs.sum_y));
   return r;
}
#endif

std::vector<std::pair<coot::residue_spec_t,float> >
map_to_model_correlation_per_residue(int imol, const std::vector<coot::residue_spec_t> &specs,
                                     unsigned short int atom_mask_mode,
                                     int imol_map) {

   float atom_radius = 1.5; // user variable?
   std::vector<std::pair<coot::residue_spec_t,float> > v;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap_reference = graphics_info_t::molecules[imol_map].xmap;
         v = coot::util::map_to_model_correlation_per_residue(mol, specs, atom_mask_mode, atom_radius, xmap_reference);
      }
   }
   return v;
}

//! \brief map to model density statistics, reported per residue
std::map<coot::residue_spec_t, coot::util::density_stats_info_t>
map_to_model_correlation_stats_per_residue(int imol,
      const std::vector<coot::residue_spec_t> &residue_specs,
      unsigned short int atom_mask_mode,
      float atom_radius,
      int imol_map) {

   std::map<coot::residue_spec_t, coot::util::density_stats_info_t> res_map;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap_reference = graphics_info_t::molecules[imol_map].xmap;
         res_map = coot::util::map_to_model_correlation_stats_per_residue(mol, residue_specs,
                                                                          atom_mask_mode, atom_radius,
                                                                          xmap_reference);
      }
   }
   return res_map;
}




#ifdef USE_GUILE
SCM
map_to_model_correlation_per_residue_scm(int imol, SCM specs_scm, unsigned short int atom_mask_mode, int imol_map) {

   SCM r = SCM_EOL;
   std::vector<coot::residue_spec_t> specs = scm_to_residue_specs(specs_scm);
   std::vector<std::pair<coot::residue_spec_t,float> >
      v = map_to_model_correlation_per_residue(imol, specs, atom_mask_mode, imol_map);
   for (unsigned int i=0; i<v.size(); i++) {
      SCM p1 = residue_spec_to_scm(v[i].first);
      SCM p2 = scm_from_double(v[i].second);
      SCM item = scm_list_2(p1, p2);
      r = scm_cons(item, r);
   }
   r = scm_reverse(r);
   return r;
}
#endif

#ifdef USE_GUILE
//! \brief map to model stats
SCM
map_to_model_correlation_stats_per_residue_scm(int imol,
          SCM specs_scm,
          unsigned short int atom_mask_mode,
          float atom_radius_for_masking,
          int imol_map) {
   SCM r = SCM_EOL;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         double map_mean = graphics_info_t::molecules[imol_map].map_mean();
         double map_sd   = graphics_info_t::molecules[imol_map].map_sigma();
         std::vector<coot::residue_spec_t> specs = scm_to_residue_specs(specs_scm);
         std::map<coot::residue_spec_t, coot::util::density_stats_info_t> res_map;
         res_map = map_to_model_correlation_stats_per_residue(imol, specs, atom_mask_mode, atom_radius_for_masking, imol_map);
         std::map<coot::residue_spec_t, coot::util::density_stats_info_t>::const_iterator it;
         for (it=res_map.begin(); it!=res_map.end(); it++) {
            const coot::residue_spec_t &spec = it->first;
            const coot::util::density_stats_info_t &dcs = it->second;
            SCM residue_spec_scm = residue_spec_to_scm(spec);

            double mean = dcs.sum/dcs.sum_weight;
            double var = mean * mean - dcs.sum_sq/dcs.sum_weight;
            std::pair<double, double> mv = dcs.mean_and_variance();
            mean = mv.first;
            var = mv.second;
            SCM mean_scm = scm_from_double(mean);
            SCM variance_scm = scm_from_double(var);
            SCM dcs_scm = scm_list_2(mean_scm, variance_scm);

            SCM item_scm = scm_list_2(residue_spec_scm, dcs_scm);
            r = scm_cons(item_scm, r);
         }
      }
   }
   r = scm_reverse(r);
   return r;

}
#endif


#ifdef USE_GUILE

SCM qq_plot_map_and_model_scm(int imol,
                              SCM residue_specs_scm,
                              SCM neighb_residue_specs_scm,
                              unsigned short int atom_mask_mode,
                              int imol_map) {

   SCM r = SCM_BOOL_F;

#ifdef HAVE_GOOCANVAS

   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         std::vector<coot::residue_spec_t> specs = scm_to_residue_specs(residue_specs_scm);
         std::vector<coot::residue_spec_t> nb_residues = scm_to_residue_specs(neighb_residue_specs_scm);
         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
         if (mol) {
            std::vector<std::pair<double, double> > v =
               coot::util::qq_plot_for_map_over_model(mol, specs, nb_residues, atom_mask_mode, xmap);

            if (0)
               for (unsigned int i=0; i<v.size(); i++)
                  std::cout << "     " << v[i].first << "    " << v[i].second << std::endl;

            // Goograph
            coot::goograph *g = new coot::goograph;
            int trace = g->trace_new();
            g->set_plot_title("Difference Map QQ Plot");
            g->set_axis_label(coot::goograph::X_AXIS, "Reference Normal Quantile");
            g->set_axis_label(coot::goograph::Y_AXIS, "Difference Map Quantile");
            g->set_trace_type(trace, coot::graph_trace_info_t::PLOT_TYPE_SCATTER);
            g->set_data(trace, v);
            std::pair<double, double> mmx = g->min_max_x();
            std::pair<double, double> mmy = g->min_max_y();
            lig_build::pos_t p1(mmx.first,  mmx.first);
            lig_build::pos_t p2(mmy.second, mmy.second);
            if (mmy.first  > mmx.first)  p1 = lig_build::pos_t(mmy.first,  mmy.first);
            if (mmy.second > mmx.second) p2 = lig_build::pos_t(mmy.second, mmy.second);
            g->add_annotation_line(p1, p2, "#444444", 1, false, false, false);
            g->show_dialog();
         }
      }
   }
#endif // HAVE_GOOCANVAS
   return r;
}
#endif


#ifdef USE_PYTHON
PyObject *
map_to_model_correlation_per_residue_py(int imol, PyObject *specs_py, unsigned short int atom_mask_mode, int imol_map) {

   std::vector<coot::residue_spec_t> specs = py_to_residue_specs(specs_py);
   std::vector<std::pair<coot::residue_spec_t,float> >
      v = map_to_model_correlation_per_residue(imol, specs, atom_mask_mode, imol_map);

   PyObject *r = PyList_New(v.size());
   for (unsigned int i=0; i<v.size(); i++) {
      PyObject *p0 = residue_spec_to_py(v[i].first);
      PyObject *p1 = PyFloat_FromDouble(v[i].second);
      PyObject *item = PyList_New(2);
      PyList_SetItem(item, 0, p0);
      PyList_SetItem(item, 1, p1);
      PyList_SetItem(r, i, item);
   }
   return r;
}

PyObject *qq_plot_map_and_model_py(int imol,
                                   PyObject *residue_specs_py,
                                   PyObject *neighb_residue_specs_py,
                                   unsigned short int atom_mask_mode,
                                   int imol_map) {

   PyObject *r = Py_False;
#ifdef HAVE_GOOCANVAS
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         std::vector<coot::residue_spec_t> specs = py_to_residue_specs(residue_specs_py);
         std::vector<coot::residue_spec_t> nb_residues = py_to_residue_specs(neighb_residue_specs_py);
         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
         if (mol) {
            std::vector<std::pair<double, double> > v =
               coot::util::qq_plot_for_map_over_model(mol, specs, nb_residues, atom_mask_mode, xmap);

            if (0)
               for (unsigned int i=0; i<v.size(); i++)
                  std::cout << "     " << v[i].first << "    " << v[i].second << std::endl;

            // Goograph
            coot::goograph *g = new coot::goograph;
            int trace = g->trace_new();
            g->set_plot_title("Difference Map QQ Plot");
            g->set_axis_label(coot::goograph::X_AXIS, "Reference Normal Quantile");
            g->set_axis_label(coot::goograph::Y_AXIS, "Difference Map Quantile");
            g->set_trace_type(trace, coot::graph_trace_info_t::PLOT_TYPE_SCATTER);
            g->set_data(trace, v);
            std::pair<double, double> mmx = g->min_max_x();
            std::pair<double, double> mmy = g->min_max_y();
            lig_build::pos_t p1(mmx.first,  mmx.first);
            lig_build::pos_t p2(mmy.second, mmy.second);
            if (mmy.first  > mmx.first)  p1 = lig_build::pos_t(mmy.first,  mmy.first);
            if (mmy.second > mmx.second) p2 = lig_build::pos_t(mmy.second, mmy.second);
            g->add_annotation_line(p1, p2, "#444444", 1, false, false, false);
            g->show_dialog();
         }
      }
   }

#endif // HAVE_GOOCANVAS
   // BL wonders:: always return false?,
   // PE:: Yes for the moment, nothing interesting is returned.
   if (PyBool_Check(r)) {
      Py_XINCREF(r);
   }
   return r;
}

#endif // USE_PYTHON

#ifdef USE_PYTHON
//! \brief return two lists: a list of vertices and a list of indices for connection
PyObject *map_contours(int imol, float contour_level) {

   PyObject *r = Py_False;

   if (is_valid_map_molecule(imol)) {
      graphics_info_t g;
      coot::Cartesian centre = g.RotationCentre();
      float radius = graphics_info_t::box_radius_xray;
      std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > contours =
         graphics_info_t::molecules[imol].get_contours(contour_level, radius, centre);

      std::cout << "got -------------------- " << contours.size() << " lines " << std::endl;
      r = PyList_New(contours.size());

      for (unsigned int i=0; i<contours.size(); i++) {
         PyObject *point_pair_py = PyList_New(2);
         PyObject *p1 = PyList_New(3);
         PyObject *p2 = PyList_New(3);
         PyList_SetItem(p1, 0, PyFloat_FromDouble(contours[i].first.x()));
         PyList_SetItem(p1, 1, PyFloat_FromDouble(contours[i].first.y()));
         PyList_SetItem(p1, 2, PyFloat_FromDouble(contours[i].first.z()));
         PyList_SetItem(p2, 0, PyFloat_FromDouble(contours[i].second.x()));
         PyList_SetItem(p2, 1, PyFloat_FromDouble(contours[i].second.y()));
         PyList_SetItem(p2, 2, PyFloat_FromDouble(contours[i].second.z()));

         PyList_SetItem(point_pair_py, 0, p1);
         PyList_SetItem(point_pair_py, 1, p2);

         PyList_SetItem(r, i, point_pair_py);
      }
   }
   if (PyBool_Check(r))
      Py_XINCREF(r);
   return r;
}
// \}
#endif // USE_PYTHON

//! \brief return two lists: a list of vertices and a list of indices for connection
PyObject *map_contours_as_triangles(int imol, float contour_level) {

   PyObject *r_py = Py_False;

   if (is_valid_map_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].update_map_internal();
      std::vector<glm::vec3> vertices = g.molecules[imol].map_as_mesh.just_vertices();
      const std::vector<g_triangle> &tris    =  g.molecules[imol].map_as_mesh.triangles;

      std::cout << "verticies size " << vertices.size() << std::endl;
      std::cout << "tris size " << tris.size() << std::endl;

      r_py = PyList_New(2);
      PyObject *vertices_py = PyList_New(vertices.size());
      PyObject *tris_py     = PyList_New(tris.size());
      //#if 0
      for (unsigned int i=0; i<vertices.size(); i++) {
         PyObject *vert_py = PyList_New(3);
         PyList_SetItem(vert_py, 0, PyFloat_FromDouble(vertices[i][0]));
         PyList_SetItem(vert_py, 1, PyFloat_FromDouble(vertices[i][1]));
         PyList_SetItem(vert_py, 2, PyFloat_FromDouble(vertices[i][2]));
         PyList_SetItem(vertices_py, i, vert_py);
      }
      for (unsigned int i=0; i<tris.size(); i++) {
         PyObject *tri_py = PyList_New(3);
         PyList_SetItem(tri_py, 0, PyLong_FromLong(tris[i][0]));
         PyList_SetItem(tri_py, 1, PyLong_FromLong(tris[i][1]));
         PyList_SetItem(tri_py, 2, PyLong_FromLong(tris[i][2]));
         PyList_SetItem(tris_py, i, tri_py);
      }
      // #endif
      PyList_SetItem(r_py, 0, vertices_py);
      PyList_SetItem(r_py, 1, tris_py);
   }

   if (PyBool_Check(r_py))
      Py_XINCREF(r_py);
   return r_py;

}


int sharpen_blur_map(int imol_map, float b_factor) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t g;
      imol_new = g.create_molecule();
      clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
      clipper::Xmap<float> xmap_new = coot::util::sharpen_blur_map(xmap, b_factor);
      std::string map_name = g.molecules[imol_map].name_; // use get_name() when it arrives
      if (b_factor < 0)
         map_name += " Sharpen ";
      else
         map_name += " Blur ";
      map_name += coot::util::float_to_string(b_factor);
      bool is_em_flag = graphics_info_t::molecules[imol_map].is_EM_map();
      g.molecules[imol_new].install_new_map(xmap_new, map_name, is_em_flag);
      float contour_level = graphics_info_t::molecules[imol_map].get_contour_level();
      graphics_info_t::molecules[imol_new].set_contour_level(contour_level);
      float cl = 5.0; // rmsd
      graphics_info_t::molecules[imol_new].set_contour_level_by_sigma(cl);
      graphics_draw();
   }
   return imol_new;
}

int sharpen_blur_map_with_resampling(int imol_map, float b_factor, float resample_factor) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t g;
      imol_new = g.create_molecule();
      clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
      clipper::Xmap<float> xmap_new = coot::util::sharpen_blur_map_with_resample(xmap, b_factor, resample_factor);
      std::string map_name = g.molecules[imol_map].name_; // use get_name() when it arrives
      if (b_factor < 0)
	 map_name += " Sharpen ";
      else
	 map_name += " Blur ";
      map_name += coot::util::float_to_string(b_factor);
      bool is_em_map_flag = g.molecules[imol_map].is_EM_map();
      g.molecules[imol_new].install_new_map(xmap_new, map_name, is_em_map_flag);
      float contour_level = g.molecules[imol_map].get_contour_level();
      g.molecules[imol_new].set_contour_level(contour_level);
      graphics_draw();
   }
   return imol_new;
}


void sharpen_blur_map_with_resampling_threaded_version(int imol_map, float b_factor, float resample_factor) {

   auto sharpen_blur_inner = +[] (std::promise<clipper::Xmap<float>> return_value,
                                  const clipper::Xmap<float> xmap,
                                  float b_factor,
                                  float resample_factor) {

      return_value.set_value(coot::util::sharpen_blur_map_with_resample(xmap, b_factor, resample_factor));
   };

   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t g;
      clipper::Xmap<float> xmap = g.molecules[imol_map].xmap;
      // make a name for the new map
      std::string map_name = g.molecules[imol_map].name_; // use get_name() when it arrives
      if (b_factor < 0)
         map_name += " Sharpen ";
      else
         map_name += " Blur ";
      map_name += coot::util::float_to_string(b_factor);
      bool is_em_map_flag = g.molecules[imol_map].is_EM_map();
      float contour_level = g.molecules[imol_map].get_contour_level();
      std::promise<clipper::Xmap<float>> computation_result_promise;

      struct sbr_callback_data_t {
         sbr_callback_data_t(const std::string &n, bool f, float cl, ProgressBarPopUp&& pp) : new_map_name(n), is_em_map_flag(f), contour_level(cl), popup(std::move(pp)) {}
         std::string new_map_name;
         bool is_em_map_flag;
         float contour_level;
         std::future<clipper::Xmap<float>> computation_result;
         ProgressBarPopUp popup;
      };
      sbr_callback_data_t *sbrcd_p = new sbr_callback_data_t(map_name, is_em_map_flag, contour_level, ProgressBarPopUp("Sharpen Blur", "Computing..."));
      sbrcd_p->computation_result = computation_result_promise.get_future();

      std::thread thread(sharpen_blur_inner, std::move(computation_result_promise), std::move(xmap), b_factor, resample_factor);
      thread.detach();

      auto check_it = +[] (gpointer data) {
         if(!data) {
            return FALSE;
         }
         sbr_callback_data_t *sbrcd_p = static_cast<sbr_callback_data_t *>(data);
         if (sbrcd_p->computation_result.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
            graphics_info_t g;
            int imol_new = g.create_molecule();
            auto result = sbrcd_p->computation_result.get();
            g.molecules[imol_new].install_new_map(result, sbrcd_p->new_map_name, sbrcd_p->is_em_map_flag);
            g.molecules[imol_new].set_contour_level(sbrcd_p->contour_level);
            g.set_imol_refinement_map(imol_new);
            graphics_draw();
            // hides the progress bar popup automatically
            delete sbrcd_p;
            return FALSE;
         } else {
            sbrcd_p->popup.pulse();
         }
         return TRUE;
      };

      GSourceFunc f = GSourceFunc(check_it);
      g_timeout_add(50, f, sbrcd_p);

   }

}


#ifdef USE_GUILE
//! \brief make many sharpened or blurred maps
//!
//! blurred maps are generated by using a positive value of b_factor.
//!
void multi_sharpen_blur_map_scm(int imol_map, SCM b_factors_list_scm) {

   if (is_valid_map_molecule(imol_map)) {

      std::vector<float> b_factors;
      SCM l_scm = scm_length(b_factors_list_scm);
      int l = scm_to_int(l_scm);
      for (int i=0; i<l; i++) {
         float f = scm_to_double(scm_list_ref(b_factors_list_scm, scm_from_int(i)));
         b_factors.push_back(f);
      }

      try {
	 graphics_info_t g;
	 const clipper::Xmap<float> &xmap_orig = g.molecules[imol_map].xmap;
	 std::vector<clipper::Xmap<float> > xmaps(b_factors.size());
	 coot::util::multi_sharpen_blur_map(xmap_orig, b_factors, &xmaps);
	 float contour_level = g.molecules[imol_map].get_contour_level();
	 for (std::size_t i=0; i<b_factors.size(); i++) {
	    const clipper::Xmap<float> &xmap_new = xmaps[i];
	    float b_factor = b_factors[i];
	    int imol_new = graphics_info_t::create_molecule();
	    std::string map_name = "Map";
	    if (b_factor < 0)
	       map_name += " Sharpen ";
	    else
	       map_name += " Blur ";
	    map_name += coot::util::float_to_string(b_factor);
	    bool is_em_map_flag = graphics_info_t::molecules[imol_map].is_EM_map();
	    g.molecules[imol_new].install_new_map(xmap_new, map_name, is_em_map_flag);
	    graphics_info_t::molecules[imol_new].set_contour_level(contour_level*exp(-0.02*b_factor));
	 }

      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }

}
#endif

#ifdef USE_PYTHON
void multi_sharpen_blur_map_py(int imol_map, PyObject *b_factors_list_py) {

   if (is_valid_map_molecule(imol_map)) {

      std::vector<float> b_factors;
      int l = PyObject_Length(b_factors_list_py);
      for (int i=0; i<l; i++) {
         PyObject *o = PyList_GetItem(b_factors_list_py, i);
         double b = PyFloat_AsDouble(o);
         b_factors.push_back(b);
      }

      try {
	 graphics_info_t g;
	 const clipper::Xmap<float> &xmap_orig = g.molecules[imol_map].xmap;
	 std::vector<clipper::Xmap<float> > xmaps(b_factors.size());
	 coot::util::multi_sharpen_blur_map(xmap_orig, b_factors, &xmaps);
	 float contour_level = g.molecules[imol_map].get_contour_level();
	 for (std::size_t i=0; i<b_factors.size(); i++) {
	    const clipper::Xmap<float> &xmap_new = xmaps[i];
	    float b_factor = b_factors[i];
	    int imol_new = graphics_info_t::create_molecule();
	    std::string map_name = "Map";
	    if (b_factor < 0)
	       map_name += " Sharpen ";
	    else
	       map_name += " Blur ";
	    map_name += coot::util::float_to_string(b_factor);
	    bool is_em_map_flag = graphics_info_t::molecules[imol_map].is_EM_map();
	    g.molecules[imol_new].install_new_map(xmap_new, map_name, is_em_map_flag);
	    graphics_info_t::molecules[imol_new].set_contour_level(contour_level*exp(-0.02*b_factor));
	 }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }

}
#endif


#ifdef USE_GUILE
SCM amplitude_vs_resolution_scm(int imol_map) {

   // return a list of (list sum count reso_average_recip)

   SCM r = SCM_EOL;
   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t g;
      const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
      // amplitude_vs_resolution decides the number of bins
      std::vector<coot::amplitude_vs_resolution_point> data = coot::util::amplitude_vs_resolution(xmap);
      std::cout << "amplitude_vs_resolution_scm() with data.size() " << data.size() << std::endl;
      for (std::size_t i=0; i<data.size(); i++) {
         SCM n = scm_list_3(scm_from_double(data[i].get_average_fsqrd()),
                            scm_from_int(data[i].count),
                            scm_from_double(data[i].get_invresolsq()));
         r = scm_cons(n, r);
      }

      std::pair<bool, float> l1(true,  0.05); // 4.5A
      std::pair<bool, float> l2(false, 0.29); // 1.8A
      float b = coot::util::b_factor(data, l1, l2);
      std::cout << "### b-factor: " << b << std::endl;
   }
   r = scm_reverse(r);
   return r;
}
#endif

float
b_factor_from_map(int imol_map) {

   float b_factor = -1;
   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t g;
      const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
      // amplitude_vs_resolution decides the number of bins
      std::vector<coot::amplitude_vs_resolution_point> data = coot::util::amplitude_vs_resolution(xmap);
      std::cout << "b_factor_from_map() with data.size() " << data.size() << std::endl;
      std::pair<bool, float> l1(true,  0.05); // 4.5A
      std::pair<bool, float> l2(false, 0.29); // 1.8A
      float b = coot::util::b_factor(data, l1, l2);
      std::cout << "### b-factor: " << b << std::endl;
   }
   return b_factor;
}

#ifdef USE_PYTHON
PyObject *amplitude_vs_resolution_py(int imol_map) {

   // return a list of [sum_fsqrd count reso_average_recip]

   PyObject *r = Py_False;

   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t g;
      clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
      std::vector<coot::amplitude_vs_resolution_point> data = coot::util::amplitude_vs_resolution(xmap);
      r = PyList_New(data.size());
      for (std::size_t i=0; i<data.size(); i++) {
         PyObject *o = PyList_New(3);
         // std::cout << "set o " << data[i].get_average_fsqrd() << std::endl;
         PyList_SetItem(o, 0, PyFloat_FromDouble(data[i].get_invresolsq()));
         PyList_SetItem(o, 1, PyLong_FromLong(data[i].count));
         PyList_SetItem(o, 2, PyFloat_FromDouble(data[i].get_average_fsqrd()));
         PyList_SetItem(r, i, o);
      }
   }

   if (PyBool_Check(r))
      Py_XINCREF(r);

   return r;
}
#endif


/*! \brief Calculate structure factors from the model and update the given difference
           map accordingly */
void sfcalc_genmap(int imol_model, int imol_map_with_data_attached, int imol_updating_difference_map) {

   graphics_info_t g;
   g.sfcalc_genmap(imol_model, imol_map_with_data_attached, imol_updating_difference_map);

}

/*! \brief As above, calculate structure factors from the model and update the given difference
           map accordingly - but difference map gets updated automatically on modification of
           the imol_model molecule */
void set_auto_updating_sfcalc_genmap(int imol_model,
                                     int imol_map_with_data_attached,
                                     int imol_updating_difference_map) {

   std::cout << "::::::::: set_auto_updating_sfcalc_genmap() --- start " << imol_model
             << " " << imol_map_with_data_attached << " " << imol_updating_difference_map
             << std::endl;

   // we need a notification that imol_model has been modified. Hmm.
   // Maybe the way to do that is that make_bonds type checked looks to see
   // if there is an "update_map" attached/related to this model - and then
   // update it. Hmm.

   // or maybe set up a timeout function that looks at the generation number of
   // imol_model (the backup number?) and if it is different to (more than) current,
   // then update the difference map - I like this plan more.

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_with_data_attached)) {
         if (is_valid_map_molecule(imol_updating_difference_map)) {
            if (map_is_difference_map(imol_updating_difference_map)) {

               if (false)
                  std::cout << "DEBUG:: in set_auto_updating_sfcalc_genmap() making a uump " << imol_model
                            << " " << imol_map_with_data_attached << " " << imol_updating_difference_map
                            << std::endl;

               updating_model_molecule_parameters_t ummp(imol_model, imol_map_with_data_attached, imol_updating_difference_map);
               updating_model_molecule_parameters_t *u = new updating_model_molecule_parameters_t(ummp);
               GSourceFunc f = GSourceFunc(graphics_info_t::molecules[imol_updating_difference_map].watch_coordinates_updates);
               graphics_info_t g;
               if (g.updating_maps_timeout_function_idx == UPDATING_MAPS_TIMEOUT_FUNCTION_IDX_UNSET)
                  g.updating_maps_timeout_function_idx = g_timeout_add(400, f, u);
               else
                  info_dialog("WARNING:: No can do.\nAn updating maps has already been started");
            }
         }
      }
   }
}

/*! \brief As above, calculate structure factors from the model and update the given difference
           map accordingly - but difference map gets updated automatically on modification of
           the imol_model molecule */
void set_auto_updating_sfcalc_genmaps(int imol_model, int imol_map_with_data_attached, int imol_updating_2fofc_map, int imol_updating_fofc_map) {

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_with_data_attached)) {
         if (is_valid_map_molecule(imol_updating_fofc_map)) {
            if (map_is_difference_map(imol_updating_fofc_map)) {
               if (is_valid_map_molecule(imol_updating_fofc_map)) {

                  updating_model_molecule_parameters_t ummp(imol_model, imol_map_with_data_attached,
                                                            imol_updating_2fofc_map, imol_updating_fofc_map);
                  updating_model_molecule_parameters_t *u = new updating_model_molecule_parameters_t(ummp);
                  // notice that the trigger in this case is on the *model* (not the difference map as above)
                  GSourceFunc f = GSourceFunc(graphics_info_t::molecules[imol_model].updating_coordinates_updates_genmaps);
                  g_timeout_add(700, f, u);
               }
            }
         }
      }
   }
}




//! \brief Go to the centre of the molecule - for Cryo-EM Molecules
//!
//!        and recontour at a sensible value.
void go_to_map_molecule_centre(int imol_map) {
   if (is_valid_map_molecule(imol_map)) {
      clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
      coot::util::map_molecule_centre_info_t mmci = coot::util::map_molecule_centre(xmap);
      if (mmci.success) {
         graphics_info_t::molecules[imol_map].set_contour_level(mmci.suggested_contour_level);
         clipper::Coord_orth nc = mmci.updated_centre;
         set_rotation_centre(nc.x(), nc.y(), nc.z());
      }
   }
}

//! \brief enable radial map colouring
void set_radial_map_colouring_enabled(int imol, int state) {

   if (is_valid_map_molecule(imol))
      graphics_info_t::molecules[imol].set_radial_map_colouring_do_radial_colouring(state);

   graphics_draw();
}


//! \brief radial map colouring centre
void set_radial_map_colouring_centre(int imol, float x, float y, float z) {
   if (is_valid_map_molecule(imol))
      graphics_info_t::molecules[imol].set_radial_map_colouring_centre(x,y,z);
}

//! \brief radial map colouring min
void set_radial_map_colouring_min_radius(int imol, float r) {
   if (is_valid_map_molecule(imol))
      graphics_info_t::molecules[imol].set_radial_map_colouring_min_radius(r);
}

//! \brief radial map colouring max
void set_radial_map_colouring_max_radius(int imol, float r) {
   if (is_valid_map_molecule(imol))
      graphics_info_t::molecules[imol].set_radial_map_colouring_max_radius(r);
}

//! \brief radial map colouring inverted colour map
void set_radial_map_colouring_invert(int imol, int invert_state) {
   if (is_valid_map_molecule(imol))
      graphics_info_t::molecules[imol].set_radial_map_colouring_invert(invert_state);
}

//! \brief radial map colouring saturation
//!
//! saturation is a number between 0 and 1, typically 0.5
void set_radial_map_colouring_saturation(int imol, float saturation) {
   if (is_valid_map_molecule(imol))
      graphics_info_t::molecules[imol].set_radial_map_colouring_saturation(saturation);
}

int flip_hand(int imol) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol)) {
      clipper::Xmap<float> xmap = graphics_info_t::molecules[imol].xmap;
      coot::util::flip_hand(&xmap);
      imol_new = graphics_info_t::create_molecule();
      std::string name = "Map ";
      name += coot::util::int_to_string(imol);
      name += " Flipped Hand";
      float contour_level = graphics_info_t::molecules[imol].get_contour_level();
      bool is_em_flag = graphics_info_t::molecules[imol].is_EM_map();
      graphics_info_t::molecules[imol_new].install_new_map(xmap, name, is_em_flag);
      graphics_info_t::molecules[imol_new].set_contour_level(contour_level);
      graphics_draw();
   }
   return imol_new;

}

void
add_density_map_cap() {

   int imol_map = imol_refinement_map();
   if (is_valid_map_molecule(imol_map)) {

      graphics_info_t g;
      clipper::Coord_orth base_point = g.get_rotation_centre_co();
      base_point -= clipper::Coord_orth(10, 10, 0);
      double x_axis_step_size = 0.5;
      double y_axis_step_size = 0.5;

      float z = -0.999; // screen z, front clipping plane
      glm::vec3 base        = graphics_info_t::unproject_to_world_coordinates(glm::vec3(-1.0f,-1.0f, z));
      glm::vec3 plus_x_axis = graphics_info_t::unproject_to_world_coordinates(glm::vec3(-1.0f, 1.0f, z));
      glm::vec3 plus_y_axis = graphics_info_t::unproject_to_world_coordinates(glm::vec3( 1.0f,-1.0f, z));

      clipper::Coord_orth base_co(base.x, base.y, base.z);
      clipper::Coord_orth plus_x_axis_co(plus_x_axis.x, plus_x_axis.y, plus_x_axis.z);
      clipper::Coord_orth plus_y_axis_co(plus_y_axis.x, plus_y_axis.y, plus_y_axis.z);
      clipper::Coord_orth delta_x_co = plus_x_axis_co - base_co;
      clipper::Coord_orth delta_y_co = plus_y_axis_co - base_co;

      double l = std::sqrt(delta_x_co.lengthsq());
      unsigned int n_x_axis_points = static_cast<int>(l/x_axis_step_size + 1);
      unsigned int n_y_axis_points = n_x_axis_points;

      std::cout << "debug:: base " << glm::to_string(base) << " x-axis " << glm::to_string(plus_x_axis)
                << std::endl;
      std::cout << "debug:: l " << l << " n_x_axis_points " << n_x_axis_points << std::endl;

      // clipper::Coord_orth x_axis_uv(1, 0, 0);
      // clipper::Coord_orth y_axis_uv(0, 1, 0);
      clipper::Coord_orth x_axis_uv(delta_x_co.unit());
      clipper::Coord_orth y_axis_uv(delta_y_co.unit());

      g.molecules[imol_map].setup_map_cap(&graphics_info_t::shader_for_map_caps,
                                          base_co, x_axis_uv, y_axis_uv,
                                          x_axis_step_size, y_axis_step_size,
                                          n_x_axis_points, n_y_axis_points);

      graphics_draw();

   }
}


//! \brief colour meshes (e.g. Ribbon diagrams) by map
//!
//! scale might be 2 and offset 1 (for example)
void recolour_mesh_by_map(int imol_model, int imol_map, float scale_factor, float offset) {

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         graphics_info_t g;
         const clipper::Xmap<float> &xmap(g.molecules[imol_map].xmap);
         g.molecules[imol_model].recolour_ribbon_by_map(xmap, scale_factor, offset);
         graphics_draw();
      }
   }
}

//! \brief test function for analysis of multiple map
int analyse_map_point_density_change(const std::vector<int> &map_number_list, int imol_map_mask) {

   std::vector<std::pair<clipper::Xmap<float> *, float> > xmaps;
   for (const auto &i : map_number_list) {
      if (graphics_info_t::is_valid_map_molecule(i)) {
         float rmsd = graphics_info_t::molecules[i].map_sigma();
         xmaps.push_back(std::make_pair(&graphics_info_t::molecules[i].xmap, rmsd));
      }
   }

   clipper::Xmap<float> xmap_for_mask;
   if (is_valid_map_molecule(imol_map_mask)) {
      xmap_for_mask = graphics_info_t::molecules[imol_map_mask].xmap;
   }

   std::cout << "DEBUG:: in analyse_map_point_density_change() with xmaps size " << xmaps.size() << std::endl;
   if (! xmaps.empty()) {

      // clipper::Xmap<float> linear_fit_map = coot::util::analyse_map_point_density_change(xmaps);
      clipper::Xmap<float> zde = coot::util::zero_dose_extrapolation(xmaps, xmap_for_mask);

      int new_molecule_number = graphics_info_t::create_molecule();
      bool is_EM_flag = true;
      std::string label = "negative linear_fit_of_decay";
      label = "zde";
      graphics_info_t::molecules[new_molecule_number].install_new_map(zde, label, is_EM_flag);
      // graphics_info_t::molecules[new_molecule_number].set_map_is_difference_map(true);
      return new_molecule_number;
   } else {
      return -1;
   }
}

#ifdef USE_PYTHON
int analyse_map_point_density_change_py(PyObject *map_number_list_py, int imol_map_mask) {

   std::vector<int> mnl;
   if (PyList_Check(map_number_list_py)) {
      int n = PyObject_Length(map_number_list_py);
      for (int i=0; i<n; i++) {
         PyObject *o = PyList_GetItem(map_number_list_py, i);
         if (PyLong_Check(o)) {  // this will need to be changed for Python3
            int imol = PyLong_AsLong(o);
            mnl.push_back(imol);
         }
      }
   }
   if (!mnl.empty()) {
      return analyse_map_point_density_change(mnl, imol_map_mask);
   } else {
      return -1;
   }
}
#endif

// use (or not) vertex gradients for the specified map
void set_use_vertex_gradients_for_map_normals(int imol, int state) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].set_use_vertex_gradients_for_map_normals(state);
      graphics_info_t::graphics_draw();
   }
}

//! the map should be displayed and not a difference map
void set_use_vertex_gradients_for_map_normals_for_latest_map() {
   use_vertex_gradients_for_map_normals_for_latest_map();
}

//! the map should be displayed and not a difference map
void use_vertex_gradients_for_map_normals_for_latest_map() {

   std::cout << "----------- use_vertex_gradients_for_map_normals_for_latest_map() ------ " << std::endl;

   int imol_map = -1;
   int n_mol = graphics_info_t::n_molecules();

   for (int i=(n_mol-1); i>=0; i--) {
      if (is_valid_map_molecule(i)) {
         if (graphics_info_t::molecules[i].is_displayed_p()) {
            if (graphics_info_t::molecules[i].is_difference_map_p())
               continue;
            imol_map = i;
            break;
         }
      }
   }

   if (imol_map != -1) {
      bool state = true;
      std::cout << "debug:: calling set_use_vertex_gradients_for_map_normals() for imol " << imol_map << std::endl;
      graphics_info_t::molecules[imol_map].set_use_vertex_gradients_for_map_normals(state);
      graphics_draw();
   }

}
