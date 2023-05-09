
#ifdef USE_PYTHON
#include <Python.h>
#endif

#include "graphics-info.h"

/*! \brief shiftfield B-factor refinement */
void 
graphics_info_t::shiftfield_b_factor_refinement(int imol) {

   int imol_map = Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_map)) {
      try {
         molecules[imol_map].fill_fobs_sigfobs(); // caches
         const clipper::HKL_data<clipper::data32::F_sigF> *fobs_data = molecules[imol_map].get_original_fobs_sigfobs();
         const clipper::HKL_data<clipper::data32::Flag> *free_flag   = molecules[imol_map].get_original_rfree_flags();
         if (fobs_data && free_flag) {
            molecules[imol].shiftfield_b_factor_refinement(*fobs_data, *free_flag);
         } else {
            std::cout << "ERROR:: null pointer in function " << __FUNCTION__ << std::endl;
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }

}

/*! \brief shiftfield xyz refinement */
void
graphics_info_t::shiftfield_xyz_factor_refinement(int imol) {

   std::cout << "Not implemented." << std::endl;

}


//static
bool
graphics_info_t::showing_intermediate_atoms_from_refinement() {

   if (use_graphics_interface_flag) {
      if (moving_atoms_asc) {
         if (last_restraints) {
            return true;
         }
      }
   }
   return false;
}


int
graphics_info_t::rail_points_total() const { // the sum of all the rail points accumulated
   return api::rail_points_t::total(rail_point_history);
}

int
graphics_info_t::calculate_new_rail_points(const updating_model_molecule_parameters_t &ummp) {

   // std::cout << ":::::::::: calculate_new_rail_points() history size: " << rail_point_history.size() << std::endl;

   if (is_valid_map_molecule(ummp.imol_fofc_map)) {

      float rmsd = molecules[ummp.imol_fofc_map].get_map_sigma_current();
      if (rail_point_history.empty()) {
         api::rail_points_t prev = api::rail_points_t(rmsd);
         if (false)
            std::cout << "...... prev and rmsd " << prev.rmsd_of_difference_map << " " << rmsd
                      << " 2Cooooot-Points " << 100000 * (prev.rmsd_of_difference_map - rmsd)
                      << std::endl;
         api::rail_points_t new_points(rmsd, prev);
         rail_point_history.push_back(new_points);
         return new_points.map_rail_points_delta;
      } else {
         // std::cout << ":::::::::: calculate_new_rail_points() B "<< std::endl;
         const api::rail_points_t &prev = rail_point_history.back();
         api::rail_points_t new_points(rmsd, prev);
         rail_point_history.push_back(new_points);
         return new_points.map_rail_points_delta;
      }
   } else {
      return 0;
   }
}

void
graphics_info_t::updating_maps_update_the_coot_points_overlay() {

   GtkWidget *label_1 = get_widget_from_builder("coot-points-frame-points-label");
   GtkWidget *label_2 = get_widget_from_builder("coot-points-frame-r-factor-label"); // now total
   GtkWidget *label_3 = get_widget_from_builder("coot-points-frame-free-r-factor-label"); // R-factors

   if (rail_point_history.empty()) {
      // std::cout << "------------- update the overlay! A" << std::endl;
      gtk_label_set_text(GTK_LABEL(label_1), "-----");
      gtk_label_set_text(GTK_LABEL(label_2), "-----");
      gtk_label_set_text(GTK_LABEL(label_3), "-----");
   } else {
      int d = rail_point_history.back().map_rail_points_delta;
      // std::cout << "------------- update the overlay! B" << std::endl;
      std::string plus;
      if (d > 0) plus = "+";
      std::string colour = "#dddddd";
      if (d < 0) colour = "#ff3333";
      if (d > 0) colour = "#33ff33";
      std::string l_1 = std::string("<span foreground='");
      l_1 += colour;
      l_1 += std::string("'>");
      l_1 += "New Coot Points:   " + plus + std::to_string(d);
      l_1 += std::string("</span>");
      // std::cout << l_1 << std::endl;
      std::string l_2 = "Total Coot Points: " + std::to_string(api::rail_points_t::total(rail_point_history));
      std::string l_3 = "R-factors: ";
      l_3 += coot::util::float_to_string_using_dec_pl(100.0f * latest_sfcalc_stats.r_factor, 2);
      l_3 += "%, ";
      l_3 += coot::util::float_to_string_using_dec_pl(100.0f * latest_sfcalc_stats.free_r_factor, 2);
      l_3 += "%";
      gtk_label_set_markup(GTK_LABEL(label_1), l_1.c_str());
      gtk_label_set_text(GTK_LABEL(label_2), l_2.c_str());
      gtk_label_set_text(GTK_LABEL(label_3), l_3.c_str());
   }

   auto coot_points_frame_callback = +[] (gpointer user_data) {
      GtkWidget *frame = get_widget_from_builder("coot-points-frame");
      if (frame) {
         gtk_widget_set_visible(frame, FALSE);
      }
      return FALSE;
   };

   GtkWidget *frame = get_widget_from_builder("coot-points-frame");
   if (frame) {
      gtk_widget_set_visible(frame, TRUE);
      GSourceFunc cb = G_SOURCE_FUNC(coot_points_frame_callback);
      g_timeout_add(3000, cb, nullptr);
   }

}
