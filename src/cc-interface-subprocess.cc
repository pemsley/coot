
#include "utils/subprocess.hpp"
#include "cc-interface.hh"
#include "graphics-info.h" // static state information in here

#include "utils/logging.hh"
extern logging logger;

std::pair<bool, std::string> graphics_info_t::acedrg_link = std::pair(false, "");

void
run_acedrg_link_generation(const std::string &acedrg_link_command) {

   auto run_acedrg_func = [] (const std::string &acedrg_link_command) {
      std::string link_command_file_name = "acedrg-link-in.txt";
      xdg_t xdg;
      std::ofstream ofs(link_command_file_name);
      ofs << acedrg_link_command << std::endl;
      ofs.close();
      std::cout << "DEBUG:: in run_acedrg_func: link-info acedrg-input file_name is " << link_command_file_name << std::endl;
      std::string cif_link_file_name_stub = "acedrg-link-from-coot";
      std::string acedrg_link_from_coot_file_name = cif_link_file_name_stub + "_link.cif";
      std::vector<std::string> cmd_list = {"acedrg", "-L", link_command_file_name, "-o", cif_link_file_name_stub};
      logger.log(log_t::INFO, "link CIF (acedrg ouput) file name will be", acedrg_link_from_coot_file_name);
      try {
         subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
         graphics_info_t::acedrg_link.second = acedrg_link_from_coot_file_name;
      }
      catch (const subprocess::CalledProcessError &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
      graphics_info_t::acedrg_link.first = true; // done
   };

   auto check_it = +[] (G_GNUC_UNUSED gpointer data) {

      // we are inside because this is a +[] lambda function
      auto file_name_to_string = +[] (const std::string &file_name) {
         std::string s;
         std::ifstream f(file_name.c_str(), std::ios::binary);
         if (!f) {
            std::cout << "WARNING:: Failed to open " << file_name << std::endl;
         } else {
            std::ostringstream ostrm;
            ostrm << f.rdbuf();
            s = ostrm.str();
         }
         return s;
      };

      if (graphics_info_t::acedrg_link.first) {
         std::string file_name = graphics_info_t::acedrg_link.second;
         if (! file_name.empty()) {
            // we can't do this (yet?) It's something about where the static lives.
            // graphics_info_t::log.log(logging::INFO, "read dictionary", file_name);
            // std::cout << "INFO:: read dictionary " << file_name << std::endl;
            logger.log(log_t::INFO, "read dictionary", file_name);
            add_status_bar_text("INFO:: read dictionary " + file_name);
            std::cout << "DEBUG:: about to read cif dictionary " << file_name << std::endl;
            read_cif_dictionary(file_name);
         } else {
            std::cout << "WARNING:: failed to make dictionary \"" << file_name << "\"" << std::endl;
            std::string err_info_log = "AcedrgOut_errorInfo.txt";
            if (coot::file_exists(err_info_log)) {
               std::string s = file_name_to_string(err_info_log);
               std::string ss = std::string("WARNING:: ") + s;
               graphics_info_t g;
               g.info_dialog(ss, false);
            } else {
               // std::cout << "INFO:: " << err_info_log << " does not exist" << std::endl;
               logger.log(log_t::INFO, err_info_log, "does not exist");
            }
            // graphics_info_t::log.log(logging::WARNING, "failed to make link dictionary", file_name);
         }
         graphics_info_t::acedrg_link.first = false; // reset
         GtkWidget *w = widget_from_builder("acedrg_running_frame");
         if (w) gtk_widget_set_visible(w, FALSE);
         return (gboolean)false; // remove the timeout
      } else {
         return (gboolean)true; // keep the timeout
      }
   };

   graphics_info_t::acedrg_link.first = false;
   std::thread thread(run_acedrg_func, acedrg_link_command);
   std::string err_info_log = "AcedrgOut_errorInfo.txt";
   if (coot::file_exists(err_info_log)) {
      std::filesystem::path p(err_info_log);
      std::filesystem::path new_p = p.parent_path() / "AcedrgOut_errorInfo.txt-previous";
      std::filesystem::rename(p, new_p);
   }
   thread.detach();
   GSourceFunc f = GSourceFunc(check_it);
   g_timeout_add(400, f, nullptr);
   GtkWidget *w = widget_from_builder("acedrg_running_frame");
   if (w) gtk_widget_set_visible(w, FALSE);

}
