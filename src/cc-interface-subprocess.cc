
#include "utils/subprocess.hpp"
#include "cc-interface.hh"
#include "graphics-info.h" // static state information in here

#include "utils/logging.hh"
extern logging logger;

std::pair<bool, std::string> graphics_info_t::acedrg_link = std::pair(false, "");
std::string graphics_info_t::acedrg_link_command_output = "";

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

      // Capture acedrg's own stdout and stderr (merged) ourselves, so that we always
      // have something concrete to show the user on failure, regardless of whether
      // acedrg happens to write its own AcedrgOut_errorInfo.txt for this kind of error.
      std::string acedrg_output_log_file_name = "acedrg-link-generation-output.log";
      graphics_info_t::acedrg_link_command_output.clear();
      try {
         subprocess::Popen p(cmd_list, subprocess::output{subprocess::PIPE}, subprocess::error{subprocess::STDOUT});
         auto res = p.communicate();
         std::string captured_output(res.first.buf.data(), res.first.length);
         std::ofstream log_ofs(acedrg_output_log_file_name);
         log_ofs << captured_output;
         log_ofs.close();
         if (p.retcode() == 0) {
            graphics_info_t::acedrg_link.second = acedrg_link_from_coot_file_name;
         } else {
            graphics_info_t::acedrg_link_command_output = captured_output;
            std::cout << "WARNING:: acedrg exited with retcode " << p.retcode() << std::endl;
         }
      }
      catch (const subprocess::CalledProcessError &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
         graphics_info_t::acedrg_link_command_output = e.what();
      }
      catch (const subprocess::OSError &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
         graphics_info_t::acedrg_link_command_output = e.what();
      }
      graphics_info_t::acedrg_link.first = true; // done
   };

   auto check_it = +[] (G_GNUC_UNUSED gpointer data) {

      // we are inside because this is a +[] lambda function
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
            std::string notification_text = "Link generated: " + file_name + "\nCurrent link applied to the model.";
            graphics_info_t::info_dialog(notification_text);
         } else {
            std::cout << "WARNING:: failed to make dictionary \"" << file_name << "\"" << std::endl;

            std::string ss = "WARNING:: Acedrg failed to generate the link.\n";

            const std::string &captured = graphics_info_t::acedrg_link_command_output;
            std::string acedrg_output_log = "acedrg-link-generation-output.log";

            if (captured.find("execve failed") != std::string::npos) {
               // acedrg could not even be started - almost always means it's not
               // installed, or not on PATH for this instance of Coot.
               ss += "\nAcedrg could not be started: " + captured +
                  "\n\nMake sure acedrg is installed and available on your PATH " +
                  "(e.g. source your CCP4 setup script) before starting Coot, then try again.";
            } else if (! coot::file_exists(acedrg_output_log) && ! captured.empty()) {
               // some other failure before/while spawning acedrg - no output log was
               // ever written, so show what we captured directly.
               ss += "\n" + captured;
            }

            // Our own capture of acedrg's stdout/stderr - this always exists when
            // acedrg was actually run (as opposed to failing to start at all).
            if (coot::file_exists(acedrg_output_log)) {
               ss += "\nCheck the log: " + acedrg_output_log;
            }

            // acedrg also sometimes writes its own error-specific log for chemistry problems.
            std::string err_info_log = "AcedrgOut_errorInfo.txt";
            if (coot::file_exists(err_info_log)) {
               std::string err_info_log_full_path = std::filesystem::absolute(err_info_log).string();
               ss += "\nSee also: " + err_info_log_full_path;
            }

            graphics_info_t g;
            g.info_dialog(ss, false);
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
   if (w) gtk_widget_set_visible(w, TRUE);

}
