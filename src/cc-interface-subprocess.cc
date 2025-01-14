
#include "utils/subprocess.hpp"
#include "cc-interface.hh"
#include "graphics-info.h" // static state information in here

std::pair<bool, std::string> graphics_info_t::acedrg_link = std::pair(false, "");

void
run_acedrg_link_generation(const std::string &acedrg_link_command) {

   auto run_acedrg_func = [] (const std::string &acedrg_link_command) {
      std::string link_command_file_name = "acedrg-link-in.txt";
      xdg_t xdg;
      std::filesystem::path runtime_dir = xdg.get_runtime_dir();
      std::string file_name = (runtime_dir / link_command_file_name).string();
      std::ofstream ofs(file_name);
      ofs << acedrg_link_command << std::endl;
      ofs.close();
      std::vector<std::string> cmd_list = {"acedrg", "-L", file_name};
      subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
      graphics_info_t::acedrg_link.second = "dictionary.cif";
   };

   auto check_it = +[] (G_GNUC_UNUSED gpointer data) {
      if (graphics_info_t::acedrg_link.first) {
         std::string file_name = graphics_info_t::acedrg_link.second;
         if (! file_name.empty()) {
            // graphics_info_t::log.log(logging::INFO, "read dictionary", file_name);
            std::cout << "INFO:: read dictionary " << file_name << std::endl;
            read_cif_dictionary(file_name);
         } else {
            std::cout << "WARNING:: failed to make dictionary " << file_name << std::endl;
            // graphics_info_t::log.log(logging::WARNING, "failed to make link dictionary", file_name);
         }
      }
   };

   graphics_info_t::acedrg_link.first = false;
   std::thread thread(run_acedrg_func, acedrg_link_command);
   thread.detach();
   GSourceFunc f = GSourceFunc(check_it);
   g_timeout_add(400, f, nullptr);

}
