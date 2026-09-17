
#include "utils/subprocess.hpp"
#include "coot-utils/json.hpp"
#include "graphics-info.h"
#include "read-molecule.hh"

int
graphics_info_t::servalcat_refine_xray_with_keywords(int imol, int imol_map, const std::string &output_prefix,
                                                     const std::string &keyword_pairs_json) {

   // ---------------- blocking! ----------------

   auto fill_key_map = [] (const std::string &keyword_pairs_json) {
      // keyword_pairs_json is a simple list of keyword key:data pair items
      //
      // We accept either a JSON object {"weight": "0.5", "ncycle": "10"}
      // or a JSON array of pairs [["weight", "0.5"], ["ncycle", "10"]].

      using json = nlohmann::json;

      // convert a json value to a string (values may arrive as numbers or bools)
      auto value_as_string = [] (const json &v) -> std::string {
         if (v.is_string()) return v.get<std::string>();
         return v.dump(); // numbers, bools, etc. -> their textual form
      };

      std::map<std::string, std::string> kvm;
      if (keyword_pairs_json.empty())
         return kvm;

      try {
         json j = json::parse(keyword_pairs_json);
         if (j.is_object()) {
            for (auto it = j.begin(); it != j.end(); ++it)
               kvm[it.key()] = value_as_string(it.value());
         } else if (j.is_array()) {
            for (const auto &item : j) {
               if (item.is_array() && item.size() == 2)
                  kvm[item[0].get<std::string>()] = value_as_string(item[1]);
            }
         } else {
            std::cout << "WARNING:: " << __FUNCTION__ << "(): unexpected JSON top-level type in "
                      << keyword_pairs_json << std::endl;
         }
      }
      catch (const json::exception &e) {
         std::cout << "WARNING:: " << __FUNCTION__ << "(): JSON parse error " << e.what()
                   << " for " << keyword_pairs_json << std::endl;
      }
      return kvm;
   };

   int imol_refined_model = -1;

   // 20260703-PE  this is copied from servalcat_refine_xray_internal() in molecules_container_t.
   // Maybe it should be refactored into a core/shared function ideal (say) - non-trivial.
   //
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {

         std::map<std::string, std::string> kvm;
         // fill kvm from keyword_pairs_json
         kvm = fill_key_map(keyword_pairs_json);

         bool set_weight = false;
         std::string weight_str;
         if (! kvm.empty()) {
            for (const auto &kv : kvm) {
               if (kv.first == "weight") {
                  set_weight = true;
                  weight_str = kv.second;
               }
            }
         }

         bool clibd_mon_is_set = false;
         char *e = getenv("CLIBD_MON");
         if (e) {
            std::string env(e);
            if (std::filesystem::exists(env))
               clibd_mon_is_set = true;
         }
         if (clibd_mon_is_set) {

            std::string mtz_file    = molecules[imol_map].Refmac_mtz_filename();
            std::string fobs_col    = molecules[imol_map].Refmac_fobs_col();
            std::string sigfobs_col = molecules[imol_map].Refmac_sigfobs_col();
            std::string r_free_col  = molecules[imol_map].Refmac_r_free_col();

            if (false) {
               std::cout << "debug:: mtz_file "    << mtz_file    << std::endl;
               std::cout << "debug:: fobs_col "    << fobs_col    << std::endl;
               std::cout << "debug:: sigfobs_col " << sigfobs_col << std::endl;
               std::cout << "debug:: r_free_col "  << r_free_col  << std::endl;
            }

            if (! mtz_file.empty()) {

               bool read_pdb_output = false; // this gets set to true if the output pdb is sane

               std::string c(",");
               std::string labin = fobs_col + c + sigfobs_col + c + r_free_col;

               std::string dir_1 = "coot-servalcat";
               coot::util::create_directory(dir_1);
               std::string prefix = coot::util::append_dir_file(dir_1, output_prefix);
               std::string  input_pdb_file_name = prefix + std::string("-in.pdb");
               std::string output_pdb_file_name = prefix + std::string(".pdb"); // named by servalcat
               int status = molecules[imol].export_coordinates(input_pdb_file_name);
               // see https://www.ebi.ac.uk/pdbe/docs/cldoc/object/cl_obj_rdwr.html#CMMDBManager::WritePDBASCII
               if (status == 0) {
                  bool output_pdb_file_name_exists = false;
                  std::filesystem::file_time_type output_pdb_file_name_time;
                  std::filesystem::path p(output_pdb_file_name);
                  if (std::filesystem::exists(p)) {
                     output_pdb_file_name_exists = true;
                     output_pdb_file_name_time = std::filesystem::last_write_time(p);
                  }
                  std::vector<std::string> cmd_list = {"servalcat", "refine_xtal_norefmac",
                                                       "-s", "xray", "--model", input_pdb_file_name,
                                                       "--hklin", mtz_file, "--labin", labin,
                                                       "-o", prefix};
                  if (set_weight) {
                     cmd_list.push_back("--weight");
                     cmd_list.push_back(weight_str);
                  }

                  if (true) {
                     std::cout << "commandline: ";
                     for (unsigned int i=0; i<cmd_list.size(); i++) std::cout << " " << cmd_list[i];
                     std::cout << "\n";
                  }
                  std::cout << "running servalcat..." << std::endl;
                  subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
                  if (std::filesystem::exists(p)) {
                     if (output_pdb_file_name_exists) {
                        std::filesystem::file_time_type new_output_pdb_file_name_time = std::filesystem::last_write_time(p);
                        auto t1 =     output_pdb_file_name_time.time_since_epoch();
                        auto t2 = new_output_pdb_file_name_time.time_since_epoch();
                        auto tt1 = std::chrono::duration_cast<std::chrono::seconds>(t1).count();
                        auto tt2 = std::chrono::duration_cast<std::chrono::seconds>(t2).count();
                        auto d = tt2 - tt1;
                        if (d > 0)
                           read_pdb_output = true;
                     } else {
                        read_pdb_output = true;
                     }
                     if (read_pdb_output) {
                        imol_refined_model = read_coordinates(output_pdb_file_name);
                     }
                  } else {
                     std::cout << "WARNING:: " << __FUNCTION__ << "(): path does not exist " << p << std::endl;
                  }
               } else {
                  std::cout << "WARNING::" << __FUNCTION__ << "(): bad status on writing servalcat input file" << std::endl;
               }
            } else {
               std::cout << "WARNING::" << __FUNCTION__ << "(): mtz file_name was empty" << std::endl;
            }
         } else {
            std::cout << "WARNING::" << __FUNCTION__ << "(): CLIBD_MON was not set correctly" << std::endl;
         }
      } else {
         std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

   return imol_refined_model;

}

