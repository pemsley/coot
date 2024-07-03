
#include "subprocess.hpp"

// return the output mtz file name
std::string servalcat_fofc(const std::string &half_map_1, const std::string &half_map_2,
                           const std::string &pdb_file_name,
                           const std::string &prefix,
                           float resolution) {

   std::string output_fn = prefix + std::string(".mtz");
   std::vector<std::string> cmd_list = {"servalcat", "fofc",
      "--halfmaps", half_map_1, half_map_2,
      "--trim", "--trim_mtz", "--resolution", std::to_string(resolution),
      "--model", pdb_file_name, "-o", prefix};
   subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
   if (false) {
      std::cout << "Data : " << obuf.buf.data() << std::endl;
      std::cout << "Data len: " << obuf.length << std::endl;
   }
   return output_fn;
}

int main(int argc, char **argv) {

   int status = 0;

   std::thread thread(servalcat_fofc, "emd_32143_half_map_1.map", "emd_32143_half_map_2.map",
                      "pdb7vvl.ent", "test-7vvl-diff-map", 2.8);
   std::cout << "... waiting on join..." << std::endl;
   thread.join();
   std::cout << "... joined." << std::endl;
   return status;
}
