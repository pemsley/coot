
#include <filesystem>
#include <fstream>
#include <iostream>

#include "utils/coot-utils.hh"
#include "inchikey-store.hh"

std::optional<std::pair<std::string, std::string> > coot::inchikey_store_t::get_match(const std::string &key) const {

   std::map<std::string, std::pair<std::string, std::string> >::const_iterator it = store.find(key);
   if (it == store.end()) {
      // failed to match
      return std::nullopt;
   } else {
      // return the match
      return it->second;
   }
}

void
coot::inchikey_store_t::parse_from_data_components_file() {

   bool debug = true;
   std::string pkg_data_dir = coot::package_data_dir();
   std::filesystem::path dir(pkg_data_dir);
   std::filesystem::path fn = dir / "data" / "Components-inchikey.ich";
   if (debug) {
      std::cout << "::: parse_from_data_components_file() dir: " << dir << std::endl;
      std::cout << "::: parse_from_data_components_file() pkg_data_dir: " << pkg_data_dir << std::endl;
      std::cout << "::: parse_from_data_components_file() " << fn << std::endl;
   }

   if (std::filesystem::exists(fn)) {
      std::ifstream f(fn);
      if (f) {
         std::string line;
         while(std::getline(f, line)) {
            std::vector<std::string> parts = util::split_string_no_blanks(line, "\t");
            if (parts.size() == 3) {
                const std::string &k = parts[0];
                const std::string &c = parts[1];
                const std::string &n = parts[2];
                store[k] = std::make_pair(c,n);
            } else {
               std::cout << "Failed to parse: " << parts.size() << " " << line << std::endl;
            }
         }
      }
   } else {
      std::cout << "WARNING:: File does not exist: " << fn << std::endl;
   }
   if (debug)
      std::cout << "DEBUG:: inchikey store size " << store.size() << std::endl;
}
