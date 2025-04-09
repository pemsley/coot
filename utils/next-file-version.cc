
#include <filesystem>
#include <iostream>
#include <string>

#include "next-file-version.hh"

namespace fs = std::filesystem;

std::optional<int>
coot::parse_version(const std::string &filename,
	      const std::string &stub, const std::string &extension,
	      const std::string &directory) {

   if (filename.length() < stub.length() + extension.length() + 1) {
      return std::nullopt; // Filename too short
   }

   if (filename == stub + "." + extension) {
      return 0;
   }

   if (filename.rfind("." + extension) == filename.length() - extension.length() - 1) {
      std::string base = filename.substr(0, filename.length() - extension.length() - 1);
      if (base.length() > stub.length() && base.substr(0, stub.length()) == stub && base[stub.length()] == '-') {
	 std::string version_str = base.substr(stub.length() + 1);
	 if (!version_str.empty() && std::all_of(version_str.begin(), version_str.end(), ::isdigit)) {
	    try {
	       return std::stoi(version_str);
	    } catch (const std::out_of_range& e) {
	       // Handle potential out-of-range for very large numbers
	       return std::nullopt;
	    } catch (const std::invalid_argument& e) {
	       // Should not happen due to ::isdigit check
	       return std::nullopt;
	    }
	 }
      }
   }

   return std::nullopt;
}

int
coot::next_file_version(const std::string& stub, const std::string& extension, const std::string& directory) {

   int max_version = -1;
   fs::path dir_path(directory);

   if (!fs::is_directory(dir_path)) {
      std::cout << "this path!" << std::endl;
      return 0; // Directory doesn't exist, so no existing files
   }

   for (const auto &entry : fs::directory_iterator(dir_path)) {
      if (entry.is_regular_file()) {
	 std::string filename = entry.path().filename().string();
	 if (auto version = parse_version(filename, stub, extension, directory)) {
	    max_version = std::max(max_version, *version);
	 }
      }
   }

  return max_version + 1;
}


std::string
coot::get_versioned_file_name(const std::string& stub, const std::string& extension, const std::string& directory) {

   int nfv = next_file_version(stub, extension, directory);
   std::string fn = stub + "." + extension;
   if (nfv > 0) {
      fn = stub + "-" + std::to_string(nfv) + "." + extension;
   }
   return fn;
}
