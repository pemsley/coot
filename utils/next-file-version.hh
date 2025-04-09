
namespace coot {

  std::optional<int>
  parse_version(const std::string &filename,
		const std::string &stub, const std::string &extension,
		const std::string &directory);
  
  int next_file_version(const std::string& stub, const std::string& extension, const std::string& directory);

  // don't add the directory - let the calling function do that if it wants to
  std::string get_versioned_file_name(const std::string& stub, const std::string& extension, const std::string& directory);

}
