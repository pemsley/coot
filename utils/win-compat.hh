
#include <string>

namespace coot {

   std::string get_fixed_font();
   bool is_regular_file(const std::string &file_name);
   bool is_dir_or_link( const std::string &file_name);
   std::string uri_to_file_name(const std::string &file_name);
   int rename_win(const char *old_filename, const char *new_filename);

}
