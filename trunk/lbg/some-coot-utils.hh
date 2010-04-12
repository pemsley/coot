
#include <string>

namespace coot {

   bool is_directory_p(const std::string &dir);

   namespace util { 
      std::string int_to_string(int i);
      std::string float_to_string(float f);
   }

}
