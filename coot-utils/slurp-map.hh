

#include <string>

#include <clipper/ccp4/ccp4_map_io.h>
#include <sys/stat.h>

namespace coot {
   namespace util {

      bool is_basic_em_map_file(const std::string &file_name);

      // inf check_only is true, then just read the header, check that it is sane
      // and return that status (don't touch the xmap). Otherwise, fill the xmap.
      //
      bool slurp_fill_xmap_from_map_file(const std::string &file_name,
                                         clipper::Xmap<float> *xmap_p,
                                         bool check_only=false);

      bool slurp_parse_xmap_data(char *data, clipper::Xmap<float> *xmap_p, bool check_only=false);

   }

}
