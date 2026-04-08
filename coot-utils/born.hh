
#include <clipper/contrib/sfscale.h>

namespace coot {

   clipper::HKL_data<clipper::data32::F_phi> make_bohr(const clipper::HKL_data<clipper::data32::F_phi> &data_in);

   void born(const std::string &mtz_file_name, const std::string &f_col, const std::string &phi_col, const std::string &map_name);

}
