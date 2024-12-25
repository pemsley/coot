
#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/map.h>
#include <nanobind/ndarray.h>

#undef MAKE_ENHANCED_LIGAND_TOOLS
#include "molecules-container.hh"

#ifdef NB_VERSION_MAJOR
//! External refinement using servalcat, using data that has already been associated.
//!
//! @return the imol of the refined model.
int
molecules_container_t::servalcat_refine_xray_with_keywords(int imol, int imol_map, const std::string &output_prefix,
                                                           const nanobind::dict &key_value_pairs) {

   std::map<std::string, std::string> kvm;
   for (const auto& item : key_value_pairs) {
      try {
         nanobind::handle f = item.first;
         nanobind::handle s = item.second;
         std::string key   = nanobind::cast<std::string>(f);
         std::string value = nanobind::cast<std::string>(s);
         kvm[key] = value;
      }
      catch (const nanobind::cast_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
   }

   return servalcat_refine_xray_internal(imol, imol_map, output_prefix, kvm);
}
#endif
