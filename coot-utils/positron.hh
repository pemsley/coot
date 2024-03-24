#ifndef COOT_POSITRON_HH
#define COOT_POSITRON_HH

#include <string>
#include <vector>

namespace coot {
   class positron_metadata_t {
   public:
      float x;
      float y;
      std::vector<float> params;
      positron_metadata_t(float xx, float yy) : x(xx), y(yy) { params.resize(6,0); }
   };

   // should the following be part of a container class? If you like, yes.

   void read_positron_metadata(std::vector<positron_metadata_t> *data, const std::string &z_data, const std::string &s_data);

   // 20240323-PE Should we also pass a "must-be-closer-than-limit?"
   // Currently we use 0.1 in both x and y.
   // If we fail to find a close point, then return -1.
   int get_closest_positron_metadata_point(const std::vector<positron_metadata_t> &positron_metadata, float x, float y);

}

#endif // COOT_POSITRON_HH
