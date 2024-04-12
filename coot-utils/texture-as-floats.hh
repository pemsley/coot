#ifndef COOT_TEXTURE_AS_FLOATS
#define COOT_TEXTURE_AS_FLOATS

#include <clipper/core/xmap.h>

class texture_as_floats_t {
   public:
   texture_as_floats_t() : width(0), height(0), x_size(0), y_size(0), z_position(0) {}
   texture_as_floats_t(const clipper::Xmap<float> &xmap, int section_index);
   int width;
   int height;
   float x_size;
   float y_size;
   float z_position;
   std::vector<float> image_data;
   bool empty() { return image_data.empty(); }
};

#endif // COOT_TEXTURE_AS_FLOATS
