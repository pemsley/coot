#ifndef COOT_SRC_MINI_TEXTURE_HH
#define COOT_SRC_MINI_TEXTURE_HH

#include <clipper/core/xmap.h>

class mini_texture_t {
   public:
   mini_texture_t(int w, int h, unsigned char *d) :
      width(w), height(h), x_size(0), y_size(0), z_position(0), image_data(d),
      data_value_for_top_of_range(0), data_value_for_bottom_of_range(0) {}
   mini_texture_t(const clipper::Xmap<float> &xmap, int section_index,
                  float data_value_for_bottom, float data_value_for_top);
   ~mini_texture_t();
   void clear(); // delets image_data
   int width;
   int height;
   float x_size;
   float y_size;
   float z_position;
   unsigned char *image_data;
   float data_value_for_top_of_range;
   float data_value_for_bottom_of_range;
   bool empty() { return image_data == nullptr; }
};

#endif // COOT_SRC_MINI_TEXTURE_HH
