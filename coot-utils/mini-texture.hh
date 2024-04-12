#ifndef COOT_SRC_MINI_TEXTURE_HH
#define COOT_SRC_MINI_TEXTURE_HH

#include <clipper/core/xmap.h>

class mini_texture_t {
   public:
   mini_texture_t(int w, int h, unsigned char *d) : width(w), height(h), image_data(d) {}
   mini_texture_t(const clipper::Xmap<float> &xmap, int section_index);
   ~mini_texture_t();
   void clear();
   int width;
   int height;
   unsigned char *image_data;
};

#endif // COOT_SRC_MINI_TEXTURE_HH
