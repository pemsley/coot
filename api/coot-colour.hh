
#ifndef COOT_COLOUR_HH
#define COOT_COLOUR_HH

#include <vector>
#include <glm/glm.hpp>

#include "utils/colour-holder.hh"  // funny old thing. Consolidate the colour class one day.

namespace coot { 
   class colour_t {
      void init(float r, float g, float b) {
	 col.resize(3);
	 col[0] = r;
	 col[1] = g;
	 col[2] = b;
      }
      std::vector<float> convert_to_hsv() const;
      void convert_from_hsv(const std::vector<float> &hsv);
   public:
      std::vector<float> col;
      colour_t() { init(0.5, 0.5, 0.5); }
      colour_t(float r, float g, float b) { init(r,g,b); }
      void set(float r, float g, float b) { init(r,g,b); }
      float &operator[](const unsigned int &idx) { return col[idx]; }
      const float &operator[](const unsigned int &idx) const { return col[idx]; }
      void rotate(float f);
      void average(const colour_t &other) {
	 for (unsigned int idx=0; idx<3; idx++)
	    col[idx] = 0.5 * (col[idx] + other[idx]);
      }
      void brighter(float f) {
	 for (unsigned int idx=0; idx<3; idx++)
	    col[idx] *= f;
	 for (unsigned int idx=0; idx<3; idx++)
	    if (col[idx] > 1.0)
	       col[idx] = 1.0;
      }
      glm::vec4 to_glm() const {
         return glm::vec4(col[0], col[1], col[2], 1.0f);
      }
      colour_holder to_colour_holder() const {
         return colour_holder(col[0], col[1], col[2]);
      }
   };
   std::ostream& operator<<(std::ostream &s, colour_t col);
}
#endif
