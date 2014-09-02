
#ifndef COOT_COLOUR_HH
#define COOT_COLOUR_HH

#include <GL/gl.h>

namespace coot { 
   class colour_t {
      void init(float r, float g, float b) {
	 col.resize(3);
	 col[0] = r;
	 col[1] = g;
	 col[2] = b;
      }
   public:
      std::vector<float> col;
      colour_t() { init(0.5, 0.5, 0.5); }
      colour_t(float r, float g, float b) { init(r,g,b); }
      void set(float r, float g, float b) { init(r,g,b); }
      void glcolor() const { glColor3f(col[0], col[1], col[2]); }
   };
   std::ostream& operator<<(std::ostream &s, colour_t col);
}
#endif
