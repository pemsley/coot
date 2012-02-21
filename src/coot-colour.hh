
#ifndef COOT_COLOUR_HH
#define COOT_COLOUR_HH

namespace coot { 
   class colour_t {
   public:
      std::vector<float> col;
      colour_t(float r, float g, float b) {
	 col.resize(3);
	 col[0] = r;
	 col[1] = g;
	 col[2] = b;
      }
      colour_t() {
	 col.resize(3);
	 col[0] = 0.5;
	 col[1] = 0.5;
	 col[2] = 0.5;
      }
   };
   std::ostream& operator<<(std::ostream &s, colour_t col);
}
#endif
