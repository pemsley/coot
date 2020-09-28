
#ifndef COLOUR_HOLDER_HH
#define COLOUR_HOLDER_HH

#include <vector>
#include <string>

namespace coot {
   class colour_holder {
   public:
      // values between 0 and 1.0
      float red;
      float green;
      float blue;
      colour_holder(const float &r, const float &g, const float &b) {
	 red = r;
	 green = g;
	 blue = b;
      }
      // needed because it's in a vector.
      colour_holder() {
	 red = 0.5;
	 green = 0.5;
	 blue = 0.5;
      }
      explicit colour_holder(const std::vector<float> &c_in) {
	 if (c_in.size() == 3) {
	    red   = c_in[0];
	    green = c_in[1];
	    blue  = c_in[2];
	 }
      }
      explicit colour_holder(const std::string &hex_colour_string);
      colour_holder(double value, double min, double max,
                    bool use_deuteranomaly_mode,
		    const std::string &dum); // somewhere between green and red
      std::string hex() const;
      void scale_intensity(float scale);
      void brighten(float brighten_amount); // presumes at full saturation, so no  simple scale of intensity
      friend std::ostream& operator<< (std::ostream& s, const colour_holder &ch);
   };
   std::ostream& operator<< (std::ostream& s, const colour_holder &ch);

   // colour conversion
   std::vector<float> convert_hsv_to_rgb(const std::vector<float> &hsv);
   std::vector<float> convert_rgb_to_hsv(const std::vector<float> &in_vals);
   colour_holder hsv_to_colour(const std::vector<float> &hsv);


}

#endif // COLOUR_HOLDER_HH
