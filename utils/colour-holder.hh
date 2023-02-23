
#ifndef COLOUR_HOLDER_HH
#define COLOUR_HOLDER_HH

#include <vector>
#include <string>

namespace coot {

   //! a container for colours
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
      // this constructor is needed because colour_holder can be used in a vector.
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
      void pastelize(float degree);
      void make_pale(float degree);
      std::string hex() const;
      void rotate_by(float angle); // fractions of a circle
      void scale_intensity(float scale);
      void brighten(float brighten_amount); // presumes at full saturation, so no  simple scale of intensity
      friend std::ostream& operator<< (std::ostream& s, const colour_holder &ch);
   };
   std::ostream& operator<< (std::ostream& s, const colour_holder &ch);

   // colour conversion
   std::vector<float> convert_hsv_to_rgb(const std::vector<float> &hsv);
   std::vector<float> convert_rgb_to_hsv(const std::vector<float> &in_vals);
   colour_holder hsv_to_colour(const std::vector<float> &hsv);

   //! the constructor above uses a hash colour, the argument here is things like "red", "orange", "blue";
   //! if a hashed colour hex is give to this function, then the above constructor
   //! is used to generate the return value
   colour_holder colour_holder_from_colour_name(const std::string &colour_name);

}

#endif // COLOUR_HOLDER_HH
