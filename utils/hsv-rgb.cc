
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "coot-utils.hh"

std::vector<float>
coot::convert_rgb_to_hsv(const std::vector<float> &in_vals) {
   
   std::vector<float> cols(3);

   // convert to hsv
   float maxc = -1.0;
   float minc = 9.0;

   for (int i=0; i<3; i++) {
      if (maxc < in_vals[i]) maxc = in_vals[i];
      if (minc > in_vals[i]) minc = in_vals[i];
   }
   cols[2] = maxc;

   if (minc == maxc) {
      cols[0] = 0.0;
      cols[1] = 0.0;
      cols[2] = maxc;
   } else { 

      float range = maxc - minc; 
      cols[1] = range/maxc;
      float rc = (maxc - in_vals[0]) / range;
      float gc = (maxc - in_vals[1]) / range;
      float bc = (maxc - in_vals[2]) / range;
      if (in_vals[0] == maxc){ 
	 cols[0] = bc-gc;
      } else {
	 if (in_vals[1]==maxc) {
	    cols[0] = 2.0+rc-bc;
	 } else {
	    cols[0] = 4.0 + gc-rc;
	 }
      }
      cols[0] = cols[0]/6.0 - floorf(cols[0]/6.0);
   }
   return cols; 
}

std::vector<float>
coot::convert_hsv_to_rgb(const std::vector<float> &hsv)  {

   std::vector<float> rgb(3);

   if (hsv[1] == 0.0) {
      rgb[0] = hsv[2]; 
      rgb[1] = hsv[2]; 
      rgb[2] = hsv[2];
   } else {
      float fi = floorf(hsv[0]*6.0);
      float f  = (hsv[0]*6.0) - fi;
      float p = hsv[2]*(1.0 - hsv[1]);
      float q = hsv[2]*(1.0 - hsv[1]*f);
      float t = hsv[2]*(1.0 - hsv[1]*(1.0-f));

      int i = int(fi);
      switch (i) {

      case 0:
      case 6:
	 rgb[0] = hsv[2]; 
	 rgb[1] = t; 
	 rgb[2] = p;
	 break;

      case 1:
	 rgb[0] = q;
	 rgb[1] = hsv[2]; 
	 rgb[2] = p;
	 break;

      case 2:
	 rgb[0] = p;
	 rgb[1] = hsv[2]; 
	 rgb[2] = t;
	 break;

      case 3:
	 rgb[0] = p;
	 rgb[1] = q; 
	 rgb[2] = hsv[2];
	 break;

      case 4:
	 rgb[0] = t;
	 rgb[1] = p; 
	 rgb[2] = hsv[2];
	 break;

      case 5:
	 rgb[0] = hsv[2];
	 rgb[1] = p; 
	 rgb[2] = q;
	 break;
      }
   }
   return rgb; 
}


coot::colour_holder
coot::hsv_to_colour(const std::vector<float> &hsv) {

   std::vector<float> v = coot::convert_hsv_to_rgb(hsv);
   return colour_holder(v[0],v[1],v[2]);
}

// 
coot::colour_holder::colour_holder(const std::string &hex_colour_string) { 

   // fallback
   red = 0.5;
   green = 0.5;
   blue = 0.5;

   if (hex_colour_string.length() == 7) {
      if (hex_colour_string[0] == '#') {
	 std::string p_1 = hex_colour_string.substr(1,2);
	 std::string p_2 = hex_colour_string.substr(3,2);
	 std::string p_3 = hex_colour_string.substr(5,2);
	 int i_1, i_2, i_3;   
	 std::stringstream ss1;
	 std::stringstream ss2;
	 std::stringstream ss3;
	 ss1 << std::hex << p_1;
	 ss1 >> i_1;
	 ss2 << std::hex << p_2;
	 ss2 >> i_2;
	 ss3 << std::hex << p_3;
	 ss3 >> i_3;
	 red   = float(i_1)/255;
	 green = float(i_2)/255;
	 blue  = float(i_3)/255;
// debug	 
// 	 std::cout << "colour_holder hexstring " << hex_colour_string
// 		   << "  p_1  :" << p_1 << ": "
// 		   << "  p_2  :" << p_2 << ": "
// 		   << "  p_3  :" << p_3 << ": "
// 		   << " -> "
// 		   << i_1 << " " << i_2 << " " << i_3 << std::endl;
      } 
   } 
}

// // dum is a holder for a colour map selection.
// // 
coot::colour_holder::colour_holder(double value, double min_z, double max_z,
                                   bool use_deuteranomaly_mode,
				   const std::string &dum) {

   // Given a min, max range of 0,1
   // If value ~0, we want ~green
   // if value ~1, we want ~red

   float this_z = value;
   float range = max_z - min_z;
   float f = (this_z-min_z)/range;
   if (f > 1.0) f = 1.0;
   if (f < 0.0) f = 0.0;

   blue = 0.25 - (f-0.5)*(f-0.5);
   red = powf(f, 0.2);
   green = powf(1.0-f, 0.2);

   if (use_deuteranomaly_mode) {
      blue = f;
   }

}

void
coot::colour_holder::scale_intensity(float f) {

   red   *= f;
   green *= f;
   blue  *= f;

}

void
coot::colour_holder::rotate_by(float angle) {

   auto convert_to_hsv = [] (float red, float green, float blue) {
                            std::vector<float> v = { red, green, blue};
                            return convert_rgb_to_hsv(v);
                         };

   // references to member functions
   auto convert_from_hsv = [] (const std::vector<float> &v, float &red, float &green, float &blue) {
                              std::vector<float> o = convert_hsv_to_rgb(v);
                              red   = o[0];
                              green = o[1];
                              blue  = o[2];
                           };

   std::vector<float> hsv = convert_to_hsv(red, green, blue);
   hsv[0] += angle;
   while (hsv[0] > 1.0)
      hsv[0] -= 1.0;

   // not sure that this does any good. I need to test what convert_rgb_to_hsv() returns
   // for some sane and non-sane input values.
   while (hsv[0] < 0.0)
      hsv[0] += 1.0;

   convert_from_hsv(hsv, red, green, blue); // modify red, green, blue
}



std::ostream&
coot::operator<< (std::ostream& s, const coot::colour_holder &ch) {

   s << "colour{" << std::fixed << std::setprecision(3) << ch.red << " " << ch.green << " " << ch.blue
     << "}";
   return s;
} 

std::string
coot::colour_holder::hex() const {

   std::stringstream ss1;
   std::stringstream ss2;
   std::stringstream ss3;

   float c_red   = red;
   float c_green = green;
   float c_blue  = blue;

   if (c_red   < 0) c_red = 0; 
   if (c_red   > 1) c_red = 1; 
   if (c_green < 0) c_green = 0; 
   if (c_green > 1) c_green = 1; 
   if (c_blue  < 0) c_blue = 0; 
   if (c_blue  > 1) c_blue = 1;

   std::string hexstring = "#";

   ss1 << std::hex << std::setw(2) << std::setfill('0') << int(c_red  * 255);
   ss2 << std::hex << std::setw(2) << std::setfill('0') << int(c_green* 255);
   ss3 << std::hex << std::setw(2) << std::setfill('0') << int(c_blue * 255);
   hexstring += ss1.str();
   hexstring += ss2.str();
   hexstring += ss3.str();
   return hexstring;
} 

void
coot::colour_holder::brighten(float amount) {

   red   += amount;
   green += amount;
   blue  += amount;
   if (red   > 1.0) red   = 1.0;
   if (green > 1.0) green = 1.0;
   if (blue  > 1.0) blue  = 1.0;
}
