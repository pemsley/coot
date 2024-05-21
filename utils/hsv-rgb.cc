/*
 * utils/hsv-rgb.hh
 *
 * Copyright 2009 by University of Oxford
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


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

void
coot::colour_holder::pastelize(float degree) {

   float col[3] = {red, green, blue};
   for (unsigned int i=0; i<3; i++) {
      const float &cc = col[i];
      float r = 1.0f - cc;
      col[i] += r * degree;
      col[i] *= (1.0f - 0.5f * degree); // I don't want bright pastel
   }
   red   = col[0];
   green = col[1];
   blue  = col[2];

}

void
coot::colour_holder::make_pale(float degree) {

   float col[3] = {red, green, blue};
   for (unsigned int i=0; i<3; i++) {
      const float &cc = col[i];
      float delta_to_max = 1.0 - cc;
      float f = degree * delta_to_max;
      col[i] += f;
      if (col[i] < 0.0) col[i] = 0.0;
      if (col[i] > 1.0) col[i] = 1.0;
   }
   red   = col[0];
   green = col[1];
   blue  = col[2];
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

   s << "colour{" << std::fixed << std::setprecision(3) << ch.red << " " << ch.green << " " << ch.blue << " " << ch.alpha
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
