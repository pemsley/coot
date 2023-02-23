
#include "utils/colour-holder.hh"
#include "colour-functions.hh"

// amount is not in degrees, it is in fractions of a circle, e.g. 10/360.
//
std::vector<float> rotate_rgb(std::vector<float> &rgb, float amount) {

#if 0
   if (true) { // print a rotation to colour table
      int n_cols = 100; // either side
      for (int i = -n_cols; i<n_cols; i++) {
         float rotation_size = 0.01f * static_cast<float>(i);
         std::vector<float> orig_colours = { 0.0f,  0.8f, 0.0f };
         std::vector<float> rgb_new = rotate_rgb(orig_colours, rotation_size);
         std::cout << "debug colours::" << rgb_new[0] << " " << rgb_new[1] << " " << rgb_new[2]
                   << " using rotation_size " << rotation_size << std::endl;
      }

      // Result:
      //
      // if we start at solid green then rotation_size for "no rotation" is 0.0
      //                                 rotation_size for full rotation is -0.33 (solid red)
   }
#endif

   std::vector<float> hsv = coot::convert_rgb_to_hsv(rgb);

   // add 20 degrees to hue (or whatever)
   // std::cout << "adding " << amount << " to hue " << std::endl;
   hsv[0] += amount;
   while  (hsv[0] > 1.0) {
      hsv[0] -= 1.0;
   }

   std::vector<float> r = coot::convert_hsv_to_rgb(hsv);
   //     std::cout << "rotate from ("
   // << rgb[0] << " " << rgb[1] << " " << rgb[2] << ")\n"
   //           << "         to ("
   // << rgb[0] << " " << rgb[1] << " " << rgb[2] << ")\n";
   return r;
   // return convert_hsv_to_rgb(hsv);

}

