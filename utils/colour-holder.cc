/*
 * utils/colour-holder.cc
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


// 20230119-PE Move other colour_holder functions here not that this file exists

#include <cmath>
#include <sstream>
#include "colour-holder.hh"

//
coot::colour_holder::colour_holder(const std::string &hex_colour_string) {

   // fallback
   red = 0.5;
   green = 0.5;
   blue = 0.5;

   if (hex_colour_string.length() == 7 || hex_colour_string.length() == 9) {
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
         red   = static_cast<float>(i_1)/255.0f;
         green = static_cast<float>(i_2)/255.0f;
         blue  = static_cast<float>(i_3)/255.0f;
         alpha = 1.0f;

         if (hex_colour_string.length() == 9) {
            std::stringstream ss4;
            int i_a;
            std::string p_4 = hex_colour_string.substr(7,2);
            ss4 << std::hex << p_4;
            ss1 >> i_a;
            alpha = static_cast<float>(i_3)/255.0f;
         }
// debug
//          std::cout << "colour_holder hexstring " << hex_colour_string
//                    << "  p_1  :" << p_1 << ": "
//                    << "  p_2  :" << p_2 << ": "
//                    << "  p_3  :" << p_3 << ": "
//                    << " -> "
//                    << i_1 << " " << i_2 << " " << i_3 << std::endl;
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

   blue  = 0.25 - (f-0.5)*(f-0.5);
   red   = powf(f, 0.2);
   green = powf(1.0-f, 0.2);
   alpha = 1.0f;

   if (use_deuteranomaly_mode) {
      blue = f;
   }
}

coot::colour_holder
coot::colour_holder_from_colour_name(const std::string &c) {

   coot::colour_holder colour;
   colour.red = 0.4;
   colour.green = 0.4;
   colour.blue = 0.4;

   if (c.length() == 7) {
      if (c[0] == '#') {
         return coot::colour_holder(c); // hex colour string
      }
   }

   if (c.length() == 9) {
      if (c[0] == '#') {
         return coot::colour_holder(c); // hex colour string
      }
   }

   if (c == "blue") {
      colour.red = 0.1; colour.green = 0.1;
      colour.blue = 0.8;
   } else {
      if (c == "sky") {
         colour.red = 0.53 * 0.6;
         colour.green = 0.81 * 0.6;
         colour.blue = 0.92 * 0.6;
      } else {
         if (c == "green") {
            colour.red   = 0.05;
            colour.green = 0.8;
            colour.blue  = 0.05;
         } else {
            if (c == "greentint") {  // old Hydrogen-bond colour
               colour.red = 0.08;
               colour.green = 0.30;
               colour.blue = 0.08;
            } else {
               if (c == "darkpurple") { // new Hydrogen-bond colour
                  colour.red   = 0.48;
                  colour.green = 0.05;
                  colour.blue  = 0.5;
               } else {
                  if (c == "sea") {
                     colour.red = 0.1;
                     colour.green = 0.6;
                     colour.blue = 0.6;
                  } else {
                     if (c == "yellow") {
                        colour.red = 0.8;
                        colour.green = 0.8;
                        colour.blue = 0.0;
                     } else {
                        if (c == "yellowtint") {
                           colour.red = 0.65;
                           colour.green = 0.65;
                           colour.blue = 0.4;
                        } else {
                           if (c == "orange") {
                              colour.red = 0.9;
                              colour.green = 0.6;
                              colour.blue = 0.1;
                           } else {
                              if (c == "red") {
                                 colour.red = 0.9;
                                 colour.green = 0.1;
                                 colour.blue = 0.1;
                              } else {
                                 if (c == "hotpink") {
                                    colour.red = 0.9;
                                    colour.green = 0.2;
                                    colour.blue = 0.6;
                                 } else {
                                    if (c == "pink") {
                                       colour.red = 0.9;
                                       colour.green = 0.3;
                                       colour.blue = 0.3;
                                    } else {
                                       if (c == "cyan") {
                                          colour.red = 0.1;
                                          colour.green = 0.7;
                                          colour.blue = 0.7;
                                       } else {
                                          if (c == "aquamarine") {
                                             colour.red = 0.1;
                                             colour.green = 0.8;
                                             colour.blue = 0.6;
                                          } else {
                                             if (c == "forestgreen") {
                                                colour.red   = 0.6;
                                                colour.green = 0.8;
                                                colour.blue  = 0.1;
                                             } else {
                                                if (c == "yellowgreen") {
                                                   colour.red   = 0.6;
                                                   colour.green = 0.8;
                                                   colour.blue  = 0.2;
                                                } else {
                                                   if (c == "goldenrod") {
                                                      colour.red   = 0.85;
                                                      colour.green = 0.65;
                                                      colour.blue  = 0.12;
                                                   } else {
                                                      if (c == "orangered") {
                                                         colour.red   = 0.9;
                                                         colour.green = 0.27;
                                                         colour.blue  = 0.0;
                                                      } else {
                                                         if (c == "magenta") {
                                                            colour.red   = 0.7;
                                                            colour.green = 0.2;
                                                            colour.blue  = 0.7;
                                                         } else {
                                                            if (c == "cornflower") {
                                                               colour.red   = 0.38;
                                                               colour.green = 0.58;
                                                               colour.blue  = 0.93;
                                                            } else {
                                                               if (c == "royalblue") {
                                                                  colour.red   = 0.25;
                                                                  colour.green = 0.41;
                                                                  colour.blue  = 0.88;
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

//    std::cout << "debug:: in colour_values_from_colour_name from colour " << c
//              << " we assign colour values "
//              << colour[0] << " "
//              << colour[1] << " "
//              << colour[2] << "\n";
   return colour;
}

