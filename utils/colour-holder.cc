
// 20230119-PE Move other colour_holder functions here not that this file exists



#include "colour-holder.hh"

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
// 	     << " we assign colour values "
// 	     << colour[0] << " "
// 	     << colour[1] << " "
// 	     << colour[2] << "\n";
   return colour;
}

