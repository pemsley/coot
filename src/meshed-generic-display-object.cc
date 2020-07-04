
#include <mmdb2/mmdb_manager.h>
#include "Python.h"

// Perhaps this (needed to compile meshed-generic-display-object.hh) should go in Mesh.hh.
#include <epoxy/gl.h>
#ifndef __APPLE__
#include <epoxy/glx.h>
#endif

#include "meshed-generic-display-object.hh"
#include "oct.hh"
#include "graphics-info.h"
#include "cylinder.hh"

glm::vec3 coord_orth_to_glm(const clipper::Coord_orth &co) {
   return glm::vec3(co.x(), co.y(), co.z());
}

void
meshed_generic_display_object::add_line(const coot::colour_holder &colour,
                                        const std::string &colour_name, int line_width,
                                        const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords) {
}


void
meshed_generic_display_object::add_point(const coot::colour_holder &colour_in,
                                         const std::string &colour_name,
                                         const int &size_in,
                                         const clipper::Coord_orth &coords_in) {

   unsigned int num_subdivisions = 3;
   float radius = 0.03 * size_in; // changing the scaling is fun
   glm::vec4 col(colour_in.red, colour_in.green, colour_in.blue, 1.0);
   glm::vec3 position = coord_orth_to_glm(coords_in);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
      oct = make_octasphere(num_subdivisions, position, radius, col);
   mesh.import(oct);

}

void
meshed_generic_display_object::add_cylinder(const std::pair<glm::vec3, glm::vec3 > &start_end,
                                            coot::colour_holder &col, float line_radius,
                                            unsigned int n_slices,
                                            bool cap_start, bool cap_end) {

   float h = glm::distance(start_end.first, start_end.second);
   cylinder c(start_end, line_radius, line_radius, h, n_slices, 2);
   glm::vec4 colour(col.red, col.green, col.blue, 1.0f);
   if (false)
      std::cout << "add_cylinder: " << glm::to_string(start_end.first) << " " << glm::to_string(start_end.second) << " "
                << c.vertices.size() << " " << c.triangle_indices_vec.size() << " with height " << h << std::endl;
   // c.add_flat_start_cap();
   // c.add_flat_end_cap();
   c.add_octahemisphere_start_cap();
   c.add_octahemisphere_end_cap();
   for (unsigned int i=0; i<c.vertices.size(); i++)
      c.vertices[i].color = colour;
   mesh.import(c.vertices, c.triangle_indices_vec);

}


void
meshed_generic_display_object::add_dodecahedron(const coot::colour_holder &colour_in,
                                                const std::string &colour_name,
                                                double radius, const clipper::Coord_orth &pos) {
}

void
meshed_generic_display_object::add_pentakis_dodecahedron(const coot::colour_holder &colour_in,
                                                         const std::string &colour_name,
                                                         double stellation_factor,
                                                         double radius,
                                                         const clipper::Coord_orth &pos) {

}
void
meshed_generic_display_object::add_arc(const arc_t &arc) {

}


coot::colour_holder
colour_values_from_colour_name(const std::string &c) {

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
            if (c == "greentint") {
               colour.red = 0.45;
               colour.green = 0.63;
               colour.blue = 0.45;
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

//    std::cout << "debug:: in colour_values_from_colour_name from colour " << c
// 	     << " we assign colour values "
// 	     << colour[0] << " "
// 	     << colour[1] << " "
// 	     << colour[2] << "\n";
   return colour;
}

