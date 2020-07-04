
// Perhaps this should go in Mesh.hh.
#include <epoxy/gl.h>
#ifndef __APPLE__
#include <epoxy/glx.h>
#endif

#include "meshed-generic-display-object.hh"

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

