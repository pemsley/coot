
#include <vector>
#include "generic-vertex.hh"
#include "g_triangle.hh"


// construct these around the origin - the instanced matrices will put them
// in the right place.
std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
make_spherical_surface_circular_patch(float atom_radius,
                                      float solid_theta,
                                      float h_scale,
                                      float v_scale,
                                      unsigned int n_slices);


std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
make_curved_bar_patch(float atom_radius, float theta_start, float theta_end, float height,
                      unsigned int n_slices, float curve);
