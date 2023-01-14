#ifndef VNC_VERTEX_TO_GENERIC_VERTEX_HH
#define VNC_VERTEX_TO_GENERIC_VERTEX_HH

#include "coot-utils/vertex.hh"
#include "generic-vertex.hh"

s_generic_vertex vnc_vertex_to_generic_vertex(const coot::api::vnc_vertex &v);
s_generic_vertex vn_vertex_to_generic_vertex(const coot::api::vn_vertex &v);

#endif // VNC_VERTEX_TO_GENERIC_VERTEX_HH
