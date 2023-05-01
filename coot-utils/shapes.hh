#ifndef SHAPES_HH
#define SHAPES_HH

#include "shape-types.hh"
#include "simple-mesh.hh"

namespace coot {

   simple_mesh_t arrow_mesh(shapes::arrow_t a);

   simple_mesh_t cone_mesh(shapes::cone_t c);

   simple_mesh_t sphere_mesh(shapes::sphere_t s);

   simple_mesh_t torus_mesh(shapes::torus_t t);

   simple_mesh_t arc_mesh(shapes::arc_t a);

   // add dodec, pentakis-dodec at some stage

}


#endif // SHAPES_HH
