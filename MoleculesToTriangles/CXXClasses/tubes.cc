
#include <iostream>

#include "MoleculesToTriangles/CXXClasses/NRStuff.h"
#include "coot-utils/cylinder.hh"

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include "geometry/residue-and-atom-specs.hh"
#include "tubes.hh"

// Put this function somewhere outside of api - so that Coot can use it without getting tangled with api
//
coot::simple_mesh_t
make_tubes_representation(mmdb::Manager *mol,
                          const std::string &atom_selection_str,
                          const std::string &colour_scheme,
                          int secondaryStructureUsageFlag) {

   std::cout << "******************************** my make_tubes_representation() " << std::endl;
   coot::simple_mesh_t m_all;
   return m_all;
}
