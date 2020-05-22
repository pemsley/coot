
#ifndef COOT_CYLINDER_HH
#define COOT_CYLINDER_HH

#include "generic-vertex.hh"
#include "coords/Cartesian.h"

class tri_indices {
public:
   unsigned int idx[3];
   tri_indices() {}
   tri_indices(const unsigned int &i0, const unsigned int &i1, const unsigned int &i2) {
      idx[0] = i0;
      idx[1] = i1;
      idx[2] = i2;
   }
};

class cylinder {
public:
   std::vector<generic_vertex> vertices;
   std::vector<tri_indices> triangle_indices_vec;
   cylinder() {} // needs to be filled
   cylinder(const coot::CartesianPair &cp,
            float base_radius, float top_radius, float height,
            unsigned int n_slices, unsigned int n_stacks);
};

#endif
