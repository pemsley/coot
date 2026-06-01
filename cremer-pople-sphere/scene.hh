//
// Created by Jordan Dialpuri on 04/04/2026.
//

#ifndef COOT_SCENE_H
#define COOT_SCENE_H

#include <string>
#include "coot-utils/simple-mesh.hh"
#include "data.h"


coot::simple_mesh_t make_mesh();

// Given Cremer-Pople angles (theta in [0,π], phi in [0,2π]), produce a ring
// mesh placed on the sphere surface and return it.  Suitable for calling when
// a viewer reports a click on the sphere at a known (theta, phi).
coot::simple_mesh_t make_ring_at_sphere_click(float theta, float phi);

// Place a pin (thin shaft + sphere head) on the sphere surface at (theta, phi).
coot::simple_mesh_t make_pinpoint_at_sphere_point(float theta, float phi);

// Create the Cremer-Pople sphere mesh based on the provided Cremer-Pople data
coot::simple_mesh_t create_cremer_pople_sphere(const std::vector<CremerPopleData>& cremer_pople_data);

#endif  // COOT_SCENE_H
