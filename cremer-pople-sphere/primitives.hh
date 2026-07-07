//
// Created by Jordan Dialpuri on 04/04/2026.
//

#ifndef COOT_PRIMITIVES_H
#define COOT_PRIMITIVES_H

#include <glm/glm.hpp>
#include "coot-utils/simple-mesh.hh"

coot::simple_mesh_t make_stick(
    const glm::vec3 &a,
    const glm::vec3 &b,
    float radius,
    const glm::vec4 &colour,
    int slices = 8
);

coot::simple_mesh_t make_atom_sphere(
    const glm::vec3 &pos,
    float radius,
    const glm::vec4 &colour
);

#endif  // COOT_PRIMITIVES_H
