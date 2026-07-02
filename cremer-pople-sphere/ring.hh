//
// Created by Jordan Dialpuri on 04/04/2026.
//

#ifndef COOT_RING_H
#define COOT_RING_H

#include <glm/glm.hpp>
#include <array>
#include "coot-utils/simple-mesh.hh"

std::array<glm::vec3, 6> ring_atoms_equal_bonds(
    float theta, float phi,
    float Q, float radius
);

glm::mat3 outward_frame(const glm::vec3 &dir);

coot::simple_mesh_t make_ring(
    float cp_theta, float cp_phi,
    const glm::vec3 &center,
    const glm::vec3 &dir,
    float bond_length,
    const glm::vec4 &atom_colour,
    const glm::vec4 &bond_colour
);

#endif  // COOT_RING_H
