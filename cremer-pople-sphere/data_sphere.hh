//
// Created by Jordan Dialpuri on 04/04/2026.
//

#ifndef COOT_DATA_SPHERE_H
#define COOT_DATA_SPHERE_H

#include <glm/glm.hpp>
#include "coot-utils/simple-mesh.hh"

float get_energy(float theta, float phi);
glm::vec4 energy_to_colour(float t);

coot::simple_mesh_t make_data_sphere(
    float radius,
    int lat_seg,
    int lon_seg
);

#endif  // COOT_DATA_SPHERE_H
