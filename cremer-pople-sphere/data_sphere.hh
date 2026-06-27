//
// Created by Jordan Dialpuri on 04/04/2026.
//

#ifndef COOT_DATA_SPHERE_H
#define COOT_DATA_SPHERE_H

#include <vector>
#include <glm/glm.hpp>
#include "coot-utils/simple-mesh.hh"

float get_energy(float theta, float phi);
glm::vec4 energy_to_colour(float t);
glm::vec4 density_to_colour(float t);

// Energy-based coloring (no data required)
// coot::simple_mesh_t make_data_sphere(
//     float radius,
//     int lat_seg,
//     int lon_seg
// );

// Density-based coloring from a precomputed lat x lon grid (see density_map.hh)
coot::simple_mesh_t make_data_sphere(
    float radius,
    int lat_seg,
    int lon_seg,
    const std::vector<float>& density_grid
);

#endif  // COOT_DATA_SPHERE_H
