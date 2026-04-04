//
// Created by Jordan Dialpuri on 04/04/2026.
//


#include "primitives.hh"
#include "coot-utils/oct.hh"

#include <iostream>

coot::simple_mesh_t make_stick(
    const glm::vec3 &a,
    const glm::vec3 &b,
    float radius,
    const glm::vec4 &colour,
    int slices)
{
     cylinder cyl({a, b}, radius, radius,
                glm::length(b-a), colour, slices);
   return coot::simple_mesh_t(cyl);
}

coot::simple_mesh_t make_atom_sphere(
    const glm::vec3 &pos,
    float radius,
    const glm::vec4 &colour)
{
   auto oct = make_octasphere(1, pos, radius, colour);

   coot::simple_mesh_t m;
   m.vertices = oct.first;
   m.triangles = oct.second;
   return m;
}