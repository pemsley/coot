//
// Created by Jordan Dialpuri on 04/04/2026.
//


#include "primitives.hh"
#include "coot-utils/oct.hh"

#include <cmath>
#include <utility>
#include <iostream>

coot::simple_mesh_t make_stick(
    const glm::vec3 &a,
    const glm::vec3 &b,
    float radius,
    const glm::vec4 &colour,
    int slices)
{
   glm::vec3 start = a;
   glm::vec3 finish = b;

   // coot's cylinder mis-orients an axis pointing along +z (it treats it as the
   // pathological -z case and flips it). A cylinder is a symmetric tube, so for
   // an (almost) +z axis we build it from finish->start instead: that gives a
   // -z axis, which is oriented correctly, and the tube lands in the same place.
   const glm::vec3 n = glm::normalize(finish - start);
   if (std::fabs(n.x) < 1e-5f && std::fabs(n.y) < 1e-5f && n.z > 0.98f)
      std::swap(start, finish);

   cylinder cyl({start, finish}, radius, radius,
                glm::length(finish - start), colour, slices);
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