//
// Created by Jordan Dialpuri on 04/04/2026.
//

#include "ring.hh"
#include "primitives.hh"
#include <cmath>

std::array<glm::vec3, 6> ring_atoms_equal_bonds(
    float theta, float phi,
    float Q, float radius)
{
   const float pi = (float)M_PI;
   std::array<glm::vec3,6> coords;

   float q2 = Q * sin(theta);
   float q3 = Q * cos(theta);

   for (int j = 0; j < 6; j++) {
      float angle = 2*pi*j/6;

      float x = radius*cos(angle);
      float y = radius*sin(angle);

      float z = sqrt(2.0f/6.0f) *
          (q2*cos(2*pi*j/3 + phi) + q3*cos(pi*j));

      coords[j] = {x,y,z};
   }

   return coords;
}

glm::mat3 outward_frame(const glm::vec3 &dir) {
   glm::vec3 lz = glm::normalize(dir);
   glm::vec3 ref = (fabs(lz.z) < 0.9f) ? glm::vec3(0,0,1) : glm::vec3(0,1,0);
   glm::vec3 lx = glm::normalize(glm::cross(ref, lz));
   glm::vec3 ly = glm::cross(lz, lx);
   return {lx, ly, lz};
}

coot::simple_mesh_t make_ring(
    float cp_theta, float cp_phi,
    const glm::vec3 &center,
    const glm::vec3 &dir,
    float bond_length,
    const glm::vec4 &atom_colour,
    const glm::vec4 &bond_colour)
{
   coot::simple_mesh_t mesh;
   const glm::vec4  &oxygen_colour = {0.85f, 0.1f, 0.1f, 1.0f};

   auto frame = outward_frame(dir);
   auto local = ring_atoms_equal_bonds(cp_theta, cp_phi, 0.3f, 0.6f);

   std::array<glm::vec3,6> world;
   for (int i = 0; i < 6; i++) {
      world[i] = center + frame * local[i];
   }

   for (int i = 0; i < 6; i++) {
      mesh.add_submesh(make_stick(world[i], world[(i+1)%6],
                                 bond_length*0.08f, bond_colour));
   }

   for (int i = 0; i < 6; i++) {

      const glm::vec4 col = (i == 0) ? oxygen_colour : atom_colour;
      mesh.add_submesh(make_atom_sphere(world[i], bond_length * 0.12f, col));
   }

   return mesh;
}