//
// Created by Jordan Dialpuri on 04/04/2026.
//


#include "scene.hh"
#include "data_sphere.hh"
#include "primitives.hh"
#include "ring.hh"
#include "coot-utils/shape-types.hh"
#include "coot-utils/shapes.hh"

coot::colour_holder to_colour(const glm::vec4 &c) {
   return coot::colour_holder(c.r, c.g, c.b, c.a);
}

coot::simple_mesh_t make_cage() {
   coot::simple_mesh_t mesh;
   mesh.set_name("Cremer-Pople Cage");

   constexpr float sphere_r = 4.0f;
   constexpr float tube_r    = 0.025f;

   constexpr int n_lat = 7;
   constexpr int n_lon = 6;
   constexpr int n_ring_atoms = 64;

   const coot::colour_holder color = {1.0f, 1.0f, 1.0f, 1.0f};

   for (int i = 0; i < n_lat; i++) {
      const float theta = M_PI * (static_cast<float>(i + 1) /(n_lat+1) - 0.5f);

      float z = sphere_r * sin(theta);
      const float R = sphere_r * cos(theta);

      shapes::torus_t torus(
          {0,0,z},
          {0,0,1},
          R,
          tube_r
      );

      torus.col = color;
      torus.n_ring_atoms = n_ring_atoms;

      mesh.add_submesh(coot::torus_mesh(torus));
   }

   for (int i = 0; i < n_lon; i++) {
      const float phi = M_PI * static_cast<float>(i) / n_lon;
      shapes::torus_t torus(
          clipper::Coord_orth(0,0,0),
          clipper::Coord_orth(std::sin(phi), -std::cos(phi),  0.0f),
          sphere_r,
          tube_r
      );

      torus.col = color;
      torus.n_ring_atoms = n_ring_atoms;

      mesh.add_submesh(coot::torus_mesh(torus));
   }

   mesh.add_submesh(make_data_sphere(3.99f, 128, 128));

   return mesh;
}


coot::simple_mesh_t make_radial_conformations() {
   coot::simple_mesh_t mesh;
   struct RingSpec { glm::vec3 dir; float cp_theta; float cp_phi; };
   const float pi            = static_cast<float>(M_PI);
   const float placement_r   = 5.0f;   // ring centres placed here
   const float cage_outer_r  = 4.00f;  // connector starts just outside cage
   const glm::vec4 atom_col({0.95f, 0.85f, 0.20f, 1.0f});  // gold atoms
   const glm::vec4 bond_col({0.75f, 0.65f, 0.15f, 1.0f});  // dark gold bonds
   const glm::vec4 line_color({0.75f, 0.75f, 0.75f, 1.0f});  // grey connectors

   const std::vector<RingSpec> placements = {
      // poles
      { { 0.0f,  0.0f,  1.0f }, 0.0f,       0.0f          },  // 4C1 chair
      { { 0.0f,  0.0f, -1.0f }, pi,          0.0f          },  // 1C4 chair
      // equatorial — boats
      { { 1.0f,  0.0f,  0.0f }, pi/2.0f,  0.0f             },  // B1,4
      { {-0.5f,  0.866f, 0.0f}, pi/2.0f,  2.0f*pi/3.0f     },  // B2,5
      { {-0.5f, -0.866f, 0.0f}, pi/2.0f,  4.0f*pi/3.0f     },  // B3,6
      // equatorial — twist-boats
      { { 0.5f,  0.866f, 0.0f}, pi/2.0f,  pi/3.0f          },  // 1S3
      { {-1.0f,  0.0f,  0.0f }, pi/2.0f,  pi               },  // 2S6
      { { 0.5f, -0.866f, 0.0f}, pi/2.0f,  5.0f*pi/3.0f     },  // 3S1
   };

   for (const auto &s : placements) {
      glm::vec3 ring_center = s.dir * placement_r;
      mesh.add_submesh(make_ring(s.cp_theta, s.cp_phi, ring_center, s.dir,
                                  0.8, atom_col, bond_col));

      glm::vec3 start = s.dir * cage_outer_r;
      glm::vec3 end   = ring_center;
      constexpr float radius = 0.02f;
      constexpr int slices = 16;
      mesh.add_submesh(make_stick(start, end, radius, line_color, slices));
   }

   return mesh;

}

coot::simple_mesh_t make_mesh() {
   coot::simple_mesh_t mesh;

   mesh.add_submesh(make_cage());
   mesh.add_submesh(make_radial_conformations());
   return mesh;
}