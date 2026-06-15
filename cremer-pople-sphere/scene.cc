//
// Created by Jordan Dialpuri on 04/04/2026.
//


#include "scene.hh"
#include "data_sphere.hh"
#include "density_map.hh"
#include "primitives.hh"
#include "ring.hh"
#include "coot-utils/shape-types.hh"
#include "coot-utils/shapes.hh"

constexpr float SPHERE_R = 100.0f;

coot::colour_holder to_colour(const glm::vec4 &c) {
   return coot::colour_holder(c.r, c.g, c.b, c.a);
}

std::string find_bin_file() {
   return "/data/coot/data/carbohydrate/n-glycan-validated-rscc.bin";
}

coot::simple_mesh_t make_cage() {
   coot::simple_mesh_t mesh;
   mesh.set_name("Cremer-Pople Cage");

   constexpr float sphere_r = SPHERE_R;
   constexpr float tube_r    = 0.006f * SPHERE_R;

   constexpr int n_lat = 7;
   constexpr int n_lon = 6;
   constexpr int n_ring_atoms = 2;

   const coot::colour_holder color = {1.0f, 1.0f, 1.0f, 0.5f};

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

   constexpr int default_lat = 128;
   constexpr int default_lon = 128;

   std::string density_path = find_bin_file();
   if (!density_path.empty()) {
      std::vector<float> density;
      int used_lat = default_lat;
      int used_lon = default_lon;

      density = load_density_grid_binary(density_path, used_lat, used_lon);
      
      mesh.add_submesh(make_data_sphere(0.99f * SPHERE_R, used_lat, used_lon, density));
      // mesh.add_submesh(density.empty()
         // ? make_data_sphere(3.99f, default_lat, default_lon)
         // : make_data_sphere(3.99f, used_lat, used_lon, density));
   } else {
      // mesh.add_submesh(make_data_sphere(3.99f, default_lat, default_lon));
   }

   return mesh;
}


coot::simple_mesh_t make_radial_conformations() {
   coot::simple_mesh_t mesh;
   struct RingSpec { glm::vec3 dir; float cp_theta; float cp_phi; };
   const float pi            = static_cast<float>(M_PI);
   const float placement_r   = 1.50f * SPHERE_R;  // ring centres placed here
   const float cage_outer_r  = 1.00f * SPHERE_R;  // connector starts just outside cage
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
                                  0.25f * SPHERE_R, atom_col, bond_col));

      glm::vec3 start = s.dir * cage_outer_r;
      glm::vec3 end   = ring_center;
      constexpr float radius = 0.005f * SPHERE_R;
      constexpr int slices = 16;
      mesh.add_submesh(make_stick(start, end, radius, line_color, slices));
   }

   return mesh;

}

/**
 * Constructs a composite mesh by combining submeshes generated from internal utility functions.
 *
 * The composite mesh is built by adding two submeshes: one representing a cage structure
 * and another representing radial conformations.
 *
 * @return A `coot::simple_mesh_t` object containing the combined mesh data.
 */
coot::simple_mesh_t make_mesh(bool add_radial_conformations) {
   coot::simple_mesh_t mesh;

   mesh.add_submesh(make_cage());
   if (add_radial_conformations)
      mesh.add_submesh(make_radial_conformations());
   return mesh;
}

/**
 * Creates a ring mesh centered at a point on a sphere, corresponding to the specified spherical coordinates.
 *
 * @param theta Polar angle in radians, representing the angle from the positive z-axis.
 * @param phi Azimuthal angle in radians, representing the angle from the positive x-axis in the xy-plane.
 *
 * @return A `coot::simple_mesh_t` object representing the generated ring mesh.
 */
coot::simple_mesh_t make_ring_at_sphere_click(float theta, float phi) {
   // Direction on the unit sphere corresponding to (theta, phi)
   const glm::vec3 dir(
       std::sin(theta) * std::cos(phi),
       std::sin(theta) * std::sin(phi),
       std::cos(theta)
   );

   // Place the ring slightly outside the cage surface
   constexpr float placement_r = 1.50f * SPHERE_R;
   const float cage_outer_r  = 1.00f * SPHERE_R;  // connector starts just outside cage

   const glm::vec3 center = dir * placement_r;

   const glm::vec4 atom_col(0.95f, 0.85f, 0.20f, 1.0f);  // gold
   const glm::vec4 bond_col(0.75f, 0.65f, 0.15f, 1.0f);

   coot::simple_mesh_t mesh;
   auto ring_mesh =  make_ring(theta, phi, center, dir, 0.2f * SPHERE_R, atom_col, bond_col);
   mesh.add_submesh(ring_mesh);

   glm::vec3 start = dir * cage_outer_r;
   glm::vec3 end   = dir * placement_r;
   constexpr float radius = 0.005f * SPHERE_R;
   constexpr int slices = 16;
   mesh.add_submesh(make_stick(start, end, radius, atom_col, slices));

   return mesh;
}



coot::simple_mesh_t make_pinpoint_at_sphere_point(float theta, float phi) {
   const glm::vec3 dir(
       std::sin(theta) * std::cos(phi),
       std::sin(theta) * std::sin(phi),
       std::cos(theta)
   );

   constexpr float base_r  = 1.00f * SPHERE_R;  // sphere surface
   constexpr float tip_r   = 1.10f * SPHERE_R;  // top of shaft
   constexpr float head_r  = 1.12f * SPHERE_R;  // centre of pin head

   const glm::vec4 shaft_col(0.85f, 0.85f, 0.85f, 1.0f);  // light grey shaft
   const glm::vec4 head_col (0.90f, 0.15f, 0.15f, 1.0f);  // red head

   coot::simple_mesh_t mesh;
   mesh.add_submesh(make_stick(dir * base_r, dir * tip_r, 0.005f * SPHERE_R, shaft_col, 8));
   mesh.add_submesh(make_atom_sphere(dir * head_r, 0.02f * SPHERE_R, head_col));
   return mesh;
}

coot::simple_mesh_t create_cremer_pople_sphere(const std::vector<CremerPopleData>& cremer_pople_data,
                                               bool add_radial_conformations) {
   coot::simple_mesh_t main_mesh;
   auto mesh = make_mesh(add_radial_conformations);
   main_mesh.add_submesh(mesh);

   for (const auto point : cremer_pople_data) {
      main_mesh.add_submesh(make_pinpoint_at_sphere_point(point.theta, point.phi));
   }

   return main_mesh;
}