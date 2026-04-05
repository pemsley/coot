//
// Created by Jordan Dialpuri on 04/04/2026.
//

#include "data.h"

#include <cmath>
#include <iostream>
#include <random>
#include "scene.hh"

int main(int argc, char **argv) {
   std::string output      = "output.glb";

   if (argc > 1) output      = argv[1];

   coot::simple_mesh_t main_mesh;
   auto mesh = make_mesh();
   main_mesh.add_submesh(mesh);

   constexpr float theta = 45 * static_cast<float>(M_PI) / 180;
   constexpr float phi = 0;
   auto current_selection_mesh = make_ring_at_sphere_click(theta, phi);
   main_mesh.add_submesh(current_selection_mesh);

   // Scatter random pinpoints uniformly across the sphere surface.
   // Uniform sampling requires theta = acos(U[-1,1]); plain U[0,π] clusters at the poles.
   std::mt19937 rng(42);
   std::uniform_real_distribution<float> dist_cos(-1.0f, 1.0f);
   std::uniform_real_distribution<float> dist_phi(0.0f, 2.0f * static_cast<float>(M_PI));

   constexpr int n_pinpoints = 50;
   for (int i = 0; i < n_pinpoints; ++i) {
      const float pin_theta = std::acos(dist_cos(rng));
      const float pin_phi   = dist_phi(rng);
      main_mesh.add_submesh(make_pinpoint_at_sphere_point(pin_theta, pin_phi));
   }

   main_mesh.export_to_gltf(output, 0.6f, 0.2f, true);
}


coot::simple_mesh_t create_cremer_pople_sphere(const std::vector<CremerPopleData>& cremer_pople_data) {
   coot::simple_mesh_t main_mesh;
   auto mesh = make_mesh();
   main_mesh.add_submesh(mesh);

   for (const auto point : cremer_pople_data) {
      main_mesh.add_submesh(make_pinpoint_at_sphere_point(point.theta, point.phi));
   }

   return main_mesh;
}