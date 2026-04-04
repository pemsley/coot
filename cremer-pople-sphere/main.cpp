//
// Created by Jordan Dialpuri on 04/04/2026.
//

#include <iostream>
#include "scene.hh"

int main(int argc, char **argv) {
   std::string output      = "output.glb";
   std::string density_csv = "";

   if (argc > 1) output      = argv[1];
   if (argc > 2) density_csv = argv[2];

   coot::simple_mesh_t main_mesh;
   auto mesh = make_mesh(density_csv);
   main_mesh.add_submesh(mesh);

   constexpr float theta = 45 * static_cast<float>(M_PI) / 180;
   constexpr float phi = 0;
   auto current_selection_mesh = make_ring_at_sphere_click(theta, phi);
   main_mesh.add_submesh(current_selection_mesh);
   main_mesh.export_to_gltf(output, 0.6f, 0.2f, true);
}