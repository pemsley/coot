//
// Created by Jordan Dialpuri on 04/04/2026.
//

#include <iostream>
#include "scene.hh"

int main(int argc, char **argv) {
   std::string output = "output.glb";
   if (argc > 1) output = argv[1];

   auto mesh = make_mesh();

   std::cout << "Mesh: " << mesh.vandt() << "\n";
   mesh.export_to_gltf(output, 0.6f, 0.2f, true);

   std::cout << "Written " << output << "\n";
}