//
// Created by Jordan Dialpuri on 04/04/2026.
//

#include "density_map.hh"
#include <fstream>
#include <iostream>
#include <sstream>

std::vector<float> load_density_grid_binary(const std::string& filename,
                                             int& lat_out, int& lon_out) {
   std::ifstream f(filename, std::ios::binary);
   if (!f.is_open()) {
      std::cerr << "density_map: cannot open \"" << filename << "\"\n";
      return {};
   }

   // Magic
   char magic[8];
   f.read(magic, 8);
   if (!f || std::string(magic, 8) != std::string("CPDENS\0\0", 8)) {
      std::cerr << "density_map: bad magic in \"" << filename << "\"\n";
      return {};
   }

   // Version
   uint32_t version = 0;
   f.read(reinterpret_cast<char*>(&version), sizeof(version));
   if (!f || version != 1u) {
      std::cerr << "density_map: unsupported version " << version << "\n";
      return {};
   }

   // Grid dimensions
   int32_t lat = 0, lon = 0;
   f.read(reinterpret_cast<char*>(&lat), sizeof(lat));
   f.read(reinterpret_cast<char*>(&lon), sizeof(lon));
   if (!f || lat <= 0 || lon <= 0 || lat > 4096 || lon > 4096) {
      std::cerr << "density_map: implausible grid dimensions "
                << lat << "x" << lon << "\n";
      return {};
   }

   // Float payload
   std::vector<float> grid(lat * lon);
   f.read(reinterpret_cast<char*>(grid.data()),
          static_cast<std::streamsize>(lat * lon * sizeof(float)));
   if (!f) {
      std::cerr << "density_map: truncated data in \"" << filename << "\"\n";
      return {};
   }

   std::cout << "density_map: loaded precomputed " << lat << "x" << lon
             << " grid from \"" << filename << "\"\n";

   lat_out = static_cast<int>(lat);
   lon_out = static_cast<int>(lon);
   return grid;
}

