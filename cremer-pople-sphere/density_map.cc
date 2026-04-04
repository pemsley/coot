//
// Created by Jordan Dialpuri on 04/04/2026.
//

#include "density_map.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>

// ---------------------------------------------------------------------------
// Binary loader
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// CSV loader (slow path / fallback)
// ---------------------------------------------------------------------------

std::vector<float> load_density_grid_csv(const std::string& filename, int lat, int lon) {
   std::vector<float> counts(lat * lon, 0.0f);

   std::ifstream file(filename);
   if (!file.is_open()) {
      std::cerr << "density_map: cannot open \"" << filename << "\"\n";
      return {};
   }

   // Skip header line
   std::string line;
   std::getline(file, line);

   int n_loaded = 0, n_skipped = 0;

   while (std::getline(file, line)) {
      std::istringstream ss(line);
      std::string sq, sphi, stheta;
      if (!std::getline(ss, sq,     ',') ||
          !std::getline(ss, sphi,   ',') ||
          !std::getline(ss, stheta, ','))
         continue;

      float phi_deg   = std::stof(sphi);
      float theta_deg = std::stof(stheta);

      if (theta_deg < 0.0f || phi_deg < 0.0f) { ++n_skipped; continue; }

      int i = static_cast<int>(theta_deg / 180.0f * static_cast<float>(lat));
      int j = static_cast<int>(phi_deg   / 360.0f * static_cast<float>(lon));
      i = std::clamp(i, 0, lat - 1);
      j = std::clamp(j, 0, lon - 1);

      counts[i * lon + j] += 1.0f;
      ++n_loaded;
   }

   std::cout << "density_map: " << n_loaded << " conformations loaded from CSV, "
             << n_skipped << " skipped (invalid)\n";

   float max_count = *std::max_element(counts.begin(), counts.end());
   if (max_count > 0.0f) {
      const float log_max = std::log1p(max_count);
      for (auto& c : counts)
         c = std::log1p(c) / log_max;
   }

   return counts;
}
