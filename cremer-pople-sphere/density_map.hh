//
// Created by Jordan Dialpuri on 04/04/2026.
//

#ifndef COOT_DENSITY_MAP_H
#define COOT_DENSITY_MAP_H

#include <string>
#include <vector>

// Load a precomputed binary written by precompute_density.py.
// The grid dimensions are read from the file header and written to lat_out / lon_out.
// Returns an empty vector on any error.
std::vector<float> load_density_grid_binary(const std::string& filename,
                                             int& lat_out, int& lon_out);

// Parse a raw Cremer-Pople CSV (columns: q, phi_deg, theta_deg) and bin the
// conformations into a lat x lon grid.  Values are log1p-normalised to [0,1].
// Returns an empty vector if the file cannot be opened.
std::vector<float> load_density_grid_csv(const std::string& filename, int lat, int lon);

#endif  // COOT_DENSITY_MAP_H
