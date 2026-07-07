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



#endif  // COOT_DENSITY_MAP_H
