/*
 * src/colour-functions.hh
 *
 * Copyright 2022 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef COLOUR_FUNCTIONS_HH
#define COLOUR_FUNCTIONS_HH

#include <vector>

// colour helper function
// double combine_colour(double v, int col_part_index); 

std::vector<float> rotate_rgb(std::vector<float> &in_vals, float amount); 
std::vector<float> convert_rgb_to_hsv(const std::vector<float> &rgb);
std::vector<float> convert_hsv_to_rgb(const std::vector<float> &hsv);


#endif // COLOR_FUNCTIONS_HH
