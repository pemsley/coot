#ifndef COLOUR_FUNCTIONS_HH
#define COLOUR_FUNCTIONS_HH

#include <vector>

// colour helper function
// double combine_colour(double v, int col_part_index); 

std::vector<float> rotate_rgb(std::vector<float> &in_vals, float amount); 
std::vector<float> convert_rgb_to_hsv(const std::vector<float> &rgb);
std::vector<float> convert_hsv_to_rgb(const std::vector<float> &hsv);


#endif // COLOR_FUNCTIONS_HH
