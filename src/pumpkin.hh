#ifndef PUMPKIN_HH
#define PUMPKIN_HH

#include <vector>
#include "generic-vertex.hh"
#include "coot-utils/g_triangle.hh"

std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > pumpkin();
std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > pumpkin_stalk();


#endif
