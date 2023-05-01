
#ifndef PRIDEOUT_OCTASPHERE_HH
#define PRIDEOUT_OCTASPHERE_HH

#include <vector>

#include <glm/glm.hpp>

#include "coot-utils/g_triangle.hh"

std::pair<std::vector<glm::vec3>, std::vector<g_triangle> >
tessellate_octasphere_patch(unsigned int num_subdivisions);



#endif // PRIDEOUT_OCTASPHERE_HH
