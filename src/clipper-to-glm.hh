#ifndef CLIPPER_TO_GLM_HH
#define CLIPPER_TO_GLM_HH

#include <glm/glm.hpp>

#include <clipper/core/coords.h>

glm::vec3 clipper_to_glm(const clipper::Coord_orth &co);

#endif // CLIPPER_TO_GLM_HH
