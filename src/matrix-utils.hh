
#ifndef COOT_MATRIX_UTILS_HH
#define COOT_MATRIX_UTILS_HH

#include <clipper/core/coords.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "coot-utils/coot-coord-utils.hh"

clipper::Mat33<double> glm_to_mat33(const glm::mat4 &quat_mat);

glm::quat coot_quaternion_to_glm(const coot::util::quaternion &q);

coot::util::quaternion glm_to_coot_quaternion(const glm::quat &q);

#endif // COOT_MATRIX_UTILS_HH
