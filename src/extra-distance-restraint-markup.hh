#ifndef EXTRA_DISTANCE_RESTRAINT_MARKUP_HH
#define EXTRA_DISTANCE_RESTRAINT_MARKUP_HH

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

// this is converted from the bonds vector in extra-restraints-representation.hh

class extra_distance_restraint_markup_instancing_data_t {
public:
   float width;
   float length;
   glm::vec3 position; // position of the base (which will be at an atom)
   glm::mat3 orientation;
   glm::vec4 colour;
};

#endif // EXTRA_DISTANCE_RESTRAINT_MARKUP_HH
