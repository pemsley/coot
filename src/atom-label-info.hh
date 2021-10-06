
#ifndef ATOM_LABEL_INFO_HH
#define ATOM_LABEL_INFO_HH

#include <string>
#include <glm/glm.hpp>

// not really an atom label.

class atom_label_info_t {
public:
   std::string label;
   glm::vec3 position;
   glm::vec4 colour;
   atom_label_info_t(const std::string l, const glm::vec3 &p, const glm::vec4 c) : label(l), position(p), colour(c) {}
};


#endif // ATOM_LABEL_INFO_HH

