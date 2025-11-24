
#include "cylinder.hh"

class dashed_cylinders_info_t {
public:
   dashed_cylinders_info_t() {}
   cylinder c;
   std::vector<std::pair<glm::mat4, std::vector<glm::vec3> > > oris_and_offsets;
   glm::vec3 scales;
};

// make a "dashed" line from mini-cylinders between two points
dashed_cylinders_info_t get_dashed_cylinders(const std::vector<std::pair<glm::vec3, glm::vec3> > &positions,
                                             unsigned int n_dashes);
