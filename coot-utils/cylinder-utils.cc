
#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "cylinder.hh"
#include "cylinder-utils.hh"


// make a "dashed" line from mini-cylinders between two points
dashed_cylinders_info_t get_dashed_cylinders(const std::vector<std::pair<glm::vec3, glm::vec3> > &positions,
                                             unsigned int n_dashes) {

   dashed_cylinders_info_t dci;
   dci.c = cylinder();

   for (unsigned int i=0; i<positions.size(); i++) {

      const auto &p1 = positions[i].first;
      const auto &p2 = positions[i].second;

      glm::vec3 delta = p2 - p1;
      glm::vec3 delta_uv = glm::normalize(delta);
      glm::vec3 z_dir(0,0,1);
      float d = glm::distance(p1,p2);

      cylinder c;
      c.init_unit(26);
      c.add_flat_end_cap();
      c.add_flat_start_cap();
      glm::mat4 ori = glm::orientation(delta_uv, z_dir);

      float length_per_dash = d / (2.0 * static_cast<float>(n_dashes));

      std::vector<glm::vec3> offsets;
      for (unsigned int i=0; i<n_dashes; i++) {
         glm::vec3 offset = p1 + delta_uv * 2.0f * length_per_dash * static_cast<float>(i);
         offsets.push_back(offset);
      }

      std::pair<glm::mat4, std::vector<glm::vec3> > o_o(ori, offsets);
      dci.oris_and_offsets.push_back(o_o);

   }

   float radius = 0.05;
   dci.scales = glm::vec3(radius, radius, 0.04);

   return dci;
}
