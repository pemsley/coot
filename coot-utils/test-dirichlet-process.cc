
#include <iostream>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

#include "utils/coot-utils.hh"
#include "dirichlet-process.hh"

std::vector<glm::vec3> make_test_points() {

   std::vector<glm::vec3> v;
   float disp = 5.0;
   for (unsigned int i=0; i<20; i++) {
      float x = 7.0f + coot::util::random_f() * disp;
      float y = 1.0f + coot::util::random_f() * disp;
      float z = 2.0f + coot::util::random_f() * disp;
      glm::vec3 p(x,y,z);
      v.push_back(p);
   }
   for (unsigned int i=0; i<20; i++) {
      float x = 2.0f + coot::util::random_f() * disp;
      float y = 7.0f + coot::util::random_f() * disp;
      float z = 2.0f + coot::util::random_f() * disp;
      glm::vec3 p(x,y,z);
      v.push_back(p);
   }
   for (unsigned int i=0; i<20; i++) {
      float x = 2.0f + coot::util::random_f() * disp;
      float y = 3.0f + coot::util::random_f() * disp;
      float z = 7.0f + coot::util::random_f() * disp;
      glm::vec3 p(x,y,z);
      v.push_back(p);
   }
   for (unsigned int i=0; i<4; i++) {
      float x = 2.0f + coot::util::random_f() * disp;
      float y = 2.0f + coot::util::random_f() * disp;
      float z = 2.0f + coot::util::random_f() * disp;
      glm::vec3 p(x,y,z);
      v.push_back(p);
   }

   return v;
}

int main(int argc, char **argv) {

   int status = 0;
   double alpha = 2.0;
   double beta  = 0.03;
   std::vector<glm::vec3> v = make_test_points();
   DirichletProcessClustering dpc(alpha, beta);
   std::vector<unsigned int> clustered_points = dpc.fit(v);

   for (unsigned int i=0; i<clustered_points.size(); i++) {
      const auto &idx = clustered_points[i];
      std::cout << " " << idx << " at "
                << v[i].x << " " << v[i].y << " " << v[i].z << std::endl;
   }

   return status;
}
