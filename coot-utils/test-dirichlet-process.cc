
#include <iostream>
#include <set>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

#include "utils/coot-utils.hh"
#include "dirichlet-process.hh"

std::vector<glm::vec3> make_test_points() {

   std::vector<glm::vec3> v;
   float disp = 2.8; // 3.5 looks "real"
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
   for (unsigned int i=0; i<8; i++) {
      float x = 2.0f + coot::util::random_f() * disp;
      float y = 2.0f + coot::util::random_f() * disp;
      float z = 2.0f + coot::util::random_f() * disp;
      glm::vec3 p(x,y,z);
      v.push_back(p);
   }

   return v;
}

unsigned int do_simple(const std::vector<glm::vec3> &test_points, double alpha, double beta) {

   // double alpha = 0.1;
   // double beta  = 0.01;

   if (true) {
      for (unsigned int i=0; i<test_points.size(); i++) {
         std::cout << "input point " << i << " " << glm::to_string(test_points[i]) << std::endl;
      }
   }

   DirichletProcessClustering dpc(alpha, beta);
   std::vector<unsigned int> clustered_points = dpc.fit(test_points);

   for (unsigned int i=0; i<clustered_points.size(); i++) {
      const auto &idx = clustered_points[i];
      std::cout << " output cluster " << idx << " at "
                << test_points[i].x << " " << test_points[i].y << " " << test_points[i].z
                << std::endl;
   }

   std::set<unsigned int> cluster_set;
   for (unsigned int i=0; i<clustered_points.size(); i++) {
      const auto &idx = clustered_points[i];
      cluster_set.insert(idx);
   }

   return cluster_set.size();

}

void do_ranges() {

   std::vector<glm::vec3> test_points = make_test_points();

   std::vector<double> exp_values;
   for (int i=-200; i<200; i++) {
      double v = static_cast<double>(i)/2.0;
      exp_values.push_back(v);
   }

   for (double aa : exp_values) {
      for (double bb : exp_values ) {
         double alpha = std::pow(1.1, aa);
         double beta  = std::pow(1.1, bb);
         unsigned int n_clusters = do_simple(test_points, alpha, beta);
         std::cout << "result: alpha " << alpha << " beta " << beta
                   << " aa " << aa << " bb " << bb
                   << " n-clusters: " << n_clusters << std::endl;
      }
   }
}

int main(int argc, char **argv) {

   int status = 0;

   double alpha = 189.0;   // these values are in the "middle" of the "robust" zone
   double beta  = 0.085;

   std::vector<glm::vec3> test_points = make_test_points();

   bool output_ranges = false;
   if (argc > 1) {
      std::string argv1 = argv[1];
      if (argv1 == "ranges")
         output_ranges = true;
   }

   if (output_ranges)
      do_ranges();
   else
      do_simple(test_points, alpha, beta);

   return status;
}
