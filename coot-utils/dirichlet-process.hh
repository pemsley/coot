/* coot-utils/dirichlet-process.hh
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef COOT_UTILS_DIRICHLET_PROCESS
#define COOT_UTILS_DIRICHLET_PROCESS

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <limits>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/vec3.hpp>
#include <glm/gtx/norm.hpp>

class DirichletProcessClustering {

public:
   DirichletProcessClustering(double alpha, double beta)
      : alpha(alpha), beta(beta), generator(std::random_device{}()) {}

   std::vector<unsigned int> fit(const std::vector<glm::vec3>& data) {

      size_t n = data.size();
      std::vector<unsigned int> cluster_assignments(n, -1);
      std::vector<glm::vec3> cluster_means;
      std::vector<int> cluster_counts;

      std::cout << "size n " << n << std::endl;

      for (size_t i = 0; i < n; i++) {

         double total_probability = 0.0;
         std::vector<double> probabilities;

         // Calculate probabilities for existing clusters
         for (size_t k = 0; k < cluster_means.size(); k++) {
            double l = likelihood(data[i], cluster_means[k]);
            double prob = cluster_counts[k] * l;
            probabilities.push_back(prob);
            total_probability += prob;
            if (false)
               std::cout << "cluster-means-index-k " << k << " l " << l << " prob " << prob
                         << " running-total " << total_probability << std::endl;
         }

         // Calculate probability for a new cluster
         double l = likelihood_new_cluster(data[i]);
         double new_cluster_prob = alpha * l;
         probabilities.push_back(new_cluster_prob);
         total_probability += new_cluster_prob;

         // Normalize probabilities
         for (auto &prob : probabilities)
            prob /= total_probability;

         if (false)
            std::cout << "point-index-i " << i << " l " << l << " new_cluster_prob " << new_cluster_prob
                   << " total_probability " << total_probability << std::endl;

         if (false) {
            for (unsigned int ii=0; ii<probabilities.size(); ii++) {
               const auto &prob = probabilities[ii];
               std::cout << "   for point " << i << " normalized-prob " << ii << ": " << prob << std::endl;
            }
         }

         // Sample a cluster based on probabilities
         unsigned int sampled_cluster = sample_from_distribution(probabilities);
         // std::cout << "sample_from_distribution returned " << sampled_cluster << std::endl;
         if (sampled_cluster == cluster_means.size()) {
            // New cluster
            // std::cout << "new cluster " << std::endl;
            cluster_means.push_back(data[i]);
            cluster_counts.push_back(1);
         } else {
            // std::cout << "existing cluster " << sampled_cluster << std::endl;
            // Existing cluster
            glm::vec3 new_mean = update_mean(cluster_means[sampled_cluster], data[i], cluster_counts[sampled_cluster]);
            cluster_means[sampled_cluster] = new_mean;
            if (false)
               std::cout << "mean of cluster " << sampled_cluster << " updated to "
                         << new_mean.x << " " << new_mean.y << " " << new_mean.z << std::endl;
            ++cluster_counts[sampled_cluster];
         }

         cluster_assignments[i] = sampled_cluster;
      }

      return cluster_assignments;
   }

private:
   double alpha; // Concentration parameter
   double beta;  // Variance of the Gaussian likelihood
   std::default_random_engine generator;

   double likelihood(const glm::vec3& point, const glm::vec3& mean) {
      double dist_squared = glm::distance2(point, mean);
      double a = dist_squared / (2.0 * beta);
      const double lowest_double = std::numeric_limits<double>::denorm_min();
      if (a < 708.0)
         return std::exp(-a) / std::sqrt(2.0 * M_PI * beta);
      else
         return lowest_double;
   }

   double likelihood_new_cluster(const glm::vec3& point) {
      // Assume a zero-mean Gaussian prior for new clusters
      glm::vec3 zero_mean(0.0f);
      return likelihood(point, zero_mean);
   }

   unsigned int sample_from_distribution(const std::vector<double>& probabilities) {
      std::discrete_distribution<unsigned int> distribution(probabilities.begin(), probabilities.end());
      return distribution(generator);
   }

   glm::vec3 update_mean(const glm::vec3& current_mean, const glm::vec3& new_point, int current_count) {
      return (current_mean * static_cast<float>(current_count) + new_point) / static_cast<float>(current_count + 1);
   }
};

#endif // COOT_UTILS_DIRICHLET_PROCESS
