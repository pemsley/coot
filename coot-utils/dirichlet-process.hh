
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

   std::vector<int> fit(const std::vector<glm::vec3>& data) {
      size_t n = data.size();
      std::vector<int> cluster_assignments(n, -1);
      std::vector<glm::vec3> cluster_means;
      std::vector<int> cluster_counts;

      for (size_t i = 0; i < n; ++i) {
	 double total_probability = 0.0;
	 std::vector<double> probabilities;

	 // Calculate probabilities for existing clusters
	 for (size_t k = 0; k < cluster_means.size(); ++k) {
	    double prob = cluster_counts[k] * likelihood(data[i], cluster_means[k]);
	    probabilities.push_back(prob);
	    total_probability += prob;
	 }

	 // Calculate probability for a new cluster
	 double new_cluster_prob = alpha * likelihood_new_cluster(data[i]);
	 probabilities.push_back(new_cluster_prob);
	 total_probability += new_cluster_prob;

	 // Normalize probabilities
	 for (auto& prob : probabilities) {
	    prob /= total_probability;
	 }

	 // Sample a cluster based on probabilities
	 int sampled_cluster = sample_from_distribution(probabilities);
	 if (sampled_cluster == cluster_means.size()) {
	    // New cluster
	    cluster_means.push_back(data[i]);
	    cluster_counts.push_back(1);
	 } else {
	    // Existing cluster
	    cluster_means[sampled_cluster] = update_mean(cluster_means[sampled_cluster], data[i], cluster_counts[sampled_cluster]);
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
      return std::exp(-dist_squared / (2.0 * beta)) / std::sqrt(2.0 * M_PI * beta);
   }

   double likelihood_new_cluster(const glm::vec3& point) {
      // Assume a zero-mean Gaussian prior for new clusters
      glm::vec3 zero_mean(0.0f);
      return likelihood(point, zero_mean);
   }

   int sample_from_distribution(const std::vector<double>& probabilities) {
      std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
      return distribution(generator);
   }

   glm::vec3 update_mean(const glm::vec3& current_mean, const glm::vec3& new_point, int current_count) {
      return (current_mean * static_cast<float>(current_count) + new_point) / static_cast<float>(current_count + 1);
   }
};

