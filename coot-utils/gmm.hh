#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <random>
#include <algorithm> // shuffle
#include <limits>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp> // For printing vectors and matrices
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtx/norm.hpp>

// A small constant to prevent division by zero or log of zero
const double EPSILON = 1e-6;

// A single Gaussian component in the mixture model.
struct Gaussian {
   double weight;
   glm::vec3 mean;
   glm::mat3 covariance;
};

// The Gaussian Mixture Model class that performs the EM algorithm.
class GMM {
public:
   // Constructor.
   GMM(int k) : k_(k) {
       components_.resize(k_);
   }

   // The main function to fit the model to the data.
   void fit(const std::vector<glm::vec3>& data, int max_iterations = 100) {
      if (data.empty()) {
          std::cerr << "Error: Data vector is empty." << std::endl;
          return;
      }

      initialize(data);

      for (int i = 0; i < max_iterations; ++i) {

         // std::cout << "Iteration " << i + 1 << std::endl;

         // E-step: Calculate responsibilities
         std::vector<std::vector<double>> responsibilities = expectationStep(data);
         // M-step: Update parameters based on responsibilities
         maximizationStep(data, responsibilities);
      }
   }

   // Compute the log-likelihood of the data given the current model
   double computeLogLikelihood(const std::vector<glm::vec3>& data) const {

      double log_likelihood = 0.0;
      for (const auto& point : data) {
         double point_likelihood = 0.0;
         for (int j = 0; j < k_; ++j) {
            double prob = gaussianPDF(point, components_[j].mean, components_[j].covariance);
            point_likelihood += components_[j].weight * prob;
         }

         if (point_likelihood > EPSILON) {
            log_likelihood += std::log(point_likelihood);
         } else {
            // Prevent log(0) by using a small value
            log_likelihood += std::log(EPSILON);
         }
      }

      return log_likelihood;
   }

   // Compute the Bayesian Information Criterion (BIC)
   // BIC = -2 * log(L) + k * log(n)
   // where L is the likelihood, k is the number of parameters, n is the number of data points
   // Lower BIC is better
   double computeBIC(const std::vector<glm::vec3>& data) const {
      double log_likelihood = computeLogLikelihood(data);
      int n = data.size();

      // Number of free parameters:
      // For each component:
      //   - 1 weight (but weights sum to 1, so k-1 free parameters)
      //   - 3 mean parameters
      //   - 6 unique covariance parameters (3x3 symmetric matrix has 6 unique values)
      // Total: (k-1) + k*3 + k*6 = k*9 + k - 1 = 10k - 1
      int num_parameters = 10 * k_ - 1;
      double bic = -2.0 * log_likelihood + num_parameters * std::log(n);
      return bic;
   }

   // Assign each data point to its most likely cluster
   std::vector<int> predict(const std::vector<glm::vec3>& data) const {

      std::vector<int> assignments(data.size());
      for (size_t i = 0; i < data.size(); ++i) {
         double max_responsibility = -1.0;
         int best_cluster = 0;

         // Compute responsibility for each component
         double responsibility_sum = 0.0;
         std::vector<double> responsibilities(k_);
         for (int j = 0; j < k_; ++j) {
            double prob = gaussianPDF(data[i], components_[j].mean, components_[j].covariance);
            responsibilities[j] = components_[j].weight * prob;
            responsibility_sum += responsibilities[j];
         }
         // Normalize and find the maximum
         for (int j = 0; j < k_; ++j) {
            if (responsibility_sum > EPSILON) {
               responsibilities[j] /= responsibility_sum;
            }
            if (responsibilities[j] > max_responsibility) {
               max_responsibility = responsibilities[j];
               best_cluster = j;
            }
         }
         assignments[i] = best_cluster;
      }
      return assignments;
   }

   // Function to print the learned parameters of the GMM.
   void printParameters(const std::string &lab) const {

      std::cout << "\n--- " << lab << " GMM Parameters ---" << std::endl;
      for (int i = 0; i < k_; ++i) {
         std::cout << "Component " << i + 1 << ":" << std::endl;
         std::cout << "  Weight: " << components_[i].weight << std::endl;
         std::cout << "  Mean: " << glm::to_string(components_[i].mean) << std::endl;
         std::cout << "  Covariance: " << glm::to_string(components_[i].covariance) << std::endl;
      }
   }

   int getK() const { return k_; }

private:
    int k_;
    std::vector<Gaussian> components_;
    std::mt19937 rng_{std::random_device{}()};

    // Initializes the GMM parameters.
    void initialize(const std::vector<glm::vec3>& data) {
        // std::cout << "Initializing GMM with k=" << k_ << "..." << std::endl;
        // Simple initialization:
        // - Randomly pick `k` data points as initial means.
        // - Set covariance to identity matrix.
        // - Set weights equally.
        std::shuffle(components_.begin(), components_.end(), rng_);

        std::vector<int> random_indices(data.size());
        std::iota(random_indices.begin(), random_indices.end(), 0);
        std::shuffle(random_indices.begin(), random_indices.end(), rng_);

        for (int i = 0; i < k_; ++i) {
            components_[i].weight = 1.0 / k_;
            components_[i].mean = data[random_indices[i]];
            components_[i].covariance = glm::mat3(1.0f); // Identity matrix
        }
    }

    // Calculates the probability density function (PDF) for a multivariate Gaussian.
    double gaussianPDF(const glm::vec3 &x, const glm::vec3 &mean, const glm::mat3 &cov) const {

       // PDF formula: (1 / ((2*pi)^(d/2) * |Cov|^(1/2))) * exp(-1/2 * (x-mean)^T * Cov^-1 * (x-mean))
        // Here, d=3.
        try {
            glm::mat3 inv_cov = glm::inverse(cov);
            double det = glm::determinant(cov);
            if (det <= 0) {
                // If the determinant is non-positive, something went wrong.
               // std::cout << "something went wrong:      det: " << det << std::endl;
               // std::cout << "                           cov: " << glm::to_string(cov) << std::endl;
               // std::cout << "                      inv_cov:  " << glm::to_string(inv_cov) << std::endl;
               return EPSILON; // Return small value instead of 0
            }

            double exponent = -0.5 * glm::dot(x - mean, inv_cov * (x - mean));
            double normalization = 1.0 / (std::pow(2.0 * M_PI, 1.5) * std::sqrt(det));

            return normalization * std::exp(exponent);
        } catch (const std::exception& e) {
            std::cerr << "PDF calculation error: " << e.what() << std::endl;
            return EPSILON;
        }
    }

   // The Expectation step: calculates the responsibilities of each component.
   std::vector<std::vector<double>> expectationStep(const std::vector<glm::vec3>& data) {

      std::vector<std::vector<double>> responsibilities(data.size(), std::vector<double>(k_));

      for (size_t i = 0; i < data.size(); ++i) {
         double responsibility_sum = 0.0;
         for (int j = 0; j < k_; ++j) {
            double prob = gaussianPDF(data[i], components_[j].mean, components_[j].covariance);
            responsibilities[i][j] = components_[j].weight * prob;
            responsibility_sum += responsibilities[i][j];
         }

         // std::cout << "responsibility_sum " << responsibility_sum << std::endl;
         // Normalize responsibilities
         if (responsibility_sum > EPSILON) {
            for (int j = 0; j < k_; ++j) {
               responsibilities[i][j] /= responsibility_sum;
            }
         }
      }
      return responsibilities;
   }

   // The Maximization step: updates the parameters of each component.
   void maximizationStep(const std::vector<glm::vec3> &data, const std::vector<std::vector<double>>  &responsibilities) {

      for (int j = 0; j < k_; ++j) {
         double responsibility_sum_j = 0.0;
         glm::vec3 new_mean = glm::vec3(0.0f);
         for (size_t i = 0; i < data.size(); ++i) {
             responsibility_sum_j += responsibilities[i][j];
             double r_ij = responsibilities[i][j];
             glm::vec3 rd = data[i] * glm::vec3(r_ij, r_ij, r_ij);
             new_mean += rd;
         }

         if (responsibility_sum_j > EPSILON) {
             // Update weight
             components_[j].weight = responsibility_sum_j / data.size();

             // Update mean
             components_[j].mean = new_mean / (float)responsibility_sum_j;

             // Update covariance
             glm::mat3 new_cov(0.0f);
             for (size_t i = 0; i < data.size(); ++i) {
                 glm::vec3 diff = data[i] - components_[j].mean;
                 float r_ij = responsibilities[i][j]; // casting
                 glm::mat3 op = glm::outerProduct(diff, diff);
                 new_cov += r_ij * op;
             }

             // Final covariance update: Divide by sum of responsibilities
             glm::mat3 final_cov = new_cov / static_cast<float>(responsibility_sum_j);

             // Add regularization to the diagonal (prevents singularities)
             final_cov += glm::mat3(EPSILON);

             components_[j].covariance = final_cov;
         }
      }
   }
};

// Helper function to select optimal number of clusters using BIC
// Tests k values from k_min to k_max and returns the k with lowest BIC
struct GMMFitResult {
   int optimal_k;
   double best_bic;
   GMM best_model;
   std::vector<int> cluster_assignments;

   GMMFitResult(int k) : optimal_k(k), best_bic(std::numeric_limits<double>::max()), best_model(k) {}
};

inline GMMFitResult fitGMMWithBIC(const std::vector<glm::vec3>& data,
                                   int k_min = 1,
                                   int k_max = 10,
                                   int max_iterations = 100,
                                   bool verbose = true) {

   if (data.empty()) {
      std::cerr << "Error: Cannot fit GMM to empty data" << std::endl;
      return GMMFitResult(1);
   }

   // Ensure k_max doesn't exceed the number of data points
   k_max = std::min(k_max, static_cast<int>(data.size()));
   k_min = std::max(1, k_min);

   if (k_min > k_max) {
      std::cerr << "Error: k_min > k_max" << std::endl;
      return GMMFitResult(1);
   }

   GMMFitResult result(k_min);

   if (verbose) {
      std::cout << "\n=== BIC Model Selection ===" << std::endl;
      std::cout << "Testing k from " << k_min << " to " << k_max << std::endl;
   }

   for (int k = k_min; k <= k_max; ++k) {
      GMM gmm(k);
      gmm.fit(data, max_iterations);

      double bic = gmm.computeBIC(data);
      double log_likelihood = gmm.computeLogLikelihood(data);

      if (verbose) {
         std::cout << "k=" << k
                   << " | BIC=" << std::fixed << bic
                   << " | Log-Likelihood=" << log_likelihood << std::endl;
      }

      if (bic < result.best_bic) {
         result.best_bic = bic;
         result.optimal_k = k;
         result.best_model = gmm; // This will copy the model
         result.cluster_assignments = gmm.predict(data);
      }
   }

   if (verbose) {
      std::cout << "\nOptimal number of clusters: " << result.optimal_k
                << " (BIC=" << std::fixed << result.best_bic << ")" << std::endl;
      result.best_model.printParameters("Optimal Model");
   }

   return result;
}
