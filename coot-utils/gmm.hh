
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <random>
#include <algorithm> // shuffle

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

private:
    int k_;
    std::vector<Gaussian> components_;
    std::mt19937 rng_{std::random_device{}()};

    // Initializes the GMM parameters.
    void initialize(const std::vector<glm::vec3>& data) {
        std::cout << "Initializing GMM..." << std::endl;
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
    double gaussianPDF(const glm::vec3 &x, const glm::vec3 &mean, const glm::mat3 &cov) {

       // PDF formula: (1 / ((2*pi)^(d/2) * |Cov|^(1/2))) * exp(-1/2 * (x-mean)^T * Cov^-1 * (x-mean))
        // Here, d=3.
        try {
            glm::mat3 inv_cov = glm::inverse(cov);
            double det = glm::determinant(cov);
            if (det <= 0) {
                // If the determinant is non-positive, something went wrong.
               std::cout << "something went wrong:      det: " << det << std::endl;
               std::cout << "                           cov: " << glm::to_string(cov) << std::endl;
               std::cout << "                      inv_cov:  " << glm::to_string(inv_cov) << std::endl;
               return 0.0;
            }

            double exponent = -0.5 * glm::dot(x - mean, inv_cov * (x - mean));
            double normalization = 1.0 / (std::pow(2.0 * M_PI, 1.5) * std::sqrt(det));

            return normalization * std::exp(exponent);
        } catch (const std::exception& e) {
            std::cerr << "PDF calculation error: " << e.what() << std::endl;
            return 0.0;
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

