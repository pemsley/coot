
#include <iostream>
#include <vector>
#include <random>

// Including the GLM headers for vec3 and mat3, and necessary functions.
// Note: GLM is a header-only library, so you typically just need to include the headers.
// You might need to add -I<path_to_glm> to your compiler command.
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp> // For printing vectors and matrices
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtx/norm.hpp>

#include "gmm.hh"

int main(int argc, char **argv) {

   // We'll generate some synthetic data from a known GMM to test the algorithm.
    const int num_points = 200;
    const int num_clusters = 3;
    std::vector<glm::vec3> data;
    std::mt19937 generator(42); // Seeded for reproducibility
    std::normal_distribution<double> dist_x(0, 1.0);
    std::normal_distribution<double> dist_y(0, 1.0);
    std::normal_distribution<double> dist_z(0, 1.0);

    std::vector<glm::vec3> true_means = {
        glm::vec3( 5,  5,  5),
        glm::vec3(-5,  0, -5),
        glm::vec3( 0, -5,  5)
    };

    for (int i = 0; i < num_points; ++i) {
        int cluster_idx = i % num_clusters;
        glm::vec3 point;
        point.x = dist_x(generator) + true_means[cluster_idx].x;
        point.y = dist_y(generator) + true_means[cluster_idx].y;
        point.z = dist_z(generator) + true_means[cluster_idx].z;
        std::cout << "point " << i << " " << point.x << " " << point.y << " " << point.z << std::endl;
        data.push_back(point);
    }

    std::cout << "Generated " << data.size() << " data points." << std::endl;
    std::cout << "True means are: "
              << glm::to_string(true_means[0]) << ", "
              << glm::to_string(true_means[1]) << ", "
              << glm::to_string(true_means[2]) << std::endl;

    // Create and fit the GMM.
    GMM gmm(num_clusters);
    gmm.fit(data, 200);

    // Print the learned parameters.
    gmm.printParameters("Final");

    return 0;
}

