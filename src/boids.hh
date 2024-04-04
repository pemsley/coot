/*
 * src/boids.hh
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <vector>
#include <glm/glm.hpp>

namespace fun {

   class boid {
      glm::vec3 velocity_delta_alignment(const std::vector<boid> &bonds,
                                         const std::vector<unsigned int> &seeable_other_boids) const;
      glm::vec3 velocity_delta_cohesion(const std::vector<boid> &bonds,
                                        const std::vector<unsigned int> &seeable_other_boids) const;
      glm::vec3 velocity_delta_no_bumps(const std::vector<boid> &bonds,
                                        const std::vector<unsigned int> &seeable_other_boids) const;
      glm::vec3 velocity_delta_no_bumps_in_objects(float box_lim) const;

   public:
      boid() { index = -1; }
      boid(int idx) { index = idx; }
      boid(int idx, const glm::vec3 &position, const glm::vec3 &velocity, const glm::vec3 &c) :
         index(idx), position(position), velocity(velocity), colour(c) {}
      unsigned int index;
      glm::vec3 position;
      glm::vec3 velocity;
      glm::vec3 colour;
      // wings are spread during turning, decrease velocity and (slow) gliding
      // If you are going to change the wing spreadness, it needs to be done
      // with instancing. (just changing 3 coordinate (maybe 2?) for each boid
      float wing_spreadness; // in radians

      glm::vec3 calc_velocity_delta(const std::vector<boid> &bonds,
                                    const std::vector<unsigned int> &seeable_other_boids,
                                    float box_lim) const;
      void apply_velocity_delta(const glm::vec3 &velocity_delta, float time_step);
      glm::mat4 make_mat() const;
   };

   class boids_container_t {
      std::vector<boid> boids;
   public:
      boids_container_t() { boids_box_limit = 30.0; }
      float boids_box_limit;
      void make_boids(unsigned int n_boids);
      std::vector<unsigned int> get_seeable_other_boids(unsigned int idx_this_boid);
      unsigned int size() const { return boids.size(); }
      const boid &operator[](unsigned int i) { return boids[i]; }
      void update();
   };

}

