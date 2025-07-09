/*
 * src/boids.cc
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

#include <iostream>
#include <iomanip>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/rotate_vector.hpp>

#include "boids.hh"


glm::vec3
fun::boid::velocity_delta_alignment(const std::vector<boid> &boids,
                                    const std::vector<unsigned int> &seeable_other_boids) const {

   glm::vec3 r(0,0,0);

   if (! seeable_other_boids.empty()) {
      for (unsigned int i=0; i<seeable_other_boids.size(); i++)
         r += boids[seeable_other_boids[i]].velocity;
      r *= 1.0f/static_cast<float>(seeable_other_boids.size());
   }
   return 0.5f * r;
}

glm::vec3
fun::boid::velocity_delta_cohesion(const std::vector<boid> &boids,
                                   const std::vector<unsigned int> &seeable_other_boids) const {

   glm::vec3 mid(0,0,0);
   // std::cout << "cohesion seeable others: " << seeable_other_boids.size() << std::endl;
   if (! seeable_other_boids.empty()) {
      float sf = 1.0f/static_cast<float>(seeable_other_boids.size());
      for (unsigned int i=0; i<seeable_other_boids.size(); i++) {
         mid += boids[seeable_other_boids[i]].position;
      }
      mid *= sf;
   }
   glm::vec3 delta_position = mid - position;
   glm::vec3 r = 2.2f * delta_position;
   std::cout << "cohesion-delta: " << r.x << " " << r.y << " " << r.z << std::endl;
   return r;
}


glm::vec3
fun::boid::velocity_delta_no_bumps(const std::vector<boid> &boids,
                                   const std::vector<unsigned int> &seeable_other_boids) const {

   const float too_close_limit = 5.0;
   glm::vec3 delta_v(0,0,0);

   if (! seeable_other_boids.empty()) {
      for (unsigned int i=0; i<seeable_other_boids.size(); i++) {
         const glm::vec3 &position_for_other_boid = boids[seeable_other_boids[i]].position;
         glm::vec3 delta_position = position_for_other_boid - position;
         float dd = glm::distance2(delta_position, glm::vec3(0,0,0));
         // std::cout << "dd: " << dd << std::endl;
         if (dd < too_close_limit * too_close_limit) {
            float weight = too_close_limit - sqrt(dd);
            delta_v += - weight * delta_position;
         }
      }
   }

   // std::cout << "velocity_delta_no_bumps returns " << glm::to_string(delta_v) << std::endl;
   return delta_v * 0.501f;
}

glm::vec3
fun::boid::velocity_delta_no_bumps_in_objects(float box_lim) const {

   glm::vec3 delta_v(0,0,0);
   for (unsigned int i=0; i<3; i++) {
      float pp = position[i]; // x y or z
      float dp = pp - box_lim;
      float dm = pp + box_lim;
      float sf = -0.1;
      delta_v[i] += sf * dp;
      delta_v[i] += sf * dm;
   }
   return delta_v;
}

#include "utils/coot-utils.hh"

glm::vec3
fun::boid::calc_velocity_delta(const std::vector<boid> &boids,
                               const std::vector<unsigned int> &seeable_other_boids,
                               float box_lim) const {

   float p0 = 2.0 * coot::util::random_f() - 1.0;
   float p1 = 2.0 * coot::util::random_f() - 1.0;
   float p2 = 2.0 * coot::util::random_f() - 1.0;

   glm::vec3 random_v(p0, p1, p2);
   glm::vec3 vd1 = velocity_delta_no_bumps_in_objects(box_lim);
   glm::vec3 vd2 = velocity_delta_no_bumps(boids,  seeable_other_boids);
   glm::vec3 vd3 = velocity_delta_cohesion(boids,  seeable_other_boids);
   glm::vec3 vd4 = velocity_delta_alignment(boids, seeable_other_boids);

   glm::vec3 delta_v = vd1 + vd3 + vd4;

   std::cout << "compare"
	     << " " << vd1.x << " " << vd1.y << " " << vd1.z
	     << " " << vd3.x << " " << vd3.y << " " << vd3.z
	     << std::endl;

   // float delta_v_dd = glm::distance2(delta_v, glm::vec3(0,0,0));
   // float max_d = 100000.0;
   // float max_dd = max_d * max_d;
   // if (delta_v_dd > max_dd) {
   //    glm::vec3 delta_v_uv = glm::normalize(delta_v);
   //    delta_v = max_d * delta_v_uv;
   // }

   float sff = 1.1f;
   // delta_v = glm::vec3(0.0f,0.0f,0.0f);
   return delta_v * sff;

}


void
fun::boids_container_t::make_boids(unsigned int n_boids) {

   for (unsigned int i=0; i<n_boids; i++) {
      glm::vec3 colour(0.4, 0.4, 0.6);
      float p0 = 2.0 * coot::util::random_f() - 1.0;
      float p1 = 2.0 * coot::util::random_f() - 1.0;
      float p2 = 2.0 * coot::util::random_f() - 1.0;
      float v0 = 2.0 * coot::util::random_f() - 1.0;
      float v1 = 2.0 * coot::util::random_f() - 1.0;
      float v2 = 2.0 * coot::util::random_f() - 1.0;
      glm::vec3 position(p0, p1, p2);
      glm::vec3 velocity(v0, v1, v2);
      boid b(i, 10.0f * position, 10.0f * velocity, colour);
      boids.push_back(b);
   }

   if (false) {
      for (unsigned int i=0; i<boids.size(); i++) {
         boid &boid = boids[i];
         std::cout << " make_boids() " << i << " "
                   << glm::to_string(boid.position) <<  " "
                   << glm::to_string(boid.velocity) <<  " "
                   << glm::to_string(boid.colour) << std::endl;
      }
   }
   
}

void
fun::boid::apply_velocity_delta(unsigned int idx, const glm::vec3 &velocity_delta, float time_step) {

   velocity += 0.01f * velocity_delta;

   // velocity normalization
   if (false) {
      float speed_sqrd = glm::distance2(velocity, glm::vec3(0,0,0));
      if (speed_sqrd < 1.0)
         velocity = glm::normalize(velocity);
      if (speed_sqrd > 25.0)
         velocity = 5.0f * glm::normalize(velocity);
   }

   // glm::vec3 delta_position = 0.005f * velocity;
   float sf = 0.01;
   glm::vec3 delta_position(velocity.x * sf, velocity.y * sf, velocity.z * sf);
   if (false)
      std::cout << "position-delta " << idx << " "
		<< delta_position.x << " "
		<< delta_position.y << " "
		<< delta_position.z << std::endl;
   position += delta_position;
   // std::cout << "position: " << idx << " " << position.x << " " << position.y << " " << position.z
   // << std::endl;

}


std::vector<unsigned int>
fun::boids_container_t::get_seeable_other_boids(unsigned int idx_this_boid) {

   float max_distance_to_another_boid = 400.0;
   float dd_max = max_distance_to_another_boid * max_distance_to_another_boid;

   const boid &this_boid = boids[idx_this_boid];

   std::vector<unsigned int> others;
   others.reserve(10);
   for (unsigned int i=0; i<boids.size(); i++) {
      if (i == idx_this_boid)
         continue;
      glm::vec3 delta_position = boids[i].position - this_boid.position;
      float dd_distance_to_other = glm::distance2(boids[i].position, this_boid.position);
      // std::cout << "distance to other " << sqrt(dd_distance_to_other) << std::endl;
      if (dd_distance_to_other < dd_max) {
         glm::vec3 velocity_normalized       = glm::normalize(this_boid.velocity);
         glm::vec3 delta_position_normalized = glm::normalize(delta_position);
         float dp_vel_delta_pos = glm::dot(delta_position, velocity_normalized);
         if (dp_vel_delta_pos > -0.5)
            others.push_back(i);
      }
   }

   //    std::cout << "for boid " << idx_this_boid
   // << " returning seeable_other_boids size " << others.size() << std::endl;

   // debugging
   others.clear();
   for (unsigned int i=0; i<boids.size(); i++) {
      if (i == idx_this_boid)
         continue;
      others.push_back(i);
   }

   return others;
}

void
fun::boids_container_t::update() {

   std::cout << "boids update()" << std::endl;

   float timestep = 0.02; // seconds

   if (false) {
      for (unsigned int i=0; i<boids.size(); i++) {
         const boid &boid = boids[i];
	 std::cout << i << " "
		   << glm::to_string(boid.position) <<  " "
		   << glm::to_string(boid.velocity) <<  " "
		   << glm::to_string(boid.colour) << std::endl;
      }
   }

   std::vector<glm::vec3> velocity_deltas(boids.size());
   std::vector<glm::vec3> velocities(boids.size());
   for (unsigned int i=0; i<boids.size(); i++) {
      boid &boid = boids[i];
      std::vector<unsigned int> other_boids = get_seeable_other_boids(i);
      glm::vec3 velocity_delta = boid.calc_velocity_delta(boids, other_boids, boids_box_limit);
      velocity_deltas[i] = velocity_delta;
      velocities[i] = boid.velocity;
   }

   // // sanitize the velocities
   // float sum_vel = 0;
   // for (unsigned int i=0; i<boids.size(); i++) {
   //    float d = glm::distance(velocities[i], glm::vec3(0,0,0));
   //    sum_vel += d;
   // }
   // float av_vel = sum_vel/static_cast<float>(boids.size());
   // float target_average = 1.0;
   // float vel_sf = target_average/av_vel;
   
   
   for (unsigned int i=0; i<boids.size(); i++) {
      boids[i].apply_velocity_delta(i, velocity_deltas[i], timestep);
   }

   // for (unsigned int i=0; i<boids.size(); i++) {
   //    boids[i].velocity *= vel_sf;
   // }

   if (false) {
      for (unsigned int i=0; i<boids.size(); i++) {
	 float dd = glm::distance2(velocity_deltas[i], glm::vec3(0,0,0));
	 std::cout << "boid accel " << std::setw(2) << i << " " << sqrt(dd)
		   << " " << std::setw(12) << velocity_deltas[i].x
		   << " " << std::setw(12) << velocity_deltas[i].y
		   << " " << std::setw(12) << velocity_deltas[i].z
		   << std::endl;
      }
   }

}

glm::mat4
fun::boid::make_mat() const {

   glm::vec3 normalized_vel = glm::normalize(velocity);
   glm::mat4 m = glm::transpose(glm::orientation(normalized_vel, glm::vec3(0.0, 0.0, 1.0)));
   // m = glm::mat4(1.0f); // remove orientation
   glm::mat4 mt = glm::translate(m, position);
   mt = m;
   return mt;
}
