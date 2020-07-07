
#include <algorithm>
#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "Particle.hh"

void
Particle::update() {
   float v_scale = 0.009;
   glm::vec3 delta = v_scale * velocity;
   position += delta;
   velocity *= 0.994f;
   // colour.w *= 0.99;
   colour.g += 0.01;
   colour.r -= 0.002;
   colour.b += 0.004;
   life -= 0.18;
}

void
particle_container_t::remove_old_particles() {

   auto remover = [] (const Particle &p) {
                     return (p.life <= 0.0);
                  };

   particles.erase(std::remove_if(particles.begin(), particles.end(), remover), particles.end());
}

float
particle_container_t::random() const {

   long int r = ::random();
   return static_cast<float>(r)/static_cast<float>(RAND_MAX);
}

void
particle_container_t::make_particles(unsigned int n_particles, const glm::vec3 &rotation_centre) {

   std::cout << "making particles with rotation_centre " << glm::to_string(rotation_centre) << std::endl;

   for (unsigned int i=0; i<n_particles; i++) {
      float p0 = 2.0 * random() - 1.0;
      float p1 = 2.0 * random() - 1.0;
      float p2 = 2.0 * random() - 1.0;
      while((p0*p0 + p1*p1) > 1.0) { // don't be "boxy"
          p0 = 2.0 * random() - 1.0;
          p1 = 2.0 * random() - 1.0;
      }
      float s_pos = 1.92;
      float s_vel = 11.1;
      glm::vec3 pos(s_pos * p0, s_pos * p1, s_pos * p2);
      pos += rotation_centre;
      float v0 = p0;
      float v1 = p1;
      float v2 = p2;
      glm::vec4 col(0.66, 0.66, 0.2, 1.0);
      glm::vec3 vel(s_vel * v0, s_vel * v1, s_vel * v2);
      std::cout << "Particle " << i << " " << glm::to_string(pos) << "\tvelocity " << glm::to_string(vel) << " \t"
                << glm::to_string(col) << std::endl;
      Particle p(pos, vel, col, 10.0f - 9.0f * random());
      particles.push_back(p);
   }
}

// pass the time?
void
particle_container_t::update_particles() {

   glm::vec3 prev;
   for (unsigned int i=0; i<particles.size(); i++) {
      particles[i].update();
   }
   remove_old_particles();
}
