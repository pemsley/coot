
#include <algorithm>
#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "Particle.hh"

void
Particle::update() {
   float v_scale = 0.009;
   v_scale = 0.015;
   glm::vec3 delta = v_scale * velocity;
   position += delta;
   velocity *= 0.992f;
   // colour.w *= 0.99;
   colour.g += 0.02;
   colour.r -= 0.02;
   colour.b += 0.001;
   life -= 0.18;
   float r = 0.5 * (1.0 + random());
   rotation += 0.01 * r; // currently a uniform is used, not this (which means they all spin at the same rate)
}

void
particle_container_t::remove_old_particles() {

   auto remover = [] (const Particle &p) {
                     return (p.life <= 0.0);
                  };

   particles.erase(std::remove_if(particles.begin(), particles.end(), remover), particles.end());

   if (! particles.empty()) {
      if (particles[0].life <= 0.0)
         particles.clear();
   }
      
}

float
particle_container_t::random() const {

   const float d = static_cast<float>(RAND_MAX);
   long int r = ::random();
   return static_cast<float>(r)/d;
}

void
particle_container_t::make_particles(unsigned int n_particles) {

   particles.clear();
   particles.reserve(n_particles);
   for (unsigned int i=0; i<n_particles; i++) {
      float p0 = 2.0 * random() - 1.0;
      float p1 = 2.0 * random() - 1.0;
      float p2 = 2.0 * random() - 1.0;
      while((p0*p0 + p1*p1) > 1.0) { // don't be "boxy"
          p0 = 2.0 * random() - 1.0;
          p1 = 2.0 * random() - 1.0;
      }
      float sc_pos = 0.1f;
      float sc_vel = 12.1f;
      glm::vec3 pos = sc_pos * glm::vec3(p0, p1, p2);
      glm::vec3 vel = sc_vel * glm::vec3(p0, p1, p2);
      glm::vec4 col(0.96, 0.26, 0.4, 1.0);
      if (false)
         std::cout << "Particle " << i << " " << glm::to_string(pos) << "\tvelocity "
                   << glm::to_string(vel) << " \t" << glm::to_string(col) << std::endl;
      Particle p(pos, vel, col, 10.0f - 9.0f * random());
      particles.push_back(p);
   }
}

// pass the time?
void
particle_container_t::update_particles() {

   for (unsigned int i=0; i<particles.size(); i++)
      particles[i].update();

   remove_old_particles();
}
