
#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <vector>
#include <glm/ext.hpp>

class Particle {
public:
   glm::vec3 position;
   glm::vec3 velocity;
   glm::vec4 colour;
   float life;
   float rotation; // radians
   float opacity;
   float colour_change_rate;
   Particle(const glm::vec3 &p, const glm::vec3 &v, const glm::vec4 &c, float l) :
      position(p), velocity(v), colour(c), life(l), colour_change_rate(1.0) {
      opacity = 1.0f; rotation = 0.0f; }
   // update the position, velocity, colour and life
   void update();
};

class particle_container_t {
   float random() const;
public:
   std::vector<Particle> particles;
   void make_particles(unsigned int n_particles, const std::vector<glm::vec3> &positions);
   void make_gone_diego_particles(const std::vector<glm::vec3> &positions); // usually just 1 or 2
   void update_particles();
   void remove_old_particles();
   unsigned int size() const { return particles.size(); }
   bool empty() const { return (particles.empty()); }

};

#endif
