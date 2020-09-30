
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
   Particle(const glm::vec3 &p, const glm::vec3 &v, const glm::vec4 &c, float l) :
      position(p), velocity(v), colour(c), life(l) { opacity = 1.0f; rotation = 0.0f; }
   // update the position, velocity, colour and life
   void update();
};

class particle_container_t {
   float random() const;
public:
   std::vector<Particle> particles;
   void make_particles(unsigned int n_particles);
   void update_particles();
   void remove_old_particles();
   void create_particle(const glm::vec3 &p);
   unsigned int size() const { return particles.size(); }
   bool empty() const { return (particles.size() == 0); }

};

#endif
