
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
   void update_gone_diego_particle();
   void update_gone_diff_map_particle();
};

class particle_container_t {
   float random() const;
public:
   std::vector<Particle> particles;
   void make_particles(unsigned int n_particles, const std::vector<glm::vec3> &positions);
   void make_gone_diego_particles(unsigned int n_particles_per_burst,
                                  const std::vector<glm::vec3> &positions,
                                  const glm::vec3 &screen_x_uv,
                                  const glm::vec3 &screen_y_uv); // usually just 1 or 2
   void make_gone_diff_map_peaks_particles(unsigned int n_particles_per_burst,
                                  const std::vector<std::pair<glm::vec3, float> > &positions,
                                  const glm::vec3 &screen_x_uv,
                                  const glm::vec3 &screen_y_uv); // usually just 1 or 2
   void update_particles();
   void update_gone_diego_particles();
   void update_gone_diff_map_particles();
   void remove_old_particles();
   unsigned int size() const { return particles.size(); }
   bool empty() const { return (particles.empty()); }
   bool have_particles_with_life() const;
   void clear() { particles.clear(); }

};

#endif
