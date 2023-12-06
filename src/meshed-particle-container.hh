#ifndef MESHED_PARTICLE_CONTAINER_HH
#define MESHED_PARTICLE_CONTAINER_HH

// for the cases where each particle container has its own mesh.
// Like gone-diego particles

#include "Mesh.hh"

class meshed_particle_container_t {
public:
   meshed_particle_container_t(const Mesh &mesh, const particle_container_t &pc) : mesh(mesh), particle_container(pc) {}
   Mesh mesh;
   particle_container_t particle_container;
   void removed_expired();
};

#endif // MESHED_PARTICLE_CONTAINER_HH
