/*
 * src/meshed-particle-container.hh
 *
 * Copyright 2023 by Medical Research Council
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
