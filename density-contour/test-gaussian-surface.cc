/*
 * density-contour/test-gaussian-surface.cc
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */

#include <iostream>
#include <string>
#include "gaussian-surface.hh"

int main(int argc, char **argv) {

   if (argc > 1) {
      std::string file_name = argv[1];

      mmdb::Manager *mol = new mmdb::Manager;
      mol->ReadCoorFile(file_name.c_str());

      coot::gaussian_surface_t gauss_surf(mol, "A");
      coot::simple_mesh_t mesh = gauss_surf.get_surface();

      std::cout << "test-gaussian-surface got " << mesh.vertices.size() << " vertices and " << mesh.triangles.size() << std::endl;
   }
   return 0;
}
