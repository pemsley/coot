/* coot-utils/read-cube.cc
 *
 * Copyright 2026 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>

#include "read-cube.hh"
#include "utils/coot-utils.hh"

// Bohr to Angstrom
#define BOHR_TO_ANGSTROM 0.52917721067

// Read a Gaussian/ORCA cube file.
//
// Format (all on their own lines):
//   1: comment
//   2: comment (often the property name)
//   3: natoms ox oy oz         (natoms < 0 => an MO cube, with an extra line
//                               after the atoms giving the orbital indices)
//   4: n1 v1x v1y v1z          (n1 < 0 => lengths in Angstrom, else Bohr)
//   5: n2 v2x v2y v2z
//   6: n3 v3x v3y v3z
//   |natoms| lines: Z charge x y z
//   (mo cube only) 1 line: nmo idx1 idx2 ... idxnmo
//   then the volumetric data, axis-2 fastest (then the MO index for an mo cube)
//
coot::util::cube_t
coot::util::read_cube_file(const std::string &file_name, bool read_grid_data) {

   cube_t cube;

   std::ifstream f(file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: read_cube_file(): cannot open " << file_name << std::endl;
      return cube;
   }

   std::string line1, line2;
   std::getline(f, line1);
   std::getline(f, line2);
   cube.title = line1 + " " + line2;

   // line 3: natoms and origin
   int natoms = 0;
   double ox = 0, oy = 0, oz = 0;
   {
      std::string line;
      std::getline(f, line);
      std::istringstream iss(line);
      if (! (iss >> natoms >> ox >> oy >> oz)) {
         std::cout << "WARNING:: read_cube_file(): bad atom-count/origin line in "
                   << file_name << std::endl;
         return cube;
      }
   }
   if (natoms < 0) {
      cube.is_mo_cube = true;
      natoms = -natoms;
   }

   // lines 4-6: voxel counts and axis step vectors. The sign of the first
   // voxel count flags the units for the whole grid (and the atoms).
   int raw_n[3] = {0, 0, 0};
   double vec[3][3];
   for (int i=0; i<3; i++) {
      std::string line;
      std::getline(f, line);
      std::istringstream iss(line);
      if (! (iss >> raw_n[i] >> vec[i][0] >> vec[i][1] >> vec[i][2])) {
         std::cout << "WARNING:: read_cube_file(): bad grid-axis line in "
                   << file_name << std::endl;
         return cube;
      }
   }

   bool in_bohr = (raw_n[0] >= 0); // negative voxel count => already Angstrom
   double to_ang = in_bohr ? BOHR_TO_ANGSTROM : 1.0;

   for (int i=0; i<3; i++)
      cube.n_points[i] = std::abs(raw_n[i]);

   cube.origin = clipper::Coord_orth(ox * to_ang, oy * to_ang, oz * to_ang);
   for (int i=0; i<3; i++)
      cube.axis_step[i] = clipper::Coord_orth(vec[i][0] * to_ang,
                                              vec[i][1] * to_ang,
                                              vec[i][2] * to_ang);

   // atom lines
   cube.atoms.reserve(natoms);
   for (int iat=0; iat<natoms; iat++) {
      std::string line;
      if (! std::getline(f, line)) {
         std::cout << "WARNING:: read_cube_file(): file ended while reading atoms in "
                   << file_name << std::endl;
         return cube;
      }
      std::istringstream iss(line);
      int z = 0;
      double charge = 0, x = 0, y = 0, z_coord = 0;
      if (! (iss >> z >> charge >> x >> y >> z_coord)) {
         std::cout << "WARNING:: read_cube_file(): bad atom line in "
                   << file_name << std::endl;
         return cube;
      }
      clipper::Coord_orth pos(x * to_ang, y * to_ang, z_coord * to_ang);
      cube.atoms.push_back(cube_atom_t(z, charge, pos));
   }

   // mo cube: the orbital-index line
   if (cube.is_mo_cube) {
      std::string line;
      if (! std::getline(f, line)) {
         std::cout << "WARNING:: read_cube_file(): file ended before the MO index line in "
                   << file_name << std::endl;
         return cube;
      }
      std::istringstream iss(line);
      int nmo = 0;
      iss >> nmo;
      for (int i=0; i<nmo; i++) {
         int idx = 0;
         if (iss >> idx)
            cube.mo_indices.push_back(idx);
      }
   }

   if (read_grid_data) {
      std::size_t n_expected =
         static_cast<std::size_t>(cube.n_points[0]) *
         static_cast<std::size_t>(cube.n_points[1]) *
         static_cast<std::size_t>(cube.n_points[2]);
      if (cube.is_mo_cube && ! cube.mo_indices.empty())
         n_expected *= cube.mo_indices.size();
      cube.data.reserve(n_expected);
      float v;
      while (f >> v)
         cube.data.push_back(v);
      if (cube.data.size() != n_expected) {
         std::cout << "WARNING:: read_cube_file(): expected " << n_expected
                   << " grid values but read " << cube.data.size() << " from "
                   << file_name << std::endl;
      }
   }

   cube.read_success = true;
   return cube;
}


mmdb::Manager *
coot::util::cube_t::make_mmdb_manager() const {

   if (atoms.empty())
      return nullptr;

   // atomic number -> element symbol
   std::map<int, std::string> z_to_symbol;
   {
      std::vector<std::pair<std::string, int> > l = coot::util::atomic_number_atom_list();
      for (const auto &p : l)
         z_to_symbol[p.second] = p.first;
   }

   mmdb::Manager *mol = new mmdb::Manager;
   mmdb::Model *model_p = new mmdb::Model;
   mmdb::Chain *chain_p = new mmdb::Chain;
   mmdb::Residue *residue_p = new mmdb::Residue;
   chain_p->SetChainID("A");
   residue_p->SetResName("MOL");
   residue_p->seqNum = 1;

   std::map<int, int> element_counts; // per-element serial, for atom names

   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      const cube_atom_t &ca = atoms[iat];
      std::string symbol = "C"; // fallback
      std::map<int, std::string>::const_iterator it = z_to_symbol.find(ca.atomic_number);
      if (it != z_to_symbol.end())
         symbol = it->second;

      int count = ++element_counts[ca.atomic_number];
      std::string atom_name = symbol + coot::util::int_to_string(count);

      // right-justified 2-char element name (e.g. " C", "NA")
      std::string ele = symbol;
      if (ele.length() == 1)
         ele = " " + ele;

      mmdb::Atom *at = new mmdb::Atom;
      at->SetCoordinates(ca.position.x(), ca.position.y(), ca.position.z(), 1.0, 20.0);
      at->SetAtomName(atom_name.c_str());
      at->SetElementName(ele.c_str());
      at->Het = 1;
      residue_p->AddAtom(at);
   }

   chain_p->AddResidue(residue_p);
   model_p->AddChain(chain_p);
   mol->AddModel(model_p);
   mol->FinishStructEdit();
   return mol;
}


bool
coot::util::cube_t::data_has_negative_lobe() const {

   if (data.empty()) return false;
   float dmin = *std::min_element(data.begin(), data.end());
   float dmax = *std::max_element(data.begin(), data.end());
   // a real negative lobe (not just numerical noise near a positive density)
   return (dmin < -0.01f * std::fabs(dmax));
}


bool
coot::util::cube_t::make_xmap(clipper::Xmap<float> *xmap_p) const {

   if (data.empty()) {
      std::cout << "WARNING:: make_xmap(): no grid data" << std::endl;
      return false;
   }
   if (n_points[0] <= 0 || n_points[1] <= 0 || n_points[2] <= 0) {
      std::cout << "WARNING:: make_xmap(): bad grid dimensions" << std::endl;
      return false;
   }

   // This path handles the usual ORCA/Gaussian cube: an axis-aligned,
   // orthogonal grid. Check that the step vectors point along x, y, z.
   double dx = axis_step[0].x();
   double dy = axis_step[1].y();
   double dz = axis_step[2].z();
   auto is_off_axis = [] (double off, double along) {
      return std::fabs(off) > 1e-4 * (std::fabs(along) + 1e-9);
   };
   bool not_aligned =
      is_off_axis(axis_step[0].y(), dx) || is_off_axis(axis_step[0].z(), dx) ||
      is_off_axis(axis_step[1].x(), dy) || is_off_axis(axis_step[1].z(), dy) ||
      is_off_axis(axis_step[2].x(), dz) || is_off_axis(axis_step[2].y(), dz);
   if (not_aligned) {
      std::cout << "WARNING:: make_xmap(): cube grid is not axis-aligned/orthogonal - "
                << "cannot build an Xmap (this would need an NXmap)" << std::endl;
      return false;
   }
   if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
      std::cout << "WARNING:: make_xmap(): non-positive grid spacing" << std::endl;
      return false;
   }

   int nu = n_points[0];
   int nv = n_points[1];
   int nw = n_points[2];

   double a = nu * dx;
   double b = nv * dy;
   double c = nw * dz;

   clipper::Cell cell(clipper::Cell_descr(a, b, c, 90.0, 90.0, 90.0));
   clipper::Spacegroup spgr(clipper::Spgr_descr(1)); // P1
   clipper::Grid_sampling gs(nu, nv, nw);
   xmap_p->init(spgr, cell, gs);

   // The cube origin is (usually) not an integer number of voxels from the
   // cell origin, so round it. Any residual is a sub-voxel (< half a voxel)
   // registration error - invisible for orbital display.
   int nx_start = static_cast<int>(std::lround(origin.x() / dx));
   int ny_start = static_cast<int>(std::lround(origin.y() / dy));
   int nz_start = static_cast<int>(std::lround(origin.z() / dz));

   // cube data order: axis-2 (w) fastest, then v, then u; for an mo cube the
   // orbital index is innermost - we use the first orbital only.
   int n_mo = (is_mo_cube && ! mo_indices.empty()) ? static_cast<int>(mo_indices.size()) : 1;

   clipper::Xmap<float>::Map_reference_coord mrc(*xmap_p);
   std::size_t idx = 0;
   for (int i=0; i<nu; i++) {
      for (int j=0; j<nv; j++) {
         for (int k=0; k<nw; k++) {
            float val = data[idx];
            idx += n_mo; // step over the other orbitals at this voxel
            mrc.set_coord(clipper::Coord_grid(i + nx_start, j + ny_start, k + nz_start));
            (*xmap_p)[mrc] = val;
         }
      }
   }
   return true;
}
