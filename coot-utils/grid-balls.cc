
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <tuple>
#include <functional>

#include "grid-balls.hh"

std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::grid_balls_t::get_extents(mmdb::Manager *mol) const {

   float mol_x_min =  1.0e30;
   float mol_y_min =  1.0e30;
   float mol_z_min =  1.0e30;
   float mol_x_max = -1.0e30;
   float mol_y_max = -1.0e30;
   float mol_z_max = -1.0e30;

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        if (at->x < mol_x_min) mol_x_min = at->x;
                        if (at->y < mol_y_min) mol_y_min = at->y;
                        if (at->z < mol_z_min) mol_z_min = at->z;
                        if (at->x > mol_x_max) mol_x_max = at->x;
                        if (at->y > mol_y_max) mol_y_max = at->y;
                        if (at->z > mol_z_max) mol_z_max = at->z;
                     }
                  }
               }
            }
         }
      }
   }

   // test here for sane min and max values

   return std::make_pair(clipper::Coord_orth(mol_x_min, mol_y_min, mol_z_min),
                         clipper::Coord_orth(mol_x_max, mol_y_max, mol_z_max));

}

int
coot::grid_balls_t::grid_index(const triple_index_t &t) const {
   return t.ix + t.iy * nx + t.iz * nx * ny;
}

coot::grid_balls_t::triple_index_t
coot::grid_balls_t::deindex(int i) const {
   triple_index_t t;
   t.ix = i % nx;
   t.iy = (i / nx) % ny;
   t.iz = i / (nx * ny);
   return t;
}

coot::grid_balls_t::grid_balls_t(int imol_in, mmdb::Manager *mol, const coot::protein_geometry *geom_p_in,
                                 float small_ball_radius, float big_ball_radius) {

   imol = imol_in;
   this->mol = mol;
   geom_p = geom_p_in;
   probe_in_radius  = small_ball_radius;
   probe_out_radius = big_ball_radius;

   n_grids_per_angstrom = 1.0; //testing

   // get extents of molecule
   std::pair<clipper::Coord_orth, clipper::Coord_orth> ee = get_extents(mol);

   // calculate grid size
   float mol_x_min = ee.first.x();
   float mol_y_min = ee.first.y();
   float mol_z_min = ee.first.z();
   float mol_x_max = ee.second.x();
   float mol_y_max = ee.second.y();
   float mol_z_max = ee.second.z();

   float extra_extents = 5.0; // angstroms

   grid_min_x = mol_x_min - extra_extents;
   grid_min_y = mol_y_min - extra_extents;
   grid_min_z = mol_z_min - extra_extents;
   grid_max_x = mol_x_max + extra_extents;
   grid_max_y = mol_y_max + extra_extents;
   grid_max_z = mol_z_max + extra_extents;

   nx = int((grid_max_x - grid_min_x) * n_grids_per_angstrom) + 1;
   ny = int((grid_max_y - grid_min_y) * n_grids_per_angstrom) + 1;
   nz = int((grid_max_z - grid_min_z) * n_grids_per_angstrom) + 1;

   grid.resize(nx * ny * nz);

   test_grid();
   test_index_deindex();

   brick_the_model(mol);
   find_cavities();

}

coot::grid_balls_t::triple_index_t
coot::grid_balls_t::mol_space_to_grid_point(const point_3d_t &p) const {

   triple_index_t t;
   float x = p.x - grid_min_x;
   float y = p.y - grid_min_y;
   float z = p.z - grid_min_z;

   t.ix = static_cast<int>(std::round(x * n_grids_per_angstrom));
   t.iy = static_cast<int>(std::round(y * n_grids_per_angstrom));
   t.iz = static_cast<int>(std::round(z * n_grids_per_angstrom));

   if (false) {
      // std::cout << "grid_min_x " << grid_min_x << " grid_min_y " << grid_min_y << " grid_min_z " << grid_min_z << std::endl;
      std::cout << "p.x " << p.x << " p.y " << p.y << " p.z " << p.z << std::endl;
      std::cout << "x " << x << " y " << y << " z " << z << std::endl;
      std::cout << "t.ix " << t.ix << " t.iy " << t.iy << " t.iz " << t.iz << std::endl;
   }

   return t;

}

coot::grid_balls_t::point_3d_t
coot::grid_balls_t::grid_point_to_mol_space(const triple_index_t &t) const {

   point_3d_t p;
   p.x = float(t.ix+0) / n_grids_per_angstrom + grid_min_x;
   p.y = float(t.iy+0) / n_grids_per_angstrom + grid_min_y;
   p.z = float(t.iz+0) / n_grids_per_angstrom + grid_min_z;

   return p;
}


void
coot::grid_balls_t::test_grid() const {

   std::cout << "testing grid to space..." << std::endl;
   int n_correct = 0;
   int n_wrong = 0;

   for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
         for (int iz=0; iz<nz; iz++) {
            point_3d_t p = grid_point_to_mol_space(triple_index_t(ix, iy, iz));
            triple_index_t t = mol_space_to_grid_point(p);
            int tix_as_int =  int(round(t.ix));
            int tiy_as_int =  int(round(t.iy));
            int tiz_as_int =  int(round(t.iz));
            if (ix != tix_as_int) {
               std::cout << "Error in grid indexing X: input: " << ix << " " << iy << " " << iz
                         << " as_int: " << tix_as_int << " " << tiy_as_int << " " << tiz_as_int
                         << " result: " << t.ix << " " << t.iy << " " << t.iz
                         << "\n";
               n_wrong++;
            }
            if (iy != tiy_as_int) {
               std::cout << "Error in grid indexing Y: input: " << ix << " " << iy << " " << iz
                         << " as_int: " << tix_as_int << " " << tiy_as_int << " " << tiz_as_int
                         << " result: " << t.ix << " " << t.iy << " " << t.iz
                         << "\n";
               n_wrong++;
            }
            if (iz != tiz_as_int) {
               std::cout << "Error in grid indexing Z: input " << ix << " " << iy << " " << iz
                         << " as_int: " << tix_as_int << " " << tiy_as_int << " " << tiz_as_int
                         << " result: " << t.ix << " " << t.iy << " " << t.iz
                         << "\n";
               n_wrong++;
            }


            if (t.ix != ix || t.iy != iy || t.iz != iz) {
            } else {
               n_correct++;
            }
         }
      }
   }
   int n_total = n_correct + n_wrong;
   std::cout << "testing done. n_correct: " << n_correct << " n_wrong " << n_wrong
             << "  " << 100.0 * static_cast<float>(n_wrong)/static_cast<float>(n_total) << " %" << std::endl;

}


void
coot::grid_balls_t::test_index_deindex() const {

   std::cout << "testing index/deindex..." << std::endl;
   int n_correct = 0;
   int n_wrong = 0;

   for (int i=0; i<nx*ny*nz; i++) {
      triple_index_t t = deindex(i);
      int i_as_int = grid_index(t);
      if (i != i_as_int) {
         std::cout << "Error in index/deindex: input: " << i << " as_int: " << i_as_int
                   << " result: " << t.ix << " " << t.iy << " " << t.iz
                   << "\n";
         n_wrong++;
      } else {
         n_correct++;
      }
   }
   int n_total = n_correct + n_wrong;
   std::cout << "testing for index/deindex done: n_correct: " << n_correct << " n_wrong " << n_wrong
             << "  " << 100.0 * static_cast<float>(n_wrong)/static_cast<float>(n_total) << " %" << std::endl;

}


namespace {
   // Fallback when the dictionary has no VdW radius for this atom
   // (get_vdw_radius() returns -1.1). Bondi radii - matching the intent of
   // the defaults in coot-utils/hole.cc (but without its " 0"/oxygen typo).
   float element_to_vdw_radius(const std::string &element) {
      std::string e = element;
      // mmdb elements are right-justified in a 2-char field, e.g. " C"
      std::string::size_type f = e.find_first_not_of(' ');
      if (f != std::string::npos) e = e.substr(f);
      if (e == "H")  return 1.20f;
      if (e == "C")  return 1.70f;
      if (e == "N")  return 1.55f;
      if (e == "O")  return 1.52f;
      if (e == "S")  return 1.80f;
      if (e == "P")  return 1.80f;
      return 1.70f; // reasonable default
   }

   // The VdW radius for an atom: from the dictionary if available, else the
   // element-based fallback above.
   float atom_vdw_radius(const coot::protein_geometry *geom_p, int imol,
                         mmdb::Atom *at, const std::string &residue_name) {
      double radius = -1.1;
      if (geom_p)
         radius = geom_p->get_vdw_radius(std::string(at->name), residue_name,
                                         imol, false); // heavy-atom VdW, as in KVFinder
      if (radius <= 0.0)
         radius = element_to_vdw_radius(at->element);
      return static_cast<float>(radius);
   }
}

// Mark every grid point that lies within an atom's VdW radius as occupied,
// and record which atom(s) cover it. This is the KVFinder "occupancy grid"
// step: the empty (unmarked) points are what the probes later roll through.
void
coot::grid_balls_t::brick_the_model(mmdb::Manager *mol) {

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  std::string residue_name = residue_p->GetResName();
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (at->isTer()) continue;

                     float radius = atom_vdw_radius(geom_p, imol, at, residue_name);

                     // rasterise the VdW sphere into the grid
                     point_3d_t atom_pos(at->x, at->y, at->z);
                     triple_index_t centre = mol_space_to_grid_point(atom_pos);
                     int ir = static_cast<int>(std::ceil(radius * n_grids_per_angstrom));
                     float radius_sqrd = radius * radius;

                     for (int dz=-ir; dz<=ir; dz++) {
                        int gz = centre.iz + dz;
                        if (gz < 0 || gz >= nz) continue;
                        for (int dy=-ir; dy<=ir; dy++) {
                           int gy = centre.iy + dy;
                           if (gy < 0 || gy >= ny) continue;
                           for (int dx=-ir; dx<=ir; dx++) {
                              int gx = centre.ix + dx;
                              if (gx < 0 || gx >= nx) continue;

                              triple_index_t t(gx, gy, gz);
                              point_3d_t gp = grid_point_to_mol_space(t);
                              float ddx = gp.x - atom_pos.x;
                              float ddy = gp.y - atom_pos.y;
                              float ddz = gp.z - atom_pos.z;
                              float d_sqrd = ddx*ddx + ddy*ddy + ddz*ddz;
                              if (d_sqrd <= radius_sqrd) {
                                 int idx = grid_index(t);
                                 grid[idx].v = 1.0f;         // occupied
                                 grid[idx].atoms.push_back(at);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

}


// Rasterise every atom, inflated to (VdW radius + probe_radius), into a grid.
// A "blocked" point is one where a probe of this radius, centred there, would
// clash with the biomolecule - so the *free* (unblocked) points are the valid
// probe-centre positions.
std::vector<unsigned char>
coot::grid_balls_t::make_blocked_grid(float probe_radius) const {

   std::vector<unsigned char> blocked(nx * ny * nz, 0);
   if (! mol) return blocked;

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (! model_p) continue;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (! residue_p) continue;
            std::string residue_name = residue_p->GetResName();
            int n_atoms = residue_p->GetNumberOfAtoms();
            for (int iat=0; iat<n_atoms; iat++) {
               mmdb::Atom *at = residue_p->GetAtom(iat);
               if (at->isTer()) continue;

               float radius = atom_vdw_radius(geom_p, imol, at, residue_name) + probe_radius;

               point_3d_t atom_pos(at->x, at->y, at->z);
               triple_index_t centre = mol_space_to_grid_point(atom_pos);
               int ir = static_cast<int>(std::ceil(radius * n_grids_per_angstrom));
               float radius_sqrd = radius * radius;

               for (int dz=-ir; dz<=ir; dz++) {
                  int gz = centre.iz + dz;
                  if (gz < 0 || gz >= nz) continue;
                  for (int dy=-ir; dy<=ir; dy++) {
                     int gy = centre.iy + dy;
                     if (gy < 0 || gy >= ny) continue;
                     for (int dx=-ir; dx<=ir; dx++) {
                        int gx = centre.ix + dx;
                        if (gx < 0 || gx >= nx) continue;
                        triple_index_t t(gx, gy, gz);
                        point_3d_t gp = grid_point_to_mol_space(t);
                        float ddx = gp.x - atom_pos.x;
                        float ddy = gp.y - atom_pos.y;
                        float ddz = gp.z - atom_pos.z;
                        if (ddx*ddx + ddy*ddy + ddz*ddz <= radius_sqrd)
                           blocked[grid_index(t)] = 1;
                     }
                  }
               }
            }
         }
      }
   }
   return blocked;
}

// Flood-fill (DFS) the free (unblocked) probe-centre positions that are
// connected to the bulk exterior, seeded from a grid corner - which sits in the
// padding, so it is free unless the probe is large enough to have swallowed it
// (in which case we fall back to the first free cell).
std::vector<unsigned char>
coot::grid_balls_t::flood_fill_exterior(const std::vector<unsigned char> &blocked) const {

   std::vector<unsigned char> reachable(nx * ny * nz, 0);

   int seed = -1;
   int corner = grid_index(triple_index_t(0, 0, 0));
   if (! blocked[corner]) {
      seed = corner;
   } else {
      for (int i=0; i<nx*ny*nz; i++) {
         if (! blocked[i]) { seed = i; break; }
      }
   }
   if (seed < 0) return reachable; // wholly blocked - should not happen

   const int dirs[6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
   std::vector<int> stack;
   stack.push_back(seed);
   reachable[seed] = 1;
   while (! stack.empty()) {
      int idx = stack.back();
      stack.pop_back();
      triple_index_t t = deindex(idx);
      for (int d=0; d<6; d++) {
         int gx = t.ix + dirs[d][0];
         int gy = t.iy + dirs[d][1];
         int gz = t.iz + dirs[d][2];
         if (gx < 0 || gx >= nx) continue;
         if (gy < 0 || gy >= ny) continue;
         if (gz < 0 || gz >= nz) continue;
         int nidx = grid_index(triple_index_t(gx, gy, gz));
         if (! blocked[nidx] && ! reachable[nidx]) {
            reachable[nidx] = 1;
            stack.push_back(nidx);
         }
      }
   }
   return reachable;
}

// A grid point is "accessible" to a probe of this radius if an exterior-connected
// free probe-centre lies within probe_radius of it - i.e. the probe surface,
// rolled in from the bulk, reaches it.
std::vector<unsigned char>
coot::grid_balls_t::make_accessible_grid(float probe_radius) const {

   std::vector<unsigned char> blocked   = make_blocked_grid(probe_radius);
   std::vector<unsigned char> reachable = flood_fill_exterior(blocked);
   std::vector<unsigned char> accessible(nx * ny * nz, 0);

   int ir = static_cast<int>(std::ceil(probe_radius * n_grids_per_angstrom));
   float thr = probe_radius * n_grids_per_angstrom; // in grid units
   float thr_sqrd = thr * thr;

   for (int idx=0; idx<nx*ny*nz; idx++) {
      if (reachable[idx]) { accessible[idx] = 1; continue; } // a free centre itself
      triple_index_t t = deindex(idx);
      bool found = false;
      for (int dz=-ir; dz<=ir && !found; dz++) {
         int gz = t.iz + dz;
         if (gz < 0 || gz >= nz) continue;
         for (int dy=-ir; dy<=ir && !found; dy++) {
            int gy = t.iy + dy;
            if (gy < 0 || gy >= ny) continue;
            for (int dx=-ir; dx<=ir; dx++) {
               int gx = t.ix + dx;
               if (gx < 0 || gx >= nx) continue;
               if (dx*dx + dy*dy + dz*dz > thr_sqrd) continue;
               if (reachable[grid_index(triple_index_t(gx, gy, gz))]) { found = true; break; }
            }
         }
      }
      if (found) accessible[idx] = 1;
   }
   return accessible;
}

// KVFinder cavity extraction: a point is a cavity point if the small (Probe In)
// surface reaches it but the big (Probe Out) surface does not.
void
coot::grid_balls_t::find_cavities() {

   std::vector<unsigned char> accessible_in  = make_accessible_grid(probe_in_radius);
   std::vector<unsigned char> accessible_out = make_accessible_grid(probe_out_radius);

   cavity_grid.assign(nx * ny * nz, 0);
   unsigned int n_cavity_points = 0;
   for (int i=0; i<nx*ny*nz; i++) {
      if (accessible_in[i] && ! accessible_out[i]) {
         cavity_grid[i] = 1;
         n_cavity_points++;
      }
   }
   std::cout << "find_cavities(): " << n_cavity_points << " cavity grid points of "
             << (nx*ny*nz) << std::endl;

   cluster_cavities();
   compute_lining_residues();
   float min_volume_to_segment = 100.0;
   float min_prominence = 0.9;  // keep necks deeper than this (spectrum tops ~1.0 A)
   int min_seed_depth = 1;      // cavities are thin shells (max depth ~2.7 A)
   float merge_dist = 4.0;
   segment_cavities(min_volume_to_segment, min_seed_depth, merge_dist, min_prominence);
}

// Cluster the raw cavity points into distinct cavities (connected components,
// 6-connected DFS), compute each cavity's volume, drop those below volume_min,
// and rank the survivors by volume - largest first, as KVFinder does.
void
coot::grid_balls_t::cluster_cavities(float volume_min) {

   cavities.clear();

   const int n = nx * ny * nz;
   std::vector<int> labels(n, 0); // 0 = unlabelled
   const int dirs[6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
   float voxel_volume = 1.0f / (n_grids_per_angstrom * n_grids_per_angstrom * n_grids_per_angstrom);

   int n_labels = 0;
   for (int start=0; start<n; start++) {
      if (! cavity_grid[start]) continue;
      if (labels[start]) continue;

      n_labels++;
      cavity_t cav;
      cav.id = n_labels;

      std::vector<int> stack;
      stack.push_back(start);
      labels[start] = n_labels;
      while (! stack.empty()) {
         int idx = stack.back();
         stack.pop_back();
         cav.grid_indices.push_back(idx);
         triple_index_t t = deindex(idx);
         for (int d=0; d<6; d++) {
            int gx = t.ix + dirs[d][0];
            int gy = t.iy + dirs[d][1];
            int gz = t.iz + dirs[d][2];
            if (gx < 0 || gx >= nx) continue;
            if (gy < 0 || gy >= ny) continue;
            if (gz < 0 || gz >= nz) continue;
            int nidx = grid_index(triple_index_t(gx, gy, gz));
            if (cavity_grid[nidx] && ! labels[nidx]) {
               labels[nidx] = n_labels;
               stack.push_back(nidx);
            }
         }
      }

      cav.volume = static_cast<float>(cav.grid_indices.size()) * voxel_volume;

      // centroid (mol space) and surface area (exposed voxel faces). A face is
      // exposed when its neighbour is not a cavity point (or off-grid). Distinct
      // cavities can never be face-adjacent (they'd be one component), so testing
      // cavity_grid is enough.
      float face_area = 1.0f / (n_grids_per_angstrom * n_grids_per_angstrom);
      unsigned int n_exposed_faces = 0;
      double sx = 0.0, sy = 0.0, sz = 0.0;
      for (unsigned int i=0; i<cav.grid_indices.size(); i++) {
         triple_index_t t = deindex(cav.grid_indices[i]);
         point_3d_t p = grid_point_to_mol_space(t);
         sx += p.x; sy += p.y; sz += p.z;
         for (int d=0; d<6; d++) {
            int gx = t.ix + dirs[d][0];
            int gy = t.iy + dirs[d][1];
            int gz = t.iz + dirs[d][2];
            if (gx < 0 || gx >= nx || gy < 0 || gy >= ny || gz < 0 || gz >= nz) {
               n_exposed_faces++; // grid edge - treat as exposed
            } else {
               if (! cavity_grid[grid_index(triple_index_t(gx, gy, gz))])
                  n_exposed_faces++;
            }
         }
      }
      cav.surface_area = static_cast<float>(n_exposed_faces) * face_area;
      double np = static_cast<double>(cav.grid_indices.size());
      if (np > 0.0)
         cav.centre = point_3d_t(static_cast<float>(sx/np),
                                 static_cast<float>(sy/np),
                                 static_cast<float>(sz/np));

      cavities.push_back(cav);
   }

   // volume filter
   if (volume_min > 0.0f) {
      std::vector<cavity_t> kept;
      for (unsigned int i=0; i<cavities.size(); i++)
         if (cavities[i].volume >= volume_min)
            kept.push_back(cavities[i]);
      cavities = kept;
   }

   // rank by volume, largest first, then renumber
   std::sort(cavities.begin(), cavities.end(),
             [](const cavity_t &a, const cavity_t &b) { return a.volume > b.volume; });
   for (unsigned int i=0; i<cavities.size(); i++)
      cavities[i].id = i + 1;

   std::cout << "cluster_cavities(): " << cavities.size() << " cavities" << std::endl;
   for (unsigned int i=0; i<cavities.size(); i++) {
      const cavity_t &c = cavities[i];
      std::cout << "   cavity " << c.id << ": " << c.grid_indices.size()
                << " points, volume " << c.volume << " A^3, area "
                << c.surface_area << " A^2, centre ("
                << c.centre.x << " " << c.centre.y << " " << c.centre.z << ")"
                << std::endl;
   }
}

// Write a plain table of cavity point positions: one "cavity_id x y z" line per
// grid point, preceded by a per-cavity comment header (points/volume/area/centre).
// A script can read this and add each cavity's points to a generic display object
// (one colour per cavity) via to_generic_object_add_points().
void
coot::grid_balls_t::write_cavity_points(const std::string &file_name) const {

   std::ofstream f(file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: write_cavity_points() cannot open " << file_name << std::endl;
      return;
   }

   f << std::fixed << std::setprecision(3);
   f << "# grid-balls cavity points\n";
   f << "# " << cavities.size() << " cavities\n";
   f << "# columns: cavity_id x y z\n";

   int n_points = 0;
   for (unsigned int ic=0; ic<cavities.size(); ic++) {
      const cavity_t &c = cavities[ic];
      f << "# cavity " << c.id
        << "  points " << c.grid_indices.size()
        << "  volume " << c.volume
        << "  area "   << c.surface_area
        << "  centre " << c.centre.x << " " << c.centre.y << " " << c.centre.z << "\n";
      for (unsigned int i=0; i<c.grid_indices.size(); i++) {
         point_3d_t p = grid_point_to_mol_space(deindex(c.grid_indices[i]));
         f << c.id << " " << p.x << " " << p.y << " " << p.z << "\n";
         n_points++;
      }
   }
   f.close();
   std::cout << "write_cavity_points(): wrote " << n_points << " points for "
             << cavities.size() << " cavities to " << file_name << std::endl;
}

// Fill cavity_t::lining_residues for each cavity in cavs: any residue with an atom
// within contact_dist of one of the cavity's grid points. The cavities are used as a
// spatial index (a label grid keyed by position in cavs), then the atoms are walked
// once. Works for top-level cavities and for (spatially disjoint) subpockets alike.
void
coot::grid_balls_t::compute_lining_residues_for(std::vector<cavity_t> &cavs, float contact_dist) {

   if (! mol) return;
   if (cavs.empty()) return;

   const int n = nx * ny * nz;

   // label grid: 0 = none, else (position-in-cavs + 1)
   std::vector<int> label(n, 0);
   for (unsigned int ic=0; ic<cavs.size(); ic++)
      for (unsigned int i=0; i<cavs[ic].grid_indices.size(); i++)
         label[cavs[ic].grid_indices[i]] = ic + 1;

   std::vector<std::set<mmdb::Residue *> > res_sets(cavs.size());

   int   R   = static_cast<int>(std::ceil(contact_dist * n_grids_per_angstrom));
   float cd2 = contact_dist * contact_dist;

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (! model_p) continue;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (! residue_p) continue;
            int n_atoms = residue_p->GetNumberOfAtoms();
            for (int iat=0; iat<n_atoms; iat++) {
               mmdb::Atom *at = residue_p->GetAtom(iat);
               if (at->isTer()) continue;

               point_3d_t ap(at->x, at->y, at->z);
               triple_index_t centre = mol_space_to_grid_point(ap);

               std::set<int> hit; // labels this atom is in contact with
               for (int dz=-R; dz<=R; dz++) {
                  int gz = centre.iz + dz;
                  if (gz < 0 || gz >= nz) continue;
                  for (int dy=-R; dy<=R; dy++) {
                     int gy = centre.iy + dy;
                     if (gy < 0 || gy >= ny) continue;
                     for (int dx=-R; dx<=R; dx++) {
                        int gx = centre.ix + dx;
                        if (gx < 0 || gx >= nx) continue;
                        triple_index_t t(gx, gy, gz);
                        int L = label[grid_index(t)];
                        if (L == 0) continue;
                        if (hit.find(L) != hit.end()) continue; // already have it
                        point_3d_t gp = grid_point_to_mol_space(t);
                        float ddx = gp.x - ap.x;
                        float ddy = gp.y - ap.y;
                        float ddz = gp.z - ap.z;
                        if (ddx*ddx + ddy*ddy + ddz*ddz <= cd2)
                           hit.insert(L);
                     }
                  }
               }
               for (std::set<int>::const_iterator it=hit.begin(); it!=hit.end(); ++it)
                  res_sets[*it - 1].insert(residue_p);
            }
         }
      }
   }

   for (unsigned int ic=0; ic<cavs.size(); ic++) {
      const std::set<mmdb::Residue *> &s = res_sets[ic];
      cavs[ic].lining_residues.assign(s.begin(), s.end());
      std::sort(cavs[ic].lining_residues.begin(), cavs[ic].lining_residues.end(),
                [](mmdb::Residue *a, mmdb::Residue *b) {
                   int cc = std::string(a->GetChainID()).compare(b->GetChainID());
                   if (cc != 0) return cc < 0;
                   return a->GetSeqNum() < b->GetSeqNum();
                });
   }
}

namespace {
   void print_lining_residues(const coot::grid_balls_t::cavity_t &c, const std::string &tag) {
      std::cout << "   " << tag << " (" << c.lining_residues.size() << " lining residues):";
      for (unsigned int i=0; i<c.lining_residues.size(); i++) {
         mmdb::Residue *r = c.lining_residues[i];
         std::cout << " " << r->GetChainID() << r->GetSeqNum() << r->GetInsCode()
                   << "(" << r->GetResName() << ")";
      }
      std::cout << std::endl;
   }
}

void
coot::grid_balls_t::compute_lining_residues(float contact_dist) {

   compute_lining_residues_for(cavities, contact_dist);

   std::cout << "compute_lining_residues(): (contact " << contact_dist << " A)" << std::endl;
   for (unsigned int ic=0; ic<cavities.size(); ic++)
      if (! cavities[ic].lining_residues.empty())
         print_lining_residues(cavities[ic], "cavity " + std::to_string(cavities[ic].id));
}

// Split one cavity into subpockets by distance-transform watershed:
//  1. D = distance-to-boundary of each cavity voxel (multi-source BFS from the wall)
//  2. seeds = deep local maxima of D (>= min_seed_depth), merged if within merge_dist
//  3. priority-flood watershed from the seeds (highest D first) - the basin boundary
//     falls at the saddle (the neck between chambers).
// Returns the subpockets (with volume/area/centre filled, parent_id = c.id), or an
// empty vector if the cavity has a single chamber (does not split).
std::vector<coot::grid_balls_t::cavity_t>
coot::grid_balls_t::segment_one_cavity(const cavity_t &c, int min_seed_depth,
                                       float merge_dist, float min_prominence) const {

   const bool debug = false; // set true for per-cavity depth/seed/merge diagnostics
   std::vector<cavity_t> subs;
   const int n = nx * ny * nz;
   const int d6[6][3] = { {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1} };

   // cavity membership
   std::vector<char> in_cav(n, 0);
   for (unsigned int i=0; i<c.grid_indices.size(); i++) in_cav[c.grid_indices[i]] = 1;

   // distance-to-boundary D, via a chamfer (3-4-5) distance transform computed by
   // Dijkstra from the cavity wall. D is ~CH_SCALE x the true distance in grid cells,
   // so a wide chamber centre is genuinely deeper than a thin neck (unlike the coarse
   // city-block BFS this replaced) - which spreads out the merge prominences so that
   // min_prominence becomes a usable dial. (The depth/prominence thresholds below
   // carry the CH_SCALE factor to stay in Angstroms.)
   const int CH_SCALE = 3; // chamfer face-step weight; D units = CH_SCALE per grid cell
   const int D_INF = 1 << 29;
   std::vector<int> D(n, 0);
   for (unsigned int i=0; i<c.grid_indices.size(); i++) D[c.grid_indices[i]] = D_INF;

   // 26 neighbour offsets with chamfer weights (face 3, edge 4, corner 5)
   int off26[26][3];
   int wt26[26];
   {
      int m26 = 0;
      for (int dz=-1; dz<=1; dz++)
         for (int dy=-1; dy<=1; dy++)
            for (int dx=-1; dx<=1; dx++) {
               if (dx==0 && dy==0 && dz==0) continue;
               off26[m26][0]=dx; off26[m26][1]=dy; off26[m26][2]=dz;
               int ax=(dx<0?-dx:dx), ay=(dy<0?-dy:dy), az=(dz<0?-dz:dz);
               int l1 = ax+ay+az;
               wt26[m26] = (l1==1) ? 3 : (l1==2 ? 4 : 5);
               m26++;
            }
   }

   // seed: cavity cells touching the wall (a non-cavity neighbour or the grid edge)
   // start at their smallest wall-step distance; then relax inward (Dijkstra).
   std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int> >,
                       std::greater<std::pair<int,int> > > pqd; // (dist, idx)
   for (unsigned int i=0; i<c.grid_indices.size(); i++) {
      int idx = c.grid_indices[i];
      triple_index_t t = deindex(idx);
      int best = D_INF;
      for (int kk=0; kk<26; kk++) {
         int gx=t.ix+off26[kk][0], gy=t.iy+off26[kk][1], gz=t.iz+off26[kk][2];
         bool wall = (gx<0||gx>=nx||gy<0||gy>=ny||gz<0||gz>=nz);
         if (! wall) wall = ! in_cav[grid_index(triple_index_t(gx,gy,gz))];
         if (wall && wt26[kk] < best) best = wt26[kk];
      }
      if (best < D_INF) { D[idx] = best; pqd.push(std::make_pair(best, idx)); }
   }
   while (! pqd.empty()) {
      std::pair<int,int> topp = pqd.top(); pqd.pop();
      int dcur = topp.first, idx = topp.second;
      if (dcur > D[idx]) continue;
      triple_index_t t = deindex(idx);
      for (int kk=0; kk<26; kk++) {
         int gx=t.ix+off26[kk][0], gy=t.iy+off26[kk][1], gz=t.iz+off26[kk][2];
         if (gx<0||gx>=nx||gy<0||gy>=ny||gz<0||gz>=nz) continue;
         int nidx = grid_index(triple_index_t(gx,gy,gz));
         if (! in_cav[nidx]) continue;
         int nd = dcur + wt26[kk];
         if (nd < D[nidx]) { D[nidx] = nd; pqd.push(std::make_pair(nd, nidx)); }
      }
   }

   if (debug) { // dynamic range of the distance field (thin cavities have small max)
      int d_max = 0;
      for (unsigned int i=0; i<c.grid_indices.size(); i++)
         if (D[c.grid_indices[i]] > d_max) d_max = D[c.grid_indices[i]];
      std::cout << "   [segment cavity " << c.id << "] " << c.grid_indices.size()
                << " cells, max depth " << (float(d_max)/(n_grids_per_angstrom*CH_SCALE))
                << " A (min_seed_depth " << min_seed_depth << " A)" << std::endl;
   }

   // seed candidates: 26-neighbour local maxima of D at least min_seed_depth (A) deep.
   // D is in chamfer units, so scale the threshold by CH_SCALE*gpa.
   int min_seed_units = static_cast<int>(min_seed_depth * n_grids_per_angstrom * CH_SCALE + 0.5f);
   std::vector<int> cand;
   for (unsigned int i=0; i<c.grid_indices.size(); i++) {
      int idx = c.grid_indices[i];
      int d0 = D[idx];
      if (d0 < min_seed_units) continue;
      triple_index_t t = deindex(idx);
      bool is_max = true;
      for (int dz=-1; dz<=1 && is_max; dz++)
         for (int dy=-1; dy<=1 && is_max; dy++)
            for (int dx=-1; dx<=1; dx++) {
               if (dx==0 && dy==0 && dz==0) continue;
               int gx=t.ix+dx, gy=t.iy+dy, gz=t.iz+dz;
               if (gx<0||gx>=nx||gy<0||gy>=ny||gz<0||gz>=nz) continue;
               int nidx = grid_index(triple_index_t(gx,gy,gz));
               if (in_cav[nidx] && D[nidx] > d0) { is_max = false; break; }
            }
      if (is_max) cand.push_back(idx);
   }

   // deepest-first, then greedily keep candidates that are >= merge_dist apart
   std::sort(cand.begin(), cand.end(), [&D](int a, int b){ return D[a] > D[b]; });
   float md_grid = merge_dist * n_grids_per_angstrom;
   float md2 = md_grid * md_grid;
   std::vector<int> seeds;
   for (unsigned int i=0; i<cand.size(); i++) {
      triple_index_t p = deindex(cand[i]);
      bool near = false;
      for (unsigned int k=0; k<seeds.size(); k++) {
         triple_index_t q = deindex(seeds[k]);
         float dx=p.ix-q.ix, dy=p.iy-q.iy, dz=p.iz-q.iz;
         if (dx*dx+dy*dy+dz*dz <= md2) { near = true; break; }
      }
      if (! near) seeds.push_back(cand[i]);
   }

   if (debug)
      std::cout << "   [segment cavity " << c.id << "] " << cand.size()
                << " candidate maxima -> " << seeds.size() << " seeds" << std::endl;
   if (seeds.size() <= 1) return subs; // single chamber - no split

   // priority-flood watershed from the seeds (highest D first). Seeds are pre-claimed
   // so a lower-D seed cannot be swallowed by a neighbour's flood before it is popped
   // (which used to leave empty basins). Each boundary cell joins the basin of its
   // highest-D neighbour (assign-on-expansion).
   std::vector<int> lab(n, 0);
   std::priority_queue<std::tuple<int,int,int> > pq; // (D, idx, label)
   for (unsigned int i=0; i<seeds.size(); i++) {
      lab[seeds[i]] = static_cast<int>(i+1);
      pq.push(std::make_tuple(D[seeds[i]], seeds[i], static_cast<int>(i+1)));
   }
   while (! pq.empty()) {
      std::tuple<int,int,int> tp = pq.top(); pq.pop();
      int idx = std::get<1>(tp);
      int L   = std::get<2>(tp);
      triple_index_t t = deindex(idx);
      for (int d=0; d<6; d++) {
         int gx=t.ix+d6[d][0], gy=t.iy+d6[d][1], gz=t.iz+d6[d][2];
         if (gx<0||gx>=nx||gy<0||gy>=ny||gz<0||gz>=nz) continue;
         int nidx = grid_index(triple_index_t(gx,gy,gz));
         if (in_cav[nidx] && lab[nidx]==0) {
            lab[nidx] = L;
            pq.push(std::make_tuple(D[nidx], nidx, L));
         }
      }
   }

   int k = static_cast<int>(seeds.size());

   // --- persistence-based basin merging (fight watershed over-segmentation) ---
   // Merge two basins when the neck between them (their saddle: the deepest point
   // at which they connect) is nearly as deep as the shallower basin's peak - i.e.
   // low prominence. A thin neck between two real chambers gives a low saddle and
   // high prominence, so those survive.
   std::vector<int> peak(k+1, 0);
   for (int i=0; i<k; i++) peak[i+1] = D[seeds[i]];

   std::vector<std::vector<int> > sad(k+1, std::vector<int>(k+1, -1));
   for (unsigned int i=0; i<c.grid_indices.size(); i++) {
      int idx = c.grid_indices[i];
      int a = lab[idx];
      if (a < 1) continue;
      triple_index_t t = deindex(idx);
      for (int d=0; d<6; d++) {
         int gx=t.ix+d6[d][0], gy=t.iy+d6[d][1], gz=t.iz+d6[d][2];
         if (gx<0||gx>=nx||gy<0||gy>=ny||gz<0||gz>=nz) continue;
         int nidx = grid_index(triple_index_t(gx,gy,gz));
         if (! in_cav[nidx]) continue;
         int b = lab[nidx];
         if (b >= 1 && b != a) {
            int s = std::min(D[idx], D[nidx]);
            if (s > sad[a][b]) { sad[a][b] = s; sad[b][a] = s; }
         }
      }
   }

   std::vector<int> uf(k+1);
   for (int i=0; i<=k; i++) uf[i] = i;
   float min_prom_grid = min_prominence * n_grids_per_angstrom * CH_SCALE; // A -> chamfer units
   int n_basins = k;
   while (true) {
      std::vector<int> peak_root(k+1, 0);
      for (int b=1; b<=k; b++) {
         int r=b; while (uf[r]!=r) { uf[r]=uf[uf[r]]; r=uf[r]; }
         if (peak[b] > peak_root[r]) peak_root[r] = peak[b];
      }
      std::map<std::pair<int,int>, int> root_sad; // (root_a,root_b) -> highest saddle
      for (int a=1; a<=k; a++) {
         int ra=a; while (uf[ra]!=ra) { uf[ra]=uf[uf[ra]]; ra=uf[ra]; }
         for (int b=a+1; b<=k; b++) {
            if (sad[a][b] < 0) continue;
            int rb=b; while (uf[rb]!=rb) { uf[rb]=uf[uf[rb]]; rb=uf[rb]; }
            if (ra==rb) continue;
            std::pair<int,int> key(std::min(ra,rb), std::max(ra,rb));
            std::map<std::pair<int,int>,int>::iterator it = root_sad.find(key);
            if (it==root_sad.end() || sad[a][b] > it->second) root_sad[key] = sad[a][b];
         }
      }
      float best_prom = 1.0e9;
      std::pair<int,int> best_key(-1,-1);
      for (std::map<std::pair<int,int>,int>::iterator it=root_sad.begin(); it!=root_sad.end(); ++it) {
         int ra=it->first.first, rb=it->first.second;
         float prom = static_cast<float>(std::min(peak_root[ra], peak_root[rb]) - it->second);
         if (prom < best_prom) { best_prom = prom; best_key = it->first; }
      }
      if (best_key.first < 0) break;           // nothing adjacent left
      if (best_prom >= min_prom_grid) break;   // all remaining basins prominent enough
      int ra = best_key.first, rb = best_key.second;
      if (debug)
         std::cout << "      merge at prominence "
                   << (best_prom/(n_grids_per_angstrom*CH_SCALE)) << " A  ("
                   << n_basins << " -> " << (n_basins-1) << ")" << std::endl;
      n_basins--;
      std::vector<int> pr(k+1, 0);             // peaks again for the tie decision
      for (int b=1; b<=k; b++) { int r=b; while(uf[r]!=r){uf[r]=uf[uf[r]];r=uf[r];} if(peak[b]>pr[r]) pr[r]=peak[b]; }
      if (pr[ra] >= pr[rb]) uf[rb] = ra; else uf[ra] = rb; // shallower joins deeper
   }

   // compact surviving roots to 1..m and relabel the voxels
   std::vector<int> root_to_new(k+1, 0);
   int m = 0;
   for (int b=1; b<=k; b++) {
      int r=b; while (uf[r]!=r) { uf[r]=uf[uf[r]]; r=uf[r]; }
      if (root_to_new[r]==0) root_to_new[r] = ++m;
   }
   for (unsigned int i=0; i<c.grid_indices.size(); i++) {
      int idx = c.grid_indices[i];
      int L = lab[idx];
      if (L < 1) continue;
      int r=L; while (uf[r]!=r) { uf[r]=uf[uf[r]]; r=uf[r]; }
      lab[idx] = root_to_new[r];
   }
   if (m <= 1) return subs; // merged back to a single chamber - no split

   // build subpockets from the merged labels
   subs.resize(m);
   for (int i=0; i<m; i++) { subs[i].parent_id = c.id; subs[i].id = i+1; }
   for (unsigned int i=0; i<c.grid_indices.size(); i++) {
      int idx = c.grid_indices[i];
      int L = lab[idx];
      if (L >= 1 && L <= m) subs[L-1].grid_indices.push_back(idx);
   }

   // per-subpocket volume, centre, surface area (a face is exposed if its neighbour
   // is not in the *same* subpocket - so the neck cut counts as surface).
   float voxel_volume = 1.0f / (n_grids_per_angstrom*n_grids_per_angstrom*n_grids_per_angstrom);
   float face_area    = 1.0f / (n_grids_per_angstrom*n_grids_per_angstrom);
   for (int i=0; i<m; i++) {
      cavity_t &s = subs[i];
      int my = i+1;
      s.volume = static_cast<float>(s.grid_indices.size()) * voxel_volume;
      double sx=0.0, sy=0.0, sz=0.0;
      unsigned int faces=0;
      for (unsigned int j=0; j<s.grid_indices.size(); j++) {
         triple_index_t t = deindex(s.grid_indices[j]);
         point_3d_t p = grid_point_to_mol_space(t);
         sx+=p.x; sy+=p.y; sz+=p.z;
         for (int d=0; d<6; d++) {
            int gx=t.ix+d6[d][0], gy=t.iy+d6[d][1], gz=t.iz+d6[d][2];
            if (gx<0||gx>=nx||gy<0||gy>=ny||gz<0||gz>=nz) { faces++; continue; }
            if (lab[grid_index(triple_index_t(gx,gy,gz))] != my) faces++;
         }
      }
      double np = static_cast<double>(s.grid_indices.size());
      if (np > 0.0) s.centre = point_3d_t(float(sx/np), float(sy/np), float(sz/np));
      s.surface_area = static_cast<float>(faces) * face_area;
   }

   return subs;
}

// Segment every cavity above min_volume_to_segment into subpockets, collect them in
// the subpockets member, compute their lining residues, and report.
void
coot::grid_balls_t::segment_cavities(float min_volume_to_segment, int min_seed_depth,
                                     float merge_dist, float min_prominence) {

   subpockets.clear();
   for (unsigned int ic=0; ic<cavities.size(); ic++) {
      const cavity_t &c = cavities[ic];
      if (c.volume < min_volume_to_segment) continue;
      std::vector<cavity_t> subs = segment_one_cavity(c, min_seed_depth, merge_dist, min_prominence);
      for (unsigned int i=0; i<subs.size(); i++)
         subpockets.push_back(subs[i]);
   }

   compute_lining_residues_for(subpockets, 5.0);

   std::cout << "segment_cavities(): " << subpockets.size() << " subpockets" << std::endl;
   for (unsigned int i=0; i<subpockets.size(); i++) {
      const cavity_t &s = subpockets[i];
      std::cout << "   subpocket " << s.parent_id << "." << s.id
                << ": " << s.grid_indices.size() << " points, volume " << s.volume
                << " A^3, centre (" << s.centre.x << " " << s.centre.y << " " << s.centre.z << ")"
                << std::endl;
      print_lining_residues(s, "subpocket " + std::to_string(s.parent_id) + "." + std::to_string(s.id));
   }
}

// Write a plain table of subpocket point positions (colour_id x y z), same format as
// write_cavity_points, with a running colour id over all subpockets.
void
coot::grid_balls_t::write_subpocket_points(const std::string &file_name) const {

   std::ofstream f(file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: write_subpocket_points() cannot open " << file_name << std::endl;
      return;
   }

   f << std::fixed << std::setprecision(3);
   f << "# grid-balls subpocket points\n";
   f << "# " << subpockets.size() << " subpockets\n";
   f << "# columns: colour_id x y z  (colour_id runs over all subpockets)\n";

   int n_points = 0;
   for (unsigned int i=0; i<subpockets.size(); i++) {
      const cavity_t &s = subpockets[i];
      int colour_id = i + 1;
      f << "# subpocket " << colour_id << "  parent " << s.parent_id << "  sub " << s.id
        << "  points " << s.grid_indices.size() << "  volume " << s.volume
        << "  centre " << s.centre.x << " " << s.centre.y << " " << s.centre.z << "\n";
      for (unsigned int j=0; j<s.grid_indices.size(); j++) {
         point_3d_t p = grid_point_to_mol_space(deindex(s.grid_indices[j]));
         f << colour_id << " " << p.x << " " << p.y << " " << p.z << "\n";
         n_points++;
      }
   }
   f.close();
   std::cout << "write_subpocket_points(): wrote " << n_points << " points for "
             << subpockets.size() << " subpockets to " << file_name << std::endl;
}
