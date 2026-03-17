/*
 * coot-utils/assembly-interface.cc
 *
 * Copyright 2026 by Medical Research Council
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

#include <iostream>
#include <sstream>
#include <cmath>
#include <set>
#include <algorithm>
#include <iomanip>

#include <clipper/core/clipper_util.h>

#include "fast-eigens.hh"
#include "assembly-interface.hh"

// -----------------------------------------------------------------------
// ASP (Atomic Solvation Parameter) values
// From Table 1 of Krissinel & Henrick (2007), column A
// Units: cal/(mol * A^2)
// -----------------------------------------------------------------------

double
coot::asp_value(asp_atom_type_t type) {
   switch (type) {
      case asp_atom_type_t::C:       return  16.0;
      case asp_atom_type_t::S:       return  41.0;
      case asp_atom_type_t::O_N:     return -11.0;
      case asp_atom_type_t::N_PLUS:  return -37.0;
      case asp_atom_type_t::O_MINUS: return -17.0;
   }
   return 0.0;
}

// -----------------------------------------------------------------------
// Classify an atom into one of 5 ASP atom types
// Following the Eisenberg & McLachlan scheme:
//   C      - carbon atoms
//   S      - sulfur atoms
//   O/N    - neutral oxygen and nitrogen
//   N+     - charged nitrogen (Arg, Lys, His N-terminal)
//   O-     - charged oxygen (Asp, Glu C-terminal)
// -----------------------------------------------------------------------

coot::asp_atom_type_t
coot::classify_asp_atom_type(mmdb::Atom *atom_p) {

   std::string element(atom_p->element);
   // trim whitespace
   while (!element.empty() && element[0] == ' ') element.erase(0, 1);
   while (!element.empty() && element.back() == ' ') element.pop_back();

   if (element == "C") return asp_atom_type_t::C;
   if (element == "S") return asp_atom_type_t::S;

   if (element == "N") {
      // Check if this is a charged nitrogen
      std::string resname(atom_p->GetResName());
      std::string atname(atom_p->GetAtomName());
      // trim
      while (!atname.empty() && atname[0] == ' ') atname.erase(0, 1);
      while (!atname.empty() && atname.back() == ' ') atname.pop_back();

      // Arg: NH1, NH2 (guanidinium, pKa ~12.5)
      if (resname == "ARG" && (atname == "NH1" || atname == "NH2" || atname == "NE"))
         return asp_atom_type_t::N_PLUS;
      // Lys: NZ (amino, pKa ~10.5)
      if (resname == "LYS" && atname == "NZ")
         return asp_atom_type_t::N_PLUS;
      // His: ND1, NE2 (imidazole, pKa ~6.0, partially charged at physiological pH)
      if (resname == "HIS" && (atname == "ND1" || atname == "NE2"))
         return asp_atom_type_t::N_PLUS;
      // N-terminus
      if (atname == "N") {
         // Could check if it's the first residue, but for simplicity
         // treat backbone N as neutral
      }
      return asp_atom_type_t::O_N;
   }

   if (element == "O") {
      std::string resname(atom_p->GetResName());
      std::string atname(atom_p->GetAtomName());
      while (!atname.empty() && atname[0] == ' ') atname.erase(0, 1);
      while (!atname.empty() && atname.back() == ' ') atname.pop_back();

      // Asp: OD1, OD2 (carboxylate, pKa ~3.9)
      if (resname == "ASP" && (atname == "OD1" || atname == "OD2"))
         return asp_atom_type_t::O_MINUS;
      // Glu: OE1, OE2 (carboxylate, pKa ~4.1)
      if (resname == "GLU" && (atname == "OE1" || atname == "OE2"))
         return asp_atom_type_t::O_MINUS;
      // C-terminal OXT
      if (atname == "OXT")
         return asp_atom_type_t::O_MINUS;
      return asp_atom_type_t::O_N;
   }

   // Default: treat as carbon (e.g. Se, metals, etc.)
   return asp_atom_type_t::C;
}

// -----------------------------------------------------------------------
// VDW radii for ASA calculation
// -----------------------------------------------------------------------

double
coot::vdw_radius(const std::string &element) {
   // Standard radii used in ASA calculations
   if (element == "C")  return 1.70;
   if (element == "N")  return 1.55;
   if (element == "O")  return 1.52;
   if (element == "S")  return 1.80;
   if (element == "P")  return 1.80;
   if (element == "H")  return 1.20;
   if (element == "SE") return 1.90;
   if (element == "FE") return 1.47;
   if (element == "ZN") return 1.39;
   if (element == "MG") return 1.73;
   if (element == "CA") return 1.97;
   if (element == "MN") return 1.39;
   if (element == "CU") return 1.40;
   if (element == "CO") return 1.26;
   if (element == "NI") return 1.63;
   return 1.70;  // default
}

// -----------------------------------------------------------------------
// monomer_id_t methods
// -----------------------------------------------------------------------

bool
coot::monomer_id_t::operator<(const monomer_id_t &o) const {
   if (chain_id != o.chain_id) return chain_id < o.chain_id;
   if (symop_no != o.symop_no) return symop_no < o.symop_no;
   if (x_shift != o.x_shift) return x_shift < o.x_shift;
   if (y_shift != o.y_shift) return y_shift < o.y_shift;
   return z_shift < o.z_shift;
}

std::string
coot::monomer_id_t::label() const {
   std::ostringstream oss;
   oss << chain_id << "[" << symop_no;
   if (x_shift != 0 || y_shift != 0 || z_shift != 0)
      oss << "," << x_shift << "," << y_shift << "," << z_shift;
   oss << "]";
   return oss.str();
}

// -----------------------------------------------------------------------
// Fibonacci sphere - generate uniformly distributed points on unit sphere
// -----------------------------------------------------------------------

std::vector<clipper::Coord_orth>
coot::fibonacci_sphere_points(int n_points) {

   std::vector<clipper::Coord_orth> points(n_points);
   double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
   double angle_increment = 2.0 * M_PI / golden_ratio;

   for (int i=0; i<n_points; i++) {
      double y = 1.0 - (2.0 * i + 1.0) / n_points;
      double radius = std::sqrt(1.0 - y * y);
      double theta = angle_increment * i;
      double x = std::cos(theta) * radius;
      double z = std::sin(theta) * radius;
      points[i] = clipper::Coord_orth(x, y, z);
   }
   return points;
}

// -----------------------------------------------------------------------
// Shrake-Rupley accessible surface area
// For each atom, count the fraction of test points on its solvent sphere
// that are not buried inside any neighbouring atom's solvent sphere.
// -----------------------------------------------------------------------

std::vector<double>
coot::shrake_rupley_asa(const std::vector<clipper::Coord_orth> &positions,
                        const std::vector<double> &radii,
                        double probe_radius,
                        int n_sphere_points) {

   int n_atoms = positions.size();
   std::vector<double> asa(n_atoms, 0.0);
   if (n_atoms == 0) return asa;

   std::vector<clipper::Coord_orth> sphere_pts = fibonacci_sphere_points(n_sphere_points);
   double area_per_point_factor = 4.0 * M_PI / n_sphere_points;

   // Build a simple spatial grid for neighbour lookup
   // Find the maximum extended radius for the neighbour search
   double max_radius = 0.0;
   for (int i=0; i<n_atoms; i++) {
      double r = radii[i] + probe_radius;
      if (r > max_radius) max_radius = r;
   }
   double neighbour_cutoff = 2.0 * max_radius;
   double neighbour_cutoff_sq = neighbour_cutoff * neighbour_cutoff;

   for (int i=0; i<n_atoms; i++) {
      double ri = radii[i] + probe_radius;
      double ri_sq_area = ri * ri * area_per_point_factor;

      // Collect neighbours of atom i
      std::vector<int> neighbours;
      for (int j=0; j<n_atoms; j++) {
         if (j == i) continue;
         double dx = positions[j].x() - positions[i].x();
         double dy = positions[j].y() - positions[i].y();
         double dz = positions[j].z() - positions[i].z();
         double d_sq = dx*dx + dy*dy + dz*dz;
         if (d_sq < neighbour_cutoff_sq)
            neighbours.push_back(j);
      }

      int n_accessible = 0;
      for (int s=0; s<n_sphere_points; s++) {
         double px = positions[i].x() + ri * sphere_pts[s].x();
         double py = positions[i].y() + ri * sphere_pts[s].y();
         double pz = positions[i].z() + ri * sphere_pts[s].z();

         bool buried = false;
         for (unsigned int jn=0; jn<neighbours.size(); jn++) {
            int j = neighbours[jn];
            double rj = radii[j] + probe_radius;
            double dx = px - positions[j].x();
            double dy = py - positions[j].y();
            double dz = pz - positions[j].z();
            if (dx*dx + dy*dy + dz*dz < rj * rj) {
               buried = true;
               break;
            }
         }
         if (!buried) n_accessible++;
      }
      asa[i] = n_accessible * ri_sq_area;
   }
   return asa;
}

// -----------------------------------------------------------------------
// Molecular mass calculation
// -----------------------------------------------------------------------

double
coot::molecular_mass(mmdb::Manager *mol_p, int selHnd) {

   mmdb::PPAtom atoms = nullptr;
   int n_atoms = 0;
   mol_p->GetSelIndex(selHnd, atoms, n_atoms);

   double mass = 0.0;
   for (int i=0; i<n_atoms; i++) {
      std::string el(atoms[i]->element);
      while (!el.empty() && el[0] == ' ') el.erase(0, 1);
      while (!el.empty() && el.back() == ' ') el.pop_back();
      if      (el == "C")  mass += 12.011;
      else if (el == "N")  mass += 14.007;
      else if (el == "O")  mass += 15.999;
      else if (el == "S")  mass += 32.065;
      else if (el == "P")  mass += 30.974;
      else if (el == "H")  mass +=  1.008;
      else if (el == "SE") mass += 78.96;
      else if (el == "FE") mass += 55.845;
      else if (el == "ZN") mass += 65.38;
      else if (el == "MG") mass += 24.305;
      else if (el == "CA") mass += 40.078;
      else if (el == "MN") mass += 54.938;
      else if (el == "CU") mass += 63.546;
      else if (el == "NA") mass += 22.990;
      else if (el == "K")  mass += 39.098;
      else if (el == "CL") mass += 35.453;
      else                  mass += 12.0;  // default approximate
   }
   return mass;
}

// -----------------------------------------------------------------------
// Principal moments of inertia
// -----------------------------------------------------------------------

std::tuple<double, double, double>
coot::principal_moments_of_inertia(const std::vector<clipper::Coord_orth> &positions,
                                   const std::vector<double> &masses) {

   int n = positions.size();
   if (n == 0) return std::make_tuple(0.0, 0.0, 0.0);

   // Centre of mass
   double total_mass = 0.0;
   double cx = 0.0, cy = 0.0, cz = 0.0;
   for (int i=0; i<n; i++) {
      cx += masses[i] * positions[i].x();
      cy += masses[i] * positions[i].y();
      cz += masses[i] * positions[i].z();
      total_mass += masses[i];
   }
   if (total_mass < 1e-10) return std::make_tuple(0.0, 0.0, 0.0);
   cx /= total_mass;
   cy /= total_mass;
   cz /= total_mass;

   // Inertia tensor
   double Ixx = 0.0, Iyy = 0.0, Izz = 0.0;
   double Ixy = 0.0, Ixz = 0.0, Iyz = 0.0;
   for (int i=0; i<n; i++) {
      double dx = positions[i].x() - cx;
      double dy = positions[i].y() - cy;
      double dz = positions[i].z() - cz;
      double m = masses[i];
      Ixx += m * (dy*dy + dz*dz);
      Iyy += m * (dx*dx + dz*dz);
      Izz += m * (dx*dx + dy*dy);
      Ixy -= m * dx * dy;
      Ixz -= m * dx * dz;
      Iyz -= m * dy * dz;
   }

   // Diagonalise using Jacobi method (from fast-eigens)
   clipper::Matrix<double> mat(3, 3);
   mat(0,0) = Ixx;  mat(0,1) = Ixy;  mat(0,2) = Ixz;
   mat(1,0) = Ixy;  mat(1,1) = Iyy;  mat(1,2) = Iyz;
   mat(2,0) = Ixz;  mat(2,1) = Iyz;  mat(2,2) = Izz;

   auto [e1, e2, e3] = fast_eigens(mat, true);
   return std::make_tuple(e1, e2, e3);
}

// -----------------------------------------------------------------------
// assembly_detector_t constructor
// -----------------------------------------------------------------------

coot::assembly_detector_t::assembly_detector_t(mmdb::Manager *mol_in, double contact_dist_in)
   : mol_p(mol_in), contact_dist(contact_dist_in) {
}

// -----------------------------------------------------------------------
// Generate symmetry-related chains within contact distance of the ASU
// -----------------------------------------------------------------------

void
coot::assembly_detector_t::generate_symmetry_chains() {

   symmetry_chains.clear();

   if (!mol_p) return;

   // First, get the ASU chains as the reference set
   int imod = 1;
   mmdb::Model *model_p = mol_p->GetModel(imod);
   if (!model_p) return;

   int n_chains = model_p->GetNumberOfChains();
   if (n_chains == 0) return;

   // Calculate the extent of the ASU for neighbour search
   double x_min = 1e30, x_max = -1e30;
   double y_min = 1e30, y_max = -1e30;
   double z_min = 1e30, z_max = -1e30;

   int all_sel = mol_p->NewSelection();
   mol_p->SelectAtoms(all_sel, 1, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                      "*", "*", "*", "*");
   mmdb::PPAtom all_atoms = nullptr;
   int n_all = 0;
   mol_p->GetSelIndex(all_sel, all_atoms, n_all);

   for (int i=0; i<n_all; i++) {
      if (all_atoms[i]->x < x_min) x_min = all_atoms[i]->x;
      if (all_atoms[i]->x > x_max) x_max = all_atoms[i]->x;
      if (all_atoms[i]->y < y_min) y_min = all_atoms[i]->y;
      if (all_atoms[i]->y > y_max) y_max = all_atoms[i]->y;
      if (all_atoms[i]->z < z_min) z_min = all_atoms[i]->z;
      if (all_atoms[i]->z > z_max) z_max = all_atoms[i]->z;
   }

   // Add padding for contact distance
   double pad = contact_dist + 5.0;
   x_min -= pad; x_max += pad;
   y_min -= pad; y_max += pad;
   z_min -= pad; z_max += pad;

   // First store the ASU chains themselves
   for (int ic=0; ic<n_chains; ic++) {
      mmdb::Chain *chain_p = model_p->GetChain(ic);
      std::string chain_id(chain_p->GetChainID());

      symmetry_chain_t sc;
      sc.id = monomer_id_t(chain_id, 0, 0, 0, 0);

      int chain_sel = mol_p->NewSelection();
      mol_p->SelectAtoms(chain_sel, 1, chain_id.c_str(),
                         mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                         "*", "*", "*", "*");
      mmdb::PPAtom chain_atoms = nullptr;
      int n_chain_atoms = 0;
      mol_p->GetSelIndex(chain_sel, chain_atoms, n_chain_atoms);

      sc.positions.resize(n_chain_atoms);
      sc.source_atoms.resize(n_chain_atoms);
      for (int i=0; i<n_chain_atoms; i++) {
         sc.positions[i] = clipper::Coord_orth(chain_atoms[i]->x,
                                                chain_atoms[i]->y,
                                                chain_atoms[i]->z);
         sc.source_atoms[i] = chain_atoms[i];
      }
      // Identity transform
      for (int r=0; r<4; r++)
         for (int c=0; c<4; c++)
            sc.transform[r][c] = (r == c) ? 1.0 : 0.0;

      symmetry_chains.push_back(sc);
      mol_p->DeleteSelection(chain_sel);
   }

   // Now generate symmetry mates
   int n_symops = mol_p->GetNumberOfSymOps();
   if (n_symops <= 0) {
      mol_p->DeleteSelection(all_sel);
      return;  // no symmetry information available
   }

   // Loop over symmetry operations and cell translations
   for (int isym=0; isym<n_symops; isym++) {
      for (int ix=-1; ix<=1; ix++) {
         for (int iy=-1; iy<=1; iy++) {
            for (int iz=-1; iz<=1; iz++) {
               // Skip identity
               if (isym == 0 && ix == 0 && iy == 0 && iz == 0) continue;

               mmdb::mat44 tmat;
               int err = mol_p->GetTMatrix(tmat, isym, ix, iy, iz);
               if (err != 0) continue;

               // For each chain in the ASU
               for (int ic=0; ic<n_chains; ic++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ic);
                  std::string chain_id(chain_p->GetChainID());

                  int chain_sel = mol_p->NewSelection();
                  mol_p->SelectAtoms(chain_sel, 1, chain_id.c_str(),
                                     mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                                     "*", "*", "*", "*");
                  mmdb::PPAtom chain_atoms = nullptr;
                  int n_chain_atoms = 0;
                  mol_p->GetSelIndex(chain_sel, chain_atoms, n_chain_atoms);

                  if (n_chain_atoms == 0) {
                     mol_p->DeleteSelection(chain_sel);
                     continue;
                  }

                  // Apply transformation and check if any atom is within range
                  std::vector<clipper::Coord_orth> transformed(n_chain_atoms);
                  bool any_in_range = false;
                  for (int i=0; i<n_chain_atoms; i++) {
                     double tx = tmat[0][0]*chain_atoms[i]->x +
                                 tmat[0][1]*chain_atoms[i]->y +
                                 tmat[0][2]*chain_atoms[i]->z + tmat[0][3];
                     double ty = tmat[1][0]*chain_atoms[i]->x +
                                 tmat[1][1]*chain_atoms[i]->y +
                                 tmat[1][2]*chain_atoms[i]->z + tmat[1][3];
                     double tz = tmat[2][0]*chain_atoms[i]->x +
                                 tmat[2][1]*chain_atoms[i]->y +
                                 tmat[2][2]*chain_atoms[i]->z + tmat[2][3];
                     transformed[i] = clipper::Coord_orth(tx, ty, tz);

                     if (tx >= x_min && tx <= x_max &&
                         ty >= y_min && ty <= y_max &&
                         tz >= z_min && tz <= z_max)
                        any_in_range = true;
                  }

                  if (any_in_range) {
                     symmetry_chain_t sc;
                     sc.id = monomer_id_t(chain_id, isym, ix, iy, iz);
                     sc.positions = transformed;
                     sc.source_atoms.resize(n_chain_atoms);
                     for (int i=0; i<n_chain_atoms; i++)
                        sc.source_atoms[i] = chain_atoms[i];
                     for (int r=0; r<4; r++)
                        for (int c=0; c<4; c++)
                           sc.transform[r][c] = tmat[r][c];
                     symmetry_chains.push_back(sc);
                  }

                  mol_p->DeleteSelection(chain_sel);
               }
            }
         }
      }
   }
   mol_p->DeleteSelection(all_sel);

   std::cout << "INFO:: assembly_detector: generated " << symmetry_chains.size()
             << " chains (" << n_chains << " ASU + "
             << (symmetry_chains.size() - n_chains) << " symmetry mates)"
             << std::endl;
}

// -----------------------------------------------------------------------
// Find interfaces between all pairs of chains
// -----------------------------------------------------------------------

void
coot::assembly_detector_t::find_interfaces() {

   interfaces.clear();
   double contact_dist_sq = contact_dist * contact_dist;

   int n_chains = symmetry_chains.size();
   for (int i=0; i<n_chains; i++) {
      for (int j=i+1; j<n_chains; j++) {
         const symmetry_chain_t &ci = symmetry_chains[i];
         const symmetry_chain_t &cj = symmetry_chains[j];

         // Quick check: do any atoms come within contact distance?
         bool have_contact = false;
         for (unsigned int ai=0; ai<ci.positions.size() && !have_contact; ai++) {
            for (unsigned int aj=0; aj<cj.positions.size() && !have_contact; aj++) {
               double dx = ci.positions[ai].x() - cj.positions[aj].x();
               double dy = ci.positions[ai].y() - cj.positions[aj].y();
               double dz = ci.positions[ai].z() - cj.positions[aj].z();
               if (dx*dx + dy*dy + dz*dz < contact_dist_sq)
                  have_contact = true;
            }
         }

         if (have_contact) {
            assembly_interface_t iface;
            iface.monomer_1 = ci.id;
            iface.monomer_2 = cj.id;

            calculate_interface_bsa(&iface, ci, cj);
            count_interface_contacts(&iface, ci, cj);
            calculate_interface_energetics(&iface);

            // Only keep interfaces with non-trivial BSA
            if (iface.buried_surface_area > 10.0) {
               interfaces.push_back(iface);
            }
         }
      }
   }

   // Sort by BSA descending
   std::sort(interfaces.begin(), interfaces.end(),
             [](const assembly_interface_t &a, const assembly_interface_t &b) {
                return a.buried_surface_area > b.buried_surface_area;
             });

   std::cout << "INFO:: assembly_detector: found " << interfaces.size()
             << " interfaces with BSA > 10 A^2" << std::endl;
}

// -----------------------------------------------------------------------
// Calculate buried surface area for an interface
// BSA = ASA(chain_i) + ASA(chain_j) - ASA(chain_i + chain_j)
// Also compute per-atom-type BSA for solvation energy calculation
// -----------------------------------------------------------------------

void
coot::assembly_detector_t::calculate_interface_bsa(assembly_interface_t *iface_p,
                                                   const symmetry_chain_t &chain_1,
                                                   const symmetry_chain_t &chain_2) {

   int n1 = chain_1.positions.size();
   int n2 = chain_2.positions.size();
   int n_total = n1 + n2;

   // Collect positions and radii for both chains
   std::vector<clipper::Coord_orth> all_pos(n_total);
   std::vector<double> all_radii(n_total);

   for (int i=0; i<n1; i++) {
      all_pos[i] = chain_1.positions[i];
      std::string el(chain_1.source_atoms[i]->element);
      while (!el.empty() && el[0] == ' ') el.erase(0, 1);
      while (!el.empty() && el.back() == ' ') el.pop_back();
      all_radii[i] = vdw_radius(el);
   }
   for (int i=0; i<n2; i++) {
      all_pos[n1 + i] = chain_2.positions[i];
      std::string el(chain_2.source_atoms[i]->element);
      while (!el.empty() && el[0] == ' ') el.erase(0, 1);
      while (!el.empty() && el.back() == ' ') el.pop_back();
      all_radii[n1 + i] = vdw_radius(el);
   }

   // ASA of chains separately
   std::vector<clipper::Coord_orth> pos1(chain_1.positions.begin(), chain_1.positions.end());
   std::vector<double> rad1(all_radii.begin(), all_radii.begin() + n1);
   std::vector<double> asa1 = shrake_rupley_asa(pos1, rad1);

   std::vector<clipper::Coord_orth> pos2(chain_2.positions.begin(), chain_2.positions.end());
   std::vector<double> rad2(all_radii.begin() + n1, all_radii.end());
   std::vector<double> asa2 = shrake_rupley_asa(pos2, rad2);

   // ASA of complex
   std::vector<double> asa_complex = shrake_rupley_asa(all_pos, all_radii);

   // BSA = ASA_separate - ASA_complex (per atom)
   double total_bsa = 0.0;
   std::map<asp_atom_type_t, double> bsa_by_type;
   bsa_by_type[asp_atom_type_t::C]       = 0.0;
   bsa_by_type[asp_atom_type_t::S]       = 0.0;
   bsa_by_type[asp_atom_type_t::O_N]     = 0.0;
   bsa_by_type[asp_atom_type_t::N_PLUS]  = 0.0;
   bsa_by_type[asp_atom_type_t::O_MINUS] = 0.0;

   clipper::Coord_orth centre_sum(0, 0, 0);
   int n_interface_atoms = 0;
   std::set<mmdb::Residue *> iface_residues_1;
   std::set<mmdb::Residue *> iface_residues_2;

   // Chain 1 atoms
   for (int i=0; i<n1; i++) {
      double delta = asa1[i] - asa_complex[i];
      if (delta > 0.1) {
         total_bsa += delta;
         asp_atom_type_t atype = classify_asp_atom_type(chain_1.source_atoms[i]);
         bsa_by_type[atype] += delta;
         centre_sum = clipper::Coord_orth(centre_sum.x() + chain_1.positions[i].x(),
                                          centre_sum.y() + chain_1.positions[i].y(),
                                          centre_sum.z() + chain_1.positions[i].z());
         n_interface_atoms++;
         iface_residues_1.insert(chain_1.source_atoms[i]->GetResidue());
      }
   }

   // Chain 2 atoms
   for (int i=0; i<n2; i++) {
      double delta = asa2[i] - asa_complex[n1 + i];
      if (delta > 0.1) {
         total_bsa += delta;
         asp_atom_type_t atype = classify_asp_atom_type(chain_2.source_atoms[i]);
         bsa_by_type[atype] += delta;
         centre_sum = clipper::Coord_orth(centre_sum.x() + chain_2.positions[i].x(),
                                          centre_sum.y() + chain_2.positions[i].y(),
                                          centre_sum.z() + chain_2.positions[i].z());
         n_interface_atoms++;
         iface_residues_2.insert(chain_2.source_atoms[i]->GetResidue());
      }
   }

   iface_p->buried_surface_area = total_bsa;
   iface_p->buried_area_by_type = bsa_by_type;
   iface_p->n_interface_atoms_1 = 0;
   iface_p->n_interface_atoms_2 = 0;
   iface_p->n_interface_residues_1 = iface_residues_1.size();
   iface_p->n_interface_residues_2 = iface_residues_2.size();

   if (n_interface_atoms > 0) {
      iface_p->centre = clipper::Coord_orth(centre_sum.x() / n_interface_atoms,
                                             centre_sum.y() / n_interface_atoms,
                                             centre_sum.z() / n_interface_atoms);
   }
}

// -----------------------------------------------------------------------
// Count hydrogen bonds, salt bridges, and disulfide bonds at interface
// -----------------------------------------------------------------------

void
coot::assembly_detector_t::count_interface_contacts(assembly_interface_t *iface_p,
                                                    const symmetry_chain_t &chain_1,
                                                    const symmetry_chain_t &chain_2) {

   int n_hb = 0;
   int n_sb = 0;
   int n_ss = 0;

   int n1 = chain_1.positions.size();
   int n2 = chain_2.positions.size();

   // H-bond criteria: donor-acceptor distance 2.5-3.5 A
   // Salt bridge: charged N...charged O distance < 4.0 A
   // Disulfide: S...S distance 1.8-2.5 A

   for (int i=0; i<n1; i++) {
      mmdb::Atom *at1 = chain_1.source_atoms[i];
      std::string el1(at1->element);
      while (!el1.empty() && el1[0] == ' ') el1.erase(0, 1);
      while (!el1.empty() && el1.back() == ' ') el1.pop_back();

      for (int j=0; j<n2; j++) {
         mmdb::Atom *at2 = chain_2.source_atoms[j];
         std::string el2(at2->element);
         while (!el2.empty() && el2[0] == ' ') el2.erase(0, 1);
         while (!el2.empty() && el2.back() == ' ') el2.pop_back();

         double dx = chain_1.positions[i].x() - chain_2.positions[j].x();
         double dy = chain_1.positions[i].y() - chain_2.positions[j].y();
         double dz = chain_1.positions[i].z() - chain_2.positions[j].z();
         double dist_sq = dx*dx + dy*dy + dz*dz;

         // Disulfide bond: S-S, 1.8-2.5 A
         if (el1 == "S" && el2 == "S") {
            std::string an1(at1->GetAtomName());
            std::string an2(at2->GetAtomName());
            // trim
            while (!an1.empty() && an1[0] == ' ') an1.erase(0, 1);
            while (!an1.empty() && an1.back() == ' ') an1.pop_back();
            while (!an2.empty() && an2[0] == ' ') an2.erase(0, 1);
            while (!an2.empty() && an2.back() == ' ') an2.pop_back();
            if (an1 == "SG" && an2 == "SG") {
               double dist = std::sqrt(dist_sq);
               if (dist > 1.8 && dist < 2.5) {
                  n_ss++;
                  continue;
               }
            }
         }

         // Salt bridge: charged N...charged O, < 4.0 A
         if (dist_sq < 16.0) {
            asp_atom_type_t t1 = classify_asp_atom_type(at1);
            asp_atom_type_t t2 = classify_asp_atom_type(at2);
            if ((t1 == asp_atom_type_t::N_PLUS && t2 == asp_atom_type_t::O_MINUS) ||
                (t1 == asp_atom_type_t::O_MINUS && t2 == asp_atom_type_t::N_PLUS)) {
               n_sb++;
               continue;
            }
         }

         // Hydrogen bond: N/O...N/O, 2.5-3.5 A (donor-acceptor)
         if (dist_sq < 12.25 && dist_sq > 6.25) {  // 3.5^2 and 2.5^2
            if ((el1 == "N" || el1 == "O") && (el2 == "N" || el2 == "O")) {
               n_hb++;
            }
         }
      }
   }

   iface_p->n_hydrogen_bonds = n_hb;
   iface_p->n_salt_bridges = n_sb;
   iface_p->n_disulfide_bonds = n_ss;
}

// -----------------------------------------------------------------------
// Calculate energetics for an interface
// Eq. 9: delta_G_int = sum over interfaces of
//   (sum_k omega_k * delta_sigma_ij^k + N_hb*E_hb + N_sb*E_sb + N_db*E_db)
// -----------------------------------------------------------------------

void
coot::assembly_detector_t::calculate_interface_energetics(assembly_interface_t *iface_p) {

   // Solvation energy: sum of ASP * buried area per atom type
   // Convert from cal/mol to kcal/mol
   double delta_g_solv = 0.0;
   for (auto &pair : iface_p->buried_area_by_type) {
      delta_g_solv += asp_value(pair.first) * pair.second;
   }
   delta_g_solv /= 1000.0;  // cal -> kcal

   // Contact energy from H-bonds, salt bridges, disulfide bonds
   // Values from Table 1 of Krissinel & Henrick (2007)
   double E_hb = 0.44;   // kcal/mol per H-bond
   double E_sb = 0.15;   // kcal/mol per salt bridge
   double E_db = 4.00;   // kcal/mol per disulfide bond

   double delta_g_cont = iface_p->n_hydrogen_bonds * E_hb
                       + iface_p->n_salt_bridges * E_sb
                       + iface_p->n_disulfide_bonds * E_db;

   // Total binding energy (negative = favourable)
   // Note: solvation term is already negative for hydrophobic burial
   // The binding energy is: delta_G_int = delta_G_solv + delta_G_cont
   // But in PISA convention, negative delta_G_int means favorable binding
   iface_p->delta_g_solv = delta_g_solv;
   iface_p->delta_g_cont = delta_g_cont;
   iface_p->delta_g_int = delta_g_solv + delta_g_cont;
}

// -----------------------------------------------------------------------
// Calculate properties of each monomeric unit
// -----------------------------------------------------------------------

void
coot::assembly_detector_t::calculate_monomer_properties() {

   monomer_props.clear();

   for (auto &sc : symmetry_chains) {
      monomer_properties_t props;
      props.id = sc.id;
      props.n_atoms = sc.positions.size();

      // Count residues
      std::set<mmdb::Residue *> residues;
      for (auto at : sc.source_atoms)
         residues.insert(at->GetResidue());
      props.n_residues = residues.size();

      // Mass
      std::vector<double> atom_masses(props.n_atoms);
      double total_mass = 0.0;
      for (int i=0; i<props.n_atoms; i++) {
         std::string el(sc.source_atoms[i]->element);
         while (!el.empty() && el[0] == ' ') el.erase(0, 1);
         while (!el.empty() && el.back() == ' ') el.pop_back();
         double m = 12.0;
         if      (el == "C")  m = 12.011;
         else if (el == "N")  m = 14.007;
         else if (el == "O")  m = 15.999;
         else if (el == "S")  m = 32.065;
         else if (el == "H")  m =  1.008;
         atom_masses[i] = m;
         total_mass += m;
      }
      props.mass = total_mass;

      // ASA
      std::vector<double> radii(props.n_atoms);
      for (int i=0; i<props.n_atoms; i++) {
         std::string el(sc.source_atoms[i]->element);
         while (!el.empty() && el[0] == ' ') el.erase(0, 1);
         while (!el.empty() && el.back() == ' ') el.pop_back();
         radii[i] = vdw_radius(el);
      }
      std::vector<double> asa = shrake_rupley_asa(sc.positions, radii);
      double total_asa = 0.0;
      for (auto a : asa) total_asa += a;
      props.surface_area = total_asa;

      // Principal moments of inertia
      auto [J1, J2, J3] = principal_moments_of_inertia(sc.positions, atom_masses);
      props.J1 = J1;
      props.J2 = J2;
      props.J3 = J3;

      // Symmetry number: default 1 (no internal symmetry)
      // Proper detection would require structural alignment, which is
      // expensive. For a first implementation, use 1.
      props.symmetry_number = 1;

      monomer_props[sc.id] = props;
   }
}

// -----------------------------------------------------------------------
// Main entry point: detect all interfaces
// -----------------------------------------------------------------------

void
coot::assembly_detector_t::detect_interfaces() {

   generate_symmetry_chains();
   find_interfaces();
   calculate_monomer_properties();
}

// -----------------------------------------------------------------------
// Summary output
// -----------------------------------------------------------------------

std::string
coot::assembly_detector_t::interfaces_summary() const {

   std::ostringstream oss;
   oss << std::fixed << std::setprecision(1);

   oss << "Assembly Interface Analysis" << std::endl;
   oss << "==========================" << std::endl;
   oss << "Total symmetry chains: " << symmetry_chains.size() << std::endl;
   oss << "Total interfaces found: " << interfaces.size() << std::endl;
   oss << std::endl;

   oss << std::setw(20) << "Interface"
       << std::setw(10) << "BSA(A^2)"
       << std::setw(8)  << "H-bonds"
       << std::setw(8)  << "SaltBr"
       << std::setw(8)  << "SS"
       << std::setw(12) << "dG_solv"
       << std::setw(12) << "dG_cont"
       << std::setw(12) << "dG_int"
       << std::endl;
   oss << std::string(90, '-') << std::endl;

   for (const auto &iface : interfaces) {
      std::string label = iface.monomer_1.label() + " - " + iface.monomer_2.label();
      oss << std::setw(20) << label
          << std::setw(10) << iface.buried_surface_area
          << std::setw(8)  << iface.n_hydrogen_bonds
          << std::setw(8)  << iface.n_salt_bridges
          << std::setw(8)  << iface.n_disulfide_bonds
          << std::setw(12) << std::setprecision(2) << iface.delta_g_solv
          << std::setw(12) << iface.delta_g_cont
          << std::setw(12) << iface.delta_g_int
          << std::endl;
   }

   oss << std::endl;
   oss << "Monomer properties:" << std::endl;
   oss << std::setw(15) << "Monomer"
       << std::setw(10) << "Mass(Da)"
       << std::setw(10) << "ASA(A^2)"
       << std::setw(8)  << "Atoms"
       << std::setw(8)  << "Res"
       << std::endl;
   oss << std::string(51, '-') << std::endl;

   for (const auto &pair : monomer_props) {
      if (pair.first.is_asu()) {
         oss << std::setw(15) << pair.first.label()
             << std::setw(10) << std::setprecision(0) << pair.second.mass
             << std::setw(10) << std::setprecision(1) << pair.second.surface_area
             << std::setw(8)  << pair.second.n_atoms
             << std::setw(8)  << pair.second.n_residues
             << std::endl;
      }
   }

   return oss.str();
}
