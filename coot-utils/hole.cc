/* coot-utils/hole.cc
 * 
 * Copyright 2011, 2012 by The University of Oxford
 * Copyright 2015, 2016 by Medical Research Council
 * Author: Paul Emsley
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
#include "clipper/ccp4/ccp4_map_io.h"
// #include "clipper/contrib/skeleton.h" // needed?

#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-hole.hh"



coot::hole::hole(mmdb::Manager *mol_in, 
                 const clipper::Coord_orth &from_pt_in,
                 const clipper::Coord_orth &to_pt_in,
                 const coot::protein_geometry &geom) {

   mol = mol_in;
   from_pt = from_pt_in;
   to_pt = to_pt_in;
   assign_vdw_radii(geom); 
   colour_map_multiplier = 1.0;
   colour_map_offset = 0.0;
}

void
coot::hole::assign_vdw_radii(const coot::protein_geometry &geom) {

   int imol = 0; // dummy
   bool use_vdwH_flag = 0; // extended atoms

   std::map<std::pair<std::string, std::string>, double> cached_radii;
   std::map<std::pair<std::string, std::string>, double>::const_iterator it;

   double radius; 
   
   radius_handle = mol->RegisterUDReal(mmdb::UDR_ATOM, "atom_radius");
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         mmdb::Residue *residue_p;
         mmdb::Atom *at;
         for (int ires=0; ires<nres; ires++) { 
            residue_p = chain_p->GetResidue(ires);
            int n_atoms = residue_p->GetNumberOfAtoms();
            std::string residue_name = residue_p->GetResName();
            for (int iat=0; iat<n_atoms; iat++) {
               at = residue_p->GetAtom(iat);

               std::string atom_name = at->name;

               // try cache first
               std::pair<std::string, std::string> p(atom_name, residue_name);

               if (cached_radii.find(p) != cached_radii.end()) {
                  radius = it->second;
               } else {
                  radius = geom.get_vdw_radius(atom_name, residue_name, imol, use_vdwH_flag);
               }
               if (radius > 0) {
                  at->PutUDData(radius_handle, radius);
               } else {
                  std::string ele = at->element;
                  // make a reasonable default
                  mmdb::realtype radius = 1.7;
                   if (ele == " N")
                      radius = 1.55;
                   if (ele == " 0")
                      radius = 1.52;
                   if (ele == " H")
                      radius = 1.2;
                  at->PutUDData(radius_handle, radius); 
               }
            }
         }
      }
   }
}

void
coot::hole::debug_atom_radii() const {

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         mmdb::Residue *residue_p;
         mmdb::Atom *at;
         for (int ires=0; ires<nres; ires++) { 
            residue_p = chain_p->GetResidue(ires);
            int n_atoms = residue_p->GetNumberOfAtoms();
            std::string residue_name = residue_p->GetResName();
            for (int iat=0; iat<n_atoms; iat++) {
               at = residue_p->GetAtom(iat);
               mmdb::realtype radius;
               at->GetUDData(radius_handle, radius);
               std::cout << "   " << coot::atom_spec_t(at) << " with radius " << radius<< std::endl;
            }
         }
      }
   }
} 

void 
coot::hole::make_atom_selection(int selhnd, const clipper::Coord_orth &pt,
                                double radius_prev) const {

   // what protein atoms are within radius_prev of pt?

   mol->SelectAtoms(selhnd, 0, "*",
                    mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                    "!HOH", // residue names
                    "*", // atom names
                    "H,C,O,N,S", // elements
                    "*", "*", "*", // altlocs, segments, charges
                    -1, -1, // occupancy (no limits)
                    pt.x(), pt.y(), pt.z(), radius_prev, mmdb::SKEY_NEW);

   if (0) {
      mmdb::PPAtom atom_selection = 0;
      int n_selected_atoms;
      mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
      std::cout << "in make_atom_selection() selected " << n_selected_atoms
                << " atoms in sphere selection around  " << pt.format() << std::endl;
   }
}


std::pair<std::vector<std::pair<clipper::Coord_orth, double> >, std::vector<coot::hole_surface_point_t> >
coot::hole::generate() {

   int nsteps = 200;
   double radius_prev = 7; // max prob size + 1.7 (max radius) + a bit.

   // position and sphere radius
   // 
   std::vector<std::pair<clipper::Coord_orth, double> > probe_path;
   
   clipper::Coord_orth diff = to_pt - from_pt;
   v_hat = clipper::Coord_orth(diff.unit());

   clipper::Coord_orth prev_point = from_pt; // initally.
   clipper::Coord_orth frag = 1/double(nsteps) * diff;
   
   for (int istep=0; istep<=nsteps; istep++) {
      clipper::Coord_orth l(double(istep) * frag);
      clipper::Coord_orth pt_linear(from_pt + l);

      // We are considering where to place the start point for the new
      // level (pt).  It should be
      // 1) just down a bit (diff) from the previous level
      // and/or
      // 2) a step along the vector from_pt -> to_pt.
      //
      // If/when the point previous new_pt has wandered off into the
      // unknown (because it was so high that it managed to escape
      // being in the channel) then what do we do? - because using 1)
      // from there would mean that we would miss then channel
      // completely.
      //
      // But following the 2) rule would mean that we would not be
      // able to follow kinky paths and that the precision of the
      // (start and) finish points would be important.
      //
      // So let's choose rule 1 except where that means we have gone
      // way off (5A?)
      // 
      clipper::Coord_orth pt = prev_point + frag;
      double dist_diff = clipper::Coord_orth::length(pt_linear, pt);

      if (0)
         std::cout << "dist_diff: " << dist_diff << " " << pt_linear.format()
                   << " " << pt.format() << " frag: " << frag.format() << " l:"
                   << l.format() << " i_step: " << istep << std::endl;

      if (dist_diff > 7)
         pt = pt_linear;
      

      int selhnd = mol->NewSelection();
      make_atom_selection(selhnd, pt, radius_prev);

      // if that didn't select atoms, then we increase the radius, try this for up to 10 rounds:
      int n_rounds = 10;
      int i_round = 0;
      mmdb::PPAtom atom_selection = 0;
      int n_selected_atoms;
      mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
      while ((n_selected_atoms == 0) && (i_round < n_rounds)) { 
         mol->DeleteSelection(selhnd);// not sure if needed
         selhnd = mol->NewSelection();
         make_atom_selection(selhnd, pt, radius_prev*(i_round+1+0.2)*(i_round+1+0.2));
         i_round++;
      }

      
      std::pair<clipper::Coord_orth, double> new_pt = optimize_point(pt, selhnd);
      probe_path.push_back(new_pt);
      if (false) { 
         double d = sqrt(l.lengthsq());
         std::cout << "istep: " << istep << " l: " << d
                   << " ss: " << new_pt.second << "        "
                   << new_pt.first.x() << " "
                   << new_pt.first.y() << " "
                   << new_pt.first.z() << " "
                   << std::endl;
      }
      mol->DeleteSelection(selhnd);
      // next round
      prev_point = new_pt.first;
   }

   std::vector<coot::hole_surface_point_t> surface_points = get_surface_points(probe_path);
   // write_probe_path_using_spheres(surface_points, "probe-surface-points");

   return std::pair<std::vector<std::pair<clipper::Coord_orth, double> >, std::vector<coot::hole_surface_point_t> > (probe_path, surface_points);
}

// By moving about in the plane we try to maximise the distance
// between the point and the atoms around it.
// 
// Return the position and the max dist.
// 
std::pair<clipper::Coord_orth, double>
coot::hole::optimize_point(const clipper::Coord_orth &pt, int selhnd) {

   int max_moves_without_improvement = 80;
   int move_count = 0; // increment if we try to make a move that does not improve
   double biggest_allowed_distance_to_protein = 5.0;
   
   mmdb::PPAtom atom_selection = 0;
   int n_selected_atoms;
   mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
   // std::cout << "optimize_point(): n_selected_atoms: " << n_selected_atoms
   // << std::endl;
   
   clipper::Coord_orth current_pt = pt;

   double D_max = 0.1;
   float rmi = 1.0/float(RAND_MAX);

   double current_ss = sphere_size(pt, selhnd);
   while (move_count < max_moves_without_improvement) {
      double d1 = 2.0 * coot::util::random() * rmi -1.0;
      double d2 = 2.0 * coot::util::random() * rmi -1.0;
      double d3 = 2.0 * coot::util::random() * rmi -1.0;
      clipper::Coord_orth y_rand(d1, d2, d3);
      clipper::Coord_orth y_hat_rand(y_rand.unit());
      clipper::Coord_orth y_prime_rand =
         y_hat_rand - (clipper::Coord_orth::dot(v_hat, y_hat_rand)) * v_hat;
      
      clipper::Coord_orth trial_pt = current_pt + D_max * y_prime_rand;

      // sphere_size() can return negative (no atoms selected).
      // 
      double ss = sphere_size(trial_pt, selhnd);
      if (ss > current_ss) {
            current_pt = trial_pt;
            current_ss = ss;
//             std::cout << "   updated current_ss to " << current_ss
//                       << " for start point " << pt.format() << " and trial point "
//                       << trial_pt.format() << "\n";
            move_count = 0; // reset
            if (ss > biggest_allowed_distance_to_protein) {
               // std::cout << "breaking..." << std::endl;
               break; // that's enough, the point has escaped presumably
            }
      } else {
         move_count++;
      }
   }
   return std::pair<clipper::Coord_orth, double> (current_pt, current_ss);
}


double
coot::hole::sphere_size(const clipper::Coord_orth &pt, int selhnd) const {

   double r = -1;
   mmdb::PPAtom atom_selection = 0;
   int n_selected_atoms;
   mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

   // std::cout << "sphere_size: n_selected_atoms: " << n_selected_atoms
   // << std::endl;
   
   double largest_possible_sphere = 99999;
   bool was_set = 0;
   mmdb::realtype atom_vdw_radius;

   for (int iat=0; iat<n_selected_atoms; iat++) {
      clipper::Coord_orth atom_pos(atom_selection[iat]->x,
                                   atom_selection[iat]->y,
                                   atom_selection[iat]->z);
      double r_1 = clipper::Coord_orth::length(atom_pos, pt);
      atom_selection[iat]->GetUDData(radius_handle, atom_vdw_radius);
      double r = r_1 - atom_vdw_radius;
      if (r < largest_possible_sphere) {
         largest_possible_sphere = r;
         if (! was_set)
            was_set = 1;
      }
   }
   
   if (was_set)
      return largest_possible_sphere;
   else
      return -1; // fail
}

// point on the wall of the probe path.
// 
void
coot::hole::write_probe_path(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path) const {

   std::string file_name = "probe-surface-points";

   unsigned int n_theta_points = 20;
   double inv_n_theta = 2 * M_PI/double(n_theta_points);

   std::ofstream render_stream(file_name.c_str());

   if (!render_stream) {
      std::cout << "failed to open " << file_name << std::endl;
   } else { 
  
      clipper::Coord_orth arb(1.1, 1.2, 1.3);
      clipper::Coord_orth plane_vect(clipper::Coord_orth::cross(v_hat, arb));
      clipper::Coord_orth unit_plane_vect(plane_vect.unit());

      for (unsigned int i=0; i<probe_path.size(); i++) {

         unsigned int max_upstream_check_limit = i-5;
         unsigned int max_downstream_check_limit = i+5;

         if (max_downstream_check_limit>=probe_path.size())
            max_downstream_check_limit=probe_path.size()-1;

         render_stream << probe_path[i].first.x() << " "
                       << probe_path[i].first.y() << " "
                       << probe_path[i].first.z() << " \"red\"\n";
         
         std::string colour = "blue";

         if (probe_path[i].second < 2.1)
            colour = "blue";
         if (probe_path[i].second < 1.9)
            colour = "cyan";
         if (probe_path[i].second < 1.7)
            colour = "green";
         if (probe_path[i].second < 1.5)
            colour = "greentint";
         if (probe_path[i].second < 1.3)
            colour = "sea";
         if (probe_path[i].second < 1.1)
            colour = "yellow";
         if (probe_path[i].second < 0.9)
            colour = "yelllowtint";
         if (probe_path[i].second < 0.7)
            colour = "orange";
         if (probe_path[i].second < 0.5)
            colour = "red";
         if (probe_path[i].second < 0.3)
            colour = "hotpink";
         

         for (unsigned int itheta=0; itheta<n_theta_points; itheta++) {
            double theta = inv_n_theta * double(itheta);
            if (coot::util::even_p(i))
               theta += inv_n_theta * 0.5;

            // now let's rotate unit_plane_vect around the circle.
            clipper::Coord_orth pos(probe_path[i].second * unit_plane_vect);
            clipper::Coord_orth circle_point =
               coot::util::rotate_around_vector(v_hat,
                                                pos,
                                                clipper::Coord_orth(0,0,0),
                                                theta);
            clipper::Coord_orth surface_point = probe_path[i].first + circle_point;

            render_stream << surface_point.x() << " "
                          << surface_point.y() << " "
                          << surface_point.z() << " "
                          << "\"" << colour << "\""
                          << "\n";
         }
      }
   }

}

void
coot::hole::write_probe_path_using_spheres(const std::vector<coot::hole_surface_point_t> &surface_points,
                                           const std::string &file_name) const {

   std::ofstream render_stream(file_name.c_str());
   if (!render_stream) {
      std::cout << "failed to open " << file_name << std::endl;
   } else {
      for (unsigned int i=0; i<surface_points.size(); i++) { 
         render_stream << surface_points[i].position.format() << " "
                       << surface_points[i].normal.format() << " "
                       << surface_points[i].colour << "\n";
      }
   } 
}


std::vector<coot::hole_surface_point_t> 
coot::hole::get_surface_points(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path) const {

   std::vector<coot::hole_surface_point_t> surface_points; // return this.
   
   double dot_density = 0.35;
   // 
   double phi_step = 5.0 * (M_PI/180.0);
   double theta_step = 5.0 * (M_PI/180.0);

   for (unsigned int i=0; i<probe_path.size(); i++) {

      // int because we max_upstream_check_limit might go negative (briefly).
      // 
      int max_upstream_check_limit = i-10;
      int max_downstream_check_limit = i+10;
      const clipper::Coord_orth &this_centre = probe_path[i].first;
      const double &r = probe_path[i].second;

      if (max_upstream_check_limit<0)
         max_upstream_check_limit = 0;
      if (max_downstream_check_limit>=int(probe_path.size()))
         max_downstream_check_limit=probe_path.size()-1;

      coot::colour_holder vertex_colour = probe_size_to_colour(probe_path[i].second);
      // std::cout << "         " << vertex_colour << std::endl;
      bool even = 1;

      // the begin and end point we don't want to draw the tops and bottoms of the "bulbs".
      // Leave them open.
      //
      bool check_ends = 0;
      double dot_prod_multiplier = 1;
      if (i == 0) {
         check_ends = 1;
      }
      if (i == (probe_path.size() -1)) {
         check_ends = 1;
         dot_prod_multiplier = -1;
      }
            
      for (double theta=0; theta<M_PI; theta+=theta_step) {
         double phi_step_inner = phi_step + 0.1 * pow(theta-0.5*M_PI, 2);
         double sin_theta = sin(theta);
         double cos_theta = cos(theta);
         for (double phi=0; phi<2*M_PI; phi+=phi_step_inner) {
            if (even) { // even
               clipper::Coord_orth origin_based_sphere_point(r*cos(phi)*sin_theta,
                                                             r*sin(phi)*sin_theta,
                                                             r*cos_theta);
               clipper::Coord_orth surface_point = this_centre + origin_based_sphere_point;
               clipper::Coord_orth normal(origin_based_sphere_point.unit());
               bool reject = 0;
               for (int iprev=max_upstream_check_limit; iprev<=max_downstream_check_limit; iprev++) {
                  if (iprev != int(i)) {

                     clipper::Coord_orth diff = probe_path[iprev].first - surface_point;
                     double d_sqrd = diff.lengthsq();
                     double probe_prev_sqrd = probe_path[iprev].second * probe_path[iprev].second;
                     if (d_sqrd < probe_prev_sqrd) {
                        reject = 1;
                        break;
                     }
                  }
               }
               if (! reject) {

                  bool do_it = 1;
                     
                  // turn off do_it so that we have open cups at the end points.
                  // 
                  if (check_ends == 1) {
                     double dp = clipper::Coord_orth::dot(v_hat, origin_based_sphere_point);
                     if (dp * dot_prod_multiplier < 0)
                        do_it = 0;
                  }
                     
                  if (do_it) { 
                     coot::hole_surface_point_t sp(surface_point, normal, vertex_colour);
                     surface_points.push_back(sp);
                  }
               }
            }
            // next round
            even = 1 - even;
         }
      }
   }
   return surface_points;
} 



// for generating a map, return (min_x, min_y_min_z), (max_x,
// max_y, max_z) for the points in the path
// 
std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::hole::get_min_and_max(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path) const {

   clipper::Coord_orth min_pos(0,0,0);
   clipper::Coord_orth max_pos(0,0,0);
   
   std::pair<clipper::Coord_orth, clipper::Coord_orth> p(min_pos, max_pos);
   if (probe_path.size() > 0) {
      min_pos = clipper::Coord_orth( 1e20,  1e20,  1e20);
      max_pos = clipper::Coord_orth(-1e20, -1e20, -1e20);
      for (unsigned int i=0; i<probe_path.size(); i++) {
         const clipper::Coord_orth &pt = probe_path[i].first;
         if (pt[0] < min_pos[0]) min_pos[0] = pt[0];
         if (pt[1] < min_pos[1]) min_pos[1] = pt[1];
         if (pt[2] < min_pos[2]) min_pos[2] = pt[2];
         if (pt[0] > max_pos[0]) max_pos[0] = pt[0];
         if (pt[1] > max_pos[1]) max_pos[1] = pt[1];
         if (pt[2] > max_pos[2]) max_pos[2] = pt[2];
      }
      p = std::pair<clipper::Coord_orth, clipper::Coord_orth>(min_pos, max_pos);
   } 
   return p;
} 

clipper::Xmap<float>
coot::hole::carve_a_map(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path,
                        const clipper::Xmap<float> &xmap_ref,
                        const std::string &file_name) const {

   // get min_max coords of probe path
   //
   std::pair<clipper::Coord_orth, clipper::Coord_orth> min_max = get_min_and_max(probe_path);
   clipper::Coord_orth middle(0.5*(min_max.first[0]+ min_max.second[0]),
                              0.5*(min_max.first[1]+ min_max.second[1]),
                              0.5*(min_max.first[2]+ min_max.second[2]));

   float border = 50;
   clipper::Cell_descr cell_descr(min_max.second[0]-min_max.first[0] + border,
                                  min_max.second[1]-min_max.first[1] + border,
                                  min_max.second[2]-min_max.first[2] + border,
                                  M_PI_2, M_PI_2, M_PI_2);
   clipper::Cell cell(cell_descr);

   clipper::ftype sampling = 2.0;
   clipper::Resolution reso(2.0);
   clipper::Grid_sampling grid(clipper::Spacegroup::p1(), cell, reso, sampling);

   float radius = clipper::Coord_orth::length(middle, min_max.first);

   // get grid range
   // gr0: a grid range of the correct size (at the origin)
   // gr1: a grid range of the correct size (around the correct place, comg)
   clipper::Grid_range gr0(cell, grid, radius);
   clipper::Grid_range gr1(gr0.min() + middle.coord_frac(cell).coord_grid(grid),
                           gr0.max() + middle.coord_frac(cell).coord_grid(grid));

   std::cout << "Here with cell "           << cell.format() << std::endl;
   std::cout << "Here with gr1 "            << gr1.format() << std::endl;
   std::cout << "Here with middle "         << middle.format() << std::endl;
   std::cout << "Here with min_max.first "  << min_max.first.format() << std::endl;
   std::cout << "Here with min_max.second " << min_max.second.format() << std::endl;

   clipper::Xmap<float> xmap = xmap_ref;
   
   // clipper::Coord_grid offset =
   // xmap.coord_map(nxmap.coord_orth(clipper::Coord_map(0.0,0.0,0.0))).coord_grid();

   std::cout << "put stuff in nxmap " << std::endl;
   typedef clipper::NXmap<float>::Map_reference_index NRI;
   // for (NRI inx = nxmap.first(); !inx.last(); inx.next()) { nxmap[inx] = 1; }
   
   // for (NRI inx = nxmap.first(); !inx.last(); inx.next()) { }
   
   // clipper::Skeleton_basic::Neighbours neighb(nxmap);
   
   for (unsigned int i=0; i<probe_path.size(); i++) { 
      const clipper::Coord_orth &pt = probe_path[i].first;
      const double &r = probe_path[i].second;
      clipper::Coord_frac cf = pt.coord_frac(xmap_ref.cell());
      clipper::Coord_map  cm = cf.coord_map(xmap_ref.grid_sampling());
      clipper::Coord_grid cg = cm.coord_grid();
      // xmap.set_data(cg, 2);
      mask_around_coord(pt, r, &xmap);
      // std::cout << "info::" << cg.format() << " set to 2"  << std::endl;
   }

   clipper::CCP4MAPfile mapout;
   mapout.open_write(file_name);
   mapout.set_cell(cell);
   mapout.export_xmap(xmap);
   mapout.close_write();
   std::cout << "wrote map " << file_name << std::endl;

   return xmap;
}

void
coot::hole::mask_around_coord(const clipper::Coord_orth &co, float atom_radius,
                              clipper::Xmap<float> *xmap) const {
   
   clipper::Coord_frac cf = co.coord_frac(xmap->cell());

   clipper::Coord_frac box0(
                            cf.u() - atom_radius/xmap->cell().descr().a(),
                            cf.v() - atom_radius/xmap->cell().descr().b(),
                            cf.w() - atom_radius/xmap->cell().descr().c());

   clipper::Coord_frac box1(
                            cf.u() + atom_radius/xmap->cell().descr().a(),
                            cf.v() + atom_radius/xmap->cell().descr().b(),
                            cf.w() + atom_radius/xmap->cell().descr().c());

   clipper::Grid_map grid(box0.coord_grid(xmap->grid_sampling()),
                          box1.coord_grid(xmap->grid_sampling()));

   float atom_radius_sq = atom_radius * atom_radius;
   int nhit = 0;
   int nmiss = 0;

   clipper::Xmap_base::Map_reference_coord ix( *xmap, grid.min() ), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) { 
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) { 
         for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
            if ( (iw.coord().coord_frac(xmap->grid_sampling()).coord_orth(xmap->cell()) - co).lengthsq() < atom_radius_sq) {

               float masked_map_val = 4.9;

               if (false)
                  std::cout << "masked " << masked_map_val << " point at " 
                            << iw.coord().coord_frac(xmap->grid_sampling()).coord_orth(xmap->cell()).format()
                            << " centre point: " << co.format() << " " 
                            << (iw.coord().coord_frac(xmap->grid_sampling()).coord_orth(xmap->cell()) - co).lengthsq()
                            << "\n";
               
               (*xmap)[iw] = masked_map_val;
               nhit++;
            } else {
               nmiss++;
            }
         }
      }
   }
   std::cout << "nhit " << nhit << " nmiss " << nmiss << std::endl;
}



clipper::NXmap<float>
coot::hole::carve_a_map(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path,
                        const std::string &file_name) const {

   std::pair<clipper::Coord_orth, clipper::Coord_orth> min_max = get_min_and_max(probe_path);
   clipper::Coord_orth middle(0.5*(min_max.first[0]+ min_max.second[0]),
                              0.5*(min_max.first[1]+ min_max.second[1]),
                              0.5*(min_max.first[2]+ min_max.second[2]));

   float border = 200;
   clipper::Cell_descr cell_descr(min_max.second[0]-min_max.first[0] + border,
                                  min_max.second[1]-min_max.first[1] + border,
                                  min_max.second[2]-min_max.first[2] + border,
                                  M_PI_2, M_PI_2, M_PI_2);
   clipper::Cell cell(cell_descr);

   clipper::ftype sampling = 2.0;
   clipper::Resolution reso(2.0);
   clipper::Grid_sampling grid(clipper::Spacegroup::p1(), cell, reso, sampling);

   float radius = clipper::Coord_orth::length(middle, min_max.first);

   // get grid range
   // gr0: a grid range of the correct size (at the origin)
   // gr1: a grid range of the correct size (around the correct place, comg)
   clipper::Grid_range gr0(cell, grid, radius);
   clipper::Grid_range gr1(gr0.min() + middle.coord_frac(cell).coord_grid(grid),
                           gr0.max() + middle.coord_frac(cell).coord_grid(grid));

   std::cout << "Here with NX min_max.first  "  << min_max.first.format() << std::endl;
   std::cout << "Here with NX min_max.second "  << min_max.second.format() << std::endl;
   std::cout << "Here with NX middle "           << middle.format() << std::endl;
   std::cout << "max - min                   = " << (min_max.second - min_max.first).format()
             << std::endl;
   std::cout << "Here with NX cell           = " << cell.format() << std::endl;
   std::cout << "Here with NX gr0             min: "  << gr0.min().format() << " max: "
             << gr0.max().format() << std::endl;
   std::cout << "Here with NX gr1             min: "  << gr1.min().format() << " max: "
             << gr1.max().format() << std::endl;

   // init nxmap
   clipper::NXmap<float> nxmap(cell, grid, gr1);
   
   std::cout << "created nxmap " << std::endl;

   clipper::Coord_grid offset = -min_max.first.coord_frac(cell).coord_map(grid).coord_grid();

   // debug
   clipper::Coord_grid offset_orig_grid = offset;
   clipper::Coord_map  offset_orig_map  = offset_orig_grid.coord_map();
   clipper::Coord_frac offset_orig_frac = offset_orig_map.coord_frac(grid);
   clipper::Coord_orth offset_orig_orth = offset_orig_frac.coord_orth(cell);
   std::cout << "Here with offset orig grid " << offset_orig_grid.format() << std::endl;
   std::cout << "Here with offset orig frac " << offset_orig_frac.format() << std::endl;
   std::cout << "Here with offset orig orth " << offset_orig_orth.format() << std::endl;

   // clipper::Coord_grid offset_grid(12, 23, 4); // good for 50
   // clipper::Coord_grid offset_grid(11, 22, 3); // pretty good for 40
   // clipper::Coord_grid offset_grid(11, 23, 3); // good for 100 (needs 11.5 for x)
   // clipper::Coord_grid offset_grid(11, 23, 3);    // good for 200.

   // try with molecule centred around 100,100,100
   clipper::Coord_grid offset_grid(19, 36, 3);
                                                  // So offset_grid is practically indepdendent
                                                  // of border.  I guess though, that
                                                  // it depends on resolution and sampling.
   
   clipper::Coord_map  offset_map  = offset_grid.coord_map();
   clipper::Coord_frac offset_frac = offset_map.coord_frac(grid);
   clipper::Coord_orth offset_orth = offset_frac.coord_orth(cell);

   std::cout << "Here with extra offset grid " << offset_grid.format() << std::endl;
   std::cout << "Here with extra offset frac " << offset_frac.format() << std::endl;
   std::cout << "Here with extra offset orth " << offset_orth.format() << std::endl;

   std::cout << "put stuff in nxmap " << std::endl;
   typedef clipper::NXmap<float>::Map_reference_index NRI;
   for (NRI inx = nxmap.first(); !inx.last(); inx.next()) { nxmap[inx] = 1; }
   
   offset += offset_grid;

   for (unsigned int i=0; i<probe_path.size(); i++) { 
      const clipper::Coord_orth &pt = probe_path[i].first;
      clipper::Coord_frac cf = pt.coord_frac(cell);
      clipper::Coord_map  cm = cf.coord_map(grid);
      clipper::Coord_grid cg = cm.coord_grid();
      if (false)
         std::cout << "info:: at " << cg.format() << " + " << offset.format()
                   << " setting to 0"  << std::endl;
      nxmap.set_data(cg+offset, 0);
   }

   clipper::CCP4MAPfile mapout;
   mapout.open_write(file_name);
   mapout.set_cell(cell);
   mapout.export_nxmap(nxmap);
   mapout.close_write();
   std::cout << "wrote map " << file_name << std::endl;

   return nxmap;

} 
