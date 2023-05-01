/* ligand/coot-map-heavy.cc
 *
 * Copyright 2004, 2005, 2006 The University of York
 * Copyright 2013 by Medical Research Council
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
 * 02110-1301, USA.
 */


#include <thread>
#include <iomanip>

#include "clipper/core/map_interp.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/rotation.h"

#include "utils/coot-utils.hh"
#include "coot-map-heavy.hh"
#include "coot-map-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh" // torsionable-bonds

#ifdef HAVE_GSL
// return status success (something moved = 1) or error/failure (0)
//
int
coot::util::fit_to_map_by_simplex_rigid(mmdb::PPAtom atom_selection,
					int n_selected_atoms,
					const clipper::Xmap<float> &xmap) {

   int i_r = 0;

   const gsl_multimin_fminimizer_type *T =
      gsl_multimin_fminimizer_nmsimplex;
   gsl_multimin_fminimizer *s = NULL;
   gsl_vector *ss, *x;
   gsl_multimin_function minex_func;

   size_t iter = 0;
   int rval = GSL_CONTINUE;
   int status = GSL_SUCCESS;
   double ssval;
   coot::util::simplex_param_t par;

   // Set the par parameters:
   par.n_atoms = n_selected_atoms;
   par.orig_atoms = atom_selection;
   clipper::Coord_orth co(0.0, 0.0, 0.0);
   for (int i=0; i<n_selected_atoms; i++)
      co += clipper::Coord_orth(atom_selection[i]->x,
				atom_selection[i]->y,
				atom_selection[i]->z);
   co = 1/float(n_selected_atoms) * co;
   par.atoms_centre = co;
   par.xmap = &xmap;


   int np = 3 * n_selected_atoms;

   ss = gsl_vector_alloc(np);

   if (ss == NULL) {
      GSL_ERROR_VAL("failed to allocate space for ss", GSL_ENOMEM, 0);
   }

   x = gsl_vector_alloc(3*n_selected_atoms);
   gsl_vector_set_all (ss, 1.0);
   gsl_vector_set_all (x,  0.01); // no rotations or translations initially.

   // setup_simplex_x_internal(x, atom_selection, n_selected_atoms);

   minex_func.f = &coot::util::my_f_simplex_rigid_internal;
   minex_func.n = np;
   minex_func.params = (void *)&par;

   s = gsl_multimin_fminimizer_alloc (T, np);
   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

   while (rval == GSL_CONTINUE) {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
	 break;
      rval = gsl_multimin_test_size (gsl_multimin_fminimizer_size(s),
				     1e-3);
      ssval = gsl_multimin_fminimizer_size(s);

      if (rval == GSL_SUCCESS) {
	 std::cout << "converged at minimum\n";
	 i_r = 1;
	 simplex_apply_shifts_rigid_internal(s->x, par);
      }

//       std::cout << iter << " "
// 		<< gsl_vector_get(s->x, 0) << " "
// 		<< gsl_vector_get(s->x, 1) << " "
// 		<< gsl_vector_get(s->x, 2) << " "
// 		<< gsl_vector_get(s->x, 3) << " "
// 		<< gsl_vector_get(s->x, 4) << " "
// 		<< gsl_vector_get(s->x, 5) << " ";
//       std::cout << "f() " << s->fval << " ssize " << ssval << "\n";
   }


   gsl_vector_free(x);
   gsl_vector_free(ss);
   gsl_multimin_fminimizer_free(s);

   return i_r;
}

// internal simplex setup function:
void
coot::util::setup_simplex_x_internal_rigid(gsl_vector *x,
					   mmdb::PPAtom atom_selection,
					   int n_selected_atoms) {

   // unnecessary, I think.  All rots and trans angles are set to 0.
}

// Simplex Refinement is finished, force-feed the results (new
// coordinates) back into the atom_selection.
//
void
coot::util::simplex_apply_shifts_rigid_internal(gsl_vector *s,
						coot::util::simplex_param_t &par) {


   double sin_t;
   double cos_t;

   sin_t = sin(-clipper::Util::d2rad(gsl_vector_get(s, 3)));
   cos_t = cos(-clipper::Util::d2rad(gsl_vector_get(s, 3)));
   clipper::Mat33<double> x_mat(1,0,0, 0,cos_t,-sin_t, 0,sin_t,cos_t);

   sin_t = sin(-clipper::Util::d2rad(gsl_vector_get(s, 4)));
   cos_t = cos(-clipper::Util::d2rad(gsl_vector_get(s, 4)));
   clipper::Mat33<double> y_mat(cos_t,0,sin_t, 0,1,0, -sin_t,0,cos_t);

   sin_t = sin(-clipper::Util::d2rad(gsl_vector_get(s, 5)));
   cos_t = cos(-clipper::Util::d2rad(gsl_vector_get(s, 5)));
   clipper::Mat33<double> z_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);

   clipper::Mat33<double> angle_mat = x_mat * y_mat * z_mat;

   clipper::Coord_orth trans(gsl_vector_get(s, 0),
			     gsl_vector_get(s, 1),
			     gsl_vector_get(s, 2));

   clipper::RTop_orth rtop(angle_mat, trans);
   clipper::Coord_orth point;

   for (int i=0; i<par.n_atoms; i++) {

      clipper::Coord_orth orig_p(par.orig_atoms[i]->x,
				 par.orig_atoms[i]->y,
				 par.orig_atoms[i]->z);

      clipper::Coord_orth point = orig_p.transform(rtop);
      point = par.atoms_centre + (orig_p - par.atoms_centre).transform(rtop);

      par.orig_atoms[i]->x = point.x();
      par.orig_atoms[i]->y = point.y();
      par.orig_atoms[i]->z = point.z();
   }
}


// The function that returns the value:
//
// Remember that v are the rotation and translation variables,
// which are applied to the (original) coordinates (in params)
//
double
coot::util::my_f_simplex_rigid_internal (const gsl_vector *v,
					 void *params) {

   coot::util::simplex_param_t *p = (coot::util::simplex_param_t *) params;

   double score = 0;

   double sin_t;
   double cos_t;

   sin_t = sin(-clipper::Util::d2rad(gsl_vector_get(v, 3)));
   cos_t = cos(-clipper::Util::d2rad(gsl_vector_get(v, 3)));
   clipper::Mat33<double> x_mat(1,0,0, 0,cos_t,-sin_t, 0,sin_t,cos_t);

   sin_t = sin(-clipper::Util::d2rad(gsl_vector_get(v, 4)));
   cos_t = cos(-clipper::Util::d2rad(gsl_vector_get(v, 4)));
   clipper::Mat33<double> y_mat(cos_t,0,sin_t, 0,1,0, -sin_t,0,cos_t);

   sin_t = sin(-clipper::Util::d2rad(gsl_vector_get(v, 5)));
   cos_t = cos(-clipper::Util::d2rad(gsl_vector_get(v, 5)));
   clipper::Mat33<double> z_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);

   clipper::Mat33<double> angle_mat = x_mat * y_mat * z_mat;

   clipper::Coord_orth trans(gsl_vector_get(v, 0),
			     gsl_vector_get(v, 1),
			     gsl_vector_get(v, 2));

   clipper::RTop_orth rtop(angle_mat, trans);
   clipper::Coord_orth point;

   for (int i=0; i<p->n_atoms; i++) {

      clipper::Coord_orth orig_p(p->orig_atoms[i]->x,
				 p->orig_atoms[i]->y,
				 p->orig_atoms[i]->z);
      point = p->atoms_centre + (orig_p - p->atoms_centre).transform(rtop);

      // we are trying to minimize, don't forget:
      score -= density_at_point(*(p->xmap), point);
   }

   return score;
}

#endif // HAVE_GSL

float
coot::util::z_weighted_density_at_point(const clipper::Coord_orth &pt,
					const std::string &ele,
					const std::vector<std::pair<std::string, int> > &atom_number_list,
					const clipper::Xmap<float> &map_in) {

   float d = coot::util::density_at_point(map_in, pt);
   float z = coot::util::atomic_number(ele, atom_number_list);
   if (z< 0.0)
      z = 6; // carbon, say
   return d*z;
}

float
coot::util::z_weighted_density_at_point_linear_interp(const clipper::Coord_orth &pt,
						      const std::string &ele,
						      const std::vector<std::pair<std::string, int> > &atom_number_list,
						      const clipper::Xmap<float> &map_in) {

   float d = coot::util::density_at_point_by_linear_interpolation(map_in, pt);
   float z = coot::util::atomic_number(ele, atom_number_list);
   if (z< 0.0)
      z = 6; // carbon, say
   return d*z;
}

float
coot::util::z_weighted_density_at_nearest_grid(const clipper::Coord_orth &pt,
                                               const std::string &ele,
                                               const std::vector<std::pair<std::string, int> > &atom_number_list,
                                               const clipper::Xmap<float> &map_in) {

   float d = density_at_point_by_nearest_grid(map_in, pt);
   float z = atomic_number(ele, atom_number_list);
   if (z< 0.0)
      z = 6; // carbon, say
   return d*z;

}

float
coot::util::z_weighted_density_score(const std::vector<mmdb::Atom *> &atoms,
				     const std::vector<std::pair<std::string, int> > &atom_number_list,
				     const clipper::Xmap<float> &map) {
   float sum_d = 0;
   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      clipper::Coord_orth co(atoms[iat]->x, atoms[iat]->y, atoms[iat]->z);
      float d = z_weighted_density_at_point(co, atoms[iat]->element, atom_number_list, map);
      sum_d += d;
   }
   return  sum_d;
}

float
coot::util::z_weighted_density_score(const minimol::molecule &mol,
				     const std::vector<std::pair<std::string, int> > &atom_number_list,
				     const clipper::Xmap<float> &map) {
   float sum_d = 0;
   std::vector<coot::minimol::atom *> atoms = mol.select_atoms_serial();
   for (unsigned int i=0; i<atoms.size(); i++) {
      float d = z_weighted_density_at_point(atoms[i]->pos, atoms[i]->element, atom_number_list, map);
      sum_d += d;
   }
   return sum_d;
}

float
coot::util::z_weighted_density_score_linear_interp(const minimol::molecule &mol,
						   const std::vector<std::pair<std::string, int> > &atom_number_list,
						   const clipper::Xmap<float> &map) {
   float sum_d = 0;
   std::vector<coot::minimol::atom *> atoms = mol.select_atoms_serial();
   for (unsigned int i=0; i<atoms.size(); i++) {
      float d = z_weighted_density_at_point_linear_interp(atoms[i]->pos, atoms[i]->element, atom_number_list, map);
      sum_d += d;
   }
   return sum_d;
}

float
coot::util::z_weighted_density_score_nearest(const minimol::molecule &mol,
                                             const std::vector<std::pair<std::string, int> > &atom_number_list,
                                             const clipper::Xmap<float> &map) {

   float sum_d = 0;
   std::vector<coot::minimol::atom *> atoms = mol.select_atoms_serial();
   for (unsigned int i=0; i<atoms.size(); i++) {
      float d = z_weighted_density_at_nearest_grid(atoms[i]->pos, atoms[i]->element, atom_number_list, map);
      sum_d += d;
   }
   return sum_d;
}

float
coot::util::z_weighted_density_score_new(const std::vector<std::pair<mmdb::Atom *, float> > &atom_atom_number_pairs,
					 const clipper::Xmap<float> &map) {

   float sum_d = 0;
   for (unsigned int iat=0; iat<atom_atom_number_pairs.size(); iat++) {
      const mmdb::Atom *at = atom_atom_number_pairs[iat].first;
      clipper::Coord_orth co(at->x, at->y, at->z);
      float d = coot::util::density_at_point(map, co) * atom_atom_number_pairs[iat].second;
      sum_d += d;
   }
   return sum_d;
}



// High density fitting is down-weighted and negative density fitting
// is upweighted.
//
float
coot::util::biased_z_weighted_density_score(const minimol::molecule &mol,
					    const std::vector<std::pair<std::string, int> > &atom_number_list,
					    const clipper::Xmap<float> &map) {

   float sum_d = 0;
   std::vector<coot::minimol::atom *> atoms = mol.select_atoms_serial();
   for (unsigned int i=0; i<atoms.size(); i++) {
      float d = coot::util::z_weighted_density_at_point(atoms[i]->pos, atoms[i]->element,
							atom_number_list, map);
      float d_v = -exp(-d)+1.0;
      sum_d += d_v;
   }
   return sum_d;
}

std::vector<std::pair<std::string, float> >
coot::util::score_atoms(const minimol::residue &residue,
			const clipper::Xmap<float> &xmap) {

   std::vector<std::pair<std::string, float> > v;
   for (unsigned int i=0; i<residue.n_atoms(); i++) {
      const clipper::Coord_orth &pt = residue[i].pos;
      float d = density_at_point(xmap, pt);
      std::pair<std::string, float> p(residue[i].name, d);
      v.push_back(p);
   }
   return v;
}




// annealing_factor is default arg 1.0, but (now) between 1.0 and 0.1 (inclusive)
// for 10 threads.
//
std::pair<clipper::RTop_orth, std::vector<mmdb::Atom> >
coot::util::jiggle_atoms(const std::vector<mmdb::Atom *> &atoms,
   const clipper::Coord_orth &centre_pt,
   float jiggle_trans_scale_factor,
   float annealing_factor) {

      if (annealing_factor <= 0)
      annealing_factor = 1.0;

      clipper::RTop_orth rtop = make_rtop_orth_for_jiggle_atoms(jiggle_trans_scale_factor,
                                                                annealing_factor);
         std::vector<mmdb::Atom> new_atoms(atoms.size());
         for (unsigned int i=0; i<atoms.size(); i++) {
            new_atoms[i].Copy(atoms[i]);
         }
         // now apply rtop to atoms (shift the atoms relative to the
         // centre_pt before doing the wiggle
         for (unsigned int i=0; i<atoms.size(); i++) {
            clipper::Coord_orth pt_rel(atoms[i]->x - centre_pt.x(),
            atoms[i]->y - centre_pt.y(),
            atoms[i]->z - centre_pt.z());
            clipper::Coord_orth new_pt = pt_rel.transform(rtop);
            new_pt += centre_pt;
            new_atoms[i].x = new_pt.x();
            new_atoms[i].y = new_pt.y();
            new_atoms[i].z = new_pt.z();
         }
   return std::pair<clipper::RTop_orth, std::vector<mmdb::Atom> > (rtop, new_atoms);
}

// We want to pass atoms of the same type as we get back (so that they
// can be dynamically updated, rather than jiggling the original
// atoms)
//
// annealing_factor is default arg 1.0
//
std::pair<clipper::RTop_orth, std::vector<mmdb::Atom> >
coot::util::jiggle_atoms(const std::vector<mmdb::Atom > &atoms,
			 const clipper::Coord_orth &centre_pt,
			 float jiggle_trans_scale_factor,
			 float annealing_factor) {

   std::vector<mmdb::Atom> new_atoms(atoms.size());
   // now apply rtop to atoms (shift the atoms relative to the
   // centre_pt before doing the wiggle
   clipper::RTop_orth rtop = make_rtop_orth_for_jiggle_atoms(jiggle_trans_scale_factor, annealing_factor);
   for (unsigned int i=0; i<atoms.size(); i++) {
      clipper::Coord_orth pt_rel(atoms[i].x - centre_pt.x(),
				 atoms[i].y - centre_pt.y(),
				 atoms[i].z - centre_pt.z());
      clipper::Coord_orth new_pt = pt_rel.transform(rtop);
      new_pt += centre_pt;
      new_atoms[i].x = new_pt.x();
      new_atoms[i].y = new_pt.y();
      new_atoms[i].z = new_pt.z();
   }
   return std::pair<clipper::RTop_orth, std::vector<mmdb::Atom> > (rtop, new_atoms);
}

clipper::RTop_orth
coot::util::make_rtop_orth_for_jiggle_atoms(float jiggle_trans_scale_factor,
                                            float annealing_factor) {

   // Read this:
   // https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices
   // make a quaternion where the 4 q values are sampled from a normal distribution
   // normalize
   // convert to 3x3 matrix

   float rmi = 1.0/float(RAND_MAX);
   float small_frac = coot::util::random() * rmi;

   double scaling_factor = 1.0;
   bool small_rot_flag = false;
   if (small_frac < 0.1) {
      scaling_factor = 0.01;
      small_rot_flag = true;
   }

   // If these angles are small, then we get small rotations of the model
   //
   double rand_ang_1 = 2.0 * M_PI * coot::util::random() * rmi;
   double rand_ang_2 = 2.0 * M_PI * coot::util::random() * rmi;
   double rand_ang_3 = 2.0 * M_PI * coot::util::random() * rmi;

   if (small_rot_flag) {
      rand_ang_1 -= M_PI;
      rand_ang_2 -= M_PI;
      rand_ang_3 -= M_PI;
      rand_ang_1 *= scaling_factor;
      rand_ang_2 *= scaling_factor;
      rand_ang_3 *= scaling_factor;
   }

   double r1 = (2.0 * coot::util::random() -1.0) * rmi;
   double r2 = (2.0 * coot::util::random() -1.0) * rmi;
   double r3 = (2.0 * coot::util::random() -1.0) * rmi;

   double rand_pos_1 = r1 * 1.1 * jiggle_trans_scale_factor * annealing_factor * scaling_factor;
   double rand_pos_2 = r2 * 1.1 * jiggle_trans_scale_factor * annealing_factor * scaling_factor;
   double rand_pos_3 = r3 * 1.1 * jiggle_trans_scale_factor * annealing_factor * scaling_factor;

   clipper::Euler<clipper::Rotation::EulerXYZr> e(rand_ang_1, rand_ang_2, rand_ang_3);
   clipper::Mat33<double> r = e.rotation().matrix();
   clipper::Coord_orth shift(rand_pos_1, rand_pos_2, rand_pos_3);
   clipper::RTop_orth rtop(r, shift);
   return rtop;

}

typedef clipper::NXmap<float>::Map_reference_index NRI;
typedef clipper::Xmap<float>::Map_reference_index  MRI;

std::vector<std::pair<clipper::Xmap<float>::Map_reference_index,
		      clipper::Xmap<float>::Map_reference_index> >
coot::make_map_reference_index_start_stops(const clipper::Xmap<float> &xmap, int n_threads_in) {

   unsigned int n_threads = n_threads_in;
   std::vector<std::pair<MRI, MRI> > map_ref_start_stops;
   bool debug = true;
   unsigned int count = 0;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next())
      count++;
   unsigned int n_per_thread = count/n_threads;
   if (n_per_thread * n_threads < count)
      n_per_thread += 1;
   unsigned int i_count = 0;
   MRI start_ix = xmap.first();
   MRI stop_ix;
   for (ix = xmap.first(); !ix.last(); ix.next()) {
      if (i_count == n_per_thread) {
	 stop_ix =  ix; // first out of range
	 std::pair<MRI, MRI> p(start_ix, stop_ix);
	 map_ref_start_stops.push_back(p);
	 start_ix = ix;
	 i_count = 0;
      } else {
	 i_count++;
      }
   }
   if (! map_ref_start_stops.back().second.last()) {
      std::pair<MRI, MRI> p(start_ix, ix);
      map_ref_start_stops.push_back(p);
   }

   return map_ref_start_stops;
}

std::vector<std::pair<NRI, NRI> >
coot::make_map_reference_index_start_stops(const clipper::NXmap<float> &nxmap, int n_threads) {

   std::vector<std::pair<NRI, NRI> > map_ref_start_stops(n_threads);
   bool debug = false;

   int nu = nxmap.grid().nu();
   int nv = nxmap.grid().nv();
   int nw = nxmap.grid().nw();
   int nu_step = nu/n_threads;
   // round that up if not exact
   if ((nu_step * n_threads) < nu)
      nu_step++;
   if (debug) {
      std::cout << "nxmap size() " << nxmap.grid().size() << std::endl;
      std::cout << "nxmap grid() " << nxmap.grid().format() << std::endl;
      std::cout << "nu " << nu << std::endl;
      std::cout << "nu_step " << nu_step << std::endl;
   }
   NRI grid_end(nxmap, clipper::Coord_grid(nu-1, nv-1, nw)); // first out-of-grid index
   for (int i=0; i<n_threads; i++) {
      NRI layer_start = NRI(nxmap, clipper::Coord_grid(nu_step*i,     0, 0));
      NRI layer_end   = NRI(nxmap, clipper::Coord_grid(nu_step*(i+1), 0, 0));
      if (layer_end.index() > nxmap.grid().size())
         layer_end = grid_end;
      map_ref_start_stops[i] = std::pair<NRI, NRI> (layer_start, layer_end);
      if (debug)
       	 std::cout << "debug::" << layer_start.index() << " " << layer_end.index() << std::endl;
   }
   return map_ref_start_stops;
}

void
xmap_to_nxmap_workpackage(const clipper::Xmap<float> &xmap,
			  clipper::NXmap<float> *nxmap_p,
			  const std::pair<NRI, NRI> &start_stop) {

   // std::cout << "starting workpackage" << std::endl;

   clipper::Coord_grid offset =
      xmap.coord_map(nxmap_p->coord_orth(clipper::Coord_map(0.0, 0.0, 0.0))).coord_grid();

   // std::cout << "debug:: " << start_stop.first.index() << " " << start_stop.second.index() << std::endl;

   clipper::Xmap<float>::Map_reference_coord ix(xmap);
   for (NRI inx = start_stop.first; inx.index() != start_stop.second.index(); inx.next()) {
      ix.set_coord(inx.coord() + offset);
      (*nxmap_p)[inx] = xmap[ix];
   }
}



clipper::NXmap<float>
coot::util::make_nxmap(const clipper::Xmap<float> &xmap, mmdb::Manager *mol, int SelectionHandle, float border) {

   // if we are given an EM map then we presume that the molecule is nicely situated in
   // that map. It might not be, but we will fix that later.

   if (is_EM_map(xmap))
      return make_nxmap_from_EM_P1_map(xmap);
   else
      return make_nxmap_from_xmap(xmap, mol, SelectionHandle, border);

}

clipper::NXmap<float>
coot::util::make_nxmap_from_EM_P1_map(const clipper::Xmap<float> &xmap) {

   clipper::Cell cell = xmap.cell();
   clipper::Grid_sampling grid_sampling = xmap.grid_sampling();
   clipper::Grid_range gr_asu = xmap.grid_asu();
   clipper::NXmap<float> nxmap(cell, grid_sampling, gr_asu);

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next() )  { // iterator index.
      clipper::Coord_grid cg = ix.coord();
      nxmap.set_data(cg, xmap[ix]);
   }
   return nxmap;
}


clipper::NXmap<float>
coot::util::make_nxmap_from_xmap(const clipper::Xmap<float> &xmap, mmdb::Manager *mol, int SelectionHandle, float border) {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   bool debug = false;

   std::pair<clipper::Coord_orth, clipper::Coord_orth> p = util::extents(mol, SelectionHandle);

   clipper::Coord_orth p1 = p.first;
   clipper::Coord_orth p2 = p.second;

   // std::cout << "debug:: make_nxmap() extents " << p.first.format() << " " << p.second.format() << std::endl;

   p1 -= clipper::Coord_orth(border,border,border);
   p2 += clipper::Coord_orth(border,border,border);

   std::pair<clipper::Coord_orth, clipper::Coord_orth> pp(p1, p2);

   std::pair<clipper::Coord_frac, clipper::Coord_frac> fc_pair =
      find_struct_fragment_coord_fracs_v2(pp, xmap.cell());

   clipper::Cell cell = xmap.cell();
   clipper::Grid_sampling grid_sampling = xmap.grid_sampling();

   float radius = sqrt((p2-p1).lengthsq());
   clipper::Coord_orth comg(0.5*(p1.x()+p2.x()),
			    0.5*(p1.y()+p2.y()),
			    0.5*(p1.z()+p2.z()));

   // make a non-cube, ideally - needs different extents to do that.

   // get grid range
   // gr0: a grid range of the correct size (at the origin, going + and - in small box)
   // gr1: a grid range of the correct size (around the correct place, comg)
   clipper::Grid_range gr0(cell, grid_sampling, radius);
   clipper::Grid_range gr1(gr0.min() + comg.coord_frac(cell).coord_grid(grid_sampling),
			   gr0.max() + comg.coord_frac(cell).coord_grid(grid_sampling));

   // Here I need to update the grid range, gr1 to get a "good" radix (radices?)

   // this constructor will fail throwing a std::length_error if the map is too big.
   //
   // init nxmap
   clipper::NXmap<float> nxmap(cell, grid_sampling, gr1);
   clipper::Xmap<float>::Map_reference_coord ix(xmap);
   clipper::Coord_grid offset =
      xmap.coord_map(nxmap.coord_orth(clipper::Coord_map(0.0,0.0,0.0))).coord_grid();
   int n_threads = 4; // the number of layers
   std::vector<std::pair<NRI, NRI> > map_ref_start_stops(n_threads);
   int nu = nxmap.grid().nu();
   int nv = nxmap.grid().nv();
   int nw = nxmap.grid().nw();
   int nu_step = nu/n_threads;
   // round that up if not exact
   if ((nu_step * n_threads) < nu)
      nu_step++;
   if (debug) {
      std::cout << "nxmap size() " << nxmap.grid().size() << std::endl;
      std::cout << "nxmap grid() " << nxmap.grid().format() << std::endl;
      std::cout << "nu " << nu << " nv " << nv << " nw " << nw << std::endl;
      std::cout << "nu_step " << nu_step << std::endl;
   }

   // MaxInt (32): 2,147,483,647

   // Is there a more elegant way to get an NRI for the end of the grid?
   NRI grid_end(nxmap, clipper::Coord_grid(nu-1, nv-1, nw)); // first bad index
   if (debug)
      std::cout << "grid_end: " << grid_end.index() << " " << grid_end.coord().format() << std::endl;
   for (int i=0; i<n_threads; i++) {
      NRI layer_start = NRI(nxmap, clipper::Coord_grid(nu_step*i,     0, 0));
      NRI layer_end   = NRI(nxmap, clipper::Coord_grid(nu_step*(i+1), 0, 0));
      if (layer_end.index() > nxmap.grid().size())
	 layer_end = grid_end;
      map_ref_start_stops[i] = std::pair<NRI, NRI> (layer_start, layer_end);
      if (debug)
	 std::cout << "debug::" << layer_start.index() << " " << layer_end.index() << std::endl;
   }

   bool single_thread_method = false;
   if (single_thread_method) {
      for (NRI inx = nxmap.first(); !inx.last(); inx.next()) {
	 ix.set_coord(inx.coord() + offset);
	 nxmap[inx] = xmap[ix];
      }
   } else {
      std::vector<std::thread> threads;
      for (int i=0; i<n_threads; i++) {
	 threads.push_back(std::thread(xmap_to_nxmap_workpackage,
				       std::cref(xmap), &nxmap, std::cref(map_ref_start_stops[i])));
	 std::this_thread::sleep_for(std::chrono::microseconds(1));
      }

      for (int i=0; i<n_threads; i++)
	 threads[i].join();
   }
   if (debug)
      std::cout << "returning from make_nxmap() " << nxmap.grid().format() << std::endl;
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   if (false)
      std::cout << "------------------- make_nxmap() timing: " << d10 << " milliseconds " << std::endl;
   return nxmap;
}


// border is default argument, value 3.
//
clipper::NXmap<float>
coot::util::make_nxmap(const clipper::Xmap<float> &xmap, atom_selection_container_t asc, float border) {

   return make_nxmap(xmap, asc.mol, asc.SelectionHandle, border);
}

#include "clipper/contrib/edcalc.h"

//
clipper::NXmap<float>
coot::util::make_edcalc_map(const clipper::NXmap<float>& map_ref,  // for metrics
			    mmdb::Manager *mol, int atom_selection_handle) {

   // init nxmap
   clipper::NXmap<float> nxmap(map_ref.grid(), map_ref.operator_orth_grid());

   clipper::EDcalc_iso<float> edc(3.0); // radius
   mmdb::Atom **sel_atoms = 0;
   int n_sel_atoms;
   std::vector<clipper::Atom> l;
   mol->GetSelIndex(atom_selection_handle, sel_atoms, n_sel_atoms);
   for (int ii=0; ii<n_sel_atoms; ii++) {
      mmdb::Atom *at = sel_atoms[ii];
      std::string ele(at->element);
      clipper::Coord_orth pt(at->x, at->y, at->z);
      clipper::Atom cat;
      cat.set_element(ele);
      cat.set_coord_orth(pt);
      cat.set_u_iso(at->tempFactor * 0.0125);
      // cat.set_u_iso(0.1);
      cat.set_occupancy(1.0);
      l.push_back(cat);
   }

   clipper::Atom_list al(l);
   edc(nxmap, al);

   return nxmap;

}
