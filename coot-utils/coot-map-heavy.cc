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




// annealing_factor is default arg 1.0
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

   // If these angles are small, then we get small rotations of the model
   // 
   double rand_ang_1 = 2* M_PI * coot::util::random() * rmi * annealing_factor;
   double rand_ang_2 = 2* M_PI * coot::util::random() * rmi * annealing_factor;
   double rand_ang_3 = 2* M_PI * coot::util::random() * rmi * annealing_factor;

   double rand_pos_1 = coot::util::random() * rmi * 0.1 * jiggle_trans_scale_factor * annealing_factor;
   double rand_pos_2 = coot::util::random() * rmi * 0.1 * jiggle_trans_scale_factor * annealing_factor;
   double rand_pos_3 = coot::util::random() * rmi * 0.1 * jiggle_trans_scale_factor * annealing_factor;
   clipper::Euler<clipper::Rotation::EulerXYZr> e(rand_ang_1, rand_ang_2, rand_ang_3);
   clipper::Mat33<double> r = e.rotation().matrix();
   clipper::Coord_orth shift(rand_pos_1, rand_pos_2, rand_pos_3);
   clipper::RTop_orth rtop(r, shift);
   return rtop;

}

// needed?
#include "atom-selection-container.hh"
#include "coords/mmdb-crystal.h"

clipper::NXmap<float>
coot::util::make_nxmap(const clipper::Xmap<float> &xmap, mmdb::Manager *mol, int SelectionHandle, float border) {

   std::pair<clipper::Coord_orth, clipper::Coord_orth> p = util::extents(mol, SelectionHandle);

   clipper::Coord_orth p1 = p.first;
   clipper::Coord_orth p2 = p.second;

   std::cout << "debug:: make_nxmap() extents " << p.first.format() << " " << p.second.format() << std::endl;

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
   typedef clipper::NXmap<float>::Map_reference_index NRI;
   for (NRI inx = nxmap.first(); !inx.last(); inx.next()) {
      ix.set_coord(inx.coord() + offset);
      nxmap[inx] = xmap[ix];
   }
   return nxmap;
}


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
