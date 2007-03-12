/* ligand/coot-map-heavy.cc
 * 
 * Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */



#include "clipper/core/map_interp.h"
#include "clipper/core/hkl_compute.h"

#include "coot-map-heavy.hh"
#include "coot-map-utils.hh"
#include "coot-coord-utils.hh"


#ifdef HAVE_GSL
// return status success (something moved = 1) or error/failure (0)
//
int
coot::util::fit_to_map_by_simplex(PPCAtom atom_selection,
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
   
   minex_func.f = &coot::util::my_f_simplex_internal;
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
	 simplex_apply_shifts_internal(s->x, par);
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
coot::util::setup_simplex_x_internal(gsl_vector *x,
				     PPCAtom atom_selection,
				     int n_selected_atoms) {

   // unnecessary, I think.  All rots and trans angles are set to 0.
}

// Simplex Refinement is finished, force-feed the results (new
// coordinates) back into the atom_selection.
// 
void
coot::util::simplex_apply_shifts_internal(gsl_vector *s,
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
coot::util::my_f_simplex_internal (const gsl_vector *v,
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

