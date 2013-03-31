/* ideal/extra-restraints-kk.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2011 by Kevin Keating
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


// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL

#include "simple-restraint.hh"

double
coot::distortion_score_start_pos(const coot::simple_restraint &start_pos_restraint,
			    void *params,
                            const gsl_vector *v) {
   
   // first extract the object from params 
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 

   if (int(v->size) != int(restraints->init_positions_size()) ) {
      cout << "very worry. A bug. " << v->size << " "
	   << restraints->init_positions_size() << endl;
      return 0.0;
   }

   // get the current coordinates
   int idx = 3*(start_pos_restraint.atom_index_1 - 0);
   clipper::Coord_orth current_pos(gsl_vector_get(v,idx), 
                                   gsl_vector_get(v,idx+1), 
                                   gsl_vector_get(v,idx+2));
   
   // get the original coordinates
   clipper::Coord_orth start_pos(restraints->initial_position(idx),
                                 restraints->initial_position(idx+1),
                                 restraints->initial_position(idx+2));
   
   double weight = 1.0/(start_pos_restraint.sigma * start_pos_restraint.sigma);
   double dist = clipper::Coord_orth::length(current_pos, start_pos);   
   return weight * dist * dist;
}


void coot::my_df_start_pos (const gsl_vector *v, 
		      void *params, 
		      gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   
   if (int(v->size) != int(restraints->init_positions_size()) ) {
      cout << "very worry. A bug. " << v->size << " "
	   << restraints->init_positions_size() << endl;
      return;
   }
   
   double val;
   int idx;
   for (int i=0; i<restraints->size(); i++) {
      
      if ( (*restraints)[i].restraint_type == coot::START_POS_RESTRAINT) {
         
         idx = 3*((*restraints)[i].atom_index_1); 
         
         double sigma = (*restraints)[i].sigma;
         double constant_part = 2.0 / (sigma * sigma);
         
         double dist_x = gsl_vector_get(v, idx)   - restraints->initial_position(idx);
         double dist_y = gsl_vector_get(v, idx+1) - restraints->initial_position(idx+1);
         double dist_z = gsl_vector_get(v, idx+2) - restraints->initial_position(idx+2);
         double squared_dist = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
         
         val = constant_part * dist_x;
         gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + val);
      
         val = constant_part * dist_y;
         gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + val);
      
         val = constant_part * dist_z;
         gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + val);
         
      }
   }
}

#endif // HAVE_GSL
