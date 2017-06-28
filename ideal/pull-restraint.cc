
#include "simple-restraint.hh"


mmdb::Atom *
coot::restraints_container_t::add_atom_pull_restraint(atom_spec_t spec, clipper::Coord_orth pos) {

   mmdb::Atom *at = 0;

   clear_atom_pull_restraint();  // clear old ones
   for (int iat=0; iat<n_atoms; iat++) { 
      atom_spec_t atom_spec(atom[iat]);
      if (atom_spec == spec) {
	 if (! fixed_check(iat)) { 
	    add_target_position_restraint(iat, pos);
	    at = atom[iat];
	 }
	 break;
      }
   }
   return at;
}

#include <algorithm>


// clear them all - but at the moment the user can only set one of them.
void
coot::restraints_container_t::clear_atom_pull_restraint() {

   unsigned int pre_size = restraints_vec.size();
   if (pre_size > 0) {
      restraints_vec.erase(std::remove_if(restraints_vec.begin(),
					  restraints_vec.end(),
					  target_position_eraser),
			   restraints_vec.end());
      unsigned int post_size = restraints_vec.size();
      // std::cout << "debug:: clear_atom_pull_restraint() pre size: " << pre_size << " post size: "
      // << post_size << std::endl;
   }
}

bool coot::target_position_eraser(const simple_restraint &r) {
   return (r.restraint_type == restraint_type_t(TARGET_POS_RESTRANT));
}



double
coot::distortion_score_target_pos(const coot::simple_restraint &rest,
				  void *params,
				  const gsl_vector *v) {

   bool harmonic_restraint = false; // 20170628 let's try the (f) ln(cosh) and (df) hyperbolic tan
                                    // style restraint

   // get the current coordinates
   int idx = 3*(rest.atom_index_1);
   clipper::Coord_orth current_pos(gsl_vector_get(v,idx), 
                                   gsl_vector_get(v,idx+1), 
                                   gsl_vector_get(v,idx+2));


   double dist = clipper::Coord_orth::length(current_pos, rest.atom_pull_target_pos);

   if (harmonic_restraint) {
      double sigma = 0.04; // (slightly refined) guess, copy below
      double weight = 1.0/(sigma*sigma);
      // std::cout << "distortion_score_target_pos() returning " << weight * dist * dist << std::endl;
      return weight * dist * dist;
   } else {
      // f = some_scale * ln (cosh(z))
      // where z = dist/top_out_dist
      // Note that f is practically max (1) when z = 2.

      // duplicate these settings in the my_df_target_pos().
      //
      double top_out_dist = 4.0; // Angstroms, needs tweaking?
      double scale = 30000.0;       // needs tweaking

      double z = dist/top_out_dist;
      double e_z = exp(z);
      double cosh_z = 0.5 * (e_z + 1.0/e_z);
      double ln_cosh_z = log(cosh_z);
      return scale * ln_cosh_z;
   }
}

void coot::my_df_target_pos(const gsl_vector *v, 
			    void *params, 
			    gsl_vector *df) {

   bool harmonic_restraint = false; // 20170628 let's try the (f) ln(cosh) and (df) hyperbolic tan
                                    // style restraint

   restraints_container_t *restraints_p = static_cast<restraints_container_t *>(params);

   // patch madness
   // std::cout << "my_df_target_pos "
   // << restraints_p->restraints_limits_target_pos.first << " "
   // << restraints_p->restraints_limits_target_pos.second<< std::endl;

   int restraints_size = restraints_p->size();
   // 20170628 investigate using restraints_limits_target_pos
   for (int i=0; i<restraints_size; i++) {
      const simple_restraint &rest = (*restraints_p)[i];
      if (rest.restraint_type == TARGET_POS_RESTRANT) {
	 double sigma = 0.04; // change as above in distortion score
	 int idx = 3*(rest.atom_index_1);

// 	 clipper::Coord_orth current_pos(gsl_vector_get(v,idx),
// 					 gsl_vector_get(v,idx+1),
// 					 gsl_vector_get(v,idx+2));
	 
	 if (harmonic_restraint) {
	    double constant_part = 2.0 / (sigma * sigma);

	    double dist_x = gsl_vector_get(v, idx)   - rest.atom_pull_target_pos[0];
	    double dist_y = gsl_vector_get(v, idx+1) - rest.atom_pull_target_pos[1];
	    double dist_z = gsl_vector_get(v, idx+2) - rest.atom_pull_target_pos[2];
	    // double squared_dist = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;

	    *gsl_vector_ptr(df, idx  ) += constant_part * dist_x;
	    *gsl_vector_ptr(df, idx+1) += constant_part * dist_y;
	    *gsl_vector_ptr(df, idx+2) += constant_part * dist_z;

	 } else {

	    double scale = 30000.0;       // needs tweaking
	    double top_out_dist = 4.0;   // Angstroms, needs tweaking?
	    double k = 1.0 / top_out_dist;

	    clipper::Coord_orth current_pos(gsl_vector_get(v,idx),
					    gsl_vector_get(v,idx+1),
					    gsl_vector_get(v,idx+2));
	    double dist_x = gsl_vector_get(v, idx)   - rest.atom_pull_target_pos[0];
	    double dist_y = gsl_vector_get(v, idx+1) - rest.atom_pull_target_pos[1];
	    double dist_z = gsl_vector_get(v, idx+2) - rest.atom_pull_target_pos[2];

	    double dist = clipper::Coord_orth::length(current_pos, rest.atom_pull_target_pos);
	    double z = dist/top_out_dist;
	    double e_2z = exp(2.0 * z);
	    double tanh_z = (e_2z - 1) / (e_2z + 1);
	    double constant_part = scale * k * tanh_z / dist;

	    *gsl_vector_ptr(df, idx  ) += constant_part * dist_x;
	    *gsl_vector_ptr(df, idx+1) += constant_part * dist_y;
	    *gsl_vector_ptr(df, idx+2) += constant_part * dist_z;

	 }
      }
   }
   
   
}


bool
coot::restraints_container_t::turn_off_when_close_target_position_restraint() {

   bool status = false;
   unsigned int pre_size = restraints_vec.size();
   restraints_vec.erase(std::remove_if(restraints_vec.begin(),
				       restraints_vec.end(),
				       turn_off_when_close_target_position_restraint_eraser(atom, n_atoms)),
			restraints_vec.end());
   unsigned int post_size = restraints_vec.size();
   if (post_size < pre_size) status = true;
   return status;
}


