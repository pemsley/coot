
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

   // get the current coordinates
   int idx = 3*(rest.atom_index_1);
   clipper::Coord_orth current_pos(gsl_vector_get(v,idx), 
                                   gsl_vector_get(v,idx+1), 
                                   gsl_vector_get(v,idx+2));

   double sigma = 0.04; // (slightly refined) guess, copy below
   double weight = 1.0/(sigma*sigma);
   double dist = clipper::Coord_orth::length(current_pos, rest.atom_pull_target_pos);
   // std::cout << "distortion_score_target_pos() returning " << weight * dist * dist << std::endl;
   return weight * dist * dist;
}

void coot::my_df_target_pos(const gsl_vector *v, 
			    void *params, 
			    gsl_vector *df) {

   restraints_container_t *restraints = (restraints_container_t *)params;
   for (int i=0; i<restraints->size(); i++) {
      if ( (*restraints)[i].restraint_type == TARGET_POS_RESTRANT) {
	 const simple_restraint &rest = (*restraints)[i];
	 double sigma = 0.04; // change as above in distortion score
	 int idx = 3*(rest.atom_index_1);

// 	 clipper::Coord_orth current_pos(gsl_vector_get(v,idx),
// 					 gsl_vector_get(v,idx+1),
// 					 gsl_vector_get(v,idx+2));
	 
         double constant_part = 2.0 / (sigma * sigma);
         
         double dist_x = gsl_vector_get(v, idx)   - rest.atom_pull_target_pos[0];
         double dist_y = gsl_vector_get(v, idx+1) - rest.atom_pull_target_pos[1];
         double dist_z = gsl_vector_get(v, idx+2) - rest.atom_pull_target_pos[2];
         // double squared_dist = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;

	 *gsl_vector_ptr(df, idx  ) += constant_part * dist_x;
	 *gsl_vector_ptr(df, idx+1) += constant_part * dist_y;
	 *gsl_vector_ptr(df, idx+2) += constant_part * dist_z;
	 
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


