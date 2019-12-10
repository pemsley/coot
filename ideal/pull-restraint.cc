#include <algorithm>

#include "simple-restraint.hh"


mmdb::Atom *
coot::restraints_container_t::add_atom_pull_restraint(const atom_spec_t &spec, clipper::Coord_orth pos) {

   mmdb::Atom *at = 0;

   // 20180217 now we replace the target position if we can, rather than delete and add.

   std::vector<simple_restraint>::iterator it;
   for (it=restraints_vec.begin(); it!=restraints_vec.end(); it++) {
      if (it->restraint_type == restraint_type_t(TARGET_POS_RESTRAINT)) {
	 if (it->atom_spec == spec) {
	    at = atom[it->atom_index_1];

	    // wait until you get the lock
	    bool unlocked = false;
	    while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
	       std::this_thread::sleep_for(std::chrono::nanoseconds(10));
	       unlocked = false;
	    }
            bool is_different = true;
            double d2 = (it->atom_pull_target_pos-pos).lengthsq();
            if (d2 < 0.0001)
               is_different = false;
	    it->atom_pull_target_pos = pos;
            it->is_closed = false; // if it had been closed before, it's active again
            if (is_different)
               needs_reset = true;
	    restraints_lock = false; // unlocked
            if (false) // debugging
               std::cout << "add_atom_pull_restraint() update position for " << it->atom_index_1 << " "
                         << atom_spec_t(at) << " " << pos.format() << "\n";
	    break;
	 }
      }
   }

   if (! at) {
      for (int iat=0; iat<n_atoms; iat++) { 
	 atom_spec_t atom_spec(atom[iat]);
	 if (atom_spec == spec) {
	    if (! fixed_check(iat)) { 
	       add_target_position_restraint(iat, spec, pos);
	       at = atom[iat];
	    }
	    break;
	 }
      }
   }

   // needs_reset = true; // always true makes the refinement smoother for some reason.
   return at;
}


void
coot::restraints_container_t::add_target_position_restraint(int idx, const atom_spec_t &spec,
							    clipper::Coord_orth &target_pos) {

   simple_restraint r(TARGET_POS_RESTRAINT, idx, spec, target_pos);

#ifdef HAVE_CXX_THREAD
   // wait until you get the lock
   bool unlocked = false;
   while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
      // std::cout << "waiting in add_target_position_restraint()" << std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
      unlocked = false;
   }

   if (false) // debug
      std::cout << "add_target_position_restraint() for " << idx << " "
                << spec << target_pos.format() << "\n";

   restraints_vec.push_back(r); // push_back_restraint()
   post_add_new_restraint(); // adds new restraint to one of the vectors of the restraint indices
   restraints_lock = false; // unlock
   needs_reset = true;
#endif // HAVE_CXX_THREAD
}



void
coot::restraints_container_t::clear_atom_pull_restraint(const coot::atom_spec_t &spec) {

   if (false)
      std::cout << "restraints_container_t clear_atom_pull_restraint for " << spec
		<< " called " << std::endl;

   unsigned int pre_size = size();
   if (pre_size > 0) {
      // extra loop needed.
      /*
      restraints_vec.erase(std::remove_if(restraints_vec.begin(),
					  restraints_vec.end(),
					  target_position_for_atom_eraser(spec)),
			   restraints_vec.end());
      */
      for (unsigned int i=0; i<restraints_vec.size(); i++) {
         if (restraints_vec[i].restraint_type == TARGET_POS_RESTRAINT)
            if (restraints_vec[i].atom_spec == spec)
               restraints_vec[i].close();
      }
      unsigned int post_size = size();
      if (false)
	 std::cout << "debug:: clear_atom_pull_restraint() pre size: " << pre_size << " post size: "
		   << post_size << std::endl;
   }
}


// clear them all
void
coot::restraints_container_t::clear_all_atom_pull_restraints() {

   bool unlocked = false;
   while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::nanoseconds(100));
      unlocked = false;
   }
   unsigned int pre_size = size();
   if (pre_size > 0) {
      /*
      restraints_vec.erase(std::remove_if(restraints_vec.begin(),
					  restraints_vec.end(),
					  target_position_eraser),
			   restraints_vec.end());
      */
      for (unsigned int i=0; i<restraints_vec.size(); i++) {
         if (restraints_vec[i].restraint_type == TARGET_POS_RESTRAINT)
            restraints_vec[i].close();
      }
      unsigned int post_size = size();
      if (false)
	 std::cout << "debug:: clear_atom_pull_restraint() pre size: " << pre_size << " post size: "
		   << post_size << std::endl;
   }
   restraints_lock = false;
   needs_reset = true; // probably
}

bool coot::target_position_eraser(const simple_restraint &r) {
   return (r.restraint_type == restraint_type_t(TARGET_POS_RESTRAINT));
}



double
coot::distortion_score_target_pos(const coot::simple_restraint &rest,
				  double scale_factor,
				  const gsl_vector *v) {

   if (rest.is_closed) return 0.0;

   // void *params, no longer passed

   // 20170628 let's try the (f) ln(cosh) and (df) hyperbolic tan
   // style restraint
   //
   // 20180602, let's try harmonic again
   //
   bool harmonic_restraint = true;

   // get the current coordinates
   int idx = 3*(rest.atom_index_1);
   clipper::Coord_orth current_pos(gsl_vector_get(v,idx), 
                                   gsl_vector_get(v,idx+1), 
                                   gsl_vector_get(v,idx+2));

   // restraints_container_t *restraints_p = static_cast<restraints_container_t *>(params);

   // we don't need to calculate dist
   //double dist = clipper::Coord_orth::length(current_pos, rest.atom_pull_target_pos);
   double dist_sq = (current_pos-rest.atom_pull_target_pos).lengthsq();

   if (harmonic_restraint) {
      double sigma = 0.03; // (slightly refined) guess, copy below
      double weight = 1.0/(sigma*sigma);
      // return weight * dist * dist;
      return weight * dist_sq;
   } else {
      // f = some_scale * ln (cosh(z))
      // where z = dist/top_out_dist
      // Note that f is practically max (1) when z = 2.

      // duplicate these settings in the my_df_target_pos().
      //
      double top_out_dist = 4.0; // Angstroms, needs tweaking?

      // double scale = log_cosh_target_distance_scale_factor;       // needs tweaking
      double scale = scale_factor;

      double dist = sqrt(dist_sq);
      double z = dist/top_out_dist;
      double e_z = exp(z);
      double cosh_z = 0.5 * (e_z + 1.0/e_z);
      double ln_cosh_z = log(cosh_z);
      // std::cout << " " << scale * ln_cosh_z << std::endl;
      return scale * ln_cosh_z;
   }
}

void coot::my_df_target_pos(const gsl_vector *v, 
			    void *params, 
			    gsl_vector *df) {

   // 20170628 let's try the (f) ln(cosh) and (df) hyperbolic tan
   // style restraint
   // 20180602-PE let's try harmonic again

   bool harmonic_restraint = true;
   
   restraints_container_t *restraints_p = static_cast<restraints_container_t *>(params);

   // patch madness
   // std::cout << "my_df_target_pos "
   // << restraints_p->restraints_limits_target_pos.first << " "
   // << restraints_p->restraints_limits_target_pos.second<< std::endl;

   int restraints_size = restraints_p->size();
   // 20170628 investigate using restraints_limits_target_pos
   for (int i=0; i<restraints_size; i++) {
      const simple_restraint &rest = (*restraints_p)[i];
      if (rest.restraint_type == TARGET_POS_RESTRAINT) {
         if (rest.is_closed) continue;
	 double sigma = 0.03; // change as above in distortion score
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

	    double scale = restraints_p->log_cosh_target_distance_scale_factor;
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


std::vector<coot::atom_spec_t>
coot::restraints_container_t::turn_off_atom_pull_restraints_when_close_to_target_position(const atom_spec_t &dragged_atom_spec) {

   // If, as we drag and the refinement proceeds, some atoms get close to their target positions, then
   // we should drop the restraint - but not for the atom being dragged - that shouldn't lose its restraint.
   // We want a list of atoms that are now close to their target position (mostly is empty).

   std::vector<atom_spec_t> v; // a list of specs for removed restraints (they no longer need to be drawn)

   double close_dist = 0.6;  // was 0.5; // was 0.4 [when there was just 1 pull atom restraint]
                             // sync with above function, or make a data member

   // we don't want to do this at the same time as clear_all_atom_pull_retraints

   bool unlocked = false;
   while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::nanoseconds(100));
      unlocked = false;
   }
   // set the return value:
   //
   std::vector<simple_restraint>::iterator it;
   for(it=restraints_vec.begin(); it!=restraints_vec.end(); it++) {
      if (it->restraint_type == restraint_type_t(TARGET_POS_RESTRAINT)) {
	 mmdb::Atom *at = atom[it->atom_index_1];
	 if (atom_spec_t(at) != dragged_atom_spec) {
	    clipper::Coord_orth pos(at->x, at->y, at->z);
	    double d = sqrt((pos - it->atom_pull_target_pos).lengthsq());
	    if (d < close_dist) {
               it->close();
	       v.push_back(it->atom_spec);
	    }
	 }
      }
   }

   // no need to remove the pull restraint! Just turn it off.

   // now do the remove
//    restraints_vec.erase(std::remove_if(restraints_vec.begin(),
// 				       restraints_vec.end(),
// 				       turn_off_when_close_target_position_restraint_eraser(close_dist, atom, n_atoms, dragged_atom_spec)),
// 			restraints_vec.end());

   restraints_lock = false;

   if (v.size() > 0)
      needs_reset = true;

   return v;
}
