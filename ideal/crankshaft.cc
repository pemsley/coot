
#include "crankshaft.hh"

#include "geometry/main-chain.hh"
#include "coot-utils/coot-coord-utils.hh"

coot::crankshaft_set::crankshaft_set(mmdb::Residue *res_0,
				     mmdb::Residue *res_1,
				     mmdb::Residue *res_2,
				     mmdb::Residue *res_3) {

   if (! res_0) throw(std::runtime_error("Null residue 0"));
   if (! res_1) throw(std::runtime_error("Null residue 1"));
   if (! res_2) throw(std::runtime_error("Null residue 2"));
   if (! res_3) throw(std::runtime_error("Null residue 3"));

   v.resize(8, 0);
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms;

   ca_1 = 0; // initally unset
   ca_2 = 0;

   mmdb::Residue *residue_p = res_0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      std::string at_name = at->name;
      if (at_name == " C  ") { // PDBv3 fixme
	 v[0] = at;
	 break;
      }
   }
   residue_p = res_1;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      std::string at_name = at->name;
      if (at_name == " N  ") { // PDBv3 fixme
	 v[1] = at;
      }
      if (at_name == " C  ") { // PDBv3 fixme
	 v[2] = at;
      }
      if (at_name == " O  ") { // PDBv3 fixme
	 v[3] = at;
      }
      if (at_name == " CA ") { // PDBv3 fixme
	 ca_1 = at;
      }
   }
   residue_atoms = 0;
   n_residue_atoms = 0;
   residue_p = res_2;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      std::string at_name = at->name;
      if (at_name == " N  ") { // PDBv3 fixme
	 v[4] = at;
      }
      if (at_name == " H  ") { // PDBv3 fixme
	 v[5] = at;
      }
      if (at_name == " C  ") { // PDBv3 fixme
	 v[6] = at;
      }
      if (at_name == " CA ") { // PDBv3 fixme
	 ca_2 = at;
      }
   }
   residue_atoms = 0;
   n_residue_atoms = 0;
   residue_p = res_3;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      std::string at_name = at->name;
      if (at_name == " N  ") { // PDBv3 fixme
	 v[7] = at;
      }
   }

   if (! ca_1) throw(std::runtime_error("missing ca_1"));
   if (! ca_2) throw(std::runtime_error("missing ca_2"));

   // check that we have all the atoms - except the HN
   int n_atoms = 0;
   for (std::size_t i=0; i<v.size(); i++) {
      if (v[i])
	 n_atoms++;
      else
	 if (i==5)
	    n_atoms++;
   }
   if (n_atoms == 8) {
      // happy path
   } else {
      throw(std::runtime_error("missing a mainchain atom"));
   }
}

// n_samples = 60 default arg
std::vector<std::pair<float, float> >
coot::crankshaft::spin_search(const coot::residue_spec_t &spec,
			      const zo::rama_table_set &zorts,
			      int n_samples) const {

   std::vector<std::pair<float, float> > v;

   std::pair<mmdb::Residue *, mmdb::Residue *> rs = util::get_this_and_next_residues(spec, mol);
   if (rs.first) {
      if (rs.second) {
	 mmdb::Residue *res_1 = rs.first;
	 mmdb::Residue *res_2 = rs.second;
	 mmdb::Atom *ca_1 = get_atom(res_1, " CA ");
	 mmdb::Atom *ca_2 = get_atom(res_2, " CA ");
	 if (ca_1 && ca_2) {
	    // when we crankshaft the C, O and N of res_1(n) and res_2(n+1),
	    // we need to know the C of residue n-1 and the N of residue n+2
	    // because the phi,psi of n and n+1 depend on them.
	    mmdb::Residue *res_0 = util::get_previous_residue(spec, mol);
	    if (res_0) {
	       // res_3 is needed for phi,psi of res_2
	       residue_spec_t spec_2(res_2);
	       mmdb::Residue *res_3 = util::get_following_residue(spec_2, mol);
	       if (res_3) {

		  std::string rt1 = "ALL!nP"; // hack
		  std::string rt2 = "ALL!nP";

		  float div = 1/float(n_samples);
		  crankshaft_set cs(res_0, res_1, res_2, res_3);
		  for (int i=0; i<n_samples; i++) {
		     float a = float(i) * div * 2*M_PI;
		     std::pair<phi_psi_t, phi_psi_t> ppp = cs.phi_psis(a);
		     std::pair<float,float> pr = probability_of_spin_orientation(ppp, rt1, rt2, zorts);
		     if (false)
			std::cout << "score: " << spec << " " << i << " "
				  << pr.first << " " << pr.second << std::endl;
		  }
	       }
	    }
	 } else {
	    std::cout << "missing mainchain atom(s) for " << spec << std::endl;
	 }
      } else {
	 std::cout << "missing second residue for " << spec << std::endl;
      }
   } else {
      std::cout << "missing first residue: " << spec << std::endl;
   }
   return v;
}

void
coot::crankshaft_set::move_the_atoms(float ang) {

   // move atom indices 2 3 4 5 (if not null)

   if (v.size() > 5) {
      int indices[] = { 2, 3, 4, 5};
      clipper::Coord_orth p_ca_1 = co(ca_1);
      clipper::Coord_orth p_ca_2 = co(ca_2);
      clipper::Coord_orth dir = p_ca_2 - p_ca_1;
      for (std::size_t i=0; i<4; i++) {
	 mmdb::Atom *at = v[indices[i]];
	 if (at) {
	    clipper::Coord_orth at_pos = co(v[indices[i]]);
	    clipper::Coord_orth at_pos_new = util::rotate_around_vector(dir, at_pos, p_ca_1, ang);
	    at->x = at_pos_new.x();
	    at->y = at_pos_new.y();
	    at->z = at_pos_new.z();
	 }
      }
   }
}

// can throw a runtime_error
//
// we need a rama_table_set to get the residue types
// (maybe that function should be a static of a rama_table_set)
//
coot::triple_crankshaft_set::triple_crankshaft_set(const residue_spec_t &spec_first_residue,
						   const zo::rama_table_set &zorts,
						   mmdb::Manager *mol) {

   std::pair<mmdb::Residue *, mmdb::Residue *> rs =
      util::get_this_and_next_residues(spec_first_residue, mol);
   if (rs.first) {
      if (rs.second) {
	 mmdb::Residue *res_1 = rs.first;
	 mmdb::Residue *res_2 = rs.second;

	 // when we crankshaft the C, O and N of res_1(n) and res_2(n+1),
	 // we need to know the C of residue n-1 and the N of residue n+2
	 // because the phi,psi of n and n+1 depend on them.
	 mmdb::Residue *res_0 = util::get_previous_residue(spec_first_residue, mol);
	 if (res_0) {
	    // res_3 is needed for phi,psi of res_2
	    residue_spec_t spec_2(res_2);
	    mmdb::Residue *res_3 = util::get_following_residue(spec_2, mol);
	    if (res_3) {
	       residue_spec_t spec_3(res_3);
	       mmdb::Residue *res_4 = util::get_following_residue(spec_3, mol);
	       if (res_4) {
		  residue_spec_t spec_4(res_4);
		  mmdb::Residue *res_5 = util::get_following_residue(spec_4, mol);
		  if (res_5) {
		     std::string this_residue_name = res_1->GetResName();
		     std::string next_residue_name = res_2->GetResName();
		     std::string res_name_3 = res_3->GetResName();
		     std::string res_name_4 = res_4->GetResName();
		     std::string res_name_5 = res_5->GetResName();
		     // I think that get_residue_type() can be a static
		     std::string rt1 = zorts.get_residue_type(this_residue_name, next_residue_name);
		     std::string rt2 = zorts.get_residue_type(next_residue_name, res_name_3);
		     std::string rt3 = zorts.get_residue_type(res_name_3, res_name_4);
		     std::string rt4 = zorts.get_residue_type(res_name_4, res_name_5);
		     residue_types.resize(5); // the type of 0th residue not relevant
		     residue_types[1] = rt1;
		     residue_types[2] = rt2;
		     residue_types[3] = rt3;
		     residue_types[4] = rt4;
		     if (false)
			std::cout << "making cs 0 from "
				  << residue_spec_t(res_0) << " "
				  << residue_spec_t(res_1) << " "
				  << residue_spec_t(res_2) << " "
				  << residue_spec_t(res_3) << " "
				  << std::endl;
		     if (false)
			std::cout << "making cs 1 from "
				  << residue_spec_t(res_1) << " "
				  << residue_spec_t(res_2) << " "
				  << residue_spec_t(res_3) << " "
				  << residue_spec_t(res_4) << " "
				  << std::endl;
		     if (false)
			std::cout << "making cs 2 from "
				  << residue_spec_t(res_2) << " "
				  << residue_spec_t(res_3) << " "
				  << residue_spec_t(res_4) << " "
				  << residue_spec_t(res_5) << " "
				  << std::endl;
		     cs[0] = crankshaft_set(res_0, res_1, res_2, res_3);
		     cs[1] = crankshaft_set(res_1, res_2, res_3, res_4);
		     cs[2] = crankshaft_set(res_2, res_3, res_4, res_5);
		  }
	       }
	    }
	 }
      }
   }
}


void
coot::triple_crankshaft_set::move_the_atoms(float best_angles[]) {

   std::cout << "move the atoms with peptide rotations "
	     << clipper::Util::rad2d(best_angles[0]) << " "
	     << clipper::Util::rad2d(best_angles[1]) << " "
	     << clipper::Util::rad2d(best_angles[2]) << " "
	     << std::endl;

   for (unsigned int i=0; i<3; i++)
      cs[i].move_the_atoms(best_angles[i]);

}

//
// not const because we can change the atoms of mol if apply_best_angles_flag is set
void
coot::crankshaft::triple_spin_search(const residue_spec_t &spec_first_residue,
				     const zo::rama_table_set &zorts,
				     bool apply_best_angles_flag,
				     int n_samples) {

   // 0.1-2-3-4.5  (. for ref, - is spin)

   float div = 1.0/float(n_samples);

   try {
      triple_crankshaft_set tcs(spec_first_residue, zorts, mol);

      float best_angles[3];
      best_angles[0] = -10; best_angles[1] = -10; best_angles[2] = -10;
      float best_prob = 0;

      for (int i=0; i<n_samples; i++) {
	 float a0 = float(i) * div * 2*M_PI;
	 phi_psi_t pp_1 = tcs.phi_psi(0, a0);
	 float pr_1 = zorts.value(pp_1, tcs.residue_type(1));
	 for (int j=0; j<n_samples; j++) {
	    float a1 = float(j) * div * 2*M_PI;
	    phi_psi_t pp_2 = tcs.phi_psi(1, a1);
	    float pr_2 = zorts.value(pp_2, tcs.residue_type(1));
	    for (int k=0; k<n_samples; k++) {
	       float a2 = float(k) * div * 2*M_PI;
	       // the last peptides gives us the probabilities for 2 phi,psi pairs
	       std::pair<phi_psi_t, phi_psi_t> pp_3 = tcs.phi_psis_last(a2);
	       std::pair<float,float> pr_3_pair =
		  probability_of_spin_orientation(pp_3, tcs.residue_type(3), tcs.residue_type(4), zorts);
	       float log_prob = pr_1 + pr_2 + pr_3_pair.first + pr_3_pair.second;
	       if (log_prob > best_prob) {
		  best_prob = log_prob;
		  best_angles[0] = a0;
		  best_angles[1] = a1;
		  best_angles[2] = a2;
	       }
	    }
	 }
      }

      std::cout << "INFO::  best log prob: " << best_prob << "  angles: "
		<< clipper::Util::rad2d(best_angles[0]) << " "
		<< clipper::Util::rad2d(best_angles[1]) << " "
		<< clipper::Util::rad2d(best_angles[2]) << " "
		<< std::endl;

      if (apply_best_angles_flag)
	 tcs.move_the_atoms(best_angles);
   }

   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
}

// not const because we can change the atoms of mol if apply_best_angles_flag is set
//
// change the return type
void
coot::crankshaft::find_maxima(const residue_spec_t &spec_first_residue,
			      const zo::rama_table_set &zorts,
			      unsigned int n_samples) {

   n_samples = 1; // hack for testing

   // 0.1-2-3-4.5  (. for ref, - is spin)

   try {
      triple_crankshaft_set tcs(spec_first_residue, zorts, mol);

      float div = 1.0/float(n_samples);
      for (std::size_t i_sample=0; i_sample<n_samples; i_sample++) {

	 // make these random numbers 0 -> 2pi
	 float start_angles[] = { 10.0, 20.0, 30.0 }; // in radians
	 run_optimizer(start_angles, tcs, zorts);
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }

}

// The function and the derivatives need to be turned upside down at some
// stage, because actually, we want to maximize, not minimize.
//
// static
double
coot::crankshaft::optimize_a_triple::f(const gsl_vector *v, void *params) {

   // convert from spin angles (in v) to phi,psi and then log probability

   // extract the table refs
   const param_holder_t *param_holder = reinterpret_cast<param_holder_t *> (params);
   const zo::rama_table_set &zorts = param_holder->zorts;

   if (false)
      std::cout << "in f() with tcs has residue types "
		<< param_holder->tcs.residue_type(1) << " "
		<< param_holder->tcs.residue_type(2) << std::endl;

   float log_prob_sum = 0;
   for (unsigned int i=0; i<3; i++) { // i is cs index, not residue_type index
      phi_psi_t pp = param_holder->tcs[i].phi_psi(gsl_vector_get(v, i));
      float pr = param_holder->tcs.log_prob(pp, i+1, zorts); // index i+1 is for residue type
      log_prob_sum += pr;
   }

   std::cout << "f() with angles "
	     << gsl_vector_get(v, 0) << " " <<  gsl_vector_get(v, 1) << " " << gsl_vector_get(v, 2) << " "
	     << " returns " << log_prob_sum << std::endl;
   return log_prob_sum;
}


/* The gradient of f, df/da */
void
coot::crankshaft::optimize_a_triple::df(const gsl_vector *v,
					void *params,
					gsl_vector *df_vec) {

   // convert from spin angles (in v) to phi,psi and then log probability


   // When it comes to refinement, moving the atoms:
   //
   // L = log_prob_sum
   // a is peptide rotation angle
   //
   // dL/da = dphi/da  * dL/dphi + dpsi/da * dL/dpsi
   //
   // zorts df gives us dL/dphi and dL/dpsi
   //
   // dphi/da = dphi/dx * dx/da + dphi/dy * dy/da + dphi/dz * dz/da
   //
   // torsion restraint derivatives give us dphi/d{x,y,z}
   //
   // what is dx/da?
   //
   // dx/da = rcos(a), dy/da = -rsin(a), dz/da = 0 in the coordinate system
   // of the peptide - where CA-CA is along (0,0,1).
   // X,Y plane orientation is arbitrary (use the position of the N atom?)
   // You will need to work out r.
   // You will need to rotate this into the world coordinates frame.

   bool do_numerical = true; // debugging

   // extract the table refs
   const param_holder_t *param_holder = reinterpret_cast<param_holder_t *> (params);
   const zo::rama_table_set &zorts = param_holder->zorts;
   for (unsigned int i=0; i<3; i++) {
      phi_psi_t pp = param_holder->tcs[i].phi_psi(gsl_vector_get(v, i));
      const std::string &rt = param_holder->tcs.residue_type(i+1);
      std::cout << "debug:: in df() for peptide index " << i << " rt " << rt << std::endl;
      std::pair<mmdb::realtype,mmdb::realtype> grads = zorts.df(pp, rt);
      std::cout << "debug:: in df() for peptide index " << i << " pp " << pp.phi << " " << pp.psi
		<< " gradients " << grads.first << " " << grads.second << std::endl;
      gsl_vector_set(df_vec, i, grads.first);

      if (do_numerical) {
	 float delta = 0.01;
	 phi_psi_t pp_1 = param_holder->tcs[i].phi_psi(gsl_vector_get(v, i) + delta);
	 phi_psi_t pp_2 = param_holder->tcs[i].phi_psi(gsl_vector_get(v, i) - delta);
	 float pr_1 = param_holder->tcs.log_prob(pp_1, i+1, zorts);
	 float pr_2 = param_holder->tcs.log_prob(pp_2, i+1, zorts);
	 float numerical_grad = (pr_1 - pr_2) / (2 * delta);
	 std::cout << "debug:: in df() for peptide index " << i << " pp "
		   << pp.phi << " " << pp.psi
		   << " gradients " << grads.first << " " << grads.second
		   << " numerical " << numerical_grad << std::endl;
	 gsl_vector_set(df_vec, i, numerical_grad);
      }
   }
}

/* Compute both f and df together. */
void
coot::crankshaft::optimize_a_triple::fdf(const gsl_vector *x, void *params,
					    double *f_in, gsl_vector *df_in) {
  *f_in = f(x, params);
  df(x, params, df_in);
}

void
coot::crankshaft::run_optimizer(float start_angles[],
				const coot::triple_crankshaft_set &tcs,
				const zo::rama_table_set &zorts) {

  size_t iter = 0;
  int status;

  /* Position of the minimum (1,2), scale factors 10,20, height 30. */
  // double par[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };

  std::vector<std::string> residue_types; // where does this come from?
  param_holder_t param_holder(zorts, tcs);

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = 3;
  my_func.f = &optimize_a_triple::f;
  my_func.df = &optimize_a_triple::df;
  my_func.fdf = &optimize_a_triple::fdf;
  my_func.params = reinterpret_cast<void *> (&param_holder);

  /* Starting point, x = (5,7) */

  x = gsl_vector_alloc(3);
  gsl_vector_set(x, 0, start_angles[0]); // refactor
  gsl_vector_set(x, 1, start_angles[1]);
  gsl_vector_set(x, 2, start_angles[2]);

  const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_pr;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, 3);

  gsl_multimin_fdfminimizer_set(s, &my_func, x, 2.0, 1.0); // within 1 degree is good enough

  do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        break;

      status = gsl_multimin_test_gradient(s->gradient, 1e-3);

      if (status == GSL_SUCCESS)
        printf ("Minimum found at:\n");

      std::cout << iter << " "
		<< gsl_vector_get(s->x, 0) << " "
		<< gsl_vector_get(s->x, 1) << " "
		<< gsl_vector_get(s->x, 2) << " "
		<< s->f << "\n";

  } while (status == GSL_CONTINUE && iter < 100);

  if (status == GSL_ENOPROG) {
     std::cout << "No progress" << std::endl;
  }

  std::cout << "iter: " << iter << std::endl;

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

}

mmdb::Atom *
coot::crankshaft::get_atom(mmdb::Residue *res_1, const std::string &atom_name_in) const {

   mmdb::Atom *r = 0;
   mmdb::Atom **residue_atoms_1 = 0;
   int n_residue_atoms_1;
   res_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   for (int iat=0; iat<n_residue_atoms_1; iat++) {
      mmdb::Atom *at = residue_atoms_1[iat];
      std::string atom_name = at->name;
      if (atom_name == atom_name_in) {
	 r = at;
	 break;
      }
   }
   return r;
}

coot::phi_psi_t
coot::crankshaft_set::phi_psi(const clipper::Coord_orth &C_pos,
			      const clipper::Coord_orth &N_pos) const {

   clipper::Coord_orth p0 = co(v[0]);
   clipper::Coord_orth p1 = co(v[1]);
   clipper::Coord_orth p_ca_1 = co(ca_1);
   clipper::Coord_orth p_ca_2 = co(ca_2);
   double torsion_phi_1 = clipper::Coord_orth::torsion(p0, p1, p_ca_1, C_pos);
   double torsion_psi_1 = clipper::Coord_orth::torsion(p1, p_ca_1, C_pos, N_pos);

   return phi_psi_t(torsion_phi_1, torsion_psi_1);
}

coot::phi_psi_t
coot::crankshaft_set::phi_psi(float ang) const {

   if (false) {
      std::cout << "here in phi_psi with ang " << ang << std::endl;
      std::cout << "here in phi_psi with ca_1 " << ca_1 << " " << atom_spec_t(ca_1) << std::endl;
      std::cout << "here in phi_psi with ca_2 " << ca_2 << " " << atom_spec_t(ca_2) << std::endl;
   }

   if (! ca_1)
      std::cout << "ERROR:: ca_1 is unset in crankshaft_set::phi_psi()" << std::endl;
   if (! ca_2)
      std::cout << "ERROR:: ca_2 is unset in crankshaft_set::phi_psi()" << std::endl;

   clipper::Coord_orth p_ca_1 = co(ca_1);
   clipper::Coord_orth p_ca_2 = co(ca_2);

   clipper::Coord_orth dir = p_ca_2 - p_ca_1;
   clipper::Coord_orth C_pos = co(v[2]);
   clipper::Coord_orth N_pos = co(v[4]);

   clipper::Coord_orth C_new = util::rotate_around_vector(dir, C_pos, p_ca_1, ang);
   clipper::Coord_orth N_new = util::rotate_around_vector(dir, N_pos, p_ca_1, ang);

   return phi_psi(C_new, N_new);
   
}



// both phi,psi pairs - for the 2 CAs involved
//
std::pair<coot::phi_psi_t, coot::phi_psi_t>
coot::crankshaft_set::get_phi_phis(const clipper::Coord_orth &C_pos,
				   const clipper::Coord_orth &N_pos) const {
   // atom-idx-2 is res-1-C
   // atom-idx-4 is res-2-N

   // phi-1: 0 1 CA_1 2
   // psi-1: 1 CA_1 2 4
   // phi-2: 2 4 CA_2 6
   // psi-2: 4 CA_2 6 7

   clipper::Coord_orth p0 = co(v[0]);
   clipper::Coord_orth p1 = co(v[1]);
   clipper::Coord_orth p2 = co(v[2]);
   clipper::Coord_orth p3 = co(v[3]);
   clipper::Coord_orth p4 = co(v[4]);
   // clipper::Coord_orth p5 = co(v[5]); might be null
   clipper::Coord_orth p6 = co(v[6]);
   clipper::Coord_orth p7 = co(v[7]);

   clipper::Coord_orth p_ca_1 = co(ca_1);
   clipper::Coord_orth p_ca_2 = co(ca_2);

   clipper::Coord_orth dir = p_ca_2 - p_ca_1;

   double torsion_phi_1 = clipper::Coord_orth::torsion(p0, p1, p_ca_1, C_pos);
   double torsion_psi_1 = clipper::Coord_orth::torsion(p1, p_ca_1, C_pos, N_pos);
   double torsion_phi_2 = clipper::Coord_orth::torsion(C_pos, N_pos, p_ca_2, p6);
   double torsion_psi_2 = clipper::Coord_orth::torsion(N_pos, p_ca_2, p6, p7);

   phi_psi_t pp1(torsion_phi_1, torsion_psi_1, 0);
   phi_psi_t pp2(torsion_phi_2, torsion_psi_2, 0);

   return std::pair<coot::phi_psi_t, coot::phi_psi_t> (pp1, pp2);
}


// what are the phi,psis when rotated by ang radians?
//
// ang in radians
std::pair<coot::phi_psi_t, coot::phi_psi_t>
coot::crankshaft_set::phi_psis(float ang) const {

   // phi-1: 0 1 CA_1 2
   // psi-1: 1 CA_1 2 4
   // phi-2: 2 4 CA_2 6
   // psi-2: 4 CA_2 6 7

   // return what the torsion would be if it has been rotated around CA-1 CA-2 by ang radian
   //

   clipper::Coord_orth p_ca_1 = co(ca_1);
   clipper::Coord_orth p_ca_2 = co(ca_2);

   clipper::Coord_orth dir = p_ca_2 - p_ca_1;
   clipper::Coord_orth C_pos = co(v[2]);
   clipper::Coord_orth N_pos = co(v[4]);

   clipper::Coord_orth C_new = util::rotate_around_vector(dir, C_pos, p_ca_1, ang);
   clipper::Coord_orth N_new = util::rotate_around_vector(dir, N_pos, p_ca_1, ang);

   return get_phi_phis(C_new, N_new);

}



// caller ensures that res_1 and res_2 are valid
std::vector<mmdb::Atom *>
coot::crankshaft::get_mainchain_atoms(mmdb::Residue *res_1, mmdb::Residue *res_2) const {

   std::vector<mmdb::Atom *> v;
   mmdb::Atom **residue_atoms_1 = 0;
   int n_residue_atoms_1;
   res_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   mmdb::Atom **residue_atoms_2 = 0;
   int n_residue_atoms_2;
   res_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

   if (n_residue_atoms_1 > 0) {
      if (n_residue_atoms_2 > 0) {
	 for (int iat=0; iat<n_residue_atoms_1; iat++) {
	    mmdb::Atom *at = residue_atoms_1[iat];
	    if (is_main_chain_p(at)) {
	       v.push_back(at);
	    }
	 }
      }
   }
   return v;
}

void
coot::crankshaft::test() const {

   zo::rama_table_set zorts;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 std::cout << "chain " << chain_p << std::endl;
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) {
	    mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	    residue_spec_t rs(residue_p);
	    std::cout << "residue " << rs << std::endl;
	    std::vector<std::pair<float, float> > r = spin_search(rs, zorts);
	    if (r.size()) {
	       std::cout << "Residue " << rs << std::endl;
	       for (std::size_t i=0; i<r.size(); i++) {
		  std::cout << i << "   " << r[i].first << " " << r[i].second << std::endl;
	       }
	    }
	 }
      }
   }
}

float
coot::crankshaft::probability_of_spin_orientation(const coot::phi_psi_t&pp,
						  const std::string &res_table_type_1,
						  const zo::rama_table_set &zorts) const {

   return zorts.value(pp, res_table_type_1);
}


std::pair<float, float>
coot::crankshaft::probability_of_spin_orientation(const std::pair<phi_psi_t, phi_psi_t> &ppp,
						  const std::string &res_table_type_1,
						  const std::string &res_table_type_2,
						  const zo::rama_table_set &zorts) const {

   return std::pair<float, float> (zorts.value(res_table_type_1, ppp.first.phi,  ppp.first.psi),
				   zorts.value(res_table_type_2, ppp.second.phi, ppp.second.psi));
   
}

coot::phi_psi_t
coot::triple_crankshaft_set::phi_psi(unsigned int peptide_idx, float angle) const {
   return cs[peptide_idx].phi_psi(angle);
}

std::pair<coot::phi_psi_t, coot::phi_psi_t>
coot::triple_crankshaft_set::phi_psis_last(float ang_third) const {
   return cs[2].phi_psis(ang_third);
}



coot::triple_crankshaft_set::triple_crankshaft_set(mmdb::Residue *res_0,
						   mmdb::Residue *res_1,
						   mmdb::Residue *res_2,
						   mmdb::Residue *res_3,
						   mmdb::Residue *res_4,
						   mmdb::Residue *res_5,
						   const std::vector<std::string> &residue_types_in) {

   std::cout << "constructing a crankshaft_set for 0" << std::endl;
   cs[0] = crankshaft_set(res_0, res_1, res_2, res_3);
   std::cout << "constructing a crankshaft_set for 1" << std::endl;
   cs[1] = crankshaft_set(res_1, res_2, res_3, res_4);
   std::cout << "constructing a crankshaft_set for 2" << std::endl;
   cs[2] = crankshaft_set(res_2, res_3, res_4, res_5);
   residue_types = residue_types_in;
   
}

