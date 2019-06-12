

#include "helix-analysis.hh"


void
coot::helix_params_container_t::make(mmdb::Manager *mol_in, const std::string chain_id,
				     int resno_helix_start, int resno_helix_end) {

   mol = mol_in;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      if (chain_id == chain_p->GetChainID()) { 
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 mmdb::Atom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int resno = residue_p->GetSeqNum();

	    if (resno >= resno_helix_start) {
	       if (resno <= (resno_helix_end - 2)) { // that was a bit hacky!
		  std::cout << "chain_id: " << residue_p->GetChainID()
			    << " resno: " << resno << std::endl;
		  atom_quad quad = get_quad(" CA ", chain_p, ires);

		  if (quad.filled_p()) {
		     double tors = quad.torsion();
		     helix_params_t t(resno, quad, tors);
		  }
	       }
	    }
	 }
      }
   }
}


// get the set of atoms starting from the given residue serial number.
coot::atom_quad
coot::helix_params_container_t::get_quad(const std::string &atom_name,
					 mmdb::Chain *chain_p,
					 int res_serial_no) {

   atom_quad quad;
   int nres = chain_p->GetNumberOfResidues();
   for (int ires=0; ires<4; ires++) {
      int isr = res_serial_no + ires;
      if (isr < nres) {
	 mmdb::Residue *res_p = chain_p->GetResidue(isr);
	 if (res_p) {
	    mmdb::Atom *at = res_p->GetAtom(atom_name.c_str());
	    if (at) {
	       if (ires == 0) quad.atom_1 = at;
	       if (ires == 1) quad.atom_2 = at;
	       if (ires == 2) quad.atom_3 = at;
	       if (ires == 3) quad.atom_4 = at;
	    } 
	 } 
      } 
   }
   return quad;

} 

void
coot::helix_params_t::calc_A() {

   // Atom-3 is the first one we can calculate a internal torsion for (i.e. tau_23).
   // 
   // phi is the angle at Atom-3 (2-3-4).
   // 
   try {
      double phi_3  = clipper::Util::d2rad(quad.angle_3());
      double tau_23 = clipper::Util::d2rad(quad.torsion());
      double cp = cos(phi_3);
      double sp = sin(phi_3);
      double ct = cos(tau_23);
      double st = sin(tau_23);
      A = clipper::Mat33<double> (-cp,  -sp, 0,
				  sp*ct, -cp*ct, -st,
				  sp*st, -cp*st, ct);
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING::" << rte.what() << std::endl;
   } 

} 

void
coot::helix_params_t::calc_B() {

   if (quad.atom_2 && quad.atom_3) { 
      clipper::Coord_orth pt_2(quad.atom_2->x, quad.atom_2->y, quad.atom_2->z);
      clipper::Coord_orth pt_3(quad.atom_3->x, quad.atom_3->y, quad.atom_3->z);
      double d = clipper::Coord_orth::length(pt_2, pt_3);
      B = clipper::Coord_orth(d, 0, 0);
      clipper::RTop_orth A_rtop(A, clipper::Coord_orth(0,0,0));
      clipper::RTop_orth A_transpose_rtop(A.transpose(), clipper::Coord_orth(0,0,0));

      clipper::Coord_orth B_pr  = B.transform(A_transpose_rtop);
      clipper::Coord_orth B_dpr = B.transform(A_rtop);

      clipper::Coord_orth C    = B_pr - B;
      clipper::Coord_orth C_pr = B - B_dpr;
      double len_C = sqrt(C.lengthsq());
      double cos_theta = clipper::Coord_orth::dot(C,C_pr)/(len_C*len_C);
      double theta = acos(cos_theta);
      double sin_theta = sin(theta);
      std::cout << "theta: " << clipper::Util::rad2d(theta) << " degrees "
		<< " cos(theta): " << cos_theta << " "
		<< " sin(theta): " << sin_theta << std::endl;

      clipper::Coord_orth e_xi(C.unit());
      clipper::Coord_orth c_crossed(clipper::Coord_orth::cross(C, C_pr));

      // Let's try another method to get to sin theta:
      // d sin(theta) = B.(CxC')/C^2
      // 
      // double d_sin_theta = clipper::Coord_orth::dot(B, c_crossed)/(len_C*len_C);
      // Oh, what is d?  Fail...
      
      
      clipper::Coord_orth e_zeta(c_crossed.x() * sin_theta/(len_C*len_C),
				 c_crossed.y() * sin_theta/(len_C*len_C),
				 c_crossed.z() * sin_theta/(len_C*len_C));
      clipper::Coord_orth e_nu(clipper::Coord_orth::cross(e_zeta, e_xi));

       std::cout << "C: " << C.format() << " "
		 << "B_pr " << B_pr.format() << " "
		 << "B " << B.format()
		 << std::endl;

      std::cout << "e_nu: "   << e_nu.format()   << " "
		<< "e_xi: "   << e_xi.format()   << " "
		<< "e_zeta: " << e_zeta.format() << " "
		<< std::endl;

      std::cout << "e lengths: "
		<< " e_nu: "   << sqrt(e_nu.lengthsq()) << "  "
		<< " e_xi: "   << sqrt(e_xi.lengthsq()) << "  "
		<< " e_zeta: " << sqrt(e_zeta.lengthsq()) << std::endl;

      clipper::Coord_orth a_bits(A(2,1) - A(1,2),
				 A(0,2) - A(2,0),
				 A(1,0) - A(0,1));
      std::cout << "abits: " << a_bits.format() << " length: " << sqrt(a_bits.lengthsq())
		<< std::endl;

      clipper::Coord_orth a_bits_2(A(1,2) - A(2,1),
				   A(2,0) - A(0,2),
				   A(0,1) - A(1,0));
      std::cout << "abits: " << a_bits_2.format() << " length: " << sqrt(a_bits_2.lengthsq())
		<< std::endl;

      // half sin theta
      double hst = 0.5 * sin (theta);
      clipper::Coord_orth e_zeta_2(a_bits.x() * hst,
				   a_bits.y() * hst,
				   a_bits.z() * hst);

      std::cout << "e_zeta_2: " << e_zeta_2.format() << " length: "
		<< sqrt(e_zeta_2.lengthsq()) << std::endl;
   }
}
