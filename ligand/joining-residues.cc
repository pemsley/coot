
#include "residue_by_phi_psi.hh"


// psi_conditional is meaningful when offset == 1
// phi_conditional is meaningful when offset == -1 (otherwise not, of course)
//
coot::minimol::fragment
coot::residue_by_phi_psi::make_2_res_joining_frag_new(const std::string &chain_id,
						      const connecting_atoms_t &current_atom_positions,
						      const double &phi_conditional,
						      const double &psi_conditional,
						      const phi_psi_t &pp1,
						      const phi_psi_t &pp2,
						      int seqnum,
						      int offset // + or - 1
						      ) const {

   if (offset == 1) {
      return make_2_res_joining_frag_new_building_forwards(chain_id, current_atom_positions,
							   psi_conditional,
							   pp1, pp2, seqnum);
   } else {
      return make_2_res_joining_frag_new_building_backwards(chain_id, current_atom_positions,
							    phi_conditional,
							    pp1, pp2, seqnum);
   }
}

coot::minimol::fragment
coot::residue_by_phi_psi::make_2_res_joining_frag_new_building_forwards(const std::string &chain_id,
									const connecting_atoms_t &current_res_pos_in,
									const double &psi_conditional,
									const phi_psi_t &pp1,
									const phi_psi_t &pp2,
									int seqnum) const {
   bool debug = false;
   minimol::fragment frag(chain_id);

   // we need edit the postions of the reference atoms to add a bit of (useful) jitter
   connecting_atoms_t current_res_pos = current_res_pos_in;

   // refs!
   clipper::Coord_orth &current_n  = current_res_pos.N_pos;
   clipper::Coord_orth &current_ca = current_res_pos.CA_pos;

   // add a bit of jitter
   double rand_lim = 0.25; // 0.1 -> 6.2, 0.3 -> 6.4 ; 0.5 -> 5.9; 0.25 -> 6.4 // I wonder what these numbers are.
   current_n  += clipper::Coord_orth(rand_lim * (util::random()/float (RAND_MAX) - 0.5),
				     rand_lim * (util::random()/float (RAND_MAX) - 0.5),
				     rand_lim * (util::random()/float (RAND_MAX) - 0.5));
   current_ca += clipper::Coord_orth(rand_lim * (util::random()/float (RAND_MAX) - 0.5),
				     rand_lim * (util::random()/float (RAND_MAX) - 0.5),
				     rand_lim * (util::random()/float (RAND_MAX) - 0.5));

   // in radians
   std::pair<bool, double> phi_current = current_res_pos.get_phi();

   // Check if we have 2 residues (this one and the one upstream) on our starting residue, which
   // gives us a phi, from which we can make a condtional psi.
   // 
   // If we don't have 2 residues, then do a normal selection of random phi,psi
   
   if (phi_current.first) {

      minimol::residue res1 = construct_next_res_from_rama_angles(pp1.phi, psi_conditional, pp1.tau,
								  seqnum + 1, current_res_pos);

      connecting_atoms_t just_built_res(res1[" N  "].pos, res1[" CA "].pos, res1[" C  "].pos);
      just_built_res.set_upstream_C(current_res_pos.C_pos);

      minimol::residue res2 = construct_next_res_from_rama_angles(pp2.phi, pp1.psi, pp2.tau,
								  seqnum + 2, just_built_res);

      // now set set the occupancy of res2 to 0.5 or so, because we care more that the
      // first residue is in density
      for (unsigned int iat=0; iat<res2.atoms.size(); iat++)
	 res2.atoms[iat].occupancy = 0.5;

      // now let's correct the position of the O atom of residue 1,
      // now that we have residue 2.
      //
      double angle   = clipper::Util::d2rad(123.0); // N-C-O
      double torsion = clipper::Util::d2rad(0.0);
      clipper::Coord_orth o_pos(res2[" CA "].pos, res2[" N  "].pos, res1[" C  "].pos, 1.231, angle, torsion);
      for (unsigned int iat=0; iat<res1.atoms.size(); iat++)
	 if (res1[iat].name == " O  ")
	    res1[iat].pos = o_pos;

      try {
	 frag.addresidue(res1, 0);
	 frag.addresidue(res2, 0);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: make_2_res_joining_frag_new_building_forwards() "
		   << rte.what() << std::endl;
      }

   } else {

      // singleton

      phi_psi_t pp0 = get_phi_psi_by_random(rama, rama_max, false); // make one up for the starting residue
      minimol::residue res1 = construct_next_res_from_rama_angles(pp1.phi, pp0.psi, pp1.tau,
								  seqnum + 1, current_res_pos);
      connecting_atoms_t just_built_res(res1[" N  "].pos, res1[" CA "].pos, res1[" C  "].pos);
      just_built_res.set_upstream_C(current_res_pos.C_pos);
      minimol::residue res2 = construct_next_res_from_rama_angles(pp2.phi, pp1.psi, pp2.tau,
								  seqnum + 2, just_built_res);
      for (unsigned int iat=0; iat<res2.atoms.size(); iat++)
	 res2.atoms[iat].occupancy = 0.5;
      // now let's correct the position of the O atom of residue 1,
      // now that we have residue 2.
      //
      double angle   = clipper::Util::d2rad(123.0); // N-C-O
      double torsion = clipper::Util::d2rad(0.0);
      clipper::Coord_orth o_pos(res2[" CA "].pos, res2[" N  "].pos, res1[" C  "].pos, 1.231, angle, torsion);
      for (unsigned int iat=0; iat<res1.atoms.size(); iat++)
	 if (res1[iat].name == " O  ")
	    res1[iat].pos = o_pos;
      try {
	 frag.addresidue(res1, 0);
	 frag.addresidue(res2, 0);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: make_2_res_joining_frag_new_building_forwards() "
		   << rte.what() << std::endl;
      }
   }

   return frag;
}

coot::minimol::fragment
coot::residue_by_phi_psi::make_2_res_joining_frag_new_building_backwards(const std::string &chain_id,
									 const connecting_atoms_t &current_res_pos,
									 const double &phi_conditional,
									 const phi_psi_t &pp1,
									 const phi_psi_t &pp2,
									 int seqnum) const {

   bool debug = false;
   minimol::fragment frag(chain_id);

   clipper::Coord_orth current_n  = current_res_pos.N_pos;
   clipper::Coord_orth current_ca = current_res_pos.CA_pos;

   // in radians
   std::pair<bool, double> psi_current = current_res_pos.get_psi();

   if (psi_current.first) {

      // add a bit of jitter
      double rand_lim = 0.25;
      current_n  += clipper::Coord_orth(rand_lim * (util::random()/float (RAND_MAX) - 0.5),
					rand_lim * (util::random()/float (RAND_MAX) - 0.5),
					rand_lim * (util::random()/float (RAND_MAX) - 0.5));
      current_ca += clipper::Coord_orth(rand_lim * (util::random()/float (RAND_MAX) - 0.5),
					rand_lim * (util::random()/float (RAND_MAX) - 0.5),
					rand_lim * (util::random()/float (RAND_MAX) - 0.5));
      
      double phi_conditional_deg = clipper::Util::rad2d(get_phi_by_random_given_psi(psi_current.second, rama));

      minimol::residue res1 = construct_prev_res_from_rama_angles(phi_conditional_deg, pp1.psi, pp1.tau,
								  seqnum - 1, current_res_pos);

      connecting_atoms_t just_built_res(res1[" N  "].pos, res1[" CA "].pos, res1[" C  "].pos);
      just_built_res.set_downstream_N(current_res_pos.N_pos);
      minimol::residue res2 = construct_prev_res_from_rama_angles(pp1.phi, pp2.psi, pp2.tau,
								  seqnum - 2, just_built_res);
      for (unsigned int iat=0; iat<res2.atoms.size(); iat++)
	 res2.atoms[iat].occupancy = 0.5;
      try {
	 frag.addresidue(res2, 0);
	 frag.addresidue(res1, 0);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: make_2_res_joining_frag_new_building_forwards() "
		   << rte.what() << std::endl;
      }

   } else {

      // std::cout << "------------------ singleton " << std::endl;
      // singleton
      phi_psi_t pp0 = get_phi_psi_by_random(rama, rama_max, false); // make one up for the starting residue
      minimol::residue res1 = construct_prev_res_from_rama_angles(pp0.phi, pp1.psi, pp1.tau,
								  seqnum - 1, current_res_pos);
      connecting_atoms_t just_built_res(res1[" N  "].pos, res1[" CA "].pos, res1[" C  "].pos);
      just_built_res.set_downstream_N(current_res_pos.N_pos);
      minimol::residue res2 = construct_prev_res_from_rama_angles(pp1.phi, pp2.psi, pp2.tau,
								  seqnum - 2, just_built_res);
      for (unsigned int iat=0; iat<res2.atoms.size(); iat++)
	 res2.atoms[iat].occupancy = 0.5;
      try {
	 frag.addresidue(res2, 0);
	 frag.addresidue(res1, 0);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: make_2_res_joining_frag_new_building_forwards() "
		   << rte.what() << std::endl;
      }

   }

   return frag;
}
