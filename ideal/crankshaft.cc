/*
 * ideal/crankshaft.cc
 *
 * Copyright 2017 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#include <algorithm> // for sort
#include <iomanip>

#include "crankshaft.hh"

#include "geometry/main-chain.hh"
#include "coot-utils/coot-coord-utils.hh"

// for refinement (put refinement in new file?)
#include "simple-restraint.hh"
#include "coot-utils/coot-map-utils.hh"
#include "utils/logging.hh"
extern logging logger;

std::ostream &
coot::operator<<(std::ostream &s, const coot::crankshaft::scored_triple_angle_set_t &as) {

   s << "minus-log-prob: " << as.minus_log_prob << " from angles";
   for (std::size_t i=0; i<as.angles.size(); i++)
      s << std::setw(9) << clipper::Util::rad2d(as.angles[i]) << " ";
   return s;
}

std::ostream &
coot::operator<<(std::ostream &s, const coot::crankshaft::scored_nmer_angle_set_t &as) {

   s << "minus_log_prob: " << std::setw(8) << as.minus_log_prob
     << " combi-score: " << std::right << std::setprecision(3) << std::fixed << as.combi_score
     << " from angles";
   for (std::size_t i=0; i<as.angles.size(); i++)
      s << std::setw(9) << clipper::Util::rad2d(as.angles[i]) << " ";
   return s;
}

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

   if (false) {
      std::cout << "in crankshaft_set constructor for res_1 " << residue_spec_t(res_1)
		<< " ca_1 was set to " << ca_1 << std::endl;
      std::cout << "in crankshaft_set constructor for res_1 " << residue_spec_t(res_1)
		<< " ca_2 was set to " << ca_2 << std::endl;
   }

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
      make_trans_from_non_pro_cis_if_needed();
   } else {
      throw(std::runtime_error("missing a mainchain atom"));
   }
}

// move the atoms of of res_1 and res_2 in mol
// Find if this is cis and move the N, (H), C and O to make it more like a trans
// N (and H) will move most of course.
// C and O will move a bit in the direction of new N position
//
void
coot::crankshaft_set::make_trans_from_non_pro_cis_if_needed() {

   if (ca_1) {
      if (ca_2) {
	 const std::string res_name_2 = ca_2->GetResName();
	 if (res_name_2 != "PRO") {
	    if (is_cis()) {
	       mmdb::Atom *N_at = v[4];
	       if (N_at) {
		  clipper::Coord_orth N_pos = co(N_at);
		  clipper::Coord_orth C_pos = co(v[2]);
		  clipper::Coord_orth O_pos = co(v[3]);
		  clipper::Coord_orth p_ca_1 = co(ca_1);
		  clipper::Coord_orth p_ca_2 = co(ca_2);
		  clipper::Coord_orth dir = p_ca_2 - p_ca_1;

		  clipper::Coord_orth N_pos_new = util::rotate_around_vector(dir, N_pos, p_ca_1, M_PI);
		  update_position(N_at, N_pos_new);

		  // Move the N and C closer together (as needed) along the N-C vector.
		  // Ideal distance 1.33A

		  clipper::Coord_orth n_c_vec = C_pos - N_pos_new;
		  clipper::Coord_orth n_c_unit(n_c_vec.unit());
		  double bl = sqrt(n_c_vec.lengthsq());
		  double bl_delta = bl - 1.33;
		  // std::cout << "   cis N-C correction bl_delta " << bl_delta << std::endl;

		  clipper::Coord_orth pos_corr = bl_delta * 0.5 * n_c_unit;

		  N_pos_new += pos_corr;
		  C_pos     -= pos_corr;
		  O_pos     -= pos_corr;

		  update_position(N_at, N_pos_new);
		  update_position(v[2], C_pos);
		  update_position(v[3], O_pos);

		  if (v[5]) { // H, null probably
		     clipper::Coord_orth H_pos = co(v[5]);
		     clipper::Coord_orth at_pos_new = util::rotate_around_vector(dir, H_pos, p_ca_1, M_PI);
		     update_position(v[5], at_pos_new);
		  }
	       }
	    }
	 }
      }
   }
}

bool
coot::crankshaft_set::is_cis() const {

   bool cis = false;
   if (ca_1) {
      if (ca_2) {
	 clipper::Coord_orth N_pos = co(v[4]);
	 clipper::Coord_orth C_pos = co(v[2]);
	 clipper::Coord_orth p_ca_1 = co(ca_1);
	 clipper::Coord_orth p_ca_2 = co(ca_2);
	 clipper::Coord_orth dir = p_ca_2 - p_ca_1;
	 double tors = clipper::Coord_orth::torsion(p_ca_2, N_pos, C_pos, p_ca_1);
	 if (false)
	    std::cout << "   debug:: torsion for ca_1 " << atom_spec_t(ca_1) << " "
		      << clipper::Util::rad2d(tors) << " degrees" << std::endl;
	 if (std::abs(tors) < M_PI_2)
	    cis = true;
      }
   }
   // std::cout << "   debug:: is cis for ca_1 " << atom_spec_t(ca_1) << " " << cis << std::endl;
   return cis;

}


std::vector<mmdb::Residue *>
coot::nmer_crankshaft_set::residues() const {

   std::vector<mmdb::Residue *> v;
   for (std::size_t i=0; i<size(); i++) {
      mmdb::Atom *at_1 = cs_vec[i].ca_1;
      mmdb::Atom *at_2 = cs_vec[i].ca_2;
      mmdb::Residue *r_1 = at_1->residue;
      mmdb::Residue *r_2 = at_2->residue;
      if (std::find(v.begin(), v.end(), r_1) == v.end()) v.push_back(r_1);
      if (std::find(v.begin(), v.end(), r_2) == v.end()) v.push_back(r_2);
   }
   return v;
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

		     // debug
		     // mol->WritePDBASCII("uncis.pdb");
		  }
	       }
	    }
	 }
      }
   }
}

coot::nmer_crankshaft_set::nmer_crankshaft_set(const residue_spec_t &spec_mid_residue,
					       unsigned int n_peptides,
					       const zo::rama_table_set &zorts,
					       mmdb::Manager *mol) {

#ifdef HAVE_CXX11

   // we first need to find what is the first residue (the N-terminal residue of the
   // first rotating peptide). res_0 is the residue before the first residue.

   // This will not do the right thing if there are insertion codes. In that case
   // more cleverness is required.
   int first_from_mid_delta = std::round((n_peptides-1)/2);
   residue_spec_t spec_first_residue(spec_mid_residue.chain_id,
				     spec_mid_residue.res_no - first_from_mid_delta);

   // if n_peptides is 1: then we will need to fill cs[0] with:
   // res_0: residue before the first spec
   // res_1: this residue (which has ca_1) of which we move the C,O
   // res_2: next residue (which has ca_2) or which we move the N
   // res_3: needed for phi,psi of res_2

   // set mess if we fail to find residues of a peptide cranshaft set
   //
   std::string mess; // if this get set, throw an exception.

   mmdb::Residue *residue_p = util::get_residue(spec_first_residue, mol);
   if (residue_p) {
      // when we crankshaft the C, O and N of res_1(n) and res_2(n+1),
      // we need to know the C of residue n-1 and the N of residue n+2
      // because the phi,psi of n and n+1 depend on them.
      mmdb::Residue *res_0 = util::get_previous_residue(spec_first_residue, mol);
      if (res_0) {

	 // a base_residue_spec is res_1 of a crankshaft set (not res_0)
	 //
	 residue_spec_t current_base_residue_spec = spec_first_residue;
	 mmdb::Residue *current_base_residue = residue_p;

	 // maybe it's silly to have "off-by-one" residue type index?
	 residue_types.resize(n_peptides+1);

	 for (std::size_t i_pep=0; i_pep<n_peptides; i_pep++) {
	    mmdb::Residue *pep_res_0 = util::get_previous_residue (current_base_residue_spec, mol);
	    mmdb::Residue *pep_res_2 = util::get_following_residue(current_base_residue_spec, mol);
	    if (pep_res_2) {
	       residue_spec_t current_res_2_spec(pep_res_2);
	       mmdb::Residue *pep_res_3 = util::get_following_residue(current_res_2_spec, mol);
	       if (pep_res_3) {

		  // this can throw a runtime_error
		  crankshaft_set cs(pep_res_0, current_base_residue, pep_res_2, pep_res_3);
		  cs_vec.push_back(cs);

		  // maybe it's silly to have "off-by-one" residue type index?
		  std::string rn_1 = current_base_residue->GetResName();
		  std::string rn_2 = pep_res_2->GetResName();
		  std::string rt = zorts.get_residue_type(rn_1, rn_2);
		  residue_types[i_pep+1] = rt;

		  // setup for next
		  current_base_residue = pep_res_2;
		  current_base_residue_spec = residue_spec_t(pep_res_2);
	       } else {
		  mess = "Fail to get residue 3 for peptide index ";
		  mess += util::int_to_string(i_pep);
		  break;
	       }
	    } else {
	       mess = "Fail to get residue 2 for peptide index ";
	       mess += util::int_to_string(i_pep);
	       break;
	    }
	 }
      } else {
	 mess = "Fail to get residue 0 for starting residue ";
      }
   }

   if (! mess.empty()) throw std::runtime_error(mess);

#endif // HAVE_CXX11

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

      logger.log(log_t::INFO, logging::ltw(" best log prob: "), logging::ltw(best_prob), logging::ltw("  angles: "), logging::ltw(clipper::Util::rad2d(best_angles[0])), logging::ltw(" " + std::to_string(clipper::Util::rad2d(best_angles[1])) + " "), logging::ltw(clipper::Util::rad2d(best_angles[2])));
      // std::cout << "INFO::  best log prob: " << best_prob << "  angles: "
      //           << clipper::Util::rad2d(best_angles[0]) << " "
      //           << clipper::Util::rad2d(best_angles[1]) << " "
      //           << clipper::Util::rad2d(best_angles[2]) << " "
      //           << std::endl;

      if (apply_best_angles_flag)
	 tcs.move_the_atoms(best_angles);
   }

   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
}

// not const because we can change the atoms of mol if apply_best_angles_flag is set
//
std::vector<coot::crankshaft::scored_triple_angle_set_t>
coot::crankshaft::find_maxima_from_triples(const residue_spec_t &spec_first_residue,
					   const zo::rama_table_set &zorts,
					   unsigned int n_samples) {

   // n_samples = 1; // hack for testing

   // 0.1-2-3-4.5  (. for ref, - is spin)

   std::vector<scored_triple_angle_set_t> results;

   try {
      triple_crankshaft_set tcs(spec_first_residue, zorts, mol);

      float div = 1.0/float(n_samples);
      for (std::size_t i_sample=0; i_sample<n_samples; i_sample++) {

	 // make these random numbers 0 -> 2pi
	 float start_angles[] = { 1.0, 2.0, 3.0 }; // in radians

	 for (std::size_t i=0; i<3; i++)
	    start_angles[i] = 2 * M_PI * float(coot::util::random())/float(RAND_MAX);

	 // testing wild solutions
	 // start_angles[0] = clipper::Util::d2rad(264.849);
	 // start_angles[1] = clipper::Util::d2rad(122.957);
	 // start_angles[2] = clipper::Util::d2rad(300.472);
	 //
	 // start_angles[0] = clipper::Util::d2rad();
	 // start_angles[1] = clipper::Util::d2rad();
	 // start_angles[2] = clipper::Util::d2rad();
	 //
	 // start_angles[0] = clipper::Util::d2rad(351.391);
	 // start_angles[1] = clipper::Util::d2rad(287.775);
	 // start_angles[2] = clipper::Util::d2rad(49.7706);

	 scored_triple_angle_set_t as = run_optimizer(start_angles, tcs, zorts);

	 if (as.filled()) {
	    if (false)
	       std::cout << "sample " << i_sample << " " << as.minus_log_prob << " from angles "
			 << std::setw(9) << clipper::Util::rad2d(as.angles[0]) << " "
			 << std::setw(9) << clipper::Util::rad2d(as.angles[1]) << " "
			 << std::setw(9) << clipper::Util::rad2d(as.angles[2]) << std::endl;

	    // add it if not already something similar
	    //
	    bool add_it = true;

	    if (false) { // what's close?
	       float closest_delta = 999.9;
	       for (std::size_t i=0; i<results.size(); i++) {
		  float diff = results[i].angles[0] - as.angles[0];
		  float this_delta_2 = diff * diff;
		  diff = results[i].angles[1] - as.angles[1];
		  this_delta_2 += diff * diff;
		  diff = results[i].angles[2] - as.angles[2];
		  this_delta_2 += diff * diff;
		  float this_delta = sqrt(this_delta_2);
		  if (this_delta < closest_delta) {
		     closest_delta = this_delta;
		  }
	       }
	       std::cout << "closest-delta: " << closest_delta << std::endl;
	    }

	    for (std::size_t i=0; i<results.size(); i++) {
	       if (results[i].is_close(as)) {
		  add_it = false;
		  break;
	       }
	    }
	    if (add_it) {
	       if (false)
		  std::cout << "as: " << as << " from starting angles "
			    << clipper::Util::rad2d(start_angles[0]) << " "
			    << clipper::Util::rad2d(start_angles[1]) << " "
			    << clipper::Util::rad2d(start_angles[2]) << " "
			    << std::endl;
	       results.push_back(as);
	    } else {
	       // std::cout << "we already have that " << std::endl;
	    }
	 }
      }
      std::cout << "From " << n_samples << " samples, found  " << results.size()
		<< " unique solutions" << std::endl;
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }

   std::sort(results.begin(), results.end());

   if (false) {
      std::ofstream f("solutions");
      for (std::size_t i=0; i<results.size(); i++) {
	 f << "   value "
	   << results[i].minus_log_prob << " at (degrees) "
	   << clipper::Util::rad2d(results[i].angles[0]) << " "
	   << clipper::Util::rad2d(results[i].angles[1]) << " "
	   << clipper::Util::rad2d(results[i].angles[2]) << "\n";
      }
   }
   return results;
}

// not const because we can change the atoms of mol if apply_best_angles_flag is set
//
// "scored" because the probability has been calculated (not used in refinement)
//
// only allow solutions that are within log_prob_filter_n_sigma sds of the top solution
// i.e. 2.0 will allow more solutions than 1.0.
//
std::vector<coot::crankshaft::scored_nmer_angle_set_t>
coot::crankshaft::find_maxima(const residue_spec_t &spec_mid_residue,
			      unsigned int n_peptides, // the length of the nmer
			      const zo::rama_table_set &zorts,
			      float log_prob_filter_n_sigma,
			      unsigned int n_samples) {

   // n_samples = 1; // hack for testing

   // 0.1-2-3-4.5  (. for ref, - is spin)

   std::vector<scored_nmer_angle_set_t> results;

#ifdef HAVE_CXX_THREAD

   try {
      // if the atoms are not there, this can throw a std::runtime_error.
      nmer_crankshaft_set cs(spec_mid_residue, n_peptides, zorts, mol);
      if (cs.size() == 0) return results;

      float div = 1.0/float(n_samples);

      unsigned int n_threads = coot::get_max_number_of_threads();

      if (n_threads > 0) {

	 // first split the samples into batches for each thread
	 std::size_t current_thread = 0;
	 std::vector<std::vector<std::size_t> > samples_for_thread(n_threads);
	 for (std::size_t i_sample=0; i_sample<n_samples; i_sample++) {
	    samples_for_thread[current_thread].push_back(i_sample);
	    ++current_thread;
	    if (current_thread == n_threads)
	       current_thread = 0;
	 }
	 results.resize(n_samples);
	 std::vector<std::thread> threads;
	 std::vector<scored_nmer_angle_set_t> sas_vec(n_threads);
	 for (std::size_t i_thread=0; i_thread<n_threads; i_thread++) {
	    if (samples_for_thread[i_thread].size() > 0) {
	       threads.push_back(std::thread(run_optimizer_in_thread,
					     std::cref(samples_for_thread[i_thread]),
					     std::cref(cs), std::cref(zorts), &results));
	    }
	 }
	 for (unsigned int i_thread=0; i_thread<threads.size(); i_thread++)
	    threads.at(i_thread).join();

      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }


   if (results.size() > 0) {

      // first remove failed solutions (no angles)
      // std::cout << "before null-eraser: " << results.size() << std::endl;
      results.erase(std::remove_if(results.begin(), results.end(), null_eraser), results.end());
      // std::cout << "after null-eraser: " << results.size() << std::endl;

      // next remove solutions that are similar to each other
      // pre-sort so that std::unique does the adjacent values comparison correctly
      std::sort(results.begin(), results.end());
      // std::cout << "before similar-eraser: " << results.size() << std::endl;
      results.erase(std::unique(results.begin(), results.end(), scored_solution_comparer), results.end());
      // std::cout << "after  similar-eraser: " << results.size() << std::endl;
      std::sort(results.begin(), results.end()); // resort after std::unique

      // now filter out those with "low" log_prob values
      //
      std::vector<float> v(results.size());
      for (std::size_t i=0; i<results.size(); i++)
	 v[i] = results[i].minus_log_prob;

      float top_score = results[0].minus_log_prob; // highest (most negative)

      // delete those results that are "a long way" below the top_score
      // or are below the mean.
      util::stats_data sd(v);
      float at_least = top_score + log_prob_filter_n_sigma*sd.sd;
      float mean = sd.mean;

      // std::cout << "before log_prob filter we have " << results.size() << " solutions" << std::endl;
      results.erase(std::remove_if(results.begin(), results.end(), eraser(at_least, mean)), results.end());
      // std::cout << "after  log_prob filter we have " << results.size() << " solutions" << std::endl;

      if (false) {
	 std::ofstream f("solutions");
	 for (std::size_t i=0; i<results.size(); i++) {
	    f << "   value "
	      << results[i].minus_log_prob << " at (degrees) ";
	    for (std::size_t iang=0; iang<results[i].angles.size(); iang++)
	       f << clipper::Util::rad2d(results[i].angles[iang]) << " ";
	    f << "\n";
	 }
      }
   }

#endif // HAVE_CXX_THREAD

   return results;
}


// The function and the derivatives need to be turned upside down at some
// stage, because actually, we want to maximize, not minimize.
//
// static
double
coot::crankshaft::optimize_a_triple::f(const gsl_vector *v, void *params) {

   // convert from spin angles (in v) to phi,psi and then log probability

   // extract the table refs
   const triple_set_param_holder_t *param_holder = reinterpret_cast<triple_set_param_holder_t *> (params);
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

   if (false)
      std::cout << "f() with angles "
		<< gsl_vector_get(v, 0) << " " <<  gsl_vector_get(v, 1) << " " << gsl_vector_get(v, 2) << " "
		<< " returns " << log_prob_sum << std::endl;
   return -log_prob_sum; // we want maxima
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
   const triple_set_param_holder_t *param_holder = reinterpret_cast<triple_set_param_holder_t *> (params);
   const zo::rama_table_set &zorts = param_holder->zorts;

   for (unsigned int i=0; i<3; i++) {
      phi_psi_t pp = param_holder->tcs[i].phi_psi(gsl_vector_get(v, i));
      const std::string &rt = param_holder->tcs.residue_type(i+1);
      // std::cout << "debug:: in df() for peptide index " << i << " rt " << rt << std::endl;
      std::pair<mmdb::realtype,mmdb::realtype> grads = zorts.df(pp, rt);
      if (false)
	 std::cout << "debug:: in df() analyticals for peptide index " << i
		   << " pp " << pp.phi << " " << pp.psi
		   << " gradients " << grads.first << " " << grads.second << std::endl;

      if (do_numerical) {
	 float delta = 0.002;
	 float current_val = gsl_vector_get(v, i);
	 phi_psi_t pp_1 = param_holder->tcs[i].phi_psi(current_val + delta);
	 phi_psi_t pp_2 = param_holder->tcs[i].phi_psi(current_val - delta);
	 float pr_1 = param_holder->tcs.log_prob(pp_1, i+1, zorts);
	 float pr_2 = param_holder->tcs.log_prob(pp_2, i+1, zorts);
	 float numerical_grad = - (pr_1 - pr_2) / (2 * delta);
	 if (false)
	    std::cout << "debug:: in df() for peptide index " << i
		      << " angles "
		      << clipper::Util::rad2d(gsl_vector_get(v, i)) << " "
		      << " pp " << pp.phi << " " << pp.psi
		      << " numerical gradient from " << pr_1 << " - " << pr_2 << " / 2 * " << delta
		      << " = " << numerical_grad << std::endl;
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

// static
double
coot::crankshaft::optimize_an_nmer::f(const gsl_vector *v, void *params) {

   // convert from spin angles (in v) to phi,psi and then log probability

   // extract the table refs
   const nmer_set_param_holder_t *param_holder = reinterpret_cast<nmer_set_param_holder_t *> (params);
   const zo::rama_table_set &zorts = param_holder->zorts;

   // debug zorts

   // debug nmer_crankshaft_set &cs

   float log_prob_sum = 0;
   for (unsigned int i=0; i<param_holder->cs.size(); i++) { // i is cs index, not residue_type index
      phi_psi_t pp = param_holder->cs[i].phi_psi(gsl_vector_get(v, i));
      float pr = param_holder->cs.log_prob(pp, i+1, zorts); // index i+1 is for residue type
      log_prob_sum += pr;
   }
   return -log_prob_sum; // we want maxima
}


/* The gradient of f, df/da */
void
coot::crankshaft::optimize_an_nmer::df(const gsl_vector *v,
				       void *params,
				       gsl_vector *df_vec) {

   // convert from spin angles (in v) to phi,psi and then log probability
   // using numerical gradients! :-/

   // extract the table refs
   const nmer_set_param_holder_t *param_holder = reinterpret_cast<nmer_set_param_holder_t *> (params);
   const zo::rama_table_set &zorts = param_holder->zorts;

   for (unsigned int i=0; i<param_holder->cs.size(); i++) {

      if (true) {
	 float delta = 0.002;
	 float current_val = gsl_vector_get(v, i);
	 phi_psi_t pp_1 = param_holder->cs[i].phi_psi(current_val + delta);
	 phi_psi_t pp_2 = param_holder->cs[i].phi_psi(current_val - delta);
	 float pr_1 = param_holder->cs.log_prob(pp_1, i+1, zorts);
	 float pr_2 = param_holder->cs.log_prob(pp_2, i+1, zorts);
	 float numerical_grad = - (pr_1 - pr_2) / (2 * delta);
	 gsl_vector_set(df_vec, i, numerical_grad);
      }
   }
}

/* Compute both f and df together. */
void
coot::crankshaft::optimize_an_nmer::fdf(const gsl_vector *x, void *params,
					double *f_in, gsl_vector *df_in) {
  *f_in = f(x, params);
  df(x, params, df_in);
}

coot::crankshaft::scored_triple_angle_set_t
coot::crankshaft::run_optimizer(float start_angles[],
				const coot::triple_crankshaft_set &tcs,
				const zo::rama_table_set &zorts) {

  size_t iter = 0;
  int status;

  triple_set_param_holder_t param_holder(zorts, tcs);

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = 3;
  my_func.f = &optimize_a_triple::f;
  my_func.df = &optimize_a_triple::df;
  my_func.fdf = &optimize_a_triple::fdf;
  my_func.params = reinterpret_cast<void *> (&param_holder);

  x = gsl_vector_alloc(3);
  gsl_vector_set(x, 0, start_angles[0]); // refactor
  gsl_vector_set(x, 1, start_angles[1]);
  gsl_vector_set(x, 2, start_angles[2]);

  const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_pr;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, 3);

  // small moves Ellie, small moves...
  gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1.0); // within 1 degree is good enough

  do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        break;

      status = gsl_multimin_test_gradient(s->gradient, 5e-2);

      if (status == GSL_SUCCESS)
	 if (false)
	    std::cout << "Minimum found at: "
		      << iter << " "
		      << gsl_vector_get(s->x, 0) << " "
		      << gsl_vector_get(s->x, 1) << " "
		      << gsl_vector_get(s->x, 2) << " value "
		      << s->f << "\n";

  } while (status == GSL_CONTINUE && iter < 1000);

  scored_triple_angle_set_t sas; // empty

  if (status == GSL_ENOPROG) {
     // std::cout << "No progress" << std::endl;
  } else {

     // sometimes we get wild solutions - not clear why (they seem to have sensibly log probs)
     // so trash those here.
     bool sane = true;
     for (std::size_t i=0; i<3; i++) {
	if (gsl_vector_get(s->x, i) < -M_PI_2) {
	   sane = false;
	   break;
	}
	if (gsl_vector_get(s->x, i) > (2.5 * M_PI)) {
	   sane = false;
	   break;
	}
     }

     if (sane) {
	std::vector<float> sv(3);
	for (std::size_t i=0; i<3; i++) {
	   sv[i] = gsl_vector_get(s->x, i);
	   if (sv[i] < 0) {
	      sv[i] += 2*M_PI;
	   }
	   if (sv[i] > 6.283) {
	      sv[i] -= 2*M_PI;
	   }
	}
	float score = s->f;
	sas = scored_triple_angle_set_t(tcs, sv, score);
     }
  }

  gsl_multimin_fdfminimizer_free(s);
  s = nullptr;
  gsl_vector_free(x);
  x = nullptr;
  return sas;

}

void
coot::crankshaft::run_optimizer_in_thread(const std::vector<std::size_t> &samples_for_thread,
					  const coot::nmer_crankshaft_set &cs,
					  const zo::rama_table_set &zorts,
					  std::vector<coot::crankshaft::scored_nmer_angle_set_t> *results) {
   if (cs.size() == 0) return;

   for (std::size_t j=0; j<samples_for_thread.size(); j++) {
      std::vector<float> start_angles(cs.n_peptides(), 0);
      for (std::size_t i=0; i<cs.n_peptides(); i++)
	 start_angles[i] = 2 * M_PI * float(coot::util::random())/float(RAND_MAX);

      coot::crankshaft::scored_nmer_angle_set_t sas = run_optimizer(start_angles, cs, zorts);
      run_optimizer(start_angles, cs, zorts);
      results->at(samples_for_thread[j]) = sas;

      // coot::crankshaft::scored_nmer_angle_set_t sas;
   }
}


coot::crankshaft::scored_nmer_angle_set_t
coot::crankshaft::run_optimizer(const std::vector<float> &start_angles,
				const coot::nmer_crankshaft_set &cs,
				const zo::rama_table_set &zorts) {

   scored_nmer_angle_set_t sas; // empty

   if (cs.size() == 0) { // Jude Short crash: crashes when cs.size() is 0 (it seems).
                         // This fixes the proximal cause, at least (gsl_vector_alloc(0)).
      std::cout << "ERROR:: empty crankshaft set " << std::endl;
      return sas;
   }

   size_t iter = 0;
   int status;

   nmer_set_param_holder_t param_holder(zorts, cs);

   gsl_vector *x;
   gsl_multimin_function_fdf my_func;

   my_func.n = cs.size();
   my_func.f = &optimize_an_nmer::f;
   my_func.df = &optimize_an_nmer::df;
   my_func.fdf = &optimize_an_nmer::fdf;
   my_func.params = reinterpret_cast<void *> (&param_holder);

   x = gsl_vector_alloc(cs.size());
   for (std::size_t i=0; i<start_angles.size(); i++)
      gsl_vector_set(x, i, start_angles[i]);

   const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_pr;
   gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, start_angles.size());

   // small moves Ellie, small moves...
   gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1.0); // within 1 degree is good enough

   do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);

      if (status)
	 break;

      status = gsl_multimin_test_gradient(s->gradient, 5e-2);

      if (status == GSL_SUCCESS)
	 if (false)
	    std::cout << "Minimum found at: "
		      << iter << " "
		      << gsl_vector_get(s->x, 0) << " "
		      << gsl_vector_get(s->x, 1) << " "
		      << gsl_vector_get(s->x, 2) << " value "
		      << s->f << "\n";

   } while (status == GSL_CONTINUE && iter < 1000);

   if (status == GSL_ENOPROG) {
      // std::cout << "No progress" << std::endl;
   } else {

      // sometimes we get wild solutions - not clear why (they seem to have sensibly log probs)
      // so trash those here.
      bool sane = true;
      for (std::size_t i=0; i<cs.size(); i++) {
	 if (gsl_vector_get(s->x, i) < -M_PI_2) {
	    sane = false;
	    break;
	 }
	 if (gsl_vector_get(s->x, i) > (2.5 * M_PI)) {
	    sane = false;
	    break;
	 }
      }

      if (sane) {
	 std::vector<float> sv(cs.size());
	 for (std::size_t i=0; i<cs.size(); i++) {
	    sv[i] = gsl_vector_get(s->x, i);
	    if (sv[i] < 0) {
	       sv[i] += 2*M_PI;
	    }
	    if (sv[i] > 6.283) {
	       sv[i] -= 2*M_PI;
	    }
	 }
	 float score = s->f;
	 sas = scored_nmer_angle_set_t(cs, sv, score);
      }
   }

   gsl_multimin_fdfminimizer_free(s);
   gsl_vector_free(x);
   return sas;

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

   if (v.size() == 0) {
      throw(std::runtime_error("unset crankshaft_set in phi_psi()"));
   }
   clipper::Coord_orth p0 = co(v[0]);
   clipper::Coord_orth p1 = co(v[1]);
   clipper::Coord_orth p_ca_1 = co(ca_1);
   clipper::Coord_orth p_ca_2 = co(ca_2);
   double torsion_phi = clipper::Coord_orth::torsion(p0, p1, p_ca_1, C_pos);
   double torsion_psi = clipper::Coord_orth::torsion(p1, p_ca_1, C_pos, N_pos);

   return phi_psi_t(torsion_phi, torsion_psi);
}

coot::phi_psi_t
coot::crankshaft_set::phi_psi(float ang) const {

   if (false) {
      std::cout << "here in phi_psi with ang " << ang << std::endl;
      std::cout << "here in phi_psi with ca_1 " << ca_1 << " " << atom_spec_t(ca_1) << std::endl;
      std::cout << "here in phi_psi with ca_2 " << ca_2 << " " << atom_spec_t(ca_2) << std::endl;
   }

   if (! ca_1) {
      std::string m = "ERROR:: ca_1 is unset in crankshaft_set::phi_psi()";
      if (true) {
	 std::cout << m << std::endl;
	 std::cout << "Here are the " << v.size() << " v atoms of the crankshaft_set " << std::endl;
	 for (std::size_t i=0; i<v.size(); i++)
	    std::cout << "   " << atom_spec_t(v[i]) << std::endl;
      }
      throw(std::runtime_error(m));
   }
   if (! ca_2) {
      std::string m = "ERROR:: ca_2 is unset in crankshaft_set::phi_psi()";
      if (true) {
	 std::cout << m << std::endl;
	 std::cout << "Here are the " << v.size() << " v atoms of the crankshaft_set " << std::endl;
	 for (std::size_t i=0; i<v.size(); i++)
	    std::cout << "   " << atom_spec_t(v[i]) << std::endl;
      }
      throw(std::runtime_error(m));
   }

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

   cs[0] = crankshaft_set(res_0, res_1, res_2, res_3);
   cs[1] = crankshaft_set(res_1, res_2, res_3, res_4);
   cs[2] = crankshaft_set(res_2, res_3, res_4, res_5);
   residue_types = residue_types_in;
}

// restores the atom positions in mol after write
void
coot::crankshaft::move_the_atoms_write_and_restore(scored_triple_angle_set_t sas,
						   const std::string &pdb_file_name) {

   std::map<mmdb::Atom *, clipper::Coord_orth> original_positions;
   std::map<mmdb::Atom *, clipper::Coord_orth>::const_iterator it;

   // 
   int indices[] = { 2, 3, 4, 5};
   for (std::size_t i=0; i<3; i++) {
      // maybe this should be done in a crankshaft_set
      for (std::size_t iat=0; iat<4; iat++) {
	 mmdb::Atom *at = sas[i].v[indices[iat]];
	 if (at) {
	    clipper::Coord_orth pos = co(at);
	    original_positions[at] = pos;
	 }
      }
      sas[i].move_the_atoms(sas.angles[i]);
      // std::cout << "writing " << pdb_file_name << std::endl;
      mol->WritePDBASCII(pdb_file_name.c_str());
   }

   // move the atoms back
   // std::cout << "repositioning " << original_positions.size() << " atoms " << std::endl;
   for (it=original_positions.begin(); it!=original_positions.end(); it++)
      update_position(it->first, it->second);

}

// move the atoms, create a copy of mol, restore the atom positions
mmdb::Manager *
coot::crankshaft::new_mol_with_moved_atoms(scored_triple_angle_set_t sas) {

   // needs thread protection
   // auto tp_1 = std::chrono::high_resolution_clock::now();
   std::map<mmdb::Atom *, clipper::Coord_orth> original_positions;
   std::map<mmdb::Atom *, clipper::Coord_orth>::const_iterator it;

   // 
   int indices[] = { 2, 3, 4, 5};
   for (std::size_t i=0; i<3; i++) {
      // Maybe this should be done in a crankshaft_set.
      // But no. Repositioning the atoms after the molecule copy
      // (which is ugly now) would be super-hideous.
      for (std::size_t iat=0; iat<4; iat++) {
	 mmdb::Atom *at = sas[i].v[indices[iat]];
	 if (at) {
	    clipper::Coord_orth pos = co(at);
	    original_positions[at] = pos;
	 }
      }
      sas[i].move_the_atoms(sas.angles[i]);
   }

   // auto tp_2 = std::chrono::high_resolution_clock::now();
   // ~ 1ms.
   mmdb::Manager *mol_new = new mmdb::Manager;
   mol_new->Copy(mol, mmdb::MMDBFCM_All);
   // auto tp_3 = std::chrono::high_resolution_clock::now();

   // move the atoms back
   // std::cout << "repositioning " << original_positions.size() << " atoms " << std::endl;
   for (it=original_positions.begin(); it!=original_positions.end(); it++)
      update_position(it->first, it->second);

   // auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
   // auto d32 = std::chrono::duration_cast<std::chrono::microseconds>(tp_3 - tp_2).count();
   // std::cout << "d21  " << d21 << " d32 " << d32 << " microseconds\n";

   return mol_new;
}

// move the atoms, create a copy of mol, restore the atom positions
mmdb::Manager *
coot::crankshaft::new_mol_with_moved_atoms(scored_nmer_angle_set_t sas) {

   std::map<mmdb::Atom *, clipper::Coord_orth> original_positions;
   std::map<mmdb::Atom *, clipper::Coord_orth>::const_iterator it;
   int indices[] = { 2, 3, 4, 5}; // v atom indices of the moving atoms in a crankshaft_set
   for (std::size_t i=0; i<sas.size(); i++) {
      for (std::size_t iat=0; iat<4; iat++) {
	 mmdb::Atom *at = sas[i].v[indices[iat]];
	 if (at) {
	    clipper::Coord_orth pos = co(at);
	    original_positions[at] = pos;
	 }
      }
      sas[i].move_the_atoms(sas.angles[i]);
   }

   mmdb::Manager *mol_new = new mmdb::Manager;
   mol_new->Copy(mol, mmdb::MMDBFCM_All);

   for (it=original_positions.begin(); it!=original_positions.end(); it++)
      update_position(it->first, it->second);

   return mol_new;
}



// -----------------------------------------------------------------------------
//   refine (in threads)
// -----------------------------------------------------------------------------

// static
coot::crankshaft::molecule_score_t
coot::crankshaft::refine_and_score_mol(mmdb::Manager *mol,
				       const std::vector<coot::residue_spec_t > &refine_residue_specs,
				       const std::vector<coot::residue_spec_t> &residue_specs_for_scoring,
				       const coot::protein_geometry &geom,
				       const clipper::Xmap<float> &xmap,
				       float map_weight,
				       const std::string &output_pdb_file_name,
				       ctpl::thread_pool *thread_pool_p, int n_threads_in) {

   molecule_score_t score;

   if (mol) {
      float radius = 2.2; // for the moment
      std::vector<std::pair<bool, mmdb::Residue *> > residues;
      std::vector<mmdb::Link> links;
      std::vector<coot::atom_spec_t> fixed_atom_specs;
      coot::restraint_usage_Flags flags = BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;

      flags = TYPICAL_RESTRAINTS;
      int imol = 0; // dummy
      int nsteps_max = 4000;
      bool make_trans_peptide_restraints = true;
      short int print_chi_sq_flag = 1;
      bool do_rama_plot_restraints = true;
      pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
      int restraints_rama_type = restraints_container_t::RAMA_TYPE_ZO;

      std::vector<std::pair<bool, mmdb::Residue *> > refine_residues;
      for (std::size_t ires=0; ires<refine_residue_specs.size(); ires++) {
	 mmdb::Residue *residue_p = util::get_residue(refine_residue_specs[ires], mol);
	 if (residue_p) {
	    std::pair<bool, mmdb::Residue *> p(false, residue_p);
	    refine_residues.push_back(p);
	 }
      }


      auto tp_0 = std::chrono::high_resolution_clock::now();

      restraints_container_t restraints(refine_residues, links, geom, mol, fixed_atom_specs, &xmap);
      restraints.thread_pool(thread_pool_p, n_threads_in);
      restraints.set_quiet_reporting();
      restraints.add_map(map_weight);
      restraints.set_rama_type(restraints_rama_type);
      restraints.set_rama_plot_weight(1);
      restraints.make_restraints(imol, geom, flags, 1, make_trans_peptide_restraints,
				 1.0, do_rama_plot_restraints, true, true, false, pseudos);
      restraints.minimize(flags, nsteps_max, print_chi_sq_flag);
      if (! output_pdb_file_name.empty())
	 restraints.write_new_atoms(output_pdb_file_name);

      coot::geometry_distortion_info_container_t gdic = restraints.geometric_distortions();
      for (std::size_t id=0; id<gdic.geometry_distortion.size(); id++) {
	 // std::cout << "   " << gdic.geometry_distortion[id] << std::endl;
      }


      auto tp_1 = std::chrono::high_resolution_clock::now();

      bool mainchain_only_flag = true;
      float score_map = coot::util::map_score_by_residue_specs(mol, residue_specs_for_scoring,
							       xmap, mainchain_only_flag);

      auto tp_2 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
      if (false) {
	 std::cout << "refine mol : d10  " << d10 << " microseconds\n";
	 std::cout << "map_score_by_residue_specs: d21  " << d21 << " microseconds\n";
	 std::cout << "scores map: " << score_map << " distortion " << gdic.distortion()
		   << " " << output_pdb_file_name << std::endl;
      }

      score.density_score = score_map;
      score.model_score = gdic.distortion();
   }
   return score;
}



// for use the threads:
// static
void
coot::crankshaft::refine_and_score_mols(std::vector<mmdb::Manager *> mols,
					const std::vector<unsigned int> &mols_thread_vec,
					const std::vector<coot::residue_spec_t> &residue_specs_for_refining,
					const std::vector<coot::residue_spec_t> &residue_specs_for_scoring,
					const coot::protein_geometry &geom,
					const clipper::Xmap<float> &xmap,
					float map_weight,
					std::vector<molecule_score_t> *mol_scores,
					ctpl::thread_pool *thread_pool_p, int n_threads) {

   for (std::size_t i=0; i<mols_thread_vec.size(); i++) {

      molecule_score_t ms = refine_and_score_mol(mols[mols_thread_vec[i]],
						 residue_specs_for_refining,
						 residue_specs_for_scoring,
						 geom, xmap, map_weight, "", thread_pool_p, n_threads);
      mol_scores->at(mols_thread_vec[i]) = ms;
   }
}

// return at max n_solutions
//
// static
std::vector<mmdb::Manager *>
coot::crankshaft::crank_refine_and_score(const coot::residue_spec_t &rs, // mid-residue
					 unsigned int n_peptides,
					 const clipper::Xmap<float> &xmap,
					 // atom_selection_container_t asc,
					 mmdb::Manager *mol_in,
					 float map_weight,
					 int n_samples, int n_solutions,
					 ctpl::thread_pool *thread_pool_p, int n_threads_in) {

   std::vector<mmdb::Manager *> solution_molecules;

   // promote this to an argument for this function
   float log_prob_filter_n_sigma = 1.0; // only the top few.

   mmdb::Residue *prev_res = util::get_previous_residue(rs, mol_in);
   if (! prev_res) {
      std::cout << "WARNING:: No residue previous to " << rs << std::endl;
   } else {
      crankshaft cs(mol_in);
      zo::rama_table_set zorts;

      residue_spec_t prev_residue_spec(prev_res);
      logger.log(log_t::INFO, logging::ltw("using residue specifier: "), logging::ltw(rs.format()));
      // std::cout << "INFO:: using residue specifier: " << rs << std::endl;

      if (n_samples == -1) {
	 // choose for me (based on the number of threads and the number of peptides

	 unsigned int n_threads = coot::get_max_number_of_threads();

	 // the more peptides, the less samples
	 if (n_threads > 0)
	    n_samples = std::lround(100 * n_threads/pow(n_peptides, 0.7));
	 else
	    n_samples = 5;
	 logger.log(log_t::INFO, logging::ltw("Chose n_samples "), logging::ltw(n_samples), logging::ltw(" for "), logging::ltw(n_peptides), logging::ltw(" peptides  and using "), logging::ltw(static_cast<unsigned int>(n_threads)));
	 // std::cout << "INFO:: Chose n_samples " << n_samples << " for " << n_peptides
	 //           << " peptides " << " and using " << n_threads << " threads" << std::endl;
      }

      // changed the API to find_maxima to be the middle of the n_peptides (3)
      // residues
      //

      std::vector<coot::crankshaft::scored_nmer_angle_set_t> sas =
	 cs.find_maxima(rs, n_peptides, zorts, log_prob_filter_n_sigma, n_samples);

      logger.log(log_t::INFO, logging::ltw("Will refine "), logging::ltw(static_cast<unsigned int>(sas.size())), logging::ltw(" crankshaft solutions"));
      // std::cout << "INFO:: Will refine " << sas.size() << " crankshaft solutions" << std::endl;
      if (false)
	 for (std::size_t i=0; i<sas.size(); i++)
	 std::cout << "   " << sas[i] << std::endl;

      coot::protein_geometry geom;
      geom.set_verbose(0);
      geom.init_standard();

      if (false) { // debug models
	 for (std::size_t i=0; i<sas.size(); i++) {
	    std::string output_pdb_file_name = "crankshaft-pre-ref-" +
	       coot::util::int_to_string(i) + ".pdb";
	    mmdb::Manager *mol = cs.new_mol_with_moved_atoms(sas[i]); // delete mol
	    mol->WritePDBASCII(output_pdb_file_name.c_str());
	    delete mol;
	 }
      }

      // which residues shall we refine?
      std::vector<coot::residue_spec_t> local_residue_specs;
      if (sas.size() > 0) {
	 std::vector<mmdb::Residue *> cs_residues = sas[0].residues();
	 for (std::size_t ires=0; ires<cs_residues.size(); ires++)
	    local_residue_specs.push_back(coot::residue_spec_t(cs_residues[ires]));
	 // maybe add some near neighbours?
      }

      logger.log(log_t::INFO, "refine these residue specs: ");
      // std::cout << "INFO:: refine these residue specs: " << std::endl;
      for (std::size_t ilr=0; ilr<local_residue_specs.size(); ilr++)
	 std::cout << "   " << local_residue_specs[ilr];
      std::cout << std::endl;

      std::vector<mmdb::Manager *> mols(sas.size(), 0);
      std::vector<crankshaft::molecule_score_t> mol_scores(sas.size());
      for (std::size_t i=0; i<sas.size(); i++) {
	 // delete mol when finished with it
	 mmdb::Manager *mol = cs.new_mol_with_moved_atoms(sas[i]); // d
	 mols[i] = mol;
      }


      unsigned int n_threads = coot::get_max_number_of_threads();

      if (n_threads > 0) {

	 std::vector<std::thread> threads;
	 unsigned int thread_idx = 0;
	 // std::cout << "debug:: making mols_thread_vec " << n_threads << std::endl;
	 std::vector<std::vector<unsigned int> > mols_thread_vec(n_threads);
	 // std::cout << "debug:: done making mols_thread_vec " << std::endl;
	 // put the mols in mols_thread_vec
	 for (unsigned int ii=0; ii<mols.size(); ii++) {
	    // std::cout << "adding mol " << ii << " to thread " << thread_idx << std::endl;
	    mols_thread_vec[thread_idx].push_back(ii);
	    thread_idx++;
	    if (thread_idx == n_threads) thread_idx = 0;
	 }

	 if (false) {
	    std::cout << "debug:: the mols_thread_vecs are: " << std::endl;
	    for (std::size_t j=0; j<mols_thread_vec.size(); j++) {
	       std::cout << "thread " << j << " of " << mols_thread_vec.size() << " threads " << std::endl;
	       for (std::size_t k=0; k<mols_thread_vec[j].size(); k++) {
		  std::cout << " " << mols_thread_vec[j][k];
	       }
	       std::cout << std::endl;
	    }
	 }

	 for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {

	    auto tp_0 = std::chrono::high_resolution_clock::now();
	    const std::vector<unsigned int> &i_mols = mols_thread_vec[i_thread];

	    // Don't bother to push a thread if there are no molecules to refine
	    // (this happens when there are more threads than there
	    // are crankshaft probability maxima/solutions).
	    //
	    if (mols_thread_vec[i_thread].size() > 0) {
	       threads.push_back(std::thread(refine_and_score_mols, mols, mols_thread_vec[i_thread],
					     std::cref(local_residue_specs), std::cref(local_residue_specs),
					     std::cref(geom), std::cref(xmap), map_weight,
					     &mol_scores, thread_pool_p, n_threads));
	       // threads.push_back(std::thread(dummy_func, mols)); // testing
	       auto tp_1 = std::chrono::high_resolution_clock::now();
	       auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
	       if (false)
		  std::cout << "time to push thread: " << d10 << " microseconds\n";
	    }
	 }

	 auto tp_2 = std::chrono::high_resolution_clock::now();
	 for (unsigned int i_thread=0; i_thread<threads.size(); i_thread++)
	    threads.at(i_thread).join();
	 auto tp_3 = std::chrono::high_resolution_clock::now();
	 auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
	 if (false)
	    std::cout << "debug:: time to join threads: " << d32 << " milliseconds\n";

      } else {
	 std::cout << "ERROR:: No threads" << std::endl;
      }

      for (std::size_t i=0; i<sas.size(); i++) {

	 // the higher the combi_score the better
	 const molecule_score_t &ms = mol_scores[i];
	 float combi_score = 0.01 * map_weight * ms.density_score - ms.model_score - sas[i].minus_log_prob;
	 sas[i].set_combi_score(ms.density_score, map_weight, ms.model_score);
	 if (false) {
	    std::cout << "scores: " << i << " minus-log-prob "
		      << std::setw(9) << sas[i].minus_log_prob << " combi-score "
		      << std::setw(9) << sas[i].combi_score << " "
		      << std::setw(9) << ms.density_score << " "
		      << std::setw(9) << ms.model_score << " combi-score "
		      << std::setw(9) << std::right << std::setprecision(3) << std::fixed
		      << sas[i].combi_score << std::endl;
	    std::cout.setf(std::ios::fixed, std::ios::floatfield); // reset from std::fixed
	 }
      }

      // Now sort the molecules using the sas combi-scores. To do that we
      // need the sas and the molecule to be attached
      std::vector<std::pair<scored_nmer_angle_set_t, mmdb::Manager *> > sas_mol_pairs(sas.size());
      for (std::size_t i=0; i<sas.size(); i++) {
	 sas_mol_pairs[i] = std::pair<scored_nmer_angle_set_t, mmdb::Manager *> (sas[i], mols[i]);
      }

      std::sort(sas_mol_pairs.begin(), sas_mol_pairs.end(),
		scored_nmer_angle_set_t::sorter_by_combi_score);

      for (std::size_t i=0; i<sas_mol_pairs.size(); i++) {
	 std::cout << "scores: " << i << " minus-log-prob "
		   << std::setw(9) << sas_mol_pairs[i].first.minus_log_prob << " combi-score "
		   << std::setw(9) << std::right << std::setprecision(3) << std::fixed
		   << sas_mol_pairs[i].first.combi_score << std::endl;
	 std::cout.setf(std::ios::fixed, std::ios::floatfield); // reset from std::fixed
      }

      for (std::size_t isol=0; isol<sas_mol_pairs.size(); isol++)
	 if (static_cast<int>(isol) < n_solutions)
	    solution_molecules.push_back(sas_mol_pairs[isol].second);
      
      // write all solutions - for debugging
      //
      if (false) {
	 for (std::size_t i=0; i<sas.size(); i++) {
	    std::string file_name = "crankshaft-post-refinement-";
	    file_name += coot::util::int_to_string(i);
	    file_name += ".pdb";
	    mols[i]->WritePDBASCII(file_name.c_str());
	 }
      }

      // delete the molecules that we are not returning
      for (std::size_t i=0; i<sas.size(); i++)
	 if (std::find(solution_molecules.begin(),
		       solution_molecules.end(), mols[i]) == solution_molecules.end())
	    delete mols[i];
   }

   return solution_molecules;
}
