
#include "crankshaft.hh"

#include "geometry/main-chain.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "zo-rama.hh"

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

std::vector<std::pair<float, float> >
coot::crankshaft::spin_search(const coot::residue_spec_t &spec) const {

   std::vector<std::pair<float, float> > v;

   std::pair<mmdb::Residue *, mmdb::Residue *> rs = util::get_this_and_next_residues(spec, mol);
   if (rs.first) {
      if (rs.second) {
	 mmdb::Residue *res_1 = rs.first;
	 mmdb::Residue *res_2 = rs.second;
	 std::vector<mmdb::Atom *> mc_atoms = get_mainchain_atoms(res_1, res_2);
	 if (mc_atoms.size()) {
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
		     int n_samples = 60;
		     float div = 1/float(n_samples);
		     crankshaft_set cs(res_0, res_1, res_2, res_3);
		     for (int i=0; i<n_samples; i++) {
			float a = float(i) * div * 2*M_PI;
			std::pair<float, float> pr = cs.probability_of_spin_orientation(a);
			std::cout << "score: " << i << " " << pr.first << " " << pr.second
				  << std::endl;
		     }
		  }
	       }
	    } else {
	       std::cout << "missing mainchain atom(s) for " << spec << std::endl;
	    }
	 }
      } else {
	 std::cout << "missing second residue " << spec << std::endl;
      }
   } else {
      std::cout << "missing first residue " << spec << std::endl;
   }
   return v;
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



// ang in radians
std::pair<float, float>
coot::crankshaft_set::probability_of_spin_orientation(float ang) const {

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

   std::pair<phi_psi_t, phi_psi_t> ppp = get_phi_phis(C_new, N_new);

   std::string res_type_1 = "ALL!nP";
   std::string res_type_2 = "ALL!nP";

   zo::rama_table_set rts;

   float v1 = rts.value(res_type_1, ppp.first.phi,  ppp.first.psi);
   float v2 = rts.value(res_type_2, ppp.second.phi, ppp.second.psi);

   return std::pair<float, float> (v1, v2);
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
	    std::vector<std::pair<float, float> > r = spin_search(rs);
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
