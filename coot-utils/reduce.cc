
#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"
#include "reduce.hh"

   // Spin-search: OH, SH, NH3, MET methyls

   // OH: SER, TYR, THR

   // look at the Cartesian function position_by_torsion().
   
   // main-chain
   //
   // placeable_main[" H  "] // bisect C(n-1)-N(n) and CA(n)-N(n)

   // CA is placed by tetrahedron - c.f. push_chiral_hydrogen, but use unit vectors to
   //                                    the neighbours

   // Hs on CB: HB2 and HB3:
   // for CYS, ASP, GLU, PHE, HIS, LYS, LEU, MET, ASN, PRO, GLU, ARG, 
   //     SER, TYR, VAL, TRP, TYR
   //    bisect CA-CB and CG-CB, call that b
   //           CA-CG, unit, call that c,
   //           a fragment of bxc and of b gives the delta from CB for HB2 and HB3

   // ALA: CB Hs: HB1, HB2, HB3, BL-A-T: average positions from N-CA-CB-HBx C-CA-CB-HBx

   // CYS: SG H : HB BL-A-T: CA-CB-SG-HG (180)

   // GLU: CG Hs: Bisect as CB Hs, but use CB-CG-CD

   // GLN: CG Hs: Bisect as CB Hs, but use CB-CG-CD

   // PHE Ring Hs: average BL-A-T of ring torsions (180)

   // HIS: average of  BL-A-T of ring torsions (180)

   // ILE: CB HB: on CB from tetrahedron of CA, CG1, CG2.
   //           HG1 and H2: bisect CB, CG, CD1
   //           HD1, HD2, HD3 BL-A-T from CB, CG1, CD1
   //           HG21 HG22, HG33: BL-A-T from CA-CB-CG2-HG2x

   // LYS: CG Hs: Bisect as CB Hs, but use CB-CG-CD
   //      CD Hs: Bisect as CB Hs, but use CG-CD-CE
   //      CE Hs: Bisect as CB Hs, but use CD-CE-NZ
   //      NZ Hs: HZ1,2,3: BL-A-T from  CD-CE-NZ-HZx

   // LEU:

   // MET: CG Hs: Bisect as CB Hs, but use CB-CG-SD
   //         HE1,2,3 BL-A-T from CG, SD, CG

   // PRO: +++

   // GLU: CG Hs: Bisect as CB Hs, but use CB-CG-CD

   // ARG: CG Hs: Bisect as CB Hs, but use CB-CG-CD
   //      CD Hs: Bisect as CB Hs, but use CG-CD-CE
   //      CE Hs: Bisect as CB Hs, but use CD-CE-NZ
   //      +++

   // SER: BL-A-T: CA, CG, OG.

   // THR: CG2 Hs: HG2[1,2,3] BL-A-T: CA, CG, CG.
   //         HG1: BL-A-T on CA, CB, OG (180)  # spin-search

   // VAL: HGxy:  BL-A-T on CA, CB, CG[1,2], HGxy

   // TRP: HD1: average BL-A-T: CD2, CG, CD1 and CE2, NE1, CD1
   //      HE1: average BL-A-T: CG, CD1, NE1 and CD2, CE2, NE1
   //      HZ2: average BL-A-T: CD2, CE2, CZ2 and CZ3, CH2, CZ2
   //      HH2: average BL-A-T: CE2, CZ2, CH2 and CE3, CZ3, CH2
   //      HZ3: average BL-A-T: CZ2, CH2, CZ3 and CD2, CE3, CZ3
   //      HE3: average BL-A-T: CE2, CD2, CE3 and CH2, CZ3, CE3

   // TYR Ring Hs: average BL-A-T of ring torsions (180)
   //          HH: BL-A-T: CE1, CZ, OH # spin-search


clipper::Coord_orth
coot::reduce::position_by_bond_length_angle_torsion(mmdb::Atom *at_1,  // CA
						    mmdb::Atom *at_2,  // CB
						    mmdb::Atom *at_3,  // CG
						    double bl,
						    double angle_rad,
						    double torsion_rad) const {
   
   clipper::Coord_orth at_1_pos = co(at_1);
   clipper::Coord_orth at_2_pos = co(at_2);
   clipper::Coord_orth at_3_pos = co(at_3);

   clipper::Coord_orth new_pos(at_1_pos, at_2_pos, at_3_pos, bl, angle_rad, torsion_rad);

   return new_pos;
}

clipper::Coord_orth
coot::reduce::position_by_bisection(mmdb::Atom *at_1,
				    mmdb::Atom *at_2,
				    mmdb::Atom *at_3,
				    double bl) const {

   // bisect the normals
   clipper::Coord_orth at_1_pos = co(at_1);
   clipper::Coord_orth at_2_pos = co(at_2);
   clipper::Coord_orth at_3_pos = co(at_3);

   clipper::Coord_orth vec_1_to_2(at_2_pos-at_1_pos);
   clipper::Coord_orth vec_3_to_2(at_2_pos-at_3_pos);
   clipper::Coord_orth vec_1_to_2_uv(vec_1_to_2.unit());
   clipper::Coord_orth vec_3_to_2_uv(vec_3_to_2.unit());

   // vector from point to mid-atom
   clipper::Coord_orth bisect_delta = 0.5 * (vec_1_to_2_uv + vec_3_to_2_uv);
   clipper::Coord_orth bisect_delta_uv(bisect_delta.unit());

   clipper::Coord_orth Hp1 = at_2_pos + bl * bisect_delta_uv;

   return Hp1;
}

std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::reduce::position_pair_by_bisection(mmdb::Atom *at_1,  // CA
					 mmdb::Atom *at_2,  // CB
					 mmdb::Atom *at_3,  // CG
					 double bl,
					 double alpha // angle btwn the two H atoms
					 ) const {

   // bisect the normals
   clipper::Coord_orth at_1_pos = co(at_1);
   clipper::Coord_orth at_2_pos = co(at_2);
   clipper::Coord_orth at_3_pos = co(at_3);

   clipper::Coord_orth vec_1_to_2(at_2_pos-at_1_pos);
   clipper::Coord_orth vec_3_to_2(at_2_pos-at_3_pos);
   clipper::Coord_orth vec_1_to_2_uv(vec_1_to_2.unit());
   clipper::Coord_orth vec_3_to_2_uv(vec_3_to_2.unit());

   // vector from point to mid-atom
   clipper::Coord_orth bisect_delta = 0.5 * (vec_1_to_2_uv + vec_3_to_2_uv);
   clipper::Coord_orth bisect_delta_uv(bisect_delta.unit());

   // vector from CA->CG
   clipper::Coord_orth vec_1_to_3 = at_3_pos - at_1_pos;
   clipper::Coord_orth vec_1_to_3_uv(vec_1_to_3.unit());

   // cpu, the vector between the H atoms is along this vector
   clipper::Coord_orth cpu(clipper::Coord_orth::cross(vec_1_to_3_uv, bisect_delta_uv));

   double scale_fac_bisector = bl * sin(0.5*(M_PI-alpha));
   double scale_fac_cpu      = bl * cos(0.5*(M_PI-alpha));
   
   clipper::Coord_orth Hp1 = at_2_pos + scale_fac_bisector * bisect_delta_uv - scale_fac_cpu * cpu;
   clipper::Coord_orth Hp2 = at_2_pos + scale_fac_bisector * bisect_delta_uv + scale_fac_cpu * cpu;

   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (Hp1, Hp2);
}

clipper::Coord_orth
coot::reduce::position_by_tetrahedron(mmdb::Atom *at_central,
				      mmdb::Atom *at_n_1,
				      mmdb::Atom *at_n_2,
				      mmdb::Atom *at_n_3,
				      double bl) const {

   clipper::Coord_orth at_c_pos = co(at_central);
   clipper::Coord_orth at_n1_pos = co(at_n_1);
   clipper::Coord_orth at_n2_pos = co(at_n_2);
   clipper::Coord_orth at_n3_pos = co(at_n_3);
   clipper::Coord_orth vec_1_to_c(at_c_pos-at_n1_pos);
   clipper::Coord_orth vec_2_to_c(at_c_pos-at_n2_pos);
   clipper::Coord_orth vec_3_to_c(at_c_pos-at_n3_pos);
   clipper::Coord_orth vec_1_to_c_uv(vec_1_to_c.unit());
   clipper::Coord_orth vec_2_to_c_uv(vec_2_to_c.unit());
   clipper::Coord_orth vec_3_to_c_uv(vec_3_to_c.unit());

   clipper::Coord_orth under_pos = 0.3333333333 *
      (vec_1_to_c_uv + vec_2_to_c_uv + vec_3_to_c_uv);
   clipper::Coord_orth under_pos_uv(under_pos.unit());

   clipper::Coord_orth H_pos = at_c_pos + bl * under_pos_uv;
   return H_pos;
}

void
coot::reduce::add_hydrogen_atoms() {

   if (mol) {
      add_riding_hydrogens();
      mol->FinishStructEdit();
   }
}

void
coot::reduce::add_riding_hydrogens() {

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 mmdb::Residue *residue_prev_p;
	 mmdb::Atom *at;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p      = chain_p->GetResidue(ires);
	    if (ires > 0)
	       residue_prev_p = chain_p->GetResidue(ires-1);
	    else
	       residue_prev_p = 0;
	    // what about strange missing residues - where we can place the CA HA
	    // (but not the N's H).
	    add_riding_hydrogens(residue_p, residue_prev_p);
	 }
      }
   }
}

// only call this with a valid residue_p
void
coot::reduce::add_riding_hydrogens(mmdb::Residue *residue_p, mmdb::Residue *residue_prev_p) {

   std::string res_name = residue_p->GetResName();
   double bl = 0.97;
   double bl_arom  = 0.93;
   double bl_amino = 0.86;
   double bl_oh    = 0.84;
   double bl_sh    = 1.2;
   if (res_name == "ALA") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      torsion_info_t torsion_1(" N  ", " CA ", " CB ", bl, 109, 180);
      add_methyl_Hs(" HB1", " HB2", " HB3", torsion_1, residue_p);
   }
   if (res_name == "CYS") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " SG ", bl, 107, residue_p);
      add_SH_H(" HG ", " SG ", " CB ", " CA ", bl_sh, 109.5, 180, residue_p);
   }
   if (res_name == "ASP") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
   }
   if (res_name == "GLU") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
   }
   if (res_name == "PHE") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_aromatic_hydrogen(" HD1", " CG ", " CD1", " CE1", bl_arom, residue_p);
      add_aromatic_hydrogen(" HE1", " CD1", " CE1", " CZ ", bl_arom, residue_p);
      add_aromatic_hydrogen(" HD2", " CG ", " CD2", " CE2", bl_arom, residue_p);
      add_aromatic_hydrogen(" HE2", " CD2", " CE2", " CZ ", bl_arom, residue_p);
      add_aromatic_hydrogen(" HZ ", " CE1", " CZ ", " CE2", bl_arom, residue_p);
   }
   if (res_name == "GLY") {
      add_main_chain_hydrogens(residue_p, residue_prev_p, true);
   }
   if (res_name == "HIS") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_his_ring_C_Hs(residue_p);
   }
   if (res_name == "ILE") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens("HG12", "HG13", " CB ", " CG1", " CD1", bl, 107, residue_p);
      torsion_info_t ti(" CB ", " CG1", " CD1", bl, 109, 180);
      add_methyl_Hs("HD11", "HD12", "HD13", ti, residue_p);
      torsion_info_t t2(" CA ", " CB ", " CG2", bl, 109, 180);
      add_methyl_Hs("HG21", "HG22", "HG23", t2, residue_p);
      add_tetrahedral_hydrogen(" HB ", " CB ", " CA ", " CG1", " CG2", bl, residue_p);
   }
   if (res_name == "LYS") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HD2", " HD3", " CG ", " CD ", " CE ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HE2", " HE3", " CD ", " CE ", " NZ ", bl, 107, residue_p);
      torsion_info_t ti(" CD ", " CE ", " NZ ", bl, 109, 180);
      add_methyl_Hs(" HZ1", " HZ2", " HZ3", ti, residue_p);
   }
   if (res_name == "LEU") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      torsion_info_t t1(" CB ", " CG ", " CD1", bl, 109, 180);
      torsion_info_t t2(" CB ", " CG ", " CD2", bl, 109, 180);
      add_methyl_Hs("HD11", "HD12", "HD13", t1, residue_p);
      add_methyl_Hs("HD21", "HD22", "HD23", t2, residue_p);
      add_tetrahedral_hydrogen(" HG ", " CG ", " CB ", " CD1", " CD2", bl, residue_p);
   }
   if (res_name == "MET") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " SD ", bl, 107, residue_p);
      torsion_info_t t1(" CG ", " SD ", " CE", bl, 109, 180);
      add_methyl_Hs(" HE1", " HE2", " HE3", t1, residue_p);
   }
   if (res_name == "MSE") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB1", " HB2", " CA ", " CB ", " CG ", bl, 107, residue_p);
   }
   if (res_name == "ASN") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_amino_hydrogens("HD21", "HD22", " ND2", " CG ", " OD1", bl_amino, residue_p);
   }
   if (res_name == "PRO") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HD2", " HD3", " CG ", " CD ", " N  ", bl, 107, residue_p);
   }
   if (res_name == "GLN") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      add_amino_hydrogens("HE21", "HE22", " NE2", " CD ", " OE1", bl_amino, residue_p);
   }
   if (res_name == "ARG") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HD2", " HD3", " CG ", " CD ", " NE ", bl, 107, residue_p);
      add_guanidinium_hydrogens(residue_p);
   }
   if (res_name == "SER") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " OG ", bl, 107, residue_p);
      add_OH_H(" HG ", " OG ", " CB ", " CA ", bl_oh, 109.5, 180, residue_p);
   }
   if (res_name == "THR") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_tetrahedral_hydrogen(" HB ", " CB ", " CA ", " OG1", " CG2", bl, residue_p);
      torsion_info_t ti(" CA ", " CB ", " CG2", bl, 109, 180);
      add_methyl_Hs("HG21", "HG22", "HG23", ti, residue_p);
      add_OH_H(" HG1", " OG1", " CB ", " CA ", bl_oh, 109.5, 180, residue_p);
   }
   if (res_name == "VAL") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      torsion_info_t t1(" CA ", " CB ", " CG1 ", bl, 109, 180);
      torsion_info_t t2(" CA ", " CB ", " CG2 ", bl, 109, 180);
      add_methyl_Hs("HG11", "HG12", "HG13", t1, residue_p);
      add_methyl_Hs("HG21", "HG22", "HG23", t2, residue_p);
      add_tetrahedral_hydrogen(" HB ", " CB ", " CA ", " CG1", " CG2", bl, residue_p);
   }
   if (res_name == "TRP") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_trp_indole_hydrogens(residue_p);
   }
   if (res_name == "TYR") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_aromatic_hydrogen(" HD1", " CG ", " CD1", " CE1", bl_arom, residue_p);
      add_aromatic_hydrogen(" HE1", " CD1", " CE1", " CZ ", bl_arom, residue_p);
      add_aromatic_hydrogen(" HD2", " CG ", " CD2", " CE2", bl_arom, residue_p);
      add_aromatic_hydrogen(" HE2", " CD2", " CE2", " CZ ", bl_arom, residue_p);
      add_OH_H(" HH ", " OH ", " CZ ", " CE1", bl_oh, 109.5, 180, residue_p);
   }
}

// is_gly is false by default
void
coot::reduce::add_main_chain_hydrogens(mmdb::Residue *residue_p, mmdb::Residue *residue_prev_p,
				       bool is_gly) {

   if (is_gly) {
      double bl = 0.97;
      add_2_sp3_hydrogens(" HA2", " HA3", " N  ", " CA ", " C  ", bl, 107, residue_p);
      add_main_chain_H(residue_p, residue_prev_p);
   } else {
      add_main_chain_HA(residue_p);
      std::string residue_name(residue_p->GetResName());
      if (util::is_standard_amino_acid_name(residue_name))
	 if (residue_name != "PRO")
	    add_main_chain_H(residue_p, residue_prev_p);
   }
}


// The H on the N
// This function can be (is) called for the first residue in a chain, that is the N-terminus
// and doesn't have a previous residue. That is checked for.
void
coot::reduce::add_main_chain_H(mmdb::Residue *residue_p, mmdb::Residue *residue_prev_p) {

   if (residue_prev_p) {
      // Try position by torsion based on O-C-N-H
      double bl = 0.86;
      if (residue_p->isNTerminus()) {
	 // NH3+ - needs spin search - these are not riding
      } else {
	 std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
	 for (unsigned int i=0; i<alt_confs.size(); i++) {
	    mmdb::Atom *at_ca     = residue_p->GetAtom(" CA ", 0, alt_confs[i].c_str());
	    mmdb::Atom *at_n      = residue_p->GetAtom(" N  ", 0, alt_confs[i].c_str());
	    mmdb::Atom *at_c_prev = residue_prev_p->GetAtom(" C  ", 0, alt_confs[i].c_str());
	    mmdb::Atom *at_o_prev = residue_prev_p->GetAtom(" O  ", 0, alt_confs[i].c_str());
	    if (at_ca && at_n && at_c_prev && at_o_prev) {

	       clipper::Coord_orth at_c_pos  = co(at_c_prev);
	       clipper::Coord_orth at_o_pos  = co(at_o_prev);
	       clipper::Coord_orth at_n_pos  = co(at_n);
	       clipper::Coord_orth at_ca_pos = co(at_ca);
	       double angle = clipper::Util::d2rad(125.0);
	       clipper::Coord_orth H_pos(at_ca_pos, at_c_pos, at_n_pos, bl, angle, M_PI);
	       mmdb::realtype bf = at_n->tempFactor;
	       add_hydrogen_atom(" H  ", H_pos, bf, residue_p);
	    }
	 }
      }
   }
}

// The H on the CA
void
coot::reduce::add_main_chain_HA(mmdb::Residue *residue_p) {

   double bl = 0.97;

   // PDBv3 FIXME
   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      mmdb::Atom *at_ca = residue_p->GetAtom(" CA ", 0, alt_confs[i].c_str());
      mmdb::Atom *at_n1 = residue_p->GetAtom(" C  ", 0, alt_confs[i].c_str());
      mmdb::Atom *at_n2 = residue_p->GetAtom(" N  ", 0, alt_confs[i].c_str());
      mmdb::Atom *at_n3 = residue_p->GetAtom(" CB ", 0, alt_confs[i].c_str());
      if (at_ca && at_n1 && at_n2 && at_n3) {
	 clipper::Coord_orth pos = position_by_tetrahedron(at_ca, at_n1, at_n2, at_n3, bl);
	 mmdb::realtype bf = at_ca->tempFactor;
	 add_hydrogen_atom(" HA ", pos, bf, residue_p);
      }
   }
}

void
coot::reduce::add_hydrogen_atom(std::string atom_name, clipper::Coord_orth &pos,
				mmdb::realtype bf, mmdb::Residue *residue_p) {

   mmdb::Atom *new_H = new mmdb::Atom;
   new_H->SetAtomName(atom_name.c_str());
   new_H->SetElementName(" H"); // PDBv3 FIXME
   new_H->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0, bf);
   residue_p->AddAtom(new_H);
}


void 
coot::reduce::add_methyl_Hs(const std::string &at_name_1,  // HB1 (for example)
			    const std::string &at_name_2,  // HB2 + 120 degress
			    const std::string &at_name_3,  // HB3 - 120 degree
			    torsion_info_t torsion_1, torsion_info_t torsion_2,
			    mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      clipper::Coord_orth p11;
      clipper::Coord_orth p12;
      clipper::Coord_orth p13;
      clipper::Coord_orth p21;
      clipper::Coord_orth p22;
      clipper::Coord_orth p23;
      bool have_1 = false;
      bool have_2 = false;
      mmdb::Atom *at_1 = residue_p->GetAtom(torsion_1.at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_2 = residue_p->GetAtom(torsion_1.at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_3 = residue_p->GetAtom(torsion_1.at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_1 && at_2 && at_3) {
	 have_1 = true;
	 p11 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_1.bond_length,
						     clipper::Util::d2rad(torsion_1.angle_deg),
						     clipper::Util::d2rad(torsion_1.torsion_deg));
	 p12 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_1.bond_length,
						     clipper::Util::d2rad(torsion_1.angle_deg),
						     clipper::Util::d2rad(torsion_1.torsion_deg + 120));
	 p13 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_1.bond_length,
						     clipper::Util::d2rad(torsion_1.angle_deg),
						     clipper::Util::d2rad(torsion_1.torsion_deg - 120));
      }
      at_1 = residue_p->GetAtom(torsion_2.at_name_1.c_str(), 0, alt_confs[i].c_str());
      at_2 = residue_p->GetAtom(torsion_2.at_name_2.c_str(), 0, alt_confs[i].c_str());
      at_3 = residue_p->GetAtom(torsion_2.at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_1 && at_2 && at_3) {
	 have_2 = true;
	 p21 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_2.bond_length,
						    clipper::Util::d2rad(torsion_2.angle_deg),
						    clipper::Util::d2rad(torsion_2.torsion_deg));
	 p22 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_2.bond_length,
						     clipper::Util::d2rad(torsion_2.angle_deg),
						     clipper::Util::d2rad(torsion_2.torsion_deg + 120));
	 p23 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_2.bond_length,
						     clipper::Util::d2rad(torsion_2.angle_deg),
						     clipper::Util::d2rad(torsion_2.torsion_deg - 120));
      }

      // can construct result differently if we only have either p1 or p2.  THis will do for now.
      // 
      if (have_1) {
	 // this may make them too short. Hmm.
// 	 clipper::Coord_orth pav_1 = 0.5 * (p11 + p21);
// 	 clipper::Coord_orth pav_2 = 0.5 * (p12 + p22);
// 	 clipper::Coord_orth pav_3 = 0.5 * (p13 + p23);
	 
	 clipper::Coord_orth pav_1 = p11;
	 clipper::Coord_orth pav_2 = p12;
	 clipper::Coord_orth pav_3 = p13;
	 mmdb::realtype bf = at_3->tempFactor;
	 add_hydrogen_atom(at_name_1, pav_1, bf, residue_p);
	 add_hydrogen_atom(at_name_2, pav_2, bf, residue_p);
	 add_hydrogen_atom(at_name_3, pav_3, bf, residue_p);
      }
   }
}
void 
coot::reduce::add_methyl_Hs(const std::string &at_name_1,  // HB1 (for example)
			    const std::string &at_name_2,  // HB2 + 120 degress
			    const std::string &at_name_3,  // HB3 - 120 degree
			    torsion_info_t torsion_1,
			    mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      clipper::Coord_orth p11;
      clipper::Coord_orth p12;
      clipper::Coord_orth p13;
      bool have_1 = false;
      mmdb::Atom *at_1 = residue_p->GetAtom(torsion_1.at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_2 = residue_p->GetAtom(torsion_1.at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_3 = residue_p->GetAtom(torsion_1.at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_1 && at_2 && at_3) {
	 have_1 = true;
	 p11 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_1.bond_length,
						     clipper::Util::d2rad(torsion_1.angle_deg),
						     clipper::Util::d2rad(torsion_1.torsion_deg));
	 p12 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_1.bond_length,
						     clipper::Util::d2rad(torsion_1.angle_deg),
						     clipper::Util::d2rad(torsion_1.torsion_deg + 120));
	 p13 = position_by_bond_length_angle_torsion(at_1, at_2, at_3,
						     torsion_1.bond_length,
						     clipper::Util::d2rad(torsion_1.angle_deg),
						     clipper::Util::d2rad(torsion_1.torsion_deg - 120));
      }
      if (have_1) {
	 
	 clipper::Coord_orth pav_1 = p11;
	 clipper::Coord_orth pav_2 = p12;
	 clipper::Coord_orth pav_3 = p13;
	 mmdb::realtype bf = at_3->tempFactor;
	 add_hydrogen_atom(at_name_1, pav_1, bf, residue_p);
	 add_hydrogen_atom(at_name_2, pav_2, bf, residue_p);
	 add_hydrogen_atom(at_name_3, pav_3, bf, residue_p);
      }
   }
}


void
coot::reduce::add_2_sp3_hydrogens(const std::string &H_at_name_1,
				  const std::string &H_at_name_2,
				  const std::string &at_name_1,
				  const std::string &at_name_2,
				  const std::string &at_name_3,
				  double bond_length,
				  double angle_between_Hs, // in degrees
				  mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      
      mmdb::Atom *at_1 = residue_p->GetAtom(at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_2 = residue_p->GetAtom(at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_3 = residue_p->GetAtom(at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_1 && at_2 && at_3) {
	 std::pair<clipper::Coord_orth, clipper::Coord_orth> Hs =
	    position_pair_by_bisection(at_1, at_2, at_3, bond_length,
				       clipper::Util::d2rad(angle_between_Hs));
	 mmdb::realtype bf = at_2->tempFactor;
	 add_hydrogen_atom(H_at_name_1, Hs.first,  bf, residue_p);
	 add_hydrogen_atom(H_at_name_2, Hs.second, bf, residue_p);
      } else {
	 if (!alt_confs[i].empty()) {
	    std::cout << "Residue " << residue_spec_t(residue_p) << " " << residue_p->GetResName()
		      << " alt-conf \"" << alt_confs[i] << "\"" << std::endl;
	    std::cout << "Fail to add " << H_at_name_1 << " " << H_at_name_2 << " at_1: " << at_1 << std::endl;
	    std::cout << "            " << H_at_name_1 << " " << H_at_name_2 << " at_2: " << at_2 << std::endl;
	    std::cout << "            " << H_at_name_1 << " " << H_at_name_2 << " at_3: " << at_3 << std::endl;
	 }
      }
   }
}


void
coot::reduce::add_tetrahedral_hydrogen(const std::string &H_at_name,
				       const std::string &at_central_name,
				       const std::string &neighb_at_name_1,
				       const std::string &neighb_at_name_2,
				       const std::string &neighb_at_name_3,
				       double bond_length,
				       mmdb::Residue *residue_p) {

   
   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      mmdb::Atom *at_central = residue_p->GetAtom(at_central_name.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_1 = residue_p->GetAtom(neighb_at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_2 = residue_p->GetAtom(neighb_at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_3 = residue_p->GetAtom(neighb_at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_central && at_n_1 && at_n_2 && at_n_3) {
	 clipper::Coord_orth H_pos = position_by_tetrahedron(at_central, at_n_1, at_n_2, at_n_3,
							     bond_length);
	 mmdb::realtype bf = at_central->tempFactor;
	 add_hydrogen_atom(H_at_name, H_pos,  bf, residue_p);
      }
   }
}


void
coot::reduce::add_aromatic_hydrogen(const std::string &H_at_name,
				    const std::string &neighb_at_name_1,
				    const std::string &neighb_at_name_2, // add to this
				    const std::string &neighb_at_name_3,
				    double bl, mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      mmdb::Atom *at_n_1 = residue_p->GetAtom(neighb_at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_2 = residue_p->GetAtom(neighb_at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_3 = residue_p->GetAtom(neighb_at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_n_1 && at_n_2 && at_n_3) {
	 mmdb::realtype bf = at_n_2->tempFactor;
	 clipper::Coord_orth H_pos = position_by_bisection(at_n_1, at_n_2, at_n_3, bl);
	 add_hydrogen_atom(H_at_name, H_pos, bf, residue_p);
      } else {
	 std::cout << "Fail Residue " << residue_spec_t(residue_p) << " " << residue_p->GetResName()
		   << " alt-conf \"" << alt_confs[i] << "\""
		   << " failed in add_aromatic_hydrogen " << std::endl;
	 std::cout << "Fail to add " << neighb_at_name_1 << " at_1: " << at_n_1 << std::endl;
	 std::cout << "            " << neighb_at_name_2 << " at_2: " << at_n_2 << std::endl;
	 std::cout << "            " << neighb_at_name_3 << " at_3: " << at_n_3 << std::endl;
      }
   }
}

void
coot::reduce::add_amino_hydrogens(const std::string &H_at_name_1,
				  const std::string &H_at_name_2,
				  const std::string &at_name_1,   // NE2
				  const std::string &at_name_2,   // CG
				  const std::string &at_name_3,   // OE1
				  double bl_amino,
				  mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      mmdb::Atom *at_n_1 = residue_p->GetAtom(at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_2 = residue_p->GetAtom(at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_3 = residue_p->GetAtom(at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_n_1 && at_n_2 && at_n_3) {
	 clipper::Coord_orth Hp1 = position_by_bond_length_angle_torsion(at_n_3, at_n_2, at_n_1,
									 bl_amino,
									 clipper::Util::d2rad(120),
									 clipper::Util::d2rad(180));
	 clipper::Coord_orth Hp2 = position_by_bond_length_angle_torsion(at_n_3, at_n_2, at_n_1,
									 bl_amino,
									 clipper::Util::d2rad(120),
									 clipper::Util::d2rad(0));
	 mmdb::realtype bf = at_n_1->tempFactor;
	 add_hydrogen_atom(H_at_name_1, Hp1, bf, residue_p);
	 add_hydrogen_atom(H_at_name_2, Hp2, bf, residue_p);
      } else {
	 std::cout << "Fail Residue " << residue_spec_t(residue_p) << " " << residue_p->GetResName()
		   << " alt-conf \"" << alt_confs[i] << "\""
		   << " failed in add_amino_hydrogens" << std::endl;
	 std::cout << "Fail to add " << at_name_1 << " at_1: " << at_name_1 << " " << at_n_1 << std::endl;
	 std::cout << "            " << at_name_2 << " at_2: " << at_name_2 << " " << at_n_2 << std::endl;
	 std::cout << "            " << at_name_3 << " at_3: " << at_name_3 << " " << at_n_3 << std::endl;
      }
   }
}

void
coot::reduce::add_guanidinium_hydrogens(mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {

      // HE
      std::string H_at_name = " HE ";
      double bl = 0.86;
      mmdb::Atom *at_n_1 = residue_p->GetAtom(" CD ", 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_2 = residue_p->GetAtom(" NE ", 0, alt_confs[i].c_str());
      mmdb::Atom *at_n_3 = residue_p->GetAtom(" CZ ", 0, alt_confs[i].c_str());
      if (at_n_1 && at_n_2 && at_n_3) {
	 mmdb::realtype bf = at_n_2->tempFactor;
	 clipper::Coord_orth H_pos = position_by_bisection(at_n_1, at_n_2, at_n_3, bl);
	 add_hydrogen_atom(H_at_name, H_pos, bf, residue_p);
      } else {
	 std::cout << "Fail Residue " << residue_spec_t(residue_p) << " " << residue_p->GetResName()
		   << " alt-conf \"" << alt_confs[i] << "\""
		   << " failed in add_aromatic_hydrogen " << std::endl;
	 std::cout << "Fail to add guanidinium-H " << " CD " << " at_1: " << at_n_1 << std::endl;
	 std::cout << "                          " << " NE " << " at_2: " << at_n_2 << std::endl;
	 std::cout << "                          " << " CZ " << " at_3: " << at_n_3 << std::endl;
      }

      // HH[12][12]
      //
      at_n_1 = residue_p->GetAtom(" NE", 0, alt_confs[i].c_str());
      at_n_2 = residue_p->GetAtom(" CZ ", 0, alt_confs[i].c_str());
      mmdb::Atom *at_nh1 = residue_p->GetAtom(" NH1", 0, alt_confs[i].c_str());
      mmdb::Atom *at_nh2 = residue_p->GetAtom(" NH2", 0, alt_confs[i].c_str());
      if (at_n_1 && at_n_2 && at_n_3) {
	 double bf_nh1 = at_nh2->tempFactor;
	 double bf_nh2 = at_nh2->tempFactor;
	 double a = clipper::Util::d2rad(120);
	 double t = clipper::Util::d2rad(180);
	 clipper::Coord_orth hh11 = position_by_bond_length_angle_torsion(at_n_1, at_n_2, at_nh1, bl, a, 0);
	 clipper::Coord_orth hh12 = position_by_bond_length_angle_torsion(at_n_1, at_n_2, at_nh1, bl, a, t);
	 clipper::Coord_orth hh21 = position_by_bond_length_angle_torsion(at_n_1, at_n_2, at_nh2, bl, a, 0);
	 clipper::Coord_orth hh22 = position_by_bond_length_angle_torsion(at_n_1, at_n_2, at_nh2, bl, a, t);
	 add_hydrogen_atom("HH11", hh11, bf_nh1, residue_p);
	 add_hydrogen_atom("HH12", hh12, bf_nh2, residue_p);
	 add_hydrogen_atom("HH21", hh21, bf_nh2, residue_p);
	 add_hydrogen_atom("HH22", hh22, bf_nh2, residue_p);
      }
   }
}

void
coot::reduce::add_trp_indole_hydrogens(mmdb::Residue *residue_p) {

   double bl = 0.86; // H on N
   double bl_arom = 0.93;
   add_trp_indole_hydrogen(" HD1", " CG ", " CD1", " NE1", bl, residue_p);
   add_trp_indole_hydrogen(" HE1", " CD1", " NE1", " CE2", bl_arom, residue_p);
   add_trp_indole_hydrogen(" HE3", " CD2", " CE3", " CZ3", bl_arom, residue_p);
   add_trp_indole_hydrogen(" HZ3", " CE3", " CZ3", " CH2", bl_arom, residue_p);
   add_trp_indole_hydrogen(" HH2", " CZ3", " CH2", " CZ2", bl_arom, residue_p);
   add_trp_indole_hydrogen(" HZ2", " CH2", " CZ2", " CE2", bl_arom, residue_p);
}

void
coot::reduce::add_trp_indole_hydrogen(const std::string &H_name,
				      const std::string &at_name_1,
				      const std::string &at_name_2,
				      const std::string &at_name_3,
				      double bl,
				      mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      mmdb::Atom *at_1 = residue_p->GetAtom(at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_2 = residue_p->GetAtom(at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_3 = residue_p->GetAtom(at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_1 && at_2 && at_3) {
	 clipper::Coord_orth H_pos = position_by_bisection(at_1, at_2, at_3, bl);
	 double bf = at_2->tempFactor;
	 add_hydrogen_atom(H_name, H_pos, bf, residue_p);
      }
   }
}

// this will need a spin-search
void
coot::reduce::add_OH_H(const std::string &H_name,
		       const std::string &at_name_1,  // OG
		       const std::string &at_name_2,  // CB
		       const std::string &at_name_3,  // CA
		       double bl,
		       double angle,      // deg
		       double tor_inital, // deg
		       mmdb::Residue *residue_p) {

   add_xH_H(H_name, at_name_1, at_name_2, at_name_3, bl, angle, tor_inital, residue_p);

}

// this will need a spin-search
void
coot::reduce::add_SH_H(const std::string &H_name,
		       const std::string &at_name_1,  // OG
		       const std::string &at_name_2,  // CB
		       const std::string &at_name_3,  // CA
		       double bl,
		       double angle,      // deg
		       double tor_inital, // deg
		       mmdb::Residue *residue_p) {

   add_xH_H(H_name, at_name_1, at_name_2, at_name_3, bl, angle, tor_inital, residue_p);

}


// this will need a spin-search
void
coot::reduce::add_xH_H(const std::string &H_name,
		       const std::string &at_name_1,  // OG
		       const std::string &at_name_2,  // CB
		       const std::string &at_name_3,  // CA
		       double bl,
		       double angle,      // deg
		       double tor_inital, // deg
		       mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      mmdb::Atom *at_1 = residue_p->GetAtom(at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_2 = residue_p->GetAtom(at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_3 = residue_p->GetAtom(at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_1 && at_2 && at_3) {
	 clipper::Coord_orth H_pos = position_by_bond_length_angle_torsion(at_3, at_2, at_1, bl,
									   clipper::Util::d2rad(angle),
									   clipper::Util::d2rad(tor_inital));
	 double bf = at_2->tempFactor;
	 add_hydrogen_atom(H_name, H_pos, bf, residue_p);
      }
   }
}

void
coot::reduce::add_his_ring_C_Hs(mmdb::Residue *residue_p) {

   double bl_arom = 0.93;
   add_his_ring_H(" HD2", " CG ", " CD2", "NE2", bl_arom, residue_p);
   add_his_ring_H(" HE2", " ND1", " CE1", "NE2", bl_arom, residue_p);

}

void
coot::reduce::add_his_ring_H(const std::string &H_name,
			     const std::string &at_name_1,
			     const std::string &at_name_2,
			     const std::string &at_name_3,
			     double bl_arom,
			     mmdb::Residue *residue_p) {

   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      mmdb::Atom *at_1 = residue_p->GetAtom(at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_2 = residue_p->GetAtom(at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_3 = residue_p->GetAtom(at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_1 && at_2 && at_3) {
	 clipper::Coord_orth H_pos = position_by_bisection(at_1, at_2, at_3, bl_arom);
	 double bf = at_2->tempFactor;
	 add_hydrogen_atom(H_name, H_pos, bf, residue_p);
      }
   }

}


#include "atom-overlaps.hh"

void
coot::reduce::find_best_his_protonation_orientation(mmdb::Residue *residue_p) {

   // Put the H on the ND1 then on the NE2 - see which gives the best
   // score - choose that.
   //
   // Do orientation hypotheses later.

   if (geom_p) {
      std::string res_name = residue_p->GetResName();
      if (res_name == "HIS") {
	 double bl = 0.86;
	 add_his_ring_H(" HE2", " CE1", "NE2", " CD2", bl, residue_p);
	 std::vector<mmdb::Residue *> neighbs = coot::residues_near_residue(residue_p, mol, 5);
	 atom_overlaps_container_t ao(residue_p, neighbs, mol, geom_p, 0.5);
	 atom_overlaps_dots_container_t aod = ao.contact_dots_for_ligand();
      }
   } else {
      std::cout << "No geometry" << std::endl;
   }

}

void
coot::reduce::delete_atom_by_name(const std::string &at_name, mmdb::Residue *residue_p) {

   bool an_atom_was_deleted = true; // so we can start the while loop
   while (an_atom_was_deleted) {
      an_atom_was_deleted = false;
      int n_atoms = residue_p->GetNumberOfAtoms();
      for (int iat=0; iat<n_atoms; iat++) {
	 mmdb::Atom *at = residue_p->GetAtom(iat);
	 std::string ele(at->element);
	 if (ele == " H" || ele == " D") {
	    residue_p->DeleteAtom(iat);
	    an_atom_was_deleted = true;
	    break;
	 }
      }
   }
}

void
coot::reduce::delete_hydrogen_atoms() {

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
	 mmdb::Chain *chain_p;
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::Residue *residue_p;
	    mmdb::Atom *at;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       int n_atoms = residue_p->GetNumberOfAtoms();
	       bool an_atom_was_deleted = true; // so we can start the while loop
	       while (an_atom_was_deleted) {
		  an_atom_was_deleted = false;
		  for (int iat=0; iat<n_atoms; iat++) {
		     at = residue_p->GetAtom(iat);
		     std::string ele(at->element);
		     if (ele == " H" || ele == " D") {
			residue_p->DeleteAtom(iat);
			an_atom_was_deleted = true;
			break;
		     }
		  }
	       }
	    }
	 }
      }
   }
}
