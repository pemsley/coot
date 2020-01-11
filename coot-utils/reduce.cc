
#include <algorithm>
#include <iomanip>
#include <string.h>
#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"
#include "reduce.hh"
#include "atom-overlaps.hh"

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
	    residue_p = chain_p->GetResidue(ires);
	    if (ires > 0)
	       residue_prev_p = chain_p->GetResidue(ires-1);
	    else
	       residue_prev_p = 0;
	    // what about strange missing residues - where we can place the CA HA
	    // (but not the N's H).
	    bool done = add_riding_hydrogens(residue_p, residue_prev_p);
	    if (! done) {
	       hydrogen_placement_by_dictionary(residue_p);
	    } else {
	       // if this was a conventional residue, then if this was the N-terminus, we
	       // want to ad NH3+ hydrogens too.
	       if (ires==0) {
		  double bl_amino = 0.86; // add 0.03 (0.89) to match richardson reduce length. Curious
		  torsion_info_t ti(" C  ", " CA ", " N  ", bl_amino, 109, 180);
		  add_methyl_Hs(" H1 ", " H2 ", " H3 ", ti, residue_p); // not methyl
	       }
	    }
	 }
      }
   }

   mol->FinishStructEdit();

   // spin methyls and hydroxyls - work on cliques

   spinables.cliquize();
   // std::vector<std::vector<atom_with_attached_Hs> > cliques = spinables.get_cliques();

   // debug
   if (true) {
      std::cout << "--------------------------- " << spinables.cliques.size() << " cliques ----------------"
		<< std::endl;
      for (std::size_t icl=0; icl<spinables.cliques.size(); icl++) {
	 if (spinables.cliques[icl].size() > 1) {
	    std::cout << "    " << icl << " --- " << std::endl;
	    for (std::size_t j=0; j<spinables.cliques[icl].size(); j++) {
	       const atom_with_attached_Hs &typed_atom = spinables.cliques[icl][j];
	       std::cout << "      " << atom_spec_t(typed_atom.at) << std::endl;
	    }
	 }
      }
   }

   // Now, which N has the H on the HISs?
   
   model_p = mol->GetModel(imod);
   if (model_p) {
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (chain_p) {
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       if (residue_p) {
		  std::string res_name = residue_p->GetResName();
		  if (res_name == "HIS") {
		     find_best_his_protonation_orientation(residue_p);
		  }
	       }
	    }
	 }
      }
   }
}

// only call this with a valid residue_p
//
// return a status indication that this residue type was handled
//
bool
coot::reduce::add_riding_hydrogens(mmdb::Residue *residue_p, mmdb::Residue *residue_prev_p) {

   bool done = false;
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
      done = true;
   }
   if (res_name == "CYS") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " SG ", bl, 107, residue_p);
      add_SH_H(" HG ", " SG ", " CB ", " CA ", bl_sh, 109.5, 180, residue_p);
      done = true;
   }
   if (res_name == "ASP") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      done = true;
   }
   if (res_name == "GLU") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      done = true;
   }
   if (res_name == "PHE") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_aromatic_hydrogen(" HD1", " CG ", " CD1", " CE1", bl_arom, residue_p);
      add_aromatic_hydrogen(" HE1", " CD1", " CE1", " CZ ", bl_arom, residue_p);
      add_aromatic_hydrogen(" HD2", " CG ", " CD2", " CE2", bl_arom, residue_p);
      add_aromatic_hydrogen(" HE2", " CD2", " CE2", " CZ ", bl_arom, residue_p);
      add_aromatic_hydrogen(" HZ ", " CE1", " CZ ", " CE2", bl_arom, residue_p);
      done = true;
   }
   if (res_name == "GLY") {
      add_main_chain_hydrogens(residue_p, residue_prev_p, true);
      done = true;
   }
   if (res_name == "HIS") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_his_ring_C_Hs(residue_p);
      done = true;
   }
   if (res_name == "ILE") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens("HG12", "HG13", " CB ", " CG1", " CD1", bl, 107, residue_p);
      torsion_info_t ti(" CB ", " CG1", " CD1", bl, 109, 180);
      add_methyl_Hs("HD11", "HD12", "HD13", ti, residue_p);
      torsion_info_t t2(" CA ", " CB ", " CG2", bl, 109, 180);
      add_methyl_Hs("HG21", "HG22", "HG23", t2, residue_p);
      add_tetrahedral_hydrogen(" HB ", " CB ", " CA ", " CG1", " CG2", bl, residue_p);
      done = true;
   }
   if (res_name == "LYS") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HD2", " HD3", " CG ", " CD ", " CE ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HE2", " HE3", " CD ", " CE ", " NZ ", bl, 107, residue_p);
      torsion_info_t ti(" CD ", " CE ", " NZ ", bl, 109, 180);
      add_methyl_Hs(" HZ1", " HZ2", " HZ3", ti, residue_p);
      done = true;
   }
   if (res_name == "LEU") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      torsion_info_t t1(" CB ", " CG ", " CD1", bl, 109, 180);
      torsion_info_t t2(" CB ", " CG ", " CD2", bl, 109, 180);
      add_methyl_Hs("HD11", "HD12", "HD13", t1, residue_p);
      add_methyl_Hs("HD21", "HD22", "HD23", t2, residue_p);
      add_tetrahedral_hydrogen(" HG ", " CG ", " CB ", " CD1", " CD2", bl, residue_p);
      done = true;
   }
   if (res_name == "MET") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " SD ", bl, 107, residue_p);
      torsion_info_t t1(" CG ", " SD ", " CE", bl, 109, 180);
      add_methyl_Hs(" HE1", " HE2", " HE3", t1, residue_p);
      done = true;
   }
   if (res_name == "MSE") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " SE ", bl, 107, residue_p);
      torsion_info_t t1(" CG ", " SE ", " CE", bl, 109, 180);
      add_methyl_Hs(" HE1", " HE2", " HE3", t1, residue_p);
      done = true;
   }
   if (res_name == "ASN") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_amino_hydrogens("HD21", "HD22", " ND2", " CG ", " OD1", bl_amino, residue_p);
      done = true;
   }
   if (res_name == "PRO") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HD2", " HD3", " CG ", " CD ", " N  ", bl, 107, residue_p);
      done = true;
   }
   if (res_name == "GLN") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      add_amino_hydrogens("HE21", "HE22", " NE2", " CD ", " OE1", bl_amino, residue_p);
      done = true;
   }
   if (res_name == "ARG") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HG2", " HG3", " CB ", " CG ", " CD ", bl, 107, residue_p);
      add_2_sp3_hydrogens(" HD2", " HD3", " CG ", " CD ", " NE ", bl, 107, residue_p);
      add_guanidinium_hydrogens(residue_p);
      done = true;
   }
   if (res_name == "SER") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " OG ", bl, 107, residue_p);
      add_OH_H(" HG ", " OG ", " CB ", " CA ", bl_oh, 109.5, 180, residue_p);
      done = true;
   }
   if (res_name == "THR") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_tetrahedral_hydrogen(" HB ", " CB ", " CA ", " OG1", " CG2", bl, residue_p);
      torsion_info_t ti(" CA ", " CB ", " CG2", bl, 109, 180);
      add_methyl_Hs("HG21", "HG22", "HG23", ti, residue_p);
      add_OH_H(" HG1", " OG1", " CB ", " CA ", bl_oh, 109.5, 180, residue_p);
      done = true;
   }
   if (res_name == "VAL") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      torsion_info_t t1(" CA ", " CB ", " CG1 ", bl, 109, 180);
      torsion_info_t t2(" CA ", " CB ", " CG2 ", bl, 109, 180);
      add_methyl_Hs("HG11", "HG12", "HG13", t1, residue_p);
      add_methyl_Hs("HG21", "HG22", "HG23", t2, residue_p);
      add_tetrahedral_hydrogen(" HB ", " CB ", " CA ", " CG1", " CG2", bl, residue_p);
      done = true;
   }
   if (res_name == "TRP") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_trp_indole_hydrogens(residue_p);
      done = true;
   }
   if (res_name == "TYR") {
      add_main_chain_hydrogens(residue_p, residue_prev_p);
      add_2_sp3_hydrogens(" HB2", " HB3", " CA ", " CB ", " CG ", bl, 107, residue_p);
      add_aromatic_hydrogen(" HD1", " CG ", " CD1", " CE1", bl_arom, residue_p);
      add_aromatic_hydrogen(" HE1", " CD1", " CE1", " CZ ", bl_arom, residue_p);
      add_aromatic_hydrogen(" HD2", " CG ", " CD2", " CE2", bl_arom, residue_p);
      add_aromatic_hydrogen(" HE2", " CD2", " CE2", " CZ ", bl_arom, residue_p);
      add_OH_H(" HH ", " OH ", " CZ ", " CE1", bl_oh, 109.5, 180, residue_p);
      done = true;
   }
   return done;
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
	       add_hydrogen_atom(" H  ", H_pos, bf, alt_confs[i], residue_p);
	    }
	 }
      }
   }
}

// like above but for ligands (both second neighbours come from ligand residue)
void
coot::reduce::add_amino_single_H(const std::string &H_at_name,
				 const std::string &first_neighb,
				 const std::vector<std::string> &second_neighb_vec,
				 double bl,
				 mmdb::Residue *residue_p) {

   if (second_neighb_vec.size() == 2) {
      add_amino_single_H(H_at_name, second_neighb_vec[0], first_neighb, second_neighb_vec[1], bl, residue_p);
   }
}

// add H to second atom by bisection
void
coot::reduce::add_amino_single_H(const std::string H_at_name,
				 const std::string &at_name_1,
				 const std::string &at_name_2,
				 const std::string &at_name_3,
				 double bl,
				 mmdb::Residue *residue_p) {

   add_trp_indole_hydrogen(H_at_name, at_name_1, at_name_2, at_name_3, bl, residue_p);
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
	 add_hydrogen_atom(" HA ", pos, bf, alt_confs[i], residue_p);
      }
   }
}

mmdb::Atom *
coot::reduce::add_hydrogen_atom(std::string atom_name, clipper::Coord_orth &pos,
				mmdb::realtype bf,
				const std::string &altconf,
				mmdb::Residue *residue_p) {

   mmdb::Atom *new_H = new mmdb::Atom;
   new_H->SetAtomName(atom_name.c_str());
   new_H->SetElementName(" H"); // PDBv3 FIXME
   new_H->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0, bf);
   if (! altconf.empty())
      strncpy(new_H->altLoc, altconf.c_str(), 18); // 19 is mmdb limit, I think

   // now test if the atom is there already.
   //
   // It it was, then modify the coords and if not, then add it as before

   int n_atoms = residue_p->GetNumberOfAtoms();
   mmdb::Atom **residue_atoms = 0;
   mmdb::Atom *at = 0;
   bool already_exits = 0;
   residue_p->GetAtomTable(residue_atoms, n_atoms);
   for (int i=0; i<n_atoms; i++) {
      std::string residue_atom_name = residue_atoms[i]->name;
      std::string residue_atom_alt_conf  = residue_atoms[i]->altLoc;
      if (residue_atom_name == atom_name) {
	 if (residue_atom_alt_conf == altconf) {
	    already_exits = true;
	    at = residue_atoms[i];
	    break;
	 }
      }
   }

   if (! already_exits) {
      residue_p->AddAtom(new_H);
      return new_H;
   } else {
      delete new_H;
      at->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0, bf);
      return at;
   }
}

// this is also called for the hydrogens on a LYS NZ.
//
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
	 mmdb::Atom *at_0 = add_hydrogen_atom(at_name_1, pav_1, bf, alt_confs[i], residue_p);
	 mmdb::Atom *at_1 = add_hydrogen_atom(at_name_2, pav_2, bf, alt_confs[i], residue_p);
	 mmdb::Atom *at_2 = add_hydrogen_atom(at_name_3, pav_3, bf, alt_confs[i], residue_p);
	 std::vector<mmdb::Atom *> h_atoms(3);
	 h_atoms[0] = at_0;
	 h_atoms[1] = at_1;
	 h_atoms[2] = at_2;
	 spinables.add(at_3, atom_with_attached_Hs::METHYL, h_atoms);
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
	 mmdb::Atom *at_0 = add_hydrogen_atom(at_name_1, pav_1, bf, alt_confs[i], residue_p);
	 mmdb::Atom *at_1 = add_hydrogen_atom(at_name_2, pav_2, bf, alt_confs[i], residue_p);
	 mmdb::Atom *at_2 = add_hydrogen_atom(at_name_3, pav_3, bf, alt_confs[i], residue_p);
	 std::vector<mmdb::Atom *> h_atoms(3);
	 h_atoms[0] = at_0;
	 h_atoms[1] = at_1;
	 h_atoms[2] = at_2;
	 spinables.add(at_3, atom_with_attached_Hs::METHYL, h_atoms);
      }
   }
}


// choose_only_farthest_position is default false.
// if choose_only_farthest_position is true, only add the H_at_name_1 and
// add it in the position that is furthese from the averge of the second neighbours
// e.g. the single H on the N of piperidine.
//
void
coot::reduce::add_2_sp3_hydrogens(const std::string &H_at_name_1,
				  const std::string &H_at_name_2,
				  const std::string &at_name_1,
				  const std::string &at_name_2,
				  const std::string &at_name_3,
				  double bond_length,
				  double angle_between_Hs, // in degrees
				  mmdb::Residue *residue_p,
				  bool choose_only_farthest_position) {

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
	 if (! choose_only_farthest_position) {
	    add_hydrogen_atom(H_at_name_1, Hs.first,  bf, alt_confs[i], residue_p);
	    add_hydrogen_atom(H_at_name_2, Hs.second, bf, alt_confs[i], residue_p);
	 } else {
	    clipper::Coord_orth at_pos_1 = co(at_1);
	    clipper::Coord_orth at_pos_3 = co(at_3);
	    clipper::Coord_orth mp(0.5 * (at_pos_1 + at_pos_3));
	    double d1 = (Hs.first).lengthsq();
	    double d2 = (Hs.second).lengthsq();
	    if (d1 > d2)
	       add_hydrogen_atom(H_at_name_1, Hs.first,  bf, alt_confs[i], residue_p);
	    else
	       add_hydrogen_atom(H_at_name_1, Hs.second, bf, alt_confs[i], residue_p);
	 }
      } else {
	 if (!alt_confs[i].empty()) {
            if (verbose_output) {
	       std::cout << "Residue " << residue_spec_t(residue_p) << " " << residue_p->GetResName()
		         << " alt-conf \"" << alt_confs[i] << "\"" << std::endl;
	       std::cout << "Fail to add " << H_at_name_1 << " " << H_at_name_2 << " at_1: " << at_1 << std::endl;
	       std::cout << "            " << H_at_name_1 << " " << H_at_name_2 << " at_2: " << at_2 << std::endl;
	       std::cout << "            " << H_at_name_1 << " " << H_at_name_2 << " at_3: " << at_3 << std::endl;
            }
	 }
      }
   }
}

// choose_only_farthest_position is default false.
// if choose_only_farthest_position is true, only add the H_at_name_1 and
// add it in the position that is furthese from the averge of the second neighbours
//
void
coot::reduce::add_2_sp3_hydrogens(const std::string &H_at_name_1,
				  const std::string &H_at_name_2,
				  const std::string &first_neighb,
				  const std::vector<std::string> &second_neighb_vec,
				  double bond_length,
				  double angle_between_Hs, // in degrees
				  mmdb::Residue *residue_p,
				  bool choose_only_farthest_position) {

   if (second_neighb_vec.size() == 2) {
      const std::string &second_1 = second_neighb_vec[0];
      const std::string &second_2 = second_neighb_vec[1];
      add_2_sp3_hydrogens(H_at_name_1, H_at_name_2, second_1, first_neighb, second_2, bond_length,
			  angle_between_Hs, residue_p, choose_only_farthest_position);
   } else {
      std::cout << "WARNING:: in add_2_sp3_hydrogens() second_neighb_vec.size() is "
		<< second_neighb_vec.size() << std::endl;
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

   // Like on a CB of VAL
   //
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
	 add_hydrogen_atom(H_at_name, H_pos,  bf, alt_confs[i], residue_p);
      }
   }
}

void
coot::reduce::add_tetrahedral_hydrogen(const std::string &H_at_name,
				       const std::string &first_neighb,
				       const std::vector<std::string> &second_neighb_vec,
				       double bl, mmdb::Residue *residue_p) {

   if (false) {
      std::cout << "atom " << first_neighb << " has " << second_neighb_vec.size()
		<< " neighbours" << std::endl;
      for (unsigned int i=0; i<second_neighb_vec.size(); i++)
	 std::cout << "   " << second_neighb_vec[i] << std::endl;
   }

   if (second_neighb_vec.size() == 3)
      add_tetrahedral_hydrogen(H_at_name, first_neighb,
			       second_neighb_vec[0],
			       second_neighb_vec[1],
			       second_neighb_vec[2],
			       bl, residue_p);
   else
      std::cout << "WARNING:: atom " << first_neighb << " had " << second_neighb_vec.size()
		<< " neighbours  (not 3)" << std::endl;
}


// This does alphiphatic hydrogens just as well, also aldehyde hydrogens
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
	 add_hydrogen_atom(H_at_name, H_pos, bf, alt_confs[i], residue_p);
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
	 add_hydrogen_atom(H_at_name_1, Hp1, bf, alt_confs[i], residue_p);
	 add_hydrogen_atom(H_at_name_2, Hp2, bf, alt_confs[i], residue_p);
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
coot::reduce::add_amino_hydrogens(const std::string &H_at_name_1,
				  const std::string &H_at_name_2,
				  const std::string &first_neighb,
				  const std::vector<std::string> &second_neighb_vec,
				  const std::map<std::string, std::vector<std::string> > &third_neighb_map,
				  double bl_amino,
				  mmdb::Residue *residue_p) {

   if (second_neighb_vec.size() > 0) {
      std::string second = second_neighb_vec[0];
      std::map<std::string, std::vector<std::string> >::const_iterator it;
      it = third_neighb_map.find(second);
      if (it != third_neighb_map.end()) {
	 std::vector<std::string> thirds = it->second;
	 if (thirds.size() > 0) {
	    const std::string &third = thirds[0];
	    add_amino_hydrogens(H_at_name_1, H_at_name_2, first_neighb, second, third, bl_amino, residue_p);
	 }
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
	 add_hydrogen_atom(H_at_name, H_pos, bf, alt_confs[i], residue_p);
      } else {
	 std::cout << "Fail Residue " << residue_spec_t(residue_p) << " " << residue_p->GetResName()
		   << " alt-conf \"" << alt_confs[i] << "\""
		   << " failed in add_guanidinium_hydrogens " << std::endl;
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
      if (at_n_1 && at_n_2 && at_nh1 && at_nh2) {
	 double bf_nh1 = at_nh2->tempFactor;
	 double bf_nh2 = at_nh2->tempFactor;
	 double a = clipper::Util::d2rad(120);
	 double t = clipper::Util::d2rad(180);
	 clipper::Coord_orth hh11 = position_by_bond_length_angle_torsion(at_n_1, at_n_2, at_nh1, bl, a, 0);
	 clipper::Coord_orth hh12 = position_by_bond_length_angle_torsion(at_n_1, at_n_2, at_nh1, bl, a, t);
	 clipper::Coord_orth hh21 = position_by_bond_length_angle_torsion(at_n_1, at_n_2, at_nh2, bl, a, 0);
	 clipper::Coord_orth hh22 = position_by_bond_length_angle_torsion(at_n_1, at_n_2, at_nh2, bl, a, t);
	 add_hydrogen_atom("HH11", hh11, bf_nh1, alt_confs[i], residue_p);
	 add_hydrogen_atom("HH12", hh12, bf_nh2, alt_confs[i], residue_p);
	 add_hydrogen_atom("HH21", hh21, bf_nh2, alt_confs[i], residue_p);
	 add_hydrogen_atom("HH22", hh22, bf_nh2, alt_confs[i], residue_p);
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
	 add_hydrogen_atom(H_name, H_pos, bf, alt_confs[i], residue_p);
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

      
void
coot::reduce::add_OH_H(const std::string &H_at_name,
		       const std::string &first_neighb,
		       const std::vector<std::string> &second_neighb_vec,
		       const std::map<std::string, std::vector<std::string> > &third_neighb_map,
		       double bond_length,
		       double ang_deg,
		       double torsion_deg,
		       mmdb::Residue *residue_p) {

   if (second_neighb_vec.size() > 0) {
      std::string second = second_neighb_vec[0];
      std::map<std::string, std::vector<std::string> >::const_iterator it;
      it = third_neighb_map.find(second);
      if (it != third_neighb_map.end()) {
	 std::vector<std::string> thirds = it->second;
	 if (thirds.size() > 0) {
	    std::string third = thirds[0];
	    add_OH_H(H_at_name, first_neighb, second_neighb_vec[0], third,
		     bond_length, ang_deg, torsion_deg, residue_p);
	 }
      } else {
	 std::cout << "failed to find key " << second << " in thirds map" << std::endl;
      }
   }
}


// this will need a spin-search
std::vector<mmdb::Atom *>
coot::reduce::add_SH_H(const std::string &H_name,
		       const std::string &at_name_1,  // OG
		       const std::string &at_name_2,  // CB
		       const std::string &at_name_3,  // CA
		       double bl,
		       double angle,      // deg
		       double tor_inital, // deg
		       mmdb::Residue *residue_p) {

   if (is_ss_bonded(residue_p)) {
      // don't add an H on the S
      std::vector<mmdb::Atom *> empty;
      return empty;
   } else {
      return add_xH_H(H_name, at_name_1, at_name_2, at_name_3, bl, angle, tor_inital, residue_p);
   }

}


// this will need a spin-search
std::vector<mmdb::Atom *>
coot::reduce::add_xH_H(const std::string &H_name,
		       const std::string &at_name_1,  // OG
		       const std::string &at_name_2,  // CB
		       const std::string &at_name_3,  // CA
		       double bl,
		       double angle,      // deg
		       double tor_inital, // deg
		       mmdb::Residue *residue_p) {

   std::vector<mmdb::Atom *> r;
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
	 mmdb::Atom *at = add_hydrogen_atom(H_name, H_pos, bf, alt_confs[i], residue_p);
	 r.push_back(at);
	 spinables.add(at_1, atom_with_attached_Hs::HYDROXYL, at); // maybe need SULFHYDRYL separate?
      } else {
         std::cout << " a lookup fail for " << at_name_1 << " " << at_name_2 << " " << at_name_3 << " placing " << H_name << std::endl;
      }
   }
   return r;
}

void
coot::reduce::add_his_ring_C_Hs(mmdb::Residue *residue_p) {

   double bl_arom = 0.93;
   add_his_ring_H(" HD2", " CG ", " CD2", "NE2", bl_arom, residue_p);
   add_his_ring_H(" HE1", " ND1", " CE1", "NE2", bl_arom, residue_p);

}

std::vector<mmdb::Atom *>
coot::reduce::add_his_ring_H(const std::string &H_name,
			     const std::string &at_name_1,
			     const std::string &at_name_2,
			     const std::string &at_name_3,
			     double bl_arom,
			     mmdb::Residue *residue_p) {

   std::vector<mmdb::Atom *> r;
   std::vector<std::string> alt_confs = util::get_residue_alt_confs(residue_p);
   for (unsigned int i=0; i<alt_confs.size(); i++) {
      mmdb::Atom *at_1 = residue_p->GetAtom(at_name_1.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_2 = residue_p->GetAtom(at_name_2.c_str(), 0, alt_confs[i].c_str());
      mmdb::Atom *at_3 = residue_p->GetAtom(at_name_3.c_str(), 0, alt_confs[i].c_str());
      if (at_1 && at_2 && at_3) {
	 clipper::Coord_orth H_pos = position_by_bisection(at_1, at_2, at_3, bl_arom);
	 double bf = at_2->tempFactor;
	 mmdb::Atom *at = add_hydrogen_atom(H_name, H_pos, bf, alt_confs[i], residue_p);
	 r.push_back(at);
      }
   }
   return r;
}

void
coot::reduce::add_his_ring_H(const std::string &H_at_name,
			     const std::string &first_neigh,
			     const std::vector<std::string> second_neighb_vec,
			     double bl,
			     mmdb::Residue *residue_p) {

   if (second_neighb_vec.size() == 2) {
      add_his_ring_H(H_at_name, second_neighb_vec[0], first_neigh, second_neighb_vec[1], bl, residue_p);
   }
}

void
coot::reduce::add_aromatic_hydrogen(const std::string &H_at_name,
				    const std::string &first_neigh,
				    const std::vector<std::string> second_neighb_vec,
				    double bl,
				    mmdb::Residue *residue_p) {

   if (second_neighb_vec.size() == 2) {
      add_aromatic_hydrogen(H_at_name, second_neighb_vec[0], first_neigh, second_neighb_vec[1],
			    bl, residue_p);
   }
}

void
coot::reduce::find_best_his_protonation_orientation(mmdb::Residue *residue_p) {

   // Put the H on the ND1 then on the NE2 - see which gives the best
   // score - choose that.
   //
   // Do orientation hypotheses later - maybe.

   if (geom_p) {
      std::string res_name = residue_p->GetResName();
      if (res_name == "HIS") {
	 double bl = 0.86;
	 std::vector<mmdb::Atom *> v = add_his_ring_H(" HE2", " CE1", "NE2", " CD2", bl, residue_p);
	 std::vector<mmdb::Residue *> neighbs = coot::residues_near_residue(residue_p, mol, 5);
	 atom_overlaps_container_t ao_1(residue_p, neighbs, mol, geom_p, 0.5);
	 atom_overlaps_dots_container_t aod_1 = ao_1.contact_dots_for_ligand();
	 double s1 = aod_1.score();
	 // this only does the first alt conf. It can get messy with alt confs.

	 if (v.size() > 0) {
	    mmdb::Atom *at = v[0];
	    delete at;
	    mol->FinishStructEdit();
	 }

 	 v = add_his_ring_H(" HD1", " CG ", "ND1", " CE1", bl, residue_p);
 	 atom_overlaps_container_t ao_2(residue_p, neighbs, mol, geom_p, 0.5);
 	 atom_overlaps_dots_container_t aod_2 = ao_2.contact_dots_for_ligand();
 	 double s2 = aod_2.score();
 	 std::cout << "DEBUG:: HIS protonation scores: " << residue_spec_t(residue_p)
		   << " NE2: "
		   << std::right << std::setprecision(1) << std::fixed << s1
		   << " vs ND1: "
		   << std::right << std::setprecision(1) << std::fixed << s2 << std::endl;
	 if (v.size() > 0) { // sanity check
	    if (s1 > s2) {
	       // delete HD1
	       delete v[0];
	       // add the first one back
	       add_his_ring_H(" HE2", " CE1", "NE2", " CD2", bl, residue_p);
	       mol->FinishStructEdit();
	    }
	 }
      }
   } else {
      std::cout << "No geometry" << std::endl;
   }
}

void
coot::reduce::switch_his_protonation(mmdb::Residue *residue_p,
				     mmdb::Atom *current_H_atom) {

   if (current_H_atom) {
      std::string atom_name = current_H_atom->name;
      std::string new_atom_name;
      if (atom_name == " HD1") new_atom_name = " HE2";
      if (atom_name == " HE2") new_atom_name = " HD1";
      if (! new_atom_name.empty()) {

	 // ND1 -> HD1   and   NE2 -> HE2 by dictionary
	 //
	 // New atom HD1 first
	 std::cout << "switch_his_protonation() " << 1 << std::endl;
	 std::string at_name_1 = " CG ";
	 std::string at_name_2 = " ND1";
	 std::string at_name_3 = " CE1";
	 if (new_atom_name == " HE2") {
	    at_name_1 = " CE1";
	    at_name_2 = " NE2";
	    at_name_3 = " CD2";
	 }
	 std::string alt_conf = current_H_atom->altLoc;
	 mmdb::Atom *at_1 = residue_p->GetAtom(at_name_1.c_str(), 0, alt_conf.c_str());
	 mmdb::Atom *at_2 = residue_p->GetAtom(at_name_2.c_str(), 0, alt_conf.c_str());
	 mmdb::Atom *at_3 = residue_p->GetAtom(at_name_3.c_str(), 0, alt_conf.c_str());
	 if (at_1 && at_2 && at_3) {
	    std::cout << "switch_his_protonation() " << 2 << " " << new_atom_name << std::endl;
	    current_H_atom->SetAtomName(new_atom_name.c_str());
	    double bl_arom = 0.93;
	    clipper::Coord_orth pos = position_by_bisection(at_1, at_2, at_3, bl_arom);
	    double bf = current_H_atom->tempFactor;
	    current_H_atom->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0, bf);
	 }
      }
   }
}

bool
coot::reduce::hack_ss_bond_test(mmdb::Residue *CYS_1_residue_p, mmdb::Model *model_p) const {

   bool status = false;

   mmdb::Atom *CYS_1_SG = 0;
   int n_atoms = CYS_1_residue_p->GetNumberOfAtoms();
   for (int iat=0; iat<n_atoms; iat++) {
      mmdb::Atom *at = CYS_1_residue_p->GetAtom(iat);
      std::string atom_name = at->GetAtomName();
      if (atom_name == " SG ") {
         CYS_1_SG = at;
         break;
      }
   }
   if (! CYS_1_SG) return false;

   clipper::Coord_orth pt_1 = co(CYS_1_SG);
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) {
         mmdb::Residue *residue_p = chain_p->GetResidue(ires);
         std::string res_name = residue_p->GetResName();
         if (residue_p != CYS_1_residue_p) {
            if (res_name == "CYS") {
               int n_atoms = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms; iat++) {
                  mmdb::Atom *at = residue_p->GetAtom(iat);
                  std::string atom_name = at->GetAtomName();
                  if (atom_name == " SG ") {
                     clipper::Coord_orth pt_2 = co(at);
                     double dd = (pt_2-pt_1).lengthsq();
                     if (dd < 3.0 * 3.0) {
                        status = true;
                        break;
                     }
                  }
               }
            }
         }
         if (status) break;
      }
   }

   return status;
}

bool
coot::reduce::is_ss_bonded(mmdb::Residue *residue_p) const {

   bool status = false;
   if (residue_p) {
      std::string res_name = residue_p->GetResName();
      if (res_name == "CYS") {
	 int imod = 1;
	 mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
	    // check SS bonds here
            //
            // Oh dear, we can't interograte SSBONDs.
            bool i = hack_ss_bond_test(residue_p, model_p);
            if (i) status = true;
         }
      }
   }

   return status;
}

bool
coot::reduce::is_linked(const std::string &atom_name, mmdb::Residue *residue_p) const {

   bool status = false;

   std::string chain_id = residue_p->GetChainID();
   std::string ins_code = residue_p->GetInsCode();
   int res_no = residue_p->GetSeqNum();

   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
	 unsigned int n_links = model_p->GetNumberOfLinks();
	 for (unsigned int i=1; i<=n_links; i++) {
	    mmdb::Link *link = model_p->GetLink(i);
	    if (! link) {
	       std::cout << "ERROR:: null link " << i << " in ref" << std::endl;
	    } else {
	       int link_res_no_1 = link->seqNum1;
	       int link_res_no_2 = link->seqNum2;
	       std::string link_ins_code_1 = link->insCode1;
	       std::string link_ins_code_2 = link->insCode2;
	       std::string link_chain_id_1 = link->chainID1;
	       std::string link_chain_id_2 = link->chainID2;
	       std::string link_atom_name_1 = link->atName1;
	       std::string link_atom_name_2 = link->atName2;
	       if (chain_id == link_chain_id_1) {
		  if (res_no == link_res_no_1) {
		     if (ins_code == link_ins_code_1) {
			if (atom_name == link_atom_name_1) {
			   status = true;
			   break;
			}
		     }
		  }
	       }
	       if (chain_id == link_chain_id_2) {
		  if (res_no == link_res_no_2) {
		     if (ins_code == link_ins_code_2) {
			if (atom_name == link_atom_name_2) {
			   status = true;
			   break;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (false)
      std::cout << "debug:: is_linked " << coot::residue_spec_t(residue_p) << " "
		<< atom_name << " " << status << std::endl;
   return status;
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

// skip waters, HOH
void
coot::reduce::hydrogen_placement_by_dictionary(mmdb::Residue *residue_p) {

   std::string res_name = residue_p->GetResName();
   if (res_name != "HOH") {
      if (geom_p) {
	 std::pair<bool, dictionary_residue_restraints_t> p =
	    geom_p->get_monomer_restraints(res_name, imol);
	 if (p.first) {
	    hydrogen_placement_by_dictionary(p.second, residue_p);
	 }
      }
   }
}

void
coot::reduce::hydrogen_placement_by_dictionary(const dictionary_residue_restraints_t &rest,
					       mmdb::Residue *residue_p) {

   std::vector<std::string> done_atom_name_list; // so that we don't add some atoms twice
   for (unsigned int iat=0; iat<rest.atom_info.size(); iat++) {
      if (rest.atom_info[iat].is_hydrogen()) {
	 const std::string &H_at_name = rest.atom_info[iat].atom_id_4c;
	 // if we haven't done it already...
	 if (std::find(done_atom_name_list.begin(), done_atom_name_list.end(), H_at_name) == done_atom_name_list.end()) {
	    // skip the HO3' on RNA and DNA. I could instead test for presence/position of next
	    // residue, but this easier and will be correct for most cases.
	    if ((rest.residue_info.group == "DNA" || rest.residue_info.group == "RNA") &&
		H_at_name == "HO3'") {
	       continue;
	    } else {
	       // to which atom is this hydrogen connected?
	       std::vector<unsigned int> neighbs = rest.neighbours(iat, false);
	       if (neighbs.size() == 1) {
		  // what else would it be?
		  const unsigned int &iat_neighb = neighbs[0];
		  const std::string &energy_type = rest.atom_info[iat_neighb].type_energy;
		  const std::string &first_neigh = rest.atom_info[iat_neighb].atom_id_4c;
		  if (! is_linked(first_neigh, residue_p)) {
		     if (! energy_type.empty()) {
			std::vector<std::string> v =
			   place_hydrogen_by_connected_atom_energy_type(iat, iat_neighb, rest, residue_p);
			done_atom_name_list.insert(done_atom_name_list.end(), v.begin(), v.end());
		     } else {
			place_hydrogen_by_connected_2nd_neighbours(iat, iat_neighb, rest, residue_p);
		     }
		  }
	       }
	    }
	 }
      }
   }

}

// return a list of placed atoms (sometimes, e.g. NH2) placing one atom means placing two of them
//
std::vector<std::string>
coot::reduce::place_hydrogen_by_connected_atom_energy_type(unsigned int iat,
							   unsigned int iat_neighb,
							   const dictionary_residue_restraints_t &rest,
							   mmdb::Residue *residue_p) {

   std::vector<std::string> v;
   const std::string &energy_type = rest.atom_info[iat_neighb].type_energy;
   return place_hydrogen_by_connected_atom_energy_type(energy_type, iat, iat_neighb, rest, residue_p);

}


// return a list of placed atoms (sometimes, e.g. NH2) placing one atom means placing two of them
//
std::vector<std::string>
coot::reduce::place_hydrogen_by_connected_atom_energy_type(const std::string &energy_type,
							   unsigned int iat,
							   unsigned int iat_neighb,
							   const dictionary_residue_restraints_t &rest,
							   mmdb::Residue *residue_p) {

   std::vector<std::string> v;
   const std::string &H_at_name = rest.atom_info[iat].atom_id_4c;
   if (false)
      std::cout << " hydrogen atom " << H_at_name << " in residue type " << rest.residue_info.comp_id
		<< " by_energy_type: " << energy_type << std::endl;
   const std::string &first_neighb = rest.atom_info[iat_neighb].atom_id_4c;
   std::vector<std::string> second_neighb_vec = rest.neighbours(first_neighb, false);

   if (false) {
      std::cout << " place_hydrogen_by_connected_atom_energy_type second_neighb_vec.size() "
		<< second_neighb_vec.size() << std::endl;
      for (unsigned int ii=0; ii<second_neighb_vec.size(); ii++)
	 std::cout << " place_hydrogen_by_connected_atom_energy_type second_neighb_vec: "
		   << ii << " " << second_neighb_vec[ii] << std::endl;
   }

   if (energy_type == "CR16" || energy_type == "CR15") {
      double bl = 1.08;
      add_aromatic_hydrogen(H_at_name, first_neighb, second_neighb_vec, bl, residue_p);
      v.push_back(H_at_name);
   }

   if (energy_type == "NR15" || energy_type == "NR16") {
      double bl = 0.86;
      add_his_ring_H(H_at_name, first_neighb, second_neighb_vec, bl, residue_p);
      v.push_back(H_at_name);
   }

   if (energy_type == "NH2" || energy_type == "C2") {
      double bl_amino = 1.01; // makes refmac mon lib
      double bl = 0.97;
      if (energy_type == "NH2") bl = bl_amino;
      std::string H_at_other = get_other_H_name(first_neighb, H_at_name, rest);
      if (! H_at_other.empty()) {
	 std::map<std::string, std::vector<std::string> >
	    tnm = third_neighbour_map(first_neighb, second_neighb_vec, rest);
	 add_amino_hydrogens(H_at_name, H_at_other, first_neighb, second_neighb_vec, tnm,
			     bl, residue_p);
	 v.push_back(H_at_name);
	 v.push_back(H_at_other);
      }
   }

   if (energy_type == "CH1") {
      // sp3 carbon with one H atom
      double bl = 0.97;
      add_tetrahedral_hydrogen(H_at_name, first_neighb, second_neighb_vec, bl, residue_p);
      v.push_back(H_at_name);
   }

   if (energy_type == "CH2") {
      // sp3 carbon with 2 H atoms
      double bl = 0.97;
      std::string H_at_other = get_other_H_name(first_neighb, H_at_name, rest);
      if (! H_at_other.empty()) {
	 double angle_between_Hs = 107.0; // we could get this from the dictionary
	 add_2_sp3_hydrogens(H_at_name, H_at_other, first_neighb, second_neighb_vec,
			     bl, angle_between_Hs, residue_p);
	 v.push_back(H_at_name);
	 v.push_back(H_at_other);
      }
   }

   if (energy_type == "NT1") {

      bool choose_only_farthest_position = true; // we find two H positions, choose
	                                         // the one farthest from the 2nd neighbs.

      // H atom on sp3 nitrogen connected to 2 sp3 carbon atoms, e.g. in piperidine
      // choose equitorial orientation, not axial.  Which means, build both hydrogen positions
      // and chose the one that is furthest from the average position of the second neighbours.
      //
      double angle_between_Hs = 107.0; // we could get this from the dictionary
      double bl = 0.86;
      std::string H_at_other = "Hdum";
      add_2_sp3_hydrogens(H_at_name, H_at_other, first_neighb, second_neighb_vec,
			  bl, angle_between_Hs, residue_p, choose_only_farthest_position);
      v.push_back(H_at_name);
   }

   if (energy_type == "NT2") {
      double angle_between_Hs = 107.0; // we could get this from the dictionary
      double bl = 0.86;
      std::string H_at_other = get_other_H_name(first_neighb, H_at_name, rest);
      add_2_sp3_hydrogens(H_at_name, H_at_other, first_neighb, second_neighb_vec,
			  bl, angle_between_Hs, residue_p);
      v.push_back(H_at_name);
      v.push_back(H_at_other);
   }

   if (energy_type == "OH1") {
      // needs a spin-search
      double bl_oh = 0.84;
      std::map<std::string, std::vector<std::string> > tnm =
	 third_neighbour_map(first_neighb, second_neighb_vec, rest);
      add_OH_H(H_at_name, first_neighb, second_neighb_vec, tnm,
	       bl_oh, 109.5, 180, residue_p);
      v.push_back(H_at_name);
   }

   if (energy_type == "SH1") {
      // needs a spin-search
      double bl_oh = 1.2; // 1.33 from ener lib.
      std::map<std::string, std::vector<std::string> > tnm =
	 third_neighbour_map(first_neighb, second_neighb_vec, rest);
      add_OH_H(H_at_name, first_neighb, second_neighb_vec, tnm,
	       bl_oh, 109.5, 180, residue_p);
      v.push_back(H_at_name);
   }

   if (energy_type == "CH3" || energy_type == "NT3") {
      std::vector<std::string> H_other = get_other_H_names(H_at_name, first_neighb, rest);
      if (H_other.size() == 2) {
	 std::map<std::string, std::vector<std::string> > tnm =
	    third_neighbour_map(first_neighb, second_neighb_vec, rest);
	 double bl = 0.97;
	 double ang = 109;
	 double tor = 180;
	 torsion_info_t tor_info(first_neighb, second_neighb_vec, tnm, bl, ang, tor); // fn reverses order
	 add_methyl_Hs(H_at_name, H_other[0], H_other[1], tor_info, residue_p);
	 v.push_back(H_at_name);
	 v.push_back(H_other[0]);
	 v.push_back(H_other[1]);
      }
   }

   if (energy_type == "NH1") {
      // e.g. H on the N in an amino acid peptide
      double bl = 0.86;
      add_amino_single_H(H_at_name, first_neighb, second_neighb_vec, bl, residue_p);
      v.push_back(H_at_name);
   }

   if (energy_type == "C1") {
      // aldehyde Hydrogen atom
      if (second_neighb_vec.size() == 2) {
	 double bl = 0.97;
	 add_aromatic_hydrogen(H_at_name, second_neighb_vec[0], first_neighb, second_neighb_vec[1],
			       bl, residue_p);
	 v.push_back(H_at_name);
      }
   }

   if (v.empty()) {
      std::cout << "FAIL: -------- place_hydrogen_by_connected_atom_energy_type()"
		<< " H_at_name \"" << H_at_name << "\" neighb: \""
		<< first_neighb << "\" energy_type "
		<< energy_type << " for comp_id " << rest.residue_info.comp_id
		<< std::endl;
   }
   return v;
}

// atoms that are connected to the second neighbour that aren't the first neighbour.
//
std::map<std::string, std::vector<std::string> > 
coot::reduce::third_neighbour_map(const std::string &first_neighb,
				  const std::vector<std::string> second_neighb_vec,
				  const coot::dictionary_residue_restraints_t &rest) const {

   std::map<std::string, std::vector<std::string> > m;
   for (unsigned int i=0; i<second_neighb_vec.size(); i++) {
      std::vector<std::string> nn = rest.neighbours(second_neighb_vec[i], false);
      for (unsigned int j=0; j<nn.size(); j++) {
	 if (nn[j] != first_neighb)
	    m[second_neighb_vec[i]].push_back(nn[j]);
      }
   }
   return m;
}

void
coot::reduce::place_hydrogen_by_connected_2nd_neighbours(unsigned int iat,
							 unsigned int iat_neighb,
							 const dictionary_residue_restraints_t &rest,
							 mmdb::Residue *residue_p) {
   
   std::vector<unsigned int> neighbs = rest.neighbours(iat_neighb, false);
   std::string ele_neighb = rest.atom_info[iat_neighb].type_symbol;

   // something complicated
   for (unsigned int iat=0; iat<neighbs.size(); iat++) {
   }
}

// return a null string on failure.
//
std::string
coot::reduce::get_other_H_name(const std::string &first_neighb,
			       const std::string &H_at_name,
			       const dictionary_residue_restraints_t &dict) const {

   return dict.get_other_H_name(H_at_name);
}

// return a empty vector on failure
//
std::vector<std::string>
coot::reduce::get_other_H_names(const std::string &H_at_name,
				const std::string &first_neighb,
				const dictionary_residue_restraints_t &dict) const {

   return dict.get_other_H_names(H_at_name);
}

void
coot::reduce::atoms_with_spinnable_Hs::add(mmdb::Atom *at,
					   atom_with_attached_Hs::hydrogen_t type,
					   const std::vector<mmdb::Atom *> &attached_hydrogen_atoms) {

   std::string alt_loc(at->altLoc);
   atom_with_attached_Hs awaH(at, type, attached_hydrogen_atoms);
   typed_atoms[alt_loc].push_back(awaH);
}

void
coot::reduce::atoms_with_spinnable_Hs::add(mmdb::Atom *at,
					   atom_with_attached_Hs::hydrogen_t type,
					   mmdb::Atom *attahed_hydrogen_atom) {

   std::string alt_loc(at->altLoc);
   std::vector<mmdb::Atom *> v;
   v.push_back(attahed_hydrogen_atom);
   atom_with_attached_Hs awaH(at, type, v);
   typed_atoms[alt_loc].push_back(awaH);
}

void
coot::reduce::atoms_with_spinnable_Hs::cliquize() {

   // maybe the cliques need to be split on alt conf
   // std::map<std::string, std::vector<std::pair<hydrogen_t, mmdb::Atom *> > >::const_iterator it;
   std::map<std::string, std::vector<atom_with_attached_Hs> >::const_iterator it;

   if (true) {
      for(it=typed_atoms.begin(); it!=typed_atoms.end(); it++) {
	 const std::string &key = it->first;
	 std::cout << "cliquize " << typed_atoms[key].size() << " spinables for altconf "
		   << key << std::endl;
      }
   }

   const double d_crit = 4.0;

   // make this a member data item
   //
   // maybe this loop is slow for big proteins
   // rethink using the results of SeekContacts()?
   //
   // atoms in the same side-chain are not in the same clique
   //
   for(it=typed_atoms.begin(); it!=typed_atoms.end(); it++) {
      const std::string &key = it->first;
      const std::vector<atom_with_attached_Hs> &atoms = it->second;
      for (std::size_t iat=0; iat<atoms.size(); iat++) {
	 bool done = false;
	 clipper::Coord_orth at_pos = co(atoms[iat].at);
	 for (unsigned int icl=0; icl<cliques.size(); icl++) {
	    for (unsigned int j=0; j<cliques[icl].size(); j++) {
	       if (atoms[iat].at->residue != cliques[icl][j].at->residue) {
		  clipper::Coord_orth pos = co(cliques[icl][j].at);
		  clipper::Coord_orth diff(at_pos - pos);
		  double dv_sqrd = diff.lengthsq();
		  if (dv_sqrd < d_crit * d_crit) {
		     cliques[icl].push_back(atoms[iat]);
		     done = true;
		     break;
		  }
	       }
	    }
	    if (done)
	       break;
	 }
	 if (! done) {
	    // start a new clique;
	    std::vector<atom_with_attached_Hs> new_clique;
	    new_clique.push_back(atoms[iat]);
	    cliques.push_back(new_clique);
	 }
      }
   }
}


void
coot::reduce::atoms_with_spinnable_Hs::resolve_clashing_clique(const std::vector<atom_with_attached_Hs> &clique) {



}
