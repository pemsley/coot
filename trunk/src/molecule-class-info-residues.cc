/* src/molecule-class-info-residues.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009 by The University of Oxford
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

// Must include these headers to get molecule_class_info_t.h to parse.
//

#include <string>
#include <algorithm>
#include <stdexcept>

#include <GL/glut.h>

// #ifdef MAKE_ENTERPRISE_TOOLS
// // includes order important, otherwise we get dcgettext() problems.
// #include "rdkit-interface.hh" // needed for add_hydrogens()
// #endif

#include <mmdb/mmdb_manager.h>
#include "mmdb-extras.h"
#include "Cartesian.h"
#include "mmdb-crystal.h"
// 
#include "rotamer.hh" // in ligand

#include "base-pairing.hh"

#include "molecule-class-info.h"
#include "coot-utils.hh"
#include "coot-hydrogens.hh"

#include "rama_plot.hh"

// 1: success
// 0: failure
// 
bool
molecule_class_info_t::progressive_residues_in_chain_check_by_chain(const char *chain_id) const {

   bool r = 0;
   
   if (atom_sel.n_selected_atoms > 0) { 
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 CChain *chain_p = atom_sel.mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) { 
	    r = coot::progressive_residues_in_chain_check(chain_p);
	    break;
	 }
      }
   }

   return r; 
} 


// Only apply charges if the molecule contains lots of hydrogens.
//
// so return a flag, whether or not the charges were applied.
// 
bool
molecule_class_info_t::apply_charges(const coot::protein_geometry &geom) {

   // More than 15% of the atoms have to be hydrogens for use to set
   // the charges on all the atoms (to something other than
   // CXX_UNSET_CHARGE).

   bool charges_applied_flag = 0; // unset initially.
   
   float fraction_hydrogens = 0.15;

   if (atom_sel.n_selected_atoms > 0) { 
      CMMDBManager *mol = atom_sel.mol;

      int n_H = 0;
      int n_all = 0; 
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 std::string ele(atom_sel.atom_selection[i]->element);
	 if (ele == " H" || ele == " D") {
	    n_H++; 
	 }
	 n_all++;
      }
      
      if ( (float(n_H)/float(n_all) > fraction_hydrogens) || n_all < 100) {

	 // first set all atom charges to unset:
	 for (int i=0; i<atom_sel.n_selected_atoms; i++)
	    atom_sel.atom_selection[i]->charge = CXX_UNSET_CHARGE - 0.1;

	 charges_applied_flag = 1;

	 // Now add real charges from the dictionary
	 // 
	 int imod = 1;
	 CModel *model_p = mol->GetModel(imod);
	 CChain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       std::string res_type = residue_p->GetResName();
	       std::pair<short int, coot::dictionary_residue_restraints_t> rp = 
		  geom.get_monomer_restraints(res_type);
	       if (rp.first) {
		  try { 
		     coot::dipole p(rp.second, residue_p);
		     p.fill_charged_atoms(residue_p, rp.second);
		  }
		  catch (std::runtime_error mess) {
		     std::cout << mess.what() << std::endl;
		  }
	       }
	    }
	 }
      }
   }
   return charges_applied_flag;
} 

int
molecule_class_info_t::assign_hetatms() {

   int r = 0;
   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      int imod = 1;
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    r += coot::hetify_residue_atoms_as_needed(residue_p);
	 }
      }
   }
   return r;
}


std::pair<bool, std::string>
molecule_class_info_t::sprout_hydrogens(const std::string &chain_id,
					int res_no,
					const std::string &ins_code,
					const coot::protein_geometry &geom) {

   std::pair<bool, std::string> r(0, "");

   make_backup();
   CResidue *residue_p = get_residue(chain_id, res_no, ins_code);
   std::vector<coot::atom_spec_t> fixed_atoms;
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (unsigned int i=0; i<n_residue_atoms; i++)
      if (std::string(residue_atoms[i]->element) != " H")
	 fixed_atoms.push_back(residue_atoms[i]);
    
   if (residue_p) {
      std::string residue_type = residue_p->GetResName();
      std::pair<bool, coot::dictionary_residue_restraints_t> p = 
	 geom.get_monomer_restraints_at_least_minimal(residue_type);
      if (! p.first) {
	 std::cout << "WARNING:: sprout_hydrogens(): No restraints for residue type "
		   << residue_type << std::endl;
      } else {
	 r = coot::add_hydrogens(residue_p, p.second);

	 std::cout << "DEBUG:: coot::add_hydrogens() returns " << r.first
		   << " " << r.second << std::endl;
	 
	 if (r.first) {

	    std::string residue_name = residue_p->GetResName();

	    std::pair<bool, coot::dictionary_residue_restraints_t> rp = 
	       geom.get_monomer_restraints(residue_name);

	    if (rp.first) { 
	    
	       // those Hs were just attached with non-good geometry, we
	       // need to minimise.  Keep all atoms fixed except all hydrogens.
	       std::vector<std::pair<bool,CResidue *> > residues;
	       std::pair<bool, CResidue *> p(0, residue_p);
	       residues.push_back(p);
	       coot::restraints_container_t restraints(residues, geom,
						       atom_sel.mol, fixed_atoms);
	       bool do_torsions = 0;

	       // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED; // fail
	       coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
	       // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES; // test
	       
	       // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES; // pass
	       // coot::restraint_usage_Flags flags = coot::BONDS_AND_PLANES; // pass
	       // coot::restraint_usage_Flags flags = coot::BONDS_AND_NON_BONDED; // fail
	       // coot::restraint_usage_Flags flags = coot::BONDS_AND_PLANES; // pass
	       // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES; // pass
	       // coot::restraint_usage_Flags flags = coot::BONDS_AND_NON_BONDED;
	       
	       int n_restraints = restraints.make_restraints(geom, flags, do_torsions,
							     0, 0, coot::NO_PSEUDO_BONDS);
	       restraints.minimize(flags);

	       // analyse results
	       PPCAtom residue_atoms = 0;
	       int n_residue_atoms;
	       residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	       
	       // 
	       std::vector<coot::dict_chiral_restraint_t> cr = rp.second.chiral_restraint;

	       for (unsigned int icr=0; icr<cr.size(); icr++) {
		  std::string centre_atom = cr[icr].atom_id_c_4c();
		  std::vector<std::pair<short int, coot::atom_spec_t> > v = 
		     coot::is_inverted_chiral_atom_p(cr[icr], residue_p);
		  if (v.size() ) {
		     for (unsigned int i=0; i<v.size(); i++) {
			if (v[i].first) {
			   std::cout << "fix this bad chiral centre "
				     << v[i].first << " "
				     << v[i].second << std::endl;
			   std::vector<std::string> attached_Hs =
			      rp.second.get_attached_H_names(v[i].second.atom_name);
			   if (attached_Hs.size() > 1) {

			      coot::atom_spec_t spec_1 = v[i].second;
			      coot::atom_spec_t spec_2 = v[i].second;
			      spec_1.atom_name = attached_Hs[0];
			      spec_2.atom_name = attached_Hs[1];
			      CAtom *at_1 = get_atom(spec_1);
			      CAtom *at_2 = get_atom(spec_2);

			      if (! at_1) {
				 std::cout << " failed to get atom with spec " << spec_1
					   << std::endl;
			      } else {
				 if (! at_2) {
				    std::cout << " failed to get atom with spec " << spec_2
					      << std::endl;
				 } else {
				    clipper::Coord_orth t(at_1->x, at_1->y, at_1->z);
				    at_1->x = at_2->x;
				    at_1->y = at_2->y;
				    at_1->z = at_2->z;
				    at_2->x = t.x();
				    at_2->y = t.y();
				    at_2->z = t.z();
				 }
			      }
			   }
			}
		     }
		  }
	       }

	       flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
	       n_restraints = restraints.make_restraints(geom, flags, do_torsions,
							 0, 0, coot::NO_PSEUDO_BONDS);
	       restraints.minimize(flags);
	    }
	 }
      }
   }
   if (r.first) {
      have_unsaved_changes_flag = 1;
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked();
   }
   return r;
}

std::vector<std::string>
molecule_class_info_t::no_dictionary_for_residue_type_as_yet(const coot::protein_geometry &geom) const {

   std::vector<std::string> v;
   int imod = 1;
   if (atom_sel.mol) { 
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    std::string residue_name = residue_p->GetResName();
	    if (! geom.have_at_least_minimal_dictionary_for_residue_type(residue_name)) {

	       // add it to v, if it is not already there:
	       if (std::find(v.begin(), v.end(), residue_name) == v.end())
		  v.push_back(residue_name);
	    } 
	 }
      }
   }
   return v;
}


// ------------- helper function to orient_view() --------------------
// can throw an std::runtime exception;
clipper::Coord_orth
molecule_class_info_t::get_vector(const coot::residue_spec_t &central_residue_spec, // ligand typically
				  const coot::residue_spec_t &neighbour_residue_spec) const {

   clipper::Coord_orth r(0,0,0);

   CResidue *central_residue = get_residue(central_residue_spec);
   CResidue *neighbour_residue = get_residue(neighbour_residue_spec);

   if (! central_residue) {
      std::string message = "Missing residue ";
      message += central_residue_spec.chain;
      message += " "; 
      message += central_residue_spec.resno;
      throw std::runtime_error(message);
   } else { 
      if (! neighbour_residue) {
	 std::string message = "Missing residue ";
	 message += neighbour_residue_spec.chain;
	 message += " ";
	 message += neighbour_residue_spec.resno;
	 throw std::runtime_error(message);
      } else {
	 // OK! go.

	 double min_dist = 42e32;
	 clipper::Coord_orth shortest_dist(0,0,0); // "best"
	 PPCAtom c_residue_atoms = 0;
	 int c_n_residue_atoms;
	 PPCAtom n_residue_atoms = 0;
	 int n_n_residue_atoms;
	 central_residue->GetAtomTable  (c_residue_atoms, c_n_residue_atoms);
	 neighbour_residue->GetAtomTable(n_residue_atoms, n_n_residue_atoms);
	 clipper::Coord_orth sum_1(0,0,0);
	 clipper::Coord_orth sum_2(0,0,0);
	 for (unsigned int iat=0; iat<c_n_residue_atoms; iat++) {
	    if (! c_residue_atoms[iat]->isTer()) { 
	       clipper::Coord_orth pt_1(c_residue_atoms[iat]->x,
					c_residue_atoms[iat]->y,
					c_residue_atoms[iat]->z);
	       sum_1 += pt_1;
	    }
	 }
	 for (unsigned int jat=0; jat<n_n_residue_atoms; jat++) {
	    if (! n_residue_atoms[jat]->isTer()) { 
	       clipper::Coord_orth pt_2(n_residue_atoms[jat]->x,
					n_residue_atoms[jat]->y,
					n_residue_atoms[jat]->z);
	       sum_2 += pt_2;
	    }
	 }

	 if (sum_1.lengthsq() > 0.0001) {
	    if (sum_2.lengthsq() > 0.0001) {
	       double frac_1 = 1.0/double(c_n_residue_atoms);
	       double frac_2 = 1.0/double(n_n_residue_atoms);
	       r = clipper::Coord_orth(sum_2 * frac_2 - sum_1 * frac_1);
	    } else {
	       throw std::runtime_error("No atoms in residue?");
	    }
	 } else {
	    throw std::runtime_error("No atoms in residue?");
	 } 
      }
   }

   return r;
} 

void
molecule_class_info_t::match_ligand_atom_names(const std::string &chain_id, int res_no, const std::string &ins_code,
					       CResidue *res_ref) {

   CResidue *res_mov = get_residue(chain_id, res_no, ins_code);

   if (! res_mov) {
      std::cout << "No residue for moving atom names:  " << chain_id << " " << res_no << " "  << ins_code
		<< std::endl;
      
   } else {
      bool match_hydrogens_also = 1;
      coot::graph_match_info_t gm = coot::graph_match(res_mov, res_ref, 0, match_hydrogens_also);
      gm.match_names(res_mov);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   } 
} 


coot::rama_score_t
molecule_class_info_t::get_all_molecule_rama_score() const {

   coot::rama_score_t rs;

   coot::rama_plot rp;
   rp.generate_phi_psis(atom_sel.mol, true);

   unsigned int n = rp.n_phi_psi_model_sets();

   clipper::Ramachandran r_gly, r_pro, r_non_gly_pro;

   // Lovell et al. 2003, 50, 437 Protein Structure, Function and Genetics values:
   ftype level_prefered = 0.02; 
   ftype level_allowed = 0.002;
   
   //clipper defaults: 0.01 0.0005

   r_gly.init(clipper::Ramachandran::Gly);
   r_gly.set_thresholds(level_prefered, level_allowed);
   //
   r_pro.init(clipper::Ramachandran::Pro);
   r_pro.set_thresholds(level_prefered, level_allowed);
   // 
   r_non_gly_pro.init(clipper::Ramachandran::NonGlyPro);
   r_non_gly_pro.set_thresholds(level_prefered, level_allowed);

   double zero_cut_off = 1e-6;

   double log_p_sum = 0.0;
   double log_p_non_sec_str_sum = 0.0;

   for (unsigned int imod=1; imod<n; imod++) {
      if (imod<=atom_sel.mol->GetNumberOfModels()) {
	 
	 coot::phi_psis_for_model_t pp = rp.get_phi_psis_for_model(imod);
	 atom_sel.mol->GetModel(pp.model_number)->CalcSecStructure();
	 std::map<coot::residue_spec_t, coot::util::phi_psi_t>::const_iterator it;
	 for (it=pp.phi_psi.begin(); it!=pp.phi_psi.end(); it++) {
	    CResidue *residue_p = get_residue(it->first);
	    if (residue_p) { 
	       bool do_it = true; // unset for secondary structure
	       int sse = residue_p->SSE;
	       // std::cout << "residue->SSE is " << sse << " vs " << SSE_Strand << " and " << SSE_Helix
	       // << std::endl;
	       switch (residue_p->SSE)  {
	       case SSE_Strand:
		  do_it = 0;
		  break;
	       case SSE_Helix:
		  do_it = 0;
		  break;
	       }
	    
	       clipper::Ramachandran &rama = r_non_gly_pro;
	       if (it->second.residue_name() == "GLY")
		  rama = r_gly;
	       if (it->second.residue_name() == "PRO")
		  rama = r_pro;
	       ftype p = rama.probability(clipper::Util::d2rad(it->second.phi()),
					  clipper::Util::d2rad(it->second.psi()));
	       if (p < zero_cut_off) {
		  // std::cout << "........ was a zero" << std::endl;
		  rs.n_zeros++;
	       } else {
		  std::pair<coot::residue_spec_t, double> pair(it->first, p);
		  rs.scores.push_back(pair);
		  log_p_sum += log(p);
		  if (do_it) {
		     rs.scores_non_sec_str.push_back(pair);
		     log_p_non_sec_str_sum += log(p);
		  } 
	       }
	    }
	 }
      }
   }

   rs.score = log_p_sum;
   rs.score_non_sec_str = log_p_non_sec_str_sum;

   return rs;
}


coot::rotamer_score_t
molecule_class_info_t::get_all_molecule_rotamer_score(const coot::rotamer_probability_tables &rot_prob_tables) const {

   coot::rotamer_score_t rs;
   double sum_log_p = 0.0;

   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p) {
	       std::string residue_name = residue_p->GetResName();
	       if (coot::util::is_standard_residue_name(residue_name)) { 
		  try { 
		     std::vector<coot::rotamer_probability_info_t> v =
			rot_prob_tables.probability_this_rotamer(residue_p);
		     if (v.size() > 0) {
			if (v[0].state == coot::rotamer_probability_info_t::OK) { 
			   double p = v[0].probability;
			   if (p > 0) { 
			      // std::cout << "   " << coot::residue_spec_t(residue_p) << "  " << p << std::endl;
			      rs.add(coot::residue_spec_t(residue_p), p);
			      sum_log_p += log(p*0.01); // probability was in percents.
			   } else {
			      rs.n_pass++;
			   } 
			} else {
			   rs.n_pass++;
			} 
		     } else {
			rs.n_pass++;
		     } 
		  }
		  catch (std::runtime_error rte) {
		     std::cout << "Error:: " << rte.what() << std::endl;
		     rs.n_pass++;
		  }
	       }
	    } 
	 }
      }
   }

   rs.score = sum_log_p;
   
   return rs;
}

// --------- HETATMs ------------
int
molecule_class_info_t::residue_has_hetatms(const std::string &chain_id,
					   int resno,
					   const std::string &ins_code) const {

   // return 1 if any of the atoms are HETATMs.
   // 
   int r = -1;
   CResidue *residue_p = get_residue(chain_id, resno, ins_code);
   if (residue_p) {
      r = coot::util::residue_has_hetatms(residue_p);
   }
   return r;
}


int
molecule_class_info_t::hetify_residue_atoms(const std::string &chain_id,
					    int resno,
					    const std::string &ins_code) {

   int r = -1;
   CResidue *residue_p = get_residue(chain_id, resno, ins_code);
   if (residue_p) {
      make_backup();
      int n_atoms = coot::hetify_residue_atoms_as_needed(residue_p);
      if (n_atoms > 0)
	 r = 1;
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   }
   return r;
}

//////////////////////////////////////////////////////////////////////////////
////////////// animated ligands  /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void 
coot::animated_ligand_interactions_t::draw(CMMDBManager *mol,
					   const gl_context_info_t &gl_info,
					   const long &start_time) const {

   CAtom *at_1 = coot::util::get_atom(ligand_atom_spec, mol);
   CAtom *at_2 = coot::util::get_atom(interacting_residue_atom_spec, mol);

   if (at_1 && at_2) {
      
      clipper::Coord_orth pt_1 = coot::util::get_coords(at_1);
      clipper::Coord_orth pt_2 = coot::util::get_coords(at_2);
//       glBegin(GL_LINES);
//       glVertex3f(pt_1.x(), pt_1.y(),pt_1.z());
//       glVertex3f(pt_2.x(), pt_2.y(),pt_2.z());
//       glEnd();

      if (1) {
	 glEnable (GL_BLEND);
	 glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	 long now_time = glutGet(GLUT_ELAPSED_TIME);

	 double add_fac = 1.4;
	 double m_fac = 1/(1+add_fac) - 0.01;
	 
	 float opacity = (sin(double(now_time - start_time) * 0.002)+add_fac)*m_fac;
	 glColor4f(0.3, 0.3, 0.8, opacity);

	 GLfloat  ambientLight[] = { 0.2f, 0.2f, 0.2f, 0.5f };
	 GLfloat  diffuseLight[] = { 0.2f, 0.2f, 0.2f, 0.5f };
	 GLfloat specularLight[] = { 0.2f, 0.2f, 0.2f, 0.5f };
   
	 // Assign created components to GL_LIGHT2
	 // glLightfv(GL_LIGHT2, GL_AMBIENT, ambientLight);
	 // glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuseLight);
	 // glLightfv(GL_LIGHT2, GL_SPECULAR, specularLight);

	 glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);
	 glDisable(GL_COLOR_MATERIAL);

	 GLfloat  mat_specular[]  = {0.18, 0.18, 0.80,  opacity};
	 GLfloat  mat_ambient[]   = {0.20, 0.20, 0.20, opacity};
	 GLfloat  mat_diffuse[]   = {0.20, 0.20, 0.90, opacity};
	 GLfloat  mat_shininess[] = {100.0};

	 if (bond_type == coot::fle_ligand_bond_t::H_BOND_ACCEPTOR_SIDECHAIN ||
	     bond_type == coot::fle_ligand_bond_t::H_BOND_DONOR_SIDECHAIN) {

	    mat_specular[1] = 0.80; // g
	    mat_specular[2] = 0.10; // b
	    mat_diffuse[1]  = 0.80; // g 
	    mat_diffuse[2]  = 0.10; // b
	 } 

	 glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
	 glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	 glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_ambient);
	 glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);


	 double top =  0.07;
	 double base = 0.07;
	 int slices  = 12;
	 int stacks  = 2;
	 int n_parts = 25;
	 
	 clipper::Coord_orth bond_vec = pt_1-pt_2;
	 clipper::Coord_orth bond_frag(bond_vec * (1.0/double(n_parts)));
	 
	 for (unsigned int ipart=0; ipart<n_parts; ipart++) {

	    if (coot::util::even_p(ipart)) { 
	       glPushMatrix();
	       
	       double height = sqrt(bond_frag.lengthsq());

	       clipper::Coord_orth base_point = pt_2 + clipper::Coord_orth(bond_frag * double(ipart));

	       glTranslatef(base_point.x(), base_point.y(), base_point.z());


	       // 	    This code from ccp4mg's cprimitive.cc (but modified)
	       //  	    ----- 
	       double ax;
	       double rx = 0; 
	       double ry = 0;
	       double length = height;
	       double vz = bond_frag.z();
	 
	       bool rot_x = false;
	       if(fabs(vz)>1e-7){
		  ax = 180.0/M_PI*acos(vz/length);
		  if(vz<0.0) ax = -ax;
		  rx = -bond_frag.y()*vz;
		  ry = bond_frag.x()*vz;
	       }else{
		  double vx = bond_frag.x();
		  double vy = bond_frag.y();
		  ax = 180.0/M_PI*acos(vx/length);
		  if(vy<0) ax = -ax;
		  rot_x = true;
	       }
	 
	       if (rot_x) { 
		  glRotated(90.0, 0.0, 1.0, 0.0);
		  glRotated(ax,  -1.0, 0.0, 0.0);
	       } else {
		  glRotated(ax, rx, ry, 0.0);
	       }
	       // 	    --------

	       GLUquadric* quad = gluNewQuadric();
	       // glScalef(1.0, 1.0, -1.0); // account for mg maths :-)
	       gluCylinder(quad, base, top, height, slices, stacks);
	       gluDeleteQuadric(quad);
	       glPopMatrix();
	    }
	 }
      }
   } 
}

// add more variables :-)
void molecule_class_info_t::add_animated_ligand_interaction(const coot::fle_ligand_bond_t &lb) {

   animated_ligand_interactions_vec.push_back(lb);
}

void
molecule_class_info_t::draw_animated_ligand_interactions(const gl_context_info_t &gl_info,
							 const long &start_time) const {

   if (draw_animated_ligand_interactions_flag)
      for (unsigned int i=0; i<animated_ligand_interactions_vec.size(); i++) {
	 if (0) 
	    std::cout << "---- interaction " << i << " of "
		      << animated_ligand_interactions_vec.size() << std::endl;
	 animated_ligand_interactions_vec[i].draw(atom_sel.mol, gl_info, start_time);
      } 

} 
							    

std::pair<bool, clipper::Coord_orth>
molecule_class_info_t::residue_centre(const std::string &chain_id, int resno, const std::string &ins_code) const {

   CResidue *residue_p = get_residue(chain_id, resno, ins_code);
   if (residue_p) {
      return residue_centre(residue_p);
   } else {
      std::cout << "Residue not found "
		<< coot::residue_spec_t(chain_id, resno, ins_code)
		<< std::endl;
      bool r = 0;
      clipper::Coord_orth pos(0,0,0);
      return std::pair<bool, clipper::Coord_orth> (r, pos);
   }
}


std::pair<bool, clipper::Coord_orth>
molecule_class_info_t::residue_centre(CResidue *residue_p) const {

   bool r = 0;
   clipper::Coord_orth pos(0,0,0);
   int n_atoms = 0;

   if (residue_p) {
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 if (! residue_atoms[iat]->isTer()) { 
	    clipper::Coord_orth p(residue_atoms[iat]->x,
				  residue_atoms[iat]->y,
				  residue_atoms[iat]->z);
	    pos += p;
	    n_atoms++;
	 }
      }
      if (n_atoms) {
	 pos = clipper::Coord_orth(pos.x()/double(n_atoms),
				   pos.y()/double(n_atoms),
				   pos.z()/double(n_atoms));
	 r = 1; 
      } 
   } 
   return std::pair<bool, clipper::Coord_orth> (r, pos);
}



// ------------------- ligand centre ---------------------
// we want a button that goes to the ligand when we click it.
// (if we are already on a ligand, go to the next one).
// 
// the first value flags if this contains a useful return value.
// 1: normal case, go ahead and use the coord
// 0:  No ligands found
// -1: No movement because we are at the (single) ligand already
//
coot::new_centre_info_t
molecule_class_info_t::new_ligand_centre(const clipper::Coord_orth &current_centre,
					 int n_atoms_min) const {

   clipper::Coord_orth pos(0,0,0);
   coot::new_ligand_position_type status = coot::NO_LIGANDS;
   coot::residue_spec_t residue_spec;

   std::vector<std::pair<clipper::Coord_orth, coot::residue_spec_t> > ligand_centres;
   
   if (atom_sel.n_selected_atoms > 0) { 
      int imod = 1;
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    if (n_atoms >= n_atoms_min) { 
	       bool is_het = false;
	       for (int iat=0; iat<n_atoms; iat++) {
		  at = residue_p->GetAtom(iat);
		  if (at->Het) {
		     std::string res_name = residue_p->GetResName();
		     if (res_name != "HOH")
			if (res_name != "WAT")
			   is_het = true;
		     break;
		  } 
	       }
	       if (is_het) {
		  std::pair<bool, clipper::Coord_orth> res_centre = residue_centre(residue_p);
		  if (res_centre.first) {
		     std::pair<clipper::Coord_orth, coot::residue_spec_t> p(res_centre.second, residue_p);
		     ligand_centres.push_back(p);
		  }
	       }
	    } 
	 }
      }
   }

   int current_centre_index = -1; // unset
   if (ligand_centres.size()) {
      for (unsigned int ilig=0; ilig<ligand_centres.size(); ilig++) { 
	 double d = clipper::Coord_orth::length(current_centre, ligand_centres[ilig].first);
	 if (d < 5) {
	    current_centre_index = ilig;
	    break;
	 } 
      }

      if (current_centre_index == -1) {
	 // go to the first
	 status = coot::NORMAL_CASE;
	 pos = ligand_centres[0].first;
	 residue_spec = ligand_centres[0].second;
      } else {
	 if (ligand_centres.size() > 1) {
	    int next_ligand_centre_index = 0;
	    // we are on the last one
	    if (current_centre_index == (ligand_centres.size()-1)) {
	       next_ligand_centre_index = 0;
	    } else {
	       // normal case, go to the next one
	       next_ligand_centre_index = current_centre_index + 1;
	    }
	    status = coot::NORMAL_CASE;
	    pos =          ligand_centres[next_ligand_centre_index].first;
	    residue_spec = ligand_centres[next_ligand_centre_index].second;
	 } else {
	    // there was only one ligand and we are centred on it.
	    // Do nothing.
	    status = coot::SINGLE_LIGAND_NO_MOVEMENT;
	    residue_spec = ligand_centres[0].second;
	 }
      }
   }

   coot::new_centre_info_t nci(status, pos, residue_spec);
   return nci;
} 


int
molecule_class_info_t::watson_crick_pair_for_residue_range(const std::string &chain_id,
							   int resno_start, int resno_end,
							   CMMDBManager *standard_residues_mol) {
   int status = 0;
   std::vector<CResidue *> new_residues;
   CModel *model_p = NULL;
   if (resno_end < resno_start)
      std::swap(resno_start, resno_end);

   for (int ires=resno_start; ires<=resno_end; ires++) { 
      CResidue *res = get_residue(chain_id, ires, "");
      if (!res) {
	 std::cout << "Residue not found in  " << chain_id << " " << ires
		   << std::endl;
      } else {
	 model_p = res->GetModel(); // set it to the model that contains this chain.
	 // these residues are (simple) deep copied
	 CResidue *res_wc =
	    coot::watson_crick_partner(res, standard_residues_mol);
	 if (res_wc) {
	    new_residues.push_back(res_wc);
	 }
      } 
   }

   if (new_residues.size()) {
      make_backup();
      CChain *chain_p = new CChain;
      // set the chain id
      std::pair<short int, std::string> u = unused_chain_id();
      if (u.first) { 
	 chain_p->SetChainID(u.second.c_str());
	 for (unsigned int ires=0; ires<new_residues.size(); ires++) { 
	    chain_p->AddResidue(new_residues[ires]);
	    new_residues[ires]->seqNum = new_residues.size() - ires;
	    status = 1;
	 }
	 model_p->AddChain(chain_p);
	 
	 atom_sel.mol->FinishStructEdit();
	 update_molecule_after_additions();
      }
   } 
   return status;
} 

// Try to add a LINK record if LINK is not blank.
//
// Passed residue new_res does not have a useful residue_number.
// 
std::pair<bool, CResidue *>
molecule_class_info_t::add_residue(CResidue *new_res,
				   const std::string &chain_id_in) {

   bool status = false;
   CResidue *res_copied = NULL;
   int imod = 1;
   if (new_res) { 
      if (atom_sel.n_selected_atoms > 0) { 
	 CModel *model_p = atom_sel.mol->GetModel(imod);
	 CChain *chain_p;
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    std::string chain_id(chain_p->GetChainID());
	    if (chain_id == chain_id_in) {
	       res_copied = copy_and_add_residue_to_chain(chain_p, new_res);
	       status = true;
	       have_unsaved_changes_flag = 1;
	       atom_sel.mol->FinishStructEdit();
	       update_molecule_after_additions();
	       break;
	    }
	 }
      }
   }
   return std::pair<bool, CResidue *> (status, res_copied);
}


// Add a LINK record if link_type is not blank (link_type is for example "NAG-ASN")
// 
coot::residue_spec_t
molecule_class_info_t::add_linked_residue(const coot::residue_spec_t &spec_in,
					  const std::string &new_residue_comp_id,
					  const std::string &link_type,
					  coot::protein_geometry *geom_p) {

   coot::residue_spec_t new_residue_spec;   
   CResidue *residue_ref = get_residue(spec_in);
   if (residue_ref) {
      coot::beam_in_linked_residue lr(residue_ref, link_type, new_residue_comp_id, geom_p);
      CResidue *result = lr.get_residue();
      // get_residue() can (and often does) modify residue_ref (for
      // deleting link mod atom, for example). So we need a FinishStructEdit() here
      atom_sel.mol->FinishStructEdit();
      
      std::pair<bool, CResidue *> status_pair = add_residue(result, spec_in.chain);

      if (status_pair.first) {

	 new_residue_spec = coot::residue_spec_t(status_pair.second);
	 
	 try { 
	    coot::dict_link_info_t link_info(residue_ref, status_pair.second,
					     link_type, *geom_p);
	    make_link(link_info.spec_ref, link_info.spec_new, link_type, link_info.dist);
	 }
	 catch (std::runtime_error rte) {
	    // we didn't find the info to make the link specs, oh well...
	    std::cout << "WARNING:: add_linked_residue() catches exception \""
		      << rte.what() << "\"" << std::endl;
	 } 
      }
   }
   return new_residue_spec;
} 

// this can throw a std::runtime_error
coot::dict_link_info_t::dict_link_info_t(CResidue *residue_ref,
					 CResidue *residue_new,
					 const std::string &link_type,
					 const coot::protein_geometry &geom) {

   if (0) { 
      std::cout << "dict_link_info_t() constructor start with link_type: "
		<< link_type << std::endl;
      std::cout << "dict_link_info_t() constructor start with residue_ref: "
		<< coot::residue_spec_t(residue_ref) << std::endl;
      std::cout << "dict_link_info_t() constructor start with residue_new: "
		<< coot::residue_spec_t(residue_new) << std::endl;
   }
   
   bool ifound = false;
   if (! residue_ref) {
      throw (std::runtime_error("Null residue_ref"));
   } else { 
      if (! residue_ref) {
	 throw (std::runtime_error("Null residue_new"));
      } else {
	 coot::dictionary_residue_link_restraints_t rr = geom.link(link_type);
	 if (rr.link_id == "") {
	    throw (std::runtime_error("Link not found in dictionary"));
	 } else {

	    bool order_switch_flag = check_for_order_switch(residue_ref,
							    residue_new,
							    link_type, geom);

	    CResidue *res_1 = residue_ref;
	    CResidue *res_2 = residue_new;

	    if (order_switch_flag) { 
	       std::swap(res_1, res_2);
	    }
	    
	    // we found it (i.e. not null object)
	    coot::residue_spec_t res_spec_ref(res_1);
	    coot::residue_spec_t res_spec_new(res_2);
	    for (unsigned int ibond=0; ibond<rr.link_bond_restraint.size(); ibond++) {
	       PPCAtom residue_atoms_1 = 0;
	       int n_residue_atoms_1;
	       res_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
	       for (unsigned int iat1=0; iat1<n_residue_atoms_1; iat1++) { 
		  std::string atom_name_1(residue_atoms_1[iat1]->name);
		  if (atom_name_1 == rr.link_bond_restraint[ibond].atom_id_1_4c()) {
		     // OK so the first atom matched
		     PPCAtom residue_atoms_2 = 0;
		     int n_residue_atoms_2;
		     res_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
		     for (unsigned int iat2=0; iat2<n_residue_atoms_2; iat2++) { 
			std::string atom_name_2(residue_atoms_2[iat2]->name);
			if (atom_name_2 == rr.link_bond_restraint[ibond].atom_id_2_4c()) {
			   ifound = 1;
			   spec_ref = coot::atom_spec_t(res_spec_ref.chain,
							res_spec_ref.resno,
							res_spec_ref.insertion_code,
							atom_name_1, "");
			   spec_new = coot::atom_spec_t(res_spec_new.chain,
							res_spec_new.resno,
							res_spec_new.insertion_code,
							atom_name_2, "");
			   dist = coot::distance(residue_atoms_1[iat1],
						 residue_atoms_2[iat2]);
			   break;
			}
			if (ifound)
			   break;
		     }
		  }
		  if (ifound)
		     break;
	       }
	       if (ifound)
		  break;
	    }

	    if (! ifound)
	       throw std::runtime_error("Dictionary links atom not found in link residues");
	 }
      }
   }
}

bool
coot::dict_link_info_t::check_for_order_switch(CResidue *residue_ref,
					       CResidue *residue_new,
					       const std::string &link_type,
					       const coot::protein_geometry &geom) const {

   bool order_switch_flag = false;
   std::string comp_id_ref = residue_ref->GetResName();
   std::string comp_id_new = residue_new->GetResName();

   try {
      std::string group_ref = geom.get_group(residue_ref);
      // std::cout << "got group_ref: " << group_ref << std::endl;
      std::string group_new = geom.get_group(residue_new);
      // std::cout << "got group_new:o " << group_new << std::endl;
      std::vector<std::pair<coot::chem_link, bool> > link_infos =
	 geom.matching_chem_link(comp_id_ref, group_ref, comp_id_new, group_new);
      if (0)
	 std::cout << "DEBUG:: in check_for_order_switch() found " << link_infos.size()
		   << " link infos " << std::endl;
      for (unsigned int ilink=0; ilink<link_infos.size(); ilink++) {
	 // std::cout << "   chem_link: " << ilink << " " << link_infos[ilink].first
	 // << " " << link_infos[ilink].second << std::endl;
	 if (link_infos[ilink].first.Id() == link_type) { 
	    order_switch_flag = link_infos[ilink].second;
	    break;
	 }
      }
   }
   catch (std::runtime_error rte) {
      std::cout << "WARNING:: check_for_order_switch() exception: " << rte.what() << std::endl;
   } 
   return order_switch_flag;
}

// multi-residue torsion map fitting interface
void
molecule_class_info_t::multi_residue_torsion_fit(const std::vector<coot::residue_spec_t> &residue_specs,
						 const clipper::Xmap<float> &xmap,
						 coot::protein_geometry *geom_p) {

   CMMDBManager *moving_mol =
      coot::util::create_mmdbmanager_from_residue_specs(residue_specs, atom_sel.mol);
   
   // do we need to send over the base atom too?  Or just say
   // that it's the first atom in moving_mol?
   // 
   coot::multi_residue_torsion_fit_map(moving_mol, xmap, geom_p);

   atom_selection_container_t moving_atoms_asc = make_asc(moving_mol);
   replace_coords(moving_atoms_asc, 1, 1);
   
}


