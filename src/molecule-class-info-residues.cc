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
#include <stdexcept>


// #ifdef MAKE_ENTERPRISE_TOOLS
// // includes order important, otherwise we get dcgettext() problems.
// #include "rdkit-interface.hh" // needed for add_hydrogens()
// #endif

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "Cartesian.h"
#include "mmdb-crystal.h"
// 
#include "rotamer.hh" // in ligand

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


// only apply charges if the molecule contains lots of hydrogens.
void
molecule_class_info_t::apply_charges(const coot::protein_geometry &geom) {

   // More than 20% of the atoms have to be hydrogens for use to set
   // the charges on all the atoms (to something other than
   // CXX_UNSET_CHARGE).
   
   float fraction_hydrogens = 0.2;

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
      
      if ( (float(n_H)/float(n_all) > fraction_hydrogens) || n_all < 40) {

	 // first set all atom charges to unset:
	 for (int i=0; i<atom_sel.n_selected_atoms; i++)
	    atom_sel.atom_selection[i]->charge = CXX_UNSET_CHARGE - 0.1;
      

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


bool
molecule_class_info_t::sprout_hydrogens(const std::string &chain_id,
					int res_no,
					const std::string &ins_code,
					const coot::protein_geometry &geom) {

   bool r = 0;

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
	 std::cout << "DEBUG:: coot::add_hydrogens() returns " << r << std::endl;
	 if (r) {


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
		     coot::is_bad_chiral_atom_p(cr[icr], residue_p);
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
   if (r) {
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
	    if (! c_residue_atoms[jat]->isTer()) { 
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
