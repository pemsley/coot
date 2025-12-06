/* src/molecule-class-info-residues.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009 by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

// Must include these headers to get molecule_class_info_t.h to parse.
//

#include "compat/coot-sysdep.h"

#include <string.h>

// #ifdef MAKE_ENHANCED_LIGAND_TOOLS
// // includes order important, otherwise we get dcgettext() problems.
// #include "rdkit-interface.hh" // needed for add_hydrogens()
// #endif

#include <mmdb2/mmdb_manager.h>
#include <mmdb2/mmdb_tables.h>  // for mmdb::Get1LetterCode()

#include "ligand/rotamer.hh" // in ligand

#include "coot-utils/glyco-tree.hh"

#include "ligand/base-pairing.hh"
#include "ideal/torsion-bonds.hh"

#include "geometry/mol-utils.hh"
#include "utils/coot-utils.hh"
#include "ideal/add-linked-cho.hh"

#include "molecule-class-info.h"
#include "coot-hydrogens.hh"
#include "graphics-info.h"


// 1: success
// 0: failure
//
bool
molecule_class_info_t::progressive_residues_in_chain_check_by_chain(const char *chain_id) const {

   bool r = 0;

   if (atom_sel.n_selected_atoms > 0) {
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 mmdb::Chain *chain_p = atom_sel.mol->GetChain(1,ichain);
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
      mmdb::Manager *mol = atom_sel.mol;

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
	 mmdb::Model *model_p = mol->GetModel(imod);
	 if (model_p) {
	    mmdb::Chain *chain_p;
	    int nchains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::PResidue residue_p;
	       for (int ires=0; ires<nres; ires++) {
		  residue_p = chain_p->GetResidue(ires);
		  std::string res_type = residue_p->GetResName();
		  std::pair<short int, coot::dictionary_residue_restraints_t> rp =
		     geom.get_monomer_restraints(res_type, imol_no);
		  if (rp.first) {
		     try {
			coot::dipole p(rp.second, residue_p);
			p.fill_charged_atoms(residue_p, rp.second);
		     }
		     catch (const std::runtime_error &mess) {
			std::cout << mess.what() << std::endl;
		     }
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
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (! model_p) continue;
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::PResidue residue_p;
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
   mmdb::Residue *residue_p = get_residue(chain_id, res_no, ins_code);
   std::vector<coot::atom_spec_t> fixed_atoms;
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++)
      if (std::string(residue_atoms[i]->element) != " H")
	 fixed_atoms.push_back(coot::atom_spec_t(residue_atoms[i]));

   if (residue_p) {
      std::string residue_type = residue_p->GetResName();
      std::pair<bool, coot::dictionary_residue_restraints_t> p =
	 geom.get_monomer_restraints_at_least_minimal(residue_type, imol_no);
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
	       geom.get_monomer_restraints(residue_name, imol_no);

	    if (rp.first) {

	       std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs(residue_p);
	       for (unsigned int iac=0; iac<alt_confs.size(); iac++) {
		  const std::string &alt_conf = alt_confs[iac];


		  mmdb::Manager *residue_mol = coot::util::create_mmdbmanager_from_residue(residue_p);
		  mmdb::Residue *residue_cp_p = coot::util::get_first_residue(residue_mol);

		  coot::util::delete_alt_confs_except(residue_cp_p, alt_conf);
		  residue_mol->FinishStructEdit();

		  if (false) { // -------- debug ----------
		     int imod = 1;
		     mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
		     mmdb::Chain *chain_p;
		     int n_chains = model_p->GetNumberOfChains();
		     for (int ichain=0; ichain<n_chains; ichain++) {
			chain_p = model_p->GetChain(ichain);
			int nres = chain_p->GetNumberOfResidues();
			mmdb::Residue *residue_inner_p;
			for (int ires=0; ires<nres; ires++) {
			   residue_inner_p = chain_p->GetResidue(ires);
			   std::cout << "copied mol residue " << coot::residue_spec_t(residue_inner_p)
				     << " " << residue_inner_p << " vs " << residue_cp_p
				     << std::endl;
			}
		     }
		  } // ------ end debug ----------

		  // those Hs were just attached with non-good geometry, we
		  // need to minimise.  Keep all atoms fixed except all hydrogens.
		  std::vector<std::pair<bool,mmdb::Residue *> > residues;
		  std::pair<bool, mmdb::Residue *> pp(0, residue_cp_p);
		  residues.push_back(pp);
		  clipper::Xmap<float> dummy_xmap;

		  coot::restraints_container_t restraints(residues,
							  atom_sel.links,
							  geom, residue_mol, fixed_atoms,
							  &dummy_xmap);
                  int n_threads = coot::get_max_number_of_threads() - 1;
                  if (n_threads < 1) n_threads = 1;
                  restraints.thread_pool(&graphics_info_t::static_thread_pool, n_threads);
		  bool do_torsions = 0;

		  coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
		  bool do_trans_peptide_restraints = false;
		  int n_restraints = restraints.make_restraints(imol_no, geom, flags, do_torsions,
								do_trans_peptide_restraints,
								0, 0, false, false, false, coot::NO_PSEUDO_BONDS);
		  restraints.minimize(flags);
		  residue_mol->FinishStructEdit();

		  // do we need to optimize some more because the atom names have been swapped?
		  bool needs_more = sprout_hydrogens_correct_chirals_maybe(residue_cp_p, alt_conf, rp.second);

		  if (needs_more) {
		     coot::restraints_container_t restraints_2(residues,
							       atom_sel.links,
							       geom, residue_mol, fixed_atoms, &dummy_xmap);
		     flags = coot::CHIRAL_VOLUMES;
		     flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
		     n_restraints = restraints_2.make_restraints(imol_no, geom,
								 flags, do_torsions,
								 do_trans_peptide_restraints,
								 0, 0, false, false, false, coot::NO_PSEUDO_BONDS);
		     restraints_2.minimize(flags);
		  }

		  // from_res, to_res
		  sprout_hydrogens_transfer_hydrogen_positions(residue_cp_p, residue_p, alt_conf);
	       }
	    }
	 }
      }
   }
   if (r.first) {
      have_unsaved_changes_flag = 1;
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked(__FUNCTION__);
   }
   return r;
}

// if the original residue has alt confs, residue_cp_p should have
// only one alt conf when this is called (e.g. only "A"s).
//
bool
molecule_class_info_t::sprout_hydrogens_correct_chirals_maybe(mmdb::Residue *residue_cp_p,
							      const std::string &alt_conf,
							      const coot::dictionary_residue_restraints_t &rp) {

   bool correct_chiral = false;

   // analyse results
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_cp_p->GetAtomTable(residue_atoms, n_residue_atoms);

   //
   std::vector<coot::dict_chiral_restraint_t> cr = rp.chiral_restraint;

   for (unsigned int icr=0; icr<cr.size(); icr++) {
      std::string centre_atom = cr[icr].atom_id_c_4c();
      std::vector<std::pair<short int, coot::atom_spec_t> > v =
	 coot::is_inverted_chiral_atom_p(cr[icr], residue_cp_p);
      if (v.size() ) {
	 for (unsigned int i=0; i<v.size(); i++) {
	    if (v[i].first) {
	       std::cout << "fix this bad chiral centre "
			 << v[i].first << " "
			 << v[i].second << std::endl;
	       std::vector<std::string> attached_Hs =
		  rp.get_attached_H_names(v[i].second.atom_name);
	       if (attached_Hs.size() > 1) {

		  coot::atom_spec_t spec_1 = v[i].second;
		  coot::atom_spec_t spec_2 = v[i].second;
		  spec_1.atom_name = attached_Hs[0];
		  spec_2.atom_name = attached_Hs[1];
		  mmdb::Atom *at_1 = get_atom(spec_1);
		  mmdb::Atom *at_2 = get_atom(spec_2);

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
			correct_chiral = true;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return correct_chiral;
}

void
molecule_class_info_t::sprout_hydrogens_transfer_hydrogen_positions(mmdb::Residue *from_res_p, mmdb::Residue *to_res_p,
								    const std::string &alt_conf) {

   mmdb::PPAtom residue_atoms_1 = 0;
   mmdb::PPAtom residue_atoms_2 = 0;
   int n_residue_atoms_1;
   int n_residue_atoms_2;
   from_res_p->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   to_res_p->GetAtomTable  (residue_atoms_2, n_residue_atoms_2);

   for (int iat_1=0; iat_1<n_residue_atoms_1; iat_1++) {
      mmdb::Atom *at_1 = residue_atoms_1[iat_1];
      std::string ele_1 = at_1->element;
      if (ele_1 == " H") {                // PDBv3 FIXME
	 std::string name_1(at_1->name);
	 std::string altc_1(at_1->altLoc);
	 if (altc_1 == alt_conf) {
	    for (int iat_2=0; iat_2<n_residue_atoms_2; iat_2++) {
	       mmdb::Atom *at_2 = residue_atoms_2[iat_2];
	       std::string ele_2 = at_2->element;
	       if (ele_2 == " H") {       // PDBv3 FIXME
		  std::string name_2(at_2->name);
		  std::string altc_2(at_2->altLoc);

		  if (name_1 == name_2) {
		     if (altc_1 == altc_2) {
			at_2->x = at_1->x;
			at_2->y = at_1->y;
			at_2->z = at_1->z;
		     }
		  }
	       }
	    }
	 }
      }
   }
}



std::vector<std::string>
molecule_class_info_t::no_dictionary_for_residue_type_as_yet(const coot::protein_geometry &geom) const {

   std::vector<std::string> v;
   if (has_model()) {
      int imod = 1;
      if (atom_sel.mol) {
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 if (model_p) {
	    mmdb::Chain *chain_p;
	    int nchains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::Residue *residue_p;
	       for (int ires=0; ires<nres; ires++) {
		  residue_p = chain_p->GetResidue(ires);
		  std::string residue_name = residue_p->GetResName();
		  if (! geom.have_at_least_minimal_dictionary_for_residue_type(residue_name, imol_no)) {

// 		  // now check the HETSYNs
// 		  for (unsigned int ihet=0; ihet<model_p->HetCompounds.nHets; ihet++) {
// 		     std::string het_id = model_p->HetCompounds[ihet].hetID;
// 		     std::cout << "het_id: " << het_id << std::endl;
// 		     if (het_id == residue_name) {
// 			int ns = model_p->HetCompounds[ihet].nSynonyms;
// 			for (unsigned int i_syn=0; i_syn<ns; i_syn++) {

// 			}
// 		     }
// 		  }


		     // add it to v, if it is not already there:
		     if (std::find(v.begin(), v.end(), residue_name) == v.end())
			v.push_back(residue_name);
		  }
	       }
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

   mmdb::Residue *central_residue = get_residue(central_residue_spec);
   mmdb::Residue *neighbour_residue = get_residue(neighbour_residue_spec);

   if (! central_residue) {
      std::string message = "Missing residue ";
      message += central_residue_spec.chain_id;
      message += " ";
      message += central_residue_spec.res_no;
      throw std::runtime_error(message);
   } else {
      if (! neighbour_residue) {
	 std::string message = "Missing residue ";
	 message += neighbour_residue_spec.chain_id;
	 message += " ";
	 message += neighbour_residue_spec.res_no;
	 throw std::runtime_error(message);
      } else {
	 // OK! go.

	 double min_dist = 42e32;
	 clipper::Coord_orth shortest_dist(0,0,0); // "best"
	 mmdb::PPAtom c_residue_atoms = 0;
	 int c_n_residue_atoms;
	 mmdb::PPAtom n_residue_atoms = 0;
	 int n_n_residue_atoms;
	 central_residue->GetAtomTable  (c_residue_atoms, c_n_residue_atoms);
	 neighbour_residue->GetAtomTable(n_residue_atoms, n_n_residue_atoms);
	 clipper::Coord_orth sum_1(0,0,0);
	 clipper::Coord_orth sum_2(0,0,0);
	 for (int iat=0; iat<c_n_residue_atoms; iat++) {
	    if (! c_residue_atoms[iat]->isTer()) {
	       clipper::Coord_orth pt_1(c_residue_atoms[iat]->x,
					c_residue_atoms[iat]->y,
					c_residue_atoms[iat]->z);
	       sum_1 += pt_1;
	    }
	 }
	 for (int jat=0; jat<n_n_residue_atoms; jat++) {
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

#include "coot-phi-psi.hh" // 20230611-PE new header

void
molecule_class_info_t::match_ligand_atom_names(const std::string &chain_id, int res_no, const std::string &ins_code,
					       mmdb::Residue *res_ref) {

   mmdb::Residue *res_mov = get_residue(chain_id, res_no, ins_code);

   if (! res_mov) {
      std::cout << "No residue for moving atom names:  " << chain_id << " " << res_no << " "  << ins_code
		<< std::endl;

   } else {
      bool match_hydrogens_also = 1;
      coot::graph_match_info_t gm = coot::graph_match(res_mov, res_ref, 0, match_hydrogens_also);
      gm.match_names(res_mov);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
}

coot::rama_score_t
molecule_class_info_t::get_all_molecule_rama_score() const {

   G_GNUC_UNUSED auto debug_the_probabilities = [] () {
      clipper::Ramachandran rama;
      ftype level_prefered = graphics_info_t::rama_level_prefered;
      ftype level_allowed  = graphics_info_t::rama_level_allowed;
      rama.init(clipper::Ramachandran::All2); // needs adjustment
      rama.set_thresholds(level_prefered, level_allowed);

      std::ofstream f_out("debug-rama.table");

      for (int i_phi= -180; i_phi<=180; i_phi++) {
         for (int i_psi= -180; i_psi<=180; i_psi++) {
            double phi = i_phi;
            double psi = i_psi;
            ftype p = rama.probability(clipper::Util::d2rad(phi), clipper::Util::d2rad(psi));
            f_out << "phi " << phi << " psi " << psi << " pr " << p << "\n";
         }
      }
      f_out.close();
   };
   
   // debug_the_probabilities();

   coot::rama_score_t rs;

   ftype level_prefered = graphics_info_t::rama_level_prefered;
   ftype level_allowed  = graphics_info_t::rama_level_allowed;

   // 20230608-PE Clipper does have TOP8000 by now - I am not going to test for it.
   //
   clipper::Ramachandran r_gly, r_pro, r_non_gly_pro, r_all;
   clipper::Ramachandran r_ileval, r_pre_pro, r_non_gly_pro_pre_pro_ileval;
   r_gly.init(clipper::Ramachandran::Gly2);
   r_gly.set_thresholds(level_prefered, level_allowed);
   //
   r_pro.init(clipper::Ramachandran::Pro2);
   r_pro.set_thresholds(level_prefered, level_allowed);
   //
   // as usual we make a first approximation and add some new ones...
   r_non_gly_pro.init(clipper::Ramachandran::NoGPIVpreP2);
   r_non_gly_pro.set_thresholds(level_prefered, level_allowed);
   //
   r_ileval.init(clipper::Ramachandran::IleVal2);
   r_ileval.set_thresholds(level_prefered, level_allowed);
   //
   r_pre_pro.init(clipper::Ramachandran::PrePro2);
   r_pre_pro.set_thresholds(level_prefered, level_allowed);
   //
   r_non_gly_pro_pre_pro_ileval.init(clipper::Ramachandran::NoGPIVpreP2);
   r_non_gly_pro_pre_pro_ileval.set_thresholds(level_prefered, level_allowed);
   //
   r_all.init(clipper::Ramachandran::All2);
   r_all.set_thresholds(level_prefered, level_allowed);

   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) {

      coot::phi_psis_for_model_t pp(atom_sel.mol, imod);
      std::map<coot::residue_spec_t, coot::util::phi_psi_with_residues_t>::const_iterator it;
      for (it=pp.phi_psi.begin(); it!=pp.phi_psi.end(); ++it) {
         const auto &residue_spec = it->first;
         const auto &ppr = it->second;
         mmdb::Residue *residue_p = get_residue(residue_spec);
         if (residue_p) {
            // there seems to be a problem here of using the correct rama...
            clipper::Ramachandran *rama_p = &r_non_gly_pro;
            if (ppr.is_pre_pro()) rama_p = &r_pre_pro;
            if (ppr.residue_name() == "GLY") rama_p = &r_gly;
            if (ppr.residue_name() == "PRO") rama_p = &r_pro;
            if (ppr.residue_name() == "ILE") rama_p = &r_ileval;
            if (ppr.residue_name() == "VAL") rama_p = &r_ileval;

            double phi = ppr.phi();
            double psi = ppr.psi();
            ftype p = rama_p->probability(clipper::Util::d2rad(phi), clipper::Util::d2rad(psi));
            coot::rama_score_t::scored_phi_psi_t scored_phi_psi(residue_spec, p, ppr);
            scored_phi_psi.set_residues(it->second);
            rs.scores.push_back(scored_phi_psi);
            if (false)
               std::cout << "debug:: get_all_molecule_rama_score() added " << residue_spec
                         << " phi " << phi << " psi " << psi << " "
                         << ppr.residue_name() << " " << ppr.is_pre_pro() << " pr: " << p << std::endl;
         }
      }
   }

   return rs;
}


coot::rama_score_t
molecule_class_info_t::get_all_molecule_rama_score_old() const {

   coot::rama_score_t rs;

#ifdef HAVE_GOOCANVAS
   coot::rama_plot rp;
   rp.generate_phi_psis(atom_sel.mol, true);

   int n = rp.n_phi_psi_model_sets();

   clipper::Ramachandran r_gly, r_pro, r_non_gly_pro, r_all;

   // Lovell et al. 2003, 50, 437 Protein Structure, Function and Genetics values:
   // BL says:: should these be the one set?!
   //ftype level_prefered = 0.02;
   //ftype level_allowed = 0.002;

   ftype level_prefered = graphics_info_t::rama_level_prefered;
   ftype level_allowed = graphics_info_t::rama_level_allowed;

   // clipper defaults: 0.01 0.0005

#ifdef CLIPPER_HAS_TOP8000
   clipper::Ramachandran r_ileval, r_pre_pro, r_non_gly_pro_pre_pro_ileval;
   r_gly.init(clipper::Ramachandran::Gly2);
   r_gly.set_thresholds(level_prefered, level_allowed);
   //
   r_pro.init(clipper::Ramachandran::Pro2);
   r_pro.set_thresholds(level_prefered, level_allowed);
   //
   // as usual we make a first approximation and add some new ones...
   r_non_gly_pro.init(clipper::Ramachandran::NoGPIVpreP2);
   r_non_gly_pro.set_thresholds(level_prefered, level_allowed);
   //
   r_ileval.init(clipper::Ramachandran::IleVal2);
   r_ileval.set_thresholds(level_prefered, level_allowed);
   //
   r_pre_pro.init(clipper::Ramachandran::PrePro2);
   r_pre_pro.set_thresholds(level_prefered, level_allowed);
   //
   r_non_gly_pro_pre_pro_ileval.init(clipper::Ramachandran::NoGPIVpreP2);
   r_non_gly_pro_pre_pro_ileval.set_thresholds(level_prefered, level_allowed);
   //
   r_all.init(clipper::Ramachandran::All2);
   r_all.set_thresholds(level_prefered, level_allowed);
#else
   r_gly.init(clipper::Ramachandran::Gly5);
   r_gly.set_thresholds(level_prefered, level_allowed);
   //
   r_pro.init(clipper::Ramachandran::Pro5);
   r_pro.set_thresholds(level_prefered, level_allowed);
   //
   r_non_gly_pro.init(clipper::Ramachandran::NonGlyPro5);
   r_non_gly_pro.set_thresholds(level_prefered, level_allowed);

   r_all.init(clipper::Ramachandran::All5);
   r_all.set_thresholds(level_prefered, level_allowed);
#endif

   double zero_cut_off = 1e-6;

   double log_p_sum = 0.0;
   double log_p_non_sec_str_sum = 0.0;
   int region;

   for (int imod=1; imod<n; imod++) {
      if (imod<=atom_sel.mol->GetNumberOfModels()) {

         coot::phi_psis_for_model_t pp = rp.get_phi_psis_for_model(imod);
         atom_sel.mol->GetModel(pp.model_number)->CalcSecStructure();
         std::map<coot::residue_spec_t, coot::util::phi_psi_with_residues_t>::const_iterator it;
         for (it=pp.phi_psi.begin(); it!=pp.phi_psi.end(); it++) {
            mmdb::Residue *residue_p = get_residue(it->first);
            if (residue_p) {
               bool do_it = true; // unset for secondary structure
               int sse = residue_p->SSE;

               switch (residue_p->SSE)  {
               case mmdb::SSE_Strand:
                  do_it = false;
                  break;
               case mmdb::SSE_Helix:
                  do_it = false;
                  break;
               }

               clipper::Ramachandran rama = r_non_gly_pro;
               if (it->second.residue_name() == "GLY")
                   rama = r_gly;
               if (it->second.residue_name() == "PRO")
                 rama = r_pro;
               double phi;
               double psi;
               phi = it->second.phi();
               psi = it->second.psi();
               ftype p = rama.probability(clipper::Util::d2rad(phi),
                                          clipper::Util::d2rad(psi));
               std::pair<coot::residue_spec_t, int> rama_pair;
               if (rama.allowed(clipper::Util::d2rad(phi),
                                clipper::Util::d2rad(psi))) {
                  region = coot::rama_plot::RAMA_ALLOWED;
                  if (rama.favored(clipper::Util::d2rad(phi),
                                   clipper::Util::d2rad(psi))) {
                     region = coot::rama_plot::RAMA_PREFERRED;
                  }
               } else {
                  region = coot::rama_plot::RAMA_OUTLIER;
               }
               rama_pair = std::make_pair(it->first, region);
               rs.region.push_back(rama_pair);

               if (p < zero_cut_off) {
                  // std::cout << "........ was a zero" << std::endl;
                  rs.n_zeros++;
               } else {
		  // old
                  // std::pair<coot::residue_spec_t, double> pair(it->first, p);
                  // rs.scores.push_back(pair);

		  coot::rama_score_t::scored_phi_psi_t scored_phi_psi(it->first, p, it->second);
		  scored_phi_psi.set_residues(it->second);

                  rs.scores.push_back(scored_phi_psi);
                  log_p_sum += log(p);
                  if (do_it) {
                     rs.scores_non_sec_str.push_back(scored_phi_psi);
                     log_p_non_sec_str_sum += log(p);
                  }
               }
            }
         }
      }
   }

   rs.score = log_p_sum;
   rs.score_non_sec_str = log_p_non_sec_str_sum;
#endif
   return rs;
}


coot::rotamer_score_t
molecule_class_info_t::get_all_molecule_rotamer_score(const coot::rotamer_probability_tables &rot_prob_tables) const {

   coot::rotamer_score_t rs;
   double sum_log_p = 0.0;

   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (! model_p) break;
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
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
		  catch (const std::runtime_error &rte) {
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
   mmdb::Residue *residue_p = get_residue(chain_id, resno, ins_code);
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
   mmdb::Residue *residue_p = get_residue(chain_id, resno, ins_code);
   if (residue_p) {
      make_backup();
      int n_atoms = coot::hetify_residue_atoms_as_needed(residue_p);
      if (n_atoms > 0)
	 r = 1;
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
   return r;
}

bool
molecule_class_info_t::has_residue_with_name(const std::string comp_id) const {

   bool r = false;

   if (has_model()) {
      for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 if (! model_p) break;
	 mmdb::Chain *chain_p;
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::Residue *residue_p;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       std::string residue_name = residue_p->GetResName();
	       if (residue_name == comp_id) {
		  r = true;
		  break;
	       }
	    }
	    if (r)
	       break;
	 }
	 if (r)
	    break;
      }
   }
   return r;
}


//////////////////////////////////////////////////////////////////////////////
////////////// animated ligands  /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void
coot::animated_ligand_interactions_t::draw(mmdb::Manager *mol,
					   const gl_context_info_t &gl_info,
					   const long &start_time) const {

#if 0 // I wonder when I will get this working again... many months 8 May 2019.

   mmdb::Atom *at_1 = coot::util::get_atom(ligand_atom_spec, mol);
   mmdb::Atom *at_2 = coot::util::get_atom(interacting_residue_atom_spec, mol);

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

	 long now_time = 0; // glutGet(GLUT_ELAPSED_TIME);
   std::cout << "Fix animated ligand interactions\n";
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

	 GLfloat  mat_specular[]  = {0.18, 0.18, 0.80, opacity};
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

	 for (int ipart=0; ipart<n_parts; ipart++) {

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
	       gluCylinder(quad, base, top, height, slices, stacks);
	       gluDeleteQuadric(quad);
	       glPopMatrix();
	    }
	 }
      }
   }
#endif
}

// add more variables :-)
void molecule_class_info_t::add_animated_ligand_interaction(const pli::fle_ligand_bond_t &lb) {

   animated_ligand_interactions_vec.push_back(lb);
}

void
molecule_class_info_t::draw_animated_ligand_interactions(const gl_context_info_t &gl_info,
							 const long &start_time) const {

   if (draw_animated_ligand_interactions_flag)
      for (unsigned int i=0; i<animated_ligand_interactions_vec.size(); i++) {
	 if (true)
	    std::cout << "---- interaction " << i << " of "
		      << animated_ligand_interactions_vec.size() << std::endl;
	 animated_ligand_interactions_vec[i].draw(atom_sel.mol, gl_info, start_time);
      }

}


std::pair<bool, clipper::Coord_orth>
molecule_class_info_t::residue_centre(const std::string &chain_id, int resno, const std::string &ins_code) const {

   mmdb::Residue *residue_p = get_residue(chain_id, resno, ins_code);
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
molecule_class_info_t::residue_centre(const coot::residue_spec_t &spec) const {

   clipper::Coord_orth p(0,0,0);
   std::pair<bool, clipper::Coord_orth> r(false, p);

   mmdb::Residue *residue_p = get_residue(spec);
   if (residue_p) {
      r = residue_centre(residue_p);
   }
   return r;
}



std::pair<bool, clipper::Coord_orth>
molecule_class_info_t::residue_centre(mmdb::Residue *residue_p) const {

   bool r = 0;
   clipper::Coord_orth pos(0,0,0);
   int n_atoms = 0;

   if (residue_p) {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
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

// return a negative number if not valid
float
molecule_class_info_t::distance_between_residues(mmdb::Residue *r1, mmdb::Residue *r2) const {

   float dist = -1;
   std::pair<bool, clipper::Coord_orth> c1 = residue_centre(r1);
   std::pair<bool, clipper::Coord_orth> c2 = residue_centre(r2);

   if (c1.first && c2.first) {
      dist = clipper::Coord_orth::length(c1.second, c2.second);
   }

   return dist;
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

   std::vector<std::string> ignored_go_to_ligand_residue_types;
   ignored_go_to_ligand_residue_types.push_back("HOH");
   ignored_go_to_ligand_residue_types.push_back("WAT");
   ignored_go_to_ligand_residue_types.push_back("MSE");
   ignored_go_to_ligand_residue_types.push_back("ACE");
   ignored_go_to_ligand_residue_types.push_back("PCA");
   ignored_go_to_ligand_residue_types.push_back("GOL");
   ignored_go_to_ligand_residue_types.push_back("SO4");

   std::vector<std::pair<clipper::Coord_orth, coot::residue_spec_t> > ligand_centres;

   if (atom_sel.n_selected_atoms > 0) {
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
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
	       if (n_atoms >= n_atoms_min) {
		  bool is_het = false;
		  for (int iat=0; iat<n_atoms; iat++) {
		     at = residue_p->GetAtom(iat);
		     if (at->Het) {
			std::string res_name = residue_p->GetResName();

			if (std::find(ignored_go_to_ligand_residue_types.begin(),
				      ignored_go_to_ligand_residue_types.end(),
				      res_name) == ignored_go_to_ligand_residue_types.end()) {
			   is_het = true;
			   break;
			}
		     }
		  }
		  if (is_het) {
		     std::pair<bool, clipper::Coord_orth> res_centre = residue_centre(residue_p);
		     if (res_centre.first) {
              std::pair<clipper::Coord_orth, coot::residue_spec_t> p(res_centre.second, coot::residue_spec_t(residue_p));
              ligand_centres.push_back(p);
		     }
		  }
	       }
	    }
	 }
      }
   }

   int current_centre_index = -1; // unset
   double closest_middle = 9999.9;
   if (ligand_centres.size()) {
      for (unsigned int ilig=0; ilig<ligand_centres.size(); ilig++) {
	 double d = clipper::Coord_orth::length(current_centre, ligand_centres[ilig].first);
	 if (d < 5) {
	    if (d < closest_middle) {
	       current_centre_index = ilig;
	       closest_middle = d;
	    }
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
	    if (current_centre_index == (int(ligand_centres.size())-1)) {
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
							   mmdb::Manager *standard_residues_mol) {
   int status = 0;
   std::vector<mmdb::Residue *> new_residues;
   mmdb::Model *model_p = NULL;
   if (resno_end < resno_start)
      std::swap(resno_start, resno_end);

   for (int ires=resno_start; ires<=resno_end; ires++) {
      mmdb::Residue *res = get_residue(chain_id, ires, "");
      if (!res) {
	 std::cout << "Residue not found in  " << chain_id << " " << ires
		   << std::endl;
      } else {
	 model_p = res->GetModel(); // set it to the model that contains this chain.
	 // these residues are (simple) deep copied
	 mmdb::Residue *res_wc =
	    coot::watson_crick_partner(res, standard_residues_mol);
	 if (res_wc) {
	    new_residues.push_back(res_wc);
	 }
      }
   }

   if (model_p) {
      if (new_residues.size()) {
	 make_backup();
	 mmdb::Chain *chain_p = new mmdb::Chain;
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
	 } else {
	    delete chain_p;
	 }
      }
   }
   return status;
}

// Try to add a LINK record if LINK is not blank.
//
// Passed residue new_res does not have a useful residue_number.
//
std::pair<bool, mmdb::Residue *>
molecule_class_info_t::add_residue(mmdb::Residue *new_res,
				   const std::string &chain_id_in) {

   bool status = false;
   mmdb::Residue *res_copied = NULL;
   int imod = 1;
   if (new_res) {
      if (atom_sel.n_selected_atoms > 0) {
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 if (model_p) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       std::string chain_id(chain_p->GetChainID());
	       if (chain_id == chain_id_in) {
		  make_backup();
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
   }
   return std::pair<bool, mmdb::Residue *> (status, res_copied);
}


// Add a LINK record if link_type is not blank (link_type is for example "NAG-ASN")
//
coot::residue_spec_t
molecule_class_info_t::add_linked_residue_by_beam_in(const coot::residue_spec_t &spec_in,
						     const std::string &new_residue_comp_id,
						     const std::string &link_type,
						     coot::protein_geometry *geom_p) {

   coot::residue_spec_t new_residue_spec;
   mmdb::Residue *residue_ref = get_residue(spec_in);
   if (residue_ref) {
      try {
	 coot::beam_in_linked_residue lr(residue_ref, link_type, new_residue_comp_id, geom_p);
	 mmdb::Residue *result = lr.get_residue();
	 // lr.get_residue() can (and often does) modify residue_ref (by
	 // deleting a link mod atom, for example). So we need a FinishStructEdit() here
	 atom_sel.mol->FinishStructEdit();

	 std::pair<bool, mmdb::Residue *> status_pair = add_residue(result, spec_in.chain_id);

	 if (status_pair.first) {
	    new_residue_spec = coot::residue_spec_t(status_pair.second);
	    coot::dict_link_info_t link_info(residue_ref, status_pair.second,
					     link_type, *geom_p);
	    make_link(link_info.spec_ref, link_info.spec_new, link_type, link_info.dist, *geom_p);
	 }

	 delete result; // no longer needed, we've made (and added) a copy of it.
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }
   return new_residue_spec;
}

#include "coot-utils/glyco-torsions.hh"

// 20140429
//
// Add a LINK record if link_type is not blank (link_type is for example "NAG-ASN")
//
coot::residue_spec_t
molecule_class_info_t::add_linked_residue_by_atom_torsions(const coot::residue_spec_t &spec_in,
							   const std::string &new_residue_comp_id,
							   const std::string &link_type,
							   coot::protein_geometry *geom_p,
							   float default_b_factor_new_atoms) {
   coot::residue_spec_t new_residue_spec;
   mmdb::Residue *residue_ref = get_residue(spec_in);
   if (residue_ref) {
      try {
	 coot::link_by_torsion_t l(link_type, new_residue_comp_id);
	 l.set_temperature_factor(default_b_factor_new_atoms);
	 mmdb::Residue *result = l.make_residue(residue_ref);

	 atom_sel.mol->FinishStructEdit();
	 std::pair<bool, mmdb::Residue *> status_pair = add_residue(result, spec_in.chain_id);
	 if (status_pair.first) {
	    new_residue_spec = coot::residue_spec_t(status_pair.second);
	    coot::dict_link_info_t link_info(residue_ref, status_pair.second, link_type, *geom_p); // exception caught
	    make_link(link_info.spec_ref, link_info.spec_new, link_type, link_info.dist, *geom_p);
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }
   return new_residue_spec;


}


// multi-residue torsion map fitting interface
void
molecule_class_info_t::multi_residue_torsion_fit(const std::vector<coot::residue_spec_t> &residue_specs,
						 const clipper::Xmap<float> &xmap,
						 int n_trials,
						 coot::protein_geometry *geom_p) {

   mmdb::Manager *moving_mol =
      coot::util::create_mmdbmanager_from_residue_specs(residue_specs, atom_sel.mol);

   // do we need to send over the base atom too?  Or just say
   // that it's the first atom in moving_mol?
   //
   std::vector<coot::residue_spec_t> neighbour_specs;
   for (unsigned int i=0; i<residue_specs.size(); i++) {
      std::vector<coot::residue_spec_t> this_res_neighbs = residues_near_residue(residue_specs[i], 8);
      for (unsigned int j=0; j<this_res_neighbs.size(); j++) {
	 if (std::find(neighbour_specs.begin(), neighbour_specs.end(),
		       this_res_neighbs[j]) == neighbour_specs.end()) {
	    if (std::find(residue_specs.begin(), residue_specs.end(), this_res_neighbs[j]) == residue_specs.end()) {
	       neighbour_specs.push_back(this_res_neighbs[j]);
	    }
	 }
      }
   }

   // 20170613 we don't want to include the ASN to which the first NAG is attached into
   // the atoms to be avoided. So (crudely) remove ASN residues from env neighbours
   std::vector<std::pair<bool, clipper::Coord_orth> > avoid_these_atoms;
   for (unsigned int i=0; i<neighbour_specs.size(); i++) {
      mmdb::Residue *residue_p = get_residue(neighbour_specs[i]);
      if (residue_p) {
	 std::string res_name = residue_p->GetResName();
	 if (res_name != "ASN") {
	    mmdb::Atom **residue_atoms = 0;
	    int n_residue_atoms;
	    residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       mmdb::Atom *at = residue_atoms[iat];
	       clipper::Coord_orth pt = coot::co(at);

	       // if an atom of the residue specs is right now close to a neighbour residue
	       // atom, then that neighbour is likely to be bonded (or angle-related)
	       // so don't add such atoms to the avoid_these_atoms - yes, this is a
	       // hack, but it's better than it was.
	       //
	       bool already_close = false;
	       //
	       mmdb::Model *model_p = moving_mol->GetModel(1);
	       if (model_p) {
		  int n_chains = model_p->GetNumberOfChains();
		  for (int ichain=0; ichain<n_chains; ichain++) {
		     mmdb::Chain *chain_p = model_p->GetChain(ichain);
		     int nres = chain_p->GetNumberOfResidues();
		     for (int ires=0; ires<nres; ires++) {
			mmdb::Residue *moving_residue_p = chain_p->GetResidue(ires);
			int n_atoms = moving_residue_p->GetNumberOfAtoms();
			for (int iat=0; iat<n_atoms; iat++) {
			   mmdb::Atom *moving_at = moving_residue_p->GetAtom(iat);
			   clipper::Coord_orth pt_moving = coot::co(moving_at);
			   if ((pt - pt_moving).lengthsq() < 2.8*2.8) {
			      already_close = true;
			      break;
			   }
			}
			if (already_close)
			   break;
		     }
		  }
	       }

	       if (! already_close) {
		  std::string ele = at->element;
		  if (ele != " H") { // PDBv3 fixme
		     bool t = false;
		     std::string res_name = at->GetResName();
		     if (res_name == "HOH") t = true;
		     std::pair<bool, clipper::Coord_orth> p(t, pt);
		     avoid_these_atoms.push_back(p);
		  }
	       }
	    }
	 }
      }
   }

   coot::multi_residue_torsion_fit_map(imol_no, moving_mol, xmap, avoid_these_atoms, n_trials, geom_p);

   atom_selection_container_t moving_atoms_asc = make_asc(moving_mol);
   replace_coords(moving_atoms_asc, 1, 1);

}



coot::residue_spec_t
molecule_class_info_t::get_residue_by_type(const std::string &residue_type) const {

   coot::residue_spec_t spec;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    std::string residue_name(residue_p->GetResName());
	    if (residue_name == residue_type) {
	       spec = coot::residue_spec_t(residue_p);
	       break;
	    }
	 }
	 if (! spec.unset_p())
	    break;
      }
   }

   return spec;
}

std::vector<coot::residue_spec_t>
molecule_class_info_t::get_residues_by_type(const std::string &residue_type) const {

   std::vector<coot::residue_spec_t> v;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    std::string residue_name(residue_p->GetResName());
	    if (residue_name == residue_type) {
	       coot::residue_spec_t spec(residue_p);
	       v.push_back(spec);
	    }
	 }
      }
   }
   return v;
}



void
molecule_class_info_t::split_water(std::string chain_id, int res_no, std::string ins_code,
				   const clipper::Xmap<float> &xmap,
				   float map_sigma) {

   bool modified = false;
   coot::residue_spec_t rs(chain_id, res_no, ins_code);
   mmdb::Residue *residue_p = get_residue(rs);
   if (residue_p) {
      int n_atoms = residue_p->GetNumberOfAtoms();
      if (n_atoms == 1) {
	 mmdb::Atom *at = residue_p->GetAtom(" O  "); // PDBv3
	 if (at) {
	    make_backup();
        float old_occ = at->occupancy;
	    mmdb::Atom *new_at = new mmdb::Atom;
	    new_at->Copy(at);
	    // force down a new altconf
	    residue_p->AddAtom(new_at);
	    strncpy(    at->altLoc, "A", 2);
	    strncpy(new_at->altLoc, "B", 2);
	    at->x     -= 0.5;
	    new_at->x += 0.5;
	    at->occupancy     = old_occ * 0.5;
	    new_at->occupancy = old_occ * 0.5;
	    modified = true;
	    atom_sel.mol->FinishStructEdit();
	    update_molecule_after_additions(); // needed to update the atom_index user data
	                                       // (used in full_atom_spec_to_atom_index()
	                                       // in fit_to_map_by_random_jiggle().
	 }
      }
   }

   if (modified) {

      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<mmdb::Chain *> chains; // empty
      fit_to_map_by_random_jiggle(residue_atoms, n_residue_atoms, xmap, map_sigma, 10, 1, false, chains);

      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
   }
}

std::vector<coot::residue_spec_t>
molecule_class_info_t::all_residues() const {

   std::vector<coot::residue_spec_t> r;

   if (!atom_sel.mol) return r;
   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
	 mmdb::Chain *chain_p;
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::Residue *residue_p;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       if (residue_p) {
		  coot::residue_spec_t spec(residue_p);
		  spec.int_user_data = ires;
		  r.push_back(spec);
	       }
	    }
	 }
      }
   }
   return r;
}

std::vector<mmdb::Residue *>
molecule_class_info_t::get_all_protein_residues() const {

   std::vector<mmdb::Residue *> v;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               v.push_back(residue_p);
            }
         }
      }
   }

   return v;
}



std::vector<coot::residue_spec_t>
molecule_class_info_t::het_groups() const {

   std::vector<coot::residue_spec_t> r;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    std::string rn = residue_p->GetResName();
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    for (int iat=0; iat<n_atoms; iat++) {
	       mmdb::Atom *at = residue_p->GetAtom(iat);
	       if (at->Het) {
		  if (rn != "HOH") {
		     r.push_back(coot::residue_spec_t(residue_p));
		     break;
		  }
	       }
	    }
	 }
      }
   }
   return r;
}

std::string
molecule_class_info_t::get_sequence_as_block(const std::string &chain_id) const {

   bool with_spaces = false; // block spaced output is easier to read
   std::string seq;
   if (atom_sel.mol) {
      mmdb::Manager *mol = atom_sel.mol;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (std::string(chain_p->GetChainID()) == chain_id) { 
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::PResidue residue_p;
	    int residue_count_block = 0;
	    int residue_count_line = 0;
	    if (nres > 0 ) {
	       residue_count_block = chain_p->GetResidue(0)->GetSeqNum();
	       residue_count_line  = residue_count_block;
	       if (residue_count_block > 0)
		  while (residue_count_block > 10)
		     residue_count_block -= 10;
	       if (residue_count_line > 0)
		  while (residue_count_line > 50)
		     residue_count_line -= 50;
	    }
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       seq += coot::util::three_letter_to_one_letter(residue_p->GetResName());
	       if (residue_count_block == 10) {
		  if (with_spaces)
		     seq += " ";
		  residue_count_block = 0;
	       }
	       if (residue_count_line == 50) {
		  seq += "\n";
		  residue_count_line = 0;
	       }
	       residue_count_block++;
	       residue_count_line++;
	    }
	 }
      }
   }
   return seq;

}


// the length of the string is guaranteed to the the length of the vector
std::pair<std::string, std::vector<mmdb::Residue *> >
molecule_class_info_t::sequence_from_chain(mmdb::Chain *chain_p) const {

   mmdb::PResidue *residues = 0;
   int n_residues;
   chain_p->GetResidueTable(residues, n_residues);
   std::string s;
   std::vector<mmdb::Residue *> v;
   char r[10];

   if (residues) {
      for (int i=0; i<n_residues; i++) {
	 std::string this_residue = "X";
	 mmdb::pstr rn = residues[i]->GetResName();
	 std::string residue_name(residues[i]->GetResName());
	 mmdb::Get1LetterCode(rn, r); // mmdb
	 this_residue = r[0];
	 if (residue_name != "HOH") {
	    s += this_residue;
	    v.push_back(residues[i]);
	 }
      }
   }
   return std::pair<std::string, std::vector<mmdb::Residue *> >(s,v);
}


// return null on failure.  seq_trip is something like "ACE".
mmdb::Atom *
molecule_class_info_t::get_centre_atom_from_sequence_triplet(const std::string &seq_trip) const {

   mmdb::Atom *at = 0;
   std::string useq_trip = coot::util::upcase(seq_trip);

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (! model_p) return at;
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   int n_same_triplet = 0;
   std::map<mmdb::Chain *, std::vector<coot::residue_spec_t> > triplet_map;
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      std::pair<std::string, std::vector<mmdb::Residue *> > seq = sequence_from_chain(chain_p);
      std::string::size_type ifound = 0; // initial non-npos value
      std::string::size_type ifound_prev = 0; // initial non-npos value
      std::string::size_type i_offset = 0;

      while (ifound != std::string::npos) {
	 ifound = seq.first.substr(ifound_prev+i_offset).find(useq_trip);
	 if (ifound != std::string::npos)
	    i_offset++;
	 // 	 std::cout << "finding " << useq_trip << " in " << seq.first.substr(ifound_prev+i_offset)
	 // 		   << " yields " << ifound << std::endl;
	 if (ifound != std::string::npos) {
	    // 	    std::cout << "... which is not npos " << std::endl;
	    int idx = ifound_prev+ifound+i_offset;  // middle residue
	    if (idx < int(seq.second.size())) {
	       // this should always be so
	       mmdb::Residue *r = seq.second[idx];
	       int iat = intelligent_this_residue_atom(r);
	       if (iat != -1) {
		  // std::cout << "adding spec " << coot::residue_spec_t(atom_sel.atom_selection[iat])
		  // << std::endl;
		  triplet_map[chain_p].push_back(coot::residue_spec_t(coot::atom_spec_t(atom_sel.atom_selection[iat])));
		  if (at == 0)
		     at = atom_sel.atom_selection[iat];
	       }
	    } else {
	       std::cout << "ERROR:: !!! .... out of range index " << idx << " not < " << seq.second.size()
			 << std::endl;
	    }
	 }
	 ifound_prev += ifound;
      }
   }

   if (at) {
      coot::residue_spec_t at_rs(at->residue);
      int n_ncs = 0;
      std::map<mmdb::Chain *, std::vector<coot::residue_spec_t> >::const_iterator it;
      for (it=triplet_map.begin(); it!=triplet_map.end(); it++) {
	 for (unsigned int i=0; i<it->second.size(); i++) {
	    const coot::residue_spec_t vi_spec(it->second[i]);
	    // std::cout << "comparing " << it->second[i] << " and " << vi_spec << std::endl;
	    if (vi_spec == at_rs) {
	       // this is the same atom, skip it.
	    } else {
	       n_same_triplet++;
	    }
	 }
      }
      std::cout << "INFO:: " << "centred on first occurance of triplet " << useq_trip << ", there are "
		<< n_same_triplet+1 << " in total" << std::endl;
      // if there are triplets that are not NCS related, write out all triplets
      if ((n_same_triplet) > 0) {
	 std::cout << "------- (middle) residues matching " << useq_trip << " --------" << std::endl;
	 for (it=triplet_map.begin(); it!=triplet_map.end(); it++) {
	    for (unsigned int i=0; i<it->second.size(); i++) {
	       const coot::residue_spec_t vi_spec(it->second[i]);
	       std::cout << "   " << vi_spec<< std::endl;
	    }
	 }
      }
   }
   return at;
}


void
molecule_class_info_t::rotate_residue(const coot::residue_spec_t &rs,
				      const clipper::Coord_orth &around_vec,
				      const clipper::Coord_orth &origin_offset,
				      double angle) {

   mmdb::Residue *residue_p = get_residue(rs);
   if (residue_p) {
      make_backup();
      coot::util::rotate_residue(residue_p, around_vec, origin_offset, angle);
      have_unsaved_changes_flag = 1;
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked(__FUNCTION__);
   }
}

std::vector<coot::fragment_info_t>
molecule_class_info_t::get_fragment_info(bool screen_output_also) const {

   std::vector<coot::fragment_info_t> v;

   if (! atom_sel.mol) return v;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (! model_p) return v;

   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      if (nres > 0) {
	 coot::fragment_info_t fi(chain_p->GetChainID());
	 mmdb::Residue *residue_p_prev  = 0;
	 mmdb::Residue *residue_p_start = 0;
	 mmdb::Residue *residue_p_this;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p_this = chain_p->GetResidue(ires);
	    // if this was the second residue or further along...
	    if (residue_p_prev) {
	       if ((residue_p_prev->GetSeqNum() + 1) != residue_p_this->GetSeqNum()) {
                  coot::residue_spec_t spec_1(residue_p_start);
                  coot::residue_spec_t spec_2(residue_p_prev);
		  coot::fragment_info_t::fragment_range_t fr(spec_1, spec_2);
		  fi.add_range(fr);
		  residue_p_start = residue_p_this; // start a new fragment (part-way through a chain)
	       }
	    } else {
	       // we are starting a new fragment (at the beginning of a chain)
	       residue_p_start = residue_p_this;
	    }
	    residue_p_prev = residue_p_this;
	 }

	 // and the last fragment of the chain
	 if (residue_p_start) {
            coot::residue_spec_t spec_1(residue_p_start);
            coot::residue_spec_t spec_2(residue_p_this);
	    coot::fragment_info_t::fragment_range_t r(spec_1, spec_2);
	    fi.add_range(r);
	 }

	 // done this fragment info:
	 v.push_back(fi);
      }
   }

   if (screen_output_also) {
      std::cout << "------------------" << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
	 std::cout << "   chain-id: " << v[i].chain_id << std::endl;
	 for (unsigned int j=0; j<v[i].ranges.size(); j++) {
	    std::cout << "      "
		      << v[i].ranges[j].start_res.res_no << " "
		      << v[i].ranges[j].end_res.res_no << std::endl;
	 }
      }
      std::cout << "------------------" << std::endl;
   }
   return v;
}

// return a new residue and dictionary in due course.
// return a null residue on failure
std::pair<mmdb::Residue *, coot::dictionary_residue_restraints_t>
molecule_class_info_t::invert_chiral_centre(const std::string &chain_id, int res_no,
					    const std::string &ins_code,
					    const std::string &atom_name,
					    const coot::protein_geometry &geom) {

   mmdb::Atom *chiral_atom = NULL;
   mmdb::Residue *residue_p = get_residue(chain_id, res_no, ins_code);

   coot::dictionary_residue_restraints_t new_restraints;
   mmdb::Residue *new_residue = 0;

   if (residue_p) {

      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 if (std::string(at->name) == atom_name) {
	    chiral_atom = at;
	    break;
	 }
      }

      if (chiral_atom) {
	 std::string comp_id = residue_p->GetResName();
	 std::pair<bool, coot::dictionary_residue_restraints_t> rest =
	    geom.get_monomer_restraints(comp_id, imol_no);
	 if (rest.first) {
	    std::vector<coot::dict_chiral_restraint_t> cr = rest.second.chiral_restraint;
	    for (unsigned int ir=0; ir<cr.size(); ir++) {
	       if (cr[ir].atom_id_c_4c() == atom_name) {
		  if (cr[ir].volume_sign == coot::dict_chiral_restraint_t::CHIRAL_RESTRAINT_POSITIVE ||
		      cr[ir].volume_sign == coot::dict_chiral_restraint_t::CHIRAL_RESTRAINT_NEGATIVE) {
		     std::string new_name = coot::suggest_new_comp_id(rest.second.residue_info.comp_id);
		     if (! new_name.empty()) {
			new_restraints = rest.second;
			new_restraints.residue_info.comp_id = new_name;
			new_restraints.residue_info.three_letter_code = new_name;
			new_restraints.residue_info.name = coot::util::remove_trailing_whitespace(new_restraints.residue_info.name);
			new_restraints.residue_info.name += " with inverted ";
			new_restraints.residue_info.name += atom_name;
			new_restraints.residue_info.name += " chiral centre";
			new_restraints.chiral_restraint[ir].invert_target_volume();
			new_residue = coot::util::deep_copy_this_residue(residue_p);
			new_residue->SetResName(new_name.c_str());
			break;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return std::pair<mmdb::Residue *, coot::dictionary_residue_restraints_t> (new_residue,
									new_restraints);

}

void
molecule_class_info_t::residue_partial_alt_locs_split_residue(coot::residue_spec_t spec,
							      int i_bond,
							      double theta, // degrees
							      bool wag_the_dog,
							      coot::protein_geometry *geom_p) {

   mmdb::Residue *residue_p = get_residue(spec);
   if (residue_p) {
      bool found = false;
      std::string residue_type = residue_p->GetResName();
      std::pair<short int, coot::dictionary_residue_restraints_t> r =
	 geom_p->get_monomer_restraints(residue_type, imol_no);

      if (r.first) {
	 bool find_hydrogen_torsions_flag = false;
	 std::vector <coot::dict_torsion_restraint_t> torsion_restraints =
	    r.second.get_non_const_torsions(find_hydrogen_torsions_flag);

	 if (i_bond >= 0 && i_bond < int(torsion_restraints.size())) {

	    std::pair<std::string, std::string> atom_names;
	    atom_names.first  = torsion_restraints[i_bond].atom_id_2_4c();
	    atom_names.second = torsion_restraints[i_bond].atom_id_3_4c();
	    if ((atom_names.first != "") && (atom_names.second != "")) {
	       mmdb::PPAtom residue_atoms;
	       int nResidueAtoms;
	       residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
	       if (nResidueAtoms > 0) {
		  for (int iat1=0; iat1<nResidueAtoms; iat1++) {
		     std::string ra1=residue_atoms[iat1]->name;
		     if (ra1 == atom_names.first) {
			std::string alt_conf = residue_atoms[iat1]->altLoc;

			if (alt_conf.empty()) { // don't split an already-split residue

			   for (int iat2=0; iat2<nResidueAtoms; iat2++) {
			      std::string ra2=residue_atoms[iat2]->name;
			      if (ra2 == atom_names.second) {

				 // OK we have both atoms.
				 // now create new atoms and spin them around the bond
				 found = true;

				 std::string monomer_type = residue_p->GetResName();
				 atom_selection_container_t residue_asc; // refactor: make a constructor
				 residue_asc.n_selected_atoms = nResidueAtoms;
				 residue_asc.atom_selection = residue_atoms;
				 residue_asc.mol = 0;
				 coot::contact_info contact = coot::getcontacts(residue_asc, monomer_type, imol_no, geom_p);
				 // contact.print(); // debug
				 // std::vector<std::vector<int> > contact_indices = contact.get_contact_indices();
				 // or this?
				 std::vector<std::vector<int> > contact_indices =
				    contact.get_contact_indices_with_reverse_contacts();

				 if (false) { // debug - we needed reverse contacts
				    std::cout << "-------- contact_indices ------" << std::endl;
				    for (unsigned int ii=0; ii<contact_indices.size(); ii++) {
				       for (unsigned int jj=0; jj<contact_indices[ii].size(); jj++) {
					  std::cout << ii << "   " << contact_indices[ii][jj] << " ";
				       }
				       std::cout << std::endl;
				    }
				 }

				 int base_atom_index = 0;

				 try {
				    //
				    coot::atom_tree_t tree(contact_indices, base_atom_index, residue_p, alt_conf);

				    bool reverse_fragment = false;
				    std::pair<unsigned int, unsigned int> sizes = tree.fragment_sizes(ra1, ra2, false);

				    // rotate the smallest fragment by default
				    if (sizes.first > sizes.second)
				       reverse_fragment = true;

				    // now consider user input:
				    if (wag_the_dog)
				       reverse_fragment = !reverse_fragment;


				    // we want to make a copy of the atoms before they were moved so that we can
				    // set the alt conf and position of them when we make a copy
				    std::vector<std::pair<mmdb::Atom *, clipper::Coord_orth> > atom_store;
				    std::vector<int> moving_atoms =
				       tree.get_moving_atom_indices(ra1, ra2, reverse_fragment);
				    for (unsigned int ii=0; ii<moving_atoms.size(); ii++) {
				       mmdb::Atom *at = residue_atoms[moving_atoms[ii]];
				       clipper::Coord_orth co = coot::co(at);
				       std::pair<mmdb::Atom *, clipper::Coord_orth> p(at, co);
				       atom_store.push_back(p);
				    }
				    tree.rotate_about(ra1, ra2, theta, reverse_fragment);

				    if (atom_store.size()) {
				       make_backup();
				       std::cout << "These are the moving atoms: " << std::endl;
				       for (unsigned int ii=0; ii<atom_store.size(); ii++) {
					  std::cout << "here 0 " << ii << " of " << atom_store.size() << std::endl;
					  mmdb::Atom *at = atom_store[ii].first;
					  const clipper::Coord_orth &pos = atom_store[ii].second;
					  if (true)
					     std::cout << "   " << coot::atom_spec_t(at) << std::endl;
					  strcpy(at->altLoc, "A");
					  at->occupancy = 0.5;
					  mmdb::Atom *at_B = new mmdb::Atom;
					  *at_B = *at;
					  at_B->SetCoordinates(pos.x(), pos.y(), pos.z(), 0.5, at->tempFactor);
					  strcpy(at_B->altLoc, "B");
					  std::cout << "adding atom " << at_B << std::endl;
					  residue_p->AddAtom(at_B);
					  std::cout << "done adding atom " << at_B << std::endl;
				       }

				       have_unsaved_changes_flag = 1;
				       atom_sel.mol->FinishStructEdit();
				       atom_sel = make_asc(atom_sel.mol);
				       make_bonds_type_checked(__FUNCTION__);
				    }
				 }

				 catch (const std::runtime_error &rte) {
				    std::cout << "ERROR:: runtime_error: " << rte.what() << std::endl;
				 }
				 catch (const std::exception &rte) {
				    std::cout << "ERROR:: exception: " << rte.what() << std::endl;
				 }
			      }
			      if (found)
				 break;
			   }
			   if (found)
			      break;
			}
		     }
		     if (found)
			break;
		  }
	       }
	    }
	 }
      }
   }
}


int
molecule_class_info_t::delete_chain(const std::string &chain_id) {

   int done = false;
   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    if (chain_p) {
	       std::string this_chain_id(chain_p->GetChainID());
	       if (this_chain_id == chain_id) {
                  make_backup();
		  model_p->DeleteChain(ichain);
		  done = true;
	       }
	    }
	 }
      }
   }

   if (done) {
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
   }

   return done;

}

int
molecule_class_info_t::delete_sidechains_for_chain(const std::string &chain_id) {

   int done = false;

   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    if (chain_p) {
	       std::string this_chain_id = chain_p->GetChainID();
	       if (this_chain_id == chain_id) {

		  make_backup();

		  int nres = chain_p->GetNumberOfResidues();
		  // delete the specific atoms of eacho of the residues:
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::PResidue residue_p = chain_p->GetResidue(ires);
		     if (residue_p) {
			mmdb::PPAtom atoms = 0;
			int n_atoms = 0;
			bool was_deleted = false;
			residue_p->GetAtomTable(atoms, n_atoms);
			for (int i=0; i<n_atoms; i++) {
			   if (! (coot::is_main_chain_or_cb_p(atoms[i]))) {
			      residue_p->DeleteAtom(i);
			      was_deleted = true;
			   }
			}
			if (was_deleted)
			   residue_p->TrimAtomTable();
		     }
		  }
		  done = true;
	       }
	    }
	 }
      }
   }

   if (done) {
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
   }
   return done;
}

bool
molecule_class_info_t::delete_sidechain(mmdb::Residue *residue_p) {

   bool was_deleted = false;

   if (residue_p) {
      // make_backup();
      mmdb::PPAtom atoms = 0;
      int n_atoms = 0;
      bool was_deleted = false;
      residue_p->GetAtomTable(atoms, n_atoms);
      for (int i=0; i<n_atoms; i++) {
         if (! (coot::is_main_chain_or_cb_p(atoms[i]))) {
            residue_p->DeleteAtom(i);
            was_deleted = true;
         }
      }
      if (was_deleted)
         residue_p->TrimAtomTable();
   }
   return was_deleted;
}

int
molecule_class_info_t::delete_sidechain_range(const coot::residue_spec_t &res_1,
					      const coot::residue_spec_t &res_2) {

   int done = false;

   int resno_min = res_1.res_no;
   int resno_max = res_2.res_no;

   if (resno_max < resno_min) {
      resno_max = res_1.res_no;
      resno_min = res_2.res_no;
   }

   std::string chain_id = res_1.chain_id;
   if (res_2.chain_id != res_1.chain_id)
      return done; // nothing doing.

   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    if (chain_p) {
	       std::string this_chain_id = chain_p->GetChainID();
	       if (this_chain_id == chain_id) {

		  make_backup();

		  int nres = chain_p->GetNumberOfResidues();
		  // delete the specific atoms of eacho of the residues:
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::PResidue residue_p = chain_p->GetResidue(ires);
		     if (residue_p) {
			mmdb::PPAtom atoms = 0;
			int n_atoms = 0;
			bool was_deleted = false;
			int resno_this = residue_p->GetSeqNum();
			if (resno_this >= resno_min) {
			   if (resno_this <= resno_max) {
			      residue_p->GetAtomTable(atoms, n_atoms);
			      for (int i=0; i<n_atoms; i++) {
				 if (! (coot::is_main_chain_or_cb_p(atoms[i]))) {
				    residue_p->DeleteAtom(i);
				    was_deleted = true;
				 }
			      }
			      if (was_deleted)
				 residue_p->TrimAtomTable();
			   }
			}
		     }
		  }
		  done = true;
	       }
	    }
	 }
      }
   }

   if (done) {
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
   }
   return done;
}

void
molecule_class_info_t::delete_all_carbohydrate() {

   if (atom_sel.mol) {
      make_backup();
      coot::util::delete_all_carbohydrate(atom_sel.mol);
      make_bonds_type_checked();
   }
}



// carbohydrate validation tools
//
// should pass the cif read number pointer
//
// should this be const?
void
molecule_class_info_t::glyco_tree_internal_distances_fn(const coot::residue_spec_t &residue_spec,
							coot::protein_geometry *geom_p,
							const std::string &file_name) {

   if (atom_sel.mol) {
      mmdb::Manager *mol = atom_sel.mol;
      mmdb::Residue *residue_p = get_residue(residue_spec);
      if (residue_p) {
	 int mmcif_read_number = 51;
	 std::vector<std::string> types_with_no_dictionary =
	    no_dictionary_for_residue_type_as_yet(*geom_p);
	 for (unsigned int i=0; i<types_with_no_dictionary.size(); i++)
	    geom_p->try_dynamic_add(types_with_no_dictionary[i], mmcif_read_number++);
	 coot::glyco_tree_t t(residue_p, mol, geom_p);
	 double dist_lim = 20;
	 t.internal_distances(dist_lim, file_name);
      }
   }
}

#include "coot-utils/secondary-structure-headers.hh"

void
molecule_class_info_t::add_secondary_structure_header_records(bool overwrite) {

   // if there is secondary structure already, don't overwrite it.
   bool do_it = true;
   if (atom_sel.mol) {
      if (! overwrite) {
	 mmdb::Model *model_p = atom_sel.mol->GetModel(1);
	 int nhelix = model_p->GetNumberOfHelices();
	 int nsheet = model_p->GetNumberOfSheets();
	 if ((nhelix > 0) || (nsheet > 0)) {
	    do_it = false;
	 }
      }
      if (do_it) {
	 int n_models = atom_sel.mol->GetNumberOfModels();
	 for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
	    mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	    int ss_status = model_p->CalcSecStructure(1);
	    coot::secondary_structure_header_records ssr(atom_sel.mol, false);
	    if (ss_status == mmdb::SSERC_Ok) {
	       std::cout << "INFO:: SSE status was OK\n";
	    } else {
	       std::cout << "INFO:: SSE status was not OK\n";
	    }
	 }
      }
   }
}



// Spin N and the sidechain around the CA-C vector
//
// useful for add terminal residue N-terminal addition extension
//
// angle in degrees.
void molecule_class_info_t::spin_N(const coot::residue_spec_t &residue_spec, float angle) {

   mmdb::Residue *residue_p = get_residue(residue_spec);
   if (residue_p) {

      double a = clipper::Util::d2rad(angle);
      // PDBv3 FIXME
      coot::atom_spec_t ca_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " CA ", "");
      coot::atom_spec_t  c_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " C  ", "");
      coot::atom_spec_t  o_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " O  ", "");
      mmdb::Atom *ca = coot::util::get_atom(ca_spec, residue_p);
      mmdb::Atom *c  = coot::util::get_atom( c_spec, residue_p);
      mmdb::Atom *o  = coot::util::get_atom( o_spec, residue_p);
      if (ca && c && o) {
	 make_backup();
	 clipper::Coord_orth ca_pos = coot::co(ca);
	 clipper::Coord_orth  c_pos = coot::co(c);
	 clipper::Coord_orth dir = ca_pos - c_pos;
	 clipper::Coord_orth origin_shift = c_pos;
	 mmdb::Atom **residue_atoms = 0;
	 int n_residue_atoms;
	 residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    mmdb::Atom *at = residue_atoms[iat];
	    if (at) {
	       if (at != ca && at!=c && at!=o) {
		  clipper::Coord_orth pos = coot::co(at);
		  clipper::Coord_orth pt = coot::util::rotate_around_vector(dir, pos, origin_shift, a);
		  coot::update_position(at, pt);
	       }
	    }
	 }

	 have_unsaved_changes_flag = 1;
	 atom_sel.mol->FinishStructEdit();
	 atom_sel = make_asc(atom_sel.mol);
	 make_bonds_type_checked(__FUNCTION__);

      }
   }
}


// place the O (because we have added a new residue)
bool
molecule_class_info_t::move_atom(const std::string &atom_name_in, mmdb::Residue *res_p, const clipper::Coord_orth &new_O_pos) {

   // Hmm! we are passed res_p - that's unusual.

   // just change the position of the first atom that matches atom_name_in
   bool done = false;

   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms;
   res_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      mmdb::Atom *at = residue_atoms[i];
      std::string atom_name(at->name);
      if (atom_name == atom_name_in) {
	 at->x = new_O_pos.x();
	 at->y = new_O_pos.y();
	 at->z = new_O_pos.z();
	 have_unsaved_changes_flag = 1;
	 atom_sel.mol->FinishStructEdit();
	 atom_sel = make_asc(atom_sel.mol);
	 make_bonds_type_checked(__FUNCTION__);
	 done = true;
	 break;
      }
   }
   return done;
}

#include "coot-utils/merge-atom-selections.hh"

// merge change/fragments of this molecule
// return 1 if a merge was done;
int
molecule_class_info_t::merge_fragments() {

   int status = 1;

   make_backup();

   coot::merge_atom_selections(atom_sel.mol); // doesn't return a value, should it?

   have_unsaved_changes_flag = 1;
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   make_bonds_type_checked(__FUNCTION__);

   return status;

}


std::vector<mmdb::Atom *>
molecule_class_info_t::closest_atoms_in_neighbour_residues(mmdb::Residue *residue_ref_p, float radius) const {

   std::vector<mmdb::Atom *> v;

   // we don't want to include main chain atoms that are peptide bonded to this one (residue_ref_p)

   if (residue_ref_p) {
      std::vector<mmdb::Residue *> rv = coot::residues_near_residue(residue_ref_p, atom_sel.mol, radius);
      for (unsigned int i=0; i<rv.size(); i++) {
         mmdb::Residue *r = rv[i];
         int n_residue_atoms = 0;
         mmdb::Atom **residue_atoms;
         float closest_dist = 99999.9;
         mmdb::Atom *closest_atom = 0;
         r->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int j=0; j<n_residue_atoms; j++) {
            mmdb::Atom *at = residue_atoms[j];
            if (! at->isTer()) {
               clipper::Coord_orth co = coot::co(at);
               mmdb::Atom **residue_ref_atoms = 0;
               int n_residue_ref_atoms = 0;
               residue_ref_p->GetAtomTable(residue_ref_atoms, n_residue_ref_atoms);
               for (int k=0; k<n_residue_ref_atoms; k++) {
                  mmdb::Atom *at_ref = residue_ref_atoms[k];
                  if (! at_ref->isTer()) {
                     clipper::Coord_orth co_ref = coot::co(at_ref);
                     double dd = (co-co_ref).lengthsq();
                     double d = sqrt(dd);
                     if (d < closest_dist) {
                        bool exclude = false;
                        if (coot::is_main_chain_p(at)) {
                           if (residue_ref_p->chain == r->chain) {
                              int res_no_delta = r->GetSeqNum() - residue_ref_p->GetSeqNum();
                              if (abs(res_no_delta) < 2)
                                 exclude = true;
                           }
                        }
                        if (! exclude) {
                           closest_dist = d;
                           closest_atom = at;
                        }
                     }
                  }
               }
            }
         }
         if (closest_atom) {
            v.push_back(closest_atom);
         }
      }
   }

   std::cout << "debug:: got " << v.size() << " closest atoms " << std::endl;

   return v;
}

void
molecule_class_info_t::label_closest_atoms_in_neighbour_atoms(coot::residue_spec_t residue_spec, float radius) {

   mmdb::Residue *residue_ref_p = get_residue(residue_spec);
   if (residue_ref_p) {
      int atom_index_handle = atom_sel.UDDAtomIndexHandle; // or SelectionHandle?
      std::vector<mmdb::Atom *>v = closest_atoms_in_neighbour_residues(residue_ref_p, radius);
      for (unsigned int i=0; i<v.size(); i++) {
         mmdb::Atom *at = v[i];
         int atom_index = -1;
         at->GetUDData(atom_index_handle, atom_index);
         if (atom_index >= 0)
            if (atom_index < atom_sel.n_selected_atoms)
               labelled_atom_index_list.push_back(atom_index);
      }
   }
}


#include "coot-utils/atom-overlaps.hh"

// this is something that you do after alignment.
void
molecule_class_info_t::resolve_clashing_sidechains_by_deletion(const coot::protein_geometry *geom_p) {

   auto matches_existing_pair = [] (const std::pair<mmdb::Residue *, mmdb::Residue *> &p,
                                    const std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > &v) {
                                   for (unsigned int i=0; i<v.size(); i++) {
                                      if (p.first == v[i].first)
                                         if (p.second == v[i].second)
                                            return true;
                                      if (p.first == v[i].second)
                                         if (p.second == v[i].first)
                                            return true;
                                   }
                                   return false;
                                };

   const std::vector<std::string> sized_residues = {
                                                    "GLY", "ALA", "SER", "CYS", "PRO", "THR",
                                                    "VAL", "ASP", "ASN", "GLU", "GLN", "HIS",
                                                    "ILE", "LEU", "LYS", "MET", "ARG", "PHE",
                                                    "TYR", "TRP", "MSE" };
   auto get_index = [] (const std::string &rt, const std::vector<std::string> &v) {
                       int idx = -1;
                       int vs = v.size();
                       for (int i=0; i<vs; i++) {
                          if (rt == v[i])
                             return i;
                       }
                       return idx;
                    };

      auto residue_is_bigger = [sized_residues, get_index] (mmdb::Residue *r1, mmdb::Residue *r2) {
                               bool status = false;
                               std::string rt_1 = r1->GetResName();
                               std::string rt_2 = r2->GetResName();
                               int idx_1 = get_index(rt_1, sized_residues);
                               int idx_2 = get_index(rt_2, sized_residues);
                               if (idx_1 > idx_2) status = true;
                               return status;
                            };

   bool ignore_waters = true;
   mmdb::Manager *mol = atom_sel.mol;
   coot::atom_overlaps_container_t overlaps(mol, geom_p, ignore_waters, 0.5, 0.25);
   std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > overlapping_residues;
   overlaps.make_all_atom_overlaps();
   std::vector<coot::atom_overlap_t> olv = overlaps.overlaps;
   for (std::size_t ii=0; ii<olv.size(); ii++) {
      const coot::atom_overlap_t &o = olv[ii];
      if (o.overlap_volume > 2.0) {
         mmdb::Residue *r1 = o.atom_1->residue;
         mmdb::Residue *r2 = o.atom_2->residue;
         std::pair<mmdb::Residue *, mmdb::Residue *>  p(r1, r2);
         if (! matches_existing_pair(p, overlapping_residues))
            overlapping_residues.push_back(p);
      }
   }

   bool things_have_changed = false;
   for (unsigned int ii=0; ii<overlapping_residues.size(); ii++) {
      mmdb::Residue *residue_p = overlapping_residues[ii].first;
      if (residue_is_bigger(overlapping_residues[ii].second, residue_p))
         residue_p = overlapping_residues[ii].second;
      delete_sidechain(residue_p);
      things_have_changed = true;
   }
   if (things_have_changed) {
      have_unsaved_changes_flag = 1;
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked(__FUNCTION__);
   }
}


void
molecule_class_info_t::resolve_clashing_sidechains_by_rebuilding(const coot::protein_geometry *geom_p,
                                                                 int imol_refinement_map) {

   auto matches_existing_pair = [] (const std::pair<mmdb::Residue *, mmdb::Residue *> &p,
                                    const std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > &v) {
                                   for (unsigned int i=0; i<v.size(); i++) {
                                      if (p.first == v[i].first)
                                         if (p.second == v[i].second)
                                            return true;
                                      if (p.first == v[i].second)
                                         if (p.second == v[i].first)
                                            return true;
                                   }
                                   return false;
                                };

   bool ignore_waters = true;
   mmdb::Manager *mol = atom_sel.mol;
   coot::atom_overlaps_container_t overlaps(mol, geom_p, ignore_waters, 0.5, 0.25);
   std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > overlapping_residues;
   overlaps.make_all_atom_overlaps();
   std::vector<coot::atom_overlap_t> olv = overlaps.overlaps;
   for (std::size_t ii=0; ii<olv.size(); ii++) {
      const coot::atom_overlap_t &o = olv[ii];
      if (o.overlap_volume > 2.0) {
         mmdb::Residue *r1 = o.atom_1->residue;
         mmdb::Residue *r2 = o.atom_2->residue;
         std::pair<mmdb::Residue *, mmdb::Residue *>  p(r1, r2);
         if (! matches_existing_pair(p, overlapping_residues))
            overlapping_residues.push_back(p);
      }
   }

   bool things_have_changed = false;
   for (unsigned int ii=0; ii<overlapping_residues.size(); ii++) {
      mmdb::Residue *residue_1_p = overlapping_residues[ii].first;
      mmdb::Residue *residue_2_p = overlapping_residues[ii].second;

      delete_sidechain(residue_1_p);
      delete_sidechain(residue_2_p);

      fill_partial_residue(coot::residue_spec_t(residue_1_p), geom_p, imol_refinement_map);
      fill_partial_residue(coot::residue_spec_t(residue_2_p), geom_p, imol_refinement_map);
      things_have_changed = true;
   }

   if (things_have_changed) {
      have_unsaved_changes_flag = 1;
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked(__FUNCTION__);
   }
}

// carbohydrate building - WTA for the moment
void
molecule_class_info_t::add_named_glyco_tree(const std::string &glycosylation_type,
                                            coot::protein_geometry *geom_p,
                                            const coot::residue_spec_t &asn_res_spec,
                                            const clipper::Xmap<float> &xmap) {

   float mt = coot::util::median_temperature_factor(atom_sel.atom_selection, atom_sel.n_selected_atoms, 2.0, 2222.2, false, false);
   float b_factor_for_new_atoms = 1.55 * mt;
   make_backup();
   coot::cho::add_named_glyco_tree(glycosylation_type, &atom_sel, imol_no, b_factor_for_new_atoms, xmap,
                                   geom_p, asn_res_spec.chain_id, asn_res_spec.res_no);
   have_unsaved_changes_flag = true;
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   make_bonds_type_checked(__FUNCTION__);

}

