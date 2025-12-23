/* src/molecule-class-info-build.cc
 * 
 * Copyright 2011 The University of Oxford
 * Copyright 2013 by Medical Research Council
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

#include <algorithm>

#include <string.h>

#include "compat/coot-sysdep.h"

// include files needed to include molecule-class-info.h correctly. Useful.
#include <mmdb2/mmdb_manager.h>

#include "coords/Cartesian.hh"
#include "coords/mmdb-extras.hh"
#include "coords/mmdb-crystal.hh"

#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/reduce.hh"
#include "coot-utils/reduce.hh"

#include "molecule-class-info.h"


std::vector<ProteinDB::Chain>
molecule_class_info_t::protein_db_loops(const std::vector<coot::residue_spec_t> &residue_specs,
					int nfrags,
					const clipper::Xmap<float> &xmap) {

   std::string pkg_data_dir = coot::package_data_dir();
   std::string dir = coot::util::append_dir_dir(pkg_data_dir, "protein_db");
   std::string file_name = coot::util::append_dir_file(dir, "protein.db");

   std::vector<clipper::Coord_orth> clash_coords;

   ProteinDB::Chain chain = make_fragment_chain(residue_specs);

   ProteinDB::ProteinDBSearch protein_db_search(file_name);
   std::vector<ProteinDB::Chain> chains = protein_db_search.search(chain, nfrags, xmap, clash_coords);

   return chains;
}


ProteinDB::Chain
molecule_class_info_t::make_fragment_chain(const std::vector<coot::residue_spec_t> &residue_specs) const {

   ProteinDB::Chain chain;

   std::vector<coot::residue_spec_t> local_specs = residue_specs;
   std::map<std::string, int> chain_id_map;

   std::vector<coot::residue_spec_t>::const_iterator it;
   for (it=local_specs.begin(); it!= local_specs.end(); it++)
      chain_id_map[it->chain_id]++;

   if (chain_id_map.size() != 1) {
      std::cout << "WARNING:: all residues need to be in the same chain. Aborted loop selection"
		<< std::endl;
   } else {

      std::sort(local_specs.begin(), local_specs.end());

      int i_loop_res=0;
      int n_loop_residues = local_specs.back().res_no - local_specs.begin()->res_no + 1;
      
      for (int i_loop_res=0; i_loop_res<n_loop_residues; i_loop_res++) {

	 coot::residue_spec_t test_spec(*local_specs.begin());
	 test_spec.res_no += i_loop_res;

	 if (std::find(local_specs.begin(), local_specs.end(), test_spec) == local_specs.end()) {
	    // Add a null residue
	    std::cout << "Added a null for " << test_spec << std::endl;
	    ProteinDB::Residue residue;
	    chain.add_residue(residue);
	 } else { 
	    mmdb::Residue *residue_p = get_residue(test_spec);
	    if (! residue_p) {
	       std::cout << "oops - missing residue " << test_spec << std::endl;
	       std::cout << "Added a null for " << test_spec << std::endl;
	       ProteinDB::Residue residue;
	       chain.add_residue(residue);
	    } else {
	       std::string type = residue_p->GetResName();
	       clipper::Coord_orth  n_pos;
	       clipper::Coord_orth ca_pos;
	       clipper::Coord_orth  c_pos;

	       mmdb::Atom *at_n  = residue_p->GetAtom(" N  ");
	       mmdb::Atom *at_ca = residue_p->GetAtom(" CA ");
	       mmdb::Atom *at_c  = residue_p->GetAtom(" C  ");
	       if (at_n)
		  n_pos = clipper::Coord_orth( at_n->x,  at_n->y,  at_n->z);
	       if (at_ca)
		  ca_pos = clipper::Coord_orth(at_ca->x, at_ca->y, at_ca->z);
	       if (at_c)
		  n_pos = clipper::Coord_orth( at_c->x,  at_c->y,  at_c->z);


	       if (at_n && at_ca && at_c) { 
		  ProteinDB::Residue residue(n_pos, ca_pos, c_pos, type);
		  chain.add_residue(residue);
		  // std::cout << "Added a real residue for " << test_spec << std::endl;
	       }
	    }
	 }
      }
   }
   return chain;
} 

void
molecule_class_info_t::add_hydrogens_from_file(const std::string &reduce_pdb_out) {

   std::cout << "adding hydrogens from PDB file " << reduce_pdb_out << std::endl;

   make_backup(__FUNCTION__);
   bool added = 0;
   atom_selection_container_t asc = get_atom_selection(reduce_pdb_out, false, true, true);
   if (asc.read_success) { 
      int imod = 1;
      mmdb::Model *new_model_p = asc.mol->GetModel(imod);
      mmdb::Chain *new_chain_p;
      int n_new_chains = new_model_p->GetNumberOfChains();
      for (int i_new_chain=0; i_new_chain<n_new_chains; i_new_chain++) {
	 new_chain_p = new_model_p->GetChain(i_new_chain);
	 int n_new_res = new_chain_p->GetNumberOfResidues();
	 mmdb::Residue *new_residue_p;
	 mmdb::Atom *new_at;
	 for (int i_new_res=0; i_new_res<n_new_res; i_new_res++) { 
	    new_residue_p = new_chain_p->GetResidue(i_new_res);
	    int n_new_atoms = new_residue_p->GetNumberOfAtoms();
	    for (int i_new_at=0; i_new_at<n_new_atoms; i_new_at++) {
	       new_at = new_residue_p->GetAtom(i_new_at);
	       std::string ele = new_at->element;
	       if (ele == " H" || ele == "H") { // PDB v2 and v3.

		  const char *chain_id  = new_at->GetChainID();
		  const char *atom_name = new_at->GetAtomName();
		  int resno             = new_at->GetSeqNum();
		  const char *ins_code  = new_at->GetInsCode();

	       
		  int selHnd = atom_sel.mol->NewSelection(); // d 
		  int nSelResidues;
		  mmdb::PPResidue SelResidues;
		  atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 1,
				       chain_id, 
				       resno, ins_code,
				       resno, ins_code,
				       "*",  // residue name
				       "*",  // Residue must contain this atom name?
				       "*",  // Residue must contain this Element?
				       "*",  // altLocs
				       mmdb::SKEY_NEW // selection key
				       );
		  atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
	       
		  if (nSelResidues != 1) {
		     std::cout << "Ooops in add_hydrogens_from_file() - expected 1 residue but found "
			       << nSelResidues << " residues with attributes \"" << chain_id << "\" "
			       << resno << " \"" << ins_code << "\""
			       << std::endl;
		  } else {

		     // normal case

		     // if the atom exists, update the coordinates, else, add an atom.
		     mmdb::Atom *at_in_residue = SelResidues[0]->GetAtom(atom_name);
		     
		     if (at_in_residue) {
			at_in_residue->x = new_at->x;
			at_in_residue->y = new_at->y;
			at_in_residue->z = new_at->z;
		     } else { 
		     
			mmdb::Atom *at_copy = new mmdb::Atom;
			at_copy->Copy(new_at);
			SelResidues[0]->AddAtom(at_copy);
			added = 1;
		     }
		  }
		  atom_sel.mol->DeleteSelection(selHnd);
	       } 
	    }
	 }
      }
   }
   have_unsaved_changes_flag = 1; // because we do a backup whatever...
   if (added) {
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
   }
}

void
molecule_class_info_t::add_hydrogen_atoms_to_residue(const coot::residue_spec_t &rs) {

   make_backup(__FUNCTION__);
   mmdb::Residue *residue_this_p = get_residue(rs);
   mmdb::Residue *residue_prev_p = coot::util::get_previous_residue(rs, atom_sel.mol);

   bool go_nuclear = false;
   coot::reduce r(atom_sel.mol, imol_no);
   r.add_hydrogens_to_residue(residue_this_p, residue_prev_p, go_nuclear);

   have_unsaved_changes_flag = 1; // because we do a backup whatever...
   atom_sel.mol->FinishStructEdit();
   update_molecule_after_additions();
}


// --------- LINKs ---------------
void
molecule_class_info_t::make_link(const coot::atom_spec_t &spec_1, const coot::atom_spec_t &spec_2,
				 const std::string &link_name, float length,
				 const coot::protein_geometry &geom) {

   // 2014: link_name and length are not part curently of a mmdb::Link.
   // Perhaps they should not be passed then?

   mmdb::Atom *at_1 = get_atom(spec_1);
   mmdb::Atom *at_2 = get_atom(spec_2);

   if (! at_1) {
      std::cout << "WARNING:: atom " << spec_1 << " not found - abandoning LINK addition " << std::endl;
   } else {
      if (! at_2) {
	 std::cout << "WARNING:: atom " << spec_1 << " not found - abandoning LINK addition " << std::endl;
      } else {

	 mmdb::Model *model_1 = at_1->GetModel();
	 mmdb::Model *model_2 = at_1->GetModel();

	 if (model_1 != model_2) {

	    std::cout << "WARNING:: specified atoms have mismatching models - abandoning LINK addition"
		      << std::endl;

	 } else {

	    make_backup(__FUNCTION__);

	    mmdb::Manager *mol = atom_sel.mol;

            mmdb::Link *link = new mmdb::Link; // sym ids default to 1555 1555

	    strncpy(link->atName1,  at_1->GetAtomName(), 20);
	    strncpy(link->aloc1,    at_1->altLoc, 20);
	    strncpy(link->resName1, at_1->GetResName(), 19);
	    strncpy(link->chainID1, at_1->GetChainID(), 9);
	    strncpy(link->insCode1, at_1->GetInsCode(), 9);
	    link->seqNum1         = at_1->GetSeqNum();

	    strncpy(link->atName2,  at_2->GetAtomName(), 20);
	    strncpy(link->aloc2,    at_2->altLoc, 20);
	    strncpy(link->resName2, at_2->GetResName(), 19);
	    strncpy(link->chainID2, at_2->GetChainID(), 9);
	    strncpy(link->insCode2, at_2->GetInsCode(), 9);
	    link->seqNum2         = at_2->GetSeqNum();

	    model_1->AddLink(link);
	    have_unsaved_changes_flag = 1;
	    atom_sel.mol->FinishStructEdit();

	    // now, do we need to do any deletions to the model that
	    // are defined in the dictionary?
	    //
	    std::vector<std::pair<bool, mmdb::Residue *> > residues(2);
	    residues[0] = std::pair<bool, mmdb::Residue *> (0, at_1->residue);
	    residues[1] = std::pair<bool, mmdb::Residue *> (0, at_2->residue);
	    std::vector<coot::atom_spec_t> dummy_fixed_atom_specs;

	    // convert to restraints_container_t interface
	    mmdb::Link local_link;
	    local_link.Copy(link);
	    std::vector<mmdb::Link> links(1);
	    links[0] = local_link;
	    clipper::Xmap<float> dummy_xmap;

	    coot::restraints_container_t restraints(residues,
						    links,
						    geom, mol, dummy_fixed_atom_specs, &dummy_xmap);
	    coot::bonded_pair_container_t bpc = restraints.bonded_residues_from_res_vec(geom);

	    asn_hydrogen_position_swap(residues); // HD21 and HD22 (HD22 will be deleted)
	    bpc.apply_chem_mods(geom);
	    atom_sel.mol->FinishStructEdit();

	    update_molecule_after_additions();
	 }
      }
   }
}

void
molecule_class_info_t::update_any_link_containing_residue(const coot::residue_spec_t &old_spec,
							  const coot::residue_spec_t &new_spec) {

   if (atom_sel.mol) {
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 if ((old_spec.model_number == imod) || (old_spec.model_number == mmdb::MinInt4)) {
	    mmdb::LinkContainer *links = model_p->GetLinks();
	    int n_links = model_p->GetNumberOfLinks();
	    for (int ilink=1; ilink<=n_links; ilink++) {
	       mmdb::Link *link_p = model_p->GetLink(ilink);
	       // mmdb::Link *link = static_cast<mmdb::Link *>(links->GetContainerClass(ilink));

	       if (link_p) {

		  // must pass a valid link_p
		  std::pair<coot::atom_spec_t, coot::atom_spec_t> link_atoms = coot::link_atoms(link_p, model_p);
		  coot::residue_spec_t res_1(link_atoms.first);
		  coot::residue_spec_t res_2(link_atoms.second);

		  // the atom names and the alt confs will stay the same
		  // (this function is invoked after renumber residues).

		  if (res_1 == old_spec) {
		     if (res_1.chain_id != new_spec.chain_id)
			strncpy(link_p->chainID1, new_spec.chain_id.c_str(), 9);
		     if (res_1.res_no != new_spec.res_no)
			link_p->seqNum1 = new_spec.res_no;
		     mmdb::Residue *new_res = get_residue(new_spec);
		     if (new_res) {
			std::string res_name_new(new_res->GetResName());
			std::string current_name_1(link_p->resName1);
			if (res_name_new != current_name_1) {
			   strncpy(link_p->resName1, res_name_new.c_str(), 19);
			}
		     }
		  }
		  if (res_2 == old_spec) {
		     if (res_2.chain_id != new_spec.chain_id)
			strncpy(link_p->chainID2, new_spec.chain_id.c_str(), 9);
		     if (res_2.res_no != new_spec.res_no)
			link_p->seqNum2 = new_spec.res_no;
		     mmdb::Residue *new_res = get_residue(new_spec);
		     if (new_res) {
			std::string res_name_new(new_res->GetResName());
			std::string current_name_2(link_p->resName2);
			if (res_name_new != current_name_2) {
			   strncpy(link_p->resName1, res_name_new.c_str(), 19);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
}


void
molecule_class_info_t::delete_any_link_containing_residue(const coot::residue_spec_t &res_spec) {

   if (atom_sel.mol) {
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 if ((res_spec.model_number == imod) || (res_spec.model_number == mmdb::MinInt4)) {
	    mmdb::LinkContainer *links = model_p->GetLinks();
	    int n_links = model_p->GetNumberOfLinks();
	    for (int ilink=1; ilink<=n_links; ilink++) {
	       mmdb::Link *link_p = model_p->GetLink(ilink);
	       // mmdb::Link *link = static_cast<mmdb::Link *>(links->GetContainerClass(ilink));

	       if (link_p) {

		  // must pass a valid link_p
		  std::pair<coot::atom_spec_t, coot::atom_spec_t> link_atoms = coot::link_atoms(link_p, model_p);
		  coot::residue_spec_t res_1(link_atoms.first);
		  coot::residue_spec_t res_2(link_atoms.second);
		  // std::cout << "found link " << res_1 << " to "  << res_2 << std::endl;
		  if (res_spec == res_1) { 
		     delete_link(link_p, model_p);
		  } 
		  if (res_spec == res_2) { 
		     delete_link(link_p, model_p);
		  }
	       } else {
		  std::cout << "ERROR:: Null link_p for link " << ilink << " of " << n_links << std::endl;
	       } 
	    }
	 }
      }
   }
}

void
molecule_class_info_t::delete_link(mmdb::Link *link, mmdb::Model *model_p) {

   // Copy out the links, delete all links and add the saved links back

   std::vector<mmdb::Link *> saved_links;
   int n_links = model_p->GetNumberOfLinks();
   for (int ilink=1; ilink<=n_links; ilink++) {
      mmdb::Link *model_link = model_p->GetLink(ilink);
      if (model_link != link) { 
	 mmdb::Link *copy_link = new mmdb::Link(*model_link);
	 saved_links.push_back(copy_link);
      }
   }

   model_p->RemoveLinks();
   for (unsigned int i=0; i<saved_links.size(); i++) { 
      model_p->AddLink(saved_links[i]);
   }
}

// there is a copy of this in add-linked-cho now.
void
molecule_class_info_t::asn_hydrogen_position_swap(std::vector<std::pair<bool, mmdb::Residue *> > residues) {

   if (residues[0].second) {
      if (residues[1].second) {
	 std::string rn0(residues[0].second->GetResName());
	 std::string rn1(residues[1].second->GetResName());
	 mmdb::Residue *r_0 = 0;
	 mmdb::Residue *r_1 = 0;
	 if (rn0 == "ASN") {
	    if (rn1 == "NAG") {
	       r_0 = residues[0].second;
	       r_1 = residues[1].second;
	    }
	 }
	 if (rn1 == "ASN") {
	    if (rn0 == "NAG") {
	       r_1 = residues[0].second;
	       r_0 = residues[1].second;
	    }
	 }

	 if (r_1 && r_0) {
	    mmdb::Atom *at_hd21 = 0;
	    mmdb::Atom *at_hd22 = 0;
	    mmdb::Atom **residue_atoms_0 = 0;
	    int n_residue_atoms_0;
	    r_0->GetAtomTable(residue_atoms_0, n_residue_atoms_0);
	    for (int iat=0; iat<n_residue_atoms_0; iat++) {
	       mmdb::Atom *at = residue_atoms_0[iat];
	       std::string atom_name(at->GetAtomName());
	       if (atom_name == "HD21") at_hd21 = at;
	       if (atom_name == "HD22") at_hd22 = at;
	    }
	    if (at_hd21 && at_hd22) {
	       clipper::Coord_orth co21 = coot::co(at_hd21);
	       clipper::Coord_orth co22 = coot::co(at_hd22);
	       at_hd21->x = co22.x();
	       at_hd21->y = co22.y();
	       at_hd21->z = co22.z();
	       at_hd22->x = co21.x(); // this atom will be deleted.
	       at_hd22->y = co21.y();
	       at_hd22->z = co21.z();
	    }
	 }
      }
   }
}


void
molecule_class_info_t::move_reference_chain_to_symm_chain_position(coot::Symm_Atom_Pick_Info_t naii) {

   if (naii.success) {

      make_backup(__FUNCTION__);

      mmdb::mat44 my_matt;
      mmdb::mat44 pre_shift_matt;
      
      int err_1 = atom_sel.mol->GetTMatrix(my_matt,
					   naii.symm_trans.isym(),
					   naii.symm_trans.x(),
					   naii.symm_trans.y(),
					   naii.symm_trans.z());
      
      int err_2 = atom_sel.mol->GetTMatrix(pre_shift_matt, 0,
					   -naii.pre_shift_to_origin.us,
					   -naii.pre_shift_to_origin.vs,
					   -naii.pre_shift_to_origin.ws);

      if (err_1 == 0) { 
	 if (err_2 == 0) {

	    // This used to work - but no longer if pre_shift_matt and my_matt
	    // are not references - strange
	    
	    mmdb::Chain *moving_chain = atom_sel.atom_selection[naii.atom_index]->residue->chain;
	    for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
	       mmdb::Atom *at = atom_sel.atom_selection[iat];
	       if (at->residue->chain == moving_chain) { 
		  // at->Transform(&pre_shift_matt);
		  // at->Transform(&my_matt);
	       }
	    }

	    coot::util::transform_chain(atom_sel.mol, moving_chain,
					atom_sel.n_selected_atoms, atom_sel.atom_selection,
					pre_shift_matt);
	    coot::util::transform_chain(atom_sel.mol, moving_chain,
					atom_sel.n_selected_atoms, atom_sel.atom_selection,
					my_matt);


	    have_unsaved_changes_flag = 1; // because we do a backup whatever...
	    atom_sel.mol->FinishStructEdit();
	    update_molecule_after_additions();
	    // update NCS ghosts, if they were present
	    if (ncs_ghosts.size() > 0) {
	       fill_ghost_info(true, 0.7);
	    }
	    update_symmetry();
	 }
      }
   }
}

#include "cootaneer/buccaneer-prot.h"

void
molecule_class_info_t::globularize() {


   mmdb::Manager *mol = atom_sel.mol;

   if (mol) { 
      make_backup(__FUNCTION__);

      bool nucleotides = false;

      // now check if we have nucleotides
      std::pair<unsigned int, unsigned int> number_of_typed_residues(0,0);
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (chain_p) {
               std::pair<unsigned int, unsigned int> nric = coot::util::get_number_of_protein_or_nucleotides(chain_p);
               number_of_typed_residues.first  = nric.first;
               number_of_typed_residues.second = nric.second;
            }
         }
      }
      if (number_of_typed_residues.second > number_of_typed_residues.first)
         nucleotides = true;

      clipper::MiniMol mm;
      clipper::MMDBfile* mmdbfile = static_cast<clipper::MMDBfile*>(mol);
      mmdbfile->import_minimol(mm);
      bool r = ProteinTools::globularise(mm, nucleotides);

      mmdbfile->export_minimol(mm);

      have_unsaved_changes_flag = 1; // because we do a backup whatever...
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
      update_symmetry();

   }

} 


void
molecule_class_info_t::reduce(coot::protein_geometry *geom_p) {

   bool go_nuclear = false; // pass this, I think

   make_backup(__FUNCTION__);
   mmdb::Manager *mol = atom_sel.mol;
   coot::reduce r(mol, imol_no);
   r.add_geometry(geom_p);
   r.add_hydrogen_atoms(go_nuclear);
   update_molecule_after_additions();

   if (false) { // debug
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 mmdb::Atom *at = atom_sel.atom_selection[i];
	 int idx;
	 if (at->GetUDData(atom_sel.UDDAtomIndexHandle, idx) == mmdb::UDDATA_Ok) {
	    std::string fail;
	    if (idx != i)
	       fail = "FAIL";
	    std::cout << "atom " << coot::atom_spec_t(at) << " has index " << i << " and udd " << idx << " " << fail << "\n";
	 } else {
	    std::cout << "bad GetUDData() for atom " << coot::atom_spec_t(at) << std::endl;
	 }
      }
   }

   update_symmetry();
}

void
molecule_class_info_t::switch_HIS_protonation(coot::residue_spec_t res_spec) {

   mmdb::Residue *residue_p = get_residue(res_spec);
   if (residue_p) {

      mmdb::Atom *at = 0;
      mmdb::Atom *atom_1 = 0;
      mmdb::Atom *atom_2 = 0;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 mmdb::Atom *res_at = residue_atoms[i];
	 std::string atom_name = res_at->name;
	 if (atom_name == " HD1")
	    atom_1 = res_at;
	 if (atom_name == " HE2")
	    atom_2 = res_at;
      }
      
      if (atom_1 && ! atom_2) at = atom_1;
      if (atom_2 && ! atom_1) at = atom_2;

      if (at) {
	 make_backup(__FUNCTION__);
	 coot::reduce r(atom_sel.mol, imol_no);
	 r.switch_his_protonation(residue_p, at);
	 update_molecule_after_additions();
	 update_symmetry();
      }
   }
}



#include "ideal/crankshaft.hh"

void
molecule_class_info_t::crankshaft_peptide_rotation_optimization(const coot::residue_spec_t &rs,
								unsigned int n_peptides,
								const clipper::Xmap<float> &xmap,
								float map_weight,
								int n_samples,
								ctpl::thread_pool *thread_pool_p,
								int n_threads) {
   if (n_threads < 1) n_threads = 1;
   int n_solutions = 1;

   std::vector<mmdb::Manager *> mols =
      coot::crankshaft::crank_refine_and_score(rs, n_peptides, xmap, atom_sel.mol, map_weight,
					       n_samples, n_solutions, thread_pool_p, n_threads);

   if (mols.size() == 1) {
      make_backup(__FUNCTION__);
      std::cout << "DEBUG:: crankshaft updated " << std::endl;
      // what do we do with the old atom selection and mol?
      // should we delete them here?
      //
      atom_sel = make_asc(mols[0]);
      have_unsaved_changes_flag = 1;
      update_molecule_after_additions();
      update_symmetry();
   }
}

int
molecule_class_info_t::add_residue_with_atoms(const coot::residue_spec_t &residue_spec, const std::string &res_name, const std::vector<coot::minimol::atom> &list_of_atoms) {

   std::cout << "start add_residue_with_atoms()" << std::endl;
   int r_added = 0;

   mmdb::Residue *residue_p = get_residue(residue_spec);
   if (! residue_p) {
      // make one then
      mmdb::Chain *chain_p_for_residue = 0; // set this
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id(chain_p->GetChainID());
            if (chain_id == residue_spec.chain_id)
               chain_p_for_residue = chain_p;
         }
         if (chain_p_for_residue) {
            residue_p = new mmdb::Residue(chain_p_for_residue, res_name.c_str(), residue_spec.res_no, residue_spec.ins_code.c_str());
         } else {
            chain_p_for_residue = new mmdb::Chain;
            chain_p_for_residue->SetChain(residue_spec.chain_id.c_str());
            model_p->AddChain(chain_p_for_residue);
            residue_p = new mmdb::Residue(chain_p_for_residue, res_name.c_str(), residue_spec.res_no, residue_spec.ins_code.c_str());
         }
      } else {
         std::cout << "ERROR:: in add_residue_with_atoms() null model_p " << imod << std::endl;
      }
   }
   if (residue_p) {
      for (unsigned int i=0; i<list_of_atoms.size(); i++) {
         const coot::minimol::atom &atom = list_of_atoms[i];
         mmdb::Atom *at = atom.make_atom();
         if (at) {
            residue_p->AddAtom(at);
         }
      }
      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle); // correct?
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      update_molecule_after_additions();
      update_symmetry();
   }
   return r_added;
}



int
molecule_class_info_t::trim_molecule_by_b_factor(float limit, bool keep_higher_flag) {

   // delete residues if the B-factor is higher (or lower)

   int status = 1;
   make_backup(__FUNCTION__);
   mmdb::Manager *mol = atom_sel.mol;
   std::set<mmdb::Residue *> deletable_residues;
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        float b_factor = at->tempFactor;
                        if (keep_higher_flag) {
                           if (b_factor < limit) {
                              deletable_residues.insert(residue_p);
                           }
                        } else {
                           if (b_factor > limit) {
                              deletable_residues.insert(residue_p);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
#if 0
   if (! deletable_residues.empty()) {
      have_unsaved_changes_flag = true;
      std::set<mmdb::Residue *>::const_iterator it;
      unsigned int n_residues = deletable_residues.size();
      unsigned int ith_residue = 0;
      for (it=deletable_residues.begin(); it!=deletable_residues.end(); ++it) {
         mmdb::Residue *r = *it;
         if (true)
            std::cout << "deleting " << ith_residue << " residue of " << n_residues << " " << r << " " << coot::residue_spec_t(r) << std::endl;
         delete r;
         ith_residue++;
      }
      mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      mol->FinishStructEdit();
      coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
      make_bonds_type_checked();
   }
#endif

   if (! deletable_residues.empty()) {
      std::vector<coot::residue_spec_t> residues;
      std::set<mmdb::Residue *>::const_iterator it;
      for (it=deletable_residues.begin(); it!=deletable_residues.end(); ++it) {
         mmdb::Residue *r = *it;
         coot::residue_spec_t spec(r);
         residues.push_back(spec);
      }
      delete_residues(residues);
   }
   return status;
}


void
molecule_class_info_t::pLDDT_to_b_factor() {

   auto converter_function = [] (float b_in) {
      float b_out = 2.0 * (100.0 - b_in);
      if (b_out < 2.0) b_out = 2.0;
      return b_out;
   };

   float m_b_factor_pre = coot::util::average_temperature_factor(atom_sel.atom_selection, atom_sel.n_selected_atoms, 0, 1000, 0, 0);
   make_backup(__FUNCTION__);
   mmdb::Manager *mol = atom_sel.mol;
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        float b_factor = at->tempFactor;
                        float new_b_factor = converter_function(b_factor);
                        at->tempFactor = new_b_factor;
                        if (true) {
                           std::string atom_name = at->name;
                           if (atom_name == " CA ") {
                              std::cout << "converted b-factor " << b_factor << " " << new_b_factor << std::endl;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   float m_b_factor_post = coot::util::average_temperature_factor(atom_sel.atom_selection, atom_sel.n_selected_atoms, 0, 1000, 0, 0);

   std::cout << "INFO:: average b-factor-pre: " << m_b_factor_pre << " post: " << m_b_factor_post << std::endl;

   have_unsaved_changes_flag = true;
   make_bonds_type_checked(); // we might be looking at B-factor representation

}
