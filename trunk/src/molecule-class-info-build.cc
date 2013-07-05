/* src/molecule-class-info-build.cc
 * 
 * Copyright 2011 The University of Oxford
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

#include <algorithm>

#include <string.h>

// include files needed to include molecule-class-info.h correctly. Useful.
#include <mmdb/mmdb_manager.h> 
#include "coords/Cartesian.h"
#include "coords/mmdb-extras.h"
#include "coords/mmdb-crystal.h"

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
      chain_id_map[it->chain]++;

   if (chain_id_map.size() != 1) {
      std::cout << "WARNING:: all residues need to be in the same chain. Aborted loop selection"
		<< std::endl;
   } else {

      std::sort(local_specs.begin(), local_specs.end());

      int i_loop_res=0;
      int n_loop_residues = local_specs.back().resno - local_specs.begin()->resno + 1;
      
      for (int i_loop_res=0; i_loop_res<n_loop_residues; i_loop_res++) {

	 coot::residue_spec_t test_spec(*local_specs.begin());
	 test_spec.resno += i_loop_res;

	 if (std::find(local_specs.begin(), local_specs.end(), test_spec) == local_specs.end()) {
	    // Add a null residue
	    std::cout << "Added a null for " << test_spec << std::endl;
	    ProteinDB::Residue residue;
	    chain.add_residue(residue);
	 } else { 
	    CResidue *residue_p = get_residue(test_spec);
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

	       CAtom *at_n  = residue_p->GetAtom(" N  ");
	       CAtom *at_ca = residue_p->GetAtom(" CA ");
	       CAtom *at_c  = residue_p->GetAtom(" C  ");
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

   make_backup();
   bool added = 0;
   atom_selection_container_t asc = get_atom_selection(reduce_pdb_out, 0);
   if (asc.read_success) { 
      int imod = 1;
      CModel *new_model_p = asc.mol->GetModel(imod);
      CChain *new_chain_p;
      int n_new_chains = new_model_p->GetNumberOfChains();
      for (int i_new_chain=0; i_new_chain<n_new_chains; i_new_chain++) {
	 new_chain_p = new_model_p->GetChain(i_new_chain);
	 int n_new_res = new_chain_p->GetNumberOfResidues();
	 CResidue *new_residue_p;
	 CAtom *new_at;
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
		  PPCResidue SelResidues;
		  atom_sel.mol->Select(selHnd, STYPE_RESIDUE, 1,
				       chain_id, 
				       resno, ins_code,
				       resno, ins_code,
				       "*",  // residue name
				       "*",  // Residue must contain this atom name?
				       "*",  // Residue must contain this Element?
				       "*",  // altLocs
				       SKEY_NEW // selection key
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
		     CAtom *at_in_residue = SelResidues[0]->GetAtom(atom_name);
		     
		     if (at_in_residue) {
			at_in_residue->x = new_at->x;
			at_in_residue->y = new_at->y;
			at_in_residue->z = new_at->z;
		     } else { 
		     
			CAtom *at_copy = new CAtom;
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

// --------- LINKs ---------------
void
molecule_class_info_t::make_link(const coot::atom_spec_t &spec_1, const coot::atom_spec_t &spec_2,
				 const std::string &link_name, float length,
				 const coot::protein_geometry &geom) {

   // I don't see how link_name and length are part of a CLink.
   // Perhaps they should not be passed then?
   
   CAtom *at_1 = get_atom(spec_1);
   CAtom *at_2 = get_atom(spec_2);

   if (! at_1) {
      std::cout << "WARNING:: atom " << spec_1 << " not found - abandoning LINK addition " << std::endl;
   } else { 
      if (! at_2) {
	 std::cout << "WARNING:: atom " << spec_1 << " not found - abandoning LINK addition " << std::endl;
      } else {

	 CModel *model_1 = at_1->GetModel();
	 CModel *model_2 = at_1->GetModel();

	 if (model_1 != model_2) {

	    std::cout << "WARNING:: specified atoms have mismatching models - abandoning LINK addition"
		      << std::endl;

	 } else {

	    make_backup();

	    CMMDBManager *mol = atom_sel.mol;
	 
	    CLink *link = new CLink; // sym ids default to 1555 1555

	    strncpy(link->atName1,  at_1->GetAtomName(), 19);
	    strncpy(link->aloc1,    at_1->altLoc, 9);
	    strncpy(link->resName1, at_1->GetResName(), 19);
	    strncpy(link->chainID1, at_1->GetChainID(), 9);
	    strncpy(link->insCode1, at_1->GetInsCode(), 9);
	    link->seqNum1         = at_1->GetSeqNum();

	    strncpy(link->atName2,  at_2->GetAtomName(), 19);
	    strncpy(link->aloc2,    at_2->altLoc, 9);
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
	    std::vector<std::pair<bool, CResidue *> > residues(2);
	    residues[0] = std::pair<bool, CResidue *> (0, at_1->residue);
	    residues[1] = std::pair<bool, CResidue *> (0, at_2->residue);
	    std::vector<coot::atom_spec_t> dummy_fixed_atom_specs;
	    coot::restraints_container_t restraints(residues, geom, mol, dummy_fixed_atom_specs);
	    coot::bonded_pair_container_t bpc = restraints.bonded_residues_from_res_vec(geom);
	    bpc.apply_chem_mods(geom);
	    atom_sel.mol->FinishStructEdit();
	    
	    update_molecule_after_additions();
	 }
      } 
   }
}

void
molecule_class_info_t::delete_any_link_containing_residue(const coot::residue_spec_t &res_spec) {

   if (atom_sel.mol) {
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (unsigned int imod=1; imod<=n_models; imod++) {
	 CModel *model_p = atom_sel.mol->GetModel(imod);
	 if ((res_spec.model_number == imod) || (res_spec.model_number == MinInt4)) {
	    CLinkContainer *links = model_p->GetLinks();
	    int n_links = model_p->GetNumberOfLinks();
	    for (unsigned int ilink=1; ilink<=n_links; ilink++) { 
	       CLink *link_p = model_p->GetLink(ilink);
	       // CLink *link = static_cast<CLink *>(links->GetContainerClass(ilink));

	       if (link_p) { 

		  // must pass a valid link_p
		  std::pair<coot::atom_spec_t, coot::atom_spec_t> link_atoms = coot::link_atoms(link_p);
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
molecule_class_info_t::delete_link(CLink *link, CModel *model_p) {

   // Copy out the links, delete all links and add the saved links back
   
   std::vector<CLink *> saved_links;
   int n_links = model_p->GetNumberOfLinks();
   for (unsigned int ilink=1; ilink<=n_links; ilink++) {
      CLink *model_link = model_p->GetLink(ilink);
      if (model_link != link) { 
	 CLink *copy_link = new CLink(*model_link);
	 saved_links.push_back(copy_link);
      }
   }
   model_p->RemoveLinks();
   for (unsigned int i=0; i<saved_links.size(); i++) { 
      model_p->AddLink(saved_links[i]);
   }
}


void
molecule_class_info_t::move_reference_chain_to_symm_chain_position(coot::Symm_Atom_Pick_Info_t naii) {

   if (naii.success) {

      make_backup();

      mat44 my_matt;
      mat44 pre_shift_matt;
      
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
	    CChain *moving_chain = atom_sel.atom_selection[naii.atom_index]->residue->chain;
	    for (unsigned int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
	       CAtom *at = atom_sel.atom_selection[iat];
	       if (at->residue->chain == moving_chain) { 
		  at->Transform(pre_shift_matt);
		  at->Transform(my_matt);
	       }
	    }


	    have_unsaved_changes_flag = 1; // because we do a backup whatever...
	    atom_sel.mol->FinishStructEdit();
	    update_molecule_after_additions();
	    update_symmetry();
	 }
      }
   } 
}
