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

// include files needed to include molecule-class-info.h correctly. Useful.
#include "Cartesian.h"
#include "mmdb_manager.h" 
#include "mmdb-extras.h"
#include "mmdb-crystal.h"

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
	    if (ele == " H") {

	       const char *chain_id  = new_at->GetChainID();
	       const char *atom_name = new_at->GetAtomName();
	       int resno             = new_at->GetSeqNum();
	       const char *ins_code  = new_at->GetInsCode();

	       
	       int selHnd = atom_sel.mol->NewSelection();
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
		  std::cout << "Ooops in add_hydrogens_from_file() " << std::endl;
	       } else {

		  // normal case
		  CAtom *at_copy = new CAtom;
		  at_copy->Copy(new_at);
		  SelResidues[0]->AddAtom(at_copy);
		  added = 1;
	       }
	       atom_sel.mol->DeleteSelection(selHnd);
	    } 
	 }
      }
   }
   if (added) { 
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
   }
} 
