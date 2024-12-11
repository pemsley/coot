
#include <iostream>
#include <vector>
#include <string>

#include "coot-utils/atom-selection-container.hh"
#include "ghost-molecule-display.hh"
#include "ncs-ghosts.hh"

// I copied this function over without reading it - I now think that this might not be the function that
// I want.

// fill ncs ghosts:  detect exact ncs automatically.
//
std::vector<coot::ghost_molecule_display_t>
fill_ghost_info(atom_selection_container_t *atom_sel, float homology_lev, bool is_from_shelx_ins_flag) {

   if (true)
      std::cout << "DEBUG::   --------------- in fill_ghost_info ------- with homology_lev "
                << homology_lev << std::endl;

   bool do_rtops_flag = true;
   std::vector<std::string> chain_ids;
   std::vector<std::vector<std::pair<std::string, int> > > residue_types;
   std::vector<int> chain_atom_selection_handles;
   std::vector<short int> first_chain_of_this_type;

   bool allow_offset_flag = 0;
   if (is_from_shelx_ins_flag) allow_offset_flag = true;

   std::vector<coot::ghost_molecule_display_t> ncs_ghosts;

   if (atom_sel->n_selected_atoms > 0) {

      int n_models = atom_sel->mol->GetNumberOfModels();
      if (n_models > 0) {
	 int imod = 1; // otherwise madness

	 mmdb::Model *model_p = atom_sel->mol->GetModel(imod);
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();

	 // std::cout << "DEUBG:: in fill_ghost_info() nchains is " << nchains << std::endl;

	 if (nchains > 0) {
	    chain_ids.resize(nchains);
	    residue_types.resize(nchains);
	    chain_atom_selection_handles.resize(nchains);
	    first_chain_of_this_type.resize(nchains, 1);
	    for (int ichain=0; ichain<nchains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
	       if (! chain_p->isSolventChain()) {
		  chain_ids[ichain] = chain_p->GetChainID();
		  int iselhnd = atom_sel->mol->NewSelection();
		  mmdb::PAtom *atom_selection = NULL;
		  int nSelAtoms;
		  atom_sel->mol->SelectAtoms(iselhnd, imod,
					    chain_p->GetChainID(),
					    mmdb::ANY_RES, "*",
					    mmdb::ANY_RES, "*",
					    "*", "*", "*", "*");
		  atom_sel->mol->GetSelIndex(iselhnd, atom_selection, nSelAtoms);
		  chain_atom_selection_handles[ichain] = iselhnd;

		  // debugging the atom selection
		  if (false) {
		     mmdb::PPAtom selatoms_1 = NULL;
		     int n_sel_atoms_1 = 0;
		     atom_sel->mol->GetSelIndex(iselhnd, selatoms_1, n_sel_atoms_1);
 		     std::cout << "DEBUG:: fill_ghost_info: first atom of " << n_sel_atoms_1
 			       << " in " << chain_p->GetChainID() << "  " << iselhnd
 			       << "  selection " << selatoms_1[0] << std::endl;
		  }

		  int nres = chain_p->GetNumberOfResidues();
		  residue_types[ichain].resize(nres);
// 		  std::cout << "INFO:: residues_types[" << ichain << "] resized to "
// 			    << residue_types[ichain].size() << std::endl;
		  for (int ires=0; ires<nres; ires++) {
                     mmdb::PResidue residue_p = chain_p->GetResidue(ires);
		     std::string resname(residue_p->name);
		     residue_types[ichain][ires] = std::pair<std::string, int> (resname, residue_p->seqNum);
		  }
	       }
	    }
	 }
      }

      // add_ncs_ghosts_no_explicit_master(chain_ids, residue_types, first_chain_of_this_type,
      //   				chain_atom_selection_handles, do_rtops_flag, homology_lev,
      //   				allow_offset_flag);

   }
   return ncs_ghosts;
}
