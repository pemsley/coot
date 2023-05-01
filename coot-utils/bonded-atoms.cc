
#include "compat/coot-sysdep.h"
#include "coot-coord-utils.hh"
#include "bonded-atoms.hh"

// called for atoms in different residues - we already know
// that at_2 will be in a different residue to at_1.
//
//
bool
coot::are_polymer_bonded(mmdb::Atom *at_1, mmdb::Atom *at_2) {

   bool state = false;

   std::string at_1_name(at_1->GetAtomName());
   std::string at_2_name(at_2->GetAtomName());

   if (at_1_name == " C  ")
      if (at_2_name == " N  ")
	 state = true;

   if (at_1_name == " O3'")
      if (at_2_name == " P  ")
	 state = true;

   if (! state) {
      if (at_1_name == " O4 " || at_1_name == " O3 " || at_1_name == " O2 ") {
	 if (at_2_name == " C1 ") {
	    // Let's add a distance check here
	    clipper::Coord_orth pt_1 = co(at_1);
	    clipper::Coord_orth pt_2 = co(at_2);
	    float bond_dist_crit = 1.86; // catches SG-CG
	    if ((pt_1-pt_2).lengthsq() < bond_dist_crit*bond_dist_crit) {
	       state = true;
	    }
	 }
      }
   }

   return state;

}

// std::vector<std::set> > ?
std::vector<std::vector<unsigned int> >
coot::make_bonds(mmdb::Manager *mol, int n_selected_atoms, int udd_atom_index_handle) {

   std::vector<std::vector<unsigned int> > v(n_selected_atoms);
   for (int i=0; i<n_selected_atoms; i++)
      v[i].reserve(4);

   bool calc_only = true; // doesn't change anything

   mol->MakeBonds(calc_only); // doesn't do a good job for alt confs, it seems

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {

	 // make a bond tree for this model

	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    if (nres > 1) {
	       for (int ires=0; ires<(nres-1); ires++) {
		  mmdb::Residue *residue_1_p = chain_p->GetResidue(ires);
		  mmdb::Residue *residue_2_p = chain_p->GetResidue(ires+1);
		  int n_atoms_res_1 = residue_1_p->GetNumberOfAtoms();
		  int n_atoms_res_2 = residue_2_p->GetNumberOfAtoms();
		  for (int iat=0; iat<n_atoms_res_1; iat++) {
		     mmdb::Atom *at_i = residue_1_p->GetAtom(iat);
		     for (int jat=0; jat<n_atoms_res_2; jat++) {
			mmdb::Atom *at_j = residue_2_p->GetAtom(jat);
			// std::cout << "iat " << iat << " jat " << jat << std::endl;
			bool b = are_polymer_bonded(at_i, at_j);
			if (b) {
			   int idx_1;
			   int idx_2;
			   at_i->GetUDData(udd_atom_index_handle, idx_1);
			   at_j->GetUDData(udd_atom_index_handle, idx_2);
			   if ((idx_1 < 0) || (idx_1 >= n_selected_atoms)) {
			      std::cout << "atom index problem " << idx_1 << " "
					<< n_selected_atoms << std::endl;
			   } else {
			      if ((idx_2 < 0) || (idx_2 >= n_selected_atoms)) {
				 std::cout << "atom index problem " << idx_2 << " "
					   << n_selected_atoms << std::endl;
			      } else {

				 // Happy path

				 v[idx_1].push_back(idx_2);
				 v[idx_2].push_back(idx_1);
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }

	 // --- make residue internal bonds

	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       int n_atoms = residue_p->GetNumberOfAtoms();
	       for (int iat=0; iat<n_atoms; iat++) {
		  mmdb::Atom *at = residue_p->GetAtom(iat);
		  int nb = at->GetNBonds();
		  mmdb::AtomBond *atom_bonds = 0;
		  int n_atom_bonds = 0;
		  at->GetBonds(atom_bonds, n_atom_bonds);
		  if (false)
		     std::cout << "at " << atom_spec_t(at) << " has " << n_atom_bonds
			       << " bonds " << std::endl;
		  if (n_atom_bonds > 0) {
		     for (int ib=0; ib<n_atom_bonds; ib++) {
			mmdb::Atom *at_bonded = atom_bonds[ib].atom;
			int idx_1;
			int idx_2;
			at->GetUDData(udd_atom_index_handle, idx_1);
			at_bonded->GetUDData(udd_atom_index_handle, idx_2);
			if ((idx_1 < 0) || (idx_1 >= n_selected_atoms)) {
			   std::cout << "internal bonds error idx_1 " << idx_1 << std::endl;
			} else {
			   if ((idx_2 < 0) || (idx_2 >= n_selected_atoms)) {
			      std::cout << "internal bonds error idx_2 " << idx_2 << std::endl;
			   } else {
			      v[idx_1].push_back(idx_2);
			   }
			}
			// v[idx_2].push_back(idx_1); // MakeBonds does both ways
		     }
		  }
	       }
	    }
	 }
      }
   }

   mol->RemoveBonds();

   return v;
}

std::vector<std::vector<unsigned int> >
coot::find_1_4_connections(const std::vector<std::vector<unsigned int> > &bonds_vec) {

   std::vector<std::vector<unsigned int> > v(bonds_vec.size());
   for (std::size_t i=0; i<bonds_vec.size(); i++)
      v[i].reserve(4);

   // needs backtracking check

   for (std::size_t i=0; i<bonds_vec.size(); i++) {
      // std::cout << "i " << i << std::endl;
      const std::vector<unsigned int> &v1 = bonds_vec[i];
      for (std::size_t j=0; j<v1.size(); j++) {
	 const std::vector<unsigned int> &v2 = bonds_vec[v1[j]];
	 // std::cout << "j " << v1[j] << std::endl;
	 for (std::size_t k=0; k<v2.size(); k++) {
	    if (v2[k] != i) {
	       // std::cout << "k " << v2[k] << std::endl;
	       const std::vector<unsigned int> &v3 = bonds_vec[v2[k]];
	       for (std::size_t l=0; l<v3.size(); l++) {
		  if (v3[l] != i) {
		     if (v3[l] != v1[j]) {
			// std::cout << "1-4: " << i << " " << v3[l] << std::endl;
			v[i].push_back(v3[l]);
		     }
		  }
	       }
	    }
	 }
      }
   }

   return v;

}
