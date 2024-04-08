/*
 * analysis/typed-distances.cc
 *
 * Copyright 2017 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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


#include "coot-utils/coot-coord-utils.hh"
#include "typed-distances.hh"

void
coot::typed_distances::init() {

   if (mol) {
      int SelHnd = mol->NewSelection(); // d

      // select only from model 1, not hydrogens, FIXME
      mol->SelectAtoms(SelHnd, 1, "*",
		       mmdb::ANY_RES, "*",
		       mmdb::ANY_RES, "*",
		       "*","*","!H","*", mmdb::SKEY_NEW);

      generate(SelHnd);

      // done

      mol->DeleteSelection(SelHnd);
   }
}

void
coot::typed_distances::generate(int SelHnd) {

   float max_dist = 5.7; // should be passed

   mmdb::PPAtom atom_selection = 0;
   int n_atoms;
   mol->GetSelIndex(SelHnd, atom_selection, n_atoms);
   if (n_atoms) {
      mmdb::Contact *pscontact = NULL;
      int n_contacts;
      float min_dist = 0.01;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
	 for (int j=0; j<4; j++)
	    my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mol->SeekContacts(atom_selection, n_atoms,
			atom_selection, n_atoms,
			0, max_dist,
			0, // in same residue
			pscontact, n_contacts,
			0, &my_matt, i_contact_group); // makes reverses also
      if (n_contacts > 0) {
	 // std::cout << "found " << n_contacts << " contacts " << std::endl;
	 unsigned int half_ww = 4; // residues either side (and self)
	 find_residues_within_window(half_ww); // puts residues in residues_within_window
	 if (pscontact) {
	    for (int i=0; i<n_contacts; i++) {

	       mmdb::Atom *central_at = atom_selection[pscontact[i].id1];
	       mmdb::Atom *neighb_at  = atom_selection[pscontact[i].id2];

	       if (! in_self_or_bonded_residue(central_at, neighb_at)) {
		  mmdb::Residue *residue_p = central_at->GetResidue();
		  atom_type_t tc = get_type(central_at);
		  atom_type_t tn = get_type(neighb_at);
		  float dist = distance(central_at, neighb_at);
		  int atom_pair_bin_id = get_atom_pair_bin_id(tc, tn);
		  if (atom_pair_bin_id == -1) {
		     std::cout << " bad atom-pair-bin-id "
			       << atom_spec_t(central_at) << " "
			       << atom_spec_t(neighb_at) << std::endl;
		  } else {
		     // happy path

		     std::map<mmdb::Residue *, std::map<int, std::vector<float> > >::iterator it =
			residue_distances_map.find(residue_p);

		     if (it != residue_distances_map.end()) {
			std::map<int, std::vector<float> >::iterator it_inner =
			   it->second.find(atom_pair_bin_id);
			if (it_inner != it->second.end()) {
			   it_inner->second.push_back(dist);
			} else {
			   residue_distances_map[residue_p][atom_pair_bin_id].push_back(dist);
			}
		     } else {
			residue_distances_map[residue_p][atom_pair_bin_id].push_back(dist);
		     }
		  }
	       }
	    }
	 }
      }
   }
}

void
coot::typed_distances::output() const {

   // using a map in a function that is const leads to unclear error messages from the
   // compiler
   //

   bool do_windowed_residue_table = true;
   std::map<mmdb::Residue *, std::map<int, std::vector<float> > >::const_iterator it;
   std::map<mmdb::Residue *, std::map<int, std::vector<float> > >::const_iterator it_r;

   if (do_windowed_residue_table) {
      std::cout << "--------------- windowed residue table -----------" << std::endl;

      std::cout << "residue_distances_map size " << residue_distances_map.size() << std::endl;

      for(it=residue_distances_map.begin(); it!=residue_distances_map.end(); it++) {
	 std::string res_name = it->first->GetResName();
	 if (res_name != "HOH") {
	    const std::vector<mmdb::Residue *> &riw = residues_within_window.at(it->first);
	    std::vector<float> sums(6, 0.0);
	    std::vector<unsigned int> n_distances(6, 0);
	    for (int j=0; j<6; j++) {
	       for (std::size_t ir=0; ir<riw.size(); ir++) {
		  mmdb::Residue *res_inner = riw[ir];
		  it_r = residue_distances_map.find(res_inner);
		  if (it_r == residue_distances_map.end()) {
		     std::cout << "Oopps - bad residue for map " << std::endl;
		  } else {
		     // happy path
		     const std::map<int, std::vector<float> > &m = it_r->second;
		     std::map<int, std::vector<float> >::const_iterator it_inner = m.find(j);
		     if (it_inner == m.end()) {
			if (false)
			   std::cout << "Oopps failed to find bin-index " << j << " "
				     << residue_spec_t(riw[ir]) << std::endl;
		     } else {
			const std::vector<float> &v = it_inner->second;
			n_distances[j] += v.size();
			for (std::size_t i=0; i<v.size(); i++) {
			   // std::cout << "print j " << j << std::endl;
			   sums[j] += v[i];
			}
		     }
		  }
	       }
	    }

	    int total_n_distances = 0;
	    std::vector<float> fractions(6, 0.0);
	    std::vector<float> means(6, 0.0);
	    for (int j=0; j<6; j++)
	       total_n_distances += n_distances[j];

	    for (int j=0; j<6; j++) {
	       fractions[j] = static_cast<float>(n_distances[j])/static_cast<float>(total_n_distances);
	    }

	    std::cout << "Sums: " << residue_spec_t(it->first);
	    for (int j=0; j<6; j++)
	       std::cout << " type " << j << ":  " << n_distances[j] << " " << fractions[j] << " ";
	    std::cout << "\n";
	 }
      }
   }
}

coot::typed_distances::atom_type_t
coot::typed_distances::get_type(mmdb::Atom *at) const {

   atom_type_t t1(NONE);
   std::string ele(at->element);
   if (ele == " C") t1 = atom_type_t(C);
   if (ele == " O") t1 = atom_type_t(O);
   if (ele == " S") t1 = atom_type_t(O);
   if (ele == " N") t1 = atom_type_t(N);

   return t1;
}


// coot::typed_distances::bin_type_t
int
coot::typed_distances::get_atom_pair_bin_id(const atom_type_t &t1, const atom_type_t &t2) const {

   int bin_id = -1; // use int, because I can't iterate through enum vars

   if (t1 == atom_type_t(C)) {
      if (t2 == atom_type_t(C)) bin_id = 0;
      if (t2 == atom_type_t(O)) bin_id = 1;
      if (t2 == atom_type_t(N)) bin_id = 2;
   }
   if (t1 == atom_type_t(O)) {
      if (t2 == atom_type_t(C)) bin_id = 1;
      if (t2 == atom_type_t(O)) bin_id = 3;
      if (t2 == atom_type_t(N)) bin_id = 4;
   }
   if (t1 == atom_type_t(N)) {
      if (t2 == atom_type_t(C)) bin_id = 2;
      if (t2 == atom_type_t(O)) bin_id = 4;
      if (t2 == atom_type_t(N)) bin_id = 5;
   }

   // std::cout << "bind-id  " << bin_id << " from " << t1 << " " << t2 << std::endl;
      
   return bin_id;
}

// fills residues_within_window
void
coot::typed_distances::find_residues_within_window(int half_wl) {

   // fills std::map<mmdb::Residue *, std::vector<mmdb::Residue *> residues_within_window;
   // captures self also

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) {
	    mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	    for(int window_res=ires-half_wl; window_res<=ires+half_wl; window_res++) {
	       if (window_res>=0) {
		  if (window_res<nres) {
		     mmdb::Residue *r = chain_p->GetResidue(window_res);
		     if (r) {
			residues_within_window[residue_p].push_back(r);
		     }
		  }
	       }
	    }
	 }
      }
   }
   

}
