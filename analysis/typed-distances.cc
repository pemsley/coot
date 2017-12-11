
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

   mmdb::PPAtom atom_selection = 0;
   int n_atoms;
   mol->GetSelIndex(SelHnd, atom_selection, n_atoms);
   if (n_atoms) {
      float max_dist = 5;
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
	 setup_bin_distances(); // uses mol
	 unsigned int half_ww = 4; // residues either side (and self)
	 find_residues_within_window(half_ww);
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

		     unsigned int dist_bin_id = get_dist_bin_id(dist);
		     bin_distances[residue_p][atom_pair_bin_id][dist_bin_id]++;
		  }
	       }
	    }
	 }
      }
   }
}

void
coot::typed_distances::output() const {

   std::map<mmdb::Residue *, std::vector<std::vector<unsigned int> > >::const_iterator it;

   bool do_per_single_residue_table = false;
   bool do_windowed_residue_table = true;

   if (do_per_single_residue_table) {
      for(it=bin_distances.begin(); it!=bin_distances.end(); it++) {

	 std::cout << " " << residue_spec_t(it->first);
	 for (std::size_t j=0; j<6; j++) {
	    std::cout << " type-idx " << j;
	    for (unsigned int idx_dist=0; idx_dist<n_distance_bins; idx_dist++)
	       std::cout << " " << it->second[j][idx_dist];
	 }
	 std::cout << "\n";
      }
   }

   if (do_windowed_residue_table) {
      std::cout << "--------------- windowed reidue table -----------" << std::endl;
      for(it=bin_distances.begin(); it!=bin_distances.end(); it++) {
	 std::cout << " " << residue_spec_t(it->first);
	 const std::vector<mmdb::Residue *> &riw = residues_within_window.at(it->first);
	 for (std::size_t idx_at_pair=0; idx_at_pair<6; idx_at_pair++) {
	    std::cout << " type-idx " << idx_at_pair;
	    for (unsigned int idx_dist=0; idx_dist<n_distance_bins; idx_dist++) {
	       unsigned int sum = 0;
	       for (std::size_t i=0; i<riw.size(); i++)
		  sum += bin_distances.at(riw[i])[idx_at_pair][idx_dist];
	       std::cout << " " << sum;
	    }
	 }
	 std::cout << "\n";
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

void
coot::typed_distances::setup_bin_distances() {

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) {
	    mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	    bin_distances[residue_p].resize(6);
	    for (std::size_t idx_atom_pair_type=0; idx_atom_pair_type<6; idx_atom_pair_type++) {
	       bin_distances[residue_p][idx_atom_pair_type].resize(n_distance_bins);
	       for (std::size_t idx_dist=0; idx_dist<n_distance_bins; idx_dist++) {
		  bin_distances[residue_p][idx_atom_pair_type][idx_dist] = 0;
	       }
	    }
	 }
      }
   }
}

int
coot::typed_distances::get_atom_pair_bin_id(const atom_type_t &t1, const atom_type_t &t2) const {

   int bin_id = -1;

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
