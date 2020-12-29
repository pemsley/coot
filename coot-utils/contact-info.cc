
#include <algorithm>

#include "coot-coord-utils.hh"
#include "contact-info.hh"

// This contact_info constructor does not take the alt conf(s) into
// account.  That is becuase (in the current scenario) the alt conf
// selection has already taken place before we get here.  If you want
// to account for alt confs, then you'll have to write a new
// constructor.
//
coot::contact_info::contact_info(const atom_selection_container_t &asc, 
				 const std::string &monomer_type,
				 int imol,
				 coot::protein_geometry *geom_p) {

   // before messing about here, are you sure that you are looking at
   // the correct cif file for this residue?

   std::pair<bool, coot::dictionary_residue_restraints_t> r = 
      geom_p->get_monomer_restraints(monomer_type, imol);

   if (r.first) {
      std::map<std::string, map_index_t> name_map;
      for (int i=0; i<asc.n_selected_atoms; i++) {
	 std::string atom_name(asc.atom_selection[i]->name);
	 name_map[atom_name] = map_index_t(i);
      }

      for (unsigned int ib=0; ib<r.second.bond_restraint.size(); ib++) {
	 std::string n_1 = r.second.bond_restraint[ib].atom_id_1_4c();
	 std::string n_2 = r.second.bond_restraint[ib].atom_id_2_4c();
	 coot::map_index_t ind_1 = name_map[n_1];
	 coot::map_index_t ind_2 = name_map[n_2];
	 if (ind_1.is_assigned() && ind_2.is_assigned()) { 
	    contacts_pair p(ind_1.index(), ind_2.index());
	    contacts.push_back(p);
	 }
      }
   }
}

// we can throw an exeption if any restraints are not found.
//
// The atom selection here has already sifted out the unwanted alt confs.
// 
coot::contact_info::contact_info(const atom_selection_container_t &asc, int imol,
				 coot::protein_geometry *geom_p,
				 const coot::bonded_pair_container_t &bonded_pairs) {

   std::vector<mmdb::Residue *> residues;
   std::map<mmdb::Residue *, std::vector<int> > atoms_in_residue;

   // fill residues and atoms_in_residue
   for (int i=0; i<asc.n_selected_atoms; i++) {
      mmdb::Residue *r = asc.atom_selection[i]->residue;
      if (std::find(residues.begin(), residues.end(), r) == residues.end())
	 residues.push_back(r);
      atoms_in_residue[r].push_back(i);
   }
   std::map<mmdb::Residue *, coot::dictionary_residue_restraints_t> res_restraints;
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      std::string rn = residues[ires]->GetResName();
      std::pair<bool, coot::dictionary_residue_restraints_t> rest = 
	 geom_p->get_monomer_restraints(rn, imol);
      if (! rest.first) {
	 std::string m = "Restraints not found for type ";
	 m += rn;
	 throw std::runtime_error(m);
      }
      res_restraints[residues[ires]] = rest.second;
   }

   contacts_from_monomer_restraints(asc, res_restraints);

   // now handle the bonded_pairs (they have residue_1, residue_2 and a link name)
   for (unsigned int ib=0; ib<bonded_pairs.bonded_residues.size(); ib++) {
      mmdb::Residue *res_1 = bonded_pairs.bonded_residues[ib].res_1;
      mmdb::Residue *res_2 = bonded_pairs.bonded_residues[ib].res_2;
      if (std::find(residues.begin(), residues.end(), res_1) != residues.end()) { 
	 if (std::find(residues.begin(), residues.end(), res_2) != residues.end()) {
	    // OK, both residues were in the atom selection (as it should be)
	    std::string comp_id_1 = res_1->GetResName();
	    std::string comp_id_2 = res_2->GetResName();
	    std::string group_1 = geom_p->get_group(res_1);
	    std::string group_2 = geom_p->get_group(res_2);
	    std::vector<std::pair<coot::chem_link, bool> > mcl = 
	       geom_p->matching_chem_link(comp_id_1, group_1,
					  comp_id_2, group_2);
	    std::cout << "debug:: found " << mcl.size() << " matching chem links"
		      << std::endl;
	    // there should be just one mcl of course, but ... by the book...
	    for (unsigned int ilink=0; ilink<mcl.size(); ilink++) {
	       bool order_switch = mcl[ilink].second;
	       dictionary_residue_link_restraints_t lr = 
		  geom_p->link(mcl[ilink].first.Id()); // or is it chem_link_name?
	       if (lr.link_id != "") {
		  // non-empty link, i.e. it was looked up OK.
		  for (unsigned int ilr=0; ilr<lr.link_bond_restraint.size(); ilr++) { 
		     std::string link_bond_atom_name_1 = lr.link_bond_restraint[ilr].atom_id_1_4c();
		     std::string link_bond_atom_name_2 = lr.link_bond_restraint[ilr].atom_id_2_4c();

		     if (order_switch == false) { 
			for (unsigned int iat_1=0; iat_1<atoms_in_residue[res_1].size(); iat_1++) {
			   std::string atom_name_1 = asc.atom_selection[iat_1]->name;
			   if (link_bond_atom_name_1 == atom_name_1) {
			      for (unsigned int iat_2=0; iat_2<atoms_in_residue[res_2].size(); iat_2++) {
				 std::string atom_name_2 = asc.atom_selection[iat_2]->name;
				 if (link_bond_atom_name_2 == atom_name_2) {
				    contacts_pair p(iat_1, iat_2);
				    contacts.push_back(p);
				 }
			      }
			   }
			}
		     } else {

			// order switch 
			for (unsigned int iat_1=0; iat_1<atoms_in_residue[res_1].size(); iat_1++) {
			   std::string atom_name_1 = asc.atom_selection[iat_1]->name;
			   if (link_bond_atom_name_2 == atom_name_1) {
			      for (unsigned int iat_2=0; iat_2<atoms_in_residue[res_2].size(); iat_2++) {
				 std::string atom_name_2 = asc.atom_selection[iat_2]->name;
				 if (link_bond_atom_name_1 == atom_name_2) {
				    contacts_pair p(iat_1, iat_2);
				    contacts.push_back(p);
				 }
			      }
			   }
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
coot::contact_info::contacts_from_monomer_restraints(const atom_selection_container_t asc,
				       std::map<mmdb::Residue *, coot::dictionary_residue_restraints_t> &res_restraints) {
   
   // make contacts from monomer restraints
   // 
   for (int iat=0; iat<asc.n_selected_atoms; iat++) {
      mmdb::Atom *at_1 = asc.atom_selection[iat];
      std::string atom_name_1 = at_1->name;
      for (int jat=0; jat<asc.n_selected_atoms; jat++) {
	 // if they are in the same residue...
	 mmdb::Atom *at_2 = asc.atom_selection[jat];
	 if (at_1->residue == at_2->residue) { 
	    std::string atom_name_2 = at_2->name;
	    // was there a bond between them?
	    const std::vector<coot::dict_bond_restraint_t> &bond_restraints =
	       res_restraints[at_1->residue].bond_restraint;
	    for (unsigned int ibond=0; ibond<bond_restraints.size(); ibond++) {
	       if (bond_restraints[ibond].atom_id_1_4c() == atom_name_1) { 
		  if (bond_restraints[ibond].atom_id_2_4c() == atom_name_2) {
		     contacts_pair p(iat, jat);
		     contacts.push_back(p);
		     break;
		  }
	       }
	       // and the reverse indexing of that...
	       if (bond_restraints[ibond].atom_id_1_4c() == atom_name_2) { 
		  if (bond_restraints[ibond].atom_id_2_4c() == atom_name_1) {
		     contacts_pair p(jat, iat);
		     contacts.push_back(p);
		     break;
		  }
	       }
	    }
	 }
      } 
   }
}


void
coot::contact_info::setup_from_monomer_restraints(const atom_selection_container_t &asc,
						  int imol,
						  coot::protein_geometry *geom_p) {

   std::vector<mmdb::Residue *> residues;
   std::map<mmdb::Residue *, std::vector<int> > atoms_in_residue;

   // fill residues and atoms_in_residue
   for (int i=0; i<asc.n_selected_atoms; i++) {
      mmdb::Residue *r = asc.atom_selection[i]->residue;
      if (std::find(residues.begin(), residues.end(), r) == residues.end())
	 residues.push_back(r);
      atoms_in_residue[r].push_back(i);
   }
   std::map<mmdb::Residue *, coot::dictionary_residue_restraints_t> res_restraints;
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      std::string rn = residues[ires]->GetResName();
      std::pair<bool, coot::dictionary_residue_restraints_t> rest = 
	 geom_p->get_monomer_restraints(rn, imol);
      if (! rest.first) {
	 std::string m = "Restraints not found for type ";
	 m += rn;
	 throw std::runtime_error(m);
      }
      res_restraints[residues[ires]] = rest.second;
   }
   contacts_from_monomer_restraints(asc, res_restraints);
}

coot::contact_info::contact_info(const atom_selection_container_t &asc,
				 int imol,
				 coot::protein_geometry *geom_p, 
				 const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &link_bond_atoms) {


   setup_from_monomer_restraints(asc, imol, geom_p);
   
   // now the link bond restraints
   for (unsigned int ilb=0; ilb<link_bond_atoms.size(); ilb++) {
      bool ifound = 0;
      for (int i=0; i<asc.n_selected_atoms; i++) {
	 if (asc.atom_selection[i] == link_bond_atoms[ilb].first) {
	    for (int j=0; j<asc.n_selected_atoms; j++) {
	       if (asc.atom_selection[j] == link_bond_atoms[ilb].second) {
		  contacts_pair p(j, i);
		  // std::cout << "---- added link bond contact " << i << " " << j << std::endl;
		  contacts.push_back(p);
		  ifound = 1;
		  break;
	       }
	    }
	 }
	 if (ifound)
	    break;
      }
   }
}

template <class T>
coot::contact_info::contact_info(mmdb::Manager *mol, int imol,
				 int selhnd,
				 const std::vector<T> &link_torsions,
				 coot::protein_geometry *geom_p) {

   atom_selection_container_t asc(mol, selhnd);
   setup_from_monomer_restraints(asc, imol, geom_p);
   // now the bond between monomers (middle atoms must be in different residues).
   for (unsigned int itor=0; itor<link_torsions.size(); itor++) { 
      bool ifound = false;
      mmdb::Residue *r1 = link_torsions[itor].atom_2->residue;
      mmdb::Residue *r2 = link_torsions[itor].atom_3->residue;
      if (r1 != r2) {
	 for (int i=0; i<asc.n_selected_atoms; i++) {
	    if (asc.atom_selection[i] == link_torsions[itor].atom_2) {
	       for (int j=0; j<asc.n_selected_atoms; j++) {
		  if (asc.atom_selection[j] == link_torsions[itor].atom_3) {
		     contacts_pair p(j, i);
		     std::cout << "---- contact_info() constructor added link bond contact "
			       << i << " " << j << std::endl;
		     contacts.push_back(p);
		     ifound = true;
		     break;
		  }
	       }
	    }
	    if (ifound)
	       break;
	 }
      }
   }
}


// instantiate that:
template coot::contact_info::contact_info(mmdb::Manager *mol,
					  int imol,
					  int selhnd,
					  const std::vector<coot::torsion_atom_quad> &link_torsions,
					  coot::protein_geometry *geom_p);

// try to get the bonds/contacts from the dictionary.  If there are no
// bonds, then fall back to the distance based search.
coot::contact_info
coot::getcontacts(const atom_selection_container_t &asc,
		  const std::string &monomer_type,
		  int imol,
		  coot::protein_geometry *geom_p) {

   coot::contact_info ci(asc, monomer_type, imol, geom_p);
   if (ci.n_contacts() == 0)
      return coot::getcontacts(asc);
   return ci;
   
} 



coot::contact_info
coot::getcontacts(const atom_selection_container_t &asc) {

   // back here again, eh? :)
   // 
   // Yep 20100518 :)
   // Yep 20100702
   
   // std::cout << "DEBUG:: getcontacts() in mmdb-extras" << std::endl;

   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   float min_dist = 0.1;
   float max_dist = 2.4; // long!  Filtered later.
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   asc.mol->SeekContacts(asc.atom_selection, asc.n_selected_atoms,
			 asc.atom_selection, asc.n_selected_atoms,
			 min_dist, max_dist, // min, max distances
			 0,        // seqDist 0 -> in same res also
			 pscontact, n_contacts,
			 0, &my_matt, i_contact_group);

   coot::contact_info ci(asc.atom_selection, pscontact, n_contacts);

   // Now, do we need to handle MSE extra bonds?
   // 
   if (std::string(asc.atom_selection[0]->GetResName()) == "MSE") { 
      ci.add_MSE_Se_bonds(asc);
   }
   
   delete [] pscontact;
   return ci;
}

// one way only: 0: 1 2 3
// but 1: does not have 0 as a member index.
// 
std::vector<std::vector<int> >
coot::contact_info::get_contact_indices() const {
   
   std::vector<std::vector<int> > v;
   int max_index = 0; 
   for (unsigned int i=0; i<contacts.size(); i++) {
      if (contacts[i].id1 > max_index)
	 max_index = contacts[i].id1;
      if (contacts[i].id2 > max_index)
	 max_index = contacts[i].id2;
   }
   if (max_index > 0) { 
      v.resize(max_index+1);
      for (unsigned int i=0; i<contacts.size(); i++) {
	 v[contacts[i].id1].push_back(contacts[i].id2);
      } 
   }
   return v;
}

// with reverses, e.g. 0->1 and 1->0 too
//
std::vector<std::vector<int> >
coot::contact_info::get_contact_indices_with_reverse_contacts() const {

   std::vector<std::vector<int> > v = get_contact_indices();

   for (unsigned int i=0; i<v.size(); i++) { 
      for (unsigned int j=0; j<v[i].size(); j++) {

	 // so we have i -> j, now add j -> i (if i is not already in
	 // the list of j)

	 std::vector<int>::const_iterator it = std::find(v[v[i][j]].begin(), v[v[i][j]].end(), i);
	 if (it == v[v[i][j]].end()) // not found
	    v[v[i][j]].push_back(i);
      }
   }

   return v;
} 



void
coot::contact_info::add_MSE_Se_bonds(const atom_selection_container_t &asc) {

   int SE_index = -1;
   int CE_index = -1;
   int CG_index = -1;
   for (int i=0; i<asc.n_selected_atoms; i++) {
      std::string atom_name = asc.atom_selection[i]->name;
      if (atom_name == "SE  ") SE_index = i;
      if (atom_name == " CE ") CE_index = i;
      if (atom_name == " CG ") CG_index = i;
   }
   if (SE_index != -1) { 
      if (CE_index != -1) { 
	 if (CG_index != -1) {
	    contacts.push_back(coot::contact_info::contacts_pair(CG_index, SE_index));
	    contacts.push_back(coot::contact_info::contacts_pair(SE_index, CE_index));
	 }
      }
   }
}

void
coot::contact_info::setup_atom_radii() {

   atom_radii.resize(23);
   atom_radii[ 0] = std::pair<std::string, mmdb::realtype> (" C", 0.77);
   atom_radii[ 1] = std::pair<std::string, mmdb::realtype> (" N", 0.65);
   atom_radii[ 2] = std::pair<std::string, mmdb::realtype> (" O", 0.6);
   atom_radii[ 3] = std::pair<std::string, mmdb::realtype> (" H", 0.35);
   // atom_radii[ 4] = std::pair<std::string, mmdb::realtype> (" S", 0.9);
   atom_radii[ 4] = std::pair<std::string, mmdb::realtype> (" S", 1.1); // S-S bonds 2.16A?
   atom_radii[ 5] = std::pair<std::string, mmdb::realtype> (" P", 1.0);
   atom_radii[ 6] = std::pair<std::string, mmdb::realtype> ("SE", 1.15);
   atom_radii[ 7] = std::pair<std::string, mmdb::realtype> ("BR", 1.15);
   atom_radii[ 8] = std::pair<std::string, mmdb::realtype> ("CL", 1.0);
   atom_radii[ 9] = std::pair<std::string, mmdb::realtype> (" I", 1.4);
   atom_radii[10] = std::pair<std::string, mmdb::realtype> (" F", 0.5);
   atom_radii[11] = std::pair<std::string, mmdb::realtype> (" K", 2.2);
   atom_radii[12] = std::pair<std::string, mmdb::realtype> ("AS", 1.3);
   atom_radii[13] = std::pair<std::string, mmdb::realtype> ("NA", 1.8);
   atom_radii[14] = std::pair<std::string, mmdb::realtype> ("MG", 1.5);
   atom_radii[15] = std::pair<std::string, mmdb::realtype> ("AU", 1.4);
   atom_radii[16] = std::pair<std::string, mmdb::realtype> ("BE", 1.05);
   atom_radii[17] = std::pair<std::string, mmdb::realtype> ("FE", 1.4);
   atom_radii[18] = std::pair<std::string, mmdb::realtype> ("ZN", 1.35);
   atom_radii[19] = std::pair<std::string, mmdb::realtype> ("PD", 1.6);
   atom_radii[20] = std::pair<std::string, mmdb::realtype> ("PB", 1.46);
   atom_radii[21] = std::pair<std::string, mmdb::realtype> ("PT", 1.46);
   atom_radii[22] = std::pair<std::string, mmdb::realtype> ("AG", 1.36);
}

mmdb::realtype
coot::contact_info::get_radius(const std::string &element) const {

   mmdb::realtype r = 0.9;
   for (unsigned int i=0; i<atom_radii.size(); i++) {
      if (atom_radii[i].first == element) {
	 r = atom_radii[i].second;
	 break;
      } 
   } 
   return r;
}

void
coot::contact_info::print() const {

   std::vector<std::vector<int> > v = get_contact_indices();
   std::cout << " ===================================== " << std::endl;
   std::cout << " ======= size: " << v.size() << " ======== " << std::endl;
   std::cout << " ===================================== " << std::endl;

   for (unsigned int ic1=0; ic1<v.size(); ic1++) {
      std::cout << "  index " << ic1 << " : ";
      for (unsigned int ic2=0; ic2<v[ic1].size(); ic2++)
	 std::cout << v[ic1][ic2] << " ";
      std::cout << std::endl;
   }
   std::cout << "===" << std::endl;

}

