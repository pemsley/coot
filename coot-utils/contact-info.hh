
#ifndef CONTACT_INFO_HH
#define CONTACT_INFO_HH

#include "atom-selection-container.hh"

namespace coot {

   class contact_info {

    class contacts_pair { 
    public:
      int id1;
      int id2;
      contacts_pair(int id1_in, int id2_in) { 
	id1 = id1_in; 
	id2 = id2_in;
      }
    };

    std::vector<std::pair<std::string, mmdb::realtype> > atom_radii;
    void setup_atom_radii();
    mmdb::realtype get_radius(const std::string &element) const;

    void contacts_from_monomer_restraints(const atom_selection_container_t asc, 
			    std::map<mmdb::Residue *, dictionary_residue_restraints_t> &res_restraints); // non-const for map [] usage

    void setup_from_monomer_restraints(const atom_selection_container_t &asc,
				       int imol,
				       protein_geometry *geom_p);

  public:
    std::vector<contacts_pair> contacts;
    contact_info(mmdb::Contact *con_in, int nc) {
      for (int i=0; i<nc; i++) { 
	contacts.push_back(contacts_pair(con_in[i].id1, con_in[i].id2));
      }
    }
    contact_info(mmdb::PPAtom atom_selection, mmdb::Contact *con_in, int nc) {
      setup_atom_radii();
      for (int i=0; i<nc; i++) { 
	mmdb::Atom *at_1 = atom_selection[con_in[i].id1];
	mmdb::Atom *at_2 = atom_selection[con_in[i].id2];
	std::string ele_1 = at_1->element;
	std::string ele_2 = at_2->element;
	mmdb::realtype dx = at_1->x - at_2->x;
	mmdb::realtype dy = at_1->y - at_2->y;
	mmdb::realtype dz = at_1->z - at_2->z;
	mmdb::realtype dist_2 = dx*dx + dy*dy + dz*dz;
	mmdb::realtype dist = sqrt(dist_2);
	mmdb::realtype r1 = get_radius(ele_1);
	mmdb::realtype r2 = get_radius(ele_2);
	if (dist < (r1 + r2 + 0.1))
	  contacts.push_back(contacts_pair(con_in[i].id1, con_in[i].id2));
      }
    }
    // This contact_info constructor does not take the alt conf(s) into
    // account.  That is becuase (in the current scenario) the alt conf
    // selection has already taken place before we get here.  If you want
    // to account for alt confs, then you'll have to write a new
    // constructor.
    //
    contact_info(const atom_selection_container_t &asc,
		 const std::string &monomer_type,
		 int imol,
		 protein_geometry *geom_p);

    // Here we look up the contacts for each monomer in the atom
    // selection. We also allow descriptions of bonds between
    // monomers.   
    // 
    // Can throw a std::runtime_error.
    // 
    contact_info(const atom_selection_container_t &asc,
		 int imol,
		 protein_geometry *geom_p, 
		 const bonded_pair_container_t &bonded_pairs);

    contact_info(const atom_selection_container_t &asc,
		 int imol,
		 protein_geometry *geom_p, 
		 const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &link_bond_atoms);

    // like above, but we have link atom quads (selhnd is a selection
    // handle - usually all atoms in mol, but not necessarily).
    // 
    template<class T>
    contact_info(mmdb::Manager *mol,
		 int imol,
		 int selhnd,
		 const std::vector<T> &link_torsions,
		 protein_geometry *geom_p);

    void add_MSE_Se_bonds(const atom_selection_container_t &asc);
    int n_contacts() const { return contacts.size(); } 

    // one way only
    std::vector<std::vector<int> > get_contact_indices() const;

    // with reverses, e.g. 0->1 and 1->0 too 
    std::vector<std::vector<int> > get_contact_indices_with_reverse_contacts() const;

    void print() const; // debug info
   };

  contact_info getcontacts(const atom_selection_container_t &asc); 
  contact_info getcontacts(const atom_selection_container_t &asc, 
			   const std::string &monomer_type, int imol,
			   protein_geometry *geom_p); 

   
}

#endif // CONTACT_INFO_HH
