
#ifndef COOT_H_BONDS_HH
#define COOT_H_BONDS_HH

#include <vector>
#include <map>
#include <algorithm>

#include <mmdb/mmdb_manager.h>

#include "protein-geometry.hh"
#include "coot-coord-utils.hh"

namespace coot {

   class h_bond {
   public:
      CAtom *hb_hydrogen; // McDonald and Thornton H-bond algorithm
      CAtom *donor;
      CAtom *acceptor;
      CAtom *donor_neigh;
      CAtom *acceptor_neigh;
      double angle_1;  // degrees
      double angle_2;
      double angle_3;
      double dist;  // H-bond length
      bool ligand_atom_is_donor; // for use when hb_hydrogen is NULL -
				 // no hydrogens in H-bond analysis.
      bool hydrogen_is_ligand_atom;
      bool bond_has_hydrogen_flag;
      
      h_bond() {
	 hb_hydrogen = NULL;
	 donor = NULL;
	 acceptor = NULL;
	 donor_neigh = NULL;
	 acceptor_neigh = NULL;
	 ligand_atom_is_donor = 1;
	 angle_1 = -1; 
	 angle_2 = -1; 
	 angle_3 = -1;
	 hydrogen_is_ligand_atom = 0; // no hydrogen
	 bond_has_hydrogen_flag = 0;
      }
      h_bond(CAtom *d, CAtom *a) {
	 hb_hydrogen = NULL;
	 donor = d;
	 acceptor = a;
	 donor_neigh = NULL;
	 acceptor_neigh = NULL;
	 ligand_atom_is_donor = 0;
	 dist = -1;
	 angle_1 = -1; 
	 angle_2 = -1; 
	 angle_3 = -1; 
	 hydrogen_is_ligand_atom = 0; // no hydrogen
	 bond_has_hydrogen_flag = 0;
      }

      // for McDonald and Thornton H-bonds, where the Hs are explicit.
      // One first of these atoms is the hydrogen, the other is the
      // acceptor.
      //
      // pass ligand_atom_is_H_flag as 1 when ligand atom is the H.
      // 
      h_bond(CAtom *h, CAtom *a, bool ligand_atom_is_H_flag) {
	 hb_hydrogen = h;
	 bond_has_hydrogen_flag = 1;
	 acceptor = a;
	 donor = NULL;
	 donor_neigh = NULL;
	 acceptor_neigh = NULL;
	 hydrogen_is_ligand_atom = ligand_atom_is_H_flag;
	 ligand_atom_is_donor = ligand_atom_is_H_flag;
	 dist = -1;
	 angle_1 = -1; 
	 angle_2 = -1; 
	 angle_3 = -1; 
      }
      bool operator<(const h_bond &hb_2) const {
	 return (residue_spec_t(donor) < residue_spec_t(hb_2.donor));
      }
      bool operator==(const h_bond &hb_2) const {
	 atom_spec_t sd1(donor);
	 atom_spec_t sa1(acceptor);
	 atom_spec_t sd2(hb_2.donor);
	 atom_spec_t sa2(hb_2.acceptor);
	 return ((sd1 == sd2) && (sa1 == sa2));
      }
      bool has_hydrogen() const { return bond_has_hydrogen_flag; }
      bool ligand_atom_is_H() const { return hydrogen_is_ligand_atom; }
      friend std::ostream & operator<<(std::ostream &s, h_bond hb);
   };
   std::ostream & operator<<(std::ostream &s, h_bond hb);


   class h_bonds {
      // return the udd handle (donor, acceptors and HB hydrogens, that is)
      int  mark_donors_and_acceptors(int donor_selHnd, int acceptor_selHnd, CMMDBManager *mol,
				     const protein_geometry &geom);
      
      // What is the nearest neighbour of the atoms in mol?
      // 
      std::map<CAtom *,  std::vector<std::pair<CAtom *, float> > >
      make_neighbour_map(int selHnd_1, int selHnd_2, CMMDBManager *mol);
      
      std::pair<bool, h_bond> 
      make_h_bond_from_ligand_hydrogen(CAtom *at_1, // H on ligand
				       CAtom *at_2, // acceptor on residue
				       const std::vector<std::pair<CAtom *, float> > &nb_1,
				       const std::vector<std::pair<CAtom *, float> > &nb_2) const;
      std::pair<bool, h_bond> 
      make_h_bond_from_environment_residue_hydrogen(CAtom *at_1, // acceptor on ligand
						    CAtom *at_2, // H on residue
						    const std::vector<std::pair<CAtom *, float> > &nb_1,
						    const std::vector<std::pair<CAtom *, float> > &nb_2) const;
      
   public:
      h_bonds() {}
      
      std::vector<h_bond>
      get(int selHnd_1, int selHnd_2, CMMDBManager *mol, const protein_geometry &geom);

      std::vector<h_bond>
      get_mcdonald_and_thornton(int selHnd_1, int selHnd_2, CMMDBManager *mol,
				const protein_geometry &geom,
				realtype max_dist =3.9);
      
      class atom_sorter {
	 CAtom *at;
	 coot::residue_spec_t at_res_spec;
      public:
	 atom_sorter(CAtom *at_in) {
	    at = at_in;
	    at_res_spec = coot::residue_spec_t(at);
	 }
	 bool operator()(const std::pair<CAtom *, float> &p1,
			 const std::pair<CAtom *, float> &p2) const {
	    
	    coot::residue_spec_t n1_res_spec(p1.first);
	    coot::residue_spec_t n2_res_spec(p2.first);
	    
	    // 
	    if ((n1_res_spec == at_res_spec) && (n2_res_spec == at_res_spec)) {
	       // they are both in the same residue as the donor
	       return (p1.second < p2.second);
	    } else {
	       if (n1_res_spec == at_res_spec) { 
		  // n_1 is in (and n_2 isn't), therefore n_1 is "shorter"
		  return 1;
	       } else {
		  if (n2_res_spec == at_res_spec) { 
		     // n_2 is in (and n_1 isn't), therefore n_2 is "shorter"
		     return 0;
		  } else {
		     // neither of them are in the residue
		     return (p1.second < p2.second);
		  } 
	       } 
	    } 
	 }
      };
   };
}

#endif // COOT_H_BONDS_HH

