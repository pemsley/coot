#ifndef MOORHEN_H_BONDS_HH
#define MOORHEN_H_BONDS_HH

#include <string>
#include <vector>

namespace moorhen {

   class h_bond_atom {
   public:
      int serial;
      float x, y, z;
      float charge;
      float occ;
      float b_iso;
      std::string element;
      std::string name;
      int model;
      std::string chain;
      int res_no;
      std::string residue_name;
      std::string altLoc;
      h_bond_atom() { serial = -1; x = 0; y = 0; z = 0; charge = 0; occ = 0;
         b_iso = 0; model = -1; res_no = -1; }
   };

   class h_bond {
   public:

      h_bond_atom hb_hydrogen; // McDonald and Thornton H-bond algorithm
      h_bond_atom donor;
      h_bond_atom acceptor;
      h_bond_atom donor_neigh;
      h_bond_atom acceptor_neigh;

      double angle_1;  // degrees
      double angle_2;
      double angle_3;
      double dist;  // H-bond length
      bool ligand_atom_is_donor; // for use when hb_hydrogen is NULL -
                                 // no hydrogens in H-bond analysis.
      bool hydrogen_is_ligand_atom;
      bool bond_has_hydrogen_flag;
      h_bond() { angle_1 = 0; angle_2 = 0; angle_3 = 0; dist = 0;
         ligand_atom_is_donor = false;
         hydrogen_is_ligand_atom = false;
         bond_has_hydrogen_flag = false; }
   };

}

#endif // MOORHEN_H_BONDS_HH
