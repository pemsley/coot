
#include <string>
#include "lidia-core/cod-atom-types.hh"
#include "utils/coot-utils.hh"
#include "read-sm-cif.hh"
#include "coot-coord-utils.hh"

void find_the_ribose(const std::string &cif_file_name,
                     const coot::dictionary_residue_restraints_t &restraints,
                     const RDKit::RWMol &rdkm, const std::vector<cod::atom_type_t> &atom_types) {

   // a ribose ring connected to a purine (atom names from ZUWNUA.DePu.cif
   //
   // C6: C[5]-3_3_0_0:H-3:N[5]-3_3_3:O[5]-3_3 474 type-3: C[5](N[5]C[5,6]C[5])(C[5]C[5]HH)(O[5]C[5])(H) type-4: C[5](N[5]C[5,6]C[5])(C[5]C[5]HH)(O[5]C[5])(H){1|C<4>,1|N<2>,1|O<2>,2|C<3>,3|H<1>}
   // C9: C[5]-3_3_3_0:C[5]-3_3_3_0:H-3:H-3 474 type-3: C[5](C[5]N[5]O[5]H)(C[5]C[5]HO)(H)2 type-4: C[5](C[5]N[5]O[5]H)(C[5]C[5]HO)(H)2{1|C<4>,1|H<1>,2|C<3>}
   // C8: C[5]-3_3_3_0:C[5]-3_3_0_0:H-3:O-3_0 474 type-3: C[5](C[5]C[5]HH)(C[5]O[5]CH)(OH)(H) type-4: C[5](C[5]C[5]HH)(C[5]O[5]CH)(OH)(H){1|H<1>,1|N<3>}
   // C7: C-3_3_0_0:C[5]-3_3_3_0:H-3:O[5]-3_3 474 type-3: C[5](C[5]C[5]HO)(O[5]C[5])(CHHO)(H) type-4: C[5](C[5]C[5]HO)(O[5]C[5])(CHHO)(H){1|N<3>,3|H<1>}
   // O1: C[5]-3_3_3_0:C[5]-3_3_3_0 358 type-3: O[5](C[5]C[5]N[5]H)(C[5]C[5]CH) type-4: O[5](C[5]C[5]N[5]H)(C[5]C[5]CH){1|O<2>,2|C<3>,3|H<1>}

   // level-2 names:
   //
   // 0: C[5]-3_3_0_0:H-3:N[5]-3_3_3:O[5]-3_3
   // 1: C[5]-3_3_3_0:C[5]-3_3_3_0:H-3:H-3
   // 2: C[5]-3_3_3_0:C[5]-3_3_0_0:H-3:O-3_0
   // 3: C-3_3_0_0:C[5]-3_3_3_0:H-3:O[5]-3_3
   // 4: C[5]-3_3_3_0:C[5]-3_3_3_0
   //
   // bonds go 0-1, 1-2, 2-3, 3-4, 4-0 or backwards

   // types from ZUWNUA (before sorter fix)
   // C6: C27 C[5]-3_3_0_0:H-3:N[5]-3_3_3:O[5]-3_3
   // C9: C28 C[5]-3_3_3_0:C[5]-3_3_3_0:H-3:H-3
   // C8: C29 C[5]-3_3_0_0:C[5]-3_3_3_0:H-3:O-3_0  C8: C[5]-3_3_3_0:C[5]-3_3_0_0:H-3:O-3_0 
   // C7: C30 C-3_3_0_0:C[5]-3_3_3_0:H-3:O[5]-3_3
   // O1: O7 C[5]-3_3_3_0:C[5]-3_3_3_0

   // types from FIKHAI
   // C6: C1 C[5]-3_3_0_0:H-3:N[5]-3_3_3:O[5]-3_3
   // C9: C2 C[5]-3_3_3_0:C[5]-3_3_3_0:H-3:H-3
   // C8: C3 C[5]-3_3_3_0:C[5]-3_3_0_0:H-3:O-3_3
   // C7: C4 C-3_3_0_0:C[5]-3_3_3_0:H-3:O[5]-3_3
   // O1: O2 C[5]-3_3_3_0:C[5]-3_3_3_0

   std::vector<std::string> ribose_fragment_types;
   ribose_fragment_types.push_back("C[5]-3_3_0_0:H-3:N[5]-3_3_3:O[5]-3_3");
   ribose_fragment_types.push_back("C[5]-3_3_3_0:C[5]-3_3_3_0:H-3:H-3");
   ribose_fragment_types.push_back("C[5]-3_3_3_0:C[5]-3_3_0_0:H-3:O-3_0");
   ribose_fragment_types.push_back("C-3_3_0_0:C[5]-3_3_3_0:H-3:O[5]-3_3");
   ribose_fragment_types.push_back("C[5]-3_3_3_0:C[5]-3_3_3_0");

   std::vector<std::pair<unsigned int, unsigned int> > fragment_atom_index_for_bond;
   fragment_atom_index_for_bond.push_back(std::make_pair(0, 1));
   fragment_atom_index_for_bond.push_back(std::make_pair(1, 2));
   fragment_atom_index_for_bond.push_back(std::make_pair(2, 3));
   fragment_atom_index_for_bond.push_back(std::make_pair(3, 4));
   fragment_atom_index_for_bond.push_back(std::make_pair(4, 0));

   std::vector<std::pair<std::string, std::string> > found_bonds;

   // name -> type
   std::map<std::string, std::string> atom_type_2s_map;

   // diagnostics - count the number of Hydrogen atoms
   //
   unsigned int n_Hydrogen_atoms = 0;
   for (unsigned int i=0; i<rdkm.getNumAtoms(); i++) {
      auto at_p = rdkm[i];
      int n = at_p->getAtomicNum();
      if (n == 1)
         n_Hydrogen_atoms++;
   }

   std::cout << "Molecule has " << n_Hydrogen_atoms << " Hydrogen atoms "<< std::endl;

   for (unsigned int i=0; i<atom_types.size(); i++) {
      auto at_p = rdkm[i];
      std::string name;
      try {
         at_p->getProp("name", name);
      }
      catch (const KeyErrorException &kee) {
         std::cout << kee.what() << " missing atom name for atom with index " << i << std::endl;
      }
      const auto &type = atom_types[i];
      if (true)
         std::cout << name << " level-2.str: " << type.level_2.string()
                   << " hash: " << type.hash_value << " type-3: "
                   << type.level_3 << " type-4: " << type.level_4 << std::endl;
      atom_type_2s_map[name] = type.level_2.string();
   }

   for (const auto &bond : restraints.bond_restraint) {
      std::string atom_name_1 = coot::util::remove_whitespace(bond.atom_id_1());
      std::string atom_name_2 = coot::util::remove_whitespace(bond.atom_id_2());
      std::string type_1 = atom_type_2s_map[atom_name_1];
      std::string type_2 = atom_type_2s_map[atom_name_2];
      auto it_1 = std::find(ribose_fragment_types.begin(), ribose_fragment_types.end(), type_1);
      auto it_2 = std::find(ribose_fragment_types.begin(), ribose_fragment_types.end(), type_2);
      if (it_1 != ribose_fragment_types.end()) {
         if (it_2 != ribose_fragment_types.end()) {
            std::cout << "Yey! " << atom_name_1 << " " << atom_name_2 << std::endl;
            // now are they a bonded pair in the interesting fragment?
            for (const auto &atom_index_pair_for_bond : fragment_atom_index_for_bond) {
               const std::string t_1 = ribose_fragment_types[atom_index_pair_for_bond.first];
               const std::string t_2 = ribose_fragment_types[atom_index_pair_for_bond.second];
               bool order_swap = false;
               bool found = false;
               if (type_1 == t_1) {
                  if (type_2 == t_2) {
                     found = true;
                  }
               }
               if (type_1 == t_2) {
                  if (type_2 == t_1) {
                     order_swap = true;
                     found = true;
                  }
               }
               if (found) {
                  std::pair<std::string, std::string> bp(atom_name_1, atom_name_2);
                  if (order_swap) bp = std::pair<std::string, std::string>(atom_name_2, atom_name_1);
                  found_bonds.push_back(bp);
               }
            }
         }
      }
   }

   if (found_bonds.size() == fragment_atom_index_for_bond.size()) {
      std::cout << "we found the ribose" << std::endl;
   } else {
      std::cout << "didn't find the ribose in " << cif_file_name << " "
                << found_bonds.size() << " bonds were found "
                << "with " << n_Hydrogen_atoms << " Hydrogen atoms" << std::endl;
   }

}

int main(int argc, char **argv) {

   if (argc > 1) {

      std::string cif_file_name = argv[1];

      coot::smcif smcif;
      mmdb::Manager *mol = smcif.read_sm_cif(cif_file_name);

      mmdb::Residue *residue = coot::util::get_first_residue(mol);

      coot::dictionary_residue_restraints_t restraints(residue);

      if (true) {
         std::cout << "restraints info: " << std::endl;
         std::cout << "            " << restraints.atom_info.size() << " atoms" << std::endl;
         std::cout << "            " << restraints.bond_restraint.size() << " bonds" << std::endl;

         if (false) {
            std::vector<coot::dict_atom> ai = restraints.atom_info;
            for (const auto &atom : ai) {
               std::cout << "   " << atom << std::endl;
            }
            for (const auto &bond : restraints.bond_restraint) {
               std::cout << "   " << bond << std::endl;
            }
         }
      }

      try {
         RDKit::RWMol rdkm = coot::rdkit_mol(restraints);

         // now do something with rdkm
         cod::atom_types_t types;
         std::vector<cod::atom_type_t> atom_types = types.get_cod_atom_types(rdkm);

         if (atom_types.size() == rdkm.getNumAtoms()) {
            find_the_ribose(cif_file_name, restraints, rdkm, atom_types);
         }
      }
      catch (const std::exception &e) {
         std::cout << "Exception " << e.what() << " when parsing " << cif_file_name << std::endl;
      }
   }
}
