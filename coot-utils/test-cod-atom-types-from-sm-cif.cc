
#include <string>
#include "lidia-core/cod-atom-types.hh"
#include "utils/coot-utils.hh"
#include "read-sm-cif.hh"
#include "coot-coord-utils.hh"

mmdb::math::Graph make_mmdb_graph_for_search_fragment() {

   mmdb::math::Graph graph;

   mmdb::math::Vertex **V;
   mmdb::math::Edge   **E;

   mmdb::Atom *at_1 = new::mmdb::Atom;
   mmdb::Atom *at_2 = new::mmdb::Atom;
   mmdb::Atom *at_3 = new::mmdb::Atom;
   mmdb::Atom *at_4 = new::mmdb::Atom;
   mmdb::Atom *at_5 = new::mmdb::Atom;

   return graph;
}

void find_the_ribose(const std::string &cif_file_name,
                     const coot::dictionary_residue_restraints_t &restraints,
                     const RDKit::RWMol &rdkm, const std::vector<cod::atom_type_t> &atom_types) {

   unsigned int n_Hydrogen_atoms = 0;
   for (unsigned int i=0; i<rdkm.getNumAtoms(); i++) {
      auto at_p = rdkm[i];
      int n = at_p->getAtomicNum();
      if (n == 1)
         n_Hydrogen_atoms++;
   }

   std::cout << "Molecule has " << n_Hydrogen_atoms << " Hydrogen atoms "<< std::endl;
   std::map<std::string, int> atom_name_hash_code_map;

   for (unsigned int i=0; i<rdkm.getNumAtoms(); i++) {
      auto at_p = rdkm[i];
      cod::atom_type_t full_type = atom_types[i];
      std::string name;
      try {
         at_p->getProp("name", name);
         atom_name_hash_code_map[name] = full_type.hash_value;
      }
      catch (const KeyErrorException &kee) {
         std::cout << kee.what() << " missing atom name for atom with index " << i << std::endl;
      }
   }

   // Ring carbon is 474, oxygen is 358
   std::vector<int> hashes_for_fragment = { 474, 474, 474, 474, 358};
   std::vector<std::pair<unsigned int, unsigned int> > fragment_atom_index_for_bond;
   fragment_atom_index_for_bond.push_back(std::make_pair(0, 1));
   fragment_atom_index_for_bond.push_back(std::make_pair(1, 2));
   fragment_atom_index_for_bond.push_back(std::make_pair(2, 3));
   fragment_atom_index_for_bond.push_back(std::make_pair(3, 4));
   fragment_atom_index_for_bond.push_back(std::make_pair(4, 0));
   
   mmdb::math::Graph graph = make_mmdb_graph_for_search_fragment();
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
