
#include <string>
#include "lidia-core/cod-atom-types.hh"
#include "utils/coot-utils.hh"
#include "read-sm-cif.hh"
#include "coot-coord-utils.hh"
#include "atom-selection-container.hh"

mmdb::Residue *get_search_fragment_residue(const std::string &file_name) {

   mmdb::Residue *residue_p = 0;
   atom_selection_container_t asc = get_atom_selection(file_name, true, true, false);
   if (asc.read_success) {
      residue_p = coot::util::get_first_residue(asc.mol);
   }
   return residue_p;
}


std::pair<bool, mmdb::math::Graph> make_mmdb_graph_for_search_fragment() {

   // 20220126-PE get the search fragment from a PDB file for now.
   // Use the first residue

   bool status = false;
   mmdb::math::Graph    graph;
   return std::pair<bool, mmdb::math::Graph>(status, graph);
}

void find_the_ribose(const std::string &cif_file_name,
                     const coot::dictionary_residue_restraints_t &restraints,
                     const RDKit::RWMol &rdkm,
                     const std::vector<cod::atom_type_t> &atom_types) {

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

   std::string search_fragment_pdb_file_name = "search-fragment.pdb";
   mmdb::Residue *residue_for_search_fragment = get_search_fragment_residue(search_fragment_pdb_file_name);

   if (residue_for_search_fragment) {

      // and now get the residue to be search from the cif_file_name
      coot::smcif smcif;
      mmdb::Manager *mol = smcif.read_sm_cif(cif_file_name);
      mmdb::Residue *residue_p = coot::util::get_first_residue(mol);

      if (residue_p) {
         coot::graph_match_info_t gmi = coot::graph_match(residue_p, residue_for_search_fragment, false, false);
      } else {
         std::cout << "No residue from this cif file " << cif_file_name << std::endl;
      }
   } else {
      std::cout << "No search residue from search fragment pdb file " << search_fragment_pdb_file_name << std::endl;
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
