
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/reduce.hh"
#include "flev.hh"

int main(int argc, char **argv) {

   int status = 0;

   coot::protein_geometry geom;
   geom.set_verbose(false);
   geom.init_standard();

   std::string pdb_file_name = "1x8b.pdb";
   std::string ligand_code = "824";
   std::string chain_id = "A";
   std::string ins_code = "";
   int res_no = 901;

   if (argc == 2) {
      std::string code = argv[1];
      if (code == "2vtq") {
         pdb_file_name = "2vtq.pdb";
         ligand_code = "LZA";
         res_no = 1299;
      }
      if (code == "2cmf") {
         pdb_file_name = "2cmf.pdb";
         ligand_code = "F11";
         res_no = 1536;
      }
   }


   // Also test 5a3h, 2wot

   bool use_gemmi = false;
   atom_selection_container_t atom_sel = get_atom_selection(pdb_file_name, use_gemmi, true, false);
   if (atom_sel.read_success) {

      int cif_read_number =  50;
      geom.try_dynamic_add("LZA", cif_read_number);
      cif_read_number++;
      geom.try_dynamic_add("824", cif_read_number);
      int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
      bool try_autoload_if_needed = false;
      bool rs = geom.have_dictionary_for_residue_type("LZA", imol_enc, try_autoload_if_needed);
      std::cout << "have_dictionary_for_residue_type LZA: " << rs << std::endl;
      rs = geom.have_dictionary_for_residue_type("824", imol_enc, try_autoload_if_needed);
      std::cout << "have_dictionary_for_residue_type 824: " << rs << std::endl;
      rs = geom.have_dictionary_for_residue_type("F11", imol_enc, try_autoload_if_needed);
      std::cout << "have_dictionary_for_residue_type F11: " << rs << std::endl;

      std::string output_file_name = "test.svg";
      float radius = 4.2; // was 4.4; // was 4.8;
      int imol = coot::protein_geometry::IMOL_ENC_ANY;
      coot::reduce r(atom_sel.mol, imol);
      r.add_geometry(&geom);
      bool nuclear_flag = false;
      r.add_hydrogen_atoms(nuclear_flag);
      mmdb::Manager *ligand_mol = geom.mol_from_dictionary(ligand_code, imol, false);

      if (true) {
         atom_sel.mol->WritePDBASCII("test-flev-with-H.pdb");
      }

      if (ligand_mol)
         pli::fle_view_with_rdkit_internal(atom_sel.mol, imol, &geom,
                                           chain_id, res_no, ins_code, radius, "svg", output_file_name);

   }
   return status;
}
