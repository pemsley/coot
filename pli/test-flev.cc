
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/reduce.hh"
#include "flev.hh"

int main(int argc, char **argv) {

   int status = 0;

   coot::protein_geometry geom;
   geom.set_verbose(false);
   geom.init_standard();
   std::string pdb_file_name = "2vtq.pdb";
   std::string ligand_code = "LZA";
   std::string chain_id = "A";
   std::string ins_code = "";
   int res_no = 1299;

   pdb_file_name = "1x8b.pdb";
   ligand_code = "824";
   res_no = 901;

   bool use_gemmi = false;
   atom_selection_container_t atom_sel = get_atom_selection(pdb_file_name, use_gemmi, true, false);
   if (atom_sel.read_success) {

      std::string output_file_name = "test.svg";
      float radius = 4.4; // was 4.8;
      int imol = coot::protein_geometry::IMOL_ENC_ANY;
      coot::reduce r(atom_sel.mol, imol);
      bool nuclear_flag = false;
      r.add_hydrogen_atoms(nuclear_flag);
      mmdb::Manager *ligand_mol = geom.mol_from_dictionary(ligand_code, imol, false);
      if (ligand_mol)
         pli::fle_view_with_rdkit_internal(atom_sel.mol, imol, &geom,
                                           chain_id, res_no, ins_code, radius, "svg", output_file_name);

   }
   return status;
}
