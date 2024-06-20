
#include "coot-utils/atom-selection-container.hh"
#include "flev.hh"

int main(int argc, char **argv) {

   int status = 0;

   coot::protein_geometry geom;
   geom.set_verbose(false);
   geom.init_standard();
   std::string pdb_file_name = "test.pdb";
   bool use_gemmi = false;
   atom_selection_container_t atom_sel = get_atom_selection(pdb_file_name, use_gemmi, true, false);
   if (atom_sel.read_success) {

      std::string output_file_name = "test.svg";
      float radius = 4.8;
      std::string chain_id = "A";
      int res_no = 32;
      std::string ins_code = "";
      int imol = 0;
      pli::fle_view_with_rdkit_internal(atom_sel.mol, imol, &geom,
                                        chain_id, res_no, ins_code, radius, "svg", output_file_name);

   }
   return status;
}
