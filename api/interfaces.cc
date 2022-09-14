
#include "interfaces.hh"

#include "coot-utils/atom-selection-container.hh"
#include "ideal/pepflip.hh"

std::string flipPeptide(const std::string &pdb_file_name_in, const coot::residue_spec_t &rs,
                        const std::string &pdb_file_name_out) {

   std::string s;
   std::string alt_conf; // this should be an argument

   atom_selection_container_t asc = get_atom_selection(pdb_file_name_in);
   int result = coot::pepflip(asc.mol, rs.chain_id, rs.res_no, rs.ins_code, alt_conf);
   if (result != 0) {
      asc.mol->WritePDBASCII(pdb_file_name_out.c_str());
      s = pdb_file_name_out;
   }

   return s;
}

int flipPeptide_mmdb(mmdb::Manager *mol, const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = coot::pepflip(mol, rs.chain_id, rs.res_no, rs.ins_code, alt_conf);
   return result;

}
