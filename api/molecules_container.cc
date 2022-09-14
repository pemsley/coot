
#include "molecules_container.hh"
#include "ideal/pepflip.hh"

int
molecules_container_t::flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf) {
   
   int result = coot::pepflip(atom_sel.mol, rs.chain_id, rs.res_no, rs.ins_code, alt_conf);
   return result;
}
