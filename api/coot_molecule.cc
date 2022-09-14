
#include "coot_molecule.hh"
#include "ideal/pepflip.hh"

bool
coot_molecule_t::is_valid_model_molecule() {

   bool status = false;
   if (atom_sel.mol)
      status = true;
   return status;

}

bool
coot_molecule_t::is_valid_map_molecule() {

   bool status = false;
   if (! xmap.is_null()) {
      status = true;
   }
   return status;

}

int coot_molecule_t::flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = coot::pepflip(atom_sel.mol, rs.chain_id, rs.res_no, rs.ins_code, alt_conf);
   return result;

}
