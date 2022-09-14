
#include "molecules_container.hh"
#include "ideal/pepflip.hh"

int
molecules_container_t::flipPeptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      result = molecules[imol].flipPeptide(rs, alt_conf);
   }
   return result;
}
