
#include "molecules_container.hh"
#include "ideal/pepflip.hh"

bool
molecules_container_t::is_valid_model_molecule(int imol) {
   bool status = false;
   if (imol >= 0) {
      int ms = molecules.size();
      if (imol < ms) {
         status = molecules[imol].is_valid_model_molecule();
      }
   }
   return status;
}

bool
molecules_container_t::is_valid_map_molecule(int imol) {
   bool status = false;
   if (imol >= 0) {
      int ms = molecules.size();
      if (imol < ms) {
         status = molecules[imol].is_valid_map_molecule();
      }
   }
   return status;
}

int
molecules_container_t::flipPeptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      result = molecules[imol].flipPeptide(rs, alt_conf);
   }
   return result;
}


int
molecules_container_t::read_pdb(const std::string &file_name) {

   int status = -1;
   atom_selection_container_t asc = get_atom_selection(file_name);
   if (asc.read_success) {
      molecules.push_back(coot_molecule_t(asc));
      status = molecules.size() -1;
   }
   return status;
}
