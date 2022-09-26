
#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <vector>
#include "coot_molecule.hh"
#include "validation-information.hh"
#include "simple-mesh.hh"
#include "coords/Cartesian.h"
#include "coot-utils/coot-rama.hh"

class molecules_container_t {

   std::vector<coot::molecule_t> molecules;

public:
   molecules_container_t() {}

   bool is_valid_model_molecule(int);
   bool is_valid_map_molecule(int);

   int flipPeptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf);
   int flipPeptide(int imol, const std::string &cid, const std::string &alt_conf);
   int read_pdb(const std::string &file_name);
   int read_mtz(const std::string &file_name, const std::string &f, const std::string &phi, const std::string &weight,
                bool use_weight, bool is_a_difference_map);
   coot::validation_information_t density_fit_analysis(int imol_model, int imol_map);
   coot::simple_mesh_t ramachandran_validation_markup_mesh(int imol);
   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > ramachandran_validation(int imol);
   coot::molecule_t & operator[] (unsigned int imol) {
      // maybe this should throw an exception on out-of-range?
      return molecules[imol];
   }
   mmdb::Manager *mol(unsigned int imol) {
      if (is_valid_model_molecule(imol)) {
         return molecules[imol].atom_sel.mol;
      } else {
         return nullptr;
      }
   }

};

#endif // MOLECULES_CONTAINER_HH
