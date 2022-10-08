
#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <vector>
#include "coot_molecule.hh"
#include "validation-information.hh"
#include "simple-mesh.hh"
#include "coords/Cartesian.h"
#include "coords/ramachandran-container.hh"
#include "coot-utils/coot-rama.hh"

class molecules_container_t {

   std::vector<coot::molecule_t> molecules;
   coot::protein_geometry geom;
   ramachandrans_container_t ramachandrans_container;

public:
   molecules_container_t() : ramachandrans_container(ramachandrans_container_t()) {}

   bool is_valid_model_molecule(int) const;
   bool is_valid_map_molecule(int) const;

   int flipPeptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf);
   int flipPeptide(int imol, const std::string &cid, const std::string &alt_conf);
   int read_pdb(const std::string &file_name);
   int read_mtz(const std::string &file_name, const std::string &f, const std::string &phi, const std::string &weight,
                bool use_weight, bool is_a_difference_map);
   coot::validation_information_t density_fit_analysis(int imol_model, int imol_map);
   coot::simple_mesh_t test_origin_cube() const;
   coot::simple_mesh_t ramachandran_validation_markup_mesh(int imol) const;
   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > ramachandran_validation(int imol) const;
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
   bool contains_unsaved_models() const {
      for (const auto &m : molecules) {
         if (m.have_unsaved_changes()) return true;
      }
      return false;
   }
   void save_unsaved_model_changes() {
      for (const auto &m : molecules) {
         if (m.have_unsaved_changes()) {
            // something fun here.
         }
      }
   }
   // returns either the specified atom or null if not found
   mmdb::Atom *get_atom(int imol, const coot::atom_spec_t &atom_spec) const;
   // returns either the specified residue or null if not found
   mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &residue_spec) const;

   int writeMap(int imol, const std::string &file_name) const;

   // not const because the internal state of a coot_molecule is changed
   coot::simple_mesh_t get_map_contours_mesh(int imol, clipper::Coord_orth position, float radius, float contour_level);


};

#endif // MOLECULES_CONTAINER_HH
