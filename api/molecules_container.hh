
#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <vector>
#include "coot_molecule.hh"
#include "validation-information.hh"
#include "simple-mesh.hh"
#include "coords/Cartesian.h"
#include "coords/ramachandran-container.hh"
#include "coot-utils/coot-rama.hh"
#include "utils/coot-utils.hh"

class molecules_container_t {

   std::vector<coot::molecule_t> molecules;
   coot::protein_geometry geom;
   coot::rotamer_probability_tables rot_prob_tables;
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
   void fill_rotamer_probability_tables() {
      if (! rot_prob_tables.tried_and_failed()) {

         std::string tables_dir = coot::package_data_dir();
         char *data_dir = getenv("COOT_DATA_DIR");
         if (data_dir) {
            tables_dir = data_dir;
         }
         tables_dir += "/rama-data";
         rot_prob_tables.set_tables_dir(tables_dir);
         rot_prob_tables.fill_tables();
      }
   }


   // returns either the specified atom or null if not found
   mmdb::Atom *get_atom(int imol, const coot::atom_spec_t &atom_spec) const;
   // returns either the specified residue or null if not found
   mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &residue_spec) const;

   int writeMap(int imol, const std::string &file_name) const;

   // not const because the internal state of a coot_molecule is changed
   coot::simple_mesh_t get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level);

   // get the rotamer dodecs for the model, not const because it regenerates the bonds.
   coot::simple_mesh_t get_rotamer_dodecs(int imol);

   int auto_fit_rotamer(int imol, const std::string &chain_id, int res_no, const std::string &ins_code, const std::string &alt_conf,
                        int imol_map);

   int delete_atom(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                   const std::string &atom_name, const std::string &alt_conf);

   int delete_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   int delete_residue_atoms_with_alt_conf(int imol, const std::string &chain_id, int res_no, const std::string &ins_code, const std::string &alt_conf);

   // add these
   //
   // delete residue
   // delete residues
   // add_terminal_residue
   // update_map

};

#endif // MOLECULES_CONTAINER_HH
