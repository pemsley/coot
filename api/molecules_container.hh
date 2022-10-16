
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
   float starting_gru_score_map;
   float starting_gru_score_model;
   std::vector<std::pair<float, float> > gru_point_history;
   void init() {
      imol_refinement_map = -1;
      imol_difference_map = -1;
      starting_gru_score_map   = 0;
      starting_gru_score_model = 0;
      geometry_init_standard(); // do this by default now
   }

public:

   molecules_container_t() : ramachandrans_container(ramachandrans_container_t()) {init();}
   int imol_refinement_map; // direct access
   int imol_difference_map; // direct access

   coot::atom_spec_t atom_cid_to_atom_spec(const std::string &cid) const;
   coot::residue_spec_t residue_cid_to_residue_spec(const std::string &cid) const;

   coot::simple_mesh_t test_origin_cube() const;
   coot::molecule_t & operator[] (unsigned int imol) {
      // maybe this should throw an exception on out-of-range?
      return molecules[imol];
   }
   mmdb::Manager *mol(unsigned int imol) const {
      if (is_valid_model_molecule(imol)) {
         return molecules[imol].atom_sel.mol;
      } else {
         return nullptr;
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

   // -------------------------------- backup and saving -----------------------------------

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
   // -------------------------------- generic utils -----------------------------------

   bool is_valid_model_molecule(int) const;
   bool is_valid_map_molecule(int) const;

   // -------------------------------- geometry/dictionaries --------------------------------

   void geometry_init_standard();
   int load_dictionary_file(const std::string &monomer_cif_file_name);

   // -------------------------------- coordinates utils -----------------------------------

   int read_pdb(const std::string &file_name);

   // returns either the specified atom or null if not found
   mmdb::Atom *get_atom(int imol, const coot::atom_spec_t &atom_spec) const;
   // returns either the specified residue or null if not found
   mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &residue_spec) const;

   int undo(int imol);

   int redo(int imol);

   // -------------------------------- map utils -------------------------------------------

   int read_mtz(const std::string &file_name, const std::string &f, const std::string &phi, const std::string &weight,
                bool use_weight, bool is_a_difference_map);
   int writeMap(int imol, const std::string &file_name) const;
   float get_map_rmsd_approx(int imol_map) const;

   // not const because the internal state of a coot_molecule is changed
   coot::simple_mesh_t get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level);

   std::vector<coot::molecule_t::difference_map_peaks_info_t> difference_map_peaks(int imol_map, int imol_protein, float n_rmsd) const;


   // -------------------------------- coordinates modelling -------------------------------

   int auto_fit_rotamer(int imol, const std::string &chain_id, int res_no, const std::string &ins_code, const std::string &alt_conf,
                        int imol_map);

   int delete_atom(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                   const std::string &atom_name, const std::string &alt_conf);
   int delete_atom_using_cid(int imol, const std::string &cid);

   int delete_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   int delete_residue_using_cid(int imol, const std::string &cid);

   int delete_residue_atoms_with_alt_conf(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                                          const std::string &alt_conf);
   int delete_residue_atoms_using_cid(int imol, const std::string &cid);

   // return a useful message if the addition did not work
   std::pair<int, std::string> add_terminal_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   // updates imol_model (of course)
   int add_waters(int imol_model, int imol_map);
   std::vector<std::string> non_standard_residue_types_in_model(int imol) const;
   int delete_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   int fill_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   int mutate_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code, const std::string &res_type);
   int flip_peptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf);
   int flip_peptide(int imol, const std::string &cid, const std::string &alt_conf);


   // -------------------------------- coordinates validation ------------------------------

   // get the rotamer dodecs for the model, not const because it regenerates the bonds.
   coot::simple_mesh_t get_rotamer_dodecs(int imol);
   coot::simple_mesh_t ramachandran_validation_markup_mesh(int imol) const;
   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > ramachandran_validation(int imol) const;


   // -------------------------------- coordinates and map validation ----------------------

   coot::validation_information_t density_fit_analysis(int imol_model, int imol_map);


   // -------------------------------- Gru Points ------------------------------------------

   // calling this adds to the gru_points history
   float calculate_current_gru_points(int imol_model, int imol_map) const;

   // so that you can check by how much the gru_points have improved in the recent move.
   float get_previous_gru_points() const;

   // do these need to be public?
   void reset_the_gru_points(int imol);
   float get_gru_points_starting_model_geometry_score(int model) const;
   float get_gru_points_starting_map_rmsd(int imol_map) const;


   // -------------------------------- Updating Maps ---------------------------------------

   // reset the gru_points (calls reset_the_gru_points()), updates the maps (using internal/clipper SFC)
   // so, update your contour lines meshes after calling this function.
   int connect_updating_maps(int imol_model, int imol_map_2fofc, int imol_map_fofc);

   // call this before calling connect_updating_maps()
   void associate_data_mtz_file_with_model(int imol, const std::string &data_mtz_file_name,
                                           const std::string &f_col, const std::string &sigf_col);

   // add these
   //
   // delete residue
   // delete residues
   // add_terminal_residue
   // update_map

};

#endif // MOLECULES_CONTAINER_HH
