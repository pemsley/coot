
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
   static std::atomic<bool> on_going_updating_map_lock;

   class gru_points_t {
   public:
      int model_gru_points_delta; // for the latest change, I mean
      int   map_gru_points_delta;
      float rmsd_of_difference_map;
      gru_points_t(float rmsd) {
         model_gru_points_delta = 0;
         map_gru_points_delta = 0;
         rmsd_of_difference_map = rmsd;
      }
      gru_points_t(float rmsd_diff_map_current, const gru_points_t &gru_points_prev) {
         model_gru_points_delta = 0;
         rmsd_of_difference_map = rmsd_diff_map_current;
         map_gru_points_delta = gru_points_delta(gru_points_prev);
      }
      int gru_points_delta(const gru_points_t &prev) {
         return int(10000.0 * (prev.rmsd_of_difference_map - rmsd_of_difference_map));
      }
      static int total(const std::vector<gru_points_t> &gru_point_history) {
         int sum = 0;
         for (const auto &item : gru_point_history) {
            sum += item.map_gru_points_delta;
         }
         return sum;
      }
   };
   std::vector<gru_points_t> gru_point_history; // map and model (model currently not used)

   void init() {
      imol_refinement_map = -1;
      imol_difference_map = -1;
      geometry_init_standard(); // do this by default now
   }

public:

   molecules_container_t() : ramachandrans_container(ramachandrans_container_t()) {init();}
   int imol_refinement_map; // direct access
   int imol_difference_map; // direct access

   // the test for these failing is spec.empty()
   coot::atom_spec_t atom_cid_to_atom_spec(int imol, const std::string &cid) const;
   coot::residue_spec_t residue_cid_to_residue_spec(int imol, const std::string &cid) const;

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
         bool ignore_lys_and_arg_flag = true; // 20221018-PE remove this flag when rotamer probabiity
                                              // tables are read from a binary file (and is fast enough
                                              // to include lys and arg).
         rot_prob_tables.fill_tables(ignore_lys_and_arg_flag);
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
            // something fun here. - whatever it is though, don't put it in this header.
         }
      }
   }
   // -------------------------------- generic utils -----------------------------------

   void display_molecule_names_table() const;
   bool is_valid_model_molecule(int) const;
   bool is_valid_map_molecule(int) const;

   // -------------------------------- geometry/dictionaries --------------------------------

   void geometry_init_standard();
   int load_dictionary_file(const std::string &monomer_cif_file_name);

   // -------------------------------- coordinates utils -----------------------------------

   int read_pdb(const std::string &file_name);
   int write_coordinates(int imol, const std::string &file_name) const;

   // returns either the specified atom or null if not found
   mmdb::Atom *get_atom(int imol, const coot::atom_spec_t &atom_spec) const;
   // returns either the specified residue or null if not found
   mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &residue_spec) const;
   // returns either the specified atom or null if not found
   mmdb::Atom *get_atom_using_cid(int imol, const std::string &cid) const;
   // returns either the specified residue or null if not found
   mmdb::Residue *get_residue_using_cid(int imol, const std::string &cid) const;

   int undo(int imol); // 20221016-PE not working yet

   int redo(int imol); // 20221016-PE not working yet

   // -------------------------------- map utils -------------------------------------------

   // return the imol for the new molecule
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
   int flip_peptide_using_cid(int imol, const std::string &cid, const std::string &alt_conf);


   // -------------------------------- coordinates validation ------------------------------

   // get the rotamer dodecs for the model, not const because it regenerates the bonds.
   coot::simple_mesh_t get_rotamer_dodecs(int imol);
   coot::simple_mesh_t ramachandran_validation_markup_mesh(int imol) const;
   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > ramachandran_validation(int imol) const;


   // -------------------------------- Coordinates and map validation ----------------------

   coot::validation_information_t density_fit_analysis(int imol_model, int imol_map);


   // -------------------------------- Gru Points ------------------------------------------

   // calling this adds to the gru_points history. Make this pairs when we add model scoring.
   //
   int calculate_new_gru_points(int imol_diff_map);

   int gru_points_total() const; // the sum of all the gru ponts accumulated

   // reset the gru_points (calls reset_the_gru_points()), updates the maps (using internal/clipper SFC)
   // so, update your contour lines meshes after calling this function.
   int connect_updating_maps(int imol_model, int imol_map_2fofc, int imol_map_fofc);
   // call this before calling connect_updating_maps(). Perhaps this should be associated with the model?
   // (currently we use a map because that is what Coot used before).
   void associate_data_mtz_file_with_map(int imol, const std::string &data_mtz_file_name,
                                         const std::string &f_col, const std::string &sigf_col,
                                         const std::string &free_r_col);


   // -------------------------------- Updating Maps ---------------------------------------

   void sfcalc_genmap(int imol_model,
                      int imol_map_with_data_attached,
                      int imol_updating_difference_map);

   coot::util::sfcalc_genmap_stats_t
   sfcalc_genmaps_using_bulk_solvent(int imol_model,
                                     int imol_map_with_data_attached,
                                     int imol_updating_difference_map);

   // add these
   //
   // delete residue
   // delete residues
   // add_terminal_residue
   // update_map

};

#endif // MOLECULES_CONTAINER_HH
