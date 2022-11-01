
#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <vector>
#include "coords/Cartesian.h"
#include "coords/ramachandran-container.hh"
#include "coot_molecule.hh"
#include "coot-utils/coot-rama.hh"
#include "utils/coot-utils.hh"
#include "ideal/simple-restraint.hh"
#include "atom-pull.hh"
#include "validation-information.hh"
#include "simple-mesh.hh"

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
      explicit gru_points_t(float rmsd) {
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

   // --------------------- refinement --------------------------

   // 201803004:
   // refinement now uses references to Xmaps.
   // A dummy_map is created and a reference to that is created. Then
   // the reference is reset to a real xmap in a molecule (imol_for_map).
   // But, for a reason I don't understand, the refinement crashes when I do that.
   // When the initial dummy_xmap doesn't go out of scope, then the refinement is OK.
   // So this static dummy map is the map that doesn't go out of scope.
   // We only need one of it, so it goes here, rather than get created every
   // time we do a refinement. It may need to be public in future.
   //
   // 20221018-PE:
   // Now that we are in api, then I am now no longer sure that this should be static
   // or what static means in WebAssembly.
   static clipper::Xmap<float> *dummy_xmap;
   float map_weight;

   static ctpl::thread_pool static_thread_pool; // does this need to be static?
   bool show_timings;

   coot::restraints_container_t *last_restraints;
   bool continue_threaded_refinement_loop;
   bool refinement_is_quiet;
   int cif_dictionary_read_number;
   // return the state of having found restraints.
   std::string adjust_refinement_residue_name(const std::string &resname) const;
   bool make_last_restraints(const std::vector<std::pair<bool,mmdb::Residue *> > &local_resiudes,
			     const std::vector<mmdb::Link> &links,
			     const coot::protein_geometry &geom,
			     mmdb::Manager *mol_for_residue_selection,
			     const std::vector<coot::atom_spec_t> &fixed_atom_specs,
			     coot::restraint_usage_Flags flags,
			     bool use_map_flag,
			     const clipper::Xmap<float> *xmap_p);
   coot::refinement_results_t refine_residues_vec(int imol,
                                                  const std::vector<mmdb::Residue *> &residues,
                                                  const std::string &alt_conf,
                                                  mmdb::Manager *mol);

   int find_serial_number_for_insert(int seqnum_new,
                                     const std::string &ins_code_for_new,
                                     mmdb::Chain *chain_p) const;
   atom_selection_container_t make_moving_atoms_asc(mmdb::Manager *residues_mol,
                                                    const std::vector<mmdb::Residue *> &residues) const;
   // return 0 if any of the residues in selection don't have (at least) bond
   // restraints.  Try to auto-load the dictionary cifs and try again.
   // The vector is a list of residues for which no restraints could be found.
   std::pair<int, std::vector<std::string> >
     check_dictionary_for_residue_restraints(int imol, mmdb::PResidue *SelResidues, int nSelResidues);
   std::pair<int, std::vector<std::string> >
     check_dictionary_for_residue_restraints(int imol, const std::vector<mmdb::Residue *> &residues);
   std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> >
   create_mmdbmanager_from_res_vector(const std::vector<mmdb::Residue *> &residues,
                                      int imol,
                                      mmdb::Manager *mol_in,
                                      std::string alt_conf);
   // simple mmdb::Residue * interface to refinement.  20081216
   coot::refinement_results_t
   generate_molecule_and_refine(int imol,  // needed for UDD Atom handle transfer
                                const std::vector<mmdb::Residue *> &residues,
                                const std::string &alt_conf,
                                mmdb::Manager *mol,
                                bool use_map_flag=true);
   bool refinement_immediate_replacement_flag = true;
   int imol_moving_atoms;
   enum moving_atoms_asc_t {
      NEW_COORDS_UNSET = 0,       // moving_atoms_asc_type values
      NEW_COORDS_ADD = 1,                 // not used?
      NEW_COORDS_REPLACE = 2,
      NEW_COORDS_REPLACE_CHANGE_ALTCONF = 3,
      NEW_COORDS_INSERT = 4,
      NEW_COORDS_INSERT_CHANGE_ALTCONF = 5};
   short int moving_atoms_asc_type;
   static void thread_for_refinement_loop_threaded();

   static std::atomic<bool> restraints_lock;
   static void get_restraints_lock(const std::string &calling_function_name);
   static void release_restraints_lock(const std::string &calling_function_name);
   static std::string restraints_locking_function_name; //  static because it is set by above
   
   bool particles_have_been_shown_already_for_this_round_flag;

   static std::vector<atom_pull_info_t> atom_pulls;
   static void all_atom_pulls_off();
   static void atom_pull_off(const coot::atom_spec_t &spec);
   static void atom_pulls_off(const std::vector<coot::atom_spec_t> &specs);
   std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > make_rotamer_torsions(const std::vector<std::pair<bool, mmdb::Residue *> > &local_residues) const;
   // this is like mini-rsr:
   int refine_direct(int imol, std::vector<mmdb::Residue *> rv, const std::string &alt_loc,
                     mmdb::Manager *mol);


   // --------------------- init --------------------------

   void init() {
      imol_refinement_map = -1;
      imol_difference_map = -1;
      geometry_init_standard(); // do this by default now
      refinement_immediate_replacement_flag = true; // 20221018-PE for WebAssembly for the moment
      imol_moving_atoms = -1;
      refinement_is_quiet = true;
      show_timings = true;
      cif_dictionary_read_number = 40;
      // refinement
      continue_threaded_refinement_loop = false;
      particles_have_been_shown_already_for_this_round_flag = false;
      map_weight = 50.0;
      map_sampling_rate = 2.2;
   }

public:

   molecules_container_t() : ramachandrans_container(ramachandrans_container_t()) {init();}

   int imol_refinement_map; // direct access
   int imol_difference_map; // direct access
   void set_imol_refinement_map(int i) { imol_refinement_map = i; }
   void set_map_weight(float w) { map_weight = w; }
   float get_map_weight() const { return map_weight; }

   // the test for these failing is spec.empty()
   coot::atom_spec_t atom_cid_to_atom_spec(int imol, const std::string &cid) const;
   coot::residue_spec_t residue_cid_to_residue_spec(int imol, const std::string &cid) const;

   coot::simple_mesh_t test_origin_cube() const;
   //! set the show_timings flag
   void set_show_timings(bool s) { show_timings = s; }

#ifdef SWIG
#else
   coot::molecule_t & operator[] (unsigned int imol) {
      // maybe this should throw an exception on out-of-range?
      return molecules[imol];
   }
#endif
   mmdb::Manager *get_mol(unsigned int imol) const { // 20221018-PE function name change
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

   std::string get_molecule_name(int imol) const;
   void display_molecule_names_table() const;
   bool is_valid_model_molecule(int) const;
   bool is_valid_map_molecule(int) const;
   int close_molecule(int imol);

   // -------------------------------- geometry/dictionaries --------------------------------

   void geometry_init_standard();
   int load_dictionary_file(const std::string &monomer_cif_file_name);

   // -------------------------------- coordinates utils -----------------------------------

   int read_pdb(const std::string &file_name);
   int get_monomer(const std::string &monomer_name);
   int get_monomer_from_dictionary(const std::string &comp_id, bool idealised_flag);
   // 20221030-PE nice to have one day:
   // int get_monomer_molecule_by_network_and_dict_gen(const std::string &text);

   int write_coordinates(int imol, const std::string &file_name) const;

   // Mode is "COLOUR-BY-CHAIN-AND-DICTIONARY"
   coot::simple_mesh_t get_bonds_mesh(int imol, const std::string &mode);

   // returns either the specified atom or null if not found
   mmdb::Atom *get_atom(int imol, const coot::atom_spec_t &atom_spec) const;
   // returns either the specified residue or null if not found
   mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &residue_spec) const;
   // returns either the specified atom or null if not found
   mmdb::Atom *get_atom_using_cid(int imol, const std::string &cid) const;
   // returns either the specified residue or null if not found
   mmdb::Residue *get_residue_using_cid(int imol, const std::string &cid) const;

   int undo(int imol);

   int redo(int imol);

   // -------------------------------- map utils -------------------------------------------

   // return the imol for the new molecule
   float map_sampling_rate;
   void set_map_sampling_rate(float msr) { map_sampling_rate = msr; }
   // return the new molecule number or -1 on failure
   int read_mtz(const std::string &file_name, const std::string &f, const std::string &phi, const std::string &weight,
                bool use_weight, bool is_a_difference_map);
   // return the new molecule number or -1 on failure
   int read_ccp4_map(const std::string &file_name, bool is_a_difference_map);
   int writeMap(int imol, const std::string &file_name) const;
   float get_map_rmsd_approx(int imol_map) const;

   // not const because the internal state of a coot_molecule is changed
   coot::simple_mesh_t get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level);

   std::vector<coot::molecule_t::difference_map_peaks_info_t> difference_map_peaks(int imol_map, int imol_protein, float n_rmsd) const;


   // -------------------------------- coordinates modelling -------------------------------

   int auto_fit_rotamer(int imol, const std::string &chain_id, int res_no, const std::string &ins_code, const std::string &alt_conf,
                        int imol_map);

   //where scope in ["ATOM","WATER","RESIDUE","CHAIN","MOLECULE"]
   int delete_using_cid(int imol, const std::string &cid, const std::string &scope);

   int delete_atom(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                   const std::string &atom_name, const std::string &alt_conf);
   int delete_atom_using_cid(int imol, const std::string &cid);

   int delete_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   int delete_residue_using_cid(int imol, const std::string &cid);

   int delete_residue_atoms_with_alt_conf(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                                          const std::string &alt_conf);
   int delete_residue_atoms_using_cid(int imol, const std::string &cid);

   int delete_chain_using_cid(int imol, const std::string &cid);

   //! @return a useful message if the addition did not work
   std::pair<int, std::string> add_terminal_residue_directly(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   //! @return a useful message if the addition did not work
   // std::pair<int, std::string> add_terminal_residue_directly_using_cid(int imol, const std::string &cid);
   // get rid of the pair as a return, so that I can compile the binding
   int add_terminal_residue_directly_using_cid(int imol, const std::string &cid);

   // updates imol_model (of course)
   int add_waters(int imol_model, int imol_map);
   std::vector<std::string> non_standard_residue_types_in_model(int imol) const;
   int delete_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   int fill_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   int flip_peptide(int imol, const coot::atom_spec_t &atom_spec, const std::string &alt_conf);
   int flip_peptide_using_cid(int imol, const std::string &atom_cid, const std::string &alt_conf);

   int mutate(int imol, const std::string &cid, const std::string &new_residue_type);

   int side_chain_180(int imol, const std::string &atom_cid);

   int move_molecule_to_new_centre(int imol, float x, float y, float z);
   coot::Cartesian get_molecule_centre(int imol) const;

   //! return the new molecule number (or -1 on no atoms selected)
   int copy_fragment_using_cid(int imol, const std::string &cid);
   //! return the new molecule number (or -1 on no atoms selected)
   int copy_fragment_using_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end);

   // -------------------------------- coordinates refinement ------------------------------

   // mode {SINGLE, TRIPLE, QUINTUPLE, HEPTUPLE, SPHERE, BIG_SPHERE, CHAIN, ALL};
   //
   int refine_residues_using_atom_cid(int imol, const std::string &cid, const std::string &mode);
   int refine_residues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                       const std::string &alt_conf, const std::string &mode);
   int refine_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end);

   void set_refinement_is_verbose() { refinement_is_quiet = false; }

   // -------------------------------- coordinates validation ------------------------------

   // get the rotamer dodecs for the model, not const because it regenerates the bonds.
   coot::simple_mesh_t get_rotamer_dodecs(int imol);
   coot::simple_mesh_t ramachandran_validation_markup_mesh(int imol) const;
   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > ramachandran_validation(int imol) const;

   void coot_all_atom_contact_dots_instanced(mmdb::Manager *mol, int imol);

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

   //! Given a point on the front clipping plane (x1, y1, z1) and a point on the back clipping plane (x2, y2, z2)
   //! this function searches imol_refinement_map (if set) to find a the centre of a blob above the contour level.
   //! Blobs at the "front" are selected in preference to blobs at the back.
   //! If no blob is found, then the first of the pair is false.
   //! If it is found, then the second is (obviously) the centre of the blob.
   //! 20221022-PE in future, this function should/will be provided with a list of displayed maps and their
   //! contour levels - but for now, it uses (only) imol_refinement_map.
   //! blob_under_pointer_to_screen_centre().
   std::pair<bool, clipper::Coord_orth> go_to_blob(float x1, float y1, float z1, float x2, float y2, float z2,
                                                   float contour_level);

#ifdef SWIG
   PyObject *simple_mesh_to_pythonic_mesh(const coot::simple_mesh_t &mesh);
   PyObject *get_pythonic_bonds_mesh(int imol);
   PyObject *get_pythonic_model_mesh(int imol, unsigned int mesh_index);
   PyObject *get_pythonic_map_mesh(int imol, float x, float y, float z, float radius, float contour_level);
#endif

};

#endif // MOLECULES_CONTAINER_HH
