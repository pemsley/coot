
#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <vector>

#ifdef SWIG
#include "Python.h"
#endif

#include "coords/Cartesian.h"
#include "coords/ramachandran-container.hh"
#include "coot_molecule.hh"
#include "coot-utils/coot-rama.hh"
#include "coot-utils/coot-coord-extras.hh" // the missing atoms type
#include "utils/coot-utils.hh"
#include "ideal/simple-restraint.hh" // needed?
#include "atom-pull.hh"
#include "validation-information.hh"
#include "coot-utils/simple-mesh.hh"
#include "phi-psi-prob.hh"
#include "instancing.hh"
#include "coot-utils/coot-map-utils.hh"

//! the container of molecules. The class for all **libcootapi** functions.
class molecules_container_t {

   std::vector<coot::molecule_t> molecules;
   coot::protein_geometry geom;
   coot::rotamer_probability_tables rot_prob_tables;
   ramachandrans_container_t ramachandrans_container;
   static std::atomic<bool> on_going_updating_map_lock;
   bool draw_missing_residue_loops_flag;

   class rail_points_t {
   public:
      int model_rail_points_delta; // for the latest change, I mean
      int   map_rail_points_delta;
      float rmsd_of_difference_map;
      explicit rail_points_t(float rmsd) {
         model_rail_points_delta = 0;
         map_rail_points_delta = 0;
         rmsd_of_difference_map = rmsd;
      }
      rail_points_t(float rmsd_diff_map_current, const rail_points_t &rail_points_prev) {
         model_rail_points_delta = 0;
         rmsd_of_difference_map = rmsd_diff_map_current;
         map_rail_points_delta = rail_points_delta(rail_points_prev);
      }
      int rail_points_delta(const rail_points_t &prev) {
         float fudge = 2.4; // 20230117-PE makes 1000 rail points equal ~1% in R-factor for the tutorial data
         return int(100000.0 * fudge * (prev.rmsd_of_difference_map - rmsd_of_difference_map));
      }
      static int total(const std::vector<rail_points_t> &rail_point_history) {
         int sum = 0;
         for (const auto &item : rail_point_history) {
            sum += item.map_rail_points_delta;
         }
         return sum;
      }
   };
   std::vector<rail_points_t> rail_point_history; // map and model (model currently not used)

   class updating_maps_info_f {
   public:
      bool maps_need_an_update;
      int imol_model;
      int imol_2fofc;
      int imol_fofc;
      int imol_with_data_info_attached;
      updating_maps_info_f() {
         maps_need_an_update = false;
         imol_model = -1;
         imol_2fofc = -1;
         imol_fofc = -1;
         imol_with_data_info_attached = -1;
      }
   };
   updating_maps_info_f updating_maps_info;
   void set_updating_maps_need_an_update(int imol); // if model imol was changed, let's update the map when
                                                    // the next contouring mesh is requested.
                                                    // Checks the above information before acting, of course.
                                                    // No action if imol is the the model for updating maps.
   void update_updating_maps(int imol); // called from the get_map_contours_mesh() function

   coot::util::sfcalc_genmap_stats_t latest_sfcalc_stats;


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
   bool make_last_restraints(const std::vector<std::pair<bool,mmdb::Residue *> > &local_residues,
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
   int refine_direct(int imol, std::vector<mmdb::Residue *> rv, const std::string &alt_loc);

   double phi_psi_probability(const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) const;

   void read_standard_residues();
   atom_selection_container_t standard_residues_asc;

   // for auto-read mtz
   int valid_labels(const std::string &mtz_file_name, const std::string &f_col, const std::string &phi_col,
                    const std::string &weight_col, int use_weights) const;
   
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
      map_sampling_rate = 1.8;
      draw_missing_residue_loops_flag = true;
      read_standard_residues();
      make_backups_flag = true;
      interrupt_long_term_job = false;
   }

public:

   //! the one and only constructor
   explicit molecules_container_t(bool verbose=true) : ramachandrans_container(ramachandrans_container_t()) {
      if (! verbose) geom.set_verbose(false);
      init();
   }

   //! the refinement map - direct access. When refinement is performed, this is the map
   //! that will be used. Many (indeed most) of thesee functions explicity take a map. If the map
   //! is not known by the calling function then this map can be used as the map molecule index
   int imol_refinement_map; // direct access
   //! the difference map - direct access
   //!
   //! I am not sure that this is needed - or will ever be.
   int imol_difference_map; // direct access

   // -------------------------------- Basic Utilities -----------------------------------
   //! \name Basic Utilities

   //! Allow the user to disable/enable backups (`state` is `true` for "enable"). The default is `true`.
   void set_make_backups(bool state) { make_backups_flag = state; }
   //! @return the backup-enabled state
   bool get_make_backups() const { return make_backups_flag; }
   //! the backup-enable state (raw public if needed/prefered)
   static bool make_backups_flag; // does this need to be static?

   //! set the map used for refinement and fitting
   void set_imol_refinement_map(int i) { imol_refinement_map = i; }
   //! set the map weight
   void set_map_weight(float w) { map_weight = w; }
   //! @return the map weight
   float get_map_weight() const { return map_weight; }

   //! Convert atom cid string to a coot atom specifier.
   //! The test for these failing is `spec.empty()`
   coot::atom_spec_t atom_cid_to_atom_spec(int imol, const std::string &cid) const;

   //! Convert residue cid string to a coot residue specifier.
   //! @return the residues spec.  `spec.empty()` is true on failure.
   coot::residue_spec_t residue_cid_to_residue_spec(int imol, const std::string &cid) const;

   //! this set the show_timings flag. Various (not all) functions in this class can calculate how long
   //! they took to run. Setting this will write the time to taken (in milliseconds) to stdout.
   //! The default is `true`.
   void set_show_timings(bool s) { show_timings = s; }

   // -------------------------------- generic utils -----------------------------------
   //! \name Generic Utils

   //! @return the name of the molecule
   std::string get_molecule_name(int imol) const;
   //! debugging function: display the table of molecule and names
   void display_molecule_names_table() const;
   //! @return is this a valid model?
   bool is_valid_model_molecule(int imol) const;
   //! @return is this a valid map?
   bool is_valid_map_molecule(int imol_map) const;
   //! @return is this a difference map?
   bool is_a_difference_map(int imol_map) const;
   //! close the molecule (and delete dynamically allocated memory)
   //! @return 1 on successful closure and 0 on failure to close
   int close_molecule(int imol);

   //! @return the mesh of a unit solid cube at the origin
   coot::simple_mesh_t test_origin_cube() const;

#ifdef SWIG
#else
#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! don't use this in emscript
   coot::molecule_t & operator[] (unsigned int imol) {
      // maybe this should throw an exception on out-of-range?
      return molecules[imol];
   }
#endif
#endif

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! don't use this in emscript
   mmdb::Manager *get_mol(unsigned int imol) const { // 20221018-PE function name change
      if (is_valid_model_molecule(imol)) {
         return molecules[imol].atom_sel.mol;
      } else {
         return nullptr;
      }
   }
#endif

   //! fill the rotamer probability tables (currently not ARG and LYS)
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
   //! \name Backup and Saving

   //! @return a flag of unsaved models state - i.e. if any of them are unsaved, then this returns true.
   bool contains_unsaved_models() const {
      for (const auto &m : molecules) {
         if (m.have_unsaved_changes()) return true;
      }
      return false;
   }

   //! Save the unsaved model - this function has not yet been written!
   void save_unsaved_model_changes() {
      for (const auto &m : molecules) {
         if (m.have_unsaved_changes()) {
            // something fun here. - whatever it is though, don't put it in this header.
         }
      }
   }

   // -------------------------------- geometry/dictionaries --------------------------------
   //! \name Geometry and Dictionaries

   //! read the stardard list of residues
   void geometry_init_standard();

   //! @return a vector of non-standard residues (so that they can be used for auxiliary dictionary import)
   std::vector<std::string> non_standard_residue_types_in_model(int imol) const;

   // -------------------------------- coordinates utils -----------------------------------
   //! \name Coordinates Utils

   //! read a PDB file (or mmcif coordinates file, despite the name)
   //! @return the new molecule index on success and -1 on failure
   int read_pdb(const std::string &file_name);

   //! get the active atom given the screen centre
   //!
   //! ``displayed_model_molecules_list`` is a colon-separated list of molecules, *e.g.* "2:3:4"
   //! @return the molecule index and the atom cid. On failure (no molecules with atoms in them, say) then
   //! return -1 and a blank string.
   std::pair<int, std::string> get_active_atom(float x, float y, float z, const std::string &displayed_model_molecules_list) const;

   //! import a dictionary cif - imol_enc is used to specify to which molecule this dictionary should apply.
   //! Use IMOL_ENC_ANY to mean "it applies to all molecules."
   //!
   //! IMOL_ENC_ANY = -999999
   //! @return 1 on success and 0 on failure
   int import_cif_dictionary(const std::string &cif_file_name, int imol_enc);
   //! get a monomer
   //! @return the new molecule index on success and -1 on failure
   int get_monomer(const std::string &monomer_name);
   //! get a monomer for a particular molecule - use -999999 (IMOL_ENC_ANY) if no molecule-specific dictionary is needed.
   //! @return the new molecule index on success and -1 on failure
   int get_monomer_from_dictionary(const std::string &comp_id, int imol, bool idealised_flag);
   //! get monomer and place it at the given position for a particular molecule - use -999999 if no molecule-specific dictionary is needed
   //! @return the new molecule index on success and -1 on failure
   int get_monomer_and_position_at(const std::string &comp_id, int imol, float x, float y, float z);

   // 20221030-PE nice to have one day:
   // int get_monomer_molecule_by_network_and_dict_gen(const std::string &text);

   //! write the coordinate to the give file name
   //! @return 1 on success and 0 on failure
   int write_coordinates(int imol, const std::string &file_name) const;

   //! By default missing loops are drawn. This function allows missing loops to not be
   //! drawn. Sometimes that can clarify the representation. This is a lightweight function
   //! that sets a flag that is used by subsequent calls to ``get_bonds_mesh()``.
   void set_draw_missing_residue_loops(bool state);

   //! get the bonds mesh.
   //!
   //! ``mode`` is "COLOUR-BY-CHAIN-AND-DICTIONARY" - more modes to follow
   //!
   //! ``against_a_dark_background`` allows the bond colours to be relevant for the background.
   //! When the background is dark, the colours should (as a rule) be bright and pastelly.
   //! When the background is light/white, the colour darker and more saturated.
   //!
   //! ``smoothness_factor`` controls the number of triangles used to make the bond cylinders
   //! and spheres for the atoms - it rises in powers of 4. 1 is the smallest ``smoothness_factor``,
   //! 2 looks nice (but maybe is slower to transfer) and 3 is best.
   //!
   //! ``bond_width`` is the bond width in Angstroms. 0.12 is a reasonable default value.
   //!
   //! ``atom_radius_to_bond_width_ratio`` allows the representation of "ball and stick". To do so use a value
   //! between (say) 1.5 and 3.0. The ratio for "liquorice" representation is 1.0 (of course).
   //!
   //! @return a ``coot::simple_mesh_t``
   coot::simple_mesh_t get_bonds_mesh(int imol, const std::string &mode,
                                      bool against_a_dark_background,
                                      float bond_width, float atom_radius_to_bond_width_ratio,
                                      int smoothness_factor);

   //! get the instanced bonds mesh.
   //!
   //! The arguments are as above:
   //!
   //! ``mode`` is "COLOUR-BY-CHAIN-AND-DICTIONARY" - more modes to follow
   //!
   //! ``against_a_dark_background`` allows the bond colours to be relevant for the background.
   //! When the background is dark, the colours should (as a rule) be bright and pastelly.
   //! When the background is light/white, the colour darker and more saturated.
   //!
   //! ``smoothness_factor`` controls the number of triangles used to make the bond cylinders
   //! and spheres for the atoms - it rises in powers of 4. 1 is the smallest ``smoothness_factor``,
   //! 2 looks nice and 3 is best. Instancing may mean that smoothness factor 3 should
   //! be used by default.
   //!
   //! ``bond_width`` is the bond width in Angstroms. 0.12 is a reasonable default value.
   //!
   //! ``atom_radius_to_bond_width_ratio`` allows the representation of "ball and stick". To do so use a value
   //! between (say) 1.5 and 3.0. The ratio for "liquorice" representation is 1.0 (of course). 1.7 or 1.8
   //! looks nice.
   //!
   //! @return a ``coot::instanced_mesh_t``
   coot::instanced_mesh_t get_bonds_mesh_instanced(int imol, const std::string &mode,
                                                   bool against_a_dark_background,
                                                   float bond_width, float atom_radius_to_bond_width_ratio,
                                                   int smoothness_factor);

   //! If any colour rule has been set for this molecule, then we will use those. Otherwise, colorChainsScheme()
   //! will be called (and that his its internal colour-by-chain colouring scheme).
   //!
   void add_colour_rule(int imol, const std::string &selection_cid, const std::string &colour);

   //! delete the colour rules for the given molecule
   void delete_colour_rules(int imol);

   //! get the colour rules
   std::vector<std::pair<std::string, std::string> > get_colour_rules(int imol) const;

   //! print the colour rules
   void print_colour_rules(int imol) const;

   //! Update float parameter for MoleculesToTriangles molecular mesh
   void M2T_updateFloatParameter(int imol, const std::string &param_name, float value);

   //! Update int parameter for MoleculesToTriangles molecular mesh
   void M2T_updateIntParameter(int imol, const std::string &param_name, int value);

   //! get ribbon and surface representation
   coot::simple_mesh_t get_molecular_representation_mesh(int imol, const std::string &cid, const std::string &colour_scheme,
                                                         const std::string &style);

   //! get a Gaussian surface representation
   //!
   //! These values seem to give a reasonable quite smooth surface:
   //!
   //! sigma = 4.4
   //!
   //! contour_level = 4.0
   //!
   //! box_radius = 5.0
   //!
   //! grid_scale = 0.7
   //! @return a simple mesh composed of a number of Gaussian surfaces (one for each chain)
   coot::simple_mesh_t get_gaussian_surface(int imol, float sigma, float contour_level,
                                            float box_radius, float grid_scale) const;

   //! get chemical feaatures for the specified residue
   coot::simple_mesh_t get_chemical_features_mesh(int imol, const std::string &cid) const;

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! @returns either the specified atom or null if not found - don't use this in emscript
   mmdb::Atom *get_atom(int imol, const coot::atom_spec_t &atom_spec) const;
   //! @returns either the specified residue or null if not found - don't use this in emscript
   mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &residue_spec) const;
   //! @returns either the specified atom or null if not found - don't use this in emscript
   mmdb::Atom *get_atom_using_cid(int imol, const std::string &cid) const;
   //! @returns either the specified residue or null if not found - don't use this in emscript
   mmdb::Residue *get_residue_using_cid(int imol, const std::string &cid) const;
   //! get the atom position - don't use this in emscript
   std::pair<bool, coot::Cartesian> get_atom_position(int imol, coot::atom_spec_t &atom_spec);
#endif

   //! @return the number of atoms in the specified model, or 0 on error
   unsigned int get_number_of_atoms(int imol) const;

   //! @return vector of chain-ids for the given molecule
   std::vector<std::string> get_chains_in_model(int imol) const;
   //! @return vector of single letter codes - in a pair with the given residue spec
   std::vector<std::pair<coot::residue_spec_t, std::string> > get_single_letter_codes_for_chain(int imol, const std::string &chain_id) const;

   //! @return a list of residue that don't have a dictionary
   std::vector<std::string> get_residue_names_with_no_dictionary(int imol) const;

   //! @return an object that has information about residues without dictionaries and residues with missing atom
   //! in the the specified molecule
   std::vector<coot::residue_spec_t> residues_with_missing_atoms(int imol);

   //! Ths function is not const because missing_atoms() takes a non-const pointer to the geometry
   // (20230117-PE I should fix that)
   //!
   //! @return an object that has information about residues without dictionaries and residues with missing atom
   //! in the the specified molecule
   coot::util::missing_atom_info missing_atoms_info_raw(int imol);

   //! undo
   //! @return 1 on successful undo, return 0 on failure
   int undo(int imol);

   //! redo
   //! @return 1 on successful redo, return 0 on failure
   int redo(int imol);

   // -------------------------------- map utils -------------------------------------------
   //! \name Map Utils

   //! @return the map sampling rate (default is 1.8)
   float map_sampling_rate;

   //! set the map sampling rate (default is 1.8). Higher numbers mean smoother maps, but they take
   //! longer to generate, longer to transfer, longer to parse and longer to draw
   void set_map_sampling_rate(float msr) { map_sampling_rate = msr; }
   //! Read the given mtz file.
   //! @return the new molecule number or -1 on failure
   int read_mtz(const std::string &file_name, const std::string &f, const std::string &phi, const std::string &weight,
                bool use_weight, bool is_a_difference_map);
   //! Read the given mtz file.
   //! @return a vector of the maps created from reading the file
   std::vector<int> auto_read_mtz(const std::string &file_name);
   //! @return the new molecule number or -1 on failure
   int read_ccp4_map(const std::string &file_name, bool is_a_difference_map);
   //! write a map. This function was be renamed from ``writeMap``
   //! @return 1 on a successful write, return 0 on failure.
   int write_map(int imol, const std::string &file_name) const;
   //! @return the map rmsd (epsilon testing is not used). -1 is returned if `imol_map` is not a map molecule index.
   float get_map_rmsd_approx(int imol_map) const;

   //! get the mesh for the map contours.
   //!
   //! This function is not **const** because the internal state of a `coot_molecule_t` is changed.
   //! @return a `simple_mesh_t` for the map contours of the specified map
   coot::simple_mesh_t get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level);

   coot::util::sfcalc_genmap_stats_t get_latest_sfcalc_stats() const { return latest_sfcalc_stats; }

   class r_factor_stats {
      public:
      float r_factor;  // 0 to 1
      float free_r_factor;
      int rail_points_total;
      int rail_points_new;
   };

   r_factor_stats get_r_factor_stats();

   // -------------------------------- coordinates modelling -------------------------------
   //! \name Coordinates Modelling

   //! auto-fit rotamer
   //! @return 1 on successful modification, return 0 on failure
   int auto_fit_rotamer(int imol, const std::string &chain_id, int res_no, const std::string &ins_code, const std::string &alt_conf,
                        int imol_map);

   //! delete item
   //!
   //! where scope is one of the strings: ["ATOM","WATER","RESIDUE","CHAIN","MOLECULE", "LITERAL"]
   //! @return 1 on successful modification, return 0 on failure
   std::pair<int, unsigned int> delete_using_cid(int imol, const std::string &cid, const std::string &scope);

   //! delete atom
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_atom(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                   const std::string &atom_name, const std::string &alt_conf);
   //! delete atom using atom cid
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_atom_using_cid(int imol, const std::string &cid);

   //! delete residue
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   //! delete residue using cid
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_residue_using_cid(int imol, const std::string &cid);

   //! delete residue atoms using alt_conf
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_residue_atoms_with_alt_conf(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                                          const std::string &alt_conf);
   //! delete residue atoms using cid
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_residue_atoms_using_cid(int imol, const std::string &cid);

   //! delete side chain
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   //! delete side chain
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_side_chain_using_cid(int imol, const std::string &cid);

   //! delete chain.
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_chain_using_cid(int imol, const std::string &cid);

   //! delete the atoms specified in the CID selection
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_literal_using_cid(int imol, const std::string &cid);

   //! add a residue onto the end of the chain by fitting to density
   //! @return a first of 1 on success. Return a useful message in second if the addition did not work
   std::pair<int, std::string> add_terminal_residue_directly(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);
   //! @return a useful message if the addition did not work
   // std::pair<int, std::string> add_terminal_residue_directly_using_cid(int imol, const std::string &cid);
   //! This used to return a pair, but I removed it so that I could compile the binding
   int add_terminal_residue_directly_using_cid(int imol, const std::string &cid);

   //! add waters, updating imol_model (of course)
   //! @return 1 on a success, 0 on failure.
   int add_waters(int imol_model, int imol_map);

   //! add hydrogen atoms, updating imol_model (of course)
   //! @return 1 on success, 0 on failure.
   int add_hydrogen_atoms(int imol_model);

   //! delete hydrogen atoms, updating imol_model (of course)
   //! @return 1 on a successful deletion, 0 on failure.
   int delete_hydrogen_atoms(int imol_model);

   //! add an alternative conformation for the specified residue
   //! @return 1 on a successful addition, 0 on failure.
   int add_alternative_conformation(int imol_model, const std::string &cid);

   //! fill the specified residue
   //! @return 1 on a successful fill, 0 on failure.
   int fill_partial_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   //! fill the specified residue
   //! @return 1 on a successful fill, 0 on failure.
   int fill_partial_residue_using_cid(int imol, const std::string &cid);

   //! fill all the the partially-filled residues in the molecule
   //! @return 1 on a successful fill, 0 on failure.
   int fill_partial_residues(int imol);

   //! flip peptide
   //! @return 1 on a successful flip
   int flip_peptide(int imol, const coot::atom_spec_t &atom_spec, const std::string &alt_conf);
   //! flip peptide using an atom CID
   //! @return 1 on a successful flip
   int flip_peptide_using_cid(int imol, const std::string &atom_cid, const std::string &alt_conf);

   //! eigen-flip ligand
   void eigen_flip_ligand(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   //! eigen-flip ligand using CID
   void eigen_flip_ligand_using_cid(int imol, const std::string &residue_cid);

   //! mutate residue
   //! @return 1 on a successful move, 0 on failure.
   int mutate(int imol, const std::string &cid, const std::string &new_residue_type);

   //! rotate last chi angle of the side chain by 180 degrees
   //! @return 1 on a successful move, 0 on failure.
   int side_chain_180(int imol, const std::string &atom_cid);

   //! JED-Flip the ligand (or residue) at the specified atom.
   //! @return a non-blank message if there is a problem
   std::string jed_flip(int imol, const std::string &atom_cid, bool invert_selection);

   //! move the molecule to the given centre
   //! @return 1 on a successful move, 0 on failure.
   int move_molecule_to_new_centre(int imol, float x, float y, float z);

   //! get molecule centre
   //! @return the molecule centre
   coot::Cartesian get_molecule_centre(int imol) const;

   //! copy a fragment
   //! @return the new molecule number (or -1 on no atoms selected)
   int copy_fragment_using_cid(int imol, const std::string &cid);
   //! copy a residue-range fragment
   //! @return the new molecule number (or -1 on no atoms selected)
   int copy_fragment_using_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end);

   //! apply transformation to atom selection in the given molecule.
   //! @return the number of atoms moved.
   int apply_transformation_to_atom_selection(int imol, const std::string &atoms_selection_cid,
                                              int n_atoms, // for validation of the atom selection, (int because mmdb atom type)
                                              float m00, float m01, float m02,
                                              float m10, float m11, float m12,
                                              float m20, float m21, float m22,
                                              float c0, float c1, float c2, // the centre of the rotation
                                              float t0, float t1, float t2); // translation

   //! update the positions of the atoms in the residue
   int new_positions_for_residue_atoms(int imol, const std::string &residue_cid, std::vector<coot::molecule_t::moved_atom_t> &moved_atoms);

   //! update the positions of the atoms in the residues
   int new_positions_for_atoms_in_residues(int imol, const std::vector<coot::molecule_t::moved_residue_t> &moved_residues);

   //! ``list_of_other_molecules`` is a colon-separated list of molecules, *e.g.* "2:3:4"
   //! @return the first is a flag set to 1 if a merge occurred (and 0 if it did not)
   //! the second is a vector of merge results, i.e. if you merged a ligand, what is the new
   //! residue spec of the ligand, and if you merged a (polymer) chain, what is the new chain-id of
   //! that chain.
   std::pair<int, std::vector<merge_molecule_results_info_t> >
   merge_molecules(int imol, const std::string &list_of_other_molecules);

   //! this is called by the above function and is useful for other non-api functions (such as add_compound()).
   std::pair<int, std::vector<merge_molecule_results_info_t> >
   merge_molecules(int imol, std::vector<mmdb::Manager *> mols);

   //! Convert a cis peptide to a trans or vice versa.
   //! @return 1 on a successful conversion.
   int cis_trans_convert(int imol, const std::string &atom_cid);

   //! replace a fragment
   //! 
   //! _i.e._ replace the atoms of ``imol_base`` by those of the atom selection ``atom_selection`` in ``imol_reference``
   //! (``imol_base`` is the molecule that is modified).
   //! 
   //! @return the success status
   int replace_fragment(int imol_base, int imol_reference, const std::string &atom_selection);

   // -------------------------------- Coordinates Refinement ------------------------------
   //! \name Coordinates Refinement

   //! refine the residues
   //
   //! ``mode`` is one of {SINGLE, TRIPLE, QUINTUPLE, HEPTUPLE, SPHERE, BIG_SPHERE, CHAIN, ALL};
   //! @returns a value of 1 if the refinement was performed and 0 if it was not.
   int refine_residues_using_atom_cid(int imol, const std::string &cid, const std::string &mode);
   //! refine the residues
   //! @returns a value of 1 if the refinement was performed and 0 if it was not.
   int refine_residues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                       const std::string &alt_conf, const std::string &mode);
   //! refine residue range
   //! @returns a value of 1 if the refinement was performed and 0 if it was not.
   int refine_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end);

   //! for debugging the refinement - write out some diagnositics - some might be useful
   void set_refinement_is_verbose() { refinement_is_quiet = false; }

   //! Recently (20230202) the smoothness factor has been added as an extra argument
   //! `smoothness_factor` is 1, 2 or 3 (3 is the most smooth).
   //! @return the instanced mesh for the specified ligand
   coot::instanced_mesh_t contact_dots_for_ligand(int imol, const std::string &cid, unsigned int smoothness_factor) const;

   //! Recently (20230202) the smoothness factor has been added as an extra argument
   //! `smoothness_factor` is 1, 2 or 3 (3 is the most smooth).
   //! @return the instanced mesh for the specified molecule
   coot::instanced_mesh_t all_molecule_contact_dots(int imol, unsigned int smoothness_factor) const;

   // -------------------------------- Coordinates validation ------------------------------
   //! \name Coordinates Validation

   //! get the rotamer dodecs for the model, not const because it regenerates the bonds.
   //! @return a `coot::simple_mesh_t`
   coot::simple_mesh_t get_rotamer_dodecs(int imol);

   //! get the rotamer dodecs for the model, not const because it regenerates the bonds.
   //! @return an `instanced_mesh_t`
   coot::instanced_mesh_t get_rotamer_dodecs_instanced(int imol);

   //! get the ramachandran validation markup mesh
   //!
   //! 20221126-PE: the function was renamed from ``ramachandran_validation_markup_mesh()``.
   //! @return a `coot::simple_mesh_t`
   coot::simple_mesh_t get_ramachandran_validation_markup_mesh(int imol) const;
   //! get the data for Ramachandran validation, which importantly contains probability information
   //! @return a vector of `phi_psi_prob_t`
   std::vector<coot::phi_psi_prob_t> ramachandran_validation(int imol) const;

   // -------------------------------- Coordinates and map validation ----------------------
   //! \name Coordinates and Map Validation

   //! density fit validation information
   //! @returns a `coot::validation_information_t`
   coot::validation_information_t density_fit_analysis(int imol_model, int imol_map) const;

   //! density correlation validation information
   //! @returns a `coot::validation_information_t`
   coot::validation_information_t density_correlation_analysis(int imol_model, int imol_map) const;

   //! rotamer validation information
   //! @returns a `coot::validation_information_t`
   coot::validation_information_t rotamer_analysis(int imol_model) const;

   //! ramachandran validation information (formatted for a graph, not 3d)
   //! @returns a `coot::validation_information_t`
   coot::validation_information_t ramachandran_analysis(int imol_model) const;

   //! ramachandran validation information (formatted for a graph, not 3d) for a given chain in a given molecule
   //! 20230127-PE This function does not exist yet.
   //!
   //! @returns a `coot::validation_information_t`
   coot::validation_information_t ramachandran_analysis_for_chain(int imol_model, const std::string &chain_id) const;

   //! peptide omega validation information
   //! @returns a `coot::validation_information_t`
   coot::validation_information_t peptide_omega_analysis(int imol_model) const;

   //! get interesting places (does not work yet)
   //! @return a vector of `coot::validation_information_t`
   std::vector<coot::molecule_t::interesting_place_t> get_interesting_places(int imol, const std::string &mode) const;

   //! get difference map peaks
   //! @return a vector of `coot::validation_information_t`
   std::vector<coot::molecule_t::interesting_place_t> difference_map_peaks(int imol_map, int imol_protein, float n_rmsd) const;

   //! get pepflips based on the difference map
   //! @return a vector of `coot::validation_information_t`
   std::vector<coot::molecule_t::interesting_place_t> pepflips_using_difference_map(int imol_coords, int imol_difference_map, float n_sigma) const;

   //! calculate the MMRRCC for the residues in the chain
   //! Multi Masked Residue Range Corellation Coefficient
   std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
             std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
   mmrrcc(int imol, const std::string &chain_id, int imol_map) const;

   //! calculate the MMRRCC for the residues in the chain
   //! Multi Masked Residue Range Corellation Coefficient
   std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
             std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
   mmrrcc_internal(const atom_selection_container_t &asc,
                   const std::string &chain_id,
                   const clipper::Xmap<float> &xmap) const;
   // -------------------------------- Rail Points ------------------------------------------
   //! \name Rail Points!

   //! calling this adds to the rail_points history. Make this pairs when we add model scoring.
   //! @returns the new rail points (since last modification)
   int calculate_new_rail_points();

   //! the total rail points
   //! @returns the sum of all rail points accumulated since the maps were connected.
   int rail_points_total() const;

   // -------------------------------- Updating Maps ---------------------------------------
   //! \name Updating Maps

   //! associate a data mtz file with a molecule
   //!
   //! call this before calling connect_updating_maps().
   void associate_data_mtz_file_with_map(int imol, const std::string &data_mtz_file_name,
                                         const std::string &f_col, const std::string &sigf_col,
                                         const std::string &free_r_col);

   //! reset the rail_points (calls reset_the_rail_points()), updates the maps (using internal/clipper SFC)
   //! so, update your contour lines meshes after calling this function.
   //! @return 1 if the connection was successful.
   int connect_updating_maps(int imol_model, int imol_with_data_info_attached, int imol_map_2fofc, int imol_map_fofc);

   //! sfcalc and re-generate maps. This is a low-level function - generally one would use the updating maps
   //! method rather than this
   void sfcalc_genmap(int imol_model,
                      int imol_map_with_data_attached,
                      int imol_updating_difference_map);

   //! sfcalc and re-generate maps (low-level function). This functions uses bulk solvent.
   //!
   //! Call this function after connecting maps for updating maps to set the initial R-factor
   //! and store the initial map flatness.
   //!
   //! @return a class of interesting statistics
   coot::util::sfcalc_genmap_stats_t
   sfcalc_genmaps_using_bulk_solvent(int imol_model,
                                     int imol_2fofc_map,
                                     int imol_updating_difference_map,
                                     int imol_map_with_data_attached);

   // -------------------------------- Go To Blob ---------------------------------------
   //! \name Go to Blob

   //! Given a point on the front clipping plane (x1, y1, z1) and a point on the back clipping plane (x2, y2, z2)
   //! this function searches imol_refinement_map (if set) to find a the centre of a blob above the contour level.
   //! Blobs at the "front" are selected in preference to blobs at the back.
   //! If no blob is found, then the first of the pair is false.
   //! If it is found, then the second is (obviously) the centre of the blob.
   //!
   //! 20221022-PE: in future, this function should/will be provided with a list of displayed maps and their
   //! contour levels - but for now, it uses (only) imol_refinement_map.
   //! Use a string to pass the map information (map index and contour level), something like "0 0.45:1 1.2:2 0.1"
   std::pair<bool, clipper::Coord_orth> go_to_blob(float x1, float y1, float z1, float x2, float y2, float z2,
                                                   float contour_level);


   // -------------------------------- Ligand Functions ---------------------------------------
   //! \name Ligand Functions

   //! Ligand Fitting
   //!
   //! I am not yet clear what extra cut-offs and flags need to be added here.
   //! You can expect this to take about 20 seconds.
   //!
   //! For trivial (i.e non-flexible) ligands you should instead use the jiggle-fit algorithm, which
   //! takes a fraction of a second. (That is the algorithm used for "Add Other Solvent Molecules" in Coot.)
   //!
   //! @return a vector indices of molecules for the best fitting ligands to this blob.
   std::vector<int> fit_ligand_right_here(int imol_protein, int imol_map, int imol_ligand, float x, float y, float z,
                                          float n_rmsd, bool use_conformers, unsigned int n_conformers);

   //! "Jiggle-Fit Ligand"
   //! if n_trials is 0, then a sensible default value will be used.
   //! if translation_scale_factor is negative then a sensible default value will be used.
   //! @return a value less than -99.9 on failure to fit.
   float fit_to_map_by_random_jiggle(int imol, const coot::residue_spec_t &res_spec, int n_trials, float translation_scale_factor);

   //! "Jiggle-Fit Ligand" with a different interface
   //! as above, if n_trials is 0, then a sensible default value will be used.
   //! if translation_scale_factor is negative then a sensible default value will be used.
   //! @return a value less than -99.9 on failure to fit.
   float fit_to_map_by_random_jiggle_using_cid(int imol, const std::string &cid, int n_trials, float translation_scale_factor);

   //! This is a ligand function, not really a ligand-fitting function.
   //!
   //! It won't work unless the dictionary for that ligand has been imported.
   //! But more importantly than that, it doesn't work yet.
   //! Except for unusual cases, ``imol`` will be IMOL_ENC_ANY (-666666)
   //! @return the string for the SVG representation.
   std::string get_svg_for_residue_type(int imol, const std::string &comp_id) const;

   //! This function is for adding compounds/molecules like buffer agents and precipitants or anions and cations.
   //! _i.e._ those ligands that can be positioned without need for internal torsion angle manipulation.
   //!
   //! ``tlc`` is the three-letter-code/compound-id
   //!
   //! ``imol_dict`` is the molecule to which the ligand is attached (if any). Typically this will be IMOL_ENC_ANY (-666666).
   //!
   //! ``imol_map`` is the molecule number of the map that will be used for fitting.
   //!
   //! @return the success status, 1 or good, 0 for not good.
   int add_compound(int imol, const std::string &tlc, int imol_dict, int imol_map, float x, float y, float z);

   //! @return a vector of residue specifiers for the ligand residues - the residue name is encoded
   //! in the `string_user_data` data item of the residue specifier
   std::vector<coot::residue_spec_t> get_non_standard_residues_in_molecule(int imol) const;

   // -------------------------------- Testing -------------------------------------
   //! \name Testing functions

   class ltj_stats_t {
   public:
      unsigned int count;
      float function_value;
      std::chrono::time_point<std::chrono::high_resolution_clock> timer_start;
      std::chrono::time_point<std::chrono::high_resolution_clock> timer;
      ltj_stats_t() : timer(std::chrono::high_resolution_clock::now()) {
         count = 0;
         function_value = 0;
         timer_start = timer;
      }
      //! This function is called by the long-term job, udating the timer and count
      void update_timer() {
         timer = std::chrono::high_resolution_clock::now();
         count += 1;
      }
      //! This function is called by the interrogation function - and is to help
      //! the uer know how the job is going.
      //!
      //! 20230127-PE A nice graph of the change of the function value seems like a good idea
      double time_difference() {
         timer = std::chrono::high_resolution_clock::now();
         auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(timer - timer_start).count();
         return d10;
      }
   };

   //! long term job
   ltj_stats_t long_term_job_stats;

   // not for user-control
   bool interrupt_long_term_job;

   //! start a long-term job.
   //!
   //! if `n_seconds` is 0, then run forever (or until interrupted)
   //!
   void testing_start_long_term_job(unsigned int n_seconds);

   //! stop the long-term job runnning (testing function)
   void testing_stop_long_term_job();

   //! get the stats for the long-term job (testing function)
   ltj_stats_t testing_interrogate_long_term_job() { return long_term_job_stats; }


   // -------------------------------- Other ---------------------------------------

#ifdef SWIG
   //! \name Python functions

   PyObject *simple_mesh_to_pythonic_mesh(const coot::simple_mesh_t &mesh);
   PyObject *get_pythonic_bonds_mesh(int imol);
   PyObject *get_pythonic_model_mesh(int imol, unsigned int mesh_index);
   PyObject *get_pythonic_map_mesh(int imol, float x, float y, float z, float radius, float contour_level);
   PyObject *get_pythonic_molecular_representation_mesh(int imol, const std::string &atom_selection,
                                                        const std::string &colour_sheme,
                                                        const std::string &style);
#endif

};

#endif // MOLECULES_CONTAINER_HH
