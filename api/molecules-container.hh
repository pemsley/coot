
#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <memory>
#ifdef SWIG
#include "Python.h"
#endif

#include <vector>

#ifdef HAVE_SSMLIB
#include <ssm/ssm_align.h>
#endif

#if NB_VERSION_MAJOR // for flychecking
#include <nanobind/nanobind.h>
#endif

#include "compat/coot-sysdep.h"

#include "coords/Cartesian.hh"
#include "coords/ramachandran-container.hh"
#include "coot-molecule.hh"
#include "coot-utils/coot-rama.hh"
#include "coot-utils/coot-coord-extras.hh" // the missing atoms type
#include "coot-utils/coot-map-utils.hh"
#include "utils/coot-utils.hh"
#include "utils/setup-syminfo.hh"
#include "ideal/simple-restraint.hh" // needed?
#include "atom-pull.hh"
#include "validation-information.hh"
#include "superpose-results.hh"
#include "lsq-results.hh"
#include "coot-utils/simple-mesh.hh"
#include "coot-utils/texture-as-floats.hh"
#include "phi-psi-prob.hh"
#include "instancing.hh"
#include "coot-colour.hh" // put this in utils
#include "saved-strand-info.hh"
#include "svg-store-key.hh"
#include "moorhen-h-bonds.hh"
#include "header-info.hh"
#include "positioned-atom-spec.hh"
#include "user-defined-colour-table.hh"

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
#ifdef SKIP_FOR_PYTHON_DOXYGEN
#else
   //! Set updating maps need an update (private)
   //!
   //! If model imol was changed, let's update the map when the next contouring mesh is requested
   //!
   //! @param imol is the model molecule index
   updating_maps_info_f updating_maps_info;
   void set_updating_maps_need_an_update(int imol);
                                                    // Checks the above information before acting, of course.
                                                    // No action if imol is the the model for updating maps.

   //! Update the updating maps without generating a mesh (private)
   //!
   //! @param imol is the model molecule index
   void update_updating_maps(int imol); // called from the get_map_contours_mesh() function

   coot::util::sfcalc_genmap_stats_t latest_sfcalc_stats;
#endif
   // --------------------- superposition --------------------------

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
   float geman_mcclure_alpha;

   bool use_rama_plot_restraints;
   float rama_plot_restraints_weight;

   bool use_torsion_restraints;
   float torsion_restraints_weight;

   ctpl::thread_pool thread_pool;
   bool show_timings;

   coot::restraints_container_t *last_restraints;
   bool continue_threaded_refinement_loop;
   bool refinement_is_quiet;
   int cif_dictionary_read_number;

   //! @param resname is the 3 letter code for the residue, e.g. "ALA" for alanine
   //!
   //! @return the state of having found restraints
   std::string adjust_refinement_residue_name(const std::string &resname) const;
#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
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

   atom_selection_container_t make_moving_atoms_asc(mmdb::Manager *residues_mol,
                                                    const std::vector<mmdb::Residue *> &residues) const;

   int find_serial_number_for_insert(int seqnum_new,
                                     const std::string &ins_code_for_new,
                                     mmdb::Chain *chain_p) const;

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

   static void thread_for_refinement_loop_threaded();

   static std::atomic<bool> restraints_lock;
   static void get_restraints_lock(const std::string &calling_function_name);
   static void release_restraints_lock(const std::string &calling_function_name);
   static std::string restraints_locking_function_name; //  static because it is set by above

   bool particles_have_been_shown_already_for_this_round_flag;

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   int servalcat_refine_xray_internal(int imol, int imol_map, const std::string &output_prefix,
                                      const std::map<std::string, std::string> &key_value_pairs);
#endif


#ifdef SKIP_FOR_PYTHON_DOXYGEN
#else
   //! Get LSQ matrix internal (private)
   //!
   //! @param imol_ref the reference model molecule index
   //! @param imol_mov the moving model molecule index
   std::pair<short int, clipper::RTop_orth> get_lsq_matrix_internal(int imol_ref, int imol_mov, bool summary_to_screen) const;
#endif

   coot::validation_information_t
   get_q_score_validation_information(mmdb::Manager *mol, int udd_q_score, bool do_per_atom) const;

#endif


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

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   static void all_atom_pulls_off();
   static std::vector<atom_pull_info_t> atom_pulls;
   // nanobinds doesn't have a atom_spec_t, does it?
   static void atom_pull_off(const coot::atom_spec_t &spec);
   static void atom_pulls_off(const std::vector<coot::atom_spec_t> &specs);
#endif

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > make_rotamer_torsions(const std::vector<std::pair<bool, mmdb::Residue *> > &local_residues) const;
#endif

   //! Real space refinement.
   //!
   //! the `n_cycles` parameter allows partial refinement - so for an animated representation one would call this
   //! with a small number (10, 20, 100?) and call it again if the refine status is still yet to reach completion
   //! GSL_CONTINUE (-2). And then make a call to get the bonds mesh (or other molecular representation).
   //! If n_cycles is negative, this means "refine to completion."
   //!
   //! @return success/progress status

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! Refine direct
   int refine_direct(int imol, std::vector<mmdb::Residue *> rv, const std::string &alt_loc, int n_cycles);

   //! get phi,psi probability
   double phi_psi_probability(const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) const;
#endif

#ifdef SKIP_FOR_PYTHON_DOXYGEN
#else
   //! Read the standard protein, RNA, and DNA dictionaries (private)
   void read_standard_residues();

   std::map<svg_store_key_t, std::string> ligand_svg_store;

   atom_selection_container_t standard_residues_asc;
#endif

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else

   coot::graph_match_info_t overlap_ligands_internal(int imol_ligand, int imol_ref, const std::string &chain_id_ref,
                                                     int resno_ref, bool apply_rtop_flag);

   int install_model(const coot::molecule_t &m);

   superpose_results_t
   superpose_with_atom_selection(atom_selection_container_t asc_ref,
                                 atom_selection_container_t asc_mov,
                                 int imol_mov,
                                 std::string moving_mol_name,
                                 std::string reference_mol_name,
                                 bool move_copy_of_imol2_flag);
#endif

#ifdef HAVE_SSMLIB

   void print_ssm_sequence_alignment(ssm::Align *SSMAlign,
				     atom_selection_container_t asc_ref,
				     atom_selection_container_t asc_mov,
				     mmdb::PAtom *atom_selection1,
				     mmdb::PAtom *atom_selection2,
				     int n_selected_atoms_1, int n_selected_atoms_2,
				     bool move_copy_of_imol2_flag);

   coot::validation_information_t
   make_ssm_sequence_alignment_as_validation_information(ssm::Align *SSMAlign,
                                                    atom_selection_container_t asc_ref,
                                                    atom_selection_container_t asc_mov,
                                                    mmdb::PAtom *atom_selection1, mmdb::PAtom *atom_selection2,
                                                    int n_selected_atoms_1, int n_selected_atoms_2,
                                                    bool move_copy_of_imol2_flag);

   void make_and_print_horizontal_ssm_sequence_alignment(ssm::Align *SSMAlign,
							 atom_selection_container_t asc_ref,
							 atom_selection_container_t asc_mov,
							 mmdb::PAtom *atom_selection1,
							 mmdb::PAtom *atom_selection2,
							 int n_selected_atoms_1, int n_selected_atoms_2) const;

   void map_secondary_structure_headers(ssm::Align *SSMAlign,
					atom_selection_container_t asc_ref,
					atom_selection_container_t asc_mov,
					mmdb::PAtom *atom_selection1,
					mmdb::PAtom *atom_selection2,
					int n_selected_atoms_1, int n_selected_atoms_2) const;

   std::pair<std::string, std::string>
      get_horizontal_ssm_sequence_alignment(ssm::Align *SSMAlign,
					   atom_selection_container_t asc_ref,
					   atom_selection_container_t asc_mov,
					   mmdb::PAtom *atom_selection1, mmdb::PAtom *atom_selection2,
					   int n_selected_atoms_1, int n_selected_atoms_2) const;

   // for gesampt this will be vector of vector
   std::vector<std::pair<coot::residue_validation_information_t, coot::residue_validation_information_t> >
   get_pairs(ssm::Align *SSMAlign,
             atom_selection_container_t asc_ref,
             atom_selection_container_t asc_mov,
             mmdb::PAtom *atom_selection1, mmdb::PAtom *atom_selection2,
             int n_selected_atoms_1, int n_selected_atoms_2) const;


   void print_horizontal_ssm_sequence_alignment(std::pair<std::string, std::string> aligned_sequences) const;

   std::string generate_horizontal_ssm_sequence_alignment_string(const std::pair<std::string, std::string> &aligned_sequences) const;

#endif  // HAVE_SSMLIB

#ifdef SKIP_FOR_PYTHON_DOXYGEN
#else
   //! Check valid labels for auto-read mtz function (private)
   int valid_labels(const std::string &mtz_file_name, const std::string &f_col, const std::string &phi_col,
                    const std::string &weight_col, int use_weights) const;

   // water fitting settings
   float ligand_water_to_protein_distance_lim_max;
   float ligand_water_to_protein_distance_lim_min;
   float ligand_water_variance_limit;
   float ligand_water_sigma_cut_off;
#endif

   unsigned int max_number_of_simple_mesh_vertices;

   // --------------------- init --------------------------
#ifdef SKIP_FOR_PYTHON_DOXYGEN
#else
   void init(); // private
#endif

#ifdef SKIP_FOR_PYTHON_DOXYGEN
#else
   //! Write some debugging info to standard out
   void debug() const;

   bool map_is_contoured_using_thread_pool_flag;
   double contouring_time;
#endif

public:

   //! the one and only constructor
   explicit molecules_container_t(bool verbose=true);

   ~molecules_container_t();

   //! the refinement map - direct access. When refinement is performed, this is the map
   //! that will be used. Many (indeed most) of thesee functions explicity take a map. If the map
   //! is not known by the calling function then this map can be used as the map molecule index
   int imol_refinement_map; // direct access
   //! the difference map - direct access
   //!
   //! I am not sure that this is needed - or will ever be.
   int imol_difference_map; // direct access

   bool use_gemmi; // for mmcif and PDB parsing. 20240112-PE set to true by default in init()

   // -------------------------------- Basic Utilities -----------------------------------
   //! \name Basic Utilities

   //! Get the package version
   //!
   //! @return the package version, e.g. "1.1.11" - if this is a not yet a release version
   //! the version will end in a "+", such as "1.1.11+"
   std::string package_version() const;

   //! Set the state of using GEMMI for coordinates parsing
   //!
   //! @param state is True to mean that it is enabled. The default is True.
   void set_use_gemmi(bool state) { use_gemmi = state; }

   //! Get the state of using GEMMI for coordinates parsing
   bool get_use_gemmi() { return use_gemmi; }


   //! Allow the user to disable/enable backups
   //!
   //! @param state is True to mean that it is enabled. The default is True.
   void set_make_backups(bool state) { make_backups_flag = state; }

   //! Get the state of the backups
   //!
   //! @return the backup-enabled state
   bool get_make_backups() const { return make_backups_flag; }

   //! the backup-enable state (raw public if needed/preferred)
   bool make_backups_flag;

   //! File name to string
   //!
   //! @param file_name is the name of the given file
   //!
   //! @return the string of the contents of the given file-name.
   std::string file_name_to_string(const std::string &file_name) const;

   //! Get the number of molecules (either map or model)
   //!
   //! @return the number of molecules
   unsigned int get_number_of_molecules() const { return molecules.size(); }

   //! Add a number of empty molecules to the internal vector/list of molecules
   //!
   //! Note this is not like STL `reserve` as it will increase the molecule index
   //! of the next added molecule by `n_empty`.
   //!
   //! @param n_empty the number of empty molecules to create
   void create_empty_molecules(unsigned int n_empty);

   //! Set the map used for refinement and fitting
   //!
   //! @param i the map molecule index used for refinement and fitting
   void set_imol_refinement_map(int i) { imol_refinement_map = i; }

   //! Set the map weight
   //!
   //! @param w the map weight to be used for refinement, e.g. 50.0
   void set_map_weight(float w) { map_weight = w; }

   //! Get the map weight
   //!
   //! @return the map weight
   float get_map_weight() const { return map_weight; }

   //! Scale map
   //!
   //! @param imol is the model molecule index
   //! @param scale_factor is the scale factor
   void scale_map(int imol_map, float scale_factor);

   //! Convert atom cid string to a coot atom specifier
   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)
   //!
   //! @return the atom spec.,  `spec.empty()` is true on failure.
   coot::atom_spec_t atom_cid_to_atom_spec(int imol, const std::string &cid) const;

   //! Convert residue cid string to a coot residue specifier
   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //!
   //! @return the residues spec.,  `spec.empty()` is true on failure.
   coot::residue_spec_t residue_cid_to_residue_spec(int imol, const std::string &cid) const;

   //! Set the show_timings flag
   //!
   //! Various (not all) functions in this class can calculate how long
   //! they took to run. Setting this will write the time to taken (in milliseconds) to stdout.
   //!
   //! @param s is True to mean that it is enabled. The default is True.
   void set_show_timings(bool s) { show_timings = s; }

   // duplicate?
   // coot::protein_geometry & get_geom() { return geom; }

   //! Get header info
   //!
   //! (the header info is sparce at the moment)
   //!
   //! @param imol is the model molecule index
   //!
   //! @return an object with header info.
   moorhen::header_info_t get_header_info(int imol) const;

   //! Get imol_enc_any (enc: encoded)
   //!
   //! imol_enc_any refers to the molecule number for dictionary that can be used with any molecule
   //!
   //! @return the value of `imol_enc_any`
   int get_imol_enc_any() const;

   // -------------------------------- generic utils -----------------------------------
   //! \name Generic Utils

   //! Get the molecule name
   //!
   //! @param imol is the model molecule index
   //!
   //! @return the name of the molecule
   std::string get_molecule_name(int imol) const;

   //! Set the molecule name
   //!
   //! @param imol is the model or map molecule index
   //! @param new_name is the new name of the model or map
   void set_molecule_name(int imol, const std::string &new_name);

   //! Debugging function: display the table of molecule and names
   void display_molecule_names_table() const;

   //! Check if the model index is valid
   //!
   //! e.g. if the molecule is a map you will have an invalid model
   //!
   //! @param imol is the model molecule index
   //!
   //! @return True or False
   bool is_valid_model_molecule(int imol) const;

   //! Check if the map index is valid
   //!
   //! e.g. if the map is a model you will have an invalid map
   //!
   //! @param imol_map is the map molecule index
   //!
   //! @return True or False
   bool is_valid_map_molecule(int imol_map) const;

   //! Check if it the map is a difference map
   //!
   //! @param imol_map is the map molecule index
   //!
   //! @return True or False
   bool is_a_difference_map(int imol_map) const;

   //! Create an empty molecule
   //!
   //! @return the index of the new molecule
   int new_molecule(const std::string &name);

   //! Close the molecule (and delete dynamically allocated memory)
   //!
   //! @param imol is the model molecule index
   //!
   //! @return 1 on successful closure and 0 on failure to close
   int close_molecule(int imol);

   //! Delete the most recent/last closed molecule in the molecule vector, until the first
   //! non-closed molecule is found (working from the end)
   void end_delete_closed_molecules();

   //! Delete the most recent/last molecule in the molecule vector
   void pop_back();

   //! Delete all molecules
   void clear();

   //! Get the eigenvalues of the specified residue
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //!
   //! @return the eigenvalues of the atoms in the specified residue
   std::vector<double> get_eigenvalues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   //! Get a simple test mesh
   //!
   //! @return the mesh of a unit solid cube at the origin
   coot::simple_mesh_t test_origin_cube() const;

#ifdef SWIG
#else
#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! don't use this in emscript
   coot::molecule_t & operator[] (unsigned int imol);
#endif
#endif

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! don't use this in ecmascript
   mmdb::Manager *get_mol(unsigned int imol) const;
#endif

   //! Fill the rotamer probability tables (currently not ARG and LYS)
   void fill_rotamer_probability_tables();

   //! Access to a compressed file that contains the rotamer probabilities
   //!
   //! libcootapi will fill the rotamer probabilities tables from this compressed data stream
   //! (placeholder only)
   void accept_rotamer_probability_tables_compressed_data(const std::string &data_stream);

   // -------------------------------- backup and saving -----------------------------------
   //! \name Backup and Saving

   //! Check if there are unsaved changes for this model
   //!
   //! e.g. as yet not written to disk
   //!
   //! @return a flag of unsaved models state - e.g. if any of them are unsaved, then this returns True.
   bool contains_unsaved_models() const;

   //! Save the unsaved model - this function has not yet been written!
   void save_unsaved_model_changes();

   // -------------------------------- geometry/dictionaries --------------------------------
   //! \name Geometry and Dictionaries

   //! Read the standard list of residues
   void geometry_init_standard();

   //! Get a list of non-standard residues in the given molecule
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a vector/list of non-standard residues
   std::vector<std::string> non_standard_residue_types_in_model(int imol) const;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   //! Extract ligand restraints from the dictionary store and make an rdkit molecule
   //! Result to be eaten by C++ only.
   //!
   //! @param residue_name the residue name
   //! @param imol_enc the molecule for the ligand (typically is imol_enc_any)
   //! @return an RDKit RDMol.
   // RDKit::RWMol get_rdkit_mol(const std::string &residue_name, int imol_enc);

   // std::shared_ptr<RDKit::RWMol> get_rdkit_mol_shared(const std::string &residue_name, int imol_enc);

   //! get the 64base-encoded pickled string that represents the given residue/ligand name
   //!
   //! @param residue_name the residue name
   //! @param imol_enc the molecule for the ligand (typically is imol_enc_any)
   //! @return a pickle string, return an empty string on failure.
   std::string get_rdkit_mol_pickle_base64(const std::string &residue_name, int imol_enc);
#endif

   // -------------------------------- coordinates utils -----------------------------------
   //! \name Coordinates Utils

   //! Read a coordinates file (mmcif or PDB)
   //!
   //! @param file_name is the name of the coordinates file
   //!
   //! @return the new molecule index on success and -1 on failure
   int read_coordinates(const std::string &file_name);

   //! Read a PDB file (or mmcif coordinates file, despite the name)
   //!
   //! It does the same job as `read_coordinates` but has (perhaps) a more familiar name
   //!
   //! @param file_name is the name of the coordinates file
   //!
   //! @return the new molecule index on success and -1 on failure
   int read_pdb(const std::string &file_name);

   //! Read a small molecule CIF file
   //!
   //! @param file_name is the cif file-name
   //!
   //! @return the new molecule index on success and -1 on failure
   int read_small_molecule_cif(const std::string &file_name);

   //! Print the secondary structure information
   //!
   //! @param imol is the model molecule index
   void print_secondary_structure_info(int imol) const;

   //! Read a PDB file (or mmcif coordinates file, despite the name) to
   //! replace the current molecule
   //!
   //! This will only work if the molecules is already a model molecule
   //!
   //! @param imol is the model molecule index
   //! @param file_name is the name of the coordinates file
   void replace_molecule_by_model_from_file(int imol, const std::string &pdb_file_name);

   //! Split an NMR model into multiple models
   //!
   //! @param imol is the model molecule index
   //!
   //! @return the vector/list of new molecule indices
   std::vector<int> split_multi_model_molecule(int imol);

   //! Make a multi-model molecule given the input molecules
   //!
   //! @param model_molecules_list is a colon-separated list of molecules, (e.g. "2:3:4")
   //!
   //! @return the new molecule index, -1 if no models were found in the model molecules list
   int make_ensemble(const std::string &model_molecules_list);

   //! Get the molecule as a PDB string
   //!
   //! @param imol is the model molecule index
   //!
   //! @return the model molecule imol as a string. Return emtpy string on error
   std::string molecule_to_PDB_string(int imol) const;

   //! Get the molecule as an mmCIF string
   //!
   //! @param imol is the model molecule index
   //!
   //! @return the model molecule imol as a string. Return emtpy string on error
   std::string molecule_to_mmCIF_string(int imol) const;

   //! Get the active atom given the screen centre
   //!
   //! @param x is x position of the centre of the screen
   //! @param y is y position of the centre of the screen
   //! @param z is z position of the centre of the screen
   //! @param displayed_model_molecules_list is a colon-separated list of molecules, (e.g. "2:3:4")
   //!
   //! @return the molecule index and the atom cid. On failure (no molecules with atoms in them) then
   //! return -1 and a blank string.
   std::pair<int, std::string> get_active_atom(float x, float y, float z, const std::string &displayed_model_molecules_list) const;

   //! Import a dictionary cif
   //!
   //! @param imol_enc is used to specify to which molecule this dictionary should apply. Use IMOL_ENC_ANY to mean "it applies to all molecules."
   //! IMOL_ENC_ANY = -999999
   //!
   //! @return 1 on success and 0 on failure
   int import_cif_dictionary(const std::string &cif_file_name, int imol_enc);

   //! Get the cif file name
   //!
   //! @param comp_id is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine
   //! @param imol_enc is the molecule index for the residue type/compound_id
   //!
   //! @return the dictionary read for the given residue type, return an empty string on failure
   std::string get_cif_file_name(const std::string &comp_id, int imol_enc) const;

   //! Get the cif restraints as a string
   //!
   //! @param comp_id is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine
   //! @param imol_enc is the molecule index for the residue type/compound_id
   //!
   //! @return a string that is the contents of a dictionary cif file
   std::string get_cif_restraints_as_string(const std::string &comp_id, int imol_enc) const;

   //! Copy the dictionary that is specific for `imol_current` so that it can be used with a new molecule
   //!
   //! @param monomer_name is the 3 letter code of the monomer in the dictionary, e.g. "ALA" for alanine
   //! @param imol_current is the model molecule index with the dictionary to be copied from
   //! @param imol_new is the model molecule index the dictionary will be copied into
   bool copy_dictionary(const std::string &monomer_name, int imol_current, int imol_new);

   //! Get a monomer
   //!
   //! @param monomer_name is the 3-letter code of the monomer in the dictionary, e.g. "ALA" for alanine
   //!
   //! @return the new molecule index on success and -1 on failure
   int get_monomer(const std::string &monomer_name);

   //! Get a monomer for a particular molecule
   //!
   //! @param comp_id is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine
   //! @param imol is the model molecule index, use -999999 (IMOL_ENC_ANY) if no molecule-specific dictionary is needed
   //! @param idealised_flag means that the coordinates have been minimised with a molecular modelling minimisation algo,
   //!        usually the value is True
   //!
   //! @return the new molecule index on success and -1 on failure
   int get_monomer_from_dictionary(const std::string &comp_id, int imol, bool idealised_flag);

   //! Get monomer and place it at the given position for a particular molecule
   //!
   //! @param comp_id is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine
   //! @param imol is the model molecule index, use -999999 (IMOL_ENC_ANY) if no molecule-specific dictionary is needed
   //! @param x is the x value of the target position
   //! @param y is the y value of the target position
   //! @param z is the z value of the target position
   //!
   //! @return the new molecule index on success and -1 on failure
   int get_monomer_and_position_at(const std::string &comp_id, int imol, float x, float y, float z);

   //! Match atom between 2 dictionaries
   //!
   //! @param comp_id_1 is the 3-letter code for the residue/ligand in the first model, e.g. "ALA" for alanine
   //! @param imol_1 is the model molecule index of the first model
   //! @param comp_id_2 is the 3-letter code for the residue/ligand in the second model, e.g. "ALA" for alanine
   //! @param imol_2 is the model molecule index of the second model
   //!
   //! @return the atom name match on superposing the atoms of the given dictionaries
   std::map<std::string, std::string>
   dictionary_atom_name_map(const std::string &comp_id_1, int imol_1, const std::string &comp_id_2, int imol_2);

   //! get types
   std::vector<std::string> get_types_in_molecule(int imol) const;

   // 20221030-PE nice to have one day:
   // int get_monomer_molecule_by_network_and_dict_gen(const std::string &text);

   //! Get the groups for a vector/list of monomers
   //!
   //! e.g. "NON POLYMER", "PEPTIDE", etc
   //!
   //! @param residue_names is a list of residue names, e.g. ["ALA", "TRP"]
   //!
   //! @return the group for the given list of residue names
   std::vector<std::string> get_groups_for_monomers(const std::vector<std::string> &residue_names) const;

   //! Get the group for a particular monomer
   //!
   //! e.g. "NON POLYMER", "PEPTIDE", etc
   //! @param residue_name is the is the 3-letter code for the residue, e.g. "ALA" for alanine
   //!
   //! @return the group for the given residue name
   std::string get_group_for_monomer(const std::string &residue_name) const;

   //! Get the hydrogen bond type of a particular atom in a given residue type
   //!
   //! @param compound_id is the 3-letter code for the residue/ligand in the first model, e.g. "TYR" for tyrosine
   //! @param imol_enc is the molecule index for the residue type/compound_id
   //! @param atom_name is the name of the atom, e.g. "OH"
   //!
   //! @return the hb_type for the given atom. On failure return an empty string.
   //! Valid types are: "HB_UNASSIGNED" ,"HB_NEITHER", "HB_DONOR", "HB_ACCEPTOR", "HB_BOTH", "HB_HYDROGEN".
   std::string get_hb_type(const std::string &compound_id, int imol_enc, const std::string &atom_name) const;

   //! Get the GPhL extra restraint information (from the input cif file)
   //!
   //! @param compound_id is the 3-letter code for the residue/ligand in the first model, e.g. "TYR" for tyrosine
   //! @param imol_enc is the molecule index for the residue type/compound_id
   //!
   //! @return a vector/list of string pairs that were part of a gphl_chem_comp_info.
   //! return an empty vector/list on failure to find any such info.
   std::vector<std::pair<std::string, std::string> > get_gphl_chem_comp_info(const std::string &compound_id, int imol_enc);

   //! Get a list of atom names and their associated AceDRG atom types
   //!
   //! @param compound_id is the 3-letter code for the residue/ligand in the first model, e.g. "TYR" for tyrosine
   //! @param imol_enc is the molecule index for the residue type/compound_id
   //!
   //! @return a list of atom names and their associated AceDRG atom types, return an empty list
   //! on failure (e.g. when atoms types are not in the dictionary)
   std::vector<std::pair<std::string, std::string> > get_acedrg_atom_types(const std::string &compound_id, int imol_enc) const;

   //! Get AceDRG atom types for ligand bonds
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the atom selection CID e.g "//A/15" (residue 15 of chain A)
   //!
   //! @return a `coot::acedrg_types_for_residue_t` - which contains a vector/list of bond descriptions.
   coot::acedrg_types_for_residue_t get_acedrg_atom_types_for_ligand(int imol, const std::string &residue_cid) const;

   //! Set the occupancy for the given atom selection
   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)
   //! @param occ_new is the new occupancy
   void set_occupancy(int imol, const std::string &cid, float occ_new);

   //! Get atom selection as json
   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)
   std::string get_molecule_selection_as_json(int imol, const std::string &cid) const;

   //! Write a PNG for the given compound_id.
   //!
   //! Currently this function does nothing (drawing is done with the not-allowed cairo)
   //!
   //! @param compound_id is the 3-letter code for the residue/ligand in the first model, e.g. "TYR" for tyrosine
   //! @param imol is the model molecule index
   void write_png(const std::string &compound_id, int imol, const std::string &file_name) const;

   //! Write the coordinates to the given file name
   //!
   //! @param imol is the model molecule index
   //! @param file_name is the name of the new model
   //!
   //! @return 1 on success and 0 on failure
   int write_coordinates(int imol, const std::string &file_name) const;

   //! Set the state for drawing missing loops
   //!
   //! By default missing loops are drawn. This function allows missing loops to not be
   //! drawn. Sometimes that can clarify the representation. This is a lightweight function
   //! that sets a flag that is used by subsequent calls to get_bonds_mesh()
   //!
   //! @param state is True to mean it is enabled
   void set_draw_missing_residue_loops(bool state);

   //! Get the bonds mesh
   //!
   //! @param mode is "COLOUR-BY-CHAIN-AND-DICTIONARY", "CA+LIGANDS" or "VDW-BALLS"
   //! @param against_a_dark_background allows the bond colours to be relevant for the background.
   //! When the background is dark, the colours should (as a rule) be bright and pastelly.
   //! When the background is light/white, the colour are darker and more saturated.
   //! @param smoothness_factor controls the number of triangles used to make the bond cylinders
   //! and spheres for the atoms - it rises in powers of 4. 1 is the smallest smoothness_factor,
   //! 2 looks nice and 3 is best.
   //! @param bond_width is the bond width in Angstroms. 0.12 is a reasonable default value.
   //! @param atom_radius_to_bond_width_ratio allows the representation of "ball and stick". To do so use a value
   //! between 1.5 and 3.0. The ratio for "liquorice" representation is 1.0.
   //!
   //! @return a `simple_mesh_t`
   coot::simple_mesh_t get_bonds_mesh(int imol, const std::string &mode,
                                      bool against_a_dark_background,
                                      float bond_width, float atom_radius_to_bond_width_ratio,
                                      int smoothness_factor);

   //! Get the instanced bonds mesh.
   //!
   //! @param mode is "COLOUR-BY-CHAIN-AND-DICTIONARY" - more modes to follow
   //! @param against_a_dark_background allows the bond colours to be relevant for the background.
   //! When the background is dark, the colours should (as a rule) be bright and pastelly.
   //! When the background is light/white, the colour are darker and more saturated.
   //! @param bond_width is the bond width in Angstroms. 0.12 is a reasonable default value.
   //! @param atom_radius_to_bond_width_ratio allows the representation of "ball and stick". To do so use a value
   //! between 1.5 and 3.0. The ratio for "liquorice" representation is 1.0.
   //! @param show_atoms_as_aniso_flag if true, if possible, show the atoms with thermal ellipsoids.
   //! @param show_aniso_atoms_as_ortep_flag if true, show any anisotropic atoms with ortep style.
   //! @param draw_hydrogen_atoms_flag if true, bonds to hydrogen atoms should be added.
   //! @param smoothness_factor controls the number of triangles used to make the bond cylinders
   //! and spheres for the atoms - it rises in powers of 4. 1 is the smallest smoothness_factor,
   //! 2 looks nice and 3 is best. Instancing may mean that smoothness factor 3 should
   //! be used by default.
   //! @return a `instanced_mesh_t`
   coot::instanced_mesh_t get_bonds_mesh_instanced(int imol, const std::string &mode,
                                                   bool against_a_dark_background,
                                                   float bond_width, float atom_radius_to_bond_width_ratio,
                                                   bool show_atoms_as_aniso_flag,
                                                   bool show_aniso_atoms_as_ortep_flag,
                                                   bool draw_hydrogen_atoms_flag,
                                                   int smoothness_factor);

   //! As `get_bonds_mesh_instanced` above, but only return the bonds for the atom selection.
   //! Typically one would call this with a wider bond width than one would use for standards atoms (all molecule)
   //!
   //! @param atom_selection_cid e.g. "//A/15" (all the atoms in residue 15 of chain A)
   //! @return a `instanced_mesh_t`
   coot::instanced_mesh_t get_bonds_mesh_for_selection_instanced(int imol, const std::string &atom_selection_cid,
                                                                 const std::string &mode,
                                                                 bool against_a_dark_background,
                                                                 float bond_width, float atom_radius_to_bond_width_ratio,
                                                                 bool show_atoms_as_aniso_flag,
                                                                 bool show_aniso_atoms_as_ortep_flag,
                                                                 bool draw_hydrogen_atoms_flag,
                                                                 int smoothness_factor);

   //! Get the Goodsell style mesh
   //!
   //! @param colour_wheel_rotation_step' the amount, in degrees, that the colour wheel advances between different chains.
   //!        97 degrees is a reasonable starting value
   //! @param saturation' a number between 0 and 1, where 0 is grey and 1 is "lego-like" colour scheme. 0.5 is a nice middle value
   //! @param goodselliness is the degree to which the C-atoms are desaturated, 0.3 is a reasonable value
   //!
   //! @return a `instanced_mesh_t`
   coot::instanced_mesh_t get_goodsell_style_mesh_instanced(int imol, float colour_wheel_rotation_step,
                                                            float saturation, float goodselliness);

   //! Export map molecule as glTF file
   //!
   //! glTF files can be imported into Blender or other 3D graphics applications
   //!
   //! @param imol is the model molecule index
   //! @param x x of the center of the screen
   //! @param y y of the center of the screen
   //! @param z z of the center of the screen
   //! @param radius e.g. 12.0 for X-ray map and 100.0 for Cryo-EM map
   //! @param contour_level e.g. 1.5 sd for X-ray, 5.0 sd for cryo-EM
   //! @param file_name extension should be .glb
   void export_map_molecule_as_gltf(int imol, float pos_x, float pos_y, float pos_z, float radius, float contour_level,
                                    const std::string &file_name);

   //! Export model molecule as glTF file
   //!
   //! glTF files can be imported into Blender or other 3D graphics applications
   //!
   //! Same parameters as the `get_bonds_mesh` function.
   //! `draw_hydrogen_atoms_flag` and `draw_missing_residue_loops` are typically False.
   //! This API will change - we want to specify surfaces and ribbons too.
   void export_model_molecule_as_gltf(int imol,
                                      const std::string &selection_cid,
                                      const std::string &mode,
                                      bool against_a_dark_background,
                                      float bonds_width, float atom_radius_to_bond_width_ratio, int smoothness_factor,
                                      bool draw_hydrogen_atoms_flag, bool draw_missing_residue_loops,
                                      const std::string &file_name);

   //! Export molecular representation as glTF file
   //!
   //! glTF files can be imported into Blender
   //!
   //! @param imol is the model molecule index
   //! @param atom_selection_cid : e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param colour_scheme is one of "colorRampChainsScheme", "colorBySecondaryScheme", "Chain"
   //! @param style "Ribbon" or "MolecularSurface"
   //! @param secondary_structure_usage_flag  0 (USE_HEADER), 1 (DONT_USE) or 2 (CALC_SECONDARY_STRUCTURE)
   //! @param file_name of the glTF (the file will be compressed, so choose ".glb" as the extension)
   void export_molecular_representation_as_gltf(int imol, const std::string &atom_selection_cid,
                                               const std::string &colour_scheme, const std::string &style,
                                               int secondary_structure_usage_flag,
                                               const std::string &file_name);

   //! export chemical features for the specified residue
   //!
   void export_chemical_features_as_gltf(int imol, const std::string &cid,
                                         const std::string &file_name) const;

   //! set the gltf PBR roughness factor
   //!
   //! @param imol is the model molecule index
   //! @param roughness_factor is the factor for the roughness (0.0 to 1.0)
   void set_gltf_pbr_roughness_factor(int imol, float roughness_factor);

   //! set the gltf PBR metalicity factor
   //!
   //! @param imol is the model molecule index
   //! @param metalicity is the factor for the roughness (0.0 to 1.0)
   void set_gltf_pbr_metalicity_factor(int imol, float metalicity);

   //! Get colour table (for testing)
   //!
   //! @param imol is the model molecule index
   //! @param against_a_dark_background allows the bond colours to be relevant for the background.
   //!
   //! @return the colour table
   std::vector<glm::vec4> get_colour_table(int imol, bool against_a_dark_background) const;

   //! Set the colour wheel rotation base for the specified molecule (in degrees)
   //!
   //! @param imol is the model molecule index
   //! @param r is the rotation angle in degrees
   void set_colour_wheel_rotation_base(int imol, float r);

   //! Set the base colour - to be used as a base for colour wheel rotation
   //!
   //! RGB colour codes,
   //! e.g. green is r:0, g: 255, b:0
   //!
   //! @param imol is the model molecule index
   void set_base_colour_for_bonds(int imol, float r, float g, float b);

   //! Add an atom selection cid for atoms and bonds not to be drawn
   //!
   //! @param imol is the model molecule index
   //! @param atom_selection_cid: e.g "//A/15" (all the atoms in residue 15 of chain A)
   void add_to_non_drawn_bonds(int imol, const std::string &atom_selection_cid);

   //! Clear the set of non-drawn atoms (so that they can be displayed again)
   //!
   //! @param imol is the model molecule index
   void clear_non_drawn_bonds(int imol);

   //! Print non-drawn bonds
   //!
   //! @param imol is the model molecule index
   void print_non_drawn_bonds(int imol) const;

   //! User-defined colour-index to colour
   void set_user_defined_bond_colours(int imol, const std::map<unsigned int, std::array<float, 4> > &colour_map);

   //! Set the user-defined residue selections (CIDs) to colour index
   //!
   //! @param imol is the model molecule index
   void set_user_defined_atom_colour_by_selection(int imol, const std::vector<std::pair<std::string, unsigned int> > &indexed_residues_cids,
                                                  bool colour_applies_to_non_carbon_atoms_also);

   //! Add a colour rule for M2T representations
   void add_colour_rule(int imol, const std::string &selection_cid, const std::string &colour);

   //! Add multiple colour rules
   //!
   //! @param selections_and_colours_combo_string e.g. "//A/1^#cc0000|//A/2^#cb0002|//A/3^#c00007",
   //! where "|" is the separator for each rule
   //! and "^" is the separator for the selection string and the colour string
   void add_colour_rules_multi(int imol, const std::string &selections_and_colours_combo_string);

   //! Delete the colour rules for the given molecule
   //!
   //! @param imol is the model molecule index
   void delete_colour_rules(int imol);

   //! Get the colour rules
   //!
   //! @param imol is the model molecule index
   std::vector<std::pair<std::string, std::string> > get_colour_rules(int imol) const;

   //! Print the colour rules
   //!
   //! @param imol is the model molecule index
   void print_colour_rules(int imol) const;

   //! Use bespoke carbon atom colour
   //!
   //! @param imol is the model molecule index
   void set_use_bespoke_carbon_atom_colour(int imol, bool state);

   //! Set bespoke carbon atom colour
   //!
   //! @param imol is the model molecule index
   void set_bespoke_carbon_atom_colour(int imol, const coot::colour_t &col);

   //! Update float parameter for MoleculesToTriangles molecular mesh
   void M2T_updateFloatParameter(int imol, const std::string &param_name, float value);

   //! Update int parameter for MoleculesToTriangles molecular mesh
   void M2T_updateIntParameter(int imol, const std::string &param_name, int value);

   //! Get ribbon and molecular surface representation
   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param colour_scheme should be one of "colorRampChainsScheme", "colorBySecondaryScheme",
   //! "Chain"
   //! @param style "Ribbon" or "MolecularSurface"
   //! @param secondary_structure_usage_flag 0 (USE_HEADER), 1 (DONT_USE) or 2 (CALC_SECONDARY_STRUCTURE).
   //!
   //! @return a `simple_mesh_t`
   coot::simple_mesh_t get_molecular_representation_mesh(int imol, const std::string &cid, const std::string &colour_scheme,
                                                         const std::string &style, int secondary_structure_usage_flag);

   //! Get a Gaussian surface representation
   //!
   //! @param imol is the model molecule index
   //! @param sigma default 4.4
   //! @param contour_level default 4.0
   //! @param box_radius default 5.0
   //! @param grid_scale default 0.7
   //! @param b_factor default 100.0 (use 0.0 for no FFT-B-factor smoothing)
   //!
   //! @return a `simple_mesh_t` composed of a number of Gaussian surfaces (one for each chain)
   coot::simple_mesh_t get_gaussian_surface(int imol, float sigma, float contour_level,
                                            float box_radius, float grid_scale, float b_factor) const;

   //! Get chemical features for the specified residue
   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //!
   //! @return a `simple_mesh_t`
   coot::simple_mesh_t get_chemical_features_mesh(int imol, const std::string &cid) const;

   //! get an (mmdb-style) atom
   //!
   //! If more than one atom is selected by the selection cid, then the first
   //! atom is returned.
   //!
   //! Don't use this in emscript.
   //!
   //! @param imol is the model molecule index
   //! @param cid is the coordinate-id for the atom.
   //! @returns either the specified atom or nullopt (None) if not found
   mmdb::Atom *get_atom_using_cid(int imol, const std::string &cid) const;

   //! get an (mmdb-style) residue
   //!
   //! If more than one residue is selected by the selection cid, then the first
   //! residue is returned.
   //!
   //! Don't use this in emscript.
   //!
   //! @param imol is the model molecule index
   //! @param cid is the coordinate-id for the residue
   //! @returns either the specified residue or nullopt (None) if not found
   mmdb::Residue *get_residue_using_cid(int imol, const std::string &cid) const;

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! get atom - internal (C++) usage only
   //!
   //! @returns either the specified atom or null if not found - don't use this in emscript
   mmdb::Atom *get_atom(int imol, const coot::atom_spec_t &atom_spec) const;
   //! get residue - internal (C++) usage only
   //!
   //! @returns either the specified residue or null if not found - don't use this in emscript
   mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &residue_spec) const;
   //! get the atom position - don't use this in emscript
   std::pair<bool, coot::Cartesian> get_atom_position(int imol, coot::atom_spec_t &atom_spec);
#endif

   //! Residue is nucleic acid?
   //!
   //! Every residue in the selection is checked
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //!
   //! @return a bool
   bool residue_is_nucleic_acid(int imol, const std::string &cid) const;

   //! Get the residue CA position
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //!
   //! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
   std::vector<double> get_residue_CA_position(int imol, const std::string &cid) const;

   //! Get the average residue position
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //!
   //! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
   std::vector<double> get_residue_average_position(int imol, const std::string &cid) const;

   //! Get the average residue side-chain position
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //!
   //! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
   std::vector<double> get_residue_sidechain_average_position(int imol, const std::string &cid) const;

   //! Get number of atoms
   //!
   //! @param imol is the model molecule index
   //!
   //! @return the number of atoms in the specified model, or 0 on error
   unsigned int get_number_of_atoms(int imol) const;

   //! Get molecule diameter
   //!
   //! @param imol is the model molecule index
   //!
   //! @return an estimate of the diameter of the model molecule (-1 on failure)
   float get_molecule_diameter(int imol) const;

   //! Get number of hydrogen atoms
   //!
   //! @param imol is the model molecule index
   //!
   //! @return the number of hydrogen atoms in the specified model, or -1 on error
   int get_number_of_hydrogen_atoms(int imol) const;

   //! Get the chain IDs in the given molecule
   //!
   //! @param imol is the model molecule index
   //!
   //! @return vector/list of chain-ids for the given molecule
   std::vector<std::string> get_chains_in_model(int imol) const;

   //! Get the chains that are related by NCS or molecular symmetry
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a vector/list of vector/list of chain ids, e.g. [[A,C], [B,D]] (for hemoglobin).
   std::vector<std::vector<std::string> > get_ncs_related_chains(int imol) const;

   //! Get the single letter codes for the residues in the specified chain
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //!
   //! @return vector/list of single letter codes - in a pair with the given residue spec
   std::vector<std::pair<coot::residue_spec_t, std::string> > get_single_letter_codes_for_chain(int imol, const std::string &chain_id) const;

   //! Get a list of residues that don't have a dictionary
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a list of residue names (compound_ids) for which there is no dictionary in the geometry store
   std::vector<std::string> get_residue_names_with_no_dictionary(int imol) const;

   //! Get residue name
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //!
   //! @return the residue name, return a blank string on residue not found.
   std::string get_residue_name(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) const;

   //! Get the residue type
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/16" (residue 16 of chain A)
   //! @return a string. Return an empty string on failure
   std::string get_residue_type(int imol, const std::string &cid) const;

   //! Get the SMILES string for the give residue type
   //!
   //! @param residue 3 letter-code/name of the compound-id
   //! @param imol_enc is the molecule index for the residue type/compound_id
   //!
   //! @return the SMILES string if the residue type can be found in the dictionary store
   //! or the empty string on a failure.
   std::string get_SMILES_for_residue_type(const std::string &residue_name, int imol_enc) const;

   //! Get residues with missing atoms
   //!
   //! @param imol is the model molecule index
   //!
   //! @return an object that has information about residues without dictionaries and residues with missing atom
   //! in the the specified molecule
   std::vector<coot::residue_spec_t> residues_with_missing_atoms(int imol);

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! Ths function is not const because missing_atoms() takes a non-const pointer to the geometry
   // (20230117-PE I should fix that)
   //!
   //! @param imol is the model molecule index
   //!
   //! @return an object that has information about residues without dictionaries and residues with missing atom
   //! in the the specified molecule
   coot::util::missing_atom_info missing_atoms_info_raw(int imol);
#endif

   //! Get missing residue ranges
   //!
   //! @param imol is the model molecule index
   //! @return missing residue ranges
   std::vector<coot::residue_range_t> get_missing_residue_ranges(int imol) const;

   //! Get a list of residues specs that have atoms within distance of the atoms of the specified residue
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param dist is the distance in Angstrom
   //!
   //! @return a list of residue specs
   std::vector<coot::residue_spec_t> get_residues_near_residue(int imol, const std::string &residue_cid, float dist) const;

   //! Get atom distances
   //!
   //! @param imol is the model molecule index
   //! @param cid_res_1 is the first atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)
   //! @param cid_res_2 is the second atom selection CID e.g "//A/17/NH" (atom NH in residue 17 of chain A)
   //! @param dist is the distance in Angstrom
  std::vector<coot::atom_distance_t>
  get_distances_between_atoms_of_residues(int imol, const std::string &cid_res_1, const std::string &cid_res_2,
					  float dist_max) const;

   //! Superposition (using SSM)
   //!
   //! The specified chain of the moving molecule is superposed onto the chain in the reference molecule (if possible).
   //!
   //! @param imol_ref the reference model molecule index
   //! @param chain_id_ref the chain ID for the reference chain
   //! @param imol_mov the moving model molecule index
   //! @param chain_id_mov the chain ID for the moving chain
   superpose_results_t SSM_superpose(int imol_ref, const std::string &chain_id_ref,
                                     int imol_mov, const std::string &chain_id_mov);

   //! Superpose using LSQ - setup the matches
   //!
   //! @param chain_id_ref the chain ID for the reference chain
   //! @param res_no_ref_start the starting residue number in the reference chain
   //! @param res_no_ref_end the ending residue number in the reference chain
   //! @param chain_id_mov the chain ID for the moving chain
   //! @param res_no_mov_start the starting residue number in the moving chain
   //! @param res_no_mov_end the ending residue number in the moving chain
   //! @param match_type 0: all, 1: main, 2: CAs, 3: N, CA, C, 4: N, CA, CB, C
   void add_lsq_superpose_match(const std::string &chain_id_ref, int res_no_ref_start, int res_no_ref_end,
                                const std::string &chain_id_mov, int res_no_mov_start, int res_no_mov_end,
                                int match_type);

   //! Superpose using LSQ for a scpecific atom - setup the matches
   //!
   //! @param chain_id_ref the chain ID for the reference chain
   //! @param res_no_ref the residue number in the reference chain
   //! @param atom_name_ref the name of the reference atom
   //! @param chain_id_mov the chain ID for the moving chain
   //! @param res_no_mov the residue number in the moving chain
   //! @param atom_name_mov the name of the moving atom
   void add_lsq_superpose_atom_match(const std::string &chain_id_ref, int res_no_ref, const std::string &atom_name_ref,
                                     const std::string &chain_id_mov, int res_no_mov, const std::string &atom_name_mov);

   //! Clear any existing lsq matchers
   void clear_lsq_matches();

   std::vector<coot::lsq_range_match_info_t> lsq_matchers;

   //! Apply the superposition using LSQ
   //!
   //! @param imol_ref the reference model molecule index
   //! @param imol_mov the moving model molecule index
   //! @return the success status, i.e. whether or not there were enough atoms to superpose
   bool lsq_superpose(int imol_ref, int imol_mov);

   //! Transform a map and create a new map
   //!
   //! @param imol_map map molecule index
   //! @param lsq_matrix is an object of type lsq_results_t, is the object returned by `get_lsq_matrix()`
   //! @param x is the point in the map about which the map is transformed
   //! @param y is the point in the map about which the map is transformed
   //! @param z is the point in the map about which the map is transformed
   //! @param radius the radius of the transformed map, typically between 10 and 100 A
   //!
   //! @return the molecule index of the new map, -1 for failure
   int transform_map_using_lsq_matrix(int imol_map, lsq_results_t lsq_matrix, float x, float y, float z, float radius);
   //! Get LSQ matrix
   //!
   //! don't apply it to the coordinates
   //!
   //! @param imol_ref the reference model molecule index
   //! @param imol_mov the moving model molecule index
   //! @param summary_to_screen if True, write a summary of the statistics to the output
   //!
   //! @return the transformation matrix in a simple class
   lsq_results_t get_lsq_matrix(int imol_ref, int imol_mov, bool summary_to_screen) const;

   //! Get symmetry
   //!
   //! now comes in a simple container that also includes the cell
   coot::symmetry_info_t
   get_symmetry(int imol, float symmetry_search_radius, float centre_x, float centre_y, float centre_z) const;

   //! Get the cell
   //!
   //! Check that `is_set` is true before use.
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a `cell_t`
   ::api::cell_t get_cell(int imol)  const;

   //! Get the middle of the "molecule blob" in cryo-EM reconstruction maps
   //!
   //! @param imol is the map molecule index
   //!
   //! @return a `coot::util::map_molecule_centre_info_t`
   coot::util::map_molecule_centre_info_t get_map_molecule_centre(int imol) const;

   //! Undo
   //!
   //! @param imol is the model molecule index
   //!
   //! @return 1 on successful undo, return 0 on failure
   int undo(int imol);

   //! Redo
   //!
   //! @param imol is the model molecule index
   //!
   //! @return 1 on successful redo, return 0 on failure
   int redo(int imol);

   //! Get the torsion of the specified atom in the specified residue
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID, e.g. //A/15 (residue 15 in chain A)
   //! @param atom_names is a list of atom names, e.g. [" CA ", " CB ", " CG ", " CD1"]
   //!
   //! @return a pair, the first of which is a succes status (1 success, 0 failure), the second is the torsion in degrees
   std::pair<int, double> get_torsion(int imol, const std::string &cid, const std::vector<std::string> &atom_names);

   //! Change the B factors
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID, e.g. //A/15 (residue 15 in chain A)
   //! @param temp_fact is the isotropic ADP/temperature factor, e.g.,  22
   void set_temperature_factors_using_cid(int imol, const std::string &cid, float temp_fact);

   // -------------------------------- map utils -------------------------------------------
   //! \name Map Utils

#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! @return the map sampling rate (default is 1.8)
   float map_sampling_rate;
#endif

   //! Get map sampling rate
   //!
   //! @return the map sampling rate, the default is 1.8
   float get_map_sampling_rate() { return map_sampling_rate; }

   //! Set the map sampling rate
   //!
   //! Higher numbers mean smoother maps, but they take
   //! longer to generate, longer to transfer, longer to parse and longer to draw
   //!
   //! @param msr is the map sampling rate to set, the default is 1.8
   void set_map_sampling_rate(float msr) { map_sampling_rate = msr; }

   //! Read the given mtz file
   //!
   //! @param file_name is the name of the MTZ file
   //! @param f F column, "FWT"
   //! @param phi phi column, "PHWT"
   //! @param weight weight column, "W"
   //! @param use_weight flag for weights usage, False
   //! @param is_a_difference_map False
   //!
   //! @return the new molecule number or -1 on failure
   int read_mtz(const std::string &file_name, const std::string &f, const std::string &phi, const std::string &weight,
                bool use_weight, bool is_a_difference_map);

   //! Replace map by mtz
   //!
   //! @param imol is the map molecule index
   //! @param f F column, "FWT"
   //! @param phi phi column, "PHWT"
   //! @param weight weight column, "W"
   //! @param use_weight flag for weights usage, False
   int replace_map_by_mtz_from_file(int imol, const std::string &file_name, const std::string &f, const std::string &phi,
                                    const std::string &weight, bool use_weight);

   //! class for the information about columns extracted from auto-reading the given mtz file
   class auto_read_mtz_info_t {
   public:
      //! molecule index
      int idx;
      //! F column
      std::string F;
      //! phi column
      std::string phi;
      //! weights column
      std::string w;
      //! flag for weights usage
      bool weights_used;
      //! F_obs column. There were not available if the return value is empty
      std::string F_obs;
      //! sigF_obs column
      std::string sigF_obs;
      //! R-Free column. There were not available if the return value is empty
      std::string Rfree;
      // constructor
      auto_read_mtz_info_t() {idx = -1; weights_used = false; }
      // constructor
      auto_read_mtz_info_t(int index, const std::string &F_in, const std::string &phi_in) :
         idx(index), F(F_in), phi(phi_in), weights_used(false) {}
      // set Fobs sigFobs column labels
      void set_fobs_sigfobs(const std::string &f, const std::string &s) {
         F_obs = f;
         sigF_obs = s;
      }
   };

   //! Auto read the given mtz file
   //!
   //! @param file_name is the name of the MTZ file
   //!
   //! @return a vector of the maps created from reading the file
   std::vector<auto_read_mtz_info_t> auto_read_mtz(const std::string &file_name);

   //! Read a CCP4 (or MRC) map
   //!
   //! There is currently a size limit of 1000 pixels per edge.
   //!
   //! @param file_name is the name of the map file
   //! @param is_a_difference_map False
   //!
   //! @return the new molecule number or -1 on failure
   int read_ccp4_map(const std::string &file_name, bool is_a_difference_map);

   //! Write a map
   //!
   //! @param imol is the map molecule index
   //! @param file_name is the name of the new map file
   //!
   //! @return 1 on a successful write, return 0 on failure.
   int write_map(int imol, const std::string &file_name) const;

   //! Get map mean
   //!
   //! @param imol is the map molecule index
   //!
   //! @return the mean of the map or -1 if imol is not a valid map molecule index
   float get_map_mean(int imol) const;

   //! Get map rmsd approx
   //!
   //! the function is approximate because the epsilon factor is not taken into account
   //!
   //! @param imol_map is the map molecule index
   //!
   //! @return the map rmsd. -1 is returned if imol_map is not a valid map molecule index.
   float get_map_rmsd_approx(int imol_map) const;

   //! Get map histogram
   //!
   //! @param imol is the map molecule index
   //! @param n_bins is the number of bins - 200 is a reasonable default.
   //! @param zoom_factor (reduces the range by the given factor)
   //! centred around the median (typically 1.0 but usefully can vary until ~20.0).
   //!
   //! @return the map histogram
   coot::molecule_t::histogram_info_t get_map_histogram(int imol, unsigned int n_bins, float zoom_factor) const;


   //! @param imol is the map molecule index
   //!
   //! @return the suggested initial contour level. Return -1 on not-a-map
   float get_suggested_initial_contour_level(int imol) const;

   //! Check if a map is an EM map or not
   //!
   //! @param imol is the map molecule index
   //!
   //! @return the "EM" status of this molecule. Return false on not-a-map.
   bool is_EM_map(int imol) const;

   //! Create a new map that is blurred/sharpened
   //!
   //! @param imol is the map molecule index
   //! @param b_factor e.g. 100.0
   //! @param in_place_flag True if you want to replace the current map, False if you want to create a new map
   //!
   //! @return the molecule index of the new map or -1 on failure or if `in_place_flag` was true.
   int sharpen_blur_map(int imol_map, float b_factor, bool in_place_flag);

   //! Create a new map that is blurred/sharpened and resampling
   //!
   //! Note that resampling can be slow, a resample_factor of 1.5 is about the limit of the trade of of prettiness for speed.
   //!
   //! @param imol_map is the map molecule index
   //! @param b_factor e.g. 100.0
   //! @param resample_factor e.g. 1.4
   //! @param in_place_flag True if you want to replace the current map, False if you want to create a new map
   //!
   //! @return the molecule index of the new map or -1 on failure or if `in_place_flag` was true.
   int sharpen_blur_map_with_resample(int imol_map, float b_factor, float resample_factor, bool in_place_flag);

   //! Mask map by atom selection
   //!
   //! (note the argument order is reversed compared to the coot api).
   //!
   //! @param imol_coords is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param atom_radius is the atom radius. Use a negative number to mean "default".
   //! @param invert_flag changes the parts of the map that are masked, so to highlight the density
   //! for a ligand one would pass the cid for the ligand and invert_flag as true, so that the
   //! parts of the map that are not the ligand are suppressed.
   //!
   //! @return the index of the new map or -1 on failure
   int mask_map_by_atom_selection(int imol_coords, int imol_map, const std::string &cid, float atom_radius, bool invert_flag);

   //! Partition the input map
   //!
   //! Each voxel in the map is assigned to the chain
   //! to which it is nearest. Unlike masking, the generated maps are not restricted to be
   //! "close" to the atoms in the atom selection.
   //!
   //! e.g. maskChains for ChimeraX - JiangLab
   //!
   //! @param imol_map is the map molecule index
   //! @param imol_model is the model molecule index
   //!
   //! @return a vector/list of the molecules indices of the newly created maps
   std::vector<int> partition_map_by_chain(int imol_map, int imol_model);

   //! Make a masked map
   //!
   //! @param imol_map_ref is the map molecule index
   //! @param imol_model is the model molecule index
   //! @param atom_selection_cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param radius is the radius of the map, e.g. 12.0 for X-ray map and 100.0 for Cryo-EM map
   //!
   //! @return the index of the newly created mask. Return -1 on failure.
   int make_mask(int imol_map_ref, int imol_model, const std::string &atom_selection_cid, float radius);

   //! Generate a new map which is the hand-flipped version of the input map
   //!
   //! @param imol_map is the map molecule index
   //!
   //! @return the molecule index of the new map, or -1 on failure.
   int flip_hand(int imol_map);

   //! Make a vector/list of maps that are split by chain-id of the input imol
   //!
   //! @param imol is the model molecule index
   //! @param imol_map is the map molecule index
   //!
   //! @return a vector/list of the map molecule indices.
   std::vector<int> make_masked_maps_split_by_chain(int imol, int imol_map);

   //! dedust map
   //!
   //! @param imol_map the map molecule index
   //!
   //! @return the map molecule index of the dedusted map or -1 on failure
   int dedust_map(int imol);

   //! Set the map colour
   //!
   //! The next time a map mesh is requested, it will have this colour.
   //! This does not apply to/affect the colour of the difference maps.
   //!
   //! RGB colour codes,
   //! e.g. green is r:0, g: 255, b:0
   //!
   //! @param imol is the model molecule index
   //!
   void set_map_colour(int imol, float r, float g, float b);

   //! Set the state of the mode of the threading in map contouring
   //!
   //! @param state is True to mean that it is enabled. The default is True.
   void set_map_is_contoured_with_thread_pool(bool state);

   //! Get the mesh for the map contours
   //!
   //! This function is not const because the internal state of a `coot_molecule_t` is changed.
   //!
   //! @param imol is the model molecule index
   //! @param position_x is the x coordinate of the target position
   //! @param position_y is the y coordinate of the target position
   //! @param position_z is the z coordinate of the target position
   //! @param radius is the radius of the map, e.g. 12.0 for X-ray map and 100.0 for Cryo-EM map
   //! @param contour_level e.g. 1.5 sd for X-ray, 5.0 sd for cryo-EM
   //!
   //! @return a `simple_mesh_t` for the map contours of the specified map
   coot::simple_mesh_t get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level);

   //! Get the mesh for the map contours using another map for colouring
   //!
   //! @param imol_ref is the reference map index
   //! @param imol_map_for_coloring is the target map index
   //! @param position_x is the x coordinate of the target position
   //! @param position_y is the y coordinate of the target position
   //! @param position_z is the z coordinate of the target position
   //! @param radius is the radius of the map, e.g. 12.0 for X-ray map and 100.0 for Cryo-EM map
   //! @param contour_level e.g. 1.5 sd for X-ray, 5.0 sd for cryo-EM
   //! @param other_map_for_colouring_min_value e.g. -1.0 in the case of correlation map
   //! @param other_map_for_colouring_max_value e.g. 1.0 in the case of correlation map
   //! @param invert_colour_ramp e.g. red to blue rather than blue to red
   //!
   //! @return a `simple_mesh_t` for the map contours of the specified map
   coot::simple_mesh_t get_map_contours_mesh_using_other_map_for_colours(int imol_ref, int imol_map_for_colouring,
                                                                         double position_x, double position_y, double position_z,
                                                                         float radius, float contour_level,
                                                                         float other_map_for_colouring_min_value,
                                                                         float other_map_for_colouring_max_value,
                                                                         bool invert_colour_ramp);
   //! Set the map saturation
   //!
   //! @param s is the map saturation, e.g. a number between 0 and 1, where 0 is grey and 1 is "lego-like" colour scheme.
   //!        0.5 is a nice middle value
   void set_map_colour_saturation(int imol, float s);
   void set_colour_map_for_map_coloured_by_other_map(std::vector<std::pair<double, std::vector<double> > > colour_table );

   user_defined_colour_table_t colour_map_by_other_map_user_defined_table;

   //! Get map vertices histogram
   //!
   //! Note not const because get_map_contours_mesh() is not const
   //!
   //! @param imol is the map molecule index
   //! @param n_bins is the number of bins - 40 is a reasonable default.
   //!
   //! @return the map vertices histogram
   coot::molecule_t::histogram_info_t get_map_vertices_histogram(int imol, int imol_map_for_sampling,
								 double position_x, double position_y, double position_z,
								 float radius, float contour_level,
								 unsigned int n_bins);



   //! Get the latest sfcalc stats
   //!
   //! @return a sfcalc_genmap_stats_t object
   coot::util::sfcalc_genmap_stats_t get_latest_sfcalc_stats() const { return latest_sfcalc_stats; }

   class r_factor_stats {
      public:
      float r_factor;  // 0 to 1
      float free_r_factor;
      int rail_points_total;
      int rail_points_new;
   };

   //! Get the R-factors
   //!
   //! @return an object `r_factor_stats`
   r_factor_stats get_r_factor_stats();

   //! Get the R factor stats as a string
   //!
   //! @param rfs is the name of the string
   //!
   //! @return a string with the R-factor stats
   std::string r_factor_stats_as_string(const r_factor_stats &rfs) const;

   //! Get the average map
   //!
   //! This function does no normalization of the scales,
   //! presuming that they are pre-normalized.
   //!
   //! @param imol_maps is a colon-separated list of map indices e.g. "2:3:4"
   //! @param scales is the list of weights corresponding to the list of maps.
   //! The number of scales factors should match the number of maps
   //!
   //! @return the index of the new map, or -1 on failure.
   int average_map(const std::string &imol_maps, std::vector<float> &scales);

   //! This function does no normalisation of the scales, presuming that they are pre-normalized.
   //!
   //! The number of maps in imol_maps should match the size of the scale vector. If not, nothing
   //! will happen and False will be returned
   //!
   //! @param imol_map is the map molecule index
   //! @param imol_maps is a colon-separated list of map indices e.g. "2:3:4"
   //! @param scales is the list of weights corresponding to the list of maps
   //!
   //! @return the success status
   bool regen_map(int imol_map, const std::string &imol_maps, const std::vector<float> &scales);

   // -------------------------------- coordinates modelling -------------------------------
   //! \name Coordinates Modelling

   //! Auto-fit rotamer
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //! @param imol_map is the map molecule index
   //!
   //! @return 1 on successful modification, return 0 on failure
   int auto_fit_rotamer(int imol, const std::string &chain_id, int res_no, const std::string &ins_code, const std::string &alt_conf,
                        int imol_map);

   //! Change to the next rotamer (rotamer cycling is implicit if needed)
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //!
   //! @return the change information.
   coot::molecule_t::rotamer_change_info_t change_to_next_rotamer(int imol, const std::string &residue_cid, const std::string &alt_conf);

   //! Change to the previous rotamer (rotamer cycling is implicit if needed)
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //!
   //! @return the change information.
   coot::molecule_t::rotamer_change_info_t change_to_previous_rotamer(int imol, const std::string &residue_cid, const std::string &alt_conf);

   //! Change to the first (0th) rotamer
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //!
   //! @return the change information.
   coot::molecule_t::rotamer_change_info_t change_to_first_rotamer(int imol, const std::string &residue_cid, const std::string &alt_conf);

   //! Change to the nth rotamer
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //!
   //! @return the state of the change.
   int set_residue_to_rotamer_number(int imol, const std::string &residue_cid, const std::string &alt_conf, int rotamer_number);

   //! Delete item
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //! @param scope is one of the strings: ["ATOM", "WATER", "RESIDUE"," CHAIN"," MOLECULE", "LITERAL"]
   //!
   //! @return 1 on successful deletion, return 0 on failure
   std::pair<int, unsigned int> delete_using_cid(int imol, const std::string &cid, const std::string &scope);

   //! Delete atom
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //! @param atom_name is the name of the atom, e.g. "OH"
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_atom(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                   const std::string &atom_name, const std::string &alt_conf);

   //! Delete atom using cid
   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_atom_using_cid(int imol, const std::string &cid);

   //! Delete residue
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   //! Delete residue using cid
   //!
   //! @param imol is the model molecule index
   //! @param cid is the residue selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_residue_using_cid(int imol, const std::string &cid);

   //! Delete residue atoms using alt_conf
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_residue_atoms_with_alt_conf(int imol, const std::string &chain_id, int res_no,
                                                                   const std::string &ins_code, const std::string &alt_conf);
   //! Delete residue atoms using cid
   //!
   //! This is the same as `delete_atom_using_cid`. It will be deleted in the future
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_residue_atoms_using_cid(int imol, const std::string &cid);

   //! Delete side chain
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   //! Delete side chain using cid
   //!
   //! @param imol is the model molecule index
   //! @param cid is the residue selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_side_chain_using_cid(int imol, const std::string &cid);

   //! Delete chain using chain cid
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A" (chain A), "//*" (all chains)
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_chain_using_cid(int imol, const std::string &cid);

   //! Delete the atoms specified in the cid selection
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15/OH" (atom OH in residue 15)
   //!
   //! @return 1 on successful deletion, return 0 on failure to delete.
   std::pair<int, unsigned int> delete_literal_using_cid(int imol, const std::string &cid);

   //! delete all carbohydrate
   //!
   //! @param imol is the model molecule index
   //!
   //! @return true on successful deletion, return false on no deletion.
   bool delete_all_carbohydrate(int imol);

   // (I should have) change(d) that stupid (alt) loc (I should have made you leave your key)
   //
   //! Change alternate conformation
   //!
   //! Note that this function only deals with (swaps) alt confs "A" and "B" - any
   //! alt-conf other than that is ignored.
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 in chain A)
   //! @param change_mode is either "residue", "main-chain", "side-chain" or a comma-separated atom-name
   //! pairs (e.g "N,CA") - you can (of course) specify just one atom, e.g.: "N".
   //!
   //! @return the success status (1 is done, 0 means failed to do)
   int change_alt_locs(int imol, const std::string &cid, const std::string &change_mode);

   //! Add a residue onto the end of the chain by fitting to density
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //!
   //! @return first: 1 on success, second is failure message
   std::pair<int, std::string> add_terminal_residue_directly(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   // std::pair<int, std::string> add_terminal_residue_directly_using_cid(int imol, const std::string &cid);
   //
   //! Add a residue onto the end of the chain by fitting to density using cid
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15/OH" (atom OH in residue 15)
   //! @return success status (1 for good, 0 for not done)
   int add_terminal_residue_directly_using_cid(int imol, const std::string &cid);

   //! Add a residue onto the end of the chain by fitting to density using Buccaneer building and cid
   //!
   //! This function has been removed - is is now a noop.
   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15)
   int add_terminal_residue_directly_using_bucca_ml_growing_using_cid(int imol, const std::string &cid);

   //! Add a residue onto the end of the chain by fitting to density using Buccaneer building
   //!
   //! @param imol is the model molecule index
   //! @param spec is the residue specifier, residue_spec_t("A", 10, "")
   int add_terminal_residue_directly_using_bucca_ml_growing(int imol, const coot::residue_spec_t &spec);

   //! Parameter for `add_waters` function
   //!
   //! @param d is the min distance, default  2.4
   void set_add_waters_water_to_protein_distance_lim_min(float d) {
      ligand_water_to_protein_distance_lim_min = d;
   }

   //! Parameter for `add_waters` function
   //!
   //! @param d is the max distance, default  3.4
   void set_add_waters_water_to_protein_distance_lim_max(float d) {
      ligand_water_to_protein_distance_lim_max = d;
   }

   //! Parameter for `add_waters` function
   //!
   //! @param d is the variance limit, default is 0.1
   void set_add_waters_variance_limit(float d) {
      ligand_water_variance_limit = d;
   }

   //! Parameter for `add_waters` function
   //!
   //! @param d is the sigma cutoff, default is 1.75
   void set_add_waters_sigma_cutoff(float d) {
      ligand_water_sigma_cut_off = d;
   }

   //! Add waters
   //!
   //! @param imol is the model molecule index
   //! @param imol_map is the map molecule index
   //!
   //! @return the number of waters added on a success, -1 on failure.
   int add_waters(int imol_model, int imol_map);

   //! Flood with dummy atoms
   //!
   //! @param imol is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param n_rmsd e.g., 4.0
   //!
   //! @return the number of waters added on a success, -1 on failure.
   int flood(int imol_model, int imol_map, float n_rmsd);

   //! Add hydrogen atoms
   //!
   //! @param imol_model is the model molecule index
   //!
   //! @return 1 on success, 0 on failure.
   int add_hydrogen_atoms(int imol_model);

   //! Delete hydrogen atoms
   //!
   //! @param imol_model is the model molecule index
   //!
   //! @return 1 on a successful deletion, 0 on failure.
   int delete_hydrogen_atoms(int imol_model);

   //! Add an alternative conformation for the specified residue
   //!
   //! @param imol_model is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 in chain A)
   //!
   //! @return 1 on a successful addition, 0 on failure.
   int add_alternative_conformation(int imol_model, const std::string &cid);

   //! Fill the specified residue
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //!
   //! @return 1 on a successful fill, 0 on failure.
   int fill_partial_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   //! Fill the specified residue using cid
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 in chain A)
   //!
   //! @return 1 on a successful fill, 0 on failure.
   int fill_partial_residue_using_cid(int imol, const std::string &cid);

   //! Fill all the the partially-filled residues in the molecule
   //!
   //! @param imol is the model molecule index
   //!
   //! @return 1 on a successful fill, 0 on failure.
   int fill_partial_residues(int imol);

   //! Add N-linked glycosylation
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param glycosylation_name is the type of glycosylation, one of:
   //!       "NAG-NAG-BMA" or "high-mannose" or "hybrid" or "mammalian-biantennary" or "plant-biantennary"
   //! @param asn_chain_id is the chain-id of the ASN to which the carbohydrate is to be added
   //! @param asn_res_no is the residue number of the ASN to which the carbohydrate is to be added
   void add_named_glyco_tree(int imol_model, int imol_map, const std::string &glycosylation_name,
                             const std::string &asn_chain_id, int asn_res_no);

#if NB_VERSION_MAJOR
#else
   //! Flip peptide
   //!
   //! @param imol is the model molecule index
   //! @param atom_spec is the atom specifier, atom_spec_t("A", 10, "", " CA ", "")
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //!
   //! @return 1 on a successful flip
   int flip_peptide(int imol, const coot::atom_spec_t &atom_spec, const std::string &alt_conf);
#endif

   //! Flip peptide using cid
   //!
   //! @param imol is the model molecule index
   //! @param atom_cid is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //!
   //! @return 1 on a successful flip
   int flip_peptide_using_cid(int imol, const std::string &atom_cid, const std::string &alt_conf);

   //! Eigen-flip the specified ligand
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   void eigen_flip_ligand(int imol, const std::string &chain_id, int res_no, const std::string &ins_code);

   //! Eigen-flip ligand using cid
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the residue selection CID e.g "//A/15" (residue 15 of chain A)
   void eigen_flip_ligand_using_cid(int imol, const std::string &residue_cid);

   //! Mutate residue
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the residue selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param new_residue_type is the 3-letter code of the new residue, e.g. "TYR" for tyrosine
   //!
   //! @return 1 on a successful move, 0 on failure.
   int mutate(int imol, const std::string &cid, const std::string &new_residue_type);

   //! Rotate last chi angle of the side chain by 180 degrees
   //!
   //! @param imol is the model molecule index
   //! @param atom_cid is the atom selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A)
   //!
   //! @return 1 on a successful move, 0 on failure.
   int side_chain_180(int imol, const std::string &atom_cid);

   //! JED-Flip the ligand (or residue) at the specified atom
   //!
   //! @param imol is the model molecule index
   //! @param atom_cid is the residue selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param invert_selection is True if you want to use the larger fragment
   //!
   //! @return a non-blank message if there is a problem
   std::string jed_flip(int imol, const std::string &atom_cid, bool invert_selection);

   //! Move the molecule to the given centre
   //!
   //! @param imol is the model molecule index
   //! @param x is the x coordinate of the new centre of the screen
   //! @param y is the y coordinate of the new centre of the screen
   //! @param z is the z coordinate of the new centre of the screen
   //!
   //! @return 1 on a successful move, 0 on failure.
   int move_molecule_to_new_centre(int imol, float x, float y, float z);

   //! Interactive B-factor refinement
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param factor might typically be 0.9 or 1.1
   void multiply_residue_temperature_factors(int imol, const std::string &cid, float factor);

   //! Get molecule centre
   //!
   //! @param imol is the model molecule index
   //!
   //! @return the molecule centre
   coot::Cartesian get_molecule_centre(int imol) const;

   //! Get Radius of Gyration
   //!
   //! @param imol is the model molecule index
   //!
   //! @return the molecule centre. If the number is less than zero, there
   //! was a problem finding the molecule or atoms.
   double get_radius_of_gyration(int imol) const;

   //! Copy the molecule
   //!
   //! @param imol the specified molecule
   //!
   //! @return the new molecule number
   int copy_molecule(int imol);

   //! Copy a fragment given the multi_cid selection string
   //!
   //! @param imol is the model molecule index
   //! @param multi_cids is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||/B/56-66"
   //!
   //! @return the new molecule number (or -1 on no atoms selected)
   int copy_fragment_using_cid(int imol, const std::string &multi_cid);

   //! Copy a fragment given the multi_cid selection string for refinement
   //!
   //! Use this in preference to `copy_fragment_using_cid` when copying
   //! a molecule fragment to make a molten zone for refinement.
   //! That is because this version quietly also copies the residues near the residues of the selection,
   //! so that those residues can be used for links and non-bonded contact restraints.
   //!
   //! @param imol is the model molecule index
   //! @param multi_cids is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||//B/56-66"
   //!
   //! @return the new molecule number (or -1 on no atoms selected)
   int copy_fragment_for_refinement_using_cid(int imol, const std::string &multi_cid);

   //! Copy a residue-range fragment
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A"
   //! @param res_no_start the starting residue number
   //! @param res_no_ref_end the ending residue number
   //!
   //! @return the new molecule number (or -1 on no atoms selected)
   int copy_fragment_using_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end);

   //! Apply transformation to atom selection in the given molecule
   //!
   //! @return the number of atoms moved.
   int apply_transformation_to_atom_selection(int imol, const std::string &atoms_selection_cid,
                                              int n_atoms, // for validation of the atom selection, (int because mmdb atom type)
                                              float m00, float m01, float m02,
                                              float m10, float m11, float m12,
                                              float m20, float m21, float m22,
                                              float c0, float c1, float c2, // the centre of the rotation
                                              float t0, float t1, float t2); // translation

   //! Update the positions of the atoms in the residue
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the residue selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param moved_atoms is a list of the atoms moved in the specified residue, e.g. moved_atom_t(" CA ", 1, 2, 3)
   int new_positions_for_residue_atoms(int imol, const std::string &residue_cid, std::vector<coot::api::moved_atom_t> &moved_atoms);

   //! Update the positions of the atoms in the residues
   //!
   //! @param imol is the model molecule index
   //! @param moved_residue is a list of the residues with moved atoms, e.g. moved_residue_t("A", 10, "")
   int new_positions_for_atoms_in_residues(int imol, const std::vector<coot::api::moved_residue_t> &moved_residues);

   //! Merge molecules
   //!
   //! @param imol is the model molecule index
   //! @param list_of_other_molecules is a colon-separated list of molecules, e.g. "2:3:4"
   //!
   //! @return the first is a flag set to 1 if a merge occurred (and 0 if it did not)
   //! the second is a vector of merge results, i.e. if you merged a ligand, what is the new
   //! residue spec of the ligand, and if you merged a (polymer) chain, what is the new chain-id of
   //! that chain.
   std::pair<int, std::vector<merge_molecule_results_info_t> >
   merge_molecules(int imol, const std::string &list_of_other_molecules);

   // this is called by the above function and is useful for other non-api functions (such as add_compound()).

#ifdef SKIP_FOR_PYTHON_DOXYGEN
//#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   std::pair<int, std::vector<merge_molecule_results_info_t> >
   merge_molecules(int imol, std::vector<mmdb::Manager *> mols);
#endif

   //! Convert a cis peptide to a trans or vice versa
   //!
   //! @param imol is the model molecule index
   //! @param atom_cid is the atom selection CID e.g "//A/15/OH" (atom OH residue 15 of chain A)
   //!
   //! @return 1 on a successful conversion.
   int cis_trans_convert(int imol, const std::string &atom_cid);

   //! Replace a residue
   //!
   //! This has a different meaning of "replace" to replace_fragment(). In this function
   //! the atoms are not merely moved/"slotted in to place", but the residue type is
   //! changed - new atoms are introduce and others are deleted (typically).
   //!
   //! Change the type of a residue (for example, "TYR" to "CYS")
   //! The algorithm will superpose the mainchain CA, C and N and try to set matching torsion
   //! to the angle that they were in the reference structure.
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the residue selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param new_residue_type is the 3-letter code of the new residue, e.g "CYS"
   //! @param imol_enc is the molecule index for the residue type/compound_id
   void replace_residue(int imol, const std::string &residue_cid, const std::string &new_residue_type, int imol_enc);

   //! Replace a fragment
   //!
   //! @param imol_base is the base model index
   //! @param imol_reference is the reference model index
   //! @param atom_selection is the selection CID e.g "//A/15-17" (residue 15, 16 and 17 of chain A)
   //!
   //! @return the success status
   int replace_fragment(int imol_base, int imol_reference, const std::string &atom_selection);

   //! Rigid-body fitting
   //!
   //! @param imol is the model molecule index
   //! @param multi_cids is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||/B/56-66"
   //! @param imol_map is the map molecule index
   int rigid_body_fit(int imol, const std::string &multi_cid, int imol_map);

   //! Rotate atoms around torsion
   //!
   //! the bond is presumed to be between atom-2 and atom-3. Atom-1 and atom-4 are
   //! used to define the absolute torsion angle.
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the residue selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param atom_name_1 e.g. " CA "
   //! @param atom_name_2 e.g. " CB "
   //! @param atom_name_3 e.g. " CG "
   //! @param atom_name_4 e.g. " CD1"
   //! @param torsion_angle e.g. 12.3 degrees
   //!
   //! @return status 1 if successful, 0 if not.
   int rotate_around_bond(int imol, const std::string &residue_cid,
                          const std::string &atom_name_1,
                          const std::string &atom_name_2,
                          const std::string &atom_name_3,
                          const std::string &atom_name_4,
                          double torsion_angle);

   //! Change the chain ID
   //!
   //! @param imol is the model molecule index
   //! @param from_chain_id e.g. "A"
   //! @param to_chain_id e.g. "C"
   //! @param use_resno_range use residue number range, typically True
   //! @param start_resno the starting residue number of the range
   //! @param end_resno the ending residue number of the range
   //!
   //! @return -1 on a conflict, 1 on good, 0 on did nothing,
   //! return also an information/error message
   std::pair<int, std::string> change_chain_id(int imol, const std::string &from_chain_id,
                                               const std::string &to_chain_id,
                                               bool use_resno_range,
                                               int start_resno, int end_resno);

   //! Split a residue into alt-confs
   //!
   //! do nothing if the residue already has alt-confs.
   //!
   //! @param imol the modified model
   //! @param residue_cid the modified residue, the residue selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param imol_diff_map is the difference map that is used to determine the residue split
   //!
   //! @return split success status
   int split_residue_using_map(int imol, const std::string &residue_cid, int imol_diff_map);

   //! Associate a sequence with a molecule
   //!
   //! @param imol is the model molecule index
   //! @param name_or_chain_id e.g. "A"
   //! @param sequence is the model sequence
   void associate_sequence(int imol, const std::string &name_or_chain_id, const std::string &sequence);

   //! Assign a sequence to a molecule
   //!
   //! Often one might copy out a fragment from a more complete
   //! molecule (and then copy it back after the sequence has been added). This runs
   //! `backrub_rotamer()` on the newly assigned residues
   //!
   //! @param imol is the model molecule index
   //! @param imol_map is the map molecule index
   void assign_sequence(int imol_model, int imol_map);

   //! Get the sequence information
   //!
   //! @param imol is the molecule index
   //! @return the sequence information
   std::vector<std::pair<std::string, std::string> > get_sequence_info(int imol) const;

   //! Get mutation information
   //!
   //! The reference sequece is that which has been provided using the
   //! `associate_sequence()` function
   //!
   //! @param imol is the model molecule index
   //! @return the mismatches/mutations as insertions, deletions or mutations
   coot::chain_mutation_info_container_t get_mutation_info(int imol) const;

   // -------------------------------- Coordinates Refinement ------------------------------
   //! \name Coordinates Refinement

   //! Refine the residues using cid
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param mode is the mode of real space refinement e.g. "SINGLE", "TRIPLE", "QUINTUPLE", "HEPTUPLE",
   //!              "SPHERE", "BIG_SPHERE", "CHAIN", "ALL"
   //! @param n_cycles is the number of refinement cycles
   //!
   //! @return a value of 1 if the refinement was performed and 0 if it was not.
   int refine_residues_using_atom_cid(int imol, const std::string &cid, const std::string &mode, int n_cycles);

   //! Refine the residues
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no is the residue number, e.g. 12
   //! @param ins_code is the insertion code, e.g. "A"
   //! @param alt_conf is the alternate conformation, e.g. "A" or "B"
   //! @param mode is the mode of real space refinement e.g. "SINGLE", "TRIPLE", "QUINTUPLE", "HEPTUPLE",
   //!               "SPHERE", "BIG_SPHERE", "CHAIN", "ALL"
   //! @param n_cycles is the number of refinement cycles
   //!
   //! @return a value of 1 if the refinement was performed and 0 if it was not.
   int refine_residues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                       const std::string &alt_conf, const std::string &mode, int n_cycles);

   //! Refine residue range
   //!
   //! @param imol is the model molecule index
   //! @param chain_id e.g. "A" for chain A
   //! @param res_no_start the starting residue number
   //! @param res_no_ref_end the ending residue number
   //! @param n_cycles is the number of refinement cycles
   //!
   //! @returns a value of 1 if the refinement was performed and 0 if it was not.
   int refine_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end, int n_cycles);

   //! Minimise/optimise the geometry of the specified residue(s)
   //!
   //! The use of "energy" should not be taken literally here
   //!
   //! @param imol is the model molecule index
   //! @param atom_selection_cid is the selection CID e.g. "//A/15" (residue 15 of chain A)
   //! @param n_cycles is the number of refinement cycles. If you pass n_cycles = 100 (or some such) then you can
   //!         get the mesh for the partially optimized ligand/residues
   //! @param do_rama_plot_restraints is the flag for the usage of Ramachandran plot restraints
   //! @param rama_plot_weight is the flag to set the Ramachandran plot restraints weight
   //! @param do_torsion_restraints is the flag for the usage of torsion restraints
   //! @param torsion_weight is the flag to set the torsion restraints weight
   //! @param refinement_is_quiet is used to reduce the amount of diagnostic text written to the output
   //!
   //! @return the success status 1 if the minimization was performed and 0 if it was not.
   std::pair<int, coot::instanced_mesh_t>
   minimize_energy(int imol, const std::string &atom_selection_cid,
                   int n_cycles,
                   bool do_rama_plot_restraints, float rama_plot_weight,
                   bool do_torsion_restraints, float torsion_weight, bool refinement_is_quiet);

   //! @param imol is the model molecule index
   //! @param atom_selection_cid is the selection CID e.g. "//A/15" (residue 15 of chain A)
   //! @param n_cycles is the number of refinement cycles. If you pass n_cycles = 100 (or some such) then you can
   //!         get the mesh for the partially optimized ligand/residues
   //! @param do_rama_plot_restraints is the flag for the usage of Ramachandran plot restraints
   //! @param rama_plot_weight is the flag to set the Ramachandran plot restraints weight
   //! @param do_torsion_restraints is the flag for the usage of torsion restraints
   //! @param torsion_weight is the flag to set the torsion restraints weight
   //! @param refinement_is_quiet is used to reduce the amount of diagnostic text written to the output
   //!
   //! @return the function value at termination
   float
   minimize(int imol, const std::string &atom_selection_cid,
            int n_cycles,
            bool do_rama_plot_restraints, float rama_plot_weight,
            bool do_torsion_restraints, float torsion_weight, bool refinement_is_quiet);

   //! Fix atoms during refinement
   //!
   //! Does nothing at the moment
   //!
   //! @param imol is the model molecule index
   //! @param atom_selection_cid is the selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A)
   void fix_atom_selection_during_refinement(int imol, const std::string &atom_selection_cid);

   //! Add or update restraint (if it has a pull restraint already)
   //!
   //! @param imol is the model molecule index
   //! @param atom_cid is the selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A)
   //! @param pos_x is the x coordinate of the target position of the specified atom
   //! @param pos_y is the y coordinate of the target position of the specified atom
   //! @param pos_z is the z coordinate of the target position of the specified atom
   void add_target_position_restraint(int imol, const std::string &atom_cid, float pos_x, float pos_y, float pos_z);

   //! Clear target_position restraint
   //!
   //! @param imol is the model molecule index
   //! @param atom_cid is the selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A)
   void clear_target_position_restraint(int imol, const std::string &atom_cid);

   //! Clear target_position restraint if it is (or they are) close to their target position
   //!
   //! @param imol is the model molecule index
   void turn_off_when_close_target_position_restraint(int imol);

   //! Control the logging
   //!
   //! @param level is the logging level, level is either "LOW" or "HIGH" or "DEBUGGING"
   void set_logging_level(const std::string &level);

   //! make the logging output go to a file
   //!
   //! @param file_name the looging file name
   void set_logging_file(const std::string &file_name);

   //! Turn on or off ramachandran restraints
   //!
   //! @param state is True to mean that it is enabled
   void set_use_rama_plot_restraints(bool state) { use_rama_plot_restraints = state; }

   //! Get the state of the rama plot restraints usage in refinement
   //!
   //! @return the state
   bool get_use_rama_plot_restraints() const { return use_rama_plot_restraints; }

   //! Set the Ramachandran plot restraints weight
   //!
   //! @param f is the weight to set, default 1.0
   void set_rama_plot_restraints_weight(float f) { rama_plot_restraints_weight = f; }

   //! Get the Ramachandran plot restraints weight
   //!
   //! @return the Ramachandran plot restraints weight
   float get_rama_plot_restraints_weight() const { return rama_plot_restraints_weight; }

   //! Turn on or off torsion restraints
   //!
   //! @param state is True to mean that it is enabled
   void set_use_torsion_restraints(bool state) { use_torsion_restraints = state; }

   //! Get the state of the rama plot restraints usage in refinement
   //!
   //! @return the state
   bool get_use_torsion_restraints() const { return use_torsion_restraints; }

   //! Set the torsiont restraints weight
   //!
   //! @param f is the weight to set, default value is 1.0
   void set_torsion_restraints_weight(float f) { torsion_restraints_weight = f; }

   //! Get the torsion restraints weight
   //!
   //! @return the torsion restraints weight
   float get_torsion_restraints_weight() const { return torsion_restraints_weight; }

   //! Initialise the refinement of (all of) molecule imol_frag
   //!
   //! @param imol_frag is the model molecule index of the fragment
   //! @param imol_ref is the model molecule index of the reference
   //! @param imol_map is the map molecule index
   void init_refinement_of_molecule_as_fragment_based_on_reference(int imol_frag, int imol_ref, int imol_map);

   //! Run some cycles of refinement and return a mesh
   //!
   //! That way we can see the molecule animate as it refines
   //!
   //! @param imol is the model molecule index
   //! @param n_cycles is the number of refinement cycles
   //!
   //! @return a pair: the first of which is the status of the refinement: GSL_CONTINUE, GSL_SUCCESS, GSL_ENOPROG (no progress).
   //! i.e. don't call this function again unless the status is GSL_CONTINUE (-2);
   //! The second is a coot::instanced_mesh_t
   std::pair<int, coot::instanced_mesh_t> refine(int imol, int n_cycles);

   //! Create a new position for the given atom and create a new bonds mesh based on that
   //!
   //! This is currently "heavyweight" as the bonds mesh is calculated from scratch (it is not (yet) merely a distortion
   //! of an internally-stored mesh).
   //!
   //! @param imol is the model molecule index
   //! @param atom_cid is the selection CID e.g "//A/15/OH" (atom OH of residue 15 of chain A)
   //! @param pos_x is the x coordinate of the target position of the specified atom
   //! @param pos_y is the y coordinate of the target position of the specified atom
   //! @param pos_z is the z coordinate of the target position of the specified atom
   //! @param n_cycles specifies the number of refinement cyles to run after the target position of the atom has been applied.
   //! If n_cycles is -1 then, no cycles are done and the mesh is bonds merely calculated.
   //!
   //! @return a `instanced_mesh_t`
   coot::instanced_mesh_t add_target_position_restraint_and_refine(int imol, const std::string &atom_cid,
                                                                   float pos_x, float pos_y, float pos_z,
                                                                   int n_cycles);
   //! Clear any and all drag-atom target position restraints
   //!
   //! @param imol is the model molecule index
   void clear_target_position_restraints(int imol);

   //! Call this after molecule refinement has finished (say when the molecule molecule is accepted into the
   //! original molecule)
   //!
   //! @param imol is the model molecule index
   void clear_refinement(int imol);

   //! For debugging the refinement - write out some diagnositics - some might be useful.
   //!
   //! API change 20240226 - this function now takes a boolean argument
   //!
   //! @param state is True to mean that it is enabled
   void set_refinement_is_verbose(bool state) { refinement_is_quiet = !state; }

   //! Set the refinement Geman-McClure alpha
   //!
   //! @param a is the Geman-McClure alpha, e.g. 0.01
   void set_refinement_geman_mcclure_alpha(float a) { geman_mcclure_alpha = a; }

   //! Get the refinement Geman-McClure alpha
   //!
   //! @return the Geman-McClure alpha
   float get_geman_mcclure_alpha() const { return geman_mcclure_alpha; }

   //! Generate GM self restraints for the whole molecule
   //!
   //! @param imol is the model molecule index
   //! @param local_dist_max is the maximum distance, e.g. 4.6
   int generate_self_restraints(int imol, float local_dist_max);

   //! Generate GM self restraints for the given chain
   //!
   //! @param imol is the model molecule index
   //! @param local_dist_max is the maximum distance, e.g. 4.6
   //! @param chain_id e.g. "A" for chain A
   void generate_chain_self_restraints(int imol, float local_dist_max,
                                       const std::string &chain_id);

   //! Generate GM self restraints for the given residues
   //!
   //! @param imol is the model molecule index
   //! @param local_dist_max is the maximum distance, e.g. 4.6
   //! @param residue_cids is a "||"-separated list of residues, e.g. "//A/12||//A/14||//B/56"
   void generate_local_self_restraints(int imol, float local_dist_max,
                                       const std::string &residue_cids);

   //! Generate parallel plane restraints (for RNA and DNA)
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid_1 is the selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param residue_cid_2 is the selection CID e.g "//A/17" (residue 17 of chain A)
   void add_parallel_plane_restraint(int imol,
                                     const std::string &residue_cid_1,
                                     const std::string &residue_cid_2);

   //! Get the mesh for extra restraints
   //!
   //! @param imol is the model molecule index
   //! @param mode is currently unused
   coot::instanced_mesh_t get_extra_restraints_mesh(int imol, int mode);

   //! Read extra restraints (e.g. from ProSMART)
   //!
   //! @param imol is the model molecule index
   int read_extra_restraints(int imol, const std::string &file_name);

   //! Clear the extra restraints
   //!
   //! @param imol is the model molecule index
   void clear_extra_restraints(int imol);

   //! External refinement using servalcat, using data that has already been associated.
   //!
   //! @param imol is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param output_prefix is the prefix of the output filename, e.g. "ref-1"
   //!
   //! @return the imol of the refined model.
   int servalcat_refine_xray(int imol, int imol_map, const std::string &output_prefix);

#if NB_VERSION_MAJOR
   //! Use servalcat keywords
   //!
   //! @param imol is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param output_prefix is the prefix of the output filename, e.g. "ref-1"
   //! @param key_value_pairs is a dictionary of key-value pairs for the servalcat keywords, e.g. resolution: 2.05
   //!
   //! @return the imol of the refined model.
   int servalcat_refine_xray_with_keywords(int imol, int imol_map, const std::string &output_prefix,
                                           const nanobind::dict &key_value_pairs);
#endif

   // -------------------------------- Coordinates validation ------------------------------
   //! \name Coordinates Validation

   //! Get the rotamer dodecs for the model
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a `simple_mesh_t`
   coot::simple_mesh_t get_rotamer_dodecs(int imol);

   //! Get the rotamer dodecs for the model instanced
   //!
   //! @param imol is the model molecule index
   //!
   //! @return an `instanced_mesh_t`
   coot::instanced_mesh_t get_rotamer_dodecs_instanced(int imol);

   //! Get the Ramachandran validation markup mesh
   //!
   //! 20221126-PE: the function was renamed from ramachandran_validation_markup_mesh().
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a `simple_mesh_t`
   coot::simple_mesh_t get_ramachandran_validation_markup_mesh(int imol) const;

   //! Get the data for Ramachandran validation, which importantly contains probability information
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a vector/list of `phi_psi_prob_t`
   std::vector<coot::phi_psi_prob_t> ramachandran_validation(int imol) const;

   //! Contact dots for ligand
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param smoothness_factor is 1, 2 or 3 (3 is the most smooth). Recently added (20230202)
   //!
   //! @return the instanced mesh for the specified ligand
   coot::instanced_mesh_t contact_dots_for_ligand(int imol, const std::string &cid, unsigned int smoothness_factor) const;

   //! Contact dots for the whole molecule/model
   //!
   //! @param imol is the model molecule index
   //! @param smoothness_factor is 1, 2 or 3 (3 is the most smooth). Recently added (20230202)
   //!
   //! @return the instanced mesh for the specified molecule.
   coot::instanced_mesh_t all_molecule_contact_dots(int imol, unsigned int smoothness_factor) const;

   //! Get a simple molecule
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param draw_hydrogen_atoms_flag is the flag for drawing H atoms
   //!
   //! @return a simple::molecule_t for the specified residue.
   // Note this function is not const because we pass a pointer to the protein_geometry geom.
   coot::simple::molecule_t get_simple_molecule(int imol, const std::string &residue_cid, bool draw_hydrogen_atoms_flag);

   //! @param imol is the model molecule index
   //! @param spec is the residue specifier, e.g. residue_spec_t("A", 10, "")
   //! @param max_dist specifies the maximum distance of the interaction, typically 3.8
   //!
   //! @return a vector of lines for non-bonded contacts and hydrogen bonds
   generic_3d_lines_bonds_box_t
   make_exportable_environment_bond_box(int imol, coot::residue_spec_t &spec, float max_dist);

   //! Get hydrogen bonds
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //! @param mcdonald_and_thornton_mode turns on the McDonald & Thornton algorithm - using explicit hydrogen atoms
   //!
   //! @return a vector of hydrogen bonds around the specified residue (typically a ligand)
   std::vector<moorhen::h_bond> get_h_bonds(int imol, const std::string &cid_str, bool mcdonald_and_thornton_mode) const;

   //! Get the mesh for ligand validation vs dictionary, coloured by badness
   //!
   //! Greater then 3 standard deviations is fully red.
   //! Less than 0.5 standard deviations is fully green
   //!
   //! @param imol is the model molecule index
   //! @param ligand_cid is the ligand selection CID e.g "//A/15" (ligand 15 of chain A)
   // Function is not const because it might change the protein_geometry geom.
   coot::simple_mesh_t get_mesh_for_ligand_validation_vs_dictionary(int imol, const std::string &ligand_cid);

   //! Ligand validation
   //!
   //! @param imol is the model molecule index
   //! @param ligand_cid is the ligand selection CID e.g "//A/15" (ligand 15 of chain A)
   //! @param include_non_bonded_contacts is the flag to include non bonded contacts
   //!
   //! @return a vector/list of interesting geometry - one for each chain involved
   std::vector<coot::geometry_distortion_info_container_t>
   get_ligand_validation_vs_dictionary(int imol, const std::string &ligand_cid, bool include_non_bonded_contacts);

   //! General fragment distortion analysis
   //!
   //! @param imol is the model molecule index
   //! @param selection_cid is the selection CID e.g "//A/15-23"
   //! @param include_non_bonded_contacts is the flag to include non bonded contacts
   //!
   //! @return a vector/list of interesting geometry - one for each chain involved
   std::vector<coot::geometry_distortion_info_container_t>
   get_validation_vs_dictionary_for_selection(int imol, const std::string &selection_cid, bool include_non_bonded_contacts);

   //! Get ligand distortion
   //!
   //! a more simple interface to the above
   //!
   //! @param imol is the model molecule index
   //! @param selection_cid is the selection CID e.g "//A/15-23"
   //! @param include_non_bonded_contacts is the flag to include non bonded contacts
   //!
   //! @return a pair: the first is the status (1 for OK, 0 for failed to determine the distortion)
   std::pair<int, double> get_ligand_distortion(int imol, const std::string &ligand_cid, bool include_non_bonded_contacts);

   //! Match ligand torsions
   //!
   //! @param imol_ligand is the ligand molecule index
   //! @param imol_ref is the reference model molecule index
   //! @param chain_id_ref is the reference chain, e.g. "A"
   //! @param resno_ref is the reference residue number, e.g. 12
   //!
   //! @return the success status
   bool match_ligand_torsions(int imol_ligand, int imol_ref, const std::string &chain_id_ref, int resno_ref);

   //! Match ligand positions
   //!
   //! i.e. do a least-squares superposition of the atoms that match in the graphs of the
   //! two specified ligands - typically one would use this function after matching ligand torsions.
   //!
   //! @param imol_ligand is the ligand molecule index
   //! @param imol_ref is the reference model molecule index
   //! @param chain_id_ref is the reference chain, e.g. "A"
   //! @param resno_ref is the reference residue number, e.g. 12
   //!
   //! @return the success status
   bool match_ligand_position(int imol_ligand, int imol_ref, const std::string &chain_id_ref, int resno_ref);

   //! Match ligand torsions and positions
   //!
   //! @param imol_ligand is the ligand molecule index
   //! @param imol_ref is the reference model molecule index
   //! @param chain_id_ref is the reference chain, e.g. "A"
   //! @param resno_ref is the reference residue number, e.g. 12
   //!
   //! @return the success status.
   bool match_ligand_torsions_and_position(int imol_ligand, int imol_ref, const std::string &chain_id_ref, int resno_ref);

   //! Match ligand torsions and positions using cid
   //!
   //! @param imol_ligand is the ligand molecule index
   //! @param imol_ref is the reference model molecule index
   //! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   bool match_ligand_torsions_and_position_using_cid(int imol_ligand, int imol_ref, const std::string &cid);

   // not const because it can dynamically add dictionaries
   //! @param imol is the model molecule index
   coot::atom_overlaps_dots_container_t get_overlap_dots(int imol);

   //! This function not const because it can dynamically add dictionaries
   //!
   //! @param imol is the model molecule index
   //! @param cid_ligand is the ligand selection CID e.g "//A/15" (ligand 15 of chain A)
   coot::atom_overlaps_dots_container_t get_overlap_dots_for_ligand(int imol, const std::string &cid_ligand);

   //! Get Atom Overlaps
   // not const because it can dynamically add dictionaries
   //! This function used to be called get_overlaps()
   //!
   //! @param imol is the model molecule index
   //! @return a vector of atom overlap objects
   std::vector<coot::plain_atom_overlap_t> get_atom_overlaps(int imol);

   //! Get the atom overlap score
   //!
   //! @param imol the model molecule index
   //! @return the overlap score - a negative number indicates failure
   float get_atom_overlap_score(int imol);

   //! Gat Atom Overlaps for a ligand or residue
   // not const because it can dynamically add dictionaries
   //! @param imol is the model molecule index
   //! @param cid_ligand is the ligand selection CID e.g "//A/15" (ligand 15 of chain A)
   //! @return a vector of atom overlap objects
   std::vector<coot::plain_atom_overlap_t> get_overlaps_for_ligand(int imol, const std::string &cid_ligand);

   //! Get the atom differences between two molecules
   //! typically after refinement
   //!
   //! @param imol1 is the first model molecule index
   //! @param imol2 is the second model molecule index
   //!
   //! @return a vector/list of `positioned_atom_spec_t`
   std::vector <positioned_atom_spec_t>
   get_atom_differences(int imol1, int imol2);

   //! get pucker info
   //!
   //! @param imol is the model molecule index
   //! @return a json string or an empty string on failure
   std::string get_pucker_analysis_info(int imol) const;


   // -------------------------------- Coordinates and map validation ----------------------
   //! \name Coordinates and Map Validation

   //! Density fit validation information.
   //!
   //! This function returns the sum of the densiy of the atoms in the residue
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map is the map molecule index
   //!
   //! @returns an object `validation_information_t`
   coot::validation_information_t density_fit_analysis(int imol_model, int imol_map) const;

   //! @return the sum of the density of the given atoms in the specified CID
   //!  return -1001 on failure to find the residue or any atoms in the residue or if imol_map is not a map
   double get_sum_density_for_atoms_in_residue(int imol, const std::string &cid,
                                               const std::vector<std::string> &atom_names,
                                               int imol_map);

   //! get the number of atoms in a given residue
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the selection CID e.g "//A/15" (residue 15 of chain A)
   //! @return the number of atoms in the residue, or -1 on failure
   int get_number_of_atoms_in_residue(int imol, const std::string &residue_cid) const;

   //! Get the density correlation validation information
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map is the map molecule index
   //!
   //! @returns an object `validation_information_t`
   coot::validation_information_t density_correlation_analysis(int imol_model, int imol_map) const;

   //! Get the rotamer validation information
   //!
   //! @param imol_model is the model molecule index
   //!
   //! @returns an object `validation_information_t`
   coot::validation_information_t rotamer_analysis(int imol_model) const;

   //! Get the ramachandran validation information (formatted for a graph, not 3D)
   //!
   //! @param imol_model is the model molecule index
   //!
   //! @returns an object `validation_information_t`
   coot::validation_information_t ramachandran_analysis(int imol_model) const;

   //! Get the ramachandran validation information (formatted for a graph, not 3D) for a given chain in a given molecule
   //!
   //! This function does not exist yet (20230127-PE)
   //!
   //! @param imol_model is the model molecule index
   //! @param chain_id e.g. "A"
   //!
   //! @returns an object `validation_information_t`
   coot::validation_information_t ramachandran_analysis_for_chain(int imol_model, const std::string &chain_id) const;

   //! Peptide omega validation information
   //!
   //! @param imol_model is the model molecule index
   //!
   //! @returns an object `validation_information_t`
   coot::validation_information_t peptide_omega_analysis(int imol_model) const;

   //! Get the median temperature factor for the model
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a negative number on failure
   float get_median_temperature_factor(int imol) const;

   //! Get the atom temperature factor
   //!
   //! @param imol is the model molecule index
   //! @param atom_cid is the selection cid for the atom
   //!
   //! @return a negative number on failure, otherwise the temperature factor
   float get_temperature_factor_of_atom(int imol, const std::string &atom_cid) const;

   //! Get interesting places
   //!
   //! This function does not work yet
   //!
   //! @return a vector/list of `validation_information_t`
   std::vector<coot::molecule_t::interesting_place_t> get_interesting_places(int imol, const std::string &mode) const;

   //! Get difference map peaks
   //!
   //! @param imol_map is the map molecule index
   //! @param imol_protein is the model molecule index
   //! @param n_rmsd number of sd, e.g. 4.8
   //!
   //! @return a vector/list of `validation_information_t`
   std::vector<coot::molecule_t::interesting_place_t> difference_map_peaks(int imol_map, int imol_protein, float n_rmsd) const;

   //! Get pepflips based on the difference map
   //!
   //! @param imol_coords is the model molecule index
   //! @param imol_difference_map is the difference map molecule index
   //! @param n_sigma number of sd, e.g. 4.8
   //!
   //! @return a vector/list of `validation_information_t`
   std::vector<coot::molecule_t::interesting_place_t> pepflips_using_difference_map(int imol_coords, int imol_difference_map, float n_sigma) const;

   //! Unmodelled blobs
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param rmsd_cut_off is the low map limit for cluster generation
   //!        1.4 is a reasonable value.
   //!
   //! @return a vector/list of `validation_information_t`
   std::vector<coot::molecule_t::interesting_place_t> unmodelled_blobs(int imol_model, int imol_map,
                                                                       float rmsd_cut_off) const;

   //! Check waters, using implicit logical OR
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param b_factor_lim typical value is 60.0
   //! @param outlier_sigma_level typical value is 0.8
   //! @param min_dist typical value is 2.3
   //! @param max_dist typical value is 3.5
   //! @param ignore_part_occ_contact_flag typical value is False
   //! @param ignore_zero_occ_flag typical value is False
   //!
   //! @return a vector/list of atom specifiers
   // Use the string_user_data of the spec for the button label
   std::vector <coot::atom_spec_t>
   find_water_baddies(int imol_model, int imol_map,
                      float b_factor_lim,
                      float outlier_sigma_level,
                      float min_dist, float max_dist,
                      bool ignore_part_occ_contact_flag,
                      bool ignore_zero_occ_flag);

   //! Get HOLE
   //!
   //! HOLE is a program for the analysis of the pore dimesions of ion channels. See Smart et al., 1996.
   //!
   //! @return a list of spheres on the surface of the pore
   coot::instanced_mesh_t get_HOLE(int imol,
                                   float start_pos_x, float start_pos_y, float start_pos_z,
                                   float end_pos_x, float end_pos_y, float end_pos_z) const;

   //! Calculate the MMRRCC for the residues in the chain
   //!
   //! Multi Masked Residue Range Corellation Coefficient
   //!
   //! @param imol is the model molecule index
   //! @param chain_id is the model chain_id
   //! @param n_residue_per_residue_range is the number of residues in the residue range. 11
   //!        is a reasonable number for a smooth plot
   //! @param imol_map is the map molecule index
#ifdef SWIG
#else
   std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
             std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
   get_mmrrcc(int imol, const std::string &chain_id, unsigned int n_residue_per_residue_range, int imol_map) const;

   //! This is a wrapper for get_mmrrcc(), using 11 for the `n_residue_per_residue_range`.
   //!
   //! Multi Masked Residue Range Corellation Coefficient
   //!
   //! @param imol is the model molecule index
   //! @param chain_id is the model chain_id
   //! @param imol_map is the map molecule index
   std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
             std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
   mmrrcc(int imol, const std::string &chain_id, int imol_map) const;
#endif

   // calculate the MMRRCC for the residues in the chain
   // Multi Masked Residue Range Corellation Coefficient
#ifdef SWIG
#else
   std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
             std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
   mmrrcc_internal(const atom_selection_container_t &asc,
                   const std::string &chain_id,
                   unsigned int n_residue_per_residue_range,
                   const clipper::Xmap<float> &xmap) const;
#endif

   //! Fourier Shell Correlation (FSC) between maps
   //!
   //! @param imol_map_1 is the first map molecule index
   //! @param imol_map_2 is the second map molecule index
   //!
   //! @return a vector/list or pairs of graph points (resolution, correlation). The resolution is in
   //! inverse Angstroms squared.
   //! An empty list is returned on failure
   std::vector<std::pair<double, double> > fourier_shell_correlation(int imol_map_1, int imol_map_2) const;

   //! Make a FSC-scaled map
   //!
   //! @param imol_ref is the reference map molecule index
   //! @param imol_map_for_scaling is the second map molecule index
   //!
   //! @return the molecule index of the new map
   int make_power_scaled_map(int imol_ref, int imol_map_for_scaling);

   //! Get the Q Score (Pintilie et al.)
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map is the map molecule index
   //!
   //! @return a `validation_information_t` object
   coot::validation_information_t get_q_score(int imol_model, int imol_map) const;

   //! Get the Pintile et al. Q Score for a particular residue (typically a ligand)
   //!
   //! @param cid If the `cid` matches more than one residue the score will be returned for all of the
   //! residues covered in the `cid`. Typically, of course the `cid` will be something like
   //! "//A/301".
   //!
   //! @return a `validation_information_t` object
   coot::validation_information_t get_q_score_for_cid(int imol_model, const std::string &cid, int imol_map) const;

   //! get mean and variance of map at non-waters
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map is the map molecule index
   //!
   //! @return the mean and variance or a negative number on failure
   std::pair<float,float> get_mean_and_variance_of_density_for_non_water_atoms(int imol_coords, int imol_map) const;

   //! Get spherical variance - typically for water atoms
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map is the map molecule index
   //!
   //! @return the variance or a negative number on failure
   float get_spherical_variance(int imol_map, int imol_model, const std::string &atom_cid, float mean_density_other_atoms) const;

   // -------------------------------- Rail Points ------------------------------------------
   //! \name Rail Points!

   //! Calling this adds to the rail_points history. Make this pairs when we add model scoring
   //!
   //! @returns the new rail points (since last modification)
   int calculate_new_rail_points();

   //! The total rail points
   //!
   //! @returns the sum of all rail points accumulated since the maps were connected.
   int rail_points_total() const;

   // -------------------------------- Updating Maps ---------------------------------------
   //! \name Updating Maps

   //! Associate a data mtz file with a molecule
   //!
   //! This function is called before calling "connect_updating_maps()"
   //!
   //! @param imol is the map molecule index
   //! @param data_mtz_file_name is the name of the mtz file
   //! @param f_col is the F column, e.g. "FOBS"
   //! @param sigf_col e.g. "SIGFOBS"
   //! @param free_r_col e.g. "RFREE"
   void associate_data_mtz_file_with_map(int imol, const std::string &data_mtz_file_name,
                                         const std::string &f_col, const std::string &sigf_col,
                                         const std::string &free_r_col);

   //! Connect updating maps
   //!
   //! Reset the rail_points (calls "reset_the_rail_points()"), updates the maps (using internal/clipper SFC).
   //! Update your contour lines meshes after calling this function.
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_with_data_info_attached is the map index with the data have been attached by the previous function (associate_data_mtz_file_with_map)
   //! @param imol_map_2fofc is the map molecule index of the 2FO-FC map
   //! @param imol_map_fofc is the map molecule index of the FO-FC map
   //!
   //! @return 1 if the connection was successful
   int connect_updating_maps(int imol_model, int imol_with_data_info_attached, int imol_map_2fofc, int imol_map_fofc);

   //! Calculate SF and re-generate maps
   //!
   //! This is a low-level function - generally one would use the updating maps method rather than this
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_map_with_data_attached is the map index with the data have been attached by the previous function (associate_data_mtz_file_with_map)
   //! @param imol_updating_difference_map is the index of the difference map that you want to update when the model updates
   void sfcalc_genmap(int imol_model,
                      int imol_map_with_data_attached,
                      int imol_updating_difference_map);

   //! Calculate SF and re-generate maps using bulk solvent
   //!
   //! This is a low-level function. Call this function after connecting maps for updating maps to set the initial R-factor
   //! and store the initial map flatness.
   //!
   //! @param imol_model is the model molecule index
   //! @param imol_2fofc_map is the map molecule index of the 2FO-FC map
   //! @param imol_updating_difference_map is the index of the difference map that you want to update when the model updates
   //! @param imol_map_with_data_attached is the map index with the data have been attached by the previous function
   //!        (associate_data_mtz_file_with_map)
   //!
   //! @return a class of interesting statistics.
   //! On failure to calculate SFs and generate the maps the returned r_factor
   //! in the returned stats will be set to -1.
   coot::util::sfcalc_genmap_stats_t
   sfcalc_genmaps_using_bulk_solvent(int imol_model,
                                     int imol_2fofc_map,
                                     int imol_updating_difference_map,
                                     int imol_map_with_data_attached);

   //! Shift-field B-factor refinement
   //!
   //! This function presumes that the Fobs, sigFobs
   //! and RFree data have been filled in the `imol_map_with_data_attached` molecule
   //!
   //! @param imol is the model molecule index
   //! @param imol_map_with_data_attached is the map index with the data have been attached by the previous function
   //!        (associate_data_mtz_file_with_map)
   //!
   //! @return success status
   bool shift_field_b_factor_refinement(int imol, int imol_with_data_attached);

   //! Get density at position
   //! @param imol_map is the map molecule index
   //! @param x is the x coordinate of the target position
   //! @param y is the y coordinate of the target position
   //! @param z is the z coordinate of the target position
   //!
   //! @return the density value
   float get_density_at_position(int imol_map, float x, float y, float z) const;

   //! @param imol_diff_map is the map molecule index of the difference map
   //! @param screen_centre_x is the position x of the center of the screen
   //! @param screen_centre_y is the position y of the center of the screen
   //! @param screen_centre_z is the position z of the center of the screen
   //!
   //! @return a vector of the position where the difference map has been flattened.
   //! The associated float value is the amount that the map has been flattened.
   // This is a light-weight fetch, the values have already been computed, here
   // were are merely copying them.
   std::vector<std::pair<clipper::Coord_orth, float> > get_diff_diff_map_peaks(int imol_diff_map,
                                                                               float screen_centre_x,
                                                                               float screen_centre_y,
                                                                               float screen_centre_z) const;


   //! Get the stored data set file name
   //!
   //! @param imol is the model molecule index
   std::string get_data_set_file_name(int imol) const;

   // -------------------------------- Go To Blob ---------------------------------------
   //! \name Go to Blob

   //! Given a point on the front clipping plane (x1, y1, z1) and a point on the back clipping plane (x2, y2, z2)
   //! this function searches imol_refinement_map (if set) to find a the centre of a blob above the contour level.
   //! Blobs at the "front" are selected in preference to blobs at the back.
   //! If no blob is found, then the first of the pair is false.
   //! If it is found, then the second is the centre of the blob.
   //!
   //! In future, this function should/will be provided with a list of displayed maps and their
   //! contour levels - but for now, it uses (only) imol_refinement_map.
   //!
   //! Use a string to pass the map information (map index and contour level), something like "0 0.45:1 1.2:2 0.1"
   //! @param x1 is the x point of the front clipping plane
   //! @param y1 is the y point of the front clipping plane
   //! @param z1 is the z point of the front clipping plane
   //! @param x2 is the x point of the back clipping plane
   //! @param y2 is the y point of the back clipping plane
   //! @param z2 is the z point of the back clipping plane
   std::pair<bool, clipper::Coord_orth> go_to_blob(float x1, float y1, float z1, float x2, float y2, float z2,
                                                   float contour_level);


   // -------------------------------- Ligand Functions ---------------------------------------
   //! \name Ligand Functions

   //! Fit the ligand at specified position.
   //!
   // I am not yet clear what extra cut-offs and flags need to be added here.
   //! You can expect this to take about 20 seconds.
   //! For trivial (i.e non-flexible) ligands you should instead use the jiggle-fit algorithm, which
   //! takes a fraction of a second. (That is the algorithm used for "Add Other Solvent Molecules" in Coot.)
   //!
   //! @param imol_protein is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param imol_ligand is the ligand molecule index
   //! @param x is the x position of the blob
   //! @param y is the y position of the blob
   //! @param z is the z position of the blob
   //! @param n_rmsd number of sd, e.g. 4.8
   //! @param use_conformers is True for flexible ligands
   //! @param n_conformers set the number of conformers
   //!
   //! @return a vector/list of indices of molecules for the best fitting ligands to this blob.
   std::vector<int> fit_ligand_right_here(int imol_protein, int imol_map, int imol_ligand, float x, float y, float z,
                                          float n_rmsd, bool use_conformers, unsigned int n_conformers);

   //! Ligand Fitting
   //!
   //! @return a vector of indices of molecules for the best fitting ligands each of the "possible ligand" blobs.

   class fit_ligand_info_t {
   public:
      int imol; // the imol of the fitted ligand
      int cluster_idx;  // the index of the cluster
      int ligand_idx;  // the ligand idx for a given cluster
      float fitting_score;
      float cluster_volume;
      fit_ligand_info_t(int i, int c, int l) : imol(i), cluster_idx(c), ligand_idx(l) { fitting_score = -1.0; cluster_volume = -1.0;}
      fit_ligand_info_t() { imol = -1; cluster_idx = -1; ligand_idx = -1; fitting_score = 1.0; cluster_volume = -1.0; }
      //! @return the fitting score
      float get_fitting_score() const { return fitting_score; }
      float get_cluster_volume() const { return cluster_volume; }
   };

   //! Fit ligand
   //!
   //! @param imol_protein is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param imol_ligand is the ligand molecule index
   //! @param n_rmsd the number of sd used as a cut-off for the map level when finding clusters, e.g. 1.2
   //! @param use_conformers is True for flexible ligands
   //! @param n_conformers set the number of conformers
   //!
   //! @return a vector/list of interesting information about the fitted ligands
   std::vector<fit_ligand_info_t> fit_ligand(int imol_protein, int imol_map, int imol_ligand,
                                             float n_rmsd, bool use_conformers, unsigned int n_conformers);

   //! Fit multiple ligands (place-holder)
   //!
   //! @param imol_protein is the model molecule index
   //! @param imol_map is the map molecule index
   //! @param multi_ligand_molecule_number_list is a colon-separated list of molecules, e.g. "2:3:4"
   //! @param n_rmsd is number of sd, e.g. 4.8
   //! @param use_conformers is True for flexible ligands
   //! @param n_conformers set the number of conformers
   //!
   //! @return an empty vector (at the moment)
   std::vector<fit_ligand_info_t> fit_ligand_multi_ligand(int imol_protein, int imol_map, const std::string &multi_ligand_molecule_number_list,
                                                          float n_rmsd, bool use_conformers, unsigned int n_conformers);

   //! Jiggle-Fit Ligand
   //!
   //! @param imol is the model molecule index
   //! @param res_spec is the residue specifier, e.g. residue_spec_t("A", 10, "")
   //! @param n_trials is the number of trials, if n_trials is 0, then a sensible default value will be used.
   //! @param translation_scale_factor if is negative then a sensible default value will be used.
   //!
   //! @return a value less than -99.9 on failure to fit.
   float fit_to_map_by_random_jiggle(int imol, const coot::residue_spec_t &res_spec, int n_trials, float translation_scale_factor);

   //! Jiggle-Fit Ligand using cid
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID, e.g "//A/15" (ligand 15 of chain A)
   //! @param n_trials is the number of trials, if n_trials is 0, then a sensible default value will be used.
   //! @param translation_scale_factor if is negative then a sensible default value will be used.
   //!
   //! @return a value less than -99.9 on failure to fit.
   float fit_to_map_by_random_jiggle_using_cid(int imol, const std::string &cid, int n_trials, float translation_scale_factor);

   //! Jiggle-Fit an atom selection, typically a whole molecule or a chain
   //!
   //! @param imol is the model molecule index
   //! @param cid is the selection CID, e.g. "//A" (chain A)
   //! @param b_factor e.g. 100.0
   //! @param n_trials is the number of trials, if n_trials is 0, then a sensible default value will be used.
   //! @param translation_scale_factor if is negative then a sensible default value will be used.
   //!
   //! @return a value less than -99.9 on failure to fit.
   float fit_to_map_by_random_jiggle_with_blur_using_cid(int imol, int imol_map, const std::string &cid, float b_factor,
                                                         int n_trials, float translation_scale_factor);

   //! This function is for adding compounds/molecules like buffer agents and precipitants or anions and cations.
   //! e.g. those ligands that can be positioned without need for internal torsion angle manipulation.
   //!
   //! @param imol is the model molecule index
   //! @param tlc is the 3-letter-code/compound-id
   //! @param imol_dict is the molecule to which the ligand is attached (if any). Typically this will be IMOL_ENC_ANY (-999999).
   //! @param imol_map is the map molecule index
   //! @param x is the x position
   //! @param y is the y position
   //! @param z is the z position
   //!
   //! @return the success status, 1 for good, 0 for not good.

   int add_compound(int imol, const std::string &tlc, int imol_dict, int imol_map, float x, float y, float z);
   // This is a ligand function, not really a ligand-fitting function.
   //!
   //! Get svg for residue type
   //!
   //! It won't work unless the dictionary for that ligand has been imported.
   //! The native output renderings are not very good at the moment.
   //! (The RDKit renderings are pretty good).
   //!
   //! @param imol is the model molecule index, except for unusual cases, it will be IMOL_ENC_ANY (-999999)
   //! @param comp_id is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine
   //! @param use_rdkit_svg is the flag for using the rdkit svg renderer
   //! @param background_type is one of:
   //!  - "light-bonds/transparent-bg"
   //!  - "light-bonds/opaque-bg"
   //!  - "dark-bonds/transparent-bg"
   //!  - "dark-bonds/opaque-bg"
   //!
   //! If you want to load them into another image, you'd typicaly want "dark-bonds/transparent-bg"
   //! If you want to see ligands, e.g. in a grid or list, you'd typically want "dark-bonds/opaque-bg"
   //! which will give you a white rectangle behind the ligand figure.
   //!
   //! This function is not const because it caches the svgs.
   //!
   //! @return the string for the SVG representation.
   std::string get_svg_for_residue_type(int imol, const std::string &comp_id,
                                        bool use_rdkit_svg,
                                        const std::string &background_type);

   //! Get SVG for 2d ligand environment view (FLEV)
   //!
   //! The caller should make sure that the dictionary for the ligand has been loaded - this
   //! function won't do that. It will add hydrogen atoms if needed.
   //!
   //! From time to time (depending on the ligand) this function will fail to produce a
   //! result.
   //!
   //! Not const because get_monomer_restraints_at_least_minimal() is called. Hmm.
   //!
   //! @param imol is the model molecule index
   //! @param residue_cid is the cid for the residue
   //! @param add_key should a key be added to the figure?
   //! @return an svg string of the representation. On failure, return an empty string.
   std::string get_svg_for_2d_ligand_environment_view(int imol, const std::string &residue_cid, bool add_key);

   //! Get non-standard residues in a model
   //!
   //! @param imol is the model molecule index
   //!
   //! @return a vector/list of residue specifiers - the residue name is encoded
   //! in the `string_user_data` data item of the residue specifier
   std::vector<coot::residue_spec_t> get_non_standard_residues_in_molecule(int imol) const;

   //! Try to read the dictionaries for any residue type in imol that as yet does not have
   //! a dictionary
   //!
   //! @param imol is the model molecule index
   //! @return true if there were no dictionary for new types that couldn't be read.
   bool try_read_dictionaries_for_new_residue_types(int imol);

   //! Get the conformers that can be generated by variation around rotatable bonds as described in the dictionary.
   //!
   //! Torsions that are marked as "const" are excluded from the variation, as are pyranose ring torsions
   //! and torsions that rotate hydrogen atoms
   //!
   //! @param comp_id is the 3-letter code for the residue/ligand, e.g. "ALA" for alanine
   //! @param imol_enc is the molecule index for the residue type/compound_id
   //! @param remove_internal_clash_conformers is the flag for removing internal clash
   //!
   //! @return a vector/list of indices of the new molecules
   std::vector<int> get_dictionary_conformers(const std::string &comp_id, int imol_enc, bool remove_internal_clash_conformers);

   //! @param imol is the map molecule index
   //! @param section_id e.g. 2
   //! @param axis e.g. 0 for X-axis, 1 for Y-axis, 2 for Z-axis
   //!
   //!
   //! The new arguments, `data_value_for_bottom`, `data_value_for_top` should be pre-calculated (don't
   //! calculate them for every call to this function).
   //!
   //! @return a texture_as_floats_t object for the given section.
   //! On failure, the image_data vector is empty.
   texture_as_floats_t get_map_section_texture(int imol, int section_id, int axis,
                                               float data_value_for_bottom, float data_value_for_top) const;

   //! @param imol_map is the map molecule index
   //! @param axis_id is 0 for X-axis, 1 for Y-axis, 2 for Z-axis
   //!
   //! @return the number of sections in the map along the given axis, -1 on failure.
   int get_number_of_map_sections(int imol_map, int axis_id) const;

   // -------------------------------- Others -------------------------------------
   //! \name Other Features

   //! @param file_name is the glTF file name
   //!
   //! @return a `simple_mesh_t` from the given file.
   coot::simple_mesh_t make_mesh_from_gltf_file(const std::string &file_name);

   //! Get octahemisphere
   //!
   //! @param n_divisions is a number divisible by 2, at least 4 (typically 16)
   //!
   //! @return a unit-vector end-cap octohemisphere mesh
   coot::simple_mesh_t get_octahemisphere(unsigned int n_divisions) const;

   unsigned int get_max_number_of_simple_mesh_vertices() const;
   void set_max_number_of_simple_mesh_vertices(unsigned int n);

   //! Predicted alignment error (AlphaFold)
   //! @return a string of a png
   std::string pae_png(const std::string &pae_file_name) const;

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

   //! Testing function
   //!
   //! start a long-term job.
   //!
   //! @param n_seconds is the number of seconds, if is 0, then run forever (or until interrupted)
   void testing_start_long_term_job(unsigned int n_seconds);

   //! Testing function
   //!
   //! stop the long-term job runnning
   void testing_stop_long_term_job();

   //! Testing function
   //!
   //! get the stats for the long-term job
   ltj_stats_t testing_interrogate_long_term_job() { return long_term_job_stats; }

   //! Testing function
   //!
   //! get the time for contouring in milliseconds
   double get_contouring_time() const { return contouring_time; }

   //! Testing function
   //!
   //! set the maximum number of threads for both the thread pool and the vector of threads
   //!
   //! @param n_threads is the number of threads
   void set_max_number_of_threads(unsigned int n_threads);

   //! Testing function
   //!
   //! Deprecated name for the "set_max_number_of_threads()" function
   void set_max_number_of_threads_in_thread_pool(unsigned int n_threads);

   //! Testing function
   //!
   //! get the time to run a test function in milliseconds
   //!
   //! @param n_threads is the number of threads
   double test_the_threading(int n_threads);

   //! Testing function
   //!
   //! @param n_threads_per_batch is the number of threads per batch
   //! @param n_batches is the number batches
   //!
   //! @return the time per batch in microseconds
   double test_launching_threads(unsigned int n_threads_per_batch, unsigned int n_batches) const;

   //! Testing function
   //!
   //! @param n_threads is the number of threads
   //!
   //! @return time in microseconds
   double test_thread_pool_threads(unsigned int n_threads);

   //! Testing function
   //!
   //! a test for mmdb/gemmi/mmcif functionality
   //!
   //! @param last_test_only is True to mean that only that last test should be run.
   //! The default is False.
   //! This is useful to set to True while a test is being developed.
   //!
   //! @return the success status: 1 means that all the tests passed.
   int mmcif_tests(bool last_test_only);

   // I want this function in the C++ documentation, but not the Python API documentation.
   // Hmm.
#ifdef DOXYGEN_SHOULD_PARSE_THIS
#else
   //! get access to protein geometry
   coot::protein_geometry & get_geometry() {
      return geom;
   }
#endif

   // -------------------------------- Blender Interface ---------------------------------------

   //! \name Functions for Blender Interface

   //! Function for Blender interface
   void make_mesh_for_map_contours_for_blender(int imol, float x, float y, float z, float level, float radius);
   //! Function for Blender interface
   void make_mesh_for_bonds_for_blender(int imol, const std::string &mode, bool against_a_dark_background,
                                      float bond_width, float atom_radius_to_bond_width_ratio,
                                      int smoothness_factor);
   //! Function for Blender interface
   //!
   //! Make an (internal) mesh
   //!
   //! This function doesn't return a value, instead it stores a `blender_mesh_t` blender_mesh
   //! in this model. One then (shortly later) uses get_triangles_for_blender(imol) (etc)
   //! to import this mesh into blender.
   //!
   //! @modifies internal state to fill the internal `blender_mesh` object
   //!
   //!
   void make_mesh_for_molecular_representation_for_blender(int imol,
                                                           const std::string &cid,
                                                           const std::string &colour_scheme,
                                                           const std::string &style,
                                                           int secondary_structure_usage_flag);
   //! Function for Blender interface
   void make_mesh_for_gaussian_surface_for_blender(int imol, float sigma, float contour_level, float box_radius, float grid_scale, float b_factor);
   //! blender

  //! Function for Blender interface
   void make_mesh_for_goodsell_style_for_blender(int imol, float colour_wheel_rotation_step,
                                                 float saturation, float goodselliness);

   //! Function for Blender interface
   std::vector<float> get_colour_table_for_blender(int imol);
   //! Function for Blender interface
   std::vector<float> get_vertices_for_blender(int imol);
   //! Function for Blender interface
   std::vector<int>   get_triangles_for_blender(int imol);

   // -------------------------------- Other ---------------------------------------

   void test_function(const std::string &s);

#if NB_VERSION_MAJOR
   // skip this (old) block for nanobinds
#else
#ifdef DOXYGEN_SHOULD_PARSE_THIS

   //! \name Old Python functions

   //! old mesh mode: do not use with nanobind
   enum mesh_mode_t { UNKNOWN, SINGLE_COLOUR, MULTI_COLOUR };
   //! old function: do not use with nanobind
   PyObject *simple_mesh_to_pythonic_mesh(const coot::simple_mesh_t &mesh, int mesh_mode);
   //! old function: do not use with nanobind
   PyObject *get_pythonic_bonds_mesh(int imol, const std::string &mode, bool against_a_dark_background,
                                     float bond_width, float atom_radius_to_bond_width_ratio,
                                     int smoothness_factor);
   //! old function: do not use with nanobind
   PyObject *get_pythonic_map_mesh(int imol, float x, float y, float z, float radius, float contour_level);
   //! old function: do not use with nanobind
   PyObject *get_pythonic_molecular_representation_mesh(int imol, const std::string &atom_selection,
                                                        const std::string &colour_sheme,
                                                        const std::string &style,
                                                        int secondary_structure_usage_flag);
   //! old function: do not use with nanobind get Gaussion surface mesh
   PyObject *get_pythonic_gaussian_surface_mesh(int imol, float sigma, float contour_level,
                                                float box_radius, float grid_scale, float fft_b_factor);
   //! old function: do not use with nanobind: get a pythonic mesh of the molecule (bonds)
   //!
   //! @return a pair - the first of which (index 0) is the list of atoms, the second (index 1) is the list of bonds.
   //! An atom is a list:
   //!
   //! 0: atom-name
   //!
   //! 1: atom-element
   //!
   //! 2: position (a list of 3 floats)
   //!
   //! 3: formal charge (an integer)
   //!
   //! 4: aromaticity flag (boolean)
   //!
   //! make a "proper" simple  molecule python class one day.
   PyObject *get_pythonic_simple_molecule(int imol, const std::string &cid, bool include_hydrogen_atoms_flag);

#endif
#endif

};

#endif // MOLECULES_CONTAINER_HH
