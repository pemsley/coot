#ifndef COOT_MOLECULE_HH
#define COOT_MOLECULE_HH

#include <utility>
#include <atomic>

#include <clipper/core/xmap.h>
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/coot-rama.hh"
#include "coot-utils/sfcalc-genmap.hh"
#include "coot-utils/atom-tree.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"
#include "ideal/simple-restraint.hh"
#include "ideal/extra-restraints.hh"
#include "coot-utils/simple-mesh.hh"
#include "ghost-molecule-display.hh"

#include "density-contour/CIsoSurface.h"
#include "gensurf.hh"
#include "coot-utils/coot-shelx.hh"

#include "coot-colour.hh" // put this in utils

#include "coords/mmdb-extras.h"
#include "merge-molecule-results-info-t.hh"
#include "phi-psi-prob.hh"

#include "coot-utils/atom-overlaps.hh"

#include "instancing.hh"

namespace coot {

   // give this a type
   enum { UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2,
          COLOUR_BY_CHAIN_BONDS=3,
          CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
          BONDS_NO_HYDROGENS=15,
          CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
          CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
          CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS=17,
          COLOUR_BY_MOLECULE_BONDS=8,
          COLOUR_BY_RAINBOW_BONDS=9,
          COLOUR_BY_B_FACTOR_BONDS=10,
          COLOUR_BY_OCCUPANCY_BONDS=11,
          COLOUR_BY_USER_DEFINED_COLOURS____BONDS=12,
          COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS=13 };

   class molecule_t {

      class molecule_save_info_t {
      public:
         // ints because we subtract from them.
         std::pair<time_t, int> last_saved;
         int modification_index;
         int max_modification_index;
         molecule_save_info_t() : last_saved(std::make_pair(0,0)), modification_index(0),
                                  max_modification_index(0) {}
         void new_modification(const std::string &mod_string) {
            modification_index++;
            std::cout << "new_modification! moved on to " << modification_index << " by " << mod_string << std::endl;
            if (modification_index > max_modification_index)
               max_modification_index = modification_index;
         }
         void made_a_save() {
            // this is called when the server says that it has saved the file
            last_saved.first = time(nullptr);
            last_saved.second = modification_index;
         }
         bool have_unsaved_changes() const {
            return modification_index > last_saved.second;
         }
         std::string index_string() const {
            return std::to_string(modification_index);
         }
         void set_modification_index(int idx) {
            // undo() and redo() (acutally restore_from_backup()) use this.
            // the index is the idx number given to restore_from_backup() by
            // the functions below
            modification_index = idx;
         }
         int get_previous_modification_index() const {
            return modification_index - 1;
         }
         int get_next_modification_index() const {
            return modification_index + 1;
         }
      };

      molecule_save_info_t save_info;

      int imol_no; // this molecule's index in the container vector
      int ligand_flip_number;
      std::string name;
      bool is_from_shelx_ins_flag;
      ShelxIns shelxins;

      // private
      void makebonds(coot::protein_geometry *geom, coot::rotamer_probability_tables *rotamer_tables_p, std::set<int> &no_bonds_to_these_atoms,
                     bool draw_hydrogen_atoms_flag, bool draw_missing_loops_flag);

#if defined __has_builtin
#if __has_builtin (__builtin_FUNCTION)
      void make_bonds_type_checked(coot::protein_geometry *geom, coot::rotamer_probability_tables *rot_prob_tables_p, bool draw_hydrogen_atoms_flag, bool draw_missing_loops_flag, const char *s = __builtin_FUNCTION());
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, bool draw_hydrogen_atoms_flag, bool draw_missing_loops_flag, const char *s = __builtin_FUNCTION());
#else
      void make_bonds_type_checked(coot::protein_geometry *geom, const char *s = 0);
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, bool draw_hydrogen_atoms_flag, bool draw_missing_loops_flag, const char *s =0);
#endif
#else // repeat above
      void make_bonds_type_checked(coot::protein_geometry *geom, const char *s = 0);
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, bool draw_hydrogen_atoms_flag, bool draw_missing_loops_flag, const char *s =0);
#endif

      int bonds_box_type; // public accessable via get_bonds_box_type(); // wass Bonds_box_type()
      graphical_bonds_container bonds_box;
      int get_bonds_box_type() const { return bonds_box_type; }

      // this is the bond dictionary also mode.
      // 20221011-PE force_rebonding arg is not currently used.
      void make_colour_by_chain_bonds(coot::protein_geometry *geom,
                                      const std::set<int> &no_bonds_to_these_atoms,
                                      bool change_c_only_flag,
                                      bool goodsell_mode,
                                      bool draw_hydrogen_atoms_flag,
                                      bool draw_missing_loops_flag,
                                      bool do_rota_markup=false,
                                      coot::rotamer_probability_tables *rotamer_tables_p = nullptr,
                                      bool force_rebonding=true);
      void make_ca_bonds();
      // just a copy of the version in src
      float bonds_colour_map_rotation;
      std::vector<glm::vec4> make_colour_table(bool against_a_dark_background) const;
      glm::vec4 get_bond_colour_by_colour_wheel_position(int icol, int bonds_box_type) const;
      coot::colour_t get_bond_colour_by_mol_no(int colour_index, bool against_a_dark_background) const;
      coot::colour_t get_bond_colour_basic(int colour_index, bool against_a_dark_background) const;
      bool use_bespoke_grey_colour_for_carbon_atoms;
      coot::colour_t bespoke_carbon_atoms_colour;

      void update_map_triangles(float radius, coot::Cartesian centre, float contour_level);

      bool is_EM_map() const;
      short int is_em_map_cached_flag; // -1 mean unset (so set it, 0 means no, 1 means yes)
      short int is_em_map_cached_state(); // set is_em_map_cached_flag if not set
      coot::ghost_molecule_display_t map_ghost_info;

      bool xmap_is_diff_map;
      bool has_xmap() const { return is_valid_map_molecule(); }

      // save the data used for the fourier, so that we can use it to
      // sharpen the map:
      // uncommenting the following line causes a crash in the multi-molecule
      // (expand molecule space) test.
      bool original_fphis_filled;
      bool original_fobs_sigfobs_filled;
      bool original_fobs_sigfobs_fill_tried_and_failed;
      clipper::HKL_data< clipper::datatypes::F_phi<float> >  *original_fphis_p;
      clipper::HKL_data< clipper::datatypes::F_sigF<float> > *original_fobs_sigfobs_p;
      clipper::HKL_data< clipper::data32::Flag> *original_r_free_flags_p;

      void clear_draw_vecs();
      void clear_diff_map_draw_vecs();
      std::vector<coot::density_contour_triangles_container_t> draw_vector_sets;
      std::vector<coot::density_contour_triangles_container_t> draw_diff_map_vector_sets;
      std::vector<std::pair<int, TRIANGLE> > map_triangle_centres; // with associated mid-points and indices

      // This function no longer does a backup or updates the save_info!
      // The calling function should do that.
      void replace_coords(const atom_selection_container_t &asc,
                          bool change_altconf_occs_flag,
                          bool replace_coords_with_zero_occ_flag);
      // helper function for above function
      bool movable_atom(mmdb::Atom *mol_atom, bool replace_coords_with_zero_occ_flag) const;
      bool moving_atom_matches(mmdb::Atom *at, int this_mol_index_maybe) const;
      void adjust_occupancy_other_residue_atoms(mmdb::Atom *at,
                                                mmdb::Residue *residue,
                                                short int force_sum_1_flag);

      // return -1 on failure
      int full_atom_spec_to_atom_index(const coot::atom_spec_t &atom_spec) const;
      // return -1 on no atom found.
      int full_atom_spec_to_atom_index(const std::string &chain,
                                       int resno,
                                       const std::string &insertion_code,
                                       const std::string &atom_name,
                                       const std::string &alt_conf) const;

      std::string name_for_display_manager() const;
      std::string dotted_chopped_name() const;
      std::string get_save_molecule_filename(const std::string &dir);
      std::string make_backup(); // returns warning message, otherwise empty string
      void save_history_file_name(const std::string &file);
      std::vector<std::string> history_filename_vec;
      std::string save_time_string;
      void restore_from_backup(int mod_index, const std::string &cwd);

      std::vector<coot::atom_spec_t> fixed_atom_specs;

      std::pair<int, mmdb::Residue *>
      find_serial_number_for_insert(int seqnum_for_new,
                                    const std::string &ins_code_for_new,
                                    const std::string &chain_id) const;

      // remove TER record from residue
      //
      void remove_TER_internal(mmdb::Residue *res_p);
      void remove_TER_on_last_residue(mmdb::Chain *chain_p);
      std::pair<bool, std::string> unused_chain_id() const;
      int append_to_molecule(const coot::minimol::molecule &water_mol);

      glm::vec4 colour_holder_to_glm(const coot::colour_holder &ch) const;

      std::pair<bool, coot::Cartesian> get_HA_unit_vector(mmdb::Residue *r) const;

      //! modify the im
      void setup_cylinder_clashes(instanced_mesh_t &im, const atom_overlaps_dots_container_t &c,
                                  float ball_size, unsigned int num_subdivisions,
                                  const std::string &molecule_name_stub) const;

      //! modify the im
      void setup_dots(instanced_mesh_t &im,
                      const atom_overlaps_dots_container_t &c,
                      float ball_size, unsigned int num_subdivisions,
                      const std::string &molecule_name_stub) const;

      // ====================== SHELX stuff ======================================

      std::pair<int, std::string> write_shelx_ins_file(const std::string &filename) const;
      int read_shelx_ins_file(const std::string &filename);
      // return the success status, 0 for fail, 1 for good.
      int add_shelx_string_to_molecule(const std::string &str);
      bool is_from_shelx_ins() const { return is_from_shelx_ins_flag; }

      void trim_atom_label_table();
      void delete_ghost_selections();
      void update_symmetry();
      bool show_symmetry;
      void delete_any_link_containing_residue(const coot::residue_spec_t &res_spec);
      void delete_link(mmdb::Link *link, mmdb::Model *model_p);

      bool sanity_check_atoms(mmdb::Manager *mol) const; // sfcalc_genmap crashes after merge of ligand.
      // Why? Something wrong with the atoms after merge?
      // Let's diagnose.... Return false on non-sane.

      // ====================== Jiggle-Fit (internal) ================================

      float fit_to_map_by_random_jiggle(mmdb::PPAtom atom_selection,
                                        int n_atoms,
                                        const clipper::Xmap<float> &xmap,
                                        float map_sigma,
                                        int n_trials,
                                        float jiggle_scale_factor,
                                        bool use_biased_density_scoring,
                                        std::vector<mmdb::Chain *> chains_for_moving);

      coot::minimol::molecule rigid_body_fit(const coot::minimol::molecule &mol_in,
                                             const clipper::Xmap<float> &xmap,
                                             float map_sigma) const;

      // ====================== init ======================================

      void init() {
         imol_no = -1; // unset
         ligand_flip_number = 0;
         bonds_box_type = UNSET_TYPE;
         is_em_map_cached_flag = false;
         xmap_is_diff_map = false;
         is_from_shelx_ins_flag = false;
         show_symmetry = false;
         default_temperature_factor_for_new_atoms = 20.0;
         original_fphis_filled = false;
         original_fobs_sigfobs_filled = false;
         original_fobs_sigfobs_fill_tried_and_failed = false;
         original_fphis_p = nullptr;
         original_fobs_sigfobs_p = nullptr;
         original_r_free_flags_p = nullptr;
         refmac_r_free_flag_sensible = false;
         use_bespoke_grey_colour_for_carbon_atoms = false;

         float rotate_colour_map_on_read_pdb = 0.24;
         bonds_colour_map_rotation = (imol_no + 1) * rotate_colour_map_on_read_pdb;
         while (bonds_colour_map_rotation > 360.0)
            bonds_colour_map_rotation -= 360.0;

      }

   public:

      // enum refine_residues_mode {SINGLE, TRIPLE, QUINTUPLE, HEPTUPLE, SPHERE, BIG_SPHERE, CHAIN, ALL};

      atom_selection_container_t atom_sel;
      // set this on reading a pdb file
      float default_temperature_factor_for_new_atoms; // direct access

      molecule_t(const std::string &name_in, int mol_no_in) : name(name_in) {init(); imol_no = mol_no_in; }
      explicit molecule_t(atom_selection_container_t asc, int imol_no_in, const std::string &name_in) : name(name_in), atom_sel(asc) {
         init();
         imol_no = imol_no_in;
         default_temperature_factor_for_new_atoms =
            util::median_temperature_factor(atom_sel.atom_selection,
                                            atom_sel.n_selected_atoms,
                                            99999.9, 0.0, false, false);
      }

      // ------------------------ close

      int close_yourself();

      // ----------------------- structure factor stuff ------------------------------------------------------

      void fill_fobs_sigfobs(); // re-reads MTZ file (currently 20210816-PE)
      // used to be a const ref. Now return the whole thing!. Caller must call
      // fill_fobs_sigfobs() directly before using this function - meh, not a good API.
      // Return a *pointer* to the data so that we don't get this hideous non-reproducable
      // crash when we access this data item after the moelcule vector has been resized
      // 20210816-PE.
      clipper::HKL_data<clipper::data32::F_sigF> *get_original_fobs_sigfobs() const {
         if (!original_fobs_sigfobs_filled) {
            std::string m("Original Fobs/sigFobs is not filled");
            throw(std::runtime_error(m));
         }
         return original_fobs_sigfobs_p;
      }

      clipper::HKL_data<clipper::data32::Flag> *get_original_rfree_flags() const {
         if (!original_fobs_sigfobs_filled) {
            std::string m("Original Fobs/sigFobs is not filled - so no RFree flags");
            throw(std::runtime_error(m));
         }
         return original_r_free_flags_p;
      }
      // use this molecules mol and the passed data to make a map for some other
      // molecule
      int sfcalc_genmap(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                        const clipper::HKL_data<clipper::data32::Flag> &free,
                        clipper::Xmap<float> *xmap_p);
      coot::util::sfcalc_genmap_stats_t
      sfcalc_genmaps_using_bulk_solvent(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                        const clipper::HKL_data<clipper::data32::Flag> &free,
                                        clipper::Xmap<float> *xmap_2fofc_p,
                                        clipper::Xmap<float> *xmap_fofc_p);
      // String munging helper function (for reading mtz files).
      // Return a pair.first string of length 0 on error to construct dataname(s).
      std::pair<std::string, std::string> make_import_datanames(const std::string &fcol,
                                                                const std::string &phi_col,
                                                                const std::string &weight_col,
                                                                int use_weights) const;
      // make these private?
      std::string refmac_fobs_col;
      std::string refmac_sigfobs_col;
      std::string refmac_mtz_filename;
      std::string refmac_r_free_col;
      std::string Refmac_fobs_col() const { return refmac_fobs_col; }
      std::string Refmac_sigfobs_col() const { return refmac_sigfobs_col; }
      std::string Refmac_mtz_filename() const { return refmac_mtz_filename; }
      bool refmac_r_free_flag_sensible;

      void associate_data_mtz_file_with_map(const std::string &data_mtz_file_name,
                                            const std::string &f_col, const std::string &sigf_col,
                                            const std::string &r_free_col);

      // ----------------------- xmap

      clipper::Xmap<float> xmap; // public because the filling function needs access

      // public access to the lock (from threads)
      static std::atomic<bool> draw_vector_sets_lock;

      // ----------------------- utils

      std::string get_name() const { return name; }
      int get_molecule_index() const { return imol_no; }
      // void set_molecule_index(int idx) { imol_no = idx; } // 20221011-PE needed?
      bool is_valid_model_molecule() const;
      bool is_valid_map_molecule() const;
      unsigned int get_number_of_atoms() const;
      mmdb::Residue *cid_to_residue(const std::string &cid) const;
      mmdb::Atom *cid_to_atom(const std::string &cid) const;
      std::pair<bool, coot::residue_spec_t> cid_to_residue_spec(const std::string &cid) const;
      std::pair<bool, coot::atom_spec_t> cid_to_atom_spec(const std::string &cid) const;
      std::vector<std::string> get_residue_names_with_no_dictionary(const coot::protein_geometry &geom) const;
      int insert_waters_into_molecule(const coot::minimol::molecule &water_mol);

      // ----------------------- model utils

      // public
      void make_bonds(protein_geometry *geom, coot::rotamer_probability_tables *rot_prob_tables_p,
                      bool draw_hydrogen_atoms_flag, bool draw_missing_loops_flag);
      // returns either the specified atom or null if not found
      mmdb::Atom *get_atom(const atom_spec_t &atom_spec) const;
      // returns either the specified residue or null if not found
      mmdb::Residue *get_residue(const coot::residue_spec_t &residue_spec) const;

      bool have_unsaved_changes() const { return save_info.have_unsaved_changes(); }
      int undo(); // 20221018-PE return status not yet useful
      int redo(); // likewise
      int write_coordinates(const std::string &file_name) const; // return 0 on OK, 1 on failure
      std::vector<coot::atom_spec_t> get_fixed_atoms() const;

      std::vector<std::string> chains_in_model() const;
      std::vector<std::pair<coot::residue_spec_t, std::string> > get_single_letter_codes_for_chain(const std::string &chain_id) const;

      residue_spec_t get_residue_closest_to(mmdb::Manager *mol, const clipper::Coord_orth &co) const;

      std::vector<std::string> get_chain_ids() const;

      // ----------------------- model bonds

      simple_mesh_t get_bonds_mesh(const std::string &mode, coot::protein_geometry *geom,
                                   bool against_a_dark_background, float bonds_width, float atom_radius_to_bond_width_ratio,
                                   int smoothness_factor,
                                   bool draw_hydrogen_atoms_flag,
                                   bool draw_missing_residue_loops);

      instanced_mesh_t get_bonds_mesh_instanced(const std::string &mode, coot::protein_geometry *geom,
                                                bool against_a_dark_background, float bonds_width, float atom_radius_to_bond_width_ratio,
                                                int smoothness_factor,
                                                bool draw_hydrogen_atoms_flag,
                                                bool draw_missing_residue_loops);

      //! If any colour rule has been set for this molecule, then we will use those. Otherwise, colorChainsScheme() will be called
      //! (and that his its internal colour-by-chain colouring scheme).
      //!
      //! the `colour_rules` is a vector of things like: ("//A", "red")
      std::vector<std::pair<std::string, std::string> > colour_rules;

      //! Add a colour rule: eg. ("//A", "red")
      void add_colour_rule(const std::string &selection, const std::string &colour_name);

      //! delete all the colour rules
      void delete_colour_rules();

      void print_colour_rules() const;

      simple_mesh_t get_molecular_representation_mesh(const std::string &cid,
                                                      const std::string &colour_scheme,
                                                      const std::string &style) const;

      simple_mesh_t get_gaussian_surface() const;

      simple_mesh_t get_chemical_features_mesh(const std::string &cid, const coot::protein_geometry &geom) const;

      bool hydrogen_atom_should_be_drawn() const { return false; } // 20221018-PE for now.
      void set_use_bespoke_carbon_atom_colour(bool state) {
         use_bespoke_grey_colour_for_carbon_atoms = state;
         // make_bonds_type_checked("set_use_bespoke_carbon_atom_colour");
      }
      void set_bespoke_carbon_atom_colour(const coot::colour_t &col) {
         bespoke_carbon_atoms_colour = col;
         // make_bonds_type_checked("set_bespoke_carbon_atom_colour");
      }

      // ----------------------- model analysis functions

      std::vector<std::string> non_standard_residue_types_in_model() const;
      std::vector<phi_psi_prob_t> ramachandran_validation(const ramachandrans_container_t &rc) const;
      // not const because it recalculates the bonds.
      simple_mesh_t get_rotamer_dodecs(protein_geometry *geom_p, rotamer_probability_tables *rpt);

      instanced_mesh_t get_rotamer_dodecs_instanced(protein_geometry *geom_p, rotamer_probability_tables *rpt);

      omega_distortion_info_container_t peptide_omega_analysis(const protein_geometry &geom,
                                                               const std::string &chain_id,
                                                               bool mark_cis_peptides_as_bad_flag) const;

      std::vector<coot::residue_spec_t> get_non_standard_residues_in_molecule() const;

      //! @return the instanced mesh for the specified ligand
      coot::instanced_mesh_t contact_dots_for_ligand(const std::string &cid, const coot::protein_geometry &geom,
                                                     unsigned int num_subdivisions) const;

      //! @return the instanced mesh for the specified molecule
      coot::instanced_mesh_t all_molecule_contact_dots(const coot::protein_geometry &geom,
                                                       unsigned int num_subdivisions) const;

      // ------------------------ model-changing functions

      int move_molecule_to_new_centre(const coot::Cartesian &new_centre);
      coot::Cartesian get_molecule_centre() const;

      int flip_peptide(const atom_spec_t &rs, const std::string &alt_conf);
      int auto_fit_rotamer(const std::string &chain_id, int res_no, const std::string &ins_code,
                           const std::string &alt_conf,
                           const clipper::Xmap<float> &xmap, const coot::protein_geometry &pg);

      std::pair<bool,float> backrub_rotamer(const std::string &chain_id, int res_no,
                                            const std::string &ins_code, const std::string &alt_conf,
                                            const clipper::Xmap<float> &xmap,
                                            const coot::protein_geometry &pg);

      // return the number of deleted atoms
      int delete_atoms(const std::vector<atom_spec_t> &atoms);
      int delete_atom(atom_spec_t &atom_spec);
      int delete_residue(residue_spec_t &residue_spec);
      int delete_residue_atoms_with_alt_conf(coot::residue_spec_t &residue_spec, const std::string &alt_conf);
      int delete_chain_using_atom_cid(const std::string &cid);
      int delete_literal_using_cid(const std::string &cid); // cid is an atom selection, e.g. containing a residue range

      std::pair<int, std::string> add_terminal_residue_directly(const residue_spec_t &spec,
                                                                const std::string &new_res_type,
                                                                const protein_geometry &geom,
                                                                const clipper::Xmap<float> &xmap);

      int add_compound(const dictionary_residue_restraints_t &monomer_restraints, const Cartesian &position,
                       const clipper::Xmap<float> &xmap, float map_rmsd);


      //! add an alternative conformation for the specified residue
      int add_alternative_conformation(const std::string &cid);

      //! add atoms to a partially-filled side chaain
      int fill_partial_residue(const residue_spec_t &res_spec, const std::string &alt_conf,
                               const clipper::Xmap<float> &xmap, const protein_geometry &geom);

      //! add atoms to a partially-filled side chaain
      int fill_partial_residues(const clipper::Xmap<float> &xmap, protein_geometry *geom);

      int mutate(const residue_spec_t &spec, const std::string &new_res_type);

      int side_chain_180(const residue_spec_t &residue_spec, const std::string &alt_conf,
                         coot::protein_geometry *geom_p); // sub functions are non-const

      int delete_side_chain(const residue_spec_t &residue_spec);

      std::string jed_flip(coot::residue_spec_t &spec, const std::string &atom_name, const std::string &alt_conf,
                           bool invert_selection, protein_geometry *geom);

      // move this up
      std::string jed_flip_internal(coot::atom_tree_t &tree,
                                    const std::vector<coot::dict_torsion_restraint_t> &interesting_torsions,
                                    const std::string &atom_name,
                                    bool invert_selection);

      // return a non-null string on a problem
      std::string jed_flip_internal(coot::atom_tree_t &tree,
                                    const dict_torsion_restraint_t &torsion,
                                    const std::string &atom_name,
                                    bool invert_selection);

      coot::minimol::molecule eigen_flip_residue(const residue_spec_t &residue_spec);

      int apply_transformation_to_atom_selection(const std::string &atom_selection_cid,
                                                 int n_atoms_in_selection,
                                                 clipper::Coord_orth &rotation_centre,
                                                 clipper::RTop_orth &rtop);

      class moved_atom_t {
      public:
         std::string atom_name;
         std::string alt_conf;
         float x, y, z;
         int index; // for fast lookup. -1 is used for "unknown"
         moved_atom_t(const std::string &a, const std::string &alt, float x_in, float y_in, float z_in) :
            atom_name(a), alt_conf(alt), x(x_in), y(y_in), z(z_in), index(-1) {}
         moved_atom_t(const std::string &a, const std::string &alt, float x_in, float y_in, float z_in, int idx) :
            atom_name(a), alt_conf(alt), x(x_in), y(y_in), z(z_in), index(idx) {}
      };

      class moved_residue_t {
      public:
         std::string chain_id;
         int res_no;
         std::string ins_code;
         std::vector<moved_atom_t> moved_atoms;
         moved_residue_t(const std::string &c, int rn, const std::string &i) : chain_id(c), res_no(rn), ins_code(i) {}
         void add_atom(const moved_atom_t &mva) {moved_atoms.push_back(mva); }
      };

      //! set new positions for the atoms in the specified residue
      int new_positions_for_residue_atoms(const std::string &residue_cid, const std::vector<moved_atom_t> &moved_atoms);

      //! set new positions for the atoms of the specified residues
      int new_positions_for_atoms_in_residues(const std::vector<moved_residue_t> &moved_residues);

      //! not for wrapping (should be private)
      int new_positions_for_residue_atoms(mmdb::Residue *residue_p, const std::vector<moved_atom_t> &moved_atoms);

      //! merge molecules - copy the atom of mols into this molecule
      //! @return the number of atoms added.
      int merge_molecules(const std::vector<mmdb::Manager *> &mols);

      //! My ligands don't jiggle-jiggle...
      //!
      //! Hey, what do you know, they actually do.
      float fit_to_map_by_random_jiggle(const residue_spec_t &res_spec, const clipper::Xmap<float> &xmap, float map_rmsd,
                                        int n_trials, float translation_scale_factor);

      int cis_trans_conversion(const std::string &atom_cid, mmdb::Manager *standard_residues_mol);

      //! @return the success status
      int replace_fragment(atom_selection_container_t asc);

      // ----------------------- merge molecules

      // merge molecules helper functions

      bool is_het_residue(mmdb::Residue *residue_p) const;
      // return state, max_resno + 1, or 0, 1 of no residues in chain.
      //
      std::pair<short int, int> next_residue_number_in_chain(mmdb::Chain *w,
                                                             bool new_res_no_by_hundreds=false) const;

      mmdb::Residue *copy_and_add_residue_to_chain(mmdb::Chain *this_model_chain,
                                                   mmdb::Residue *add_model_residue,
                                                   bool new_resno_by_hundreds_flag=true);
      void copy_and_add_chain_residues_to_chain(mmdb::Chain *new_chain, mmdb::Chain *this_molecule_chain);
      std::vector<std::string> map_chains_to_new_chains(const std::vector<std::string> &adding_model_chains,
                                                        const std::vector<std::string> &this_model_chains) const;
      // that's too complicated for try_add_by_consolidation(), we just want this:
      std::string suggest_new_chain_id(const std::string &current_chain_id) const;
      std::pair<bool, std::vector<std::string> > try_add_by_consolidation(mmdb::Manager *adding_mol);
      bool merge_molecules_just_one_residue_homogeneous(atom_selection_container_t molecule_to_add);
      bool merge_molecules_just_one_residue_at_given_spec(atom_selection_container_t molecule_to_add,
                                                          residue_spec_t target_spec);

      // return success status and spec if new residue if possible.
      std::pair<bool, coot::residue_spec_t> merge_ligand_to_near_chain(mmdb::Manager *mol);

      // return success status and spec if new residue if possible.
      std::pair<int, std::vector<merge_molecule_results_info_t> >
      merge_molecules(const std::vector<atom_selection_container_t> &add_molecules);

      // ----------------------- refinement

      coot::extra_restraints_t extra_restraints;
      //! refinement tool
      std::vector<mmdb::Residue *> select_residues(const residue_spec_t &spec, const std::string &mode) const;
      //! resno_start and resno_end are inclusive
      std::vector<mmdb::Residue *> select_residues(const std::string &chain_id, int resno_start, int resno_end) const;

      int refine_direct(std::vector<mmdb::Residue *> rv, const std::string &alt_loc, const clipper::Xmap<float> &xmap,
                        float map_weight, const coot::protein_geometry &geom, bool refinement_is_quiet);

      // ----------------------- map functions

      // return -1.1 on not-a-map
      float get_map_rmsd_approx() const;
      int write_map(const std::string &file_name) const;
      void set_map_is_difference_map(bool flag);
      bool is_difference_map_p() const;

      // changes the internal map mesh holder (hence not const)
      coot::simple_mesh_t get_map_contours_mesh(clipper::Coord_orth position, float radius, float contour_level);

      //! The container class for an interesting place.
      //!
      //! This documentation doesn't work and I don't know why.
      class interesting_place_t {
      public:
         //! Feature
         std::string feature_type;
         //! Residue specifier
         residue_spec_t residue_spec; // use this for sorting a combination of interesting_place_t types.
         //! Position
         float x, y, z;
         //! button label
         std::string button_label;
         //! actual value of the feature
         float feature_value; // e.g. peak-height (not all "interesting" feature values can be captured by a float of course)
         //! synthetic badness (for "score by badness")
         float badness; // a nubmer between 100.0 and 0.0 (inclusive) if it's negative then it's not set.
         //! constructor
         interesting_place_t() {}
         //! constructor
         interesting_place_t(const std::string &ft, const residue_spec_t &rs, const clipper::Coord_orth &pt, const std::string &bl) :
            feature_type(ft), residue_spec(rs), button_label(bl) {
            x = pt.x(); y = pt.y(); z = pt.z();
            feature_value = -1; // something "unset"
            badness = -1.1; // "unset"
         }
         //! constructor
         interesting_place_t(const std::string &ft, const clipper::Coord_orth &pt, const std::string &bl) : feature_type(ft), button_label(bl) {
            x = pt.x(); y = pt.y(); z = pt.z();
            feature_value = -1; // something "unset"
            badness = -1.1; // "unset"
         }
         //! internal to libcootapi function to set the values
         void set_feature_value(const float &f) { feature_value = f; }
         void set_badness_value(const float &b) { badness = b; }
      };

      //! difference maps peaks class
      class difference_map_peaks_info_t {
      public:
         clipper::Coord_orth pos;
         float peak_height; // nrmsd
         // maybe other useful stuff here in future
         difference_map_peaks_info_t(const clipper::Coord_orth &p, float ph) : pos(p), peak_height(ph) {}
      };

      // the molecule is passed so that the peaks are placed around the protein
      std::vector<interesting_place_t> difference_map_peaks(mmdb::Manager *mol, float n_rmsd) const;

   };
}


#endif // COOT_MOLECULE_HH
