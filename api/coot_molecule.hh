#ifndef COOT_MOLECULE_HH
#define COOT_MOLECULE_HH

#include <utility>
#include <atomic>

#include <clipper/core/xmap.h>
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/coot-rama.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"
#include "simple-mesh.hh"
#include "ghost-molecule-display.hh"

#include "density-contour/CIsoSurface.h"
#include "gensurf.hh"
#include "coot-utils/coot-shelx.hh"


namespace coot {

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
         std::pair<time_t, unsigned int> last_saved;
         unsigned int modification_index;
         molecule_save_info_t() : last_saved(std::make_pair(0,0)), modification_index(0) {}
         void new_modification() {
            modification_index++;
         }
         void made_a_save() {
            // this is called when the server says that it has saved the file
            last_saved.first = time(nullptr);
            last_saved.second = modification_index;
         }
         bool have_unsaved_changes() const {
            return modification_index > last_saved.second;
         }
      };

      molecule_save_info_t save_info;

      int imol_no; // this molecule's index in the container vector
      std::string name;
      bool is_from_shelx_ins_flag;
      ShelxIns shelxins;

      // private
      void makebonds(coot::protein_geometry *geom, coot::rotamer_probability_tables *rotamer_tables_p, std::set<int> &no_bonds_to_these_atoms);

#if defined __has_builtin
#if __has_builtin (__builtin_FUNCTION)
      void make_bonds_type_checked(coot::protein_geometry *geom, coot::rotamer_probability_tables *rot_prob_tables_p, const char *s = __builtin_FUNCTION());
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, const char *s = __builtin_FUNCTION());
#else
      void make_bonds_type_checked(coot::protein_geometry *geom, const char *s = 0);
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, const char *s =0);
#endif
#else // repeat above
      void make_bonds_type_checked(coot::protein_geometry *geom, const char *s = 0);
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, const char *s =0);
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
                                      bool do_rota_markup=false,
                                      coot::rotamer_probability_tables *rotamer_tables_p = nullptr,
                                      bool force_rebonding=true);
      void make_ca_bonds();

      void update_map_triangles(float radius, coot::Cartesian centre, float contour_level);

      bool is_EM_map() const;
      short int is_em_map_cached_flag; // -1 mean unset (so set it, 0 means no, 1 means yes)
      short int is_em_map_cached_state(); // set is_em_map_cached_flag if not set
      coot::ghost_molecule_display_t map_ghost_info;

      bool xmap_is_diff_map;
      bool has_xmap() const { return is_valid_map_molecule(); }

      void clear_draw_vecs();
      void clear_diff_map_draw_vecs();
      std::vector<coot::density_contour_triangles_container_t> draw_vector_sets;
      std::vector<coot::density_contour_triangles_container_t> draw_diff_map_vector_sets;
      std::vector<std::pair<int, TRIANGLE> > map_triangle_centres; // with associated mid-points and indices

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

      void make_backup();

      void trim_atom_label_table();
      void delete_ghost_selections();
      void update_symmetry();
      bool show_symmetry;

   public:

      atom_selection_container_t atom_sel;
      molecule_t() {}
      explicit molecule_t(atom_selection_container_t asc, int imol_no_in) : atom_sel(asc) {
         init();
         imol_no = imol_no_in;
      }

      void init() { // add imol_no here?
         bonds_box_type = UNSET_TYPE;
         is_em_map_cached_flag = false;
         xmap_is_diff_map = false;
         is_from_shelx_ins_flag = false;
         show_symmetry = false;
      }

      clipper::Xmap<float> xmap; // public because the filling function needs access

      // public access to the lock (from threads)
      static std::atomic<bool> draw_vector_sets_lock;

      // utils

      std::string get_name() const { return name; }
      int get_molecule_index() const { return imol_no; }
      // void set_molecule_index(int idx) { imol_no = idx; } // 20221011-PE needed?
      bool is_valid_model_molecule() const;
      bool is_valid_map_molecule() const;
      std::pair<bool, coot::residue_spec_t> cid_to_residue_spec(const std::string &cid);

      // model utils

      // public
      void make_bonds(coot::protein_geometry *geom, coot::rotamer_probability_tables *rot_prob_tables_p);
      // returns either the specified atom or null if not found
      mmdb::Atom *get_atom(const coot::atom_spec_t &atom_spec) const;
      // returns either the specified residue or null if not found
      mmdb::Residue *get_residue(const coot::residue_spec_t &residue_spec) const;

      bool have_unsaved_changes() const { return save_info.have_unsaved_changes(); }

      // model analysis functions

      std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > ramachandran_validation() const;
      // not const because it recalculates the bonds.
      coot::simple_mesh_t get_rotamer_dodecs(coot::protein_geometry *geom_p,
                                             coot::rotamer_probability_tables *rpt);

      // model-changing functions

      int flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf);
      int auto_fit_rotamer(const std::string &chain_id, int res_no, const std::string &ins_code,
                           const std::string &alt_conf,
                           const clipper::Xmap<float> &xmap, const coot::protein_geometry &pg);

      std::pair<bool,float> backrub_rotamer(const std::string &chain_id, int res_no,
                                            const std::string &ins_code, const std::string &alt_conf,
                                            const coot::protein_geometry &pg);
      // return the number of deleted atoms
      int delete_atoms(const std::vector<coot::atom_spec_t> &atoms);
      int delete_atom(coot::atom_spec_t &atom_spec);
      int delete_residue(coot::residue_spec_t &residue_spec);
      int delete_residue_atoms_with_alt_conf(coot::residue_spec_t &residue_spec, const std::string &alt_conf);

      // map functions

      int writeMap(const std::string &file_name) const;

      // changes the internal map mesh holder (hence not const)
      coot::simple_mesh_t get_map_contours_mesh(clipper::Coord_orth position, float radius, float contour_level);

   };
}


#endif // COOT_MOLECULE_HH
