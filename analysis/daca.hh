/*
 * analysis/data.hh
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include <vector>
#include <utility>
#include <string>
#include <unordered_set>

#include <clipper/core/coords.h>
#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"


namespace coot {

   class daca {

   public:

      class box_index_t {
      public:
         int idx_x;
         int idx_y;
         int idx_z;
         float box_width;
         box_index_t(int ii, int jj, int kk) : idx_x(ii), idx_y(jj), idx_z(kk) {
            box_width = 1.0;
         }
         box_index_t(const clipper::Coord_orth &pos);
         bool operator<(const box_index_t &other) const;
         clipper::Coord_orth coord_orth() const;
         // friend std::ostream &operator<<(std::ostream &s, const box_index_t &bi);
         float d() const;
         float d_squared() const;
      };


      double get_radius(const std::string &ele) const;
      // std::ostream &operator<<(std::ostream &s, const box_index_t &bi);
      std::vector<std::pair<mmdb::Atom *, float> >
      solvent_exposure_old_version(int selhnd, mmdb::Manager *mol) const;
      std::vector<std::pair<mmdb::Residue *, float> >
      solvent_exposure_old_version_v2(mmdb::Manager *mol, bool side_chain_only=true) const;
      std::vector<std::pair<mmdb::Residue *, float> >
      solvent_exposure(mmdb::Manager *mol, bool side_chain_only=true) const;

   private:
      enum mode_t {REFERENCE, ANALYSIS};
      // the atoms in the reference fragments have a particular order and (for each fragmemt)
      // are centred on the origin.
      std::map<std::string, std::vector<std::vector<clipper::Coord_orth> > > reference_fragments;
      void fill_reference_fragments();
      bool boxes_have_been_resized = false;
      std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > > boxes; // reference
      std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > > boxes_for_testing; // this pdb
      std::vector<mmdb::Residue *> helical_residues; // fill this before calling calculate_daca()
      std::map<float, float> envelope_distance_map; // the key is squared values
      struct quality_index_model_t {
         float intercept;
         float slope;
         float stddev_residual;
      };
      std::map<std::string, quality_index_model_t> quality_index_models; // key: "ALA-helix" etc.
      void fill_helix_flags(mmdb::Model *model_p, mmdb::Manager *mol);
      std::pair<bool, clipper::RTop_orth>  // return status also, as this can fail.
      get_frag_to_reference_rtop(const std::string &res_name,
                                 const unsigned int &frag_idx,
                                 const std::vector<mmdb::Atom *> &fragment_atoms) const;
      void add_to_box(mode_t mode,
                      const std::string &res_name,
                      bool residue_is_helical_flag,
                      unsigned int frag_index,
                      const box_index_t &box_index,
                      const std::string &atom_type,
                      unsigned int counts=1);
      int get_reference_counts(const std::string &res_name,
                               bool residue_is_helical_flag,
                               unsigned int frag_index,
                               const box_index_t &box_index,
                               const std::string &atom_type) const;
      std::vector<std::vector<mmdb::Atom *> > get_daca_fragments(mmdb::Residue *reference_residue_p) const;
      std::vector<std::pair<mmdb::Atom *, std::string> >
         make_typed_atoms(mmdb::Model *model_p, const protein_geometry &geom) const;
      std::vector<std::pair<mmdb::Atom *, std::string> >
         make_symmetry_typed_atoms(mmdb::Manager *mol,
                                   mmdb::Model *model_p,
                                   const protein_geometry &geom,
                                   float expansion_radius,
                                   std::vector<mmdb::Atom *> *symm_atom_store_p) const;

      float calculate_daca(mmdb::Residue *reference_residue_p,
                         const std::vector<std::pair<mmdb::Atom *, std::string> > &typed_atoms,
                         mode_t mode);
      bool atom_is_close_to_a_residue_atom(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const;
      bool atom_is_neighbour_mainchain(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const;
      void debug_boxes(const std::string &debug_prefix="") const;
      void compare_boxes() const;
      void presize_boxes(mode_t mode=REFERENCE);
      void normalize();
      void normalize_v2();
      void envelope();
      void smooth();

   public:
      daca() { fill_reference_fragments(); boxes_have_been_resized = false; }
      std::vector<std::vector<std::string> > atom_names_for_fragments(const std::string &res_name) const;
      float gompertz_scale(const float &dist);
      void write_tables_using_reference_structures_from_dir(const std::string &input_pdb_files_dir_name,
                                                            const std::string &output_tables_dir);
      void read_tables(const std::string &dir);
      void read_many_tables(const std::vector<std::string> &dirs);
      void write_tables(const std::string &dir_name) const;
      void read_quality_index_table(const std::string &file_name);
      void score_molecule(const std::string &pdb_file_name);
      void make_data_for_figure_2(const std::string &pdb_dir);
      void make_quality_index_table(const std::string &input_table_file,
                                    const std::string &output_file_name);
      void cook();
      //! Run REFERENCE and ANALYSIS passes on the same PDB file.
      //! Returns a pair: (n_total_contacts, n_misses).
      //! For a correct implementation n_misses should be 0.
      std::pair<int, int> self_test(const std::string &pdb_file_name);
      //! Build the reference from the original PDB, then perturb atom positions
      //! by a random offset in [-perturbation, +perturbation] along each axis
      //! and score the perturbed model. Returns (n_total_contacts, n_misses).
      std::pair<int, int> self_test_perturbed(const std::string &pdb_file_name,
                                              float perturbation);

   };


}
