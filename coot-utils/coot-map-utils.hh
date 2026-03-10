/* coot-utils/coot-map-utils.hh
 *
 * Copyright 2004, 2005, 2006, 2007 The University of York
 * Copyright 2015 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
 * 02110-1301, USA.
 */

#ifndef COOT_MAP_UTILS_HH
#define COOT_MAP_UTILS_HH

#include <map>

#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include <clipper/core/hkl_data.h>
#include <clipper/contrib/sfcalc_obs.h>
#include "coot-coord-utils.hh"
#include "coot-density-stats.hh"
#include <mmdb2/mmdb_manager.h>
#include "amp-reso.hh"

namespace coot {

   namespace util {

      clipper::RTop_orth make_rtop_orth_from(mmdb::mat44 mat);

      // cubic interpolation
      float density_at_point(const clipper::Xmap<float> &map_in,
                             const clipper::Coord_orth &co);
      float density_at_point_by_cubic_interp(const clipper::NXmap<float> &map_in,
                                             const clipper::Coord_map &cm);
      // linear interpolation (faster) use for jiggle-fit of chains and the like
      float density_at_point_by_linear_interpolation(const clipper::Xmap<float> &map_in,
                                                     const clipper::Coord_orth &co);
      // nearest grid point - faster yet
      float density_at_point_by_nearest_grid(const clipper::Xmap<float> &map_in,
                                             const clipper::Coord_orth &co);
      // NXmap versions
      float density_at_point_by_nearest_grid(const clipper::NXmap<float> &nxmap,
                                             const clipper::Coord_orth &co);
      float density_at_point_by_linear_interp(const clipper::NXmap<float> &nxmap,
                                              const clipper::Coord_orth &co);

      float density_at_map_point(const clipper::Xmap<float> &map_in,
                                 const clipper::Coord_map &cg);

      clipper::Grad_orth<double> gradient_at_point(const clipper::Xmap<float> &map_in,
                                                   const clipper::Coord_orth &co);

      // return a variance of -1 on error.
      std::pair<float, float> mean_and_variance(const clipper::Xmap<float> &xmap);

      density_stats_info_t density_around_point(const clipper::Coord_orth &point,
                                                const clipper::Xmap<float> &xmap,
                                                float d);

      // This is a console/testing function.  Should not be used in a
      // real graphics program.  Use instead import_map_from() with a
      // precalculated map.
      //
      // return 1 if map is filled, 0 if not (e.g. mtz file not found).
      //
      bool map_fill_from_mtz(clipper::Xmap<float> *xmap,
                             std::string mtz_file_name,
                             std::string f_col,
                             std::string phi_col,
                             std::string weight_col,
                             short int use_weights,
                             float sampling_rate=1.5);

      bool map_fill_from_mtz(clipper::Xmap<float> *xmap,
                             std::string mtz_file_name,
                             std::string f_col,
                             std::string phi_col,
                             std::string weight_col,
                             short int use_weights,
                             float reso_limit_high,
                             short int use_reso_limit_high,
                             float sampling_rate=1.5);

      // needed by above:
      void filter_by_resolution(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata,
                                const float &reso_low,
                                const float &reso_high);


      // return the max gridding of the map, e.g. 0.5 for a 1A map.
      //
      // return the maximum Angstrom/grid of the given map.
      //
      float max_gridding(const clipper::Xmap<float> &xmap);


      // The sum of the density at the atom centres, optionally
      // weighted by atomic number.
      //
      float map_score(mmdb::PPAtom atom_selection,
                      int n_selected_atoms,
                      const clipper::Xmap<float> &xmap,
                      short int with_atomic_weighting);

      // The sum of the density at the atom centres, weighted by occupancy
      //
      float map_score(std::vector<mmdb::Atom *> atoms, const clipper::Xmap<float> &xmap);

      float map_score_atom(mmdb::Atom *atom, const clipper::Xmap<float> &xmap);

      float map_score_by_residue_specs(mmdb::Manager *mol,
                                       const std::vector<residue_spec_t> &res_specs,
                                       const clipper::Xmap<float> &xmap,
                                       bool main_chain_only_flag = false);

      clipper::Xmap<float> sharpen_blur_map(const clipper::Xmap<float> &xmap_in, float b_factor);

      // sharpen/blur self
      void sharpen_blur_map(clipper::Xmap<float> *xmap, float b_factor);

      clipper::Xmap<float> sharpen_blur_map_with_resample(const clipper::Xmap<float> &xmap_in, float b_factor,
                                                          float resample_factor);

      clipper::Xmap<float> sharpen_blur_map_with_reduced_sampling(const clipper::Xmap<float> &xmap_in, float b_factor,
                                                                  float resample_factor);

      // pass a pointer to a vector of maps that has the same size as the number of B-factors
      //
      void multi_sharpen_blur_map(const clipper::Xmap<float> &xmap_in,
                                  const std::vector<float> &b_factors,
                                  std::vector<clipper::Xmap<float> > *maps_p);

      // not sure if this works ATM
      clipper::Xmap<float> sharpen_map(const clipper::Xmap<float> &xmap_in,
                                       float sharpen_factor);

      //! map molecule centre
      class map_molecule_centre_info_t {
      public:
         //! success flag
         bool success;
         //! new centre
         clipper::Coord_orth updated_centre;
         //! suggested contour level
         float suggested_contour_level;
         //! the suggested radius
         float suggested_radius;
         //! sum of densities - for whatever use that may be.
         double sum_of_densities; // for scoring origins
         map_molecule_centre_info_t() {
            success = false;
            sum_of_densities = -1;
            suggested_contour_level = 0.0;
            suggested_radius = -1;
         }
      };

      map_molecule_centre_info_t map_molecule_centre(const clipper::Xmap<float> &xmap);

      map_molecule_centre_info_t map_molecule_recentre_from_position(const clipper::Xmap<float> &xmap,
                                                                     const clipper::Coord_orth &current_centre);

      // if n_bins is -1, let the function decide how many bins
      //
      // actually, we return bins of amplitude squares.
      std::vector<amplitude_vs_resolution_point>
      amplitude_vs_resolution(const clipper::Xmap<float> &xmap_in, int n_bins = -1);

      // pass a flag for the resolution limit saying if this limit should be used.
      // rule of thumb: low resolution limit 0.12
      float b_factor(const std::vector<amplitude_vs_resolution_point> &fsqrd_data,
                     std::pair<bool, float> reso_low_invresolsq  = std::pair<bool, float>(false, -1),
                     std::pair<bool, float> reso_higy_invresolsq = std::pair<bool, float>(false, -1));

      clipper::Xmap<float> transform_map(const clipper::Xmap<float> &xmap_in,
                                         const clipper::Spacegroup &new_space_group,
                                         const clipper::Cell &new_cell,
                                         const clipper::RTop_orth &rtop,
                                         const clipper::Coord_orth &about_pt,
                                         float box_size);

      clipper::Grid_sampling suggested_grid_sampling(const clipper::Grid_sampling &orig_sampling,
                                                     const clipper::Cell &orig_cell,
                                                     const clipper::Spacegroup &new_space_group,
                                                     const clipper::Cell &new_cell);

      clipper::Xmap<float> laplacian_transform(const clipper::Xmap<float> &xmap_in);

      std::vector<float> density_map_points_in_sphere(clipper::Coord_orth pt, float radius,
                                                      const clipper::Xmap<float> &xmap_in);

      // pass a negative atom_selection to build an atom map for the whole molecule
      //
      clipper::Xmap<float> calc_atom_map(mmdb::Manager *mol,
                                         int atom_selection_handle,
                                         const clipper::Cell &cell,
                                         const clipper::Spacegroup &space_group,
                                         const clipper::Grid_sampling &sampling);

      clipper::Xmap<float> mask_map(const clipper::Xmap<float> &xmap_in,
                                    const std::vector<mmdb::Residue *> &neighbs);

      clipper::Xmap<float> make_map_mask(const clipper::Spacegroup &space_group,
                                         const clipper::Cell &cell,
                                         const clipper::Grid_sampling &grid_sampling,
                                         mmdb::Manager *mol,
                                         int atom_selection_handle,
                                         float radius,
                                         float smooth);

      // return a number less than -1 on badness
      // (perhaps this should return the atom map and the mask map)
      //
      // 0: all-atoms
      // 1: main-chain atoms if is standard amino-acid, else all atoms
      // 2: side-chain atoms if is standard amino-acid, else all atoms
      // 3: side-chain atoms-exclusing CB if is standard amino-acid, else all atoms
      //
      float map_to_model_correlation(mmdb::Manager *mol,
                                     const std::vector<residue_spec_t> &specs_for_correl,
                                     const std::vector<residue_spec_t> &specs_for_masking_neighbs,
                                     unsigned short int atom_mask_mode,
                                     float atom_radius, // for masking
                                     const clipper::Xmap<float> &xmap_from_sfs);

      density_correlation_stats_info_t
      map_to_model_correlation_stats(mmdb::Manager *mol,
                                     const std::vector<residue_spec_t> &specs_for_correl,
                                     const std::vector<residue_spec_t> &specs_for_masking_neighbs,
                                     unsigned short int atom_mask_mode,
                                     float atom_radius, // for masking
                                     const clipper::Xmap<float> &xmap_from_sfs,
                                     map_stats_t map_stats_flag);

      // the second of the pair contains the correlation for the given residue spec.
      //
      std::vector<std::pair<residue_spec_t, float> >
      map_to_model_correlation_per_residue(mmdb::Manager *mol,
                                           const std::vector<residue_spec_t> &specs,
                                           unsigned short int atom_mask_mode,
                                           float atom_radius, // for masking
                                           const clipper::Xmap<float> &xmap_from_sfs);

      std::map<coot::residue_spec_t, density_stats_info_t>
      map_to_model_correlation_stats_per_residue(mmdb::Manager *mol,
                                                 const std::vector<residue_spec_t> &specs,
                                                 unsigned short int atom_mask_mode,
                                                 float atom_radius, // for masking
                                                 const clipper::Xmap<float> &xmap);

      // n_residues_per_run should be an odd number more than 2 (say 11)
      //
      std::pair<std::map<coot::residue_spec_t, density_correlation_stats_info_t>, std::map<coot::residue_spec_t, density_correlation_stats_info_t> >
      map_to_model_correlation_stats_per_residue_run(mmdb::Manager *mol,
                                                     const std::string &chain_id,
                                                     const clipper::Xmap<float> &xmap,
                                                     unsigned int n_residues_per_run,
                                                     bool exclude_CON,
                                                     float atom_mask_radius=2.8,
                                                     float NOC_mask_radius=1.8); // optimized on strepavidin

      // helper
      std::pair<clipper::Coord_frac, clipper::Coord_frac>
      find_struct_fragment_coord_fracs_v2(const std::pair<clipper::Coord_orth, clipper::Coord_orth> &selection_extents,
                                          const clipper::Cell &cell);


      // which uses these for map statistics aggregation
      class map_stats_holder_helper_t {
      public:
         double sum_x;
         double sum_x_squared;
         double sum_y;
         double sum_y_squared;
         double sum_xy;
         int n;
         map_stats_holder_helper_t() {
            sum_x = 0;
            sum_x_squared = 0;
            sum_y = 0;
            sum_y_squared = 0;
            sum_xy = 0;
            n = 0;
         }
         void add_xy(const double &x, const double &y) {
            sum_x += x;
            sum_y += y;
            sum_x_squared += x*x;
            sum_y_squared += y*y;
            sum_xy += x*y;
            n++;
         }
      };


      // should this be here, or is it heavy?
      std::vector<std::pair<double, double> >
      qq_plot_for_map_over_model(mmdb::Manager *mol,
                                 const std::vector<coot::residue_spec_t> &specs,
                                 const std::vector<coot::residue_spec_t> &nb_residues,
                                 int atom_mask_mode,
                                 const clipper::Xmap<float> &xmap);

      // return a map and its standard deviation.  scale is applied to
      // map_in_2 before substraction.
      std::pair<clipper::Xmap<float>, float>
      difference_map(const clipper::Xmap<float> &xmap_in_1,
                     const clipper::Xmap<float> &xmap_in_2,
                     float map_scale);

      // like above, but average
      clipper::Xmap<float>
      average_map(const std::vector<std::pair<clipper::Xmap<float>, float> > &maps_and_scales_vec);


      // like above, but variance
      clipper::Xmap<float> variance_map(const std::vector<std::pair<clipper::Xmap<float>, float> > &maps_and_scales_vec);

      // Similar to average_map() but modify the map, don't return a new one
      // Also this function presumes that the maps have the same gridding (which makes it faster)
      //
      void
      regen_weighted_map(clipper::Xmap<float> *xmap_in,
                         const std::vector<std::pair<clipper::Xmap<float> *, float> > &maps_and_scales_vec);

      // Spin the torsioned atom round the rotatable bond and find the
      // orientation (in degrees) from the current position that is in
      // the highest density.
      //
      // return a torsion
      std::pair<float, float> spin_search(const clipper::Xmap<float> &xmap, mmdb::Residue *res, coot::torsion tors);

      // Return a map that is a copy of the given map with interpolation,
      // with grid spacing at most 0.5A (by generated by integer scaling
      // factor of the input map)
      //
      clipper::Xmap<float> reinterp_map_fine_gridding(const clipper::Xmap<float> &xmap);

      // return a map that is a copy of the given map with interpolation
      // e.g. a sampling_multiplier of 2 will double the number of grid points (in each direction).
      //
      clipper::Xmap<float> reinterp_map(const clipper::Xmap<float> &xmap_in, float sampling_multiplier);


      // make a copy of map_in, but in the cell, spacegroup and gridding of reference_map
      clipper::Xmap<float> reinterp_map(const clipper::Xmap<float> &xmap_in,
                                        const clipper::Xmap<float> &reference_xmap);

      // negative becomes positive and positive becomes negative.
      // Apply an offset so that most of the map is above zero.
      //
      void reverse_map(clipper::Xmap<float> *xmap_p);

      class map_fragment_info_t {
         // sans recentre at origin
         void init(const clipper::Xmap<float> &xmap,
                   const clipper::Coord_orth &centre,
                   float radius);
         //
         void init_making_map_centred_at_origin(const clipper::Xmap<float> &xmap,
                                                const clipper::Coord_orth &centre,
                                                float radius);
         float box_radius_a_internal;
         float box_radius_b_internal;
         float box_radius_c_internal;
      public:
         map_fragment_info_t(const clipper::Xmap<float> &xmap,
                             const clipper::Coord_orth &centre,
                             float radius, bool centre_at_origin = false);
         clipper::Xmap<float> xmap;
         clipper::Coord_grid offset;
         // transfer xmap (small, at origin) into xmap_p (big)
         void unshift(clipper::Xmap<float> *xmap_p, const clipper::Coord_orth &centre);
         void simple_origin_shift(const clipper::Xmap<float> &ip_xmap,
                                  const clipper::Coord_orth &centre,
                                  float radius);
         clipper::Grid_map make_grid_map(const clipper::Xmap<float> &input_xmap,
                                         const clipper::Coord_orth &centre) const;
      };


      //
      class simple_residue_triple_t {
      public:
         mmdb::Residue *this_residue;
         mmdb::Residue *next_residue;
         mmdb::Residue *prev_residue;
         std::string alt_conf;
         simple_residue_triple_t() {
            this_residue = 0;
            prev_residue = 0;
            next_residue = 0;
         }
         simple_residue_triple_t(mmdb::Residue *this_residue_in,
                                 mmdb::Residue *prev_residue_in,
                                 mmdb::Residue *next_residue_in,
                                 const std::string &alt_conf_in) : alt_conf(alt_conf_in) {
            this_residue = this_residue_in;
            prev_residue = prev_residue_in;
            next_residue = next_residue_in;
         }
      };

      //
      class residue_triple_t {
      public:
         mmdb::Residue *this_residue;
         mmdb::Residue *next_residue;
         mmdb::Residue *prev_residue;
         std::string alt_conf;
         residue_triple_t() {
            this_residue = 0;
            prev_residue = 0;
            next_residue = 0;
         }
         residue_triple_t(mmdb::Residue *this_residue_in,
                          mmdb::Residue *prev_residue_in,
                          mmdb::Residue *next_residue_in,
                          const std::string &alt_conf_in) : alt_conf(alt_conf_in) {
            this_residue = this_residue_in;
            prev_residue = prev_residue_in;
            next_residue = next_residue_in;
         }
         ~residue_triple_t() {
            delete this_residue;
            delete next_residue;
            delete prev_residue;
         }
         residue_triple_t deep_copy() {
            mmdb::Residue *this_residue_cp = deep_copy_this_residue(this_residue);
            mmdb::Residue *prev_residue_cp = deep_copy_this_residue(this_residue);
            mmdb::Residue *next_residue_cp = deep_copy_this_residue(this_residue);
            return residue_triple_t(this_residue_cp,
                                    prev_residue_cp,
                                    next_residue_cp,
                                    alt_conf);
         }
      };

      class backrub_residue_triple_t : public residue_triple_t {

      public:
         //  Note to self, the residue copy may need the deep_copy
         //  that does the atom index transfer too.

         backrub_residue_triple_t(mmdb::Residue *this_residue_in,
                                  mmdb::Residue *prev_residue_in,
                                  mmdb::Residue *next_residue_in,
                                  const std::string &alt_conf_in) : residue_triple_t(this_residue_in,
                                                                                     prev_residue_in,
                                                                                     next_residue_in,
                                                                                     alt_conf_in) {
            trim_this_residue_atoms();
            trim_prev_residue_atoms();
            trim_next_residue_atoms();
         }

         // Delete atoms that don't have this alt conf (or "" alt conf).
         void trim_this_residue_atoms();
         // As above, and also delete all atoms that are not C or O
         void trim_prev_residue_atoms();
         // As trim_this_residue_atoms, and also delete all atoms that are not N or H.
         void trim_next_residue_atoms();
         void trim_residue_atoms_generic(mmdb::Residue *residue_p,
                                         std::vector<std::string> keep_atom_vector,
                                         bool use_keep_atom_vector);

      };


      class map_ref_triple_t {
      public:
         double dist_sq;
         clipper::Xmap<float>::Map_reference_coord iw;
         float density;
         map_ref_triple_t(const double &d_in,
                          const clipper::Xmap<float>::Map_reference_coord &iw_in,
                          const float &den_in) {
            dist_sq = d_in;
            iw = iw_in;
            density = den_in;
         }
         map_ref_triple_t() {}
         bool operator<(const map_ref_triple_t &mrt) const {
            return (mrt.dist_sq < dist_sq);
         }
      };


      class segment_map {
         enum {UNASSIGNED = -1, TOO_LOW = -2 };
         // sorting function used by above
         static bool compare_density_values_map_refs(const std::pair<clipper::Xmap_base::Map_reference_index, float> &v1,
                                                     const std::pair<clipper::Xmap_base::Map_reference_index, float> &v2);

         // change values in segmented_map
         //
         void flood_fill_segmented_map(clipper::Xmap<std::pair<bool, int> > *segmented_map,
                                       const clipper::Xmap<float> &xmap,  // for Neighbours
                                       const clipper::Coord_grid &seed_point,
                                       int from_val, int to_val);

         // return a vector of grid points as we trace the route to the
         // local peak from start_point by steepest ascent.
         //
         std::vector<clipper::Coord_grid> path_to_peak(const clipper::Coord_grid &start_point,
                                                       const clipper::Xmap<float> &xmap_new);
         static bool sort_segment_vec(const std::pair<int, int> &a,
                                      const std::pair<int, int> &b);
         int find_biggest_segment(const std::map<int, std::vector<clipper::Coord_grid> > &segment_id_map,
                                  const std::map<int, int> &segment_id_counter_map) const;
         // test function
         int find_smallest_segment(const std::map<int, std::vector<clipper::Coord_grid> > &segment_id_map,
                                   const std::map<int, int> &segment_id_counter_map) const;
         void resegment_watershed_points(clipper::Xmap<int> *xmap_int,
                                         const clipper::Xmap<float> &xmap) const;

      public:
         segment_map() {};
         // Return the number of segments and the segmented map.
         //
         // -1 means no segment. low_level is the level below which
         // segmentation should not occur (don't make blobs in the
         // noise).
         //
         // This is Pintilie flooding (with extra watershed remapping)
         //
         std::pair<int, clipper::Xmap<int> > segment(const clipper::Xmap<float> &xmap_in, float low_level);

         // This is the flood-down method
         //
         std::pair<int, clipper::Xmap<int> > segment_emsley_flood(const clipper::Xmap<float> &xmap_in,
                                                                  float low_level);

         // multi-scale segmentation.  Return a segmented map.
         //
         std::pair<int, clipper::Xmap<int> > segment(const clipper::Xmap<float> &xmap_in,
                                                     float low_level,
                                                     float b_factor_increment, // per round
                                                     int n_rounds);
      };

      class soi_variance {
         const clipper::Xmap<float> &xmap;
         clipper::Xmap<float> make_variance_map() const;
         clipper::Xmap<float> solvent_treatment_map() const;
         clipper::Xmap<float> protein_treatment_map() const;
         static void apply_variance_values(clipper::Xmap<float> &variance_map,
                                           const clipper::Xmap<float> &xmap,
                                           const std::vector<clipper::Coord_grid> &soi_gps,
                                           const std::vector<clipper::Xmap_base::Map_reference_index> &grid_indices);
      public:
         soi_variance(const clipper::Xmap<float> &xmap_in) : xmap(xmap_in) { }
         void proc(float solvent_content_frac);
         static bool mri_var_pair_sorter(const std::pair<clipper::Xmap_base::Map_reference_index, float> &p1,
                                         const std::pair<clipper::Xmap_base::Map_reference_index, float> &p2);
      };

      // attach the chain-id to each returned map
      //
      // you can set an informative "state" message
      std::vector<std::pair<std::string, clipper::Xmap<float> > >
      partition_map_by_chain(const clipper::Xmap<float> &xmap, mmdb::Manager *mol,
                             std::string *state_string_p);

      bool is_EM_map(const clipper::Xmap<float> &xmap);


      typedef std::pair<double, double> phitheta;

      std::vector<phitheta> make_phi_thetas(unsigned int n_pts);
      float average_of_sample_map_at_sphere_points(clipper::Coord_orth &centre,
                                                   float radius,
                                                   const std::vector<phitheta> &phi_thetas,
                                                   clipper::Xmap<float> &xmap);

      std::vector<std::pair<clipper::Resolution, double> >
      fsc(const clipper::Xmap<float> &xmap_1, const clipper::Xmap<float> &xmap_2);

     // scale map_for_scaling using relion-like resolution binning
     clipper::Xmap<float>
     power_scale(const clipper::Xmap<float> &xmap_ref, const clipper::Xmap<float> &xmap_for_scaling);

      void
      compare_structure_factors(const clipper::Xmap<float> &xmap_1, const clipper::Xmap<float> &xmap_2);

      void flip_hand(clipper::Xmap<float> *xmap_p);

      // input is map and its rmsd
      clipper::Xmap<float>
      analyse_map_point_density_change(const std::vector<std::pair<clipper::Xmap<float> *, float> > &xmaps,
                                       const clipper::Xmap<float> &xmap_for_mask);

      clipper::Xmap<float> zero_dose_extrapolation(const std::vector<std::pair<clipper::Xmap<float> *, float> > &xmaps,
                                                   const clipper::Xmap<float> &xmap_for_mask);

      clipper::Xmap<float> real_space_zero_dose_extrapolation(const std::vector<clipper::Xmap<float> *> &xmaps,
                                                   const clipper::Xmap<float> &xmap_for_mask);

      int split_residue_using_map(mmdb::Residue *residue_p, mmdb::Manager *mol, const clipper::Xmap<float> &xmap);

      std::vector<std::vector<float> >
      get_density_on_cylinder(const clipper::Coord_orth &pt_1, const clipper::Coord_orth &pt_2,
                              const clipper::Coord_orth &pt_ref, const clipper::Xmap<float> &xmap,
                              double radius, unsigned int n_length, unsigned int n_ring);

   }
}

#endif // COOT_MAP_UTILS_HH
