#ifndef SIDE_CHAIN_DENSITIES_HH
#define SIDE_CHAIN_DENSITIES_HH


#include <atomic>
// This file is here because we need to know the rotamer name
//
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/fragment-container.hh"
#include "utils/coot-fasta.hh"

namespace coot {

   class density_box_t {
   public:
      // density_box is a copy of a filled pointer
      density_box_t(float *density_box, mmdb::Residue *residue_p, int n_steps);
      density_box_t() { init(); } // needed because it's unsed in a map
      float *density_box;
      mmdb::Residue *residue_p;
      double mean;
      double mean_of_positives;
      double var;

      double mean_around_ca;
      double mean_of_positives_around_ca;
      double var_around_ca;

      bool is_weird;

      int n_steps; // either side of the middle

      void init() {
         density_box = 0; residue_p = 0; n_steps = 0; mean=0; var = -1;
         mean_around_ca = 0; mean_of_positives_around_ca = 0;
         var_around_ca = -1;
         is_weird = false;
         mean_of_positives = 0;
      }
      void scale_by(float scale_factor) {
         if (n_steps > 0) {
            int n = 2 * n_steps + 1;
            int n3 = n * n * n;
            for (int i=0; i<n3; i++)
               if (density_box[i] > -1000)
                  density_box[i] *= scale_factor;
         }
      }
      int nnn() const {
         int n = 2 * n_steps + 1;
         return n * n * n;
      }
      std::pair<float, float> mean_and_variance() const;
      void self_normalize();
      void clear() {
         delete [] density_box;
         density_box = 0;
      }
      // caller checks for valid index
      float operator[](const unsigned int &idx) const {
         return density_box[idx];
      }
      bool empty() const { return (n_steps == 0); }
      void set_stats(const double &mean_in, const double &var_in, const double &mean_of_positives_in) {
         mean = mean_in;
         mean_of_positives = mean_of_positives_in;
         var  = var_in;
      }
      void set_around_ca_stats(const double &mean_in, const double &var_in, const double &mean_of_positives_in) {
         mean_around_ca = mean_in;
         var_around_ca = var_in;
         mean_of_positives_around_ca = mean_of_positives_in;
      }

      void normalize_using_ca_stats();
   };


   class side_chain_densities {

      // SAMPLE_FOR_RESIDUE means that we are generating a box to test a map
      // where we don't know the residue type
      //
      enum mode_t { GEN_USABLE_POINTS, SAMPLE_FOR_DB, SAMPLE_FOR_RESIDUE };

      int n_steps;
      float grid_box_radius;
      std::string data_dir; // we can we find side-chain-data?
      std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > get_residue_axes(mmdb::Residue *res) const;
      std::vector<clipper::Coord_orth>
      make_axes(const clipper::Coord_orth &pt_ca_this,
                const clipper::Coord_orth &pt_cb_this,
                const clipper::Coord_orth &pt_c_this,
                const clipper::Coord_orth &pt_n_this) const;
      std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> >
      get_residue_axes_type_GLY(mmdb::Residue *this_residue) const;
      std::tuple<int, int, int> grid_index_to_grid(int idx) const;

      void proc_mol(const std::string &id, mmdb::Manager *mol, const clipper::Xmap<float> &xmap);
      void proc_chain(const std::string &id, mmdb::Chain *chain, const clipper::Xmap<float> &xmap);
      void proc_residue(mmdb::Residue *residue, const clipper::Xmap<float> &xmap);
      // residue_next_p can be null - residue_next_p needs to be set when
      // generating usable grid points
      // User should check that the densty box is ok using empty()
      density_box_t sample_map(mmdb::Residue *residue_this_p,
                               mmdb::Residue *residue_next_p,
                               mode_t mode,
                               const clipper::Coord_orth &cb_pt,
                               const std::vector<clipper::Coord_orth> &axes,
                               const clipper::Xmap<float> &xmap,
                               std::string gen_pts_file_name = "") const;
      std::string get_rotamer_name(mmdb::Residue *r) const;
      bool is_close_to_atoms(const std::vector<std::pair<double, clipper::Coord_orth> > &atom_positions,
                             const clipper::Coord_orth &test_position) const;
      void write_density_box(float *density_box, int n_steps, const std::string &id,
                             mmdb::Residue * residue_p) const;
      void write_density_box(const density_box_t &db, const std::string &id) const;

      void store_density_box(const density_box_t &density_box) {
         density_boxes.push_back(density_box);
      }
      std::vector<density_box_t> density_boxes;
      void normalize_density_boxes(const std::string &id);
      void normalize_density_boxes_v1(const std::string &id);
      void normalize_density_boxes_v2(const std::string &id);
      void normalize_density_boxes_v3(const std::string &id);
      void add_mean_and_variance_to_individual_density_blocks();
      void write_density_boxes() const;
      double get_log_likelihood(const unsigned int &grid_idx,
                                const density_box_t &block,
                                const double &mean,
                                const double &variance,
                                const double &skew) const;
      double get_log_likelihood_ratio(const unsigned int &grid_idx,
                                      const density_box_t &block,
                                      const std::string &rotamer_dir, // for debugging
                                      const double &step_size,
                                      const double &mean,
                                      const double &variance,
                                      const double &skew) const;
      bool get_test_map_is_above_model_mean(const unsigned int &grid_idx,
                                            const density_box_t &block,
                                            const double &mean) const;
      std::map<std::string, double>
      compare_block_vs_all_rotamers(density_box_t block,
                                    mmdb::Residue *residue_p, // for debugging
                                    const std::string &data_dir,
                                    const std::pair<bool, std::vector<std::pair<std::string, std::string> > > &rotamer_limits,
                                    const clipper::Xmap<float> &xmap);
      // fills the rotamer grid cache, non-const
      std::pair<bool, double>
      compare_block_vs_rotamer(density_box_t block,
                               mmdb::Residue *residue_p,
                               const std::string &rotamer_dir,
                               const clipper::Xmap<float> &xmap);
      bool in_sphere(int grid_idx, const int &n_steps) const; // manhattan - not a good test
      bool in_sphere(const clipper::Coord_orth &pt, // Euclidean
                     const clipper::Coord_orth &cb,
                     const double &d_max) const;
      std::set<int> useable_grid_points;
      void fill_useable_grid_points_vector(const std::string &file_name);
      double get_grid_point_distance_from_grid_centre(const unsigned int &idx,
                                                      const double &step_size) const;

      clipper::Coord_orth make_pt_in_grid(int ix, int iy, int iz, const float &step_size,
                                          const std::vector<clipper::Coord_orth> &axes) const;

      std::string dir_to_key(const std::string &str) const;
      std::pair<std::string, std::string> map_key_to_residue_and_rotamer_names(const std::string &key) const;

      double null_hypothesis_scale;
      double null_hypothesis_sigma;

      // cache variable
      std::map<std::string, std::map<unsigned int, std::tuple<double, double, double> > > rotamer_dir_grid_stats_map_cache;
      // cache variable (for better normalization of user/test map/model)
      std::map<mmdb::Residue *, density_box_t> density_block_map_cache;
      // called by fill_residue_blocks()
      void normalize_density_blocks();
      // use the above cache
      density_box_t get_block(mmdb::Residue *residue_p) const;
      // and the cache of likelihoods for the best rotamer of each residue type at each position:
      // this gets added to using the results lock
      std::map<mmdb::Residue *, std::map<std::string, double> > best_score_for_res_type_cache;

      std::map<int, std::string> make_sequence_for_chain(mmdb::Chain *chain_p) const;

      bool like_the_others(const std::map<int, std::string> &chain,
                           const std::vector<std::map<int, std::string> > &other_chains) const;

      // return negative values on failure
      std::tuple<double, double, double>
      get_stats_around_ca(mmdb::Residue *residue_this,
                          const std::vector<clipper::Coord_orth> &axes,
                          float step_size,
                          const clipper::Xmap<float> &xmap) const;

      std::map<std::string, double> relabun;
      double get_relabun(const std::string &res_name);

   public:

      std::string id;
      std::atomic<bool> results_addition_lock;
      side_chain_densities(const std::string &id_in, mmdb::Manager *mol,
                           int n_steps_in, float grid_box_radius_in,
                           const clipper::Xmap<float> &xmap,
                           const std::string &file_name) : id(id_in) {
         n_steps = n_steps_in;
         grid_box_radius = grid_box_radius_in;
         fill_useable_grid_points_vector(file_name);
         proc_mol(id, mol, xmap);
         null_hypothesis_scale = 1.0;
         null_hypothesis_sigma = 1.0;
         set_default_magic_numbers(); // probably not needed
         results_addition_lock = false;
      }

      // constructor for testing residues vs the "database" - using
      // get_rotamer_likelihoods()
      //
      // also constructor for generating the useable grid points file(s)
      //
      // maybe this should be a different class?
      //
      side_chain_densities(int n_steps_in,
                           float grid_box_radius_in,
                           const std::string &useable_grid_points_file_name) {
         init(n_steps_in, grid_box_radius_in, useable_grid_points_file_name);
      }

      // "wrapper" for above using default values
      side_chain_densities();

      void init(int n_steps_in, float grid_box_radius_in, const std::string &useable_grid_points_file_name) {
         n_steps = n_steps_in;
         grid_box_radius = grid_box_radius_in;
         fill_useable_grid_points_vector(useable_grid_points_file_name);
         null_hypothesis_scale = 1.0;
         null_hypothesis_sigma = 1.0;
         set_default_magic_numbers();
         results_addition_lock = false;
      }

      // magic numbers
      double mn_log_likelihood_ratio_difference_min;
      double mn_scale_for_normalized_density;
      double mn_density_block_sample_x_max;

      class results_t {
      public:
         int offset;
         double sum_score;
         unsigned int n_scored_residues;
         std::string sequence;
         std::vector<std::pair<std::string, std::string> > sequence_residue_type_and_rotamer_name;
         std::string sequence_name;
         std::string true_sequence; // for testing/analysis
         results_t() {
            offset = -1; // unset
            n_scored_residues = 0;
            sum_score = 0;
         }
         results_t(const int &offset_in, const float &f, const unsigned int &n_scored_residues_in,
                   const std::string &running_sequence_in,
                   const std::string &gene_name_in,
                   const std::string &true_sequence_in) : offset(offset_in),
                                                          sum_score(f),
                                                          n_scored_residues(n_scored_residues_in),
                                                          sequence(running_sequence_in),
                                                          sequence_name(gene_name_in),
                                                          true_sequence(true_sequence_in) {}
      };
      std::map<std::string, std::vector<results_t> > results_container;
      void get_results_addition_lock();
      void release_results_addition_lock();

      void set_default_magic_numbers() {
         // magic numbers
         mn_log_likelihood_ratio_difference_min = -18.0;
         mn_scale_for_normalized_density = 1.0;
         mn_density_block_sample_x_max = 13.0;
      }

      std::vector<mmdb::Residue *> make_a_run_of_residues(mmdb::Manager *mol, const std::string &chain_id,
                                                          int resno_start, int resno_end) const;

      void set_magic_number(const std::string &mn_name, double val);

      // Confusing output if this is not called after constructor. So maybe
      // the data_dir should be an argument to the constructor.
      void set_data_dir(const std::string &dir) { data_dir = dir; }

      std::map<std::string, double>
      likelihood_of_each_rotamer_at_this_residue(mmdb::Residue *residue_p,
                                                 const clipper::Xmap<float> &xmap,
                                                 bool limit_to_correct_rotamers_only=false,
                                                 bool verbose_output_mode = false);

      void set_null_hypothesis_scale_and_sigma(double scale, double sigma) {
         null_hypothesis_scale = scale;
         null_hypothesis_sigma = sigma;
      }

      // a function to density block map cache
      void fill_residue_blocks(const std::vector<mmdb::Residue *> &residues,
                               const clipper::Xmap<float> &xmap);
      // above is called by
      void fill_residue_blocks(mmdb::Manager *mol, const std::string &chain_id,
                               int resno_start, int resno_end,
                               const clipper::Xmap<float> &xmap);

      // we want to find the probability distribution from all the sample of that type
      // of rotamer for that particular residue type.
      //
      static void combine_directory(const std::string &dir, int n_steps,
                                    double mn_unreliable_minimum_counts,
                                    double mn_unreliable_minimum_counts_for_low_variance,
                                    double mn_unreliable_minimum_variance,
                                    double mn_use_this_variance_for_unreliable);

      std::map<std::string, double>
      get_rotamer_likelihoods(mmdb::Residue *residue_p,
                              const clipper::Xmap<float> &xmap,
                              bool limit_to_correct_rotamers_only = false,
                              bool verbose_output_mode = true);

      void gen_useable_grid_points(mmdb::Residue *residue_this_p,
                                   mmdb::Residue *residue_next_p,
                                   int n_steps, float grid_box_radius,
                                   const std::string &file_name) const;
      void check_stats(mmdb::Residue *residue_p,
                       const std::string &res_name,
                       const std::string &rot_name) const;

      void check_useable_grid_points(mmdb::Residue *residue_p,
                                     const std::string &useable_grid_points_mapped_to_residue_file_name) const;

      // return an error message (if any of the residues didn't have a mainchain and CB) or a vector of residues
      //
      std::pair<std::string, std::vector<mmdb::Residue *> >
      setup_test_sequence(mmdb::Manager *mol, const std::string &chain_id, int resno_start, int resno_end,
                               const clipper::Xmap<float> &xmap);

      void test_sequence(const std::vector<mmdb::Residue *> &a_run_of_residues,
                         const clipper::Xmap<float> &xmap,
                         const std::string &sequence_name,    // from fasta file
                         const std::string &sequence);

      void setup_likelihood_of_each_rotamer_at_every_residue(const std::vector<mmdb::Residue *> &a_run_of_residues,
                                                             const clipper::Xmap<float> &xmap);

      // find the best result stored by the above function.
      results_t get_result(bool only_return_result_if_probably_correct, bool print_sequencing_solutions_flag=false) const;

      // return the "guessed" sequence
      std::string
      guess_the_sequence(mmdb::Manager *mol,
                         const std::string &chain_id,
                         int resno_start, int resno_end,
                         const clipper::Xmap<float> &xmap,
                         bool verbose_output_mode = false);

      // Have a guess at the sequence - choose the best fitting residue at every position
      // and turn that into a string.
      //
      std::string best_guess(mmdb::Manager *mol,
                             const std::string &chain_id,
                             int resno_start, int resno_end,
                             const clipper::Xmap<float> &xmap) const;

      bool test_grid_point_to_coords_interconversion() const;
   };

   std::vector<std::pair<fragment_container_t::fragment_range_t, std::vector<side_chain_densities::results_t> > >
   get_fragment_sequence_scores(mmdb::Manager *mol,
                                const fasta_multi &fam,
                                const clipper::Xmap<float> &xmap);

}



#endif // SIDE_CHAIN_DENSITIES_HH
