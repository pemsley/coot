
// This file is here because we need to know the rotamer name
//
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"

namespace coot {

   class density_box_t {
   public:
      // density_box is a copy of a filled pointer
      density_box_t(float *density_box, mmdb::Residue *residue_p, int n_steps);
      float *density_box;
      mmdb::Residue *residue_p;
      int n_steps; // either side of the middle
      void scale_by(float scale_factor) {
	 int n = 2 * n_steps + 1;
	 int nnn = n * n * n;
	 for (int i=0; i<nnn; i++)
	    if (density_box[i] > -1000)
	       density_box[i] *= scale_factor;
      }
      int nnn() const {
	 int n = 2 * n_steps + 1;
	 return n * n * n;
      }
      std::pair<float, float> mean_and_variance() const;
      void self_normalize();
      void clear() {
	 delete [] density_box;
      }
      // caller checks for valid index
      float operator[](const unsigned int &idx) const {
	 return density_box[idx];
      }
      bool empty() const { return (n_steps == 0); }
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
      void store_density_box(const density_box_t &density_box) {
	 density_boxes.push_back(density_box);
      }
      std::vector<density_box_t> density_boxes;
      void normalize_density_boxes(const std::string &id);
      void write_density_boxes() const;
      double get_log_likelihood(const unsigned int &grid_idx,
				const density_box_t &block,
				const double &mean,
				const double &variance,
				const double &skew) const;
      std::map<std::string, double>
      compare_block_vs_all_rotamers(density_box_t block,
				    const std::string &data_dir,
				    const clipper::Xmap<float> &xmap);
      // fills the rotamer grid cache, non-const
      std::pair<bool, double>
      compare_block_vs_rotamer(density_box_t block,
			       const std::string &rotamer_dir,
			       const clipper::Xmap<float> &xmap);
      bool in_sphere(int grid_idx, const int &n_steps) const; // manhattan - not a good test
      bool in_sphere(const clipper::Coord_orth &pt, // Euclidean
		     const clipper::Coord_orth &cb,
		     const double &d_max) const;
      std::set<int> useable_grid_points;
      void fill_useable_grid_points_vector(const std::string &file_name);

      clipper::Coord_orth make_pt_in_grid(int ix, int iy, int iz, const float &step_size,
					  const std::vector<clipper::Coord_orth> &axes) const;

      std::map<std::string, double>
      likelihood_of_each_rotamer_at_this_residue(mmdb::Residue *residue_p,
						 const clipper::Xmap<float> &xmap);

      std::string dir_to_key(const std::string &str) const;
      std::pair<std::string, std::string> map_key_to_residue_and_rotamer_names(const std::string &key) const;

      // why is this not in utils or something?
      char single_letter_code(const std::string &res_name) const;

      // cache variable
      std::map<std::string, std::map<unsigned int, std::tuple<double, double, double> > > rotamer_dir_grid_stats_map_cache;

   public:

      std::string id;
      side_chain_densities(const std::string &id_in, mmdb::Manager *mol,
			   int n_steps_in, float grid_box_radius_in,
			   const clipper::Xmap<float> &xmap,
			   const std::string &file_name) {
	 n_steps = n_steps_in;
	 grid_box_radius = grid_box_radius_in;
	 id = id_in;
	 fill_useable_grid_points_vector(file_name);
	 proc_mol(id, mol, xmap);
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
			   const std::string &useable_grid_points_file_name) : n_steps(n_steps_in), grid_box_radius(grid_box_radius_in) { fill_useable_grid_points_vector(useable_grid_points_file_name); }

      void set_data_dir(const std::string &dir) { data_dir = dir; }

      // we want to find the probability distribution from all the sample of that type
      // of rotamer for that particular residue type.
      //
      static void combine_directory(const std::string &dir, int n_steps);

      std::map<std::string, double>
      get_rotamer_likelihoods(mmdb::Residue *residue_p, const clipper::Xmap<float> &xmap);

      void gen_useable_grid_points(mmdb::Residue *residue_this_p,
				   mmdb::Residue *residue_next_p,
				   int n_steps, float grid_box_radius,
				   const std::string &file_name) const;
      void check_stats(mmdb::Residue *residue_p,
		       const std::string &res_name,
		       const std::string &rot_name) const;

      void check_useable_grid_points(mmdb::Residue *residue_p,
				     const std::string &useable_grid_points_mapped_to_residue_file_name) const;
      void test_sequence(mmdb::Manager *mol, const std::string &chain_id, int resno_start, int resno_end,
			 const clipper::Xmap<float> &xmap, const std::string &sequence);


      void probability_of_each_rotamer_at_each_residue(mmdb::Manager *mol,
						       const std::string &chain_id,
						       int resno_start, int resno_end,
						       const clipper::Xmap<float> &xmap);

      // Have a guess at the sequence - choose the best fitting residue at every position
      // and turn that into a string.
      //
      std::string best_guess(mmdb::Manager *mol,
			     const std::string &chain_id,
			     int resno_start, int resno_end,
			     const clipper::Xmap<float> &xmap) const;

   };
}

