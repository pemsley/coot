
#include <atomic>

#include <clipper/core/xmap.h>
#include "mini-mol/mini-mol.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

   // Do I want to make this a class?

   minimol::fragment
   multi_build_N_terminal_ALA(mmdb::Residue *res_p,
			      mmdb::Residue *upstream_neighbour_p,
			      const std::string &chain_id,
			      float b_factor_in,
			      int n_trials,
			      const protein_geometry &geom,
			      const clipper::Xmap<float> &xmap_in,
			      std::pair<float, float> mean_and_variance,
			      bool debug_trials);

   minimol::fragment
   multi_build_C_terminal_ALA(mmdb::Residue *res_p,
			      mmdb::Residue *downstream_neighbour_p,
			      const std::string &chain_id,
			      float b_factor_in,
			      int n_trials,
			      const protein_geometry &geom,
			      const clipper::Xmap<float> &xmap_in,
			      std::pair<float, float> mean_and_variance,
			      bool debug_trials);

   minimol::fragment
   multi_build_terminal_ALA(int offset, // direction
			    mmdb::Residue *res_p,
			    mmdb::Residue *up_or_down_stream_neighbour_p, // depends on offset
			    const std::string &chain_id,
			    float b_factor_in,
			    int n_trials,
			    const protein_geometry &geom,
			    const clipper::Xmap<float> &xmap_in,
			    std::pair<float, float> mean_and_variance,
			    bool debug_trials);

   class stored_fragment_t {
   public:
      class position_triple_t {
	 void fill_residue_atom_positions(const minimol::residue &res);
      public:
	 position_triple_t(const minimol::residue &res) { fill_residue_atom_positions(res); }
	 clipper::Coord_orth positions[3];
      };
   private:
      static mmdb::Residue *get_standard_residue_instance(const std::string &residue_type_in,
							  mmdb::Manager *standard_residues_mol);
      static void apply_sequence(stored_fragment_t &frag_to_be_modified,
				 mmdb::Manager *mol, const std::string &best_seq, int rnoffset,
				 mmdb::Manager *standard_residues_mol,
				 std::atomic<unsigned int> &store_lock);

      static bool matches_position(const position_triple_t &p1,
				   const position_triple_t &p2,
				   const std::vector<clipper::RTop_orth> &symms,
				   double d_crit_sqrd);
      std::vector<std::pair<int, position_triple_t> > residue_atom_positions;
      void fill_residue_atom_positions();
   public:
      stored_fragment_t() { standard_residues_mol = 0; }
      stored_fragment_t(const minimol::fragment &f, int build_dir_in,
			bool with_sidechains_in, mmdb::Manager *m) {
	 frag = f;
	 with_sidechains = with_sidechains_in;
	 sidechains_tried = false;
	 standard_residues_mol = m;
	 fill_residue_atom_positions();
	 build_dir = build_dir_in;
      }
      minimol::fragment frag;
      int build_dir;
      bool with_sidechains;
      bool sidechains_tried;
      static bool try_assign_sidechains(coot::stored_fragment_t &stored_frag,
					std::atomic<unsigned int> &locked,
					const clipper::Xmap<float> &xmap,
					const std::vector<std::pair<std::string, std::string> > &sequences,
					mmdb::Manager *standard_residues_mol);
      mmdb::Manager *standard_residues_mol;
      bool matches_position_in_fragment(const stored_fragment_t::position_triple_t &res_triple_for_testing,
					const std::vector<clipper::RTop_orth> &symms) const;
   };

   class stored_fragment_container_t {
   public:
      stored_fragment_container_t() {
	 all_fragments_stored = false;
      }
      std::vector<stored_fragment_t> stored_fragments;
      bool all_fragments_stored;
      void add(const stored_fragment_t &sf, std::atomic<unsigned int> &store_lock);
      std::size_t size() const { return stored_fragments.size(); }
      std::size_t n_sequenced() const {
	 std::size_t n = 0;
	 for (std::size_t i=0; i<stored_fragments.size(); i++) {
	    if (stored_fragments[i].with_sidechains)
	       n++;
	 }
	 return n;
      }
   };

   class multi_build_terminal_residue_addition {

      void update_O_position_in_prev_residue(mmdb::Residue *res_p,
					     minimol::fragment *many_residues,
					     const minimol::residue &res);

      std::pair<bool, minimol::residue>
      try_to_recover_from_bad_fit_forwards(mmdb::Residue *res_p, mmdb::Residue *res_prev_p,
					   const std::string &chain_id, int n_trials,
					   const protein_geometry &geom,
					   const clipper::Xmap<float> &xmap,
					   std::pair<float, float> mv);

      bool crashing_into_self(const minimol::fragment &many_residues, int seqnum, int offset);

      // make seeds, store them and grow them, sequence them, compare them, consolidate them
      //
      void start_from_map(const protein_geometry &geom);
      const clipper::Xmap<float> &xmap;
      std::pair<float, float> mv; // mean_and_variance

      stored_fragment_container_t fragment_store;
      std::atomic<unsigned int> store_lock;

      void add_to_fragment_store(const minimol::fragment &new_fragment,
				 int build_dir,
				 bool with_sidechains=false);

      // this can't be started with an std::async because it isn't guaranteed to run.
      // So now we start it with a std::thread, and use thread.join().
      //
      static void store_manager(stored_fragment_container_t &fragment_store,
				std::atomic<unsigned int> &store_lock, // as above
				const clipper::Xmap<float> &xmap,
				const std::vector<std::pair<std::string, std::string> > &sequences);

      mmdb::Manager *standard_residues_mol;
      void setup_standard_residues_mol();
      std::vector<clipper::RTop_orth> symms;
      void setup_symms(); // fill symms from xmap
      bool was_this_already_built_p(minimol::residue &res,
				    unsigned int seed_number,
				    int build_dir,
				    std::atomic<unsigned int> &store_lock) const;
      std::vector<std::pair<std::string, std::string> > sequences;
      void init_no_go();
      clipper::Xmap<unsigned int> no_go; // true means "already traced"
      void mask_no_go_map(const minimol::fragment &frag);

      // return true if the atoms of the residue are in the no-go map (i.e. it
      // had been masked because it has been built before).
      //
      bool is_in_no_go_map(minimol::residue &res) const;

   public:

      multi_build_terminal_residue_addition(const protein_geometry &geom,
					    const clipper::Xmap<float> &xmap_in,
					    std::pair<float, float> mv_in,
					    const std::vector<std::pair<std::string, std::string> > &sequences_in) :
	 xmap(xmap_in) {
	 mv = mv_in;
	 setup_standard_residues_mol();
	 setup_symms();
	 sequences = sequences_in;
	 init_no_go();

	 start_from_map(geom);
      }
      minimol::fragment
      forwards_2018(unsigned int iseed,
		    mmdb::Residue *res_p,
		    mmdb::Residue *upstream_neighbour_p,
		    const std::string &chain_id,
		    float b_factor_in,
		    int n_trials,
		    const protein_geometry &geom,
		    const clipper::Xmap<float> &xmap,
		    std::pair<float, float> mv,
		    bool debug_trials);

      minimol::fragment
      backwards_2018(unsigned int iseed,
		     mmdb::Residue *res_p,
		     mmdb::Residue *upstream_neighbour_p,
		     const std::string &chain_id,
		     float b_factor_in,
		     int n_trials,
		     const protein_geometry &geom,
		     const clipper::Xmap<float> &xmap,
		     std::pair<float, float> mv,
		     bool debug_trials);

      // because kludge, this needs to be public
      void refine_end(coot::minimol::fragment *many_residues,
		      int seqnum, int offset,
		      const protein_geometry &geom,
		      const clipper::Xmap<float> &xmap_in);

      // because kludge, this needs to be public
      bool does_residue_fit(const coot::minimol::residue &res,
			    const clipper::Xmap<float> &xmap,
			    std::pair<float, float> mv);

   };

}
