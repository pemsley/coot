
#include <vector>
#include <string>

// clipper coords
#include "clipper/core/coords.h"
#include "mini-mol.hh"

// for directory reading
#include <sys/types.h>   
#include <dirent.h>     
// #include <unistd.h> not needed?

// #include <string.h> // for string matching // needed?

namespace coot { 

   class coord_fragment_wrapper_t {
      
   public:
      std::vector<clipper::Coord_orth> fragment; 
      std::string segment_id; 
      int size() { return fragment.size(); }
   };

   class main_fragment_t {
      int i_start_res_;
   public:
      int ilength; 
      int molecule_number; 
      std::vector<float> sqrt_eigen_values; 
      std::string segment_id; 
      main_fragment_t(int i_start_res, int mol_no,
		      std::vector<float> eigns_sqrt, std::string seg_id,
		      int ilen); 
      int i_start_res() const {return i_start_res_;}
   };

   class db_fitting_result {
   public:
      clipper::RTop_orth rtop;
      int db_frag_index;
      float deviance;
      int ilength;
      int istart_res_of_ca_target;
      db_fitting_result(const clipper::RTop_orth rtop_in, int db_frag_i, float devi, int ist_res, int ilen) { 
	 rtop = rtop_in; 
	 db_frag_index = db_frag_i; 
	 deviance = devi;
	 ilength = ilen;
	 istart_res_of_ca_target = ist_res; 
      }
   };

   class weighted_residue : public minimol::residue {
   public:
      float weight_sum;
      float weight_sum_cb;
      short int have_cb_flag;
      int cb_index;

      weighted_residue() {
	 have_cb_flag = 0;
 	 weight_sum = 0.0;
 	 weight_sum_cb = 0.0;
      }

      // Adding residues with the weight_target_devi.
      // 
      // Given a residue with the weight_target_devi.  We apply
      // weight_pos_in_frag... 
      //
      // Note that we implicitly presume that all in_res has CA, C, N and O,
      // but only some have CB.  So we need to keep a separate track of the
      // sum of the CB weights (which is different from the sum of all
      // weights).
      //
      // Notice that when we start, the vector atoms is empty, so we do a
      // check that only fails for the first time, where we do an addatom,
      // rather than adding to the position of already esisting atoms.
      //
      // For the first time, we check that CA, C, N, O exist after we have
      // added in_res.
      // 
      void add_residue_pos(const minimol::residue &in_res,
			   const clipper::RTop_orth &rtop,
			   float weight);

      minimol::residue pull_residue() const; 
   };

   class db_main {
      std::vector<coot::main_fragment_t> mainchain_frag_db;
      std::vector<coot::minimol::molecule> molecule_list; 

      // For each set of 6 residues (Ca's thereof), we fit a set of
      // "similar by eigen" structures.  This is a vector of the sets
      // of (typically) 6 atoms (typically stepped by 3).
      //
      // We store also the starting residue number (of the target) of
      // the match in the db_fitting_result.
      //
      // Match target fragment fills this vector.
      // 
      std::vector <std::vector <db_fitting_result> > big_results;
      float max_devi;
      std::vector<weighted_residue> output_fragment;
      // Save these in match_target_fragment(), because they will be
      // needed in merge_fragments().
      int iresno_start;
      int iresno_end;
      std::string target_fragment_fragment_id; 

      // Given a molecule and a starting residue number, return a
      // vector of length ilength that contains the Ca positions of
      // that fragment.  We use (only) the first fragment (0) in
      // target.
      // 
      // You need to check the length of the returned vector to see if
      // this function was successful.
      // 
      std::vector<clipper::Coord_orth>
      get_target_ca_coords(int iresno, int ilength,
			   const coot::minimol::molecule &target) const;
      std::vector<db_fitting_result>
      fit_reference_structures(float max_devi,
			       std::vector<clipper::Coord_orth> target_ca,
			       int iresno_start,
			       int ilength) const; 

      short int similar_eigens(float tolerance_frac, 
			       const std::vector<float> &target, 
			       const std::vector<float> &frag) const; 

      std::vector<std::string> get_reference_pdb_list() const;
      
      coot::minimol::fragment get_fragment_ca_atoms(int istart, int ilength,
					   const coot::minimol::molecule &m) const;
      short int matches_pdb_name(std::string file_str) const;
      clipper::Matrix<float> make_cov_matrix(const std::vector<clipper::Coord_orth> &c) const;
      std::vector<clipper::Coord_orth>
      frag_to_coords(const coot::minimol::fragment &fragment) const;

      // This is an awkward function. We have to generate the fragment
      // from the molecule and then select the (ilength) Ca from the
      // right place (and segment id).
      // 
      std::vector<clipper::Coord_orth>
      mainchain_ca_coords_of_db_frag(int i, int ilength) const;
      float deviance(const std::vector<clipper::Coord_orth> &frag1, 
		     const std::vector<clipper::Coord_orth> &frag2, 
		     const clipper::RTop_orth &rtop) const; 

      // helper function for merge_fragments:
      float weight_pos_in_frag(int ipos, int ilength) const;
      // 
      coot::minimol::residue pull_db_residue(const coot::db_fitting_result &res,
					     int ipos) const; 

   public:
      db_main() { max_devi = 2.0; }

      // Fill the molecule_list with (minimol) molecules and
      // mainchain_frag_db with fragment info.
      // 
      int fill_with_fragments(int ilength);

      // So, we have a molecule, we want to chop it into fragments,
      // starting at residue istart.  We don't want to keep the actual
      // fragments, we will just store just the eigen values of the
      // fragment that we pull out (and the discard), along with
      // enough other information to get the fragment back if/when we
      // want to use it in recombination.
      //
      // 
      // Fill big results:
      // 
      void match_target_fragment(const coot::minimol::molecule &target_cas,
				 int iresno_target_start, // typically 1
				 int iresno_target_end, 
				 int ilength);

      // Molecular Graphics usage: We only want to load the graphics
      // once.  So let's do that and store it in a static.
      // 
      short int is_empty() const; 

      // This is the big complex one:
      //
      void merge_fragments();

      // And we return the results in a minimol::fragment:
      // 
      minimol::fragment mainchain_fragment() const;

      // clear results
      void clear_results(); 
   };

} // namespace coot
