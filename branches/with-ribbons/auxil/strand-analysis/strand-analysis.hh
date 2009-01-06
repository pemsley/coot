
#include "sec-str-dist-check.hh" // in parallel directory
#include "clipper/core/coords.h"

std::vector<clipper::Coord_orth> z_control_points(int nres);

std::pair<bool, clipper::RTop_orth>
orient_strand_on_z(int handle, CMMDBManager *mol);

// fiddle with mol
void apply_rtop_to_strand(int SelHnd, CMMDBManager *mol,
			  const clipper::RTop_orth &rtop);


void read_dir(const std::string &dir_name, int strand_length);

bool
matches_pdb_name(const std::string &file_str);

std::vector<std::string>
get_reference_pdb_list(const std::string &dir_name);



namespace coot {

   class stats_t {
   public:
      float mean;
      float var;
      stats_t(float mean_in, float var_in) {
	 mean = mean_in;
	 var = var_in;
      } 
   };

   stats_t get_rtop_and_apply(const std::vector<clipper::Coord_orth> &found_atoms_strand_1,
			      const std::vector<clipper::Coord_orth> &found_atoms_strand_2,
			      PPCResidue SelResidues2, int nSelResidues2);

   class strand_info_t {
   public: 
      std::string name;
      CStrand strand;
      int length;
      CMMDBManager *mol;
      int SelectionHandle;
   }; 

   class strands_t {

      std::vector<strand_info_t> strand_infos;
      float residual_between_strands(int istrand, int jstrand, int this_length) const; 
      std::vector<std::pair<std::string,std::vector<float> > >
      filter_bad_fits(const std::vector<std::pair<std::string,std::vector<float> > > &dist_arr) const;
      std::vector<std::pair<float, float> >
      get_stats(const std::vector<std::pair<std::string, std::vector<float> > > &dist_arr) const;
      std::vector<std::pair<std::string, std::vector<float> > >
      clear_out_column(int i, const std::vector<std::pair<std::string, std::vector<float> > > &d) const;

   public:
      void analyse_pdb_file(const std::string &filename);
      void add_strand(const std::string &filename,
		      CMMDBManager *mol,
		      CStrand *strand,
		      int SelectionHandle);
      void strand_analysis(CModel *model_p,
			   CMMDBManager *mol,
			   const std::string &filename);

      void post_read_analysis(int this_length) const;
   };

}
