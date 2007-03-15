
#include "sec-str-dist-check.hh" // in parallel directory
#include "clipper/core/coords.h"

std::vector<clipper::Coord_orth> z_control_points(int nres);

std::pair<bool, clipper::RTop_orth>
orient_strand_on_z(int handle, CMMDBManager *mol);

// fiddle with mol
void apply_rtop_to_strand(int SelHnd, CMMDBManager *mol,
			  const clipper::RTop_orth &rtop);


void read_dir(const std::string &dir_name);

bool
matches_pdb_name(const std::string &file_str);

std::vector<std::string>
get_reference_pdb_list(const std::string &dir_name);



namespace coot {

   class strand_info_t {
   public: 
      std::string filename;
      CStrand strand;
      int length;

   }; 

   class strands_t {

      std::vector<strand_info_t> strand_infos;
      float residual_between_strands(int istrand, int jstrand, int this_length) const; 
   public:
      void analyse_pdb_file(const std::string &filename);
      void add_strand(const std::string &filename,
		      CMMDBManager *mol,
		      CStrand *strand,
		      int SelectionHandle);
      void strand_analysis(CModel *model_p,
			   CMMDBManager *mol,
			   const std::string &filename);

      void post_read_analysis() const;
   };

}
