
#include "sec-str-dist-check.hh" // in parallel directory
#include "clipper/core/coords.h"

void strand_analysis(CModel *model_p, CMMDBManager *mol,
		     const std::string &filename);

std::vector<clipper::Coord_orth> z_control_points(int nres);

std::pair<bool, clipper::RTop_orth>
orient_strand_on_z(int handle, CMMDBManager *mol);

// fiddle with mol
void apply_rtop_to_strand(int SelHnd, CMMDBManager *mol,
			  const clipper::RTop_orth &rtop);


void analyse_pdb_file(const std::string &filename);
void read_dir(const std::string &dir_name);

bool
matches_pdb_name(const std::string &file_str);

std::vector<std::string>
get_reference_pdb_list(const std::string &dir_name);


