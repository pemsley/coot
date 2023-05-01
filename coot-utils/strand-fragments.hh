
#include <vector>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

namespace coot {
	std::pair<std::vector<std::vector<clipper::Coord_orth> >,
                  std::vector<std::vector<clipper::Coord_orth> > >
	   mol_to_5_residue_strand_fragments(mmdb::Manager *mol);
}
