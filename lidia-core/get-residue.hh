
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   namespace lidia_utils {
      mmdb::Residue *get_residue(const residue_spec_t &res_spec, mmdb::Manager *mol);

      clipper::Coord_orth co(mmdb::Atom *at);

   }
}
