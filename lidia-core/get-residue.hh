
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   // 20180617-PE not yet canonical - might be dangerous to have this here

   mmdb::Residue *get_residue(const residue_spec_t &res_spec, mmdb::Manager *mol);

   clipper::Coord_orth co(mmdb::Atom *at);
}
