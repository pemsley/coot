
#include <vector>
#include <string>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

namespace coot {

   // return a true in the first if the second is a valid angle (in degress)
   //
   // vector will be typically the difference in positoin between a ring atom and
   // an atom attached to that ring atom.  This function is used to distinugish
   // alpha and beta linking in carbohydrates
   //
   std::pair<bool, double> angle_betwen_plane_and_vector(mmdb::Residue *residue_p,
							 const std::vector<std::string> &ring_atom_names,
							 const std::string &altconf,
							 const clipper::Coord_orth &vector);

   std::pair<bool, double> angle_betwen_plane_and_vector(mmdb::Residue *residue_p,
							 mmdb::Atom *atom_in_ring,
							 mmdb::Atom *bonding_atom);

}
