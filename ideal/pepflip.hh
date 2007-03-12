

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include "mmdb_manager.h"
#include "clipper/core/coords.h"

namespace coot {

   // mmdb-style interface.
   //
   // Given a mol and a residue number and chain of the first residue
   // in the peptide (i.e. the residue with the C and O atoms), flip
   // the C and O atoms of this peptide and the N of the next one (if
   // exists) round a line joining this Ca to the next one.
   //
   // Return status is 0 if the flip did not happen (because, for
   // example, either or both of the Ca's could not be found).
   //
   // Typically, one would copy one's mol (and save it) before calling
   // this.
   // 
   int pepflip(CMMDBManager *mol, int resno, 
	       const std::string &altconf,
	       const std::string &chain_id);

   // You are advised against using this externally.
   // 
   std::vector<clipper::Coord_orth> 
   flip_internal(const std::vector<clipper::Coord_orth> &ca,
		 const std::vector<CAtom *> &atoms); 
   

}
