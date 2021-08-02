#ifndef MERGE_C_AND_N_TERMINII_HH
#define MERGE_C_AND_N_TERMINII_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>

namespace coot {

   // not by finding overlapping fragments, try to merge using close N and C terminii and fitting
   // a possible missing residue or two between the N and C terminii and using symmetry
   //
   // use_symmetry, if set, will cause Coot to try to link symmetry-related fragment
   // but this is not coded up yet.
   //
   void merge_C_and_N_terminii(mmdb::Manager *mol,
                               const clipper::Xmap<float> &xmap,
                               bool use_symmetry=true, bool using_missing_loop_fit=true);

   void merge_C_and_N_terminii_0_gap(mmdb::Manager *mol);

}


#endif // MERGE_C_AND_N_TERMINII_HH
