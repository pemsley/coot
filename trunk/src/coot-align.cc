
#include <iostream>
#include <string>
#include <vector>

#include "mmdb_align.h"
#include "mmdb_tables.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "coot-align.hh"

void
coot::chain_mutation_info_container_t::rationalize_insertions() {

   // What is the plan here?
   //
   // We have a list of single insertions.  We want to move the single
   // insertions to a set of insertion ranges.

   if (single_insertions.size() > 0) {

      int min_resno =  9999;
      int max_resno = -9999;
      for (int i=0; i<single_insertions.size(); i++) {
	 // std::cout << "single insertion residue number: "
	 // << single_insertions[i].first.resno << std::endl;
	 if (single_insertions[i].first.resno < min_resno)
	    min_resno = single_insertions[i].first.resno;
	 if (single_insertions[i].first.resno > max_resno)
	    max_resno = single_insertions[i].first.resno;
      }
      //       std::cout << "DEBUG:: rational insertion min and max: "
      // << min_resno << " " << max_resno << std::endl;

      std::vector<std::pair<int, std::string> > ins((max_resno - min_resno + 1), (std::pair<int, std::string>(0, "")));
      for (int i=0; i<single_insertions.size(); i++) {
	 ins[single_insertions[i].first.resno-min_resno].first++;
	 ins[single_insertions[i].first.resno-min_resno].second =
	    single_insertions[i].second;
      }
      short int in_insertion_range = 0;
      int ires_start_of_this_insertion_range = 0;
      int running_end = 0;
      std::vector<std::string> running_types;
      for (int ires=min_resno; ires<=max_resno; ires++) {
	 // std::cout << ires << " " << ins[ires-min_resno].first << std::endl;
	 if (ins[ires-min_resno].first == 1) {
	    if (in_insertion_range == 0) {
	       ires_start_of_this_insertion_range = ires;
	       running_types.resize(0);
	    }
	    in_insertion_range = 1;
	    running_types.push_back(ins[ires-min_resno].second);
	    running_end = ires;
	 } else {
	    // we finish a insertion range maybe:
	    if (in_insertion_range == 1) {
// 	       std::cout << " insertion range from " << ires_start_of_this_insertion_range
// 			 << " to " << running_end << std::endl;
	       insertions.push_back(coot::mutate_insertion_range_info_t(ires_start_of_this_insertion_range, running_types));
	    }
	    in_insertion_range = 0;
	 }
      }
      // were we in an insertion range when we finished?
      if (in_insertion_range) {
// 	 std::cout << " insertion range from " << ires_start_of_this_insertion_range
// 		   << " to " << max_resno << std::endl;
	 insertions.push_back(coot::mutate_insertion_range_info_t(ires_start_of_this_insertion_range, running_types));
      }

      std::cout << "INFO:: There were " << insertions.size() << " insertion ranges\n";
      for (int irange=0; irange<insertions.size(); irange++) {
	 std::cout << "From " << insertions[irange].start_resno << " to "
		   << insertions[irange].end_resno() << std::endl;
	 for (int it=0; it<insertions[irange].types.size(); it++) {
	    std::cout << "    " << insertions[irange].types[it] << std::endl;
	 }
      }
	 
   }
}
