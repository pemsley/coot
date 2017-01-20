
#include "coot-utils.hh"

// unsigned int coot::coot_n_threads = 0;

unsigned int coot::get_max_number_of_threads() {

   unsigned int n_threads_default = 4;

   // only do this once:
   // (if coot_n_threads is unset, then set it)
   //
   // is this initially set to 0 by the compiler?
   //
   if (coot_n_threads == 0) {
      const char *e = getenv("COOT_N_THREADS");
      if (e) {
	 try {
	    // can throw an runtime_error exception on unable to convert
	    unsigned int t = util::string_to_int(e);
	    if (t > 1024) t = n_threads_default;
	    coot_n_threads = t;
	 }
	 catch (const std::runtime_error &e) {
	    coot_n_threads = n_threads_default;
	 }
      } else {
	 coot_n_threads = n_threads_default;
      }
   }
   return coot_n_threads;
}
