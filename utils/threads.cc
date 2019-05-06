
#include <stdexcept>
#include <iostream>
#include "coot-utils.hh"

#ifdef  _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

// can return -1 if name is invalid
//
long coot::get_number_of_threads_by_system_call()  {

#ifdef _WIN32
SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  return sysinfo.dwNumberOfProcessors;
#else
  return sysconf(_SC_NPROCESSORS_CONF);
#endif

}

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
	    if (t >= 1024) t = n_threads_default;
	    coot_n_threads = t;
	 }
	 catch (const std::runtime_error &e) {
	    coot_n_threads = 1;
	 }
      } else {
	 // no environment variable, use system call.
	 coot_n_threads = n_threads_default;
	 long n = get_number_of_threads_by_system_call();
	 if (n > 0)
	    coot_n_threads = n;
      }
   }
   return coot_n_threads;
}
