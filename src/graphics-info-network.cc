
#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#include "graphics-info.h"

#ifdef USE_LIBCURL
// define the static
std::vector<coot::simple_curl_handler_t> graphics_info_t::curl_handlers;
#endif

#ifdef USE_LIBCURL
bool 
graphics_info_t::add_curl_handle_and_file_name(std::pair<CURL *, std::string> p) {

   while (curl_handlers_lock == 1) {
      coot::usleep(int(100*float(coot::util::random())/float(RAND_MAX)));
   }
   bool done = add_curl_handle_and_file_name_inner(p);
   if (! done)
      add_curl_handle_and_file_name(p); // stack overflow?
   return 1;
}
#endif

#ifdef USE_LIBCURL
bool 
graphics_info_t::add_curl_handle_and_file_name_inner(std::pair<CURL *, std::string> p) {

   bool done=0;
   if (curl_handlers_lock == 0) {
      curl_handlers_lock = 1;
      coot::simple_curl_handler_t sch(p.first, p.second);
      curl_handlers.push_back(sch); // curl_handlers is a static vector
      curl_handlers_lock = 0;
      done = 1;
   }
   return done;
} 
#endif



#ifdef USE_LIBCURL
// remove all handles from curl_handlers that have the filename file_name.
bool
graphics_info_t::remove_curl_handle_with_file_name(std::string file_name) {

   while (curl_handlers_lock == 1) {
      coot::usleep(int(100*float(coot::util::random())/float(RAND_MAX)));
   }
   bool done = remove_curl_handle_with_file_name_inner(file_name);
   if (! done)
      remove_curl_handle_with_file_name(file_name); // stack overflow?
   return 1;
}
#endif

#ifdef USE_LIBCURL
// remove all handles from curl_handlers that have the filename file_name.
bool
graphics_info_t::remove_curl_handle_with_file_name_inner(std::string file_name) {

   bool done = 0;
   if (curl_handlers_lock == 0) {
      curl_handlers_lock = 1;
      bool all_erased = 0;
      std::vector<coot::simple_curl_handler_t>::iterator it;
      while (! all_erased) {
	 all_erased = 1;
	 for (it=curl_handlers.begin(); it!=curl_handlers.end(); it++) {
	    if (it->file_name == file_name) { 
	       curl_handlers.erase(it);
	       all_erased = 0;
	       break;
	    }
	 }
      }
      curl_handlers_lock = 0;
      done = 1;
   }
   return done;
}
#endif


#ifdef USE_LIBCURL
CURL *
graphics_info_t::get_curl_handle_for_file_name(const std::string &filename) const {

   while (curl_handlers_lock == 1) {
      coot::usleep(int(100*float(coot::util::random())/float(RAND_MAX)));
   }
   return get_curl_handle_for_file_name_inner(filename);
}
#endif

#ifdef USE_LIBCURL
CURL *
graphics_info_t::get_curl_handle_for_file_name_inner(const std::string &filename) const {

   CURL *c = NULL;
   if (curl_handlers_lock == 0) {
      curl_handlers_lock = 1;
      for (unsigned int i=0; i<curl_handlers.size(); i++) {
	 if (curl_handlers[i].file_name == filename) {
	    c = curl_handlers[i].c;
	    break;
	 }
      }
      curl_handlers_lock = 0;
   }
   return c;
}
#endif

#ifdef USE_LIBCURL
// static
bool
graphics_info_t::curl_handler_stop_it_flag_set(CURL *c) {

   while (curl_handlers_lock == 1) {
      coot::usleep(int(100*float(coot::util::random())/float(RAND_MAX)));
   }
   return curl_handler_stop_it_flag_set_inner(c);
} 
#endif

#ifdef USE_LIBCURL
// static
bool
graphics_info_t::curl_handler_stop_it_flag_set_inner(CURL *c) {

   // 20100422: check the stop status of the curl handler c, not all of them.

   bool r = 0;
   std::vector<coot::simple_curl_handler_t>::iterator it;
   for (it=curl_handlers.begin(); it!=curl_handlers.end(); it++) {
      if (it->c == c) { 
	 if (it->stop_is_set()) {
	    r = 1;
	    break;
	 }
      }
   }
   // std::cout << "curl_handler_stop_it_flag_set() returns " << r << std::endl;
   return r;
} 
#endif

#ifdef USE_LIBCURL
// static 
void
graphics_info_t::set_stop_curl_download_flag(const std::string &file_name) {

   while (curl_handlers_lock == 1) {
      coot::usleep(int(100*float(coot::util::random())/float(RAND_MAX)));
   }
   set_stop_curl_download_flag_inner(file_name);
}
#endif

#ifdef USE_LIBCURL
// static 
void
graphics_info_t::set_stop_curl_download_flag_inner(const std::string &file_name) {

   for (unsigned int i=0; i<curl_handlers.size(); i++) {
      if (curl_handlers[i].file_name == file_name) {
	 curl_handlers[i].set_stop();
	 break;
      }
   }
}
#endif
