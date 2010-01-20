
#include "graphics-info.h"

#ifdef USE_LIBCURL
// define the static
std::vector<coot::simple_curl_handler_t> graphics_info_t::curl_handlers;
#endif

#ifdef USE_LIBCURL
void
graphics_info_t::add_curl_handle_and_file_name(std::pair<CURL *, std::string> p) {

   coot::simple_curl_handler_t sch(p.first, p.second);
   curl_handlers.push_back(sch);

}
#endif


#ifdef USE_LIBCURL
// remove all handles from curl_handlers that have the filename file_name.
void
graphics_info_t::remove_curl_handle_with_file_name(std::string file_name) {

   bool done = 0;
   std::vector<coot::simple_curl_handler_t>::iterator it;
   while (!done) {
      done = 1;
      for (it=curl_handlers.begin(); it!=curl_handlers.end(); it++) {
	 if (it->file_name == file_name) { 
	    curl_handlers.erase(it);
	    done = 0;
	    break;
	 }
      }
   }

}
#endif


#ifdef USE_LIBCURL
CURL *
graphics_info_t::get_curl_handle_for_file_name(const std::string &filename) const {

   CURL *c = NULL;
   for (unsigned int i=0; i<curl_handlers.size(); i++) {
      if (curl_handlers[i].file_name == filename) {
	 return curl_handlers[i].c;
      }
   }
   return c;
}
#endif

#ifdef USE_LIBCURL
// static
bool
graphics_info_t::curl_handler_stop_it_flag_set(CURL *c) {

   bool r = 0;
   std::vector<coot::simple_curl_handler_t>::iterator it;
   for (it=curl_handlers.begin(); it!=curl_handlers.end(); it++) {
      if (it->stop_is_set()) {
	 r = 1;
	 break;
      } else {
      }
   }
   return r;
} 
#endif

#ifdef USE_LIBCURL
// static 
void
graphics_info_t::set_stop_curl_download_flag(const std::string &file_name) {

   for (unsigned int i=0; i<curl_handlers.size(); i++) {
      if (curl_handlers[i].file_name == file_name) {
	 curl_handlers[i].set_stop();
	 break;
      }
   }
}
#endif
