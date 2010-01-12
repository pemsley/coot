
#include "graphics-info.h"

#ifdef USE_LIBCURL
// define the static
std::vector<std::pair<CURL *, std::string> > graphics_info_t::curl_handlers;
#endif

#ifdef USE_LIBCURL
void
graphics_info_t::add_curl_handle_and_file_name(std::pair<CURL *, std::string> p) {

   curl_handlers.push_back(p);

}
#endif


#ifdef USE_LIBCURL
// remove all handles from curl_handlers that have the filename file_name.
void
graphics_info_t::remove_curl_handle_with_file_name(std::string file_name) {

   bool done = 0;
   std::vector<std::pair<CURL *, std::string> >::iterator it;
   while (!done) {
      done = 1;
      for (it=curl_handlers.begin(); it!=curl_handlers.end(); it++) {
	 if (it->second == file_name) { 
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
      if (curl_handlers[i].second == filename) {
	 return curl_handlers[i].first;
      }
   }
   return c;
}
#endif
