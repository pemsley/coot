/* src/main.cc
 * 
 * Copyright 2009 by The University of Oxford
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <iostream>
#include <stdio.h>

#include <gtk/gtk.h>

#ifdef USE_LIBCURL
#include <curl/curl.h>
#endif

#include "guile-fixups.h"
#include "cc-interface.hh"

// return 0 on success
#ifdef USE_LIBCURL
int coot_get_url(const char *url, const char *file_name) {

   FILE *f = fopen(file_name, "w");
   if (f) { 
      CURL *c = curl_easy_init();
      curl_easy_setopt(c, CURLOPT_URL, url);
      curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_coot_curl_data_to_file);
      curl_easy_setopt(c, CURLOPT_WRITEDATA, f);
      CURLcode success = curl_easy_perform(c);
      curl_easy_cleanup(c);
      fclose(f);
      return success;
   } else {
      return 2; // file not opened.
   } 
}
#endif /* USE_LIBCURL */

#ifdef USE_LIBCURL
std::string coot_get_url_as_string_internal(const char *url) {
   std::string s;
   CURL *c = curl_easy_init();
   curl_easy_setopt(c, CURLOPT_URL, url);
   curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_coot_curl_data);
   curl_easy_setopt(c, CURLOPT_WRITEDATA, &s);
   CURLcode success = curl_easy_perform(c);
   curl_easy_cleanup(c);
   return s;
}
#endif // USE_LIBCURL

#ifdef USE_LIBCURL
#ifdef USE_GUILE
SCM coot_get_url_as_string(const char *url) {
   SCM r = SCM_BOOL_F;
   std::string s = coot_get_url_as_string_internal(url);
   r = scm_from_locale_string(s.c_str());
   return r;
}
#endif /* USE_GUILE */
#endif /* USE_LIBCURL */


#ifdef USE_LIBCURL
size_t
write_coot_curl_data(void *buffer, size_t size, size_t nmemb, void *userp) {

   // std::cout << "size: " << size << " nmeb: " << nmemb;
   if (buffer) {
      char *s = static_cast<char *> (buffer);
      std::string res(s);
      // std::cout << res << std::endl;
      std::string *sp = static_cast<std::string *>(userp);
      *sp += res;
   } else {
      std::cout << std::endl;
   }
   return nmemb; // slightly naughty, we should return the size of the
		 // data that we actually delt with.
}
#endif /* USE_LIBCURL */


#ifdef USE_LIBCURL
size_t
write_coot_curl_data_to_file(void *buffer, size_t size, size_t nmemb, void *userp) {

   if (buffer) {
      FILE *f = static_cast<FILE *> (userp);
      const char *s = static_cast<const char *> (buffer);
      for (size_t i=0; i<nmemb; i++) 
	 fputc(s[i], f);
   }
   return nmemb; // slightly naughty, we should return the size of the
		 // data that we actually delt with.
}
#endif /* USE_LIBCURL */


