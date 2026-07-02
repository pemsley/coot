/* layla/url_fetch.cpp
 *
 * Copyright 2026 by Medical Research Council
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

#include <fstream>
#include <iostream>
#include <cstdio>     // std::rename, std::remove

#include "url_fetch.hpp"

#ifdef USE_LIBCURL
#ifndef HAVE_CURL_H
#define HAVE_CURL_H
#endif
#include <curl/curl.h>
#endif // USE_LIBCURL

#ifndef VERSION
#define VERSION "unknown"
#endif

namespace {

#ifdef USE_LIBCURL
   // Accumulate the response into the std::string passed as userp.
   size_t write_to_string(void *buffer, size_t size, size_t nmemb, void *userp) {
      if (buffer && userp) {
         std::string *sp = static_cast<std::string *>(userp);
         sp->append(static_cast<char *>(buffer), size * nmemb);
      }
      return size * nmemb;
   }

   std::string user_agent_string() {
      std::string ua = "coot ";
      ua += VERSION;
      ua += " https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/";
      return ua;
   }
#endif // USE_LIBCURL
}

std::string
coot::layla::get_url_as_string(const std::string &url) {

   std::string s;

#ifdef USE_LIBCURL
   char errbuf[CURL_ERROR_SIZE];
   errbuf[0] = '\0';
   std::string ua = user_agent_string();

   CURL *c = curl_easy_init();
   if (! c) return s;
   curl_easy_setopt(c, CURLOPT_URL, url.c_str());
   curl_easy_setopt(c, CURLOPT_NOSIGNAL, 1L);
   curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 6L);
   curl_easy_setopt(c, CURLOPT_TIMEOUT, 30L);
   curl_easy_setopt(c, CURLOPT_SSL_VERIFYPEER, 0L);
   curl_easy_setopt(c, CURLOPT_USERAGENT, ua.c_str());
   curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);
   curl_easy_setopt(c, CURLOPT_ERRORBUFFER, errbuf);
   curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_to_string);
   curl_easy_setopt(c, CURLOPT_WRITEDATA, &s);

   CURLcode cc = curl_easy_perform(c);
   long http_code = 0;
   curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &http_code);
   curl_easy_cleanup(c);

   if (cc != CURLE_OK) {
      std::cout << "WARNING:: layla::get_url_as_string() failed for " << url
                << ": " << errbuf << std::endl;
      s.clear();
   } else if (http_code >= 400) {
      std::cout << "WARNING:: layla::get_url_as_string() HTTP " << http_code
                << " for " << url << std::endl;
      s.clear();
   }
#else
   std::cout << "WARNING:: layla::get_url_as_string() - no libcurl support" << std::endl;
#endif // USE_LIBCURL

   return s;
}

int
coot::layla::get_url(const std::string &url, const std::string &output_file_name,
                     long min_bytes) {

   std::string s = get_url_as_string(url);
   if (static_cast<long>(s.size()) < min_bytes)
      return 1; // failure (also covers empty result from a failed fetch)

   std::string tmp_file_name = output_file_name + ".tmp";
   {
      std::ofstream f(tmp_file_name.c_str(), std::ios::binary);
      if (! f) {
         std::cout << "WARNING:: layla::get_url() cannot write " << tmp_file_name << std::endl;
         return 1;
      }
      f.write(s.data(), s.size());
      f.close();
      if (! f) {
         std::cout << "WARNING:: layla::get_url() write error for " << tmp_file_name << std::endl;
         std::remove(tmp_file_name.c_str());
         return 1;
      }
   }

   if (std::rename(tmp_file_name.c_str(), output_file_name.c_str()) != 0) {
      std::cout << "WARNING:: layla::get_url() cannot rename " << tmp_file_name
                << " to " << output_file_name << std::endl;
      std::remove(tmp_file_name.c_str());
      return 1;
   }
   return 0;
}
