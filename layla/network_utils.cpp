/* layla/network_utils.cpp
 *
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>     // std::rename, std::remove

#include "utils/coot-utils.hh"
#include "network_utils.hpp"

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

std::string
coot::layla::get_drug_via_wikipedia_and_chembl_curl(const std::string &drugname_in) {

   auto write_file = [] (const std::string &s, const std::string &file_name) {
      std::ofstream f(file_name.c_str());
      f << s;
      f.close();
   };

   std::string r;
   std::string drug_name = coot::util::downcase(drugname_in);
   std::string url_pre  = "https://en.wikipedia.org/w/api.php?format=xml&action=query&titles=";
   std::string url_post = "&prop=revisions&rvprop=content";
   std::string url = url_pre + drug_name + url_post;
   std::string result = get_url_as_string(url);

   if (result.length() > 10) {
      std::string::size_type p1 = result.find(">#REDIRECT ");
      if (p1 != std::string::npos) {
         // Redirected. find the new drug name
         std::string s1 = result.substr(p1+13);
         std::string::size_type p2 = s1.find("]]");
         if (p2 != std::string::npos) {
            std::string s2 = s1.substr(0,p2);
            drug_name = s2;
            url = url_pre + drug_name + url_post;
            result = get_url_as_string(url);
         }
      }

      // Use ChEMBL (DrugBank is now blocked by Cloudflare)
      auto get_chembl_code_from_html = [] (const std::string &s) {
         std::string r;
         std::stringstream ss(s);
         std::string line;
         while (std::getline(ss, line)) {
            if (line.find(" ChEMBL ") != std::string::npos) {
               if (line.length() < 80) {
                  std::vector<std::string> parts = coot::util::split_string_no_blanks(line);
                  // format: "| ChEMBL = 25"
                  if (parts.size() == 4) {
                     r = parts[3];
                  }
               }
            }
         }
         return r;
      };

      // extract the molfile string from the ChEMBL JSON response.
      // The molfile is the value of "molfile" key inside "molecule_structures".
      // It is a multi-line string with \n escaped as literal characters in JSON.
      auto extract_molfile_from_chembl_json = [] (const std::string &json_str) {
         std::string r;
         std::string key = "\"molfile\": \"";
         std::string::size_type p1 = json_str.find(key);
         if (p1 != std::string::npos) {
            std::string::size_type start = p1 + key.length();
            // find the closing quote - the molfile value ends with "
            // but contains escaped newlines \\n, so we look for an unescaped quote
            std::string mol;
            bool in_escape = false;
            for (std::string::size_type i=start; i<json_str.length(); i++) {
               char c = json_str[i];
               if (in_escape) {
                  if (c == 'n') mol += '\n';
                  else mol += c;
                  in_escape = false;
               } else {
                  if (c == '\\') {
                     in_escape = true;
                  } else if (c == '"') {
                     break;
                  } else {
                     mol += c;
                  }
               }
            }
            if (mol.length() > 10)
               r = mol;
         }
         return r;
      };

      std::string chembl_code = get_chembl_code_from_html(result);
      if (! chembl_code.empty()) {
         std::cout << "INFO:: ChEMBL code: " << chembl_code << std::endl;
         std::string chembl_url = "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL"
                                + chembl_code + ".json";
         std::string chembl_json = get_url_as_string(chembl_url);
         if (chembl_json.length() > 10) {
            std::string molfile = extract_molfile_from_chembl_json(chembl_json);
            if (! molfile.empty()) {
               std::string mol_file_name = "CHEMBL" + chembl_code + ".mol";
               write_file(molfile, mol_file_name);
               r = mol_file_name;
            } else {
               std::cout << "WARNING:: Failed to extract molfile from ChEMBL JSON" << std::endl;
            }
         } else {
            std::cout << "WARNING:: Failed to fetch ChEMBL data for CHEMBL" << chembl_code << std::endl;
         }
      } else {
         std::cout << "WARNING:: No ChEMBL code found in Wikipedia page for " << drug_name << std::endl;
      }

   }
   return r;
}
