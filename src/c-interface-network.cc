/* src/c-interface-network.cc
 * 
 * Copyright 2009, 2012 by The University of Oxford
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

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#include <iostream>
#include <stdexcept>

#include <stdio.h>

#include <gtk/gtk.h>

// this is defined in graphics-info.h, so omit here?!
//#ifdef USE_LIBCURL
//#ifndef HAVE_CURL_H
//#define HAVE_CURL_H
//#include <curl/curl.h>
//#endif // HAVE_CURL_H
//#endif

#include "guile-fixups.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"
#include "cc-interface-network.hh"

// #include "c-interface.h" // for coot_version().  Do we really need all of that..?

#include "graphics-info.h" // because that is where the curl handlers and filenames vector is stored

// return 0 on success
#ifdef USE_LIBCURL
int coot_get_url(const char *url, const char *file_name) {
   return coot_get_url_and_activate_curl_hook(url, file_name, 0);
}
#endif /* USE_LIBCURL */


#ifdef USE_LIBCURL
int coot_get_url_and_activate_curl_hook(const char *url, const char *file_name,
					short int activate_curl_hook_flag) {

   // This can take a while to download files - i.e. can block.  If
   // this is called in a guile thread, then that is bad because it
   // stops the syncing of guile threads (1.8 behavour). See "Blocking
   // in Guile Mode" in the Guile Reference Manual.  Apparently 1.9
   // does not have this issue.  Anyway, with 1.8 we should not run a
   // long lived thread in guile-mode, i.e. we need to exit guile-mode
   // and that is done with scm_without_guile(), and we pass the
   // long-lived function to that.
   //
   // Thanks to Andy Wingo for help here.

   // write binary
   FILE *f = fopen(file_name, "wb");

   if (f) {

      // If you see a crash in curl_easy_init() then the problem is
      // somewhere else.  Curl_open() doesn't do anything except a bit
      // of mallocing.  So the memory is messed up elsewhere and beforehand.
      CURL *c = curl_easy_init();
      long int no_signal = 1;

      std::pair<FILE *, CURL *> p_for_write(f,c);
      curl_easy_setopt(c, CURLOPT_URL, url);
      curl_easy_setopt(c, CURLOPT_NOSIGNAL, no_signal);
      curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 6);
      curl_easy_setopt(c, CURLOPT_SSL_VERIFYPEER, FALSE);
      std::string user_agent_str = "Coot-";
      user_agent_str += VERSION;
      user_agent_str += " https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/";
      curl_easy_setopt(c, CURLOPT_USERAGENT, user_agent_str.c_str());
      curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_coot_curl_data_to_file);
      curl_easy_setopt(c, CURLOPT_WRITEDATA, &p_for_write);
      std::pair <CURL *, std::string> p(c,file_name);
      CURLcode success = CURLcode(-1);
      if (activate_curl_hook_flag) {
	 graphics_info_t g;
	 g.add_curl_handle_and_file_name(p);
#ifdef USE_GUILE
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
	 // good return values (same as the non-wrapped function call)
	 success = CURLcode(GPOINTER_TO_INT(scm_without_guile(wrapped_curl_easy_perform, c)));
#else
	 std::cout << "Can't do this with this old guile" << std::endl;
#endif
#else
#ifdef USE_PYTHON
         Py_BEGIN_ALLOW_THREADS;
	 success = curl_easy_perform(c);
         Py_END_ALLOW_THREADS;
#else
	 success = curl_easy_perform(c);
#endif /* USE_PYTHON */
#endif /* USE_GUILE	*/
	 g.remove_curl_handle_with_file_name(file_name);
      } else {
	 success = curl_easy_perform(c);
         // std::cout << "coot_get_url_and_activate_curl_hook() here with status " << success << std::endl;
      }

      fclose(f);
      curl_easy_cleanup(c);
      return success;
   } else {
      return 2; // file not opened.
   }
}
#endif /* USE_LIBCURL */

#ifdef USE_LIBCURL
void *wrapped_curl_easy_perform(void *data) {

   CURL *c = (CURL *) data;
   int i = curl_easy_perform(c);
   return GINT_TO_POINTER(i);
}
#endif

#ifdef USE_LIBCURL
std::string coot_get_url_as_string_internal(const char *url) {
   std::string s;

   // Wikipedia wants a user-agent
   // 
   std::string user_agent = PACKAGE;
   user_agent += " ";
   user_agent += VERSION;
   user_agent += " http://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/";

   long int no_signal = 1; 
   CURL *c = curl_easy_init();
   curl_easy_setopt(c, CURLOPT_URL, url);
   curl_easy_setopt(c, CURLOPT_NOSIGNAL, no_signal);
   curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 10);
   curl_easy_setopt(c, CURLOPT_USERAGENT, user_agent.c_str());
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

#ifdef USE_PYTHON
PyObject *coot_get_url_as_string_py(const char *url) {
   PyObject *r  = Py_False;
   std::string s = coot_get_url_as_string_internal(url);
   r = PyString_FromString(s.c_str());
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif /* USE_PYTHON */
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

   // Wow, userp is complicated.  We need to know if the Cancel
   // button of the download GUI has been pressed.  Doing so sets the
   // stop boolean for the approriate curl_handler.
   // 
   // Here we lookup the curl handler and see if the stop flag has
   // been set.  If so, don't do anything and return 0.

   std::pair<FILE *, CURL *> *p = static_cast<std::pair<FILE *, CURL *> *> (userp);

   if (graphics_info_t::curl_handler_stop_it_flag_set(p->second)) {
      // signal an error to the libcurl so that it will abort the transfer.
      return size_t(0);
   } else { 
      if (buffer) {
	 const char *s = static_cast<const char *> (buffer);
	 for (size_t i=0; i<nmemb; i++) 
	    fputc(s[i], p->first);
	 return nmemb; // slightly naughty, we should return the size of the
	 // data that we actually delt with.
      } else {
	 return size_t(0);
      }
   }
}
#endif /* USE_LIBCURL */

// Return the file name of the mdl mol file.
// 
// can throw a std::runtime_error
//
std::string
get_drug_mdl_via_wikipedia_and_drugbank(std::string drugname) {

   try {
      if (graphics_info_t::prefer_python) {
	 
#ifdef USE_PYTHON
	 std::string s = get_drug_via_wikipedia_and_drugbank_py(drugname);
	 if (s.empty()) {
	    std::cout << "INFO:: get_drug_via_wikipedia result-not-a-string" << std::endl;
	 }
	 return s;
#endif

      } else {

#ifdef USE_GUILE
	 std::string s = get_drug_via_wikipedia_and_drugbank_scm(drugname);
	 if (s.empty()) {
	    std::cout << "INFO:: get_drug_via_wikipedia_scm result-not-a-string" << std::endl;
	 }
	 return s;
#endif
      }
   }
   catch (...) {
      std::cout << "WARNING:: trapped an exception in get_drug_mdl_via_wikipedia_and_drugbank() "
		<< std::endl;
   }
   return std::string("");
}

#ifdef USE_GUILE

#include "c-interface-scm.hh" // for debugging

// Return the file name of the mdl mol file.
// 
// Return empty string on error
//
std::string
get_drug_via_wikipedia_and_drugbank_scm(const std::string &drugname) {

   std::string s;
   std::string command = "(get-drug-via-wikipedia ";
   command += single_quote(drugname);
   command += ")";
   SCM r = safe_scheme_command(command);
   std::string sr = scm_to_locale_string(display_scm(r));
   if (scm_is_true(scm_string_p(r)))
      s = scm_to_locale_string(r);
   return s;
}
#endif

#ifdef USE_PYTHON
// Return the file name of the mdl mol file.
// 
// Return empty string on error
//
std::string
get_drug_via_wikipedia_and_drugbank_py(const std::string &drugname) {

   std::string s;
   std::string command = "get_drug_via_wikipedia(";
   command += single_quote(drugname);
   command += ")";
   PyObject *r = safe_python_command_with_return(command);
   if (PyString_Check(r))
      s = PyString_AsString(r);
   Py_XDECREF(r);
   return s;
}
#endif



#ifdef USE_LIBCURL
#ifdef USE_GUILE
SCM curl_progress_info(const char *file_name) {

   SCM r = SCM_BOOL_F;
   graphics_info_t g;
   std::string f(file_name);
   CURLINFO info;
   double dv;
   long int liv;
   CURL *c = g.get_curl_handle_for_file_name(f);

   if (c) { 

      r = SCM_EOL;
      info = CURLINFO_CONTENT_LENGTH_DOWNLOAD;
      CURLcode status = curl_easy_getinfo(c, info, &dv);
      if (status == CURLE_OK) {
	 SCM scm_v   = scm_double2num(dv);
	 SCM scm_sym = scm_str2symbol("content-length-download");
	 r = scm_cons(scm_cons(scm_sym, scm_v), r);
      }

      info = CURLINFO_SIZE_DOWNLOAD;
      status = curl_easy_getinfo(c, info, &dv);
      if (status == CURLE_OK) {
	 SCM scm_v   = scm_double2num(dv);
	 SCM scm_sym = scm_str2symbol("size-download");
	 r = scm_cons(scm_cons(scm_sym, scm_v), r);
      }
   } else {
      // std::cout << "Found no CURL handle for  " << f << std::endl;
   } 
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
/* this shall be a dictionary or False*/
PyObject *curl_progress_info_py(const char *file_name) {

   PyObject *r = Py_False;
   graphics_info_t g;
   std::string f(file_name);
   CURLINFO info;
   double dv;
   long int liv;
   CURL *c = g.get_curl_handle_for_file_name(f);

   if (c) { 

     r = PyDict_New();
     info = CURLINFO_CONTENT_LENGTH_DOWNLOAD;
     CURLcode status = curl_easy_getinfo(c, info, &dv);
     if (status == CURLE_OK) {
       PyObject *py_v   = PyFloat_FromDouble(dv);
       PyObject *py_sym = PyString_FromString("content-length-download");
       PyDict_SetItem(r, py_sym, py_v);
     }

     info = CURLINFO_SIZE_DOWNLOAD;
     status = curl_easy_getinfo(c, info, &dv);
     if (status == CURLE_OK) {
       PyObject *py_v   = PyFloat_FromDouble(dv);
       PyObject *py_sym = PyString_FromString("size-download");
       PyDict_SetItem(r, py_sym, py_v);
     }
   } else {
      // std::cout << "Found no CURL handle for  " << f << std::endl;
   } 
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON
#endif /* USE_LIBCURL */


#ifdef USE_LIBCURL
void stop_curl_download(const char *file_name) {  // stop curling the to file_name;

   graphics_info_t g;
   g.set_stop_curl_download_flag(file_name);

} 
#endif // USE_LIBCURL



#ifdef USE_LIBCURL
void curl_test_make_a_post() {

#ifdef NO_ANOMYMOUS_DATA
   // don't do anything
#else
   
   CURL *c = curl_easy_init();
   long int no_signal = 1; // for multi-threading

   std::string url = "http://localhost/test/cootpost.py/slurp";
   std::string post_string = "name=";
   post_string += "#1234";
   post_string += "&version=";
   post_string += VERSION;
   post_string += "&sys_build_type=";
   post_string += COOT_SYS_BUILD_TYPE;

   std::cout << "posting " << post_string << std::endl;
   std::cout << "posting to  " << url << std::endl;
   curl_easy_setopt(c, CURLOPT_NOSIGNAL, no_signal);
   curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 10);
   curl_easy_setopt(c, CURLOPT_URL, url.c_str());
   curl_easy_setopt(c, CURLOPT_POSTFIELDS, post_string.c_str());

   CURLcode status = curl_easy_perform(c);
   if (status != CURLE_OK)
      std::cout << "curl_make_a_post() failed " << curl_easy_strerror(status)
		<< std::endl;
   curl_easy_cleanup(c);
#endif
}
#endif // USE_LIBCURL



#ifdef USE_LIBCURL
void
curl_post(const std::string &url, const std::string &post_string) {
   
   CURL *c = curl_easy_init();
   long int no_signal = 1; // for multi-threading
   curl_easy_setopt(c, CURLOPT_NOSIGNAL, no_signal);
   curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 10);
   curl_easy_setopt(c, CURLOPT_URL, url.c_str());
   curl_easy_setopt(c, CURLOPT_POSTFIELDS, post_string.c_str());

   CURLcode status = curl_easy_perform(c);
   if (status != CURLE_OK)
      std::cout << "curl_post() failed " << curl_easy_strerror(status)
		<< std::endl;
   curl_easy_cleanup(c);
}
#endif // USE_LIBCURL
