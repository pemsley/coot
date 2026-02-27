
/* src/c-interface-network.cc
 * 
 * Copyright 2009, 2012 by The University of Oxford
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"


#include <iostream>
#include <stdexcept>
#include <filesystem>

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

#include "read-molecule.hh" // now with std::string args
#include <thread>
#include <chrono>
#include <future>
#include "gtk-utils.hh"
#include "curl-utils.hh"
#include <zlib.h>

#include "coot-utils/pae.hh"

#include "utils/logging.hh"
extern logging logger;


// return 0 on success
#ifdef USE_LIBCURL
int coot_get_url(const std::string &url, const std::string  &file_name) {

   std::optional<ProgressNotifier> notifier = std::nullopt;
   return coot_get_url_and_activate_curl_hook(url, file_name, 0, notifier);
}

int coot_get_url_with_notifier(const std::string &url,
                               const std::string  &file_name,
                               std::optional<ProgressNotifier> notifier) {

   return coot_get_url_and_activate_curl_hook(url, file_name, 0, notifier);
}
#endif /* USE_LIBCURL */


#ifdef USE_LIBCURL
int coot_curl_progress_callback(void *clientp,
                                curl_off_t dltotal,
                                curl_off_t dlnow,
                                curl_off_t ultotal,
                                curl_off_t ulnow) {

   ProgressNotifier* notifier_ptr = (ProgressNotifier*)(clientp);
   if(dltotal == 0) {
      dltotal++;
   }
   // g_debug("Inside coot_curl_progress_callback(); dlnow=%li, dltotal=%li", dlnow, dltotal);
   notifier_ptr->update_progress((float)dlnow/(float)dltotal);
   return 0;
}

int coot_get_url_and_activate_curl_hook(const std::string &url,
                                        const std::string &file_name,
					short int activate_curl_hook_flag,
                                        std::optional<ProgressNotifier> notifier) {

   std::cout << "DEBUG:: in coot_get_url_and_activate_curl_hook():\n url:       "
	     << url << "\n file-name: " << file_name << std::endl;

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


   auto is_html = [] (const std::string &file_name) {
	 std::ifstream f(file_name.c_str());
	 std::string line;
	 std::getline(f, line);
         if (line.find("html") != std::string::npos)
            return true;
         return false;
   };

   // 20250526-PE
   // github returns a file with contents "404: Not Found" (no newline) when the file is not
   // on the server.
   // Gnarly.
   auto is_error_message = [] (const std::string &file_name) {

#if 0 // I don't understand why this doesn't work (size is 0)
      std::ifstream f(file_name.c_str());
      if (f) {
         std::filebuf *frdbuf = f.rdbuf();
         std::size_t size = frdbuf->pubseekoff(0, f.end, f.in);
         frdbuf->pubseekpos(0, f.in);
         std::cout << "............. size " << size << std::endl;
         std::string file_contents;
         if (file_contents.find("404: ") != std::string::npos) {
            std::cout << "............... found the 404" << std::endl;
            return true;
         } else {
            std::cout << "................... no find the 404" << std::endl;
         }
         std::cout << "................. file_contents: " << file_contents << "\n" << std::endl;
      }
      return false;
#endif

      bool status = false;
      std::filesystem::path p(file_name);
      std::uintmax_t file_size = std::filesystem::file_size(p);
      if (file_size < 10)
         if (file_size < 0)
            status = true;
      return status;

   };


   auto file_to_string = [] (const std::string &file_name) {
      std::string s;
      std::string line;
      std::ifstream f(file_name);
      if (!f) {
         std::cout << "Failed to open " << file_name << std::endl;
      } else {
         std::cout << "file_to_string() A " << std::endl;
         while (std::getline(f, line)) {
            std::cout << "file_to_string() B " << line << std::endl;
            s += line;
            s += "\n";
         }
         std::cout << "file_to_string() C (while terminated) " << std::endl;
      }
      return s;
   };

   // write binary
   FILE *f = fopen(file_name.c_str(), "wb");

   if (f) {

      // If you see a crash in curl_easy_init() then the problem is
      // somewhere else.  Curl_open() doesn't do anything except a bit
      // of mallocing.  So the memory is messed up elsewhere and beforehand.
      CURL *c = curl_easy_init();
      long int no_signal = 1;
      int to = 301; // we can tolerate longer waits with sub-thread
      std::string ext = coot::util::file_name_extension(file_name);
      if (ext == ".gz") {
         std::string ext_2 = coot::util::file_name_extension(coot::util::name_sans_extension(file_name));
         if (ext_2 == ".xml")
            to = 6; // bored of a waiting for slow EBI server
      }

      bool debug = false;
      std::pair<FILE *, CURL *> p_for_write(f,c);
      curl_easy_setopt(c, CURLOPT_URL, url.c_str());
      curl_easy_setopt(c, CURLOPT_NOSIGNAL, no_signal);
      curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 6);
      // 20231004-PE remove the timeout
      // curl_easy_setopt(c, CURLOPT_TIMEOUT, to); // maximum time the request is allowed to take
      curl_easy_setopt(c, CURLOPT_SSL_VERIFYPEER, FALSE);
      std::string user_agent_str = "Coot-";
      user_agent_str += VERSION;
      user_agent_str += " https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/";
      curl_easy_setopt(c, CURLOPT_USERAGENT, user_agent_str.c_str());
      curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_coot_curl_data_to_file);
      curl_easy_setopt(c, CURLOPT_WRITEDATA, &p_for_write);
      if(notifier.has_value()) {
         auto* notifier_ptr = &notifier.value();
         curl_easy_setopt(c, CURLOPT_XFERINFOFUNCTION, coot_curl_progress_callback);
         curl_easy_setopt(c, CURLOPT_XFERINFODATA, notifier_ptr);
         curl_easy_setopt(c, CURLOPT_NOPROGRESS, 0);
      }
      std::pair <CURL *, std::string> p(c,file_name);
      CURLcode success = CURLcode(-1);

      if (debug)
         std::cout << "debug:: here in coot_get_url_and_activate_curl_hook() with activate_curl_hook_flag "
                   << activate_curl_hook_flag << std::endl;

      if (activate_curl_hook_flag) {
	 graphics_info_t g;
	 g.add_curl_handle_and_file_name(p);
#ifdef USE_GUILE
	 success = CURLcode(GPOINTER_TO_INT(scm_without_guile(wrapped_curl_easy_perform, c)));
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
         // 20250526-PE unfortunately curl_easy_perform() is returning 0 even
         // if the file was not on the server.
	 success = curl_easy_perform(c);

         if (coot::file_exists(file_name)) {
            std::cout << "DEBUG:: file_exists " << file_name << std::endl;
            std::filesystem::path fp(file_name);
            std::uintmax_t size = std::filesystem::file_size(fp);
            std::cout << "DEBUG:: file_size " << size << std::endl;
            if (is_html(file_name)) {
               std::cout << "DEBUG:: file is html " << file_name << std::endl;
               success = CURLcode(23); // CURL write error (say)
               int rm_status = remove(file_name.c_str()); // Ciao Bella
               if (rm_status == 0)
                  logger.log(log_t::INFO, logging::function_name_t("coot_get_url_and_activate_curl_hook"),
                             file_name, "removed (failed download)");
            } else {
               if (debug)
                  std::cout << "::::::::::::::: file is not html " << file_name << std::endl;
               if (is_error_message(file_name)) {
                  if (debug)
                     std::cout << "::::::::::::::: file is error message " << file_name << std::endl;
                  success = CURLcode(23); // CURL write error (say)
                  int rm_status = remove(file_name.c_str()); // Ciao Bella
                  if (rm_status == 0)
                     logger.log(log_t::INFO, logging::function_name_t("coot_get_url_and_activate_curl_hook"),
                                file_name, "removed (failed download)");
               } else {
                  if (debug)
                     std::cout << "::::::::::::::: file is not error message " << file_name << std::endl;
               }
            }
         } else {
            success = CURLcode(23); // CURL write error (say)
         }
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
   user_agent += " https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/";
   char buff[1024];

   long int no_signal = 1;
   CURL *c = curl_easy_init();
   curl_easy_setopt(c, CURLOPT_URL, url);
   curl_easy_setopt(c, CURLOPT_NOSIGNAL, no_signal);
   curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 4);
   curl_easy_setopt(c, CURLOPT_SSL_VERIFYPEER, FALSE);
   curl_easy_setopt(c, CURLOPT_USERAGENT, user_agent.c_str());
   CURLcode cc = curl_easy_setopt(c, CURLOPT_ERRORBUFFER, buff);
   curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_coot_curl_data);
   curl_easy_setopt(c, CURLOPT_WRITEDATA, &s);
   CURLcode success = curl_easy_perform(c);
   if (success != 0) {
      std::cout << "WARNING:: coot_get_url_as_string_internal with arg " << url << " failed" << std::endl;
      std::cout << "ERROR: " << buff << std::endl;
   }
   curl_easy_cleanup(c);
   // std::cout << "DEBUG:: coot_get_url_as_string_internal() returns this string\n---start---\n" << s << "\n---end---\n"<< std::endl;
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

   PyObject *r = Py_False;
   std::string s = coot_get_url_as_string_internal(url);
   r = PyUnicode_FromString(s.c_str());
   if (r) {
      if (PyUnicode_Check(r)) {
         // all good
      } else {
         std::cout << "DEBUG:: coot_get_url_as_string_py(): A response was not unicode" << std::endl;
         PyErr_Clear(); // clear the above error
         r = PyBytes_FromStringAndSize(s.c_str(), s.size());
      }
      if (r) {
         if (PyBool_Check(r)) {
           Py_INCREF(r);
         }
      }
   } else {
      std::cout << "DEBUG:: coot_get_url_as_string_py(): B response was not unicode" << std::endl;
      PyErr_Clear(); // clear the above error
      r = PyBytes_FromStringAndSize(s.c_str(), s.size());
   }
   return r;
}
#endif /* USE_PYTHON */
#endif /* USE_LIBCURL */


#ifdef USE_LIBCURL
size_t
write_coot_curl_data(void *buffer, size_t size, size_t nmemb, void *user_data) {

   // don't stop on null if we have binary data
   size_t realsize = size * nmemb;
   std::string *response = static_cast<std::string*>(user_data);
   char *ptr = static_cast<char *>(buffer);
   response->append(ptr, realsize);
   return realsize;  // continue reading
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
	    // std::cout << "INFO:: get_drug_via_wikipedia result-not-a-string" << std::endl;
	    logger.log(log_t::INFO, "get_drug_via_wikipedia result-not-a-string");
	 }
	 return s;
#endif

      } else {

#ifdef USE_GUILE
	 std::string s = get_drug_via_wikipedia_and_drugbank_scm(drugname);
	 if (s.empty()) {
	    // std::cout << "INFO:: get_drug_via_wikipedia_scm result-not-a-string" << std::endl;
	    logger.log(log_t::INFO, "get_drug_via_wikipedia_scm result-not-a-string");
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
   std::string command = "coot_utils.fetch_drug_via_wikipedia(";
   command += single_quote(drugname);
   command += ")";
   PyObject *r = safe_python_command_with_return(command);
   if(!r) {
      std::cout<<"fixme: Call to Python get_drug_via_wikipedia('"<<drugname<<"') returned a null pointer.\n";
      return s;
   }
   if (PyUnicode_Check(r))
     s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(r));
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
	 SCM scm_v   = scm_from_double(dv);
	 SCM scm_sym = scm_string_to_symbol(scm_from_locale_string("content-length-download"));
	 r = scm_cons(scm_cons(scm_sym, scm_v), r);
      }

      info = CURLINFO_SIZE_DOWNLOAD;
      status = curl_easy_getinfo(c, info, &dv);
      if (status == CURLE_OK) {
	 SCM scm_v   = scm_from_double(dv);
	 SCM scm_sym = scm_string_to_symbol(scm_from_locale_string("size-download"));
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
       PyObject *py_sym = myPyString_FromString("content-length-download");
       PyDict_SetItem(r, py_sym, py_v);
     }

     info = CURLINFO_SIZE_DOWNLOAD;
     status = curl_easy_getinfo(c, info, &dv);
     if (status == CURLE_OK) {
       PyObject *py_v   = PyFloat_FromDouble(dv);
       PyObject *py_sym = myPyString_FromString("size-download");
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

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"

void fetch_and_superpose_alphafold_models_using_active_molecule() {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = active_atom_spec();
   if (aa.first) {
      int imol = aa.second.first;
      fetch_and_superpose_alphafold_models(imol);
   }
}

void fetch_and_superpose_alphafold_models(int imol) {

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {
         bool found_a_uniprot_dbref = false;
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_refs = chain_p->GetNumberOfDBRefs();
               std::string chain_id = chain_p->GetChainID();
               for (int ref_no=0; ref_no<n_refs; ref_no++) {
                  mmdb::DBReference  *ref = chain_p->GetDBRef(ref_no);  // 0..nDBRefs-1
                  std::string db = ref->database;
                  std::string db_accession = ref->dbAccession;
                  // std::cout << "INFO:: DBREF Chain " << chain_id << " " << db << " " << db_accession << std::endl;
                  logger.log(log_t::INFO, "DBREF Chain", chain_id, db, db_accession);
                  if (db == "UNP") {  // uniprot
                     found_a_uniprot_dbref = true;
                     int imol_af = fetch_alphafold_model_for_uniprot_id(db_accession);
                     if (is_valid_model_molecule(imol_af)) {
                        superpose_with_chain_selection(imol, imol_af, chain_id.c_str(), "A", 1, 0, 0);
                        graphics_info_t g;
                        g.molecules[imol_af].ca_plus_ligands_representation(g.Geom_p(), true);
                        graphics_draw();
                     }
                  }
               }
            }
         }
         if (!found_a_uniprot_dbref) {
            std::string m = "WARNING:: no DBREF found in molecule header";
            info_dialog(m.c_str());
            std::cout << m << " " << imol << std::endl;
         }
      }
   }
}

void display_png_from_string_in_a_dialog(const std::string &string, const std::string &title);

int
fetch_alphafold_model_for_uniprot_id(const std::string &uniprot_id) {

   auto join = [] (const std::string &d, const std::string &f) {
      return d + std::string("/") + f;
   };

   auto write_string_to_file = [] (const std::string &s, const std::string &fn) {

      std::cout << "write string to file! " << fn << std::endl;
      std::ofstream f(fn);
      f << s;
      f << "\n";
      f.close();
   };

   int imol = -1; // return this
   // https://alphafold.ebi.ac.uk/files/AF-Q7N8I7-F1-model_v4.pdb
   // https://alphafold.ebi.ac.uk/files/AF-Q7N8I7-F1-predicted_aligned_error_v4.json
   std::string fn_tail_pdb = std::string("AF-") + uniprot_id + std::string("-F1-model_v6.pdb");
   std::string fn_tail_pae = std::string("AF-") + uniprot_id + std::string("-F1-predicted_aligned_error_v6.json");

   xdg_t xdg;
   std::string download_dir = join(xdg.get_cache_home().string(), "coot-download");
   // make coot-download if needed
   std::string dld = coot::get_directory(download_dir);
   if (! dld.empty()) {
      download_dir = dld;
      std::string fn_pdb = coot::util::append_dir_file(download_dir, fn_tail_pdb);
      std::string url = std::string("https://alphafold.ebi.ac.uk/files/") + fn_tail_pdb;
      bool needs_downloading = true;
      if (coot::file_exists_and_non_tiny(fn_pdb, 500)) needs_downloading = false;
      if (needs_downloading) {
         coot_get_url(url.c_str(), fn_pdb.c_str());
         if (coot::file_exists_and_non_tiny(fn_pdb, 500)) {
            imol = handle_read_draw_molecule_and_move_molecule_here(fn_pdb);
         } else {
            std::string m = "WARNING:: UniProt ID " + uniprot_id + std::string(" not found");
            info_dialog(m.c_str());
         }
      }

      // display_pae_from_file_in_a_dialog() needs the model to exist
      // (so that it knows which model to manipulate)
      imol = handle_read_draw_molecule_and_move_molecule_here(fn_pdb);
      graphics_draw();

      // now let's do the PAE file
      std::string fn_pae = coot::util::append_dir_file(download_dir, fn_tail_pae);
      url = std::string("https://alphafold.ebi.ac.uk/files/") + fn_tail_pae;
      needs_downloading = true;
      if (coot::file_exists_and_non_tiny(fn_pae, 500)) needs_downloading = false;
      if (needs_downloading) {
         coot_get_url(url.c_str(), fn_pae.c_str());
         if (coot::file_exists_and_non_tiny(fn_pae, 500)) {
            display_pae_from_file_in_a_dialog(imol, fn_pae);
         }
      } else {
         display_pae_from_file_in_a_dialog(imol, fn_pae);
      }

   }
   return imol;
}



#ifdef USE_LIBCURL
void stop_curl_download(const char *file_name) {  // stop curling the to file_name;

   graphics_info_t g;
   g.set_stop_curl_download_flag(file_name);

}
#endif // USE_LIBCURL

#ifdef USE_LIBCURL
void fetch_emdb_map(const std::string &emd_accession_code) {

   std::string map_gz_url = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-" +
      emd_accession_code + "/map/emd_" + emd_accession_code + ".map.gz";
   std::string download_dir = "coot-download";
   download_dir = coot::get_directory(download_dir.c_str());
   std::string gz_fnl = "emd_" + emd_accession_code + ".map.gz";
   std::string fnl    = "emd_" + emd_accession_code + ".map";
   std::string gz_fn = coot::util::append_dir_file(download_dir, gz_fnl);
   std::string fn    = coot::util::append_dir_file(download_dir,    fnl);

   if (coot::file_exists_and_non_tiny(fn)) {
      read_ccp4_map(fn, false);
      g_info("Reading CCP4 map from cached downloads...");
      return;
   }

   std::string label_str = "Downloading map for " + emd_accession_code + " from EMDB...";
   ProgressBarPopUp popup("Coot: Downloading Map", label_str);
   std::thread worker([=](ProgressBarPopUp&& pp){

      std::shared_ptr<ProgressBarPopUp> popup = std::make_shared<ProgressBarPopUp>(std::move(pp));
      int status = coot_get_url_with_notifier(map_gz_url, gz_fn, ProgressNotifier(popup));

      if (status != 0) { // if it's bad
         g_warning("Download failed. Status=%i", status);
         return;
      }

      std::string gzipedBytes;

      std::ifstream file(gz_fn);
      std::stringstream ss;

      while (!file.eof()) {
         gzipedBytes += (char) file.get();
      }
      file.close();
      if (gzipedBytes.size() == 0) {
         g_warning("The downloaded file (%s) is empty or could not be read.", gz_fn.c_str());
         return;
      }

      unsigned int full_length   = gzipedBytes.size();
      unsigned int half_length   = gzipedBytes.size()/2;
      unsigned int uncomp_length = full_length;
      Bytef* uncomp = new Bytef[uncomp_length];
      z_stream strm;
      strm.next_in   = (Bytef *) gzipedBytes.c_str();
      strm.avail_in  = full_length;
      strm.total_out = 0;
      strm.zalloc    = Z_NULL;
      strm.zfree     = Z_NULL;
      bool done = false;

      int err = inflateInit2(&strm, (16 + MAX_WBITS));
      if (err != Z_OK) {
         delete [] uncomp;
         g_warning("The downloaded file (%s) could not be decompressed (zlib error: %i) [1].", gz_fn.c_str(), err);
         return;
      }

      while (!done) {
         // If the output buffer is too small
         if (strm.total_out >= uncomp_length) {
            // Increase size of output buffer
            // Bytef* uncomp2 = static_cast<Bytef *>(calloc(sizeof(Bytef), uncomp_length + half_length));
            Bytef* uncomp2 = new Bytef[uncomp_length + half_length];
            memset(uncomp2, 0, uncomp_length + half_length);
            memcpy(uncomp2, uncomp, uncomp_length);
            uncomp_length += half_length;
            delete [] uncomp;
            uncomp = uncomp2;
         }

         strm.next_out  = static_cast<Bytef *>(uncomp + strm.total_out);
         strm.avail_out = uncomp_length - strm.total_out;

         // keep inflating...
         int err = inflate (&strm, Z_SYNC_FLUSH);
         if (err == Z_STREAM_END) {
            done = true;
         } else if (err != Z_OK) {
            g_warning("The downloaded file (%s) could not be decompressed (zlib error: %i) [2].", gz_fn.c_str(), err);
            break;
         }
      }

      err = inflateEnd (&strm);
      if (err != Z_OK) {
         delete [] uncomp;
         g_warning("The downloaded file (%s) could not be decompressed (zlib error: %i) [3].", gz_fn.c_str(), err);
         return;
      }

      for (size_t i = 0; i < strm.total_out; ++i) {
         ss << uncomp[i];
      }
      delete [] uncomp;

      g_info("The downloaded file has been successfully decompressed. Writing it down...");
      std::ofstream out(fn);
      out << ss.str();
      out.close();
      g_info("Deleting the downloaded archive...");
      remove(gz_fn.c_str());

      struct callback_data {
         std::string fn;
         std::shared_ptr<ProgressBarPopUp> popup;
      };
      callback_data* cbd = new callback_data{fn, std::move(popup)};
      g_idle_add(+[](gpointer user_data) {
         callback_data* cbd = (callback_data*) user_data;
         g_info("Reading CCP4 map from downloaded file...");
         int imol = read_ccp4_map(cbd->fn, false);
         go_to_map_molecule_centre(imol);
         delete cbd;
         return FALSE; // only one
      }, cbd);

   }, std::move(popup));
   worker.detach();


}
#endif // USE_LIBCURL

#ifdef USE_LIBCURL
int fetch_cod_entry(const std::string &cod_code) {

   int imol = -1;
   std::string url = "https://www.crystallography.net/cod/" + cod_code + ".cif";
   std::cout << "url: " << url << std::endl;

   xdg_t xdg;
   std::filesystem::path ch = xdg.get_cache_home();
   if (! std::filesystem::exists(ch))
      std::filesystem::create_directories(ch);
   std::filesystem::path download_dir = ch / "coot-download";
   if (! std::filesystem::exists(download_dir))
      std::filesystem::create_directories(download_dir);

   std::string fn_tail = cod_code + std::string(".cif");
   std::filesystem::path fn = download_dir / fn_tail;
   if (std::filesystem::exists(fn)) {
      imol = read_small_molecule_cif(fn.c_str());
   } else {
      coot_get_url(url.c_str(), fn.c_str());
      if (coot::file_exists_and_non_tiny(fn.string())) {
         imol = read_small_molecule_cif(fn.c_str());
      } else {
         std::cout << "DEBUG:: failed to download " << url << std::endl;
      }
   }
   return imol;
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
   curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 6);
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
   curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 6);
   curl_easy_setopt(c, CURLOPT_URL, url.c_str());
   curl_easy_setopt(c, CURLOPT_POSTFIELDS, post_string.c_str());

   CURLcode status = curl_easy_perform(c);
   if (status != CURLE_OK)
      std::cout << "curl_post() failed " << curl_easy_strerror(status)
		<< std::endl;
   curl_easy_cleanup(c);
}
#endif // USE_LIBCURL




/*  ----------------------------------------------------------------------- */
/*                  Client/Server                                           */
/*  ----------------------------------------------------------------------- */
//! \brief client/server functions
#ifdef USE_PYTHON
void set_python_draw_function(const std::string &s) {

   graphics_info_t g;
   g.set_python_draw_function(s);

}
#endif // USE_PYTHON

#include "get-monomer.hh"

void get_monomer_dictionary_in_subthread(const std::string &comp_id,
					 bool run_get_monomer_post_fetch_flag) {

   struct cif_data {
      int fetch_done;
      bool run_get_monomer_post_fetch_flag;
      std::string comp_id;
   };

   auto delete_cif_file_if_exists = [] (std::filesystem::path cif_file_path) {

      if (std::filesystem::exists(cif_file_path)) {
         std::filesystem::remove(cif_file_path);
      }
   };

   auto get_dict = [delete_cif_file_if_exists] (const std::string &comp_id, struct cif_data *c) {

      // https://raw.githubusercontent.com/MonomerLibrary/monomers/refs/heads/master/l/L36.cif
      // put it in xdg cache monomers/l/L36.cif

      const char rs = comp_id[0];
      const char v = tolower(rs); // get the sub directory name
      std::string letter(1, v);
      std::string cif_file_name = comp_id + ".cif";

      // Bernhard says not to use std::filesystem for URLs - OK.
      std::string github = "https://raw.githubusercontent.com/MonomerLibrary/monomers/refs/heads/master";
      std::string url_string = github + "/" +  letter + "/" + cif_file_name;

      xdg_t xdg;

      std::filesystem::path ch = xdg.get_cache_home();

      if (! std::filesystem::exists(ch))
	 std::filesystem::create_directories(ch);

      if (! std::filesystem::exists(ch))
	 std::filesystem::create_directories(ch);
      std::filesystem::path monomers_path = ch / "monomers";
      if (!std::filesystem::exists(monomers_path))
	 std::filesystem::create_directories(monomers_path);
      std::filesystem::path letter_dir_path = monomers_path / letter;
      if (!std::filesystem::exists(letter_dir_path))
	 std::filesystem::create_directories(letter_dir_path);

      if (std::filesystem::exists(ch)) {
	 std::filesystem::path monomers_path = ch / "monomers";
	 if (std::filesystem::exists(monomers_path)) {
	    std::filesystem::path letter_dir_path = monomers_path / letter;
	    if (std::filesystem::exists(letter_dir_path)) {
	       std::filesystem::path cif_file_path = letter_dir_path / cif_file_name;
               // return 0 on success
	       int status = coot_get_url(url_string, cif_file_path.string());
               logger.log(log_t::DEBUG, logging::function_name_t("get_monomer_dictionary_in_subthread()"),
                          "coot_get_url() returned status ", status);
               if (status != 0)
                  delete_cif_file_if_exists(cif_file_path);
	    } else {
	       // std::cout << "INFO:: " << letter_dir_path << " does not exist " << std::endl;
	       logger.log(log_t::INFO, letter_dir_path.string(), "does not exist");
	    }
	 } else {
	    //std::cout << "INFO:: " << monomers_path << " does not exist " << std::endl;
	    logger.log(log_t::INFO, monomers_path.string(), "does not exist");
	 }
      } else {
	 // std::cout << "INFO:: " << ch << " does not exist " << std::endl;
	 logger.log(log_t::INFO, ch.string(), "does not exist");
      }
      c->fetch_done = 1;
   };


   if (! comp_id.empty()) {
      // std::cout << "INFO:: Get " << comp_id << " in a sub-thread" << std::endl;
      logger.log(log_t::INFO, "Get", comp_id, "in a sub-thread");

      struct cif_data *c = new cif_data;
      c->fetch_done = 0;
      c->comp_id = comp_id;
      c->run_get_monomer_post_fetch_flag = run_get_monomer_post_fetch_flag;

      std::thread thread(get_dict, comp_id, c);
      thread.detach();

      auto timeout_func = +[] (gpointer data) {
         int s = 1;
	 struct cif_data *cif_data = static_cast<struct cif_data *>(data);
	 if (cif_data->fetch_done == 1) s = 0;
         if (! cif_data->comp_id.empty()) {
            const char rs = cif_data->comp_id[0];
            const char v = tolower(rs); // get the sub directory name
            std::string letter(1, v);
            std::string cif_file_name = cif_data->comp_id + ".cif";
            xdg_t xdg;
            std::filesystem::path ch = xdg.get_cache_home();
            std::filesystem::path monomers_path = ch / "monomers";
            std::filesystem::path letter_dir_path = monomers_path / letter;
            std::filesystem::path cif_file_path = letter_dir_path / cif_file_name;
            // std::cout << "INFO:: call read_cif_dictionary() on this cif "
            // << cif_file_path << std::endl;
            if (std::filesystem::exists(cif_file_path)) {
               logger.log(log_t::INFO, "call read_cif_dictionary() on this cif", cif_file_path.string());
               int read_status = read_cif_dictionary(cif_file_path.string());
               if (read_status > 0) {
                  if (cif_data->run_get_monomer_post_fetch_flag)
                     get_monomer(cif_data->comp_id);
               } else {
                  logger.log(log_t::WARNING, logging::function_name_t("get_monomer_dictionary_in_subthread"),
                             "Failed to read", cif_file_path);
               }
            } else {
               // failed to download:
               logger.log(log_t::WARNING, logging::function_name_t("get_monomer_dictionary_in_subthread"),
                          "File does not exist", cif_file_path);
               std::string message = std::string("Failed to download ") + cif_file_path.string() +
                  std::string(" from github monomers");
               info_dialog(message.c_str());
            }
         }
	 return gboolean(s);
      };
      g_timeout_add(500, GSourceFunc(timeout_func), c);
   }
}
