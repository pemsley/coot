/* layla/python-utils.cpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
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
#include "utils/coot-utils.hh"
#include "python_utils.hpp"


// Include curl
#ifdef USE_LIBCURL
#ifndef HAVE_CURL_H
#define HAVE_CURL_H
// defined in new python!?
#ifdef socklen_t
#undef socklen_t
#endif
#include <curl/curl.h>
#endif // HAVE_CURL_H
#endif

// Let's try again just using curl. This is not python so should not be in a python file
// However, it is namespace coot::layla, and there is nowhere else for those functions
// to go at the moment. Just rename the file if/when this works, I guess
// 
std::string
coot::layla::get_drug_via_wikipedia_and_drugbank_curl(const std::string &drugname_in) {

   auto write_coot_curl_data = +[] (void *buffer, size_t size, size_t nmemb, void *userp) {

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
   };

   auto get_url_as_string = [write_coot_curl_data] (const std::string &url) {

      std::string user_agent = "coot";
      user_agent += " ";
      user_agent += VERSION;
      user_agent += " https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/";
      char buff[1024];
      std::string s;

      long int no_signal = 1; 
      CURL *c = curl_easy_init();
      curl_easy_setopt(c, CURLOPT_URL, url.c_str());
      curl_easy_setopt(c, CURLOPT_NOSIGNAL, no_signal);
      curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 4);
      curl_easy_setopt(c, CURLOPT_SSL_VERIFYPEER, 0);
      curl_easy_setopt(c, CURLOPT_USERAGENT, user_agent.c_str());
      curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, true); // 20230919-PE new, so that fetch from DrugBank works.
      CURLcode cc = curl_easy_setopt(c, CURLOPT_ERRORBUFFER, buff);
      curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_coot_curl_data);
      curl_easy_setopt(c, CURLOPT_WRITEDATA, &s);
      CURLcode success = curl_easy_perform(c);
      // std::cout << "DEBUG:: curl_easy_perform() for " << url << " success " << success << std::endl;
      if (success != 0) {
         std::cout << "WARNING:: coot_get_url_as_string with arg " << url << " failed" << std::endl;
         std::cout << "ERROR: " << buff << std::endl;
      }
      curl_easy_cleanup(c);
      // std::cout << "DEBUG:: get_url_as_string() result size " << s.length() << std::endl;
      return s;
   };

   auto write_file = [] (const std::string &s, const std::string &file_name) {
      std::ofstream f(file_name.c_str());
      f << s;
      f.close();
   };

   auto get_drugbank_code_from_html = [] (const std::string &s) {
      // hacky function - use an xml parser instead/as well

      std::string r;
      // std::string::size_type p1 = s.find("| DrugBank ");

      std::stringstream ss(s);
      std::string line;
      while (std::getline(ss, line)) {
         if (line.find(" DrugBank ") != std::string::npos) {
            if (line.length() < 80) {
               std::vector<std::string> parts = coot::util::split_string_no_blanks(line);
               if (parts.size() == 4) {
                  r = parts[3];
               }
            }
         }
      }
      return r;
   };

   std::string r;
   std::string drug_name = coot::util::downcase(drugname_in);
   std::string url_pre  = "https://en.wikipedia.org/w/api.php?format=xml&action=query&titles=";
   std::string url_post = "&prop=revisions&rvprop=content";
   std::string url = url_pre + drug_name + url_post;
   // std::cout << "debug:: url:: " << url << std::endl;
   std::string result = get_url_as_string(url);

   if (result.length() > 10) {
      std::string::size_type p1 = result.find(">#REDIRECT ");
      if (p1 != std::string::npos) {
         // Redirected. find the new drug name
         std::string s1 = result.substr(p1+13);
         // std::cout << "s1: \"" << s1 << "\"" << std::endl;
         std::string::size_type p2 = s1.find("]]");
         if (p2 != std::string::npos) {
            std::string s2 = s1.substr(0,p2);
            // std::cout << "s2: \"" << s2 << "\"" << std::endl;
            drug_name = s2;
            url = url_pre + drug_name + url_post;
            result = get_url_as_string(url);
         }
      }
      // std::cout << "result: \"" << result << "\"" << std::endl;

      std::string db_code = get_drugbank_code_from_html(result);
      std::cout << "DEBUG:: db_code: " << db_code << std::endl;

      std::string db_mol_url_pre = "https://www.drugbank.ca/structures/small_molecule_drugs/";
      std::string db_mol_url_post = ".mol";
      std::string db_mol_url = db_mol_url_pre + db_code + db_mol_url_post;
      std::string db_result = get_url_as_string(db_mol_url);

      // std::cout << "DEBUG:: db_mol_url: " << db_mol_url << std::endl;
      // std::cout << "DEBUG:: db_result \"" << db_result << "\"" << std::endl;
      if (db_result.length() > 10) {
         std::string db_file_name = db_code + ".mol";
         write_file(db_result, db_file_name);
         r = db_file_name;
      }
   }

   return r;
}



#if 0 // Python - bleugh.

void coot::layla::setup_python_basic(int argc, char **argv) {

#ifdef USE_PYTHON

#ifdef USE_PYMAC_INIT

                   asdlkfasdf
  //  (on Mac OS, call PyMac_Initialize() instead)
  //http://www.python.org/doc/current/ext/embedding.html
  //
  PyMac_Initialize();

#else

   wchar_t** _argv = static_cast<wchar_t **>(PyMem_Malloc(sizeof(wchar_t*)*argc));
   for (int i=0; i<argc; i++) {
      wchar_t* arg = Py_DecodeLocale(argv[i], NULL);
      _argv[i] = arg;
   }
   Py_InitializeEx(0);
   PySys_SetArgv(argc, _argv);

   // We expect these to be null because we are outside a python script.
   PyObject *globals = PyEval_GetGlobals();
   // std::cout << "in setup_python_basic() globals " << globals << std::endl;
   PyObject *locals  = PyEval_GetLocals();
   // std::cout << "in setup_python_basic() locals " << locals << std::endl;
   

#endif // USE_PYMAC_INIT

   auto get_pythondir = [] () {
                           std::string p = coot::prefix_dir();
                           std::string dp   = coot::util::append_dir_dir(p,   "lib");
                           std::string python_version = "python";
                           python_version += coot::util::int_to_string(PY_MAJOR_VERSION);
                           python_version += ".";
                           python_version += coot::util::int_to_string(PY_MINOR_VERSION);
                           std::string ddp  = coot::util::append_dir_dir(dp,  python_version);
                           std::string dddp = coot::util::append_dir_dir(ddp, "site-packages");
                           return dddp;
                        };
   auto get_pkgpythondir = [get_pythondir] () {
                              std::string d = get_pythondir();
                              std::string dp   = coot::util::append_dir_dir(d, "coot");
                              return dp;
                           };

   // std::string pkgpydirectory = PKGPYTHONDIR;
   // std::string pydirectory = PYTHONDIR;
   // use ${prefix}/lib/python3.9/site-package for PYTHONDIR
   // use ${pythondir}/coot' for PKGPYTHONDIR (i.e. PYTHONDIR + "/coot")

   std::string pkgpydirectory = get_pkgpythondir();
   std::string    pydirectory = get_pythondir();

   if (true) {
      std::cout << "debug:: in setup_python()    pydirectory is " << pydirectory << std::endl;
      std::cout << "debug:: in setup_python() pkgpydirectory is " << pkgpydirectory << std::endl;
   }

   PyObject *sys_path = PySys_GetObject("path");
   PyList_Append(sys_path, PyUnicode_FromString(pydirectory.c_str()));
   PyList_Append(sys_path, PyUnicode_FromString(pkgpydirectory.c_str()));

   // int err = PyRun_SimpleString("import coot");

#endif // USE_PYTHON

}

void coot::layla::setup_python_coot_module() {

   PyObject *coot = PyImport_ImportModule("coot");
   if (! coot) {
      std::cout << "ERROR:: setup_python_coot_module() Null coot" << std::endl;
   } else {
      if (false)
         std::cout << "INFO:: setup_python_coot_module() good coot module" << std::endl;
   }
}


void coot::layla::setup_python_module(const std::string &module_name) {

   PyObject *coot = PyImport_ImportModule(module_name.c_str());
   if (! coot) {
      std::cout << "ERROR:: setup_python_coot_module() Null module for " << module_name << std::endl;
   } else {
      if (true)
         std::cout << "INFO:: setup_python_coot_module() good module " << module_name << std::endl;
   }
}


PyObject *coot::layla::safe_python_command_with_return(const std::string &python_cmd) {

   std::cout << "--------------- start layla_safe_python_command_with_return(): " << python_cmd << std::endl;

   // 20220330-PE I think that this is super ricketty now!
   // Does it only find things in dynamic_atom_overlaps_and_other_outliers module?
   // this function was empty before today, returning NULL.

   // std::cout << "in safe_python_command_with_return() A " << python_cmd << std::endl;

   // command = "import coot; " + python_cmd;
   std::string command = python_cmd;

   PyObject* result = nullptr;
   PyObject *am = PyImport_AddModule("__main__");

   if (am) {
      PyObject* d = PyModule_GetDict(am);

      std::cout << "running command: " << command << std::endl;
      PyObject* source_code = Py_CompileString(command.c_str(), "adhoc", Py_eval_input);
      PyObject* func = PyFunction_New(source_code, d);
      result = PyObject_CallObject(func, PyTuple_New(0));
      std::cout << "--------------- in safe_python_command_with_return() result at: " << result << std::endl;
      if (result) {
         if(!PyUnicode_Check(result)) {
             std::cout << "--------------- in safe_python_command_with_return() result is probably not a string." << std::endl;
         }
      } else {
         std::cout << "--------------- in safe_python_command_with_return() result was null" << std::endl;
         if(PyErr_Occurred()) {
            std::cout << "--------------- in safe_python_command_with_return() Printing Python exception:" << std::endl;
            PyErr_Print();
         }
      }

      // debugging
      // PyRun_String("import coot; print(dir(coot))", Py_file_input, d, d);
      Py_XDECREF(func);
      Py_XDECREF(source_code);
   } else {
      std::cout << "ERROR:: Hopeless failure: module for __main__ is null" << std::endl;
   }
   return result;
}


std::string coot::layla::get_drug_via_wikipedia_and_drugbank_py(const std::string &drugname) {

   std::string s;

   setup_python_basic(0,0);
   setup_python_module("coot");
   setup_python_module("coot_utils");

   std::string command = "coot_utils.fetch_drug_via_wikipedia(";
   command += coot::util::single_quote(drugname);
   command += ")";
   PyObject *r = safe_python_command_with_return(command);
   if(!r) {
      std::cout<<"FIXME: Call to Python get_drug_via_wikipedia('" << drugname << "') returned a null pointer.\n";
      return s;
   }
   if (PyUnicode_Check(r))
     s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(r));
   Py_XDECREF(r);
   return s;
}
#endif // Python stuff
