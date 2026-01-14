
#include <Python.h>

#include <string>
#include <vector>
#include <iostream>
#include <thread>

#include <gtk/gtk.h>
#include <netinet/in.h>
#include <fcntl.h>
#include <coot-utils/json.hpp>

#include <utils/coot-utils.hh>

#include "python-3-interface.hh"
#include "python-results-container.hh"
#include "graphics-info.h"
#include "c-interface.h"

using json = nlohmann::json;

struct param_doc {
   std::string name;
   long kind;
   std::string annotation;
};

struct func_doc {
   std::string function_name;
   std::string documentation;
   std::vector<param_doc> params;
};



static int server_fd = -1;
static int client_fd = -1;
gint coot_socket_listener_idle_func(gpointer data);

void init_coot_socket_listener() {

   int port = graphics_info_t::remote_control_port_number;
   server_fd = socket(AF_INET, SOCK_STREAM, 0);
   if (server_fd < 0) {
      std::cerr << "Error: Unable to create socket\n";
      return;
   }

   sockaddr_in addr;
   std::memset(&addr, 0, sizeof(addr));
   addr.sin_family = AF_INET;
   addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK); // localhost only
   addr.sin_port = htons(port);

   int optval = 1;
   setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR, &optval, sizeof(optval));

   if (bind(server_fd, (sockaddr*)&addr, sizeof(addr)) < 0) {
      std::cerr << "Error: Unable to bind socket\n";
      close(server_fd);
      server_fd = -1;
      return;
   }
   if (listen(server_fd, 1) < 0) {
      std::cerr << "Error: Unable to listen\n";
      close(server_fd);
      server_fd = -1;
      return;
   }

   // Make server socket non-blocking
   int flags = fcntl(server_fd, F_GETFL, 0);
   fcntl(server_fd, F_SETFL, flags | O_NONBLOCK);

   // log this
   std::cout << "INFO:: Socket listener initialized on port " << port << std::endl;
}

// called by c_inner_main() if we have guile
void make_socket_listener_maybe() {

   if (graphics_info_t::try_port_listener) {
      if (graphics_info_t::coot_socket_listener_idle_function_token == -1) {
         // if (graphics_info_t::listener_socket_have_good_socket_state) {
         // 2025-12-07-PE I am not sure that that is useful in this new
         // listener.
         if (true) {
            init_coot_socket_listener();
            int timeout_interval = 1000;
            graphics_info_t::coot_socket_listener_idle_function_token =
              g_timeout_add(timeout_interval, coot_socket_listener_idle_func, nullptr);
         }
      }
   }
}

void set_coot_listener_socket_state_internal(int sock_state) {
   graphics_info_t::listener_socket_have_good_socket_state = sock_state;
}

void set_remote_control_port(int port_number) {
  graphics_info_t::remote_control_port_number = port_number;
}

int get_remote_control_port_number() {
  return graphics_info_t::remote_control_port_number;
}


gint coot_socket_listener_idle_func(gpointer data) {

   auto get_param_docs = [] (PyObject *attr) {

      // std::cout << "get_param_docs(): --- start --- " << attr << std::endl;
      std::vector<param_doc> param_doc_vec; // return this

      if (attr && PyCallable_Check(attr)) {
         PyObject *inspect_module_name = myPyString_FromString("inspect");
         PyObject *inspect = PyImport_Import(inspect_module_name);
         PyObject* signature = PyObject_CallMethod(inspect, "signature", "O", attr);

         if (!signature) {
            std::cout << "WARNING:: coot_socket_listener_idle_func(): get_param_docs(): null signature" << std::endl;
            PyErr_Print();
            PyErr_Clear();
         } else {

            PyObject* params = PyObject_GetAttrString(signature, "parameters");
            PyObject *key, *value;
            Py_ssize_t pos = 0;

            // std::cout << "DEBUG:: params type: " << params->ob_type->tp_name << std::endl;
            PyObject* empty = PyObject_GetAttrString(inspect, "_empty");

            if (params) {
                // Get keys from the mappingproxy
                PyObject* keys = PyMapping_Keys(params);
                if (keys) {

                    // PyObject *dp = display_python(keys);
                    // if (dp)
                    //    std::cout << "DEBUG:: keys: " << PyBytes_AS_STRING(PyUnicode_AsUTF8(dp)) << std::endl;
                    // else
                    //    std::cout << "DEBUG:: null dp" << std::endl;

                    Py_ssize_t size = PyList_Size(keys);

                    for (Py_ssize_t i = 0; i < size; i++) {
                        PyObject* key = PyList_GetItem(keys, i);
                        const char* param_name = PyUnicode_AsUTF8(key);

                        // PyObject* value = PyObject_GetItem(params, key);
                        PyObject* value = PyMapping_GetItemString(params, param_name);  // ← Use this instead

                        /* 
                        // std::cout << "DEBUG: processing parameter: " << param_name << std::endl;
                        {
                           if (value) {
                              // Get the type name to see what we actually have
                              const char* type_name = value->ob_type->tp_name;
                              // type_name is "Parameter"
                              // std::cout << "DEBUG:: key " << i << " name: " << param_name << " type: " << type_name << std::endl;

                              // Try to get the 'kind' attribute to verify it's a Parameter object
                              PyObject *kindObj = PyObject_GetAttrString(value, "kind");
                              if (kindObj) {
                                  long kind = PyLong_AsLong(kindObj);
                                  std::cout << "DEBUG::   -> has kind attribute: " << kind << std::endl;
                                  Py_DECREF(kindObj);
                              } else {
                                  std::cout << "DEBUG::   -> NO kind attribute (not a Parameter object!)" << std::endl;
                                  PyErr_Clear();
                              }
                           }
                        }
                        */

                        // Now you can access the Parameter object
                        PyObject *kindObj    = PyObject_GetAttrString(value, "kind");
                        PyObject *defaultObj = PyObject_GetAttrString(value, "default");
                        PyObject *annotObj   = PyObject_GetAttrString(value, "annotation");

                        long kind = PyLong_AsLong(kindObj);
                        bool has_annot = (annotObj != empty);
                        std::string annotation;
                        if (has_annot) {
                           PyObject* s = PyObject_Str(annotObj);
                           annotation = PyUnicode_AsUTF8(s);
                           Py_DECREF(s);
                        }

                        param_doc dp;
                        dp.name = param_name;
                        dp.kind = PyLong_AsLong(kindObj);
                        dp.annotation = annotation;

                        param_doc_vec.push_back(dp);

                        Py_DECREF(value);
                    }
                    Py_DECREF(keys);
                }
            }

            /* -------------- put this in another function ---------- 
            PyObject* retAnn = PyObject_GetAttrString(sig, "return_annotation");
            bool has_ret_annot = (retAnn != empty);

            std::string return_type;
            if (has_ret_annot) {
                PyObject* s = PyObject_Str(retAnn);
                return_type = PyUnicode_AsUTF8(s);
                Py_DECREF(s);
            }
            Py_DECREF(retAnn);
            */

         }
      }
      // std::cout << "DEBUG:: ----------- returning " << param_doc_vec.size() << " param docs" << std::endl;
      return param_doc_vec;
   };

   auto handle_list_tools = [get_param_docs] (int block_index) {

      // the first is the module name, e.g. coot, coot_utils and the
      // second is the list of functions in that module
      // return functions
      std::vector<std::pair<std::string, std::vector<func_doc> > > functions;
      int n_blocks = 400;

      PyObject* inspect = PyImport_ImportModule("inspect");

      unsigned int count = 0;
      std::vector<std::string> module_names = {"coot", "coot_utils"};
      for (const std::string &module_name : module_names) {
         std::vector<func_doc> func_doc_vec; // functions that are in this module
         PyObject* pModule = PyImport_ImportModule(module_name.c_str());

         if (!pModule) {
            PyErr_Print();
            std::cerr << "Failed to import module " << module_name << "\n";
            continue;
         }

         // Get list of attributes (like dir(coot))
         PyObject* pDir = PyObject_Dir(pModule);
         if (!pDir) {
             PyErr_Print();
             Py_DECREF(pModule);
             continue;
         }

         // Iterate over list
         Py_ssize_t size = PyList_Size(pDir);
         for (Py_ssize_t i=0; i<size; ++i) {

            PyObject* pName = PyList_GetItem(pDir, i); // borrowed ref
            const char* name = PyUnicode_AsUTF8(pName);

            PyObject* attr = PyObject_GetAttrString(pModule, name);
            if (PyCallable_Check(attr)) {
               func_doc fd;
               fd.function_name = name;
               PyObject* docObj = PyObject_GetAttrString(attr, "__doc__");
               if (docObj) {
                  if (docObj == Py_None) {
                     if (false)
                        std::cout << "debug:: name " << name << " docobj: pynone" << std::endl;
                  } else {
                     std::string doc = PyBytes_AS_STRING(PyUnicode_AsUTF8String(docObj));
                     if (false)
                        std::cout << "debug:: name " << name << " docobj: " << doc << std::endl;
                     if (! doc.empty())
                        fd.documentation = doc;
                  }
               }
               bool attr_is_a_real_function = false;
               if (PyCFunction_Check(attr)) attr_is_a_real_function = true;
               if (PyFunction_Check(attr))  attr_is_a_real_function = true;
               if (PyMethod_Check(attr))    attr_is_a_real_function = true;
               if (attr_is_a_real_function) {
                  std::vector<param_doc> param_docs = get_param_docs(attr);
                  if (! param_docs.empty())
                     fd.params = param_docs;

                  if (false)
                     std::cout << "DEBUG:: in handle_list_tools() n_blocks is " << n_blocks
                               << " count/n_blocks is " << count/n_blocks << " block_index is "
                               << block_index << std::endl;

                  if (count/n_blocks == block_index || block_index == -1)
                     func_doc_vec.push_back(fd);
                  Py_XDECREF(attr);
                  count++;
               }
            }
         }

         if (! func_doc_vec.empty())
            functions.push_back(std::make_pair(module_name, func_doc_vec));

         Py_DECREF(pDir);
         Py_DECREF(pModule);
      }
      return functions;
   };

   auto handle_list_tools_with_search_pattern = [handle_list_tools] (const std::string &raw_pattern) {

      // note that raw_pattern gets split into words - and each word needs to pass the filter.

      std::cout << "DEBUG:: in handle_list_tools_with_search_pattern(): " << raw_pattern << std::endl;

      std::vector<std::string> search_words = coot::util::split_string(raw_pattern, " ");

      std::vector<std::pair<std::string, std::vector<func_doc> > > functions = handle_list_tools(-1);
      std::vector<std::pair<std::string, std::vector<func_doc> > > filtered_functions;
      for (const auto &module_pair : functions) {
         const auto &module = module_pair.first;
         const auto &funcs = module_pair.second;
         std::vector<func_doc> matchers;
         for (const auto &func : funcs) {

            bool found = true;
            for (const auto &word : search_words) {
               if (func.function_name.find(word) == std::string::npos) {
                  if (func.documentation.find(word) == std::string::npos)
                     found = false;
               }
               if (! found) break;
            }
            if (found)
               matchers.push_back(func);
         }
         if (! matchers.empty())
            filtered_functions.push_back(std::make_pair(module, matchers));
      }
      return filtered_functions;
   };

   auto make_response_from_functions = [] (const std::vector<std::pair<std::string, std::vector<func_doc> > > &functions, int id, int n_max, bool add_python_exec) {

      json j_functions;
      json j_item;
      if (add_python_exec) {
         json j_item;
         j_item["name"] = "python.exec";
         j_item["description"] = "Execute Python Code";
         j_item["params"] = json::array({"code"});
         j_functions.push_back(j_item);
         json j_item_m;
         j_item_m["name"] = "python.exec_multiline";
         j_item_m["description"] = "Execute Multi-line Python Code";
         j_item_m["params"] = json::array({"code"});
         j_functions.push_back(j_item_m);
      }

      unsigned int n = 0;
      for (const auto &module_pair : functions) {
         const auto &module = module_pair.first;
         const auto &funcs = module_pair.second;
         for (const auto &func : funcs) {
            if (n >= n_max)
               if (n_max != -1)
                  continue; // tmp limit
            json j_item;
            j_item["name"] = module + "." + func.function_name;
            // std::cout << "debug:: j_item[name] is " << module + "." + func.function_name << std::endl;
            if (!func.documentation.empty()) {
               j_item["description"] = func.documentation;
               // std::cout << "DEBUG:: documentation for " << func.function_name << ":\n" << func.documentation << std::endl;
            }
            json j_params = json::array();
            if (! func.params.empty()) {
               for (unsigned int i=0; i<func.params.size(); i++) {
                  json jp;
                  const auto &param = func.params[i];
                  json j_param_name = param.name;
                  jp["name"] = j_param_name;
                  jp["kind"] = param.kind;
                  if (! param.annotation.empty()) {
                     jp["annotation"] = param.annotation;
                  }
                  j_params.push_back(jp);
               }
            }
            j_item["params"] = j_params;
            j_functions.push_back(j_item);
            n++;
         }
      }
      return j_functions;
   };

   auto handle_execute_python_result = [] (execute_python_results_container_t rrr, int id) {

      std::string response_str;
      if (rrr.result)
         if (PyBool_Check(rrr.result) || rrr.result == Py_None)
            Py_INCREF(rrr.result);
      const char *mess = "%s";
      PyObject *dest = myPyString_FromString(mess);
      PyObject *o = PyUnicode_Format(dest, rrr.result);
      if (o) {
         std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(o));
         json j_response;
         j_response["jsonrpc"] = "2.0";
         j_response["id"] = std::to_string(id);
         json j_result;
         j_result["value"] = s;
         if (! rrr.stdout.empty())
            j_result["stdout"] = rrr.stdout;
         j_response["result"] = j_result;
         if (! rrr.error_message.empty()) {
            json j_error;
            j_error["code"] = -32001;
            j_error["message"] = rrr.error_message;
            j_response["error"] = j_error;
         }
         response_str = j_response.dump();
      } 
      if (! rrr.error_message.empty()) {
         json j_response;
         j_response["jsonrpc"] = "2.0";
         j_response["id"] = std::to_string(id);
         json j_error;
         j_error["code"] = -32001;
         j_error["message"] = rrr.error_message;
         j_response["error"] = j_error;
         response_str = j_response.dump();
         std::cout << "DEBUG:: here in handle_string_as_json(): B is the full dump:::::::::::::::\n" << response_str << std::endl;
      }
      return response_str;
   };

   // run return result wrapped in json.
   //
   auto handle_string_as_json = [handle_list_tools, handle_list_tools_with_search_pattern,
                                 make_response_from_functions, handle_execute_python_result] (const std::string &buf_str) {

      // std::cout << "handle this string: " << buf_str << ":" << std::endl;

      std::string response_str; // can return blank if we don't get a method that we  know

      int id = -1; // unset/unfound
      if (! buf_str.empty()) {
         json req = json::parse(buf_str);
         json::const_iterator j_id = req.find("id");
         if (j_id != req.end()) {
            id = j_id.value();
         } else {
            std::cout << "handle_string_as_json(): id not found - sad" << std::endl;
         }
         try {

            if (req["method"] == "python.exec") {
               std::string code = req["params"]["code"];
               execute_python_results_container_t rrr = execute_python_code_with_result_internal(code);
               response_str = handle_execute_python_result(rrr, id);
            }

            if (req["method"] == "python.exec_multiline") {
               std::string code = req["params"]["code"];
               execute_python_results_container_t rrr = execute_python_multiline_code_with_result_internal(code);
               response_str = handle_execute_python_result(rrr, id);
            }

            if (req["method"] == "mcp.list_tools") {
               int block_index = -1; // means all tools
               std::vector<std::pair<std::string, std::vector<func_doc> > > functions = handle_list_tools(block_index);
               bool add_python_exec = true;
               int n_max = 405;
               json response_funcs = make_response_from_functions(functions, id, n_max, add_python_exec);
               json j_response;
               j_response["jsonrpc"] = "2.0";
               j_response["id"] = std::to_string(id);
               j_response["result"] = response_funcs;
               response_str = j_response.dump();
            }

            if (req["method"] == "mcp.list_tools_block") {
               int block_index = 0;
               json::const_iterator it_params = req.find("params");
               if (it_params != req.end()) {
                  const json &j_params = *it_params;
                  // std::cout << "DEBUG:: ::::::::::::::: here in list tools block with j_params " << j_params << std::endl;
                  json::const_iterator it_block_index = j_params.find("block_index");
                  if (it_block_index != j_params.end()) {
                     block_index = it_block_index->get<int>();
                     std::cout << "DEBUG:: block_index: " << block_index << std::endl;
                  } else {
                     std::cout << "DEBUG:: block_index not found" << std::endl;
                  }
               }
               int n_max = 405;
               bool add_python_exec = true; // false maybe.
               std::vector<std::pair<std::string, std::vector<func_doc> > > functions = handle_list_tools(block_index);
               json response_funcs = make_response_from_functions(functions, id, n_max, add_python_exec);
               json j_response;
               j_response["jsonrpc"] = "2.0";
               j_response["id"] = std::to_string(id);
               j_response["result"] = response_funcs;
               response_str = j_response.dump();
            }

            if (req["method"] == "mcp.search") {
               json::const_iterator it_params = req.find("params");
               if (it_params != req.end()) {
                  const json &j_params = *it_params;
                  json::const_iterator it_pattern = j_params.find("pattern");
                  if (it_pattern != j_params.end()) {
                     std::string pattern = it_pattern->get<std::string>();
                     std::vector<std::pair<std::string, std::vector<func_doc> > > functions = handle_list_tools_with_search_pattern(pattern);
                     int n_max = 405;
                     bool add_python_exec = false;
                     json response_funcs = make_response_from_functions(functions, id, n_max, add_python_exec);
                     json j_response;
                     j_response["jsonrpc"] = "2.0";
                     j_response["id"] = std::to_string(id);
                     j_response["result"] = response_funcs;
                     response_str = j_response.dump();
                  } else {
                     std::cout << "DEBUG:: :::::: failed to find pattern in j_params" << std::endl;
                  }
               }
            }

         }
         catch (const std::exception &e) {
            std::cout << "WARNING:: coot_socket_listener_idle_func(): catch handle_string_as_json fail "
                      << buf_str << std::endl;
         }
      }
      return response_str;

   };

   auto write_all = [] (int fd, const std::string &s) {

      const char* p = static_cast<const char*>(s.c_str());
      size_t length = s.length();
      while (length > 0) {
         ssize_t n = write(fd, p, length);
         if (n < 0) {
            std::cout << "DEBUG:: write() error: " << errno << " (" << strerror(errno) << ")" << std::endl;
            switch(errno) {
               case (EINTR):
                  std::cout << "DEBUG:: EINTR" << std::endl;
                  // try again
                  break;
               case (EAGAIN):
                  std::cout << "DEBUG:: EAGAIN" << std::endl;
                  // output buffer was filled - write() needs to be called again to finish sending s
                  std::this_thread::sleep_for(std::chrono::microseconds(200));
                  break;
               default:
                  std::cout << "DEBUG:: neither EINTR not EAGAIN "
                            << EINTR << " " << EAGAIN << " errno " << errno << std::endl;
                  return false;  // real error
            }
         }
         if (n == 0) {
            return false;  // connection closed
         }
         if (n > 0) {
            p += n;
            length -= n;
         }
      }
      return true;
   };

   std::cout << "listening..." << std::endl;

   // If no active client connection, accept one non-blockingly
   if (client_fd < 0 && server_fd >= 0) {
      client_fd = accept(server_fd, nullptr, nullptr);
      if (client_fd >= 0) {
         // Set client socket non-blocking too
         int flags = fcntl(client_fd, F_GETFL, 0);
         fcntl(client_fd, F_SETFL, flags | O_NONBLOCK);
         std::cout << "Accepted new client socket connection.\n";
      }
   }
   // If we have a client, read any incoming data non-blockingly
   if (client_fd >= 0) {
      char buffer[4096];
      ssize_t n_read = read(client_fd, buffer, sizeof(buffer) - 1);
      // std::cout << "debug:: n_read: " << n_read << std::endl;
      if (n_read > 4) {
         int n_sent = int(buffer[3]) + 256 * int(buffer[2]) + 256 * 256 * int(buffer[1]) + 256 * 26 * 256 * int(buffer[0]);
         // std::cout << "debug:: n_sent: " << n_sent << std::endl;
         buffer[n_read] = '\0';
         std::cout << "<- Received: " << buffer+4 << std::endl;
         std::string buf_as_string(buffer+4, n_read);
         std::string r = handle_string_as_json(buf_as_string);

         const std::string &response_str = r;
         int32_t len = response_str.size();

         if (true) {
            std::cout << "-> Final-response:" << std::endl;
            std::cout << response_str << std::endl;
            std::cout << "DEBUG:: response_str len " << len << std::endl;
         }

         char header[4];
         header[0] = (len >> 24) & 0xFF;
         header[1] = (len >> 16) & 0xFF;
         header[2] = (len >> 8)  & 0xFF;
         header[3] = (len)       & 0xFF;

         std::string framed;
         framed.reserve(4 + response_str.size());
         framed.append(header, 4);        // 4-byte header
         framed.append(response_str);     // JSON body
         bool write_statue = write_all(client_fd, framed);
         // std::cout << "DEBUG:: write_status: " << write_statue << std::endl;

      } else if (n_read == 0) {
         // Client disconnected
         std::cout << "Client disconnected.\n";
         close(client_fd);
         client_fd = -1;
      } else if (n_read < 0 && errno != EWOULDBLOCK && errno != EAGAIN) {
         std::cerr << "Socket read error\n";
         close(client_fd);
         client_fd = -1;
      }
      // else: nothing to read right now
   }
   // Always return 1 (TRUE) to keep idle handler running
   return 1;

}

