/* layla/python-utils.hpp
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

#ifndef LAYLA_PYTHON_UTILS_HPP
#define LAYLA_PYTHON_UTILS_HPP
#include <string>
#include "Python.h"

namespace coot::layla {

      std::string get_drug_via_wikipedia_and_drugbank_curl(const std::string &drugname);
#if 0
      void setup_python_basic(int argc, char **argv);
      void setup_python_coot_module();
      void setup_python_module(const std::string &module_name);
      PyObject *safe_python_command_with_return(const std::string &python_cmd);
      std::string get_drug_via_wikipedia_and_drugbank_py(const std::string &drugname);
#endif

}


#endif // LAYLA_PYTHON_UTILS_HPP
