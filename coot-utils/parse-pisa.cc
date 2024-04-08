/*
 * coot-utils/parse-pisa.cc
 *
 * Copyright 2023 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include "pugixml.hpp"
#include <iostream>

#include "parse-pisa.hh"
#include "pugiconfig.hpp"

int parse_pisa(const std::string &file_name) {

   pugi::xml_document doc;
   pugi::xml_parse_result result = doc.load_file(file_name.c_str());
   if (!result)
      return -1;
        
   for (pugi::xml_node tool: doc.child("Profile").child("Tools").children("Tool")) {
      int timeout = tool.attribute("Timeout").as_int();
      if (timeout > 0)
         std::cout << "Tool " << tool.attribute("Filename").value() << " has timeout " << timeout << "\n";
   }
   return 0;
}


