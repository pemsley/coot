/* layla/inchi_key_database.cpp
 * 
 * Copyright 2025 by Global Phasing Ltd.
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

#include "inchi_key_database.hpp"
#include <string>
#include <sstream>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <utility>

coot::layla::InchiKeyDatabase coot::layla::parseInchikeyDatabase(std::istream& rawData) {
   InchiKeyDatabase inchiMap;
   std::string line;

   while (std::getline(rawData, line)) {
      if (line.empty() || line.starts_with('#'))
         continue;  // Skip empty or comment lines

      std::vector<std::string> segments;
      std::istringstream lineStream(line);
      std::string segment;

      while (std::getline(lineStream, segment, '\t')) {
         segments.push_back(segment);
      }

      if (segments.size() != 3) {
         throw std::runtime_error("Unexpected number of segments in InchiKey database line: " + line);
      }

      const auto& inchiKey = segments[0];
      const auto& monomer_id = segments[1];
      const auto& chem_name = segments[2];

      inchiMap[inchiKey] = {monomer_id, chem_name};
   }

   return inchiMap;
}
