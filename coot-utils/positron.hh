/*
 * coot-utils/positron.hh
 *
 * Copyright 2024 by Medical Research Council
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
#ifndef COOT_POSITRON_HH
#define COOT_POSITRON_HH

#include <string>
#include <vector>

namespace coot {
   class positron_metadata_t {
   public:
      float x;
      float y;
      std::vector<float> params;
      positron_metadata_t(float xx, float yy) : x(xx), y(yy) { params.resize(6,0); }
   };

   // should the following be part of a container class? If you like, yes.

   void read_positron_metadata(std::vector<positron_metadata_t> *data, const std::string &z_data, const std::string &s_data);

   // 20240323-PE Should we also pass a "must-be-closer-than-limit?"
   // Currently we use 0.1 in both x and y.
   // If we fail to find a close point, then return -1.
   int get_closest_positron_metadata_point(const std::vector<positron_metadata_t> &positron_metadata, float x, float y);

   class positron_metadata_container_t {
   public:
      positron_metadata_container_t() {}
      explicit positron_metadata_container_t(const std::vector<positron_metadata_t> &positron_metadata) : metadata(positron_metadata) {}
      std::vector<positron_metadata_t> metadata;
      // Currently we use limit of 0.1 in both x and y.
      // If we fail to find a close point, then return -1.
      int get_closest_positron_metadata_point(const std::pair<float, float> &z) const {
         return coot::get_closest_positron_metadata_point(metadata, z.first, z.second);
      }
      size_t size() const { return metadata.size(); }
   };

}

#endif // COOT_POSITRON_HH
