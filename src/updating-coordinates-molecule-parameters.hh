/*
 * src/updating-coordinates-molecule-parameters.hh
 *
 * Copyright 2019 by Medical Research Council
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

#ifndef UPDATING_COORDINATES_MOLECULE_PARAMETERS_T
#define UPDATING_COORDINATES_MOLECULE_PARAMETERS_T

#include<string>
#include <sys/types.h>// stat
#include <sys/stat.h>

// This class is for reading the output of Refmac
//
class updating_coordinates_molecule_parameters_t {
public:
   int imol;
   std::string pdb_file_name;
   timespec ctime;
   updating_coordinates_molecule_parameters_t() {
      imol = -1;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
   }
   explicit updating_coordinates_molecule_parameters_t(const std::string &file_name) : pdb_file_name(file_name) {
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
      imol = -1; // is this used?
   }
   void update_from_stat_info(const struct stat &s) {

#ifdef WINDOWS_MINGW
      ctime.tv_sec = s.st_ctime;
      ctime.tv_nsec = 0.; // not available!? Lets hope not necessary
#else
#ifndef _POSIX_SOURCE
      ctime = s.st_ctimespec; // Mac OS X?
#else
      ctime = s.st_ctim;
#endif
#endif
   }
};

// This class is for updating the difference map when the model has changed.
//
class updating_model_molecule_parameters_t {
public:
   int imol_coords;
   int imol_map_with_data_attached;
   int imol_2fofc_map;  // sigmaA weighted, that is
   int imol_fofc_map; // ditto
   updating_model_molecule_parameters_t() {
      imol_coords = -1;
      imol_map_with_data_attached = -1;
      imol_2fofc_map = -1;
      imol_fofc_map = -1;
   }
   updating_model_molecule_parameters_t(int imol_coords_in, int imol_d, int imol_map_in) :
      imol_coords(imol_coords_in), imol_map_with_data_attached(imol_d), imol_fofc_map(imol_map_in) { imol_2fofc_map = -1;
   }
   updating_model_molecule_parameters_t(int imol_coords_in, int imol_data, int imol_map_2fofc_in, int imol_map_fofc_in) :
      imol_coords(imol_coords_in), imol_map_with_data_attached(imol_data), imol_2fofc_map(imol_map_2fofc_in), imol_fofc_map(imol_map_fofc_in) {}
   std::string format() const {
      std::string ss;
      ss += "imol_coords: ";
      ss += std::to_string(imol_coords);
      ss += std::string(" ");

      ss += "imol_map_with_data_attached: ";
      ss += std::to_string(imol_map_with_data_attached);
      ss += std::string(" ");

      ss += "imol_2fofc_map: ";
      ss += std::to_string(imol_2fofc_map);
      ss += std::string(" ");

      ss += "imol_fofc_map: ";
      ss += std::to_string(imol_fofc_map);

      return ss;
   }
};

#endif
