/*
 * coot-utils/read-amber-trajectory.hh
 *
 * Copyright 2025 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef READ_AMBER_TRAJECTORY_HH
#define READ_AMBER_TRAJECTORY_HH

#include <string>
#include <vector>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   //! Read an Amber NetCDF trajectory file and create a multi-model molecule
   //!
   //! @param topology_mol The template molecule providing atom names, residues, connectivity
   //! @param trajectory_file_name Path to the Amber NetCDF trajectory file (.nc)
   //! @param start_frame First frame to read (0-indexed), -1 for first frame
   //! @param end_frame Last frame to read (0-indexed), -1 for last frame
   //! @param stride Read every nth frame (1 = all frames)
   //! @return A new mmdb::Manager with each frame as a separate model, or nullptr on failure
   mmdb::Manager *read_amber_trajectory(mmdb::Manager *topology_mol,
                                        const std::string &trajectory_file_name,
                                        int start_frame = -1,
                                        int end_frame = -1,
                                        int stride = 1);

   //! Get information about an Amber NetCDF trajectory file
   //!
   //! @param trajectory_file_name Path to the Amber NetCDF trajectory file
   //! @param n_frames Output: number of frames in trajectory
   //! @param n_atoms Output: number of atoms per frame
   //! @param has_box Output: whether periodic box information is present
   //! @return true if file was read successfully
   bool get_amber_trajectory_info(const std::string &trajectory_file_name,
                                  int &n_frames,
                                  int &n_atoms,
                                  bool &has_box);

} // namespace coot

#endif // READ_AMBER_TRAJECTORY_HH
