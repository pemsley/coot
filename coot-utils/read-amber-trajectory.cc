/*
 * coot-utils/read-amber-trajectory.cc
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

#include <iostream>
#include <vector>
#include <cstring>

#include "read-amber-trajectory.hh"

#ifdef HAVE_NETCDF
#include <netcdf.h>
#endif

namespace coot {

#ifdef HAVE_NETCDF

   // Helper function to check NetCDF errors
   static bool check_nc_error(int status, const std::string &operation) {
      if (status != NC_NOERR) {
         std::cout << "ERROR:: NetCDF " << operation << ": " << nc_strerror(status) << std::endl;
         return false;
      }
      return true;
   }

   bool get_amber_trajectory_info(const std::string &trajectory_file_name,
                                  int &n_frames,
                                  int &n_atoms,
                                  bool &has_box) {
      n_frames = 0;
      n_atoms = 0;
      has_box = false;

      int ncid;
      int status = nc_open(trajectory_file_name.c_str(), NC_NOWRITE, &ncid);
      if (!check_nc_error(status, "opening file " + trajectory_file_name))
         return false;

      // Get dimensions
      int frame_dim_id, atom_dim_id;
      size_t frame_len, atom_len;

      status = nc_inq_dimid(ncid, "frame", &frame_dim_id);
      if (!check_nc_error(status, "getting frame dimension")) {
         nc_close(ncid);
         return false;
      }
      nc_inq_dimlen(ncid, frame_dim_id, &frame_len);
      n_frames = static_cast<int>(frame_len);

      status = nc_inq_dimid(ncid, "atom", &atom_dim_id);
      if (!check_nc_error(status, "getting atom dimension")) {
         nc_close(ncid);
         return false;
      }
      nc_inq_dimlen(ncid, atom_dim_id, &atom_len);
      n_atoms = static_cast<int>(atom_len);

      // Check for box information
      int cell_lengths_id;
      status = nc_inq_varid(ncid, "cell_lengths", &cell_lengths_id);
      has_box = (status == NC_NOERR);

      nc_close(ncid);
      return true;
   }

   mmdb::Manager *read_amber_trajectory(mmdb::Manager *topology_mol,
                                        const std::string &trajectory_file_name,
                                        int start_frame,
                                        int end_frame,
                                        int stride) {

      std::cout << "INFO:: read_amber_trajectory: starting..." << std::endl;

      if (!topology_mol) {
         std::cout << "ERROR:: read_amber_trajectory: null topology molecule" << std::endl;
         return nullptr;
      }

      if (stride < 1) stride = 1;

      // Get trajectory info
      std::cout << "INFO:: read_amber_trajectory: getting trajectory info..." << std::endl;
      int n_frames_total, n_atoms_traj;
      bool has_box;
      if (!get_amber_trajectory_info(trajectory_file_name, n_frames_total, n_atoms_traj, has_box)) {
         return nullptr;
      }
      std::cout << "INFO:: read_amber_trajectory: trajectory has " << n_frames_total
                << " frames, " << n_atoms_traj << " atoms" << std::endl;

      // Count atoms in topology
      int n_atoms_topology = 0;
      for (int imod = 1; imod <= topology_mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model = topology_mol->GetModel(imod);
         if (model) {
            int n_chains = model->GetNumberOfChains();
            for (int ic = 0; ic < n_chains; ic++) {
               mmdb::Chain *chain = model->GetChain(ic);
               if (chain) {
                  int n_res = chain->GetNumberOfResidues();
                  for (int ir = 0; ir < n_res; ir++) {
                     mmdb::Residue *res = chain->GetResidue(ir);
                     if (res) {
                        n_atoms_topology += res->GetNumberOfAtoms();
                     }
                  }
               }
            }
         }
         break; // Only count first model
      }
      std::cout << "INFO:: read_amber_trajectory: topology has " << n_atoms_topology << " atoms" << std::endl;

      if (n_atoms_topology != n_atoms_traj) {
         std::cout << "ERROR:: read_amber_trajectory: atom count mismatch - "
                   << "topology has " << n_atoms_topology << " atoms, "
                   << "trajectory has " << n_atoms_traj << " atoms" << std::endl;
         return nullptr;
      }

      // Adjust frame range
      if (start_frame < 0) start_frame = 0;
      if (end_frame < 0 || end_frame >= n_frames_total) end_frame = n_frames_total - 1;
      if (start_frame > end_frame) {
         std::cout << "ERROR:: read_amber_trajectory: invalid frame range" << std::endl;
         return nullptr;
      }

      int n_frames_to_load = (end_frame - start_frame) / stride + 1;
      std::cout << "INFO:: read_amber_trajectory: will load " << n_frames_to_load
                << " frames (from " << start_frame << " to " << end_frame
                << ", stride " << stride << ")" << std::endl;

      // Open NetCDF file
      std::cout << "INFO:: read_amber_trajectory: opening NetCDF file..." << std::endl;
      int ncid;
      int status = nc_open(trajectory_file_name.c_str(), NC_NOWRITE, &ncid);
      if (!check_nc_error(status, "opening file"))
         return nullptr;

      // Get coordinates variable
      int coords_var_id;
      status = nc_inq_varid(ncid, "coordinates", &coords_var_id);
      if (!check_nc_error(status, "getting coordinates variable")) {
         nc_close(ncid);
         return nullptr;
      }

      // Use the original topology molecule's model as template
      mmdb::Model *template_model = topology_mol->GetModel(1);
      if (!template_model) {
         std::cout << "ERROR:: read_amber_trajectory: no model in topology" << std::endl;
         nc_close(ncid);
         return nullptr;
      }

      // Get chain count for iteration
      int n_chains = template_model->GetNumberOfChains();

      // Create new molecule for trajectory (empty - we'll add models)
      std::cout << "INFO:: read_amber_trajectory: creating molecule..." << std::endl;
      mmdb::Manager *traj_mol = new mmdb::Manager();

      // Buffer for one frame of coordinates
      std::vector<float> coords(n_atoms_traj * 3);

      std::cout << "INFO:: read_amber_trajectory: loading frames..." << std::endl;

      int model_number = 1;
      int frames_loaded = 0;
      for (int iframe = start_frame; iframe <= end_frame; iframe += stride) {
         // Progress output every 50 frames
         if (frames_loaded % 50 == 0) {
            std::cout << "INFO:: read_amber_trajectory: loading frame " << iframe
                      << " (" << frames_loaded << "/" << n_frames_to_load << ")" << std::endl;
         }

         // Read coordinates for this frame
         size_t start[3] = {static_cast<size_t>(iframe), 0, 0};
         size_t count[3] = {1, static_cast<size_t>(n_atoms_traj), 3};

         status = nc_get_vara_float(ncid, coords_var_id, start, count, coords.data());
         if (!check_nc_error(status, "reading coordinates for frame " + std::to_string(iframe))) {
            continue; // Skip this frame but try others
         }

         // Create a new model for this frame
         mmdb::Model *new_model = new mmdb::Model();
         new_model->SetMMDBManager(traj_mol, model_number);

         // Copy structure from template and update coordinates
         // atom_idx tracks position in the coordinate array (same order as template_atoms)
         int atom_idx = 0;
         for (int ic = 0; ic < n_chains; ic++) {
            mmdb::Chain *template_chain = template_model->GetChain(ic);
            if (!template_chain) continue;

            mmdb::Chain *new_chain = new mmdb::Chain();
            new_chain->SetChainID(template_chain->GetChainID());
            new_model->AddChain(new_chain);

            int n_res = template_chain->GetNumberOfResidues();
            for (int ir = 0; ir < n_res; ir++) {
               mmdb::Residue *template_res = template_chain->GetResidue(ir);
               if (!template_res) continue;

               mmdb::Residue *new_res = new mmdb::Residue();
               new_res->SetResID(template_res->GetResName(),
                                 template_res->GetSeqNum(),
                                 template_res->GetInsCode());
               new_chain->AddResidue(new_res);

               int n_at = template_res->GetNumberOfAtoms();
               for (int ia = 0; ia < n_at; ia++) {
                  mmdb::Atom *template_atom = template_res->GetAtom(ia);
                  if (!template_atom) continue;

                  mmdb::Atom *new_atom = new mmdb::Atom();
                  new_atom->Copy(template_atom);

                  // Use running atom index - atoms are in same order as trajectory
                  if (atom_idx < n_atoms_traj) {
                     new_atom->x = coords[atom_idx * 3 + 0];
                     new_atom->y = coords[atom_idx * 3 + 1];
                     new_atom->z = coords[atom_idx * 3 + 2];
                  }
                  atom_idx++;

                  new_res->AddAtom(new_atom);
               }
            }
         }

         traj_mol->AddModel(new_model);
         model_number++;
         frames_loaded++;
      }

      std::cout << "INFO:: read_amber_trajectory: finished reading frames, finalizing..." << std::endl;

      nc_close(ncid);

      // Finalize
      std::cout << "INFO:: read_amber_trajectory: PDBCleanup..." << std::endl;
      traj_mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL | mmdb::PDBCLEAN_INDEX);
      std::cout << "INFO:: read_amber_trajectory: FinishStructEdit..." << std::endl;
      traj_mol->FinishStructEdit();

      std::cout << "INFO:: read_amber_trajectory: loaded " << (model_number - 1)
                << " frames from " << trajectory_file_name << std::endl;

      return traj_mol;
   }

#else // HAVE_NETCDF not defined

   bool get_amber_trajectory_info(const std::string &trajectory_file_name,
                                  int &n_frames,
                                  int &n_atoms,
                                  bool &has_box) {
      std::cout << "WARNING:: get_amber_trajectory_info: NetCDF support not compiled" << std::endl;
      n_frames = 0;
      n_atoms = 0;
      has_box = false;
      return false;
   }

   mmdb::Manager *read_amber_trajectory(mmdb::Manager *topology_mol,
                                        const std::string &trajectory_file_name,
                                        int start_frame,
                                        int end_frame,
                                        int stride) {
      std::cout << "WARNING:: read_amber_trajectory: NetCDF support not compiled" << std::endl;
      return nullptr;
   }

#endif // HAVE_NETCDF

} // namespace coot
