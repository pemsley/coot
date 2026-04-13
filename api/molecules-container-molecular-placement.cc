/*
 * api/molecules-container-molecular-replacement.cc
 *
 * Copyright 2025 by Medical Research Council
 * Author: Paul Emsley
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
#include <iomanip>
#include <sstream>

#include "molecules-container.hh"
#include "ligand/molecular-replacement.hh"
#include "coot-utils/coot-map-utils.hh"

std::vector<molecules_container_t::mr_solution_t>
molecules_container_t::molecular_placement_fit(int imol_map, int imol_model,
                                               float x, float y, float z,
                                               int n_rotation_solutions,
                                               int n_translation_solutions) {

   std::vector<mr_solution_t> all_solutions;

   if (! is_valid_map_molecule(imol_map)) {
      std::cout << "ERROR:: " << __FUNCTION__ << " invalid map molecule " << imol_map << std::endl;
      return all_solutions;
   }
   if (! is_valid_model_molecule(imol_model)) {
      std::cout << "ERROR:: " << __FUNCTION__ << " invalid model molecule " << imol_model << std::endl;
      return all_solutions;
   }

   clipper::Coord_orth target_centre(x, y, z);
   // Note: do not hold a reference to molecules[imol_map].xmap across push_back
   // calls — vector reallocation would invalidate it. Copy the xmap pointer index
   // and re-fetch as needed.
   mmdb::Manager *mol_model = molecules[imol_model].atom_sel.mol;

   // --- Run the core MR search (in ligand/) ---

   std::vector<coot::mr_solution_t> core_solutions =
      coot::molecular_replacement_search(molecules[imol_map].xmap, mol_model, target_centre,
                                         n_rotation_solutions, n_translation_solutions);

   // --- Wrap each core result into an API-level solution ---

   for (unsigned int i=0; i<core_solutions.size(); i++) {

      auto &cs = core_solutions[i];
      if (! cs.placed_mol) continue;

      unsigned int sol_idx = i;

      int imol_placed = molecules.size();
      atom_selection_container_t asc = make_asc(cs.placed_mol);
      std::string name = "MR-solution-" + std::to_string(sol_idx);
      coot::molecule_t new_mol(asc, imol_placed, name);
      molecules.push_back(new_mol);

      // Compute mean density at CA (protein) and C1' (nucleic acid) positions
      float mean_density_ca = 0.0f;
      {
         int sel = cs.placed_mol->NewSelection();
         cs.placed_mol->SelectAtoms(sel, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                                    "*", " CA , C1'", "*", "*");
         mmdb::PPAtom atoms = nullptr;
         int n_atoms = 0;
         cs.placed_mol->GetSelIndex(sel, atoms, n_atoms);
         if (n_atoms > 0) {
            float sum = 0.0f;
            for (int j=0; j<n_atoms; j++) {
               clipper::Coord_orth pt(atoms[j]->x, atoms[j]->y, atoms[j]->z);
               sum += coot::util::density_at_point(molecules[imol_map].xmap, pt);
            }
            mean_density_ca = sum / static_cast<float>(n_atoms);
         }
         cs.placed_mol->DeleteSelection(sel);
      }

      std::ostringstream pdb_name;
      pdb_name << "mr-solution-" << std::setfill('0') << std::setw(3) << sol_idx << ".pdb";
      write_coordinates(imol_placed, pdb_name.str());

      mr_solution_t sol;
      sol.rotation = cs.rotation;
      sol.rotation_score = cs.rotation_score;
      sol.translation = cs.translation;
      sol.translation_score = cs.translation_score;
      sol.mean_density_at_ca = mean_density_ca;
      sol.imol = imol_placed;
      sol.pdb_filename = pdb_name.str();
      all_solutions.push_back(sol);

      // make_asc stores the pointer (does not deep copy), so placed_mol
      // is now owned by the molecule_t in molecules[]. Do not delete it.
   }

   // --- Summary table ---

   std::cout << "INFO:: MR solution summary:" << std::endl;
   std::cout << "INFO:: " << std::setw(5) << "Index"
             << std::setw(10) << "RF_score"
             << std::setw(10) << "TF_sigma"
             << std::setw(12) << "CA/C1'_dens"
             << "   " << "PDB file" << std::endl;
   std::cout << "INFO:: " << std::string(65, '-') << std::endl;
   for (unsigned int i=0; i<all_solutions.size(); i++) {
      const auto &s = all_solutions[i];
      std::cout << "INFO:: " << std::setw(5) << i
                << std::fixed << std::setprecision(2)
                << std::setw(10) << s.rotation_score
                << std::setprecision(1)
                << std::setw(10) << s.translation_score
                << std::setprecision(3)
                << std::setw(12) << s.mean_density_at_ca
                << "   " << s.pdb_filename << std::endl;
   }

   return all_solutions;
}
