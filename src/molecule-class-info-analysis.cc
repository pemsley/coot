/*
 * src/molecule-class-info-analysis.cc
 *
 * Copyright 2020 by Medical Research Council
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


#include "molecule-class-info.h"


coot::model_composition_stats_t
molecule_class_info_t::get_model_composition_statistics() const {

   coot::model_composition_stats_t m;

   return m;

}

std::vector<std::string>
molecule_class_info_t::get_types_in_molecule() const {

   // 20250513-PE move this into coot-utils

   std::vector<std::string> v;
   std::set<std::string> s;
   mmdb::Manager *mol = atom_sel.mol;
   if (mol) {
      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
	 mmdb::Model *model_p = mol->GetModel(imod);
	 if (model_p) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(ichain);
	       int n_res = chain_p->GetNumberOfResidues();
	       for (int ires=0; ires<n_res; ires++) {
		  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		  if (residue_p) {
		     std::string type = residue_p->GetResName();
		     s.insert(type);
		  }
	       }
	    }
	 }
      }
   }
   for (const auto &item : s) {
      v.push_back(item);
   }
   return v;
}


#include "coot-utils/grid-balls.hh"
#include "density-contour/gaussian-surface.hh"
#include "colour-holder-to-glm.hh"

std::vector<coot::simple_mesh_t>
molecule_class_info_t::show_cavities(const coot::protein_geometry *geom_p) {

   std::vector<coot::simple_mesh_t> cavity_meshes;

   // geom.init_refmac_mon_lib("PC1.cif", read_number); need PC1 for testing.
   // but do it in python

   int imol = 0;
   coot::grid_balls_t gb(imol, atom_sel.mol, geom_p, 1.4, 4.5);
   gb.write_cavity_points("cavity-points.table");
   gb.write_subpocket_points("subpocket-points.table");

   // A Gaussian surface for each (non-trivial) cavity: build an mmdb molecule from
   // the cavity's grid points (create_mmdbmanager_from_points() uses chain "A"), then
   // contour it. Skip tiny cavities that can't make a meaningful surface.
   const std::string chain_id = "A";
   // 20260722-PE These have been tuned
   float sigma         = 1.2;
   float contour_level = 3.2;       // was 3
   float box_radius    = 3.0;
   float grid_scale    = 2.0;
   float b_factor      = 2.0;
   unsigned int min_points_for_mesh = 70;

   for (unsigned int ic=0; ic<gb.cavities.size(); ic++) {
      const coot::grid_balls_t::cavity_t &cav = gb.cavities[ic];
      if (cav.grid_indices.size() < min_points_for_mesh) continue;

      std::vector<clipper::Coord_orth> pts;
      pts.reserve(cav.grid_indices.size());
      for (unsigned int i=0; i<cav.grid_indices.size(); i++) {
         coot::grid_balls_t::point_3d_t p =
            gb.grid_point_to_mol_space(gb.deindex(cav.grid_indices[i]));
         pts.push_back(clipper::Coord_orth(p.x, p.y, p.z));
      }

      mmdb::Manager *cav_mol = coot::util::create_mmdbmanager_from_points(pts, b_factor);
      cav_mol->FinishStructEdit();
      coot::gaussian_surface_t gauss_surf(cav_mol, chain_id, sigma, contour_level,
                                          box_radius, grid_scale, b_factor);
      coot::simple_mesh_t surface_mesh = gauss_surf.get_surface();
      surface_mesh.name = "Cavity Surface " + std::to_string(ic);
      cavity_meshes.push_back(surface_mesh);
      delete cav_mol;
   }

   // give each cavity a distinct hue (the caller displays these meshes)
   for (unsigned int i=0; i<cavity_meshes.size(); i++) {
      coot::colour_holder ch(0.66, 0.44, 0.44, 0.5);
      ch.rotate_by(0.22 * static_cast<float>(i));
      glm::vec4 col = colour_holder_to_glm(ch);
      coot::simple_mesh_t &sm = cavity_meshes[i];
      for (unsigned int j=0; j<sm.vertices.size(); j++)
         sm.vertices[j].color = col;
   }

   return cavity_meshes;
}
