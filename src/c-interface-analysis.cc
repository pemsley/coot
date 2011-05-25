/* src/c-interface.h
 * 
 * Copyright 2011 by The University of Oxford
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include "graphics-info.h"
#include "cc-interface.hh"
#include "c-interface.h"

#include "coot-hole.hh"

/* ------------------------------------------------------------------------- */
/*                      HOLE                                                 */
/* ------------------------------------------------------------------------- */
/*! \name Coot's Hole implementation */
void hole(int imol, float start_x, float start_y, float start_z, float end_x, float end_y, float end_z,
	  float colour_map_multiplier, float colour_map_offset,
	  int n_runs, bool show_probe_radius_graph_flag) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      CMMDBManager *mol = g.molecules[imol].atom_sel.mol;
      clipper::Coord_orth p_1(start_x, start_y, start_z);
      clipper::Coord_orth p_2(  end_x,   end_y,   end_z);

      coot::hole hole(mol, p_1, p_2, *g.Geom_p());
      hole.set_colour_shift(colour_map_multiplier, colour_map_offset);
      std::pair<std::vector<std::pair<clipper::Coord_orth, double> >, std::vector<coot::hole_surface_point_t> >
	 hole_path_and_surface = hole.generate();

      int obj_path    = new_generic_object_number("Probe path");
      int obj_surface = new_generic_object_number("Probe surface");
   
      for (unsigned int i=0; i<hole_path_and_surface.first.size(); i++) {
	 to_generic_object_add_point(obj_path, "red", 3,
				     hole_path_and_surface.first[i].first.x(),
				     hole_path_and_surface.first[i].first.y(),
				     hole_path_and_surface.first[i].first.z());
      }

      for (unsigned int i=0; i<hole_path_and_surface.second.size(); i++) { 
	 to_generic_object_add_point(obj_surface,
				     hole_path_and_surface.second[i].colour.hex().c_str(),
				     1, // pixel
				     hole_path_and_surface.second[i].position.x(),
				     hole_path_and_surface.second[i].position.y(),
				     hole_path_and_surface.second[i].position.z());
      }

      set_display_generic_object(obj_path,    1);
      set_display_generic_object(obj_surface, 1);

   }

} 
