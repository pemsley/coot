/* src/molecule-class-info-maps.cc
 * 
 * Copyright 2011 by The University of York
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
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

// Having to set up the include files like this so that
// molecule-class-info.h can be parsed, is silly.

// For stat, mkdir:
#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <mmdb/mmdb_manager.h>
#include "mmdb-extras.h"
#include "Cartesian.h"
#include "mmdb-crystal.h"
#include "molecule-class-info.h"

#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/cns/cns_map_io.h"



// export map fragment (.ext)
//
void
molecule_class_info_t::export_map_fragment(float radius,
					   clipper::Coord_orth comg,
					   const std::string &file_name) const {

   if (has_map()) {
      clipper::Xmap<float> &xmap = xmap_list[0];
      clipper::Cell cell = xmap.cell();
      clipper::Grid_sampling grid = xmap.grid_sampling();

      // get grid range
      clipper::Grid_range gr0(cell, grid, radius);
      clipper::Grid_range gr1(gr0.min() + comg.coord_frac(cell).coord_grid(grid),
			      gr0.max() + comg.coord_frac(cell).coord_grid(grid));

      // init nxmap
      clipper::NXmap<float> nxmap(cell, grid, gr1);
      clipper::Xmap<float>::Map_reference_coord ix(xmap);
      clipper::Coord_grid offset =
	 xmap.coord_map(nxmap.coord_orth(clipper::Coord_map(0.0,0.0,0.0))).coord_grid();
      typedef clipper::NXmap<float>::Map_reference_index NRI;
      for (NRI inx = nxmap.first(); !inx.last(); inx.next()) {
	 ix.set_coord(inx.coord() + offset);
	 nxmap[inx] = xmap[ix];
      }
      clipper::CCP4MAPfile mapout;
      mapout.open_write(file_name);
      mapout.set_cell(cell);
      mapout.export_nxmap(nxmap);
      mapout.close_write();
   }
} 
