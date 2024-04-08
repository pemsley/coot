/* src/molecule-class-info-maps.cc
 * 
 * Copyright 2011 by The University of York
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <string>
#include <iomanip>
#include <fstream>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/Cartesian.h"
#include "coords/mmdb-crystal.h"
#include "molecule-class-info.h"

#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/cns/cns_map_io.h"



// export map fragment (.ext)
//
void
molecule_class_info_t::export_map_fragment(float radius,
					   clipper::Coord_orth comg,
					   const std::string &file_name) const {

   if (has_xmap()) {
      clipper::Cell cell = xmap.cell();
      clipper::Grid_sampling grid = xmap.grid_sampling();

      // get grid range
      // gr0: a grid range of the correct size (at the origin, going + and - in small box)
      // gr1: a grid range of the correct size (around the correct place, comg)
      clipper::Grid_range gr0(cell, grid, radius);
      clipper::Grid_range gr1(gr0.min() + comg.coord_frac(cell).coord_grid(grid),
			      gr0.max() + comg.coord_frac(cell).coord_grid(grid));

      // init nxmap
      clipper::NXmap<float> nxmap_local(cell, grid, gr1);
      clipper::Xmap<float>::Map_reference_coord ix(xmap);
      clipper::Coord_grid offset =
	 xmap.coord_map(nxmap_local.coord_orth(clipper::Coord_map(0.0,0.0,0.0))).coord_grid();
      typedef clipper::NXmap<float>::Map_reference_index NRI;
      for (NRI inx = nxmap_local.first(); !inx.last(); inx.next()) {
	 ix.set_coord(inx.coord() + offset);
	 nxmap_local[inx] = xmap[ix];
      }
      try { 
	 clipper::CCP4MAPfile mapout;
	 mapout.open_write(file_name);
	 mapout.set_cell(cell);
	 mapout.export_nxmap(nxmap_local);
	 mapout.close_write();
      }
      catch (const clipper::Message_base &except) {
	 std::cout << "WARNING:: Failed to write " << file_name  << std::endl;
      } 
   }
}


void
molecule_class_info_t::export_map_fragment_to_plain_file(float radius,
							 clipper::Coord_orth centre,
							 const std::string &filename) const {

   if (has_xmap()) {
      clipper::Cell cell = xmap.cell();
      clipper::Grid_sampling grid = xmap.grid_sampling();

      // get grid range
      // gr0: a grid range of the correct size (at the origin, going + and - in small box)
      // gr1: a grid range of the correct size (around the correct place, centre)
      clipper::Grid_range gr0(cell, grid, radius);
      clipper::Grid_range gr1(gr0.min() + centre.coord_frac(cell).coord_grid(grid),
			      gr0.max() + centre.coord_frac(cell).coord_grid(grid));

      // Example map metadata
      //
      // Number of columns, rows, sections ...............   19   19   17
      // Map mode ........................................    2
      // Start and stop points on columns, rows, sections    33   51   48   66   46   62
      // Grid sampling on x, y, z ........................  256  256  360
      // Cell dimensions .................................  142.2 142.2 200.41 90.0 90.0 120.0
      // Fast, medium, slow axes .........................    Y    X    Z
      // Minimum density .................................    -0.08186
      // Maximum density .................................     0.04684
      // Mean density ....................................    -0.02189
      // Rms deviation from mean density .................     0.02270
      // Space-group .....................................    1

      // init nxmap
      clipper::NXmap<float> nxmap_local(cell, grid, gr1);
      clipper::Xmap<float>::Map_reference_coord ix(xmap);
      clipper::Coord_grid offset =
	 xmap.coord_map(nxmap_local.coord_orth(clipper::Coord_map(0.0,0.0,0.0))).coord_grid();

      std::cout << " xmap cell:   " << cell.format() << std::endl;
      std::cout << " xmap grid_sampling: " << grid.format() << std::endl;
      std::cout << " xmap gr0:    " << gr0.format() << std::endl;
      std::cout << " xmap gr0.min:    " << gr0.min().format() << std::endl;
      std::cout << " xmap gr0.max:    " << gr0.max().format() << std::endl;
      std::cout << " xmap gr1:    " << gr1.format() << std::endl;
      std::cout << "nxmap Offset: " << offset.format() << std::endl; // start grid points

      std::ofstream f(filename.c_str());
      if (f) {
	 typedef clipper::NXmap<float>::Map_reference_index NRI;
	 int count = 0;
	 for (NRI inx = nxmap_local.first(); !inx.last(); inx.next()) {
	    ix.set_coord(inx.coord() + offset);
	    f << std::setw(12) << xmap[ix] << " ";
	    count++;
	    if (count%6==0) {
	       f << std::endl;
	       count = 0;
	    }
	 }
      }
      f << std::endl;
      f.close();
   }

}
