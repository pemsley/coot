/* coot-utils/read-cube.hh
 *
 * Copyright 2026 by Medical Research Council
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

#ifndef COOT_UTILS_READ_CUBE_HH
#define COOT_UTILS_READ_CUBE_HH

#include <string>
#include <vector>

#include <clipper/core/coords.h>

namespace mmdb { class Manager; }
namespace clipper { template<class T> class Xmap; }

namespace coot {
   namespace util {

      // A single atom as read from a Gaussian/ORCA "cube" file. The atomic
      // number and position come straight from the file; position has been
      // converted to Angstroms.
      //
      class cube_atom_t {
      public:
         int atomic_number;
         double nuclear_charge;          // the second column in the atom line
         clipper::Coord_orth position;   // Angstroms
         cube_atom_t(int z, double c, const clipper::Coord_orth &p) :
            atomic_number(z), nuclear_charge(c), position(p) {}
      };

      // The parsed contents of a Gaussian/ORCA cube file.
      //
      // The grid is described by an origin and three axis step vectors (as in
      // the file). For a standard ORCA/Gaussian cube the step vectors are
      // axis-aligned and mutually orthogonal. All lengths are in Angstroms
      // (cube files are usually in Bohr; the reader converts).
      //
      class cube_t {
      public:
         bool read_success;
         std::string title;                    // the two comment lines, joined
         clipper::Coord_orth origin;           // Angstroms
         clipper::Coord_orth axis_step[3];     // per-axis voxel step vector, Angstroms
         int n_points[3];                       // number of voxels along each axis
         bool is_mo_cube;                       // true if this holds molecular orbital data
         std::vector<int> mo_indices;           // orbital indices (mo cube only)
         std::vector<cube_atom_t> atoms;
         // Volumetric data, fastest-varying index last (axis 2, then the MO
         // index for an mo cube). Size = n_points[0]*n_points[1]*n_points[2]
         // (* mo_indices.size() for an mo cube). May be empty if only the
         // header and atoms were wanted.
         std::vector<float> data;

         cube_t() : read_success(false), is_mo_cube(false) {
            n_points[0] = 0; n_points[1] = 0; n_points[2] = 0;
         }

         // Build a (new) mmdb::Manager holding the cube's atoms, as a single
         // residue. The caller owns the returned pointer. Returns nullptr if
         // there are no atoms.
         mmdb::Manager *make_mmdb_manager() const;

         // Fill a clipper Xmap from the grid data (requires that the data was
         // read). The grid is wrapped in a P1 cell sized (n * step) with the
         // origin carried as an integer grid offset (the EM-map precedent), so
         // the map registers with the atoms in true coordinates. Returns false
         // if there is no grid data or the grid is not axis-aligned/orthogonal
         // (which would need a general cell or an NXmap). For an mo cube only
         // the first orbital is used.
         bool make_xmap(clipper::Xmap<float> *xmap_p) const;

         // True if the grid data has a real negative lobe (orbital-like), so it
         // should be displayed as a difference map (+/- contours).
         bool data_has_negative_lobe() const;
      };

      // Read a Gaussian/ORCA cube file. On failure the returned cube_t has
      // read_success == false.
      //
      // If read_grid_data is false, only the header and atoms are parsed (the
      // data vector is left empty) - useful when only the model is wanted.
      //
      cube_t read_cube_file(const std::string &file_name, bool read_grid_data=true);

   }
}

#endif // COOT_UTILS_READ_CUBE_HH
