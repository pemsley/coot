
#ifndef REBOX_MAP_HH
#define REBOX_MAP_HH

#include <clipper/core/xmap.h>
#include <clipper/core/coords.h>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   namespace util {

      struct reboxed_map_t {
         clipper::Xmap<float> xmap;
         clipper::Coord_orth offset; // translation applied to the map origin
         bool success;
         std::string message;
         reboxed_map_t() : offset(0,0,0), success(false) {}
      };

      //! Rebox a map to a cubic box around the atoms in the selection.
      //!
      //! The new map is a P1 cubic Xmap with orthogonal axes whose grid spacing
      //! matches the original map. The cube edge length is the maximum extent
      //! of the atom selection plus border, expanded to fit n_pixels_per_edge
      //! grid points if that value is positive.
      //!
      //! @param xmap_in        the input map
      //! @param mol            the molecule
      //! @param SelectionHandle the mmdb atom selection handle
      //! @param border         padding in Angstroms beyond the atom extents
      //! @param n_pixels_per_edge number of grid points along each edge (0 = auto)
      //! @return a reboxed_map_t containing the new map and the translation offset
      reboxed_map_t rebox_map(const clipper::Xmap<float> &xmap_in,
                              mmdb::Manager *mol,
                              int SelectionHandle,
                              float border,
                              int n_pixels_per_edge);

   }
}

#endif // REBOX_MAP_HH
