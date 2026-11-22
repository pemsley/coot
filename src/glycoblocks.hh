/*
 * src/glycoblocks.hh
 *
 * 3D-SNFG "glycoblocks" for carbohydrates, ported from Moorhen (wasm_src/glycoblocks.cc).
 * See glycoblocks.cc for the copyright and licence.
 */

#ifndef GLYCOBLOCKS_HH
#define GLYCOBLOCKS_HH

#include <string>
#include <mmdb2/mmdb_manager.h>
#include "api/instancing.hh"

namespace coot {
   namespace glycoblocks {
      //! the glycoblocks (SNFG shapes and the linkages between them) for the sugars in the
      //! residues of the given selection. three_d: the 3D-SNFG style - cubes and spheres (and
      //! thick prisms) rather than the flat blocks of Moorhen
      coot::instanced_mesh_t DrawSugarBlocks(mmdb::Manager *molHnd, const std::string &cid_str, bool three_d = false,
                                             float scale_factor = 1.0f,
                                             const glm::vec4 &block_colour = glm::vec4(-1.0f));
   }
}

#endif // GLYCOBLOCKS_HH
