/* layla/similar_ligands.hpp
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
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

#ifndef LAYLA_SIMILAR_LIGANDS_HPP
#define LAYLA_SIMILAR_LIGANDS_HPP

#include <string>
#include <gtk/gtk.h>

#include "ligand_editor_canvas.hpp" // CootLigandEditorCanvas

namespace coot { class protein_geometry; }
namespace RDKit { class ROMol; }

namespace coot {
   namespace layla {

      // Resolve the shipped COD-types similarity database file:
      // $COOT_COD_TYPES_DB else <package-data-dir>/cod-types-db-normalized.bin
      std::string cod_types_db_path();

      // Type the query molecule by its acedrg/COD atom types, search the database
      // for monomers with Jaccard similarity > min_jaccard, and present them in a
      // window of lazily-rendered 2D thumbnails. geom_p is used to fetch/draw the
      // matched monomers (local monomer library, else GitHub monomers, else PDBe
      // CCD). parent is the transient-for window. canvas + mol_idx identify the
      // sketched molecule, so that "Assign atom names from this ligand" can write
      // names back onto it.
      void search_for_similar_ligands(GtkWindow *parent,
                                      CootLigandEditorCanvas *canvas,
                                      unsigned int mol_idx,
                                      protein_geometry *geom_p,
                                      const RDKit::ROMol &query_mol,
                                      double min_score = 0.35);
   }
}

#endif // LAYLA_SIMILAR_LIGANDS_HPP
