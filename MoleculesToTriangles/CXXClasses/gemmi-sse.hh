/*
 * MoleculesToTriangles/CXXClasses/gemmi-sse.hh
 *
 * gemmi-native secondary-structure assignment (replaces mmdb
 * model->CalcSecStructure / residue->SSE) using gemmi's DSSP (dssp.hpp).
 *
 * The per-residue result is stamped onto gemmi::Residue::flag:
 *   'H' = helix, 'E' = strand, 'L' = loop. Non-polypeptide / non-polymer
 *   residues are left untouched ('\0', treated as "no secondary structure").
 *
 * (Future: a lighter port of mmdb's CalcSecStructure could replace the DSSP
 *  call here without changing this interface.)
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef MOLECULESTOTRIANGLES_GEMMI_SSE_HH
#define MOLECULESTOTRIANGLES_GEMMI_SSE_HH

#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      // Run DSSP over model 0 and stamp Residue::flag with 'H'/'E'/'L'.
      void assign_secondary_structure(gemmi::Structure &st);

      // Convenience accessors (read the stamped flag).
      inline bool residue_is_helix (const gemmi::Residue &res) { return res.flag == 'H'; }
      inline bool residue_is_strand(const gemmi::Residue &res) { return res.flag == 'E'; }
   }
}

#endif // MOLECULESTOTRIANGLES_GEMMI_SSE_HH
