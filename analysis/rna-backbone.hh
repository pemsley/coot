/*
 * analysis/rna-backbone.hh
 *
 * Copyright 2023 by Medical Research Council
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

#ifndef COOT_UTILS_RNA_BACKBONE_HH
#define COOT_UTILS_RNA_BACKBONE_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include "analysis/stats.hh"

namespace coot {

   class rna_backbone_t {
   public:
      stats::single base_samples;
      stats::single backbone_samples;

      rna_backbone_t(mmdb::Manager *mol, const clipper::Xmap<float> &xmap);
      void scan_all();
   };

}


#endif // COOT_UTILS_RNA_BACKBONE_HH
