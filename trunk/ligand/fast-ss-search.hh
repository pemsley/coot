/* ligand/fast-ss-search.cc
 * 
 * Copyright 2008 The University of York
 * Author: Kevin Cowtan
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#include "clipper/core/coords.h"
#include "clipper/core/xmap.h"
#include "mini-mol.hh"

namespace coot {


   /*! Result class */
   class SearchResult {
   public:
     float score; int rot; int trn;
     bool operator <( SearchResult other ) const { return score < other.score; }
   };


   /*! Fast template search class.
     Suitable for interactive secondary structure searches. */
   class SSfind {
   public:
     enum SSTYPE { ALPHA2, ALPHA3, ALPHA4, BETA2, BETA3, BETA4 };
     typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;

     void prep_xmap( const clipper::Xmap<float>& xmap, const double radius, const double rhocut = 0.0 );
     void prep_target( SSfind::SSTYPE type, int num_residues );
     void set_target_score( const double score );
     const std::vector<SearchResult>& search( const std::vector<clipper::RTop_orth>& ops, const double rhocut, const double frccut = 0.0 );
   
     const std::vector<Pair_coord>          target_coords() { return target_cs; }
     const std::vector<clipper::Coord_orth> calpha_coords() { return calpha_cs; }
     const std::vector<SearchResult>& results() const { return rslts; }
     const float score_cutoff() const { return bestcut; }

   private:
     std::vector<Pair_coord>          target_cs;
     std::vector<clipper::Coord_orth> calpha_cs;
     std::vector<float> mapbox;
     float bestcut;
     clipper::Grid grid;
     clipper::Grid_range mxgr;
     clipper::Mat33<> grrot;
     std::vector<SearchResult> rslts;
   };


   class fast_secondary_structure_search : public SSfind {
   public:
     void operator()( const clipper::Xmap<float>& xmap,
		      const clipper::Coord_orth& centre,
		      double radius, int num_residues, SSTYPE type );
     minimol::molecule mol;
     bool success;
   };
}
