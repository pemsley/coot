/* ligand/fast-ss-search.cc
 * 
 * Copyright 2008 The University of York
 * Author: Kevin Cowtan
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

#include "clipper/core/coords.h"
#include "clipper/core/xmap.h"
#include "mini-mol/mini-mol.hh"

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
     enum SSTYPE { ALPHA2, ALPHA3, ALPHA4, ALPHA3S,
		   BETA2,  BETA3,  BETA4,  BETA3S  };
     typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;

     class Target {
     public:
       Target( SSfind::SSTYPE type, int num_residues );
       const std::vector<Pair_coord>          target_coords() { return target_cs; }
       const std::vector<clipper::Coord_orth> calpha_coords() { return calpha_cs; }
     private:
       std::vector<Pair_coord>          target_cs;
       std::vector<clipper::Coord_orth> calpha_cs;
     };

     void prep_xmap( const clipper::Xmap<float>& xmap, const double radius );
     void prep_search( const clipper::Xmap<float>& xmap );
     void prep_search( const clipper::Xmap<float>& xmap, const double rhocut, const double radcut, const clipper::Coord_orth centre );
     std::vector<SearchResult> search( const std::vector<Pair_coord>& target_cs, const std::vector<clipper::RTop_orth>& ops, const double rhocut, const double frccut = 0.0 ) const;

   private:
     std::vector<float> mapbox;
     std::vector<int> srctrn;
     clipper::Grid grid;
     clipper::Grid_range mxgr;
     clipper::Mat33<> grrot;
   };


   class fast_secondary_structure_search : public SSfind {
   public:
     static int join_offset( const std::vector<clipper::Coord_orth>& frag1, const std::vector<clipper::Coord_orth>& frag2 );
     static double join_score( const std::vector<clipper::Coord_orth>& frag1, const std::vector<clipper::Coord_orth>& frag2 );
     static std::vector<std::vector<clipper::Coord_orth> > join( std::vector<std::vector<clipper::Coord_orth> > frags, int extn=0 );
     void operator()( const clipper::Xmap<float>& xmap,
		      const clipper::Coord_orth& centre,
		      double radius, std::vector<SSfind::Target> targets );
     minimol::molecule mol;
     bool success;
   };
}
