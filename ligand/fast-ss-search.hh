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

   /*! Fast template search class.
     Suitable for interactive secondary structure searches. */
   class SSfind {
   public:
     enum SSTYPE { ALPHA2, ALPHA3, ALPHA4, BETA2, BETA3, BETA4 };

     //! constructor: takes search angle step (degrees or radians)
     SSfind( const double step = 15.0 );
     //! perform secondary structure search using specified SS target
     std::vector<std::vector<clipper::Coord_orth> > operator() ( const clipper::Xmap<float>& xmap, int num_residues, SSTYPE type );
     //! perform fast template search using arbitrary target
     std::vector<clipper::RTop_orth> operator() ( const clipper::Xmap<float>& xmap, std::vector<std::pair<clipper::Coord_orth,clipper::Coord_orth> >& all_co );
     //! return list of scores to go with the results from the search functions
     const std::vector<float>& scores() { return scores_; };
   private:
     double step_;
     std::vector<float> scores_;
   };


   class fast_secondary_structure_search : public SSfind {
   public:
     fast_secondary_structure_search( const clipper::Xmap<float> &xmap,
				      int num_residues, SSTYPE type );
     minimol::molecule mol;
     bool success;
   };
}
