/* ideal/regularize-minimol.hh
 * 
 * Copyright 2004, 2005 The University of York
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

#ifndef REGULARIZE_MINIMOL_HH
#define REGULARIZE_MINIMOL_HH

#ifdef HAVE_GSL

#include "mini-mol/mini-mol.hh"
#include "simple-restraint.hh"


namespace coot {

   minimol::molecule
   regularize_minimol_molecule(const minimol::molecule &molin,
			       const protein_geometry &geom);

   // update the coordinates of frag_in internal to this function
   void refine_minimol_fragment(minimol::fragment &frag_in,
                                const protein_geometry &geom,
                                const clipper::Xmap<float> &xmap,
                                float weight = 60.0,
                                bool do_GM = false,
                                float GM_distance_lim = 4.0f,
                                float GM_alpha = 0.02f);

}

#endif // HAVE_GSL
#endif // REGULARIZE_MINIMOL_HH

