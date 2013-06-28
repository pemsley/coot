/* ligand/rigid-body.hh
 * 
 * Copyright 2005 by Paul Emsley, The University of York
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
 * 02110-1301, USA.
 */

#ifndef COOT_RIGID_BODY_HH
#define COOT_RIGID_BODY_HH

#include "mini-mol.hh"
#include "clipper/core/xmap.h"

// move the atoms of mol
namespace coot { 
   void rigid_body_fit(minimol::molecule *mol, const clipper::Xmap<float> &xmap);
   clipper::Vec3<double>
   get_rigid_body_angle_components(const std::vector<minimol::atom *> &atoms_p,
				   const clipper::Coord_orth &mean_pos,
				   const std::vector<clipper::Grad_orth<float> > &grad_vec,
				   double gradient_scale);
   void apply_angles_to_molecule(const clipper::Vec3<double> &angles,
				 const std::vector<minimol::atom *> *atoms_p,
				 const clipper::Coord_orth &mean_pos);
   float score_molecule(const coot::minimol::molecule &m,
			const clipper::Xmap<float> &xmap);

   std::pair<bool, clipper::RTop_orth>
   get_rigid_body_fit_rtop(minimol::molecule *mol,
			   const clipper::Xmap<float> &xmap);
   // as above but make a local RTop, that is, remove local_centre
   // from coordinates before calculating the RTop.
   std::pair<bool, clipper::RTop_orth>
   get_rigid_body_fit_rtop(minimol::molecule *mol,
			   const clipper::Coord_orth &local_centre,
			   const clipper::Xmap<float> &xmap);

}

#endif // COOT_RIGID_BODY_HH
