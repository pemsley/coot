/* ligand/ligand.hh
 * 
 * Copyright 2002, 2003, 2004, 2005, 2007 by The University of York
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
 * 02110-1301, USA.
 */

#include "rigid-body.hh"
#include "clipper/core/map_interp.h"

#include "coot-utils.hh"
#include "coot-map-utils.hh" // score_molecule


// move the atoms of m, don't use atoms with near zero occupancy for
// calculation of shifts.
//
void
coot::rigid_body_fit(coot::minimol::molecule *m, const clipper::Xmap<float> &xmap) {

   int round_max = 200;
   int iround = 0;
   double move_by_length = 1.0; // just to past test initially
   int n_atoms = 0;  // the number of non-zero occupancy atoms.

   double t = (xmap.cell().a()/xmap.grid_sampling().nu() +
	       xmap.cell().c()/xmap.grid_sampling().nv() +
	       xmap.cell().b()/xmap.grid_sampling().nw())/3.0;

   double gradient_scale = 0.3*t*t;  // somewhat arbitrary.

   while ( (iround < round_max) && (move_by_length > 0.0002) ) {
      clipper::Coord_orth midpoint(0,0,0);

      n_atoms = 0;
      for (unsigned int ifrag=0; ifrag<m->fragments.size(); ifrag++) {
	 for (unsigned int ires=(*m)[ifrag].min_res_no();
	      ires<(*m)[ifrag].max_residue_number(); ires++) {
	    for (unsigned int iat=0; iat<(*m)[ifrag][ires].atoms.size(); iat++) {
	       if (fabs((*m)[ifrag][ires][iat].occupancy) > 0.001) {
		  midpoint += (*m)[ifrag][ires][iat].pos;
		  n_atoms++;
	       }
	    }
	 }
      }
      if (n_atoms > 0) {

	 // assign the atom weights
	 std::vector<minimol::atom *> atoms = m->select_atoms_serial();
	 std::vector<double> atom_z_weights(atoms.size());
	 std::vector<std::pair<std::string, int> > atom_list = coot::util::atomic_number_atom_list();
	 for (int i=0; i<atoms.size(); i++) {
	    double z = coot::util::atomic_number(atoms[i]->element, atom_list);
	    if (z < 0.0) {
	       std::cout << "Unknown element :" << atoms[i]->element << ": " << std::endl;
	       z = 6.0; // as for carbon
	    } 
	    atom_z_weights[i] = z;
	 }
	 
	 double one_over = 1/double(n_atoms);
	 clipper::Coord_orth mean_pos = one_over * midpoint;

	 clipper::Grad_map<float> grad;
	 clipper::Grad_frac<float> grad_frac;
	 clipper::Grad_orth<float> grad_orth;
	 float dv, sum_dx = 0, sum_dy = 0, sum_dz = 0;
	 float sum_occ_fac = 0.0;
	 std::vector<clipper::Grad_orth<float> > grad_vec(atoms.size());
	 // only apply shifts for atoms with non-zero occupancy
	 for (int ii=0; ii<atoms.size(); ii++) {
	    if (fabs(atoms[ii]->occupancy) > 0.001) {
	       float occ_fac = atoms[ii]->occupancy;
	       if (occ_fac > 10.0)
		  occ_fac = 1.0;
	       else 
		  if (occ_fac < 0.0)
		  occ_fac = 1.0;
		     
	       clipper::Coord_orth atom_pos = atoms[ii]->pos;
	       clipper::Coord_frac atom_pos_frc = atom_pos.coord_frac(xmap.cell());
	       clipper::Coord_map  atom_pos_map = atom_pos_frc.coord_map(xmap.grid_sampling());
	       clipper::Interp_cubic::interp_grad(xmap, atom_pos_map, dv, grad);
	       grad_frac = grad.grad_frac(xmap.grid_sampling());
	       grad_orth = grad_frac.grad_orth(xmap.cell());
	       sum_dx += grad_orth.dx() * occ_fac * atom_z_weights[ii];
	       sum_dy += grad_orth.dy() * occ_fac * atom_z_weights[ii];
	       sum_dz += grad_orth.dz() * occ_fac * atom_z_weights[ii];
	       grad_vec[ii] = grad_orth;
	       sum_occ_fac += occ_fac;
	    } else {
	       grad_vec[ii] = clipper::Grad_orth<float>(0,0,0);
	    }
	 }
	 double datfrac = 1.0/double (sum_occ_fac);
	 clipper::Grad_orth<float> av_grad(sum_dx * datfrac,
					   sum_dy * datfrac,
					   sum_dz * datfrac);
	 
	 clipper::Coord_orth moved_by =
	    clipper::Coord_orth(gradient_scale*av_grad.dx(),
				gradient_scale*av_grad.dy(),
				gradient_scale*av_grad.dz());
	 clipper::Vec3<double> angles; 
	 if (n_atoms > 1)
	    angles = get_rigid_body_angle_components(atoms, mean_pos,
						     grad_vec, gradient_scale);
	 
	 move_by_length = sqrt(moved_by.lengthsq());
	 // std::cout << "move_by_length " << move_by_length << std::endl;
	 mean_pos += moved_by;
	 
	 for (unsigned int ifrag=0; ifrag<m->fragments.size(); ifrag++) {
	    for (unsigned int ires=(*m)[ifrag].min_res_no();
		 ires<=(*m)[ifrag].max_residue_number();
		 ires++) {
	       for (unsigned int iat=0; iat<(*m)[ifrag][ires].atoms.size(); iat++) {
		  (*m)[ifrag][ires][iat].pos += moved_by;
	       }
	    }
	 }
	 
	 if (n_atoms > 1) 
	    apply_angles_to_molecule(angles, &atoms, mean_pos); // modify atoms
	 
	 iround++;
      }
   }
}


clipper::Vec3<double>
coot::get_rigid_body_angle_components(const std::vector<minimol::atom *> &atoms,
				      const clipper::Coord_orth &mean_pos,
				      const std::vector<clipper::Grad_orth<float> > &grad_vec,
				      double gradient_scale) {


   // For each atom, we multiply the vector between the atom position
   // and the centre of rotation (V) by rotation_component to get the
   // vector pependicular to this vector in the appropriate plane
   // (e.g. XY plane for rotation round Z) - let's call that new vector
   // Vp.  We get the unit vector of Vp: Vpu.  We form the dot product
   // of the vector of gradients with Vpu.  This gives us a distance
   // (+ve or -ve) d_Vpu.
   //
   // We want an angle however, and that is atan2(d_Vpu,|V|).
   // Average that angle for contributions from all the atoms.
   // Do this for \alpha_x, \alpha_y, \alpha_z.

   // To get the angle contributions, we need to find the vector
   // perpendicular to the vector that joins the rotation centre to
   // the atom position.  And we need that, in the plane corresponding
   // to the axis about which we are rotating, (i.e. we need the
   // vector in the XY plane for rotations about z (i.e. the z
   // component should be zero).
   //
   clipper::RTop_orth rotation_component[3];
   rotation_component[0] = clipper::RTop_orth(clipper::Mat33<double>(0,0,0,0,0,1,0,-1,0),
					      clipper::Coord_frac(0,0,0));
   rotation_component[1] = clipper::RTop_orth(clipper::Mat33<double>(0,0,-1,0,0,0,1,0,0),
					      clipper::Coord_frac(0,0,0));
   rotation_component[2] = clipper::RTop_orth(clipper::Mat33<double>(0,1,0,-1,0,0,0,0,0),
					      clipper::Coord_frac(0,0,0));


   double sum_alpha[3];
   double sum_cos_theta[3];
   double sum_grad[3];
   double Vp_rms_sum[3];

   for (int ir=0; ir<3; ir++) { 
      sum_grad[ir] = 0.0;
      sum_alpha[ir] = 0.0;
      sum_cos_theta[ir] = 0.0; // debugging
      Vp_rms_sum[ir] = 0.0;
   }

   clipper::Coord_orth V, Vp[3];
   double dot_prod[3];
   int n_atoms = 0;
   for (unsigned int ii=0; ii<atoms.size(); ii++) {
      if (fabs(atoms[ii]->occupancy) > 0.001) { 
	 V = atoms[ii]->pos - mean_pos; // vector from centre to atom
	 clipper::Coord_orth grad(grad_vec[ii].dx(),
				  grad_vec[ii].dy(),
				  grad_vec[ii].dz());
	 for (int ir=0; ir<3; ir++) {  // rotation axis
	    Vp[ir] = V.transform(rotation_component[ir]);
	    Vp_rms_sum[ir] += Vp[ir].lengthsq();
	    dot_prod[ir] = clipper::Coord_orth::dot(grad,Vp[ir]);
	    sum_grad[ir] += dot_prod[ir];
	 }
	 n_atoms++; 
      }
   }

   double Vp_av_len[3];
   for (int ir=0; ir<3; ir++) // rotation axis
      Vp_av_len[ir] = 0.0;

   if (n_atoms > 0) { 
      for (int ir=0; ir<3; ir++) {  // rotation axis
	 Vp_av_len[ir] = sqrt(Vp_rms_sum[ir]/double(n_atoms));
      }
   }

   clipper::Vec3<double> a;
   for (int ir=0; ir<3; ir++) {  
      a[ir] = gradient_scale * 0.01 * sum_grad[ir]/Vp_av_len[ir];
      // Very useful diagnostic      
//       std::cout << "Vp_av_len[" << ir << "] is " << Vp_av_len[ir]; 
//       std::cout << "  a[" << ir << "] is " << a[ir]*57.3 << " degrees " << std::endl;
   }
   return a;
}


void 
coot::apply_angles_to_molecule(const clipper::Vec3<double> &angles ,
			       const std::vector<minimol::atom *> *atoms_p,
			       const clipper::Coord_orth &mean_pos) {


   double sin_t;
   double cos_t;

   sin_t = sin(-angles[0]);
   cos_t = cos(-angles[0]);
   clipper::Mat33<double> x_mat(1,0,0, 0,cos_t,-sin_t, 0,sin_t,cos_t);

   sin_t = sin(-angles[1]);
   cos_t = cos(-angles[1]);
   clipper::Mat33<double> y_mat(cos_t,0,sin_t, 0,1,0, -sin_t,0,cos_t);

   sin_t = sin(-angles[2]);
   cos_t = cos(-angles[2]);
   clipper::Mat33<double> z_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);

   clipper::Mat33<double> angle_mat = x_mat * y_mat * z_mat;
   clipper::RTop_orth angle_op(angle_mat, clipper::Coord_frac(0,0,0));

   // std::cout << angle_op.format() << std::endl;
   // std::cout << "mean pos in apply_angles: " << mean_pos.format() << std::endl; 
   // std::cout << (*atoms_p)[0]->pos.format() << std::endl;
   for (unsigned int ii=0; ii<atoms_p->size(); ii++) {
      (*atoms_p)[ii]->pos -= mean_pos;
      (*atoms_p)[ii]->pos = (*atoms_p)[ii]->pos.transform(angle_op);
      (*atoms_p)[ii]->pos += mean_pos;
   }
   // std::cout << (*atoms_p)[0]->pos.format() << std::endl;

}


float
coot::score_molecule(const coot::minimol::molecule &m,
		     const clipper::Xmap<float> &xmap) {

   float score = 0.0;
   for (unsigned int ifrag=0; ifrag<m.fragments.size(); ifrag++) {
      for (unsigned int ires=m[ifrag].min_res_no(); ires<m[ifrag].max_residue_number();
	   ires++) {
	 for (unsigned int iat=0; iat<m[ifrag][ires].atoms.size(); iat++) {
	    score += coot::util::density_at_point(xmap, m[ifrag][ires][iat].pos);
	 }
      }
   }
   return score; 
}
