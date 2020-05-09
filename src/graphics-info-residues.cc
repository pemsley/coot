/* src/graphics-info-residues.cc
 * 
 * Copyright 2011 by The University of Oxford.
 * Copyright 2015, 2016 by Medical Research Council
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#include <algorithm>
#include <iomanip>

#include "graphics-info.h"
#include "interface.h" // for create_multi_residue_torsion_dialog()
#include "ideal/torsion-bonds.hh"

#include "trackball.h"
void
graphics_info_t::multi_torsion_residues(int imol, const std::vector<coot::residue_spec_t> &v) {

   if (is_valid_model_molecule(imol)) {

      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (! mol)
	 std::cout << "no (reference) mol" << std::endl;

      // we probably don't need residues, but we do need to check that
      // we have restraints for each of the residue types.
      // 
      std::vector<std::string> residue_types;
      for (unsigned int i=0; i<v.size(); i++) { 
	 mmdb::Residue *r= molecules[imol].get_residue(v[i]);
	 if (r) {
	    std::string comp_id = r->GetResName();
	    if (std::find(residue_types.begin(), residue_types.end(), comp_id) ==
		residue_types.end()) {
	       residue_types.push_back(comp_id);
	    }
	 } else {
	    std::cout << "WARNING:: residues " << v[i] << " not found in molecule "
		      << imol << std::endl;
	 }
      }

      // uses dynamic add
      bool dict_status = geom_p->have_dictionary_for_residue_types(residue_types, imol,
								   cif_dictionary_read_number);

      mmdb::Manager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);
      if (! moving_mol) {
	 std::cout << "WARNING:: multi_torsion_residues() no moving mol" << std::endl;
      } else {
	 imol_moving_atoms = imol;
	 // now select everything in moving mol
	 int selhnd = moving_mol->NewSelection();
	 mmdb::PPAtom atom_selection = 0;
	 int n_selected_atoms;
	 moving_mol->SelectAtoms(selhnd, 0, "*",
				 mmdb::ANY_RES, "*",
				 mmdb::ANY_RES, "*",
				 "*", // residue name
				 "*",
				 "*", 
				 "*"); // alt-loc
	 moving_mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

	 try { 
	    std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > pairs = 
	       coot::torsionable_bonds(imol, mol, atom_selection, n_selected_atoms, Geom_p());
	    
	    GtkWidget *w = wrapped_create_multi_residue_torsion_dialog(pairs);
	    gtk_widget_show(w);

	    moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
	    make_moving_atoms_graphics_object(imol, make_asc(moving_mol));
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: " << rte.what() << std::endl;
	 } 

	 moving_mol->DeleteSelection(selhnd);
      }
   }
} 


// bottom left flat ligand view:
// 
void
graphics_info_t::setup_graphics_ligand_view_aa() {

   if (show_graphics_ligand_view_flag) { // user control
      std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
      if (active_atom.first) {
	 mmdb::Residue *residue_p = molecules[active_atom.second.first].get_residue(active_atom.second.second);
	 if (0)
	    std::cout << "-------------------------- setup_graphics_ligand_view_aa() imol: "
		      << active_atom.second.first
		      << " residue: " << coot::residue_spec_t(residue_p) << std::endl;
	 setup_graphics_ligand_view(active_atom.second.first, residue_p, active_atom.second.second.alt_conf);
      }
   }
}

// bottom left flat ligand view:
// 
void
graphics_info_t::setup_graphics_ligand_view_aa(int imol) {

   if (show_graphics_ligand_view_flag) { // user control
      std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec(imol);
      if (active_atom.first) {
	 mmdb::Residue *residue_p = molecules[active_atom.second.first].get_residue(active_atom.second.second);
	 if (0)
	    std::cout << "-------------------------- setup_graphics_ligand_view_aa() imol: "
		      << active_atom.second.first
		      << " residue: " << coot::residue_spec_t(residue_p) << std::endl;
	 setup_graphics_ligand_view(active_atom.second.first, residue_p, active_atom.second.second.alt_conf);
      }
   }
}

void
graphics_info_t::setup_graphics_ligand_view(int imol, mmdb::Residue *residue_p, const std::string &alt_conf) {

   if (show_graphics_ligand_view_flag) { // user control
      if (!use_graphics_interface_flag) {
	 graphics_ligand_view_flag = false;
      } else { 
	 if (!residue_p) {
	    graphics_ligand_view_flag = false;
	 } else {
	    if (coot::util::residue_has_hetatms(residue_p) != 1) {
	       graphics_ligand_view_flag = false;
	    } else {
	       if (residue_p->GetNumberOfAtoms() > 1) { 
		  if (0)
		     std::cout << "   setup_graphics_ligand() on residue "
			       << coot::residue_spec_t(residue_p) << std::endl;
		  graphics_ligand_view_flag =
		     graphics_ligand_mol.setup_from(imol, residue_p, alt_conf, Geom_p(), background_is_black_p());

		  // This overwrites the atom info in the status bar - I don't like that.
		  // There needs to be a different mechanism to report the residue type.
		  // 
		  if (false) {
		     std::string res_name = residue_p->GetResName();
		     std::pair<bool, coot::dictionary_residue_restraints_t> p = 
			Geom_p()->get_monomer_restraints_at_least_minimal(res_name, imol);
		     if (! p.first) {
			// 
		     } else {
			const coot::dictionary_residue_restraints_t &restraints = p.second;
			add_status_bar_text(restraints.residue_info.name);
		     }
		  }
	       }
	    }
	 }
      }
   }
}



void
graphics_info_t::graphics_ligand_view() {

   if (! show_graphics_ligand_view_flag)  // user control
      return;

   if (! is_valid_model_molecule(graphics_ligand_mol.imol))
      return;

   if (! molecules[graphics_ligand_mol.imol].is_displayed_p())
      return;
   
   if (graphics_ligand_view_flag) { 

      try {
	 
	 graphics_info_t g;
         GtkAllocation allocation = get_glarea_allocation();

	 std::pair<lig_build::pos_t, lig_build::pos_t> ext = 
	    g.graphics_ligand_mol.ligand_extents();

	 float sc = 28;
	 float h = float(allocation.height);
	 float w = float(allocation.width);
	 float ar = h/w;
	 glPushMatrix();
	 glLoadIdentity();
	 glScalef(sc/w, sc/h, sc/100); // stops the view ligand
				       // changing size as the window
				       // is reshaped

	 // std::cout << "extents: top_left " << ext.first
	 // << "   bottom right" << ext.second << std::endl;
	 // 
	 // If the offset scale factor (now 0.9) is 1.0, then when we
	 // have big molecules, they sit too much towards the centre
	 // of the screen (i.e. the offset correction is too much).
	 // 
	 // glTranslatef(-20.5-0.8*ext.first.x, -21.0+0.8*ext.second.y, 0);
	 // glTranslatef(-20.5-0.8*ext.first.x, -20.5-0.8*ext.first.y, 0);

	 double screen_bottom_left_x_pos = -1 * w/sc;
	 double screen_bottom_left_y_pos = -1 * h/sc;
	 double x_trans = screen_bottom_left_x_pos -0.8*ext.first.x + 3;
	 double y_trans = screen_bottom_left_y_pos -0.8*ext.first.y + 2;
	 
	 glTranslatef(x_trans, y_trans, 0);

	 glMatrixMode(GL_PROJECTION);
	 glPushMatrix();
	 glLoadIdentity();
   
	 GLfloat col[3] = { 0.6, 0.6, 0.6 };
	 glColor3fv(col);

	 glEnable(GL_LINE_SMOOTH);
	 glEnable(GL_BLEND);
	 if (g.background_is_black_p())
	    glBlendFunc(GL_SRC_ALPHA,GL_ZERO);

	 if (0) {
	    GLfloat col[3] = { 0.9, 0.1, 0.1 };
	    glColor3fv(col);
	    double x_pos = ext.first.x;
	    double y_pos = ext.first.y;
	    glLineWidth(1.0);
	    glBegin(GL_LINES);
	    glVertex3d(x_pos-1, y_pos,   0);
	    glVertex3d(x_pos+1, y_pos,   0);
	    glVertex3d(x_pos,   y_pos-1, 0);
	    glVertex3d(x_pos,   y_pos+1, 0);
	    glEnd();
	 } 
	 
	 glLineWidth(2.0);
	 g.graphics_ligand_mol.render();

	 if (0) {
	    // box at corner of screen;
	    double x_pos = -1 * w/sc - x_trans;
	    double y_pos = -1 * h/sc - y_trans;
	    glColor3f(0.4, 0.7, 0.9);
	    glBegin(GL_POLYGON);
	    glVertex3d(x_pos-3, y_pos-3, 0);
	    glVertex3d(x_pos-3, y_pos+3, 0);
	    glVertex3d(x_pos+3, y_pos+3, 0);
	    glVertex3d(x_pos+3, y_pos-3, 0);
	    glEnd();
	 } 

	 // debug box, ligand space
	 if (0) { 
	    // the lines ------------------
	    glBegin(GL_LINES);
   
	    glVertex3f( 0.0,  0.0, 0.0);
	    glVertex3f( 1.0,  0.0, 0.0);

	    glVertex3f( 1.0,  0.0, 0.0);
	    glVertex3f( 1.0,  1.0, 0.0);
   
	    glVertex3f( 1.0,  1.0, 0.0);
	    glVertex3f( 0.0,  1.0, 0.0);
   
	    glVertex3f( 0.0,  1.0, 0.0);
	    glVertex3f( 0.0,  0.0, 0.0);
   
	    glEnd();
	    // end of the lines ------------
	 }

	 glDisable(GL_LINE_SMOOTH);
	 glDisable(GL_BLEND);
   
	 glPopMatrix();
	 glMatrixMode(GL_MODELVIEW);
	 glPopMatrix();
      }

      catch(const std::runtime_error &rte) {
	 // PE20130605
	 // This is useful for debugging, but not for production
	 // (e.g. where the comp-id is ZN we get here currently and it
	 // is printed on every frame)
	 // 
	 // std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }
}


void
graphics_info_t::perpendicular_ligand_view(int imol, const coot::residue_spec_t &residue_spec) {

   mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(residue_spec);
   if (residue_p) {

      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

      if (n_residue_atoms == 0) {
         std::cout << "ERROR: Residue " << residue_spec << "has 0 atoms " << std::endl;
      } else {
         clipper::Coord_orth running_centre(0.0, 0.0, 0.0);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            running_centre += coot::co(at);
         }
         double scale = 1/double(n_residue_atoms);
         clipper::Coord_orth mean_pos(running_centre.x() * scale,
                                      running_centre.y() * scale,
                                      running_centre.z() * scale);
         clipper::Matrix<double> mat(3,3);
         for (int ii=0; ii<3; ii++)
            for (int jj=0; jj<3; jj++)
               mat(ii,jj) = 0.0;

         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            clipper::Coord_orth co = coot::co(at);
            mat(0,0) += (co.x() - mean_pos.x()) * (co.x() - mean_pos.x());
            mat(0,1) += (co.x() - mean_pos.x()) * (co.y() - mean_pos.y());
            mat(0,2) += (co.x() - mean_pos.x()) * (co.z() - mean_pos.z());
            mat(1,0) += (co.y() - mean_pos.y()) * (co.x() - mean_pos.x());
            mat(1,1) += (co.y() - mean_pos.y()) * (co.y() - mean_pos.y());
            mat(1,2) += (co.y() - mean_pos.y()) * (co.z() - mean_pos.z());
            mat(2,0) += (co.z() - mean_pos.z()) * (co.x() - mean_pos.x());
            mat(2,1) += (co.z() - mean_pos.z()) * (co.y() - mean_pos.y());
            mat(2,2) += (co.z() - mean_pos.z()) * (co.z() - mean_pos.z());
         }

         graphics_info_t g;
         std::pair<bool, clipper::Coord_orth> residue_centre = g.molecules[imol].residue_centre(residue_p);
         if (residue_centre.first) {

            if (false) // new style ligand navigation 20160720
               g.setRotationCentre(residue_centre.second);

            if (false) { // debug matrices/eigen

               std::cout << "raw mat" << std::endl;
               for (int ii=0; ii<3; ii++) {
                  std::cout << "(";
                  for (int jj=0; jj<3; jj++) {
                     std::cout <<  "  " << mat(ii,jj);
                  }
                  std::cout << ")\n";
               }
            }

            std::vector<double> eigens = mat.eigen(false);

            if (false) { // debug matrices
               std::cout << "eigens:\n     " << sqrt(eigens[0]) << " " << sqrt(eigens[1]) << " "
                         << sqrt(eigens[2]) << std::endl;
               std::cout << "post mat" << std::endl;
               for (int ii=0; ii<3; ii++) {
                  std::cout << "(";
                  for (int jj=0; jj<3; jj++) {
                     std::cout << std::right << std::setprecision(4) << std::fixed <<  "   " << mat(ii,jj);
                  }
                  std::cout << ")\n";
               }
               std::cout.setf(std::ios::fixed, std::ios::floatfield);
            }

            clipper::Mat33<double> md(mat(0,0), mat(0,1), mat(0,2),
                                      mat(1,0), mat(1,1), mat(1,2),
                                      mat(2,0), mat(2,1), mat(2,2));

            // OK, so md gives us an unsorted eigenvector matrix,
            // which, when we apply it to the view orients the ligands
            // along one of screen X, Y or Z.  Ideally, we'd like the
            // long axis of the ligand to be along screen X (and the
            // thinnest direction of the ligand along screen Z).  It
            // might be possible to rotate md by rotating the columns
            // (or rows).  But I don't do that.  I now test for the
            // longest eigenvector and that gives us the axis about
            // which I need to rotate the view by 90 degrees to put
            // the long axis of the ligand along screen z (we don't
            // need a rotation if eigen[0] is the longest vector of
            // course (unless we consider the tweak to rotate the
            // second longest eigenvector along screen y).  The order
            // the rotations are applied is important.  The flags are
            // set so that screen-x rotation happens last.

            bool need_x_rotate = false;
            bool need_y_rotate = false;
            bool need_z_rotate = false;

            coot::util::quaternion vq(md);

            // convert a coot::util::quaternion to a float *.

            // ugly - we should have equivalent of this in the quaternion class.
            //
            float vqf[4];
            vqf[0] = vq.q0; vqf[1] = vq.q1; vqf[2] = vq.q2; vqf[3] = vq.q3;
            glm::quat target_quat(vq.q0, vq.q1, vq.q2, vq.q3);

            // test: initally along screen Y?
            if (eigens[1] > eigens[0]) {
               if (eigens[1] > eigens[2]) {
                  need_z_rotate = true;
                  if (eigens[2] > eigens[0]) // minor correction
                     need_x_rotate = true;
               }
            }

            // test: initally along screen Z?
            if (eigens[2] > eigens[0]) {
               if (eigens[2] > eigens[1]) {
                  need_y_rotate = true;
                  if (eigens[0] > eigens[1]) // minor correction
                     need_x_rotate = true;
               }
            }

            // test: initally along screen X? (only minor correction may be needed)
            if (eigens[0] > eigens[1])
               if (eigens[0] > eigens[2])
                  if (eigens[2] > eigens[1])
                     need_x_rotate = true;

            if (need_z_rotate) {
               float spin_quat[4];
               trackball(spin_quat, 1.0, 1.0, 1.0, 1.0 + 0.0174533*180, 0.4);
               add_quats(spin_quat, vqf, vqf);
            }
            if (need_y_rotate) {
               float spin_quat[4];
               trackball(spin_quat, 0, 0, 0.0174533*45, 0.000, 0.8);
               add_quats(spin_quat, vqf, vqf);
            }
            if (need_x_rotate) {
               float spin_quat[4];
               trackball(spin_quat, 0, 0, 0.0, 0.0174533*45, 0.8);
               add_quats(spin_quat, vqf, vqf);
            }

            // interpolate between current view and new view
            //

            float nice_zoom = 23.1; // "ligand-sized"
            const clipper::Coord_orth &rc = residue_centre.second;
            coot::Cartesian res_centre(rc.x(), rc.y(), rc.z());
            coot::Cartesian rot_centre = RotationCentre();
            coot::view_info_t view1(g.glm_quat,  rot_centre,      zoom, "current");
            coot::view_info_t view2(target_quat, res_centre, nice_zoom, "ligand-perp");
            int nsteps = 60;
            coot::view_info_t::interpolate(view1, view2, nsteps);

         }
      }
   }
}


void
graphics_info_t::eigen_flip_active_residue() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();

   if (aa_spec_pair.first) {
      int imol = aa_spec_pair.second.first;
      mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
      mmdb::Residue *residue_p = at->GetResidue();
      if (residue_p) {
         std::string chain_id = residue_p->GetChainID();
         int res_no = residue_p->GetSeqNum();
         molecules[imol].eigen_flip_residue(chain_id, res_no);
         graphics_draw();
      }
   }
}
