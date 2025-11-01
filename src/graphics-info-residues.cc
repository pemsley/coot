/* src/graphics-info-residues.cc
 *
 * Copyright 2011 by The University of Oxford.
 * Copyright 2015, 2016 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "compat/coot-sysdep.h"

#include <algorithm>
#include <iomanip>

#include "graphics-info.h"
#include "interface.h" // for create_multi_residue_torsion_dialog()
#include "ideal/torsion-bonds.hh"

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
      bool dict_status = geom_p->have_restraints_dictionary_for_residue_types(residue_types, imol,
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
	    gtk_widget_set_visible(w, TRUE);

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
graphics_info_t::setup_graphics_ligand_view_using_active_atom() {

   if (show_graphics_ligand_view_flag) { // user control
      std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
      if (active_atom.first) {
         coot::atom_spec_t spec(active_atom.second.second);
	 mmdb::Residue *residue_p = molecules[active_atom.second.first].get_residue(coot::residue_spec_t(spec));
	 if (false)
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
graphics_info_t::setup_graphics_ligand_view_using_active_atom(int only_in_imol) {

   if (show_graphics_ligand_view_flag) { // user control
      std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec(only_in_imol);
      if (active_atom.first) {
         coot::atom_spec_t spec(active_atom.second.second);
	 mmdb::Residue *residue_p = molecules[active_atom.second.first].get_residue(coot::residue_spec_t(spec));
	 if (false)
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
		  if (false)
		     std::cout << "debug:: setup_graphics_ligand() on residue "
			       << coot::residue_spec_t(residue_p) << std::endl;

                  gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));

                  graphics_ligand_view_imol = imol;

		  graphics_ligand_view_flag =
		     graphics_ligand_mesh_molecule.setup_from(imol, residue_p, alt_conf, Geom_p());

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
	       } else {
                  graphics_ligand_view_flag = false; // because no atoms
               }
	    }
	 }
      }
   }
}


#if 0 // delete this function when new version works - keep this for reference for now.
void
graphics_info_t::graphics_ligand_view() {

   std::cout << "old code graphics_ligand_view() " << std::endl;

#if 0
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
#endif
}

#endif

#include <glm/gtx/rotate_vector.hpp>
#include "matrix-utils.hh"

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
               g.setRotationCentre(coot::Cartesian(residue_centre.second));

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

            if (true) { // debug matrices
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

            coot::util::quaternion vq(md);
            glm::quat target_quat(0,0,0,1);
            // target_quat *= glm::rotate(target_quat, 1.57f, glm::vec3(0,1,0));

            // wrong
            // target_quat *= glm::conjugate(glm::quat(vq.q1, vq.q2, vq.q3, vq.q0));
            // target_quat *=   glm::inverse(glm::quat(vq.q1, vq.q2, vq.q3, vq.q0));
            // target_quat *= glm::quat(vq.q1, vq.q2, vq.q3, vq.q0);
            // target_quat *= glm::quat(vq.q0, vq.q1, vq.q2, vq.q3);
            // target_quat *= glm::inverse(glm::quat(vq.q0, vq.q1, vq.q2, vq.q3));
            // target_quat *= glm::conjugate(glm::quat(vq.q0, vq.q1, vq.q2, vq.q3));
            // target_quat = glm::conjugate(glm::quat(vq.q1, vq.q2, vq.q3, vq.q0)) * target_quat;
            // target_quat = glm::quat(vq.q1, vq.q2, vq.q3, vq.q0) * target_quat;
            // target_quat = glm::inverse(glm::quat(vq.q1, vq.q2, vq.q3, vq.q0)) * target_quat;
            // target_quat = glm::conjugate(glm::quat(vq.q0, vq.q1, vq.q2, vq.q3)) * target_quat;

            target_quat = glm::conjugate(glm::quat(vq.q1, vq.q2, vq.q3, vq.q0)) * target_quat;

            std::cout << "target_quat: " << glm::to_string(target_quat) << std::endl;

            // interpolate between current view and new view
            //

            float nice_zoom = 23.1; // "ligand-sized"
            const clipper::Coord_orth &rc = residue_centre.second;
            coot::Cartesian res_centre(rc.x(), rc.y(), rc.z());
            coot::Cartesian rot_centre = RotationCentre();
            coot::view_info_t view1(g.view_quaternion,  rot_centre,      zoom, "current");
            // target_quat in view2 has been wrongly calculated
            std::cout << "INFO:: perpendicular_ligand_view() calling interpolate with view2 position "
                      << res_centre << std::endl;
            coot::view_info_t view2(target_quat, res_centre, nice_zoom, "ligand-perp");
            int nsteps = go_to_ligand_animate_view_n_steps; // 20220303-PE No API to change this yet
            // std::cout << "perpendicular_ligand_view() calling interpolate()" << std::endl;
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


// static
void
graphics_info_t::add_terminal_residue_using_active_atom() {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();

   if (aa_spec_pair.first) {
      int imol = aa_spec_pair.second.first;
      coot::atom_spec_t &spec = aa_spec_pair.second.second;
      mmdb::Atom *at = molecules[imol].get_atom(spec);
      mmdb::Residue *residue_p = at->GetResidue();
      if (residue_p) {
         int atom_indx = -1;
         if (at->GetUDData(molecules[imol].atom_sel.UDDAtomIndexHandle, atom_indx) == mmdb::UDDATA_Ok) {
            std::string terminus_type = g.molecules[imol].get_term_type(atom_indx);
            std::string res_type = "ALA"; // elsewhere we are more clever than this.
            std::string chain_id = spec.chain_id;
            execute_add_terminal_residue(imol, terminus_type, residue_p, chain_id, res_type, true);
         }
      }
   }
}




// do it if have intermediate atoms and ctrl is pressed.
//
// axis: 0 for Z, 1 for X.
//
short int
graphics_info_t::rotate_intermediate_atoms_maybe(unsigned int width, unsigned int height) {

   // for rotation, use trackball rotation about the centre of the fragment.

   short int handled_flag = 0;

   if (moving_atoms_asc) {
      if (! last_restraints) { // we don't want this to happen during refinement
         if (moving_atoms_asc->n_selected_atoms > 0) {
            if (control_is_pressed) {
               handled_flag = true;
               bool sane_deltas = true;
               if (mouse_current_x - GetMouseBeginX() >  80) sane_deltas = false;
               if (mouse_current_x - GetMouseBeginX() < -80) sane_deltas = false;
               if (mouse_current_y - GetMouseBeginY() >  80) sane_deltas = false;
               if (mouse_current_y - GetMouseBeginY() < -80) sane_deltas = false;

               clipper::RTop_orth screen_matrix;
               clipper::RTop_orth screen_matrix_transpose;

#if 0
               {
                  coot::Cartesian centre = unproject_xyz(0, 0, 0.5);
                  coot::Cartesian front  = unproject_xyz(0, 0, 0.0);
                  coot::Cartesian right  = unproject_xyz(1, 0, 0.5);
                  coot::Cartesian top    = unproject_xyz(0, 1, 0.5);

                  coot::Cartesian screen_x = (right - centre);
                  // coot::Cartesian screen_y = (top   - centre);
                  coot::Cartesian screen_y = (centre - top);
                  coot::Cartesian screen_z = (front - centre);

                  screen_x.unit_vector_yourself();
                  screen_y.unit_vector_yourself();
                  screen_z.unit_vector_yourself();

                  float mat[16];
                  mat[ 0] = screen_x.x(); mat[ 1] = screen_x.y(); mat[ 2] = screen_x.z();
                  mat[ 4] = screen_y.x(); mat[ 5] = screen_y.y(); mat[ 6] = screen_y.z();
                  mat[ 8] = screen_z.x(); mat[ 9] = screen_z.y(); mat[10] = screen_z.z();

                  mat[ 3] = 0; mat[ 7] = 0; mat[11] = 0;
                  mat[12] = 0; mat[13] = 0; mat[14] = 0; mat[15] = 1;

                  screen_matrix.rot()(0,0) = screen_x.x();
                  screen_matrix.rot()(0,1) = screen_x.y();
                  screen_matrix.rot()(0,2) = screen_x.z();
                  screen_matrix.rot()(1,0) = screen_y.x();
                  screen_matrix.rot()(1,1) = screen_y.y();
                  screen_matrix.rot()(1,2) = screen_y.z();
                  screen_matrix.rot()(2,0) = screen_z.x();
                  screen_matrix.rot()(2,1) = screen_z.y();
                  screen_matrix.rot()(2,2) = screen_z.z();
                  screen_matrix.trn()[0] = 0.0;
                  screen_matrix.trn()[1] = 0.0;
                  screen_matrix.trn()[2] = 0.0;
                  screen_matrix_transpose.rot() = screen_matrix.rot().transpose();
               }

               if (sane_deltas) {
                  clipper::Coord_orth mac = moving_atoms_centre();
                  float spin_quat[4];
                  trackball(spin_quat,
                            (2.0*mouse_current_x - width)  /width,
                            (height -  2.0*mouse_current_y)/height,
                            (2.0*GetMouseBeginX() - width) /width,
                            (height - 2.0*GetMouseBeginY())/height,
                            get_trackball_size() );
                  GL_matrix m;
                  m.from_quaternion(spin_quat);
                  clipper::Mat33<double> m33 = m.to_clipper_mat();
                  clipper::Coord_orth zero(0,0,0);
                  for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
                     mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
                     clipper::Coord_orth origin_based_pos = coot::co(at) - mac;
                     clipper::RTop_orth rtop(m33, zero);
                     clipper::RTop_orth combined(screen_matrix_transpose * rtop * screen_matrix);
                     clipper::Coord_orth new_pos = origin_based_pos.transform(combined);
                     new_pos += mac;
                     at->x = new_pos.x(); at->y = new_pos.y(); at->z = new_pos.z();
                  }
                  Bond_lines_container bonds(*moving_atoms_asc, true);
                  regularize_object_bonds_box.clear_up();
                  regularize_object_bonds_box = bonds.make_graphical_bonds();
                  graphics_draw();
               }
#endif
            }
         }
      }
   }
   return handled_flag;
}


void
graphics_info_t::add_distance_labels_for_environment_distances() {

   bool background_is_black_flag = background_is_black_p();
   for (int i=0; i<environment_object_bonds_box.num_colours; i++) {
      float dark_bg_cor = 0.0;
      if (! background_is_black_flag)
         dark_bg_cor = 0.29;
      float fidx = static_cast<float>(i);
      coot::colour_holder col(0.96-dark_bg_cor, 0.96-0.4*fidx-dark_bg_cor, 0.5+0.5*fidx-dark_bg_cor);
      graphical_bonds_lines_list<graphics_line_t> &ll = environment_object_bonds_box.bonds_[i];
      for (int j=0; j< ll.num_lines; j++) {
         glm::vec3 s = cartesian_to_glm(ll.pair_list[j].positions.getStart());
         glm::vec3 e = cartesian_to_glm(ll.pair_list[j].positions.getFinish());
         float d = glm::distance(s,e);
         glm::vec3 mid_point = 0.5f * (s + e);
         glm::vec3 offset_mid_point = mid_point + glm::vec3(0.15, 0.05, 0.05);
         std::ostringstream ss;
         ss << std::setprecision(2) << std::fixed << d;
         // ss << std::fixed << d;
         std::string text(ss.str());

         glm::vec4 ccol(col.red, col.green, col.blue, 1.0);
         add_label(text, offset_mid_point, 0.8f * ccol);

      }
   }
}

