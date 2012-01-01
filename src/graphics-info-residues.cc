/* src/graphics-info-residues.cc
 * 
 * Copyright 2011 by The University of Oxford.
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
 * Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <algorithm>

#include "graphics-info.h"
#include "interface.h" // for create_multi_residue_torsion_dialog()

void
graphics_info_t:: multi_torsion_residues(int imol, const std::vector<coot::residue_spec_t> &v) {

   if (is_valid_model_molecule(imol)) {

      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (! mol)
	 std::cout << "no (reference) mol" << std::endl;

      // we probably don't need residues, but we do need to check that
      // we have restraints for each of the residue types.
      // 
      std::vector<std::string> residue_types;
      for (unsigned int i=0; i<v.size(); i++) { 
	 CResidue *r= molecules[imol].get_residue(v[i]);
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
      bool dict_status = geom_p->have_dictionary_for_residue_types(residue_types);

      CMMDBManager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);
      if (! moving_mol) {
	 std::cout << "WARNING:: multi_torsion_residues() no moving mol" << std::endl;
      } else {
	 imol_moving_atoms = imol;
	 // now select everything in moving mol
	 int selhnd = moving_mol->NewSelection();
	 PPCAtom atom_selection = 0;
	 int n_selected_atoms;
	 moving_mol->SelectAtoms(selhnd, 0, "*",
				 ANY_RES, "*",
				 ANY_RES, "*",
				 "*", // residue name
				 "*",
				 "*", 
				 "*"); // alt-loc
	 moving_mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

	 try { 
	    std::vector<std::pair<CAtom *, CAtom *> > pairs = 
	       coot::torsionable_bonds(mol, atom_selection, n_selected_atoms, Geom_p());
	    
	    GtkWidget *w = wrapped_create_multi_residue_torsion_dialog(pairs);
	    gtk_widget_show(w);

	    moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
	    make_moving_atoms_graphics_object(make_asc(moving_mol));
	 }
	 catch (std::runtime_error rte) {
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

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
   if (active_atom.first) {
      CResidue *residue_p = molecules[active_atom.second.first].get_residue(active_atom.second.second);
      setup_graphics_ligand_view(residue_p);
   }
}

void
graphics_info_t::setup_graphics_ligand_view(CResidue *residue_p) {

   if (!residue_p) {
      graphics_ligand_view_flag = false;
   } else {
      if (coot::util::residue_has_hetatms(residue_p) != 1) {
	 std::cout << "   setup_graphics_ligand() clear for "
		   << coot::residue_spec_t(residue_p) << std::endl;
	 graphics_ligand_view_flag = false;
      } else { 
	 std::cout << "   setup_graphics_ligand() on residue "
		   << coot::residue_spec_t(residue_p) << std::endl;
	 graphics_ligand_view_flag =
	    graphics_ligand_mol.setup_from(residue_p, Geom_p(), background_is_black_p());
      }
   }
}



void
graphics_info_t::graphics_ligand_view() {

   if (graphics_ligand_view_flag) { 

      try {
	 
	 graphics_info_t g;
	 std::pair<lig_build::pos_t, lig_build::pos_t> ext = 
	    g.graphics_ligand_mol.ligand_extents();

	 // float sc = 0.1;
	 float sc = 0.043;
	 glPushMatrix();
	 glLoadIdentity();
	 glScalef(sc, sc, sc);
	 //       glTranslatef(-8, -8., 0); // unscaled molecule

	 // std::cout << "extents: top_left " << ext.first
	 // << "   bottom right" << ext.second << std::endl;
	 // 
	 // If the offset scale factor (now 0.9) is 1.0, then when we
	 // have big molecules, they sit too much towards the centre
	 // of the screen (i.e. the offset correction is too much).
	 // 
	 // glTranslatef(-20.5-0.8*ext.first.x, -21.0+0.8*ext.second.y, 0);
	 glTranslatef(-20.5-0.8*ext.first.x, -20.5-0.8*ext.first.y, 0);
	 
	 glMatrixMode(GL_PROJECTION);
	 glPushMatrix();
	 glLoadIdentity();
   
	 GLfloat col[3] = { 0.6, 0.6, 0.6 };
	 glColor3fv(col);
	 glLineWidth(2.0);

	 glEnable(GL_LINE_SMOOTH);
	 glEnable(GL_BLEND);
	 if (g.background_is_black_p())
	    glBlendFunc(GL_SRC_ALPHA,GL_ZERO);

	 g.graphics_ligand_mol.render();

	 // debug box
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

      catch (std::runtime_error rte) {
	 std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }
} 

