/* ideal/torsion-bonds.hh
 * 
 * Copyright 2015 by Medical Research Council
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

// Created: 20170721-PE
// 
// There should be more functions here.

#ifndef IDEAL_TORSION_BONDS_HH
#define IDEAL_TORSION_BONDS_HH

namespace coot {

   // this can throw an exception
   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
   torsionable_bonds(int imol,
		     mmdb::Manager *mol, mmdb::PPAtom atom_selection, int n_selected_atoms,
		     protein_geometry *geom);
   // not sure this needs to public
   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
   torsionable_link_bonds(std::vector<mmdb::Residue *> residues_in,
			  mmdb::Manager *mol,
			  protein_geometry *geom);

   // And the atom_quad version of that (for setting link torsions)
   // 
   std::vector<torsion_atom_quad>
   torsionable_quads(int imol,
		     mmdb::Manager *mol, mmdb::PPAtom atom_selection, int n_selected_atoms,
		     protein_geometry *geom);
   std::vector<torsion_atom_quad>
   torsionable_link_quads(int imol,
			  std::vector<mmdb::Residue *> residues_in,
			  mmdb::Manager *mol, protein_geometry *geom_p);
   
   // this can throw an std::runtime exception
   void multi_residue_torsion_fit_map(int imol,
				      mmdb::Manager *mol,
				      const clipper::Xmap<float> &xmap,
				      const std::vector<std::pair<bool, clipper::Coord_orth> > &avoid_these_atoms, // flag is is-water?b
				      int n_trials,
				      protein_geometry *geom_p); 
   // which calls 
   double get_rand_angle(double current_angle, const torsion_atom_quad &quad, int itrial,
			 int n_trials,
			 bool allow_conformer_switch,
			 bool small_torsion_changes);
   
   // Does this model bang into itself?
   // Don't give atoms that are both in a quad a bang score
   // 
   double get_self_clash_score(mmdb::Manager *mol,
			       mmdb::PPAtom atom_selection,
			       int n_selected_atoms,
			       const std::vector<torsion_atom_quad> &quads);

   double get_environment_clash_score(mmdb::Manager *mol,
				      mmdb::PPAtom atom_selection,
				      int n_selected_atoms,
				      const std::vector<std::pair<bool, clipper::Coord_orth> > &avoid_these_atoms); // flag is is-water?

   bool both_in_a_torsion_p(mmdb::Atom *at_1,
			    mmdb::Atom *at_2,
			    const std::vector<torsion_atom_quad> &quads);
} 

#endif // IDEAL_TORSION_BONDS_HH
