/* ligand/backrub-rotamer.hh
 * 
 * Copyright 2009 by The University of Oxford.
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

#ifndef BACKRUB_ROTAMER_HH
#define BACKRUB_ROTAMER_HH

#include <string>

#include <clipper/core/xmap.h>
#include "mini-mol/mini-mol.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "coot-utils/atom-selection-container.hh"

namespace coot {

   // not just for backrubbing.
   // 20201030-PE now return the clashing waters too.
        std::pair<float, std::vector<mmdb::Atom *> > get_clash_score(const minimol::molecule &a_rotamer,
                         atom_selection_container_t asc, int water_interaction_mode);

   class backrub {
      mmdb::Residue *orig_this_residue;
      mmdb::Residue *orig_prev_residue;
      mmdb::Residue *orig_next_residue;
      std::string alt_conf;
      clipper::Coord_orth ca_prev;
      clipper::Coord_orth ca_next;
      clipper::Coord_orth ca_this;
      // thow an exception on failure to find CA of prev or next.
      void setup_this_and_prev_next_ca_positions();
      minimol::fragment make_test_fragment(mmdb::Residue *r, double rotation_angle) const;
      std::string chain_id;
      float score_fragment(minimol::fragment &frag) const;
      // clipper::Xmap<float> xmap;
      const clipper::Xmap<float> *xmap_p;
      mmdb::Manager *stored_mol;
      minimol::residue
      make_residue_include_only(mmdb::Residue *orig_prev_residue,
				const std::vector<std::string> &prev_res_atoms) const;
      clipper::Coord_orth rotamer_residue_centre() const;
      float residue_radius(const clipper::Coord_orth &rc);

      // do a check of the residue numbers and chaid id so that "same
      // residue" clashes are not counted.
      // 20201030-PE now return the clashing waters too.
      std::pair<float, std::vector<mmdb::Atom*> >
      get_clash_score(const minimol::molecule &a_rotamer,
		                mmdb::PPAtom sphere_atoms, int n_sphere_atoms,
                     int water_interaction_mode) const;

      void rotate_individual_peptide(mmdb::Residue *r, double rotation_angle,
				     minimol::fragment *f) const;
      // an analysis function
      // 
      double sample_individual_peptide(mmdb::Residue *r, double rotation_angle,
				       minimol::fragment *f,
				       mmdb::Residue *residue_front,
				       mmdb::Residue *residue_back,
				       bool is_leading_peptide_flag) const;
      // A modelling function.
      // Fiddle with the atom positions of fragment f.
      // 
      void rotate_individual_peptides_back_best(mmdb::Residue *r, double rotation_angle,
						minimol::fragment *f) const;

      // fiddle with the peptide position in f
      // 
      void apply_back_rotation(minimol::fragment *f,				   
			       bool is_leading_peptide_flag,
			       double best_back_rotation_angle) const;

      std::vector<mmdb::Atom *> clashing_waters; // these are turned into waters_for_deletion

   public:

      // Throw an exception on failure to construct the backrub internals.
      // Throw an exception if this, previous or next residues are null
      backrub(const std::string &chain_id_in,
	      mmdb::Residue *this_r,
	      mmdb::Residue *prev_r,
	      mmdb::Residue *next_r,
	      const std::string &alt_conf_in,
	      mmdb::Manager *mol_in,
	      const clipper::Xmap<float> *xmap_in_p) : alt_conf(alt_conf_in) {
            orig_this_residue = this_r;
            orig_prev_residue = prev_r;
            orig_next_residue = next_r;
            setup_this_and_prev_next_ca_positions();
            chain_id = chain_id_in;
            xmap_p = xmap_in_p;
            stored_mol = mol_in;
         }

      // throw an exception on failure to get a good search result.
      std::pair<coot::minimol::molecule, float> search(const dictionary_residue_restraints_t &rest);

      // maybe we need to delete a water or two to get a good fit for the side chain?
      std::vector<atom_spec_t> waters_for_deletion() const;
   };

   void backrub_molecule(mmdb::Manager *mol, const clipper::Xmap<float> *xmap_p,
                        const coot::protein_geometry &pg);

}


#endif // BACKRUB_ROTAMER_HH
