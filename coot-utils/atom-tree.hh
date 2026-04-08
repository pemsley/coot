/*
 * coot-utils/atom-tree.hh
 *
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
//
#ifndef ATOM_TREE_HH
#define ATOM_TREE_HH

#include "geometry/protein-geometry.hh"
#include "atom-vertex.hh"
#include "map-index.hh"
#include "mini-mol/mini-mol.hh"

namespace coot {

  class atom_tree_t {

  protected:

      // either we create this class with a residue - in which case
      // residue is set and atom_selection is not,
      //
      // or we create this class with an atom selection, in which case
      // residue is not set.
      //
      // This difference is used in rotate_internal, where we decide
      // which atoms to move.
      //
      //
      mmdb::Residue *residue;
      mmdb::PPAtom atom_selection;   // for the multi-residue interface
      int     n_selected_atoms; // (we can't do residue_p->GetAtomTable())
      bool made_from_minimol_residue_flag;
      std::vector<std::pair<int, int> > bonds;
      std::vector<atom_vertex> atom_vertex_vec;
      void construct_internal(const dictionary_residue_restraints_t &rest,
                              mmdb::Residue *res,
                              const std::string &altconf);
      std::map<std::string, map_index_t, std::less<std::string> > name_to_index;
   private:
      bool fill_atom_vertex_vec(const dictionary_residue_restraints_t &rest, mmdb::Residue *res,
                                const std::string &altconf);
      bool fill_atom_vertex_vec_using_contacts(const std::vector<std::vector<int> > &contact_indices,
                                               int base_atom_index);
      bool fill_torsions(const dictionary_residue_restraints_t &rest, mmdb::Residue *res,
                         const std::string &altconf);
      void fill_name_map(const std::string &altconf);

      // return bool of 0 on not able to fill (not an exception).
      //
      std::pair<bool, atom_index_quad>
      get_atom_index_quad(const dict_torsion_restraint_t &tr,
                          mmdb::Residue *res, const std::string &altconf) const;

      std::vector<map_index_t> get_back_atoms(const map_index_t &index2) const;

      // don't add base index to forward atoms (of base_index).
      std::pair<int, std::vector<map_index_t> > get_forward_atoms(const map_index_t &base_index,
                                                                  const map_index_t &index2) const;
      std::vector<map_index_t>
      uniquify_atom_indices(const std::vector<map_index_t> &vin) const;

      std::vector<map_index_t>
      get_unique_moving_atom_indices(const std::string &atom1,
                                     const std::string &atom2,
                                     bool reversed_flag);


      // Return the complementary indices c.f. the moving atom set,
      // but do not include index2 or index3 in the returned set (they
      // do not move even with the reverse flag (of course)).
      std::vector<map_index_t>
      complementary_indices(const std::vector<map_index_t> &moving_atom_indices,
                            const map_index_t &index2,
                            const map_index_t &index3) const;

      // add forward_atom_index as a forward atom of this_index - but
      // only if forward_atom_index is not already a forward atom of
      // this_index.
      void add_unique_forward_atom(int this_index, int forward_atom_index);

      // so now we have a set of moving and non-moving atoms:
      //
      // Note: the angle is in radians.
      void rotate_internal(std::vector<map_index_t> moving_atom_indices,
                           const clipper::Coord_orth &dir,
                           const clipper::Coord_orth &base_atom_pos,
                           double angle);
      double quad_to_torsion(const map_index_t &index2) const;

      map_index_t get_index(mmdb::Atom *atom) const;

      bool in_forward_atoms(const map_index_t &bond_atom_index,
                            const map_index_t &fixed) const;


      // factored out:
      bool fill_atom_vertex_vec_using_contacts_by_atom_selection(const std::vector<std::vector<int> > &contact_indices,
                                                                 mmdb::PPAtom residue_atoms,
                                                                 int n_residue_atoms,
                                                                 int base_atom_index);


   public:

      // The angles are in degrees.
      //
      class tree_dihedral_info_t {
      public:
         atom_name_quad quad;
         double dihedral_angle;
         tree_dihedral_info_t(const atom_name_quad &quad_in, double ang_in) :
            quad(quad_in), dihedral_angle(ang_in) {}
         tree_dihedral_info_t() {}
         friend std::ostream& operator<<(std::ostream &o, tree_dihedral_info_t t);
      };

      class tree_dihedral_quad_info_t {
      public:
         atom_quad quad;
         double dihedral_angle;
         map_index_t fixed_atom_index; // from which the reverse flag gets derived;
         tree_dihedral_quad_info_t() {}
         // fixed can be unset, then we don't care if this is reversed or not.
         tree_dihedral_quad_info_t(const atom_quad &quad_in, double angle_in, const map_index_t &fixed) :
            quad(quad_in), fixed_atom_index(fixed) {
            dihedral_angle = angle_in;
         }
      };

      // the constructor throws an exception if there is no tree in
      // the restraints.
      //
      // FIXME.  In fact, we can handle there being no tree.  We
      // should fall back to using the bonds in the restraint. And
      // throw an exception if there are no bonds.
      //
      atom_tree_t(const dictionary_residue_restraints_t &rest, mmdb::Residue *res,
                  const std::string &altconf);

      // the constructor, given a list of bonds and a base atom index.
      // Used perhaps as the fallback when the above raises an
      // exception.
      //
      atom_tree_t(const std::vector<std::vector<int> > &contact_indices,
                  int base_atom_index,
                  mmdb::Residue *res,
                  const std::string &alconf);

      // constructor can throw an exception
      //
      // the constructor should not throw an exception if there are no
      // tree in the restraints.  It should instead try the bonds in
      // the restraints.  If there are no bonds then it throws an
      // exception.
      //
      atom_tree_t(const dictionary_residue_restraints_t &rest,
                  const minimol::residue &res_in,
                  const std::string &altconf);

      atom_tree_t(const dictionary_residue_restraints_t &rest,
                  const std::vector<std::vector<int> > &contact_indices,
                  int base_atom_index,
                  const minimol::residue &res_in,
                  const std::string &altconf);

      // constructor, given a list of bonds and a base atom index.
      //
      // Used for multi-residue tree generation.
      //
      atom_tree_t(const std::vector<std::vector<int> > &contact_indices,
                  int base_atom_index,
                  mmdb::Manager *mol,
                  int selection_handle);


      ~atom_tree_t() {
         if (made_from_minimol_residue_flag) {
            // question: does this delete the atoms in the residue too?
            delete residue;
            residue = 0;
         }
      }

      // when moving around the given atom pair, what are the fragment
      // sizes?  We ask this because generally, we want to flip the
      // small fragment.
      //
      std::pair<unsigned int, unsigned int>
      fragment_sizes(const std::string &atom1,
                     const std::string &atom2,
                     bool reversed_flag);

      // return a unique set
      std::vector<int>
      get_moving_atom_indices(const std::string &atom1,
                              const std::string &atom2,
                              bool reversed_flag);

      // Rotate round the 2 middle atoms of the torsion by angle (in
      // degress).  This is a relative rotation - not setting the
      // torsion angle.  The atoms of the mmdb::Residue residue are
      // manipulated.
      //
      // This can be used either by the residue or atom_selection
      // method of tree generation.
      //
      // If the reversed flag is true, this allows the rotation of the
      // "base" atoms, rather than the branch (dog wags rather than
      // tail).  atom1 and atom2 can be passed in either way round,
      // this function will sort out the position in the tree.  The
      // fragment rotation is reversed by setting the reversed_flag
      // (not the atom order).
      //
      // Return the new torsion angle (use the embedded torsion on
      // index2 if you can) Otherwise return -1.0;
      //
      // this can throw a std::runtime_error exception
      //
      double rotate_about(const std::string &atom1, const std::string &atom2,
                          double angle,
                          bool reversed_flag);


      // Return the new torsion angle (use the embedded torsion on
      // index2 if you can) Otherwise return -1.0;
      //
      // this can throw a std::runtime_error exception
      //
      // for use with multi-residue (atom-selection) (we can't use
      // atom names, because the may occur in different residues).
      //
      double rotate_about(int index_1, int index_2, double angle, bool reversed_flag);

      // this can throw an exception
      //
      // input angle in degrees
      //
      // return the dihedral angle (in degrees)
      //
      double set_dihedral(const std::string &atom1, const std::string &atom2,
                          const std::string &atom3, const std::string &atom4,
                          double angle);

      double set_dihedral(const map_index_t &i1,
                          const map_index_t &i2,
                          const map_index_t &i3,
                          const map_index_t &i4,
                          double angle);

      // reverse_flag to 1 to rotate "base" atoms rather than "tail" atoms.
      //
      double set_dihedral(const atom_quad &quad, double angle, bool reverse_flag);

      // this can throw an exception
      //
      // return the set of angles - should be the same that they were
      // set to (for validation).
      //
      // angle in degrees
      std::vector<double> set_dihedral_multi(const std::vector<tree_dihedral_info_t> &di);

      // this can throw an exception.
      //
      // angle in the double part of the pair, in degrees.
      //
      // return the set of angles - should be the same that they were
      // set to (for validation).

      std::vector<double> set_dihedral_multi(const std::vector<tree_dihedral_quad_info_t> &quads);


      minimol::residue GetResidue() const; // for use with above
                                           // function, where the
                                           // class constructor is
                                           // made with a minimol::residue.


   };
   std::ostream& operator<<(std::ostream &o, atom_tree_t::tree_dihedral_info_t t);

}


#endif // ATOM_TREE_HH
