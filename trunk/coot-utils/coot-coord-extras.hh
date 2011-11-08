/* coot-utils/coot-coord-extras.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2009 by the University of Oxford
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

#ifndef HAVE_COOT_COORD_EXTRAS_HH
#define HAVE_COOT_COORD_EXTRAS_HH

#include <iostream>
#include <map>

#include "coot-coord-utils.hh"
#include "protein-geometry.hh"
#include "atom-quads.hh"
#include "mini-mol.hh"
#include "bonded-pairs.hh"


namespace coot {

   namespace util { 
      // geom_p gets updated to include the residue restraints if necessary
      // 
      std::pair<int, std::vector<std::string> >
      check_dictionary_for_residues(PCResidue *SelResidues, int nSelResidues,
				    protein_geometry *geom_p, int read_number);

      // We also now pass regular_residue_flag so that the indexing of the
      // contacts is inverted in the case of not regular residue.  I don't
      // know why this is necessary, but I have stared at it for hours, this
      // is a quick (ugly hack) fix that works.  I suspect that there is
      // some atom order dependency in mgtree that I don't understand.
      // Please fix (remove the necessity of depending on
      // regular_residue_flag) if you know how.
      // 
      std::vector<std::vector<int> >
      get_contact_indices_from_restraints(CResidue *residue,
					  protein_geometry *geom_p,
					  bool regular_residue_flag,
					  bool add_reverse_contacts);

      std::vector<std::vector<int> >
      get_contact_indices_for_PRO_residue(PPCAtom residue_atom,
					  int nResidueAtoms, 
					  coot::protein_geometry *geom_p);

      class missing_atom_info {
      public:
	 std::vector<std::string> residues_with_no_dictionary;
	 std::vector<CResidue *>  residues_with_missing_atoms;
	 std::vector<std::pair<CResidue *, std::vector<CAtom *> > > atoms_in_coords_but_not_in_dict;
	 missing_atom_info() {}
	 missing_atom_info(const std::vector<std::string> &residues_with_no_dictionary_in,
			   const std::vector<CResidue *>  &residues_with_missing_atoms_in,
			   const std::vector<std::pair<CResidue *, std::vector<CAtom *> > > &atoms_in_coords_but_not_in_dict_in) {
	    residues_with_no_dictionary = residues_with_no_dictionary_in;
	    residues_with_missing_atoms = residues_with_missing_atoms_in;
	    atoms_in_coords_but_not_in_dict = atoms_in_coords_but_not_in_dict_in;
	 }
      };

      class dict_atom_info_t {
      public:
	 std::string name;
	 short int is_Hydrogen_flag;
	 dict_atom_info_t(const std::string &name_in, short int is_Hydrogen_flag_in) {
	    name = name_in;
	    is_Hydrogen_flag = is_Hydrogen_flag_in;
	 }
      };

      // a trivial helper class
      class dict_residue_atom_info_t {
      public:
	 std::string residue_name;
	 std::vector<dict_atom_info_t> atom_info;
	 dict_residue_atom_info_t(const std::string &residue_name_in,
				  const std::vector<dict_atom_info_t> &atom_info_in) {
	    residue_name = residue_name_in;
	    atom_info = atom_info_in;
	 }
	 // Here is the clever stuff, get the restraints info for a
	 // particular residue and from that set the atom_info.
	 dict_residue_atom_info_t(const std::string &residue_name,
				  protein_geometry *geom_p);
	 bool is_empty_p() const {
	    return (atom_info.size() == 0);
	 }
      };

      // do we need to pass read number to this function too?
      bool is_nucleotide_by_dict_dynamic_add(CResidue *residue_p, coot::protein_geometry *geom_p);
      bool is_nucleotide_by_dict(CResidue *residue_p, const coot::protein_geometry &geom);
   }

   CResidue *GetResidue(const minimol::residue &r); // For use with wiggly
					   // ligands, constructed from
					   // a minimol residue, the
					   // get_contact_indices_from_restraints()
					   // needs a CResidue *.  Caller disposes.
      


   class recursive_forwards_container_t {
   public:
      bool done;
      std::vector<int> forwards;
      recursive_forwards_container_t(bool done_in, const std::vector<int> &forwards_in) {
	 done = done_in;
	 forwards = forwards_in;
      }
      recursive_forwards_container_t() {
	 done = 0;
      } 
   };

   class atom_vertex {

   public:
      enum connection_type_t { START, END, STANDARD, NONE };
      connection_type_t connection_type;
      std::vector<int> forward;
      std::vector<int> backward;
      std::pair<bool,atom_index_quad> torsion_quad;
      atom_vertex() {
	 connection_type = NONE;
	 torsion_quad.first = 0;
      }
   }; 

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
      CResidue *residue;
      PPCAtom atom_selection;   // for the multi-residue interface
      int     n_selected_atoms; // (we can't do residue_p->GetAtomTable())
      bool made_from_minimol_residue_flag; 
      std::vector<std::pair<int, int> > bonds;
      std::vector<atom_vertex> atom_vertex_vec;
      void construct_internal(const dictionary_residue_restraints_t &rest,
			      CResidue *res,
			      const std::string &altconf);
      std::map<std::string, map_index_t, std::less<std::string> > name_to_index;
   private: 
      bool fill_atom_vertex_vec(const dictionary_residue_restraints_t &rest, CResidue *res,
				const std::string &altconf);
      bool fill_atom_vertex_vec_using_contacts(const std::vector<std::vector<int> > &contact_indices,
					       int base_atom_index);
      bool fill_torsions(const dictionary_residue_restraints_t &rest, CResidue *res,
			 const std::string &altconf);
      void fill_name_map(const std::string &altconf);

      // return bool of 0 on not able to fill (not an exception).
      // 
      std::pair<bool, atom_index_quad>
      get_atom_index_quad(const dict_torsion_restraint_t &tr,
			  CResidue *res, const std::string &altconf) const;
      
      std::vector<map_index_t> get_back_atoms(const map_index_t &index2) const;

      // don't add base index to forward atoms (of base_index).
      std::pair<int, std::vector<map_index_t> > get_forward_atoms(const map_index_t &base_index,
								  const map_index_t &index2) const;
      std::vector<map_index_t>
      uniquify_atom_indices(const std::vector<map_index_t> &vin) const;

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

      // factored out:
      bool fill_atom_vertex_vec_using_contacts_by_atom_selection(const std::vector<std::vector<int> > &contact_indices,
								 PPCAtom residue_atoms,
								 int n_residue_atoms,
								 int base_atom_index);
      
      
   public:

      // The angles are in degrees (they get converted to radians in
      // set_dihedral()).
      // 
      class tree_dihedral_info_t {
      public:
	 atom_name_quad quad;
	 double dihedral_angle;
	 tree_dihedral_info_t(const atom_name_quad quad_in, double ang_in) {
	    quad = quad_in;
	    dihedral_angle = ang_in;
	 }
	 tree_dihedral_info_t() {}
	 friend std::ostream& operator<<(std::ostream &o, tree_dihedral_info_t t);
      };
      
      // the constructor throws an exception if there is no tree in
      // the restraints.
      // 
      // FIXME.  In fact, we can handle there being no tree.  We
      // should fall back to using the bonds in the restraint. And
      // throw an exception if there are no bonds.
      //
      atom_tree_t(const dictionary_residue_restraints_t &rest, CResidue *res,
		  const std::string &altconf);

      // the constructor, given a list of bonds and a base atom index.
      // Used perhaps as the fallback when the above raises an
      // exception.
      // 
      atom_tree_t(const std::vector<std::vector<int> > &contact_indices,
		  int base_atom_index, 
		  CResidue *res,
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
		  CMMDBManager *mol,
		  int selection_handle);
      

      ~atom_tree_t() {
	 if (made_from_minimol_residue_flag)
	    // question: does this delete the atoms in the residue too?
	    delete residue; 
      }


      // Rotate round the 2 middle atoms of the torsion by angle (in
      // degress).  This is a relative rotation - not setting the
      // torsion angle.  The atoms of the CResidue residue are
      // manipulated.
      //
      // The reversed flag (being true) allows the rotation of the
      // base, rather than the branch (dog wags rather than tail).
      // atom1 and atom2 can be passed in either way round, this
      // function will sort out the position in the tree.  The
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


      // this can throw an exception
      // 
      // return the set of angles - should be the same that they were
      // set to (for validation).
      //
      // angle in degrees
      std::vector<double> set_dihedral_multi(const std::vector<tree_dihedral_info_t> &di);
      //
      minimol::residue GetResidue() const; // for use with above
					   // function, where the
					   // class constructor is
					   // made with a minimol::residue.


   };
   std::ostream& operator<<(std::ostream &o, atom_tree_t::tree_dihedral_info_t t);

   // Move the atoms in res_moving.
   // 
   // Return the number of rotated torsions.
   //
   class match_torsions {
      enum { REFERENCE_TORSION, MOVING_TORSION};
      CResidue *res_moving;
      CResidue *res_ref;
      dictionary_residue_restraints_t moving_residue_restraints;
      std::pair<bool, double> apply_torsion(const atom_name_quad &quad_moving,
					    const atom_name_quad &reference,
					    const std::string &alt_conf);
      // return in radians
      std::pair<bool, double> get_torsion(int torsion_type, const coot::atom_name_quad &quad) const;
      // return in radians
      std::pair<bool, double> get_torsion(CResidue *res, const coot::atom_name_quad &quad) const;
      
   public: 
      match_torsions(CResidue *res_moving, CResidue *res_ref,
		     const dictionary_residue_restraints_t &moving_residue_restraints_in);
      int match (const std::vector<dict_torsion_restraint_t>  &tr_moving,
		 const std::vector<dict_torsion_restraint_t>  &tr_ref);

   };


   // which uses:
   //
   // (just the torsionable bonds (middle atom pairs)) by looking at
   // the monomer restraints
   // 
   std::vector<std::pair<CAtom *, CAtom *> >
   torsionable_bonds_monomer_internal(CResidue *residue_p,
				      PPCAtom atom_selection, int n_selected_atoms,
				      bool include_pyranose_ring_torsions_flag,
				      protein_geometry *geom_p);

   // only uses the LINKR records in the first model.
   //
   bonded_pair_container_t    
   linkrs_in_atom_selection(CMMDBManager *mol, PPCAtom atom_selection, int n_selected_atoms,
			    protein_geometry *geom_p);
      


}

#endif // HAVE_COOT_COORD_EXTRAS_HH
