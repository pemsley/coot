/* coot-utils/coot-coord-extras.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 by The University of York
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

#include "protein-geometry.hh"

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
					  short int regular_residue_flag);

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
      short int is_nucleotide_by_dict_dynamic_add(CResidue *residue_p, coot::protein_geometry *geom_p);
      short int is_nucleotide_by_dict(CResidue *residue_p, const coot::protein_geometry &geom);
   }

}

#endif // HAVE_COOT_COORD_EXTRAS_HH
