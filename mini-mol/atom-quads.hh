/* mini-mol/atom-quads.hh
 * 
 * Copyright  2009 The University of Oxford
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


#ifndef LIGAND_ATOM_QUAD_HH
#define LIGAND_ATOM_QUAD_HH

#include <vector>
#include <stdexcept>

#ifndef HAVE_STRING
#include <string>
#define HAVE_STRING
#endif

#include <mmdb2/mmdb_manager.h>

namespace coot { 

   // 4 atom names for a torsion:
   class atom_name_quad {
      std::string atom_name_[4];
   public:

      // an atom quad can (now) refer to atoms in different residues.
      // So make placeholders for the residue index (1 or 2)
      int atom_residue_index[4]; 
				 
      atom_name_quad() {} // sigh.  needed for torsioned_atoms_info_t,
			  // don't know why.
      atom_name_quad(const std::string &atom_name_0,
                     const std::string &atom_name_1,
                     const std::string &atom_name_2,
                     const std::string &atom_name_3) {

	 atom_name_[0] = atom_name_0;
	 atom_name_[1] = atom_name_1;
	 atom_name_[2] = atom_name_2;
	 atom_name_[3] = atom_name_3;
	 
	 atom_residue_index[0] = 1;
	 atom_residue_index[1] = 1;
	 atom_residue_index[2] = 1;
	 atom_residue_index[3] = 1;
      }
      bool all_non_blank() const {
	 if ((atom_name_[0] == "") || (atom_name_[1] == "") || (atom_name_[2] == "") || (atom_name_[3] == ""))
	    return 0;
	 else
	    return 1;
      } 
      const std::string &atom_name(int i) const {
	 if (i>=0 && i<4) { 
	    return atom_name_[i];
	 } else {
	    throw std::runtime_error("out of bounds index on atom_name_quad::atom_name()");
	 }
      }
      void set_atom_residue_index(int atom_index, int residue_index) {
	 if (atom_index>=0 && atom_index<4) {
	    if (residue_index == 1 || residue_index == 2) {
	       atom_residue_index[atom_index] = residue_index;
	    }
	 } 
      } 
      double torsion(mmdb::Residue *residue) const;
      friend std::ostream& operator<<(std::ostream &o, const atom_name_quad &q);
   };
   std::ostream& operator<<(std::ostream &o, const atom_name_quad &q);

   class atom_name_torsion_quad : public atom_name_quad {
   public:
      std::string id;
      double torsion;
      atom_name_torsion_quad(const std::string &id_in,
			     const std::string &at_name_1_in,
			     const std::string &at_name_2_in,
			     const std::string &at_name_3_in,
			     const std::string &at_name_4_in,
			     double tors_in) : atom_name_quad(at_name_1_in,
							      at_name_2_in,
							      at_name_3_in,
							      at_name_4_in) {
	 torsion = tors_in;
	 id = id_in;
      }
   };

   class atom_index_quad {
   public:
      int index1;
      int index2;
      int index3;
      int index4;
      atom_index_quad() { // unassigned pair
         index1 = -1;
         index2 = -1;
         index3 = -1;
         index4 = -1;
      }
      atom_index_quad(int i1, int i2, int i3, int i4) {
         index1 = i1;
         index2 = i2;
         index3 = i3;
         index4 = i4;
      }

      // Return the torsion in degrees.  Use the indices to index into
      // residue res to find the atoms for the torsion.
      // 
      // Can throw an exception on not all indices found
      // 
      double torsion(mmdb::Residue *res) const;

      // as above, but we reference atoms in a selection rather than a
      // residue.
      // 
      double torsion(mmdb::PPAtom atom_selection, int n_selected_atoms) const;
      
      friend std::ostream& operator<<(std::ostream &o, const atom_index_quad &q);
   };
   std::ostream& operator<<(std::ostream &o, const atom_index_quad &q);

   class atom_quad {
      atom_quad setup_chiral_quad(mmdb::Residue *residue_with_O, mmdb::Residue *residue_with_chiral_centre,
				  const std::string &O_name,
				  const std::vector<std::string> &chiral_atom_names) const;
      
   public:
      mmdb::Atom *atom_1;
      mmdb::Atom *atom_2;
      mmdb::Atom *atom_3;
      mmdb::Atom *atom_4;
      std::string name;
      atom_quad(mmdb::Atom *atom_1_in, mmdb::Atom *atom_2_in, mmdb::Atom *atom_3_in, mmdb::Atom *atom_4_in) {
	 atom_1 = atom_1_in;
	 atom_2 = atom_2_in;
	 atom_3 = atom_3_in;
	 atom_4 = atom_4_in;
	 name = "";
      }
      atom_quad() {
	 atom_1 = NULL;
	 atom_2 = NULL;
	 atom_3 = NULL;
	 atom_4 = NULL;
	 name = "";
      }
      atom_quad(mmdb::Residue *first, mmdb::Residue *second, const std::string &link);
      bool filled_p() const; // ! were there any nulls?
      
      friend std::ostream& operator<<(std::ostream &o, const atom_quad &q);
      // Can throw a std::runtime_error if any of the atoms are null.
      double angle_2() const; // angle 1-2-3 in degrees
      double angle_3() const; // angle 2-3-4 in degrees
      // Can throw a std::runtime_error if any of the atoms are null.
      double torsion() const; // in degrees
      // Can throw a std::runtime_error if any of the atoms are null.
      double chiral_volume() const;
      atom_quad reverse() const {
	 return atom_quad(atom_4, atom_3, atom_2, atom_1);
      }
      // Can throw a std::runtime_error if any of the atoms are null.
      atom_name_quad get_atom_name_quad() const;
   }; 
   std::ostream& operator<<(std::ostream &o, const atom_quad &q);

   class torsion_atom_quad : public atom_quad {
   public:
      int period;
      double angle;
      double angle_esd;
      std::string residue_name; // set in the case of monomer torsions.
      torsion_atom_quad() : atom_quad() {}
      torsion_atom_quad(mmdb::Atom *atom_1_in,
			mmdb::Atom *atom_2_in,
			mmdb::Atom *atom_3_in,
			mmdb::Atom *atom_4_in,
			double angle_in, double angle_esd_in, int period_in) :
	 atom_quad(atom_1_in, atom_2_in, atom_3_in, atom_4_in) {
	 period = period_in;
	 angle = angle_in;
	 angle_esd = angle_esd_in;
      }
   };

}

#endif // LIGAND_ATOM_QUAD_HH

