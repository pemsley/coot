/* geometry/protein-geometry.cc
 * 
 * Copyright 2004, 2005 The University of York
 * Copyright 2008, 2009, 2011, 2012 The University of Oxford
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

#ifndef PROTEIN_GEOMETRY_HH
#define PROTEIN_GEOMETRY_HH

#define MONOMER_DIR_STR "COOT_SBASE_DIR"

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#include <map>
#include <stdexcept>

#include <mmdb/mmdb_mmcif.h>

#ifndef  __MMDB_SBase__
#include <mmdb/mmdb_sbase.h>
#endif

#include "clipper/core/coords.h"

#include "atom-quads.hh"

namespace coot {

   std::string atom_id_mmdb_expand(const std::string &atomname);
   std::string atom_id_mmdb_expand(const std::string &atomname, const std::string &element);

   class pdbx_chem_comp_descriptor_item {
   public:
      std::string type;
      std::string program;
      std::string program_version;
      std::string descriptor;
      pdbx_chem_comp_descriptor_item(const std::string &type_in,
				     const std::string &program_in,
				     const std::string &program_version_in,
				     const std::string &descriptor_in) {
	 type = type_in;
	 program = program_in;
	 program_version = program_version_in;
	 descriptor = descriptor_in;
      }
   };
   
   class pdbx_chem_comp_descriptor_container_t {
   public:
      std::vector<pdbx_chem_comp_descriptor_item> descriptors;
   };
   
   
   class dict_chem_comp_t {
      void setup_internal(const std::string &comp_id_in,
			  const std::string &three_letter_code_in,
			  const std::string &name_in,
			  const std::string &group_in,
			  int number_atoms_all_in,
			  int number_atoms_nh_in,
			  const std::string &description_level_in) {
	 comp_id = comp_id_in;
	 three_letter_code = three_letter_code_in;
	 name = name_in;
	 group = group_in;
	 number_atoms_all = number_atoms_all_in;
	 number_atoms_nh = number_atoms_nh_in;
	 description_level = description_level_in;
      }
   public:
      std::string comp_id;
      std::string three_letter_code;
      std::string name;
      std::string group; // e.g. "L-peptide"
      int number_atoms_all;
      int number_atoms_nh;
      std::string description_level;
      dict_chem_comp_t(const std::string &comp_id_in,
		       const std::string &three_letter_code_in,
		       const std::string &name_in,
		       const std::string &group_in,
		       int number_atoms_all_in,
		       int number_atoms_nh_in,
		       const std::string &description_level_in) {
	 setup_internal(comp_id_in,
			three_letter_code_in, name_in, group_in,
			number_atoms_all_in, number_atoms_nh_in,
			description_level_in);
      }
      dict_chem_comp_t() {
 	 setup_internal("", "", "", "", 0, 0, "");
      }
      dict_chem_comp_t(const dict_chem_comp_t &din) {
	 setup_internal(din.comp_id, din.three_letter_code, din.name,
			din.group, din.number_atoms_all, din.number_atoms_nh,
			din.description_level);
      } 
   };


   class basic_dict_restraint_t {
      std::string atom_id_1_;
      std::string atom_id_2_;

   public:
      basic_dict_restraint_t() {} // for planes
      basic_dict_restraint_t(const std::string &at1,
			     const std::string &at2);
      std::string atom_id_1() const { return atom_id_1_;}
      std::string atom_id_2() const { return atom_id_2_;}
      std::string atom_id_1_4c() const {  // 4 character return;
	 return atom_id_mmdb_expand(atom_id_1_);
      }
      std::string atom_id_2_4c() const {
	 return atom_id_mmdb_expand(atom_id_2_);
      }
      
   }; 

   class dict_bond_restraint_t : public basic_dict_restraint_t {
      std::string type_;
      double dist_;
      double dist_esd_;
      bool have_target_values; 
   
   public:
      // dict_bond_restraint_t() {};
      dict_bond_restraint_t(std::string atom_id_1_in,
			    std::string atom_id_2_in,
			    std::string type,
			    double dist_in,
			    double dist_esd_in) :
      basic_dict_restraint_t(atom_id_1_in, atom_id_2_in) {

	 dist_ = dist_in;
	 dist_esd_ = dist_esd_in;
	 type_ = type;
	 have_target_values = 1;
      }


      dict_bond_restraint_t(std::string atom_id_1_in,
			    std::string atom_id_2_in,
			    std::string type) :
	 basic_dict_restraint_t(atom_id_1_in, atom_id_2_in) {
	 have_target_values = 0;
	 type_ = type;
      }
      
      std::string type() const { return type_; }
      // can throw a std::runtime_error exception (if target values not set)
      double value_dist() const {
	 if (have_target_values)
	    return dist_;
	 else
	    throw std::runtime_error("unset target distance geometry");
      }
      // can throw a std::runtime_error exception
      double value_esd () const {
	 if (have_target_values)
	    return dist_esd_;
	 else
	    throw std::runtime_error("unset target distance geometry");
      }
   };

   class dict_angle_restraint_t : public basic_dict_restraint_t {
      std::string atom_id_3_;
      double angle_;
      double angle_esd_;
   public:
      // dict_angle_restraint_t() {};
      dict_angle_restraint_t(std::string atom_id_1,
			     std::string atom_id_2,
			     std::string atom_id_3,
			     double angle,
			     double angle_esd) :
      basic_dict_restraint_t(atom_id_1, atom_id_2) {
	 atom_id_3_ = atom_id_3;
	 angle_ = angle;
	 angle_esd_ = angle_esd; 
      };
      
      std::string atom_id_3_4c() const { return atom_id_mmdb_expand(atom_id_3_);}
      double angle() const { return angle_; }
      double esd ()  const { return angle_esd_;}
      friend std::ostream& operator<<(std::ostream &s, const dict_angle_restraint_t &rest);
   };
   std::ostream& operator<<(std::ostream &s, const dict_angle_restraint_t &rest);

   // Note hydrogen torsions can only be detected at the container
   // (protein_geometry) level, because we don't have acces to the
   // elements here (only the atom names).
   // 
   class dict_torsion_restraint_t : public basic_dict_restraint_t {
      std::string id_;
      std::string atom_id_3_;
      std::string atom_id_4_;
      double angle_;
      double angle_esd_;
      int period;
   public:
      
      // dict_torsion_restraint_t() {}; 
      dict_torsion_restraint_t(std::string id_in,
			       std::string atom_id_1,
			       std::string atom_id_2,
			       std::string atom_id_3,
			       std::string atom_id_4,
			       double angle,
			       double angle_esd,
			       int period_in) :
      basic_dict_restraint_t(atom_id_1, atom_id_2)
      {
	 id_ = id_in;
	 atom_id_3_ = atom_id_3;
	 atom_id_4_ = atom_id_4;
	 angle_ = angle;
	 angle_esd_ = angle_esd;
	 period = period_in;
      };
      std::string atom_id_3_4c() const { return atom_id_mmdb_expand(atom_id_3_);}
      std::string atom_id_4_4c() const { return atom_id_mmdb_expand(atom_id_4_);}
      std::string atom_id_3() const { return atom_id_3_;}
      std::string atom_id_4() const { return atom_id_4_;}
      std::string id() const { return id_;}
      bool is_const() const; // is the id const?  (Don't consider angle or sd).
      int periodicity() const { return period; }
      double angle() const { return angle_; }
      double esd ()  const { return angle_esd_;}
      friend std::ostream& operator<<(std::ostream &s, const dict_torsion_restraint_t &rest);
      bool is_pyranose_ring_torsion() const;
      bool is_ring_torsion(const std::vector<std::vector<std::string> > &ring_atoms_sets) const;
      // hack for mac, ostream problems
      std::string format() const;
   };
   std::ostream& operator<<(std::ostream &s, const dict_torsion_restraint_t &rest); 

   // ------------------------------------------------------------------------
   // class dict_plane_restraint_t 
   // ------------------------------------------------------------------------
   // 
   class dict_plane_restraint_t : public basic_dict_restraint_t {
      std::vector<std::string> atom_ids;
      double dist_esd_;  // despite separate entries for each atom in
			// the dictionary, I decide that it is
			// sensible to say that all atoms in a plane
			// have the same esd deviation from the plane.
   public:
      dict_plane_restraint_t() {};
      dict_plane_restraint_t(const std::string &plane_id_in,
			     const std::string &atom_id_in,
			     double dist_esd_in) {
	 plane_id = plane_id_in;
	 dist_esd_ = dist_esd_in;
	 atom_ids.push_back(atom_id_in);
      };
      dict_plane_restraint_t(const std::string &plane_id_in,
			     const std::vector<std::string> &plane_atom_ids,
			     double dist_esd_in) {
	 plane_id = plane_id_in;
	 dist_esd_ = dist_esd_in;
	 atom_ids = plane_atom_ids;
      };
      std::string plane_id; // or int plane_id number 1, 2 3.
      double dist_esd() const { return dist_esd_; }
      std::string atom_id(int i) const { return atom_id_mmdb_expand(atom_ids[i]); }
      int n_atoms() const { return atom_ids.size(); }
      const std::string &operator[](int i) const { return atom_ids[i];}
      void push_back_atom(const std::string &at) { atom_ids.push_back(at); }
      friend std::ostream&  operator<<(std::ostream &s, dict_plane_restraint_t rest);
   };

   std::ostream&  operator<<(std::ostream &s, dict_plane_restraint_t rest);
   
   // ------------------------------------------------------------------------
   // class dict_chiral_restraint_t
   // ------------------------------------------------------------------------
   // 
   class dict_chiral_restraint_t : public basic_dict_restraint_t {
      bool is_both_flag; 
      std::string chiral_id; // or int chiral_id 1, 2, 3.
      std::string local_atom_id_centre; // CA usually.
      std::string local_atom_id_1;
      std::string local_atom_id_2;
      std::string local_atom_id_3;
      double target_volume_;
      double volume_sigma_;
      // angles in radians.
      double assign_chiral_volume_target_internal(double a, double b, double c,
						  double alpha, double beta, double gamma);
   public:
      int volume_sign;  // +/- 1, checked by is_bad_chiral_atom_p and
			// set by nomenclature checking
      enum { CHIRAL_RESTRAINT_BOTH = -2,
	     CHIRAL_VOLUME_RESTRAINT_VOLUME_SIGN_UNASSIGNED = -3 };
      
      dict_chiral_restraint_t() {};
      dict_chiral_restraint_t(const std::string &chiral_id_in,
			      const std::string &atom_id_centre_in,
			      const std::string &atom_id_1_in,
			      const std::string &atom_id_2_in,
			      const std::string &atom_id_3_in,
			      int volume_sign_in) {

	 chiral_id = chiral_id_in;
	 local_atom_id_centre = atom_id_centre_in;
	 local_atom_id_1 = atom_id_1_in;
	 local_atom_id_2 = atom_id_2_in;
	 local_atom_id_3 = atom_id_3_in;
	 volume_sign = volume_sign_in;
	 target_volume_ = -999.9;  // unassigned
	 volume_sigma_  = -999.9;
	 is_both_flag = 0;
	 if (volume_sign_in == CHIRAL_RESTRAINT_BOTH) { 
	    is_both_flag = 1;
	    volume_sigma_ = 1.0 ; // mark as assigned
	 } 
      }
      std::string Chiral_Id() const { return chiral_id; }
      std::string atom_id_1_4c() const { return atom_id_mmdb_expand(local_atom_id_1);}
      std::string atom_id_2_4c() const { return atom_id_mmdb_expand(local_atom_id_2);}
      std::string atom_id_3_4c() const { return atom_id_mmdb_expand(local_atom_id_3);}
      std::string atom_id_c_4c() const { return atom_id_mmdb_expand(local_atom_id_centre);}
      double assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds,
					 const std::vector <dict_angle_restraint_t> &angles);
      double target_volume() const { return target_volume_;}
      double volume_sigma()  const { return volume_sigma_;}
      bool is_a_both_restraint() const { return is_both_flag; } 
      bool has_unassigned_chiral_volume() const {
	 return (volume_sigma_ < 0.0) ? 1 : 0;
      }

   };


   // ------------------------------------------------------------------------
   // class dict_atom
   // ------------------------------------------------------------------------
   // 
   // one of these for each atom in a dictionary_residue_restraints_t
   // (i.e. each atom in a residue/comp_id).
   // 
   class dict_atom {
   public:
      enum { IDEAL_MODEL_POS, REAL_MODEL_POS}; 
      std::string atom_id;
      std::string atom_id_4c;
      std::string type_symbol;
      std::string type_energy;
      std::pair<bool, float> partial_charge;
      std::pair<bool, clipper::Coord_orth> pdbx_model_Cartn_ideal;
      std::pair<bool, clipper::Coord_orth> model_Cartn;
      dict_atom(const std::string &atom_id_in,
		const std::string &atom_id_4c_in,
		const std::string &type_symbol_in,
		const std::string &type_energy_in,
		std::pair<bool, float> partial_charge_in) {
	 atom_id = atom_id_in;
	 atom_id_4c = atom_id_4c_in;
	 type_symbol = type_symbol_in;
	 type_energy = type_energy_in;
	 partial_charge = partial_charge_in;
      }
      dict_atom() {}; // for resize(0);
      void add_pos(int pos_type, const std::pair<bool, clipper::Coord_orth> &model_pos_ideal);
   };

   // ------------------------------------------------------------------------
   // class dict_chem_comp_tree_t
   // ------------------------------------------------------------------------
   // 
   class dict_chem_comp_tree_t  : public basic_dict_restraint_t {
   public:
      std::string atom_id;
      std::string atom_back;
      std::string atom_forward;
      std::string connect_type;
      dict_chem_comp_tree_t(const std::string &atom_id_in, const std::string &atom_back_in,
			    const std::string &atom_forward_in, const std::string &connect_type_in) {
	 atom_id = atom_id_mmdb_expand(atom_id_in);
	 atom_back = atom_id_mmdb_expand(atom_back_in);
	 atom_forward = atom_id_mmdb_expand(atom_forward_in);
	 connect_type = connect_type_in;
      }
      dict_chem_comp_tree_t() {}
   };

   // ------------------------------------------------------------------------
   // class dictionary_residue_restraints_t
   // ------------------------------------------------------------------------
   // 
   class dictionary_residue_restraints_t {
      bool has_partial_charges_flag;
      bool filled_with_sbase_data_flag; // if set, this means that
					// there is only bond orders
					// (at the moment) and atom
					// names.
   public:
      dictionary_residue_restraints_t(std::string comp_id_in,
				      int read_number_in) {
	 has_partial_charges_flag = 0;
	 // comp_id = comp_id_in;
	 read_number = read_number_in;
	 filled_with_sbase_data_flag = 0;
      }
      dictionary_residue_restraints_t() {
	 // comp_id = ""; /* things are unset */
	 filled_with_sbase_data_flag = 0;
	 has_partial_charges_flag = 0;
	 read_number = -1;
      }
      dictionary_residue_restraints_t(bool constructor_for_sbase_restraints) {
	 // comp_id = ""; /* things are unset */
	 filled_with_sbase_data_flag = 1;
	 has_partial_charges_flag = 0;
	 read_number = -1;
      }
      std::string cif_file_name;
      void clear_dictionary_residue();
      dict_chem_comp_t residue_info;
      std::vector <dict_atom> atom_info;
      // std::string comp_id; // i.e. residue type name
      std::string comp_id() const { return residue_info.comp_id; }
      std::vector <dict_chem_comp_tree_t> tree;
      int read_number; 
      std::vector    <dict_bond_restraint_t>    bond_restraint;
      std::vector   <dict_angle_restraint_t>   angle_restraint;
      std::vector <dict_torsion_restraint_t> torsion_restraint;
      std::vector  <dict_chiral_restraint_t>  chiral_restraint;
      std::vector   <dict_plane_restraint_t>   plane_restraint;
      pdbx_chem_comp_descriptor_container_t descriptors;
      // Return 1 for hydrogen or deuterium, 0 for not found or not a hydrogen.
      bool is_hydrogen(const std::string &atom_name) const;
      int assign_chiral_volume_targets(); // return the number of targets made.
      bool has_unassigned_chiral_volumes() const;
      bool has_partial_charges() const;
      void set_has_partial_charges(bool state) {
	 has_partial_charges_flag = state;
      }
      std::vector<dict_torsion_restraint_t> get_non_const_torsions(bool include_hydrogen_torsions_flag) const;

      // compares atoms of torsion_restraint vs the ring atoms.
      // bool is_ring_torsion(const dict_torsion_restraint_t &torsion_restraint) const;
      bool is_ring_torsion(const atom_name_quad &quad) const;
      
      void write_cif(const std::string &filename) const;
      // look up the atom id in the atom_info (dict_atom vector)
      std::string atom_name_for_tree_4c(const std::string &atom_id) const;

      // look up the atom id in the atom_info (dict_atom vector).
      // Return "" on no atom found with name atom_name;
      // 
      std::string element(const std::string &atom_name) const;

      // likewise look up the energy type.  Return "" on no atom found
      // with that atom_name.
      // 
      std::string type_energy(const std::string &atom_name) const;

      std::vector<std::vector<std::string> > get_ligand_ring_list() const;

      std::vector<std::vector<std::string> > get_ligand_aromatic_ring_list() const;

      std::vector<std::string> get_attached_H_names(const std::string &atom_name) const;

      bool is_from_sbase_data() const { return filled_with_sbase_data_flag; }

      std::vector<std::string> neighbours(const std::string &atom_name, bool allow_hydrogen_neighbours_flag) const;

      // return "" on not found
      std::string get_bond_type(const std::string &name_1, const std::string &name_2) const;

      // replace the restraints that we have with new_restraints,
      // keeping restraints that in the current set bu not in
      // new_restraints
      void conservatively_replace_with(const dictionary_residue_restraints_t &new_restraints);
      void conservatively_replace_with_bonds (const dictionary_residue_restraints_t &new_restraints);
      void conservatively_replace_with_angles(const dictionary_residue_restraints_t &new_restraints);

      //
      void remove_redundant_plane_resetraints();
      bool is_redundant_plane_resetraints(std::vector<dict_plane_restraint_t>::iterator it);
      
   };

   // ------------------------------------------------------------------------
   // class dict_link_bond_restraint_t
   // ------------------------------------------------------------------------
   // 
   class dict_link_bond_restraint_t : public basic_dict_restraint_t {
      double value_dist;
      double value_dist_esd;
   public:
      dict_link_bond_restraint_t (int atom_1_comp_id_in,
				  int atom_2_comp_id_in,
				  const std::string &atom_id_1,
				  const std::string &atom_id_2,
				  realtype value_dist_in,
				  realtype value_dist_esd_in) :
	 basic_dict_restraint_t(atom_id_1, atom_id_2) {
	 atom_1_comp_id = atom_1_comp_id_in;
	 atom_2_comp_id = atom_2_comp_id_in;
	 value_dist = value_dist_in;
	 value_dist_esd = value_dist_esd_in;
      }
      int atom_1_comp_id, atom_2_comp_id;
      double dist() const { return value_dist; }
      double esd()  const { return value_dist_esd; }
   }; 

   class dict_link_angle_restraint_t : public basic_dict_restraint_t {
      double angle_;
      double angle_esd_;
      std::string atom_id_3_;
   public:
      dict_link_angle_restraint_t (int atom_1_comp_id_in,
				   int atom_2_comp_id_in,
				   int atom_3_comp_id_in,
				   const std::string &atom_id_1,
				   const std::string &atom_id_2,
				   const std::string &atom_id_3_in,
				   realtype value_angle_in,
				   realtype value_angle_esd_in) :
	 basic_dict_restraint_t(atom_id_1, atom_id_2) {
	 atom_id_3_ = atom_id_3_in;
	 atom_1_comp_id = atom_1_comp_id_in;
	 atom_2_comp_id = atom_2_comp_id_in;
	 atom_3_comp_id = atom_3_comp_id_in;
	 angle_ = value_angle_in;
	 angle_esd_ = value_angle_esd_in;
      }
      int atom_1_comp_id, atom_2_comp_id, atom_3_comp_id;
      double angle() const { return angle_; }
      double angle_esd() const { return angle_esd_;}
      std::string atom_id_3_4c() const { return atom_id_mmdb_expand(atom_id_3_);}
   }; 

   class dict_link_torsion_restraint_t : public basic_dict_restraint_t {
      double angle_;
      double angle_esd_;
      std::string atom_id_3_;
      std::string atom_id_4_;
      std::string id_;
      int period_;

   public:
         
      dict_link_torsion_restraint_t (int atom_1_comp_id_in,
				     int atom_2_comp_id_in,
				     int atom_3_comp_id_in,
				     int atom_4_comp_id_in,
				     const std::string &atom_id_1_in,
				     const std::string &atom_id_2_in,
				     const std::string &atom_id_3_in,
				     const std::string &atom_id_4_in,
				     double value_angle_in,
				     double value_angle_esd,
				     int period,
				     const std::string &id) :
	 basic_dict_restraint_t(atom_id_1_in, atom_id_2_in)  {

	 atom_id_3_ = atom_id_3_in;
	 atom_id_4_ = atom_id_4_in;
	 atom_1_comp_id = atom_1_comp_id_in;
	 atom_2_comp_id = atom_2_comp_id_in;
	 atom_3_comp_id = atom_3_comp_id_in;
	 atom_4_comp_id = atom_4_comp_id_in;
	 angle_ = value_angle_in;
	 angle_esd_ = value_angle_esd;
	 period_ = period;
	 id_ = id;
      }
      int atom_1_comp_id, atom_2_comp_id, atom_3_comp_id, atom_4_comp_id;
      double angle() const { return angle_; }
      double angle_esd() const { return angle_esd_;}
      std::string atom_id_3_4c() const { return atom_id_mmdb_expand(atom_id_3_);}
      std::string atom_id_4_4c() const { return atom_id_mmdb_expand(atom_id_4_);}
      int period() const { return period_; }
      std::string id() const { return id_;}
      bool is_pyranose_ring_torsion() const; 
   }; 

   class dict_link_chiral_restraint_t : public basic_dict_restraint_t {
   public:
      int atom_1_comp_id, atom_2_comp_id, atom_3_comp_id, atom_c_comp_id;
   private:
      std::string atom_id_c_;
      std::string atom_id_3_;
      std::string id_;
      int volume_sign;
      double target_volume_;
      double target_sigma_;
      std::string chiral_id;
   public:
      dict_link_chiral_restraint_t(const std::string &chiral_id_in,
				   int atom_c_comp_id_in,
				   int atom_1_comp_id_in,
				   int atom_2_comp_id_in,
				   int atom_3_comp_id_in,
				   const std::string &atom_id_c_in,
				   const std::string &atom_id_1_in,
				   const std::string &atom_id_2_in,
				   const std::string &atom_id_3_in,
				   int volume_sign_in
				   ) : basic_dict_restraint_t(atom_id_1_in, atom_id_2_in) {
	 atom_c_comp_id = atom_c_comp_id_in;
	 atom_1_comp_id = atom_1_comp_id_in;
	 atom_2_comp_id = atom_2_comp_id_in;
	 atom_3_comp_id = atom_3_comp_id_in;
	 atom_id_c_ = atom_id_c_in;
	 atom_id_3_ = atom_id_3_in;
	 volume_sign = volume_sign_in;
	 target_volume_ = -999.9;  // unassigned
	 target_sigma_  = -999.9;
	 chiral_id = chiral_id_in;
      }
      std::string Chiral_Id() const { return chiral_id; }
      double assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds_1,
					 const std::vector <dict_angle_restraint_t> &angles_1,
					 const std::vector <dict_bond_restraint_t> &bonds_2,
					 const std::vector <dict_angle_restraint_t> &angles_2,
					 const std::vector <dict_link_bond_restraint_t> &link_bonds,
					 const std::vector <dict_link_angle_restraint_t> &link_angles);
      double target_volume() const { return target_volume_; }
      double target_sigma() const  { return target_sigma_; }
      bool has_unassigned_chiral_volume() const {
	 return (target_sigma_ < 0.0);
      }
   };


   class dict_link_plane_restraint_t : public basic_dict_restraint_t {
      double dist_esd_;
   public:
      dict_link_plane_restraint_t(const std::string &atom_id,
				  const std::string &plane_id_in,
				  int atom_comp_id,
				  double dist_esd) {
	 plane_id = plane_id_in;
	 dist_esd_ = dist_esd;
	 atom_ids.push_back(atom_id);
	 atom_comp_ids.push_back(atom_comp_id);
      }
      std::string plane_id; 
      std::vector<std::string> atom_ids;
      std::vector<int> atom_comp_ids;
      int n_atoms() const { return atom_ids.size(); }
      double dist_esd() const { return dist_esd_; }
      std::string atom_id(int i) const { return atom_id_mmdb_expand(atom_ids[i]); }
   };

   // for link_ids such as TRANS, PTRANS (proline), CIS etc.
   // 
   class dictionary_residue_link_restraints_t {
   public:
      dictionary_residue_link_restraints_t(const std::string link_id_in) {
	 link_id = link_id_in;}
      dictionary_residue_link_restraints_t() { link_id = ""; /* things unset */ }
      std::string link_id;
      std::vector <dict_link_bond_restraint_t>    link_bond_restraint;
      std::vector <dict_link_angle_restraint_t>   link_angle_restraint;
      std::vector <dict_link_torsion_restraint_t> link_torsion_restraint;
      std::vector <dict_link_plane_restraint_t>   link_plane_restraint;
      std::vector <dict_link_chiral_restraint_t>  link_chiral_restraint;
      int assign_link_chiral_volume_targets(); // return the number of link targets made.
      bool has_unassigned_chiral_volumes() const;
   };

   class simple_cif_reader {
      std::vector<std::string> names;
      std::vector<std::string> three_letter_codes;
   public:
      simple_cif_reader(const std::string &cif_dictionary_file_name);
      bool has_restraints_for(const std::string &res_type);
   };

   class chem_link {
      std::string id;
      std::string chem_link_comp_id_1;
      std::string chem_link_mod_id_1;
      std::string chem_link_group_comp_1;
      std::string chem_link_comp_id_2;
      std::string chem_link_mod_id_2;
      std::string chem_link_group_comp_2;
      std::string chem_link_name;
   public: 
      chem_link(const std::string &id_in,
		const std::string &chem_link_comp_id_1_in,
		const std::string &chem_link_mod_id_1_in,
		const std::string &chem_link_group_comp_1_in,
		const std::string &chem_link_comp_id_2_in,
		const std::string &chem_link_mod_id_2_in,
		const std::string &chem_link_group_comp_2_in,
		const std::string &chem_link_name_in) {

       id = id_in;
       chem_link_comp_id_1 = chem_link_comp_id_1_in;
       chem_link_mod_id_1 = chem_link_mod_id_1_in;
       chem_link_group_comp_1 = chem_link_group_comp_1_in;
       chem_link_comp_id_2 = chem_link_comp_id_2_in;
       chem_link_mod_id_2 = chem_link_mod_id_2_in; 
       chem_link_group_comp_2 = chem_link_group_comp_2_in;
       chem_link_name = chem_link_name_in;
      }
      friend std::ostream& operator<<(std::ostream &s, chem_link lnk);
      // pair: matches need-order-switch-flag
      std::pair<bool, bool> matches_comp_ids_and_groups(const std::string &comp_id_1,
							const std::string &group_1,
							const std::string &comp_id_2,
							const std::string &group_2) const;
      std::string Id() const { return id; }
      bool is_peptide_link_p() const {
	 if (id == "TRANS" || id == "PTRANS" || id == "NMTRANS" ||
	     id == "CIS"   || id == "PCIS"   || id == "NMCIS")
	    return 1;
	 else
	    return 0;
      }
      std::pair<std::string, std::string> chem_mod_names() const {
	 return std::pair<std::string, std::string> (chem_link_mod_id_1, chem_link_mod_id_2);
      } 
   };
   std::ostream& operator<<(std::ostream &s, chem_link lnk);

   // ------------------------------------------------------------------------
   //                  chem_mods
   // ------------------------------------------------------------------------
   
   // a container for the data_mod_list chem_mods (used to be simply
   // chem_mod, but that name is now used as a class that contains the
   // actual chem mods (with lists of atoms, bonds, angles and so on).
   // 
   class list_chem_mod {
   public:
      std::string name;
      std::string id;
      std::string group_id;
      std::string comp_id;
      list_chem_mod(const std::string &id_in,
		    const std::string &name_in,
		    const std::string &comp_id_in,
		    const std::string &group_id_in) {
	 id = id_in;
	 name = name_in;
	 comp_id = comp_id_in;
	 group_id = group_id_in;
      } 
      friend std::ostream& operator<<(std::ostream &s, list_chem_mod mod);
   };
   std::ostream& operator<<(std::ostream &s, list_chem_mod mod);

   enum chem_mod_function_t { CHEM_MOD_FUNCTION_UNSET,
			      CHEM_MOD_FUNCTION_ADD,
			      CHEM_MOD_FUNCTION_CHANGE,
			      CHEM_MOD_FUNCTION_DELETE };
   
   class chem_mod_atom {
   public:
      chem_mod_atom(const std::string &function_in,
		    const std::string &atom_id_in,
		    const std::string &new_atom_id_in,
		    const std::string &new_type_symbol_in,
		    const std::string &new_type_energy_in,
		    realtype new_partial_charge_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id = atom_id_in;
	 new_atom_id = new_atom_id_in;
	 new_type_symbol = new_type_symbol_in;
	 new_type_energy = new_type_energy_in;
	 new_partial_charge = new_partial_charge_in;
      }
      chem_mod_function_t function;
      std::string atom_id;
      std::string new_atom_id;
      std::string new_type_symbol;
      std::string new_type_energy;
      realtype new_partial_charge;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_atom &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_atom &a);

   class chem_mod_tree {
      chem_mod_function_t function;
      std::string atom_id;
      std::string atom_back;
      std::string back_type;
      std::string atom_forward;
      std::string connect_type;
   public:
      chem_mod_tree (const std::string &function_in,
		     const std::string &atom_id_in,
		     const std::string &atom_back_in,
		     const std::string &back_type_in,
		     const std::string &atom_forward_in,
		     const std::string &connect_type_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id = atom_id_in;
	 atom_back = atom_back_in;
	 back_type = back_type_in;
	 atom_forward = atom_forward_in;
	 connect_type = connect_type_in;
      }
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_tree &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_tree &a);

   class chem_mod_bond {
   public:
      chem_mod_bond(const std::string &function_in,
		    const std::string &atom_id_1_in,
		    const std::string &atom_id_2_in,
		    const std::string &new_type_in,
		    realtype new_value_dist_in,
		    realtype new_value_dist_esd_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id_1 = atom_id_1_in;
	 atom_id_2 = atom_id_2_in;
	 new_type = new_type_in;
	 new_value_dist = new_value_dist_in;
	 new_value_dist_esd = new_value_dist_esd_in;
      }
      chem_mod_function_t function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string new_type;
      realtype new_value_dist;
      realtype new_value_dist_esd;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_bond &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_bond &a);

   class chem_mod_angle {
   public:
      chem_mod_angle(const std::string &function_in,
		     const std::string &atom_id_1_in,
		     const std::string &atom_id_2_in,
		     const std::string &atom_id_3_in,
		     realtype new_value_angle_in,
		     realtype new_value_angle_esd_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id_1 = atom_id_1_in;
	 atom_id_2 = atom_id_2_in;
	 atom_id_3 = atom_id_3_in;
	 new_value_angle = new_value_angle_in;
	 new_value_angle_esd = new_value_angle_esd_in;
      }
      chem_mod_function_t function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      std::string new_type;
      realtype new_value_angle;
      realtype new_value_angle_esd;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_angle &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_angle &a);

   class chem_mod_tor {
   public:
      chem_mod_tor(const std::string &function_in,
		   const std::string &atom_id_1_in,
		   const std::string &atom_id_2_in,
		   const std::string &atom_id_3_in,
		   const std::string &atom_id_4_in,
		   realtype new_value_angle_in,
		   realtype new_value_angle_esd_in,
		   int new_period_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id_1 = atom_id_1_in;
	 atom_id_2 = atom_id_2_in;
	 atom_id_3 = atom_id_3_in;
	 atom_id_4 = atom_id_4_in;
	 new_value_angle = new_value_angle_in;
	 new_value_angle_esd = new_value_angle_esd_in;
	 new_period = new_period_in;
      }
      chem_mod_function_t function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      std::string atom_id_4;
      std::string new_type;
      realtype new_value_angle;
      realtype new_value_angle_esd;
      int new_period;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_tor &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_tor &a);

   class chem_mod_plane {
   public:
      chem_mod_plane(const std::string &plane_id_in,
		     const std::string &function_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 plane_id = plane_id_in;
      }
      chem_mod_function_t function;
      std::string plane_id;
      std::vector<std::pair<std::string, realtype> > atom_id_esd;
      void add_atom(const std::string &atom_id, realtype esd) {
	 std::pair<std::string, realtype> p(atom_id, esd);
	 atom_id_esd.push_back(p);
      }
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_plane &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_plane &a);

   class chem_mod_chir {
   public:
      chem_mod_chir(const std::string &function_in,
		    const std::string &atom_id_centre_in,
		    const std::string &atom_id_1_in,
		    const std::string &atom_id_2_in,
		    const std::string &atom_id_3_in,
		    int new_volume_sign_in) {
	 function = CHEM_MOD_FUNCTION_UNSET;
	 if (function_in == "add")
	    function = CHEM_MOD_FUNCTION_ADD;
	 if (function_in == "delete")
	    function = CHEM_MOD_FUNCTION_DELETE;
	 if (function_in == "change")
	    function = CHEM_MOD_FUNCTION_CHANGE;
	 atom_id_centre = atom_id_centre_in;
	 atom_id_1 = atom_id_1_in;
	 atom_id_2 = atom_id_2_in;
	 atom_id_3 = atom_id_3_in;
	 new_volume_sign = new_volume_sign_in;
      }
      chem_mod_function_t function;
      std::string atom_id_centre;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      int new_volume_sign;
      friend std::ostream& operator<<(std::ostream &s, const chem_mod_chir &a);
   };
   std::ostream& operator<<(std::ostream &s, const chem_mod_chir &a);

   


   // ------------------------------------------------------------------------
   //                  energy lib
   // ------------------------------------------------------------------------
   
   class energy_lib_atom {
   public:
      enum { HB_UNASSIGNED=-1, HB_NEITHER, HB_DONOR, HB_ACCEPTOR, HB_BOTH, HB_HYDROGEN };
      std::string type;
      realtype weight;
      int hb_type;
      // radii are negative if not assigned.
      realtype vdw_radius;
      realtype vdwh_radius;
      realtype ion_radius;
      std::string element;
      int valency; // negative if unset
      int sp_hybridisation; // negative if unset
      energy_lib_atom(const std::string &type_in,
		      int hb_type_in,
		      float weight_in,
		      float vdw_radius_in,
		      float vdwh_radius_in,
		      float ion_radius_in,
		      const std::string &element_in,
		      int valency_in,
		      int sp_hybridisation_in) {
	 type = type_in;
	 hb_type = hb_type_in;
	 weight = weight_in;
	 vdw_radius  = vdw_radius_in;
	 vdwh_radius = vdwh_radius_in;
	 ion_radius = ion_radius_in;
	 element = element_in;
	 valency = valency_in;
	 sp_hybridisation = sp_hybridisation_in;
      }
      // for the map
      energy_lib_atom() {
	 type = "";
	 hb_type = HB_UNASSIGNED;
	 weight = -1;
	 vdw_radius = -1;
	 ion_radius = -1;
	 element = "";
	 valency = -1;
	 sp_hybridisation = -1;
      } 
      friend std::ostream& operator<<(std::ostream &s, const energy_lib_atom &at);
   };
   std::ostream& operator<<(std::ostream &s, const energy_lib_atom &at);

   

   class energy_lib_bond {
   public:
      std::string atom_type_1;
      std::string atom_type_2;
      std::string type; // single, double, aromatic etc
      float spring_constant; // for energetics
      float length;
      float esd;
      energy_lib_bond() {
	 type = "unset";
	 length = 0;
	 esd = 0;
      }
      energy_lib_bond(const std::string &atom_type_1_in,
		      const std::string &atom_type_2_in,
		      const std::string &type_in,
		      float spring_constant_in,
		      float length_in,
		      float esd_in) {
	 atom_type_1 = atom_type_1_in;
	 atom_type_2 = atom_type_2_in;
	 type = type_in;
	 spring_constant = spring_constant_in;
	 length = length_in;
	 esd = esd_in;
      }
      // Order-dependent.  Call twice - or more.
      bool matches(const std::string &type_1, const std::string &type_2,
		   const std::string &bond_type_in,
		   bool permissive_type) const {
	 
	 bool r = false;
	 if (type == bond_type_in) { 
	    if (atom_type_1 == type_1) {
	       if (atom_type_2 == "") { 
		  if (permissive_type) 
		     r = true;
	       } else {
		  if (atom_type_2 == type_2)
		     r = true;
	       }
	    }
	 }
	 return r;
      } 
      friend std::ostream& operator<<(std::ostream &s, const energy_lib_bond &bond);
   };
   std::ostream& operator<<(std::ostream &s, const energy_lib_bond &bond);

   class energy_lib_angle {
   public:
      std::string atom_type_1;
      std::string atom_type_2;
      std::string atom_type_3;
      float spring_constant; // for energetics
      float angle;
      float angle_esd;
      energy_lib_angle() {
	 angle = 120;
	 angle_esd = 6;
	 spring_constant = 450;
      }
      energy_lib_angle(const std::string &atom_type_1_in,
		       const std::string &atom_type_2_in,
		       const std::string &atom_type_3_in,
		       float spring_constant_in,
		       float value_in,
		       float value_esd_in) {
	 
	 atom_type_1 = atom_type_1_in;
	 atom_type_2 = atom_type_2_in;
	 atom_type_3 = atom_type_3_in;
	 spring_constant = spring_constant_in;
	 angle = value_in;
	 angle_esd = value_esd_in;
      }
      bool matches(const std::string &type_1,
		   const std::string &type_2,
		   const std::string &type_3,
		   bool permissive_1, bool permissive_3) const {

	 bool r = false;
	 // must match the middle atom at least.
	 if (atom_type_2 == type_2) {
	 
	    // first atom matches
	    if (atom_type_1 == type_1) {
	       if (atom_type_3 == type_3)
		  r = true;
	       if (atom_type_3 == "")
		  if (permissive_3)
		     r = true;
	    }

	    // 3rd atom  match
	    if (atom_type_3 == type_3) {
	       if (atom_type_1 == "")
		  if (permissive_1)
		     r = true;
	    }

	    // permissive 1 and 3
	    if (permissive_1 && permissive_3) {
// 	       std::cout << "looking at \"" << atom_type_1 << "\" and \"" << atom_type_3
// 			 << "\"" << std::endl;
	       if (atom_type_1 == "")
		  if (atom_type_3 == "")
		     r = true;
	    }
	 }
	 return r;
      }
   };

   class energy_lib_torsion {
   public:
      std::string atom_type_1;
      std::string atom_type_2;
      std::string atom_type_3;
      std::string atom_type_4;
      std::string label;
      float spring_constant; // for energetics
      float angle;
      int period;
      energy_lib_torsion() {
	 spring_constant = 0.0;
	 angle = 0.0;
	 period = 0;
      }
      energy_lib_torsion(const std::string &atom_type_1_in,
			 const std::string &atom_type_2_in,
			 const std::string &atom_type_3_in,
			 const std::string &atom_type_4_in,
			 float spring_constant_in,
			 float value_in,
			 int period_in) {
	 
	 atom_type_1 = atom_type_1_in;
	 atom_type_2 = atom_type_2_in;
	 atom_type_3 = atom_type_3_in;
	 atom_type_4 = atom_type_4_in;
	 spring_constant = spring_constant_in;
	 angle = value_in;
	 period = period_in;
      }
      // order dependent.  Call twice.
      bool matches(const std::string &type_2, const std::string &type_3) const {
	 bool r = false;
	 if (atom_type_2 == type_2)
	    if (atom_type_3 == type_3)
	       r = true;
	 return r;
      }
      friend std::ostream& operator<<(std::ostream &s, const energy_lib_torsion &torsion);
   };
   std::ostream& operator<<(std::ostream &s, const energy_lib_torsion &torsion);

   // --------------------------
   // energy container
   // --------------------------
   // 
   class energy_lib_t {

      // so that we can return an angle, a status and a message.  We
      // don't want to keep calling get_angle when the first time we
      // get a ENERGY_TYPES_NOT_FOUND.
      // 
      class energy_angle_info_t {
      public:
	 enum { OK, ANGLE_NOT_FOUND, ENERGY_TYPES_NOT_FOUND};
	 short int status;
	 energy_lib_angle angle;
	 std::string message;
	 energy_angle_info_t() {
	    status = ANGLE_NOT_FOUND;
	 }
	 energy_angle_info_t(short int status, const energy_lib_angle &angle, std::string message);
      };

      class energy_torsion_info_t {
      public:
	 enum { OK, TORSION_NOT_FOUND, ENERGY_TYPES_NOT_FOUND};
	 short int status;
	 energy_lib_angle angle;
	 std::string message;
	 energy_torsion_info_t() {
	    status = TORSION_NOT_FOUND;
	 }
	 energy_torsion_info_t(short int status,
			       const energy_lib_torsion &torsion,
			       std::string message);
      };
      
      // if permissive is true, allow the bond to be matched by
      // default/"" energy type.  Order dependent.
      energy_lib_bond get_bond(const std::string &atom_type_1,
			       const std::string &atom_type_2,
			       const std::string &bond_type, // refmac energy lib format
			       bool permissive) const;

      // if permissive is true, allow the bond to be matched by
      // default/"" energy type.  Order dependent.
      // 
      energy_angle_info_t get_angle(const std::string &atom_type_1,
				    const std::string &atom_type_2,
				    const std::string &atom_type_3,
				    bool permissive_atom_2,
				    bool permissive_atom_3) const;

      
   public:
      std::map<std::string, energy_lib_atom> atom_map;
      std::vector<energy_lib_bond> bonds;
      std::vector<energy_lib_angle> angles;
      std::vector<energy_lib_torsion> torsions;
      
      energy_lib_t() {}
      energy_lib_t(const std::string &file_name) { read(file_name); }

      // Will throw an std::runtime_error if not found.
      // 
      energy_lib_bond get_bond(const std::string &atom_type_1,
			       const std::string &atom_type_2,
			       const std::string &bond_type) const; // refmac energy lib format 
      // 
      energy_lib_angle get_angle(const std::string &atom_type_1,
				 const std::string &atom_type_2,
				 const std::string &atom_type_3) const;

      // types of the 2 middle atoms
      // Will throw an std::runtime_error if not found.
      energy_lib_torsion get_torsion(const std::string &atom_type_2,
				     const std::string &atom_type_3) const;
      
      void read(const std::string &file_name);
      void add_energy_lib_atom(    const energy_lib_atom    &atom);
      void add_energy_lib_bond(    const energy_lib_bond    &bond);
      void add_energy_lib_angle(   const energy_lib_angle   &angle);
      void add_energy_lib_torsion(const energy_lib_torsion &torsion);
      void add_energy_lib_atoms( PCMMCIFLoop mmCIFLoop);
      void add_energy_lib_bonds( PCMMCIFLoop mmCIFLoop);
      void add_energy_lib_angles(PCMMCIFLoop mmCIFLoop);
      void add_energy_lib_torsions(PCMMCIFLoop mmCIFLoop);
   };


   // ---------------------------------------------------------------
   // helper classes for linkage selection
   // ---------------------------------------------------------------

   class glycosidic_distance {
   public:
      double distance;
      CAtom *at1;
      CAtom *at2;
      glycosidic_distance(CAtom *at1_in, CAtom *at2_in, double d) {
	 at1 = at1_in;
	 at2 = at2_in;
	 distance = d;
      }
      bool operator<(const glycosidic_distance &d1) const {
	 if (d1.distance < distance)
	    return 1;
	 else
	    return 0;
      }
   };

   

   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------
   // class protein_geometry     the container class
   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------
   // 
   // consider molecule_geometry
   class protein_geometry {
      
      PCSBase SBase; 

      enum { MON_LIB_LIST_CIF = -999}; // A negative number special
				       // flag to tell the reader that
				       // this is not a normal
				       // residue's restraints and
				       // that the chem_comp should
				       // not be added to the list
				       // (currently).

      enum { UNSET_NUMBER = -1 };  // An unset number, for example the
      // number of atoms.

      bool close_float_p (const realtype &f1, const realtype &f2) { //testing func
	 float d = fabsf(f1-f2);
	 if (d < 0.001)
	    return true;
	 else
	    return false;
      }
      
      // std::vector<simple_residue_t> residue; 
      std::vector<std::string> residue_codes;
      bool verbose_mode;

      std::vector<dictionary_residue_restraints_t> dict_res_restraints;
      std::vector<dictionary_residue_link_restraints_t> dict_link_res_restraints;
      std::vector<chem_link> chem_link_vec;
      std::vector<list_chem_mod>  chem_mod_vec;

      // the monomer data in list/mon_lib_list.cif, not the
      // restraints, just id, 3-letter-code, name, group,
      // number-of-atoms, description_level.
      // Added to by the simple_mon_lib* functions.
      //
      // Use a map for faster lookups.  the key is the comp_id;
      // 
      std::map<std::string,dictionary_residue_restraints_t> simple_monomer_descriptions;

      int  comp_atom   (PCMMCIFLoop mmCIFLoop); 
      std::string comp_atom_pad_atom_name(const std::string &atom_id, const std::string &type_symbol) const;
      // return the comp_id
      std::string chem_comp   (PCMMCIFLoop mmCIFLoop);
      void comp_tree   (PCMMCIFLoop mmCIFLoop); 
      int  comp_bond   (PCMMCIFLoop mmCIFLoop); 
      void comp_angle  (PCMMCIFLoop mmCIFLoop); 
      void comp_torsion(PCMMCIFLoop mmCIFLoop); 
      void comp_plane  (PCMMCIFLoop mmCIFLoop); 
      std::pair<int, std::vector<std::string> >
      comp_chiral (PCMMCIFLoop mmCIFLoop);  // return the number of chirals and a vector
						 // of monomer names that have had
						 // chirals added (almost certainly just
						 // one of them, of course).

      void add_chem_links (PCMMCIFLoop mmCIFLoop); // references to the modifications
                                                // to the link groups (the modifications
                                                // themselves are in data_mod_list)
      int  link_bond   (PCMMCIFLoop mmCIFLoop); 
      void link_angle  (PCMMCIFLoop mmCIFLoop); 
      void link_torsion(PCMMCIFLoop mmCIFLoop); 
      void link_plane  (PCMMCIFLoop mmCIFLoop);
      int  link_chiral  (PCMMCIFLoop mmCIFLoop); // return number of new chirals
      void pdbx_chem_comp_descriptor(PCMMCIFLoop mmCIFLoop); 

      // return the comp id (so that later we can associate the file name with the comp_id).
      // 
      std::string chem_comp_component(PCMMCIFStruct structure);
      // non-looping (single) tor
      void chem_comp_tor_structure(PCMMCIFStruct structure);
      // non-looping (single) chir
      void chem_comp_chir_structure(PCMMCIFStruct structure);


      void mon_lib_add_chem_comp(const std::string &comp_id,
				 const std::string &three_letter_code,
				 const std::string &name,
				 const std::string &group,
				 int number_atoms_all, int number_atoms_nh,
				 const std::string &description_level);

      void mon_lib_add_atom(const std::string &comp_id,
			    const std::string &atom_id,
			    const std::string &atom_id_4c,
			    const std::string &type_symbol,
			    const std::string &type_energy,
			    const std::pair<bool, realtype> &partial_charge,
			    const std::pair<bool, clipper::Coord_orth> &model_pos,
			    const std::pair<bool, clipper::Coord_orth> &model_pos_ideal);

      // called because they were all at origin, for example.
      void delete_atom_positions(const std::string &comp_id, int pos_type);
			    
      void mon_lib_add_tree(std::string comp_id,
			    std::string atom_id,
			    std::string atom_back,
			    std::string atom_forward,
			    std::string connect_type);
   
      void mon_lib_add_bond(std::string comp_id,
			    std::string atom_id_1, std::string atom_id_2,
			    std::string type,
			    realtype value_dist, realtype value_dist_esd);

      void mon_lib_add_bond_no_target_geom(std::string comp_id,
					   std::string atom_id_1, std::string atom_id_2,
					   std::string type);
      
      void mon_lib_add_angle(std::string comp_id,
			     std::string atom_id_1,
			     std::string atom_id_2,
			     std::string atom_id_3,
			     realtype value_angle, realtype value_angle_esd);

      void mon_lib_add_torsion(std::string comp_id,
			       std::string torsion_id,
			       std::string atom_id_1,
			       std::string atom_id_2,
			       std::string atom_id_3,
			       std::string atom_id_4,
			       realtype value_angle, realtype value_angle_esd,
			       int period);

      void mon_lib_add_chiral(std::string comp_id,
			      std::string id,
			      std::string atom_id_centre,
			      std::string atom_id_1,
			      std::string atom_id_2,
			      std::string atom_id_3,
			      std::string volume_sign); 

      // Add a plane atom, we look through restraints trying to find a
      // comp_id, and then a plane_id, if we find it, simply push back
      // the atom name, if not, we create a restraint.
      // 
      void mon_lib_add_plane(const std::string &comp_id,
			     const std::string &plane_id,
			     const std::string &atom_id,
			     const realtype &dist_esd);

      void add_restraint(std::string comp_id, const dict_bond_restraint_t &restr);
      void add_restraint(std::string comp_id, const dict_angle_restraint_t &restr);
      void add_restraint(std::string comp_id, const dict_torsion_restraint_t &restr);
      void add_restraint(std::string comp_id, const dict_chiral_restraint_t &rest);

      void add_pdbx_descriptor(const std::string &comp_id,
			       pdbx_chem_comp_descriptor_item &descr);

      // comp_tree need to convert unquoted atom_ids to 4char atom
      // names.  So we look them up in the atom table.  If name not
      // found, return the input string
      //
      std::string atom_name_for_tree_4c(const std::string &comp_id, const std::string &atom_id) const;

      // for simple monomer descriptions:

      // return the comp id (so that later we can associate the file name with the comp_id).
      // 
      std::string simple_mon_lib_chem_comp   (PCMMCIFLoop mmCIFLoop);
      // add to simple_monomer_descriptions not dict_res_restraints.


      void add_cif_file_name(const std::string &cif_filename,
			     const std::string &comp_id1,
			     const std::string &comp_id2);

      void simple_mon_lib_add_chem_comp(const std::string &comp_id,
					const std::string &three_letter_code,
					const std::string &name,
					const std::string &group,
					int number_atoms_all, int number_atoms_nh,
					const std::string &description_level);

      // mod stuff (references by chem links)

      int add_chem_mods(PCMMCIFData data);
      int add_chem_mod(PCMMCIFLoop mmCIFLoop);
      int add_mod(PCMMCIFData data);

      int add_mods(PCMMCIFData data);
      int add_chem_mods(PCMMCIFLoop mmCIFLoop);
      
      // which calls:
      void add_chem_mod_atom( PCMMCIFLoop mmCIFLoop);
      void add_chem_mod_bond( PCMMCIFLoop mmCIFLoop);
      void add_chem_mod_tree( PCMMCIFLoop mmCIFLoop);
      void add_chem_mod_angle(PCMMCIFLoop mmCIFLoop);
      void add_chem_mod_tor(  PCMMCIFLoop mmCIFLoop);
      void add_chem_mod_chir( PCMMCIFLoop mmCIFLoop);
      void add_chem_mod_plane(PCMMCIFLoop mmCIFLoop);


      // synonyms (for RNA/DNA)
      // 
      void add_synonyms(PCMMCIFData data);
      void add_chem_comp_synonym(PCMMCIFLoop mmCIFLoop);
      class residue_name_synonym {
      public:
	 residue_name_synonym(std::string &comp_id_in,
			      std::string &comp_alternative_id_in,
			      std::string &mod_id_in) {
	    comp_id = comp_id_in;
	    comp_alternative_id = comp_alternative_id_in;
	    mod_id = mod_id_in;
	 }
	 std::string comp_id;
	 std::string comp_alternative_id;
	 std::string mod_id;
      };
      std::vector<residue_name_synonym> residue_name_synonyms;

      // link stuff
      int init_links(PCMMCIFData data);

      void link_add_bond(const std::string &link_id,
			 int atom_1_comp_id,
			 int atom_2_comp_id,
			 const std::string &atom_id_1,
			 const std::string &atom_id_2,
			 realtype value_dist,
			 realtype value_dist_esd);
      
      void link_add_angle(const std::string &link_id,
			  int atom_1_comp_id,
			  int atom_2_comp_id,
			  int atom_3_comp_id,
			  const std::string &atom_id_1,
			  const std::string &atom_id_2,
			  const std::string &atom_id_3,
			  realtype value_dist,
			  realtype value_dist_esd);

      // we want to allow synthetic/programatic addition of link
      // torsion restraints (so we can then know the rotatable bonds
      // in a link) so this should be public?

      void link_add_torsion(const std::string &link_id,
			    int atom_1_comp_id,
			    int atom_2_comp_id,
			    int atom_3_comp_id,
			    int atom_4_comp_id,
			    const std::string &atom_id_1,
			    const std::string &atom_id_2,
			    const std::string &atom_id_3,
			    const std::string &atom_id_4,
			    realtype value_dist,
			    realtype value_dist_esd,
			    int period,
			    const std::string &id); // psi, phi (or carbo link id)

      void link_add_chiral(const std::string &link_id,
			   int atom_c_comp_id,
			   int atom_1_comp_id,
			   int atom_2_comp_id,
			   int atom_3_comp_id,
			   const std::string &atom_id_c,
			   const std::string &atom_id_1,
			   const std::string &atom_id_2,
			   const std::string &atom_id_3,
			   int volume_sign);

      void link_add_plane(const std::string &link_id,
			  const std::string &atom_id,
			  const std::string &plane_id,
			  int atom_comp_id,
			  double dist_esd);

      void assign_chiral_volume_targets();

      // the chiral restraint for this comp_id(s) may need filtering
      // (i.e. removing some of them if they are not real chiral centres
      // (e.g. from prodrg restraints)).
      // 
      void filter_chiral_centres(const std::vector<std::string> & comp_id_for_filtering);

      // Return a filtered list, that is don't include chiral centers that
      // are connected to more than one hydrogen.
      // 
      std::vector<dict_chiral_restraint_t> filter_chiral_centres(const dictionary_residue_restraints_t &restraints);
      
      void assign_link_chiral_volume_targets();
      int read_number;

      std::vector <std::string> standard_protein_monomer_files() const; 
      // a wrapper to init_refmac_mon_lib
      int refmac_monomer(const std::string &s, // dir
			 const std::string &protein_mono); // extra path to file

      // not const because we can do a dynamic add.
      int get_monomer_type_index(const std::string &monomer_type);
      std::string get_padded_name(const std::string &atom_id, const int &comp_id_index) const;

      // return a list of torsion restraints that are unique for atoms 2 and
      // 3 (input restraints vector can potentially have many restraints
      // that have the same atoms 2 and 3).
      //
      std::vector <dict_torsion_restraint_t>
      filter_torsion_restraints(const std::vector <dict_torsion_restraint_t> &restraints_in) const;
      static bool torsion_restraints_comparer(const dict_torsion_restraint_t &a, const dict_torsion_restraint_t &b);

      // bool is the need-order-switch-flag
      std::vector<std::pair<chem_link, bool> > matching_chem_link(const std::string &comp_id_1,
								  const std::string &group_1,
								  const std::string &comp_id_2,
								  const std::string &group_2,
								  bool allow_peptide_link_flag) const;

      int init_refmac_mon_lib_inner(std::string filename, int read_number_in);


      energy_lib_t energy_lib;

      std::pair<bool, dictionary_residue_restraints_t>
      get_monomer_restraints_internal(const std::string &monomer_type, bool allow_minimal_flag) const;

      // created for looking up energy types cheaply (hopefully).
      // return -1 on restraints not found.
      // 
      int get_monomer_restraints_index(const std::string &monomer_type, bool allow_minimal_flag) const;

      std::vector<std::string> non_auto_load_residue_names;
      bool is_non_auto_load_ligand(const std::string resname) const;
      void fill_default_non_auto_load_residue_names(); // build-it defaults

   public:

      protein_geometry() {
	 read_number = 0;
	 set_verbose(1);
	 SBase = NULL;
	 fill_default_non_auto_load_residue_names();
      }
      
      // SBase things

      // init_sbase() should return SBase_OK.  If not, trouble.
      // 
      int init_sbase(const std::string &sbase_monomer_dir); // inits SBase
      void read_sbase_residues();
      // return NULL on unable to make residue
      CResidue *get_sbase_residue(const std::string &res_name) const;

      // return a vector of compound ids and compound names for
      // compounds that contain compound_name_substring.
      //
      std::vector<std::pair<std::string, std::string> >
      matching_sbase_residues_names(const std::string &compound_name_substring) const;

      // and fills these chem mod classes, simple container class
      // indexed with map on the mod_id


      class chem_mod {

      public:
	 chem_mod() {};
	 std::vector<chem_mod_atom>  atom_mods;
	 std::vector<chem_mod_tree>  tree_mods;
	 std::vector<chem_mod_bond>  bond_mods;
	 std::vector<chem_mod_angle> angle_mods;
	 std::vector<chem_mod_tor>   tor_mods;
	 std::vector<chem_mod_plane> plane_mods;
	 std::vector<chem_mod_chir>  chir_mods;
	 void add_mod_atom(const chem_mod_atom &chem_atom) {
	    atom_mods.push_back(chem_atom);
	 }
	 void add_mod_tree(const chem_mod_tree &chem_tree) {
	    tree_mods.push_back(chem_tree);
	 }
	 void add_mod_bond(const chem_mod_bond &chem_bond) {
	    bond_mods.push_back(chem_bond);
	 }
	 void add_mod_angle(const chem_mod_angle &chem_angle) {
	    angle_mods.push_back(chem_angle);
	 }
	 void add_mod_tor(const chem_mod_tor &chem_tor) {
	    tor_mods.push_back(chem_tor);
	 }
	 void add_mod_plane(const chem_mod_plane &chem_plane) {
	    plane_mods.push_back(chem_plane);
	 }
	 void add_mod_chir(const chem_mod_chir &chem_chir) {
	    chir_mods.push_back(chem_chir);
	 }
	 chem_mod_plane &operator[](const chem_mod_plane &plane) {
	    for (unsigned int iplane=0; iplane<plane_mods.size(); iplane++) {
	       if (plane_mods[iplane].plane_id == plane.plane_id) {
		  return plane_mods[iplane];
	       }
	    }
	    plane_mods.push_back(plane);
	    return plane_mods.back();
	 }
      };
      
      
      std::map<std::string, chem_mod> mods;
      // can throw an std::runtime exception if there is no chem mod
      // for the link_id (and that's fine)
      std::pair<chem_mod, chem_mod> get_chem_mods_for_link(const std::string &link_id) const;
      void debug_mods() const;
      
      
      // Refmac monomer lib things
      // 
      // Return the number of bond restraints
      // 
      int init_refmac_mon_lib(std::string filename, int read_number_in);

      int size() const { return dict_res_restraints.size(); }
      const dictionary_residue_restraints_t & operator[](int i) const {
	 return dict_res_restraints[i]; }
      const dictionary_residue_link_restraints_t & link(int i) const {
	 return dict_link_res_restraints[i]; }
      dictionary_residue_link_restraints_t link(const std::string &id_in) const {
	 dictionary_residue_link_restraints_t r;
	 for (unsigned int id=0; id<dict_link_res_restraints.size(); id++) {
	    if (dict_link_res_restraints[id].link_id == id_in) {
	       r = dict_link_res_restraints[id];
	       break;
	    }
	 } 
	 return r;
      }

      // return "" on comp_id not found, else return the file name.
      // 
      std::string get_cif_file_name(const std::string &comp_id) const;
      
      int link_size() const { return dict_link_res_restraints.size(); }
      void info() const;
      std::string three_letter_code(const unsigned int &i) const;

      void set_verbose(bool verbose_mode_in);

      int init_standard(); // standard protein residues and links.
       			   // Return the current read_number
      
      // Return 0 on failure to do a dynamic add, otherwise return the
      // number of atoms read.
      // 
      int try_dynamic_add(const std::string &resname, int read_number);  // return success status?
                                              // this is not const if we use dynamic add.

      // return true on having deleted;
      bool delete_mon_lib(std::string comp_id); // delete comp_id from dict_res_restraints
						// (if it exists there).



      // return a pair, the first is status (1 if the name was found, 0 if not)
      // 
      std::pair<bool, std::string> get_monomer_name(const std::string &comp_id) const;
      
      std::vector <dict_torsion_restraint_t>
      get_monomer_torsions_from_geometry(const std::string &monomer_type);
      std::vector <dict_chiral_restraint_t>
      get_monomer_chiral_volumes(const std::string monomer_type) const;

      // as above, except filter out of the returned vectors torsions
      // that move (or are based on) hydrogens.
      // 
      std::vector <dict_torsion_restraint_t>
      get_monomer_torsions_from_geometry(const std::string &monomer_type, 
					 bool find_hydrogen_torsions) const;

      // Return success status in first (0 is fail) and the second is
      // a whole residue's restraints so that we can use it to test if
      // an atom is a hydrogen.  Must have a full entry (not from
      // sbase) for the first of the returned pair to be true.
      // 
      std::pair<bool, dictionary_residue_restraints_t>
      get_monomer_restraints(const std::string &monomer_type) const;

      // Return success status in first (0 is fail) and the second is
      // a whole residue's restraints so that we can use it to test if
      // an atom is a hydrogen.
      //
      // the dictionary_residue_restraints_t is returned even we have
      // a minimal restraints description (e.g. from sbase).
      // 
      std::pair<bool, dictionary_residue_restraints_t>
      get_monomer_restraints_at_least_minimal(const std::string &monomer_type) const;

      

      // non-const because we try to read in stuff from SBase when
      // it's not in the dictionary yet.  SBase gives us bond orders.
      //
      // This relies on SBase being setup before we get to make this
      // call (init_sbase()).
      // 
      std::pair<bool, dictionary_residue_restraints_t>
      get_bond_orders(const std::string &monomer_type);

      // used to created data from sbase to put into protein_geometry
      // object.
      //
      // return a flag to signify success.
      //
      bool fill_using_sbase(const std::string &monomer_type);

      // If monomer_type is not in dict_res_restraints, then add a new
      // item to the dict_res_restraints and add mon_res_in.  Return 1
      // for replaced 0 for added.
      // 
      bool replace_monomer_restraints(std::string monomer_type,
				      const dictionary_residue_restraints_t &mon_res_in);

      // Keep everything that we have already, replace only those
      // parts that are in mon_res_in.  If there is not already
      // something there then do nothing (because there are tree and
      // atom and comp_id info) missing from Mogul information.

      // Used to update bond and angle restraints from Mogul.
      //
      // status returned was if there was something already there.
      // 
      bool replace_monomer_restraints_conservatively(std::string monomer_type,
						     const dictionary_residue_restraints_t &mon_res_in);
      void replace_monomer_restraints_conservatively_bonds(int irest,
							   const dictionary_residue_restraints_t &mon_res);
      void replace_monomer_restraints_conservatively_angles(int irest,
							    const dictionary_residue_restraints_t &mon_res);


      // this function is no longer const because it can run try_dynamic_add
      //
      bool have_dictionary_for_residue_type(const std::string &monomer_type,
					    int read_number);

      // this is const because there is no dynamic add.
      //
      // if there is just an sbase entry, then this returns false.
      // 
      bool have_dictionary_for_residue_type_no_dynamic_add(const std::string &monomer_type) const;
      
      // this is const because there is no dynamic add.
      //
      // if there is (even ) an sbase entry, then this returns true.
      // 
      bool have_at_least_minimal_dictionary_for_residue_type(const std::string &monomer_type) const;
      
      // likewise not const
      bool have_dictionary_for_residue_types(const std::vector<std::string> &residue_types);


      // return a pair, overall status, and pair of residue names and
      // atom names that dont't match.
      //
      std::pair<bool, std::vector<std::pair<std::string, std::vector<std::string> > > >
      atoms_match_dictionary(const std::vector<CResidue *> residues,
			     bool check_hydrogens_too_flag) const;

      std::pair<bool, std::vector<std::string> >
      atoms_match_dictionary(CResidue *res, bool check_hydrogens_too_flag) const;

      // return a pair: a status, yes/no atoms match and a vector of
      // atoms whose names do not match.
      // 
      std::pair<bool, std::vector<std::string> >
      atoms_match_dictionary(CResidue *res,
			     bool check_hydrogens_too_flag,
			     const dictionary_residue_restraints_t &restraints) const;

      // add "synthetic" 5 atom planar peptide restraint
      void add_planar_peptide_restraint();
      void remove_planar_peptide_restraint();
      bool planar_peptide_restraint_state() const;
      
      // restraints for omega for both CIS and TRANS links (and
      // PTRANS)
      void add_omega_peptide_restraints();
      void remove_omega_peptide_restraints();

      // a list of comp_ids that match the string in the chem_comp
      // name using the simple_monomer_descriptions, return the
      // comp_id and name:
      std::vector<std::pair<std::string, std::string> > matching_names(const std::string &test_string, short int allow_minimal_descriptions) const;

      // make a connect file specifying the bonds to Hydrogens
      bool hydrogens_connect_file(const std::string &resname,
				  const std::string &filename) const;
      
      std::vector<std::string> monomer_types() const;

      // calls try_dynamic_add if needed.
      // make HETATMs if non-standard residue name.
      CMMDBManager *mol_from_dictionary(const std::string &three_letter_code,
					bool idealised_flag);
      
      // Used by above (or maybe you just want a residue?)
      // (Can return NULL).
      // 
      // Something strange happens with internal-to-a-CResidue atom
      // indexing when I tried to use this.  The problem was resoloved
      // by using mol_from_dictionary() above and getting the first
      // residue.  Something to do with atom indexing on checking
      // in...?
      // 
      CResidue *get_residue(const std::string &comp_id, bool idealised_flag);

      // Thow a std::runtime_error exception if we can't get the group of r
      std::string get_group(CResidue *r) const;
      // and the string version of this
      std::string get_group(const std::string &res_name) const;

      // bool is the need-order-switch-flag
      std::vector<std::pair<chem_link, bool> >
      matching_chem_link(const std::string &comp_id_1,
			 const std::string &group_1,
			 const std::string &comp_id_2,
			 const std::string &group_2) const;
      
      // Try to find a link that is not a peptide link (because that
      // fails on a distance check).  This is the method to find
      // isopeptide links (which again need to be distance checked in
      // find_link_type_rigourous()).
      // 
      // bool the need-order-switch-flag
      std::vector<std::pair<chem_link, bool> >
      matching_chem_link_non_peptide(const std::string &comp_id_1,
				     const std::string &group_1,
				     const std::string &comp_id_2,
				     const std::string &group_2) const;

      // return "" on failure.
      // no order switch is considered.
      // 
      std::string find_glycosidic_linkage_type(CResidue *first, CResidue *second) const;

      std::pair<std::string, bool>
      find_glycosidic_linkage_type_with_order_switch(CResidue *first, CResidue *second) const;

      // can throw a std::runtime_error
      chem_link get_chem_link(const std::string &link_id);

      void print_chem_links() const;
      static int chiral_volume_string_to_chiral_sign(const std::string &chiral_vol_string);
      static std::string make_chiral_volume_string(int chiral_sign);

      bool linkable_residue_types_p(const std::string &this_res_type,
				    const std::string &env_residue_res_type);

      bool OXT_in_residue_restraints_p(const std::string &residue_type) const;

      void read_energy_lib(const std::string &file_name);

      // return HB_UNASSIGNED when not found
      // 
      int get_h_bond_type(const std::string &atom_name, const std::string &monomer_name) const;

      // Find the bonded neighbours of the given atoms - throw an
      // exception if residue name is not in dictionary.
      // 
      std::vector<std::string> get_bonded_neighbours(const std::string &residue_name,
						     const std::string &atom_name_1,
						     const std::string &atom_name_2,
						     bool also_2nd_order_neighbs_flag=0) const;

      // add a monomer restraints description.
      //
      void add(const dictionary_residue_restraints_t &rest) {
	 dict_res_restraints.push_back(rest);
      }

      // a new pdb file has been read in (say).  The residue types
      // have been compared to the dictionary.  These (comp_ids) are
      // the types that are not in the dictionary.  Try to load an
      // sbase description at least so that we can draw their bonds
      // correctly.  Use fill_using_sbase().
      // 
      bool try_load_sbase_description(const std::vector<std::string> &comp_ids);

      // expand to 4c, the atom_id, give that it should match an atom of the type comp_id.
      // Used in chem mods, when we don't know the comp_id until residue modification time.
      // 
      std::string atom_id_expand(const std::string &atom_id,
				 const std::string &comp_id) const;

      // return "" if not found, else return the energy type found in ener_lib.cif
      // 
      std::string get_type_energy(const std::string &atom_name,
				  const std::string &residue_name) const;

      // return -1.1 on failure to look up.
      // 
      double get_vdw_radius(const std::string &atom_name,
			    const std::string &residue_name,
			    bool use_vdwH_flag) const;


      // Find the non-bonded contact distance
      // 
      // Return a pair, if not found the first is 0.  Look up in the energy_lib.
      // 
      std::pair<bool, double> get_nbc_dist(const std::string &energy_type_1,
					   const std::string &energy_type_2) const;


      // This uses have_dictionary_for_residue_type() (and thus
      // try_dynamic_add() if needed).
      // 
      // Return -1 if residue type not found.
      // 
      int n_hydrogens(const std::string &residue_type);

      // Add XXX or whatever to non-auto residue names
      // 
      void add_non_auto_load_residue_name(const std::string &res_name);
      void remove_non_auto_load_residue_name(const std::string &res_name);

      std::vector<std::string> monomer_restraints_comp_ids() const;

      // can throw a std::runtime_error
      std::string Get_SMILES_for_comp_id(const std::string &comp_id) const;

   };

} // namespace coot

#endif //  PROTEIN_GEOMETRY_HH
