/* geometry/protein-geometry.cc
 * 
 * Copyright 2004, 2005 The University of York
 * Copyright 2008, 2009 The University of Oxford
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

#define MONOMER_DIR_STR "SBASE_DIR"

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#include "mmdb_mmcif.h"

#define USE_SBASE

#ifdef USE_SBASE 
// needs #include "mmdb_sbase.h"
#ifndef  __MMDB_SBase__
#include "mmdb_sbase.h"
#endif
#endif // USE_SBASE

#include "clipper/core/coords.h"

// #include "db-main.h"

namespace coot {

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
   };


   class basic_dict_restraint_t {
      std::string atom_id_1_;
      std::string atom_id_2_;

   protected:
      std::string atom_id_mmdb_expand(const std::string &atomname) const {
	 std::string r;
	 int ilen = atomname.length();
	 
	 if (ilen == 4) return atomname;
	 
	 if (ilen == 1) {
	    r = " ";
	    r += atomname;
	    r += "  ";
	 } else {
	    if (ilen == 2) { 
	       r = " ";
	       r += atomname;
	       r += " ";
	    } else {
	       if (ilen == 3) {
		  r = " ";
		  r += atomname;
	       } else {
		  r = atomname;
	       }
	    }
	 }
	 return r;
      }

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
      }
      std::string type() const { return type_; }
      double dist() const { return dist_; }
      double esd () const { return dist_esd_;}
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
   };

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
      enum { CHIRAL_RESTRAINT_BOTH = -2};
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
   public:
      dictionary_residue_restraints_t(std::string comp_id_in,
				      int read_number_in) {
	 has_partial_charges_flag = 0;
	 comp_id = comp_id_in;
	 read_number = read_number_in;
      }
      void clear_dictionary_residue();
      dict_chem_comp_t residue_info;
      std::vector <dict_atom> atom_info;
      std::string comp_id; // i.e. residue type name
      std::vector <dict_chem_comp_tree_t> tree;
      int read_number; 
      std::vector    <dict_bond_restraint_t>    bond_restraint;
      std::vector   <dict_angle_restraint_t>   angle_restraint;
      std::vector <dict_torsion_restraint_t> torsion_restraint;
      std::vector  <dict_chiral_restraint_t>  chiral_restraint;
      std::vector   <dict_plane_restraint_t>   plane_restraint;
      // Return 1 for hydrogen or deuterium, 0 for not found or not a hydrogen.
      bool is_hydrogen(const std::string &atom_name) const;
      int assign_chiral_volume_targets(); // return the number of targets made.
      bool has_unassigned_chiral_volumes() const;
      bool has_partial_charges() const;
      void set_has_partial_charges(bool state) {
	 has_partial_charges_flag = state;
      }
      std::vector<dict_torsion_restraint_t> get_non_const_torsions(bool include_hydrogen_torsions_flag) const;
      void write_cif(const std::string &filename) const;
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
   };
   std::ostream& operator<<(std::ostream &s, chem_link lnk);

   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------
   // class protein_geometry     the container class
   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------
   // 
   // consider molecule_geometry
   class protein_geometry {
      
#ifdef USE_SBASE   
      PCSBase SBase; 
#endif // USE_SBASE

      enum { MON_LIB_LIST_CIF = -999}; // A negative number special
				       // flag to tell the reader that
				       // this is not a normal
				       // residue's restraints and
				       // that the chem_comp should
				       // not be added to the list
				       // (currently).

      enum { UNSET_NUMBER = -1 };  // An unset number, for example the
				  // number of atoms.
      
      // std::vector<simple_residue_t> residue; 
      std::vector<std::string> residue_codes;
      bool verbose_mode;

      std::vector<dictionary_residue_restraints_t> dict_res_restraints;
      std::vector<dictionary_residue_link_restraints_t> dict_link_res_restraints;
      std::vector<chem_link> chem_link_vec;

      // the monomer data in list/mon_lib_list.cif, not the
      // restraints, just id, 3-letter-code, name, group,
      // number-of-atoms, description_level.
      // Added to by the simple_mon_lib* functions.
      std::vector<dictionary_residue_restraints_t> simple_monomer_descriptions;

      int  comp_atom   (PCMMCIFLoop mmCIFLoop); 
      std::string comp_atom_pad_atom_name(const std::string &atom_id, const std::string &type_symbol) const;
      void chem_comp   (PCMMCIFLoop mmCIFLoop);
      void comp_tree   (PCMMCIFLoop mmCIFLoop); 
      int  comp_bond   (PCMMCIFLoop mmCIFLoop); 
      void comp_angle  (PCMMCIFLoop mmCIFLoop); 
      void comp_torsion(PCMMCIFLoop mmCIFLoop); 
      int  comp_chiral (PCMMCIFLoop mmCIFLoop);  // return number of chirals.
      void comp_plane  (PCMMCIFLoop mmCIFLoop); 

      void add_chem_links (PCMMCIFLoop mmCIFLoop); // references to the modifications
                                                // to the link groups (the modifications
                                                // themselves are in data_mod_list)
      int  link_bond   (PCMMCIFLoop mmCIFLoop); 
      void link_angle  (PCMMCIFLoop mmCIFLoop); 
      void link_torsion(PCMMCIFLoop mmCIFLoop); 
      void link_plane  (PCMMCIFLoop mmCIFLoop);
      int  link_chiral  (PCMMCIFLoop mmCIFLoop); // return number of new chirals

      void chem_comp_component(PCMMCIFStruct structure);


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
			    

      void mon_lib_add_tree(std::string comp_id,
			    std::string atom_id,
			    std::string atom_back,
			    std::string atom_forward,
			    std::string connect_type);
   
      void mon_lib_add_bond(std::string comp_id,
			    std::string atom_id_1, std::string atom_id_2,
			    std::string type,
			    realtype value_dist, realtype value_dist_esd);

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

      // for simple monomer descriptions:
      void simple_mon_lib_chem_comp   (PCMMCIFLoop mmCIFLoop);
      // add to simple_monomer_descriptions not dict_res_restraints.
      void simple_mon_lib_add_chem_comp(const std::string &comp_id,
					const std::string &three_letter_code,
					const std::string &name,
					const std::string &group,
					int number_atoms_all, int number_atoms_nh,
					const std::string &description_level);

      // mod stuff (references by chem links)
      int add_mods(PCMMCIFData data);
      int add_chem_mods(PCMMCIFLoop mmCIFLoop); 

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
      void assign_link_chiral_volume_targets();
      int read_number;

      void delete_mon_lib(std::string comp_id); // delete comp_id from dict_res_restraints
						// (if it exists there).

      std::vector <std::string> standard_protein_monomer_files() const; 
      // a wrapper to init_refmac_mon_lib
      int refmac_monomer(const std::string &s, // dir
			 const std::string &protein_mono); // extra path to file

      // not const because we can do a dynamic add.
      int get_monomer_type_index(const std::string &monomer_type);
      std::string get_padded_name(const std::string &atom_id, const int &comp_id_index) const;
      std::vector <coot::dict_torsion_restraint_t>
      filter_torsion_restraints(const std::vector <coot::dict_torsion_restraint_t> &restraints_in) const;
      static bool torsion_restraints_comparer(const coot::dict_torsion_restraint_t &a, const coot::dict_torsion_restraint_t &b);

      // bool is the need-order-switch-flag
      std::vector<std::pair<chem_link, bool> > matching_chem_link(const std::string &comp_id_1,
								  const std::string &group_1,
								  const std::string &comp_id_2,
								  const std::string &group_2,
								  bool allow_peptide_link_flag) const;


   public:

      protein_geometry() { read_number = 0; set_verbose(1); }
#ifdef USE_SBASE   
      // SBase things
      int init_sbase(const std::string &sbase_monomer_dir); // inits SBase
      void read_sbase_residues();
#endif // USE_SBASE   
      
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
      int link_size() const { return dict_link_res_restraints.size(); }
      void info() const;
      std::string three_letter_code(const unsigned int &i) const;

      void set_verbose(bool verbose_mode_in);

      int init_standard(); // standard protein residues and links.
       			   // Return the current read_number
      
      // Return 0 on failure to do a dynamic add 
      // 
      int try_dynamic_add(const std::string &resname, int read_number);  // return success status?
      // this is not const if we use dynamic add.
      std::vector <coot::dict_torsion_restraint_t>
      get_monomer_torsions_from_geometry(const std::string &monomer_type);
      std::vector <coot::dict_chiral_restraint_t>
      get_monomer_chiral_volumes(const std::string monomer_type) const;

      // as above, except filter out of the returned vectors torsions
      // that move (or are based on) hydrogens.
      // 
      std::vector <coot::dict_torsion_restraint_t>
      get_monomer_torsions_from_geometry(const std::string &monomer_type, 
					 short int find_hydrogen_torsions) const;

      // Return success status in first (0 is fail) and the second is
      // a whole residue's restraints so that we can use it to test if
      // an atom is a hydrogen.
      // 
      std::pair<short int, dictionary_residue_restraints_t>
      get_monomer_restraints(const std::string &monomer_type) const;

      // If monomer_type is not in dict_res_restraints, then add a new
      // item to the dict_res_restraints and add mon_res_in.  Return 1
      // for replaced 0 for added.
      // 
      bool replace_monomer_restraints(std::string monomer_type,
				      const dictionary_residue_restraints_t &mon_res_in);

      // this function is no longer const because it can run try_dynamic_add
      //
      int have_dictionary_for_residue_type(const std::string &monomer_type,
					   int read_number);
      // likewise not const
      bool have_dictionary_for_residue_types(const std::vector<std::string> &residue_types);

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
      CMMDBManager *mol_from_dictionary(const std::string &three_letter_code,
					bool idealised_flag) const;

      // Thow an exception if we can't get the group of r
      std::string get_group(CResidue *r) const;

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

      void print_chem_links() const;
      static int chiral_volume_string_to_chiral_sign(const std::string &chiral_vol_string);
      static std::string make_chiral_volume_string(int chiral_sign);

      bool linkable_residue_types_p(const std::string &this_res_type,
				    const std::string &env_residue_res_type);

   };

} // namespace coot

#endif //  PROTEIN_GEOMETRY_HH
