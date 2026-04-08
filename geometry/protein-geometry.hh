/* geometry/protein-geometry.cc
 * 
 * Copyright 2004, 2005 The University of York
 * Copyright 2008, 2009, 2011, 2012 The University of Oxford
 * Copyright 2013, 2014, 2015 by Medical Research Council
 * Author: Paul Emsley
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifndef PROTEIN_GEOMETRY_HH
#define PROTEIN_GEOMETRY_HH

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <set>

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#include <map>
#include <stdexcept>
#include <algorithm>

#include <mmdb2/mmdb_utils.h>
#include <mmdb2/mmdb_math_graph.h>

#ifdef HAVE_CCP4SRS
#ifndef CCP4SRS_BASE_H
#include <ccp4srs/ccp4srs_manager.h>
#endif
#endif

#include "chem_mods.hh"

#include "match-results.hh"

#include "energy-lib.hh"

#include "clipper/core/coords.h"

#include "hb-types.hh"
#include "mini-mol/atom-quads.hh"

#include "metal-ligand.hh"

#include "pdbe-chem-comp-atom-depiction.hh"
#include "gphl-chem-comp-info.hh"
#include "pdbx-chem-comp-description-generator.hh"

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
                                     const std::string &descriptor_in) :
         type(type_in), program(program_in), program_version(program_version_in),
         descriptor(descriptor_in) { }
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
      std::string group; // e.g. "L-peptide" (pdbx calls this "type")
      std::string pdbx_type; // e.g. HETAIN
      std::string formula;
      std::string mon_nstd_parent_comp_id;
      std::string pdbx_synonyms;
      std::string pdbx_initial_date;
      std::string pdbx_modified_date;
      std::string pdbx_ambiguous_flag;
      std::string pdbx_release_status;
      std::string pdbx_replaced_by;
      std::string pdbx_replaces;
      std::string formula_weight;
      std::string one_letter_code;
      std::string pdbx_model_coordinates_missing_flag;
      std::string pdbx_ideal_coordinates_missing_flag;
      std::string pdbx_ideal_coordinates_details; // e.g. Corina
      std::string pdbx_model_coordinates_db_code;
      std::string pdbx_subcomponent_list;
      std::string pdbx_processing_site; // e.g. RCSB
      int number_atoms_all;
      int number_atoms_nh;
      std::string description_level;
      int pdb_formal_charge;
      dict_chem_comp_t(const std::string &comp_id_in,
                       const std::string &three_letter_code_in,
                       const std::string &name_in,
                       const std::string &group_in,
                       int number_atoms_all_in,
                       int number_atoms_nh_in,
                       const std::string &description_level_in) :
         comp_id(""), three_letter_code(three_letter_code_in), name(name_in), group(group_in),
         pdb_formal_charge(0) {
         setup_internal(comp_id_in,
                        three_letter_code_in, name_in, group_in,
                        number_atoms_all_in, number_atoms_nh_in,
                        description_level_in);
      }
      dict_chem_comp_t() : pdb_formal_charge(0) {
          setup_internal("", "", "", "", 0, 0, "");
      }
      dict_chem_comp_t(const dict_chem_comp_t &din) : pdb_formal_charge(0) {
         setup_internal(din.comp_id, din.three_letter_code, din.name,
                        din.group, din.number_atoms_all, din.number_atoms_nh,
                        din.description_level);
      }
      friend std::ostream& operator<<(std::ostream &s, const dict_chem_comp_t &rest);
   };
   std::ostream& operator<<(std::ostream &s, const dict_chem_comp_t &rest);


   class basic_dict_restraint_t {
      std::string atom_id_1_;
      std::string atom_id_2_;
      std::string atom_id_1_4c_;
      std::string atom_id_2_4c_;

   public:
      basic_dict_restraint_t() {} // for planes
      basic_dict_restraint_t(const std::string &at1,
                             const std::string &at2) {
         set_atom_id_1(at1);
         set_atom_id_2(at2);
      }
      std::string atom_id_1() const { return atom_id_1_;}
      std::string atom_id_2() const { return atom_id_2_;}
      std::string atom_id_1_4c() const {  // 4 character return;
         return atom_id_1_4c_;
      }
      std::string atom_id_2_4c() const {
         return atom_id_2_4c_;
      }
      void set_atom_id_1(const std::string &id) {
         atom_id_1_ = id;
         atom_id_1_4c_ = atom_id_mmdb_expand(id);
      }
      void set_atom_id_2(const std::string &id) {
         atom_id_2_ = id;
         atom_id_2_4c_ = atom_id_mmdb_expand(id);
      }
   };

   class dict_bond_restraint_t : public basic_dict_restraint_t {
      std::string type_;  // bond order
      double dist_;
      double dist_esd_;
      bool have_target_values;
      double dist_nuclear_;
      double dist_nuclear_esd_;
      bool dist_nuclear_was_set;
      bool atom_has_only_this_non_hydrogen_bond_first;
      bool atom_has_only_this_non_hydrogen_bond_second;

   public:
      enum bond_length_type_t { UNKNOWN, NUCLEAR_POSITION, ELECTONS_POSITION };
      enum aromaticity_t { NON_AROMATIC, AROMATIC, UNASSIGNED };
      aromaticity_t aromaticity;
      bond_length_type_t bond_length_type_;
      // dict_bond_restraint_t() {};
      dict_bond_restraint_t(const std::string &atom_id_1_in,
                            const std::string &atom_id_2_in,
                            const std::string &type,
                            double dist_in,
                            double dist_esd_in,
                            double dist_nuclear_in,
                            double dist_nuclear_esd_in,
                            bool dist_nuclear_was_set_in,
                            aromaticity_t arom_in=UNASSIGNED,
                            bond_length_type_t blt=UNKNOWN) :
         basic_dict_restraint_t(atom_id_1_in, atom_id_2_in), type_(type) {

         dist_ = dist_in;
         dist_esd_ = dist_esd_in;
         have_target_values = 1;
         aromaticity = arom_in;
         bond_length_type_ = blt;
	 atom_has_only_this_non_hydrogen_bond_first  = false;
	 atom_has_only_this_non_hydrogen_bond_second = false;
         dist_nuclear_ = -1.0;
         dist_nuclear_esd_ = -1.0;
         dist_nuclear_was_set = dist_nuclear_was_set_in;
         if (dist_nuclear_was_set_in) {
            dist_nuclear_     = dist_nuclear_in;
            dist_nuclear_esd_ = dist_nuclear_esd_in;
         }
      }

      dict_bond_restraint_t(const std::string &atom_id_1_in,
                            const std::string &atom_id_2_in,
                            const std::string &type,
                            aromaticity_t arom_in=UNASSIGNED) :
         basic_dict_restraint_t(atom_id_1_in, atom_id_2_in), type_(type), dist_(0.0), dist_esd_(0.0) {
         have_target_values = 0;
         aromaticity = arom_in;
         bond_length_type_ = UNKNOWN;
	 atom_has_only_this_non_hydrogen_bond_first  = false;
	 atom_has_only_this_non_hydrogen_bond_second = false;
         dist_nuclear_ = -1.0;
         dist_nuclear_esd_ = -1.0;
         dist_nuclear_was_set = false;
      }
      dict_bond_restraint_t() {} // boost::python needs this

      std::string type() const { return type_; }
      bool is_bond_to_hydrogen_atom() const;
      int mmdb_bond_type() const; // for mmdb::math::Graph mmdb::math::Edge usage
      // can throw a std::runtime_error exception (if target values not set)
      double value_dist() const {
         if (have_target_values)
            return dist_;
         else
            throw std::runtime_error("value_dist(): unset target distance");
      }
      // can throw a std::runtime_error exception
      double value_esd () const {
         if (have_target_values)
            return dist_esd_;
         else
            throw std::runtime_error("value_esd(): unset target-distance");
      }
      bool matches_names(const dict_bond_restraint_t &r) const {
         if (atom_id_1() == r.atom_id_1())
            if (atom_id_2() == r.atom_id_2())
               return true;
         if (atom_id_1() == r.atom_id_2())
            if (atom_id_2() == r.atom_id_1())
               return true;
         return false;
      }
      void set_atom_1_atom_id(const std::string &id) { set_atom_id_1(id); }
      void set_atom_2_atom_id(const std::string &id) { set_atom_id_2(id); }
      // this function called by dictionary parser after all the bonds have been read
      void set_only_bond(const std::string &pos, bool b) {
         if (pos == "first")  atom_has_only_this_non_hydrogen_bond_first  = b;
         if (pos == "second") atom_has_only_this_non_hydrogen_bond_second = b;
      }
      friend std::ostream& operator<<(std::ostream &s, const dict_bond_restraint_t &rest);
   };
   std::ostream& operator<<(std::ostream &s, const dict_bond_restraint_t &rest);

   class dict_angle_restraint_t : public basic_dict_restraint_t {
      std::string atom_id_3_;
      std::string atom_id_3_4c_;
      double angle_;
      double angle_esd_;
   public:
      // dict_angle_restraint_t() {};
      dict_angle_restraint_t(const std::string &atom_id_1,
                             const std::string &atom_id_2,
                             const std::string &atom_id_3,
                             double angle,
                             double angle_esd) :
         basic_dict_restraint_t(atom_id_1, atom_id_2), atom_id_3_(atom_id_3) {
         atom_id_3_4c_ = atom_id_mmdb_expand(atom_id_3_);
         angle_ = angle;
         angle_esd_ = angle_esd;
      };
      dict_angle_restraint_t() {} // boost::python needs this

      std::string atom_id_3() const { return atom_id_3_;}
      std::string atom_id_3_4c() const { return atom_id_3_4c_; }
      double angle() const { return angle_; }
      double esd ()  const { return angle_esd_;}

      bool matches_names(const dict_angle_restraint_t &r) const {
         if (atom_id_1() == r.atom_id_1())
            if (atom_id_2() == r.atom_id_2())
               if (atom_id_3() == r.atom_id_3())
                  return true;
         if (atom_id_1() == r.atom_id_3())
            if (atom_id_2() == r.atom_id_2())
               if (atom_id_3() == r.atom_id_1())
                  return true;
         return false;
      }
      void set_atom_1_atom_id(const std::string &id) { set_atom_id_1(id); }
      void set_atom_2_atom_id(const std::string &id) { set_atom_id_2(id); }
      void set_atom_3_atom_id(const std::string &id) { atom_id_3_ = id; }

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
      std::string atom_id_3_4c_;
      std::string atom_id_4_4c_;
      double angle_;
      double angle_esd_;
      int period;
   public:

      dict_torsion_restraint_t(const std::string &id_in,
                               const std::string &atom_id_1,
                               const std::string &atom_id_2,
                               const std::string &atom_id_3,
                               const std::string &atom_id_4,
                               double angle,
                               double angle_esd,
                               int period_in) :
         basic_dict_restraint_t(atom_id_1, atom_id_2), id_(id_in), atom_id_3_(atom_id_3), atom_id_4_(atom_id_4)
      {
         atom_id_3_4c_ = atom_id_mmdb_expand(atom_id_3_);
         atom_id_4_4c_ = atom_id_mmdb_expand(atom_id_4_);
         angle_ = angle;
         angle_esd_ = angle_esd;
         period = period_in;
      };
      std::string atom_id_3_4c() const { return atom_id_3_4c_; }
      std::string atom_id_4_4c() const { return atom_id_4_4c_; }
      std::string atom_id_3() const { return atom_id_3_;}
      std::string atom_id_4() const { return atom_id_4_;}
      std::string id() const { return id_;}
      bool is_const() const; // is the id const?  (Don't consider angle or sd).
      int periodicity() const { return period; }
      double angle() const { return angle_; }
      double esd ()  const { return angle_esd_;}
      friend std::ostream& operator<<(std::ostream &s, const dict_torsion_restraint_t &rest);
      bool is_pyranose_ring_torsion(const std::string &comp_id) const;
      bool is_ring_torsion(const std::vector<std::vector<std::string> > &ring_atoms_sets) const;
      // hack for mac, ostream problems
      bool is_peptide_torsion() const;
      std::string format() const;
      void set_atom_1_atom_id(const std::string &id) { set_atom_id_1(id); }
      void set_atom_2_atom_id(const std::string &id) { set_atom_id_2(id); }
      void set_atom_3_atom_id(const std::string &id) { atom_id_3_ = id; }
      void set_atom_4_atom_id(const std::string &id) { atom_id_4_ = id; }
   };
   std::ostream& operator<<(std::ostream &s, const dict_torsion_restraint_t &rest);

   // ------------------------------------------------------------------------
   // class dict_plane_restraint_t
   // ------------------------------------------------------------------------
   //
   class dict_plane_restraint_t : public basic_dict_restraint_t {
      std::vector<std::pair<std::string, double> > atom_ids;
      double dist_esd_;  // despite separate entries for each atom in
                        // the dictionary, I decide that it is
                        // sensible to say that all atoms in a plane
                        // have the same esd deviation from the plane.
   public:
      dict_plane_restraint_t() {};
      dict_plane_restraint_t(const std::string &plane_id_in,
                             const std::vector<std::pair<std::string, double> > &atom_esds_ins) :
         atom_ids(atom_esds_ins),
         plane_id(plane_id_in) { dist_esd_ = 0.02; }
      dict_plane_restraint_t(const std::string &plane_id_in,
                             const std::string &atom_id_in,
                             double dist_esd_in) :
         plane_id(plane_id_in) {
         atom_ids.push_back(std::pair<std::string, float> (atom_id_in, dist_esd_in));
      };
      dict_plane_restraint_t(const std::string &plane_id_in,
                             const std::vector<std::string> &plane_atom_ids,
                             double dist_esd_in) : plane_id(plane_id_in) {

         atom_ids.resize(plane_atom_ids.size());
         for (unsigned int i=0; i<plane_atom_ids.size(); i++)
            atom_ids[i] = std::pair<std::string, float> (plane_atom_ids[i], dist_esd_in);

      };
      std::string plane_id; // or int plane_id number 1, 2 3.
      double dist_esd(int i) const { return atom_ids[i].second; }
      std::string atom_id(int i) const { return atom_id_mmdb_expand(atom_ids[i].first); }
      int n_atoms() const { return atom_ids.size(); }
      bool empty() const { return (atom_ids.size() == 0); }
      std::pair<std::string, double> operator[](unsigned int i) const { return atom_ids[i];}
      bool matches_names(const dict_plane_restraint_t &r) const;
      void push_back_atom(const std::string &at, float esd) {
         atom_ids.push_back(std::pair<std::string, float> (at, esd));
      }
      friend std::ostream&  operator<<(std::ostream &s, dict_plane_restraint_t rest);
      void set_atom_ids(const std::vector<std::pair<int, std::string> > &from_tos) {
         for (unsigned int i=0; i<from_tos.size(); i++)
            atom_ids[from_tos[i].first].first = from_tos[i].second;
      }
      void set_dist_esd(unsigned int idx, double sigma_in) {
         if (idx < atom_ids.size())
            atom_ids[idx].second = sigma_in;
      }
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
             CHIRAL_VOLUME_RESTRAINT_VOLUME_SIGN_UNASSIGNED = -3,
             CHIRAL_RESTRAINT_POSITIVE = 1,
             CHIRAL_RESTRAINT_NEGATIVE = -1};
      
      dict_chiral_restraint_t() {};
      dict_chiral_restraint_t(const std::string &chiral_id_in,
                              const std::string &atom_id_centre_in,
                              const std::string &atom_id_1_in,
                              const std::string &atom_id_2_in,
                              const std::string &atom_id_3_in,
                              int volume_sign_in) :
         chiral_id(chiral_id_in), local_atom_id_centre(atom_id_centre_in),
         local_atom_id_1(atom_id_1_in),
         local_atom_id_2(atom_id_2_in),
         local_atom_id_3(atom_id_3_in)
      {
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
      std::string get_atom_id_centre() const { return local_atom_id_centre; }
      void set_atom_1_atom_id(const std::string &id) { local_atom_id_1 = id; }
      void set_atom_2_atom_id(const std::string &id) { local_atom_id_2 = id; }
      void set_atom_3_atom_id(const std::string &id) { local_atom_id_3 = id; }
      void set_atom_c_atom_id(const std::string &id) { local_atom_id_centre = id; }
      
      double assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds,
                                         const std::vector <dict_angle_restraint_t> &angles);
      double target_volume() const { return target_volume_;}
      double volume_sigma()  const { return volume_sigma_;}
      bool is_a_both_restraint() const { return is_both_flag; } 
      bool has_unassigned_chiral_volume() const {
         return (volume_sigma_ < 0.0) ? 1 : 0;
      }
      void invert_target_volume() {
         if (volume_sign == CHIRAL_RESTRAINT_POSITIVE) {
            volume_sign = CHIRAL_RESTRAINT_NEGATIVE;
            target_volume_ = -target_volume_;
         } else { 
            if (volume_sign == CHIRAL_RESTRAINT_NEGATIVE) { 
               volume_sign = CHIRAL_RESTRAINT_POSITIVE;
               target_volume_ = -target_volume_;
            }
         }
      }
      bool matches_names(const dict_chiral_restraint_t &r) const {
         if (atom_id_c_4c() != r.atom_id_c_4c()) { 
            return false;
         } else {
            if (atom_id_1_4c() == r.atom_id_1_4c())
               if (atom_id_2_4c() == r.atom_id_2_4c())
                  if (atom_id_3_4c() == r.atom_id_3_4c())
                     return true;
            if (atom_id_1_4c() == r.atom_id_2_4c())
               if (atom_id_2_4c() == r.atom_id_3_4c())
                  if (atom_id_3_4c() == r.atom_id_1_4c())
                     return true;
            if (atom_id_1_4c() == r.atom_id_3_4c())
               if (atom_id_2_4c() == r.atom_id_1_4c())
                  if (atom_id_3_4c() == r.atom_id_2_4c())
                     return true;
         }
         return false;
      }
      friend std::ostream& operator<<(std::ostream &s, const dict_chiral_restraint_t &rest);
   };
   std::ostream& operator<<(std::ostream &s, const dict_chiral_restraint_t &rest);

   // ------------------------------------------------------------------------
   // class dict_improper_dihedral_restraint_t
   // ------------------------------------------------------------------------
   //
   class dict_improper_dihedral_restraint_t : public basic_dict_restraint_t {

      // c.f. a chiral restraint, atom id 2 is the "chiral" atom (the
      // "base" of the tetrahedron).
      std::string local_atom_id_1;
      std::string local_atom_id_2;
      std::string local_atom_id_3;
      std::string local_atom_id_4;

   public:

      dict_improper_dihedral_restraint_t() {};
      dict_improper_dihedral_restraint_t(const std::string &atom_id_1_in,
                                         const std::string &atom_id_2_in,
                                         const std::string &atom_id_3_in,
                                         const std::string &atom_id_4_in) :
         local_atom_id_1(atom_id_1_in),
         local_atom_id_2(atom_id_2_in),
         local_atom_id_3(atom_id_3_in),
         local_atom_id_4(atom_id_4_in) {
         sigma  = 0.005; // will need fixing
      }
      double sigma;
      std::string atom_id_1_4c() const { return atom_id_mmdb_expand(local_atom_id_1);}
      std::string atom_id_2_4c() const { return atom_id_mmdb_expand(local_atom_id_2);}
      std::string atom_id_3_4c() const { return atom_id_mmdb_expand(local_atom_id_3);}
      std::string atom_id_4_4c() const { return atom_id_mmdb_expand(local_atom_id_4);}
      std::string get_atom_id_centre() const { return local_atom_id_2; }
      friend std::ostream& operator<<(std::ostream &s, const dict_improper_dihedral_restraint_t &rest);
   };
   std::ostream& operator<<(std::ostream &s, const dict_improper_dihedral_restraint_t &rest);


   // ------------------------------------------------------------------------
   // class dict_atom
   // ------------------------------------------------------------------------
   // 
   // one of these for each atom in a dictionary_residue_restraints_t
   // (i.e. each atom in a residue/comp_id).
   // 
   class dict_atom {
      void init() {
         aromaticity = UNASSIGNED;
         ordinal_id = -1;
         is_hydrogen_flag = false;
      }
   public:
      enum aromaticity_t { NON_AROMATIC, AROMATIC, UNASSIGNED };
      enum { IDEAL_MODEL_POS, REAL_MODEL_POS};
      std::string atom_id;
      std::string atom_id_4c;
      std::string type_symbol;
      std::string type_energy;
      std::string acedrg_atom_type;
      aromaticity_t aromaticity;
      bool is_hydrogen_flag;
      std::pair<bool, float> partial_charge;
      std::pair<bool, int> formal_charge;
      std::pair<bool, std::string> pdbx_stereo_config;
      std::pair<bool, clipper::Coord_orth> pdbx_model_Cartn_ideal;
      std::pair<bool, clipper::Coord_orth> model_Cartn;
      int ordinal_id;
      dict_atom(const std::string &atom_id_in,
                const std::string &atom_id_4c_in,
                const std::string &type_symbol_in,
                const std::string &type_energy_in,
                std::pair<bool, float> partial_charge_in) :
         atom_id(atom_id_in), atom_id_4c(atom_id_4c_in), type_symbol(type_symbol_in),
         type_energy(type_energy_in), partial_charge(partial_charge_in)
      {
         init();
         if (type_energy == "H" ) is_hydrogen_flag = true;
         if (type_symbol == "H" ) is_hydrogen_flag = true;
         if (type_symbol == " H") is_hydrogen_flag = true;
         if (type_symbol == " D") is_hydrogen_flag = true;
      }
      dict_atom() { init(); }; // for resize(0);
      void add_pos(int pos_type, const std::pair<bool, clipper::Coord_orth> &model_pos_ideal);
      void add_ordinal_id(int ordinal_id_in) { ordinal_id = ordinal_id_in; }
      bool is_hydrogen() const { return is_hydrogen_flag; }
      friend std::ostream& operator<<(std::ostream &s, const dict_atom &at);
   };
   std::ostream& operator<<(std::ostream &s, const dict_atom &at);

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
                            const std::string &atom_forward_in, const std::string &connect_type_in) :
         atom_id(atom_id_mmdb_expand(atom_id_in)),
         atom_back(atom_id_mmdb_expand(atom_back_in)),
         atom_forward(atom_id_mmdb_expand(atom_forward_in)),
         connect_type(connect_type_in) {}
      dict_chem_comp_tree_t() {}
   };


   // ------------------------------------------------------------------------
   // class dictionary_match_info_t
   // ------------------------------------------------------------------------
   //

   class dictionary_match_info_t; // yes.

   // ------------------------------------------------------------------------
   // class dictionary_residue_restraints_t
   // ------------------------------------------------------------------------
   //

   class dictionary_residue_restraints_t {

      class eraser {
      public:
         std::vector<std::string> baddies;
         explicit eraser(const std::vector<std::string> &baddies_in) : baddies(baddies_in) {}
         bool operator()(const dict_atom &at) const {
            return (std::find(baddies.begin(), baddies.end(), at.atom_id_4c) != baddies.end());
         }
         bool operator()(const dict_bond_restraint_t &br) {
            if (std::find(baddies.begin(), baddies.end(), br.atom_id_1_4c()) != baddies.end())
               return true;
            if (std::find(baddies.begin(), baddies.end(), br.atom_id_2_4c()) != baddies.end())
               return true;
            return false;
         }
         bool operator()(const dict_angle_restraint_t &ar) {
            if (std::find(baddies.begin(), baddies.end(), ar.atom_id_1_4c()) != baddies.end())
               return true;
            if (std::find(baddies.begin(), baddies.end(), ar.atom_id_2_4c()) != baddies.end())
               return true;
            if (std::find(baddies.begin(), baddies.end(), ar.atom_id_3_4c()) != baddies.end())
               return true;
            return false;
         }
      };

      class atom_pair_t {
      public:
         mmdb::Atom *at_1;
         mmdb::Atom *at_2;
         atom_pair_t() {
            at_1 = 0;
            at_2 = 0;
         }
         atom_pair_t(mmdb::Atom *a1, mmdb::Atom *a2) { at_1 = a1; at_2 = a2; }
         bool operator==(const std::pair<mmdb::Atom *, mmdb::Atom *> &t) {
            if (t.first == at_1) {
               if (t.second == at_2) {
                  return true;
               } else {
                  return false;
               }
            } else {
               return false;
            }
         }
         // return null if they both match.
         mmdb::Atom *shared_atom(const atom_pair_t &pair_in) {
            mmdb::Atom *shared_atom = NULL;
            if (pair_in.at_1 == at_1) {
               if (pair_in.at_2 != at_2) {
                  shared_atom = at_1;
               }
            } else {
               if (pair_in.at_2 == at_2) {
                  shared_atom = at_2;
               }
            }

            // now with swapped indices
            if (pair_in.at_1 == at_2) {
               if (pair_in.at_2 != at_1) {
                  shared_atom = at_2;
               }
            } else {
               if (pair_in.at_2 == at_1) {
                  shared_atom = at_1;
               }
            }
            return shared_atom;
         }
      };
      void init(mmdb::Residue *r);
      bool has_partial_charges_flag;
      bool nuclear_distances_flag; // oops. I didn't know that bond_length_type_t NUCLEAR_POSITION was a thing.
                                   // this needs to be consolidated.
      bool filled_with_bond_order_data_only_flag; // if set, this means that
                                        // there is only bond orders
                                        // (at the moment) and atom
                                        // names.
      std::vector<std::pair<std::string, std::string> >
      extra_name_swaps_from_name_clash(const std::vector<std::pair<std::string, std::string> > &change_name) const;
      std::string invent_new_name(const std::string &ele,
                                  const std::vector<std::string> &other_invented_names) const;


      void write_cif_pdbx_chem_comp_descriptor(mmdb::mmcif::Data *data) const;

      // imol_enc can be the model molecule number or
      // -1 for all
      // -2 for auto
      // -3 for unset
      // int imol_enc;

      void delete_atoms_from_restraints(const std::vector<std::string> &H_atoms_to_be_deleted);

   public:
      dictionary_residue_restraints_t(const std::string &comp_id_in,
                                      int read_number_in) {
         // a dictionary_residue_restraints_t no longer has a comp_id
         residue_info.comp_id = comp_id_in;
         has_partial_charges_flag = 0;
         read_number = read_number_in;
         filled_with_bond_order_data_only_flag = 0;
         nuclear_distances_flag = false;
      }
      dictionary_residue_restraints_t() {
         filled_with_bond_order_data_only_flag = 0;
         has_partial_charges_flag = 0;
         read_number = -1;
         nuclear_distances_flag = false;
      }
      explicit dictionary_residue_restraints_t(bool constructor_for_srs_restraints) {
         filled_with_bond_order_data_only_flag = 1;
         has_partial_charges_flag = 0;
         read_number = -1;
         nuclear_distances_flag = false;
         if (constructor_for_srs_restraints) {
         }
      }
      // fake a dictionary (bond and angle restraints) from the
      // coordinates in the residue.  Fake up some bond and angle
      // esds.
      explicit dictionary_residue_restraints_t(mmdb::Residue *residue_p);
      explicit dictionary_residue_restraints_t(mmdb::Manager *mol); // mol contains one residue in a hierarchy

      std::string cif_file_name;
      void clear_dictionary_residue();
      bool is_filled() const {
         if (bond_restraint.size() > 0)
            if (atom_info.size() > 0)
               return true;
         return false;
      }
      dict_chem_comp_t residue_info;
      std::vector <dict_atom> atom_info;
      unsigned int number_of_atoms() const { return atom_info.size(); }
      unsigned int number_of_non_hydrogen_atoms() const;
      // if the name matches a from (first), change it to the second
      void atom_id_swap(const std::vector<std::pair<std::string, std::string> > &from_tos);
      // std::string comp_id; // i.e. residue type name
      std::string comp_id() const { return residue_info.comp_id; }
      std::vector <dict_chem_comp_tree_t> tree;
      int read_number;
      std::vector    <dict_bond_restraint_t>    bond_restraint;
      std::vector   <dict_angle_restraint_t>   angle_restraint;
      std::vector <dict_torsion_restraint_t> torsion_restraint;
      std::vector  <dict_chiral_restraint_t>  chiral_restraint;
      std::vector   <dict_plane_restraint_t>   plane_restraint;
      std::vector   <dict_improper_dihedral_restraint_t>   improper_dihedral_restraint;
      pdbx_chem_comp_descriptor_container_t descriptors;
      pdbx_chem_comp_description_generator_t description_generation;
      chem_comp_atom_depiction_t depiction;
      gphl_chem_comp_info_t gphl_chem_comp_info;

      // Return 1 for hydrogen or deuterium, 0 for not found or not a hydrogen.
      bool is_hydrogen(const std::string &atom_name) const;
      bool is_hydrogen(unsigned int ) const; // the index of an atom in atom_info is passed.
      int assign_chiral_volume_targets(); // return the number of targets made.
      // only call the following when passing an H
      std::string get_bonded_atom(const std::string &H_atom_name) const;
      bool has_unassigned_chiral_volumes() const;
      bool has_partial_charges() const;
      void set_has_partial_charges(bool state) {
         has_partial_charges_flag = state;
      }
      std::vector<dict_torsion_restraint_t> get_non_const_torsions(bool include_hydrogen_torsions_flag) const;

      // compares atoms of torsion_restraint vs the ring atoms.
      // bool is_ring_torsion(const dict_torsion_restraint_t &torsion_restraint) const;
      bool is_ring_torsion(const atom_name_quad &quad) const;

      // I think that 20 improper dihedrals can be calculated faster than a plane restraint
      // on the graphics card (and possibly CPU):
      // So if this is a plane restraint, what are the sets of 4 atoms for the
      // improper dihedrals? Atoms to be constructed in "known" order.
      // Call this after the restraint has been filled (well, the angle restraints, at least).
      // Maybe instead of returning, I want to fill a
      // std::vector<improper_dihedral_restraints_t> improper_restraint
      // It will look like a chiral restraint.
      std::vector<atom_name_quad> plane_restraint_to_improper_dihedrals(unsigned int i) const;

      void set_use_nuclear_distances(bool state) { nuclear_distances_flag = state; }
      void write_cif(const std::string &filename) const;
      // look up the atom id in the atom_info (dict_atom vector)
      std::string atom_name_for_tree_4c(const std::string &atom_id) const;
      // quote atom name as needed - i.e. CA -> CA, CA' -> "CA'"
      std::string quoted_atom_name(const std::string &an) const;

      // look up the atom id in the atom_info (dict_atom vector).
      // Return "" on no atom found with name atom_name;
      //
      std::string element(const std::string &atom_name) const;

      bool is_bond_to_hydrogen_atom(const dict_bond_restraint_t &br) const;

      // likewise look up the energy type.  Return "" on no atom found
      // with that atom_name.
      //
      std::string type_energy(const std::string &atom_name) const;

      std::vector<std::vector<std::string> > get_ligand_ring_list() const;

      // return null on failure
      mmdb::Residue *GetResidue(bool idealize_flag, float b_factor) const;

      // This is very slow if you call it a number of times.
      // Better to extract the ring info with get_ligand_ring_list()
      // and test for atom_name_1 and atom_name_2 using that.
      bool in_same_ring(const std::string &atom_name_1, const std::string &atom_name_2) const;

      // Here for convenience, but it doesn't rely on class functions or data items
      // (could/should be static?)
      bool in_same_ring(const std::string &atom_name_1, const std::string &atom_name_2,
                        const std::vector<std::vector<std::string> > &ring_list) const;

      void add_pyranose_pseudo_ring_plane_restraints(const std::string &plane_id,
                                                     std::vector<std::string> &atom_name_vec,
                                                     double esd);

      bool ligand_has_aromatic_bonds_p() const;

      std::vector<std::vector<std::string> > get_ligand_aromatic_ring_list() const;

      std::vector<std::string> get_attached_H_names(const std::string &atom_name) const;

      bool is_bond_order_data_only() const { return filled_with_bond_order_data_only_flag; }

      std::vector<std::string> neighbours(const std::string &atom_name, bool allow_hydrogen_neighbours_flag) const;
      // same thing with indexing into the atom_info vector.  No protection
      // for out of bounds atom_idx value (i.e. atom_idx must be valid).
      std::vector<unsigned int> neighbours(unsigned atom_idx,
                                           bool allow_hydrogen_neighbours_flag) const;

      // return "" on not found
      std::string get_bond_type(const std::string &name_1, const std::string &name_2) const;

      // replace the restraints that we have with new_restraints,
      // keeping restraints that in the current set but not in
      // new_restraints
      void conservatively_replace_with(const dictionary_residue_restraints_t &new_restraints);
      void conservatively_replace_with_bonds (const dictionary_residue_restraints_t &new_restraints);
      void conservatively_replace_with_angles(const dictionary_residue_restraints_t &new_restraints);
      void replace_coordinates(const dictionary_residue_restraints_t &mon_res_in);

      //
      void remove_redundant_plane_restraints();
      bool is_redundant_plane_restraint(std::vector<dict_plane_restraint_t>::iterator it) const;
      bool is_redundant_plane_restraint(unsigned int idx) const;
      void reweight_subplanes(); // if an atom is in more than one plane restraint, then
                                 // reduce its esd.

      // quiet means don't tell me about matches
      bool compare(const dictionary_residue_restraints_t &new_restraints,
                   double bond_length_tolerance,
                   double bond_esd_tolerance,
                   double angle_tolerance,
                   double angle_esd_tolerance,
                   bool compare_hydrogens=false,
                   bool output_energy_types=false,
                   bool quiet=false) const;

      // return a dictionary that is a copy of this dictionary, but
      // trying to match the names of the atoms of ref.  Do graph
      // matching to find the set of atom names that match/need to be
      // changed.
      // 

      // If new_comp_id is "auto", suggest_new_comp_id() is called to
      // generate a comp_id string.
      //
      // If residue_p is not null, then change the atom names
      // in the residue.
      //
      dictionary_match_info_t
      match_to_reference(const dictionary_residue_restraints_t &ref,
                         mmdb::Residue *residue_p,
                         const std::string &new_comp_id,
                         const std::string &new_compound_name) const;

      // change the atom names and the residue type of the passed residue.
      bool change_names(mmdb::Residue *residue_p,
                        const std::vector<std::pair<std::string, std::string> > &change_name,
                        const std::string &new_comp_id) const;

      // make a mmdb::math::Graph from the atom_info and bond restraints.
      //
      // Caller disposes of the memory with a delete().
      mmdb::math::Graph *make_graph(bool use_hydrogen) const;

      // Are the atoms only of elements C,N.O,H,F,Cl,I,Br,P,S?
      bool comprised_of_organic_set() const;

      // are the number of atoms of each element the same ie. they have the same chemical formula?
      //
      bool composition_matches(const dictionary_residue_restraints_t &other) const;

      // for hydrogens
      bool is_connected_to_donor(const std::string &H_at_name_4c,
                                 const energy_lib_t &energy_lib) const;

      // return an empty string on failure
      std::string get_other_H_name(const std::string &H_at_name) const;
      // return an empty vector on failure
      std::vector<std::string> get_other_H_names(const std::string &H_at_name) const;

      void remove_phosphate_hydrogens();
      void remove_sulphate_hydrogens();
      void remove_PO4_SO4_hydrogens(const std::string &P_or_S);
      void remove_carboxylate_hydrogens();
      void move_3GP_atoms();

      friend std::ostream& operator<<(std::ostream &s, const dictionary_residue_restraints_t &rest);

#ifdef HAVE_CCP4SRS
      bool fill_using_ccp4srs(ccp4srs::Manager *srs_manager, const std::string &monomer_type);
#endif // HAVE_CCP4SRS
      
   };
   std::ostream& operator<<(std::ostream &s, const dictionary_residue_restraints_t &rest);

   class dictionary_match_info_t {
   public:
      unsigned int n_matches;
      dictionary_residue_restraints_t dict;
      std::vector<std::pair<std::string, std::string> > name_swaps;
      std::vector<std::string> same_names;
      std::string new_comp_id;
      dictionary_match_info_t() {
         n_matches = 0;
      }
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
                                  mmdb::realtype value_dist_in,
                                  mmdb::realtype value_dist_esd_in) :
         basic_dict_restraint_t(atom_id_1, atom_id_2) {
         atom_1_comp_id = atom_1_comp_id_in;
         atom_2_comp_id = atom_2_comp_id_in;
         value_dist = value_dist_in;
         value_dist_esd = value_dist_esd_in;
      }
      int atom_1_comp_id, atom_2_comp_id;
      double dist() const { return value_dist; }
      double esd()  const { return value_dist_esd; }
      bool matches(const std::string &at_name_1, const std::string &at_name_2) const {
         bool status = true;
         if (at_name_1 != atom_id_1()) status = false;
         if (at_name_2 != atom_id_2()) status = false;
         return status;
      }
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
                                   mmdb::realtype value_angle_in,
                                   mmdb::realtype value_angle_esd_in) :
         basic_dict_restraint_t(atom_id_1, atom_id_2), atom_id_3_(atom_id_3_in) {
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
                                     const std::string &id_in) :
         basic_dict_restraint_t(atom_id_1_in, atom_id_2_in),
         atom_id_3_(atom_id_3_in),
         atom_id_4_(atom_id_4_in), id_(id_in) {

         atom_1_comp_id = atom_1_comp_id_in;
         atom_2_comp_id = atom_2_comp_id_in;
         atom_3_comp_id = atom_3_comp_id_in;
         atom_4_comp_id = atom_4_comp_id_in;
         angle_ = value_angle_in;
         angle_esd_ = value_angle_esd;
         period_ = period;
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
                                   ) : basic_dict_restraint_t(atom_id_1_in, atom_id_2_in),
                                       atom_1_comp_id(atom_1_comp_id_in),
                                       atom_2_comp_id(atom_2_comp_id_in),
                                       atom_3_comp_id(atom_3_comp_id_in),
                                       atom_c_comp_id(atom_c_comp_id_in),
                                       atom_id_c_(atom_id_c_in),
                                       atom_id_3_(atom_id_3_in),
                                       chiral_id(chiral_id_in) {
         volume_sign = volume_sign_in;
         target_volume_ = -999.9;  // unassigned
         target_sigma_  = -999.9;
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
                                  double dist_esd) : plane_id(plane_id_in) {
         dist_esd_ = dist_esd;
         atom_ids.push_back(atom_id);
         atom_comp_ids.push_back(atom_comp_id);
      }
      std::string plane_id; 
      std::vector<std::string> atom_ids;
      std::vector<int> atom_comp_ids;
      unsigned int n_atoms() const { return atom_ids.size(); }
      double dist_esd() const { return dist_esd_; }
      void set_dist_esd(double d) { dist_esd_ = d; }
      std::string atom_id(int i) const { return atom_id_mmdb_expand(atom_ids[i]); }
   };

   // for link_ids such as TRANS, PTRANS (proline), CIS etc.
   // 
   class dictionary_residue_link_restraints_t {
   public:
      explicit dictionary_residue_link_restraints_t(const std::string &link_id_in) : link_id(link_id_in) {}
      dictionary_residue_link_restraints_t() : link_id("") {}
      std::string link_id;
      std::vector <dict_link_bond_restraint_t>    link_bond_restraint;
      std::vector <dict_link_angle_restraint_t>   link_angle_restraint;
      std::vector <dict_link_torsion_restraint_t> link_torsion_restraint;
      std::vector <dict_link_plane_restraint_t>   link_plane_restraint;
      std::vector <dict_link_chiral_restraint_t>  link_chiral_restraint;
      int assign_link_chiral_volume_targets(); // return the number of link targets made.
      bool has_unassigned_chiral_volumes() const;
      bool empty() const { return link_id.empty(); }
   };

   class simple_cif_reader {
      std::vector<std::string> names;
      std::vector<std::string> three_letter_codes;
   public:
      explicit simple_cif_reader(const std::string &cif_dictionary_file_name);
      bool has_restraints_for(const std::string &res_type);
   };

   class chem_link {
   public:
      std::string id;
      std::string chem_link_comp_id_1;
      std::string chem_link_mod_id_1;
      std::string chem_link_group_comp_1;
      std::string chem_link_comp_id_2;
      std::string chem_link_mod_id_2;
      std::string chem_link_group_comp_2;
      std::string chem_link_name;
      unsigned int hash_code;
      chem_link() { hash_code = 0; }
      chem_link(const std::string &id_in,
                const std::string &chem_link_comp_id_1_in,
                const std::string &chem_link_mod_id_1_in,
                const std::string &chem_link_group_comp_1_in,
                const std::string &chem_link_comp_id_2_in,
                const std::string &chem_link_mod_id_2_in,
                const std::string &chem_link_group_comp_2_in,
                const std::string &chem_link_name_in) :
         id(id_in),
         chem_link_comp_id_1(chem_link_comp_id_1_in),
         chem_link_mod_id_1(chem_link_mod_id_1_in),
         chem_link_group_comp_1(chem_link_group_comp_1_in),
         chem_link_comp_id_2(chem_link_comp_id_2_in),
         chem_link_mod_id_2(chem_link_mod_id_2_in),
         chem_link_group_comp_2(chem_link_group_comp_2_in),
         chem_link_name(chem_link_name_in)
      {
         hash_code = make_hash_code(chem_link_comp_id_1_in, chem_link_comp_id_2,
                                    chem_link_group_comp_1_in, chem_link_group_comp_2_in);
      }
      friend std::ostream& operator<<(std::ostream &s, chem_link lnk);
      static unsigned int make_hash_code(const std::string &comp_id_1, const std::string &comp_id_2, const std::string &group_1, const std::string &group_2);
      unsigned int get_hash_code() const { return hash_code; }
      bool operator==(const chem_link &cl) const {
         // cheap test - relies on the dictionary having unique entries
         return (cl.Id() == id);
      }
      bool operator<(const chem_link &cl) const {
         if (cl.Id() < id) {
            return true;
         } else {
            if (cl.Id() > id) {
               return false;
            } else {
              // Ids equal
              return (cl.get_hash_code() < hash_code);
            }
         }
      }
      // pair: matches need-order-switch-flag
      std::pair<bool, bool> matches_comp_ids_and_groups_hashed(const std::string &comp_id_1,
                                                        const std::string &group_1,
                                                        const std::string &comp_id_2,
                                                        const std::string &group_2) const;
      // caller should handle the reversing or the order - because they may match be *both* ways round
      bool matches_comp_ids_and_groups(const std::string &comp_id_1,
                                       const std::string &group_1,
                                       const std::string &comp_id_2,
                                       const std::string &group_2) const;
      std::pair<bool, bool> matches_comp_ids_and_groups_hashed(unsigned int hash_code_forward,
                                                               unsigned int hash_code_backwards) const;
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


   // ---------------------------------------------------------------
   // helper classes for linkage selection
   // ---------------------------------------------------------------

   class glycosidic_distance {
   public:
      double distance;
      mmdb::Atom *at1;
      mmdb::Atom *at2;
      glycosidic_distance(mmdb::Atom *at1_in, mmdb::Atom *at2_in, double d) {
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

   class read_refmac_mon_lib_info_t {
   public:
      unsigned int n_atoms;
      unsigned int n_bonds;
      unsigned int n_links;
      int monomer_idx;
      std::vector<std::string> error_messages;
      bool success;
      std::string comp_id; // the first if there are many.  Test for blank when used.
      read_refmac_mon_lib_info_t() {
         n_atoms = 0;
         n_bonds = 0;
         n_links = 0;
         success = true;
         monomer_idx = -1;
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

      enum { MON_LIB_LIST_CIF = -999}; // A negative number special
                                       // flag to tell the reader that
                                       // this is not a normal
                                       // residue's restraints and
                                       // that the chem_comp should
                                       // not be added to the list
                                       // (currently).

      enum { UNSET_NUMBER = -1 };  // An unset number, for example the
      // number of atoms.

      class restraint_eraser {
      public:
         std::vector<std::string> names;
         // the constructor, can be information that needs to be used
         // internally in the operator() function.  This is run once
         // 
         explicit restraint_eraser(const std::vector<std::string> &names_in) : names(names_in) {}

         explicit restraint_eraser(const std::set<std::string> &names_in) {
            std::set<std::string>::const_iterator it;
            for (it=names_in.begin(); it != names_in.end(); ++it)
               names.push_back(*it);
         }

         // return true for deletion
         bool operator()(const dict_torsion_restraint_t &r) const {
            int n_match = 0;
            for (unsigned int i=0; i<names.size(); i++) {
               if (r.atom_id_1_4c() == names[i]) n_match++;
               if (r.atom_id_2_4c() == names[i]) n_match++;
               if (r.atom_id_3_4c() == names[i]) n_match++;
               if (r.atom_id_4_4c() == names[i]) n_match++;
            }
            return (n_match == 4);
         }
      };

      //testing func
      bool close_float_p(const mmdb::realtype &f1, const mmdb::realtype &f2) const {
         float d = fabs(f1-f2);
         if (d < 0.001)
            return true;
         else
            return false;
      }
      
      // std::vector<simple_residue_t> residue; 
      std::vector<std::string> residue_codes;
      bool verbose_mode;

      // int imol_enc_current; // given the current read of a dictionary
      // set the molecule_number for those entries to imol_enc
      // the first number of the pair is the imol_enc
      //
      std::vector<std::pair<int, dictionary_residue_restraints_t> > dict_res_restraints;
      std::vector<dictionary_residue_link_restraints_t> dict_link_res_restraints;
      std::map<unsigned int, std::vector<chem_link> > chem_link_map;
      std::vector<list_chem_mod>  chem_mod_vec;

      // the monomer data in list/mon_lib_list.cif, not the
      // restraints, just id, 3-letter-code, name, group,
      // number-of-atoms, description_level.
      // Added to by the simple_mon_lib* functions.
      //
      // Use a map for faster lookups.  the key is the comp_id;
      // 
      std::map<std::string,dictionary_residue_restraints_t> simple_monomer_descriptions;

      int  comp_atom(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc, bool is_from_pdbx_model_atom=false); 
      std::string comp_atom_pad_atom_name(const std::string &atom_id, const std::string &type_symbol) const;
      // return the comp_id
      std::string chem_comp(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);
      void comp_tree   (mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);
      int  comp_bond   (mmdb::mmcif::PLoop mmCIFLoop, int imol_enc, bool is_from_pdbx_model_bond=false);
      void comp_angle  (mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);
      void comp_torsion(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);
      void comp_plane  (mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);
      std::pair<int, std::vector<std::string> >
      comp_chiral(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);
                                                 // return the number of chirals and a vector
                                                 // of monomer names that have had
                                                 // chirals added (almost certainly just
                                                 // one of them, of course).

      void chem_comp_acedrg(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);

      void add_chem_links (mmdb::mmcif::PLoop mmCIFLoop); // references to the modifications
                                                // to the link groups (the modifications
                                                // themselves are in data_mod_list)
      int  link_bond   (mmdb::mmcif::PLoop mmCIFLoop); 
      void link_angle  (mmdb::mmcif::PLoop mmCIFLoop); 
      void link_torsion(mmdb::mmcif::PLoop mmCIFLoop); 
      void link_plane  (mmdb::mmcif::PLoop mmCIFLoop);
      int  link_chiral  (mmdb::mmcif::PLoop mmCIFLoop); // return number of new chirals
      void pdbx_chem_comp_descriptor(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);

      void pdbe_chem_comp_atom_depiction(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);

      void pdbx_chem_comp_description_generator(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);

      void gphl_chem_comp_info(mmdb::mmcif::PStruct structure, int imol_enc);

      // return the comp id (so that later we can associate the file name with the comp_id).
      // 
      std::string chem_comp_component( mmdb::mmcif::PStruct structure, int imol_enc);
      std::string pdbx_chem_comp_model(mmdb::mmcif::PStruct structure, int imol_enc);
      // non-looping (single) tor
      void chem_comp_tor_structure(mmdb::mmcif::PStruct structure, int imol_enc);
      // non-looping (single) chir
      void chem_comp_chir_structure(mmdb::mmcif::PStruct structure, int imol_enc);

      void parse_lib_info(mmdb::mmcif::PStruct structure);

      void mon_lib_add_chem_comp(const std::string &comp_id,
                                 int imol_enc,
                                 const std::string &three_letter_code,
                                 const std::string &name,
                                 const std::string &group,
                                 int number_atoms_all, int number_atoms_nh,
                                 const std::string &description_level);

      // old
      void mon_lib_add_atom(const std::string &comp_id,
                            int imol_enc,
                            const std::string &atom_id,
                            const std::string &atom_id_4c,
                            const std::string &type_symbol,
                            const std::string &type_energy,
                            const std::pair<bool, mmdb::realtype> &partial_charge,
                            const std::pair<bool, int> &formal_charge,
                            dict_atom::aromaticity_t arom_in,
                            const std::pair<bool, clipper::Coord_orth> &model_pos,
                            const std::pair<bool, clipper::Coord_orth> &model_pos_ideal);

      void mon_lib_add_atom(const std::string &comp_id,
                            int imol_enc,
                            const dict_atom &atom_info);

      void mon_lib_add_acedrg_atom_type(const std::string &comp_id, int imol_enc,
                                        const std::string &atom_id,
                                        const std::string &atom_type);

      // called because they were all at origin, for example.
      void delete_atom_positions(const std::string &comp_id, int imol_enc, int pos_type);
                            
      void mon_lib_add_tree(std::string comp_id,
                            int imol_enc,
                            std::string atom_id,
                            std::string atom_back,
                            std::string atom_forward,
                            std::string connect_type);

      // If value_dist_nuclear_esd is 0 or less than 0, the don't use
      // the value_dist_nuclear or value_dist_nuclear_esd
      void mon_lib_add_bond(std::string comp_id,
                            int imol_enc,
                            std::string atom_id_1, std::string atom_id_2,
                            std::string type,
                            mmdb::realtype value_dist, mmdb::realtype value_dist_esd,
                            mmdb::realtype value_dist_nuclear, mmdb::realtype value_dist_nuclear_esd,
                            dict_bond_restraint_t::aromaticity_t arom_in,
                            dict_bond_restraint_t::bond_length_type_t lbt_in);

      void mon_lib_add_bond_no_target_geom(std::string comp_id,
                                           int imol_enc,
                                           std::string atom_id_1, std::string atom_id_2,
                                           std::string type,
                                           dict_bond_restraint_t::aromaticity_t arom_in);
      
      void mon_lib_add_angle(std::string comp_id,
                             int imol_enc,
                             std::string atom_id_1,
                             std::string atom_id_2,
                             std::string atom_id_3,
                             mmdb::realtype value_angle, mmdb::realtype value_angle_esd);

      void mon_lib_add_torsion(std::string comp_id,
                               int imol_enc,
                               std::string torsion_id,
                               std::string atom_id_1,
                               std::string atom_id_2,
                               std::string atom_id_3,
                               std::string atom_id_4,
                               mmdb::realtype value_angle, mmdb::realtype value_angle_esd,
                               int period);

      void mon_lib_add_chiral(std::string comp_id,
                              int imol_enc,
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
                             int imol_enc,
                             const std::string &plane_id,
                             const std::string &atom_id,
                             const mmdb::realtype &dist_esd);

      void add_restraint(std::string comp_id, int imol_enc, const dict_bond_restraint_t &restr);
      void add_restraint(std::string comp_id, int imol_enc, const dict_angle_restraint_t &restr);
      void add_restraint(std::string comp_id, int imol_enc, const dict_torsion_restraint_t &restr);
      void add_restraint(std::string comp_id, int imol_enc, const dict_chiral_restraint_t &rest);

      void add_pdbx_descriptor(const std::string &comp_id,
                               int imol_enc,
                               pdbx_chem_comp_descriptor_item &descr);

      // comp_tree need to convert unquoted atom_ids to 4char atom
      // names.  So we look them up in the atom table.  If name not
      // found, return the input string
      //
      std::string atom_name_for_tree_4c(const std::string &comp_id, const std::string &atom_id) const;

      // for simple monomer descriptions:

      // return the comp id (so that later we can associate the file name with the comp_id).
      // 
      std::string simple_mon_lib_chem_comp(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc);
      // add to simple_monomer_descriptions not dict_res_restraints.

      void add_cif_file_name(const std::string &cif_filename,
                             const std::string &comp_id1,
                             const std::string &comp_id2,
                             int imol_enc);

      void simple_mon_lib_add_chem_comp(const std::string &comp_id,
                                        int imol_enc,
                                        const std::string &three_letter_code,
                                        const std::string &name,
                                        const std::string &group,
                                        int number_atoms_all, int number_atoms_nh,
                                        const std::string &description_level);

      // mod stuff (references by chem links)

      int add_chem_mods(mmdb::mmcif::PData data);
      int add_chem_mod(mmdb::mmcif::PLoop mmCIFLoop);
      int add_mod(mmdb::mmcif::PData data);

      int add_mods(mmdb::mmcif::PData data);
      int add_chem_mods(mmdb::mmcif::PLoop mmCIFLoop);
      
      // which calls:
      void add_chem_mod_atom( mmdb::mmcif::PLoop mmCIFLoop);
      void add_chem_mod_bond( mmdb::mmcif::PLoop mmCIFLoop);
      void add_chem_mod_tree( mmdb::mmcif::PLoop mmCIFLoop);
      void add_chem_mod_angle(mmdb::mmcif::PLoop mmCIFLoop);
      void add_chem_mod_tor(  mmdb::mmcif::PLoop mmCIFLoop);
      void add_chem_mod_chir( mmdb::mmcif::PLoop mmCIFLoop);
      void add_chem_mod_plane(mmdb::mmcif::PLoop mmCIFLoop);


      // synonyms (for RNA/DNA)
      // 
      void add_synonyms(mmdb::mmcif::PData data);
      void add_chem_comp_synonym(mmdb::mmcif::PLoop mmCIFLoop);
      class residue_name_synonym {
      public:
         residue_name_synonym(std::string &comp_id_in,
                              std::string &comp_alternative_id_in,
                              std::string &mod_id_in) :
            comp_id(comp_id_in),
            comp_alternative_id(comp_alternative_id_in),
            mod_id(mod_id_in) {}
         std::string comp_id;
         std::string comp_alternative_id;
         std::string mod_id;
      };
      std::vector<residue_name_synonym> residue_name_synonyms;

      // link stuff
      int init_links(mmdb::mmcif::PData data);

      void link_add_bond(const std::string &link_id,
                         int atom_1_comp_id,
                         int atom_2_comp_id,
                         const std::string &atom_id_1,
                         const std::string &atom_id_2,
                         mmdb::realtype value_dist,
                         mmdb::realtype value_dist_esd);
      
      void link_add_angle(const std::string &link_id,
                          int atom_1_comp_id,
                          int atom_2_comp_id,
                          int atom_3_comp_id,
                          const std::string &atom_id_1,
                          const std::string &atom_id_2,
                          const std::string &atom_id_3,
                          mmdb::realtype value_dist,
                          mmdb::realtype value_dist_esd);

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
                            mmdb::realtype value_dist,
                            mmdb::realtype value_dist_esd,
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
      void filter_chiral_centres(int imol, const std::vector<std::string> & comp_id_for_filtering);

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
      // 20240110-PE don't use this function - use get_monomer_restraints_index().
      int get_monomer_type_index(const std::string &monomer_type);
      std::string get_padded_name(const std::string &atom_id, const int &comp_id_index) const;

      // return a list of torsion restraints that are unique for atoms 2 and
      // 3 (input restraints vector can potentially have many restraints
      // that have the same atoms 2 and 3).
      //
      std::vector <dict_torsion_restraint_t>
      filter_torsion_restraints(const std::vector <dict_torsion_restraint_t> &restraints_in) const;
      static bool torsion_restraints_comparer(const dict_torsion_restraint_t &a, const dict_torsion_restraint_t &b);

      std::vector<chem_link> matching_chem_links(const std::string &comp_id_1,
                                                 const std::string &group_1,
                                                 const std::string &comp_id_2,
                                                 const std::string &group_2,
                                                 bool allow_peptide_link_flag) const;

      energy_lib_t energy_lib;

      std::pair<bool, dictionary_residue_restraints_t>
      get_monomer_restraints_internal(const std::string &monomer_type,
                                      int imol_enc,
                                      bool allow_minimal_flag) const;

      std::vector<std::string> non_auto_load_residue_names;
      void fill_default_non_auto_load_residue_names(); // build-it defaults

      // return empty file name on failure.
      std::string comp_id_to_file_name(const std::string &comp_id) const;

#ifdef HAVE_CCP4SRS
      ccp4srs::Manager *ccp4srs;
#endif

      void add_molecule_number_to_entries(const std::vector<std::string> &comp_ids, int imol_enc);

      std::vector<atom_name_torsion_quad>
      get_reference_monomodal_torsion_quads(const std::string &res_name) const;

   public:

      protein_geometry() {
         read_number = 0;
         set_verbose(1);
         parse_metal_NOS_distance_tables();
#if HAVE_CCP4SRS
         ccp4srs = NULL;
#endif
         fill_default_non_auto_load_residue_names();
      }

      enum { IMOL_ENC_ANY = -999999, IMOL_ENC_AUTO = -999998, IMOL_ENC_UNSET = -999997 };

      // CCP4 SRS things

      int init_ccp4srs(const std::string &ccp4srs_dir); // inits CCP4SRS
      void read_ccp4srs_residues();
      // return NULL on unable to make residue
      mmdb::Residue *get_ccp4srs_residue(const std::string &res_name) const;

      // return a vector of compound ids and compound names for
      // compounds that contain compound_name_substring.
      //
      std::vector<std::pair<std::string, std::string> >
      matching_ccp4srs_residues_names(const std::string &compound_name_substring) const;

      // and fills these chem mod classes, simple container class
      // indexed with map on the mod_id

      std::map<std::string, chem_mod> mods;
      // can throw an std::runtime exception if there is no chem mod
      // for the link_id (and that's fine)
      std::pair<chem_mod, chem_mod> get_chem_mods_for_link(const std::string &link_id) const;
      void debug_mods() const;

      // when making the bonds mesh, we would like to have a mesh that includes a hemisphere
      // where there are atoms that have only one non-Hydrogen bond. So the dictionary
      // needs to know where such atoms here. We set them in this function.
      void set_only_bonds(int dict_idx);

      // Refmac monomer lib things
      //
      // Return the number of bond restraints
      //
      read_refmac_mon_lib_info_t
      init_refmac_mon_lib(std::string filename, int read_number_in,
                          int imol_enc = IMOL_ENC_ANY);

      unsigned int size() const { return dict_res_restraints.size(); }
      const std::pair<int, dictionary_residue_restraints_t> & operator[](int i) const {
         // debugging SRS compilation
         // std::cout << "const operator[] for a geom " << i << " of size "
         //           << dict_res_restraints.size() << std::endl;
         return dict_res_restraints[i]; }
      const dictionary_residue_link_restraints_t & link(int i) const {
         return dict_link_res_restraints[i]; }

      dictionary_residue_link_restraints_t link(const std::string &id_in) const;

      // return "" on comp_id not found, else return the file name.
      //
      std::string get_cif_file_name(const std::string &comp_id, int imol_enc) const;

      int link_size() const { return dict_link_res_restraints.size(); }
      void info() const;
      std::string three_letter_code(const unsigned int &i) const;

      void set_verbose(bool verbose_mode_in);

      int init_standard(); // standard protein residues and links.
                                  // Return the current read_number

      // if the dictionary is already in the store, then do nothing, otherwise
      // try_dynamic_add().
      //
      // Internally this function wraps have_dictionary_for_residue_type() which
      // is not as good a name as check_and_try_dynamic_add() for what the function
      // does.
      //
      int check_and_try_dynamic_add(const std::string &resname, int imol_enc, int read_number);  // return success status?

      // Return 0 on failure to do a dynamic add, otherwise return the
      // number of atoms read.
      //
      int try_dynamic_add(const std::string &resname, int read_number);  // return success status?
      // this is not const if we use dynamic add.

      bool matches_imol(int imol_dict, int imol_enc) const;

      // return true on having deleted;
      bool delete_mon_lib(const std::string &comp_id, int imol_enc); // delete comp_id from dict_res_restraints
                                                                        // (if it exists there).
                                                                     // 20161004 we also need to match
                                                                     // imols before deletion occurs

      void print_dictionary_store() const;

      // return a pair, the first is status (1 if the name was found, 0 if not)
      //
      std::pair<bool, std::string> get_monomer_name(const std::string &comp_id, int imol_enc) const;

      // return 2-3 filtered torsions
      std::vector <dict_torsion_restraint_t>
      get_monomer_torsions_from_geometry(const std::string &monomer_type,
                                         int imol_enc);
      std::vector <dict_chiral_restraint_t>
      get_monomer_chiral_volumes(const std::string &monomer_type,
                                 int imol_enc) const;

      // as above, except filter out of the returned vectors torsions
      // that move (or are based on) hydrogens.
      // return 2-3 filtered torsions
      //
      std::vector <dict_torsion_restraint_t>
      get_monomer_torsions_from_geometry(const std::string &monomer_type,
                                         int imol_enc,
                                         bool find_hydrogen_torsions) const;

      std::pair<bool, dict_atom> get_monomer_atom_info(const std::string &monomer_name,
                                                       const std::string &atom_name,
                                                       int imol_enc) const;

      std::vector<std::pair<int, std::string> > get_monomer_names() const;

      bool copy_monomer_restraints(const std::string &monomer_type, int imol_enc_current, int imol_enc_new);

      // Return success status in first (0 is fail) and the second is
      // a whole residue's restraints so that we can use it to test if
      // an atom is a hydrogen.  Must have a full entry (not minimal
      // for the first of the returned pair to be true.
      // 
      std::pair<bool, dictionary_residue_restraints_t>
      get_monomer_restraints(const std::string &monomer_type,
                             int imol_enc) const;

      // Return success status in first (0 is fail) and the second is
      // a whole residue's restraints so that we can use it to test if
      // an atom is a hydrogen.
      //
      // the dictionary_residue_restraints_t is returned even we have
      // a minimal restraints description (e.g. from ccp4srs).
      // 
      std::pair<bool, dictionary_residue_restraints_t>
      get_monomer_restraints_at_least_minimal(const std::string &monomer_type,
                                              int imol_enc) const;

      // caller ensures that idx is valid
      const dictionary_residue_restraints_t &
      get_monomer_restraints(unsigned int idx) const {
         return dict_res_restraints[idx].second;
      }

      // created for looking up energy types cheaply (hopefully).
      // 20161007 Now used by refinement monomer restraints lookup.
      // return -1 on restraints not found.
      // 
      int get_monomer_restraints_index(const std::string &monomer_type,
                                       int imol_enc,
                                       bool allow_minimal_flag) const;

      // 20250124-PE and the reverse
      //
      // return the second blank on lookup failure
      std::pair<int, std::string> get_monomer_name(int monomer_index) const;

      // non-const because we try to read in stuff from ccp4srs when
      // it's not in the dictionary yet.  ccp4srs gives us bond orders.
      //
      // This relies on ccp4srs being setup before we get to make this
      // call (init_ccp4srs()).
      //
      std::pair<bool, dictionary_residue_restraints_t>
      get_bond_orders(const std::string &monomer_type);

      // used to created data from ccp4srs to put into protein_geometry
      // object.
      //
      // return a flag to signify success.
      //
      bool fill_using_ccp4srs(const std::string &monomer_type);

      // If monomer_type is not in dict_res_restraints, then add a new
      // item to the dict_res_restraints and add mon_res_in.  Return 1
      // for replaced 0 for added.
      //
      bool replace_monomer_restraints(std::string monomer_type,
                                      int imol_enc,
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
                                            int imol_enc,
                                            int read_number,
                                            bool try_autoload_if_needed=true);

      // 20240110-PE I need to know if the residues have (bond) restraints - i.e. they
      // can be refined. The above function is for drawing things.
      // return false if there are no bond restraints
      bool have_restraints_dictionary_for_residue_type(const std::string &monomer_type,
                                            int imol_enc,
                                            int read_number,
                                            bool try_autoload_if_needed=true);

      // this is const because there is no dynamic add.
      //
      // if there is just an ccp4srs entry, then this returns false.
      //
      bool have_dictionary_for_residue_type_no_dynamic_add(const std::string &monomer_type, int imol) const;

      // this is const because there is no dynamic add.
      //
      // if there is (even) a ccp4srs entry, then this returns true.
      //
      bool have_at_least_minimal_dictionary_for_residue_type(const std::string &monomer_type,
                                                             int imol) const;

      // likewise not const
      bool have_dictionary_for_residue_types(const std::vector<std::string> &residue_types,
                                             int imol_enc,
                                             int read_number);

      // likewise not const.
      // Return false if there are no bond restraints
      bool have_restraints_dictionary_for_residue_types(const std::vector<std::string> &residue_types,
                                                        int imol_enc,
                                                        int read_number);

      // return a pair, overall status, and pair of residue names and
      // atom names that dont't match.
      //
      std::pair<bool, std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > >
      atoms_match_dictionary(int imol, const std::vector<mmdb::Residue *> &residues,
                             bool check_hydrogens_too_flag,
                             bool apply_bond_distance_check) const;

      std::pair<bool, std::vector<std::string> >
      atoms_match_dictionary(int imol,
                             mmdb::Residue *res,
                             bool check_hydrogens_too_flag,
                             bool apply_bond_distance_check) const;

      // return a pair: a status, yes/no atoms match and a vector of
      // atoms whose names do not match.
      //
      std::pair<bool, std::vector<std::string> >
      atoms_match_dictionary(mmdb::Residue *res,
                             bool check_hydrogens_too_flag,
                             bool apply_bond_distance_check,
                             const dictionary_residue_restraints_t &restraints) const;

      bool
      atoms_match_dictionary_bond_distance_check(mmdb::Residue *residue_p,
                                                 bool check_hydrogens_too_flag,
                                                 const dictionary_residue_restraints_t &restraints) const;

      // add "synthetic" 5 atom planar peptide restraint
      void add_planar_peptide_restraint();
      bool make_tight_planar_peptide_restraint();
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
      mmdb::Manager *mol_from_dictionary(const std::string &three_letter_code,
                                         int imol_enc,
                                         bool idealised_flag);
      
      // make HETATMs if non-standard residue name.
      mmdb::Manager *mol_from_dictionary(int monomer_index,
                                         int imol_enc,
                                         bool idealised_flag);

      // find the missing names (if any)
      // (and then call try_dynamic_add would be the typical usage.)
      std::vector<std::string> residue_names_with_no_dictionary(mmdb::Manager *mol, int imol) const;

      bool read_extra_dictionaries_for_molecule(mmdb::Manager *mol, int imol, int *read_number_p);

      // Used by above (or maybe you just want a residue?)
      // (Can return NULL).
      // 
      // Something strange happens with internal-to-a-mmdb::Residue atom
      // indexing when I tried to use this.  The problem was resoloved
      // by using mol_from_dictionary() above and getting the first
      // residue.  Something to do with atom indexing on checking
      // in...?
      // 
      mmdb::Residue *get_residue(const std::string &comp_id,
                                 int imol_enc,
                                 bool idealised_flag,
                                 bool try_autoload_if_needed=true, float b_factor=20.0);

      // Thow a std::runtime_error exception if we can't get the group of r
      std::string get_group(mmdb::Residue *r) const;
      // and the string version of this
      std::string get_group(const std::string &res_name) const;

      // bool is the need-order-switch-flag
      std::vector<chem_link>
      matching_chem_links(const std::string &comp_id_1,
                          const std::string &group_1,
                          const std::string &comp_id_2,
                          const std::string &group_2) const;

      // Try to find a link that is not a peptide link (because that
      // fails on a distance check).  This is the method to find
      // isopeptide links (which again need to be distance checked in
      // find_link_type_rigourous()).
      //
      // bool the need-order-switch-flag
      std::vector<chem_link>
      matching_chem_links_non_peptide(const std::string &comp_id_1,
                                      const std::string &group_1,
                                      const std::string &comp_id_2,
                                      const std::string &group_2) const;
      // In this version the mol is passed so that we can find links
      // that match header LINKs or SSBonds
      std::vector<chem_link>
      matching_chem_links_non_peptide(const std::string &comp_id_1,
                                      const std::string &group_1,
                                      const std::string &comp_id_2,
                                      const std::string &group_2,
                                      mmdb::Manager *mol) const;

      // return "" on failure.
      // no order switch is considered.
      //
      std::string find_glycosidic_linkage_type_by_distance(mmdb::Residue *first, mmdb::Residue *second) const;
      std::string find_glycosidic_linkage_type(mmdb::Residue *first, mmdb::Residue *second,
                                               mmdb::Manager *mol) const;
      bool are_linked_in_order(mmdb::Residue *first,
                               mmdb::Residue *second,
                               mmdb::Link *link) const;

      std::pair<std::string, bool>
      find_glycosidic_linkage_type_with_order_switch(mmdb::Residue *first, mmdb::Residue *second) const;

      // can throw a std::runtime_error
      chem_link get_chem_link(const std::string &link_id) const;

      void print_chem_links() const;
      static int chiral_volume_string_to_chiral_sign(const std::string &chiral_vol_string);
      static std::string make_chiral_volume_string(int chiral_sign);

      bool linkable_residue_types_p(const std::string &this_res_type,
                                    const std::string &env_residue_res_type);

      bool OXT_in_residue_restraints_p(const std::string &residue_type) const;

      void read_energy_lib(const std::string &file_name);

      // return HB_UNASSIGNED when not found
      //
      hb_t get_h_bond_type(const std::string &atom_name,
                           const std::string &monomer_name,
                           int imol_enc) const;

      // return HB_UNASSIGNED when not found
      //
      hb_t get_h_bond_type(const std::string &type_energy) const;

      // Find the bonded neighbours of the given atoms - throw an
      // exception if residue name is not in dictionary.
      // 
      std::vector<std::string> get_bonded_neighbours(const std::string &comp_id, int imol_enc,
                                                     const std::string &atom_name_1,
                                                     const std::string &atom_name_2,
                                                     bool also_2nd_order_neighbs_flag=0) const;

      std::vector<std::pair<std::string, std::string> >
      get_bonded_and_1_3_angles(const std::string &comp_id,
                                int imol_enc) const;

      // add a monomer restraints description.
      //
      void add(int imol_enc, const dictionary_residue_restraints_t &rest) {
         std::pair<int, dictionary_residue_restraints_t> restp(imol_enc, rest);
         dict_res_restraints.push_back(restp);
      }

      // a new pdb file has been read in (say).  The residue types
      // have been compared to the dictionary.  These (comp_ids) are
      // the types that are not in the dictionary.  Try to load an
      // ccp4srs description at least so that we can draw their bonds
      // correctly.  Use fill_using_ccp4srs().
      // 
      bool try_load_ccp4srs_description(const std::vector<std::string> &comp_ids);

      // This is made public so that we can decide if we want to set IMOL_ENC_ANY for this cif file
      // 
      bool is_non_auto_load_ligand(const std::string &resname) const;

      // expand to 4c, the atom_id, give that it should match an atom of the type comp_id.
      // Used in chem mods, when we don't know the comp_id until residue modification time.
      // 
      std::string atom_id_expand(const std::string &atom_id,
                                 const std::string &comp_id,
                                 int imol_enc) const;

      // return "" if not found, else return the energy type found in ener_lib.cif
      // 
      std::string get_type_energy(const std::string &atom_name,
                                  const std::string &residue_name,
                                  int imol) const;

      // return -1.1 on failure to look up.
      // 
      double get_vdw_radius(const std::string &atom_name,
                            const std::string &residue_name,
                            int imol,
                            bool use_vdwH_flag) const;

      // calculated once and then stored
      bool atom_is_metal(mmdb::Atom *atom) const;

      // So that we can use semi-sensible distance restraints for metals to waters, Zn-HIS
      // for example
      std::map<std::string, double> metal_O_map;
      std::map<std::string, double> metal_N_map;
      std::map<std::string, double> metal_S_map;
      bool parse_metal_NOS_distance_tables(); // fill the above sets
      // extract values from these sets - return 0.0 on failure
      double get_metal_O_distance(const std::string &metal) const;
      double get_metal_N_distance(const std::string &metal) const;

      // 20230517-PE
      // New style metal distances from Keitaro
      //
      void read_metal_distances(const std::string &file_name);
      // Here is where those distances are stored:
      std::map<std::string, std::vector<metal_ligand_t> > metals_store;
      // debugging:
      void print_metal_store() const; 

      // Find the non-bonded contact distance
      // 
      // Return a pair, if not found the first is 0.  Look up in the energy_lib.
      // 
      // 20151126 We need to know more than just the energy types. We need to know 
      //          if the atoms are in the same ring (if so then we apply 1-4 distance 
      //          corrections). 
      //
      //          We pass in_same_residue_flag because 1-4 nbc-distance shortening can
      //          be applied within a residue, but we leave them full between
      //          residues.  This is a bit of a hack currently. We really need to find
      //          if the atoms are 1-4 related.
      //
      std::pair<bool, double> get_nbc_dist(const std::string &energy_type_1,
                                           const std::string &energy_type_2,
                                           bool in_same_residue_flag = true,
                                           bool in_same_ring_flag = true) const;

      // faster, when the caller has cached the metal state
      std::pair<bool, double> get_nbc_dist_v2(const std::string &energy_type_1,
                                              const std::string &energy_type_2,
                                              const std::string &element_1,
                                              const std::string &element_2,
                                              bool is_metal_atom_1,
                                              bool is_metal_atom_2,
                                              bool extended_atoms_mode, // if the model does not have Hydrogen atoms
                                              bool in_same_residue_flag = true,
                                              bool in_same_ring_flag = true) const;

      energy_lib_atom get_energy_lib_atom(const std::string &ener_type) const;
      
      // Return -1 if residue type not found.
      // 
      int n_non_hydrogen_atoms(const std::string &residue_type);

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
      std::string Get_SMILES_for_comp_id(const std::string &comp_id,  int imol_enc) const;

      // debug
      void debug() const;

      class dreiding_torsion_energy_t {
      public:
         double Vjk;
         double phi0_jk;
         double n_jk;
         dreiding_torsion_energy_t(double Vjk_in, double phi0_jk_in, double n_jk_in) {
            Vjk = Vjk_in;
            phi0_jk = phi0_jk_in;
            n_jk = n_jk_in;
         }
         double E(double phi) const { return 0.5*Vjk*(1.0-cos(n_jk*(phi-phi0_jk))); }
      };

      // can thow a std::runtime_error
      double dreiding_torsion_energy(const std::string &comp_id,
                                     int imol_enc,
                                     mmdb::Atom *atom_1,
                                     mmdb::Atom *atom_2,
                                     mmdb::Atom *atom_3,
                                     mmdb::Atom *atom_4) const;
      double dreiding_torsion_energy(const std::string &comp_id,
                                     int imol_enc,
                                     const atom_quad &quad) const;
      double dreiding_torsion_energy(double phi, int sp_a1, int sp_a2,
                                     const std::string &bond_order,
                                     bool at_1_deloc_or_arom,
                                     bool at_2_deloc_or_arom) const;
      double dreiding_torsion_energy(double phi,
                                     double Vjk, double phi0_jk, double n_jk) const;
      double dreiding_torsion_energy(double phi, dreiding_torsion_energy_t dr) const;
      dreiding_torsion_energy_t dreiding_torsion_energy_params(const std::string &comp_id,
                                                               int imol_enc,
                                                               const atom_quad &quad) const;
      dreiding_torsion_energy_t dreiding_torsion_energy_params(double phi, int sp_a1, int sp_a2,
                                                               const std::string &bond_order,
                                                               bool at_1_deloc_or_arom,
                                                               bool at_2_deloc_or_arom) const;
      // use auto-load if not present
      void use_unimodal_ring_torsion_restraints(int imol, const std::string &res_name, int mmcif_read_number);

      // pass the atom names and the desired torsion value - sigma is not specified
      // by the user.
      void use_unimodal_ring_torsion_restraints(int imol, const std::string &res_name,
                                                const std::vector<atom_name_torsion_quad> &tors_info_vec,
                                                int mmcif_read_number);

      // (list "pseudo-ring-1" (list " C1 " " C2 " " C4 " " C5 ") 0.01)
      // (list "pseudo-ring-2" (list " C2 " " C3 " " C5 " " O5 ") 0.01)
      // (list "pseudo-ring-3" (list " C3 " " C4 " " O5 " " C1 ") 0.01)
      void add_pyranose_pseudo_ring_plane_restraints(const std::string &comp_id, int imol_enc,
                                                     const std::string &plane_id,
                                                     std::vector<std::string> &atom_name_vec,
                                                     double esd);

      // to use improper dihedrals rather than plane restraints, call this
      // after reading restraints. Fills the improper dihedral restraints
      // If you read new restraints, you will need to call this function again
      // and then delete_plane_restraints().
      //
      void all_plane_restraints_to_improper_dihedrals();
      void delete_plane_restraints();

      std::vector<std::pair<std::string, std::string> > get_acedrg_atom_types(const std::string &comp_id,
                                                                              int imol_enc) const;

#ifdef HAVE_CCP4SRS
      match_results_t residue_from_best_match(mmdb::math::Graph &graph1, mmdb::math::Graph &graph2,
                                              mmdb::math::GraphMatch &match, int n_match,
                                              ccp4srs::Monomer *monomer_p) const;
      std::vector<match_results_t>
      compare_vs_ccp4srs(mmdb::math::Graph *graph_1, float similarity, int n_vertices,
                         int srs_idx_start, int srs_idx_end,
                         bool fill_graph_matches) const;
      int ccp4_srs_n_entries() const;

      // return empty string if not available.
      std::vector<std::string> get_available_ligand_comp_id(const std::string &hoped_for_head,
                                                            unsigned int n_top=10) const;

#endif // HAVE_CCP4SRS

   };

} // namespace coot

#endif //  PROTEIN_GEOMETRY_HH
