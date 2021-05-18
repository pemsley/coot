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


#ifndef ENERGY_LIB_HH
#define ENERGY_LIB_HH

#include <string>
#include <vector>
#include <map>

#include <mmdb2/mmdb_manager.h>
#include "hb-types.hh"

namespace coot {

   // ------------------------------------------------------------------------
   //                  energy lib
   // ------------------------------------------------------------------------

   class energy_lib_atom {
   public:
      std::string type;
      mmdb::realtype weight;
      hb_t hb_type;
      // radii are negative if not assigned.
      mmdb::realtype vdw_radius;
      mmdb::realtype vdwh_radius;
      mmdb::realtype ion_radius;
      std::string element;
      int valency; // negative if unset
      int sp_hybridisation; // negative if unset
      energy_lib_atom(const std::string &type_in,
                      hb_t hb_type_in,
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
      bool needed_permissive;
      energy_lib_bond() {
         type = "unset";
         length = 0;
         esd = 0;
         needed_permissive = false;
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
         needed_permissive = false;
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
      bool filled() const {
         return (type != "unset");
      }
      void set_needed_permissive() {
         needed_permissive = true;
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
//                std::cout << "looking at \"" << atom_type_1 << "\" and \"" << atom_type_3
//                          << "\"" << std::endl;
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

      void read(const std::string &file_name,
                bool print_info_message_flag=false);
      void add_energy_lib_atom(    const energy_lib_atom    &atom);
      void add_energy_lib_bond(    const energy_lib_bond    &bond);
      void add_energy_lib_angle(   const energy_lib_angle   &angle);
      void add_energy_lib_torsion(const energy_lib_torsion &torsion);
      void add_energy_lib_atoms( mmdb::mmcif::PLoop mmCIFLoop);
      void add_energy_lib_bonds( mmdb::mmcif::PLoop mmCIFLoop);
      void add_energy_lib_angles(mmdb::mmcif::PLoop mmCIFLoop);
      void add_energy_lib_torsions(mmdb::mmcif::PLoop mmCIFLoop);
   };
}


#endif // ENERGY_LIB_HH

