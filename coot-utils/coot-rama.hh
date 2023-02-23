/* coot-utils/coot-coord-rama.hh
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2008 by The University of Oxford
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


#ifndef COOT_RAMA_HH
#define COOT_RAMA_HH

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include <mmdb2/mmdb_manager.h>


namespace coot {

   namespace util {

      class phi_psi_t {

         //
         double phi_;
         double psi_;
         std::string lab;
         //    enum residue_type {ALA,GLY,SER,THR,ASP,GLU,ASN,GLN,
         //                       LEU,ILE,PHE,TYR,HIS,CYS,MET,TRP,
         //                       ARG,LYS,PRO,VAL,CYH};
         std::string residue_name_;
         bool is_filled_;
         bool is_pre_pro_;
         bool is_allowed_is_set_;
         bool is_allowed_;

      public:
         // torsion angles in degrees
         phi_psi_t(double phi_in, double psi_in,
                   const std::string &res_name,
                   const std::string &residue_label,
                   int resno,
                   const std::string &ins_code_in,
                   const std::string &chain_id_in) :
            lab(residue_label), residue_name_(res_name), chain_id(chain_id_in), ins_code(ins_code_in) {
            phi_ = phi_in;
            psi_ = psi_in;
            residue_number = resno;
            is_filled_ = 1;
            is_pre_pro_ = false;
            is_allowed_is_set_ = false;
            is_allowed_ = false;
         }
         phi_psi_t(double a, double b,
                   const std::string &res_name,
                   const std::string &residue_label,
                   int resno,
                   const std::string &ins_code_in,
                   const std::string &chain_id_in,
                   bool is_pre_pro) :
            lab(residue_label), residue_name_(res_name), chain_id(chain_id_in), ins_code(ins_code_in) {
            phi_ = a;
            psi_ = b;
            residue_number = resno;
            is_filled_ = 1;
            is_pre_pro_ = is_pre_pro;
            is_allowed_is_set_ = false;
            is_allowed_ = false;
         }

         phi_psi_t() {
            is_filled_ = false;
            is_pre_pro_ = false;
            phi_ = 0.0;
            psi_ = 0.0;
            residue_number = -1;
            is_allowed_is_set_ = false;
            is_allowed_ = false;
         };
         // this can throw an exception (e.g. bonding atoms too far
         // apart).  Uses get_phi_psi() below
         phi_psi_t(mmdb::Residue *prev, mmdb::Residue *this_res, mmdb::Residue *next);

         double phi() const {return phi_;}
         double psi() const {return psi_;}
         std::string label() const {return lab;}
         std::string residue_name() const { return residue_name_; }
         // why not use a residue spec here?
         std::string chain_id;
         int residue_number;  // rename this res_no
         std::string ins_code;
         bool is_filled() const { return is_filled_; }
         bool is_allowed_is_set() const { return is_allowed_is_set_; }
         bool is_allowed() const { return is_allowed_; }
         bool is_pre_pro() const { return is_pre_pro_; }
         void set_is_allowed(bool state) {
            is_allowed_is_set_ = true;
            is_allowed_ = state;
         }
         friend std::ostream& operator<<(std::ostream &a, phi_psi_t v);
      };

      std::ostream& operator<<(std::ostream &s, phi_psi_t v);

      // throw an exception on failure to get angles or nSelResidues is not 3.
      phi_psi_t ramachandran_angles(mmdb::PResidue *SelResidues, int nSelResidues);

      class phi_psi_with_residues_t : public phi_psi_t {
      public:
         mmdb::Residue *residue_prev;
         mmdb::Residue *residue_this;
         mmdb::Residue *residue_next;
         phi_psi_with_residues_t(mmdb::Residue *r_1,
                                 mmdb::Residue *r_2,
                                 mmdb::Residue *r_3) : phi_psi_t(r_1, r_2, r_3) {
            residue_prev = r_1;
            residue_this = r_2;
            residue_next = r_3;
         }
         explicit phi_psi_with_residues_t(const phi_psi_t &pp) : phi_psi_t(pp) {
            residue_prev = 0;
            residue_this = 0;
            residue_next = 0;
         }
         phi_psi_with_residues_t() : phi_psi_t() {
            residue_prev = 0;
            residue_this = 0;
            residue_next = 0;
         }
      };

      // used by ramachandran_angles:
      std::pair<bool, phi_psi_with_residues_t> get_phi_psi(mmdb::PResidue *SelResidue);
      std::pair<bool, phi_psi_with_residues_t> get_phi_psi(mmdb::Residue *residue_0,
                                                           mmdb::Residue *residue_1,
                                                           mmdb::Residue *residue_2);

      class phi_psi_pair_helper_t {
      public:
         phi_psi_t first;
         phi_psi_t second;
         bool is_valid_pair_flag;
         phi_psi_pair_helper_t(const phi_psi_t &f, const phi_psi_t &s, bool valid_flag) : first(f), second(s) {
            is_valid_pair_flag = valid_flag;
         }
         phi_psi_pair_helper_t() {
            is_valid_pair_flag = 0;
         }
      };
   }
}

#endif // COOT_RAMA_HH

