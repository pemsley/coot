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

#include <mmdb/mmdb_manager.h>


namespace coot {

   namespace util { 
   
      class phi_psi_t {
      
	 //
	 double phi_;
	 double psi_;
	 std::string lab;
	 int residue_number;
	 std::string ins_code;
	 //    enum residue_type {ALA,GLY,SER,THR,ASP,GLU,ASN,GLN, 
	 // 		      LEU,ILE,PHE,TYR,HIS,CYS,MET,TRP,
	 // 		      ARG,LYS,PRO,VAL,CYH};
	 std::string residue_name_;
	 bool is_filled_;
      
      public:
	 phi_psi_t(double a, double b,
		   const std::string &res_name,
		   const std::string &residue_label,
		   int resno,
		   const std::string ins_code_in,
		   std::string chainid) {
	    phi_ = a;
	    psi_ = b;
	    lab = residue_label;
	    residue_name_ = res_name;
	    residue_number = resno;
	    ins_code = ins_code_in;
	    chain_id = chainid;
	    is_filled_ = 1;
	 }
	 phi_psi_t() {
	    is_filled_ = 0;
	 };
	 // this can throw an exception (e.g. bonding atoms too far
	 // apart).  Uses get_phi_psi() below
	 phi_psi_t(CResidue *prev, CResidue *this_res, CResidue *next); 
      
	 double phi() const {return phi_;}
	 double psi() const {return psi_;}
	 std::string label() const {return lab;}
	 std::string residue_name() const { return residue_name_; }
	 std::string chain_id;
	 bool is_filled() const {
	    return is_filled_;
	 }
	 friend std::ostream& operator<<(std::ostream &a, phi_psi_t v);
      };

      std::ostream& operator<<(std::ostream &s, phi_psi_t v);

      // throw an exception on failure to get angles or nSelResidues is not 3.
      phi_psi_t ramachandran_angles(PCResidue *SelResidues, int nSelResidues);

      // used by ramachandran_angles:
      std::pair<bool, phi_psi_t> get_phi_psi(PCResidue *SelResidue);
      std::pair<bool, phi_psi_t> get_phi_psi(CResidue *residue_0,
					     CResidue *residue_1,
					     CResidue *residue_2);

      class phi_psi_pair_helper_t {
      public:
	 phi_psi_t first;
	 phi_psi_t second;
	 bool is_valid_pair_flag;
	 phi_psi_pair_helper_t(const phi_psi_t &f, const phi_psi_t &s, bool valid_flag) {
	    first = f;
	    second = s;
	    is_valid_pair_flag = valid_flag;
	 }
	 phi_psi_pair_helper_t() {
	    is_valid_pair_flag = 0;
	 }
      };
   }
}

#endif // COOT_RAMA_HH

