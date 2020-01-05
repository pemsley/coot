/* coot-utils/coot-coord-rama.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2008 by The University of Oxford
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

#include <stdexcept>

#include "utils/coot-utils.hh"
#include "coot-rama.hh"
#include "clipper/core/coords.h"
#include "clipper/core/ramachandran.h"

#include "coot-coord-utils.hh" // for residue_spec_t

std::ostream& coot::util::operator<<(std::ostream &s, coot::util::phi_psi_t v) {

   s << v.label() << " phi=" << v.phi() << ", psi=" << v.psi();
   return s; 

}

coot::util::phi_psi_t
coot::util::ramachandran_angles(mmdb::PResidue *SelResidues, int nSelResidues) {

   if (nSelResidues != 3) {
      std::string mess = "EXCEPTION: ramachandran_angles was given ";
      mess += coot::util::int_to_string(nSelResidues);
      mess += " residue";
      if (nSelResidues != 1)
	 mess += "s";
      mess += ", not 3";
      throw std::runtime_error(mess);
   }

   std::pair<bool, coot::util::phi_psi_t> p = coot::util::get_phi_psi(SelResidues);
   if (p.first == 0) {
      std::string mess = "EXCEPTION: failed to get atoms for phi psis.";
      throw std::runtime_error(mess);
   }
   return p.second;
} 

// SelResidue is guaranteed to have 3 residues (there is no protection
// for that in this function).
std::pair<bool, coot::util::phi_psi_with_residues_t>
coot::util::get_phi_psi(mmdb::PResidue *SelResidue) {
   return get_phi_psi(SelResidue[0], SelResidue[1], SelResidue[2]);
} 

// each residue needs to be non-null (there is no protection
// for that in this function).
std::pair<bool, coot::util::phi_psi_with_residues_t>
coot::util::get_phi_psi(mmdb::Residue *residue_0, mmdb::Residue *residue_1, mmdb::Residue *residue_2) {

   bool is_valid_flag = 0;
   bool is_pre_pro = 0;
   coot::util::phi_psi_t phi_psi; // part of the returned value
   int nResidueAtoms;
   mmdb::PPAtom res_selection;
   int natom = 0;
   int ires = residue_1->GetSeqNum();
   clipper::Coord_orth c_prev, n_this, ca_this, c_this, n_next;

   residue_0->GetAtomTable(res_selection, nResidueAtoms);
   if (nResidueAtoms > 0) {
      for (int j=0; j<nResidueAtoms; j++) {
	 std::string atom_name = res_selection[j]->name;
	 if (atom_name == " C  ") {
	    c_prev = clipper::Coord_orth(res_selection[j]->x,
					 res_selection[j]->y,
					 res_selection[j]->z);
	    natom++;
	 }
      }
   }
   residue_1->GetAtomTable(res_selection, nResidueAtoms);
   if (nResidueAtoms > 0) {
      for (int j=0; j<nResidueAtoms; j++) {
	 std::string atom_name = res_selection[j]->name;
	 if (atom_name == " C  ") {
	    c_this = clipper::Coord_orth(res_selection[j]->x,
					 res_selection[j]->y,
					 res_selection[j]->z);
	    natom++;
	 }
	 if (atom_name == " CA ") {
	    ca_this = clipper::Coord_orth(res_selection[j]->x,
					  res_selection[j]->y,
					  res_selection[j]->z);
	    natom++;
	 }
	 if (atom_name == " N  ") {
	    n_this = clipper::Coord_orth(res_selection[j]->x,
					 res_selection[j]->y,
					 res_selection[j]->z);
	    natom++;
	 }
      }
   }

   residue_2->GetAtomTable(res_selection, nResidueAtoms);
   if (std::string(residue_2->GetResName()) == "PRO")
      is_pre_pro = 1;
   if (nResidueAtoms > 0) {
      for (int j=0; j<nResidueAtoms; j++) {
	 std::string atom_name = res_selection[j]->name;
	 if (atom_name == " N  ") {
	    n_next = clipper::Coord_orth(res_selection[j]->x,
					 res_selection[j]->y,
					 res_selection[j]->z);
	    natom++;
	 }
      }
   }

   if (natom == 5) {
      char num[30];
      snprintf(num,20,"%d",ires); 
      std::string label(num);
      std::string segid = residue_1->GetChainID();
      std::string inscode = residue_1->GetInsCode();
      label += inscode;
      label += " ";
      label += segid;
      label += " ";
      label += residue_1->name;
      
      double phi   = clipper::Util::rad2d(ca_this.torsion(c_prev, n_this, ca_this, c_this));
      double psi   = clipper::Util::rad2d(ca_this.torsion(n_this, ca_this, c_this, n_next));
      
      phi_psi = coot::util::phi_psi_t(phi, psi,
				      residue_1->name,
				      label.c_str(),
				      ires,
				      inscode,
                                      segid,
                                      is_pre_pro);
      // peptide bonding atoms have to be within 2.0A, or this is not
      // a valid peptide.
      // 
      double dist_1 = clipper::Coord_orth::length(c_prev, n_this);
      double dist_2 = clipper::Coord_orth::length(c_this, n_next);

      if (dist_1 < 2.0) 
	 if (dist_2 < 2.0) 
	    is_valid_flag = 1;
      
   } else {
      // std::cout << "only found " << natom << " atoms " << std::endl;
   }

   coot::util::phi_psi_with_residues_t phi_psi_with_residues(phi_psi);
   phi_psi_with_residues.residue_prev = residue_0;
   phi_psi_with_residues.residue_this = residue_1;
   phi_psi_with_residues.residue_next = residue_2;

   return std::pair<bool, coot::util::phi_psi_with_residues_t> (is_valid_flag, phi_psi_with_residues);
}


// this can throw an exception (e.g. bonding atoms too far
// apart).
coot::util::phi_psi_t::phi_psi_t(mmdb::Residue *prev, mmdb::Residue *this_res, mmdb::Residue *next) {

   std::pair<bool, coot::util::phi_psi_with_residues_t> bpp = coot::util::get_phi_psi(prev, this_res, next);

   if (! bpp.first) {
      std::string mess = "bad residues for phi,psi calculation";
      throw std::runtime_error(mess);
   } else { 
      *this = bpp.second;
   }
}
