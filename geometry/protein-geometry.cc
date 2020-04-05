/* geometry/protein-geometry.cc
 * 
 * Copyright 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
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

#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>  // needed for sort? Yes.
#include <stdexcept>  // Thow execption.

#include "utils/win-compat.hh"
#include "mini-mol/atom-quads.hh"
#include "geometry/protein-geometry.hh"
#include "utils/coot-utils.hh"

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define DATADIR "C:/coot/share"
#define PKGDATADIR DATADIR
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

#include "clipper/core/clipper_util.h"

#include "compat/coot-sysdep.h"

#include "lbg-graph.hh"

// std::string 
// coot::basic_dict_restraint_t::atom_id_1_4c() const {
//    return atom_id_mmdb_expand(atom_id_1_); 
// }

// std::string 
// coot::basic_dict_restraint_t::atom_id_2_4c() const {
//    return atom_id_mmdb_expand(atom_id_2_);
// }

// std::string 
// coot::basic_dict_restraint_t::atom_id_mmdb_expand(const std::string &atomname) const {

//    std::string r;
//    int ilen = atomname.length();

//    if (ilen == 4) return atomname;
   
//    if (ilen == 1) {
//       r = " ";
//       r += atomname;
//       r += "  ";
//    } else {
//       if (ilen == 2) { 
// 	 r = " ";
// 	 r += atomname;
// 	 r += " ";
//       } else {
// 	 if (ilen == 3) {
// 	    r = " ";
// 	    r += atomname;
// 	 } else {
// 	    r = atomname;
// 	 }
//       }
//    }
//    return r;
// }

std::string
coot::atom_id_mmdb_expand(const std::string &atomname) { 
   std::string r;
   int ilen = atomname.length();
      
   if (ilen == 4) return atomname;
      
   if (ilen == 1) {
      r = " ";
      r += atomname;
      r += "  ";
   } else {
      if (ilen == 2) {

	 // 20180512-PE we have to be more clever here for metals.
	 // But what about CA! - argh! we shouldn't be using this function.
	 // We need to know the residue name to pad correctly.
	 //
	 bool done = false;
	 if (atomname == "MG" || atomname == "NA" || atomname == "LI" || atomname == "LI" || atomname == "AL" || atomname == "SI" ||
	     atomname == "CL" || atomname == "SC" || atomname == "TI" || atomname == "CR" || atomname == "MN" || atomname == "FE" ||
	     atomname == "CO" || atomname == "NI" || atomname == "CU" || atomname == "ZN" || atomname == "GA" || atomname == "AS" ||
	     atomname == "SE" || atomname == "BR" || atomname == "RB" || atomname == "SR" || atomname == "RE" || atomname == "OS" ||
	     atomname == "IR" || atomname == "PT" || atomname == "AU" || atomname == "HG" || atomname == "PB" || atomname == "BI") {
	    r += atomname;
	    r += "  ";
	 } else {
	    r = " ";
	    r += atomname;
	    r += " ";
	 }
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

std::string
coot::atom_id_mmdb_expand(const std::string &atomname, const std::string &element) {

   std::string r = coot::atom_id_mmdb_expand(atomname);

   if (element.length() == 2 && element[0] != ' ') {
      if (atomname.length() == 1) { // unlikely
	 r = " ";
	 r += atomname;
	 r += "  ";
      } else {
	 if (atomname.length() == 2) {
	    r = atomname;
	    r += "  ";
	 } else {
	    if (atomname.length() == 3) {
	       r = atomname;
	       r += " ";
	    } else {
	       r = atomname;
	    }
	 }
      }
   }
   if (0)  // debug
      std::cout << "Given :" << atomname << ": and element :" <<
	 element << ": returning :" << r << ":" << std::endl;
   return r;
}


bool
coot::protein_geometry::matches_imol(int imol_dict, int imol_enc) const {

   if (imol_dict == IMOL_ENC_ANY)
      return true;
   else
      return (imol_dict == imol_enc);
}


// for mmdb::math::Graph mmdb::math::Edge usage
//
// if the bond type is "deloc" then return a single bond.  This
// doesn't matter (at the moment) because the code using this doesn't
// care about the bond order, it only cares about connectivity (or so
// I think).
// 
int
coot::dict_bond_restraint_t::mmdb_bond_type() const {
   int bt = 1;
   if (type_ == "double")
      bt = 2;
   if (type_ == "triple")
      bt = 3;
   return bt;
}


void
coot::protein_geometry::assign_chiral_volume_targets() {

   for (unsigned int idict=0; idict<dict_res_restraints.size(); idict++) {
      if (dict_res_restraints[idict].second.has_unassigned_chiral_volumes()) {
	 if (0) 
	    std::cout << "DEBUG:: assign_chiral_volume_targets for dict_res_restraints entry: "
		      << idict << " "
		      << dict_res_restraints[idict].second.residue_info.comp_id
		      << " has unassigned chiral volumes" << std::endl;
	 dict_res_restraints[idict].second.assign_chiral_volume_targets();
      }
   }

   assign_link_chiral_volume_targets();
}

// the chiral restraint for this comp_id(s) may need filtering
// (i.e. removing some of them if they are not real chiral centres
// (e.g. from prodrg restraints)).
void
coot::protein_geometry::filter_chiral_centres(int imol, const std::vector<std::string> &comp_ids_for_filtering) {

   for (unsigned int ichir=0; ichir<comp_ids_for_filtering.size(); ichir++) {
      bool minimal = false;
      int idx = get_monomer_restraints_index(comp_ids_for_filtering[ichir], imol, minimal);
      if (idx != -1) {
	 const coot::dictionary_residue_restraints_t &restraints =
	    dict_res_restraints[idx].second;
	 std::vector<coot::dict_chiral_restraint_t> new_chirals =
	    filter_chiral_centres(restraints);
	 dict_res_restraints[idx].second.chiral_restraint = new_chirals;
      }
   }
}

// Return a filtered list, that is don't include chiral centers that
// are connected to more than one hydrogen.
// 
std::vector<coot::dict_chiral_restraint_t>
coot::protein_geometry::filter_chiral_centres(const dictionary_residue_restraints_t &restraints) {

   std::vector<coot::dict_chiral_restraint_t> v;
   for (unsigned int ichir=0; ichir<restraints.chiral_restraint.size(); ichir++) { 
      int n_H=0;
      for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
	 if (restraints.bond_restraint[ib].atom_id_1_4c() ==
	     restraints.chiral_restraint[ichir].atom_id_c_4c()) {
	    if (restraints.element(restraints.bond_restraint[ib].atom_id_2_4c()) == " H")
	       n_H++;
	 }
	 if (restraints.bond_restraint[ib].atom_id_2_4c() ==
	     restraints.chiral_restraint[ichir].atom_id_c_4c()) {
	    if (restraints.element(restraints.bond_restraint[ib].atom_id_1_4c()) == " H")
	       n_H++;
	 }
      }
      if (n_H <= 1)
	 v.push_back(restraints.chiral_restraint[ichir]);
   }
   return v;
} 

void
coot::protein_geometry::assign_link_chiral_volume_targets() {

   for (unsigned int idict=0; idict<dict_link_res_restraints.size(); idict++) {
      if (dict_link_res_restraints[idict].has_unassigned_chiral_volumes()) {
	 dict_link_res_restraints[idict].assign_link_chiral_volume_targets();
      }
   }
}

bool
coot::dict_plane_restraint_t::matches_names(const coot::dict_plane_restraint_t &r) const {

   bool status = true;
   unsigned int n_found = 0;
   if (atom_ids.size() != r.atom_ids.size())
      return false;
   if (atom_ids.size() > 0)
      status = false; // initial setting.
   
   for (unsigned int i=0; i<atom_ids.size(); i++) {
      const std::string &ref_atom = atom_ids[i].first;
      for (unsigned int j=0; j<r.atom_ids.size(); j++) { 
	 if (r.atom_ids[j].first == ref_atom) {
	    n_found++;
	    break;
	 }
      }
   }
   if (n_found == atom_ids.size())
      status = true;
   return status;
}


// constructor.  Caller should make sure that there are no bonds
// before constructing this (mol->RemoveBonds());
//
coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t(mmdb::Residue *residue_p) {

   init(residue_p);
}

void
coot::dictionary_residue_restraints_t::init(mmdb::Residue *residue_p) {

   filled_with_bond_order_data_only_flag = false;
   if (residue_p) {
      mmdb::PModel   model;
      mmdb::PChain   chain;
      // mmdb::PResidue res = 0;
      mmdb::math::Graph    graph;
      mmdb::math::PPVertex V;
      mmdb::math::PPEdge   E;
      int       i, im,ic,ir, nV,nE, k1,k2;

      graph.MakeGraph   ( residue_p,NULL );
      graph.GetVertices ( V,nV );
      graph.GetEdges    ( E,nE );

      mmdb::PPAtom residue_atoms = 0;
      int nResidueAtoms;
      residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
	 
      if (0) { //debug
	 for (int iv=0; iv<nV; iv++) { 
	    std::cout << "vertex " << iv << " of " << nV << " " << V[iv] << std::endl;
	 }
			      
	 for (int ie=0; ie<nE; ie++) { 
	    std::cout << "edge " << ie << " of " << nE << " " << E[ie] << std::endl;
	    std::cout << "edge " << ie << " of " << nE << " with vertex1 "
		      << E[ie]->GetVertex1() << " and vertex2 "
		      << E[ie]->GetVertex2() << std::endl;
	 }
      }
			   
      for (i=0;i<nE;i++)  {

	 // mmdb 1.25.3 on pc offset needs to be -1
	 int idx_offset = -1;

	 bool debug = false;
	 
	 if (debug) { 
	    std::cout << "V index for k1 " << E[i]->GetVertex1() << std::endl;
	    std::cout << "V index for k2 " << E[i]->GetVertex2() << std::endl;
	 }
	 // this indexing needs testing.
	 k1 = V[E[i]->GetVertex1()+idx_offset]->GetUserID();
	 k2 = V[E[i]->GetVertex2()+idx_offset]->GetUserID();
	 if (debug) 
	    std::cout << "1 adding bond to atom  " << k1 << " of " << nResidueAtoms << std::endl;
	 residue_p->atom[k1]->AddBond ( residue_p->atom[k2],E[i]->GetType() );
	 if (debug)
	    std::cout << "2 adding bond to atom  " << k2 << " of " << nResidueAtoms << std::endl;
	 residue_p->atom[k2]->AddBond ( residue_p->atom[k1],E[i]->GetType() );
	 if (debug) 
	    std::cout << "added bond of type " << E[i]->GetType() << " to "
		      << k1 << " and " << k2 << std::endl;
      }

	 
      // bool calc_only = true;
      // mol->MakeBonds(calc_only);  // crash, hence above hack.

      std::string comp_id = residue_p->GetResName();
      std::string group("monomer");
      std::string desc_level(".");
      int n_all = nResidueAtoms;
      int n_non_H = 0;
      for (int iat=0; iat<nResidueAtoms; iat++) {
	 std::string ele(residue_atoms[iat]->element);
	 if (ele != "H" && ele != " H" && ele != "D" && ele != " D")
	    n_non_H++;
      }
	    
      residue_info = dict_chem_comp_t(comp_id, comp_id, comp_id, group,
				      n_all, n_non_H, desc_level);
      // also fill atom_info with dict_atom objects
      for (int iat=0; iat<nResidueAtoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 std::string ele(residue_atoms[iat]->element);
	 dict_atom da(at->name, at->name, ele, "", std::pair<bool, float> (false, 0));
	 atom_info.push_back(da);
      }

      std::vector<atom_pair_t> bond_pairs;
      for (int iat=0; iat<nResidueAtoms; iat++) { 
	 mmdb::Atom *at_1 = residue_atoms[iat];
	 if (at_1) { 
	    int n_bonds_1 = at_1->GetNBonds();
	    mmdb::AtomBond *AtomBonds = NULL;
	    int n_bonds_2; 
	    at_1->GetBonds(AtomBonds, n_bonds_2);
	    for (int ibond=0; ibond<n_bonds_2; ibond++) { 
	       mmdb::Atom *at_2 = AtomBonds[ibond].atom;
	       if (at_2) { 
		  if (at_1 < at_2) { // pointer comparison
		     std::string at_name_1(at_1->name);
		     std::string at_name_2(at_2->name);
		     std::string type = "single";
		     if (AtomBonds[ibond].order == 2)
			type = "double";
		     if (AtomBonds[ibond].order == 3)
			type = "triple";
		     clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
		     clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
		     double dist = sqrt((pt_1-pt_2).lengthsq());
		     double dist_esd = 0.02;
		     dict_bond_restraint_t br(at_name_1, at_name_2, type, dist, dist_esd);
		     bond_restraint.push_back(br);
		     bond_pairs.push_back(atom_pair_t(at_1, at_2));
		  }
	       }
	    }
	 }
      }

      // now find angle by finding bond-pair pairs that that share
      // an atom
      // 
      for (unsigned int ibp=0; ibp<bond_pairs.size(); ibp++) { 
	 for (unsigned int jbp=ibp; jbp<bond_pairs.size(); jbp++) {
	    if (ibp != jbp) {
	       mmdb::Atom *shared_atom = bond_pairs[ibp].shared_atom(bond_pairs[jbp]);
	       if (shared_atom) {
		  mmdb::Atom *at_1 = bond_pairs[ibp].at_1;
		  mmdb::Atom *at_2 = bond_pairs[ibp].at_2;
		  mmdb::Atom *at_3 = bond_pairs[jbp].at_1;
		  if (at_1 == shared_atom) {
		     at_1 = bond_pairs[ibp].at_2;
		     at_2 = bond_pairs[ibp].at_1; // shared atom
		  } 
		  if (at_3 == shared_atom)
		     at_3 = bond_pairs[jbp].at_2;

		  // this test should not be needed (because bond_pairs vector is only made
		  // of atoms that are non-null) - but is here to keep the static analyzer quiet.
		  // 
		  if (at_1 && at_2 && at_3) { 

		     clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
		     clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
		     clipper::Coord_orth pt_3(at_3->x, at_3->y, at_3->z);
		     // doesn't exist (mmdb problem)?
		     // double angle = BondAngle(at_1, at_2, at_3);
		     double angle = clipper::Util::rad2d(clipper::Coord_orth::angle(pt_1, pt_2, pt_3));
		     // std::cout << "angle: " << angle << std::endl;
		     if (angle > 0.001) {
			std::string at_name_1(at_1->name);
			std::string at_name_2(at_2->name);
			std::string at_name_3(at_3->name);
			double angle_esd = 3;
			dict_angle_restraint_t ar(at_name_1, at_name_2, at_name_3, angle, angle_esd);
			angle_restraint.push_back(ar);
			if (0)
			   std::cout << "found an angle restraint "
				     << at_name_1 << " " << at_name_2 << " " << at_name_3 << std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
}

// mol contains one residue in a hierarchy, the residue from which the
// dictionary should be constructed.
// 
coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t(mmdb::Manager *mol) {

   mmdb::Residue *residue_p = NULL;
   filled_with_bond_order_data_only_flag = true; // it has nothing initially
   
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) {
	 residue_p = chain_p->GetResidue(ires);
	 if (residue_p)
	    break;
      }
      if (residue_p)
	 break;
   }
   

   if (residue_p) {
      mol->RemoveBonds();
      init(residue_p);
   }
} 



void
coot::dictionary_residue_restraints_t::clear_dictionary_residue() {

   residue_info = coot::dict_chem_comp_t("", "", "", "", 0, 0, "");
   has_partial_charges_flag = 0;

   // need different constructors.
//    atom_info.resize(0);
//    bond_restraint.resize(0);
//    angle_restraint.resize(0);
//    torsion_restraint.resize(0);
//    chiral_restraint.resize(0);
//    plane_restraint.resize(0);
}


// Because they were all at origin, for example.
// 
void
coot::protein_geometry::delete_atom_positions(const std::string &comp_id,
					      int imol_enc,
					      int pos_type) {
   
   // for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
   // if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {

   int idr = get_monomer_restraints_index(comp_id, imol_enc, false);
   if (idr >= 0) {
      for (unsigned int iat=0; iat<dict_res_restraints[idr].second.atom_info.size(); iat++) { 
	 if (pos_type == dict_atom::IDEAL_MODEL_POS)
	    dict_res_restraints[idr].second.atom_info[iat].pdbx_model_Cartn_ideal.first = false;
	 if (pos_type == dict_atom::REAL_MODEL_POS)
	    dict_res_restraints[idr].second.atom_info[iat].model_Cartn.first = false;
      }
   }
}

// the input from_tos should be 4 character atom names. PDBv3 FIXME
void
coot::dictionary_residue_restraints_t::atom_id_swap(const std::vector< std::pair<std::string, std::string> > &from_tos) {

   std::vector<std::pair<int, std::string> > alter_idx;
   // find the set of difference to apply (by filling alter_idx), then apply them

   // -------------- Atoms ---------------------------------
   // 
   for (unsigned int iat=0; iat<atom_info.size(); iat++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (atom_info[iat].atom_id_4c == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(iat, from_tos[j].second));
	 }
      }
   }

   for (unsigned int idx=0; idx<alter_idx.size(); idx++) {
      std::string was = atom_info[alter_idx[idx].first].atom_id_4c;
      atom_info[alter_idx[idx].first].atom_id_4c = alter_idx[idx].second;
      atom_info[alter_idx[idx].first].atom_id    = util::remove_whitespace(alter_idx[idx].second);
   } 
   alter_idx.clear();

   // -------------- Atom Tree ------------------------------
   // 
   for (unsigned int iat=0; iat<tree.size(); iat++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (tree[iat].atom_id == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(iat, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      tree[alter_idx[idx].first].atom_id = alter_idx[idx].second;
   alter_idx.clear();
   
   
   // -------------- Bonds 1st atom --------------------------
   //
   for (unsigned int ib=0; ib<bond_restraint.size(); ib++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (bond_restraint[ib].atom_id_1_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ib, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      bond_restraint[alter_idx[idx].first].set_atom_1_atom_id(util::remove_whitespace(alter_idx[idx].second));
   alter_idx.clear();
   
   // -------------- Bonds 2nd atom --------------------------
   //
   for (unsigned int ib=0; ib<bond_restraint.size(); ib++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (bond_restraint[ib].atom_id_2_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ib, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      bond_restraint[alter_idx[idx].first].set_atom_2_atom_id(util::remove_whitespace(alter_idx[idx].second));
   alter_idx.clear();

   
   // -------------- Angles 1st atom ------------------------
   //
   for (unsigned int ia=0; ia<angle_restraint.size(); ia++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (angle_restraint[ia].atom_id_1_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ia, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      angle_restraint[alter_idx[idx].first].set_atom_1_atom_id(alter_idx[idx].second);
   alter_idx.clear();

   // -------------- Angles 2nd atom ------------------------
   //
   for (unsigned int ia=0; ia<angle_restraint.size(); ia++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (angle_restraint[ia].atom_id_2_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ia, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      angle_restraint[alter_idx[idx].first].set_atom_2_atom_id(alter_idx[idx].second);
   alter_idx.clear();

   // -------------- Angles 3rd atom ------------------------
   //
   for (unsigned int ia=0; ia<angle_restraint.size(); ia++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (angle_restraint[ia].atom_id_3_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ia, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      angle_restraint[alter_idx[idx].first].set_atom_3_atom_id(alter_idx[idx].second);
   alter_idx.clear();

   // -------------- Torsion 1st atom ------------------------
   //
   for (unsigned int it=0; it<torsion_restraint.size(); it++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (torsion_restraint[it].atom_id_1_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(it, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      torsion_restraint[alter_idx[idx].first].set_atom_1_atom_id(alter_idx[idx].second);
   alter_idx.clear();
   

   // -------------- Torsion 2nd atom ------------------------
   //
   for (unsigned int it=0; it<torsion_restraint.size(); it++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (torsion_restraint[it].atom_id_2_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(it, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      torsion_restraint[alter_idx[idx].first].set_atom_2_atom_id(alter_idx[idx].second);
   alter_idx.clear();

   // -------------- Torsion 3rd atom ------------------------
   //
   for (unsigned int it=0; it<torsion_restraint.size(); it++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (torsion_restraint[it].atom_id_3_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(it, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      torsion_restraint[alter_idx[idx].first].set_atom_3_atom_id(alter_idx[idx].second);
   alter_idx.clear();

   // -------------- Torsion 4th atom ------------------------
   //
   for (unsigned int it=0; it<torsion_restraint.size(); it++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (torsion_restraint[it].atom_id_4_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(it, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      torsion_restraint[alter_idx[idx].first].set_atom_4_atom_id(alter_idx[idx].second);
   alter_idx.clear();

   // -------------- Chiral Centre atom ------------------------
   //
   for (unsigned int ic=0; ic<chiral_restraint.size(); ic++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (chiral_restraint[ic].atom_id_c_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ic, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      chiral_restraint[alter_idx[idx].first].set_atom_c_atom_id(alter_idx[idx].second);
   alter_idx.clear();

   // -------------- Chiral 1st atom ------------------------
   //
   for (unsigned int ic=0; ic<chiral_restraint.size(); ic++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (chiral_restraint[ic].atom_id_1_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ic, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      chiral_restraint[alter_idx[idx].first].set_atom_1_atom_id(alter_idx[idx].second);
   alter_idx.clear();
   
   // -------------- Chiral 2nd atom ------------------------
   //
   for (unsigned int ic=0; ic<chiral_restraint.size(); ic++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (chiral_restraint[ic].atom_id_2_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ic, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      chiral_restraint[alter_idx[idx].first].set_atom_2_atom_id(alter_idx[idx].second);
   alter_idx.clear();
   
   // -------------- Chiral 3rd atom ------------------------
   //
   for (unsigned int ic=0; ic<chiral_restraint.size(); ic++) { 
      for (unsigned int j=0; j<from_tos.size(); j++) { 
	 if (chiral_restraint[ic].atom_id_3_4c() == from_tos[j].first) {
	    alter_idx.push_back(std::pair<int, std::string>(ic, from_tos[j].second));
	 }
      }
   }
   for (unsigned int idx=0; idx<alter_idx.size(); idx++)
      chiral_restraint[alter_idx[idx].first].set_atom_3_atom_id(alter_idx[idx].second);
   alter_idx.clear();
   
   // -------------- Planes ------------------------
   // 
   for (unsigned int ip=0; ip<plane_restraint.size(); ip++) { 
      dict_plane_restraint_t &pr = plane_restraint[ip];
      alter_idx.clear();
      for (int iat=0; iat<pr.n_atoms(); iat++) { 
	 for (unsigned int j=0; j<from_tos.size(); j++) {
	    if (pr.atom_id(iat) == from_tos[j].first)
	       alter_idx.push_back(std::pair<int, std::string>(iat, from_tos[j].second));
	 }
      }
      if (alter_idx.size())
	 pr.set_atom_ids(alter_idx);
   }
}


// If include_hydrogen_torsions_flag is set, we check the neighbours
// of the atom-2 and atom-3 to see if this is really a pure hydrogen
// torsion.
// 
std::vector<coot::dict_torsion_restraint_t>
coot::dictionary_residue_restraints_t::get_non_const_torsions(bool include_hydrogen_torsions_flag) const {

   std::vector<coot::dict_torsion_restraint_t> v;
   for (unsigned int i=0; i<torsion_restraint.size(); i++) {
      if (! torsion_restraint[i].is_const()) {
	 if (include_hydrogen_torsions_flag) { 
	    v.push_back(torsion_restraint[i]);
	 } else {
	    // only add this torsion if neither of end atoms of the torsion are hydrogen.
	    if (!is_hydrogen(torsion_restraint[i].atom_id_1())) { 
	       if (!is_hydrogen(torsion_restraint[i].atom_id_4())) { 
		  v.push_back(torsion_restraint[i]);
	       } else {

		  // OK, so atom-4 was a hydrogen in the dictionary,
		  // but is there an atom attached to atom_id_3 that
		  // is not a hydrogen and not atom_id_3?  (Is so,
		  // then this is not a pure hydrogen torsion, and we
		  // can add it to the list).
		  //

		  std::vector<std::string> v_n = neighbours(torsion_restraint[i].atom_id_3(), false);
		  for (unsigned int i_neighb=0; i_neighb<v_n.size(); i_neighb++) { 
		     if (v_n[i_neighb] != torsion_restraint[i].atom_id_4()) { 
			if (v_n[i_neighb] != torsion_restraint[i].atom_id_2()) { 
			   if (v_n[i_neighb] != torsion_restraint[i].atom_id_1()) {
			      if (! is_hydrogen(v_n[i_neighb])) {
				 v.push_back(torsion_restraint[i]);
			      }
			   }
			}
		     }
		  }
	       }
	    } else {

	       // OK, so atom-1 was a hydrogen in the dictionary,
	       // but is there an atom attached to atom_id_3 that
	       // is not a hydrogen and not atom_id_3?  (Is so,
	       // then this is not a pure hydrogen torsion, and we
	       // can add it to the list).
	       //
	       std::vector<std::string> v_n = neighbours(torsion_restraint[i].atom_id_2(), false);
	       for (unsigned int i_neighb=0; i_neighb<v_n.size(); i_neighb++) { 
		  if (v_n[i_neighb] != torsion_restraint[i].atom_id_1()) { 
		     if (v_n[i_neighb] != torsion_restraint[i].atom_id_3()) { 
			if (v_n[i_neighb] != torsion_restraint[i].atom_id_4()) {
			   if (! is_hydrogen(v_n[i_neighb])) {
			      v.push_back(torsion_restraint[i]);
			   }
			}
		     }
		  }
	       }
	    } 
	 }
      }
   }
   return v;
}

bool
coot::dict_torsion_restraint_t::is_pyranose_ring_torsion(const std::string &comp_id) const {

   // Needs fixup for PDBv3
   bool status = false;
   std::string ring_atoms[6] = { " C1 ", " C2 ", " C3 ", " C4 ", " C5 ", " O5 " };
   if (comp_id == "XYP")
      for (unsigned int i=0; i<6; i++)
	 ring_atoms[i][3] = 'B'; // danger on PDBv3 fixup.

   int n_matches = 0;
   for (unsigned int i=0; i<6; i++) { 
      if (atom_id_2_4c() == ring_atoms[i])
	 n_matches++;
      if (atom_id_3_4c() == ring_atoms[i])
	 n_matches++;
   }
   if (n_matches == 2)
      status = true;
   return status;
}

bool
coot::dict_link_torsion_restraint_t::is_pyranose_ring_torsion() const {

   // Needs fixup for PDBv3
   bool status = false;
   std::string ring_atoms[6] = { " C1 ", " C2 ", " C3 ", " C4 ", " C5 ", " O5 " };

   int n_matches = 0;
   for (unsigned int i=0; i<6; i++) { 
      if (atom_id_2_4c() == ring_atoms[i])
	 n_matches++;
      if (atom_id_3_4c() == ring_atoms[i])
	 n_matches++;
   }
   if (n_matches == 2)
      status = true;
   return status;
} 

bool
coot::dict_torsion_restraint_t::is_ring_torsion(const std::vector<std::vector<std::string> > &ring_atoms_sets) const {

   bool match = false; 
   std::vector<std::string> torsion_restraint_atom_names(2);
   torsion_restraint_atom_names[0] = atom_id_2_4c();
   torsion_restraint_atom_names[1] = atom_id_3_4c();
   
   for (unsigned int iring=0; iring<ring_atoms_sets.size(); iring++) { 
      const std::vector<std::string> &ring_atom_names = ring_atoms_sets[iring];

      int n_match = 0;
      for (unsigned int iname_1=0; iname_1<ring_atom_names.size(); iname_1++) {
	 for (unsigned int iname_2=0; iname_2<torsion_restraint_atom_names.size(); iname_2++) { 
	    if (ring_atom_names[iname_1] == torsion_restraint_atom_names[iname_2])
	       n_match++;
	 }
      }
      if (n_match == 2) {
	 match = true;
	 break;
      }
   }
   return match;
} 


bool
coot::dict_torsion_restraint_t::is_const() const {

   bool const_flag = 0;
   if (id_.length() > 5) {
      std::string bit = id_.substr(0,5);
      if (bit == "CONST")
	 const_flag = 1;
      if (bit == "const")
	 const_flag = 1;
   }
   return const_flag;
}


std::string
coot::dictionary_residue_restraints_t::atom_name_for_tree_4c(const std::string &atom_id) const {

   std::string r = atom_id;
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {
      if (atom_info[iat].atom_id == atom_id) {
	 r = atom_info[iat].atom_id_4c;
      }
   }
   return r;
}


// look up the atom id in the atom_info (dict_atom vector).
// Return "" on no atom found with name atom_name;
// 
std::string
coot::dictionary_residue_restraints_t::element(const std::string &atom_name) const {

   std::string r = "";
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {
      if (atom_info[iat].atom_id_4c == atom_name) {
	 r = atom_info[iat].type_symbol;
	 break;
      }
   }
   if (r.length() == 1)
      r = " " + r;
   
   // std::cout << " dictionary_residue_restraints_t::element()"
   // 	     << " on atom name :" << atom_name << ": returns :" << r << ":" << std::endl;
   return r;
}

// likewise look up the energy type.  Return "" on no atom found
// with that atom_name.
// 
std::string
coot::dictionary_residue_restraints_t::type_energy(const std::string &atom_name) const {

   std::string r = "";

   // If you are reading this, then you are looking in a dictionary looked up from an index
   // that is out of bounds.
   // 
   // std::cout << "dictionary_has " << atom_info.size() << " atoms" << std::endl;
   
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {
      if (false)
	 std::cout << "comparing :" << atom_name << ": with :" << atom_info[iat].atom_id_4c
		   << ":" << std::endl;
      if (atom_info[iat].atom_id_4c == atom_name) { // PDBv3 FIXME
	 r = atom_info[iat].type_energy;
	 break;
      }
   }
   return r;
}


std::vector<std::string>
coot::dictionary_residue_restraints_t::neighbours(const std::string &atom_name, bool allow_hydrogen_neighbours_flag) const {

   std::vector<std::string> n;
   for (unsigned int i=0; i<bond_restraint.size(); i++) { 
      if (bond_restraint[i].atom_id_1() == atom_name) {
	 if (allow_hydrogen_neighbours_flag || ! is_hydrogen(bond_restraint[i].atom_id_2())) {
	    n.push_back(bond_restraint[i].atom_id_2());
	 }
      }
      if (bond_restraint[i].atom_id_2() == atom_name) {
	 if (allow_hydrogen_neighbours_flag || ! is_hydrogen(bond_restraint[i].atom_id_1())) {
	    n.push_back(bond_restraint[i].atom_id_1());
	 }
      }
   }
   return n;
}

std::vector<unsigned int>
coot::dictionary_residue_restraints_t::neighbours(unsigned int idx, bool allow_hydrogen_neighbours_flag) const {

   std::vector<unsigned int> v;
   std::string atom_name = atom_info[idx].atom_id_4c;
   for (unsigned int i=0; i<bond_restraint.size(); i++) { 
      if (bond_restraint[i].atom_id_1() == atom_name) {
	 const std::string &other_atom_name = bond_restraint[i].atom_id_2();
	 if (allow_hydrogen_neighbours_flag || ! is_hydrogen(other_atom_name)) {
	    // what is the index of atom_id_2?  - bleugh.  This shows
	    // that the dictionary store should work with atom indices, not
	    // atom names (i.e. do it like RDKit does it).
	    for (unsigned int iat=0; iat<atom_info.size(); iat++) {
	       if (atom_info[iat].atom_id_4c == other_atom_name) {
		  v.push_back(iat);
		  break;
	       }
	    }
	 }
      }
      if (bond_restraint[i].atom_id_2() == atom_name) {
	 const std::string &other_atom_name = bond_restraint[i].atom_id_1();
	 if (allow_hydrogen_neighbours_flag || ! is_hydrogen(other_atom_name)) {
	    // what is the index of atom_id_2?  - bleugh.  This shows
	    // that the dictionary store should work with atom indices, not
	    // atom names (i.e. do it like RDKit does it).
	    for (unsigned int iat=0; iat<atom_info.size(); iat++) {
	       if (atom_info[iat].atom_id_4c == other_atom_name) {
		  v.push_back(iat);
		  break;
	       }
	    }
	 }
      }
   }
   return v;
}





std::vector<std::string>
coot::dictionary_residue_restraints_t::get_attached_H_names(const std::string &atom_name) const {

   std::vector<std::string> v;
   for (unsigned int i=0; i<bond_restraint.size(); i++) { 
      if (bond_restraint[i].atom_id_1() == atom_name) {
	 // if (element(bond_restraint[i].atom_id_2()) == " H")
	 if (is_hydrogen(bond_restraint[i].atom_id_2()))
	    v.push_back(bond_restraint[i].atom_id_2());
      }
      if (bond_restraint[i].atom_id_2() == atom_name) {
	 // if (element(bond_restraint[i].atom_id_1()) == " H")
	 if (is_hydrogen(bond_restraint[i].atom_id_1()))
	    v.push_back(bond_restraint[i].atom_id_1());
      }
   }
   return v;
}


// return an empty string on failure
std::string
coot::dictionary_residue_restraints_t::get_other_H_name(const std::string &H_at_name) const {

   std::string r;

   // if it's a hydrogen atom name as input then the neighbour of that won't be
   // a hydrogen
   //
   std::vector<std::string> neighbs = neighbours(H_at_name, false);

   if (neighbs.size() == 1) {
      const std::string &n = neighbs[0];
      for (unsigned int i=0; i<bond_restraint.size(); i++) {
	 if (bond_restraint[i].atom_id_1() == n) {
	    if (bond_restraint[i].atom_id_2() != H_at_name) {
	       if (false)
		  std::cout << "here 1 with br " << bond_restraint[i] << " ele :"
			    << element(bond_restraint[i].atom_id_2_4c()) << ":" << std::endl;
	       if (element(bond_restraint[i].atom_id_2_4c()) == " H") {
		  r = bond_restraint[i].atom_id_2_4c();
		  break;
	       }
	    }
	 }
	 if (bond_restraint[i].atom_id_2() == n) {
	    if (bond_restraint[i].atom_id_1() != H_at_name) {
	       if (false)
		  std::cout << "here 2 with br " << bond_restraint[i] << " ele :"
			    << element(bond_restraint[i].atom_id_1_4c()) << ":" << std::endl;
	       if (element(bond_restraint[i].atom_id_1_4c()) == " H") {
		  r = bond_restraint[i].atom_id_1_4c();
		  break;
	       }
	    }
	 }
      }
   }
   return r;
}

// return an empty vector on failure
std::vector<std::string>
coot::dictionary_residue_restraints_t::get_other_H_names(const std::string &H_at_name) const {

   std::vector<std::string> v;

   // if it's a hydrogen atom name as input then the neighbour of that won't be
   // a hydrogen
   //
   std::vector<std::string> neighbs = neighbours(H_at_name, false);

   if (neighbs.size() == 1) {
      const std::string &n = neighbs[0];
      for (unsigned int i=0; i<bond_restraint.size(); i++) {
	 if (bond_restraint[i].atom_id_1() == n) {
	    if (bond_restraint[i].atom_id_2() != H_at_name) {
	       if (false)
		  std::cout << "here 1 with br " << bond_restraint[i] << " ele :"
			    << element(bond_restraint[i].atom_id_2_4c()) << ":" << std::endl;
	       if (element(bond_restraint[i].atom_id_2_4c()) == " H") {
		  v.push_back(bond_restraint[i].atom_id_2_4c());
	       }
	    }
	 }
	 if (bond_restraint[i].atom_id_2() == n) {
	    if (bond_restraint[i].atom_id_1() != H_at_name) {
	       if (false)
		  std::cout << "here 2 with br " << bond_restraint[i] << " ele :"
			    << element(bond_restraint[i].atom_id_1_4c()) << ":" << std::endl;
	       if (element(bond_restraint[i].atom_id_1_4c()) == " H") {
		  v.push_back(bond_restraint[i].atom_id_1_4c());
	       }
	    }
	 }
      }
   }
   return v;
}




void 
coot::protein_geometry::add_restraint(std::string comp_id,
				      int imol_enc,
				      const dict_torsion_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   bool ifound = false;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    ifound = true;
	    dict_res_restraints[i].second.torsion_restraint.push_back(restr); 
	    break;
	 }
      }
   }

   // it was not there
   if (! ifound) {
      dictionary_residue_restraints_t drr(comp_id, read_number);
      drr.torsion_restraint.push_back(restr);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, drr);
      dict_res_restraints.push_back(p);
      // add the angle to the newly created dictionary_residue_restraints_t

      // I don't think that we need to do this - they can be pre-added.
      // use back? or reverse iterator? FIXME
      // dict_res_restraints[dict_res_restraints.size()-1].second.torsion_restraint.push_back(restr);
   }
}

void 
coot::protein_geometry::add_restraint(std::string comp_id,
				      int imol_enc,
				      const dict_chiral_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   bool ifound = false;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    ifound = true;
	    dict_res_restraints[i].second.chiral_restraint.push_back(restr); 
	    break;
	 }
      }
   }


   // it was not there
   if (! ifound) {
      std::cout << "---------------------------- oops missing in add_restraint() chiral " << std::endl;
      dictionary_residue_restraints_t drr(comp_id, read_number);
      drr.chiral_restraint.push_back(restr);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, drr);
      dict_res_restraints.push_back(p);
      // add the angle to the newly created dictionary_residue_restraints_t
      // dict_res_restraints[dict_res_restraints.size()-1].
   }
}


// static int
int
coot::protein_geometry::chiral_volume_string_to_chiral_sign(const std::string &volume_sign) {

   int volume_sign_int = coot::dict_chiral_restraint_t::CHIRAL_VOLUME_RESTRAINT_VOLUME_SIGN_UNASSIGNED;
   if (volume_sign.length() > 3) { 

       if (volume_sign.substr(0,3) == "pos") { 
	  volume_sign_int = 1;
       }
       if (volume_sign.substr(0,3) == "neg") { 
	  volume_sign_int = -1;
       }
       if (volume_sign.substr(0,3) == "POS") { 
	  volume_sign_int = 1;
       }
       if (volume_sign.substr(0,3) == "NEG") { 
	  volume_sign_int = -1;
       }
       if (volume_sign == "both" || volume_sign == "BOTH") { 
	  volume_sign_int = coot::dict_chiral_restraint_t::CHIRAL_RESTRAINT_BOTH;
       }
    }
    return volume_sign_int;
}


// and the reverse 
//
// static 
std::string 
coot::protein_geometry::make_chiral_volume_string(int chiral_sign) { 

  std::string s;

  if (chiral_sign == -1) 
    s = "negative";
  if (chiral_sign == 1) 
    s = "positive";
  if (chiral_sign == coot::dict_chiral_restraint_t::CHIRAL_RESTRAINT_BOTH)
    s = "both";
  return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const dict_atom &at) {

   s << "dict_atom: "
     << "atom_id :" << at.atom_id << ":  "
     << "atom-id-4c :" << at.atom_id_4c << ":  "
     << "type-symbol :" << at.type_symbol << ":  "
     << "pdbx_stereo_config: " << at.pdbx_stereo_config.first
     << " \"" << at.pdbx_stereo_config.second << "\" ";
   if (at.formal_charge.first)
      s << "formal-charge " << at.formal_charge.second << " ";
   else 
      s << "no-formal-charge ";
   if (at.partial_charge.first)
      s << "partial-charge " << at.partial_charge.second << " ";
   else 
      s << "no-partial-charge ";
   s << "model-pos " << at.model_Cartn.first << " ";
   if (at.model_Cartn.first)
      s << at.model_Cartn.second.format() << " ";
   s << "ideal-pos " << at.pdbx_model_Cartn_ideal.first << " ";
   if (at.pdbx_model_Cartn_ideal.first)
      s << at.pdbx_model_Cartn_ideal.second.format();
   return s;
} 


std::ostream&
coot::operator<<(std::ostream &s, const dict_bond_restraint_t &rest) {

   s << "[bond-restraint: " 
     << rest.atom_id_1_4c() << " "
     << rest.atom_id_2_4c() << " "
     << rest.type() << " " << std::setw(7) << rest.value_dist() << " " << rest.value_esd()
     <<  "]";
   return s;
}

std::ostream& coot::operator<<(std::ostream &s, const coot::dict_chem_comp_t &rest) {

   s << "[dict_chem_comp comp_id: \"" << rest.comp_id << "\" 3-letter-code: \""
     << rest.three_letter_code << "\" name: \"" << rest.name << "\" group: \"" << rest.group
     << "\" descr-level: \"" << rest.description_level << "\" "  << rest.number_atoms_all
     << " " << rest.number_atoms_nh << "]";
   return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const dict_chiral_restraint_t &rest) {

   s << "[chiral: " << rest.Chiral_Id() << " "
     << rest.atom_id_c_4c() << " "
     << rest.atom_id_1_4c() << " "
     << rest.atom_id_2_4c() << " "
     << rest.atom_id_3_4c() << " "
     << rest.volume_sign << "]";
   return s;
}


std::ostream& coot::operator<<(std::ostream&s, coot::dict_plane_restraint_t rest) {

   s << "[plane-restraint: " << rest.plane_id << " " << " {"
     << rest.n_atoms() << " atoms} ";
   for (int iatom=0; iatom<rest.n_atoms(); iatom++) {
      s << ":" << rest[iatom].first << " " << rest[iatom].second << ": ";
   }
   s << "]";
   return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const dict_angle_restraint_t &rest) {

   s << "[angle-restraint: " 
     << rest.atom_id_1_4c() << " "
     << rest.atom_id_2_4c() << " "
     << rest.atom_id_3_4c() << " "
     << rest.angle_ << " " << rest.angle_esd_
     <<  "]";
   return s;
} 

std::ostream&
coot::operator<<(std::ostream &s, const coot::dict_torsion_restraint_t &rest) {
   s << "[torsion-restraint: " << rest.id() << " "
     << "\"" << rest.atom_id_1_4c() << "\"" << " "
     << "\"" << rest.atom_id_2_4c() << "\"" << " "
     << "\"" << rest.atom_id_3_4c() << "\"" << " "
     << "\"" << rest.atom_id_4_4c() << "\"" << " "
     << rest.angle() << " " 
     << rest.esd() << " " 
     << rest.periodicity();
   if (rest.is_const())
      s << " CONST ";
   s << "]";
   return s;
}

// hack for mac, ostream problems
std::string
coot::dict_torsion_restraint_t::format() const {

   std::string s = "[torsion-restraint: ";
   s +=  id();
   s += " ";
   s +=  atom_id_1_4c();
   s +=  " ";
   s +=  atom_id_2_4c();
   s +=   " ";
   s +=  atom_id_3_4c();
   s +=   " ";
   s +=  atom_id_4_4c();
   s +=   " ";
   s +=  coot::util::float_to_string(angle());
   s +=  " " ;
   s +=  coot::util::float_to_string(esd());
   s +=  " ";
   s +=  coot::util::int_to_string(periodicity());
   if (is_const())
      s +=  " CONST ";
   s +=  "]";
   return s;
}


std::ostream &
coot::operator<<(std::ostream &s, const coot::dictionary_residue_restraints_t &rest) {

   std::cout << "--- dict " << rest.residue_info.comp_id << std::endl;
   std::cout << "    " << rest.atom_info.size() << " atoms" << std::endl;
   for (unsigned int iat=0; iat<rest.atom_info.size(); iat++)
      std::cout << "   " << rest.atom_info[iat] << std::endl;
   std::cout << "    " << rest.bond_restraint.size() << " bonds" << std::endl;
   for (unsigned int ibond=0; ibond<rest.bond_restraint.size(); ibond++)
      std::cout << "   " << rest.bond_restraint[ibond] << std::endl;

   return s;
}





std::string
coot::protein_geometry::get_padded_name(const std::string &atom_id,
					const int &comp_id_index) const {
   std::string s;
   if (comp_id_index < 0) {
      std::cout << "ERROR:: disaster in get_padded_name for comp_id_index "
		<< comp_id_index << " and atom name \"" << atom_id << "\"" << std::endl;
      return s;
   } else {
      for (unsigned int iat=0; iat<dict_res_restraints[comp_id_index].second.atom_info.size(); iat++) {
	 if (dict_res_restraints[comp_id_index].second.atom_info[iat].atom_id == atom_id) {
	    s = dict_res_restraints[comp_id_index].second.atom_info[iat].atom_id_4c;
	    break;
	 }
      }
   }
   return s;
}




void
coot::protein_geometry::info() const {


   std::cout << "::::: MONOMER GEOMETRY:" << std::endl;
   for (unsigned int idr=0; idr<size(); idr++) {
      // ejd says that "restraints" has an "n" in it.  Fixed.
      std::cout << dict_res_restraints[idr].second.residue_info.comp_id << std::endl;
      std::cout << "   " << dict_res_restraints[idr].second.bond_restraint.size()
		<< " bond restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].second.angle_restraint.size()
		<< " angle restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].second.torsion_restraint.size()
		<< " torsion restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].second.plane_restraint.size()
		<< " plane restraints " << std::endl;
//       for (int i=0; i<dict_res_restraints[idr].plane_restraint.size(); i++) {
// 	 for (int j=0; j<dict_res_restraints[idr].plane_restraint[i].atom_ids.size(); j++) {
// 	    std::cout << "      " << dict_res_restraints[idr].plane_restraint[i].plane_id
// 		      << " " << dict_res_restraints[idr].plane_restraint[i].atom_ids[j]
// 		      << std::endl;
// 	 }
//       }
   }

   std::cout << "::::: LINK GEOMETRY:" << std::endl;
   for (unsigned int idr=0; idr<dict_link_res_restraints.size(); idr++) {
      std::cout << dict_link_res_restraints[idr].link_id << std::endl;
      std::cout << "   " << dict_link_res_restraints[idr].link_bond_restraint.size()
		<< " link bond restraits " << std::endl;
      std::cout << "   " << dict_link_res_restraints[idr].link_angle_restraint.size()
		<< " link angle restraits " << std::endl;
      std::cout << "   " << dict_link_res_restraints[idr].link_torsion_restraint.size()
		<< " link torsion restraits " << std::endl;
      std::cout << "   " << dict_link_res_restraints[idr].link_plane_restraint.size()
		<< " link plane restraits " << std::endl;
   }
   
} 

// static
unsigned int
coot::chem_link::make_hash_code(const std::string &comp_id_1, const std::string &comp_id_2, const std::string &group_1, const std::string &group_2) {

   unsigned int hash = 0;
   unsigned int hash_c1 = 0;
   unsigned int hash_c2 = 0;
   unsigned int hash_g1 = 0;
   unsigned int hash_g2 = 0;

   std::string local_group_1 = group_1;
   std::string local_group_2 = group_2;

   if (local_group_1 == "L-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "L-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "P-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "P-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "M-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "M-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "D-pyranose")   local_group_1 = "pyranose";
   if (local_group_2 == "D-pyranose")   local_group_2 = "pyranose";
   if (local_group_1 == "D-SACCHARIDE") local_group_1 = "pyranose";  // CCD annotation for MAN, etc
   if (local_group_2 == "D-SACCHARIDE") local_group_2 = "pyranose";
   if (local_group_1 == "SACCHARIDE")   local_group_1 = "pyranose";  // CCD annotation for FUC
   if (local_group_2 == "SACCHARIDE")   local_group_2 = "pyranose";

   if (false)
      std::cout << "debug:: in make_hash_code " << group_1 << " " << group_2 << " -> "
		<< local_group_1 << " " << local_group_2 << std::endl;

   for (unsigned int i = 0; i < comp_id_1.length(); i++) {
     unsigned int chr = comp_id_1[i];
     hash_c1  = ((hash_c1 << 5) - hash_c1) + chr;
     hash_c1 |= 0; // Convert to 32bit integer
   }
   for (unsigned int i = 0; i < comp_id_2.length(); i++) {
     unsigned int chr = comp_id_2[i];
     hash_c2  = ((hash_c2 << 5) - hash_c2) + chr;
     hash_c2 |= 0;
   }
   for (unsigned int i = 0; i < local_group_1.length(); i++) {
     unsigned int chr = local_group_1[i];
     hash_g1  = ((hash_g1 << 5) - hash_g1) + chr;
     hash_g1 |= 0;
   }
   for (unsigned int i = 0; i < local_group_2.length(); i++) {
     unsigned int chr = local_group_2[i];
     hash_g2  = ((hash_g2 << 5) - hash_g2) + chr;
     hash_g2 |= 0;
   }

   // hash = hash_c1 + 2 * hash_c2 + 4 * hash_g1 + hash_g2;

   hash = hash_g1 + 8 * hash_g2;

   // hash |= 0;

   return hash;
}

std::pair<bool, bool>
coot::chem_link::matches_comp_ids_and_groups_hashed(const std::string &comp_id_1,
					     const std::string &group_1,
					     const std::string &comp_id_2,
					     const std::string &group_2) const {

   bool match = false; // initially
   bool order_switch = false;

   std::string local_group_1 = group_1;
   std::string local_group_2 = group_2;

   // chem_links specify "peptide" or "pyranose", but comp_groups are "L-peptide"/"D-pyranose".
   // So allow them to match.
   // 201201013 (Friday) allow M-peptides to match too.
   if (local_group_1 == "L-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "L-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "M-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "M-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "D-pyranose") local_group_1 = "pyranose";
   if (local_group_2 == "D-pyranose") local_group_2 = "pyranose";
   if (local_group_1 == "D-SACCHARIDE") local_group_1 = "pyranose";  // CCD annotation for MAN, etc
   if (local_group_1 == "SACCHARIDE") local_group_1 = "pyranose";    // CCD annotation for FUC
   if (local_group_2 == "D-SACCHARIDE") local_group_2 = "pyranose";
   if (local_group_2 == "SACCHARIDE") local_group_2 = "pyranose";

   unsigned int hash_test = make_hash_code(comp_id_1, comp_id_2, local_group_1, local_group_2);

   if (hash_test == hash_code) {
      match = true;
      order_switch = false;
   }

   return std::pair<bool, bool>(match, order_switch);
}


std::pair<bool, bool>
coot::chem_link::matches_comp_ids_and_groups_hashed(unsigned int hash_test_forwards,
                                                    unsigned int hash_test_backwards) const {

   bool match = false; // initially
   bool order_switch = false;

   if (hash_test_forwards == hash_code)
      match = true;

   if (hash_test_backwards == hash_code) {
      match = true;
      order_switch = true;
   }

   return std::pair<bool, bool>(match, order_switch);
}




std::pair<bool, bool>
coot::chem_link::matches_comp_ids_and_groups(const std::string &comp_id_1,
					     const std::string &group_1,
					     const std::string &comp_id_2,
					     const std::string &group_2) const {

   bool debug = false;

   if (debug) {
      std::cout << "   ------ DEBUG:: in matches_comp_ids_and_groups() "
		<< id << " chem_link_name " << chem_link_name << ": input comp_ids "
		<< comp_id_1 << " and " << comp_id_2 << " vs ref link-comp_id-1 :"
		<< chem_link_comp_id_1 << ": ref link-comp_id-2:"
		<< chem_link_comp_id_2 << ":" << std::endl; 
      std::cout << "         for chem_link_comp_name " << chem_link_name << ": input groups "
		<< group_1 << " and " << group_2 << " vs ref link-group-1 :"
		<< chem_link_group_comp_1 << ": ref link-group-2 :"
		<< chem_link_group_comp_2 << ":" << std::endl;
   }

   unsigned int hash_test = make_hash_code(comp_id_1, comp_id_2, group_1, group_2);

   if (debug)
      std::cout << "   checking:: hash comparison this: " << hash_code << " test: "
		<< hash_test << " " << *this << "\n";

   bool match = false; // initially
   bool order_switch = false;

   std::string local_group_1 = group_1; 
   std::string local_group_2 = group_2;

   // chem_links specify "peptide" or "pyranose", but comp_groups are "L-peptide"/"D-pyranose".
   // So allow them to match.
   // 201201013 (Friday) allow M-peptides to match too.
   // If you are thinkg of adding/changing this, change the one in make_hash_code() too.
   if (local_group_1 == "L-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "L-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "P-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "P-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "M-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "M-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "D-pyranose") local_group_1 = "pyranose";
   if (local_group_2 == "D-pyranose") local_group_2 = "pyranose";
   if (local_group_1 == "D-SACCHARIDE") local_group_1 = "pyranose";  // CCD annotation for MAN, etc
   if (local_group_2 == "D-SACCHARIDE") local_group_2 = "pyranose";
   if (local_group_1 == "SACCHARIDE")   local_group_1 = "pyranose";    // CCD annotation for FUC
   if (local_group_2 == "SACCHARIDE")   local_group_2 = "pyranose";

   if (debug)
      std::cout << "     sigh... check match: dict \""
		<< chem_link_group_comp_1 << "\" vs model \"" << local_group_1 << "\" -and- dict \""
		<< chem_link_group_comp_2 << "\" vs model \"" << local_group_2 << "\"\n";

   if (local_group_2 == "SACCHARIDE") local_group_2 = "pyranose";
   
   if (debug) {
      std::cout << "   ------ DEBUG:: in matches_comp_ids_and_groups() "
		<< id << " chem_link_name " << chem_link_name << ": input comp_ids "
		<< comp_id_1 << " and " << comp_id_2 << " vs ref link-comp_id-1 :"
		<< chem_link_comp_id_1 << ": ref link-comp_id-2:"
		<< chem_link_comp_id_2 << ":" << std::endl;
      std::cout << "         for chem_link_comp_name " << chem_link_name << ": input groups "
		<< local_group_1 << " and " << local_group_2 << " vs ref link-group-1 :"
		<< chem_link_group_comp_1 << ": ref link-group-2 :"
		<< chem_link_group_comp_2 << ":" << std::endl;
   }

   if (((chem_link_group_comp_1 == "") || (chem_link_group_comp_1 == local_group_1)) &&
       ((chem_link_group_comp_2 == "") || (chem_link_group_comp_2 == local_group_2))) {
      if (((chem_link_comp_id_1 == "") || (chem_link_comp_id_1 == comp_id_1)) &&
	  ((chem_link_comp_id_2 == "") || (chem_link_comp_id_2 == comp_id_2))) {
	 match = true;
      }
   }

   if (((chem_link_group_comp_1 == "DNA/RNA") && (local_group_1 == "RNA") && 
	(chem_link_group_comp_1 == "DNA/RNA") && (local_group_2 == "RNA")))
      match = true;

   if (((chem_link_group_comp_1 == "DNA/RNA") && (local_group_1 == "DNA") && 
	(chem_link_group_comp_1 == "DNA/RNA") && (local_group_2 == "DNA")))
      match = true;

   if (match) {
      
      // OK, nothing more to do
      
   } else { 
      
      // And what about if the residues come here backward? We should
      // report a match and that they should be reversed to the calling
      // function?  

      // reverse index 
      if (((chem_link_group_comp_1 == "") || (chem_link_group_comp_1 == local_group_2)) &&
	  ((chem_link_group_comp_2 == "") || (chem_link_group_comp_2 == local_group_1))) { 
	 if (((chem_link_comp_id_1 == "") || (chem_link_comp_id_1 == comp_id_2)) &&
	     ((chem_link_comp_id_2 == "") || (chem_link_comp_id_2 == comp_id_1))) {
	    // std::cout << "debug:: matched with order switch " << std::endl;
	    match = true;
	    order_switch = true;
	 }
      }
   }
   
   if (debug)
      if (match)
	 std::cout << "!!! matches_comp_ids_and_groups() passed comp_id_1: \""
		   << comp_id_1 << "\" and group_1: \"" << group_1
		   << "\" and comp_id_2: \"" 
		   << comp_id_2 << "\" and group_2: \"" << group_2
		   << "\" and this: " << *this
		   << "\" returns " << match << std::endl;
   
   return std::pair<bool, bool>(match, order_switch);
}

std::ostream& coot::operator<<(std::ostream &s, coot::chem_link lnk) {

   s << "[chem_link: id: " << lnk.id
     << " [comp_id1: \"" << lnk.chem_link_comp_id_1 << "\" group_1: \"" << lnk.chem_link_group_comp_1
     << "\" mod_1: \"" << lnk.chem_link_mod_id_1 << "\"] to "
     << " [comp_id2: \"" << lnk.chem_link_comp_id_2 << "\" group_2: \"" << lnk.chem_link_group_comp_2
     << "\" mod_2: \"" << lnk.chem_link_mod_id_2 << "\"] " << lnk.chem_link_name << "]";
   return s; 
}

std::ostream& coot::operator<<(std::ostream &s, coot::list_chem_mod mod) {

   s << "[list_chem_mod: id: " << mod.id << " " 
     << "name: " << mod.name 
     << " comp_id :" << mod.comp_id 
     << ": group_id: " << mod.group_id
     << "]";
   return s;
} 


// return "" on failure.
// no order switch is considered.
// 
std::string
coot::protein_geometry::find_glycosidic_linkage_type(mmdb::Residue *first, mmdb::Residue *second) const {

   // Fixup needed for PDBv3

   bool debug = false;
   double critical_dist = 2.4; // A, less than that and Coot should
			       // try to make the bond.
                               // 20170505: changed to 2.4, was 3.0.
                               // Needed to stop the beta1-6 linked MAN
                               // on a BMA linkinking to the NAG to which
                               // the BMA is bonded (A605-602 in 3whe).
                              
   mmdb::PPAtom res_selection_1 = NULL;
   mmdb::PPAtom res_selection_2 = NULL;
   int i_no_res_atoms_1;
   int i_no_res_atoms_2;
   double d;
   std::vector<coot::glycosidic_distance> close;
   
   first->GetAtomTable( res_selection_1, i_no_res_atoms_1);
   second->GetAtomTable(res_selection_2, i_no_res_atoms_2);

   for (int i1=0; i1<i_no_res_atoms_1; i1++) { 
      clipper::Coord_orth a1(res_selection_1[i1]->x,
			     res_selection_1[i1]->y,
			     res_selection_1[i1]->z);
      for (int i2=0; i2<i_no_res_atoms_2; i2++) {
	 clipper::Coord_orth a2(res_selection_2[i2]->x,
				res_selection_2[i2]->y,
				res_selection_2[i2]->z);
	 d = (a1-a2).lengthsq();
	 if (d < critical_dist*critical_dist) {
	    close.push_back(coot::glycosidic_distance(res_selection_1[i1],
						      res_selection_2[i2],
						      sqrt(d)));
	 }
      }
   }

   std::sort(close.begin(), close.end());

   // if you consider to uncomment this to debug a repulsion instead
   // the forming of a glycosidic bond, consider the residue numbering
   // of the residues involved: the "residue 1" should have the O4 and
   // the "residue 2" (+1 residue number) should have the C1.
   // 
   if (debug) {
      std::cout << "DEBUG:: find_glycosidic_linkage_type() "
		<< "number of sorted distances: "
		<< close.size() << std::endl;
      for (unsigned int i=0; i<close.size(); i++) {
	 std::cout << "#### glyco close: " << close[i].distance << "  "
		   << close[i].at1->GetChainID() << " " 
		   << close[i].at1->GetSeqNum() << " " 
		   << close[i].at1->GetAtomName() << " " 
		   << " to "
		   << close[i].at2->GetChainID() << " " 
		   << close[i].at2->GetSeqNum() << " " 
		   << close[i].at2->GetAtomName() << " " 
		   << std::endl;
      }
   }

   std::string link_type("");

   // glyco_chiral constructor can throw an exception
   try { 
   
      float smallest_link_dist = 99999.9;
      for (unsigned int i=0; i<close.size(); i++) {
	 std::string name_1(close[i].at1->name);
	 std::string name_2(close[i].at2->name);


	 // First test the NAG-ASN link (that order - as per dictionary)
	 //
	 if (name_1 == " C1 ")
	    if (name_2 == " ND2")
	       if (close[i].distance < smallest_link_dist) {
		  smallest_link_dist = close[i].distance;
		  link_type = "NAG-ASN";
	       }
      
	 if (name_1 == " O4 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-4");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-4";
		  }
	       }
      
	 if (name_1 == " O2 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-2");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-2";
		  }
	       }

	 if (name_1 == " O3 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-3");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-3";
		  }
	       }
      

	 // The BETA2-3 link should never happen :-)
	 // There are no biosynthetic pathways to make an BETA2-3 link for a SIA.
	 // (SIA BETA2-3 would be axial if it existed)
	 //
	 if (name_1 == " C2 " )
	    if (name_2 == " O3 ")
	       if (std::string(close[i].at1->GetResName()) == "SIA") {
		  if (close[i].distance < smallest_link_dist) {
		     coot::atom_quad glyco_chiral_quad(first, second, "ALPHA2-3");
		     std::cout << "   glyco_chiral ALPHA2-3 "
			       << close[i].at1->GetResName() << " "
			       << close[i].at2->GetResName() << " "
			       << glyco_chiral_quad.chiral_volume() << std::endl;
		     if (glyco_chiral_quad.chiral_volume() > 0.0) {
			smallest_link_dist = close[i].distance;
			link_type = "ALPHA2-3";
		     }
		  }
	       }

	 // 20180111 Add ALPHA2-6 links for SIA
	 if (name_1 == " C2 " )
	    if (name_2 == " O6 ")
	       if (std::string(close[i].at1->GetResName()) == "SIA") {
		  if (close[i].distance < smallest_link_dist) {
		     coot::atom_quad glyco_chiral_quad(first, second, "ALPHA2-6");
		     std::cout << "   glyco_chiral ALPHA2-6 "
			       << close[i].at1->GetResName() << " "
			       << close[i].at2->GetResName() << " "
			       << glyco_chiral_quad.chiral_volume() << std::endl;
		     if (glyco_chiral_quad.chiral_volume() > 0.0) {
			smallest_link_dist = close[i].distance;
			link_type = "ALPHA2-6";
		     }
		  }
	       }

	 if (name_1 == " O6 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-6");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-6";
		  }
	       }
	       
	 if (name_1 == " O2 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-2");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-2";
		  }
	       }
	       
	 if (name_1 == " O3 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-3");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-3";
		  }
	       }

	 if (name_1 == " C2 " )
	    if (name_2 == " O3 ")
	       if (std::string(close[i].at1->GetResName()) == "SIA") { 
		  if (close[i].distance < smallest_link_dist) {
		     coot::atom_quad glyco_chiral_quad(first, second, "ALPHA2-3");
		     std::cout << "   glyco_chiral ALPHA2-3 "
			       << close[i].at1->GetResName() << " "
			       << close[i].at2->GetResName() << " "
			       << glyco_chiral_quad.chiral_volume() << std::endl;
		     if (glyco_chiral_quad.chiral_volume() < 0.0) { 
			smallest_link_dist = close[i].distance;
			link_type = "ALPHA2-3";
		     }
		  }
	       }
      
	 if (name_1 == " O4 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-4");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-4";
		  }
	       }
      
	 if (name_1 == " O6 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-6");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-6";
		  }
	       }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING::" << rte.what() << std::endl;
   }

   if (false) 
      std::cout << "   debug:: find_glycosidic_linkage_type() for "
		<< first->GetChainID() << " " << first->GetSeqNum() << " " << first->GetInsCode()
		<< first->GetResName() << ","
		<< second->GetChainID() << " " << second->GetSeqNum() << " " << second->GetInsCode()
		<< second->GetResName() 
		<< " returns \"" << link_type << "\""
		<< std::endl;
   
   return link_type;
}

std::string
coot::protein_geometry::find_glycosidic_linkage_type(mmdb::Residue *first,
						     mmdb::Residue *second,
						     mmdb::Manager *mol) const {
   
   // First check that the residues are LINK - or are sequential.
   // If so, then we can find the link as above.
   
   std::string link_type;
   bool are_linked = false;
   bool are_sequential = false;

   // Test for sequential/tandem
   //
   std::string chain_id_1 =  first->GetChainID();
   std::string chain_id_2 = second->GetChainID();
   int resno_1 =  first->GetSeqNum();
   int resno_2 = second->GetSeqNum();
   if (chain_id_1 == chain_id_2) {
      if (resno_1 == (resno_2+1)) are_sequential = true;
      if (resno_2 == (resno_1+1)) are_sequential = true;
   } 

   if (! are_sequential) {
      // Test for link in molecule
      std::string ins_code_1 =  first->GetInsCode();
      std::string ins_code_2 = second->GetInsCode();
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) { 
	 mmdb::LinkContainer *links = model_p->GetLinks();
	 int n_links = model_p->GetNumberOfLinks();
	 for (int ilink=1; ilink<=n_links; ilink++) { 
	    mmdb::Link *link = model_p->GetLink(ilink);
	    if (link) {
	       are_linked = are_linked_in_order(first, second, link);
	       if (! are_linked)
		  are_linked = are_linked_in_order(first, second, link);
	       if (are_linked)
		  break;
	    }
	 }
      }
   } 


   if (are_linked || are_sequential)
      link_type = find_glycosidic_linkage_type(first, second);

   return link_type;
}


bool
coot::protein_geometry::are_linked_in_order(mmdb::Residue *first,
					    mmdb::Residue *second,
					    mmdb::Link *link) const {

   bool linked = false;
   // c.f link_atoms() - but we can't use that because that's a
   // coot-utils function and geometry is more primitive.
   std::string link_chain_id_1(link->chainID1);
   std::string link_chain_id_2(link->chainID2);
   std::string chain_id_1 =  first->GetChainID();
   std::string chain_id_2 = second->GetChainID();
   int resno_1 =  first->GetSeqNum();
   int resno_2 = second->GetSeqNum();
   if (link_chain_id_1 == chain_id_1) { 
      if (link_chain_id_2 == chain_id_2) {
	 int link_reso_1 = link->seqNum1;
	 int link_reso_2 = link->seqNum2;
	 if (link_reso_1 == resno_1) { 
	    if (link_reso_2 == resno_2) {
	       std::string link_ins_code_1 = link->insCode1;
	       std::string link_ins_code_2 = link->insCode2;
	       std::string ins_code_1 =  first->GetInsCode();
	       std::string ins_code_2 = second->GetInsCode();
	       if (link_ins_code_1 == ins_code_1) { 
		  if (link_ins_code_2 == ins_code_2) {
		     linked = true;
		  }
	       }
	    }
	 }
      }
   }

   return linked;
} 
					    

std::pair<std::string, bool>
coot::protein_geometry::find_glycosidic_linkage_type_with_order_switch(mmdb::Residue *first, mmdb::Residue *second) const {

   std::pair<std::string, bool> r("", false);

   std::string l = find_glycosidic_linkage_type(first, second);

   if (l == "") { 
      l = find_glycosidic_linkage_type(second,first);
      if (l != "") {
	 r.first = l;
	 r.second = true;
      } 
   } else {
      r.first = l;
      r.second = false;
   } 
   return r;
} 



      
void
coot::protein_geometry::set_verbose(bool verbose_flag) {
   verbose_mode = verbose_flag;
}


// return empty file name on failure.
std::string
coot::protein_geometry::comp_id_to_file_name(const std::string &comp_id) const {

   std::string file_name;

   if (comp_id.length() > 0) { 
      char *cmld = getenv("COOT_MONOMER_LIB_DIR");
      std::string d;
      if (! cmld) {
	 d = PKGDATADIR;
	 d = util::append_dir_dir(d, "lib");
	 d = util::append_dir_dir(d, "data");
	 d = util::append_dir_dir(d, "monomers");
      } else {
	 d = cmld;
      }

      if (! d.empty()) {
	 std::string c0(1,comp_id[0]);
	 std::string lc0 = util::downcase(c0);
	 d = util::append_dir_dir(d, lc0);
	 std::string lfn = comp_id + ".cif";
	 file_name = util::append_dir_file(d, lfn);
      }
   }
   return file_name;
} 


// Return 0 on failure to do a dynamic add (actually, the number of
// atoms read).
//
// try_dynamic_add() will add with an imol of IMOL_ENC_ANY
// 
int
coot::protein_geometry::try_dynamic_add(const std::string &resname, int read_number) {

   int success = 0;  // fail initially and is set to the number of
		     // atoms read from the mmcif dictionary file in
		     // init_refmac_mon_lib().

   // If this is INH, DRG etc, don't try to auto-add
   // 
   if (is_non_auto_load_ligand(resname)) {
      std::cout << "INFO:: comp-id: " << resname
		<< " is marked for non-autoloading - stopping now "
		<< std::endl;
      return success;
   }

   // So what is happening here?
   //
   // It is a heirachy of setting
   //
   // The highest priority is COOT_REFMAC_LIB_DIR, if that is set we use it.
   // If COOT_REFMAC_LIB_DIR is then try CLIB (a CCP4 setting).
   // If that is not set, then we fall back to the default directory:
   // $prefix/share/coot onto which we tag a "lib" dir.
   // 

   char *s  = getenv("COOT_REFMAC_LIB_DIR");
   char *cmld = getenv("COOT_MONOMER_LIB_DIR");

   if (s) {

   } else {
      s  = getenv("CLIB");

      if (! s) {
	 std::cout << "DEBUG:: try_dynamic_add() using package_data_dir(): " << package_data_dir()
		   << std::endl;
	 std::string tmp_string = package_data_dir();
	 tmp_string = util::append_dir_dir(tmp_string, "lib");
	 s = new char[tmp_string.length() + 1];
	 strcpy(s, tmp_string.c_str());
      } else {
	 if (verbose_mode)
	    std::cout << "INFO:: using standard CCP4 Refmac dictionary"
		      << " to search for \"" << resname << "\"" << std::endl; 
      }
   }

   if (!s) {

      std::cout << "debug:: try_dynamic_add() bad news - null s" << std::endl;

   } else {
      std::string filename(s);
      std::string beta_anomer_name;
      std::string alpha_anomer_name;
      std::string alt_beta_anomer_name;
      std::string alt_alpha_anomer_name;

      if (cmld) { 
	 filename = cmld;
	 filename += "/";
      } else {
	 filename += "/data/monomers/";
      }

      filename = coot::util::intelligent_debackslash(filename);
      if (resname.length() > 0) {
	 const char rs = resname[0];
	 const char v = tolower(rs); // get the sub directory name
         char v1[2];
	 v1[0] = v;
	 v1[1] = '\0';
	 std::string letter(v1);
	 filename += letter;
	 filename += "/";
	 std::string upcased_resname_filename = filename;
	 if (resname.length() > 2 ) { 
	    if (resname[2] != ' ') { 
	       filename += resname;
	       upcased_resname_filename += coot::util::upcase(resname);
	    } else {
	       filename += resname.substr(0,2);
	       upcased_resname_filename += coot::util::upcase(resname.substr(0,2));
	    }
	 } else {
	    filename += resname;
	    upcased_resname_filename += coot::util::upcase(resname);	    
	 }
	 beta_anomer_name = filename;
	 beta_anomer_name += "-b-D.cif";
	 alt_beta_anomer_name = filename;
	 alt_beta_anomer_name += "-B-D.cif";
	 alpha_anomer_name  = filename;
	 alpha_anomer_name += "-a-L.cif";
	 alt_alpha_anomer_name  = filename;
	 alt_alpha_anomer_name += "-A-L.cif";
	 filename += ".cif";
	 upcased_resname_filename += ".cif";
	 
	 struct stat buf;
	 int istat = stat(filename.c_str(), &buf);
	 if (istat == 0) {
	    if (coot::is_regular_file(filename)) {

	       coot::read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(filename, read_number);
	       success = rmit.success;
	    } else {

	       // error/warning
	       if (! coot::is_dir_or_link(filename)) {
		  std::cout << "WARNING: " << filename << ": no such file (or directory)\n";
	       } else {
		  std::cout << "ERROR: dictionary " << filename << " is not a regular file" << std::endl;
	       }
	    }
	 } else { 
	    
	    // Regular file doesn't exist,
	    
	    // Try the upcased filename
	    // 
	    istat = stat(upcased_resname_filename.c_str(), &buf);
	    if (istat == 0) {
	       read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(upcased_resname_filename, read_number);
	       success = rmit.success;
	    } else { 

	       // try the beta anomer version
	       istat = stat(beta_anomer_name.c_str(), &buf);
	       if (istat == 0) {
		  read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(beta_anomer_name, read_number);
		  success = rmit.success;
	       } else {
		  // try the upcased file name e.g. xxx/NAG-B-D.cif
		  istat = stat(alt_beta_anomer_name.c_str(), &buf);
		  if (istat == 0) {
		     read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(alt_beta_anomer_name, read_number);
		     success = rmit.success;
		  } else {
		     // alpha?
		     istat = stat(alpha_anomer_name.c_str(), &buf);
		     if (istat == 0) {
			read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(alpha_anomer_name, read_number);
			success = rmit.success;
		     } else {
			istat = stat(alt_alpha_anomer_name.c_str(), &buf);
			if (istat == 0) {
			   read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(alt_alpha_anomer_name, read_number);
			   success = rmit.success;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   // now, did we get something with minimal description? If so,
   // delete it, it was a fail.
   // 
   std::pair<bool, dictionary_residue_restraints_t> p = get_monomer_restraints(resname, IMOL_ENC_ANY);
   if (! p.first) {
      success = 0;
   } else {
      // elide minimal description restraints.
      if (p.second.residue_info.description_level == "M") { 
	 success = 0;
	 delete_mon_lib(resname, IMOL_ENC_ANY);
      }
   }
   return success;
}


bool
coot::protein_geometry::is_non_auto_load_ligand(const std::string resname) const {

   bool r = false;
   std::vector<std::string>::const_iterator it;
   for (it=non_auto_load_residue_names.begin(); it!=non_auto_load_residue_names.end(); it++) {
      if (*it == resname) {
	 r = true;
	 break;
      } 
   }
   return r;
}

void
coot::protein_geometry::add_non_auto_load_residue_name(const std::string &res_name) {

   bool found = false;
   std::vector<std::string>::const_iterator it;
   for (it=non_auto_load_residue_names.begin(); it!=non_auto_load_residue_names.end(); it++) {
      if (*it == res_name) {
	 found = true;
	 break;
      }
      if (found)
	 break;
   }
   if (! found)
      non_auto_load_residue_names.push_back(res_name);
}

void
coot::protein_geometry::remove_non_auto_load_residue_name(const std::string &res_name) {

   std::vector<std::string>::iterator it;
   for (it=non_auto_load_residue_names.begin(); it!=non_auto_load_residue_names.end(); it++) {
      if (*it == res_name) {
	 non_auto_load_residue_names.erase(it);
	 break;
      }
   }
} 



void
coot::protein_geometry::fill_default_non_auto_load_residue_names() { // build-it default
   non_auto_load_residue_names.push_back("XXX");
   non_auto_load_residue_names.push_back("LIG");
   non_auto_load_residue_names.push_back("DRG");
   non_auto_load_residue_names.push_back("INH");
   non_auto_load_residue_names.push_back("LG0");
   non_auto_load_residue_names.push_back("LG1");
   non_auto_load_residue_names.push_back("LG2");
   non_auto_load_residue_names.push_back("LG3");
   non_auto_load_residue_names.push_back("LG4");
   non_auto_load_residue_names.push_back("LG5");
   non_auto_load_residue_names.push_back("LG6");
   non_auto_load_residue_names.push_back("LG7");
   non_auto_load_residue_names.push_back("LG8");
   non_auto_load_residue_names.push_back("LG9");
}


std::vector <coot::dict_torsion_restraint_t>
coot::protein_geometry::get_monomer_torsions_from_geometry(const std::string &monomer_type,
							   int imol_enc) {

   bool ifound = false;
   std::vector <coot::dict_torsion_restraint_t> rv;

   int idr = get_monomer_restraints_index(monomer_type, imol_enc, false);

   if (idr >= 0 ) {
      ifound = true;
      return dict_res_restraints[idr].second.torsion_restraint;
   }

   // check synonyms before 3-letter codes.  Maybe needs matches_imol() test
   // 
   if (! ifound) { 
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) { 
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) { // is this the right test?
	       int ndict = dict_res_restraints.size();
	       for (int i=0; i<ndict; i++) {
		  if (dict_res_restraints[i].second.residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		     ifound = 1;
		     rv = dict_res_restraints[i].second.torsion_restraint;
		     break;
		  }
	       }
	    }
	 }
	 if (ifound)
	    break;
      }
   }

   if (ifound == 0) {
      int read_number = 40; // just a filler, FIXME.
      int nbonds = try_dynamic_add(monomer_type, read_number);
      if (nbonds > 0) {
	 // OK, we got it, what are the torsions?
	 for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
	    if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	       ifound = 1;
	       rv = dict_res_restraints[i].second.torsion_restraint;
	       break;
	    }
	 }
      }
   }

   if (ifound == 0) { // which it should be if we get here with that
		      // return in the loop
      std::cout << "WARNING: residue type " << monomer_type << " not found "
		<< "in restraints dictionary (torsion)" << std::endl;
   }
   rv = filter_torsion_restraints(rv);
   return rv;
}  

std::vector <coot::dict_torsion_restraint_t>
coot::protein_geometry::get_monomer_torsions_from_geometry(const std::string &monomer_type,
							   int imol_enc,
							   bool find_hydrogen_torsions_flag) const {

   bool ifound = 0;
   int ii = -1; // unset, set when torsion found, ii is the index in
		// dict_res_restraints for the torsion restraints,
		// needed so that we can ask about hydrogens.
   std::vector <coot::dict_torsion_restraint_t> unfiltered_torsion_restraints;
   std::vector <coot::dict_torsion_restraint_t> filtered_torsion_restraints;

   int idr = get_monomer_restraints_index(monomer_type, imol_enc, false);

   if (idr >= 0) {
      ifound = 1;
      unfiltered_torsion_restraints = dict_res_restraints[idr].second.torsion_restraint;
      ii = idr; // used for filtering out hydrogens
   }

   if (! ifound) { // which it should be if we get here with that
		      // return in the loop
      // try dynamic add?
      std::cout << "WARNING: residue type " << monomer_type << " not found "
		<< "in restraints dictionary (in get_monomer_torsions_from_geometry(mon, hy)" << std::endl;
   } else {
      if (find_hydrogen_torsions_flag) {
	 filtered_torsion_restraints = unfiltered_torsion_restraints;
      } else {
	 // we don't want torsions that move Hydrogens
	 int nt = dict_res_restraints[ii].second.torsion_restraint.size();
	 for (int it=0; it<nt; it++) { 
	    if (!dict_res_restraints[ii].second.is_hydrogen(dict_res_restraints[ii].second.torsion_restraint[it].atom_id_1())) {
	       if (!dict_res_restraints[ii].second.is_hydrogen(dict_res_restraints[ii].second.torsion_restraint[it].atom_id_4())) {
		  filtered_torsion_restraints.push_back(dict_res_restraints[ii].second.torsion_restraint[it]);
	       }
	    }
	 }
      }
   }

   // more filtering (only one version of a torsion_restraint that
   // have the same middle atoms).
   // 
   filtered_torsion_restraints = filter_torsion_restraints(filtered_torsion_restraints);
   return filtered_torsion_restraints;
}

// Return only one version of a torsions restraint that have the same
// middle atoms.
// 
std::vector <coot::dict_torsion_restraint_t>
coot::protein_geometry::filter_torsion_restraints(const std::vector <coot::dict_torsion_restraint_t> &restraints_in) const {

   std::vector <coot::dict_torsion_restraint_t> r;

   for (unsigned int i=0; i<restraints_in.size(); i++) {
      std::string a2 = restraints_in[i].atom_id_2_4c();
      std::string a3 = restraints_in[i].atom_id_3_4c();
      bool match = false;
      for (unsigned int j=0; j<r.size(); j++) {
	 if (r[j].atom_id_2_4c() == a2)
	    if (r[j].atom_id_3_4c() == a3)
	       match = true;
	 if (r[j].atom_id_2_4c() == a3)
	    if (r[j].atom_id_3_4c() == a2)
	       match = true;
      }
      if (! match) {
	 const dict_torsion_restraint_t &tr = restraints_in[i];
	 r.push_back(tr);
      }
   }


   // why do we crash (on the mac) in the following sort?
   //
   // for (unsigned int i=0; i<restraints_in.size(); i++)
   //    std::cout << " filter_torsion_restraints(): " << restraints_in[i] << std::endl;
   // 20180227 still crashing.
   // Because it is running the sort even if r is of length 0? No, not that.
   // std::cout << " r.size() " << r.size() << std::endl;
   try {
      std::sort(r.begin(), r.end(), torsion_restraints_comparer);
   }
   catch (const std::bad_alloc &e) {
      std::cout << "ERROR caught bad alloc when sorting in filter_torsion_restraints()" << std::endl;
   }
   return r;
}

// static
bool
coot::protein_geometry::torsion_restraints_comparer(const coot::dict_torsion_restraint_t &a, const coot::dict_torsion_restraint_t &b) {

//    const std::string &a2 = a.atom_id_2_4c();
//    const std::string &a3 = a.atom_id_3_4c();
//    const std::string &b2 = b.atom_id_2_4c();
//    const std::string &b3 = b.atom_id_3_4c();

//    const std::string &a2 = a.atom_id_2();
//    const std::string &a3 = a.atom_id_3();
//    const std::string &b2 = b.atom_id_2();
//    const std::string &b3 = b.atom_id_3();

   // we leave these (debugging/test) copies in for now (was testing if reference was the problem
   // (it wasn't).
   //
   std::string a2 = a.atom_id_2();
   std::string a3 = a.atom_id_3();
   std::string b2 = b.atom_id_2();
   std::string b3 = b.atom_id_3();

   if (a2 < b2) {
      return false;
   } else {
      if (a2 > b2) {
 	 return true;
      } else {
	 // a2 and b2 are equal

	 if (a3 <= b3) { // 20180227-PE changed to <=, stops crash

	    // std::cout << "a3 " << a3 << std::endl;
	    // std::cout << "b3 " << b3 << std::endl;

	    return false;

	 }
      }
   }

   //       else
   // 	 if (a3 < b3)
   // 	    return false;

   return true;
}


std::vector <coot::dict_chiral_restraint_t>
coot::protein_geometry::get_monomer_chiral_volumes(const std::string monomer_type,
						   int imol_enc) const { 

   bool ifound = 0;
   std::vector <coot::dict_chiral_restraint_t> rv;
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].first == imol_enc) {
	 if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	    ifound = true;
	    rv = dict_res_restraints[i].second.chiral_restraint;
	    break;
	 }
      }
   }
   if (! ifound) {
      for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
	 if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	    if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	       ifound = true;
	       rv = dict_res_restraints[i].second.chiral_restraint;
	       break;
	    }
	 }
      }
   }
   // OK so the monomer_type did not match the comp_id.  Perhaps the
   // comp_id was not the same as the three letter code, so let's
   // check the monomer_type against the three_letter_codes.
   // 
   if (ifound == 0) {
      for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
	 if (dict_res_restraints[i].second.residue_info.three_letter_code == monomer_type) {
	    ifound = 1;
	    rv = dict_res_restraints[i].second.chiral_restraint;
	    break;
	 }
      }
   } 
   if (ifound == 0) {
      // try dynamic add?
      std::cout << "WARNING: residue type " << monomer_type << " not found "
		<< "in restraints dictionary (chiral)" << std::endl;
   }
   return rv;
}

std::pair<bool, coot::dict_atom>
coot::protein_geometry::get_monomer_atom_info(const std::string &monomer_name,
					      const std::string &atom_name,
					      int imol_enc) const {

   bool status = false;
   dict_atom da;

   std::pair<bool, dictionary_residue_restraints_t> dr = get_monomer_restraints(monomer_name, imol_enc);
   if (dr.first) {
      const std::vector<dict_atom> &atom_info = dr.second.atom_info;
      for (std::size_t i=0; i<atom_info.size(); i++) {
	 dict_atom a = atom_info[i];
	 // std::cout << "comparing atom names \"" << atom_name << "\" and \""
	 // << a.atom_id_4c << "\"" << std::endl;
	 if (atom_name == a.atom_id_4c) {
	    da = a;
	    status = true;
	 }
      }
   }

   return std::pair<bool, dict_atom> (status, da);
}




std::vector <std::string>
coot::protein_geometry::standard_protein_monomer_files() const {

   std::vector <std::string> s;

   s.push_back("a/ALA.cif");
   s.push_back("a/ASP.cif");
   s.push_back("a/ASN.cif");
   s.push_back("c/CYS.cif");
   s.push_back("g/GLN.cif");
   s.push_back("g/GLY.cif");
   s.push_back("g/GLU.cif");
   s.push_back("p/PHE.cif");
   s.push_back("h/HIS.cif");
   s.push_back("i/ILE.cif");
   s.push_back("l/LYS.cif");
   s.push_back("l/LEU.cif");
   s.push_back("m/MET.cif");
   s.push_back("m/MSE.cif");
   s.push_back("p/PRO.cif");
   s.push_back("a/ARG.cif");
   s.push_back("s/SER.cif");
   s.push_back("t/THR.cif");
   s.push_back("v/VAL.cif");
   s.push_back("t/TRP.cif");
   s.push_back("t/TYR.cif");

   s.push_back("p/PO4.cif");
   s.push_back("s/SO4.cif");
   s.push_back("g/GOL.cif");
   s.push_back("c/CIT.cif");
   s.push_back("e/EDO.cif");
   // s.push_back("e/ETH.cif");

   // s.push_back("a/AR.cif"); // argon
   // s.push_back("c/CR.cif"); // Chromium
   // s.push_back("g/GR.cif"); old

   // s.push_back("a/AD.cif");
   // s.push_back("c/CD.cif");
   // s.push_back("g/GD.cif");
   // s.push_back("t/TD.cif");
   // s.push_back("u/UR.cif");

   // new-style (CCP4 6.2) RNA names
   s.push_back("a/A.cif");
   s.push_back("c/C.cif");
   s.push_back("g/G.cif");
   s.push_back("u/U.cif");

   // new-style (CCP4 6.2) DNA names
   s.push_back("d/DA.cif");
   s.push_back("d/DC.cif");
   s.push_back("d/DG.cif");
   s.push_back("d/DT.cif");
   
   s.push_back("h/HOH.cif");
   s.push_back("n/NA.cif"); // for quiet

   return s;
}


// optional arg: bool try_autoload_if_needed=true
bool
coot::protein_geometry::have_dictionary_for_residue_type(const std::string &monomer_type,
							 int imol_enc,
							 int read_number_in,
							 bool try_autoload_if_needed) {

   bool ifound = 0;
   int ndict = dict_res_restraints.size();
   read_number = read_number_in;

   // ---------------- FIXME ----------------------------------------------------
   // ---------------- FIXME ----------------------------------------------------
   // ---------------- FIXME ----------------------------------------------------
   int idr = get_monomer_restraints_index(monomer_type, imol_enc, true);
   if (idr >= 0) {
      ifound = true;
   }

   // check synonyms before checking three-letter-codes
   
   if (! ifound) {
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) { 
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    for (int j=0; j<ndict; j++) {
	       if (dict_res_restraints[j].second.residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		  ifound = true;
		  break;
	       }
	    }
	 }
	 if (ifound)
	    break;
      }
   }
	 
   // OK so the monomer_type did not match the comp_id.  Perhaps the
   // comp_id was not the same as the three letter code, so let's
   // check the monomer_type against the three_letter_codes.
   // 
   if (ifound == 0) {
      for (int i=0; i<ndict; i++) {
	 if (dict_res_restraints[i].second.residue_info.three_letter_code == monomer_type) {
	    if (! dict_res_restraints[i].second.is_bond_order_data_only()) {
	       ifound = 1;
	       break;
	    }
	 }
      }
   }

   if (! ifound) {
      if (try_autoload_if_needed) {
	 ifound = try_dynamic_add(monomer_type, read_number);
	 // std::cout << "here in have_dictionary_for_residue_type() try_dynamic_add returned "
	 // << ifound << std::endl;
      }
   }

   return ifound;
}

// this is const because there is no dynamic add
// This is a test for "shall we add a new dictionary?"
bool
coot::protein_geometry::have_dictionary_for_residue_type_no_dynamic_add(const std::string &monomer_type) const {

   bool ifound = 0;

   int ndict = dict_res_restraints.size();
   for (int i=0; i<ndict; i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 if (matches_imol(dict_res_restraints[i].first, IMOL_ENC_ANY)) {
	    if (! dict_res_restraints[i].second.is_bond_order_data_only()) {
	       ifound = 1;
	       break;
	    }
	 }
      }
   }
   return ifound;
} 


bool
coot::protein_geometry::have_dictionary_for_residue_types(const std::vector<std::string> &residue_types,
							  int imol_enc,
							  int read_number) {

   bool have_all = 1;
   for (unsigned int i=0; i<residue_types.size(); i++) {
      int ifound = have_dictionary_for_residue_type(residue_types[i], imol_enc, read_number);
      if (ifound == 0) {
	 have_all = 0;
      } 
      read_number++;
   }
   return have_all;
}

// this is const because there is no dynamic add
//
// maybe imol should be passed.
bool
coot::protein_geometry::have_at_least_minimal_dictionary_for_residue_type(const std::string &monomer_type,
									  int imol) const {

   //std::cout << "debug:: in have_at_least_minimal_dictionary_for_residue_type() " << monomer_type
   // << " " << imol << std::endl;

   bool ifound = false;
   std::size_t ndict = dict_res_restraints.size();

   for (std::size_t i=0; i<ndict; i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 if (matches_imol(dict_res_restraints[i].first, IMOL_ENC_ANY)) {
	    ifound = true;
	    break;
	 }
	 if (matches_imol(dict_res_restraints[i].first, imol)) {
	    ifound = true;
	    break;
	 }
      }
   }
   return ifound;
} 



// Check that the atom names in the residue match the atom names in
// the dictionary.  Atom " OXT" is treated as a special case (it does
// not cause a failure when " OXT" is not in the dictionary for the
// residue).
//
// Similarly, we can not return problematic status (0) if the
// hydrogens do not match, if caller wishes.
//
// This does only monomer by monomer testing.
//
// There is no test of DEL-O1 (for example) atoms when making links
// between monomers.
// 
// Return in pair.first the state of the match and in .second, the
// list of atoms that do not match the dictionary.
// 
std::pair<bool, std::vector<std::string> >
coot::protein_geometry::atoms_match_dictionary(mmdb::Residue *residue_p,
					       bool check_hydrogens_too_flag,
					       bool apply_bond_distance_check,
					       const coot::dictionary_residue_restraints_t &restraints) const {

   std::vector<std::string> atom_name_vec;
   bool status = 1; // nothing fails to match (so far).

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   bool debug = false;
   if (debug) {
      std::cout << "=== atoms_match_dictionary() with these residue atom names ======= " << std::endl;
      for (int i=0; i<n_residue_atoms; i++) {
	 std::cout << i << "  :" << residue_atoms[i]->name << ": and ele "
		   << residue_atoms[i]->element << std::endl;
      } 
      std::cout << "=== atoms_match_dictionary() with these residue atom names ======= " << std::endl;
      for (unsigned int irat=0; irat<restraints.atom_info.size(); irat++) {
	 std::cout << irat << "  :" << restraints.atom_info[irat].atom_id_4c
		   << ":" << std::endl;
      } 
   }
   
   for (int i=0; i<n_residue_atoms; i++) {

      if (! residue_atoms[i]->isTer()) { 
	 std::string residue_atom_name(residue_atoms[i]->name);
	 std::string ele(residue_atoms[i]->element);

	 bool found = 0;
	 // PDBv3 FIXME
	 if (ele == " H" || ele == " D")
	    if (check_hydrogens_too_flag == false)
	       found = true;

	 if (! found) {
	    for (unsigned int irestraint_atom_name=0; irestraint_atom_name<restraints.atom_info.size(); irestraint_atom_name++) {
	       if (restraints.atom_info[irestraint_atom_name].atom_id_4c == residue_atom_name) {
		  found = 1;
		  break;
	       }
	    }
	 }
	 if (! found) {
	    if (residue_atom_name != " OXT") { 
	       atom_name_vec.push_back(residue_atom_name);
	       status = 0;
	    }
	 }
      }
   }

   // We can finally fail to match because we have some very long
   // bonds (but no atom name mismatches, of course)
   // 
   if (status && apply_bond_distance_check)
      status = atoms_match_dictionary_bond_distance_check(residue_p, check_hydrogens_too_flag, restraints);

   return std::pair<bool, std::vector<std::string> > (status, atom_name_vec);
}


std::pair<bool, std::vector<std::string> >
coot::protein_geometry::atoms_match_dictionary(int imol,
					       mmdb::Residue *residue_p,
					       bool check_hydrogens_too_flag,
					       bool apply_bond_distance_check) const {

   std::string res_name(residue_p->GetResName());
   std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
      get_monomer_restraints(res_name, imol);

   if (restraints.first) {
      return atoms_match_dictionary(residue_p, check_hydrogens_too_flag,
				    apply_bond_distance_check, restraints.second);
   } else { 
      std::vector<std::string> atom_name_vec;
      return std::pair<bool, std::vector<std::string> > (false, atom_name_vec);
   }
}


// return a pair, overall status, and vector of pairs of residue names and
// atom names that dont't match.
//
std::pair<bool, std::vector<std::pair<std::string, std::vector<std::string> > > >
coot::protein_geometry::atoms_match_dictionary(int imol,
					       const std::vector<mmdb::Residue *> &residues,
					       bool check_hydrogens_too_flag,
					       bool apply_bond_distance_check) const {

   bool status = true;
   std::vector<std::pair<std::string, std::vector<std::string> > > p;
   
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      std::string res_name(residues[ires]->GetResName());
      std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
	 get_monomer_restraints(res_name, imol);
      if (restraints.first) {
	 std::pair<bool, std::vector<std::string> > r =
	    atoms_match_dictionary(residues[ires],
				   check_hydrogens_too_flag,
				   apply_bond_distance_check,
				   restraints.second);
	 if (r.first == 0) {
	    std::pair<std::string, std::vector<std::string> > p_bad(res_name, r.second);
	    p.push_back(p_bad);
	    status = 0;
	 }
      } else {
	 std::cout << "ERROR:: atoms_match_dictionary() --- no restraints" << std::endl;
      }
   }
   return std::pair<bool, std::vector<std::pair<std::string, std::vector<std::string> > > > (status, p);
}

bool
coot::protein_geometry::atoms_match_dictionary_bond_distance_check(mmdb::Residue *residue_p,
								   bool check_hydrogens_too_flag,
								   const coot::dictionary_residue_restraints_t &restraints) const { 

   bool status = true; // good
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   if (n_residue_atoms > 2) { 
      for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
	 for (int iat=0; iat<(n_residue_atoms-1); iat++) {
	    const mmdb::Atom *at_1 = residue_atoms[iat];
	    std::string atom_name_1(at_1->name);
	    if (restraints.bond_restraint[ibond].atom_id_1_4c() == atom_name_1) { 
	       for (int jat=iat+1; jat<n_residue_atoms; jat++) {
		  const mmdb::Atom *at_2 = residue_atoms[jat];
		  std::string atom_name_2(at_2->name);
		  if (restraints.bond_restraint[ibond].atom_id_2_4c() == atom_name_2) {
		     std::string alt_conf_1(at_1->altLoc);
		     std::string alt_conf_2(at_2->altLoc);
		     if (alt_conf_1 == alt_conf_2) {
			clipper::Coord_orth pt1(at_1->x, at_1->y, at_1->z);
			clipper::Coord_orth pt2(at_2->x, at_2->y, at_2->z);
			double d = (pt1-pt2).lengthsq();
			if (d > 10) {
			   status = false;
			   break;
			}
		     }
		  }
		  if (status == false)
		     break;
	       }
	    }
	    if (status == false)
	       break;
	 }
	 if (status == false)
	    break;
      }
   }
   return status;
}

// return a pair, the first is status (1 if the name was found, 0 if not)
// 
std::pair<bool, std::string>
coot::protein_geometry::get_monomer_name(const std::string &comp_id, int imol_enc) const {

   std::pair<bool, std::string> r(false,"");

   bool allow_minimal_flag = true;
   std::pair<bool, dictionary_residue_restraints_t> rest =
      get_monomer_restraints_internal(comp_id, imol_enc, allow_minimal_flag);

   if (rest.first) { 
      r.first = true;
      std::string s = rest.second.residue_info.name;
      r.second = coot::util::remove_trailing_whitespace(s);
   }
   return r;
} 




// Try comparing vs the comp_id first, if that fails compare the
// three_letter_code to the monomer_type.
//
// In future, try to come here only with the monomer_type adjusted to
// the comp_id, for example, monomer_type should be "NAG-b-D", not
// "NAG".
// 
std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints(const std::string &monomer_type,
					       int imol_enc) const {

   return get_monomer_restraints_internal(monomer_type, imol_enc, 0);
   
}

std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints_at_least_minimal(const std::string &monomer_type,
								int imol_enc) const {

   return get_monomer_restraints_internal(monomer_type, imol_enc, 1);
}

std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints_internal(const std::string &monomer_type,
							int imol_enc,
							bool allow_minimal_flag) const {

   // 20161028
   // Compiling with SRS causes a crash when we access dict_res_restraints.
   // Needs more testing.

   coot::dictionary_residue_restraints_t t(std::string("(null)"), 0);
   std::pair<bool, dictionary_residue_restraints_t> r(0,t);
   unsigned int nrest = dict_res_restraints.size();

   // This is how it used to be - starting from the beginning.
   // 20161004
   // Now we want to start from the end.
   // 
//    for (unsigned int i=0; i<nrest; i++) {
//       if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
// 	 if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
// 	    r.second = dict_res_restraints[i].second;
// 	    r.first = true;
// 	    break;
// 	 }
//       }
//    }

   // Compiler error.  Not sure why
   // std::vector<std::pair<int, dictionary_residue_restraints_t> >::reverse_iterator rit;
   // for (rit=dict_res_restraints.rbegin(); rit!=dict_res_restraints.rend(); rit--) {
   // }

   // Try to match most recent exact match by model number
   //
   for (int i=(nrest-1); i>=0; i--) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
 	 if (dict_res_restraints[i].first == imol_enc) {
 	    r.second = dict_res_restraints[i].second;
 	    r.first = true;
 	    break;
 	 }
      }
   }

   // Try to match any dictionary (this is what will match in most cases)
   //
   if (!r.first) {
      for (int i=(nrest-1); i>=0; i--) {
	 // std::cout << "Matching :" << dict_res_restraints[i].second.residue_info.comp_id
	 // << ": :" << monomer_type << ":" << std::endl;
	 if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	    // std::cout << "   comparing " << dict_res_restraints[i].first << " " << imol_enc << std::endl;
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       r.second = dict_res_restraints[i].second;
	       r.first = true;
	       // std::cout << "found " << std::endl;
	       break;
	    }
	 }
      }
   }
   
   if (!r.first) {
      // OK, that failed, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) { 
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    int ndict = dict_res_restraints.size();
	    for (int j=0; j<ndict; j++) {
	       if (dict_res_restraints[j].second.residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		  r.first = true;
		  r.second = dict_res_restraints[j].second;
		  break;
	       }
	    }
	 }
	 if (r.first)
	    break;
      }
   }

   if (!r.first) {
      for (unsigned int i=0; i<nrest; i++) {
	 if (dict_res_restraints[i].second.residue_info.three_letter_code == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       if (allow_minimal_flag || (! dict_res_restraints[i].second.is_bond_order_data_only())) { 
		  r.second = dict_res_restraints[i].second;
		  r.first = true;
		  break;
	       }
	    }
	 }
      }
   }
   return r;
}

// return -1 on monomer not found.
int
coot::protein_geometry::get_monomer_restraints_index(const std::string &monomer_type,
						     int imol_enc,
						     bool allow_minimal_flag) const {

   int r = -1;
   bool debug = false;

   unsigned int nrest = dict_res_restraints.size();
   for (unsigned int i=0; i<nrest; i++) {
      if (debug)
	 std::cout << "in get_monomer_restraints_index() comparing \""
		   << dict_res_restraints[i].second.residue_info.comp_id << "\" vs \"" << monomer_type
		   << "\" and " << dict_res_restraints[i].first << " " <<  imol_enc
		   << " with allow_minimal_flag " << allow_minimal_flag << std::endl;
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	    // if (dict_res_restraints[i].first == imol_enc) {
	    if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].second.is_bond_order_data_only())) { 
	       r = i;
	       break;
	    }
	 }
      }
   }

   if (r == -1) {
      for (unsigned int i=0; i<nrest; i++) {
	 if (dict_res_restraints[i].second.residue_info.comp_id  == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].second.is_bond_order_data_only())) { 
		  r = i;
		  break;
	       }
	    }
	 }
      }
   }
   
   if (r == -1) {
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) { 
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       int ndict = dict_res_restraints.size();
	       for (int j=0; j<ndict; j++) {
		  if (dict_res_restraints[j].second.residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		     r = j;
		     break;
		  }
	       }
	    }
	    if (r != -1)
	       break;
	 }
      }
   }
   
   if (r == -1) {
      for (unsigned int i=0; i<nrest; i++) {
	 if (dict_res_restraints[i].second.residue_info.three_letter_code  == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].second.is_bond_order_data_only())) { 
		  r = i;
		  break;
	       }
	    }
	 }
      }
   }
   
   return r;
}


std::string
coot::protein_geometry::get_type_energy(const std::string &atom_name,
					const std::string &residue_name,
					int imol) const { 
   // return "" if not found, else return the energy type found in ener_lib.cif
   //
   std::string r;
   int indx = get_monomer_restraints_index(residue_name, imol, 1); // in get_type_energy()
   if (indx != -1) {
      const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[indx].second;
      r = restraints.type_energy(atom_name);
   }
   return r;
}

// return -1.1 on failure to look up.
// 
double
coot::protein_geometry::get_vdw_radius(const std::string &atom_name,
				       const std::string &residue_name,
				       int imol,
				       bool use_vdwH_flag) const {
   double r = -1.1;
   int indx = get_monomer_restraints_index(residue_name, imol, 1); // in get_vdw_radius
   if (indx != -1) {
      const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[indx].second;
      std::string et = restraints.type_energy(atom_name);
      if (et != "") {
	 std::map<std::string, energy_lib_atom>::const_iterator it =
	    energy_lib.atom_map.find(et);
	 if (it != energy_lib.atom_map.end()) {
	    if (use_vdwH_flag)
	       r = it->second.vdwh_radius;
	    else
	       r = it->second.vdw_radius;
	 }
      }
   } else {
      std::cout << "  no restraints for type " << residue_name << std::endl;
   }
   return r;
}

// calculated once and then stored
bool
coot::protein_geometry::atom_is_metal(mmdb::Atom *atom) const {

   // PDBv3 FIXME

   bool status = false;
   std::string atom_name(atom->GetAtomName());
   if (atom_name == "NA" || atom_name == "CA" || atom_name == "LI") {
      return true;
   } else {
      if (atom_name == "BE" || atom_name == "K" || atom_name == "RB") {
	 return true;
      } else {
	 if (atom_name == "SR" || atom_name == "CS" || atom_name == "BA") {
	    return true;
	 } else {
	    if (atom_name == "SC" || atom_name == "TI" || atom_name == "V" || atom_name == "CR") {
	       return true;
	    } else {
	       if (atom_name == "MN" || atom_name == "FE" || atom_name == "CO" || atom_name == "NI") {
		  return true;
	       } else {
		  if (atom_name == "CU" || atom_name == "ZN" || atom_name == "ZR" || atom_name == "MO") {
		     return true;
		  } else {
		     if (atom_name == "AG" || atom_name == "AU" || atom_name == "PT" || atom_name == "HG") {
			return true;
		     } else {
			if (atom_name == "OS" || atom_name == "PB" || atom_name == " K" || atom_name == " W") {
			   return true;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   return status;
}


// expand to 4c, the atom_id, give that it should match an atom of the type comp_id.
// Used in chem mods, when we don't know the comp_id until residue modification time.
// 
std::string
coot::protein_geometry::atom_id_expand(const std::string &atom_id,
				       const std::string &comp_id,
				       int imol_enc) const {

   std::string s = atom_id;
   bool allow_minimal = true;
   int idx = get_monomer_restraints_index(comp_id, imol_enc, allow_minimal);
   if (idx != -1) {
      const coot::dictionary_residue_restraints_t &restraints =
	 dict_res_restraints[idx].second;
      const std::vector<dict_atom> &atoms = restraints.atom_info;
      for (unsigned int iat=0; iat<atoms.size(); iat++) { 
	 if (atoms[iat].atom_id == atom_id) { 
	    s = atoms[iat].atom_id_4c;
	    break;
	 }
      }
   }
//    std::cout << "atom_id_expand() \"" << atom_id << "\" expanded to \""
// 	     << s << "\"" << std::endl;
   return s; 
} 






// Hmmm... empty function, needs examining.
// 
// Use dynamic add if necessary.
// 
// Return -1 if the thing was not found or added.
int
coot::protein_geometry::get_monomer_type_index(const std::string &monomer_type) { 

   int i = -1;

   return i;
}

// Return 1 for hydrogen or deuterium, 0 for not found or not a hydrogen.
//
bool
coot::dictionary_residue_restraints_t::is_hydrogen(const std::string &atom_name) const {

   bool r = false;
   for (unsigned int i=0; i<atom_info.size(); i++) {
      if (false)
	 std::cout << "in is_hydrogen() comparing \"" << atom_info[i].atom_id_4c << "\" with \""
		   << atom_name << "\"" << std::endl;
      if (atom_info[i].atom_id_4c == atom_name) {
	 if (false)
	    std::cout << "in is_hydrogen found atom name " << atom_name << " and atom has type_symbol \""
		      << atom_info[i].type_symbol << "\"" << std::endl;
	 if (atom_info[i].type_symbol == "H" ||
	     atom_info[i].type_symbol == " H" ||
	     atom_info[i].type_symbol == "D") {
	    r = true;
	    break;
	 }
      }
   }
   return r;
}

bool
coot::dict_atom::is_hydrogen() const {

   bool r = false;
   if (type_symbol == "H" ||
       type_symbol == " H" ||
       type_symbol == "D")
      r = true;
   return r;
}

bool
coot::dictionary_residue_restraints_t::is_hydrogen(unsigned int idx) const {

   bool r = false;
   if (idx < atom_info.size()) {
      const std::string &ele = atom_info[idx].type_symbol;
      if (ele == " H" || ele == "H" || ele == "D")
	 r = true;
   } 
   return r;
}

unsigned int
coot::dictionary_residue_restraints_t::number_of_non_hydrogen_atoms() const {

   unsigned int r = 0;
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {
      if (! is_hydrogen(iat))
	 r++;
   }
   return r;
}

std::string
coot::dictionary_residue_restraints_t::get_bonded_atom(const std::string &H_atom_name) const {

  std::string r;

  for (unsigned int i=0; i<bond_restraint.size(); i++) {
    if (bond_restraint[i].atom_id_1_4c() == H_atom_name) {
      r = bond_restraint[i].atom_id_2_4c();
      break;
    }
    if (bond_restraint[i].atom_id_2_4c() == H_atom_name) {
      r = bond_restraint[i].atom_id_1_4c();
      break;
    }
  }
  return r;
}



// c.f. dict_torsion_restraint_t::is_ring_torsion()
bool
coot::dictionary_residue_restraints_t::is_ring_torsion(const coot::atom_name_quad &quad) const {

   bool match = false; 
   std::vector<std::string> torsion_atom_names(2);
   torsion_atom_names[0] = quad.atom_name(1);
   torsion_atom_names[1] = quad.atom_name(2);
   std::vector<std::vector<std::string> > ring_atoms_sets = get_ligand_ring_list();

   for (unsigned int iring=0; iring<ring_atoms_sets.size(); iring++) { 
      const std::vector<std::string> &ring_atom_names = ring_atoms_sets[iring];
      int n_match = 0;
      for (unsigned int iname_1=0; iname_1<ring_atom_names.size(); iname_1++) {
	 for (unsigned int iname_2=0; iname_2<torsion_atom_names.size(); iname_2++) { 
	    if (ring_atom_names[iname_1] == torsion_atom_names[iname_2])
	       n_match++;
	 }
      }
      if (n_match == 2) {
	 match = true;
	 break;
      }
   }

   return match;
} 

bool
coot::dictionary_residue_restraints_t::has_unassigned_chiral_volumes() const {
   bool r = 0;
   for (unsigned int ic=0; ic<chiral_restraint.size(); ic++) {
      if (chiral_restraint[ic].has_unassigned_chiral_volume()) {
	 r = 1;
	 break;
      }
   }
   return r;
}

bool
coot::dictionary_residue_link_restraints_t::has_unassigned_chiral_volumes() const {
   bool r = 0;
   for (unsigned int ic=0; ic<link_chiral_restraint.size(); ic++) {
      if (link_chiral_restraint[ic].has_unassigned_chiral_volume()) {
	 r = 1;
	 break;
      }
   }
   return r;
}


int
coot::dictionary_residue_restraints_t::assign_chiral_volume_targets() {

   int ich = 0;

   if (0) 
      std::cout << "DEBUG:: in dictionary_residue_restraints_t::assign_chiral_volume_targets "
		<< "there are " << chiral_restraint.size() << " chiral restraints for "
		<< residue_info.comp_id << " \n";
    
   for (unsigned int i=0; i<chiral_restraint.size(); i++) {
      chiral_restraint[i].assign_chiral_volume_target(bond_restraint, angle_restraint);
      ich++;
   }
   return ich;
}

int
coot::dictionary_residue_link_restraints_t::assign_link_chiral_volume_targets() {

   int ic = 0;
   for (unsigned int i=0; i<link_chiral_restraint.size(); i++) {
      std::vector <coot::dict_bond_restraint_t> bond_restraints_1;
      std::vector <coot::dict_bond_restraint_t> bond_restraints_2;
      std::vector <coot::dict_angle_restraint_t> angle_restraints_1;
      std::vector <coot::dict_angle_restraint_t> angle_restraints_2;
      std::vector <coot::dict_link_bond_restraint_t> link_bonds;
      std::vector <coot::dict_link_angle_restraint_t> link_angles;
      
      link_chiral_restraint[i].assign_chiral_volume_target(bond_restraints_1,
							   angle_restraints_1,
							   bond_restraints_2,
							   angle_restraints_2,
							   link_bonds,
							   link_angles);
      ic++;
   }
   return ic;
}

double
coot::dict_chiral_restraint_t::assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds,
							   const std::vector <dict_angle_restraint_t> &angles) {

   double vol = -1;
   double a = -1, b = -1, c = -1;
   double alpha = -1, beta = -1, gamma = -1;
   std::string mmdb_centre_atom =  atom_id_mmdb_expand(local_atom_id_centre);
   std::string mmdb_local_atom_id_1 = atom_id_mmdb_expand(local_atom_id_1);
   std::string mmdb_local_atom_id_2 = atom_id_mmdb_expand(local_atom_id_2);
   std::string mmdb_local_atom_id_3 = atom_id_mmdb_expand(local_atom_id_3);
   
   // local_atom_id_centre to local_atom_id_1 bond length
   for (unsigned int i=0; i<bonds.size(); i++) {
      try { 
	 if (bonds[i].atom_id_1_4c() == mmdb_centre_atom) {
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_1)) { 
	       a = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_2)) { 
	       b = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_3)) { 
	       c = bonds[i].value_dist();
	    }
	 }
	 if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_centre)) {
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_1)) { 
	       a = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_2)) { 
	       b = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_3)) { 
	       c = bonds[i].value_dist();
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 // do nothing, it's not really an error if the dictionary
	 // doesn't have target geometry (the bonding description came
	 // from a Chemical Component Dictionary entry for example).
      } 
   }

   for (unsigned int i=0; i<angles.size(); i++) {
      if (angles[i].atom_id_2_4c() == mmdb_centre_atom) {
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_2 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_3) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_2 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_3))  {
	    alpha = clipper::Util::d2rad(angles[i].angle());
	 }
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_3) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_3))  {
	    beta = clipper::Util::d2rad(angles[i].angle());
	 }
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_2) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_2))  {
	    gamma = clipper::Util::d2rad(angles[i].angle());
	 }
      }
   }

   
   if (a > 0 && b > 0 && c > 0) {
//       std::cout << "DEBUG:: found all distances in chiral restraint" << std::endl;
      if (alpha > 0 && beta > 0 && gamma > 0) {
// 	 std::cout << "DEBUG:: found all angles in chiral restraint" << std::endl;
	 vol = assign_chiral_volume_target_internal(a, b, c, alpha, beta, gamma);
      } else {
// 	 std::cout << "DEBUG:: failed to find all angles in chiral restraint"
// 		   << alpha << " " << beta << " " << gamma << std::endl;
      }
   } else {
//       std::cout << "DEBUG:: failed to find all distances in chiral restraint"
// 		<< a << " " << b << " " << c << std::endl;
   }
   return vol;
}

double
coot::dict_link_chiral_restraint_t::assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds_1,
								const std::vector <dict_angle_restraint_t> &angles_1,
								const std::vector <dict_bond_restraint_t> &bonds_2,
								const std::vector <dict_angle_restraint_t> &angles_2,
								const std::vector <dict_link_bond_restraint_t> &link_bonds,
								const std::vector <dict_link_angle_restraint_t> &link_angles) {

   double d = 0;

   return d;
} 


// angles in radians.
double
coot::dict_chiral_restraint_t::assign_chiral_volume_target_internal(double a, double b, double c,
								    double alpha, double beta, double gamma) {

   // set target_volume
   // from: abc ( 1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma) + 2(cos(alpha) + cos(beta) + cos(gamma)))^0.5

   double cos_alpha = cos(alpha);
   double cos_beta  = cos(beta);
   double cos_gamma = cos(gamma);
   
   double cos_2_alpha = cos_alpha * cos_alpha;
   double cos_2_beta  = cos_beta  * cos_beta;
   double cos_2_gamma = cos_gamma * cos_gamma;


   double param = 1-cos_2_alpha-cos_2_beta-cos_2_gamma + 2*cos_alpha*cos_beta*cos_gamma;
   if (false) {
      std::cout << "assign_chiral_volume_target_internal() input a, b, c, alpha, beta, gamma " << a << " "
		<< b << " " << c << " "
		<< clipper::Util::rad2d(alpha) << " "
		<< clipper::Util::rad2d(beta) << " "
		<< clipper::Util::rad2d(gamma) << " " << std::endl;

      double a_bit = 1-cos_2_alpha-cos_2_beta-cos_2_gamma;
      double b_bit = 2 * cos_alpha * cos_beta * cos_gamma;
      double c_bit = a_bit + b_bit;
      std::cout << "    bits: " << a_bit << " " << b_bit << " " << c_bit << std::endl;
      std::cout << "    param: " << param << std::endl;
   }
   if (param < 0) param = 0;
   target_volume_ = volume_sign * a*b*c*sqrt(param);

   // volume_sigma_ = 0.2;  // seems reasonable give target voluemes of about 2.6
   // volume_sigma_ = 0.3; // test
   volume_sigma_ = 0.18; // meh, 0.3 seemed to give "too many" chiral errors when
                         // refining large numbers of atoms in poor/EM maps.

   if (false)
      std::cout << "DEBUG:: assign_chiral_volume_target_internal() target_volume chiral: "
		<< target_volume_ << std::endl;

   // give the appropriate dictionary, this can return a target_volume_ of nan
   return target_volume_;
}


std::string
coot::protein_geometry::three_letter_code(const unsigned int &i) const {

   std::string r = dict_res_restraints[i].second.residue_info.three_letter_code;
   if (r == "")
      r = dict_res_restraints[i].second.residue_info.comp_id;
   return r;
}

// add "synthetic" 5 atom planar peptide restraint
void
coot::protein_geometry::add_planar_peptide_restraint() {

   std::string plane_id = "plane-5-atoms";
   mmdb::realtype dist_esd = 0.08; // was 0.11 (why?)

   std::string atom_id; 
   std::vector<std::pair<int, std::string> > v;
   v.push_back(std::pair<int, std::string> (1, "CA"));
   v.push_back(std::pair<int, std::string> (1, "C"));
   v.push_back(std::pair<int, std::string> (1, "O"));
   v.push_back(std::pair<int, std::string> (2, "N"));
   v.push_back(std::pair<int, std::string> (2, "CA"));

   for (unsigned int i=0; i<v.size(); i++) {
      link_add_plane("TRANS",  v[i].second, plane_id, v[i].first, dist_esd);
      link_add_plane("PTRANS", v[i].second, plane_id, v[i].first, dist_esd);
   }
}

bool
coot::protein_geometry::make_tight_planar_peptide_restraint() {

   std::string link_id("TRANS");
   std::string plane_id("plane-5-atoms");
   bool ifound = false;

   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) { // e.g "TRANS"
	 std::vector<coot::dict_link_plane_restraint_t>::iterator it;
	 for (it = dict_link_res_restraints[i].link_plane_restraint.begin();
	      it != dict_link_res_restraints[i].link_plane_restraint.end(); it++) {
	    if (it->plane_id == plane_id) {
	       it->set_dist_esd(0.03); // guess value
	       ifound = true;
	       break;
	    }
	 }
      }
   }
   return ifound;
}


void
coot::protein_geometry::remove_planar_peptide_restraint() {

   std::string link_id = "TRANS";
   std::string plane_id = "plane-5-atoms";
   unsigned int ifound = 0;

   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if ((dict_link_res_restraints[i].link_id ==  "TRANS")  ||
          (dict_link_res_restraints[i].link_id == "PTRANS")) {

	 std::vector<coot::dict_link_plane_restraint_t>::iterator it;
	 for (it = dict_link_res_restraints[i].link_plane_restraint.begin();
	      it != dict_link_res_restraints[i].link_plane_restraint.end(); it++) {
	    if (it->plane_id == plane_id) {
	       ifound++;
	       if (0)
		  std::cout << "INFO:: before removal of plane3 TRANS has " 
			    << dict_link_res_restraints[i].link_plane_restraint.size()
			    << " plane restraints\n";
	       
	       // let's remove it
 	       dict_link_res_restraints[i].link_plane_restraint.erase(it);
	       
	       if (0)
		  std::cout << "INFO::  after removal of plane3 TRANS has " 
			    << dict_link_res_restraints[i].link_plane_restraint.size()
			    << " plane restraints\n";
	       break;
	    }
	 }
      }
      if (ifound == 2)
	 break;
   }
}

// Do the link restraints contain a planar peptide restraint?
bool
coot::protein_geometry::planar_peptide_restraint_state() const {

   std::string link_id = "TRANS";
   std::string plane_id = "plane-5-atoms";
   bool ifound = false;

   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) { // e.g "TRANS"
	 
	 std::vector<coot::dict_link_plane_restraint_t>::const_iterator it;
	 for (it = dict_link_res_restraints[i].link_plane_restraint.begin();
	      it != dict_link_res_restraints[i].link_plane_restraint.end(); it++) {
	    if (it->plane_id == plane_id) {
	       ifound = true;
	       break;
	    }
	 }
      }
   }
   return ifound;
}


// restraints for omega for both CIS and TRANS links (and
// PTRANS)
void
coot::protein_geometry::add_omega_peptide_restraints() {

   double esd = 5.0; // perhaps this should be passed?
   std::vector<std::pair<std::string, double> > v;
   v.push_back(std::pair<std::string, double> ("TRANS",  180.0));
   v.push_back(std::pair<std::string, double> ("PTRANS", 180.0));
   v.push_back(std::pair<std::string, double> ("CIS",    0.0));
   v.push_back(std::pair<std::string, double> ("PCIS",   0.0));

   for (unsigned int iv=0; iv<v.size(); iv++) {
      std::string link_id = v[iv].first;
      // period is 0 (like the dictionary).  A good thing?
      link_add_torsion(link_id, 1, 1, 2, 2, "CA", "C", "N", "CA", v[iv].second, esd, 0, "omega");
   }

}



void
coot::protein_geometry::remove_omega_peptide_restraints() {

   std::vector<std::string> v;
   v.push_back("TRANS");
   v.push_back("PTRANS");
   v.push_back("CIS");
   v.push_back("PCIS");

   bool ifound = 0;
   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == v[i]) { // is TRANS, say

	 std::vector<coot::dict_link_torsion_restraint_t>::iterator it;
	 for (it = dict_link_res_restraints[i].link_torsion_restraint.begin();
	      it != dict_link_res_restraints[i].link_torsion_restraint.end(); it++) {
	    if (it->id() == "omega") {
	       ifound = 1;
 	       dict_link_res_restraints[i].link_torsion_restraint.erase(it);
	       break;
	    }
	 }
      }
   }
}



// a list of three-letter-codes (should that be comp_ids?) that match
// the string in the chem_comp name using the simple_monomer_descriptions
// 
std::vector<std::pair<std::string, std::string> >
coot::protein_geometry::matching_names(const std::string &test_string,
				       short int allow_minimal_descriptions) const {

   std::vector<std::pair<std::string, std::string> > v;
   std::map<std::string,coot::dictionary_residue_restraints_t>::const_iterator it;

   std::vector<std::string> test_name_fragments =
      util::split_string(test_string, " ");

   for (it=simple_monomer_descriptions.begin();
	it!=simple_monomer_descriptions.end();
	it++) {
      std::string name_dc = util::downcase(it->second.residue_info.name);

      std::string::size_type ifound = std::string::npos;      

      for (unsigned int i=0; i<test_name_fragments.size(); i++) {

	 const std::string &test_string = test_name_fragments[i];
	 std::string test_string_dc = util::downcase(test_string);

	 ifound = name_dc.find(test_string_dc);
	 if (ifound == std::string::npos)
	    break;
      }
      if (ifound != std::string::npos) {

	 // 	 std::cout << "test_string :" << test_string << ": matched :"
	 // 		   << it->second.residue_info.comp_id << ": :"
	 // 		   << it->second.residue_info.name
	 // 		   << ":" << std::endl;

	 std::pair<std::string, std::string> p(it->second.residue_info.comp_id,
					       it->second.residue_info.name);
	 v.push_back(p);
      }
   }
   return v;
}

void
coot::dictionary_residue_restraints_t::remove_redundant_plane_restraints() {


   bool match = true; // synthetic first value

   // erase_if usage would be more elegant here.
   //
   // This might be better done with indices, then we can remove the
   // higher planes (rather than the lower ones)
   // 
   while (match) {

      match = false;
      std::vector<dict_plane_restraint_t>::iterator it;
      for (it=plane_restraint.begin(); it!=plane_restraint.end(); it++) { 
	 if (is_redundant_plane_restraint(it)) {
	    // std::cout << "   erase plane " << it->plane_id << std::endl;
	    plane_restraint.erase(it);
	    match = true;
	    break;
	 }
      }
   }
} 

// is the plane restraint of it_ref redundant? (i.e. has an exact copy
// (i.e. do all the atom names match? another plane in the list?
// 
bool
coot::dictionary_residue_restraints_t::is_redundant_plane_restraint(std::vector<dict_plane_restraint_t>::iterator it_ref) const {

   bool match = false;
   std::vector<dict_plane_restraint_t>::const_iterator it_this;
   for (it_this=plane_restraint.begin(); it_this!=it_ref; it_this++) {
      
      if (it_this->n_atoms() >= it_ref->n_atoms()) {
	 
	 // do all of the atoms in this_rest have matchers in ref_rest?
	 //
	 int n_match = 0;
	 for (int i=0; i<it_ref->n_atoms(); i++) {
	    for (int j=0; j<it_this->n_atoms(); j++) {
	       if (it_this->atom_id(j) == it_ref->atom_id(i)) {
		  n_match++;
		  break;
	       }
	    }
	 }
	 if (n_match == it_ref->n_atoms()) {
	    if (0) { 
	       std::cout << "test plane     " << it_ref->plane_id << " matches list plane id "
			 << it_this->plane_id << " ref plane: ";
	       for (int iat=0; iat<it_ref->n_atoms(); iat++)
		  std::cout << " " << it_ref->atom_id(iat);
	       std::cout << " vs list-plane ";
	       for (int iat=0; iat<it_this->n_atoms(); iat++)
		  std::cout << " " << it_this->atom_id(iat);

	       std::cout << std::endl;
	    } 
	    match = true;
	    break;
	 }
      }
   }
   return match;
}

// if an atom is in more than one plane restraint, then
// reduce its esd.
void coot::dictionary_residue_restraints_t::reweight_subplanes() {

   std::map<std::string, int> name_map;
   std::vector<dict_plane_restraint_t>::iterator it;
   for (it=plane_restraint.begin(); it!=plane_restraint.end(); it++) {
      for (int i=0; i<it->n_atoms(); i++) {
	 const std::string &atom_name = it->atom_id(i);
	 name_map[atom_name]++;
      }
   }

//    std::map<std::string, int>::iterator it_names;
//    for(it_names=name_map.begin(); it_names!=name_map.end(); it_names++) {
//       if (it_names->second != 1) {
// 	 // reweight.
// 	 double w_multiplier = sqrt(double(it_names->second));
// 	 const std::string &atom_name = it_names->first;
// 	 for (it=plane_restraint.begin(); it!=plane_restraint.end(); it++) {
// 	    for (unsigned int i=0; i<it->n_atoms(); i++) {
// 	       if (it->atom_id(i) == atom_name) {
// 		  it->set_dist_esd(i, it->dist_esd(i) * w_multiplier);
// 	       }
// 	    }
// 	 }
//       }
//    }


   for (unsigned int i=0; i<plane_restraint.size(); i++) {
      dict_plane_restraint_t &rest_1 = plane_restraint[i];
      for (unsigned int j=0; j<plane_restraint.size(); j++) {
	 dict_plane_restraint_t &rest_2 = plane_restraint[j];
	 if (i != j) {
	    std::vector<int> matchers_1;
	    for (int ii=0; ii<rest_1.n_atoms(); ii++) {
	       for (int jj=0; jj<rest_2.n_atoms(); jj++) {
		  if (rest_1.atom_id(ii) == rest_2.atom_id(jj)) {
		     matchers_1.push_back(ii);
		  }
	       }
	    }
	    if (matchers_1.size() > 3) {
	       for (unsigned int im=0; im<matchers_1.size(); im++) {
		  rest_1.set_dist_esd(matchers_1[im], rest_1.dist_esd(matchers_1[im]) * 1.4142);
	       }
	    }
	 }
      }
   }
}




bool
coot::simple_cif_reader::has_restraints_for(const std::string &res_type) {

   bool r = 0;
   for (unsigned int i=0; i<three_letter_codes.size(); i++) {
      if (three_letter_codes[i] == res_type) {
	 r = 1;
	 break;
      }
   }
   return r;
}

// replace (return 1)  or add (if not replacable) (return 0).
// 
bool
coot::protein_geometry::replace_monomer_restraints(std::string monomer_type,
						   int imol_enc,
						   const coot::dictionary_residue_restraints_t &mon_res_in) {
   bool s = false;

   dictionary_residue_restraints_t mon_res = mon_res_in;
   mon_res.assign_chiral_volume_targets();
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    dict_res_restraints[i].second = mon_res;
	    s = true;
	    break;
	 }
      }
   }

   if (!s) {
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, mon_res_in);
      dict_res_restraints.push_back(p);
   }
   return s;
}


// Keep everything that we have already, replace only those
// parts that are in mon_res_in.
// Used to update bond and angle restraints from Mogul.
// 
bool
coot::protein_geometry::replace_monomer_restraints_conservatively(std::string monomer_type, 
								  const dictionary_residue_restraints_t &mon_res) {

   bool s = false;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 replace_monomer_restraints_conservatively_bonds( i, mon_res);
	 replace_monomer_restraints_conservatively_angles(i, mon_res);
	 s = true;
	 break;
      }
   }
   return s;
}

void
coot::protein_geometry::replace_monomer_restraints_conservatively_bonds(int irest,
									const dictionary_residue_restraints_t &mon_res) {

   // replace bonds
   // 
   for (unsigned int ibond=0; ibond<dict_res_restraints[irest].second.bond_restraint.size(); ibond++) { 
      for (unsigned int jbond=0; jbond<mon_res.bond_restraint.size(); jbond++) {

	 if (dict_res_restraints[irest].second.bond_restraint[ibond].atom_id_1_4c() ==
	     mon_res.bond_restraint[jbond].atom_id_1_4c()) {
	    if (dict_res_restraints[irest].second.bond_restraint[ibond].atom_id_2_4c() ==
		mon_res.bond_restraint[jbond].atom_id_2_4c()) {
	       dict_res_restraints[irest].second.bond_restraint[ibond] =
		  mon_res.bond_restraint[jbond];
	       break;
	    }
	 }

	 if (dict_res_restraints[irest].second.bond_restraint[ibond].atom_id_1_4c() ==
	     mon_res.bond_restraint[jbond].atom_id_2_4c()) {
	    if (dict_res_restraints[irest].second.bond_restraint[ibond].atom_id_2_4c() ==
		mon_res.bond_restraint[jbond].atom_id_1_4c()) {
	       dict_res_restraints[irest].second.bond_restraint[ibond] =
		  mon_res.bond_restraint[jbond];
	       break;
	    }
	 }
      }
   }
}




void
coot::protein_geometry::replace_monomer_restraints_conservatively_angles(int irest,
									 const dictionary_residue_restraints_t &mon_res) {

   for (unsigned int iangle=0; iangle<dict_res_restraints[irest].second.angle_restraint.size(); iangle++) { 
      for (unsigned int jangle=0; jangle<mon_res.angle_restraint.size(); jangle++) {

	 // middle atom the same
	 // 
	 if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_2_4c() ==
	     mon_res.angle_restraint[jangle].atom_id_2_4c()) {

	    // check for either way round of 1 and 3:

	    if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_1_4c() ==
		mon_res.angle_restraint[jangle].atom_id_1_4c()) {
	       if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_3_4c() ==
		   mon_res.angle_restraint[jangle].atom_id_3_4c()) {
		  dict_res_restraints[irest].second.angle_restraint[iangle] =
		     mon_res.angle_restraint[jangle];
	       }
	    }

	    if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_1_4c() ==
		mon_res.angle_restraint[jangle].atom_id_3_4c()) {
	       if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_3_4c() ==
		   mon_res.angle_restraint[jangle].atom_id_1_4c()) {
		  dict_res_restraints[irest].second.angle_restraint[iangle] =
		     mon_res.angle_restraint[jangle];
	       }
	    }
	 }
      }
   }
}

// replace the restraints that we have with new_restraints,
// keeping restraints that in the current set bu not in
// new_restraints
void
coot::dictionary_residue_restraints_t::conservatively_replace_with(const dictionary_residue_restraints_t &new_restraints) {

   conservatively_replace_with_bonds(new_restraints);
   conservatively_replace_with_angles(new_restraints);

}

void
coot::dictionary_residue_restraints_t::conservatively_replace_with_bonds (const dictionary_residue_restraints_t &new_restraints) {

   for (unsigned int ibond=0; ibond<bond_restraint.size(); ibond++) { 
      for (unsigned int jbond=0; jbond<new_restraints.bond_restraint.size(); jbond++) {
	 if (bond_restraint[ibond].atom_id_1_4c() ==
	     new_restraints.bond_restraint[jbond].atom_id_1_4c()) {
	    if (bond_restraint[ibond].atom_id_2_4c() ==
		new_restraints.bond_restraint[jbond].atom_id_2_4c()) {
	       bond_restraint[ibond] = new_restraints.bond_restraint[jbond];
	       break;
	    }
	 } 

	 if (bond_restraint[ibond].atom_id_1_4c() ==
	     new_restraints.bond_restraint[jbond].atom_id_2_4c()) {
	    if (bond_restraint[ibond].atom_id_2_4c() ==
		new_restraints.bond_restraint[jbond].atom_id_1_4c()) {
	       bond_restraint[ibond] = new_restraints.bond_restraint[jbond];
	       break;
	    }
	 } 
      }
   }
}

void
coot::dictionary_residue_restraints_t::conservatively_replace_with_angles(const dictionary_residue_restraints_t &new_restraints) {

   for (unsigned int iangle=0; iangle<angle_restraint.size(); iangle++) { 
      for (unsigned int jangle=0; jangle<new_restraints.angle_restraint.size(); jangle++) {
	 
	 if (angle_restraint[iangle].atom_id_1_4c() ==
	     new_restraints.angle_restraint[jangle].atom_id_1_4c()) {
	    if (angle_restraint[iangle].atom_id_2_4c() ==
		new_restraints.angle_restraint[jangle].atom_id_2_4c()) {
	       if (angle_restraint[iangle].atom_id_3_4c() ==
		   new_restraints.angle_restraint[jangle].atom_id_3_4c()) {
		  angle_restraint[iangle] = new_restraints.angle_restraint[jangle];
	       }
	    }
	 }
	 if (angle_restraint[iangle].atom_id_3_4c() ==
	     new_restraints.angle_restraint[jangle].atom_id_1_4c()) {
	    if (angle_restraint[iangle].atom_id_2_4c() ==
		new_restraints.angle_restraint[jangle].atom_id_2_4c()) {
	       if (angle_restraint[iangle].atom_id_1_4c() ==
		   new_restraints.angle_restraint[jangle].atom_id_3_4c()) {
		  angle_restraint[iangle] = new_restraints.angle_restraint[jangle];
	       }
	    }
	 }
      }
   }
}

void
coot::dictionary_residue_restraints_t::replace_coordinates(const dictionary_residue_restraints_t &mon_res_in) {

   for (unsigned int iat=0; iat<atom_info.size(); iat++) { 
      dict_atom &at = atom_info[iat];
      
      for (unsigned int iat=0; iat<mon_res_in.atom_info.size(); iat++) { 
	 const dict_atom &at_ref = mon_res_in.atom_info[iat];

	 if (at_ref.atom_id_4c == at.atom_id_4c) {
	    at.pdbx_model_Cartn_ideal = at_ref.pdbx_model_Cartn_ideal;
	    at.model_Cartn            = at_ref.model_Cartn;
	 } 
      }
   }
} 




std::vector<std::string>
coot::protein_geometry::monomer_types() const {
   std::vector<std::string> v;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      v.push_back(dict_res_restraints[i].second.residue_info.comp_id);
   }
   return v;
}

// Thow an exception if we can't get the group of r.
//
// 20100420 First compare against the three_letter_code and if not
// found try testing against the comp_id (added 20100420).  This is
// needed so that this function returns the group from a Ur/U (that's
// the comp_id/tlc) [and similarly for other bases].  This is needed
// so that find_link_type_rigourous() works for link type "p" (RNA/DNA
// stuff).  (Problem found when trying to sphere refine a on RNA -
// 3l0u).
// 
std::string
coot::protein_geometry::get_group(mmdb::Residue *r) const {

   std::string res_name = r->GetResName();
   return get_group(res_name);
}


std::string
coot::protein_geometry::get_group(const std::string &res_name_in) const {
   
   bool found = 0;
   std::string group;
   std::string res_name = res_name_in;
   if (res_name.length() > 3)
      res_name = res_name.substr(0,2);
   for (unsigned int i=0; i<size(); i++) {
      if (three_letter_code(i) == res_name) {
	 found = 1;
	 group = (*this)[i].second.residue_info.group;
	 break;
      }
   }

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].second.residue_info.comp_id == res_name) {
	 found = 1;
	 group = dict_res_restraints[i].second.residue_info.group;
	 break;
      }
   }

   if (! found) {
      std::string s = "No dictionary group found for residue type :";
      s += res_name;
      s += ":";
      throw std::runtime_error(s);
   }
   return group;
}


std::vector<std::string>
coot::protein_geometry::residue_names_with_no_dictionary(mmdb::Manager *mol, int imol_no) const {

   std::vector<std::string> v;

   if (mol) {
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       if (residue_p) {
		  std::string residue_name = residue_p->GetResName();
		  if (! have_dictionary_for_residue_type_no_dynamic_add(residue_name))
		     v.push_back(residue_name);
	       }
	    }
	 }
      }
   }
   return v;
}

bool
coot::protein_geometry::read_extra_dictionaries_for_molecule(mmdb::Manager *mol, int imol, int *read_number_p) {

   std::vector<std::string> v = residue_names_with_no_dictionary(mol, imol);
   bool success = true;
   for (std::size_t i=0; i<v.size(); i++) {
      const std::string &rn = v[i];
      int success_for_residue = try_dynamic_add(rn, *read_number_p);
      if (success_for_residue == 0)
	 success = false;
      *read_number_p += 1;
   }

   return success;
}



// optional arg: bool try_autoload_if_needed=true.
// optional arg: float b_factor.
mmdb::Residue *
coot::protein_geometry::get_residue(const std::string &comp_id, int imol_enc,
				    bool idealised_flag,
				    bool try_autoload_if_needed, float b_factor) {

   // If the coordinates for the model are (0,0,0) then this function
   // returns a null.
   
   mmdb::Residue *residue_p = NULL;

   // might use try_dynamic_add (if needed).
   bool r = have_dictionary_for_residue_type(comp_id, imol_enc, try_autoload_if_needed);
   if (r) {
      for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
	 if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	    residue_p = dict_res_restraints[i].second.GetResidue(idealised_flag, b_factor);
	    break;
	 }
      }
   }
   return residue_p;
} 



mmdb::Manager *
coot::protein_geometry::mol_from_dictionary(const std::string &three_letter_code,
					    int imol_enc,
					    bool idealised_flag) {

   mmdb::Manager *mol = NULL;
   mmdb::Residue *residue_p = get_residue(three_letter_code, imol_enc, idealised_flag);
   if (residue_p) { 
      mmdb::Chain *chain_p = new mmdb::Chain;
      chain_p->SetChainID("A");
      chain_p->AddResidue(residue_p);
      mmdb::Model *model_p = new mmdb::Model;
      model_p->AddChain(chain_p);
      mol = new mmdb::Manager;
      mol->AddModel(model_p);
   } else {
      std::cout << "WARNING:: Null residue in mol_from_dictionary() for " << three_letter_code << std::endl;
   }
   return mol;
}

mmdb::Manager *
coot::protein_geometry::mol_from_dictionary(int monomer_index,
					    int imol_enc,
					    bool idealised_flag) {

   mmdb::Manager *mol = NULL;
   mmdb::Residue *residue_p = 0;
   float b_factor = 30.0;
   int r_size = dict_res_restraints.size();

   if (monomer_index >= 0)
      if (monomer_index < r_size)
	 residue_p = dict_res_restraints[monomer_index].second.GetResidue(idealised_flag, b_factor);

   if (residue_p) { 
      mmdb::Chain *chain_p = new mmdb::Chain;
      chain_p->SetChainID("A");
      chain_p->AddResidue(residue_p);
      mmdb::Model *model_p = new mmdb::Model;
      model_p->AddChain(chain_p);
      mol = new mmdb::Manager;
      mol->AddModel(model_p);
   } else {
      std::cout << "WARNING:: Null residue in mol_from_dictionary() for idx "
		<< monomer_index << std::endl;
   }

   std::cout << "DEBUG:: mol_from_dictionary() returns " << mol << std::endl;
   return mol;
}


void
coot::protein_geometry::print_chem_links() const {

   // simple vector
   // for (unsigned int i_chem_link=0; i_chem_link<chem_link_vec.size(); i_chem_link++)
   // std::cout<< i_chem_link << " " << chem_link_vec[i_chem_link] << "\n";

   std::map<unsigned int, std::vector<chem_link> >::const_iterator it;

   for (it=chem_link_map.begin(); it!=chem_link_map.end(); it++) {
      const std::vector<chem_link> &v = it->second;
      std::vector<chem_link>::const_iterator itv;
      for (itv=v.begin(); itv!=v.end(); itv++) {
	 const chem_link &cl = *itv;
	 std::cout << "   " << it->first << " " << cl << "\n";
      }
   }

} 

// delete comp_id from dict_res_restraints (if it exists there).
bool
coot::protein_geometry::delete_mon_lib(const std::string &comp_id, int imol_enc) {

   bool deleted = false;
   std::vector<std::pair<int, coot::dictionary_residue_restraints_t> >::iterator it;
   for (it=dict_res_restraints.begin(); it!=dict_res_restraints.end(); it++) {
      if (it->second.residue_info.comp_id == comp_id) {
	 if (it->first == imol_enc) {
	    dict_res_restraints.erase(it);
	    deleted = 1;
	    break;
	 }
      }
   }
   
   return deleted;
} 

bool
coot::protein_geometry::linkable_residue_types_p(const std::string &this_res_type,
						 const std::string &env_res_type) {

   std::pair<short int, coot::dictionary_residue_restraints_t> r1 = get_monomer_restraints(this_res_type, IMOL_ENC_ANY);
   std::pair<short int, coot::dictionary_residue_restraints_t> r2 = get_monomer_restraints(env_res_type, IMOL_ENC_ANY);

   bool r = 0;
   if (r1.first) {
      if (r1.second.residue_info.group != "non-polymer")
	 r = 1;
   }
   if (r2.first) {
      if (r2.second.residue_info.group != "non-polymer")
	 r = 1;
   }
   return r;
} 

bool
coot::protein_geometry::OXT_in_residue_restraints_p(const std::string &residue_type) const {

   bool r = 0;
   std::pair<bool, coot::dictionary_residue_restraints_t> p =
      get_monomer_restraints(residue_type, IMOL_ENC_ANY);
   if (p.first) {
      for (unsigned int i=0; i<p.second.atom_info.size(); i++) {
	 if (p.second.atom_info[i].atom_id_4c == " OXT") {
	    r = 1;
	    break;
	 }
      }
   } else {
      if (0) 
	 std::cout << "INFO:: residue type :" << residue_type << ": not found in dictionary"
		   << std::endl;
   } 
   return r;
}

mmdb::Residue *
coot::dictionary_residue_restraints_t::GetResidue(bool idealised_flag, float b_factor) const {

   mmdb::Residue *residue_p = NULL;
   std::vector<mmdb::Atom *> atoms;

   bool make_hetatoms = ! coot::util::is_standard_residue_name(residue_info.comp_id);
   int atom_index = 0;
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {

      clipper::Coord_orth p(0,0,0);
      bool flag_and_have_coords = false;

      if (idealised_flag && atom_info[iat].pdbx_model_Cartn_ideal.first) {
	 p = atom_info[iat].pdbx_model_Cartn_ideal.second;
	 flag_and_have_coords = true;
      }

      if (! flag_and_have_coords) {
	 // OK, try model_Cartn (and that is idealised if the dictionary was refmac)
	 // (better than nothing).
	 // 
	 if (atom_info[iat].model_Cartn.first) {
	    p = atom_info[iat].model_Cartn.second;
	    flag_and_have_coords = true;
	 }
      }

      if (flag_and_have_coords) { 
	 mmdb::Atom *atom = new mmdb::Atom;
	 mmdb::realtype occ = 1.0;
	 mmdb::realtype b = b_factor;
	 std::string ele = atom_info[iat].type_symbol; // element
	 atom->SetCoordinates(p.x(), p.y(), p.z(), occ, b);
	 atom->SetAtomName(atom_info[iat].atom_id_4c.c_str());

	 //	       Strange things happen when I use this...
	 // 	       atom->SetAtomName(atom_index, -1,
	 // 				 atom_info[iat].atom_id_4c.c_str(),
	 // 				 "", "", ele.c_str());
	       
	 atom->SetElementName(ele.c_str());
	 if (make_hetatoms)
	    atom->Het = 1;
	 atoms.push_back(atom);
	 atom_index++; // for next round
      }
   }

   if (atoms.size() > 0) {
      residue_p = new mmdb::Residue;
      residue_p->SetResID(residue_info.comp_id.c_str(), 1, "");
      for (unsigned int iat=0; iat<atoms.size(); iat++) 
	 residue_p->AddAtom(atoms[iat]);
   } 
   return residue_p; 
} 

std::vector<std::vector<std::string> >
coot::dictionary_residue_restraints_t::get_ligand_ring_list() const {

   // get a list of bonds, so that they can be used to find rings.
   // 
   std::vector<std::pair<std::string, std::string> > bonds;
   for (unsigned int irest=0; irest<bond_restraint.size(); irest++) {
      std::pair<std::string, std::string> p(bond_restraint[irest].atom_id_1_4c(),
					    bond_restraint[irest].atom_id_2_4c());
      bonds.push_back(p);
   }

   // used in non-necessarily-aromatic way...
   coot::aromatic_graph_t bond_list(bonds);
   std::vector<std::vector<std::string> > ring_list = bond_list.ring_list();

   if (0) {
      std::cout << "----------- " << ring_list.size() << " rings ---------- " << std::endl;
      for (unsigned int i=0; i<ring_list.size(); i++) {
	 std::cout << "ring " << i << "\n   ";
	 for (unsigned int j=0; j<ring_list[i].size(); j++) { 
	    std::cout << ring_list[i][j] << "  ";
	 }
	 std::cout << std::endl;
      }
   }
   return ring_list;
}

bool
coot::dictionary_residue_restraints_t::in_same_ring(const std::string &atom_name_1, const std::string &atom_name_2) const { 

   bool match = false;
   std::vector<std::vector<std::string> > ring_list = get_ligand_ring_list();

   for (unsigned int i=0; i<ring_list.size(); i++) {
      unsigned int n_match = 0;
      for (unsigned int j=0; j<ring_list[i].size(); j++) {
	 if (ring_list[i][j] == atom_name_1)
	    n_match++;
	 if (ring_list[i][j] == atom_name_2)
	    n_match++;
      }
      if (n_match == 2) {
	 match = true;
	 break;
      }
   }
   return match;
}

bool
coot::dictionary_residue_restraints_t::in_same_ring(const std::string &atom_name_1, const std::string &atom_name_2,
						    const std::vector<std::vector<std::string> > &ring_list) const {

   bool match = false;

   for (unsigned int i=0; i<ring_list.size(); i++) {
      unsigned int n_match = 0;
      for (unsigned int j=0; j<ring_list[i].size(); j++) {
	 if (ring_list[i][j] == atom_name_1)
	    n_match++;
	 if (ring_list[i][j] == atom_name_2)
	    n_match++;
      }
      if (n_match == 2) {
	 match = true;
	 break;
      }
   }
   return match;
}



bool
coot::dictionary_residue_restraints_t::ligand_has_aromatic_bonds_p() const {

   for (unsigned int irest=0; irest<bond_restraint.size(); irest++)
      if (bond_restraint[irest].type() == "aromatic")
	 return true;
   return false;
}


std::vector<std::vector<std::string> >
coot::dictionary_residue_restraints_t::get_ligand_aromatic_ring_list() const {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   // 
   std::vector<std::pair<std::string, std::string> > bonds;
   for (unsigned int irest=0; irest<bond_restraint.size(); irest++) {
      if (bond_restraint[irest].type() == "aromatic") {
	 std::pair<std::string, std::string> p(bond_restraint[irest].atom_id_1_4c(),
					       bond_restraint[irest].atom_id_2_4c());
	 bonds.push_back(p);
      }
   }
   
   coot::aromatic_graph_t arom(bonds);
   std::vector<std::vector<std::string> > ring_list = arom.ring_list();

   if (0) {
      std::cout << "----------- " << ring_list.size() << " rings ---------- " << std::endl;
      for (unsigned int i=0; i<ring_list.size(); i++) {
	 std::cout << "ring " << i << "\n   ";
	 for (unsigned int j=0; j<ring_list[i].size(); j++) { 
	    std::cout << ring_list[i][j] << "  ";
	 }
	 std::cout << std::endl;
      }
   }
   return ring_list;
}



// Find the bonded neighbours of the given atoms - throw an
// exception if residue name is not in dictionary.
// 
std::vector<std::string>
coot::protein_geometry::get_bonded_neighbours(const std::string &residue_name,
					      int imol_enc,
					      const std::string &atom_name_1,
					      const std::string &atom_name_2, 
					      bool also_2nd_order_neighbs_flag) const {

   std::vector<std::string> v;

   std::vector<std::string> v_2nd_order; // only filled for triple bonds

   std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
      get_monomer_restraints_at_least_minimal(residue_name, imol_enc);

   if (restraints.first) { 
      for (unsigned int i=0; i<restraints.second.bond_restraint.size(); i++) {
// 	 std::cout << "-- comparing " << atom_name_1 << " " << atom_name_2 << " -- to -- "
// 		   << restraints.second.bond_restraint[i].atom_id_1_4c() << " "
// 		   << restraints.second.bond_restraint[i].atom_id_2_4c() << std::endl;
	 if (restraints.second.bond_restraint[i].atom_id_1_4c() == atom_name_1)
	    if (restraints.second.bond_restraint[i].atom_id_2_4c() != atom_name_2) {
	       // std::cout << " adding a " << restraints.second.bond_restraint[i].atom_id_2_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_2_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name,  imol_enc,
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
	 if (restraints.second.bond_restraint[i].atom_id_1_4c() == atom_name_2)
	    if (restraints.second.bond_restraint[i].atom_id_2_4c() != atom_name_1) { 
	       // std::cout << " adding b " << restraints.second.bond_restraint[i].atom_id_2_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_2_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name, imol_enc,
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
	 if (restraints.second.bond_restraint[i].atom_id_2_4c() == atom_name_1)
	    if (restraints.second.bond_restraint[i].atom_id_1_4c() != atom_name_2) {
	       // std::cout << " adding c " << restraints.second.bond_restraint[i].atom_id_1_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_1_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name, imol_enc,
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
	 if (restraints.second.bond_restraint[i].atom_id_2_4c() == atom_name_2) {
	    if (restraints.second.bond_restraint[i].atom_id_1_4c() != atom_name_1) {
	       // std::cout << " adding d " << restraints.second.bond_restraint[i].atom_id_1_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_1_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name, imol_enc,
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
	 }
      }

      // add the neighbour neighbours to v (if needed):
      if (also_2nd_order_neighbs_flag) { 
	 for(unsigned int in=0; in<v_2nd_order.size(); in++)
	    if (std::find(v.begin(), v.end(), v_2nd_order[in]) == v.end())
	       v.push_back(v_2nd_order[in]);
      } 

      if (v.size()) {
	 // add the initial atom names if they are not already there
	 if (std::find(v.begin(), v.end(), atom_name_1) == v.end())
	    v.push_back(atom_name_1);
	 if (std::find(v.begin(), v.end(), atom_name_2) == v.end())
	    v.push_back(atom_name_2);
      } 
   } else {
      std::string m = "No dictionary for ";
      m += residue_name;
      throw std::runtime_error(m);
   } 
   return v;
}

std::vector<std::pair<std::string, std::string> >
coot::protein_geometry::get_bonded_and_1_3_angles(const std::string &comp_id, int imol) const {

   std::vector<std::pair<std::string, std::string> > v;
   int idx = get_monomer_restraints_index(comp_id, imol, true);
   if (idx != -1) {
      const dictionary_residue_restraints_t &rest = dict_res_restraints[idx].second;
      for (unsigned int i=0; i<rest.bond_restraint.size(); i++) {
	 std::pair<std::string, std::string> p(rest.bond_restraint[i].atom_id_1_4c(),
					       rest.bond_restraint[i].atom_id_2_4c());
	 v.push_back(p);
      }
      for (unsigned int i=0; i<rest.angle_restraint.size(); i++) {
	 std::pair<std::string, std::string> p(rest.angle_restraint[i].atom_id_1_4c(),
					       rest.angle_restraint[i].atom_id_3_4c());
	 v.push_back(p);
      }
   }
   return v;
} 



// return "" on not found
std::string
coot::dictionary_residue_restraints_t::get_bond_type(const std::string &name_1,
						     const std::string &name_2) const {

   std::string r("unknown");
   for (unsigned int i=0; i<bond_restraint.size(); i++) {
      if (bond_restraint[i].atom_id_1_4c() == name_1) { 
	 if (bond_restraint[i].atom_id_2_4c() == name_2) { 
	    r = bond_restraint[i].type();
	    break;
	 }
      }
      if (bond_restraint[i].atom_id_1_4c() == name_2) { 
	 if (bond_restraint[i].atom_id_2_4c() == name_1) { 
	    r = bond_restraint[i].type();
	    break;
	 }
      }
   }
   return r;
}

// should this return a set - or filter the results for duplicates?
//
std::vector<std::string>
coot::protein_geometry::monomer_restraints_comp_ids() const {

   std::vector<std::string> v;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++)
      v.push_back(dict_res_restraints[i].second.residue_info.comp_id);
   return v;
} 


// can throw a std::runtime_error
std::string
coot::protein_geometry::Get_SMILES_for_comp_id(const std::string &comp_id) const {

   bool found = false;
   std::string s; 
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {

      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {

	 unsigned int nd = dict_res_restraints[i].second.descriptors.descriptors.size();
	 for (unsigned int idesc=0; idesc<nd; idesc++) { 
	    if (dict_res_restraints[i].second.descriptors.descriptors[idesc].type == "SMILES_CANONICAL") { 
	       s = dict_res_restraints[i].second.descriptors.descriptors[idesc].descriptor;
	       found = true;
	       break;
	    }
	 }
      }
      if (found)
	 break;
   }

   if (! found){
      // check non-canonical
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {

      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {

	 unsigned int nd = dict_res_restraints[i].second.descriptors.descriptors.size();
	 for (unsigned int idesc=0; idesc<nd; idesc++) { 
	    if (dict_res_restraints[i].second.descriptors.descriptors[idesc].type == "SMILES") {
	       s = dict_res_restraints[i].second.descriptors.descriptors[idesc].descriptor;
	       found = true;
	       break;
	    }
	 }
      }
      if (found)
	 break;
   }
}

   if (! found){
      std::string mess = "No SMILES in dictionary for ";
      mess += comp_id;
      throw (std::runtime_error(mess));
   }
   return s;
} 

      



// This uses have_dictionary_for_residue_type() (and thus
// try_dynamic_add() if needed).
// 
// Return -1 if residue type not found.
// 
int
coot::protein_geometry::n_hydrogens(const std::string &residue_type) {

   int n_hydrogens = -1;
   
   std::pair<bool, dictionary_residue_restraints_t> r =
      get_monomer_restraints(residue_type, IMOL_ENC_ANY);

   if (r.first) {
      n_hydrogens = 0; // not not-found
      for (unsigned int iat=0; iat<r.second.atom_info.size(); iat++) { 
	 if (r.second.atom_info[iat].type_symbol == " H")
	    n_hydrogens++;
	 else 
	 if (r.second.atom_info[iat].type_symbol == "H")
	    n_hydrogens++;
      }
   }
   return n_hydrogens;
} 

// This uses have_dictionary_for_residue_type() (and thus
// try_dynamic_add() if needed).
// 
// Return -1 if residue type not found.
// 
int
coot::protein_geometry::n_non_hydrogen_atoms(const std::string &residue_type) {

   int n_atoms = -1;
   
   std::pair<bool, dictionary_residue_restraints_t> r =
      get_monomer_restraints(residue_type, IMOL_ENC_ANY);

   if (r.first) {
      n_atoms = 0; // not not-found
      for (unsigned int iat=0; iat<r.second.atom_info.size(); iat++) { 
	 if (r.second.atom_info[iat].type_symbol != " H")
	    if (r.second.atom_info[iat].type_symbol != "H")
	       n_atoms++;
      }
   }
   return n_atoms;
}




coot::energy_lib_atom
coot::protein_geometry::get_energy_lib_atom(const std::string &ener_type) const {

   coot::energy_lib_atom at;
   
   std::map<std::string, energy_lib_atom>::const_iterator it =
      energy_lib.atom_map.find(ener_type);

   if (it != energy_lib.atom_map.end())
      at = it->second;
   
   return at;
   
}

// use auto-load if not present
void
coot::protein_geometry::use_unimodal_ring_torsion_restraints(int imol, const std::string &res_name,
							     int mmcif_read_number) {

   bool minimal = false;
   int idx = get_monomer_restraints_index(res_name, imol, minimal);
   if (idx == -1) {
      try_dynamic_add(res_name, mmcif_read_number);
      idx = get_monomer_restraints_index(res_name, imol, minimal);
   }

   if (idx != -1) {
      // continue

      // (first is the imol)
      std::vector <dict_torsion_restraint_t> &torsion_restraints =
	 dict_res_restraints[idx].second.torsion_restraint;

      // We don't want to clear them all, just the ring torsions
      // orig_torsion_restraints.clear();
      std::vector<std::string> ring_atom_names;
      ring_atom_names.push_back(" C1 "); ring_atom_names.push_back(" C2 ");
      ring_atom_names.push_back(" C3 "); ring_atom_names.push_back(" C4 ");
      ring_atom_names.push_back(" C5 "); ring_atom_names.push_back(" O5 ");
      if (res_name == "XYP")
	 for (unsigned int i=0; i<ring_atom_names.size(); i++)
	    ring_atom_names[i][3] = 'B';

      if (false) {
	 for (unsigned int i=0; i<torsion_restraints.size(); i++)
	    std::cout << " torsion restraint: " << i << " " << torsion_restraints[i] << std::endl;
	 std::cout << "...............  pre-delete size: " << torsion_restraints.size()
		   << " for " << res_name << std::endl;
      }

      torsion_restraints.erase(std::remove_if(torsion_restraints.begin(),
					      torsion_restraints.end(),
					      restraint_eraser(ring_atom_names)), torsion_restraints.end());

      if (false)
	 std::cout << "............... post-delete size: " << torsion_restraints.size()
		   << " for " << res_name << std::endl;

      std::vector<atom_name_torsion_quad> quads = get_reference_monomodal_torsion_quads(res_name);
      for (unsigned int i=0; i<quads.size(); i++) {
	 const atom_name_torsion_quad &quad = quads[i];
	 dict_torsion_restraint_t tors(quad.id,
				       quad.atom_name(0), quad.atom_name(1), quad.atom_name(2), quad.atom_name(3),
				       quad.torsion, 4.0, 1);
	 torsion_restraints.push_back(tors);
      }
   }
}

// pass the atom names and the desired torsion value - sigma is not specified
// by the user.
void
coot::protein_geometry::use_unimodal_ring_torsion_restraints(int imol, const std::string &res_name,
							     const std::vector<coot::atom_name_torsion_quad> &tors_info_vec,
							     int mmcif_read_number) {

   bool minimal = false;
   int idx = get_monomer_restraints_index(res_name, imol, minimal);
   if (idx == -1) {
      try_dynamic_add(res_name, mmcif_read_number);
      idx = get_monomer_restraints_index(res_name, imol, minimal);
   }

   if (idx != -1) {
      // continue

      // (first is the imol)
      std::vector <dict_torsion_restraint_t> &torsion_restraints =
	 dict_res_restraints[idx].second.torsion_restraint;

      // We don't want to clear all torsion restraints, just the ring torsions
      //
      std::set<std::string> ring_atom_names;
      // we generate the ring atoms names from the atoms in the passed quads
      for (std::size_t i=0; i<tors_info_vec.size(); i++) {
	 const atom_name_torsion_quad &q = tors_info_vec[i];
	 ring_atom_names.insert(q.atom_name(0));
	 ring_atom_names.insert(q.atom_name(1));
	 ring_atom_names.insert(q.atom_name(2));
	 ring_atom_names.insert(q.atom_name(3));
      }

      if (false) {
	 std::cout << "...............  pre-delete size: " << torsion_restraints.size()
		   << " for " << res_name << std::endl;
      }

      torsion_restraints.erase(std::remove_if(torsion_restraints.begin(),
					      torsion_restraints.end(),
					      restraint_eraser(ring_atom_names)), torsion_restraints.end());

      if (false)
	 std::cout << "............... post-delete size: " << torsion_restraints.size()
		   << " for " << res_name << std::endl;

      for (unsigned int i=0; i<tors_info_vec.size(); i++) {
	 const atom_name_torsion_quad &quad = tors_info_vec[i];
	 // this could have a nicer interface - using the quad here.
	 dict_torsion_restraint_t tors(quad.id,
				       quad.atom_name(0), quad.atom_name(1), quad.atom_name(2), quad.atom_name(3),
				       quad.torsion, 4.0, 1);
	 torsion_restraints.push_back(tors);
      }
   }
}


std::vector<coot::atom_name_torsion_quad>
coot::protein_geometry::get_reference_monomodal_torsion_quads(const std::string &res_name) const {

   std::vector<coot::atom_name_torsion_quad> v;
   if (res_name == "MAN" || res_name == "GLC" || res_name == "GAL") {
      v.push_back(coot::atom_name_torsion_quad("var_4",  "O5", "C5", "C4", "C3", -61.35));
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1", "O5", "C5", "C4",  67.67));
      v.push_back(coot::atom_name_torsion_quad("var_11", "C5", "O5", "C1", "C2", -67.59));
      v.push_back(coot::atom_name_torsion_quad("var_9",  "C3", "C2", "C1", "O5",  61.20));
      v.push_back(coot::atom_name_torsion_quad("var_6",  "C5", "C4", "C3", "C2",  53.75));
   }
   if (res_name == "BMA") {
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1", "O5", "C5", "C4",  67.57));
      v.push_back(coot::atom_name_torsion_quad("var_4",  "O5", "C5", "C4", "C3", -61.31));
      v.push_back(coot::atom_name_torsion_quad("var_6",  "C5", "C4", "C3", "C2",  53.80));
      v.push_back(coot::atom_name_torsion_quad("var_9",  "C3", "C2", "C1", "O5",  61.28));
      v.push_back(coot::atom_name_torsion_quad("var_11", "C5", "O5", "C1", "C2", -67.52));
   }
   if (res_name == "NAG") {
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1", "O5", "C5", "C4",  61.18));
      v.push_back(coot::atom_name_torsion_quad("var_5",  "O5", "C5", "C4", "C3", -57.63));
      v.push_back(coot::atom_name_torsion_quad("var_7",  "C5", "C4", "C3", "C2",  57.01));
      v.push_back(coot::atom_name_torsion_quad("var_10", "C3", "C2", "C1", "O5",  57.60));
      v.push_back(coot::atom_name_torsion_quad("var_1",  "C5", "O5", "C1", "C2", -61.19));
   }
   if (res_name == "FUC") {
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1", "O5", "C5", "C4", -67.61));
      v.push_back(coot::atom_name_torsion_quad("var_5",  "C5", "C4", "C3", "C2", -53.87));
      v.push_back(coot::atom_name_torsion_quad("var_11", "C5", "O5", "C1", "C2",  67.595));
      v.push_back(coot::atom_name_torsion_quad("var_9",  "C3", "C2", "C1", "O5", -61.26));
      v.push_back(coot::atom_name_torsion_quad("var_4",  "O5", "C5", "C4", "C3",  61.32));
   }
   if (res_name == "XYP") {
      v.push_back(coot::atom_name_torsion_quad("var_4",  "O5B", "C5B", "C4B", "C3B", -61.35));
      v.push_back(coot::atom_name_torsion_quad("var_2",  "C1B", "O5B", "C5B", "C4B",  67.67));
      v.push_back(coot::atom_name_torsion_quad("var_11", "C5B", "O5B", "C1B", "C2B", -67.59));
      v.push_back(coot::atom_name_torsion_quad("var_9",  "C3B", "C2B", "C1B", "O5B",  61.20));
      v.push_back(coot::atom_name_torsion_quad("var_6",  "C5B", "C4B", "C3B", "C2B",  53.75));
   }

   return v;
}


// debug
void
coot::protein_geometry::debug() const {

   std::cout << "### debug(): " << dict_res_restraints.size() << " entries " << std::endl;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      int imol = dict_res_restraints[i].first;
      std::string imol_str = std::string("          ") + util::int_to_string(imol);
      if (imol == IMOL_ENC_ANY)
	 imol_str = "IMOL_ENC_ANY";
      if (imol == IMOL_ENC_AUTO)
	 imol_str = "IMOL_ENC_AUTO";
      if (imol == IMOL_ENC_UNSET)
	 imol_str = "IMOL_ENC_UNSET";
      std::cout << "     " << i << " imol: " << imol_str << " \""
		<< dict_res_restraints[i].second.residue_info << "\"" << std::endl;
   }
}

void
coot::dictionary_residue_restraints_t::remove_phosphate_hydrogens() {

   remove_PO4_SO4_hydrogens(" P");

}

void
coot::dictionary_residue_restraints_t::remove_sulphate_hydrogens() {

   remove_PO4_SO4_hydrogens(" S");
}

void
coot::dictionary_residue_restraints_t::remove_PO4_SO4_hydrogens(const std::string &P_ele) {

   // This is needed to make the dictionary match the modifications to the RDKit molecule
   // when removing Hydrogen atoms from Oxygen atoms on phosphates (and we need to do that
   // so that the atom types match the Acedrg tables (no bonds to Hydrogens in phosphates).

   // find Ps
   // find Os connected to Ps
   // find Hs connected to Os.
   // delete bond, angle, etc, restraints that contain those Hs
   std::vector<std::string> H_atoms_to_be_deleted;

   unsigned int n_atoms = atom_info.size();
   for (unsigned int i=0; i<n_atoms; i++) {
      if (element(atom_info[i].atom_id_4c) == P_ele) {

	 // this block needs reverse indexing check also - e.g. the P can be the second atom
	 //

	 std::vector<std::string> oxygen_list;
	 // is there a bond from an O to this P?
	 unsigned int n_bonds = bond_restraint.size();
	 for (unsigned int j=0; j<n_bonds; j++) {
	    const dict_bond_restraint_t &br = bond_restraint[j];
	    // is an atom of this bond the phosphate atom?
	    if (br.atom_id_1_4c() == atom_info[i].atom_id_4c) {
	       // yes it is.  Is there an oxygen atom bonded to this phosphate atom?
	       for (unsigned int k=0; k<n_bonds; k++) {
		  if (j != k) {
		     const dict_bond_restraint_t &br_inner = bond_restraint[k];
		     if (br_inner.atom_id_1_4c() == atom_info[i].atom_id_4c) {
			// is this an oxygen on the other side of the bond?
			if (element(br_inner.atom_id_2_4c()) == " O") {
			   // add it if it was not already in the list:
			   if (std::find(oxygen_list.begin(),
					 oxygen_list.end(),
					 br_inner.atom_id_2_4c()) == oxygen_list.end()) {
			      if (false)
				 std::cout << "adding " << util::single_quote(br_inner.atom_id_2_4c())
					   << std::endl;
			      oxygen_list.push_back(br_inner.atom_id_2_4c());
			   }
			}
		     }
		  }
	       }
	    }
	 }

	 if (oxygen_list.size() > 1) { // 20170603 unsure about this test
	    if (false)
	       std::cout << "found " << oxygen_list.size() << " oxygen atoms attached to "
			 << util::single_quote(atom_info[i].atom_id_4c) << std::endl;
	    // All these oxygen atoms are bonded to the same phosphate atoms
	    // delete all hydrogen atoms attached to all these oxygen atoms
	    //
	    std::vector<std::string> hydrogen_atom_delete_list;
	    std::vector<std::string> oxygen_atom_charge_list;
	    for (unsigned int j=0; j<n_bonds; j++) {
	       const dict_bond_restraint_t &br = bond_restraint[j];
	       if (std::find(oxygen_list.begin(),
			     oxygen_list.end(),
			     br.atom_id_1_4c()) != oxygen_list.end()) {
		  if (element(br.atom_id_2_4c()) == " H") {
		     hydrogen_atom_delete_list.push_back(br.atom_id_2_4c());
		     oxygen_atom_charge_list.push_back(br.atom_id_1_4c());
		  }
	       }
	       // reverse indexing
	       if (std::find(oxygen_list.begin(),
			     oxygen_list.end(),
			     br.atom_id_2_4c()) != oxygen_list.end()) {
		  if (element(br.atom_id_1_4c()) == " H") {
		     hydrogen_atom_delete_list.push_back(br.atom_id_1_4c());
		     oxygen_atom_charge_list.push_back(br.atom_id_2_4c());
		  }
	       }
	    }
	    if (false) {
	       std::cout << "Delete these " << hydrogen_atom_delete_list.size() << " hydrogen atoms"
			 << std::endl;
	       std::cout << "Charge these " << oxygen_atom_charge_list.size() << " oxygen atoms"
			 << std::endl;
	    }

	    delete_atoms_from_restraints(hydrogen_atom_delete_list);

	    for (unsigned int j=0; j<oxygen_atom_charge_list.size(); j++) {
	       for (unsigned int k=0; k<atom_info.size(); k++) {
		  if (atom_info[k].atom_id_4c == oxygen_atom_charge_list[j]) {
		     atom_info[k].formal_charge.first = true;
		     atom_info[k].formal_charge.first = -1;
		  }
	       }
	    }
	 }
      }
   }
   
}

void
coot::dictionary_residue_restraints_t::remove_carboxylate_hydrogens() {

   // This is a pain.

   std::vector<std::string> H_atoms_to_be_deleted;
   std::vector<std::string> oxygen_atom_charge_list;

   unsigned int n_atoms = atom_info.size();
   for (unsigned int i=0; i<n_atoms; i++) {

      if (element(atom_info[i].atom_id_4c) == " C") { // PDBv3 FIXME
	 std::vector<std::string> oxygen_list;
	 int n_bonds_to_C = 0;
	 int C_O_bond_idx = -1; // initially unset
	 std::string O_with_possibly_H;
	 int n_single_bonds = 0;
	 int n_double_bonds = 0;
	 unsigned int n_bonds = bond_restraint.size();
	 for (unsigned int j=0; j<n_bonds; j++) {
	    const dict_bond_restraint_t &br = bond_restraint[j];

	    //
	    // note to self - do reverse indexing also
	    //

	    // is an atom of this bond the C atom?
	    //
	    if (br.atom_id_1_4c() == atom_info[i].atom_id_4c) {
	       // yes it is.
	       n_bonds_to_C++;

	       std::string atom_id_O = br.atom_id_2_4c();
	       if (element(atom_id_O) == " O") {
		  if (br.type() == "single") {
		     C_O_bond_idx = j;
		     O_with_possibly_H = atom_id_O;
		     oxygen_atom_charge_list.push_back(O_with_possibly_H);
		     n_single_bonds++;
		  }
		  if (br.type() == "double") {
		     n_double_bonds++;
		  }
	       }
	    }
	    if (br.atom_id_2_4c() == atom_info[i].atom_id_4c) {
	       // yes it is.
	       n_bonds_to_C++;

	       std::string atom_id_O = br.atom_id_1_4c();
	       if (element(atom_id_O) == " O") {
		  if (br.type() == "single") {
		     C_O_bond_idx = j;
		     O_with_possibly_H = atom_id_O;
		     oxygen_atom_charge_list.push_back(O_with_possibly_H);
		     n_single_bonds++;
		  }
		  if (br.type() == "double") {
		     n_double_bonds++;
		  }
	       }
	    }
	 }

	 if (n_bonds_to_C == 3) {
	    if (n_single_bonds == 1) {
	       if (n_double_bonds == 1) {
		  if (! O_with_possibly_H.empty()) {
		     // go through the bond restraints looking for both the named Oxygen atom
		     // and a hydrogen.
		     std::string delete_H_atom;
		     for (unsigned int k=0; k<n_bonds; k++) {
			const dict_bond_restraint_t &br = bond_restraint[k];
			if (br.atom_id_1_4c() == O_with_possibly_H) {
			   if (element(br.atom_id_2_4c()) == " H") {
			      // delete this hydrogen and this bond
			      H_atoms_to_be_deleted.push_back(br.atom_id_2_4c());
			   }
			}
			if (br.atom_id_2_4c() == O_with_possibly_H) {
			   if (element(br.atom_id_1_4c()) == " H") {
			      // delete this hydrogen and this bond
			      H_atoms_to_be_deleted.push_back(br.atom_id_1_4c());
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   std::cout << "Here with H_atoms_to_be_deleted size() " << H_atoms_to_be_deleted.size() << std::endl;
   if (H_atoms_to_be_deleted.size() > 0) {
      delete_atoms_from_restraints(H_atoms_to_be_deleted);
      for (unsigned int j=0; j<oxygen_atom_charge_list.size(); j++) {
	 for (unsigned int k=0; k<atom_info.size(); k++) {
	    if (atom_info[k].atom_id_4c == oxygen_atom_charge_list[j]) {
	       atom_info[k].formal_charge.first = true;
	       atom_info[k].formal_charge.first = -1;
	    }
	 }
      }
   }
}


void
coot::dictionary_residue_restraints_t::delete_atoms_from_restraints(const std::vector<std::string> &hydrogen_atom_delete_list) {

   if (hydrogen_atom_delete_list.size()) {
      for (unsigned int j=0; j<hydrogen_atom_delete_list.size(); j++) {
	 atom_info.erase(std::remove_if(atom_info.begin(), atom_info.end(), eraser(hydrogen_atom_delete_list)),
			 atom_info.end());
	 bond_restraint.erase(std::remove_if(bond_restraint.begin(),
					     bond_restraint.end(),
					     eraser(hydrogen_atom_delete_list)),
			      bond_restraint.end());
	 angle_restraint.erase(std::remove_if(angle_restraint.begin(),
					      angle_restraint.end(),
					      eraser(hydrogen_atom_delete_list)),
			       angle_restraint.end());
      }
   }
}


// function here: convert_plane_restraints_to_improper_dihedrals()
// called by the user after having read the dictionary
// (and before refinement starts).

std::vector<coot::atom_name_quad>
coot::dictionary_residue_restraints_t::plane_restraint_to_improper_dihedrals(unsigned int idx) const {

   std::vector<atom_name_quad> q;
   std::map<std::string, bool> done_quads;

   q.reserve(4); // may need tweaking

   if (idx < plane_restraint.size()) {
      const dict_plane_restraint_t &pr = plane_restraint[idx];
      std::set<std::string> plane_rest_atom_name_map;

      for (int i=0; i<pr.n_atoms(); i++)
	 plane_rest_atom_name_map.insert(pr[i].first);

      int n_angle_restraints = angle_restraint.size();

      // now makes sets for all of the angle restraints

      std::vector<std::set<std::string> > angle_restraint_atom_names_map_vector(n_angle_restraints);

      if (false)
         std::cout << "plane_rest_atom_name_map has " << plane_rest_atom_name_map.size()
		   << " atoms " << std::endl;
      std::set<std::string>::const_iterator it;
      for (it=plane_rest_atom_name_map.begin(); it!=plane_rest_atom_name_map.end(); it++)
	 std::cout << "   " << *it;
      std::cout << std::endl;

      for (int i=0; i<n_angle_restraints; i++) {
	 const dict_angle_restraint_t &ar = angle_restraint[i];
	 for (int j=0; j<n_angle_restraints; j++) {
	    angle_restraint_atom_names_map_vector[i].insert(ar.atom_id_1_4c());
	    angle_restraint_atom_names_map_vector[i].insert(ar.atom_id_2_4c());
	    angle_restraint_atom_names_map_vector[i].insert(ar.atom_id_3_4c());
	 }
      }

      for (int i=0; i<n_angle_restraints; i++) {
	 const dict_angle_restraint_t &ar_1 = angle_restraint[i];

	 // are the atoms of that angle restraint in the plane restraint?
	 bool success = true;
	 for(it =angle_restraint_atom_names_map_vector[i].begin();
	     it!=angle_restraint_atom_names_map_vector[i].end();
	     it++) {
	    if (plane_rest_atom_name_map.find(*it) == plane_rest_atom_name_map.end()) {
	       success = false;
	       break;
	    } else {
	       // std::cout << "found i " << *it << " in plane rest atom map " << std::endl;
	    }
	 }
	 if (! success) continue;

	 for (int j=(i+1); j<n_angle_restraints; j++) {
	    const dict_angle_restraint_t &ar_2 = angle_restraint[j];

	    // are the atoms of this second angle_restraint in the plane restraint?
	    bool success_inner = true;
	    for(it =angle_restraint_atom_names_map_vector[j].begin();
		it!=angle_restraint_atom_names_map_vector[j].end();
		it++) {
	       if (plane_rest_atom_name_map.find(*it) == plane_rest_atom_name_map.end()) {
		  success_inner = false;
		  break;
	       }
	    }
	    if (! success_inner) continue;

	    // OK, do ar_1 and ar_2 share 2 atoms?

	    int n_match = 0;
	    for (it= angle_restraint_atom_names_map_vector[j].begin();
		 it!=angle_restraint_atom_names_map_vector[j].end();
		 it++) {
	       if (angle_restraint_atom_names_map_vector[i].find(*it) !=
		   angle_restraint_atom_names_map_vector[i].end()) {
		  n_match++;
	       }
	    }
	    if (n_match == 2) {
	       // i and j are have atoms that are in pr

	       // std::cout << "make a quad from " << angle_restraint[i] << " " << angle_restraint[j]
	       // << std::endl;

	       std::set<std::string> quad_atom_names;
	       for (it= angle_restraint_atom_names_map_vector[i].begin();
		    it!=angle_restraint_atom_names_map_vector[i].end();
		    it++) {
		  // std::string a = *it;
                  // std::cout << "a " << a << std::endl;
		  quad_atom_names.insert(*it);
	       }
	       for (it= angle_restraint_atom_names_map_vector[j].begin();
		    it!=angle_restraint_atom_names_map_vector[j].end();
		    it++) {
		  // const std::string &b = *it;
                  // std::cout << "b " << b << std::endl;
		  quad_atom_names.insert(*it);
	       }

	       if (quad_atom_names.size() == 4) {

                  // we need to put the atoms in to the atom name quad in the correct order
                  // which is not the order of ti quad_atom_names set.

                  std::string atom_name_1 = angle_restraint[i].atom_id_1_4c();
                  std::string atom_name_2 = angle_restraint[i].atom_id_2_4c();
                  std::string atom_name_3 = angle_restraint[i].atom_id_3_4c();

                  // what is the atom in the second angle restraint that is not in the first?
                  std::string other_atom;
	          for (it= angle_restraint_atom_names_map_vector[j].begin();
		       it!=angle_restraint_atom_names_map_vector[j].end();
		       it++) {
		     const std::string &b = *it;
                     if (angle_restraint_atom_names_map_vector[i].find(b) == angle_restraint_atom_names_map_vector[i].end()) {
                        other_atom = b;
                        if (false) { // debugging
                           std::cout << "Other atom \"" << b << "\" was not in the set ";
                           std::set<std::string>::const_iterator it_2;
                           for (it_2=angle_restraint_atom_names_map_vector[i].begin(); it_2!=angle_restraint_atom_names_map_vector[i].end(); it_2++)
                               std::cout << " \"" << *it_2 << "\"";
                           std::cout << std::endl;
                        }
                     }
                  }
                  if (!other_atom.empty()) {
                     std::string atom_name_4 = other_atom;

                     // make a key for the done quads from the sorted atom names
                     std::set<std::string> set_for_key;
                     set_for_key.insert(atom_name_1);
                     set_for_key.insert(atom_name_2);
                     set_for_key.insert(atom_name_3);
                     set_for_key.insert(atom_name_4);
                     std::set<std::string>::const_iterator it_s = set_for_key.begin();
                     std::string key = *it_s;
                     for(int ii=0; ii<3; ii++) {
                        it_s++;
                        key += "+" + *it_s;
                     }

		     if (done_quads.find(key) == done_quads.end()) {
		        atom_name_quad qn(atom_name_1, atom_name_2, atom_name_3, atom_name_4);
                        done_quads[key] = true;
			q.push_back(qn);
			// std::cout << "......... added this quad " << qn << std::endl;
		     }
		  }
	       } else {
		  std::cout << "error: weird quad " << std::endl;
	       }
	    }
	 }
      }
   }
   return q;
}

void coot::protein_geometry::all_plane_restraints_to_improper_dihedrals() {

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      dict_res_restraints[i].second.improper_dihedral_restraint.clear();
      for (unsigned int j=0; j<dict_res_restraints[i].second.plane_restraint.size(); j++) {
         std::vector<atom_name_quad> qs = dict_res_restraints[i].second.plane_restraint_to_improper_dihedrals(j);
         for(unsigned int k=0; k<qs.size(); k++) {
            dict_improper_dihedral_restraint_t r(qs[k].atom_name(0), qs[k].atom_name(1), qs[k].atom_name(2), qs[k].atom_name(3));
            dict_res_restraints[i].second.improper_dihedral_restraint.push_back(r);
         }
      }
   }

}

void coot::protein_geometry::delete_plane_restraints() {

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      dict_res_restraints[i].second.plane_restraint.clear();
   }

}

