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

bool
coot::dictionary_residue_restraints_t::is_bond_to_hydrogen_atom(const coot::dict_bond_restraint_t &br) const {
   bool is_H = false;
   std::string ele_1 = element(br.atom_id_1_4c());
   std::string ele_2 = element(br.atom_id_2_4c());
   if (ele_1 == " H") {
      is_H = true;
   } else {
      if (ele_2 == " H") {
         is_H = true;
      }
   }
   return is_H;
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
   nuclear_distances_flag = false;

   if (residue_p) {
      // mmdb::PResidue res = 0;
      mmdb::math::Graph    graph;
      mmdb::math::PPVertex V;
      mmdb::math::PPEdge   E;
      int       i, nV,nE, k1,k2;

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

   if (false)
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




mmdb::Residue *
coot::dictionary_residue_restraints_t::GetResidue(bool idealised_flag, float b_factor) const {

   // std::cout << "in GetResidue() idealised_flag is " << idealised_flag << std::endl;

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

      if (false) {  // debug
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for(int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::cout << "debug:: GetResidue() " << iat << " " << at->GetAtomName()
                         << at->x << " " << at->y << " " << at->z << std::endl;
            }
         }

      }
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
      for (it=plane_rest_atom_name_map.begin(); it!=plane_rest_atom_name_map.end(); ++it)
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
	     ++it) {
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
		++it) {
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
                           for (it_2=angle_restraint_atom_names_map_vector[i].begin(); it_2!=angle_restraint_atom_names_map_vector[i].end(); ++it_2)
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
                        ++it_s;
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
