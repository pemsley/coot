/* lidia-core/rdkit-interface.cc
 * 
 * Copyright 2010, 2011, 2012 by The University of Oxford
 * Copyright 2012, 2013, 2014, 2015, 2016 by Medical Research Council
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
#include <string>
#include "clipper/core/coords.h"
#include "clipper/minimol/minimol.h"
#include "geometry/protein-geometry.hh"
#include "lidia-core/use-rdkit.hh"
#include "mmdb2/mmdb_mattype.h"
#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <cstring>  // Fixes ::strchr complaints on 4.4.7 (hal)
#include <queue>
#include "utils/coot-utils.hh"
#include "rdkit-interface.hh"
#include <GraphMol/Chirality.h>  // for CIP ranks
#include "geometry/residue-and-atom-specs.hh"
#include "geometry/mol-utils.hh"

// what does this give us?
// #include "coot-utils/coot-coord-utils.hh" // after rdkit-interface.hh to avoid ::strchr problems
#include "neighbour-sorter.hh"


// This can throw an runtime_error exception (residue not in
// dictionary).
//
// Return an RDKit::RWMol with one conformer (perhaps we should do a
// conformer for each alt conf).
//
// alt_conf is an optional arg.
// 
RDKit::RWMol
coot::rdkit_mol(mmdb::Residue *residue_p, int imol_enc, const coot::protein_geometry &geom) {

   // std::cout << "--------------- starting rdkit_mol() ----- " << std::endl;

   if (! residue_p) {
      throw std::runtime_error("Null residue in coot::rdkit_mol()");
   } else {

      std::string res_name = residue_p->GetResName();
      if (false)
         std::cout << "====================  here in rdkit_mol() with geometry with res_name \""
                   << res_name << "\"" << std::endl;

      std::pair<bool, coot::dictionary_residue_restraints_t> p = 
         geom.get_monomer_restraints_at_least_minimal(res_name, imol_enc);
      if (! p.first) {

         std::string m = "rdkit_mol(): residue type ";
         m += res_name;
         m += " not in dictionary";
         throw(std::runtime_error(m));

      } else {
         if (false)
            std::cout << "......... calling rdkit_mol() with restraints that have "
                      << p.second.bond_restraint.size() << " bond restraints"
                      << std::endl;
         return rdkit_mol(residue_p, p.second);
      }
   }
}

// alt_conf is an optional argument, default "". If non-empty, it means "only select
// conformers with alt-confs that match `alt_conf`.
// undelocalize is an optional argument, default false
RDKit::RWMol
coot::rdkit_mol(mmdb::Residue *residue_p,
                const coot::dictionary_residue_restraints_t &restraints,
                const std::string &alt_conf,
                bool do_undelocalize) {

   auto get_alt_confs_in_residue = [] (mmdb::Residue *residue_p) {
      std::vector<std::string> v;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string alt_conf = at->altLoc;
            if (std::find(v.begin(), v.end(), alt_conf) == v.end()) {
               v.push_back(alt_conf);
            }
         }
      }
      return v;
   };

   // std::cout << "--------------- starting second rdkit_mol() ----- " << std::endl;

   bool debug = false;

   if (debug)
      std::cout << "=========== in rdkit_mol() with restraints that have "
                << restraints.atom_info.size() << " atoms, "
                << restraints.bond_restraint.size() << " bond restraints with do_undelocalize "
                << do_undelocalize
                << " for residue " << residue_spec_t(residue_p) << " and alt conf "
                << "\"" << alt_conf << "\"" << std::endl;

   if (debug)
      for (unsigned int ii=0; ii<restraints.atom_info.size(); ii++)
         std::cout << ii << "   " << restraints.atom_info[ii] << std::endl;

   RDKit::RWMol m;

   std::string n = coot::util::remove_trailing_whitespace(restraints.residue_info.name);
   m.setProp("_Name", n);
   m.setProp("ResName", std::string(residue_p->GetResName()));
   m.setProp("ResNumber", residue_p->GetSeqNum());
   m.setProp("ChainID", std::string(residue_p->GetChainID()));
   m.setProp("alt_id", alt_conf);

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   // this is so that we don't add multiple copies of an atom with
   // the same name (that is, only add the first atom of a given
   // atom name in a residue with alt confs).
   // 
   // We also need added_atoms so that that is what we iterate over
   // when handling bond restraints (we don't want to consider B
   // conformer for restraints if they are ignored when we are adding
   // atoms).
   // 
   std::vector<std::string> added_atom_names;
   std::vector<mmdb::Atom *>     added_atoms; // gets added to as added_atom_names gets added to.
   std::map<std::string, int> atom_index;
   int current_atom_id = 0;
   std::vector<std::pair<int, int> > bonded_atoms; // vector of the indices of the atoms that we will
                                       // add to the rdkit molecule. The first index is into the
                                       // residue_atoms and the second into restraints.atom_info.
                                       // We don't want to add atoms that are
                                       // not bonded to anything (e.g. hydrogens with mismatching
                                       // names).
  
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat_1=0; iat_1<n_residue_atoms; iat_1++) {
      mmdb::Atom *at_1 = residue_atoms[iat_1];
      if (! at_1->Ter) {
         std::string atom_name_1(at_1->name);
         std::string atom_alt_conf(at_1->altLoc);
         if (debug)
            std::cout << "rdkit_mol() handling atom " << iat_1 << " of " << n_residue_atoms
                      << " with mmdb::Residue atom name " << atom_name_1
                      << " alt-conf \"" << atom_alt_conf << "\""<< std::endl;
         if (true) {
            bool found_a_bonded_atom = false;
            for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
               if (restraints.bond_restraint[ib].atom_id_1_4c() == atom_name_1) {
                  // is the atom to which atom_name_1 is bonded in the
                  // atoms of the residue?
                  for (int iat_2=0; iat_2<n_residue_atoms; iat_2++) {
                     mmdb::Atom *at_2 = residue_atoms[iat_2];
                     std::string atom_name_2 = at_2->name;
                     if (atom_name_2 == restraints.bond_restraint[ib].atom_id_2_4c()) {
                        found_a_bonded_atom = true;
                        break;
                     } 
                  }
               } 
               if (restraints.bond_restraint[ib].atom_id_2_4c() == atom_name_1) {
                  // is the atom to wich atom_name_1 is bonded in the
                  // atoms of the residue?
                  for (int iat_2=0; iat_2<n_residue_atoms; iat_2++) {
                     mmdb::Atom *at_2 = residue_atoms[iat_2];
                     std::string atom_name_2 = at_2->name;
                     if (atom_name_2 == restraints.bond_restraint[ib].atom_id_1_4c()) {
                        found_a_bonded_atom = true;
                        break;
                     }
                  }
               }
            }

            // what is the index of this atom name in the atom_info? (we need this so that we
            // can apply the charge to the rdkit atom later.)
            int ai_idx = -1;
            for (unsigned int jj=0; jj<restraints.atom_info.size(); jj++) {
               if (restraints.atom_info[jj].atom_id_4c == atom_name_1) {
                  ai_idx = jj;
                  break;
               }
            }

            if (found_a_bonded_atom) {
               if (ai_idx > -1) {
                  std::pair<int, int> p(iat_1, ai_idx);
                  bool already_there = false;
                  for (unsigned int iba=0; iba<bonded_atoms.size(); iba++) {
                     if (bonded_atoms[iba].second == ai_idx) {
                        already_there = true;
                        break;
                     }
                  }
                  if (! already_there) {
                     bonded_atoms.push_back(p);
                  } else {
                     std::cout << "WARNING:: rdkit_mol() atom index already there " << ai_idx << std::endl;
                  }
               }
            }
         }
      }
   }

   if (debug) { 
      std::cout << "DEBUG:: number of bonded atoms with alt conf \"" << alt_conf << "\" found: "
                << bonded_atoms.size() << std::endl;
   }

   if (! bonded_atoms.empty()) {

      for (unsigned int iat=0; iat<bonded_atoms.size(); iat++) {

         mmdb::Atom *at = residue_atoms[bonded_atoms[iat].first];
         std::string atom_name(at->name);
         if (debug)
            std::cout << "   handling atom " << iat << " of " << n_residue_atoms << " bonded_atoms " 
                      << atom_name << " ";
         
         // only add the atom if the atom_name is not in the list of
         // already-added atom names.
         if (std::find(added_atom_names.begin(), added_atom_names.end(), atom_name) != added_atom_names.end()) {
            std::cout << "!!!! Problem? atom name \"" << atom_name
                      << "\" was already added" << std::endl;

         } else {
            RDKit::Atom *rdkit_at = new RDKit::Atom;
            try {
               std::string ele_capped =
                  coot::util::capitalise(coot::util::remove_leading_spaces(at->element));
               int atomic_number = tbl->getAtomicNumber(ele_capped);
               rdkit_at->setAtomicNum(atomic_number);
               // rdkit_at->setMass(tbl->getAtomicWeight(atomic_number));
               rdkit_at->setIsotope(0);
               rdkit_at->setProp("name", atom_name);

               // formal charge
               const coot::dict_atom &atom_info = restraints.atom_info[bonded_atoms[iat].second];
               if (false) {
                  std::cout << "in rdkit_mol() using atom_info " << atom_info << std::endl;
               }
               if (atom_info.formal_charge.first)
                  rdkit_at->setFormalCharge(atom_info.formal_charge.second);

               // set the valence from they type energy.  Abstract?
               //
               std::string type_energy = restraints.type_energy(at->name);
               if (type_energy != "") {
                  if (type_energy == "NT") {
                     bool charge_it = true;
                     // but don't charge it if there are 3 non-hydrogen bonds.

                     // Actually, not charging is good for layout, but is it
                     // good for pyrogen? Hmmm.
                     //
                     // A better test would be to check any of the bonds to this
                     // atom for aromaticity (in which case, don't charge)
                     //
                     bool include_H_neighb_bonds = false;
                     if (restraints.neighbours(atom_name, include_H_neighb_bonds).size() == 3)
                        charge_it = false;

                     if (charge_it)
                        rdkit_at->setFormalCharge(1);
                  }
                  
                  // other NT*s will drop hydrogens in RDKit, no need to
                  // fix up formal charge (unless there is a hydrogen! Hmm).
               
                  // now that phosphates are undelocalized, we don't
                  // need to make this hack:
                  // 
                  // if (type_energy == "P") {
                  //    rdkit_at->setFormalCharge(1);
                  // }
                  
               }

               set_atom_chirality(rdkit_at, at, residue_p, restraints);

               if (false) {
                  RDKit::Atom::ChiralType ct = rdkit_at->getChiralTag();
                  std::string cts = "!";
                  if (ct == RDKit::Atom::CHI_UNSPECIFIED)     cts = "-";
                  if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CW)  cts = " CW";
                  if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CCW) cts = "CCW";
                  if (ct == RDKit::Atom::CHI_OTHER)           cts = "Oth";
                  // std::cout << "############# After chiral set: atom name " << atom_name
                  // << " Chir: " << ct << " " << cts << std::endl;
               }

               m.addAtom(rdkit_at);
               
               if (debug)
                  std::cout << "     adding atom with name \"" << atom_name
                            << "\" to added_atom_names which is currently of size "
                            << added_atom_names.size();

               added_atom_names.push_back(atom_name);
               added_atoms.push_back(residue_atoms[iat]);
               atom_index[atom_name] = current_atom_id;
               current_atom_id++; // for next round
            }
            catch (const std::exception &rte) {
               std::cout << rte.what() << std::endl;
            }
         }
         if (debug)
            std::cout << std::endl;
      }

      if (debug) {
         std::cout << "DEBUG:: number of atoms in rdkit mol: " << m.getNumAtoms() << std::endl;
      } 

      // Doing wedge bonds before we set the chirality doesn't make sense.
      // So this code needs to be moved down.
   
      for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
         if (false)
            std::cout << "   handling bond " << ib << " of " << restraints.bond_restraint.size()
                      << " :" << restraints.bond_restraint[ib].atom_id_1_4c() << ": " 
                      << " :" << restraints.bond_restraint[ib].atom_id_2_4c() << ": " 
                      << restraints.bond_restraint[ib].type()
                      << std::endl;
         RDKit::Bond::BondType type = convert_bond_type(restraints.bond_restraint[ib].type());
         RDKit::Bond *bond = new RDKit::Bond(type);

         // if the restraints said that this bond was delocalized - we
         // want to know that when we create an mdl file for mogul.
         //
         if (restraints.bond_restraint[ib].type() == "deloc") {
            std::string prop_type("restraints-type");
            std::string bond_type("deloc");
            bond->setProp(prop_type, bond_type);
         }
      
         std::string atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
         std::string atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();
         std::string ele_1 = restraints.element(atom_name_1);
         std::string ele_2 = restraints.element(atom_name_2);
         int idx_1 = -1; // unset initially
         int idx_2 = -1; // unset

         // this block sets idx_1 and idx_2
         //
         for (unsigned int iat=0; iat<added_atom_names.size(); iat++) {
            if (added_atom_names[iat] == atom_name_1) { 
               idx_1 = iat;
               break;
            }
         }
         for (unsigned int iat=0; iat<added_atom_names.size(); iat++) { 
            if (added_atom_names[iat] == atom_name_2) { 
               idx_2 = iat;
               break;
            }
         }

         if (idx_1 != -1) { 
            if (idx_2 != -1) {

               // wedge bonds should have the chiral centre as the first atom.
               //
               bool swap_order = false;
               if (restraints.chiral_restraint.size()) {
                  swap_order = chiral_check_order_swap(m[idx_1], m[idx_2], restraints.chiral_restraint);
               } else {
                  // use the atoms rdkit chiral status
                  // 20170810: this makes sense for atoms that are marked as S/R in the pdbx_stereo_config
                  //           but not for other ones that RDKit conjures up for as yet unknown reason
                  //           e.g. CT in FAK
                  swap_order = chiral_check_order_swap(m[idx_1], m[idx_2]);
               }

               // Finally, heuristic: so that atoms with more than 1 bond are end atom if the
               // other atom has just 1 bond (this one).
               if (! swap_order) {
                  swap_order = chiral_check_order_swap_singleton(m[idx_1], m[idx_2], restraints);
               }

               if (! swap_order) {  // normal
                  bond->setBeginAtomIdx(idx_1);
                  bond->setEndAtomIdx(  idx_2);
               } else {
                  bond->setBeginAtomIdx(idx_2);
                  bond->setEndAtomIdx(  idx_1);
               } 

               if (type == RDKit::Bond::AROMATIC) {
                  bond->setIsAromatic(true);
                  m[idx_1]->setIsAromatic(true);
                  m[idx_2]->setIsAromatic(true);
               }
               // Are you here again?
               // You're here because the dictionary for this ligand was double-read
               // and added to, not replaced, the previous dictionary.  There are two
               // sets of atoms and bonds. Use have_dictionary_for_residue_type_no_dynamic_add()
               // to check before adding the restraints. Or pre-trash the current restraints
               // Or increment the read number.
               //
               //
               m.addBond(bond); // by default, this does a copy of bond. It can take the ownership

            } else {
               if (ele_2 != " H") {

                  if (atom_name_2 == " OXT" ||
                      atom_name_1 == " O1 ") {
                     // shut up about linked carbohdyrates and modified residues
                  } else { 
                     std::cout << "WARNING:: oops, bonding in rdkit_mol() "
                               << "failed to get atom index idx_2 for atom name: "
                               << atom_name_2 << " ele :" << ele_2 << ":" << std::endl;
                     if (false) {
                        std::cout << "Here's the atoms we have:\n";
                        for (unsigned int iat=0; iat<added_atom_names.size(); iat++) 
                           std::cout << std::setw(2) << iat << " :" << added_atom_names[iat] << ":\n";
                     }
                  }
                  // give up trying to construct this thing then.
                  std::string message = "Failed to get atom index for atom name \"";
                  message += atom_name_2;
                  message += "\" in residue of type ";
                  message += residue_p->GetResName();
                  message += " ";
                  message += residue_spec_t(residue_p).format();
                  throw std::runtime_error(message);
               }
            }
         } else {
            if (ele_1 != " H") { 
               if (atom_name_2 == " OXT" ||
                   atom_name_1 == " O1 ") {
                  // shut up about linked carbohdyrates and modified residues
               } else { 
                  std::cout << "WARNING:: oops, bonding in rdkit_mol() "
                            << "failed to get atom index idx_1 for atom name: \""
                            << atom_name_1 << "\" ele :" << ele_1 << ":" << std::endl;
                  if (false) {
                     std::cout << "Here's the atoms we have:\n";
                     for (unsigned int iat=0; iat<added_atom_names.size(); iat++) 
                        std::cout << std::setw(2) << iat << " :" << added_atom_names[iat] << ":";
                     // give up trying to construct this thing then. (Come back with a full dictionary molecule)
                     std::cout << std::endl;
                  }
               }
               std::string message = "Failed to get atom index for atom name \"";
               message += atom_name_1;
               message += "\"";
               message += "\" in residue of type ";
               message += residue_p->GetResName();
               message += " ";
               message += residue_spec_t(residue_p).format();
               throw std::runtime_error(message);
            }
         }
         delete bond;
      }
 
      if (debug) { 
         std::cout << "DEBUG:: rdkit_mol() number of bond restraints:    "
                   << restraints.bond_restraint.size() << std::endl;
         std::cout << "------- post construction of atoms ------" << std::endl;
         debug_rdkit_molecule(&m);
      }

      // Now all the bonds are in place.  We can now try to add an extra H to the ring N that
      // need them.  We do this after, because all the molecule bonds have to be in place for
      // us to know if the H is there already (tested in add_H_to_ring_N_as_needed()).
   
      std::vector<int> Hs_added_list; // a list of atoms to which extra Hs have been added
                                      // (this is needed, because we only want to add a
                                      // one H to an aromatic N and as we run though the
                                      // bond list, we typically find N-C and C-N bonds,
                                      // which would result in the H being added twice -
                                      // not good.

      for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) { 

         RDKit::Bond::BondType type = convert_bond_type(restraints.bond_restraint[ib].type());
         if (type == RDKit::Bond::AROMATIC) {
            std::string atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
            std::string atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();
            std::string ele_1 = restraints.element(atom_name_1);
            std::string ele_2 = restraints.element(atom_name_2);
            int idx_1 = -1; // unset
            int idx_2 = -1; // unset

            // we can't run through n_residue_atoms because the atom in the
            // mmdb::Residue may not have been added to the atom in the rdkit
            // molecule (as is the case for an alt conf).

            for (unsigned int iat=0; iat<m.getNumAtoms(); iat++) {
               try {
                  std::string name;
                  RDKit::Atom *at_p = m[iat];
                  at_p->getProp("name", name);
                  if (name == atom_name_1)
                     idx_1 = iat;
                  if (name == atom_name_2)
                     idx_2 = iat;
               }
               catch (const KeyErrorException &err) {
                  // this happens for alt conf
                  // std::cout << "caught no-name exception in rdkit_mol H-block" << std::endl;
               }
            }

            if (debug) { 
               std::cout << "idx_1 " << idx_1 << " " << n_residue_atoms << std::endl;
               std::cout << "idx_2 " << idx_2 << " " << n_residue_atoms << std::endl;
            }
      
            if (idx_1 != -1) { 
               if (idx_2 != -1) {         
   
                  // special edge case for aromatic ring N that may need an H attached for
                  // kekulization (depending on energy type).

                  if (m[idx_1]->getAtomicNum() == 7) {
                     if (std::find(Hs_added_list.begin(), Hs_added_list.end(), idx_1) == Hs_added_list.end()) { 
                        std::string n = add_H_to_ring_N_as_needed(&m, idx_1, atom_name_1, restraints);
                        if (0)
                           std::cout << "testing 1 idx_1 " << idx_1 << " idx_2 " << idx_2
                                     << " n was :" << n << ":" << std::endl;
                        if (n != "")
                           added_atom_names.push_back(n);
                        Hs_added_list.push_back(idx_1);
                     }
                  }
                  if (m[idx_2]->getAtomicNum() == 7) {
                     if (std::find(Hs_added_list.begin(), Hs_added_list.end(), idx_2) == Hs_added_list.end()) { 
                        std::string n = add_H_to_ring_N_as_needed(&m, idx_2, atom_name_2, restraints);
                        // std::cout << "testing 2 idx_1 " << idx_1 << " idx_2 " << idx_2
                        // << " n was :" << n << ":" << std::endl;
                        if (n != "")
                           added_atom_names.push_back(n);
                        Hs_added_list.push_back(idx_2);
                     }
                  }
               }
            }
         }
      }

      if (debug)
         debug_rdkit_molecule(&m);

      if (do_undelocalize) { 
         if (debug)
            std::cout << "=============== calling undelocalise() " << &m << std::endl;
         coot::undelocalise(&m);
      }

      if (debug)
         std::cout << "---------------------- calling assign_formal_charges() -----------"
                   << std::endl;
      coot::assign_formal_charges(&m); // those not in the cif file, that is

      if (debug)
         std::cout << "---------------------- getting ring info findSSSR() -----------"
                   << std::endl;
      std::vector<std::vector<int> > ring_info;
      RDKit::MolOps::findSSSR(m, ring_info);

      if (debug) {
         // what's the ring info then?
         RDKit::RingInfo* ring_info_p = m.getRingInfo();
         unsigned int n_rings = ring_info_p->numRings();
         std::cout << "found " << n_rings << " rings" << std::endl;
      }

      if (debug)
         std::cout << "---------------------- calling cleanUp() -----------" << std::endl;
      RDKit::MolOps::cleanUp(m);

      // OK, so cleanUp() doesn't fix the N charge problem our prodrg molecule
      // 
      if (false) { // debug, formal charges
         std::cout << "::::::::::::::::::::::::::: after cleanup :::::::::::::::::"
                   << std::endl;
         int n_mol_atoms = m.getNumAtoms();
         for (int iat=0; iat<n_mol_atoms; iat++) {
       RDKit::Atom* at_p = m[iat];
       std::string name = "";
            try {
               at_p->getProp("name", name);
            }
            catch (const KeyErrorException &kee) {
               std::cout << "caught no-name for atom exception in rdkit_mol(): "
                         <<  kee.what() << std::endl;
            }
            int formal_charge = at_p->getFormalCharge();
            std::cout << name << " formal_charge " << formal_charge << std::endl;
         }
         std::cout << "::::::::::: done " << std::endl;
      }

      // 2016014-PE needs investigating.
      // 
      // do we need to do this to fix up chirality?
      // rdkit_mol_sanitize(m);
   
      if (debug)
         std::cout << "DEBUG:: sanitizeMol() " << std::endl;
      RDKit::MolOps::sanitizeMol(m);

      // Now all the atoms have been added. If we try to run assignAtomCIPRanks() too early
      // we get:
      // 
      // ****
      // Pre-condition Violation
      // getNumImplicitHs() called without preceding call to calcImplicitValence()
      // Violation occurred on line 166 in file xxx/rdkit-Release_2015_03_1/Code/GraphMol/Atom.cpp
      // Failed Expression: d_implicitValence>-1
      // ****

      // 
      // Either from 3d or via R/S.  We need to calculate the CIP ranks.
      // 
      // Let's try R/S first.
      // 
      // Which means that we assign R/S according the the dictionary - and
      // ignore the positions of the atoms.
      //
      // My understanding of the CHI_CCW and CHI_CW based on atom positions is this:
      // Take the first 3 non-hydrogen atom neighbours of the chiral centre atom (CCat) in order,
      // calculate the vectors to these atoms from the chiral centre atom: p1, p2, p3
      // calculate p1.p2xp3: if it's positive then CW, if it's negative CHI_CCW.
      // This is not the problem though.  We don't have positions - we want to set the chirality
      // so that, when the postions are calculated, the chirality is correct in 3D.
      // 
      // OK, so the 3 non-hydrogen neighbours of CCat are not necessarily the top-ranked CIP atoms.
      // If one of the substituents is a H atom, then they will be - this is the easy case.
      // If (n_neighbs == 3) then if they come (in the bonds list for this atom) in the same order
      // as the CIP ranks them - or are a circular permutation of that then if R then CW
      // if S the CCW. Similar and reversed logic for non-circular permutation.
      // 
      // OK, more tricky, we have 4 non-hydrogen neighbours.  The CCW only relates to the
      // first 3 neighbours.
      //
      // First test are the first 3 neighbours the same as the highest 3 CIP ranked neighbours?
      // If yes, then same case as above with 3 neighbours.
      //
      // If not, then R for the CIP ranked neighbours will be CCW for the first 3 (when you swap
      // an atom, and the atom has a CIP rank of less than both or more than both the others, then
      // the direction changes).  If the "inserted" atom has a CIP-rank between the other two then
      // the rotation direction is preserved.
      // 
      // Likewise S will be CW for the first 3.
      // 
      // Needs testing.
      //
      // Lots of things that work:
      //
      // Fails: B7H?   C17 is not a chiral centre - but the wedgebondmol picture suggest that it is
      //               (this code doesn't think that it is).
      //

      unsigned int n_atoms = m.getNumAtoms();

      if (n_atoms > 0) {
         RDKit::UINT_VECT ranks(m.getNumAtoms(), -1);
         RDKit::Chirality::assignAtomCIPRanks(m, ranks);
         //
         for (unsigned int iat=0; iat<bonded_atoms.size(); iat++) {
            const coot::dict_atom &atom_info = restraints.atom_info[bonded_atoms[iat].second];
      
            if (atom_info.pdbx_stereo_config.first) {
               if (atom_info.pdbx_stereo_config.second == "R" ||
                   atom_info.pdbx_stereo_config.second == "S") {

                  // accumulate neigbs of the chiral atom here:
                  std::vector<std::pair<const RDKit::Atom *, unsigned int> > neighbs;
            
                  // what are the atoms bonded to this rdkit atom?
                  unsigned int idx_iat = bonded_atoms[iat].first;
                  if (idx_iat >= n_atoms) {
                     // bad!
                     std::cout << "ERROR:: rdkit_mol() chiral-check: trying to get atom with "
                               << "index  " << idx_iat << " but molecule has " << n_atoms
                               << " atoms" << std::endl;
                  } else {
                     // happy path
                     RDKit::Atom *rdkit_at = m[idx_iat];  // probably - or always?

                     RDKit::ROMol::OEDGE_ITER beg,end;
                     boost::tie(beg,end) = m.getAtomBonds(rdkit_at);
                     while(beg != end){
                        const RDKit::Bond *bond=m[*beg];
                        ++beg;
                        const RDKit::Atom *nbr=bond->getOtherAtom(rdkit_at);
                        unsigned int cip_rank = 0;
                        nbr->getProp(RDKit::common_properties::_CIPRank, cip_rank);
                        if (false)
                           std::cout << "debug:: in rdkit_mol(residue *version) iat bonded: "
                                     << iat << " cip_rank neighb: " << cip_rank << std::endl;
                        std::pair<const RDKit::Atom *, unsigned int> p(nbr, cip_rank);
                        neighbs.push_back(p);
                     }

                     if (false) {
                        std::cout << "in rdkit_mol(residue, ..) atom " << rdkit_at << " has stereoconfig "
                                  << atom_info.pdbx_stereo_config.second << " and "
                                  << neighbs.size() << " non-H neighbours " << std::endl;
                        std::cout << "---------- unsorted neighbs: " << std::endl;
                        for (unsigned int jj=0; jj<neighbs.size(); jj++) {
                           std::cout << neighbs[jj].first << " " << neighbs[jj].second << std::endl;
                        }
                     }

                     std::vector<std::pair<const RDKit::Atom *, unsigned int> > sorted_neighbs = neighbs;
                     std::sort(sorted_neighbs.begin(), sorted_neighbs.end(), cip_rank_sorter);

                     if (false) {
                        std::cout << "---------- CIP sorted neighbs: " << std::endl;
                        for (unsigned int jj=0; jj<sorted_neighbs.size(); jj++) { 
                           std::cout << jj << " " << sorted_neighbs[jj].first << " "
                                     << sorted_neighbs[jj].second << std::endl;
                        }
                     }

                     bool inverted = true; // set this using cleverness

                     if (neighbs.size() == 3) {

                        neighbs.resize(3);
                        sorted_neighbs.resize(3);

                        if (neighbs[0] == sorted_neighbs[0])
                           if (neighbs[1] == sorted_neighbs[1])
                              if (neighbs[2] == sorted_neighbs[2])
                                 inverted = false;
                     
                        if (neighbs[0] == sorted_neighbs[1])
                           if (neighbs[1] == sorted_neighbs[2])
                              if (neighbs[2] == sorted_neighbs[0])
                                 inverted = false;

                        if (neighbs[0] == sorted_neighbs[2])
                           if (neighbs[1] == sorted_neighbs[0])
                              if (neighbs[2] == sorted_neighbs[1])
                                 inverted = false;

                     } else {

                        if (neighbs.size() == 4) { // what else can it be?

                           // are the first 3 atoms of neighbour list the three atoms of highest CIP rank?
                  
                           bool atom_sets_match = false;
                           std::vector<const RDKit::Atom *> needed_atoms(3);
                           needed_atoms[0] = sorted_neighbs[1].first;
                           needed_atoms[1] = sorted_neighbs[2].first;
                           needed_atoms[2] = sorted_neighbs[3].first;

                           unsigned int n_found = 0;
                           for (unsigned int jj=0; jj<3; jj++) {
                              for (unsigned int ii=0; ii<3; ii++) {
                                 if (needed_atoms[ii] == neighbs[jj].first)
                                    n_found += 1;
                              }
                           }

                           if (n_found == 3) {

                              // as above

                              if (neighbs[0] == sorted_neighbs[0])
                                 if (neighbs[1] == sorted_neighbs[1])
                                    if (neighbs[2] == sorted_neighbs[2])
                                       inverted = false;

                              if (neighbs[0] == sorted_neighbs[1])
                                 if (neighbs[1] == sorted_neighbs[2])
                                    if (neighbs[2] == sorted_neighbs[0])
                                       inverted = false;

                              if (neighbs[0] == sorted_neighbs[2])
                                 if (neighbs[1] == sorted_neighbs[0])
                                    if (neighbs[2] == sorted_neighbs[1])
                                       inverted = false;

                              if (atom_info.pdbx_stereo_config.second == "R") {
                                 if (inverted)
                                    rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
                                 else
                                    rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
                              }
               
                              if (atom_info.pdbx_stereo_config.second == "S") {
                                 if (inverted)
                                    rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
                                 else
                                    rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
                              }
                     
                           } else {

                              // tricky case: the high CIP ranked atoms are not the first 3 neighbours of rdkit_at
                              unsigned int idx_cip_rank_lowest = 0;
                              unsigned int cip_rank_lowest = 99999;
                              for (unsigned int jj=0; jj<4; jj++) {
                                 if (neighbs[jj].second < cip_rank_lowest) {
                                    cip_rank_lowest = neighbs[jj].second;
                                    idx_cip_rank_lowest = jj;
                                 }
                              }

                              // idx_cip_rank_lowest should be something other than 0 now
                              //
                              // This part needs testing
                              //
                              if (false)
                                 std::cout << "debug idx_cip_rank_lowest " << idx_cip_rank_lowest << std::endl;
                              //
                              // these need checking
                              if (idx_cip_rank_lowest == 1)
                                 inverted = true;
                              if (idx_cip_rank_lowest == 2)
                                 inverted = false;
                              if (idx_cip_rank_lowest == 3)
                                 inverted = true;

                           }
                        } else {
                           std::cout << "WARNING:: crazy atom - too many connections " << atom_info << std::endl;
                        }

                     }

                     if (atom_info.pdbx_stereo_config.second == "R") {
                        if (inverted)
                           rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
                        else
                           rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
                     }
               
                     if (atom_info.pdbx_stereo_config.second == "S") {
                        if (inverted)
                           rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
                        else
                           rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
                     }
                  }
               }
            }
         }

         if (debug)
            std::cout << "DEBUG:: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "
                      << "in constructing rdkit molecule, now adding a conf " 
                      << "number of atoms comparison added_atom names size: " 
                      << added_atom_names.size() << " vs m.getNumAtoms() " 
                      << m.getNumAtoms() << std::endl;

         std::vector<std::string> alt_confs_in_residue = get_alt_confs_in_residue(residue_p);

         for (unsigned int iconf=0; iconf<alt_confs_in_residue.size(); iconf++) {

            std::string alt_conf = alt_confs_in_residue[iconf];

            RDKit::Conformer *conf = new RDKit::Conformer(m.getNumAtoms());
            conf->set3D(true);

            // Add positions to the conformer (only the first instance of an
            // atom with a particular atom name).
            //
            for (int iat=0; iat<n_residue_atoms; iat++) {
               std::string atom_name(residue_atoms[iat]->name);
               std::string atom_alt_conf(residue_atoms[iat]->altLoc);
               if (true) { // was alt-conf test
                  std::map<std::string, int>::const_iterator it = atom_index.find(atom_name);
                  if (it != atom_index.end()) {
                     RDGeom::Point3D pos(residue_atoms[iat]->x,
                                         residue_atoms[iat]->y,
                                         residue_atoms[iat]->z);
                     conf->setAtomPos(it->second, pos);
                     if (debug)
                        std::cout << "in construction of rdkit mol: making a conformer atom "
                                  << iat << " " << it->second << " " << atom_name << " at pos "
                                  << pos << std::endl;
                  }
               }
            }

            int conf_id = m.addConformer(conf);

            // RDKit::MolOps::assignStereochemistry(m, false, true, true); // this does not assign 
            // stereochemistry on m
   

            // 20161013 do we need to this these days?  It adds chirality to the SD in 2ZC
            // Needs more consideration.
            //
            RDKit::MolOps::assignChiralTypesFrom3D(m, conf_id, true);

         }
   
         if (debug) 
            std::cout << "ending construction of rdkit mol: n_atoms " << m.getNumAtoms()
                      << std::endl;

         set_energy_lib_atom_types(&m);

         // debugging
         // RDKit::MolToMolFile(m, "rdkit.mol");
      
      } // number of atoms test

   } // number of bonded atoms test
   
   return m;
}

// Fill a conformer with the 2D depiction if you can with the first at the conformer id, else first is -1
std::pair<int, RDKit::RWMol>
coot::rdkit_mol_with_2d_depiction(const dictionary_residue_restraints_t &restraints) {

   RDKit::RWMol mol = rdkit_mol(restraints);
   if (restraints.depiction.empty()) {
      return std::make_pair(-1, mol);
   } else {
      // happy path
      RDKit::MolOps::removeHs(mol);
      // Do we want a call to "PrepreMolForDrawing()" here?
      std::cout << "atom number compare " << mol.getNumAtoms() << " " << restraints.depiction.atoms.size()
                << std::endl;
      if (mol.getNumAtoms() == restraints.depiction.atoms.size()) {
         RDKit::Conformer *conf = new RDKit::Conformer(mol.getNumAtoms());
         conf->set3D(false);
         for (unsigned int i=0; i<restraints.depiction.atoms.size(); i++) {
            const auto &atom(restraints.depiction.atoms[i]);
            RDGeom::Point3D pos(atom.x, atom.y, 0.0);
            conf->setAtomPos(i, pos);
         }
         int iconf = mol.addConformer(conf, true);
         std::cout << "debug:: Happy return iconf: " << iconf << std::endl;
         return std::make_pair(iconf, mol);
      } else {
         return std::make_pair(-1, mol);
      }
   }
}


void
coot::set_atom_chirality(RDKit::Atom *rdkit_at,
                         mmdb::Atom *at,
                         mmdb::Residue *residue_p,
                         const coot::dictionary_residue_restraints_t &restraints) {

   // set the chirality
   // (if this atom is has restraints-style chiral info)
   //
   bool done_chiral = false;

   std::string atom_name = at->name;

   for (unsigned int ichi=0; ichi<restraints.chiral_restraint.size(); ichi++) {
      const dict_chiral_restraint_t &cr = restraints.chiral_restraint[ichi];
      if (cr.atom_id_c_4c() == atom_name) {
         done_chiral = true;
         if (!cr.has_unassigned_chiral_volume()) {
            rdkit_at->setProp("mmcif_chiral_N1", util::remove_whitespace(cr.atom_id_1_4c()));
            rdkit_at->setProp("mmcif_chiral_N2", util::remove_whitespace(cr.atom_id_2_4c()));
            rdkit_at->setProp("mmcif_chiral_N3", util::remove_whitespace(cr.atom_id_3_4c()));
            if (!cr.is_a_both_restraint()) {
               // e.g. RDKit::Atom::CHI_TETRAHEDRAL_CCW;
               RDKit::Atom::ChiralType chiral_tag = get_chiral_tag(residue_p, restraints, at);
               rdkit_at->setChiralTag(chiral_tag);

               std::string bc("positive");
               if (cr.volume_sign == dict_chiral_restraint_t::CHIRAL_RESTRAINT_NEGATIVE)
                  bc = "negative";

               rdkit_at->setProp("mmcif_chiral_volume_sign", bc);

            } else {
               std::string bc("both");
               rdkit_at->setProp("mmcif_chiral_volume_sign", bc);
            } 
         }
      }
   }

   // set chirality
   // (if this atom has Chemical Component Dictionary style chirality (R/S pdbx_stereo_config_flag)
   //
   if (! done_chiral) {
      for (unsigned int i=0; i<restraints.atom_info.size(); i++) { 
         if (restraints.atom_info[i].atom_id_4c == atom_name) {
            set_atom_chirality(rdkit_at, restraints.atom_info[i]);
         }
      }
   }
}

void
coot::set_atom_chirality(RDKit::Atom *rdkit_at, const coot::dict_atom &dict_atom) {

   bool debug = false;

   if (dict_atom.pdbx_stereo_config.first) {
      if (dict_atom.pdbx_stereo_config.second == "R") {

         // "work it out later" using rdkit sanitize doesn't seem to work
         //
         // RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_UNSPECIFIED;
         RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;

         if (debug)
            std::cout << "   pdbx_stereo_config: " << dict_atom.atom_id << " R -> CW " << std::endl;
         rdkit_at->setChiralTag(chiral_tag);
         std::string cip = "R";
         rdkit_at->setProp("_CIPCode", cip);
      }
      if (dict_atom.pdbx_stereo_config.second == "S") {
         RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
         std::string cip = "S";
         rdkit_at->setProp("_CIPCode", cip);
         rdkit_at->setChiralTag(chiral_tag);
         if (debug)
            std::cout << "   pdbx_stereo_config: " << dict_atom.atom_id << " S -> CCW " << std::endl;
      }
      if (dict_atom.pdbx_stereo_config.second == "N") {
         if (false) // otherwise too noisy
            std::cout << "pdbx_stereo_config says N for " << dict_atom.atom_id << std::endl;
      }
   } else {
      if (false)
         std::cout << "No pdbx_stereoconfig for atom " << dict_atom.atom_id << std::endl;
   }
}

// sorts atoms so that the smallest ranks are at the top (close to index 0).
//
bool
coot::cip_rank_sorter(const std::pair<const RDKit::Atom *, unsigned int> &at_1,
                      const std::pair<const RDKit::Atom *, unsigned int> &at_2) {

   return (at_1.second < at_2.second);
}


#if (RDKIT_VERSION >= RDKIT_VERSION_CHECK(2018, 3, 1))
bool
coot::chiral_check_order_swap(RDKit::Atom* at_1, RDKit::Atom* at_2,
                              const std::vector<dict_chiral_restraint_t>  &chiral_restraints) {
#else
bool
coot::chiral_check_order_swap(RDKit::Atom *at_1, RDKit::Atom *at_2,
                              const std::vector<dict_chiral_restraint_t>  &chiral_restraints) {
#endif
   bool status = false;

   // set status to true (only) if the first atom is not chiral and the
   // second one is.
   
   try {
      bool second_is_chiral = false;
      std::string name_1;
      std::string name_2;
      at_1->getProp("name", name_1);
      at_2->getProp("name", name_2);

      for (unsigned int ich=0; ich<chiral_restraints.size(); ich++) {
         if (chiral_restraints[ich].atom_id_c_4c() == name_2) {
            second_is_chiral = true;
            break;
         }
      }

      if (second_is_chiral) {
         bool first_is_chiral = false;
         for (unsigned int ich=0; ich<chiral_restraints.size(); ich++) {
            if (chiral_restraints[ich].atom_id_c_4c() == name_1) {
               first_is_chiral = true;
               break;
            }
         }

         if (!first_is_chiral)
            status = true;
      }
   }
   catch (...) {
      // this should not catch anything, the names should be set as properties
   }

   return status;
}

#if (RDKIT_VERSION >= RDKIT_VERSION_CHECK(2018, 3, 1))
bool
coot::chiral_check_order_swap(RDKit::Atom* at_1, RDKit::Atom* at_2) {
#else
bool
coot::chiral_check_order_swap(RDKit::Atom *at_1, RDKit::Atom *at_2) {
#endif
   bool status = false;

   RDKit::Atom::ChiralType chiral_tag_1 = at_1->getChiralTag();
   if ( (chiral_tag_1 != RDKit::Atom::CHI_TETRAHEDRAL_CW) &&
        (chiral_tag_1 != RDKit::Atom::CHI_TETRAHEDRAL_CCW)) {
      RDKit::Atom::ChiralType chiral_tag_2 = at_2->getChiralTag();
      if (chiral_tag_2 == RDKit::Atom::CHI_TETRAHEDRAL_CW)
         status = true;
      if (chiral_tag_2 == RDKit::Atom::CHI_TETRAHEDRAL_CCW)
         status = true;
   }

   return status;
}

bool
coot::chiral_check_order_swap_singleton(RDKit::Atom *at_1, RDKit::Atom *at_2,
                                        const dictionary_residue_restraints_t  &restraints) {

   // This improves things, but needs further improvement because it doesn't correct the
   // atom order of a bond of a CH3 connected to a C that is connected to 3 other CH3s -
   // e.g. GN5
   //

   bool status = false;

   try {
      std::string name_1;
      std::string name_2;
      at_1->getProp("name", name_1);
      at_2->getProp("name", name_2);

      bool allow_H = true;
      std::vector<std::string> n_1 = restraints.neighbours(name_1, allow_H);
      std::vector<std::string> n_2 = restraints.neighbours(name_2, allow_H);
      if (n_1.size() == 1)
         if (n_2.size() > 1)
            status = true;

      // perhaps it was an OH?
      if (! status) {
         allow_H = false;
         n_1 = restraints.neighbours(name_1, allow_H);
         n_2 = restraints.neighbours(name_2, allow_H);
         if (n_1.size() == 1)
            if (n_2.size() > 1)
               status = true;

         // actually, let's turn around everything where the second atom has more
         // non-H bonds than the first
         if (n_2.size() > n_1.size())
            status = true;
      }
   }
   catch (...) {
      // this should not catch anything, the names should be set as properties
   }
   return status;
}



// fill the coords from the dictionary if you can.
// 
RDKit::RWMol
coot::rdkit_mol(const coot::dictionary_residue_restraints_t &r) {

   bool debug = false;
   RDKit::RWMol m;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

   std::map<std::string, int> added_atoms; // atom name to rdkit atom index

   // ------------------------------------ Atoms -----------------------------

   for (unsigned int iat=0; iat<r.atom_info.size(); iat++) { 

      try {
         RDKit::Atom *at = new RDKit::Atom;
         std::string atom_name = r.atom_info[iat].atom_id_4c;
         std::string ele_capped =
            coot::util::capitalise(coot::util::remove_leading_spaces(r.atom_info[iat].type_symbol));
         int atomic_number = tbl->getAtomicNumber(ele_capped);
         at->setAtomicNum(atomic_number);
         // at->setMass(tbl->getAtomicWeight(atomic_number));
         at->setIsotope(0);
         at->setProp("name", atom_name);

         if (false)
            std::cout << "debug:: in rdkit_mol(restraints) iat " << iat << " name: " << atom_name
                      << " formal charge set "   << r.atom_info[iat].formal_charge.first
                      << " formal charge value " << r.atom_info[iat].formal_charge.second << std::endl;

         // formal charge
         if (r.atom_info[iat].formal_charge.first)
            at->setFormalCharge(r.atom_info[iat].formal_charge.second);

         // set the chirality (if this atom is a chiral centre of a chiral restraint).
         //
         bool done_chiral = false;
         for (unsigned int ichi=0; ichi<r.chiral_restraint.size(); ichi++) {
            if (r.chiral_restraint[ichi].atom_id_c_4c() == r.atom_info[iat].atom_id_4c) {
               if (!r.chiral_restraint[ichi].has_unassigned_chiral_volume()) {
                  if (!r.chiral_restraint[ichi].is_a_both_restraint()) {
                     RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
                     if (r.chiral_restraint[ichi].volume_sign ==
                         dict_chiral_restraint_t::CHIRAL_RESTRAINT_NEGATIVE)
                        chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
                     // std::cout << ".... seting chiral tag to " << chiral_tag << std::endl;
                     at->setChiralTag(chiral_tag);
                     done_chiral = true;
                  }
               }
            }
         }

         // need to try to get chiral info using atom_info[iat].pdbx_stereo_config
         //
         if (! done_chiral)
            set_atom_chirality(at, r.atom_info[iat]);

         if (false) {
            RDKit::Atom::ChiralType ct = at->getChiralTag();
            std::string cts = "!";
            if (ct == RDKit::Atom::CHI_UNSPECIFIED)     cts = "-";
            if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CW)  cts = " CW";
            if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CCW) cts = "CCW";
            if (ct == RDKit::Atom::CHI_OTHER)           cts = "Oth";
            // std::cout << "############# After chiral set: atom name " << atom_name
            // << " Chir: " << ct << " " << cts << std::endl;
         }
         
         int idx = m.addAtom(at);
         std::string key = coot::util::remove_whitespace(r.atom_info[iat].atom_id_4c);
         added_atoms[key] = idx; // for making bonds.
      }
      catch (const std::exception &rte) {
         std::cout << rte.what() << std::endl;
      }
   }

   // 20160702 add atoms to a conformer
   RDKit::Conformer *conf = new RDKit::Conformer(m.getNumAtoms());
   conf->set3D(true);
   for (unsigned int iat=0; iat<r.atom_info.size(); iat++) {
      if (false)
         std::cout << "rdkit_mol(): atom info loop iat " << iat << " "
                   << r.atom_info[iat].pdbx_model_Cartn_ideal.first
                   << " " << r.atom_info[iat].model_Cartn.first << std::endl;
      try {
         if (r.atom_info[iat].pdbx_model_Cartn_ideal.first) {
            RDGeom::Point3D pos(r.atom_info[iat].pdbx_model_Cartn_ideal.second.x(),
                                r.atom_info[iat].pdbx_model_Cartn_ideal.second.y(),
                                r.atom_info[iat].pdbx_model_Cartn_ideal.second.z());
            conf->setAtomPos(iat, pos);
         } else {
            if (r.atom_info[iat].model_Cartn.first) {
               RDGeom::Point3D pos(r.atom_info[iat].model_Cartn.second.x(),
                                   r.atom_info[iat].model_Cartn.second.y(),
                                   r.atom_info[iat].model_Cartn.second.z());
               conf->setAtomPos(iat, pos);
            }
         }
      }
      catch (const std::exception &rte) {
         std::cout << rte.what() << std::endl;
      }
   }
   m.addConformer(conf);

   
   // ------------------------------------ Bonds -----------------------------

   if (false) { // for checking spaces in atom names
      // the check in the below loop uses atom names with whitespace stripped.
      for (const auto &item : added_atoms) {
         std::cout << "debug added_atom \"" << item.first << "\" " << item.second << std::endl;
      }
   }

   int n_atoms = m.getNumAtoms();
   std::map<std::string, int>::const_iterator it_1;
   std::map<std::string, int>::const_iterator it_2;
   int idx_1, idx_2;
   for (unsigned int ib=0; ib<r.bond_restraint.size(); ib++) {
      const dict_bond_restraint_t &br = r.bond_restraint[ib];
      std::string at_1 = coot::util::remove_whitespace(br.atom_id_1());
      std::string at_2 = coot::util::remove_whitespace(br.atom_id_2());
      it_1 = added_atoms.find(at_1);
      it_2 = added_atoms.find(at_2);
      if (it_1 != added_atoms.end()) { 
         if (it_2 != added_atoms.end()) {
            idx_1 = it_1->second;
            idx_2 = it_2->second;
            RDKit::Bond::BondType type = convert_bond_type(br.type());
            RDKit::Bond *bond = new RDKit::Bond(type);

            if (idx_1 < n_atoms) {
               if (idx_2 < n_atoms) {
            
                  // wedge bonds should have the chiral centre as the first atom.
                  //
                  bool swap_order = false;

                  if (r.chiral_restraint.size()) {
                     swap_order = chiral_check_order_swap(m[idx_1], m[idx_2], r.chiral_restraint);
                  } else {
                     // use the atoms rdkit chiral status
                     swap_order = chiral_check_order_swap(m[idx_1], m[idx_2]);
                  }
                  if (! swap_order) { // normal
                     bond->setBeginAtomIdx(idx_1);
                     bond->setEndAtomIdx(  idx_2);
                  } else {
                     bond->setBeginAtomIdx(idx_2);
                     bond->setEndAtomIdx(  idx_1);
                  } 
            
                  if (type == RDKit::Bond::AROMATIC) {
                     bond->setIsAromatic(true);
                     m[idx_1]->setIsAromatic(true);
                     m[idx_2]->setIsAromatic(true);
                  }
                  m.addBond(bond); // does a copy (it can take ownership with extra arg)
               } else {
                  std::cout << "ERROR:: atom indexing problem " << idx_2 << " " << n_atoms << std::endl;
               }
            } else {
               std::cout << "ERROR:: atom indexing problem " << idx_1 << " " << n_atoms << std::endl;
            }
            delete bond;
         }
      }
   }

   if (debug)
      std::cout << "DEBUG:: rdkit_mol(): numbonds " << m.getNumBonds() << std::endl;

   set_3d_conformer_state(&m);

   if (debug)
      std::cout << "---------------------- calling assign_formal_charges() -----------"
                << std::endl;
   coot::assign_formal_charges(&m);

   if (debug)
      std::cout << "---------------------- getting ring info findSSSR() -----------"
                << std::endl;

   std::vector<std::vector<int> > ring_info;
   RDKit::MolOps::findSSSR(m, ring_info);

   if (debug) {
      // what's the ring info then?
      RDKit::RingInfo* ring_info_p = m.getRingInfo();
      unsigned int n_rings = ring_info_p->numRings();
      std::cout << "INFO:: ring-info: found " << n_rings << " rings" << std::endl;
   }

   if (debug)
      std::cout << "---------------------- calling cleanUp() -----------" << std::endl;
   RDKit::MolOps::cleanUp(m);
   if (debug)
      std::cout << "---------------------- calling sanitizeMol() -----------" << std::endl;
   RDKit::MolOps::sanitizeMol(m); // doesn't seem to do chirality assignement
                                  // if chiral centres are set to CHI_UNSPECIFIED
                                  // (I thought that it should - needs more digging)
                                  // Now chiral centres are set translating
                                  // pdbx_stereo_config R to CW and S to CCW.
                                  // which presumes that the pdbx CIP codes are the
                                  // same as RDKit's.

   set_energy_lib_atom_types(&m); // Refmac types used for H-bonding
   return m;
}



// can throw a std::runtime_error or std::exception.
//
// should kekulize flag be an argmuent?
// 
RDKit::RWMol
coot::rdkit_mol_sanitized(mmdb::Residue *residue_p, int imol_enc, const protein_geometry &geom) {

   RDKit::RWMol mol = rdkit_mol(residue_p, imol_enc, geom);
   rdkit_mol_sanitize(mol);
   return mol;
}

void
coot::rdkit_mol_sanitize(RDKit::RWMol &mol) {
   
   // clear out any cached properties
   mol.clearComputedProps();
   // clean up things like nitro groups
   RDKit::MolOps::cleanUp(mol);
   // update computed properties on atoms and bonds:
   mol.updatePropertyCache();
   RDKit::MolOps::Kekulize(mol);
   RDKit::MolOps::assignRadicals(mol);
            
   // then do aromaticity perception
   RDKit::MolOps::setAromaticity(mol);

   // std::cout << "-------- coot::rdkit_mol_sanitize() mol: " << std::endl;
   // debug_rdkit_molecule(&mol);
    
   // set conjugation
   RDKit::MolOps::setConjugation(mol);
               
   // set hybridization
   RDKit::MolOps::setHybridization(mol); // non-linear ester bonds.

   // remove bogus chirality specs:
   RDKit::MolOps::cleanupChirality(mol);

}

// hack the setting of 3D state, seems not to
// be done for mdl files when zs are 0.
void
coot::set_3d_conformer_state(RDKit::RWMol *mol) {

   if (mol) {
      for (unsigned int iconf=0; iconf<mol->getNumConformers(); iconf++) { 
         RDKit::Conformer &conf = mol->getConformer(iconf);
         int n_atoms = conf.getNumAtoms();
         bool all_zero_z = true;
         for (int iat=0; iat<n_atoms; iat++) { 
            RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
            if ((r_pos.z < -0.01) || (r_pos.z > 0.01)) {
               all_zero_z = false;
               break;
            }
         }
         if (all_zero_z) {
            // std::cout << "conformer " << iconf << " set 3d false" << std::endl;
            conf.set3D(false);
         } else {
            // std::cout << "conformer " << iconf << " set 3d true" << std::endl;
            conf.set3D(true);
         }
         // std::cout << "conformer " << iconf << " is3D(): " << conf.is3D() << std::endl;
      }
   } else {
      std::cout << "WARNING:: in set_3d_conformer_state() null mol " << std::endl;
   }
}

// e.g. reading from a MolFile
bool
coot::has_zero_coords(RDKit::RWMol *mol, unsigned int iconf) {

   bool zero_z = true;
   
   if (mol) {
      if (iconf<mol->getNumConformers()) {
         const RDKit::Conformer &conf = mol->getConformer(iconf);
         int n_atoms = conf.getNumAtoms();
         for (int iat=0; iat<n_atoms; iat++) { 
            const RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
            if (r_pos.lengthSq() > 0.1) {
               zero_z = false;
               break;
            }
         }
      }
   }
   return zero_z;
}




RDKit::Bond::BondType
coot::convert_bond_type(const std::string &t) {

   // It was a mistake to store the bond order as a string, wasn't it?
   
   RDKit::Bond::BondType bt = RDKit::Bond::UNSPECIFIED;
   if (t == "single")
      bt = RDKit::Bond::SINGLE;
   if (t == "double")
      bt = RDKit::Bond::DOUBLE;
   if (t == "triple")
      bt = RDKit::Bond::TRIPLE;
   if (t == "coval")
      bt = RDKit::Bond::SINGLE;
   if (t == "deloc")
      bt = RDKit::Bond::ONEANDAHALF;
   if (t == "aromatic")
      bt = RDKit::Bond::AROMATIC;
   if (t == "arom")
      bt = RDKit::Bond::AROMATIC;
   if (t == "aromat")
      bt = RDKit::Bond::AROMATIC;
   
   if (0) { // debug
      std::cout << "created RDKit bond type " << bt;
      if (bt == RDKit::Bond::AROMATIC)
         std::cout << " (aromatic)";
      std::cout << std::endl;
   }
   
   return bt;
}

// used in the rdkit_mol() "constructor", e.g. in thumbnails
// 
RDKit::Atom::ChiralType
coot::get_chiral_tag(mmdb::Residue *residue_p,
                     const dictionary_residue_restraints_t &restraints,
                     mmdb::Atom *atom_p) {

   RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_UNSPECIFIED; // as yet

   if (! residue_p) return chiral_tag;
   
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   std::string atom_name = atom_p->name;
   bool debug = false;

   // To make RDKit/SMILES chiral tags, we consider the order of the 3
   // atoms in the atom list that are after the chiral centre.
   // 
   // But, the order of atoms from a residue in a PDB is can be
   // arbitrary, so maybe chiral atom appears after all of the atoms
   // to which it is bonded... What to do?  Let's take the last
   // 3. That makes sense, the first-appearing atom is the "from" atom
   // in SMILES encoding of chiralilty.
   
   // does the order of the restraints match the order of the atoms?
   //
   for (unsigned int ichi=0; ichi<restraints.chiral_restraint.size(); ichi++) { 
      if (restraints.chiral_restraint[ichi].atom_id_c_4c() == atom_name) {
         const coot::dict_chiral_restraint_t &chiral_restraint = restraints.chiral_restraint[ichi];

         int n_neigbours_found = 0;
         std::vector<int> ni(4, -1); // neighbour_indices: gap for c atom at 0.
         bool atom_orders_match = false;

         for (int iat=0; iat<n_residue_atoms; iat++) {
            std::string atom_name_local = residue_atoms[iat]->name;
            if (atom_name_local == chiral_restraint.atom_id_1_4c()) {
               ni[1] = iat;
               n_neigbours_found++;
            }
            if (atom_name_local == chiral_restraint.atom_id_2_4c()) {
               ni[2] = iat;
               n_neigbours_found++;
            }
            if (atom_name_local == chiral_restraint.atom_id_3_4c()) {
               ni[3] = iat;
               n_neigbours_found++;
            }
         }

         if (n_neigbours_found == 3) {

            // 3 2 1
            if ((ni[3] > ni[2]) && (ni[2] > ni[1])) { 
               atom_orders_match = true;
               // std::cout << "match by method A " << std::endl;
            } 
            // circular permutation, 1 3 2 
            if ((ni[1] > ni[3]) && (ni[3] > ni[2])) { 
               atom_orders_match = true;
               // std::cout << "match by method B " << std::endl;
            } 
            // circular permutation, 2 1 3
            if ((ni[2] > ni[1]) && (ni[1] > ni[3])) {
               // std::cout << "match by method C " << std::endl;
               atom_orders_match = true;
            } 
            
            // This bit needs checking
            // 
            if (atom_orders_match) {
               if (debug) 
                  std::cout << "atom orders match:     true, vol sign "
                            << std::setw(2) << chiral_restraint.volume_sign
                            << " -> CCW "
                            << "for atom \"" << atom_name 
                            << "\" neighbs "
                            << std::setw(2) << ni[1] << " "
                            << std::setw(2) << ni[2] << " "
                            << std::setw(2) << ni[3] << " \""
                            << chiral_restraint.atom_id_1_4c() << "\" \""
                            << chiral_restraint.atom_id_2_4c() << "\" \""
                            << chiral_restraint.atom_id_3_4c() << "\"" << std::endl;
               if (chiral_restraint.volume_sign == 1)
                  chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
               if (chiral_restraint.volume_sign == -1)
                  chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
            } else {
               if (debug) 
                  std::cout << "atom orders NOT match: true, vol sign "
                            << std::setw(2) << chiral_restraint.volume_sign
                            << " ->  CW "
                            << "for atom \"" << atom_name 
                            << "\" neighbs "
                            << std::setw(2) << ni[1] << " "
                            << std::setw(2) << ni[2] << " "
                            << std::setw(2) << ni[3] << " \""
                            << chiral_restraint.atom_id_1_4c() << "\" \""
                            << chiral_restraint.atom_id_2_4c() << "\" \"" 
                            << chiral_restraint.atom_id_3_4c() << "\"" << std::endl;
               if (chiral_restraint.volume_sign == 1)
                  chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
               if (chiral_restraint.volume_sign == -1)
                  chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
            }
         }
         break; // because we found the restraint that matches the passed atom
      }
   }
   // CHI_UNSPECIFIED:     0
   // CHI_TETRAHEDRAL_CW:  1
   // CHI_TETRAHEDRAL_CCW: 2
   // CHI_OTHER:           3
   // std::cout << "returning chiral_tag " << chiral_tag << std::endl;
   return chiral_tag;
}


// used in the rdkit_mol() "constructor", e.g. in thumbnails
// 
RDKit::Atom::ChiralType
coot::get_chiral_tag_v2(mmdb::Residue *residue_p,
                        const dictionary_residue_restraints_t &restraints,
                        mmdb::Atom *atom_p) {

   RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_UNSPECIFIED; // as yet
   
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   std::string atom_name = atom_p->name;

   std::cout << "Called get_chiral_tag_v2() whti atom name " << atom_name << std::endl;

   // do we have 4 atoms bonded to the chiral centre atom?

   for (unsigned int ich=0; ich<restraints.chiral_restraint.size(); ich++) { 
      const dict_chiral_restraint_t &cr = restraints.chiral_restraint[ich];

      if (cr.atom_id_c_4c() == atom_name) { 
         std::vector<chiral_neighbour_info_t> neighbs;

         // loop over residue atoms to find chiral centre atom
         //
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               mmdb::Atom *chiral_atom = 0;
               std::string atom_name_local(at->name);
               if (atom_name_local == cr.atom_id_c_4c()) {
                  chiral_atom = at;
               }

               if (atom_name_local == cr.atom_id_1_4c()) neighbs.push_back(chiral_neighbour_info_t(at, iat, 1));
               if (atom_name_local == cr.atom_id_2_4c()) neighbs.push_back(chiral_neighbour_info_t(at, iat, 2));
               if (atom_name_local == cr.atom_id_3_4c()) neighbs.push_back(chiral_neighbour_info_t(at, iat, 3));
            }
         }

         if (neighbs.size() != 3) {
            std::cout << "Errg.  Not all chiral neighbours found " << neighbs.size() << " "
                      << atom_name << std::endl;
         } else {
            // now we get the 4th atom by looking at the atoms bonded to atom_p
            // (the we don't aleady have)

            for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
               std::string other_atom;
               const dict_bond_restraint_t &br = restraints.bond_restraint[ib];
               if (br.atom_id_1_4c() == atom_name)
                  other_atom = br.atom_id_2_4c();
               if (br.atom_id_2_4c() == atom_name)
                  other_atom = br.atom_id_1_4c();
            
               if (! other_atom.empty()) {
                  for (int iat=0; iat<n_residue_atoms; iat++) { 
                     mmdb::Atom *at = residue_atoms[iat];
                     std::string atom_name_local(at->name);
                     if (0) 
                        std::cout << iat << " comparing :" << atom_name_local << ": :"
                                  << other_atom << ":" << std::endl;
                     
                     if (atom_name_local == other_atom) {

                        // is at in the neighbs already?
                        std::vector<chiral_neighbour_info_t>::const_iterator nit;
                        bool found = false;
                        for (nit=neighbs.begin(); nit!=neighbs.end(); nit++) {
                           if (nit->at == at) {
                              found = true;
                              break;
                           }
                        }

                        if (!found) {
                           std::cout << atom_name_local << " was not found in neighbs vec" << std::endl;
                           chiral_neighbour_info_t cni(at, iat, 0);
                           neighbs.push_back(cni);
                           std::cout << "neighbs now of size() " << neighbs.size() << std::endl;
                           break;
                        } else {
                           std::cout << atom_name_local << " was already in in neighbs vec" << std::endl;
                        } 
                     }
                  }
               }
            }

            if (neighbs.size() != 4) {
               std::cout << "WARNING:: Errgh.  Not we don't have 4 chiral-centre neighbours "
                         << neighbs.size() << std::endl;
            } else {

               
               std::sort(neighbs.begin(), neighbs.end(), chiral_neighbour_info_t::neighbour_sorter);
               std::vector<chiral_neighbour_info_t> back_neighbs;

               back_neighbs.push_back(neighbs[1]);
               back_neighbs.push_back(neighbs[2]);
               back_neighbs.push_back(neighbs[3]);

               std::cout << "back_neighbs: (sorted) "
                         << back_neighbs[0].idx_mmcif << " "
                         << back_neighbs[1].idx_mmcif << " "
                         << back_neighbs[2].idx_mmcif << " "
                         << std::endl;
               std::cout << "back_neighbs:          "
                         << back_neighbs[0].idx_atom_list << " "
                         << back_neighbs[1].idx_atom_list << " "
                         << back_neighbs[2].idx_atom_list << " "
                         << std::endl;

               //                // 3 2 1
               //                if ((ni[3] > ni[2]) && (ni[2] > ni[1])) { 
               //                   atom_orders_match = true;
               //                   // std::cout << "match by method A " << std::endl;
               //                } 
               //                // circular permutation, 1 3 2 
               //                if ((ni[1] > ni[3]) && (ni[3] > ni[2])) { 
               //                   atom_orders_match = true;
               //                   // std::cout << "match by method B " << std::endl;
               //                } 
               //                // circular permutation, 2 1 3
               //                if ((ni[2] > ni[1]) && (ni[1] > ni[3])) {
               //                   // std::cout << "match by method C " << std::endl;
               //                   atom_orders_match = true;
               //                } 
               
               
               bool atom_orders_match = false;
               // 2 1 0
               if ((back_neighbs[2].idx_atom_list > back_neighbs[1].idx_atom_list) &&
                   (back_neighbs[1].idx_atom_list > back_neighbs[0].idx_atom_list)) {
                  atom_orders_match = true;
               }
               // 0 2 1 
               if ((back_neighbs[0].idx_atom_list > back_neighbs[2].idx_atom_list) &&
                   (back_neighbs[2].idx_atom_list > back_neighbs[1].idx_atom_list)) {
                  atom_orders_match = true;
               }
               // 1 0 2
               if ((back_neighbs[1].idx_atom_list > back_neighbs[0].idx_atom_list) &&
                   (back_neighbs[0].idx_atom_list > back_neighbs[2].idx_atom_list)) {
                  atom_orders_match = true;
               }

               if (atom_orders_match) {
                  if (cr.volume_sign == 1)
                     chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
                  else 
                     chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
               } else {
                  if (cr.volume_sign == -1)
                     chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
                  else 
                     chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
               } 
            }
         }
      } 
   }

   return chiral_tag;

} 

// static
bool
coot::chiral_neighbour_info_t::neighbour_sorter(const coot::chiral_neighbour_info_t &v1,
                                                const coot::chiral_neighbour_info_t &v2) {

   return (v1.idx_mmcif < v2.idx_mmcif);
}



// tweaking function used by rdkit mol construction function.
// (change mol maybe).
//
// Don't try to add an H if there is already an H on this atom (typically, this ligand had
// hydrogens (everywhere) already).
//
// Try to find the name of the Hydrogen from the bond restraints.
// If not found, add an atom called "-".
// 
// @return the added hydrogen name - or "" if nothing was added.
// 
std::string
coot::add_H_to_ring_N_as_needed(RDKit::RWMol *mol,
                                int idx, const std::string &atom_name,
                                const coot::dictionary_residue_restraints_t &restraints) {


   std::string r = "";
   std::string energy_type = restraints.type_energy(atom_name);
   
   if ((energy_type == "NR15") || (energy_type == "NR16")) {

      // first, is there an H bonded to this atom already?
      bool already_there = 0;
      unsigned int n_bonds = mol->getNumBonds();
      for (unsigned int ib=0; ib<n_bonds; ib++) {
         const RDKit::Bond *bond_p = mol->getBondWithIdx(ib);
         int idx_1 = bond_p->getBeginAtomIdx();
         int idx_2 = bond_p->getEndAtomIdx();
         if (idx_1 == idx)
            if ((*mol)[idx_2]->getAtomicNum() == 1) {
               already_there = 1;
               break;
            }
         if (idx_2 == idx)
            if ((*mol)[idx_1]->getAtomicNum() == 1) {
               already_there = 1;
               break;
            }
      }

      if (! already_there) { 
      
         // -------------  add an H atom --------------------------

         //
         // Probably is better if we do as Greg Landrum suggests: just
         // N_at->setNumExplicitHs(1), then we don't need to make an H
         // and a bond.
         // 
         RDKit::Atom *at = new RDKit::Atom;
         at->setAtomicNum(1);
         int idx_for_H = mol->addAtom(at);

         // std::string name_H = "-";
         // what is the Name of this H?

         std::string name_H = "-";
         for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) { 
            if (restraints.bond_restraint[ib].atom_id_1_4c() == atom_name) {
               if (restraints.element(restraints.bond_restraint[ib].atom_id_2_4c()) == " H") { 
                  name_H = restraints.bond_restraint[ib].atom_id_2_4c();
                  break;
               }
            } 
            if (restraints.bond_restraint[ib].atom_id_2_4c() == atom_name) {
               if (restraints.element(restraints.bond_restraint[ib].atom_id_1_4c()) == " H") { 
                  name_H = restraints.bond_restraint[ib].atom_id_1_4c();
                  break;
               }
            }
         }
      
      
         if (name_H != "-") {
            at->setProp("name", name_H);
         }
         r = name_H;

         // -------------  now add a bond --------------------------

         RDKit::Bond *bond = new RDKit::Bond(RDKit::Bond::SINGLE);
         bond->setBeginAtomIdx(idx);
         bond->setEndAtomIdx(idx_for_H);
         bool take_ownership = true;
         mol->addBond(bond, take_ownership);
      }
   }
   return r;
}

void
coot::mogulify_mol(RDKit::RWMol &mol) {

   // Mogul doesn't care about charges
   // charge_guanidinos(&mol);

   // unfused aromatic rings should have aromatic bonds, not single
   // and double bonds (fused rings should be kekulized).
   //
   // I've checked this - it seems that mogul can cope with single and
   // double bonds in the mol file for benzene (the target geometry is
   // symmetric).

   // aromatic rings that are pi-bonded to metals should have aromatic
   // bonds in the ring.

   // perchlorate: 3 double bonds and 1 single bond.

   // phosphate: phosphonate, phosphinate: use one double bond and 3 single bonds.

   // sulfone and sulfonamide use 2 double S=O bonds

   // sulfoxide: use a double S=O bond
   
   mogulify_nitro_groups(&mol);
   
   // std::cout << "---------------------- before ----------------" << std::endl;
   // debug_rdkit_molecule(&mol);
   // RDKit::MolOps::Kekulize(mol);
   // std::cout << "---------------------- after  ----------------" << std::endl;
   // debug_rdkit_molecule(&mol);
   
}

void
coot::charge_guanidinos(RDKit::RWMol *rdkm) {

   // std::cout << "-------------------- trying to delocalize_guanidinos() " << std::endl;

   // find an sp2 hybridized C, of degree 3 connected to 3 Ns.  Deloc
   // those three bonds.
   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {

      if ((*ai)->getAtomicNum() == 6) {
         RDKit::Atom *C_at = *ai;
         int idx_c = C_at->getIdx();
         unsigned int degree = rdkm->getAtomDegree(C_at);
         if (degree == 3) { 
            std::vector<RDKit::Bond *> CN_bonds;
            RDKit::Bond *C_N_double_bond = NULL;
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(C_at);
            while(nbrIdx != endNbrs) {
               if (rdkm->getAtomWithIdx(*nbrIdx)->getAtomicNum() == 7) { 
                  RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_c, *nbrIdx);
                  if (bond) {
                     CN_bonds.push_back(bond);
                     if (!C_N_double_bond) {
                        if (bond->getBondType() == RDKit::Bond::DOUBLE)
                           C_N_double_bond = bond;
                     } else {
                        C_N_double_bond = NULL; // !! something strange
                     } 
                  }
               }
               ++nbrIdx;
            }
            // std::cout << "found " << CN_bonds.size() << " N bonds to this C " << std::endl;
            if (CN_bonds.size() == 3) {
               if (C_N_double_bond) { 
                  int idx_n = C_N_double_bond->getOtherAtomIdx(idx_c);
                  // (*rdkm)[idx_n]->setFormalCharge(+1);
               }
            }
         }
      }
   }
}

void
coot::mogulify_nitro_groups(RDKit::RWMol *rdkm) {


   // std::cout << "--------------------- mogulify nitros --------------- " << std::endl;
   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {
      // do we have a N connected to an O via a double and an O via a single?
      // 
      // (if the molecule is sanitized, then the N should have a
      // charge too - that is ignored for now)

      if ((*ai)->getAtomicNum() == 7) {
         RDKit::Atom *N_at = *ai;
         int idx_c = N_at->getIdx();
         unsigned int degree = rdkm->getAtomDegree(N_at);
         if (degree == 3) {
            // fill these if you can
            RDKit::Bond *double_bond = NULL;
            RDKit::Bond *single_bond = NULL;
            
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(N_at);
            while(nbrIdx != endNbrs) {
               if (rdkm->getAtomWithIdx(*nbrIdx)->getAtomicNum() == 8) { 
                  RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_c, *nbrIdx);
                  if (bond) {
                     if (bond->getBondType() == RDKit::Bond::DOUBLE)
                        double_bond = bond;
                     if (bond->getBondType() == RDKit::Bond::SINGLE)
                        single_bond = bond;
                  }
               }
               ++nbrIdx;
            }

            if (double_bond && single_bond) {
               single_bond->setBondType(RDKit::Bond::DOUBLE);
            }
         }
      }
   }
}





lig_build::molfile_molecule_t
coot::make_molfile_molecule(const RDKit::ROMol &rdkm, int iconf) {

   lig_build::molfile_molecule_t mol;
   int n_conf  = rdkm.getNumConformers();

   if (n_conf) {
      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

      RDKit::Conformer conf = rdkm.getConformer(iconf);
      int n_mol_atoms = rdkm.getNumAtoms();

      for (int iat=0; iat<n_mol_atoms; iat++) {
         const RDKit::Atom* at_p = rdkm[iat];
         RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
         std::string name = "ZZZZ";
         try {
            at_p->getProp("name", name);
         }
         catch (const KeyErrorException &kee) {
            // std::cout << "caught no-name for atom exception in make_molfile_molecule(): "
            // <<  kee.what() << std::endl;

            try {
               at_p->getProp("_Name", name);  // RDKit's version of an atom name.
            }
            catch  (const KeyErrorException &kee) {
               // std::cout << "no _Name for " << at_p << " " << kee.what() << std::endl;
               try {
                  at_p->getProp("_TriposAtomName", name);  // RDKit's version of an atom name from a Mol2 File.
               }
               catch  (const KeyErrorException &kee) {
                  if (0)
                     std::cout << "no name, _Name or _TriposAtomName for "
                               << at_p << " " << kee.what() << std::endl;
               }
            }
         }
         clipper::Coord_orth pos(r_pos.x, r_pos.y, r_pos.z);
         int n = at_p->getAtomicNum();
         std::string element = tbl->getElementSymbol(n);
         int charge = at_p->getFormalCharge();
         lig_build::molfile_atom_t mol_atom(pos, element, name);

         mol_atom.formal_charge = charge;
         RDKit::Atom::ChiralType ct = at_p->getChiralTag();
         if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CW)
            mol_atom.chiral = RDKit::Atom::CHI_TETRAHEDRAL_CW;
         if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CCW)
            mol_atom.chiral = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
         mol_atom.aromatic = at_p->getIsAromatic();
         if (false)
            std::cout << "added atom " << mol_atom << std::endl;
         mol.add_atom(mol_atom);
      }

      unsigned int n_bonds = rdkm.getNumBonds();
      for (unsigned int ib=0; ib<n_bonds; ib++) {
         const RDKit::Bond *bond_p = rdkm.getBondWithIdx(ib);
         int idx_1 = bond_p->getBeginAtomIdx();
         int idx_2 = bond_p->getEndAtomIdx();
         lig_build::bond_t::bond_type_t bt = convert_bond_type(bond_p->getBondType());
         if (false)
            std::cout << "   make_molfile_molecule() " << idx_1 << " " << idx_2 << " from type  "
                      << bond_p->getBondType() << " to " << bt << std::endl;
         lig_build::molfile_bond_t mol_bond(idx_1, idx_2, bt);
         RDKit::Bond::BondDir bond_dir = bond_p->getBondDir();
         if (bond_dir != RDKit::Bond::NONE) {
            if (bond_dir == RDKit::Bond::BEGINWEDGE) { 
               mol_bond.bond_type = lig_build::bond_t::OUT_BOND;
               // std::cout << "found a BEGINWEDGE between atoms "
               // << mol.atoms[idx_1] << " and " << mol.atoms[idx_2]
               // << std::endl;
            } 
            if (bond_dir == RDKit::Bond::BEGINDASH) { 
               mol_bond.bond_type = lig_build::bond_t::IN_BOND;
            }
         }
         mol.add_bond(mol_bond);
      }
   }
   return mol;
}

// returns NULL on fail. Caller deletes.
//
mmdb::Residue *
coot::make_residue(const RDKit::ROMol &rdkm, int iconf, const std::string &res_name) {

   // replace this function by making a residue directly instead of via a molfile.
   // If there are no atom names, make them from the element and atom number

   mmdb::Residue *residue_p = NULL;
   lig_build::molfile_molecule_t mol = coot::make_molfile_molecule(rdkm, iconf);

   // now convert mol to a mmdb::Residue *
   // 
   if (mol.atoms.size()) {
      residue_p = new mmdb::Residue;
      residue_p->seqNum = 1;
      residue_p->SetResName(res_name.c_str());
      mmdb::Chain *chain_p = new mmdb::Chain;
      chain_p->SetChainID("");
      chain_p->AddResidue(residue_p);
      for (unsigned int iat=0; iat<mol.atoms.size(); iat++) { 
         mmdb::Atom *at = new mmdb::Atom;
         std::string atom_name = mol.atoms[iat].name; // overridden hopefully
         at->SetAtomName(atom_name.c_str());
         at->SetElementName(mol.atoms[iat].element.c_str());
         at->SetCoordinates(mol.atoms[iat].atom_position.x(),
                            mol.atoms[iat].atom_position.y(),
                            mol.atoms[iat].atom_position.z(),
                            1.0, 30.0);
         at->Het = 1;
         residue_p->AddAtom(at);
      }
   }

   return residue_p;
}

// coot::dictionary_residue_restraints_t
// coot::make_dictionary(const RDKit::ROMol &rdkm, int iconf, const std::string &res_name) {

//    coot::dictionary_residue_restraints_t d;

//    return d;

// } 


lig_build::bond_t::bond_type_t
coot::convert_bond_type(const RDKit::Bond::BondType &type) {

   lig_build::bond_t::bond_type_t t = lig_build::bond_t::SINGLE_BOND; // Hmmm..
   
   if (type == RDKit::Bond::SINGLE)
      t = lig_build::bond_t::SINGLE_BOND;
   if (type == RDKit::Bond::DOUBLE)
      t = lig_build::bond_t::DOUBLE_BOND;
   if (type == RDKit::Bond::TRIPLE)
      t = lig_build::bond_t::TRIPLE_BOND;
   if (type == RDKit::Bond::AROMATIC)
      t = lig_build::bond_t::AROMATIC_BOND;
    
    // PUTS_IT_BACK?
    // enabling this makes phosphate bonds disconnected in the
    // thumbnail and Residue->2D
    // 
    //     if (type == RDKit::Bond::ONEANDAHALF)
    //        t = lig_build::bond_t::DELOC_ONE_AND_A_HALF;

   return t;
}

// fiddle with rdkm.  This can throw a std::exception
// 
void
coot::remove_non_polar_Hs(RDKit::RWMol *rdkm) {

   unsigned int n_bonds = rdkm->getNumBonds();
   // std::vector<RDKit::ATOM_SPTR> atoms_to_be_deleted;
   std::vector<RDKit::Atom *> atoms_to_be_deleted;
   for (unsigned int ib=0; ib<n_bonds; ib++) {
      RDKit::Bond *bond_p = rdkm->getBondWithIdx(ib);
      int idx_1 = bond_p->getBeginAtomIdx();
      int idx_2 = bond_p->getEndAtomIdx();

      RDKit::Atom* at_p_1 = (*rdkm)[idx_1];
      RDKit::Atom* at_p_2 = (*rdkm)[idx_2];
      // If this was a bond for a hydrogen attached to a carbon, delete it.
      if ((at_p_1->getAtomicNum() == 1) && (at_p_2->getAtomicNum() == 6)) {
            // rdkm->removeBond(idx_1, idx_2);
            atoms_to_be_deleted.push_back(at_p_1);
         }
      if ((at_p_2->getAtomicNum() == 1) && (at_p_1->getAtomicNum() == 6)) {
            // rdkm->removeBond(idx_1, idx_2);
            atoms_to_be_deleted.push_back(at_p_2);
         }
      if ((at_p_1->getAtomicNum() == 1) && (at_p_2->getAtomicNum() == 5)) {
            // rdkm->removeBond(idx_1, idx_2);
            atoms_to_be_deleted.push_back(at_p_1);
         }
      if ((at_p_2->getAtomicNum() == 1) && (at_p_1->getAtomicNum() == 5)) {
            // rdkm->removeBond(idx_1, idx_2);
            atoms_to_be_deleted.push_back(at_p_2);
         }
   }
   for (unsigned int i=0; i<atoms_to_be_deleted.size(); i++) { 
      rdkm->removeAtom(atoms_to_be_deleted[i]);
   }

   /// std::cout << "DEBUG:: remove_non_polar_Hs() clearComputedProps() " << std::endl;
   rdkm->clearComputedProps();
   // clean up things like nitro groups
   // std::cout << "DEBUG:: remove_non_polar_Hs() cleanUp() " << std::endl;
   RDKit::MolOps::cleanUp(*rdkm);
   // update computed properties on atoms and bonds:

   coot::assign_formal_charges(rdkm);

   // std::cout << "DEBUG:: remove_non_polar_Hs() updatePropertyCache() " << std::endl;
   rdkm->updatePropertyCache();

   // std::cout << "DEBUG:: remove_non_polar_Hs() kekulize() " << std::endl;
   RDKit::MolOps::Kekulize(*rdkm);
   // std::cout << "DEBUG:: remove_non_polar_Hs() assignRadicals() " << std::endl;
   RDKit::MolOps::assignRadicals(*rdkm);
   // set conjugation
   // std::cout << "DEBUG:: remove_non_polar_Hs() setConjugation() " << std::endl;
   RDKit::MolOps::setConjugation(*rdkm);
   // set hybridization
   // std::cout << "DEBUG:: remove_non_polar_Hs() setHybridization() " << std::endl;
   RDKit::MolOps::setHybridization(*rdkm);
   // remove bogus chirality specs:
   // std::cout << "DEBUG:: remove_non_polar_Hs() cleanupChirality() " << std::endl;
   RDKit::MolOps::cleanupChirality(*rdkm);

}


// Delete a hydrogen (if possible) from a N with valence 4.
void
coot::delete_excessive_hydrogens(RDKit::RWMol *rdkm) {

   unsigned int n_mol_atoms = rdkm->getNumAtoms();
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {

#if (RDKIT_VERSION >= RDKIT_VERSION_CHECK(2018, 3, 1))
      RDKit::Atom* at_p = (*rdkm)[iat];
#else
         RDKit::Atom *at_p = (*rdkm)[iat];
#endif
      if (at_p->getAtomicNum() == 7) {

#if (RDKIT_VERSION >= RDKIT_VERSION_CHECK(2017, 3, 1))
         RDKit::Atom::ValenceType which = RDKit::Atom::ValenceType::EXPLICIT;
         int e_valence = at_p->getValence(which);
#else
         int e_valence = at_p->getExplicitValence();
#endif

         // std::cout << " atom N has explicit valence: " << e_valence << std::endl;

         if (e_valence == 4) {

            RDKit::ROMol::OEDGE_ITER current, end;
            boost::tie(current, end) = rdkm->getAtomBonds(at_p);
            RDKit::Atom *last_hydrogen_p = NULL;
            while (current != end) {
               RDKit::Bond* bond=(*rdkm)[*current];
               // is this a bond to a hydrogen?
               int idx = bond->getOtherAtomIdx(iat);
               RDKit::Atom* at_other_p = (*rdkm)[iat];
               if (at_other_p->getAtomicNum() == 1)
                  last_hydrogen_p = at_other_p;
               current++;
            }

            if (last_hydrogen_p) {
               // delete it then
               rdkm->removeAtom(last_hydrogen_p);
            }
         }
      }
   }
}

// Put a +1 charge on Ns with 4 bonds,
void
coot::assign_formal_charges(RDKit::RWMol *rdkm) {

   bool debug = false;
   int n_mol_atoms = rdkm->getNumAtoms();
   if (debug)
      std::cout << "---------------------- in assign_formal_charges() with " << n_mol_atoms
                << " atoms -----------" << std::endl;

   for (int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::Atom* at_p = (*rdkm)[iat];
      // debug
      if (0)
         std::cout << "in assign_formal_charges() calcExplicitValence on atom "
                   << iat << "/" << n_mol_atoms
                   << "  " << at_p->getAtomicNum() << std::endl;
      at_p->calcExplicitValence(false);
   }

   for (int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::Atom* at_p = (*rdkm)[iat];
      if (debug) {

#if (RDKIT_VERSION >= RDKIT_VERSION_CHECK(2017, 3, 1))
         RDKit::Atom::ValenceType which = RDKit::Atom::ValenceType::EXPLICIT;
         int v = at_p->getValence(which);
#else
         int v = at_p->getExplicitValence();
#endif
         std::cout << "atom " << iat << "/" << n_mol_atoms << "  " << at_p->getAtomicNum()
                   << " with valence " << v << std::endl;
      }
      if (at_p->getAtomicNum() == 7) { // N
         if (debug)
            std::cout << " incoming atom N has charge: " << at_p->getFormalCharge() << std::endl;

#if (RDKIT_VERSION >= RDKIT_VERSION_CHECK(2017, 3, 1))
         RDKit::Atom::ValenceType which = RDKit::Atom::ValenceType::EXPLICIT;
         int e_valence = at_p->getValence(which);
#else
         int e_valence = at_p->getExplicitValence();
#endif
         if (debug)
            std::cout << " atom N has explicit valence: " << e_valence << std::endl;
         if (e_valence == 4) {
            if (debug)
               std::cout << ".......... assign_formal_charges: found a N with valence 4..."
                         << at_p << std::endl;
            at_p->setFormalCharge(1);
         }
         if (debug)
            std::cout << " atom N has charge: " << at_p->getFormalCharge() << std::endl;
      }
      if (at_p->getAtomicNum() == 12) { // Mg
         at_p->setFormalCharge(2);
      }
   }

   charge_phosphates(rdkm);

   if (debug)
      std::cout << "----------- normal completion of assign_formal_charges()" << std::endl;
}

// a wrapper for the above, matching hydrogens names to the
// dictionary.  Add atoms to residue_p, return success status.
//
// This calls undelocalise.  Is that what we want to do?
//
std::pair<bool, std::string>
coot::add_hydrogens_with_rdkit(mmdb::Residue *residue_p,
                              const coot::dictionary_residue_restraints_t &restraints) {

   bool r = 0;
   std::string error_message;
   if (residue_p) {
      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
      try {

         // first save the existing hydrogen names. We don't want to
         // add a hydrogen with the same name as an atom we already
         // have.
         //
         std::vector<std::string> existing_H_names;
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            if (! residue_atoms[iat]->isTer()) {
               std::string ele = residue_atoms[iat]->element;
               if (ele == " H")
                  existing_H_names.push_back(residue_atoms[iat]->name);
            }
         }

         // get out now if thre are existing hydrogens,
         // MolOps::addHs() throws an exception (not clear to me why)
         // if there are hydrogens on the molecule already.
         if (existing_H_names.size()) {
            return std::pair<bool, std::string> (0, "Ligand contains (some) hydrogens already");
         }


         std::vector<std::string> residue_alt_confs = util::get_residue_alt_confs(residue_p);

         for (unsigned int i_alt_conf=0; i_alt_conf<residue_alt_confs.size(); i_alt_conf++) { 

            std::string alt_conf = residue_alt_confs[i_alt_conf];
            RDKit::RWMol m_no_Hs = rdkit_mol(residue_p, restraints, residue_alt_confs[i_alt_conf]);
            unsigned int n_mol_atoms = m_no_Hs.getNumAtoms();

            coot::undelocalise(&m_no_Hs);
            for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
               RDKit::Atom* at_p = m_no_Hs[iat];
               at_p->calcImplicitValence(true);
            }

            bool explicit_only = false;
            bool add_coords = true;
            RDKit::ROMol m_no_Hs_ro(m_no_Hs);
            RDKit::ROMol *m_pre = RDKit::MolOps::addHs(m_no_Hs_ro, explicit_only, add_coords);
            RDKit::RWMol m(*m_pre);
            // I think m_pre should be deleted.
            delete m_pre;
            unsigned int n_atoms_new = m.getNumAtoms();
            unsigned int n_conf = m.getNumConformers();

            double vdwThresh=10.0;
            int confId = 0;
            bool ignoreInterfragInteractions=true;
            int maxIters = 500;

            ForceFields::ForceField *ff =
               RDKit::UFF::constructForceField(m, vdwThresh, confId,
                                               ignoreInterfragInteractions);

            for (unsigned int iat=0; iat<n_mol_atoms; iat++)
               ff->fixedPoints().push_back(iat);

            ff->initialize();
            int res=ff->minimize(maxIters);
            delete ff;
         

            if (! n_conf) {
               std::cout << "ERROR:: mol with Hs: no conformers" << std::endl;
            } else { 
               RDKit::Conformer conf = m.getConformer(0);
               std::vector<std::string> H_names_already_added;
         
               for (unsigned int iat=0; iat<n_atoms_new; iat++) {
                  RDKit::Atom* at_p = m[iat];
                  RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
                  std::string name = "";
                  try {
                     at_p->getProp("name", name);
                     mmdb::Atom *res_atom = residue_p->GetAtom(name.c_str());
                     if (res_atom) {
                        std::cout << "setting heavy atom " << name << " to "
                                  << r_pos << std::endl;
                        res_atom->x = r_pos.x;
                        res_atom->y = r_pos.y;
                        res_atom->z = r_pos.z;
                     }
                     
                  }
                  catch (const KeyErrorException &kee) {

                     // OK...
                     //
                     // typically when we get here, that's because the
                     // atom is a new one, generated by RDKit.
                  
                     std::string name_i = coot::infer_H_name(iat, at_p, &m, restraints,
                                                           H_names_already_added);
                     if (! name_i.empty()) {

                        // add atom if the name is not already there:
                        // 
                        if (std::find(existing_H_names.begin(),
                                      existing_H_names.end(),
                                      name_i) == existing_H_names.end()) {
                        
                           H_names_already_added.push_back(name_i);
                        
                           int n = at_p->getAtomicNum();
                           std::string element = tbl->getElementSymbol(n);
                        
                           mmdb::Atom *at = new mmdb::Atom;
                           at->SetAtomName(name_i.c_str());
                           // at->SetElementName(element.c_str()); // FIXME?
                           at->SetElementName(" H");  // PDBv3 FIXME
                           at->SetCoordinates(r_pos.x, r_pos.y, r_pos.z, 1.0, 30.0);
                           at->Het = 1;
                           if (alt_conf != "") {
                              strncpy(at->altLoc, alt_conf.c_str(), alt_conf.length()+1);
                           }
                           residue_p->AddAtom(at);
                           r = 1;
                        }
                     }
                  }
               }
            }
         }
         
         // delete m;
      }
      catch (const std::runtime_error &e) {
         std::cout << e.what() << std::endl;
      }
      catch (const std::exception &rdkit_error) {
         std::cout << rdkit_error.what() << std::endl;
      }
   }
   return std::pair<bool, std::string> (r, error_message);
}

// atom_p is a hydrogen atom we presume, of degree 1.  This is tested
// before the restraints are checked.
// 
std::string
coot::infer_H_name(int iat,
         RDKit::Atom* atom_p,
                   const RDKit::ROMol *mol,
                   const dictionary_residue_restraints_t &restraints,
                   const std::vector<std::string> &H_names_already_added) {

   std::string r = "";

   unsigned int deg = mol->getAtomDegree(atom_p);
   if (deg == 1) {
      RDKit::ROMol::OEDGE_ITER current, end;
      boost::tie(current, end) = mol->getAtomBonds(atom_p);
      while (current != end) {
            const RDKit::Bond* bond=(*mol)[*current];
            int idx = bond->getOtherAtomIdx(iat);
            const RDKit::Atom* other_atom_p = (*mol)[idx];
         std::string bonding_atom_name;
         try {
            other_atom_p->getProp("name", bonding_atom_name);
            // in the restraints, what is the name of the hydrogen
            // bonded to atom with name bonding_atom_name? (if any).
            std::vector<std::string> nv = 
               restraints.get_attached_H_names(bonding_atom_name);
            for (unsigned int i=0; i<nv.size(); i++) { 
               if (std::find(H_names_already_added.begin(),
                             H_names_already_added.end(), nv[i]) ==
                   H_names_already_added.end()) {
                  r = nv[i];
                  break;
               }
            }
         }
         catch (const KeyErrorException &kee) {
            // this should not happen, there should be no way we get
            // here where we have a hydrogen (check in calling
            // function) with no name attached to another atom with no
            // name.
            std::cout << "ERROR:: in infer_H_name() bonding atom with no name "
                      << std::endl;
         }
         current++;
      }
   }
   // std::cout << "returning infered H name :" << r << ":" << std::endl;
   return r;
}


// return a kekulized molecule
RDKit::RWMol
coot::remove_Hs_and_clean(const RDKit::ROMol &rdkm, bool set_aromaticity) {

   RDKit::ROMol *rdk_mol_with_no_Hs_ro = RDKit::MolOps::removeHs(rdkm);
   RDKit::RWMol rdk_mol_with_no_Hs = *rdk_mol_with_no_Hs_ro;

   // clear out any cached properties
   rdk_mol_with_no_Hs.clearComputedProps();
   // clean up things like nitro groups
   RDKit::MolOps::cleanUp(rdk_mol_with_no_Hs);
   rdk_mol_with_no_Hs.updatePropertyCache();

   RDKit::MolOps::assignRadicals(rdk_mol_with_no_Hs);

   if (set_aromaticity)
      RDKit::MolOps::setAromaticity(rdk_mol_with_no_Hs);

   // set conjugation
   RDKit::MolOps::setConjugation(rdk_mol_with_no_Hs);

   // set hybridization
   RDKit::MolOps::setHybridization(rdk_mol_with_no_Hs); // non-linear ester bonds, yay!

   // remove bogus chirality specs:
   RDKit::MolOps::cleanupChirality(rdk_mol_with_no_Hs);

   // When I add this some of the above might be redundant.
   //
   unsigned int ops = RDKit::MolOps::SANITIZE_KEKULIZE;
   RDKit::MolOps::sanitizeMol(rdk_mol_with_no_Hs, ops); //sets ringinfo

   RDKit::MolOps::Kekulize(rdk_mol_with_no_Hs);

   return rdk_mol_with_no_Hs;
} 



// tweak rdkmol.
// 
// Don't do anything if rdk_mol is not 3d (return -1).
// 
int
coot::add_2d_conformer(RDKit::ROMol *rdk_mol, double weight_for_3d_distances) {

   bool debug = false;

   int icurrent_conf = 0; // the conformer number from which the
                          // distance matrix is generated.  Should this
                          // be passed?

   unsigned int n_conf  = rdk_mol->getNumConformers();
   if (n_conf == 0) {
      std::cout << "WARNING:: no conformers in add_2d_conformer() - aborting"
                << std::endl;
      return -1;
   }

   if (! rdk_mol->getConformer(icurrent_conf).is3D())
      return -1;


   // OK, go........
   
   unsigned int n_mol_atoms = rdk_mol->getNumAtoms();   

   if (n_mol_atoms == 0) return -1;

   if (debug) 
      std::cout << "::::: add_2d_conformer before compute2DCoords n_atoms: "
                << rdk_mol->getConformer(0).getNumAtoms()
                << " n_bonds " << rdk_mol->getNumBonds() << std::endl;

   // We must call calcImplicitValence() before getNumImplictHs()
   // [that is to say that compute2DCoords() calls getNumImplictHs()
   // without calling calcImplicitValence() and that results in problems].
   //
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::Atom* at_p = (*rdk_mol)[iat];
      at_p->calcImplicitValence(true);
   }

   int n_items = n_mol_atoms * (n_mol_atoms - 1)/2;

   double *cData = new double[n_items]; // handled by smart pointer, I think.
   for (int i=0; i<n_items; i++)
      cData[i] = -1;

   // fill cData with distances (don't include hydrogen distance
   // metrics (makes layout nicer?)).
   // 
   RDKit::Conformer conf = rdk_mol->getConformer(icurrent_conf);
   int ic_index = 0;
   for (unsigned int iat=1; iat<n_mol_atoms; iat++) {
      RDKit::Atom* iat_p = (*rdk_mol)[iat];
      if (iat_p->getAtomicNum() != 1) { 
         RDGeom::Point3D &pos_1 = conf.getAtomPos(iat);
         // std::cout << "   in 3d conformer: pos " << iat << " is " << pos_1 << std::endl;
         for (unsigned int jat=0; jat<iat; jat++) {
          RDKit::Atom* jat_p = (*rdk_mol)[jat];
          if (jat_p->getAtomicNum() != 1) { 
               RDGeom::Point3D &pos_2 = conf.getAtomPos(jat);
               RDGeom::Point3D diff = pos_1 - pos_2;

               // (thanks to JED for useful discussions)
               ic_index = iat*(iat - 1)/2 + jat;

                if (iat < jat)
                   ic_index = jat*(jat -1)/2 + iat;

               if (ic_index >= n_items)
                  std::cout << "indexing problem! " << ic_index << " but limit "
                            << n_items << std::endl;
               if (false)
                  std::cout << "mimic: atoms " << iat << " " << jat
                            << " ic_index " << ic_index << " for max " << n_items
                            << " dist " << diff.length() << std::endl;
               cData[ic_index] = diff.length();
            }
         }
      }
   }

   RDDepict::DOUBLE_SMART_PTR dmat(cData);
   
   // AFAICS, the number of conformers before and after calling
   // compute2DCoords() is the same (1).
   // 
   // Does it clear first? Yes, optional parameter clearConfs=true.
   // all atoms, long along x-axis, flip 2 rotatable bonds, try it for
   // 20 samples
   // 
   // int iconf = RDDepict::compute2DCoords(*rdk_mol, NULL, 1, 0, 2, 20);
   //

   int nRB = RDKit::Descriptors::calcNumRotatableBonds(*rdk_mol);

   // The can screw up the coordinates of rdk_mol.  I don't know
   // why. It is not due to distance matrix (dmat) I am pretty sure of
   // that.  other confs are cleared, so this should return 0.
    int iconf =
       RDDepict::compute2DCoordsMimicDistMat(*rdk_mol, &dmat, true, true,
                                             weight_for_3d_distances, nRB, 200);

   conf = rdk_mol->getConformer(iconf);
   RDKit::WedgeMolBonds(*rdk_mol, &conf);

   if (debug) { // .................... debug ...................
      std::cout << "::::: add_2d_conformer after  compute2DCoords n_atoms: "
                << rdk_mol->getConformer(0).getNumAtoms()
                << " n_bonds " << rdk_mol->getNumBonds() << std::endl;
   
      std::cout << ":::::: in add_2d_conformer here are the coords: "
                << std::endl;
      conf = rdk_mol->getConformer(iconf);
      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
         RDKit::Atom* at_p = (*rdk_mol)[iat];
         RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
         std::string name = "";
         try {
            at_p->getProp("name", name);
         }
         catch (const KeyErrorException &kee) {
         }
         std::cout << "   " << iat << " " << name << "  " << r_pos << std::endl;
      }
   }

   return iconf;
}

// undelocalise if you can, i.e. there is are 2 deloc bonds to a atom
// (a carbon), make one of these a double (the second, non-N bound)
// and the first (N-C) a single bond.
// 
void
coot::undelocalise(RDKit::RWMol *rdkm) {

   charge_undelocalized_guanidinos(rdkm);
   undelocalise_aminos(rdkm);
   undelocalise_nitros(rdkm);
   undelocalise_methyl_carboxylates(rdkm);
   undelocalise_carboxylates(rdkm); // after above
   undelocalise_phosphates(rdkm);
   undelocalise_sulphates(rdkm);
   charge_metals(rdkm);
   charge_sp3_borons(rdkm);
}


void
coot::undelocalise_aminos(RDKit::RWMol *rdkm) {


   // Plan:
   //
   // 1) Identify a bond with ONEANDAHALF bond order (bond_1)
   // 2) find a nitrogen in the bond
   // 3)   if yes, then is the other atom a carbon
   // 4)      if yes, then find another bond to this carbon that
   //            is deloc (bond_2)
   // 5)         if the other atom of this other bond an oxygen?
   // 6)            if yes, then make bond_1 single, bond_2 double

   bool debug = false;

   int n_bonds = rdkm->getNumBonds();
   RDKit::ROMol::BondIterator bondIt;
   RDKit::ROMol::BondIterator bondIt_inner;
   for(bondIt=rdkm->beginBonds(); bondIt!=rdkm->endBonds(); ++bondIt) {
      if ((*bondIt)->getBondType() == RDKit::Bond::ONEANDAHALF) {
         // was one of these atoms a Nitrogen?
         RDKit::Atom *atom_1 = (*bondIt)->getBeginAtom();
         RDKit::Atom *atom_2 = (*bondIt)->getEndAtom();
         bool do_it = 0;
         if (atom_1->getAtomicNum() == 7) {
            if (atom_2->getAtomicNum() == 6) {
               do_it = 1;
            }
         }
         if (atom_2->getAtomicNum() == 7) {
            if (atom_1->getAtomicNum() == 6) {
               do_it = 1;
               std::swap(atom_1, atom_2);
            }
         }

         if (do_it) {

               // atom_1 is a Nitrogen, atom_2 is a Carbon.  Does the
            // carbon have a bond (not this one) to an oxygen that is
            // delocalised?
            //
            RDKit::Atom *N_at = atom_1;
            RDKit::Atom *C_at = atom_2;
            
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(C_at);
            while(nbrIdx != endNbrs) {
               const RDKit::Atom* at = (*rdkm)[*nbrIdx];
               if (at->getAtomicNum() == 8) { 
                  RDKit::Bond *bond_inner = rdkm->getBondBetweenAtoms(C_at->getIdx(), *nbrIdx);
                  if (bond_inner) {
                     if (bond_inner->getBondType() == RDKit::Bond::ONEANDAHALF) {
                        (*bondIt)->setBondType(RDKit::Bond::SINGLE);
                        bond_inner->setBondType(RDKit::Bond::DOUBLE);
                     }
                  }
               }
               ++nbrIdx;
            }
         }
      }
   }
}

void
coot::undelocalise_nitros(RDKit::RWMol *rdkm) {

   // std::cout << "--------- undelocalise_nitros() " << std::endl;

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {
      if ((*ai)->getAtomicNum() == 7) {
         RDKit::Atom *N_at = *ai;
         int idx_n = N_at->getIdx();
         unsigned int degree = rdkm->getAtomDegree(N_at);
         if (degree == 3) {
            // fill these if you can
            std::vector<RDKit::Bond *> deloc_bonds;
            
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(N_at);
            while(nbrIdx != endNbrs) {
             const RDKit::Atom* at = (*rdkm)[*nbrIdx];
               if (rdkm->getAtomWithIdx(*nbrIdx)->getAtomicNum() == 8) { 
                  RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_n, *nbrIdx);
                  if (bond) {
                     if (bond->getBondType() == RDKit::Bond::ONEANDAHALF) {
                        deloc_bonds.push_back(bond);
                     }
                  }
               }
               ++nbrIdx;
            }

            if (deloc_bonds.size() == 2) {
               deloc_bonds[0]->setBondType(RDKit::Bond::DOUBLE);
               deloc_bonds[1]->setBondType(RDKit::Bond::SINGLE);
               int idx_O = deloc_bonds[1]->getOtherAtomIdx(idx_n);
               // mogul ignores these, I think
               (*rdkm)[idx_O]->setFormalCharge(-1);
               N_at->setFormalCharge(+1);
            } 
         }
      }
   }
}



// run this after undelocalise_methyl_carboxylates.
// 
void
coot::undelocalise_carboxylates(RDKit::RWMol *rdkm) {

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {

      // Is there a carbon that is deloc attached to 2 oxygens.  (
      if ((*ai)->getAtomicNum() == 6) {
         RDKit::Atom *C_at = *ai;
         int idx_c = C_at->getIdx();
         std::vector<RDKit::Bond *> deloc_O_bonds;
         RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
         boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(C_at);
         while(nbrIdx != endNbrs) {
            const RDKit::Atom* at = (*rdkm)[*nbrIdx];
            RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_c, *nbrIdx);
            if (bond) {
               if (bond->getBondType() == RDKit::Bond::ONEANDAHALF)
                  deloc_O_bonds.push_back(bond);
            }
            ++nbrIdx;
         }

         if (deloc_O_bonds.size() == 2) {
            deloc_O_bonds[0]->setBondType(RDKit::Bond::SINGLE);
            deloc_O_bonds[1]->setBondType(RDKit::Bond::DOUBLE);
            int idx_o = deloc_O_bonds[0]->getOtherAtomIdx(idx_c);
            (*rdkm)[idx_o]->setFormalCharge(-1);
         }
      }
   }
}


void
coot::undelocalise_methyl_carboxylates(RDKit::RWMol *rdkm) {

   // The valence of 2 1/2 on a methyl on one of the oxygens of a carboxylate
   // 
   RDKit::ROMol::BondIterator bondIt;
   RDKit::ROMol::BondIterator bondIt_inner;
   for(bondIt=rdkm->beginBonds(); bondIt!=rdkm->endBonds(); ++bondIt) {
      if ((*bondIt)->getBondType() == RDKit::Bond::ONEANDAHALF) {
         RDKit::Atom *atom_1 = (*bondIt)->getBeginAtom();
         RDKit::Atom *atom_2 = (*bondIt)->getEndAtom();

         if (atom_1->getAtomicNum() == 6) {
            if (atom_2->getAtomicNum() == 8) {

               // rename for clarity
               RDKit::Atom *central_C = atom_1;
               RDKit::Atom *O1 = atom_2;

               for(bondIt_inner=rdkm->beginBonds(); bondIt_inner!=rdkm->endBonds(); ++bondIt_inner) {
                  if ((*bondIt_inner)->getBondType() == RDKit::Bond::ONEANDAHALF) {
                     RDKit::Atom *atom_1_in = (*bondIt_inner)->getBeginAtom();
                     RDKit::Atom *atom_2_in = (*bondIt_inner)->getEndAtom();
                     if (atom_1_in == central_C) {
                        if (atom_2_in != O1) {
                           if (atom_2_in->getAtomicNum() == 8) {

                              // OK, we have a carbon (atom_1) bonded to two Os via delocs -
                              // the oxygens are atom_2 and atom_2_in
                              // 
                              // rename for clarity
                              //
                              RDKit::Atom *O2 = atom_2_in;

                              // bondIt and bondIt_inner are the bonds that we will ultimately modify
                              // 
                              deloc_O_check_inner(rdkm, central_C, O1, O2, *bondIt, *bondIt_inner);

                           }
                        }
                     }

                     // The central carbon was the other atom?
                     if (atom_2_in == central_C) {
                        if (atom_1_in != O1) {
                           if (atom_1_in->getAtomicNum() == 8) {

                              // OK, we have a carbon (atom_1) bonded to two Os via delocs -
                              // the oxygens are atom_2 and atom_2_in
                              // 
                              // rename for clarity
                              //
                              RDKit::Atom *O2 = atom_1_in;
                              deloc_O_check_inner(rdkm, central_C, O1, O2, *bondIt, *bondIt_inner);
                           } 
                        } 
                     }
                  }
               }
            }
         }
            
         if (atom_1->getAtomicNum() == 8) {
            if (atom_2->getAtomicNum() == 6) {
               // rename for clarity
               RDKit::Atom *central_C = atom_2;
               RDKit::Atom *O1 = atom_1;
               
               for(bondIt_inner=rdkm->beginBonds(); bondIt_inner!=rdkm->endBonds(); ++bondIt_inner) {
                  if ((*bondIt_inner)->getBondType() == RDKit::Bond::ONEANDAHALF) {
                     RDKit::Atom *atom_1_in = (*bondIt_inner)->getBeginAtom();
                     RDKit::Atom *atom_2_in = (*bondIt_inner)->getEndAtom();
                     if (atom_1_in == central_C) {
                        if (atom_2_in != O1) {
                           if (atom_2_in->getAtomicNum() == 8) {

                              // again, we have detected a central carbon bonded to two
                              // oxygens with deloc bonds.
                              RDKit::Atom *O2 = atom_2_in;
                              deloc_O_check_inner(rdkm, central_C, O1, O2, *bondIt, *bondIt_inner);
                           }
                        }
                     }
                     // The central carbon was the other atom?
                     if (atom_2_in == central_C) {
                        if (atom_1_in != O1) {
                           if (atom_1_in->getAtomicNum() == 8) {
                              RDKit::Atom *O2 = atom_1_in;
                              deloc_O_check_inner(rdkm, central_C, O1, O2, *bondIt, *bondIt_inner);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}



// fiddle with the bonds in rdkm as needed.
void
coot::deloc_O_check_inner(RDKit::RWMol *rdkm, RDKit::Atom *central_C,
                          RDKit::Atom *O1, RDKit::Atom *O2,
                          RDKit::Bond *b1, RDKit::Bond *b2) {

   // std::cout << "debug:: deloc_O_check_inner: " << rdkm << " " << central_C << std::endl;

   RDKit::ROMol::BondIterator bondIt_in_in;
   // OK, so was there something attached to either of the Oxygens?
   //
   for(bondIt_in_in=rdkm->beginBonds(); bondIt_in_in!=rdkm->endBonds(); ++bondIt_in_in) {
      if ((*bondIt_in_in)->getBondType() == RDKit::Bond::SINGLE) {
         RDKit::Atom *atom_1_in_in = (*bondIt_in_in)->getBeginAtom();
         RDKit::Atom *atom_2_in_in = (*bondIt_in_in)->getEndAtom();

         // check atom_1_in_in vs the first oxygen
         if (atom_1_in_in == O1) {
            if (atom_2_in_in != central_C) {

               // OK, so O1 was bonded to something else
               //
               b1->setBondType(RDKit::Bond::SINGLE);
               b2->setBondType(RDKit::Bond::DOUBLE);
            }
         }

         // check vs the second oxygen
         if (atom_1_in_in == O2) {
            if (atom_2_in_in != central_C) {
               // OK, so O2 was bonded to something else
               b1->setBondType(RDKit::Bond::DOUBLE);
               b2->setBondType(RDKit::Bond::SINGLE);
            }
         }


         // check atom_2_in_in vs the first oxygen
         if (atom_2_in_in == O1) {
            if (atom_1_in_in != central_C) {
               // OK, so O1 was bonded to something else
               //
               b1->setBondType(RDKit::Bond::SINGLE);
               b2->setBondType(RDKit::Bond::DOUBLE);
            }
         }

         // check atom_2_in_in vs the second oxygen
         if (atom_2_in_in == O2) {
            if (atom_1_in_in != central_C) {
               // OK, so O2 was bonded to something else
               b1->setBondType(RDKit::Bond::DOUBLE);
               b2->setBondType(RDKit::Bond::SINGLE);
            }
         }
      }
   }
}

void
coot::undelocalise_phosphates(RDKit::ROMol *rdkm) {

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {

      if ((*ai)->getAtomicNum() == 15) {
         RDKit::Atom *P_at = *ai;
         int idx_1 = P_at->getIdx();
         std::vector<RDKit::Bond *> deloc_O_bonds;

         RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
         boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(P_at);
         while(nbrIdx != endNbrs) {
          const RDKit::Atom* at = (*rdkm)[*nbrIdx];
          RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_1, *nbrIdx);
            if (bond) {
               if (bond->getBondType() == RDKit::Bond::ONEANDAHALF)
                  deloc_O_bonds.push_back(bond);
            }
            ++nbrIdx;
         }

         if (deloc_O_bonds.size() == 4) {
            // PO4 monomer, PO4(-3)
            // make 1 single, one double.
            deloc_O_bonds[0]->setBondType(RDKit::Bond::DOUBLE);
            deloc_O_bonds[1]->setBondType(RDKit::Bond::SINGLE);
            deloc_O_bonds[2]->setBondType(RDKit::Bond::SINGLE);
            deloc_O_bonds[3]->setBondType(RDKit::Bond::SINGLE);
            // Handle formal charge too.
            int idx_o_1 = deloc_O_bonds[1]->getOtherAtomIdx(idx_1);
            int idx_o_2 = deloc_O_bonds[2]->getOtherAtomIdx(idx_1);
            int idx_o_3 = deloc_O_bonds[3]->getOtherAtomIdx(idx_1);
            RDKit::Atom* at_p_1 = (*rdkm)[idx_o_1];
            RDKit::Atom* at_p_2 = (*rdkm)[idx_o_2];
            RDKit::Atom* at_p_3 = (*rdkm)[idx_o_3];
            at_p_1->setFormalCharge(-1);
            at_p_2->setFormalCharge(-1);
            at_p_3->setFormalCharge(-1);
         }

         if (deloc_O_bonds.size() == 3) {
            // make 2 single and one double.  Handle formal charge too.
            deloc_O_bonds[0]->setBondType(RDKit::Bond::SINGLE);
            deloc_O_bonds[1]->setBondType(RDKit::Bond::SINGLE);
            deloc_O_bonds[2]->setBondType(RDKit::Bond::DOUBLE);
            int idx_o_0 = deloc_O_bonds[0]->getOtherAtomIdx(idx_1);
            int idx_o_1 = deloc_O_bonds[1]->getOtherAtomIdx(idx_1);

            RDKit::Atom* at_p_0 = (*rdkm)[idx_o_0];
            RDKit::Atom* at_p_1 = (*rdkm)[idx_o_1];
            at_p_0->setFormalCharge(-1);
            at_p_1->setFormalCharge(-1);
         }

         if (deloc_O_bonds.size() == 2) {
            // make 1 single, one double. Handle formal charge too.
            deloc_O_bonds[0]->setBondType(RDKit::Bond::SINGLE);
            deloc_O_bonds[1]->setBondType(RDKit::Bond::DOUBLE);
            int idx_o_0 = deloc_O_bonds[0]->getOtherAtomIdx(idx_1);
            RDKit::Atom* at_p_0 = (*rdkm)[idx_o_0];
            at_p_0->setFormalCharge(-1);
         }
      }
   }
}


void
coot::undelocalise_sulphates(RDKit::ROMol *rdkm) {

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {

      if ((*ai)->getAtomicNum() == 16) {
         RDKit::Atom *S_at = *ai;
         int idx_1 = S_at->getIdx();
         std::vector<RDKit::Bond *> deloc_O_bonds;

         RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
         boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(S_at);
         while(nbrIdx != endNbrs) {
            const RDKit::Atom* at = (*rdkm)[*nbrIdx];
            RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_1, *nbrIdx);
            if (bond) {
               if (bond->getBondType() == RDKit::Bond::ONEANDAHALF)
                  deloc_O_bonds.push_back(bond);
            }
            ++nbrIdx;
         }

         if (deloc_O_bonds.size() >= 3) {
            // SO4 monomer, SO4(-2)
            // make 1 single, two double.
            deloc_O_bonds[0]->setBondType(RDKit::Bond::DOUBLE);
            deloc_O_bonds[1]->setBondType(RDKit::Bond::DOUBLE);
            deloc_O_bonds[2]->setBondType(RDKit::Bond::SINGLE);
            if (deloc_O_bonds.size() == 4) {
               deloc_O_bonds[3]->setBondType(RDKit::Bond::SINGLE);
               // Handle formal charge too.
               int idx_o_2 = deloc_O_bonds[2]->getOtherAtomIdx(idx_1);
               int idx_o_3 = deloc_O_bonds[3]->getOtherAtomIdx(idx_1);
               RDKit::Atom* at_p_2 = (*rdkm)[idx_o_2];
               RDKit::Atom* at_p_3 = (*rdkm)[idx_o_3];
               at_p_2->setFormalCharge(-1);
               at_p_3->setFormalCharge(-1);
            } else {
               int idx_o_2 = deloc_O_bonds[2]->getOtherAtomIdx(idx_1); // this single-bonded O
               RDKit::Atom* at_p_2 = (*rdkm)[idx_o_2];
               at_p_2->setFormalCharge(-1);
            }
         }
      }
   }
}


void
coot::charge_sp3_borons(RDKit::RWMol *rdkm) {

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {
      if ((*ai)->getAtomicNum() == 5) {
         unsigned int degree = rdkm->getAtomDegree(*ai);
         if (degree == 4)
            (*ai)->setFormalCharge(-1);
      }
   }
}


void
coot::charge_metals(RDKit::RWMol *rdkm) {

   // hackety hack code.  Needs improvement/thinking about

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {
      if ((*ai)->getAtomicNum() == 11) { // Na
         (*ai)->setFormalCharge(+1);
      }
      if ((*ai)->getAtomicNum() == 12) { // Mg
         // std::cout << "............................... charging Mg" << std::endl;
         (*ai)->setFormalCharge(+2);
      }
      if ((*ai)->getAtomicNum() == 20) { // Ca
         (*ai)->setFormalCharge(+2);
      }
   }
}

void
coot::charge_undelocalized_guanidinos(RDKit::RWMol *rdkm) {

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {

      // Find a C with 3 deloc bonds to N.  If found, charge the C.
      if ((*ai)->getAtomicNum() == 6) {
         RDKit::Atom *C_at = *ai;
         int idx_c = C_at->getIdx();
         unsigned int degree = rdkm->getAtomDegree(C_at);
         int deloc_bond_count = 0;
         if (degree == 3) {
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(C_at);
            while(nbrIdx != endNbrs) {
               const RDKit::Atom* at = (*rdkm)[*nbrIdx];
               if (rdkm->getAtomWithIdx(*nbrIdx)->getAtomicNum() == 7) {
                  RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_c, *nbrIdx);
                  // std::cout << ".... found a C-N bond " << bond->getBondType() << std::endl;
                  if (bond->getBondType() == RDKit::Bond::ONEANDAHALF) {
                     // std::cout << "... it was ONEANDAHALF" << std::endl;
                     deloc_bond_count++;
                  }
               }
               ++nbrIdx;
            }
         }

         // std::cout << "...... deloc C-N bond count:" << deloc_bond_count << std::endl;
         if (deloc_bond_count == 3) {
            // std::cout << ".... charging the C" << std::endl;
            C_at->setFormalCharge(+1);
         }
      }
   }
}

// when using a refmac cif dictionary, we construct a molecule with
// deloc bonds (e.g. on a phosphate)
// valence on P: (1 1/2) * 3 + 1 -> 6 => problem.
//
// So, in that case, +1 charge the P.  This might be a hack.
//
// return the number of deleted atoms
void
coot::charge_phosphates(RDKit::RWMol *rdkm) {

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {

      if ((*ai)->getAtomicNum() == 15) {
         RDKit::Atom *P_at = *ai;
         int idx_1 = P_at->getIdx();
         std::vector<RDKit::Bond *> deloc_O_bonds;

         RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
         boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(P_at);
         while(nbrIdx != endNbrs) {
          const RDKit::Atom* at = (*rdkm)[*nbrIdx];
            RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_1, *nbrIdx);
            if (bond) {
               if (bond->getBondType() == RDKit::Bond::ONEANDAHALF)
                  deloc_O_bonds.push_back(bond);
            }
            ++nbrIdx;
         }

         if (deloc_O_bonds.size() == 3) {

            // a typical terminal phosphate (e.g. AMP)
            //
            // (intermediate phosphates in ATP have 1 + 1 + 11/2 + 11/2, which is OK)
            //

            P_at->setFormalCharge(1);
         }
      }
   }
}

// return the number of atoms added (e.g. -2)
int
coot::remove_phosphate_hydrogens(RDKit::RWMol *m, bool deloc_bonds) { 

   return remove_PO4_SO4_hydrogens(m, 15, deloc_bonds);

}

// return the number of atoms added (e.g. -1)
int
coot::remove_sulphate_hydrogens(RDKit::RWMol *m, bool deloc_bonds) {

   return remove_PO4_SO4_hydrogens(m, 16, deloc_bonds);

}


// return the number of atoms added (e.g. -1)
int
coot::remove_PO4_SO4_hydrogens(RDKit::RWMol *m, 
                               unsigned int atomic_num, 
                               bool deloc_bonds) {

   int n_added = 0;
   bool debug = false;
   RDKit::ROMol::AtomIterator ai;
   std::vector<RDKit::Atom *> H_atoms_to_be_deleted;
   for(ai=m->beginAtoms(); ai!=m->endAtoms(); ai++) {

      unsigned int this_atomic_num = (*ai)->getAtomicNum(); // convert int to unsigned int
      if (this_atomic_num == atomic_num) {
         RDKit::Atom *P_at = *ai;
         int idx_1 = P_at->getIdx();
         // std::cout << "new thingate centre " << P_at << " " << idx_1 << std::endl;
         std::vector<RDKit::Bond *> single_PO_bonds; // with a hydrogen attached
         std::vector<RDKit::Bond *> double_PO_bonds;
         std::vector<RDKit::Atom *> O_atoms_for_charging;
         std::vector<RDKit::Atom *> probable_phosphate_hydrogens;

         RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
         boost::tie(nbrIdx, endNbrs) = m->getAtomNeighbors(P_at);
         while (nbrIdx != endNbrs) {
            RDKit::Atom* at = (*m)[*nbrIdx];
            RDKit::Bond *bond = m->getBondBetweenAtoms(idx_1, *nbrIdx);
            if (bond) {

               if (at->getAtomicNum() == 8) {
                  if (bond->getBondType() == RDKit::Bond::SINGLE) { 
                     const int &idx_O = *nbrIdx;
                     const RDKit::Atom* O_at = at;

                     // Is the other atom of the O a hydrogen?
                     RDKit::ROMol::OEDGE_ITER current, end;
                     boost::tie(current, end) = m->getAtomBonds(O_at);
                     while (current != end) {

                        RDKit::Bond* o_bond=(*m)[*current];
                        // is this a bond to a hydrogen?
                        int idx_H = o_bond->getOtherAtomIdx(idx_O);
                        RDKit::Atom* at_other_p = (*m)[idx_H];
                        if (at_other_p->getAtomicNum() == 1) {
                           single_PO_bonds.push_back(bond);
                           O_atoms_for_charging.push_back(at);
                           probable_phosphate_hydrogens.push_back(at_other_p);
                        } else {
                           // std::cout << at_other_p << " was not a hydrogen" << std::endl;
                        }
                        current++;
                     }

                  }
                  if (bond->getBondType() == RDKit::Bond::DOUBLE) {
                     double_PO_bonds.push_back(bond);
                     // 20171217 surely we can't mean to charge an O with a double bond?
                     // O_atoms_for_charging.push_back(at.get());
                  }
               }
            }
            ++nbrIdx;
         }

         if (debug) {
            std::string ele = (atomic_num == 15) ? "P" : "S";
            std::cout << "DEBUG:: atom_idx: " << idx_1 << " Found "
                      << single_PO_bonds.size() << " single " << ele << "O bonds and "
                      << double_PO_bonds.size() << " double " << ele << "O bonds "
                      << std::endl;
         }

         bool do_strip_Hs = false;
         if (atomic_num == 15)
            if (single_PO_bonds.size() == 2 || single_PO_bonds.size() == 1) // terminal and mid PO4s
               if (double_PO_bonds.size() == 1)
                  do_strip_Hs = true;

         if (atomic_num == 16)
            if (single_PO_bonds.size() == 1)  // SO bonds of course in this case
               if (double_PO_bonds.size() == 2)
                  do_strip_Hs = true;

         if (do_strip_Hs) {

            if (debug) {
               std::string thingate = (atomic_num== 16) ? "sulphate" : "phosphate";
               std::cout << " :::::: found a " << thingate << " :::::::::::::" << std::endl;
            }

            for (unsigned int ip=0; ip<probable_phosphate_hydrogens.size(); ip++)
               H_atoms_to_be_deleted.push_back(probable_phosphate_hydrogens[ip]);

            // 20150622-PE
            // If we charge the Os, then (for AMP from PDBe-AMP.cif) we end up
            // with
            // "Explicit valence for atom # 1 O, 3, is greater than permitted"
            // when we call sanitizeMol() from hydrogen_transformations()
            // (directly after this function is called).
            //

            // 20170608 let's try to charge the O atoms:
            for (unsigned int ii=0; ii<O_atoms_for_charging.size(); ii++)
                O_atoms_for_charging[ii]->setFormalCharge(-1);

            if (deloc_bonds) { 

               if (O_atoms_for_charging.size() == 3) {

                  for (unsigned int ii=0; ii<single_PO_bonds.size(); ii++)
                     single_PO_bonds[ii]->setBondType(RDKit::Bond::ONEANDAHALF);
                  for (unsigned int ii=0; ii<double_PO_bonds.size(); ii++)
                     double_PO_bonds[ii]->setBondType(RDKit::Bond::ONEANDAHALF);
                  P_at->setFormalCharge(+1);

               }
            }
         }
      }
   }

   if (debug)
      std::cout << " in rdkit remove_PO4_SO4_hydrogens() remove these "
                << H_atoms_to_be_deleted.size() << " Hydrogen atoms for atomic number "
                << atomic_num << std::endl;

   for (unsigned int idel=0; idel<H_atoms_to_be_deleted.size(); idel++) {

      // remove bonds for these atoms then delete the atom
      //
      RDKit::ROMol::OEDGE_ITER current, end;
      boost::tie(current, end) = m->getAtomBonds(H_atoms_to_be_deleted[idel]);
      while (current != end) {
         RDKit::Bond* bond= (*m)[*current];
         int idx = H_atoms_to_be_deleted[idel]->getIdx();
         int idx_other = bond->getOtherAtomIdx(idx);
         if (debug) { // debug
            std::string name_1;
            std::string name_2;

            RDKit::Atom *other_at = (*m)[idx_other];
            H_atoms_to_be_deleted[idel]->getProp("name", name_1);
            other_at->getProp("name", name_2);
            std::cout << "----- removeBond between " << idx << " " << idx_other << " " << name_1 << " " << name_2
                      << std::endl;
         }
         m->removeBond(idx, idx_other);
         current++;
      }

      if (debug) {
         std::string atom_name;
         H_atoms_to_be_deleted[idel]->getProp("name", atom_name);
         std::cout << "-------- delete atom " << H_atoms_to_be_deleted[idel]
                                    << " " << atom_name << std::endl;
      }
      m->removeAtom(H_atoms_to_be_deleted[idel]);
      n_added--;
   }
   if (H_atoms_to_be_deleted.size() > 0) {

      std::string s = (H_atoms_to_be_deleted.size() > 1) ? "s" : "";
      std::string thingate = (atomic_num== 16) ? "sulphate" : "phosphate";
      std::cout << "INFO:: Deleted " << H_atoms_to_be_deleted.size()
                << " " << thingate << " hydrogen atom" << s  << std::endl;
   }
   // return the number of atoms added (e.g. -1)
   return n_added;
}

int
coot::remove_carboxylate_hydrogens(RDKit::RWMol *m, bool deloc_bonds) {

   // No HO2 on O2 in BEZ (benzoic acid)

   int n_added = 0;
   bool debug = false;
   RDKit::ROMol::AtomIterator ai;
   std::vector<RDKit::Atom *> H_atoms_to_be_deleted;
   for(ai=m->beginAtoms(); ai!=m->endAtoms(); ai++) {

      unsigned int this_atomic_num = (*ai)->getAtomicNum(); // convert int to unsigned int
      if (this_atomic_num == 6) {
         RDKit::Atom *C_at = *ai;
         int idx_C = C_at->getIdx();
         if (C_at->getDegree() == 3) {

            std::vector<RDKit::Bond *> single_CO_bonds; // with a hydrogen attached (presumably)
            std::vector<RDKit::Bond *> double_CO_bonds;
            std::vector<RDKit::Atom *> O_atoms_for_charging;
            std::vector<RDKit::Atom *> carboxylate_hydrogens;

            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = m->getAtomNeighbors(C_at);
            while (nbrIdx != endNbrs) {
             const RDKit::Atom* at = (*m)[*nbrIdx];
               RDKit::Bond *bond = m->getBondBetweenAtoms(idx_C, *nbrIdx);
               if (bond) {

                  if (at->getAtomicNum() == 8) {
                     if (bond->getBondType() == RDKit::Bond::SINGLE) {
                        single_CO_bonds.push_back(bond);
                     }
                     if (bond->getBondType() == RDKit::Bond::DOUBLE) {
                        double_CO_bonds.push_back(bond);
                     }
                  }
               }
               nbrIdx++;
            }
            if (single_CO_bonds.size() == 1) {
               if (double_CO_bonds.size() == 1) {
                  // was there an H atom on the other side of the single C-O bond?
                  RDKit::Bond *bond = single_CO_bonds[0];
                  RDKit::Atom *O_at = bond->getOtherAtom(C_at);
                  if (O_at->getDegree() == 2) {
                     int idx_O = O_at->getIdx();
                     RDKit::ROMol::ADJ_ITER nbrIdx_inner, endNbrs_inner;
                     boost::tie(nbrIdx_inner, endNbrs_inner) = m->getAtomNeighbors(O_at);
                     while (nbrIdx_inner != endNbrs_inner) {
                        const RDKit::Atom* at = (*m)[*nbrIdx_inner];
                        RDKit::Bond *bond_inner = m->getBondBetweenAtoms(idx_O, *nbrIdx_inner);
                        if (bond_inner) {
                           RDKit::Atom *at_H = bond_inner->getOtherAtom(O_at);
                           if (at_H->getAtomicNum() == 1) {
                              // delete this H, charge the O
                              m->removeAtom(at_H);
                              O_at->setFormalCharge(-1);
                           }
                        }
                        nbrIdx_inner++;
                     }
                  }
               }
            }
         }
      }
   }
   return n_added;

}


void
coot::debug_rdkit_molecule(const RDKit::ROMol *rdkm) {

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   std::cout << "---- Atoms: " << rdkm->getNumAtoms() << std::endl;
   for (unsigned int iat=0; iat<rdkm->getNumAtoms(); iat++) { 

      const RDKit::Atom* at_p = (*rdkm)[iat];
      std::string name;
      try {
         at_p->getProp("name", name);
      }
      catch (const KeyErrorException &err) {
      }
      catch (...) {
         // std::cout << "cannot happen" << std::endl; // yeah, it can (sigh)
      }
      int n = at_p->getAtomicNum();
      std::string element = tbl->getElementSymbol(n);
      unsigned int degree = rdkm->getAtomDegree(at_p);
      RDKit::Atom::HybridizationType ht = at_p->getHybridization();

      int f_c = at_p->getFormalCharge();
      std::cout << std::setw(3) << iat << " ele: " << std::setw(2) << std::right << element;
      if (! name.empty())
         std::cout << " name :" << name << ":";
      std::cout << " degree: " << degree;
      std::cout << " formal-charge: " << f_c << " ";
      std::cout << " hybridization: " << ht;
                                                \
      // chiral tag
      RDKit::Atom::ChiralType ct = at_p->getChiralTag();
      std::string cts = "!";
      if (ct == RDKit::Atom::CHI_UNSPECIFIED)     cts = "-";
      if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CW)  cts = " CW";
      if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CCW) cts = "CCW";
      if (ct == RDKit::Atom::CHI_OTHER)           cts = "Oth";
      std::cout << " Chir: " << cts;


      // R/S chirality
      std::string cip;
      try {
         at_p->getProp("_CIPCode", cip);
         std::cout << " CIP-Code " << cip;
      }
      catch (const KeyErrorException &err) {
         // Not an error
         // std::cout << "KeyErrorException " << err.what() << " for _CIPCode" << std::endl;
         std::cout << " CIP-Code - ";
      }
      unsigned int cip_rank;
      try {
         at_p->getProp(RDKit::common_properties::_CIPRank, cip_rank);
         std::cout << " CIP-Rank " << cip_rank;
      }
      catch (const KeyErrorException &err) {
         std::cout << " CIP-Rank - ";
      }
      catch (...) {
         std::cout << " CIP-Rank... - ";
      }
      std::cout << std::endl;
   }

   unsigned int n_bonds = rdkm->getNumBonds();
   std::cout << "---- Bonds: " << n_bonds << std::endl;
   for (unsigned int ib=0; ib<n_bonds; ib++) {
      const RDKit::Bond *bond_p = rdkm->getBondWithIdx(ib);
      int idx_1 = bond_p->getBeginAtomIdx();
      int idx_2 = bond_p->getEndAtomIdx();

      std::string n_1, n_2;

      const RDKit::Atom* at_1 = (*rdkm)[idx_1];
      const RDKit::Atom* at_2 = (*rdkm)[idx_2];
      try { 
         at_1->getProp("name", n_1);
         at_2->getProp("name", n_2);
      }
      catch (const KeyErrorException &err) {
      }
      catch (...) { }  // grr.
      std::string bond_type;
      if (bond_p->getBondType() == RDKit::Bond::SINGLE) bond_type   = "  single";
      if (bond_p->getBondType() == RDKit::Bond::DOUBLE) bond_type   = "  double";
      if (bond_p->getBondType() == RDKit::Bond::TRIPLE) bond_type   = "  triple";
      if (bond_p->getBondType() == RDKit::Bond::AROMATIC) bond_type = "aromatic";
      if (bond_p->getBondType() == RDKit::Bond::ONEANDAHALF) bond_type = "one-and-a-half";
      std::string bond_dir_str;
      RDKit::Bond::BondDir bond_dir = bond_p->getBondDir();
      if (bond_dir == RDKit::Bond::NONE)       bond_dir_str = "no-special";
      if (bond_dir == RDKit::Bond::BEGINWEDGE) bond_dir_str = "beginwedge";
      if (bond_dir == RDKit::Bond::BEGINDASH)  bond_dir_str = "begindash";
      if (bond_dir == RDKit::Bond::UNKNOWN)    bond_dir_str = "unknown";

      std::cout << "  " << std::setw(2) << ib << "th  "
                << std::setw(2) << idx_1 << " " << n_1 << " -- "
                << std::setw(2) << idx_2 << " " << n_2 << "  type "
                << std::setw(2) << bond_p->getBondType() << " " << bond_type << " bond-dir: "
                << bond_dir_str << std::endl;
   }
}




// now update the atom positions of the conformer iconf in
// rdkit_molecule using the atom positions in residue_p (perhaps this
// should be in rdkit-interface.hh/cc?)
//
// ignore alt confs.
void coot::update_coords(RDKit::RWMol *mol_p, int iconf, mmdb::Residue *residue_p) {

   int n_atoms;
   mmdb::PPAtom residue_atoms = NULL;
   residue_p->GetAtomTable(residue_atoms, n_atoms);
   RDKit::Conformer &conf = mol_p->getConformer(iconf);
   for (int iat=0; iat<n_atoms; iat++) {
      std::string residue_atom_name(residue_atoms[iat]->name);
      mmdb::Atom *r_at = residue_atoms[iat];
      for (int jat=0; jat<n_atoms; jat++) {
         RDKit::Atom* at_p = (*mol_p)[jat];
         try {
            std::string rdkit_atom_name;
            at_p->getProp("name", rdkit_atom_name);
            if (rdkit_atom_name == residue_atom_name) {
               RDGeom::Point3D r_pos(r_at->x, r_at->y, r_at->z);
               conf.setAtomPos(jat, r_pos);
            }
         }
         catch (const KeyErrorException &kee) {
            //             std::cout << "caught no-name for atom exception in update_coords(): "
            //                       <<  kee.what() << std::endl;
         }
      }
   }
}

// are all the bonds between the atoms (in the vector) all aromatic?
bool
coot::is_aromatic_ring(const std::vector<int> &ring_atom_indices,
                       RDKit::ROMol &rdkm) {

   bool arom = true;

   // for each atom in ring_atom_indices,
   //    find the bonds of that atom
   //    find the other atom index
   //    if the other atom index is a member of ring_atom_indices
   //       if the bond order is aromatic,
   //          continue

   for (unsigned int i=0; i<ring_atom_indices.size(); i++) {
      // presuming that ring_atom_indices[i] is not negative
      unsigned int idx(ring_atom_indices[i]);


      RDKit::Atom* at_p = rdkm[idx];
      RDKit::ROMol::OEDGE_ITER current, end;
      boost::tie(current, end) = rdkm.getAtomBonds(at_p);

      while (current != end) {

         RDKit::Bond* bond= rdkm[*current];
         int idx_other = bond->getOtherAtomIdx(idx);

         std::vector<int>::const_iterator it = std::find(ring_atom_indices.begin(),
                                                         ring_atom_indices.end(),
                                                         idx_other);
         if (it != ring_atom_indices.end()) {
            // this bond was in the ring
            if (bond->getBondType() == RDKit::Bond::AROMATIC) {
            } else {
               arom = false;
               break;
            }
         }
         current++;
      }
   }

   return arom;

}


// bond_index is the index of the bond to be replace.  atom_index
// is the index of the atom (i.e. one of the atoms in bond of
// bond_index), the R-group of which the atoms are marked for
// deletion and the atom_index atom is changed to "*" for later
// modification with various R-groups.
//
RDKit::ROMol *
coot::split_molecule(const RDKit::ROMol &mol, int bond_index, int atom_index) {

   RDKit::ROMol *ret_mol = NULL;
   RDKit::RWMol *working_mol = new RDKit::RWMol(mol);

   if (bond_index < 0 || bond_index >= int(working_mol->getNumBonds())) {
      std::cout << "split_molecule() bad bond index " << bond_index << " vs " << working_mol->getNumBonds()
                << "  in split_molecule()" << std::endl;
   } else {
      // happy path
      const RDKit::Bond *bond_p = working_mol->getBondWithIdx(bond_index);
      if (bond_p) {
         int idx_1 = bond_p->getBeginAtomIdx();
         int idx_2 = bond_p->getEndAtomIdx();
         if (atom_index == idx_2)
            std::swap(idx_1, idx_2);
         if (atom_index != idx_1) {
            std::cout << "bond index and atom index inconsistent - fail" << std::endl;
         } else {
            // OK, so far so good
            std::queue<int> q;
            std::vector<int> considered;
            std::vector<int> R_group_atoms;
            considered.push_back(idx_1); // not the picked atom;
            considered.push_back(idx_2); // not the first neighbour of idx_1.
            // what are the neighbours of idx_1 that is not idx_2?

            RDKit::Atom* at_p = (*working_mol)[idx_1];
            RDKit::ROMol::OEDGE_ITER current, end;
            boost::tie(current, end) = working_mol->getAtomBonds(at_p);

            // add some atoms to the queue
            while (current != end) {
               RDKit::Bond* bond= (*working_mol)[*current];
               int idx = bond->getOtherAtomIdx(idx_1);
               if (idx != idx_2) {
                  q.push(idx);
                  considered.push_back(idx);
                  R_group_atoms.push_back(idx);
               }
               current++;
            }

            while (q.size()) {
               int current_atom_idx = q.front();
               q.pop();
               RDKit::Atom* at_p = (*working_mol)[current_atom_idx];
               boost::tie(current, end) = working_mol->getAtomBonds(at_p);
               // std::cout << "current and end: " << current << " " << end << std::endl;
               while (current != end) {
                  RDKit::Bond* bond=(*working_mol)[*current];
                  int idx = bond->getOtherAtomIdx(current_atom_idx);
                  if (std::find(considered.begin(),
                                considered.end(),
                                idx) == considered.end()) {
                     q.push(idx);
                     considered.push_back(idx);
                     R_group_atoms.push_back(idx);
                  }
                  current++;
               }
            }

            std::cout << "R-group atoms has " << R_group_atoms.size() << " atoms"
                      << std::endl;
            if (1)
               for (unsigned int iat=0; iat<R_group_atoms.size(); iat++)
                  std::cout << "    " << R_group_atoms[iat] << std::endl;

            if (R_group_atoms.size()) {

               // make a list of atoms to deleted, then delete them
               // and set the returned molecule pointer to non-null
               //
               std::vector<RDKit::Atom *> atoms_to_be_deleted;
               for (unsigned int iat=0; iat<R_group_atoms.size(); iat++) {
                  // std::cout << "... deleting atom " << R_group_atoms[iat] << std::endl;
                  RDKit::Atom* at_p = (*working_mol)[R_group_atoms[iat]];
                  atoms_to_be_deleted.push_back(at_p);
               }
               for (unsigned int iat=0; iat<R_group_atoms.size(); iat++)
                  working_mol->removeAtom(atoms_to_be_deleted[iat]);

               ret_mol = working_mol;
            }
         }
      }
   }

   return ret_mol;
}


// atom_index is the atom in mol that is the "any" atom that is part of the fragment.
//
// This atoms is/becomes the '*' atom of the trial fragment (having atomic number of 0)
//
// Caller deletes the molecules of the vector.
//
std::vector<RDKit::ROMol *>
coot::join_molecules(const RDKit::ROMol &mol, int atom_index, const RDKit::ROMol &trial_fragment) {

   std::vector<RDKit::ROMol *> v;

   bool i_joining_atom_found = false;
   for (unsigned int iat=0; iat<trial_fragment.getNumAtoms(); iat++) {
      const RDKit::Atom* at_p = trial_fragment[iat];

      // Was it the '*' atom of the trial_fragment?
      //
      if (at_p->getAtomicNum() == 0) {
         i_joining_atom_found = true;
         std::cout << "in join_molecules() found '*' atom " << iat << std::endl;
         RDKit::RWMol working_mol = mol;
         // make a coordMap of the working mol atoms for later use in 2d coords generation
         //
         RDGeom::INT_POINT2D_MAP coordMap;
         try {
            int iconf = 0;
            int n_conf = working_mol.getNumConformers();
            std::cout << "------------ working mol has " << n_conf << " conformers" << std::endl;
            RDKit::Conformer conf = working_mol.getConformer(0);
            for (unsigned int iwk_at=0; iwk_at<working_mol.getNumAtoms(); iwk_at++) {
               RDGeom::Point3D p = conf.getAtomPos(iwk_at);
               coordMap[iwk_at] = RDGeom::Point2D(p.x, p.y);
            }
         }
         catch (...) {
            std::cout << "caught something on working mol conformer details" << std::endl;
         }

         std::map<unsigned int, unsigned int> atom_idx_map;
         std::cout << "working_mol has " << working_mol.getNumAtoms() << " atoms." << std::endl;

         //make copies of the atoms of trial_fragment and add them to working_mol
         for (unsigned int iat_inner=0; iat_inner<trial_fragment.getNumAtoms(); iat_inner++) {
            // we don't want to add the '*' atom of the trial_fragment
            if (iat_inner != iat) {
               at_p = trial_fragment[iat_inner];
               unsigned int new_n_atoms = working_mol.addAtom(at_p); // copies
               unsigned int working_idx = new_n_atoms; // I think
               atom_idx_map[iat_inner] = working_idx;
            }
         }

         if (false) {
            std::cout << "in join_molecules() iterating over " << trial_fragment.getNumBonds()
                      << " trial_fragment bonds " << std::endl;

            std::cout << "------ atom idx map ---------" << std::endl;
            std::map<unsigned int, unsigned int>::const_iterator it;
            for (it=atom_idx_map.begin(); it!=atom_idx_map.end(); ++it) {
               std::cout << "    atom_name " << it->first << " -> " << it->second << std::endl;
            }
            std::cout << "------ end of atom idx map ---------" << std::endl;
         }

         for (unsigned int ibond=0; ibond<trial_fragment.getNumBonds(); ibond++) {
            const RDKit::Bond *bond = trial_fragment.getBondWithIdx(ibond);
            unsigned int itb = bond->getBeginAtomIdx();
            unsigned int ite = bond->getEndAtomIdx();
            std::map<unsigned int, unsigned int>::const_iterator it_b = atom_idx_map.find(itb);
            std::map<unsigned int, unsigned int>::const_iterator it_e = atom_idx_map.find(ite);
            if (it_b != atom_idx_map.end()) {
               if (it_e != atom_idx_map.end()) {
                  unsigned int iwb = it_b->second;
                  unsigned int iwe = it_e->second;
                  working_mol.addBond(iwb, iwe, bond->getBondType());
               } else {
                  // OK, so ite could have been the '*' atom
                  if (ite == iat) {
                     unsigned int iwb = it_b->second;
                     working_mol.addBond(iwb, atom_index, bond->getBondType());
                  }
               }
            } else {
               // OK, so itb could have been the '*' atom
               if (itb == iat) {
                  if (it_e != atom_idx_map.end()) { // just for safety sake
                     unsigned int iwe = it_e->second;
                     working_mol.addBond(atom_index, iwe, bond->getBondType());
                  }
               }
            }
         }

         rdkit_mol_sanitize(working_mol);

         RDDepict::compute2DCoords(working_mol, &coordMap);
         RDKit::ROMol *romol = new RDKit::ROMol(working_mol);
         v.push_back(romol);
      }

   }
   return v;
}

mmdb::Residue *coot::residue_from_rdkit_mol(const RDKit::ROMol &mol, int conf_id, const std::string &new_comp_id) {

   std::cout << "DEBUG:: residue_from_rdkit_mol()::::::::::::::::::::::::::::::: --- start ---" << std::endl;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   mmdb::Residue *r = nullptr;
   std::vector<mmdb::Atom *> atoms;
   unsigned int n_atoms = mol.getNumAtoms();
   std::map<std::string, unsigned int> ele_count;
   if (n_atoms > 0) {
      unsigned int iconf = 0;
      const RDKit::Conformer &conf = mol.getConformer(iconf);
      for (unsigned int i=0; i<n_atoms; i++) {
         const RDKit::Atom *rat = mol.getAtomWithIdx(i);
         mmdb::Atom *at = new mmdb::Atom;
         int atomic_number = rat->getAtomicNum();
         std::string ele = tbl->getElementSymbol(atomic_number);
         ele_count[ele]++;
         try {
            std::string name;
            rat->getProp("name", name); // They must have been made by caller
            // this name extraction code is common code
            if (name.size() == 1) name = " " + name + "  ";
            if (name.size() == 2) name = " " + name + " ";
            if (name.size() == 3) name = " " + name;
            at->SetAtomName(name.c_str());
            at->SetElementName(ele.c_str());
            RDGeom::Point3D p = conf.getAtomPos(i);
            mmdb::realtype occ = 1.0;
            mmdb::realtype tFac = 30.0;
            std::cout << "DEBUG:: residue atom " << i << " \"" << name << "\" at " << p.x << " " << p.y << " " << p.z << std::endl;
            at->SetCoordinates(p.x, p.y, p.z, occ, tFac);
            atoms.push_back(at);
         }
         catch (const std::runtime_error &rte) {
            std::cout << "DEBUG:: rte " << rte.what() << std::endl;
         }
      }
      if (! atoms.empty()) {
         r = new mmdb::Residue;
         r->SetResID(new_comp_id.c_str(), 1, "");
         for (unsigned int i=0; i<atoms.size(); i++) {
            mmdb::Atom *at = atoms[i];
            r->AddAtom(at);
         }
      }
   }
   std::cout << "DEBUG:: residue_from_rdkit_mol()::::::::::::::::::::::::::::::: returning " << r << std::endl;
   return r;
}

std::optional<coot::dictionary_residue_restraints_t> coot::dictionary_from_rdkit_mol(const RDKit::ROMol &mol, int conf_id, const std::string &new_comp_id) {

   // 2025-12-28-PE - Hmmm.. I must have done something like this for pyrogen.

   auto make_bond = [] (const RDKit::ROMol &mol, const RDKit::Bond *bond_ptr,
                        const std::map<std::string, std::string> &name_map, unsigned int iconf) -> std::optional<coot::dict_bond_restraint_t> {

      const RDKit::Conformer &conf = mol.getConformer(iconf);
      unsigned int idx_1  = bond_ptr->getBeginAtomIdx();
      unsigned int idx_2  = bond_ptr->getEndAtomIdx();
      RDKit::Atom *atom_1 = bond_ptr->getBeginAtom();
      RDKit::Atom *atom_2 = bond_ptr->getEndAtom();
      std::string name_1;
      std::string name_2;
      atom_1->getProp("name", name_1);
      atom_2->getProp("name", name_2);
      std::map<std::string, std::string>::const_iterator it_1;
      std::map<std::string, std::string>::const_iterator it_2;
      it_1 = name_map.find(name_1);
      it_2 = name_map.find(name_2);
      if (it_1 != name_map.end() &&  it_2 != name_map.end()) {
         std::string type = "x-unset-x";
         if (bond_ptr->getBondType() == RDKit::Bond::SINGLE)   type = "single";
         if (bond_ptr->getBondType() == RDKit::Bond::DOUBLE)   type = "double";
         if (bond_ptr->getBondType() == RDKit::Bond::TRIPLE)   type = "triple";
         if (bond_ptr->getBondType() == RDKit::Bond::AROMATIC) type = "aromatic";
         RDGeom::Point3D p1 = conf.getAtomPos(idx_1);
         RDGeom::Point3D p2 = conf.getAtomPos(idx_2);
         RDGeom::Point3D delta = p2 - p1;
         double dd = delta.lengthSq();
         double d = std::sqrt(dd);
         double esd = 0.015;
         std::cout << "DEBUG:: dictionary_from_rdkit_mol(): bond \"" << name_1 << "\" \"" << name_2
                   << "\" type: " << type << " rdkit-type: " << bond_ptr->getBondType() << " bl: " << d << std::endl;
         dict_bond_restraint_t b(name_1, name_2, type, d, esd, 0.0, 0.0, false);
         return b;
      } else {
         std::cout << "DEBUG:: bond atom name lookup fail" << std::endl;
         return std::nullopt;
      }
   };

   std::cout << "DEBUG:: dict_from_rdkit_mol()::::::::::::::::::::::::::::::: --- start ---" << std::endl;
   std::map<std::string, std::string> name_map;
   coot::dictionary_residue_restraints_t dict;
   dict.residue_info.name    = new_comp_id;
   dict.residue_info.comp_id = new_comp_id;
   dict.residue_info.group = "non-polymer";
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   std::vector<coot::dict_atom> atoms;
   unsigned int n_atoms = mol.getNumAtoms();
   std::map<std::string, unsigned int> ele_count;
   if (n_atoms > 0) {
      unsigned int iconf = 0;
      const RDKit::Conformer &conf = mol.getConformer(iconf);
      for (unsigned int i=0; i<n_atoms; i++) {
         const RDKit::Atom *rat = mol.getAtomWithIdx(i);
         int atomic_number = rat->getAtomicNum();
         std::string ele = tbl->getElementSymbol(atomic_number);
         ele_count[ele]++;
         std::string name;
         try {
            rat->getProp("name", name); // They must have not been make by the caller
            std::string name_4c = name;
            if (name.size() == 1) name_4c = " " + name + "  ";
            if (name.size() == 2) name_4c = " " + name + " ";
            if (name.size() == 3) name_4c = " " + name;
            RDGeom::Point3D p = conf.getAtomPos(i);
            std::pair<bool, float> charge(false, 0.0);
            std::string type_energy = "missing-energy-type";  // FIXME
            rat->getProp("type_energy", type_energy);
            name_map[name] = name_4c;
            dict_atom at(name, name_4c, ele, type_energy, charge);
            int pos_type = 1;
            std::pair<bool, clipper::Coord_orth> cpos(true, clipper::Coord_orth(p.x, p.y, p.z));
            std::cout << "DEBUG:: dictionary atom " << i << " \"" << name << "\" type " << type_energy << " at " << cpos.second.format() << std::endl;
            at.add_pos(pos_type, cpos);
            at.add_ordinal_id(i);
            atoms.push_back(at);
         } catch (const std::runtime_error &rte) {
            std::cout << "DEBUG:: missing atom name for atom index " << i << std::endl;
            return std::nullopt;
         }
      }
      std::vector<dict_bond_restraint_t> bonds;
      unsigned int n_bonds = mol.getNumBonds();
      for (unsigned int i=0; i<n_bonds; i++) {
         const RDKit::Bond *bond_ptr = mol.getBondWithIdx(i);
         std::optional<coot::dict_bond_restraint_t> bond = make_bond(mol, bond_ptr, name_map, iconf);
         if (bond)
            bonds.push_back(bond.value());
      }
      dict.atom_info = atoms;
      dict.bond_restraint = bonds;
      return dict;
   } else {
      return std::nullopt;
   }
}



#endif // MAKE_ENHANCED_LIGAND_TOOLS
