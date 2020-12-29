/* pyrogen/mmff-restraints.cc
 * 
 * Copyright 2014 by Medical Research Council
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

//#include "lidia-core/use-rdkit.hh"
#include "mmff-restraints.hh"
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/BondStretch.h>
#include <ForceField/MMFF/AngleBend.h>

//#include "mmff-restraints.hh"
#include "restraints-private.hh" // for bond-order conversion

// This function can potentially alter the aromaticity of the atoms of mol
// 
// caller disposes of returned object
coot::mmff_b_a_restraints_container_t *
coot::mmff_bonds_and_angles(RDKit::ROMol &mol) {

   mmff_b_a_restraints_container_t *r = new mmff_b_a_restraints_container_t;

   RDKit::MMFF::MMFFMolProperties *mmffMolProperties = new RDKit::MMFF::MMFFMolProperties(mol);
   if (! mmffMolProperties->isValid()) {
      std::cout << "invalid properties " << std::endl;
   } else {
      // happy path

      // iterate over bonds - simple
      // 
      // ForceFields::MMFF::MMFFBondCollection *mmff_bonds =
      // ForceFields::MMFF::MMFFBondCollection::getMMFFBond();
      
      RDKit::ROMol::BondIterator bondIt;
      RDKit::ROMol::BondIterator start;
      RDKit::ROMol::BondIterator end;

      for (bondIt=mol.beginBonds(); bondIt!=mol.endBonds(); bondIt++) {
      	 unsigned int idx_1 = (*bondIt)->getBeginAtomIdx();
      	 unsigned int idx_2 = (*bondIt)->getEndAtomIdx();
      	 unsigned int iAtomType_1 = mmffMolProperties->getMMFFAtomType(idx_1);
      	 unsigned int iAtomType_2 = mmffMolProperties->getMMFFAtomType(idx_2);
      	 unsigned int bondType  = mmffMolProperties->getMMFFBondType(*bondIt);
      	 const ForceFields::MMFF::MMFFBond *mmffBondParams = 0;
         // (*mmff)(bondType, iAtomType_1, iAtomType_2);  // I don't know how to look up the bond list now.
         // FIXME
      	 if (mmffBondParams) { 
      	    double r0 = ForceFields::MMFF::Utils::calcBondRestLength(mmffBondParams);
      	    double kb = ForceFields::MMFF::Utils::calcBondForceConstant(mmffBondParams);
	    double sigma = 0.02;
	    if (kb != 0.0) // the usual case, I hope!
	       sigma = 0.04/sqrt(kb);
      	    std::string order = convert_to_energy_lib_bond_type((*bondIt)->getBondType());
      	    mmff_bond_restraint_info_t br(idx_1, idx_2, order, r0, sigma);
      	    r->bonds.push_back(br);
      	 }
      }
      
      // iterate over angles
      // 
      ForceFields::MMFF::MMFFAngleCollection *mmff_angles;
      // = ForceFields::MMFF::MMFFAngleCollection::getMMFFAngle();  // I don't know how to look up the anglelist now.
      // FIXME
      unsigned int n_atoms = mol.getNumAtoms();
      std::map<unsigned long long, bool> done_angle;
      for (unsigned int iat_1=0; iat_1<n_atoms; iat_1++) { 
	 const RDKit::Atom *at_1 = mol[iat_1];
	 RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
	 boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_1);
	 while(nbr_idx_1 != end_nbrs_1){
	    const RDKit::Atom *at_2 = mol[*nbr_idx_1];

	    RDKit::ROMol::ADJ_ITER nbr_idx_2, end_nbrs_2;
	    boost::tie(nbr_idx_2, end_nbrs_2) = mol.getAtomNeighbors(at_2);
	    while(nbr_idx_2 != end_nbrs_2){
	       const RDKit::Atom *at_3 = mol[*nbr_idx_2];
	       if (at_3 != at_1) {

		  unsigned int idx_1 = at_1->getIdx();
		  unsigned int idx_2 = at_2->getIdx();
		  unsigned int idx_3 = at_3->getIdx();

		  unsigned int m = 10000;
		  unsigned long long angle_key_1 = idx_1 * m * m + idx_2 * m + idx_3;
		  unsigned long long angle_key_2 = idx_3 * m * m + idx_2 * m + idx_1;

		  if (done_angle.find(angle_key_1) == done_angle.end() &&
		      done_angle.find(angle_key_2) == done_angle.end()) {

		     done_angle[m] = true;

		     unsigned int iAtomType_1 = mmffMolProperties->getMMFFAtomType(idx_1);
		     unsigned int iAtomType_2 = mmffMolProperties->getMMFFAtomType(idx_2);
		     unsigned int iAtomType_3 = mmffMolProperties->getMMFFAtomType(idx_3);

		     // unsigned int angle_type =
                     // mmffMolProperties->getMMFFAngleType(mol, idx_1, idx_2, idx_3);

 		     const ForceFields::MMFF::MMFFAngle *mmffAngleParams = 0;

                     // doesn't compile now
                     // (*mmff_angles)(angle_type, iAtomType_1, iAtomType_2, iAtomType_3);
		     
		     if (mmffAngleParams) {
			double a = ForceFields::MMFF::Utils::calcAngleRestValue(mmffAngleParams);
			double k = ForceFields::MMFF::Utils::calcAngleForceConstant(mmffAngleParams);

			double esd = 4.0;
			if (k != 0) // the usual case, I hope.
			   esd = 3.0/sqrt(k);

			if (false)
			   std::cout << "mmff angle restraint " << idx_1 << " " << idx_2 << " " << idx_3 << "    "
				     << a << " k: " << k << " esd: " << esd << std::endl;
			mmff_angle_restraint_info_t angle(idx_1, idx_2, idx_3, a, esd);
			r->angles.push_back(angle);
		     }
		  }
	       }
	       nbr_idx_2++;
	    }
	    nbr_idx_1++;
	 }
      }
   }
   return r;
}

// These functions can potentially alter the aromaticity of the atoms of mol
coot::dictionary_residue_restraints_t
coot::make_mmff_restraints(RDKit::ROMol &mol) {
   
   dictionary_residue_restraints_t r;
   mmff_b_a_restraints_container_t *mm_info = mmff_bonds_and_angles(mol); // d

   for (unsigned int ibond=0; ibond<mm_info->bonds_size(); ibond++) { 
      try {
	 unsigned int idx_1 = mm_info->bonds[ibond].get_idx_1();
	 unsigned int idx_2 = mm_info->bonds[ibond].get_idx_2();
	 const RDKit::Atom *at_p_1 = mol[idx_1];
	 const RDKit::Atom *at_p_2 = mol[idx_2];
	 std::string name_1 = "";
	 std::string name_2 = "";
	 at_p_1->getProp("name", name_1);
	 at_p_2->getProp("name", name_2);

	 dict_bond_restraint_t br(name_1, name_2,
				  mm_info->bonds[ibond].get_type(),
				  mm_info->bonds[ibond].get_resting_bond_length(),
				  mm_info->bonds[ibond].get_sigma());
	 r.bond_restraint.push_back(br);
      }
      catch (const KeyErrorException &kee) {
	 std::cout << "ERROR:: OOps in make_mmff_restraints(): " << kee.what() << std::endl;
      }
   }
   for (unsigned int iangle=0; iangle<mm_info->angles_size(); iangle++) { 
      try {
	 unsigned int idx_1 = mm_info->angles[iangle].get_idx_1();
	 unsigned int idx_2 = mm_info->angles[iangle].get_idx_2();
	 unsigned int idx_3 = mm_info->angles[iangle].get_idx_3();
	 const RDKit::Atom *at_p_1 = mol[idx_1];
	 const RDKit::Atom *at_p_2 = mol[idx_2];
	 const RDKit::Atom *at_p_3 = mol[idx_3];
	 std::string name_1 = "";
	 std::string name_2 = "";
	 std::string name_3 = "";
	 at_p_1->getProp("name", name_1);
	 at_p_2->getProp("name", name_2);
	 at_p_3->getProp("name", name_3);

	 dict_angle_restraint_t br(name_1, name_2, name_3,
				   mm_info->angles[iangle].get_resting_angle(),
				   mm_info->angles[iangle].get_sigma());
	 r.angle_restraint.push_back(br);
      }
      catch (const KeyErrorException &kee) {
	 std::cout << "ERROR:: OOps in make_mmff_restraints(): " << kee.what() << std::endl;
      }
   }
   


   delete mm_info;
   return r;
} 
