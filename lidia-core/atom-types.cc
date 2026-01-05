/*
 * lidia-core/atom-types.cc
 *
 * Copyright 2017 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include <cctype>
#include "utils/coot-utils.hh"
#include "rdkit-interface.hh"


// add "energy_type" properties to the atoms
void
coot::set_energy_lib_atom_types(RDKit::ROMol *mol) {

   class idx_pair {
   public:
      unsigned int idx_0;
      unsigned int idx_1;
      idx_pair(unsigned int idx_0_in,
	       unsigned int idx_1_in) :
	 idx_0(idx_0_in),
	 idx_1(idx_1_in) {}
   };
   class atom_typing {
   public:
      std::vector<unsigned int> atom_index;
      std::string energy_type;
      std::string smarts_string;
      atom_typing() {}
      atom_typing(const std::string &energy_type_in,
		  const std::string &smarts_string_in,
		  unsigned int atom_index_in) :
	 energy_type(energy_type_in),
	 smarts_string(smarts_string_in) {
	 atom_index.push_back(atom_index_in);
      }
      atom_typing(const std::string &energy_type_in,
		  const std::string &smarts_string_in,
		  idx_pair atom_indices_in) :
	 energy_type(energy_type_in),
	 smarts_string(smarts_string_in) {
	 atom_index.push_back(atom_indices_in.idx_0);
	 atom_index.push_back(atom_indices_in.idx_1);
      }
   };

   std::vector<atom_typing> smarts_list = {

        // Full coverage for C, H, O.

        // Oxygen
        atom_typing("O2",  "[OX2;H0]", 0), // ester, Os between P and C are O2, not OP
        atom_typing("OP",  "O~P",   0),
        atom_typing("OS",  "O~S",   0),
        atom_typing("OB",  "O~B",   0),
        atom_typing("OC",  "*C(=O)[OH]", idx_pair(2,3)), // carboxylic acid
        atom_typing("OC",  "*C(=O)O",    idx_pair(2,3)), // carboxylate, doesn"t match deloc bonds
        atom_typing("OH1", "[OH1]", 0), // alcohol
        atom_typing("O2",  "[oX2;H0]", 0), // ring oxygen
        atom_typing("O",   "O=*",   0), // carbonyl oxygen

        // OH2 no examples
        // OHA no examples
        // OHB no examples
        // OHC no examples
        // OC2 no exmampes
        
        // Fallback oxygen
        atom_typing("O",   "O",   0),

        // Carbon SP
        {"CSP1", "[H][C]#*",  1}, // e.g. in 2GT
        {"CSP",  "[C]#[C]",   0},
        {"CSP",  "[C]#[C]",   1},
        {"CSP",  "[C]#*",     0},

        // Carbon SP2
        {"CR56", "c12aaaac1aaa2",  0}, // works on indole
        {"CR56", "c12aaaac1aaa2",  5}, // works on indole
        {"CR56", "c12aaaan1aaa2",  0}, // same pattern as (below) N56, but catching first 56 atom
        {"CR56", "c12AAAAc1aaa2",  0}, // 6-ring doesn"t have to aromatic (8UG)
        {"CR66", "c12aaaac1aaaa2", 0},
        {"CR66", "c12aaaac1aaaa2", 5},
        {"CR6",  "c12ccccc1OCO2",  0},   // mouse, fused atoms in 6-ring not non-arom 5-ring
        {"CR6",  "c12ccccc1OCO2",  5},   // mouse, fused atoms in 6-ring not non-arom 5-ring
        {"CR66", "c12aaaac1AAAA2", 0},   // one 6-ring aromatic, other not. Needed for XXX?
        {"CR66", "c12aaaac1AAAA2", 5},   // one 6-ring aromatic, other not. Needed for XXX?
        // but makes a fail on 113.
        {"CR6",  "c12caccc1***2",  0},  // aromatic 6, (non-)aromatic 5, maybe this should be CR56?
        {"CR6",  "c12caccc1***2",  5},  // aromatic 6, (non-)aromatic 5, maybe this should be CR56?

        // note CR1  missing - can"t find example
        //      CR1H missing - can"t find example

        atom_typing("CR16", "[cr6;H1]",  0),
        atom_typing("CR6",  "[cr6;H0]",  0),
        atom_typing("CR15", "[cr5;H1]",  0),
//        ("CR5",  "C1(=O)[C,c][C,c]C=N1", 0), // carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
        atom_typing("CR5",  "[cr5;H0]",  0),
        atom_typing("CR5",  "[CR5;H0]",  0),
        atom_typing("C1",   "[CX3;H1]",    0),  // double bond, single bond and one H
        atom_typing("C2",   "[CX3;H2]=*",  0),  // double bond, and 2 H
        atom_typing("C",    "[CX3;H0;^2]", 0),
        atom_typing("C",    "[CX3]=[OX1]", 0),  // carbonyl carbon
        atom_typing("C",    "[$([CX2](=C)=C)]",   0), // bonded to 3 things not hydrogen

        // Carbon SP3
        atom_typing("CT",   "[CX4H0]", 0), // single bonded to 4 things not hydrogen
        atom_typing("CH3",  "[C;H3;^3]",   0), // bonded to something via single bond and 3 Hs
        atom_typing("CH2",  "[C;^3;H2]",   0), // standard aliphatic C.
        atom_typing("CH1",  "*[C;H1](*)*", 1), // bonded to H and 3 things 

        // sp??? needs sorting 
        atom_typing("CH2",  "[CH2]",   0), // bonded to 2 hydrogens
        
        // Carbon fallback
        atom_typing("C", "[C,c]", 0),

        // Hydrogen
        atom_typing("HCH1", "[H][CH1]",    0),
        atom_typing("HCH2", "[H][C;H2^3]", 0),
        atom_typing("HCH3", "[H][CH3]",    0),
        atom_typing("HNC1", "[H][N;H1;^2]~C(~N)~N", 0), // H of N of N=C ? 
        atom_typing("HNC2", "[H][NX3;H2;^2]", 0), // H on a NC2 (NH1 and NH2 of ARG)
        atom_typing("HNC3", "[H][NX3;H3;^2]", 0), // guess - no examples
        atom_typing("HNT1", "[H][NX4;H1;^3]", 0),
        atom_typing("HNT1", "[H][NX3;H1;^3]", 0),
        atom_typing("HNT2", "[H][NX3;H2;^3]", 0), // H connected to type NT2
        atom_typing("HNT3", "[N^3;H3][H]", 1), // NH3+ 
        atom_typing("HNH2", "[H][NH2;^2]", 0), // NH2 is sp2
        atom_typing("HNH1", "[H][NX3;H1;^2]",    0),
        atom_typing("HCR6", "[H][cr6;H1]", 0),
        atom_typing("HCR5", "[H][cr5;H1]", 0), // connected to aromatic ring C with 1 H
        atom_typing("HNR5", "[H][nr5;H1]", 0), // connected to aromatic ring C with 1 H
        atom_typing("HNR5", "[H][Nr5;H1]", 0), // guess based on above, connected to aromatic N in a 5-ring
        atom_typing("HNR6", "[H][nr6;H1]", 0), // connected to aromatic 6-ring N with 1 H
        atom_typing("HNR6", "[H][NR6;H1]", 0), // guess based on above

        // HCR missing - no examples (and how is it different to HCR1?)
        atom_typing("HCR1", "[H]c",        0),
        atom_typing("HNH1", "[H][NH1]",    0),
        atom_typing("HOH1", "[H][OH1]",    0),
        atom_typing("HOH2", "[H][OH2][H]", 0), // H of HOH - water
        atom_typing("HOH2", "[H][OH2][H]", 2), // H of HOH - water
        atom_typing("H",    "[H]",         0),

        // Nitrogen, SP3

        atom_typing("NT1", "[NX4;H1;^3]",  0),
        atom_typing("NT1", "[NX3;H1;^3]",  0),
        atom_typing("NT2", "[NX3;H2;^3]",  0), // different to mon-lib!
        atom_typing("NT3", "[NX4;H3;^3]",  0),
        atom_typing("NT",  "[NX3;H0;^3]",  0),


        // NE-CZ in ARG should be deloc (guandino) - also NC1-C
        // single (as is in ARG.cif) is not found in ener_lib!
       
        // Nitrogen, SP2
        atom_typing("NR66", "c12aaaan1aaaa2", 5), // (second) 66 atom is an N.
        atom_typing("NR56", "c12aaaan1aaa2",  5), // (second) 56 atom is an N.
        atom_typing("NR55", "c12aaan1aaa2",   4), // (second) 55 atom is an N.
        atom_typing("NC2",  "[NX3;H2^2]", 0),     // N of sp2 NH2 (as in ARG).
        atom_typing("NH2",  "[NX3^2][CX3^2]=[N^2;X3+]", 0), // amidinium (charged)... 
        atom_typing("NH2",  "[NX3^2][CX3^2]=[N^2;X3+]", 2), // amidinium (charged)... 
        atom_typing("NR15", "[nr5;X3;H1]",    0),
        atom_typing("NR5",  "[nr5;X3;H0]",    0),
        atom_typing("NR5",  "[NR;X3;H0;^2]",  0), // [NR5;X3;H0;^2] fails on 14C (also at daylight)
        atom_typing("NRD5", "[nr5;X2;H0]",    0), // guess from 071
        atom_typing("NRD5", "C1(=O)[C,c][C,c]C=N1", 5), // N bonded to carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
        atom_typing("NR16", "[nr6;H1]",    0),
        atom_typing("NRD6", "a:[nr6;X2;H0]:a",  1), // aromatic N with no H, i.e. one double one single
        atom_typing("NR6",  "[nr6]",    0),
        atom_typing("NC1",  "[H][N;H1;^2]~C(~N)~N", 1), 
        atom_typing("NC1",  "[NX3;H1;^2]C(~N)~N", 0), // N, as in NE in ARG
        atom_typing("NC1",  "[NX2;H1;^2]", 0),  // N of N=C ? 
        atom_typing("NH1",  "[NX3;H1;^2]", 0),
        atom_typing("NH2",  "[NX3;H2;^2]", 0),  // sp2, e.g. ND2 of an ASP
        atom_typing("NT",   "*n1~[o]~[o]1", 1), // guess from 16X dioxaziridine (bleugh)
        // (NT needs checking?)
        // NC2 no examples
        // NC3 no examples
        // NPA no examples
        // NPB no examples
        

        // Nitrogen SP1
        atom_typing("NS",   "[N^1]", 0),
        // NS1 no examples
        

        // fall-back nitrogen
        atom_typing("N",    "[N,n]",      0),  

        // Phosphorus
        atom_typing("P",    "P", 0),
        // Cl
        atom_typing("CL",   "[Cl]", 0),
        // F
        atom_typing("F",    "[F]",  0),
        // Br
        atom_typing("BR",    "[Br]",  0),

        // Sulfur
        atom_typing("SH1",  "[SX2H1]", 0),  // SG of CYS
        atom_typing("ST",   "[SX4]", 0),  // tetrahedral (2 single bonds, 2 double)
        atom_typing("S1",   "[S]=*", 0),
        atom_typing("S2",   "[SX2,sX2]", 0),
        atom_typing("S3",   "[SX3,sX3]", 0),
        atom_typing("S",    "[S,s]", 0),

        // Silicon
        atom_typing("SI1",  "[Si;X4]", 0), // tetragonal Si
        atom_typing("SI",   "[Si]",    0)  // Si any other

   };


   std::vector<std::string> eles = {
      "He", "Li", "Be", "B",  "Ne", "Na", "Mg", "Al", 
      "Ar", "K", "Ca",  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
      "Zn", "Ga", "Ge", "As", "Se",      "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
      "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
      "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
      "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
      "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"};


   std::vector<atom_typing> ele_smarts(eles.size());
   for (std::size_t i=0; i<eles.size(); i++)
      ele_smarts[i] = atom_typing(util::upcase(eles[i]), "["+eles[i]+"]", 0);

   smarts_list.insert(smarts_list.end(),
		      ele_smarts.begin(),
		      ele_smarts.end());

   for (std::size_t ism=0; ism<smarts_list.size(); ism++) {
      const atom_typing &smarts_type = smarts_list[ism];
      RDKit::RWMol *query = RDKit::SmartsToMol(smarts_type.smarts_string);
      std::vector<RDKit::MatchVectType> matches;
      bool recursionPossible = true;
      bool useChirality = true;
      bool uniquify = true;
      unsigned int n_matched = RDKit::SubstructMatch(*mol, *query, matches, uniquify, recursionPossible, useChirality); // 20210923-PE FIXME
      // int matched = false;
      //
      if (n_matched > 0) {
	 for (unsigned int im=0; im<matches.size(); im++) {
	    for (std::size_t j=0; j<smarts_type.atom_index.size(); j++) {
	       unsigned int match_atom_index = smarts_type.atom_index[j];
	       unsigned int idx_this_atom_1 = matches[im][match_atom_index].first;
	       unsigned int idx_this_atom_2 = matches[im][match_atom_index].second;
	       if (true)
		  std::cout << "query " << smarts_list[ism].smarts_string
			    << " matches idx pair " << idx_this_atom_1 << " " << idx_this_atom_2
			    << " " << smarts_type.energy_type
			    << std::endl;
	       RDKit::Atom *at_p = mol->getAtomWithIdx(idx_this_atom_2);
	       // if this atom has a energy_type already, pass, else set it
	       try {
		  std::string e;
		  at_p->getProp("type_energy", e);
		  std::cout << "already has type_energy \"" << e << "\""<< std::endl;
	       }
	       catch (const KeyErrorException &e) {
		  std::cout << "setting type_energy " << smarts_type.energy_type
		            << " for atom " << idx_this_atom_2 << std::endl;
		  at_p->setProp("type_energy", smarts_type.energy_type);
	       }
	    }
	 }
      }
   }
}

// dicitionaries from CCDs don't have energy atom types. We need them for ligand
// environment analysis (flev).  It is presumed that mol has energy_type atom types
// (e.g. set from the above function). mol is not modified.
//
void
coot::set_dictionary_atom_types_from_mol(dictionary_residue_restraints_t *dictionary,
					 const RDKit::ROMol *mol) {

   unsigned int n_mol_atoms = mol->getNumAtoms();

   if (false)
      std::cout << "here in set_dictionary_atom_types_from_mol() with "
		<< n_mol_atoms << " atoms" << std::endl;

   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      const RDKit::Atom *at_p = mol->getAtomWithIdx(iat);
      try {
	 std::string type;
	 std::string atom_id;
	 at_p->getProp("type_energy", type);
	 // std::cout << "debug:: here with type " << type << std::endl;
	 // at_p->getProp("atom_id", atom_id);
	 at_p->getProp("name", atom_id); // should be atom_id
	 for (std::size_t j=0; j<dictionary->atom_info.size(); j++) {
	    if (dictionary->atom_info[j].atom_id_4c == atom_id) {
	       dictionary->atom_info[j].type_energy = type;
	       // std::cout << "debug:: for atom " << atom_id << " added type " << type << std::endl;
	       break;
	    }
	 }
      }
      catch (const KeyErrorException &ke) {
	 std::cout << "WARNING:: " << ke.what() << " in set_dictionary_atom_types_from_mol() for atom index "
		   << iat << std::endl;
      }
   }
}

// make an rdkit molecule from dictionary and use the above function to the energy types
void
coot::set_dictionary_atom_types(dictionary_residue_restraints_t *dictionary) {

   RDKit::RWMol m = rdkit_mol(*dictionary);
   set_dictionary_atom_types_from_mol(dictionary, &m);

}

#endif // MAKE_ENHANCED_LIGAND_TOOLS
