
#include <vector>
#include <tuple>
#include <string>
#include <GraphMol/GraphMol.h> // Include RDKit for handling molecules and SMARTS patterns
#include <GraphMol/Substruct/SubstructMatch.h>

#include "lidia-core/rdkit-interface.hh" // something here is the right header for SmartsToMol().

   
std::vector<std::tuple<std::string, std::string, int> > get_smarts_by_element() {

   auto ele_to_smarts = [] (const std::string &ele) {
      return std::tuple<std::string, std::string, int> (ele, "[" + ele + "]", 0);
   };

   std::vector<std::string> eles = {
      "He", "Li", "Be", "B",  "Ne", "Na", "Mg", "Al",
      "Ar", "K", "Ca",  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
      "Zn", "Ga", "Ge", "As", "Se",      "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
      "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
      "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
      "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
      "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"};

   std::vector<std::tuple<std::string, std::string, int> > smarts(eles.size());
   for (unsigned int i=0; i<eles.size(); i++)
      smarts[i] = ele_to_smarts(eles[i]);
   return smarts;

}


void set_monomer_library_atom_types(const RDKit::ROMol& mol) {

   std::vector<std::tuple<std::string, std::string, int> > smarts_list = {

      {"O2",  "[OX2;H0]", 0},
      {"OP",  "O~P",   0},
      {"OS",  "O~S",   0},
      {"OB",  "O~B",   0},
      {"OC",  "*C(=O)[OH]", 2}, // carboxylic acid
      {"OC",  "*C(=O)[OH]", 3}, // carboxylic acid
      {"OC",  "*C(=O)O",    2}, // carboxylate, doesn"t match deloc bonds
      {"OC",  "*C(=O)O",    3}, // carboxylate, doesn"t match deloc bonds
      {"OH1", "[OH1]",    0}, // alcohol
      {"O2",  "[oX2;H0]", 0}, // ring oxygen
      {"O",   "O=*",      0}, // carbonyl oxygen

      // OH2 no examples
      // OHA no examples
      // OHB no examples
      // OHC no examples
      // OC2 no exmampes

      // Fallback oxygen
      {"O",   "O",   0},

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

      {"CR16", "[cr6;H1]",  0},
      {"CR6",  "[cr6;H0]",  0},
      {"CR15", "[cr5;H1]",  0},
      // ("CR5",  "C1(=O)[C,c][C,c]C=N1", 0), # carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
      {"CR5",  "[cr5;H0]",  0},
      {"CR5",  "[CR5;H0]",  0},
      {"C1",   "[CX3;H1]",    0},  // double bond, single bond and one H
      {"C2",   "[CX3;H2]=*",  0},  // double bond, and 2 H
      {"C",    "[CX3;H0;^2]", 0},
      {"C",    "[CX3]=[OX1]", 0},  // carbonyl carbon
      {"C",    "[$([CX2](=C)=C)]", 0}, // bonded to 3 things not hydrogen

      // Carbon SP3
      {"CT",   "[CX4H0]",     0}, // single bonded to 4 things not hydrogen
      {"CH3",  "[C;H3;^3]",   0}, // bonded to something via single bond and 3 Hs
      {"CH2",  "[C;^3;H2]",   0}, // standard aliphatic C.
      {"CH1",  "*[C;H1](*)*", 1}, // bonded to H and 3 things

      // sp??? needs sorting
      {"CH2",  "[CH2]",   0}, // bonded to 2 hydrogens

      // Carbon fallback
      {"C", "[C,c]", 0},

      // Hydrogen
      {"HCH1", "[H][CH1]",    0},
      {"HCH2", "[H][C;H2^3]", 0},
      {"HCH3", "[H][CH3]",    0},
      {"HNC1", "[H][N;H1;^2]~C(~N)~N", 0}, // H of N of N=C ?
      {"HNC2", "[H][NX3;H2;^2]", 0}, // H on a NC2 (NH1 and NH2 of ARG)
      {"HNC3", "[H][NX3;H3;^2]", 0}, // guess - no examples
      {"HNT1", "[H][NX4;H1;^3]", 0},
      {"HNT1", "[H][NX3;H1;^3]", 0},
      {"HNT2", "[H][NX3;H2;^3]", 0}, // H connected to type NT2
      {"HNT3", "[N^3;H3][H]", 1}, // NH3+
      {"HNH2", "[H][NH2;^2]", 0}, // NH2 is sp2
      {"HNH1", "[H][NX3;H1;^2]",    0},
      {"HCR6", "[H][cr6;H1]", 0},
      {"HCR5", "[H][cr5;H1]", 0}, // connected to aromatic ring C with 1 H
      {"HNR5", "[H][nr5;H1]", 0}, // connected to aromatic ring C with 1 H
      {"HNR5", "[H][Nr5;H1]", 0}, // guess based on above, connected to aromatic N in a 5-ring
      {"HNR6", "[H][nr6;H1]", 0}, // connected to aromatic 6-ring N with 1 H
      {"HNR6", "[H][NR6;H1]", 0}, // guess based on above

      // HCR missing - no examples (and how is it different to HCR1?)
      {"HCR1", "[H]c",        0},
      {"HNH1", "[H][NH1]",    0},
      {"HOH1", "[H][OH1]",    0},
      {"HOH2", "[H][OH2][H]", 0}, // H of HOH - water
      {"HOH2", "[H][OH2][H]", 2}, // H of HOH - water
      {"H",    "[H]",         0},

      // Nitrogen, SP3

      {"NT1", "[NX4;H1;^3]",  0},
      {"NT1", "[NX3;H1;^3]",  0},
      {"NT2", "[NX3;H2;^3]",  0}, // different to mon-lib!
      {"NT3", "[NX4;H3;^3]",  0},
      {"NT",  "[NX3;H0;^3]",  0},


      // NE-CZ in ARG should be deloc (guandino) - also NC1-C
      // single (as is in ARG.cif) is not found in ener_lib!

      // Nitrogen, SP2
      {"NR66", "c12aaaan1aaaa2", 5}, // (second} 66 atom is an N.
      {"NR56", "c12aaaan1aaa2",  5}, // (second) 56 atom is an N.
      {"NR55", "c12aaan1aaa2",   4}, // (second) 55 atom is an N.
      {"NC2",  "[NX3;H2^2]", 0},     // N of sp2 NH2 (as in ARG).
      {"NH2",  "[NX3^2][CX3^2]=[N^2;X3+]", 0}, // amidinium (charged)...
      {"NH2",  "[NX3^2][CX3^2]=[N^2;X3+]", 2}, // amidinium (charged)...
      {"NR15", "[nr5;X3;H1]",    0},
      {"NR5",  "[nr5;X3;H0]",    0},
      {"NR5",  "[NR;X3;H0;^2]",  0}, // [NR5;X3;H0;^2] fails on 14C (also at daylight)
      {"NRD5", "[nr5;X2;H0]",    0}, // guess from 071
      {"NRD5", "C1(=O)[C,c][C,c]C=N1", 5}, // N bonded to carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
      {"NR16", "[nr6;H1]",    0},
      {"NRD6", "a:[nr6;X2;H0]:a",  1}, // aromatic N with no H, i.e. one double one single
      {"NR6",  "[nr6]",    0},
      {"NC1",  "[H][N;H1;^2]~C(~N)~N", 1},
      {"NC1",  "[NX3;H1;^2]C(~N)~N", 0}, // N, as in NE in ARG
      {"NC1",  "[NX2;H1;^2]",  0},  // N of N=C ?
      {"NH1",  "[NX3;H1;^2]",  0},
      {"NH2",  "[NX3;H2;^2]",  0},  // sp2, e.g. ND2 of an ASP
      {"NT",   "*n1~[o]~[o]1", 1}, // guess from 16X dioxaziridine (bleugh)
      // (NT needs checking?)
      // NC2 no examples
      // NC3 no examples
      // NPA no examples
      // NPB no examples


      // Nitrogen SP1
      {"NS",   "[N^1]", 0},
      // NS1 no examples


      // fall-back nitrogen
      {"N",    "[N,n]",      0},

      // Phosphorus
      {"P",    "P", 0},
      // Cl
      {"CL",   "[Cl]", 0},
      // F
      {"F",    "[F]",  0},
      // Br
      {"BR",    "[Br]",  0},

      // Sulfur
      {"SH1",  "[SX2H1]", 0},  // SG of CYS
      {"ST",   "[SX4]",   0},  // tetrahedral (2 single bonds, 2 double)
      {"S1",   "[S]=*",   0},
      {"S2",   "[SX2,sX2]", 0},
      {"S3",   "[SX3,sX3]", 0},
      {"S",    "[S,s]", 0},

      // Silicon
      {"SI1",  "[Si;X4]", 0}, // tetragonal Si
      {"SI",   "[Si]",    0}  // Si any other
   };

   RDKit::RWMol new_mol(mol); // Create a modifiable copy of the input molecule

   std::vector<std::tuple<std::string, std::string, int> > smarts_by_element = get_smarts_by_element();
   smarts_list.insert(smarts_list.end(), smarts_by_element.begin(), smarts_by_element.end());

   for (const auto& smarts_entry : smarts_list) {
      const std::string& type_energy = std::get<0>(smarts_entry);
      const std::string& smarts_pattern = std::get<1>(smarts_entry);
      RDKit::RWMol* pattern = RDKit::SmartsToMol(smarts_pattern);

      if (pattern) {
         std::vector<RDKit::MatchVectType> matches;
         RDKit::SubstructMatch(new_mol, *pattern, matches);

         for (const auto& match : matches) {
            for (const auto& atom_match : match) {
               RDKit::Atom* atom = new_mol.getAtomWithIdx(atom_match.second);
               atom->setProp("type_energy", type_energy);
            }
         }
         delete pattern;
      }
   }

   // 


}
