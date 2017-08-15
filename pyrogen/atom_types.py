# Copyright 2012 by the University of Oxford
# Copyright 2015 by Medical Research Council
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

from __future__ import print_function
from rdkit import Chem

# match_atom_index can be of type int or a list - otherwise trouble.
#
# Note that atom_type properties can also have been set in hydrogen_transformations():
#
def set_atom_type(match, match_atom_index, mol, atom_type, types_type='Refmac'):
    prop_str = 'type_energy'
    if types_type == 'Parm@Frosst':
        prop_str = 'pf_atom_type'
    try:
        this_atom_idx = match[match_atom_index]
        try:
            current_type = mol.GetAtomWithIdx(this_atom_idx).GetProp(prop_str)
        except KeyError:
            try:
                mol.GetAtomWithIdx(this_atom_idx).SetProp(prop_str, atom_type)
                name = mol.GetAtomWithIdx(this_atom_idx).GetProp("name")
                if False:
                    print('   actually have set atom number %s having name %s with type %s ' %
                          (str(this_atom_idx).rjust(2), repr(name), repr(atom_type)))
            except KeyError as e:
                print('No name for atom idx', match_atom_index)

    except TypeError:
        for match_atom in match_atom_index:
            # print('   set_atom_type ', match, match_atom, mol, atom_type, types_type)
            set_atom_type(match, match_atom, mol, atom_type, types_type)

def ele_to_smart(v):
    return (v.upper(), '['+v+']', 0)

# those not handled by hand-coding
def smarts_by_element():
   eles = [
      "He", "Li", "Be", "B",  "Ne", "Na", "Mg", "Al", 
      "Ar", "K", "Ca",  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
      "Zn", "Ga", "Ge", "As", "Se",      "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
      "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
      "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
      "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
      "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"]
   return map(ele_to_smart, eles)
   
def set_atom_types(mol):
    smarts_list = [

        # Full coverage for C, H, O.

        # Oxygen
        ('O2',  "[OX2;H0]", 0), # ester, Os between P and C are O2, not OP
        ('OP',  'O~P',   0),
        ('OS',  'O~S',   0),
        ('OB',  'O~B',   0),
        ('OC',  '*C(=O)[OH]', (2,3)), # carboxylic acid
        ('OC',  '*C(=O)O',    (2,3)), # carboxylate, doesn't match deloc bonds
        ('OH1', '[OH1]', 0), # alcohol
        ('O2',  "[oX2;H0]", 0), # ring oxygen
        ('O',   'O=*',   0), # carbonyl oxygen

        # OH2 no examples
        # OHA no examples
        # OHB no examples
        # OHC no examples
        # OC2 no exmampes
        
        # Fallback oxygen
        ('O',   'O',   0),

        # Carbon SP
        ("CSP1", '[H][C]#*',  1), # e.g. in 2GT
        ("CSP",  '[C]#[C]',   (0,1)),
        ("CSP",  '[C]#*',     0),

        # Carbon SP2
        ('CR56', 'c12aaaac1aaa2',  (0,5)), # works on indole
        ('CR56', 'c12aaaan1aaa2',  0), # same pattern as (below) N56, but catching first 56 atom
        ('CR56', 'c12AAAAc1aaa2',  0), # 6-ring doesn't have to aromatic (8UG)
        ('CR66', 'c12aaaac1aaaa2', (0,5)),
        ('CR6',  'c12ccccc1OCO2',  (0,5)),   # mouse, fused atoms in 6-ring not non-arom 5-ring
        ('CR66', 'c12aaaac1AAAA2', (0,5)),   # one 6-ring aromatic, other not. Needed for XXX?
                                             # but makes a fail on 113.
        ('CR6',  'c12caccc1***2',  (0,5)),  # aromatic 6, (non-)aromatic 5, maybe this should be CR56?

        # note CR1  missing - can't find example
        #      CR1H missing - can't find example

        ('CR16', '[cr6;H1]',  0),
        ('CR6',  '[cr6;H0]',  0),
        ('CR15', '[cr5;H1]',  0),
#        ('CR5',  'C1(=O)[C,c][C,c]C=N1', 0), # carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
        ('CR5',  '[cr5;H0]',  0),
        ('CR5',  '[CR5;H0]',  0),
        ('C1',   '[CX3;H1]',    0),  # double bond, single bond and one H
        ('C2',   '[CX3;H2]=*',  0),  # double bond, and 2 H
        ('C',    '[CX3;H0;^2]', 0),
        ('C',    '[CX3]=[OX1]', 0),  # carbonyl carbon
        ('C',    '[$([CX2](=C)=C)]',   0), # bonded to 3 things not hydrogen

        # Carbon SP3
        ('CT',   '[CX4H0]', 0),     # single bonded to 4 things not hydrogen
        ('CH3',  '[C;H3;^3]',   0), # bonded to something via single bond and 3 Hs
        ('CH2',  '[C;^3;H2]',   0), # standard aliphatic C.
        ('CH1',  '*[C;H1](*)*', 1), # bonded to H and 3 things 

        # sp??? needs sorting 
        ('CH2',  '[CH2]',   0), # bonded to 2 hydrogens
        
        # Carbon fallback
        ('C', '[C,c]', 0),

        # Hydrogen
        ('HCH1', '[H][CH1]',    0),
        ('HCH2', '[H][C;H2^3]', 0),
        ('HCH3', '[H][CH3]',    0),
        ('HNC1', '[H][N;H1;^2]~C(~N)~N', 0), # H of N of N=C ? 
        ('HNC2', '[H][NX3;H2;^2]', 0), # H on a NC2 (NH1 and NH2 of ARG)
        ('HNC3', '[H][NX3;H3;^2]', 0), # guess - no examples
        ('HNT1', '[H][NX4;H1;^3]', 0),
        ('HNT1', '[H][NX3;H1;^3]', 0),
        ('HNT2', '[H][NX3;H2;^3]', 0), # H connected to type NT2
        ('HNT3', '[N^3;H3][H]', 1), # NH3+ 
        ('HNH2', '[H][NH2;^2]', 0), # NH2 is sp2
        ('HNH1', '[H][NX3;H1;^2]',    0),
        ('HCR6', '[H][cr6;H1]', 0),
        ('HCR5', '[H][cr5;H1]', 0), # connected to aromatic ring C with 1 H
        ('HNR5', '[H][nr5;H1]', 0), # connected to aromatic ring C with 1 H
        ('HNR5', '[H][Nr5;H1]', 0), # guess based on above, connected to aromatic N in a 5-ring
        ('HNR6', '[H][nr6;H1]', 0), # connected to aromatic 6-ring N with 1 H
        ('HNR6', '[H][NR6;H1]', 0), # guess based on above

        # HCR missing - no examples (and how is it different to HCR1?)
        ('HCR1', '[H]c',        0),
        ('HNH1', '[H][NH1]',    0),
        ('HOH1', '[H][OH1]',    0),
        ('HOH2', '[H][OH2][H]', (0,2)), # H of HOH - water
        ('H',    '[H]',         0),
        
        # Nitrogen, SP3

        ('NT1', '[NX4;H1;^3]',  0),
        ('NT1', '[NX3;H1;^3]',  0),
        ('NT2', '[NX3;H2;^3]',  0), # different to mon-lib!
        ('NT3', '[NX4;H3;^3]',  0),
        ('NT',  '[NX3;H0;^3]',  0),


        # NE-CZ in ARG should be deloc (guandino) - also NC1-C
        # single (as is in ARG.cif) is not found in ener_lib!
       
        # Nitrogen, SP2
        ('NR66', 'c12aaaan1aaaa2', 5), # (second) 66 atom is an N.
        ('NR56', 'c12aaaan1aaa2',  5), # (second) 56 atom is an N.
        ('NR55', 'c12aaan1aaa2',   4), # (second) 55 atom is an N.
        ('NC2',  '[NX3;H2^2]', 0),     # N of sp2 NH2 (as in ARG).
        ('NH2',  '[NX3^2][CX3^2]=[N^2;X3+]', (0,2)), # amidinium (charged)... 
        ('NR15', '[nr5;X3;H1]',    0),
        ('NR5',  '[nr5;X3;H0]',    0),
        ('NR5',  '[NR;X3;H0;^2]',  0), # [NR5;X3;H0;^2] fails on 14C (also at daylight)
        ('NRD5', '[nr5;X2;H0]',    0), # guess from 071
        ('NRD5', 'C1(=O)[C,c][C,c]C=N1', 5), # N bonded to carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
        ('NR16', '[nr6;H1]',    0),
        ('NRD6', 'a:[nr6;X2;H0]:a',  1), # aromatic N with no H, i.e. one double one single
        ('NR6',  '[nr6]',    0),
        ('NC1',  '[H][N;H1;^2]~C(~N)~N', 1), 
        ('NC1',  '[NX3;H1;^2]C(~N)~N', 0), # N, as in NE in ARG
        ('NC1',  '[NX2;H1;^2]', 0),  # N of N=C ? 
        ('NH1',  '[NX3;H1;^2]', 0),
        ('NH2',  '[NX3;H2;^2]', 0),  # sp2, e.g. ND2 of an ASP
        ('NT',   '*n1~[o]~[o]1', 1), # guess from 16X dioxaziridine (bleugh)
        # (NT needs checking?)
        # NC2 no examples
        # NC3 no examples
        # NPA no examples
        # NPB no examples
        

        # Nitrogen SP1
        ('NS',   '[N^1]', 0),
        # NS1 no examples
        

        # fall-back nitrogen
        ('N',    '[N,n]',      0),  

        # Phosphorus
        ('P',    'P', 0),
        # Cl
        ('CL',   '[Cl]', 0),
        # F
        ('F',    '[F]',  0),
        # Br
        ('BR',    '[Br]',  0),

        # Sulfur
        ('SH1',  '[SX2H1]', 0),  # SG of CYS
        ('ST',   '[SX4]', 0),  # tetrahedral (2 single bonds, 2 double)
        ('S1',   '[S]=*', 0),
        ('S2',   '[SX2,sX2]', 0),
        ('S3',   '[SX3,sX3]', 0),
        ('S',    '[S,s]', 0),

        # Silicon
        ('SI1',  '[Si;X4]', 0), # tetragonal Si
        ('SI',   '[Si]',    0)  # Si any other

        ]

    full_list = smarts_list

    for item in smarts_by_element():
       full_list.append(item)

    # print('----------------------------refmac-----------------------------------------')
    for smarts_info in full_list:
        atom_type, smarts, match_atom_index = smarts_info
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
	    if False:
		print("SMARTS ", smarts)
		print("  ", atom_type, ": ", matches)
            for match in matches:
                set_atom_type(match, match_atom_index, mol, atom_type)
        else:
            # print "SMARTS ", smarts, " --- No hits  "
            pass

    # do we return success (everything has a type) or not?
    #
    for atom in mol.GetAtoms():
       try:
          atom_type = atom.GetProp('type_energy')
       except KeyError:
          is_aromatic = atom.GetIsAromatic()
          hybrid      = atom.GetHybridization()
          print("Error:: Missing type for atom \"", atom.GetProp('name'), "\" is_aromatic: ", is_aromatic, " hybridization: ", hybrid)
          return False
    # we got to the end, good
    return True


def set_parmfrosst_atom_types(mol):

    electroneg = '[N,O,F,S,Cl,Br]' # does this also need aromtics?

    smarts_list = [

        # Hydrogen
        # initally from http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html
        ('H2', '[H][C^3;H1][O]', 0), # O is an example electroneg, more specific than H1
        ('H1', '[H][C^3][N]',    0),
        ('H1', '[H][C^3][O]',    0),
        ('H1', '[H][C^3][S]',    0),
        ('HA', '[H][C^2]',       0),
        ('H5', '[H][c^2;H1]([n])[n]', 0),
        ('H4', '[H][c^2;H1][n]', 0), # was '[H][c^2;H1][n][H]'
        ('H4', '[H][c^2;H1][s]', 0), # H4 is H-C&sp2~elneg
        ('H4', '[H][C^2;H1]O',  0),
        ('HA', '[H][cr;^2]',     0), # aromatic H, HA also bonds to [C^2]
        ('H',  '[H][NX3;H1;^2]', 0), # both this  and the one below? Checkme
        ('H',  '[H][nX3;H1;^2]', 0),
        ('HC', '[H][CH3;^3]',    0),
        ('HC', '[H][c,CH1]',     0),
        ('HC', '[H][c,CH2]',     0),
        ('HO', '[H][O;H1]',      0),
        ('HS', '[H][S]',         0),
        ('HW', '[H][O;H2]',      0),
        ('H2', '[H][NX3;H2;^2]', 0),
        ('H3', '[N^3;H3][H]',    1),
        # fallback
        ('H', '[H]',    0),

        # Carbon
        ('CT', '[CX4]', 0), # bonded to 4 things
        ('C',  '[CX3]=[O]', 0),
        ('C',  '[cX3]=[O]', 0),
        ('C',  '[c]~[cX3]=[N]', 1), # C10 in PF7
        ('CB', 'c12aaaac1aaa2', 0),
        ('CB', 'c12aaaac1aaaa2', (0,5)),
        ('CWd', '[cr5;X3;^2][c]=O', 0),
        ('CWa', 'n[cr5;R2][cr5;R2]c', 1), # C4 in PF4

        # by the book, pathological specialities
        ('CR',   '[cr5]1[nr5;H0][ar5][ar5][ar5]1', 0), # this is useful
        ('CWb',  '[cr5;R1]1[cr5][cr5][ar5][nr5;H1]1', 0),
        ('CCa',  '[cr5]1[cr5][cr5][ar5][ar5]1', 1),
        ('CC@',  '[cr5]1[cr5][ar5][ar5][ar5]1', 0),
        ('CWd',  '[cr5;H0]1@[cr5;H1]@[sr5]@[ar5]@[ar5]1', 1), # heuristic C5 in PF5
        ('CWe',  '[cr5;H0]1@[cr5;H0]([I,F,Cl,Br])@[sr5]@[ar5]@[ar5]1', 1), # heuristic C6 in PF5
        ('CC!',  '[cr5;H0]1[sr5][ar5][ar5][cr5]1', 0),

        ('CCb',  '[R]1:[c]:[R]:[R]:[R]:1', 1), # hits things that should be CR or CW
        ('CR',   'n[cr5;R2]n', 1),
        ('CRb',  'n[cr5;R1;H0]n', 1),
        ('CCc',   '[cr5][c;H0;X3;r5][cr5][nr5][cr5]', 1),

        # atom CG of a TYR has type CA
        #
        # comes before CB - bleugh
        # and CC
        ('CAa', '[a][cr6;H0;R1]([C,c,N,n])[a]', 1),
        ('CAb', 'c', 0),  # CA on fusion atoms in PF6

        ('CCd', '[cr5;X3;^2]', 0), # 5 ring unsaturated

        ('CB', '[a][cr6;H0][a]', 1), # before CA

        ('CAc', '[cr6;X3]', 0),
        ('CAd', '[cr6;H1]', 0),
        ('CAe', '[cr6]',    0), # Does this work? (PF6)
        ('CAf', 'c',        0), # this is above CC (PF6)

        ('CK', '[cr5;H1;N2]', 0),
        ('CM', 'n1cnccc1', (3,4,5)), # positions 5 or 6 in pyrimidine # or 4 presumably!
        ('CM', '[C]=[C]', (0,1)),   # by validation
        ('CM', '[C]=A', 0),   # by validation
        ('C2', '[CX2]#[CX2]', (0,1)),   # by validation
        ('C2', '[CX2]#A', 0),   # by validation
        ('C2', 'A#[CX2]', 1),   # by validation [1]
        ('CJ', 'n1n[cX3]cc1', (3,4,5)), # positions 5 or 6 in pyrimidine # or 4 presumably!
#        ('CN', '[C]'), # how is this different to CB?
        ('CQ', 'n1[cH1]nccc1', 1),
        ('CR', 'n1[cH1]ncc1', 0), # NE2 HIS
        ('CV', '[cr5;^2;H1]~n', 0),
        ('C*', '[cr5;^2;X3]', 0),
        # Carbon fallback
        ('C',  '[C,c]', 0), # sp hybrizided. Hmmm.

        # [1] needed because above [CX2]#A doesn't find the bond both ways.  I suppose
        #     that's the way SMARTS work.


        # Oxygen
        ('OS',  "[OX2;H0]", 0), # ester
        ('OS',  "[oX2;H0]", 0), # aromatic, should I add [n]?
        ('OH',  "[OH1]", 0), # alcohol
        ('OW',  "[OH2]", 0), # water
        # O2 should match sulphate (and phosphates, I think) but not
        # O=S-{non-oxygen}
        ('O2',  'OP(~O)~O',   0), # O on a phosphate
        ('O2',  'OS(~O)~O',   0), # guess based on above
        ('O2',  'O=P',   0), # O on a phosphate
        ('O2',  'O=S=O',   (0,2)), # guess based on above
        ('O',   'O=*',   0), # carbonyl oxygen
        # fallback
        ('O',  'O',   0),
        ('O',  'o',   0),


        # Nitrogen
        #
        # T
#       ('N2',   '[R][NX3;H1;^2]',   1), # {sp2 amino N} - sigh
        ('N2',   '[NX3;!R]=[R]',     0), # {sp2 amino N} - sigh PF7 - N1 is not N!
        ('N',    '[NX3;H1;^2]C=O',   0), # amide
        ('N',    '[NX3;H0;^2]C=O',   0), # PF5
        ('NA',   '[nr5;H1]',      0), # both this and the one below? Checkme
        ('NC',   '[nr6;X2;H0]',   0), # pyridine
        ('NB',  '[c][nr;X2;H0]',  0), # {e.g. HIS C=N-C} # does this catch anything? checkme
        ('NB',   '[c,s][nr5;R1;H0]c',  1),
        ('NB',   '[c;H1,s][nr5;R1;H0]n',  1), # hits N4, N5 in PF4.
        ('NB',   '[cr5;R2][nr5;R1;H0]n',  1), # hits N4, N5 in PF4.
        ('N*',   '[nr5;X2;H0]',   0),
        ('N*',   '[nr6;X2;H0]',   0),
        ('N*',   '[nr5;R2;X3;H0]', 0), # N3 in PF4
        ('NT',   '[N^3;X4]',      0),
        ('N3',   '[N^3;X3]',      0),
        ('N2',   '[NX3;H2^2]', 0),     # N of sp2 NH2 (as in ARG).
        ('N2',   'aN', 1),     # book
        ('N2',   '[NX3][C][NX3]', (0,2)),     # book
        ('N2',   '[H][NX3]CO', 1),     # book
        ('N2',   '[H][NX3]C=O', 1),     # book
        ('Nc',    '[NX3;H1;^2]',   0), # amide
        ('NC',   '[NR6;X2;H0]',   0),
        ('NL',   '[NX1]#A',   0), # triple bond N
        # fall-back nitrogen
        ('N',    '[N,n]',      0),

        ('SO', '[S](=O)[N]', 0),  # hypervalent sulfur
        ('S',  '[S,s][S,s]', (0,1)), # sulfide
        ('S',  '[s]', 0), # "sulfide" - hmm. PF5
        # sulfur
        ('Su', 'S', 0),
        ('Su', 's', 0), # yikes
        
        # P
        ('P',    'P', 0),
        # Cl
        ('CL',   '[Cl]', 0),
        # F
        ('F',    '[F]',  0),
        # Br
        ('BR',   '[Br]',  0),

        ]

    # print('----------------------------amber-----------------------------------------')
    for smarts_info in smarts_list:
        atom_type, smarts, match_atom_index = smarts_info
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
	    if False:
		print("SMARTS ", smarts)
		print("  ", atom_type, ": ", matches)
            for match in matches:
                # print('set_atom_type', match, match_atom_index, mol, atom_type)
                set_atom_type(match, match_atom_index, mol, atom_type, types_type='Parm@Frosst')
        else:
            # print "SMARTS ", smarts, " --- No hits  "
            pass

    # do we return success (everything has a type) or not?
    #
    for atom in mol.GetAtoms():
       try:
          atom_type = atom.GetProp('pf_atom_type')
       except KeyError:
          is_aromatic = atom.GetIsAromatic()
          hybrid      = atom.GetHybridization()
          print("Error:: Missing type for atom ", atom.GetProp('name'), " is_aromatic: ", is_aromatic, " hybridization: ", hybrid)
          return False
    # we got to the end, good
    return True
