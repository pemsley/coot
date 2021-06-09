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

    # Refmac and monomer_library mean the same thing.

    prop_str = 'type_energy'
    if types_type == 'Parm@Frosst':
        prop_str = 'pf_atom_type'
    try:
        # print('match {} match_atom_index {} types_type {} '.format(match, match_atom_index, types_type))
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
                if types_type != "Parm@Frosst":
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

def set_monomer_library_atom_types(mol):
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

def set_amber_atom_types(mol):
    smarts_list = [

        ('C',   '[CX3]=[OX1]', 0), # sp2 C carbonyl group
        ('CA',  'c', 0), # aromatic C
        ('CB',  'c', 0), # aromatic C of fused 5,6 membered rings
        ('CC',  '[cr5;H1]', 0), # C in HIS
        ('CD',  'C=C-C=C', 1), # CD in C=CD-CD=C
        ('CK',  'n[cr5]n', 1), # C in the 5 membered ring of purines
        ('CM',  '[cr6;H1]n', 0), # C in pyrmidines, positions 5,6
        ('CM',  '[cr6;H1]c', 0), # C in pyrmidines, positions 5,6
        ('CN',  '[cr6]n[cr5]', 0), # aromatic C 5,6 membered rings. Is that enough?
        ('CQ',  'n[cr5]n', 1), # C of 5 membered ring in purines N-C-N
        ('CR',  'c[cr5]n', 1), # as above but in HIS
        ('CT',  '[CX4]', 0), # sp3 C aliphatic
        ('CV',  'c', 0), # arom 5 membered ring with 1 N and 1 H (HIS)
        ('CW',  'c', 0), # arom 5 membered ring with 1 N-H and 1 H (HIS)
        ('C*',  'c', 0), # arom 5 membered ring with 1 subst. (TRP)
        ('C*',  'c', 0),  # arom 5 membered ring with 1 subst. (TRP)
        ('CY',  'C#N', 0), # nitrile C
        ('CZ',  'A=C=A', 1), # sp C
        ('CZ',  'C#A', 1), # sp C

        ('O',  'O=C',  0) # carbonyl oxygen
        ('OH',  'COH', 1) # alcohol oxygen
        ('OS',  'C(=O)OA', 2) # ethyl or ester oxygen
        ('OW',  'HOY', 1) # water oxygen
        ('O2',  'O=CO', 1) # carboxyl oxygen
        ('O2',  'O=PO', 1) # phosphate

        ('SI',   '[Si]',    0)  # Si any other
        ]

    full_list = smarts_list
    for smarts_info in full_list:
        atom_type, smarts, match_atom_index = smarts_info
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            if False:
                print("Amber SMARTS ", smarts)
                print("  ", atom_type, ": ", matches)
            for match in matches:
                set_atom_type(match, match_atom_index, mol, atom_type)
        else:
            # print "SMARTS ", smarts, " --- No hits  "
            pass


def set_parmfrosst_atom_types(mol):

    electroneg      = '[N,n,O,o,F,S,s,Cl,Br]'
    ring_electroneg = '[N,n,O,o,s,S]' # are the non-arom needed?
    branched_electroneg = '(~' + electroneg + ')'
    g = '[!H&!He&!Ne&!Ar&!Xe&!Kr&!Rn]'
    g_and_ar = '[c,n,s]' # I believe
    # ar = in aromatic-5-ring or aromatic-6-ring
    # r3 = in a 3 membered ring: g~g~g~@1


    # print('debug H3a: {}'.format('[H][C^3;H1]'+branched_electroneg+branched_electroneg+electroneg))
    smarts_list = [

        # Hydrogen
        # initally from http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html
        #('H2', '[H][C^3;H1][O]', 0), # O is an example electroneg, more specific than H1
        ('HP',  '[H][C^3;X4][+]',    0), # H on a CT next to atom with formal positive charge
        ('H3a', '[H][C^3;H1]'+branched_electroneg+branched_electroneg+'~'+electroneg, 0),
        ('H2a', '[H][C^3;H1]'+branched_electroneg+electroneg, 0),
        ('H2a', '[H][C^3;H2]'+branched_electroneg+electroneg, 0),
        # ('H5',  '[H][c^2;H1]([n])[n]', 0),
        ('H5',  '[H][c^2;H1]'+branched_electroneg+'~'+electroneg, 0),
        ('H5',  '[H][C^2;H1]'+branched_electroneg+'~'+electroneg, 0),
        ('H1a', '[H][C^3;H1]'+electroneg, 0),
        ('H1b', '[H][C^3][N,n]',  0), # {aliph H; 1 elneg grp}
        ('H1e', '[H][C^3]' + electroneg,  0), # {aliph H; 1 elneg grp}
        ('H4a', '[H][C^2,c^2]'+electroneg, 0), # H4 is H-C&sp2~elneg
        ('H4b', '[H][C^2,c^2]~'+electroneg, 0),
        ('H4c', '[H][C^2,c^2;H1][n]', 0), # was '[H][c^2;H1][n][H]'
        ('H4d', '[H][c^2;H1][s]', 0), # not electroneg! parm_Frosst.pcp inconsistent!
        ('HA', '[H][C^2]',       0),
        ('HA', '[H][cr;^2]',     0), # aromatic H, HA also bonds to [C^2]
        ('H',  '[H][NX3;H1;^2]', 0), # both this  and the one below? Checkme
        ('H',  '[H][nX3;H1;^2]', 0),
        ('HC', '[H][CH3;^3]',    0),
        ('HC', '[H][c,CH1]',     0),
        ('HC', '[H][c,CH2]',     0),

        # Amber weirdness:
        # H-O&x2-C&sp3-O&x2-H        > HX * * * HX        ; {acetal hydroxyl must have vdW rad}
        # H-O&x2-C&sp3-N&sp3        > HX * * *        ; {acetal hydroxyl must have vdW rad}
        ('HXa', '[H]-[OX2]-[C^3]-[OX2]-[H]',  (0,4)), # these are by the book
        ('HXb', '[H]-[OX2]-[C^3]-[N^3]',      0),
        ('HXc', '[H]-[OX2]-[C^3](N)(N)',      0), # this is to match H14 in PF-113 -
                                                  # it's not the book acetal pattern
        ('HXe', '[H]-[OX2]-[C^3]([N^3])[c]', 0),  # ;H1 on N added to distinguish
                                                  # PF-345 and PF-261 (they look very close
                                                  # to me.
        ('HXf', '[H]-[OX2]-[C^3][N;X3]', 0), # this is acetal hydroxyl
        ('HXg', '[H]-[OX2]-[C^3][n;X3]', 0), # aromatic N

        ('HO', '[H][O;H1]',      0),
        ('HS', '[H][S]',         0),
        ('HW', '[H][O;H2]',      0),
        ('H2b', '[H][C^3]'+branched_electroneg+electroneg, 0), # not sure this works

        # fallback
        ('H', '[H]',    0),

        # ------------------------------- pathological specialities ------------------------------
        #
        # fused 5,6 unsat systems
        # the v and ws go together (the same SMART) for normal 5,6 ring and "his-like" sets:
        # normal 5,6 ring
        ('CBv', 'c12aaaac1ccn2',  (0,5)), # indole, hit the C[5,6] carbons
        ('CCv', 'c12aaaac1ccn2',  6),     # indole, hits a non-fused carbon
        ('CWv', 'c12aaaac1ccn2',  7),     # indole, hits a non-fused carbon

        # 5,6 { normal 5,6-ring } CB CB CC CW *, add in the branch so that it
        # doesn't hit sat carbon in 6-ring (this may not be the right thing to catch)
        ('CW-row-2', '[c;r5;R2]([C;X3;r6])[c;r5;R2]([C;X3])[c;r5][c;r5]'+electroneg, 3),
        ('CC-row-2', '[c;r5;R2]([C;X3;r6])[c;r5;R2]([C;X3])[c;r5][c;r5]'+electroneg, 2),
        ('CB-row-2', '[c;r5;R2]([C;X3;r6])[c;r5;R2]([C;X3])[c;r5][c;r5]'+electroneg, (0,1)),
        # catch C1 pf4 (hack?)
        ('CW-row-2b', '[c;r5;R2]1[c;r5;R2][c;r5][c;r5]'+electroneg+'1', 3),
        # PF_12
        ('CB-row-2c', '[c;r5]1:[c;r5]~[#6;r5]~[#6;r5][NH]1', (0,1)),
        # ('CC-row-2c', '[c;r5;r6;R2]1:[c;r5;r6;R2][#6;r5][#6;r5]'+electroneg+'1', 2),
        # ('CW-row-2c', '[c;r5;r6;R2]1:[c;r5;r6;R2][#6;r5][#6;r5]'+electroneg+'1', 3),

        # As it stands c12aaaac1aaa2 has ~50% wrong hits
        # ('CBw', 'c12aaaac1aaa2',  (0,5)), # hit the C[5,6] carbons # which way round do we go?
        # ('CCw', 'c12aaaac1caa2',  6),     # hits a non-fused carbon
        # ('CWw', 'c12aaaac1aca2',  7),     # hits a non-fused carbon

        # 5,6 ring with "HIS-like" 5-ring (does this hit anything? second ring is 5 (or 4???))
        ('CBx', '[cr5]12aaaa[c5]2-'+electroneg+'=[cr5]-'+ring_electroneg+'1', (0,5)),
        ('CRx', '[cr5]12aaaa[c5]2-'+electroneg+'=[cr5]-'+ring_electroneg+'1', 7),

        # 5-ring of indole:  { 5-ring unsat }
        ('CR-indole', '[c;r5;H1]1n[c;R2][c;R2]c1',  0),
        # ('NB-indole', 'c1nccc1',  1),

        # PF_5
        ('CCj', '[cr5]1(N)[cr5][ar5][cr5][cr5]1', 0),
        ('CWj', '[cr5]1(N)[cr5][ar5][cr5][cr5]1', 2), # account for branch atom

        # PF_3
        ('CWi', '[cr5]1[cr5][cr5][sr5][sr5]1', 0),
        ('CCi', '[cr5]1[cr5][cr5][sr5][sr5]1', 1),

        # real histidine-like
        # note: [nH] is not the matched by [n], this HIS Nitrogen atoms are different
        # and will match differently with different routes round the ring.
        #
        ('CCy',       '[cr5]1[cr5;H]'+ring_electroneg+'[cr5]'+ring_electroneg+'1', 1), # HIS
        ('CC-testc',  '[cr5]1[cr5;H][nr5][c][nH]1', 1), # HIS
        ('CC-testcc', '[cr5]1[cr5H][nr5;H][c][nr5]1', 1), # HIS, other way
        ('CWy', '[cr5]1[cr5]'+ring_electroneg+'[cr5]'+ring_electroneg+'1', 0), # HIS
        ('CWy', '[cr5]1[cr5]'+ring_electroneg+'[cr5]'+'[nH]1', 0), # HIS
        ('CRy', '[cr5]1[cr5]'+ring_electroneg+'[cr5]'+ring_electroneg+'1', 3), # HIS # success - all hits!
        #
        # maybe we should allow CWy to match CC and CCy to match CW - 2 way rounds an imidazole

        # * * * CR * ; { 5-ring unsat }
        ('CR-5ring-A', g+'1'+g+'N=C'+g+'1', 0), # what does this catch?
        ('CR-5ring-B', '[C;H0;X3]1=NC=CC1', 0), # H0 may not be needed.
        ('NB-5ring-B', '[C;H0;X3]1=NC=CC1', 1), # H0 may not be needed.

        # 5-ring unsat C&ar5=N&ar5-g&ar5~g&ar5~g&ar5-@1        -> CR NB * * * (above the comment)
        #
        # the C=N double bond above is tested for here by X2 on the N (atom index 1)
        ('CR-5ring-a',  '[C,c;^2;R1]1~[n;R1;X2]~[r5;^2;R1;!n]~[ar5;R1]~[C,c;^2;R1]1', 0),
        ('NB-5ring-a',  '[C,c;^2;R1]1~[n;R1;X2]~[r5;^2;R1;!n]~[ar5;R1]~[C,c;^2;R1]1', 1),
        # with a sp3 - bridging C atoms can be aromatic
        ('CR-5ring-b',  '[C,c;^2]1~[N]~[r5;^3]~[r5]~[C,c]1', 0),
        ('NB-5ring-b',  '[C,c;^2]1~[N]~[r5;^3]~[r5]~[C,c]1', 1),
        # for PF4 N4
        ('CR-5ring-c',  '[C,c;^2;R1]1~[n;R1;X2]~[a;r5;^2;R1]~[ar5;R2]~[a;^2;R2]1', 0),
        ('NB-5ring-c',  '[C,c;^2;R1]1~[n;R1;X2]~[a;r5;^2;R1]~[ar5;R2]~[a;^2;R2]1', 1),
        # for PF4 N6
        ('NB-5ring-d',  '[#6;R1]1~[#7;R1]~[#6;R2]~[#7;R2]~[#7;R1]1', 1),
        # for PF4 N5
        ('NB-5ring-d',  '[#6;R2]1~[#7;R1]~[#7;R1]~[#6;R1]~[#7;R2]1', 1),
        # for PF8 N3

        # Note to reader: compare N2 in PF11 and N2 in PF_15.  Extremely similar.  But different types
        #                         NB             NA
        # maybe because of the non-H substituents on atoms idx-3 and idx-4?
        # so add H1 to make it "his-like" - testing.

        # ('NB-5ring-e',  '[c;R1]1~[n;R1]~[a;R1]~[a;R1]~[a;R1]1', 1), # idx-2 atom not a (too general)
        ('NB-5ring-e2', '[c;R1]1~[n;R1]~[n;R1]~[c;R1;H1]~[c;R1;H1]1', 1), # idx-2 atom not a
        ('NB-5ring-e3', '[c;R1]1~[n;R1]~[c;R1]~[n;R1]~[s;R1]1', 1), #
        # for PF8 N4
        ('NB-5ring-f',  '[c;R1]1[n;R1][s;R1][a;R1][a;R1]1', 1), # atom idx-2 needs to be "s" not "a"
        # for PF9 N4
        ('NB-5ring-g',  '[c;R1]1~[n;R1]~[n;R2]~[a;R2]~[c;R1]1', 1),
        # note to self - there may be more CR NBs
        ('CR-5ring-g',  '[c;R1]1~[n;R1]~[n;R2]~[c;R2]~[c;R1]1', 0),
        ('CR-5ring-h',  '[c;R1]1~[n;R1]~[n;R1]~[c;R1]~[c;R1]1', 0), # PF11
        ('NB-5ring-i',  '[c;R1]1~[n;R1;H0]~[n;R2;X3]([c;r6])~[c;R2]~[n;R1;H0]1', 1), # N3 in PF14,  but not N2 in PF_4
        ('NB-5ring-j',  '[c;R1]1~[n;R1]~[n;R1]~[c;R1]~[c;R1]1', 2), # N1 in PF15?

        # 5-ring unsat C&ar5=C&ar5-C&ar5=C&ar5-g&ar5-@1 -> CW CC CC CW
        ('CW-5ring-b',  '[cr5]1[cr5][cr5][cr5][ar5]1', (0,3)),
        ('CC-5ring-b',  '[cr5]1[cr5][cr5][cr5][ar5]1', (1,2)),

        # 5-ring unsat C&ar5=C&ar5-C&ar5=C&ar5-g&ar5-@1        > CW CC CC CW
        ('CWc',  '[cr5;R1]1[cr5][cr5][cr5][nr5]1', (0,3)),
        ('CCc',  '[cr5;R1]1[cr5][cr5][cr5][nr5]1', (1,2)),
        ('CWg',  '[cr5;R1]1[cr5][cr5][cr5][ar5]1', (0,3)),
        ('CCg',  '[cr5;R1]1[cr5][cr5][cr5][ar5]1', (1,2)),

        # 5-ring unsat CW CC * * *, Note CR NB is higher than this
        ('CW-5ring-c',  '[cr5]1[cr5][cr5][cr5][nr5]1', 0),
        ('CC-5ring-c',  '[cr5]1[cr5][cr5][cr5][nr5]1', 1),
        ('CW-5ring-d',  '[c;r5]1[c;r5][c;r5][n;r5][n;r5]1', 0),
        ('CW-5ring-d',  '[c;r5]1[c;r5][c;r5][n;r5][n;r5]1', 2), # symmetry of above
        ('CC-5ring-d',  '[c;r5]1[c;r5][c;r5][n;r5][n;r5]1', 1),

        # 5-ring unsat, NB has no H, NA has H, i.e. HIS C=N-C vs C-N(H)-C
        ('NBz',  '[nr5;R1]1[cr5;R1][cr5R1][nr5;R1;X2][ar5R1]1', 3),
        ('NAz',  '[nr5;R1]1[cr5;R1][cr5R1][nr5;R1;X2][ar5R1]1', 0), # does this find anything?
        ('CRz',  '[nr5;R1]1[cr5;R1][cr5R1][nr5;R1;X2][ar5R1]1', (1,2)),

        # ('CR-5-ring-j', '[c]1[n;H][n][c][c]1', 0), # catches C2 in PF_15 # too late
        # ('CW-5-ring-j', '[c]1[n;H][n][c][c]1', 3), # catches C4 in PF_15

        # n1cnnc1
        # ('NBzz',  '[nr5;R1]1[cr5;R1][nr5R1][nr5R1][ar5R1]1', (2,3)),
        # ('CRzz',  '[nr5;R1]1[cr5;R1][nr5R1][nr5R1][ar5R1]1', 1),
        # we can be more specific and catch the other carbon:
        # we may need to generalize again (just comment above?)
        ('NBzz',  '[nr5;R1]1[cr5;R1][nr5R1][nr5R1][cr5R1]1', (2,3)),
        ('CRzz',  '[nr5;R1]1[cr5;R1][nr5R1][nr5R1][cr5R1]1', (1,4)),

        # non-aromatic (a ring carbon atom has 2 hydrogens)

        # 5-ring unsat, NB has no H, NA has H, i.e. HIS C=N-C vs C-N(H)-C
        ('NBna',  '[N;r5;R1]1=[C;r5;R1][C;r5;R1][N;r5;R1;X2]=[A;r5;R1]1', (0,3)),
        ('CRna',  '[N;r5;R1]1=[C;r5;R1][C;r5;R1][N;r5;R1;X2]=[C;r5;R1]1', 4), # PF20, non-arom though

        # ------------------------------ back to sanity ---------------------------

        # Carbon

        # carbon order: CP {rings/fused-rings} C CR CB C* CA CM C2 CJ CT
        #

        ('CPa', '[cr6]-[cr6]', 0), # biphenyl bridge
        ('CPb', '[cr5]-[cr5]', 0), # biphenyl bridge
        ('CPe', 'c1c(cccccc)aacc1', 1), # biphenyl bridge, hits C7 in PF23

        # ('CCb',  '[R]1:[c]:[R]:[R]:[R]:1', 1), # hits things that should be CR or CW
        ('CCc',   '[cr5][c;H0;X3;r5][cr5][nr5][cr5]', 1),

        # ('CCd', '[cr5;X3;^2]', 0), # 5 ring unsaturated

        ('CRa',   '[cr5]1[nr5;H0][ar5][ar5][ar5]1', 0), # this is useful
        ('CRb',   'n[cr5;R2]n', 1),
        ('CRc',   'n[cr5;R1;H0]n', 1),
        ('CRd',   'n1[cH1]ncc1', 1), # NE2 HIS

        # ('Ca',  '[c]:[cX3]=[N]', 1), # was '[c]~[cX3]=[N]' for C10 in PF7, but hits CRs
        ('Ca',  '[C]-[CX3]=[N]', 1),
        ('Cb',  '[CX3]=[O]',     0),
        ('Cc',  '[cX3]=[O]',     0),
        ('Cd',  '[H][CR0]=[N]',  1),
        ('Ce',  '[H][CR0]=[N]',  1), # also conyl
        ('Cf',  '[!C][CR0]=[N]', 1), # also conyl

        # CB is fused aromatic bridge-head
        ('CBa', 'c12aaaac1aaa2', 0), # doesn't hit anything
        ('CBb', 'c12aaaac1aaaa2', (0,5)),
        ('CBc', 'c12aaaac1aan2',  (0,5)), # indole

        # atom CG of a TYR has type CA
        #
        ('CAa', '[a][cr6;H0;R1]([C,c,N,n])[a]', 1),
        ('CAb', 'c', 0),  # CA on fusion atoms in PF6
        ('CAc', '[cr6;X3]', 0),
        ('CAd', '[cr6;H1]', 0),
        ('CAe', '[cr6]',    0), # Does this work? (PF6)
        ('CAf', '[Cr6]',    0), #  [CR6] means in 6 rings :-)
        ('CAg', 'c1ccccc1',  (0,1,2,3,4,5)),
        ('CAh', 'c1ccccn1',  (0,1,2,3,4)),
        ('CAg', '[C^2]1~[C^2]~[C^2]~[C^2]~[C^2]~[C^2]1',  (0,1,2,3,4,5)),
        ('CAh', '[C^2]1~[C^2]~[C^2]~[C^2]~[C^2]~N1',  (0,1,2,3,4)),

        # SMARTS put an atom in a 5 ring before a 6 ring. but CA is 6-ring
        ('C*', '[cr5;^2;X3]', 0),
        ('C*', '[C;X3;^2]1C(=O)[NH]~cc1', 0), # C8 in PF_12 is not aromatic (by rdkit) but is a PF C*. Grr!
                                              # So provide the specified ring.

        ('CMa', 'n1cnccc1', (3,4,5)), # positions 5 or 6 in pyrimidine # or 4 presumably!
        ('CMb', '[C^2]=[C^2]', (0,1)),  # by validation # catches C7 in PF_0 {non-aromatic C=x}
        ('CMb', '[C;^2;R0]=[C^2]', 0),  # by validation
        ('CMc', '[C;^2;R0]=A', 0),   # by validation
        ('CMd', '[C]=a', 0),   # by validation

        ('C2', '[CX2]#[CX2]', (0,1)),   # by validation
        ('C2', '[CX2]#A', 0),   # by validation
        ('C2', 'A#[CX2]', 1),   # by validation [1]

        ('CWd', '[cr5;X3;^2][c]=O', 0),
        ('CWa', 'n[cr5;R2][cr5;R2]c', 1), # C4 in PF4

        # cyclopropyl and epoxide
        ('CJa',   '[CX4]1[CX4][CX4]1',  (0,1,2)),
        ('CJb',   '[CX4]1[CX4]A1', (0,1)), # remove unindex symmetry

        ('CT', '[CX4]', 0), # bonded to 4 things
        # Carbon fallback
        ('C-fallback',  '[C,c]', 0), # sp hybrizided. Hmmm.

        # [1] needed because above [CX2]#A doesn't find the bond both ways.  I suppose
        #     that's the way SMARTS work.


        # Oxygen
        # first is the carboxylate special:
        ('O2-carboxylate', '[C,c][CX3](=O)[O-]', (2,3)),
        # and presumably this too, but not come across it yet)
        # ('O2-carboxylic-acid', '[C,c][CX3](=O)[OH1]', (2,3)),
        # what's the point of using R0 and R1 and R2!?
        ('OSa',  "[OX2;H0;R0]", 0), # ester
        ('OSb',  "[OX2;H0;R1]", 0), # ester
        ('OSc',  "[OX2;H0;R2]", 0), # ester PF-85 - in two rings!
        ('OSd',  "[oX2;H0;r5]", 0), # aromatic, should I add [n]?
        ('OSe',  "[O;-][CX4]", 0),  # alkoxide ion
        ('OSf',  "[cX3][o][cX3](=O)c", 1), # aromatic in ring with ketone oxygen

        # this is too liberal - it matches O1 in PF 319 (10 other good matches though)
        # ('OSd6', "[oX2;H0;r6]", 0), # aromatic, should I add [n]?
        ('OSd6a', "[oX2;H0;r6]1[cR2][cR2]c(=O)cc1", 0), # a couple
        ('OSd6b', "[oX2;H0;r6]1[c](=O)[n]ccc1", 0), # PF-115
        ('OSd6c', "[oX2;H0;r6]1ccc(=O)cc1", 0), # PF-333
        ('OSd6d', "[oX2;H0;r6]1ccc(=O)cc1", 0), #

        # this matches O1 in 319 - bad!
        # ('OSf',  "[cX3]1[o][cX3;R2][cX3R2]cc1", 1), # aromatic in ring with ketone oxygen
        ('OSg',  "[cX3]1[o][cX3;R2][cX3R2]c(=O)c1", 1), # aromatic in ring with ketone oxygen
        ('OSh',  "[cX3]1[o][cX3][nX3]c(=O)c1", 1), # PF-76 ring is aromatic
        ('OSi',  "[cX3]1[o][cX3](=N)ccc1", 1), # PF7


        # should be 'O'? c.f. 319 vs 333
        #('OSg',  "[cX3]1[o][cX3;R2][cX3R2]cc1", 1), # aromatic in ring with ketone oxygen

        ('OW',  "[OH2]", 0), # water
        ('OH',  "[OH1]", 0), # alcohol

        ('O2',  "[O;-]-*", 0), # oxide ion - deloc
        # O2 should match sulphate (and phosphates, I think) but not
        # O=S-{non-oxygen}
        ('O2',  'OP(~O)~O',   0), # O on a phosphate
        ('O2',  'OS(~O)~O',   0), # guess based on above
        ('O2',  'O=P',   0),      # O on a phosphate
        ('O2',  'O=S=O',   (0,2)), # guess based on above

        # horrible hacking c.f. PF2:  O1=S1 -> O
        #                       PF56  O1-S1 -> Ou
        ('O',   'O=[C]-*',   0), # carbonyl oxygen
        ('O',   'O=[c]',     0), # weird ring iiin PF-3
        ('O',   'O=[S;X4]C',   0),
        ('O',   'O=[S;X3]([N])C',   0), # weird
        # this looks very like one of the OSs to me.
        ('O',   'c12occcc1cccc2',   1), # hit PF-319, seems pretty perverse to do this though.
                                        # why isn't it another "ester"? We have
                                        # to have ultra-specific OS types so that
                                        # we don't hit this - grump.
        # fallback
        ('Ou',  'O',   0),
        ('Ou',  'o',   0),


        # Nitrogen
        #
        # Order: NJ NL N3 NC N2 NB N N2 N NA N3 N2 N3 ND N2 N N* N3 Nu
        #
        # N(-H)-C=O,S                > N * * *        ; {explicitly reset amide N to type N}

        # N=N-C&x4-@1 > NJ NJ *        ; {aziridine}
        ('NJ',    'N=N[CX4]',      (0,1)),

        # N&x1#g > NL *                ; {univalent triple-bonded N}
        ('NL',   '[NX1]#A',   0), # triple bond N

        # N&x4 > N3                ; {+ve tetrahedral N}
        ('N3b',   '[N^3;X3]',    0),
        ('N3c',   '[NX4;+]',     0),

        # g=N&ar6-g                > * NC *        ; {pyridine N}
        ('NC',   '[nr6;X2;H0]',  0), # pyridine

        # N&x3(-H)-Acy2-C         > N2 * C        ; {protonated imine}
        ('N2b',   '[NX3;H](=C)C',   0),   # fixme - don-t forget the C here too.

        # N&x3-C!ar=N&x3                > N2 C N2        ; {amidine}
        ('N2ax',   'A[C;^2](=N)N',   (2,3)),
        ('Cax',    'A[C;^2](=N)N',   1),  # NOTE! C here

        # Nitro > N2 ; {nitro group - guess jan 1999 CIBayly}
        ('N2f', 'C[N^2](~[O;X1])~[OX1]', 1),

        # g&ar-N&ar(-H)-g > * NA * * ; {e.g. HIS C-NH-C}
        # I don't understandy why this above pattern doesn't match N3 in PF_15.  Instead that
        # falls down to N (an amino N). Grump.
        # ('NAe', 'A[nH]a', 1), # too liberal, catches N3 in PF_15
        ('NAf', '[a][n;H;r5][a]', 1), # don't match nH in 6 rings - baah.

        # Nsulfon > N3 ; {sulfonamide-type N is pyramidal. Ordered after aromatic amine}
        # N&x3-S,P(~O&x1)~O&x1
        ('N3a', '[nX3]-[S]([OX1])~[OX1]', 0),
        ('N3b', '[nX3]-[S](=[OX1])~[OX1]', 0),
        # N&x2&fcm1-S,P(~O&x1)~O&x1
        ('N3c', '[nX2;+]-[S,P]([OX1])~[OX1]', 0),
        ('N3d', '[nX2;+]-[S,P](=[OX1])~[OX1]', 0),

        # g&ar-Acy1-N(-H,g)-g&ar > * N2 * * ; {more delocalized aniline-type N oct2004 CIBayly}
        # ('N2d', '[a][n;H]a', 1), # was ('N2d', '[a][N;H]a', 1), # not this - catches N3 in PF_15 (should be N)
        # ('N2e', '[a][n](-[a,A])a', 1), # not this - catches N3 in PF_15 (should be N)

        # g&ar-Acy1-N(-!H)-!H        > * N3 * *      ; {disubst aniline-type N apr2005 CIBayly}
        ('N3d', 'a[N;H0]([!H])[!H]', 1),

        # g&sp2-N!ar-H        > * N  * ; {sp2 amino N}
        ('N', '[A;^2][N][H]', 1),
        ('N', '[c][nH][c](=O)', 1), # N3 in PF_15

        # g&ar-N!ar > * N2 ; {sp2 amino N}
        ('N2g', g_and_ar + '-[N]', 1),

        # N&ar > Nstar                ; {sp2 N}
        ('N*', 'n', 0),

        # N&sp3        > N3                ; {regular sp3 N, not charged}
        ('N3', '[N;0;^3]', 0),

        # fall-back nitrogen
        ('Nu',    '[N,n]',      0),

        ('Sua', 'N[S](=O)(O)=N', 1), # Hmm! PF56

        # many flavours of hypervalent sulfur
        #
        ('SOa', '[S](=O)[N]', 0),  # hypervalent sulfur
        ('SOb', '[C^3][S](=O)(=O)[C^2,c^2]', 1),  # hypervalent sulfur, matches PF-22, 22 is aromic C in position 4, 75 is not (although it looks close)
        ('SOc', '[c][S](=O)(=O)[c]', 1),  # hypervalent sulfur
        ('SOd', '[c]-[S;X4](=O)(C)=[N]', 1),   # hypervalent sulfur
        ('SOe', '[n]-[S;X4](=O)(=O)[C]', 1),   # hypervalent sulfur
        ('SOh', '[c]-[S;X4](=O)(=O)[c]', 1),   # should hit 699, 567
        ('SOi', '[c]-[S;X4](=O)(=O)[O]', 1),   # should hit 521
        ('SOj', '[C;^3][S;X4](=O)(=O)[O]', 1), # hits 416

        ('Sa',  '[S,s][S,s]', (0,1)), # sulfide
        ('Sb',  '[s]', 0), # "sulfide" - hmm. PF5
        ('Sc',  '[c,C]S[c]', 1), # PF15
        ('Sd',  '[c,C]SN', 1), # PF85

        # let's try something less specific than that
        ('Se',  '[S]1[CX3][AX2]=[AX3][C]1', 0), # does this hit anything?
        ('Sf',  '[S]1[CX3]=[A][A][C]1', 0), # does this hit anything?
        ('Sg',  '[C^3][SX2][n]', 1), # should hit 609
        ('Sh',  '[C^3][SX2][c]', 1), # should hit 598
        ('Si',  '[C^3;R1]1[SX2;R1][CR1]AAA[C^2]1', 1), # should hit 496

        # sulfur
        ('Sub', 'S', 0),
        ('Suc', 's', 0), # yikes

        # P
        ('P',    'P', 0),
        # Cl
        ('CL',   '[Cl]', 0),
        # F
        ('F',    '[F]',  0),
        # Br
        ('BR',   '[Br]',  0),

        ]

    for smarts_info in smarts_list:
        atom_type, smarts, match_atom_index = smarts_info
        # print("SMARTS: {}".format(smarts))
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            if False:
                print("SMARTS ", smarts)
                print("  ", atom_type, ": ", matches)
            for match in matches:
                # print('set_atom_type', match, match_atom_index, mol, pattern, atom_type)
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
