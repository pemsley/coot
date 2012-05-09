
import sys
from subprocess import call
from rdkit import Chem
from rdkit.Chem import AllChem

import coot_libs as coot
from jay_util import *

global run_mogul
run_mogul = True

def make_mogul_ins_file(mogul_ins_file_name, mogul_out_file_name, sdf_file_name):
   f = open(mogul_ins_file_name, 'w')
   if f:
     f.write('mogul molecule file ')
     f.write(sdf_file_name)
     f.write('\n')
     f.write('mogul output   file ')
     f.write(mogul_out_file_name)
     f.write('\n')
     f.write('mogul output distribution all on\n')
     f.write('bond all\n')
     f.write('angle all\n')
     f.write('torsion all\n')
     f.write('ring all\n')
     f.close()
   return f


def execute_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name):
   f = make_mogul_ins_file(mogul_ins_file_name, mogul_out_file_name, sdf_file_name)
   if f: 
      # print 'now run mogul using ins file %s' % mogul_ins_file_name
      if run_mogul:
          call(['mogul', '-ins', mogul_ins_file_name])
      return True
   else:
      return False

def atom_name_from_atomic_number_and_count(element, count):
    name = element
    name += str(count)
    return name

def add_atom_names(mol):
    nz = {}
    atom_names = []
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z in nz:
            nz[z]  = nz[z] + 1
        else:
            nz[z] = 1;
        ele = atom.GetSymbol().upper()
        name = atom_name_from_atomic_number_and_count(ele, nz[z])
        p_name = pad_atom_name(name, ele)
        atom.SetProp("name", p_name)
        atom_names.append(p_name)
    return atom_names

def convert_to_coot_bond_type(rdkit_type):
    out_type = 'single'
    if (rdkit_type == Chem.rdchem.BondType.SINGLE):
        out_type = 'single'
    if (rdkit_type == Chem.rdchem.BondType.AROMATIC):
        out_type = 'aromatic'
    if (rdkit_type == Chem.rdchem.BondType.DOUBLE):
        out_type = 'double'
    if (rdkit_type == Chem.rdchem.BondType.TRIPLE):
        out_type = 'triple'
    if (rdkit_type == Chem.rdchem.BondType.ONEANDAHALF):
        out_type = 'deloc'
    return out_type

def pad_atom_name(name, element):

    padded = name
    if (len(element) == 1):
        if (len(name) == 2):
            padded = ' ' + name + ' '
        if (len(name) == 3):
            padded = ' ' + name
    if (len(element) == 2):
        if (len(name) == 2):
            padded = name + '  '
        if (len(name) == 3):
            padded = name + ' '
    return padded

def is_smiles_files(file_name):
    bits = file_name.rsplit(".")
    if len(bits) > 1:
        return bits[1] == 'smi'
    else:
        return False
        

# return the contents of file_name
def read_file(file_name):
    f = open(file_name)
    return f.read()

        
def make_restraints_for_bond_orders(mol):
    restraints = {}
    bond_list = []
    for bond in mol.GetBonds(): 
        type = bond.GetBondType()
        coot_bond_type = convert_to_coot_bond_type(type)
        at_1 = bond.GetBeginAtom()
        at_2 = bond.GetEndAtom()
        name_1 = at_1.GetProp('name')
        name_2 = at_2.GetProp('name')
        item = [name_1, name_2, coot_bond_type, 1.0, 1.0]
        bond_list.append(item)
    restraints['_chem_comp_bond'] = bond_list
    restraints['_chem_comp'] = [mol.GetProp('comp_id'),
                                mol.GetProp('comp_id'),
                                mol.GetProp('name'),
                                'non-polymer',
                                mol.GetNumAtoms(),
                                mol.GetNumAtoms(),
                                'Partial']
    return restraints

# match_atom_index can be of type int or a list - otherwise trouble.
#
def set_atom_type(match, match_atom_index, mol, atom_type):
    # print "   trying to set ", match_atom_index, " of ", match
    try:
        this_atom = match[match_atom_index]
        # print "   this_atom:", this_atom
        try:
            current_type = mol.GetAtomWithIdx(this_atom).GetProp("atom_type")
            # print "   oops - atom ", mol.GetAtomWithIdx(this_atom).GetProp("name"), " already has type ", current_type
        except KeyError:
            mol.GetAtomWithIdx(this_atom).SetProp("atom_type", atom_type)
            print "    set atom number ", this_atom, " having name ", mol.GetAtomWithIdx(this_atom).GetProp("name"), " as ", atom_type
    except TypeError:
        for match_atom in match_atom_index:
            set_atom_type(match, match_atom, mol, atom_type)

def ele_to_smart(v):
    return (v.upper(), '['+v+']', 0)

# those not handled by hand-coding
def smarts_by_element():
   eles = [
      "He", "Li", "Be", "B",  "Ne", "Na", "Mg", "Al", "Si",
      "Ar", "K", "Ca",  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
      "Zn", "Ga", "Ge", "As", "Se",      "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
      "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
      "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
      "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
      "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"]
   return map(ele_to_smart, eles)
   

def set_atom_types(mol):
    smarts_list = [

        # Oxygen
        ('O2',  "[OX2;H0]", 0), # ester, Os between P and C are O2, not OP
        ('OP',  'O~P',   0),
        ('OS',  'O~S',   0),
        ('OB',  'O~B',   0),
        ('OC',  '*C(=O)[OH]', (2,3)), # carboxylate
        ('OH1', '[OH1]', 0), # alcohol
        ('O2',  "[oX2;H0]", 0), # ring oxygen
        ('O',   'O=*',   0), # carbonyl oxygen
        # Fallback oxygen
        ('O',   'O',   0), 
        
        # Carbon SP2
        ('CR56', 'c12aaaac1aaa2',  (0,5)), # works on indole
        ('CR66', 'c12aaaac1aaaa2', (0,5)),
        ('CR6',  'c12ccccc1OCO2',  (0,5)),   # mouse, fused atoms in 6-ring not non-arom 5-ring
        ('CR66', 'c12aaaac1AAAA2', (0,5)),   # one 6-ring aromatic, other not.
        ('CR6',  'c12caccc1***2',  (0,5)),  # aromatic 6, (non-)aromatic 5, maybe this should be CR56?

        ('CR16', '[cr6;H1]',  0),
        ('CR6',  '[cr6;H0]',  0),
        ('CR15', '[cr5;H1]',  0),
#        ('CR5',  'C1(=O)[C,c][C,c]C=N1', 0), # carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
        ('CR5',  '[cr5;H0]',  0),
        ('CR5',  '[CR5;H0]',  0),
        ('C1',   '[CX3;H1]',    0),  # double bond, single bond and one H
        ('C',    '[CX3;H0;^2]', 0),
        ('C',    '[CX3]=[OX1]', 0),  # carbonyl carbon
        ('C',    '[$([CX2](=C)=C)]',   0), # bonded to 3 things not hydrogen

        # Carbon SP3
        ('CT',   '[CX4H0]', 0), # single bonded to 4 things not hydrogen
        ('CH3',  '*[CH3^3]',    0), # bonded to 3 Hs
        ('CH2',  '[C;^3;H2]',   0), # standard aliphatic C.
        ('CH1',  '*[C](*)*',    1), # bonded to H and 3 things --- ??? mistake?

        # sp??? needs sorting 
        ('CH2',  '[CH2]',   0), # bonded to 2 hydrogens
        
        # Carbon fallback
        ('C', '[C,c]', 0),

        # Hydrogen
        ('HCH1', '[H][CH1]',    0),
        ('HCH2', '[H][C;H2^3]', 0),
        ('HCH3', '[H][CH3]',    0),
        ('HNC1', '[H][NX2;H1;^2]', 0), # H of N of N=C ? 
        ('HNT2', '[H][NX3;H2;^3]', 0), # H connected to type NT2
        ('HNH2', '[H][NH2;^2]', 0), # NH2 is sp2
        ('HNH1', '[H][NH1]',    0),
        ('HCR6', '[H][cr6;H1]', 0),
        ('HCR1', '[H]c',        0),
        ('HNH1', '[H][NH1]',    0),
        ('HOH1', '[H][OH1]',    0),
        ('HCH',  '[H]',         0),
        
        # Nitrogen, SP3
        ('NH2', '[NX3][CX3]=[NX3+]', (0,2)), # amidinium (charged)
        ('NT1', '[NX4;H1;^3]',  0),
        ('NT1', '[NX3;H1;^3]',  0),
        ('NT2', '[NX3;H2;^3]',  0), # different to mon-lib!
        ('NT',  '[NX3;H0;^3]',  0),
        
        # Nitrogen, SP2
        ('NR15', '[nr5;X3;H1]',    0),
        ('NR5',  '[nr5;X3;H0]',    0),
        ('NRD5', '[nr5;X2;H0]',    0), # guess from 071
        ('NRD5', 'C1(=O)[C,c][C,c]C=N1', 5), # N bonded to carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
        ('NR16', '[nr6;H1]',    0),
        ('NRD6', '[nr6;H0]',  0), # guess from 084
        ('NR6',  '[nr6]',    0),
        ('NC1',  '[NX2;H1;^2]', 0),  # N of N=C ? 
        ('NH1',  '[NX3;H1;^2]', 0),
        ('NH2',  '[NX3;H2;^2]', 0),  # sp2, e.g. ND2 of an ASP

        # fall-back nitrogen
        ('N',    '[N]',      0),  # how to say n or N?

        # Phosphorus
        ('P',    'P', 0),
        # Cl
        ('CL',   '[Cl]', 0),
        # F
        ('F',    '[F]',  0),
        # Br
        ('BR',    '[Br]',  0),

        # Sulfur
        ('SH1',  '[SH1]', 0),  # SG of CYS
#       ('S1',   'S=a',   0), # don't know how to specify exactly one double bond
#         ('ST',   'C[S](=O)N', 1), # guess based on 059
#         ('ST',   'c[S](=O)N', 1), # guess based on 059
        ('ST',   '[SX4]', 0), # tetrahedral (2 single bonds, 2 double)
        ('S2',   '[SX2,sX2]', 0),
        ('S3',   '[SX3,sX3]', 0),
        ('S',    '[S,s]', 0)

        ]

    full_list = smarts_list

    for item in smarts_by_element():
       full_list.append(item)

    for smarts_info in full_list:
        atom_type, smarts, match_atom_index = smarts_info
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            print "SMARTS ", smarts
            print "  ", atom_type, ": ", matches
            for match in matches:
                set_atom_type(match, match_atom_index, mol, atom_type)
        else:
            # print "SMARTS ", smarts, " --- No hits  "
            pass



def make_restraints(smiles_string, comp_id, sdf_file_name, pdb_out_file_name, mmcif_dict_name):

   mogol_exe = which('mogul')

   if (run_mogul):
      if (mogol_exe == None):
         print "mogul not found in path"
         return False

   compound_name = 'some-compound'
   m = Chem.MolFromSmiles(smiles_string)
   m_H = AllChem.AddHs(m)
   AllChem.EmbedMolecule(m_H)
   AllChem.UFFOptimizeMolecule(m_H)
   atom_names = add_atom_names(m_H)
   m_H.SetProp('comp_id', comp_id)
   m_H.SetProp('name', compound_name)
   set_atom_types(m_H)
   Chem.AllChem.ComputeGasteigerCharges(m_H)

   mb = Chem.MolToMolBlock(m_H)
   print >> file(sdf_file_name,'w'), mb
   coot.write_pdb_from_mol(m_H, comp_id, pdb_out_file_name)
   mogul_ins_file_name = 'mogul.ins'
   mogul_out_file_name = 'mogul.out'
   bor = make_restraints_for_bond_orders(m_H)
   # print "we got this bor: ", bor
   # coot.write_restraints(bor, comp_id, 'bond-orders.cif')

   # print out the set types:
   for atom in m_H.GetAtoms():
      charge = atom.GetProp('_GasteigerCharge')
      name   = atom.GetProp('name')
      try:
         atom_type   = atom.GetProp('atom_type')
         is_aromatic = atom.GetIsAromatic()
         hybrid      = atom.GetHybridization()
         print "   atom: ", name, atom.GetSymbol(), " ", is_aromatic, " ", hybrid, "   ", atom_type, "    ", charge
      except KeyError:
         print "miss", name, atom.GetSymbol(), charge

   mogul_state = execute_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name)

   if mogul_state:
      coot.mogul_out_to_mmcif_dict_by_mol(mogul_out_file_name, comp_id,
                                          compound_name, m_H, bor, mmcif_dict_name)
      return True
   return mogul_state



if __name__ == "__main__":

    smiles_string = "CC"
    sdf_file_name = 'test.sdf'
    cif_restraints_file_name = "restraints.cif"
    pdb_out_file_name = "out.pdb"
    comp_id = "XXX"
    smiles_string = 'Cc1ccccc1'
    smiles_string = 'Oc1ccccc1'
    smiles_string = 'O=C(O)c1ccc(O)cc1'
    
    if len(sys.argv) > 2:
        smiles_string = sys.argv[1]
        # if smiles_string ends in .smi, read it as if it was a file
        if is_smiles_files(smiles_string):
            smiles_string = read_file(smiles_string)
        file_stub = sys.argv[2]

        # mogul handling
        if len(sys.argv) == 4:
            if sys.argv[3] == '--no-mogul':
                run_mogul = False

        sdf_file_name = file_stub + ".sdf"
        cif_restraints_file_name = file_stub + ".cif"
        pdb_out_file_name = file_stub + ".pdb"
        comp_id = file_stub
    
    make_restraints(smiles_string, comp_id, sdf_file_name, pdb_out_file_name, cif_restraints_file_name)

