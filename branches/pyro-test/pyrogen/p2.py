
from subprocess import call
from rdkit import Chem
from rdkit.Chem import AllChem
import coot_libs as coot
from jay_util import *

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


def run_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name):
   f = make_mogul_ins_file(mogul_ins_file_name, mogul_out_file_name, sdf_file_name)
   if f: 
      # print 'now run mogul using ins file %s' % mogul_ins_file_name
      # call(['mogul', '-ins', mogul_ins_file_name])
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
        ele = atom.GetSymbol()
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
        out_type = 'arom'
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
    print "   trying to set ", match_atom_index, " of ", match
    try:
        this_atom = match[match_atom_index]
        print "   this_atom:", this_atom
        try:
            current_type = mol.GetAtomWithIdx(this_atom).GetProp("atom_type")
            print "   oops - atom ", mol.GetAtomWithIdx(this_atom).GetProp("name"), " already has type ", current_type
        except KeyError:
            mol.GetAtomWithIdx(this_atom).SetProp("atom_type", atom_type)
            print "    set ", mol.GetAtomWithIdx(this_atom).GetProp("name")
    except TypeError:
        for match_atom in match_atom_index:
            set_atom_type(match, match_atom, mol, atom_type)

def set_atom_types(mol):
    smarts_list = (

        # Oxygen
        ('O=*',    'O',   0), # carbonyl oxygen
        ('O=*O',   'OC',  0), # OD1 of ASP
        ('O~P~O',  'OP',  0),
        ('O~S~O',  'OS',  0),
        ('O~B',    'OB',  0),
        ("[O;D2]", 'O2',  0), # ester
        ('O-*',    'OH1', 0), # alcohol
        
        # Carbon
        ('[cr5;r6]', 'CR56', 0),
        ('c12ccccc1cccc2', 'CR66', (0,5)),
        ('[cr6]', 'CR6', 0),
        ('[cr5]', 'CR5', 0),
        ('[CH2]', 'C2',  0), # bonded to 2 hydrogens
        ('[$([CX2](=C)=C)]',   'C',   0), # bonded to 3 things not hydrogen
        ('*[C](*)*',   'CH1',   1), # bonded to H and 3 things
        ('[CX3]=[OX1]', 'C', 0),  # carbonyl carbon
        
        # Hydrogen
        ('[H]',   'H',   0),
        
        # Nitrogen
        ('[NX3;H2]',   'NH2',   0),
        ('[NX3;H1]',   'NH1',   0)
        )
    
    for smarts_info in smarts_list:
        smarts, atom_type, match_atom_index = smarts_info
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            print "SMARTS ", smarts
            print "  ", atom_type, matches
            for match in matches:
                set_atom_type(match, match_atom_index, mol, atom_type)




def make_restraints(smiles_string, sdf_file_name, mmcif_dict_name):

   mogol_exe = which('mogul')

   if (mogol_exe == None):
      print "mogul not found in path"
      return False
   else:
      comp_id = 'XXX'
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
      coot.write_pdb_from_mol(m_H, comp_id, "out.pdb")
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
              atom_type = atom.GetProp('atom_type')
              print "types: ", name, atom.GetSymbol(), atom_type, charge
          except KeyError:
              print "miss", name, atom.GetSymbol(), charge
      
      mogul_state = run_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name)
      
      if mogul_state:
          coot.mogul_out_to_mmcif_dict_by_mol(mogul_out_file_name, comp_id,
                                              compound_name, m_H, bor, mmcif_dict_name)
          return True
      return mogul_state


smiles_string = 'Cc1ccccc1'
smiles_string = 'Oc1ccccc1'
smiles_string = 'O=C(O)c1ccc(O)cc1'
smiles_string = 'NC(=O)c1ccc(O)cc1'
# smiles_string = 'O(c3cc4c2c1c[n+](ccc1ccc2nc4cc3)CCN5CCC(CC5)CCCC%10CCN(CC[n+]9cc8c7c6cc(OC)ccc6nc7ccc8cc9)CC%10)C'
# smiles_string = 'COc1ccc2[nH]c3ccc4cc[n+](CCN5CCC(CCCC6CCN(CC6)CC[n+]7ccc8ccc9[nH]c%10ccc(OC)cc%10c9c8c7)CC5)cc4c3c2c1'
# smiles_string = 'c1ccc2c(c1)cccc2' # napthalene, fused benze rings
# smiles_string = 'NC(=O)C(CO)N' # test for atom type C
smiles_string = 'c1[nH]c2c(c1)cccc2' # indole
sdf_file_name = 'test.sdf'
make_restraints(smiles_string, sdf_file_name, 'rest-test.cif')

