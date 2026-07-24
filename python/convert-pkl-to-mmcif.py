
"""Write a minimal CCD-style mmCIF (acedrg input) directly from an RDKit mol,
using gemmi. Atom names come from the mol's `name` prop."""
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

fn  = "thing.pkl"
out = "thing.cif"

mol = Chem.Mol(open(fn, "rb").read())

# The compound id (three-letter-code) rides along on the molecule as the
# `ResName` property (set in Coot/Layla, preserved through the AllProps pickle).
# Use it as the CCD id so acedrg reads it from _chem_comp.id and no -r is needed.
# Fall back to LIG only if the mol was never given a name.
ccd_id = mol.GetProp("ResName").strip() if mol.HasProp("ResName") else ""
if not ccd_id:
    ccd_id = "LIG"
    print("warning: mol has no ResName property; falling back to id 'LIG'")

# Capture Layla's own 2D sketch BEFORE we touch the mol. This might be the
# hand-drawn layout and is the depiction we care about most, and it is written
# out verbatim into the _pdbe_chem_comp_*_depiction loops below (rather
# than a recomputed 2D layout).
mol2d = Chem.Mol(mol)
if mol2d.GetNumConformers() == 0:
    raise RuntimeError("input mol has no conformer to use as a 2D depiction")

# The conformer coming out of Layla is a 2D sketch layout (all z = 0),
# Build a genuine 3D conformer with ETKDG. acedrg
# regenerates its own 3D anyway, but this makes the CIF's *_ideal coords sane.
mol.RemoveAllConformers()
mol = Chem.AddHs(mol)                       # add Hs *before* embedding
params = AllChem.ETKDGv3()
params.randomSeed = 1
if AllChem.EmbedMolecule(mol, params) != 0:
    # Fallback for awkward molecules: relax the strict experimental-torsion terms.
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("3D embedding failed for " + fn)
AllChem.MMFFOptimizeMolecule(mol)           # light clean-up; ignore if no MMFF params
conf = mol.GetConformer()
assert conf.Is3D()

# Match the PDB-CCD convention: aromatic ring bonds carry a KEKULIZED order
# (SING/DOUB) plus pdbx_aromatic_flag = Y — not a literal "AROM" value_order.
# clearAromaticFlags=False keeps GetIsAromatic() true so we can still set the flag.
Chem.Kekulize(mol, clearAromaticFlags=False)

# --- atom names --------------------------------------------------------------
# Preserving the Coot/Layla atom names is the whole point of this pipeline, so
# names carried on the mol take priority. Atoms without a name (the Hs added by
# AddHs above) get a generated element+counter name. CCD atom_id must be UNIQUE,
# so generated names are checked against every name already in use.
names = [None] * mol.GetNumAtoms()
used = set()

# Pass 1: take the existing (stripped) names. Guard against duplicates/blanks.
for i, a in enumerate(mol.GetAtoms()):
    if a.HasProp("name"):
        nm = a.GetProp("name").strip()
        if nm and nm not in used:
            names[i] = nm
            used.add(nm)
        elif nm in used:
            raise RuntimeError(f"duplicate atom name '{nm}' on the input mol "
                               f"(atom idx {i}); names must be unique")

# Pass 2: generate unique names for the rest (the added Hs), skipping collisions.
elem_counts = {}
for i, a in enumerate(mol.GetAtoms()):
    if names[i] is not None:
        continue
    sym = a.GetSymbol().upper()
    while True:
        elem_counts[sym] = elem_counts.get(sym, 0) + 1
        cand = f"{sym}{elem_counts[sym]}"
        if cand not in used:
            break
    names[i] = cand
    used.add(cand)

n_from_mol = sum(1 for a in mol.GetAtoms() if a.HasProp("name") and a.GetProp("name").strip())
assert len(set(names)) == len(names), "atom names are not unique"

ORDER = {Chem.BondType.SINGLE: "SING", Chem.BondType.DOUBLE: "DOUB",
         Chem.BondType.TRIPLE: "TRIP", Chem.BondType.AROMATIC: "AROM"}

# --- build the mmCIF ---------------------------------------------------------
doc = gemmi.cif.Document()
blk = doc.add_new_block(ccd_id)
blk.set_pair("_chem_comp.id", ccd_id)
blk.set_pair("_chem_comp.three_letter_code", ccd_id)
blk.set_pair("_chem_comp.name", gemmi.cif.quote(Chem.MolToSmiles(mol)))
blk.set_pair("_chem_comp.type", "non-polymer")
blk.set_pair("_chem_comp.formula", gemmi.cif.quote(Chem.rdMolDescriptors.CalcMolFormula(mol)))

atom_loop = blk.init_loop("_chem_comp_atom.", [
    "comp_id", "atom_id", "type_symbol", "charge",
    "pdbx_model_Cartn_x_ideal", "pdbx_model_Cartn_y_ideal",
    "pdbx_model_Cartn_z_ideal"])
for i, a in enumerate(mol.GetAtoms()):
    p = conf.GetAtomPosition(i)
    atom_loop.add_row([ccd_id, names[i], a.GetSymbol(), str(a.GetFormalCharge()),
                       f"{p.x:.3f}", f"{p.y:.3f}", f"{p.z:.3f}"])

bond_loop = blk.init_loop("_chem_comp_bond.", [
    "comp_id", "atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag"])
for b in mol.GetBonds():
    order = ORDER.get(b.GetBondType(), "SING")
    arom = "Y" if b.GetIsAromatic() else "N"
    bond_loop.add_row([ccd_id, names[b.GetBeginAtomIdx()],
                       names[b.GetEndAtomIdx()], order, arom])

# --- 2D depiction loops (PDBe extension) -------------------------------------
# Written from Layla's hand-drawn layout (mol2d), NOT a recomputed one. Mirrors
# pdbeccdutils' _pdbe_chem_comp_atom_depiction / _bond_depiction categories.
# mol2d holds only the drawn heavy atoms (no explicit H), which is exactly what
# a 2D depiction shows; every atom already carries its Coot/Layla name.
def dname(atom):
    return atom.GetProp("name").strip()

conf2d = mol2d.GetConformer()
atom_dep = blk.init_loop("_pdbe_chem_comp_atom_depiction.", [
    "comp_id", "atom_id", "element", "model_Cartn_x", "model_Cartn_y",
    "pdbx_ordinal"])
for i, a in enumerate(mol2d.GetAtoms()):
    p = conf2d.GetAtomPosition(i)
    atom_dep.add_row([ccd_id, dname(a), a.GetSymbol(),
                      f"{p.x:.3f}", f"{p.y:.3f}", str(i + 1)])

# Wedge/hash bond directions for the drawing, derived from the 2D layout +
# chirality (the same PrepareMolForDrawing call pdbeccdutils uses). value_order
# is the kekulized RDKit bond-type name (SINGLE/DOUBLE/...), bond_dir the wedge
# enum name (NONE/BEGINWEDGE/BEGINDASH). Bonds to H are omitted.
try:
    drawmol = rdMolDraw2D.PrepareMolForDrawing(
        mol2d, wedgeBonds=True, kekulize=True, addChiralHs=True)
except (RuntimeError, ValueError):
    try:
        drawmol = rdMolDraw2D.PrepareMolForDrawing(
            mol2d, wedgeBonds=False, kekulize=True, addChiralHs=False)
    except (RuntimeError, ValueError):
        drawmol = Chem.Mol(mol2d)
        Chem.Kekulize(drawmol, clearAromaticFlags=True)

bond_dep = blk.init_loop("_pdbe_chem_comp_bond_depiction.", [
    "comp_id", "atom_id_1", "atom_id_2", "value_order", "bond_dir",
    "pdbx_ordinal"])
dep_ordinal = 0
for b in drawmol.GetBonds():
    a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
    if a1.GetSymbol() == "H" or a2.GetSymbol() == "H":
        continue
    dep_ordinal += 1
    bond_dep.add_row([ccd_id, dname(a1), dname(a2),
                      b.GetBondType().name, b.GetBondDir().name, str(dep_ordinal)])

doc.write_file(out)
print(f"wrote {out}: {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds "
      f"({n_from_mol} names preserved from mol, "
      f"{mol.GetNumAtoms() - n_from_mol} generated); "
      f"2D depiction: {mol2d.GetNumAtoms()} atoms, {dep_ordinal} bonds")
