---
name: coot-rdkit
description: RDKit molecular manipulation and visualization within Coot's Python environment. Use when working with Coot and need to (1) Create RDKit molecules from Coot monomers, (2) Modify molecular structures (e.g., atom substitution), (3) Generate 2D chemical structure diagrams, (4) Perform cheminformatics operations on ligands or small molecules loaded in Coot.
---

# Coot-RDKit Integration

This skill provides guidance for using RDKit within Coot's Python environment for molecular manipulation and visualization.

## Key Integration Points

### Module Import
Use `coot_headless_api` (NOT `chapi`):
```python
import coot_headless_api
from rdkit import Chem
from rdkit.Chem import AllChem
```

### Creating RDKit Molecules from Coot Monomers

```python
import coot_headless_api
import base64
from rdkit import Chem

# Initialize molecules container
molecules = coot_headless_api.molecules_container_t(False)  # False = not verbose

# Get monomer from Coot's library
imol = molecules.get_monomer("AMP")  # or any other monomer code

# Get RDKit molecule as pickled base64
pickle_base64_str = molecules.get_rdkit_mol_pickle_base64("AMP", imol)

# Decode and create RDKit molecule
pickle_bytes = base64.b64decode(pickle_base64_str)
rdkit_mol = Chem.Mol(pickle_bytes)  # Use Chem.Mol(), NOT pickle.loads()
```

## Molecular Manipulation

### Atom Substitution
```python
from rdkit import Chem

# Make editable copy
mol_edit = Chem.RWMol(rdkit_mol)

# Replace atom (e.g., phosphorus to sulfur)
for atom in mol_edit.GetAtoms():
    if atom.GetSymbol() == 'P':
        atom.SetAtomicNum(16)  # 16 = sulfur
        break

# Convert back to read-only molecule
modified_mol = mol_edit.GetMol()
Chem.SanitizeMol(modified_mol)
```

## 2D Structure Visualization

### CRITICAL: Always Regenerate 2D Coordinates

When removing hydrogens or modifying structure, ALWAYS regenerate 2D coordinates:

```python
from rdkit.Chem import AllChem

# Remove hydrogens
mol_no_h = Chem.RemoveHs(mol)

# IMPORTANT: Regenerate 2D coords AFTER removing hydrogens
AllChem.Compute2DCoords(mol_no_h)

# Now generate visualization
```

### Generating SVG Diagrams

```python
from rdkit.Chem.Draw import rdMolDraw2D

drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
drawer.DrawMolecule(mol_no_h)
drawer.FinishDrawing()
svg_string = drawer.GetDrawingText()

# Save without displaying (see Data Handling below)
```

## Data Handling Best Practices

### NEVER Display Large String Data

Do NOT return or print large strings (SVG, base64, etc.) as this causes slow response times:

**BAD:**
```python
svg_data  # This displays all the text - SLOW!
```

**GOOD:**
```python
# Just save directly without displaying
# (use len() to verify if needed, though even this may not return properly)
```

### Efficient File Writing

Write files directly without displaying content:

**From Coot Python:**
```python
svg_content = drawer.GetDrawingText()
# Don't display svg_content - just reference it
```

**Then in bash or file creation:**
```bash
# Use the variable directly without echoing/catting the content
```

## Coot Python Limitations

### Single-Line Return Values Only

Coot's Python environment only returns values from single-line expressions:

**Works:**
```python
Chem.MolToSmiles(mol)  # Returns SMILES string
```

**Doesn't return properly:**
```python
x = 5
y = 10
x + y  # Won't return the value
```

**Workaround** - Define function in one call, execute in next:
```python
# Call 1: Define
def my_function():
    x = 5
    y = 10
    return x + y

# Call 2: Execute
my_function()  # Now returns 15
```

## Common Workflows

### Modify and Visualize Ligand

1. Load monomer from Coot library
2. Convert to RDKit molecule
3. Make modifications (atom substitution, etc.)
4. Remove hydrogens if desired for cleaner diagram
5. **Regenerate 2D coordinates** (critical!)
6. Generate SVG diagram
7. Save to file **without displaying**
