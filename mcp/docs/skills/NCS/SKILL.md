---
name: coot-NCS-reference-guidance
description: "NCS Reference Guidance"
---
## Using NCS/Symmetry-Related Molecules as Reference

### Overview

When building or rebuilding fragments, **use the sequence and model from NCS-related (non-crystallographically related) or molecular symmetry-related molecules** as a reference. If one copy is better built, use it to guide building of the other copies.

### Why This Helps

- Multiple copies of the same molecule in the asymmetric unit should have identical sequences
- One copy may have better density or be more complete than others
- The better-built copy provides the correct sequence and approximate coordinates

### Workflow: Compare Sequences First

Before building, compare sequences between related chains:

```python
def compare_chain_sequences(imol, chain1, chain2):
    """Compare sequences of two chains, return differences."""
    diffs = []
    for resno in range(1, 200):  # Adjust range as needed
        name1 = coot.residue_name(imol, chain1, resno, "")
        name2 = coot.residue_name(imol, chain2, resno, "")
        
        if name1 and name2 and name1 != name2:
            diffs.append((resno, name1, name2))
        elif name1 and not name2:
            diffs.append((resno, name1, "MISSING"))
        elif name2 and not name1:
            diffs.append((resno, "MISSING", name2))
    
    return diffs

# Compare chains A and B
diffs = compare_chain_sequences(0, "A", "B")
if diffs:
    print("Sequence differences:")
    for resno, a, b in diffs:
        print(f"  {resno}: A={a}, B={b}")
else:
    print("Sequences are identical")
```

### Using Better Chain as Template

If chain B is more complete, use it to guide chain A:

```python
# Chain B has residues 94-96 built, chain A doesn't
# Get sequence from chain B
target_seq = ""
for resno in range(94, 97):
    name = coot.residue_name(0, "B", resno, "")
    if name:
        # Convert 3-letter to 1-letter code
        one_letter = {'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 
                      'ILE': 'I', 'PRO': 'P', 'PHE': 'F', 'TYR': 'Y',
                      'TRP': 'W', 'SER': 'S', 'THR': 'T', 'CYS': 'C',
                      'MET': 'M', 'ASN': 'N', 'GLN': 'Q', 'ASP': 'D',
                      'GLU': 'E', 'LYS': 'K', 'ARG': 'R', 'HIS': 'H'}
        target_seq += one_letter.get(name, 'X')

print(f"Target sequence from B: {target_seq}")

# Build residues in chain A, then mutate to correct sequence
# ... add_terminal_residue calls ...

# Mutate to match chain B
coot.mutate_and_autofit_residue_range(0, "A", 94, 96, target_seq)
```

### Avoiding Common Mistakes

1. **Don't assume "auto" gives correct residue type**
   - `add_terminal_residue(imol, chain, resno, "auto", 1)` builds ALA by default
   - Always mutate to the correct sequence afterward

2. **Check for subtle differences**
   - GLY vs GLN both start with 'G' - compare 3-letter codes, not single letters
   - Use `residue_name()` not single-letter conversions for comparisons

3. **Verify after mutation**
   - Check density correlation after mutating
   - Use `auto_fit_best_rotamer()` to optimize side chain placement

### Complete Example: Build Fragment Using NCS Reference

```python
def build_fragment_from_ncs(imol, target_chain, ref_chain, start_resno, end_resno):
    """Build missing residues in target_chain using ref_chain as template."""
    
    checkpoint = coot.make_backup_checkpoint(imol, 
        f"before building {target_chain}/{start_resno}-{end_resno} from {ref_chain}")
    
    # Get reference sequence
    ref_seq = ""
    for resno in range(start_resno, end_resno + 1):
        name = coot.residue_name(imol, ref_chain, resno, "")
        if name:
            one_letter = {'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 
                          'ILE': 'I', 'PRO': 'P', 'PHE': 'F', 'TYR': 'Y',
                          'TRP': 'W', 'SER': 'S', 'THR': 'T', 'CYS': 'C',
                          'MET': 'M', 'ASN': 'N', 'GLN': 'Q', 'ASP': 'D',
                          'GLU': 'E', 'LYS': 'K', 'ARG': 'R', 'HIS': 'H'}
            ref_seq += one_letter.get(name, 'X')
    
    print(f"Reference sequence from {ref_chain}: {ref_seq}")
    
    # Find last existing residue in target chain
    last_resno = start_resno - 1
    while coot.residue_name(imol, target_chain, last_resno, ""):
        pass
    last_resno -= 1
    
    # Extend chain
    for i in range(end_resno - last_resno):
        current = last_resno + i
        result = coot.add_terminal_residue(imol, target_chain, current, "auto", 1)
        if result != 1:
            print(f"Failed to add residue after {current}")
            break
        
        # Refine as we go
        specs = [[target_chain, r, ""] for r in range(max(1, current-1), current+2)]
        coot.refine_residues_py(imol, specs)
        coot.accept_moving_atoms_py()
    
    # Mutate to correct sequence
    coot.mutate_and_autofit_residue_range(imol, target_chain, 
                                          start_resno, end_resno, ref_seq)
    
    # Final refinement
    specs = [[target_chain, r, ""] for r in range(start_resno, end_resno + 1)]
    coot.refine_residues_py(imol, specs)
    coot.accept_moving_atoms_py()
    
    # Verify
    print(f"\nBuilt residues in {target_chain}:")
    for resno in range(start_resno, end_resno + 1):
        name = coot.residue_name(imol, target_chain, resno, "")
        print(f"  {resno}: {name}")

# Usage:
build_fragment_from_ncs(0, "A", "B", 94, 96)
```

## Using NCS Matrices to Copy Fragments

For more substantial rebuilding, you can copy entire fragments from one NCS copy to another using the NCS transformation matrix.

### Step 1: Enable NCS Ghosts and Get the Matrix

```python
# Create NCS ghosts (calculates the transformation matrices)
coot.make_ncs_ghosts_maybe(imol)

# Get NCS information
ncs_chains = coot.ncs_chain_ids_py(imol)  # e.g., [['A', 'B']]
ghosts = coot.ncs_ghosts_py(imol)

# Extract the rotation matrix and translation vector
# ghosts[0] = [description, from_chain, to_chain, [rot_matrix, trans_vector], bool]
rot_trans = ghosts[0][3]
rot = rot_trans[0]    # 9 elements: m11,m12,m13,m21,m22,m23,m31,m32,m33
trans = rot_trans[1]  # 3 elements: x,y,z

print(f"NCS transforms {ghosts[0][1]} onto {ghosts[0][2]}")
```

### Step 2: Copy Fragment from Reference Chain

```python
# Delete the problem region in target chain
for resno in range(1, 11):
    coot.delete_residue(imol, "A", resno, "")

# Create a fragment molecule from the reference chain
frag_mol = coot.new_molecule_by_atom_selection(imol, "//B/1-10")
```

### Step 3: Apply NCS Transformation

```python
# Transform the fragment to the target position
coot.transform_molecule_by(frag_mol,
    rot[0], rot[1], rot[2],   # row 1 of rotation matrix
    rot[3], rot[4], rot[5],   # row 2
    rot[6], rot[7], rot[8],   # row 3
    trans[0], trans[1], trans[2])  # translation
```

### Step 4: Change Chain ID and Replace Fragment

```python
# Change chain ID from B to A
coot.change_chain_id(frag_mol, "B", "A", 0, 0, 0)

# Paste into the target molecule
coot.replace_fragment(imol, frag_mol, "//A/1-10")
```

### Step 5: Refine with Flanking Residues

```python
# Include flanking residues for better annealing at the junction
specs = [["A", r, ""] for r in range(1, 14)]  # 1-10 plus 11-13
coot.refine_residues_py(imol, specs)
coot.accept_moving_atoms_py()
```

### Complete Function: NCS Fragment Copy

```python
def copy_fragment_by_ncs(imol, from_chain, to_chain, start_resno, end_resno, flank=3):
    """Copy a fragment from one NCS chain to another using NCS matrix.
    
    Parameters
    ----------
    imol : int
        Molecule index
    from_chain : str
        Source chain ID (the better-built one)
    to_chain : str
        Target chain ID (the one to rebuild)
    start_resno : int
        First residue number to copy
    end_resno : int
        Last residue number to copy
    flank : int
        Number of flanking residues to include in refinement
    """
    
    checkpoint = coot.make_backup_checkpoint(imol, 
        f"before NCS copy {from_chain}/{start_resno}-{end_resno} to {to_chain}")
    
    # Ensure NCS ghosts are calculated
    coot.make_ncs_ghosts_maybe(imol)
    
    # Get NCS matrix
    ghosts = coot.ncs_ghosts_py(imol)
    if not ghosts or ghosts == False:
        print("No NCS ghosts found!")
        return False
    
    # Find the right ghost (from_chain -> to_chain)
    rot_trans = None
    for ghost in ghosts:
        if ghost[1] == from_chain and ghost[2] == to_chain:
            rot_trans = ghost[3]
            break
        elif ghost[1] == to_chain and ghost[2] == from_chain:
            # Need inverse - for now, skip this case
            print(f"Found inverse NCS ({to_chain} -> {from_chain}), may need inverse_rtop_py")
            rot_trans = ghost[3]
            break
    
    if not rot_trans:
        print(f"No NCS relationship found between {from_chain} and {to_chain}")
        return False
    
    rot = rot_trans[0]
    trans = rot_trans[1]
    
    # Delete target residues
    for resno in range(start_resno, end_resno + 1):
        coot.delete_residue(imol, to_chain, resno, "")
    
    # Create fragment from source chain
    atom_sel = f"//{from_chain}/{start_resno}-{end_resno}"
    frag_mol = coot.new_molecule_by_atom_selection(imol, atom_sel)
    
    if frag_mol < 0:
        print(f"Failed to create fragment from {atom_sel}")
        return False
    
    # Transform fragment
    coot.transform_molecule_by(frag_mol,
        rot[0], rot[1], rot[2],
        rot[3], rot[4], rot[5],
        rot[6], rot[7], rot[8],
        trans[0], trans[1], trans[2])
    
    # Change chain ID
    coot.change_chain_id(frag_mol, from_chain, to_chain, 0, 0, 0)
    
    # Replace fragment in target molecule
    target_sel = f"//{to_chain}/{start_resno}-{end_resno}"
    result = coot.replace_fragment(imol, frag_mol, target_sel)
    
    if result != 1:
        print("replace_fragment failed")
        coot.restore_to_backup_checkpoint(imol, checkpoint)
        return False
    
    # Refine with flanking residues
    refine_start = max(1, start_resno - flank)
    refine_end = end_resno + flank
    specs = [[to_chain, r, ""] for r in range(refine_start, refine_end + 1)
             if coot.residue_name(imol, to_chain, r, "")]
    
    coot.set_go_to_atom_chain_residue_atom_name(to_chain, start_resno, "CA")
    coot.refine_residues_py(imol, specs)
    coot.accept_moving_atoms_py()
    
    # Clean up fragment molecule
    coot.close_molecule(frag_mol)
    
    print(f"Successfully copied {from_chain}/{start_resno}-{end_resno} to {to_chain}")
    return True

# Usage:
copy_fragment_by_ncs(0, "B", "A", 1, 10, flank=3)
```

### When to Use NCS Fragment Copy

- **Large missing regions** - faster than building residue-by-residue
- **Poorly built regions** - replace with better NCS copy
- **After molecular replacement** - one copy may fit better initially
- **Symmetric assemblies** - enforce consistency across copies

### Key Takeaway

**When multiple copies of the same molecule exist, use the better-built copy as a reference for sequence and model building.** This ensures:
- Correct residue types (not default ALA)
- Consistent sequences across all copies
- Proper disulfide bonds and other features

**For larger rebuilds, use the NCS matrix to transform and paste entire fragments** - this is faster and more accurate than building from scratch.
