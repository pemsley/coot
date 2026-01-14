---
name: coot-unmodelled-blobs
description: "How to handle Unmodelled Density Blobs"
---
## Handling Unmodeled Density Blobs

### Investigating Blobs

When `find_blobs_py()` identifies unmodeled density, investigate what's actually there before acting:

```python
# Find blobs in difference map (molecule 2) at 3 sigma
blobs = coot.find_blobs_py(0, 2, 3.0)
# Returns: [[[x, y, z], score], ...]
# Higher score = larger/stronger blob

# Go to the biggest blob
if blobs:
    biggest = blobs[0]
    pos, score = biggest[0], biggest[1]
    coot.set_rotation_centre(pos[0], pos[1], pos[2])
```

### Determining What a Blob Represents

**Don't just look at which residues are "nearby" by CA distance** - this can be misleading for non-spherical blobs.

Instead, find which **atoms** are closest to the blob centre:

```python
blob_x, blob_y, blob_z = 59.92, 3.06, -4.23

import math
def dist(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

# Check all atoms near the blob
nearby_atoms = []
for chain in ['A', 'B']:
    for resno in range(1, 150):
        atoms = coot.residue_info_py(0, chain, resno, "")
        if atoms:
            res_name = coot.residue_name(0, chain, resno, "")
            for atom in atoms:
                x, y, z = atom[2]
                d = dist(blob_x, blob_y, blob_z, x, y, z)
                if d < 5.0:
                    atom_name = atom[0][0].strip()
                    nearby_atoms.append((d, chain, resno, res_name, atom_name))

nearby_atoms.sort()
for d, chain, resno, res_name, atom_name in nearby_atoms[:10]:
    print(f"{d:.1f}Å: {chain}/{resno} {res_name} {atom_name}")
```

### Chain Extension vs Terminal Atoms

**Critical distinction:**

If a blob is near the **C atom and O atom of the last residue** in a chain, it usually means:
- ❌ **NOT** a missing OXT (terminal carboxyl oxygen)
- ✅ **More residues to build** - the chain continues!

During model building, incomplete chains are common. The density beyond the last modeled residue represents **unbuilt residues**, not terminal atoms.

**Wrong approach:**
```python
# Don't do this for chain extension!
coot.add_OXT_to_residue(0, "A", 93, "")  # Wrong!
```

**Correct approach:**
```python
# Extend the chain by adding residues
checkpoint = coot.make_backup_checkpoint(0, "before chain extension")

# Add residues one at a time, refining as you go
coot.add_terminal_residue(0, "A", 93, "auto", 1)  # Adds residue 94
coot.refine_residues_py(0, [["A", 93, ""], ["A", 94, ""]])
coot.accept_moving_atoms_py()

# Check if blob is still there
blobs = coot.find_blobs_py(0, 2, 3.0)
# If blob persists (maybe smaller), add another residue

coot.add_terminal_residue(0, "A", 94, "auto", 1)  # Adds residue 95
coot.refine_residues_py(0, [["A", 94, ""], ["A", 95, ""]])
coot.accept_moving_atoms_py()

# Continue until blob is gone or no more density
```

### Iterative Chain Extension Workflow

```python
def extend_chain_into_density(imol, chain_id, last_resno, imol_map, max_residues=10):
    """Extend a chain into unmodeled density."""
    
    checkpoint = coot.make_backup_checkpoint(imol, f"before extending {chain_id}")
    
    current_resno = last_resno
    residues_added = 0
    
    for i in range(max_residues):
        # Check for remaining blob near current terminus
        blobs = coot.find_blobs_py(imol, imol_map, 3.0)
        if not blobs:
            break
            
        # Get position of current C-terminus
        atoms = coot.residue_info_py(imol, chain_id, current_resno, "")
        c_pos = None
        for atom in atoms:
            if atom[0][0].strip() == "C":
                c_pos = atom[2]
                break
        
        if not c_pos:
            break
            
        # Check if any blob is near the C-terminus
        blob_near_terminus = False
        for blob in blobs:
            pos = blob[0]
            d = ((pos[0]-c_pos[0])**2 + (pos[1]-c_pos[1])**2 + (pos[2]-c_pos[2])**2)**0.5
            if d < 6.0:  # Within 6 Angstroms
                blob_near_terminus = True
                break
        
        if not blob_near_terminus:
            break
        
        # Add next residue
        result = coot.add_terminal_residue(imol, chain_id, current_resno, "auto", 1)
        if result != 1:
            break
            
        current_resno += 1
        residues_added += 1
        
        # Refine the new region
        specs = [[chain_id, r, ""] for r in range(current_resno - 2, current_resno + 1) 
                 if r > 0]
        coot.refine_residues_py(imol, specs)
        coot.accept_moving_atoms_py()
        
        # Navigate to see progress
        coot.set_go_to_atom_chain_residue_atom_name(chain_id, current_resno, "CA")
    
    print(f"Added {residues_added} residues to chain {chain_id}")
    return residues_added

# Usage:
extend_chain_into_density(0, "A", 93, 2)
```

### When to Add OXT

Add OXT when you've finished building a fragment and there's no more density to extend into:

```python
# After extending chain A as far as the density allows
# Check there's no more blob near the terminus
blobs = coot.find_blobs_py(0, 2, 3.0)
# If no blob near C-terminus, cap the chain:
coot.add_OXT_to_residue(0, "A", 96, "")
coot.refine_residues_py(0, [["A", 95, ""], ["A", 96, ""]])
coot.accept_moving_atoms_py()
```

**The key distinction:**
- Blob near terminus → extend chain with `add_terminal_residue()`
- No blob, chain fully built → cap with `add_OXT_to_residue()`

### Summary: Blob Investigation Checklist

1. **Go to the blob** - `set_rotation_centre()`
2. **Find closest atoms** - not just closest residue CAs
3. **Check if near chain terminus** - look for C, O atoms of last residue
4. **If near terminus**: extend chain with `add_terminal_residue()`, don't add OXT
5. **If near side chain**: might be missing atoms, use `fill_partial_residue()`
6. **If isolated**: might be water, ion, or ligand
7. **Always checkpoint first** - `make_backup_checkpoint()`
8. **Refine after changes** - check if blob score decreases
