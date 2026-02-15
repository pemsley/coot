---
name: coot-figure-making
description: "Best practices for creating publication-quality molecular graphics figures in Coot using user-defined colors, ribbons, and molecular representations"
---

# Coot Figure-Making Best Practices

This skill provides guidance for creating publication-quality molecular graphics figures in Coot, with emphasis on custom coloring schemes and ribbon representations.

## Critical Rule: Use Function Wrappers

**ALWAYS wrap complex multi-line code in functions when using the MCP interface.**

### Why?

Multi-line code executed directly (without a function wrapper) can corrupt the Python interpreter if any error occurs. Function wrappers provide proper error handling and prevent interpreter crashes.

### Pattern

```python
# Step 1: Define the function with run_python_multiline()
def my_figure_function():
    # ... complex code here ...
    return result

# Step 2: Call the function with run_python()
my_figure_function()
```

### Example

```python
# ❌ BAD - Direct multi-line execution (can crash interpreter)
ss_info = coot.get_header_secondary_structure_info(0)
strands = ss_info['strands']
# ... more code ...

# ✅ GOOD - Function wrapper (safe error handling)
def setup_figure():
    ss_info = coot.get_header_secondary_structure_info(0)
    strands = ss_info['strands']
    # ... more code ...
    return "Done"

setup_figure()
```

## User-Defined Color Workflow

### Overview

Coot's user-defined color system allows custom coloring of specific selections. The workflow has strict ordering requirements:

1. Define color palette with `set_user_defined_colours_py()`
2. Assign colors to selections with `set_user_defined_atom_colour_by_selection_py()`
3. Create representation that uses those colors
4. **Critical**: Recreating a representation requires reassigning ALL colors

### Color Indices

- Indices 0-59: Reserved for Coot's internal colors
- Indices 60+: Available for user-defined colors
- Convention: Start at 60 and increment (60, 61, 62, etc.)

### Basic Pattern

```python
def setup_colors():
    # Step 1: Define colors (RGB values 0.0-1.0)
    blue = [0.3, 0.6, 1.0]
    orange = [1.0, 0.5, 0.0]
    
    # Step 2: Set color palette
    coot.set_user_defined_colours_py([
        (60, blue),
        (61, orange)
    ])
    
    # Step 3: Assign colors to selections (MMDB format)
    color_assignments = [
        ("//A/10-50", 60),   # Blue for residues 10-50
        ("//A/100-150", 61)  # Orange for residues 100-150
    ]
    coot.set_user_defined_atom_colour_by_selection_py(imol, color_assignments)
    
    # Step 4: Create representation
    coot.add_ribbon_representation_with_user_defined_colours(imol, "My Figure")
    
    return "Colors applied"
```

## Secondary Structure Coloring

### Using PDB Header Information

The most reliable way to color by secondary structure is to use the annotations in the PDB header.

```python
def color_by_secondary_structure(imol):
    # Get secondary structure from PDB header
    ss_info = coot.get_header_secondary_structure_info(imol)
    strands = ss_info['strands']
    helices = ss_info['helices']
    
    # Build MMDB selection strings for strands
    strand_selections = []
    for strand in strands:
        chain = strand['initChainID']
        start = strand['initSeqNum']
        end = strand['endSeqNum']
        selection = f"//{chain}/{start}-{end}"
        strand_selections.append(selection)
    
    # Build MMDB selection strings for helices
    helix_selections = []
    for helix in helices:
        chain = helix['initChainID']
        start = helix['initSeqNum']
        end = helix['endSeqNum']
        selection = f"//{chain}/{start}-{end}"
        helix_selections.append(selection)
    
    # Define colors
    blue = [0.3, 0.6, 1.0]        # Beta strands
    purple = [0.5, 0.3, 0.5]      # Helices
    
    coot.set_user_defined_colours_py([
        (60, blue),
        (61, purple)
    ])
    
    # Assign colors
    strand_assignments = [(sel, 60) for sel in strand_selections]
    helix_assignments = [(sel, 61) for sel in helix_selections]
    all_assignments = strand_assignments + helix_assignments
    
    coot.set_user_defined_atom_colour_by_selection_py(imol, all_assignments)
    
    # Create ribbon
    coot.add_ribbon_representation_with_user_defined_colours(imol, "Secondary Structure")
    
    return f"{len(strand_selections)} strands, {len(helix_selections)} helices"
```

### Adding Manual Secondary Structure Annotations

Sometimes secondary structure elements aren't annotated in the PDB header. You can add them manually:

```python
def add_missing_helix(imol):
    # Get existing colors
    ss_info = coot.get_header_secondary_structure_info(imol)
    strands = ss_info['strands']
    helices = ss_info['helices']
    
    # Build all selections (as before)
    strand_selections = [...]
    helix_selections = [...]
    
    # Add manually identified helix
    helix_selections.append("//A/135-139")
    
    # Reassign ALL colors and recreate ribbon
    # (must include ALL selections every time)
    coot.set_user_defined_colours_py([...])
    all_assignments = strand_assignments + helix_assignments
    coot.set_user_defined_atom_colour_by_selection_py(imol, all_assignments)
    coot.add_ribbon_representation_with_user_defined_colours(imol, "Updated")
    
    return "Helix added"
```

### Critical Rule: Reassign All Colors When Recreating

**When you recreate a ribbon representation, you MUST reassign ALL color selections, not just the new ones.**

```python
# ❌ BAD - Only assigns new helix, strands lose their color
coot.set_user_defined_atom_colour_by_selection_py(imol, [("//A/135-139", 61)])
coot.add_ribbon_representation_with_user_defined_colours(imol, "New")

# ✅ GOOD - Reassigns everything
all_assignments = strand_assignments + helix_assignments + new_helix
coot.set_user_defined_atom_colour_by_selection_py(imol, all_assignments)
coot.add_ribbon_representation_with_user_defined_colours(imol, "New")
```

## Highlighting Specific Features

### Extracting Residues to Separate Molecules

To highlight specific residues (like active site residues, ligands, chromophores), extract them to a new molecule:

```python
def highlight_feature(imol, selection, color_rgb, bond_thickness=10.0):
    # Extract to new molecule
    feature_imol = coot.new_molecule_by_atom_selection(imol, selection)
    
    if not coot.is_valid_model_molecule(feature_imol):
        return -1
    
    # Define color
    color_index = 62  # Use a different index than strands/helices
    coot.set_user_defined_colours_py([(color_index, color_rgb)])
    
    # Assign color
    coot.set_user_defined_atom_colour_by_selection_py(feature_imol, [(selection, color_index)])
    
    # Add representation with thick bonds
    coot.add_molecular_representation_py(
        feature_imol,
        selection,
        "userDefined",  # Use user-defined colors
        "Bonds"
    )
    
    # Make bonds thicker for emphasis
    coot.set_bond_thickness(feature_imol, bond_thickness)
    
    return feature_imol

# Example: Highlight chromophore in orange
highlight_feature(0, "//A/66", [1.0, 0.5, 0.0], 10.0)
```

## View Setup

### Centering and Zooming

```python
def setup_view(chain_id, resno, zoom_level=200):
    # Center on specific residue
    coot.set_go_to_atom_chain_residue_atom_name(chain_id, resno, "CA")
    
    # Set zoom level
    # 150-300: Whole molecule overview
    # 50-100: Domain level
    # 20-50: Residue detail
    coot.set_zoom(zoom_level)
```

### Hiding Bond Representation

For ribbon-only figures, hide the bond representation:

```python
# Hide bonds for molecule 0
coot.set_mol_displayed(0, 0)

# Show bonds again if needed
coot.set_mol_displayed(0, 1)
```

## Complete Example: GFP Beta Barrel Figure

This example creates a publication-quality figure showing GFP's beta barrel structure with colored secondary structure and highlighted chromophore.

```python
def make_gfp_figure():
    """
    Create a figure showing GFP with:
    - Blue beta barrel strands
    - Dark pastel helices
    - Orange chromophore with thick bonds
    """
    imol = 0  # GFP molecule
    
    # Get secondary structure
    ss_info = coot.get_header_secondary_structure_info(imol)
    strands = ss_info['strands']
    helices = ss_info['helices']
    
    # Build strand selections
    strand_selections = []
    for strand in strands:
        sel = f"//{strand['initChainID']}/{strand['initSeqNum']}-{strand['endSeqNum']}"
        strand_selections.append(sel)
    
    # Build helix selections
    helix_selections = []
    for helix in helices:
        sel = f"//{helix['initChainID']}/{helix['initSeqNum']}-{helix['endSeqNum']}"
        helix_selections.append(sel)
    
    # Add manually identified helix (not in PDB header)
    helix_selections.append("//A/135-139")
    
    # Define colors
    blue = [0.3, 0.6, 1.0]           # Beta strands
    dark_pastel = [0.5, 0.3, 0.5]    # Helices
    orange = [1.0, 0.5, 0.0]         # Chromophore
    
    coot.set_user_defined_colours_py([
        (60, blue),
        (61, dark_pastel),
        (62, orange)
    ])
    
    # Assign colors to secondary structure
    strand_assignments = [(sel, 60) for sel in strand_selections]
    helix_assignments = [(sel, 61) for sel in helix_selections]
    all_assignments = strand_assignments + helix_assignments
    coot.set_user_defined_atom_colour_by_selection_py(imol, all_assignments)
    
    # Create ribbon
    coot.add_ribbon_representation_with_user_defined_colours(imol, "GFP Barrel")
    
    # Hide bonds
    coot.set_mol_displayed(imol, 0)
    
    # Extract and highlight chromophore
    chrom_imol = coot.new_molecule_by_atom_selection(imol, "//A/66")
    coot.set_user_defined_atom_colour_by_selection_py(chrom_imol, [("//A/66", 62)])
    coot.add_molecular_representation_py(chrom_imol, "//A/66", "userDefined", "Bonds")
    coot.set_bond_thickness(chrom_imol, 10.0)
    
    # Center view
    coot.set_go_to_atom_chain_residue_atom_name("A", 100, "CA")
    coot.set_zoom(200)
    
    return f"Figure created: {len(strand_selections)} strands, {len(helix_selections)} helices, chromophore"

# To use (must be called with run_python after defining with run_python_multiline):
make_gfp_figure()
```

## Color Scheme Suggestions

### Standard Secondary Structure
- **Beta strands**: Blue `[0.3, 0.6, 1.0]`
- **Alpha helices**: Red/Purple `[0.8, 0.2, 0.4]` or `[0.5, 0.3, 0.5]`
- **Loops**: Gray (or leave uncolored)

### Highlight Schemes
- **Active site**: Bright orange `[1.0, 0.5, 0.0]`
- **Substrate binding**: Yellow `[1.0, 0.9, 0.0]`
- **Metal coordination**: Cyan `[0.0, 0.8, 0.8]`
- **Mutation sites**: Magenta `[1.0, 0.0, 1.0]`

### Domain Coloring
- **Domain 1**: Blue `[0.2, 0.4, 0.8]`
- **Domain 2**: Green `[0.2, 0.8, 0.4]`
- **Domain 3**: Orange `[0.9, 0.5, 0.2]`
- **Linker**: Gray `[0.6, 0.6, 0.6]`

## Troubleshooting

### Colors Don't Appear
**Problem**: Ribbon is gray after setting colors.

**Solution**: Make sure you call `add_ribbon_representation_with_user_defined_colours()` AFTER setting colors.

### Colors Disappear After Update
**Problem**: Added new colored region, but existing colors turned red/brown.

**Solution**: When recreating ribbon, reassign ALL color selections, not just new ones.

### Python Interpreter Crashes
**Problem**: Multi-line code causes "Failed to get __main__ module" error.

**Solution**: Always use function wrappers with `run_python_multiline()` then call with `run_python()`.

### Feature Not Visible
**Problem**: Extracted feature (ligand, chromophore) doesn't show up.

**Solution**: 
1. Check molecule is valid: `coot.is_valid_model_molecule(feature_imol)`
2. Ensure `add_molecular_representation_py()` succeeded
3. Verify feature molecule is displayed: `coot.set_mol_displayed(feature_imol, 1)`

## Graphics Quality Settings for Publication Figures

For publication-quality figures, especially for journal covers or high-impact visualizations, use these graphics settings:

### Background Color

Set an appropriate background color for your publication medium:

```python
# Light grey (80%) - excellent for print publications
coot.set_background_colour(0.8, 0.8, 0.8)

# Near-white (98%) - for very light backgrounds
coot.set_background_colour(0.98, 0.98, 0.98)

# Medium grey (50%) - good general purpose
coot.set_background_colour(0.5, 0.5, 0.5)

# White - for manuscripts requiring white backgrounds
coot.set_background_colour(1.0, 1.0, 1.0)

# Black - for dark backgrounds (presentations)
coot.set_background_colour(0.0, 0.0, 0.0)
```

### Outline Mode

Enable outline mode (also called "cel shading" or "toon shading") for a polished, professional look with dark edges around ribbons and bonds:

```python
# Enable outline mode
coot.set_use_outline(1)

# Disable outline mode
coot.set_use_outline(0)

# Query outline state
state = coot.use_outline_state()
```

### Fancy Graphics Mode

Enable advanced rendering effects for high-quality figures:

```python
def enable_fancy_graphics():
    """Enable all fancy graphics effects for publication figures"""
    
    # Ambient Occlusion (SSAO) - adds subtle shadows in crevices
    # Makes surfaces appear more 3D with depth perception
    coot.set_use_ambient_occlusion(1)
    
    # Fancy Lighting - enhanced lighting model
    # Provides better shading and highlights
    coot.set_use_fancy_lighting(1)
    
    # Depth Blur - depth of field effect
    # Blurs distant objects for focus effect
    coot.set_use_depth_blur(1)
    
    return "Fancy graphics enabled"

def disable_fancy_graphics():
    """Disable fancy graphics for faster rendering"""
    coot.set_use_ambient_occlusion(0)
    coot.set_use_fancy_lighting(0)
    coot.set_use_depth_blur(0)
    return "Fancy graphics disabled"
```

### SSAO Fine-Tuning

Ambient occlusion can be fine-tuned for different effects:

```python
# Adjust SSAO strength (default: typically around 1.0)
coot.set_ssao_strength(1.5)  # Stronger shadows

# Adjust SSAO radius (default: typically around 0.5)
coot.set_ssao_radius(0.7)  # Larger shadow radius

# Adjust SSAO bias (default: typically around 0.025)
coot.set_ssao_bias(0.03)  # Reduces shadow artifacts

# Set number of samples for SSAO (more = better quality, slower)
coot.set_ssao_kernel_n_samples(32)  # Default is often 16

# Set blur size (0, 1, or 2)
coot.set_ssao_blur_size(1)  # Smooths out SSAO shadows
```

### Shadow Settings

Coot provides real-time shadow rendering that adds depth and dimensionality to molecular structures:

```python
# Enable shadows by setting shadow strength (0 = off, higher = darker)
# Recommended range: 0.3-0.7
coot.set_shadow_strength(0.3)  # Subtle shadows (recommended)
coot.set_shadow_strength(0.5)  # Medium shadows
coot.set_shadow_strength(0.7)  # Strong shadows

# Shadow resolution (1-4, higher = sharper shadows)
# 4 is maximum quality
coot.set_shadow_resolution(4)  # Maximum resolution - sharpest shadows

# Shadow softness (1-3, higher = softer edges)
# 3 is maximum softness
coot.set_shadow_softness(3)  # Maximum softness - smoothest shadow edges

# Shadow box size (default: 66)
# Adjust if shadows are cut off
coot.set_shadow_box_size(66)
```

**Recommended shadow settings for publication:**
- **Strength**: 0.3 (subtle, doesn't overpower the structure)
- **Resolution**: 4 (maximum quality)
- **Softness**: 3 (smooth, professional appearance)

**When to use shadows:**
- Publication figures with complex 3D structures
- Presentations where depth perception is important
- Visualizations that benefit from enhanced spatial relationships
- Combined with mid-grey backgrounds for best effect

**Example:**
```python
def enable_publication_shadows():
    """Enable subtle, high-quality shadows for publication figures"""
    coot.set_shadow_strength(0.3)      # Subtle shadows
    coot.set_shadow_resolution(4)      # Maximum resolution
    coot.set_shadow_softness(3)        # Maximum softness
    return "Publication shadows enabled"

def disable_shadows():
    """Disable shadows"""
    coot.set_shadow_strength(0.0)
    return "Shadows disabled"
```

### Other Quality Settings

```python
# Anti-aliasing - smooths jagged edges
# Note: May need to restart Coot for this to take effect
coot.set_anti_aliasing(1)

# Enable fog for atmospheric depth
coot.set_use_fog(1)

# Perspective projection (more realistic depth)
coot.set_use_perspective_projection(1)
```

### Complete Publication Setup Example

```python
def setup_publication_graphics():
    """
    Configure Coot for creating publication-quality figures
    Optimized for journal covers and high-impact visualizations
    """
    
    # Background: 80% grey (excellent for print)
    coot.set_background_colour(0.8, 0.8, 0.8)
    
    # Enable outline mode for polished look
    coot.set_use_outline(1)
    
    # Enable all fancy graphics effects
    coot.set_use_ambient_occlusion(1)
    coot.set_use_fancy_lighting(1)
    coot.set_use_depth_blur(1)
    
    # Fine-tune SSAO for publication quality
    coot.set_ssao_strength(1.2)
    coot.set_ssao_radius(0.6)
    coot.set_ssao_kernel_n_samples(32)
    coot.set_ssao_blur_size(1)
    
    # Enable subtle, high-quality shadows
    coot.set_shadow_strength(0.3)
    coot.set_shadow_resolution(4)
    coot.set_shadow_softness(3)
    
    return "Publication graphics settings applied"

def setup_presentation_graphics():
    """
    Configure Coot for presentation slides (dark background)
    """
    
    # Black background for presentations
    coot.set_background_colour(0.0, 0.0, 0.0)
    
    # Enable outline mode
    coot.set_use_outline(1)
    
    # Enable fancy graphics
    coot.set_use_ambient_occlusion(1)
    coot.set_use_fancy_lighting(1)
    coot.set_use_depth_blur(1)
    
    # Enable shadows for presentations
    coot.set_shadow_strength(0.4)  # Slightly stronger for dark backgrounds
    coot.set_shadow_resolution(4)
    coot.set_shadow_softness(3)
    
    return "Presentation graphics settings applied"
```

## Summary

1. **Always use function wrappers** for complex code
2. **Set colors BEFORE creating representations**
3. **Reassign ALL colors** when recreating ribbons
4. **Use PDB header secondary structure** as the authoritative source
5. **Extract features to separate molecules** for emphasis
6. **Use consistent color schemes** for clarity
7. **Enable fancy graphics** for publication-quality figures
8. **Configure shadows** for enhanced depth (strength 0.3, resolution 4, softness 3)
9. **Choose appropriate background** for your publication medium (mid-grey works well with shadows)

Following these practices ensures reliable, publication-quality molecular graphics figures in Coot.

## Summary of Additions to coot-figure-making/SKILL.md

The following new content has been added to the figure-making skill based on today's session visualizing the 1ej6 reovirus core assembly.

---

## New Sections Added

### 1. Gaussian Surface Representations

**Key Topics Covered:**
- Creating Gaussian surfaces
- **CRITICAL RULE: Always use molecular symmetry coloring when applicable**
- Gaussian surface parameters (contour level, sigma, grid scale, box radius, B-factor)
- Working with generic display objects to selectively hide surfaces
- Downloading biological assemblies from PDBe

**Critical Best Practice:**
```python
# ALWAYS set molecular symmetry coloring for structures with symmetry
coot.set_gaussian_surface_chain_colour_mode(2)
coot.gaussian_surface(imol)
```

### 2. Downloading Biological Assemblies from PDBe

**Correct URL Format:**
- Use **hyphen** not underscore: `1ej6-assembly1` NOT `1ej6_assembly-1`
- File is **gzipped**: `.cif.gz` extension
- Base URL: `https://www.ebi.ac.uk/pdbe/static/entry/download/`

**Complete example code provided** for downloading, decompressing, and loading assemblies.

### 3. Advanced Graphics Settings for Gaussian Surfaces

**Publication-Quality SSAO Settings:**
- **256-512 kernel samples** for screenshots (not default 32)
- **Large radius (20-25)** for large assemblies
- Complete optimized settings for publication figures

**CRITICAL Discovery: Outline and Depth Blur Are Mutually Exclusive**
- Despite API allowing both, shader implementation makes them exclusive
- **Recommendation:** Use outline mode for Gaussian surfaces
- Clear documentation with correct/incorrect examples

### 4. Complete Publication Setup

**Optimized settings for Gaussian surface figures:**
- Dark grey background (0.2, 0.2, 0.2) for dramatic contrast
- Outline mode for polished cel-shaded look
- High-quality SSAO (256 samples, radius 25)
- Fancy lighting
- Complete working code example

---

## Updated Summary Section

Reorganized into three categories:

1. **Ribbon Representations** (existing best practices)
2. **Gaussian Surface Representations** (NEW - 6 best practices)
3. **General Publication Graphics** (existing, refined)

---

## Key Lessons Learned

### Image-Making Rules
1. **Always use molecular symmetry coloring** for structures with symmetry
2. **Use 256-512 SSAO kernel samples** for publication screenshots
3. **Outline and depth blur are mutually exclusive** in the shader
4. **Large SSAO radius (20-25)** needed for large assemblies
5. **Correct PDBe assembly URL format** with hyphen and .cif.gz

### Why This Matters
- Molecular symmetry coloring **immediately reveals biological organization**
- High kernel samples produce **smooth, professional shadows** without noise
- Large assemblies need different parameters than small molecules
- Biological assemblies show the **true functional form** of the structure

---

## Code Examples Added

1. Download biological assembly from PDBe
2. Hide first half of generic display objects
3. Setup publication-quality SSAO for Gaussian surfaces
4. Complete publication graphics configuration
5. All with proper error handling and documentation

---

## Session Context

These updates came from successfully:
1. Fetching PDB 1ej6 (reovirus core, 300-chain assembly)
2. Downloading biological assembly via PDBe API
3. Creating Gaussian surface with molecular symmetry coloring
4. Optimizing graphics settings for publication quality
5. Discovering shader limitations (outline/depth blur exclusivity)
6. Finding optimal SSAO parameters through experimentation

The result: A beautiful visualization of icosahedral viral symmetry with 5 colors showing the distribution of 5 unique chain types across 300 chains in the biological assembly.

