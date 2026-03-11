---
name: coot-validation
description: "Comprehensive structure validation combining model-to-map analysis and unmodeled density detection"
---
# Coot Structure Validation Best Practices

## Overview

When performing structure validation in Coot with both a model and a map, you need to analyze the structure from three complementary perspectives:

1. **Model-to-Map Validation**: How well does the existing model fit the density?
2. **Map-to-Model Validation**: Where is there significant density that is NOT explained by the model?
3. **Atom Overlap Validation**: Are there steric clashes between atoms in the model?

All three perspectives are essential for comprehensive validation.

## The Three Types of Validation

### Model-to-Map: Finding Problems in Your Model

These functions analyze how well your current model fits the density. They lead you to **places in the model** that need attention:

- Poor density correlation
- Ramachandran outliers
- Rotamer outliers
- Geometry violations

### Understanding Rotamer Validation: Two Different Scores

**CRITICAL: `rotamer_graphs_py()` and `score_rotamers_py()` report different things**

There are two functions that report rotamer information, and they measure fundamentally different aspects:

#### 1. `rotamer_graphs_py(imol)` - Continuous Probability Density
Returns the **probability density** at the exact chi angles of the current conformation.
- Measures: "How likely is THIS specific chi1, chi2, chi3... combination?"
- Scale: 0-100%, where 100% = peak of the probability distribution
- Use for: **Primary validation metric** - this is what you should check

```python
rotamers = coot.rotamer_graphs_py(0)
# Returns: [[chain_id, resno, ins_code, score_percentage, resname], ...]
# score_percentage is the continuous probability density at the actual chi angles
```

**Interpretation guidelines:**
- **> 50%**: Excellent - in high-density region
- **20-50%**: Good - acceptable conformation
- **5-20%**: Marginal - check density fit carefully
- **< 5%**: Poor - likely wrong (but check density!)
- **< 1%**: Very poor - almost certainly wrong

#### 2. `score_rotamers_py(...)` - Discrete Bin Probabilities
Returns the **discrete rotamer library** showing what % of structures have each named rotamer.
- Measures: "How common is rotamer 'm-85' vs 't80' vs 'p90' across all proteins?"
- Scale: Probabilities sum to ~100% across all discrete bins
- Use for: **Understanding alternatives** - what other conformations exist?

```python
rotamers = coot.score_rotamers_py(0, "A", 42, "", "", 1, 1, 0.001)
# Returns: [[name, probability, density_score, atom_list, richardson_score], ...]
# probability is the discrete bin frequency (e.g., m-85 appears in 43% of structures)
```

**Why the scores differ:**
- Discrete bins: "43% of TYR have m-85 rotamer" (which bin?)
- Continuous density: "83% probability at chi1=-62.3°, chi2=-85.1°" (exact angles)
- The continuous score can exceed discrete bin probabilities!

#### Rotamer Complexity and Expected Scores

**Number of chi angles determines score expectations:**

| Residue Type | Chi Angles | Typical # Rotamers | Best Rotamer % | Good Score |
|--------------|------------|-------------------|----------------|------------|
| VAL, THR, SER | 1 | 3 | 70-75% | > 40% |
| PHE, TYR, ASP, ASN | 2 | 4-9 | 35-45% | > 20% |
| GLU, GLN, MET, ILE, LEU | 3 | 9-15 | 20-35% | > 10% |
| LYS, ARG | 4 | 30-35 | 9-10% | > 5% |

**Key insight:** More chi angles = more rotamers = lower individual probabilities

**Examples from validation:**
- **LEU (3 chi)**: Best possible = 59%, so 30% is good, 5% is poor
- **ARG (4 chi)**: Best possible = 9%, so 5% is good, 0.5% is poor
- **VAL (1 chi)**: Best possible = 73%, so 40% is good, 5% is terrible

#### Proper Rotamer Validation Workflow

**DON'T:** Use arbitrary cutoffs like "< 10% is bad"

**DO:** Use context-aware validation:

```python
# Step 1: Get current rotamer scores
rotamers = coot.rotamer_graphs_py(0)

# Step 2: For flagged residues, get alternatives
for chain, resno, inscode, score, resname in rotamers:
    if score < 20:  # Preliminary flag
        # Get all possible rotamers to understand context
        alternatives = coot.score_rotamers_py(0, chain, resno, "", "", 1, 1, 0.001)
        
        if len(alternatives) == 0:
            continue  # GLY, ALA - no rotamers
        
        # Sort by density fit
        sorted_alts = sorted(alternatives, key=lambda x: x[2], reverse=True)
        best_density = sorted_alts[0][2]
        current_density = sorted_alts[0][2]  # Approximate
        
        # Check how many rotamers exist
        n_rotamers = len(alternatives)
        
        # Decision logic
        if n_rotamers < 5 and score < 10:
            # Few rotamers (VAL, PHE, etc.) and low score = likely wrong
            print(f"PROBLEM: {chain}/{resno} {resname}: {score:.1f}% (few alternatives)")
        elif n_rotamers > 20 and score < 2:
            # Many rotamers (LYS, ARG) and very low score = likely wrong
            print(f"PROBLEM: {chain}/{resno} {resname}: {score:.1f}% (many alternatives)")
        elif best_density - current_density > 3.0:
            # Alternative has much better density fit
            print(f"PROBLEM: {chain}/{resno} {resname}: better rotamer available")
```

#### Combined Validation: Rotamers + Density + Clashes

**A residue needs fixing if:**
1. **Rotamer score is low FOR THAT RESIDUE TYPE** (bottom 10% of possibilities), AND
2. **Density correlation is poor** (< 0.7), AND/OR
3. **Causes steric clashes** (> 2.0 Å³ overlap)

**A low rotamer score alone is NOT sufficient** - always check:
- Is this low for this residue type? (compare to alternatives)
- Does the density support this conformation?
- Are there clashes that would be resolved by changing rotamers?

### Atom Overlaps: Finding Steric Clashes

Atom overlap detection identifies clashes between atoms that may not be caught by local geometry validation. These reveal **packing problems** such as:

- Clashes between distant residues
- Side-chain/side-chain clashes
- Backbone/side-chain clashes
- Clashes with symmetry mates

**Critical insight**: Ramachandran and rotamer validation catch local geometry problems (within a residue or its immediate neighbors), while atom overlap detection catches global packing problems between any atoms in the structure.

### Map-to-Model: Finding Missing Features

The blob-finding function identifies regions of significant density that are not explained by your current model. It leads you to **places in the map** where you might be missing:

- Waters
- Ligands
- Alternative conformations
- Metal ions
- Other small molecules
- Missing residues or loops

## Critical Function: find_blobs_py()

**Always include blob detection when performing structure validation with a map.**

```python
blobs = coot.find_blobs_py(
    imol_model=0,              # your protein model
    imol_map=1,                # the map to search (often difference map)
    cut_off_density_level=3.0  # sigma threshold (typically 2.5-4.0)
)

# Returns: list of (position, score) tuples
# [(clipper::Coord_orth, float), ...]
```

### Parameters

- `imol_model`: The model molecule - density explained by this model will be excluded
- `imol_map`: The map to search for blobs (usually a difference map, but can be regular map)
- `cut_off_density_level`: Sigma threshold for blob detection
  - **3.0 sigma**: Standard threshold for significant features
  - **2.5 sigma**: More sensitive, finds weaker features
  - **4.0 sigma**: Conservative, only strong features

### Understanding the Results

```python
for position, score in blobs:
    x = position.x()
    y = position.y()
    z = position.z()
    print(f"Blob at ({x:.2f}, {y:.2f}, {z:.2f}) - score: {score:.2f}")
```

The **score** represents the strength/volume of the unmodeled density. Higher scores indicate more significant features that should be investigated.

### Recentering the View

The user likes to see what you are considering and how you change the model, so, if you
can, try to use coot.set_rotation_centre() or coot.set_go_to_atom_chain_residue_atom_name()
or some such to bring the currently interesting issue to the centre of the screen.

## Complete Validation Workflow

### 1. Model-to-Map Validation

```python
# Ramachandran outliers
rama_outliers = coot.all_molecule_ramachandran_score_py(0)

# Rotamer outliers  
rotamer_outliers = coot.rotamer_graphs_py(0)

# Per-residue density correlation
correlation_stats = coot.map_to_model_correlation_stats_per_residue_range_py(
    0,      # imol_model
    "A",    # chain_id
    1,      # start_resno
    100,    # end_resno
    1       # imol_map
)

# Geometry validation
chiral = coot.chiral_volume_errors_py(0)
```

### 2. Atom Overlap Validation

```python
# Get worst 30 atom overlaps
overlaps = coot.molecule_atom_overlaps_py(0, 30)

# Check for severe clashes
severe_clashes = [o for o in overlaps if o['overlap-volume'] > 5.0]
if severe_clashes:
    print(f"WARNING: {len(severe_clashes)} severe clashes found!")
    
# For full analysis (caution: can be very large!)
# all_overlaps = coot.molecule_atom_overlaps_py(0, -1)
```

### 3. Map-to-Model Validation (Blobs)

```python
# Find unmodeled density in difference map
diff_map_blobs = coot.find_blobs_py(
    imol_model=0,
    imol_map=2,  # difference map
    cut_off_density_level=3.0
)

# Find features in regular map (alternative approach)
regular_map_blobs = coot.find_blobs_py(
    imol_model=0,
    imol_map=1,  # 2mFo-DFc map
    cut_off_density_level=1.0  # Lower threshold for fitted map
)
```

### 3. Comprehensive Validation Report

```python
def comprehensive_validation(imol_model, imol_map, imol_diff_map=None):
    """
    Perform complete structure validation combining model and map analysis.
    
    Returns dictionary with all validation metrics.
    """
    results = {}
    
    # Model-to-map validation
    results['ramachandran'] = coot.all_molecule_ramachandran_score_py(imol_model)
    results['rotamers'] = coot.rotamer_graphs_py(imol_model)
    
    # Atom overlap validation
    results['atom_overlaps'] = coot.molecule_atom_overlaps_py(imol_model, 30)
    severe_clashes = [o for o in results['atom_overlaps'] if o['overlap-volume'] > 5.0]
    results['severe_clash_count'] = len(severe_clashes)
    
    # Per-residue correlation (requires chain info)
    import coot_utils
    chains = coot_utils.chain_ids(imol_model)
    results['correlation_by_chain'] = {}
    
    for chain in chains:
        n_residues = coot.chain_n_residues(chain, imol_model)
        if n_residues > 0:
            stats = coot.map_to_model_correlation_stats_per_residue_range_py(
                imol_model, chain, 1, 9999, imol_map
            )
            results['correlation_by_chain'][chain] = stats
    
    # Map-to-model validation (blobs)
    if imol_diff_map is not None:
        results['diff_map_blobs'] = coot.find_blobs_py(
            imol_model, imol_diff_map, 3.0
        )
    
    results['map_blobs'] = coot.find_blobs_py(
        imol_model, imol_map, 1.0
    )
    
    return results

# Usage
validation = comprehensive_validation(
    imol_model=0,
    imol_map=1,
    imol_diff_map=2
)
```

## Interpreting Blob Results

### What Different Maps Tell You

**Difference Map (mFo-DFc) Blobs:**
- **Positive blobs (>3σ)**: Missing atoms/features - something should be added here
- **Negative blobs (<-3σ)**: Incorrectly modeled atoms - something should be removed/moved
- Most reliable for finding genuine missing features

**Regular Map (2mFo-DFc) Blobs:**
- Less sensitive to model bias
- Good for finding larger missing features (domains, ligands)
- Use lower sigma threshold (0.5-1.5σ)

### Common Blob Interpretations

```python
blobs = coot.find_blobs_py(0, 2, 3.0)  # diff map, 3 sigma

# Large score (>50): Likely missing ligand, metal, or several waters
# Medium score (10-50): Likely 1-3 waters or alternative conformation  
# Small score (3-10): Likely single water or weak alternative conformation

for position, score in blobs:
    if score > 50:
        print(f"Large feature at {position} - investigate for ligand/metal")
    elif score > 10:
        print(f"Medium feature at {position} - likely waters")
    else:
        print(f"Small feature at {position} - check carefully")
```

## Critical Function: molecule_atom_overlaps_py()

**Always include atom overlap checking when validating structure geometry.**

```python
# Get worst 30 atom overlaps (default behavior after API update)
overlaps = coot.molecule_atom_overlaps_py(
    imol=0,
    n_pairs=30  # Number of worst overlaps to return (default: 30)
)

# Get ALL overlaps (use with caution - can be hundreds!)
all_overlaps = coot.molecule_atom_overlaps_py(
    imol=0,
    n_pairs=-1  # -1 means return all overlaps
)

# Each overlap is a dict with:
# {
#     'atom-1-spec': [imol, chain, resno, inscode, atom_name, altconf],
#     'atom-2-spec': [imol, chain, resno, inscode, atom_name, altconf],
#     'overlap-volume': float,  # in Ų
#     'radius-1': float,
#     'radius-2': float
# }
```

### Understanding Overlap Results

**Overlap volume** indicates severity:
- **>5.0 Å³**: Severe clash - atoms are deeply interpenetrating
- **2.0-5.0 Å³**: Moderate clash - needs immediate attention
- **0.5-2.0 Å³**: Minor clash - may be acceptable in some contexts
- **<0.5 Å³**: Very minor overlap - often acceptable

**Common clash patterns:**

```python
overlaps = coot.molecule_atom_overlaps_py(0, 30)

for overlap in overlaps:
    atom1 = overlap['atom-1-spec']
    atom2 = overlap['atom-2-spec']
    volume = overlap['overlap-volume']
    
    chain1, res1, atom_name1 = atom1[1], atom1[2], atom1[4]
    chain2, res2, atom_name2 = atom2[1], atom2[2], atom2[4]
    
    if volume > 5.0:
        print(f"SEVERE: {chain1}/{res1} {atom_name1} ↔ {chain2}/{res2} {atom_name2}: {volume:.2f} Ų")
    elif volume > 2.0:
        print(f"MODERATE: {chain1}/{res1} {atom_name1} ↔ {chain2}/{res2} {atom_name2}: {volume:.2f} Ų")
```

### Why Overlaps Are Essential

**Example from tutorial data:**
- Ramachandran validation found outliers at A/41-42
- Overlap validation revealed **A/41 O ↔ A/43 N: 2.07 Å³** backbone clash
- BUT also found **A/2 ↔ A/89 clashes (7.45, 6.40 Å³)** between distant residues that had PERFECT local geometry!

**Key lesson**: A model can have perfect Ramachandran and rotamer scores but catastrophic packing problems. You need
both local geometry validation (Rama/rotamer) AND global packing validation (overlaps).

## Prioritizing Validation Fixes

### Understanding Rotamer Scores in Context

**CRITICAL: Never use absolute rotamer score thresholds without considering residue type**

Before prioritizing rotamer fixes, understand what's "bad" for each residue:

```python
def assess_rotamer_severity(chain, resno, score, resname):
    """
    Determine if a rotamer score is actually problematic.
    Returns: 'critical', 'moderate', 'acceptable', or 'good'
    """
    # Get all possible rotamers to understand the distribution
    alternatives = coot.score_rotamers_py(0, chain, resno, "", "", 1, 1, 0.001)
    n_rotamers = len(alternatives)
    
    # Context-aware thresholds based on number of possible rotamers
    if n_rotamers <= 3:  # VAL, THR, SER (1 chi)
        if score < 10: return 'critical'
        elif score < 30: return 'moderate'
        else: return 'acceptable'
    elif n_rotamers <= 9:  # PHE, TYR, etc. (2 chi)
        if score < 5: return 'critical'
        elif score < 15: return 'moderate'
        else: return 'acceptable'
    elif n_rotamers <= 15:  # GLU, GLN, MET (3 chi)
        if score < 3: return 'critical'
        elif score < 10: return 'moderate'
        else: return 'acceptable'
    else:  # LYS, ARG (4 chi, 30+ rotamers)
        if score < 1: return 'critical'
        elif score < 5: return 'moderate'
        else: return 'acceptable'
```

### 1. Address High-Confidence Issues First

**Priority 1: Combined problems (multiple red flags)**
1. **Severe atom overlaps (>5 Å³)** - atoms deeply interpenetrating
2. **Poor rotamer + poor density + clashes** - triple failure
   - Example: Score < 5% for 2-chi residue, correlation < 0.5, clashes > 2 Å³
3. **Ramachandran outliers with poor density correlation** - backbone is wrong
4. **Large difference map blobs (>4σ)** - definitely missing features

**Priority 2: Single severe issues**
1. **Context-aware rotamer outliers with poor density**:
   - VAL/THR/SER < 10% AND correlation < 0.7
   - PHE/TYR/ASP < 5% AND correlation < 0.7
   - GLU/GLN/MET < 3% AND correlation < 0.7
   - LYS/ARG < 1% AND correlation < 0.7
2. **Moderate atom overlaps (2-5 Å³) between distant residues**

**Important:** A low rotamer score with GOOD density correlation (>0.8) may be correct - it could be a genuine unusual but real conformation. Don't "fix" it unless there's supporting evidence (clashes, poor density, chemical implausibility).

### 2. Investigate Moderate Issues

1. **Medium difference map blobs (3-4σ)** - probably real features
2. **Context-appropriate moderate rotamer scores with marginal density**:
   - Check if alternative rotamer has much better density fit
   - Compare alternatives with `score_rotamers_py()`
3. **Minor atom overlaps (0.5-2 Å³)** - may need adjustment
4. **Moderate geometry outliers** - may need refinement

### 3. Review Low-Priority Items

1. **Small blobs near model** - might be noise or minor adjustments
2. **Very minor overlaps (<0.5 Å³)** - often acceptable
3. **Isolated geometry outliers with good density** - may be genuine
4. **Borderline Ramachandran outliers** - check context

## Automated Validation Example

```python
def validate_and_fix_chain(imol_model, chain_id, imol_map, imol_diff_map):
    """
    Automated validation and suggested fixes for a chain.
    """
    issues = []
    
    # 1. Check for atom overlaps
    overlaps = coot.molecule_atom_overlaps_py(imol_model, 50)
    for overlap in overlaps:
        atom1 = overlap['atom-1-spec']
        atom2 = overlap['atom-2-spec']
        volume = overlap['overlap-volume']
        
        # Only report if at least one atom is in this chain
        if atom1[1] == chain_id or atom2[1] == chain_id:
            severity = 'high' if volume > 5.0 else ('medium' if volume > 2.0 else 'low')
            issues.append({
                'type': 'atom_overlap',
                'atom1': f"{atom1[1]}/{atom1[2]} {atom1[4]}",
                'atom2': f"{atom2[1]}/{atom2[2]} {atom2[4]}",
                'severity': severity,
                'value': volume
            })
    
    # 2. Check correlation for each residue
    stats = coot.map_to_model_correlation_stats_per_residue_range_py(
        imol_model, chain_id, 1, 9999, imol_map
    )
    
    for residue_spec, correlation in stats:
        if correlation < 0.7:  # Poor fit threshold
            issues.append({
                'type': 'poor_correlation',
                'residue': residue_spec,
                'severity': 'high',
                'value': correlation
            })
    
    # 3. Find nearby blobs that might explain poor correlation
    blobs = coot.find_blobs_py(imol_model, imol_diff_map, 3.0)
    
    for position, score in blobs:
        issues.append({
            'type': 'unmodeled_density',
            'position': (position.x(), position.y(), position.z()),
            'severity': 'high' if score > 50 else 'medium',
            'score': score
        })
    
    # 4. Check Ramachandran
    rama = coot.all_molecule_ramachandran_score_py(imol_model)
    for outlier in rama:
        if outlier[4] == 'OUTLIER':  # Ramachandran region
            issues.append({
                'type': 'ramachandran_outlier',
                'residue': outlier[0:3],  # chain, resno, inscode
                'severity': 'high'
            })
    
    return sorted(issues, key=lambda x: {'high': 0, 'medium': 1, 'low': 2}[x['severity']])

# Usage
issues = validate_and_fix_chain(0, "A", 1, 2)
for issue in issues[:10]:  # Top 10 issues
    print(f"{issue['type']}: {issue}")
```

## Common Patterns

### Water Placement from Blobs

```python
# Find blobs in difference map
blobs = coot.find_blobs_py(0, 2, 3.0)

# Add waters at blob positions
for position, score in blobs:
    if 5 < score < 30:  # Typical water blob size
        # Check if appropriate for water
        x, y, z = position.x(), position.y(), position.z()
        # Add water at this position
        coot.place_typed_atom_at_pointer("HOH")
```

### Missing Residue Detection

```python
# Look for large blobs that might be missing residues
blobs = coot.find_blobs_py(0, 2, 3.0)

missing_residue_candidates = [
    (pos, score) for pos, score in blobs 
    if score > 100  # Large feature
]

for position, score in missing_residue_candidates:
    print(f"Large unmodeled density at {position} - check for missing residues")
```

## Key Takeaways

1. **Always check atom overlaps** - local geometry can be perfect while global packing is catastrophic
2. **Always run blob detection** when you have both model and map
3. **Use difference maps** (mFo-DFc) for most sensitive blob detection
4. **Combine all three validation types** (model-to-map, overlaps, map-to-model) for complete picture
5. **Context-aware rotamer validation** - never use absolute thresholds:
   - Understand the difference between `rotamer_graphs_py()` (continuous density) and `score_rotamers_py()` (discrete bins)
   - Use `rotamer_graphs_py()` scores for validation (measures actual chi angles)
   - Consider number of possible rotamers: 5% is terrible for VAL but good for ARG
   - Always check density correlation and available alternatives before "fixing"
6. **Prioritize by severity AND confidence** - fix issues with multiple red flags first
7. **Iterate** - fixing one issue may reveal others
8. **Document** - keep track of what you fixed and why

## Function Reference

### Essential Validation Functions

```python
# Rotamer validation - primary metric (continuous probability density)
rotamers = coot.rotamer_graphs_py(imol)
# Returns: [[chain_id, resno, ins_code, score_percentage, resname], ...]

# Rotamer alternatives - for understanding context
alternatives = coot.score_rotamers_py(imol, chain, resno, "", "", imol_map, 1, 0.001)
# Returns: [[name, probability, density_score, atom_list, richardson_score], ...]

# Atom overlap detection
overlaps = coot.molecule_atom_overlaps_py(imol, n_pairs=30)  # Default: 30 worst
all_overlaps = coot.molecule_atom_overlaps_py(imol, n_pairs=-1)  # All overlaps

# Blob detection (map-to-model)
blobs = coot.find_blobs_py(imol_model, imol_map, sigma_cutoff)

# Ramachandran validation
rama = coot.all_molecule_ramachandran_score_py(imol)

# Density correlation (model-to-map)
corr = coot.map_to_model_correlation_stats_per_residue_range_py(
    imol, chain, imol_map, n_per_range, exclude_NOC_flag
)

# Geometry validation
chiral = coot.chiral_volume_errors_py(imol)
```

Remember: 
- **Model-to-map tells you what's wrong with your model**
- **Atom overlaps tell you about packing problems** 
- **Map-to-model tells you what you're missing**
- **Rotamer scores must be interpreted in context of residue type**

