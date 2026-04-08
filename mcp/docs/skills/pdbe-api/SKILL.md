---
name: pdbe-api
description: "Query the PDBe (Protein Data Bank in Europe) REST API and Solr search API from within Coot to access structure metadata, validation data, revision history, search capabilities, and download coordinate files"
---

# PDBe API Access from Coot

## Overview

The PDBe (Protein Data Bank in Europe) provides comprehensive REST and Solr-based APIs for programmatic access to structure data, validation reports, compound information, revision history, and search capabilities. Coot can access these APIs directly using the `coot_get_url_as_string_py()` function, which now supports both text and binary data.

## Core Function

**`coot.coot_get_url_as_string_py(url)`** - Fetch URL content

Returns:
- Python `str` for text/JSON content (valid UTF-8)
- Python `bytes` for binary content (gzipped files, images, etc.)

```python
import json

# Example 1: Get structure summary (returns string)
result = coot.coot_get_url_as_string_py("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/4wa9")
data = json.loads(result)

# Example 2: Download gzipped coordinates (returns bytes)
import gzip
compressed = coot.coot_get_url_as_string_py("https://files.rcsb.org/download/4wa9.cif.gz")
decompressed = gzip.decompress(compressed)
imol = coot.read_coordinates_as_string(decompressed.decode('utf-8'), "4wa9")
```

## Main API Endpoints

### Entry-based API
Base URL: `https://www.ebi.ac.uk/pdbe/api/`

Documentation: https://www.ebi.ac.uk/pdbe/api/doc/

### Aggregated API
Base URL: `https://www.ebi.ac.uk/pdbe/graph-api/`

Documentation: https://pdbe.org/graph-api

### Search API (Solr)
Base URL: `https://www.ebi.ac.uk/pdbe/search/pdb/select?`

Documentation: https://www.ebi.ac.uk/pdbe/api/doc/search.html

## Common Query Patterns

### 1. Structure Summary and Metadata

Get basic information about a structure including deposition date, revision date, authors, and experimental method:

```python
import json

pdb_id = "4wa9"
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)

# Handle both string and bytes responses
if isinstance(result, bytes):
   result = result.decode('utf-8')

data = json.loads(result)

# Extract key information
entry = data[pdb_id][0]
print(f"Title: {entry['title']}")
print(f"Release date: {entry['release_date']}")
print(f"Revision date: {entry['revision_date']}")
print(f"Method: {entry['experimental_method']}")
print(f"Authors: {entry['entry_authors']}")
```

**Key fields in response:**
- `title` - Structure title
- `release_date` - Original deposition date (YYYYMMDD)
- `revision_date` - Most recent revision date (YYYYMMDD)
- `experimental_method` - List of experimental methods
- `entry_authors` - List of authors
- `number_of_entities` - Count of different entity types (protein, ligand, water, etc.)

### 2. Molecule/Entity Information

Get organism and molecule details:

```python
import json

pdb_id = "4wa9"
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

data = json.loads(result)

if pdb_id in data:
   for entity in data[pdb_id]:
      molecule_name = entity.get('molecule_name', ['N/A'])[0] if isinstance(entity.get('molecule_name', []), list) else entity.get('molecule_name', 'N/A')
      source = entity.get('source', [{}])[0] if isinstance(entity.get('source', []), list) else entity.get('source', {})
      organism = source.get('organism_scientific_name', 'N/A')
      expression_host = source.get('expression_host_scientific_name', 'N/A')
      
      print(f"Molecule: {molecule_name}")
      print(f"  Source organism: {organism}")
      if expression_host != 'N/A' and expression_host != organism:
         print(f"  Expression host: {expression_host}")
```

### 3. Compound/Ligand Information

Get detailed information about a specific compound including formula, SMILES, InChI, and revision history:

```python
import json

comp_id = "AXI"  # 3-letter code
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{comp_id}"
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

data = json.loads(result)

compound = data[comp_id][0]
print(f"Name: {compound['name']}")
print(f"Formula: {compound['formula']}")
print(f"Weight: {compound['weight']}")
print(f"Creation date: {compound['creation_date']}")
print(f"Revision date: {compound['revision_date']}")
print(f"InChI: {compound['inchi']}")
print(f"SMILES: {compound['smiles'][0]['name']}")
```

**Use case:** Check if a ligand definition was recently revised, which might explain geometry changes.

### 4. Downloading Coordinate Files

Download and load PDB/mmCIF files directly:

```python
import gzip

# Download current version from RCSB
pdb_id = "4wa9"
url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"

print(f"Downloading {pdb_id}...")
compressed = coot.coot_get_url_as_string_py(url)
print(f"Downloaded {len(compressed)} bytes (compressed)")

# Decompress
decompressed = gzip.decompress(compressed)
print(f"Decompressed to {len(decompressed)} bytes")

# Load into Coot
imol = coot.read_coordinates_as_string(decompressed.decode('utf-8'), f"{pdb_id}")
print(f"Loaded as molecule {imol}")
```

**Note:** The wwPDB versioned archive exists but is not currently accessible via HTTPS through this API. Use the current version from RCSB or PDBe.

### 5. Validation Reports

Get residue-wise outliers including clashes, geometry outliers, and density fit issues:

```python
import json

pdb_id = "4wa9"
url = f"https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

# Note: Response may have multiple JSON objects, parse carefully
# Access specific chain/residue validation data
```

**Available validation endpoints:**
- `/validation/residuewise_outlier_summary/entry/{pdb_id}` - Residue-level outliers
- `/validation/rama_sidechain_listing/entry/{pdb_id}` - Ramachandran and rotamer outliers
- `/validation/global_percentiles/entry/{pdb_id}` - Overall quality metrics

### 6. Structure Status and Revision History

Check if a structure has been superseded or revised:

```python
import json

pdb_id = "4wa9"
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/status/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

data = json.loads(result)

status = data[pdb_id][0]
print(f"Status: {status['status_code']}")  # REL = released, OBS = obsolete
print(f"Since: {status['since']}")
print(f"Superseded by: {status['superceded_by']}")
print(f"Obsoletes: {status['obsoletes']}")
```

### 7. Ligand Binding Sites

Get information about ligand binding sites and interactions:

```python
import json

pdb_id = "4wa9"
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

data = json.loads(result)

# Access ligand information per chain
for entity in data[pdb_id]:
    print(f"Chain: {entity['chain_id']}")
    for ligand in entity.get('ligands', []):
        print(f"  Ligand: {ligand['chem_comp_id']}")
        print(f"  Residue: {ligand['author_residue_number']}")
```

### 8. Assembly Information

Get biological assembly information:

```python
import json

pdb_id = "2hyy"
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

data = json.loads(result)

for assembly in data[pdb_id]:
    print(f"Assembly {assembly['assembly_id']}: {assembly['name']}")
    print(f"  Form: {assembly['form']}")
    print(f"  Preferred: {assembly['preferred']}")
```

## Solr Search API

The Solr search API allows complex queries across the entire PDB. However, it has important limitations.

### What Solr Search CAN Do (Well-Indexed Fields)

✅ **Metadata searches:**
- By release/deposition date: `release_year:2025`
- By experimental method: `experimental_method:"X-ray diffraction"`
- By resolution: `resolution:[* TO 2.0]` or `resolution:[1.5 TO 2.5]`
- By organism: `organism_scientific_name:"Homo sapiens"`

✅ **Presence/absence queries:**
- Has protein: `number_of_protein_chains:[1 TO *]`
- Has carbohydrate: `has_carb_polymer:Y`
- Has bound molecules: `has_bound_molecule:Y`
- Has modified residues: `has_modified_residues:Y`

✅ **Component searches:**
- Specific ligand: `chem_comp_id:ATP`
- Ligand name: `ligand_name:imatinib`
- Molecule name: `molecule_name:*kinase*`

✅ **Author/citation:**
- By author: `entry_authors:"Smith J"`
- By UniProt: `uniprot_accession:P12345`

✅ **Combined queries:**
```python
# Example: Human kinases with resolution < 2Å from 2024
query = 'release_year:2024 AND organism_scientific_name:"Homo sapiens" AND molecule_name:*kinase* AND resolution:[* TO 2.0]'
```

### What Solr Search CANNOT Do

❌ **Detailed connectivity:** Cannot search for "THR covalently bonded to NAG" or other specific atom-level connections

❌ **Geometry queries:** Cannot search for "bonds longer than X" or "angles outside range Y"

❌ **Spatial relationships:** Cannot search for "atoms within 5Å of ligand"

❌ **Sequence motifs:** Cannot search for "structures with GXGXXG motif"

❌ **Complex structural features:** Cannot search for "beta-barrel with 8 strands"

❌ **Validation specifics:** Cannot search for "residues with Ramachandran outliers at position X"

**The Pattern:** Solr indexes **metadata and simple categorical data**, not **structural details or relationships**.

For analyses requiring connectivity or geometry (like finding O-glycosylated threonines), you must:
1. Use Solr to find candidates (e.g., structures with NAG + resolution < 2.5Å)
2. Download those structures
3. Parse mmCIF connectivity tables locally
4. Extract geometric parameters

### Basic Search Syntax

```python
import json

# Simple search for high-resolution X-ray structures from 2024
query = "release_year:2024 AND experimental_method:\"X-ray diffraction\" AND resolution:[* TO 1.5]"
url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={query}&wt=json&rows=10&fl=pdb_id,title,resolution"

result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

data = json.loads(result)

print(f"Found {data['response']['numFound']} structures")
for doc in data['response']['docs']:
    print(f"  {doc['pdb_id']}: {doc.get('resolution', 'N/A')}Å")
    print(f"    {doc.get('title', 'N/A')[:70]}")
```

### Common Solr Search Fields

**Identifiers & Metadata:**
- `pdb_id` - PDB entry ID
- `molecule_name` - Molecule name
- `molecule_type` - Entity type: use `"Protein"` (capital P) — NOT `polypeptide(L)` or `protein`
- `molecule_sequence` - One-letter sequence string (stored but **not full-text indexed** — wildcard search like `*C*C*` returns 0 results; fetch and filter in Python instead)
- `polymer_length` - Length of the polymer entity in residues (supports range queries: `[5 TO 30]`)
- `number_of_polymer_residues` - Total residues across all chains in the entry
- `number_of_protein_chains` - Number of protein chains
- `organism_scientific_name` - Source organism
- `experimental_method` - Experimental method (e.g., `"X-ray diffraction"`)
- `resolution` - Structure resolution
- `ligand_name` - Ligand/compound name
- `citation_title` - Publication title
- `deposition_date` - Deposition date
- `revision_date` - Revision date

**Experimental Details:**
- `experimental_method` - Method (e.g., "X-ray diffraction", "Electron Microscopy", "Solution NMR")
- `resolution` - Structure resolution (numeric, use ranges like `[1.0 TO 2.0]`)
- `em_resolution` - EM-specific resolution
- `data_quality` - Overall quality metric

**Molecular Content:**
- `molecule_name` - Molecule name (supports wildcards: `*kinase*`)
- `molecule_type` - Type (Protein, DNA, RNA, etc.)
- `organism_scientific_name` - Source organism
- `organism_synonyms` - Alternative organism names
- `genus` - Organism genus
- `expression_host_scientific_name` - Expression system

**Ligands & Modifications:**
- `chem_comp_id` - Chemical component 3-letter code
- `ligand_name` - Ligand name
- `has_bound_molecule` - Y/N
- `has_carb_polymer` - Y/N (has carbohydrate)
- `has_modified_residues` - Y/N
- `number_of_bound_molecules` - Count

**Authors & Citations:**
- `entry_authors` - Entry authors
- `citation_authors` - Publication authors
- `citation_title` - Paper title
- `citation_year` - Publication year
- `pubmed_id` - PubMed ID

**Protein Details:**
- `uniprot_accession` - UniProt accession
- `uniprot_id` - UniProt ID
- `gene_name` - Gene name
- `go_id` - Gene Ontology ID

**Structure Properties:**
- `number_of_protein_chains` - Count
- `number_of_polymer_entities` - Count
- `assembly_composition` - Assembly type
- `symmetry_group` - Symmetry

### Practical Search Examples

**Example 1: High-resolution X-ray structures from 2024**
```python
import json

url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?q=release_year:2024%20AND%20experimental_method:\"X-ray%20diffraction\"%20AND%20resolution:[*%20TO%201.5]&wt=json&rows=5&fl=pdb_id,title,resolution"
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

Sorting by `sort=molecular_weight+asc` gives a 400 error. Use `polymer_length` instead
as a proxy for size, or use `number_of_polymer_residues` for entry-level size.

**5. Use `fq` (filter query) for range constraints**

Range filtering on numeric fields works well as a filter query:
```python
# Filter to entities with 5-30 residues:
url = "...&q=molecule_type:Protein&fq=polymer_length:[5+TO+30]&sort=polymer_length+asc..."
```

**Example 2: Human kinase structures**
```python
query = "organism_scientific_name:\"Homo sapiens\" AND molecule_name:*kinase*"
url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={query}&wt=json&rows=10&fl=pdb_id,title,resolution"
```

**Example 3: Cryo-EM structures better than 3Å from 2025**
```python
query = "release_year:2025 AND experimental_method:\"Electron Microscopy\" AND resolution:[* TO 3.0]"
url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={query}&wt=json&rows=10&fl=pdb_id,title,em_resolution"
```

**Example 4: Structures with carbohydrates**
```python
query = "has_carb_polymer:Y"
url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={query}&wt=json&rows=10&fl=pdb_id,title"
```

**Example 5: Structures of a specific protein from different species**
```python
# Find ABL1 structures from different mammals
query = "molecule_name:*ABL1* OR molecule_name:*ABL*kinase*"
url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={query}&wt=json&rows=100&fl=pdb_id,organism_scientific_name,title"
```

**6. Discover available fields by fetching a sample document with `fl=*`**

When you don't know what fields are in the index:
```python
url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?q=*:*&wt=json&rows=1&fl=*"
data = json.loads(coot.coot_get_url_as_string_py(url))
for k in sorted(data['response']['docs'][0].keys()):
    print(k)
```
Fields prefixed `q_` and `t_` are query/text variants of the base fields — ignore them
when exploring the schema.

**7. There is no `disulfide` or `bond_types` field in the Solr index**

To find structures with disulfide bonds, you must:
- Use Solr to find small proteins with ≥2 Cys in their sequence (fetch + filter in Python), then
- Use the PDBe REST API (`/pdb/entry/molecules/{pdb_id}`) to confirm the sequence and structure.

### Advanced Search Examples

**Find structures with specific ligand:**
```python
query = "ligand_name:axitinib"
```

**Find high-resolution kinase structures:**
```python
query = "molecule_name:kinase AND resolution:[0 TO 2.0]"
```

**Find structures revised in 2024:**
```python
query = "revision_date:[20240101 TO 20241231]"
```

## Practical Workflows

### Detecting Structure Revisions

Check if a structure has been significantly revised since release:

```python
import json
from datetime import datetime

def check_structure_revision(pdb_id):
    """Check if structure was revised and when"""
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}"
    result = coot.coot_get_url_as_string_py(url)
    
    if isinstance(result, bytes):
        result = result.decode('utf-8')
    
    data = json.loads(result)
    
    entry = data[pdb_id][0]
    release = entry['release_date']
    revision = entry['revision_date']
    
    # Convert to datetime for comparison
    release_dt = datetime.strptime(release, "%Y%m%d")
    revision_dt = datetime.strptime(revision, "%Y%m%d")
    
    days_diff = (revision_dt - release_dt).days
    years_diff = days_diff / 365.25
    
    print(f"PDB {pdb_id}:")
    print(f"  Released: {release}")
    print(f"  Revised: {revision}")
    print(f"  Time since release: {years_diff:.1f} years")
    
    if days_diff > 30:
        print(f"  WARNING: Structure revised {days_diff} days after release")
        return True
    return False

# Example usage
check_structure_revision("4wa9")
```

### Checking Ligand Revisions

Determine if a ligand definition was updated, which might explain geometry changes:

```python
import json

def check_ligand_revision(comp_id):
    """Check when a ligand was last revised"""
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{comp_id}"
    result = coot.coot_get_url_as_string_py(url)
    
    if isinstance(result, bytes):
        result = result.decode('utf-8')
    
    data = json.loads(result)
    
    compound = data[comp_id][0]
    print(f"Compound {comp_id} ({compound['name']}):")
    print(f"  Created: {compound['creation_date']}")
    print(f"  Revised: {compound['revision_date']}")
    
    if compound['creation_date'] != compound['revision_date']:
        print(f"  WARNING: Ligand definition was revised")
        return True
    return False

# Example usage
check_ligand_revision("AXI")
```

### Downloading and Comparing Structures

Download structures from different species and compare them:

```python
import gzip
import json

def download_and_load_structure(pdb_id):
    """Download and load a structure from RCSB"""
    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    
    print(f"Downloading {pdb_id}...")
    compressed = coot.coot_get_url_as_string_py(url)
    
    # Check if it's an error response (HTML)
    if isinstance(compressed, str) and compressed.startswith("<!DOCTYPE"):
        print(f"ERROR: Could not download {pdb_id}")
        return None
    
    decompressed = gzip.decompress(compressed)
    imol = coot.read_coordinates_as_string(decompressed.decode('utf-8'), pdb_id)
    print(f"Loaded as molecule {imol}")
    
    return imol

def compare_species_structures(pdb_id1, pdb_id2):
    """Download two structures and superpose them"""
    # Download both structures
    imol1 = download_and_load_structure(pdb_id1)
    imol2 = download_and_load_structure(pdb_id2)
    
    if imol1 is None or imol2 is None:
        print("Failed to download one or both structures")
        return
    
    # Superpose (using CA atoms from chain A, residues 240-400)
    print(f"\nSuperposing {pdb_id2} onto {pdb_id1}...")
    sel1 = "//A/240-400/CA"
    sel2 = "//A/240-400/CA"
    result = coot.superpose_with_atom_selection(imol1, imol2, sel1, sel2, 0)
    
    if result >= 0:
        print(f"Success! Structures superposed.")
    else:
        print("Superposition failed!")
    
    return imol1, imol2

# Example: Compare human and mouse ABL1
# First find structures using Solr
query = "molecule_name:*ABL1*"
url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={query}&wt=json&rows=100&fl=pdb_id,organism_scientific_name"
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
   result = result.decode('utf-8')

data = json.loads(result)

# Find human and mouse structures
human_pdbs = []
mouse_pdbs = []
for doc in data['response']['docs']:
   org = doc.get('organism_scientific_name', ['Unknown'])[0]
   if 'Homo sapiens' in org:
      human_pdbs.append(doc['pdb_id'])
   elif 'Mus musculus' in org:
      mouse_pdbs.append(doc['pdb_id'])

print(f"Human ABL1 structures: {len(human_pdbs)}")
print(f"Mouse ABL1 structures: {len(mouse_pdbs)}")

# Compare first human and mouse structures
if human_pdbs and mouse_pdbs:
   compare_species_structures(human_pdbs[0], mouse_pdbs[0])
```

### Finding Related Structures

Search for structures with the same ligand and protein:

```python
import json

def find_related_structures(protein_name, ligand_name=None):
    """Find structures containing specific protein-ligand combination"""
    if ligand_name:
        query = f'molecule_name:*{protein_name}* AND chem_comp_id:{ligand_name}'
    else:
        query = f'molecule_name:*{protein_name}*'
    
    url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={query}&wt=json&rows=50&fl=pdb_id,title,resolution,organism_scientific_name"
    
    result = coot.coot_get_url_as_string_py(url)
    
    if isinstance(result, bytes):
        result = result.decode('utf-8')
    
    data = json.loads(result)
    
    print(f"Found {data['response']['numFound']} structures")
    for doc in data['response']['docs']:
        org = doc.get('organism_scientific_name', ['N/A'])
        if isinstance(org, list):
            org = org[0] if org else 'N/A'
        
        print(f"  {doc['pdb_id']}: {doc.get('title', 'N/A')[:60]}")
        print(f"    Resolution: {doc.get('resolution', 'N/A')} Å")
        print(f"    Organism: {org}")

# Example usage
find_related_structures("ABL1", "STI")  # ABL1 with imatinib
```

## Error Handling

Always wrap API calls in try/except blocks and handle both string and bytes responses:

```python
import json

def safe_pdbe_query(url):
    """Safely query PDBe API with error handling"""
    try:
        result = coot.coot_get_url_as_string_py(url)
        if not result or result == "":
            print(f"Empty response from {url}")
            return None
        
        # Handle bytes response
        if isinstance(result, bytes):
            result = result.decode('utf-8')
        
        # Check for HTML error pages
        if result.startswith("<!DOCTYPE") or result.startswith("<html"):
            print(f"Received HTML error page instead of JSON")
            print(result[:200])
            return None
        
        data = json.loads(result)
        return data
    
    except json.JSONDecodeError as e:
        print(f"JSON parsing error: {e}")
        print(f"Response was: {result[:200]}...")
        return None
    
    except Exception as e:
        print(f"Error querying PDBe API: {e}")
        return None

# Example usage
data = safe_pdbe_query("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/4wa9")
if data:
    print("Success!")
```

## Common Issues and Solutions

### Issue: Binary vs Text Data

The function returns `bytes` for binary data (gzipped files) and `str` for text (JSON). Always check the type:

```python
result = coot.coot_get_url_as_string_py(url)

if isinstance(result, bytes):
    # Binary data - might be gzipped
    if result.startswith(b'\x1f\x8b'):  # gzip magic bytes
        import gzip
        decompressed = gzip.decompress(result)
        content = decompressed.decode('utf-8')
    else:
        content = result.decode('utf-8')
else:
    # Already a string
    content = result
```

### Issue: JSON parsing errors with validation endpoints

Some validation endpoints return multiple JSON objects or malformed responses. Handle carefully:

```python
# Instead of json.loads(), parse line by line or handle errors
try:
    data = json.loads(result)
except json.JSONDecodeError:
    # Try alternative parsing or just display raw result
    print("Could not parse JSON, raw response:")
    print(result[:1000])
```

### Issue: Unicode decode errors

If you get UnicodeDecodeError, the response might contain non-UTF-8 bytes. This should be handled automatically by the function now, but if you encounter issues:

```python
try:
    result = coot.coot_get_url_as_string_py(url)
except Exception as e:
    print(f"Error fetching URL: {e}")
```

### Issue: URL encoding for complex queries

Always encode special characters in Solr queries:

```python
import urllib.parse

query = "molecule_name:\"Protein kinase\" AND resolution:[0 TO 2.0]"
encoded = urllib.parse.quote(query)
url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={encoded}&wt=json"
```

### Issue: Rate limiting

The PDBe API may rate limit excessive requests. Add delays between batch queries:

```python
import time

pdb_ids = ["4wa9", "2hyy", "1iep"]
for pdb_id in pdb_ids:
    data = safe_pdbe_query(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}")
    # Process data...
    time.sleep(0.5)  # Wait 500ms between requests
```

## Quick Reference

### Most Useful Endpoints

| Purpose | Endpoint |
|---------|----------|
| Structure summary | `/pdb/entry/summary/{pdb_id}` |
| Molecule/organism info | `/pdb/entry/molecules/{pdb_id}` |
| Compound info | `/pdb/compound/summary/{comp_id}` |
| Validation outliers | `/validation/residuewise_outlier_summary/entry/{pdb_id}` |
| Structure status | `/pdb/entry/status/{pdb_id}` |
| Ligand binding sites | `/pdb/entry/ligand_monomers/{pdb_id}` |
| Search structures | `/search/pdb/select?q={query}` |
| Download coordinates | `https://files.rcsb.org/download/{pdb_id}.cif.gz` |

### Common Solr Query Patterns

| Query | Purpose |
|-------|---------|
| `pdb_id:4wa9` | Specific PDB entry |
| `molecule_name:*kinase*` | By protein name (wildcards) |
| `chem_comp_id:ATP` | Structures with specific ligand |
| `resolution:[0 TO 2.0]` | High resolution structures |
| `release_year:2025` | Structures from 2025 |
| `revision_year:2024` | Recently revised structures |
| `experimental_method:"X-ray diffraction"` | By experimental method |
| `organism_scientific_name:"Homo sapiens"` | By organism |
| `has_carb_polymer:Y` | Has carbohydrate |
| `has_bound_molecule:Y` | Has ligands |

### Combining Queries with AND/OR

```python
# Human kinases with resolution < 2Å from 2024
query = 'release_year:2024 AND organism_scientific_name:"Homo sapiens" AND molecule_name:*kinase* AND resolution:[* TO 2.0]'

# ABL1 from human OR mouse
query = 'molecule_name:*ABL1* AND (organism_scientific_name:"Homo sapiens" OR organism_scientific_name:"Mus musculus")'
```

## Integration with Coot Workflows

### Example: Automated Structure Quality Check

```python
import json

def structure_quality_report(imol):
    """Generate quality report using PDBe API data"""
    
    # Get PDB ID from molecule
    pdb_file = coot.molecule_name(imol)
    # Extract PDB ID from filename (assumes format like "pdb4wa9.ent" or "4wa9")
    import re
    match = re.search(r'(\d\w{3})', pdb_file.lower())
    if not match:
        print("Could not extract PDB ID from filename")
        return
    
    pdb_id = match.group(1)
    
    # Get structure info
    data = safe_pdbe_query(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}")
    if not data:
        return
    
    entry = data[pdb_id][0]
    
    print("=" * 60)
    print(f"STRUCTURE QUALITY REPORT: {pdb_id.upper()}")
    print("=" * 60)
    print(f"Title: {entry['title']}")
    print(f"Method: {entry['experimental_method']}")
    print(f"Released: {entry['release_date']}")
    print(f"Revised: {entry['revision_date']}")
    
    # Check for significant revisions
    if entry['revision_date'] != entry['release_date']:
        from datetime import datetime
        release = datetime.strptime(entry['release_date'], "%Y%m%d")
        revision = datetime.strptime(entry['revision_date'], "%Y%m%d")
        days = (revision - release).days
        print(f"\n⚠️  STRUCTURE REVISED {days} days after release")
        print("   Check PDBe for revision details")
    
    print("=" * 60)

# Usage: structure_quality_report(0)
```

## Resources

- **PDBe API Documentation**: https://www.ebi.ac.uk/pdbe/api/doc/
- **Aggregated API**: https://pdbe.org/graph-api
- **Search API**: https://www.ebi.ac.uk/pdbe/api/doc/search.html
- **Mailing List**: pdbe-api-users@ebi.ac.uk
- **GitHub Examples**: https://github.com/PDBeurope/pdbe-api-training
- **RCSB Downloads**: https://files.rcsb.org/download/

## Summary

The PDBe API provides rich programmatic access to structure metadata, validation data, and search capabilities. Using `coot.coot_get_url_as_string_py()`, you can:

1. **Download coordinate files** - Get structures in mmCIF/PDB format (gzipped)
2. **Check revision history** - Detect structures and ligands that have been revised
3. **Access validation reports** - Get quality metrics and outlier information
4. **Search across the PDB** - Find related structures, compare organisms, filter by properties
5. **Get compound information** - Access chemical details, SMILES, InChI
6. **Verify structure status** - Check for supersession or obsolescence
7. **Integrate external data** - Bring PDB metadata into Coot workflows

**Key Capabilities:**
- Binary data support (download gzipped files)
- Comprehensive metadata access
- Powerful search with well-understood limitations
- Cross-species structure comparison
- Revision tracking and provenance checking

**Key Limitations:**
- Solr search cannot query detailed connectivity or geometry
- Versioned coordinates not accessible via HTTPS (use current versions)
- For analyses requiring atom-level connectivity, download and parse structures locally

This enables powerful automated quality checks, cross-species structure comparison, data-driven validation, and integration of PDB metadata into Coot-based structural biology workflows.
