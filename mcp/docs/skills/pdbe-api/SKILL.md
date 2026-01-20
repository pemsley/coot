---
name: pdbe-api
description: "Query the PDBe (Protein Data Bank in Europe) REST API and Solr search API from within Coot to access structure metadata, validation data, revision history, and search capabilities"
---

# PDBe API Access from Coot

## Overview

The PDBe (Protein Data Bank in Europe) provides comprehensive REST and Solr-based APIs for programmatic access to structure data, validation reports, compound information, revision history, and search capabilities. Coot can access these APIs directly using the `coot_get_url_as_string_py()` function.

## Core Function

**`coot.coot_get_url_as_string_py(url)`** - Fetch URL content as string

```python
# Example: Get structure summary
result = coot.coot_get_url_as_string_py("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/4wa9")

# Parse JSON response
import json
data = json.loads(result)
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

### 2. Compound/Ligand Information

Get detailed information about a specific compound including formula, SMILES, InChI, and revision history:

```python
import json

comp_id = "AXI"  # 3-letter code
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{comp_id}"
result = coot.coot_get_url_as_string_py(url)
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

### 3. Validation Reports

Get residue-wise outliers including clashes, geometry outliers, and density fit issues:

```python
import json

pdb_id = "4wa9"
url = f"https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)

# Note: Response may have multiple JSON objects, parse carefully
# Access specific chain/residue validation data
```

**Available validation endpoints:**
- `/validation/residuewise_outlier_summary/entry/{pdb_id}` - Residue-level outliers
- `/validation/rama_sidechain_listing/entry/{pdb_id}` - Ramachandran and rotamer outliers
- `/validation/global_percentiles/entry/{pdb_id}` - Overall quality metrics

### 4. Structure Status and Revision History

Check if a structure has been superseded or revised:

```python
import json

pdb_id = "4wa9"
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/status/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)
data = json.loads(result)

status = data[pdb_id][0]
print(f"Status: {status['status_code']}")  # REL = released, OBS = obsolete
print(f"Since: {status['since']}")
print(f"Superseded by: {status['superceded_by']}")
print(f"Obsoletes: {status['obsoletes']}")
```

### 5. Ligand Binding Sites

Get information about ligand binding sites and interactions:

```python
import json

pdb_id = "4wa9"
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)
data = json.loads(result)

# Access ligand information per chain
for entity in data[pdb_id]:
    print(f"Chain: {entity['chain_id']}")
    for ligand in entity.get('ligands', []):
        print(f"  Ligand: {ligand['chem_comp_id']}")
        print(f"  Residue: {ligand['author_residue_number']}")
```

### 6. Assembly Information

Get biological assembly information:

```python
import json

pdb_id = "2hyy"
url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/{pdb_id}"
result = coot.coot_get_url_as_string_py(url)
data = json.loads(result)

for assembly in data[pdb_id]:
    print(f"Assembly {assembly['assembly_id']}: {assembly['name']}")
    print(f"  Form: {assembly['form']}")
    print(f"  Preferred: {assembly['preferred']}")
```

## Solr Search API

The Solr search API allows complex queries across the entire PDB.

### Basic Search Syntax

```python
import json
import urllib.parse

# Simple search
query = "molecule_name:Dihydrofolate AND organism_scientific_name:Human"
encoded_query = urllib.parse.quote(query)
url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={encoded_query}&wt=json&rows=10"

result = coot.coot_get_url_as_string_py(url)
data = json.loads(result)

# Access results
for doc in data['response']['docs']:
    print(f"PDB ID: {doc['pdb_id']}")
    print(f"Title: {doc.get('title', 'N/A')}")
```

### Common Solr Search Fields

- `pdb_id` - PDB entry ID
- `molecule_name` - Molecule name
- `organism_scientific_name` - Source organism
- `experimental_method` - Experimental method (e.g., "X-ray diffraction")
- `resolution` - Structure resolution
- `ligand_name` - Ligand/compound name
- `citation_title` - Publication title
- `deposition_date` - Deposition date
- `revision_date` - Revision date

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

### Finding Related Structures

Search for structures with the same ligand and protein:

```python
import json
import urllib.parse

def find_related_structures(protein_name, ligand_name):
    """Find structures containing specific protein-ligand combination"""
    query = f'molecule_name:"{protein_name}" AND ligand_name:"{ligand_name}"'
    encoded = urllib.parse.quote(query)
    url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q={encoded}&wt=json&rows=50"
    
    result = coot.coot_get_url_as_string_py(url)
    data = json.loads(result)
    
    print(f"Found {data['response']['numFound']} structures")
    for doc in data['response']['docs']:
        print(f"  {doc['pdb_id']}: {doc.get('title', 'N/A')}")
        print(f"    Resolution: {doc.get('resolution', 'N/A')} Å")

# Example usage
find_related_structures("ABL1", "imatinib")
```

## Error Handling

Always wrap API calls in try/except blocks:

```python
import json

def safe_pdbe_query(url):
    """Safely query PDBe API with error handling"""
    try:
        result = coot.coot_get_url_as_string_py(url)
        if not result or result == "":
            print(f"Empty response from {url}")
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
| Compound info | `/pdb/compound/summary/{comp_id}` |
| Validation outliers | `/validation/residuewise_outlier_summary/entry/{pdb_id}` |
| Structure status | `/pdb/entry/status/{pdb_id}` |
| Ligand binding sites | `/pdb/entry/ligand_monomers/{pdb_id}` |
| Search structures | `/search/pdb/select?q={query}` |

### Common Solr Query Patterns

| Query | Purpose |
|-------|---------|
| `pdb_id:4wa9` | Specific PDB entry |
| `molecule_name:kinase` | By protein name |
| `ligand_name:imatinib` | By ligand name |
| `resolution:[0 TO 2.0]` | High resolution structures |
| `revision_date:[20240101 TO 20241231]` | Recently revised |
| `experimental_method:"X-ray diffraction"` | By experimental method |

## Integration with Coot Workflows

### Example: Automated Structure Quality Check

```python
import json

def structure_quality_report(imol):
    """Generate quality report using PDBe API data"""
    
    # Get PDB ID from molecule
    pdb_file = coot.molecule_name(imol)
    # Extract PDB ID from filename (assumes format like "pdb4wa9.ent")
    import re
    match = re.search(r'pdb(\w{4})', pdb_file.lower())
    if not match:
        print("Could not extract PDB ID from filename")
        return
    
    pdb_id = match.group(1)
    
    # Get structure info
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}"
    result = coot.coot_get_url_as_string_py(url)
    data = json.loads(result)
    
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

## Summary

The PDBe API provides rich programmatic access to structure metadata, validation data, and search capabilities. Using `coot.coot_get_url_as_string_py()`, you can:

1. Check structure and ligand revision history
2. Access validation reports and quality metrics
3. Search for related structures
4. Get compound information and chemical details
5. Verify structure status and supersession
6. Integrate external data into Coot workflows

This enables powerful automated quality checks, structure comparison workflows, and data-driven validation within Coot.
