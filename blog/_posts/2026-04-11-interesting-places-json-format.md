---
layout: post
title: "Coot Interesting Places JSON Format"
date: Sat 11 Apr 2026 06:41:20 BST
---

# Coot Interesting Places JSON Format

Coot can display a navigable dialog of "interesting places" in a structure. For example,
one might want to see the positions of validation outliers such as clashes, or indeed any
noteworthy positions you might want to flag.
The data are provided as a JSON file and loaded with:

```python
coot.read_interesting_places_json_file("my_outliers.json")
```

The dialog is applied to whichever molecule is currently active.

## Top-Level Structure

```json
{
    "title": "My Validation Results",
    "sections": [ ... ]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `title` | string | No | Title shown in the dialog window. Defaults to `"<Title>"` if omitted. |
| `sections` | array | Yes | Array of section objects (see below). |

## Sections

Each section becomes a labelled group in the dialog. A section contains a title and an array of items:

```json
{
    "title": "Van der Waals Clashes",
    "items": [ ... ]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `title` | string | No | Heading for this group of items. |
| `items` | array | Yes | Array of item objects. |

## Items

Each item represents a single interesting place. When the user clicks it, Coot navigates to the corresponding position in the structure.

Every item must have a `position-type` field that determines how Coot
locates the place. The supported position types are described below.

### Common Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `position-type` | string | Yes | One of: `"by-atom-spec"`, `"by-atom-spec-pair"`, `"by-residue-spec"`, `"by-coordinates"`. |
| `label` | string | Yes | Text shown on the button/entry in the dialog. |
| `badness` | number | No | A numeric severity score (between 0 and 1, where 1 is fully bad). When present, items are colour-coded by badness in the dialog. |

### Position Type: `by-atom-spec`

Navigate to a single atom.

```json
{
    "position-type": "by-atom-spec",
    "label": "Rotamer outlier: A PHE 42",
    "atom-spec": ["A", 42, "", "CA", ""],
    "badness": 3.5
}
```

An atom-spec is a 5-element array:

| Index | Type | Description |
|-------|------|-------------|
| 0 | string | Chain ID |
| 1 | integer | Residue number |
| 2 | string | Insertion code (empty string if none) |
| 3 | string | Atom name |
| 4 | string | Alternate conformation ID (empty string if none) |

### Position Type: `by-residue-spec`

Navigate to a residue (centres on the first atom found).

```json
{
    "position-type": "by-residue-spec",
    "label": "Missing sidechain: A ARG 105",
    "residue-spec": ["A", 105, ""],
    "badness": 2.0
}
```

A residue-spec is a 3-element array:

| Index | Type | Description |
|-------|------|-------------|
| 0 | string | Chain ID |
| 1 | integer | Residue number |
| 2 | string | Insertion code (empty string if none) |

Note: a 5-element array (matching atom-spec format) is also accepted --- the
last two elements are ignored.

### Position Type: `by-atom-spec-pair`

Navigate to a pair (or more) of atoms. Used for clashes, bond/angle/torsion
outliers, or any item that involves an interaction between atoms.

```json
{
    "position-type": "by-atom-spec-pair",
    "label": "Clash: A GLN 323 O - A ILE 324 HG12, z=-6.2",
    "atom-1-spec": ["A", 323, "", "O", ""],
    "atom-2-spec": ["A", 324, "", "HG12", ""]
}
```

For angle outliers (3 atoms) and torsion outliers (4 atoms), additional
atom specs can be included:

```json
{
    "position-type": "by-atom-spec-pair",
    "label": "Torsion outlier: omega, A ALA 55 - A LYS 56, z=-13.9",
    "atom-1-spec": ["A", 55, "", "CA", ""],
    "atom-2-spec": ["A", 55, "", "C", ""],
    "atom-3-spec": ["A", 56, "", "N", ""],
    "atom-4-spec": ["A", 56, "", "CA", ""]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `atom-1-spec` | atom-spec | Yes | First atom. |
| `atom-2-spec` | atom-spec | Yes | Second atom. |
| `atom-3-spec` | atom-spec | No | Third atom (for angles/torsions). |
| `atom-4-spec` | atom-spec | No | Fourth atom (for torsions). |

### Position Type: `by-coordinates`

Navigate to an arbitrary point in space (no atom or residue required).

```json
{
    "position-type": "by-coordinates",
    "label": "Unmodelled density blob",
    "position": [12.5, 34.2, -7.8]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `position` | array of 3 numbers | Yes | Cartesian coordinates [x, y, z] in Angstroms. |

## Complete Example

```json
{
    "title": "Servalcat Validation",
    "sections": [
        {
            "title": "Van der Waals Clashes",
            "items": [
                {
                    "position-type": "by-atom-spec-pair",
                    "label": "Clash: A GLN 323 O - A ILE 324 HG12, z=-6.2",
                    "atom-1-spec": ["A", 323, "", "O", ""],
                    "atom-2-spec": ["A", 324, "", "HG12", ""]
                },
                {
                    "position-type": "by-atom-spec-pair",
                    "label": "Clash: A LYS 343 CB - A GLU 344 HA, z=-5.6",
                    "atom-1-spec": ["A", 343, "", "CB", ""],
                    "atom-2-spec": ["A", 344, "", "HA", ""]
                }
            ]
        },
        {
            "title": "Rotamer Outliers",
            "items": [
                {
                    "position-type": "by-residue-spec",
                    "label": "Rotamer outlier: A PHE 42",
                    "residue-spec": ["A", 42, ""],
                    "badness": 3.5
                }
            ]
        },
        {
            "title": "Unmodelled Blobs",
            "items": [
                {
                    "position-type": "by-coordinates",
                    "label": "Blob near A ASP 200 (12.3 sigma)",
                    "position": [45.2, 12.8, -3.1]
                }
            ]
        }
    ]
}
```

## Notes

- The dialog operates on the currently active molecule. Make sure the
  relevant model is centered before loading the JSON file.
- Sections can be empty (`"items": []`) but will still appear as a group
  heading.
- The `badness` field is optional. When present, it affects the colour of
  the entry in the dialog (higher values = worse).
- Atom and residue specs must exactly match atoms in the loaded model ---
  if an atom is not found, that item will be silently skipped.
- Insertion codes and alternate conformations should be empty strings
  (`""`) when not applicable, not omitted from the array.
