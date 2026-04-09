#!/usr/bin/env python3
"""
Quick-and-dirty molecule comparison using Acedrg atom type tables.

Each molecule is represented by its sorted set of unique non-H Acedrg atom types
extracted from acedrg-generated CIF files. Comparison uses Jaccard similarity
on these sets, with comm-3 style detailed output.
"""

import os
import re
import sys
import json
import struct
import time
from pathlib import Path


def _extract_elements(s):
    """Extract element symbols from a string like 'C[6a]C[6a]H' → ['C','C','H'].

    Skips over [...] ring annotations.
    """
    elements = []
    i = 0
    while i < len(s):
        if s[i] == '[':
            # Skip ring annotation
            j = s.index(']', i)
            i = j + 1
        elif s[i].isupper():
            # Start of element symbol
            elem = s[i]
            i += 1
            if i < len(s) and s[i].islower() and s[i] != '[' and not s[i:i+1] == '[':
                # Check it's actually a lowercase letter part of element (e.g. Cl, Br)
                # not the start of a [...] annotation
                if i < len(s) and s[i].islower():
                    elem += s[i]
                    i += 1
            elements.append(elem)
        else:
            i += 1
    return elements


def _parse_neighbor_groups(atype):
    """Parse an atom type into (central_atom, [(sorted_elements, count), ...]).

    E.g. 'C[6a](C[6a]C[6a]H)2(OC){1|C<3>,2|H<1>}'
    → ('C[6a]', [('CCH', 2), ('CO', 1)])
    """
    # Strip environment
    idx = atype.find('{')
    if idx >= 0:
        atype = atype[:idx]

    # Extract central atom (element + optional ring annotation)
    m = re.match(r'^([A-Z][a-z]?(?:\[[^\]]+\])?)', atype)
    if not m:
        return atype, []
    central = m.group(1)
    rest = atype[m.end():]

    # Parse parenthesized groups with optional repeat count
    groups = []
    i = 0
    while i < len(rest):
        if rest[i] == '(':
            # Find matching close paren
            depth = 1
            j = i + 1
            while j < len(rest) and depth > 0:
                if rest[j] == '(':
                    depth += 1
                elif rest[j] == ')':
                    depth -= 1
                j += 1
            group_content = rest[i+1:j-1]
            # Check for repeat count after ')'
            count = 1
            k = j
            while k < len(rest) and rest[k].isdigit():
                k += 1
            if k > j:
                count = int(rest[j:k])
            # Extract and sort elements from group
            elements = _extract_elements(group_content)
            sorted_elems = ''.join(sorted(elements))
            groups.append((sorted_elems, count))
            i = k
        else:
            i += 1

    return central, groups


def _normalize_atom_type(atype):
    """Normalize: keep central atom+ring, sort elements within each neighbor group,
    sort groups, keep repeat counts. Strip environment.

    E.g. 'C[6a](C[6a]C[6a]H)2(OC){1|C<3>,2|H<1>}'
    → 'C[6a](CCH)2(CO)'
    """
    central, groups = _parse_neighbor_groups(atype)
    if not groups:
        return central
    # Sort groups by (elements, count) for canonical form
    groups.sort()
    parts = []
    for elems, count in groups:
        if count > 1:
            parts.append(f"({elems}){count}")
        else:
            parts.append(f"({elems})")
    return central + ''.join(parts)


def simplify_atom_type(atype, level='full'):
    """Simplify an atom type string to a given level of detail.

    Levels:
      'full'       — no change (e.g. C[6a](C[6a]C[6a]H)2(OC){1|C<3>,2|H<1>})
      'no_env'     — strip {…} environment  (e.g. C[6a](C[6a]C[6a]H)2(OC))
      'normalized' — strip env + ring from neighbors, sort groups
                     (e.g. C[6a](CCH)2(CO))
      'core'       — element + ring only    (e.g. C[6a])
    """
    if level == 'full':
        return atype
    if level == 'no_env':
        idx = atype.find('{')
        if idx >= 0:
            return atype[:idx]
        return atype
    if level == 'normalized':
        return _normalize_atom_type(atype)
    if level == 'core':
        m = re.match(r'^([A-Z][a-z]?(?:\[[^\]]+\])?)', atype)
        return m.group(1) if m else atype
    return atype


def parse_cif_acedrg_types(path, level='full'):
    """Extract sorted unique non-H atom types from an acedrg CIF file.

    Parses the _chem_comp_acedrg loop and returns (comp_id, sorted_types).
    Returns (None, []) if no acedrg table found.
    level: 'full', 'no_env', or 'core' — controls atom type simplification.
    """
    comp_id = None
    atom_types = set()
    in_acedrg_loop = False
    acedrg_columns = []
    atom_type_col = -1
    comp_id_col = -1

    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            stripped = line.strip()

            # Grab comp_id from data_ block
            if stripped.startswith('data_comp_') and not stripped.startswith('data_comp_list'):
                comp_id = stripped[len('data_comp_'):]

            # Detect start of the acedrg loop
            if stripped == 'loop_':
                in_acedrg_loop = False
                acedrg_columns = []
                atom_type_col = -1
                comp_id_col = -1
                # Peek: the next _chem_comp_acedrg lines will tell us
                continue

            if stripped.startswith('_chem_comp_acedrg.'):
                in_acedrg_loop = True
                col_name = stripped.split('.')[1]
                acedrg_columns.append(col_name)
                if col_name == 'atom_type':
                    atom_type_col = len(acedrg_columns) - 1
                if col_name == 'comp_id':
                    comp_id_col = len(acedrg_columns) - 1
                continue

            if in_acedrg_loop:
                if stripped == '' or stripped.startswith('#') or stripped.startswith('loop_') or stripped.startswith('_'):
                    in_acedrg_loop = False
                    continue

                if atom_type_col < 0:
                    continue

                # Parse the data line — fields are space-separated,
                # but atom_type contains spaces inside its notation.
                # The atom_type is always the last column (col index 2 of 3 columns:
                # comp_id, atom_id, atom_type).
                # The atom_type runs from after the second whitespace-delimited token to EOL.
                parts = stripped.split(None, len(acedrg_columns) - 1)
                if len(parts) >= len(acedrg_columns):
                    atype = parts[atom_type_col].strip()
                    # Skip hydrogen types
                    if not atype.startswith('H(') and atype != 'H':
                        atom_types.add(simplify_atom_type(atype, level))

    return comp_id, sorted(atom_types)


def scan_monomers(directory, level='full', verbose=True):
    """Walk a monomers directory and parse all CIFs with acedrg tables.

    Returns dict: {comp_id: sorted_types_list}
    """
    db = {}
    n_files = 0
    n_with_acedrg = 0

    monomers_path = Path(directory)
    # Monomers are typically organized as monomers/A/ALA.cif, monomers/B/... etc.
    cif_files = sorted(monomers_path.rglob('*.cif'))

    if verbose:
        print(f"Found {len(cif_files)} CIF files in {directory}")

    for cif_path in cif_files:
        n_files += 1
        comp_id, types = parse_cif_acedrg_types(str(cif_path), level)
        if types:
            n_with_acedrg += 1
            if comp_id:
                db[comp_id] = types

        if verbose and n_files % 2000 == 0:
            print(f"  processed {n_files} files, {n_with_acedrg} with acedrg types...")

    if verbose:
        print(f"Done: {n_files} files, {n_with_acedrg} with acedrg tables, {len(db)} unique comp_ids")

    return db


def comm3_compare(set1, set2):
    """Compare two sorted lists like comm -3.

    Returns (only_in_1, only_in_2, common) as sorted lists.
    """
    s1 = set(set1)
    s2 = set(set2)
    common = sorted(s1 & s2)
    only1 = sorted(s1 - s2)
    only2 = sorted(s2 - s1)
    return only1, only2, common


def jaccard(set1, set2):
    """Jaccard similarity: |intersection| / |union|."""
    s1 = set(set1)
    s2 = set(set2)
    if not s1 and not s2:
        return 1.0
    intersection = len(s1 & s2)
    union = len(s1 | s2)
    return intersection / union if union > 0 else 0.0


def search(query_types, db, top_n=20):
    """Compare query against all molecules in db, return top-N by Jaccard.

    Returns list of (comp_id, jaccard_score, n_common, n_only_query, n_only_ref).
    """
    query_set = set(query_types)
    results = []
    for comp_id, ref_types in db.items():
        ref_set = set(ref_types)
        common = len(query_set & ref_set)
        union = len(query_set | ref_set)
        j = common / union if union > 0 else 0.0
        results.append((comp_id, j, common, len(query_set - ref_set), len(ref_set - query_set)))

    results.sort(key=lambda x: -x[1])
    return results[:top_n]


def build_dictionary(db):
    """Build a global dictionary of all unique atom type strings.

    Returns (type_to_id, id_to_type) where IDs are 0-based integers.
    """
    all_types = set()
    for types in db.values():
        all_types.update(types)
    sorted_types = sorted(all_types)
    type_to_id = {t: i for i, t in enumerate(sorted_types)}
    return type_to_id, sorted_types


def encode_db(db, type_to_id):
    """Encode the database as compact integer arrays.

    Returns dict: {comp_id: sorted_list_of_type_ids}
    """
    encoded = {}
    for comp_id, types in db.items():
        encoded[comp_id] = sorted(type_to_id[t] for t in types)
    return encoded


def save_db_json(db, type_to_id, id_to_type, path):
    """Save the database as JSON for prototype use."""
    data = {
        'dictionary': id_to_type,
        'molecules': {comp_id: sorted(type_to_id[t] for t in types)
                      for comp_id, types in db.items()}
    }
    with open(path, 'w') as f:
        json.dump(data, f, separators=(',', ':'))
    return os.path.getsize(path)


def save_db_binary(db, type_to_id, id_to_type, path):
    """Save the database in a compact binary format.

    Format:
      Header: magic(4) + n_types(4) + n_molecules(4)
      Dictionary: for each type: len(2) + utf8_bytes
      Molecules: for each molecule: comp_id_len(1) + comp_id + n_types(2) + type_ids(2 each)
    """
    with open(path, 'wb') as f:
        n_types = len(id_to_type)
        n_molecules = len(db)
        # Header
        f.write(b'ACDR')
        f.write(struct.pack('<II', n_types, n_molecules))

        # Dictionary
        for t in id_to_type:
            encoded = t.encode('utf-8')
            f.write(struct.pack('<H', len(encoded)))
            f.write(encoded)

        # Molecules
        for comp_id, types in sorted(db.items()):
            cid = comp_id.encode('utf-8')
            f.write(struct.pack('<B', len(cid)))
            f.write(cid)
            ids = sorted(type_to_id[t] for t in types)
            f.write(struct.pack('<H', len(ids)))
            for tid in ids:
                f.write(struct.pack('<H', tid))

    return os.path.getsize(path)


def load_db_binary(path):
    """Load a binary database.

    Returns (db, id_to_type) where db = {comp_id: [type_strings]}.
    """
    with open(path, 'rb') as f:
        magic = f.read(4)
        assert magic == b'ACDR', f"Bad magic: {magic}"
        n_types, n_molecules = struct.unpack('<II', f.read(8))

        id_to_type = []
        for _ in range(n_types):
            slen = struct.unpack('<H', f.read(2))[0]
            id_to_type.append(f.read(slen).decode('utf-8'))

        db = {}
        for _ in range(n_molecules):
            cid_len = struct.unpack('<B', f.read(1))[0]
            comp_id = f.read(cid_len).decode('utf-8')
            nt = struct.unpack('<H', f.read(2))[0]
            ids = [struct.unpack('<H', f.read(2))[0] for _ in range(nt)]
            db[comp_id] = [id_to_type[i] for i in ids]

    return db, id_to_type


def report_sizes(db, type_to_id, id_to_type):
    """Print data size analysis."""
    n_mol = len(db)
    n_types = len(id_to_type)

    # Raw string size
    raw_bytes = 0
    total_type_entries = 0
    for comp_id, types in db.items():
        raw_bytes += len(comp_id)
        for t in types:
            raw_bytes += len(t)
        total_type_entries += len(types)

    # Dictionary size
    dict_bytes = sum(2 + len(t.encode('utf-8')) for t in id_to_type)

    # Integer array size (2 bytes per ID + 3 bytes per molecule header)
    int_array_bytes = n_mol * 3 + total_type_entries * 2

    avg_types = total_type_entries / n_mol if n_mol > 0 else 0

    print(f"\n--- Size Analysis ---")
    print(f"Molecules:           {n_mol}")
    print(f"Unique atom types:   {n_types}")
    print(f"Avg types/molecule:  {avg_types:.1f}")
    print(f"Raw string data:     {raw_bytes:,} bytes ({raw_bytes/1024:.1f} KB)")
    print(f"Dictionary:          {dict_bytes:,} bytes ({dict_bytes/1024:.1f} KB)")
    print(f"Integer arrays:      {int_array_bytes:,} bytes ({int_array_bytes/1024:.1f} KB)")
    print(f"Total (binary):      {dict_bytes + int_array_bytes + 12:,} bytes "
          f"({(dict_bytes + int_array_bytes + 12)/1024:.1f} KB)")

    if n_mol > 0:
        scale = 50000 / n_mol
        est = (dict_bytes + int_array_bytes * scale + 12)
        print(f"\nExtrapolated to 50k molecules:")
        print(f"  Dictionary:        {dict_bytes:,} bytes (same)")
        print(f"  Integer arrays:    {int(int_array_bytes * scale):,} bytes "
              f"({int_array_bytes * scale / 1024 / 1024:.1f} MB)")
        print(f"  Total:             {int(est):,} bytes ({est/1024/1024:.1f} MB)")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Acedrg atom type molecule similarity')
    sub = parser.add_subparsers(dest='command')

    # Common argument for atom type simplification
    level_help = ('Atom type detail level: '
                  'full (default), no_env (strip {…} environment), core (element+ring only)')

    # parse: parse a single CIF
    p_parse = sub.add_parser('parse', help='Parse a single CIF file')
    p_parse.add_argument('cif', help='Path to acedrg CIF file')
    p_parse.add_argument('--level', default='full', choices=['full', 'no_env', 'normalized', 'core'], help=level_help)

    # build: scan monomers and build DB
    p_build = sub.add_parser('build', help='Build reference DB from monomers directory')
    p_build.add_argument('directory', help='Path to monomers directory')
    p_build.add_argument('-o', '--output', default='acedrg-refdb.bin', help='Output binary DB file')
    p_build.add_argument('--json', help='Also save as JSON (for inspection)')
    p_build.add_argument('--level', default='full', choices=['full', 'no_env', 'normalized', 'core'], help=level_help)

    # search: query against a DB
    p_search = sub.add_parser('search', help='Search for similar molecules')
    p_search.add_argument('query_cif', help='Query CIF file')
    p_search.add_argument('db', help='Reference DB file (.bin)')
    p_search.add_argument('-n', '--top', type=int, default=20, help='Number of top hits')
    p_search.add_argument('-v', '--verbose', action='store_true', help='Show comm-3 details for top hits')
    p_search.add_argument('--level', default='full', choices=['full', 'no_env', 'normalized', 'core'], help=level_help)

    # compare: quick comparison across all 3 levels for one query
    p_cmp = sub.add_parser('compare-levels', help='Compare similarity at all 3 detail levels')
    p_cmp.add_argument('query_cif', help='Query CIF file')
    p_cmp.add_argument('directory', help='Monomers directory')
    p_cmp.add_argument('-n', '--top', type=int, default=10, help='Number of top hits')

    args = parser.parse_args()

    if args.command == 'parse':
        comp_id, types = parse_cif_acedrg_types(args.cif, args.level)
        print(f"Comp ID: {comp_id}  (level={args.level})")
        print(f"Non-H atom types ({len(types)}):")
        for t in types:
            print(f"  {t}")

    elif args.command == 'build':
        t0 = time.time()
        db = scan_monomers(args.directory, level=args.level)
        t1 = time.time()
        print(f"Scan took {t1-t0:.1f}s")

        type_to_id, id_to_type = build_dictionary(db)
        report_sizes(db, type_to_id, id_to_type)

        size = save_db_binary(db, type_to_id, id_to_type, args.output)
        print(f"\nSaved binary DB: {args.output} ({size:,} bytes, {size/1024/1024:.2f} MB)")

        if args.json:
            jsize = save_db_json(db, type_to_id, id_to_type, args.json)
            print(f"Saved JSON DB:   {args.json} ({jsize:,} bytes, {jsize/1024/1024:.2f} MB)")

    elif args.command == 'search':
        # Load query
        q_comp_id, q_types = parse_cif_acedrg_types(args.query_cif, args.level)
        if not q_types:
            print(f"No acedrg types found in {args.query_cif}")
            sys.exit(1)
        print(f"Query: {q_comp_id} ({len(q_types)} non-H atom types)")

        # Load DB
        db, id_to_type = load_db_binary(args.db)
        print(f"Reference DB: {len(db)} molecules, {len(id_to_type)} unique types\n")

        # Search
        results = search(q_types, db, top_n=args.top)

        print(f"Top {len(results)} matches:")
        print(f"{'Rank':>4}  {'Comp ID':<8}  {'Jaccard':>7}  {'Common':>6}  {'Query-only':>10}  {'Ref-only':>8}")
        print("-" * 60)
        for i, (comp_id, j, common, only_q, only_r) in enumerate(results, 1):
            print(f"{i:>4}  {comp_id:<8}  {j:>7.4f}  {common:>6}  {only_q:>10}  {only_r:>8}")

        if args.verbose and results:
            print("\n--- Detailed comparison for top 5 ---")
            for comp_id, j, _, _, _ in results[:5]:
                ref_types = db[comp_id]
                only_q, only_r, common = comm3_compare(q_types, ref_types)
                print(f"\n{q_comp_id} vs {comp_id} (Jaccard={j:.4f}):")
                if common:
                    print(f"  Common ({len(common)}):")
                    for t in common:
                        print(f"    {t}")
                if only_q:
                    print(f"  Only in {q_comp_id} ({len(only_q)}):")
                    for t in only_q:
                        print(f"    {t}")
                if only_r:
                    print(f"  Only in {comp_id} ({len(only_r)}):")
                    for t in only_r:
                        print(f"    {t}")

    elif args.command == 'compare-levels':
        print(f"Comparing similarity at all 3 detail levels for query from {args.query_cif}")
        print(f"Scanning {args.directory}...\n")

        for level in ['full', 'no_env', 'normalized', 'core']:
            q_comp_id, q_types = parse_cif_acedrg_types(args.query_cif, level)
            if not q_types:
                print(f"[{level}] No acedrg types found")
                continue

            db = scan_monomers(args.directory, level=level, verbose=False)
            type_to_id, id_to_type = build_dictionary(db)

            print(f"=== Level: {level} ===")
            print(f"Query {q_comp_id}: {len(q_types)} types | "
                  f"DB: {len(db)} molecules, {len(id_to_type)} unique types")

            results = search(q_types, db, top_n=args.top)
            print(f"{'Rank':>4}  {'Comp ID':<8}  {'Jaccard':>7}  {'Common':>6}  "
                  f"{'Q-only':>6}  {'R-only':>6}")
            print("-" * 50)
            for i, (comp_id, j, common, only_q, only_r) in enumerate(results, 1):
                print(f"{i:>4}  {comp_id:<8}  {j:>7.4f}  {common:>6}  "
                      f"{only_q:>6}  {only_r:>6}")
            print()

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
