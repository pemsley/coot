#!/usr/bin/env python3
"""
Compute and search COD/acedrg atom types for CCD monomer entries using coot_headless_api.

Subcommands:
   build   — Batch-compute types for all monomers, produce JSON + binary DB
   search  — Compute types for query molecules, search against a pre-built DB

Usage:
   python3 compute-cod-types.py build /path/to/monomers [-o cod-types-db] [--level normalized]
   python3 compute-cod-types.py search "TYR,PHE" --db cod-types-db-normalized.json --monomers /path/to/monomers
   python3 compute-cod-types.py search "04*" --db cod-types-db-normalized.json --monomers /path/to/monomers
"""

import os
import re
import sys
import json
import time
import struct
import fnmatch
import argparse
from pathlib import Path


# ---------------------------------------------------------------------------
# Atom type simplification
# ---------------------------------------------------------------------------

def _extract_elements(s):
   """Extract element symbols from a string like 'C[6a]C[6a]H' → ['C','C','H']."""
   elements = []
   i = 0
   while i < len(s):
      if s[i] == '[':
         j = s.index(']', i)
         i = j + 1
      elif s[i].isupper():
         elem = s[i]
         i += 1
         if i < len(s) and s[i].islower() and s[i] != '[':
            elem += s[i]
            i += 1
         elements.append(elem)
      else:
         i += 1
   return elements


def _parse_neighbor_groups(atype):
   """Parse an atom type into (central_atom, [(sorted_elements, count), ...])."""
   idx = atype.find('{')
   if idx >= 0:
      atype = atype[:idx]

   m = re.match(r'^([A-Z][a-z]?(?:\[[^\]]+\])?)', atype)
   if not m:
      return atype, []
   central = m.group(1)
   rest = atype[m.end():]

   groups = []
   i = 0
   while i < len(rest):
      if rest[i] == '(':
         depth = 1
         j = i + 1
         while j < len(rest) and depth > 0:
            if rest[j] == '(':
               depth += 1
            elif rest[j] == ')':
               depth -= 1
            j += 1
         group_content = rest[i+1:j-1]
         count = 1
         k = j
         while k < len(rest) and rest[k].isdigit():
            k += 1
         if k > j:
            count = int(rest[j:k])
         elements = _extract_elements(group_content)
         sorted_elems = ''.join(sorted(elements))
         groups.append((sorted_elems, count))
         i = k
      else:
         i += 1
   return central, groups


def _normalize_atom_type(atype):
   """Normalize: keep central atom+ring, sort elements within neighbor groups,
   sort groups, keep repeat counts. Strip environment."""
   central, groups = _parse_neighbor_groups(atype)
   if not groups:
      return central
   groups.sort()
   parts = []
   for elems, count in groups:
      if count > 1:
         parts.append(f"({elems}){count}")
      else:
         parts.append(f"({elems})")
   return central + ''.join(parts)


def simplify_atom_type(atype, level='full'):
   """Simplify an atom type string to a given level of detail."""
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


# ---------------------------------------------------------------------------
# Binary format
# ---------------------------------------------------------------------------

def save_db_binary(db, path):
   """Save reference DB in compact binary format.

   Format: magic(4) + n_types(4) + n_molecules(4)
           + dictionary entries + molecule entries
   """
   all_types = set()
   for types in db.values():
      all_types.update(types)
   type_list = sorted(all_types)
   type_to_id = {t: i for i, t in enumerate(type_list)}

   with open(path, 'wb') as f:
      f.write(b'ACDR')
      f.write(struct.pack('<II', len(type_list), len(db)))
      for t in type_list:
         t_bytes = t.encode('utf-8')
         f.write(struct.pack('<H', len(t_bytes)))
         f.write(t_bytes)
      for comp_id, types in sorted(db.items()):
         cid_bytes = comp_id.encode('utf-8')
         f.write(struct.pack('<B', len(cid_bytes)))
         f.write(cid_bytes)
         ids = sorted(type_to_id[t] for t in types)
         f.write(struct.pack('<H', len(ids)))
         for tid in ids:
            f.write(struct.pack('<H', tid))

   return len(type_list)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_all_cif_files(monomers_dir):
   """Walk the monomers directory and return all .cif file paths."""
   cif_files = []
   monomers_path = Path(monomers_dir)
   for subdir in sorted(monomers_path.iterdir()):
      if subdir.is_dir():
         for cif_file in sorted(subdir.glob('*.cif')):
            cif_files.append(cif_file)
   return cif_files


def find_cif_for_comp_id(monomers_dir, comp_id):
   """Find the CIF file for a compound ID in the monomers directory.

   Convention: monomers/<first-char>/<comp_id>.cif
   """
   first_char = comp_id[0].upper()
   cif_path = Path(monomers_dir) / first_char / f"{comp_id}.cif"
   if cif_path.exists():
      return cif_path
   return None


def compute_types_for_compound(mc, comp_id, cif_path, imol_enc, level):
   """Import a CIF and compute simplified non-H atom types for one compound.

   Returns sorted list of type strings, or None on failure.
   """
   status = mc.import_cif_dictionary(str(cif_path), imol_enc)
   if status == 0:
      return None

   try:
      pairs = mc.get_computed_acedrg_atom_types(comp_id, imol_enc)
   except RuntimeError as e:
      print(f"  FAIL {comp_id}: {e}", file=sys.stderr)
      return None

   if len(pairs) == 0:
      return None

   types = set()
   for atom_name, cod_type in pairs:
      if cod_type.startswith('H(') or cod_type == 'H':
         continue
      types.add(simplify_atom_type(cod_type, level))

   return sorted(types) if types else None


def resolve_query(query_str, db_comp_ids):
   """Resolve a query string to a list of comp_ids.

   query_str can be:
     - single comp_id: "TYR"
     - comma-separated: "TYR,PHE,ALA"
     - glob pattern: "04*", "A?E"
   """
   if ',' in query_str:
      return [c.strip() for c in query_str.split(',')]

   # Try as glob pattern against known comp_ids
   matches = sorted(cid for cid in db_comp_ids if fnmatch.fnmatch(cid, query_str))
   if matches:
      return matches

   # Single comp_id (may or may not be in DB)
   return [query_str]


# ---------------------------------------------------------------------------
# build subcommand
# ---------------------------------------------------------------------------

def cmd_build(args):
   print(f"Importing coot_headless_api...")
   import coot_headless_api as chapi

   print(f"Creating molecules container...")
   mc = chapi.molecules_container_t(False)

   imol_enc = -999999  # IMOL_ENC_ANY
   level = args.level

   print(f"Scanning {args.monomers_dir} for CIF files...")
   cif_files = find_all_cif_files(args.monomers_dir)
   print(f"Found {len(cif_files)} CIF files")

   if args.limit > 0:
      cif_files = cif_files[:args.limit]
      print(f"Processing first {args.limit} entries only")

   db = {}
   n_success = 0
   n_fail = 0
   n_empty = 0
   t_start = time.time()

   for i, cif_path in enumerate(cif_files):
      comp_id = cif_path.stem

      if (i + 1) % 500 == 0:
         elapsed = time.time() - t_start
         rate = (i + 1) / elapsed
         remaining = (len(cif_files) - i - 1) / rate
         print(f"  [{i+1}/{len(cif_files)}] {comp_id} "
               f"({n_success} ok, {n_fail} fail, {n_empty} empty) "
               f"[{elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining]")

      types = compute_types_for_compound(mc, comp_id, cif_path, imol_enc, level)
      if types is None:
         status = mc.import_cif_dictionary(str(cif_path), imol_enc)
         if status == 0:
            n_fail += 1
         else:
            n_empty += 1
         continue

      if types:
         db[comp_id] = types
         n_success += 1
         if args.verbose:
            print(f"  OK {comp_id}: {len(types)} non-H types")
      else:
         n_empty += 1

   elapsed = time.time() - t_start
   print(f"\nDone in {elapsed:.1f}s")
   print(f"  Success: {n_success}")
   print(f"  Failed:  {n_fail}")
   print(f"  Empty:   {n_empty}")

   # Save JSON
   json_path = f"{args.output}-{level}.json"
   with open(json_path, 'w') as f:
      json.dump(db, f, indent=1, sort_keys=True)
   json_size = os.path.getsize(json_path)
   print(f"\nJSON saved to {json_path} ({json_size/1024/1024:.2f} MB)")

   # Save binary
   bin_path = f"{args.output}-{level}.bin"
   n_unique = save_db_binary(db, bin_path)
   bin_size = os.path.getsize(bin_path)
   print(f"Binary saved to {bin_path} ({bin_size/1024/1024:.2f} MB)")
   print(f"Unique atom types: {n_unique}")

   if n_success > 0:
      avg_types = sum(len(v) for v in db.values()) / n_success
      print(f"\nAverage types per molecule: {avg_types:.1f}")
      projected_50k = bin_size * (50000 / n_success)
      print(f"Projected size for 50k molecules: {projected_50k/1024/1024:.1f} MB")


# ---------------------------------------------------------------------------
# search subcommand
# ---------------------------------------------------------------------------

def cmd_search(args):
   # Load reference database
   with open(args.db) as f:
      db = json.load(f)
   print(f"Reference DB: {len(db)} molecules")

   level = args.level
   query_ids = resolve_query(args.query, db.keys())
   print(f"Query: {len(query_ids)} compound(s)")

   # For each query, get its types: from DB if available, else compute
   need_compute = [qid for qid in query_ids if qid not in db]
   query_types = {}

   # Use pre-existing types from DB where available
   for qid in query_ids:
      if qid in db:
         query_types[qid] = db[qid]

   # Compute types for compounds not in the DB
   if need_compute:
      if not args.monomers:
         print(f"ERROR: {len(need_compute)} compound(s) not in DB "
               f"and no --monomers directory specified: {need_compute}",
               file=sys.stderr)
         sys.exit(1)

      print(f"Computing types for {len(need_compute)} compound(s) not in DB...")
      import coot_headless_api as chapi
      mc = chapi.molecules_container_t(False)
      imol_enc = -999999

      for comp_id in need_compute:
         cif_path = find_cif_for_comp_id(args.monomers, comp_id)
         if cif_path is None:
            print(f"  WARNING: no CIF found for {comp_id}", file=sys.stderr)
            continue
         types = compute_types_for_compound(mc, comp_id, cif_path, imol_enc, level)
         if types:
            query_types[comp_id] = types
         else:
            print(f"  WARNING: no types computed for {comp_id}", file=sys.stderr)

   if not query_types:
      print("No query compounds with valid types.", file=sys.stderr)
      sys.exit(1)

   # Search each query against the full database
   top_n = args.top

   for qid in sorted(query_types.keys()):
      q_set = set(query_types[qid])
      scores = []
      for rid, ref_types in db.items():
         if rid == qid:
            continue
         r_set = set(ref_types)
         common = len(q_set & r_set)
         if common == 0:
            continue
         union = len(q_set | r_set)
         scores.append((common / union, rid, common, len(q_set - r_set), len(r_set - q_set)))

      scores.sort(reverse=True)
      hits = scores[:top_n]

      print(f"\n{qid} ({len(q_set)} types):")
      print(f"  {'Rank':>4}  {'Hit':<8}  {'Jaccard':>7}  {'Common':>6}  {'Q-only':>6}  {'R-only':>6}")
      print(f"  {'-'*50}")
      for rank, (j, rid, common, only_q, only_r) in enumerate(hits, 1):
         print(f"  {rank:>4}  {rid:<8}  {j:>7.4f}  {common:>6}  {only_q:>6}  {only_r:>6}")

   # Also write space-separated output file
   out_path = args.output
   if out_path:
      with open(out_path, 'w') as f:
         for qid in sorted(query_types.keys()):
            q_set = set(query_types[qid])
            scores = []
            for rid, ref_types in db.items():
               if rid == qid:
                  continue
               r_set = set(ref_types)
               common = len(q_set & r_set)
               if common == 0:
                  continue
               union = len(q_set | r_set)
               scores.append((common / union, rid))
            scores.sort(reverse=True)
            for j, rid in scores[:top_n]:
               f.write(f"{qid} {rid} {j:.4f}\n")
      print(f"\nResults written to {out_path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
   parser = argparse.ArgumentParser(
      description='Compute and search COD/acedrg atom types for CCD monomers')
   sub = parser.add_subparsers(dest='command')

   level_help = 'Atom type simplification level (default: normalized)'
   level_choices = ['full', 'no_env', 'normalized', 'core']

   # build
   p_build = sub.add_parser('build',
      help='Batch-compute types for all monomers, produce JSON + binary DB')
   p_build.add_argument('monomers_dir',
      help='Path to monomers directory')
   p_build.add_argument('--output', '-o', default='cod-types-db',
      help='Output filename prefix (default: cod-types-db)')
   p_build.add_argument('--level', '-l', default='normalized',
      choices=level_choices, help=level_help)
   p_build.add_argument('--limit', type=int, default=0,
      help='Process only first N entries (0 = all)')
   p_build.add_argument('--verbose', '-v', action='store_true',
      help='Print progress for each compound')

   # search
   p_search = sub.add_parser('search',
      help='Search for similar molecules in a pre-built DB')
   p_search.add_argument('query',
      help='Comma-separated comp_ids (e.g. "TYR,PHE") or glob pattern (e.g. "04*")')
   p_search.add_argument('--db', required=True,
      help='Reference DB JSON file')
   p_search.add_argument('--monomers',
      help='Monomers directory (needed for compounds not in DB)')
   p_search.add_argument('-n', '--top', type=int, default=20,
      help='Number of top hits per query (default: 20)')
   p_search.add_argument('--level', '-l', default='normalized',
      choices=level_choices, help=level_help)
   p_search.add_argument('--output', '-o', default=None,
      help='Write space-separated results to file')

   args = parser.parse_args()

   if args.command == 'build':
      cmd_build(args)
   elif args.command == 'search':
      cmd_search(args)
   else:
      parser.print_help()


if __name__ == '__main__':
   main()
