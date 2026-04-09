#!/usr/bin/env python3
"""Merge one or more cod-types JSON databases into a single binary file.

Later files override earlier ones for overlapping comp_ids.

Usage: python3 json-to-bin.py -o output.bin file1.json [file2.json ...]
"""

import json
import os
import struct
import sys
import argparse


def main():
   parser = argparse.ArgumentParser(
      description='Merge JSON cod-types databases into a single binary file')
   parser.add_argument('json_files', nargs='+', help='Input JSON files')
   parser.add_argument('-o', '--output', required=True, help='Output binary file')
   args = parser.parse_args()

   # Merge: later files win
   db = {}
   for json_path in args.json_files:
      with open(json_path) as f:
         data = json.load(f)
      n_before = len(db)
      n_new = len(set(data.keys()) - set(db.keys()))
      n_override = len(set(data.keys()) & set(db.keys()))
      db.update(data)
      print(f"{json_path}: {len(data)} molecules ({n_new} new, {n_override} overridden)")

   # Build dictionary
   all_types = set()
   for types in db.values():
      all_types.update(types)
   type_list = sorted(all_types)
   type_to_id = {t: i for i, t in enumerate(type_list)}

   # Write binary
   with open(args.output, 'wb') as f:
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

   bin_size = os.path.getsize(args.output)
   print(f"\nMerged: {len(db)} molecules, {len(type_list)} unique types")
   print(f"Written to {args.output} ({bin_size:,} bytes, {bin_size/1024/1024:.2f} MB)")


if __name__ == '__main__':
   main()
