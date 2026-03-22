#!/usr/bin/env python3
"""
Speed test: just import all monomer CIF dictionaries without computing COD types.
Compare timing with compute-cod-types.py to isolate where the bottleneck is.
"""

import sys
import time
from pathlib import Path

import coot_headless_api as chapi

def find_all_cif_files(monomers_dir):
   cif_files = []
   monomers_path = Path(monomers_dir)
   for subdir in sorted(monomers_path.iterdir()):
      if subdir.is_dir():
         for cif_file in sorted(subdir.glob('*.cif')):
            cif_files.append(cif_file)
   return cif_files

def main():
   monomers_dir = sys.argv[1] if len(sys.argv) > 1 else '/Users/pemsley/Projects/monomers/monomers'

   print(f"Creating molecules container...")
   mc = chapi.molecules_container_t(False)
   imol_enc = -999999

   print(f"Scanning {monomers_dir} for CIF files...")
   cif_files = find_all_cif_files(monomers_dir)
   print(f"Found {len(cif_files)} CIF files")

   n_ok = 0
   n_fail = 0
   t_start = time.time()

   for i, cif_path in enumerate(cif_files):
      status = mc.import_cif_dictionary(str(cif_path), imol_enc)
      if status:
         n_ok += 1
      else:
         n_fail += 1

      if (i + 1) % 500 == 0:
         elapsed = time.time() - t_start
         rate = (i + 1) / elapsed
         remaining = (len(cif_files) - i - 1) / rate
         print(f"  [{i+1}/{len(cif_files)}] "
               f"({n_ok} ok, {n_fail} fail) "
               f"{rate:.1f}/s "
               f"[{elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining]")

   elapsed = time.time() - t_start
   print(f"\nDone in {elapsed:.1f}s ({len(cif_files)/elapsed:.1f}/s)")
   print(f"  OK: {n_ok}  Fail: {n_fail}")

if __name__ == '__main__':
   main()
