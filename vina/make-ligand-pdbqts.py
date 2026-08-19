#!/usr/bin/env python3
"""Batch-prepare flexible-ligand PDBQT files from acedrg restraints dictionaries.

Each input .cif is an acedrg restraints dictionary (they may all share the same
three-letter code, e.g. "LIG"). For each one the ligand is built from the
dictionary and written as <cif-basename>.pdbqt in the output directory. Outputs
are named from the cif file name (not the comp_id) so that dictionaries that
share a comp_id do not clash.

Usage:
    python make-ligand-pdbqts.py <out_dir> <ligand.cif> [more.cif ...]
    python make-ligand-pdbqts.py <out_dir> "dir/*.cif"

Set COOT_HEADLESS_API_DIR to the build directory that has coot_headless_api*.so,
or make sure it is otherwise importable.
"""

import os
import sys
import glob


def import_chapi():
    d = os.environ.get("COOT_HEADLESS_API_DIR")
    if d and os.path.isdir(d):
        sys.path.insert(0, d)
    try:
        import coot_headless_api
        return coot_headless_api
    except ImportError as e:
        sys.exit("cannot import coot_headless_api (%s); set COOT_HEADLESS_API_DIR" % e)


def comp_id_from_cif(path):
    """The ligand three-letter code, from the non-list data_comp_ block."""
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("data_comp_") and line != "data_comp_list":
                    return line[len("data_comp_"):]
    except IOError:
        pass
    return "LIG"


def expand(args):
    cifs = []
    for a in args:
        if any(c in a for c in "*?["):
            cifs.extend(sorted(glob.glob(a)))
        else:
            cifs.append(a)
    return cifs


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: make-ligand-pdbqts.py <out_dir> <ligand.cif> [more.cif ...]")

    out_dir = sys.argv[1]
    cifs = expand(sys.argv[2:])
    if not cifs:
        sys.exit("no cif files found")
    os.makedirs(out_dir, exist_ok=True)

    chapi = import_chapi()
    mc = chapi.molecules_container_t(True)
    mc.set_use_gemmi(False)
    enc = mc.get_imol_enc_any()

    n_ok = 0
    failures = []
    for cif in cifs:
        comp_id = comp_id_from_cif(cif)
        base = os.path.splitext(os.path.basename(cif))[0]
        out = os.path.join(out_dir, base + ".pdbqt")
        try:
            mc.import_cif_dictionary(cif, enc)   # re-imports comp_id for this ligand
            imol = mc.get_monomer(comp_id)
            if imol < 0:
                failures.append((base, "get_monomer('%s') failed" % comp_id))
                continue
            n = mc.export_ligand_as_pdbqt(imol, "/*/*/*", out)
            if n > 0:
                n_ok += 1
                print("[ok]   %-20s -> %s (%d atoms)" % (base, out, n))
            else:
                failures.append((base, "export wrote 0 atoms"))
        except Exception as e:
            failures.append((base, str(e)))

    print("\ndone: %d/%d ligands prepared into %s" % (n_ok, len(cifs), out_dir))
    if failures:
        print("failures (%d):" % len(failures))
        for base, why in failures:
            print("  %-20s %s" % (base, why))


if __name__ == "__main__":
    main()
