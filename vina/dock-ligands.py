#!/usr/bin/env python3
"""Dock a set of ligand PDBQT files against a receptor with AutoDock Vina.

Runs vina once per ligand, saves the docked poses, and writes a CSV of the best
affinity per ligand, sorted best (most negative) first.

Usage:
    python dock-ligands.py --receptor rec.pdbqt --ligands "ligs/*.pdbqt" \
        --out results_dir  (--config box.txt | --center X Y Z --size SX SY SZ)

Options:
    --vina PATH          vina executable (default: "vina" on PATH)
    --exhaustiveness N   vina search exhaustiveness (default 8)
    --cpu N              cpus per vina process (default 1)
    --jobs N             vina processes to run concurrently (default 1)
    --seed N             random seed (optional, for reproducibility)

Notes:
  - A box centre with negative coordinates must not be passed on the vina command
    line (vina mis-parses the leading '-'), so this script always writes the box to
    a small config file and passes it with --config.
"""

import os
import sys
import glob
import csv
import argparse
import tempfile
import subprocess
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed


def parse_args():
    p = argparse.ArgumentParser(description="Batch Vina docking of ligand PDBQT files.")
    p.add_argument("--receptor", required=True, help="receptor .pdbqt")
    p.add_argument("--ligands", required=True, help="ligand .pdbqt files (glob, quote it)")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--config", help="vina box config file (center_x/size_x/... lines)")
    p.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"))
    p.add_argument("--size", nargs=3, type=float, metavar=("SX", "SY", "SZ"))
    p.add_argument("--vina", default=shutil.which("vina") or "vina")
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--cpu", type=int, default=1)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--seed", type=int, default=None)
    return p.parse_args()


def resolve_box(args, out_dir):
    """Return a path to a vina box config file (writing one if needed)."""
    if args.config:
        return args.config
    if args.center and args.size:
        path = os.path.join(out_dir, "_box.txt")
        with open(path, "w") as f:
            f.write("center_x = %.3f\ncenter_y = %.3f\ncenter_z = %.3f\n" % tuple(args.center))
            f.write("size_x = %.3f\nsize_y = %.3f\nsize_z = %.3f\n" % tuple(args.size))
        return path
    sys.exit("supply a box: --config <box.txt> or --center X Y Z --size SX SY SZ")


def best_affinity(out_pdbqt):
    """Best (mode 1) affinity from the first REMARK VINA RESULT line, or None."""
    try:
        with open(out_pdbqt) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except (IOError, ValueError, IndexError):
        pass
    return None


def dock_one(args, config_path, ligand):
    base = os.path.splitext(os.path.basename(ligand))[0]
    out_pdbqt = os.path.join(args.out, base + "_out.pdbqt")
    cmd = [args.vina, "--receptor", args.receptor, "--ligand", ligand,
           "--config", config_path, "--out", out_pdbqt,
           "--exhaustiveness", str(args.exhaustiveness), "--cpu", str(args.cpu)]
    if args.seed is not None:
        cmd += ["--seed", str(args.seed)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        return (base, None, "timeout")
    if r.returncode != 0:
        err = ""
        for ln in r.stdout.splitlines():
            if "error" in ln.lower() or "Error" in ln:
                err = ln.strip()
        return (base, None, err or "vina returned %d" % r.returncode)
    aff = best_affinity(out_pdbqt)
    if aff is None:
        return (base, None, "no affinity in output")
    return (base, aff, None)


def main():
    args = parse_args()
    ligands = sorted(glob.glob(args.ligands))
    if not ligands:
        sys.exit("no ligand pdbqt files matched: " + args.ligands)
    if not os.path.exists(args.vina):
        sys.exit("vina executable not found: %s (use --vina)" % args.vina)
    os.makedirs(args.out, exist_ok=True)
    config_path = resolve_box(args, args.out)

    print("docking %d ligands against %s (vina exhaustiveness %d, %d job(s) x %d cpu)"
          % (len(ligands), args.receptor, args.exhaustiveness, args.jobs, args.cpu))

    results = []
    done = 0
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(dock_one, args, config_path, lig): lig for lig in ligands}
        for fut in as_completed(futs):
            base, aff, err = fut.result()
            results.append((base, aff, err))
            done += 1
            if err:
                print("[%3d/%d] %-20s FAILED: %s" % (done, len(ligands), base, err))
            else:
                print("[%3d/%d] %-20s %8.2f kcal/mol" % (done, len(ligands), base, aff))

    # CSV, sorted best (most negative) first; failures at the end
    csv_path = os.path.join(args.out, "docking_scores.csv")
    ranked = sorted(results, key=lambda r: (r[1] is None, r[1] if r[1] is not None else 0.0))
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "ligand", "affinity_kcal_per_mol", "error"])
        for i, (base, aff, err) in enumerate(ranked, start=1):
            w.writerow([i, base, "" if aff is None else "%.2f" % aff, err or ""])

    n_ok = sum(1 for _, aff, _ in results if aff is not None)
    print("\ndone: %d/%d docked. scores in %s" % (n_ok, len(ligands), csv_path))
    top = [r for r in ranked if r[1] is not None][:10]
    if top:
        print("top hits:")
        for i, (base, aff, _) in enumerate(top, start=1):
            print("  %2d. %-20s %8.2f" % (i, base, aff))


if __name__ == "__main__":
    main()
