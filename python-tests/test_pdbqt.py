# test_pdbqt.py
# Copyright 2026 by Global Phasing Ltd.
# Author: Paul Emsley
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

"""Tests for PDBQT (AutoDock/Vina) export via molecules_container_t:

  - export_molecule_as_pdbqt(imol, file_name)      rigid receptor
  - export_ligand_as_pdbqt(imol, cid, file_name)   flexible ligand + torsion tree

The tests validate that the written files are well-formed PDBQT: correct
fixed-column layout, valid AutoDock atom types, parseable partial charges,
non-polar hydrogens merged (receptor), and a balanced ROOT/BRANCH/ENDBRANCH/
TORSDOF torsion tree (ligand).
"""

import os
import sys

try:
    import pytest
except ImportError:
    pytest = None


# ---------------------------------------------------------------------------
# locate and import the headless-api module (skip cleanly if not available)
# ---------------------------------------------------------------------------

def _import_chapi():
    # Prefer an explicit build directory over any installed copy: a stale
    # system/Homebrew coot_headless_api can shadow a freshly built one.
    here = os.path.dirname(os.path.abspath(__file__))
    dirs = []
    env = os.environ.get("COOT_HEADLESS_API_DIR")
    if env:
        dirs.append(env)
    for rel in ("../../build-chapi-with-homebrew-deps", "../../build", "../build", "build"):
        dirs.append(os.path.normpath(os.path.join(here, rel)))
    # insert in reverse so the first existing candidate ends up at the front of sys.path
    for d in reversed(dirs):
        if os.path.isdir(d):
            sys.path.insert(0, d)
    try:
        import coot_headless_api
        return coot_headless_api
    except ImportError:
        return None


def _find_data_file():
    here = os.path.dirname(os.path.abspath(__file__))
    for rel in ("../data/tutorial-modern.pdb", "data/tutorial-modern.pdb",
                "../api/tutorial-modern.pdb"):
        p = os.path.normpath(os.path.join(here, rel))
        if os.path.exists(p):
            return p
    env = os.environ.get("COOT_TEST_DATA_DIR")
    if env:
        p = os.path.join(env, "tutorial-modern.pdb")
        if os.path.exists(p):
            return p
    return None


chapi = _import_chapi()
DATA = _find_data_file()

if pytest is not None:
    pytestmark = pytest.mark.skipif(chapi is None or DATA is None,
                                    reason="coot_headless_api or tutorial-modern.pdb not available")

# The AutoDock (AD4) atom types this exporter can emit.
AD_TYPES = {"C", "A", "N", "NA", "OA", "S", "SA", "HD", "H", "P",
            "F", "Cl", "Br", "I", "B", "Mg", "Ca", "Mn", "Fe", "Cu", "Zn", "Se"}


# ---------------------------------------------------------------------------
# structural validators
# ---------------------------------------------------------------------------

def _parse_atoms(lines):
    """Return list of (serial, ad_type) and assert per-atom PDBQT column layout."""
    atoms = []
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")):
            assert len(ln) >= 79, "PDBQT atom line shorter than 79 columns: %r" % ln
            serial = int(ln[6:11])
            ad_type = ln[77:79].strip()
            assert ad_type in AD_TYPES, "unexpected AutoDock type %r in %r" % (ad_type, ln)
            float(ln[70:76])  # partial charge must parse (raises on failure)
            atoms.append((serial, ad_type))
    return atoms


def _assert_contiguous_serials(atoms):
    for i, (serial, _) in enumerate(atoms):
        assert serial == i + 1, "atom serials not contiguous at %d (got %d)" % (i + 1, serial)


def _assert_torsion_tree(lines):
    """Assert a well-formed torsion tree and return TORSDOF."""
    stack = []
    nbranch = 0
    torsdof = None
    have_root = have_endroot = False
    defined = set()
    for ln in lines:
        rec = ln.split()[0] if ln.split() else ""
        if ln.startswith(("ATOM", "HETATM")):
            defined.add(int(ln[6:11]))
        elif rec == "ROOT":
            have_root = True
        elif rec == "ENDROOT":
            have_endroot = True
        elif rec == "BRANCH":
            a, b = int(ln.split()[1]), int(ln.split()[2])
            assert a in defined, "BRANCH parent atom %d not yet defined" % a
            stack.append((a, b))
            nbranch += 1
        elif rec == "ENDBRANCH":
            a, b = int(ln.split()[1]), int(ln.split()[2])
            assert stack, "ENDBRANCH without matching BRANCH"
            assert stack.pop() == (a, b), "ENDBRANCH %d %d does not match its BRANCH" % (a, b)
        elif rec == "TORSDOF":
            torsdof = int(ln.split()[1])
    assert have_root and have_endroot, "missing ROOT/ENDROOT"
    assert not stack, "unclosed BRANCH records"
    assert torsdof is not None, "missing TORSDOF"
    assert torsdof == nbranch, "TORSDOF %d != number of BRANCH records %d" % (torsdof, nbranch)
    return torsdof


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

def _tmp(name):
    import tempfile
    return os.path.join(tempfile.gettempdir(), name)


def test_receptor_pdbqt():
    """Rigid receptor: valid columns/types, non-polar H merged, no torsion tree."""
    mc = chapi.molecules_container_t(True)
    mc.set_use_gemmi(False)
    imol = mc.read_pdb(DATA)
    assert imol >= 0
    mc.add_hydrogen_atoms(imol)              # Coot's own reduce; export needs H present

    out = _tmp("test_pdbqt_receptor.pdbqt")
    n = mc.export_molecule_as_pdbqt(imol, out)
    assert n > 0, "no atoms written"

    lines = open(out).read().splitlines()
    atoms = _parse_atoms(lines)
    assert len(atoms) == n
    _assert_contiguous_serials(atoms)

    # rigid receptor: no torsion-tree records
    assert not any(ln.startswith(("ROOT", "BRANCH", "TORSDOF")) for ln in lines)

    types = set(t for _, t in atoms)
    assert "HD" in types, "expected polar hydrogens (HD) in the receptor"
    assert "H" not in types, "non-polar hydrogens should be merged, not written as H"
    os.remove(out)


def test_ligand_pdbqt():
    """Flexible ligand: valid file plus a balanced torsion tree with TORSDOF."""
    mc = chapi.molecules_container_t(True)
    mc.set_use_gemmi(False)
    imol = mc.get_monomer("TYR")             # a residue with a rotatable side chain
    assert imol >= 0

    out = _tmp("test_pdbqt_ligand.pdbqt")
    n = mc.export_ligand_as_pdbqt(imol, "/*/*/*", out)
    assert n > 0, "no atoms written"

    lines = open(out).read().splitlines()
    atoms = _parse_atoms(lines)
    assert len(atoms) == n
    _assert_contiguous_serials(atoms)

    torsdof = _assert_torsion_tree(lines)
    assert torsdof >= 1, "expected at least one rotatable bond for TYR"

    types = set(t for _, t in atoms)
    assert "A" in types, "expected aromatic ring carbons (A) for TYR"
    os.remove(out)


# ---------------------------------------------------------------------------
# allow running directly: python test_pdbqt.py
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if chapi is None:
        sys.exit("coot_headless_api not importable (set PYTHONPATH or COOT_HEADLESS_API_DIR)")
    if DATA is None:
        sys.exit("tutorial-modern.pdb not found (set COOT_TEST_DATA_DIR)")
    failures = 0
    for t in (test_receptor_pdbqt, test_ligand_pdbqt):
        try:
            t()
            print("PASS:", t.__name__)
        except AssertionError as e:
            failures += 1
            print("FAIL:", t.__name__, "-", e)
    sys.exit(1 if failures else 0)
