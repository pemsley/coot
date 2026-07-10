# coot_commands/commands/validation.py
#
# Copyright 2026 Jordan Dialpuri, Medical Research Council Laboratory of Molecular Biology
#
# This file is part of Coot
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

"""Validation commands that summarise model outliers as text.

Rather than opening Coot's (C++/GTK) validation graphs - which are not
scriptable from the embedded interpreter - these commands call the
low-level ``coot.*`` scoring functions directly (the same ones the
quick-test validation uses) and report a short text summary.
"""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (resolve_model, resolve_map, as_float, as_int,
                                 first_present, ArgType, CommandError,
                                 ACTIVE_MODEL_NOTE)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Validation"

# Coordinates of the blobs found by the most recent "find blobs" command,
# largest first, so a follow-up "go to blob N" can recentre on one.  Each
# entry is an (x, y, z) tuple; index 0 is blob 1.
_last_blobs: list[tuple[float, float, float]] = []


def _rama_improbables(imol: int) -> list:
    """Residues with Ramachandran probability below 0.02."""
    try:
        rs = coot.all_molecule_ramachandran_score_py(imol)
        scored = rs[5]
        return [item for item in scored if item[2] < 0.02]
    except Exception:
        return []


def _atom_overlaps(imol: int) -> list:
    """Atom overlaps with a clash volume above 2.0 A^3."""
    try:
        overlaps = coot.molecule_atom_overlaps_py(imol)
        return [o for o in overlaps
                if isinstance(o, dict) and o.get("overlap-volume", 0) > 2.0]
    except Exception:
        return []


def _as_list(value) -> list:
    """Coerce a maybe-list binding result to a list."""
    return list(value) if isinstance(value, (list, tuple)) else []


@command(r"validate anomalies(?: (?:of |in |for )?model (?P<model>\S+))?|"
         r"(?:validate|check|find) (?:outliers|anomalies)"
         r"(?: (?:of |in |for )?model (?P<model2>\S+))?",
         examples=["validate anomalies", "find outliers"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "model2": ArgType.MODEL},
         notes="Summarises geometry outliers - Ramachandran improbables, "
               "atom overlaps (clashes), C-beta deviations and chiral volume "
               "errors - as text. " + ACTIVE_MODEL_NOTE)
def validate_anomalies(model: Optional[str] = None,
                       model2: Optional[str] = None) -> str:
    """Summarise model outliers (Ramachandran, clashes, C-beta, chirals)."""
    imol = resolve_model(first_present(model, model2))
    if coot is None:
        return f"Validated model {imol}"

    rama = _rama_improbables(imol)
    overlaps = _atom_overlaps(imol)
    c_beta = _as_list(coot.c_beta_deviations_py(imol)) \
        if hasattr(coot, "c_beta_deviations_py") else []
    chirals = _as_list(coot.chiral_volume_errors_py(imol)) \
        if hasattr(coot, "chiral_volume_errors_py") else []

    total = len(rama) + len(overlaps) + len(c_beta) + len(chirals)
    if total == 0:
        return f"No outliers found in model {imol}"

    parts = [f"Model {imol}: {total} outlier(s) -",
             f"{len(rama)} Ramachandran, {len(overlaps)} clash(es), "
             f"{len(c_beta)} C-beta, {len(chirals)} chiral."]
    if rama:
        worst = sorted(rama, key=lambda item: item[2])[:3]
        worst_str = ", ".join(f"{item[1]} ({item[2]:.3f})" for item in worst)
        parts.append(f"Worst Ramachandran: {worst_str}")
    return " ".join(parts)




@command(r"(?:validate|check|find|list) (?:unmodell?ed |unbuilt |missing )?blobs"
         r"(?: (?:of |in |for )?model (?P<model>\S+))?"
         r"(?: (?:against |using |with |in )?map (?P<map>\S+))?"
         r"(?: (?:above |over |at )(?P<sigma>[\d.]+)\s*(?:sigma|rmsd)?)?",
         examples=["find blobs", "check unmodelled blobs",
                   "find blobs above 1.5 sigma"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "map": ArgType.MAP},
         notes="Lists peaks of density not accounted for by the model - "
               "candidate sites for waters, ligands or unbuilt residues. "
               "Uses the refinement map and a 1.0 sigma cut-off by default; "
               "add e.g. 'above 1.5 sigma' to raise the threshold.")
def validate_unmodelled_blobs(model: Optional[str] = None,
                              map: Optional[str] = None,
                              sigma: Optional[str] = None) -> str:
    """Summarise unmodelled density blobs (candidate build sites)."""
    imol = resolve_model(model)
    cut_off = as_float(sigma, "sigma cut-off") if sigma is not None else 1.0
    if coot is None:
        return f"Validated model {imol}"

    imol_map = resolve_map(map)
    blobs = _as_list(coot.find_blobs_py(imol, imol_map, cut_off))
    global _last_blobs
    if not blobs:
        _last_blobs = []
        return (f"No unmodelled blobs found in model {imol} "
                f"above {cut_off:g} sigma")

    # find_blobs_py returns [[x, y, z], volume] pairs; report the largest and
    # remember their positions so "go to blob N" can recentre on one.
    blobs.sort(key=lambda b: b[1], reverse=True)
    _last_blobs = [(b[0][0], b[0][1], b[0][2]) for b in blobs]
    worst = "; ".join(
        f"blob {i} at ({b[0][0]:.1f}, {b[0][1]:.1f}, {b[0][2]:.1f}) "
        f"vol {b[1]:.1f}"
        for i, b in enumerate(blobs[:3], start=1))
    return (f"Model {imol}: {len(blobs)} unmodelled blob(s) above "
            f"{cut_off:g} sigma. Largest: {worst}")


@command(r"(?:go to|centre on|center on|goto) blob (?P<n>\d+)",
         examples=["go to blob 1", "centre on blob 2"],
         category=CATEGORY,
         notes="Centres the view on one of the blobs from the most recent "
               "'find blobs' command, numbered from 1 (largest first). Run "
               "'find blobs' first to populate the list.")
def go_to_blob(n: str) -> str:
    """Centre the view on a blob from the last 'find blobs' result."""
    index = as_int(n, "blob number")
    if not _last_blobs:
        raise CommandError("no blobs to go to - run 'find blobs' first")
    if index < 1 or index > len(_last_blobs):
        raise CommandError(
            f"blob {index} does not exist - there are {len(_last_blobs)} "
            f"(numbered 1 to {len(_last_blobs)})")
    x, y, z = _last_blobs[index - 1]
    if coot is not None:
        coot.set_rotation_centre(x, y, z)
    return f"Centred on blob {index} at ({x:.1f}, {y:.1f}, {z:.1f})"


