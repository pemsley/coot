# coot_commands/commands/refine.py
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

"""Commands that run real-space refinement against the refinement map."""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import resolve_model, as_int, CommandError

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Refinement"


def _refine_zone(imol: int, chain_id: str, res1: int, res2: int) -> str:
    """Real-space refine a residue range, accepting the result immediately.

    Mirrors Coot's own ``refine_active_residue`` helper: a refinement map
    must be set, and we flip on immediate replacement so the refinement is
    applied without the interactive accept/reject dialog, restoring the
    previous state afterwards.
    """
    if coot is None:
        return f"Refined {chain_id}/{res1}-{res2} of model {imol}"
    if coot.imol_refinement_map() < 0:
        raise CommandError("no map set for refinement - open a map first")
    lo, hi = (res1, res2) if res1 <= res2 else (res2, res1)
    replacement_state = coot.refinement_immediate_replacement_state()
    coot.set_refinement_immediate_replacement(1)
    try:
        coot.refine_zone(imol, chain_id, lo, hi, "")
        coot.accept_regularizement()
    finally:
        if replacement_state == 0:
            coot.set_refinement_immediate_replacement(0)
    if lo == hi:
        return f"Refined {chain_id}/{lo} of model {imol}"
    return f"Refined {chain_id}/{lo}-{hi} of model {imol}"


@command(r"(?:real[- ]space )?refine "
         r"(?:model (?P<model>\S+) )?"
         r"(?:residues? )?(?P<chain>[A-Za-z0-9])[ /](?P<res1>-?\d+)"
         r"(?: to |\s*-\s*)(?P<res2>-?\d+)",
         examples=["refine A 45 to 50", "refine A/45-50",
                   "real space refine A 100 to 105"],
         category=CATEGORY,
         notes="Real-space refines the given residue range against the "
               "refinement map. With no model number, uses the active model.")
def refine_range(chain: str, res1: str, res2: str,
                 model: Optional[str] = None) -> str:
    """Real-space refine a range of residues."""
    imol = resolve_model(model)
    return _refine_zone(imol, chain.upper(),
                        as_int(res1, "residue number"),
                        as_int(res2, "residue number"))


@command(r"(?:real[- ]space )?refine "
         r"(?:model (?P<model>\S+) )?"
         r"(?:residue |res )?(?P<chain>[A-Za-z0-9])[ /](?P<resno>-?\d+)",
         examples=["refine A 45", "refine A/45", "refine residue B 12"],
         category=CATEGORY,
         notes="Real-space refines a single residue against the refinement "
               "map. With no model number, uses the active model.")
def refine_residue(chain: str, resno: str,
                   model: Optional[str] = None) -> str:
    """Real-space refine a single residue."""
    imol = resolve_model(model)
    res = as_int(resno, "residue number")
    return _refine_zone(imol, chain.upper(), res, res)


@command(r"(?:real[- ]space )?refine(?: (?:the )?active(?: residue)?"
         r"| this(?: residue)?| here)?$",
         examples=["refine", "refine active residue", "refine here"],
         category=CATEGORY,
         notes="Real-space refines the active residue (the last one clicked) "
               "against the refinement map.")
def refine_active(**_: Optional[str]) -> str:
    """Real-space refine the active residue."""
    if coot is None:
        return "Refined the active residue"
    active = coot.active_residue_py()
    if not active:
        return "No active residue - click a model first"
    imol, chain, resno = active[0], active[1], active[2]
    return _refine_zone(imol, chain, resno, resno)
