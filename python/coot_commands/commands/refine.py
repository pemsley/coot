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
from coot_commands.types import (RES_SPEC, resolve_model, resolve_residue,
                                 as_int, ArgType, CommandError,
                                 ACTIVE_MODEL_NOTE, ACTIVE_RESIDUE_NOTE)

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


def _chain_residue_range(imol: int, chain_id: str) -> tuple[int, int]:
    """First and last residue numbers of a chain, from its serial numbers."""
    n = coot.chain_n_residues(chain_id, imol)
    if n <= 0:
        raise CommandError(f"no chain {chain_id} in model {imol}")
    first = coot.seqnum_from_serial_number(imol, chain_id, 0)
    last = coot.seqnum_from_serial_number(imol, chain_id, n - 1)
    return first, last

def _refine_sphere(imol: int, chain_id: str, res: int, ins: str = "") -> str:
    """Real-space refine the sphere of atoms around a residue.

    ``rsr_sphere_refine`` acts on the *active* atom rather than taking a
    residue argument, so we first set the go-to atom to the target residue
    (which also makes it active), then refine, using the same immediate-
    replacement save/restore dance as :func:`_refine_zone`.
    """
    if coot is None:
        return f"Refined the sphere around {chain_id}/{res} of model {imol}"
    if coot.imol_refinement_map() < 0:
        raise CommandError("no map set for refinement - open a map first")
    coot.set_go_to_atom_molecule(imol)
    if coot.set_go_to_atom_from_res_spec_py([chain_id, res, ins]) <= 0:
        return f"No residue {chain_id}/{res} in model {imol}"
    replacement_state = coot.refinement_immediate_replacement_state()
    coot.set_refinement_immediate_replacement(1)
    try:
        coot.rsr_sphere_refine()
        coot.accept_regularizement()
    finally:
        if replacement_state == 0:
            coot.set_refinement_immediate_replacement(0)
    return f"Refined the sphere around {chain_id}/{res} of model {imol}"

@command(r"(?:real[- ]space )?refine "
         r"(?:model (?P<model>\S+) )?"
         r"(?:residues? )?(?P<chain>[A-Za-z0-9])[ /](?P<res1>-?\d+)"
         r"(?: to |\s*-\s*)(?P<res2>-?\d+)",
         examples=["refine A 45 to 50", "refine A/45-50",
                   "real space refine A 100 to 105"],
         category=CATEGORY,
         notes="Real-space refines the given residue range against the "
               "refinement map. " + ACTIVE_MODEL_NOTE)
def refine_range(chain: str, res1: str, res2: str,
                 model: Optional[str] = None) -> str:
    """Real-space refine a range of residues."""
    imol = resolve_model(model)
    return _refine_zone(imol, chain.upper(),
                        as_int(res1, "residue number"),
                        as_int(res2, "residue number"))


@command(r"(?:real[- ]space )?refine "
         r"(?:model (?P<model>\S+) )?"
         r"chain (?P<chain>[A-Za-z0-9])",
         examples=["refine chain A", "refine chain B"],
         category=CATEGORY,
         notes="Real-space refines a whole chain (its full residue range) "
               "against the refinement map. " + ACTIVE_MODEL_NOTE)
def refine_chain(chain: str, model: Optional[str] = None) -> str:
    """Real-space refine an entire chain."""
    imol = resolve_model(model)
    chain_id = chain.upper()
    if coot is None:
        return f"Refined chain {chain_id} of model {imol}"
    first, last = _chain_residue_range(imol, chain_id)
    _refine_zone(imol, chain_id, first, last)
    return f"Refined chain {chain_id} ({first}-{last}) of model {imol}"


@command(r"(?:real[- ]space )?refine "
         r"(?:model (?P<model>\S+) )?"
         r"(?:residue |res )?(?P<chain>[A-Za-z0-9])[ /](?P<resno>-?\d+)",
         examples=["refine A 45", "refine A/45", "refine residue B 12"],
         category=CATEGORY,
         notes="Real-space refines a single residue against the refinement "
               "map. " + ACTIVE_MODEL_NOTE)
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
         notes="Real-space refines the active residue (the one at the centre "
               "of the screen) against the refinement map.")
def refine_active(**_: Optional[str]) -> str:
    """Real-space refine the active residue."""
    if coot is None:
        return "There has been an issue"
    active = coot.active_residue_py()
    if not active:
        return "No active residue - centre on a model first"
    imol, chain, resno = active[0], active[1], active[2]
    return _refine_zone(imol, chain, resno, resno)


@command(r"(?:real[- ]space )?refine sphere"
         r"(?: model (?P<model>\S+))?"
         r"(?: (?:residue |res )?" + RES_SPEC + r")?",
         examples=["refine sphere A/89", "refine sphere A 89", "refine sphere"],
         category=CATEGORY,
         notes="Real-space refines the sphere of atoms around the named "
               "residue against the refinement map. " + ACTIVE_RESIDUE_NOTE)
def refine_sphere(chain: Optional[str] = None, resno: Optional[str] = None,
                  model: Optional[str] = None) -> str:
    """Real-space refine the sphere around a residue."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    return _refine_sphere(imol, chain_id, res, ins)


@command(r"(?:shiftfield )?refine (?:the )?"
         r"(?:b[- ]?factors?|b[- ]?values?|adps?|temperature factors?|"
         r"atomic displacement parameters?)"
         r"(?: (?:of |in |for )?model (?P<model>\S+))?",
         examples=["refine b factors", "refine b-factors of model 0",
                   "refine adps"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Refines the atomic B-factors (ADPs) of the whole model using "
               "the shiftfield method, against the reflection data. Needs a "
               "map with reflection data set as the refinement map (unlike "
               "real-space refinement, this is a whole-molecule reciprocal-"
               "space step). " + ACTIVE_MODEL_NOTE)
def refine_b_factors(model: Optional[str] = None) -> str:
    """Refine the atomic B-factors (ADPs) of a model."""
    imol = resolve_model(model)
    if coot is None:
        return f"Refined B-factors of model {imol}"
    if coot.imol_refinement_map() < 0:
        raise CommandError(
            "no map set for refinement - open a map with reflection data first")
    coot.shiftfield_b_factor_refinement(imol)
    return f"Refined B-factors of model {imol}"
