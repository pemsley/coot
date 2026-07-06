# coot_commands/commands/navigation.py
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

"""Commands that move the view centre around the structure."""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import resolve_model, as_int, ArgType

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Navigation"


@command(r"(?:go to|centre on|center on|goto) "
         r"(?:model (?P<model>\S+) )?"
         r"(?:residue |res )?(?P<chain>[A-Za-z0-9])[ /](?P<resno>-?\d+)",
         examples=["go to A 45", "centre on A/45", "go to model 0 B 12"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Centres on the given chain/residue, picking a sensible atom "
               "(CA for protein, P/C1' for nucleotides, and so on) so it "
               "works for waters, ligands and nucleic acids too. With no "
               "model number, uses the active model. Chain and residue may "
               "be separated by a space or a slash.")
def go_to_residue(chain: str, resno: str, model: Optional[str] = None) -> str:
    """Centre the view on a chain/residue."""
    imol = resolve_model(model)
    res = as_int(resno, "residue number")
    chain_id = chain.upper()
    if coot is not None:
        coot.set_go_to_atom_molecule(imol)
        success = coot.set_go_to_atom_from_res_spec_py([chain_id, res, ""])
        if success <= 0:
            return f"No residue {chain_id}/{res} in model {imol}"
    return f"Centred on {chain_id}/{res} of model {imol}"


def _step_residue(direction: str) -> str:
    """Centre on the residue next to the active one, in *direction*.

    ``direction`` is "next" or "previous".  We ask Coot for the active
    atom, let ``goto_(next|prev)_atom_maybe_py`` work out the neighbouring
    residue's intelligent atom (the same logic the space bar uses), then
    centre on it - so stepping skips to the next residue, not the next
    atom within the current one.
    """
    if coot is None:
        return f"Stepped to the {direction} residue"
    active = coot.active_residue_py()
    if not active:
        return "No active residue - click a model first"
    imol, chain, resno, ins_code, atom_name = active[:5]
    coot.set_go_to_atom_molecule(imol)
    finder = (coot.goto_next_atom_maybe_py if direction == "next"
              else coot.goto_prev_atom_maybe_py)
    nxt = finder(chain, resno, ins_code, atom_name)
    if not nxt:
        return f"No {direction} residue"
    n_chain, n_resno, _n_ins, n_atom = nxt
    coot.set_go_to_atom_chain_residue_atom_name(n_chain, n_resno, n_atom)
    return f"Centred on {n_chain}/{n_resno} of model {imol}"


@command(r"(?:go to |goto )?(?:next|forward)(?: residue| res)?",
         examples=["next residue", "next", "forward residue"],
         category=CATEGORY,
         notes="Centres on the next residue after the active one (like the "
               "space bar). Click a model first to set the active residue.")
def next_residue(**_: Optional[str]) -> str:
    """Centre on the next residue."""
    return _step_residue("next")


@command(r"(?:go to |goto )?(?:previous|prev|back)(?: residue| res)?",
         examples=["previous residue", "prev", "back residue"],
         category=CATEGORY,
         notes="Centres on the residue before the active one (like "
               "shift-space). Click a model first to set the active residue.")
def previous_residue(**_: Optional[str]) -> str:
    """Centre on the previous residue."""
    return _step_residue("previous")


@command(r"centre at (?P<x>-?[\d.]+) (?P<y>-?[\d.]+) (?P<z>-?[\d.]+)",
         examples=["centre at 12.0 4.5 -3.2"],
         category=CATEGORY,
         notes="Centres the view on the given orthogonal coordinates.")
def centre_at_xyz(x: str, y: str, z: str) -> str:
    """Centre the view on an x, y, z position."""
    from coot_commands.types import as_float
    xf, yf, zf = (as_float(v, "coordinate") for v in (x, y, z))
    if coot is not None:
        coot.set_rotation_centre(xf, yf, zf)
    return f"Centred at ({xf:g}, {yf:g}, {zf:g})"
