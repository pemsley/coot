# coot_commands/commands/navigation.py
#
# Copyright 2026 by Medical Research Council
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

from coot_commands.registry import command
from coot_commands.types import resolve_model, as_int

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
         notes="Centres on the CA of the given chain/residue. With no model "
               "number, uses the active model. Chain and residue may be "
               "separated by a space or a slash.")
def go_to_residue(chain, resno, model=None):
    """Centre the view on a chain/residue."""
    imol = resolve_model(model)
    res = as_int(resno, "residue number")
    if coot is not None:
        coot.set_go_to_atom_molecule(imol)
        coot.set_go_to_atom_chain_residue_atom_name(chain.upper(), res, " CA ")
    return f"Centred on {chain.upper()}/{res} of model {imol}"


@command(r"centre at (?P<x>-?[\d.]+) (?P<y>-?[\d.]+) (?P<z>-?[\d.]+)",
         examples=["centre at 12.0 4.5 -3.2"],
         category=CATEGORY,
         notes="Centres the view on the given orthogonal coordinates.")
def centre_at_xyz(x, y, z):
    """Centre the view on an x, y, z position."""
    from coot_commands.types import as_float
    xf, yf, zf = (as_float(v, "coordinate") for v in (x, y, z))
    if coot is not None:
        coot.set_rotation_centre(xf, yf, zf)
    return f"Centred at ({xf:g}, {yf:g}, {zf:g})"
