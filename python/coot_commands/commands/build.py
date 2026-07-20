# coot_commands/commands/build.py
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

"""Commands that add atoms and monomers to the model at the view centre."""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Building"


@command(r"add water(?:s)?(?: here| at (?:the )?(?:centre|center|pointer))?$",
         examples=["add water here", "add water"],
         category=CATEGORY,
         notes="Places a water atom at the screen centre (the pointer "
               "position), like the 'Place Atom At Pointer' tool.")
def add_water(**_: Optional[str]) -> str:
    """Add a water at the screen centre."""
    if coot is not None:
        coot.place_typed_atom_at_pointer("Water")
    return "Added a water at the centre"


@command(r"add (?:solvent|ligand|monomer) (?P<code>\w+)"
         r"(?: here| at (?:the )?(?:centre|center))?$",
         examples=["add solvent SO4 here", "add ligand GOL"],
         category=CATEGORY,
         notes="Fetches the named monomer from the dictionary and moves it to "
               "the screen centre. Use \"add water\" for a single water atom.")
def add_monomer(code: str) -> str:
    """Add a monomer (by 3-letter code) at the screen centre."""
    comp_id = code.upper()
    if coot is None:
        return f"Added {comp_id} at the centre"
    imol = coot.get_monomer(comp_id)
    if not isinstance(imol, int) or imol < 0:
        return f"Could not find monomer {comp_id} in the dictionary"
    coot.move_molecule_to_screen_centre_internal(imol)
    return f"Added {comp_id} as model {imol} at the centre"
