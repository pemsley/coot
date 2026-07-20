# coot_commands/commands/display.py
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

"""Commands that show and hide models and maps."""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (RES_SPEC, resolve_model, resolve_map,
                                 resolve_residue, ArgType,
                                 ACTIVE_MODEL_NOTE, ACTIVE_MAP_NOTE)

try:
    import coot
except ImportError:
    coot = None  


CATEGORY = "Display"


@command(r"(?:show|display) model(?: (?P<model>\S+))?",
         examples=["show model 0", "display model"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes=ACTIVE_MODEL_NOTE)
def show_model(model: Optional[str] = None) -> str:
    """Show (display) a model."""
    imol = resolve_model(model)
    if coot is not None:
        coot.set_mol_displayed(imol, 1)
    return f"Showing model {imol}"


@command(r"(?:hide|undisplay) model(?: (?P<model>\S+))?",
         examples=["hide model 0", "undisplay model"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes=ACTIVE_MODEL_NOTE)
def hide_model(model: Optional[str] = None) -> str:
    """Hide (undisplay) a model."""
    imol = resolve_model(model)
    if coot is not None:
        coot.set_mol_displayed(imol, 0)
    return f"Hiding model {imol}"


@command(r"(?:show|display) map(?: (?P<map>\S+))?",
         examples=["show map 1", "display map"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP},
         notes=ACTIVE_MAP_NOTE)
def show_map(map: Optional[str] = None) -> str:
    """Show (display) a map."""
    imol = resolve_map(map)
    if coot is not None:
        coot.set_map_displayed(imol, 1)
    return f"Showing map {imol}"


@command(r"(?:hide|undisplay) map(?: (?P<map>\S+))?",
         examples=["hide map 1", "undisplay map"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP},
         notes=ACTIVE_MAP_NOTE)
def hide_map(map: Optional[str] = None) -> str:
    """Hide (undisplay) a map."""
    imol = resolve_map(map)
    if coot is not None:
        coot.set_map_displayed(imol, 0)
    return f"Hiding map {imol}"


@command(r"(?:show|display) only(?: (?:the )?)?active",
         examples=["show only active"],
         category=CATEGORY,
         notes="Hides every model except the active one.")
def display_only_active(**_: Optional[str]) -> str:
    """Show only the active model, hiding the others."""
    if coot is not None:
        coot.display_only_active()
    return "Showing only the active model"


@command(r"(?:show|display) (?:residue )?environment(?: " + RES_SPEC + r")?",
         examples=["show residue environment", "show environment A/89"],
         category=CATEGORY,
         notes="Shows the environment distances (contacts and H-bonds) around "
               "a residue, centring on it. With no residue named, uses the "
               "active residue (the one at the centre of the screen).")
def show_environment(chain: Optional[str] = None,
                     resno: Optional[str] = None) -> str:
    """Show the environment distances around a residue."""
    if coot is None:
        return "Showing the residue environment"
    imol, chain_id, res, ins = resolve_residue(chain, resno)
    coot.set_go_to_atom_molecule(imol)
    coot.set_go_to_atom_from_res_spec_py([chain_id, res, ins])
    coot.set_show_environment_distances(1)
    return f"Showing the environment of {chain_id}/{res} of model {imol}"


@command(r"(?:hide|undisplay) (?:residue )?environment",
         examples=["hide residue environment", "hide environment"],
         category=CATEGORY,
         notes="Hides the residue environment distances.")
def hide_environment(**_: Optional[str]) -> str:
    """Hide the residue environment distances."""
    if coot is not None:
        coot.set_show_environment_distances(0)
    return "Hiding the residue environment"