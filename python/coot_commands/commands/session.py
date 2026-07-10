# coot_commands/commands/session.py
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

"""Session-level commands: list molecules, load the tutorial, view sequences.

Whole-model operations (merge, superpose) live in
:mod:`coot_commands.commands.models`.
"""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (resolve_model, molecule_name, ArgType,
                                 ACTIVE_MODEL_NOTE, loaded_models, loaded_maps)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Session"


def _name(imol: int) -> str:
    """A short, path-less molecule name, falling back to 'molecule N'."""
    return molecule_name(imol) or f"molecule {imol}"


@command(r"list (?:all )?(?:models|molecules)",
         examples=["list models", "list molecules"],
         category=CATEGORY,
         notes="Lists the loaded models with their molecule number and name.")
def list_models(**_: Optional[str]) -> str:
    """List the loaded models."""
    if coot is None:
        return "No models loaded"
    imols = [int(i) for i in loaded_models()]
    if not imols:
        return "No models loaded"
    lines = [f"{len(imols)} model(s):"]
    lines += [f"  {imol}: {_name(imol)}" for imol in imols]
    return "\n".join(lines)


@command(r"list (?:all )?maps",
         examples=["list maps"],
         category=CATEGORY,
         notes="Lists the loaded maps with their molecule number and name, "
               "marking which are difference maps.")
def list_maps(**_: Optional[str]) -> str:
    """List the loaded maps."""
    if coot is None:
        return "No maps loaded"
    imols = [int(i) for i in loaded_maps()]
    if not imols:
        return "No maps loaded"
    lines = [f"{len(imols)} map(s):"]
    for imol in imols:
        try:
            is_diff = coot.map_is_difference_map(imol) == 1
        except Exception:
            is_diff = False
        suffix = " (difference map)" if is_diff else ""
        lines.append(f"  {imol}: {_name(imol)}{suffix}")
    return "\n".join(lines)


@command(r"load tutorial(?: (?:model(?: and data)?|data))?",
         examples=["load tutorial", "load tutorial model and data"],
         category=CATEGORY,
         notes="Loads the bundled tutorial model and its data (map "
               "coefficients), the same as File > Open Tutorial.")
def load_tutorial(**_: Optional[str]) -> str:
    """Load the tutorial model and data."""
    if coot is not None:
        coot.load_tutorial_model_and_data()
    return "Loaded the tutorial model and data"


@command(r"(?:show|display|view) sequence(?: (?:of |for )?model (?P<model>\S+))?",
         examples=["show sequence", "show sequence of model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Opens the sequence view for the model. " + ACTIVE_MODEL_NOTE)
def show_sequence(model: Optional[str] = None) -> str:
    """Show the sequence view for a model."""
    imol = resolve_model(model)
    if coot is not None:
        coot.sequence_view(imol)
    return f"Showing the sequence of model {imol}"
