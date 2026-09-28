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

from coot_commands.ansi import swatch
from coot_commands.registry import command
from coot_commands.types import (resolve_model, molecule_name, ArgType,
                                 CommandError, ACTIVE_MODEL_NOTE,
                                 loaded_models, loaded_maps)

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
        lines.append(f"  {imol}: {_swatch(imol)}{_name(imol)}{suffix}")
    return "\n".join(lines)


def _swatch(imol: int) -> str:
    """A coloured swatch (plus trailing space) for a map, or '' if unknown."""
    try:
        colours = coot.get_map_colour_py(imol)
    except Exception:
        colours = None
    # get_map_colour_py returns [[r,g,b], [r,g,b]] (the map colour and the
    # negative-level colour) as 0-1 floats, or False for a non-map molecule.
    if not colours:
        return ""
    return swatch(colours[0]) + " "


def _describe_new(models: list[str], maps: list[str]) -> str:
    """Describe the molecules a load added, e.g. 'model 0, 2Fo-Fc map 1, ...'.

    Only the molecules that appeared are named, so loading the tutorial into a
    session that already had things open does not claim credit for them.
    """
    parts = [f"model {imol} ({_name(int(imol))})" for imol in models]
    for imol in maps:
        try:
            is_diff = coot.map_is_difference_map(int(imol)) == 1
        except Exception:  # noqa: BLE001 - fall back to an unlabelled map
            is_diff = False
        kind = "Fo-Fc difference map" if is_diff else "2Fo-Fc map"
        parts.append(f"{kind} {imol}")
    return ", ".join(parts)


@command(r"load tutorial(?: (?:model(?: and data)?|data))?",
         examples=["load tutorial", "load tutorial model and data"],
         category=CATEGORY,
         notes="Loads the bundled tutorial dataset, the same as File > Open "
               "Tutorial. Both the model AND its data are always loaded - "
               "there is no way to get one without the other, and all the "
               "phrasings ('load tutorial', 'load tutorial model', 'load "
               "tutorial data') do exactly the same thing. You end up with "
               "three molecules: the model (tutorial-modern), a 2Fo-Fc map "
               "and an Fo-Fc difference map, both computed from the bundled "
               "MTZ. The 2Fo-Fc map carries the reflection data and becomes "
               "the refinement map, so refinement and rotamer fitting work "
               "straight away. Auto-updating maps are then switched on, so "
               "both maps recompute as you edit the model; turn that off "
               "again with 'set updating maps off'.")
def load_tutorial(**_: Optional[str]) -> str:
    """Load the tutorial model and data, with updating maps switched on."""
    if coot is None:
        return "Loaded the tutorial model and data"

    # Diff the molecule list around the load so the summary names what this
    # command actually added rather than everything that happens to be open.
    before_models, before_maps = set(loaded_models()), set(loaded_maps())
    coot.load_tutorial_model_and_data()
    new_models = sorted(set(loaded_models()) - before_models, key=int)
    new_maps = sorted(set(loaded_maps()) - before_maps, key=int)

    loaded = _describe_new(new_models, new_maps)
    summary = f"Loaded the tutorial: {loaded}" if loaded else \
        "Loaded the tutorial model and data"

    # Updating maps is what makes the tutorial behave the way the tutorial
    # describes - the density follows the atoms as you edit. It needs the
    # refinement map and a difference map, which the load has just provided;
    # if something is nonetheless missing (another dataset already open makes
    # the refinement map ambiguous), the data is still loaded, so report the
    # reason rather than failing the whole command.
    from coot_commands.commands.settings import start_updating_maps
    try:
        return f"{summary}. {start_updating_maps()}"
    except CommandError as e:
        return (f"{summary}. Could not turn on updating maps ({e}) - turn "
                f"them on yourself with 'set updating maps on'")


@command(r"(?:open|display|view) sequence(?: (?:of |for )?model (?P<model>\S+))?",
         examples=["open sequence", "open sequence of model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Opens the sequence view for the model. Close it again with "
               "'close sequence'. " + ACTIVE_MODEL_NOTE)
def open_sequence(model: Optional[str] = None) -> str:
    """Open the sequence view for a model."""
    imol = resolve_model(model)
    if coot is not None:
        coot.sequence_view(imol)
    return f"Opened the sequence of model {imol}"


@command(r"close sequence(?: (?:of |for )?model (?P<model>\S+))?",
         examples=["close sequence", "close sequence of model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Closes the sequence view opened by 'open sequence'. "
               + ACTIVE_MODEL_NOTE)
def close_sequence(model: Optional[str] = None) -> str:
    """Close the sequence view for a model."""
    imol = resolve_model(model)
    if coot is not None:
        coot.remove_sequence_view_from_sequence_view_box(imol)
    return f"Closed the sequence of model {imol}"
