# coot_commands/context_tools.py
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

"""Custom context tools for the assistant: live queries, not actions.

These are :func:`~coot_commands.tools.custom_tool` handlers - always offered to
the model, unlike the retrieval-filtered commands - that let it read Coot's
current state so it can resolve what the user *means* before acting.  The
motivating case: the user says "refine here" or "what's this residue?"; the
model calls :func:`get_active_residue` to learn the residue at the centre of
the screen, then acts on it.

Importing this module registers the tools (it is imported by the package
``__init__``).  Add a query here and it becomes available to the assistant with
no other wiring.  They run wherever ``execute_tool`` runs - in-process, or
inside a live Coot over the socket - so the values are always current.
"""

from __future__ import annotations

from coot_commands.tools import custom_tool

try:
    import coot
except ImportError:
    coot = None


@custom_tool(
    "get_active_residue",
    "Return the residue currently at the centre of the screen (the \"active\" "
    "residue) as its chain, residue number and model number. Call this when the "
    "user refers to \"here\", \"this residue\", \"the current residue\" or "
    "similar, to find out which residue and model they mean before acting.")
def get_active_residue() -> str:
    """Report the residue at the centre of the screen."""
    if coot is None:
        return "the Coot API is not available"
    active = coot.active_residue_py()
    if not active:
        return "No active residue - centre on a model first"
    imol, chain, resno = active[0], active[1], active[2]
    ins_code = active[3] if len(active) > 3 else ""
    spec = f"{chain}/{resno}" + (f" (insertion code '{ins_code}')" if ins_code else "")
    return f"Active residue: {spec} of model {imol}"


@custom_tool(
    "get_active_map",
    "Return the map molecule number currently set for refinement (the map that "
    "refinement commands use). Call this to find out which map is active.")
def get_active_map() -> str:
    """Report the molecule number of the map set for refinement."""
    if coot is None:
        return "the Coot API is not available"
    imol = coot.imol_refinement_map()
    if imol is None or imol < 0:
        return "No map is set for refinement - open a map first"
    return f"Refinement map: molecule {imol}"
