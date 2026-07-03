# coot_commands/commands/display.py
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

"""Commands that show and hide models and maps."""

from coot_commands.registry import command
from coot_commands.types import resolve_model, resolve_map

try:
    import coot
except ImportError:
    coot = None  # allows this module to be imported/tested without Coot


CATEGORY = "Display"


@command(r"(?:show|display) model(?: (?P<model>\S+))?",
         examples=["show model 0", "display model"],
         category=CATEGORY,
         notes="With no number, acts on the active model.")
def show_model(model=None):
    """Show (display) a model."""
    imol = resolve_model(model)
    if coot is not None:
        coot.set_mol_displayed(imol, 1)
    return f"Showing model {imol}"


@command(r"(?:hide|undisplay) model(?: (?P<model>\S+))?",
         examples=["hide model 0", "undisplay model"],
         category=CATEGORY,
         notes="With no number, acts on the active model.")
def hide_model(model=None):
    """Hide (undisplay) a model."""
    imol = resolve_model(model)
    if coot is not None:
        coot.set_mol_displayed(imol, 0)
    return f"Hiding model {imol}"


@command(r"(?:show|display) map(?: (?P<map>\S+))?",
         examples=["show map 1", "display map"],
         category=CATEGORY,
         notes="With no number, acts on the active (refinement) map.")
def show_map(map=None):
    """Show (display) a map."""
    imol = resolve_map(map)
    if coot is not None:
        coot.set_map_displayed(imol, 1)
    return f"Showing map {imol}"


@command(r"(?:hide|undisplay) map(?: (?P<map>\S+))?",
         examples=["hide map 1", "undisplay map"],
         category=CATEGORY,
         notes="With no number, acts on the active (refinement) map.")
def hide_map(map=None):
    """Hide (undisplay) a map."""
    imol = resolve_map(map)
    if coot is not None:
        coot.set_map_displayed(imol, 0)
    return f"Hiding map {imol}"


@command(r"(?:show|display) only(?: (?:the )?)?active",
         examples=["show only active"],
         category=CATEGORY,
         notes="Hides every model except the active one.")
def display_only_active(**_):
    """Show only the active model, hiding the others."""
    if coot is not None:
        coot.display_only_active()
    return "Showing only the active model"
