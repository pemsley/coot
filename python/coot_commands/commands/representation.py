# coot_commands/commands/representation.py
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

"""Commands that change how models are drawn: carbon colours, symmetry."""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (resolve_model, first_present, ArgType,
                                 ACTIVE_MODEL_NOTE)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Representation"


@command(r"(?:colou?r )?(?:carbons?|carbon colou?rs?) (?:of model (?P<model>\S+) )?grey|"
         r"grey carbons?(?: (?:for|of) model (?P<model2>\S+))?",
         examples=["colour carbons grey", "grey carbons", "grey carbons for model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "model2": ArgType.MODEL},
         notes="Uses grey for carbon atoms. " + ACTIVE_MODEL_NOTE)
def grey_carbons(model: Optional[str] = None, model2: Optional[str] = None) -> str:
    """Use grey carbon colours for a model."""
    imol = resolve_model(first_present(model, model2))
    if coot is not None:
        coot.set_use_grey_carbons_for_molecule(imol, 1)
    return f"Using grey carbons for model {imol}"


@command(r"(?:colou?r )?(?:carbons?|carbon colou?rs?) (?:of model (?P<model>\S+) )?coloured|"
         r"colou?red carbons?(?: (?:for|of) model (?P<model2>\S+))?",
         examples=["colour carbons coloured", "coloured carbons",
                   "coloured carbons for model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "model2": ArgType.MODEL},
         notes="Uses per-element carbon colouring. " + ACTIVE_MODEL_NOTE)
def coloured_carbons(model: Optional[str] = None, model2: Optional[str] = None) -> str:
    """Use element-coloured carbons for a model."""
    imol = resolve_model(first_present(model, model2))
    if coot is not None:
        coot.set_use_grey_carbons_for_molecule(imol, 0)
    return f"Using coloured carbons for model {imol}"


# A ribbon/surface representation is added to a model as one or more extra
# "meshes" (the colour-by-chain ribbon adds one mesh per chain).  Coot's
# add_molecular_representation() returns 0 regardless of what it added, and its
# remove_molecular_representation() does not touch these meshes, so we cannot
# hide by a "representation index".  Instead we record the *mesh indices* that
# each show appended and toggle their draw state with set_draw_mesh(), which is
# the call that actually shows/hides a mesh.  Keyed by (imol, style) -> mesh
# indices.
_REPS: dict[tuple[int, str], list[int]] = {}

# secondary_structure_usage: 2 == CALC_SECONDARY_STRUCTURE, so ribbons get
# sensible helices/strands even when the model has no header SS records.
_CALC_SECONDARY_STRUCTURE = 2


def _mesh_count(imol: int) -> int:
    """Number of extra meshes attached to a model.

    draw_mesh_state() returns -1 for an out-of-range index, so we probe
    upward until it does - there is no direct "how many meshes" call.
    """
    n = 0
    while coot.draw_mesh_state(imol, n) != -1:
        n += 1
    return n


def _show_rep(imol: int, style: str, colour_scheme: str, label: str) -> str:
    """Show a ribbon/surface representation, tracking the meshes it adds.

    If we have shown this representation before, re-enable those meshes
    rather than building a duplicate set.
    """
    if coot is None:
        return f"Showing {label} for model {imol}"
    tracked = _REPS.get((imol, style))
    if tracked:
        for idx in tracked:
            coot.set_draw_mesh(imol, idx, 1)
        return f"Showing {label} for model {imol}"
    before = _mesh_count(imol)
    coot.add_molecular_representation_py(
        imol, "//", colour_scheme, style, _CALC_SECONDARY_STRUCTURE)
    after = _mesh_count(imol)
    if after <= before:
        return f"Could not add {label} to model {imol}"
    _REPS[(imol, style)] = list(range(before, after))
    return f"Showing {label} for model {imol}"


def _hide_rep(imol: int, style: str, label: str) -> str:
    """Hide a previously shown ribbon/surface representation, if any.

    The mesh indices stay tracked so a later "show" re-enables the same
    meshes instead of rebuilding them.
    """
    idxs = _REPS.get((imol, style))
    if not idxs:
        return f"No {label} to hide for model {imol}"
    if coot is not None:
        for idx in idxs:
            coot.set_draw_mesh(imol, idx, 0)
    return f"Hiding {label} for model {imol}"


@command(r"show ribbons?(?: (?:for|of) model (?P<model>\S+))?",
         examples=["show ribbons", "show ribbons for model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Draws a ribbon (cartoon) representation of the whole model, "
               "coloured by chain. " + ACTIVE_MODEL_NOTE)
def show_ribbons(model: Optional[str] = None) -> str:
    """Show a ribbon representation of a model."""
    imol = resolve_model(model)
    return _show_rep(imol, "Ribbon", "colorRampChainsScheme", "ribbons")


@command(r"hide ribbons?(?: (?:for|of) model (?P<model>\S+))?",
         examples=["hide ribbons", "hide ribbons for model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Removes the ribbon representation added by 'show ribbons'. "
               + ACTIVE_MODEL_NOTE)
def hide_ribbons(model: Optional[str] = None) -> str:
    """Hide the ribbon representation of a model."""
    imol = resolve_model(model)
    return _hide_rep(imol, "Ribbon", "ribbons")


@command(r"show (?:molecular )?surface(?: (?:for|of) model (?P<model>\S+))?",
         examples=["show surface", "show surface for model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Draws a molecular surface for the whole model, coloured by "
               "chain. " + ACTIVE_MODEL_NOTE)
def show_surface(model: Optional[str] = None) -> str:
    """Show a molecular surface of a model."""
    imol = resolve_model(model)
    return _show_rep(imol, "MolecularSurface", "Chain", "molecular surface")


@command(r"hide (?:molecular )?surface(?: (?:for|of) model (?P<model>\S+))?",
         examples=["hide surface", "hide surface for model 0"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Removes the molecular surface added by 'show surface'. "
               + ACTIVE_MODEL_NOTE)
def hide_surface(model: Optional[str] = None) -> str:
    """Hide the molecular surface of a model."""
    imol = resolve_model(model)
    return _hide_rep(imol, "MolecularSurface", "molecular surface")


@command(r"show symmetry",
         examples=["show symmetry"],
         category=CATEGORY)
def show_symmetry(**_: Optional[str]) -> str:
    """Show symmetry-related molecules."""
    if coot is not None:
        coot.set_show_symmetry_master(1)
    return "Showing symmetry"


@command(r"hide symmetry",
         examples=["hide symmetry"],
         category=CATEGORY)
def hide_symmetry(**_: Optional[str]) -> str:
    """Hide symmetry-related molecules."""
    if coot is not None:
        coot.set_show_symmetry_master(0)
    return "Hiding symmetry"
