# coot_commands/types.py
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

"""Shared argument coercion and entity resolution for commands.

Command handlers receive their regex named groups as strings (or ``None``
for an omitted optional group).  The helpers here turn those into the
values the Coot API wants - integers, molecule numbers resolved to the
active model/map, RGB colours - in one place, so handlers stay short and
behave consistently.

When a molecule is not named, we fall back to the *active* molecule (the
one the GUI menus would act on): ``resolve_model`` uses the active
residue's molecule, and ``resolve_map`` uses the refinement map.
"""

from __future__ import annotations

import enum
from typing import Any, List, Optional

try:
    import coot
except ImportError:
    coot = None  # allows this module to be imported/tested without Coot


class CommandError(Exception):
    """Raised when an argument cannot be resolved.

    ``dispatch`` lets this propagate; the entry point turns it into a
    readable message for the user, so handlers can validate freely.
    """


def _require_coot() -> Any:
    if coot is None:
        raise CommandError("the Coot API is not available")
    return coot


def as_int(value: Optional[str], what: str = "value") -> int:
    """Coerce a captured string to an int, or raise CommandError."""
    try:
        return int(value)
    except (TypeError, ValueError):
        raise CommandError(f"expected a whole number for {what}, got {value!r}")


def as_float(value: Optional[str], what: str = "value") -> float:
    """Coerce a captured string to a float, or raise CommandError."""
    try:
        return float(value)
    except (TypeError, ValueError):
        raise CommandError(f"expected a number for {what}, got {value!r}")


def active_model() -> int:
    """Return the molecule number of the active model, or raise.

    The active model is the one under the pointer / last picked, matching
    the "active molecule" the GUI menus use.  The Python binding is
    ``active_residue_py``; we fall back to ``active_residue`` and then to
    the first coordinate molecule for robustness across builds.
    """
    c = _require_coot()
    getter = getattr(c, "active_residue_py", None) or getattr(c, "active_residue", None)
    if getter is not None:
        active = getter()
        if active:
            return active[0]
    first = getattr(c, "first_coords_imol", None)
    if first is not None:
        imol = first()
        if imol is not None and imol >= 0:
            return imol
    raise CommandError("no active model - load or click a model first")


def active_map() -> int:
    """Return the molecule number of the active (refinement) map, or raise."""
    c = _require_coot()
    imol = c.imol_refinement_map()
    if imol is not None and imol >= 0:
        return imol
    raise CommandError("no active map - open a map first")


def resolve_model(value: Optional[str] = None) -> int:
    """Resolve a model reference to a molecule number (imol).

    Accepts a plain integer, or ``None``/"active"/"this"/"current" to mean
    the active model.  This is also the hook for richer references later
    ("the ligand", a model by name); keeping it here means every command
    that takes a model gains them at once.
    """
    if value is None or str(value).lower() in ("active", "this", "current", "it"):
        return active_model()
    return as_int(value, "model number")


def resolve_map(value: Optional[str] = None) -> int:
    """Resolve a map reference to a molecule number, defaulting to the active map."""
    if value is None or str(value).lower() in ("active", "this", "current", "it"):
        return active_map()
    return as_int(value, "map number")


# Named colours -> (r, g, b) floats in 0..1, for map/background colour commands.
COLOURS: dict[str, tuple[float, float, float]] = {
    "black":   (0.0, 0.0, 0.0),
    "white":   (1.0, 1.0, 1.0),
    "red":     (1.0, 0.0, 0.0),
    "green":   (0.0, 1.0, 0.0),
    "blue":    (0.0, 0.0, 1.0),
    "yellow":  (1.0, 1.0, 0.0),
    "cyan":    (0.0, 1.0, 1.0),
    "magenta": (1.0, 0.0, 1.0),
    "orange":  (1.0, 0.5, 0.0),
    "purple":  (0.5, 0.0, 0.5),
    "grey":    (0.5, 0.5, 0.5),
    "gray":    (0.5, 0.5, 0.5),
    "pink":    (1.0, 0.6, 0.8),
    "sky":     (0.5, 0.7, 1.0),
    "salmon":  (1.0, 0.6, 0.5),
}


def resolve_colour(name: str) -> tuple[float, float, float]:
    """Resolve a colour name to an (r, g, b) tuple, or raise CommandError."""
    key = str(name).strip().lower()
    if key in COLOURS:
        return COLOURS[key]
    raise CommandError(
        f"unknown colour {name!r} (known: {', '.join(sorted(set(COLOURS)))})")


# ---------------------------------------------------------------------------
#   Argument types (for tab completion)
# ---------------------------------------------------------------------------
#
# A command declares the *type* of each of its regex named groups via the
# ``arg_types`` argument to the ``@command`` decorator, e.g.
# ``arg_types={"model": ArgType.MODEL}``.  The completion engine
# (:mod:`coot_commands.completion`) asks the type for its current candidate
# values, so pressing Tab after "show model " lists the loaded models.  The
# type is metadata only - it does not affect how a command runs.


def loaded_models() -> List[str]:
    """Molecule numbers of the currently loaded models, as strings."""
    if coot is None:
        return []
    try:
        return [str(i) for i in range(coot.graphics_n_molecules())
                if coot.is_valid_model_molecule(i) == 1]
    except Exception:
        return []


def loaded_maps() -> List[str]:
    """Molecule numbers of the currently loaded maps, as strings."""
    if coot is None:
        return []
    try:
        return [str(i) for i in range(coot.graphics_n_molecules())
                if coot.is_valid_map_molecule(i) == 1]
    except Exception:
        return []


def colour_names() -> List[str]:
    """The named colours accepted by the colour/background commands."""
    return sorted(set(COLOURS))


class ArgType(enum.Enum):
    """The kind of value a command argument accepts.

    Each member knows how to enumerate its current candidate values via
    :meth:`candidates`, which the completion engine calls when the user
    presses Tab at that argument's position.  Add a member here (and a
    branch in :meth:`candidates`) to teach completion about a new kind of
    argument; commands then opt in with ``arg_types={"group": ArgType.X}``.
    """

    MODEL = "model"
    MAP = "map"
    COLOUR = "colour"

    def candidates(self) -> List[str]:
        """Return the current completion candidates for this type."""
        if self is ArgType.MODEL:
            return loaded_models()
        if self is ArgType.MAP:
            return loaded_maps()
        if self is ArgType.COLOUR:
            return colour_names()
        return []
