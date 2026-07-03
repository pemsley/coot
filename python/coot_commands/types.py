# coot_commands/types.py
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

try:
    import coot
except ImportError:
    coot = None  # allows this module to be imported/tested without Coot


class CommandError(Exception):
    """Raised when an argument cannot be resolved.

    ``dispatch`` lets this propagate; the entry point turns it into a
    readable message for the user, so handlers can validate freely.
    """


def _require_coot():
    if coot is None:
        raise CommandError("the Coot API is not available")
    return coot


def as_int(value, what="value"):
    """Coerce a captured string to an int, or raise CommandError."""
    try:
        return int(value)
    except (TypeError, ValueError):
        raise CommandError(f"expected a whole number for {what}, got {value!r}")


def as_float(value, what="value"):
    """Coerce a captured string to a float, or raise CommandError."""
    try:
        return float(value)
    except (TypeError, ValueError):
        raise CommandError(f"expected a number for {what}, got {value!r}")


def active_model():
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


def active_map():
    """Return the molecule number of the active (refinement) map, or raise."""
    c = _require_coot()
    imol = c.imol_refinement_map()
    if imol is not None and imol >= 0:
        return imol
    raise CommandError("no active map - open a map first")


def resolve_model(value=None):
    """Resolve a model reference to a molecule number (imol).

    Accepts a plain integer, or ``None``/"active"/"this"/"current" to mean
    the active model.  This is also the hook for richer references later
    ("the ligand", a model by name); keeping it here means every command
    that takes a model gains them at once.
    """
    if value is None or str(value).lower() in ("active", "this", "current", "it"):
        return active_model()
    return as_int(value, "model number")


def resolve_map(value=None):
    """Resolve a map reference to a molecule number, defaulting to the active map."""
    if value is None or str(value).lower() in ("active", "this", "current", "it"):
        return active_map()
    return as_int(value, "map number")


# Named colours -> (r, g, b) floats in 0..1, for map/background colour commands.
COLOURS = {
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


def resolve_colour(name):
    """Resolve a colour name to an (r, g, b) tuple, or raise CommandError."""
    key = str(name).strip().lower()
    if key in COLOURS:
        return COLOURS[key]
    raise CommandError(
        f"unknown colour {name!r} (known: {', '.join(sorted(set(COLOURS)))})")
