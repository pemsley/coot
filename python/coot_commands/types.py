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


def first_present(*values: Optional[str]) -> Optional[str]:
    """Return the first non-``None`` value, or ``None`` if all are ``None``.

    Handy where regex alternation captures the same argument under two group
    names (e.g. ``model`` and ``model2``): ``resolve_model(first_present(
    model, model2))`` replaces the ``model if model is not None else model2``
    ternary that would otherwise appear in every such handler.
    """
    for value in values:
        if value is not None:
            return value
    return None


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


# Canonical wording for the "argument is optional and defaults to the active
# molecule" note.  Commands whose model/map group is optional append one of
# these to their own ``notes`` (concatenating with a space) so the sentence
# reads identically everywhere, rather than each command paraphrasing it.
ACTIVE_MODEL_NOTE = "With no model number, acts on the active model."
ACTIVE_MAP_NOTE = "With no map number, acts on the active map."
ACTIVE_RESIDUE_NOTE = ("With no residue named, acts on the active residue "
                       "(the one at the centre of the screen).")


# A residue reference like "A/89" or "A 89" - a single-character chain id and a
# (possibly negative) residue number, separated by a space or a slash.  Many
# commands take one, so the regex fragment and its parser live here rather than
# being copy-pasted into every module.  Embed the fragment in a pattern and
# pass the captured groups straight to :func:`parse_res`.
RES_SPEC = r"(?P<chain>[A-Za-z0-9])[ /](?P<resno>-?\d+)"

# The same reference made optional, with its own leading space: when the user
# omits it the ``chain``/``resno`` groups arrive as ``None`` and the handler
# falls back to the active residue via :func:`resolve_residue`.  Embed it
# straight after the command word, e.g. ``rf"backrub{OPT_RES_SPEC}"``.  This is
# what lets every per-residue command act on the residue at the screen centre
# when no residue is named - and, just as importantly, keeps the command
# *matching* in that case rather than falling through to a shorter command.
OPT_RES_SPEC = r"(?: " + RES_SPEC + r")?"


def parse_res(chain: Optional[str], resno: Optional[str]) -> tuple[str, int]:
    """Turn captured ``chain``/``resno`` groups into (upper chain id, int)."""
    return str(chain).upper(), as_int(resno, "residue number")


def active_model() -> int:
    """Return the molecule number of the active model, or raise.

    The active model is the one holding the *active residue* - the residue
    at the centre of the screen (see :func:`active_residue`) - matching the
    "active molecule" the GUI menus act on.  The Python binding is
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
    raise CommandError("no active model - load a model first")


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


def require_model(imol: int) -> int:
    """Return *imol* if it is a loaded model molecule, else raise CommandError.

    Use this for commands that name an explicit model number (rather than
    defaulting to the active one) so a wrong or stale number - or a map number
    - is reported clearly instead of being passed to the Coot API, where it
    would silently misbehave or crash.
    """
    if coot is None:
        return imol
    if coot.is_valid_model_molecule(imol) == 1:
        return imol
    if coot.is_valid_map_molecule(imol) == 1:
        raise CommandError(f"molecule {imol} is a map, not a model")
    raise CommandError(f"no model molecule {imol}")


def require_map(imol: int) -> int:
    """Return *imol* if it is a loaded map molecule, else raise CommandError.

    The map analogue of :func:`require_model`.
    """
    if coot is None:
        return imol
    if coot.is_valid_map_molecule(imol) == 1:
        return imol
    if coot.is_valid_model_molecule(imol) == 1:
        raise CommandError(f"molecule {imol} is a model, not a map")
    raise CommandError(f"no map molecule {imol}")


def active_residue() -> tuple[int, str, int, str]:
    """Return (imol, chain_id, resno, ins_code) of the active residue, or raise.

    The active residue is the one at the **centre of the screen** - Coot
    picks the displayed atom closest to the rotation centre (with CA
    substitution), so it tracks wherever the view is centred, not the last
    click.  The Python binding is ``active_residue_py`` (falling back to
    ``active_residue``), returning ``[imol, chain, resno, ins_code,
    atom_name, ...]``.
    """
    c = _require_coot()
    getter = getattr(c, "active_residue_py", None) or getattr(c, "active_residue", None)
    active = getter() if getter is not None else None
    if not active:
        raise CommandError("no active residue - centre on a residue first")
    return active[0], active[1], active[2], active[3]


def resolve_residue(chain: Optional[str] = None, resno: Optional[str] = None,
                    model: Optional[str] = None) -> tuple[int, str, int, str]:
    """Resolve a residue reference to (imol, chain_id, resno, ins_code).

    With *chain* and *resno* given, parses them (via :func:`parse_res`)
    against the resolved model.  With neither given, falls back to the
    *active* residue - the one at the centre of the screen (see
    :func:`active_residue`) - which also fixes the model.  This is the
    residue analogue of :func:`resolve_model`, so any command that names a
    residue gains the "act on the residue at the screen centre" default just
    by making its spec optional.
    """
    if chain is not None and resno is not None:
        imol = resolve_model(model)
        chain_id, res = parse_res(chain, resno)
        return imol, chain_id, res, ""
    return active_residue()


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


def molecule_name(imol: int) -> str:
    """Short, path-less name of a molecule, or "" if unavailable.

    Prefers ``molecule_name_stub_py`` (name without directory or file
    extension) and falls back to the full ``molecule_name``.
    """
    if coot is None:
        return ""
    try:
        name = coot.molecule_name_stub_py(imol, 0)
        if name:
            return name
    except Exception:
        pass
    try:
        return coot.molecule_name(imol) or ""
    except Exception:
        return ""


def describe_molecule(imol: int) -> str:
    """A short ``'map N'`` / ``'model N'`` label for a molecule number.

    Names the *kind* of molecule (not its filename), for messages like
    "Deleted map 1". Falls back to ``"molecule N"`` when the number is neither
    a valid map nor a valid model, or when Coot is unavailable.
    """
    if coot is None:
        return f"molecule {imol}"
    try:
        if coot.is_valid_map_molecule(imol):
            return f"map {imol}"
        if coot.is_valid_model_molecule(imol):
            return f"model {imol}"
    except Exception:
        pass
    return f"molecule {imol}"


def molecule_label(value: str) -> str:
    """Decorate a molecule-number string with its name, for display lists.

    Returns e.g. ``"0 (tutorial-modern)"`` for completion menus, while the
    value inserted into the command line stays the bare number.  Falls back
    to the bare value when the number or name cannot be resolved.
    """
    try:
        name = molecule_name(int(value))
    except (TypeError, ValueError):
        return value
    return f"{value} ({name})" if name else value


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

    def label(self, value: str) -> str:
        """A display label for a candidate value, shown in Tab completion.

        The value inserted into the command line is unchanged; this only
        affects what the user sees in the option list.  Model and map
        numbers gain their molecule name (e.g. "0 (tutorial-modern)").
        """
        if self is ArgType.MODEL or self is ArgType.MAP:
            return molecule_label(value)
        return value
