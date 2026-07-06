# coot_commands/commands/delete.py
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

"""Commands that close (delete) molecules - a single map/model, or all."""

from __future__ import annotations

from typing import List, Optional

from coot_commands.registry import command
from coot_commands.types import as_int, CommandError

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Delete"


def _all_molecules() -> List[int]:
    """Molecule numbers of every loaded model and map."""
    if coot is None:
        return []
    return [i for i in range(coot.graphics_n_molecules())
            if coot.is_valid_model_molecule(i) or coot.is_valid_map_molecule(i)]


def _describe(imol: int) -> str:
    """A short 'map 1'/'model 0' label for a molecule number."""
    if coot is None:
        return f"molecule {imol}"
    if coot.is_valid_map_molecule(imol):
        return f"map {imol}"
    if coot.is_valid_model_molecule(imol):
        return f"model {imol}"
    return f"molecule {imol}"


def _close_one(imol: int, kind: Optional[str]) -> str:
    """Close one molecule, checking it is of the expected *kind* if given.

    *kind* is "map", "model" or None (either).  We check before closing so
    "delete map 0" on a model reports the mismatch rather than silently
    closing the wrong thing.
    """
    if coot is None:
        return f"Deleted {kind or 'molecule'} {imol}"
    is_map = bool(coot.is_valid_map_molecule(imol))
    is_model = bool(coot.is_valid_model_molecule(imol))
    if not (is_map or is_model):
        raise CommandError(f"no molecule {imol} to delete")
    if kind == "map" and not is_map:
        raise CommandError(f"molecule {imol} is not a map")
    if kind == "model" and not is_model:
        raise CommandError(f"molecule {imol} is not a model")
    label = _describe(imol)
    coot.close_molecule(imol)
    return f"Deleted {label}"


@command(r"(?:delete|close|remove) map (?P<imol>\d+)",
         examples=["delete map 1", "close map 1"],
         category=CATEGORY,
         notes="Closes the given map, freeing its molecule number.")
def delete_map(imol: str) -> str:
    """Delete (close) a map."""
    return _close_one(as_int(imol, "map number"), "map")


@command(r"(?:delete|close|remove) model (?P<imol>\d+)",
         examples=["delete model 0", "close model 0"],
         category=CATEGORY,
         notes="Closes the given model, freeing its molecule number.")
def delete_model(imol: str) -> str:
    """Delete (close) a model."""
    return _close_one(as_int(imol, "model number"), "model")


@command(r"(?:delete|close|remove) (?:molecule |mol )(?P<imol>\d+)",
         examples=["delete molecule 2", "close mol 2"],
         category=CATEGORY,
         notes="Closes the given molecule (map or model) by number.")
def delete_molecule(imol: str) -> str:
    """Delete (close) a molecule of either kind."""
    return _close_one(as_int(imol, "molecule number"), None)


# Registered before the plain "delete all" so the confirmed form matches first.
@command(r"(?:delete|close|remove) (?:all|everything) confirm",
         examples=["delete all confirm"],
         category=CATEGORY,
         notes="Confirms and carries out 'delete all', closing every loaded "
               "map and model.")
def delete_all_confirm(**_: Optional[str]) -> str:
    """Close every loaded map and model (confirmed)."""
    molecules = _all_molecules()
    if not molecules:
        return "Nothing to delete"
    if coot is not None:
        # Close from the highest number down: closing renumbers nothing here,
        # but iterating a stable snapshot keeps it robust either way.
        for imol in sorted(molecules, reverse=True):
            coot.close_molecule(imol)
    return f"Deleted all {len(molecules)} molecules"


@command(r"(?:delete|close|remove) (?:all|everything)",
         examples=["delete all", "close everything"],
         category=CATEGORY,
         notes="Asks for confirmation before closing every loaded map and "
               "model. Type 'delete all confirm' to go ahead.")
def delete_all(**_: Optional[str]) -> str:
    """Ask to confirm closing every loaded map and model."""
    molecules = _all_molecules()
    if not molecules:
        return "Nothing to delete"
    return (f"This will close all {len(molecules)} molecules. "
            "Type 'delete all confirm' to proceed.")
