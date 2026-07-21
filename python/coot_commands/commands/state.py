# coot_commands/commands/state.py
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

"""Session state commands: undo (and, in future, redo)."""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "State"


@command(r"undo",
         examples=["undo"],
         category=CATEGORY,
         notes="Undoes the last modification.")
def undo(**_: Optional[str]) -> str:
    """Undo the last action."""
    if coot is not None:
        coot.apply_undo()
    return "Undid the last action"

