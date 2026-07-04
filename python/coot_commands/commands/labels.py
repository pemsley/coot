# coot_commands/commands/labels.py
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

"""Commands for atom labels."""

from coot_commands.registry import command

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Labels"


@command(r"(?:clear|remove) (?:all )?labels",
         examples=["clear labels", "remove all labels"],
         category=CATEGORY)
def clear_labels(**_):
    """Remove all atom labels."""
    if coot is not None:
        coot.remove_all_atom_labels()
    return "Cleared all atom labels"
