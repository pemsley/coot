# coot_commands/commands/help.py
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

"""The built-in ``help`` command.

This lists every registered command, grouped by category, using the
metadata each one carries - so the help text stays correct automatically
as commands are added.
"""

from __future__ import annotations

from coot_commands.registry import command, all_commands


@command(r"(?:help|commands|what can i do)\??",
         examples=["help"],
         category="Help")
def show_help() -> str:
    """List the available commands, grouped by category."""
    by_category = {}
    for cmd in all_commands():
        by_category.setdefault(cmd.category, []).append(cmd)

    lines = ["Available commands (type any of the examples):", ""]
    for category in sorted(by_category):
        lines.append(f"{category}:")
        for cmd in by_category[category]:
            example = cmd.examples[0] if cmd.examples else cmd.name
            summary = cmd.help_text or ""
            lines.append(f"  {example:<26} {summary}")
        lines.append("")
    return "\n".join(lines).rstrip()
