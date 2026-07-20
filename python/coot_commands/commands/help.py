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

from coot_commands.registry import command


@command(r"(?:help|commands|what can i do)\??",
         examples=["help"],
         category="Help")
def show_help() -> str:
    """List the available commands, grouped by category."""
    # Reuse the docs generator's grouping so help and the Markdown reference
    # order commands identically.  Imported lazily to avoid an import cycle
    # during command auto-discovery.
    from coot_commands.docs import commands_by_category

    lines = ["Available commands (type any of the examples):", ""]
    for category, cmds in commands_by_category().items():
        lines.append(f"{category}:")
        for cmd in cmds:
            example = cmd.examples[0] if cmd.examples else cmd.name
            summary = cmd.help_text or ""
            lines.append(f"  {example:<26} {summary}")
        lines.append("")
    return "\n".join(lines).rstrip()
