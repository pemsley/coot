# coot_commands/docs.py
#
# Copyright 2026 Jordan Dialpuri, Medical Research Council Laboratory of Molecular Biology

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

"""Generate documentation from the registered commands.

Because every command carries its own help text, examples, category and
notes (see :func:`coot_commands.registry.command`), the reference
documentation can be generated straight from the registry - so it never
drifts from the code.  Run this module as a script to (re)write the
Markdown reference::

    python3 -m coot_commands.docs > doc/command-reference.md

or from inside Coot::

    import coot_commands.docs as d
    print(d.to_markdown())
"""

from __future__ import annotations

from collections import OrderedDict

import coot_commands  # noqa: F401  - ensures all command modules are imported
from coot_commands.registry import Command, all_commands


def commands_by_category() -> "OrderedDict[str, list[Command]]":
    """Return an OrderedDict of category -> list of Command, sorted."""
    grouped: OrderedDict[str, list[Command]] = OrderedDict()
    for cmd in sorted(all_commands(), key=lambda c: (c.category, c.name)):
        grouped.setdefault(cmd.category, []).append(cmd)
    return grouped


def to_markdown() -> str:
    """Render the full command reference as a Markdown string."""
    lines = [
        "# Coot command reference",
        "",
        "The commands below are typed into the **Command** tab of the "
        "Python/AI terminal. Text is matched case-insensitively and extra "
        "whitespace is ignored. Where a command takes a model or map number, "
        "omitting it acts on the *active* molecule.",
        "",
        "> This file is generated from the command definitions "
        "(`python/coot_commands/`). Do not edit by hand - run "
        "`python3 -m coot_commands.docs` to regenerate.",
        "",
    ]
    for category, commands in commands_by_category().items():
        lines.append(f"## {category}")
        lines.append("")
        for cmd in commands:
            example = cmd.examples[0] if cmd.examples else cmd.name
            lines.append(f"### `{example}`")
            lines.append("")
            if cmd.help_text:
                lines.append(cmd.help_text)
                lines.append("")
            if len(cmd.examples) > 1:
                lines.append("Examples:")
                lines.append("")
                for ex in cmd.examples:
                    lines.append(f"- `{ex}`")
                lines.append("")
            if cmd.notes:
                lines.append(cmd.notes)
                lines.append("")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


if __name__ == "__main__":
    import sys
    sys.stdout.write(to_markdown())
