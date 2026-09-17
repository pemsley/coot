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


def slug(category: str) -> str:
    """The GitHub-style anchor for a ``## category`` heading.

    Headings become anchors by lower-casing, dropping anything that is not
    alphanumeric, a space or a hyphen, and turning spaces into hyphens - so
    "Model editing" links as "#model-editing".
    """
    kept = [c for c in category.lower() if c.isalnum() or c in " -"]
    return "".join(kept).strip().replace(" ", "-")


# The prose that opens the reference.  It explains what the command
# interface is for and the handful of conventions that hold across every
# command, so each entry below needs only to describe itself.
_INTRO = """\
The commands below are typed into the **Command** tab of the Python/AI
terminal. They exist so that the things you do most often while building a
model are a short phrase rather than a hunt through the menus or a remembered
API call: `refine chain A`, `fetch 4hhb`, `colour map 1 blue`, `go to A/89`.

The guiding idea is that you should be able to type roughly what you mean.
Each command is a pattern rather than a fixed string, so the obvious synonyms
and phrasings reach the same place - `show`/`display`, `hide`/`undisplay`,
`centre on`/`go to` - and both British and American spellings of "colour" are
accepted. Matching ignores case and collapses runs of whitespace, and dictated
input is normalised first, so "model zero" spoken into the field works exactly
like "model 0" typed.

Arguments follow the same principle. Where a command takes a model or map
number, omitting it acts on the *active* molecule, which is usually the one
you mean; where it takes a residue you can name it (`A/89`) or leave it out
and act on the residue at the centre of the screen. Tab completion fills in
what it can and offers the live candidates - the molecules actually loaded,
the colours actually available, the files actually on disk.

Every command is a small Python function that carries its own help text,
examples, category and notes alongside it (see
`python/coot_commands/registry.py`). That metadata is the single source of
truth: the in-app `help` command and this reference are both generated from
it, so neither can drift from what the code really accepts. Adding a command
means writing the function and filling in the decorator; the documentation
follows for free."""


def to_markdown() -> str:
    """Render the full command reference as a Markdown string."""
    grouped = commands_by_category()
    lines = [
        "# Coot command reference",
        "",
        _INTRO,
        "",
        "> This file is generated from the command definitions "
        "(`python/coot_commands/`). Do not edit by hand - run "
        "`python3 -m coot_commands.docs` to regenerate.",
        "",
        "## Contents",
        "",
    ]
    for category, commands in grouped.items():
        count = len(commands)
        plural = "" if count == 1 else "s"
        lines.append(f"- [{category}](#{slug(category)}) - "
                     f"{count} command{plural}")
    lines.append("")

    for category, commands in grouped.items():
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
    from coot_commands.registry import unmatched_examples

    # Warn (but still generate) if any advertised example can't be dispatched
    # by its own command - such examples fail when typed and mislead tab
    # completion. See registry.unmatched_examples.
    for name, example in unmatched_examples():
        sys.stderr.write(
            f"warning: example {example!r} does not match command {name}\n")

    sys.stdout.write(to_markdown())
