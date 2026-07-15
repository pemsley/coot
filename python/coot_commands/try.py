# coot_commands/try.py
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

"""Dry-run a command line against the registry without running Coot.

Answers the two questions you actually have while writing a pattern:

* **Which command wins?**  ``dispatch`` runs the *first* registered
  pattern that matches, so this shows that command, the arguments it would
  be called with, and where it lives.
* **Is anything shadowing it?**  It also lists every *other* command whose
  pattern matches the same input.  Those are unreachable for this input
  (the first one wins) - a classic cause of "my new command never fires".

It never calls the handler, so it is safe to run outside Coot::

    python3 -m coot_commands.try "refine chain A"

With no argument it reads lines from stdin, so you can pipe or type several
inputs to probe.
"""

from __future__ import annotations

import inspect
import sys
from typing import List, Optional, Tuple

import coot_commands  # noqa: F401  - imports all command modules so they register
from coot_commands.registry import Command, all_commands, normalise


def _location(cmd: Command) -> str:
    """`file.py:line` for a command's handler, for a clickable pointer."""
    try:
        path = inspect.getsourcefile(cmd.handler) or "?"
        line = cmd.handler.__code__.co_firstlineno
        return f"{path.rsplit('/', 1)[-1]}:{line}"
    except (TypeError, OSError):
        return "?"


def _format_call(cmd: Command, groups: dict) -> str:
    """Render the handler call that would be made, e.g. ``f(chain='A')``."""
    args = ", ".join(f"{k}={v!r}" for k, v in groups.items())
    return f"{cmd.name}({args})"


def matches(text: str) -> List[Tuple[Command, dict]]:
    """Every command whose pattern matches *text*, in registration order.

    The first entry is the one ``dispatch`` would run; the rest are
    shadowed for this input.  Each is paired with the captured groups
    (``match.groupdict()``) it would pass to its handler.
    """
    norm = normalise(text)
    found = []
    for cmd in all_commands():
        m = cmd.regex.match(norm)
        if m:
            found.append((cmd, m.groupdict()))
    return found


def explain(text: str) -> str:
    """Human-readable dry-run report for a single input line."""
    norm = normalise(text)
    found = matches(text)
    lines = [f"input (normalised): {norm!r}", ""]
    if not found:
        lines.append("no command matched - nothing would run.")
        return "\n".join(lines)

    winner, groups = found[0]
    lines.append(f"MATCH  {winner.name}  [{winner.category}]  {_location(winner)}")
    if groups:
        for key, value in groups.items():
            lines.append(f"    {key} = {value!r}")
    else:
        lines.append("    (no arguments captured)")
    lines.append(f"    would call: {_format_call(winner, groups)}")

    if len(found) > 1:
        lines.append("")
        lines.append("also matched (shadowed - the first match above wins):")
        for cmd, _ in found[1:]:
            lines.append(f"    {cmd.name}  [{cmd.category}]  {_location(cmd)}")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv:
        print(explain(" ".join(argv)))
        return 0
    # No argument: treat each stdin line as an input to probe.
    for raw in sys.stdin:
        line = raw.strip()
        if not line:
            continue
        print(explain(line))
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
