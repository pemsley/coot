# coot_command_interface.py
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

"""Entry point for Coot's natural-language Command tab.

The C++ Command terminal calls :func:`run_command` with each line the
user types.  This module is deliberately thin: it delegates to the
:mod:`coot_commands` package, where the registry and the individual
command modules live.  Keeping this shim stable means the C++ side never
has to change as commands are added.
"""

from __future__ import annotations

from typing import Optional


def run_command(text: Optional[str]) -> str:
    """Parse *text* and execute the matching Coot command.

    Returns a string describing the outcome, suitable for display in the
    Command terminal.  Never raises: import, parsing, or execution errors
    are reported in the returned string instead.
    """
    if text is None or not text.strip():
        return ""

    try:
        from coot_commands.registry import dispatch
    except Exception as e:
        return f"Error: could not load command modules: {e}"

    try:
        result = dispatch(text)
    except Exception as e:
        import traceback
        detail = str(e) or repr(e)
        tb = traceback.format_exc().rstrip()
        return (f"Error running {text.strip()!r}: {type(e).__name__}: {detail}\n"
                f"{tb}")

    if result is None:
        return f"Unrecognised command: {text.strip()!r}  (try \"help\")"

    return result


def complete_command(text: Optional[str]) -> str:
    """Tab-complete the command line *text* for the Command tab.

    Returns a string the C++ side interprets: the first line is the new
    entry text (the input with the current word completed as far as is
    unambiguous); any following line lists the candidate words to display
    when the completion is ambiguous.  An empty string means "nothing to
    complete".  Never raises.
    """
    try:
        from coot_commands.completion import complete
        replacement, options = complete(text or "")
    except Exception:
        return ""

    if not replacement and not options:
        return ""
    if options:
        return replacement + "\n" + "  ".join(options)
    return replacement
