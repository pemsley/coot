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


def _read_clipboard() -> Optional[str]:
    """Return the system clipboard text, or ``None`` if it can't be read.

    Tries the platform clipboard readers directly (macOS ``pbpaste`` first,
    then the common Linux tools) rather than a GTK clipboard read, which is
    asynchronous in GTK4 and awkward to use from a synchronous command.
    """
    import shutil
    import subprocess

    for argv in (["pbpaste"], ["wl-paste", "--no-newline"],
                 ["xclip", "-selection", "clipboard", "-o"], ["xsel", "-b"]):
        if shutil.which(argv[0]) is None:
            continue
        try:
            out = subprocess.run(argv, capture_output=True, text=True, timeout=5)
        except Exception:
            continue
        if out.returncode == 0:
            return out.stdout
    return None


def run_from_clipboard() -> str:
    """Run the clipboard contents as a command, and return the result.

    A bridge for voice input: macOS Dictation types into native text fields
    but not into Coot's VTE terminal, so dictate a command into any native
    field (Spotlight, Notes, ...), copy it, then call this to run it.  The
    text goes through :func:`run_command`, so spoken wording ("model zero")
    is rewritten to canonical form just like typed input.

    Bind it to a key or a toolbar button for a keyboard-free flow; called
    with no argument it is short enough to type as the one keystroke voice
    can't supply.
    """
    text = _read_clipboard()
    if text is None:
        result = ("Could not read the clipboard (no pbpaste/xclip/xsel/"
                  "wl-paste found)")
    elif not text.strip():
        result = "Clipboard is empty - dictate and copy a command first"
    else:
        result = run_command(text)
    print(result)
    return result
