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

"""Natural-language command interface for Coot.

The C++ "Command" terminal tab passes each line the user types to
``run_command()``.  This module parses that text into a Coot Python API
call and executes it.  Because all of the parsing lives here (pure
Python), new commands can be added by editing this file and
re-installing - no C++ recompile is needed.

To add a command, write a handler that takes the regex ``match`` object
and returns a human-readable string, then register it in ``COMMANDS``
with a regex pattern.  Patterns are matched against the lower-cased,
whitespace-normalised input.
"""

import re

try:
    import coot
except ImportError:
    coot = None


def _set_model_displayed(match, state):
    imol = int(match.group("imol"))
    if coot is not None:
        coot.set_mol_displayed(imol, state)
    return ("Showing" if state else "Hiding") + " model %d" % imol


def _cmd_show_model(match):
    return _set_model_displayed(match, 1)


def _cmd_hide_model(match):
    return _set_model_displayed(match, 0)


# Each entry is (compiled_pattern, handler).  Handlers receive the regex
# match object and return a message string.  The first pattern that
# matches wins, so order from most to least specific.
COMMANDS = [
    (re.compile(r"^(?:show|display)\s+model\s+(?P<imol>\d+)$"),     _cmd_show_model),
    (re.compile(r"^(?:hide|undisplay)\s+model\s+(?P<imol>\d+)$"),   _cmd_hide_model),
]


def run_command(text):
    """Parse *text* and execute the matching Coot command.

    Returns a string describing the outcome, suitable for display in the
    Command terminal.  Never raises: parsing or execution errors are
    reported in the returned string instead.
    """
    if text is None:
        return ""
    # Normalise: collapse internal whitespace and lower-case for matching.
    stripped = text.strip()
    if not stripped:
        return ""
    normalised = re.sub(r"\s+", " ", stripped).lower()

    for pattern, handler in COMMANDS:
        match = pattern.match(normalised)
        if match:
            try:
                return handler(match)
            except Exception as e:
                return "Error running %r: %s" % (stripped, e)

    return "Unrecognised command: %r" % stripped
