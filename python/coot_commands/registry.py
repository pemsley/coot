# coot_commands/registry.py
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

"""Command registry and matching engine for the Coot command interface.

Command modules register handlers with the :func:`command` decorator,
pairing a regular expression (with named groups for arguments) to a
function.  :func:`dispatch` normalises the user's text, finds the first
matching pattern, and calls its handler with the captured groups as
keyword arguments.

All matching logic (whitespace normalisation, case-insensitivity, and in
future fuzzy "did you mean" suggestions) lives here, so every command
benefits without having to think about it.
"""

from __future__ import annotations

import re
from typing import Callable, Iterable, Optional

# A command handler takes keyword arguments (the regex named groups, each a
# str or None) and returns a message string.
Handler = Callable[..., str]


class Command:
    """A single registered command: a pattern, a handler, and metadata.

    The metadata (``help_text``, ``examples``, ``category``, ``notes`` and
    the handler docstring) is what both the in-app ``help`` command and the
    Markdown documentation generator read, so documenting a command is just
    a matter of filling these in on the decorator.
    """

    def __init__(self, pattern: str, handler: Handler, help_text: str,
                 examples: Iterable[str], category: str,
                 notes: Optional[str]) -> None:
        self.regex = re.compile(pattern, re.IGNORECASE)
        self.pattern = pattern
        self.handler = handler
        self.help_text = help_text
        self.examples = tuple(examples)
        self.category = category
        self.notes = notes

    @property
    def name(self) -> str:
        return self.handler.__name__

    @property
    def description(self) -> str:
        """The long-form description: the handler docstring, if any."""
        return (self.handler.__doc__ or "").strip()


# The ordered list of all registered commands.  The first pattern that
# matches wins, so more specific patterns should be registered first
# (within a module, that means declaring them earlier).
_COMMANDS = []


def command(pattern: str, help: Optional[str] = None,
            examples: Iterable[str] = (), category: str = "General",
            notes: Optional[str] = None) -> Callable[[Handler], Handler]:
    """Decorator registering *handler* for inputs matching *pattern*.

    *pattern* is a regular expression matched (case-insensitively) against
    the whitespace-normalised input.  Named groups become keyword
    arguments passed to the handler, e.g. ``(?P<model>\\S+)`` calls
    ``handler(model="0")``.  Make an argument optional (e.g. a trailing
    ``(?: (?P<model>\\S+))?``) and it arrives as ``None`` - the shared
    resolvers in :mod:`coot_commands.types` then fall back to the active
    molecule.

    Documentation metadata:

    * *help* - one-line summary (defaults to the handler docstring's first
      line) shown by the ``help`` command and as the row summary in docs.
    * *examples* - example input strings; the first is treated as the
      canonical form. Use British spelling as the primary example.
    * *category* - groups the command in ``help`` output and the generated
      Markdown (e.g. "View", "Maps", "Navigation").
    * *notes* - optional extra prose for the generated docs (caveats,
      argument details) that would be too long for *help*.
    """
    def decorator(handler: Handler) -> Handler:
        doc_first_line = (handler.__doc__ or "").strip().splitlines()
        default_help = doc_first_line[0] if doc_first_line else ""
        help_text = help if help is not None else default_help
        _COMMANDS.append(Command(pattern, handler, help_text, examples, category, notes))
        return handler
    return decorator


def all_commands() -> list[Command]:
    """Return the list of registered :class:`Command` objects."""
    return list(_COMMANDS)


def normalise(text: str) -> str:
    """Collapse runs of whitespace and strip the ends."""
    return re.sub(r"\s+", " ", text.strip())


def dispatch(text: str) -> Optional[str]:
    """Find and run the command matching *text*.

    Returns the handler's result string, or ``None`` if nothing matched
    (the caller decides how to report an unrecognised command).  Handler
    exceptions propagate to the caller.
    """
    norm = normalise(text)
    for cmd in _COMMANDS:
        match = cmd.regex.match(norm)
        if match:
            return cmd.handler(**match.groupdict())
    return None
