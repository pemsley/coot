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
from typing import Callable, Dict, Iterable, Optional

# A command handler takes keyword arguments (the regex named groups, each a
# str or None) and returns a message string.
Handler = Callable[..., str]

# Maps a regex named group to the kind of value it accepts.  Used only by the
# completion engine; the value is typically an ``ArgType`` member (see
# :mod:`coot_commands.types`) but any object with a ``candidates()`` method, a
# zero-argument callable, or a plain iterable of strings works too.
ArgTypes = Dict[str, object]


class Command:
    """A single registered command: a pattern, a handler, and metadata.

    The metadata (``help_text``, ``examples``, ``category``, ``notes`` and
    the handler docstring) is what both the in-app ``help`` command and the
    Markdown documentation generator read, so documenting a command is just
    a matter of filling these in on the decorator.
    """

    def __init__(self, pattern: str, handler: Handler, help_text: str,
                 examples: Iterable[str], category: str,
                 notes: Optional[str],
                 arg_types: Optional[ArgTypes] = None) -> None:
        self.regex = re.compile(pattern, re.IGNORECASE)
        self.pattern = pattern
        self.handler = handler
        self.help_text = help_text
        self.examples = tuple(examples)
        self.category = category
        self.notes = notes
        self.arg_types = dict(arg_types) if arg_types else {}

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
            notes: Optional[str] = None,
            arg_types: Optional[ArgTypes] = None) -> Callable[[Handler], Handler]:
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
    * *arg_types* - maps a named group to the kind of value it accepts (an
      :class:`~coot_commands.types.ArgType`), so tab completion can offer
      live candidates, e.g. ``arg_types={"model": ArgType.MODEL}``.
    """
    def decorator(handler: Handler) -> Handler:
        doc_first_line = (handler.__doc__ or "").strip().splitlines()
        default_help = doc_first_line[0] if doc_first_line else ""
        help_text = help if help is not None else default_help
        _COMMANDS.append(Command(pattern, handler, help_text, examples,
                                 category, notes, arg_types))
        return handler
    return decorator


def all_commands() -> list[Command]:
    """Return the list of registered :class:`Command` objects."""
    return list(_COMMANDS)


def normalise(text: str) -> str:
    """Collapse runs of whitespace and strip the ends."""
    return re.sub(r"\s+", " ", text.strip())


def unmatched_examples() -> list[tuple[str, str]]:
    """Return ``(command_name, example)`` pairs that don't match their command.

    Every example should be dispatchable by its own command - otherwise the
    phrasing advertised in ``help`` and the reference silently fails, and
    tab completion (which learns a command's shape by matching examples
    against the regex) mis-classifies the example's words.  Tests and the
    docs generator call this to catch such drift at authoring time.
    """
    problems = []
    for cmd in _COMMANDS:
        for example in cmd.examples:
            if not cmd.regex.match(normalise(example)):
                problems.append((cmd.name, example))
    return problems


def dispatch(text: str) -> Optional[str]:
    """Find and run the command matching *text*.

    Returns the handler's result string, or ``None`` if nothing matched
    (the caller decides how to report an unrecognised command).  Handler
    exceptions propagate to the caller.

    The input is first passed through :func:`coot_commands.speech.from_speech`,
    which rewrites dictated forms ("model zero" -> "model 0") to canonical
    text, so spoken commands (e.g. macOS Dictation typed into the Command tab)
    match the same patterns as typed ones.  It is a no-op on already-canonical
    input.
    """
    from coot_commands.speech import from_speech
    norm = normalise(from_speech(text))
    for cmd in _COMMANDS:
        match = cmd.regex.match(norm)
        if match:
            return cmd.handler(**match.groupdict())
    return None
