# coot_commands/completion.py
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

"""Tab completion for the Coot command interface.

The Command tab calls :func:`complete` with whatever the user has typed so
far and completes the current word, the way a shell does.  There are two
kinds of completion, and both fall out of the command definitions with no
per-command bookkeeping:

* *Keywords* - the fixed words of a command ("show", "model", "map").
  These come from each command's ``examples``.  We match an example
  against the command's own regex to learn which of its words are literal
  keywords and which are argument values, so only keywords are offered as
  keyword completions.

* *Argument values* - the live options for an argument, e.g. the loaded
  model numbers after "show model ".  A command declares the type of each
  argument with ``arg_types`` (see :class:`coot_commands.types.ArgType`);
  the completion engine asks that type for its current candidates.

Completion works one word at a time.  Given the words typed so far, we
collect every command word that could come next, then return the longest
common prefix (so a single Tab fills in as much as is unambiguous) plus the
full list of options to show the user when there is more than one.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

from coot_commands.registry import all_commands, normalise

# Leading verbs that are interchangeable across the validation family, so
# typing any of them at the start of a line offers the same commands (the
# regexes already accept them all - this just makes them equally
# discoverable by Tab).  "show", "list" and "find" are deliberately left out:
# they head their own curated families ("show model", "list maps",
# "find ligand") and must not be broadened with validation nouns.
_INTERCHANGEABLE_LEAD_VERBS = frozenset({"check", "validate"})


def _spec_candidates(spec: object) -> List[str]:
    """Resolve an ``arg_types`` entry to its current candidate strings.

    Accepts an :class:`~coot_commands.types.ArgType` (anything with a
    ``candidates()`` method), a zero-argument callable, or a plain iterable
    of strings.  Any failure yields no candidates rather than raising.
    """
    if spec is None:
        return []
    getter = getattr(spec, "candidates", None)
    if callable(getter):
        try:
            return [str(v) for v in getter()]
        except Exception:
            return []
    if callable(spec):
        try:
            return [str(v) for v in spec()]
        except Exception:
            return []
    if isinstance(spec, str):
        return [spec]
    try:
        return [str(v) for v in spec]
    except TypeError:
        return []


def _spec_label(spec: object, value: str) -> str:
    """Display label for an argument *value*, or the value itself.

    A type may define ``label(value)`` to decorate the value shown in the
    completion menu (e.g. adding a molecule name) without changing the value
    inserted into the command line.  Any failure falls back to the value.
    """
    labeller = getattr(spec, "label", None)
    if callable(labeller):
        try:
            return str(labeller(value))
        except Exception:
            return value
    return value


def _example_words(cmd, example: str) -> Tuple[List[str], List[bool], List]:
    """Split *example* into words, marking which are argument values.

    Returns ``(words, is_arg, arg_name)`` parallel lists.  A word is an
    argument value when it matches one of the command regex's named groups;
    the group name is recorded in *arg_name* so the caller can look the
    argument's type up in ``cmd.arg_types``.
    """
    ex = normalise(example)
    words = ex.split(" ") if ex else []
    is_arg = [False] * len(words)
    arg_name: List = [None] * len(words)

    match = cmd.regex.match(ex)
    if match:
        # Map each captured value back to its group name; the first group to
        # claim a value wins (values are distinct in practice).
        value_group = {}
        for name, value in match.groupdict().items():
            if value:
                value_group.setdefault(value, name)
        for i, word in enumerate(words):
            if word in value_group:
                is_arg[i] = True
                arg_name[i] = value_group[word]
    return words, is_arg, arg_name


def _prefix_matches(cmd, words_ex, is_arg, arg_name, fixed) -> bool:
    """Do the already-typed *fixed* words fit this example's earlier words?

    A keyword word normally matches exactly (case-insensitively).  The one
    exception is an interchangeable leading verb: if the example reads
    "check ramachandran" but the user typed "validate" (both in
    :data:`_INTERCHANGEABLE_LEAD_VERBS`), we accept it - provided swapping the
    verb into the example still matches the command's own regex - so the
    validation verbs surface each other's completions.  Other verbs stay
    literal, keeping the curated "show"/"list" menus intact.

    An argument word matches a typed value only when that value is one of the
    argument's current candidates - otherwise "perspective view" would let
    "view" be offered after any first word.  When an argument has no
    candidates to check against (a free-form argument, or a type whose
    values are momentarily unavailable) we accept the typed value.
    """
    for i, typed in enumerate(fixed):
        if is_arg[i]:
            options = _spec_candidates(cmd.arg_types.get(arg_name[i]))
            if options and typed.lower() not in {o.lower() for o in options}:
                return False
        elif words_ex[i].lower() != typed.lower():
            if not _is_verb_synonym(cmd, words_ex, i, typed):
                return False
    return True


def _is_verb_synonym(cmd, words_ex, i, typed) -> bool:
    """Is *typed* an accepted stand-in for this example's leading verb?

    Only the first word (``i == 0``) counts, and both the typed word and the
    example's word must be interchangeable validation verbs.  We then confirm
    the command's regex still matches the example with the verb swapped in,
    so a verb is only offered where the command genuinely accepts it.
    """
    if i != 0:
        return False
    if typed.lower() not in _INTERCHANGEABLE_LEAD_VERBS:
        return False
    if words_ex[i].lower() not in _INTERCHANGEABLE_LEAD_VERBS:
        return False
    swapped = [typed, *words_ex[1:]]
    return bool(cmd.regex.match(normalise(" ".join(swapped))))


def _sorted(candidates: Sequence[str]) -> List[str]:
    """Sort candidates numerically when they are all numbers, else by name."""
    items = list(candidates)
    if items and all(c.isdigit() for c in items):
        return sorted(items, key=int)
    return sorted(items)


def _common_prefix(strings: Sequence[str]) -> str:
    """Longest common (case-sensitive) prefix of *strings*."""
    if not strings:
        return ""
    lo = min(strings)
    hi = max(strings)
    n = 0
    for a, b in zip(lo, hi):
        if a != b:
            break
        n += 1
    return lo[:n]


def _gather(fixed: Sequence[str], partial: str) -> dict:
    """Candidate next words after the *fixed* words, filtered by *partial*.

    Returns a dict mapping each candidate value (what would be inserted) to
    a ``(label, is_arg)`` pair.  *label* is what the user sees in an
    ambiguous option list - it differs from the value only for argument
    types that decorate it (e.g. a model number gaining its molecule name).
    *is_arg* marks live argument values (model/map numbers etc.) apart from
    the fixed keyword vocabulary, so callers can treat the two differently.
    """
    pos = len(fixed)
    partial_low = partial.lower()
    found: dict = {}
    for cmd in all_commands():
        for example in cmd.examples:
            words_ex, is_arg, arg_name = _example_words(cmd, example)
            if len(words_ex) <= pos:
                continue
            if not _prefix_matches(cmd, words_ex, is_arg, arg_name, fixed):
                continue
            if is_arg[pos]:
                spec = cmd.arg_types.get(arg_name[pos])
                for value in _spec_candidates(spec):
                    if value.lower().startswith(partial_low):
                        found[value] = (_spec_label(spec, value), True)
            else:
                keyword = words_ex[pos]
                if keyword.lower().startswith(partial_low):
                    found.setdefault(keyword, (keyword, False))
    return found


def _is_complete(words: Sequence[str]) -> bool:
    """Do the joined *words* already form a complete, dispatchable command?

    Used to stop greedy multi-word completion at a valid command: we only
    pull in a further forced word when the line so far is not yet a command
    in its own right - so "check unmodelled" (invalid alone) extends to
    "check unmodelled blobs", but "ramachandran" (valid alone) is left as is.
    """
    text = normalise(" ".join(words))
    if not text:
        return False
    return any(m and m.end() == len(text)
               for m in (cmd.regex.match(text) for cmd in all_commands()))


def _extend_forced(fixed: List[str]) -> List[str]:
    """Pull in any further words that are forced, for multi-word completion.

    After an unambiguous keyword completion, keep appending words while the
    line is not yet a complete command and the next word is the single
    possible keyword - so a mandatory multi-word tail like "unmodelled
    blobs" fills in on one Tab.  Stops at a complete command, an ambiguous
    choice, or an argument value (which the user must supply themselves).
    """
    while not _is_complete(fixed):
        nxt = _gather(fixed, "")
        if len(nxt) != 1:
            break
        value, (_, is_arg) = next(iter(nxt.items()))
        if is_arg:
            break
        fixed = fixed + [value]
    return fixed


def complete(text: str) -> Tuple[str, List[str]]:
    """Complete the command line *text*.

    Returns ``(replacement, options)``:

    * *replacement* is the new command-line text - the input with the
      current word extended as far as is unambiguous.  It is empty when
      there is nothing to complete or no progress can be made.
    * *options* is the list of candidate words to show the user, non-empty
      only when the completion is ambiguous (more than one candidate).
    """
    raw = text or ""
    norm = normalise(raw)
    trailing_space = bool(raw) and raw[-1].isspace()
    words = norm.split(" ") if norm else []

    # Are we completing a fresh word (after a space, or on an empty line) or
    # extending the word under the cursor?
    if trailing_space or not words:
        fixed, partial = words, ""
    else:
        fixed, partial = words[:-1], words[-1]

    found = _gather(fixed, partial)
    if not found:
        return "", []

    ordered = _sorted(list(found))
    if len(ordered) == 1:
        # Unambiguous: fill it in, then pull in any forced multi-word tail,
        # and add a space ready for the next word.
        value = ordered[0]
        chosen = fixed + [value]
        if not found[value][1]:  # a keyword, so a tail may follow
            chosen = _extend_forced(chosen)
        return " ".join(chosen) + " ", []

    # Ambiguous: extend the word by the shared prefix and list the options.
    # The prefix is computed from the values (what gets inserted), while the
    # options shown to the user carry their display labels.
    prefix = _common_prefix(ordered)
    word = prefix if len(prefix) >= len(partial) else partial
    replacement = " ".join(fixed + [word]) if (fixed or word) else ""
    return replacement, [found[v][0] for v in ordered]
