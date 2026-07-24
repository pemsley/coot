# coot_commands/speech.py
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

"""Turn dictated speech into the canonical text the command patterns expect.

macOS Dictation (and other speech-to-text) types its transcription straight
into the focused field, so a spoken command arrives as ordinary text - but
worded the way people speak rather than the way the command regexes are
written.  :func:`from_speech` rewrites those spoken forms into the canonical
tokens, and :func:`~coot_commands.registry.dispatch` runs it on every input,
so "superpose model zero onto model one" reaches the same handler as
"superpose model 0 onto model 1".

It only ever rewrites number words, "point"/"minus" and spoken separators; a
command that is already typed with digits passes through unchanged, so this is
safe to run on all input, not just dictated input.  No command keyword is a
number word, and the number words themselves are spelled distinctly from their
homophones ("two" not "to", "four" not "for"), so the rewrite does not clash
with the vocabulary.
"""

from __future__ import annotations

import re
from typing import List, Tuple

_UNITS = {
    "zero": 0, "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
    "six": 6, "seven": 7, "eight": 8, "nine": 9,
}
_TEENS = {
    "ten": 10, "eleven": 11, "twelve": 12, "thirteen": 13, "fourteen": 14,
    "fifteen": 15, "sixteen": 16, "seventeen": 17, "eighteen": 18,
    "nineteen": 19,
}
_TENS = {
    "twenty": 20, "thirty": 30, "forty": 40, "fifty": 50, "sixty": 60,
    "seventy": 70, "eighty": 80, "ninety": 90,
}
_NUMBER_WORDS = {**_UNITS, **_TEENS, **_TENS}


def _is_number_word(word: str) -> bool:
    w = word.lower()
    return w in _NUMBER_WORDS or w in ("hundred", "thousand")


def _parse_number(tokens: List[str], start: int) -> Tuple[int, int]:
    """Parse a run of spoken number words at *start*.

    Returns ``(value, count)`` where *count* is how many tokens were consumed
    (0 if *start* is not a number word).  Follows ordinary English number
    grammar - "forty five" -> 45, "one hundred twenty" -> 120 - but stops
    rather than merging two bare units ("four five" -> 4, then 5), so a number
    read out digit by digit is not silently turned into its sum.
    """
    total = 0        # sum of scale-completed groups (… thousand)
    current = 0      # the group being built
    unit_open = False  # a units/teens value has been placed in this group
    count = 0
    saw = False
    n = len(tokens)

    while start + count < n:
        w = tokens[start + count].lower()
        if w in _UNITS:
            # A unit fits after a tens word ("forty" "five") or a hundred
            # boundary, but never straight after another unit/teen - so a
            # number read out digit by digit ("four" "five") is not merged.
            if unit_open:
                break
            current += _NUMBER_WORDS[w]
            unit_open = True
        elif w in _TEENS:
            # A teen (10-19) only starts a group or follows a hundred.
            if unit_open or (current % 100) != 0:
                break
            current += _NUMBER_WORDS[w]
            unit_open = True
        elif w in _TENS:
            if unit_open or (current % 100) != 0:
                break
            current += _NUMBER_WORDS[w]
        elif w == "hundred" and saw:
            current = (current or 1) * 100
            unit_open = False
        elif w == "thousand" and saw:
            total += (current or 1) * 1000
            current = 0
            unit_open = False
        elif w == "and" and saw and start + count + 1 < n \
                and _is_number_word(tokens[start + count + 1]):
            pass  # spoken connector, e.g. "one hundred and five"
        else:
            break
        saw = True
        count += 1

    if not saw:
        return (0, 0)
    return (total + current, count)


_POINT = re.compile(r"(\d) (?:point|dot) (\d)", re.IGNORECASE)
_NEGATIVE = re.compile(r"\b(?:minus|negative|dash) (\d)", re.IGNORECASE)
_SPACED_SLASH = re.compile(r" ?/ ?")


def from_speech(text: str) -> str:
    """Rewrite dictated *text* into canonical command text.

    Idempotent on already-canonical (digit) input.
    """
    if not text:
        return text or ""

    collapsed = re.sub(r"\s+", " ", text.strip())
    if not collapsed:
        return ""

    tokens = collapsed.split(" ")
    out: List[str] = []
    i = 0
    while i < len(tokens):
        value, consumed = _parse_number(tokens, i)
        if consumed:
            out.append(str(value))
            i += consumed
        else:
            out.append(tokens[i])
            i += 1
    result = " ".join(out)

    # "one point five" -> "1 point 5" -> "1.5"
    result = _POINT.sub(r"\1.\2", result)
    # "minus 5" / "negative 5" -> "-5" (for negative residue numbers)
    result = _NEGATIVE.sub(r"-\1", result)
    # spoken or spaced chain/residue separator -> a bare slash ("A / 45" -> "A/45")
    result = result.replace(" slash ", "/").replace(" stroke ", "/")
    result = _SPACED_SLASH.sub("/", result)
    return result
