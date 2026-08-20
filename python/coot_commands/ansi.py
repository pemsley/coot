# coot_commands/ansi.py
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

"""Small helpers for colouring command output with ANSI escape codes.

The Command tab renders a subset of ANSI SGR codes (24-bit truecolour, the
eight standard colours and bold - see ``command_output_append`` in
``src/vte.cc``), and the Assistant's VTE terminal renders them natively.  A
command handler can therefore wrap text in these helpers to add colour; a
sink that does not understand the codes still shows the plain text, so this
is safe to use unconditionally.
"""

from __future__ import annotations

from typing import Sequence, Union

RESET = "\033[0m"

# Either 0-1 floats (as Coot's colour APIs return) or 0-255 ints.
RGB = Sequence[Union[int, float]]


def _to_255(component: Union[int, float]) -> int:
    """Coerce one colour component to an int in 0-255, clamped."""
    value = round(component * 255) if component <= 1.0 else round(component)
    return max(0, min(255, int(value)))


def truecolour(rgb: RGB) -> str:
    """The SGR sequence that sets the 24-bit foreground colour."""
    r, g, b = (_to_255(c) for c in tuple(rgb)[:3])
    return f"\033[38;2;{r};{g};{b}m"


def colourise(text: str, rgb: RGB) -> str:
    """Wrap *text* so it renders in the given colour, then resets."""
    return f"{truecolour(rgb)}{text}{RESET}"


def swatch(rgb: RGB, glyph: str = "██") -> str:
    """A small coloured block (default two full-block characters)."""
    return colourise(glyph, rgb)
