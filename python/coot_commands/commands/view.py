# coot_commands/commands/view.py
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

"""Commands controlling the view: background, projection, motion, zoom.

These mirror the Draw > View menu.  Several are toggles (spin, rock,
fullscreen) - issuing the command again turns the effect off.
"""

from coot_commands.registry import command
from coot_commands.types import resolve_colour, as_float, CommandError

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "View"


@command(r"(?:set )?background (?:colou?r )?(?:to )?(?P<colour>\S+)",
         examples=["background black", "set background colour to white"],
         category=CATEGORY,
         notes="Colour names: black, white, grey, and the other named "
               "colours accepted by colour commands.")
def set_background(colour):
    """Set the background colour."""
    r, g, b = resolve_colour(colour)
    if coot is not None:
        coot.set_background_colour(r, g, b)
    return f"Background set to {colour.lower()}"


@command(r"spin(?: view)?",
         examples=["spin"],
         category=CATEGORY,
         notes="Toggles idle spinning; issue again to stop.")
def spin_view(**_):
    """Toggle spinning the view."""
    if coot is not None:
        coot.toggle_idle_spin_function()
    return "Toggled spin"


@command(r"rock(?: view)?",
         examples=["rock"],
         category=CATEGORY,
         notes="Toggles idle rocking; issue again to stop.")
def rock_view(**_):
    """Toggle rocking the view."""
    if coot is not None:
        coot.toggle_idle_rock_function()
    return "Toggled rock"


@command(r"(?:use )?(?P<mode>orthographic|perspective)(?: (?:view|projection))?",
         examples=["orthographic", "perspective view"],
         category=CATEGORY)
def set_projection(mode):
    """Switch between orthographic and perspective projection."""
    perspective = 1 if mode.lower() == "perspective" else 0
    if coot is not None:
        coot.set_use_perspective_projection(perspective)
    return f"Using {mode.lower()} projection"


@command(r"(?:toggle )?fullscreen",
         examples=["fullscreen"],
         category=CATEGORY,
         notes="Toggles fullscreen; issue again to leave fullscreen.")
def toggle_fullscreen(**_):
    """Toggle fullscreen mode."""
    if coot is not None:
        coot.fullscreen()
    return "Toggled fullscreen"


@command(r"(?:set )?zoom(?: (?:to|=))? (?P<factor>[\d.]+)",
         examples=["zoom to 30", "set zoom 50"],
         category=CATEGORY,
         notes="Larger numbers zoom out. Typical range ~10-100.")
def set_zoom(factor):
    """Set the view zoom factor."""
    f = as_float(factor, "zoom factor")
    if coot is not None:
        coot.set_zoom(f)
        coot.graphics_draw()
    return f"Zoom set to {f:g}"
