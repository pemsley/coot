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

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (RES_SPEC, resolve_colour, resolve_model,
                                 resolve_residue, as_float, first_present,
                                 CommandError, ArgType, ACTIVE_MODEL_NOTE)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "View"


@command(r"(?:set )?background (?:colou?r )?(?:to )?(?P<colour>\S+)",
         examples=["background black", "set background colour to white"],
         category=CATEGORY,
         arg_types={"colour": ArgType.COLOUR},
         notes="Colour names: black, white, grey, and the other named "
               "colours accepted by colour commands.")
def set_background(colour: str) -> str:
    """Set the background colour."""
    r, g, b = resolve_colour(colour)
    if coot is not None:
        coot.set_background_colour(r, g, b)
    return f"Set background to {colour.lower()}"


@command(r"spin(?: view)?",
         examples=["spin"],
         category=CATEGORY,
         notes="Toggles idle spinning; issue again to stop.")
def spin_view(**_: Optional[str]) -> str:
    """Toggle spinning the view."""
    if coot is not None:
        coot.toggle_idle_spin_function()
    return "Toggled spin"


@command(r"rock(?: view)?",
         examples=["rock"],
         category=CATEGORY,
         notes="Toggles idle rocking; issue again to stop.")
def rock_view(**_: Optional[str]) -> str:
    """Toggle rocking the view."""
    if coot is not None:
        coot.toggle_idle_rock_function()
    return "Toggled rock"


@command(r"(?:use )?(?P<mode>orthographic|perspective)(?: (?:view|projection))?",
         examples=["orthographic", "perspective view"],
         category=CATEGORY,
         arg_types={"mode": ("orthographic", "perspective")})
def set_projection(mode: str) -> str:
    """Switch between orthographic and perspective projection."""
    perspective = 1 if mode.lower() == "perspective" else 0
    if coot is not None:
        coot.set_use_perspective_projection(perspective)
    return f"Using {mode.lower()} projection"


@command(r"(?:toggle )?fullscreen",
         examples=["fullscreen"],
         category=CATEGORY,
         notes="Toggles fullscreen; issue again to leave fullscreen.")
def toggle_fullscreen(**_: Optional[str]) -> str:
    """Toggle fullscreen mode."""
    if coot is not None:
        coot.fullscreen()
    return "Toggled fullscreen"


@command(r"(?:set )?zoom(?: (?:to|=))? (?P<factor>[\d.]+)",
         examples=["zoom to 30", "set zoom 50"],
         category=CATEGORY,
         notes="Larger numbers zoom out. Typical range ~10-100.")
def set_zoom(factor: str) -> str:
    """Set the view zoom factor."""
    f = as_float(factor, "zoom factor")
    if coot is not None:
        coot.set_zoom(f)
        coot.graphics_draw()
    return f"Set zoom to {f:g}"


# ---------------------------------------------------------------------------
#   Zooming to fit a selection
# ---------------------------------------------------------------------------
#
# "zoom to 30" above sets the number; these set the *view*, working out the
# number from how big the thing is (coot_commands.focus). They are ordered
# residue, chain, molecule so the most specific pattern is tried first.


@command(r"(?:zoom|view|frame) (?:to |on )?(?:fit )?"
         r"(?:residue |res )?" + RES_SPEC +
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["zoom to A/45", "zoom to residue A 45", "zoom on A/45"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Centres on the residue and picks a zoom that fits it with "
               "enough of its surroundings to judge it. " + ACTIVE_MODEL_NOTE)
def zoom_to_residue(chain: str, resno: str, model: Optional[str] = None) -> str:
    """Zoom to fit one residue."""
    from coot_commands import focus
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    zoom = focus.show_residue(imol, chain_id, res, ins)
    if zoom is None:
        return f"Could not zoom to {chain_id}/{res} of model {imol}"
    return f"Zoomed to {chain_id}/{res} of model {imol} (zoom {zoom:.0f})"


@command(r"(?:zoom|view|frame) (?:to |on )?(?:fit )?chain (?P<chain>[A-Za-z0-9])"
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["zoom to chain A", "zoom on chain B"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Centres on the chain and zooms out far enough to see all of "
               "it. " + ACTIVE_MODEL_NOTE)
def zoom_to_chain(chain: str, model: Optional[str] = None) -> str:
    """Zoom to fit a whole chain."""
    from coot_commands import focus
    imol = resolve_model(model)
    chain_id = chain.upper()
    zoom = focus.show_chain(imol, chain_id)
    if zoom is None:
        return f"No chain {chain_id} in model {imol}"
    return f"Zoomed to chain {chain_id} of model {imol} (zoom {zoom:.0f})"


@command(r"(?:zoom|view|frame) (?:to |on )?fit"
         r"(?: (?:of |in |on )?model (?P<model>\S+))?$"
         r"|(?:zoom|view|frame) (?:to |on )?(?:the )?"
         r"(?:whole |entire )?model(?: (?P<model2>\S+))?$",
         examples=["zoom to fit", "zoom to fit model 0", "zoom to model 1"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "model2": ArgType.MODEL},
         notes="Zooms out until the whole model is in view, like a "
               "'fit to screen'. " + ACTIVE_MODEL_NOTE)
def zoom_to_fit(model: Optional[str] = None,
                model2: Optional[str] = None) -> str:
    """Zoom to fit a whole model."""
    from coot_commands import focus
    imol = resolve_model(first_present(model, model2))
    zoom = focus.show_molecule(imol)
    if zoom is None:
        return f"Could not measure model {imol} to zoom to it"
    return f"Zoomed to fit model {imol} (zoom {zoom:.0f})"


@command(r"(?:set )?(?:auto[- ]?zoom|auto[- ]?focus)(?: (?P<state>on|off))?",
         examples=["autozoom on", "autozoom off", "auto zoom on"],
         category=CATEGORY,
         arg_types={"state": ("on", "off")},
         notes="When on (the default), the view follows whatever a command "
               "acts on: centring on the residue and picking a zoom that "
               "fits it, so you can see the assistant work. Turn it off to "
               "keep the view where you put it; the explicit 'zoom to ...' "
               "commands still work either way.")
def set_autozoom(state: Optional[str] = None) -> str:
    """Turn automatic view-following on or off."""
    from coot_commands import focus
    if state is None:
        return ("Autozoom is "
                + ("on" if focus.autozoom_enabled() else "off"))
    enabled = state.lower() == "on"
    focus.set_autozoom(enabled)
    return "Turned autozoom " + ("on" if enabled else "off")
