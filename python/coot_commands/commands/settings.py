# coot_commands/commands/settings.py
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

"""Adjust and read back Coot settings.

Each setting here comes as a pair: a ``set ...`` command that changes it
and a matching ``get ...`` (also spelled "show ..." or "what is ...")
that reports the current value.  The word *set* is reserved for these
settings commands, so an action like refining or deleting never starts
with it.

Setting values are read straight from the Coot API's own getters, so a
``get`` always reflects the live state - including changes made through
the GUI rather than these commands.
"""

from __future__ import annotations

from typing import Optional, Any

from coot_commands.registry import command
from coot_commands.types import (as_int, as_float, resolve_map, require_map,
                                 resolve_model, active_map, first_difference_map,
                                 ArgType, ACTIVE_MAP_NOTE, CommandError)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Settings"

# Shared regex fragments.  ``_GET`` is the optional prefix that turns a
# setting's noun into a query: "bond thickness", "get bond thickness",
# "show the bond thickness", "what's the bond thickness" all match the same
# handler.  ``_TO`` is the connector between a setting and its new value,
# accepting "set X to N", "set X = N" and a bare "set X N" alike (mirroring
# the zoom command in view.py).
_GET = r"(?:(?:get|show|what(?:'s| is))(?: the)? )?"
_TO = r"(?: (?:to|=))? "


def _read(getter: str) -> Optional[Any]:
    """Return the current value from a Coot getter, or ``None`` if unavailable.

    Looks the getter up by name so a build that lacks it (or a test run with
    no Coot module) degrades to a clear "unavailable" message rather than an
    AttributeError.
    """
    if coot is None:
        return None
    fn = getattr(coot, getter, None)
    return fn() if fn is not None else None


def _draw() -> None:
    """Redraw the graphics if the API is present (a no-op otherwise)."""
    if coot is not None:
        coot.graphics_draw()


# ---------------------------------------------------------------------------
#   Bond thickness
# ---------------------------------------------------------------------------

@command(r"set bond thickness" + _TO + r"(?P<value>\d+)",
         examples=["set bond thickness to 4", "set bond thickness 3"],
         category=CATEGORY,
         notes="Line/stick thickness for model bonds, in pixels. "
               "A whole number; the default is 5.")
def set_bond_thickness(value: str) -> str:
    """Set the default bond (stick) thickness."""
    t = as_int(value, "bond thickness")
    if coot is not None:
        coot.set_default_bond_thickness(t)
    _draw()
    return f"Set bond thickness to {t}"


@command(_GET + r"bond thickness\??",
         examples=["get bond thickness", "what is the bond thickness"],
         category=CATEGORY)
def get_bond_thickness(**_: Optional[str]) -> str:
    """Report the current default bond thickness."""
    value = _read("get_default_bond_thickness")
    if value is None:
        return "Bond thickness is unavailable (no Coot API)"
    return f"Bond thickness is {value}"


# ---------------------------------------------------------------------------
#   Font size
# ---------------------------------------------------------------------------

@command(r"set font size" + _TO + r"(?P<value>\d+)",
         examples=["set font size to 2", "set font size 3"],
         category=CATEGORY,
         notes="Size of on-screen labels: 1 (small), 2 (medium) or 3 (large).")
def set_font_size(value: str) -> str:
    """Set the on-screen label font size."""
    size = as_int(value, "font size")
    if coot is not None:
        coot.set_font_size(size)
    _draw()
    return f"Set font size to {size}"


@command(_GET + r"font size\??",
         examples=["get font size", "what is the font size"],
         category=CATEGORY)
def get_font_size(**_: Optional[str]) -> str:
    """Report the current on-screen label font size."""
    value = _read("get_font_size")
    if value is None:
        return "Font size is unavailable (no Coot API)"
    return f"Font size is {value}"


# ---------------------------------------------------------------------------
#   Map display radius
# ---------------------------------------------------------------------------

@command(r"set map radius" + _TO + r"(?P<value>[\d.]+)",
         examples=["set map radius to 20", "set map radius 15"],
         category=CATEGORY,
         notes="Radius, in Angstroms, of the sphere of density drawn around "
               "the screen centre. Larger radii are slower to contour.")
def set_map_radius(value: str) -> str:
    """Set the map display radius (Angstroms)."""
    r = as_float(value, "map radius")
    if coot is not None:
        coot.set_map_radius(r)
    _draw()
    return f"Set map radius to {r:g} A"


@command(_GET + r"map radius\??",
         examples=["get map radius", "what is the map radius"],
         category=CATEGORY)
def get_map_radius(**_: Optional[str]) -> str:
    """Report the current map display radius."""
    value = _read("get_map_radius")
    if value is None:
        return "Map radius is unavailable (no Coot API)"
    return f"Map radius is {value:g} A"


# ---------------------------------------------------------------------------
#   Map sampling rate
# ---------------------------------------------------------------------------

@command(r"set map sampling(?: rate)?" + _TO + r"(?P<value>[\d.]+)",
         examples=["set map sampling rate to 1.8", "set map sampling 2"],
         category=CATEGORY,
         notes="Finer sampling (higher rate) makes smoother maps at the cost "
               "of memory. Applies to maps read after it is set; typical "
               "values are 1.5-2.5.")
def set_map_sampling_rate(value: str) -> str:
    """Set the map sampling rate for maps read from now on."""
    rate = as_float(value, "map sampling rate")
    if coot is not None:
        coot.set_map_sampling_rate(rate)
    return f"Set map sampling rate to {rate:g}"


@command(_GET + r"map sampling(?: rate)?\??",
         examples=["get map sampling rate", "what is the map sampling rate"],
         category=CATEGORY)
def get_map_sampling_rate(**_: Optional[str]) -> str:
    """Report the current map sampling rate."""
    value = _read("get_map_sampling_rate")
    if value is None:
        return "Map sampling rate is unavailable (no Coot API)"
    return f"Map sampling rate is {value:g}"


# ---------------------------------------------------------------------------
#   Contour step (iso-level increment)
# ---------------------------------------------------------------------------

@command(r"set contour step" + _TO + r"(?P<value>[\d.]+)",
         examples=["set contour step to 0.1", "set contour step 0.05"],
         category=CATEGORY,
         notes="How much a scroll changes the contour level of a normal "
               "(non-difference) map. See 'set difference map contour step' "
               "for difference maps.")
def set_contour_step(value: str) -> str:
    """Set the contour-level scroll step for normal maps."""
    step = as_float(value, "contour step")
    if coot is not None:
        coot.set_iso_level_increment(step)
    return f"Set contour step to {step:g}"


@command(_GET + r"contour step\??",
         examples=["get contour step", "what is the contour step"],
         category=CATEGORY)
def get_contour_step(**_: Optional[str]) -> str:
    """Report the contour-level scroll step for normal maps."""
    value = _read("get_iso_level_increment")
    if value is None:
        return "Contour step is unavailable (no Coot API)"
    return f"Contour step is {value:g}"


@command(r"set (?:diff(?:erence)? map|difference) contour step" + _TO
         + r"(?P<value>[\d.]+)",
         examples=["set difference map contour step to 0.1",
                   "set diff map contour step 0.05"],
         category=CATEGORY,
         notes="How much a scroll changes the contour level of a difference "
               "map.")
def set_diff_contour_step(value: str) -> str:
    """Set the contour-level scroll step for difference maps."""
    step = as_float(value, "difference map contour step")
    if coot is not None:
        coot.set_diff_map_iso_level_increment(step)
    return f"Set difference map contour step to {step:g}"


@command(_GET + r"(?:diff(?:erence)? map|difference) contour step\??",
         examples=["get difference map contour step",
                   "what is the difference map contour step"],
         category=CATEGORY)
def get_diff_contour_step(**_: Optional[str]) -> str:
    """Report the contour-level scroll step for difference maps."""
    value = _read("get_diff_map_iso_level_increment")
    if value is None:
        return "Difference map contour step is unavailable (no Coot API)"
    return f"Difference map contour step is {value:g}"


# ---------------------------------------------------------------------------
#   Per-map contour level
# ---------------------------------------------------------------------------

_MAP_OF = r"(?:(?: of| for)? map (?P<map>\S+))?"


@command(r"set contour level" + _MAP_OF + _TO
         + r"(?P<value>-?[\d.]+)(?P<sigma> sigma)?",
         examples=["set contour level to 1.5 sigma",
                   "set contour level of map 1 to 0.3",
                   "set contour level of map 1 to 2 sigma"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP},
         notes="Sets the contour level of a map, in absolute units by "
               "default or in sigma when the value ends with 'sigma'. "
               + ACTIVE_MAP_NOTE)
def set_contour_level(value: str, map: Optional[str] = None,
                      sigma: Optional[str] = None) -> str:
    """Set the contour level of a map (absolute, or in sigma)."""
    imol = require_map(resolve_map(map))
    level = as_float(value, "contour level")
    if coot is not None:
        if sigma:
            coot.set_contour_level_in_sigma(imol, level)
        else:
            coot.set_contour_level_absolute(imol, level)
    _draw()
    unit = " sigma" if sigma else ""
    return f"Set contour level of map {imol} to {level:g}{unit}"


@command(_GET + r"contour level" + _MAP_OF + r"\??",
         examples=["get contour level", "what is the contour level of map 1"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP},
         notes=ACTIVE_MAP_NOTE)
def get_contour_level(map: Optional[str] = None) -> str:
    """Report the current contour level of a map."""
    imol = require_map(resolve_map(map))
    if coot is None:
        return "Contour level is unavailable (no Coot API)"
    absolute = coot.get_contour_level_absolute(imol)
    sigma = coot.get_contour_level_in_sigma(imol)
    return f"Contour level of map {imol} is {absolute:g} ({sigma:g} sigma)"


# ---------------------------------------------------------------------------
#   Updating (sfcalc) maps
# ---------------------------------------------------------------------------
#
# Coot's "updating maps" feature recomputes 2Fo-Fc and Fo-Fc maps from the
# reflection data as the model is edited, so density follows the atoms as you
# build.  Turning it on needs three molecules: the model, a map that carries
# the reflection data (the 2Fo-Fc map read from the MTZ - i.e. the refinement
# map), and an Fo-Fc difference map to write updates into.  We take the model
# from the active model, the data map from the refinement map, and the first
# loaded difference map, so a plain "set updating maps on" works for the usual
# one-model/one-dataset session; explain what is missing otherwise.

_UPDATING_NOTE = (
    "Coot's auto-updating maps recompute the 2Fo-Fc and difference maps from "
    "the reflection data as you edit the model, so density follows the atoms. "
    "Needs the refinement map (a map with reflection data attached) and a "
    "difference map to be loaded - 'load tutorial' or opening an MTZ provides "
    "both. Only one updating-maps session can run at a time.")


def _start_updating_maps() -> str:
    """Wire up auto-updating maps for the active model, or raise."""
    imol_model = resolve_model()
    imol_data = active_map()
    if coot is not None and coot.map_is_difference_map(imol_data) == 1:
        raise CommandError(
            f"map {imol_data} is a difference map; the refinement map must be "
            "the 2Fo-Fc map that carries the reflection data")
    imol_diff = first_difference_map()
    if imol_diff is None:
        raise CommandError(
            "no difference map loaded - updating maps needs a difference (Fo-Fc) "
            "map to update; open your data (e.g. 'load tutorial') first")
    if coot is not None:
        coot.set_auto_updating_sfcalc_genmap(imol_model, imol_data, imol_diff)
    return (f"Turned on updating maps: model {imol_model}, data from map "
            f"{imol_data}, updating difference map {imol_diff}")


def _stop_updating_maps() -> str:
    """Stop watching the active model for updating maps."""
    imol_model = resolve_model()
    if coot is not None and hasattr(coot, "stop_updating_molecule"):
        coot.stop_updating_molecule(imol_model)
    return f"Turned off updating maps for model {imol_model}"


@command(r"(?:set )?updating maps?(?: (?P<state>on|off))?",
         examples=["set updating maps on", "set updating maps off"],
         category=CATEGORY,
         arg_types={"state": ("on", "off")},
         notes="Turns Coot's auto-updating maps on or off. " + _UPDATING_NOTE)
def set_updating_maps(state: Optional[str] = None) -> str:
    """Turn Coot's auto-updating (sfcalc) maps on or off."""
    if state is not None and state.lower() == "off":
        return _stop_updating_maps()
    return _start_updating_maps()


@command(r"stop updating maps?",
         examples=["stop updating maps"],
         category=CATEGORY,
         notes="Turns off auto-updating maps (same as 'set updating maps "
               "off'). " + _UPDATING_NOTE)
def stop_updating_maps(**_: Optional[str]) -> str:
    """Stop Coot's auto-updating (sfcalc) maps."""
    return _stop_updating_maps()
