# coot_commands/commands/maps.py
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

"""Commands for map contouring and colouring."""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (resolve_map, require_map, resolve_colour,
                                 as_float, ArgType, ACTIVE_MAP_NOTE,
                                 CommandError)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Maps"


@command(r"contour (?:map (?P<map>\S+) )?(?:to|at) (?P<level>[\d.]+) ?(?:sigma|rmsd|rms)",
         examples=["contour map 1 to 1.5 sigma", "contour to 1.2 rmsd"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP},
         notes="Sets the contour level in sigma (map RMSD units). "
               + ACTIVE_MAP_NOTE)
def contour_sigma(level: str, map: Optional[str] = None) -> str:
    """Set a map's contour level in sigma."""
    imol = resolve_map(map)
    lvl = as_float(level, "contour level")
    if coot is not None:
        coot.set_contour_level_in_sigma(imol, lvl)
    return f"Contoured map {imol} at {lvl:g} sigma"


@command(r"contour (?:map (?P<map>\S+) )?(?:to|at) (?P<level>[\d.]+)",
         examples=["contour map 1 to 0.35", "contour to 0.3"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP},
         notes="Sets the absolute contour level (map units). " + ACTIVE_MAP_NOTE
               + " Add \"sigma\" to contour in RMSD units instead.")
def contour_absolute(level: str, map: Optional[str] = None) -> str:
    """Set a map's absolute contour level."""
    imol = resolve_map(map)
    lvl = as_float(level, "contour level")
    if coot is not None:
        coot.set_contour_level_absolute(imol, lvl)
    return f"Contoured map {imol} at {lvl:g}"


@command(r"(?:(?:make )?(?:a )?smoother copy(?: of)?|smooth(?:en)?)"
         r"(?: map)?(?: (?P<map>(?!by\b|x\b)\S+))?"
         r"(?:(?: by| x) (?P<factor>[\d.]+))?$",
         examples=["make a smoother copy of map 1", "smooth map 1 by 2",
                   "smooth map"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP},
         notes="Interpolates the map onto a finer grid, making a new map "
               "molecule and leaving the original alone. The factor "
               "multiplies the existing grid sampling (the default is 1.25, "
               "as in Calculate -> Map Tools); larger factors are smoother "
               "and use more memory. Unlike 'set map sampling rate' this "
               "acts on a map that is already loaded, but it interpolates "
               "the existing grid rather than re-transforming the structure "
               "factors, so it adds no real detail. " + ACTIVE_MAP_NOTE)
def smooth_map(map: Optional[str] = None, factor: Optional[str] = None) -> str:
    """Make a smoother (finer-sampled) copy of a map."""
    imol = require_map(resolve_map(map))
    mult = 1.25 if factor is None else as_float(factor, "smoothing factor")
    if mult <= 1.0:
        raise CommandError(f"the smoothing factor must be more than 1, got {mult:g}")
    if coot is None:
        return f"Would smooth map {imol} by a factor of {mult:g}"
    imol_new = coot.smooth_map(imol, mult)
    if imol_new < 0:
        raise CommandError(f"could not make a smoother copy of map {imol}")
    return f"Made map {imol_new}: map {imol} smoothed by a factor of {mult:g}"


@command(r"colou?r map(?: (?P<map>\S+))? (?P<colour>\S+)",
         examples=["colour map 1 blue", "colour map cyan"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP, "colour": ArgType.COLOUR},
         notes=ACTIVE_MAP_NOTE)
def colour_map(colour: str, map: Optional[str] = None) -> str:
    """Set a map's colour."""
    imol = resolve_map(map)
    r, g, b = resolve_colour(colour)
    if coot is not None:
        coot.set_map_colour(imol, r, g, b)
    return f"Coloured map {imol} {colour.lower()}"


@command(r"(?:make (?:map )?)?(?:map (?P<map>\S+) )?(?:is )?(?:a )?difference map$",
         examples=["map 1 is a difference map", "make map 2 a difference map"],
         category=CATEGORY,
         arg_types={"map": ArgType.MAP},
         notes="Marks the map as a difference map (green/red, contoured "
               "either side of zero). " + ACTIVE_MAP_NOTE)
def set_difference_map(map: Optional[str] = None) -> str:
    """Mark a map as a difference map."""
    imol = resolve_map(map)
    if coot is not None:
        coot.set_map_is_difference_map(imol, 1)
    return f"Map {imol} marked as a difference map"
