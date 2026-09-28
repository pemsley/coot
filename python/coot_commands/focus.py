# coot_commands/focus.py
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

"""Put what is being worked on on the screen, at a sensible size.

Centring is not enough.  Coot opens at a zoom that shows a whole domain, so
a command that centres on one residue leaves it a speck in the middle of the
view; and a zoom tight enough for a side chain is useless for a chain.  The
right zoom depends on how big the thing being worked on is, so this module
computes it from the selection's actual extent rather than guessing.

The arithmetic is exact, not a fudge factor.  Coot's projection is

    glOrtho(-0.3*zoom*aspect, 0.3*zoom*aspect, -0.3*zoom, 0.3*zoom, ...)

(``draw.cc``), so the view is :data:`VIEW_HEIGHT_PER_ZOOM` * ``zoom``
angstroms tall.  Fitting something *S* angstroms across therefore needs a
zoom of ``S / 0.6``, and everything below is that one relationship plus a
padding factor so the selection does not touch the edges.

Who calls this
--------------
:func:`follow_residue` is called from
:func:`coot_commands.types.resolve_residue` - the one function every
per-residue command uses to work out what it is acting on.  Hooking there
means the view follows the assistant automatically as it works, and a new
per-residue command gets the behaviour without having to remember to ask for
it.  Turn it off with the ``autozoom`` command or ``COOT_AUTOZOOM=0``; the
explicit ``zoom to ...`` commands still work when it is off, since those are
a request rather than a side effect.
"""

from __future__ import annotations

import os
from typing import List, Optional, Sequence, Tuple

from coot_commands import scoring

try:
    import coot
except ImportError:
    coot = None


# Angstroms of view height per unit of zoom - from the glOrtho call in draw.cc,
# whose vertical half-extent is 0.3 * zoom.
VIEW_HEIGHT_PER_ZOOM = 0.6

# How much bigger than the selection the view should be. This is not one
# number, because the padding is doing two different jobs. Framing a residue,
# the extra space IS the point - you judge a side chain by its neighbours, its
# density and what it might be clashing with, none of which are inside its own
# bounding box. Framing a whole chain or molecule, there is no wider context
# to show and the padding is just a margin so the edges are not flush.
PADDING = 1.8            # a residue: show the neighbourhood
PADDING_SPAN = 1.4       # a residue range: a little context either end
PADDING_WHOLE = 1.15     # a chain or molecule: just a margin

# Nothing is ever fitted more tightly than this many angstroms across. A
# glycine is barely 3 A wide, and a view that tight shows the residue and
# nothing to judge it against - you want its neighbours, the density, and
# whatever it might be clashing with.
MIN_EXTENT = 8.0

# Zoom limits. Below ~8 the near clipping plane starts cutting into the
# residue. The upper limit has to clear the largest thing anyone would ask to
# fit - a 500 A capsid needs ~960 - so it is a guard against a degenerate
# extent rather than a working range.
MIN_ZOOM = 8.0
MAX_ZOOM = 1000.0

# How many residues to sample when measuring a span or a chain. Reading every
# residue's atoms to measure a 300-residue chain costs far more than the
# answer is worth, and a dozen samples plus PADDING bounds it well enough.
MAX_SAMPLES = 12

Point = Tuple[float, float, float]


# ---------------------------------------------------------------------------
#   The setting
# ---------------------------------------------------------------------------

_autozoom: Optional[bool] = None


def autozoom_enabled() -> bool:
    """Is the view following what commands act on?

    Defaults to on - the point of the feature is that you do not have to ask -
    unless ``COOT_AUTOZOOM`` says otherwise (``0``, ``off`` or ``false``).
    """
    global _autozoom
    if _autozoom is None:
        raw = os.environ.get("COOT_AUTOZOOM", "").strip().lower()
        _autozoom = raw not in ("0", "off", "false", "no")
    return _autozoom


def set_autozoom(enabled: bool) -> None:
    """Turn view-following on or off for the rest of the session."""
    global _autozoom
    _autozoom = bool(enabled)


# ---------------------------------------------------------------------------
#   Zoom from extent
# ---------------------------------------------------------------------------


def points_extent(points: Sequence[Point]) -> float:
    """The largest span of *points* along any axis, in angstroms.

    Axis-aligned rather than a true diameter: the view is axis-aligned too,
    and the padding covers the difference.
    """
    if not points:
        return 0.0
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    return max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def centre_of(points: Sequence[Point]) -> Optional[Point]:
    """The midpoint of the bounding box of *points*, or ``None`` if empty.

    The bounding-box centre, not the centroid: we are framing a view, so what
    matters is that the extremes are equally far from the middle, not where
    the mass happens to sit.
    """
    if not points:
        return None
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    return ((max(xs) + min(xs)) / 2.0,
            (max(ys) + min(ys)) / 2.0,
            (max(zs) + min(zs)) / 2.0)


def zoom_for_extent(extent: float, padding: float = PADDING,
                    min_extent: float = MIN_EXTENT) -> float:
    """The zoom factor that fits *extent* angstroms in the view, with padding."""
    height = max(extent, min_extent) * padding
    return max(MIN_ZOOM, min(MAX_ZOOM, height / VIEW_HEIGHT_PER_ZOOM))


# ---------------------------------------------------------------------------
#   Reading coordinates
# ---------------------------------------------------------------------------


def residue_points(imol: int, chain_id: str, resno: int,
                   ins_code: str = "") -> List[Point]:
    """Every atom position in one residue."""
    points: List[Point] = []
    for atom in scoring.residue_atoms(imol, chain_id, resno, ins_code):
        try:
            points.append((float(atom[2][0]), float(atom[2][1]),
                           float(atom[2][2])))
        except (IndexError, TypeError, ValueError):
            continue
    return points


def _sampled(values: Sequence[int], n: int = MAX_SAMPLES) -> List[int]:
    """At most *n* of *values*, evenly spread, always keeping both ends."""
    items = list(values)
    if len(items) <= n:
        return items
    step = (len(items) - 1) / float(n - 1)
    picked = {int(round(i * step)) for i in range(n)}
    picked.add(len(items) - 1)
    return [items[i] for i in sorted(picked)]


def span_points(imol: int, chain_id: str, res1: int, res2: int) -> List[Point]:
    """Atom positions sampled across a residue range."""
    lo, hi = (res1, res2) if res1 <= res2 else (res2, res1)
    points: List[Point] = []
    for resno in _sampled(range(lo, hi + 1)):
        points.extend(residue_points(imol, chain_id, resno))
    return points


def chain_residue_numbers(imol: int, chain_id: str) -> List[int]:
    """The residue numbers of a chain, in order, or ``[]``."""
    if coot is None:
        return []
    try:
        n = coot.chain_n_residues(chain_id, imol)
        return [coot.seqnum_from_serial_number(imol, chain_id, i)
                for i in range(n)]
    except Exception:  # noqa: BLE001 - a missing chain must not raise here
        return []


def chain_points(imol: int, chain_id: str) -> List[Point]:
    """Atom positions sampled along a whole chain."""
    points: List[Point] = []
    for resno in _sampled(chain_residue_numbers(imol, chain_id)):
        points.extend(residue_points(imol, chain_id, resno))
    return points


def molecule_points(imol: int) -> List[Point]:
    """Atom positions sampled across every chain of a molecule."""
    if coot is None:
        return []
    points: List[Point] = []
    try:
        n_chains = coot.n_chains(imol)
    except Exception:  # noqa: BLE001
        return points
    for i in range(n_chains):
        try:
            chain_id = coot.chain_id_py(imol, i)
        except Exception:  # noqa: BLE001
            continue
        if isinstance(chain_id, str):
            points.extend(chain_points(imol, chain_id))
    return points


# ---------------------------------------------------------------------------
#   Moving the view
# ---------------------------------------------------------------------------


def _apply(points: Sequence[Point], centre: bool = True,
           padding: float = PADDING) -> Optional[float]:
    """Zoom to fit *points* (and centre on them), returning the zoom set."""
    if coot is None or not points:
        return None
    zoom = zoom_for_extent(points_extent(points), padding)
    if centre:
        middle = centre_of(points)
        if middle is not None:
            coot.set_rotation_centre(*middle)
    coot.set_zoom(zoom)
    coot.graphics_draw()
    return zoom


def show_residue(imol: int, chain_id: str, resno: int,
                 ins_code: str = "") -> Optional[float]:
    """Centre on a residue and zoom to fit it, returning the zoom set.

    The go-to-atom call sets Coot's go-to state (and checks the residue is
    really there), but the view is then centred on the residue's *bounding-box
    middle* rather than left on the intelligent atom.  That matters: the CA of
    a long side chain sits at one end of the residue, so centring there while
    sizing the zoom from the whole bounding box frames it off-centre, pushing
    the far tip - and anything it is clashing with - to the edge or off screen
    entirely.  Centring on the middle makes the usable radius symmetric, which
    is what lets a contact partner stay in view (see :func:`show_contact`).

    The residue is still the active one afterwards: Coot picks the active
    residue as the displayed atom nearest the rotation centre, and that is now
    this residue's own middle.
    """
    if coot is None:
        return None
    points = residue_points(imol, chain_id, resno, ins_code)
    try:
        coot.set_go_to_atom_molecule(imol)
        if coot.set_go_to_atom_from_res_spec_py([chain_id, resno, ins_code]) <= 0:
            return None
    except Exception:  # noqa: BLE001 - a view failure must not fail a command
        return None
    return _apply(points)


def atom_point(imol: int, chain_id: str, resno: int, ins_code: str,
               atom_name: str) -> Optional[Point]:
    """The position of one named atom, or ``None`` if it is not there."""
    wanted = atom_name.strip()
    for atom in scoring.residue_atoms(imol, chain_id, resno, ins_code):
        try:
            if str(atom[0][0]).strip() == wanted:
                return (float(atom[2][0]), float(atom[2][1]),
                        float(atom[2][2]))
        except (IndexError, TypeError, ValueError):
            continue
    return None


def show_points(points: Sequence[Point],
                padding: float = PADDING) -> Optional[float]:
    """Centre on arbitrary positions and zoom to fit them all."""
    return _apply(points, padding=padding)


def show_contact(imol: int, spec_1, spec_2) -> Optional[float]:
    """Frame both partners of a contact - a clash, an H-bond - together.

    Takes the two atom specs Coot reports for the pair (see
    :func:`coot_commands.scoring.parse_atom_spec`) and frames *both residues*,
    not just the two atoms.  Two clashing atoms are by definition ~3 A apart,
    so fitting the atoms alone would zoom in past the point of usefulness:
    what you need to see is the two side chains and which way they would have
    to move.  Falls back to the atoms themselves if the residues cannot be
    read, and returns ``None`` if neither can.
    """
    points: List[Point] = []
    for spec in (spec_1, spec_2):
        parsed = scoring.parse_atom_spec(spec)
        if parsed is None:
            continue
        chain_id, resno, ins_code, atom_name = parsed
        residue = residue_points(imol, chain_id, resno, ins_code)
        if residue:
            points.extend(residue)
        else:
            point = atom_point(imol, chain_id, resno, ins_code, atom_name)
            if point is not None:
                points.append(point)
    return _apply(points) if points else None


def show_span(imol: int, chain_id: str, res1: int, res2: int) -> Optional[float]:
    """Centre on a residue range and zoom to fit all of it."""
    return _apply(span_points(imol, chain_id, res1, res2),
                  padding=PADDING_SPAN)


def show_chain(imol: int, chain_id: str) -> Optional[float]:
    """Centre on a chain and zoom to fit it."""
    return _apply(chain_points(imol, chain_id), padding=PADDING_WHOLE)


def show_molecule(imol: int) -> Optional[float]:
    """Centre on a whole molecule and zoom to fit it."""
    return _apply(molecule_points(imol), padding=PADDING_WHOLE)


# ---------------------------------------------------------------------------
#   The automatic hook
# ---------------------------------------------------------------------------


def follow_residue(imol: int, chain_id: str, resno: int,
                   ins_code: str = "") -> None:
    """Bring a residue into view, if autozoom is on.

    Called from :func:`coot_commands.types.resolve_residue`, so it runs once
    per per-residue command, *before* the command acts - which is what you
    want when the assistant is working: you see the residue, then you see it
    change.  Never raises: failing to move the view must not fail the command
    that was going to do the real work.
    """
    if not autozoom_enabled() or coot is None:
        return
    try:
        show_residue(imol, chain_id, resno, ins_code)
    except Exception:  # noqa: BLE001 - the view is never worth an exception
        pass
