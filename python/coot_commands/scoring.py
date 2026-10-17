# coot_commands/scoring.py
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

"""Measure one residue, so "did that help?" has an answer.

Everything else in this package *does* something; this module *measures*.
It exists because the assistant cannot follow a procedure like "auto-fit
the rotamer, and if that does not improve things, try a backrub" without a
number to compare before and after.  Coot's validation commands score a
whole model and report the worst few residues, which is the wrong
granularity for that decision.

Three quantities describe a residue well enough to choose between fits:

* **rotamer probability** - how usual the side-chain conformation is,
  from Coot's Richardson rotamer tables via ``coot.rotamer_score``.  It is
  a *percentage* (0-100), not a fraction, and what counts as bad depends on
  the residue type - see :func:`rotamer_band_for`.  Glycine and alanine have
  no rotamer, so it is ``None``.
* **density fit** - the sum of the map value at each atom, weighted by
  occupancy (``coot.density_score_residue_py``).  Bigger is better, but it
  scales with the number of atoms, so it is comparable *for the same
  residue before and after a change* and meaningless between residues.
  :attr:`ResidueScore.density_per_atom` is the comparable-ish version.
* **clash volume** - the summed overlap volume of this residue's atoms
  with everything else.  Coot only exposes whole-molecule overlaps, so
  this costs a full overlap calculation; callers that measure repeatedly
  in a loop should pass ``clashes=False`` for the intermediate steps.

:func:`capture_positions` and :func:`restore_positions` bracket a trial
edit so a fit that turns out worse can be undone precisely, without
disturbing Coot's own undo stack.

The module imports and runs without Coot (every reader degrades to
``None``), so the comparison logic is testable standalone.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import coot
except ImportError:
    coot = None


# Rotamer probability bands, as percentages, following MolProbity's
# convention - and Coot's own: its rotamer graphs saturate the distortion bar
# at a probability of 0.3 (graphics-info-graphs.cc), which is the outlier line.
# What counts as a bad rotamer depends on the residue type, so there is no one
# threshold.  The score is a probability density at the actual chi angles, and
# the more chi angles a residue has the more rotamers it is spread across and
# the lower every individual score becomes: arginine's *best possible* score is
# about 9%, where valine's is about 73%.  A flat cut-off therefore calls a
# badly-placed valine "favoured" and flags every arginine in the structure.
#
# The bands below come from the guidance in mcp/docs/skills/validation, which
# gives the "good" figure per residue class; the outlier line is set at half of
# it, which reproduces the worked examples given there (valine at 5% terrible,
# leucine at 5% poor, arginine at 5% good but 0.5% poor).
#
# Grouping is by how many rotamers a residue type has, which is what drives the
# score - not strictly by chi count, which is why leucine and isoleucine sit
# with the three-chi group.
_ROTAMER_BANDS = (
    (("VAL", "THR", "SER", "CYS", "PRO"),               20.0, 40.0),
    (("PHE", "TYR", "ASP", "ASN", "HIS", "TRP"),        10.0, 20.0),
    (("GLU", "GLN", "MET", "ILE", "LEU"),                5.0, 10.0),
    (("LYS", "ARG"),                                     2.5,  5.0),
)

# For a residue type not listed above (a modified residue, say): the middle of
# the range, so it is neither flagged nor excused on no evidence.
DEFAULT_ROTAMER_BAND = (5.0, 10.0)


def rotamer_band_for(res_name: str) -> Tuple[float, float]:
    """``(outlier below, favoured at or above)`` for a residue type, in percent."""
    name = (res_name or "").strip().upper()
    for names, outlier, favoured in _ROTAMER_BANDS:
        if name in names:
            return outlier, favoured
    return DEFAULT_ROTAMER_BAND

# Residues with no side-chain torsions, so no rotamer to score or fit.
NO_ROTAMER = frozenset(("GLY", "ALA"))

# Fractional change in density fit below which two fits count as equally good.
# Real-space fits jitter by well under a percent, so a candidate has to beat
# the incumbent by more than this to be preferred - otherwise the assistant
# would keep "improving" a residue by chasing map noise.
DENSITY_MARGIN = 0.02

# Water and the like: scoring a rotamer for these is meaningless.
_NON_PROTEIN = frozenset(("HOH", "WAT", "DOD"))


class ResidueScore:
    """A snapshot of how good one residue looks, at one moment.

    Any field may be ``None`` when the underlying measurement was not
    available (no map loaded, no rotamer for this residue type, Coot not
    present); :meth:`describe` says so rather than inventing a number.
    """

    def __init__(self, imol: int, chain_id: str, resno: int, ins_code: str,
                 res_name: str, rotamer: Optional[float],
                 density: Optional[float], n_atoms: int,
                 clash_volume: Optional[float] = None,
                 n_clashes: Optional[int] = None) -> None:
        self.imol = imol
        self.chain_id = chain_id
        self.resno = resno
        self.ins_code = ins_code
        self.res_name = res_name
        self.rotamer = rotamer
        self.density = density
        self.n_atoms = n_atoms
        self.clash_volume = clash_volume
        self.n_clashes = n_clashes

    @property
    def spec(self) -> str:
        """The residue as ``A/45``, with any insertion code appended."""
        return f"{self.chain_id}/{self.resno}{self.ins_code.strip()}"

    @property
    def density_per_atom(self) -> Optional[float]:
        """Density fit divided by atom count - the size-independent version."""
        if self.density is None or self.n_atoms <= 0:
            return None
        return self.density / self.n_atoms

    @property
    def thresholds(self) -> Tuple[float, float]:
        """This residue type's ``(outlier below, favoured at or above)``."""
        return rotamer_band_for(self.res_name)

    @property
    def is_outlier(self) -> bool:
        """True when the side chain is in a rotamer outlier conformation.

        Judged against this residue type's own band, not a flat cut-off: an
        arginine at 3% is unremarkable where a valine at 3% is badly wrong.
        A residue with no rotamer (glycine, alanine) is never an outlier; nor
        is one we could not score, since we should not claim a problem we did
        not measure.
        """
        return self.rotamer is not None and self.rotamer < self.thresholds[0]

    @property
    def rotamer_band(self) -> str:
        """``"outlier"`` / ``"allowed"`` / ``"favoured"`` for the rotamer."""
        if self.rotamer is None:
            return "none"
        outlier, favoured = self.thresholds
        if self.rotamer < outlier:
            return "outlier"
        if self.rotamer < favoured:
            return "allowed"
        return "favoured"

    def describe(self) -> str:
        """A one-line summary, phrased for the assistant to read and compare.

        Every number carries its units and its threshold, because the model
        reads this string with no other context and would otherwise have to
        guess whether "0.21" is good.
        """
        parts = [f"{self.spec} {self.res_name} of model {self.imol}:"]
        parts.append(self._rotamer_phrase())
        parts.append(self._density_phrase())
        if self.clash_volume is not None:
            if self.clash_volume <= 0.0:
                parts.append("no clashes")
            else:
                parts.append(f"clashes {self.clash_volume:.1f} A^3 over "
                             f"{self.n_clashes} atom pair(s)")
        return " ".join(parts)

    def _rotamer_phrase(self) -> str:
        if self.res_name in NO_ROTAMER:
            return f"rotamer n/a ({self.res_name} has no side-chain rotamer);"
        if self.rotamer is None:
            return "rotamer not scored;"
        if self.rotamer <= 0.0:
            return ("rotamer 0% - matches no library rotamer (an outlier, or "
                    "the side chain has missing atoms);")
        outlier, favoured = self.thresholds
        return (f"rotamer {self.rotamer:.2f}% ({self.rotamer_band}; for "
                f"{self.res_name} outlier is below {outlier:g}% and favoured "
                f"is {favoured:g}% or more - these depend on the residue "
                f"type, so do not compare them across types);")

    def _density_phrase(self) -> str:
        if self.density is None:
            return "density fit unavailable (no refinement map set);"
        per_atom = self.density_per_atom
        return (f"density fit {self.density:.2f} over {self.n_atoms} atoms "
                f"({per_atom:.2f}/atom, higher is better);")


# ---------------------------------------------------------------------------
#   Readers
# ---------------------------------------------------------------------------


def residue_name(imol: int, chain_id: str, resno: int, ins_code: str = "") -> str:
    """The residue's three-letter type, or ``""`` if it cannot be read.

    ``residue_name`` returns a ``std::string`` straight through SWIG;
    ``residue_name_py`` is the older binding, kept as a fallback so this works
    across builds.  Neither raises for a residue that is not there.
    """
    if coot is None:
        return ""
    for getter in ("residue_name", "residue_name_py"):
        fn = getattr(coot, getter, None)
        if fn is None:
            continue
        try:
            name = fn(imol, chain_id, resno, ins_code)
        except Exception:  # noqa: BLE001 - a missing residue must not raise here
            continue
        if isinstance(name, str) and name.strip():
            return name.strip()
    return ""


def residue_atoms(imol: int, chain_id: str, resno: int,
                  ins_code: str = "") -> List:
    """The residue's atoms as ``residue_info_py`` returns them, or ``[]``.

    Each entry is ``[[name, alt_conf], [occ, b, element, seg_id], [x, y, z],
    index]``.
    """
    if coot is None:
        return []
    try:
        info = coot.residue_info_py(imol, chain_id, resno, ins_code)
    except Exception:  # noqa: BLE001
        return []
    return list(info) if isinstance(info, (list, tuple)) else []


def rotamer_percent(imol: int, chain_id: str, resno: int, ins_code: str = "",
                    alt_conf: str = "") -> Optional[float]:
    """Rotamer probability as a percentage, or ``None`` if it has no rotamer.

    ``coot.rotamer_score`` returns 0.0 for every failure mode it knows -
    glycine/alanine, missing atoms, no matching library rotamer - so we screen
    out the residue types that legitimately have no rotamer first, and let a
    remaining 0.0 through as the meaningful signal it is (see
    :meth:`ResidueScore._rotamer_phrase`).
    """
    if coot is None:
        return None
    name = residue_name(imol, chain_id, resno, ins_code)
    if name in NO_ROTAMER or name in _NON_PROTEIN:
        return None
    try:
        return float(coot.rotamer_score(imol, chain_id, resno, ins_code,
                                        alt_conf))
    except Exception:  # noqa: BLE001
        return None


def density_fit(imol: int, chain_id: str, resno: int, ins_code: str = "",
                imol_map: Optional[int] = None) -> Optional[float]:
    """Summed map value over the residue's atoms, or ``None`` with no map."""
    if coot is None:
        return None
    if imol_map is None:
        try:
            imol_map = coot.imol_refinement_map()
        except Exception:  # noqa: BLE001
            return None
    if imol_map is None or imol_map < 0:
        return None
    try:
        return float(coot.density_score_residue_py(
            imol, [chain_id, resno, ins_code], imol_map))
    except Exception:  # noqa: BLE001
        return None


# ---------------------------------------------------------------------------
#   Specs
# ---------------------------------------------------------------------------
#
# Coot's validation bindings answer with specs rather than text, in two shapes:
#
#   atom spec     [user_data, chain_id, res_no, ins_code, atom_name, alt_conf]
#                 (graphics_info_t::atom_spec_to_py)
#   residue spec  [chain_id, res_no, ins_code]
#                 (residue_spec_to_py, c-interface.cc)
#
# A command that reports "8 clashes" without saying where they are leaves the
# reader - and the assistant - with nothing to act on, so the parsers below
# turn those specs into "A/45" and "A/45 CB". The residue spec grew its current
# shape recently (it used to carry a leading status flag), so parsing tolerates
# an extra leading element rather than trusting the length.


def parse_residue_spec(spec) -> Optional[Tuple[str, int, str]]:
    """``(chain_id, resno, ins_code)`` from a residue spec, or ``None``.

    Accepts the three-element form and the older four-element one that
    prefixed a status flag, telling them apart by which element is the chain
    id (a string) followed by the residue number (an int).
    """
    if not isinstance(spec, (list, tuple)):
        return None
    items = list(spec)
    # Drop a leading flag/index so both spec generations parse the same way.
    if items and not isinstance(items[0], str):
        items = items[1:]
    if len(items) < 2 or not isinstance(items[0], str):
        return None
    try:
        resno = int(items[1])
    except (TypeError, ValueError):
        return None
    ins_code = str(items[2]) if len(items) > 2 and isinstance(items[2], str) else ""
    return items[0].strip(), resno, ins_code.strip()


def parse_atom_spec(spec) -> Optional[Tuple[str, int, str, str]]:
    """``(chain_id, resno, ins_code, atom_name)`` from an atom spec, or ``None``."""
    if not isinstance(spec, (list, tuple)) or len(spec) < 5:
        return None
    try:
        resno = int(spec[2])
    except (TypeError, ValueError):
        return None
    return (str(spec[1]).strip(), resno, str(spec[3]).strip(),
            str(spec[4]).strip())


def format_residue_spec(spec) -> str:
    """A residue spec as ``A/45``, or ``"?"`` if it cannot be read."""
    parsed = parse_residue_spec(spec)
    if parsed is None:
        return "?"
    chain_id, resno, ins_code = parsed
    return f"{chain_id}/{resno}{ins_code}"


def format_atom_spec(spec) -> str:
    """An atom spec as ``A/45 CB``, or ``"?"`` if it cannot be read."""
    parsed = parse_atom_spec(spec)
    if parsed is None:
        return "?"
    chain_id, resno, ins_code, atom_name = parsed
    residue = f"{chain_id}/{resno}{ins_code}"
    return f"{residue} {atom_name}" if atom_name else residue


def _spec_matches(spec, chain_id: str, resno: int, ins_code: str = "") -> bool:
    """Does an atom spec name this residue?"""
    parsed = parse_atom_spec(spec)
    if parsed is None:
        return False
    return (parsed[0] == chain_id and parsed[1] == resno
            and parsed[2] == ins_code.strip())


def atom_overlaps(imol: int, n_max: int = -1) -> List[Dict]:
    """Every atom overlap in a molecule, largest first, or ``[]``.

    ``molecule_atom_overlaps_py`` takes a maximum pair count as its second
    argument; ``-1`` asks for all of them, which is what a per-residue
    question needs (the list is sorted globally, so a bounded request could
    truncate away the overlaps belonging to the residue we care about).
    Older builds took only the molecule number, hence the fallback.
    """
    if coot is None:
        return []
    try:
        overlaps = coot.molecule_atom_overlaps_py(imol, n_max)
    except TypeError:
        try:
            overlaps = coot.molecule_atom_overlaps_py(imol)
        except Exception:  # noqa: BLE001
            return []
    except Exception:  # noqa: BLE001
        return []
    if not isinstance(overlaps, (list, tuple)):
        return []  # the binding returns False on failure
    return [o for o in overlaps if isinstance(o, dict)]


def residue_clashes(imol: int, chain_id: str, resno: int,
                    ins_code: str = "") -> Tuple[Optional[float], Optional[int]]:
    """``(total overlap volume, number of overlapping pairs)`` for a residue.

    Returns ``(None, None)`` when overlaps cannot be computed at all, which is
    different from ``(0.0, 0)`` - "measured, and clash-free".
    """
    if coot is None:
        return None, None
    overlaps = atom_overlaps(imol)
    total, count = 0.0, 0
    for o in overlaps:
        here = (_spec_matches(o.get("atom-1-spec"), chain_id, resno, ins_code)
                or _spec_matches(o.get("atom-2-spec"), chain_id, resno,
                                 ins_code))
        if here:
            total += float(o.get("overlap-volume", 0.0))
            count += 1
    return total, count


def score_residue(imol: int, chain_id: str, resno: int, ins_code: str = "",
                  imol_map: Optional[int] = None,
                  clashes: bool = True) -> ResidueScore:
    """Measure one residue: rotamer, density fit and (optionally) clashes.

    *clashes* costs a whole-molecule overlap calculation, so pass ``False``
    for the intermediate measurements inside a fitting loop and ``True`` for
    the ones the user sees.
    """
    name = residue_name(imol, chain_id, resno, ins_code)
    atoms = residue_atoms(imol, chain_id, resno, ins_code)
    clash_volume, n_clashes = (residue_clashes(imol, chain_id, resno, ins_code)
                               if clashes else (None, None))
    return ResidueScore(
        imol=imol, chain_id=chain_id, resno=resno, ins_code=ins_code,
        res_name=name or "?",
        rotamer=rotamer_percent(imol, chain_id, resno, ins_code),
        density=density_fit(imol, chain_id, resno, ins_code, imol_map),
        n_atoms=len(atoms),
        clash_volume=clash_volume, n_clashes=n_clashes)


# ---------------------------------------------------------------------------
#   Comparison
# ---------------------------------------------------------------------------


def is_better(candidate: ResidueScore, incumbent: ResidueScore) -> bool:
    """Is *candidate* a better fit for this residue than *incumbent*?

    The order of preference matches how one would judge it by eye:

    1. Getting out of rotamer-outlier territory beats everything - an
       unusual side-chain conformation is the problem we set out to fix.
    2. Otherwise the better density fit wins, but only by more than
       :data:`DENSITY_MARGIN`, so noise-level differences do not count.
    3. With the density a wash, prefer the more probable rotamer.
    """
    if candidate.is_outlier != incumbent.is_outlier:
        return not candidate.is_outlier
    if candidate.density is not None and incumbent.density is not None:
        margin = DENSITY_MARGIN * abs(incumbent.density)
        if candidate.density > incumbent.density + margin:
            return True
        if candidate.density < incumbent.density - margin:
            return False
    if candidate.rotamer is not None and incumbent.rotamer is not None:
        return candidate.rotamer > incumbent.rotamer
    return False


def compare(before: ResidueScore, after: ResidueScore) -> str:
    """A short "x -> y" phrase describing what a fit changed."""
    bits = []
    if before.rotamer is not None and after.rotamer is not None:
        bits.append(f"rotamer {before.rotamer:.2f}% -> {after.rotamer:.2f}%")
    if before.density is not None and after.density is not None:
        bits.append(f"density {before.density:.2f} -> {after.density:.2f}")
    return ", ".join(bits) if bits else "no measurable change"


# ---------------------------------------------------------------------------
#   Trial edits
# ---------------------------------------------------------------------------
#
# A composite command tries a fit, measures it, and keeps it only if it helped.
# Coot's undo would serve for the "put it back" step, but it is a single shared
# stack that the user also drives from the toolbar; rewinding it from a command
# would silently eat whatever they had queued there.  Writing the coordinates
# back directly is both narrower and reversible in exactly one step.


def capture_positions(imol: int, chain_id: str, resnos: Sequence[int],
                      ins_code: str = "") -> List:
    """Record the coordinates of every atom in the given residues.

    Pass the neighbours as well as the residue itself where a fit can move
    them: a backrub rotates the flanking peptides, so restoring residue *i*
    alone would leave *i-1* and *i+1* displaced.
    """
    saved = []
    for resno in resnos:
        for atom in residue_atoms(imol, chain_id, resno, ins_code):
            try:
                name, alt_conf = atom[0][0], atom[0][1]
                x, y, z = atom[2][0], atom[2][1], atom[2][2]
            except (IndexError, TypeError):
                continue
            saved.append((chain_id, resno, ins_code, name, alt_conf, x, y, z))
    return saved


def restore_positions(imol: int, saved: Iterable) -> int:
    """Put coordinates captured by :func:`capture_positions` back.

    Returns the number of atoms restored.  Uses the batched
    ``set_atom_attributes_py`` so the whole residue moves in one structure
    edit and one redraw.
    """
    settings = []
    for (chain_id, resno, ins_code, name, alt_conf, x, y, z) in saved:
        for attribute, value in (("x", x), ("y", y), ("z", z)):
            settings.append([imol, chain_id, resno, ins_code, name, alt_conf,
                             attribute, float(value)])
    if not settings or coot is None:
        return 0
    coot.set_atom_attributes_py(settings)
    return len(settings) // 3
