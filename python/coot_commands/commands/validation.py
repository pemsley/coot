# coot_commands/commands/validation.py
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

"""Validation commands that summarise model outliers as text.

Rather than opening Coot's (C++/GTK) validation graphs - which are not
scriptable from the embedded interpreter - these commands call the
low-level ``coot.*`` scoring functions directly (the same ones the
quick-test validation uses) and report a short text summary.
"""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (resolve_model, resolve_map, as_float, as_int,
                                 first_present, first_difference_map, ArgType,
                                 CommandError, ACTIVE_MODEL_NOTE, ACTIVE_MAP_NOTE)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Validation"

# Coordinates of the blobs found by the most recent "find blobs" command,
# largest first, so a follow-up "go to blob N" can recentre on one.  Each
# entry is an (x, y, z) tuple; index 0 is blob 1.
_last_blobs: list[tuple[float, float, float]] = []


def _rama_improbables(imol: int) -> list:
    """Residues with Ramachandran probability below 0.02."""
    try:
        rs = coot.all_molecule_ramachandran_score_py(imol)
        scored = rs[5]
        return [item for item in scored if item[2] < 0.02]
    except Exception:
        return []


def _atom_overlaps(imol: int) -> list:
    """Atom overlaps with a clash volume above 2.0 A^3."""
    try:
        overlaps = coot.molecule_atom_overlaps_py(imol)
        return [o for o in overlaps
                if isinstance(o, dict) and o.get("overlap-volume", 0) > 2.0]
    except Exception:
        return []


def _as_list(value) -> list:
    """Coerce a maybe-list binding result to a list."""
    return list(value) if isinstance(value, (list, tuple)) else []


@command(r"validate anomalies(?: (?:of |in |for )?model (?P<model>\S+))?|"
         r"(?:validate|check|find) (?:outliers|anomalies)"
         r"(?: (?:of |in |for )?model (?P<model2>\S+))?",
         examples=["validate anomalies", "find outliers"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "model2": ArgType.MODEL},
         notes="Summarises geometry outliers - Ramachandran improbables, "
               "atom overlaps (clashes), C-beta deviations and chiral volume "
               "errors - as text. " + ACTIVE_MODEL_NOTE)
def validate_anomalies(model: Optional[str] = None,
                       model2: Optional[str] = None) -> str:
    """Summarise model outliers (Ramachandran, clashes, C-beta, chirals)."""
    imol = resolve_model(first_present(model, model2))
    if coot is None:
        return f"Validated model {imol}"

    rama = _rama_improbables(imol)
    overlaps = _atom_overlaps(imol)
    c_beta = _as_list(coot.c_beta_deviations_py(imol)) \
        if hasattr(coot, "c_beta_deviations_py") else []
    chirals = _as_list(coot.chiral_volume_errors_py(imol)) \
        if hasattr(coot, "chiral_volume_errors_py") else []

    total = len(rama) + len(overlaps) + len(c_beta) + len(chirals)
    if total == 0:
        return f"No outliers found in model {imol}"

    parts = [f"Model {imol}: {total} outlier(s) -",
             f"{len(rama)} Ramachandran, {len(overlaps)} clash(es), "
             f"{len(c_beta)} C-beta, {len(chirals)} chiral."]
    if rama:
        worst = sorted(rama, key=lambda item: item[2])[:3]
        worst_str = ", ".join(f"{item[1]} ({item[2]:.3f})" for item in worst)
        parts.append(f"Worst Ramachandran: {worst_str}")
    return " ".join(parts)




@command(r"(?:validate|check|find|list) (?:unmodell?ed |unbuilt |missing )?blobs"
         r"(?: (?:of |in |for )?model (?P<model>\S+))?"
         r"(?: (?:against |using |with |in )?map (?P<map>\S+))?"
         r"(?: (?:above |over |at )(?P<sigma>[\d.]+)\s*(?:sigma|rmsd)?)?",
         examples=["find blobs", "check unmodelled blobs",
                   "find blobs above 2 sigma"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "map": ArgType.MAP},
         notes="Lists peaks of density not accounted for by the model - "
               "candidate sites for waters, ligands or unbuilt residues. "
               "The search masks out density within 1.9 A of the model and "
               "finds the peaks left over, so it belongs on a difference "
               "(mFo-DFc) map: by default it uses the loaded difference map "
               "at a 3.0 sigma cut-off. On a 2mFo-DFc map almost everything "
               "is above a low cut-off, so blobs would appear all over - name "
               "a map with 'using map N' only if you mean to. Add e.g. "
               "'above 2 sigma' to change the threshold.")
def validate_unmodelled_blobs(model: Optional[str] = None,
                              map: Optional[str] = None,
                              sigma: Optional[str] = None) -> str:
    """Summarise unmodelled density blobs (candidate build sites)."""
    imol = resolve_model(model)
    cut_off = as_float(sigma, "sigma cut-off") if sigma is not None else 3.0
    if coot is None:
        return f"Validated model {imol}"

    # Default to the difference map - the peaks that survive masking the model
    # there are genuinely unmodelled density.  Fall back to the active map only
    # if no difference map is loaded (and warn, since that is noisy).
    warning = ""
    if map is not None:
        imol_map = resolve_map(map)
    else:
        imol_map = first_difference_map()
        if imol_map is None:
            imol_map = resolve_map(None)
            warning = (" (no difference map loaded, so this used the active "
                       "map, where blobs are unreliable)")

    blobs = _as_list(coot.find_blobs_py(imol, imol_map, cut_off))
    global _last_blobs
    if not blobs:
        _last_blobs = []
        return (f"No unmodelled blobs found in model {imol} on map {imol_map} "
                f"above {cut_off:g} sigma{warning}")

    # find_blobs_py returns [[x, y, z], volume] pairs; report the largest and
    # remember their positions so "go to blob N" can recentre on one.
    blobs.sort(key=lambda b: b[1], reverse=True)
    _last_blobs = [(b[0][0], b[0][1], b[0][2]) for b in blobs]
    worst = "; ".join(
        f"blob {i} at ({b[0][0]:.1f}, {b[0][1]:.1f}, {b[0][2]:.1f}) "
        f"vol {b[1]:.1f}"
        for i, b in enumerate(blobs[:3], start=1))
    return (f"Model {imol}: {len(blobs)} unmodelled blob(s) on map {imol_map} "
            f"above {cut_off:g} sigma{warning}. Largest: {worst}")


@command(r"(?:go to|centre on|center on|goto) blob (?P<n>\d+)",
         examples=["go to blob 1", "centre on blob 2"],
         category=CATEGORY,
         notes="Centres the view on one of the blobs from the most recent "
               "'find blobs' command, numbered from 1 (largest first). Run "
               "'find blobs' first to populate the list.")
def go_to_blob(n: str) -> str:
    """Centre the view on a blob from the last 'find blobs' result."""
    index = as_int(n, "blob number")
    if not _last_blobs:
        raise CommandError("no blobs to go to - run 'find blobs' first")
    if index < 1 or index > len(_last_blobs):
        raise CommandError(
            f"blob {index} does not exist - there are {len(_last_blobs)} "
            f"(numbered 1 to {len(_last_blobs)})")
    x, y, z = _last_blobs[index - 1]
    if coot is not None:
        coot.set_rotation_centre(x, y, z)
    return f"Centred on blob {index} at ({x:.1f}, {y:.1f}, {z:.1f})"


# ---------------------------------------------------------------------------
#   Per-category text validation
# ---------------------------------------------------------------------------
#
# These report a single validation category each (the granular counterpart to
# "validate anomalies", which bundles several).  They call the scriptable
# ``coot.*_py`` scorers and summarise the result as text, so they work in the
# terminal without opening a GTK graph.  A trailing "of model N" targets a
# particular model; otherwise they act on the active one.

# The optional "of/in/for model N" suffix, shared by every command here.  Two
# spellings of the capture group let a pattern offer it on both sides of an
# alternation without clashing on the group name.
_MODEL = r"(?: (?:of |in |for )?model (?P<model>\S+))?"
_MODEL2 = r"(?: (?:of |in |for )?model (?P<model2>\S+))?"


@command(r"(?:(?:check|validate|list|show) )?ramachandran(?: outliers)?" + _MODEL,
         examples=["check ramachandran", "ramachandran outliers"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Lists residues in improbable regions of the Ramachandran plot "
               "(probability below 0.02). " + ACTIVE_MODEL_NOTE)
def check_ramachandran(model: Optional[str] = None) -> str:
    """List Ramachandran outliers for a model."""
    imol = resolve_model(model)
    if coot is None:
        return f"Checked Ramachandran for model {imol}"
    rama = _rama_improbables(imol)
    if not rama:
        return f"No Ramachandran outliers in model {imol}"
    worst = sorted(rama, key=lambda item: item[2])[:5]
    detail = ", ".join(f"{item[1]} ({item[2]:.3f})" for item in worst)
    return f"Model {imol}: {len(rama)} Ramachandran outlier(s). Worst: {detail}"


@command(r"(?:(?:check|validate|list|show) )?rotamers?(?: outliers)?" + _MODEL,
         examples=["check rotamers", "rotamer outliers"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Reports the least probable side-chain rotamers (low "
               "probability = unusual). " + ACTIVE_MODEL_NOTE)
def check_rotamers(model: Optional[str] = None) -> str:
    """List the least probable rotamers for a model."""
    imol = resolve_model(model)
    if coot is None:
        return f"Checked rotamers for model {imol}"
    try:
        rots = _as_list(coot.rotamer_graphs_py(imol))
    except Exception:
        rots = []
    # Each entry is [chain, resno, ins_code, probability, name]; sort worst first.
    scored = [r for r in rots if isinstance(r, (list, tuple)) and len(r) >= 5]
    if not scored:
        return f"No rotamers scored for model {imol}"
    scored.sort(key=lambda r: r[3])
    worst = scored[:5]
    detail = ", ".join(f"{r[0]}/{r[1]} {r[4]} ({r[3]:g})" for r in worst)
    return (f"Model {imol}: {len(scored)} rotamer(s) scored. "
            f"Least probable: {detail}")


@command(r"(?:(?:check|validate|list|show|find) )?"
         r"(?:clashes|atom overlaps|overlaps)" + _MODEL,
         examples=["check clashes", "atom overlaps"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Lists atom overlaps (steric clashes) with a clash volume "
               "above 2 A^3. " + ACTIVE_MODEL_NOTE)
def check_clashes(model: Optional[str] = None) -> str:
    """List steric clashes (atom overlaps) for a model."""
    imol = resolve_model(model)
    if coot is None:
        return f"Checked clashes for model {imol}"
    overlaps = _atom_overlaps(imol)
    if not overlaps:
        return f"No significant clashes in model {imol}"
    worst = sorted(overlaps, key=lambda o: o.get("overlap-volume", 0),
                   reverse=True)[:3]
    detail = ", ".join(f"{o.get('overlap-volume', 0):.1f} A^3" for o in worst)
    return f"Model {imol}: {len(overlaps)} clash(es). Largest: {detail}"


@command(r"(?:(?:check|validate|list|show|find) )?cis[- ]?peptides?" + _MODEL,
         examples=["check cis peptides", "cis peptides"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Counts cis peptide bonds. Cis peptides before proline are "
               "common and usually fine; other cis peptides are worth "
               "checking. " + ACTIVE_MODEL_NOTE)
def check_cis_peptides(model: Optional[str] = None) -> str:
    """Count cis peptide bonds in a model."""
    imol = resolve_model(model)
    if coot is None:
        return f"Checked cis peptides for model {imol}"
    try:
        cis = _as_list(coot.cis_peptides_py(imol))
    except Exception:
        cis = []
    if not cis:
        return f"No cis peptides in model {imol}"
    return f"Model {imol}: {len(cis)} cis peptide(s)"


@command(r"(?:(?:check|validate|list|show) )?"
         r"(?:non[- ]?standard|unusual) residues?" + _MODEL,
         examples=["check non-standard residues", "non-standard residues"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Lists residue types that are not standard amino acids or "
               "water - ligands, modified residues and the like. "
               + ACTIVE_MODEL_NOTE)
def check_non_standard_residues(model: Optional[str] = None) -> str:
    """List non-standard residue types in a model."""
    imol = resolve_model(model)
    if coot is None:
        return f"Checked non-standard residues for model {imol}"
    try:
        names = [str(n) for n in _as_list(coot.non_standard_residue_names_py(imol))]
    except Exception:
        names = []
    if not names:
        return f"No non-standard residues in model {imol}"
    unique = sorted(set(names))
    return (f"Model {imol}: {len(unique)} non-standard residue type(s): "
            f"{', '.join(unique)}")


@command(r"(?:(?:check|validate|list|show|find) )?missing atoms?" + _MODEL,
         examples=["check missing atoms", "missing atoms"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Lists residues that are missing modelled atoms. "
               + ACTIVE_MODEL_NOTE)
def check_missing_atoms(model: Optional[str] = None) -> str:
    """List residues with missing atoms in a model."""
    imol = resolve_model(model)
    if coot is None:
        return f"Checked missing atoms for model {imol}"
    try:
        info = _as_list(coot.missing_atom_info_py(imol))
    except Exception:
        info = []
    if not info:
        return f"No residues with missing atoms in model {imol}"
    specs = [f"{e[0]}/{e[1]}" for e in info
             if isinstance(e, (list, tuple)) and len(e) >= 2]
    shown = ", ".join(specs[:5])
    more = "" if len(specs) <= 5 else f" (+{len(specs) - 5} more)"
    return (f"Model {imol}: {len(info)} residue(s) with missing atoms: "
            f"{shown}{more}")


@command(r"(?:(?:check|validate|list|show) )?highly[- ]?coordinated waters?"
         + _MODEL + r"|(?:check|validate|list|show) waters?" + _MODEL2,
         examples=["check waters", "highly coordinated waters"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "model2": ArgType.MODEL},
         notes="Lists waters with an unusually high number of close contacts "
               "(coordination number 5 or more within 3.2 A), which may be "
               "misassigned ions. " + ACTIVE_MODEL_NOTE)
def check_waters(model: Optional[str] = None,
                 model2: Optional[str] = None) -> str:
    """List highly-coordinated waters (possible ions) in a model."""
    imol = resolve_model(first_present(model, model2))
    if coot is None:
        return f"Checked waters for model {imol}"
    try:
        waters = _as_list(coot.highly_coordinated_waters_py(imol, 5, 3.2))
    except Exception:
        waters = []
    if not waters:
        return f"No highly-coordinated waters in model {imol}"
    return (f"Model {imol}: {len(waters)} highly-coordinated water(s) "
            f"(possible ions)")


# ---------------------------------------------------------------------------
#   GUI validation tools
# ---------------------------------------------------------------------------
#
# These open Coot's interactive validation windows rather than returning text.
# The overlay (dynamic_validation_internal) is the modern all-in-one panel;
# the others open a single classic tool.

@command(r"(?:open|show|display) validation(?: overlay| panel| graphs?)?"
         + _MODEL + r"|validation(?: overlay| panel| graphs?)" + _MODEL2,
         examples=["open validation", "validation overlay"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "model2": ArgType.MODEL},
         notes="Opens Coot's interactive validation overlay (Ramachandran, "
               "rotamers, density fit and more) for a model against the "
               "active map. " + ACTIVE_MODEL_NOTE + " " + ACTIVE_MAP_NOTE)
def open_validation(model: Optional[str] = None,
                    model2: Optional[str] = None) -> str:
    """Open the interactive validation overlay for a model."""
    imol = resolve_model(first_present(model, model2))
    if coot is None:
        return f"Opened validation overlay for model {imol}"
    imol_map = resolve_map(None)
    coot.dynamic_validation_internal(imol, imol_map)
    return f"Opened validation overlay for model {imol} against map {imol_map}"


@command(r"(?:(?:show|find|generate|check) )?difference map peaks"
         r"(?: (?:above|over|at) (?P<sigma>[\d.]+)\s*(?:sigma|rmsd)?)?" + _MODEL,
         examples=["difference map peaks", "difference map peaks above 4 sigma"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Marks peaks in the Fo-Fc difference map - candidate sites for "
               "missing atoms, waters or ligands, and for parts of the model "
               "in wrong density. Uses a 4 sigma cut-off by default. "
               "Needs a difference map to be loaded. " + ACTIVE_MODEL_NOTE)
def show_difference_map_peaks(model: Optional[str] = None,
                             sigma: Optional[str] = None) -> str:
    """Mark peaks in the difference map (missing/wrong density)."""
    imol = resolve_model(model)
    level = as_float(sigma, "sigma cut-off") if sigma is not None else 4.0
    imol_diff = first_difference_map()
    if imol_diff is None:
        raise CommandError(
            "no difference map loaded - open your data (e.g. 'load tutorial') "
            "to get an Fo-Fc difference map first")
    if coot is not None:
        # difference_map_peaks(diff_map, model, level, max_closeness,
        #                      positive, negative, around_model_only)
        coot.difference_map_peaks(imol_diff, imol, level, 2.0, 1, 1, 1)
    return (f"Marked difference map peaks for model {imol} on map "
            f"{imol_diff} above {level:g} sigma")


@command(r"(?:(?:check|validate|show|find) )?"
         r"(?:gln[/ ]?(?:and )?asn|glutamine and asparagine)"
         r"(?: (?:b[- ]?factor )?outliers)?" + _MODEL,
         examples=["check gln and asn", "gln asn outliers"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Opens the Gln/Asn B-factor outlier analysis, flagging "
               "glutamine and asparagine side chains that may need a 180 "
               "degree flip. " + ACTIVE_MODEL_NOTE)
def check_gln_asn(model: Optional[str] = None) -> str:
    """Open the Gln/Asn side-chain flip (B-factor outlier) analysis."""
    imol = resolve_model(model)
    if coot is not None:
        coot.gln_asn_b_factor_outliers(imol)
    return f"Opened Gln/Asn outlier analysis for model {imol}"


