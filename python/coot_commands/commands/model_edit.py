# coot_commands/commands/model_edit.py
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

"""Commands that edit a single residue: delete, mutate, rotamer fitting.

All of these act on one residue named by a spec like ``A/72`` (or ``A 72``);
with no model number they use the active model.  The rotamer commands
(``autofit``) need a refinement map, mirroring :mod:`coot_commands.commands.refine`.
"""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (RES_SPEC, OPT_RES_SPEC, resolve_residue,
                                 ArgType, CommandError, centre_on_residue,
                                 ACTIVE_RESIDUE_NOTE)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Model editing"


@command(r"delete residue" + OPT_RES_SPEC + r"(?: (?:of |from )?model (?P<model>\S+))?",
         examples=["delete residue A/72", "delete residue A 72", "delete residue"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Deletes the named residue. " + ACTIVE_RESIDUE_NOTE)
def delete_residue(chain: Optional[str] = None, resno: Optional[str] = None,
                   model: Optional[str] = None) -> str:
    """Delete a single residue."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    if coot is not None:
        coot.delete_residue(imol, chain_id, res, ins)
    return f"Deleted {chain_id}/{res} of model {imol}"


@command(r"(?:replace|mutate) residue" + OPT_RES_SPEC +
         r"(?: (?:with|to|into) (?P<restype>[A-Za-z]{1,3}))"
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["replace residue A/72 with ALA", "mutate residue A 72 to GLY",
                   "mutate residue to ALA"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Mutates the residue to the given (1- or 3-letter) type. "
               + ACTIVE_RESIDUE_NOTE)
def replace_residue(restype: str, chain: Optional[str] = None,
                    resno: Optional[str] = None,
                    model: Optional[str] = None) -> str:
    """Mutate a residue to another type."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    target = restype.upper()
    if coot is not None:
        coot.mutate(imol, chain_id, res, ins, target)
    return f"Mutated {chain_id}/{res} of model {imol} to {target}"


@command(r"(?:add|build) alt[- ]?conf(?:ormer|ormation)?" + OPT_RES_SPEC +
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["add altconf A/72", "add alt conf A 72", "add altconf"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Adds an alternate conformer to the named residue. "
               + ACTIVE_RESIDUE_NOTE)
def add_altconf(chain: Optional[str] = None, resno: Optional[str] = None,
                model: Optional[str] = None) -> str:
    """Add an alternate conformation to a residue."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    if coot is not None:
        coot.add_alt_conf_py(imol, chain_id, res, ins, "", 0)
    return f"Added an alternate conformer to {chain_id}/{res} of model {imol}"


@command(r"(?:auto[- ]?fit|autofit)(?: rotamer)?" + OPT_RES_SPEC +
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["autofit A/72", "auto-fit rotamer A 72", "autofit"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Auto-fits the best rotamer for the residue against the "
               "refinement map (open a map first). " + ACTIVE_RESIDUE_NOTE)
def autofit_rotamer(chain: Optional[str] = None, resno: Optional[str] = None,
                    model: Optional[str] = None) -> str:
    """Auto-fit the best rotamer for a residue."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    if coot is None:
        return f"Auto-fitted {chain_id}/{res} of model {imol}"
    imol_map = coot.imol_refinement_map()
    if imol_map < 0:
        raise CommandError("no map set for fitting - open a map first")
    score = coot.auto_fit_best_rotamer(imol, chain_id, res, ins, "",
                                       imol_map, 1, 0.01)
    return (f"Auto-fitted {chain_id}/{res} of model {imol} "
            f"(rotamer score {score:.3g})")


@command(r"backrub(?: rotamer)?" + OPT_RES_SPEC +
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["backrub rotamer A/89", "backrub A 89", "backrub rotamer"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Applies a backrub rotamer fit to the named residue. "
               + ACTIVE_RESIDUE_NOTE)
def backrub_residue(chain: Optional[str] = None, resno: Optional[str] = None,
                    model: Optional[str] = None) -> str:
    """Apply a backrub rotamer fit to a residue."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    if coot is not None:
        coot.backrub_rotamer(imol, chain_id, res, ins, "")
    return f"Backrub-fitted {chain_id}/{res} of model {imol}"


@command(r"add ?oxt"
         r"(?: to| on)?"
         r"(?: (?:the )?active(?: residue)?| here)?"
         r"(?: (?:residue )?" + RES_SPEC + r")?"
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["add OXT to A/89", "add OXT to A 89", "add OXT"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Adds a terminal OXT oxygen to a residue (usually a chain's "
               "C-terminus). With no residue named, acts on the active "
               "residue (the one at the centre of the screen); with no model "
               "number, the active model.")
def add_oxt(chain: Optional[str] = None, resno: Optional[str] = None,
            model: Optional[str] = None) -> str:
    """Add a terminal OXT atom to a residue."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    if coot is not None:
        coot.add_OXT_to_residue(imol, chain_id, res, ins)
    return f"Added OXT to {chain_id}/{res} of model {imol}"


@command(r"add terminal residue"
         r"(?: to| on| onto)?"
         r"(?: (?:the )?active(?: residue)?| here)?"
         r"(?: (?:residue )?" + RES_SPEC + r")?"
         r"(?: (?:of |in )?model (?P<model>\S+))?"
         r"(?: (?:as|type) (?P<restype>[A-Za-z]{1,3}))?",
         examples=["add terminal residue to A/89", "add terminal residue",
                   "add terminal residue to A 89 as ALA"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Builds a new residue onto the end of a chain, attached to the "
               "named terminal residue and fitted against the refinement map. "
               "The type defaults to 'auto' (guessed from any sequence); add "
               "'as ALA' to force one. With no residue named, acts on the "
               "active residue; with no model number, the active model.")
def add_terminal_residue(chain: Optional[str] = None,
                         resno: Optional[str] = None,
                         model: Optional[str] = None,
                         restype: Optional[str] = None) -> str:
    """Add a terminal residue onto the end of a chain."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    rtype = restype.upper() if restype else "auto"
    suffix = f" ({rtype})" if rtype != "auto" else ""
    if coot is None:
        return f"Added a terminal residue{suffix} to {chain_id}/{res} of model {imol}"
    status = coot.add_terminal_residue(imol, chain_id, res, rtype, 1)
    if not status:
        return (f"Could not add a terminal residue to {chain_id}/{res} "
                f"of model {imol}")
    return f"Added a terminal residue{suffix} to {chain_id}/{res} of model {imol}"


@command(r"pep[- ]?flip" + OPT_RES_SPEC +
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["pepflip A/89", "pep flip A 89", "pepflip"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Flips the peptide bond following the named residue by 180. "
               + ACTIVE_RESIDUE_NOTE)
def pepflip_residue(chain: Optional[str] = None, resno: Optional[str] = None,
                    model: Optional[str] = None) -> str:
    """Flip the peptide following a residue."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    if coot is not None:
        centre_on_residue(imol, chain_id, res, ins)  # show it before flipping
        coot.pepflip(imol, chain_id, res, ins, "")
    return f"Flipped the peptide at {chain_id}/{res} of model {imol}"
