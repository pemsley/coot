# coot_commands/commands/fix.py
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

"""Commands that try a fix, measure it, and keep it only if it helped.

Everywhere else in this package a command is one Coot operation.  The
commands here are *procedures*: they run several operations, score the
residue after each (see :mod:`coot_commands.scoring`) and keep whichever
result is best, reverting the rest.

They exist for the assistant's sake.  Asking a small local model to run
that loop itself - fit, measure, compare, decide, undo - takes half a dozen
tool-calling rounds, and every round is a chance to mis-read a number or
forget to check.  Encoding the loop in Python makes it deterministic: the
model's only decision is *which* procedure to run, and the returned string
narrates every step so it can report honestly on what happened.

That narration is the whole interface.  A handler here must never return a
bare "Done" - the string it returns is all the model will ever know about
the outcome, so it carries the before and after numbers and says plainly
when nothing worked.
"""

from __future__ import annotations

from typing import List, Optional

from coot_commands import scoring
from coot_commands.registry import command
from coot_commands.scoring import NO_ROTAMER, ResidueScore, compare, is_better
from coot_commands.types import (OPT_RES_SPEC, resolve_residue, ArgType,
                                 CommandError, ACTIVE_RESIDUE_NOTE)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Fixing"


def _require_refinement_map() -> int:
    """The refinement map, or a clear error - every fit here needs one."""
    if coot is None:
        return -1
    imol_map = coot.imol_refinement_map()
    if imol_map < 0:
        raise CommandError("no map set for fitting - open a map first")
    return imol_map


class _Attempt:
    """One candidate conformation: how we got there, and how good it is."""

    def __init__(self, label: str, score: ResidueScore, positions: List) -> None:
        self.label = label
        self.score = score
        self.positions = positions


@command(r"(?:fix|sort out|improve) (?:the )?rotamer" + OPT_RES_SPEC +
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["fix rotamer A/45", "fix rotamer A 45", "fix rotamer"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Tries the standard rotamer repair sequence and keeps whichever "
               "result is best: score the residue, auto-fit the best rotamer, "
               "score again, then apply a backrub and score once more. The "
               "conformation with the best score is kept and the others are "
               "reverted, so this never leaves the residue worse than it "
               "started. 'Best' means: out of rotamer-outlier territory first, "
               "then the highest density fit. Needs a refinement map. "
               + ACTIVE_RESIDUE_NOTE)
def fix_rotamer(chain: Optional[str] = None, resno: Optional[str] = None,
                model: Optional[str] = None) -> str:
    """Auto-fit then backrub a side chain, keeping whichever scores best."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    imol_map = _require_refinement_map()
    if coot is None:
        return f"Fixed the rotamer of {chain_id}/{res} of model {imol}"

    name = scoring.residue_name(imol, chain_id, res, ins)
    if not name:
        raise CommandError(f"no residue {chain_id}/{res} in model {imol}")
    if name in NO_ROTAMER:
        return (f"{chain_id}/{res} is {name}, which has no side-chain rotamer "
                f"- nothing to fit. Check the backbone instead (Ramachandran, "
                f"peptide flips).")

    def look(clashes: bool = False) -> ResidueScore:
        return scoring.score_residue(imol, chain_id, res, ins,
                                     imol_map=imol_map, clashes=clashes)

    # A backrub rotates the flanking peptides, so a faithful revert has to
    # restore the neighbours too, not just this residue.
    window = [res - 1, res, res + 1]

    def snapshot() -> List:
        return scoring.capture_positions(imol, chain_id, window, ins)

    start = _Attempt("as it was", look(clashes=True), snapshot())
    attempts = [start]
    steps = [f"Started at {start.score.describe()}"]

    # 1. Auto-fit: search the rotamer library for the best fit to the density.
    try:
        fit_score = coot.auto_fit_best_rotamer(imol, chain_id, res, ins, "",
                                               imol_map, 1, 0.01)
        after_autofit = _Attempt("auto-fit", look(), snapshot())
        attempts.append(after_autofit)
        steps.append(f"Auto-fit (search score {fit_score:.3g}): "
                     f"{compare(start.score, after_autofit.score)}")
    except Exception as e:  # noqa: BLE001 - report, then still try the backrub
        steps.append(f"Auto-fit failed ({e})")

    # 2. Backrub: a small backbone adjustment that often rescues a side chain
    #    the rotamer search could not place.  Applied on top of the auto-fit,
    #    which is how one would do it by hand.
    previous = attempts[-1]
    try:
        coot.backrub_rotamer(imol, chain_id, res, ins, "")
        after_backrub = _Attempt("backrub", look(), snapshot())
        attempts.append(after_backrub)
        steps.append(f"Backrub: {compare(previous.score, after_backrub.score)}")
    except Exception as e:  # noqa: BLE001
        steps.append(f"Backrub failed ({e})")

    # 3. Keep the best of them.  Ties go to the earlier attempt, so a fit that
    #    changes nothing measurable leaves the model as it was.
    best = start
    for attempt in attempts[1:]:
        if is_better(attempt.score, best.score):
            best = attempt
    if best is not attempts[-1]:
        scoring.restore_positions(imol, best.positions)

    final = look(clashes=True)
    if best is start:
        verdict = (f"Neither improved it, so {chain_id}/{res} was left as it "
                   f"was.")
    else:
        verdict = f"Kept the {best.label} result."
    return (f"{chain_id}/{res} {name} of model {imol}: " + " ".join(steps)
            + f". {verdict} Now at {final.describe()}")


@command(r"(?:score|measure|check) residue" + OPT_RES_SPEC +
         r"(?: (?:of |in )?model (?P<model>\S+))?",
         examples=["score residue A/45", "score residue A 45", "score residue"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Reports how good one residue looks right now: its rotamer "
               "probability (with the outlier threshold), its density fit "
               "against the refinement map, and its clash volume. Use it "
               "before and after a fit to see whether the fit helped. "
               + ACTIVE_RESIDUE_NOTE)
def score_residue(chain: Optional[str] = None, resno: Optional[str] = None,
                  model: Optional[str] = None) -> str:
    """Report the rotamer, density fit and clashes of one residue."""
    imol, chain_id, res, ins = resolve_residue(chain, resno, model)
    if coot is None:
        return f"Scored {chain_id}/{res} of model {imol}"
    if not scoring.residue_name(imol, chain_id, res, ins):
        raise CommandError(f"no residue {chain_id}/{res} in model {imol}")
    return scoring.score_residue(imol, chain_id, res, ins).describe()
