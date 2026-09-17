# coot_commands/commands/models.py
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

"""Whole-model operations: merging and superposing models.

These act on entire molecules at once (unlike the per-residue edits in
:mod:`coot_commands.commands.model_edit`), so they always name their models
explicitly rather than defaulting to the active one.
"""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import as_int, require_model, ArgType, CommandError

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Models"


@command(r"merge model[s]? (?P<first>\d+) (?:and|into|with|\+) model[s]? (?P<second>\d+)",
         examples=["merge model 1 and model 2"],
         category=CATEGORY,
         notes="Merges the second model into the first, so the combined model "
               "keeps the first model's number.")
def merge_models(first: str, second: str) -> str:
    """Merge one model into another."""
    target = require_model(as_int(first, "model number"))
    source = require_model(as_int(second, "model number"))
    if target == source:
        raise CommandError("cannot merge a model into itself")
    if coot is not None:
        coot.merge_molecules_py([source], target)
    return f"Merged model {source} into model {target}"


@command(r"superpose (?:model )?(?P<source>\S+) (?:onto|on|to) "
         r"(?:model )?(?P<target>\S+)"
         r"(?P<inplace> in[- ]?place| moving| no[- ]?copy)?",
         examples=["superpose model 0 onto model 1",
                   "superpose model 0 onto model 1 in place"],
         category=CATEGORY,
         arg_types={"source": ArgType.MODEL, "target": ArgType.MODEL},
         notes="Superposes the source model onto the target (reference) model "
               "by secondary-structure matching (SSM). By default a superposed "
               "copy is made, leaving the source model untouched; add 'in "
               "place' to move the source model itself instead.")
def superpose_models(source: str, target: str,
                     inplace: Optional[str] = None) -> str:
    """Superpose one model onto another (SSM)."""
    imol_source = require_model(as_int(source, "model number"))
    imol_target = require_model(as_int(target, "model number"))
    if imol_source == imol_target:
        raise CommandError("cannot superpose a model onto itself")
    copy = inplace is None
    if coot is not None:
        # superpose(reference, moving, flag): flag 1 copies the moving
        # molecule to a new one, 0 moves it in place.
        coot.superpose(imol_target, imol_source, 1 if copy else 0)
    if copy:
        return (f"Superposed a copy of model {imol_source} onto "
                f"model {imol_target}")
    return f"Superposed model {imol_source} onto model {imol_target} (moved)"
