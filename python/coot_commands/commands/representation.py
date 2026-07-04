# coot_commands/commands/representation.py
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

"""Commands that change how models are drawn: carbon colours, symmetry."""

from coot_commands.registry import command
from coot_commands.types import resolve_model

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Representation"


@command(r"(?:colou?r )?(?:carbons?|carbon colou?rs?) (?:of model (?P<model>\S+) )?grey|"
         r"grey carbons?(?: (?:for|of) model (?P<model2>\S+))?",
         examples=["colour carbons grey", "grey carbons"],
         category=CATEGORY,
         notes="Uses grey for carbon atoms. With no model number, acts on "
               "the active model.")
def grey_carbons(model=None, model2=None):
    """Use grey carbon colours for a model."""
    imol = resolve_model(model if model is not None else model2)
    if coot is not None:
        coot.set_use_grey_carbons_for_molecule(imol, 1)
    return f"Model {imol}: grey carbons"


@command(r"(?:colou?r )?(?:carbons?|carbon colou?rs?) (?:of model (?P<model>\S+) )?coloured|"
         r"colou?red carbons?(?: (?:for|of) model (?P<model2>\S+))?",
         examples=["colour carbons coloured", "coloured carbons"],
         category=CATEGORY,
         notes="Uses per-element carbon colouring. With no model number, "
               "acts on the active model.")
def coloured_carbons(model=None, model2=None):
    """Use element-coloured carbons for a model."""
    imol = resolve_model(model if model is not None else model2)
    if coot is not None:
        coot.set_use_grey_carbons_for_molecule(imol, 0)
    return f"Model {imol}: coloured carbons"


@command(r"show symmetry",
         examples=["show symmetry"],
         category=CATEGORY)
def show_symmetry(**_):
    """Show symmetry-related molecules."""
    if coot is not None:
        coot.set_show_symmetry_master(1)
    return "Showing symmetry"


@command(r"hide symmetry",
         examples=["hide symmetry"],
         category=CATEGORY)
def hide_symmetry(**_):
    """Hide symmetry-related molecules."""
    if coot is not None:
        coot.set_show_symmetry_master(0)
    return "Hiding symmetry"
