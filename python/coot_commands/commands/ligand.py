# coot_commands/commands/ligand.py
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

"""Commands for fitting ligands into density.

The workflow mirrors Coot's Ligand > Find Ligands dialog: fetch the
ligand's dictionary as a molecule with ``get_monomer``, register the
protein, map and ligand for the search, then run ``execute_ligand_search``
which returns the molecule numbers of the fitted solutions.
"""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import resolve_model, resolve_map, ArgType, CommandError

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Ligand"


def _fit_ligand(comp_id: str, imol_protein: int, imol_map: int,
                here: bool) -> str:
    """Fit the ligand ``comp_id`` into ``imol_map`` and report the result.

    ``here`` restricts the search to the current view centre (Ligand >
    Find Ligands "Here"); otherwise the whole map is clustered and the top
    sites searched.  We fetch the ligand's dictionary as a molecule via
    ``get_monomer`` (which pulls the standard restraints for known codes),
    clear any ligands left over from a previous search, then run it.
    """
    if coot is None:
        return f"Fitted ligand {comp_id}"

    imol_ligand = coot.get_monomer(comp_id)
    if imol_ligand < 0:
        raise CommandError(
            f"could not get a dictionary for ligand {comp_id!r} - "
            f"check the three-letter code, or import its dictionary first")

    coot.add_ligand_clear_ligands()
    coot.set_ligand_search_protein_molecule(imol_protein)
    coot.set_ligand_search_map_molecule(imol_map)
    coot.add_ligand_search_ligand_molecule(imol_ligand)
    coot.set_find_ligand_here_cluster(1 if here else 0)

    solutions = coot.execute_ligand_search_py()
    solutions = list(solutions) if isinstance(solutions, (list, tuple)) else []

    # The input ligand molecule was only needed as a search template; hide
    # it so it does not clutter the view alongside the fitted copies.
    coot.set_mol_displayed(imol_ligand, 0)

    where = "at the view centre" if here else f"across map {imol_map}"
    if not solutions:
        return f"No fit found for ligand {comp_id} {where}"
    return (f"Fitted ligand {comp_id} {where}: {len(solutions)} solution(s) "
            f"in molecule(s) {', '.join(str(s) for s in solutions)}")


@command(r"(?:fit|find|add) ligand (?P<comp>[A-Za-z0-9]{1,5})"
         r"(?: (?:into |to |in )?model (?P<model>\S+))?"
         r"(?: (?:into |to |in |against |using )?map (?P<map>\S+))?"
         r"(?P<here> here)?",
         examples=["fit ligand LIG", "fit ligand ATP here",
                   "find ligand NAG into model 0 map 1"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL, "map": ArgType.MAP},
         notes="Fits a ligand (by three-letter code) into density, the same "
               "as Ligand > Find Ligands. Uses the active model and the "
               "refinement map by default. Add 'here' to search only at the "
               "current view centre - handy after 'go to blob N'. Returns the "
               "molecule numbers of the fitted solutions.")
def fit_ligand(comp: str, model: Optional[str] = None,
               map: Optional[str] = None, here: Optional[str] = None) -> str:
    """Fit a ligand into density by its three-letter code."""
    imol_protein = resolve_model(model)
    imol_map = resolve_map(map)
    return _fit_ligand(comp.upper(), imol_protein, imol_map, bool(here))
