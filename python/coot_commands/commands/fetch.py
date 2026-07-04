# coot_commands/commands/fetch.py
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

"""Commands that fetch structures from online databases.

These wrap the fetch helpers in :mod:`get_ebi` (the same code the File >
Fetch menu uses): the PDBe for coordinates and EDS maps, and PDB-REDO for
re-refined models and maps.
"""

from __future__ import annotations

from typing import Any, Optional

from coot_commands.registry import command
from coot_commands.types import CommandError

CATEGORY = "Fetch"


def _get_ebi() -> Any:
    """Import the get_ebi helper module, or raise a readable error."""
    try:
        import get_ebi
        return get_ebi
    except ImportError as e:
        raise CommandError(f"fetch is unavailable ({e})")


def _report(code: str, result: Any, what: str) -> str:
    """Turn a get_ebi return value into a user message.

    The helpers return a molecule number, a list of molecule numbers, or
    ``False`` on failure.
    """
    code = code.upper()
    if result is False or result is None:
        return f"Could not fetch {code} ({what})"
    if isinstance(result, (list, tuple)):
        loaded = [m for m in result if isinstance(m, int) and m >= 0]
        return f"Fetched {code} ({what}) into molecules {loaded}"
    if isinstance(result, int) and result >= 0:
        return f"Fetched {code} ({what}) as model {result}"
    return f"Fetched {code} ({what})"


@command(r"fetch (?:pdb[- ]?redo )?(?P<code>\w{4,})(?: from)? pdb[- ]?redo|"
         r"fetch pdb[- ]?redo (?P<code2>\w{4,})",
         examples=["fetch 1abc from pdb-redo", "fetch pdb-redo 1abc"],
         category=CATEGORY,
         notes="Fetches the re-refined model and maps from PDB-REDO.")
def fetch_pdb_redo(code: Optional[str] = None, code2: Optional[str] = None) -> str:
    """Fetch a re-refined model and maps from PDB-REDO."""
    acc = code if code is not None else code2
    result = _get_ebi().get_pdb_redo(acc.lower())
    return _report(acc, result, "PDB-REDO")


@command(r"fetch (?:eds )?maps? (?:for )?(?P<code>\w{4,})|"
         r"fetch (?:pdb )?(?P<code2>\w{4,}) (?:and|with|\+) maps?",
         examples=["fetch 1abc and map", "fetch map for 1abc"],
         category=CATEGORY,
         notes="Fetches coordinates and the 2Fo-Fc / Fo-Fc maps from the "
               "Electron Density Server (EDS).")
def fetch_pdb_and_map(code: Optional[str] = None, code2: Optional[str] = None) -> str:
    """Fetch coordinates and maps (via EDS)."""
    acc = code if code is not None else code2
    result = _get_ebi().get_eds_pdb_and_mtz(acc.lower())
    return _report(acc, result, "PDB + map")


@command(r"fetch (?:pdb )?(?P<code>\w{4,})",
         examples=["fetch 1abc", "fetch pdb 4hhb"],
         category=CATEGORY,
         notes="Fetches coordinates from the PDBe. Add \"and map\" to also "
               "fetch maps, or \"from pdb-redo\" for the re-refined version.")
def fetch_pdb(code: str) -> str:
    """Fetch coordinates from the PDBe."""
    result = _get_ebi().get_ebi_pdb(code.lower())
    return _report(code, result, "PDB")
