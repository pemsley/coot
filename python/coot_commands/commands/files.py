# coot_commands/commands/files.py
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

"""Load files from disk by path.

Reads a coordinate file (PDB/mmCIF), a map (CCP4/MRC), reflection data
(MTZ) or a restraints dictionary straight from a file path, dispatching on
the file extension - the same readers used by the File > Open menu.  A path
containing a shell wildcard (``*``, ``?`` or ``[...]``) is expanded and every
match is loaded.  For fetching structures from online databases by accession
code, see :mod:`coot_commands.commands.fetch`.
"""

from __future__ import annotations

import glob
import os

from coot_commands.registry import command
from coot_commands.types import ArgType, CommandError, molecule_name

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Files"

# Extensions we recognise, grouped by the reader they map to.  The check is
# on the lower-cased name, and a trailing ".gz" is stripped first, so
# "model.pdb.gz" is treated as ".pdb".
_COORD_EXTS = (".pdb", ".ent", ".cif", ".mmcif", ".pdbx")
_MAP_EXTS = (".map", ".ccp4", ".mrc")
_MTZ_EXTS = (".mtz",)

# Characters that make a path a shell glob rather than a single file.
_GLOB_CHARS = "*?["


def _clean_path(path: str) -> str:
    """Strip surrounding quotes/whitespace and expand '~' and env vars."""
    path = path.strip().strip('"').strip("'")
    return os.path.expanduser(os.path.expandvars(path))


def _resolve_paths(raw: str) -> list[str]:
    """Turn a raw path argument into a list of existing files.

    A wildcard (``*``, ``?`` or ``[...]``) is expanded with :mod:`glob`, so
    "load /data/*.pdb" loads every matching file; a plain path yields a
    single-element list.  Raises :class:`CommandError` if nothing matches.
    """
    path = _clean_path(raw)
    if any(c in path for c in _GLOB_CHARS):
        matches = sorted(m for m in glob.glob(path) if os.path.isfile(m))
        if not matches:
            raise CommandError(f"no files match: {path}")
        return matches
    if not os.path.isfile(path):
        raise CommandError(f"no such file: {path}")
    return [path]


def _extension(path: str) -> str:
    """The lower-case extension, ignoring a trailing '.gz'."""
    lower = path.lower()
    if lower.endswith(".gz"):
        lower = lower[:-3]
    return os.path.splitext(lower)[1]


def _name(imol: int) -> str:
    """A short, path-less molecule name, falling back to 'molecule N'."""
    return molecule_name(imol) or f"molecule {imol}"


def _load_one(path: str) -> str:
    """Read a single coordinate/map/MTZ file, returning a report line.

    Raises :class:`CommandError` if the extension is unknown or the reader
    reports failure.  Assumes ``coot`` is available (callers handle the
    headless case).
    """
    ext = _extension(path)
    base = os.path.basename(path)

    if ext in _COORD_EXTS:
        imol = coot.handle_read_draw_molecule_with_recentre(path, 1)
        if imol < 0:
            raise CommandError(f"could not read coordinates from {path}")
        return f"Loaded {base} as model {imol} ({_name(imol)})"

    if ext in _MAP_EXTS:
        imol = coot.handle_read_ccp4_map(path, 0)
        if imol < 0:
            raise CommandError(f"could not read map from {path}")
        return f"Loaded {base} as map {imol}"

    if ext in _MTZ_EXTS:
        imols = coot.auto_read_make_and_draw_maps(path)
        loaded = [m for m in imols if isinstance(m, int) and m >= 0]
        if not loaded:
            raise CommandError(f"could not read maps from {path}")
        return f"Loaded maps from {base} as {loaded}"

    raise CommandError(
        f"don't know how to load '{ext or base}' - expected coordinates "
        "(.pdb/.cif), a map (.map/.ccp4/.mrc) or reflections (.mtz)")


# Restraints dictionary - matched before the generic loader below so an
# explicit "load dictionary foo.cif" doesn't get read as coordinates.
@command(r"(?:load|open|read|import) "
         r"(?:cif )?(?:dict(?:ionary)?|restraints?|monomer(?: library)?|geometry) "
         r"(?P<path>.+)",
         examples=["load dictionary /path/to/LIG.cif",
                   "load dictionary ~/dict/*.cif"],
         category=CATEGORY,
         arg_types={"path": ArgType.PATH},
         notes="Reads one or more restraints (monomer) dictionaries, the same "
               "as File > Import CIF dictionary. A wildcard loads every match.")
def load_dictionary(path: str) -> str:
    """Load a restraints (monomer) dictionary from a path."""
    paths = _resolve_paths(path)
    if coot is None:
        return "Would load dictionary " + ", ".join(paths)
    for p in paths:
        coot.read_cif_dictionary(p)
    names = ", ".join(os.path.basename(p) for p in paths)
    return f"Loaded dictionary {names}"


# Only fire when the argument looks like a path - the first token contains a
# '/', '~', '.' or a wildcard, or ends in a known extension. This keeps the
# verbs "open" and "read" from shadowing word commands like "open sequence" /
# "open validation", and lets "load tutorial" fall through to its own handler.
@command(r"(?:load|open|read) "
         r"(?=\S*[/~.*?\[]|\S+\.(?:pdb|ent|cif|mmcif|pdbx|map|ccp4|mrc|mtz)\b)"
         r"(?P<path>.+)",
         examples=["load /path/to/model.cif",
                   "load ~/structures/4hhb.pdb",
                   "load /path/to/*.pdb",
                   "load /path/to/data.mtz"],
         category=CATEGORY,
         arg_types={"path": ArgType.PATH},
         notes="Loads a file from disk by path, choosing the reader from the "
               "extension: coordinates (.pdb, .ent, .cif, .mmcif, .pdbx), a "
               "map (.map, .ccp4, .mrc) or reflection data (.mtz). A trailing "
               ".gz is handled, and a wildcard (e.g. '*.pdb') loads every "
               "matching file. For a restraints dictionary use 'load "
               "dictionary'; to fetch from an online database by accession "
               "code use 'fetch'.")
def load_file(path: str) -> str:
    """Load coordinate, map or MTZ file(s) from a path or wildcard."""
    paths = _resolve_paths(path)
    if coot is None:
        return "Would load " + ", ".join(paths)

    results: list[str] = []
    errors: list[str] = []
    for p in paths:
        try:
            results.append(_load_one(p))
        except CommandError as e:
            errors.append(str(e))

    if not results:
        raise CommandError("; ".join(errors) or "nothing loaded")
    message = "\n".join(results)
    if errors:
        message += "\n" + "\n".join(f"(skipped: {e})" for e in errors)
    return message
