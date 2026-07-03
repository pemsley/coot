# coot_commands/__init__.py
#
# Copyright 2026 by Medical Research Council
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

"""Natural-language command package for Coot.

Importing this package auto-discovers and imports every module in the
``coot_commands.commands`` sub-package, which is what causes each
command's :func:`~coot_commands.registry.command` decorator to run and
register the command.  To add commands, drop a new module into
``coot_commands/commands/`` (and list it in ``python/Makefile.am``) - no
other file needs editing.

The public entry point is :func:`~coot_commands.registry.dispatch`, used
by the ``coot_command_interface`` shim that the C++ Command tab calls.
"""

import importlib
import pkgutil

from . import commands


def _discover_command_modules():
    """Import every module under coot_commands.commands so they register."""
    for _finder, name, _ispkg in pkgutil.iter_modules(commands.__path__):
        importlib.import_module(f"{__name__}.commands.{name}")


_discover_command_modules()
