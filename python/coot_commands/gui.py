# coot_commands/gui.py
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

"""GUI feedback for commands, via Coot's C++ (GTK4) dialogs.

Commands must NOT use the legacy Python GUI helpers (``coot_gui`` /
``interesting_things_gui``): those import ``gi.repository.Gtk``, which fails
to initialise in Coot's embedded interpreter. The functions here call the
C++ GTK4 dialogs exposed on the ``coot`` module instead, which work. Every
call is best-effort and never raises, so GUI feedback can't break a command.
"""

from __future__ import annotations

try:
    import coot
except ImportError:
    coot = None


def status(text: str) -> None:
    """Show a one-line message in Coot's status bar."""
    if coot is not None and text:
        try:
            coot.add_status_bar_text(text)
        except Exception:
            pass


def dialog(text: str) -> None:
    """Pop up a Coot information dialog (use sparingly)."""
    if coot is not None and text:
        try:
            coot.info_dialog(text)
        except Exception:
            pass
