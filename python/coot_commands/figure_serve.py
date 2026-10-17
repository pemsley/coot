# coot_commands/figure_serve.py
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

"""Run the figure agent as a subprocess that Coot's Figure tab drives.

This mirrors :mod:`coot_commands.agent_serve`: Coot spawns this process for
the Figure tab, writes one request per line to ``stdin`` and reads
newline-delimited JSON *events* from ``stdout``.  The slow multimodal model
calls happen here; each appearance change and each render is executed back in
the live Coot over its JSON-RPC socket, on Coot's main thread.

Protocol
--------
Requests (one JSON object per line on stdin)::

    {"goal": "a clean figure of the active model on a dark background"}
    {"max_iters": 8}          # optional, sets the iteration cap for later goals

(a bare line of text is also accepted and treated as the ``goal``).

Events (one JSON object per line on stdout)::

    {"type": "ready"}
    {"type": "status", "backend": "...", "model": "...", "rpc": true, ...}
    {"type": "tools",  "names": [...]}
    {"type": "render", "iteration": 0, "image": "/path/figure_iter_00.png"}
    {"type": "critique", "iteration": 0, "text": "..."}
    {"type": "step",   "tool": "...", "args": {...}, "result": "..."}
    {"type": "final",  "text": "...", "iterations": n}
    {"type": "error",  "message": "..."}
    {"type": "done"}                     # one per request, always last

The ``image`` paths point at PNGs the loop saves in a shared work directory,
so the GUI can load and display the current render each iteration.  The port
to reach Coot comes from ``COOT_RPC_PORT``; the vision backend/model/endpoint
come from the ``COOT_VISION_*`` / ``ANTHROPIC_API_KEY`` environment (see
:mod:`coot_commands.figure`).
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from typing import Any, Dict, Optional, TextIO

from coot_commands.figure import (make_socket_renderer, make_vision_chat,
                                   run_figure_agent)
from coot_commands.socket_client import CootSocketClient, make_socket_executor

DEFAULT_MAX_ITERS = 6


def _emit(out: TextIO, event: Dict[str, Any]) -> None:
    """Write one JSON event as a line and flush (the GUI reads line by line)."""
    out.write(json.dumps(event) + "\n")
    out.flush()


def _parse_request(line: str):
    """Classify a stdin *line*: ``("goal", str)``, ``("max_iters", int)`` or None."""
    line = line.strip()
    if not line:
        return None
    try:
        parsed = json.loads(line)
    except json.JSONDecodeError:
        return ("goal", line)  # tolerate a bare, unquoted goal line
    if isinstance(parsed, dict):
        if parsed.get("reset"):
            return ("reset", None)
        if "max_iters" in parsed:
            try:
                return ("max_iters", int(parsed["max_iters"]))
            except (TypeError, ValueError):
                return None
        goal = parsed.get("goal") or parsed.get("text")
        if isinstance(goal, str) and goal.strip():
            return ("goal", goal)
        return None
    if isinstance(parsed, str):
        return ("goal", parsed) if parsed.strip() else None
    return None


def _startup_status(client: CootSocketClient) -> Dict[str, Any]:
    """Probe the RPC socket and report the configured vision backend/model."""
    from coot_commands.agent import DEFAULT_MODEL
    backend = os.environ.get("COOT_VISION_BACKEND", "ollama").lower()
    model = os.environ.get("COOT_VISION_MODEL", DEFAULT_MODEL)
    rpc_ok, rpc_detail = True, ""
    try:
        client.connect()  # also warms the connection reused for tool calls
    except Exception as e:  # noqa: BLE001
        rpc_ok, rpc_detail = False, str(e)
    return {"type": "status", "backend": backend, "model": model,
            "rpc": rpc_ok, "rpc_detail": rpc_detail}


def serve(stdin: TextIO, stdout: TextIO, *,
          client: Optional[CootSocketClient] = None,
          save_dir: Optional[str] = None,
          startup_status: bool = True) -> None:
    """Read goals from *stdin*, stream figure-agent events to *stdout*, until EOF.

    One *conversation* is threaded across goals for the life of the process, so
    the agent remembers what it already changed (a ``{"reset": true}`` line
    starts a new figure).  *client* and *save_dir* are injectable for testing;
    by default a real socket client and a temp directory shared with Coot are used.
    """
    client = client or CootSocketClient()
    if save_dir is None:
        save_dir = tempfile.mkdtemp(prefix="coot_figure_")
    render = make_socket_renderer(client, save_dir)
    execute = make_socket_executor(client)
    chat = make_vision_chat()
    max_iters = DEFAULT_MAX_ITERS
    conversation: list = []
    saved = 0  # cumulative render count, so iteration PNGs get unique names

    _emit(stdout, {"type": "ready"})
    if startup_status:
        _emit(stdout, _startup_status(client))

    for line in stdin:
        request = _parse_request(line)
        if request is None:
            continue
        kind, value = request
        if kind == "reset":
            conversation = []
            _emit(stdout, {"type": "reset"})
            continue
        if kind == "max_iters":
            max_iters = max(1, value)
            _emit(stdout, {"type": "config", "max_iters": max_iters})
            continue

        try:
            run_figure_agent(value, render=render, chat=chat, execute=execute,
                             messages=conversation, save_dir=save_dir,
                             save_index=saved, max_iters=max_iters,
                             on_event=lambda e: _emit(stdout, e), verbose=False)
        except Exception as e:  # noqa: BLE001 - report any failure to the GUI
            _emit(stdout, {"type": "error", "message": str(e)})
        saved += max_iters  # reserve a filename range for the renders just done
        _emit(stdout, {"type": "done"})

    client.close()


def main() -> int:
    serve(sys.stdin, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
