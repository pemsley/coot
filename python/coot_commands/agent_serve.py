# coot_commands/agent_serve.py
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

"""Run the agent as a subprocess that Coot's Assistant tab drives.

This is the process Coot spawns for the Assistant tab.  It reads one request
per line from ``stdin`` and writes newline-delimited JSON *events* to
``stdout``, so the GUI can stream progress without ever blocking on the slow
model calls (those happen here, in this separate process).  Tool calls are
executed back in the live Coot over its socket
(:mod:`coot_commands.socket_client`), which runs them on Coot's main thread.

Protocol
--------
Requests (one JSON object per line on stdin)::

    {"text": "load the tutorial data and refine A 89"}

(a bare line of text is also accepted and treated as the ``text``).

Events (one JSON object per line on stdout)::

    {"type": "ready"}
    {"type": "tools",  "names": [...]}
    {"type": "step",   "tool": "...", "args": {...}, "result": "..."}
    {"type": "final",  "text": "..."}
    {"type": "error",  "message": "..."}
    {"type": "done"}                     # one per request, always last

The port to reach Coot comes from ``COOT_RPC_PORT`` (Coot sets it when it
spawns us); the model and endpoints come from the same ``COOT_AGENT_*`` /
``COOT_EMBED_*`` environment as the standalone CLI.
"""

from __future__ import annotations

import json
import os
import sys
import urllib.error
import urllib.request
from typing import Any, Dict, Optional, TextIO
from urllib.parse import urlparse

from coot_commands.agent import run_agent
from coot_commands.socket_client import CootSocketClient, make_socket_executor


def _emit(out: TextIO, event: Dict[str, Any]) -> None:
    """Write one JSON event as a line and flush (the GUI reads line by line)."""
    out.write(json.dumps(event) + "\n")
    out.flush()


def _parse_request(line: str):
    """Classify a stdin *line*: ``("reset", None)``, ``("text", str)`` or None.

    ``{"reset": true}`` starts a new conversation; ``{"text": "..."}`` (or a
    bare, unquoted line) is a request; anything else is skipped.
    """
    line = line.strip()
    if not line:
        return None
    try:
        parsed = json.loads(line)
    except json.JSONDecodeError:
        return ("text", line)  # tolerate a bare, unquoted request line
    if isinstance(parsed, dict):
        if parsed.get("reset"):
            return ("reset", None)
        text = parsed.get("text")
        if isinstance(text, str) and text.strip():
            return ("text", text)
        return None
    if isinstance(parsed, str):
        return ("text", parsed) if parsed.strip() else None
    return None


def _probe_ollama(timeout: float = 2.0):
    """Check the model server is reachable; return ``(ok, detail)``.

    A GET to the server root is enough - Ollama answers it, and any HTTP
    response (even an error status) proves the server is up.
    """
    from coot_commands.agent import DEFAULT_URL, _normalise_chat_url
    url = _normalise_chat_url(os.environ.get("COOT_AGENT_URL", DEFAULT_URL))
    parts = urlparse(url)
    base = f"{parts.scheme}://{parts.netloc}"
    try:
        with urllib.request.urlopen(base, timeout=timeout) as resp:
            resp.read(64)
        return True, ""
    except urllib.error.HTTPError as e:
        return True, f"HTTP {e.code}"  # the server responded, so it is reachable
    except Exception as e:  # noqa: BLE001 - any failure means unreachable
        return False, str(e)


def _startup_status(client: CootSocketClient) -> Dict[str, Any]:
    """Probe the RPC socket and the model server for a GUI readiness indicator."""
    from coot_commands.agent import DEFAULT_MODEL
    model = os.environ.get("COOT_AGENT_MODEL", DEFAULT_MODEL)
    rpc_ok, rpc_detail = True, ""
    try:
        client.connect()  # also warms the connection reused for tool calls
    except Exception as e:  # noqa: BLE001
        rpc_ok, rpc_detail = False, str(e)
    ollama_ok, ollama_detail = _probe_ollama()
    return {"type": "status", "model": model,
            "rpc": rpc_ok, "rpc_detail": rpc_detail,
            "ollama": ollama_ok, "ollama_detail": ollama_detail}


def _context_stats(conversation: list) -> Dict[str, Any]:
    """Approximate how much context the running conversation is using.

    We have no tokenizer here, so tokens are estimated at ~4 characters each
    over the serialised messages - enough for a GUI "how full is the context"
    indicator, labelled as approximate.
    """
    return {"messages": len(conversation),
            "approx_tokens": max(0, len(json.dumps(conversation)) // 4)}


def serve(stdin: TextIO, stdout: TextIO, *,
          client: Optional[CootSocketClient] = None,
          startup_status: bool = True) -> None:
    """Read requests from *stdin*, stream events to *stdout*, until EOF.

    A single *conversation* is threaded across requests for the life of the
    process, so the agent remembers earlier turns (a ``{"reset": true}`` line
    starts a new one).  *client* is injectable for testing; by default a real
    socket client to Coot is created (connecting lazily on the first tool call).
    On start it emits a ``status`` event (RPC + model reachability) for the GUI
    readiness indicator, unless *startup_status* is false.
    """
    client = client or CootSocketClient()
    execute = make_socket_executor(client)
    conversation: list = []
    _emit(stdout, {"type": "ready"})
    if startup_status:
        _emit(stdout, _startup_status(client))

    for line in stdin:
        request = _parse_request(line)
        if request is None:
            continue
        kind, text = request
        if kind == "reset":
            conversation = []
            _emit(stdout, {"type": "reset"})
            continue
        try:
            run_agent(text, messages=conversation, execute=execute,
                      on_event=lambda e: _emit(stdout, e), verbose=False)
        except Exception as e:  # noqa: BLE001 - report any failure to the GUI
            _emit(stdout, {"type": "error", "message": str(e)})
        _emit(stdout, {"type": "context", **_context_stats(conversation)})
        _emit(stdout, {"type": "done"})

    client.close()


def main() -> int:
    serve(sys.stdin, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
