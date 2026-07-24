# coot_commands/agent.py
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

"""Drive Coot from natural language with a small local model.

This is the loop that ties a language model to the command registry: it
sends the user's request plus the command *tools* (see
:mod:`coot_commands.tools`) to a local, OpenAI-compatible chat endpoint,
runs whatever tool calls the model emits, feeds the results back, and
repeats until the model answers in plain text.

The default endpoint is Ollama's OpenAI-compatible API
(``http://localhost:11434/v1/chat/completions``); any server speaking that
protocol works.  A 7-8B instruct model with solid tool-calling - Qwen2.5-7B
or Qwen3-8B at 4-bit fits comfortably in 16 GB - is the target.  Override
the model and URL with the ``COOT_AGENT_MODEL`` and ``COOT_AGENT_URL``
environment variables, or the keyword arguments.

Only the Python standard library is used, so there is no new build
dependency.  The chat transport is injectable (the *chat* argument), which
keeps the loop itself testable without a running model server - and lets a
future in-process, GTK-thread-aware transport slot in without touching the
loop.

Usage inside Coot's Python tab::

    import coot_commands.agent as agent
    print(agent.run_agent("go to residue A 45 and colour it by chain"))

or standalone (drives Ollama; handlers no-op without Coot, so this also
exercises the loop end to end)::

    python3 -m coot_commands.agent "add a water near A 45"

Note: :func:`run_agent` is synchronous and blocks on the model call.  Called
from Coot's GUI thread it will freeze the display until it returns; wiring it
into the Command tab without blocking (a worker thread that marshals each tool
call back onto the main loop) is deliberately left as a follow-up.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import Any, Callable, Dict, List, Optional

from coot_commands.tools import command_tools, custom_tools, execute_tool

DEFAULT_URL = "http://localhost:11434/v1/chat/completions"
DEFAULT_MODEL = "gemma4"
# How many commands to expose per request when retrieval is on.  A small model
# chooses far better from ~12 tools than from all ~90; see coot_commands.retrieval.
DEFAULT_TOP_K = 12

# A message-producing transport: given the running message list and the tool
# definitions, return the assistant's reply message (the OpenAI
# ``choices[0].message`` dict, with an optional ``tool_calls`` list).
ChatFn = Callable[[List[Dict[str, Any]], List[Dict[str, Any]]], Dict[str, Any]]

SYSTEM_PROMPT = (
    "You are the assistant inside Coot, a program for building and refining "
    "macromolecular models (proteins, nucleic acids, ligands) into experimental "
    "density from X-ray crystallography or cryo-EM. You act by calling the "
    "provided tools; each tool is a Coot command. Work in small steps: call one "
    "tool at a time, use each result to decide the next, and finish with a short "
    "plain-text summary of what you did.\n"
    "\n"
    "Molecules: every model and map has an integer molecule number, shared "
    "across models and maps (e.g. model 0, map 1). When the user does not name a "
    "molecule, omit that argument so the command acts on the active molecule. "
    "Residues are referenced by chain and residue number, written 'A/45' or "
    "'A 45'. When the user says 'here', 'this residue', 'the current residue' or "
    "similar, call get_active_residue to find out which residue and model they "
    "mean before acting.\n"
    "\n"
    "Structural-biology terms - map the user's shorthand to the right command. "
    "Real-space refinement (RSR, 'refine') locally optimises atoms into the "
    "density. A rotamer is a side-chain conformation; fitting or fixing a rotamer "
    "picks the best-fitting one. A Ramachandran outlier is a residue with an "
    "unusual backbone phi/psi combination. A peptide flip ('pepflip') rotates a "
    "peptide bond by ~180 degrees to correct the backbone; a backrub is a small "
    "local backbone adjustment. A clash is atoms too close together; C-beta "
    "deviations and chiral-volume errors are geometry problems. ADPs (B-factors) "
    "describe atomic displacement; occupancy is the fraction of an atom present; "
    "an alt conf is an alternate conformation; OXT is the C-terminal oxygen. "
    "Waters are ordered solvent; a ligand or monomer is a bound small molecule. "
    "The refinement map is the map refinement uses; a 2Fo-Fc map shows density, "
    "while a difference (Fo-Fc) map shows model-vs-data disagreement - green "
    "(positive) peaks suggest missing atoms, red (negative) peaks suggest atoms "
    "that should not be there. Validation flags these problems so you can fix "
    "them.\n"
    "\n"
    "The conversation may span several requests: remember what you did earlier "
    "(for example, a residue you just refined) and use it as context. Only call "
    "tools that are provided, with the arguments they define - never invent a "
    "tool or argument. If a request is ambiguous or very large in scope, do the "
    "most sensible part and state what you assumed, or ask one brief clarifying "
    "question."
)


def _normalise_chat_url(url: str) -> str:
    """Accept a full endpoint or just a base, and return the chat endpoint.

    POSTing to the Ollama base URL (``http://localhost:11434``) returns 405
    Method Not Allowed, so we tolerate a base or a ``.../v1`` root and append
    the ``/v1/chat/completions`` path, and strip a trailing slash (which
    otherwise redirects).
    """
    url = url.rstrip("/")
    if url.endswith("/chat/completions"):
        return url
    if url.endswith("/v1"):
        return url + "/chat/completions"
    return url + "/v1/chat/completions"


def _ollama_chat(model: str, url: str, timeout: float,
                 messages: List[Dict[str, Any]],
                 tools: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Default transport: one round-trip to an OpenAI-compatible endpoint."""
    url = _normalise_chat_url(url)
    payload = {
        "model": model,
        "messages": messages,
        "tools": tools,
        "tool_choice": "auto",
        "stream": False,
        # Low temperature: we want deterministic tool selection, not prose.
        "temperature": 0.0,
    }
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=data, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        # The server's body carries the real reason (e.g. "model 'x' not
        # found"), which urllib otherwise hides behind a bare status code.
        detail = e.read().decode("utf-8", "replace").strip()
        raise RuntimeError(
            f"chat request to {url} with model {model!r} failed "
            f"(HTTP {e.code}): {detail}") from None
    return body["choices"][0]["message"]


# Executes a tool by name with a dict of arguments, returning a result string.
# The default runs commands in-process; the GUI/agent_serve path injects one
# that runs them over the socket into a live Coot (see coot_commands.socket_client).
ExecuteFn = Callable[[str, Dict[str, Any]], str]

# Receives structured progress events so a consumer (the GUI transcript, a test)
# can observe the run without parsing printed text.  Event shapes:
#   {"type": "tools",  "names": [...]}
#   {"type": "step",   "tool": name, "args": {...}, "result": "..."}
#   {"type": "final",  "text": "..."}
#   {"type": "stopped","steps": n}
EventFn = Callable[[Dict[str, Any]], None]


def _run_tool_calls(tool_calls: List[Dict[str, Any]],
                    execute: ExecuteFn,
                    emit: EventFn) -> List[Dict[str, Any]]:
    """Execute each tool call, returning the ``role: tool`` reply messages."""
    replies = []
    for call in tool_calls:
        function = call.get("function", {})
        name = function.get("name", "")
        raw_args = function.get("arguments") or "{}"
        try:
            args = json.loads(raw_args) if isinstance(raw_args, str) else raw_args
        except json.JSONDecodeError:
            args = {}
        result = execute(name, args)
        emit({"type": "step", "tool": name, "args": args, "result": result})
        replies.append({
            "role": "tool",
            "tool_call_id": call.get("id", ""),
            "content": result,
        })
    return replies


def _retrieved_tools(user_text: str, top_k: int,
                     emit: EventFn) -> List[Dict[str, Any]]:
    """Tools for the commands most relevant to *user_text*, top_k of them.

    Falls back to the full command set if retrieval fails (e.g. the embedding
    model is not pulled or the server is down) so a request never breaks just
    because the optional embeddings are unavailable.
    """
    from coot_commands import retrieval
    try:
        names = retrieval.select_tools(user_text, top_k)
        emit({"type": "tools", "names": names})
        return command_tools(names)
    except Exception as e:  # noqa: BLE001 - retrieval is best-effort
        emit({"type": "tools", "names": [], "error": str(e)})
        return command_tools()


def _make_emit(on_event: Optional[EventFn], verbose: bool) -> EventFn:
    """Build the event sink: fans out to *on_event* and/or a human-readable print."""
    def emit(event: Dict[str, Any]) -> None:
        if on_event is not None:
            on_event(event)
        if verbose:
            print(_format_event(event))
    return emit


def _format_event(event: Dict[str, Any]) -> str:
    """Render an event as the one-line form the CLI/verbose mode prints."""
    kind = event.get("type")
    if kind == "tools":
        names = event.get("names") or []
        if event.get("error"):
            return f"  (retrieval unavailable: {event['error']}; using all tools)"
        return f"  (retrieved {len(names)} tools: {', '.join(names)})"
    if kind == "step":
        args = ", ".join(f"{k}={v!r}" for k, v in (event.get("args") or {}).items())
        return f"  -> {event.get('tool')}({args}): {event.get('result')}"
    if kind == "final":
        return event.get("text", "")
    if kind == "stopped":
        return f"Stopped after {event.get('steps')} tool-calling rounds without a final answer."
    return json.dumps(event)


def run_agent(user_text: str, *,
              model: Optional[str] = None,
              url: Optional[str] = None,
              tools: Optional[List[Dict[str, Any]]] = None,
              chat: Optional[ChatFn] = None,
              execute: Optional[ExecuteFn] = None,
              on_event: Optional[EventFn] = None,
              messages: Optional[List[Dict[str, Any]]] = None,
              top_k: Optional[int] = DEFAULT_TOP_K,
              max_steps: int = 8,
              timeout: float = 120.0,
              verbose: bool = True) -> str:
    """Fulfil *user_text* by letting the model call Coot commands.

    Returns the model's final plain-text reply.  *chat* overrides the transport
    (used by the tests); by default a fresh Ollama transport is built from
    *model*/*url* (falling back to ``COOT_AGENT_MODEL``/``COOT_AGENT_URL`` then
    the module defaults).  *execute* overrides how a tool call is run (default:
    in-process :func:`coot_commands.tools.execute_tool`; the GUI injects a
    socket-backed executor into a live Coot).  *on_event* receives structured
    progress events (see :data:`EventFn`), for a GUI transcript or tests.
    *messages* is the running conversation: pass the same list across calls to
    give the agent memory of earlier requests (it is seeded with the system
    prompt if empty and appended to in place); omit it for a one-shot call.
    *tools* overrides the exposed command set; when it is ``None`` and *top_k*
    is set, embedding retrieval narrows the ~90 commands to the *top_k* most
    relevant (pass ``top_k=None`` to expose them all).  *max_steps* caps the
    tool-calling rounds so a confused model cannot loop forever.
    """
    model = model or os.environ.get("COOT_AGENT_MODEL", DEFAULT_MODEL)
    url = url or os.environ.get("COOT_AGENT_URL", DEFAULT_URL)
    execute = execute or execute_tool
    emit = _make_emit(on_event, verbose)
    if tools is None:
        commands = _retrieved_tools(user_text, top_k, emit) if top_k else command_tools()
        # Custom context tools (e.g. get_active_residue) are always available,
        # so "here"/"this residue" can be resolved whatever the request says.
        tools = custom_tools() + commands
    if chat is None:
        def chat(messages, tools):
            return _ollama_chat(model, url, timeout, messages, tools)

    # Seed a fresh conversation, or continue a caller-supplied one (giving the
    # agent memory across requests); either way append this request's turn.
    if messages is None:
        messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    elif not messages:
        messages.append({"role": "system", "content": SYSTEM_PROMPT})
    messages.append({"role": "user", "content": user_text})

    for _step in range(max_steps):
        message = chat(messages, tools)
        messages.append(message)
        tool_calls = message.get("tool_calls")
        if not tool_calls:
            text = (message.get("content") or "").strip()
            emit({"type": "final", "text": text})
            return text
        messages.extend(_run_tool_calls(tool_calls, execute, emit))

    emit({"type": "stopped", "steps": max_steps})
    return ("Stopped after {} tool-calling rounds without a final answer."
            .format(max_steps))


def main(argv: Optional[List[str]] = None) -> int:
    import sys
    args = sys.argv[1:] if argv is None else argv
    if not args:
        sys.stderr.write('usage: python3 -m coot_commands.agent "<request>"\n')
        return 2
    # verbose=True already prints the step and final lines as they happen.
    run_agent(" ".join(args), verbose=True)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
