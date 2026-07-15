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

from coot_commands.tools import command_tools, execute_tool

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
    "macromolecular models into experimental density. Carry out the user's "
    "request by calling the provided tools; each tool is a Coot command. "
    "Molecules are referred to by integer number (models and maps share the "
    "numbering). When the user does not name a molecule, omit the argument and "
    "the command acts on the active one. Call one tool at a time, use the "
    "result of each call to decide the next, and when the task is done reply "
    "with a short plain-text summary of what you did. Do not invent tools or "
    "arguments that were not provided."
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


def _run_tool_calls(tool_calls: List[Dict[str, Any]],
                    verbose: bool) -> List[Dict[str, Any]]:
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
        result = execute_tool(name, args)
        if verbose:
            shown = ", ".join(f"{k}={v!r}" for k, v in args.items())
            print(f"  -> {name}({shown}): {result}")
        replies.append({
            "role": "tool",
            "tool_call_id": call.get("id", ""),
            "content": result,
        })
    return replies


def _retrieved_tools(user_text: str, top_k: int,
                     verbose: bool) -> List[Dict[str, Any]]:
    """Tools for the commands most relevant to *user_text*, top_k of them.

    Falls back to the full command set if retrieval fails (e.g. the embedding
    model is not pulled or the server is down) so a request never breaks just
    because the optional embeddings are unavailable.
    """
    from coot_commands import retrieval
    try:
        names = retrieval.select_tools(user_text, top_k)
        if verbose:
            print(f"  (retrieved {len(names)} tools: {', '.join(names)})")
        return command_tools(names)
    except Exception as e:  # noqa: BLE001 - retrieval is best-effort
        if verbose:
            print(f"  (retrieval unavailable: {e}; using all tools)")
        return command_tools()


def run_agent(user_text: str, *,
              model: Optional[str] = None,
              url: Optional[str] = None,
              tools: Optional[List[Dict[str, Any]]] = None,
              chat: Optional[ChatFn] = None,
              top_k: Optional[int] = DEFAULT_TOP_K,
              max_steps: int = 8,
              timeout: float = 120.0,
              verbose: bool = True) -> str:
    """Fulfil *user_text* by letting the model call Coot commands.

    Returns the model's final plain-text reply.  *chat* overrides the transport
    (used by the tests); by default a fresh Ollama transport is built from
    *model*/*url* (falling back to ``COOT_AGENT_MODEL``/``COOT_AGENT_URL`` then
    the module defaults).  *tools* overrides the exposed command set; when it is
    ``None`` and *top_k* is set, embedding retrieval narrows the ~90 commands to
    the *top_k* most relevant (pass ``top_k=None`` to expose them all).
    *max_steps* caps the tool-calling rounds so a confused model cannot loop
    forever.
    """
    model = model or os.environ.get("COOT_AGENT_MODEL", DEFAULT_MODEL)
    url = url or os.environ.get("COOT_AGENT_URL", DEFAULT_URL)
    if tools is None:
        tools = _retrieved_tools(user_text, top_k, verbose) if top_k else command_tools()
    if chat is None:
        def chat(messages, tools):
            return _ollama_chat(model, url, timeout, messages, tools)

    messages: List[Dict[str, Any]] = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": user_text},
    ]

    for _step in range(max_steps):
        message = chat(messages, tools)
        messages.append(message)
        tool_calls = message.get("tool_calls")
        if not tool_calls:
            return (message.get("content") or "").strip()
        messages.extend(_run_tool_calls(tool_calls, verbose))

    return ("Stopped after {} tool-calling rounds without a final answer."
            .format(max_steps))


def main(argv: Optional[List[str]] = None) -> int:
    import sys
    args = sys.argv[1:] if argv is None else argv
    if not args:
        sys.stderr.write('usage: python3 -m coot_commands.agent "<request>"\n')
        return 2
    print(run_agent(" ".join(args)))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
