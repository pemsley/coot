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
import re
import urllib.error
import urllib.request
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

from coot_commands.tools import command_tools, custom_tools, execute_tool

# Thinking models expose their reasoning in one of two ways: a dedicated field
# on the message (``reasoning`` / ``reasoning_content`` / ``thinking``, used by
# Ollama, DeepSeek-R1, ...) or an inline ``<think>...</think>`` block in the
# content (Qwen and friends). split_thinking() normalises both so callers can
# show the reasoning separately from the visible answer.
_THINK_RE = re.compile(r"<think>(.*?)</think>", re.DOTALL | re.IGNORECASE)


def split_thinking(message: Dict[str, Any]) -> Tuple[str, str]:
    """Return ``(thinking, visible_content)`` for a chat reply message.

    Pulls any reasoning field and any ``<think>`` block out of the content, so
    the reasoning can be surfaced on its own and the visible answer is clean.
    """
    thinking = (message.get("reasoning") or message.get("reasoning_content")
                or message.get("thinking") or "")
    content = message.get("content") or ""
    if isinstance(content, list):  # a multimodal content array: keep text parts
        content = " ".join(p.get("text", "") for p in content
                           if isinstance(p, dict) and p.get("type") == "text")
    blocks = _THINK_RE.findall(content)
    if blocks:
        inline = "\n".join(b.strip() for b in blocks)
        thinking = f"{thinking}\n{inline}".strip() if thinking else inline
        content = _THINK_RE.sub("", content)
    return thinking.strip(), content.strip()

DEFAULT_URL = "http://localhost:11434/v1/chat/completions"
DEFAULT_MODEL = "gemma4-long"
# How many commands to expose per request when retrieval is on.  A small model
# chooses far better from a handful of tools than from all ~90; see
# coot_commands.retrieval.  16 keeps the surface small while leaving margin so a
# common intent phrased with filler ("Oh no! go to A 89") still clears the cut.
DEFAULT_TOP_K = 16
# How many procedures (coot_commands.playbooks) to inject per request.  These
# are what tell the model what to *do* about a situation - fit a rotamer, then
# measure, then try a backrub - as opposed to what each command does.  Two is
# deliberate: one is usually the right procedure and the second covers a
# near-miss, while more would crowd out the conversation on a small model.
DEFAULT_PLAYBOOK_K = 2

# Commands exposed on every request, whatever retrieval picked.  Measuring a
# residue is relevant to any request that changes one - and a procedure that
# says "score it before and after" is worthless if score_residue was not among
# the retrieved tools.  refine_sphere earns its place for the same reason from
# the other end: settling the surroundings is the last step of almost every
# repair, and without it the model gets as far as "this now needs a local
# refinement" and has to stop.  Keep this list very short; it is a tax on
# every call.
PINNED_COMMANDS = ("score_residue", "refine_sphere")

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
    "Act rather than plan. Decide the single next call and make it - do not "
    "list the calls you intend to make, and do not restate a plan you have "
    "already worked out. If you find yourself going over the same decision a "
    "second time, that means it is time to call the tool, not to think again. "
    "When a step names two candidates, pick one, say which, and act on it.\n"
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
    "Deciding what to do. A request often names a problem rather than a "
    "command ('this rotamer is wrong', 'what should I do next'). When a "
    "procedure for the situation is provided, follow its steps rather than "
    "improvising. Measure before you act and again afterwards - score_residue "
    "reports a residue's rotamer probability and density fit - so you can tell "
    "whether a fix helped instead of assuming it did. Never report a fix as "
    "successful without a number showing it worked, and when something does "
    "not improve, say so and stop rather than trying the same move again.\n"
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


def _models_url(url: str) -> str:
    """Turn a chat endpoint (or a bare base) into Ollama's model-list URL.

    Ollama serves its catalogue from ``/api/tags`` on the server root, which
    sits alongside - not under - the OpenAI-compatible ``/v1`` prefix.
    """
    base = _normalise_chat_url(url)
    base = base[: -len("/v1/chat/completions")]
    return base + "/api/tags"


def list_models(url: Optional[str] = None, timeout: float = 2.0) -> List[str]:
    """The models the local Ollama server currently has, name-sorted.

    Returns ``[]`` when the server is unreachable or answers with something
    unexpected, so a caller can fall back without special-casing errors. Used
    by the Assistant tab to fill its model dropdown; keep the timeout short
    because Coot's GUI thread blocks on this call.
    """
    url = url or os.environ.get("COOT_AGENT_URL", DEFAULT_URL)
    try:
        with urllib.request.urlopen(_models_url(url), timeout=timeout) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, OSError, ValueError):
        return []
    models = payload.get("models") if isinstance(payload, dict) else None
    if not isinstance(models, list):
        return []
    names = [m.get("name") for m in models
             if isinstance(m, dict) and isinstance(m.get("name"), str)]
    # Drop embedding-only models: they are pulled for coot_commands.retrieval,
    # not for chat, and picking one would only produce a confusing failure.
    # /api/tags gives us no capability flag, so this is a name heuristic.
    names = [n for n in names if "embed" not in n.lower()]
    return sorted(set(names))


def _merge_tool_call_deltas(accumulated: Dict[int, Dict[str, Any]],
                            deltas: List[Dict[str, Any]]) -> None:
    """Fold streamed ``tool_calls`` fragments into *accumulated*, keyed by index.

    Ollama happens to send a tool call complete in a single delta, but the
    OpenAI streaming protocol allows it to arrive in pieces - the id and name
    in one chunk and the JSON arguments a few characters at a time after it -
    and other servers do exactly that.  Merging by ``index`` and concatenating
    the argument text handles both, treating Ollama's whole-in-one form as the
    degenerate case rather than assuming it.
    """
    for delta in deltas:
        if not isinstance(delta, dict):
            continue
        index = delta.get("index", 0)
        call = accumulated.setdefault(
            index, {"id": "", "type": "function",
                    "function": {"name": "", "arguments": ""}})
        if delta.get("id"):
            call["id"] = delta["id"]
        if delta.get("type"):
            call["type"] = delta["type"]
        function = delta.get("function") or {}
        if function.get("name"):
            call["function"]["name"] = function["name"]
        if function.get("arguments"):
            call["function"]["arguments"] += function["arguments"]


def _assemble_stream(lines: Iterable[bytes],
                     on_delta: Optional[EventFn] = None) -> Dict[str, Any]:
    """Build one reply message from an OpenAI server-sent-event stream.

    Separated from the HTTP call so it can be tested against a canned stream,
    and so the two transports assemble their result identically: what this
    returns has the same shape as the non-streaming ``choices[0].message``,
    and the loop above it cannot tell which produced it.

    *on_delta* receives a ``thinking_delta`` event per reasoning fragment, so
    the GUI can show the model working instead of a long silence.  Only the
    reasoning is streamed, not the answer: the answer is rendered as Markdown
    when it arrives (see ``text_view_append_markdown`` in vte.cc), which needs
    whole lines, and it is short enough that it appears in one go anyway.
    """
    content, reasoning = [], []
    tool_calls: Dict[int, Dict[str, Any]] = {}
    for raw in lines:
        line = raw.decode("utf-8", "replace").strip()
        if not line.startswith("data:"):
            continue                       # blank separators and comments
        body = line[len("data:"):].strip()
        if body == "[DONE]":
            break
        try:
            chunk = json.loads(body)
        except json.JSONDecodeError:
            continue                       # a truncated frame is not fatal
        choices = chunk.get("choices") or []
        if not choices:
            continue
        delta = choices[0].get("delta") or {}
        piece = delta.get("reasoning") or delta.get("reasoning_content") or ""
        if piece:
            reasoning.append(piece)
            if on_delta is not None:
                on_delta({"type": "thinking_delta", "text": piece})
        if delta.get("content"):
            content.append(delta["content"])
        if delta.get("tool_calls"):
            _merge_tool_call_deltas(tool_calls, delta["tool_calls"])

    message: Dict[str, Any] = {"role": "assistant", "content": "".join(content)}
    if reasoning:
        message["reasoning"] = "".join(reasoning)
    if tool_calls:
        message["tool_calls"] = [tool_calls[i] for i in sorted(tool_calls)]
    return message


def _ollama_chat(model: str, url: str, timeout: float,
                 messages: List[Dict[str, Any]],
                 tools: List[Dict[str, Any]],
                 on_delta: Optional[EventFn] = None) -> Dict[str, Any]:
    """Default transport: one round-trip to an OpenAI-compatible endpoint.

    Streams when *on_delta* is given, purely so the reasoning can be shown as
    it is produced; the assembled reply is identical either way.
    """
    url = _normalise_chat_url(url)
    streaming = on_delta is not None
    payload = {
        "model": model,
        "messages": messages,
        "tools": tools,
        "tool_choice": "auto",
        "stream": streaming,
        # Low temperature: we want deterministic tool selection, not prose.
        "temperature": 0.0,
        # Ollama's OpenAI-compatible endpoint maps this to the native
        # "think" level for models that support variable thinking effort.
        "reasoning_effort": "medium",
    }
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=data, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if streaming:
                # Iterating the response yields lines as they arrive, which is
                # the whole point - reading it all first would reintroduce the
                # silence we are removing.
                return _assemble_stream(resp, on_delta)
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
#   {"type": "tools",     "names": [...]}
#   {"type": "playbooks", "names": [...]}
#   {"type": "step",      "tool": name, "args": {...}, "result": "..."}
#   {"type": "thinking",  "text": "..."}
#   {"type": "final",     "text": "..."}
#   {"type": "stopped",   "steps": n}
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


def _retrieved_tools(user_text: str, top_k: int, emit: EventFn,
                     pinned: Iterable[str] = ()) -> List[Dict[str, Any]]:
    """Tools for the commands most relevant to *user_text*, top_k of them.

    *pinned* names are added to the retrieved set whatever the ranking said -
    :data:`PINNED_COMMANDS` plus the commands the selected playbooks refer to.
    Unknown names are dropped by :func:`command_tools`, so a playbook that
    names a command that no longer exists degrades quietly rather than failing.

    Falls back to the full command set if retrieval fails (e.g. the embedding
    model is not pulled or the server is down) so a request never breaks just
    because the optional embeddings are unavailable.
    """
    from coot_commands import retrieval
    try:
        names = retrieval.select_tools(user_text, top_k)
    except Exception as e:  # noqa: BLE001 - retrieval is best-effort
        emit({"type": "tools", "names": [], "error": str(e)})
        return command_tools()
    for name in pinned:
        if name not in names:
            names.append(name)
    emit({"type": "tools", "names": names})
    return command_tools(names)


def _retrieved_playbooks(user_text: str, k: int,
                         emit: EventFn) -> Tuple[str, List[str]]:
    """``(guidance text, tool names to pin)`` for the procedures that fit.

    Best-effort, exactly like tool retrieval: the assistant works without
    playbooks, so an unreachable embedding server costs guidance rather than
    the whole request.
    """
    from coot_commands import playbooks
    try:
        names = playbooks.select_playbooks(user_text, k)
    except Exception as e:  # noqa: BLE001 - guidance is best-effort
        emit({"type": "playbooks", "names": [], "error": str(e)})
        return "", []
    emit({"type": "playbooks", "names": names})
    return playbooks.render(names), playbooks.tools_for(names)


# How many earlier user turns join the current one to form the retrieval query.
RETRIEVAL_CONTEXT_TURNS = 2


def _retrieval_query(user_text: str, messages: Optional[List[Dict[str, Any]]],
                     n_turns: int = RETRIEVAL_CONTEXT_TURNS) -> str:
    """The text that tools and procedures are ranked against.

    Not just this request.  Mid-conversation a request is often meaningless on
    its own - "oh no!! let's get 32 sorted ASAP" carries no hint that this is
    about rotamers, so ranking it alone retrieves noise (font size, zoom) and
    no procedure clears the relevance floor.  The task was established a turn
    or two earlier, so the last few user turns join the current one and the
    query describes what the user is actually working on.

    Only user turns are used: tool results and the model's own replies are
    long, repeat the vocabulary of whatever it just did, and would pull the
    ranking towards the last action rather than the current intent.
    """
    said = [m.get("content") for m in (messages or [])
            if m.get("role") == "user"]
    recent = [c for c in said[-n_turns:] if isinstance(c, str) and c.strip()]
    return " ".join(recent + [user_text])


def _set_playbook_message(messages: List[Dict[str, Any]], text: str) -> None:
    """Rewrite the leading system message to carry the current guidance.

    The procedures that matter change from request to request, so each turn
    replaces the last turn's guidance rather than adding to it - otherwise a
    small model's context fills with advice about problems already dealt with.

    The guidance has to live *inside* the first system message, not in a second
    one before the user turn.  Several chat templates - Qwen3's among them -
    raise "System message must be at the beginning" when a system message
    appears after the conversation has started, which is exactly where a
    per-request message lands from the second request onwards.  Folding it into
    message 0 keeps one system message at the front, which every template
    accepts, and puts the procedures where a small model reads instructions
    best anyway.
    """
    system = {"role": "system", "content": SYSTEM_PROMPT}
    if messages and messages[0].get("role") == "system":
        messages[0] = system
    else:
        messages.insert(0, system)
    if text:
        messages[0]["content"] = SYSTEM_PROMPT + "\n\n" + text


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
    if kind == "playbooks":
        names = event.get("names") or []
        if event.get("error"):
            return f"  (no procedures retrieved: {event['error']})"
        if not names:
            return "  (no procedures matched)"
        return f"  (procedures: {', '.join(names)})"
    if kind == "step":
        args = ", ".join(f"{k}={v!r}" for k, v in (event.get("args") or {}).items())
        return f"  -> {event.get('tool')}({args}): {event.get('result')}"
    if kind == "warning":
        return f"  (warning: {event.get('message')})"
    if kind == "final":
        return event.get("text", "")
    if kind == "stopped":
        return f"Stopped after {event.get('steps')} tool-calling rounds without a final answer."
    return json.dumps(event)


EMPTY_ANSWER_NOTE = (
    "The model finished without writing an answer. This usually means its "
    "context window filled up: Ollama loads a model with a 4096-token context "
    "by default, which the system prompt, the tool definitions and a few tool "
    "results together exceed. Restarting the Ollama server with "
    "OLLAMA_CONTEXT_LENGTH=16384 (or more) fixes it. What was done is listed "
    "in the steps above."
)


def _recover_empty_answer(messages: List[Dict[str, Any]], thinking: str,
                          chat: ChatFn, emit: EventFn) -> str:
    """Salvage a reply that has reasoning, or nothing at all, but no answer.

    A model that runs out of context stops mid-thought and returns empty
    content, which would otherwise surface as a blank answer after a run that
    actually did the work.  Ask once more, without tools and with an explicit
    instruction to write plain text - a short request that often fits where
    the previous one did not.  Failing that, hand back the reasoning, and
    failing that, say what has most likely gone wrong.
    """
    emit({"type": "warning", "message": "empty answer; asking again"})
    nudge = list(messages) + [{
        "role": "user",
        "content": "Write your answer now, as plain text: what did you do, "
                   "and what were the before and after numbers?"}]
    try:
        retry = chat(nudge, [])
        _, content = split_thinking(retry)
        if content:
            return content
    except Exception as e:  # noqa: BLE001 - the fallbacks below still apply
        emit({"type": "warning", "message": f"retry failed: {e}"})
    if thinking:
        return thinking
    return EMPTY_ANSWER_NOTE


def run_agent(user_text: str, *,
              model: Optional[str] = None,
              url: Optional[str] = None,
              tools: Optional[List[Dict[str, Any]]] = None,
              chat: Optional[ChatFn] = None,
              execute: Optional[ExecuteFn] = None,
              on_event: Optional[EventFn] = None,
              messages: Optional[List[Dict[str, Any]]] = None,
              top_k: Optional[int] = DEFAULT_TOP_K,
              playbook_k: Optional[int] = DEFAULT_PLAYBOOK_K,
              max_steps: int = 16,
              timeout: float = 300.0,
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
    *timeout* is urllib's, which applies per socket operation rather than to
    the request as a whole.  That distinction matters: without streaming the
    entire reply arrives in one read at the end, so any generation longer than
    the timeout fails outright, whereas a streamed reply resets it on every
    chunk and only a genuine stall trips it.
    *tools* overrides the exposed command set; when it is ``None`` and *top_k*
    is set, embedding retrieval narrows the ~90 commands to the *top_k* most
    relevant (pass ``top_k=None`` to expose them all).  *playbook_k* is how
    many procedures (:mod:`coot_commands.playbooks`) are retrieved and injected
    as guidance for the request, and whose commands are pinned into the tool
    set; pass ``0`` or ``None`` to send none.  *max_steps* caps the
    tool-calling rounds so a confused model cannot loop forever.
    """
    model = model or os.environ.get("COOT_AGENT_MODEL", DEFAULT_MODEL)
    url = url or os.environ.get("COOT_AGENT_URL", DEFAULT_URL)
    execute = execute or execute_tool
    emit = _make_emit(on_event, verbose)

    # Rank against the request *in context* - see _retrieval_query.
    query = _retrieval_query(user_text, messages)

    # Pick the procedures first: they decide which commands must be exposed,
    # since a procedure that names a command the model was not given is worse
    # than no procedure at all.
    guidance, playbook_commands = (
        _retrieved_playbooks(query, playbook_k, emit) if playbook_k
        else ("", []))
    if tools is None:
        pinned = list(PINNED_COMMANDS) + playbook_commands
        commands = (_retrieved_tools(query, top_k, emit, pinned) if top_k
                    else command_tools())
        # Custom context tools (e.g. get_active_residue) are always available,
        # so "here"/"this residue" can be resolved whatever the request says.
        tools = custom_tools() + commands
    # Whether this turn's reasoning already reached the caller a fragment at a
    # time, so the aggregate "thinking" event below can be skipped rather than
    # repeating the whole thing under it.  Reset per chat call, since only the
    # streaming transport sets it and an injected one never will.
    streamed = {"this_turn": False}
    if chat is None:
        def chat(messages, tools):
            streamed["this_turn"] = False

            def on_delta(event: Dict[str, Any]) -> None:
                streamed["this_turn"] = True
                emit(event)

            return _ollama_chat(model, url, timeout, messages, tools,
                                on_delta=on_delta)

    # Seed a fresh conversation, or continue a caller-supplied one (giving the
    # agent memory across requests); either way append this request's turn.
    # _set_playbook_message writes the system message, so it also does the
    # seeding: it puts one at the front whether or not there is one already.
    if messages is None:
        messages = []
    _set_playbook_message(messages, guidance)
    messages.append({"role": "user", "content": user_text})

    for _step in range(max_steps):
        message = chat(messages, tools)
        messages.append(message)
        thinking, content = split_thinking(message)
        if thinking and not streamed["this_turn"]:
            emit({"type": "thinking", "text": thinking})
        tool_calls = message.get("tool_calls")
        if not tool_calls:
            if not content:
                content = _recover_empty_answer(messages, thinking, chat, emit)
            emit({"type": "final", "text": content})
            return content
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
