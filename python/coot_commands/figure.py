# coot_commands/figure.py
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

"""A vision-in-the-loop agent that tunes Coot's shaders into a nice figure.

Where :mod:`coot_commands.agent` drives Coot from *text*, this drives it
from what the scene actually *looks like*.  One iteration is:

    render the current view  ->  show the image to a multimodal model
    ->  it critiques the figure and calls appearance tools to improve it
    ->  apply those changes in the live Coot  ->  render again

repeating until the model is happy (it calls ``finish``) or a cap is hit.
The tool surface is the appearance/representation command registry (see
:mod:`coot_commands.commands.appearance`): SSAO, shadows, lighting, fog,
depth-of-field, background, materials, ribbons/surface, and so on.

Nothing here is Coot-specific except the *renderer* and the *executor*, both
injected: the renderer produces PNG bytes of the current view, the executor
runs a tool call.  In a live session both go over Coot's JSON-RPC socket
(see :func:`make_socket_renderer` and
:func:`coot_commands.socket_client.make_socket_executor`); in tests they are
plain callables, so the whole loop runs offline.

Vision backend
--------------
The image critique needs a *multimodal* model.  The transport is pluggable
(the ``chat`` argument), and :func:`make_vision_chat` builds a default from
the environment:

* ``COOT_VISION_BACKEND`` - ``ollama`` (default) or ``anthropic``
* ``COOT_VISION_MODEL``   - model name (default: the text agent's model)
* ``COOT_VISION_URL``     - endpoint base (default: the text agent's URL)
* ``ANTHROPIC_API_KEY`` / ``COOT_ANTHROPIC_API_KEY`` - key for the anthropic backend

The ``ollama`` backend speaks the OpenAI-compatible chat API (images as
``image_url`` data URIs); the ``anthropic`` backend talks to the Messages
API, translating the OpenAI-shaped messages/tools this loop uses to and from
Anthropic's format so the loop itself never has to care which is in use.
"""

from __future__ import annotations

import base64
import json
import os
import struct
import urllib.error
import urllib.request
import zlib
from typing import Any, Callable, Dict, List, Optional

from coot_commands.agent import split_thinking
from coot_commands.registry import all_commands
from coot_commands.tools import command_tools, execute_tool

# Categories whose commands the figure agent may call.  These are the visual
# knobs; deliberately excludes model-editing so the agent restyles, never
# rebuilds.
FIGURE_TOOL_CATEGORIES = ("Appearance", "Representation")

# A message-producing transport, same shape as coot_commands.agent.ChatFn:
# (messages, tools) -> the assistant reply message (OpenAI choices[0].message).
ChatFn = Callable[[List[Dict[str, Any]], List[Dict[str, Any]]], Dict[str, Any]]
# Produces a PNG of the current Coot view.
RenderFn = Callable[[], bytes]
# Runs a tool call by name, returning a result string (mirrors agent.ExecuteFn).
ExecuteFn = Callable[[str, Dict[str, Any]], str]
# Structured progress sink (see the event shapes emitted below).
EventFn = Callable[[Dict[str, Any]], None]


SYSTEM_PROMPT = (
    "You are a molecular-graphics art director working inside Coot, tuning how "
    "a structure is rendered so it becomes a beautiful, publication-quality "
    "figure. You are shown the CURRENT rendered image each turn and you improve "
    "it by calling the provided appearance tools (lighting, ambient occlusion, "
    "shadows, depth of field, fog, outline, background, materials, ribbon/"
    "surface representations).\n"
    "\n"
    "Work like a careful artist, one considered change at a time:\n"
    "1. Look at the image and briefly critique it in plain text - what is the "
    "subject, and what specifically holds it back (flat lighting, no depth, "
    "distracting background, too dark/bright, cluttered)?\n"
    "2. Call one or a few tools that address the biggest problem. Do not change "
    "everything at once - you must see the effect before deciding the next "
    "move.\n"
    "3. You will then be shown the new render. Reassess and continue.\n"
    "\n"
    "Principles for a strong figure: turn on fancy lighting first (most effects "
    "need it); use ambient occlusion for shape and depth; a little cast shadow "
    "grounds the model; a dark, plain background makes a bright model pop, while "
    "white suits line-art print figures; depth of field draws the eye to the "
    "subject; keep specular highlights subtle, not plastic; avoid over-dark "
    "sooty ambient occlusion or crushed gamma. Prefer restraint - a clean, "
    "well-lit figure beats a busy one.\n"
    "\n"
    "Keep your reasoning brief - a few sentences at most - then act; long "
    "deliberation wastes the response budget. "
    "Every turn you MUST take an action: either call at least one appearance "
    "tool to improve the figure, or call the finish tool. Never reply with only "
    "text and no tool call - a short critique is welcome, but it must accompany "
    "a tool call. When the figure genuinely looks great and further changes "
    "would not help, call finish with a one-line summary. Only call tools that "
    "are provided, with their defined arguments."
)

# Not a registry command: figure-only, handled directly by the loop so it is
# never exposed to the text Assistant.
FINISH_TOOL: Dict[str, Any] = {
    "type": "function",
    "function": {
        "name": "finish",
        "description": ("Call this when the figure looks great and no further "
                        "changes are needed."),
        "parameters": {
            "type": "object",
            "properties": {
                "summary": {"type": "string",
                            "description": "One sentence describing the final look."},
            },
            "required": ["summary"],
        },
    },
}


# ---------------------------------------------------------------------------
#   TGA -> PNG  (pure standard library; no Pillow dependency)
# ---------------------------------------------------------------------------
#
# Coot's screendump writes an uncompressed 32-bpp BGRA Targa with a bottom-left
# origin (see src/screendump-tga.cc: TGAhead = {0,2,0,0,0,0,w,h,32}, GL_BGRA).
# Vision APIs want PNG/JPEG, so we decode that specific TGA and re-encode PNG
# here - keeping the whole loop dependency-free, like the text agent.

def _png_chunk(tag: bytes, data: bytes) -> bytes:
    """One PNG chunk: length, tag, data, CRC over tag+data."""
    return (struct.pack(">I", len(data)) + tag + data +
            struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))


def _encode_png(width: int, height: int, rgb: bytes) -> bytes:
    """Encode top-to-bottom, 8-bit RGB bytes as a PNG."""
    raw = bytearray()
    stride = width * 3
    for y in range(height):
        raw.append(0)  # filter type 0 (none) for this scanline
        raw += rgb[y * stride:(y + 1) * stride]
    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    return (b"\x89PNG\r\n\x1a\n" +
            _png_chunk(b"IHDR", ihdr) +
            _png_chunk(b"IDAT", zlib.compress(bytes(raw), 6)) +
            _png_chunk(b"IEND", b""))


# Default cap on the long edge of the image sent to the model. A retina
# screendump is huge; at full size it floods the context window and slows the
# model (and can truncate a thinking model before it emits its tool call). ~1024
# keeps ample detail for an aesthetic critique. Override per-renderer, or via
# COOT_FIGURE_MAX_PX for the socket renderer.
DEFAULT_MAX_PX = 1024


def _tga_to_rgb(tga: bytes):
    """Decode a Coot screendump TGA to ``(rgb_bytes, width, height)``.

    Handles the 24- and 32-bpp uncompressed true-colour TGAs Coot writes,
    honouring the origin bit so the image is never upside down.
    """
    if len(tga) < 18:
        raise ValueError("TGA too short to contain a header")
    id_len = tga[0]
    image_type = tga[2]
    width = tga[12] | (tga[13] << 8)
    height = tga[14] | (tga[15] << 8)
    bpp = tga[16]
    top_origin = bool(tga[17] & 0x20)
    if image_type != 2:
        raise ValueError(f"unsupported TGA image type {image_type} "
                         "(expected 2, uncompressed true-colour)")
    if bpp not in (24, 32):
        raise ValueError(f"unsupported TGA depth {bpp} (expected 24 or 32)")
    nbytes = bpp // 8
    start = 18 + id_len
    pix = tga[start:start + width * height * nbytes]
    if len(pix) < width * height * nbytes:
        raise ValueError("TGA pixel data truncated")

    src_stride = width * nbytes
    rgb = bytearray(width * height * 3)
    dst_stride = width * 3
    for y in range(height):
        # File rows run bottom-to-top unless the origin bit is set.
        src_y = y if top_origin else (height - 1 - y)
        row = memoryview(pix)[src_y * src_stride:(src_y + 1) * src_stride]
        out = memoryview(rgb)[y * dst_stride:(y + 1) * dst_stride]
        out[0::3] = row[2::nbytes]  # R  (source is BGR(A))
        out[1::3] = row[1::nbytes]  # G
        out[2::3] = row[0::nbytes]  # B
    return bytes(rgb), width, height


def _downscale_rgb(rgb: bytes, width: int, height: int, max_px: int):
    """Box-average downscale RGB bytes so the long edge is <= *max_px*.

    Uses an integer shrink factor and averages each factor x factor block, which
    (unlike nearest-neighbour subsampling) preserves thin features such as bonds.
    Returns ``(rgb, width, height)`` unchanged when already small enough.
    """
    long_edge = max(width, height)
    if not max_px or long_edge <= max_px:
        return rgb, width, height
    factor = (long_edge + max_px - 1) // max_px  # ceil -> smallest shrink that fits
    nw, nh = width // factor, height // factor
    out = bytearray(nw * nh * 3)
    row_bytes = width * 3
    inv = 1.0 / (factor * factor)
    oi = 0
    for oy in range(nh):
        y0 = oy * factor
        for ox in range(nw):
            base = y0 * row_bytes + ox * factor * 3
            r = g = b = 0
            for dy in range(factor):
                p = base + dy * row_bytes
                for _ in range(factor):
                    r += rgb[p]; g += rgb[p + 1]; b += rgb[p + 2]
                    p += 3
            out[oi] = int(r * inv)
            out[oi + 1] = int(g * inv)
            out[oi + 2] = int(b * inv)
            oi += 3
    return bytes(out), nw, nh


def tga_to_png(tga: bytes, max_px: Optional[int] = None) -> bytes:
    """Convert a Coot screendump TGA to PNG bytes, optionally downscaled.

    When *max_px* is given, the image is box-average downscaled so its long edge
    is at most *max_px* pixels before encoding.
    """
    rgb, width, height = _tga_to_rgb(tga)
    if max_px:
        rgb, width, height = _downscale_rgb(rgb, width, height, max_px)
    return _encode_png(width, height, rgb)


def _data_uri(png: bytes) -> str:
    """A ``data:image/png;base64,...`` URI for a PNG byte string."""
    return "data:image/png;base64," + base64.b64encode(png).decode("ascii")


# ---------------------------------------------------------------------------
#   Rendering the live view over the socket
# ---------------------------------------------------------------------------

def make_socket_renderer(client, workdir: str,
                         tga_name: str = "coot_figure.tga",
                         max_px: Optional[int] = None) -> RenderFn:
    """Return a ``render()`` that screendumps the live Coot and returns PNG bytes.

    *client* is a connected :class:`~coot_commands.socket_client.CootSocketClient`.
    The screendump is written to *workdir* (which Coot and this process share on
    the same machine); ``screendump_image`` runs on Coot's main thread, so by
    the time the RPC call returns the file is complete.  The PNG is downscaled so
    its long edge is at most *max_px* pixels (default: ``COOT_FIGURE_MAX_PX`` or
    :data:`DEFAULT_MAX_PX`) - a full retina screendump is too large to send.
    """
    tga_path = os.path.join(workdir, tga_name)
    if max_px is None:
        try:
            max_px = int(os.environ.get("COOT_FIGURE_MAX_PX", DEFAULT_MAX_PX))
        except ValueError:
            max_px = DEFAULT_MAX_PX

    def render() -> bytes:
        # screendump_image forces a .tga extension, so pass one and it is used
        # verbatim; the call redraws twice before capture.
        client.exec_python("__import__('coot').screendump_image({!r})".format(tga_path))
        with open(tga_path, "rb") as fh:
            return tga_to_png(fh.read(), max_px=max_px)

    return render


# ---------------------------------------------------------------------------
#   Vision transports
# ---------------------------------------------------------------------------

def _http_json(url: str, payload: Dict[str, Any], headers: Dict[str, str],
               timeout: float) -> Dict[str, Any]:
    """POST JSON and return the parsed JSON reply, surfacing the server's error."""
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data,
                                 headers={"Content-Type": "application/json", **headers})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        detail = e.read().decode("utf-8", "replace").strip()
        raise RuntimeError(f"vision request to {url} failed (HTTP {e.code}): "
                           f"{detail}") from None


def _ollama_vision_chat(model: str, url: str, timeout: float,
                        messages: List[Dict[str, Any]],
                        tools: List[Dict[str, Any]]) -> Dict[str, Any]:
    """OpenAI-compatible multimodal chat (Ollama, llama.cpp server, ...)."""
    from coot_commands.agent import _normalise_chat_url
    payload = {"model": model, "messages": messages, "tools": tools,
               "tool_choice": "auto", "stream": False, "temperature": 0.3,
               "max_tokens": 81920}
    body = _http_json(_normalise_chat_url(url), payload, {}, timeout)
    return body["choices"][0]["message"]


def _to_anthropic(messages: List[Dict[str, Any]], tools: List[Dict[str, Any]]):
    """Translate OpenAI-shaped messages/tools into Anthropic Messages format.

    Returns ``(system, anthropic_messages, anthropic_tools)``.  Text and image
    content, assistant ``tool_calls`` and ``role: tool`` results are mapped to
    Anthropic content blocks; the system message becomes the top-level system.
    """
    system = ""
    out: List[Dict[str, Any]] = []
    for msg in messages:
        role = msg.get("role")
        if role == "system":
            system = msg.get("content") or ""
            continue
        if role == "tool":
            out.append({"role": "user", "content": [{
                "type": "tool_result",
                "tool_use_id": msg.get("tool_call_id", ""),
                "content": msg.get("content", ""),
            }]})
            continue
        content = msg.get("content")
        blocks: List[Dict[str, Any]] = []
        if isinstance(content, str) and content:
            blocks.append({"type": "text", "text": content})
        elif isinstance(content, list):
            for part in content:
                if part.get("type") == "text":
                    blocks.append({"type": "text", "text": part.get("text", "")})
                elif part.get("type") == "image_url":
                    uri = part["image_url"]["url"]
                    b64 = uri.split(",", 1)[1]
                    blocks.append({"type": "image", "source": {
                        "type": "base64", "media_type": "image/png", "data": b64}})
        for call in msg.get("tool_calls") or []:
            fn = call.get("function", {})
            try:
                args = json.loads(fn.get("arguments") or "{}")
            except json.JSONDecodeError:
                args = {}
            blocks.append({"type": "tool_use", "id": call.get("id", ""),
                           "name": fn.get("name", ""), "input": args})
        if blocks:
            out.append({"role": "assistant" if role == "assistant" else "user",
                        "content": blocks})
    anth_tools = [{"name": t["function"]["name"],
                   "description": t["function"].get("description", ""),
                   "input_schema": t["function"].get("parameters",
                                                     {"type": "object", "properties": {}})}
                  for t in tools]
    return system, out, anth_tools


def _from_anthropic(body: Dict[str, Any]) -> Dict[str, Any]:
    """Translate an Anthropic Messages reply back to an OpenAI-style message."""
    text_parts: List[str] = []
    thinking_parts: List[str] = []
    tool_calls: List[Dict[str, Any]] = []
    for block in body.get("content", []):
        kind = block.get("type")
        if kind == "text":
            text_parts.append(block.get("text", ""))
        elif kind == "thinking":
            thinking_parts.append(block.get("thinking", ""))
        elif kind == "tool_use":
            tool_calls.append({
                "id": block.get("id", ""),
                "type": "function",
                "function": {"name": block.get("name", ""),
                             "arguments": json.dumps(block.get("input", {}))},
            })
    message: Dict[str, Any] = {"role": "assistant",
                               "content": "\n".join(p for p in text_parts if p)}
    if thinking_parts:
        message["reasoning"] = "\n".join(p for p in thinking_parts if p)
    if tool_calls:
        message["tool_calls"] = tool_calls
    return message


def _anthropic_vision_chat(model: str, url: str, key: str, timeout: float,
                           messages: List[Dict[str, Any]],
                           tools: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Anthropic Messages API transport, in OpenAI-message clothing."""
    system, anth_messages, anth_tools = _to_anthropic(messages, tools)
    # A generous cap: a thinking model needs room for its reasoning AND the
    # answer/tool call - too small and it is cut off before it acts.
    payload: Dict[str, Any] = {"model": model, "max_tokens": 4096,
                               "messages": anth_messages, "tools": anth_tools}
    if system:
        payload["system"] = system
    base = url.rstrip("/")
    if base.endswith("/v1/messages"):
        endpoint = base
    elif base.endswith("/v1"):
        endpoint = base + "/messages"
    else:
        endpoint = base + "/v1/messages"
    headers = {"x-api-key": key, "anthropic-version": "2023-06-01"}
    return _from_anthropic(_http_json(endpoint, payload, headers, timeout))


def make_vision_chat(timeout: float = 120.0) -> ChatFn:
    """Build the default vision transport from the environment (see module docs)."""
    from coot_commands.agent import DEFAULT_MODEL, DEFAULT_URL
    backend = os.environ.get("COOT_VISION_BACKEND", "ollama").lower()
    model = os.environ.get("COOT_VISION_MODEL", DEFAULT_MODEL)
    if backend == "anthropic":
        key = (os.environ.get("COOT_ANTHROPIC_API_KEY")
               or os.environ.get("ANTHROPIC_API_KEY", ""))
        url = os.environ.get("COOT_VISION_URL", "https://api.anthropic.com")

        def chat(messages, tools):
            return _anthropic_vision_chat(model, url, key, timeout, messages, tools)
        return chat

    url = os.environ.get("COOT_VISION_URL") or os.environ.get("COOT_AGENT_URL", DEFAULT_URL)

    def chat(messages, tools):
        return _ollama_vision_chat(model, url, timeout, messages, tools)
    return chat


# ---------------------------------------------------------------------------
#   The loop
# ---------------------------------------------------------------------------

def figure_tool_names() -> List[str]:
    """Registry command names the figure agent may call (visual categories)."""
    return [c.name for c in all_commands() if c.category in FIGURE_TOOL_CATEGORIES]


def _strip_old_images(messages: List[Dict[str, Any]]) -> None:
    """Replace images in all but the most recent user turn with a text note.

    Local multimodal models cope best with a single image in context, and it is
    the *latest* render that matters; older ones only cost tokens.  We rewrite
    earlier image turns to a short placeholder in place.
    """
    image_turns = [m for m in messages
                   if m.get("role") == "user" and isinstance(m.get("content"), list)]
    for msg in image_turns[:-1]:
        texts = [p.get("text", "") for p in msg["content"] if p.get("type") == "text"]
        msg["content"] = (" ".join(t for t in texts if t)
                          + " [earlier render omitted]").strip()


def _run_tool_calls(tool_calls, execute, emit) -> tuple:
    """Execute non-finish tool calls; return (tool_reply_messages, finish_summary)."""
    replies: List[Dict[str, Any]] = []
    finish_summary: Optional[str] = None
    for call in tool_calls:
        fn = call.get("function", {})
        name = fn.get("name", "")
        raw = fn.get("arguments") or "{}"
        try:
            args = json.loads(raw) if isinstance(raw, str) else raw
        except json.JSONDecodeError:
            args = {}
        if name == "finish":
            finish_summary = args.get("summary", "Looks good.")
            replies.append({"role": "tool", "tool_call_id": call.get("id", ""),
                            "content": "finished"})
            continue
        result = execute(name, args)
        emit({"type": "step", "tool": name, "args": args, "result": result})
        replies.append({"role": "tool", "tool_call_id": call.get("id", ""),
                        "content": result})
    return replies, finish_summary


def run_figure_agent(goal: str, *,
                     render: RenderFn,
                     chat: Optional[ChatFn] = None,
                     execute: Optional[ExecuteFn] = None,
                     on_event: Optional[EventFn] = None,
                     tool_names: Optional[List[str]] = None,
                     messages: Optional[List[Dict[str, Any]]] = None,
                     save_dir: Optional[str] = None,
                     save_index: int = 0,
                     max_iters: int = 6,
                     verbose: bool = False) -> str:
    """Iteratively tune Coot's appearance toward *goal*, judging by the render.

    *render* returns PNG bytes of the current view; *chat* is the multimodal
    transport (default: :func:`make_vision_chat`); *execute* runs a tool call
    (default: in-process :func:`~coot_commands.tools.execute_tool`; the live
    session injects a socket executor).  Each iteration is reported through
    *on_event* (and printed when *verbose*), including a ``render`` event whose
    ``image`` is the path of the saved PNG when *save_dir* is given.  Returns the
    model's closing summary.

    *messages* is the running conversation: pass the same list across goals to
    give the agent memory of what it already changed (seeded with the system
    prompt if empty, appended to in place); omit it for a one-shot run.  Only the
    model's *visible* answer and its tool calls/results are kept in history - the
    verbose reasoning is dropped so the context stays small.  *save_index* is the
    starting number for saved iteration PNGs, so filenames stay unique across
    goals in one session.
    """
    chat = chat or make_vision_chat()
    execute = execute or execute_tool
    names = tool_names if tool_names is not None else figure_tool_names()
    tools = command_tools(names) + [FINISH_TOOL]

    def emit(event: Dict[str, Any]) -> None:
        if on_event is not None:
            on_event(event)
        if verbose:
            print(event.get("type"), {k: v for k, v in event.items()
                                      if k not in ("type", "image")})

    emit({"type": "tools", "names": names})
    if messages is None:
        messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    elif not messages:
        messages.append({"role": "system", "content": SYSTEM_PROMPT})

    # A bare-text reply (no tool call, no finish) is NOT treated as "done" - a
    # chatty model would otherwise quit after one change. We nudge it to act or
    # finish, and only give up after two such no-ops in a row so a model that
    # simply won't call tools still terminates.
    noops = 0

    for i in range(max_iters):
        png = render()
        image_path = ""
        if save_dir:
            image_path = os.path.join(save_dir, f"figure_iter_{save_index + i:03d}.png")
            with open(image_path, "wb") as fh:
                fh.write(png)
        emit({"type": "render", "iteration": i, "image": image_path,
              "bytes": len(png)})

        if i == 0:
            prompt = (f"Goal: {goal}\n\nThis is the current Coot render. Critique "
                      "it, then call one or more appearance tools to improve it.")
        elif noops:
            prompt = ("You replied without changing anything. Either call an "
                      "appearance tool now to make a concrete improvement toward "
                      "the goal, or call finish if the figure is genuinely done.")
        else:
            prompt = ("Here is the updated render after your changes. Reassess "
                      "and call more tools to keep improving it, or call finish "
                      "if it now looks great.")
        messages.append({"role": "user", "content": [
            {"type": "text", "text": prompt},
            {"type": "image_url", "image_url": {"url": _data_uri(png)}},
        ]})
        _strip_old_images(messages)

        message = chat(messages, tools)
        thinking, critique = split_thinking(message)
        tool_calls = message.get("tool_calls")
        # Store a CLEANED reply in history: the visible answer and any tool calls,
        # but NOT the (often huge) reasoning - keeping it would bloat the context
        # every turn and starve the next reply of room, truncating it.
        cleaned: Dict[str, Any] = {"role": "assistant", "content": critique}
        if tool_calls:
            cleaned["tool_calls"] = tool_calls
        messages.append(cleaned)
        if thinking:
            emit({"type": "thinking", "iteration": i, "text": thinking})
        if critique:
            emit({"type": "critique", "iteration": i, "text": critique})

        if not tool_calls:
            noops += 1
            if noops >= 2:
                emit({"type": "final", "text": critique or "Done.",
                      "iterations": i + 1})
                return critique or "Done."
            continue  # nudge and re-render next iteration
        noops = 0

        replies, finish_summary = _run_tool_calls(tool_calls, execute, emit)
        messages.extend(replies)
        if finish_summary is not None:
            emit({"type": "final", "text": finish_summary, "iterations": i + 1})
            return finish_summary

    emit({"type": "stopped", "iterations": max_iters})
    return f"Stopped after {max_iters} iterations without finishing."


def main(argv: Optional[List[str]] = None) -> int:
    """CLI: drive a live Coot over its socket. Requires the RPC listener running.

    Usage::  python3 -m coot_commands.figure "a clean figure of the active model"
    """
    import sys
    import tempfile
    from coot_commands.socket_client import CootSocketClient, make_socket_executor

    args = sys.argv[1:] if argv is None else argv
    if not args:
        sys.stderr.write('usage: python3 -m coot_commands.figure "<goal>"\n')
        return 2
    workdir = tempfile.mkdtemp(prefix="coot_figure_")
    client = CootSocketClient()
    client.connect()
    render = make_socket_renderer(client, workdir)
    run_figure_agent(" ".join(args), render=render,
                     execute=make_socket_executor(client),
                     save_dir=workdir, verbose=True)
    client.close()
    print(f"(iteration images saved in {workdir})")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
