# test_coot_figure.py
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

"""Standalone tests for the figure agent (no Coot, no model server).

Run from the python/ directory:

    python3 test_coot_figure.py

Covers the pure-Python TGA->PNG converter, the appearance command tools, the
vision-in-the-loop agent (with injected render/chat/execute), the Anthropic
message adapter, and the figure_serve stdin/stdout bridge.  Everything runs
offline; nothing here needs Coot or a model server.  Also discoverable by
pytest.
"""

import io
import json
import struct
import zlib

import coot_commands  # noqa: F401  - triggers command discovery/registration
from coot_commands import figure, figure_serve
from coot_commands.registry import all_commands, unmatched_examples


# ---------------------------------------------------------------------------
#   TGA -> PNG
# ---------------------------------------------------------------------------

def _coot_tga(width, height, rows_top_to_bottom):
    """Build a Coot-style TGA (32bpp BGRA, bottom-origin) from top-to-bottom rows.

    *rows_top_to_bottom* is a list of rows, each a list of (r, g, b) pixels.
    """
    head = struct.pack("<9h", 0, 2, 0, 0, 0, 0, width, height, 32)
    body = bytearray()
    for row in reversed(rows_top_to_bottom):  # file stores bottom row first
        for (r, g, b) in row:
            body += bytes([b, g, r, 255])
    return head + bytes(body)


def _png_dims_and_rows(png):
    """Decode a (small) PNG we wrote: return (w, h, rows) of (r, g, b) pixels."""
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    w, h = struct.unpack(">II", png[16:24])
    off, idat = 8, b""
    while off < len(png):
        ln = struct.unpack(">I", png[off:off + 4])[0]
        if png[off + 4:off + 8] == b"IDAT":
            idat += png[off + 8:off + 8 + ln]
        off += 12 + ln
    raw = zlib.decompress(idat)
    stride = w * 3
    rows = []
    for y in range(h):
        line = raw[y * (stride + 1) + 1:y * (stride + 1) + 1 + stride]
        rows.append([tuple(line[x * 3:x * 3 + 3]) for x in range(w)])
    return w, h, rows


def test_tga_to_png_colours_and_origin():
    tga = _coot_tga(2, 2, [[(255, 0, 0), (0, 255, 0)],
                           [(0, 0, 255), (255, 255, 255)]])
    w, h, rows = _png_dims_and_rows(figure.tga_to_png(tga))
    assert (w, h) == (2, 2)
    # Origin is flipped correctly: the top row we supplied stays on top.
    assert rows[0] == [(255, 0, 0), (0, 255, 0)]
    assert rows[1] == [(0, 0, 255), (255, 255, 255)]


def test_tga_to_png_handles_24bpp():
    head = struct.pack("<9h", 0, 2, 0, 0, 0, 0, 1, 1, 24)
    tga = head + bytes([30, 20, 10])  # BGR
    w, h, rows = _png_dims_and_rows(figure.tga_to_png(tga))
    assert (w, h) == (1, 1)
    assert rows[0][0] == (10, 20, 30)


def test_tga_to_png_rejects_bad_input():
    for bad in (b"", b"\x00" * 5):
        try:
            figure.tga_to_png(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("expected ValueError for malformed TGA")


def test_downscale_caps_long_edge_and_averages():
    # 4x2 image; cap long edge at 2 -> factor 2 -> 2x1, each output = block mean.
    rows = [[(0, 0, 0), (100, 100, 100), (10, 10, 10), (30, 30, 30)],
            [(200, 200, 200), (100, 100, 100), (50, 50, 50), (90, 90, 90)]]
    tga = _coot_tga(4, 2, rows)
    w, h, out = _png_dims_and_rows(figure.tga_to_png(tga, max_px=2))
    assert (w, h) == (2, 1)
    # left block mean = (0+100+200+100)/4 = 100; right = (10+30+50+90)/4 = 45
    assert out[0][0] == (100, 100, 100)
    assert out[0][1] == (45, 45, 45)


def test_downscale_noop_when_small_enough():
    tga = _coot_tga(2, 2, [[(1, 2, 3), (4, 5, 6)], [(7, 8, 9), (10, 11, 12)]])
    w, h, _ = _png_dims_and_rows(figure.tga_to_png(tga, max_px=1024))
    assert (w, h) == (2, 2)  # unchanged


def test_data_uri_round_trips():
    png = figure.tga_to_png(_coot_tga(1, 1, [[(1, 2, 3)]]))
    uri = figure._data_uri(png)
    assert uri.startswith("data:image/png;base64,")
    import base64
    assert base64.b64decode(uri.split(",", 1)[1]) == png


# ---------------------------------------------------------------------------
#   Appearance command tools
# ---------------------------------------------------------------------------

def test_appearance_commands_registered():
    names = {c.name for c in all_commands() if c.category == "Appearance"}
    for expected in ("set_fancy_lighting", "set_ambient_occlusion",
                     "set_shadow_strength", "set_background", "set_fog",
                     "set_depth_of_field", "set_model_specular"):
        assert expected in names, expected


def test_appearance_examples_all_match():
    appearance = {c.name for c in all_commands() if c.category == "Appearance"}
    problems = [(n, e) for (n, e) in unmatched_examples() if n in appearance]
    assert problems == [], problems


def test_figure_tool_names_are_visual_only():
    names = set(figure.figure_tool_names())
    assert "set_ambient_occlusion" in names   # Appearance
    assert "show_ribbons" in names             # Representation
    # Never edits the model.
    assert not any(n.startswith(("delete", "mutate", "refine")) for n in names)


# ---------------------------------------------------------------------------
#   The loop
# ---------------------------------------------------------------------------

def _one_pixel_png():
    return figure.tga_to_png(_coot_tga(1, 1, [[(10, 20, 30)]]))


def test_run_figure_agent_applies_tools_then_finishes():
    calls = {"n": 0}

    def render():
        return _one_pixel_png()

    def chat(messages, tools):
        # The latest user turn must carry exactly one image (older ones stripped).
        img_turns = [m for m in messages if m.get("role") == "user"
                     and isinstance(m.get("content"), list)]
        assert len(img_turns) == 1
        # finish must be offered as a tool.
        assert any(t["function"]["name"] == "finish" for t in tools)
        calls["n"] += 1
        if calls["n"] == 1:
            return {"role": "assistant", "content": "Flat; enabling AO.",
                    "tool_calls": [{"id": "a", "type": "function", "function": {
                        "name": "set_ambient_occlusion",
                        "arguments": '{"state": "on"}'}}]}
        return {"role": "assistant", "content": "Great now.",
                "tool_calls": [{"id": "b", "type": "function", "function": {
                    "name": "finish", "arguments": '{"summary": "Clean figure."}'}}]}

    executed = []

    def execute(name, args):
        executed.append((name, args))
        return f"{name} ok"

    events = []
    out = figure.run_figure_agent("nice figure", render=render, chat=chat,
                                  execute=execute, on_event=events.append,
                                  max_iters=6)
    assert out == "Clean figure."
    # finish is handled by the loop, never routed to execute.
    assert executed == [("set_ambient_occlusion", {"state": "on"})]
    types = [e["type"] for e in events]
    assert types.count("render") == 2
    assert types[-1] == "final"


def test_history_is_threaded_and_reasoning_stripped():
    # A shared conversation across two goals gives the agent memory of what it
    # did, and the stored history must NOT contain the verbose <think> reasoning.
    def chat(messages, tools):
        return {"role": "assistant",
                "content": "<think>lots of private reasoning here</think>Enabling AO.",
                "tool_calls": [{"id": "a", "type": "function", "function": {
                    "name": "set_ambient_occlusion", "arguments": '{"state": "on"}'}}]}

    convo = []
    figure.run_figure_agent("goal one", render=_one_pixel_png, chat=chat,
                            execute=lambda n, a: "ok", messages=convo, max_iters=1)
    # System prompt + user(goal) + cleaned assistant + tool reply.
    assert convo[0]["role"] == "system"
    assistant_msgs = [m for m in convo if m.get("role") == "assistant"]
    assert assistant_msgs and assistant_msgs[0]["content"] == "Enabling AO."
    assert all("reasoning" not in str(m.get("content", "")) for m in assistant_msgs)

    # Second goal continues the SAME conversation (memory retained).
    before = len(convo)
    figure.run_figure_agent("goal two", render=_one_pixel_png, chat=chat,
                            execute=lambda n, a: "ok", messages=convo, max_iters=1)
    assert len(convo) > before
    assert convo.count({"role": "system", "content": figure.SYSTEM_PROMPT}) == 1
    goals = [m for m in convo if m.get("role") == "user"]
    assert len(goals) == 2  # both goals present in one thread


def test_run_figure_agent_two_noops_ends():
    # A chatty model that never calls a tool: one bare-text reply is nudged, a
    # second in a row ends the loop (so it can't quit after a single no-op).
    seen = []

    def chat(messages, tools):
        seen.append(1)
        return {"role": "assistant", "content": "It already looks perfect."}

    events = []
    out = figure.run_figure_agent("x", render=_one_pixel_png, chat=chat,
                                  execute=lambda n, a: "ok", max_iters=5,
                                  on_event=events.append)
    assert out == "It already looks perfect."
    assert len(seen) == 2  # nudged once, then terminated on the second no-op
    assert [e["type"] for e in events].count("render") == 2


def test_run_figure_agent_noop_then_tool_continues():
    # A no-op followed by a real tool call must NOT end the loop early.
    turns = {"n": 0}

    def chat(messages, tools):
        turns["n"] += 1
        if turns["n"] == 1:
            return {"role": "assistant", "content": "Let me think..."}
        if turns["n"] == 2:
            return {"role": "assistant", "content": "Enabling AO.",
                    "tool_calls": [{"id": "a", "type": "function", "function": {
                        "name": "set_ambient_occlusion",
                        "arguments": '{"state": "on"}'}}]}
        return {"role": "assistant", "content": "Done.",
                "tool_calls": [{"id": "b", "type": "function", "function": {
                    "name": "finish", "arguments": '{"summary": "ok"}'}}]}

    executed = []
    out = figure.run_figure_agent("x", render=_one_pixel_png, chat=chat,
                                  execute=lambda n, a: executed.append(n) or "ok",
                                  max_iters=6)
    assert out == "ok"
    assert executed == ["set_ambient_occlusion"]


def test_run_figure_agent_stops_at_cap():
    def chat(messages, tools):
        return {"role": "assistant", "content": "one more tweak",
                "tool_calls": [{"id": "a", "type": "function", "function": {
                    "name": "set_fog", "arguments": '{"state": "on"}'}}]}

    events = []
    out = figure.run_figure_agent("x", render=_one_pixel_png, chat=chat,
                                  execute=lambda n, a: "ok", max_iters=2,
                                  on_event=events.append)
    assert "Stopped after 2" in out
    assert [e["type"] for e in events].count("render") == 2


def test_strip_old_images_keeps_only_latest():
    messages = [
        {"role": "system", "content": "s"},
        {"role": "user", "content": [{"type": "text", "text": "first"},
                                     {"type": "image_url",
                                      "image_url": {"url": "data:image/png;base64,AAA"}}]},
        {"role": "assistant", "content": "ok"},
        {"role": "user", "content": [{"type": "text", "text": "second"},
                                     {"type": "image_url",
                                      "image_url": {"url": "data:image/png;base64,BBB"}}]},
    ]
    figure._strip_old_images(messages)
    assert isinstance(messages[1]["content"], str)  # first render collapsed
    assert "first" in messages[1]["content"]
    assert isinstance(messages[3]["content"], list)  # latest render kept
    assert any(p.get("type") == "image_url" for p in messages[3]["content"])


# ---------------------------------------------------------------------------
#   Anthropic adapter
# ---------------------------------------------------------------------------

def test_split_thinking_from_field_and_tags():
    from coot_commands.agent import split_thinking
    # dedicated reasoning field
    t, c = split_thinking({"reasoning": "I should enable AO.", "content": "Enabling AO."})
    assert t == "I should enable AO." and c == "Enabling AO."
    # inline <think> block (Qwen-style)
    t, c = split_thinking({"content": "<think>black scene, need lighting</think>Turning on lighting."})
    assert t == "black scene, need lighting" and c == "Turning on lighting."
    # neither
    t, c = split_thinking({"content": "just an answer"})
    assert t == "" and c == "just an answer"


def test_figure_loop_emits_thinking():
    def chat(messages, tools):
        return {"role": "assistant",
                "content": "<think>flat and dark</think>Enabling lighting.",
                "tool_calls": [{"id": "a", "type": "function", "function": {
                    "name": "finish", "arguments": '{"summary": "done"}'}}]}
    events = []
    figure.run_figure_agent("x", render=_one_pixel_png, chat=chat,
                            execute=lambda n, a: "ok", on_event=events.append)
    thinking = [e for e in events if e["type"] == "thinking"]
    critique = [e for e in events if e["type"] == "critique"]
    assert thinking and thinking[0]["text"] == "flat and dark"
    assert critique and critique[0]["text"] == "Enabling lighting."


def test_from_anthropic_captures_thinking():
    msg = figure._from_anthropic({"content": [
        {"type": "thinking", "thinking": "reasoning here"},
        {"type": "text", "text": "answer"}]})
    assert msg["reasoning"] == "reasoning here"
    assert msg["content"] == "answer"


def test_anthropic_translation_round_trip():
    messages = [
        {"role": "system", "content": "SYS"},
        {"role": "user", "content": [
            {"type": "text", "text": "improve"},
            {"type": "image_url", "image_url": {"url": "data:image/png;base64,QUJD"}}]},
        {"role": "assistant", "content": "ok", "tool_calls": [
            {"id": "x", "type": "function",
             "function": {"name": "set_fog", "arguments": '{"state": "on"}'}}]},
        {"role": "tool", "tool_call_id": "x", "content": "Fog on"},
    ]
    tools = figure.command_tools(["set_fog"]) + [figure.FINISH_TOOL]
    system, am, at = figure._to_anthropic(messages, tools)
    assert system == "SYS"
    blocks = [b for m in am for b in m["content"]]
    kinds = {b["type"] for b in blocks}
    assert {"image", "tool_use", "tool_result"} <= kinds
    assert {"set_fog", "finish"} == {t["name"] for t in at}

    back = figure._from_anthropic({"content": [
        {"type": "text", "text": "great"},
        {"type": "tool_use", "id": "y", "name": "finish",
         "input": {"summary": "done"}}]})
    assert back["tool_calls"][0]["function"]["name"] == "finish"
    assert json.loads(back["tool_calls"][0]["function"]["arguments"]) == {"summary": "done"}


# ---------------------------------------------------------------------------
#   figure_serve bridge
# ---------------------------------------------------------------------------

def test_parse_request_forms():
    assert figure_serve._parse_request("hello") == ("goal", "hello")
    assert figure_serve._parse_request('{"goal": "x"}') == ("goal", "x")
    assert figure_serve._parse_request('{"text": "y"}') == ("goal", "y")
    assert figure_serve._parse_request('{"max_iters": 9}') == ("max_iters", 9)
    assert figure_serve._parse_request("   ") is None


def test_serve_streams_events(monkeypatch=None):
    # Stub the Coot/model-facing pieces so serve runs fully offline.
    orig = (figure_serve.make_socket_renderer, figure_serve.make_socket_executor,
            figure_serve.make_vision_chat, figure_serve.run_figure_agent)
    try:
        figure_serve.make_socket_renderer = lambda c, d: (lambda: b"PNG")
        figure_serve.make_socket_executor = lambda c: (lambda n, a: "ok")
        figure_serve.make_vision_chat = lambda: "CHAT"

        seen = {"goals": [], "convo_ids": set()}

        def fake_run(goal, *, render, chat, execute, messages, save_dir,
                     save_index, max_iters, on_event, verbose):
            seen["goals"].append(goal)
            seen["convo_ids"].add(id(messages))  # same list threaded across goals
            messages.append({"role": "user", "content": goal})
            on_event({"type": "render", "iteration": 0,
                      "image": save_dir + "/figure_iter_00.png"})
            on_event({"type": "final", "text": "done " + goal, "iterations": 1})
        figure_serve.run_figure_agent = fake_run

        class FakeClient:
            def connect(self): pass
            def close(self): pass

        stdin = io.StringIO('{"max_iters": 8}\nmake it pretty\n'
                            '{"reset": true}\nnow different\n')
        stdout = io.StringIO()
        figure_serve.serve(stdin, stdout, client=FakeClient(),
                           save_dir="/tmp/figtest")
        assert seen["goals"] == ["make it pretty", "now different"]
        assert len(seen["convo_ids"]) == 2  # reset swapped in a new conversation
    finally:
        (figure_serve.make_socket_renderer, figure_serve.make_socket_executor,
         figure_serve.make_vision_chat, figure_serve.run_figure_agent) = orig

    events = [json.loads(l) for l in stdout.getvalue().splitlines()]
    types = [e["type"] for e in events]
    assert types[0] == "ready"
    assert any(e["type"] == "status" for e in events)
    assert any(e["type"] == "config" and e["max_iters"] == 8 for e in events)
    assert any(e["type"] == "render" for e in events)
    assert any(e["type"] == "reset" for e in events)
    assert any(e["type"] == "final" and e["text"] == "done make it pretty"
               for e in events)
    assert types[-1] == "done"


def _run():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    failures = 0
    for test in tests:
        try:
            test()
            print(f"PASS {test.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"FAIL {test.__name__}: {e}")
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    return failures == 0


if __name__ == "__main__":
    import sys
    sys.exit(0 if _run() else 1)
