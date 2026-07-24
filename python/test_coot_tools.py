# test_coot_tools.py
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

"""Standalone tests for the language-model bridge (no Coot, no model server).

Run from the python/ directory:

    python3 test_coot_tools.py

Covers coot_commands.tools (registry -> tool schemas, tool call -> handler)
and coot_commands.agent's loop with an injected fake transport, so nothing
here needs Coot or Ollama.  Also discoverable by pytest.
"""

import io
import json
import socket
import struct
import threading

import coot_commands  # noqa: F401  - triggers command discovery/registration
from coot_commands import agent
from coot_commands.registry import all_commands
from coot_commands.tools import command_to_tool, command_tools, execute_tool


def _tool_named(name):
    for tool in command_tools():
        if tool["function"]["name"] == name:
            return tool
    return None


def test_every_command_becomes_a_valid_tool():
    tools = command_tools()
    assert len(tools) == len({c.name for c in all_commands()})
    for tool in tools:
        assert tool["type"] == "function"
        fn = tool["function"]
        assert isinstance(fn["name"], str) and fn["name"]
        assert isinstance(fn["description"], str) and fn["description"]
        params = fn["parameters"]
        assert params["type"] == "object"
        assert isinstance(params["properties"], dict)


def test_tool_names_are_unique():
    names = [t["function"]["name"] for t in command_tools()]
    assert len(names) == len(set(names))


def test_required_and_optional_params_from_signature():
    # go_to_residue(chain, resno, model=None): chain/resno required, model not.
    tool = _tool_named("go_to_residue")
    assert tool is not None
    params = tool["function"]["parameters"]
    assert set(params["properties"]) == {"chain", "resno", "model"}
    assert set(params["required"]) == {"chain", "resno"}


def test_arg_type_enriches_description_and_omit_hint():
    tool = _tool_named("go_to_residue")
    model_desc = tool["function"]["parameters"]["properties"]["model"]["description"]
    assert "Model" in model_desc
    assert "active" in model_desc  # optional -> omit hint appended


def test_examples_included_in_description():
    tool = _tool_named("go_to_residue")
    assert "Example phrasings" in tool["function"]["description"]


def test_kwargs_only_handler_has_no_params():
    # next_residue(**_) takes no real arguments.
    tool = _tool_named("next_residue")
    assert tool is not None
    assert tool["function"]["parameters"]["properties"] == {}
    assert "required" not in tool["function"]["parameters"]


def test_execute_tool_dispatches_to_handler():
    # centre_at_xyz runs without Coot (the coot-is-None branch returns a string).
    out = execute_tool("centre_at_xyz", {"x": "12.0", "y": "4.5", "z": "-3.2"})
    assert out == "Centred at (12, 4.5, -3.2)"


def test_execute_tool_coerces_numeric_arguments_to_strings():
    # A model may emit JSON numbers, not strings; handlers expect strings.
    out = execute_tool("centre_at_xyz", {"x": 1, "y": 2, "z": 3})
    assert out == "Centred at (1, 2, 3)"


def test_execute_unknown_tool_reports_error():
    out = execute_tool("no_such_command", {})
    assert out.startswith("Error: unknown command")


def test_execute_tool_reports_command_error():
    # A bad coordinate raises CommandError inside the handler; execute_tool
    # turns it into a readable string rather than propagating.
    out = execute_tool("centre_at_xyz", {"x": "not-a-number", "y": "0", "z": "0"})
    assert out.startswith("Error:")


def test_command_tools_subset_by_name():
    tools = command_tools(names=["go_to_residue", "no_such_command"])
    names = [t["function"]["name"] for t in tools]
    assert names == ["go_to_residue"]  # unknown names are dropped


def test_command_to_tool_matches_registry_entry():
    cmd = next(c for c in all_commands() if c.name == "centre_at_xyz")
    tool = command_to_tool(cmd)
    assert tool["function"]["name"] == "centre_at_xyz"
    assert set(tool["function"]["parameters"]["properties"]) == {"x", "y", "z"}


# --- agent loop (fake transport) --------------------------------------------

def _fake_chat_script(*replies):
    """Return a chat transport that yields *replies* in order, recording calls."""
    state = {"i": 0, "seen": []}

    def chat(messages, tools):
        state["seen"].append((list(messages), tools))
        reply = replies[state["i"]]
        state["i"] += 1
        return reply

    chat.state = state
    return chat


def test_agent_runs_a_tool_call_then_returns_final_text():
    chat = _fake_chat_script(
        {"role": "assistant", "content": None, "tool_calls": [
            {"id": "c1", "function": {
                "name": "centre_at_xyz",
                "arguments": json.dumps({"x": "1", "y": "2", "z": "3"})}}]},
        {"role": "assistant", "content": "Done - centred the view."},
    )
    out = agent.run_agent("centre at 1 2 3", chat=chat, top_k=None, verbose=False)
    assert out == "Done - centred the view."
    # The tool result must be fed back before the final turn.
    final_messages = chat.state["seen"][-1][0]
    tool_msgs = [m for m in final_messages if m.get("role") == "tool"]
    assert tool_msgs and tool_msgs[0]["content"] == "Centred at (1, 2, 3)"
    assert tool_msgs[0]["tool_call_id"] == "c1"


def test_agent_handles_reply_with_no_tool_calls():
    chat = _fake_chat_script({"role": "assistant", "content": "Hello!"})
    assert agent.run_agent("hi", chat=chat, top_k=None, verbose=False) == "Hello!"


def test_agent_stops_after_max_steps():
    loop_reply = {"role": "assistant", "content": None, "tool_calls": [
        {"id": "c", "function": {"name": "centre_at_xyz",
                                 "arguments": "{\"x\":\"0\",\"y\":\"0\",\"z\":\"0\"}"}}]}
    chat = _fake_chat_script(*([loop_reply] * 10))
    out = agent.run_agent("spin", chat=chat, max_steps=3, top_k=None, verbose=False)
    assert "Stopped after 3" in out
    assert chat.state["i"] == 3


# --- retrieval (fake embeddings) --------------------------------------------

def _bag_of_words_embed(vocab):
    """A fake embedder: each text -> a count vector over *vocab*.

    Deterministic and network-free, so retrieval ranking can be asserted:
    documents sharing more query words score higher under cosine.
    """
    def embed(texts):
        vectors = []
        for text in texts:
            words = text.lower().split()
            vectors.append([float(words.count(term)) for term in vocab])
        return vectors
    return embed


def test_retriever_ranks_by_similarity():
    from coot_commands.retrieval import ToolRetriever
    vocab = ["water", "add", "refine", "residue", "centre", "colour"]
    documents = {
        "add_water": "add a water molecule",
        "refine_residue": "refine a residue",
        "set_colour": "set the colour",
    }
    retriever = ToolRetriever(documents, _bag_of_words_embed(vocab))
    assert retriever.select("add a water please", k=1) == ["add_water"]
    assert retriever.select("refine this residue", k=1) == ["refine_residue"]


def test_retriever_k_limits_results_and_orders_them():
    from coot_commands.retrieval import ToolRetriever
    vocab = ["water", "refine", "colour"]
    documents = {"add_water": "water", "refine_residue": "refine",
                 "set_colour": "colour"}
    retriever = ToolRetriever(documents, _bag_of_words_embed(vocab))
    top = retriever.select("water refine colour", k=2)
    assert len(top) == 2


def test_command_documents_cover_every_command():
    from coot_commands.retrieval import command_documents
    docs = command_documents()
    assert len(docs) == len({c.name for c in all_commands()})
    assert all(text.strip() for text in docs.values())


def test_agent_uses_retrieved_subset_when_top_k_set():
    from coot_commands import retrieval
    captured = {}

    def chat(messages, tools):
        captured["tools"] = tools
        return {"role": "assistant", "content": "ok"}

    fake = retrieval.ToolRetriever(
        {"add_water": "add water", "refine_residue": "refine"},
        _bag_of_words_embed(["water", "refine", "add"]))
    orig = retrieval._default_retriever
    retrieval._default_retriever = fake
    try:
        agent.run_agent("add a water", chat=chat, top_k=1, verbose=False)
    finally:
        retrieval._default_retriever = orig
    names = [t["function"]["name"] for t in captured["tools"]]
    assert "add_water" in names                       # the retrieved command
    assert "get_active_residue" in names              # custom tools are pinned


def test_agent_falls_back_to_all_tools_when_retrieval_fails():
    from coot_commands import retrieval

    def boom(texts):
        raise RuntimeError("no embedding server")

    captured = {}

    def chat(messages, tools):
        captured["tools"] = tools
        return {"role": "assistant", "content": "ok"}

    fake = retrieval.ToolRetriever({"add_water": "add water"}, boom)
    orig = retrieval._default_retriever
    retrieval._default_retriever = fake
    try:
        agent.run_agent("do something", chat=chat, top_k=5, verbose=False)
    finally:
        retrieval._default_retriever = orig
    # Fallback exposes the full command set (plus the pinned custom tools).
    from coot_commands.tools import custom_tools
    assert len(captured["tools"]) == (
        len({c.name for c in all_commands()}) + len(custom_tools()))


def test_chat_url_normalisation():
    from coot_commands.agent import _normalise_chat_url as n
    full = "http://127.0.0.1:11435/v1/chat/completions"
    assert n("http://127.0.0.1:11435") == full           # bare base (the 405 case)
    assert n("http://127.0.0.1:11435/") == full          # trailing slash
    assert n("http://127.0.0.1:11435/v1") == full        # v1 root
    assert n(full) == full                               # already full: unchanged
    assert n(full + "/") == full                         # full with trailing slash


def test_embed_url_normalisation():
    from coot_commands.retrieval import _normalise_embed_url as n
    full = "http://127.0.0.1:11435/api/embed"
    assert n("http://127.0.0.1:11435") == full
    assert n("http://127.0.0.1:11435/") == full
    assert n("http://127.0.0.1:11435/api") == full
    assert n(full) == full
    assert n(full + "/") == full


# --- pluggable executor + events -------------------------------------------

def test_run_agent_uses_injected_executor_and_emits_events():
    calls = []

    def execute(name, args):
        calls.append((name, args))
        return "did " + name

    events = []
    chat = _fake_chat_script(
        {"role": "assistant", "content": None, "tool_calls": [
            {"id": "c1", "function": {
                "name": "add_water", "arguments": "{}"}}]},
        {"role": "assistant", "content": "Added a water."},
    )
    out = agent.run_agent("add water", chat=chat, execute=execute,
                          on_event=events.append, top_k=None, verbose=False)
    assert out == "Added a water."
    assert calls == [("add_water", {})]                      # our executor ran
    kinds = [e["type"] for e in events]
    assert "step" in kinds and kinds[-1] == "final"
    step = next(e for e in events if e["type"] == "step")
    assert step["tool"] == "add_water" and step["result"] == "did add_water"


# --- socket client (loopback fake server) -----------------------------------

def _recv_exactly(conn, n):
    buf = b""
    while len(buf) < n:
        chunk = conn.recv(n - len(buf))
        if not chunk:
            break
        buf += chunk
    return buf


def _fake_coot_server(responder):
    """A one-shot loopback server framing like json-rpc.cc; returns its port."""
    srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    srv.bind(("127.0.0.1", 0))
    srv.listen(1)
    port = srv.getsockname()[1]

    def run():
        conn, _ = srv.accept()
        (length,) = struct.unpack(">I", _recv_exactly(conn, 4))
        request = json.loads(_recv_exactly(conn, length).decode("utf-8"))
        payload = json.dumps(responder(request)).encode("utf-8")
        conn.sendall(struct.pack(">I", len(payload)) + payload)
        conn.close()
        srv.close()

    threading.Thread(target=run, daemon=True).start()
    return port


def test_socket_client_exec_python_returns_value():
    from coot_commands.socket_client import CootSocketClient
    seen = {}

    def responder(req):
        seen["req"] = req
        return {"jsonrpc": "2.0", "id": str(req["id"]),
                "result": {"value": "Centred on A/45"}}

    port = _fake_coot_server(responder)
    client = CootSocketClient(port=port)
    assert client.exec_python("1 + 1") == "Centred on A/45"
    assert seen["req"]["method"] == "python.exec"
    assert seen["req"]["params"]["code"] == "1 + 1"
    client.close()


def test_socket_client_raises_on_server_error():
    from coot_commands.socket_client import CootSocketClient, CootSocketError
    port = _fake_coot_server(
        lambda req: {"jsonrpc": "2.0", "id": str(req["id"]),
                     "error": {"code": -32001, "message": "boom"}})
    client = CootSocketClient(port=port)
    try:
        client.exec_python("bad")
        assert False, "expected CootSocketError"
    except CootSocketError as e:
        assert "boom" in str(e)
    client.close()


def test_socket_executor_builds_execute_tool_call():
    from coot_commands.socket_client import CootSocketClient, make_socket_executor
    seen = {}

    def responder(req):
        seen["code"] = req["params"]["code"]
        return {"jsonrpc": "2.0", "id": str(req["id"]),
                "result": {"value": "Centred on A/45 of model 0"}}

    port = _fake_coot_server(responder)
    execute = make_socket_executor(CootSocketClient(port=port))
    result = execute("go_to_residue", {"chain": "A", "resno": "45"})
    assert result == "Centred on A/45 of model 0"
    # The generated code invokes execute_tool with the name and JSON args.
    assert "execute_tool" in seen["code"]
    assert "go_to_residue" in seen["code"]
    assert '"chain": "A"' in seen["code"] or "'chain': 'A'" in seen["code"]


# --- agent_serve event loop -------------------------------------------------

def test_agent_serve_streams_events_per_request(monkeypatch=None):
    from coot_commands import agent_serve

    # Stub run_agent so the serve loop is tested without a model: it drives one
    # tool call through the injected executor and emits step/final events.
    def fake_run_agent(text, *, messages, execute, on_event, verbose):
        on_event({"type": "step", "tool": "load_tutorial", "args": {},
                  "result": execute("load_tutorial", {})})
        on_event({"type": "final", "text": "done: " + text})

    orig = agent_serve.run_agent
    agent_serve.run_agent = fake_run_agent

    class FakeClient:
        def exec_python(self, code):
            return "Loaded the tutorial model and data"

        def close(self):
            pass

    try:
        stdin = io.StringIO('{"text": "load tutorial"}\n')
        stdout = io.StringIO()
        agent_serve.serve(stdin, stdout, client=FakeClient(), startup_status=False)
    finally:
        agent_serve.run_agent = orig

    events = [json.loads(line) for line in stdout.getvalue().splitlines()]
    kinds = [e["type"] for e in events]
    assert kinds[0] == "ready"
    assert "step" in kinds
    assert any(e["type"] == "final" and e["text"] == "done: load tutorial"
               for e in events)
    assert any(e["type"] == "context" and "approx_tokens" in e for e in events)
    assert kinds[-1] == "done"


# --- custom (context) tools -------------------------------------------------

def test_custom_tools_are_registered_and_schematised():
    from coot_commands.tools import custom_tools
    names = [t["function"]["name"] for t in custom_tools()]
    assert "get_active_residue" in names
    for tool in custom_tools():
        assert tool["type"] == "function"
        assert tool["function"]["description"]


def test_execute_tool_dispatches_to_custom_tool():
    from coot_commands import tools

    @tools.custom_tool("unit_probe", "test probe",
                       parameters={"type": "object",
                                   "properties": {"x": {"type": "string"}}})
    def _probe(x=None):
        return f"probe:{x}"

    try:
        assert tools.execute_tool("unit_probe", {"x": "7"}) == "probe:7"
        # Unknown/extra args are filtered out, not passed through as a TypeError.
        assert tools.execute_tool("unit_probe", {"x": "7", "bogus": "9"}) == "probe:7"
    finally:
        tools._CUSTOM_TOOLS.pop("unit_probe", None)


def test_get_active_residue_without_coot():
    # Standalone (no coot) it reports the API is unavailable rather than raising.
    from coot_commands.tools import execute_tool
    out = execute_tool("get_active_residue", {})
    assert "Coot API is not available" in out


def test_agent_pins_custom_tools_even_with_no_commands():
    captured = {}

    def chat(messages, tools):
        captured["tools"] = tools
        return {"role": "assistant", "content": "ok"}

    from coot_commands.tools import custom_tools
    agent.run_agent("do nothing", chat=chat, tools=None, top_k=None, verbose=False)
    names = [t["function"]["name"] for t in captured["tools"]]
    for custom in custom_tools():
        assert custom["function"]["name"] in names


def test_run_agent_threads_conversation_across_calls():
    conversation = []
    chat1 = _fake_chat_script({"role": "assistant", "content": "refined A 42"})
    agent.run_agent("refine the worst residue", chat=chat1,
                    messages=conversation, top_k=None, verbose=False)
    # The running conversation retains system + this turn.
    assert [m["role"] for m in conversation] == ["system", "user", "assistant"]

    chat2 = _fake_chat_script({"role": "assistant", "content": "ok"})
    agent.run_agent("focus on the residue you just refined", chat=chat2,
                    messages=conversation, top_k=None, verbose=False)
    # The second call sees the first turn's history (that's the memory).
    seen = chat2.state["seen"][0][0]
    contents = [m.get("content") for m in seen]
    assert "refine the worst residue" in contents
    assert "refined A 42" in contents
    assert "focus on the residue you just refined" in contents


def test_agent_serve_threads_conversation_and_reset():
    from coot_commands import agent_serve
    snapshots = []

    def fake_run_agent(text, *, messages, execute, on_event, verbose):
        snapshots.append(list(messages))          # history coming into this turn
        messages.append({"role": "user", "content": text})
        messages.append({"role": "assistant", "content": "ok:" + text})
        on_event({"type": "final", "text": "ok:" + text})

    class FakeClient:
        def exec_python(self, code):
            return "x"

        def close(self):
            pass

    orig = agent_serve.run_agent
    agent_serve.run_agent = fake_run_agent
    try:
        stdin = io.StringIO('{"text": "first"}\n{"text": "second"}\n'
                            '{"reset": true}\n{"text": "third"}\n')
        stdout = io.StringIO()
        agent_serve.serve(stdin, stdout, client=FakeClient(), startup_status=False)
    finally:
        agent_serve.run_agent = orig

    assert snapshots[0] == []                                   # first: no history
    assert any(m.get("content") == "first" for m in snapshots[1])  # second sees first
    assert snapshots[2] == []                                   # after reset: cleared
    events = [json.loads(line) for line in stdout.getvalue().splitlines()]
    assert any(e["type"] == "reset" for e in events)


def test_agent_serve_emits_startup_status():
    from coot_commands import agent_serve

    class FakeClient:
        def connect(self):
            pass          # RPC reachable

        def exec_python(self, code):
            return "x"

        def close(self):
            pass

    orig_probe = agent_serve._probe_ollama
    agent_serve._probe_ollama = lambda timeout=2.0: (True, "")
    try:
        stdout = io.StringIO()
        agent_serve.serve(io.StringIO(""), stdout, client=FakeClient(),
                          startup_status=True)
    finally:
        agent_serve._probe_ollama = orig_probe

    events = [json.loads(line) for line in stdout.getvalue().splitlines()]
    status = [e for e in events if e["type"] == "status"]
    assert status, "expected a status event"
    assert status[0]["rpc"] is True
    assert status[0]["ollama"] is True
    assert "model" in status[0]


def test_agent_serve_startup_status_reports_rpc_failure():
    from coot_commands import agent_serve

    class DeadClient:
        def connect(self):
            raise RuntimeError("connection refused")

        def close(self):
            pass

    orig_probe = agent_serve._probe_ollama
    agent_serve._probe_ollama = lambda timeout=2.0: (False, "no server")
    try:
        stdout = io.StringIO()
        agent_serve.serve(io.StringIO(""), stdout, client=DeadClient(),
                          startup_status=True)
    finally:
        agent_serve._probe_ollama = orig_probe

    status = [json.loads(l) for l in stdout.getvalue().splitlines()
              if json.loads(l)["type"] == "status"][0]
    assert status["rpc"] is False and "refused" in status["rpc_detail"]
    assert status["ollama"] is False


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
