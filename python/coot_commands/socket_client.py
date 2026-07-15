# coot_commands/socket_client.py
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

"""Talk to a running Coot from another process, over its JSON-RPC socket.

The agent runs as a separate process (so its slow model calls never block
Coot's GUI), but its tool calls must reach the live Coot.  Coot already
serves a length-prefixed JSON-RPC socket on localhost (see ``src/json-rpc.cc``,
started via ``make_socket_listener_maybe``); this client speaks that wire
format and runs a ``python.exec`` there.  Because the server processes the
request on Coot's GTK idle function, the executed code runs on the main
thread - exactly where the Coot API is safe to call.

The frame format is a 4-byte big-endian length prefix followed by the JSON
payload, in both directions.  :meth:`CootSocketClient.exec_python` sends one
request and reads one response (a single client, used sequentially).

:func:`make_socket_executor` adapts the client into the ``execute`` callback
that :func:`coot_commands.agent.run_agent` expects: it runs a command by
calling :func:`coot_commands.tools.execute_tool` *inside* Coot, so the same
registry and argument handling apply whether a command runs in-process or
over the socket.
"""

from __future__ import annotations

import json
import os
import socket
import struct
import time
from typing import Any, Dict, Optional

DEFAULT_HOST = "127.0.0.1"
# Coot's default remote-control port (graphics_info_t::remote_control_port_number;
# vte.cc falls back to 9090 when it is unset).
DEFAULT_PORT = 9090


class CootSocketError(RuntimeError):
    """A transport-level failure talking to Coot (connection or protocol)."""


class CootSocketClient:
    """A client for Coot's length-prefixed JSON-RPC socket."""

    def __init__(self, host: str = DEFAULT_HOST, port: Optional[int] = None,
                 timeout: float = 30.0) -> None:
        self.host = host
        self.port = port if port is not None else int(
            os.environ.get("COOT_RPC_PORT", DEFAULT_PORT))
        self.timeout = timeout
        self._sock: Optional[socket.socket] = None
        self._next_id = 1

    def connect(self, retries: int = 15, delay: float = 0.2) -> None:
        """Connect to Coot, retrying briefly so a startup race can't fail us.

        Coot brings the listener up and spawns this process at nearly the same
        moment, so the first connect can land a hair too early. We retry for
        ~*retries* x *delay* seconds before giving up with the last error.
        """
        if self._sock is not None:
            return
        last_error: Optional[OSError] = None
        for attempt in range(max(1, retries)):
            try:
                self._sock = socket.create_connection(
                    (self.host, self.port), timeout=self.timeout)
                return
            except OSError as e:
                last_error = e
                if attempt < retries - 1:
                    time.sleep(delay)
        raise CootSocketError(
            f"cannot connect to Coot at {self.host}:{self.port} after "
            f"{retries} attempts: {last_error}") from None

    def close(self) -> None:
        if self._sock is not None:
            try:
                self._sock.close()
            finally:
                self._sock = None

    def _send_frame(self, payload: bytes) -> None:
        assert self._sock is not None
        self._sock.sendall(struct.pack(">I", len(payload)) + payload)

    def _recv_exactly(self, n: int) -> bytes:
        assert self._sock is not None
        chunks = []
        remaining = n
        while remaining > 0:
            chunk = self._sock.recv(remaining)
            if not chunk:
                raise CootSocketError("Coot closed the connection")
            chunks.append(chunk)
            remaining -= len(chunk)
        return b"".join(chunks)

    def _recv_frame(self) -> bytes:
        (length,) = struct.unpack(">I", self._recv_exactly(4))
        return self._recv_exactly(length)

    def exec_python(self, code: str) -> str:
        """Evaluate *code* (a single expression) in Coot; return its value string.

        Raises :class:`CootSocketError` on a transport failure or if the server
        reports an error.
        """
        self.connect()
        request_id = self._next_id
        self._next_id += 1
        request = {
            "jsonrpc": "2.0",
            "id": request_id,
            "method": "python.exec",
            "params": {"code": code},
        }
        self._send_frame(json.dumps(request).encode("utf-8"))
        response = json.loads(self._recv_frame().decode("utf-8"))
        if "error" in response:
            message = response["error"].get("message", "unknown error")
            raise CootSocketError(f"Coot error: {message}")
        result = response.get("result") or {}
        return result.get("value", "")


def make_socket_executor(client: CootSocketClient):
    """Return an ``execute(name, args)`` that runs a command inside Coot.

    The command runs via :func:`coot_commands.tools.execute_tool` on the Coot
    side, as a single ``__import__`` expression so no separate import statement
    is needed, mirroring how the Command tab evaluates its Python.
    """
    def execute(name: str, args: Dict[str, Any]) -> str:
        args_json = json.dumps(args)
        code = (
            "__import__('coot_commands.tools', fromlist=['execute_tool'])"
            ".execute_tool({name!r}, __import__('json').loads({args!r}))"
        ).format(name=name, args=args_json)
        return client.exec_python(code)
    return execute
