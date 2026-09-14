# coot_commands/tools.py
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

"""Expose the command registry to a tool-calling language model.

This is the bridge that lets a small local model (Gemma, Qwen, ...) drive
Coot: it turns each :func:`~coot_commands.registry.command` into an
OpenAI-style *tool* (a JSON schema of name, description and parameters)
and runs a tool call the model emits back through the registered handler.

Nothing here talks to a model or a network - see :mod:`coot_commands.agent`
for the loop that does.  Keeping the bridge separate means it is pure
Python and unit-testable without Coot or a model server: the command
modules ``try: import coot / except ImportError: coot = None``, so the
schemas build and (side-effect-free) handlers run standalone.

Why not hand the model Coot's ~400 raw API functions?  A 7-8B model
degrades badly past a few dozen tools.  The curated ``@command`` set - each
with a one-line description and example phrasings - is a far better tool
surface, and :func:`command_tools` accepts a *names* subset so a future
retrieval step can narrow it further per request.

The handler signature is the source of truth for a command's parameters:
by convention it mirrors the regex's named groups (e.g.
``go_to_residue(chain, resno, model=None)``), so :func:`inspect.signature`
yields both the parameter list and which are required (no default) versus
optional (default ``None``, resolved to the active molecule).
"""

from __future__ import annotations

import inspect
from typing import Any, Dict, Iterable, List, Optional

from coot_commands.registry import Command, all_commands
from coot_commands.types import ArgType, CommandError

# JSON-schema parameter descriptions per argument kind.  The value the model
# supplies is always coerced to a string before the handler sees it (handlers
# expect the same strings the regex would capture), so every parameter is typed
# "string"; the ArgType only enriches the human-readable description.
_ARG_DESCRIPTIONS = {
    ArgType.MODEL: "Model (molecule) number, e.g. \"0\".",
    ArgType.MAP: "Map (molecule) number, e.g. \"1\".",
    ArgType.COLOUR: "A colour name, e.g. \"red\" or \"sky blue\".",
}

# What to add for an optional argument of a given kind: the fallback the shared
# resolvers apply when the model omits it (see coot_commands.types).
_ARG_OMIT_HINTS = {
    ArgType.MODEL: " Omit to act on the active model.",
    ArgType.MAP: " Omit to use the map set for refinement.",
}


def _handler_params(cmd: Command) -> List[inspect.Parameter]:
    """The command's real arguments: named handler params, minus ``**kwargs``.

    Some handlers accept ``**_`` (they take no arguments but must swallow the
    named groups ``dispatch`` would pass); those contribute no tool parameters.
    """
    sig = inspect.signature(cmd.handler)
    return [p for p in sig.parameters.values()
            if p.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD,
                          inspect.Parameter.KEYWORD_ONLY)]


def _param_schema(cmd: Command, param: inspect.Parameter) -> Dict[str, Any]:
    """JSON schema for one parameter, described using its :class:`ArgType`."""
    arg_type = cmd.arg_types.get(param.name)
    required = param.default is inspect.Parameter.empty
    description = _ARG_DESCRIPTIONS.get(arg_type, "")
    if not required:
        description += _ARG_OMIT_HINTS.get(arg_type, "")
    schema: Dict[str, Any] = {"type": "string"}
    if description:
        schema["description"] = description.strip()
    return schema


def _description(cmd: Command) -> str:
    """The tool description: the command's help plus its example phrasings.

    Example phrasings matter a lot for a small model - they show the natural
    language that maps to this command, so the model can pick the right tool
    from an ambiguous request.
    """
    text = (cmd.help_text or cmd.description or cmd.name).strip()
    if cmd.examples:
        text += "\n\nExample phrasings: " + "; ".join(
            f'"{ex}"' for ex in cmd.examples[:4])
    return text


def command_to_tool(cmd: Command) -> Dict[str, Any]:
    """Render one :class:`Command` as an OpenAI-style tool definition."""
    properties: Dict[str, Any] = {}
    required: List[str] = []
    for param in _handler_params(cmd):
        properties[param.name] = _param_schema(cmd, param)
        if param.default is inspect.Parameter.empty:
            required.append(param.name)
    parameters: Dict[str, Any] = {"type": "object", "properties": properties}
    if required:
        parameters["required"] = required
    return {
        "type": "function",
        "function": {
            "name": cmd.name,
            "description": _description(cmd),
            "parameters": parameters,
        },
    }


def _commands_by_name() -> Dict[str, Command]:
    """Map tool name -> command, keeping the first on a name clash.

    Tool names must be unique; two commands sharing a handler ``__name__``
    (different modules, same function name) would otherwise collide, so we keep
    the first-registered and skip the rest.
    """
    by_name: Dict[str, Command] = {}
    for cmd in all_commands():
        by_name.setdefault(cmd.name, cmd)
    return by_name


def command_tools(names: Optional[Iterable[str]] = None) -> List[Dict[str, Any]]:
    """Return tool definitions for the registered commands.

    Pass *names* to expose only a subset (e.g. the output of a retrieval step
    that picked the commands relevant to a request); by default every command
    is exposed.
    """
    by_name = _commands_by_name()
    selected = list(by_name) if names is None else [n for n in names if n in by_name]
    return [command_to_tool(by_name[n]) for n in selected]


# Custom (context/query) tools: agent-only tools that are NOT @command regex
# handlers. A command is an ACTION triggered by typed or spoken language; a
# custom tool is typically a QUERY returning live context - e.g. the residue at
# the centre of the screen - so the model can resolve deictic references like
# "here", "this residue" or "the current position" before it acts. Register one
# with @custom_tool (see coot_commands.context_tools). They are always exposed
# to the model (never dropped by retrieval), since such context is relevant
# regardless of how a request is worded.
_CUSTOM_TOOLS: Dict[str, Dict[str, Any]] = {}


def custom_tool(name: str, description: str,
                parameters: Optional[Dict[str, Any]] = None) -> Callable:
    """Register *handler* as an agent tool named *name*.

    *parameters* is a JSON-schema object for the arguments (default: none).  The
    handler returns a result string (like a command handler); it receives the
    model-supplied arguments as keyword arguments.
    """
    schema_params = parameters or {"type": "object", "properties": {}}

    def decorator(handler: Callable[..., str]) -> Callable[..., str]:
        _CUSTOM_TOOLS[name] = {
            "schema": {
                "type": "function",
                "function": {
                    "name": name,
                    "description": description,
                    "parameters": schema_params,
                },
            },
            "handler": handler,
        }
        return handler

    return decorator


def custom_tools() -> List[Dict[str, Any]]:
    """Tool definitions for the always-available custom (context) tools."""
    return [entry["schema"] for entry in _CUSTOM_TOOLS.values()]


def _custom_kwargs(handler: Callable, arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Filter *arguments* to those the custom *handler* actually accepts."""
    sig = inspect.signature(handler)
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
        return dict(arguments)
    allowed = {name for name, p in sig.parameters.items()
               if p.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD,
                             inspect.Parameter.KEYWORD_ONLY)}
    return {k: v for k, v in arguments.items() if k in allowed}


def execute_tool(name: str, arguments: Optional[Dict[str, Any]] = None) -> str:
    """Run the tool *name* (custom tool or command) with *arguments*.

    Returns the handler's result string, or a readable error string (never
    raises) so the agent loop can feed failures straight back to the model.

    Command arguments are coerced to strings and passed by the handler's
    parameter name, mirroring exactly what ``registry.dispatch`` passes from a
    regex ``groupdict`` - a missing argument arrives as ``None`` and the shared
    resolvers fall back to the active molecule (or raise a clear
    :class:`CommandError`).
    """
    arguments = arguments or {}

    custom = _CUSTOM_TOOLS.get(name)
    if custom is not None:
        try:
            return custom["handler"](**_custom_kwargs(custom["handler"], arguments))
        except CommandError as e:
            return f"Error: {e}"
        except Exception as e:  # noqa: BLE001 - report any failure to the model
            return f"Error running '{name}': {e}"

    cmd = _commands_by_name().get(name)
    if cmd is None:
        return f"Error: unknown command '{name}'"
    kwargs = {}
    for param in _handler_params(cmd):
        value = arguments.get(param.name)
        kwargs[param.name] = None if value is None else str(value)
    try:
        return cmd.handler(**kwargs)
    except CommandError as e:
        return f"Error: {e}"
    except Exception as e:  # noqa: BLE001 - report any handler failure to the model
        return f"Error running '{name}': {e}"
