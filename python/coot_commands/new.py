# coot_commands/new.py
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

"""Scaffold a new command from a pattern and examples.

Writing the boilerplate by hand is the tedious part: the licence header, the
``import`` lines, and a handler signature whose parameters exactly match the
pattern's named groups (with the optional ones defaulting to ``None``).  This
helper does that for you.

It derives the signature *from the examples*: it matches each example against
the pattern, and a named group that is ``None`` in any example is treated as
optional (``= None``), the rest as required.  So give at least one example that
omits each optional argument and the signature comes out right.

Run it interactively::

    python3 -m coot_commands.new

It prints a ready-to-paste ``@command`` block.  If you name a command *file*
that does not exist yet, it offers to create it (header + imports + the stub)
and reminds you to add it to ``python/Makefile.am``.  It never edits an
existing file - appending blindly would put a specific pattern *after* the
general ones and get it shadowed (see ordering in ``doc/writing-commands.md``),
so for an existing file it prints the block for you to place by hand.
"""

from __future__ import annotations

import os
import re
import sys
from typing import List

from coot_commands.registry import normalise

_GROUP_RE = re.compile(r"\(\?P<([A-Za-z_]\w*)>")

_HEADER = '''\
# coot_commands/commands/{module}.py
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

"""{summary}"""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import resolve_model, resolve_map, as_int, as_float, CommandError

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "{category}"
'''


def group_names(pattern: str) -> List[str]:
    """The named groups in *pattern*, in order of appearance."""
    seen: List[str] = []
    for name in _GROUP_RE.findall(pattern):
        if name not in seen:
            seen.append(name)
    return seen


def signature(pattern: str, examples: List[str]) -> str:
    """Build the handler parameter list for *pattern* from its *examples*.

    A named group is optional (``name=None``) if it fails to participate in
    at least one example match; otherwise it is required.  Groups that never
    appear in any successful match are assumed optional, to be safe.
    """
    names = group_names(pattern)
    if not names:
        return ""
    regex = re.compile(pattern, re.IGNORECASE)
    optional = {name: False for name in names}
    matched_any = False
    for example in examples:
        m = regex.match(normalise(example))
        if not m:
            continue
        matched_any = True
        groups = m.groupdict()
        for name in names:
            if groups.get(name) is None:
                optional[name] = True
    if not matched_any:
        optional = {name: True for name in names}
    required = [n for n in names if not optional[n]]
    opt = [n for n in names if optional[n]]
    parts = [f"{n}: str" for n in required]
    parts += [f"{n}: Optional[str] = None" for n in opt]
    return ", ".join(parts)


def stub(name: str, pattern: str, examples: List[str],
         help_text: str, notes: str) -> str:
    """Render the ``@command`` decorator + handler stub as source text."""
    ex_lines = ", ".join(repr(e) for e in examples) or repr(pattern)
    deco = [f'@command(r"{pattern}",',
            f'         examples=[{ex_lines}],',
            f'         category=CATEGORY,']
    if notes:
        deco.append(f'         notes={notes!r},')
    # Drop trailing comma on the last kwarg, close the call.
    deco[-1] = deco[-1].rstrip(",") + ")"

    sig = signature(pattern, examples)
    doc = help_text or "TODO: one-line description."
    body = ["    \"\"\"" + doc + "\"\"\"",
            "    # TODO: resolve arguments and call the coot.* API, then return a",
            "    # short status string. Coerce captured strings via types.py helpers.",
            "    raise CommandError(\"not implemented yet\")"]
    return "\n".join(deco + [f"def {name}({sig}) -> str:"] + body) + "\n"


def _prompt(label: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    try:
        value = input(f"{label}{suffix}: ").strip()
    except EOFError:
        value = ""
    return value or default


def _prompt_examples() -> List[str]:
    print("Examples (one per line, blank to finish; the first is canonical):")
    examples: List[str] = []
    while True:
        try:
            line = input("  example> ").strip()
        except EOFError:
            break
        if not line:
            break
        examples.append(line)
    return examples


def interactive() -> int:
    print("Scaffold a new Coot command. Ctrl-C to abort.\n")
    name = _prompt("Handler function name (e.g. refine_chain)")
    if not name.isidentifier():
        sys.stderr.write(f"error: {name!r} is not a valid function name\n")
        return 1
    category = _prompt("Category", "General")
    pattern = _prompt("Pattern (regex; named groups become arguments)")
    if not pattern:
        sys.stderr.write("error: a pattern is required\n")
        return 1
    try:
        re.compile(pattern)
    except re.error as exc:
        sys.stderr.write(f"error: pattern is not a valid regex: {exc}\n")
        return 1
    examples = _prompt_examples()
    help_text = _prompt("One-line help")
    notes = _prompt("Notes (optional, longer prose for the docs)")

    # Validate the examples up front - the same check the test suite enforces.
    regex = re.compile(pattern, re.IGNORECASE)
    bad = [e for e in examples if not regex.match(normalise(e))]
    if bad:
        sys.stderr.write("\nwarning: these examples do NOT match the pattern "
                         "(fix the pattern or the example):\n")
        for e in bad:
            sys.stderr.write(f"  {e!r}\n")

    block = stub(name, pattern, examples, help_text, notes)

    module = _prompt("\nTarget command module (file stem under commands/, "
                     "e.g. refine)")
    print()
    if not module:
        print("# Paste this into a file in coot_commands/commands/:\n")
        print(block)
        return 0

    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "commands", f"{module}.py")
    if os.path.exists(path):
        print(f"# {module}.py already exists - not editing it (ordering matters:")
        print("# specific patterns must come before general ones). Paste this in,")
        print("# placing it above any more-general pattern:\n")
        print(block)
        return 0

    summary = f"Commands for {category.lower()}."
    contents = _HEADER.format(module=module, summary=summary,
                              category=category) + "\n\n" + block
    with open(path, "w") as fh:
        fh.write(contents)
    print(f"wrote {path}")
    print("\nNext:")
    print(f"  1. add 'coot_commands/commands/{module}.py' to "
          "python/Makefile.am (nobase_dist_pkgpython_PYTHON)")
    print("  2. implement the handler body (it raises NotImplemented for now)")
    print(f"  3. dry-run it:  python3 -m coot_commands.try "
          f"{examples[0]!r}" if examples else
          "  3. dry-run it with:  python3 -m coot_commands.try '<your input>'")
    return 0


def main() -> int:
    try:
        return interactive()
    except KeyboardInterrupt:
        sys.stderr.write("\naborted\n")
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
