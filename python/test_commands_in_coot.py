# test_commands_in_coot.py
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

"""Integration tests for the command interface, run *inside* Coot.

The companion ``test_coot_commands.py`` stubs out ``coot`` and checks parsing,
routing and completion without a running program.  This module is the other
half: it drives the *real* command interface (``coot_command_interface.
run_command``) against a live Coot session with the tutorial loaded, and
verifies the observable effect of each command through the ``coot`` API - a
map really is contoured, a residue really is deleted, the view really recentres.

Run it from inside Coot - at the Python prompt (Calculate > Scripting > Python,
or the Command tab in Python mode)::

    import test_commands_in_coot as t
    t.run_all()

or as a start-up script::

    coot --script test_commands_in_coot.py          # with graphics
    coot --no-graphics --script test_commands_in_coot.py

It prints ``PASS``/``FAIL`` per check and a summary.  When run as a script it
sets the process exit code (0 = everything passed) so it can gate CI.

It is destructive to the *session* (it loads the tutorial, edits it, and closes
every molecule at the end) but touches nothing on disk.  Network fetches,
dictionary look-ups and whole-model merges are skipped by default (see
``_SWEEP_SKIP``) so the run stays fast and offline.
"""

from __future__ import annotations

import sys

try:
    import coot
except ImportError:
    coot = None

import coot_command_interface as cli
from coot_commands.registry import all_commands, normalise
from coot_commands.types import loaded_models, loaded_maps


# ---------------------------------------------------------------------------
#   Tiny test harness
# ---------------------------------------------------------------------------

_passes = 0
_failures: list[str] = []


def check(cond: bool, label: str, detail: str = "") -> bool:
    """Record one assertion; print PASS/FAIL and return the condition."""
    global _passes
    if cond:
        _passes += 1
        print(f"  PASS {label}")
    else:
        msg = f"{label}" + (f" - {detail}" if detail else "")
        _failures.append(msg)
        print(f"  FAIL {msg}")
    return bool(cond)


def run(text: str) -> str:
    """Run a command line through the real interface, echoing the result."""
    out = cli.run_command(text)
    first = out.splitlines()[0] if out else ""
    print(f"    $ {text}\n      -> {first}")
    return out


def _looks_like_error(msg: str) -> bool:
    """True if a ``run_command`` result is a dispatch/execution failure.

    ``run_command`` never raises; it turns an unmatched command or a handler
    exception into a message with one of these prefixes.  A handler that ran
    fine but reports "no such residue" is *not* an error here - the command
    dispatched and executed, which is what the smoke sweep checks.
    """
    return msg.startswith(("Error running", "Unrecognised command", "Error:"))


def ran_ok(text: str) -> bool:
    """Run *text* and return True if it dispatched and executed without error."""
    return not _looks_like_error(run(text))


def approx(a: float, b: float, tol: float = 1e-3) -> bool:
    return abs(a - b) <= tol


# ---------------------------------------------------------------------------
#   Session fixture: a loaded tutorial with known model / map / residues
# ---------------------------------------------------------------------------

class Session:
    """The loaded tutorial, with a chain and some valid residue numbers."""

    def __init__(self) -> None:
        self.model = int(loaded_models()[0])
        maps = [int(m) for m in loaded_maps()]
        # The tutorial loads the 2Fo-Fc map first, then the Fo-Fc difference map.
        self.map = maps[0]
        self.diff_map = maps[-1]
        self.chain = coot.chain_id_py(self.model, 0)
        n = coot.chain_n_residues(self.chain, self.model)
        # Work from a quarter of the way in, away from the (possibly ragged)
        # termini, and grab a run of consecutive residues to act on.
        base = max(1, n // 4)
        self.resnos = [coot.seqnum_from_serial_number(self.model, self.chain,
                                                      base + k)
                       for k in range(12)]
        self.last_resno = coot.seqnum_from_serial_number(self.model, self.chain,
                                                         n - 1)

    def res(self, k: int) -> int:
        return self.resnos[k]


def setup_session() -> "Session":
    """Close everything, load the tutorial, and set the refinement map."""
    run("delete all confirm")
    run("load tutorial")
    if not loaded_models() or not loaded_maps():
        raise RuntimeError("tutorial did not load - cannot run integration tests")
    sess = Session()
    coot.set_imol_refinement_map(sess.map)
    # Centre somewhere definite so there is an active residue.
    run(f"go to {sess.chain} {sess.res(0)}")
    return sess


# ---------------------------------------------------------------------------
#   Assertive tests, grouped by category
# ---------------------------------------------------------------------------

def test_display(s: Session) -> None:
    print("[Display]")
    run(f"hide model {s.model}")
    check(coot.mol_is_displayed(s.model) == 0, "hide model", "still displayed")
    run(f"show model {s.model}")
    check(coot.mol_is_displayed(s.model) == 1, "show model", "not displayed")
    run(f"hide map {s.map}")
    check(coot.map_is_displayed(s.map) == 0, "hide map", "still displayed")
    run(f"show map {s.map}")
    check(coot.map_is_displayed(s.map) == 1, "show map", "not displayed")
    check(ran_ok("show environment"), "show environment")
    check(ran_ok("hide environment"), "hide environment")


def test_maps(s: Session) -> None:
    print("[Maps]")
    run(f"contour map {s.map} to 0.42")
    check(approx(coot.get_contour_level_absolute(s.map), 0.42, 1e-2),
          "contour absolute", f"level={coot.get_contour_level_absolute(s.map)}")
    run(f"contour map {s.map} to 1.3 sigma")
    check(approx(coot.get_contour_level_in_sigma(s.map), 1.3, 1e-2),
          "contour sigma", f"sigma={coot.get_contour_level_in_sigma(s.map)}")
    check(ran_ok(f"colour map {s.map} blue"), "colour map")
    run(f"map {s.diff_map} is a difference map")
    check(coot.map_is_difference_map(s.diff_map) == 1, "difference map")


def test_view(s: Session) -> None:
    print("[View]")
    check(ran_ok("background black"), "background")
    check(ran_ok("set zoom to 30"), "zoom")
    check(ran_ok("orthographic"), "orthographic")
    check(ran_ok("perspective view"), "perspective")
    # Toggles: flip on then back off so we leave the session as we found it.
    check(ran_ok("spin") and ran_ok("spin"), "spin toggle")
    check(ran_ok("rock") and ran_ok("rock"), "rock toggle")


def test_navigation(s: Session) -> None:
    print("[Navigation]")
    run(f"go to {s.chain} {s.res(2)}")
    active = coot.active_residue_py()
    check(bool(active) and active[1] == s.chain and active[2] == s.res(2),
          "go to residue", f"active={active[:3] if active else None}")
    run("centre at 10 20 30")
    check(approx(coot.rotation_centre_position(0), 10.0) and
          approx(coot.rotation_centre_position(1), 20.0) and
          approx(coot.rotation_centre_position(2), 30.0), "centre at xyz")
    run(f"go to {s.chain} {s.res(2)}")
    check(ran_ok("next residue"), "next residue")
    check(ran_ok("previous residue"), "previous residue")


def test_refinement(s: Session) -> None:
    print("[Refinement]")
    coot.set_imol_refinement_map(s.map)
    check(ran_ok(f"refine {s.chain} {s.res(0)} to {s.res(2)}"), "refine range")
    check(ran_ok(f"refine {s.chain} {s.res(3)}"), "refine residue")
    check(ran_ok(f"backrub {s.chain} {s.res(4)}"), "backrub")
    check(ran_ok(f"refine sphere {s.chain} {s.res(5)}"), "refine sphere")
    # Bare form must fall back to the active residue, not error out.
    run(f"go to {s.chain} {s.res(3)}")
    check(ran_ok("backrub rotamer"), "backrub (active residue)")


def test_model_editing(s: Session) -> None:
    print("[Model editing]")
    coot.set_imol_refinement_map(s.map)
    check(ran_ok(f"autofit {s.chain} {s.res(2)}"), "autofit")
    r_mut = s.res(6)
    run(f"mutate residue {s.chain} {r_mut} to ALA")
    check(coot.does_residue_exist_p(s.model, s.chain, r_mut, "") == 1,
          "mutate (residue kept)")
    check(ran_ok(f"add altconf {s.chain} {s.res(7)}"), "add altconf")
    check(ran_ok(f"pepflip {s.chain} {s.res(1)}"), "pepflip")
    check(ran_ok(f"add OXT to {s.chain}/{s.last_resno}"), "add OXT")
    r_del = s.res(9)
    exists_before = coot.does_residue_exist_p(s.model, s.chain, r_del, "") == 1
    run(f"delete residue {s.chain} {r_del}")
    gone = coot.does_residue_exist_p(s.model, s.chain, r_del, "") == 0
    check(exists_before and gone, "delete residue",
          f"before={exists_before} gone={gone}")


def test_building(s: Session) -> None:
    print("[Building]")
    check(ran_ok("add water"), "add water")


def test_validation(s: Session) -> None:
    print("[Validation]")
    check(ran_ok("validate anomalies"), "validate anomalies")
    out = run("find blobs")
    check(not _looks_like_error(out), "find blobs")
    if "blob(s)" in out or "blob 1" in out:
        check(ran_ok("go to blob 1"), "go to blob")


def test_representation(s: Session) -> None:
    print("[Representation]")
    for cmd in ("show ribbons", "hide ribbons", "show surface", "hide surface",
                "grey carbons", "coloured carbons",
                "show symmetry", "hide symmetry"):
        check(ran_ok(cmd), cmd)


def test_labels_state_session_help(s: Session) -> None:
    print("[Labels / State / Session / Help]")
    check(ran_ok("clear labels"), "clear labels")
    check(ran_ok("undo"), "undo")
    check("model" in run("list models"), "list models")
    check("map" in run("list maps"), "list maps")
    check("Available commands" in run("help"), "help")


# ---------------------------------------------------------------------------
#   Generic smoke sweep: every registered command's first example dispatches
# ---------------------------------------------------------------------------

# Handlers skipped by the automatic sweep: network fetches, dictionary
# look-ups, whole-model ops needing a second model, UI-disrupting toggles, and
# the destructive close-everything commands (covered explicitly at the end).
_SWEEP_SKIP = {
    "fetch_pdb", "fetch_pdb_and_map", "fetch_pdb_redo", "fetch_monomer",
    "add_monomer", "fit_ligand",
    "merge_models", "superpose_models",
    "toggle_fullscreen",
    "go_to_blob",  # needs a prior "find blobs"; covered in test_validation
    "delete_map", "delete_model", "delete_molecule",
    "delete_all", "delete_all_confirm",
}


def test_command_sweep(s: Session) -> None:
    """Run every non-skipped command's canonical example; expect a clean run.

    This keeps coverage complete automatically: a newly added command is
    exercised here the moment it is registered, without editing this file.
    """
    print("[Sweep: every command's first example]")
    for cmd in all_commands():
        if cmd.name in _SWEEP_SKIP or not cmd.examples:
            continue
        example = cmd.examples[0]
        out = cli.run_command(example)
        ok = not _looks_like_error(out)
        first = out.splitlines()[0] if out else "(empty)"
        check(ok, f"{cmd.name}: {example!r}", first)


# ---------------------------------------------------------------------------
#   Destructive teardown tests: run last, they close molecules
# ---------------------------------------------------------------------------

def test_delete(s: Session) -> None:
    print("[Delete]  (tears the session down)")
    run(f"delete map {s.diff_map}")
    check(coot.is_valid_map_molecule(s.diff_map) == 0, "delete map")
    run(f"delete model {s.model}")
    check(coot.is_valid_model_molecule(s.model) == 0, "delete model")
    run("delete all confirm")
    remaining = [i for i in range(coot.graphics_n_molecules())
                 if coot.is_valid_model_molecule(i) or coot.is_valid_map_molecule(i)]
    check(remaining == [], "delete all", f"left: {remaining}")


# ---------------------------------------------------------------------------
#   Runner
# ---------------------------------------------------------------------------

def run_all() -> bool:
    """Run the whole suite; return True if every check passed."""
    global _passes, _failures
    _passes, _failures = 0, []

    if coot is None:
        print("This suite must run inside Coot (import coot failed).")
        return False

    # Sanity-check the registry matches its own advertised examples first; if
    # this fails the sweep below is meaningless.
    from coot_commands.registry import unmatched_examples
    bad = unmatched_examples()
    check(bad == [], "all examples match their own command", str(bad))

    session = setup_session()

    # Non-destructive assertive tests.
    for test in (test_display, test_maps, test_view, test_navigation,
                 test_refinement, test_model_editing, test_building,
                 test_validation, test_representation,
                 test_labels_state_session_help):
        try:
            test(session)
        except Exception as e:
            check(False, f"{test.__name__} raised", f"{type(e).__name__}: {e}")

    # Generic coverage sweep (still non-destructive: destructive handlers are
    # in _SWEEP_SKIP), then reset the fixture in case the sweep left a mess.
    try:
        test_command_sweep(session)
    except Exception as e:
        check(False, "test_command_sweep raised", f"{type(e).__name__}: {e}")
    session = setup_session()

    # Destructive teardown last.
    try:
        test_delete(session)
    except Exception as e:
        check(False, "test_delete raised", f"{type(e).__name__}: {e}")

    total = _passes + len(_failures)
    print(f"\n{_passes}/{total} checks passed")
    if _failures:
        print("Failures:")
        for f in _failures:
            print(f"  - {f}")
    return not _failures


if __name__ == "__main__":
    sys.exit(0 if run_all() else 1)
