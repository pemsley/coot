# test_coot_commands.py
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

"""Standalone tests for the coot_commands package (no Coot required).

Run from the python/ directory:

    python3 test_coot_commands.py

The command modules do ``try: import coot / except ImportError: coot = None``,
so they load without Coot; the argument-value tests inject a small fake ``coot``
so that model/map completion is deterministic.  The suite is also discoverable
by pytest (the ``test_*`` functions) if it is ever available.
"""

import contextlib

import coot_commands  # noqa: F401  - triggers command discovery/registration
from coot_commands import types
from coot_commands.commands import ligand as ligand_mod
from coot_commands.commands import model_edit as model_edit_mod
from coot_commands.commands import models as models_mod
from coot_commands.commands import refine as refine_mod
from coot_commands.commands import session as session_mod
from coot_commands.commands import settings as settings_mod
from coot_commands.commands import validation as validation_mod
from coot_commands.completion import complete
from coot_commands.registry import unmatched_examples
import coot_command_interface as cli


class _FakeCoot:
    """Just enough of the coot API for the completion value providers.

    Models live at molecule numbers 0, 1, 2; a map lives at 3.
    """

    def graphics_n_molecules(self):
        return 4

    def is_valid_model_molecule(self, i):
        return 1 if i in (0, 1, 2) else 0

    def is_valid_map_molecule(self, i):
        return 1 if i == 3 else 0


class _FakeCootFull(_FakeCoot):
    """A richer fake for exercising command *behaviour* via dispatch.

    Adds molecule names, an active model/map, and the handful of scoring
    and search calls the ligand/blob/list commands make, recording their
    arguments so tests can assert what the command asked Coot to do.
    """

    NAMES = {0: "tutorial-modern", 1: "1hr2_final",
             2: "ligand-blob", 3: "tutorial-map"}

    def __init__(self):
        self.centre = None
        self.search = {}
        self.displayed = {}

    # --- names -----------------------------------------------------------
    def molecule_name_stub_py(self, imol, include_path_flag):
        return self.NAMES.get(imol, "")

    def molecule_name(self, imol):
        return self.NAMES.get(imol, "")

    # --- active model / map ---------------------------------------------
    def active_residue_py(self):
        return [0, "A", 10, "", " CA "]

    def imol_refinement_map(self):
        return 3

    def map_is_difference_map(self, imol):
        return 0

    # --- unmodelled blobs -----------------------------------------------
    def find_blobs_py(self, imol_model, imol_map, cut_off):
        # Deliberately not volume-sorted, so the command's own sort is tested.
        return [[[1.0, 2.0, 3.0], 5.0], [[4.0, 5.0, 6.0], 9.0]]

    def set_rotation_centre(self, x, y, z):
        self.centre = (x, y, z)

    # --- residue edits (OXT / terminal residue) --------------------------
    def add_OXT_to_residue(self, imol, chain_id, resno, ins_code):
        self.search["oxt"] = (imol, chain_id, resno, ins_code)
        return 1

    def add_terminal_residue(self, imol, chain_id, resno, restype, immediate):
        self.search["terminal"] = (imol, chain_id, resno, restype, immediate)
        return 1

    def backrub_rotamer(self, imol, chain_id, resno, ins_code, alt_conf):
        self.search["backrub"] = (imol, chain_id, resno, ins_code)

    # --- superposition ---------------------------------------------------
    def superpose(self, imol1, imol2, move_flag):
        self.search["superpose"] = (imol1, imol2, move_flag)

    # --- ligand fitting --------------------------------------------------
    def get_monomer(self, comp_id):
        self.search["comp_id"] = comp_id
        return 4

    def add_ligand_clear_ligands(self):
        self.search["cleared"] = True

    def set_ligand_search_protein_molecule(self, imol):
        self.search["protein"] = imol

    def set_ligand_search_map_molecule(self, imol):
        self.search["map"] = imol

    def add_ligand_search_ligand_molecule(self, imol):
        self.search["ligand"] = imol

    def set_find_ligand_here_cluster(self, state):
        self.search["here"] = state

    def execute_ligand_search_py(self):
        return [5]

    def set_mol_displayed(self, imol, state):
        self.displayed[imol] = state


class _FakeSettings(_FakeCootFull):
    """Records settings writes and serves them back through the getters.

    Each ``set_*`` stores into ``self.settings``; the matching ``get_*``
    reads it back, so a set-then-get round trip returns the value the
    command wrote - just as the real getters read Coot's live state.
    """

    def __init__(self):
        super().__init__()
        self.settings = {}
        # Per-map contour levels, keyed by molecule number.
        self.absolute = {}
        self.sigma = {}
        self.drawn = 0

    def graphics_draw(self):
        self.drawn += 1

    # --- global scalar settings -----------------------------------------
    def set_default_bond_thickness(self, t):
        self.settings["bond_thickness"] = t

    def get_default_bond_thickness(self):
        return self.settings.get("bond_thickness", 5)

    def set_font_size(self, i):
        self.settings["font_size"] = i

    def get_font_size(self):
        return self.settings.get("font_size", 2)

    def set_map_radius(self, f):
        self.settings["map_radius"] = f

    def get_map_radius(self):
        return self.settings.get("map_radius", 10.0)

    def set_map_sampling_rate(self, r):
        self.settings["sampling"] = r

    def get_map_sampling_rate(self):
        return self.settings.get("sampling", 1.8)

    def set_iso_level_increment(self, v):
        self.settings["iso_step"] = v

    def get_iso_level_increment(self):
        return self.settings.get("iso_step", 0.05)

    def set_diff_map_iso_level_increment(self, v):
        self.settings["diff_step"] = v

    def get_diff_map_iso_level_increment(self):
        return self.settings.get("diff_step", 0.05)

    # --- per-map contour level ------------------------------------------
    def set_contour_level_absolute(self, imol, level):
        self.absolute[imol] = level

    def get_contour_level_absolute(self, imol):
        return self.absolute.get(imol, 0.0)

    def set_contour_level_in_sigma(self, imol, level):
        self.sigma[imol] = level

    def get_contour_level_in_sigma(self, imol):
        return self.sigma.get(imol, 0.0)


class _FakeValidation(_FakeCootFull):
    """Records validation/updating-map calls and serves sample scorer data.

    Molecules: models 0,1,2; a 2Fo-Fc (data) map at 3 and an Fo-Fc
    difference map at 4.  The scriptable scorers return one outlier each so
    the summary commands have something to report.
    """

    def __init__(self):
        super().__init__()
        self.updating = None
        self.stopped = None
        self.overlay = None
        self.peaks = None
        self.gln_asn = None
        self.blobs_call = None

    # --- molecule inventory (add the difference map at 4) ----------------
    def graphics_n_molecules(self):
        return 5

    def is_valid_map_molecule(self, i):
        return 1 if i in (3, 4) else 0

    def map_is_difference_map(self, imol):
        return 1 if imol == 4 else 0

    # --- updating maps ---------------------------------------------------
    def set_auto_updating_sfcalc_genmap(self, imol_model, imol_data, imol_diff):
        self.updating = (imol_model, imol_data, imol_diff)

    def stop_updating_molecule(self, imol):
        self.stopped = imol

    # --- GUI validation openers -----------------------------------------
    def dynamic_validation_internal(self, imol, imol_map):
        self.overlay = (imol, imol_map)

    def difference_map_peaks(self, imol_map, imol_coords, level, closeness,
                             pos, neg, around):
        self.peaks = (imol_map, imol_coords, level)

    def gln_asn_b_factor_outliers(self, imol):
        self.gln_asn = imol

    def find_blobs_py(self, imol_model, imol_map, cut_off):
        self.blobs_call = (imol_model, imol_map, cut_off)
        return [[[1.0, 2.0, 3.0], 5.0], [[4.0, 5.0, 6.0], 9.0]]

    # --- scriptable scorers (one outlier each) --------------------------
    def all_molecule_ramachandran_score_py(self, imol):
        # index 5 holds the per-residue scores: [_, name, probability, ...].
        return [0, 0, 0, 0, 0, [[0, "A 45", 0.01], [0, "A 46", 0.50]]]

    def molecule_atom_overlaps_py(self, imol):
        return [{"overlap-volume": 3.0}, {"overlap-volume": 1.0}]

    def rotamer_graphs_py(self, imol):
        return [["A", 45, "", 1.2, "mt"], ["A", 46, "", 88.0, "m"]]

    def cis_peptides_py(self, imol):
        return [["cis-a"], ["cis-b"]]

    def non_standard_residue_names_py(self, imol):
        return ["LIG", "LIG", "NAG"]

    def missing_atom_info_py(self, imol):
        return [["A", 89], ["B", 12]]

    def highly_coordinated_waters_py(self, imol, coordination, dist):
        return [["A", 201]]


@contextlib.contextmanager
def _use_coot(fake, *modules):
    """Install *fake* as ``coot`` in ``types`` and each command *module*.

    Command handlers read their module-global ``coot`` at call time, so a
    dispatch-level test must patch the module that owns the handler as well
    as ``types`` (which resolves the active model/map).  State is restored
    on exit even if the test raises.
    """
    mods = [types, *modules]
    saved = [m.coot for m in mods]
    for m in mods:
        m.coot = fake
    try:
        yield fake
    finally:
        for m, s in zip(mods, saved):
            m.coot = s


def test_examples_match_their_own_pattern():
    """Every example must be dispatchable by its own command."""
    problems = unmatched_examples()
    assert problems == [], f"examples that don't match their pattern: {problems}"


def test_keyword_completion_fills_unambiguous_word():
    assert complete("sh") == ("show ", [])
    assert complete("show mo") == ("show model ", [])


def test_keyword_completion_lists_options_and_shared_prefix():
    replacement, options = complete("show ")
    assert options == ["environment", "map", "model", "only", "residue",
                       "ribbons", "surface", "symmetry"]
    assert replacement == "show "
    # "view" must NOT leak here from the "perspective view" example.
    assert "view" not in options


def test_keyword_completion_extends_shared_prefix():
    # Both "map" and "model" start with "m", so Tab fills the shared "m".
    replacement, options = complete("show m")
    assert options == ["map", "model"]
    assert replacement == "show m"


def test_colour_value_completion():
    _, options = complete("background ")
    assert options == types.colour_names()
    assert complete("colour map 1 bl") == ("colour map 1 bl", ["black", "blue"])


def test_projection_literal_completion():
    assert complete("ort") == ("orthographic ", [])
    assert complete("per") == ("perspective ", [])


def test_no_completion_for_unknown_input():
    assert complete("wibble") == ("", [])
    assert cli.complete_command("wibble") == ""


def test_model_value_completion_with_fake_coot():
    saved = types.coot
    types.coot = _FakeCoot()
    try:
        # Numbers sort numerically, not lexically.
        replacement, options = complete("show model ")
        assert options == ["0", "1", "2"], options
        assert replacement == "show model "
        # A single matching value is filled in.
        assert complete("show model 1") == ("show model 1 ", [])
        # Maps and models are distinct sets.
        assert complete("show map ") == ("show map 3 ", [])
    finally:
        types.coot = saved


def test_interface_protocol_string():
    # First line is the replacement; a second line lists ambiguous options.
    assert cli.complete_command("sh") == "show "
    assert cli.complete_command("show ") == \
        "show \nenvironment  map  model  only  residue  ribbons  " \
        "surface  symmetry"


def test_model_completion_shows_molecule_names():
    # Tab-listing a model/map argument decorates each number with its name,
    # while the value inserted into the line stays the bare number.
    with _use_coot(_FakeCootFull()):
        replacement, options = complete("show model ")
        assert options == ["0 (tutorial-modern)", "1 (1hr2_final)",
                           "2 (ligand-blob)"], options
        assert replacement == "show model "
        # Filling in a single unambiguous value uses the number, not the label.
        assert complete("show model 1") == ("show model 1 ", [])


def test_new_command_keyword_completion():
    # list -> the two list commands' nouns.
    assert complete("list ")[1] == ["maps", "models", "molecules"]
    # fit -> only "ligand", so it fills in.
    assert complete("fit ") == ("fit ligand ", [])
    # find and "go to" offer the new nouns among their options.
    assert "blobs" in complete("find ")[1]
    assert "blob" in complete("go to ")[1]


def test_list_models_and_maps_report_names():
    with _use_coot(_FakeCootFull(), session_mod):
        models = cli.run_command("list models")
        assert "3 model(s):" in models
        assert "0: tutorial-modern" in models
        maps = cli.run_command("list maps")
        assert "1 map(s):" in maps
        assert "3: tutorial-map" in maps


def test_list_commands_without_coot():
    # Handlers must degrade gracefully when Coot is unavailable.
    assert cli.run_command("list models") == "No models loaded"
    assert cli.run_command("list maps") == "No maps loaded"


def test_find_blobs_then_go_to_blob():
    fake = _FakeCootFull()
    with _use_coot(fake, validation_mod):
        out = cli.run_command("find blobs")
        assert "2 unmodelled blob(s)" in out
        # Largest blob (volume 9) is numbered 1.
        assert "blob 1 at (4.0, 5.0, 6.0)" in out
        go = cli.run_command("go to blob 1")
        assert fake.centre == (4.0, 5.0, 6.0)
        assert "Centred on blob 1" in go
        # Out-of-range blob is a clean error, not a crash.
        assert "does not exist" in cli.run_command("go to blob 9")


def test_find_blobs_defaults_to_difference_map():
    # With a difference map loaded (molecule 4), a bare "find blobs" must use
    # it at the 3 sigma default - not the 2Fo-Fc map at a low cut-off, which
    # is what produced blobs "in the middle of nowhere".
    fake = _FakeValidation()
    with _use_coot(fake, validation_mod):
        out = cli.run_command("find blobs")
        assert fake.blobs_call == (0, 4, 3.0)  # model 0, difference map 4, 3 sigma
        assert "on map 4" in out
        # An explicit cut-off and map still win.
        cli.run_command("find blobs using map 3 above 2 sigma")
        assert fake.blobs_call == (0, 3, 2.0)


def test_find_blobs_warns_without_difference_map():
    # No difference map: fall back to the active map, but say so.
    fake = _FakeCootFull()
    fake.find_blobs_py = lambda m, mp, c: []
    with _use_coot(fake, validation_mod):
        out = cli.run_command("find blobs")
        assert "no difference map loaded" in out


def test_go_to_blob_without_blobs_errors():
    validation_mod._last_blobs = []
    out = cli.run_command("go to blob 1")
    assert "find blobs" in out


def test_fit_ligand_runs_search():
    fake = _FakeCootFull()
    with _use_coot(fake, ligand_mod):
        out = cli.run_command("fit ligand LIG")
        assert "Fitted ligand LIG" in out
        assert "5" in out  # the fitted solution's molecule number
        assert fake.search["comp_id"] == "LIG"
        assert fake.search["protein"] == 0    # active model
        assert fake.search["map"] == 3        # refinement map
        assert fake.search["here"] == 0       # whole-map search
        # The "here" variant restricts the search to the view centre.
        here = cli.run_command("fit ligand ATP here")
        assert fake.search["here"] == 1
        assert "at the view centre" in here


def test_add_oxt_to_named_residue():
    fake = _FakeCootFull()
    with _use_coot(fake, model_edit_mod):
        out = cli.run_command("add OXT to A/89")
        assert "Added OXT to A/89" in out
        # Model defaults to the active model (0); ins code blank.
        assert fake.search["oxt"] == (0, "A", 89, "")


def test_add_oxt_defaults_to_active_residue():
    # With no residue named, the command must act on the active residue,
    # which the fake reports as A/10 of model 0.
    fake = _FakeCootFull()
    with _use_coot(fake, model_edit_mod):
        out = cli.run_command("add OXT")
        assert fake.search["oxt"] == (0, "A", 10, "")
        assert "Added OXT to A/10" in out


def test_add_terminal_residue_named_and_typed():
    fake = _FakeCootFull()
    with _use_coot(fake, model_edit_mod):
        # Default type is "auto".
        cli.run_command("add terminal residue to A/89")
        assert fake.search["terminal"] == (0, "A", 89, "auto", 1)
        # An explicit "as ALA" forces the residue type.
        out = cli.run_command("add terminal residue to A 89 as ALA")
        assert fake.search["terminal"] == (0, "A", 89, "ALA", 1)
        assert "(ALA)" in out


def test_add_terminal_residue_defaults_to_active():
    fake = _FakeCootFull()
    with _use_coot(fake, model_edit_mod):
        out = cli.run_command("add terminal residue")
        # Active residue is A/10 of model 0.
        assert fake.search["terminal"] == (0, "A", 10, "auto", 1)
        assert "Added a terminal residue to A/10" in out


def test_backrub_without_spec_uses_active_residue():
    # Regression: a per-residue command with no residue spec must still match
    # its own command and act on the active residue (A/10 of model 0) - not
    # fall through to the "back" navigation command.
    fake = _FakeCootFull()
    with _use_coot(fake, model_edit_mod):
        out = cli.run_command("backrub rotamer")
        assert fake.search["backrub"] == (0, "A", 10, "")
        assert "Backrub-fitted A/10" in out


def test_superpose_copies_by_default():
    fake = _FakeCootFull()
    with _use_coot(fake, models_mod):
        out = cli.run_command("superpose model 0 onto model 1")
        # superpose(reference=1, moving=0, flag=1 -> copy).
        assert fake.search["superpose"] == (1, 0, 1)
        assert "copy of model 0 onto model 1" in out
        # "in place" moves the source instead of copying (flag 0).
        moved = cli.run_command("superpose model 0 onto model 1 in place")
        assert fake.search["superpose"] == (1, 0, 0)
        assert "(moved)" in moved


def test_superpose_onto_itself_errors():
    assert "onto itself" in cli.run_command("superpose model 2 onto model 2")


def test_fit_ligand_unknown_code_errors():
    class _NoDict(_FakeCootFull):
        def get_monomer(self, *_):
            return -1
    with _use_coot(_NoDict(), ligand_mod):
        out = cli.run_command("fit ligand ZZZ")
        assert "could not get a dictionary" in out


def test_scalar_settings_round_trip():
    # Each 'set X' writes the value; the matching 'get X' reads it back.
    fake = _FakeSettings()
    with _use_coot(fake, settings_mod):
        assert "Set bond thickness to 4" in cli.run_command("set bond thickness to 4")
        assert fake.settings["bond_thickness"] == 4
        assert cli.run_command("get bond thickness") == "Bond thickness is 4"
        # "what is the ..." and a bare noun are the same query.
        assert cli.run_command("what is the bond thickness") == "Bond thickness is 4"
        assert cli.run_command("bond thickness") == "Bond thickness is 4"

        cli.run_command("set font size 3")
        assert fake.settings["font_size"] == 3
        assert cli.run_command("get font size") == "Font size is 3"

        cli.run_command("set map radius to 20")
        assert fake.settings["map_radius"] == 20.0
        assert cli.run_command("get map radius") == "Map radius is 20 A"

        cli.run_command("set map sampling rate to 2.5")
        assert fake.settings["sampling"] == 2.5
        assert cli.run_command("get map sampling rate") == "Map sampling rate is 2.5"

        cli.run_command("set contour step to 0.1")
        assert fake.settings["iso_step"] == 0.1
        assert cli.run_command("get contour step") == "Contour step is 0.1"

        cli.run_command("set difference map contour step to 0.2")
        assert fake.settings["diff_step"] == 0.2
        assert cli.run_command("get difference map contour step") == \
            "Difference map contour step is 0.2"


def test_set_settings_rejects_non_numeric_value():
    # The value must be a number: the regex requires a digit, so a worded
    # value simply doesn't match and falls through to "unrecognised".
    assert "Unrecognised command" in cli.run_command("set font size big")


def test_contour_level_absolute_and_sigma():
    fake = _FakeSettings()
    with _use_coot(fake, settings_mod):
        # Named map, absolute units.
        out = cli.run_command("set contour level of map 3 to 0.3")
        assert fake.absolute[3] == 0.3
        assert "Set contour level of map 3 to 0.3" in out
        assert fake.drawn > 0  # setting the level redraws
        # No map named -> the active (refinement) map, in sigma.
        out = cli.run_command("set contour level to 1.5 sigma")
        assert fake.sigma[3] == 1.5
        assert "Set contour level of map 3 to 1.5 sigma" in out
        # get reports both the absolute level and its sigma equivalent.
        assert cli.run_command("get contour level") == \
            "Contour level of map 3 is 0.3 (1.5 sigma)"


def test_contour_level_on_a_model_is_rejected():
    # Molecule 0 is a model, not a map: a clear error rather than a bad call.
    fake = _FakeSettings()
    with _use_coot(fake, settings_mod):
        out = cli.run_command("set contour level of map 0 to 0.3")
        assert "is a model, not a map" in out
        assert 0 not in fake.absolute


def test_get_setting_without_coot():
    # With no Coot API, a getter degrades to a clear message, not an error.
    assert cli.run_command("get bond thickness") == \
        "Bond thickness is unavailable (no Coot API)"


def test_open_and_close_sequence():
    fake = _FakeCootFull()
    calls = {}
    fake.sequence_view = lambda imol: calls.setdefault("open", imol)
    fake.remove_sequence_view_from_sequence_view_box = \
        lambda imol: calls.setdefault("close", imol)
    with _use_coot(fake, session_mod):
        assert "Opened the sequence of model 0" in cli.run_command("open sequence")
        assert calls["open"] == 0
        assert "Closed the sequence of model 0" in cli.run_command("close sequence")
        assert calls["close"] == 0
        # "show sequence" is no longer the verb; it should not open the view.
        assert "Unrecognised command" in cli.run_command("show sequence")


def test_updating_maps_on_uses_data_and_difference_maps():
    fake = _FakeValidation()
    with _use_coot(fake, settings_mod):
        out = cli.run_command("set updating maps on")
        # (active model 0, refinement/data map 3, difference map 4).
        assert fake.updating == (0, 3, 4)
        assert "Turned on updating maps" in out
        # Both "off" and the bare "stop" turn it off again.
        assert "Turned off" in cli.run_command("set updating maps off")
        assert fake.stopped == 0
        fake.stopped = None
        assert "Turned off" in cli.run_command("stop updating maps")
        assert fake.stopped == 0


def test_updating_maps_needs_a_difference_map():
    # _FakeCootFull has only a non-difference map (3), so there is nothing to
    # update into: a clear error, and no sfcalc wired up.
    fake = _FakeValidation()
    fake.map_is_difference_map = lambda imol: 0  # hide the difference map
    with _use_coot(fake, settings_mod):
        out = cli.run_command("set updating maps on")
        assert "no difference map" in out
        assert fake.updating is None


def test_text_validation_summaries():
    fake = _FakeValidation()
    with _use_coot(fake, validation_mod):
        assert "1 Ramachandran outlier" in cli.run_command("check ramachandran")
        assert "A/45" in cli.run_command("check rotamers")
        assert "1 clash" in cli.run_command("check clashes")
        assert "2 cis peptide" in cli.run_command("check cis peptides")
        nonstd = cli.run_command("check non-standard residues")
        assert "LIG" in nonstd and "NAG" in nonstd
        assert "A/89" in cli.run_command("check missing atoms")
        assert "highly-coordinated water" in cli.run_command("check waters")


def test_gui_validation_openers():
    fake = _FakeValidation()
    with _use_coot(fake, validation_mod):
        out = cli.run_command("open validation")
        assert fake.overlay == (0, 3)  # active model, active map
        assert "validation overlay" in out
        peaks = cli.run_command("difference map peaks above 4 sigma")
        assert fake.peaks == (4, 0, 4.0)  # diff map, model, level
        assert "difference map peaks" in peaks
        assert "Opened Gln/Asn" in cli.run_command("check gln and asn")
        assert fake.gln_asn == 0


def test_refine_b_factors_runs_shiftfield():
    fake = _FakeCootFull()
    calls = {}
    fake.imol_refinement_map = lambda: 3
    fake.shiftfield_b_factor_refinement = lambda imol: calls.__setitem__("bf", imol)
    with _use_coot(fake, refine_mod):
        out = cli.run_command("refine b factors")
        assert calls["bf"] == 0          # active model
        assert "Refined B-factors of model 0" in out
        # Alternative spellings reach the same command.
        cli.run_command("refine adps of model 1")
        assert calls["bf"] == 1


def test_refine_residue_centres_on_the_target_first():
    fake = _FakeCootFull()
    calls = {}
    fake.imol_refinement_map = lambda: 3
    fake.set_go_to_atom_molecule = lambda imol: calls.__setitem__("goto_mol", imol)
    fake.set_go_to_atom_from_res_spec_py = lambda spec: calls.setdefault("goto_spec", spec) or 1
    fake.refinement_immediate_replacement_state = lambda: 0
    fake.set_refinement_immediate_replacement = lambda s: None
    fake.refine_zone = lambda imol, ch, lo, hi, ins: calls.__setitem__("zone", (imol, ch, lo, hi))
    fake.accept_regularizement = lambda: None
    with _use_coot(fake, refine_mod):
        out = cli.run_command("refine A 45")
    # The view was centred on A/45 before the refinement ran.
    assert calls["goto_spec"] == ["A", 45, ""]
    assert calls["zone"] == (0, "A", 45, 45)
    assert "Refined A/45" in out


def test_refine_b_factors_needs_a_map():
    fake = _FakeCootFull()
    fake.imol_refinement_map = lambda: -1   # no refinement map set
    with _use_coot(fake, refine_mod):
        out = cli.run_command("refine b-factors")
        assert "no map set for refinement" in out


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
