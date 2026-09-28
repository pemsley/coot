# coot_commands/playbooks.py
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

"""What to do next: short procedures the assistant can follow.

:mod:`coot_commands.tools` tells the model what each command *does*.  This
module tells it what to *do* - the handful of moves an experienced person
makes when they see a rotamer outlier, a clash, or a green blob in the
difference map.  That knowledge is not in any command's help text, because
it is not about one command; it is about which command to reach for, in
what order, and how to tell whether it worked.

Why not put it in the system prompt?  Because a 7-8B model follows a short
prompt far better than a long one, and only one or two of these procedures
matter to any given request.  So they are retrieved the same way commands
are (:mod:`coot_commands.retrieval` - the same embedding endpoint, the same
:class:`~coot_commands.retrieval.ToolRetriever`), and only the relevant ones
are injected, as a system message, for that request.

Each playbook also declares the *tools* its steps mention.  This matters:
tool retrieval narrows ~90 commands to a couple of dozen, and a procedure
that says "now call backrub_residue" is worse than useless if
``backrub_residue`` was not among them.  :func:`tools_for` gives the agent
the names to pin, so a selected playbook always arrives with the commands
it asks for.

Writing a playbook
------------------
Keep the steps imperative, numbered, and short, and name real commands.
Always include a measurement step - a procedure the model cannot check the
result of is one it will report as successful whatever happened.  State the
give-up condition too: "if that does not clear it, say so" is a step, and
leaving it out is how you get a model that mutates a residue three times
rather than admit the density is ambiguous.
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from coot_commands.retrieval import ToolRetriever, default_prefixes, ollama_embed


class Playbook:
    """One situation and the procedure for it.

    *situation* is the one-line "when this applies"; *triggers* is the
    vocabulary a user or a validation result might use for it (it exists
    purely to give the embedding something to match on); *steps* is the
    procedure; *tools* are the command names the steps mention, which the
    agent pins into the exposed tool set.
    """

    def __init__(self, name: str, situation: str, triggers: str, steps: str,
                 tools: Sequence[str] = ()) -> None:
        self.name = name
        self.situation = situation.strip()
        self.triggers = triggers.strip()
        self.steps = steps.strip()
        self.tools = tuple(tools)

    def document(self) -> str:
        """The text embedded to represent this playbook for retrieval."""
        return (f"{self.name.replace('-', ' ')}. {self.situation} "
                f"Symptoms: {self.triggers}. {self.steps}")

    def render(self) -> str:
        """The playbook as it is shown to the model."""
        return (f"PROCEDURE - {self.name.replace('-', ' ')}\n"
                f"When: {self.situation}\n"
                f"{self.steps}")


PLAYBOOKS: Tuple[Playbook, ...] = (

    Playbook(
        "rotamer-outlier",
        "A side chain sits in an unusual conformation - a rotamer outlier, "
        "meaning a rotamer probability below 0.3%.",
        "rotamer outlier, bad rotamer, unusual side chain, improbable "
        "rotamer, side chain in the wrong place, fix the side chain",
        "1. score_residue on the residue. Read all three numbers, not just "
        "the rotamer one: a low rotamer probability on its own is NOT enough "
        "to justify changing anything. It is worth fixing when the rotamer is "
        "an outlier AND the density fit is poor, or the residue is clashing. "
        "If the rotamer is low but the density fit is good and there are no "
        "clashes, the density is telling you the conformation is real - say "
        "so and leave it alone.\n"
        "2. The rotamer thresholds depend on the residue type and "
        "score_residue tells you which apply. Arginine and lysine spread "
        "their probability over thirty-odd rotamers so their best possible "
        "score is only about 9%; valine has three and reaches about 73%. "
        "Never compare a rotamer percentage against another residue type's.\n"
        "3. fix_rotamer on it. That auto-fits the best rotamer, then tries a "
        "backrub, and keeps whichever scores best - do not run autofit_rotamer "
        "and backrub_residue by hand unless fix_rotamer is unavailable.\n"
        "4. Read the numbers it returns. If the rotamer has reached its "
        "favoured band and the density fit went up, it worked - stop.\n"
        "5. If it is still an outlier, refine_sphere on the residue to let the "
        "surroundings relax, then score_residue again.\n"
        "6. If it still will not clear, do not keep trying. Say so, and say "
        "why it might be: weak or absent side-chain density (the side chain "
        "may be genuinely disordered), a clash with a neighbour, or the "
        "wrong residue type for the sequence.",
        tools=("score_residue", "fix_rotamer", "refine_sphere",
               "autofit_rotamer", "backrub_residue", "go_to_residue")),

    Playbook(
        "ramachandran-outlier",
        "A residue has an improbable backbone phi/psi combination - a "
        "Ramachandran outlier.",
        "ramachandran outlier, bad backbone, improbable phi psi, backbone "
        "geometry problem, rama outlier",
        "1. check_ramachandran to list the outliers and pick the worst.\n"
        "2. go_to_residue on it so you can see what you are doing.\n"
        "3. A Ramachandran outlier is usually a misplaced carbonyl, so try "
        "pepflip_residue on the residue - and if that does not help, on the "
        "residue before it, since the offending peptide may be on either "
        "side.\n"
        "4. refine_sphere afterwards to settle the surroundings.\n"
        "5. check_ramachandran again to confirm the outlier has gone. If it "
        "has not, the backbone may be genuinely misthreaded - report it "
        "rather than flipping peptides repeatedly.",
        tools=("check_ramachandran", "go_to_residue", "pepflip_residue",
               "refine_sphere", "score_residue", "undo")),

    Playbook(
        "clash",
        "Atoms are too close together - a steric clash or atom overlap.",
        "clash, clashes, atom overlap, atoms too close, bumping atoms, "
        "steric problem, bad contacts",
        "Do these in order, one tool call at a time.\n"
        "1. check_clashes to find the overlaps - unless validate_anomalies "
        "has already reported them, which runs the same search and fills the "
        "same list, in which case go straight to step 2.\n"
        "2. go_to_clash 1 to frame the worst one. Use this rather than "
        "go_to_residue: a clash is between two residues, and going to one of "
        "them can leave the other off screen.\n"
        "3. Work out which of the two residues can move, from the atom names "
        "the clash reports. N, CA, C, O and OXT are backbone atoms; every "
        "other name is a side-chain atom. A residue whose SIDE-CHAIN atom is "
        "in the clash can move - a side chain swings out of the way. A "
        "residue whose BACKBONE atom is in the clash cannot. So: if one atom "
        "is side-chain, that residue is your target; if BOTH are side-chain "
        "atoms, both residues are targets and you fix both; if both are "
        "backbone atoms, go to step 7. Say which you picked and why.\n"
        "4. score_residue on each target, one call each, before you change "
        "anything. Do this first even though fix_rotamer reports its own "
        "before-and-after: with two targets, fixing one moves atoms that the "
        "other is clashing with, so its numbers stop meaning what they meant "
        "and you need both baselines taken up front.\n"
        "5. fix_rotamer on the first target, then check_clashes. If the "
        "clash has gone, stop and report it. If there is a second target and "
        "the clash is still there, fix_rotamer on that one too, then "
        "check_clashes again.\n"
        "6. If it is still there, refine_sphere on the residue you fixed "
        "last, then check_clashes once more.\n"
        "7. Now stop - do not go back to step 5 or keep refining. Report "
        "which residues clash, their before and after scores, what you "
        "tried, and that it needs a human decision: usually two things "
        "modelled into the same density (a misplaced water, or one "
        "conformation where there are really two) or a backbone error that "
        "local fitting cannot reach.",
        tools=("check_clashes", "go_to_clash", "go_to_residue",
               "score_residue", "fix_rotamer", "refine_sphere",
               "delete_residue")),

    Playbook(
        "unmodelled-density",
        "There is density the model does not account for - a green (positive) "
        "difference-map peak, or an unmodelled blob.",
        "green blob, positive difference density, unmodelled density, "
        "something missing from the model, unexplained density, what is this "
        "blob, fo-fc peak",
        "1. validate_unmodelled_blobs to find the blobs, largest first.\n"
        "2. go_to_blob 1 to centre on the largest, then get_active_residue to "
        "find out what the blob is sitting next to. Judge by what is nearest "
        "the blob, not by which residue happens to be closest in sequence.\n"
        "3. Decide what it is. A small blob hydrogen-bonding distance from a "
        "polar atom is a water: add_water. A blob beside a residue with an "
        "incomplete side chain is that side chain: replace_residue with the "
        "same type to rebuild it, then fix_rotamer. A large blob in a pocket "
        "is a ligand: fit_ligand, once you know which one to look for.\n"
        "4. A blob at the END of a chain is the commonest case and the "
        "easiest to get wrong. Density just beyond the last residue means "
        "there are MORE RESIDUES TO BUILD - it does NOT mean a missing OXT. "
        "Do not reach for add_oxt here. Use add_terminal_residue, then "
        "refine_sphere.\n"
        "5. Then run validate_unmodelled_blobs AGAIN. If the blob is smaller "
        "but still there, the chain continues: add_terminal_residue once "
        "more and repeat. Build one residue at a time, refining and "
        "re-checking after each, until the blob has gone or has stopped "
        "shrinking. Do not add several residues before looking.\n"
        "6. add_oxt belongs at the very end, and only once there is no blob "
        "left at the terminus - capping a chain that still has density beyond "
        "it is what stops you building the residues that are really there.\n"
        "7. Never guess a ligand identity. If you cannot tell what a blob is, "
        "describe its size and surroundings and ask.",
        tools=("validate_unmodelled_blobs", "go_to_blob", "get_active_residue",
               "add_water", "add_terminal_residue", "add_oxt", "fit_ligand",
               "add_monomer", "replace_residue", "fix_rotamer",
               "refine_sphere", "show_difference_map_peaks")),

    Playbook(
        "negative-density",
        "Modelled atoms sit in density that is not there - a red (negative) "
        "difference-map peak.",
        "red blob, negative difference density, atoms in no density, red "
        "peak, atom should not be there, over-modelled",
        "1. show_difference_map_peaks to mark the peaks, and go_to_residue on "
        "the residue at the worst negative one.\n"
        "2. score_residue to confirm the density fit is genuinely poor.\n"
        "3. A whole water in negative density is usually spurious: "
        "delete_residue.\n"
        "4. A side chain in negative density is usually disordered rather "
        "than absent - the standard treatment is to keep the atoms and let "
        "the B-factors rise (refine_b_factors), not to delete them. Deleting "
        "a side chain is a modelling decision; say what you would do and ask "
        "before doing it.\n"
        "5. If only part of a residue is in negative density, try fix_rotamer "
        "first - it may simply be pointing the wrong way.",
        tools=("show_difference_map_peaks", "go_to_residue", "score_residue",
               "delete_residue", "fix_rotamer", "refine_b_factors",
               "refine_sphere")),

    Playbook(
        "missing-atoms",
        "A residue is incomplete - it has fewer atoms than its type should.",
        "missing atoms, incomplete residue, truncated side chain, stub "
        "residue, missing side chain",
        "1. check_missing_atoms to list the affected residues.\n"
        "2. go_to_residue on one and look at whether there is density for the "
        "missing part.\n"
        "3. If there is, replace_residue with the same residue type rebuilds "
        "it complete, then fix_rotamer to place the side chain and "
        "refine_sphere to settle it.\n"
        "4. score_residue to confirm the rebuilt side chain fits.\n"
        "5. If there is no density for the missing atoms, leave the residue "
        "truncated - that is a deliberate choice, not an error - and say so.",
        tools=("check_missing_atoms", "go_to_residue", "replace_residue",
               "fix_rotamer", "refine_sphere", "score_residue")),

    Playbook(
        "cis-peptide",
        "A peptide bond is in the cis conformation.",
        "cis peptide, cis bond, twisted peptide, omega outlier",
        "1. check_cis_peptides to count and locate them.\n"
        "2. A cis peptide *before a proline* is common and usually correct - "
        "leave it alone and say why.\n"
        "3. For a non-proline cis peptide, go_to_residue and look at the "
        "density: genuine non-proline cis peptides exist but are rare, so "
        "the usual cause is a modelling error.\n"
        "4. pepflip_residue then refine_sphere, and check_cis_peptides again.\n"
        "5. If the density clearly supports the cis conformation, keep it and "
        "report it as a real feature.",
        tools=("check_cis_peptides", "go_to_residue", "pepflip_residue",
               "refine_sphere")),

    Playbook(
        "suspicious-water",
        "A water looks wrong - too many close contacts, or sitting in a "
        "position no water should occupy.",
        "bad water, suspicious water, water might be an ion, highly "
        "coordinated water, too many waters, check the waters",
        "1. check_waters to find waters with unusually high coordination.\n"
        "2. A water with five or more close contacts is often a "
        "misassigned ion - sodium, magnesium and calcium all get modelled as "
        "water. Look at the coordination number and distances before "
        "deciding.\n"
        "3. go_to_residue on it and score_residue to see how well it fits.\n"
        "4. A water in weak density with no hydrogen-bonding partner should "
        "go: delete_residue.\n"
        "5. Do not swap a water for an ion on geometry alone - the "
        "identification needs the chemistry of the site and usually the "
        "B-factor behaviour too. Report the candidates and let the user "
        "decide.",
        tools=("check_waters", "go_to_residue", "score_residue",
               "delete_residue", "validate_unmodelled_blobs")),

    Playbook(
        "gln-asn-his-flip",
        "A glutamine, asparagine or histidine side chain may be flipped 180 "
        "degrees - the density cannot distinguish the amide's N from its O.",
        "gln flip, asn flip, his flip, amide flip, glutamine asparagine "
        "outlier, side chain flipped the wrong way",
        "1. check_gln_asn to find the candidates.\n"
        "2. The density looks identical either way round, so the evidence is "
        "hydrogen bonding: the correct orientation is the one whose N donates "
        "and whose O accepts to the neighbours.\n"
        "3. score_residue before and after any change - a flip should not "
        "make the density fit worse, and if it does, something else is "
        "wrong.\n"
        "4. refine_sphere afterwards.\n"
        "5. This is a judgement call about the local hydrogen-bonding "
        "network. Report which residues are flagged and what you changed, "
        "rather than flipping everything the analysis lists.",
        tools=("check_gln_asn", "go_to_residue", "score_residue",
               "refine_sphere")),

    Playbook(
        "triage-the-model",
        "The user asks what to do next, what is wrong with the model, or to "
        "improve it generally, without naming a specific problem.",
        "what should I do next, what is wrong with this model, improve the "
        "model, tidy up the structure, check my model, where should I look, "
        "fix the worst problems",
        "1. Start by measuring, not fixing: validate_anomalies for an "
        "overview, then validate_unmodelled_blobs for anything missing.\n"
        "2. Report what you found before changing anything, worst first.\n"
        "3. Work in the order that stops you undoing your own work: backbone "
        "before side chains (Ramachandran outliers and peptide flips), then "
        "side chains (rotamer outliers), then clashes, then unmodelled "
        "density, then waters.\n"
        "4. Fix a few of the worst, not all of them, and score_residue each "
        "one before and after so you can say what actually improved.\n"
        "5. Re-run validate_anomalies at the end and report the before and "
        "after counts. If a category did not improve, say so plainly.",
        tools=("validate_anomalies", "validate_unmodelled_blobs",
               "check_ramachandran", "check_rotamers", "check_clashes",
               "score_residue", "fix_rotamer", "go_to_residue",
               "open_validation")),

    Playbook(
        "poor-density-fit",
        "A residue does not fit its density - it may be shifted, pointing the "
        "wrong way, or in no density at all.",
        "residue does not fit the density, out of density, poor fit, bad "
        "density fit, residue in the wrong place, does not match the map",
        "1. score_residue to get a density fit number to compare against.\n"
        "2. refine_sphere on the residue - most poor fits are small shifts "
        "that real-space refinement pulls back in. Prefer it to refining the "
        "residue alone: a sphere refine also moves the residues that are "
        "close to this one in SPACE rather than in sequence, which is where "
        "the thing holding it out of density usually is.\n"
        "3. score_residue again. If the density fit rose, stop.\n"
        "4. If the side chain specifically is out of density, fix_rotamer.\n"
        "5. If the backbone is out of density over several residues, that is "
        "a bigger problem than these steps can fix - a register error or a "
        "misthreaded loop. Report it and where it starts, rather than "
        "refining residue by residue.",
        tools=("score_residue", "refine_sphere", "fix_rotamer",
               "refine_range", "go_to_residue", "get_active_map")),
)


# The heading that marks an injected playbook message, so the agent can find
# and replace the previous request's guidance instead of stacking a fresh copy
# on the conversation every turn.
PLAYBOOK_HEADER = (
    "Relevant procedures for this request. Follow the one that fits; if none "
    "does, ignore them. Prefer their steps to improvising, and always measure "
    "before and after so you can report what actually changed. Work through "
    "the steps one at a time: make the call the current step asks for, read "
    "what it returns, then move on. Do not restate the remaining steps or "
    "announce calls you are about to make - make the next one."
)


def _by_name() -> Dict[str, Playbook]:
    return {p.name: p for p in PLAYBOOKS}


def playbook_documents() -> Dict[str, str]:
    """Map playbook name -> the text embedded to retrieve it."""
    return {p.name: p.document() for p in PLAYBOOKS}


def tools_for(names: Iterable[str]) -> List[str]:
    """The command names the given playbooks refer to, de-duplicated in order.

    The agent pins these into the exposed tool set so a procedure never names
    a command the model has not been given.
    """
    by_name = _by_name()
    out: List[str] = []
    for name in names:
        playbook = by_name.get(name)
        if playbook is None:
            continue
        for tool in playbook.tools:
            if tool not in out:
                out.append(tool)
    return out


def render(names: Iterable[str]) -> str:
    """The system message carrying the given playbooks, or ``""`` for none."""
    by_name = _by_name()
    chosen = [by_name[n] for n in names if n in by_name]
    if not chosen:
        return ""
    return PLAYBOOK_HEADER + "\n\n" + "\n\n".join(p.render() for p in chosen)


# How similar a playbook has to be to the request before it is worth sending.
#
# Calibrated against embeddinggemma, where a request that really is about one
# of these situations scores 0.58-0.73 against it, a bare "what should I do
# next?" scores ~0.35 (short queries score low even when they are exactly on
# topic), and requests about something else entirely - "colour the map blue",
# "fetch 4hhb", "go to residue A 45" - land at 0.14-0.31.  0.32 sits in that
# gap.  A different embedding model will put the gap somewhere else; override
# with COOT_PLAYBOOK_MIN_SCORE, and note that setting it too high is the safe
# direction (no guidance, which is where this started) while setting it too
# low sends procedures to requests that have no problem to solve.
MIN_RELEVANCE = 0.32


def min_relevance() -> float:
    """The similarity floor, from ``COOT_PLAYBOOK_MIN_SCORE`` or the default."""
    raw = os.environ.get("COOT_PLAYBOOK_MIN_SCORE")
    if raw is None:
        return MIN_RELEVANCE
    try:
        return float(raw)
    except ValueError:
        return MIN_RELEVANCE


# Process-wide retriever over the playbooks, built lazily so importing this
# module never touches the network.  Separate from the command retriever
# because the two corpora are ranked independently - a request selects its own
# top few commands and its own top few procedures.
_default_retriever: Optional[ToolRetriever] = None


def default_retriever() -> ToolRetriever:
    """The shared retriever over all playbooks (Ollama embeddings)."""
    global _default_retriever
    if _default_retriever is None:
        query_prefix, document_prefix = default_prefixes()
        _default_retriever = ToolRetriever(
            playbook_documents(), ollama_embed,
            query_prefix=query_prefix, document_prefix=document_prefix)
    return _default_retriever


def select_playbooks(query: str, k: int = 2,
                     retriever: Optional[ToolRetriever] = None) -> List[str]:
    """The playbooks relevant to *query*: up to *k*, best first, possibly none.

    Unlike tool retrieval, which must always hand the model something to call,
    this returns nothing when nothing is relevant.  Ranking alone would attach
    guidance to every request - "go to residue A 45" would arrive with advice
    about missing atoms - so anything below :data:`MIN_RELEVANCE` is dropped.
    """
    if k <= 0:
        return []
    scored = (retriever or default_retriever()).select_scored(query, k)
    floor = min_relevance()
    return [name for name, score in scored if score >= floor]
