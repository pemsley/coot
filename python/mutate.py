# mutate.py
# Copyright 2005, 2006, 2009 by Bernhard Lohkamp
# Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import coot
import coot_utils

# Mutate chain-id of molecule number imol given sequence.
#
# The number of residues in chain-id must match the length of sequence.
#


def mutate_chain(imol, chain_id, sequence):

    if (len(sequence) != chain_n_iresidues(chain_id, imol)):
        print("sequence mismatch: molecule", coot.chain_n_residues(
            chain_id, imol), "new sequences:", len(sequence))
    else:
        coot.make_backup(imol)  # do backup first
        backup_mode = coot.backup_state(imol)
        # turn off backup for imol
        coot.turn_off_backup(imol)
        baddies = 0
        for ires in sequence.upper():
            res = coot.mutate_single_residue_by_serial_number(
                ires, chain_id, imol, sequence_list[ires])
            if (res != 1):
                baddies += 1
        print("multi_mutate of %s residues had %s baddies" %
              (len(sequence), baddies))

        coot.set_have_unsaved_changes(imol)
        if (backup_mode == 1):
            coot.turn_on_backup(imol)
        coot.update_go_to_atom_window_on_changed_mol(imol)
        coot.graphics_draw()

# an internal function of mutate, This presumes a protein sequence


def multi_mutate(mutate_function, imol, start_res_no, chain_id, sequence):

    if len(sequence) > 0:
        baddies = 0
        for ires in range(len(sequence)):
            result = mutate_function(
                start_res_no+ires, "", chain_id, imol, sequence[ires].upper())
            if (result != 1):
                # add a baddy if result was 0 (fail)
                baddies += 1
            print("multi_mutate of", len(sequence),
                  "residues had", baddies, "errors")
    else:
        print("BL WARNING:: no sequence or sequence of no length found!")

# The stop-res-no is inclusive, so usage e.g.
# mutate_resiude_range(0,"A",1,2,"AC")
#
# This presumes a protein sequence (not nucleic acid).
# not any more... sort of
#


def mutate_residue_range(imol, chain_id, start_res_no, stop_res_no, sequence):

    if coot_utils.is_nucleotide_chain_qm(imol, chain_id):
        mutate_nucleotide_range(
            imol, chain_id, start_res_no, stop_res_no, sequence)
    else:
        coot.make_backup(imol)
        n_residues = stop_res_no - start_res_no + 1
        if (len(sequence) != n_residues):
            print("sequence length mismatch:", len(sequence), n_residues)
        else:
            backup_mode = coot.backup_state(imol)
            coot.turn_off_backup(imol)
            multi_mutate(coot.mutate_single_residue_by_seqno, imol,
                         start_res_no, chain_id, sequence.upper())
            coot.set_have_unsaved_changes(imol)
            if backup_mode == 1:
                coot.turn_on_backup(imol)
                coot.update_go_to_atom_window_on_changed_mol(imol)
                coot.graphics_draw()

# mutate and auto fit a residue range
#
# This presumes a protein sequence (not nucleic acid).
#
# The sequence is a string of one letter codes
#


def mutate_and_autofit_residue_range(imol, chain_id, start_res_no, stop_res_no,
                                     sequence):

    # does nucleic acids now too
    mutate_residue_range(imol, chain_id, start_res_no,
                         stop_res_no, sequence.upper())
    mol_for_map = coot.imol_refinement_map()
    if (mol_for_map >= 0):
        backup_mode = coot.backup_state(imol)
        coot.turn_off_backup(imol)
        for ires in range(len(sequence)):
            clash = 1
            altloc = ""
            inscode = ""
            resno = ires + start_res_no
            print("auto-fit-best-rotamer ", resno, altloc, inscode, chain_id,
                  imol, mol_for_map, clash)
            score = coot.auto_fit_best_rotamer(resno, altloc, inscode, chain_id,
                                               imol, mol_for_map, clash, 0.5)
            print("   Best score: ", score)
#          coot_utils.number_list(start_res_no,stop_res_no)
        if (backup_mode == 1):
            coot.turn_on_backup(imol)
    else:
        print("WARNING:: no map set for refinement.  Can't fit")

# mutate and autofit whole chain
#
# This presumes a protein sequence (not nucleic acid).


def mutate_and_auto_fit(residue_number, chain_id, mol, mol_for_map,
                        residue_type):

    coot.mutate(mol, chain_id, residue_number, "", residue_type)
    coot.auto_fit_best_rotamer(
        residue_number, "", "", chain_id, mol, mol_for_map, 0, 0.5)

# a short-hand for mutate-and-auto-fit


def maf(*args):
    mutate_and_auto_fit(*args)

# return a char (well string in python), return "A" for unknown residue_type
#


def three_letter_code2single_letter(residue_type):

    dic_type_1lc = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
                    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    if (residue_type in dic_type_1lc):
        res_type_1lc = dic_type_1lc[residue_type]
    else:
        # if not standard aa then set it to Ala
        res_type_1lc = "A"
    return res_type_1lc


# Mutate a residue range of nucleotides
# where sequence is a string (for example: "atgccgta")
#
def mutate_nucleotide_range(imol, chain_id, resno_start, resno_end, sequence):

    from types import StringType
    # return False on unknown letter. We don't need to be specific about
    # RNA/DNA here because mutate-base uses canonical_base_name() to
    # set the base name according to the residue type before it was
    # mutated.
    #

    def nucleotide_letter2three_letter_code(letter):
        # not entirely sure what should happening here,
        # only upper case conversion!?
        nt = ["A", "C", "G", "T", "U"]
        up = letter.upper()
        if up in nt:
            return up
        else:
            return False

    # main line
    if (type(sequence) is not StringType):
        s = "sequence must be a string"
        coot.info_dialog(s)
        print(s)
    else:
        residue_range_length = resno_end - resno_start + 1
        seq_length = len(sequence)

        if not (residue_range_length == seq_length):
            s = "residue range must equal sequence length"
            coot.info_dialog(s)
            print(s)
        else:
            residue_types = list(map(nucleotide_letter2three_letter_code,
                                     sequence.upper()))
            for res_no, res_type in map(None,
                                        list(range(resno_start, resno_end + 1)),
                                        residue_types):
                if res_type:
                    coot.mutate_base(imol, chain_id, res_no, "", res_type)

# a wrapper for mutate_single_residue_by_seqno (which uses slightly
# inconvenient single letter code)
#
# Here residue-type is the 3-letter code
# BL says:: not working, since mutate expects a 3 letter code. Maybe
# this function should work the other way round, i.e. take a 1 letter code
# and translate to 1 letter code for "new" mutate function.
#


def mutate_residue_redundant(residue_number, chain_id, imol,
                             residue_type, insertion_code):

    ins_code = insertion_code
    dic_type_1lc = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
                    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    if (residue_type in dic_type_1lc):
        res_type_1lc = dic_type_1lc[residue_type]
    else:
        # if not standard aa then set it to Ala
        res_type_1lc = "A"
    coot.mutate(imol, chain_id, residue_number, "", res_type_1lc)

# Prompted by Tim Gruene's email to CCP4bb 20060201.
# Turn all residues (including GLY) of imol to ALA.
#
# type is an optional argument.  if type is 'SER' then build polySer,
# if type is 'GLY', build polyGly.
#
# 1 or 2 args are:
#
# 1.) imol
#
# 2.) [optional] type, if type is "SER" then build polySer,
#                      if type is "GLY", build polyGly.
#


def poly_ala(imol, res_type=False):

    if coot.is_valid_model_molecule(imol):
        coot.make_backup(imol)
        backup_mode = coot.backup_state(imol)
        imol_map = coot.imol_refinement_map()

        coot.turn_off_backup(imol)
        dic_sg = {"SER": "S", "Ser": "S", "GLY": "G", "Gly": "G"}
        if (res_type):
            # residue name given!
            if ((len(res_type) == 1) and (res_type in list(dic_sg.values()))):
                single_letter_code = res_type
            else:
                if (res_type in dic_sg):
                    single_letter_code = dic_sg[res_type]
                else:
                    print("BL WARNING:: unrecognised residue name %s given" % res_type)
                    return
        else:
            single_letter_code = "A"

        for chain_id in coot_utils.chain_ids(imol):
            n_residues = coot.chain_n_residues(chain_id, imol)
            for serial_number in range(n_residues):
                coot.mutate_single_residue_by_serial_number(serial_number,
                                                       chain_id,
                                                       imol,
                                                       single_letter_code)
        if (backup_mode == 1):
            coot.turn_on_backup(imol)

    else:
        print("BL WARNING:: no valid model molecule", imol)

# Delete (back to the CB stub) the side change in the range
# resno-start to resno-end
#


def delete_sidechain_range(imol, chain_id, resno_start, resno_end):

    if coot_utils.valid_model_molecule_qm(imol):
        coot.make_backup(imol)
        backup_mode = coot.backup_state(imol)
        coot.turn_off_backup(imol)
        for resno in coot_utils.number_list(resno_start, resno_end):
            coot.delete_residue_sidechain(imol, chain_id, resno, "", 0)
        if (backup_mode == 1):
            coot.turn_on_backup(imol)
