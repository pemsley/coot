import coot
import coot_utils

# this function is not at startup?

def deactivate_molecules_except(imol):
    for i in coot_utils.model_molecule_list():
        if (i != imol):
            coot.set_mol_active(i, 0)


# This will hydrogenate the active residue, not imol
# returns the generic_objects.probe clash score
#
def contact_score_ligand(imol, res_spec):

    deactivate_molecules_except(imol)

    chain_id = coot_utils.res_spec_to_chain_id(res_spec)
    res_no = coot_utils.res_spec_to_res_no(res_spec)
    ins_code = coot_utils.res_spec_to_ins_code(res_spec)
    ss = "//" + chain_id + "/" + str(res_no)
    imol_selection = coot.new_molecule_by_atom_selection(imol, ss)
    coot_molprobity_dir = coot_utils.get_directory("coot-molprobity")
    ligand_selection_pdb = os.path.join(coot_molprobity_dir,
                                        "tmp-selected-ligand-for-probe-" + str(imol) + ".pdb")
    protein_selection_pdb = os.path.join(coot_molprobity_dir,
                                         "tmp-protein-for-probe-" + str(imol) + ".pdb")
    dots_file_name = os.path.join(coot_molprobity_dir,
                                  "probe-" + chain_id + "-" + str(res_no) + ".dots")
    coot.set_mol_active(imol_selection, 0)
    coot.set_mol_displayed(imol_selection, 0)
    coot_utils.set_go_to_atom_molecule(imol)
    rc = residue_centre(imol, chain_id, res_no, ins_code)
    coot.set_rotation_centre(*rc)
    coot.hydrogenate_region(6)

    coot.write_pdb_file(imol_selection, ligand_selection_pdb)
    coot.write_pdb_file(imol, protein_selection_pdb)

    args = ["-q", "-unformated", "-once", # -once or -both
                                  # "-outside", # crashes generic_objects.probe
                                  # first pattern
                                  "CHAIN" + chain_id + " " + str(res_no),
                                  # second pattern
                                  "not " + str(res_no),
                                  # consider og33 (occ > 0.33)
                                  "-density50",
                                  ligand_selection_pdb, protein_selection_pdb]
    coot_utils.popen_command(probe_command, args, [], dots_file_name, False)

    # coot_utils.debugging!?
    coot.handle_read_draw_probe_dots_unformatted(dots_file_name, imol, 0)

    cs = generic_objects.probe_clash_score(dots_file_name)
    coot.graphics_draw()
    return cs


def contact_score_ligand_func():
    with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        contact_score_ligand(aa_imol, [aa_chain_id, aa_res_no, aa_ins_code])

def coot_contact_dots_ligand_func():
    with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        coot.coot_contact_dots_for_ligand_py(aa_imol, [aa_chain_id, aa_res_no, aa_ins_code])

# not ready for public yet
def coot_all_atom_contact_dots_func():
    
    with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:

        coot_probe_object_names = ['wide-contact', 'close-contact', 'small-overlap',
                                   'big-overlap', 'H-bond', 'clashes']
        n = coot.number_of_generic_objects()
        for i in range(n):
           if generic_object_name(i) in coot_probe_object_names:
              coot.close_generic_object(i)
           else:
               print(generic_object_name(i), "is not in", coot_probe_object_names)
        coot.coot_all_atom_contact_dots(aa_imol)

