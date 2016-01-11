
def deactivate_molecules_except(imol):
    for i in model_molecule_list():
        if (i != imol):
            set_molecule_active(i, 0)


# This will hydrogenate the active residue, not imol
#
def contact_score_ligand(imol, res_spec):

    deactivate_molecules_except(imol)

    chain_id = res_spec2chain_id(res_spec)
    res_no = res_spec2res_no(res_spec)
    ins_code = res_spec2ins_code(res_spec)
    ss = "//" + chain_id + "/" + str(res_no)
    imol_selection = new_molecule_by_atom_selection(imol, ss)
    ligand_selection_pdb = os.path.join("coot-molprobity",
                                        "tmp-selected-ligand-for-probe-" + str(imol) + ".pdb")
    protein_selection_pdb = os.path.join("coot-molprobity",
                                         "tmp-protein-for-probe-" + str(imol) + ".pdb")
    dots_file_name = os.path.join("coot-molprobity",
                                  "probe-" + chain_id + "-" + str(res_no) + ".dots")
    set_mol_active(imol_selection, 0)
    set_mol_displayed(imol_selection, 0)
    set_go_to_atom_molecule(imol)
    rc = residue_centre(imol, chain_id, res_no, ins_code)
    set_rotation_centre(*rc)
    hydrogenate_region(6)

    write_pdb_file(imol_selection, ligand_selection_pdb)
    write_pdb_file(imol, protein_selection_pdb)
    popen_command(probe_command, ["-q", "-u", "-once", # -once or -both
                                  # first pattern
                                  "CHAIN" + chain_id + " " + str(res_no),
                                  # second pattern
                                  "not " + str(res_no),
                                  # consider og33 (occ > 0.33)
                                  "-density30",
                                  ligand_selection_pdb, protein_selection_pdb],
                  [], dots_file_name, False)

    # debugging!?
    handle_read_draw_probe_dots_unformatted(dots_file_name, imol, 0)

    probe_clash_score(dots_file_name)
    graphics_draw()


menu = coot_menubar_menu("Ligand")

def contact_score_ligand_func():
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        contact_score_ligand(aa_imol, [aa_chain_id, aa_res_no, aa_ins_code])
        
add_simple_coot_menu_menuitem(menu, "Isolated dots for this ligand",
                              lambda func: contact_score_ligand_func())

            
