
def add_module_prosmart():

    if (have_coot_python):
        if coot_python.main_menubar():
            menu = coot_menubar_menu("Restraints")

            def generate_all_molecule_self_restraints(val):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    generate_self_restraints(aa_imol, val)

            # Generate self resraints for chains (or all molecule)
            def generate_self_restraint_func(sig, all_mol=False):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    if not all_mol:
                        # chains
                        generate_local_self_restraints(aa_imol, aa_chain_id, sig)
                    else:
                        # all molecule
                        generate_local_self_restraints(aa_imol, sig)

            # Generate self restraints for residues around sphere
            def generate_self_restraint_in_sphere_func():
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    centred_residue = [aa_chain_id, aa_res_no, aa_ins_code]
                    radius = 10
                    local_dist_max = 4.2
                    other_residues = residues_near_residue(aa_imol,
                                                           centred_residue,
                                                           radius)
                    residue_specs = (centred_residue + other_residues) if isinstance(other_residues, list) else centred_residue
                    generate_local_self_restraints_by_residues_py(aa_imol,
                                                                  residue_specs,
                                                                  local_dist_max)

            add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 3.7 for Chain",
                lambda func: generate_self_restraint_func(3.7))

            add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 4.3 for Chain",
                lambda func: generate_self_restraint_func(4.3))

            add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 5 for Chain",
                lambda func: generate_self_restraint_func(6))

            add_simple_coot_menu_menuitem(
                menu, "Generate All-Molecule Self Restraints 4.3",
                lambda func: generate_self_restraint_func(4.3, True))

            add_simple_coot_menu_menuitem(
                menu, "Generate Local Self Restraints 6",
                lambda func: generate_self_restraint_in_sphere_func())

            add_simple_coot_menu_menuitem(
                menu, "Generate All-molecule Self Restraints 4.3",
                lambda func: generate_all_molecule_self_restraints(4.3))

            add_simple_coot_menu_menuitem(
                menu, "Generate All-molecule Self Restraints 5.0",
                lambda func: generate_all_molecule_self_restraints(5.0))

            def display_extra_restraints_func(state):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    set_show_extra_restraints(aa_imol, state)

            add_simple_coot_menu_menuitem(
                menu, "Undisplay Extra Restraints",
                lambda func: display_extra_restraints_func(0))

            add_simple_coot_menu_menuitem(
                menu, "Display Extra Restraints",
                lambda func: display_extra_restraints_func(1))


            def prosmart_cut_to_func(sig_low, sig_high):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    set_extra_restraints_prosmart_sigma_limits(aa_imol,
                                                               sig_low, sig_high)
            add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 6",
                lambda func: prosmart_cut_to_func(-6, 6))

            add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 4",
                lambda func: prosmart_cut_to_func(-4, 4))

            add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 2.0",
                lambda func: prosmart_cut_to_func(-2, 2))

            add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 1.0",
                lambda func: prosmart_cut_to_func(-1, 1))

            add_simple_coot_menu_menuitem(
                menu, "Undisplay All Extra Distance Restraints",
                lambda func: prosmart_cut_to_func(0, 0))

            def restraint_to_ca_func(state):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    set_extra_restraints_representation_for_bonds_go_to_CA(aa_imol, state)

            add_simple_coot_menu_menuitem(
                menu, "Restraint Representation To CA",
                lambda func: restraint_to_ca_func(1))

            add_simple_coot_menu_menuitem(
                menu, "Restraint Representation To Home Atom",
                lambda func: restraint_to_ca_func(0))

            ## extra

            # def delete_all_extra_restraints_func():
            #     with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
            #                                aa_ins_code, aa_atom_name, aa_alt_conf]:
            #         delete_all_extra_restraints(aa_imol)

            # add_simple_coot_menu_menuitem(
            #     menu, "Delete All Extra Restraints",
            #     lambda func: delete_all_extra_restraints_func())

