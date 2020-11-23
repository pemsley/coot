
def add_module_prosmart():
    
    if True:
        if coot_gui_api.main_menubar():
            menu = coot_gui.coot_menubar_menu("ProSMART")

            def generate_self_restraint_func(sig):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    coot.generate_local_self_restraints(aa_imol, aa_chain_id, sig)

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
                    coot.generate_local_self_restraints_by_residues_py(aa_imol,
                                                                  residue_specs,
                                                                  local_dist_max)


            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 4.3 for Chain",
                lambda func: generate_self_restraint_func(4.3))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 6 for Chain",
                lambda func: generate_self_restraint_func(6))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Generate Local Self Restraints 6",
                lambda func: generate_self_restraint_in_sphere_func())

            def display_extra_restraints_func(state):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    coot.set_show_extra_restraints(aa_imol, state)

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Undisplay Extra Restraints",
                lambda func: display_extra_restraints_func(0))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Display Extra Restraints",
                lambda func: display_extra_restraints_func(1))

                        
            def prosmart_cut_to_func(sig_low, sig_high):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    coot.set_extra_restraints_prosmart_sigma_limits(aa_imol,
                                                               sig_low, sig_high)
            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 6",
                lambda func: prosmart_cut_to_func(-6, 6))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 4",
                lambda func: prosmart_cut_to_func(-4, 4))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 2.0",
                lambda func: prosmart_cut_to_func(-2, 2))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 1.0",
                lambda func: prosmart_cut_to_func(-1, 1))
            
            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Undisplay All Extra Distance Restraints",
                lambda func: prosmart_cut_to_func(0, 0))
            
            def restraint_to_ca_func(state):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    coot.set_extra_restraints_representation_for_bonds_go_to_CA(aa_imol, state)

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Restraint Representation To CA",
                lambda func: restraint_to_ca_func(1))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Restraint Representation To Home Atom",
                lambda func: restraint_to_ca_func(0))
            
            ## extra
            
            # def delete_all_extra_restraints_func():
            #     with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
            #                                aa_ins_code, aa_atom_name, aa_alt_conf]:
            #         coot.delete_all_extra_restraints(aa_imol)

            # coot_gui.add_simple_coot_menu_menuitem(
            #     menu, "Delete All Extra Restraints",
            #     lambda func: delete_all_extra_restraints_func())

