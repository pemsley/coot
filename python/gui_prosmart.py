
def add_module_prosmart():
    
    if (have_coot_python):
        if coot_python.main_menubar():
            menu = coot_menubar_menu("ProSMART")

            def prosmart_cut_to_func(sig_low, sig_high):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    set_extra_restraints_prosmart_sigma_limits(aa_imol,
                                                               sig_low, sig_high)
                    
            add_simple_coot_menu_menuitem(
                menu, "Cut to 6",
                lambda func: prosmart_cut_to_func(-6, 6))

            add_simple_coot_menu_menuitem(
                menu, "Cut to 4",
                lambda func: prosmart_cut_to_func(-4, 4))

            add_simple_coot_menu_menuitem(
                menu, "Cut to 2.5",
                lambda func: prosmart_cut_to_func(-2.5, 2.5))

            add_simple_coot_menu_menuitem(
                menu, "Cut to 2",
                lambda func: prosmart_cut_to_func(-2, 2))

            add_simple_coot_menu_menuitem(
                menu, "Cut to 1.5",
                lambda func: prosmart_cut_to_func(-1.5, 1.5))

            add_simple_coot_menu_menuitem(
                menu, "Cut to 1",
                lambda func: prosmart_cut_to_func(-1, 1))

            add_simple_coot_menu_menuitem(
                menu, "Cut to 0.5",
                lambda func: prosmart_cut_to_func(-0.5, 0.5))

            add_simple_coot_menu_menuitem(
                menu, "Cut to 0",
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

            def generate_self_restraint_func(sig):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    generate_local_self_restraints(aa_imol, aa_chain_id, sig)

            def generate_self_restraint_in_sphere_func(sig):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    current_residue_spec = [aa_chain_id, aa_res_no, aa_ins_code]
                    residue_specs = residues_near_residue(aa_imol, current_residue_spec, 6)
                    residue_specs.append(current_residue_spec)
                    generate_local_self_restraints_by_residues_py(aa_imol, residue_specs, sig)

            add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 4.3",
                lambda func: generate_self_restraint_func(4.3))

            add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 6",
                lambda func: generate_self_restraint_func(6))

            add_simple_coot_menu_menuitem(
                menu, "Generate Local Self Restraints 6",
                lambda func: generate_self_restraint_in_sphere_func(4.2))

                
            

