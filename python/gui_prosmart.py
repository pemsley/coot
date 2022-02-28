
# target is my molecule, ref is the homologous (high-res) model
#
# extra arg: include_side_chains=False
#
def run_prosmart(imol_target, imol_ref, include_side_chains=False):
  """
  target is my molecule, ref is the homologous (high-res) model

  extra arg: include_side_chains=False
  """

  dir_stub = "coot-ccp4"
  make_directory_maybe(dir_stub)
  target_pdb_file_name = os.path.join(dir_stub,
                                      molecule_name_stub(imol_target, 0).replace(" ", "_") + \
                                      "-prosmart.pdb")
  reference_pdb_file_name = os.path.join(dir_stub,
                                         molecule_name_stub(imol_ref, 0).replace(" ", "_") + \
                                         "-prosmart-ref.pdb")
  prosmart_out = os.path.join("ProSMART_Output",
                              molecule_name_stub(imol_target, 0).replace(" ", "_") + \
                              "-prosmart.txt")
  prosmart_rmax = 6.0

  write_pdb_file(imol_target, target_pdb_file_name)
  write_pdb_file(imol_ref, reference_pdb_file_name)
  prosmart_exe = find_exe("prosmart")
  if prosmart_exe:
    l = ["-p1", target_pdb_file_name,
         "-p2", reference_pdb_file_name,
         "-restrain_seqid", "30"]
    if include_side_chains:
      l += ["-side"]
    popen_command(prosmart_exe,
                  l,
                  [],
                  os.path.join(dir_stub, "prosmart.log"),
                  False)
    if (not os.path.isfile(prosmart_out)):
      print "file not found", prosmart_out
    else:
      print "Reading ProSMART restraints from", prosmart_out
      add_refmac_extra_restraints(imol_target, prosmart_out)

def add_prosmart_secondary_structure_restraints(imol, do_mc_h_bonds_also_flag):

    import os

    dir_stub = get_directory("coot-ccp4")
    stub_name = molecule_name_stub(imol, 0)

    helix_pdb_file_name_rwd = stub_name + "-helix.pdb"
    helix_pdb_file_name = os.path.join(dir_stub, helix_pdb_file_name_rwd)
    strand_pdb_file_name_rwd = stub_name + "-strand.pdb"
    strand_pdb_file_name = os.path.join(dir_stub, strand_pdb_file_name_rwd)
    h_bonds_pdb_file_name_rwd = stub_name + "-h-bonds.pdb"
    h_bonds_pdb_file_name = os.path.join(dir_stub, h_bonds_pdb_file_name_rwd)
    helix_out = os.path.join(dir_stub, "ProSMART_Output",
                             "LIB_" + stub_name + "-helix.txt")
    strand_out = os.path.join(dir_stub, "ProSMART_Output",
                              "LIB_" + stub_name + "-strand.txt")
    h_bonds_out = os.path.join(dir_stub, "ProSMART_Output",
                               stub_name + "-h-bonds.txt")  # not LIB_

    write_pdb_file(imol, helix_pdb_file_name)
    write_pdb_file(imol, strand_pdb_file_name)
    if do_mc_h_bonds_also_flag:
        write_pdb_file(imol, h_bonds_pdb_file_name)

    # ProSMART writes results in ProSMART_Output, so we change to the coot-ccp4 directory
    # so that it puts it in the place we expect it

    current_dir = os.getcwd()

    os.chdir(dir_stub)

    prosmart_exe = find_exe("prosmart")
    if prosmart_exe:
        for file_name, prosm_arg, secondary_structure in [[helix_pdb_file_name_rwd,
                                                           "-helix", "helix"],
                                                          [strand_pdb_file_name_rwd,
                                                           "-strand", "strand"]]:
            # print "BL DEBUG:: prosmart cmd", prosmart_exe, ["-p1", file_name, prosm_arg]
            popen_command(prosmart_exe,
                          ["-p1", file_name, prosm_arg],
                          [],
                          "prosmart-" + stub_name + "-" + secondary_structure + ".log",
                          True)
        if do_mc_h_bonds_also_flag:
            popen_command(prosmart_exe,
                          ["-p1", h_bonds_pdb_file_name_rwd, "-h"],
                          [],
                          "prosmart-" + stub_name + "-h-bonds.log",
                          False)
    os.chdir(current_dir)

    for fn in [helix_out, strand_out]:
        print "INFO:: reading ProSMART output file fn:", fn
        if os.path.isfile(fn):
            add_refmac_extra_restraints(imol, fn)
        else:
            s = "Missing file: " + fn
            info_dialog(s)

    if do_mc_h_bonds_also_flag:
        print "INFO:: reading ProSMART output file fn:", h_bonds_out
        if os.path.isfile(h_bonds_out):
            add_refmac_extra_restraints(imol, h_bonds_out)
        else:
            s = "Missing file: " + h_bonds_out
            info_dialog(s)


def add_module_restraints():

    if have_coot_python:
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
                menu, "Generate Self Restraints 3.7 for Chain",
                lambda func: generate_self_restraint_func(3.7))

            add_simple_coot_menu_menuitem(
                menu, "Generate Self Restraints 4.3 for Chain",
                lambda func: generate_self_restraint_func(4.3))

            add_simple_coot_menu_menuitem(
                menu, "Generate Self Restraints 5 for Chain",
                lambda func: generate_self_restraint_func(5))

            add_simple_coot_menu_menuitem(
                menu, "Generate All-molecule Self Restraints 4.3",
                lambda func: generate_all_molecule_self_restraints(4.3))

            add_simple_coot_menu_menuitem(
                menu, "Generate All-molecule Self Restraints 5.0",
                lambda func: generate_all_molecule_self_restraints(5.0))

            add_simple_coot_menu_menuitem(
                menu, "Generate All-molecule Self Restraints 6.0",
                lambda func: generate_all_molecule_self_restraints(6.0))

            add_simple_coot_menu_menuitem(
                menu, "Generate Local Self Restraints 6",
                lambda func: generate_self_restraint_in_sphere_func())

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

            add_simple_coot_menu_menuitem(
                menu, "Add Intermediate Atom Rotamer Dodecs",
                lambda func: set_show_intermediate_atoms_rota_markup(1))

            add_simple_coot_menu_menuitem(
                menu, "Add Intermediate Atom Ramachandran Spheres",
                lambda func: set_show_intermediate_atoms_rama_markup(1))

            add_coot_menu_separator(menu)
            load_from_search_load_path("user_define_restraints.py")
            add_coot_menu_separator(menu)

            def delete_all_extra_restraints_func():
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    delete_all_extra_restraints(aa_imol)
            add_simple_coot_menu_menuitem(
                menu, "Delete All Extra Restraints",
                lambda func: delete_all_extra_restraints_func())

def add_module_prosmart():

    if (have_coot_python):
        if coot_python.main_menubar():
            menu = coot_menubar_menu("ProSMART")

            def add_prosmart_ss_restraints_func(with_H_bonds=False):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    add_prosmart_secondary_structure_restraints(aa_imol, with_H_bonds)
            add_simple_coot_menu_menuitem(
                menu, "Add ProSMART Secondary Structure & H-bond Restraints",
                lambda func: add_prosmart_ss_restraints_func(True))

            add_simple_coot_menu_menuitem(
                menu, "Add ProSMART (Only) Secondary Structure Restraints",
                lambda func: add_prosmart_ss_restraints_func())

            add_simple_coot_menu_menuitem(
                menu,
                "Save as REFMAC restraints...",
                lambda func:
                generic_chooser_and_file_selector("Save REFMAC restraints for molecule",
                                                  valid_model_molecule_qm,
                                                  " Restraints file name:  ",
                                                  "refmac-restraints.txt",
                                                  lambda imol, file_name:
                                                  extra_restraints2refmac_restraints_file(imol, file_name)))


            def launch_prosmart_gui():
                def go_button_cb(*args):
                    imol_tar = get_option_menu_active_molecule(*option_menu_mol_list_pair_tar)
                    imol_ref = get_option_menu_active_molecule(*option_menu_mol_list_pair_ref)
                    do_side_chains = check_button.get_active()
                    run_prosmart(imol_tar, imol_ref, do_side_chains)
                    window.destroy()
                window = gtk.Window(gtk.WINDOW_TOPLEVEL)
                hbox = gtk.HBox(False, 0)
                vbox = gtk.VBox(False, 0)
                h_sep = gtk.HSeparator()
                chooser_hint_text_1 = " Target molecule "
                chooser_hint_text_2 = " Reference (high-res) molecule "
                go_button = gtk.Button(" ProSMART ")
                cancel_button = gtk.Button("  Cancel  ")
                check_button = gtk.CheckButton("Include Side-chains")

                option_menu_mol_list_pair_tar = generic_molecule_chooser(vbox,
                                                                         chooser_hint_text_1)
                option_menu_mol_list_pair_ref = generic_molecule_chooser(vbox,
                                                                         chooser_hint_text_2)

                vbox.pack_start(check_button, False, False, 2)
                vbox.pack_start(h_sep, False, False, 2)
                vbox.pack_start(hbox, False, False, 2)
                hbox.pack_start(go_button, False, False, 6)
                hbox.pack_start(cancel_button, False, False, 6)
                window.add(vbox)

                cancel_button.connect("clicked", lambda w: window.destroy())

                go_button.connect("clicked", go_button_cb, option_menu_mol_list_pair_tar,
                                  option_menu_mol_list_pair_ref)
                window.show_all()

            add_simple_coot_menu_menuitem(
                menu,
                "ProSMART...",
                lambda func: launch_prosmart_gui()
            )


            # def restraint_to_ca_func(state):
            #     with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
            #                                aa_ins_code, aa_atom_name, aa_alt_conf]:
            #         set_extra_restraints_representation_for_bonds_go_to_CA(aa_imol, state)

            # add_simple_coot_menu_menuitem(
            #     menu, "Restraint Representation To CA",
            #     lambda func: restraint_to_ca_func(1))

            # add_simple_coot_menu_menuitem(
            #     menu, "Restraint Representation To Home Atom",
            #     lambda func: restraint_to_ca_func(0))


