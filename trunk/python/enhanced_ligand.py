
ligand_check_refmac_mtz = False
ligand_check_refmac_fobs_col = False
ligand_check_refmac_sigfobs_col = False
ligand_check_refmac_r_free_col = False

def ligand_check_refmac_columns(f_list, sigf_list, rfree_list):

    # dummy function, not needed anywhere (as is)
    # need the above as globals...

    # happy path

    # Using the first sigf (there should only be one typically)
    # get the F label (Fx) from x/y/SIGFPx
    #
    pass
    

if enhanced_ligand_coot_p():


    def import_from_3d_generator_from_mdl():
        pass
    
    if (use_gui_qm != 2):
        menu = coot_menubar_menu("Ligand")
        
        add_simple_coot_menu_menuitem(
          menu,
          "SMILES -> 2D",
          lambda func:
          generic_single_entry("SMILES string",
                               "", " Send to 2D Viewer ",
                               lambda text: smiles_to_ligand_builder(text)
                               )
          )


        add_simple_coot_menu_menuitem(
            menu,
            "Residue -> 2D",
            lambda func:
            using_active_atom(residue_to_ligand_builder,
                              "aa_imol", "aa_chain_id", "aa_res_no",
                              "aa_ins_code", 0.015)
            )


        def flev_rdkit_func():
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                fle_view_with_rdkit(aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, 4.2)
                set_flev_idle_ligand_interactions(1)
            
        add_simple_coot_menu_menuitem(
            menu,
            "FLEV this residue",
            lambda func: flev_rdkit_func()
            )


        add_simple_coot_menu_menuitem(
            menu,
            "Toggle Ligand Interactions",
            lambda func: toggle_idle_ligand_interactions()
            )


        def go_solid_func(state):
            set_display_generic_objects_as_solid(state)
            graphics_draw()

        add_simple_coot_menu_menuitem(
            menu, "Solid Generic Objects",
            lambda func: go_solid_func(1))

        add_simple_coot_menu_menuitem(
            menu, "Unsolid Generic Objects",
            lambda func: go_solid_func(0))
        

        def show_chem_func():
            set_display_generic_objects_as_solid(1) # there may be consequences...
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                show_feats(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)
            
        add_simple_coot_menu_menuitem(
            menu,
            "Show Chemical Features",
            lambda func: show_chem_func()
            )

        # not in scheme?!
        def write_sdf_func():
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                rn = residue_name(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)
                file_name = rn + ".sdf"
                residue_to_sdf_file(aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                    file_name)

        add_simple_coot_menu_menuitem(
            menu,
            "write sdf file",
            lambda func: write_sdf_func()
            )

