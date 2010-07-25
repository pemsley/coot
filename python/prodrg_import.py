cprodrg = "cprodrg"
cprodrg = "c:/Programs/CCP4-Packages/ccp4-6.1.3/bin/cprodrg.exe"

# if there is a prodrg-xyzin set the current-time to its mtime, else False
prodrg_xyzin      = "../../coot/lbg/prodrg-in.mdl"
sbase_to_coot_tlc = "../../coot/lbg/.sbase-to-coot-comp-id"
# this is for BL win machine
prodrg_xyzin      = os.path.join(os.getenv("HOME"), "Projects/build-xp-python/lbg/prodrg-in.mdl")
sbase_to_coot_tlc = "../../build-xp-python/lbg/.sbase-to-coot-comp-id" 

def import_from_prodrg(minimize_mode):

    import operator
    
    prodrg_dir = "coot-ccp4"
    res_name = "DRG"

    make_directory_maybe(prodrg_dir)
    prodrg_xyzin  = "../../coot/lbg/prodrg-in.mdl"
    # BL win machine
    prodrg_xyzin      = os.path.join(os.getenv("HOME"), "Projects/build-xp-python/lbg/prodrg-in.mdl")
    prodrg_xyzout = os.path.join(prodrg_dir, "prodrg-" + res_name + ".pdb")
    prodrg_cif    = os.path.join(prodrg_dir, "prodrg-out.cif")
    prodrg_log    = os.path.join(prodrg_dir, "prodrg.log")
    # requires python >= 2.5 (shall we test?)
    mini_mode = "NO" if (minimize_mode == 'mini-no') else "PREP"
    # should test if cprdrg exists?
    status = popen_command(cprodrg,
                           ["XYZIN",  prodrg_xyzin,
                            "XYZOUT", prodrg_xyzout,
                            "LIBOUT", prodrg_cif],
                           ["MINI " + mini_mode, "END"],
                           prodrg_log, True)
    if operator.isNumberType(status):
        if (status == 0):
            read_cif_dictionary(prodrg_cif)
            imol = handle_read_draw_molecule_and_move_molecule_here(prodrg_xyzout)
            # FIXME
            # this wont work - have to rework using atom to allow extra args before (and after) dictionary?
            using_active_atom(match_ligand_torsions, imol, "aa_imol", "aa_chain_id", "aa_res_no")
            using_active_atom(overlap_ligands, imol, "aa_imol", "aa_chain_id", "aa_res_no")
            with_auto_accept([regularize_residues, imol, [["", 1, ""]]])
            using_active_atom(match_ligand_torsions, imol, "aa_imol", "aa_chain_id", "aa_res_no")
            using_active_atom(overlap_ligands, imol, "aa_imol", "aa_chain_id", "aa_res_no")
            additional_representation_by_attributes(imol, "", 1, 1, "", 2, 2, 0.2, 1)
            return True

if (have_coot_python):
    if coot_python.main_menubar():
        menu = coot_menubar_menu("PRODRG")
        add_simple_coot_menu_menuitem(
            menu,
            "Import (using MINI PREP)",
            lambda func:
            # run prodrg, read its output files, and run regularisation
            # on the imported PDB file
            import_from_prodrg('mini-prep')
            )

        add_simple_coot_menu_menuitem(
            menu,
            "Import (no pre-minimisation)",
            lambda func:
            # run prodrg, read its output files, [and run regularisation
            # on the imported PDB file] ; BL says:: wrong comment
            import_from_prodrg('mini-no')
            )

        add_simple_coot_menu_menuitem(
            menu,
            "Export to lbg",
            lambda func:
            using_active_atom(prodrg_flat, "aa_imol", "aa_chain_id", "aa_res_no")
            )

        add_simple_coot_menu_menuitem(
            menu,
            "FLE-View",
            lambda func:
            using_active_atom(fle_view, "aa_imol", "aa_chain_id", "aa_res_no", "aa_ins_code")
            )
        add_simple_coot_menu_menuitem(
            menu,
            "Load SBase monomer...",
            lambda func:
            generic_single_entry("Load SBase Monomer from three-letter-code: ",
                                 "",
                                 " Load ",
                                 lambda tlc:
                                 get_sbase_monomer(tlc))
            )


def get_mdl_latest_time(file_name):
    if not os.path.isfile(file_name):
        return False
    else:
        return os.stat(file_name).st_mtime

global mdl_latest_time
global sbase_transfer_latest_time
mdl_latest_time = get_mdl_latest_time(prodrg_xyzin)
sbase_transfer_latest_time = get_mdl_latest_time(sbase_to_coot_tlc)
# FIXME: this is not a proper name
def func():
    import operator
    global mdl_latest_time
    global sbase_transfer_latest_time
    
    mdl_now_time   = get_mdl_latest_time(prodrg_xyzin)
    sbase_now_time = get_mdl_latest_time(sbase_to_coot_tlc)

    # print "sbase_now_time %s    sbase_latest_time %s" %(sbase_now_time, sbase_transfer_latest_time)

    if operator.isNumberType(mdl_now_time):
        if (mdl_now_time > mdl_latest_time):
            mdl_latest_time = mdl_now_time  # globals? FIXME
            import_from_prodrg('mini-prep')

    if operator.isNumberType(sbase_now_time):
        if (sbase_now_time > sbase_transfer_latest_time):
            sbase_transfer_latest_time = sbase_now_time  # globals? FIXME
            fin = open(sbase_to_coot_tlc)
            tlc_symbol = fin.readline()  # need to read more? FIXME
            fin.close()
            imol = get_sbase_monomer(tlc_symbol)
            if not valid_model_molecule_qm(imol):
                print "failed to get SBase molecule for", tlc_symbol
            else:
                # it was read OK, do an overlap
                using_active_atom(overlap_ligands, imol, "aa_imol", "aa_chain_id", "aa_res_no")
                
    return True # return value, keep running; FIXME:: how to stop?
gobject.timeout_add(500, func)

# return False (if fail) or a list of: the molecule number of the
# selected residue, the prodrg output mol file-name, the prodrg
# output pdb file-name
#
def prodrg_flat(imol_in, chain_id_in, res_no_in):

    import operator
    import random

    selection_string = "//" + chain_id_in + "/" + str(res_no_in)
    imol = new_molecule_by_atom_selection(imol_in, selection_string)
    prodrg_input_file_name = os.path.join("coot-ccp4", "tmp-residue-for-prodrg.pdb")
    prodrg_output_mol_file = os.path.join("coot-ccp4", ".coot-to-lbg-mol")
    prodrg_output_pdb_file = os.path.join("coot-ccp4", ".coot-to-lbg-pdb")
    prodrg_output_lib_file = os.path.join("coot-ccp4", ".coot-to-lbg-lib")
    prodrg_log             = os.path.join("coot-ccp4", "tmp-prodrg-flat.log")
    set_mol_displayed(imol, 0)
    set_mol_active   (imol, 0)
    write_pdb_file(imol, prodrg_input_file_name)
    arg_list = ["XYZIN",  prodrg_input_file_name,
                "MOLOUT", prodrg_output_mol_file,
                "XYZOUT", prodrg_output_pdb_file,
                "LIBOUT", prodrg_output_lib_file]
    print "arg_list", arg_list
    status = popen_command(cprodrg,
                           arg_list,
                           ["COORDS BOTH", "MINI FLAT", "END"],
                           prodrg_log, True)
    # Does this make sense in python? status being number that is. FIXME
    if not operator.isNumberType(status):
        info_dialog("Ooops: cprodrg not found?")
        return False
    else:
        if not (status == 0):
            mess = "Something went wrong running cprodrg\n" + \
                   ("(quelle surprise)" if random.randint(0,100) <10 else "")
            info_dialog(mess)
            return False
        else:
            # normal return value (hopefully)
            return [imol,
                    prodrg_outpu_mol_file,
                    prodrg_output_pdb_file,
                    prodrg_output_lib_file]


def prodrg_plain(mode, imol_in, chain_id_in, res_no_in):

    selection_string = "//" + chain_id_in + "/" + \
                       str(res_no_in)
    imol = new_molecule_by_atom_selection(imol_in, selection_string)
    stub = os.path.join("coot-ccp4", "prodrg-tmp-" + str(os.getpid()))
    prodrg_xyzin  = stub + "-xyzin.pdb"
    prodrg_xyzout = stub + "-xyzout.pdb"
    prodrg_cif    = stub + "-dict.cif"
    prodrg_log    = stub + ".log"

    write_pdb_file(imol, prodrg_xyzin)
    result = popen_command(cprodrg,
                           ["XYZIN",  prodrg_xyzin,
                            "XYZOUT", prodrg_xyzout,
                            "LIBOUT", prodrg_cif],
                           ["MINI PREP", "END"],
                           prodrg_log, True)
    close_molecule(imol)
    return [result, prodrg_xyzout, prodrg_cif]


def fle_view(imol, chain_id, res_no, ins_code):

    import operator
    
    r_flat  = prodrg_flat (imol, chain_id, res_no)
    r_plain = prodrg_plain('mini-no', imol, chain_id, res_no)
    if (r_flat and
        operator.isNumberType(r_plain[0]) and
        r_plain[0] == 0 ):
        imol_ligand_fragment = r_flat[0]
        prodrg_output_flat_mol_file_name = r_flat[1]
        prodrg_output_flat_pdb_file_name = r_flat[2]
        prodrg_output_cif_file_name      = r_flat[3]
        prodrg_output_3d_pdb_file_name   = r_plain[1]
        # 'using_active_atom'
        active_atom = active_residue()
        imol     = active_atom[0]
        chain_id = active_atom[1]
        res_no   = active_atom[2]
        fle_view_internal(imol, chain_id, res_no, "",  # should be from active_atom!!     using_active_atom([[]])
                          imol_ligand_fragment,
                          prodrg_output_flat_mol_file_name,
                          prodrg_output_flat_pdb_file_name,
                          prodrg_output_3d_pdb_file_name,
                          prodrg_output_cif_file_name)
        # touch on Windows?? FIXME
        # either distribute touch.exe or use DOS:
        # copy /b filename.ext +,,
        popen_command("touch",
                      [os.path.join("coot-ccp4", ".coot-to-lbg-mol-ready")],
                      [],
                      "/dev/null", False)

