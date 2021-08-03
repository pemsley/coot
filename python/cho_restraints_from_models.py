
# the smaller the sigmas the more weight the restraints have
# This pushes up the weight a bit (c.f. 1.0)

global cho_bond_sigma_scale
cho_bond_sigma_scale = 3.0

def get_glyco_model_dir():
    import os
    return os.path.join(pkgdatadir(), "data", "cho-models")

def get_glyco_tree_acedrg_pyranose_dictionaries_dir():
    import os
    return os.path.join(pkgdatadir(), "data", "cho-acedrg")

def remove_duplicates(ls):
    return list(set(ls))

def read_acedrg_pyranose_dictionaries():
    import os
    cho_list = ["NAG", "MAN", "BMA", "SIA", "GLC", "FUC", "XYP"]
    for comp_id in cho_list:
        fn = os.path.join(get_glyco_tree_acedrg_pyranose_dictionaries_dir(),
                          comp_id + "-acedrg.cif")
        if os.path.isfile(fn):
            coot.read_cif_dictionary(fn)

def add_cho_restraints_for_residue(imol, residue_spec):
    if isinstance(residue_spec, list):
        id = glyco_tree_residue_id(imol, residue_spec)
        # print "BL DEBUG:: --------------------- add_cho_restraints_for_residue: glyco_tree_residue_id:", id
        # print "BL DEBUG:: --------------------- add_cho_restraints_for_residue: glyco_tree_residues:", coot.glyco_tree_residues_py(imol, residue_spec)
        add_cho_restraints_for_residue_with_id(imol, residue_spec, id)

def add_cho_restraints_for_residue_with_id(imol, residue_spec, glyco_id):

    # hacky?! 
    def pad_name(str):
        l = len(str)
        ret = "    "
        if (l > 3):
            return str
        else:
            return " " + str.ljust(3)

    # return list [atom_1, atom_2, d, esd]
    # or False
    def line2extra_bond_restraint_spec(parent_residue_spec, line):
        global cho_bond_sigma_scale
        parts = line.split()
        if (len(parts) == 7):
            at_name_1 = pad_name(parts[0])
            at_name_2 = pad_name(parts[1])
            mean_str = parts[2]
            sd_str = parts[3]
            n_str = parts[4]
            d_str = parts[5]
            mod_sarle_str = parts[6]
            # could be done in one?!
            mean = float(mean_str)
            sd = float(sd_str)
            n = float(n_str)
            d = float(d_str)
            mod_sarle = float(mod_sarle_str)

            if not (n >= 20):
                return False
            else:
                if not (d < 4.3):
                    return False
                else:
                    return [[residue_spec_to_chain_id(residue_spec),
                             res_spec_utils.residue_spec_to_res_no(residue_spec),
                             coot_utils.residue_spec_to_ins_code(residue_spec),
                             at_name_1, ""],
                            [residue_spec_to_chain_id(parent_residue_spec),
                             res_spec_utils.residue_spec_to_res_no(parent_residue_spec),
                             coot_utils.residue_spec_to_ins_code(parent_residue_spec),
                             at_name_2, ""],
                            mean, cho_bond_sigma_scale]

    def glyco_id2level_number(glyco_id):
        return glyco_id[0]

    def glyco_id2prime_arm_sym(glyco_id):
        return glyco_id[1]

    def glyco_id2residue_type(glyco_id):
        return glyco_id[2]

    def glyco_id2link_type(glyco_id):
        return glyco_id[3]

    def glyco_id2parent_residue_type(glyco_id):
        return glyco_id[4]

    def glyco_id2residue_spec(glyco_id):
        return glyco_id[5]


    # main line
    #
    if isinstance(glyco_id, list):
        dir = get_glyco_model_dir()
        level = glyco_id2level_number(glyco_id)
        prime_arm_sym = glyco_id2prime_arm_sym(glyco_id)
        res_name = glyco_id2residue_type(glyco_id)
        link_type = glyco_id2link_type(glyco_id)
        parent_res_name = glyco_id2parent_residue_type(glyco_id)
        parent_res_spec = glyco_id2residue_spec(glyco_id)
        fn_file = "model-level-" + str(level) + "-" + res_name + "-" + \
                  link_type + "-" + parent_res_name + ".tab"
        model_fn = os.path.join(dir, fn_file)
        if not os.path.isfile(model_fn):
            print("model file does not exist", model_fn)
        else:
            try:
                f = open(filename, 'r')
                lines = f.readlines()
                f.close()
            except:
                print("BL WARNING: couldnt open file ", filename)
                return False # with something?
            print("INFO:: read %s lines from file %s" %(lines, model_fn))
            new_restraints = [line2extra_bond_restraint_spec(line) for line in lines]
            # print "BL DEBUG:: ", new_restraints
            coot.add_extra_bond_restraints_py(imol, new_restraints)


def test_get_cho_restraints(imol):
    raw_carbo_tree_list = []

    for chain_id in coot_utils.chain_ids(imol):
        for res_serial in range(coot.chain_n_residues(chain_id, imol)):
            res_no = coot.seqnum_from_serial_number(imol, chain_id, res_serial)
            rn = coot.residue_name(imol, chain_id, res_no, "")
            if isinstance(rn, str):
                if (rn == "NAG"):
                    residue_spec = [chain_id, res_no, ""]
                    rl = coot.glyco_tree_residues_py(imol, residue_spec)
                    print("rl:", rl)
                    coot.print_glyco_tree(imol, chain_id, res_no, "")
                    if not isinstance(rl, list):
                        print("bad glyco-tree for residue", residue_spec)
                    else:
                        if (len(rl) > 3):
                            id = glyco_tree_residue_id(imol, residue_spec)
                            print("residue %s residue id %s" %(residue_spec, id))
                            add_cho_restraints_for_residue_with_id(imol,
                                                                   residue_spec,
                                                                   id)
                        else:
                            print("BL WARNING:: rl <=3", rl)

def correlation_coefficient_of_this_tree():
    with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                   aa_ins_code, aa_atom_name, aa_alt_conf,
                                   aa_res_spec]:
        residues = glyco_tree_residues(aa_imol, aa_res_spec)
        cc = coot.map_to_model_correlation(aa_imol, residues, [], 0,
                                      coot.imol_refinement_map())
        txt = "residues %s\ncc: %5.3f" %(residues, cc)
        print(txt)
        coot.info_dialog(txt)
