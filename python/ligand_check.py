mogul_results = mogul_results_py

def get_metrics_for_ligand(imol, chain_id, res_no, ins_code,
                           refmac_input_mtz_file_name,
                           fobs_col, sig_fobs_col, rfree_col,
                           refmac_dir):

    """
    Args:
        refmac_dir: the dir where input/output refmac files should be written.

    Return:
        [correlation, mogul, bumps, difference_map_stats, b_factor_info]
        where:
            mogul: is either a list [mogul_z_worst, mogul_out_file_name] or an
            error status
            difference_map_stats:
                the output of map_to_model_correlation_stats_py
                [correlation, var_x, var_y, n, sum_x, sum_y, D, D2 (based on mean/sd of the map at the ligand) 
                map_mean map_mean_at_ligand map_sd map_sd_at_ligand )
            b_factor_info:
                [median_ratio median_ligand mediand_-env ks_testd_-result]
    """

    def local_refmac(stub_name):

        ligand_spec = [chain_id, res_no, ins_code]
        neighbs = residues_near_residue(imol, ligand_spec, 4)

        rn = residue_name(imol, chain_id, res_no, ins_code)
        n_ligand_atoms = het_group_n_atoms(rn)

        if (not isNumber(n_ligand_atoms)):
            print "BL ERROR:: failed liagnd atoms not a number."
            return False
        else:
            refmac_out_sfs_file_name = os.path.join(refmac_dir, \
                                                    stub_name + "-with-ligand-refmac.mtz")
            with_ligand_pdb_file_name = os.path.join(refmac_dir, \
                                                     stub_name + "-with-ligand.pdb")

            make_directory_maybe(refmac_dir)
            make_directory_maybe("coot-refmac") # XYZOUT goes here
            write_pdb_file(imol, with_ligand_pdb_file_name)
            r = refmac_calc_sfs_make_mtz_with_columns(with_ligand_pdb_file_name,
                                                      refmac_input_mtz_file_name,
                                                      refmac_out_sfs_file_name,
                                                      fobs_col, sig_fobs_col, rfree_col)
            if not r:
                print "BL ERROR:: failed calculating sfs in refmac"
                return False
            else:
                # happy path
                return refmac_out_sfs_file_name

    # get_correlation at the ligand for the direct (FWT) map
    #
    def get_correlation(stub_name):

        global refmac_extra_params

        refmac_extra_params = ["refine exclude all from " + \
                               str(res_no) + " " + \
                               chain_id + " to " + \
                               str(res_no) + " " + \
                               chain_id]

        print "BL DEBUG:: in get_correlation(): refmac_extra_params:", refmac_extra_params

        ligand_spec = [chain_id, res_no, ins_code]
        refmac_out_sfs_file_name = local_refmac(stub_name)

        if not refmac_out_sfs_file_name:
            return False
        else:
            # happy path
            imol_map = make_and_draw_map(refmac_out_sfs_file_name,
                                         "FWT", "PHWT", "", 0, 0)
            neighbs = residues_near_residue(imol, ligand_spec, 4)

            c = map_to_model_correlation(imol, [ligand_spec],
                                         neighbs, 0, imol_map)
            close_molecule(imol_map)
            return c

    # return False or a list of stats.
    #
    def get_ligand_difference_map_stats(stub_name):

        global refmac_extra_params

        refmac_extra_params = False # reset - no ligand exclusion

        refmac_out_sfs_file_name = local_refmac(stub_name + "-for-ligand-diff-map")
        ligand_spec = [chain_id, res_no, ins_code]
        neighbs = residues_near_residue(imol, ligand_spec, 4)

        if not refmac_out_sfs_file_name:
            return False
        else:
            # happy path
            imol_map = make_and_draw_map(refmac_out_sfs_file_name,
                                         "DELFWT", "PHDELWT", "", 0, 1)
            # now do some stats on the map at the ligand site

            c = map_to_model_correlation_stats_py(imol, [ligand_spec],
                                                  neighbs, 10, imol_map)
            print "BL INFO:: residue %s density statistics %s!" \
                  %(ligand_spec, c)
            return c

    # Return an error status (False, i.e. not a list) or a list
    # [mogul_z_worst, mogul_out_file_name] on success.
    #
    # (PE note: we want the filename so that we can show mogul markup
    # on the ligand when we activate this function from the gui)
    def get_mogul_score(use_cache_qm):

        # if run_result is a list it contains the mogul-output file name
        # if it is False something went wrong
        #
        run_result = run_mogul("BONDS_AND_ANGLES", imol, chain_id,
                               res_no, ins_code,
                               "ligand-check", use_cache_qm)

        print "BL DEBUG:: run_results (mogul): ", chain_id, res_no, run_result

        if not run_result:
            return False
        else:
            # happy path
            if use_cache_qm:
                # return a random mogul_score
                return [123, run_result]
            mogul_results_list = mogul_results(run_result)
            if (not isinstance(mogul_results_list, list)):
                return False
            else:
                if (len(mogul_results_list) == 0):
                    mogul_score = "MOGUL_NO_STATS"
                else:
                    mogul_score = max(mogul_results_list)
                return [mogul_score, run_result]

    # return a list: n_bad_overlaps n_hydrogen_bonds n_small_overlaps n_close_contacts n_wide_contacts
    #
    def get_bump_score():
        ligand_spec = [chain_id, res_no, ins_code]
        cs = contact_score_ligand(imol, ligand_spec)
        graphics_draw()
        return cs

    # return a list of [median_ratio, median_ligand, mediand_-env, ks_test_result]
    # stub_name is only passed so that we can write diagnostics to standard-out.
    #
    def get_b_factor_distribution_metrics(stub_name):

        def median_ratio(ligand_b_factors, env_b_factors):
            # maybe the second one should include the ligand Bs?!
            return median(ligand_b_factors) / median(env_b_factors)

        def filter_out_waters(imol, env_residues):
            def is_not_water(residue_item):
                rn = residue_name(imol,
                                  residue_spec_to_chain_id(residue_item),
                                  residue_spec_to_res_no(residue_item),
                                  residue_spec_to_ins_code(residue_item))
                return (rn != "HOH" and rn != "WAT")
            return [res_item for res_item in env_residues if is_not_water(res_item)]

        # from http://stackoverflow.com/questions/24101524/finding-median-of-list-in-python
        def median(lst):
            lst = sorted(lst)
            if len(lst) < 1:
                return None
            if len(lst) %2 == 1:
                return lst[((len(lst)+1)//2)-1]
            else:
                return float(sum(lst[(len(lst)//2)-1:(len(lst)//2)+1]))/2.0

        # Return a list of length 2: 
        # The first item is a list of atoms for the residue specified by ligand-spec
        # It may be False, in which case this function failed to return a result
        # If it is not False, the 2nd is a list of atom b-factors
        #
        def ligand_environment_temperature_factors(imol, ligand_spec, radius):
            atoms = residue_info(imol,
                                 residue_spec_to_chain_id(ligand_spec),
                                 residue_spec_to_res_no(ligand_spec),
                                 residue_spec_to_ins_code(ligand_spec))
            env_residues = residues_near_residue(imol, ligand_spec, radius)
            non_water_env_residues = filter_out_waters(imol, env_residues)
            env_residues = map(lambda res_spec: residue_info(imol,
                                                             residue_spec_to_chain_id(res_spec),
                                                             residue_spec_to_res_no(res_spec),
                                                             residue_spec_to_ins_code(res_spec)),
                               non_water_env_residues)
            # this is a list of residue info not atoms, so flatten
            env_atoms = []
            map(env_atoms.extend, env_residues)
            if isinstance(atoms, list):
                r1 = map(lambda atom: atom[1][1][0] \
                         if isinstance(atom[1][1], list) \
                         else atom[1][1], atoms)
            else:
                r1 = False
            if env_atoms:
                r2 = map(lambda atom: atom[1][1][0] \
                         if isinstance(atom[1][1], list) \
                         else atom[1][1], env_atoms)
            else:
                r2 = False

            return [r1, r2]

        # main line of get_b_factor_distribution_metrics
        #
        ligand_spec = [chain_id, res_no, ins_code]
        lig_env_temp_factors = ligand_environment_temperature_factors(imol, ligand_spec, 5)

        v1 = lig_env_temp_factors[0]
        v2 = lig_env_temp_factors[1]
        print "b-factor kolmogorov-smirnov lig:", stub_name, ligand_spec, v1
        print "b-factor kolmogorov-smirnov env:", stub_name, ligand_spec, v2
        temp_factor_median_ratio = median_ratio(v1, v2)
        kolmogorov_smirnov_result = kolmogorov_smirnov(v1, v2)

        return [temp_factor_median_ratio, median(v1), median(v2),
            kolmogorov_smirnov]

    # main line of get_metrics_for_ligand
    #
    stub_name = molecule_name_stub(imol, 0)

    b_factor_info = get_b_factor_distribution_metrics(stub_name)

    # add error checking to this
    #
    cor = get_correlation(stub_name)
    if (isNumber(cor)):
        dms = get_ligand_difference_map_stats(stub_name)
        if (not isinstance(dms, list)):
            return False # error
        else:
            mog = get_mogul_score(False)  # use cache for testing only
            if mog:
                bmp = get_bump_score()
                if (isinstance(bmp, list)):
                    return [cor, mog, bmp, dms, b_factor_info]
                else:
                    return False
            else:
                return False
    else:
        return False
            
            

def run_mogul(mode, imol, chain_id, res_no, ins_code, prefix_str, use_cache_qm):
    # dummy since I cannot test mogul
    # will return bogus info for now.
    if use_cache_qm:
        # for testing return a random Z and random file name
        return [2, "dummy"]
    else:
        return False
    
