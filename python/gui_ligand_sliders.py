global cache_ligand_metrics
cache_ligand_metrics = False

# return 3 things
#
def mtz_file_name2refinement_column_labels(file_name):

    columns = coot.get_mtz_columns(file_name)
    read_success = columns.read_success

    if (not read_success == 1):

        # unhappy path
        print("Failed to read columms from file %s for map molecule %s" \
              %(file_name, imol_map()))

    else:
        
        # happy path
        f_cols = columns.f_cols
        sigf_cols = columns.sigf_cols
        rfree_cols = columns.r_free_cols

        l1 = f_cols.size()
        l2 = sigf_cols.size()
        l3 = rfree_cols.size()

        if (not all([x > 0 for x in [l1, l2 ,l3]])):

            # unhappy path
            print("Failed to find columns of the necessary types from %s: %s %s %s!" \
                  %(file_name, l1, l2, l3))
            return False, False, False
            
        else:

            # happy path
            f_col_label = f_cols[0].column_label
            sigf_col_label = sigf_cols[0].column_label
            r_free_col_label = rfree_cols[0].column_label

            return f_col_label, sigf_col_label, r_free_col_label


def ligand_validation_metrics_gui_list_wrapper_pre(
    imol, chain_id, res_no, ins_code, refmac_input_mtz_file_name,
    fp_col, sigfp_col, rfree_col, refmac_dir):

    global cache_ligand_metrics

    m = ligand_check.get_metrics_for_ligand(imol, chain_id, res_no, ins_code,
                               refmac_input_mtz_file_name,
                               fp_col, sigfp_col, rfree_col, refmac_dir)

    if cache_ligand_metrics:
        m = [0.935339450836182, 10.6290054321289,
             [0, 10, 18, 70, 143],
             [-0.0306930494488489, 0.386785203839611, 0.0146206432429434, 10865.0,
              5373.01106146351, -27.3189968762454, 0.0173681300165985, 0.0280039555722363,
              -8.86567375069092e-10, -0.0025144037621947, 0.117837198078632, 0.120915851909265],
             [0.952432048573266, 13.3150000572205, 13.9800000190735, 0.176323987538941]]

        if (not isinstance(m, list)):
            # unhappy path
            print("WARNING:: no ligand metrics found.")
        else:
            # happy path ;-)
            diff_d = 0.05
            d = m[0]
            mwz = m[1]           # Mogul Z
            contact_info = m[2]  # number of bad contacts
            bc = contact_info[0]
            low_is_good = 1
            high_id_good = 0

            # if mogul ran OK, then we can display the mogul markup
            #
            if (isinstance(mwz, list)):
                mogul_out_file_name = mwz[1]
                coot.mogul_markup(imol, chain_id, res_no, ins_code,
                             mogul_out_file_name)
                mwz = mwz[0]

            percentile_d = coot.get_ligand_percentile("density_correlation", d, high_id_good)
            percentile_diff_d = coot.get_ligand_percentile("coot_diff_map_correlation", diff_d, low_is_good)
            percentile_mwz = coot.get_ligand_percentile("mogul_z_worst", mwz, low_is_good)
            percentile_bc = coot.get_ligand_percentile("bumps_1", bc, low_is_good)

            if (percentile_d < 0):
                # just an example, but means we do not have ligands-2016.db
                txt = "BL INFO:: we dont have ligands-2016.db, so \n" + \
                      "percentiles are no available and graph meaningless!"
                coot.info_dialog(txt)
            input_to_sliders = [["Direct map density correl.", percentile_d, d],
                                [" Diff map density correl.", percentile_diff_d, diff_d],
                                ["            Mogul Z-worst", percentile_mwz, mwz],
                                ["             Bad contacts", percentile_bc, bc]]

            ligand_validation_sliders.ligand_validation_metrics_gui_list_wrapper(input_to_sliders)
                
            
        
        
#if (use_gui_qm != 2):
# disable for now since I cannot run mogul and test it in real...
if (0):
    menu = coot_gui.coot_menubar_menu("Ligand")

    def ligand_metric_slider_func():
        
        imol_map = coot.imol_refinement_map()

        if (not coot_utils.valid_map_molecule_qm(imol_map)):
            coot.add_status_bar_text("No valid refinement map molecule")
        else:
            l = refmac_parameters(imol_map)
            if l:
                refmac_input_mtz_file_name = l[0]
            else:
                refmac_input_coot.mtz_file_name = mtz_file_name(imol_map)
            refmac_dir = coot_utils.get_directory("coot-refmac")

            f_col_label, sigf_col_label, r_free_col_label = \
                         mtz_file_name2refinement_column_labels(refmac_input_mtz_file_name)
            print("    f_col_label:", f_col_label)
            print(" sigf-col-label:", sigf_col_label)
            print("rfree-col-label:", r_free_col_label)

            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                ligand_validation_sliders.ligand_validation_metrics_gui_list_wrapper_pre(
                    aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                    refmac_input_mtz_file_name,
                    f_col_label, sigf_col_label, r_free_col_label,
                    refmac_dir)

    add_simple_coot_menu_menuitem(
        menu, "Ligand Metric Sliders",
        lambda func: ligand_metric_slider_func())
        
            
        
