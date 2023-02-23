# refmac.py 
#
# Copyright 2005, 2006 by Bernhard Lohkamp
# Copyright (C) 2008 by Bernhard Lohkamp, The University of York
# Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Extra parameters can be passed to refmac using either a file
# "refmac-extra-params" or by setting the variable
# refmac-extra-params.  The LABIN line should not be part those
# extra parameters, of course - Coot takes care of that.
# refmac_extra_params should be a list, like
# refmac_extra_params = ['NCYC 0', 'WEIGHT 0.2']
#
# if refmac_extra_params is given we use that! Otherwise we check for a 
# refmac-extra-params.txt file
# also params given in .coot.py are considered very first, so check is:
#
# 1.) whatever is written here
#
# 2.) .coot.py (needs to include: global refmac_extra_params)
#
# 3.) refmac-extra-params.txt
#
# deftexi refmac_extra_params

import os
import sys
import numbers
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import GLib
import coot
import coot_utils

global refmac_extra_params
# refmac_extra_params has to be a list
# this is an example, feel free to change it (uncomment the second line
# and comment the next one out)!!
refmac_extra_params = None
#refmac_extra_params = ['NCYC 2', 'WEIGHT 0.3']

# If phase-combine-flag is 1 then we should do a phase combination
# with the coefficients that were used to make the map, specifically
# PHIB and FOM, that are the result of (say) a DM run.  It is best
# of course to use HL coefficients, if DM (say) wrote them out.
# 
# If phase-combine-flag is 0, then phib-fom-pair is ignored.  If
# phase-combine-flag is 1 phib-fom-pair is presumed to be an
# improper list of the PHIB and FOM column labels (as strings).
#
# Is there a case where you want phase recombination, but do not
# have FOMs?  Don't know, I presume not.  So we required that if
# phase-combine-flag is 1, then PHIB and FOM are present.  Else we
# don't run refmac.
#
# If ccp4i project dir is "" then we ignore it and write the refmac
# log file here.  If it is set to something, it is the prefix of the
# refmac log file.  In that case it should end in a "/".  This is
# not tested for here!  The pdb-in-filename etc are dealt with in
# the c++ execute_refmac function.
#
# BL says: all the phase stuff not properly tested, yet!
#
# if swap_map_colours_post_refmac? is not 1 then imol-mtz-molecule
# is ignored.
#


#BL says: we need to have some global parameter definitions here (for now)
global refmac_count
refmac_count = 0


global use_gui_qm

# /a/b/c -> c
def split_label(label):
    return label[label.rfind("/")+1:len(label)]

# BL says: for loggraph we use a function
# Modernisation: Actually we first try to run qtrview, but if this is not
# available we run loggraph.
# 
def run_loggraph(logfile):
    
    import os

    # Modern: try qtrview first
    if coot_utils.command_in_path_qm("logview"):
        log_view = coot_utils.find_exe("logview")
        coot_utils.run_concurrently(log_view, [logfile])
    else:
        # no qtrtview (which can be run easily on the command line

        # this is the easy code, with no exception handling, so commented out (next 3 lines):
        # loggraph_exe = os.path.join(os.environ['CCP4I_TOP'],'bin/loggraph.tcl')
        # bltwish_exe = os.path.join(os.environ['CCP4I_TCLTK'],'bltwish')
        # os.spawnl(os.P_NOWAIT, bltwish_exe , bltwish_exe , loggraph_exe , logfile)

        # now lets check properly if the executables are there:
        # 1st check TCL
        ccp4i_tcltk = os.getenv('CCP4I_TCLTK')
        if not ccp4i_tcltk:
            print("BL ERROR:: We cannot allocate $CCP4I_TCLTK so no wish available")
        else:
            wish_command = "wish"
            if coot_utils.is_windows():
                # Windows
                wish_command += '.exe'
            wish_exe = os.path.join(ccp4i_tcltk, wish_command)
            if (not os.path.isfile(wish_exe)):
                # some windows machines have CCP4I_TCLTK not including bin
                wish_exe_win_alt = os.path.join(ccp4i_tcltk,
                                                "bin",
                                                wish_command)
                if (os.path.isfile(wish_exe_win_alt)):
                    wish_exe = wish_exe_win_alt
                # maybe we want to check for bltwish instead then?!

            if (not os.path.isfile(wish_exe)):
                print("BL ERROR:: We have $CCP4I_TCLTK but we cannot find wish")
            else:
                # now that we have tcl + wish we check for loggraph.tcl
                ccp4i_top = os.getenv('CCP4I_TOP')
                if not ccp4i_top:
                    print("BL ERROR:: We cannot allocate $CCP4I_TOP so no loggraph available")
                else:
                    loggraph_exe = os.path.join(ccp4i_top, "bin", "loggraph.tcl")
                    if not os.path.isfile(loggraph_exe):
                        print("BL ERROR:: We have $CCP4I_TOP but we cannot find loggraph.tcl")
                    else:
                        # print 'BL DEBUG:: We have allocated everything to run loggraph and shall do that now'
                        #os.spawnl(os.P_NOWAIT, bltwish_exe , bltwish_exe , loggraph_exe , logfile)
                        coot_utils.run_concurrently(wish_exe, [loggraph_exe, logfile])


# make_molecules_flag is synonymous with continue after refmac run, i.e.
# read molecules, run loggraph etc., furthermore not threaded, in other words
#
# make_molecules_flag is tested for being = 0, if not 0, then this is
# the main thread and we can do graphics things, like read in a pdb
# and mtz file to make new molecules.
#
# BL says: I believe in python we have everything handled already correctly,
# so this is just a different name for now with not invoked functions.
# Not sure why we want to block the graphics? Maybe needs a new keywd arg
# to run in thread or not!?
#
def run_refmac_by_filename(pdb_in_filename, pdb_out_filename,
                           mtz_in_filename, mtz_out_filename,
                           extra_cif_lib_filename, imol_refmac_count,
                           swap_map_colours_post_refmac_p,
                           imol_mtz_molecule, show_diff_map_flag,
                           phase_combine_flag, phib_fom_pair,
                           force_n_cycles, make_molecules_flag,
                           ccp4i_project_dir, f_col, sig_f_col, r_free_col=""):

    import os
    ref_fin_fn = ".refmac-is-finished"

    def refmac_is_finished_qm():
        return os.path.isfile(ref_fin_fn)

    def refmac_run_ok_qm():
        if os.path.isfile(ref_fin_fn):
            with open(ref_fin_fn, 'r') as fn:
                s = fn.readline()
                return s == "status 0"
        else:
            return False

    # needs to return true (for keep going) or false
    def check_for_refmac_finished_then_do_stuff():
        if refmac_is_finished_qm():
            if refmac_run_ok_qm():
                set_recentre_on_read_pdb(0)
                imol = read_pdb(pdb_out_filename)
                new_map_id = make_and_draw_map(mtz_out_filename,
                                               "FWT", "PHWT", "", 0, 0)
                if (swap_map_colours_post_refmac_p == 1):
                    swap_map_colours(imol_mtz_molecule, new_map_id)
                if (show_diff_map_flag == 1):
                    make_and_draw_map(mtz_out_filename,
                                               "DELFWT", "PHDELWT", "", 0, 1)
                    # do we have anom?
                    try:
                        make_and_draw_map(mtz_out_filename, "FAN", "PHAN", "", 0, 1)
                    except:
                        print("WARNING:: no anom map")
            return False # stop
        return True # continue

    if os.path.isfile(ref_fin_fn):
        os.remove(ref_fin_fn)

 #   run_python_thread(run_refmac_by_filename_inner,
    run_refmac_by_filename_inner(
        pdb_in_filename, pdb_out_filename,
        mtz_in_filename, mtz_out_filename,
        extra_cif_lib_filename, imol_refmac_count,
        swap_map_colours_post_refmac_p,
        imol_mtz_molecule, show_diff_map_flag,
        phase_combine_flag, phib_fom_pair,
        force_n_cycles, make_molecules_flag,
        ccp4i_project_dir, f_col, sig_f_col, r_free_col="")

    # gobject.timeout_add(1000, check_for_refmac_finished_then_do_stuff)


def run_refmac_by_filename_inner(pdb_in_filename, pdb_out_filename,
                           mtz_in_filename, mtz_out_filename,
                           extra_cif_lib_filename, imol_refmac_count,
                           swap_map_colours_post_refmac_p,
                           imol_mtz_molecule, show_diff_map_flag,
                           phase_combine_flag, phib_fom_pair,
                           force_n_cycles, make_molecules_flag,
                           ccp4i_project_dir, f_col, sig_f_col, r_free_col=""):

    # Paul's scheme code is ommitting threaded print. Why? FIXME

    global refmac_count
    import os, stat, operator

    # first check if refmac exists?
    refmac_execfile = coot_utils.find_exe("refmac5", "CBIN", "CCP4_BIN", "PATH")
    if not refmac_execfile:
        print("WARNING:: run_refmac_by_filename_inner(): no refmac found")
        print("  - no new map and molecule available")
        return
    else:
        # we have refmac
        # False here is in a sub-thread, do it noiselessly; True as normal.
        to_screen_flag = False if make_molecules_flag == 0 else True

    # some additional argument jiggery-pokery: convert
    # [["/crystal/thing/R-free"]] to ["/crystal/thing/R-free"]
    if isinstance(r_free_col, list):
        if r_free_col:
            if (isinstance(r_free_col[0], list)):
                r_free_col = r_free_col[0]

    local_r_free_col = [] if not r_free_col else r_free_col[0]
    # need to check for f-col being a string or list
    labin_string = ""
    if not ((isinstance(f_col, str)) and (f_col == "")):
        if (phase_combine_flag == 3 and (len(f_col) == 2)):
            labin_string = "LABIN F+=" + split_label(f_col[0]) + \
                " SIGF+=" + split_label(sig_f_col[0]) +\
                " F-=" + split_label(f_col[1]) +\
                " SIGF-=" + split_label(sig_f_col[1])
        else:
            if (refmac_use_intensities_state()):
                labin_string = "LABIN IP=" + split_label(f_col) + \
                    " SIGIP=" + split_label(sig_f_col)
            else:
                labin_string = "LABIN FP=" + split_label(f_col) +\
                    " SIGFP=" + split_label(sig_f_col)

        if (local_r_free_col != ""):
            labin_string += " FREE=" + split_label(r_free_col)

        if (phase_combine_flag == 1):
            # we have Phi Fom pair
            if (phib_fom_pair[0] != "" and phib_fom_pair[1] !=0):
                labin_string += " - \nPHIB=" + split_label(phib_fom_pair[0]) + \
                    " FOM=" + split_label(phib_fom_pair[1])
        if (phase_combine_flag == 2):
            # we have HLs
            if (phib_fom_pair[1] == ""):
                hl_list = eval(phib_fom_pair[0])
                if (len(hl_list) == 4):
                    hl_label_ls = ["HLA", "HLB", "HLC", "HLD"]
                    labin_string += " - \n"
                    for i in range(4):
                        # shorten the label
                        hl_label = split_label(hl_list[i])
                        labin_string += " " + hl_label_ls[i] + "=" + hl_label

# BL says: command line args have to be string not list here
    command_line_args = ["XYZIN",  pdb_in_filename,
                         "XYZOUT", pdb_out_filename,
                         "HKLIN",  mtz_in_filename,
                         "HKLOUT", mtz_out_filename]

    if (extra_cif_lib_filename != ""):
        command_line_args += ["LIBIN", extra_cif_lib_filename]

    std_lines = ["MAKE HYDROGENS NO"] # Garib's suggestion 8 Sept 2003

    refinement_type = coot.get_refmac_refinement_method()
    if (refinement_type == 1):
        # do rigid body refinement
        std_lines.append("REFInement TYPE RIGID")

    imol_coords = coot.refmac_imol_coords()
    if isinstance(force_n_cycles, numbers.Number):
       if force_n_cycles >=0:
           if (refinement_type == 1):
               std_lines.append("RIGIDbody NCYCle " + str(force_n_cycles))
               if coot.get_refmac_refinement_method() == 1:
                   group_no = 1
                   if (coot.is_valid_model_molecule(imol_coords)):
                       for chain_id in coot_utils.chain_ids(imol_coords):
                           if (not coot_utils.is_solvent_chain_qm(imol_coords, chain_id)):
                               n_residues = coot.chain_n_residues(chain_id, imol_coords)
                               start_resno = coot.seqnum_from_serial_number(imol_coords ,chain_id, 0)
                               end_resno   = coot.seqnum_from_serial_number(imol_coords ,chain_id, n_residues-1)
                               rigid_line = "RIGIdbody GROUp " +  str(group_no) \
                               + " FROM " + str(start_resno) + " " + chain_id \
                               + " TO "   + str(end_resno)   + " " + chain_id 
                               std_lines.append(rigid_line)
                               group_no += 1
                   else:
                       print("WARNING:: no valid refmac imol!")
           else:
               std_lines.append("NCYC " + str(force_n_cycles))
    # TLS?
    if (refinement_type == 2):
        tls_string = "REFI TLSC 5"
        std_lines.append(tls_string)

    # TWIN?
    if False: # coot.refmac_use_twin_state():
        std_lines.append("TWIN")

    # SAD?
    if (phase_combine_flag == 3) and (labin_string == ""):
        std_lines.append("REFI SAD")

    std_lines.extend(refmac_sad_params())

    # NCS?
    if False: # coot.refmac_use_ncs_state():
        if (get_refmac_version()[1] >= 5):
            std_lines.append("NCSR LOCAL")
        else:
            coot_utils.chain_ids_from_ncs = ncs_chain_ids(imol_coords)
            if coot_utils.chain_ids_from_ncs:
                for ncs_set in coot_utils.chain_ids_from_ncs:
                    no_ncs_chains = len(ncs_set)
                    if (no_ncs_chains > 1):
                        ncs_string = "NCSRestraints NCHAins " + str(no_ncs_chains) + " CHAIns "
                        for chain in ncs_set:
                            ncs_string += " " + chain
                        std_lines.append(ncs_string)

    data_lines = add_refmac_extra_params(std_lines, force_n_cycles)

    if (not extra_params_include_weight_p(data_lines)) :
        data_lines.append("WEIGHT AUTO 5")

    data_lines.append(labin_string)

    log_file_name_disambiguator = coot_utils.strip_path(coot_utils.file_name_sans_extension(pdb_in_filename))
    # this should be a database filename:
    refmac_log_file_name = os.path.join(coot_utils.get_directory("coot-refmac"),
                                        (ccp4i_project_dir if (len(ccp4i_project_dir) > 0) else "") + \
                                        "refmac-from-coot-" + \
                                        log_file_name_disambiguator + \
                                        "-" + \
                                        str(refmac_count) + ".log")

    refmac_count = imol_refmac_count + refmac_count + 1

    print("INFO: running refmac with these command line args: ", command_line_args)
    print("INFO: running refmac with these data lines: ", data_lines)
    try:
        print("environment variable:  SYMOP: ", os.environ['SYMINFO'])
    except:
        print(" not set !")
    try:
        print("environment variable: ATOMSF: ", os.environ['ATOMSF'])
    except:
        print(" not set !")
    try:
        print("environment variable:  CLIBD: ", os.environ['CLIBD'])
    except:
        print(" not set !")
    try:
        print("environment variable:   CLIB: ", os.environ['CLIB'])
    except:
        print(" not set !")

    data_lines += ["END"]

    if coot_utils.coot_has_gobject() and sys.version_info >= (2, 4) and make_molecules_flag:

        # can spawn refmac and add button

        print("DEBUG:: run_refmac_by_filename_inner(): calling run_concurrently with ", refmac_execfile, command_line_args, data_lines, refmac_log_file_name, to_screen_flag)

        run_concurrently_results = coot_utils.run_concurrently(refmac_execfile,
                                                               command_line_args,
                                                               data_lines,
                                                               refmac_log_file_name,
                                                               to_screen_flag)
        print("DEBUG:: run_refmac_by_filename_inner(): run_concurrently_results", run_concurrently_results)
        if run_concurrently_results == False:
            print("DEBUG:: run_refmac_by_filename_inner(): failed to start running concurrently")
        else:
            global use_gui_qm
            refmac_process, logObj = run_concurrently_results
            print("HHHHHHHHHHHHHHHHHHHHHHHHere 0000000000000")
            if coot_utils.use_gui_qm:
                print("HHHHHHHHHHHHHHHHHHHHHHHHere  PATH A")
                separator   = add_coot_toolbar_separator()
                kill_button = coot_toolbar_button("Kill refmac",
                                                  "kill_process(" + str(refmac_process.pid)+ ")",
                                                  "stop.svg")
                button_tup = (kill_button, separator)
                add_status_bar_text("Running refmac")
            else:
                button_tup = None
                print("HHHHHHHHHHHHHHHHHHHHHHHHere  PATH B")

            print("HHHHHHHHHHHHHHHHHHHHHHHHere  PATH C")
            GLib.timeout_add(1000,
                             post_run_refmac,
                             imol_refmac_count,
                             imol_mtz_molecule, swap_map_colours_post_refmac_p,
                             show_diff_map_flag,
                             pdb_out_filename, mtz_out_filename, mtz_in_filename,
                             refmac_log_file_name,
                             phib_fom_pair, f_col, sig_f_col, r_free_col,
                             phase_combine_flag, make_molecules_flag,
                             refmac_process, logObj,
                             button_tup, True)
            print("HHHHHHHHHHHHHHHHHHHHHHHHere  PATH D")
    else:
        print("HHHHHHHHHHHHHHHHHHHHHHHHere  PATH OTHER")
        # no gobject and no subprocess, so run 'old' coot_utils.popen_command
        refmac_status = coot_utils.popen_command(refmac_execfile,
                                      command_line_args,
                                      data_lines,
                                      refmac_log_file_name, to_screen_flag)
        if make_molecules_flag == 0:
            # e.g. from get_recent_pdb, i.e. dont load new files, basically
            # dont thread here...
            return pdb_out_filename, mtz_out_filename  # return values
        else:
            # pass refmac_status as refmac_process!?
            post_run_refmac(imol_refmac_count,
                            imol_mtz_molecule, swap_map_colours_post_refmac_p,
                            show_diff_map_flag,
                            pdb_out_filename, mtz_out_filename, mtz_in_filename,
                            refmac_log_file_name,
                            phib_fom_pair, f_col, sig_f_col, r_free_col,
                            phase_combine_flag, make_molecules_flag,
                            refmac_status)


def post_run_refmac(imol_refmac_count,
                    imol_mtz_molecule, swap_map_colours_post_refmac_p,
                    show_diff_map_flag,
                    pdb_out_filename, mtz_out_filename, mtz_in_filename,
                    refmac_log_file_name,
                    phib_fom_pair, f_col, sig_f_col, r_free_col,
                    phase_combine_flag, make_molecules_flag,
                    refmac_process, logObj=None,
                    button = None, run_in_timer=False):

    if (run_in_timer):
        # we run refmac concurrently
        if  (refmac_process.poll() != None):
            # refmac is finished
            logObj.close()
            refmac_status = refmac_process.poll()
            if button:
                button[0].destroy()
                button[1].destroy()
            # and stop the loop - later
        else:
            # continue the loop
            return True
    else:
        # refmac run via coot_utils.popen_command
        refmac_status = refmac_process

    if (refmac_status) : # refmac ran fail...
        print("Refmac Failed with status", refmac_status)
        add_status_bar_text("Refmac failed")
        if (button):
            button[0].destroy()
            button[1].destroy()
        if (run_in_timer):
            # stop the gobject timer
            print("... or was killed")
            return False

    else : # refmac ran OK.

        # BL says:
        # popen will always write a log file, so we check the size, is 0 if failed
        #refmac_log_size = os.stat(refmac_log_file_name)[stat.ST_SIZE]
        add_status_bar_text("Refmac finished")

        # only run loggraph if refinement
        if make_molecules_flag:
            run_loggraph(refmac_log_file_name)

        # deal with R-free column
        if (r_free_col == ""):
            r_free_bit = ["",0]
        else:
            r_free_bit = [r_free_col,1]

        recentre_status = recentre_on_read_pdb()
        novalue = set_recentre_on_read_pdb(0)
        imol = read_pdb(pdb_out_filename)

        if recentre_status:
            set_recentre_on_read_pdb(1)
        else:
            set_recentre_on_read_pdb(0)
        set_refmac_counter(imol, imol_refmac_count + 1)

        if (phase_combine_flag == 3 or coot.refmac_use_twin_state()):
            # for now we have to assume the 'standard' name, let's see if we can do better...
                # this may not be necessary!?
            args = [mtz_out_filename, "FWT", "PHWT", "", 0, 0, 1, "FP", "SIGFP"] + r_free_bit
        else:
            args = [mtz_out_filename, "FWT", "PHWT", "", 0, 0, 1, f_col, sig_f_col] + r_free_bit
        new_map_id = coot.make_and_draw_map_with_refmac_params(*args)

        # set the new map as refinement map
        if coot_utils.valid_map_molecule_qm(new_map_id):
            coot.set_imol_refinement_map(new_map_id)

        # store the refmac parameters to the new map (if we used a file)
        if (coot.get_refmac_used_mtz_file_state()):
            coot.set_stored_refmac_file_mtz_filename(new_map_id, mtz_in_filename)
            if (phase_combine_flag > 0 and phase_combine_flag < 3):
                # save the phase information
                phib = ""
                fom  = ""
                hla  = ""
                hlb  = ""
                hlc  = ""
                hld  = ""
                if (phase_combine_flag == 1):
                    # we have Phi Fom pair
                    phib = phib_fom_pair[0]
                    fom  = phib_fom_pair[1]
                if (phase_combine_flag == 2):
                    # we have HLs
                    hla, hlab, hlc, hld = eval(phib_fom_pair[0])
                coot.save_refmac_phase_params_to_map(new_map_id,
                                                phib, fom,
                                                hla, hlb, hld, hlc)

        if (swap_map_colours_post_refmac_p == 1) :
            coot.swap_map_colours(imol_mtz_molecule, new_map_id)

        # difference map (flag was set):
        if (show_diff_map_flag == 1):
            if (phase_combine_flag == 3):
                # for now we have to assume the 'standard' name, let's see if we can do better...
                args = [mtz_out_filename, "DELFWT", "PHDELWT", "", 0, 1, 1, "FP", "SIGFP"] + r_free_bit
                coot.make_and_draw_map_with_refmac_params(*args)
                # make an anomalous difference map too?!
                args = [mtz_out_filename, "FAN", "PHAN", "", 0, 1, 1, "FP", "SIGFP"] + r_free_bit
                try:
                    coot.make_and_draw_map_with_refmac_params(*args)
                except:
                    print("BL INFO:: couldnt make the anomalous difference map.")
            else:
                args = [mtz_out_filename, "DELFWT", "PHDELWT", "", 0, 1, 1, f_col, sig_f_col] + r_free_bit
                coot.make_and_draw_map_with_refmac_params(*args)

        # finally run the refmac_log_reader to get interesting things,
        # if real refinement
        if make_molecules_flag:
            read_refmac_log(imol, refmac_log_file_name)

        if (run_in_timer):
            # stop the gobject timer
            return False


# Return True if the list of strings @var{params_list} contains a
# string beginning with "WEIGHT".  If not return False
#
def extra_params_include_weight_p(params_list):

   have_weight = ['WEI' in x.upper() for x in params_list]
   if (sum(have_weight) == 1):
       return True
   elif (sum(have_weight) == 0):
       return False
   else:
       print('BL WARNING:: This shouldn\'t happen, we have more than one weight defined. God knows what refmac will do now...')
       return False


# If refmac-extra-params is defined (as a list of strings), then
# add that, else read the file "refmac-extra-params.txt" and add to the
# refmac parameters
#
# Return a list a list of strings.
#
def add_refmac_extra_params(pre_lines, force_no_cycles):
  global refmac_extra_params

  post_lines = pre_lines
  extra_params = []
  if refmac_extra_params:
      # we ignore a file
      extra_params = refmac_extra_params
  else:
      # see if there is a refmac parameter file
      extras_file_name = "refmac-extra-params.txt"
      if os.path.isfile(extras_file_name):
          try:
              f = open(extras_file_name,'r')
          except IOError:
              print('BL INFO:: we dont have refmac-extra-params file')
          except:
              print("BL ERROR:: unknown error reading", extras_file_name)
          else:
              extra_params = f.readlines()
              print('BL INFO:: we have refmac-extra-params file and read lines')
              f.close()
      # else
      # no extra params file, continue

  if extra_params:
      print("refmac extra params: ", extra_params)
      # remove line with NCYC when n_cycle is 0
      if (force_no_cycles == 0):
          post_lines = [x if not (x[0:3].upper() == 'NCY') else 'NCYC 0' for x in extra_params]
      else:
          # may override standard lines, so we need to check them all
          post_lines = extra_params
      for line in pre_lines:
          # check the first 3 chars of extra param as upper
          ls = [line[0:3].upper() in x.upper() for x in extra_params]
          if not sum(ls):
              # add the param
              post_lines.append(line.upper())

  return post_lines

def refmac_ncs_params():
    ret_ls = []
    if (refmac_use_ncs_state()):
        if (get_refmac_version()[1] >= 5):
            ret_ls = ["NCSR LOCAL"]
        else:
            chain_ids_from_ncs = ncs_chain_ids(imol_coords)
            if chain_ids_from_ncs:
                for ncs_set in chain_ids_from_ncs:
                    no_ncs_chains = len(ncs_set)
                    if (no_ncs_chains > 1):
                        ncs_string = "NCSRestraints NCHAins " + str(no_ncs_chains) + " CHAIns "
                        for chain in ncs_set:
                            ncs_string += " " + chain
                        ret_ls.append(ncs_string)
    return ret_ls

def refmac_sad_params():
    ret_ls = []
    if False: # refmac_use_sad_state() == 1:
        sad_atom_ls = get_refmac_sad_atom_info()
        for sad_atom in sad_atom_ls:
            sad_string = ""
            atom_name  = sad_atom[0]
            fp         = sad_atom[1]
            fpp        = sad_atom[2]
            wavelength = sad_atom[3]
            sad_string += "ANOM FORM " + atom_name
            if (fp):
                sad_string += " " + str(fp)
            if (fpp):
                sad_string += " " + str(fpp)
            if (wavelength):
                sad_string += " " + str(wavelength)
            ret_ls.append(sad_string)

    return ret_ls
# This is not run as a sub-thread now
#
# @return the output mtz file name or False
#
def run_refmac_for_phases(imol, mtz_file_name, f_col, sig_f_col):

    import os

    if os.path.isfile(mtz_file_name):
        if coot_utils.valid_model_molecule_qm(imol):
            dir_state = coot.make_directory_maybe("coot-refmac")
            if not dir_state == 0:
                print("Failed to make coot-refmac directory\n")
            else:
                stub = os.path.join(coot_refmac_dir, "refmac-for-phases")
                pdb_in = stub + ".pdb"
                pdb_out = stub + "-tmp.pdb"
                mtz_out = stub + ".mtz"
                cif_lib_filename = ""
                coot.write_pdb_file(imol, pdb_in)
                # preserve the ncs & tls state
                ncs_state = coot.refmac_use_ncs_state()
                tls_state = coot.refmac_use_tls_state()
                #print "BL DEBUG:: states were", ncs_state, tls_state
                coot.set_refmac_use_ncs(0)
                coot.set_refmac_use_tls(0)
                run_refmac_by_filename(pdb_in, pdb_out,
                        mtz_file_name, mtz_out,
                        cif_lib_filename, 0, 0, -1,
                        1, 0, [], 0, 1, "",
                        f_col, sig_f_col, "")
                # I wish I could reset the state, but I currently cannot
                # as refmac runs only after the button has been pushed
                # (and we wont be here any more)
                coot.set_refmac_use_ncs(ncs_state)
                coot.set_refmac_use_tls(tls_state)
                #print "BL DEBUG:: reset states", ncs_state, tls_state
                # check for success???? Just check for availability of mtz outfile
                # since we return this
                if os.path.isfile(mtz_out):
                    return mtz_out

        else:
            print("BL WARNING:: no valid model molecule!")
    else:
        print("BL WARNING:: mtzfile %s not found" %mtz_file_name)

    # something went wrong somewhere...
    return False


def refmac_for_phases_and_make_map(mtz_file_name, f_col, sig_f_col):

	import os, time

	if os.path.isfile(mtz_file_name):
		molecule_chooser_gui("  Choose a molecule from which to calculate Structure factors:  ",
			lambda imol: run_refmac_for_phases(imol, mtz_file_name, f_col, sig_f_col))


# This will read a refmac log file and extract information in
# coot_gui.interesting_things_gui format, i.e. a list with either coorinates or residues
# to be flagged as interesting
#
# the list of extracted information may look like:
# [[0, "A", 43,
#   [bond_deviations, [atom, value],[atom2, value2]],
#   [vdw_deviation, [atoms, value]],
#   [more deviations]],
#   next res] 
#
# e.g. coot_gui.interesting_things_gui(
#       "Bad things by Analysis X",
#       [["Bad Chiral",0,"A",23,"","CA","A"],
#        ["Bad Density Fit",0,"B",65,"","CA",""],
#        ["Interesting blob",45.6,46.7,87.5]]
#
def read_refmac_log(imol, refmac_log_file):
    
    import os
    import re

    # return the 'next' line containing the residue information,
    # otherwise False
    # try next 5 lines
    def get_next_line(pos):
        for j in range(1,6):
            next_line = lines[pos+j]
            if "ch:" in next_line:
                return next_line
        return False
        
    
    def get_warning(i):
        this_line = lines[i]
        next_line = get_next_line(i)
        if   (("CIS" in this_line) and next_line):
            # CIS peptide found
            # we ignore multiple conformations for now...
            # e.g.            ch:AA   res:  40  ARG      -->  41  GLU
            if not "ch:" in next_line:
                next_line = lines[i+2]
            item_ls   = split_clean(next_line)
            chain     = item_ls[1]
            chain_id  = chain[-1]
            res_no1   = int(item_ls[3])
            res_name1 = item_ls[4]
            res_no2   = int(item_ls[6])
            res_name2 = item_ls[7]
            angle_str = ""
            if (len(item_ls) > 8):
                angle_str = item_ls[7]
            else:
                item_ls = this_line.split()
                angle_str = item_ls[-1]

            info_text = "pre-WARNING: CIS bond: " + chain_id + " " + str(res_no1) + " - " \
                        + str(res_no2) + "  " + angle_str

            tooltip = " ".join(this_line.split()) + "\n" + " ".join(next_line.split())

            warning_info_list.append([info_text, imol, chain_id, res_no1, "", " CA ", "",
                                      ["dummy"], "", tooltip]) # first 2 are dummies for tooltip to work

        # oh dear, how do we handle this, at the moment its for
        # (all info in 1 line), e.g.:
        #   WARNING : conn: gap.     ch:AA res: 163  THR     --> 165  LYS
        #    "gap"-link will be created
        elif ("gap" in this_line):
            # gap-link
            chain     = this_line[30:32]       # this is e.g. "AA"
            chain_id  = chain[1]
            res_no1   = int(this_line[37:41])
            res_name1 = this_line[43:46]
            res_no2   = int(this_line[54:58])
            res_name2 = this_line[60:63]
            extra_info = ""
            tooltip = " ".join(this_line.split())
            if ("-link will" in next_line):
                # we have extra info
                extra_info = "'gap'-link created"
                tooltip += "\n"
                tooltip += " ".join(next_line.split())

            info_text = "pre-WARNING: connection gap: " + chain_id + " " \
                        + str(res_no1) + " - " + str(res_no2) + extra_info

            warning_info_list.append([info_text, imol, chain_id, res_no1, "", " CA ", "",
                                      ["dummy"], "", tooltip])

        elif (("big distance" in this_line) and next_line):
            # large distance
            # e.g. "            ch:AA   res:  71  ILE     -->  72  CYS      ideal_dist=     1.329"
            item_ls = split_clean(next_line)
            chain     = item_ls[1]
            chain_id  = chain[-1]
            res_no1   = int(item_ls[3])
            res_name1 = item_ls[4]
            res_no2   = int(item_ls[6]) # ?
            res_name2 = item_ls[7]
            dist_str  = item_ls[9]

            info_text = "pre-WARNING: large distance: " + chain_id + " " \
                        + str(res_no1) + " - " + str(res_no2) 
            tooltip = " ".join(this_line.split()) + "\n" + " ".join(next_line.split())

            warning_info_list.append([info_text, imol, chain_id, res_no1, "", " CA ", "",
                                      ["dummy"], "", tooltip])

        elif (("link" in this_line) and next_line):
            # link usually SS
            # e.g. "           ch:BB   res:   7  CYS      at:SG  .->BB   res:  96  CYS      at:SG  ."
            item_ls    = split_clean(next_line)
            chain1     = item_ls[1]
            chain_id1  = chain1[-1]
            res_no1    = int(item_ls[3])
            res_name1  = item_ls[4]
            atom_name1 = item_ls[6]
            alt_conf1  = item_ls[7]
            chain2     = item_ls[9]
            chain_id2  = chain2[-1]
            res_no2    = int(item_ls[11])
            res_name2  = item_ls[12]
            atom_name2 = item_ls[14]
            alt_conf2  = item_ls[15]

            link_type  = this_line[17:26].rstrip()

            if (alt_conf1 == "."):
                alt_conf1 = ""

            info_text = "pre-WARNING: " + link_type + " link found: " \
                        + chain_id1 + " " + str(res_no1) + " " + atom_name1 + " - " \
                        + chain_id2 + " " + str(res_no2) + " " + atom_name2

            tooltip = " ".join(this_line.split()) + "\n" + " ".join(next_line.split())

            warning_info_list.append([info_text, imol, chain_id1, res_no1, "", atom_name1, alt_conf1,
                                      ["dummy"], "", tooltip])
       
            
        # elif(): not sure what other warnings there could be...
        # FIXME maybe include tooltip with some fancyness

    def get_info(i):
        this_line = lines[i]
        next_line = get_next_line(i)
        if (("link is found" in this_line) and next_line):
            # link found but not used
            # we ignore multiple conformations for now...
            # e.g.             ch:AA   res:   2  VAL      at:CA  .->ch:AA   res:  89  PHE      at:CZ  .
            item_ls = split_clean(next_line)
            chain1     = item_ls[1]
            chain_id1  = chain1[-1]
            res_no1    = int(item_ls[3])
            res_name1  = item_ls[4]
            atom_name1 = item_ls[6]
            alt_conf1  = item_ls[7]
            chain2     = item_ls[10]
            chain_id2  = chain2[-1]
            res_no2    = int(item_ls[12])
            res_name2  = item_ls[13]
            atom_name2 = item_ls[15]
            alt_conf2  = item_ls[16]
            

            if (alt_conf1 == "."):
                alt_conf1 = ""

            info_text = "pre-INFO: not used link: " \
                        + chain_id1 + " " + str(res_no1) + " " + atom_name1 + " - " \
                        + chain_id2 + " " + str(res_no2) + " " + atom_name2

            tooltip = " ".join(this_line.split()) + "\n" + " ".join(next_line.split())

            warning_info_list.append([info_text, imol, chain_id1, res_no1, "", atom_name1, alt_conf1,
                                      ["dummy"], "", tooltip])
            
        # elif ("" in this_line): not sure what other INFO there may be

    def get_bond_dist(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 5):
            # e.g.
            # "A  15 ARG C   . - A  15 ARG O   . mod.= 1.295 id.= 1.231 dev= -0.064 sig.= 0.020"
            # new
            # "A    318 THR CA  A - A    318 THR C   . mod.= 1.278 id.= 1.525 dev=  0.247 sig.= 0.021"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[0]
            res_no1   = int(item_ls[1])
            res_name1 = item_ls[2]
            atom1     = item_ls[3]
            alt_conf1 = item_ls[4]
            if (item_ls[4] == '-'):
                alt_conf1 = "."
                item_ls.insert(4, alt_conf1)

            chain_id2 = item_ls[6]
            res_no2   = int(item_ls[7])
            res_name2 = item_ls[8]
            atom2     = item_ls[9]
            alt_conf2 = item_ls[10]
            if (len(alt_conf2) > 1):
                alt_conf2 = "."
                item_ls.insert(10, alt_conf2)
            
            mod       = float(item_ls[12])
            ideal     = float(item_ls[14])
            dev       = float(item_ls[16])
            sig       = float(item_ls[18])


            if (alt_conf1 == "."): alt_conf1 = ""
            if (alt_conf2 == "."): alt_conf2 = ""

            if coot_utils.debug():
                print("BL DEBUG:: haeve item lst", item_ls)

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Bond distance",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        mod, ideal, dev, sig, line]]])
               
            i += 1

    def get_bond_angle(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 5):
            # e.g. "A  50 GLU O   B - A  51 LEU N     mod.= 100.81 id.= 123.00 dev= 22.193 sig.=  1.600"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[0]
            res_no1   = int(item_ls[1])
            res_name1 = item_ls[2]
            atom1     = item_ls[3]
            alt_conf1 = item_ls[4]
            if (item_ls[4] == '-'):
                alt_conf1 = "."
                item_ls.insert(4, alt_conf1)

            chain_id2 = item_ls[6]
            res_no2   = int(item_ls[7])
            res_name2 = item_ls[8]
            atom2     = item_ls[9]
            alt_conf2 = item_ls[10]
            if (len(alt_conf2) > 1):
                alt_conf2 = "."
                item_ls.insert(10, alt_conf2)
            
            mod       = float(item_ls[12])
            ideal     = float(item_ls[14])
            dev       = float(item_ls[16])
            sig       = float(item_ls[18])

            if (alt_conf1 == " "): alt_conf1 = ""
            if (alt_conf2 == " "): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Bond angle",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        mod, ideal, dev, sig, line]]])
               
            i += 1

    def get_torsion(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 5):
            # e.g. "A  11 ASN CA    - A  12 LEU CA    mod.=-167.65 id.= 180.00 per.= 1 dev=-12.352 sig.=  3.000"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[0]
            res_no1   = int(item_ls[1])
            res_name1 = item_ls[2]
            atom1     = item_ls[3]
            alt_conf1 = item_ls[4]
            if (item_ls[4] == '-'):
                alt_conf1 = "."
                item_ls.insert(4, alt_conf1)

            chain_id2 = item_ls[6]
            res_no2   = int(item_ls[7])
            res_name2 = item_ls[8]
            atom2     = item_ls[9]
            alt_conf2 = item_ls[10]
            if (len(alt_conf2) > 1):
                alt_conf2 = "."
                item_ls.insert(10, alt_conf2)
            
            mod       = float(item_ls[12])
            ideal     = float(item_ls[14])
            period    = int(item_ls[16])
            dev       = float(item_ls[18])
            sig       = float(item_ls[20])

            if (alt_conf1 == " "): alt_conf1 = ""
            if (alt_conf2 == " "): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Torsion angle",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        mod, ideal, dev, sig, line]]])
               
            i += 1

    def get_chiral(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 5):
            # e.g. "A  51 LEU CG    mod.=   2.79 id.=  -2.59 dev= -5.375 sig.=  0.200"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[0]
            res_no1   = int(item_ls[1])
            res_name1 = item_ls[2]
            atom1     = item_ls[3]
            alt_conf1 = item_ls[4]
            if ("mod" in item_ls[4]):
                alt_conf1 = "."
                item_ls.insert(4, alt_conf1)
            
            mod       = float(item_ls[6])
            ideal     = float(item_ls[8])
            dev       = float(item_ls[10])
            sig       = float(item_ls[12])


            if (alt_conf1 == " "): alt_conf1 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Chiral",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        -999999, "", "", "",
                                                        mod, ideal, dev, sig, line]]])
               
            i += 1

    def get_planarity(i):
        i += 3      # jump to interesting lines (only 3 here....)
        while (len(lines[i]) > 5):
            # e.g. "Atom: A  59 ASP C   B deviation=   0.31 sigma.=   0.02"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[1]
            res_no1   = int(item_ls[2])
            res_name1 = item_ls[3]
            atom1     = item_ls[4]
            alt_conf1 = item_ls[5]
            if ("dev" in alt_conf1):
                alt_conf1 = "."
                item_ls.insert(5, alt_conf1)

            dev       = float(item_ls[7])
            sig       = float(item_ls[9])

            if (alt_conf1 == "."): alt_conf1 = ""       # guess

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Planarity",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        -999999, "", "", "",
                                                        -99999., -99999., dev, sig, line]]])       
               
            i += 1

    def get_vdw(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 5):
            # e.g. "A  26 CYS SG  A - A  75 ILE CD1 . mod.= 2.812 id.= 3.820 dev= -1.008 sig.= 0.300"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[0]
            res_no1   = int(item_ls[1])
            res_name1 = item_ls[2]
            atom1     = item_ls[3]
            alt_conf1 = item_ls[4]
            if (item_ls[4] == '-'):
                alt_conf1 = "."
                item_ls.insert(4, alt_conf1)

            chain_id2 = item_ls[6]
            res_no2   = int(item_ls[7])
            res_name2 = item_ls[8]
            atom2     = item_ls[9]
            alt_conf2 = item_ls[10]
            if (len(alt_conf2) > 1):
                alt_conf2 = "."
                item_ls.insert(10, alt_conf2)
            
            mod       = float(item_ls[12])
            ideal     = float(item_ls[14])
            dev       = float(item_ls[16])
            sig       = float(item_ls[18])

            if (alt_conf1 == "."): alt_conf1 = ""
            if (alt_conf2 == "."): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["VDW",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        mod, ideal, dev, sig, line]]])
               
            i += 1

    def get_b_value(i):
        i += 3      # jump to interesting lines
        while (len(lines[i]) > 5):
            # e.g. "B   5 PHE N     - B   4 GLN C        ABS(DELTA)= 15.990   Sigma=  1.500"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[0]
            res_no1   = int(item_ls[1])
            res_name1 = item_ls[2]
            atom1     = item_ls[3]
            alt_conf1 = item_ls[4]
            if (item_ls[4] == '-'):
                alt_conf1 = "."
                item_ls.insert(4, alt_conf1)

            chain_id2 = item_ls[6]
            res_no2   = int(item_ls[7])
            res_name2 = item_ls[8]
            atom2     = item_ls[9]
            alt_conf2 = item_ls[10]
            if (len(alt_conf2) > 1):
                alt_conf2 = "."
                item_ls.insert(10, alt_conf2)
            
            dev       = float(item_ls[12])
            sig       = float(item_ls[14])

            if (alt_conf1 == "."): alt_conf1 = ""
            if (alt_conf2 == "."): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["B-value",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        -99999., -99999., dev, sig, line]]])
               
            i += 1

    def get_ncs(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 5):
            # e.g. "Positional: A  12 LEU N   . deviation = 0.544 sigma= 0.050"
            # e.g. "B-value   : B  50 ASN CA  . deviation =20.000 sigma= 1.500"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[1]
            res_no1   = int(item_ls[2])
            res_name1 = item_ls[3]
            atom1     = item_ls[4]
            alt_conf1 = item_ls[5]
            if ('dev' in alt_conf1):
                alt_conf1 = "."
                item_ls.insert(5, alt_conf1)
                
            dev       = item_ls[7]
            sig       = item_ls[9]
                
            if (alt_conf1 == "."): alt_conf1 == ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["NCS",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        -999999, "", "", "",
                                                        -99999., -99999., dev, sig, line]]])       
               
            i += 1

    def get_sphericity(i):
        i += 3      # jump to interesting lines (only 3 here....)
        while (len(lines[i]) > 5):
            # e.g.          "A  26 CYS SG  B U-value= 0.2014 0.2329 0.2399 0.0179 0.0227-0.0064 Delta= 0.051 Sigma=  0.025"
            # new, e.g.
            # "A     19 GLU OE1 B U value= 0.6897 1.3874 0.6117-0.2763-0.0869 0.2513 Delta= 0.812 Sigma=  0.063"
            # "A     38 GLU OE1   U value= 1.4840 1.8504 0.7506-0.5906 0.7687-1.0238 Delta= 2.145 Sigma=  0.190"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[0]
            res_no1   = int(item_ls[1])
            res_name1 = item_ls[2]
            atom1     = item_ls[3]
            alt_conf1 = item_ls[4]
            if ('U' in alt_conf1):
                alt_conf1 = "."
                item_ls.insert(4, alt_conf1)

            u_mat_ls  = item_ls[6:11]
            dev       = item_ls[13]
            sig       = item_ls[15]

            if (alt_conf1 == "."): alt_conf1 = ""       # check, may be '.' if none

            if coot_utils.debug():
                print("BL DEBUG:: have item_ls", item_ls)
            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Sphericity",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        -999999, "", "", "",
                                                        -99999., -99999., dev, sig, line]]])       
                                                       
            i += 1

    def get_rigid(i):
        i += 3      # jump to interesting lines, again only 3
        while (len(lines[i]) > 5):
            # e.g.    "A  12 LEU N     - A  11 ASN C      Delta  =  4.625 Sigma=  2.000"
            # new
            # "A     19 GLU OE2 B - A     19 GLU CD  B  Delta  = 70.167 Sigma=  3.000 Dist=  1.270"
            # "A     15 ALA CB    - A     15 ALA CA     Delta  = 42.418 Sigma=  3.000 Dist=  1.542"
            line = lines[i]
            item_ls = split_clean(line)
            chain_id  = item_ls[0]
            res_no1   = int(item_ls[1])
            res_name1 = item_ls[2]
            atom1     = item_ls[3]
            alt_conf1 = item_ls[4]
            if (item_ls[4] == '-'):
                alt_conf1 = "."
                item_ls.insert(4, alt_conf1)

            chain_id2 = item_ls[6]
            res_no2   = int(item_ls[7])
            res_name2 = item_ls[8]
            atom2     = item_ls[9]
            alt_conf2 = item_ls[10]
            if (len(alt_conf2) > 1):
                # i.e. dont have alt conf but "Delta" or such instead
                alt_conf2 = "."
                item_ls.insert(10, alt_conf2)
            
            if coot_utils.debug():
                print("BL DEBUG:: have item list", item_ls)
                
            dev       = float(item_ls[12])
            sig       = float(item_ls[14])

            if (alt_conf1 == "."): alt_conf1 = ""
            if (alt_conf2 == "."): alt_conf2 = ""

            if coot_utils.debug():
                print("BL DEBUG:: have item list", item_ls)

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Rigid",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        -99999., -99999., dev, sig, line]]])
               
            i += 1

    # split a line in a list containing relevant information
    # need to additionally split things like e.g. merged numbers (e.g. 0.0227-0.0064 in U-value)
    # and chain id + res no (e.g. A1115 or mod.=-1.295)
    def split_clean(line):
        import re
        ret = []
        last_item = ""
        # first replace 'dodgy' characters by ' ' to be split (e.g. :, = and >)
        item_ls = line.replace(':',' ')
        item_ls = item_ls.replace('=', ' ')
        item_ls = item_ls.replace('>', ' ')
        # now split by white space
        item_ls = item_ls.split()
        for item in item_ls:
            start = item[0]
            end   = item[-1]

            # check for '-' in string
            if (item.find("-") > 0):
                # '-' in middle
                new_ls = item.split("-")
                ret.append(new_ls[0])
                for i in range(1, len(new_ls)):
                    ret.append("-" + new_ls[i])
            # check for merged string + no (and vice versa)
            elif (start.isalpha() and end.isdigit()
                  and not (last_item == 'at')
                  and not last_item.isalpha()):
                char = True
                i = 0
                string = start
                while char:
                    char = item[i].isalpha()
                    string += item[i]
                    i += 1
                ret.append(string)
                ret.append(item[i:-1])
            
            else:
                ret.append(item)
            last_item = item

        #print "BL DEBUG:: split clean returns", ret
        return ret

    # this sorts the initial list to look like the above mentioned extracted list
    # all information is conserved, i.e. 2 residues etc.!
    def sort_and_convert_list(ls):
        sorted_ls = ls
        sorted_ls.sort()
        merged_ls = []

        def find_dev_type(res, res_dev_type):
            dev_type = False
            position = False
            for type_ls in res[3:len(res)]:
                if (res_dev_type == type_ls[0]):
                    # same type
                    position = res.index(type_ls)
                    dev_type = res_dev_type
                    
            return dev_type, position
            
        
        for i in range(len(sorted_ls)):
            res = sorted_ls[i]
            merged_flag = False
            if ((i < len(sorted_ls) - 1)):
                next_res = sorted_ls[i + 1]
                if ((res[1] == next_res[1]) and (res[2] == next_res[2])):
                    # same chain and res_no
                    res_dev_type = next_res[3][0]
                    found_dev_type, pos = find_dev_type(res, res_dev_type)
                    if (found_dev_type):
                        # same type
                        res[pos].append(next_res[3][1])
                        sorted_ls[i + 1] = res
                    else:
                        # different/new type
                        res.append(next_res[3])
                        sorted_ls[i + 1] = res
                    merged_flag = True

                else:
                    #different chain and res_no
                    merged_ls.append(res)

            else:
                #last residue
                if not (merged_flag):
                    merged_ls.append(res)

        return merged_ls

    # here we actually convert the list with all information into something we can feed to
    # interesting_thing_gui, i.e.
    # FROM
    # [0, "A", 43,
    #   [bond_deviations, [atom, value...],[atom2, value2....]],
    #   [vdw_deviation, [atoms, value...]],
    #   [more deviations]],
    #   next res] 
    #
    # TO
    # e.g.  [["A 23: Chiral Deviation CA",0,"A",23,"","CA","A","fulltext_for_tooltip"],
    #        ["B 65: Bond Angle Deviation CA - CB",0,"B",65,"","CA","","fulltext_for_tooltip],
    #        ["Multiple deviations, x bond distances, y torsion angles",,0,"A",23,"","CA","A","fulltest_for_tooltip]]
    def make_interesting_things_list(list):

        ret_ls = []
            
        for res in list:
            imol     = res[0]
            chain_id = res[1]
            res_no   = res[2]
            ref_to_res = res_no

            dev_name = chain_id + " " + str(res_no) + ": "
            tooltip = ""

            no_diff_dev = len(res) - 3

            if ((no_diff_dev > 1) or (len(res[3]) > 2)):
                res_no2   = res[3][1][4]
                if (res_no2 != -999999 and res_no2 != res_no):
                    ref_to_res = res_no2
                # multiple deviations
                for dev_type in res[3:len(res)]:
                    # iterate the different deviation types
                    no_of_this_dev = len(dev_type) - 1
                    tmp_tip = ""
                    if (no_of_this_dev > 1):
                        plural_str = "s"
                        for item in dev_type[1:len(dev_type)]:
                            tmp_ls = item[12].split()
                            tmp_tip += " ".join(tmp_ls) + "\n"

                    else:
                        plural_str = ""
                        tmp_ls = dev_type[1][12].split()
                        tmp_tip = " ".join(tmp_ls) + "\n"                       

                    if (len(res) == 4):
                        #only one type of deviation
                        tooltip = tmp_tip
                        dev_name += str(no_of_this_dev) + " " + dev_type[0] + " deviation" + plural_str
                    else:
                        # multiple types
                        tooltip += str(no_of_this_dev) + " " + dev_type[0] + " deviation" + \
                                   plural_str + ":\n" + tmp_tip
                        dev_name = chain_id + " " + str(res_no) + ": Multiple deviations"
                        
                ins_code  = ""
                atom_name = res[3][1][2]         # use the first atom to centre (maybe should use CA or next atom
                alt_conf  = res[3][1][3]
                func = [refine_zone, imol, chain_id, res_no, ref_to_res, alt_conf]
                button_name = "Refine"

                tooltip = tooltip.rstrip("\n")
                ret_ls.append([dev_name, imol, chain_id, res_no, ins_code, atom_name, alt_conf,
                               func, button_name, tooltip])  # with a fix button and tooltip
                
            else:
                # single deviation
                ins_code   = ""
                atom_name  = res[3][1][2]
                alt_conf   = res[3][1][3]
                dev_name   = dev_name + res[3][0] + " deviation: " + atom_name
                atom_name2 = res[3][1][6]
                if (atom_name2):
                    res_no2 = res[3][1][4]
                    dev_name += "-"
                    if (res_no2 == res_no):
                        # difference in same residue
                        dev_name += atom_name2
                        func = [refine_zone, imol, chain_id, res_no, res_no, alt_conf]
                    else:
                        # difference to other residue
                        dev_name = dev_name + str(res_no2) + " " + atom_name2
                        func = [refine_zone, imol, chain_id, res_no, res_no2, alt_conf]
                else:
                    func = [refine_zone, imol, chain_id, res_no, res_no, alt_conf]

                button_name = "Refine"
                tooltip = res[3][1][12]
                tmp_ls = tooltip.split()
                tooltip = " ".join(tmp_ls)
                ret_ls.append([dev_name, imol, chain_id, res_no, ins_code, atom_name, alt_conf,
                               func, button_name, tooltip])  # with a fix button and tooltip

        return ret_ls


    # main
    interesting_list = []
    refmac_all_dev_list = []
    warning_info_list = []
    lines = False
    cgmat_cycle = "CGMAT cycle number ="
    found_last = False
    last_cycle_finished = False
    ncyc = 5
    reg_ncyc      = re.compile("data line.*ncyc", re.IGNORECASE)
    reg_ncyc_only = re.compile("ncyc", re.IGNORECASE)

    # read the log file
    filename = os.path.normpath(os.path.abspath(refmac_log_file))
    if (os.path.isfile(filename)):
        try:
            fin = open(filename, 'r')
            lines = fin.readlines()
            fin.close()
        except:
            print("ERROR opening the filename ", filename)

        if (lines):
            i = 0
            no_lines = len(lines)
            while (i < no_lines):
                line = lines[i]
                
                ncycle_res = reg_ncyc.search(line)
                if (ncycle_res):     
                    # get no of cycles
                    line_ls = line.split()
                    for item in line_ls:
                        ncyc_pos = reg_ncyc_only.search(item)
                        if ('NCYC' in item.upper()):
                            ncyc = int(line_ls[line_ls.index(item) + 1])
                            break
                    #pos = ncycle_res.end()
                    #ncyc = int(line[pos: pos + 4])

                # find the Warning and Info lines at the beginning
                if   ("WARNING :" in line):
                    get_warning(i)
                elif ("INFO:" in line):
                    get_info(i)

                if (cgmat_cycle in line):
                    # we ignore all outliers before the end
                    cycle_number = int(line[(line.index(cgmat_cycle) + 20):(len(line) - 1)])
                    if (cycle_number == ncyc):
                        found_last = True

                if (found_last):
                    if ("Overall figure of merit" in line):
                        # last cycle is finished
                        last_cycle_finished = True

                if (found_last and last_cycle_finished):
                    # now extract info
                    if   ("Bond distance outliers" in line):
                        get_bond_dist(i)
                    elif ("Bond angle outliers" in line):
                        get_bond_angle(i)
                    elif ("Torsion angle outliers" in line):
                        get_torsion(i)
                    elif ("Chiral volume outliers" in line):
                        get_chiral(i)
                    elif ("Large deviation of atoms from planarity" in line):
                        get_planarity(i)
                    elif ("VDW outliers" in line):
                        get_vdw(i)
                    elif ("B-value outliers" in line):
                        get_b_value(i)
                    elif ("NCS restraint outliers" in line):
                        get_ncs(i)
                    elif ("Sphericity outliers" in line):
                        get_sphericity(i)
                    elif ("Rigid bond outliers" in line):
                        get_rigid(i)

                i += 1       # next line



        converted_refmac_all_dev_list = sort_and_convert_list(refmac_all_dev_list)

        deviation_list = make_interesting_things_list(converted_refmac_all_dev_list)

        interesting_list += warning_info_list

        interesting_list += deviation_list

        if interesting_list:
            try:
                # no threads.
                # coot_gui.run_with_gtk_threading(interesting_things_with_fix_maybe, "Refmac outliers etc.", interesting_list)
                print("No threads for interesting_things_with_fix_maybe") # fixme
                pass
            except:
                print("BL INFO:: could not show the interesting Refmac things! Probably no PyGTK!")
        else:
            print("BL INFO:: no deviations and nothing otherwise interesting found in ", filename)
        
    else:
        print("WARNING your file with name %s does not seem to exist" %filename)
        
#read_refmac_log(0, "25_refmac5.log")
#read_refmac_log(0, "refmac-from-coot-0.log")

global refmac_version_cached
refmac_version_cached = False
# returns the major refmac version, i.e. [5, 3] (for 5.3) or False if no refmac found
#
def get_refmac_version():

    global refmac_version_cached
    if refmac_version_cached:
        return refmac_version_cached
    else:

        # maybe want to cache these too!?
        refmac_execfile = coot_utils.find_exe("refmac5", "CBIN", "CCP4_BIN", "PATH")

        if (refmac_execfile):
            log_file = "refmac_version_tmp.log"
            status = coot_utils.popen_command(refmac_execfile, ["-i"], [], log_file)
            if (status == 0):
                fin = open(log_file, 'r')
                lines = fin.readlines()
                fin.close()
                os.remove(log_file)
                for line in lines:
                    if ("Program" in line):
                        version = line.split()[-1]
                        ret = list(map(int, version.split(".")))
                        refmac_version_cached = ret
                        return ret
            else:
                os.remove(log_file)
                print("INFO:: problem to get refmac version")
        else:
            return False


# Writes an input file for refmac to refine occupancies. Takes all
# residues with alt confs (each will be one group). Complete groups
# if sum occ ~1 otherwise incomplete. Takes solvent too.
#
def restraints_for_occupancy_refinement(imol, file_name="refmac_extra_params.txt",
                                        n_cycle=1):

    # simple version for now.
    # take each residue and make all alt confs to be one
    # ignore close contact ones
    res_list = coot_utils.residues_with_alt_confs(imol)
    group_id = 1
    if (not res_list):
        print("BL INFO:: no alt confs")
        return
    else:
        f = open(file_name, 'w')
        for (chain_id, res_no, ins_code) in res_list:
            atom_name = 'dummy'
            if coot_utils.is_protein_chain_qm(imol, chain_id):
                atom_name = ' CA '
            if coot_utils.is_solvent_chain_qm(imol, chain_id):
                atom_name = ' O  '
            alt_confs = coot_utils.residue_alt_confs(imol, chain_id, res_no, ins_code)
            try:
                # remove no alt confs, assums always CA is alt
                alt_confs.remove('') 
            except:
                pass # no atoms without alt confs
            group_list = []
            total_occ = 0.
            # only if we know the atom type for now
            if (atom_name != "dummy"):
                write_it = 0
                for alt_conf in alt_confs:
                    try:
                        # add CA occupancies
                        atom_spec = coot_utils.atom_specs(imol, chain_id, res_no, ins_code,
                                               atom_name, alt_conf)
                        total_occ += atom_spec[0]
                        line = "occupancy group id " + str(group_id) + \
                               " chain " + chain_id + \
                               " residue " + str(res_no) + \
                               " alt " + alt_conf + "\n"
                        f.write(line)                
                        group_list.append(group_id)
                        group_id += 1
                        write_it = 1
                    except:
                        # some problem getting atoms.
                        pass

                complete = ""
                if (abs(total_occ - 1.) > 0.1):
                    complete = "in"
                line = "occupancy group alts " + complete + "complete " + \
                       " ".join(map(str, group_list)) + "\n"
                if write_it:
                    f.write(line)
        # we could have cycle etc too
        f.write("occupancy refine ncycle "+ str(n_cycle) + "\n") 
        f.write("occupancy refine\n")
        f.close()
        
        
