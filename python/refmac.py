# refmac.py 
#
# Copyright 2005, 2006 by Bernhard Lohkamp
# Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
#

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
# 1.) .coot.py (needs to include: global refmac_extra_params)
#
# 2.) whatever is written here
#
# 3.) refmac-extra-params.txt
#
# deftexi refmac_extra_params

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


# BL says: for loggraph we use a function
def run_loggraph(logfile):
    import os
    # this is the easy code, with no exception handling, so commented out (next 3 lines):
    # loggraph_exe = os.path.join(os.environ['CCP4I_TOP'],'bin/loggraph.tcl')
    # bltwish_exe = os.path.join(os.environ['CCP4I_TCLTK'],'bltwish')
    # os.spawnl(os.P_NOWAIT, bltwish_exe , bltwish_exe , loggraph_exe , logfile)

    # now lets get some exception handling in:
    try:
        # 1st check TCL
        os.environ['CCP4I_TCLTK']
        if (os.name == 'nt'):
            # Windows
            bltwish_exe = os.path.join(os.environ['CCP4I_TCLTK'],'bltwish') + '.exe'
            #some windows machines have bltwish not in CCP4I_TCLTK but $/bin
            bltwish_exe_win_alt = os.path.join(os.environ['CCP4I_TCLTK'],'bin','bltwish') + '.exe'
            if (os.path.isfile(bltwish_exe_win_alt)):
                bltwish_exe = bltwish_exe_win_alt
        else:
            bltwish_exe = os.path.join(os.environ['CCP4I_TCLTK'],'bltwish')

        if (os.path.isfile(bltwish_exe)):
            
           # now that we have tcl + btlwish we check for loggraph.tcl
           try:
               os.environ['CCP4I_TOP']
               loggraph_exe = os.path.join(os.environ['CCP4I_TOP'],'bin/loggraph.tcl')
               if os.path.isfile(loggraph_exe):
                   print 'BL DEBUG:: We have allocated everything to run loggraph and shall do that now'
                   os.spawnl(os.P_NOWAIT, bltwish_exe , bltwish_exe , loggraph_exe , logfile)
               else:
                   print 'We have $CCP4I_TOP but we cannot find bin/loggraph.tcl'
           except KeyError:
               print 'We cannot allocate $CCP4I_TOP so no loggraph available'
        else:
           print 'We have $CCP4I_TCLTK but we cannot find bltwish.exe'
    except KeyError:
        print 'We cannot allocate $CCP4I_TCLTK so no bltwish available'


def run_refmac_by_filename(pdb_in_filename, pdb_out_filename, mtz_in_filename, mtz_out_filename, extra_cif_lib_filename, imol_refmac_count, swap_map_colours_post_refmac_p, imol_mtz_molecule, show_diff_map_flag, phase_combine_flag, phib_fom_pair, force_n_cycles, ccp4i_project_dir, f_col, sig_f_col, r_free_col):

    global refmac_count
    global refmac_extra_params
    import os, stat, operator

    refmac_execfile = find_exe("refmac5","CCP4_BIN","PATH")

    labin_string = "LABIN FP=" + f_col + " SIGFP=" + sig_f_col
    if (r_free_col != "") :
        labin_string += " FREE=" + r_free_col
# BL says: command line args have to be string not list here
    command_line_args = " XYZIN " + pdb_in_filename \
                        + " XYZOUT " + pdb_out_filename + " HKLIN " + mtz_in_filename \
                        + " HKLOUT " + mtz_out_filename

    if (extra_cif_lib_filename != "") :
        command_line_args += " LIBIN " + extra_cif_lib_filename

    std_lines = ["MAKE HYDROGENS NO"] # Garib's suggestion 8 Sept 2003

    if operator.isNumberType(force_n_cycles):
       if force_n_cycles >=0:
          std_lines.append("NCYCLES " + str(force_n_cycles))
    if (refmac_extra_params):
        extra_params = refmac_extra_params
    else:
        extra_params = get_refmac_extra_params()

    if (extra_params == None): extra_params = [""]
    data_lines = std_lines + extra_params

    if (not extra_params_include_weight_p(extra_params)) :
        data_lines.append("WEIGHT AUTO")

    data_lines.append(labin_string)

    print "refmac extra params: ", extra_params

    refmac_log_file_name = ""
    if (len(ccp4i_project_dir) > 0) :
        refmac_log_file_name += ccp4i_project_dir
    refmac_log_file_name += "refmac-from-coot-"
    refmac_log_file_name += str(refmac_count)
    refmac_log_file_name += ".log"

    refmac_count += 1

    print "INFO: running refmac with these command line args: ", command_line_args
    print "INFO: running refmac with these data lines: ", data_lines
    try:
        print "environment variable:  SYMOP: ", os.environ['SYMINFO']
    except:
        print " not set !"
    try:
        print "environment variable: ATOMSF: ", os.environ['ATOMSF'] 
    except:
        print " not set !"
    try:
        print "environment variable:  CLIBD: ", os.environ['CLIBD'] 
    except:
        print " not set !"
    try:
        print "environment variable:   CLIB: ", os.environ['CLIB'] 
    except:
        print " not set !"

    # BL says: 
    # we have to write a file with datalines in it, which we read in
    # couldnt find away around, yet
    refmac_data_file = "coot-refmac-data-lines.txt"
    input = file (refmac_data_file,'w')
    for data in data_lines:
        input.write(data + '\n')
    input.close()

    status = os.popen(refmac_execfile + command_line_args + ' < ' + refmac_data_file + ' > ' + refmac_log_file_name, 'r') # to screen too.
    refmac_finished = status.close()

    if (refmac_finished != None) : # refmac ran fail...
        print "Refmac Failed."

    else : # refmac ran OK.

        # BL says: 
        # popen will always write a log file, so we check the size, is 0 if failed
        refmac_log_size = os.stat(refmac_log_file_name)[stat.ST_SIZE]

        if (refmac_log_size > 0) :
            run_loggraph(refmac_log_file_name)

        # deal with R-free column
        if (r_free_col == ""): r_free_bit = ["",0]
        else: r_free_bit = [r_free_col,1]

        args = [mtz_out_filename,"FWT","PHWT","",0,0,1,f_col,sig_f_col] + r_free_bit

        recentre_status = recentre_on_read_pdb()
        novalue = set_recentre_on_read_pdb(0)
        print "DEBUG:: recentre status: ", recentre_status
        imol = handle_read_draw_molecule(pdb_out_filename)

        if (recentre_status) :
            set_recentre_on_read_pdb(1)
        set_refmac_counter(imol, imol_refmac_count + 1)

        new_map_id = make_and_draw_map_with_refmac_params(*args) 

	# set the new map as refinement map
	set_imol_refinement_map(new_map_id)
	# now run the refmac_log_reader to get interesting things
	read_refmac_log(imol, refmac_log_file_name)

        if (swap_map_colours_post_refmac_p == 1) :
            swap_map_colours(imol_mtz_molecule, new_map_id)

        # difference map (flag was set):
        if (show_diff_map_flag == 1) :
           args = [mtz_out_filename,"DELFWT","PHDELWT","",0,1,1,f_col,sig_f_col] + r_free_bit

           make_and_draw_map_with_refmac_params(*args)

        
            
# Return True if the list of strings @var{params_list} contains a
# string beginning with "WEIGHT".  If not return False
#
def extra_params_include_weight_p(params_list):

   have_weight = params_list.count('WEIG')
   if (have_weight == 1):
     return True
   elif (have_weight == 0):
     return False
   else:
     print 'BL WARNING:: This shouldn\'t happen, we have more than one weight defined. God knows what refmac will do now...'
     return False


# If refmac-extra-params is defined (as a list of strings), then
# return that, else read the file "refmac-extra-params.txt"
#
# Return a list a list of strings.
#
def get_refmac_extra_params():
  global refmac_extra_params

  extras_file_name = "refmac-extra-params.txt"
  try:
     f = open(extras_file_name,'r')
  except IOError:
     print 'BL WARNING:: we dont have refmac-extra-params file'
  else:
     refmac_extra_params = f.readlines()
     print 'BL INFO:: we have refmac-extra-params file and read lines'
     f.close()

#  print "BL INFO:: we added ", refmac_extra_params
     extra_params = refmac_extra_params
     return extra_params

#
def run_refmac_for_phases(imol, mtz_file_name, f_col, sig_f_col):

	import os

	if os.path.isfile(mtz_file_name):
		if valid_model_molecule_qm(imol):
			dir_state = make_directory_maybe("coot-refmac")
			if not dir_state == 0:
				print "Failed to make coot-refmac directory\n"
			else:
				stub = "coot-refmac/refmac-for-phases"
				pdb_in = stub + ".pdb"
				pdb_out = stub + "-tmp.pdb"
				mtz_out = stub + ".mtz"
				cif_lib_filename = ""
				write_pdb_file(imol, pdb_in)
				run_refmac_by_filename(pdb_in, pdb_out,
					mtz_file_name, mtz_out, 
					cif_lib_filename, 0, 0, -1,
					1, 0, [], 0, "",
					f_col, sig_f_col, "")

		else:
			print "BL WARNING:: no valid model molecule!"
	else:
		print "BL WARNING:: mtzfile %s not found" %mtz_file_name


def refmac_for_phases_and_make_map(mtz_file_name, f_col, sig_f_col):

	import os

	if os.path.isfile(mtz_file_name):
		molecule_chooser_gui("  Choose a molecule from which to calculate Structure factors:  ",
			lambda imol: run_refmac_for_phases(imol, mtz_file_name, f_col, sig_f_col))	


# This will read a refmac log file and extract information in
# interesting_things_gui format, i.e. a list with either coorinates or residues
# to be flagged as interesting
#
# the list of extracted information may look like:
# [[0, "A", 43,
#   [bond_deviations, [atom, value],[atom2, value2]],
#   [vdw_deviation, [atoms, value]],
#   [more deviations]],
#   next res] 
#
# e.g. interesting_things_gui(
#       "Bad things by Analysis X",
#       [["Bad Chiral",0,"A",23,"","CA","A"],
#        ["Bad Density Fit",0,"B",65,"","CA",""],
#        ["Interesting blob",45.6,46.7,87.5]]
#
def read_refmac_log(imol, refmac_log_file):
    
    import os
    import re
    
    def get_warning(i):
        this_line = lines[i]
        next_line = lines[i+1]
        if   ("CIS" in this_line):
            # CIS peptide found
            # we ignore multiple conformations for now...
            chain     = next_line[15:17]       # this is e.g. "AA"
            chain_id  = chain[1]
            res_no1   = int(next_line[22:26])
            res_name1 = next_line[28:31]
            res_no2   = int(next_line[40:44])
            res_name2 = next_line[46:49]
            angle_str = this_line[49:56]

            info_text = "pre-WARNING: CIS bond: " + chain_id + " " + str(res_no1) + " - " \
                        + str(res_no2) + angle_str

            warning_info_list.append([info_text, imol, chain_id, res_no1, "", " CA ", ""])
            
        elif ("gap" in this_line):
            # gap-link
            chain     = this_line[30:32]       # this is e.g. "AA"
            chain_id  = chain[1]
            res_no1   = int(this_line[37:41])
            res_name1 = next_line[43:46]
            res_no2   = int(next_line[54:58])
            res_name2 = next_line[60:63]
            extra_info = ""
            if ("-link will" in next_line):
                # we have extra info
                extra_info = "'gap'-link created"

            info_text = "pre-WARNING: connection gap: " + chain_id + " " \
                        + str(res_no1) + " - " + str(res_no2) + extra_info

            warning_info_list.append([info_text, imol, chain_id, res_no1, "", " CA ", ""])

        elif ("big distance" in this_line):
            # large distance
            chain     = next_line[15:17]       # this is e.g. "AA"
            chain_id  = chain[1]
            res_no1   = int(next_line[24:28])
            res_name1 = next_line[30:33]
            res_no2   = int(next_line[41:45])
            res_name2 = next_line[47:50]
            dist_str  = this_line[52:59]

            info_text = "pre-WARNING: large distance: " + chain_id + " " \
                        + str(res_no1) + " - " + str(res_no2) + dist_str

            warning_info_list.append([info_text, imol, chain_id, res_no1, "", " CA ", ""])
            
        # elif(): not sure what other warnings there could be...

    def get_info(i):
        this_line = lines[i]
        next_line = lines[i+1]
        if   ("link is found" in this_line):
            # link found but not used
            # we ignore multiple conformations for now...
            chain1     = next_line[15:17]       # this is e.g. "AA"
            chain_id1  = chain1[1]              # we just take the second char
            res_no1    = int(next_line[22:26])
            res_name1  = next_line[28:31]
            atom_name1 = next_line[40:43]
            alt_conf1  = next_line[44]          # guess this is alt_conf ("." here)
            chain2     = next_line[50:52]
            chain_id2  = chain2[1]
            res_no2    = int(next_line[57:61])
            res_name2  = next_line[63:66]
            atom_name2 = next_line[75:78]
            alt_conf2  = next_line[79]

            if (alt_conf1 == "."):
                alt_conf1 = ""

            info_text = "pre-INFO: not used link: " \
                        + chain_id1 + " " + str(res_no1) + " " + atom_name1 + " - " \
                        + chain_id2 + " " + str(res_no2) + " " + atom_name2

            warning_info_list.append([info_text, imol, chain_id1, res_no1, "", atom_name1, alt_conf1])
            
        # elif ("" in this_line): not sure what other INFO there may be

    def get_bond_dist(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 2):
            # e.g. "A  15 ARG C   . - A  15 ARG O   . mod.= 1.295 id.= 1.231 dev= -0.064 sig.= 0.020"
            line = lines[i]
            chain_id  = line[0]
            res_no1   = int(line[1:5])
            res_no2   = int(line[19:23])
            res_name1 = line[6:9]
            res_name2 = line[24:27]
            atom1     = line[9:13]
            atom2     = line[27:31]
            alt_conf1 = line[14]
            alt_conf2 = line[32]
            
            mod       = float(line[39:45])
            ideal     = float(line[50:56])
            dev       = float(line[61:68])
            sig       = float(line[74:80])

            if (alt_conf1 == "."): alt_conf1 = ""
            if (alt_conf2 == "."): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Bond distance",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        mod, ideal, dev, sig]]])
               
            i += 1

    def get_bond_angle(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 2):
            # e.g. "A  50 GLU O   B - A  51 LEU N     mod.= 100.81 id.= 123.00 dev= 22.193 sig.=  1.600"
            line = lines[i]
            chain_id  = line[0]
            res_no1   = int(line[1:5])
            res_no2   = int(line[19:23])
            res_name1 = line[6:9]
            res_name2 = line[24:27]
            atom1     = line[9:13]
            atom2     = line[27:31]
            alt_conf1 = line[14]
            alt_conf2 = line[32]
            
            mod       = float(line[39:46])
            ideal     = float(line[51:58])
            dev       = float(line[63:70])
            sig       = float(line[77:83])

            if (alt_conf1 == " "): alt_conf1 = ""
            if (alt_conf2 == " "): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Bond angle",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        mod, ideal, dev, sig]]])
               
            i += 1

    def get_torsion(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 2):
            # e.g. "A  11 ASN CA    - A  12 LEU CA    mod.=-167.65 id.= 180.00 per.= 1 dev=-12.352 sig.=  3.000"
            line = lines[i]
            chain_id  = line[0]
            res_no1   = int(line[1:5])
            res_no2   = int(line[19:23])
            res_name1 = line[6:9]
            res_name2 = line[24:27]
            atom1     = line[9:13]
            atom2     = line[27:31]
            alt_conf1 = line[14]            # guess
            alt_conf2 = line[32]            # guess
            
            mod       = float(line[39:46])
            ideal     = float(line[51:58])
            period    = int(line[64:66])    # ignored for now
            dev       = float(line[71:78])
            sig       = float(line[85:92])

            if (alt_conf1 == " "): alt_conf1 = ""
            if (alt_conf2 == " "): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Torsion angle",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        mod, ideal, dev, sig]]])
               
            i += 1

    def get_chiral(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 2):
            # e.g. "A  51 LEU CG    mod.=   2.79 id.=  -2.59 dev= -5.375 sig.=  0.200"
            line = lines[i]
            chain_id  = line[0]
            res_no1   = int(line[1:5])
            res_name1 = line[6:9]
            atom1     = line[9:13]
            alt_conf1 = line[14]
            
            mod       = float(line[21:28])
            ideal     = float(line[33:40])
            dev       = float(line[45:52])
            sig       = float(line[58:65])

            if (alt_conf1 == " "): alt_conf1 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Chiral",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        -999999, "", "", "",
                                                        mod, ideal, dev, sig]]])
               
            i += 1

    def get_planarity(i):
        i += 3      # jump to interesting lines (only 3 here....)
        while (len(lines[i]) > 2):
            # e.g. "Atom: A  59 ASP C   B deviation=   0.31 sigma.=   0.02"
            line = lines[i]
            chain_id  = line[6]
            res_no1   = int(line[7:11])
            res_name1 = line[12:15]
            atom1     = line[15:19]
            alt_conf1 = line[20]            # guess
            dev       = float(line[32:39])
            sig       = float(line[47:54])

            if (alt_conf1 == " "): alt_conf1 = ""       # guess

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Planarity",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        -999999, "", "", "",
                                                        -99999., -99999., dev, sig]]])       
               
            i += 1

    def get_vdw(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 2):
            # e.g. "A  26 CYS SG  A - A  75 ILE CD1 . mod.= 2.812 id.= 3.820 dev= -1.008 sig.= 0.300"
            line = lines[i]
            chain_id  = line[0]
            res_no1   = int(line[1:5])
            res_no2   = int(line[19:23])
            res_name1 = line[6:9]
            res_name2 = line[24:27]
            atom1     = line[9:13]
            atom2     = line[27:31]
            alt_conf1 = line[14]
            alt_conf2 = line[32]
            
            mod       = float(line[39:45])
            ideal     = float(line[50:56])
            dev       = float(line[61:68])
            sig       = float(line[74:79])  # there seem to be a difference between Garib's doc and refmac v. 5.2.0019
                                            # 5.2.0019 requires 74:79, Garib's docu suggests 74:80

            if (alt_conf1 == "."): alt_conf1 = ""
            if (alt_conf2 == "."): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["VDW",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        mod, ideal, dev, sig]]])
               
            i += 1

    def get_b_value(i):
        i += 3      # jump to interesting lines
        while (len(lines[i]) > 2):
            # e.g. "B   5 PHE N     - B   4 GLN C        ABS(DELTA)= 15.990   Sigma=  1.500"
            line = lines[i]
            chain_id  = line[0]
            res_no1   = int(line[1:5])
            res_no2   = int(line[19:23])
            res_name1 = line[6:9]
            res_name2 = line[24:27]
            atom1     = line[9:13]
            atom2     = line[27:31]
            alt_conf1 = line[14]
            alt_conf2 = line[32]
            
            dev       = float(line[48:55])
            sig       = float(line[64:71])

            if (alt_conf1 == " "): alt_conf1 = ""
            if (alt_conf2 == " "): alt_conf2 = ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["B-value",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        -99999., -99999., dev, sig]]])
               
            i += 1

    def get_ncs(i):
        i += 4      # jump to interesting lines
        while (len(lines[i]) > 2):
            # e.g. "Positional: A  12 LEU N   . deviation = 0.544 sigma= 0.050"
            # e.g. "B-value   : B  50 ASN CA  . deviation =20.000 sigma= 1.500"
            line = lines[i]
            ncs_type  = line[0:10]
            chain_id  = line[12]
            res_no1   = int(line[13:17])
            res_name1 = line[17:22]
            atom1     = line[22:24]
            alt_conf1 = line[26]
            
            dev       = float(line[39:45])
            sig       = float(line[52:58])

            if (alt_conf1 == "."): alt_conf1 == ""

            refmac_all_dev_list.append([imol, chain_id, res_no1, [["NCS", ncs_type],
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        -999999, "", "", "",
                                                        -99999., -99999., dev, sig]]])       
               
            i += 1

    def get_sphericity(i):
        i += 3      # jump to interesting lines (only 3 here....)
        while (len(lines[i]) > 2):
            # e.g. "A  26 CYS SG  B U-value= 0.2014 0.2329 0.2399 0.0179 0.0227-0.0064 Delta= 0.051 Sigma=  0.025"
            line = lines[i]
            chain_id  = line[0]
            res_no1   = int(line[1:5])
            res_name1 = line[6:9]
            atom1     = line[10:12]
            alt_conf1 = line[14]
            u_mat_str = line[24:66]                     # not used currently
            dev       = float(line[73:79])
            sig       = float(line[86:93])

            if (alt_conf1 == " "): alt_conf1 = ""       # check, may be '.' if none

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Sphericity",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        -999999, "", "", "",
                                                        -99999., -99999., dev, sig]]])       
                                                       
            i += 1

    def get_rigid(i):
        i += 3      # jump to interesting lines, again only 3
        while (len(lines[i]) > 2):
            # e.g. "A  12 LEU N     - A  11 ASN C      Delta  =  4.625 Sigma=  2.000"
            line = lines[i]
            chain_id  = line[0]
            res_no1   = int(line[1:5])
            res_no2   = int(line[19:23])
            res_name1 = line[6:9]
            res_name2 = line[24:27]
            atom1     = line[9:14]
            atom2     = line[27:32]
            alt_conf1 = line[14]          #check
            alt_conf2 = line[32]          #check
            dev       = float(line[43:50])
            sig       = float(line[57:63])

            refmac_all_dev_list.append([imol, chain_id, res_no1, ["Rigid",
                                                       [res_no1, res_name1, atom1, alt_conf1,
                                                        res_no2, res_name2, atom2, alt_conf2,
                                                        -99999., -99999., dev, sig]]])
               
            i += 1

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
            if ((i < len(sorted_ls) - 1)):
                merged_flag = False
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
    # e.g.  [["A 23: Chiral Deviation CA",0,"A",23,"","CA","A"],
    #        ["B 65: Bond Angle Deviation CA - CB",0,"B",65,"","CA",""],
    #        ["Multiple deviations, x bond distances, y torsion angles",,0,"A",23,"","CA","A"]]
    def make_interesting_things_list(list):

        ret_ls = []
            
        for res in list:
            imol     = res[0]
            chain_id = res[1]
            res_no   = res[2]
            ref_to_res = res_no

            dev_name = chain_id + " " + str(res_no) + ": "

            no_diff_dev = len(res) - 3

            if ((no_diff_dev > 1) or (len(res[3]) > 2)):
                res_no2   = res[3][1][4]
                if (res_no2 != -999999 and res_no2 != res_no):
                    ref_to_res = res_no2
                # multiple deviations
                for dev_type in res[3:len(res)]:
                    # iterate the different deviation types
                    no_of_this_dev = len(dev_type) - 1
                    if (no_of_this_dev > 1):
                        plural_str = "s"
                    else:
                        plural_str = ""
                    if (res.index(dev_type) < (len(res) - 1)):
                        end_str = ", "
                    else:
                        end_str = ""
                    dev_name += str(no_of_this_dev) + " " + dev_type[0] + " deviation" + plural_str + end_str
                    
                ins_code  = ""
                atom_name = res[3][1][2]         # use the first atom to centre (maybe should use CA or next atom
                alt_conf  = res[3][1][3]
                func = [refine_zone, imol, chain_id, res_no, ref_to_res, alt_conf]
                button_name = "Refine"
                ret_ls.append([dev_name, imol, chain_id, res_no, ins_code, atom_name, alt_conf, func, button_name])  # maybe a fix button
                
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
                        dev_name = dev_name + " " + str(res_no2) + atom_name2
                        func = [refine_zone, imol, chain_id, res_no, res_no2, alt_conf]

                button_name = "Refine"
                ret_ls.append([dev_name, imol, chain_id, res_no, ins_code, atom_name, alt_conf, func, button_name])  # maybe a fix button

        return ret_ls


    # main
    interesting_list = []
    refmac_all_dev_list = []
    warning_info_list = []
    lines = False
    cgmat_cycle = "CGMAT cycle number ="
    found_last = False
    last_cycle_finished = False
    ncyc = False
    reg_ncyc= re.compile("data line.*ncyc", re.IGNORECASE)

    # read the log file
    filename = os.path.normpath(os.path.abspath(refmac_log_file))
    if (os.path.isfile(filename)):
        try:
            fin = open(filename, 'r')
            lines = fin.readlines()
            fin.close()
        except:
            print "ERROR opening the filename ", filename

        if (lines):
            i = 0
            no_lines = len(lines)
            while (i < no_lines):
                line = lines[i]
                
                ncycle_res = reg_ncyc.search(line)
                if (ncycle_res):     
                    # get no of cycles
                    pos = ncycle_res.end()
                    ncyc = int(line[pos: pos + 4])

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
                    if ("ML based su of thermal parameters" in line):
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

        if (interesting_list):
            interesting_things_with_fix_maybe("Refmac outliers etc.", interesting_list)
        else:
            print "BL INFO:: no deviations and nothing otherwise interesting found in ", filename
        
    else:
        print "WARNING your file with name %s does not seem to exist" %filename
        
#read_refmac_log(0, "25_refmac5.log")
#read_refmac_log(0, "refmac-from-coot-0.log")
