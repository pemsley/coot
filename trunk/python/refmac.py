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
    command_line_args = " XYZIN " + pdb_in_filename
    + " XYZOUT " + pdb_out_filename + " HKLIN " + mtz_in_filename
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
    try: print "environment variable:  SYMOP: ", os.environ['SYMINFO']
    except: print " not set !"
    try: print "environment variable: ATOMSF: ", os.environ['ATOMSF'] 
    except: print " not set !"
    try: print "environment variable:  CLIBD: ", os.environ['CLIBD'] 
    except: print " not set !"
    try: print "environment variable:   CLIB: ", os.environ['CLIB'] 
    except: print " not set !"

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


