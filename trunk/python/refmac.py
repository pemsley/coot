# refmac wrapper for coot
#
# Copyright 2005 by Bernhard Lohkamp
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
# in contrast to the scheme refmac the "refmac-extra-params" file 
# overwrites the parameters given here
# also params given in .coot.py are considered very first, so check is:
# 1.) refmac-extra-params.txt
# 2.) .coot.py (needs to include: global refmac_extra_params)
# 3.) whatever is written here

#refmac_extra_params has to be a list
# this is an example, feel free to change it (uncomment next 2 lines)!!
#global refmac_extra_params
#refmac_extra_params = ['NCYC 0', 'WEIGHT 0.3']


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
  bltwish_exe = os.path.join(os.environ['CCP4I_TCLTK'],'bltwish')
  bltwish_exe_win = os.path.join(os.environ['CCP4I_TCLTK'],'bltwish') + '.exe'
#some windows machines have bltwish not in CCP4I_TCLTK but $/bin
  bltwish_exe_win_alt = os.path.join(os.environ['CCP4I_TCLTK'],'bin','bltwish') + '.exe'
  if (os.path.isfile(bltwish_exe) or os.path.isfile(bltwish_exe_win) or os.path.isfile(bltwish_exe_win_alt)):
    if (os.path.isfile(bltwish_exe_win_alt)): bltwish_exe = bltwish_exe_win_alt
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



def run_refmac_by_filename(pdb_in_filename,pdb_out_filename,mtz_in_filename,mtz_out_filename,extra_cif_lib_filename,imol_refmac_count,swap_map_colours_post_refmac_p,imol_mtz_molecule,show_diff_map_flag,phase_combine_flag,phib_fom_pair,ccp4i_project_dir,f_col,sig_f_col,r_free_col) :

    global refmac_count
    global refmac_extra_params
    import os
    import stat
# BL says: we assume the refmac5 to be used is in PATH
    refmac_execfile = "refmac5"

    print "got args: ",pdb_in_filename,pdb_out_filename,mtz_in_filename,mtz_out_filename,extra_cif_lib_filename,imol_refmac_count,swap_map_colours_post_refmac_p,imol_mtz_molecule,show_diff_map_flag,phase_combine_flag,phib_fom_pair,ccp4i_project_dir,f_col,sig_f_col,r_free_col


#     local_r_free_col = (r_free_col == "") "" : r_free_col
    labin_string = "LABIN FP=" + f_col + " SIGFP=" + sig_f_col
    if (r_free_col != "") :
        labin_string += " FREE=" + r_free_col
# BL says: command line args have to be string not list here
    command_line_args = " XYZIN " + pdb_in_filename + " XYZOUT " + pdb_out_filename + " HKLIN " + mtz_in_filename + " HKLOUT " + mtz_out_filename
    if (extra_cif_lib_filename != "") :
        command_line_args += " LIBIN " + extra_cif_lib_filename
#        command_line_args.append(extra_cif_lib_filename)
    std_lines = ["MAKE HYDROGENS NO"]
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
    refmac_log_file_name += "refmac-from-coot"
    refmac_log_file_name += str(refmac_count)
    refmac_log_file_name += ".log"

    refmac_count += 1

    print "INFO: running refmac with these command line args: ", command_line_args
    print "INFO: running refmac with these data lines: ", data_lines

#    refmac_all_arguments = command_line_args
#    refmac_all_arguments +=data_lines
#    refmac_all_arguments.append(refmac_log_file_name)
#    print "BL INFO: all_arguments are :", refmac_all_arguments

# BL says: we have to write a file with datalines in it, couldnt find away around, yet
    refmac_data_file = "coot-refmac-data-lines.txt"
    input = file (refmac_data_file,'w')
    for data in data_lines:
        input.write(data + '\n')
    input.close()

    status = os.popen(refmac_execfile + command_line_args + ' < ' + refmac_data_file + ' > ' + refmac_log_file_name, 'r') # to screen too.
    refmac_finished = status.close()
# BL says: popen will always write a log file, so we check the size, is 0 if failed
    refmac_log_size = os.stat(refmac_log_file_name)[stat.ST_SIZE]

    print 'BL DEBUG:: refmac_log_size is: ', refmac_log_size
    
    if (refmac_log_size > 0) :
        run_loggraph(refmac_log_file_name)

    if (refmac_finished != None) : # refmac ran fail...
        print "Refmac Failed."
# BL says: lets check if we have refmac in PATH and that's why it failed
        refmacinpath(refmac_execfile)

    else : # refmac ran OK.

#        args = [mtz_out_filename, "FWT", "PHWT", "", 0, 0, 1, f_col, sig_f_col]

        # deal with R-free column
        recentre_status = recentre_on_read_pdb()
        imol = handle_read_draw_molecule(pdb_out_filename)
        if (recentre_status) :
            set_recentre_on_read_pdb(1)
        set_refmac_counter(imol, imol_refmac_count + 1)

        new_map_id = make_and_draw_map_with_refmac_params(mtz_out_filename, "FWT", "PHWT", "", 0, 0, 1, f_col, sig_f_col, r_free_col, 0) # FIXME c.f. apply

        if (swap_map_colours_post_refmac_p == 1) :
            swap_map_colours(imol_mtz_molecule, new_map_id)

        # difference map:
        if (show_diff_map_flag == 1) :
#            args = [mtz_out_filename, "DELFWT", "PHDELWT", "", 0, 1, 1, f_col, sig_f_col]
#            if (r_free_col != "") :
#                args.append(r_free_col);

            make_and_draw_map_with_refmac_params(mtz_out_filename, "DELFWT", "PHDELWT", "", 0, 1, 1, f_col, sig_f_col, r_free_col, 0) # FIXME needs apply

        
            
# return a boolean.        
def extra_params_include_weight_p(params_list):

   have_weight = params_list.count('WEIG')
   if (have_weight == 1):
     return True
   elif (have_weight == 0):
     return False
   else:
     print 'BL WARNING:: This shouldn\'t happen, we have more than one weight defined. God knows what refmac will do now...'
     return False


# return a list
def get_refmac_extra_params():
  global refmac_extra_params

  extras_file_name = "refmac-extra-params.txt"
  try:
     f = open(extras_file_name,'r')
  except IOError:
     print 'BL INFO:: we dont have refmac-extra-params file'
  else:
     refmac_extra_params = f.readlines()
     print 'BL INFO:: we have refmac-extra-params file and read lines'
     f.close()

#  print "BL INFO:: we added ", refmac_extra_params
     extra_params = refmac_extra_params
     return extra_params

def refmacinpath(refmac_execfile):
    import os, string
     # find executable
#    mode = kw.get("mode", os.P_WAIT)
    for path in string.split(os.environ["PATH"], os.pathsep):
        file = os.path.join(path, refmac_execfile)
        file_win = os.path.join(path, refmac_execfile) + ".exe"
        try:
            if (os.path.isfile(file) or os.path.isfile(file_win)):
             print 'But we could located refmac5 in: ' , file
             return True
            else:
             pass
        except os.error:
            return False
            pass
    print "..because we cannot find refmac5"

