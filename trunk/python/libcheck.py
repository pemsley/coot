# libcheck.py 
#
# Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
# Copyright 2005, 2006 by Bernhard Lohkamp
#    Copyright (C) 2007, 2008, 2009 by Bernhard Lohkamp, The University of York
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


# BL to do:
# debug with given ccp4i_dir and lib file
# else if refmac doesnt exist
#

# this is override-able by the user in their .coot file (for example).
libcheck_exe = "libcheck"

# This files provides an interface to libcheck

# Where @var{code} is @emph{e.g.} "3GP" and @var{ccp4i-project-dir}
# is an optional arg which is the directory to which the libcheck log
# files and pdb file will go.
#
# Return -2 on @var{code} is not a string
# Return -3 on libcheck failure
# Return -4 on refmac failure
# Return @var{imol} on success
# Return @code{handle-read-draw-molecule} error code on failure to 
#       read resultant pdb file
# 
# Actually, it would be nice to know if this code represented a full
# description - or a minimal one... perhaps we can parse the log
# file and call a pop-up that will tell us.
# 
# @var{dict_cif_libin should be a string}.  If it is "" then it is 
# ignored.  If it is not "" then it is used to create input to
# libcheck (not a command line argument) so that bespoke dictionary
# libraries can produce coords using "Get Monomer".
# 
# might not do all above mentioned things in python, yet
def monomer_molecule_from_3_let_code(code, dict_cif_libin, ccp4i_project_dir = ""):

  import os, stat, shutil

  # return True if log file has "has the minimal description" in it.
  # else return False
  #
  def libcheck_minimal_qm(filename):
	
	# BL says: a quick fix for now, maybe we want to use something like
	# Paul's call_with_input_file funcn at some point
	ret = False
	try:
		f = open(filename, 'r')
		lines = f.readlines()
		for line in lines:
			if "has the minimal description" in line:
				ret = True
				break
		f.close()
	except:
		print "BL WARNING: couldnt open file ", filename
	return ret

  # move file_name to some file name with a date extension
  #
  def move_aside(file_name):
    import time
    if (os.path.isfile(file_name)):
        # we use - and not : to avoid hazzle in windows
        extension = time.strftime('%d-%m-%Y-%H-%M-%S')
        new_file_name = file_name + "-" + extension
        try:
          os.rename(file_name, new_file_name)
        except:
          print "BL WARNING:: could not rename file %s to %s" %(file_name, new_file_name)
      
  def libcheck_monomer_gui(dir_prefix, code_str, cif_file_name,
                           pdb_file_name, post_refmac_pdb_file_name):

       if (len(dict_cif_libin) == 0):
         libcheck_input = ["N", "MON " + code_str[0:3], ""]
       else:
         # BL says:: Paul's code says FILE_CIF as argument, but I beleave this should be 
         # FILE_L. FILE_CIF is for coordinate cif files and not dictionary cif file (which is
         # what we have here, I think)
         #      libcheck_input = ["N","FILE_CIF " + dict_cif_libin,"MON " + code_str,""]
         libcheck_input = ["N", "FILE_L " + dict_cif_libin,
                           "MON " + code_str, ""]

       log_file_name = os.path.join(dir_prefix, "coot-libcheck-"  + code_str + ".log")
       refmac_input = ["MODE NEWENTRY", "END"]
       
       refmac_log_file_name = os.path.join(dir_prefix, "coot-libcheck-refmac-" + code_str + ".log")
       refmac_command_line = ["LIBIN", cif_file_name, "XYZIN", pdb_file_name,
                              "XYZOUT", post_refmac_pdb_file_name]

       move_aside("libcheck.lib")
       print "passing libcheck these data lines:", libcheck_input
       
       libcheck_exe_file = find_exe(libcheck_exe, "CCP4_BIN", "PATH")
       if (libcheck_exe_file):
         libstatus = popen_command(libcheck_exe_file, [], libcheck_input, log_file_name)

         print "BL INFO:: libcheck status:", libstatus

         if (not isNumber(libstatus)):
           return -3
         else:
           if (libstatus == 0):
             # we assume libcheck run ok
             #
             # But I now find that libcheck can run OK, but
             # not produce an output file (using dict .cif
             # file from PRODRG).
             #
             # So we first need to check that the output of
             # libcheck exists.
             # 
             if (not os.path.isfile(cif_file_name)):
               print "libcheck failed to write the output cif file."
               
             else:
               # OK, now let's run refmac:
               #
               libcheck_minimal_desc_status = libcheck_minimal_qm(log_file_name)
               refmac_exe = find_exe("refmac5", "CCP4_BIN", "PATH")
               refmac_status = popen_command(refmac_exe, refmac_command_line, refmac_input, refmac_log_file_name)
               
               print "DEBUG:: libcheck-minimal? returns ", libcheck_minimal_desc_status

               if (not isNumber(refmac_status)):
                 return -4 # refmac fails
               else:
                 if (not refmac_status == 0):
                   return -4 # refmac fails elsewhere
                 else:
                   # refmac run ok
                   #
                   # if there was a minimal description,
                   # we get the real cif file in
                   # libcheck.lib.
                   libcheck_lib = "libcheck.lib"
                   if (os.path.isfile(libcheck_lib)):
                     shutil.copyfile(libcheck_lib, cif_file_name)
                   handle_libcheck_cif_and_pdb(cif_file_name,
                                               pdb_file_name,
                                               post_refmac_pdb_file_name)
                   
  def handle_libcheck_cif_and_pdb(cif_file_name, pdb_file_name, post_refmac_pdb_file_name):

    if (os.path.isfile(post_refmac_pdb_file_name) and
        os.path.isfile(cif_file_name)):
      pdb_status = handle_read_draw_molecule_with_recentre(
        post_refmac_pdb_file_name, 0)
      if (valid_model_molecule_qm(pdb_status)):
        assign_hetatms(pdb_status)
        move_molecule_here(pdb_status)
        read_cif_dictionary(cif_file_name)
        return pdb_status  # return imol of the ligand
    return -1 # on fail
    

  # main body
  if not isinstance(code, str):
    print "WARNING:: Oops code was not a string ", code
    return -2
  else:

    # do the files exist already?  If so, just read them in.
    if (ccp4i_project_dir == ""):
      dir_prefix = ccp4i_project_dir
    else:
      dir_prefix = os.path.normpath(ccp4i_project_dir)
    print "BL DEBUG:: dir prefix", dir_prefix

    code_str = str(code)

    pdb_file_name = os.path.join(dir_prefix, "libcheck_" + code_str + ".pdb")
    cif_file_name = os.path.join(dir_prefix, "libcheck_" + code_str + ".cif")
    post_refmac_pdb_file_name = os.path.join(dir_prefix, "monomer-" + code_str + ".pdb")

    if (os.path.isfile(post_refmac_pdb_file_name) and
        os.path.isfile(cif_file_name)):
      handle_libcheck_cif_and_pdb(cif_file_name, pdb_file_name, post_refmac_pdb_file_name)
      
    else:
      libcheck_exe_file = find_exe(libcheck_exe, "CCP4_BIN", "PATH")
      if (not libcheck_exe_file):
        info_dialog("You need to setup CCP4 (specifically LIBCHECK) first.")
        return -2
      else:
        libcheck_monomer_gui(dir_prefix, code_str, cif_file_name,
                             pdb_file_name, post_refmac_pdb_file_name)

#monomer_molecule_from_3_let_code("3GP","")
