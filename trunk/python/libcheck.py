# libcheck.py 
#
# Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
# Copyright 2005, 2006 by Bernhard Lohkamp
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

# BL to do:
# debug with given ccp4i_dir and lib file
# else if refmac doesnt exist
# aks paul about last line and nov
# dunno if always works, seems to have problems occasionally!
#

# This files provides an interface to libcheck

# Where @var{code} is @emph{e.g.} "3GP" and @var{ccp4i-project-dir}
# is an optional arg which is the directory to which the libcheck log
# files and pdb file will go.
#
# Return -2 on @var{code} is not a string
# Return -3 on libcheck failure
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
def monomer_molecule_from_3_let_code(code, dict_cif_libin, ccp4i_project_dir):

  import os, stat

  # return #t if log file has "has the minimal description" in it.
  # else return #f
  #
  def libcheck_minimal_qm(filename):
	
	# BL says: a quick fix for now, maybe we want to use something like
	# Paul's call_with_input_file funcn at some point
	ret = False
	try:
		f = open(filename, 'r')
		file = f.readlines()
		for line in file:
			if "has the minimal description" in line:
				ret = True
				break
		f.close()
	except:
		print "BL WARNING: couldnt open file ", filename
	return ret

  # main body

  if not isinstance(code,str):
      print "WARNING:: Oops code was not a string ", code
      return -2
  else:

   if (ccp4i_project_dir==""):
      dir_prefix = ccp4i_project_dir
   else:
      dir_prefix = ccp4i_project_dir + "/"

   code_str = str(code)

   if (len(dict_cif_libin)==0):
      libcheck_input = ["N","MON " + code_str,""]
   else:
# BL says:: Paul's code says FILE_CIF as argument, but I beleave this should be 
# FILE_L. FILE_CIF is for coordinate cif files and not dictionary cif file (which is
# what we have here, I think)
#      libcheck_input = ["N","FILE_CIF " + dict_cif_libin,"MON " + code_str,""]
      libcheck_input = ["N","FILE_L " + dict_cif_libin,"MON " + code_str,""]
   
   pdb_file_name = dir_prefix + "libcheck_" + code_str + ".pdb"
   cif_file_name = dir_prefix + "libcheck_" + code_str + ".cif"
   post_refmac_pdb_file_name = dir_prefix + "monomer-" + code_str + ".pdb"
   log_file_name = dir_prefix + "coot-libcheck-"  + code_str + ".log"
   refmac_input = ["MODE NEWENTRY","END"]
   refmac_log_file_name = dir_prefix + "coot-libcheck-refmac-" + code_str + ".log"
   refmac_command_line = " LIBIN " + cif_file_name + " XYZIN " + pdb_file_name + " XYZOUT " + post_refmac_pdb_file_name

   libcheck_exe = find_exe("libcheck", "CCP4_BIN", "PATH")

   if (libcheck_exe):
   # BL says: as with refmac we have to write a file with paramters 
   # to run libcheck. We can remove it later
      libcheck_input_file = "coot-libcheck-input.txt"
      input = file(libcheck_input_file,'w')
      for data in libcheck_input:
        input.write(data + '\n')
      input.close()
      status = os.popen(libcheck_exe + ' < ' + libcheck_input_file + ' > ' + log_file_name, 'r')
      libcheck_status = status.close()
      if (os.path.isfile(log_file_name) and not libcheck_status):
         # means we have a log file, 
         # dunno how else to check for status currently!?
         os.remove(libcheck_input_file)
         # we assume libcheck run ok
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

            # I have no clue what 'nov' means in guile! Just left it out...
            refmac_exe = find_exe("refmac5", "CCP4_BIN", "PATH")
            # BL says: as with refmac and before we have to 
            # write a file with paramters, to be removed later
            libcheck_refmac_input_file = "coot-libcheck-refmac-input.txt"
            input = file(libcheck_refmac_input_file,'w')
            for data in refmac_input:
               input.write(data + '\n')
            input.close()
#            print " run file", refmac_exe + refmac_command_line + ' < ' + libcheck_refmac_input_file + ' > ' + refmac_log_file_name
            # let's run refmac
            refmac_status = os.popen(refmac_exe + refmac_command_line + ' < ' + libcheck_refmac_input_file + ' > ' + refmac_log_file_name, 'r')
            refmac_finished = refmac_status.close()
            os.remove(libcheck_refmac_input_file)
            refmac_log_size = os.stat(refmac_log_file_name)[stat.ST_SIZE]
            # we assume if log > 0 refmac run
            if (refmac_log_size > 0 and not refmac_finished) :
               recentre_status = recentre_on_read_pdb()
               set_recentre_on_read_pdb(0)
               pdb_status = read_pdb(post_refmac_pdb_file_name)
               if (pdb_status > 0):
                  rc = rotation_centre()
                  mc = molecule_centre(pdb_status)
                  translate_molecule_by(pdb_status,rc[0]-mc[0],rc[1]-mc[1],rc[2]-mc[2])
                  set_recentre_on_read_pdb(recentre_status)
                  if libcheck_minimal_desc_status:
			read_cif_dictionary("libcheck.lib")
		  return pdb_status
            else:
              return -4     # refmac failed!?

      else:
        return -3  # libcheck failed !?
   else:
      print "BL WARNING:: Theres no libcheck hence not this function..."

#monomer_molecule_from_3_let_code("3GP","")
