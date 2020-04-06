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
#
def monomer_molecule_from_3_let_code(code, dict_cif_libin,
                                     ccp4i_project_dir = ""):

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
        print("BL WARNING: couldnt open file ", filename)
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
          print("BL WARNING:: could not rename file %s to %s" %(file_name, new_file_name))

  # the exe and log_file_name (and command_lines_args perhaps)
  # relative to dir (not the calling dir)
  # BL says:: maybe this should be in coot_utils.py
  def run_command_in_dir(run_dir, exe_name, command_line_args,
                         data_lines, log_file_name, to_screen_flag):
    import os

    print(("BL INFO:: run_command_in_dir dir:"
           "%s exe: %s cl_args: %s data_lines: %s log_file: %s"
           %(run_dir, exe_name, command_line_args, data_lines, log_file_name)))
    current = os.getcwd()
    norm_dir = os.path.normpath(run_dir)
    os.chdir(norm_dir)
    status = popen_command(exe_name, command_line_args, data_lines,
                           log_file_name, to_screen_flag)
    os.chdir(current)
    return status

  def log_file2text(log_file_name, dir):
    f = os.path.join(dir, log_file_name)
    return file2string(f)

  # This gets called in non-graphics mode too.
  #
  def libcheck_monomer_gui(dir_prefix, code_str, cif_file_name,
                           pdb_file_name, post_refmac_pdb_file_name):

    # print "================= debug: libcheck_monomer_gui: dir_prefix: %s   code_str: %s   cif_file_name: %s   pdb_file_name: %s   post_refmac_pdb_file_name: %s" %(dir_prefix, code_str, cif_file_name, pdb_file_name, post_refmac_pdb_file_name)

    if (len(dict_cif_libin) == 0):
      libcheck_input = ["N", "MON " + code_str[0:3], ""]
    else:
      # BL says:: Paul's code says FILE_CIF as argument, but I beleave this should be
      # FILE_L. FILE_CIF is for coordinate cif files and not dictionary cif file (which is
      # what we have here, I think)
      #      libcheck_input = ["N","FILE_CIF " + dict_cif_libin,"MON " + code_str,""]
      libcheck_input = ["N", "FILE_L " + dict_cif_libin,
                        "MON " + code_str, ""]

    log_file_name = "coot-libcheck-"  + code_str + ".log"
    log_file_name_in_dir = os.path.join(dir_prefix, log_file_name)
    refmac_input = ["MODE NEWENTRY", "END"]

    refmac_log_file_name = os.path.join(dir_prefix, "coot-libcheck-refmac-" + code_str + ".log")
    refmac_command_line = ["LIBIN", cif_file_name, "XYZIN", pdb_file_name,
                           "XYZOUT", post_refmac_pdb_file_name]

    move_aside(os.path.join(dir_prefix, "libcheck.lib"))
    #print "passing libcheck these data lines:", libcheck_input

    libcheck_exe_file = find_exe(libcheck_exe, "CBIN", "CCP4_BIN", "PATH")
    if (libcheck_exe_file):
      libstatus = run_command_in_dir(dir_prefix, libcheck_exe_file, [], libcheck_input, log_file_name, True)

      print("BL INFO:: libcheck status:", libstatus)

      if (not isNumber(libstatus)):

        log_text = log_file2text(log_file_name, dir_prefix)
        if (isinstance(log_text, str)):
          simple_text_dialog("Libcheck log", log_text, 200, 400)
          return -1
      else:
        if (libstatus != 0):

          log_text = log_file2text(log_file_name, dir_prefix)
          if isinstance(log_text, str):
            simple_text_dialog("Libcheck log", log_text, 400, 400)
            return -1

        else:
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
            print("libcheck failed to write the output cif file", cif_file_name)
            log_text = log_file2text(log_file_name, dir_prefix)
            if isinstance(log_text, str):
              simple_text_dialog("Libcheck log", log_text, 400, 400)
            return -1
          else:
            # OK, now let's run refmac:
            #
            libcheck_minimal_desc_status = libcheck_minimal_qm(log_file_name_in_dir)
            refmac_exe = find_exe("refmac5", "CBIN", "CCP4_BIN", "PATH")
            refmac_status = popen_command(refmac_exe, refmac_command_line, refmac_input, refmac_log_file_name)

            print("INFO:: libcheck-minimal? is ", libcheck_minimal_desc_status)

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
                libcheck_lib = os.path.join(dir_prefix, "libcheck.lib")
                #print "------------- about to copy file %s to %s in dir %s" %(libcheck_lib, cif_file_name, os.getcwd())
                if (os.path.isfile(libcheck_lib)):
                  shutil.copyfile(libcheck_lib, cif_file_name)
                return handle_libcheck_cif_and_pdb(cif_file_name,
                                                   pdb_file_name,
                                                   post_refmac_pdb_file_name)

  def handle_libcheck_cif_and_pdb(cif_file_name, pdb_file_name, post_refmac_pdb_file_name):

    #print "================= debug:: handle-libcheck-cif-and-pdb: cif-file-name: %s pdb-file-name: %s post-refmac-pdb-file-name: %s" %(cif_file_name, pdb_file_name, post_refmac_pdb_file_name)
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
    print("WARNING:: Oops code %s was not a string " %code)
    return -2
  else:

    # do the files exist already?  If so, just read them in.
    if (ccp4i_project_dir == ""):
      target_dir_name = "coot-ccp4"
      dir_name = get_directory(target_dir_name)
      dir_prefix = os.path.normpath(dir_name)
    else:
      dir_prefix = os.path.normpath(ccp4i_project_dir)
    #print "BL DEBUG:: dir prefix", dir_prefix

    code_str = str(code)

    pdb_file_name = os.path.join(dir_prefix, "libcheck_" + code_str + ".pdb")
    cif_file_name = os.path.join(dir_prefix, "libcheck_" + code_str + ".cif")
    post_refmac_pdb_file_name = os.path.join(dir_prefix, "monomer-" + code_str + ".pdb")

    if (os.path.isfile(post_refmac_pdb_file_name) and
        os.path.isfile(cif_file_name)):
      return handle_libcheck_cif_and_pdb(cif_file_name, pdb_file_name,
                                         post_refmac_pdb_file_name)

    else:
      libcheck_exe_file = find_exe(libcheck_exe, "CBIN", "CCP4_BIN", "PATH")
      if (not libcheck_exe_file):
        info_dialog("You need to setup CCP4 (specifically LIBCHECK) first.")
        return -2
      else:
        v = libcheck_monomer_gui(dir_prefix, code_str, cif_file_name,
                                 pdb_file_name, post_refmac_pdb_file_name)
        # somehow this can be none?! shouldnt
        if (not isNumber(v)):
          v = -1
        return v

#monomer_molecule_from_3_let_code("3GP","")
