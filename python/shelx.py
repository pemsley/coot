# shelx.py
# BL python adaption of shelx.scm by Paul Emsley
#
#    Copyright (C) 2005, 2006, 2007 by Bernhard Lohkamp
#    Copyright (C) 2007, 2008 by Bernhard Lohkamp, The University of York
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


# Provide functions to handle shelx fcf file (which I believe to be old-style
# cif rather than mmCIF).
# 
# So, in Coot, if it sees the extention of .fcf then it calls the
# fcf handler and the handler does the conversion, writes a new
# mmCIF file and then reads it.
# 

from types import *
import numbers

def handle_shelx_fcf_file_old(filename):
    import os

    filename = os.path.normpath(filename)

    # OK, so we have been passed a shelx fcf file.
    # We need to run this through the awk filter.
    # BL says: I guess in python it's easier to do it within python

    output_mmCIF_file_name = filename + ".cif"

    convert_shelx_fcf_to_cif(filename,output_mmCIF_file_name)

    ret = coot.auto_read_cif_data_with_phases(output_mmCIF_file_name)
    return ret
#    read_cif_data_with_phases(output_mmCIF_file_name)
#    coot.read_cif_data_with_phases_fo_fc(output_mmCIF_file_name)

def handle_shelx_fcf_file(filename):
    return coot.read_small_molecule_data_cif(filename)
    
def convert_shelx_fcf_to_cif(fcf_filename,cif_filename):

   import string, math

   fin = False
   fout = False

   intable = 0

   try:
    fin = open(fcf_filename,'r')
   except IOError:
    print("BL ERROR:: Cannot read ", fcf_filename)
   try:
    fout = file(cif_filename,"w")
   except IOError:
    print("BL ERROR:: Cannot write ", cif_filename)
   if (fin and fout):
    lines = fin.readlines()
    for line in lines:
       if "_shelx_title" in line:
          # dont write
          pass
       elif "_shelx_refln_list_code" in line:
          # dont write
          pass
       elif "_shelx_F_calc_maximum" in line:
          # dont write
          pass
       elif "_exptl_crystal_F_000" in line:
          # dont write
          pass
       elif "_reflns_d_resolution_high" in line:
          # dont write
          pass
       elif "_shelx_refln_list_code" in line:
          # dont write
          pass
       elif "_shelx_F_calc_maximum" in line:
          # dont write
          pass
       elif "_exptl_crystal_F_000" in line:
          # dont write
          pass
       elif "_reflns_d_resolution_high" in line:
          # dont write
          pass
       else:
        if "_refln_F_squared_meas" in line:
         line = line.replace("_refln_F_squared_meas"," _refln.F_meas")
        if "_refln_F_squared_sigma" in line:
         line = line.replace("_refln_F_squared_sigma"," _refln.F_meas_sigma")
        # symm equivs:
        if " _symmetry_equiv_pos_as_xyz" in line:
          line = line.replace(" _symmetry_equiv_pos_as_xyz","_symmetry_equiv.pos_as_xyz")
        if "_cell_length_a" in line:
          line = line.replace("_cell_length_a","_cell.length_a     ")
        if "_refln_F_squared_meas" in line:
         line = line.replace("_refln_F_squared_meas"," _refln.F_meas")
        if "_refln_F_squared_sigma" in line:
         line = line.replace("_refln_F_squared_sigma"," _refln.F_meas_sigma")
        # symm equivs:
        if "_symmetry_equiv_pos_as_xyz" in line:
          line = line.replace("_symmetry_equiv_pos_as_xyz","_symmetry_equiv.pos_as_xyz")
        if "_cell_length_a" in line:
          line = line.replace("_cell_length_a","_cell.length_a     ")
        if "_cell_length_b" in line:
          line = line.replace("_cell_length_b","_cell.length_b     ")
        if "_cell_length_c" in line:
          line = line.replace("_cell_length_c","_cell.length_c     ")
        if "_cell_angle_alpha" in line:
          line = line.replace("_cell_angle_alpha","_cell.angle_alpha  ")
        if "_cell_angle_beta" in line:
          line = line.replace("_cell_angle_beta","_cell.angle_beta   ")
        if "_cell_angle_gamma" in line:
          line = line.replace("_cell_angle_gamma","_cell.angle_gamma  ")
        if "_refln_index_h" in line:
          line = line.replace("_refln_index_h"," _refln.index_h")
        if "_cell_length_b" in line:
          line = line.replace("_cell_length_b","_cell.length_b     ")
        if "_cell_length_c" in line:
          line = line.replace("_cell_length_c","_cell.length_c     ")
        if "_cell_angle_alpha" in line:
          line = line.replace("_cell_angle_alpha","_cell.angle_alpha  ")
        if "_cell_angle_beta" in line:
          line = line.replace("_cell_angle_beta","_cell.angle_beta   ")
        if "_cell_angle_gamma" in line:
          line = line.replace("_cell_angle_gamma","_cell.angle_gamma  ")
        if "_refln_index_h" in line:
          line = line.replace("_refln_index_h"," _refln.index_h")
        if "_refln_index_k" in line:
          line = line.replace("_refln_index_k"," _refln.index_k")
        if "_refln_index_l" in line:
          line = line.replace("_refln_index_l"," _refln.index_l")
        if "_refln_F_calc" in line:
          line = line.replace("_refln_F_calc"," _refln.F_calc")
        if "_refln_phase_calc" in line:
          line = line.replace("_refln_phase_calc"," _refln.phase_calc")
          intable = 1
        if (intable == 1 and len(string.split(line))>5):
          col = string.split(line)
          col[3] = str(math.sqrt(float(col[3])))
          if (float(col[3]) > 0.0):
             col[4] = str(0.5 * (float(col[4]) + 0.0) / float(col[3]))
             line = string.join(col) + "\n"
        fout.write(line)

    fin.close()
    fout.close()
   else:
    print("BL ERROR:: IO error in conversion!")

global coot_shelxl_dir
coot_shelxl_dir = "coot-shelxl"

def remove_time_extensions(str):

	# e.g. remove -2007-04-27_1326.24 from 03srv164-2007-04-27_1326.24
	import re
	pattern = "-20[0-9][0-9]-[01][0-9]-[0-3][0-9]_[0-2][0-9][0-5][-0-9].[0-5][-0-9]"
	# we replace pattern with ''
	result = re.sub(pattern,"",str)
	return result

# BL says: shelxl_refine can have 3 arguments here:
#
# first: imol
#
# second: orig hkl file
#
# third: shelxh_flag , i.e. use shelxh instead of shelxl (in absenc of hkl
#        inputfile this has to be written as: shelxh_flag=YOURSETTING
#        True:  use shelxh
#        False: use shelxl (default)
#
def shelxl_refine(imol, hkl_file_in_maybe=False, shelxh_flag=False):

    func = lambda ins_file_name: coot.write_shelx_ins_file(imol, ins_file_name)
    shelxl_refine_inner(imol, hkl_file_in_maybe, func, shelxh_flag)

# shelxl refinement with input text
#
def shelxl_refine_primitive(imol, ins_text, hkl_file_in_maybe=False, shelxh_flag=False):

    func = lambda ins_file_name: coot_utils.save_string_to_file(ins_text, ins_file_name, True)
    shelxl_refine_inner(imol, hkl_file_in_maybe, func, shelxh_flag)

def shelxl_refine_inner(imol, hkl_file_in_maybe, func, shelxh_flag):

   import os
   import shutil
   global coot_shelxl_dir

   # First write out the imol-th molecule in shelx format.
   # What should the filename be?
   # Let's say that we read "abc.res"
   # We want to create "abc-coot.ins".  Ah, but what about the data file?
   # Urg.
   #
   # This, as it stands, doesn't properly deal with incrementing the
   # filenames.

   if (shelxh_flag):
       shelxl_exe = coot_utils.find_exe("shelxh", "PATH")
   else:
       shelxl_exe = coot_utils.find_exe("shelxl", "PATH")
   if (not shelxl_exe):
       print("WARNING:: can't find executable shelxl/h")
   else:
       orig_hkl_file = False
       if (hkl_file_in_maybe and os.path.isfile(hkl_file_in_maybe)):
           orig_hkl_file = hkl_file_in_maybe

       dir = os.path.join(os.path.abspath("."), coot_shelxl_dir)
       if (not os.path.isdir(dir)):
          try:
              os.mkdir(dir)
          except:
              pass
       if (not os.path.isdir(dir)):
           # We dont have shelxl dir and couldnt make it so exit!
           print("Failed to make shelxl directory ", dir)
       else:
         # run shelx
         stub = os.path.join(dir,(coot_utils.strip_path(coot_utils.strip_extension(coot.molecule_name(imol))) + "-" + coot_utils.unique_date_time_str()))
         # BL says: shelxl only excepts filename <= 80 char, so check it before
         if (len(stub) > 80):
             print("BL WARNING:: filename %s too long! Has %s characters, only 80 are allowed!" %(stub, len(stub)))
         else:
           ins_filename = stub + ".ins"
           res_filename = stub + ".res"
           lst_filename = stub + ".lst"
           fcf_filename = stub + ".fcf"
           log_filename = stub + ".log"
           hkl_filename = stub + ".hkl"

           # BL says: non-unix OS python cant make a link. grr!
           # So we have to make a copy of the passed hkl-file-in
           # if we cannot create a link
           # If we copied the file we maybe we should remove the
           # copied file after use since it's large?!
           #
           # helper function to either link or copy file
           # return True on sucess, False otherwise
           def make_symlink(src_file, dst_file):
               try:
                   print("make link from %s to %s" %(dst_file, src_file))
                   os.symlink(src_file, dst_file)
                   return True
               except:
                   # failed to create link, try so copy
                   try:
                       print("copy %s to %s" %(src_file, dst_file))
                       shutil.copyfile(src_file, dst_file)
                       return True
                   except:
                       return False

           # make a link to the passed hkl-file-in, if it was passed.
           # (if not, we presume that this file (or the link) exists
           # already. (Let's check that, and stop early if it doesn't
           # exist)

           if (orig_hkl_file):
               symlink_target = False
               if (len(orig_hkl_file) == 0):
                   # shouldnt happen!?
                   symlink_target = False
               elif (hkl_filename[0:1] == "/" or hkl_filename[0:1] == "\\"):
                   # BL nore: maybe we shoudl use coot_utils.slash_start_qm here too?!
                   symlink_target = orig_hkl_file
               else:
                   if (coot_utils.slash_start_qm(orig_hkl_file)):
                       symlink_target = orig_hkl_file
                   else:
                       symlink_target = "../" + orig_hkl_file
               if symlink_target:
                   make_symlink(symlink_target, hkl_filename)

           else:
               # no input hklin or it wasnt an accessible file
               #
               # hklin was not given, let's generate a filename
               # from the filename of the coordinates molecule
               # imol, trial names are derived from the name of
               # the coordinates molecule
               trial_file_stub     = coot_utils.strip_extension(coot.molecule_name(imol))
               trial_hkl_file_name = trial_file_stub + ".hkl"
               print("BL DEBUG:: looking for", trial_hkl_file_name)
               if (os.path.isfile(trial_hkl_file_name)):
                   print("INFO:: hkl file %s found" %trial_hkl_file_name)
                   # BL says: I dont understand the following but I gues
                   # I have tp believe Paul (at least for now)
                   # if trial-hkl-file-name and hkl-filename
                   # are in the same directory (e.g
                   # shelx-coot), then (strange as it may
                   # seem) we *don't* put the directory name
                   # in the "file that exists" (the first
                   # argument - ie. the target).  Bizarre but true.
                   hkl_filename_dir        = os.path.abspath(hkl_filename)
                   trial_hkl_file_name_dir = os.path.abspath(trial_hkl_file_name)
                   print("BL DEBUG:: comparing path", hkl_filename_dir, trial_hkl_file_name_dir)
                   if (hkl_filename_dir == trial_hkl_file_name_dir):
                       # same dirs
                       # so lets strip off the dir of the target:
                       new_target = coot_utils.strip_path(trial_hkl_file_name)
                       make_symlink(new_target, hkl_filename)
                   else:
                       # different dirs, keep the dirs:
                       make_symlink(trial_hkl_file_name, hkl_filename)
               else:
                   print("INFO:: hkl file %s does not exist!" %trial_hkl_file_name)

           if (not os.path.isfile(hkl_filename)):
              print("data (hkl) file %s not found - not running shelxl!" %hkl_filename)
           else:
              func(ins_filename)

              # running shelxl creates stub.res
              print("BL INFO:: Running shelxl as: ", shelxl_exe + " " + stub + " > " + log_filename)
              shelx_status = coot_utils.popen_command(shelxl_exe, [stub], [], log_filename)
              if (not shelx_status and
                  os.path.isfile(res_filename) and
                  os.path.isfile(fcf_filename)):
                  # it isn't a pdb file, but Coot knows what to do.
                  imol_res = coot.handle_read_draw_molecule_with_recentre(res_filename, 0)
                  handle_shelx_fcf_file(fcf_filename)
                  read_shelx_lst_file(lst_filename, imol_res)
              else:
                 print("BL WARNING:: shelxl didnt run ok! No output files found!")

#handle_shelx_fcf_file("test.fcf")
#shelxl_refine(0,"./shelx_tut/lyso-tet.hkl")

# do-shelx-lst-file
#
# ie. create a interesting-things GUI for split (and other things?)
# in a shelx .lst file.
#
def read_shelx_lst_file(file_name, imol):

    import os
    import operator

    # turn interesting-list into a GUI:
    #
    def gui_interesting_list(interesting_list):

        if (interesting_list):
           coot_gui.interesting_things_gui("Interesting Things from SHELX", interesting_list)
        else:
            print("INFO:: Nothing Interesting from LST file")
            coot.add_status_bar_text("Nothing Interesting from LST file")

    # chop off last char of string
    def chop_end(s):
        return s[:-1]

    # return last char of string:
    def last_char(s):
        return s[-1]

    # return a  chain resno inscode atom-name altconf list
    #
    def make_atom_parts(atom_id):

        list = []
        ls = atom_id.rsplit("_",1)
        if (len(ls)==2):

            # chain resno inscode atom-name altconf

            # The resno part can contain an alt conf specifier
            # e.g.: 34 or 34a. So let's try to make it a number
            # simply.  If that fails, then strip off the last char
            # and try again.  If *that* is successful, return the
            # resno and the last char is the altconf specifier.

            # BL says:: we try to make it to integer, otherwise exception!?
            resno = ls[1]
            try:
                    resno = int(resno)
                    chain_id = chain_id_for_shelxl_residue_number(imol, resno)
                    if (not chain_id):
                        print("couldn't find chain id for resno ", resno)
                        return ["", resno, "", ls[0], ""]
                    else:
                        return [chain_id, resno, "", ls[0], ""]

            except:	 	# handle the alt conf
                    s = ls[1]
                    p = chop_end(s)
                    alt_conf = last_char(s)
                    try:
                        n = int(p)
                    except:
                        return ["", 1, "", "blank", ""]	# failure
                    chain_id = chain_id_for_shelxl_residue_number(imol,n)
                    if (not chain_id):
                        print("couldn't find chain id for resno ", resno)
                        return ["", n, "", ls[0], alt_conf]
                    else:
                        return [chain_id, n, "", ls[0], alt_conf]

    # So we have a line in a table something like these:
    # 
    # "    2.7499    2.5580    0.1919    0.0400    DANG CN_2011 C_2011"
    # "                        1.5156    0.5000    FLAT O_1012 CA_1012 N_1013 CA_1013"
    # 
    # It seems that there are either 4 leading numbers or 2 leading
    # numbers.  Let's parse up the line here and return (list
    # observed target error sigma restraint) where observed and
    # target are #f if there number is missing (e.g. a "FLAT"
    # restraint).
    # 
    # There are 2 regimes we understand, either 
    # number number number number restraint-descr atom-1-desc atom-2-desc 
    # or 
    # number number restraint-descr atom-stuff [.. atom-stuff]
    # 
    # Actually here we will just make a string list of the
    # atom-stuffs (after the restraint-descr) and parse that
    # elsewhere (in do-gui) to find the atoms.  The actual
    # meaning of the atom-stuffs depends on the restraint-descr.
    # 
    # return #f on something we don't understand 
    # 
    def parse_dr_line(line):

        import string
        bits = line.split()

        if (len(bits) <=3):
            return False	# failed to understand
        else:
            try:
                n0 = float(bits[0])
            except:
                n0 = False
            try:
                n1 = float(bits[1])
            except:
                n1 = False
            try:
                n2 = float(bits[2])
            except:
                n2 = False
            try:
                n3 = float(bits[3])
            except:
                n3 = False

                if (not n0 and not n1):
                    return False	# failed to understand
                else:
                    if (n2 and n3):
                        # a 4 number line
                        if (len(bits)<=6):
                            return False 	# failed to understand
                        else:
                            lst = [n0, n1, n2, n3, bits[4], bits[5:]]
                            return lst
                    else:
                        # a 2 number line
                        if (len(bits)<=3):	# shouldnt actually happen as we checked earlier!!! Paul!?
                            return False	# failed to understand
                        else:
                            lst = [False, False, n0, n1, bits[2], bits[3:]]
                            return lst


    #
    def do_gui(disagreeable_restraints_list, interesting_list):
        #print "BL DEBUG:: DR: ", disagreeable_restraints_list

        dis_res_button_list = []
        for dr_list in disagreeable_restraints_list:

            restraint_type = dr_list[4]
            drl_dict = {"BUMP": 0, "DANG": 0, "FLAT": 0,
                        "SIMU": 1, "ISOR": 1}
            if (restraint_type in drl_dict):
                drl_index = drl_dict[restraint_type]
            else:
                drl_index = 0
            atom_parts = make_atom_parts(dr_list[5][drl_index])
            stats_string = "not sure what?"
            n2 = dr_list[2]
            n3 = dr_list[3]
            if (isinstance(n2, numbers.Number) and isinstance(n3, numbers.Number) and (n2 != 0)):
                z = n3 / abs(n2)
                stats_string = " " + str(n2) + " [Z=" + str(z) +"]"
            else:
                stats_string = ""
            if (drl_index == 0):
                rt = restraint_type + " "
            else:
                rt = restraint_type + " " + dr_list[5][0]
            button_label = "Disagreeable Restraint " + rt + " " + atom_parts[0] + " " + \
                           str(atom_parts[1]) + " " + atom_parts[3] + stats_string
            interesting_thing = [button_label, imol] + atom_parts
            dis_res_button_list.append(interesting_thing)

        gui_interesting_list(interesting_list + dis_res_button_list)


    # main body
    file_name = os.path.abspath(file_name)
    if not coot_utils.valid_model_molecule_qm(imol):
       print("WARNING:: Molecule number %i not valid" %imol)
    else:
       if not os.path.isfile(file_name):
          print("WARNING:: shelx lst file: %s does not exist" %file_name)
       else:
          fin = False
          interesting_list = []
          try:
             fin = open(file_name,'r')
          except:
             print("BL ERROR:: Cannot read ", file_name)
          if (fin):
             lines = fin.readlines()
             split_list = []
             disagreeable_restraints_list = []
             dr_count = 0
             for line in lines:
                 if ("may be split into" in line):
                   parts = line.split()
                   if len(parts)>6:
                     atom_id = parts[3]
                     # e.g; atom_id: CD1_1392
                     atom_parts = make_atom_parts(atom_id)
                     button_label = "Atom " + str(atom_parts[1]) + " " + atom_parts[3] + " may be split?"
                     split_list.append([button_label, imol] + atom_parts)
                     #split_list.append(atom_parts)
                 elif ("   Observed   Target    Error     Sigma     Restraint" in line):
                     dr_count = 1
                 elif (dr_count == 1):
                     disagreeable_restraints_list = []	# reset the disagreeable-restraints-list
                     dr_count = 2
                 elif (dr_count == 2):	# OK, in the disagreeable table
                     #print "BL DEBUG:: DR: ", line
                     dr_bits = parse_dr_line(line)
                     #print "BL DEBUG:: dr_bits: ", dr_bits,len(line)
                     if (type(dr_bits) is ListType):
                         disagreeable_restraints_list.append(dr_bits)
                     if (len(line) > 1):
                         dr_count = 2
                     else:
                         dr_count = 0
                 else:
                    pass
             fin.close()
             # print "BL DEBUG:: now run the gui with ", disagreeable_restraints_list, split_list
             do_gui(disagreeable_restraints_list, split_list)
          else:
            print("BL ERROR:: unknown error opening file ", file_name)

#read_shelx_lst_file("test.lst",0)

# Read a shelx project (i.e. the .res file, the .lst file and the
# .fcf.cif file (if it exists (else read the .fcf file (if it
# exists)))
#
# If the file-name has an extension of .lst .res .ins .fcf .log .hkl
# then strip that off before adding new extensions.
#
def read_shelx_project(file_name):
        import os

        file_name = os.path.abspath(file_name)
        extension = coot_utils.file_name_extension(file_name)
        if (extension == "lst"):
            file_stub = coot_utils.file_name_sans_extension(file_name)
        elif (extension == "ins"):
            file_stub = coot_utils.file_name_sans_extension(file_name)
        elif (extension == "insh"):
            file_stub = coot_utils.file_name_sans_extension(file_name)
        elif (extension == "log"):
            file_stub = coot_utils.file_name_sans_extension(file_name)
        elif (extension == "hkl"):
            file_stub = coot_utils.file_name_sans_extension(file_name)
        elif (extension == "res"):
            file_stub = coot_utils.file_name_sans_extension(file_name)
        elif (extension == "fcf"):
            file_stub = coot_utils.file_name_sans_extension(file_name)
        else:
            file_stub = file_name

        res_file_name = file_stub + ".res"
        lst_file_name = file_stub + ".lst"
        fcf_file_name = file_stub + ".fcf"
        fcf_cif_file_name = file_stub + ".fcf.cif"

        print(" file_name: ", file_name)
        print(" file_stub: ", file_stub)
        print(" extension: ", extension)

        if (os.path.isfile(res_file_name)):
                print("Read res file ", res_file_name)
                imol_res = coot.handle_read_draw_molecule_with_recentre(res_file_name, 0)
        else:
                print("  No res file ", res_file_name)
                imol_res = -1

        if (not coot_utils.valid_model_molecule_qm(imol_res)):
                print("WARNING:: Bad molecule from res file read.")
        else:
                if (os.path.isfile(fcf_cif_file_name)):
                    if (not os.path.isfile(fcf_file_name)):
                        print("   Read fcf-cif file ", fcf_cif_file_name)
                        coot.auto_read_cif_data_with_phases(fcf_cif_file_name)
                    else:
                        # OK both xxx.fcf and xxx.fcf.cif exist, we
                        # only want to read the xxx.fcf.cif if it is
                        # more recent than xxx.fcf (if it is not, we
                        # want to handle-shelx-fcf-file.
                        from stat import ST_MTIME
                        fcf_date = os.stat(fcf_file_name)[ST_MTIME]
                        fcf_cif_date = os.stat(fcf_cif_file_name)[ST_MTIME]
                        # print "    fcf_date %s fcf_cif_date %s" %(time.ctime(fcf_date), time.ctime(fcf_cif_date))
                        if (fcf_date < fcf_cif_date):
                            coot.auto_read_cif_data_with_phases(fcf_cif_file_name)
                        else:
                            handle_shelx_fcf_file(fcf_file_name)

                # xxx.fcf.cif did not exist:
                else:
                    if (os.path.isfile(fcf_file_name)):
                                print("   Read fcf file ", fcf_file_name)
                                handle_chelx_fcf_file(fcf_file_name)
                    else:
                        print("BL WARNING:: no fcf and or cif files found!")

                if (os.path.isfile(lst_file_name)):
                    print("  ::Read lst file ", lst_file_name)
                    read_shelx_lst_file(lst_file_name, imol_res)
                else:
                    print("   ::No lst file ", lst_file_name)
