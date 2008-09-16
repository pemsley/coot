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

# BL says: Where do we get SG from? Otherwise clipper can't cif file!
# what shall we do???
# otherwise fine, if we include SG manually, e.g
# _symmetry.space_group_name_H-M  'C 2 2 21'

from types import *

def handle_shelx_fcf_file(filename):
    import os

    filename = os.path.normpath(filename)

    # OK, so we have been passed a shelx fcf file.
    # We need to run this through the awk filter.
    # BL says: I guess in python it's easier to do it within python

    output_mmCIF_file_name = filename + ".cif"

    convert_shelx_fcf_to_cif(filename,output_mmCIF_file_name)

    ret = auto_read_cif_data_with_phases(output_mmCIF_file_name)
    return ret
#    read_cif_data_with_phases(output_mmCIF_file_name)
#    read_cif_data_with_phases_fo_fc(output_mmCIF_file_name)



def convert_shelx_fcf_to_cif(fcf_filename,cif_filename):

   import string, math

   fin = False
   fout = False

   intable = 0

   try:
    fin = open(fcf_filename,'r')
   except IOError:
    print "BL ERROR:: Cannot read ", fcf_filename
   try:
    fout = file(cif_filename,"w")
   except IOError:
    print "BL ERROR:: Cannot write ", cif_filename
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
    print "BL ERROR:: IO error in conversion!"

def remove_time_extensions(str):

	# e.g. remove -2007-04-27_1326.24 from 03srv164-2007-04-27_1326.24
	import re
	pattern = "-20[0-9][0-9]-[01][0-9]-[0-2][0-9]_[0-2][0-9][0-5][-0-9].[0-5][-0-9]"
	# we replace pattern with ''
	result = re.sub(pattern,"",str)
	return result

# BL says: shelxl_refine can have 3 arguments here:
#
# first: imol
#
# second: orig hkl file
#
# third: use shelxh instead of shelxl (can be either argument 2 or 3, depending if hkl is there)
#        is flag: 1 is on -> use shelxh
#
def shelxl_refine(*args):

   import os
   import shutil
   import operator
   shelxh_flag = False

   if (len(args)<=3 and len(args)>=0): 
    if (len(args)==1):
       orig_hkl_file = False
    elif (len(args)==2):
# BL says: we dont check if axis is string, yet at least not directly
       if (operator.isNumberType(args[1])):
          orig_hkl_file = False
          if (args[1]==1):
             shelxh_flag = True
          else:
             print "BL WARNING:: unknown second argument %s for shelxl_refine!" %args[1]
       else:
         orig_hkl_file = args[1]
    else:
       #have 3 args
       orig_hkl_file = args[1]
       if (args[2]==1):
          shelxh_flag = True
       else:
          print "BL WARNING:: unknow second argument %s for shelxl_refine!" %args[1]
    imol = args[0]

    coot_shelxl_dir = "coot-shelxl"

    # First write out the imol-th molecule in shelx format.
    # What should the filename be? 
    # Let's say that we read "abc.res"
    # We want to create "abc-coot.ins".  Ah, but what about the data file?
    # Urg.
    # 
    # This, as it stands, doesn't properly deal with incrementing the
    # filenames.

    if (shelxh_flag):
       shelxl_exe = find_exe("shelxh", "PATH")
    else:
       shelxl_exe = find_exe("shelxl", "PATH")
    if (shelxl_exe):
       dir = os.path.join(os.path.abspath("."),coot_shelxl_dir)
       if not os.path.isdir(dir):
          try: os.mkdir(dir)
          except: print "Failed to make shelxl directory ", dir
       if os.path.isdir(dir):
         # run shelx
         stub = os.path.join(dir,(strip_path(strip_extension(molecule_name(imol))) + "-" + unique_date_time_str()))
         # BL says: shelxl only excepts filename <= 80 char, so check it before
         if len(stub)<=80:
           ins_filename = stub + ".ins"
           res_filename = stub + ".res"
           fcf_filename = stub + ".fcf"
           log_filename = stub + ".log"
           hkl_filename = stub + ".hkl"
           # create a hkl masta file name based on imol name
           dir_stub, tmp = os.path.split(dir)
           print "BL DEBUG:: dir_stub, tmp, dir", dir_stub, tmp, dir
           hkl_master_file_name = os.path.join(dir_stub ,(strip_path(strip_extension(molecule_name(imol))) + ".hkl"))
           print "BL DEBUG:: stub and masta is", stub, hkl_master_file_name

           # BL says: in windows python cant make a link. grr!
           # So we have to make a copy of the passed hkl-file-in
           # we follow Paul's idea 
           # But use copy instead of link! Maybe we should remove 
           # copied file after use since it's large...
           # make a link to the passed hkl-file-in, if it was passed.
           # (if not, we presume that this file (or the link) exists
           # already. (Let's check that, and stop early if it doesn't
           # exist)

           if (orig_hkl_file):
              if len(orig_hkl_file)==0:
                 symlink_target = False
              else:
                 if os.path.isfile(orig_hkl_file):
                    symlink_target = orig_hkl_file
           else:
              if os.path.isfile(hkl_master_file_name):
                 symlink_target = hkl_master_file_name
              else:
                 symlink_target = False
           if symlink_target:
              print "copy %s to %s" %(symlink_target, hkl_filename)
              shutil.copyfile(symlink_target,hkl_filename)
           # check if masta file is there already 
           if not os.path.isfile(hkl_filename):
              print "data (hkl) file %s not found - not running shelxl!" %hkl_filename
           else:
              write_shelx_ins_file(imol,ins_filename)
              # running shelxl creates stub.res
              print "BL INFO:: Running shelxl as: ",shelxl_exe + " " + stub + " > " + log_filename
              shelx_status = popen_command(shelxl_exe, [stub], [], log_filename)
              #status = os.popen(shelxl_exe + " " + stub + " > " + log_filename,'r')
              #shelx_status = status.close()
              if (not shelx_status and os.path.isfile(res_filename) and os.path.isfile(fcf_filename)):
                 read_pdb(res_filename)
                 handle_shelx_fcf_file(fcf_filename)
              else: 
                 print "BL WARNING:: shelxl didnt run ok! No output files found!"
              # let's move the hkl file back to save space
# check this
              # shutil.move(hkl_filename,hkl_master_file_name)
         else:
            print "BL WARNING:: filename %s too long! Has %s characters, only 80 are allowed!" %(stub, len(stub))
       else:
           # We dont have shelxl dir and couldnt make it so exit!
           print "BL WARNING:: couldnt make a shelxl directory, so wont run shelxl!"
       
    else:
       print "coot warning: can't find executable shelxl"
   else:
       print "Wrong number of arguments. shelxl_refine takes 1, 2 or 3 arguments!"


#handle_shelx_fcf_file("test.fcf")
#shelxl_refine(0,"./shelx_tut/lyso-tet.hkl")

# do-shelx-lst-file
#
# ie. create a interesting-things GUI for split (and other things?)
# in a shelx .lst file.
# 
def read_shelx_lst_file(file_name,imol):

    import os

    # turn interesting-list into a GUI:
    # 
    def gui_interesting_list(interesting_list):

#        print "debug: interesting_list:", interesting_list
        if (interesting_list):
           interesting_things_gui("Interesting Things from SHELX",interesting_list)
	else:
	   print "INFO:: Nothing Interesting from LST file" 
	   add_status_bar_text("Nothing Interesting from LST file")

    # chop off last char of string
    def chop_end(s):
	res = s[:-1]
	return res

    # return last char of string:
    def last_char(s):
	res = s[-1]
	return res

    # return a  chain resno inscode atom-name altconf list
    #
    def make_atom_parts(atom_id):

	print "BL DEBUG:: in make_atom_parts"
        list = []
        ls = atom_id.rsplit("_",1)
	print "BL DEBUG:: ls", ls
        if (len(ls)==2):
	
	      # chain resno inscode atom-name altconf
	      # 
	      # The resno part can contain an alt conf specifier
	      # e.g.: 34 or 34a. So let's try to make it a number
	      # simply.  If that fails, then strip off the last char
	      # and try again.  If *that* is successful, return the
	      # resno and the last char is the altconf specifier.
	      # 
	      # BL says:: we try to make it to integer, otherwise exception!?
		resno = ls[1]
		print "BL DEBUG:: resno", resno
		try:
		   	resno = int(resno)
			print "BL DEBUG:: chain_id not avaliable yet, ignore for now"
			#chain_id = chain_id_for_shelxl_residue_number(imol, resno)
			chain_id = False
			if (not chain_id):
				print "couldn't find chain id for resno ", resno
				print "BL DEBUG:: returning", ["", resno, "", ls[0], ""]
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
				return ["",1,"","blank",""]	# failure
			chain_id = chain_id_for_shelxl_residue_number(imol,n)
			if (not chain_id):
				print "couldn't find chain id for resno ", resno
				return ["", n, "", ls[0], alt_conf]
			else:
				return [chain_id,n,"",ls[0],alt_conf]

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
					lst = [n0, n1, n2, n3, bits[4], string.join(bits[5:len(bits)])]
					return lst
			else:
				# a 2 number line
				if (len(bits)<=3):	# shouldnt actually happen as we checked earlier!!! Paul!?
					return False	# failed to understand
				else:
					lst = [False, False, n0, n1, bits[2], string.join(bits[3:len(bits)])]


    #
    def do_gui(disagreeable_restraints_list, interesting_list):
	#print "DR: ", disagreeable_restraints_list

	dis_res_button_list = []
	for dr_list in disagreeable_restraints_list:

		if (not dr_list):
			dis_res_button_list.append(interesting_list)
		else:
			restraint_type = dr_list[0][4]
			drl_dict = {"BUMP": 0, "DANG": 0, "FLAT": 0,
				    "SIMU": 1, "ISOR": 1}
			if (drl_dict.has_key(restraint_type)):
				drl_index = drl_dict[restraint_type]
			else:
				drl_index = 0
			atom_parts = make_atom_parts(dr_list[0][5][drl_index])
			stats_string = "not sure what?"
			n2 = dr_list[0][2]
			n3 = dr_list[0][3]
			if (operator.isNumberType(n2) and operator.isNumberType(n3) and (n2 <> 0)):
				z = n3 / abs(n2)
				stats_string = " " + str(n2) + " [Z=" + str(z) +"]"
			else:
				stats_string = ""
			if drl_index == 0:
				rt = retraint_type + " "
			else:
				rt = restraint_type + " " + dr_list[0][5][0]
			button_label = "Disagreeable Restraint " + rt + " " + atom_parts[0] + " " + str(atom_parts[1]) + " " + atom_parts[3] + stats_string
			interesting_thing = [button_label, [imol,atom_parts]]
			dis_res_button_list.append(interesting_thing)

	gui_interesting_list(interesting_list.append(dis_res_button_list))


    # main body
    file_name = os.path.abspath(file_name)
    if not valid_model_molecule_qm(imol):
       print "WARNING:: Molecule number %i not valid" %imol
    else:
       if not os.path.isfile(file_name):
          print "WARNING:: shelx lst file: %s does not exist" %file_name
       else:
          fin = False
          interesting_list = []
          try:
             fin = open(file_name,'r')
          except:
             print "BL ERROR:: Cannot read ", file_name
          if (fin):
             lines = fin.readlines()
	     split_list = []
 	     disagreeable_restraints_list = []
	     dr_count = 0
             for line in lines:
                 if ("may be split into" in line):
		   print "BL DEBUG:: found split",
                   parts = line.split()
		   print "BL DEBUG:: this is the parts", parts
                   if len(parts)>6:
                     atom_id = parts[3]
                     # e.g; atom_id: CD1_1392
		     print "BL DEBUG:: now make the atom parts"
                     atom_parts = make_atom_parts(atom_id)
		     buton_label = "Atom " + str(atom_parts[1]) + " " + atom_parts[3] + " may be split?"
		     print "BL DEBUG:: buton_label is :", buton_label
		     split_list.append([buton_label,imol])
		     split_list.append(atom_parts)
	         elif ("   Observed   Target    Error     Sigma     Restraint" in line):
			dr_count = 1
		 elif (dr_count == 1):
			disagreeable_retraints_list = []	# reset the disagreeable-restraints-list
			dr_count = 2
		 elif (dr_count == 2):	# OK, in the disagreeable table
			#print "DR: ", line
			dr_bits = parse_dr_line(line)
			#print "dr_bits: ", dr_bits,len(line)
			if (type(dr_bits) is ListType):
				disagreeable_restraints_list.append(dr_bits)
			if (len(line) > 1):
				dr_count = 2
			else:
				dr_count = 0
                 else:
                    pass
	     print "BL DEBUG:: closing ", fin
             fin.close()
	     #print "BL DEBUG:: now run the gui with ", disagreeable_restraints_list, split_list
	     print "BL DEBUG:: now run the gui with split", split_list
             do_gui(disagreeable_restraints_list, split_list)
          else:
            print "BL ERROR:: unknown error opening file ", file_name

#read_shelx_lst_file("test.lst",0)

# Read a shelx project (i.e. the .res file, the .lst file and the
# .fcf.cif file (if it exists (else read the .fcf file (if it
# exists)))
#
# If the file-name has an extension of .lst .res .ins .fcf .log .hkl
# then strip that off before adding new extensions.
#
def read_shelx_project(file_name):
	import os, time

	file_name = os.path.abspath(file_name)
	extension = file_name_extension(file_name)
	if (extension == "lst"):
		file_stub = file_name_sans_extension(file_name)
	elif (extension == "ins"):
		file_stub = file_name_sans_extension(file_name)
	elif (extension == "insh"):
		file_stub = file_name_sans_extension(file_name)
	elif (extension == "log"):
		file_stub = file_name_sans_extension(file_name)
	elif (extension == "hkl"):
		file_stub = file_name_sans_extension(file_name)
	elif (extension == "res"):
		file_stub = file_name_sans_extension(file_name)
	elif (extension == "fcf"):
		file_stub = file_name_sans_extension(file_name)
	else:
		file_stub = file_name

	res_file_name = file_stub + ".res"
	lst_file_name = file_stub + ".lst"
	fcf_file_name = file_stub + ".fcf"
	fcf_cif_file_name = file_stub + ".fcf.cif"

	print " file_name: ", file_name
	print " file_stub: ", file_stub
	print " extension: ", extension

	if (os.path.isfile(res_file_name)):
		print "Read res file ", res_file_name
		read_shelx_ins_file(res_file_name)
		imol_res = handle_read_draw_molecule_with-recentre(res_file_name,0)
	else:
		print "  No res file ", res_file_name
		imol_res = -1

	if (not valid_model_molecule_qm(imol_res)):
		print "WARNING:: Bad molecule from res file read."
	else:
		if (os.path.isfile(fcf_cif_file_name)):
			if (not os.path.isfile(fcf_file_name)):
				print "   Read fcf-cif file ", fcf_cif_file_name
				auto_read_cif_data_with_phases(fcf_cif_file_name)
			else:
			# OK both xxx.fcf and xxx.fcf.cif exist, we
			# only want to read the xxx.fcf.cif if it is
			# more recent than xxx.fcf (if it is not, we
			# want to handle-shelx-fcf-file.
				fcf_date = os.stat.st_mtime(fcf_file_name)
				fcf_cif_date = os.stat.st_mtime(fcf_cif_file_name)
				# print "    fcf_date %s fcf_cif_date %s" %(time.ctime(fcf_date), time.ctime(fcf_cif_date))
				if (fcf_date < fcf_cif_date):
					auto_read_cif_data_with_phases(fcf_cif_file_name)
				else:
					handle_shelx_fcf_file(fcf_file_name)

		# xxx.fcf.cif did not exist:
		else:
			if (os.path.isfile(fcf_file_name)):
				print "   Read fcf file ", fcf_file_name
				handle_chelx_fcf_file(fcf_file_name)
			else:
				print "BL WARNING:: no fcf and or cif files found!"

		if (os.path.isfile(lst_file_name)):
			print "  ::Read lst file ", lst_file_name
			read_shelx_lst_file(lst_file_name, imol_res)
		else:
			print "   ::No lst file ", lst_file_name
