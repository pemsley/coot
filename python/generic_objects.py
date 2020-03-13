# generic-objects.py
#
# Copyright 2006 by Paul Emsley, The University of York
# Copyright 2006 by Bernhard Lohkamp
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

global probe_command
global reduce_command
global reduce_molecule_updates_current
global interactive_probe_is_OK_qm

# gets set the first time we run interactive_probe.  Takes values
# 'unset' (initial value) 'yes' and 'no')
#
interactive_probe_is_OK_qm = 'unset'

# set reduce  and probe command full path here if you wish
# (or have a look in group_settings!! These overwrite whatever is here), e.g.
# reduce_command = "/home/bernhard/bin/reduce"

# BL says:: just to be consistent with Paul's names
# map to scheme names:
# deftexi generic_object_is_displayed_qm
generic_object_is_displayed_qm = generic_object_is_displayed_p

# map to scheme names:
# deftexi is_closed_generic_object_qm 
is_closed_generic_object_qm = is_closed_generic_object_p

# return a new generic object number for the given object obj-name.
# If there is no such object with name obj-name, then create a new
# one.  (Maybe you want the generic object to be cleaned if it exists
# before returning, this function does not do that).
# 
def generic_object_with_name(obj_name):

  t = generic_object_index(obj_name)
  if (t == -1):
    t = new_generic_object_number(obj_name)
  return t


# display a GUI for generic objects
#
def generic_objects_gui():

    pygtk_flag = False
    try:
      import pygtk
      pygtk.require("2.0")
      import gtk, pango
      pygtk_flag = True
    except:
      print "BL WARNING:: no pygtk2. Function wont work!!!"
	  
    if (pygtk_flag and using_gui()):
      # Now we run the gui
      def delete_event(*args):
	# BL says: first we shall close the generic objects
	# or not
	#for generic_object_number in range(n_objects):
	#	set_display_generic_object(generic_object_number, 0)
	gen_window.destroy()
	return False


      def check_button_callback(widget, generic_object_number):
	button_state = widget.get_active()
	object_state = generic_object_is_displayed_qm(generic_object_number)
	if ((button_state == True) and (object_state == 0)):
	  set_display_generic_object(generic_object_number,1)
	if ((button_state == False) and (object_state == 1)):
	  set_display_generic_object(generic_object_number,0)

      def all_check_button_callback(widget):
	show_all = widget.get_active()
	for check_button in open_generic_objects:
	  if (show_all):
	    check_button.set_active(True)
	  else:
	    check_button.set_active(False)

      n_objects = number_of_generic_objects()
	
      if (n_objects > 0):
	gen_window = gtk.Window(gtk.WINDOW_TOPLEVEL)
	gen_window.set_title("Generic objects")
	vbox = gtk.VBox(False, 0)
	open_generic_objects = []
	no_active_generic_objects = 0

	for generic_object_number in range(n_objects):

	  print "INFO:: generic object attributes: ", \
		generic_object_number, \
		generic_object_name(generic_object_number), \
		is_closed_generic_object_qm(generic_object_number)

	  if (is_closed_generic_object_qm(generic_object_number) == 0):
	    name = generic_object_name(generic_object_number)

	    if (name):
	      label = str(generic_object_number) + "  " + name
	      frame = gtk.Frame(label=None)
	      check_button = gtk.CheckButton(label)

		# this callback gets called by the
		# gtk-toggle-button-set-active just below,
		# which is why we need the state and active
		# test.
	      check_button.connect("toggled", 
		     check_button_callback, generic_object_number)
	      current_state = generic_object_is_displayed_qm(generic_object_number)
	      if (current_state == 1):
		check_button.set_active(True)
		no_active_generic_objects += 1

	      vbox.add(frame)
	      frame.add(check_button)
	      frame.show()
	      check_button.show()
	      open_generic_objects.append(check_button)

	if (len(open_generic_objects) > 1):
	  hsep = gtk.HSeparator()
	  label = "Show/hide all"
	  frame = gtk.Frame(label=None)
	  check_button = gtk.CheckButton(label)
	  check_button.connect("toggled", all_check_button_callback)
	  if (len(open_generic_objects) == no_active_generic_objects):
	    check_button.set_active(True)
	  vbox.add(hsep)
	  hsep.show()
	  vbox.add(frame)
	  frame.add(check_button)
	  frame.show()
	  check_button.show()

	gen_window.connect("delete_event", delete_event)

	vbox.show()
	gen_window.add(vbox)
	gen_window.set_border_width(10)
	gen_window.show()


# return status
#
def reduce_on_pdb_file(imol, pdb_in, pdb_out):
    return reduce_on_pdb_file_generic(imol, "build", pdb_in, pdb_out)

# return status
#
def reduce_on_pdb_file_no_flip(imol, pdb_in, pdb_out):
    return reduce_on_pdb_file_generic(imol, "no-flip", pdb_in, pdb_out)

# return status
#
def reduce_on_pdb_file_generic(imol, no_flip_or_build, pdb_in, pdb_out):

  global reduce_command
  
  print "running reduce on", pdb_in
  # could try find_exe too!
  if not command_in_path_qm(reduce_command):
    print "command for reduce %s is not found in path" %reduce_command
  else:
    # need full path to find het dict
    full_reduce_command = find_exe(reduce_command, "PATH")

    nshl = non_standard_residue_names(imol)
    ext = "-".join(nshl)
    reduce_het_dict_file_name = "coot-molprobity/reduce-het-dict-" + ext + ".txt"
    write_reduce_het_dict(imol, reduce_het_dict_file_name)
    
    dict_args = []
    if not set_reduce_het_dict():
      # now should use the and build the het dic!?
      dict_args = ["-DB", reduce_het_dict_file_name]

    mode = "-build"
    if no_flip_or_build == "no-flip":
        mode = "-NOFLIP"

    print "======= reduce_on_pdb_file: command %s args %s with pdb_out: %s" \
          %(reduce_command,
            [mode, pdb_in] + dict_args,
            pdb_out)
    status = popen_command(reduce_command,
                           [mode, pdb_in] + dict_args,
                           [],
                           pdb_out)

    print "======== status", status

    # 20130107 returning (status == 0) is problematic because sometimes/often reduce exits
    # with status 1 but it has worked OK!
    #
    # return status == 0
    return True

  
global old_pdb_style
old_pdb_style = False
  
reduce_molecule_updates_current = False
       
# run molprobity (well reduce and probe) to make generic objects (and
# display the generic objects gui)
#
def probe(imol):
  import os
  global reduce_command, probe_command
  global old_pdb_style
    
  if is_valid_model_molecule(imol):

    if not (os.path.isfile(reduce_command)):
      reduce_command = find_exe("reduce", "PATH", "CBIN", "CCP4_BIN")
    # we need to check if probe_command is defined too
    if not(os.path.isfile(probe_command)):
      probe_command = find_exe("probe", "PATH", "CBIN", "CCP4_BIN")
    make_directory_maybe("coot-molprobity")
    mol_pdb_file = "coot-molprobity/for-reduce.pdb"
    reduce_out_pdb_file = "coot-molprobity/reduced.pdb"
    reduce_het_dict_file_name = "coot-molprobity/reduce-het-dict.txt"
    write_pdb_file(imol, mol_pdb_file)
    write_reduce_het_dict(imol, reduce_het_dict_file_name)
    if not reduce_command:
      # couldnt find reduce
      print "BL WARNING:: Could not locate the program reduce!! Please check if installed!"
    else:

      dict_args = []
      if not set_reduce_het_dict():
        dict_args = ["-DB", reduce_het_dict_file_name]

      if old_pdb_style:
        # old
        arg_list = ["-build", "-oldpdb", mol_pdb_file]
      else:
        # modern 
        arg_list = ["-build", mol_pdb_file]

      arg_list += dict_args

      print "BL INFO:: running reduce: %s %s and ouptut to: %s" \
            %(reduce_command , arg_list, reduce_out_pdb_file)
      print "BL INFO:: running reduce: REDUCE_HET_DICT env var:", os.getenv('REDUCE_HET_DICT')
      
      reduce_status = popen_command(reduce_command,
                                    arg_list,
                                    [],
                                    reduce_out_pdb_file)
      
      # dont check for status as meaningless...
      if not probe_command:
        # couldnt find probe
        print "BL WARNING:: Could not locate the program probe!! Please check if installed!"
      else:
        probe_name_stub = strip_extension(strip_path(molecule_name(imol)))
        probe_pdb_in = "coot-molprobity/" + probe_name_stub + "-with-H.pdb"
        probe_out = "coot-molprobity/probe-dots.out"

        prepare_file_for_probe(reduce_out_pdb_file, probe_pdb_in)

        probe_status = popen_command(probe_command,
                                     ["-u", "-mc", "ALL", probe_pdb_in],
                                     [],
                                     probe_out)

        if (probe_status):
          print "BL WARNING:: probe failed, cannot continue!"
        else:
          # by default, we don't want to click on the
          # imol-probe molecule (I think :-)
          recentre_status = recentre_on_read_pdb()
          novalue = set_recentre_on_read_pdb(0)
          if (reduce_molecule_updates_current):
            print "======= update molecule ======="
            imol_probe = clear_and_update_model_molecule_from_file(imol, probe_pdb_in)
          else:
            print "======= read new pdb file ======="
            imol_probe = read_pdb(probe_pdb_in)

          if recentre_status == 1:
            set_recentre_on_read_pdb(1)

          # show the GUI for USER MODS
          if using_gui():
            user_mods_gui(imol_probe, reduce_out_pdb_file)

          # toggle_active_mol(imol_probe) let's not do
          # that actually.  I no longer think that the
          # new probe molecule should not be clickable
          # when it is initally displayed (that plus
          # there is some active/displayed logic problem
          # for the molecules, which means that after
          # several probes, the wrong molecule is
          # getting refined).

          handle_read_draw_probe_dots_unformatted(probe_out, imol_probe, 2)
          generic_objects_gui()
          graphics_draw()


# Write the connectivity for the non-standard (non-water) residues in
# the given molecule for which we have the dictionary.
#
# Don't return anything interesting.  
#
def write_reduce_het_dict(imol, reduce_het_dict_file_name):

  import shutil
  con_file_names = []
  for res_name in non_standard_residue_names(imol):
    f_name = "coot-molprobity/conn-" + res_name + ".txt"
    status = write_connectivity(res_name, f_name)
    if (status == 1):
      con_file_names.append(f_name)
  if con_file_names:
    fin = open(reduce_het_dict_file_name, 'w')
    for file_name in con_file_names:
      shutil.copyfileobj(open(file_name, 'rb'), fin)
    fin.close()
    
      
# Prepare file for probe, i.e. remove 'USER' from file
def prepare_file_for_probe(file_in, file_out):
    
    try:
      fin = open(file_in,'r')
    except IOError:
      print "BL WARNING:: Cannot read ", file_in
    try:
      fout = file(file_out,'w')
    except IOError:
      print "BL WARNING:: Cannot write ", file_out
    if (fin and fout):
       lines = fin.readlines()
       for line in lines:
           if (not "USER" in line):
               fout.write(line)
       fin.close()
       fout.close()

 
# run "probe" interactively, 
# which in the current implementation, means that this function 
# can run during a edit-chi angles manipulation, or after
# a real space refine zone.
# 
# Thus function presumes that there are 2 pdb files in the current 
# directory, one of which is the reference pdb file and the other
# is a pdb file containing the tmp/moving atom set.
# 
# The function takes arguments for the centre of the probe dots
# and the radius of the probe dots sphere.  The chain id and 
# residue number are also needed to pass as arguments to probe.
#
def interactive_probe(x_cen, y_cen, z_cen, radius, chain_id, res_no):

    import os, string
    global probe_command
    global interactive_probe_is_OK_qm

    probe_pdb_in_1 = "molprobity-tmp-reference-file.pdb"
    probe_pdb_in_2 = "molprobity-tmp-moving-file.pdb"
    probe_out = "coot-molprobity/molprobity-tmp-probe-dots.out"
    chain_str = ""
    if (not chain_id == ""):
      chain_str = "chain"
    atom_sel = "(file1 within " + str(radius) + " of " \
	       + str(x_cen) + ", " \
	       + str(y_cen) + ", " \
	       + str(z_cen) + ", " \
	       + "not water not (" + chain_str + string.lower(chain_id) + " " \
	       + str(res_no) + ")),file2"

    # no longer use std-bonds
    print "probe command", probe_command, \
          ["-mc", "-u", "-quiet", "-drop", "-both",
          atom_sel, "file2", probe_pdb_in_1, probe_pdb_in_2]

    # if unset, then set it.
    if (interactive_probe_is_OK_qm == 'unset'):
      if (command_in_path_qm(probe_command)):
        interactive_probe_is_OK_qm = 'yes'
      else:
        interactive_probe_is_OK_qm = 'no'
        
    if (interactive_probe_is_OK_qm == 'yes'):
       status = popen_command(probe_command,
			      ["-mc", "-u", "-quiet", "-drop", "-both",
			       atom_sel, "file2",
			       probe_pdb_in_1, probe_pdb_in_2],
			      [],
			      probe_out)

       # don't show the gui, so the imol is not needed/dummy.
       handle_read_draw_probe_dots_unformatted(probe_out, 0, 0)
       graphics_draw()

#
#
def get_probe_dots_from(pdb_file_name, point, radius):

  global probe_command
  global interactive_probe_is_OK_qm
  import os
  # if unset, then set it, try to make dir too
  if (interactive_probe_is_OK_qm == 'unset'):
    if (not command_in_path_qm(probe_command)):
      interactive_probe_is_OK_qm = 'no'
    else:
      dir_name = get_directory("coot-molprobity")
      if (dir_name):    # OK, we had it or made it
        interactive_probe_is_OK_qm = 'yes'
      else:
        interactive_probe_is_OK_qm = 'no'

  if (interactive_probe_is_OK_qm == 'yes'):
    probe_out = os.path.join("coot-molprobity", "molprobity-tmp-probe-dots.out")
    within_str = "(within " + \
                 str(radius) + \
                 " of " + \
                 str(point[0]) + \
                 ", " + \
                 str(point[1]) + \
                 ", " + \
                 str(point[2]) + \
                 ")"
    args = ["-mc", "-u", "-quiet", "-drop", "-stdbonds",
            "ALL",  # whole residues from sphere selection
                    # was needed to make this work
            # within_str  # problems with atom selection
            pdb_file_name]
    print "popen_comand on", probe_command, args
    popen_command(probe_command, args, [], probe_out, False)
    handle_read_draw_probe_dots_unformatted(probe_out, 0, 0)
    graphics_draw()

# Update the generic objects probe dots from residues within radius
# of the screen centre.
# 
# Return nothing interesting.
#
def probe_local_sphere(imol, radius):

  # We need to select more atoms than we probe because if the atom
  # selection radius and the probe radius are the same, then
  # sometimes the middle atom of a bonded angle set is missing
  # (outside the sphere) and that leads to bad clashes.
  # There are also bad clashed at the edge when alt-confed atoms
  # are not selected but non-alt-confs are which leads to missing atoms
  # in a bond angle and therefore clashes.

  pt = rotation_centre()
  imol_new = new_molecule_by_sphere_selection(imol, pt[0], pt[1], pt[2],
                                              radius, 0)
  set_mol_displayed(imol_new, 0)
  set_mol_active   (imol_new, 0)
  pdb_name = "molprobity-tmp-reference-file.pdb"
  make_directory_maybe("coot-molprobity")
  write_pdb_file_for_molprobity(imol_new, pdb_name)

  get_probe_dots_from(pdb_name, pt, radius)
  close_molecule(imol_new)


def probe_local_sphere_active_atom(radius=5.0):
  
  with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                             aa_ins_code, aa_atom_name, aa_alt_conf]:
    probe_local_sphere(aa_imol, radius)

  
# add in the conn files by concatting.
#
def write_pdb_file_for_molprobity(imol, pdb_name):

  tmp_pdb_name = file_name_sans_extension(pdb_name) + "-tmp.pdb"

  write_pdb_file(imol, tmp_pdb_name)

  # Let's add on the connectivity cards of the residues that
  # molprobity doesn't know about (which I presume are all
  # non-standard residues).  Cut out (filter) files that didn't
  # write properly.
  #
  conn_file_names = []
  for res_name in non_standard_residue_names(imol):
    f_name = os.path.join("coot-molprobity", "conn-" + res_name + ".txt")
    status = write_connectivity(res_name, f_name)
    if (status == 1):
      conn_file_names.append(f_name)

  # now, add (append) each of the con-file-names to the end of
  # pdb-name
  if (not os.path.isfile(tmp_pdb_name)):
    print "ERROR:: tmp file name not found", tmp_pdb_name
  else:
    conn_file_names.insert(0, tmp_pdb_name)
    fn_all = open(pdb_name, 'w')
    for conn_file in conn_file_names:
      if (os.path.isfile(conn_file)):
        cf = open(conn_file, 'r')
        data = cf.read()
        cf.close()
        fn_all.write(data)
    fn_all.close()

# not sure if the following ones should reside here?! Maybe rather in coot_gui?

# a toggle function for the main toolbar to switch on probe dots post refine
# and for chis/rotamers
#
def toggle_interactive_probe_dots(widget=None):
  """a toggle function for the main toolbar to switch on probe dots post refine
  and for chis/rotamers

  Keyword arguments:
  widget -- can be passed from the toolbutton
  
  """

  if widget:
    if widget.get_active():
      # the button is toggled on
      set_do_probe_dots_on_rotamers_and_chis(1)
      set_do_probe_dots_post_refine(1)
    else:
      set_do_probe_dots_on_rotamers_and_chis(0)
      set_do_probe_dots_post_refine(0)
  else:
    # no alternative for now (could just go by state and change back and forth)
    print "BL WARNING:: no widget"


def set_reduce_het_dict():
  """helper function to set env variable REDUCE_HET_DICT
  returns True if could be set (or already set), False if it couldnt be set.
  
  """

  global reduce_command

  # BL says: I think we should set REDUCE_HET_DICT
  # so let's set REDUCE_HET_DICT if not set already!
  # not sure if needed any more since we write a connectivity
  # file - let's see FIXME!!
  # if we have probe/reduce from ccp4 we should pick it up there!?
  # should be in $CCP4/share/reduce/

  red_het_dict = os.getenv('REDUCE_HET_DICT')
  red_het_file = "reduce_wwPDB_het_dict.txt"
  
  #negative path
#  if (not red_het_dict or
#      (red_het_dict and not os.path.isfile(red_het_dict))):

  # positive path (maybe should check that red_het_dict is a string?!)
  if (red_het_dict and
      (os.path.isfile(red_het_dict))):
    return True
  else:
    # not set and/or not found

    # first check if from ccp4 then it is in $CCP4/share/reduce/
    ccp4_dir = os.getenv('CCP4')  
    if (ccp4_dir and (os.path.normpath(ccp4_dir) in os.path.normpath(reduce_command))):
      # have ccp4 reduce
      red_het_dict = os.path.join(ccp4_dir, "share", "reduce",
                                 red_het_file)
      if os.path.isfile(red_het_dict):
        os.environ['REDUCE_HET_DICT'] = red_het_dict
        return True
      else:
        print "BL WARNING:: could neither find nor set REDUCE_HET_DICT !"
        return False
      
    else:
      # not from ccp4
      
      # we assume the dic is in same dir as reduce command
      full_reduce_command = find_exe(reduce_command, "PATH")
      dir, tmp = os.path.split(full_reduce_command)
      red_het_dict = os.path.join(dir, red_het_file)
      if (os.path.isfile(red_het_dict)):
        os.environ['REDUCE_HET_DICT'] = red_het_dict
        return True
      else:
        # before we give up check in share/coot
        # this is where the windows installer shall put it
        prefix_dir = os.getenv("COOT_PREFIX")
        if not prefix_dir:
          pkg_data_dir = pkgdatadir()
        else:
          pkg_data_dir = os.path.join(prefix_dir, "share", "coot")
        red_het_file = os.path.join(pkg_data_dir, red_het_file)
        if (os.path.isfile(red_het_dict)):
          os.environ['REDUCE_HET_DICT'] = red_het_dict
          return True
        else:
          # finally give up
          print "BL WARNING:: could neither find nor set REDUCE_HET_DICT !"
          return False
      
