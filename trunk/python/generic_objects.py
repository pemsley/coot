# generic-objects.py
#
# Copyright 2006 by Paul Emsley, The University of York
# Copyright 2006 by Bernhard Lohkamp
 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA

global probe_command
global reduce_command
# BL says:: put in definitions of probe and reduce command here:
# actually not used any more, I think!
probe_command_name = "probe"
reduce_command_name = "reduceV2.13.2.win"

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
	  
  if (pygtk_flag):
	  # Now we run the gui
	def delete_event(*args):
		# BL says: first we shall close the generic objects
		for generic_object_number in range(n_objects):
			set_display_generic_object(generic_object_number,0)
		gen_window.destroy()
		return False


	def check_button_callback(widget, generic_object_number):
		button_state = widget.get_active()
		object_state = generic_object_is_displayed_qm(generic_object_number)
		if ((button_state == True) and (object_state == 0)):
			set_display_generic_object(generic_object_number,1)
		if ((button_state == False) and (object_state == 1)):
			set_display_generic_object(generic_object_number,0)

	n_objects = number_of_generic_objects()
	
	if (n_objects > 0):
		gen_window = gtk.Window(gtk.WINDOW_TOPLEVEL)
		gen_window.set_title("Generic objects")
		vbox = gtk.VBox(False, 0)

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

					vbox.add(frame)
					frame.add(check_button)
					frame.show()
					check_button.show()

		gen_window.connect("delete_event", delete_event)

		vbox.show()
		gen_window.add(vbox)
		gen_window.set_border_width(10)
		gen_window.show()

       
# run molprobity (well reduce and probe) to make generic objects (and
# display the generic objects gui)
#
def probe(imol):
    import os
    global reduce_command, probe_command
    
    if is_valid_model_molecule(imol):
       #reduce seems to have a different name in win than otherwise
       if os.name == 'nt':
		if not (os.path.isfile(reduce_command)):
			reduce_command = find_exe(reduce_command_name, "PATH")
       else:
		# this is assuming reduce is called 'reduce' on other systems (which it probably isnt)
		if not (os.path.isfile(reduce_command)):
			reduce_command = find_exe("reduce", "PATH")
       # we need to check if probe_command is defined too
       if not(os.path.isfile(probe_command)):
		probe_command = find_exe(probe_command_name, "PATH")
       make_directory_maybe("coot-molprobity")
       mol_pdb_file = "coot-molprobity/for-reduce.pdb"
       reduce_out_pdb_file = "coot-molprobity/reduced.pdb"
       write_pdb_file(imol,mol_pdb_file)
       if (reduce_command):

          # BL says: I think we should set REDUCE_HET_DICT
          # so let's set REDUCE_HET_DICT if not set already!
          if (not os.getenv('REDUCE_HET_DICT')):
             # we assume the dic is in same dir as reduce command
             dir, tmp = os.path.split(reduce_command)
             connection_file = os.path.join(dir,'reduce_het_dict.txt')
             if (os.path.isfile(connection_file)):
                os.putenv('REDUCE_HET_DICT',connection_file)
             else:
                print "BL WARNING:: could neither find nor set REDUCE_HET_DICT !"
          else:
             # REDUCE_HET_DICT is defined
             pass

          print "BL INFO:: run reduce as ", reduce_command + " " + mol_pdb_file + " > " + reduce_out_pdb_file
          status = os.popen(reduce_command + " -build -oldpdb " + mol_pdb_file + " > " + reduce_out_pdb_file,'r')
#          status = os.popen(reduce_command + " " + mol_pdb_file + " > " + reduce_out_pdb_file,'r')
          reduce_status = status.close()
          if (probe_command):
             probe_pdb_in = "coot-molprobity/probe-in.pdb"
             probe_out = "coot-molprobity/probe-out.pdb"
             prepare_file_for_probe(reduce_out_pdb_file,probe_pdb_in)
             probe_command_line = probe_command + " -u -stdbonds -mc \"ALL ALL\" " + probe_pdb_in + " > " + probe_out
             print "BL INFO:: run probe as ", probe_command_line
             status = os.popen(probe_command_line,'r')
             probe_status = status.close()
                         # by default, we don't want to click on the
                         # imol-probe molecule (I think :-)
             recentre_status = recentre_on_read_pdb()
             novalue = set_recentre_on_read_pdb(0)
             imol_probe = read_pdb(probe_pdb_in)
             if recentre_status == 1: set_recentre_on_read_pdb(1)
             # toggle_active_mol(imol_probe) let's not do
             # that actually.  I no longer think that the
             # new probe molecule should not be clickable
             # when it is initally displayed (that plus
             # there is some active/displayed logic problem
             # for the molecules, which means that after
             # several probes, the wrong molecule is
             # getting refined).

	     print "BL DEBUG:: handle_dots with ", probe_out, imol_probe
             handle_read_draw_probe_dots_unformatted(probe_out,imol_probe,1)
             generic_objects_gui()
             graphics_draw()
          else:
            # couldnt find probe
            print "BL WARNING:: Could not locate the program probe!! Please check if installed!"
       else:
         # couldnt find reduce
         print "BL WARNING:: Could not locate the program reduce!! Please check if installed!"

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

    probe_pdb_in_1 = "molprobity-tmp-reference-file.pdb"
    probe_pdb_in_2 = "molprobity-tmp-moving-file.pdb"
    probe_out = "coot-molprobity/molprobity-tmp-probe-dots.out"
    atom_sel = "\"(file1 within " + str(radius) + " of " + str(x_cen) + ", " + str(y_cen) + ", " + str(z_cen) + ", " + "not water not (chain" + string.lower(chain_id) + " " + str(res_no) + ")),file2\""

    if not (os.path.isfile(probe_command)):
	probe_command = find_exe("probe", "PATH")
    if (probe_command):
       probe_line = " -mc -u -quiet -drop -stdbonds -both " + atom_sel + " file2 " + probe_pdb_in_1 + " " + probe_pdb_in_2
#       print "BL INFO:: run probe as ", probe_command + probe_line
       status = os.popen(probe_command + probe_line + " > " + probe_out,'r')
       status.close()

       # don't show the gui, so the imol is not needed/dummy.
       handle_read_draw_probe_dots_unformatted(probe_out,0,0)
       graphics_draw()


#probe(0)

