# coot-gui.py
#
# Copyright 2001 by Neil Jerram
# Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
# Copyright 2007 by Paul Emsley
# And Add: "Copyright 2007 by Bernhard Lohkamp"? Hmm...
# with adaptions from
#   Interactive Python-GTK Console
#   Copyright (C), 1998 James Henstridge <james@daa.com.au>

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


# BL adaption of coot-gui scripting window running with pygtk
global do_function

# Fire up the coot scripting gui.  This function is called from the
# main C++ code of coot.  Not much use if you don't have a gui to
# type functions in to start with.
#

import pygtk, gtk, pango
import gobject

import types
from types import *

global histpos
global history
histpos = 0
history = ['']

# Fire up the coot scripting gui.  This function is called from the
# main C++ code of coot.  Not much use if you don't have a gui to
# type functions in to start with.
#
def coot_gui():

   import sys, string
   import re
   def delete_event(*args):
       window.destroy()
       return False

#   def destroy(*args):
#       gtk.main_quit()

   def hist_func(widget, event, entry, textbuffer, text):
       if event.keyval in (gtk.keysyms.KP_Up, gtk.keysyms.Up):
          entry.emit_stop_by_name("key_press_event")
          historyUp()
          gobject.idle_add(entry.grab_focus)
       elif event.keyval in (gtk.keysyms.KP_Down, gtk.keysyms.Down):
          entry.emit_stop_by_name("key_press_event")
          historyDown()
          gobject.idle_add(entry.grab_focus)
       else:
          pass

   def historyUp():
       global histpos
       if histpos > 0:
          l = entry.get_text()
          if len(l) > 0 and l[0] == '\n': l = l[1:]
          if len(l) > 0 and l[-1] == '\n': l = l[:-1]
          history[histpos] = l
          histpos = histpos - 1
          entry.set_text(history[histpos])
#       print "BL DEBUG:: hist and pos", history, histpos
                        
   def historyDown():
       global histpos
       if histpos < len(history) - 1:
          l = entry.get_text()
          if len(l) > 0 and l[0] == '\n': l = l[1:]
          if len(l) > 0 and l[-1] == '\n': l = l[:-1]
          history[histpos] = l
          histpos = histpos + 1
          entry.set_text(history[histpos])
#       print "BL DEBUG:: hist and pos", history, histpos

   def warning_event(entry, event):
      if (event.state & gtk.gdk.CONTROL_MASK):
        # BL says:: for windows we want to exclude Ctrl so that we can copy things
        # alhough not really needed since we can use middle mouse  or drag to copy
        return
      else:
       tag1 = textbuffer.create_tag(foreground="red",weight=pango.WEIGHT_BOLD)
       insert_tag_text(tag1,"Don't type here.  Type in the white entry bar\n")
       insert_prompt()
       while gtk.events_pending():
          gtk.main_iteration(False)
       return

   def scroll_to_end():
       end = textbuffer.get_end_iter()
       textbuffer.place_cursor(end)
       text.scroll_mark_onscreen(textbuffer.get_insert())
       return gtk.FALSE
       
   def insert_normal_text(s):
       end = textbuffer.get_end_iter()
       textbuffer.insert_with_tags(end,str(s),textbuffer.create_tag(foreground="black"))
       scroll_to_end()

   def insert_tag_text(tagg,s):
       end = textbuffer.get_end_iter()
       textbuffer.insert_with_tags(end,str(s),tagg)
       scroll_to_end()

   def insert_prompt():
       end = textbuffer.get_end_iter()
       insert_normal_text("coot >> ")
       scroll_to_end()

   def do_function(widget, entry):
       global histpos
       entry_text = entry.get_text()
       print "BL INFO:: command input is: ", entry_text
       insert_tag_text(textbuffer.create_tag(foreground="red"),entry_text + "\n")
       while gtk.events_pending():
         gtk.main_iteration(False)
       his = False
       res = None
       try:
             res = eval(entry_text)
             his = True
       except SyntaxError:
         try:
             exec entry_text in globals()
             res = None
             his = True
         except:
          guile_function =  test_and_translate_guile(entry_text)
          if guile_function:
             print "BL INFO::  We have a guile command!"
             insert_normal_text("BL INFO:: Detected guile scripting!\nYou should use python commands!!\nBut I'm a nice guy and translated it for you, this time...!\n")
          else:
             insert_normal_text("BL WARNING:: Python syntax error!\n(Or you attempted to use an invalid guile command...)\n")
             type_error, error_value = sys.exc_info()[:2]
             error = str(error_value)
             insert_normal_text("Python error:\n") 
             insert_normal_text(error + "\n")

       except:
          guile_function =  test_and_translate_guile(entry_text)
          if guile_function:
             insert_normal_text("BL INFO:: Detected guile scripting!\nYou should use python commands!!\nBut I'm a nice guy and translated it for you, this time...!\n")
          else:
             insert_normal_text("BL WARNING:: Python error!\n(Or you attempted to use an invalid guile command...)\n")
             type_error, error_value = sys.exc_info()[:2]
             error = str(error_value)
             insert_normal_text("Python error:\n") 
             insert_normal_text(error + "\n")

       if res is not None:
             print "BL INFO:: result is", res
             insert_normal_text(str(res) + "\n")

       if his:
             l = entry_text + '\n'
             if len(l) > 1 and l[0] == '\n': l = l[1:]
             histpos = len(history) - 1
             # print "pedebug", l[-1]
             # can't test l[-] because python gives us a error
             # SystemError: ../Objects/longobject.c:223: bad argument to internal function
             # when we do so (on set_monomer_restriaints()).
             # if len(l) > 0 and l[-1] == '\n':
             if len(l) > 0:
                history[histpos] = l[:-1]
             else:
                history[histpos] = l

             histpos += 1
             history.append('')
       
       entry.set_text("")
       entry_text = ""
       insert_prompt()
       
   # BL says:: now we check if entry was guile command (by mistake?!)
   # if so, let's translate it
   #
   def test_and_translate_guile(py_func):
       #test for - or ( at start
       if ((string.find(py_func,"-")>0) or (string.find(py_func,"(")==0)):
          #remove ()
          tmp = string.rstrip(string.lstrip(py_func,("(")),")")
          tmp_list = string.split(tmp)
          tmp_command = string.replace(tmp_list[0],"-","_")
          tmp_list[0:1] = []
          #convert string args from list to int, float, string...
          arg = tmp_list
          for i in range(len(arg)):
              try:
                 arg[i] = int(arg[i])
              except ValueError:
                 try:
                    arg[i] = float(arg[i])
                 except ValueError:
                    # it should be a string then
                    try:
                       arg[i] = arg[i].strip('"')
                    except:
                       print "BL WARNING:: unknown type of argument!"

          python_function = tmp_command + "(*arg)"
#          print "BL DEBUG:: python func and args are", python_function, tmp_command, arg

          try:
             res = eval(python_function)
          except SyntaxError:
             exec python_function in locals()
             res = None
             return True
          except:
             return False
          if res is not None:
             print "BL DEBUG:: result is", res
             insert_normal_text(str(res) + "\n")
             return True
       else:
          print "This is not a guile command!"
          return False


   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   window.set_title("Coot Python Scripting")
   text = gtk.TextView()
   textbuffer = text.get_buffer()
   scrolled_win = gtk.ScrolledWindow()
   entry = gtk.Entry()
   completion = gtk.EntryCompletion()
   entry.set_completion(completion)
   liststore = gtk.ListStore(gobject.TYPE_STRING)
   for i in globals():
       liststore.append([i])
   completion.set_model(liststore)
   completion.set_text_column(0)
#   completion.connect("match-selected", match_cb)
   close_button = gtk.Button("  Close  ")
   vbox = gtk.VBox(False, 0)
   hbox = gtk.HBox(False, 0)
   label = gtk.Label("Command: ")
   ifont = gtk.gdk.Font("fixed")
   window.set_default_size(400,250)
   window.add(vbox)
   vbox.set_border_width(5)

   hbox.pack_start(label, False, False, 0)   # expand fill padding ???
   hbox.pack_start(entry, True, True, 0)
   vbox.pack_start(hbox, False, False, 5)
   scrolled_win.add(text)
   vbox.add(scrolled_win)
   scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
   vbox.pack_end(close_button, False, False, 5)

   entry.set_text("")
   entry.set_editable(True)
#   entry.connect("activate", do_function, entry, textbuffer, text)
   entry.connect("activate", do_function, entry)
#bl testing
   entry.connect("key_press_event", hist_func, entry, textbuffer, text)
   window.connect("delete_event", delete_event)
   close_button.connect("clicked", delete_event)

   text.set_editable(False)
   text.set_cursor_visible(False)
   text.connect("key_press_event", warning_event)
   end = textbuffer.get_end_iter()
   tag = textbuffer.create_tag(foreground="black")
   textbuffer.insert_with_tags(end,"coot >> ",tag)
#   textbuffer.set_text("coot >> ")

   text.modify_base(gtk.STATE_NORMAL, gtk.gdk.color_parse("#bfe6bf"))
   window.show_all()
#   gtk.main()


# The callback from pressing the Go button in the smiles widget, an
# interface to run libcheck.
#
def handle_smiles_go(tlc_entry,smiles_entry):

    import os, stat

    tlc_text = tlc_entry.get_text()
    smiles_text = smiles_entry.get_text()

    if (len(smiles_text) > 0):
       if ((len(tlc_text) > 0) and (len(tlc_text) < 4)): three_letter_code = tlc_text
       # next should actually not happen since I restrict the input to 3 letters!
       elif (len(tlc_text) > 0): three_letter_code = tlc_text[0:3]
       else: three_letter_code = "DUM"
          
       # ok let's run libcheck
       smiles_file = "coot-" + three_letter_code + ".smi"
       libcheck_data_lines = ["N","MON " + three_letter_code, "FILE_SMILE " + smiles_file,""]
       log_file_name = "libcheck-" + three_letter_code + ".log"
       pdb_file_name = "libcheck_" + three_letter_code + ".pdb"
       cif_file_name = "libcheck_" + three_letter_code + ".cif"
#       print "BL DEBUG:: all others: ", libcheck_data_lines, log_file_name, pdb_file_name ,cif_file_name
       #write smiles strings to a file
       smiles_input = file(smiles_file,'w')
       smiles_input.write(smiles_text)
       smiles_input.close()
       #now let's run libcheck (based on libcheck.py)
       libcheck_exe = find_exe("libcheck", "CCP4_BIN", "PATH")

       if (libcheck_exe):
          # BL says: as with refmac we have to write a file with paramters to run libcheck
          libcheck_input_file = "coot-libcheck-input.txt"
          input = file(libcheck_input_file,'w')
          for data in libcheck_data_lines:
             input.write(data + '\n')
          input.close()
          status = os.popen(libcheck_exe + ' < ' + libcheck_input_file + ' > ' + log_file_name,'r')
          libcheck_status = status.close()
          os.remove(libcheck_input_file)
          log_file_size = os.stat(log_file_name)[stat.ST_SIZE]
          if ((not libcheck_status) and (os.path.isfile(log_file_name)) and (log_file_size > 0)):
             check_libcheck_logfile(log_file_name)
             # means we have a log file, dunno how else to check for status currently!?
             if (os.path.isfile("libcheck.lib")):
                if (os.path.isfile(cif_file_name)):
                   # if we have cif_file_name already we cant move to it
                   # and I dont want to overwrite it, so we make a backup
                   try:
                      os.rename(cif_file_name,cif_file_name + ".bak")
                   except OSError:
                      # bak file exists, so let's remove and overwrite it
                      os.remove(cif_file_name + ".bak")
                      os.rename(cif_file_name,cif_file_name + ".bak")
                      print "BL INFO:: overwriting %s since same three letter code used again..." %(cif_file_name + ".bak")
                   else:
                      print "BL WARNING:: %s exists and we cant remove/overwrite it! So libcheck wont run!" %(cif_file_name + ".bak")
                os.rename("libcheck.lib",cif_file_name)
                print "BL INFO:: renamed %s to %s " %("libcheck.lib",cif_file_name)
             sc = rotation_centre()
             print "BL INFO:: reading ",pdb_file_name
             imol = handle_read_draw_molecule_with_recentre(pdb_file_name,0)
             if (is_valid_model_molecule(imol)):
                mc = molecule_centre(imol)
                sc_mc = [sc[i]-mc[i] for i in range(len(mc))]
                translate_molecule_by(imol,*sc_mc)
             read_cif_dictionary(cif_file_name)
          else: print "BL WARNING:: no log file or no length!"
       else: print " BL WARNING:: libcheck not found!"
    else:
       print "BL WARNING:: Wrong input (no smiles text)! Can't continue!"

# smiles GUI
def smiles_gui():

    def delete_event(*args):
        smiles_window.destroy()
        return False
   
    def go_button_pressed(widget, tlc_entry, smiles_entry,smiles_window):
        handle_smiles_go(tlc_entry,smiles_entry)
        smiles_window.destroy()
        delete_event()

    def smiles_connect(widget, event, tlc_entry, smiles_entry, smiles_window):
        if (event.keyval == 65293):
           handle_smiles_go(tlc_entry,smiles_entry)
           smiles_window.destroy()

    smiles_window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    smiles_window.set_title("SMILES GUI")
    vbox = gtk.VBox(False, 0)
    hbox1 = gtk.HBox(False, 0)
    hbox2 = gtk.HBox(False, 0)
    tlc_label = gtk.Label("  3-letter code ")
    tlc_entry = gtk.Entry(max=3)
    tlc_entry.set_text("")
    smiles_label = gtk.Label("SMILES string ")
    smiles_entry = gtk.Entry()
    go_button = gtk.Button("  Go  ")
    vbox.pack_start(hbox1, False, False, 0)
    vbox.pack_start(hbox2, False, False, 0)
    vbox.pack_start(go_button, False, False, 6)
    hbox1.pack_start(tlc_label, False, False, 0)
    hbox1.pack_start(tlc_entry, False, False, 0)
    hbox2.pack_start(smiles_label, False, False, 0)
    hbox2.pack_start(smiles_entry, True, True, 0)
    smiles_window.add(vbox)
    vbox.set_border_width(6)
    
    smiles_entry.connect("key-press-event", smiles_connect, tlc_entry, smiles_entry, smiles_window)
    go_button.connect("clicked", go_button_pressed,tlc_entry,smiles_entry,smiles_window)
    smiles_window.connect("delete_event",delete_event)
    
    smiles_window.show_all()

# Generic single entry widget
# 
# Pass the hint labels of the entries and a function that gets called
# when user hits "Go".  The handle-go-function accepts one argument
# that is the entry text when the go button is pressed.
#
def generic_single_entry(function_label, entry_1_default_text, go_button_label, handle_go_function):

    def delete_event(*args):
       window.destroy()
       return False

    def go_function_event(widget, smiles_entry, *args):
        handle_go_function(smiles_entry.get_text())
        delete_event()

    def key_press_event(widget, event, smiles_entry, *args):
        if (event.keyval == 65293):
           handle_go_function(smiles_entry.get_text())
           delete_event()

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    vbox = gtk.VBox(False, 0)
    hbox1 = gtk.HBox(False, 0)
    hbox2 = gtk.HBox(False, 0)
    hbox3 = gtk.HBox(False, 0)
    function_label = gtk.Label(function_label)
    smiles_entry = gtk.Entry()
    cancel_button = gtk.Button("  Cancel  ")
    go_button = gtk.Button(go_button_label)

    vbox.pack_start(hbox1, False, False, 0)
    vbox.pack_start(hbox2, False, False, 4)
    vbox.pack_start(hbox3, False, False, 4)
    hbox3.pack_start(go_button, True, True, 6)
    hbox3.pack_start(cancel_button, True, True, 6)
    hbox1.pack_start(function_label, False, 0)
    hbox2.pack_start(smiles_entry, False, 0)
    window.add(vbox)
    vbox.set_border_width(6)
 
    if isinstance(entry_1_default_text,types.StringTypes):
       smiles_entry.set_text(entry_1_default_text)
    else:
       print "BL WARNING:: entry_1_default_text was no string!!"
   
    cancel_button.connect("clicked", delete_event)

    go_button.connect("clicked", go_function_event, smiles_entry)

    smiles_entry.connect("key_press_event", key_press_event, smiles_entry)

    window.show_all()

# generic double entry widget, now with a check button
# 
# pass a the hint labels of the entries and a function 
# (handle-go-function) that gets called when user hits "Go" (which
# takes two string aguments and the active-state of the check button
# (either True of False).
# 
# if check-button-label not a string, then we don't display (or
# create, even) the check-button.  If it *is* a string, create a
# check button and add the callback handle-check-button-function
# which takes as an argument the active-state of the the checkbutton.
#
def generic_double_entry(label_1, label_2, entry_1_default_text, entry_2_default_text, check_button_label, handle_check_button_function, go_button_label, handle_go_function):

    def delete_event(*args):
       window.destroy()
       return False

    def go_function_event(*args):
	if check_button:
	        handle_go_function(tlc_entry.get_text(), smiles_entry.get_text(), check_button.get_active())
	else:
        	handle_go_function(tlc_entry.get_text(), smiles_entry.get_text())
        delete_event()

    def key_press_event(widget, event, tlc_entry, smiles_entry, check_button):
        if (event.keyval == 65293):
		if check_button:
			handle_go_function(tlc_entry.get_text(), smiles_entry.get_text(), check_button.get_active())
	        else:
			handle_go_function(tlc_entry.get_text(), smiles_entry.get_text())
           	delete_event()

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    vbox = gtk.VBox(False, 0)
    hbox1 = gtk.HBox(False, 0)
    hbox2 = gtk.HBox(False, 0)
    hbox3 = gtk.HBox(False, 0)
    tlc_label = gtk.Label(label_1)
    tlc_entry = gtk.Entry()
    smiles_label = gtk.Label(label_2)
    smiles_entry = gtk.Entry()
    h_sep = gtk.HSeparator()
    cancel_button = gtk.Button("  Cancel  ")
    go_button = gtk.Button(go_button_label)

    vbox.pack_start(hbox1, False, False, 0)
    vbox.pack_start(hbox2, False, False, 0)
    hbox3.pack_start(go_button, True, False, 6)
    hbox3.pack_start(cancel_button, True, False, 6)
    hbox1.pack_start(tlc_label, False, 0)
    hbox1.pack_start(tlc_entry, False, 0)
    hbox2.pack_start(smiles_label, False, 0)
    hbox2.pack_start(smiles_entry, False, 0)

    if type(check_button_label) is StringType:

	def check_callback(*args):
		active_state = c_button.get_active()
		handle_check_button_function(active_state)

	c_button = gtk.CheckButton(check_button_label)
	vbox.pack_start(c_button, False, False, 2)
	c_button.connect("toggled", check_callback)
	check_button = c_button
    else:
	check_button = False 	# the check-button when we don't want to see it

    vbox.pack_start(h_sep, True, False, 3)
    vbox.pack_start(hbox3, False, False, 0)
    window.add(vbox)
    vbox.set_border_width(6)

    if isinstance(entry_1_default_text,types.StringTypes):
       tlc_entry.set_text(entry_1_default_text)
    else:
       print "BL WARNING:: entry_1_default_text was no string!!"
 
    if isinstance(entry_2_default_text,types.StringTypes):
       smiles_entry.set_text(entry_2_default_text)
    else:
       print "BL WARNING:: entry_2_default_text was no string!!"

    go_button.connect("clicked", go_function_event)

    smiles_entry.connect("key-press-event", key_press_event, tlc_entry, smiles_entry, check_button)

    window.show_all()

# generic double entry widget, now with a check button
#
# OLD
#
# 
# pass a the hint labels of the entries and a function 
# (handle-go-function) that gets called when user hits "Go" (which
# takes two string aguments and the active-state of the check button
# (either True of False).
# 
# if check-button-label not a string, then we don't display (or
# create, even) the check-button.  If it *is* a string, create a
# check button and add the callback handle-check-button-function
# which takes as an argument the active-state of the the checkbutton.
#
def generic_multiple_entries_with_check_button(entry_info_list, check_button_info, go_button_label, handle_go_function):

    def delete_event(*args):
       window.destroy()
       return False

    def go_function_event(*args):
       print "Here.................. check-button is ", check_button
       if check_button:
	        handle_go_function(map(entry.get_text(), entries), check_button.get_active())
       else:
        	handle_go_function(map(entry.get_text(), entries))
       delete_event()

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    vbox = gtk.VBox(False, 0)
    hbox3 = gtk.HBox(False, 0)
    h_sep = gtk.HSeparator()
    cancel_button = gtk.Button("  Cancel  ")
    go_button = gtk.Button(go_button_label)

    # all the labelled entries
    #
    entries = []
    for entry_info in entry_info_list:

       entry_hint_text = entry_info[0]
       entry_1_default_text = entry_info[1]
       hbox1 = gtk.HBox(False, 0)

       label = gtk.Label(entry_1_hint_text)
       entry = gtk.Entry()
       entries.append(entry)

       if type(entry_1_default_text) is StringType:
          entry.set_text(entry_1_default_text)

       hbox1.pack_start(label, False, 0)
       hbox1.pack_start(entry, False, 0)
       vbox.pack_start(hbox1, False, False, 0)

    print "debug:: check-button-info: ", check_button_info
    if not (type(entry_info_list) is ListType and len(check_button_info) == 2):
       print "check_button_info failed list and length test"
       check_button = False
    else:
       if type(check_button_info[0]) is StringType:

          def check_callback(*args):
		active_state = c_button.get_active()
		handle_check_button_function(active_state)
                
          c_button = gtk.CheckButton(check_button_info[0])
          vbox.pack_start(c_button, False, False, 2)
          c_button.connect("toggled", check_callback)
          check_button = c_button
       check_button = False      # the check-button when we don't want to see it

    print "Here check button creation.................. check-button is ", check_button
    vbox.pack_start(h_sep, True, False, 3)
    vbox.pack_start(hbox3, False, False, 0)
    window.add(vbox)
    vbox.set_border_width(6)

    hbox3.pack_start(go_button, True, False, 6)
    hbox3.pack_start(cancel_button, True, False, 6)


    cancel_button.connect("clicked", delete_event)
    go_button.connect("clicked", go_function_event)

    window.show_all()

# A demo gui to move about to molecules
def molecule_centres_gui():

    def delete_event(*args):
        window.destroy()
        return False

    def callback_func(widget,molecule_number,label):
        s = "Centred on " + label
        add_status_bar_text(s)
        set_rotation_centre(*molecule_centre(molecule_number))

    # first, we create a window and a frame to be put into it.
    # 
    # we create a vbox (a vertical box container) that will contain the
    # buttons for each of the coordinate molecules
    # 
    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    frame = gtk.Frame("Molecule Centres")
    vbox = gtk.VBox(False, 3)
    
    # add the frame to the window and the vbox to the frame
    # 
    window.add(frame)
    frame.add(vbox)
    vbox.set_border_width(6)

    # for each molecule, test if this is a molecule that has a
    # molecule (it might have been closed or contain a map).  The we
    # construct a button label that is the molecule number and the
    # molecule name and create a button with that label.
    # 
    # then we attach a signal to the "clicked" event on the newly
    # created button.  In the function that subsequently happen (on a
    # click) we add a text message to the status bar and set the
    # graphics rotation centre to the centre of the molecule.  Each
    # of these buttons is packed into the vbox (the #f #f means no
    # expansion and no filling).  The padding round each button is 3
    # pixels.
    # 
    for molecule_number in molecule_number_list():
        if (is_valid_model_molecule(molecule_number)):
           name = molecule_name(molecule_number)
           label = str(molecule_number) + " " + name
           button = gtk.Button(label)
           button.connect("clicked",callback_func,molecule_number,label)
           vbox.add(button)
    window.show_all()

# A BL function to analyse the libcheck log file
#
def check_libcheck_logfile(logfile_name):
    # to check if error or warning occured in libcheck
    # maybe want to have something like that for all logfiles
    # and produce a warning for users...

    fp = open(logfile_name, 'r')
    lines = fp.readlines()
    for line in lines:
        if "WARNING: no such keyword :FILE" in line:
           print "BL WARNING:: libcheck didn't seem to run ok! Please check output carefully!!"
        else: pass
    fp.close()

# old coot test
# BL says:: probably means people will just hassle me ...
# so not implemented yet
def old_coot_qm():
# need to agree on a time format too..
#    new_release_time = 1200000000
#    current_time = time.localtime()
#    time_diff = current_time - new_release_time
    return 0

# We can either go to a place (in which case the element is a list of
# button label (string) and 3 numbers that represent the x y z
# coordinates) or an atom (in which case the element is a list of a
# button label (string) followed by the molecule-number chain-id
# residue-number ins-code atom-name altconf)
# 
# e.g. interesting_things_gui(
#       "Bad things by Analysis X",
#       [["Bad Chiral",0,"A",23,"","CA","A"],
#        ["Bad Density Fit",0,"B",65,"","CA",""],
#        ["Interesting blob",45.6,46.7,87.5]]
# 
def interesting_things_gui(title,baddie_list):

    interesting_things_with_fix_maybe(title, baddie_list)

# In this case, each baddie can have a function at the end which is
# called when the fix button is clicked.
# And extra string arguments can give the name of the button and a tooltip
# as these two are ambigious, they can only be:
# 1 extra string:  button name
# 2 extra strings: button name and tooltip
# 
def interesting_things_with_fix_maybe(title, baddie_list):

   # does this baddie have a fix at the end?.  If yes, return the
   # func, if not, return #f.
   # BL says::
   # OBS in python the fix function is a list with the 1st element the function
   # hope this will work!?
   def baddie_had_fix_qm(baddie):
       
       l = len(baddie)
       if (l < 5):
          # certainly no fix
          return False, False, False
       else:
          func_maybe1 = baddie[l-1]
          func_maybe2 = baddie[l-2]
          func_maybe3 = baddie[l-3]

          if (isinstance(func_maybe1, types.ListType) and len(func_maybe1)>0):
             # the last one is probably a funcn (no button name)
             func_maybe_strip = func_maybe1[0]
#             print "BL DEBUG:: func_maybe_strip is", func_maybe_strip
             if (callable(func_maybe_strip)):
                return func_maybe1, False, False
             else:
                return False, False, False
             
          elif (isinstance(func_maybe2, types.ListType) and len(func_maybe2)>0
                and isinstance(func_maybe1, types.StringType)):
             # the second last is function, last is button name
             func_maybe_strip = func_maybe2[0]
             button_name = func_maybe1
             if (callable(func_maybe_strip)):
                return func_maybe2, button_name, False
             else:
                return False, False, False
             
          elif (isinstance(func_maybe3, types.ListType) and len(func_maybe3)>0
                and isinstance(func_maybe2, types.StringType)
                and isinstance(func_maybe1, types.StringType)):
             # the third last is function, second last is button name, last is tooltip
             func_maybe_strip = func_maybe3[0]
             button_name = func_maybe2
             tooltip_str = func_maybe1
             if (callable(func_maybe_strip)):
                return func_maybe3, button_name, tooltip_str
             elif (func_maybe_strip == "dummy"):
                return False, False, tooltip_str
             else:
                return False, False, False
          
             
          else:
             return False, False, False

   def fix_func_call(widget, call_func):
       func_maybe_strip = call_func[0]
       func_args = call_func[1:len(call_func)]
       func_maybe_strip(*func_args)

   def delete_event(*args):
       window.destroy()
       gtk.main_quit()
       return False

   def callback_func1(widget,coords):
       set_rotation_centre(*coords)

   def callback_func2(widget,mol_no,atom_info):
       print "Attempt to go to chain: %s resno: %s atom-name: %s" %(atom_info[0],atom_info[1],atom_info[2])
       set_go_to_atom_molecule(mol_no)
       success = set_go_to_atom_chain_residue_atom_name(*atom_info)
       if success == 0:           # failed?!
          new_name = unmangle_hydrogen_name(atom_info[2])
          success2 = set_go_to_atom_chain_residue_atom_name(atom_info[0],atom_info[1],new_name)
          if success2 == 0:
             print "Failed to centre on demangled name: ", new_name
             set_go_to_atom_chain_residue_atom_name(atom_info[0],atom_info[1]," CA ")

   # main body
   # to accomodated tooltips we need to either have a gtk.Window with gtk.main()
   # or a dialog and run() it! We use Window to not block events and hope not to
   # interfere with the gtk_main() of coot itself
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   scrolled_win = gtk.ScrolledWindow()
   outside_vbox = gtk.VBox(False,2)
   inside_vbox = gtk.VBox(False,0)

   window.set_default_size(250,250)
   window.set_title(title)
   inside_vbox.set_border_width(4)

   window.add(outside_vbox)
   outside_vbox.add(scrolled_win)
   scrolled_win.add_with_viewport(inside_vbox)
   scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)
   tooltips = gtk.Tooltips()

   for baddie_items in baddie_list:

       if baddie_items == '' :
          print 'Done'
       else:
          hbox = gtk.HBox(False,0) # hbox to contain Baddie button and
                                   # and fix it button

          label = baddie_items[0]
          button = gtk.Button(label)

          inside_vbox.pack_start(hbox,False,False,2)
          hbox.pack_start(button, True, True, 1)

          # add the a button for the fix func if it exists.  Add
          # the callback.
          fix_func, button_name, tooltip_str = baddie_had_fix_qm(baddie_items)
          if (fix_func):
             if (button_name):
                fix_button_name = button_name
             else:
                fix_button_name = "  Fix  "
             fix_button = gtk.Button(fix_button_name)
             hbox.pack_end(fix_button, False, False, 2)
             fix_button.show()
             fix_button.connect("clicked", fix_func_call, fix_func)

          if (tooltip_str):
             # we have a tooltip str
             tooltips.set_tip(button, tooltip_str)


          if (len(baddie_items) == 4):               # e.g. ["blob",1,2,3]
             # we go to a place
             coords = [baddie_items[1], baddie_items[2], baddie_items[3]]
             button.connect("clicked", callback_func1, coords)

          else:
             # we go to an atom
             mol_no = baddie_items[1]
             atom_info = [baddie_items[2], baddie_items[3], baddie_items[5]]
             button.connect("clicked", callback_func2, mol_no, atom_info)

       
   outside_vbox.set_border_width(4)
   ok_button = gtk.Button("  OK  ")
   outside_vbox.pack_start(ok_button, False, False, 6)
   ok_button.connect("clicked", delete_event)

   window.connect("destroy", delete_event)

   window.show_all()
   gtk.main()

#interesting_things_gui("Bad things by Analysis X",[["Bad Chiral",0,"A",23,"","CA","A"],["Bad Density Fit",0,"B",65,"","CA",""],["Interesting blob",45.6,46.7,87.5],["Interesting blob 2",45.6,41.7,80.5]])
#interesting_things_gui("Bad things by Analysis X",[["Bad Chiral",0,"A",23,"","CA","A",[print_sequence_chain,0,'A']],["Bad Density Fit",0,"B",65,"","CA",""],["Interesting blob",45.6,46.7,87.5],["Interesting blob 2",45.6,41.7,80.5]])


# Fill an option menu with the "right type" of molecules.  If
# filter_function returns True then add it.  Typical value of
# filter_function is valid-model-molecule_qm
#
def fill_option_menu_with_mol_options(menu, filter_function):

    mols = 0
    n_molecules = graphics_n_molecules()
   
    for mol_no_ls in molecule_number_list():
        if filter_function(mol_no_ls):
           label_str = molecule_name(mol_no_ls)
           if (isinstance(label_str,types.StringTypes)):
              mlabel_str = str(mol_no_ls) + " " + label_str
              menu.append_text(mlabel_str)
              menu.set_active(0)
              mols += 1
           else:
              print "OOps molecule name for molecule %s is %s" %(mol_no_ls,label_str)
    return mols

# Fill an option menu with maps
#
def fill_option_menu_with_map_mol_options(menu):
    fill_option_menu_with_mol_options(menu, valid_map_molecule_qm)

# Helper function for molecule chooser.  Not really for users.
# 
# Return a list of models, corresponding to the menu items of the
# option menu.
# 
# The returned list will not contain references to map or closed
# molecules.
# 
def fill_option_menu_with_coordinates_mol_options(menu):
    fill_option_menu_with_mol_options(menu, valid_model_molecule_qm)

#
def fill_option_menu_with_number_options(menu, number_list, default_option_value):
    count = 0
    while count < len(number_list):
       mlabel_str = str(number_list[count])
       menu.append_text(mlabel_str)
       if (default_option_value == number_list[count]):
          menu.set_active(count)
          print "setting menu active ", default_option_value, count
       count +=1

# Helper function for molecule chooser.  Not really for users.
# 
# return the molecule number of the active item in the option menu,
# or return False if there was a problem (e.g. closed molecule)
#
# BL says:: we do it for gtk_combobox instead! option_menu is deprecated
#
def get_option_menu_active_molecule(option_menu, model_mol_list):

    model = option_menu.get_model()
    active_item = option_menu.get_active()
    # combobox has no children as such, so we just count the rows
    children = 0
    for i in model:
        children += 1

    return model[active_item][0]
#    if (children == model_mol_list):
#       return model[active_item][0]
#    else:
#       print "Failed children length test : ",children, model_mol_list
#       return False

# BL says:: we do it for gtk_combobox instead! option_menu is deprecated
# Here we return the active item in an option menu of generic items
#
def get_option_menu_active_item(option_menu, item_list):

   model = option_menu.get_model()
   active_item = option_menu.get_active()
   # combobox has no children as such, so we just count the rows
   children = 0
   for i in model:
        children += 1

   return model[active_item][0]

def molecule_chooser_gui_generic(chooser_label, callback_function, option_menu_fill_function):
 
    def delete_event(*args):
       window.destroy()
       return False

    def on_ok_clicked(*args):
        # what is the molecule number of the option menu?
        active_mol_no, junk = get_option_menu_active_molecule(option_menu,model_mol_list).split(' ',1)
        try:
           active_mol_no = int(active_mol_no)
           print "INFO: operating on molecule number ", active_mol_no
           callback_function(active_mol_no)
           delete_event()
        except:
           print "Failed to get a (molecule) number"

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    label = gtk.Label(chooser_label)
    vbox = gtk.VBox(False,6)
    hbox_buttons = gtk.HBox(False,5)
# BL says:: option menu is depricated, so we use combox instead, maybe!?!
    option_menu = gtk.combo_box_new_text()
    ok_button = gtk.Button("  OK  ")
    cancel_button = gtk.Button(" Cancel ")
    h_sep = gtk.HSeparator()
    model_mol_list = option_menu_fill_function(option_menu)

    window.set_default_size(370,100)
    window.add(vbox)
    vbox.pack_start(label,False,False,5)
    vbox.pack_start(option_menu,True,True,0)
    vbox.pack_start(h_sep,True,False,2)
    vbox.pack_start(hbox_buttons,False,False,5)
    hbox_buttons.pack_start(ok_button,True,False,5)
    hbox_buttons.pack_start(cancel_button,True,False,5)

    # button callbacks:
    ok_button.connect("clicked",on_ok_clicked,option_menu,model_mol_list)
    cancel_button.connect("clicked",delete_event)

    window.show_all()
#molecule_chooser_gui("test-gui",print_sequence_chain(0,"A"))

# Fire up a molecule chooser dialog, with a given label and on OK we
# call the call_back_fuction with an argument of the chosen molecule
# number. 
# 
# chooser-label is a directive to the user such as "Choose a Molecule"
# 
# callback-function is a function that takes a molecule number as an
# argument.
#
def molecule_chooser_gui(chooser_label,callback_function):
    molecule_chooser_gui_generic(chooser_label,callback_function,fill_option_menu_with_coordinates_mol_options)

# As above but for maps
#
def map_molecule_chooser_gui(chooser_label,callback_function):
    molecule_chooser_gui_generic(chooser_label,callback_function,fill_option_menu_with_map_mol_options)

# A pair of widgets, a molecule chooser and an entry.  The
# callback_function is a function that takes a molecule number and a
# text string.
#
def generic_chooser_and_entry(chooser_label,entry_hint_text,default_entry_text,callback_function):

    import operator

    def delete_event(*args):
       window.destroy()
       return False

    def on_ok_button_clicked(*args):
        # what is the molecule number of the option menu?
        active_mol_no, junk = get_option_menu_active_molecule(option_menu,model_mol_list).split(' ',1)

        try:
           active_mol_no = int(active_mol_no)
           print "INFO: operating on molecule number ", active_mol_no
           text = entry.get_text()
           callback_function(active_mol_no,text)
           delete_event()
        except:
           print "Failed to get a (molecule) number"

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    label = gtk.Label(chooser_label)
    vbox = gtk.VBox(False, 2)
    hbox_for_entry = gtk.HBox(False, 0)
    entry = gtk.Entry()
    entry_label = gtk.Label(entry_hint_text)
    hbox_buttons = gtk.HBox(True, 2)
    option_menu = gtk.combo_box_new_text()
    ok_button = gtk.Button("  OK  ")
    cancel_button = gtk.Button(" Cancel ")
    h_sep = gtk.HSeparator()
    model_mol_list = fill_option_menu_with_coordinates_mol_options(option_menu)
    
    window.set_default_size(400,100)
    window.add(vbox)
    vbox.pack_start(label, False, False, 5)
    vbox.pack_start(option_menu, True, True, 0)
    vbox.pack_start(hbox_for_entry, False, False, 5)
    vbox.pack_start(h_sep, True, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 5)
    hbox_buttons.pack_start(ok_button, True, False, 5)
    hbox_buttons.pack_start(cancel_button, False, False, 5)
    hbox_for_entry.pack_start(entry_label, False, False, 4)
    hbox_for_entry.pack_start(entry, True, True, 4)
    entry.set_text(default_entry_text)

    # button callbacks
    ok_button.connect("clicked",on_ok_button_clicked, entry, option_menu, callback_function)
    cancel_button.connect("clicked", delete_event)

    window.show_all()

# A pair of widgets, a chooser entry and a file selector.  The
# callback_function is a function that takes a molecule number and a
# text string (e.g. chain_id and file_name)
#
# chooser_filter is typically valid_map_molecule_qm or valid_model_molecule_qm
#
def generic_chooser_entry_and_file_selector(chooser_label, chooser_filter, entry_hint_text, default_entry_text, file_selector_hint, callback_function):

    import operator

    def delete_event(*args):
       window.destroy()
       return False

    def on_ok_button_clicked(*args):
        # what is the molecule number of the option menu?
        active_mol_no, junk = get_option_menu_active_molecule(option_menu,model_mol_list).split(' ',1)

        try:
           active_mol_no = int(active_mol_no)
           text = entry.get_text()
           file_sel_text = file_sel_entry.get_text()
           callback_function(active_mol_no, text, file_sel_text)
           delete_event()
        except:
           print "Failed to get a (molecule) number"

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    label = gtk.Label(chooser_label)
    vbox = gtk.VBox(False, 2)
    hbox_for_entry = gtk.HBox(False, 0)
    entry = gtk.Entry()
    entry_label = gtk.Label(entry_hint_text)
    hbox_buttons = gtk.HBox(True, 2)
    option_menu = gtk.combo_box_new_text()
    ok_button = gtk.Button("  OK  ")
    cancel_button = gtk.Button(" Cancel ")
    h_sep = gtk.HSeparator()
    model_mol_list = fill_option_menu_with_mol_options(option_menu, chooser_filter)
    
    window.set_default_size(400,100)
    window.add(vbox)
    vbox.pack_start(label, False, False, 5)
    vbox.pack_start(option_menu, True, True, 0)
    vbox.pack_start(hbox_for_entry, False, False, 5)
    hbox_buttons.pack_start(ok_button, True, False, 5)
    hbox_buttons.pack_start(cancel_button, False, False, 5)
    hbox_for_entry.pack_start(entry_label, False, False, 4)
    hbox_for_entry.pack_start(entry, True, True, 4)
    entry.set_text(default_entry_text)

    file_sel_entry = file_selector_entry(vbox, file_selector_hint)
    vbox.pack_start(h_sep, True, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 5)
    

    # button callbacks
    ok_button.connect("clicked",on_ok_button_clicked, entry, option_menu, callback_function)
    cancel_button.connect("clicked", delete_event)

    window.show_all()

# A pair of widgets, a molecule chooser and a file selector.  The
# callback_function is a function that takes a molecule number and a
# file_name
#
# chooser_filter is typically valid_map_molecule_qm or valid_model_molecule_qm
#
def generic_chooser_and_file_selector(chooser_label, chooser_filter, file_selector_hint, default_file_name, callback_function):

    import operator

    def delete_event(*args):
       window.destroy()
       return False

    def on_ok_button_clicked(*args):
        # what is the molecule number of the option menu?
        active_mol_no, junk = get_option_menu_active_molecule(option_menu,model_mol_list).split(' ',1)

        try:
           active_mol_no = int(active_mol_no)
           file_sel_text = file_sel_entry.get_text()
           callback_function(active_mol_no, file_sel_text)
           delete_event()
        except:
           print "Failed to get a (molecule) number"

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    label = gtk.Label(chooser_label)
    vbox = gtk.VBox(False, 2)
    hbox_for_entry = gtk.HBox(False, 0)
    hbox_buttons = gtk.HBox(True, 2)
    option_menu = gtk.combo_box_new_text()
    ok_button = gtk.Button("  OK  ")
    cancel_button = gtk.Button(" Cancel ")
    h_sep = gtk.HSeparator()
    model_mol_list = fill_option_menu_with_mol_options(option_menu, chooser_filter)
    
    window.set_default_size(400,100)
    window.add(vbox)
    vbox.pack_start(label, False, False, 5)
    vbox.pack_start(option_menu, True, True, 0)
    hbox_buttons.pack_start(ok_button, True, False, 5)
    hbox_buttons.pack_start(cancel_button, False, False, 5)

    file_sel_entry = file_selector_entry(vbox, file_selector_hint)
    vbox.pack_start(h_sep, True, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 5)
    

    # button callbacks
    ok_button.connect("clicked",on_ok_button_clicked, option_menu, callback_function)
    cancel_button.connect("clicked", delete_event)

    window.show_all()


# If a menu with label menu-label is not found in the coot main
# menubar, then create it and return it. 
# If it does exist, simply return it.
#  
def coot_menubar_menu(menu_label):

   try:
    import coot_python
    coot_main_menubar = coot_python.main_menubar()

    def menu_bar_label_list():
      ac_lab_ls = []
      ac_lab = []
      for menu_child in coot_main_menubar.get_children():
          lab = []
          # ac-lab-ls is a GtkAccelLabel in a list
          lab.append(menu_child.get_children()[0].get_text())
          lab.append(menu_child)
          ac_lab_ls.append(lab)
#      print "BL DEBUG:: list is", ac_lab_ls, ac_lab, lab
      return ac_lab_ls

    # main body
    #
    found_menu = False
    for f in menu_bar_label_list():
       if menu_label in f: 
          # print "BL DEBUG:: found menu label is ", f
          found_menu = f[1]
    if found_menu:
       return found_menu
    else:
       menu = gtk.Menu()
       menuitem = gtk.MenuItem(menu_label)
       menuitem.set_submenu(menu)
       coot_main_menubar.append(menuitem)
       menuitem.show()
       return menu
   except: print """BL WARNING:: could not import coot_python module!!\n
                    Some things, esp. extensions, may be crippled!"""

# Given that we have a menu (e.g. one called "Extensions") provide a
# cleaner interface to adding something to it:
# 
# activate-function is a thunk.
#
def add_simple_coot_menu_menuitem(menu,menu_item_label,activate_function):

    submenu = gtk.Menu()
    sub_menuitem = gtk.MenuItem(menu_item_label)

    menu.append(sub_menuitem)
    sub_menuitem.show()

    sub_menuitem.connect("activate",activate_function)


# Make an interesting things GUI for residues of molecule number
# imol that have alternate conformations.
#
def alt_confs_gui(imol):

	if valid_model_molecule_qm(imol):
		alt_conf_residues = residues_with_alt_confs(imol)
		for i in range(len(alt_conf_residues)): alt_conf_residues[i][0] = imol
		centre_atoms = map(residue_spec, alt_conf_residues)
		if centre_atoms:
			interesting_things_gui("Residues with Alt Confs",
			map(lambda alt_conf_residue_cpmd, centre_atom: 
			[alt_conf_residue_cpmd[1] + " " + str(alt_conf_residue_cpmd[2]) + " " + alt_conf_residue_cpmd[3] + " " + centre_atom[0] + " " + centre_atom[1]] 
			+ alt_conf_residue_cpmd + centre_atom
			, alt_conf_residues, centre_atoms))
		else:
			print "BL INFO:: no alternative confs detected"
	else:
		print "BL WARNING:: no valid molecule", imol

# button-list is a list of pairs (improper list) the first item of
# which is the button label text the second item is a lambda
# function, what to do when the button is pressed.
#
def generic_button_dialog(dialog_name, button_list):

	def delete_event(*args):
		window.destroy()
		return False

	# main body
	window = gtk.Window(gtk.WINDOW_TOPLEVEL)
	scrolled_win = gtk.ScrolledWindow()
	outside_vbox = gtk.VBox(False, 2)
	inside_vbox = gtk.VBox(False, 0)

	window.set_default_size(250,250)
	window.set_title(dialog_name)
	inside_vbox.set_border_width(4)

	window.add(outside_vbox)
	outside_vbox.add(scrolled_win)
	scrolled_win.add_with_viewport(inside_vbox)
	scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)

	for button_item in button_list:
		if button_item and len(button_item)==2:
			button_label = button_item[0]
			action = button_item[1]

			button = gtk.Button(button_label)
			inside_vbox.pack_start(button, False, False, 2)
			button.connect("clicked", action)
			button.show()

	outside_vbox.set_border_width(4)
	ok_button = gtk.Button("  OK  ")
	outside_vbox.pack_start(ok_button, False, False, 6)
	ok_button.connect("clicked", delete_event)

	window.show_all()


# Generic interesting things gui: user passes a function that takes 4
# args: the chain-id, resno, inscode and residue-serial-number
# (should it be needed) and returns either #f or something
# interesting (e.g. a label/value).  It is the residue-test-func of
# the residue-matching-criteria function.
# 
def generic_interesting_things(imol,gui_title_string,residue_test_func):

	if valid_model_molecule_qm(imol):

		interesting_residues = residues_matching_criteria(imol, residue_test_func)
		for i in range(len(interesting_residues)): interesting_residues[i][0] = imol
		centre_atoms = map(residue_spec, interesting_residues)
		if centre_atoms:
			# BL says:: ignoring "Atom in residue name failure" for nor
			interesting_things_gui(gui_title_string,
			map(lambda interesting_residue, centre_atom:
			[interesting_residue[0] + " " 
			+ interesting_residue[1] + " " 
			+ str(interesting_residue[2]) + " " 
			+ interesting_residue[3] + " " 
			+ centre_atom[0] + " " + centre_atom[1]]
			+ interesting_residue[1:len(interesting_residue)] + centre_atom,
			interesting_residues, centre_atoms)
			)
	else:
		print "BL WARNING:: no valid model molecule ", imol

def generic_number_chooser(number_list, default_option_value, hint_text, go_button_label, go_function):

    def delete_event(*args):
       window.destroy()
       return False

    def go_button_pressed(*args):
        active_number = int(get_option_menu_active_molecule(option_menu, number_list))
        try:
#           print "BL DEBUG:: go_function is:", go_function
#           print "BL DEBUG:: active_number is:", active_number
           go_function(active_number)
        except:
           print "Failed to get execute function"
        delete_event()


    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    vbox = gtk.VBox(False, 0)
    hbox1 = gtk.HBox(False, 0)
    hbox2 = gtk.HBox(True, 0)      # for Go and Cancel
    function_label = gtk.Label(hint_text)
    h_sep = gtk.HSeparator()
    go_button = gtk.Button(go_button_label)
    cancel_button = gtk.Button("  Cancel  ")
# BL says:: option menu is depricated, so we use combox instead, maybe!?!
    option_menu = gtk.combo_box_new_text()

    fill_option_menu_with_number_options(option_menu, number_list, default_option_value)

    vbox.pack_start(hbox1, True, False, 0)
    vbox.pack_start(function_label, False, 0)
    vbox.pack_start(option_menu, True, 0)
    vbox.pack_start(h_sep)
    vbox.pack_start(hbox2, False, False, 0)
    hbox2.pack_start(go_button, False, False, 6)
    hbox2.pack_start(cancel_button, False, False, 6)
    window.add(vbox)
    vbox.set_border_width(6)
    hbox1.set_border_width(6)
    hbox2.set_border_width(6)
    go_button.connect("clicked", go_button_pressed, option_menu, number_list, go_function)
    cancel_button.connect("clicked", delete_event)

    window.show_all()

# vbox is the vbox to which this compound widget should be added.
# button-press-func is the lambda function called on pressing return
# or the button, which takes one argument, the entry.
#
# Add this widget to the vbox here.
#
def entry_do_button(vbox, hint_text, button_label, button_press_func):

	hbox = gtk.HBox(False, 0)
	entry = gtk.Entry()
	button = gtk.Button(button_label)
	label = gtk.Label(hint_text)

	hbox.pack_start(label, False, False, 2)
	hbox.pack_start(entry, True, False, 2)
	hbox.pack_start(button, False, False, 2)
	button.connect("clicked", button_press_func, entry)

	label.show()
	button.show()
	entry.show()
	hbox.show()
	vbox.pack_start(hbox, True, False)
	return entry

# pack a hint text and a molecule chooser option menu into the given vbox.
# 
# return the option-menu and model molecule list:
def generic_molecule_chooser(hbox, hint_text):

	menu = gtk.Menu()
# BL says:: option menu is depricated, so we use combox instead, maybe!?!
	option_menu = gtk.combo_box_new_text()
	label = gtk.Label(hint_text)
	model_mol_list = fill_option_menu_with_coordinates_mol_options(option_menu)

	hbox.pack_start(label, False, False, 2)
	hbox.pack_start(option_menu, True, True, 2)

# Return an entry, insert the widget into the hbox in this function
#
def file_selector_entry(hbox, hint_text):

	vbox = gtk.VBox(False, 0)

	def file_func1(*args):
		def file_ok_sel(*args):
			t = fs_window.get_filename()
			print t
			entry.set_text(t)
			fs_window.destroy()

		fs_window = gtk.FileSelection("file selection")
		fs_window.ok_button.connect("clicked", file_ok_sel)
		fs_window.cancel_button.connect("clicked",
				lambda w: fs_window.destroy())
		fs_window.show()

	entry = entry_do_button(vbox, hint_text, "  File...  ", file_func1)

	hbox.pack_start(vbox, False, False, 2)
	vbox.show()
        return entry

# This is the same as the file_selector_entry, but using the modern FileChooser
# Return an entry, insert the widget into the hbox in this function
#
def file_chooser_entry(hbox, hint_text):

   if gtk.pygtk_version > (2,3,90):

	vbox = gtk.VBox(False, 0)

	def file_func1(*args):
		def file_ok_sel(*args):
			t = fc_window.get_filename()
			print t
			entry.set_text(t)
			fc_window.destroy()

		fc_window = gtk.FileChooserDialog("file selection",
                                                  None,
                                                  gtk.FILE_CHOOSER_ACTION_OPEN,
                                                  (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                                   gtk.STOCK_OPEN, gtk.RESPONSE_OK))
                response = fc_window.run()
                if response == gtk.RESPONSE_OK:
                   file_ok_sel(fc_window, entry)
                elif response == gtk.RESPONSE_CANCEL:
                   fc_window.destroy()

	entry = entry_do_button(vbox, hint_text, "  File...  ", file_func1)

	hbox.pack_start(vbox, False, False, 2)
	vbox.show()
        return entry
   else:
        print "PyGtk 2.3.90 or later required for this function!"
        return False

# The gui for the strand placement function
#
def place_strand_here_gui():

    generic_number_chooser(number_list(4,12), 7, 
                            " Estimated number of residues in strand",
                            "  Go  ",
                            lambda n: place_strand_here(n, 15))

# Cootaneer gui
def cootaneer_gui(imol):

	def delete_event(*args):
		window.destroy()
		return False

	def go_function_event(widget, imol):
		print "apply the sequence info here\n"
		print "then cootaneer\n"

		# no active atom won't do.  We need
		# to find the nearest atom in imol to (rotation-centre).
		#
		# if it is too far away, give a
		# warning and do't do anything.

		n_atom = closest_atom(imol)
		if n_atom:
			imol	= n_atom[0]
			chain_id= n_atom[1]
			resno	= n_atom[2]
			inscode	= n_atom[3]
			at_name	= n_atom[4]
			alt_conf= n_atom[5]
			cootaneer(imol_map, imol, [chain_id, resno, inscode, 
				at_name, alt_conf])
		else:
			print "BL WARNING:: no close atom found!"


	def add_text_to_text_box(text_box, description):
		start = text_box.get_start_iter()
		text_box.create_tag("tag", foreground="black", 
			background = "#c0e6c0")
		text_box.insert_with_tags_by_name(start, description, "tag")

	# return the (entry . textbuffer/box)
	#
	def entry_text_pair_frame(seq_info):

		frame = gtk.Frame()
		vbox = gtk.VBox(False, 3)
		entry = gtk.Entry()
		textview = gtk.TextView()
                textview.set_wrap_mode(gtk.WRAP_WORD_CHAR)
		text_box = textview.get_buffer()
		chain_id_label = gtk.Label("Chain ID")
		sequence_label = gtk.Label("Sequence")
	
		frame.add(vbox)
		vbox.pack_start(chain_id_label, False, False, 2)
		vbox.pack_start(entry, False, False, 2)
		vbox.pack_start(sequence_label, False, False, 2)
		vbox.pack_start(textview, False, False, 2)
		add_text_to_text_box(text_box, seq_info[1])
		entry.set_text(seq_info[0])
		return [frame, entry, text_box]

	# main body
	imol_map = imol_refinement_map()
	if (imol_map == -1):
           show_select_map_dialog()
           print "BL DEBUG:: probably should wait here for input!?"
	
	window = gtk.Window(gtk.WINDOW_TOPLEVEL)
	outside_vbox = gtk.VBox(False, 2)
	inside_vbox = gtk.VBox(False, 2)
	h_sep = gtk.HSeparator()
	buttons_hbox = gtk.HBox(True, 2)
	go_button = gtk.Button("  Cootaneer!  ")
	cancel_button = gtk.Button("  Cancel  ")

	seq_info_ls = sequence_info(imol)
        # print "BL DEBUG:: sequence_list and imol is", seq_info_ls, imol

	for seq_info in seq_info_ls:
		seq_widgets = entry_text_pair_frame(seq_info)
		inside_vbox.pack_start(seq_widgets[0], False, False, 2)

	outside_vbox.pack_start(inside_vbox, False, False, 2)
	outside_vbox.pack_start(h_sep, False, False, 2)
	outside_vbox.pack_start(buttons_hbox, True, False, 2)
	buttons_hbox.pack_start(go_button, True, False, 6)
	buttons_hbox.pack_start(cancel_button, True, False, 6)

	cancel_button.connect("clicked", delete_event)

	go_button.connect("clicked", go_function_event, imol)

        window.add(outside_vbox)
	window.show_all()	


# The gui for saving views
def view_saver_gui():

	def local_view_name():
		view_count = 0
		while view_count or view_count == 0:
			strr = "View"
			if view_count > 0:
				strr = strr + "-" + str(view_count)
			# now is a view already called str?
			jview = 0
			while jview or jview == 0:
				jview_name = view_name(jview)
				if jview >= n_views(): return strr
				elif jview_name == False: return strr
				elif strr == jview_name:
					view_count += 1
					break
				else:
					jview +=1
		return strr

	generic_single_entry("View Name: ", local_view_name(), " Add View ",
		lambda text: add_view_here(text))

# geometry is an improper list of ints
#
def dialog_box_of_buttons(window_name, geometry, buttons, close_button_label):

	def add_text_to_text_widget(text_box, description):
		textbuffer = text_box.get_buffer()
		start = textbuffer.get_start_iter()
		textbuffer.create_tag("tag", foreground="black", 
						background = "#c0e6c0")
		textbuffer.insert_with_tags_by_name(start, description, "tag")

	# main line
	window = gtk.Window(gtk.WINDOW_TOPLEVEL)
	scrolled_win = gtk.ScrolledWindow()
	outside_vbox = gtk.VBox(False, 2)
	inside_vbox = gtk.VBox(False, 0)

        window.set_default_size(geometry[0], geometry[1])
	window.set_title(window_name)
	inside_vbox.set_border_width(2)
	window.add(outside_vbox)
	outside_vbox.pack_start(scrolled_win, True, True, 0) # expand fill padding
	scrolled_win.add_with_viewport(inside_vbox)
	scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)

	for button_info in buttons:
		button_label = button_info[0]
		callback = button_info[1]
		if len(button_info)==2:
			description = False
		else:
			description = button_info[2]
		button = gtk.Button(button_label)

# BL says:: in python we should pass the callback as a string
		if type(callback) is StringType:
			def callback_func(button,call):
				eval(call)
			button.connect("clicked", callback_func, callback)
		else:
			button.connect("clicked", callback)

		if type(description) is StringType:
			text_box = gtk.TextView()
			text_box.set_editable(False)
			add_text_to_text_widget(text_box, description)
			inside_vbox.pack_start(text_box, False, False, 2)
			text_box.realize()
# BL says:: not working here
#			text_box.thaw()

		inside_vbox.pack_start(button, False, False, 2)

	outside_vbox.set_border_width(2)
	ok_button = gtk.Button(close_button_label)
	outside_vbox.pack_end(ok_button, False, False, 0)
	ok_button.connect("clicked", lambda w: window.destroy())
	
	window.show_all()

# geometry is an improper list of ints
# buttons is a list of: [["button_1_label, button_1_action],
#                        ["button_2_label, button_2_action]]
# The button_1_action function takes no arguments.
# The button_2_action function takes as an argument the imol
# 
def dialog_box_of_pairs_of_buttons(imol, window_name, geometry, buttons, close_button_label):

        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        scrolled_win = gtk.ScrolledWindow()
        outside_vbox = gtk.VBox(False, 2)
        inside_vbox = gtk.VBox(False, 0)

        window.set_default_size(geometry[0], geometry[1])
        window.set_title(window_name)
        inside_vbox.set_border_width(2)
        window.add(outside_vbox)
        outside_vbox.pack_start(scrolled_win, True, True, 0) # expand fill padding
        scrolled_win.add_with_viewport(inside_vbox)
        scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)

	for button_info in buttons:
		print "buttons_info ", buttons_info
		if type(button_info) is ListType:
			buttons_label_1 = button_info[0]
			callback_1 = button_info[1]

			buttons_label_2 = button_info[2]
			callback_2 = button_info[3]

			button_1 = gtk.Button(button_label_1)
			hbox = gtk.HBox(False, 2)

			print "button_label_1 ", button_label_1
			print "callback_1 ", callback_1
			print "button_label_2 ", button_label_2
			print "callback_2 ", callback_2

			# But what is imol??
			button_1.connect("clicked", callback_1)
			h_box.pack_start(button_1, False, False, 2)

			if callback_2:
				button_2 = gtk.Button(button_label_2)
				button_2.connect("clicked",
					lambda func: callback_2(imol))
			inside_vbox.pack_start(h_box, False, False, 2)

		outside_vbox.set_border_width(2)
		ok_button = gtk.Button(close_button_label)
		outside_vbox.pack_end(ok_button, False, False, 2)
		ok_button.connect("clicked", lambda w: window.destroy())
		window.show_all()

# as the dialog_box_of_buttons, butwe can put in an extra widget (extra_widget)
#
def dialog_box_of_buttons_with_widget(window_name, geometry, buttons, extra_widget, close_button_label):

	def add_text_to_text_widget(text_box, description):
		textbuffer = text_box.get_buffer()
		start = textbuffer.get_start_iter()
		textbuffer.create_tag("tag", foreground="black", 
						background = "#c0e6c0")
		textbuffer.insert_with_tags_by_name(start, description, "tag")

	# main line
	window = gtk.Window(gtk.WINDOW_TOPLEVEL)
	scrolled_win = gtk.ScrolledWindow()
	outside_vbox = gtk.VBox(False, 2)
	inside_vbox = gtk.VBox(False, 0)
        h_sep = gtk.HSeparator()

        window.set_default_size(geometry[0], geometry[1])
	window.set_title(window_name)
	inside_vbox.set_border_width(2)
	window.add(outside_vbox)
	outside_vbox.pack_start(scrolled_win, True, True, 0) # expand fill padding
	scrolled_win.add_with_viewport(inside_vbox)
	scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)

	for button_info in buttons:
		button_label = button_info[0]
		callback = button_info[1]
		if len(button_info)==2:
			description = False
		else:
			description = button_info[2]
		button = gtk.Button(button_label)

# BL says:: in python we should pass the callback as a string
		if type(callback) is StringType:
			def callback_func(button,call):
				eval(call)
			button.connect("clicked", callback_func, callback)
		else:
			button.connect("clicked", callback)

		if type(description) is StringType:
			text_box = gtk.TextView()
			text_box.set_editable(False)
			add_text_to_text_widget(text_box, description)
			inside_vbox.pack_start(text_box, False, False, 2)
			text_box.realize()
# BL says:: not working here
#			text_box.thaw()

		inside_vbox.pack_start(button, False, False, 2)

        # for the extra widget
        inside_vbox.pack_start(h_sep, False, False, 2)
        inside_vbox.pack_start(extra_widget, False, False, 2)
        
	outside_vbox.set_border_width(2)
	ok_button = gtk.Button(close_button_label)
	outside_vbox.pack_end(ok_button, False, False, 0)
	ok_button.connect("clicked", lambda w: window.destroy())
	
	window.show_all()

# A gui showing views
def views_panel_gui():

	number_of_views = n_views()
	buttons = []

	for button_number in range(number_of_views):
		button_label = view_name(button_number)
		desciption = view_description(button_number)
# BL says:: add the decisption condition!!
		buttons.append([button_label, "go_to_view_number(" + str(button_number) + ",0)", desciption])

	if len(buttons) > 1:
		def view_button_func():
			import time
			go_to_first_view(1)
			time.sleep(1)
			play_views()
		view_button = ["  Play Views ", lambda func: view_button_func()]
		buttons.insert(0,view_button)

	dialog_box_of_buttons("Views", [200,140], buttons, "  Close  ")

# nudge screen centre box.  Useful when Ctrl left-mouse has been
# taken over by another function.
#
# This is using the absolute coordinates
# 
def nudge_screen_centre_paule_gui():

	zsc = 0.02
# BL says:: in python we need some helper functn, no 'proper' lambda
# axis 0-2 = x,y,z and 0=+ and 1=-
	def nudge_screen_func(axis, operand):
		nudge = zsc * zoom_factor()
		rc = rotation_centre()
		if operand == 0:
			rc[axis] += nudge
		elif operand == 1:
			rc[axis] -= nudge
		else:
			# We should never be here
			print "BL WARNING:: something went wrong!!!"
		set_rotation_centre(*rc)

	buttons = [ 
		["Nudge +X", lambda func: nudge_screen_func(0,0)],
		["Nudge -X", lambda func: nudge_screen_func(0,1)],
		["Nudge +Y", lambda func: nudge_screen_func(1,0)],
		["Nudge -Y", lambda func: nudge_screen_func(1,1)],
		["Nudge +Z", lambda func: nudge_screen_func(2,0)],
		["Nudge -Z", lambda func: nudge_screen_func(2,1)],
	]

	dialog_box_of_buttons("Nudge Screen Centre", [200,240], buttons, "  Close ")
		

# nudge screen centre box.  Useful when Ctrl left-mouse has been
# taken over by another function.
#
# This is using the viewing coordinates
# 
def nudge_screen_centre_gui():

	zsc = 0.02
# BL says:: in python we need some helper functn, no 'proper' lambda
# axis 0-2 = x,y,z and 0=+ and 1=-
	def nudge_screen_func(axis, operand):
		nudge_factor = zsc * zoom_factor()
		rc = rotation_centre()
                mat = view_matrix_transp()
                nudge_vector = []
                for i in range(axis, axis + 7, 3):
                   nudge_vector.append(mat[i])
		if operand == 0:
                   rc = map(lambda x, y: x + y, rc, nudge_vector)
		elif operand == 1:
                   rc = map(lambda x, y: x - y, rc, nudge_vector)
		else:
			# We should never be here
			print "BL WARNING:: something went wrong!!!"
		set_rotation_centre(*rc)

	buttons = [ 
		["Nudge +X", lambda func: nudge_screen_func(0,0)],
		["Nudge -X", lambda func: nudge_screen_func(0,1)],
		["Nudge +Y", lambda func: nudge_screen_func(1,0)],
		["Nudge -Y", lambda func: nudge_screen_func(1,1)],
		["Nudge +Z", lambda func: nudge_screen_func(2,0)],
		["Nudge -Z", lambda func: nudge_screen_func(2,1)],
	]

	dialog_box_of_buttons("Nudge Screen Centre", [200,240], buttons, "  Close ")

# as nudge_screen_centre_gui but with clipping and zoom control
def nudge_screen_centre_extra_gui():

   # this is for the centre nudging
   zsc = 0.02
# BL says:: in python we need some helper functn, no 'proper' lambda
# axis 0-2 = x,y,z and 0=+ and 1=-
   def nudge_screen_func(axis, operand):
	nudge_factor = zsc * zoom_factor()
	rc = rotation_centre()
        mat = view_matrix_transp()
        nudge_vector = []
        for i in range(axis, axis + 7, 3):
           nudge_vector.append(mat[i])
        if operand == 0:
           rc = map(lambda x, y: x + y, rc, nudge_vector)
        elif operand == 1:
           rc = map(lambda x, y: x - y, rc, nudge_vector)
        else:
           # We should never be here
           print "BL WARNING:: something went wrong!!!"
        set_rotation_centre(*rc)

   buttons = [ 
		["Nudge +X", lambda func: nudge_screen_func(0,0)],
		["Nudge -X", lambda func: nudge_screen_func(0,1)],
		["Nudge +Y", lambda func: nudge_screen_func(1,0)],
		["Nudge -Y", lambda func: nudge_screen_func(1,1)],
		["Nudge +Z", lambda func: nudge_screen_func(2,0)],
		["Nudge -Z", lambda func: nudge_screen_func(2,1)],
	]

   # and this is for the clipping and zooming
   vbox = gtk.VBox(False, 0)

   def change_clipp(*args):
        set_clipping_front(clipp_adj.value)
        set_clipping_back (clipp_adj.value)

   def change_zoom(*args):
        set_zoom(zoom_adj.value)
        graphics_draw()

   # for clipping
   clipp_label = gtk.Label("Clipping")
   clipp_adj = gtk.Adjustment(0.0, -10.0, 20.0, 0.05, 4.0, 10.1)
   clipp_scale = gtk.HScale(clipp_adj)
   vbox.pack_start(clipp_label, False, False, 0)
   vbox.pack_start(clipp_scale, False, False, 0)
   clipp_label.show()
   clipp_scale.show()

   clipp_adj.connect("value_changed", change_clipp)
   
   h_sep = gtk.HSeparator()
   vbox.pack_start(h_sep, False, False, 5)

   # for zooming
   zoom = zoom_factor()
   zoom_label = gtk.Label("Zoom")
   zoom_adj = gtk.Adjustment(zoom, zoom*0.125, zoom*8, 0.01, 0.5, zoom)
   zoom_scale = gtk.HScale(zoom_adj)
   vbox.pack_start(zoom_label, False, False, 0)
   vbox.pack_start(zoom_scale, False, False, 0)
   zoom_label.show()
   zoom_scale.show()

   zoom_adj.connect("value_changed", change_zoom)

   dialog_box_of_buttons_with_widget("Nudge Screen Centre with Extras", [200,400], buttons, vbox, "  Close ")


def cis_peptides_gui(imol):

   def get_ca(atom_list):

      if (atom_list == []):
         return False
      else:
         for atom in atom_list:
            atom_name = atom[0][0]
            if (atom_name == ' CA '):
               print "BL DEBUG:: returning ", atom
               return atom     # check me

   def make_list_of_cis_peps(imol, list_of_cis_peps):

      ret = []

      for cis_pep_spec in list_of_cis_peps:
         r1 = cis_pep_spec[0]
         r2 = cis_pep_spec[1]
         omega = cis_pep_spec[2]
         atom_list_r1 = residue_info(imol, *r1[1:4])
         atom_list_r2 = residue_info(imol, *r2[1:4])
         ca_1 = get_ca(atom_list_r1)
         ca_2 = get_ca(atom_list_r2)
         chain_id = r1[1]

         if ((not ca_1) and (not ca_2)):
            ret.append(["Cis Pep (atom failure) " + r1[1] + " " + str(r1[2]),
                    imol, chain_id, r1[3], r1[2], "CA", ""])
         else:
            p1 = ca_1[2]
            p2 = ca_2[2]
            pos = map(lambda x, y: (x + y) / 2.0, p1, p2)
            tors_s1 = str(omega)
            if (len(tors_s1) < 6):
               tors_string = tors_s1
            else:
               tors_string = tors_s1[0:6]
            mess = ("Cis Pep: " + chain_id + " " +
                   str(r1[2]) + " " + 
                   residue_name(imol, *r1[1:4]) + " - " +
                   str(r2[2]) + " " +
                   residue_name(imol, *r1[1:4]) + "   " +
                   tors_string)
            ls = pos
            ls.insert(0,mess)
            ret.append(ls)
      return ret

   cis_peps = cis_peptides(imol)

   if (cis_peps == []):
      info_dialog("No Cis Peptides found")
   else:
      list_of_cis_peptides = make_list_of_cis_peps(imol, cis_peps)
      interesting_things_gui("Cis Peptides:", list_of_cis_peptides)

# Cootaneer/sequencing gui modified by BL with ideas from KC
# based on Paul's cootaneer gui and generic_chooser_entry_and_file_selector
def cootaneer_gui_bl():

	# unfortunately currently I dont see a way to avoid a global variable here
	# contains flag if file was imported and no of sequences
	global imported_sequence_file_flags
	imported_sequence_file_flags = [False, 0]

	def delete_event(*args):
		window.destroy()
		return False

	def refine_function_event(widget):
		# doesnt need to do anything
		status = refine_check_button.get_active()
		if (status):
			print "INFO:: refinement on"
		else:
			print "INFO:: refinement off"

	def go_function_event(widget):
		print "apply the sequence info here\n"
		print "then cootaneer\n"

		# no active atom won't do.  We need
		# to find the nearest atom in imol to (rotation-centre).
		#
		# if it is too far away, give a
		# warning and do't do anything.
		active_mol_no, junk = get_option_menu_active_molecule(option_menu, model_mol_list).split(' ', 1)
		imol = int(active_mol_no)
		imol_map = imol_refinement_map()

		do_it = assign_sequences_to_mol(imol)

		if (do_it):
			# now cootaneer it
			chain_ls = chain_ids(imol)
			for chain_id in chain_ls:
				res_name = resname_from_serial_number(imol, chain_id, 0)
				res_no = seqnum_from_serial_number(imol, chain_id, 0)
				ins_code = insertion_code_from_serial_number(imol, chain_id, 0)
				alt_conf = ""
				at_name = residue_spec2atom_for_centre(imol, chain_id, res_no, ins_code)[0]
				cootaneer(imol_map, imol, [chain_id, res_no, ins_code, 
							   at_name, alt_conf])
				
				refine_qm = refine_check_button.get_active()
			# refine?
			window.hide()
			if (refine_qm):
				fit_protein(imol)
					
			delete_event()


        def import_function_event(widget, selector_entry):
		# we import a sequence file and update the cootaneer table
		global imported_sequence_file_flags
		imported_sequence_file_qm = imported_sequence_file_flags[0]
		active_mol_no, junk = get_option_menu_active_molecule(option_menu, model_mol_list).split(' ', 1)
		imol = int(active_mol_no)
		
		seq_info_ls = []
		seq_file_name = selector_entry.get_text()
		if (seq_file_name):
		    # get and set sequence info
		    assign_sequence_from_file(imol, str(seq_file_name))
		    seq_info_ls = sequence_info(imol)
		    no_of_sequences = len(seq_info_ls)

		    # remove children if new file
		    if not imported_sequence_file_qm:
			    table_children = seq_table.get_children()
			    for child in table_children:
				    seq_table.remove(child)
			    widget_range = range(no_of_sequences)
		    else:
			    # we update the number of sequences
			    spin_len = int(spin_button.get_value())
			    widget_range = range(spin_len, no_of_sequences)
			    
		    # make new table
		    imported_sequence_file_flags = [True, no_of_sequences]
		    spin_button.set_value(no_of_sequences)
		    seq_table.resize(no_of_sequences, 1)
		    for i in widget_range:
			    seq_widget = entry_text_pair_frame_with_button(seq_info_ls[i])
			    seq_table.attach(seq_widget[0], 0, 1, i, i+1)
			    seq_widget[0].show_all()
		else:
			print "BL WARNING:: no filename"

	def fill_table_with_sequences(*args):
		# fills the table with sequences if they have been associated with the model imol
		# already
		global imported_sequence_file_flags
		active_mol_no, junk = get_option_menu_active_molecule(option_menu, model_mol_list).split(' ', 1)
		imol = int(active_mol_no)
		seq_info_ls = sequence_info(imol)
		if (seq_info_ls):
			# we have a sequence and fill the table
			no_of_sequences = len(seq_info_ls)
			imported_sequence_file_flags = [True, no_of_sequences]
			spin_button.set_value(no_of_sequences)
			seq_table.resize(no_of_sequences, 1)
			# remove existing children
			table_children = seq_table.get_children()
			for child in table_children:
				seq_table.remove(child)
			widget_range = range(no_of_sequences)
			for i in widget_range:
				seq_widget = entry_text_pair_frame_with_button(seq_info_ls[i])
				seq_table.attach(seq_widget[0], 0, 1, i, i+1)
				seq_widget[0].show_all()
		else:
			# no sequence information, reset the table
			clear_function_event()
			


	def add_text_to_text_buffer(text_buffer, description):
		start = text_buffer.get_start_iter()
		text_buffer.create_tag("tag", foreground="black", 
			background = "#c0e6c0")
		text_buffer.insert_with_tags_by_name(start, description, "tag")

	# return the (entry . textbuffer/box)
	#
	def entry_text_pair_frame_with_button(seq_info):
           
           def fragment_go_event(widget):
		active_mol_no, junk = get_option_menu_active_molecule(option_menu, model_mol_list).split(' ', 1)
		imol = int(active_mol_no)
		imol_map = imol_refinement_map()
		print "apply the sequence info here\n"
		print "then cootaneer\n"

		# no active atom won't do.  We need
		# to find the nearest atom in imol to (rotation-centre).
		#
		# if it is too far away, give a
		# warning and do't do anything.

		do_it = assign_sequences_to_mol(imol)
		if (do_it):

			n_atom = closest_atom(imol)
			if n_atom:
				imol	= n_atom[0]
				chain_id= n_atom[1]
				res_no	= n_atom[2]
				ins_code= n_atom[3]
				at_name	= n_atom[4]
				alt_conf= n_atom[5]
				cootaneer(imol_map, imol, [chain_id, res_no, ins_code, 
							   at_name, alt_conf])
				refine_qm = refine_check_button.get_active()
				if (chain_check_button.get_active()):
					# we try to change the chain
					from_chain_id = chain_id
					to_chain_id = entry.get_text()
					start_resno = seqnum_from_serial_number(imol, from_chain_id, 0)
					end_resno = seqnum_from_serial_number(imol, from_chain_id, (chain_n_residues(from_chain_id, imol) - 1))
					[istat, message] = change_chain_id_with_result_py(imol, from_chain_id, to_chain_id, 1, start_resno, end_resno)
					if (istat == 1):
						# renaming ok
						chain_id = to_chain_id
					else:
						# renaming failed
						if (refine_qm):
							message += "\nRefinement proceeded with old=new chain "
							message += chain_id
						info_dialog(message)
					
				if (refine_qm):
					fit_chain(imol, chain_id)
			else:
				print "BL WARNING:: no close atom found!"
		else:
			print "BL ERROR:: something went wrong assigning the sequence"

	   def chain_toggled(widget):
		   # doesnt need to do anything
		   status = chain_check_button.get_active()
		   if(status):
			   print "BL INFO:: assign chain_id too"
		   else:
			   print "BL INFO:: do not assign chain_id"
               
           frame = gtk.Frame()
           vbox = gtk.VBox(False, 3)
           hbox = gtk.HBox(False, 3)
           entry = gtk.Entry()
           textview = gtk.TextView()
           textview.set_wrap_mode(gtk.WRAP_WORD_CHAR)
	   textview.set_editable(True)
	   textview.set_size_request(300, -1)
	   textview.modify_font(pango.FontDescription("Courier 11"))
           text_buffer = textview.get_buffer()
           chain_id_label = gtk.Label("Chain ID")
           sequence_label = gtk.Label("Sequence")
	   vbox_for_buttons = gtk.VBox(False, 3)
           fragment_button = gtk.Button("  Sequence closest fragment  ")
	   chain_check_button = gtk.CheckButton("Assign Chain ID as well?")
           
           frame.add(hbox)
           vbox.pack_start(chain_id_label, False, False, 2)
           vbox.pack_start(entry, False, False, 2)
           vbox.pack_start(sequence_label, False, False, 2)
           vbox.pack_start(textview, True, False, 2)
           add_text_to_text_buffer(text_buffer, seq_info[1])
           entry.set_text(seq_info[0])
           hbox.pack_start(vbox, False, False, 2)
	   vbox_for_buttons.pack_start(chain_check_button, False, False, 6)
	   vbox_for_buttons.pack_start(fragment_button, False, False, 6)
           hbox.pack_start(vbox_for_buttons, False, False, 2)
              
           fragment_button.connect("clicked", fragment_go_event)
	   chain_check_button.connect("toggled", chain_toggled)
              
           return [frame, entry, text_buffer]

        # redraw the table when spin_button is changed
        def spin_button_changed(widget):
		global imported_sequence_file_flags
		no_of_frames = int(spin_button.get_value())
		imported_sequence_file_qm = imported_sequence_file_flags[0]
		no_of_sequences = imported_sequence_file_flags[1]

		# get table information
		table_children = seq_table.get_children()
		no_of_children = len(table_children)

		# make range for redraw
		if (imported_sequence_file_qm):
			# we only make extra ones
			redraw_range = range(no_of_sequences, no_of_frames)
			if (no_of_sequences > no_of_frames):
				child_range = []    # do not remove children
				redraw_range = []   # do not redraw
				spin_button.set_value(no_of_sequences)
				no_of_frames = no_of_sequences
			else:
				if (no_of_children == no_of_sequences):
					# do not remove any
					child_range = []
				else:
					child_range = range(0, no_of_children - no_of_sequences)
					# children seem to be added in the beginning ???????
		else:
			# we make everything new
			redraw_range = range(no_of_frames)
			child_range = range(no_of_children)
			
		# lets remove the children
		for j in child_range:
			child = table_children[j]
			seq_table.remove(child)
			
		# make new cells
		for i in redraw_range:
			make_cell(i)

		# resize the table
		seq_table.resize(no_of_frames, 1)

	# reset the table
	def clear_function_event(*args):
		global imported_sequence_file_flags
		imported_sequence_file_flags = [False, 0]
		seq_table.resize(1, 1)
		# remove children
		table_children = seq_table.get_children()
		for child in table_children:
			seq_table.remove(child)
		# make new
		make_cell(0)
		spin_button.set_value(1)
		
	# make one cell in line with deafult fill
	def make_cell(line):
		seq_widget = entry_text_pair_frame_with_button(["", "Cut and Paste Sequence to here or import a sequence file"])
		seq_table.attach(seq_widget[0], 0, 1, line, line + 1)
		seq_widget[0].show_all()
		
	# assign the in table given sequences to the model imol
	def assign_sequences_to_mol(imol):
		no_of_seq = int(spin_button.get_value())
		frame_list = seq_table.get_children()
		seq_all = []
		write_sequence = True
		for i in range(no_of_seq):
			fframe = frame_list[i]
			hhbox = fframe.get_children()[0]
			vvbox = hhbox.get_children()[0]
			child_list = vvbox.get_children()
			chain_id = child_list[1].get_text()
			seq_textbuffer = child_list[3].get_buffer()
			startiter, enditer = seq_textbuffer.get_bounds() 
			seq_in = seq_textbuffer.get_text(startiter, enditer)
			pair = [chain_id, seq_in]
			if (len(chain_id) == 0) or (
			    "Cut and Paste Sequence to here or import a sequence file") in pair:
				print "ERROR:: the input contains an invalid chain and/or sequence"
				write_sequence = False
			seq_all.append(pair)

		if (write_sequence):
			for element in seq_all:
				chain_id_new = element[0]
				seq = element[1].upper()
				# first check if chain_id is already in mol
				# if so delete it so that it can be replaced by the new sequence
				seq_info = sequence_info(imol)
				for info in seq_info:
					chain_id_old = info[0]
					if (chain_id_new == chain_id_old):
						delete_sequence_by_chain_id(imol, chain_id_old)
				
				# replace space, eol etc in sequence first
				seq = seq.replace(" ", "")
				seq = seq.replace("\r", "")       # mac?
				seq = seq.replace("\r\n", "")     # win?
				seq = seq.replace("\n", "")       # unix?
				assign_sequence_from_string(imol, chain_id_new, seq)
			return True
		else:
			add_status_bar_text("Invalid chain_id and/or sequence provided")
			return False


	# main body
	imol_map = imol_refinement_map()
	if (imol_map == -1):
		show_select_map_dialog()
	
	window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        window.set_title("Sequencing GUI")
        tooltips = gtk.Tooltips()
        label = gtk.Label("Molecule to be sequenced")
	vbox = gtk.VBox(False, 2)
        option_menu = gtk.combo_box_new_text()
        model_mol_list = fill_option_menu_with_mol_options(option_menu, valid_model_molecule_qm)
	inside_vbox = gtk.VBox(False, 2)
	seq_table = gtk.Table(1, 1, True)
        hbox_for_spin = gtk.HBox(False, 0)
        spin_label = gtk.Label("Number of Sequences:")
        spin_adj = gtk.Adjustment(1, 1, 10, 1, 4, 0)
        spin_button = gtk.SpinButton(spin_adj, 0, 0)
	refine_check_button = gtk.CheckButton("Auto-fit-refine after sequencing?")
	h_sep = gtk.HSeparator()
	h_sep2 = gtk.HSeparator()
	buttons_hbox = gtk.HBox(False, 2)
        import_button = gtk.Button("  Import and associate sequence from file  ")
	go_button = gtk.Button("  Sequence all fragments!  ")
        tooltips.set_tip(go_button, "This currently ignores all chain IDs")
	cancel_button = gtk.Button("  Cancel  ")
	clear_button = gtk.Button("  Clear all  ")

        window.set_default_size(400, 200)
        window.add(vbox)
        vbox.pack_start(label, False, False, 5)
        vbox.pack_start(option_menu, False, False, 0)

        hbox_for_spin.pack_start(spin_label, False, False, 2)
        hbox_for_spin.pack_start(spin_button, False, False, 2)
	hbox_for_spin.pack_end(refine_check_button, False, False, 2)
        vbox.pack_start(hbox_for_spin, False, False, 5)

        vbox.pack_start(inside_vbox, False, False, 2)
        inside_vbox.add(seq_table)
	make_cell(0)
	fill_table_with_sequences()
	
        vbox.pack_start(h_sep, False, False, 2)
        file_sel_entry = file_chooser_entry(vbox, "Select PIR file")
        vbox.pack_start(import_button, False, False, 6)

	buttons_hbox.pack_start(go_button, False, False, 6)
	buttons_hbox.pack_start(cancel_button, False, False, 6)
	buttons_hbox.pack_start(clear_button, False, False, 6)

        vbox.pack_start(h_sep2, False, False, 2)
        vbox.pack_start(buttons_hbox, False, False, 5)


        import_button.connect("clicked", import_function_event, file_sel_entry)

	cancel_button.connect("clicked", delete_event)

	go_button.connect("clicked", go_function_event)

	clear_button.connect("clicked", clear_function_event)
                               
	spin_adj.connect("value_changed", spin_button_changed)

	refine_check_button.connect("toggled", refine_function_event)

	option_menu.connect("changed", fill_table_with_sequences)

#        window.add(vbox)
	window.show_all()


# let the c++ part of mapview know that this file was loaded:
set_found_coot_python_gui()

