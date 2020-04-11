# coot-gui.py
#
# Copyright 2001 by Neil Jerram
# Copyright 2002-2008, 2009 by The University of York
# Copyright 2007 by Paul Emsley
# Copyright 2006, 2007, 2008, 2009 by Bernhard Lohkamp
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
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA


# Fire up the coot scripting gui.  This function is called from the
# main C++ code of coot.  Not much use if you don't have a gui to
# type functions in to start with.
#

import pygtk, gtk, pango
try:
   import gobject
except:
   print("WARNING:: no gobject available")

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
def coot_gui(own_gtk_main=False):
   
   global coot_listener_socket
   import sys, string
   import re

   if (coot_listener_socket):
      own_gtk_main = True
   
   def delete_event(*args):
       window.destroy()
       if (own_gtk_main):
          gtk.main_quit()
       return False

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

   # BL says:: now we check if entry was guile command (by mistake?!)
   # if so, let's translate it
   #
   def test_and_translate_guile(py_func):
       # test for '(' at start
       # and if a '-' before the '(' [i.e. if we have a guile command instead of a python command in the beginning,
       # we need to allow '-' after the '(']
       if ((string.find(py_func[0:py_func.find('(')],"-") > 0) or (string.find(py_func,"(") == 0)):
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
                       print("BL WARNING:: unknown type of argument!")

          python_function = tmp_command + "(*arg)"
#          print "BL DEBUG:: python func and args are", python_function, tmp_command, arg

          try:
             res = eval(python_function)
          except SyntaxError:
             exec(python_function, locals())
             res = None
             return True
          except:
             return False
          if res is not None:
             # silence is golden
             # print "BL DEBUG:: result is", res
             insert_normal_text(str(res) + "\n")
             return True
       else:
          # silence is golden
          # print "This is not a guile command!"
          return False

   def do_function(widget, entry):
       global histpos
       entry_text = entry.get_text()
       # enough of this debugging output
       # print "BL INFO:: command input is: ", entry_text
       if (entry_text != None):
          insert_tag_text(textbuffer.create_tag(foreground="darkblue"),
                          entry_text + "\n")
       while gtk.events_pending():
         gtk.main_iteration(False)
       # save to history by default now!? only in severe cases we dont
       his = True
       res = None
       try:
             res = eval(entry_text)
             his = True
       except SyntaxError:
         try:
             test_equal = entry_text.split("=")
             if ((len(test_equal) > 1)
                 and (test_equal[0].strip() in globals())
                 and (callable(eval(test_equal[0].strip())))):
                # we try to assign a value to a build-in function! Shouldnt!
                warning_text = "BL WARNING:: you tried to assign a value to the\nbuild-in function: " + test_equal[0].strip() + "!\nThis is not allowed!!\n"
                insert_normal_text(warning_text)
                res = False
                his = False
             else:
                exec(entry_text, globals())
                res = None
                his = True
         except:
          guile_function = test_and_translate_guile(entry_text)
          if guile_function:
             print("BL INFO::  We have a guile command!")
             insert_normal_text("BL INFO:: Detected guile scripting!\n\
             You should use python commands!!\n\
             But I'm a nice guy and translated it for you, this time...!\n")
          else:
             # I don't need to see this - and its often not right (not a syntax error)
             # insert_normal_text("BL WARNING:: Python syntax error!\n\
             # (Or you attempted to use an invalid guile command...)\n")
             type_error, error_value = sys.exc_info()[:2]
             error = str(error_value)
             insert_normal_text("Python error:\n") 
             insert_normal_text(error + "\n")
             insert_normal_text(str(type_error) + "\n")

       except:
          if (entry_text != None):
             guile_function =  test_and_translate_guile(entry_text)
             if guile_function:
                insert_normal_text("BL INFO:: Detected guile scripting!\nYou should use python commands!!\nBut I'm a nice guy and translated it for you, this time...!\n")
             else:
                insert_normal_text("BL WARNING:: Python error!\n(Or you attempted to use an invalid guile command...)\n")
                type_error, error_value = sys.exc_info()[:2]
                error = str(error_value)
                insert_normal_text("Python error:\n") 
                insert_normal_text(error + "\n")
                insert_normal_text(str(type_error) + "\n")
          else:
             print("BL WARNING:: None input")

       if res is not None:
             print("BL INFO:: result is", res)
             insert_normal_text(str(res) + "\n")

       if his:
             l = entry_text + '\n'
             if len(l) > 1:
                # only save if there is something to save
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
       # BL for debug on mingw
       global have_mingw
       if have_mingw:
          sys.stderr.flush()
          sys.stdout.flush()
       

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
      tmp = [i][0]
      if (not tmp[-3:len(tmp)] == '_py'):
          liststore.append([i])
   liststore.set_sort_column_id(0, gtk.SORT_ASCENDING)
   completion.set_model(liststore)
   completion.set_text_column(0)
#   completion.connect("match-selected", match_cb)
   close_button = gtk.Button("  Close  ")
   vbox = gtk.VBox(False, 0)
   hbox = gtk.HBox(False, 0)
   label = gtk.Label("Command: ")
   ifont = gtk.gdk.Font("fixed")
   window.set_default_size(550,250)
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
   # run in own gtk loop only when using sockets!
   if (own_gtk_main):
      gtk.main()


# The callback from pressing the Go button in the smiles widget, an
# interface to run libcheck.
#
def handle_smiles_go(tlc_entry, smiles_entry):

    tlc_text = tlc_entry.get_text()
    smiles_text = smiles_entry.get_text()
    use_libcheck = False
    if is_windows():
        use_libcheck = True
    new_molecule_by_smiles_string(tlc_text, smiles_text, force_libcheck=use_libcheck)

# smiles GUI
#
def smiles_gui():

   def smiles_gui_internal():
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
      if enhanced_ligand_coot_p():
          text = gtk.Label("  [SMILES interface works by using Pyrogen]  ")
      else:
          text = gtk.Label("  [SMILES interface works by using CCP4's LIBCHECK]  ")
      go_button = gtk.Button("  Go  ")
      vbox.pack_start(hbox1, False, False, 0)
      vbox.pack_start(hbox2, False, False, 4)
      vbox.pack_start(text, False, False, 2)
      vbox.pack_start(go_button, False, False, 6)
      hbox1.pack_start(tlc_label, False, False, 0)
      hbox1.pack_start(tlc_entry, False, False, 0)
      hbox2.pack_start(smiles_label, False, False, 0)
      hbox2.pack_start(smiles_entry, True, True, 4)
      smiles_window.add(vbox)
      vbox.set_border_width(6)

      smiles_entry.connect("key-press-event", smiles_connect, tlc_entry, smiles_entry, smiles_window)
      go_button.connect("clicked", go_button_pressed,tlc_entry,smiles_entry,smiles_window)
      smiles_window.connect("delete_event",delete_event)

      smiles_window.show_all()
      
   # first check that libcheck is available... if not put up and info
   # dialog.
   if enhanced_ligand_coot_p():
      smiles_gui_internal()
   else:
      if (find_exe(libcheck_exe, "CCP4_BIN", "PATH")):
         smiles_gui_internal()
      else:
         info_dialog("You need to setup CCP4 (specifically LIBCHECK) first.")
      

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
    window.set_title('Coot')
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
    hbox1.pack_start(function_label, False, False, 0)
    hbox2.pack_start(smiles_entry, True, True, 0)
    window.add(vbox)
    vbox.set_border_width(6)
 
    if isinstance(entry_1_default_text,(str,)):
       smiles_entry.set_text(entry_1_default_text)
    else:
       print("BL WARNING:: entry_1_default_text was no string!!")
   
    cancel_button.connect("clicked", delete_event)

    go_button.connect("clicked", go_function_event, smiles_entry)

    smiles_entry.connect("key_press_event", key_press_event, smiles_entry)

    window.show_all()

# generic double entry widget, now with a check button
# ...and returns the widget if requested
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
def generic_double_entry(label_1, label_2,
                         entry_1_default_text, entry_2_default_text,
                         check_button_label, handle_check_button_function,
                         go_button_label, handle_go_function,
                         return_widget = False):

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
    window.set_title('Coot')
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
    hbox3.pack_start(go_button, False, False, 6)
    hbox3.pack_start(cancel_button, True, False, 6)
    hbox1.pack_start(tlc_label, False, False, 0)
    hbox1.pack_start(tlc_entry, False, False, 0)
    hbox2.pack_start(smiles_label, False, False, 0)
    hbox2.pack_start(smiles_entry, True, True, 0)

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

    if isinstance(entry_1_default_text,(str,)):
       tlc_entry.set_text(entry_1_default_text)
    else:
       print("BL WARNING:: entry_1_default_text was no string!!")
 
    if isinstance(entry_2_default_text,(str,)):
       smiles_entry.set_text(entry_2_default_text)
    else:
       print("BL WARNING:: entry_2_default_text was no string!!")

    go_button.connect("clicked", go_function_event)
    cancel_button.connect("clicked", delete_event)

    smiles_entry.connect("key-press-event", key_press_event, tlc_entry, smiles_entry, check_button)

    window.set_default_size(400, 100)
    window.show_all()

    # return the widget
    if (return_widget):
       return window

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
       if check_button:
          handle_go_function([entry.get_text() for entry in entries], check_button.get_active())
       else:
          handle_go_function([entry.get_text() for entry in entries])
       delete_event()

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    window.set_title('Coot')
    vbox = gtk.VBox(False, 0)
    hbox3 = gtk.HBox(False, 0)
    h_sep = gtk.HSeparator()
    cancel_button = gtk.Button("  Cancel  ")
    go_button = gtk.Button(go_button_label)

    # all the labelled entries
    #
    entries = []
    for entry_info in entry_info_list:

       entry_1_hint_text = entry_info[0]
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

    # print "debug:: check-button-info: ", check_button_info
    # print "debug:: entry-info-list: ", entry_info_list

    if not (type(check_button_info) is ListType):
            # and len(check_button_info) == 2): we can't do len on a bool
       check_button = False
    else:
       if type(check_button_info[0]) is StringType:

          def check_callback(*args):
             active_state = c_button.get_active()
             check_button_info[1] = active_state

          c_button = gtk.CheckButton(check_button_info[0])
          vbox.pack_start(c_button, False, False, 2)
          c_button.connect("toggled", check_callback)
          check_button = c_button
       else:
          check_button = False      # the check-button when we don't want to see it

    # print "Here check button creation.................. check-button is ", check_button
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
           print("BL WARNING:: libcheck didn't seem to run ok! Please check output carefully!!")
        else: pass
    fp.close()

# old coot test
def old_coot_qm():

   import time, random
   # when making a new date, recall that guile times seem to be
   # differently formatted to python times (god knows why)
   # run e.g.: time.mktime((2008, 7, 24, 0,0,0,0,0,0))
   # 
   # new_release_time = 1199919600 # 10 Jan 2008
   # new_release_time = 1205622000 # 16 Mar 2008 0.3.3
   # new_release_time = 1216854000 # 24 Jul 2008 0.4
   # new_release_time = 1237244400 # 17 March 2009
   # new_release_time = 1249945200 # 11 August 2009 : 0.5
   # new_release_time = 1279926000 # 24 July 2010 : --
   # new_release_time = 1310000000 #  7 July 2011 : 0.6
   # new_release_time = 1330000000 # 23 Jan 2012
   new_release_time = 1352674800 # 12 Nov 2012
   current_time = int(time.time())
   time_diff = current_time - new_release_time
   if (time_diff > 0):
      if (time_diff > 8600):
         s = "You've got an Old Coot!\n\nIt's time to upgrade."
      else:
         if (random.randint(0,10) == 0):
            # Jorge Garcia:
            s = "(Nothing says \"patriotism\" like an Ireland shirt...)\n"
         else:
            s = "You've got an Old Coot!\n\nIt's time to upgrade."
      info_dialog(s)

#if (not coot_has_guile()):
#   old_coot_qm()

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

          if (isinstance(func_maybe1, list) and len(func_maybe1)>0):
             # the last one is probably a funcn (no button name)
             func_maybe_strip = func_maybe1[0]
#             print "BL DEBUG:: func_maybe_strip is", func_maybe_strip
             if (callable(func_maybe_strip)):
                return func_maybe1, False, False
             else:
                return False, False, False
             
          elif (isinstance(func_maybe2, list) and len(func_maybe2)>0
                and isinstance(func_maybe1, bytes)):
             # the second last is function, last is button name
             func_maybe_strip = func_maybe2[0]
             button_name = func_maybe1
             if (callable(func_maybe_strip)):
                return func_maybe2, button_name, False
             else:
                return False, False, False
             
          elif (isinstance(func_maybe3, list) and len(func_maybe3)>0
                and isinstance(func_maybe2, bytes)
                and isinstance(func_maybe1, bytes)):
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
       print("Attempt to go to chain: %s resno: %s atom-name: %s" %(atom_info[0],atom_info[1],atom_info[2]))
       set_go_to_atom_molecule(mol_no)
       success = set_go_to_atom_chain_residue_atom_name(*atom_info)
       if success == 0:           # failed?!
          new_name = unmangle_hydrogen_name(atom_info[2])
          success2 = set_go_to_atom_chain_residue_atom_name(atom_info[0],atom_info[1],new_name)
          if success2 == 0:
             print("Failed to centre on demangled name: ", new_name)
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
   #tooltips = gtk.Tooltips()

   for baddie_items in baddie_list:

       if baddie_items == '' :
          print('Done')
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
             if gtk.pygtk_version >= (2,12):
                button.set_tooltip_text(tooltip_str)
             else:
                coot_tooltips.set_tip(button, tooltip_str)


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

    mol_ls = []
    n_molecules = graphics_n_molecules()
   
    for mol_no in molecule_number_list():
        if filter_function(mol_no):
           label_str = molecule_name(mol_no)
           if (isinstance(label_str,(str,))):
              mlabel_str = str(mol_no) + " " + label_str
              menu.append_text(mlabel_str)
              menu.set_active(0)
              mol_ls.append(mol_no)
           else:
              print("OOps molecule name for molecule %s is %s" %(mol_no_ls,label_str))
    return mol_ls

# Fill an option menu with maps and return the list of maps
#
def fill_option_menu_with_map_mol_options(menu):
    return fill_option_menu_with_mol_options(menu, valid_map_molecule_qm)

# Helper function for molecule chooser.  Not really for users.
# 
# Return a list of models, corresponding to the menu items of the
# option menu.
# 
# The returned list will not contain references to map or closed
# molecules.
# 
def fill_option_menu_with_coordinates_mol_options(menu):
    return fill_option_menu_with_mol_options(menu, valid_model_molecule_qm)

#
def fill_option_menu_with_number_options(menu, number_list, default_option_value):

    for number in number_list:
       mlabel_str = str(number)
       menu.append_text(mlabel_str)
       if (default_option_value == number):
          count = number_list.index(number)
          menu.set_active(count)
          print("setting menu active ", default_option_value, count)

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
    children = len(model)

    if (children == len(model_mol_list)):
       try:
          all_model = model[active_item][0]
          imol_model, junk = all_model.split(' ', 1)
       
          return int(imol_model)
       except:
          print("INFO:: could not get active_item")
          return False
    else:
       print("Failed children length test : ",children, model_mol_list)
       return False

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

   if (children == len(item_list)):
       try:
          all_model = model[active_item][0]
          return all_model
       except:
          print("INFO:: could not get active_item")
          return False
   else:
       print("Failed children length test : ",children, item_list)
       return False

# Typically option_menu_fill_function is
# fill_option_menu_with_coordinates_mol_options
#
def molecule_chooser_gui_generic(chooser_label, callback_function, option_menu_fill_function):
 
    def delete_event(*args):
       window.destroy()
       return False

    def on_ok_clicked(*args):
        # what is the molecule number of the option menu?
        active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)
        try:
           active_mol_no = int(active_mol_no)
           print("INFO: operating on molecule number ", active_mol_no)
           try:
              callback_function(active_mol_no)
           except:
              print("BL INFO:: problem in callback_function", callback_function.__name__)
           delete_event()
        except:
           print("Failed to get a (molecule) number")

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    window.set_title('Coot')
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

# Fire up a coordinates/model molecule chooser dialog, with a given
# label and on OK we call the callback_fuction with an argument of
# the chosen molecule number. 
# 
# chooser_label is a directive to the user such as "Choose a Molecule"
# 
# callback_function is a function that takes a molecule number as an
# argument.
#
def molecule_chooser_gui(chooser_label, callback_function):
    molecule_chooser_gui_generic(chooser_label,
                                 callback_function,
                                 fill_option_menu_with_coordinates_mol_options)

# Fire up a map molecule chooser dialog, with a given label and on OK we
# call the callback_fuction with an argument of the chosen molecule
# number. 
# 
# chooser_label is a directive to the user such as "Choose a Molecule"
# 
# callback_function is a function that takes a molecule number as an
# argument.
#
def map_molecule_chooser_gui(chooser_label, callback_function):
    molecule_chooser_gui_generic(chooser_label,
                                 callback_function,
                                 fill_option_menu_with_map_mol_options)

# A pair of widgets, a molecule chooser and an entry.  The
# callback_function is a function that takes a molecule number and a
# text string.
#
def generic_chooser_and_entry(chooser_label, entry_hint_text,
                              default_entry_text, callback_function,
                              always_dismiss_on_ok_clicked=True):

   print("BL DEBUG:: --- deal with always_dissmiss...", always_dismiss_on_ok_clicked)
   # cf = lambda text, dummy: callback_function(text)
   generic_chooser_and_entry_and_check_button(chooser_label, entry_hint_text,
                                              default_entry_text, False,
                                              callback_function,
                                              always_dismiss_on_ok_clicked)

# as above , plus we also have a check-button
# [and an additional argument in the callback - actually not I think ]
# If check-button-label is false, then don't create a check-button.
def generic_chooser_and_entry_and_check_button(chooser_label, entry_hint_text,
                                              default_entry_text, check_button_label,
                                              callback_function,
                                              always_dismiss_on_ok_clicked=True):

   import operator

   def delete_event(*args):
      window.destroy()
      return False

   def on_ok_button_clicked(*args):
      # what is the molecule number of the option menu?
      active_mol_no = get_option_menu_active_molecule(option_menu,model_mol_list)

      try:
         active_mol_no = int(active_mol_no)
         print("INFO: operating on molecule number ", active_mol_no)
         text = entry.get_text()
         if check_button:
            check_button_state = check_button.get_active()
            cbf_ret = callback_function(active_mol_no, text, check_button_state)
         else:
            cbf_ret = callback_function(active_mol_no, text)
         if always_dismiss_on_ok_clicked:
            delete_event()
         else:
            if cbf_ret:
               delete_event()
            else:
               return True
      except:
         print("Failed to get a (molecule) number")

   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   window.set_title('Coot')
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
   if check_button_label:
      check_button = gtk.CheckButton(check_button_label)
      vbox.pack_start(check_button, False, False, 2)
   else:
      check_button = False
   vbox.pack_start(h_sep, True, False, 2)
   vbox.pack_start(hbox_buttons, False, False, 5)
   hbox_buttons.pack_start(ok_button, True, False, 5)
   hbox_buttons.pack_start(cancel_button, False, False, 5)
   hbox_for_entry.pack_start(entry_label, False, False, 4)
   hbox_for_entry.pack_start(entry, True, True, 4)
   entry.set_text(default_entry_text)

   # button callbacks
   ok_button.connect("clicked", on_ok_button_clicked, entry,
                     option_menu, callback_function, check_button)
   cancel_button.connect("clicked", delete_event)

   window.show_all()

# Create a window
#
# Return a pair of widgets, a chooser entry and a file selector.  The
# callback_function is a function that takes a molecule number and a
# text string (e.g. chain_id and file_name)
#
# chooser_filter is typically valid_map_molecule_qm or valid_model_molecule_qm
#
# Add option to add a checkbutton
#
def generic_chooser_entry_and_file_selector(chooser_label, chooser_filter,
                                            entry_hint_text, default_entry_text,
                                            file_selector_hint, callback_function,
                                            use_check_button=False,
                                            check_button_label="",
                                            alternative_callback_function=None):

    import operator

    def delete_event(*args):
       window.destroy()
       return False

    def on_ok_button_clicked(*args):
        # what is the molecule number of the option menu?
        active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)
        
        try:
           active_mol_no = int(active_mol_no)
           text = entry.get_text()
           file_sel_text = file_sel_entry.get_text()
           if (c_button and c_button.get_active()):
              # use alt function
              alternative_callback_function(active_mol_no,
                                            text, file_sel_text)
           else:
              callback_function(active_mol_no, text, file_sel_text)
           delete_event()
        except:
           print("Failed to get a (molecule) number")

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

    c_button = None
    if use_check_button:
       # now add a check button
       c_button = gtk.CheckButton(check_button_label)
       vbox.pack_start(c_button, False, False, 2)

    file_sel_entry = file_selector_entry(vbox, file_selector_hint)
    vbox.pack_start(h_sep, True, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 5)
    

    # button callbacks
    ok_button.connect("clicked", on_ok_button_clicked, entry, option_menu,
                      callback_function,
                      c_button, alternative_callback_function)
    cancel_button.connect("clicked", delete_event)

    window.show_all()

# Create a window.
#
# Return a pair of widgets, a molecule chooser and a file selector.  The
# callback_function is a function that takes a molecule number and a
# file_name
#
# chooser_filter is typically valid_map_molecule_qm or valid_model_molecule_qm
#
def generic_chooser_and_file_selector(chooser_label,
                                      chooser_filter,
                                      file_selector_hint,
                                      default_file_name,
                                      callback_function):

    import operator

    def delete_event(*args):
       window.destroy()
       return False

    def on_ok_button_clicked(*args):
        # what is the molecule number of the option menu?
        active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)

        try:
           active_mol_no = int(active_mol_no)
           file_sel_text = file_sel_entry.get_text()
           callback_function(active_mol_no, file_sel_text)
           delete_event()
        except:
           print("Failed to get a (molecule) number")

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

    file_sel_entry = file_selector_entry(vbox,
                                         file_selector_hint,
                                         default_file_name)
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
      for menu_child in coot_main_menubar.get_children():
          lab = []
          # lab is a GtkAccelLabel list
          lab.append(menu_child.get_children()[0].get_text())
          lab.append(menu_child)
          ac_lab_ls.append(lab)
      return ac_lab_ls

    # main body
    #
    found_menu = False
    for f in menu_bar_label_list():
       if menu_label.lstrip("_") in f: 
          #print "BL DEBUG:: found menu label is ", f
          # we shall return the submenu and not the menuitem
          found_menu = f[1].get_submenu()
    if found_menu:
       return found_menu
    else:
       menu = gtk.Menu()
       menuitem = gtk.MenuItem(menu_label)
       menuitem.set_submenu(menu)
       coot_main_menubar.append(menuitem)
       menuitem.show()
       return menu
   except: print("""BL WARNING:: could not import coot_python module!!\n
                    Some things, esp. extensions, may be crippled!""")


# Given that we have a menu (e.g. one called "Extensions") provide a
# cleaner interface to adding something to it:
# 
# activate_function is a thunk.
#
def add_simple_coot_menu_menuitem(menu, menu_item_label, activate_function):

    submenu = gtk.Menu()
    sub_menuitem = gtk.MenuItem(menu_item_label)

    menu.append(sub_menuitem)
    sub_menuitem.show()

    sub_menuitem.connect("activate", activate_function)


# Make an interesting things GUI for residues of molecule number
# imol that have alternate conformations.
#
def alt_confs_gui(imol):

   interesting_residues_gui(imol, "Residues with Alt Confs",
                            residues_with_alt_confs(imol))

# Make an interesting things GUI for residues with missing atoms
#
def missing_atoms_gui(imol):

   interesting_residues_gui(imol, "Residues with missing atoms",
                            missing_atom_info(imol))

# Make an interesting things GUI for residues with zero occupancy atoms
#
def zero_occ_atoms_gui(imol):

   atom_ls = atoms_with_zero_occ(imol)
   if atom_ls:
      interesting_things_gui("Residues with zero occupancy atoms",
                             atom_ls)
   else:
      s  = "No atoms with zero occupancy found\n"
      s += "in molecule "
      s += str(imol)
      info_dialog(s)


# Make an interesting things GUI for residues of molecule number
# imol for the given imol.   A generalization of alt-confs gui.
#
def interesting_residues_gui(imol, title, interesting_residues):

   from types import ListType
   centre_atoms = []
   if valid_model_molecule_qm(imol):
      residues = interesting_residues
      for spec in residues:
         if (type(spec) is ListType):
            centre_atoms.append(residue_spec_to_atom_for_centre(imol, *spec))
         else:
            centre_atoms.append([False])

      interesting_things_gui(
         title,
         list(map(lambda residue_cpmd, centre_atom: 
             [residue_cpmd[0] + " " +
              str(residue_cpmd[1]) + " " +
              residue_cpmd[2] + " " +
              centre_atom[0] + " " + centre_atom[1],
              imol, residue_cpmd[0], residue_cpmd[1], residue_cpmd[2],
              centre_atom[0], centre_atom[1]] if centre_atom else
              ["[oops - why did this happen?]", 0, 0, 0, 0, 0, 0],
             residues, centre_atoms)))
      
   else:
      print("BL WARNING:: no valid molecule", imol)

                
# If a toolbutton with label button_label is not found in the coot main
# toolbar, then create it and return it.
# If it does exist, the icon will be overwritten. The callback function wont!
# [NOTE: if its a Coot internal Toolbar you can
# only add a function but not overwrite its internal function!!!]
#
# return False if we cannot create the button and/or wrong no of arguments
#
# we accept 6 arguments:
#   button_label
#   callback_function   (can be a string or callable function, or list
#                        with first element fucntion, then args, i.e.
#                        [function, arg1, arg2, ...])
#   gtk-stock-item      (or icon_widget, whatever that is)
#   tooltip             a text to present as tooltip
#   toggle_button_flag  toggle button instead of normal button
#   use_button          pass the button widget onto the function
#  
def coot_toolbar_button(button_label, cb_function,
                        icon_name=False, tooltip=False,
                        toggle_button=False, use_button=False):

   import types
   # we do not exclusively use strings any more...
   #print "BL DEBUG:: create toolbutton with args!", args
   #if (not type(cb_function) is StringType):
   #   print "BL WARNING:: callback function wasn't a string! Cannot create toolbarbutton!"
   #   return False
   
   try:
      import coot_python
   except:
      print("""BL WARNING:: could not import coot_python module!!
      So we cannot make toolbar_buttons!""")
      return False
   
   coot_main_toolbar = coot_python.main_toolbar()

   # main body
   #
   found_button = False
   for f in toolbar_label_list():
      if button_label in f: 
         found_button = f[1]
   if found_button:
      # here we only try to add a new icon, we cannot overwrite the callback function!
      toolbutton = found_button
   else:
      if toggle_button:
         toolbutton = gtk.ToggleToolButton()
         toolbutton.set_label(button_label)
      else:
         toolbutton = gtk.ToolButton(icon_widget=None,
                                     label=button_label)
      coot_main_toolbar.insert(toolbutton, -1)       # insert at the end
      toolbutton.set_is_important(True)              # to display the text,
                                                     # otherwise only icon
      # tooltips?
      if tooltip:
         if gtk.pygtk_version >= (2,12):
            toolbutton.set_tooltip_text(tooltip)
         else:
            coot_tooltips.set_tip(toolbutton, tooltip)

      def cb_wrapper(widget, callback_function):
         if (type(callback_function) is bytes):
            # old style with string as function
            eval(callback_function)
         else:
            # have function as callable and maybe extra args (all in one list)
            args = []
            function = callback_function
            if (type(callback_function) is list):
               function = callback_function[0]
               args = callback_function[1:]
            # pass the widget/button as well? Maybe the cb function can
            # make use of it
            if use_button:
               args.append(widget)
            if callable(function):
               function(*args)
            else:
               print("BL ERROR:: cannot evaluate or call function", function)
      #toolbutton.connect("clicked", lambda w: eval(cb_function))
      toolbutton.connect("clicked", cb_wrapper, cb_function)
      toolbutton.show()

   if (icon_name):
      # try to add a stock item
      try:
         toolbutton.set_stock_id(icon_name)
      except:
         # try to add a icon widget
         try:
            toolbutton.set_icon_widget(icon_name)
         except:
            print("BL INFO::  icon name/widget given but could not add the icon")

   return toolbutton

# If a combobox with name is not found in the coot main
# toolbar, then create it and return it. If it is there give a warning
#
# return False if we cannot create the button and/or wrong no of arguments
#
# we accept 3 arguments:
#   label               name of the tool item (to keep track)
#   entry_list          labels for all entries of the box
#   callback_function   (can be a string or callable function, or list
#                        with first element fucntion, then args, i.e.
#                        [function, arg1, arg2, ...])
#   tooltip             a text to present as tooltip
#  
def coot_toolbar_combobox(label, entry_list, cb_function, tooltip=""):

   try:
      import coot_python
   except:
      print("""BL WARNING:: could not import coot_python module!!
      So we cannot make a toolbar combobox!""")
      return False
   
   coot_main_toolbar = coot_python.main_toolbar()

   # main body
   #
   found_combobox = False
   for f in toolbar_label_list():
      if label in f: 
         found_combobox = f[1]
   if found_combobox:
      s = "BL WARNING:: already have a comboxbox with name %s. Please " +\
          "try again!" %found_combobox
      info_dialog(s)
      return False
   else:
      toolitem = gtk.ToolItem()
      toolitem.set_name(label)
      combobox = gtk.combo_box_new_text()
      for text in entry_list:
         combobox.append_text(text)
      combobox.set_active(0)

      # tooltips?
      if tooltip:
         if gtk.pygtk_version >= (2,12):
            combobox.set_tooltip_text(tooltip)
         else:
            coot_tooltips.set_tip(combobox, tooltip)

      def cb_wrapper(widget, callback_function):
         pos = widget.get_active()
         if (type(callback_function) is bytes):
            # old style with string as function
            # does not take into account the positions pos! Or how?
            # FIXME? Maybe
            eval(callback_function)
         else:
            # have function as callable and maybe extra args (all in one list)
            args = []
            function = callback_function
            # if there is a list, then we add the pos as last arg
            if (isinstance(callback_function, list)):
               function = callback_function[0]
               args = callback_function[1:]
               args.append(pos)
            # pass the widget/button as well? Maybe the cb function can
            # make use of it
            if callable(function):
               function(*args)
            else:
               print("BL ERROR:: cannot evaluate or call function", function)

      combobox.connect("changed", cb_wrapper, cb_function)
      toolitem.add(combobox)
      coot_main_toolbar.insert(toolitem, -1)
      toolitem.show_all()

   return toolitem
   

# returns a list of existing toolbar buttons
# [[label, toolbutton],[]...]
# or False if coot_python is not available
#
def toolbar_label_list():
   try:
      import coot_python
   except:
      return False
   
   coot_main_toolbar = coot_python.main_toolbar()
   button_label_ls = []
   for toolbar_child in coot_main_toolbar.get_children():
      ls = []
      try:
         label = toolbar_child.get_label()
         ls.append(label)
         ls.append(toolbar_child)
         button_label_ls.append(ls)
      except:
        # some toolitems we have from GTK2 cannot be accessed here
        # so we pass it and dont add it to the list.
        # nothing we can do about it. Probably good as we dont want to
        # play with the GTK2 toolitems only with the PYGTK ones.
        #
        # try if it has a name instead, e.g. toolitem
        try:
           label = toolbar_child.get_name()
           ls.append(label)
           ls.append(toolbar_child)
           button_label_ls.append(ls)
        except:
           pass

   return button_label_ls


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
		centre_atoms = list(map(residue_spec, interesting_residues))
		if centre_atoms:
			# BL says:: ignoring "Atom in residue name failure" for nor
			interesting_things_gui(gui_title_string,
			list(map(lambda interesting_residue, centre_atom:
			[interesting_residue[0] + " " 
			+ interesting_residue[1] + " " 
			+ str(interesting_residue[2]) + " " 
			+ interesting_residue[3] + " " 
			+ centre_atom[0] + " " + centre_atom[1]]
			+ interesting_residue[1:len(interesting_residue)] + centre_atom,
			interesting_residues, centre_atoms))
			)
	else:
		print("BL WARNING:: no valid model molecule ", imol)

# A gui that makes a generic number chooser the go function is a
# function that takes the value of the active menu item - as a
# number.
#
def generic_number_chooser(number_list, default_option_value, hint_text,
                           go_button_label, go_function):

    def delete_event(*args):
       window.destroy()
       return False

    def go_button_pressed(*args):
        active_number = int(get_option_menu_active_item(option_menu, number_list))
        try:
#           print "BL DEBUG:: go_function is:", go_function
#           print "BL DEBUG:: active_number is:", active_number
           go_function(active_number)
        except:
           print("Failed to get execute function")
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
def entry_do_button(vbox, hint_text,
                    button_label, button_press_func,
                    entry_text=False):

   hbox = gtk.HBox(False, 0)
   entry = gtk.Entry()
   button = gtk.Button(button_label)
   label = gtk.Label(hint_text)

   hbox.pack_start(label, False, False, 2)
   hbox.pack_start(entry, True, False, 2)
   hbox.pack_start(button, False, False, 2)
   button.connect("clicked", button_press_func, entry)

   if entry_text:
      entry.set_text(entry_text)
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
    return [option_menu, model_mol_list]

# Return an entry, the widget is inserted into the hbox passed to
# this function
#
def file_selector_entry(hbox, hint_text, default_file_name = False):

   if (file_chooser_selector_state() == 0 or gtk.pygtk_version < (2,3,90)):

        vbox = gtk.VBox(False, 0)

        def file_func1(*args):
                def file_ok_sel(*args):
                        t = fs_window.get_filename()
                        print(t)
                        entry.set_text(t)
                        fs_window.destroy()

                fs_window = gtk.FileSelection("file selection")
                fs_window.ok_button.connect("clicked", file_ok_sel)
                fs_window.cancel_button.connect("clicked",
                                lambda w: fs_window.destroy())
                fs_window.show()

        entry = entry_do_button(vbox, hint_text,
                            "  File...  ", file_func1,
                            default_file_name)

        hbox.pack_start(vbox, False, False, 2)
        vbox.show()
        return entry
     
   else:
      return file_chooser_entry(hbox, hint_text, default_file_name)

# This is the same as the file_selector_entry, but using the modern FileChooser
# Return an entry, the widget is inserted into the hbox passed to
# this function
#
def file_chooser_entry(hbox, hint_text, default_file_name = False):

   if gtk.pygtk_version > (2,3,90):

        vbox = gtk.VBox(False, 0)

        def file_func1(*args):
                def file_ok_sel(*args):
                        t = fc_window.get_filename()
                        print(t)
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

        entry = entry_do_button(vbox, hint_text, "  File...  ",
                            file_func1, default_file_name)

        hbox.pack_start(vbox, False, False, 2)
        vbox.show()
        return entry
   else:
        print("PyGtk 2.3.90 or later required for this function!")
        return False

# The gui for the strand placement function
#
def place_strand_here_gui():

    generic_number_chooser(number_list(4,12), 7, 
                            " Estimated number of residues in strand",
                            "  Go  ",
                            lambda n: place_strand_here(n, 15))

# Cootaneer gui
#
def cootaneer_gui(imol):

    def delete_event(*args):
       window.destroy()
       return False

    def go_function_event(widget, imol):
           print("apply the sequence info here\n")
           print("then cootaneer\n")

           # no active atom won't do.  We need
           # to find the nearest atom in imol to (rotation-centre).
           #
           # if it is too far away, give a
           # warning and do't do anything.

           n_atom = closest_atom(imol)
           if n_atom:
              imol	= n_atom[0]
              chain_id = n_atom[1]
              resno    = n_atom[2]
              inscode  = n_atom[3]
              at_name  = n_atom[4]
              alt_conf = n_atom[5]
              cootaneer_results = cootaneer(imol_map, imol, [chain_id, resno, inscode, 
                                                             at_name, alt_conf])
              print("Cootaneering status:", cootaneer_results)
              if (cootaneer_results == 0):
                 s = "Insufficiently confident in alignment to make a fit." + \
                     "\n" + \
                     "Perhaps you could improve or extend this fragment."
                 window.destroy()
                 info_dialog(s)
           else:
              print("BL WARNING:: no close atom found!")
              window.destroy()


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
           print("BL DEBUG:: probably should wait here for input!?")

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    outside_vbox = gtk.VBox(False, 2)
    inside_vbox = gtk.VBox(False, 2)
    h_sep = gtk.HSeparator()
    buttons_hbox = gtk.HBox(True, 2)
    go_button = gtk.Button("  Cootaneer!  ")
    cancel_button = gtk.Button("  Cancel  ")

    seq_info_ls = sequence_info(imol)
    # print "BL DEBUG:: sequence_list and imol is", seq_info_ls, imol

    if not seq_info_ls:
           s = "No sequence assigned for molecule number " + str(imol)
           print(s)
           info_dialog(s)
    else:

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
#
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

   def add_view_local_func(text):
      new_view_number = add_view_here(text)
      add_view_to_views_panel(text, new_view_number)
   generic_single_entry("View Name: ", local_view_name(), " Add View ",
                        lambda text: add_view_local_func(text))

def add_view_to_views_panel(view_name, view_number):
   global views_dialog_vbox
   if (views_dialog_vbox):
      button = gtk.Button(view_name)
      button.connect("clicked", lambda func: go_to_view_number(view_number, 0))
      views_dialog_vbox.pack_start(button, False, False, 2)
      button.show()

# return a list of [h_box_buttons, window]
#
# a button is a list of [label, callback, text_description]
#
def dialog_box_of_buttons(window_name, geometry, buttons,
                          close_button_label, post_hook=None):
   return dialog_box_of_buttons_with_check_button(window_name, geometry,
                                                  buttons, close_button_label,
                                                  False, False, False, post_hook)

# geometry is an improper list of ints
#
# return a list of [h_box_buttons, window]
#
# a button is a list of [label, callback, (optional: text_description)]
# where callback is a string or list of strings to be evaluated
#
# If check_button_label is False, don't make one, otherwise create with
# the given label and "on" state
#
# Add a post hook, to execute a function after the dialog is closed.
# Default is None.
#
# Note:
#  - if label is "HSep" a horizontal separator is inserted instead of a button
#  - the description is optional
#
def dialog_box_of_buttons_with_check_button(window_name, geometry,
                                            buttons, close_button_label,
                                            check_button_label,
                                            check_button_func,
                                            check_button_is_initially_on_flag,
                                            post_close_hook=None):
   

   def add_text_to_text_widget(text_box, description):
      textbuffer = text_box.get_buffer()
      start = textbuffer.get_start_iter()
      textbuffer.create_tag("tag", foreground="black", 
                            background = "#c0e6c0")
      textbuffer.insert_with_tags_by_name(start, description, "tag")

   def close_cb_func(*args):
      if post_close_hook:
         post_close_hook()
      window.destroy()

   # main line
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   scrolled_win = gtk.ScrolledWindow()
   outside_vbox = gtk.VBox(False, 2)
   h_sep = gtk.HSeparator()
   inside_vbox = gtk.VBox(False, 0)
   
   window.set_default_size(geometry[0], geometry[1])
   window.set_title(window_name)
   inside_vbox.set_border_width(2)
   window.add(outside_vbox)

   if check_button_label:
      check_button = gtk.CheckButton(check_button_label)
      # somehow need to execute the function before we can use it in the
      # callback. This is odd to say the least. FIXME
      check_button_func(check_button, inside_vbox)
      check_button.connect("toggled", lambda func:
                           check_button_func(check_button, inside_vbox))
      if check_button_is_initially_on_flag:
         check_button.set_active(True)
      outside_vbox.pack_start(check_button, False, False, 2)
   
   outside_vbox.pack_start(scrolled_win, True, True, 0) # expand fill padding
   scrolled_win.add_with_viewport(inside_vbox)
   scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)
   
   for button_info in buttons:
      add_button_info_to_box_of_buttons_vbox(button_info, inside_vbox)
      
   outside_vbox.set_border_width(2)
   outside_vbox.pack_start(h_sep, False, False, 2)
   ok_button = gtk.Button(close_button_label)
   outside_vbox.pack_end(ok_button, False, False, 0)
   ok_button.connect("clicked", close_cb_func, window, post_close_hook)
	
   window.show_all()
   return [inside_vbox, window]

# This is exported outside of the box-of-buttons gui because the
# clear_and_add_back function (e.g. from using the check button)
# needs to add buttons - let's not duplicate that code.
#
def add_button_info_to_box_of_buttons_vbox(button_info, vbox):

   def add_text_to_text_buffer(text_buffer, description):
      start = text_buffer.get_start_iter()
      text_buffer.create_tag("tag", foreground="black", 
                          background = "#c0e6c0")
      text_buffer.insert_with_tags_by_name(start, description, "tag")
   
   button_label = button_info[0]
   if ((button_label == "HSep") and (len(button_info) == 1)):
      # insert a HSeparator rather than a button
      button = gtk.HSeparator()
   else:
      callback = button_info[1]
      if (len(button_info) == 2):
         description = False
      else:
         description = button_info[2]
      button = gtk.Button(button_label)

      # BL says:: in python we should pass the callback as a string
      if type(callback) is StringType:
         def callback_func(button, call):
            eval(call)
         button.connect("clicked", callback_func, callback)
      elif (type(callback) is ListType):
         # list can be list of strings or list of functions
         # with args
         #
         if (isinstance(callback[0], str)):
            # we have strings to evaluate
            def callback_func(button, call):
               for item in call:
                  eval(item)
            button.connect("clicked", callback_func, callback)
         else:
            def callback_func(button, call):
               for item in call:
                  item[0](*item[1:])
            button.connect("clicked", callback_func, callback)
      else:
         button.connect("clicked", callback)

      if (description):
         text_view = gtk.TextView()
         text_view.set_editable(False)
         buff = text_view.get_buffer()
         add_text_to_text_buffer(buff, description)
         vbox.pack_start(text_view, False, False, 2)
      
   vbox.pack_start(button, False, False, 2)
   button.show()


# geometry is an improper list of ints
# buttons is a list of: [[["button_1_label, button_1_action],
#                         ["button_2_label, button_2_action]], [next pair of buttons]]
# The button_1_action function is a string
# The button_2_action function is a string
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
           #print "button_info ", button_info
           if type(button_info) is ListType:
              button_label_1 = button_info[0][0]
              callback_1 = button_info[0][1]

              button_label_2 = button_info[1][0]
              callback_2 = button_info[1][1]

              button_1 = gtk.Button(button_label_1)
              h_box = gtk.HBox(False, 2)

              #print "button_label_1 ", button_label_1
              #print "callback_1 ", callback_1
              #print "button_label_2 ", button_label_2
              #print "callback_2 ", callback_2

              def callback_func(button, call):
                 eval(call)
              button_1.connect("clicked", callback_func, callback_1)
              h_box.pack_start(button_1, False, False, 2)

              if callback_2:
                 button_2 = gtk.Button(button_label_2)
                 button_2.connect("clicked", callback_func, callback_2)
                 h_box.pack_start(button_2, False, False, 2)
              inside_vbox.pack_start(h_box, False, False, 2)

        outside_vbox.set_border_width(2)
        ok_button = gtk.Button(close_button_label)
        outside_vbox.pack_end(ok_button, False, False, 2)
        ok_button.connect("clicked", lambda w: window.destroy())
        window.show_all()

# as the dialog_box_of_buttons, but we can put in an extra widget (extra_widget)
#
def dialog_box_of_buttons_with_widget(window_name, geometry,
                                      buttons, extra_widget,
                                      close_button_label):

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
                elif (type(callback) is ListType):
                   def callback_func(button, call):
                      for item in call:
                         eval(item)
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
#                        text_box.thaw()

                inside_vbox.pack_start(button, False, False, 2)

        # for the extra widget
        inside_vbox.pack_start(h_sep, False, False, 2)
        inside_vbox.pack_start(extra_widget, False, False, 2)
        
        outside_vbox.set_border_width(2)
        ok_button = gtk.Button(close_button_label)
        outside_vbox.pack_end(ok_button, False, False, 0)
        ok_button.connect("clicked", lambda w: window.destroy())
        
        window.show_all()

# A dialog box with radiobuttons e.g. to cycle thru loops
#
# the button list shall be in the form of:
# [[button_label1, "button_function1"],
#  [button_label2, "button_function2"]]
#
# function happens when toggled
# obs: button_functions are strings, but can be tuples for multiple functions
# go_function is string too!
# selected button is the button to be toggled on (default is first)
#
def dialog_box_of_radiobuttons(window_name, geometry, buttons,
                               go_button_label, go_button_function,
                               selected_button = 0,
                               cancel_button_label = "",
                               cancel_function = False):
	
   def go_function_event(widget, button_ls):
      eval(go_button_function)
      window.destroy()
      return False

   def cancel_function_cb(widget):
      #print "BL DEBUG:: just for the sake"
      if (cancel_function):
         eval(cancel_function)
      window.destroy()
      return False
	
   # main line
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   scrolled_win = gtk.ScrolledWindow()
   outside_vbox = gtk.VBox(False, 2)
   inside_vbox = gtk.VBox(False, 0)
   button_hbox = gtk.HBox(False, 0)

   window.set_default_size(geometry[0], geometry[1])
   window.set_title(window_name)
   inside_vbox.set_border_width(2)
   window.add(outside_vbox)
   outside_vbox.pack_start(scrolled_win, True, True, 0) # expand fill padding
   scrolled_win.add_with_viewport(inside_vbox)
   scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)

   button = None
   button_ls = []
   for button_info in buttons:
      button_label = button_info[0]
      callback = button_info[1]
      button = gtk.RadioButton(button, button_label)

      # BL says:: in python we should pass the callback as a string
      if type(callback) is StringType:
         def callback_func(button, call):
            eval(call)
         button.connect("toggled", callback_func, callback)
      else:
         button.connect("toggled", callback)

      inside_vbox.pack_start(button, False, False, 2)
      button_ls.append(button)

   outside_vbox.set_border_width(2)
   go_button     = gtk.Button(go_button_label)
   outside_vbox.pack_start(button_hbox, False, False, 2)
   button_hbox.pack_start(go_button, True, True, 6)
   go_button.connect("clicked", go_function_event, button_ls)
   if (cancel_button_label):
      cancel_button = gtk.Button(cancel_button_label)
      button_hbox.pack_start(cancel_button, True, True, 6)
      cancel_button.connect("clicked", cancel_function_cb)
	
   # switch on the first or selected button
   # somehow I need to emit the toggled signal too (shouldnt have to!?)
   button_ls[selected_button].set_active(True)
   button_ls[selected_button].toggled()

   window.show_all()

global views_dialog_vbox
views_dialog_vbox = False
# A gui showing views
#
def views_panel_gui():

   global views_dialog_vbox
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

   views_vbox = dialog_box_of_buttons("Views", [200,140], buttons, "  Close  ")
   views_dialog_vbox = views_vbox

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
                        print("BL WARNING:: something went wrong!!!")
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
                   rc = list(map(lambda x, y: x + y, rc, nudge_vector))
                elif operand == 1:
                   rc = list(map(lambda x, y: x - y, rc, nudge_vector))
                else:
                        # We should never be here
                        print("BL WARNING:: something went wrong!!!")
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
           rc = list(map(lambda x, y: x + y, rc, nudge_vector))
        elif operand == 1:
           rc = list(map(lambda x, y: x - y, rc, nudge_vector))
        else:
           # We should never be here
           print("BL WARNING:: something went wrong!!!")
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

   dialog_box_of_buttons_with_widget("Nudge Screen Centre with Extras",
                                     [200, 400], buttons, vbox, "  Close ")


# A gui to make a difference map (from arbitrarily gridded maps
# (that's it's advantage))
#
def make_difference_map_gui():
   
   def delete_event(*args):
      window.destroy()
      return False

   def go_function_event(widget):
      print("make diff map here\n")
      active_mol_no_ref = get_option_menu_active_molecule(option_menu_ref_mol, map_molecule_list_ref)
      active_mol_no_sec = get_option_menu_active_molecule(option_menu_sec_mol, map_molecule_list_sec)
      scale_text = scale_entry.get_text()
      scale = False
      try:
         scale = float(scale_text)
      except:
         print("can't decode scale", scale_text)
      if (scale):
         difference_map(active_mol_no_ref, active_mol_no_sec, scale)
      delete_event()

      
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   diff_map_vbox = gtk.VBox(False, 2)
   h_sep = gtk.HSeparator()
   title = gtk.Label("Make a Difference Map")
   ref_label = gtk.Label("Reference Map:")
   sec_label = gtk.Label("Substract this map:")
   second_map_hbox = gtk.HBox(False, 2)
   buttons_hbox = gtk.HBox(True, 6)
   option_menu_ref_mol = gtk.combo_box_new_text()
   option_menu_sec_mol = gtk.combo_box_new_text()
   scale_label = gtk.Label("Scale")
   scale_entry = gtk.Entry()
   ok_button = gtk.Button("   OK   ")
   cancel_button = gtk.Button(" Cancel ")

   map_molecule_list_ref = fill_option_menu_with_map_mol_options(option_menu_ref_mol)
   map_molecule_list_sec = fill_option_menu_with_map_mol_options(option_menu_sec_mol)
  
   window.add(diff_map_vbox)
   diff_map_vbox.pack_start(title, False, False, 2)
   diff_map_vbox.pack_start(ref_label, False, False, 2)
   diff_map_vbox.pack_start(option_menu_ref_mol, True, True, 2)

   diff_map_vbox.pack_start(sec_label, False, False, 2)
   diff_map_vbox.pack_start(second_map_hbox, False, False, 2)

   second_map_hbox.pack_start(option_menu_sec_mol, True, True, 2)
   second_map_hbox.pack_start(scale_label, False, False, 2)
   second_map_hbox.pack_start(scale_entry, False, False, 2)

   diff_map_vbox.pack_start(h_sep, True, False, 2)
   diff_map_vbox.pack_start(buttons_hbox, True, False, 2)
   buttons_hbox.pack_start(ok_button, True, False, 2)
   buttons_hbox.pack_start(cancel_button, True, False, 2)
   scale_entry.set_text("1.0")

   ok_button.connect("clicked", go_function_event)

   cancel_button.connect("clicked", delete_event)

   window.show_all()
   

def cis_peptides_gui(imol):

   def get_ca(atom_list):

      if (atom_list == []):
         return False
      else:
         for atom in atom_list:
            atom_name = atom[0][0]
            if (atom_name == ' CA '):
               print("BL DEBUG:: returning ", atom)
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
            pos = list(map(lambda x, y: (x + y) / 2.0, p1, p2))
            tors_s1 = str(omega)
            if (len(tors_s1) < 6):
               tors_string = tors_s1
            else:
               tors_string = tors_s1[0:6]
            mess = ("Cis Pep: " + chain_id + " " +
                   str(r1[2]) + " " + 
                   residue_name(imol, *r1[1:4]) + " - " +
                   str(r2[2]) + " " +
                   residue_name(imol, *r2[1:4]) + "   " +
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

#
def transform_map_using_lsq_matrix_gui():

  def delete_event(*args):
     window.destroy()
     return False

  def on_ok_button_clicked(*args):
     active_mol_ref = get_option_menu_active_molecule(frame_info_ref[1], frame_info_ref[2])
     active_mol_mov = get_option_menu_active_molecule(frame_info_mov[1], frame_info_mov[2])
     
     chain_id_ref     = frame_info_ref[3].get_text()
     resno_1_ref_text = frame_info_ref[4].get_text()
     resno_2_ref_text = frame_info_ref[5].get_text()

     chain_id_mov     = frame_info_mov[3].get_text()
     resno_1_mov_text = frame_info_mov[4].get_text()
     resno_2_mov_text = frame_info_mov[5].get_text()

     radius_text = radius_entry.get_text()

     imol_map = imol_refinement_map()
     cont = False
     try:
        resno_1_ref = int(resno_1_ref_text)
        resno_2_ref = int(resno_2_ref_text)
        resno_1_mov = int(resno_1_mov_text)
        resno_2_mov = int(resno_2_mov_text)
        radius = float(radius_text)
        cont = True
     except:
        print("BL ERROR:: conversion from input text to numbers failed")
        
     if (cont):
        if (not valid_map_molecule_qm(imol_map)):
           print("Must set the refinement map")
        else:
           imol_copy = copy_molecule(active_mol_mov)
           new_map_number = transform_map_using_lsq_matrix(active_mol_ref, chain_id_ref,
                                                           resno_1_ref, resno_2_ref,
                                                           imol_copy, chain_id_mov,
                                                           resno_1_mov, resno_2_mov,
                                                           imol_map, rotation_centre(), radius)
           set_molecule_name(imol_copy,
                             "Transformed copy of " + strip_path(molecule_name(active_mol_mov)))

           s =  "Transformed map: from map " + str(imol_map) + \
               " by matrix that created coords " + str(imol_copy)
           set_molecule_name(new_map_number,
                             "Transformed map: from map " + str(imol_map) + \
                             " by matrix that created coords " + str(imol_copy))
           
           set_mol_active(imol_copy, 0)
           set_mol_displayed(imol_copy, 0)

     window.destroy()
   
  # atom-sel-type is either 'Reference or 'Moving
  # 
  # return the list [frame, option_menu, model_mol_list, entries...]
  def atom_sel_frame(atom_sel_type):
     frame = gtk.Frame(atom_sel_type)
     # option_menu == combobox
     option_menu = gtk.combo_box_new_text()
     model_mol_list = fill_option_menu_with_coordinates_mol_options(option_menu)
     atom_sel_vbox = gtk.VBox(False, 2)
     atom_sel_hbox = gtk.HBox(False, 2)
     chain_id_label = gtk.Label(" Chain ID ")
     resno_1_label = gtk.Label(" Resno Start ")
     resno_2_label = gtk.Label(" Resno End ")
     chain_id_entry = gtk.Entry()
     resno_1_entry = gtk.Entry()
     resno_2_entry = gtk.Entry()

     frame.add(atom_sel_vbox)
     atom_sel_vbox.pack_start(option_menu,    False, False, 2)
     atom_sel_vbox.pack_start(atom_sel_hbox,  False, False, 2)
     atom_sel_hbox.pack_start(chain_id_label, False, False, 2)
     atom_sel_hbox.pack_start(chain_id_entry, False, False, 2)
     atom_sel_hbox.pack_start(resno_1_label,  False, False, 2)
     atom_sel_hbox.pack_start(resno_1_entry,  False, False, 2)
     atom_sel_hbox.pack_start(resno_2_label,  False, False, 2)
     atom_sel_hbox.pack_start(resno_2_entry,  False, False, 2)

     return [frame, option_menu, model_mol_list, chain_id_entry, resno_1_entry, resno_2_entry]
  
  window = gtk.Window(gtk.WINDOW_TOPLEVEL)
  dialog_name = "Map Transformation"
  main_vbox = gtk.VBox(False, 2)
  buttons_hbox = gtk.HBox(False, 2)
  cancel_button = gtk.Button("  Cancel  ")
  ok_button = gtk.Button("  Transform  ")
  usage = "Note that this will transform the current refinement map " + \
          "to around the screen centre"
  usage_label = gtk.Label(usage)
  h_sep = gtk.HSeparator()
  frame_info_ref = atom_sel_frame("Reference")
  frame_info_mov = atom_sel_frame("Moving")
  radius_hbox = gtk.HBox(False, 2)
  radius_label = gtk.Label("  Radius ")
  radius_entry = gtk.Entry()
  window.set_title(dialog_name)


  radius_hbox.pack_start(radius_label, False, False, 2)
  radius_hbox.pack_start(radius_entry, False, False, 2)

  buttons_hbox.pack_start(ok_button,     False, False, 4)
  buttons_hbox.pack_start(cancel_button, False, False, 4)

  window.add(main_vbox)
  main_vbox.pack_start(frame_info_ref[0], False, False, 2)
  main_vbox.pack_start(frame_info_mov[0], False, False, 2)
  main_vbox.pack_start(radius_hbox, False, False, 2)
  main_vbox.pack_start(usage_label, False, False, 4)
  main_vbox.pack_start(h_sep, False, False, 2)
  main_vbox.pack_start(buttons_hbox, False, False, 6)

  frame_info_ref[3].set_text("A")
  frame_info_ref[4].set_text("1")
  frame_info_ref[5].set_text("10")
  frame_info_mov[3].set_text("B")
  frame_info_mov[4].set_text("1")
  frame_info_mov[5].set_text("10")

  radius_entry.set_text("8.0")

  cancel_button.connect("clicked", delete_event)
  ok_button.connect("clicked",on_ok_button_clicked)

  window.show_all()
  if (not valid_map_molecule_qm(imol_refinement_map())):
     show_select_map_dialog()


def ncs_ligand_gui():
   
   def delete_event(*args):
      window.destroy()
      return False

   def go_button_function(widget):
      print("ncs ligand function here\n")
      active_mol_no_ref = get_option_menu_active_molecule(option_menu_ref_mol, molecule_list_ref)
      active_mol_no_lig = get_option_menu_active_molecule(option_menu_lig_mol, molecule_list_lig)
      chain_id_lig = chain_id_lig_entry.get_text()
      chain_id_ref = chain_id_ref_entry.get_text()
      resno_start = False
      resno_end = False
      try:
         resno_start = int(resno_start_entry.get_text())
      except:
         print("can't decode resno_start", resno_start_entry.get_text())
         
      resno_end_t = resno_end_entry.get_text()
      try:
         resno_end = int(resno_end_t)
      except:
         if (resno_end_t == ""):
            resno_end = resno_start
         else:
            print("can't decode resno_end", resno_end_t)

      if (resno_end and resno_start):
         make_ncs_ghosts_maybe(active_mol_no_ref)
         print("ncs ligand with", active_mol_no_ref, \
               chain_id_ref, active_mol_no_lig, chain_id_lig, \
               resno_start, resno_end)
         ncs_ligand(active_mol_no_ref,
                    chain_id_ref,
                    active_mol_no_lig,
                    chain_id_lig,
                    resno_start,
                    resno_end)

      delete_event()
      
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   ncs_ligands_vbox = gtk.VBox(False, 2)
   title = gtk.Label("Find NCS-Related Ligands")
   ref_label = gtk.Label("Protein with NCS:")
   ref_chain_hbox = gtk.HBox(False, 2)
   chain_id_ref_label = gtk.Label("NCS Master Chain")
   chain_id_ref_entry = gtk.Entry()
   lig_label = gtk.Label("Molecule containing ligand")
   specs_hbox = gtk.HBox(False, 2)
   h_sep = gtk.HSeparator()
   buttons_hbox = gtk.HBox(True, 6)
   chain_id_lig_label= gtk.Label("Chain ID: ")
   resno_start_label = gtk.Label(" Residue Number ")
   to_label = gtk.Label("  to  ")
   chain_id_lig_entry = gtk.Entry()
   resno_start_entry = gtk.Entry()
   resno_end_entry = gtk.Entry()
   ok_button = gtk.Button("   Find Candidate Positions  ")
   cancel_button = gtk.Button("    Cancel    ")
   option_menu_ref_mol = gtk.combo_box_new_text()
   option_menu_lig_mol = gtk.combo_box_new_text()

   molecule_list_ref = fill_option_menu_with_coordinates_mol_options(option_menu_ref_mol)
   molecule_list_lig = fill_option_menu_with_coordinates_mol_options(option_menu_lig_mol)

   window.add(ncs_ligands_vbox)
   ncs_ligands_vbox.pack_start(title, False, False, 6)
   ncs_ligands_vbox.pack_start(ref_label, False, False, 2)
   ncs_ligands_vbox.pack_start(option_menu_ref_mol, True, False, 2)
   ncs_ligands_vbox.pack_start(ref_chain_hbox, False, False, 2)
   ncs_ligands_vbox.pack_start(lig_label, False, False, 2)
   ncs_ligands_vbox.pack_start(option_menu_lig_mol, True, False, 2)
   ncs_ligands_vbox.pack_start(specs_hbox, False, False, 2)
   ncs_ligands_vbox.pack_start(h_sep, False, False, 2)
   ncs_ligands_vbox.pack_start(buttons_hbox, False, False, 2)

   buttons_hbox.pack_start(ok_button,     True, False, 4)
   buttons_hbox.pack_start(cancel_button, True, False, 4)

   ref_chain_hbox.pack_start(chain_id_ref_label, False, False, 2)
   ref_chain_hbox.pack_start(chain_id_ref_entry, False, False, 2)

   specs_hbox.pack_start(chain_id_lig_label, False, False, 2)
   specs_hbox.pack_start(chain_id_lig_entry, False, False, 2)
   specs_hbox.pack_start(resno_start_label, False, False, 2)
   specs_hbox.pack_start(resno_start_entry, False, False, 2)
   specs_hbox.pack_start(to_label, False, False, 2)
   specs_hbox.pack_start(resno_end_entry, False, False, 2)
   specs_hbox.pack_start(gtk.Label(" "), False, False, 2) # neatness ?!

   chain_id_lig_entry.set_size_request(32, -1)
   chain_id_ref_entry.set_size_request(32, -1)
   resno_start_entry.set_size_request(50, -1)
   resno_end_entry.set_size_request(50, -1)
   chain_id_ref_entry.set_text("A")
   chain_id_lig_entry.set_text("A")
   resno_start_entry.set_text("1")

   #tooltips = gtk.Tooltips()
   chain_tip = "'A' is a reasonable guess at the NCS master chain id.  " + \
               "If your ligand (specified below) is NOT bound to the protein's " + \
               "'A' chain, then you will need to change this chain and also " + \
               "make sure that the master molecule is specified appropriately " + \
               "in the Draw->NCS Ghost Control window."
   resno_tip = "Leave blank for a single residue"
   if gtk.pygtk_version >= (2,12):
      chain_id_ref_entry.set_tooltip_text(chain_tip)
      resno_end_entry.set_tooltip_text(resno_tip)
   else:
      coot_tooltips.set_tip(chain_id_ref_entry, chain_tip)
      coot_tooltips.set_tip(resno_end_entry, resno_tip)

   ok_button.connect("clicked", go_button_function)

   cancel_button.connect("clicked", delete_event)

   window.show_all()


# NCS jumping GUI
global ncs_jumping_time_step
ncs_jumping_time_step = 500

def ncs_jumping_gui():

   global ncs_jumping_time_step
   global timeout_function_token
   timeout_function_token = False

   def delete_event(*args):
      # first stop the jumping (if necessary)
      global timeout_function_token
      timeout_function_token = False
      
      window.destroy()
      return False

   # FIXME chekc this. Not sure if we can get a number from timeout_add or
   # if we better make a new function which returns True/False to continue/stop

   # need to return True to be called again. I return False if stop (bug in
   # gobject.source_remove prevents me from using this.
   def skip_ncs_timeout_func():
      global timeout_function_token
      skip_to_next_ncs_chain("forward")
      if timeout_function_token:
         return True
      else:
         return False
      
   def start_function_event(*args):
      global timeout_function_token
      if not isNumber(timeout_function_token):
         timeout_function_token = gobject.timeout_add(ms_step,
                                                      skip_ncs_timeout_func)
      else:
         timeout_function_token = False

         
   def stop_function_event(*args):
      global timeout_function_token
      timeout_function_token = False      

   # main body
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   outside_vbox = gtk.VBox(False, 2)
   inside_hbox = gtk.HBox(False, 2)
   cancel_hbox = gtk.HBox(False, 2) # paul says VBox?!?!
   h_sep = gtk.HSeparator()
   jump_start_button = gtk.Button("NCS Jump Start")
   jump_stop_button = gtk.Button("Stop")
   cancel_button = gtk.Button("Cancel")
   ms_step = ncs_jumping_time_step
   timeout_function_token = False

   window.set_title("Auto NCS Jumping")
   window.add(outside_vbox)
   outside_vbox.pack_start(inside_hbox, False, False, 2)
   outside_vbox.pack_start(h_sep, False, False, 2)
   outside_vbox.pack_start(cancel_hbox, False, False, 2)
   inside_hbox.pack_start(jump_start_button, False, False, 2)
   inside_hbox.pack_start(jump_stop_button, False, False, 2)
   cancel_hbox.pack_start(cancel_button, False, False, 2)

   jump_start_button.connect("clicked", start_function_event)

   jump_stop_button.connect("clicked", stop_function_event)
   
   cancel_button.connect("clicked", delete_event)

   window.show_all()

   
# GUI for ligand superpositioning by graph matching
#
def superpose_ligand_gui():
   
   def delete_event(*args):
      window.destroy()
      return False

   def go_button_function(widget):
      active_mol_no_ref_lig = get_option_menu_active_molecule(*option_menu_ref_mol_pair)
      active_mol_no_mov_lig = get_option_menu_active_molecule(*option_menu_mov_mol_pair)
      chain_id_ref = chain_id_ref_entry.get_text()
      chain_id_mov = chain_id_mov_entry.get_text()
      resno_ref = False
      resno_mov = False
      try:
         resno_ref = int(resno_ref_entry.get_text())
      except:
         print("can't decode reference resno", resno_ref_entry.get_text())
         
      try:
         resno_mov = int(resno_mov_entry.get_text())
      except:
         print("can't decode moving resno", resno_mov_entry.get_text())
         

      if (resno_ref and resno_mov):
         overlay_my_ligands(active_mol_no_mov_lig, chain_id_mov, resno_mov,
                            active_mol_no_ref_lig, chain_id_ref, resno_ref)

      delete_event()
      
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   title = gtk.Label("Superpose Ligands")
   ligands_vbox = gtk.VBox(False, 2)
   ref_chain_hbox = gtk.HBox(False, 2)
   chain_id_ref_label = gtk.Label("Ligand Chain ID: ")
   chain_id_ref_entry = gtk.Entry()
   resno_ref_label    = gtk.Label(" Residue Number ") 
   resno_ref_entry = gtk.Entry()
   
   mov_chain_hbox = gtk.HBox(False, 2)
   chain_id_mov_label= gtk.Label("Ligand Chain ID: ")
   chain_id_mov_entry = gtk.Entry()
   resno_mov_label = gtk.Label(" Residue Number ")
   resno_mov_entry = gtk.Entry()

   h_sep = gtk.HSeparator()

   buttons_hbox = gtk.HBox(True, 6)

   ok_button = gtk.Button("   Superpose 'em  ")
   cancel_button = gtk.Button("    Cancel    ")

   window.add(ligands_vbox)
   ligands_vbox.pack_start(title, False, False, 6)
   option_menu_ref_mol_pair = generic_molecule_chooser(ligands_vbox, "Model with reference ligand")
   ligands_vbox.pack_start(ref_chain_hbox, False, False, 2)

   option_menu_mov_mol_pair = generic_molecule_chooser(ligands_vbox, "Model with moving ligand")   
   ligands_vbox.pack_start(mov_chain_hbox, False, False, 2)

   ligands_vbox.pack_start(h_sep, False, False, 2)
   ligands_vbox.pack_start(buttons_hbox, False, False, 2)

   buttons_hbox.pack_start(ok_button,     True, False, 4)
   buttons_hbox.pack_start(cancel_button, True, False, 4)

   ref_chain_hbox.pack_start(chain_id_ref_label, False, False, 2)
   ref_chain_hbox.pack_start(chain_id_ref_entry, False, False, 2)
   ref_chain_hbox.pack_start(resno_ref_label, False, False, 2)
   ref_chain_hbox.pack_start(resno_ref_entry, False, False, 2)

   mov_chain_hbox.pack_start(chain_id_mov_label, False, False, 2)
   mov_chain_hbox.pack_start(chain_id_mov_entry, False, False, 2)
   mov_chain_hbox.pack_start(resno_mov_label, False, False, 2)
   mov_chain_hbox.pack_start(resno_mov_entry, False, False, 2)

#   chain_id_lig_entry.set_size_request(32, -1)
#   chain_id_ref_entry.set_size_request(32, -1)
#   resno_start_entry.set_size_request(50, -1)
#   resno_end_entry.set_size_request(50, -1)
#   chain_id_ref_entry.set_text("A")
#   chain_id_lig_entry.set_text("A")
#   resno_start_entry.set_text("1")

#   tooltips = gtk.Tooltips()
#   tooltips.set_tip(chain_id_ref_entry, "'A' is a reasonable guess at the NCS master chain id.  " +
#                    "If your ligand (specified below) is NOT bound to the protein's " +
#                    "'A' chain, then you will need to change this chain and also " +
#                    "make sure that the master molecule is specified appropriately " +
#                    "in the Draw->NCS Ghost Control window.")
#   tooltips.set_tip(resno_end_entry, "Leave blank for a single residue")

   ok_button.connect("clicked", go_button_function)

   cancel_button.connect("clicked", delete_event)

   window.show_all()


def gui_overlap_ligands(imol_ligand, imol_ref, chain_id_ref, res_no_ref):

    # we don't want to overlap-ligands if there is no dictionary
    # for the residue to be matched to.
    #
    res_name = residue_name(imol_ref, chain_id_ref, res_no_ref, "")
    restraints = monomer_restraints(res_name)
    if (not restraints):
        return False
    else:
        if (not residue_has_hetatms_qm(imol_ref, chain_id_ref, res_no_ref, "")):
            return False
        else:
            print("----------- overlap-ligands %s %s %s %s ------------" \
            %(imol_ligand, imol_ref, chain_id_ref, res_no_ref))
            # this can return the rtop operator or the False (for fail of course).
            match_ligand_torsions(imol_ligand, imol_ref, chain_id_ref, res_no_ref)
            ret = overlap_ligands(imol_ligand, imol_ref, chain_id_ref, res_no_ref)
            return ret


global std_key_bindings
std_key_bindings = [["^g", "keyboard-go-to-residue"],
                    ["^s", "quick-save-as"],
                    ["^i", "residue info"],
                    ["^z", "undo"],
                    ["^y", "redo"],
                    ["a", "refine with auto-zone"],
                    ["b", "toggle baton swivel"],
                    ["c", "toggle cross-hairs"],
                    ["d", "reduce depth of field"],
                    ["f", "increase depth of field"],
                    ["u", "undo last navigation"],
                    ["i", "toggle spin mode"],
                    ["l", "label closest atom"],
                    ["m", "zoom out"],
                    ["n", "zoom in"],
                    ["o", "other NCS chain"],
                    ["p", "update position to closest atom"],
                    ["s", "update skeleton"],
                    [".", "up in button list"],
                    [",", "down in button list"]]
   
def key_bindings_gui():

   global std_key_bindings
   def delete_event(*args):
      window.destroy()
      return False

   def box_for_binding(item, inside_vbox, buttonize_flag):
      
      binding_hbox = gtk.HBox(False, 2)
      txt = str(item[1])
      key_label = gtk.Label("   " + txt + "   ")
      name_label = gtk.Label(item[2])

      if (buttonize_flag):
         button_label = "   " + txt + "   " + item[2]
         button = gtk.Button(button_label)
         #al = gtk.Alignment(0, 0, 0, 0)
         #label = gtk.Label(button_label)
         #button.add(al)
         #al.add(label)
         binding_hbox.pack_start(button, True, True, 0)
         inside_vbox.pack_start(binding_hbox, False, False, 0)
         binding_func = item[3]
         if not (callable(binding_func)):
            s  = "Cannot call given function with button,\n"
            s += "probably a scheme function.\n"
            s += "The shortcut should still work though."
            def binding_func():
               info_dialog(s)
               print("INFO::", s)
         
         button.connect("clicked", lambda func: apply(binding_func))
         
      else:
         binding_hbox.pack_start(key_label, False, False, 2)
         binding_hbox.pack_start(name_label, False, False, 2)
         inside_vbox.pack_start(binding_hbox, False, False, 2)
      
   # main line
   #
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   scrolled_win = gtk.ScrolledWindow()
   outside_vbox = gtk.VBox(False, 2)
   inside_vbox = gtk.VBox(False, 0)
   dialog_name = "Key Bindings"
   buttons_hbox = gtk.HBox(False, 2)
   close_button = gtk.Button("  Close  ")
   std_frame = gtk.Frame("Standard Key Bindings:")
   usr_frame = gtk.Frame("User-defined Key Bindings:")
   std_frame_vbox = gtk.VBox(False, 2)
   usr_frame_vbox = gtk.VBox(False, 2)
   close_button.connect("clicked", delete_event)

   window.set_default_size(250, 350)
   window.set_title(dialog_name)
   inside_vbox.set_border_width(4)

   window.add(outside_vbox)
   outside_vbox.add(scrolled_win)
   scrolled_win.add_with_viewport(inside_vbox)
   scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)

   inside_vbox.pack_start(std_frame, False, False, 2)
   inside_vbox.pack_start(usr_frame, False, False, 2)
   
   std_frame.add(std_frame_vbox)
   usr_frame.add(usr_frame_vbox)

   py_and_scm_keybindings = key_bindings
   if (coot_has_guile()):
      scm_key_bindings = run_scheme_command("*key-bindings*")
      # check for list
      if isinstance(scm_key_bindings, list):
         # filter out doublicates
         for item in scm_key_bindings:
            scm_code, scm_key, text, tmp = item
            py_keys  = [elem[1] for elem in py_and_scm_keybindings]
            py_codes = [elem[0] for elem in py_and_scm_keybindings]
            if ((not scm_code in py_codes) and (not scm_key in py_keys)):
               item[2] = item[2] + " (scm)"
               py_and_scm_keybindings.append(item)
            else:
               item[2] = item[2] + " (scm + doublicate key)"
               py_and_scm_keybindings.append(item)

       
   for items in py_and_scm_keybindings:
      box_for_binding(items, usr_frame_vbox, True)

   for items in std_key_bindings:
      box_for_binding(["dummy"] + items, std_frame_vbox, False)

   buttons_hbox.pack_end(close_button, False, False, 6)
   outside_vbox.pack_start(buttons_hbox, False, False, 6)

   window.show_all()

# for news infos
(
   INSERT_NO_NEWS,
   INSERT_NEWS,
   STATUS,
   NEWS_IS_READY,
   GET_NEWS,
   SET_TEXT,
   STOP,
   NO_NEWS,
   ) = list(range(8))

global news_status
news_status = NO_NEWS
global news_string_1
global news_string_2
news_string_1 = False
news_string_2 = False

def coot_news_info(*args):

   import threading
   import urllib.request, urllib.parse, urllib.error
   global text_1, text_2
   url = "http:" + \
         "//www.biop.ox.ac.uk/coot/software" + \
         "/binaries/pre-releases/PRE-RELEASE-NOTES"

   def test_string():
      import time
      time.sleep(2)
      ret = "assssssssssssssssssssssssssssssssssssssss\n\n" + \
            "assssssssssssssssssssssssssssssssssssssss\n\n" + \
            "assssssssssssssssssssssssssssssssssssssss\n\n" + \
            "\n-----\n" + \
            "bill asssssssssssssssssssssssssssssssssssss\n\n" + \
            "fred asssssssssssssssssssssssssssssssssssss\n\n" + \
            "george sssssssssssssssssssssssssssssssssssss\n\n" + \
            "\n-----\n"
      return ret

   def stop():
      return

   # return [pre_release_news_string, std_release_news_string]
   def trim_news(s):
      sm_pre = s.find("-----")
      if (sm_pre == -1):
         return ["nothing", "nothing"]
      else:
         pre_news = s[0:sm_pre]
         post_pre = s[sm_pre:-1].lstrip("-")
         sm_std = post_pre.find("-----")
         if (sm_std == -1):
            return [pre_news, "nothing"]
         else:
            return [pre_news, post_pre[0:sm_std]]

   def get_news_thread():
      global news_status
      global news_string_1
      global news_string_2
      import urllib.request, urllib.parse, urllib.error
      try:
         s = urllib.request.urlopen(url).read()
         both_news = trim_news(s)
         news_string_1 = both_news[1]
         news_string_2 = both_news[0]
         news_status = NEWS_IS_READY
      except:
         pass

   def coot_news_error_handle(key, *args):
      # not used currently
      print("error: news: error in %s with args %s" %(key, args))

   def get_news():
      # no threading for now. Doesnt do the right thing
      run_python_thread(get_news_thread, ())
      
   def insert_string(s, text):
      background_colour = "#c0e6c0"
      end = text.get_end_iter()
      text.insert(end, str(s))

   def insert_news():
      global news_string_1
      global news_string_2
      insert_string(news_string_1, text_1)
      insert_string(news_string_2, text_2)
      
   def insert_no_news():
      insert_string("  No news\n", text_1)
      insert_string("  Yep - there really is no news\n", text_2)

   if (len(args) == 1):
      arg = args[0]
      if (arg == STOP):
         stop()
      if (arg == STATUS):
         return news_status
      if (arg == INSERT_NEWS):
         insert_news()
      if (arg == INSERT_NO_NEWS):
         insert_no_news()
      if (arg == GET_NEWS):
         get_news()
   if (len(args) == 3):
      if (args[0] == SET_TEXT):
         text_1 = args[1]
         text_2 = args[2]

def whats_new_dialog():
   global text_1, text_2
   text_1 = False         # the text widget
   text_2 = False
   timer_label = False
   global count
   count = 0
   ms_step = 200

   def on_close_button_clicked(*args):
      coot_news_info(STOP)
      delete_event(*args)
      
   def delete_event(*args):
      window.destroy()
      return False
   
   def check_for_new_news():
      global count
      count += 1
      timer_string = str(count * ms_step / 1000.) + "s"
      timer_label.set_alignment(0.96, 0.5)
      timer_label.set_text(timer_string)
      if (count > 100):
         timer_label.set_text("Timeout")
         coot_news_info(INSERT_NO_NEWS)
         return False  # turn off the gtk timeout function ?!
      else:
         if (coot_news_info(STATUS) == NEWS_IS_READY):
            coot_news_info(INSERT_NEWS)
            return False
         return True

   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   vbox = gtk.VBox(False, 2)
   inside_vbox = gtk.VBox(False, 2)
   scrolled_win_1 = gtk.ScrolledWindow()
   scrolled_win_2 = gtk.ScrolledWindow()
   label = gtk.Label("Lastest Coot Release Info")
   text_view_1 = gtk.TextView()
   text_view_2 = gtk.TextView()
   text_view_1.modify_base(gtk.STATE_NORMAL, gtk.gdk.color_parse("#bfe6bf"))
   text_view_2.modify_base(gtk.STATE_NORMAL, gtk.gdk.color_parse("#bfe6bf"))
   text_1 = text_view_1.get_buffer()
   text_2 = text_view_2.get_buffer()
   h_sep = gtk.HSeparator()
   close_button = gtk.Button("   Close   ")
   notebook = gtk.Notebook()
   notebook_label_pre = gtk.Label("Pre-release")
   notebook_label_std = gtk.Label("Release")
   timer_label = gtk.Label("0.0s")

   window.set_default_size(540, 400)
   vbox.pack_start(label, False, False, 10)
   vbox.pack_start(timer_label, False, False, 2)
   vbox.pack_start(notebook, True, True, 4)
   notebook.append_page(scrolled_win_1, notebook_label_std)
   notebook.append_page(scrolled_win_2, notebook_label_pre)
   vbox.pack_start(h_sep, False, False, 4)
   vbox.pack_start(close_button, False, False, 2)
   window.add(vbox)

   coot_news_info(SET_TEXT, text_1, text_2)
   coot_news_info(GET_NEWS)

   scrolled_win_1.add(text_view_1)
   scrolled_win_2.add(text_view_2)
   scrolled_win_1.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
   scrolled_win_2.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)

   close_button.connect("clicked", on_close_button_clicked)

   gobject.timeout_add(ms_step, check_for_new_news)

   window.show_all()
   

# Cootaneer/sequencing gui modified by BL with ideas from KC
# based on Paul's cootaneer gui and generic_chooser_entry_and_file_selector
#
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
         print("INFO:: refinement on")
      else:
         print("INFO:: refinement off")

   def go_function_event(widget):
      print("apply the sequence info here\n")
      print("then cootaneer\n")

      # no active atom won't do.  We need
      # to find the nearest atom in imol to (rotation-centre).
      #
      # if it is too far away, give a
      # warning and do't do anything.
      active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)
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
            at_name = residue_spec_to_atom_for_centre(imol, chain_id, res_no, ins_code)[0]
            cootaneer_result = cootaneer(imol_map, imol, [chain_id, res_no, ins_code, 
                                                          at_name, alt_conf])
            if (cootaneer_result == 0):
               s = "Insufficiently confident in alignment to make a fit." + \
                   "\n" + \
                   "Perhaps you could improve or extend this fragment."
               info_dialog(s)
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
      active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)
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
            widget_range = list(range(no_of_sequences))
         else:
            # we update the number of sequences
            spin_len = int(spin_button.get_value())
            widget_range = list(range(spin_len, no_of_sequences))

         # make new table
         imported_sequence_file_flags = [True, no_of_sequences]
         spin_button.set_value(no_of_sequences)
         seq_table.resize(no_of_sequences, 1)
         for i in widget_range:
            seq_widget = entry_text_pair_frame_with_button(seq_info_ls[i])
            seq_table.attach(seq_widget[0], 0, 1, i, i+1)
            seq_widget[0].show_all()
      else:
         print("BL WARNING:: no filename")

   def fill_table_with_sequences(*args):
      # fills the table with sequences if they have been associated with the model imol
      # already
      global imported_sequence_file_flags
      active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)
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
         widget_range = list(range(no_of_sequences))
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
         active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)
         imol = int(active_mol_no)
         imol_map = imol_refinement_map()
         print("apply the sequence info here\n")
         print("then cootaneer\n")

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
               chain_id = n_atom[1]
               res_no	= n_atom[2]
               ins_code = n_atom[3]
               at_name	= n_atom[4]
               alt_conf = n_atom[5]
               cootaneer_result = cootaneer(imol_map, imol, [chain_id, res_no, ins_code, 
                                                             at_name, alt_conf])
               if (cootaneer_result == 0):
                  s = "Insufficiently confident in alignment to make a fit." + \
                      "\n" + \
                      "Perhaps you could improve or extend this fragment."
                  info_dialog(s)
               else:
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
               print("BL WARNING:: no close atom found!")
         else:
            print("BL ERROR:: something went wrong assigning the sequence")

      def chain_toggled(widget):
         # doesnt need to do anything
         status = chain_check_button.get_active()
         if(status):
            print("BL INFO:: assign chain_id too")
         else:
            print("BL INFO:: do not assign chain_id")

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
         redraw_range = list(range(no_of_sequences, no_of_frames))
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
               child_range = list(range(0, no_of_children - no_of_sequences))
               # children seem to be added in the beginning ???????
      else:
         # we make everything new
         redraw_range = list(range(no_of_frames))
         child_range = list(range(no_of_children))

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
   def clear_function_event(widget = None, file_sel_entry = None):
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
      # reset the filechooser entry
      if file_sel_entry:
         file_sel_entry.set_text("")

   # make one cell in line with deafult fill
   def make_cell(line):
      seq_widget = entry_text_pair_frame_with_button(["",
                                                      "Cut and Paste Sequence to here or import a sequence file"])
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
            print("ERROR:: the input contains an invalid chain and/or sequence")
            write_sequence = False
         seq_all.append(pair)

      if (write_sequence):
         for element in seq_all:
            chain_id_new = element[0]
            seq = element[1].upper()
            # first check if chain_id is already in mol
            # if so delete it so that it can be replaced by the new sequence
            seq_info = sequence_info(imol)
            if seq_info:
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
   #tooltips = gtk.Tooltips()
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
   if gtk.pygtk_version >= (2,12):
      go_button.set_tooltip_text("This currently ignores all chain IDs")
   else:
      coot_tooltips.set_tip(go_button, "This currently ignores all chain IDs")
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
   file_sel_entry = file_selector_entry(vbox, "Select PIR file")
   vbox.pack_start(import_button, False, False, 6)

   buttons_hbox.pack_start(go_button, False, False, 6)
   buttons_hbox.pack_start(cancel_button, False, False, 6)
   buttons_hbox.pack_start(clear_button, False, False, 6)

   vbox.pack_start(h_sep2, False, False, 2)
   vbox.pack_start(buttons_hbox, False, False, 5)


   import_button.connect("clicked", import_function_event, file_sel_entry)

   cancel_button.connect("clicked", delete_event)

   go_button.connect("clicked", go_function_event)

   clear_button.connect("clicked", clear_function_event, file_sel_entry)

   spin_adj.connect("value_changed", spin_button_changed)

   refine_check_button.connect("toggled", refine_function_event)

   option_menu.connect("changed", fill_table_with_sequences)

#        window.add(vbox)
   window.show_all()

# a function to run a pygtk widget in a function as a thread
#
def run_with_gtk_threading(function, *args):
   import gobject
   def idle_func():
      gtk.gdk.threads_enter()
      try:
         # function(*args, **kw)
         function(*args)
         return False
      finally:
         gtk.gdk.threads_leave()
   gobject.idle_add(idle_func)


def generic_check_button(vbox, label_text, handle_check_button_function):
   def check_callback(*args):
      active_state = check_button.get_active()
      set_state = 0
      if (active_state):
         set_state = 1
      handle_check_button_function(set_state)
   check_button = gtk.CheckButton(label_text)
   vbox.pack_start(check_button, False, False, 2)
   check_button.connect("toggled", check_callback)
   return check_button

# a master gui to set all kinds of refinement parameters
#
def refinement_options_gui():

   def set_matrix_func(button, entry):
      text = entry.get_text()
      try:
         t = float(text)
         s = "Matrix set to " + text
         set_matrix(t)
         add_status_bar_text(s)
      except:
         add_status_bar_text("Failed to read a number") 

   def delete_event(*rags):
      window.destroy()
      return False

   def go_function_event(*args):
      set_matrix_func('dummy', matrix_entry)
      window.destroy()
      return False

   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   vbox = gtk.VBox(False, 0)
   hbox = gtk.HBox(False, 0)
   h_sep = gtk.HSeparator()
   h_sep2 = gtk.HSeparator()
   h_sep3 = gtk.HSeparator()
   go_button = gtk.Button("   Ok   ")
   cancel_button = gtk.Button("  Cancel  ")

   window.add(vbox)
   # add the matrix entry
   matrix_entry = entry_do_button(vbox, "set matrix: (smaller means better geometry)",
                                  "Set", set_matrix_func)
   matrix_entry.set_text(str(matrix_state()))

   vbox.pack_start(h_sep2, False, False, 2)

   # use torsion restrains?
   torsion_restraints_button = generic_check_button(vbox,
                                                 "Use Torsion Restraints?",
                                                 lambda state: set_refine_with_torsion_restraints(state))
   if (refine_with_torsion_restraints_state() == 1):
      torsion_restraints_button.set_active(True)
   else:
      torsion_restraints_button.set_active(False)

   # planar peptide restrains?
   planar_restraints_button = generic_check_button(vbox,
                                                 "Use Planar Peptide Restraints?",
                                                 lambda state: remove_planar_peptide_restraints() if state == 0 else add_planar_peptide_restraints())
   if (planar_peptide_restraints_state() == 1):
      planar_restraints_button.set_active(True)
   else:
      planar_restraints_button.set_active(False)

   # use ramachandran restrains?
   rama_restraints_button = generic_check_button(vbox,
                                                 "Use Ramachandran Restraints?",
                                                 lambda state: set_refine_ramachandran_angles(state))
   if (refine_ramachandran_angles_state() == 1):
      rama_restraints_button.set_active(True)
   else:
      rama_restraints_button.set_active(False)      

   vbox.pack_start(h_sep3, False, False, 2)

   # add rotamer check button
   rotamer_autofit_button = generic_check_button(vbox,
                                                 "Real Space Refine after Auto-fit Rotamer?",
                                                 lambda state: set_rotamer_auto_fit_do_post_refine(state))
   if (rotamer_auto_fit_do_post_refine_state() == 1):
      rotamer_autofit_button.set_active(True)
   else:
      rotamer_autofit_button.set_active(False)

   # add mutate check button
   mutate_autofit_button = generic_check_button(vbox,
                                                "Real Space Refine after Mutate and Auto-fit?",
                                                lambda state: set_mutate_auto_fit_do_post_refine(state))
   if (mutate_auto_fit_do_post_refine_state() == 1):
      mutate_autofit_button.set_active(True)
   else:
      mutate_autofit_button.set_active(False)

   # add terminal residue check button
   terminal_autofit_button = generic_check_button(vbox,
                                                  "Real Space Refine after Add Terminal Residue?",
                                                  lambda state: set_add_terminal_residue_do_post_refine(state))
   if (add_terminal_residue_do_post_refine_state() == 1):
      terminal_autofit_button.set_active(True)
   else:
      terminal_autofit_button.set_active(False)

   # add a b-factor button
   reset_b_factor_button = generic_check_button(vbox,
                                                "Reset B-Factors to Default Value after Refinement/Move?",
                                                lambda state: set_reset_b_factor_moved_atoms(state))
   if (get_reset_b_factor_moved_atoms_state() == 1):
      reset_b_factor_button.set_active(True)
   else:
      reset_b_factor_button.set_active(False)

   vbox.pack_start(h_sep, False, False, 2)
   vbox.pack_start(hbox, False, False, 0)
   hbox.pack_start(go_button, False, False, 6)
   hbox.pack_start(cancel_button, False, False, 6)

   go_button.connect("clicked", go_function_event)
   cancel_button.connect("clicked", delete_event)
   window.show_all()

# a simple window to show a progress bar
# return the window (to be destroyed elsewhere)
#
def show_progressbar(text):

   gtk.gdk.threads_init()
   def progress_run(pbar):
      pbar.pulse()
      return True

   def destroy_cb(widget, timer):
      gobject.source_remove(timer)
      timer = 0
      gtk.main_quit()
      return False
   
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   window.set_title("External Program Progress")
   window.set_border_width(0)
   window.set_default_size(300, 50)

   vbox = gtk.VBox(False, 5)
   vbox.set_border_width(10)
   window.add(vbox)

   pbar = gtk.ProgressBar()
   pbar.pulse()
   pbar.set_text(text)
   vbox.pack_start(pbar, False, False, 5)

   timer = gobject.timeout_add (100, progress_run, pbar)
   
   window.connect("destroy", destroy_cb, timer)

   window.show_all()
   global python_return
   python_return = window
   gtk.main()

import threading
# helper function to push the python threads
# this only runs python threads for 20ms every 50ms
def python_thread_sleeper():
   global python_thread_sleep
   sleep_time = python_thread_sleep
   time.sleep(sleep_time / 1000.)
   if (threading.activeCount() == 1):
     #print "BL DEBUG:: stopping timeout"
     return False
   return True

# this has locked, so that no one else can use it
global python_thread_return
python_thread_return = False
global python_thread_sleep
python_thread_sleep = 20

# function to run a python thread with function using
# args which is a tuple
# optionally pass sleep time in ms (default is 20) - usefull
# for computationally expensive threads which may have run longer
# N.B. requires gobject hence in coot_gui.py
#
def run_python_thread(function, args, sleep_time=20):

   import gobject
   
   class MyThread(threading.Thread):
      def __init__(self):
         threading.Thread.__init__(self)
      def run(self):
         global python_thread_return
         python_return_lock = threading.Lock()
         python_return_lock.acquire()
         try:
            python_thread_return = function(*args)
         finally:
            python_return_lock.release()

   global python_thread_sleep
   if (not sleep_time == 20):
      python_thread_sleep = sleep_time
   if (threading.activeCount() == 1):
      gobject.timeout_add(50, python_thread_sleeper)
   MyThread().start()


def map_sharpening_gui(imol):

   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   vbox = gtk.VBox(False, 2)
   hbox = gtk.HBox(False, 2)
   adj = gtk.Adjustment(0.0, -30, 60, 0.05, 2, 30.1)
   slider = gtk.HScale(adj)
   label = gtk.Label("\nSharpen Map:")
   lab2  = gtk.Label("Add B-factor: ")

   vbox.pack_start(label,  False, False, 2)
   vbox.pack_start(hbox,   False, False, 2)
   hbox.pack_start(lab2,   False, False, 2)
   hbox.pack_start(slider, True,  True,  2)
   window.add(vbox)
   window.set_size_request(500, 100)
   # slider.add_mark(-30, -30, 0) # not yet, needs updated pygtk

   adj.connect("value_changed", lambda func: sharpen(imol, adj.value))
   
   window.show_all()


# Associate the contents of a sequence file with a chain.
# Select file from a GUI.
# File format can be Fasta (default) or PIR
#
# added checkbox to assign for all protein chains
#
def associate_sequence_with_chain_gui(sequence_format="FASTA",
                                      do_alignment=False):

   def associate_func(imol, chain_id_in, pir_file_name,
                      use_all_chains=False):
      #print "assoc seq:", imol, chain_id, pir_file_name
      all_chains = chain_id_in
      if use_all_chains:
         all_chains = [chain for chain in chain_ids(imol) if is_protein_chain_qm (imol, chain)]
      for chain_id in all_chains:
         if (sequence_format == "FASTA"):
            associate_fasta_file(imol, chain_id, pir_file_name)
         elif (sequence_format == "PIR"):
            associate_pir_file(imol, chain_id, pir_file_name)
         else:
            info_dialog("BL INFO:: wrong sequence input format.")
            return
      if do_alignment:
         alignment_mismatches_gui(imol)

   generic_chooser_entry_and_file_selector(
      "Associate Sequence with Chain: ",
      valid_model_molecule_qm,
      "Chain ID",
      "",
      "Select " + sequence_format +" file",
      lambda imol, chain_id, pir_file_name:
      associate_func(imol, chain_id, pir_file_name),
      use_check_button=True,
      check_button_label="Assign all protein chains (ignoring chain above)",
      alternative_callback_function=(lambda imol, chain_id, pir_file_name:
                                     associate_func(imol, chain_id, pir_file_name,
                                                    use_all_chains=True)))


# Make a box-of-buttons GUI for the various modifications that need
# to be made to match the model sequence to the assigned sequence(s).
# 
# Call this when the associated sequence(s) have been read in already.
#
def alignment_mismatches_gui(imol):

   # Return CA if there is such an atom in the residue, else return
   # the first atom in the residue.
   #
   def get_sensible_atom_name(res_info):
      chain_id = res_info[2]
      res_no   = res_info[3]
      ins_code = res_info[4]
      residue_atoms = residue_info(imol, chain_id, res_no, ins_code)
      if not residue_atoms:
         return " CA "  # wont work of course
      else:
         for atoms in residue_atoms:
            if (atoms[0][0] == " CA "):
               return " CA "
         return residue_atoms[0][0][0]
      
   # main line
   am = alignment_mismatches(imol)

   if (am == []):
      info_dialog("No sequence mismatches")
   else:
      if not am:
         info_dialog("Sequence not associated - no alignment")
      else:
         #print "mutations", am[0]
         #print "deletions", am[1]
         #print "insertions", am[2]

         # a button is a list [button_label, button_action]
         def mutate_buttons():
            ret_buttons = []
            for res_info in am[0]:
               chain_id = res_info[2]
               res_no   = res_info[3]
               ins_code = res_info[4]
               button_1_label = "Mutate " + chain_id + \
                                " " + str(res_no) + \
                                " " + residue_name(imol, chain_id, res_no, ins_code) + \
                                " to " + res_info[0]
               button_1_action = ["set_go_to_atom_molecule(" + str(imol) + ")",
                                  "set_go_to_atom_chain_residue_atom_name(\'" + \
                                  chain_id + "\', " + \
                                  str(res_no) + ", " + \
                                  "\'" + get_sensible_atom_name(res_info) + "\')"]
               ret_buttons.append([button_1_label, button_1_action])
            return ret_buttons

         def delete_buttons():
            ret_buttons = []
            for res_info in am[1]:
               chain_id = res_info[2]
               res_no   = res_info[3]
               ins_code = res_info[4]
               button_1_label = "Delete " + chain_id + \
                                " " + str(res_no)
               button_1_action = ["set_go_to_atom_molecule(" + str(imol) + ")",
                                  "set_go_to_atom_chain_residue_atom_name(\'" + \
                                  chain_id + "\', " + \
                                  str(res_no) + ", " + \
                                  "\'" + get_sensible_atom_name(res_info) + "\')"]
               ret_buttons.append([button_1_label, button_1_action])
            return ret_buttons

         def insert_buttons():
            ret_buttons = []
            for res_info in am[2]:
               chain_id = res_info[2]
               res_no   = res_info[3]
               ins_code = res_info[4]
               button_1_label = "Insert " + chain_id + \
                                " " + str(res_no)
               #button_1_action = "info_dialog(" + button_1_label + ")"
               # oh dear, will that work? I dont think I can pass asignments
               # here. Need to rethink this! FIXME
               # messy but should work (without too much hazzle)
               button_1_action = ["info_dialog(\'" + button_1_label + "\')",
                                  "set_go_to_atom_chain_residue_atom_name(\'" + \
                                  chain_id + "\', " + \
                                  "nearest_residue_by_sequence(" + \
                                  str(imol) + ", \'" + \
                                  chain_id + "\', " + \
                                  str(res_no) + ", " + \
                                  "\'" + ins_code + "\')[2], " + \
                                  "\' CA \')"]
               ret_buttons.append([button_1_label, button_1_action])
            return ret_buttons

         buttons  = delete_buttons()
         buttons += mutate_buttons()
         buttons += insert_buttons()

         dialog_box_of_buttons("Residue mismatches", [300, 300],
                               buttons, "  Close  ")

# Wrapper in that we test if there have been sequence(s) assigned to
# imol before we look for the sequence mismatches
#
def wrapper_alignment_mismatches_gui(imol):

   seq_info = sequence_info(imol)
   print("BL DEBUG:: sequence_info", seq_info)
   if seq_info:
      alignment_mismatches_gui(imol)
   else:
      associate_sequence_with_chain_gui(do_alignment=True)


# Multiple residue ranges gui
# 
# Create a new top level window that contains a residue-range-vbox
# which contains a set of hboxes that contain (or allow the user to
# enter) a residue range (chain-id resno-start resno-end).
# 
# The '+' and '-' buttons on the right allow the addition of extra
# residue ranges (and remove them of course).  The residue ranges are
# added to residue-range-widgets in a somewhat ugly manner.  Note
# also that fill-residue-range-widgets-previous-data generates
# residue range widgets and adds them to residue-range-widgets.
# 
# The interesting function is make-residue-range-frame which returns
# the outside vbox (that contains the frame and + and - buttons, so
# it is not a great name for the function) and each of the entries -
# so that they can be decoded (gtk-entry-get-text) when the "Go"
# button is pressed.
#
# Using the variable saved-residue-ranges, the GUI can restore itself
# from previous (saved) residue ranges.
#
# (Notice that we are not dealing with insertion codes).
#
global saved_residue_ranges
saved_residue_ranges = []
def residue_range_gui(func, function_text, go_button_label):

   global saved_residue_ranges
   residue_range_widgets = []
   #
   def n_residue_range_vboxes(residue_range_widgets_vbox):
      ls = residue_range_widgets_vbox.get_children()
      length = len(ls)
      return length

   # Remove widget from residue-range-widgets (on '-' button
   # pressed)
   def remove_from_residue_range_widget(widget):
      for ls in residue_range_widgets:
         if (widget == ls[0]):
                 rem = ls
                 break
      residue_range_widgets.remove(rem)

   #
   def make_residue_range_frame(residue_range_vbox):

      def plus_button_cb(*args):
         # we need to add a new residue-range
         # outside-hbox into the residue-range-widgets-vbox
         rr_frame = make_residue_range_frame(residue_range_vbox)
         residue_range_vbox.pack_start(rr_frame[0], False, False, 2)
         rr_frame[0].show()
         residue_range_widgets.append(rr_frame)

      def minus_button_cb(*args):
         n = n_residue_range_vboxes(residue_range_vbox)
         if (n > 1):
            remove_from_residue_range_widget(outside_hbox)
            outside_hbox.destroy()
      
      frame = gtk.Frame()
      outside_hbox = gtk.HBox(False, 2)
      hbox = gtk.HBox(False, 2)
      text_1 = gtk.Label("  Chain-ID:")
      text_2 = gtk.Label("  Resno Start:")
      text_3 = gtk.Label("  Resno End:")
      entry_1 = gtk.Entry()
      entry_2 = gtk.Entry()
      entry_3 = gtk.Entry()
      plus_button  = gtk.Button("+")
      minus_button = gtk.Button(" - ")

      hbox.pack_start(text_1,  False, False, 0)
      hbox.pack_start(entry_1, False, False, 0)
      hbox.pack_start(text_2,  False, False, 0)
      hbox.pack_start(entry_2, False, False, 0)
      hbox.pack_start(text_3,  False, False, 0)
      hbox.pack_start(entry_3, False, False, 0)

      outside_hbox.pack_start(frame, False, False, 2)
      frame.add(hbox)
      outside_hbox.pack_start(plus_button,  False, False, 2)
      outside_hbox.pack_start(minus_button, False, False, 2)

      plus_button.connect("clicked", plus_button_cb)
      
      minus_button.connect("clicked", minus_button_cb)

      list(map(lambda x: x.show(), [frame, outside_hbox, hbox, text_1, text_2, text_3,
                               entry_1, entry_2, entry_3, plus_button, minus_button]))

      # return the thing that we need to pack and the entries we
      # need to read.
      return [outside_hbox, entry_1, entry_2, entry_3]

   #
   def make_residue_ranges():
      ls = list(map(make_residue_range, residue_range_widgets))
      return ls

   # Return a list ["A", 2, 3] or False on failure to decose entries
   #
   def make_residue_range(residue_range_widget):
      # print "make a residue range using", residue_range_widget
      entry_1 = residue_range_widget[1]
      entry_2 = residue_range_widget[2]
      entry_3 = residue_range_widget[3]
      chain_id = entry_1.get_text()
      res_no_1_txt = entry_2.get_text()
      res_no_2_txt = entry_3.get_text()
      try:
         res_no_1 = int(res_no_1_txt)
         res_no_2 = int(res_no_2_txt)
         return [chain_id, res_no_1, res_no_2]
      except:
         print("did not understand %s and %s as numbers - fail residue range" %(res_no_1_txt, res_no_2_txt))
         return False

   #
   def save_ranges(residue_range_widgets):
      global saved_residue_ranges
      residue_ranges = []
      for residue_range_widget in residue_range_widgets:
         ls = make_residue_range(residue_range_widget)
         residue_ranges.append(ls)
      saved_residue_ranges = residue_ranges

   # range_info is list [chain_id, res_no_1, res_no_2]
   #
   def fill_with_previous_range(range_info, vbox_info):
      print("fill_with_previous_range using", range_info)
      entry_1 = vbox_info[1]
      entry_2 = vbox_info[2]
      entry_3 = vbox_info[3]
      chain_id = range_info[0]
      resno_1  = range_info[1]
      resno_2  = range_info[2]

      entry_1.set_text(chain_id)
      entry_2.set_text(str(resno_1))
      entry_3.set_text(str(resno_2))
      
   #
   def fill_residue_range_widgets_previous_data(previous_ranges,
                                                first_vbox_info,
                                                residue_range_vbox):
      print("first one", previous_ranges[0])
      if previous_ranges:
         if previous_ranges[0]:
            fill_with_previous_range(previous_ranges[0], first_vbox_info)
      if (len(previous_ranges) > 1):
         for rang in previous_ranges[1:]:
            if rang:
               vbox_info = make_residue_range_frame(residue_range_vbox)
               residue_range_vbox.pack_start(vbox_info[0], False, False, 2)
               print("next one", rang)
               residue_range_widgets.append(vbox_info)
               fill_with_previous_range(rang, vbox_info)

   # main line
   #
   def cancel_button_cb(*args):
      save_ranges(residue_range_widgets)
      window.destroy()
      return False

   def go_button_cb(*args):
      from types import IntType
      save_ranges(residue_range_widgets)
      residue_ranges = make_residue_ranges()
      imol = get_option_menu_active_molecule(*mc_opt_menu_model_list)
      if (type(imol) is IntType):
         func(imol, residue_ranges)
      else:
         print("BL INFO:: couldn't get valid imol!!")
      window.destroy()
      return False
      
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   vbox = gtk.VBox(False, 0)
   residue_range_vbox = gtk.VBox(True, 2)
   residue_range_widget_info = make_residue_range_frame(residue_range_vbox)
   hbox_buttons = gtk.HBox(False, 0)
   function_label = gtk.Label(function_text)
   cancel_button = gtk.Button("  Cancel  ")
   go_button = gtk.Button(go_button_label)
   h_sep = gtk.HSeparator()
   # the first residue range
   outside_vbox_residue_range = residue_range_widget_info[0]

   residue_range_widgets = [residue_range_widget_info] #?

   # buttons
   hbox_buttons.pack_start(cancel_button, False, False, 6)
   hbox_buttons.pack_start(    go_button, False, False, 6)

   # the vbox of residue ranges
   residue_range_vbox.pack_start(outside_vbox_residue_range, False, False, 0)

   if saved_residue_ranges:
      fill_residue_range_widgets_previous_data(saved_residue_ranges,
                                               residue_range_widget_info,
                                               residue_range_vbox)
      
   # main vbox
   vbox.pack_start(function_label, False, False, 0)
   mc_opt_menu_model_list = generic_molecule_chooser(vbox,
                                                     "Molecule for Ranges:")
   vbox.pack_start(residue_range_vbox, False, False, 2)
   vbox.pack_start(h_sep, True, True, 6)
   vbox.pack_start(hbox_buttons, False, False, 0)

   window.add(vbox)
   vbox.set_border_width(6)

   cancel_button.connect("clicked", cancel_button_cb)

   go_button.connect("clicked", go_button_cb)

   window.show_all()
            

global additional_solvent_ligands
additional_solvent_ligands = []

global random_jiggle_n_trials
random_jiggle_n_trials = 50

def solvent_ligand_list():
   global additional_solvent_ligands
   return (additional_solvent_ligands +
           ["EDO", "GOL", "DMS", "ACT", "MPD", "CIT", "SO4", "PO4", "TRS",
            "TAM", "PEG", "PG4", "PE8", "EBE", "BTB"])

# add solvent molecules
#
# Change the translation jiggle-factor to 1.0, so the ligand doesn't
# move so far and get sucked into protein density (this is just a
# temporary hack, it would be better to mask the enviroment of the
# ligand by the surrounding atoms of the molecule to which the ligand
# is added - that is much harder).
#
def solvent_ligands_gui():

   global random_jiggle_n_trials
   #
   def add_ligand_func(imol, tlc):
      print("Add a %s to molecule %s here" %(tlc, imol))
      imol_ligand = get_monomer(tlc)
      if (valid_model_molecule_qm(imol_ligand)):
         # delete hydrogens from the ligand if the master molecule
         # does not have hydrogens.
         if (valid_model_molecule_qm(imol)):
            if (not molecule_has_hydrogens(imol)):
               delete_residue_hydrogens(imol_ligand, "A", 1, "", "")
         if (valid_map_molecule_qm(imol_refinement_map())):
            print("========  jiggling!  ======== ")

            merge_molecules([imol_ligand], imol)

            # We no longer do this: because we now mask the map by the neighbours
            # of the jiggling ligand and for that we need to use molecule imol
            #
            # fit_to_map_by_random_jiggle(imol_ligand, "A", 1, "",
            #                             random_jiggle_n_trials, 1.0)
            # with_auto_accept([refine_zone, imol_ligand, "A", 1, 1, ""])


            # we presume that the active residue is the jiggling residue!
            # (we'd like to know what imol_ligand: "A" 1 moves to on merging
            # but we don't)
            #
            active_atom = active_residue()
            aa_chain_id = active_atom[1]
            aa_res_no   = active_atom[2]
            fit_to_map_by_random_jiggle(imol, aa_chain_id, aa_res_no, "",
                                        random_jiggle_n_trials, 1.0)

            # if we use refine_residues, that will take note of residues
            # near this residue and make non-bonded contacts
            # (whereas refine_zone will not).
            #
            # with_auto_accept([refine_zone, imol, aa_chain_id, aa_res_no, 1, ""])
            with_auto_accept([refine_residues, imol, [[aa_chain_id, aa_res_no, ""]]])

         else:
            print("======== not jiggling - no map ======== ")
         if valid_model_molecule_qm(imol):
            set_mol_active(imol_ligand, 0)
            set_mol_displayed(imol_ligand, 0)
   
   # add a button for a 3-letter-code to the scrolled vbox that runs
   # add-ligand-func when clicked.
   #
   def add_solvent_button(comp_id, button_label, inside_vbox,
                          molecule_option_menu, model_list):
      def button_cb(*args):
         imol = get_option_menu_active_molecule(molecule_option_menu,
                                                model_list)
         add_ligand_func(imol, comp_id)

      button = gtk.Button(button_label)
      inside_vbox.pack_start(button, False, False, 1)
      button.show()
      button.connect("clicked", button_cb)

   def delete_event(*args):
      window.destroy()
      return False

   def add_new_button_cb(*args):
      global additional_solvent_ligands
      def add_button_func(txt):
         additional_solvent_ligands.append(txt)
         add_solvent_button(txt,
                            comp_id2button_label(txt), inside_vbox,
                            molecule_option_menu,
                            model_list)
      generic_single_entry("Add new 3-letter-code/comp-id",
                           "", "  Add  ",
                           lambda txt:
                              add_button_func(txt))

   def comp_id2button_label(comp_id):
      auto_load_dictionary(comp_id)
      comp_id_name = comp_id2name(comp_id)
      ret = (comp_id + ": " + comp_id_name) if comp_id_name else comp_id
      return ret
      
   # main
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   scrolled_win = gtk.ScrolledWindow()
   outside_vbox = gtk.VBox(False, 2)
   inside_vbox  = gtk.VBox(False, 2)
   label = gtk.Label("\nSolvent molecules added to molecule: ")
   frame_for_option_menu = gtk.Frame(" Choose Molecule ")
   vbox_for_option_menu = gtk.VBox(False, 6)
   molecule_option_menu = gtk.combo_box_new_text()
   model_list = fill_option_menu_with_coordinates_mol_options(molecule_option_menu)
   add_new_button = gtk.Button("  Add a new Residue Type...")
   h_sep = gtk.HSeparator()
   close_button = gtk.Button("  Close  ")

   window.set_default_size(250, 500)
   window.set_title("Solvent Ligands")
   window.set_border_width(8)
   window.add(outside_vbox)
   outside_vbox.pack_start(label, False, False, 2)
   frame_for_option_menu.add(vbox_for_option_menu)
   vbox_for_option_menu.pack_start(molecule_option_menu, False, False, 8)
   frame_for_option_menu.set_border_width(6)
   outside_vbox.pack_start(frame_for_option_menu, False, False, 2)
   outside_vbox.pack_start(scrolled_win, True, True, 0)
   scrolled_win.add_with_viewport(inside_vbox)
   scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)
   outside_vbox.pack_start(add_new_button, False, False, 6)
   outside_vbox.pack_start(h_sep, False, False, 2)
   outside_vbox.pack_start(close_button, False, False, 2)

   for comp_id in solvent_ligand_list():
      button_label = comp_id2button_label(comp_id)
      add_solvent_button(comp_id, button_label, inside_vbox,
                         molecule_option_menu, model_list)

   add_new_button.connect("clicked", add_new_button_cb)
   close_button.connect("clicked", delete_event)

   window.show_all()


# USER MODS gui
#
def user_mods_gui(imol, pdb_file_name):

   # no alt conf, no inscode
   def atom_spec_to_string(atom_spec):
      chain_id  = atom_spec[1]
      res_no    = atom_spec[2]
      ins_code  = atom_spec[3]
      atom_name = atom_spec[4]
      alt_conf  = atom_spec[5]
      return chain_id + str(res_no) + atom_name

   #
   def make_flip_buttons(flip_list):
      ret = []
      for flip in flip_list:
         atom_spec    = flip[0]
         residue_type = flip[1]
         info_string  = flip[2]
         set_string   = flip[3]
         score        = flip[4]
         chain_id  = atom_spec[1]
         res_no    = atom_spec[2]
         ins_code  = atom_spec[3]
         atom_name = atom_spec[4]
         alt_conf  = atom_spec[5]
         label = set_string + " " + \
                 chain_id + " " + str(res_no) + \
                 atom_name + " : " + \
                 info_string + " " + \
                 " score %2.2f" % score
         func = [cmd2str(set_go_to_atom_molecule, imol),
                 cmd2str(set_go_to_atom_chain_residue_atom_name,
                         chain_id, res_no, atom_name)]
         ret.append([label, func])
      return ret

   #
   def make_no_adj_buttons(no_flip_list):
      ret = []
      for no_flip_item in no_flip_list:
         specs = no_flip_item[0]
         info_string = no_flip_item[1]
         label = "No Adjustment " + \
                 " ".join(map(atom_spec_to_string, specs)) + \
                 " " + \
                 info_string
         atom_spec = specs[0]
         chain_id  = atom_spec[1]
         res_no    = atom_spec[2]
         ins_code  = atom_spec[3]
         atom_name = atom_spec[4]
         alt_conf  = atom_spec[5]
         func = [cmd2str(set_go_to_atom_molecule, imol),
                 cmd2str(set_go_to_atom_chain_residue_atom_name,
                         chain_id, res_no, atom_name)]
         ret.append([label, func])
      return ret
   
   # Return a list of buttons that are (in this case, (only) clashes,
   # unknown/problems and flips.
   #
   def filter_buttons(flip_list):
      ret = []
      for flip in flip_list:
         atom_spec    = flip[0]
         residue_type = flip[1]
         info_string  = flip[2]
         set_string   = flip[3]
         score        = flip[4]
         if len(info_string) < 3:
            pass
         else:
            # keep-letter: K is for keep, 
            # C is clash, X is "not sure"
            # F is flip.
            keep_letter = info_string[2:3]
            if keep_letter in "XFC":
               ret.append(flip)
      return ret

   # filter flag is True or False
   def clear_and_add_back(vbox, flip_list, no_adj_list, filter_flag):
      # clear
      children = vbox.get_children()
      list(map(lambda c: c.destroy(), children))
      # add back
      if not filter_flag:
         buttons = make_no_adj_buttons(no_adj_list) + \
                  make_flip_buttons(flip_list)
      else:
         # filter
         buttons = make_flip_buttons(filter_buttons(flip_list))
      list(map(lambda button_info: add_button_info_to_box_of_buttons_vbox(button_info, vbox),
          buttons))

   # main line
   #
   # user mods will return a pair of lists.
   if using_gui():
      flips = user_mods(pdb_file_name)
      flip_buttons = make_flip_buttons(flips[0])
      no_adj_buttons = make_no_adj_buttons(flips[1])
      all_buttons = no_adj_buttons + flip_buttons
      dialog_box_of_buttons_with_check_button(
         " Flips ", [300, 300], all_buttons, "  Close  ",
         "Clashes, Problems and Flips only",
         lambda check_button, vbox: clear_and_add_back(vbox, flips[0], flips[1], True)
                                    if check_button.get_active() else
                                    clear_and_add_back(vbox, flips[0], flips[1], False),
         False)

def rename_residue_gui_simple():
   active_atom = active_residue()
   if (not active_atom):
      info_dialog("No Residue Here")
   else:
      print(active_atom)
      generic_single_entry("Rename this residue", "ALA", "Rename",
                           lambda text: using_active_atom(set_residue_name,
                                                          "aa_imol", "aa_chain_id", "aa_res_no", "aa_ins_code",
                                                            text))
#                           lambda text: using_active_atom([[set_residue_name,
#                                                            ["aa_imol", "aa_chain_id", "aa_res_no", "aa_ins_code"],
#                                                            [text]]]))

def average_map_gui():
   
   # mav-widgets (mav is "map average") is something like a
   # (local-to-the-function) "static" - it can be refered to in the
   # call-backs.  Presumably, if you refactor this as a class, it can
   # be a class data item.
   #
   # BL syas:  ok, lets make a class, just for the sake of it
   class mav:
      def __init__(self):
         self.mav_widgets = []
      
      # On pressing the - button, delete the widget from the
      # mav-widgets store.  The widget is actually gtk-widget-deleted
      # elsewhere.
      #
      def remove_from_mav_widgets(self, widget):
         to_be_removed = False
         for item in self.mav_widgets:
            if (item[0] == widget):
               to_be_removed = item
               break
         if to_be_removed:
            self.mav_widgets.remove(to_be_removed)

      # Return a list of the hbox the optionmenu the model-mol-list
      # and the entry
      def add_average_molecule_widget(self, maps_vbox):

         def plus_button_cb(*args):
            mav_bits = self.add_average_molecule_widget(maps_vbox)
            mav_bits[0].show()
            self.mav_widgets.append(mav_bits)

         def minus_button_cb(*args):
            if len(self.mav_widgets) > 1:
               self.remove_from_mav_widgets(hbox)
               hbox.destroy()
            
         frame = gtk.Frame(False)
         hbox = gtk.HBox(False, 2)
         label = gtk.Label("Weight:  ")
         entry = gtk.Entry()
         optionmenu = gtk.combo_box_new_text()
         map_mol_list = fill_option_menu_with_map_mol_options(optionmenu)
         plus_button  = gtk.Button("+")
         minus_button = gtk.Button(" - ")
         hbox.pack_start(optionmenu, False, False, 2)
         hbox.pack_start(label, False, False, 2)
         hbox.pack_start(entry, False, False, 2)
         hbox.pack_start( plus_button, False, False, 2)
         hbox.pack_start(minus_button, False, False, 2)
         entry.set_size_request(40, -1)
         entry.set_text("1.0")

         plus_button.connect("clicked", plus_button_cb)

         # when the - button is clicked, delete the hbox and
         # everything in it.  Except if it was the only hbox/line,
         # and in that case, don't do anything.
         minus_button.connect("clicked", minus_button_cb)

         maps_vbox.pack_start(hbox, False, False, 2)
         # show everything we just created
         list(map(lambda x: x.show(), [frame, hbox, label, entry, optionmenu,
                                  plus_button, minus_button]))

         # print "saving map_mol_list", map_mol_list
         return [hbox, optionmenu, map_mol_list, entry]

      # main line
      #
      # create the usual outside vbox for everything, and the
      # inner-vbox which is for the map hboxes.
      #
      def main_line(self):
         
         def cancel_button_cb(*args):
            window.destroy()
            return False
         
         def ok_button_cb(*args):
            maps_to_average_list = []
            for mav_bits in self.mav_widgets:
               option_menu  = mav_bits[1]
               map_mol_list = mav_bits[2]
               entry        = mav_bits[3]
               print("map_mol_list", map_mol_list)
               map_selected = get_option_menu_active_molecule(option_menu, map_mol_list)
               text = entry.get_text()
               weight = float(text)   # try?! FIXME
               maps_to_average_list.append([map_selected, weight])
            print("maps to average", maps_to_average_list)
            average_map(maps_to_average_list)
            window.destroy()
            return False
            
         window = gtk.Window(gtk.WINDOW_TOPLEVEL)
         outer_vbox = gtk.VBox(False, 0)
         inner_vbox = gtk.VBox(False, 0)
         title = gtk.Label("Average Maps")
         h_sep = gtk.HSeparator()
         buttons_hbox = gtk.HBox(False, 2)
         mav_widget = self.add_average_molecule_widget(inner_vbox)
         cancel_button = gtk.Button("  Cancel  ")
         ok_button = gtk.Button("  Average Maps  ")

         window.add(outer_vbox)
         outer_vbox.pack_start(title, False, False, 2)
         outer_vbox.pack_start(inner_vbox, False, False, 2)
         outer_vbox.pack_start(h_sep, False, False, 2)
         outer_vbox.pack_start(buttons_hbox, False, False, 6)
         buttons_hbox.pack_end(ok_button, False, False, 6)
         buttons_hbox.pack_end(cancel_button, False, False, 6)

         # reset mav-widget when we start a new dialog (otherwise, it
         # retains the old mav-widgets, which is confusing/wrong.
         #
         self.mav_widgets = [mav_widget]

         cancel_button.connect("clicked", cancel_button_cb)

         # On clicking "Average Maps", we get a list of maps to
         # average by looking at the widgets (the option menu and the
         # entry) for each mav-bits (i.e. each widget list that is
         # returned by add-average-molecule-widget).
         #
         ok_button.connect("clicked", ok_button_cb)
         window.show_all()
         
   gui = mav()
   gui.main_line()


# simple rename residue GUI
#
def rename_residue_gui():
   active_atom = active_residue()
   if not active_atom:
      info_dialog("No Residue Here")
   else:
      aa_imol     = active_atom[0]
      aa_chain_id = active_atom[1]
      aa_res_no   = active_atom[2]
      aa_ins_code = active_atom[3]
      label = "Rename Residue [in molecule " + \
              str(aa_imol) + \
              "]: " + \
              aa_chain_id + \
              str(aa_res_no) + \
              aa_ins_code + \
              " to: "
      generic_single_entry(label, "ALA", "Rename Residue",
                           lambda text: set_residue_name(aa_imol, aa_chain_id,
                                                         aa_res_no, aa_ins_code,
                                                         text))


# molecule chooser
# 
# coordination number chooser
# 
# Max dist entry
# 
# results vbox
#    containing buttons with atom-spec labels
# 
# h-sep
#
#      close-button
#
def water_coordination_gui():

   def delete_event(*args):
      window.destroy()
      return False

   def apply_cb(*args):
      n = get_number()
      imol = get_molecule()
      d = get_distance()
      if d:
         update_water_results(imol, n, d)
      
   def key_press_event(widget, event):
      if (event.keyval == 65293):  # GDK_return
         n = get_number()
         imol = get_molecule()
         d = get_distance()
         if d:
            update_water_results(imol, n, d)
      
   def atom_spec_to_text(atom_spec):
      return " ".join(map(str, atom_spec))

   def get_ele(imol, atom_spec):
      from types import ListType
      atoms = residue_info(imol,
                           atom_spec[1],
                           atom_spec[2],
                           atom_spec[3])
      input_atom_name = atom_spec[4]
      if type(atoms) is ListType:
         for atom in atoms:
            if atom[0][0] == input_atom_name:
               return atom[1][2]
      return False
   
   # add info about bumps (to non-H-bonding atoms or so).  given a
   # water info (a central atom spec and a list of its contacts).
   def make_bump_text(imol, water_info):
      central_atom = water_info[0]
      contact_list = water_info[1]
      rv = ""
      for contact in contact_list:
         ele = get_ele(imol, contact)
         if (ele == " C"):
            rv = " Bump"
      return rv

   def is_a_metal_site_too_qm(atom_spec, metal_results):
      # print "   checking if %s is in %s... result: %s" %(atom_spec, metal_results, atom_spec in metal_results)
      return atom_spec in metal_results

   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   vbox = gtk.VBox(False, 0)
   results_vbox = gtk.VBox(False, 0)
   water_results_label = gtk.Label("Other Coordinated Waters")
   metal_results_vbox = gtk.VBox(False, 0)
   metal_results_label = gtk.Label("Potential Metals: ")
   h_sep = gtk.HSeparator()
   hint_text = "Molecule: "
   hbox_chooser = gtk.HBox(False, 0)
   hbox_max_dist = gtk.HBox(False, 0)
   hbox_number_chooser = gtk.HBox(False, 0)
   number_text = gtk.Label("Coordination Number: ")
   molecule_chooser_option_menu_and_model_list = generic_molecule_chooser(hbox_chooser, hint_text)
   molecule_chooser_option_menu = molecule_chooser_option_menu_and_model_list[0]
   model_list = molecule_chooser_option_menu_and_model_list[1]
   scrolled_win = gtk.ScrolledWindow()
   metal_results_scrolled_win = gtk.ScrolledWindow()
   number_menu = gtk.combo_box_new_text()
   number_list = list(range(3, 10))
   dist_label = gtk.Label("Max Dist: ")
   dist_entry = gtk.Entry()
   close_button = gtk.Button("  Close  ")
   apply_button = gtk.Button("  Apply  ")
   hbox_buttons = gtk.HBox(False, 6)

   def get_molecule():
      return get_option_menu_active_molecule(molecule_chooser_option_menu,
                                             model_list)
   def get_number():
      return int(get_option_menu_active_item(number_menu,
                                             number_list))
   def get_distance():
      t = dist_entry.get_text()
      try:
         ret = float(t)
      except:
         ret = False
      return ret

   def clear_previous_results():
      for this_vbox in [results_vbox, metal_results_vbox]:
         children = this_vbox.get_children()
         res = [widget.destroy() for widget in children]
         #for child in children:
         #   print "BL DEBUG:: now destroy child", child
         #   child.destroy()
         #this_vbox.hide()
         #this_vbox.show()

   def update_water_results(imol, n, d):
      results = highly_coordinated_waters(imol, n, d)
      clear_previous_results()
      import time
      time.sleep(2)
      if results:
         coordination_results = results[1]
         metal_results = results[0]
         for water_info in coordination_results:
            # a water spec is a central-atom spec and a list
            # of its neighbours.
            atom_spec = water_info[0]
            t = atom_spec_to_text(atom_spec)
            bump_text = make_bump_text(imol, water_info)
            button = gtk.Button(t + bump_text)
            if not is_a_metal_site_too_qm(atom_spec, metal_results):
               results_vbox.pack_start(button, False, False, 1)
               button.show()
               def water_func(widget, imol, water_info):
                  water_spec = water_info[0]
                  chain_id   = water_spec[1]
                  res_no     = water_spec[2]
                  atom_name  = water_spec[4]
                  set_go_to_atom_molecule(imol)
                  set_go_to_atom_chain_residue_atom_name(chain_id, res_no, atom_name)
               button.connect("clicked", water_func, imol, water_info)
         
         # now handle metal results
         for metal_site in metal_results:
            metal_text = metal_site[1]
            t = atom_spec_to_text(metal_site[0])
            button_text = t + " Potential " + metal_text
            button = gtk.Button(button_text)
            metal_results_vbox.pack_start(button, False, False, 1)
            button.show()
            def metal_func(widget, imol, metal_site):
               metal_spec = metal_site[0]
               chain_id   = metal_spec[1]
               res_no     = metal_spec[2]
               atom_name  = metal_spec[4]
               set_go_to_atom_molecule(imol)
               set_go_to_atom_chain_residue_atom_name(chain_id, res_no, atom_name)
            button.connect("clicked", metal_func, imol, metal_site)

   window.add(vbox)
   
   fill_option_menu_with_number_options(number_menu, number_list, 5)
   
   scrolled_win.add_with_viewport(results_vbox)
   scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)
   
   metal_results_scrolled_win.add_with_viewport(metal_results_vbox)
   metal_results_scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)

   hbox_max_dist.pack_start(dist_label, False, False, 2)
   hbox_max_dist.pack_start(dist_entry, False, False, 2)

   vbox.pack_start(hbox_chooser, False, False, 6)

   hbox_number_chooser.pack_start(number_text, False, False, 2)
   hbox_number_chooser.pack_start(number_menu, False, False, 2)

   vbox.pack_start(hbox_number_chooser, False, False, 6)

   vbox.pack_start(hbox_max_dist, False, False, 2)

   # metal sites
   vbox.pack_start(metal_results_label, False, False, 2)
   vbox.pack_start(metal_results_scrolled_win, True, True, 0)

   # interesting water sites
   vbox.pack_start(water_results_label, False, False, 2)
   vbox.pack_start(scrolled_win, True, True, 0)  # expand fill padding
   vbox.pack_start(h_sep, False, False, 2)
   hbox_buttons.pack_start(close_button, False, False, 2)
   hbox_buttons.pack_start(apply_button, False, False, 2)
   vbox.pack_start(hbox_buttons, False, False, 2)

   # From the Nayal and Di Cera (1996) paper, it seems that 2.7
   # and at least 4 oxygens is a good test for Na+ or other
   # interesting things
   #
   dist_entry.set_text("2.7")

   close_button.connect("clicked", delete_event)

   dist_entry.connect("key_press_event", key_press_event)

   apply_button.connect("clicked", apply_cb)

   window.show_all()

# return a list, or False (e.g. if not in same chain and molecule)
#
def min_max_residues_from_atom_specs(specs):

   print(specs)

   min_res_no = False
   max_res_no = False
   chain_id = False

   for spec in specs:
      spec_model = spec[1]
      spec_chain = spec[2]
      res_no     = spec[3]

      if isinstance(chain_id, str):
         if (chain_id != spec_chain):
            return False  # chain mismatch
      else:
         chain_id = spec_chain

      if not min_res_no:  # maybe check for number as well
         min_res_no = res_no
      else:
         if (res_no < min_res_no):
            min_res_no = res_no
      if not max_res_no:  # maybe check for number as well
         max_res_no = res_no
      else:
         if (res_no > max_res_no):
            max_res_no = res_no

   if (min_res_no and max_res_no and chain_id):
      return [min_res_no, max_res_no, chain_id]
   else:
      return False

# by default, rename loop residues to UNK.  If python True, then
# leave them as the residue names found in the database.
global db_loop_preserve_residue_names
db_loop_preserve_residue_names = False
               
def click_protein_db_loop_gui():

   global db_loop_preserve_residue_names
   def pick_loop_func(n):
      def pick_func(*atom_specs):
         residue_specs = list(map(atom_spec_to_residue_spec, atom_specs))
         imol = atom_specs[0][1]
         min_max_and_chain_id = min_max_residues_from_atom_specs(atom_specs)

         if not isinstance(min_max_and_chain_id, list):
            info_dialog("Picked atoms not in same molecule and chain")
         else:
            loop_mols = protein_db_loops(imol, residue_specs,
                                         imol_refinement_map(),
                                         10, db_loop_preserve_residue_names)
            imol_loop_orig = loop_mols[0][0]
            imol_loops_consolidated = loop_mols[0][1]
            loop_mols = loop_mols[1]
            min_resno = min_max_and_chain_id[0]
            max_resno = min_max_and_chain_id[1]
            chain_id  = min_max_and_chain_id[2]
            set_mol_active(imol_loops_consolidated, 1)

            if valid_model_molecule_qm(imol_loop_orig):
               if len(loop_mols) > 0:
                  buttons = [[str(loop_mol) + " " + molecule_name(loop_mol),
                                 lambda func: copy_residue_range(imol, chain_id,
                                                                 loop_mol, chain_id,
                                                                 min_resno, max_resno)] for loop_mol in loop_mols]
                  def toggle_func(imol):
                     toggle_active_mol(imol)
                     toggle_display_mol(imol)

                  loop_mols.append(imol_loops_consolidated)

                  dialog_box_of_buttons("Loop Candidates",
                                        [360, 200],
                                        [["Original loop", lambda func:
                                          copy_residue_range(imol, chain_id,
                                                             imol_loop_orig, chain_id,
                                                             min_resno, max_resno)
                                         ],
                                         ["Toggle All Candidate Loops", lambda func:
                                          toggle_func(imol_loops_consolidated)]
                                         ] + buttons,
                                        " Close ",
                                        lambda: [(set_mol_displayed(im, 0), set_mol_active(im,0)) for im in loop_mols])
                 
      user_defined_click(n, pick_func)
      
   generic_number_chooser(list(range(2,10)), 4,
                          "Number of residues for basis",
                          "Pick Atoms...",

                          lambda n: pick_loop_func(n))


def refmac_multi_sharpen_gui():

   def delete_event(*args):
      window.destroy()
      return False

   def sharpen_cb(widget, *args):

      # get max_band n_levels and map file name
      max_b = int(get_option_menu_active_item(option_menu_b_factor,
                                              b_factor_list))
      n_levels = int(get_option_menu_active_item(option_menu_n_levels,
                                                 n_levels_list))
      active_item_imol = get_option_menu_active_molecule(option_menu_map,
                                                         map_molecule_list)
      # There is no function to get a map file name from a molecule
      # It is not stored. So we make/guess it...
      map_file_name = molecule_name(active_item_imol)
      if (map_file_name.find(" ") > 0):
         # we have map coeffs - but then sharpen as here wont work anyway
         map_file_name = map_file_name[:map_file_name.find(" ")]
      map_file_name_stub = strip_path(file_name_sans_extension(map_file_name))
      refmac_output_mtz_file_name = "starting_map-" + map_file_name_stub + ".mtz"
      log_file_name = "refmac-multisharp-" + map_file_name_stub + ".log"
      if not os.path.isfile(map_file_name):
         info_dialog("WARNING:: file not found %s" %map_file_name)
      else:
         if not directory_is_modifiable_qm(os.getcwd()):
            m = "ERROR:: Current directory " + os.getcwd() + " is not writable"
            info_dialog(m)
            return
         print("active_item_imol", active_item_imol)
         step_size = max_b/n_levels
         numbers_string = ' '.join(str(i+1) for i in range(n_levels))
         blur_string = "SFCALC BLUR  " + numbers_string
         sharp_string = "SFCALC SHARP " + numbers_string

         cmd_line_args = ["MAPIN", map_file_name]
         data_lines = ["MODE SFCALC",
                       blur_string,
                       sharp_string,
                       "END"]
         this_dir = os.getcwd()
         if not directory_is_modifiable_qm(this_dir):
            info_dialog("WARNING:: Current directory is not writable")
         else:
            refmac_execfile = find_exe("refmac5", "CBIN", "CCP4_BIN", "PATH")
            s = popen_command(refmac_execfile,
                              cmd_line_args,
                              data_lines,
                              log_file_name,
                              False)

            try:
               if s != 0:
                  info_dialog("WARNING:: refmac5 failed")
               else:
                  # Happy path
                  print("BL DEBUG:: s", s)
                  if os.path.isfile("starting_map.mtz"):
                     os.rename("starting_map.mtz", refmac_output_mtz_file_name)
                     # offer a read-mtz dialog
                     manage_column_selector(refmac_output_mtz_file_name)

            except:
               print("BL DEBUG:: tried to rename starting-map.mtz but failed.")
               pass
         delete_event(widget)

   print("BL DEBUG:: now make a windwo")
   window = gtk.Window(gtk.WINDOW_TOPLEVEL)
   # boxes
   vbox = gtk.VBox(False, 0)
   hbox_1 = gtk.HBox(False, 0)
   hbox_2 = gtk.HBox(False, 0)
   hbox_3 = gtk.HBox(False, 0)
   # menus
   option_menu_map = gtk.combo_box_new_text()
   option_menu_b_factor = gtk.combo_box_new_text()
   option_menu_n_levels = gtk.combo_box_new_text()
   #labels
   map_label = gtk.Label("Map ")
   sb_label = gtk.Label("Sharpen & Blur in ")
   levels_label = gtk.Label(" levels up to ")
   A_label = gtk.Label(" A*A")
   # separate
   h_sep = gtk.HSeparator()
   # buttons
   ok_button = gtk.Button("   OK   ")
   cancel_button = gtk.Button(" Cancel ")
   n_levels_list = [1, 2, 3, 4, 5, 6]
   b_factor_list = [50, 100, 200, 400, 800, 2000]

   map_molecule_list = fill_option_menu_with_map_mol_options(option_menu_map)
   fill_option_menu_with_number_options(option_menu_n_levels, n_levels_list, 4)
   fill_option_menu_with_number_options(option_menu_b_factor, b_factor_list, 200)

   window.set_title("Refmac for Sharpening & Blurring")
   hbox_1.pack_start(map_label, False, False, 2)
   hbox_1.pack_start(option_menu_map, False, False, 2)
   hbox_2.pack_start(sb_label, False, False, 2)
   hbox_2.pack_start(option_menu_n_levels, False, False, 2)
   hbox_2.pack_start(levels_label, False, False, 2)
   hbox_3.pack_end(cancel_button, False, False, 12)
   hbox_3.pack_end(ok_button, False, False, 12)
   
   vbox.pack_start(hbox_1)
   vbox.pack_start(hbox_2)
   vbox.pack_start(h_sep)
   vbox.pack_start(hbox_3)
   
   cancel_button.connect("clicked", delete_event)

   ok_button.connect("clicked", sharpen_cb, option_menu_b_factor, b_factor_list,
                     option_menu_n_levels, n_levels_list,
                     option_menu_map, map_molecule_list)

   window.add(vbox)
   window.show_all()


def add_module_cryo_em():
   if coot_python.main_menubar():
      add_module_cryo_em_gui()

def add_module_ccp4():
   if coot_python.main_menubar():
      add_module_ccp4_gui()

def add_module_cryo_em_gui():
   if coot_python.main_menubar():
      menu = coot_menubar_menu("Cryo-EM")

      add_simple_coot_menu_menuitem(menu, "Sharpen/Blur...",
                                    lambda func: sharpen_blur_map_gui())

      add_simple_coot_menu_menuitem(menu, "Multi-sharpen...",
                                    lambda func: refmac_multi_sharpen_gui())

      def interactive_nudge_func():
         with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                        aa_ins_code, aa_atom_name,
                                        aa_alt_conf, aa_res_spec]:
            nudge_residues_gui(aa_imol, aa_res_spec)
      add_simple_coot_menu_menuitem(menu, "Interactive Nudge Residues...",
                                    lambda func: interactive_nudge_func())


def add_module_ccp4_gui():
   if coot_python.main_menubar():
      menu = coot_menubar_menu("CCP4")

      add_simple_coot_menu_menuitem(menu, "Make LINK via Acedrg",
                                    lambda func: acedrg_link_generation_control_window())

   
#### BL stuff
   
def scale_alt_conf_occ_gui(imol, chain_id, res_no, ins_code):

    alt_confs = residue_alt_confs(imol, chain_id, res_no, ins_code)
    try:
        # remove no alt conf, i.e. '', from list if present
        alt_confs.remove('')
    except:
        pass

    if (len(alt_confs) != 2):
        print("INFO:: no (or too many) alt confs, no gui!")
    else:
        # only do if 2 alt confs (plus no)
        res_info = residue_info(imol, chain_id, res_no, ins_code)
        alt_conf = alt_confs[0]

        def get_occ_for_alt_conf(atom_ls, alt_conf):
            # should check that they are consistent!!
            for i in range(len(atom_ls)):
                alt_conf_str = atom_ls[i][0][1]
                if (alt_conf_str == alt_conf):
                    occ = atom_ls[i][1][0]
                    return occ            

        occ_start = get_occ_for_alt_conf(res_info, alt_conf)
    
        def delete_event(*args):
            window.destroy()
            return False

        def go_function_event(widget, occ_adjust,
                              imol, chain_id, res_no, ins_code,
                              alt_confs):
            new_occ = occ_adjust.value
            print("BL INFO:: setting new occupancy ", new_occ)
            alt_conf_list = [[alt_confs[0], new_occ], [alt_confs[1], 1 - new_occ]]
            set_alt_conf_occ(imol, chain_id, res_no, ins_code, alt_conf_list)
            delete_event()

        def change_occ(*args):
            # needed?!
            # print "change occ slider?!"
            pass

        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        title = gtk.Label("Adjust alt conf occupancies")
        occ_label = gtk.Label("Occupancy")
        alt_conf_label = gtk.Label("Alt Conf: " + alt_conf)
        occ_adj = gtk.Adjustment(occ_start, 0.1, 0.99, 0.01, 0.1, 0.1)
        occ_scale = gtk.HScale(occ_adj)
        vbox = gtk.VBox(False, 0)
        scale_hbox = gtk.HBox(False, 0)
        scale_hbox.pack_start(alt_conf_label, False, False, 2)
        scale_hbox.pack_start(occ_scale, True, True, 2)
        vbox.pack_start(occ_label, False, False, 0)
        vbox.pack_start(scale_hbox, False, False, 0)

        occ_adj.connect("value_changed", change_occ)

        window.add(vbox)

        h_sep = gtk.HSeparator()
        buttons_hbox = gtk.HBox(True, 6)
        ok_button = gtk.Button("   OK   ")
        cancel_button = gtk.Button(" Cancel ")

        buttons_hbox.pack_start(ok_button, True, False, 2)
        buttons_hbox.pack_start(cancel_button, True, False, 2)

        vbox.pack_start(h_sep, True, False, 2)
        vbox.pack_start(buttons_hbox, True, False, 2)

        ok_button.connect("clicked", go_function_event, occ_adj,
                          imol, chain_id, res_no, ins_code, alt_confs)

        cancel_button.connect("clicked", delete_event)

        window.show_all()


def select_atom_alt_conf_occ_gui():

    def helper_function(*args):
        imol      = args[0][1]
        chain_id  = args[0][2]
        res_no    = args[0][3]
        ins_code  = args[0][4]
        scale_alt_conf_occ_gui(imol, chain_id, res_no, ins_code)
    user_defined_click(1, helper_function)


def toggle_backrub_rotamers(widget=None):
   """Toggle function to swtich on and off backrub rotamer fitting.
   
   Keyword arguments:
   widget -- can be passed from the toolbutton

   """

   if widget:
      if widget.get_active():
         # the button is toggled on
         set_rotamer_search_mode(ROTAMERSEARCHLOWRES)
         print("BL INFO:: Using Backrub rotamers now!")
      else:
         set_rotamer_search_mode(ROTAMERSEARCHHIGHRES)
         print("BL INFO:: NOT using Backrub rotamers any more!")

   else:
      # non graphical - but wont be able to run if this is not loaded.
      mode = rotamer_search_mode_state()
      if (mode == ROTAMERSEARCHLOWRES):
         set_rotamer_search_mode(ROTAMERSEARCHHIGHRES)
         print("BL INFO:: NOT using Backrub rotamers any more!")
      if (mode == ROTAMERSEARCHHIGHRES or
          mode == ROTAMERSEARCHAUTOMATIC):
         set_rotamer_search_mode(ROTAMERSEARCHLOWRES)
         print("BL INFO:: Using Backrub rotamers now!")
         
         
      # no alternative for now
      # need to be able to get the state of search mode.
      # easily added. FIXME
      print("BL WARNING:: no widget")
   
def toggle_hydrogen_display(widget=None):
      """Toggle function to display all hydrogens or not.

      Keyword arguments:
      widget -- can be passed from the toolbutton

      """

      if widget:
         if widget.get_active():
            # the button is toggled on
            hide_all_hydrogens()
         else:
            show_all_hydrogens()

      else:
         # non graphical - but wont be able to run if this is not loaded.
         print("BL INFO:: No display, so I dont care about the hydrogens.")
         print("BL WARNING:: no widget")


def toggle_wiimote(widget=None):
   """a toggle function to connect and disconnect from a Wiimote

   Keyword arguments:
   widget -- can be passed from the toolbutton

   """

   if widget:
      if widget.get_active():
         # the button is toggled on
         try:
            setup_wii()
         except NameError:
            print ("BL WARNING:: setup_wii not defined! " 
                   "Did you compile with WII_INTERFACE_WIIUSE?")
         except:
            print("BL WARNING:: could not set up Wii")
      else:
         try:
            stop_wii()
         except NameError:
            print ("BL WARNING:: stop_wii not defined! "
                   "Did you compile with WII_INTERFACE_WIIUSE?")
         except:
            print("BL WARNING:: could not stop wii")

   else:
      # no alternative for now (could just go by state and change back and forth)
      print("BL WARNING:: no widget")

# Simple minded dialog, search disk or not?
# Could/should be expaded to browse for exe and to select which 
# disks to search for.
#
def search_disk_dialog(program_name, path_ls):

   # graphics
   ret = False
   label_text = "Couldn't find %s in default path" %(program_name)
   for path in path_ls:
      label_text += " and "
      label_text += path
   label_text += "\n\nShall we search the whole disk?\n"

   try:
      ret = yes_no_dialog(label_text, "Search whole disk dialog")
   except:
      # no graphics
      label_text += "[y/N] >"
      result =""
      while result.lower() not in ['y', 'yes', 'n', 'no']:
         result = input(label_text)
      if result.lower() in ['y', 'yes']:
         ret = True
      else:
         ret = False

   return ret

def duplicate_range_by_atom_pick():

   """Pick two atoms and duplicate range in between."""

   def pick_range_func(*atom_specs):

      residue_specs = list(map(atom_spec_to_residue_spec, atom_specs))
      imol_1 = atom_specs[0][1]
      imol_2 = atom_specs[1][1]
      chain_id1 = atom_specs[0][2]
      chain_id2 = atom_specs[1][2]
      res_no_1 = atom_specs[0][3]
      res_no_2 = atom_specs[1][3]

      # some sanity check
      if (not imol_1 == imol_2):
         msg = (
            "BL WARNING:: not the same imols. \n"
            "imol %i and %i were selected"
            %(imol_1, imol_2))
         info_dialog_and_text(msg)
         return
      else:
         # imol ok
         if (not chain_id1 == chain_id2):
            msg = (
               "BL WARNING:: not the same chains. \n"
               "Chains %s and %s were selected"
               %(chain_id1, chain_id2))
            info_dialog_and_text(msg)
            return
         else:
            # chain ok
            
            # only allow duplication of 30 res for now (hardcode)
            # as to avoid strangeness if an atom is not selected
            # properly. And swap start, stop res if required
            res_diff = res_no_1 - res_no_2
            if (abs(res_diff) > 30):
               msg = (
                  "BL WARNING:: too many residues. \n"
                  "%i residues deExceeds the limit of 30 residues"
                  %(abs(res_diff)))
               info_dialog_and_text(msg)
               return
            elif (abs(res_diff) < 0):
               msg = (
                  "BL WARNING::  No residue selected.")
               info_dialog_and_text(msg)
               return
            if (res_no_1 > res_no_2):
               # need to swap
               res_no_1, res_no_2 = res_no_2, res_no_1
               
            # main line
            # NB: no occupancy setting here
            duplicate_residue_range(imol_1, chain_id1, res_no_1, res_no_2)

   add_status_bar_text("Pick two atoms")
   user_defined_click(2, pick_range_func)

# A simple modal dialog to accept or reject something. Pass a question
# to be asked as string.
#
def yes_no_dialog(label_text, title_text=None):
   """A simple modal dialog to accept or reject something. Pass a question
   to be asked as string.

   Args:
      label_text: the text/question to be displayed/asked in the dialog.

   Keyword args:
      title_text: title of the dialog window

   Returns:
      bool of answer.
   """
   
   dialog = gtk.Dialog(title_text, None,
                       gtk.DIALOG_MODAL | gtk.DIALOG_NO_SEPARATOR,
                       (gtk.STOCK_YES, gtk.RESPONSE_ACCEPT,
                        gtk.STOCK_NO, gtk.RESPONSE_REJECT))
   ifont = gtk.gdk.Font("fixed")
   label = gtk.Label(label_text)
   dialog.vbox.pack_end(label, True, True, 0)
   dialog.show_all()
   result = dialog.run()
   if result == gtk.RESPONSE_ACCEPT:
      ret = True
   else:
      ret = False
   dialog.destroy()
   
   return ret

   
# let the c++ part of mapview know that this file was loaded:

# print "From coot_gui.py calling set_found_coot_python_gui()" - called from rcrane loader
# too (I think)
set_found_coot_python_gui()
