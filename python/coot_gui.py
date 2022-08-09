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

# import pygtk, gtk, pango

import types
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk
import coot
import coot_utils
import coot_gui_api
import acedrg_link
import sharpen_blur
import redefine_functions as rf

# thank you ebassi!
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


try:
    import gobject
except:
    print("WARNING:: no gobject available")


global histpos
global history
histpos = 0
history = ['']


# The callback from pressing the Go button in the smiles widget, an
# interface to run libcheck.
#
def handle_smiles_go(tlc_entry, smiles_entry):

    tlc_text = tlc_entry.get_text()
    smiles_text = smiles_entry.get_text()
    use_libcheck = False
    if coot_utils.is_windows():
        use_libcheck = True
    generator_3d_import.new_molecule_by_smiles_string(
        tlc_text, smiles_text, force_libcheck=use_libcheck)

# smiles GUI
#


def smiles_gui():

    def smiles_gui_internal():
        def delete_event(*args):
            smiles_window.destroy()
            return False

        def go_button_pressed(widget, tlc_entry, smiles_entry, smiles_window):
            handle_smiles_go(tlc_entry, smiles_entry)
            smiles_window.destroy()
            delete_event()

        def smiles_connect(widget, event, tlc_entry, smiles_entry, smiles_window):
            if (event.keyval == 65293):
                handle_smiles_go(tlc_entry, smiles_entry)
                smiles_window.destroy()

        smiles_window = Gtk.Window()
        smiles_window.set_title("SMILES GUI")
        vbox = Gtk.VBox(False, 0)
        hbox1 = Gtk.HBox(False, 0)
        hbox2 = Gtk.HBox(False, 0)
        tlc_label = Gtk.Label("  3-letter code ")
        tlc_entry = Gtk.Entry(max=3)
        tlc_entry.set_text("")
        smiles_label = Gtk.Label("SMILES string ")
        smiles_entry = Gtk.Entry()
        if coot.enhanced_ligand_coot_p():
            text = Gtk.Label("  [SMILES interface works by using Pyrogen]  ")
        else:
            text = Gtk.Label(
                "  [SMILES interface works by using CCP4's LIBCHECK]  ")
        go_button = Gtk.Button("  Go  ")
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

        smiles_entry.connect("key-press-event", smiles_connect,
                             tlc_entry, smiles_entry, smiles_window)
        go_button.connect("clicked", go_button_pressed,
                          tlc_entry, smiles_entry, smiles_window)
        smiles_window.connect("delete_event", delete_event)

        smiles_window.show_all()

    # first check that libcheck is available... if not put up and info
    # dialog.
    if coot.enhanced_ligand_coot_p():
        smiles_gui_internal()
    else:
        if (coot_utils.find_exe(libcheck_exe, "CCP4_BIN", "PATH")):
            smiles_gui_internal()
        else:
            coot.info_dialog(
                "You need to setup CCP4 (specifically LIBCHECK) first.")


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

    window = Gtk.Window()
    window.set_title('Coot')
    vbox = Gtk.VBox(False, 0)
    hbox1 = Gtk.HBox(False, 0)
    hbox2 = Gtk.HBox(False, 0)
    hbox3 = Gtk.HBox(False, 0)
    function_label = Gtk.Label(function_label)
    smiles_entry = Gtk.Entry()
    cancel_button = Gtk.Button("  Cancel  ")
    go_button = Gtk.Button(go_button_label)

    vbox.pack_start(hbox1, False, False, 0)
    vbox.pack_start(hbox2, False, False, 4)
    vbox.pack_start(hbox3, False, False, 4)
    hbox3.pack_start(go_button, True, True, 6)
    hbox3.pack_start(cancel_button, True, True, 6)
    hbox1.pack_start(function_label, False, False, 0)
    hbox2.pack_start(smiles_entry, True, True, 0)
    window.add(vbox)
    vbox.set_border_width(6)

    if isinstance(entry_1_default_text, (str,)):
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
                         return_widget=False):

    def delete_event(*args):
        window.destroy()
        return False

    def go_function_event(*args):
        if check_button:
            handle_go_function(tlc_entry.get_text(),
                               smiles_entry.get_text(),
                               check_button.get_active())
        else:
            handle_go_function(tlc_entry.get_text(), smiles_entry.get_text())
            delete_event()

    def key_press_event(widget, event, tlc_entry, smiles_entry, check_button):
        if (event.keyval == 65293):
            if check_button:
                handle_go_function(tlc_entry.get_text(
                ), smiles_entry.get_text(), check_button.get_active())
            else:
                handle_go_function(tlc_entry.get_text(),
                                   smiles_entry.get_text())
                delete_event()

    window = Gtk.Window()
    window.set_title('Coot')
    vbox = Gtk.VBox(False, 0)
    hbox1 = Gtk.HBox(False, 0)
    hbox2 = Gtk.HBox(False, 0)
    hbox3 = Gtk.HBox(False, 0)
    tlc_label = Gtk.Label(label_1)
    tlc_entry = Gtk.Entry()
    smiles_label = Gtk.Label(label_2)
    smiles_entry = Gtk.Entry()
    h_sep = Gtk.HSeparator()
    cancel_button = Gtk.Button("  Cancel  ")
    go_button = Gtk.Button(go_button_label)

    vbox.pack_start(hbox1, False, False, 0)
    vbox.pack_start(hbox2, False, False, 0)
    hbox3.pack_start(go_button, False, False, 6)
    hbox3.pack_start(cancel_button, True, False, 6)
    hbox1.pack_start(tlc_label, False, False, 0)
    hbox1.pack_start(tlc_entry, False, False, 0)
    hbox2.pack_start(smiles_label, False, False, 0)
    hbox2.pack_start(smiles_entry, True, True, 0)

    try:

        def check_callback(*args):
            active_state = c_button.get_active()
            handle_check_button_function(active_state)

        c_button = Gtk.CheckButton(check_button_label)
        vbox.pack_start(c_button, False, False, 2)
        c_button.connect("toggled", check_callback)
        check_button = c_button

    except:
        check_button = False 	# the check-button when we don't want to see it

    vbox.pack_start(h_sep, True, False, 3)
    vbox.pack_start(hbox3, False, False, 0)
    window.add(vbox)
    vbox.set_border_width(6)

    if isinstance(entry_1_default_text, (str,)):
        tlc_entry.set_text(entry_1_default_text)
    else:
        print("BL WARNING:: entry_1_default_text was no string!!")

    if isinstance(entry_2_default_text, (str,)):
        smiles_entry.set_text(entry_2_default_text)
    else:
        print("BL WARNING:: entry_2_default_text was no string!!")

    go_button.connect("clicked", go_function_event)
    cancel_button.connect("clicked", delete_event)

    smiles_entry.connect("key-press-event", key_press_event,
                         tlc_entry, smiles_entry, check_button)

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
            handle_go_function([entry.get_text()
                                for entry in entries], check_button.get_active())
        else:
            handle_go_function([entry.get_text() for entry in entries])
        delete_event()

    window = Gtk.Window()
    window.set_title('Coot')
    vbox = Gtk.VBox(False, 0)
    hbox3 = Gtk.HBox(False, 0)
    h_sep = Gtk.HSeparator()
    cancel_button = Gtk.Button("  Cancel  ")
    go_button = Gtk.Button(go_button_label)

    # all the labelled entries
    #
    entries = []
    for entry_info in entry_info_list:

        entry_1_hint_text = entry_info[0]
        entry_1_default_text = entry_info[1]
        hbox1 = Gtk.HBox(False, 0)

        label = Gtk.Label(entry_1_hint_text)
        entry = Gtk.Entry()
        entries.append(entry)

        try:
            entry.set_text(entry_1_default_text)
        except:
            pass

        hbox1.pack_start(label, False, False, 0)
        hbox1.pack_start(entry, False, False, 0)
        vbox.pack_start(hbox1,  False, False, 0)

    # print "debug:: check-button-info: ", check_button_info
    # print "debug:: entry-info-list: ", entry_info_list

    try:
        def check_callback(*args):
            active_state = c_button.get_active()
            check_button_info[1] = active_state

        c_button = Gtk.CheckButton(check_button_info[0])
        vbox.pack_start(c_button, False, False, 2)
        c_button.connect("toggled", check_callback)
        check_button = c_button
    except:
        check_button = False

    print("debug:: Here check button creation........ check-button is ", check_button)
    vbox.pack_start(h_sep, True,  False, 3)
    vbox.pack_start(hbox3, False, False, 0)
    window.add(vbox)
    vbox.set_border_width(6)

    hbox3.pack_start(go_button,     True, False, 6)
    hbox3.pack_start(cancel_button, True, False, 6)

    cancel_button.connect("clicked", delete_event)
    go_button.connect("clicked", go_function_event)

    window.show_all()


def molecule_centres_gui():

    def delete_event(*args):
        window.destroy()
        return False

    def callback_func(widget, molecule_number, label):
        s = "Centred on " + label
        coot.add_status_bar_text(s)
        coot.set_rotation_centre(*molecule_centre(molecule_number))

    # first, we create a window and a frame to be put into it.
    #
    # we create a vbox (a vertical box container) that will contain the
    # buttons for each of the coordinate molecules
    #
    window = Gtk.Window()
    frame = Gtk.Frame("Molecule Centres")
    vbox = Gtk.VBox(False, 3)

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
    for molecule_number in coot_utils.molecule_number_list():
        if (coot.is_valid_model_molecule(molecule_number)):
            name = coot.molecule_name(molecule_number)
            label = str(molecule_number) + " " + name
            button = Gtk.Button(label)
            button.connect("clicked", callback_func, molecule_number, label)
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
            print(
                "BL WARNING:: libcheck didn't seem to run ok! Please check output carefully!!")
        else:
            pass
    fp.close()

# old coot test


def old_coot_qm():

    import time
    import random
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
    new_release_time = 1352674800  # 12 Nov 2012
    current_time = int(time.time())
    time_diff = current_time - new_release_time
    if (time_diff > 0):
        if (time_diff > 8600):
            s = "You've got an Old Coot!\n\nIt's time to upgrade."
        else:
            if (random.randint(0, 10) == 0):
                # Jorge Garcia:
                s = "(Nothing says \"patriotism\" like an Ireland shirt...)\n"
            else:
                s = "You've got an Old Coot!\n\nIt's time to upgrade."
        coot.info_dialog(s)

# if (not coot.coot_has_guile()):
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


def interesting_things_gui(title, baddie_list):

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

            if (isinstance(func_maybe1, list) and len(func_maybe1) > 0):
                # the last one is probably a funcn (no button name)
                func_maybe_strip = func_maybe1[0]
#             print "BL DEBUG:: func_maybe_strip is", func_maybe_strip
                if (callable(func_maybe_strip)):
                    return func_maybe1, False, False
                else:
                    return False, False, False

            elif (isinstance(func_maybe2, list) and len(func_maybe2) > 0
                  and isinstance(func_maybe1, bytes)):
                # the second last is function, last is button name
                func_maybe_strip = func_maybe2[0]
                button_name = func_maybe1
                if (callable(func_maybe_strip)):
                    return func_maybe2, button_name, False
                else:
                    return False, False, False

            elif (isinstance(func_maybe3, list) and len(func_maybe3) > 0
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
        Gtk.main_quit()
        return False

    def callback_func1(widget, coords):
        coot.set_rotation_centre(*coords)

    def callback_func2(widget, mol_no, atom_info):
        print("Attempt to go to chain: %s resno: %s atom-name: %s" %
              (atom_info[0], atom_info[1], atom_info[2]))
        coot.set_go_to_atom_molecule(mol_no)
        success = coot.set_go_to_atom_chain_residue_atom_name(*atom_info)
        if success == 0:           # failed?!
            new_name = coot.unmangle_hydrogen_name(atom_info[2])
            success2 = coot.set_go_to_atom_chain_residue_atom_name(atom_info[0], atom_info[1], new_name)
            if success2 == 0:
                print("Failed to centre on demangled name: ", new_name)
                coot.set_go_to_atom_chain_residue_atom_name(atom_info[0], atom_info[1], " CA ")

    # main body
    # to accomodated tooltips we need to either have a Gtk.Window with Gtk.main()
    # or a dialog and run() it! We use Window to not block events and hope not to
    # interfere with the gtk_main() of coot itself
    window = Gtk.Window()
    scrolled_win = Gtk.ScrolledWindow()
    outside_vbox = Gtk.VBox(False, 2)
    inside_vbox = Gtk.VBox(False, 0)

    window.set_default_size(250, 250)
    window.set_title(title)
    inside_vbox.set_border_width(4)

    window.add(outside_vbox)
    outside_vbox.add(scrolled_win)
    scrolled_win.add_with_viewport(inside_vbox)
    #tooltips = Gtk.Tooltips()

    for baddie_items in baddie_list:

        if baddie_items == '':
            print('Done')
        else:
            hbox = Gtk.HBox(False, 0)  # hbox to contain Baddie button and
            # and fix it button

            label = baddie_items[0]
            button = Gtk.Button(label)

            inside_vbox.pack_start(hbox, False, False, 2)
            hbox.pack_start(button, True, True, 1)

            # add the a button for the fix func if it exists.  Add
            # the callback.
            fix_func, button_name, tooltip_str = baddie_had_fix_qm(
                baddie_items)
            if (fix_func):
                if (button_name):
                    fix_button_name = button_name
                else:
                    fix_button_name = "  Fix  "
                fix_button = Gtk.Button(fix_button_name)
                hbox.pack_end(fix_button, False, False, 2)
                fix_button.show()
                fix_button.connect("clicked", fix_func_call, fix_func)

            if (tooltip_str):
                # we have a tooltip str
                if Gtk.pygtk_version >= (2, 12):
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
    ok_button = Gtk.Button("  OK  ")
    outside_vbox.pack_start(ok_button, False, False, 6)
    ok_button.connect("clicked", delete_event)

    window.connect("destroy", delete_event)

    window.show_all()
    Gtk.main()

#interesting_things_gui("Bad things by Analysis X",[["Bad Chiral",0,"A",23,"","CA","A"],["Bad Density Fit",0,"B",65,"","CA",""],["Interesting blob",45.6,46.7,87.5],["Interesting blob 2",45.6,41.7,80.5]])
#interesting_things_gui("Bad things by Analysis X",[["Bad Chiral",0,"A",23,"","CA","A",[print_sequence_chain,0,'A']],["Bad Density Fit",0,"B",65,"","CA",""],["Interesting blob",45.6,46.7,87.5],["Interesting blob 2",45.6,41.7,80.5]])


# Fill an option menu with the "right type" of molecules.  If
# filter_function returns True then add it.  Typical value of
# filter_function is valid-model-molecule_qm
#
def fill_option_menu_with_mol_options(menu, filter_function):

    mol_ls = []

    for mol_no in coot_utils.molecule_number_list():
        if filter_function(mol_no):
            label_str = coot.molecule_name(mol_no)
            if (isinstance(label_str, (str,))):
                mlabel_str = str(mol_no) + " " + label_str
                menu.append_text(mlabel_str)
                menu.set_active(0)
                mol_ls.append(mol_no)
            else:
                print("OOps molecule name for molecule %s is %s" %
                      (mol_no_ls, label_str))
    return mol_ls

# Fill an option menu with maps and return the list of maps
#


def fill_option_menu_with_map_mol_options(menu):
    return fill_option_menu_with_mol_options(menu, coot_utils.valid_map_molecule_qm)

def fill_option_menu_with_map_mol_with_associated_data_options(menu):
    return fill_option_menu_with_mol_options(menu, valid_map_with_associated_data_molecule_qm)

def fill_option_menu_with_difference_map_mol_options(menu):
    return fill_option_menu_with_mol_options(menu, is_difference_map_qm)

# Helper function for molecule chooser.  Not really for users.
#
# Return a list of models, corresponding to the menu items of the
# option menu.
#
# The returned list will not contain references to map or closed
# molecules.
#


def fill_option_menu_with_coordinates_mol_options(menu):
    return fill_option_menu_with_mol_options(menu, coot_utils.valid_model_molecule_qm)

#


def fill_option_menu_with_number_options(menu, number_list, default_option_value):

    print("************************* Get Rid of this! **********************")

    for number in coot_utils.number_list:
        mlabel_str = str(number)
        menu.append_text(mlabel_str)
        if (default_option_value == number):
            count = coot_utils.number_list.index(number)
            menu.set_active(count)
            print("setting menu active ", default_option_value, count)

def fill_combobox_with_number_options(combobox, number_list, active_value):

    def make_store(number_list):
        name_store = Gtk.ListStore(int, str)
        for i in number_list:
            label_str = str(i)
            name_store.append([imol, label_str])
        return name_store

    combobox_items = make_store(number_list)
    renderer_text = Gtk.CellRendererText()
    combox.set_entry_text_column(1)
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)
    for i in number_list:
        print("debug: comparing", i, active_value)
        if i == active_value:
            combobox.set_active(i)
    return combobox_items


# Helper function for molecule chooser.  Not really for users.
#
# return the molecule number of the active item in the option menu,
# or return False if there was a problem (e.g. closed molecule)
#
# BL says:: we do it for gtk_combobox instead! option_menu is deprecated
#
# Get rid of this. Search for references and replace them
#

# 20220228-PE commented out on merge - let's usee get_combobox_active_molecule() now
# def get_option_menu_active_molecule(option_menu, model_mol_list):

#     model = option_menu.get_model()
#     active_item = option_menu.get_active()
#     # combobox has no children as such, so we just count the rows
#     children = len(model)

#     if (children == len(model_mol_list)):
#         try:
#             all_model = model[active_item][0]
#             imol_model, junk = all_model.split(' ', 1)

#             return int(imol_model)
#         except:
#             print("INFO:: could not get active_item")
#             return False
#     else:
#         print("Failed children length test : ", children, model_mol_list)
#         return False
#        print "get_option_menu_active_molecule(): Failed children length test : ",children, model_mol_list
#        return False

def get_combobox_active_molecule(combobox, model_mol_list):

    model = combobox.get_model()
    active_item = combobox.get_active()
    # combobox has no children as such, so we just count the rows
    children = len(model)

    if (children == len(model_mol_list)):
       try:
          all_model = model[active_item][0]
          imol_model, junk = all_model.split(' ', 1)
       
          return int(imol_model)
       except:
          print("WARNING:: could not get active_item", combobox)
          return False
    else:
       print("WARNING:: get_combobox_active_molecule(): Failed children length test : ", children, model_mol_list)
       return False


# Get rid of this also
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
        print("Failed children length test : ", children, item_list)
        return False

# see acedrg_link for how this is used
def make_store_for_string_list_combobox(combobox, string_list, default_string):
    mol_store = Gtk.ListStore(str)
    for item in string_list:
        mol_store.append([item])
    return mol_store

def make_store_for_molecule_combobox(filter_function):
    mol_store = Gtk.ListStore(int, str)
    for imol in coot_utils.molecule_number_list():
        if filter_function(imol) == 1:
            label_str = coot.molecule_name(imol)
            m_label_str = str(imol) + ' ' + label_str
            mol_store.append([imol, m_label_str])
    return mol_store

def make_store_for_map_molecule_combobox():
    mol_store = Gtk.ListStore(int, str)
    for imol in coot_utils.molecule_number_list():
        if coot.is_valid_map_molecule(imol) == 1:
            label_str = coot.molecule_name(imol)
            m_label_str = str(imol) + ' ' + label_str
            mol_store.append([imol, m_label_str])
    return mol_store

def make_store_for_diff_map_molecule_combobox():
    mol_store = Gtk.ListStore(int, str)
    for imol in coot_utils.molecule_number_list():
        if coot.is_valid_map_molecule(imol) == 1:
            if coot.map_is_difference_map(imol) == 1:
                label_str = coot.molecule_name(imol)
                m_label_str = str(imol) + ' ' + label_str
                mol_store.append([imol, m_label_str])
    return mol_store


def make_store_for_model_molecule_combobox():
    mol_store = Gtk.ListStore(int, str)
    for imol in coot_utils.molecule_number_list():
        if coot.is_valid_model_molecule(imol) == 1:
            label_str = coot.molecule_name(imol)
            m_label_str = str(imol) + ' ' + label_str
            print("debug:: make_store_for_model_molecule_combobox appending", imol, m_label_str)
            mol_store.append([imol, m_label_str])
    return mol_store

def fill_combobox_with_molecule_options(combobox, filter_function):
    mols_ls = []
    name_store = Gtk.ListStore(int, str)
    for imol in coot_utils.molecule_number_list():
        if filter_function(imol):
            label_str = coot.molecule_name(imol)
            m_label_str = str(imol) + ' ' + label_str
            name_store.append([imol, m_label_str])
            mols_ls.append(imol)

    print("###### debug fill_combobox_with_molecule_options() returning name_store", name_store)
    return name_store


def fill_combobox_with_model_molecule_options(combobox):
    name_store = Gtk.ListStore(int, str)
    for imol in coot_utils.molecule_number_list():
        if coot.is_valid_model_molecule(imol):
            label_str = coot.molecule_name(imol)
            m_label_str = str(imol) + ' ' + label_str
            name_store.append([imol, m_label_str])

    return name_store

# Typically option_menu_fill_function is
# fill_option_menu_with_coordinates_mol_options
#


def molecule_chooser_gui_generic(chooser_label, callback_function, molecule_filter_function):

    def delete_event(*args):
        window.destroy()
        return False

    def on_ok_button_clicked(*args):
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol_model = it[0]
            callback_function(imol_model)
            window.destroy()

    def on_mol_combobox_changed(combobox):
        # this function is not useful. For this dialog, we want to do things when
        # the "OK" button is pressed
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            print("Selected: imol=%s" % it)

    window = Gtk.Window(title="Molecule Chooser")
    window.set_title('Coot: Molecule Chooser')
    label = Gtk.Label(chooser_label)
    vbox = Gtk.VBox(False, 6)
    hbox_buttons = Gtk.HBox(False, 5)

    # -------- replacing an option menu of molecules: here's how to do it --------------
    #          (also see the on_ok_button_clicked callback)

    # option_menu = Gtk.combo_box_new_text()
    combobox_items = make_store_for_molecule_combobox(molecule_filter_function)
    combobox = Gtk.ComboBox.new_with_model(combobox_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_items) > 0:
        combobox.set_active(0)
    combobox.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)
    # this line is often not needed in other cases.
    combobox.connect("changed", on_mol_combobox_changed)

    # ----------------------------------------------------------------------------------

    ok_button = Gtk.Button("  OK  ")
    cancel_button = Gtk.Button(" Cancel ")
    h_sep = Gtk.HSeparator()

    vbox.pack_start(label, False, False, 5)
    vbox.pack_start(combobox, True, True, 0)
    vbox.pack_start(h_sep, True, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 5)
    hbox_buttons.pack_start(ok_button, True, False, 5)
    hbox_buttons.pack_start(cancel_button, True, False, 5)

    window.add(vbox)

    # button callbacks:
    ok_button.connect("clicked", on_ok_button_clicked, combobox)
    cancel_button.connect("clicked", delete_event)

    window.show_all()


# molecule_chooser_gui("test-gui",print_sequence_chain(0,"A"))

# Fire up a coordinates/model molecule chooser dialog, with a given
# label and on OK we call the callback_fuction with an argument of
# the chosen molecule number.
#
# chooser_label is a directive to the user such as "Choose a Molecule"
#
# callback_function is a function that takes a molecule number as an
# argument.
#
def molecule_chooser_gui(label, callback_fn):
    molecule_chooser_gui_generic(label, callback_fn, coot.is_valid_model_molecule)

# Fire up a map molecule chooser dialog, with a given label and on OK we
# call the callback_fuction with an argument of the chosen molecule
# number.
#
# chooser_label is a directive to the user such as "Choose a Molecule"
#
# callback_function is a function that takes a molecule number as an
# argument.
#


def map_molecule_chooser_gui(label, callback_fn):
    molecule_chooser_gui_generic(label, callback_fn, coot.is_valid_map_molecule)

# A pair of widgets, a molecule chooser and an entry.  The
# callback_function is a function that takes a molecule number and a
# text string.
#


def generic_chooser_and_entry(chooser_label, entry_hint_text,
                              default_entry_text, callback_function,
                              always_dismiss_on_ok_clicked=True):

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

    def combobox_to_molecule_number(combobox):
        imol = -1
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol = it[0]
        return imol

    def on_ok_button_clicked(*args):
        # active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)
        active_mol_no = combobox_to_molecule_number(combobox)

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
        except KeyError as e:
            print("WARNING:: Failed to get a (molecule) number")
            print("WARNING:: KeyError Failed to get a (molecule) number", e)

    window = Gtk.Window()
    window.set_title('Coot')
    label = Gtk.Label(chooser_label)
    vbox = Gtk.VBox(False, 2)
    hbox_for_entry = Gtk.HBox(False, 0)
    entry = Gtk.Entry()
    entry_label = Gtk.Label(entry_hint_text)
    hbox_buttons = Gtk.HBox(True, 2)
    # option_menu = Gtk.combo_box_new_text()

    combobox_items = make_store_for_model_molecule_combobox()
    combobox = Gtk.ComboBox.new_with_model(combobox_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_items) > 0:
        combobox.set_active(0)
    combobox.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)

    ok_button = Gtk.Button("  OK  ")
    cancel_button = Gtk.Button(" Cancel ")
    h_sep = Gtk.HSeparator()

    window.set_default_size(400, 100)
    window.add(vbox)
    vbox.pack_start(label, False, False, 5)
    vbox.pack_start(combobox, True, True, 0)
    vbox.pack_start(hbox_for_entry, False, False, 5)
    if check_button_label:
        check_button = Gtk.CheckButton(check_button_label)
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
                      combobox, callback_function, check_button)
    cancel_button.connect("clicked", delete_event)

    window.show_all()

# Create a window
#
# Return a pair of widgets, a chooser entry and a file selector.  The
# callback_function is a function that takes a molecule number and a
# text string (e.g. chain_id and file_name)
#
# chooser_filter is typically coot_utils.valid_map_molecule_qm or coot_utils.valid_model_molecule_qm
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

    def get_molecule():
        tree_iter = combobox.get_active_iter()
        imol = -1
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol = it[0]
        return imol

    def on_ok_button_clicked(*args):
        # what is the molecule number of the option menu?
        # active_mol_no = get_option_menu_active_molecule(option_menu, model_mol_list)
        imol = get_molecule()
        print("on_ok_button_clicked, got imol", imol)

        try:
            active_mol_no = int(imol)
            text = entry.get_text()
            file_sel_text = file_sel_entry.get_text()
            if (c_button and c_button.get_active()):
                # use alt function
                alternative_callback_function(imol, text, file_sel_text)
            else:
                callback_function(imol, text, file_sel_text)
            delete_event()
        except KeyError as e:
            print(e)
            print("Failed to get a (molecule) number")

    window = Gtk.Window()
    window.set_title("Coot Chooser")
    label = Gtk.Label(chooser_label)
    vbox = Gtk.VBox(False, 2)
    hbox_for_entry = Gtk.HBox(False, 2)
    entry = Gtk.Entry()
    entry_label = Gtk.Label(entry_hint_text)
    hbox_buttons = Gtk.HBox(True, 2)
    ok_button = Gtk.Button("  OK  ")
    cancel_button = Gtk.Button(" Cancel ")
    h_sep = Gtk.HSeparator()
    # model_mol_list = fill_option_menu_with_mol_options(option_menu, chooser_filter)
    # name_store = fill_combobox_with_model_molecule_options(combobox)
    combobox = Gtk.ComboBox()
    combobox_mol_items = make_store_for_model_molecule_combobox()
    combobox.set_model(combobox_mol_items)
    renderer_text = Gtk.CellRendererText()
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)
    combobox.set_active(0)

    window.set_default_size(400, 100)
    window.add(vbox)
    vbox.pack_start(label, False, False, 5)
    vbox.pack_start(combobox, True, True, 0)
    vbox.pack_start(hbox_for_entry, False, False, 5)
    hbox_buttons.pack_start(ok_button, True, False, 5)
    hbox_buttons.pack_start(cancel_button, False, False, 5)
    hbox_for_entry.pack_start(entry_label, False, False, 4)
    hbox_for_entry.pack_start(entry, True, True, 4)
    entry.set_text(default_entry_text)

    c_button = None
    if use_check_button:
        # now add a check button
        c_button = Gtk.CheckButton(check_button_label)
        vbox.pack_start(c_button, False, False, 2)

    file_sel_entry = file_chooser_entry(window, vbox, file_selector_hint)
    vbox.pack_start(h_sep, True, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 5)

    # button callbacks
    ok_button.connect("clicked", on_ok_button_clicked, entry, combobox,
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
# chooser_filter is typically coot_utils.valid_map_molecule_qm or coot_utils.valid_model_molecule_qm
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

    def combobox_to_molecule_number(combobox):
        imol = -1
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol = it[0]
        return imol

    def on_ok_button_clicked(*args):
        # what is the molecule number of the option menu?
        active_mol_no = combobox_to_molecule_number(combobox)

        try:
            active_mol_no = int(active_mol_no)
            file_sel_text = file_sel_entry.get_text()
            callback_function(active_mol_no, file_sel_text)
            delete_event()
        except:
            print("Failed to get a (molecule) number")

    window = Gtk.Window()
    label = Gtk.Label(chooser_label)
    vbox = Gtk.VBox(False, 2)
    hbox_for_entry = Gtk.HBox(False, 0)
    hbox_buttons = Gtk.HBox(True, 2)
    # option_menu = Gtk.combo_box_new_text()
    combobox_items = make_store_for_model_molecule_combobox()
    combobox = Gtk.ComboBox.new_with_model(combobox_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_items) > 0:
        combobox.set_active(0)
    combobox.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)

    ok_button = Gtk.Button("  OK  ")
    cancel_button = Gtk.Button(" Cancel ")
    h_sep = Gtk.HSeparator()

    window.set_default_size(400, 100)
    window.add(vbox)
    vbox.pack_start(label, False, False, 5)
    vbox.pack_start(combobox, True, True, 0)
    hbox_buttons.pack_start(ok_button, True, False, 5)
    hbox_buttons.pack_start(cancel_button, False, False, 5)

    file_sel_entry = file_chooser_entry(window, vbox, file_selector_hint, default_file_name)
    vbox.pack_start(h_sep, True, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 5)

    # button callbacks
    ok_button.connect("clicked", on_ok_button_clicked, option_menu, callback_function)
    cancel_button.connect("clicked", delete_event)

    window.show_all()


# If a menu with label menu-label is not found in the coot main
# menubar, then create it and return it.
# If it does exist, simply return it.
#
def coot_menubar_menu(menu_label):

    try:
        coot_main_menubar = coot_gui_api.main_menubar()

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
                # print "BL DEBUG:: found menu label is ", f
                # we shall return the submenu and not the menuitem
                found_menu = f[1].get_submenu()
        if found_menu:
            return found_menu
        else:
            menu = Gtk.Menu()
            menuitem = Gtk.MenuItem(menu_label)
            menuitem.set_submenu(menu)
            coot_main_menubar.append(menuitem)
            menuitem.show()
            return menu
    except:
        print("ERROR:: coot_main_menubar() an error occurs using coot_gui_api")


# Given that we have a menu (e.g. one called "Extensions") provide a
# cleaner interface to adding something to it:
#
# activate_function is a thunk.
#
def add_simple_coot_menu_menuitem(menu, menu_item_label, activate_function):

    submenu = Gtk.Menu()  # what this for?
    sub_menuitem = Gtk.MenuItem(menu_item_label)

    menu.append(sub_menuitem)
    sub_menuitem.show()

    sub_menuitem.connect("activate", activate_function)


# Make an interesting things GUI for residues of molecule number
# imol that have alternate conformations.
#
def alt_confs_gui(imol):

    residues_list = coot_utils.residues_with_alt_confs(imol)
    interesting_residues_gui(imol, "Residues with Alt Confs", residues_list)

# Make an interesting things GUI for residues with missing atoms
#


def missing_atoms_gui(imol):

    interesting_residues_gui(imol, "Residues with missing atoms",
                             coot.missing_atom_info_py(imol))

# Make an interesting things GUI for residues with zero occupancy atoms
#


def zero_occ_atoms_gui(imol):

    atom_ls = coot_utils.atoms_with_zero_occ(imol)
    if atom_ls:
        interesting_things_gui("Residues with zero occupancy atoms",
                               atom_ls)
    else:
        s = "No atoms with zero occupancy found\n"
        s += "in molecule "
        s += str(imol)
        coot.info_dialog(s)


# Make an interesting things GUI for residues of molecule number
# imol for the given imol.   A generalization of alt-confs gui.
#
def interesting_residues_gui(imol, title, interesting_residues):

    # from types import ListType
    centre_atoms = []
    if coot_utils.valid_model_molecule_qm(imol):
        residues = interesting_residues
        for i in range(len(residues)):
            # if we get here with a 4-element spec, then fix that now
            sp = residues[i]
            if len(sp) == 4:
                residues[i] = sp[1:]

        for spec in residues:
            if spec:
                centre_atoms.append(coot_utils.residue_spec_to_atom_for_centre(imol, *spec)) # debracket spec
            else:
                centre_atoms.append([False])

        interesting_things_gui(title,
                    list(map(lambda residue_cpmd, centre_atom:
                             [residue_cpmd[0] + " " +
                              str(residue_cpmd[1]) + " " +
                              residue_cpmd[2] + " " +
                              centre_atom[0] + " " + centre_atom[1],
                              imol,
                              residue_cpmd[0],
                              residue_cpmd[1],
                              residue_cpmd[2],
                              centre_atom[0],
                              centre_atom[1]] if centre_atom else
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

    coot_main_toolbar = coot_gui_api.main_toolbar()

    # main body
    #
    found_button = False
    for f in toolbar_label_list():
        print("#### debug coot_toolbar_button", f)
    for f in toolbar_label_list():
        if button_label in f:
            found_button = f[1]
    if found_button:
        # here we only try to add a new icon, we cannot overwrite the callback function!
        toolbutton = found_button
    else:
        if toggle_button:
            toolbutton = Gtk.ToggleToolButton()
            toolbutton.set_label(button_label)
        else:
            toolbutton = Gtk.ToolButton(icon_widget=None, label=button_label)
        coot_main_toolbar.insert(toolbutton, -1)       # insert at the end
        toolbutton.set_is_important(True)              # to display the text,
        # otherwise only icon

        def cb_wrapper(widget, callback_function):
            print("DEBUG:: in cb_wrapper() for coot_toolbar_button(): widget is", widget, "and callback_function is", callback_function)
            if True:
                # have function as callable and maybe extra args (all in one list)
                args = []
                function = callback_function
                if type(callback_function) is list:
                    function = callback_function[0]
                    args = callback_function[1:]
                # pass the widget/button as well? Maybe the cb function can
                # make use of it
                if use_button:
                    args.append(widget)
                if callable(function):
                    print("DEBUG:: in cb_wrapper() was callable", function, "with args", *args)
                    function(*args)
                else:
                    print("BL ERROR:: cannot evaluate or call function", function)
        #toolbutton.connect("clicked", lambda w: eval(cb_function))
        toolbutton.connect("clicked", cb_wrapper, cb_function)
        toolbutton.show()

    if icon_name:
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

    coot_main_toolbar = coot_gui_api.main_toolbar()

    # main body
    #
    found_combobox = False
    for f in toolbar_label_list():
        if label in f:
            found_combobox = f[1]
    if found_combobox:
        s = "BL WARNING:: already have a comboxbox with name %s. Please " +\
            "try again!" % found_combobox
        coot.info_dialog(s)
        return False
    else:
        toolitem = Gtk.ToolItem()
        toolitem.set_name(label)
        combobox = Gtk.combo_box_new_text()
        for text in entry_list:
            combobox.append_text(text)
        combobox.set_active(0)

        # tooltips?
        if tooltip:
            if Gtk.pygtk_version >= (2, 12):
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

    coot_main_toolbar = coot_gui_api.main_toolbar()
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
    window = Gtk.Window()
    scrolled_win = Gtk.ScrolledWindow()
    outside_vbox = Gtk.VBox(False, 2)
    inside_vbox = Gtk.VBox(False, 0)

    window.set_default_size(250, 250)
    window.set_title(dialog_name)
    inside_vbox.set_border_width(4)

    window.add(outside_vbox)
    outside_vbox.add(scrolled_win)
    scrolled_win.add_with_viewport(inside_vbox)

    for button_item in button_list:
        if button_item and len(button_item) == 2:
            button_label = button_item[0]
            action = button_item[1]

            button = Gtk.Button(button_label)
            inside_vbox.pack_start(button, False, False, 2)
            button.connect("clicked", action)
            button.show()

    outside_vbox.set_border_width(4)
    ok_button = Gtk.Button("  OK  ")
    outside_vbox.pack_start(ok_button, False, False, 6)
    ok_button.connect("clicked", delete_event)

    window.show_all()


# Generic interesting things gui: user passes a function that takes 4
# args: the chain-id, resno, inscode and residue-serial-number
# (should it be needed) and returns either #f or something
# interesting (e.g. a label/value).  It is the residue-test-func of
# the residue-matching-criteria function.
#
def generic_interesting_things(imol, gui_title_string, residue_test_func):

    if coot_utils.valid_model_molecule_qm(imol):

        interesting_residues = coot_utils.residues_matching_criteria(
            imol, residue_test_func)
        for i in range(len(interesting_residues)):
            interesting_residues[i][0] = imol
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

# A gui that makes a generic number chooser: the go_function is a
# function that takes the value of the active menu item - as a
# number.
#
def generic_number_chooser(number_list, default_option_value, hint_text,
                           go_button_label, go_function):

    def delete_event(*args):
        window.destroy()
        return False

    def go_button_pressed(button, combobox, go_function):

        print("########### go_button_pressed: button", button)
        print("########### go_button_pressed: combobox", combobox)
        print("########### go_button_pressed: go_function", go_function)
        active_item_index = combobox.get_active()
        print("########### go_button_pressed: active_item_index", active_item_index)

        imol = -1
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol = it[0]
            pass
        print("########### go_button_pressed: imol", imol)
        go_function(imol)
        # print("Failed to get execute function")
        delete_event()

    window = Gtk.Window()
    vbox = Gtk.VBox(False, 0)
    hbox1 = Gtk.HBox(False, 0)
    hbox2 = Gtk.HBox(True, 0)      # for Go and Cancel
    function_label = Gtk.Label(hint_text)
    h_sep = Gtk.HSeparator()
    go_button = Gtk.Button(go_button_label)
    cancel_button = Gtk.Button("  Cancel  ")

    # fill_option_menu_with_number_options(option_menu, number_list, default_option_value)

    def make_store(number_list):
        name_store = Gtk.ListStore(int, str)
        for i in number_list:
            label_str = str(i)
            name_store.append([i, label_str])
        return name_store

    combobox_items = make_store(number_list)
    combobox = Gtk.ComboBox.new_with_model(combobox_items)
    renderer_text = Gtk.CellRendererText()
    combobox.set_entry_text_column(1)
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)
    for i in number_list:
        if i == default_option_value:
            combobox.set_active(i)

    vbox.pack_start(hbox1, True, False, 0)
    vbox.pack_start(function_label, True, False, 0)
    vbox.pack_start(combobox, True, False, 0)
    vbox.pack_start(h_sep, False, False, 0)
    vbox.pack_start(hbox2, False, False, 0)
    hbox2.pack_start(go_button, False, False, 6)
    hbox2.pack_start(cancel_button, False, False, 6)
    window.add(vbox)
    vbox.set_border_width(6)
    hbox1.set_border_width(6)
    hbox2.set_border_width(6)
    go_button.connect("clicked", go_button_pressed, combobox, go_function)
    cancel_button.connect("clicked", delete_event)

    window.show_all()

# vbox is the vbox to which this compound widget should be added.
# button-press-func is the lambda function called on pressing return
# or the button, which takes one argument, the entry.
#
# Add this widget to the vbox here.
#


def entry_do_button(vbox, hint_text, button_label, button_press_func, entry_text=False):

    hbox = Gtk.HBox(False, 0)
    entry = Gtk.Entry()
    button = Gtk.Button(button_label)
    label = Gtk.Label(hint_text)

    hbox.pack_start(label, False, False, 2)
    hbox.pack_start(entry, True, True, 2)
    hbox.pack_start(button, False, False, 2)
    button.connect("clicked", button_press_func, entry)

    if entry_text:
        entry.set_text(entry_text)
    label.show()
    button.show()
    entry.show()
    hbox.show()
    vbox.pack_start(hbox, True, False, 2)
    return entry

# pack a hint text and a molecule chooser option menu into the given vbox.
#
# return the option-menu and model molecule list:


def generic_molecule_chooser(hbox, hint_text):

    menu = Gtk.Menu()
    # option_menu = Gtk.combo_box_new_text()
    combobox_items = make_store_for_model_molecule_combobox()
    combobox = Gtk.ComboBox.new_with_model(combobox_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_items) > 0:
        combobox.set_active(0)
    combobox.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)

    label = Gtk.Label(hint_text)

    hbox.pack_start(label, False, False, 2)
    hbox.pack_start(combobox, True, True, 2)

    # we used to return the molecule list here also, but the combo box can
    # get its own active item using combobox_to_molecule_number()
    return combobox

# Return an entry, the widget is inserted into the hbox passed to
# this function
#


# This is the same as the file_selector_entry, but using the modern FileChooser
# Return an entry, the widget is inserted into the hbox passed to
# this function
#


def file_chooser_entry(parent_window, hbox, hint_text, default_file_name=False):

    #$ if Gtk.pygtk_version > (2, 3, 90):
    if True:

        def file_func1(*args):

            def file_ok_sel(*args):
                t = file_chooser_dialog.get_filename()
                print(t)
                entry.set_text(t)
                file_chooser_dialog.destroy()

            file_chooser_dialog = Gtk.FileChooserDialog("Please choose a file", parent_window,
                                                        Gtk.FileChooserAction.OPEN,
                                                        (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                                         Gtk.STOCK_OPEN, Gtk.ResponseType.OK))

            response = file_chooser_dialog.run()
            if response == Gtk.ResponseType.OK:
                file_ok_sel(file_chooser_dialog, entry)
            elif response == Gtk.ResponseType.CANCEL:
                file_chooser_dialog.destroy()

        vbox = Gtk.VBox(False, 0)
        entry = entry_do_button(vbox, hint_text, "  File...  ", file_func1, default_file_name)
        hbox.pack_start(vbox, False, False, 2)
        vbox.show()
        return entry
    else:
        return False

# The gui for the strand placement function
#


def place_strand_here_gui():

    generic_number_chooser(coot_utils.number_list(4, 12), 7,
                           " Estimated number of residues in strand",
                           "  Go  ",
                           lambda n: coot.place_strand_here(n, 15))

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
            imol = n_atom[0]
            chain_id = n_atom[1]
            resno = n_atom[2]
            inscode = n_atom[3]
            at_name = n_atom[4]
            alt_conf = n_atom[5]
            cootaneer_results = cootaneer(imol_map, imol, [chain_id, resno, inscode,
                                                           at_name, alt_conf])
            print("Cootaneering status:", cootaneer_results)
            if (cootaneer_results == 0):
                s = "Insufficiently confident in alignment to make a fit." + \
                    "\n" + \
                    "Perhaps you could improve or extend this fragment."
                window.destroy()
                coot.info_dialog(s)
        else:
            print("BL WARNING:: no close atom found!")
            window.destroy()

    def add_text_to_text_box(text_box, description):
        start = text_box.get_start_iter()
        text_box.create_tag("tag", foreground="black",
                            background="#c0e6c0")
        text_box.insert_with_tags_by_name(start, description, "tag")

        # return the (entry . textbuffer/box)
        #
    def entry_text_pair_frame(seq_info):

        frame = Gtk.Frame()
        vbox = Gtk.VBox(False, 3)
        entry = Gtk.Entry()
        textview = Gtk.TextView()
        textview.set_wrap_mode(Gtk.WRAP_WORD_CHAR)
        text_box = textview.get_buffer()
        chain_id_label = Gtk.Label("Chain ID")
        sequence_label = Gtk.Label("Sequence")

        frame.add(vbox)
        vbox.pack_start(chain_id_label, False, False, 2)
        vbox.pack_start(entry, False, False, 2)
        vbox.pack_start(sequence_label, False, False, 2)
        vbox.pack_start(textview, False, False, 2)
        add_text_to_text_box(text_box, seq_info[1])
        entry.set_text(seq_info[0])
        return [frame, entry, text_box]

        # main body
    imol_map = coot.imol_refinement_map()
    if (imol_map == -1):
        coot.show_select_map_dialog()
        print("BL DEBUG:: probably should wait here for input!?")

    window = Gtk.Window()
    outside_vbox = Gtk.VBox(False, 2)
    inside_vbox = Gtk.VBox(False, 2)
    h_sep = Gtk.HSeparator()
    buttons_hbox = Gtk.HBox(True, 2)
    go_button = Gtk.Button("  Cootaneer!  ")
    cancel_button = Gtk.Button("  Cancel  ")

    seq_info_ls = sequence_info(imol)
    # print "BL DEBUG:: sequence_list and imol is", seq_info_ls, imol

    if not seq_info_ls:
        s = "No sequence assigned for molecule number " + str(imol)
        print(s)
        coot.info_dialog(s)
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
                jview_name = coot.view_name_py(jview)
                if jview >= coot.n_views():
                    return strr
                elif jview_name == False:
                    return strr
                elif strr == jview_name:
                    view_count += 1
                    break
                else:
                    jview += 1
        return strr

    def add_view_local_func(text):
        new_view_number = coot.add_view_here(text)
        add_view_to_views_panel(text, new_view_number)

    generic_single_entry("View Name: ", local_view_name(), " Add View ",
                         lambda text: add_view_local_func(text))


def add_view_to_views_panel(view_name, view_number):
    global views_dialog_vbox
    if (views_dialog_vbox):
        button = Gtk.Button(view_name)
        button.connect(
            "clicked", lambda func: coot.go_to_view_number(view_number, 0))
        views_dialog_vbox.pack_start(button, False, False, 2)
        button.show()

# return a list of [h_box_buttons, window]
#
# a button is a list of [label, callback, text_description]
#


def dialog_box_of_buttons(window_name, geometry, buttons, close_button_label, post_hook=None):

    return dialog_box_of_buttons_with_check_button(window_name, geometry,
                                                   buttons, close_button_label,
                                                   False, False, False, post_hook)

# geometry is an improper list of ints
#
# return a list of [h_box_buttons, window]
#
# a button is a list of [label, callback, (optional: text_description)]
# where callback is a function that takes one argument (I think).
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
        textbuffer.create_tag("tag", foreground="black", background="#c0e6c0")
        textbuffer.insert_with_tags_by_name(start, description, "tag")

    def close_cb_func(*args):
        if post_close_hook:
            post_close_hook()
        window.hide()

    def hide_this_window(widget):
        gtk_widget_hide(widget)

    # main line
    window = Gtk.Window()
    scrolled_win = Gtk.ScrolledWindow()
    outside_vbox = Gtk.VBox(False, 2)
    h_sep = Gtk.HSeparator()
    inside_vbox = Gtk.VBox(False, 0)

    window.set_default_size(geometry[0], geometry[1])
    window.set_title(window_name)
    inside_vbox.set_border_width(2)
    window.add(outside_vbox)

    if check_button_label:
        check_button = Gtk.CheckButton(check_button_label)
        # somehow need to execute the function before we can use it in the
        # callback. This is odd to say the least. FIXME
        check_button_func(check_button, inside_vbox)
        check_button.connect("toggled", lambda func:
                             check_button_func(check_button, inside_vbox))
        if check_button_is_initially_on_flag:
            check_button.set_active(True)
        outside_vbox.pack_start(check_button, False, False, 2)

    outside_vbox.pack_start(scrolled_win, True, True, 0)  # expand fill padding
    scrolled_win.add_with_viewport(inside_vbox)

    for button_info in buttons:
        # print("DEBUG:: in dialog_box_of_buttons_with_check_button(): button_info:", button_info)
        add_button_info_to_box_of_buttons_vbox(button_info, inside_vbox)

    outside_vbox.set_border_width(2)
    outside_vbox.pack_start(h_sep, False, False, 2)
    ok_button = Gtk.Button(close_button_label)
    outside_vbox.pack_end(ok_button, False, False, 0)
    ok_button.connect("clicked", close_cb_func, window, post_close_hook)
    window.connect("delete-event", lambda widget : hide_this_window())

    window.show_all()
    return [inside_vbox, window]

def dialog_box_of_buttons_from_specs(window_name, geometry,
                                     imol, specs):

   # not needed, why here Paul?
   #func = [cmd2str(set_go_to_atom_molecule, imol),
   #        cmd2str(set_go_to_atom_chain_residue_atom_name,
   #                chain_id, res_no, atom_name)]

   buttons = []
   for spec in specs:
      label = residue_spec_to_string(spec)
      cbf = [cmd2str(set_go_to_atom_molecule, imol),
             cmd2str(set_go_to_atom_chain_residue_atom_name,
                     residue_spec_to_chain_id(spec),
                     residue_spec_to_res_no(spec), " C  ")]
      buttons.append([label, cbf])

   return dialog_box_of_buttons(window_name, geometry, buttons, " Close ")

# This is exported outside of the box-of-buttons gui because the
# clear_and_add_back function (e.g. from using the check button)
# needs to add buttons - let's not duplicate that code.
#


def add_button_info_to_box_of_buttons_vbox(button_info, vbox):

    def add_text_to_text_buffer(text_buffer, description):
        start = text_buffer.get_start_iter()
        text_buffer.create_tag("tag", foreground="black",
                               background="#c0e6c0")
        text_buffer.insert_with_tags_by_name(start, description, "tag")

    # main line

    # print("debug:: add_button_info_to_box_of_buttons_vbox() button_info:", button_info)

    button_label = button_info[0]
    # print("debug in add_button_info_to_box_of_buttons_vbox with button bits", button_label, button_info[1])
    if ((button_label == "HSep") and (len(button_info) == 1)):
        # insert a HSeparator rather than a button
        button = Gtk.HSeparator()
    else:
        callback = button_info[1]
        # print("debug:: in add_button_info_to_box_of_buttons_vbox, callback is", callback)
        if (len(button_info) == 2):
            description = False
        else:
            description = button_info[2]
        button = Gtk.Button(button_label)

        # 20210827-PE functions must be real functions, not strings
        button.connect("clicked", callback)

        if description:
            text_view = Gtk.TextView()
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

    window = Gtk.Window()
    scrolled_win = Gtk.ScrolledWindow()
    outside_vbox = Gtk.VBox(False, 2)
    inside_vbox = Gtk.VBox(False, 0)

    window.set_default_size(geometry[0], geometry[1])
    window.set_title(window_name)
    inside_vbox.set_border_width(2)
    window.add(outside_vbox)
    outside_vbox.pack_start(scrolled_win, True, True, 0)  # expand fill padding
    scrolled_win.add_with_viewport(inside_vbox)

    for button_info in buttons:
        # print "button_info ", button_info
        if type(button_info) is ListType:
            button_label_1 = button_info[0][0]
            callback_1 = button_info[0][1]

            button_label_2 = button_info[1][0]
            callback_2 = button_info[1][1]

            button_1 = Gtk.Button(button_label_1)
            h_box = Gtk.HBox(False, 2)

            # print "button_label_1 ", button_label_1
            # print "callback_1 ", callback_1
            # print "button_label_2 ", button_label_2
            # print "callback_2 ", callback_2

            def callback_func(button, call):
                eval(call)
            button_1.connect("clicked", callback_func, callback_1)
            h_box.pack_start(button_1, False, False, 2)

            if callback_2:
                button_2 = Gtk.Button(button_label_2)
                button_2.connect("clicked", callback_func, callback_2)
                h_box.pack_start(button_2, False, False, 2)
            inside_vbox.pack_start(h_box, False, False, 2)

    outside_vbox.set_border_width(2)
    ok_button = Gtk.Button(close_button_label)
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
                              background="#c0e6c0")
        textbuffer.insert_with_tags_by_name(start, description, "tag")

    # main line
    window = Gtk.Window()
    scrolled_win = Gtk.ScrolledWindow()
    outside_vbox = Gtk.VBox(False, 2)
    inside_vbox = Gtk.VBox(False, 0)
    h_sep = Gtk.HSeparator()

    window.set_default_size(geometry[0], geometry[1])
    window.set_title(window_name)
    inside_vbox.set_border_width(2)
    window.add(outside_vbox)
    outside_vbox.pack_start(scrolled_win, True, True, 0)  # expand fill padding
    scrolled_win.add_with_viewport(inside_vbox)

    for button_info in buttons:
        button_label = button_info[0]
        callback = button_info[1]
        if len(button_info) == 2:
            description = False
        else:
            description = button_info[2]
        button = Gtk.Button(button_label)

# BL says:: in python we should pass the callback as a string
        if type(callback) is StringType:
            def callback_func(button, call):
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
            text_box = Gtk.TextView()
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
    ok_button = Gtk.Button(close_button_label)
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
                               selected_button=0,
                               cancel_button_label="",
                               cancel_function=False):

    def go_function_event(widget, button_ls):
        eval(go_button_function)
        window.destroy()
        return False

    def cancel_function_cb(widget):
        # print "BL DEBUG:: just for the sake"
        if (cancel_function):
            eval(cancel_function)
        window.destroy()
        return False

    # main line
    window = Gtk.Window()
    scrolled_win = Gtk.ScrolledWindow()
    outside_vbox = Gtk.VBox(False, 2)
    inside_vbox = Gtk.VBox(False, 0)
    button_hbox = Gtk.HBox(False, 0)

    window.set_default_size(geometry[0], geometry[1])
    window.set_title(window_name)
    inside_vbox.set_border_width(2)
    window.add(outside_vbox)
    outside_vbox.pack_start(scrolled_win, True, True, 0)  # expand fill padding
    scrolled_win.add_with_viewport(inside_vbox)

    button = None
    button_ls = []
    for button_info in buttons:
        button_label = button_info[0]
        callback = button_info[1]
        button = Gtk.RadioButton(button, button_label)

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
    go_button = Gtk.Button(go_button_label)
    outside_vbox.pack_start(button_hbox, False, False, 2)
    button_hbox.pack_start(go_button, True, True, 6)
    go_button.connect("clicked", go_function_event, button_ls)
    if (cancel_button_label):
        cancel_button = Gtk.Button(cancel_button_label)
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
    number_of_views = coot.n_views()
    buttons = []

    def generator(button_number):
        func = lambda button_number_c=button_number : coot.go_to_view_number(button_number_c, 0)
        def action(arg):
            print("DEBUG:: in views_panel_gui(): action() arg is", arg, "button_number is", button_number)
            func(button_number)
        return action

    for button_number in range(number_of_views):
        button_label = coot.view_name_py(button_number)
        description = coot.view_description_py(button_number)
        # BL says:: add the description condition
        # func = "go_to_view_number(" + str(button_number) + ",0)" 20211115-PE no to string functions
        func = generator(button_number)
        buttons.append([button_label, func, description])

    if len(buttons) > 1:
        def view_button_func():
            import time
            coot.go_to_first_view(1)
            time.sleep(1)
            coot.play_views()
        view_button = ["  Play Views ", lambda func: view_button_func()]
        buttons.insert(0, view_button)

    views_vbox = dialog_box_of_buttons("Views", [200, 140], buttons, "  Close  ")
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
        nudge = zsc * coot.zoom_factor()
        rc = coot_utils.rotation_centre()
        if operand == 0:
            rc[axis] += nudge
        elif operand == 1:
            rc[axis] -= nudge
        else:
            # We should never be here
            print("BL WARNING:: something went wrong!!!")
        coot.set_rotation_centre(*rc)

    buttons = [
        ["Nudge +X", lambda func: nudge_screen_func(0, 0)],
        ["Nudge -X", lambda func: nudge_screen_func(0, 1)],
        ["Nudge +Y", lambda func: nudge_screen_func(1, 0)],
        ["Nudge -Y", lambda func: nudge_screen_func(1, 1)],
        ["Nudge +Z", lambda func: nudge_screen_func(2, 0)],
        ["Nudge -Z", lambda func: nudge_screen_func(2, 1)],
    ]

    dialog_box_of_buttons("Nudge Screen Centre", [200, 240], buttons, "  Close ")


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
        nudge_factor = zsc * coot.zoom_factor()
        rc = coot_utils.rotation_centre()
        mat = coot_utils.view_matrix_transp()
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
        coot.set_rotation_centre(*rc)

    buttons = [
        ["Nudge +X", lambda func: nudge_screen_func(0, 0)],
        ["Nudge -X", lambda func: nudge_screen_func(0, 1)],
        ["Nudge +Y", lambda func: nudge_screen_func(1, 0)],
        ["Nudge -Y", lambda func: nudge_screen_func(1, 1)],
        ["Nudge +Z", lambda func: nudge_screen_func(2, 0)],
        ["Nudge -Z", lambda func: nudge_screen_func(2, 1)],
    ]

    dialog_box_of_buttons("Nudge Screen Centre", [200, 240], buttons, "  Close ")

# as nudge_screen_centre_gui but with clipping and zoom control


def nudge_screen_centre_extra_gui():

    # this is for the centre nudging
    zsc = 0.02
# BL says:: in python we need some helper functn, no 'proper' lambda
# axis 0-2 = x,y,z and 0=+ and 1=-

    def nudge_screen_func(axis, operand):
        nudge_factor = zsc * coot.zoom_factor()
        rc = coot_utils.rotation_centre()
        mat = coot_utils.view_matrix_transp()
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
        coot.set_rotation_centre(*rc)

    buttons = [
        ["Nudge +X", lambda func: nudge_screen_func(0, 0)],
        ["Nudge -X", lambda func: nudge_screen_func(0, 1)],
        ["Nudge +Y", lambda func: nudge_screen_func(1, 0)],
        ["Nudge -Y", lambda func: nudge_screen_func(1, 1)],
        ["Nudge +Z", lambda func: nudge_screen_func(2, 0)],
        ["Nudge -Z", lambda func: nudge_screen_func(2, 1)],
    ]

    # and this is for the clipping and zooming
    vbox = Gtk.VBox(False, 0)

    def change_clipp(*args):
        coot.set_clipping_front(clipp_adj.value)
        coot.set_clipping_back(clipp_adj.value)

    def change_zoom(*args):
        coot.set_zoom(zoom_adj.value)
        coot.graphics_draw()

    # for clipping
    clipp_label = Gtk.Label("Clipping")
    clipp_adj = Gtk.Adjustment(0.0, -10.0, 20.0, 0.05, 4.0, 10.1)
    clipp_scale = Gtk.HScale(clipp_adj)
    vbox.pack_start(clipp_label, False, False, 0)
    vbox.pack_start(clipp_scale, False, False, 0)
    clipp_label.show()
    clipp_scale.show()

    clipp_adj.connect("value_changed", change_clipp)

    h_sep = Gtk.HSeparator()
    vbox.pack_start(h_sep, False, False, 5)

    # for zooming
    zoom = coot.zoom_factor()
    zoom_label = Gtk.Label("Zoom")
    zoom_adj = Gtk.Adjustment(zoom, zoom*0.125, zoom*8, 0.01, 0.5, zoom)
    zoom_scale = Gtk.HScale(zoom_adj)
    vbox.pack_start(zoom_label, False, False, 0)
    vbox.pack_start(zoom_scale, False, False, 0)
    zoom_label.show()
    zoom_scale.show()

    zoom_adj.connect("value_changed", change_zoom)

    dialog_box_of_buttons_with_widget("Nudge Screen Centre with Extras", [200, 400], buttons, vbox, "  Close ")


# A gui to make a difference map (from arbitrarily gridded maps
# (that's it's advantage))
#
def make_difference_map_gui():

    def delete_event(*args):
        window.destroy()
        return False

    def go_function(widget):
        print("make diff map here\n")
        # active_mol_no_ref = get_option_menu_active_molecule(option_menu_ref_mol, coot_utils.map_molecule_list_ref)
        # active_mol_no_sec = get_option_menu_active_molecule(option_menu_sec_mol, coot_utils.map_molecule_list_sec)

        tree_iter_ref = ref_combobox.get_active_iter()
        tree_iter_sec = sec_combobox.get_active_iter()
        model_ref = ref_combobox.get_model()
        model_sec = sec_combobox.get_model()
        it_ref = model_ref[tree_iter_ref]
        it_sec = model_sec[tree_iter_sec]
        imol_1 = it_ref[0]
        imol_2 = it_sec[0]
        scale_text = scale_entry.get_text()
        scale = False
        try:
            scale = float(scale_text)
        except:
            print("WARNING:: make_difference_map_gui(): can't decode scale", scale_text)
        if scale:
            coot.difference_map(imol_1, imol_2, scale)
        delete_event()

    window = Gtk.Window()
    diff_map_vbox = Gtk.VBox(False, 2)
    h_sep = Gtk.HSeparator()
    title = Gtk.Label("Make a Difference Map")
    ref_label = Gtk.Label("Reference Map:")
    sec_label = Gtk.Label("Substract this map:")
    second_map_hbox = Gtk.HBox(False, 2)
    buttons_hbox = Gtk.HBox(True, 6)
    # option_menu_ref_mol = Gtk.combo_box_new_text()
    # option_menu_sec_mol = Gtk.combo_box_new_text()

    scale_label = Gtk.Label("Scale")
    scale_entry = Gtk.Entry()
    ok_button = Gtk.Button("   OK   ")
    cancel_button = Gtk.Button(" Cancel ")

    # map_molecule_list_ref = fill_option_menu_with_map_mol_options(option_menu_ref_mol)
    # map_molecule_list_sec = fill_option_menu_with_map_mol_options(option_menu_sec_mol)

    # ----------------- replacing option menu -----------------------

    ref_combobox_items = make_store_for_map_molecule_combobox()
    sec_combobox_items = make_store_for_map_molecule_combobox()
    ref_combobox = Gtk.ComboBox.new_with_model(ref_combobox_items)
    sec_combobox = Gtk.ComboBox.new_with_model(sec_combobox_items)
    if len(ref_combobox_items) > 0:
        ref_combobox.set_active(0)
    if len(sec_combobox_items) > 0:
        sec_combobox.set_active(0)
    ref_combobox.set_entry_text_column(1)
    sec_combobox.set_entry_text_column(1)
    ref_renderer_text = Gtk.CellRendererText()
    sec_renderer_text = Gtk.CellRendererText()
    ref_combobox.pack_start(ref_renderer_text, True)
    sec_combobox.pack_start(sec_renderer_text, True)
    ref_combobox.add_attribute(ref_renderer_text, "text", 1)
    sec_combobox.add_attribute(sec_renderer_text, "text", 1)

    # ---------------------------------------------------------------

    window.add(diff_map_vbox)
    diff_map_vbox.pack_start(title, False, False, 2)
    diff_map_vbox.pack_start(ref_label, False, False, 2)
    # diff_map_vbox.pack_start(option_menu_ref_mol, True, True, 2)
    diff_map_vbox.pack_start(ref_combobox, True, False, 2)

    diff_map_vbox.pack_start(sec_label, False, False, 2)
    diff_map_vbox.pack_start(second_map_hbox, False, False, 2)

    # second_map_hbox.pack_start(option_menu_sec_mol, True, True, 2)
    second_map_hbox.pack_start(sec_combobox, True, False, 2)
    second_map_hbox.pack_start(scale_label, False, False, 2)
    second_map_hbox.pack_start(scale_entry, False, False, 2)

    diff_map_vbox.pack_start(h_sep, True, False, 2)
    diff_map_vbox.pack_start(buttons_hbox, True, False, 2)
    buttons_hbox.pack_start(ok_button, True, False, 2)
    buttons_hbox.pack_start(cancel_button, True, False, 2)
    scale_entry.set_text("1.0")

    ok_button.connect("clicked", go_function)

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
            atom_list_r1 = coot.residue_info_py(imol, *r1[1:4])
            atom_list_r2 = coot.residue_info_py(imol, *r2[1:4])
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
                        coot.residue_name(imol, *r1[1:4]) + " - " +
                        str(r2[2]) + " " +
                        coot.residue_name(imol, *r2[1:4]) + "   " +
                        tors_string)
                ls = pos
                ls.insert(0, mess)
                ret.append(ls)
        return ret

    cis_peps = coot.cis_peptides_py(imol)

    if (cis_peps == []):
        coot.info_dialog("No Cis Peptides found")
    else:
        list_of_cis_peptides = make_list_of_cis_peps(imol, cis_peps)
        interesting_things_gui("Cis Peptides:", list_of_cis_peptides)

#


def transform_map_using_lsq_matrix_gui():

    def delete_event(*args):
        window.destroy()
        return False

    def on_ok_button_clicked(*args):
        # active_mol_ref = get_option_menu_active_molecule(frame_info_ref[1], frame_info_ref[2])
        # active_mol_mov = get_option_menu_active_molecule(frame_info_mov[1], frame_info_mov[2])

        combobox_ref = frame_info_ref[1]
        combobox_mov = frame_info_mov[1]

        tree_iter = combobox_ref.get_active_iter()
        if tree_iter is not None:
            model = combobox_ref.get_model()
            it = model[tree_iter]
            imol_ref = it[0]

        tree_iter = combobox_mov.get_active_iter()
        if tree_iter is not None:
            model = combobox_mov.get_model()
            it = model[tree_iter]
            imol_mov = it[0]

        print("debug:: imol_ref", imol_ref)
        print("debug:: imol_mov", imol_mov)

        chain_id_ref     = frame_info_ref[2].get_text()
        resno_1_ref_text = frame_info_ref[3].get_text()
        resno_2_ref_text = frame_info_ref[4].get_text()

        chain_id_mov     = frame_info_mov[2].get_text()
        resno_1_mov_text = frame_info_mov[3].get_text()
        resno_2_mov_text = frame_info_mov[4].get_text()

        radius_text = radius_entry.get_text()

        imol_map = coot.imol_refinement_map()
        cont = False
        try:
            resno_1_ref = int(resno_1_ref_text)
            resno_2_ref = int(resno_2_ref_text)
            resno_1_mov = int(resno_1_mov_text)
            resno_2_mov = int(resno_2_mov_text)
            radius = float(radius_text)

            if not coot_utils.valid_map_molecule_qm(imol_map):
                coot.info_dialog("WARNING:: Must set the refinement map")
            else:
                imol_copy = coot.copy_molecule(imol_mov)
                new_map_number = coot_utils.transform_map_using_lsq_matrix(imol_ref, chain_id_ref,
                                                                           resno_1_ref, resno_2_ref,
                                                                           imol_copy, chain_id_mov,
                                                                           resno_1_mov, resno_2_mov,
                                                                           imol_map, coot_utils.rotation_centre(), radius)
                sp = coot_utils.strip_path(coot.molecule_name(imol_mov))
                coot.set_molecule_name(imol_copy, "Transformed copy of " + sp)

                s = "Transformed map: from map " + str(imol_map) + " by matrix that created coords " + str(imol_copy)
                coot.set_molecule_name(new_map_number, s)
                coot.set_mol_active(imol_copy, 0)
                coot.set_mol_displayed(imol_copy, 0)

        except KeyError as e:
            coot.info_dialog("ERROR:: conversion from input text to numbers failed")

        window.destroy()

    # atom-sel-type is either 'Reference or 'Moving
    #
    # return the list [frame, option_menu, model_mol_list, entries...]
    def atom_sel_frame(atom_sel_type):
        frame = Gtk.Frame()

        # ----------------- replacing option_menu ------------------------

        # option_menu = Gtk.combo_box_new_text()
        # model_mol_list = fill_option_menu_with_coordinates_mol_options(option_menu)

        combobox_items = make_store_for_model_molecule_combobox()
        combobox = Gtk.ComboBox.new_with_model(combobox_items)
        renderer_text = Gtk.CellRendererText()
        if len(combobox_items)> 0:
            combobox.set_active(0)
        combobox.set_entry_text_column(1)
        combobox.pack_start(renderer_text, True)
        combobox.add_attribute(renderer_text, "text", 1)

        atom_sel_vbox = Gtk.VBox(False, 2)
        atom_sel_hbox = Gtk.HBox(False, 2)
        chain_id_label = Gtk.Label(" Chain ID ")
        resno_1_label = Gtk.Label(" Resno Start ")
        resno_2_label = Gtk.Label(" Resno End ")
        chain_id_entry = Gtk.Entry()
        resno_1_entry = Gtk.Entry()
        resno_2_entry = Gtk.Entry()

        frame.add(atom_sel_vbox)
        atom_sel_vbox.pack_start(combobox,       False, False, 2)
        atom_sel_vbox.pack_start(atom_sel_hbox,  False, False, 2)
        atom_sel_hbox.pack_start(chain_id_label, False, False, 2)
        atom_sel_hbox.pack_start(chain_id_entry, False, False, 2)
        atom_sel_hbox.pack_start(resno_1_label,  False, False, 2)
        atom_sel_hbox.pack_start(resno_1_entry,  False, False, 2)
        atom_sel_hbox.pack_start(resno_2_label,  False, False, 2)
        atom_sel_hbox.pack_start(resno_2_entry,  False, False, 2)

        return [frame, combobox, chain_id_entry, resno_1_entry, resno_2_entry]

    window = Gtk.Window()
    dialog_name = "Map Transformation"
    main_vbox = Gtk.VBox(False, 2)
    buttons_hbox = Gtk.HBox(False, 2)
    cancel_button = Gtk.Button("  Cancel  ")
    ok_button = Gtk.Button("  Transform  ")
    usage = "Note that this will transform the current refinement map " + \
            "to around the screen centre"
    usage_label = Gtk.Label(usage)
    h_sep = Gtk.HSeparator()
    frame_info_ref = atom_sel_frame("Reference")
    frame_info_mov = atom_sel_frame("Moving")
    radius_hbox = Gtk.HBox(False, 2)
    radius_label = Gtk.Label("  Radius ")
    radius_entry = Gtk.Entry()
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

    frame_info_ref[2].set_text("A")
    frame_info_ref[3].set_text("1")
    frame_info_ref[4].set_text("10")
    frame_info_mov[2].set_text("B")
    frame_info_mov[3].set_text("1")
    frame_info_mov[4].set_text("10")

    radius_entry.set_text("8.0")

    cancel_button.connect("clicked", delete_event)
    ok_button.connect("clicked", on_ok_button_clicked)

    window.show_all()
    if (not coot_utils.valid_map_molecule_qm(coot.imol_refinement_map())):
        coot.show_select_map_dialog()


def ncs_ligand_gui():

    def delete_event(*args):
        window.destroy()
        return False

    def go_button_function(widget):
        print("ncs ligand function here\n")
        active_mol_no_ref = get_option_menu_active_molecule(
            option_menu_ref_mol, molecule_list_ref)
        active_mol_no_lig = get_option_menu_active_molecule(
            option_menu_lig_mol, molecule_list_lig)
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
            coot.make_ncs_ghosts_maybe(active_mol_no_ref)
            print("ncs ligand with", active_mol_no_ref,
                  chain_id_ref, active_mol_no_lig, chain_id_lig,
                  resno_start, resno_end)
            ncs.ncs_ligand(active_mol_no_ref,
                           chain_id_ref,
                           active_mol_no_lig,
                           chain_id_lig,
                           resno_start,
                           resno_end)

        delete_event()

    window = Gtk.Window()
    ncs.ncs_ligands_vbox = Gtk.VBox(False, 2)
    title = Gtk.Label("Find NCS-Related Ligands")
    ref_label = Gtk.Label("Protein with NCS:")
    ref_chain_hbox = Gtk.HBox(False, 2)
    chain_id_ref_label = Gtk.Label("NCS Master Chain")
    chain_id_ref_entry = Gtk.Entry()
    lig_label = Gtk.Label("Molecule containing ligand")
    specs_hbox = Gtk.HBox(False, 2)
    h_sep = Gtk.HSeparator()
    buttons_hbox = Gtk.HBox(True, 6)
    chain_id_lig_label = Gtk.Label("Chain ID: ")
    resno_start_label = Gtk.Label(" Residue Number ")
    to_label = Gtk.Label("  to  ")
    chain_id_lig_entry = Gtk.Entry()
    resno_start_entry = Gtk.Entry()
    resno_end_entry = Gtk.Entry()
    ok_button = Gtk.Button("   Find Candidate Positions  ")
    cancel_button = Gtk.Button("    Cancel    ")
    option_menu_ref_mol = Gtk.combo_box_new_text()
    option_menu_lig_mol = Gtk.combo_box_new_text()

    molecule_list_ref = fill_option_menu_with_coordinates_mol_options(
        option_menu_ref_mol)
    molecule_list_lig = fill_option_menu_with_coordinates_mol_options(
        option_menu_lig_mol)

    window.add(ncs_ligands_vbox)
    ncs.ncs_ligands_vbox.pack_start(title, False, False, 6)
    ncs.ncs_ligands_vbox.pack_start(ref_label, False, False, 2)
    ncs.ncs_ligands_vbox.pack_start(option_menu_ref_mol, True, False, 2)
    ncs.ncs_ligands_vbox.pack_start(ref_chain_hbox, False, False, 2)
    ncs.ncs_ligands_vbox.pack_start(lig_label, False, False, 2)
    ncs.ncs_ligands_vbox.pack_start(option_menu_lig_mol, True, False, 2)
    ncs.ncs_ligands_vbox.pack_start(specs_hbox, False, False, 2)
    ncs.ncs_ligands_vbox.pack_start(h_sep, False, False, 2)
    ncs.ncs_ligands_vbox.pack_start(buttons_hbox, False, False, 2)

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
    specs_hbox.pack_start(Gtk.Label(" "), False, False, 2)  # neatness ?!

    chain_id_lig_entry.set_size_request(32, -1)
    chain_id_ref_entry.set_size_request(32, -1)
    resno_start_entry.set_size_request(50, -1)
    resno_end_entry.set_size_request(50, -1)
    chain_id_ref_entry.set_text("A")
    chain_id_lig_entry.set_text("A")
    resno_start_entry.set_text("1")

    #tooltips = Gtk.Tooltips()
    chain_tip = "'A' is a reasonable guess at the NCS master chain id.  " + \
                "If your ligand (specified below) is NOT bound to the protein's " + \
                "'A' chain, then you will need to change this chain and also " + \
                "make sure that the master molecule is specified appropriately " + \
                "in the Draw->NCS Ghost Control window."
    resno_tip = "Leave blank for a single residue"
    if Gtk.pygtk_version >= (2, 12):
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
        ncs.skip_to_next_ncs_chain("forward")
        if timeout_function_token:
            return True
        else:
            return False

    def start_function_event(*args):
        global timeout_function_token
        if not coot_utils.isNumber(timeout_function_token):
            timeout_function_token = gobject.timeout_add(ms_step,
                                                         skip_ncs_timeout_func)
        else:
            timeout_function_token = False

    def stop_function_event(*args):
        global timeout_function_token
        timeout_function_token = False

    # main body
    window = Gtk.Window()
    outside_vbox = Gtk.VBox(False, 2)
    inside_hbox = Gtk.HBox(False, 2)
    cancel_hbox = Gtk.HBox(False, 2)  # paul says VBox?!?!
    h_sep = Gtk.HSeparator()
    jump_start_button = Gtk.Button("NCS Jump Start")
    jump_stop_button = Gtk.Button("Stop")
    cancel_button = Gtk.Button("Cancel")
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

    def combobox_to_molecule_number(combobox):
      imol = -1
      tree_iter = combobox.get_active_iter()
      if tree_iter is not None:
        model = combobox.get_model()
        it = model[tree_iter]
        imol = it[0]
      return imol

    def go_button_function(widget):
        active_mol_no_ref_lig = combobox_to_molecule_number(combobox_ref)
        active_mol_no_mov_lig = combobox_to_molecule_number(combobox_mov)

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
            coot_utils.overlay_my_ligands(active_mol_no_mov_lig, chain_id_mov, resno_mov,
                                          active_mol_no_ref_lig, chain_id_ref, resno_ref)

        delete_event()

    window = Gtk.Window()
    title = Gtk.Label("Superpose Ligands")
    ligands_vbox = Gtk.VBox(False, 2)
    ref_chain_hbox = Gtk.HBox(False, 2)
    chain_id_ref_label = Gtk.Label("Ligand Chain ID: ")
    chain_id_ref_entry = Gtk.Entry()
    resno_ref_label = Gtk.Label(" Residue Number ")
    resno_ref_entry = Gtk.Entry()

    mov_chain_hbox = Gtk.HBox(False, 2)
    chain_id_mov_label = Gtk.Label("Ligand Chain ID: ")
    chain_id_mov_entry = Gtk.Entry()
    resno_mov_label = Gtk.Label(" Residue Number ")
    resno_mov_entry = Gtk.Entry()

    h_sep = Gtk.HSeparator()

    buttons_hbox = Gtk.HBox(True, 6)

    ok_button = Gtk.Button("   Superpose 'em  ")
    cancel_button = Gtk.Button("    Cancel    ")

    window.add(ligands_vbox)
    ligands_vbox.pack_start(title, False, False, 6)
    combobox_ref_mol = generic_molecule_chooser(ligands_vbox, "Model with reference ligand")
    ligands_vbox.pack_start(ref_chain_hbox, False, False, 2)

    combobox_mov_mol = generic_molecule_chooser(ligands_vbox, "Model with moving ligand")
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

#   tooltips = Gtk.Tooltips()
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
    res_name = coot.residue_name(imol_ref, chain_id_ref, res_no_ref, "")
    restraints = monomer_restraints(res_name)
    if (not restraints):
        return False
    else:
        if (not coot_utils.residue_has_hetatms_qm(imol_ref, chain_id_ref, res_no_ref, "")):
            return False
        else:
            print("----------- overlap-ligands %s %s %s %s ------------"
                  % (imol_ligand, imol_ref, chain_id_ref, res_no_ref))
            # this can return the rtop operator or the False (for fail of course).
            coot.match_ligand_torsions(imol_ligand, imol_ref,
                                  chain_id_ref, res_no_ref)
            ret = overlap_ligands(imol_ligand, imol_ref,
                                  chain_id_ref, res_no_ref)
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

        binding_hbox = Gtk.HBox(False, 2)
        txt = str(item[1])
        key_label = Gtk.Label("   " + txt + "   ")
        name_label = Gtk.Label(item[2])

        if (buttonize_flag):
            button_label = "   " + txt + "   " + item[2]
            button = Gtk.Button(button_label)
            #al = Gtk.Alignment(0, 0, 0, 0)
            #label = Gtk.Label(button_label)
            # button.add(al)
            # al.add(label)
            binding_hbox.pack_start(button, True, True, 0)
            inside_vbox.pack_start(binding_hbox, False, False, 0)
            binding_func = item[3]
            if not (callable(binding_func)):
                s = "Cannot call given function with button,\n"
                s += "probably a scheme function.\n"
                s += "The shortcut should still work though."

                def binding_func():
                    coot.info_dialog(s)
                    print("INFO::", s)

            button.connect("clicked", lambda func: apply(binding_func))

        else:
            binding_hbox.pack_start(key_label, False, False, 2)
            binding_hbox.pack_start(name_label, False, False, 2)
            inside_vbox.pack_start(binding_hbox, False, False, 2)

    # main line
    #
    window = Gtk.Window()
    scrolled_win = Gtk.ScrolledWindow()
    outside_vbox = Gtk.VBox(False, 2)
    inside_vbox = Gtk.VBox(False, 0)
    dialog_name = "Key Bindings"
    buttons_hbox = Gtk.HBox(False, 2)
    close_button = Gtk.Button("  Close  ")
    std_frame = Gtk.Frame()
    usr_frame = Gtk.Frame()
    std_frame_vbox = Gtk.VBox(False, 2)
    usr_frame_vbox = Gtk.VBox(False, 2)
    close_button.connect("clicked", delete_event)

    window.set_default_size(250, 350)
    window.set_title(dialog_name)
    inside_vbox.set_border_width(4)

    window.add(outside_vbox)
    outside_vbox.add(scrolled_win)
    scrolled_win.add_with_viewport(inside_vbox)

    inside_vbox.pack_start(std_frame, False, False, 2)
    inside_vbox.pack_start(usr_frame, False, False, 2)

    std_frame.add(std_frame_vbox)
    usr_frame.add(usr_frame_vbox)

    py_and_scm_keybindings = key_bindings
    if (coot.coot_has_guile()):
        scm_key_bindings = coot.run_scheme_command("*key-bindings*")
        # check for list
        if isinstance(scm_key_bindings, list):
            # filter out doublicates
            for item in scm_key_bindings:
                scm_code, scm_key, text, tmp = item
                py_keys = [elem[1] for elem in py_and_scm_keybindings]
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
        active_mol_no = get_option_menu_active_molecule(
            option_menu, model_mol_list)
        imol = int(active_mol_no)
        imol_map = coot.imol_refinement_map()

        do_it = assign_sequences_to_mol(imol)

        if (do_it):
            # now cootaneer it
            chain_ls = coot_utils.chain_ids(imol)
            for chain_id in chain_ls:
                res_name = coot.resname_from_serial_number(imol, chain_id, 0)
                res_no = coot.seqnum_from_serial_number(imol, chain_id, 0)
                ins_code = coot.insertion_code_from_serial_number(imol, chain_id, 0)
                alt_conf = ""
                at_name = coot_utils.residue_spec_to_atom_for_centre(
                    imol, chain_id, res_no, ins_code)[0]
                cootaneer_result = cootaneer(imol_map, imol, [chain_id, res_no, ins_code,
                                                              at_name, alt_conf])
                if (cootaneer_result == 0):
                    s = "Insufficiently confident in alignment to make a fit." + \
                        "\n" + \
                        "Perhaps you could improve or extend this fragment."
                    coot.info_dialog(s)
                refine_qm = refine_check_button.get_active()
            # refine?
            window.hide()
            if (refine_qm):
                fitting.fit_protein(imol)

            delete_event()

    def import_function_event(widget, selector_entry):
        # we import a sequence file and update the cootaneer table
        global imported_sequence_file_flags
        imported_sequence_file_qm = imported_sequence_file_flags[0]
        active_mol_no = get_option_menu_active_molecule(
            option_menu, model_mol_list)
        imol = int(active_mol_no)

        seq_info_ls = []
        seq_file_name = selector_entry.get_text()
        if (seq_file_name):
            # get and set sequence info
            coot.assign_sequence_from_file(imol, str(seq_file_name))
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
        active_mol_no = get_option_menu_active_molecule(
            option_menu, model_mol_list)
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
                               background="#c0e6c0")
        text_buffer.insert_with_tags_by_name(start, description, "tag")

    # return the (entry . textbuffer/box)
    #
    def entry_text_pair_frame_with_button(seq_info):

        def fragment_go_event(widget):
            active_mol_no = get_option_menu_active_molecule(
                option_menu, model_mol_list)
            imol = int(active_mol_no)
            imol_map = coot.imol_refinement_map()
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
                    imol = n_atom[0]
                    chain_id = n_atom[1]
                    res_no = n_atom[2]
                    ins_code = n_atom[3]
                    at_name = n_atom[4]
                    alt_conf = n_atom[5]
                    cootaneer_result = cootaneer(imol_map, imol, [chain_id, res_no, ins_code,
                                                                  at_name, alt_conf])
                    if (cootaneer_result == 0):
                        s = "Insufficiently confident in alignment to make a fit." + \
                            "\n" + \
                            "Perhaps you could improve or extend this fragment."
                        coot.info_dialog(s)
                    else:
                        refine_qm = refine_check_button.get_active()
                        if (chain_check_button.get_active()):
                            # we try to change the chain
                            from_chain_id = chain_id
                            to_chain_id = entry.get_text()
                            start_resno = coot.seqnum_from_serial_number(
                                imol, from_chain_id, 0)
                            end_resno = coot.seqnum_from_serial_number(
                                imol, from_chain_id, (coot.chain_n_residues(from_chain_id, imol) - 1))
                            [istat, message] = coot.change_chain_id_with_result_py(
                                imol, from_chain_id, to_chain_id, 1, start_resno, end_resno)
                            if (istat == 1):
                                # renaming ok
                                chain_id = to_chain_id
                            else:
                                # renaming failed
                                if (refine_qm):
                                    message += "\nRefinement proceeded with old=new chain "
                                    message += chain_id
                                coot.info_dialog(message)

                        if (refine_qm):
                            fitting.fit_chain(imol, chain_id)
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

        frame = Gtk.Frame()
        vbox = Gtk.VBox(False, 3)
        hbox = Gtk.HBox(False, 3)
        entry = Gtk.Entry()
        textview = Gtk.TextView()
        textview.set_wrap_mode(Gtk.WRAP_WORD_CHAR)
        textview.set_editable(True)
        textview.set_size_request(300, -1)
        textview.modify_font(pango.FontDescription("Courier 11"))
        text_buffer = textview.get_buffer()
        chain_id_label = Gtk.Label("Chain ID")
        sequence_label = Gtk.Label("Sequence")
        vbox_for_buttons = Gtk.VBox(False, 3)
        fragment_button = Gtk.Button("  Sequence closest fragment  ")
        chain_check_button = Gtk.CheckButton("Assign Chain ID as well?")

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
                    child_range = list(
                        range(0, no_of_children - no_of_sequences))
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
    def clear_function_event(widget=None, file_sel_entry=None):
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
                            coot.delete_sequence_by_chain_id(imol, chain_id_old)

                # replace space, eol etc in sequence first
                seq = seq.replace(" ", "")
                seq = seq.replace("\r", "")       # mac?
                seq = seq.replace("\r\n", "")     # win?
                seq = seq.replace("\n", "")       # unix?
                coot.assign_sequence_from_string(imol, chain_id_new, seq)
            return True
        else:
            coot.add_status_bar_text("Invalid chain_id and/or sequence provided")
            return False

    # main body
    imol_map = coot.imol_refinement_map()
    if (imol_map == -1):
        coot.show_select_map_dialog()

    window = Gtk.Window()
    window.set_title("Sequencing GUI")
    #tooltips = Gtk.Tooltips()
    label = Gtk.Label("Molecule to be sequenced")
    vbox = Gtk.VBox(False, 2)
    option_menu = Gtk.combo_box_new_text()
    model_mol_list = fill_option_menu_with_mol_options(
        option_menu, coot_utils.valid_model_molecule_qm)
    inside_vbox = Gtk.VBox(False, 2)
    seq_table = Gtk.Table(1, 1, True)
    hbox_for_spin = Gtk.HBox(False, 0)
    spin_label = Gtk.Label("Number of Sequences:")
    spin_adj = Gtk.Adjustment(1, 1, 10, 1, 4, 0)
    spin_button = Gtk.SpinButton(spin_adj, 0, 0)
    refine_check_button = Gtk.CheckButton("Auto-fit-refine after sequencing?")
    h_sep = Gtk.HSeparator()
    h_sep2 = Gtk.HSeparator()
    buttons_hbox = Gtk.HBox(False, 2)
    import_button = Gtk.Button("  Import and associate sequence from file  ")
    go_button = Gtk.Button("  Sequence all fragments!  ")
    if Gtk.pygtk_version >= (2, 12):
        go_button.set_tooltip_text("This currently ignores all chain IDs")
    else:
        coot_tooltips.set_tip(
            go_button, "This currently ignores all chain IDs")
    cancel_button = Gtk.Button("  Cancel  ")
    clear_button = Gtk.Button("  Clear all  ")

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
    file_sel_entry = file_chooser_entry(window, vbox, "Select PIR file")
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


def generic_check_button(vbox, label_text, handle_check_button_function):
    def check_callback(*args):
        active_state = check_button.get_active()
        set_state = 0
        if (active_state):
            set_state = 1
        handle_check_button_function(set_state)
    check_button = Gtk.CheckButton(label_text)
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
            coot.set_matrix(t)
            coot.add_status_bar_text(s)
        except:
            coot.add_status_bar_text("Failed to read a number")

    def delete_event(*rags):
        window.destroy()
        return False

    def go_function_event(*args):
        set_matrix_func('dummy', matrix_entry)
        window.destroy()
        return False

    window = Gtk.Window()
    vbox = Gtk.VBox(False, 0)
    hbox = Gtk.HBox(False, 0)
    h_sep = Gtk.HSeparator()
    h_sep2 = Gtk.HSeparator()
    h_sep3 = Gtk.HSeparator()
    go_button = Gtk.Button("   Ok   ")
    cancel_button = Gtk.Button("  Cancel  ")

    window.add(vbox)
    # add the matrix entry
    matrix_entry = entry_do_button(vbox, "set matrix: (smaller means better geometry)",
                                   "Set", set_matrix_func)
    matrix_entry.set_text(str(coot.matrix_state()))

    vbox.pack_start(h_sep2, False, False, 2)

    # use torsion restrains?
    torsion_restraints_button = generic_check_button(vbox,
                                                     "Use Torsion Restraints?",
                                                     lambda state: coot.set_refine_with_torsion_restraints(state))
    if (coot.refine_with_torsion_restraints_state() == 1):
        torsion_restraints_button.set_active(True)
    else:
        torsion_restraints_button.set_active(False)

    # planar peptide restrains?
    planar_restraints_button = generic_check_button(vbox,
                                                    "Use Planar Peptide Restraints?",
                                                    lambda state: coot.remove_planar_peptide_restraints() if state == 0 else coot.add_planar_peptide_restraints())
    if (coot.planar_peptide_restraints_state() == 1):
        planar_restraints_button.set_active(True)
    else:
        planar_restraints_button.set_active(False)

    # use ramachandran restrains?
    rama_restraints_button = generic_check_button(vbox,
                                                  "Use Ramachandran Restraints?",
                                                  lambda state: coot.set_refine_ramachandran_angles(state))
    if (coot.refine_ramachandran_angles_state() == 1):
        rama_restraints_button.set_active(True)
    else:
        rama_restraints_button.set_active(False)

    vbox.pack_start(h_sep3, False, False, 2)

    # add rotamer check button
    rotamer_autofit_button = generic_check_button(vbox,
                                                  "Real Space Refine after Auto-fit Rotamer?",
                                                  lambda state: coot.set_rotamer_auto_fit_do_post_refine(state))
    if (coot.rotamer_auto_fit_do_post_refine_state() == 1):
        rotamer_autofit_button.set_active(True)
    else:
        rotamer_autofit_button.set_active(False)

    # add mutate check button
    mutate_autofit_button = generic_check_button(vbox,
                                                 "Real Space Refine after Mutate and Auto-fit?",
                                                 lambda state: coot.set_mutate_auto_fit_do_post_refine(state))
    if (coot.mutate_auto_fit_do_post_refine_state() == 1):
        mutate_autofit_button.set_active(True)
    else:
        mutate_autofit_button.set_active(False)

    # add terminal residue check button
    terminal_autofit_button = generic_check_button(vbox,
                                                   "Real Space Refine after Add Terminal Residue?",
                                                   lambda state: coot.set_add_terminal_residue_do_post_refine(state))
    if (coot.add_terminal_residue_do_post_refine_state() == 1):
        terminal_autofit_button.set_active(True)
    else:
        terminal_autofit_button.set_active(False)

    # add a b-factor button
    reset_b_factor_button = generic_check_button(vbox,
                                                 "Reset B-Factors to Default Value after Refinement/Move?",
                                                 lambda state: coot.set_reset_b_factor_moved_atoms(state))
    if (coot.get_reset_b_factor_moved_atoms_state() == 1):
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


def map_sharpening_gui(imol):

    window = Gtk.Window()
    vbox = Gtk.VBox(False, 2)
    hbox = Gtk.HBox(False, 2)
    adj = Gtk.Adjustment(0.0, -30, 60, 0.05, 2, 30.1)
    slider = Gtk.HScale(adj)
    label = Gtk.Label("\nSharpen Map:")
    lab2 = Gtk.Label("Add B-factor: ")

    vbox.pack_start(label,  False, False, 2)
    vbox.pack_start(hbox,   False, False, 2)
    hbox.pack_start(lab2,   False, False, 2)
    hbox.pack_start(slider, True,  True,  2)
    window.add(vbox)
    window.set_size_request(500, 100)
    # slider.add_mark(-30, -30, 0) # not yet, needs updated pygtk

    adj.connect("value_changed", lambda func: coot.sharpen(imol, adj.value))

    window.show_all()


# Associate the contents of a sequence file with a chain.
# Select file from a GUI.
# File format can be Fasta (default) or PIR
#
# added checkbox to assign for all protein chains
#
def associate_sequence_with_chain_gui(sequence_format="FASTA", do_alignment=False):

    def associate_func(imol, chain_id_in, pir_file_name,
                       use_all_chains=False):
        # print "assoc seq:", imol, chain_id, pir_file_name
        all_chains = chain_id_in
        if use_all_chains:
            all_chains = [chain for chain in coot_utils.chain_ids(
                imol) if coot_utils.is_protein_chain_qm(imol, chain)]
        for chain_id in all_chains:
            if (sequence_format == "FASTA"):
                coot_utils.associate_fasta_file(imol, chain_id, pir_file_name)
            elif (sequence_format == "PIR"):
                coot_utils.associate_pir_file(imol, chain_id, pir_file_name)
            else:
                coot.info_dialog("BL INFO:: wrong sequence input format.")
                return
        if do_alignment:
            alignment_mismatches_gui(imol)

    generic_chooser_entry_and_file_selector(
        "Associate Sequence with Chain: ",
        coot_utils.valid_model_molecule_qm,
        "Chain ID",
        "",
        "Select " + sequence_format + " file",
        lambda imol, chain_id, pir_file_name:
        associate_func(imol, chain_id, pir_file_name),
        use_check_button=True,
        check_button_label="Assign sequence to all protein chains (ignoring the chain specified above)",
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
         if (len(am) > 3):
            # we have alignment text for info box
            # protected for compatibiity reasons
            alignments_as_text_list = am[3]
            for alignment_text in alignments_as_text_list:
               info_dialog_with_markup(alignment_text)

         dialog_box_of_buttons("Residue mismatches", [300, 300], buttons, "  Close  ")

# Wrapper in that we test if there have been sequence(s) assigned to
# imol before we look for the sequence mismatches
#


def wrapper_alignment_mismatches_gui(imol):

    seq_info = coot.sequence_info_py(imol)
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

        frame = Gtk.Frame()
        outside_hbox = Gtk.HBox(False, 2)
        hbox = Gtk.HBox(False, 2)
        text_1 = Gtk.Label("  Chain-ID:")
        text_2 = Gtk.Label("  Resno Start:")
        text_3 = Gtk.Label("  Resno End:")
        entry_1 = Gtk.Entry()
        entry_2 = Gtk.Entry()
        entry_3 = Gtk.Entry()
        plus_button = Gtk.Button("+")
        minus_button = Gtk.Button(" - ")

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
            print("did not understand %s and %s as numbers - fail residue range" %
                  (res_no_1_txt, res_no_2_txt))
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
        resno_1 = range_info[1]
        resno_2 = range_info[2]

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
                    residue_range_vbox.pack_start(
                        vbox_info[0], False, False, 2)
                    print("next one", rang)
                    residue_range_widgets.append(vbox_info)
                    fill_with_previous_range(rang, vbox_info)

    # main line
    #
    def cancel_button_cb(*args):
        save_ranges(residue_range_widgets)
        window.destroy()
        return False

    def combobox_to_molecule_number(combobox):
        imol = -1
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol = it[0]
        return imol

    def go_button_cb(*args):

        save_ranges(residue_range_widgets)
        residue_ranges = make_residue_ranges()
        imol = combobox_to_molecule_number(combobox)
        try:
            func(imol, residue_ranges)
        except KeyError as e:
            print("ERROR:: something was amiss")
            print("ERROR::", e)
        window.destroy()
        return False

    window = Gtk.Window()
    vbox = Gtk.VBox(False, 0)
    residue_range_vbox = Gtk.VBox(True, 2)
    residue_range_widget_info = make_residue_range_frame(residue_range_vbox)
    hbox_buttons = Gtk.HBox(False, 0)
    function_label = Gtk.Label(function_text)
    cancel_button = Gtk.Button("  Cancel  ")
    go_button = Gtk.Button(go_button_label)
    h_sep = Gtk.HSeparator()
    # the first residue range
    outside_vbox_residue_range = residue_range_widget_info[0]

    residue_range_widgets = [residue_range_widget_info]  # ?

    # buttons
    hbox_buttons.pack_start(cancel_button, False, False, 6)
    hbox_buttons.pack_start(go_button, False, False, 6)

    # the vbox of residue ranges
    residue_range_vbox.pack_start(outside_vbox_residue_range, False, False, 0)

    if saved_residue_ranges:
        fill_residue_range_widgets_previous_data(saved_residue_ranges,
                                                 residue_range_widget_info,
                                                 residue_range_vbox)

    # main vbox
    vbox.pack_start(function_label, False, False, 0)
    combobox = generic_molecule_chooser(vbox, "Molecule for Ranges:")
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
# Change the coot_utils.translation jiggle-factor to 1.0, so the ligand doesn't
# move so far and get sucked into protein density (this is just a
# temporary hack, it would be better to mask the enviroment of the
# ligand by the surrounding atoms of the molecule to which the ligand
# is added - that is much harder).
#


def solvent_ligands_gui():

    global random_jiggle_n_trials
    #

    def add_ligand_func(imol, tlc):
        print("Add a %s to molecule %s here" % (tlc, imol))
        imol_ligand = coot.get_monomer(tlc)
        if (coot_utils.valid_model_molecule_qm(imol_ligand)):
            # delete hydrogens from the ligand if the master molecule
            # does not have hydrogens.
            if coot_utils.valid_model_molecule_qm(imol):
                if not coot_utils.molecule_has_hydrogens(imol):
                    coot.delete_residue_hydrogens(imol_ligand, "A", 1, "", "")
            if (coot_utils.valid_map_molecule_qm(coot.imol_refinement_map())):
                print("========  jiggling!  ======== ")

                coot.merge_molecules_py([imol_ligand], imol)

                # We no longer do this: because we now mask the map by the neighbours
                # of the jiggling ligand and for that we need to use molecule imol
                #
                # coot.fit_to_map_by_random_jiggle(imol_ligand, "A", 1, "",
                #                             random_jiggle_n_trials, 1.0)
                # coot_utils.with_auto_accept([refine_zone, imol_ligand, "A", 1, 1, ""])

                # we presume that the active residue is the jiggling residue!
                # (we'd like to know what imol_ligand: "A" 1 moves to on merging
                # but we don't)
                #
                active_atom = coot.active_residue_py()
                aa_chain_id = active_atom[1]
                aa_res_no = active_atom[2]
                coot.fit_to_map_by_random_jiggle(imol, aa_chain_id, aa_res_no, "",
                                            random_jiggle_n_trials, 1.0)

                # if we use refine_residues, that will take note of residues
                # near this residue and make non-bonded contacts
                # (whereas refine_zone will not).
                #
                # coot_utils.with_auto_accept([refine_zone, imol, aa_chain_id, aa_res_no, 1, ""])
                coot_utils.with_auto_accept(
                    [coot.refine_residues_py, imol, [[aa_chain_id, aa_res_no, ""]]])

            else:
                print("======== not jiggling - no map ======== ")
            if coot_utils.valid_model_molecule_qm(imol):
                coot.set_mol_active(imol_ligand, 0)
                coot.set_mol_displayed(imol_ligand, 0)

    # add a button for a 3-letter-code to the scrolled vbox that runs
    # add-ligand-func when clicked.
    #
    def add_solvent_button(comp_id, button_label, inside_vbox, combobox):

        def combobox_to_molecule_number(combobox):
            imol = -1
            tree_iter = combobox.get_active_iter()
            if tree_iter is not None:
                model = combobox.get_model()
                it = model[tree_iter]
                imol = it[0]
            return imol

        def button_cb(*args):
            imol = combobox_to_molecule_number(combobox)
            add_ligand_func(imol, comp_id)

        button = Gtk.Button(button_label)
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
            add_solvent_button(txt, comp_id_to_button_label(txt), inside_vbox, combobox)
        generic_single_entry("Add new 3-letter-code/comp-id", "", "  Add  ",
                             lambda txt:
                             add_button_func(txt))

    def comp_id_to_button_label(comp_id):
        coot.auto_load_dictionary(comp_id)
        comp_id_name = coot.comp_id_to_name_py(comp_id)
        label = (comp_id + ": " + comp_id_name) if comp_id_name else comp_id
        return label

    # main
    window = Gtk.Window()
    scrolled_win = Gtk.ScrolledWindow()
    outside_vbox = Gtk.VBox(False, 2)
    inside_vbox = Gtk.VBox(False, 2)
    label = Gtk.Label("\nSolvent molecules added to molecule: ")
    frame_for_option_menu = Gtk.Frame()
    vbox_for_combobox = Gtk.VBox(False, 6)
    # molecule_option_menu = Gtk.combo_box_new_text()

    combobox_items = make_store_for_model_molecule_combobox()
    combobox = Gtk.ComboBox.new_with_model(combobox_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_items) > 0:
        combobox.set_active(0)
    combobox.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)

    add_new_button = Gtk.Button("  Add a new Residue Type...")
    h_sep = Gtk.HSeparator()
    close_button = Gtk.Button("  Close  ")

    window.set_default_size(250, 500)
    window.set_title("Solvent Ligands")
    window.set_border_width(8)
    window.add(outside_vbox)
    outside_vbox.pack_start(label, False, False, 2)
    frame_for_option_menu.add(vbox_for_combobox)
    vbox_for_combobox.pack_start(combobox, False, False, 8)
    frame_for_option_menu.set_border_width(6)
    outside_vbox.pack_start(frame_for_option_menu, False, False, 2)
    outside_vbox.pack_start(scrolled_win, True, True, 0)
    scrolled_win.add_with_viewport(inside_vbox)
    # scrolled_win.set_policy(Gtk.POLICY_AUTOMATIC, Gtk.POLICY_ALWAYS)
    outside_vbox.pack_start(add_new_button, False, False, 6)
    outside_vbox.pack_start(h_sep, False, False, 2)
    outside_vbox.pack_start(close_button, False, False, 2)

    for comp_id in solvent_ligand_list():
        button_label = comp_id_to_button_label(comp_id)
        add_solvent_button(comp_id, button_label, inside_vbox, combobox)

    add_new_button.connect("clicked", add_new_button_cb)
    close_button.connect("clicked", delete_event)

    window.show_all()


# USER MODS gui
#
def user_mods_gui(imol, pdb_file_name):

    # no alt conf, no inscode
    def atom_spec_to_string(atom_spec):
        chain_id = atom_spec[1]
        res_no = atom_spec[2]
        ins_code = atom_spec[3]
        atom_name = atom_spec[4]
        alt_conf = atom_spec[5]
        return chain_id + str(res_no) + atom_name

    #
    def make_flip_buttons(flip_list):
        ret = []
        for flip in flip_list:
            atom_spec = flip[0]
            residue_type = flip[1]
            info_string = flip[2]
            set_string = flip[3]
            score = flip[4]
            chain_id = atom_spec[1]
            res_no = atom_spec[2]
            ins_code = atom_spec[3]
            atom_name = atom_spec[4]
            alt_conf = atom_spec[5]
            label = set_string + " " + \
                chain_id + " " + str(res_no) + \
                atom_name + " : " + \
                info_string + " " + \
                " score %2.2f" % score
            func = [cmd2str(coot.set_go_to_atom_molecule, imol),
                    coot_utils.cmd2str(set_go_to_atom_chain_residue_atom_name,
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
            chain_id = atom_spec[1]
            res_no = atom_spec[2]
            ins_code = atom_spec[3]
            atom_name = atom_spec[4]
            alt_conf = atom_spec[5]
            func = [cmd2str(coot.set_go_to_atom_molecule, imol),
                    coot_utils.cmd2str(set_go_to_atom_chain_residue_atom_name,
                                       chain_id, res_no, atom_name)]
            ret.append([label, func])
        return ret

    # Return a list of buttons that are (in this case, (only) clashes,
    # unknown/problems and flips.
    #
    def filter_buttons(flip_list):
        ret = []
        for flip in flip_list:
            atom_spec = flip[0]
            residue_type = flip[1]
            info_string = flip[2]
            set_string = flip[3]
            score = flip[4]
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
    if coot_utils.using_gui():
        flips = user_mods(pdb_file_name)
        flip_buttons = make_flip_buttons(flips[0])
        no_adj_buttons = make_no_adj_buttons(flips[1])
        all_buttons = no_adj_buttons + flip_buttons
        dialog_box_of_buttons_with_check_button(" Flips ", [300, 300], all_buttons, "  Close  ",
                                                "Clashes, Problems and Flips only",
                                                lambda check_button, vbox: clear_and_add_back(
                                                    vbox, flips[0], flips[1], True)
                                                if check_button.get_active() else
                                                clear_and_add_back(vbox, flips[0], flips[1], False),
                                                False)


def rename_residue_gui_simple():
    active_atom = coot.active_residue_py()
    if not active_atom:
        coot.info_dialog("INFO:: No Residue Here")
    else:
        print(active_atom)
        generic_single_entry("Rename this residue", "ALA", "Rename",
                             lambda text: coot_utils.using_active_atom(set_residue_name,
                                                                       "aa_imol", "aa_chain_id", "aa_res_no", "aa_ins_code",
                                                                       text))
#                           lambda text: coot_utils.using_active_atom([[set_residue_name,
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

            frame = Gtk.Frame()
            hbox = Gtk.HBox(False, 2)
            label = Gtk.Label("Weight:  ")
            entry = Gtk.Entry()
            combobox_items = make_store_for_molecule_combobox(coot.is_valid_map_molecule)
            combobox = Gtk.ComboBox.new_with_model(combobox_items)
            renderer_text = Gtk.CellRendererText()
            if len(combobox_items) > 0:
                combobox.set_active(0)
            combobox.set_entry_text_column(1) # Sets the model column which combo_box
                                              # should use to get strings from to be text_column
            combobox.pack_start(renderer_text, True)
            combobox.add_attribute(renderer_text, "text", 1)

            map_mol_list = coot_utils.map_molecule_list()
            plus_button = Gtk.Button("+")
            minus_button = Gtk.Button(" - ")
            hbox.pack_start(combobox, False, False, 2)
            hbox.pack_start(label, False, False, 2)
            hbox.pack_start(entry, False, False, 2)
            hbox.pack_start(plus_button, False, False, 2)
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
            list(map(lambda x: x.show(), [frame, hbox, label, entry, combobox, plus_button, minus_button]))

            # print "saving map_mol_list", map_mol_list
            return [hbox, combobox, map_mol_list, entry]

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
                    combobox = mav_bits[1]
                    map_mol_list = mav_bits[2]
                    entry = mav_bits[3]
                    # map_selected = get_option_menu_active_molecule(option_menu, map_mol_list)

                    map_selected = -1
                    tree_iter = combobox.get_active_iter()
                    if tree_iter is not None:
                        model = combobox.get_model()
                        it = model[tree_iter]
                        map_selected = it[0]

                    text = entry.get_text()
                    weight = float(text)   # try?! FIXME
                    maps_to_average_list.append([map_selected, weight])
                print("maps to average", maps_to_average_list)
                coot.average_map_py(maps_to_average_list)
                window.destroy()
                return False

            window = Gtk.Window()
            outer_vbox = Gtk.VBox(False, 0)
            inner_vbox = Gtk.VBox(False, 0)
            title = Gtk.Label("Average Maps")
            h_sep = Gtk.HSeparator()
            buttons_hbox = Gtk.HBox(False, 2)
            mav_widget = self.add_average_molecule_widget(inner_vbox)
            cancel_button = Gtk.Button("  Cancel  ")
            ok_button = Gtk.Button("  Average Maps  ")

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
    active_atom = coot.active_residue_py()
    if not active_atom:
        coot.info_dialog("No Residue Here")
    else:
        aa_imol     = active_atom[0]
        aa_chain_id = active_atom[1]
        aa_res_no   = active_atom[2]
        aa_ins_code = active_atom[3]
        label = "Rename Residue [in molecule " + str(aa_imol) + "]: " + \
                aa_chain_id + str(aa_res_no) + aa_ins_code + " to: "
        generic_single_entry(label, "ALA", "Rename Residue",
                             lambda text: coot.set_residue_name(aa_imol, aa_chain_id,
                                                           aa_res_no, aa_ins_code, text))


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
            print("now call update_water_results() with imol", imol, "coordination number", n, "imol", imol)
            update_water_results(imol, n, d)

    def key_press_event(widget, event):
        if (event.keyval == 65293):  # GDK_return
            n = get_number()
            imol = get_molecule()
            d = get_distance()
            if d:
                update_water_results(imol, n, d)

    def atom_spec_to_text(atom_spec):
        # remove the leading element (user_int, I think)
        # return " ".join(map(str, atom_spec))
        remains = atom_spec[1:]
        return " ".join(map(str, remains))

    def get_ele(imol, atom_spec):
        # print("debug:: get_ele() with imol", imol, " atom_spec", atom_spec)
        atoms = coot.residue_info_py(imol, atom_spec[1], atom_spec[2], atom_spec[3])
        # print("debug:: get_ele() with atoms", atoms)
        input_atom_name = atom_spec[4]
        try:
            for atom in atoms:
                if atom[0][0] == input_atom_name:
                    return atom[1][2]
        except KeyError as e:
            print(e)
        return False

    # add info about bumps (to non-H-bonding atoms or so).  given a
    # water info (a central atom spec and a list of its contacts).
    #
    def make_bump_text(imol, water_info):
        # print("debug:: in make_bump_text with imol", imol, "water_info", water_info)
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

    def make_store(number_list):
        name_store = Gtk.ListStore(int, str)
        for i in number_list:
            label_str = str(i)
            name_store.append([i, label_str])
        return name_store

    window = Gtk.Window()
    window.set_title("Coot Highly-Coordinated Waters")
    vbox = Gtk.VBox(False, 0)
    results_vbox = Gtk.VBox(False, 0)
    water_results_label = Gtk.Label("Other Coordinated Waters")
    metal_results_frame = Gtk.Frame()
    metal_results_vbox = Gtk.VBox(False, 0)
    metal_results_label = Gtk.Label("Potential Metals: ")
    h_sep = Gtk.HSeparator()
    hint_text = "Molecule: "
    hbox_chooser = Gtk.HBox(False, 0)
    hbox_max_dist = Gtk.HBox(False, 0)
    hbox_number_chooser = Gtk.HBox(False, 0)
    number_text = Gtk.Label("Coordination Number: ")
    # molecule_chooser_option_menu_and_model_list = generic_molecule_chooser(hbox_chooser, hint_text)
    # molecule_chooser_option_menu = molecule_chooser_option_menu_and_model_list[0]

    combobox_molecule = Gtk.ComboBox()
    combobox_mol_items = make_store_for_model_molecule_combobox()
    combobox_molecule.set_model(combobox_mol_items)

    renderer_text = Gtk.CellRendererText()
    if len(combobox_mol_items) > 0:
        combobox_molecule.set_active(0)
    combobox_molecule.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_molecule.pack_start(renderer_text, True)
    combobox_molecule.add_attribute(renderer_text, "text", 1)

    print("debug:: water_coordination_gui(): combobox_mol_items:", combobox_mol_items)
    print("debug:: water_coordination_gui(): combobox_molecule:",  combobox_molecule)
    # print("debug:: water_coordination_gui(): something:", something)

    combobox_molecule.set_active(0)

    # model_list = molecule_chooser_option_menu_and_model_list[1]
    scrolled_win = Gtk.ScrolledWindow()
    metal_results_scrolled_win = Gtk.ScrolledWindow()

    # https://python-gtk-3-tutorial.readthedocs.io/en/latest/combobox.html

    # coordination number combobox
    combobox_coordination = Gtk.ComboBoxText.new()
    for i in range(3,10):
        combobox_coordination.append_text(str(i))

    combobox_coordination.set_active(2)

    dist_label = Gtk.Label("Max Dist: ")
    dist_entry = Gtk.Entry()
    close_button = Gtk.Button("  Close  ")
    apply_button = Gtk.Button("  Apply  ")
    hbox_buttons = Gtk.HBox(False, 6)

    def get_molecule():
        tree_iter = combobox_molecule.get_active_iter()
        imol = -1
        if tree_iter is not None:
            model = combobox_molecule.get_model()
            it = model[tree_iter]
            imol = it[0]
        return imol

    def get_number():
        at = combobox_coordination.get_active_text()
        n = int(at)
        return n

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
            # for child in children:
            #   print "BL DEBUG:: now destroy child", child
            #   child.destroy()
            # this_vbox.hide()
            # this_vbox.show()

    def update_water_results(imol, n, d):
        results = coot.highly_coordinated_waters_py(imol, n, d)
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
                button = Gtk.Button(t + bump_text)
                if not is_a_metal_site_too_qm(atom_spec, metal_results):
                    results_vbox.pack_start(button, False, False, 1)
                    button.show()

                    def water_func(widget, imol, water_info):
                        water_spec = water_info[0]
                        chain_id = water_spec[1]
                        res_no = water_spec[2]
                        atom_name = water_spec[4]
                        coot.set_go_to_atom_molecule(imol)
                        coot.set_go_to_atom_chain_residue_atom_name(chain_id, res_no, atom_name)
                    button.connect("clicked", water_func, imol, water_info)

            # now handle metal results
            for metal_site in metal_results:
                metal_text = metal_site[1]
                t = atom_spec_to_text(metal_site[0])
                button_text = t + " Potential " + metal_text
                button = Gtk.Button(button_text)
                metal_results_vbox.pack_start(button, False, False, 1)
                button.show()

                def metal_func(widget, imol, metal_site):
                    metal_spec = metal_site[0]
                    chain_id = metal_spec[1]
                    res_no = metal_spec[2]
                    atom_name = metal_spec[4]
                    coot.set_go_to_atom_molecule(imol)
                    coot.set_go_to_atom_chain_residue_atom_name(
                        chain_id, res_no, atom_name)
                button.connect("clicked", metal_func, imol, metal_site)

    window.add(vbox)

    # fill_option_menu_with_number_options(number_menu, coot_utils.number_list, 5)

    scrolled_win.add_with_viewport(results_vbox)

    metal_results_scrolled_win.add_with_viewport(metal_results_frame)
    metal_results_frame.add(metal_results_vbox)
    # metal_results_scrolled_win.add_with_viewport(metal_results_vbox)
    

    vbox.pack_start(combobox_molecule, False, False, 2)

    hbox_max_dist.pack_start(dist_label, False, False, 2)
    hbox_max_dist.pack_start(dist_entry, False, False, 2)

    vbox.pack_start(hbox_chooser, False, False, 6)

    hbox_number_chooser.pack_start(number_text, False, False, 2)
    # hbox_number_chooser.pack_start(number_menu, False, False, 2)
    hbox_number_chooser.pack_start(combobox_coordination, False, False, 2)

    vbox.pack_start(hbox_number_chooser, False, False, 6)

    vbox.pack_start(hbox_max_dist, False, False, 2)

    # metal sites
    vbox.pack_start(metal_results_label, False, False, 2)
    vbox.pack_start(metal_results_scrolled_win, True, True, 0)

    # interesting water sites
    vbox.pack_start(water_results_label, False, False, 2)
    vbox.pack_start(scrolled_win, True, True, 0)  # expand fill padding
    vbox.pack_start(h_sep, False, False, 2)
    hbox_buttons.pack_start(apply_button, False, False, 2)
    hbox_buttons.pack_start(close_button, False, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 2)

    # From the Nayal and Di Cera (1996) paper, it seems that 2.7
    # and at least 4 oxygens is a good test for Na+ or other
    # interesting things
    #
    dist_entry.set_text("2.7")
    dist_entry.set_width_chars(5)

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
        res_no = spec[3]

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


# Note to self - using a python function to call a swigged function that calls a python function
# is completely painful (exceptions get lost (although not so much now that I've tried to fix that)).
# Never do this again - it was miserable trying to fix this.
#
def click_protein_db_loop_gui():

    print("DEBUG:: coot_utils:: click_protein_db_loop_gui() --- start ---")

    global db_loop_preserve_residue_names

    def pick_loop_func(n_clicks):

        print("DEBUG:: coot_utils:: click_protein_db_loop_gui() pick_loop_func() --- start ---")

        def atom_pick_func(*atom_specs):
            # residue_specs = list(map(coot_utils.atom_spec_to_residue_spec, atom_specs))
            for spec in atom_specs:
                print("debug: atom spec:", spec)
            residue_specs = [coot_utils.atom_spec_to_residue_spec(spec) for spec in atom_specs]
            imol = atom_specs[0][1]
            min_max_and_chain_id = min_max_residues_from_atom_specs(atom_specs)
            if not isinstance(min_max_and_chain_id, list):
                coot.info_dialog("Picked atoms not in same molecule and chain")
            else:
                loop_mols = coot.protein_db_loops_py(imol, residue_specs, coot.imol_refinement_map(),
                                                     10, db_loop_preserve_residue_names)
                imol_loop_orig = loop_mols[0][0]
                imol_loops_consolidated = loop_mols[0][1]
                loop_mols = loop_mols[1]
                min_resno = min_max_and_chain_id[0]
                max_resno = min_max_and_chain_id[1]
                ch_id = min_max_and_chain_id[2]

                def my_copy_function(imol, ch_id, loop_mol, min_resno, max_resno):
                    coot.copy_residue_range(imol, ch_id, loop_mol, ch_id, min_resno, max_resno)

                coot.set_mol_active(imol_loops_consolidated, 1)
                if coot_utils.valid_model_molecule_qm(imol_loop_orig):
                    if len(loop_mols) > 0:
                        # ninja python skills (thank you stackoverflow)
                        buttons = [["Insert " + str(j_mol) + " " + coot.molecule_name(j_mol), lambda button, j_mol=j_mol: my_copy_function(imol, ch_id, j_mol, min_resno, max_resno)]
                                   for j_mol in loop_mols]

                        loop_mols.append(imol_loops_consolidated)

                        def toggle_func(imol):
                            coot_utils.toggle_active_mol(imol)
                            coot_utils.toggle_display_mol(imol)

                        all_buttons = [["Original loop", lambda func: coot.copy_residue_range(imol, chain_id,
                                                                                              imol_loop_orig, chain_id,
                                                                                              min_resno, max_resno)],
                                       ["Toggle Display All Candidate Loops", lambda func: toggle_func(imol_loops_consolidated)]
                                      ] + buttons

                        dialog_box_of_buttons("Loop Candidates", [360, 200], all_buttons, " Close ",
                                               lambda: [(coot.set_mol_displayed(im, 0), coot.set_mol_active(im, 0)) for im in loop_mols])


        coot.user_defined_click_py(n_clicks, atom_pick_func)

    print("DEBUG:: coot_utils:: click_protein_db_loop_gui() --- calling generic_number_chooser() ---")
    # this is wrong! The 4th element has a label "6" (all are offset by 2)
    generic_number_chooser(list(range(2, 10)), 4,
                           "Number of residues for basis",
                           "Pick Atoms...",
                           lambda n_clicks: pick_loop_func(n_clicks))


def refmac_multi_sharpen_gui():

    def delete_event(*args):
        window.destroy()
        return False

    def sharpen_cb(widget, *args):

        # get max_band n_levels and map file name
        max_b = int(get_option_menu_active_item(option_menu_b_factor, b_factor_list))
        n_levels = int(get_option_menu_active_item(option_menu_n_levels, n_levels_list))
        active_item_imol = get_option_menu_active_molecule(option_menu_map, coot_utils.map_molecule_list)


        # There is no function to get a map file name from a molecule
        # It is not stored. So we make/guess it...
        map_file_name = coot.molecule_name(active_item_imol)
        if (map_file_name.find(" ") > 0):
            # we have map coeffs - but then sharpen as here wont work anyway
            map_file_name = map_file_name[:map_file_name.find(" ")]
        map_file_name_stub = coot_utils.strip_path(
            coot_utils.file_name_sans_extension(map_file_name))
        refmac_output_mtz_file_name = "starting_map-" + map_file_name_stub + ".mtz"
        log_file_name = "refmac-multisharp-" + map_file_name_stub + ".log"
        if not os.path.isfile(map_file_name):
            coot.info_dialog("WARNING:: file not found %s" % map_file_name)
        else:
            if not coot_utils.directory_is_modifiable_qm(os.getcwd()):
                m = "ERROR:: Current directory " + os.getcwd() + " is not writable"
                coot.info_dialog(m)
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
            if not coot_utils.directory_is_modifiable_qm(this_dir):
                coot.info_dialog("WARNING:: Current directory is not writable")
            else:
                refmac_execfile = coot_utils.find_exe(
                    "refmac5", "CBIN", "CCP4_BIN", "PATH")
                s = coot_utils.popen_command(refmac_execfile,
                                             cmd_line_args,
                                             data_lines,
                                             log_file_name,
                                             False)

                try:
                    if s != 0:
                        coot.info_dialog("WARNING:: refmac5 failed")
                    else:
                        # Happy path
                        print("BL DEBUG:: s", s)
                        if os.path.isfile("starting_map.mtz"):
                            os.rename("starting_map.mtz",
                                      refmac_output_mtz_file_name)
                            # offer a read-mtz dialog
                            coot.manage_column_selector(refmac_output_mtz_file_name)

                except:
                    print("BL DEBUG:: tried to rename starting-map.mtz but failed.")
                    pass
            delete_event(widget)

    window = Gtk.Window()
    # boxes
    vbox = Gtk.VBox(False, 0)
    hbox_1 = Gtk.HBox(False, 0)
    hbox_2 = Gtk.HBox(False, 0)
    hbox_3 = Gtk.HBox(False, 0)
    # menus
    # option_menu_map = Gtk.combo_box_new_text()
    # option_menu_b_factor = Gtk.combo_box_new_text()
    # option_menu_n_levels = Gtk.combo_box_new_text()


    # labels
    map_label = Gtk.Label("Map ")
    sb_label = Gtk.Label("Sharpen & Blur in ")
    levels_label = Gtk.Label(" levels up to ")
    A_label = Gtk.Label(" A*A")
    # separate
    h_sep = Gtk.HSeparator()
    # buttons
    ok_button = Gtk.Button("   OK   ")
    cancel_button = Gtk.Button(" Cancel ")
    n_levels_list = [1, 2, 3, 4, 5, 6]
    b_factor_list = [50, 100, 200, 400, 800, 2000]

    combobox_map_items = make_store_for_map_molecule_combobox()
    combobox_map = Gtk.ComboBox.new_with_model(combobox_map_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_map_items) > 0:
        combobox_map.set_active(0)
    combobox_map.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_map.pack_start(renderer_text, True)
    combobox_map.add_attribute(renderer_text, "text", 1)

    combobox_n_levels = Gtk.ComboBox.new_with_model(combobox_n_levels_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_map_items) > 0:
        combobox_n_levels.set_active(0)
    combobox_n_levels.set_entry_text_column(0)
    combobox_n_levels.pack_start(renderer_text, True)
    combobox_n_levels.add_attribute(renderer_text, "text", 0)

    combobox_b_factor = Gtk.ComboBox.new_with_model(combobox_b_factor_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_map_items) > 0:
        combobox_b_factor.set_active(0)
    combobox_b_factor.set_entry_text_column(0)
    combobox_b_factor.pack_start(renderer_text, True)
    combobox_b_factor.add_attribute(renderer_text, "text", 0)

    # fill_option_menu_with_number_options(option_menu_n_levels, n_levels_list, 4)
    # fill_option_menu_with_number_options(option_menu_b_factor, b_factor_list, 200)

    window.set_title("Refmac for Sharpening & Blurring")
    hbox_1.pack_start(map_label,         False, False, 2)
    hbox_1.pack_start(combobox_map,      False, False, 2)
    hbox_2.pack_start(sb_label,          False, False, 2)
    hbox_2.pack_start(combobox_n_levels, False, False, 2)
    hbox_2.pack_start(levels_label,      False, False, 2)
    hbox_3.pack_end(cancel_button,       False, False, 12)
    hbox_3.pack_end(ok_button,           False, False, 12)

    vbox.pack_start(hbox_1)
    vbox.pack_start(hbox_2)
    vbox.pack_start(h_sep)
    vbox.pack_start(hbox_3)

    cancel_button.connect("clicked", delete_event)

    ok_button.connect("clicked", sharpen_cb, combobox_b_factor, b_factor_list,
                      combobox_n_levels, n_levels_list,
                      combobox_map, coot_utils.map_molecule_list)

    window.add(vbox)
    window.show_all()


def auto_assign_sequence_from_map():

    active_atom = active_residue()
    # get these from the current fragment
    imol  = active_atom[0]
    ch_id = active_atom[1]
    res_no = active_atom[2]
    res_spec = [ch_id, res_no, ""]
    imol_map = imol_refinement_map()
    fragment_residues = linked_residues_py(res_spec, imol, 1.7)
    residue_number_list = [spec[2] for spec in fragment_residues] # bleugh. Should be spec[1]
    resno_start = min(residue_number_list)
    resno_end   = max(residue_number_list)
    new_sequence = sequence_from_map(imol, ch_id, resno_start, resno_end, imol_map)
    set_rotamer_search_mode(ROTAMERSEARCHLOWRES)
    coot.mutate_residue_range(imol, ch_id, resno_start, resno_end, new_sequence)
    backrub_rotamers_for_chain(imol, ch_id)
    refine_residues(imol, fragment_residues)

    
# ;; Associate the contents of a PIR file with a molecule.  Select file from a GUI.
# ;; 
# (define (associate-pir-with-molecule-gui do-alignment?)

#   (format #t "in associate-pir-with-molecule-gui~%") 
#   (generic-chooser-entry-and-file-selector 
#    "Associate PIR Sequence to Model: "
#    valid-model-molecule?
#    "Chain ID"
#    ""
#    "Select PIR Alignment file"
#    (lambda (imol chain-id file-name)
#      (associate-pir-file imol chain-id file-name)
#      (if do-alignment?
#         (alignment-mismatches-gui imol)))
#     *generic-chooser-entry-and-file-selector-file-entry-default-text*))


def associate_pir_wih_molecule_gui(do_alignment_flag):

    def lambda_fn(imol, ch_id, file_name):
        coot_utils.associate_pir_file(imol, ch_id, file_name)
        if do_alignment_flag:
            alignment_mismatches_gui(imol)

    generic_chooser_entry_and_file_selector("Associate PIR Sequence to Model",
                                            coot_utils.valid_model_molecule_qm, "Chain ID", "", " PIR Alignment File: ",
                                            lambda imol, ch_id, file_name: lambda_fn(imol, ch_id, file_name))


def add_module_cryo_em():
    if coot_gui_api.main_menubar():
        add_module_cryo_em_gui()


def add_module_ccp4():
    if coot_gui_api.main_menubar():
        add_module_ccp4_gui()


def add_module_pdbe():
   if coot_python.main_menubar():
      add_module_pdbe_gui()

def add_module_cryo_em_gui():

    def solidify_maps(w):
        for imol in range(coot.graphics_n_molecules()):
            if coot.is_valid_map_molecule(imol):
                coot.set_draw_solid_density_surface(imol, 1)

    def unsolidify_maps(w):
        for imol in range(coot.graphics_n_molecules()):
            if coot.is_valid_map_molecule(imol):
                coot.set_draw_solid_density_surface(imol, 0)

    def add_mol_sym_mtrix():
         with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                        aa_ins_code, aa_atom_name,
                                        aa_alt_conf, aa_res_spec]:
            add_molecular_symmetry_from_mtrix_from_self_file(aa_imol)
            set_show_symmetry_master(1)

    def go_to_box_middle():
        m_list = map_molecule_list()
        if len(m_list) > 0:
            m = m_list[-1]
            c = cell(m)
            set_rotation_centre(0.5 * c[0], 0.5 * c[1], 0.5 * c[2])

        add_simple_coot_menu_menuitem(menu, "Sharpen/Blur...",
                                      lambda func: sharpen_blur.sharpen_blur_map_gui())

    def flip_hand_local_func():
        map_molecule_chooser_gui("Select", lambda imol: flip_hand(imol))

        add_simple_coot_menu_menuitem(menu, "Add molecular symmetry using MTRIX",
                                    lambda func: add_mol_sym_mtrix())

        add_simple_coot_menu_menuitem(menu, "Sharpen/Blur...",
                                    lambda func: sharpen_blur_map_gui())

    def make_masked_maps_using_active_atom():
        active_atom = coot.active_residue_py()
        print("active_atom:", active_atom)
        if active_atom:
            imol = active_atom[0]
            coot.make_masked_maps_split_by_chain(imol, coot.imol_refinement_map())
        
    def go_to_box_middle():
        m_list = coot_utils.map_molecule_list()
        if len(m_list) > 0:
            m = m_list[-1]
            c = rf.cell(m)
            coot.set_rotation_centre(0.5 * c[0], 0.5 * c[1], 0.5 * c[2])

    if coot_gui_api.main_menubar():

        def ass_seq_assoc_seq():
            assign_sequence_to_active_fragment()

        def interactive_nudge_func():
            with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                                      aa_ins_code, aa_atom_name,
                                                      aa_alt_conf, aa_res_spec]:
                interactive_nudge_residues.nudge_residues_gui(aa_imol, aa_res_spec)
                
        menu = coot_menubar_menu("Cryo-EM")

        add_simple_coot_menu_menuitem(menu, "Multi-sharpen...",
                                      lambda func: refmac_multi_sharpen_gui())

        add_simple_coot_menu_menuitem(menu, "Sharpen/Blur...",
                                      lambda func: sharpen_blur.sharpen_blur_map_gui())

        add_simple_coot_menu_menuitem(menu, "Mask Map by Chains",
                                      lambda func: make_masked_maps_using_active_atom())

        add_simple_coot_menu_menuitem(menu, "Interactive Nudge Residues...",
                                      lambda func: interactive_nudge_func())

        add_simple_coot_menu_menuitem(menu, "Go To Map Molecule Middle",
                                      lambda func: coot.go_to_map_molecule_centre(coot.imol_refinement_map()))

        add_simple_coot_menu_menuitem(menu, "Go To Box Middle",
                                      lambda func: go_to_box_middle())

        add_simple_coot_menu_menuitem(menu, "Flip Hand of Map",
                                      lambda func: flip_hand_local_func())

        add_simple_coot_menu_menuitem(menu, "Add molecular symmetry using MTRIX",
                                    lambda func: add_mol_sym_mtrix())

        add_simple_coot_menu_menuitem(menu, "Align and Mutate using ClustalW2",
                                    lambda func:
                                    generic_chooser_entry_and_file_selector(
                                        "Align Sequence to Model: ",
                                        coot_utils.valid_model_molecule_qm,
                                        "Chain ID",
                                        "",
                                        "Select PIR Alignment file",
                                        lambda imol, chain_id, target_sequence_pif_file:
                                        run_clustalw_alignment(imol, chain_id,
                                                               target_sequence_pif_file)))

        add_simple_coot_menu_menuitem(menu, "Assign Sequence Based on Associated Sequence",
                                      lambda func: ass_seq_assoc_seq())

        add_simple_coot_menu_menuitem(menu, "Auto-assign Sequence Based on Map",
                                      lambda func: auto_assign_sequence_from_map())

        add_simple_coot_menu_menuitem(menu, "No Auto-Recontour Map Mode",
                                      lambda func: set_auto_recontour_map(0))

        add_simple_coot_menu_menuitem(menu, "Enable Auto-Recontour Map Mode",
                                      lambda func: set_auto_recontour_map(1))

        add_simple_coot_menu_menuitem(menu, "Interactive Nudge Residues...",
                                      lambda func: interactive_nudge_func())


def add_module_ccp4_gui():
    if coot_gui_api.main_menubar():
        menu = coot_menubar_menu("CCP4")

        add_simple_coot_menu_menuitem(menu, "Make LINK via Acedrg",
                                      lambda func: acedrg_link.acedrg_link_generation_control_window())


def add_module_refine():

   def chain_refine_active_atom(widget):
      active_atom = active_residue()
      if active_atom:
         aa_imol     = active_atom[0]
         aa_chain_id = active_atom[1]
         all_residues = residues_in_chain(aa_imol, aa_chain_id)
         refine_residues(aa_imol, all_residues);

   def all_atom_refine_active_atom(widget):
      active_atom = active_residue()
      if active_atom:
         aa_imol = active_atom[0]
         all_residues_in_mol = all_residues(aa_imol)
         refine_residues(aa_imol, all_residues_in_mol);

   def refine_fragment_active_atom(w):
      active_atom = active_residue()
      print("###### active_atom", active_atom)
      if active_atom:
         aa_imol = active_atom[0]
         aa_res_spec = [active_atom[1], active_atom[2], active_atom[3]] # doesn't ative_residue 
         res_list = linked_residues_py(aa_res_spec, aa_imol, 1.7)
         refine_residues(aa_imol, res_list);

   def regularize_fragment_active_atom(w):
      active_atom = active_residue()
      if active_atom:
         aa_imol = active_atom[0]
         aa_res_spec = [active_atom[1], active_atom[2], active_atom[3]] # doesn't ative_residue 
         res_list = linked_residues_py(aa_res_spec, aa_imol, 1.7)
         regularize_residues(aa_imol, res_list);

   def regularize_chain_active_atom(w):
      active_atom = active_residue()
      if active_atom:
         aa_imol = active_atom[0]
         aa_chain_id = active_atom[1]
         all_residues = residues_in_chain(aa_imol, aa_chain_id)
         regularize_residues(aa_imol, all_residues);
      
   if coot_python.main_menubar():
      menu = coot_menubar_menu("Refine")

      add_simple_coot_menu_menuitem(menu, "All-Atom Refine", all_atom_refine_active_atom)

      add_simple_coot_menu_menuitem(menu, "Chain Refine", chain_refine_active_atom)

      # they get turned on but are not active - they currently need to be turn off by the user using the Generic Display dialog
      add_simple_coot_menu_menuitem(menu, "Contact Dots On",  lambda widget: set_do_coot_probe_dots_during_refine(1))
      add_simple_coot_menu_menuitem(menu, "Contact Dots Off", lambda widget: set_do_coot_probe_dots_during_refine(0))

      add_simple_coot_menu_menuitem(menu, "Intermediate Atom Restraints On",  lambda widget: set_draw_moving_atoms_restraints(1))
      add_simple_coot_menu_menuitem(menu, "Intermediate Atom Restraints Off", lambda widget: set_draw_moving_atoms_restraints(0))

      add_simple_coot_menu_menuitem(menu, "Refine Fragment", refine_fragment_active_atom)

      add_simple_coot_menu_menuitem(menu, "Regularize Fragment", regularize_fragment_active_atom)

      add_simple_coot_menu_menuitem(menu, "Regularize Chain", regularize_chain_active_atom)

      add_simple_coot_menu_menuitem(menu, "Rama Goodness Dodecs On",  lambda w: set_show_intermediate_atoms_rota_markup(1))
      add_simple_coot_menu_menuitem(menu, "Rama Goodness Dodecs Off", lambda w: set_show_intermediate_atoms_rota_markup(0))


def scale_alt_conf_occ_gui(imol, chain_id, res_no, ins_code):

    alt_confs = coot_utils.residue_alt_confs(imol, chain_id, res_no, ins_code)
    try:
        # remove no alt conf, i.e. '', from list if present
        alt_confs.remove('')
    except:
        pass

    if (len(alt_confs) != 2):
        print("INFO:: no (or too many) alt confs, no gui!")
    else:
        # only do if 2 alt confs (plus no)
        res_info = coot.residue_info_py(imol, chain_id, res_no, ins_code)
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
            alt_conf_list = [[alt_confs[0], new_occ],
                             [alt_confs[1], 1 - new_occ]]
            coot_utils.set_alt_conf_occ(
                imol, chain_id, res_no, ins_code, alt_conf_list)
            delete_event()

        def change_occ(*args):
            # needed?!
            # print "change occ slider?!"
            pass

        window = Gtk.Window()
        title = Gtk.Label("Adjust alt conf occupancies")
        occ_label = Gtk.Label("Occupancy")
        alt_conf_label = Gtk.Label("Alt Conf: " + alt_conf)
        occ_adj = Gtk.Adjustment(occ_start, 0.1, 0.99, 0.01, 0.1, 0.1)
        occ_scale = Gtk.HScale(occ_adj)
        vbox = Gtk.VBox(False, 0)
        scale_hbox = Gtk.HBox(False, 0)
        scale_hbox.pack_start(alt_conf_label, False, False, 2)
        scale_hbox.pack_start(occ_scale, True, True, 2)
        vbox.pack_start(occ_label, False, False, 0)
        vbox.pack_start(scale_hbox, False, False, 0)

        occ_adj.connect("value_changed", change_occ)

        window.add(vbox)

        h_sep = Gtk.HSeparator()
        buttons_hbox = Gtk.HBox(True, 6)
        ok_button = Gtk.Button("   OK   ")
        cancel_button = Gtk.Button(" Cancel ")

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
        imol = args[0][1]
        chain_id = args[0][2]
        res_no = args[0][3]
        ins_code = args[0][4]
        scale_alt_conf_occ_gui(imol, chain_id, res_no, ins_code)
    coot.user_defined_click_py(1, helper_function)


def toggle_backrub_rotamers(widget=None):
    """Toggle function to swtich on and off backrub rotamer fitting.

    Keyword arguments:
    widget -- can be passed from the toolbutton

    """

    if widget:
        if widget.get_active():
            # the button is toggled on
            coot.set_rotamer_search_mode(ROTAMERSEARCHLOWRES)
            print("INFO:: Using Backrub rotamers now")
        else:
            coot.set_rotamer_search_mode(ROTAMERSEARCHHIGHRES)
            print("INFO:: No longer using Backrub rotamers")

    else:
        # non graphical - but wont be able to run if this is not loaded.
        mode = coot.rotamer_search_mode_state()
        if (mode == ROTAMERSEARCHLOWRES):
            coot.set_rotamer_search_mode(ROTAMERSEARCHHIGHRES)
            print("INFO:: No longer using Backrub rotamers")
        if (mode == ROTAMERSEARCHHIGHRES or
                mode == ROTAMERSEARCHAUTOMATIC):
            coot.set_rotamer_search_mode(ROTAMERSEARCHLOWRES)
            print("INFO:: Using Backrub rotamers")

        # no alternative for now
        # need to be able to get the state of search mode.
        # easily added. FIXME
        print("BL WARNING:: no widget")

def atom_overlaps_for_this_model():
   """Display Atom overlaps for active atom"""
   active_atom = coot.active_residue_py()
   if active_atom:
      aa_imol = active_atom[0]
      coot.coot_all_atom_contact_dots(aa_imol)


def toggle_hydrogen_display(widget=None):
    """Toggle function to display all hydrogens or not.

    Keyword arguments:
    widget -- can be passed from the toolbutton

    """

    if widget:
        if widget.get_active():
            # the button is toggled on
            coot_utils.hide_all_hydrogens()
        else:
            coot_utils.show_all_hydrogens()

    else:
        # non graphical - but wont be able to run if this is not loaded.
        print("BL INFO:: No display, so I dont care about the hydrogens.")
        print("BL WARNING:: no widget")


# Simple minded dialog, search disk or not?
# Could/should be expaded to browse for exe and to select which
# disks to search for.
#


def search_disk_dialog(program_name, path_ls):

    # graphics
    ret = False
    label_text = "Couldn't find %s in default path" % (program_name)
    for path in path_ls:
        label_text += " and "
        label_text += path
    label_text += "\n\nShall we search the whole disk?\n"

    try:
        ret = yes_no_dialog(label_text, "Search whole disk dialog")
    except:
        # no graphics
        label_text += "[y/N] >"
        result = ""
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

        residue_specs = list(
            map(coot_utils.atom_spec_to_residue_spec, coot_utils.atom_specs))
        imol_1 = coot_utils.atom_specs[0][1]
        imol_2 = coot_utils.atom_specs[1][1]
        chain_id1 = coot_utils.atom_specs[0][2]
        chain_id2 = coot_utils.atom_specs[1][2]
        res_no_1 = coot_utils.atom_specs[0][3]
        res_no_2 = coot_utils.atom_specs[1][3]

        # some sanity check
        if not (imol_1 == imol_2):
            msg = (
                "WARNING:: not the same imols. \n"
                "imol %i and %i were selected" % (imol_1, imol_2))
            coot.info_dialog_and_text(msg)
            return
        else:
            # imol ok
            if (not chain_id1 == chain_id2):
                msg = (
                    "BL WARNING:: not the same chains. \n"
                    "Chains %s and %s were selected" % (chain_id1, chain_id2))
                coot.info_dialog_and_text(msg)
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
                        % (abs(res_diff)))
                    coot.info_dialog_and_text(msg)
                    return
                elif (abs(res_diff) < 0):
                    msg = (
                        "BL WARNING::  No residue selected.")
                    coot.info_dialog_and_text(msg)
                    return
                if (res_no_1 > res_no_2):
                    # need to swap
                    res_no_1, res_no_2 = res_no_2, res_no_1

                # main line
                # NB: no occupancy setting here
                coot_utils.duplicate_residue_range(
                    imol_1, chain_id1, res_no_1, res_no_2)

    coot.add_status_bar_text("Pick two atoms")
    coot.user_defined_click_py(2, pick_range_func)

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

    dialog = Gtk.Dialog(title_text, None,
                        Gtk.DIALOG_MODAL | Gtk.DIALOG_NO_SEPARATOR,
                        (Gtk.STOCK_YES, Gtk.RESPONSE_ACCEPT,
                         Gtk.STOCK_NO, Gtk.RESPONSE_REJECT))
    ifont = Gtk.gdk.Font("fixed")
    label = Gtk.Label(label_text)
    dialog.vbox.pack_end(label, True, True, 0)
    dialog.show_all()
    result = dialog.run()
    if result == Gtk.RESPONSE_ACCEPT:
        ret = True
    else:
        ret = False
    dialog.destroy()

    return ret

# Function based on Davis et al. (2007) Molprobity: all atom contacts
# and structure validation for proteins and nucleic acids, Nucleic
# Acids Research 35, W375-W383.
#
#    "RNA sugar puckers (C3'endo or C2'endo) is strongly correlated
#    to the perpendicular distance between the following (3')
#    phosphate and either the plane of the base or the C1'-N1/9
#    glycosidic bond vector. [] .. a sugar pucker is very difficult
#    to determine directly from the electron density at resolutions
#    typical for RNAs."
#
# To paraphrase:
# The distance of the plane of the base to the following phosphate
# is highly correlated to the pucker of the ribose.
#
# An analysis of the structures in RNADB2005 shows that a critical
# distance of 3.3A provides a partition function to separate C2' from
# C3' endo puckering.  Not all ribose follow this rule.  There may be
# some errors in the models comprising RNADB2005. So we check the
# distance of the following phosphate to the plane of the ribose and
# record the riboses that are inconsitent.  We also report puckers
# that are not C2' or C3'.  The puckers are determined by the most
# out-of-plane atom of the ribose (the rms deviation of the 4 atoms
# in the plane is calculated, but not used to determine the
# puckering atom).
#
def pukka_puckers_qm(imol):

    import types

    residue_list = []
    crit_d = 3.0  # Richardson's group value to partition C2'-endo from C3'-endo

    def add_questionable(r):
        residue_list.append(r)

    def get_ribose_residue_atom_name(imol, residue_spec, pucker_atom):
        r_info = coot.residue_info_py(imol, residue_spec[0], residue_spec[1], residue_spec[2])
        t_pucker_atom = pucker_atom[0:3] + "*"
        if (pucker_atom in [at[0][0] for at in r_info]):
            return pucker_atom
        else:
            return t_pucker_atom

    # main line
    for chain_id in coot_utils.chain_ids(imol):
        if (not coot_utils.is_solvent_chain_qm(imol, chain_id)):
            n_residues = coot.chain_n_residues(chain_id, imol)

            for serial_number in range(n_residues):

                res_name = coot.resname_from_serial_number(imol, chain_id, serial_number)
                res_no = coot.seqnum_from_serial_number(imol, chain_id, serial_number)
                ins_code = coot.insertion_code_from_serial_number(imol, chain_id, serial_number)

                if (not res_name == "HOH"):

                    residue_spec = [chain_id, res_no, ins_code]
                    pi = coot.pucker_info_py(imol, residue_spec, 1)
                    if pi:
                        try:
                            if len(pi) == 4:
                                pucker_atom = pi[1]
                                if ((abs(pi[0]) > crit_d) and (pucker_atom == " C2'")):
                                    add_questionable([pucker_atom, residue_spec, "Inconsistent phosphate distance for C2' pucker"])
                                if ((abs(pi[0]) < crit_d) and (pucker_atom == " C3'")):
                                    add_questionable([pucker_atom, residue_spec, "Inconsistent phosphate distance for C3' pucker"])
                                if not ((pucker_atom == " C2'") or (pucker_atom == " C3'")):
                                    add_questionable([pucker_atom, residue_spec, "puckered atom:" + pucker_atom])
                        except NameError as e:
                            print(e)
                        except TypeError as e:
                            print(e)

    def go(imol, residue_spec, at_name):
        ch_id = residue_spec_to_chain_id(residue_spec)
        res_no = residue_spec_to_res_no(residue_spec)
        coot.set_go_to_atom_molecule(imol)
        coot.set_go_to_atom_chain_residue_atom_name(ch_id, res_no, at_name)

    def generator(imol, residue_spec, at_name):
        func = lambda imol_c = imol, residue_spec_c = residue_spec, at_name_c = at_name : go(imol_c, residue_spec_c, at_name_c)
        def action(arg):
            func(imol, residue_spec, at_name)
        return action

    if len(residue_list) == 0:
        coot.info_dialog("No bad puckers.")
    else:
        buttons = []
        for residue in residue_list:
            residue_spec = residue[1]
            pucker_atom  = residue[0]
            info_string  = residue[2]
            at_name = get_ribose_residue_atom_name(imol, residue_spec, pucker_atom)
            label = residue_spec[0] + " " + str(residue_spec[1]) + residue_spec[2] + ": " + info_string
            func = generator(imol, residue_spec, at_name)
            ls = [label, func]
            buttons.append(ls)
        dialog_box_of_buttons("Non-pukka puckers", [370, 250], buttons, "  Close  ")

# A gui to list Ramachandran outliers etc.
# May become more sophisticated at some point
#
def rama_outlier_gui():

    """
    A gui to list Ramachandran outliers etc.
    A first draft, may become more sophisticated at some point
    """

    def list_rama_outliers(imol):

        r = coot.all_molecule_ramachandran_region_py(imol)
        outliers = []
        allowed = []
        for res in r:
            if res[1] == 0:
                outliers.append(res[0])
            if res[1] == 1:
                allowed.append(res[0])

    def make_buttons_old(res_list, label_string):

         ret = []
         for res_spec in res_list:
            chain_id = res_spec[1]
            res_no = res_spec[2]
            label = label_string + ": " + \
                    chain_id + " " + str(res_no)
            func = [cmd2str(set_go_to_atom_molecule, imol),
                    cmd2str(set_go_to_atom_from_res_spec, res_spec)]
            ret.append([label, func])

         return ret

    def go_to_residue(res_spec):
        print("Here in go_to_residue() with res_spec", res_spec)

    def generator(res_spec):
        func = lambda res_spec_c=res_spec : go_to_residue(res_spec_c)
        def action(arg):
            func(res_spec)
        return action

    def make_buttons(res_list, label_prefix):

        button_list = []
        for res_spec in res_list:
            chain_id = res_spec[1]
            res_no   = res_spec[2]
            ins_code = res_spec[3]
            label = label_prefix + ": " + chain_id + " " + str(res_no)
            if ins_code:
                label += " " + ins_code
            func = generator(res_spec)
            buttons.append([label, func])

        return button_list

    outlier_buttons = make_buttons(outliers, "Outlier")
    allowed_buttons = make_buttons(allowed, "Allowed")
    all_buttons = outlier_buttons + allowed_buttons

    def clear_and_add_back(vbox, outliers_list, allowed_list, filter_flag):
         # clear
         children = vbox.get_children()
         map(lambda c: c.destroy(), children)
         # add back
         if not filter_flag:
            buttons = outliers_list + allowed_list
         else:
            # filter
            buttons = outliers_list
         map(lambda button_info: add_button_info_to_box_of_buttons_vbox(button_info, vbox),
             buttons)

    dialog_box_of_buttons_with_check_button(
          " Ramachandran issues ", [300, 300], [], "  Close  ",
          "Outliers only",
          lambda check_button, vbox: clear_and_add_back(vbox, outlier_buttons, allowed_buttons, True)
          if check_button.get_active() else
          clear_and_add_back(vbox, outlier_buttons, allowed_buttons, False),
          False)

    molecule_chooser_gui("List Rama outliers for which molecule?", lambda imol: list_rama_outliers(imol))
   

def model_map_diff_map_molecule_chooser_gui(callback_function):

    def delete_event(*args):
       window.destroy()
       return False

    def combobox_to_molecule_number(combobox):
        imol = -1
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol = it[0]
        return imol

    def on_ok_clicked(*args):
        active_mol_no_model    = combobox_to_molecule_number(combobox_model)
        active_mol_no_map      = combobox_to_molecule_number(combobox_map)
        active_mol_no_diff_map = combobox_to_molecule_number(combobox_diff_map)
        try:
           active_mol_no_model    = int(active_mol_no_model)
           active_mol_no_map      = int(active_mol_no_map)
           active_mol_no_diff_map = int(active_mol_no_diff_map)
           print("INFO: operating on molecule numbers ", active_mol_no_model, active_mol_no_map, active_mol_no_diff_map)
           try:
               auto_button_state = auto_update_checkbutton.get_active()
               callback_function(active_mol_no_model, active_mol_no_map, active_mol_no_diff_map, auto_button_state)
           except TypeError as e:
               print(e)
               print("BL INFO:: problem in callback_function", callback_function.func_name)
           delete_event()
        except TypeError as e:
            print(e)
            print( "Failed to get a (molecule) number")

    window = Gtk.Window(title="Coot: Map Molecule Chooser")
    model_chooser_label = "Model"
    map_chooser_label   = "Map (with data info attached)"
    diff_map_chooser_label   = "Difference Map"
    label_for_map      = Gtk.Label(map_chooser_label)
    label_for_diff_map = Gtk.Label(diff_map_chooser_label)
    label_for_model    = Gtk.Label(model_chooser_label)
    vbox = Gtk.VBox(False,6)
    hbox_buttons = Gtk.HBox(False,5)

    # 20220326-PE
    # combobox_model    = Gtk.combo_box_new_text()
    # combobox_diff_map = Gtk.combo_box_new_text()
    # combobox_map      = Gtk.combo_box_new_text()

    combobox_diff_map = Gtk.ComboBox()

    combobox_mol_items = make_store_for_model_molecule_combobox()
    combobox_model = Gtk.ComboBox.new_with_model(combobox_mol_items)

    renderer_text_for_model = Gtk.CellRendererText()
    if len(combobox_mol_items) > 0:
        combobox_model.set_active(0)
    combobox_model.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_model.pack_start(renderer_text_for_model, True)
    combobox_model.add_attribute(renderer_text_for_model, "text", 1)

    combobox_map_items = make_store_for_map_molecule_combobox()
    combobox_map = Gtk.ComboBox.new_with_model(combobox_map_items)
    
    renderer_text_for_map = Gtk.CellRendererText()
    if len(combobox_map_items) > 0:
        combobox_map.set_active(0)
    combobox_map.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_map.pack_start(renderer_text_for_map, True)
    combobox_map.add_attribute(renderer_text_for_map, "text", 1)

    combobox_diff_map_items = make_store_for_diff_map_molecule_combobox()
    combobox_diff_map = Gtk.ComboBox.new_with_model(combobox_diff_map_items)

    renderer_text_for_diff_map = Gtk.CellRendererText()
    if len(combobox_diff_map_items) > 0:
        combobox_diff_map.set_active(0)
    combobox_diff_map.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_diff_map.pack_start(renderer_text_for_diff_map, True)
    combobox_diff_map.add_attribute(renderer_text_for_diff_map, "text", 1)

    ok_button = Gtk.Button("  OK  ")
    cancel_button = Gtk.Button(" Cancel ")
    h_sep = Gtk.HSeparator()

    # "auto" button
    auto_update_checkbutton = Gtk.CheckButton("Auto Update")

    window.set_default_size(370,100)
    window.add(vbox)
    vbox.pack_start(label_for_model,False,False, 5)
    vbox.pack_start(combobox_model, True, True,  2)
    vbox.pack_start(label_for_map,  False,False, 5)
    vbox.pack_start(combobox_map,   True, True,  2)
    vbox.pack_start(label_for_diff_map, False, False, 5)
    vbox.pack_start(combobox_diff_map,  True,  True,  2)
    vbox.pack_start(auto_update_checkbutton, False, False, 2)
    vbox.pack_start(h_sep,True,False, 2)
    vbox.pack_start(hbox_buttons,False, False,5)
    hbox_buttons.pack_start(cancel_button, True, False, 5)
    hbox_buttons.pack_start(    ok_button, True, False, 5)

    # button callbacks:
    ok_button.connect("clicked", on_ok_clicked)
    cancel_button.connect("clicked", delete_event)

    window.show_all()

def show_updating_maps_chooser():

    def update_maps_func(toolbutton):
        # the R-factor is added to the status bar
        stats_result = coot.calculate_maps_and_stats_py(imol, imol_map, imol_map, imol_diff_map)

    def generator_update_maps_func_wrap_for_capture(imol, imol_map, imol_diff_map):
        func = lambda imol_c=imol, imol_map_c=imol_map, imol_diff_map_c=imol_diff_map : coot.calculate_maps_and_stats_py(imol, imol_map, imol_map, imol_diff_map)
        def action(arg):
           func(imol, imol_map, imol_diff_map)
        return action

    # this needs  the same generator work-around but I've run out of energy.
    def shiftfield_func():
        imol = 0
        coot.shiftfield_b_factor_refinement(imol)

    def ok_button_function(imol, imol_map, imol_diff_map, use_auto_update_mode):
        # print("imol for coords", imol, "imol_for_map", imol_map, "use_auto_update_mode", use_auto_update_mode)
        if coot.is_valid_map_molecule(imol_map):
            if use_auto_update_mode:
                coot.set_auto_updating_sfcalc_genmap(imol, imol_map, imol_diff_map)
            else:
                print("================================= calculate_maps_and_stats_py()", imol, imol_map, imol_diff_map)
                coot.calculate_maps_and_stats_py(imol, imol_map, imol_map, imol_diff_map)
            menu_bar_callback_func = generator_update_maps_func_wrap_for_capture(imol, imol_map, imol_diff_map)
            coot_toolbar_button("Update Maps", menu_bar_callback_func, use_button=True)
            coot_toolbar_button("Shiftfield B", shiftfield_func)

    model_map_diff_map_molecule_chooser_gui(ok_button_function)



# let the c++ part of mapview know that this file was loaded:

# print "From coot_gui.py calling coot.set_found_coot_python_gui()" - called from rcrane loader
# too (I think)
coot.set_found_coot_python_gui()
