# Copyright 2007, 2008 by The University of York
# Copyright 2008 by Bernhard Lohkamp
# Copyright 2007 by Paul Emsley
# Copyright 2007 by The University of Oxford
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

# we only have the menu if shelxl is in PATH (and if we have the
# python menubar)
if (have_coot_python and find_exe("shelxl", "PATH")):

    menu = coot_menubar_menu("SHELX")

    def shelx_refine_func():

       window = gtk.Window(gtk.WINDOW_TOPLEVEL)
       vbox = gtk.VBox(False, 0)
       hbox = gtk.HBox(False, 0)
       go_button = gtk.Button("  Refine  ")
       cancel_button = gtk.Button("  Cancel  ")
       entry_hint_text = "HKL data filename \n(leave blank for default)"
       chooser_hint_text = " Choose molecule for SHELX refinement  "
       h_sep = gtk.HSeparator()

       window.add(vbox)
       options_menu_mol_list_pair = generic_molecule_chooser(vbox, chooser_hint_text)
       entry = file_selector_entry(vbox, entry_hint_text)

       def shelx_delete_event(*args):
           window.destroy()
           return False

       def shelx_refine_go_funcn_event(*args):
           import operator
           txt = entry.get_text()
           imol = get_option_menu_active_molecule(option_menu_list_pair)
           if (operator.isNumberType(imol)):
               if (len(txt) == 0):
                   shelxl_refine(imol)
               else:
                   shelxl_refine(imol, txt)
           window.destroy()
           return False

       go_button.connect("clicked", shelx_refine_go_funcn_event)
       cancel_button.connect("clicked", shelx_delete_event)

       vbox.pack_start(h_sep, False, False, 2)
       vbox.pack_start(hbox, False, False, 2)
       hbox.pack_start(go_button, True, False, 0)
       hbox.pack_start(cancel_button, True, False, 0)
       window.show_all()
       
    add_simple_coot_menu_menuitem(menu, "SHELXL Refine...",
                                  lambda func:shelx_refine_func())


    def shelx_read_project_func(*args):
        def shelx_delete_event(*args):
            window.destroy()
            return False

        def shelx_read_go_funcn_event(*args):
            file_name = entry.get_text()
            read_shelx_project(file_name)
            window.destroy()
            return False

        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        hbox = gtk.HBox(False, 0)
        vbox = gtk.VBox(True, 0)
        h_sep = gtk.HSeparator()
        go_button = gtk.Button("  Read Project   ")
        cancel_button = gtk.Button("  Cancel   ")

        entry = file_selector_entry(vbox, " Project Name: ")

        cancel_button.connect("clicked", shelx_delete_event)
        go_button.connect("clicked", shelx_read_go_funcn_event)

        window.add(vbox)

        hbox.pack_start(go_button, False, False, 2)
        hbox.pack_start(cancel_button, False, False, 2)
        vbox.pack_start(h_sep, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)
        window.show_all()
         
    add_simple_coot_menu_menuitem(menu, "Read SHELX Project...",
                                  lambda func: shelx_read_project_func())


    add_simple_coot_menu_menuitem(menu, "Read LST file...",
               lambda func: generic_chooser_and_file_selector("Model Corresponding to LST file: ",
                                                              valid_model_molecule_qm,
                                                              "LST file",
                                                              "",
                                                              lambda imol, lst_file_name: read_shelx_lst_file(lst_file_name, imol)))


    add_simple_coot_menu_menuitem(menu, "Add SHELXL instruction...",
               lambda func: generic_chooser_and_entry("Add new SHELXL command to model:",
                                                      "SHELX instruction:",
                                                      "",
                                                      lambda imol, text: add_shelx_string_to_molecule(imol, text)))
