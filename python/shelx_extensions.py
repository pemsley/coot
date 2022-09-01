# shelx_extensions.py
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
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numbers
import coot
import coot_utils
import coot_gui

def add_module_shelx():
    # we only have the menu if shelxl is in PATH (and if we have the
    # python menubar)
    if not coot_utils.find_exe("shelxl", "PATH"):
        coot.info_dialog("WARNING:: Cannot find shelxl.\n\nSome SHELX plugin functions may not be working")
    if have_coot_python:
        menu = coot_gui.coot_menubar_menu("SHELX")

        def combobox_to_molecule_number(combobox):
            imol = -1
            tree_iter = combobox.get_active_iter()
            if tree_iter is not None:
                model = combobox.get_model()
                it = model[tree_iter]
                imol = it[0]
            return imol


        def shelx_refine_func():

            window = Gtk.Window()
            vbox = Gtk.VBox(False, 0)
            hbox = Gtk.HBox(False, 0)
            go_button = Gtk.Button("  Refine  ")
            cancel_button = Gtk.Button("  Cancel  ")
            entry_hint_text = "HKL data filename \n(leave blank for default)"
            chooser_hint_text = " Choose molecule for SHELX refinement  "
            h_sep = Gtk.HSeparator()

            window.add(vbox)
            # option_menu_mol_list_pair = coot_gui.generic_molecule_chooser(vbox, chooser_hint_text)

            combobox_items = make_store_for_molecule_combobox(molecule_filter_function)
            combobox = Gtk.ComboBox.new_with_model(combobox_items)
            renderer_text = Gtk.CellRendererText()
            if len(combobox_items) > 0:
                combobox.set_active(0)
                combobox.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
            combobox.pack_start(renderer_text, True)
            combobox.add_attribute(renderer_text, "text", 1)
            entry = coot_gui.file_selector_entry(vbox, entry_hint_text)

            def shelx_delete_event(*args):
                window.destroy()
                return False

            def shelx_refine_go_funcn_event(*args):
                import operator
                txt = entry.get_text()
                imol = combobox_to_molecule_number(combobox)
                if (isinstance(imol, numbers.Number)):
                    editable_shelx_gui(imol, txt)
                window.destroy()
                return False

            go_button.connect("clicked", shelx_refine_go_funcn_event)
            cancel_button.connect("clicked", shelx_delete_event)

            vbox.append(h_sep)
            vbox.append(hbox)
            hbox.append(go_button)
            hbox.append(cancel_button)
            window.show_all()

        coot_gui.add_simple_coot_menu_menuitem(
            menu, "SHELXL Refine...",
            lambda func: shelx_refine_func())


        def shelx_read_project_func(*args):
            def shelx_delete_event(*args):
                window.destroy()
                return False

            def shelx_read_go_funcn_event(*args):
                file_name = entry.get_text()
                shelx.read_shelx_project(file_name)
                window.destroy()
                return False

            window = Gtk.Window(Gtk.WINDOW_TOPLEVEL)
            hbox = Gtk.HBox(False, 0)
            vbox = Gtk.VBox(True, 0)
            h_sep = Gtk.HSeparator()
            go_button = Gtk.Button("  Read Project   ")
            cancel_button = Gtk.Button("  Cancel   ")

            entry = coot_gui.file_selector_entry(vbox, " Project Name: ")

            cancel_button.connect("clicked", shelx_delete_event)
            go_button.connect("clicked", shelx_read_go_funcn_event)

            window.add(vbox)

            hbox.append(go_button)
            hbox.append(cancel_button)
            vbox.append(h_sep)
            vbox.append(hbox)
            window.show_all()

        coot_gui.add_simple_coot_menu_menuitem(
            menu, "Read SHELX Project...",
            lambda func: shelx_read_project_func())


        coot_gui.add_simple_coot_menu_menuitem(
            menu, "Read LST file...",
            lambda func: coot_gui.generic_chooser_and_file_selector("Model Corresponding to LST file: ",
                                                           coot_utils.valid_model_molecule_qm,
                                                           "LST file",
                                                           "",
                                                           lambda imol, lst_file_name: shelx.read_shelx_lst_file(lst_file_name, imol)))


        coot_gui.add_simple_coot_menu_menuitem(
            menu, "Add SHELXL instruction...",
            lambda func: coot_gui.generic_chooser_and_entry("Add new SHELXL command to model:",
                                                   "SHELX instruction:",
                                                   "",
                                                   lambda imol, text: coot.add_shelx_string_to_molecule(imol, text)))

def shelx_ins_strings(imol):

    ins_tmp_file = "coot-tmp.ins"
    coot.write_shelx_ins_file(imol, ins_tmp_file)
    lines = []
    try:
        fin = open(ins_tmp_file, 'r')
        lines = fin.readlines()
        fin.close()
    except:
        print("INFO:: problems reading file", ins_tmp_file)
    return lines

def shelxl_refine_gui(imol, hkl_file_name_maybe=False):

    def shelx_delete_event(*args):
        window.destroy()
        return False

    def shelx_refine_go_funcn_event(*args):
        import operator
        start, end = textbuffer.get_bounds()
        txt = textbuffer.get_text(start, end)
        shelx.shelxl_refine_primitive(imol, txt, hkl_file_name_maybe)
        window.destroy()
        return False

    window = Gtk.Window()
    window.set_size_request(500, 500)
    vbox = Gtk.VBox(False, 2)
    buttons_hbox = Gtk.HBox(False, 2)
    scrolled_win = Gtk.ScrolledWindow()
    text = Gtk.TextView()
    textbuffer = text.get_buffer()
    text.set_editable(True)
    refine_button = Gtk.Button("  Refine  ")
    cancel_button = Gtk.Button("  Cancel  ")

    window.add(vbox)
    scrolled_win.add(text)
    buttons_hbox.append(refine_button)
    buttons_hbox.append(cancel_button)
    scrolled_win.set_policy(Gtk.POLICY_AUTOMATIC, Gtk.POLICY_ALWAYS)
    vbox.append(buttons_hbox)
    vbox.pack_start(scrolled_win)

    text_strings = shelx_ins_strings(imol)
    bg_col = "#c0e6c0"

    text.modify_base(Gtk.STATE_NORMAL, Gtk.gdk.color_parse(bg_col))
    text.modify_text(Gtk.STATE_NORMAL, Gtk.gdk.color_parse("black"))
    for string in text_strings:
        end = textbuffer.get_end_iter()
        textbuffer.insert(end, string)

    refine_button.connect("clicked", shelx_refine_go_funcn_event)
    cancel_button.connect("clicked", shelx_delete_event)

    window.show_all()


def editable_shelx_gui(imol, hklin_file_name):

    def shelx_delete_event(*args):
        window.destroy()
        return False

    def shelx_refine_go_funcn_event(*args):
        from types import StringType
        start, end = textbuffer.get_bounds()
        txt = textbuffer.get_text(start, end)
        if (type(txt) == StringType):
            hklin_file_info = False
            if (len(hklin_file_name) > 0):
                hklin_file_info = hklin_file_name
            shelx.shelxl_refine_primitive(imol, txt, hklin_file_info)
        window.destroy()
        return False

    window = Gtk.Window(Gtk.WINDOW_TOPLEVEL)
    text = Gtk.TextView()
    textbuffer = text.get_buffer()
    scrolled_win = Gtk.ScrolledWindow()
    cancel_button = Gtk.Button("  Cancel  ")
    run_button = Gtk.Button("  Run  ")
    vbox = Gtk.VBox(False, 0)
    buttons_hbox = Gtk.HBox(False, 0)

    window.set_size_request(450, 400)
    window.add(vbox)
    scrolled_win.add(text)
    vbox.set_border_width(5)
    buttons_hbox.append(run_button)
    buttons_hbox.append(cancel_button)
    vbox.append(buttons_hbox)
    vbox.append(scrolled_win)
    scrolled_win.set_policy(Gtk.POLICY_AUTOMATIC, Gtk.POLICY_ALWAYS)
    text.set_editable(True)

    shelx_ins_list = get_shelx_ins_list(imol)
    bg_col = "white"

    text.modify_base(Gtk.STATE_NORMAL, Gtk.gdk.color_parse(bg_col))
    text.modify_text(Gtk.STATE_NORMAL, Gtk.gdk.color_parse("black"))
    for string in shelx_ins_list:
        end = textbuffer.get_end_iter()
        textbuffer.insert(end, string)

    run_button.connect("clicked", shelx_refine_go_funcn_event)
    cancel_button.connect("clicked", shelx_delete_event)
    window.show_all()

# seems to be the same...
get_shelx_ins_list = shelx_ins_strings
