
import coot
import gi
gi.require_version('Gtk', '4.0')
import coot_gui_api
import coot_gui
import coot_utils
import fitting
import random
from gi.repository import Gtk

def sharpen_blur_map_gui():

    def on_ok_button_clicked(button, combobox, check_button_for_resample, check_button_for_refinement_map, entry_for_blur, entry_for_resample, window):
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol_map = it[0]
            state_1 = check_button_for_resample.get_active()
            state_2 = check_button_for_refinement_map.get_active()
            print("state_1", state_1)
            print("state_2", state_2)
            imol_new = -1
            t1 = entry_for_blur.get_text()
            blur_factor = float(t1)
            if state_1:
                t2 = entry_for_resample.get_text()
                resample_factor = float(t2)
                imol_new = coot.sharpen_blur_map_with_resampling(imol_map, blur_factor, resample_factor)
            else:
                imol_new = coot.sharpen_blur_map(imol_map, blur_factor)
            if state_2:
                print("set_imol_refinement_map() with ", imol_new)
                coot.set_imol_refinement_map(imol_new)

        window.destroy()

    def on_check_button_toggled(check_button, entry_2, label_2):
        print("toggled", check_button)
        if check_button.get_active():
            entry_2.set_sensitive(True)
            label_2.set_sensitive(True)
        else:
            entry_2.set_sensitive(False)
            label_2.set_sensitive(False)

    chooser_label = "Map for Sharpen/Blur"
    entry_hint_text_1 =  "Sharpen/Blur:"

    window = Gtk.Window()
    entry_hint_text_2 = "Factor:"
    default_entry_1_text = "20"
    default_entry_2_text = "1.3"

    label = Gtk.Label(chooser_label)
    vbox = Gtk.VBox(False, 2)
    hbox_for_sharpen  = Gtk.HBox(False, 0)
    hbox_for_resample = Gtk.HBox(False, 0)
    entry_1 = Gtk.Entry();
    entry_2 = Gtk.Entry();
    entry_label_1 = Gtk.Label(entry_hint_text_1)
    entry_label_2 = Gtk.Label(entry_hint_text_2)
    hbox_buttons = Gtk.HBox(True, 2)
    # combobox = Gtk.combo_box_new_text()
    ok_button    = Gtk.Button(" Make Map ")
    cancel_button = Gtk.Button(" Cancel ")
    h_sep = Gtk.HSeparator()
    # map_mol_list = fill_option_menu_with_map_mol_options(combobox)
    # map_mol_list = fill_option_menu_with_map_mol_options(combobox)
    combobox_items = coot_gui.make_store_for_map_molecule_combobox()
    combobox = Gtk.ComboBox.new_with_model(combobox_items)
    renderer_text = Gtk.CellRendererText()
    if len(combobox_items) > 0:
        combobox.set_active(0)
    combobox.set_entry_text_column(1)
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)

    # combobox.connect("changed", on_mol_combobox_changed)
    check_button_for_resample = Gtk.CheckButton("Resample")
    check_button_for_refinement_map = Gtk.CheckButton("Make the new map the Refinement Map")

    window.set_title("Coot: Sharpen/Blur Map")
    window.set_default_size(400, 100)
    window.add(vbox)

    vbox.append(label)
    vbox.append(combobox)
    vbox.append(hbox_for_sharpen)
    vbox.append(hbox_for_resample)
    vbox.append(check_button_for_refinement_map)
    vbox.append(h_sep)
    vbox.append(hbox_buttons)

    hbox_for_sharpen.append(entry_label_1)
    hbox_for_sharpen.append(entry_1)
    hbox_for_resample.append(check_button_for_resample)
    hbox_for_resample.append(entry_label_2)
    hbox_for_resample.append(entry_2)
    hbox_buttons.append(ok_button)
    hbox_buttons.append(cancel_button)
    entry_1.set_size_request(52, -1)
    entry_2.set_size_request(52, -1)
    entry_1.set_text(default_entry_1_text)
    entry_2.set_text(default_entry_2_text)

    check_button_for_refinement_map.set_active(True)
    entry_label_2.set_sensitive(False)
    entry_2.set_sensitive(False)

    cancel_button.connect("clicked", lambda button: window.destroy())
    # ok_button.connect("clicked", on_ok_button_clicked, combobox, map_mol_list,
    # entry_1, entry_2, check_button, check_button_for_refinement_map, window)

    ok_button.connect("clicked", on_ok_button_clicked, combobox, check_button_for_resample, check_button_for_refinement_map, entry_1, entry_2, window)

    check_button_for_resample.connect("toggled", on_check_button_toggled, entry_2, entry_label_2)
                                                  

    window.show_all()

if False:
    if True:
        if coot_gui_api.main_menubar():
            menu = coot_gui.coot_menubar_menu("Cryo-EM")
            if menu:
                coot_gui.add_simple_coot_menu_menuitem(menu, "Sharpen/Blur Map",
                                              lambda func: sharpen_blur_map_gui())
