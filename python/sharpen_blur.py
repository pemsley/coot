
from __future__ import print_function
import gtk

def sharpen_blur_map_gui():

    def on_ok_button_clicked(button, combobox, map_mol_list, entry_1, entry_2, check_button, check_button_for_refinement_map, window):
        print("click", button)
        imol_map = coot_gui.get_option_menu_active_molecule(combobox, map_mol_list)
        t1 = entry_1.get_text()
        blur_factor = float(t1)
        if check_button.get_active():
            t2 = entry_2.get_text()
            resample_factor = float(t2)
            imol_new = coot.sharpen_blur_map_with_resampling(imol_map, blur_factor, resample_factor)
        else:
            imol_new = coot.sharpen_blur_map(imol_map, blur_factor)
        print("check_button_for_refinement_map active", check_button_for_refinement_map.get_active())
        if check_button_for_refinement_map.get_active():
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

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    entry_hint_text_2 = "Factor:"
    default_entry_1_text = "20"
    default_entry_2_text = "1.3"

    label = gtk.Label(chooser_label)
    vbox = gtk.VBox(False, 2)
    hbox_for_sharpen  = gtk.HBox(False, 0)
    hbox_for_resample = gtk.HBox(False, 0)
    entry_1 = gtk.Entry();
    entry_2 = gtk.Entry();
    entry_label_1 = gtk.Label(entry_hint_text_1)
    entry_label_2 = gtk.Label(entry_hint_text_2)
    hbox_buttons = gtk.HBox(True, 2)
    combobox = gtk.combo_box_new_text()
    ok_button    = gtk.Button(" Make Map ")
    cancel_button = gtk.Button(" Cancel ")
    h_sep = gtk.HSeparator()
    map_mol_list = coot_gui.fill_option_menu_with_map_mol_options(combobox)
    check_button = gtk.CheckButton("Resample")
    check_button_for_refinement_map = gtk.CheckButton("Make the new map the Refinement Map")

    window.set_title("Coot: Sharpen/Blur Map")
    window.set_default_size(400, 100)
    window.add(vbox)

    vbox.pack_start(label,    False, False, 2)
    vbox.pack_start(combobox, False, False, 2)
    vbox.pack_start(hbox_for_sharpen,  False, False, 2)
    vbox.pack_start(hbox_for_resample, False, False, 2)
    vbox.pack_start(check_button_for_refinement_map, False, False, 2)
    vbox.pack_start(h_sep, False, False, 2)
    vbox.pack_start(hbox_buttons, False, False, 2)

    hbox_for_sharpen.pack_start(entry_label_1, False, False, 2)
    hbox_for_sharpen.pack_start(entry_1, False, False, 2)
    hbox_for_resample.pack_start(check_button,  False, False, 2)
    hbox_for_resample.pack_start(entry_label_2, False, False, 2)
    hbox_for_resample.pack_start(entry_2, False, False, 2)
    hbox_buttons.pack_start(ok_button, False, False, 6)
    hbox_buttons.pack_start(cancel_button, False, False, 6)
    entry_1.set_size_request(52, -1)
    entry_2.set_size_request(52, -1)
    entry_1.set_text(default_entry_1_text)
    entry_2.set_text(default_entry_2_text)

    check_button_for_refinement_map.set_active(True)
    entry_label_2.set_sensitive(False)
    entry_2.set_sensitive(False)

    cancel_button.connect("clicked", lambda button: window.destroy())
    ok_button.connect("clicked", on_ok_button_clicked, combobox, map_mol_list, entry_1, entry_2, check_button, check_button_for_refinement_map, window)
    check_button.connect("toggled", on_check_button_toggled, entry_2, entry_label_2)
                                                  

    window.show_all()

if False:
    if True:
        if coot_gui_api.main_menubar():
            menu = coot_gui.coot_menubar_menu("Cryo-EM")
            if menu:
                coot_gui.add_simple_coot_menu_menuitem(menu, "Sharpen/Blur Map",
                                              lambda func: sharpen_blur_map_gui())
