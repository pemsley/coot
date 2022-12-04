import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
import coot
import coot_utils
import coot_gui

def hole_ify():
    global start_pos, end_pos
    start_pos = False
    end_pos = False

    def status_bar_pos(position, pos_type):
        s = "Hole %s point set: (%6.2f %6.2f %6.2f)" %(pos_type,
                                                       position[0],
                                                       position[1],
                                                       position[2])
        coot.add_status_bar_text(s)

    window = Gtk.Window()
    vbox = Gtk.VBox(False, 0)
    hbox = Gtk.VBox(False, 0)
    hbox_pos_buttons = Gtk.HBox(False, 0)
    hbox_calc_cancel_buttons = Gtk.HBox(False, 0)
    start_button = Gtk.Button("Set Start Point")
    end_button   = Gtk.Button("Set End Point")
    calculate_button = Gtk.Button("Calculate")
    cancel_button = Gtk.Button("Cancel")
    # something_holder = coot_gui.generic_molecule_chooser(hbox, "HOLE-ify molecule: ")
    combobox = Gtk.ComboBox()
    # Gtk.widget_show(combobox)
    combobox_mol_items = coot_gui.make_store_for_model_molecule_combobox()
    combobox.set_model(combobox_mol_items)
    renderer_text = Gtk.CellRendererText()
    combobox.pack_start(renderer_text, True)
    combobox.add_attribute(renderer_text, "text", 1)
    combobox.set_active(0)

    window.add(vbox)
    hbox_pos_buttons.pack_start(start_button, True, True, 6)
    hbox_pos_buttons.pack_start(end_button, True, True, 6)
    hbox_calc_cancel_buttons.pack_start(calculate_button, True, True, 6)
    hbox_calc_cancel_buttons.pack_start(cancel_button, True, True, 6)
    vbox.pack_start(combobox, True, True, 6)
    vbox.pack_start(hbox, True, True, 6)
    vbox.pack_start(hbox_pos_buttons, True, True, 6)
    vbox.pack_start(hbox_calc_cancel_buttons, True, True, 6)

    def get_molecule():
        tree_iter = combobox.get_active_iter()
        imol = -1
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol = it[0]
        return imol

    def start_button_cb(*args):
        global start_pos
        start_pos = coot_utils.rotation_centre()
        print("INFD:: Start pos set to:", start_pos)
        status_bar_pos(start_pos, "start")

    start_button.connect("clicked", start_button_cb)

    def end_button_cb(*args):
        global end_pos
        end_pos = coot_utils.rotation_centre()
        print("INFO:: End pos set to:", end_pos)
        status_bar_pos(end_pos, "end")
        
    end_button.connect("clicked", end_button_cb)

    def delete_event(*args):
        window.destroy()
        return False
    
    def calculate_button_cb(*args):
        global start_pos, end_pos

        # imol = coot_gui.get_option_menu_active_molecule(*option_menu_and_model_mol_list)

        imol = get_molecule()

        if True:
            print("INFO:: Hole positions", start_pos, end_pos)
            if not isinstance(start_pos, list):
                coot.add_status_bar_text("Start position not set")
            else:
                if not isinstance(end_pos, list):
                    coot.add_status_bar_text("End position not set")
                else:
                    print("INFO:: hole", imol, start_pos, end_pos)
                    colour_map_multiplier = 1
                    colour_map_offset = 0
                    hole_args = [imol] + start_pos + end_pos + \
                                [colour_map_multiplier, colour_map_offset, 1, 1, "coot-hole-dots.table"]
                    output_file_name = "coot-hole-dots.table"
                    n_runs = 4
                    coot.hole(imol,
                              start_pos[0], start_pos[1], start_pos[2],
                              end_pos[0], end_pos[1], end_pos[2],
                              colour_map_multiplier, colour_map_offset, n_runs, True,
                              output_file_name)
                    window.destroy()

    calculate_button.connect("clicked", calculate_button_cb)

    cancel_button.connect("clicked", delete_event)

    window.show_all()

# if coot_python.main_menubar():
#     menu = coot_gui.coot_menubar_menu("Test Hole")
#     coot_gui.add_simple_coot_menu_menuitem(menu, "Test",
#             lambda func: hole_ify()
#             )
