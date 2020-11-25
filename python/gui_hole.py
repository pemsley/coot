
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

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    vbox = gtk.VBox(False, 0)
    hbox = gtk.VBox(False, 0)
    hbox_pos_buttons = gtk.HBox(False, 0)
    hbox_calc_cancel_buttons = gtk.HBox(False, 0)
    start_button = gtk.Button("  Set Start Point  ")
    end_button   = gtk.Button("  Set End Point  ")
    h_sep = gtk.HSeparator()
    calculate_button = gtk.Button("Calculate")
    cancel_button = gtk.Button("Cancel")
    hole_export_entry = gtk.Entry()
    export_text = gtk.Label("Export surface dots to File: ")
    export_hbox = gtk.HBox(False, 0)
    combobox = coot_gui.generic_molecule_chooser(hbox, "HOLE-ify molecule: ")

    window.add(vbox)
    hbox_pos_buttons.pack_start(start_button, False, False, 6)
    hbox_pos_buttons.pack_start(end_button, False, False, 6)
    hbox_calc_cancel_buttons.pack_start(calculate_button, False, False, 6)
    hbox_calc_cancel_buttons.pack_start(cancel_button, False, False, 6)
    export_hbox.pack_start(export_text, False, False, 6)
    export_hbox.pack_start(hole_export_entry, False, False, 6)
    vbox.pack_start(hbox, True, True, 6)
    vbox.pack_start(hbox_pos_buttons, True, True, 6)
    vbox.pack_start(export_hbox, True, True, 6)
    vbox.pack_start(h_sep)
    vbox.pack_start(hbox_calc_cancel_buttons, True, True, 6)

    hole_export_entry.set_text("hole_surface_dots.dat")

    def start_button_cb(*args):
        global start_pos
        start_pos = coot_utils.rotation_centre()
        print("Start pos set to:", start_pos)
        status_bar_pos(start_pos, "start")
        
    start_button.connect("clicked", start_button_cb)

    def end_button_cb(*args):
        global end_pos
        end_pos = coot_utils.rotation_centre()
        print("End pos set to:", end_pos)
        status_bar_pos(end_pos, "end")
        
    end_button.connect("clicked", end_button_cb)

    def delete_event(*args):
        window.destroy()
        return False
    
    def calculate_button_cb(*args):
        global start_pos, end_pos
        imol = combobox_to_molecule_number(combobox)
        if isinstance(imol, int):
            print(start_pos, end_pos)
            if not isinstance(start_pos, list):
                coot.add_status_bar_text("Start position not set")
            else:
                if not isinstance(end_pos, list):
                    coot.add_status_bar_text("End position not set")
                else:
                    # main?
                    print("hole", imol, start_pos, end_pos)
                    # get this from an entry ideally.
                    export_dots_file_name = hole_export_entry.get_text()
                    colour_map_multiplier = 1
                    colour_map_offset = 0
                    hole_args = [imol] + start_pos + end_pos + \
                                [colour_map_multiplier, colour_map_offset, 1,
                                 True, export_dots_file_name]
                    print("BL DEBUG:: hole_args", hole_args)
                    coot.hole(*hole_args)
                    delete_event()
                    
    calculate_button.connect("clicked", calculate_button_cb)

    cancel_button.connect("clicked", delete_event)

    window.show_all()

