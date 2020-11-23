
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
    start_button = gtk.Button("Set Start Point")
    end_button   = gtk.Button("Set End Point")
    calculate_button = gtk.Button("Calculate")
    cancel_button = gtk.Button("Cancel")
    option_menu_and_model_mol_list = coot_gui.generic_molecule_chooser(hbox,
                                                              "HOLE-ify molecule: ")

    window.add(vbox)
    hbox_pos_buttons.pack_start(start_button, True, True, 6)
    hbox_pos_buttons.pack_start(end_button, True, True, 6)
    hbox_calc_cancel_buttons.pack_start(calculate_button, True, True, 6)
    hbox_calc_cancel_buttons.pack_start(cancel_button, True, True, 6)
    vbox.pack_start(hbox, True, True, 6)
    vbox.pack_start(hbox_pos_buttons, True, True, 6)
    vbox.pack_start(hbox_calc_cancel_buttons, True, True, 6)

    def start_button_cb(*args):
        global start_pos
        start_pos = coot_utils.rotation_centre()
        print "Start pos set to:", start_pos
        status_bar_pos(start_pos, "start")
        
    start_button.connect("clicked", start_button_cb)

    def end_button_cb(*args):
        global end_pos
        end_pos = coot_utils.rotation_centre()
        print "End pos set to:", end_pos
        status_bar_pos(end_pos, "end")
        
    end_button.connect("clicked", end_button_cb)

    def delete_event(*args):
        window.destroy()
        return False
    
    def calculate_button_cb(*args):
        global start_pos, end_pos
        imol = coot_gui.get_option_menu_active_molecule(*option_menu_and_model_mol_list)
        if isinstance(imol, int):
            print start_pos, end_pos
            if not isinstance(start_pos, list):
                coot.add_status_bar_text("Start position not set")
            else:
                if not isinstance(end_pos, list):
                    coot.add_status_bar_text("End position not set")
                else:
                    # main?
                    print "hole", imol, start_pos, end_pos
                    colour_map_multiplier = 1
                    colour_map_offset = 0
                    hole_args = [imol] + start_pos + end_pos + \
                                [colour_map_multiplier, colour_map_offset, 1, 1 ]
                    coot.hole(*hole_args)
                    delete_event()
                    
    calculate_button.connect("clicked", calculate_button_cb)

    cancel_button.connect("clicked", delete_event)

    window.show_all()

if (have_coot_python):
    if coot_python.main_menubar():
        menu = coot_gui.coot_menubar_menu("Test Hole")
        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Test",
            lambda func: hole_ify()
            )
