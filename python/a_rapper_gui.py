import pygtk, gtk, pango

def rapper_process():
    #needed?
    pass

def rapper_it(imol, chain_id, start_resno, end_resno, sequence, number_of_models):

    imol_map = imol_refinement_map()
    if (not valid_map_molecule_qm(imol_map)):
        print "No valid map molecule given (possibly ambiguous)"
    else:
        str = "//" + chain_id + "/" + str(start_no) + \
              "-" + str(end_resno)
        frag_mol = new_molecule_by_atom_selection(imol, str)
        fragment_pdb = "coot-rapper-fragment-in.pdb"
        map_file = "coot-rapper.map"
        write_pdb_file(frag_mol, fragment_pdb)
        export_map(imol_map, map_file)
        print "running rapper: ", imol, chain_id, start_resno, end_resno, sequence, number_of_models

        # run rapper, not decided yet how

def stop_rapper():
    print "stopping rapper process..."

def cancel_dialog_func(widget, window):

    # do something more clever if the rapper process is still running
    window.destroy()

#
def a_rapper_gui():

    def rapper_go_func(*args):
        from types import IntType
        imol = get_option_menu_active_molecule(option_menu_pdb, model_mol_list)
        if (not (type(imol) is IntType)):
            print "bad active model"
        else:
            chain_text = entry_chain.get_text()
            start_resno_text = entry_start_resno.get_text()
            end_resno_text = entry_end_resno.get_text()
            text_sequence_buffer = text_sequence.get_buffer()
            sequence = text_sequence_buffer.get_text()
            number_of_models_text = entry_models.get_text()
            start_resno = int(start_resno_text)
            end_resno = int(end_resno_text)
            number_of_models = int(number_of_models_text)

            if (not ((type(start_resno) is IntType) and
                     (type(end_resno) is IntType) and
                     (type(number_of_models) is IntType))):
                print "Something incomprehensible: ", start_resno, end_resno, number_of_models
            else:
                rapper_it(imol, chain_text, start_resno, end_resno, sequence, number_of_models)
            

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    vbox = gtk.VBox(False, 2)
    pdb_hbox = gtk.HBox(False, 2)
    chain_hbox = gtk.HBox(False, 2)
    start_resno_hbox = gtk.HBox(False, 2)
    end_resno_hbox = gtk.HBox(False, 2)
    sequence_vbox = gtk.VBox(False, 2)
    models_hbox = gtk.HBox(False, 2)
    buttons_hbox = gtk.HBox(False, 2)
    h_sep = gtk.HSeparator()
    stop_button = gtk.Button("  Stop  ")
    go_button = gtk.Button("  Go  ")
    cancel_button = gtk.Button("  Cancel  ")

    # pdb
    label_pdb = gtk.Label("Rebuild fragment of Model: ")
    option_menu_pdb = gtk.combo_box_new_text()

    # chain
    label_chain = gtk.Label("Chain")
    entry_chain = gtk.Entry()

    # start resno
    label_start_resno = gtk.Label("Starting Res Number: ")
    entry_start_resno = gtk.Entry()

    # end resno
    label_end_resno = gtk.Label("     End Res Number: ")
    entry_end_resno = gtk.Entry()

    # sequence
    
    # sequence as it currently is:
    sequence_as_is_check_button = gtk.CheckButton("As is")

    # sequence as test
    label_sequence = gtk.Label("Sequence")
    text_sequence = gtk.TextView()

    # number of rapper models
    label_models = gtk.Label("Number of Models: ")
    entry_models = gtk.Entry()

    #
    model_mol_list = fill_option_menu_with_coordinates_mol_options(option_menu_pdb)

    pdb_hbox.pack_start(label_pdb, False, False, 2)
    pdb_hbox.pack_start(option_menu_pdb, False, False, 2)
    chain_hbox.pack_start(label_chain, False, False, 2)
    chain_hbox.pack_start(entry_chain, False, False, 2)
    start_resno_hbox.pack_start(label_start_resno, False, False, 2)
    start_resno_hbox.pack_start(entry_start_resno, False, False, 2)
    end_resno_hbox.pack_start(label_end_resno, False, False, 2)
    end_resno_hbox.pack_start(entry_end_resno, False, False, 2)

    sequence_vbox.pack_start(label_sequence, False, False, 2)
    sequence_vbox.pack_start(sequence_as_is_check_button, False, False, 2)
    sequence_vbox.pack_start(text_sequence, False, False, 2)

    models_hbox.pack_start(label_models, False, False, 2)
    models_hbox.pack_start(entry_models, False, False, 2)

    window.set_title("A Rapper GUI")
    vbox.pack_start(pdb_hbox, False, False, 2)
    vbox.pack_start(chain_hbox, False, False, 2)
    vbox.pack_start(start_resno_hbox, False, False, 2)
    vbox.pack_start(end_resno_hbox, False, False, 2)
    vbox.pack_start(sequence_vbox, False, False, 2)
    vbox.pack_start(models_hbox, False, False, 2)
    vbox.pack_start(h_sep, False, False, 2)

    buttons_hbox.pack_start(go_button, True, False, 6)
    buttons_hbox.pack_start(stop_button, True, False, 6)
    buttons_hbox.pack_start(cancel_button, True, False, 6)

    vbox.pack_start(buttons_hbox, True, False, 2)
    window.add(vbox)

    sequence_as_is_check_button.set_active(True)

    stop_button.connect("clicked", stop_rapper)

    cancel_button.connect("clicked", cancel_dialog_func, window)

    go_button.connect("clicked", rapper_go_func)

    window.show_all()

a_rapper_gui()
    
    
    
