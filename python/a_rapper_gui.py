#    a_rapper_gui.py
#    Copyright (C) 2008  Bernhard Lohkamp, The University of York
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pygtk, gtk, pango

if (os.name == 'nt'):
    home_dir = os.getenv('COOT_HOME')
else:
    home_dir = os.getenv('HOME')
rapper_dir = os.path.join(os.path.normpath(home_dir), "rappermc")
rapper_command = "rapper"

def rapper_process():
    #needed?
    pass

def rename_dir_by_date(dir_name):
    import time
    date_str = time.strftime("%Y-%m-%d--%H-%M--%S")
    new_dir_name = dir_name + date_str
    # check if exists?
    os.rename(dir_name, new_dir_name)
    return new_dir_name

# return a string or False
#
def sequence_string(imol, chain_id, resno_start, resno_end):

    def all_chars_qm(ls):
        return all(map(lambda x: x.isalpha(), ls))

    if (not coot_utils.valid_model_molecule_qm(imol)):
        return False
    else:
        single_letter_code_list = []
        for resno in range(resno_start, resno_end + 1):
            res_name = residue_name(imol, chain_id, resno, "")
            single_letter_code_list.append(three_letter_code2single_letter(res_name))
        if (not all_chars_qm):
            print "bad sequence chars ", single_letter_code_list
            return False
        else:
            return "".join(single_letter_code_list)

def rapper_it(imol, chain_id, start_resno, end_resno, sequence, number_of_models):

    imol_map = imol_refinement_map()
    if (not coot_utils.valid_map_molecule_qm(imol_map)):
        print "No valid map molecule given (possibly ambiguous)"
    else:
        str = "//" + chain_id + "/" + str(start_no) + \
              "-" + str(end_resno)
        frag_mol = new_molecule_by_atom_selection(imol, str)
        fragment_pdb = "coot-rapper-fragment-in.pdb"
        rapper_out_pdb = "rapper_out.pdb"
        rapper_mode = "model-loops-benchmark" # maybe "ca-trace" perhaps
        # rapper_mode = "ca-trace"
        length_for_rapper = 3 # test
        if (type(sequence) is StringType):
            sequence_string = sequence
        else:
            sequence_string = "AAA"
        map_file = "coot-rapper.map"
        whole_pdb_file_name = "rapper-all-atoms.pdb"
        write_pdb_file(frag_mol, fragment_pdb)
        write_pdb_file(imol, whole_pdb_file_name)
        close_molecule(frag_mol)
        export_map(imol_map, map_file)
        print "running rapper: ", imol, chain_id, start_resno, end_resno, \
              sequence, number_of_models

        if (os.path.isfile("TESTRUNS")):
            rename_dir_by_date("TESTRUNS")
        
        rapper_command_line_args = [#"params.xml",
                                    rapper_mode,
                                    "--pdb", fragment_pdb,
                                    # "--framework", whole_pdb_file_name,
                                    "--pdb-out", rapper_out_pdb,
                                    "--map", map_file,
                                    "--models", str(number_of_models),
                                    "--start", str(start_resno),
                                    "--stop", str(end_resno),
                                    # "--length", str(length_for_rapper)
                                    "--seq", sequence_string,
                                    "--mainchain-restraint-threshold", "2.0",
                                    "--sidechain-centroid-restraint-threshold", "2.0",
                                    "--sidechain-mode", "smart",
                                    "--sidechain-radius-reduction", "0.75",
                                    "--enforce-mainchain-restraints", "true",
                                    "--enforce-sidechain-centroid-restraints", "true",
                                    "--edm-fit", "true",
                                    "--rapper-dir", rapper_dir]
        rapper_status = coot_utils.popen_command(rapper_command,
                                      rapper_command_line_args,
                                      [],
                                      "rapper.log",
                                      True)
        print "rapper_status:", rapper_status
        if (os.path.isfile("TESTRUNS")):
            new_dir_name = rename_dir_by_date("TESTRUNS")
            result_pdb_file_name = os.path.join(new_dir_name, "looptest-best.pdb")
            if (os.path.isfile(result_pdb_file_name)):
                read_pdb(result_pdb_file_name)
            else:
                info_dialog("RAPPER failed - no results")

        # run rapper, not decided yet how

def stop_rapper():
    print "stopping rapper process..."

def cancel_dialog_func(widget, window):

    # do something more clever if the rapper process is still running
    window.destroy()

# loop_building_tool is either 'rapper' or 'ARP/wARP'
#
def a_rapper_gui(loop_building_tool):

    def rapper_go_func(*args):
        from types import IntType
        imol = coot_gui.get_option_menu_active_molecule(option_menu_pdb, model_mol_list)
        if (not (type(imol) is IntType)):
            print "bad active model"
        else:
            chain_text = entry_chain.get_text()
            start_resno_text = entry_start_resno.get_text()
            end_resno_text = entry_end_resno.get_text()
            text_sequence_buffer = text_sequence.get_buffer()
            sequence = text_sequence_buffer.get_text()
            print "BL DEBUG:: input sequence is", sequence
            number_of_models_text = entry_models.get_text()
            start_resno = int(start_resno_text)
            end_resno = int(end_resno_text)
            number_of_models = int(number_of_models_text)

            if (not ((type(start_resno) is IntType) and
                     (type(end_resno) is IntType) and
                     (type(number_of_models) is IntType))):
                print "Something incomprehensible: ", start_resno, end_resno, number_of_models
            else:
                if (loop_building_tool == 'rapper'):
                    seq = sequence_string(imol, chain_text, start_resno, end_resno)
                    rapper_it(imol, chain_text, start_resno, end_resno,
                              seq, number_of_models)
                    
                elif (loop_building_tool == 'ARP/wARP'):
                    # not supprted yet
                    new_start = start_resno - 1
                    new_end   = end_resno + 1
                    seq = sequence_string(imo, chain_text,
                                          new_start, new_end)
                    arp_warp_it(imol, chain_text,
                                [start_resno, end_resno],
                                [new_start, new_end],
                                seq, number_of_models)
                else:
                    print "INFO:: invalid loop building argument", loop_building_tool

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

    # sequence as it currently is:
    sequence_as_is_check_button = gtk.CheckButton("As is")

    # sequence as text
    label_sequence = gtk.Label("Sequence")
    text_sequence = gtk.TextView()

    # number of rapper models
    label_models = gtk.Label("Number of Models: ")
    entry_models = gtk.Entry()

    #
    model_mol_list = coot_gui.fill_option_menu_with_coordinates_mol_options(option_menu_pdb)

    #text_sequence.set_usize(-1, 80) # ?? using text view!?
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

    entry_models.set_text("2")

    sequence_as_is_check_button.set_active(True)

    stop_button.connect("clicked", stop_rapper)

    cancel_button.connect("clicked", cancel_dialog_func, window)

    go_button.connect("clicked", rapper_go_func)

    window.show_all()

menu = coot_gui.coot_menubar_menu("Loop")
add_simple_coot_menu_menuitem(menu, "RAPPER...",
                              lambda func: a_rapper_gui('rapper'))

add_simple_coot_menu_menuitem(menu, "ARP/wARP Loopy...",
                              lambda func: a_rapper_gui('ARP/wARP'))

a_rapper_gui()
    
    
    
