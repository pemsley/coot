
def new_molecule_with_nudged_residues(imol, residue_spec,
                                      residue_delta, nudge_by):

    imol_new = coot.copy_molecule(imol)
    chain_id = res_spec_utils.residue_spec_to_chain_id(residue_spec)
    resno_start = res_spec_utils.residue_spec_to_res_no(residue_spec) - residue_delta
    resno_end = res_spec_utils.residue_spec_to_res_no(residue_spec) + residue_delta

    if coot_utils.debug():
        print("imol:", imol)
        print("residue_spec:", residue_spec)
        print("residue_delta:", residue_delta)
        print("resno_start:", resno_start)
        print("resno_end:", resno_end)
        print("nudge_by:", nudge_by)

    status = coot.nudge_residue_sequence(imol_new, chain_id, resno_start, resno_end,
                                    nudge_by, 1)

    if (status == 0):
        # fail
        s = "Failed to nudge around " + chain_id + " " + \
            str(coot_utils.residue_spec_to_res_no(residue_spec))
        coot.add_status_bar_text(s)
        coot.close_molecule(imol_new)
        return -1  # return a bad new molecule id
    else:
        return imol_new

def nudge_residues_gui(imol, residue_spec):

    def delete_event(*args):
        window.destroy()
        return False

    def change_nudge(*args):
        v = adj.value
        v_int = int(v) # rounding!?
        rdt = entry.get_text()
        try:
            rd = int(rdt)
        except:
            print("BL WARNING:: could not convert %s to a number, set to 1 then." %rdt)
            # or shall we bail!?
            rd = 1
        residue_delta = rd

        imol_new = new_molecule_with_nudged_residues(imol, residue_spec,
                                                     residue_delta, v_int)
        # return imol_ne # why, for who?
  
    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    vbox = gtk.VBox(False, 4)
    hbox_0 = gtk.HBox(False, 4)
    hbox_1 = gtk.HBox(False, 4)
    hbox_2 = gtk.HBox(False, 4)
    label_1 = gtk.Label(" Nudge by ")
    label_2 = gtk.Label(" residues ")
    label_n = gtk.Label(" Nudge ")
    res_lab = " residues up and down from " + \
              res_spec_utils.residue_spec_to_chain_id(residue_spec) + " " + \
              str(coot_utils.residue_spec_to_res_no(residue_spec))
    label_a = gtk.Label(res_lab)
    m_lab = " Nudging residues from Molecule:\n   " + \
            str(imol) + ": " + \
            coot_utils.strip_path(coot.molecule_name(imol))
    label_m = gtk.Label(m_lab)
    entry = gtk.Entry()
    h_sep = gtk.HSeparator()
    adj = gtk.Adjustment(0., -30., 59., 0.01, 1., 29.)
    slider = gtk.HScale(adj)
    buttons_hbox = gtk.HBox(False, 4)
    cancel_button = gtk.Button(" Close ")
    residue_delta = 5 # but isnt it a variable?!

    window.set_title("Nudge Residues")
    vbox.pack_start(hbox_0, False, False, 4)
    vbox.pack_start(hbox_1, False, False, 4)
    vbox.pack_start(hbox_2, False, False, 4)
    vbox.pack_start(h_sep, False, False, 4)
    vbox.pack_start(buttons_hbox, False, False, 4)
    hbox_0.pack_start(label_m, False, False, 4)
    hbox_1.pack_start(label_n, False, False, 2)
    hbox_1.pack_start(entry, False, False, 2)
    hbox_1.pack_start(label_a, False, False, 2)
    hbox_2.pack_start(label_1, False, False, 4)
    hbox_2.pack_start(slider, True, True, 4)
    hbox_2.pack_start(label_2, False, False, 4)
    buttons_hbox.pack_end(cancel_button, False, False, 4)
    window.add(vbox)
    window.set_size_request(350, 200)
    entry.set_text("5")
    entry.set_usize(30, -1)
    slider.set_digits(0)

    cancel_button.connect("clicked", delete_event)

    adj.connect("value_changed", change_nudge)
    
    window.show_all()
    
    
