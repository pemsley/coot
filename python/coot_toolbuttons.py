# Copyright 2007 by Bernhard Lohkamp
# Copyright 2006, 2007 by The University of York

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

import pygtk
import gtk
import pango
import time
import os

  
# put all coot svgs in the default icon set
def register_coot_icons():
  import glob
  iconfactory = gtk.IconFactory()
  stock_ids   = gtk.stock_list_ids()
  pixbuf_dir = os.getenv('COOT_PIXMAPS_DIR')
  if (not pixbuf_dir):
    pixbuf_dir = os.path.join(get_pkgdatadir(), "pixmaps")
  patt = os.path.normpath(pixbuf_dir + '/*.svg')
  coot_icon_filename_ls = glob.glob(patt)
  icon_info_ls = []
  for full_name in coot_icon_filename_ls:
    name = os.path.basename(full_name)
    icon_info_ls.append([name, full_name])
  for stock_id, filename in icon_info_ls:
    # only load image files when our stock_id is not present
    if ((stock_id not in stock_ids) and not
        ('phenixed' in filename)):
      if os.path.isfile(filename):
        pixbuf = gtk.gdk.pixbuf_new_from_file(filename)
        iconset = gtk.IconSet(pixbuf)
        iconfactory.add(stock_id, iconset)
  iconfactory.add_default()

# adds a SeparatorToolItem to the coot_main_toolbar (by default at the
# last position). Return the separator or False. If there is a
# separator in the last position, dont add it
#
def add_coot_toolbar_separator():
  """adds a SeparatorToolItem to the coot_main_toolbar (by default at the
  last position). There is no return value currently.
  
  """
  
  try:
    coot_main_toolbar = coot_python.main_toolbar()
  except:
    print """BL ERROR:: coot_python module not available!!
    So we cannot make toolbar_separator!"""
    return False

  # main body
  toolbar_items = coot_main_toolbar.get_children()
  if (type(toolbar_items[-1]) != gtk.SeparatorToolItem):
    # only add if last one is not already separator
    sep = gtk.SeparatorToolItem()
    coot_main_toolbar.insert(sep, -1)
    sep.show()
    return sep
  else:
    return False

  
if (have_coot_python):

  if coot_python.main_toolbar():

    register_coot_icons()

    def activate_menuitem(widget, item):
      try:
        apply(item[1])
      except:
        print "BL INFO:: unable to execute function", item[1]

    # takes list with [["item_name", func],...]
    def create_show_pop_menu(item_n_funcn_list, event):
      menu = gtk.Menu()
      for item in item_n_funcn_list:
        item_name = item[0]
        #item_func = item[1]
        menuitem = gtk.MenuItem(item_name)
        menu.append(menuitem)
        menuitem.connect("activate", activate_menuitem, item)
      menu.show_all()
      menu.popup(None, None, None, event.button, event.time)

    def make_icons_model():
      import os, glob
      stock_icon_ls = gtk.stock_list_ids()
      # coot icons
      pixbuf_dir = os.getenv('COOT_PIXMAPS_DIR')
      if (not pixbuf_dir):
        pixbuf_dir = os.path.join(get_pkgdatadir(), "pixmaps")
      patt = os.path.normpath(pixbuf_dir + '/*.svg')
      coot_icon_filename_ls = glob.glob(patt)

      model = gtk.ListStore(gtk.gdk.Pixbuf, str, str)
      for icon_filename in coot_icon_filename_ls:
        if os.path.isfile(icon_filename):
          icon = os.path.basename(icon_filename)
          pixbuf = gtk.gdk.pixbuf_new_from_file_at_size(icon_filename, 16, 16)
          model.append([pixbuf, icon, icon_filename])

      # build in default gtk icons
      icon_theme = gtk.icon_theme_get_default()
      for icon in stock_icon_ls:
        try:
          pixbuf = icon_theme.load_icon(icon, 16, gtk.ICON_LOOKUP_USE_BUILTIN)
          model.append([pixbuf, icon, None])
        except:
          pass
      return model

    def icon_selection_combobox(model):
        
      combobox = gtk.ComboBox()
      combobox.set_wrap_width(10)
      combobox.set_model(model)
      crpx = gtk.CellRendererPixbuf()
      #crt  = gtk.CellRendererText()
      combobox.pack_start(crpx, False)
      combobox.add_attribute(crpx, 'pixbuf', 0)
      combobox.show_all()
      return combobox


    def make_toolbar_button_gui():

      def delete_event(*args):
        window.destroy()
        return False

      def go_function(*args):
        save_toolbuttons_qm = save_button.get_active()
        #print "do something, maybe, save state", save_toolbuttons_qm
        # iterate thru the widgets and get the state of check button and icon combobox
        frames = frame_vbox.get_children()
        for frame in frames:
          inner_vbox = frame.get_children()[0]
          hboxes = inner_vbox.get_children()
          for hbox in hboxes:
            left_hbox, right_hbox = hbox.get_children()
            check_button, description = left_hbox.get_children()
            label, combobox = right_hbox.get_children()
            model = combobox.get_model()
            button_label = check_button.get_label()
            if (check_button.get_active()):
              # make new button
              iter = combobox.get_active_iter()
              icon = None
              if (iter):
                icon, icon_stock, icon_filename = model.get(iter, 0, 1, 2)
                # print "icon", icon, icon_stock, icon_filename
                if (icon_stock):
                  icon = icon_stock
              toolbar_ls = list_of_toolbar_functions()
              for group in toolbar_ls:
                group_name = group[0]
                group_items = group[1:len(group)]
                for item in group_items:
                  if (len(item) >0):
                    check_button_label = item[0]
                    callback_function  = item[1]
                    description        = item[2]
                    toggle_button_item = False
                    use_button_item    = False
                    # is this foolproof? We assume an icon on i3
                    if len(item) > 4:
                      toggle_button_item = item[4]
                      if len(item) > 5:
                        use_button_item = item[5]
                    if (check_button_label == button_label):
                      new_toolbutton = coot_toolbar_button(button_label,
                                                           callback_function,
                                                           icon_name=icon,
                                                           tooltip=description,
                                                           toggle_button=toggle_button_item,
                                                           use_button=use_button_item)
                      # save
                      if save_toolbuttons_qm:
                        save_toolbar_to_init_file(button_label,
                                                  callback_function,
                                                  icon=icon,
                                                  tooltip=description,
                                                  toggle_button=toggle_button_item,
                                                  use_button=use_button_item)
                      break
            else:
              # remove an existing button?
              # check if the button is in the existing button list
              for toolbar in toolbar_label_list():
                if (button_label == toolbar[0]):
                  coot_python.main_toolbar().remove(toolbar[1])
                  # remove from save
                  remove_toolbar_from_init_file(button_label)

        window.destroy()
        return False
          
      def save_function(*args):
        print "Save me"

      window = gtk.Window(gtk.WINDOW_TOPLEVEL)
      window.set_title("Toolbar Selection")
      window.set_default_size(700, 450)
      scrolled_win = gtk.ScrolledWindow()
      scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)
      vbox = gtk.VBox(False, 2)
      frame_vbox = gtk.VBox(False, 2)
      h_sep = gtk.HSeparator()
      buttons_hbox = gtk.HBox(False, 2)
      cancel_button = gtk.Button("  Cancel  ")
      go_button = gtk.Button("   Apply   ")
      save_button = gtk.CheckButton(" Save to preferences?")
      save_button.set_active(True)
      icon_model = make_icons_model()
            
      toolbar_ls = list_of_toolbar_functions()

      for group in toolbar_ls:
        group_name = group[0]
        group_items = group[1:len(group)]
        
        inner_vbox = gtk.VBox(False, 2)
        frame = gtk.Frame(group_name)
        frame.add(inner_vbox)

        for item in group_items:
          if (len(item) >0):
            suggested_icon = False
            check_button_label = item[0]
            callback_function  = item[1]
            description        = item[2]
            if (len(item) > 3):
              suggested_icon     = item[3]
          
            def check_button_callback(*args):
              # print "check button toggled"
              # we actually dont do anything and check later for the status?
              pass

            hbox = gtk.HBox(False, 0)
            left_hbox = gtk.HBox(False, 0)
            right_hbox = gtk.HBox(False, 0)
            check_button = gtk.CheckButton(check_button_label)
            label = gtk.Label(description)
            icon_label = gtk.Label("Add icon?")
            icon_combobox = icon_selection_combobox(icon_model)
            button_icon = False

            all_toolbar_labels = map(lambda x: x[0], toolbar_label_list())
            if (check_button_label in all_toolbar_labels):
                check_button.set_active(True)
                # now try to get the icon
                try:
                  ls = toolbar_label_list()
                  button_icon = ls[all_toolbar_labels.index(check_button_label)][1].get_stock_id()
                except:
                  # no icon found
                  pass

            if (not button_icon):
              button_icon = suggested_icon

            # set the icon
            if (button_icon):
              for index in range(len(icon_model)):
                icon_iter = icon_model.get_iter(index)
                stock_id = icon_model.get_value(icon_iter, 1)
                if (stock_id == button_icon):
                  icon_combobox.set_active_iter(icon_iter)
                  break
            
            left_hbox.pack_start(check_button, True, True, 3)
            left_hbox.pack_end(label, True, True, 3)
            right_hbox.pack_start(icon_label, False, False, 3)
            right_hbox.pack_start(icon_combobox, False, False, 3)
            hbox.pack_start(left_hbox, True, True, 3)
            hbox.pack_end(right_hbox, False, False, 3)
            inner_vbox.pack_start(hbox, False, False, 3)
            check_button.connect("toggled", check_button_callback)

        frame_vbox.pack_start(frame, False, False, 3)

      frame_vbox.set_border_width(3)
      scrolled_win.add_with_viewport(frame_vbox)
      buttons_hbox.pack_start(go_button, True, False, 6)
      buttons_hbox.pack_start(cancel_button, True, False, 6)
      vbox.pack_start(scrolled_win, True, True, 0)
      vbox.pack_start(h_sep, False, False, 3)
      vbox.pack_start(save_button, False, False, 3)
      vbox.pack_start(buttons_hbox, False, False, 3)

      cancel_button.connect("clicked", delete_event)
      go_button.connect("clicked", go_function)
      save_button.connect("toggled", save_function)

      window.add(vbox)
      window.show_all()
      
    # simple GUI to remove a toolbar
    def remove_toolbar_button_gui():
      
      def remove_toolbar_button(entry_text):
        print "remove button", entry_text
        for toolbar_child in coot_main_toolbar.get_children():
          if (type(toolbar_child) == gtk.ToolButton or
              type(toolbar_child) == gtk.ToggleToolButton):
            button_label = toolbar_child.get_label()
            if (button_label == entry_text):
              coot_main_toolbar.remove(toolbar_child)
              remove_toolbar_from_init_file(button_label)
              break
      generic_single_entry("Remove toolbar button",
                           "button label",
                           "Remove",
                           remove_toolbar_button)

    # GUI to add user-defined toolbars
    # uses Assistant or a simple GUI, depending on PyGTK version
    def add_toolbar_button_gui():
      major, minor, micro = gtk.pygtk_version
      if (major >= 2 and minor >= 10):
        # have assistant
        print "BL DEBUG:: run assistant"
        add_toolbar_button_assistant()
      else:
        # no assistant available -> simple GUI
        add_toolbar_button_simple_gui()
        print "BL DEBUG:: no assistant, so simple GUI"
      

    def toolbar_hide_text():
      coot_main_toolbar.set_style(gtk.TOOLBAR_ICONS)

    def toolbar_show_text():
      coot_main_toolbar.set_style(gtk.TOOLBAR_BOTH_HORIZ)

    ####################################
    # an assistant to add a toolbutton #
    ####################################
    def add_toolbar_button_assistant():

      def cb_close(assistant):
        assi.destroy()
        return False

      def cb_name_entry_key_press(entry, event):
        assi.set_page_complete(entry.get_parent(),
                               len(entry.get_text()) > 0)

      def cb_func_entry_key_press(entry, event):
        assi.set_page_complete(entry.get_parent(),
                               len(entry.get_text()) > 0)

      def cb_entry_key_press(entry, event):
        assi.set_page_complete(entry.get_parent(),
                               len(entry.get_text()) > 0)

      def cb_apply(assistant):
        name_str = name_entry.get_text()
        func_str = func_entry.get_text()
        save_qm  = radiobutton_yes.get_active()
        icon_iter = icon_combobox.get_active_iter()
        icon = None
        if (icon_iter):
          icon, icon_stock, icon_filename = icon_model.get(icon_iter, 0, 1, 2)
          if (icon_stock):
            icon = icon_stock

        coot_toolbar_button(name_str, func_str, icon)
        if (save_qm):
          save_toolbar_to_init_file(name_str, func_str, icon)
        
      assi = gtk.Assistant()
      assi.set_title("Toolbutton Assistant")

      assi.connect('delete_event', cb_close)
      assi.connect('close', cb_close)
      assi.connect('cancel', cb_close)
      assi.connect('apply', cb_apply)

      # Construct page 0 (enter the button name)
      vbox = gtk.VBox(False, 5)
      vbox.set_border_width(5)
      vbox.show()
      assi.append_page(vbox)
      assi.set_page_title(vbox, 'Create a new Coot toolbutton')
      assi.set_page_type(vbox, gtk.ASSISTANT_PAGE_CONTENT)

      label = gtk.Label("Please enter the name for the new toolbutton...")
      label.set_line_wrap(True)
      label.show()
      vbox.pack_start(label, True, True, 0)

      name_entry = gtk.Entry()
      name_entry.connect("key-press-event", cb_name_entry_key_press)
      name_entry.show()
      vbox.pack_end(name_entry)
      # check if name exists!!!!

      # Construct page 1 (enter the function)
      vbox = gtk.VBox(False, 5)
      vbox.set_border_width(5)
      vbox.show()
      assi.append_page(vbox)
      assi.set_page_title(vbox, 'Set the function')
      assi.set_page_type(vbox, gtk.ASSISTANT_PAGE_CONTENT)

      label = gtk.Label("Please enter the python function... (e.g. refine_active_residue())")
      label.set_line_wrap(True)
      label.show()
      vbox.pack_start(label, True, True, 0)

      func_entry = gtk.Entry()
      func_entry.connect("key-press-event", cb_func_entry_key_press)
      func_entry.show()
      vbox.pack_start(func_entry, True, True, 0)
      # add advanced option? FIXME (for callable func and args!?)
      #adv_button = gtk.Button("Advanced options...")
      #adv_button.show()
      #vbox.pack_start(adv_button, False, False, 0)

      # Construct page 2 (save?)
      vbox = gtk.VBox(False, 5)
      vbox.set_border_width(5)
      vbox.show()
      assi.append_page(vbox)
      assi.set_page_title(vbox, 'Save to preferences?')
      assi.set_page_type(vbox, gtk.ASSISTANT_PAGE_CONTENT)

      label = gtk.Label("Do you want to save the button in your preferences?")
      label.set_line_wrap(True)
      label.show()
      vbox.pack_start(label, True, True, 0)

      radiobutton_yes = gtk.RadioButton(None, "Yes")
      radiobutton_no  = gtk.RadioButton(radiobutton_yes, "No")
      radiobutton_yes.set_active(True)
      radiobutton_yes.show()
      radiobutton_no.show()
      vbox.pack_start(radiobutton_yes, True, True, 0)
      vbox.pack_start(radiobutton_no,  True, True, 0)
      assi.set_page_complete(radiobutton_yes.get_parent(), True)

      # Construct page 3 (select icon)
      vbox = gtk.VBox(False, 5)
      vbox.set_border_width(5)
      vbox.show()
      assi.append_page(vbox)
      assi.set_page_title(vbox, 'Select an icon')
      assi.set_page_type(vbox, gtk.ASSISTANT_PAGE_CONTENT)

      label = gtk.Label("Please select an icon or leave as it is...")
      label.set_line_wrap(True)
      label.show()
      vbox.pack_start(label, True, True, 0)

      hbox = gtk.HBox(False, 5)
      icon_model = make_icons_model()
      icon_combobox = icon_selection_combobox(icon_model)
      icon_combobox.show()
      hbox.show()
      hbox.pack_start(icon_combobox, True, False, 2)
      vbox.pack_end(hbox)
      assi.set_page_complete(label.get_parent(), True)

      # Final page
      # As this is the last page needs to be of page_type
      # gtk.ASSISTANT_PAGE_CONFIRM or gtk.ASSISTANT_PAGE_SUMMARY
      label = gtk.Label('Thanks for using the toolbutton assistent!')
      label.set_line_wrap(True)
      label.show()
      assi.append_page(label)
      assi.set_page_title(label, 'Finish')
      assi.set_page_complete(label, True)
      assi.set_page_type(label, gtk.ASSISTANT_PAGE_CONFIRM)

      assi.show()

    def add_toolbar_button_simple_gui():

      def button_func(*args):
        # dummy dont need to do anythin here
        pass

      def do_func(*args):
        name_str = args[0]
        func_str = args[1]
        save_qm  = args[2]
        icon_iter = icon_combobox.get_active_iter()
        icon = None
        if (icon_iter):
          icon, icon_stock, icon_filename = icon_model.get(icon_iter, 0, 1, 2)
          if (icon_stock):
            icon = icon_stock

        coot_toolbar_button(name_str, func_str, icon)
        if (save_qm):
          save_toolbar_to_init_file(name_str, func_str, icon)
        
      entry_widget = generic_double_entry("Button Name:", "Python Function:",
                                          "test", "centre_of_mass(0)",
                                          "Save to Preferences?", button_func,
                                          "Create", do_func,
                                          return_widget=True)

      hbox = gtk.HBox(True, 0)
      label = gtk.Label("Icon:")
      icon_model = make_icons_model()
      icon_combobox = icon_selection_combobox(icon_model)
      hbox.pack_start(label, True, False, 0)
      hbox.pack_start(icon_combobox, True, False, 0)
      children = entry_widget.get_children()
      vbox = children[0]
      vbox.pack_start(hbox, True, False, 2)
      vbox.reorder_child(hbox, 3) 
      hbox.show_all()

    def show_pop_up_menu(widget, event):
      if (event.button == 3):
        create_show_pop_menu([["Manage Buttons (add, delete Buttons)", make_toolbar_button_gui],
                              ["Add a user-defined Button", add_toolbar_button_gui],
                              ["Remove a Button", remove_toolbar_button_gui],
                              ["Hide Text (only icons)", toolbar_hide_text],
                              ["Show Text", toolbar_show_text]],
                             event)

    coot_main_toolbar = coot_python.main_toolbar()
    coot_main_toolbar.connect("button-press-event", show_pop_up_menu)

# save a toolbar button to ~/.coot-preferences/coot_toolbuttons.py
#
def save_toolbar_to_init_file(button_label, callback_function,
                              icon=None, tooltip=None,
                              toggle_button=False, use_button=False):

  # so far only for string
  save_str = "coot_toolbar_button(\"" + button_label + "\", "
  if type(callback_function) is StringType:
    save_str += ("\"" + callback_function + "\"")
  else:
    # no lists yet. FIXME
    save_str += (callback_function.__name__)
  if icon:
    save_str += (", icon_name=\"" + str(icon) + "\"")
  if tooltip:
    save_str += (", tooltip=\"" + str(tooltip) + "\"")
  if toggle_button:
    save_str += (", toggle_button=" + str(toggle_button))
  if use_button:
    save_str += (", use_button=" + str(use_button))
  save_str += ")"
  
  home = os.getenv('HOME')
  if (not home and os.name == 'nt'):
    home = os.getenv('COOT_HOME')
  if not home:
    print "BL ERROR:: could not find a home directory"
  else:
    filename = os.path.join(home, ".coot-preferences", "coot_toolbuttons.py")
    remove_line_containing_from_file(["coot_toolbar_button", button_label],
                                     filename)
    save_string_to_file(save_str, filename)

    
# remove a toolbar from  ~/.coot-preferences/coot_toolbuttons.py
#
def remove_toolbar_from_init_file(button_label):
  home = os.getenv('HOME')
  if (not home and os.name == 'nt'):
    home = os.getenv('COOT_HOME')
  if not home:
    print "BL ERROR:: could not find a home directory"
  else:    
    filename = os.path.join(home, ".coot-preferences", "coot_toolbuttons.py")
    remove_str_ls = ["coot_toolbar_button", button_label]
    if (os.path.isfile(filename)):
      remove_line_containing_from_file(remove_str_ls, filename)


# returns a list with pre-defined toolbar-functions (stock-id is optional,
# so are a few others - almost consistent with order in coot_toolbar_button
# function)
# format:
# [[Group1,[toolbarname1, callbackfunction1, description1],
#          [toolbarname2, callbackfunction2, description2(tooltip),
#           stock-id2, toggle_button, use_button],...],
#  [Group2]....]
# OBS: callbackfunction is a string! not any more exclusively
#
def list_of_toolbar_functions():
  ls = [["Display",
         ["Stereo/Mono", "stereo_mono_toggle()", "Toggle between Stereo and Mono view", "stereo-view.svg"],
         ["Side-by-side/Mono", "side_by_side_stereo_mono_toggle()", "Toggle between Side-by-Side Stereo and Mono view", "stereo-view.svg"],
         ["Zalman Stereo/mono", "zalman_stereo_mono_toggle()", "Toggle between Zalman Stereo and Mono view", "stereo-view.svg"],
         ["Swap Stereo", "switch_stereo_sides()", "Change left and right stereo image", "undo-1.svg"],
         ["Full screen", toggle_full_screen,
          "Switch between full screen and window mode", "gtk-fullscreen",
          True, True],
          ['Hydrogens', toggle_hydrogen_display,
           "Toggle to show (or not) hydrogens",
           "delete.svg", True, True],
         ],
        ["Refinement",
         ["Sphere Refine", "sphere_refine()", "RSR around active residue", "reset-view.svg"],
         ["Sphere Refine +", "sphere_refine_plus()", "RSR around active residue +/- 1 residue", "reset-view.svg"],
         ["Repeat Refine Zone", "repeat_refine_zone()", "Repeat (so that I don't need to click the atoms again)", "rrz.svg"],
         ["Tandem Refine", "refine_tandem_residues()", "Refine a tandem 7-residue range", "something.svg"],
         ["Sphere Regularization", "sphere_regularize()", "Regularize around active residue", "reset-view.svg"],
         ["Sphere Regularization +", "sphere_regularize_plus()", "Regularize around active residue +/- 1 residue", "reset-view.svg"],
         ["Refine residue", "refine_active_residue()", "RSR active residue"],
         ["Reset B", "reset_b_factor_active_residue()", "Reset the B-Factor of active Residue"],
         ["Add Alt Conf", "altconf()", "Add alternative conformation", "add-alt-conf.svg"],
         ["Change Alt Conf Occ", "select_atom_alt_conf_occ_gui()", "Change occupancies for alternative conformations", "add-alt-conf.svg"],
         ["Edit BB", "setup_backbone_torsion_edit(1)", "Edit Backbone Torsion Angle", "flip-peptide.svg"],
         ['Torsion Gen.', "setup_torsion_general(1)", "Torsion General (after O function)", "edit-chi.svg"],
         ['Backrub Rotamers', toggle_backrub_rotamers,
          "Toggle use (or not) of backrub rotamers",
          "auto-fit-rotamer.svg", True, True],
         ['Cis<->Trans', "do_cis_trans_conversion_setup(1)", "Convert peptide: cis->trans or trans->cis", "flip-peptide.svg"],
         ["Run Refmac", "wrapped_create_run_refmac_dialog()", "Launch Refmac for Refinement", "azerbaijan.svg"]],
        ["Validation",
         ["Interactive dots", toggle_interactive_probe_dots,
          "Show dots after refinement and for chi/rotamer changes",
          "probe-clash.svg", True, True],
         ["Local probe dots", "probe_local_sphere_active_atom()",
          "Show probe dots for active atom in 4A radius.",
          "probe-clash.svg"]
         ],
        ["Building",
         ["Choose Undo Molecule", "show_set_undo_molecule_chooser()", "Choose Undo Molecule", "undo-1.svg"],
         ["Ligand Builder", "start_ligand_builder_gui()", "Ligand Builder", "go-to-ligand.svg"],
         ["Split Water", "split_active_water()", "Split water in two and refine", "add-water.svg"],
         ["Find Waters", "wrapped_create_find_waters_dialog()", "Find water molecules in map", "add-water.svg"],
         ["Build NA", "find_nucleic_acids_local(6.0)", "Find nucleic acids locally using Cootilus", "dna.svg"]
         ],
        ["NMR",[]],
        ["EM",[]],
        ["Sidechains/Alignment",[]]]
  return ls

# to make combobox
#coot_toolbar_combobox("MC restraints", ["no","aH","bS"], [set_secondary_structure_restraints_type], tooltip="Change refinement restraints for secondary structure")
  
