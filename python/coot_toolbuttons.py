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

import pygtk, gtk, pango
import time

  
have_coot_python = False
try: 
  import coot_python
  have_coot_python = True
except:
  print """BL WARNING:: could not import coot_python module!!
Some things, esp. extensions, may be crippled!"""

# put all coot svgs in the default icon set
def register_coot_icons():
  import glob
  iconfactory = gtk.IconFactory()
  stock_ids   = gtk.stock_list_ids()
  pixbuf_dir = os.getenv('COOT_PIXMAPS_DIR')
  patt = os.path.normpath(pixbuf_dir + '/*.svg')
  coot_icon_filename_ls = glob.glob(patt)
  icon_info_ls = []
  for full_name in coot_icon_filename_ls:
    name = os.path.basename(full_name)
    icon_info_ls.append([name, full_name])
  for stock_id, filename in icon_info_ls:
    # only load image files when our stock_id is not present
    if (stock_id not in stock_ids) and (not 'phenixed' in filename):
      pixbuf = gtk.gdk.pixbuf_new_from_file(filename)
      iconset = gtk.IconSet(pixbuf)
      iconfactory.add(stock_id, iconset)
  iconfactory.add_default()

if (have_coot_python):

  # an example function for the toolbutton (not sure if it should be here!)
  def stereo_mono_toggle():
    display_state = stereo_mode_state()
    if (display_state == 0):
      hardware_stereo_mode()
    else:
      mono_mode()

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
                #print "icon", icon, icon_stock, icon_filename
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
                    if (check_button_label == button_label):
                      new_toolbutton = coot_toolbar_button(button_label, callback_function, icon)
                      # save
                      save_toolbar_to_init_file(button_label, callback_function, icon)
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

      def make_icons_model():
        import os, glob
        stock_icon_ls = gtk.stock_list_ids()
        # coot icons
        pixbuf_dir = os.getenv('COOT_PIXMAPS_DIR')
        patt = os.path.normpath(pixbuf_dir + '/*.svg')
        coot_icon_filename_ls = glob.glob(patt)

        model = gtk.ListStore(gtk.gdk.Pixbuf, str, str)
        for icon_filename in coot_icon_filename_ls:
          if (not 'phenixed' in icon_filename):
            icon = os.path.basename(icon_filename)
            pixbuf = gtk.gdk.pixbuf_new_from_file(icon_filename)
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
        
      def save_function(*args):
        print "Save me"

      window = gtk.Window(gtk.WINDOW_TOPLEVEL)
      window.set_title("Toolbar Selection")
      window.set_default_size(600, 450)
      scrolled_win = gtk.ScrolledWindow()
      scrolled_win.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_ALWAYS)
      vbox = gtk.VBox(False, 2)
      frame_vbox = gtk.VBox(False, 2)
      h_sep = gtk.HSeparator()
      buttons_hbox = gtk.HBox(False, 2)
      cancel_button = gtk.Button("  Cancel  ")
      go_button = gtk.Button("   Add   ")
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
            check_button_label = item[0]
            callback_function  = item[1]
            description        = item[2]
          
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
            button_icon = None
            for toolbar_button_label in toolbar_label_list():
              if (check_button_label == toolbar_button_label[0]):
                check_button.set_active(True)
                try:
                  button_icon = toolbar_button_label[1].get_stock_id()
                  for index in range(len(icon_model)):
                    icon_iter = icon_model.get_iter(index)
                    stock_id = icon_model.get_value(icon_iter, 1)
                    if (stock_id == button_icon):
                      icon_combobox.set_active_iter(icon_iter)
                      break
                except:
                  # no button icon available
                  #print "BL DEBUG:: we dont have a button icon"
                  pass
            left_hbox.pack_start(check_button, False, False, 3)
            left_hbox.pack_end(label, False, False, 3)
            right_hbox.pack_start(icon_label, False, False, 3)
            right_hbox.pack_start(icon_combobox, False, False, 3)
            hbox.pack_start(left_hbox, False, False, 3)
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
      

    def remove_toolbar_button_gui():
      def remove_toolbar_button(entry_text):
        print "remove button", entry_text
        for toolbar_child in coot_main_toolbar.get_children():
          button_label = toolbar_child.get_label()
          if (button_label == entry_text):
            coot_main_toolbar.remove(toolbar_child)
            remove_toolbar_from_init_file(button_label)
            break
      generic_single_entry("Remove toolbar button", "button label", "Remove", remove_toolbar_button)

    def toolbar_hide_text():
      coot_main_toolbar.set_style(gtk.TOOLBAR_ICONS)

    def toolbar_show_text():
      coot_main_toolbar.set_style(gtk.TOOLBAR_BOTH_HORIZ)

    def show_pop_up_menu(widget, event):
      if (event.button == 3):
        create_show_pop_menu([["Manage Buttons/Add a new Button", make_toolbar_button_gui],
                              ["Remove a Button", remove_toolbar_button_gui],
                              ["Hide text (only icons)", toolbar_hide_text],
                              ["Show text", toolbar_show_text]],
                             event)

    coot_main_toolbar = coot_python.main_toolbar()
    coot_main_toolbar.connect("button-press-event", show_pop_up_menu)

# save a toolbar button to ~/.coot-preferences/coot_toolbuttons.py
def save_toolbar_to_init_file(button_label, callback_function, icon):

  save_str = "coot_toolbar_button(\"" + button_label + "\", \"" + callback_function
  if icon:
    save_str += ("\", \"" + icon)
  save_str += "\")"
  
  home = 'HOME'
  if (os.name == 'nt'):
    home = 'COOT_HOME'
  filename = os.path.join(os.getenv(home), ".coot-preferences", "coot_toolbuttons.py")
  remove_line_containing_from_file(["coot_toolbar_button", button_label], filename)
  save_string_to_file(save_str, filename)

# remove a toolbar from  ~/.coot-preferences/coot_toolbuttons.py
def remove_toolbar_from_init_file(button_label):
  home = 'HOME'
  if (os.name == 'nt'):
    home = 'COOT_HOME'
  filename = os.path.join(os.getenv(home), ".coot-preferences", "coot_toolbuttons.py")
  remove_str_ls = ["coot_toolbar_button", button_label]
  if (os.path.isfile(filename)):
    remove_line_containing_from_file(remove_str_ls, filename)


# returns a list with pre-defined toolbar-functions
# format:
# [[Group1,[toolbarname1, callbackfunction1, description1], [toolbarname2, callbackfunction2, description],...], [Group2]....]
# OBS: callbackfunction is a string!
def list_of_toolbar_functions():
  ls = [["Display",
         ["Stereo/Mono", "stereo_mono_toggle()", "Toggle between Stereo and Mono view"],
         ["test", "rotation_centre()", "test function"]],
        ["Refinement",
         ["Refine active residue", "refine_active_residue()", ""]],
        ["NMR",[]],
        ["EM",[]],
        ["Sidechains/Alignment",[]]]
  return ls

  
