
# This is an example script for hover previews. It is not used
# by coot.

import os
import coot
import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gdk
import coot_gui_api

def make_ligand_preview_menu(ligand_names, on_ligand_selected):
   """
   Create a toolbar menu button with ligand names that show SVG previews on hover.

   ligand_names: list of 3-letter codes, e.g. ["ATP", "GOL", "SO4", "MPD"]
   on_ligand_selected: callback function taking (ligand_name: str)
   """

   image_dir = os.path.join(os.path.expanduser("~"), ".cache", "Coot", "monomer-images")

   menu_button = Gtk.MenuButton(label="Ligands")

   popover = Gtk.Popover()
   menu_button.set_popover(popover)

   vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=2)
   vbox.set_margin_top(6)
   vbox.set_margin_bottom(6)
   vbox.set_margin_start(6)
   vbox.set_margin_end(6)
   popover.set_child(vbox)

   # The preview popover — a secondary popover attached to the main popover
   preview_popover = Gtk.Popover()
   preview_popover.set_autohide(False)
   preview_popover.set_can_focus(False)

   preview_picture = Gtk.Picture()
   preview_picture.set_size_request(200, 200)
   preview_popover.set_child(preview_picture)
   preview_popover.set_position(Gtk.PositionType.RIGHT)
   preview_popover.set_parent(popover)

   current_hover = [None]  # track which item is hovered

   def on_enter(controller, x, y, ligand_name, button):
      image_path = os.path.join(image_dir, ligand_name + ".svg")
      if not os.path.exists(image_path):
         image_path = os.path.join(image_dir, ligand_name + ".png")
      if os.path.exists(image_path):
         preview_picture.set_filename(image_path)
         # Re-parent the preview popover to the hovered button
         if preview_popover.get_parent() is not None:
            preview_popover.unparent()
         preview_popover.set_parent(button)
         preview_popover.set_position(Gtk.PositionType.RIGHT)
         preview_popover.popup()
         current_hover[0] = ligand_name
      else:
         preview_popover.popdown()
         current_hover[0] = None

   def on_leave(controller, ligand_name):
      if current_hover[0] == ligand_name:
         preview_popover.popdown()
         current_hover[0] = None

   for name in ligand_names:
      button = Gtk.Button(label=name)
      button.set_has_frame(False)
      button.add_css_class("flat")

      # Click handler
      button.connect("clicked", lambda btn, n=name: _on_clicked(btn, n, popover, on_ligand_selected))

      # Hover handlers
      motion = Gtk.EventControllerMotion()
      motion.connect("enter", lambda ctrl, x, y, n=name, b=button: on_enter(ctrl, x, y, n, b))
      motion.connect("leave", lambda ctrl, n=name: on_leave(ctrl, n))
      button.add_controller(motion)

      vbox.append(button)

   coot_gui_api.main_toolbar().append(menu_button)
   return menu_button


def _on_clicked(button, ligand_name, popover, callback):
   popover.popdown()
   callback(ligand_name)


# --- Set up the menu ---

def on_ligand_selected(ligand_name):
   print(f"Selected ligand: {ligand_name}")
   # placeholder - e.g. coot.get_monomer(ligand_name) or similar

ligands = ["ATP", "GOL", "SO4", "MPD", "CQ8"]
make_ligand_preview_menu(ligands, on_ligand_selected)

