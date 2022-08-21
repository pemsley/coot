
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk
import coot
import coot_gui_api
import coot_utils
import coot_gui

def chain_refine(imol, ch_id):
    residue_list = coot_utils.residues_in_chain(imol, ch_id)
    coot.refine_residues_py(imol, residue_list)

# this might be a coot_gui function now
def coot_menubar_menu(s):
   menu = Gtk.Menu()
   menuitem = Gtk.MenuItem(s)
   menuitem.set_submenu(menu)
   main_menubar = coot_gui_api.main_menubar()
   main_menubar.append(menuitem)
   menuitem.show()
   return menu

def add_module_refine():

    menu_refine = coot_menubar_menu("Refine")

    def all_atom_refine_func():
        active_atom = coot.active_residue_py()
        if active_atom:
            imol = active_atom[0]
            residue_list = coot_utils.all_residues(imol)
            coot.refine_residues_py(imol, residue_list)

    def chain_refine_func():
        active_atom = coot.active_residue_py()
        if active_atom:
            imol = active_atom[0]
            ch_id = active_atom[1]
            chain_refine(imol, ch_id)

    def refine_active_fragment_func():
        active_atom = coot.active_residue_py()
        if active_atom:
            imol = active_atom[0]
            ch_id = active_atom[1]
            residue_list = coot.linked_residues_py(res_spec, imol, 1.7)
            coot.refine_residues_py(imol, residue_list)
        
    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "All-atom Refine", lambda arg: all_atom_refine_func())
    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "Chain Refine",    lambda arg: chain_refine_func())
    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "Refine Fragment", lambda arg: chain_refine_func())

    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "Contact Dots On",  lambda arg: coot.set_do_coot_probe_dots_during_refine(1))
    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "Contact Dots Off", lambda arg: coot.set_do_coot_probe_dots_during_refine(0))

    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "Rotamer Markup On",   lambda arg: coot.set_show_intermediate_atoms_rota_markup(1))
    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "Rotamer Markup Off",  lambda arg: coot.set_show_intermediate_atoms_rota_markup(0))

    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "Refine Rotamers On",  lambda arg: coot.set_refine_rotamers(1))
    coot_gui.add_simple_coot_menu_menuitem(menu_refine, "Refine Rotamers Off", lambda arg: coot.set_refine_rotamers(0))
