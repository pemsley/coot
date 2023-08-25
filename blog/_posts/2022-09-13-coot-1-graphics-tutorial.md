---
layout: post
title:  "Coot Graphics Tutorials"
date: Tue 13 Sep 00:00:43 BST 2022
---

Lucrezia Catapano has made some introductory tutorials for Coot 1.

![Coot 1 Graphics Tutorial Image]({{"../../../images/Coot1-TheGraphics-FrontPageImage.png"}})

Here is the video to cover the basics:

[Coot 1 Tutorial: Basics](https://www.youtube.com/watch?v=Xhonm4K1y0c)

Here is the new (Sept 2022) graphics tutorial - map masking, Shader Preferences and scripting:

[Coot 1 Tutorial: Graphics Part 1](https://www.youtube.com/watch?v=tCpCmBTgt6s)

[Coot 1 Tutorial: Graphics Part 2](https://www.youtube.com/watch?v=OL5-aSZ-p9I)

Here is Lucrezia's script for the masking and colours:

{% highlight python %}

import coot
import coot_utils

# resize window
coot.set_graphics_window_size(880, 1000)

# read model and mapn
imol = coot.read_pdb('7p9b.pdb')
imol_map = coot.read_ccp4_map('emd_13261.map', 0)

# set map radius
coot.set_map_radius_em(99.0)

# mask map by chains
def make_masked_map_using_active_atom():
    active_atom = coot.active_residue_py()
    if active_atom:
        imol = active_atom[0]
        coot.make_masked_maps_split_by_chain(imol, coot.imol_refinement_map())
make_masked_map_using_active_atom()

# make maps as surface
coot_utils.solidify_maps()

# undisplay original map and model
coot.set_map_displayed(1, 0)
coot.set_mol_displayed(0, 0)

# set center and zoom
#coot.set_rotation_centre(166.403, 166.418, 166.409)
#coot.set_view_quaternion(1.000000, 0.000000, 0.000000, 0.000000)
coot.set_zoom(381.078)

# decrease clipping front
for i in range(12):
    coot.decrease_clipping_front()

# color maps
coot.set_map_hexcolour(2, '#A36352')
coot.set_map_hexcolour(3, '#A36352')
coot.set_map_hexcolour(4, '#D4992E')
coot.set_map_hexcolour(5, '#D4992E')
coot.set_map_hexcolour(6, '#B8A152')
coot.set_map_hexcolour(7, '#B8A152')
coot.set_map_hexcolour(8, '#B2B896' )
coot.set_map_hexcolour(9,'#B2B896' )
coot.set_map_hexcolour(10, '#CF9EB5')
coot.set_map_hexcolour(11,'#CF9EB5' )

# shader preferences
coot.set_effects_shader_brightness(1.43)
coot.set_effects_shader_gamma(0.48)
coot.set_focus_blur_z_depth(0.15)
coot.set_focus_blur_strength(1.0)

coot.graphics_draw()

{% endhighlight %}

And here is Lucrezia's script for the background colours

{% highlight python %}

import coot
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk
import coot_gui_api
import coot_utils


main_menubar = coot_gui_api.main_menubar()
menu_label = "Background"
menu = Gtk.Menu()
menuitem = Gtk.MenuItem(menu_label)
menuitem.set_submenu(menu)
main_menubar.append(menuitem)
menuitem.show()


def add_simple_coot_menu_menuitem(menu, menu_item_label, activate_function):
    submenu = Gtk.Menu()
    sub_menuitem = Gtk.MenuItem(menu_item_label)
    menu.append(sub_menuitem)
    sub_menuitem.show()
    sub_menuitem.connect("activate", activate_function)

add_simple_coot_menu_menuitem(menu, 'Red',     lambda widget: coot.set_background_colour(1,0,0))
add_simple_coot_menu_menuitem(menu, 'Blue',    lambda widget: coot.set_background_colour(0,0,1))
add_simple_coot_menu_menuitem(menu, 'Aqua',    lambda widget: coot.set_background_colour(0,1,1))
add_simple_coot_menu_menuitem(menu, 'White',   lambda widget: coot.set_background_colour(1,1,1))
add_simple_coot_menu_menuitem(menu, 'Fuchsia', lambda widget: coot.set_background_colour(1,0,1))
add_simple_coot_menu_menuitem(menu, 'Green',   lambda widget: coot.set_background_colour(0,1,0))
add_simple_coot_menu_menuitem(menu, 'Yellow',  lambda widget: coot.set_background_colour(1,1,0))

red_submenu = Gtk.Menu()
red_menuitem_sub = Gtk.MenuItem(label="Other Reds")
red_menuitem_sub.set_submenu(red_submenu)
menu.append(red_menuitem_sub)
red_menuitem_sub.show()

green_submenu = Gtk.Menu()
green_menuitem_sub = Gtk.MenuItem(label="Other Greens")
green_menuitem_sub.set_submenu(green_submenu)
menu.append(green_menuitem_sub)
green_menuitem_sub.show()

blue_submenu = Gtk.Menu()
blue_menuitem_sub = Gtk.MenuItem(label="Other Blues")
blue_menuitem_sub.set_submenu(blue_submenu)
menu.append(blue_menuitem_sub)
blue_menuitem_sub.show()

grey_submenu = Gtk.Menu()
grey_menuitem_sub = Gtk.MenuItem(label="Other Greys")
grey_menuitem_sub.set_submenu(grey_submenu)
menu.append(grey_menuitem_sub)
grey_menuitem_sub.show()

aqua_submenu = Gtk.Menu()
aqua_menuitem_sub = Gtk.MenuItem(label="Other Aquas")
aqua_menuitem_sub.set_submenu(aqua_submenu)
menu.append(aqua_menuitem_sub)
aqua_menuitem_sub.show()

yellow_submenu = Gtk.Menu()
yellow_menuitem_sub = Gtk.MenuItem(label="Other Yellows")
yellow_menuitem_sub.set_submenu(yellow_submenu)
menu.append(yellow_menuitem_sub)
yellow_menuitem_sub.show()

fuchsia_submenu = Gtk.Menu()
fuchsia_menuitem_sub = Gtk.MenuItem(label="Other Fuchsias")
fuchsia_menuitem_sub.set_submenu(fuchsia_submenu)
menu.append(fuchsia_menuitem_sub)
fuchsia_menuitem_sub.show()

for i in range(1, 11):
   f = i/10
   label = "Background colour Red " + str(f)
   add_simple_coot_menu_menuitem(red_submenu, label, lambda widget,f1=f: coot.set_background_colour(f1,0,0))

for i in range(1, 21):
   f = i/20
   label = "Background colour Blue " + str(f)
   add_simple_coot_menu_menuitem(blue_submenu, label, lambda widget,f1=f: coot.set_background_colour(0,0,f1))

for i in range(1, 21):
   f = i/20
   if i == 0:
      f = 0.03
   label = "Background colour Grey " + str(f)
   add_simple_coot_menu_menuitem(grey_submenu, label, lambda widget,f1=f: coot.set_background_colour(f1,f1,f1))
   
for i in range(1, 11):
      f = i/10
      label = "Background colour Aqua " + str(f)
      add_simple_coot_menu_menuitem(aqua_submenu, label, lambda widget,f1=f: coot.set_background_colour(0,f1,f1))

for i in range(1, 11):
      f = i/10
      label = "Background colour Fuchsia " + str(f)
      add_simple_coot_menu_menuitem(fuchsia_submenu, label, lambda widget,f1=f: coot.set_background_colour(f1,0,f1))

for i in range(1, 11):
      f = i/10
      label = "Background colour Green " + str(f)
      add_simple_coot_menu_menuitem(green_submenu, label, lambda widget,f1=f: coot.set_background_colour(0,f1,0))

for i in range(1, 11):
      f = i/10
      label = "Background colour Yellow " + str(f)
      add_simple_coot_menu_menuitem(yellow_submenu, label, lambda widget,f1=f: coot.set_background_colour(f1,f1,0))


{% endhighlight %}
