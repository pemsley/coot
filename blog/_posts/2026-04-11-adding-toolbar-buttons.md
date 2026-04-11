---
layout: post
title:  "Adding Toolbar Buttons to Coot 1"
date: Sat 11 Apr 2026 06:43:19 BST
---

# Adding a Button to the Coot Toolbar with Python

Coot's GTK4 interface exposes the main toolbar to Python scripts, making it straightforward to add custom buttons that trigger your own functions. This post walks through the pattern, drawn from working examples in the Coot source tree.

## The Minimal Pattern

Required imports:


```python
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk
import coot
import coot_gui_api

def my_callback(button):
    coot.single_model_view_next_model_number(imol)

button = Gtk.Button(label="Step")
button.connect("clicked", my_callback)
coot_gui_api.main_toolbar().append(button)
```

That's the basics. `coot_gui_api.main_toolbar()` returns the main toolbar widget, and `.append()` adds your button to the end of it. The callback receives the button widget as its argument — the standard GTK4 signal signature.

## A Complete Example: Trajectory Animation

The trajectory viewer script shows how to combine a regular button with a toggle button and a GLib timer to build a small animation controller:

```python
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib
import coot
import coot_gui_api

imol_traj = coot.read_pdb("multi-model.pdb")
coot.single_model_view_model_number(imol_traj, 1)

animation_timeout_id = None

def step_callback(button):
    coot.single_model_view_next_model_number(imol_traj)

def animate_timeout():
    """Called every 100ms during animation"""
    coot.single_model_view_next_model_number(imol_traj)
    return True  # Return True to keep the timeout running

def animate_toggled(button):
    global animation_timeout_id
    if button.get_active():
        animation_timeout_id = GLib.timeout_add(100, animate_timeout)
    else:
        if animation_timeout_id is not None:
            GLib.source_remove(animation_timeout_id)
            animation_timeout_id = None

step_button = Gtk.Button(label="Step")
step_button.connect("clicked", step_callback)
coot_gui_api.main_toolbar().append(step_button)

animate_button = Gtk.ToggleButton(label="Animate")
animate_button.connect("toggled", animate_toggled)
coot_gui_api.main_toolbar().append(animate_button)
```

Two things to note here. First, `Gtk.ToggleButton` works exactly like `Gtk.Button` for toolbar purposes — you create it, connect a signal, and append it. The difference is the signal name: `"toggled"` instead of `"clicked"`, and you query `button.get_active()` to check its state. Second, `GLib.timeout_add()` is the way to do timed callbacks in GTK4 — return `True` from the timeout function to keep it firing, or `False` to stop.

## Another Example: Alt-Conf Switcher

The alt-conf switcher combines a toolbar button with a keyboard binding. The button clears non-drawn bonds, while a key binding toggles between showing alt conf A and B:

```python
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk
import coot
import coot_gui_api

alt_conf = "A"

def switch_alt_conf():
    global alt_conf
    aa = coot.active_residue_py()
    imol = aa[0]
    cid = "//*/*/*:" + alt_conf
    coot.set_new_non_drawn_bonds(imol, cid)
    alt_conf = "B" if alt_conf == "A" else "A"

def clear_callback(button):
    active_atom = coot.active_residue_py()
    if not active_atom:
        coot.info_dialog("No Residue Here")
    else:
        coot.clear_non_drawn_bonds(active_atom[0])

coot.add_key_binding_gtk4_py(103, 0, switch_alt_conf, "Switch Alt Conf")

button = Gtk.Button(label="ClearAltConfMode")
button.connect("clicked", clear_callback)
coot_gui_api.main_toolbar().append(button)
```

The toolbar button and the key binding are independent — one doesn't need the other — but together they make a useful workflow tool. The key binding (G) toggles which alt conf is hidden; the button resets the display to show everything.

## Summary

The key points:

- **`coot_gui_api.main_toolbar()`** gives you the toolbar widget
- **`.append(widget)`** adds your button at the end
- Use **`Gtk.Button`** for simple click actions, **`Gtk.ToggleButton`** for on/off state
- Callbacks receive the button widget as their first argument
- Use **`GLib.timeout_add()`** for timed/animated operations
- These scripts run at startup when placed in the appropriate location, or can be loaded at runtime via the scripting console

