# tips-gui.py
# Copyright 2006 by Bernhard Lohkamp
# Copyright 2006 by Paul Emsley, The University of York
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

# given a number and a gtk text widget @var{text} (textbuffer!), put tip number
# @var{n} into the widget.
def show_coot_tip_from_list(n, text):
    tip = tips.tip_list()[n]
    textbuffer = text.get_buffer()
    textbuffer.set_text(tip)

# increment the tip number when the user sees a tip
def increment_coot_tip_number():
    global coot_tip_number
    if (coot_tip_number == len(tips.tip_list())-1):
       coot_tip_number = 0
    else:
       coot_tip_number += 1

# decrement the tip number when the user sees a tip
def decrease_coot_tip_number():
    global coot_tip_number
    if (coot_tip_number == 0):
       coot_tip_number = len(tips.tip_list())-1
    else:
       coot_tip_number -= 1

# run the tips gui.
def tips_gui():

    import pygtk, gtk, pango 
    import random

    global do_coot_tips_flag, coot_tip_number
    # print ":::::::::::::::::::: do_coot_tips_flag:", do_coot_tips_flag

    def delete_event(*args):
        window.destroy()
        return False

#    def destroy(*args):
#        gtk.main_quit()

    def do_next_tip(widget,data=None):
        global coot_tip_number
        increment_coot_tip_number()
        show_coot_tip_from_list(coot_tip_number,text)

    def do_prev_tip(widget,data=None):
        global coot_tip_number
        decrease_coot_tip_number()
        show_coot_tip_from_list(coot_tip_number,text)

    if do_coot_tips_flag:
       window = gtk.Window(gtk.WINDOW_TOPLEVEL)
       text = gtk.TextView()
       text.modify_base(gtk.STATE_NORMAL, gtk.gdk.color_parse("#bfe6bf"))
       text.set_editable(False)
       textbuffer = text.get_buffer()
       scrolled_win = gtk.ScrolledWindow()
       cancel_button = gtk.Button("  Close  ")
       next_tip_button = gtk.Button("  Next Tip  ")
       prev_tip_button = gtk.Button("  Previous Tip  ")
       hbox = gtk.HBox(False, 0)
       vbox = gtk.VBox(False, 0)
       window.set_default_size(490,120)
       window.add(vbox)
       window.set_title("Tips")
       vbox.set_border_width(10)
       scrolled_win.add(text)
       vbox.add(scrolled_win)
       vbox.add(hbox)
       hbox.add(next_tip_button)
       hbox.add(prev_tip_button)
       hbox.add(cancel_button)
       scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
    
       coot_tip_number = random.randint(0,len(tips.tip_list())-1)
       show_coot_tip_from_list(coot_tip_number,text)
    
       cancel_button.connect_object("clicked",delete_event, window)
       next_tip_button.connect("clicked", do_next_tip)
       prev_tip_button.connect("clicked", do_prev_tip)
       window.connect("delete_event", delete_event)
#       window.connect("destroy", delete_event)
    
       window.show_all()
#       gtk.main()

