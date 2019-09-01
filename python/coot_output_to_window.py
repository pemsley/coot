

from __future__ import print_function
import sys
import os
import pygtk
pygtk.require('2.0')

import time
import gtk
import gobject
import subprocess
import select
import pango

if False: # debugging
    print("gtk version", gtk.gtk_version)
    print("pygtk version", gtk.pygtk_version)

class Application():

    def dont_press_that(self):
        w = gtk.Window()
        w.set_title("Close Coot")
        label = gtk.Label("        You really don't want to press that        ")
        vbox = gtk.VBox(False, 6)
        hbox = gtk.HBox(True, 6)
        vbox.pack_start(label, False, False, 8)
        vbox.pack_start(hbox,  False, False, 6)
        button_yes = gtk.Button(" Yeah, I do ")
        button_no = gtk.Button("     OK      ")
        hbox.pack_start(button_yes, False, False, 8)
        hbox.pack_start(button_no,  False, False, 8)
        w.add(vbox)
        button_yes.connect("clicked", self.close_everything, w)
        button_no.connect("clicked", self.close_mini_dialog, w)
        w.show_all()

    def close_mini_dialog(self, button, mini_window):
        mini_window.destroy()

    def close_everything(self, button, mini_dialog):
        mini_dialog.destroy()
        self.window.destroy()

    def window_delete(self, something, otherthing):
        self.dont_press_that()
        return True

    def __init__(self, argv):

        path = os.getenv("PATH")
        print("PATH", path)

        self.coot_log_file = "coot-output.log"
        self.window = gtk.Window()
        self.window.set_title("Coot Output")
        self.create_widgets()
        self.window.iconify()
        self.window.connect("delete_event", self.window_delete)

        coot_arg_list = argv[1:]

        if False: # debugging
            print('----')
            print(dir(pango))
            print('----')

        buffer = self.textview.get_buffer()
        self.tag_bold    = buffer.create_tag("bold")
        self.tag_orange  = buffer.create_tag("orange")
        self.tag_brown   = buffer.create_tag("brown")
        self.tag_grey    = buffer.create_tag("grey")
        self.tag_info    = buffer.create_tag("info")
        self.tag_warning = buffer.create_tag("warning")
        self.tag_debug   = buffer.create_tag("debug")
        self.tag_error   = buffer.create_tag("error")

        self.tag_bold.set_property("foreground", "blue")
        self.tag_bold.set_property("background", "black")

        self.tag_brown.set_property("foreground", "brown")

        self.tag_info.set_property("foreground", "#002200")
        self.tag_info.set_property("weight", pango.WEIGHT_BOLD)

        self.tag_warning.set_property("foreground", "brown")
        self.tag_warning.set_property("weight", pango.WEIGHT_BOLD)

        self.tag_debug.set_property("foreground", "grey")
        self.tag_debug.set_property("weight", pango.WEIGHT_LIGHT)

        self.tag_error.set_property("foreground", "#990000")
        self.tag_error.set_property("weight", pango.WEIGHT_BOLD)

        self.tag_grey.set_property("foreground", "#444444")

        self.add_timeout_check_for_new_text()

        proc_list = ['./coot-bin', '--no-state-script']
        proc_list = ['coot']
        proc_list += coot_arg_list
        
        self.process = subprocess.Popen(proc_list, stdout=subprocess.PIPE)

        self.window.show_all()

    def create_widgets(self):
        # self.vbox = gtk.VBox(spacing=6,homogeneous=False)
        self.vbox = gtk.VBox(False, 6)
        self.vbox.set_homogeneous(False)
        self.vbox.set_border_width(6)

        self.scrolled_win = gtk.ScrolledWindow()
        self.scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)

        tw = gtk.TextView()
        self.textview = tw
        self.vbox.pack_start(self.scrolled_win)
        self.scrolled_win.add(self.textview)
    
        # self.hbox_2 = gtk.HBox(spacing=10)
        # self.button_exit = gtk.Button("Exit")
        # self.hbox_2.pack_start(self.button_exit, fill=False, expand=False)
        # self.vbox.pack_start(self.hbox_2)
        
        self.window.add(self.vbox)
        self.window.set_default_size(700,300)
    
    def callback_exit(self, widget, callback_data=None):
        gtk.main_quit()

    def add_timeout_check_for_new_text(self): 
        gobject.timeout_add(100, self.check_for_new_lines_and_show_them)
        pass

    def add_text_to_buffer_using_tag(self, text, tag):
        buffer = self.textview.get_buffer()
        prev_end = buffer.get_end_iter()
        prev_end_offset = prev_end.get_offset()
        b   = buffer.get_start_iter()
        buffer.insert(prev_end, text)
        b   = buffer.get_start_iter()
        new_end = buffer.get_end_iter()
        iter_at_prev_end = buffer.get_iter_at_offset(prev_end_offset)
        buffer.apply_tag(tag, iter_at_prev_end, new_end)

    def update_text_view_with_new_text(self, raw_text):
        text = raw_text
        if raw_text[:4] == "[1m":
            text = raw_text[4:-5] # strip escaped chars and newline
            text += "\n"
            self.add_text_to_buffer_using_tag(text, self.tag_brown)
        else:
            if text[:6] == "INFO::":
                self.add_text_to_buffer_using_tag(text, self.tag_info)
            else:
                if text[:9] == "WARNING::":
                    self.add_text_to_buffer_using_tag(text, self.tag_warning)
                else:
                    if text[:7] == "DEBUG::":
                        self.add_text_to_buffer_using_tag(text, self.tag_debug)
                    else:
                        if text[:7] == "ERROR::":
                            self.add_text_to_buffer_using_tag(text, self.tag_error)
                        else:
                            if text[:13] == "Error loading":
                                self.add_text_to_buffer_using_tag(text, self.tag_error)
                            else:
                                self.add_text_to_buffer_using_tag(text, self.tag_grey)


        adj = self.scrolled_win.get_vadjustment()
        # adj.set_value(adj.get_upper()) # needs new pygtk
        v = adj.get_value()
        new_v = v + len(text)
        adj.set_value(new_v)

    def save_to_log_file(self, text):
        f = open(self.coot_log_file, "a")
        try:
            f.write(text)
            f.close()
        except IOError as e:
            # lots of error messages?
            print(e)

    def check_for_new_lines_and_show_them(self):

       # We want to exit when the process is dead.
       # We want to return True when when we have finished reading all the
       # text in the pipeline.
       # We want to send off lines for display when we see a new-line.

       if self.process.poll() == 0:
           gtk.main_quit()
           exit()

       has_input = True
       while has_input:
          # time.sleep(0.3)
          s = select.select([self.process.stdout], [], [], 0.001)
          s_out = s[0]
          if s_out:
             running_line = self.process.stdout.readline()
             self.update_text_view_with_new_text(running_line)
             self.save_to_log_file(running_line)
          else:
             has_input = False
       return True # continue timeout

if __name__ == "__main__":
    app = Application(sys.argv)

    gtk.main()

