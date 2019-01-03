

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
    
    def __init__(self, argv):

        path = os.getenv("PATH")
        print("PATH", path)

        self.coot_log_file = "coot-output.log"
        self.window = gtk.Window()
        self.window.set_title("Coot Output")
        self.create_widgets()
        self.window.iconify()

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

       # print("debug:: check_for_new_lines_and_show_them() start")

       # We want to return False when the process is dead.
       # We want to return True when when we have finished reading all the
       # text in the pipeline.
       # In the meantime, we want to send off lines for display when we
       # see a new-line.

       # check if the process is alive!

       # set this to False when the stdout for the process dies
       continue_timeout_status = True

       has_input = True
       while has_input:
          # time.sleep(0.3)
          s = select.select([self.process.stdout], [], [], 0.001)
          s_out = s[0]
          if s_out:
             running_line = self.process.stdout.readline()
             # print("sending for display: {}".format(running_line))
             self.update_text_view_with_new_text(running_line)
             self.save_to_log_file(running_line)
          else:
             has_input = False
       # print("debug:: check_for_new_lines_and_show_them() returns", continue_timeout_status)
       return continue_timeout_status

if __name__ == "__main__":
    app = Application(sys.argv)

    gtk.main()

