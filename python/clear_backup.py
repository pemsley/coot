import numbers
#    <one line to give the program's name and a brief idea of what it does.>
#    Copyright (C) <year>  <name of author>
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

global clear_out_backup_run_n_days
global clear_out_backup_old_days

clear_out_backup_run_n_days = 7     # run every 7 days
clear_out_backup_old_days   = 7     # Files older than 7 days are condisdered
                                    # for deletion

#
def delete_coot_backup_files(action_type):

    import glob
    import time
    import operator
    global clear_out_backup_old_days

    # check global COOT_BACKUP_DIR (if exists)
    backup_env = os.getenv('COOT_BACKUP_DIR')
    dirs = []
    dir_str = ""
    files = []
    if (backup_env):
        global_backup_patt = os.path.normpath(backup_env + "/*.pdb*")
        global_files = glob.glob(global_backup_patt)
        files = global_files
        dir_str = backup_env
        dirs = [backup_env]
    # check local dir too if exists
    if (os.path.isdir("coot-backup")):
        local_backup_patt  = os.path.normpath("./coot-backup/*.pdb*")
        local_files  = glob.glob(local_backup_patt)
        files += local_files
        if (dir_str):
            dir_str += " and ./coot-backup"
        else:
            dir_str = "./coot-backup"
        dirs.append("./coot-backup")
    now = int(time.time())
    if (not isinstance(clear_out_backup_old_days, numbers.Number)):
        clear_out_backup_old_days = 7
    n_days = clear_out_backup_old_days
    last_week = (now - (n_days * 24 * 60 * 60))   # clear out more than 7 days

    def old_files_list(files):
        old_files = []
        for file in files:
            this_mtime = os.path.getmtime(file)
            #print "comparing ", last_week, this_mtime
            if (this_mtime < last_week):
                old_files.append(file)
        return old_files

    # main
    olds = old_files_list(files)
    if (action_type == "count"):
        total_size = sum(map(os.path.getsize, olds))
        mb_size = total_size / (1024. * 1024)
        #print "There are %s old files  (%s bytes) (%.1fMb) in %s" %(len(olds),
        #                                                          total_size,
        #                                                          mb_size,
        #                                                          dir_str)
        # obs: here mb_size is not a string!!!
        return [len(olds), mb_size]

    if (action_type == 'delete'):
        for file in files:
            file_name, dir_name = os.path.split(file)
            print("Deleting old file %s from %s" %(file_name, dir_name))
            os.remove(file)
        # now create a last-backup file with a time stamp:
        for directory in dirs:
            last_cleaned_file = os.path.normpath(os.path.join(directory, "last-cleaned"))
            if (os.path.isdir(directory)):
                fin = open(last_cleaned_file, 'w')
                fin.write(str(int(time.time())))
                fin.close()
        #return something?
                                                                  

# Make a GUI
#
# return True or False depending on if the GUI dialog was shown (it isn't
# shown if there are no files to delete).
#
def clear_backup_gui():

    file_stats = delete_coot_backup_files("count")
            
    # more than 1 file to possibly delete?
    if (file_stats[0] == 0):
        return False   # didn't run
    else:
        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        frame = gtk.Frame("Old Backups")
        vbox = gtk.VBox(False, 3)
        hbox = gtk.HBox(False, 10)
        ok_button = gtk.Button(label=" Clear up ")
        cancel_button = gtk.Button(label=" Stay messy ")
        h_sep = gtk.HSeparator()
        label_str = "  There are " + str(file_stats[0]) + \
                    " old backup files (%.1fMb) \n" %file_stats[1] + \
                    "   Delete Them?"
        label = gtk.Label(label=label_str)
        ok_button.connect("clicked", lambda w: list(map(eval, ["delete_coot_backup_files('delete')", "coot_real_exit(0)"])))
        cancel_button.connect("clicked", lambda w: coot.coot_real_exit(0))

        ok_text = " Consider yourself patted on the back! "
        cancel_text = "A less pejorative label here might be \"Keep\" or \"Cancel\" " + \
                      "but seeing as (for the moment) I like my intestines where they are " + \
                      "and not used as hosiery fastenings for Systems Adminstrators then " + \
                      "we get this rather nannying label..."
        if gtk.pygtk_version >= (2,12):
            ok_button.set_tooltip_text(ok_text)
            cancel_button.set_tooltip_text(cancel_text)
        else:
            coot_tooltips.set_tip(ok_button, ok_text)
            coot_tooltips.set_tip(cancel_button, cancel_text)

        window.add(frame)
        frame.set_border_width(6)
        frame.add(vbox)
        vbox.append(label)
        vbox.pack_start(h_sep)
        vbox.pack_start(hbox)
        hbox.set_border_width(10)
        vbox.set_border_width(6)
        hbox.append(ok_button)
        hbox.append(cancel_button)

        ok_button.set_flags(gtk.CAN_DEFAULT)
        cancel_button.set_flags(gtk.CAN_DEFAULT)
        ok_button.grab_default()

        window.show_all()
        return True  #it ran


# return a status, True or False, did the gui run?
#
# Note that clear-backup-gui returns either true or False too.
#
# If this function returns False, then coot.coot_real_exit() just exits with
# coot.coot_real_exit().  Otherwise we wait for the GUI.
#
def clear_backups_maybe():

    import time
    import operator
    global clear_out_backup_run_n_days
    
    now = int(time.time())
    if (not isinstance(clear_out_backup_run_n_days, numbers.Number)):
        clear_out_backup_run_n_days = 7
    n_days = clear_out_backup_run_n_days
    last_week = (now - (n_days * 24 * 60 * 60))   # clear out every 7 days

    # check global COOT_BACKUP_DIR if exists
    backup_env = os.getenv('COOT_BACKUP_DIR')
    dirs = []
    if (backup_env):
        dirs = [backup_env]
    # check local dir too if exists
    if (os.path.isdir("coot-backup")):
        dirs.append(os.path.normpath("./coot-backup"))

    all_last_cleaned_files = []
    for directory in dirs:
        last_cleaned_file = os.path.join(directory, "last-cleaned")
        if (os.path.isfile(last_cleaned_file)):
            all_last_cleaned_files.append(last_cleaned_file)
            
    if (len(all_last_cleaned_files) == 0):
        ret = clear_backup_gui()
        return ret
    else:
        for last_cleaned_file in all_last_cleaned_files:
            fin = open(last_cleaned_file, 'r')
            val = fin.read()
            fin.close()
            try:
                ival = int(val)
                if (ival < last_week):
                    ret = clear_backup_gui()
                    return ret
                else:
                    print("INFO:: backup clearout done %s days ago" %((now - val) * 24 * 60 * 60))
                    return False
            except:
                return False
            

    
