# The file name of the phone home script, executed with python
# FIXME - need to be installed in xxx/etc (or some such).

def phone_home_cmd():
    if (os.name == 'nt'):
        home = os.getenv('COOT_HOME')
    else:
        home = os.getenv('HOME')
    file_name = "Projects/coot/update-binary/phone_home.py"
    full_name = os.path.join(home, os.path.normpath(file_name))
    return full_name

def get_revision_from_string(str):
    # e.g. str is "coot-0.6-pre-1-revision-2060" (with a newline at the
    # end too).  We want to return 2060 (a number) from here (or False).
    s, tmp = str.rsplit("\n")
    ls = s.split("-")
    try:
        return int(ls[-1])
    except:
        return False

# Is this true?
def get_stable_release_from_server_string(stri):
    return stri

# needs testing!
def get_stable_release_from_coot_version():
    s = coot_version()
    ls = s.split(" ")
    return ls[0]

# return True or False
def new_version_on_server(str, is_pre_release):

    from types import StringType
    if is_pre_release:

        # pre-releases are handle with svn revision numbers (as ints).
        #
        server_rev = get_revision_from_string(str)
        if server_rev:
            return server_rev > svn_revision
    else:
        # For a stable release: 
        #
        # The Coot release number (e.g. "0.6.4") is a string that can
        # be compared lexographically.
        #
        server_version = get_stable_release_from_server_string(str)
        this_build_version = get_stable_release_from_coot_version()
        if (server_version is StringType):
            return server_version > this_build_version

def notify_of_new_version(str):
    ls  = str.split("c") # ??? FIXME
    ls2 = ls[-1].split()    # ??? FIXME

    download_binary_dialog(ls2[0])

# version_string is something like: "coot-0.6-pre-1-revision-2060"
def download_binary_dialog(version_string):

    import os
    is_windows = False
    if (os.name == 'nt'):
        is_windows = True

    def delete_event(*args):
        window.destroy()
        return False

    def ok_button_event(*args):
        args = [download_bin_cmd, "binary", coot_sys_build_type()]
        st = coot_version
        if "-pre-" in st:
            args.append("pre-release")
        args.extend(["version", version_string])
        print "BL DEBUG:: args are", args
        popen_commmand("python", args, "tmp_coot_download_cmd.log", False)


    s = "New revision available " + \
        "for this binary type:\n" + \
        coot_sys_build_type()     + \
        "\n"                      + \
        version_string

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    dialog_name = "Download "
    if is_windows:
        dialog_name += "installer"
    else:
        dialog_name += "binary"
    main_vbox = gtk.VBox(False, 6)
    cancel_button = gtk.Button("  Cancel  ")
    ok_button = gtk.Button("  Download  ")
    buttons_hbox = gtk.HBox(False, 6)
    h_sep = gtk.HSeparator()
    info_string = gtk.Label(s)

    buttons_hbox.pack_start(ok_button, True, False, 6)
    buttons_hbox.pack_start(cancel_button, True, False, 6)

    main_vbox.pack_start(info_string,  True, False, 6)
    main_vbox.pack_start(hsep,         True, False, 6)
    main_vbox.pack_start(buttons_hbox, True, False, 6)

    window.add(main_vbox)

    cancel_button.connect("clicked", delete_event)
    ok_button.connect("clicked", ok_button_event)

    window.show_all()


def check_for_updates():

    from types import StringType
    import time
    global count
    server_info_status = 'pending'
    def update_info_thread():
        print "get updates info thread"
        # here we construct args to popen_command,
        # adding in "pre-release" if this binary is a
        # pre-release.
        # args ends up as something like:
        # ["xxx/phone_home.py", "pre-release" 
        # "binary", "Linux-1386-fedora-10-python-gtk2"
        # "command-line", "/home/xx/coot/bin/coot"]
        update_coot_log = "tmp-update-coot.log"
        args = ["binary", coot_sys_build_type(),
                "command-line", "somthing"]
        if ("-pre-" in coot_version()):
            args.append("pre-release")
        args.insert(0, phone_home_cmd())
        if (os.path.isfile(update_coot_log)):
            os.remove(update_coot_log)
        print "BL DEBUG:: args are", args
        status = popen_command("python", args, [], update_coot_log, False)
        print "BL DEBUG:: status is", status
        if (os.path.isfile(update_coot_log)):
            # OK, so the server said something
            # Set the status here, so that the
            # function that looks to see wether
            # or not the server responded is
            # notified.
            #
            fin = open(update_coot_log, 'r')
            line = fin.readline()
            print "got line", line
            server_info_status = line
            # not now for testing! FIXME
            # os.remove(update_coot_log)
        # thread function ends here
        
    run_python_thread(update_info_thread, ())
    is_pre_release = "-pre-" in coot_version()
    count = 0
    def timeout_count():
        global count
        if (count > 2000):  # try for 20 seconds, otherwise timeout.
            # fail_with_timeout
            print "final fail: server_info_status:", server_info_status
            return False  # stop running this idle function
        elif (server_info_status is StringType):
            if (new_version_on_server(server_info_status, is_pre_release)):
                notify_of_new_version(server_info_status)
            else:
                s = "No version newer than this revision (" + \
                    str(svn_revision())                     + \
                    ")."
                info_dialog(s)
            return False  # stop running idle function
        else:
            time.sleep(0.010)
            # print "server_info_status", server_info_status
            count += 1
            return True
            
    gobject.idle_add(timeout_count)
            
    


menu = coot_menubar_menu("Updates")
add_simple_coot_menu_menuitem(menu,
                              "Check for updates...",
                              lambda func: check_for_updates())
