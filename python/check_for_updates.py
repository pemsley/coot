import numbers


# Is this true? Dont understand this
def get_stable_release_from_server_string(stri):
    return stri

# needs testing!
def get_stable_release_from_coot_version():
    s = coot_version()
    ls = s.split(" ")
    return ls[0]

# return True or False (None on error, maybe should return some Error!?)
# although if no return, then var becomes None anyway.. FIXME).
def new_version_on_server(stri, is_pre_release):

    from types import StringType
    if is_pre_release:

        # pre-releases are handle with svn revision numbers (as ints).
        #
        server_rev = coot_utils.get_revision_from_string(stri)
        if server_rev:
            return server_rev > svn_revision()
        else:
            return None

    else:
        # For a stable release: 
        #
        # The Coot release number (e.g. "0.6.4") is a string that can
        # be compared lexographically.
        #
        server_version = get_stable_release_from_server_string(stri)
        this_build_version = get_stable_release_from_coot_version()
        if (server_version is StringType):
            return server_version > this_build_version
    return None

# show the dialog
# 
def notify_of_new_version(stri, use_curl=False):
    #print "notify_of_new_version given string:", stri

    # Maybe protect from being called again (download is running)
    download_binary_dialog(coot_split_version_string(stri), use_curl)

# version_string is something like: "coot-0.6-pre-1-revision-2060"
def download_binary_dialog(version_string, use_curl=False):

    import os
    import operator
    global pending_install_in_place
    global stop_download
    global file_name_for_progress_bar
    stop_download = False
    file_name_for_progress_bar = False
    if (pending_install_in_place ==  "cancelled"):
        pending_install_in_place = False
    
    #print "BL DEBUG:: in download binary with version", version_string
    #print "BL DEBUG:: download_binary_dialog with", version_string, use_curl

    def do_download_dialog(use_curl=False):
        #
        # print "running download-binary-dialog with version-string arg:", version_string
        #
        def update_progress_bar(progress_bar):
            if file_name_for_progress_bar:
                curl_info = curl_progress_info(file_name_for_progress_bar)
                if curl_info:
                    v1 = curl_info['content-length-download']
                    v2 = curl_info['size-download']                        
                    if isinstance(v1, numbers.Number):
                        if isinstance(v2, numbers.Number):
                            f = v2 / v1
                            #print "count %s, active_count %s, f: %s" %(count, active_count, f)
                            progress_bar.set_fraction(f)

        def set_progress_bar_full(progress_bar):
            progress_bar.set_fraction(1)

        def set_file_name_func(file_name):
            global file_name_for_progress_bar
            file_name_for_progress_bar = file_name

        def pending_install_in_place_func(val):
            global pending_install_in_place
            pending_install_in_place = val

        def delete_event(*args):
            # only close window (BL says:: maybe should be hide and a global window
            # in case we want to go back to it!?)
            window.destroy()
            return False

        def cancel_event(*args):
            # stop the download process(es) if
            # they were running (they get set to
            # False when they are not running)
            # 
            global pending_install_in_place
            if file_name_for_progress_bar:
                stop_curl_download(file_name_for_progress_bar)
            pending_install_in_place = "cancelled"  # not fail!
            return False

        def ok_button_event(*args):
            global pending_install_in_place
            import operator

            # only start when ok button is pressed
            gobject.timeout_add(500, idle_func) # update every 500ms is god enough?!
            if not isinstance(revision, numbers.Number):
                info_dialog("Failed to communicate with server")
            else:
                # BL says:: check if we have a download available already?! (from before)
                if (version_string in get_latest_pending_version()):
                    # installer already here (win only?!)
                    pending_install_in_place = True
                else:
                    # download
                    ok_button.set_sensitive(False)
                    #align.show()
                    progress_bar.show()
                    cancel_button.show()
                    def threaded_func():
                        global pending_install_in_place
                        ret = run_download_binary_curl(revision, version_string,
                                                       pending_install_in_place_func,
                                                       set_file_name_func,
                                                       progress_bar,
                                                       use_curl)

                        if ((not ret) and
                            (not pending_install_in_place == "cancelled")):
                            print("run_download_binary_curl failed")
                            pending_install_in_place = "fail"

                    coot_gui.run_python_thread(threaded_func, [])

        #  and a timeout checking the progress of the download:
        def idle_func():
            global pending_install_in_place
            if (pending_install_in_place == "fail"):
                window.destroy()
                info_dialog("Failure to download and install binary")
                return False  # stop idle func
            if (pending_install_in_place == "cancelled"):
                window.destroy()
                return False  # stop idle func when download cancelled
            if (pending_install_in_place == "full"):  # dont go yet!? Why?
                set_progress_bar_full(progress_bar)
            if pending_install_in_place:
                window.destroy()
                restart_dialog()
                return False  # stop idle func
            else:
                update_progress_bar(progress_bar)
                return True   # continue running

        s = "   New revision available " + \
            "for this binary type:   \n" + \
            coot_sys_build_type()     + \
            "\n"                      + \
            "\n"                      + \
            version_string
        revision = coot_utils.get_revision_from_string(version_string)

        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        dialog_name = "Download "
        if coot_utils.is_windows():
            dialog_name += "installer"
        else:
            dialog_name += "binary"
        main_vbox = gtk.VBox(False, 6)
        cancel_button = gtk.Button("  Cancel  ")
        close_button = gtk.Button("  Close  ")
        ok_button = gtk.Button("  Download and Pre-install ")
        buttons_hbox = gtk.HBox(False, 6)
        progress_bar = gtk.ProgressBar()
        h_sep = gtk.HSeparator()
        info_string = gtk.Label(s)

        # needed?
        #align = gtk.Alignment(0.5, 0.5, 0, 0)
        #align.add(pbar)

        window.set_title(dialog_name)
        buttons_hbox.pack_start(ok_button, True, False, 6)
        buttons_hbox.pack_start(cancel_button, True, False, 6)
        buttons_hbox.pack_start(close_button, True, False, 6)

        main_vbox.pack_start(info_string,  True, False, 6)
        main_vbox.pack_start(h_sep,        True, False, 6)
        main_vbox.pack_start(buttons_hbox, True, False, 6)
        main_vbox.pack_start(progress_bar, True, False, 6)
        #main_vbox.pack_start(align,        True, False, 6)
        main_vbox.set_border_width(6)

        window.add(main_vbox)
        window.set_border_width(6)

        cancel_button.connect("clicked", cancel_event)
        close_button.connect("clicked", delete_event)
        ok_button.connect("clicked", ok_button_event)
        window.show_all()
        progress_bar.hide()
        #align.hide()
        cancel_button.hide()
        idle_func()  # call once for to get existing updates

    # main line
    coot_prefix = os.getenv("COOT_PREFIX")
    coot_prefix = os.path.normpath(coot_prefix) # needed for WinCoot?! FIXME
    if not coot_utils.directory_is_modifiable_qm(coot_prefix):
        if not coot_prefix:
            info_dialog("COOT_PREFIX is not set.  Download not started.")
        else:
            info_dialog("Directory " + coot_prefix + \
                        " is not modifiable.\n" + \
                        "Download/install not started.")
    else:
        do_download_dialog(use_curl)
        if (pending_install_in_place == "cancelled"):
            pending_install_in_place = False  # reset for next run
        
def restart_dialog(extra_text=""):

    def delete_event(*args):
        window.destroy()
        return False

    def ok_button_event(*args):
        coot_command = "coot"
        coot_args =[]
        if coot_utils.is_windows():
            # we need the name of the installer, should be in pending-install
            coot_command = get_latest_pending_version()  # shall be full path name
            prefix_dir = os.getenv("COOT_PREFIX")

            coot_args =["/instdir=" + prefix_dir,
                        "/startdir=" + os.getcwd(),
                        "/autoupdate"]
            
        # create a coot background subprocess
        coot_utils.run_concurrently(coot_command, coot_args)
            
        window.destroy()
        coot_real_exit(0)

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    dialog_name = "Restart Required to complete installation"
    main_vbox = gtk.VBox(False, 6)
    cancel_button = gtk.Button("  Later  ")
    ok_button = gtk.Button("  Restart Now  ")
    label = gtk.Label("  Restart required to complete install " + extra_text)
    buttons_hbox = gtk.HBox(False, 6)
    h_sep = gtk.HSeparator()

    window.set_title(dialog_name)
    window.add(main_vbox)
    window.set_border_width(6)
    main_vbox.pack_start(label,        False, False, 6)
    main_vbox.pack_start(h_sep,        False, False, 5)
    main_vbox.pack_start(buttons_hbox, False, False, 6)
    buttons_hbox.pack_start(ok_button, False, False, 6)
    buttons_hbox.pack_start(cancel_button, False, False, 6)

    cancel_button.connect("clicked", delete_event)
    ok_button.connect("clicked", ok_button_event)

    window.show_all()


# Paul's test binary?
# http://www.biop.ox.ac.uk/coot/software/binaries/pre-releases/coot-0.6-pre-1-revision-2535-binary-Linux-i386-centos-4-gtk2.tar.gz


def check_for_updates_gui(use_curl=False):

    from types import StringType
    import time
    global server_info_status

    server_info_status = False

    def handle_latest_version_server_response(txt_from_server):
        global server_info_status
        # OK, so the server said something.
        # Set the status here, so that the
        # function that looks to see whether
        # or not the server responded is
        # notified.
        #
        if ("The requested URL" and "was not found on this server" in txt_from_server):
            server_info_status = "file-not-found"
        else:
            server_info_status = txt_from_server
    
    def get_server_info_status_thread():
        url = coot_utils.make_latest_version_url()
        print("INFO:: get URL", url)
        if use_curl:
            # non pythonic
            #x=get_url_as_string(url) # FIXME to trick the firewall
            latest_version_server_response = coot_get_url_as_string(url)
        else:
            # pythonic version
            import urllib.request, urllib.parse, urllib.error
            try:
                latest_version_server_response = urllib.request.urlopen(url).read()
            except:
                info_dialog("Could not establish connection to get latest version on server")
                return
        try:
            handle_latest_version_server_response(latest_version_server_response)
        except:
            print("BL INFO:: problem getting server response from for url", url)            # for now we give file-not-found, there should be some other form
            # of error
            handle_latest_version_server_response("The requested URL was not found on this server")

    # main line
    #
    coot_gui.run_python_thread(get_server_info_status_thread, ())
    is_pre_release = coot_utils.pre_release_qm()  # BL: why?
    global count
    count = 0
    def timeout_count():
        global count
        global server_info_status
        if (count > 2000):  # try for 20 seconds, otherwise timeout.
            # fail_with_timeout
            print("final fail: server_info_status:", server_info_status)
            # maybe some info here too?!!? 
            return False  # stop running this idle function
        elif (server_info_status == "file-not-found"):
            s = "No " + \
                ("pre-release" if is_pre_release else "release") + \
                " binary for this system (" + \
                coot_sys_build_type() + \
                ") on the binary server."
            info_dialog(s)
            return False  # stop running idle function
        elif (server_info_status == ""):
            # BL says:: supposedly - havent checked
            # this happens when coot can't get past the proxy
            # 
            info_dialog("Can't communicate with server")
        elif server_info_status:
            if (new_version_on_server(server_info_status, is_pre_release)):
                notify_of_new_version(server_info_status, use_curl)
            else:
                s = "No version newer than this revision (" + \
                    str(svn_revision())                     + \
                    ")."
                info_dialog(s)
            return False # stop running idle function
        else:
            time.sleep(0.010)
            # print "server_info_status", server_info_status
            count += 1
            return True
            
    #run_with_gtk_threading(timeout_count)
    gobject.idle_add(timeout_count)

# returns latest installer version in pending-install or empty string ""
# if no installer for windows use only!
def get_latest_pending_version():
    if coot_utils.is_windows():
        import glob
        prefix_dir = os.getenv("COOT_PREFIX")
        pending_dir = os.path.join(prefix_dir, "pending-install")
        if os.path.isdir(pending_dir):
            installer_glob = os.path.join(pending_dir, "WinCoot*.exe")
            installer_filenames = glob.glob(installer_glob)
            installer_name = ""
            if installer_filenames:
                installer_name = installer_filenames[0]
            if (len(installer_filenames) > 1):
                # shall use the newest?!
                tmp_time = 0
                for filename in installer_filenames:
                    mtime = os.path.getmtime(filename)
                    if mtime > tmp_time:
                        tmp_time = mtime
                        installer_name = filename
                # now delete the other ones
                for filename in installer_filenames:
                    if (filename != installer_name):
                        os.remove(filename)
            return installer_name
    return ""

# hack to check if a pending install for WinCoot?!
# if so, show the restart dialog?!
# not sure yet how to do it otherwise (without chaning runwincoot.bat)
if coot_utils.is_windows():
    pending_version_full = get_latest_pending_version()
    if pending_version_full:
        pending_version = os.path.basename(pending_version_full)
        restart_dialog("\n\nFound %s to install\n" %pending_version + \
                       "Please \"Restart Now\"")
    
