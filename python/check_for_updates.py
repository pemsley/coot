
def get_revision_from_string(stri):
    # e.g. str is "coot-0.6-pre-1-revision-2060" (with a newline at the
    # end too).  We want to return 2060 (a number) from here (or False).
    lss = stri.rsplit("\n")
    ls = lss[0].split("-")
    try:
        return int(ls[-1])
    except:
        return False

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
        server_rev = get_revision_from_string(stri)
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

# First generate a string with no trailing newline.
# 
def notify_of_new_version(stri):
    #print "notify_of_new_version given string:", stri

    import os

    if not is_windows():
        ls2 = stri[stri.find("c"):-1]  # for 'coot' start
    else:
        ls2 = stri[stri.find("W"):-1]  # for 'WinCoot' start

    # Maybe protect from being called again (download is running)
    download_binary_dialog(ls2)

global pending_install_in_place
pending_install_in_place = False

# version_string is something like: "coot-0.6-pre-1-revision-2060"
def download_binary_dialog(version_string):

    import os
    global pending_install_in_place
    global stop_download
    stop_download = False

    #rint "BL DEBUG:: in download binary with version", version_string
    #rint "BL DEBUG:: pending_install_in_place", pending_install_in_place

    def delete_event(*args):
        # only close window (BL says:: maybe should be hide and a global window
        # in case we want to go back to it!?)
        window.destroy()
        return False

    def cancel_event(*args):
        # This cancels download and closes window
        # 
        global stop_download
        global pending_install_in_place
        stop_download = True
        pending_install_in_place = "stop"
        window.destroy()
        return False

    def ok_button_event(*args):
        global pending_install_in_place
        # only start when ok button is pressed
        gobject.timeout_add(100, idle_func) # update every 100ms is god enough?!
        # BL says:: check if we have a download available already?! (from before)
        if (version_string in get_latest_pending_version()):
            # installer already here (win only)
            pending_install_in_place = True
        else:
            # download
            ok_button.set_sensitive(False)
            align.show()
            cancel_button.show()
            def threaded_func():
                global stop_download
                ret = run_download_binary_curl(revision, version_string)
                if (not ret) and (not stop_download):
                    print "run_download_binary_curl failed"
                    pending_install_in_place = "fail"

            run_python_thread(threaded_func, [])


    # get the binary, the action that happens when the download button is pressed.
    # This is using python/urllib by default, set use_curl to True to use curl
    #
    def run_download_binary_curl(revision, version_string, use_curl = False):
        print "::::: run_download_binary_curl.... with revision %s with version_string %s" \
              %(revision, version_string)
        prefix = os.getenv("COOT_PREFIX")
        prefix = os.path.normpath(prefix)  # FIXME do we need this?
        print "::::: run-download-binary-curl.... prefix is", prefix
        if not prefix:  # do we need to check if prefix is string?
            print "OOps! Can't find COOT_PREFIX"
        else:
            pre_release_flag = "-pre" in coot_version()
            ys = "www.ysbl.york.ac.uk/~emsley/software/binaries/"
            binary_type = coot_sys_build_type()
            if (binary_type == "Linux-i386-fedora-3") or \
               (binary_type == "Linux-i386-fedora-3-python") or \
               (binary_type == "Linux-i386-fedora-8-python-gtk2") or \
               (binary_type == "Linux-i386-fedora-8-gtk2") or \
               (binary_type == "Linux-i386-fedora-10-python-gtk2") or \
               (binary_type == "Linux-i386-fedora-10-gtk2") or \
               (binary_type == "Linux-i686-ubuntu-8.04.3") or \
               (binary_type == "Linux-i686-ubuntu-8.04.3-python"):
                host_dir = ys
            else:
                host_dir = "www.biop.ox.ac.uk/coot/software/binaries/"
            
            if is_windows():
                host_dir = "www.ysbl.york.ac.uk/~lohkamp/software/binaries/"
            
            tar_file_name = version_string
            if is_windows():
                tar_file_name += ".exe"
            else:
                tar_file_name += "-binary-" + binary_type + ".tar.gz"

            release_dir = "releases/"
            if ("ysbl.york.ac.uk" in host_dir):  # includes windows
                if pre_release_flag:
                    release_dir = "nightlies/pre-release/"
                else:
                    release_dir = "stable/"
            else:
                if pre_release_flag:
                    release_dir = "pre-releases/"
                else:
                    release_dir = "releases/"

            url = "http://" + host_dir + release_dir + tar_file_name
            md5_url = url + ".md5sum"
            md5_tar_file_name = tar_file_name +".md5sum"

            print "md5sum url for curl:", md5_url
            print "url for curl:", url

            if use_curl:
                # this is for curl
                curl_exe = find_exe("curl", os.path.join(prefix, "bin"), "PATH")

                popen_command(curl_exe,
                              [md5_url, "-o", md5_tar_file_name],
                              [], "tmp-coot-get-md5-url.log", False)
                popen_command(curl_exe,
                              [url, "-o", tar_file_name],
                              [], "tmp-coot-get-url.log", False)
            else:
                # this is for pythonic urllib
                import urllib
                def progress_function(count, blockSize, totalSize):
                    global stop_download
                    percent = int(count*blockSize*100/totalSize)
                    sys.stdout.write("\rDownloading " + tar_file_name + "...  %d%%" % percent)  # FIXME (format!!)
                    sys.stdout.flush()
                    pbar.set_text("Downloading %s %%" %percent)
                    pbar.set_fraction(percent/100.)
                    if stop_download:
                        # Brute force exit of thread!
                        print "\n\nBL INFO:: stopping download"
                        sys.exit()
                
                try:
                    md5_url_local_file_name, md5_url_info =  urllib.urlretrieve(md5_url, md5_tar_file_name)
                except:
                    print "BL ERROR:: could not download", md5_url
                    return False
                try:
                    print "\n"
                    url_local_file_name, url_info =  urllib.urlretrieve(url, tar_file_name, progress_function)
                    print "\n"
                except:
                    print "BL ERROR:: could not download", url
                    return False

            if not os.path.isfile(tar_file_name):
                print "Ooops: %s does not exist after attempted download" %tar_file_name
                return False
            else:
                if not os.path.isfile(md5_tar_file_name):
                    print "Ooops: %s does not exist after attempted download" %md5_tar_file_name
                    return False
                else:
                    if not match_md5sums(tar_file_name, md5_tar_file_name):
                        return False
                    else:
                        success = install_coot_tar_file(tar_file_name)
                        if success:
                            global pending_install_in_place
                            pending_install_in_place = True
                            return True
                        return False

    #
    # print "running download-binary-dialog with version-string arg:", version_string
    #

    s = "   New revision available " + \
        "for this binary type:   \n" + \
        coot_sys_build_type()     + \
        "\n"                      + \
        "\n"                      + \
        version_string
    revision = get_revision_from_string(version_string)
 
    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    dialog_name = "Download "
    if is_windows():
        dialog_name += "installer"
    else:
        dialog_name += "binary"
    main_vbox = gtk.VBox(False, 6)
    cancel_button = gtk.Button("  Cancel  ")
    close_button = gtk.Button("  Close  ")
    ok_button = gtk.Button("  Download and Pre-install ")
    buttons_hbox = gtk.HBox(False, 6)
    h_sep = gtk.HSeparator()
    info_string = gtk.Label(s)

    align = gtk.Alignment(0.5, 0.5, 0, 0)
    pbar = gtk.ProgressBar()
    align.add(pbar)

    window.set_title(dialog_name)
    buttons_hbox.pack_start(ok_button, True, False, 6)
    buttons_hbox.pack_start(cancel_button, True, False, 6)
    buttons_hbox.pack_start(close_button, True, False, 6)

    main_vbox.pack_start(info_string,  True, False, 6)
    main_vbox.pack_start(h_sep,        True, False, 6)
    main_vbox.pack_start(buttons_hbox, True, False, 6)
    main_vbox.pack_start(align,        True, False, 6)
    main_vbox.set_border_width(6)

    window.add(main_vbox)
    window.set_border_width(6)

    cancel_button.connect("clicked", cancel_event)
    close_button.connect("clicked", delete_event)
    ok_button.connect("clicked", ok_button_event)

    #  now the timer that waits for the binary...
    def idle_func():
        #import time  # not needed since we use timeout_add
        #time.sleep(0.01)
        global pending_install_in_place
        if (pending_install_in_place == "fail"):
            window.destroy()
            info_dialog = "Failure to download and install binary"
            return False  # stop idle func
        if (pending_install_in_place == "stop"):
            pending_install_in_place = False  # reset for next run
            window.destroy()
            return False  # stop idle func when download canceled
        if pending_install_in_place:
            window.destroy()
            restart_dialog()
            return False  # stop idle func
        else:
            return True   # continue running

    idle_func()  # call once for to get existing updates
    window.show_all()
    align.hide()
    cancel_button.hide()


# return success status as a boolean
#
def install_coot_tar_file(tar_file_name, use_tar = False):

    prefix_dir = os.getenv("COOT_PREFIX")
    prefix_dir = os.path.normpath(prefix_dir)  # FIXME do we need this?
    if not prefix_dir:
        print "OOps could not get COOT_PREFIX"
        return False
    if not directory_is_modifiable_qm(prefix_dir):
        print "OOps directory %s is not modifiable" %prefix_dir
        return False
    else:
        pending_dir = os.path.join(prefix_dir, "pending-install")
        if not os.path.isdir(pending_dir):
            os.mkdir(pending_dir)
        if not os.path.isdir(pending_dir):
            print "OOps could not create", pending_dir
            return False
        else:
            a_tar_file_name = os.path.abspath(tar_file_name)
            # with working dir !?
            current_dir = os.getcwd()
            os.chdir(pending_dir)
            print "now current dir is", os.getcwd()
            if use_tar:
                # non-pythonic
                popen_command("tar", ["xzf", a_tar_file_name], [],
                              "untar.log", False)
            else:
                # pythonic
                if not is_windows():
                    import tarfile
                    tar = tarfile.open(a_tar_file_name)
                    tar.extractall()
                    tar.close()
                else:
                    if os.path.isfile(tar_file_name):
                        # needs to be removed on WIN32 first
                        os.remove(tar_file_name)  
                    os.rename(a_tar_file_name, tar_file_name)
                    # shall remove md5sum at some point!? FIXME                
                    
            os.chdir(current_dir)
            print "final current dir is", os.getcwd()
            return True # ?

#
def directory_is_modifiable_qm(dir_in):
    ret = False
    # check existence:
    ret = os.access(dir_in, os.F_OK)
    if ret:
        # check readability
        ret = os.access(dir_in, os.R_OK)
        if ret:
            # check writability
            ret = os.access(dir_in, os.W_OK)
            if ret:
                # check executability (needed?!)
                ret = os.access(dir_in, os.X_OK)
    return ret

# return as a string, or False
def get_target_md5_string(file_name):
    if not os.path.isfile(file_name):
        return False
    else:
        fin = open(file_name, 'r')
        lines = fin.readlines()
        fin.close()
        first_line = lines[0]
        return first_line[0:first_line.find(" ")]

# return a string
def get_md5sum_string(file_name):
    if not os.path.isfile(file_name):
        return False
    else:
        import hashlib
        fin = open(file_name, 'rb')
        md5sum = hashlib.md5(fin.read()).hexdigest()
        fin.close()
        return md5sum

def match_md5sums(tar_file_name, target_md5sum_file_name):
    # necessary to check for files? done already above?
    if not os.path.isfile(tar_file_name):
        return False
    else:
        if not os.path.isfile(target_md5sum_file_name):
            print "OOps! %s does not exist" %target_md5sum_file_name
            return False
        else:
            target_md5_string = get_target_md5_string(target_md5sum_file_name)
            # remove the md5sum file (not needed any more)
            os.remove(target_md5sum_file_name)
            md5_string = get_md5sum_string(tar_file_name)
            if not target_md5_string:    # need to test if string?
                print "OOps %s is not a string" %target_md5_string
                return False
            else:
                if not md5_string:       # as above
                    print "OOps %s is not a string" %md5_string
                    return False
                else:
                    if not (target_md5_string == md5_string):
                        print "Oops: md5sums do not match %s %s.  Doing nothing" \
                              %(target_md5_string, md5_string)
                        return False
                    else:
                        return True

        
def restart_dialog(extra_text=""):

    def delete_event(*args):
        window.destroy()
        return False

    def ok_button_event(*args):
        coot_command = "coot"
        coot_args =[]
        if is_windows():
            # we need the name of the installer, should be in pending-install
            coot_command = get_latest_pending_version()  # shall be full path name
            prefix_dir = os.getenv("COOT_PREFIX")

            coot_args =["/instdir=" + prefix_dir,
                        "/startdir=" + os.getcwd(),
                        "/autoupdate"]
            
        # create a coot background subprocess
        run_concurrently(coot_command, coot_args)
            
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


def check_for_updates_gui():

    from types import StringType
    import time
    global count
    global server_info_status

    server_info_status = False

    # return a boolean
    def pre_release_qm():
        return "-pre" in coot_version()

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
    
    # here we construct args to popen_command,
    # adding in "pre-release" if this binary is a
    # pre-release.
    # args ends up as something like:
    # ["xxx/phone_home.py", "pre-release" 
    # "binary", "Linux-1386-fedora-10-python-gtk2"
    # "command-line", "/home/xx/coot/bin/coot"]
    #
    def make_latest_version_url():
        build_type = coot_sys_build_type()
        # FIXME this is only for biop files !!!
        # what about versions in York?
        host = "http://www.biop.ox.ac.uk/coot/software/binaries/"
        pre_dir = "pre-releases" if pre_release_qm else "releases"
        if is_windows():
            host = "http://www.ysbl.york.ac.uk/~lohkamp/software/binaries/"
            pre_dir = "nightlies/pre-release" if pre_release_qm else "stable"
        url = host + \
              pre_dir + \
              "/" + \
              "type-binary-" + \
              build_type + \
              "-latest.txt"
        return url

    def get_server_info_status_thread():
        url = make_latest_version_url()
        print "INFO:: get URL", url
        # non pythonic
        #latest_version_server_response = coot_get_url_as_string(url)
        # pythonic version
        import urllib
        latest_version_server_response = urllib.urlopen(url).read()
        handle_latest_version_server_response(latest_version_server_response)

    # main line
    #
    run_python_thread(get_server_info_status_thread, ())
    is_pre_release = pre_release_qm()  # BL: why?
    count = 0
    def timeout_count():
        global count
        global server_info_status
        if (count > 2000):  # try for 20 seconds, otherwise timeout.
            # fail_with_timeout
            print "final fail: server_info_status:", server_info_status
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
        elif server_info_status:
            if (new_version_on_server(server_info_status, is_pre_release)):
                notify_of_new_version(server_info_status)
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
    if is_windows():
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
            return installer_name
    return ""

# hack to check if a pending install for WinCoot?!
# if so, show the restart dialog?!
# not sure yet how to do it otherwise (without chaning runwincoot.bat)
if is_windows():
    pending_version_full = get_latest_pending_version()
    if pending_version_full:
        pending_version = os.path.basename(pending_version_full)
        restart_dialog("\n\nFound %s to install\n" %pending_version + \
                       "Please \"Restart Now\"")
    
