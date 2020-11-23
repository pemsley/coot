
# This happens when user clicks on the "Launch JLigand" button.
# It starts a jligand and puts it in the background.
#
def launch_jligand_function():

    global jligand_jar
    global jligand_home_env
    global java_command
    
    jligand.start_jligand_listener()
    # maybe this should rather check PATH or similar!? FIXME
    if not os.path.isfile(jligand_jar):

        # Boo.  Give us a warning dialog
        #
        s = "jligand java jar file: " + jligand_jar + " not found"

        # make an extra message telling us that JLIGAND_HOME is
        # not set if it is not set.
        env_message = "Environment variable JLIGAND_HOME not set\n\n" \
                      if not jligand_home_env else ""
        info_dialog(env_message + s)

    else:
        # OK, it does exist - run it!
        #
        java_exe = coot_utils.find_exe(java_command)
        if not java_exe:
            print("BL INFO:: no java found")
        else:
            # first check if we can run it with coot, i.e. is '-version'
            # a valid command line arg
            jligand_version = ["-jar", jligand_jar, "-version"]
            cmd = java_exe + " " + \
                  coot_utils.string_append_with_spaces(jligand_version)
            res = coot_utils.shell_command_to_string(cmd)
            if (not res):
                message = "Sorry, your JLigand:\n\n " + jligand_jar + "\n\n" + \
                          "is not new enough to work with Coot!\n" + \
                          "Please download a new one!"
                info_dialog(message)
            else:
                coot_utils.run_concurrently(java_exe, jligand_args)
                # beam in a new menu to the menu bar:
                if True:
                    if coot_gui_api.main_menubar():
                        jligand_menu = coot_gui.coot_menubar_menu("JLigand")
                        coot_gui.add_simple_coot_menu_menuitem(
                            jligand_menu, "Send Link to JLigand (click 2 monomers)",
                            lambda func: click_select_residues_for_jligand()
                            )

# This happens when user clicks on the "Select Residues for JLigand"
# (or some such) button.  It expects the user to click on atoms of
# the two residues involved in the link.
#
def click_select_residues_for_jligand():

    global imol_jligand_link
    
    def link_em(*args):
        print("we received these clicks", args)
        if (len(args) == 2):
            click_1 = args[0]
            click_2 = args[1]
            print("click_1:", click_1)
            print("click_2:", click_2)
            if ((len(click_1) == 7)
                and (len(click_2) ==7)):
                resname_1 = residue_name(click_1[1],
                                         click_1[2],
                                         click_1[3],
                                         click_1[4])
                resname_2 = residue_name(click_2[1],
                                         click_2[2],
                                         click_2[3],
                                         click_2[4])
                imol_click_1 = click_1[1]
                imol_click_2 = click_2[1]
                chain_click_1 = click_1[2]
                chain_click_2 = click_2[2]
                resno_click_1 = click_1[3]
                resno_click_2 = click_2[3]
                if not (isinstance(resname_1, str) and
                        isinstance(resname_2, str)):
                    print("Bad resnames: %s and %s" %(resname_1, resname_2))
                else:
                    if not (imol_click_1 == imol_click_2):
                        msg = "Two different molecules %s and %s selected.\n" \
                              %(imol_click_1, imol_click_2) + \
                              "Make sure to select residues in the same molecule."
                        info_dialog(msg)
                        imol_jligand_link = False
                    elif (chain_click_1 == chain_click_2 and
                          resno_click_1 == resno_click_2):
                        msg = "Same residue %s %s selected.\n" \
                              %(chain_click_1, resno_click_1) + \
                              "Make sure to select different residues."
                        info_dialog(msg)
                        imol_jligand_link = False
                    else:
                        # happy path
                        imol_jligand_link = imol_click_1
                        jligand.write_file_for_jligand(jligand.click2res_spec(click_1), resname_1,
                                               jligand.click2res_spec(click_2), resname_2)

    user_defined_click(2, link_em)
