
# Define an environment variable for the place where jligand.jar
# resides in a canonical distribution
# 

# should be global?
to_jligand_secret_file_name = ".coot-to-jligand-8lcs"
global from_jligand_secret_file_name
from_jligand_secret_file_name = ".jligand-to-coot"

jligand_home = os.getenv("JLIGAND_HOME") if os.getenv("JLIGAND_HOME") else "."

java_command = "java"  # may have to absolute FIXME
jligand_jar = os.path.join(os.path.normpath(jligand_home), "jligand.jar")
jligand_args = ["-jar", jligand_jar, "-coot"]

def write_file_for_jligand(resname_1, resname_2):

    fin = open(to_jligand_secret_file_name, 'w')
    int_time = time.time()
    fin.write("CODE " + resname_1 + "\n")
    fin.write("CODE " + resname_2 + "\n")
    fin.write("TIME " + str(int_time) + "\n")
    fin.close()

if True:
  if coot_gui_api.main_menubar():
      
      menu = coot_gui.coot_menubar_menu("JLigand")

      def launch_jligand():
          if not os.path.isfile(jligand_jar):
              s = "jligand java jar file: " + jligand_jar + " not found"
              coot.info_dialog(s)
          else:
              java_exe = coot_utils.find_exe(java_command)
              if java_exe:
                  coot_utils.run_concurrently(java_exe, jligand_args)
              else:
                  print "BL INFO:: no java found"
          
      coot_gui.add_simple_coot_menu_menuitem(
          menu, "Launch Jligand",
          lambda func: launch_jligand()
          )


      def make_jligand_link():
          def link_em(*args):
              print "we received these clicks", args
              if (len(args) == 2):
                  click_1 = args[0]
                  click_2 = args[1]
                  if ((len(click_1) == 7)
                      and (len(click_2) ==7)):
                      resname_1 = coot.residue_name(click_1[1],
                                               click_1[2],
                                               click_1[3],
                                               click_1[4])
                      resname_2 = coot.residue_name(click_2[1],
                                               click_2[2],
                                               click_2[3],
                                               click_2[4])
                      if not (isinstance(resname_1, str) and
                              isinstance(resname_2, str)):
                          print "Bad resnames: %s and %s" %(resname_1, resname_2)
                      else:
                          write_file_for_jligand(resname_1, resname_2)
                          
          user_defined_click(2, link_em)
          
      coot_gui.add_simple_coot_menu_menuitem(
          menu, "Send Link to JLigand (click 2 monomers)",
          lambda func: make_jligand_link()
          )

def handle_read_from_jligand_file():
    global from_jligand_secret_file_name
    if (os.path.isfile(from_jligand_secret_file_name)):
        coot.read_cif_dictionary(from_jligand_secret_file_name)

global startup_mtime
startup_mtime = prodrg_import.get_file_latest_time(from_jligand_secret_file_name)

def jligand_timeout_func():
    
    global startup_mtime
    global from_jligand_secret_file_name

    import operator
    
    now_time = prodrg_import.get_file_latest_time(from_jligand_secret_file_name)
    #print "startup-mtime: %s   now-time: %s" %(startup_mtime, now_time)
    if operator.isNumberType(now_time):
        if not operator.isNumberType(startup_mtime):
            startup_mtime = now_time
            print "just set startup-mtime to (first)", startup_mtime
            handle_read_from_jligand_file()
        else:
            if (now_time > startup_mtime):
                startup_mtime = now_time
                print "just set startup-mtime to", startup_mtime
                handle_read_from_jligand_file()
                
    return True # we never expire...
                
gobject.timeout_add(500, jligand_timeout_func)
