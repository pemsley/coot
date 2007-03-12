# BL python script for rendering povray images in coot

def povray_args():
    return " +FN16 +A"
# BL says: dont know how usefull this function is/will be....

def povray_image():
    import os, webbrowser

    coot_povray_file_name = "coot.pov"
    coot_povray_png_file_name = "coot-povray.png"

    povray(coot_povray_file_name)
    print "calling povray with args: ", povray_args()
    extra_args = "+H600 +W600"
    args = coot_povray_file_name + " " + povray_args() + " " + extra_args
    if os.name == "nt":
       render_exe = "pvengine.exe"
    else:
# BL says: dunno what povray exe is called on other systems, just assume for now
       render_exe = "pvengine"
    povray_exe = find_render_exe("POV*",render_exe)
    if (povray_exe):
      povray_call = povray_exe + " /EXIT /RENDER " + args + " +o" + coot_povray_png_file_name
      print "BL DEBUG:: povray command line", povray_call
      os.system(povray_call)
      print "calling display..."
      try:
         webbrowser.open(coot_povray_png_file_name)
      except OSError:
         print "BL WARNING:: We can't find rendered file ",coot_povray_png_file_name

#BL testing
#povray_image()
 
