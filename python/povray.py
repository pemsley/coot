# Copyright 2004, 2005 by The University of York
# Copyright 2005 by Bernhard Lohkamp
 
#;;; This program is free software; you can redistribute it and/or modify
#;;; it under the terms of the GNU General Public License as published by
#;;; the Free Software Foundation; either version 2 of the License, or (at
#;;; your option) any later version.
 
#;;; This program is distributed in the hope that it will be useful, but
#;;; WITHOUT ANY WARRANTY; without even the implied warranty of
#;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#;;; General Public License for more details.
 
#;;; You should have received a copy of the GNU General Public License
#;;; along with this program; if not, write to the Free Software
#;;; Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#;; direction is either 'forwards or 'backwards
#;; 
#;; start-resno is higher than stop-resno if we are building backwards
#;; 
#;; (fit-gap 0 "A"  23 26)   ; we'll build forwards
#;; (fit-gap 0 "A"  26 23)   ; we'll build backwards
#;;

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
 
