# povray.py
#
# Copyright 2005, 2006 by Paul Emsley, The University of York
# Copyright 2005, 2006 by Bernhard Lohkamp
# Copyright 2007, 2008 by Bernhard Lohkamp, The University of York
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


# BL python script for rendering povray images in coot

# define some names:
coot_povray_file_name = "coot.pov"
coot_povray_png_file_name = "coot-povray.png"
coot_povray_log_file_name = "coot-povray.log"


# define povray exe name here
# think is differnt under non windows!
import os
if (os.name == 'nt'):
	import platform
	bits, junk = platform.architecture()
	if (bits == '64bit'):
		povray_command_name = "pvengine64"
	else:
		povray_command_name = "pvengine"
else:
	povray_command_name = "povray"

# args not including the output filename
#def raster3d.povray_args():
#    return " +FN16 +A"
# BL says: dont know how usefull this function is/will be....
# BL says:: this is actually defined in raster3d.py

# Run provray using current displayed image and write .pov file to
# default filename
#
def povray_image():
    import os
    import webbrowser
    import sys

    povray(coot_povray_file_name)
    print("calling povray with args: ", raster3d.povray_args())
    extra_args = "-UV +H600 +W600"
    if (os.name == 'nt'):
        args = " /EXIT /RENDER "
    else:
        args = " "
    args = args + coot_povray_file_name + " " + raster3d.povray_args() + " " + extra_args
    # BL says: dunno what povray exe is called on other systems,
    # just assume is same for now
    povray_exe = coot_utils.find_exe(povray_command_name, "PATH")
    if (povray_exe):
      povray_call = povray_exe + args + " +o" + coot_povray_png_file_name
      print("BL DEBUG:: povray_call", povray_call)
      major, minor, micro, releaselevel, serial = sys.version_info
      if (major >= 2 and minor >=4):
          # new style
          import subprocess
          status = subprocess.call(povray_call, shell=True)
          if status:
              # something went wrong with raster3d
              # maybe same for system call?!?
              print("BL WARNING:: some error in povray")
              return
      else:
          os.system(povray_call)
      print("INFO:: displaying...")
      try:
         webbrowser.open(coot_povray_png_file_name,1,1)
      except OSError:
         print("BL WARNING:: We can't find rendered file ",coot_povray_png_file_name)

#povray_image()
