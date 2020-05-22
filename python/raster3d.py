# raster3d.py
#
# Copyright 2005, 2006 by Bernhard Lohkamp
# Copyright 2007 by Bernhard Lohkamp, The University of York
# Copyright 2005, 2006 by Paul Emsley, The University of York
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


# BL pyton script for rendering images in COOT

# args not including the output filename
def povray_args():
    return ["+FN16", "+A"]
# BL says: dont know how usefull this function is/will be....

# run raster3d
#
def render_image():
    import os
    import webbrowser
    import sys

    coot_r3d_file_name="coot.r3d"
    version = raster3d_version()
    if (os.name == 'nt' and version[0] == 2):
        coot_image_file_name = "coot.tif"
        image_format = " -tiff "
    else:
        coot_image_file_name = "coot.png"
        image_format = " -png "
    raster3d(coot_r3d_file_name)
    r3d_exe = find_exe("render", "PATH")
    if (r3d_exe):
       r3d_dir = os.path.dirname(r3d_exe)
       os.environ['R3D_LIB'] = r3d_dir + "/materials"
       r3d_call = r3d_exe + image_format + " -labels " + coot_image_file_name + " < " + coot_r3d_file_name
       print("BL DEBUG:: r3d_call is ", r3d_call)
       print("calling render...")

       major, minor, micro, releaselevel, serial = sys.version_info
       if (major >= 2 and minor >=4):
           # new style
           import subprocess
           status = subprocess.call(r3d_call, shell=True)
           if status:
               # something went wrong with raster3d
               # maybe same for system call?!?
               print("BL WARNING:: some error in raster3d")
               return
       else:
           status = os.system(r3d_call)
       print("calling display...")
       try:
         webbrowser.open(coot_image_file_name,1,1)
       except OSError:
         print("BL WARNING:: We can't find rendered file ",coot_image_file_name)

# Run either raster3d or povray
#
# @var{image_type} is either 'raster3d' or 'povray'
#
def raytrace(image_type, source_file_name, image_file_name, x_size, y_size):
   import os
   import webbrowser
   import sys

   if (image_type == "raster3d"):
    # BL says: this is for windows, since tif file
    version = raster3d_version()
    do_tif = False
    if version:
        if version[0] == 2:
            do_tif = True
    if (os.name == 'nt' and do_tif):
        render_exe = "render.exe"
        image_file_name += ".tif"
        image_format = " -tiff "
    else:
        render_exe = "render"
        image_format = " -png "
    r3d_exe = find_exe("render", "CCP4_BIN", "PATH")
    if (r3d_exe):
       r3d_dir = os.path.dirname(r3d_exe)
       os.environ['R3D_LIB'] = r3d_dir + "/materials"
       #we have to check filenames for spaces for dodgy windows path
       image_file_name_mod, source_file_name_mod, space_flag = \
       check_file_names_for_space_and_move(image_file_name, source_file_name)
       r3d_call = r3d_exe + image_format + " -labels" + image_file_name_mod + " < " + source_file_name_mod
       print("BL DEBUG:: r3d_call is ", r3d_call)
       print("calling render...")

       major, minor, micro, releaselevel, serial = sys.version_info
       if (major >= 2 and minor >=4):
           import subprocess
           status = subprocess.call(r3d_call, shell=True)
           if status:
               # something went wrong with raster3d
               # maybe same for system call?!?
               print("BL WARNING:: some error in raster3d")
               return
       else:
           status = os.system(r3d_call)
       # now we have to copy files back if necessary!!
       if (space_flag):
          import shutil
          shutil.move(image_file_name_mod,image_file_name)
          shutil.move(source_file_name_mod,source_file_name)

       print("calling display...")
       try:
         webbrowser.open(image_file_name,1,1)
       except OSError:
         print("BL WARNING:: We can't find rendered file ",image_file_name)

   elif (image_type == "povray"):
    image_file_name += ".png"
    #BL says: conversion of filename, again! Windows is driving me crazy...
    image_file_name = os.path.normpath(image_file_name)
    # again, we just assume povray exe is pvengine on all systems
    povray_exe = find_exe(povray_command_name, "PATH")
    if (povray_exe):
      image_file_name_mod, source_file_name_mod, space_flag = \
        check_file_names_for_space_and_move(image_file_name, source_file_name)
      if (os.name == 'nt'):
          args = ["/EXIT", "/RENDER"]
      else:
          args = " "
      args += [source_file_name_mod] + povray_args() + ["-UV" , "+W" + str(x_size) , "+H" + str(y_size)]
      print("BL INFO:: run povray with args: ", args)
      povray_call = [povray_exe] +  args + ["+o" + image_file_name_mod]
      print("BL DEBUG:: povray command line", povray_call)
      major, minor, micro, releaselevel, serial = sys.version_info
      if (major >= 2 and minor >=4):
          import subprocess
          status = subprocess.call(povray_call)
          if status:
              # something went wrong with raster3d
              # maybe same for system call?!?
              print("BL WARNING:: some error in povray")
              return
      else:
          os.system(povray_call)
      # now we have to copy files back if necessary!!
      if (space_flag):
         import shutil
         shutil.move(image_file_name_mod,image_file_name)
         shutil.move(source_file_name_mod,source_file_name)
      else: pass
      print("calling display...")
      try:
         webbrowser.open(image_file_name,1,1)
      except OSError:
         print("BL WARNING:: We can't find rendered file ",image_file_name)

   else:
     print("Image type ", image_type, " unknown!")

# Converts a ppm file to a bmp file (for windows users) and opens in
# browser or viewer
# actually dont use ppm any more but png, so not needed as such any more
# keep for backwards compatibility
#
def ppm2bmp(ppm_file_name):
    import os
    import webbrowser
    import sys

    bmp_file_name, extension = os.path.splitext(ppm_file_name)
    # FIXME: horrible hacks!!!
    if (extension == ".ppm"):
        bmp_file_name += ".bmp"

        # BL says: need to do some wired stuff to make sure cmd/system works with
        # space in file name , e.g. ' & " thingys
        cmd = 'ppm2bmp "'
        ppm2bmp_call = cmd + ppm_file_name + '"'
        major, minor, micro, releaselevel, serial = sys.version_info
        if (major >= 2 and minor >=4):
            import subprocess
            status = subprocess.call(ppm2bmp_call, shell=True)
            if status:
                # something went wrong with raster3d
                # maybe same for system call?!?
                print("BL WARNING:: some error in ppm2bmp")
                return
        else:
            os.system(ppm2bmp_call)
    if (extension == ".png"):
        bmp_file_name = ppm_file_name
    if (not os.path.isfile(bmp_file_name)):
        print("BL WARNING:: Cannot find png/bmp file ", bmp_file_name)
    else:
        print("calling display...")
        try:
            webbrowser.open(bmp_file_name,1,1)
        except OSError:
            print("BL WARNING:: We can't open screendump file ",bmp_file_name)

# Tests file names for spaces. there is certainly on problem on windows
# not sure about other OS, yet
def check_file_names_for_space_and_move(image_file_name,source_file_name):

    import string, os, shutil

    space_flag = False

    if (((image_file_name.find(" ") > -1) or (source_file_name.find(" ") > -1)) and (os.name == 'nt')):
       # we have spaces, so tmp copy src to C: and run there, then
       # copy tmp back to where it should be
       image_file_name_mod = "C:\\" + os.path.basename(image_file_name)
       source_file_name_mod = "C:\\" + os.path.basename(source_file_name)
       shutil.move(source_file_name,source_file_name_mod)
       space_flag = True
       return image_file_name_mod, source_file_name_mod, space_flag
    else:
       return image_file_name, source_file_name, space_flag

#raytrace("povray","coot.pov","coot.png",600,600)

# return version as a tuple of 3 numbers (or 2 plus letter)
# or False if not found or problems running raster3d
#
def raster3d_version():

    import os

    raster3d_exe = find_exe("render", "CCP4_BIN", "PATH")

    if not raster3d_exe:
        return False
    else:
        log_file = "raster3_version.tmp"
        # BL note: currently -i is a bogus switch, so gives info
        status = popen_command(raster3d_exe, ["-i"], [], log_file)
        if (status != 0):
            return False
        else:
            fin = open(log_file, 'r')
            lines = fin.readlines()
            fin.close()
            os.remove(log_file)
            for line in lines:
                if ("Raster3D" in line):
                    version_string = line[9:]
                    tmp = version_string.split(".")
                    try:
                        major_version = int(tmp[0])
                    except:
                        print("BL INFO:: problem extracting major version " + \
                              "from raster3d")
                        return False
                    if ("-" in tmp[1]):
                        # have new style version
                        tmp_min = tmp[1].split("-")
                        try:
                            minor = int(tmp_min[0])
                            micro = int(tmp_min[1])
                        except:
                            print("BL INFO:: problem extracting minor version " + \
                                  "from raster3d")

                            return False
                        return [major, minor, micro]
                    else:
                        # old style
                        if (len(tmp[1]) != 2):
                            print("BL INFO:: cannot deal with this version.")
                            return False
                        else:
                            try:
                                minor = int(tmp[1][0])
                            except:
                                print("BL INFO:: problem extracting minor " + \
                                      "version of raster3d (old style).")
                                return False
                            micro = tmp[1][1]
                            return [major, minor, micro]
            return False
