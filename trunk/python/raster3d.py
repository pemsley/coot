#BL pyton script for rendering images in COOT

def povray_args():
    return " +FN16 +A"
# BL says: dont know how usefull this function is/will be....

def render_image():
    import os, webbrowser

    coot_r3d_file_name="coot.r3d"
    coot_tif_file_name="coot.tif"
    raster3d(coot_r3d_file_name)
    if os.name == "nt":
       render_exe = "render.exe"
    else:
       render_exe = "render"
    r3d_exe = find_render_exe("raster3d",render_exe)
    r3d_dir = os.path.dirname(r3d_exe)
    if (r3d_exe):
       os.environ['R3D_LIB'] = r3d_dir + "/materials"
       r3d_call = r3d_exe + " -tiff " + coot_tif_file_name + " < " + coot_r3d_file_name
       print "calling render..."
       
# we assume for now that render is in path
       status = os.system(r3d_call)
       print "calling display..."
       try:
         webbrowser.open(coot_tif_file_name)
       except OSError:
         print "BL WARNING:: We can't find rendered file ",coot_tif_file_name

def raytrace(image_type, source_file_name, image_file_name, x_size, y_size):
   import os, webbrowser

   if (image_type == "raster3d"):
# BL says: this is for windows, since tif file
    if os.name == "nt":
       render_exe = "render.exe"
       image_file_name += ".tif"
    else:
       render_exe = "render"
       image_file_name += ".png"
    r3d_exe = find_render_exe("raster3d",render_exe)
    r3d_dir = os.path.dirname(r3d_exe)
    if (r3d_exe):
       os.environ['R3D_LIB'] = r3d_dir + "/materials"
       r3d_call = r3d_exe + " -tiff " + image_file_name + " < " + source_file_name
       print "calling render..."

# we assume for now that render is in path
       status = os.system(r3d_call)
#       print "Status ",status

       print "calling display..."
       try:
         webbrowser.open(image_file_name)
       except OSError:
         print "BL WARNING:: We can't find rendered file ",image_file_name

   elif (image_type == "povray"):
    image_file_name += ".png"
#BL says: conversion of filename, again! Windows is driving me crazy...
    image_file_name = os.path.normpath(image_file_name)
    args = source_file_name + povray_args() + " +W" + str(x_size) + " +H" + str(y_size)
    print "INFO:: povray with args: ", args
# BL says: this should work for all POV-Ray installations, as long as it 
# has POV in name
    if os.name == "nt":
       render_exe = "pvengine.exe"
    else:
# BL says: dunno what povray exe is called on other systems, just assume for now
       render_exe = "pvengine"
    povray_exe = find_render_exe("POV*",render_exe)
    if (povray_exe):
      povray_call = povray_exe + " /EXIT /RENDER " + args + " +o" + image_file_name
      print "BL DEBUG:: povray command line", povray_call
      os.system(povray_call)
      print "calling display..."
      try:
         webbrowser.open(image_file_name)
      except OSError:
         print "BL WARNING:: We can't find rendered file ",image_file_name

   else:
     print "Image type ", image_type, " unknown!"


def find_render_exe(render_program,render_program_exe):
    import os, string, fnmatch
    render_exe = None
    render_dir = None

    if (os.name=="nt"):
# for now we assume everything will be in C. Could be more flexible in windows but would require win32api
# which I dont wanna use since it's not official python!
# change C to whatever other name if your raster3d is not ther or put it in your path!!
      drive = "C:\\"
    else:
      drive = "/"
# first try to find it in PATH
    for path in string.split(os.environ["PATH"], os.pathsep):
        if (fnmatch.fnmatch(os.path.basename(path),render_program)):
           render_dir = path
           if (render_program=="raster3d"):
              render_exe = os.path.join(render_dir,render_program_exe)
              if (os.path.isfile(render_exe)):
                 return render_exe
                 break
           elif (render_program=="POV*"):
              render_exe = os.path.join(render_dir,"bin",render_program_exe)
              if (os.path.isfile(render_exe)):
                 return render_exe
                 break
           else:
              print "BL INFO:: We dont dont know ", render_program
    if (render_dir):
        print "We found ", render_program , " in ($PATH): ",render_dir
        pass
    else:
#if not there search everywhere
     for root, dir, files in os.walk(drive):
        if (fnmatch.fnmatch(os.path.basename(root),render_program)):
           render_dir = root
           if (render_program=="raster3d"):
              render_exe = os.path.join(render_dir,render_program_exe)
              if (os.path.isfile(render_exe)):
                 return render_exe
                 break
           elif (render_program=="POV*"):
              render_exe = os.path.join(render_dir,"bin",render_program_exe)
              if (os.path.isfile(render_exe)):
                 return render_exe
                 break
           else:
              print "BL INFO:: We dont dont know ", render_program
    if (render_dir):
       print "We found ", render_program, " not in $PATH, but in: ",render_dir
       pass
#we havent found raster3d dir at all:
    else:
       print "BL WARNING:: We havent found ", render_program, "!!"
       return render_exe

def ppm2bmp(ppm_file_name):
    import os, webbrowser

    bmp_file_name ,extension = os.path.splitext(ppm_file_name)
    bmp_file_name += ".bmp"
#    ppm2bmp_call = "ppm2bmp " + ppm_file_name + " " + bmp_file_name
    ppm2bmp_call = "ppm2bmp " + ppm_file_name 
    os.system(ppm2bmp_call)
    if (os.path.isfile(bmp_file_name)):
      print "calling display..."
      try:
         webbrowser.open(bmp_file_name)
      except OSError:
         print "BL WARNING:: We can't find screendump file ",bmp_file_name
    else:
      print "BL INFO:: Cannot find bmp file ",bmp_file_name

#testing:
#raytrace("povray","coot.pov","coot.png",600,600)
