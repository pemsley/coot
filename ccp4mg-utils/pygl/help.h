/*
     pygl/help.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2011 University of York
     Copyright (C) 2012 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/
#if defined (_WIN32) && not defined (WINDOWS_MINGW)
#include <windows.h>
#if !defined (__GNUC__)
#define snprintf _snprintf
#endif
#endif


#ifndef __GL_HELP__

#define __GL_HELP__
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#ifdef USE_GLX
#include <GL/glx.h>
#endif


#include "ppmutil.h"
#include <vector>

#include "matrix.h"
#include "cartesian.h"
#include "plane.h"
#include "volume.h"

#include <vector>

class OffScreenBuffer{
  unsigned width;
  unsigned height;
#if defined (USE_GLX)
  Pixmap pixmap;
  GLXPixmap glxpixmap;
  Display *pixmapdpy;
  GLXContext offcontext;
#elif defined (_WIN32)
  HDC hDC;
  HGLRC hGLRC;
  HBITMAP hBitmap, hOldBitmap;
  void setupPixelFormat();
  void setupDIB();
#endif
#if defined (__APPLE_CC__) // These not used at moment.
/*
CGLContextObj ctx;
CGLContextObj glut_ctx;
unsigned char *pixels;
*/
#endif
 public:
  OffScreenBuffer(unsigned w, unsigned h);
  ~OffScreenBuffer();
  void MakeCurrent();
  void Write(const char *filename, int width=-1, int height=-1, int trans=0);
  image_info get_pixdata(int trans=0);
  unsigned GetWidth() const {return width;};
  unsigned GetHeight() const {return height;};
};

class OnScreenBuffer{
#if defined (USE_GLX)
  unsigned long readable;
  unsigned long drawable;
  Display *dpy;
  GLXContext context;
#elif defined (_WIN32)
  HDC hDC;
  HGLRC hGLRC;
#endif
 public:
  OnScreenBuffer();
  ~OnScreenBuffer(){};
  void MakeCurrent();
};

void SetStereoAvailable(int stereo_available_in);

#if defined (_WIN32) && not defined (WINDOWS_MINGW)
#define EXAMPLE_DLL __declspec(dllexport)
#endif
class MGFramebufferObject {
 protected:
  GLuint _texture;
  GLuint _id;
  int _width;
  int _height;
 public:
#if defined (_WIN32) && not defined (WINDOWS_MINGW)
  __stdcall EXAMPLE_DLL MGFramebufferObject();
  __stdcall EXAMPLE_DLL ~MGFramebufferObject();
  __stdcall EXAMPLE_DLL MGFramebufferObject(int theWidth, int theHeight);
  unsigned __stdcall EXAMPLE_DLL texture();
  unsigned __stdcall EXAMPLE_DLL handle();
  int __stdcall EXAMPLE_DLL width();
  int __stdcall EXAMPLE_DLL height();
  void __stdcall EXAMPLE_DLL bind();
  void __stdcall EXAMPLE_DLL release();
#else
  MGFramebufferObject();
  ~MGFramebufferObject();
  MGFramebufferObject(int theWidth, int theHeight);
  unsigned texture();
  unsigned handle();
  int width();
  int height();
  void bind();
  void release();
#endif
};

class MGDepthFramebufferObject : public MGFramebufferObject {
 public:
#if defined (_WIN32) && not defined (WINDOWS_MINGW)
   __stdcall EXAMPLE_DLL MGDepthFramebufferObject(int theWidth, int theHeight);
   __stdcall EXAMPLE_DLL ~MGDepthFramebufferObject();
#else
   MGDepthFramebufferObject(int theWidth, int theHeight);
   ~MGDepthFramebufferObject();
#endif
};

class MGDepthCompareFramebufferObject : public MGFramebufferObject {
 public:
#if defined (_WIN32) && not defined (WINDOWS_MINGW)
   __stdcall EXAMPLE_DLL MGDepthCompareFramebufferObject(int theWidth, int theHeight);
   __stdcall EXAMPLE_DLL ~MGDepthCompareFramebufferObject();
#else
   MGDepthCompareFramebufferObject(int theWidth, int theHeight);
   ~MGDepthCompareFramebufferObject();
#endif
};

#endif
