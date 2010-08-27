/*
     pygl/help.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC

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
#if defined (_WIN32)
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

// These are  helper functions which don't depend on an OPENGL context.
image_info get_pixdata(int trans=0);
void write_pixdata(const char *filename, int width=-1, int height=-1, int trans=0);
int findprimc(const std::vector<Cartesian> &xyzbox, const std::vector<Cartesian> &primorigin, const Cartesian &origin, const matrix &objrotmat);
int findprimc(double *xyzmpc, double *xyzmmc, double *xyzpmc, double *xyzppc, std::vector<Cartesian> primorigin, Cartesian origin, matrix objrotmat);
int findprimc_main(const Cartesian &xyzmpf, const Cartesian &xyzmpb, const Cartesian &xyzmmf, const Cartesian &xyzmmb, const Cartesian &xyzpmf, const Cartesian &xyzpmb, const Cartesian &xyzppf, const Cartesian &xyzppb, const std::vector<Cartesian> &primorigin,const Cartesian &origin,const matrix &objrotmat);
std::vector<Cartesian>getxyzc(double x,double y);
const double *GLf2f(const GLfloat *in, int size);
const GLfloat *f2GLf(double *in, int size);
const GLfloat *buildrotmatrix_from_c(matrix a);
const GLfloat *buildrotmatrix(
GLfloat a0, GLfloat a1, GLfloat a2, GLfloat a3, 
GLfloat a4, GLfloat a5, GLfloat a6, GLfloat a7, 
GLfloat a8, GLfloat a9, GLfloat a10, GLfloat a11, 
GLfloat a12, GLfloat a13, GLfloat a14, GLfloat a15);

int CheckIfStereoAvailable(void);
int CheckIfAlphaAvailable(void);

image_info_yuv_t get_yuvdata(int trans=0);

Volume GetClippingPlanes();
Volume GetFrontAndBackClippingPlanes();
bool isPointInClippingVolume(const Cartesian &p, const Volume &v);

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

#endif
