/*
     pygl/help.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

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

#ifdef USE_GLX
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#endif

#include <math.h>

#define GLUT_DISABLE_ATEXIT_HACK

#ifdef __APPLE_CC__
#include <OpenGL/OpenGL.h>
#include <OpenGL/glu.h>
#include <ApplicationServices/ApplicationServices.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#ifdef _USE_GLUT_
#include <GL/glut.h>
#endif

#ifdef USE_GLX
#include <GL/glx.h>
#ifndef sgi
#include <GL/glext.h>
#endif
#endif

/*
#ifdef __APPLE_CC__
#include <OpenGL/OpenGL.h>
#include <OpenGL/CGLTypes.h>
#include <OpenGL/CGLContext.h>
#endif
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include "help.h"
#include "ppmutil.h"
#include "plane.h"
#include "cartesian.h"
#include "matrix.h"
#include "geomutil.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#define PIBY2 (M_PI * 2)

using namespace std;

#ifdef _WIN32
void Win32Error(void){
    DWORD dw = GetLastError();
    LPVOID lpMsgBuf;
    FormatMessage( 
    FORMAT_MESSAGE_ALLOCATE_BUFFER | 
    FORMAT_MESSAGE_FROM_SYSTEM | 
    FORMAT_MESSAGE_IGNORE_INSERTS,
    NULL,
    dw,
    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
    (LPTSTR) &lpMsgBuf,
    0,
    NULL );
    printf("Error: %ld %s\n",dw,(LPTSTR)lpMsgBuf);
}

void OffScreenBuffer::setupDIB(){
    BITMAPINFO *bmInfo;
    BITMAPINFOHEADER *bmHeader;
    UINT usage;
    VOID *base;
    int bmiSize;
    int bitsPerPixel;

    bmiSize = sizeof(*bmInfo);
    bitsPerPixel = GetDeviceCaps(hDC, BITSPIXEL);

    switch (bitsPerPixel) {
    case 8:
	/* bmiColors is 256 WORD palette indices */
	bmiSize += (256 * sizeof(WORD)) - sizeof(RGBQUAD);
	break;
    case 16:
	/* bmiColors is 3 WORD component masks */
	bmiSize += (3 * sizeof(DWORD)) - sizeof(RGBQUAD);
	break;
    case 24:
    case 32:
    default:
	/* bmiColors not used */
	break;
    }

    bmInfo = (BITMAPINFO *) calloc(1, bmiSize);
    bmHeader = &bmInfo->bmiHeader;

    bmHeader->biSize = sizeof(*bmHeader);
    bmHeader->biWidth = width;
    bmHeader->biHeight = height;
    bmHeader->biPlanes = 1;			/* must be 1 */
    bmHeader->biBitCount = bitsPerPixel;
    bmHeader->biXPelsPerMeter = 0;
    bmHeader->biYPelsPerMeter = 0;
    bmHeader->biClrUsed = 0;			/* all are used */
    bmHeader->biClrImportant = 0;		/* all are important */

    switch (bitsPerPixel) {
    case 8:
	bmHeader->biCompression = BI_RGB;
	bmHeader->biSizeImage = 0;
	usage = DIB_PAL_COLORS;
	/* bmiColors is 256 WORD palette indices */
	{
	    WORD *palIndex = (WORD *) &bmInfo->bmiColors[0];
	    int i;

	    for (i=0; i<256; i++) {
		palIndex[i] = i;
	    }
	}
	break;
    case 16:
	bmHeader->biCompression = BI_RGB;
	bmHeader->biSizeImage = 0;
	usage = DIB_RGB_COLORS;
	/* bmiColors is 3 WORD component masks */
	{
	    DWORD *compMask = (DWORD *) &bmInfo->bmiColors[0];

	    compMask[0] = 0xF800;
	    compMask[1] = 0x07E0;
	    compMask[2] = 0x001F;
	}
	break;
    case 24:
    case 32:
    default:
	bmHeader->biCompression = BI_RGB;
	bmHeader->biSizeImage = 0;
	usage = DIB_RGB_COLORS;
	/* bmiColors not used */
	break;
    }

    hBitmap = CreateDIBSection(hDC, bmInfo, usage, &base, NULL, 0);
    if (hBitmap == NULL) {
	(void) MessageBox(WindowFromDC(hDC),
		"Failed to create DIBSection.",
		"OpenGL application error",
		MB_ICONERROR | MB_OK);
	exit(1);
    }

    hOldBitmap = (HBITMAP)SelectObject(hDC, hBitmap);

    free(bmInfo);
}

void OffScreenBuffer::setupPixelFormat()
{
    PIXELFORMATDESCRIPTOR pfd = {
	sizeof(PIXELFORMATDESCRIPTOR),	/* size of this pfd */
	1,				/* version num */
	PFD_SUPPORT_OPENGL,		/* support OpenGL */
	0,				/* pixel type */
	24,				/* 8-bit color depth */
	0, 0, 0, 0, 0, 0,		/* color bits (ignored) */
	0,				/* no alpha buffer */
	0,				/* alpha bits (ignored) */
	0,				/* no accumulation buffer */
	0, 0, 0, 0,			/* accum bits (ignored) */
	32,				/* depth buffer */
	0,				/* no stencil buffer */
	0,				/* no auxiliary buffers */
	PFD_MAIN_PLANE,			/* main layer */
	0,				/* reserved */
	0, 0, 0,			/* no layer, visible, damage masks */
    };
    int SelectedPixelFormat;
    BOOL retVal;

    pfd.cColorBits = GetDeviceCaps(hDC, BITSPIXEL);

    pfd.iPixelType = PFD_TYPE_RGBA;

    pfd.dwFlags |= PFD_DRAW_TO_BITMAP;

    SelectedPixelFormat = ChoosePixelFormat(hDC, &pfd);
    if (SelectedPixelFormat == 0) {
	(void) MessageBox(WindowFromDC(hDC),
		"Failed to find acceptable pixel format.",
		"OpenGL application error",
		MB_ICONERROR | MB_OK);
	return;
    }

    retVal = SetPixelFormat(hDC, SelectedPixelFormat, &pfd);
    if (retVal != TRUE) {
	(void) MessageBox(WindowFromDC(hDC),
		"Failed to set pixel format.",
		"OpenGL application error",
		MB_ICONERROR | MB_OK);
	return;
    }
}
#endif

#if defined (USE_GLX) 
#include <GL/glx.h>


image_info GetDrawableImage(Drawable win, Display *dpy, unsigned width, unsigned height){

  // XY or Z determines speed. Not sure how to determine which to use though.
  XImage *image = XGetImage (dpy, win, 0, 0, width, height, AllPlanes, ZPixmap);

  XColor *rgbs = new XColor[width];

  Colormap colormap = DefaultColormap(dpy,DefaultScreen(dpy));
  unsigned char *bitmap = new unsigned char[width*height*3];

  for (unsigned i = 0; i < height; i++){
    for (unsigned j = 0; j < width; j++){
       rgbs[j].pixel =  XGetPixel(image, j, i);
    }
    XQueryColors(dpy,colormap,rgbs,width);
    for (unsigned j = 0; j < width; j++){
      bitmap[(i*width+j)*3]   = rgbs[j].red;
      bitmap[(i*width+j)*3+1] = rgbs[j].green;
      bitmap[(i*width+j)*3+2] = rgbs[j].blue;
    }
  }

  image_info iinfo(width, height, bitmap);

  return iinfo;
}
#endif

void OffScreenBuffer::MakeCurrent(void){
#if defined (USE_GLX)
  glXMakeCurrent(pixmapdpy,glxpixmap,offcontext);
#elif defined (_WIN32)
  wglMakeCurrent( hDC, hGLRC );
#endif
}

image_info OffScreenBuffer::get_pixdata(int trans){
  image_info iinfo;
#if defined (USE_GLX)
  if(!pixmap||!glxpixmap)
    return iinfo;
#endif
  iinfo = ::get_pixdata(trans);
  return iinfo;
}

void OffScreenBuffer::Write(const char *filename, int width, int height,int trans){
  image_info iinfo = get_pixdata(trans);
  iinfo.invert();
  if(width>0&&height>0) iinfo.ScaleImage(width,height);
  iinfo.write(filename);
}

OffScreenBuffer::~OffScreenBuffer(){
#if defined (USE_GLX)
  glXDestroyContext(pixmapdpy,offcontext);
  glXDestroyGLXPixmap(pixmapdpy,glxpixmap);
  XFreePixmap(pixmapdpy, pixmap);
  XCloseDisplay(pixmapdpy);
#endif
}

OffScreenBuffer::OffScreenBuffer(unsigned w, unsigned h){
#if defined (USE_GLX) //&& !defined (__APPLE_CC__)

  width = w;
  height = h;

  int attribSingle_acc[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ACCUM_RED_SIZE, 1,
      GLX_ACCUM_GREEN_SIZE, 1,
      GLX_ACCUM_BLUE_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      None };
  int attribDouble_acc[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ACCUM_RED_SIZE, 1,
      GLX_ACCUM_GREEN_SIZE, 1,
      GLX_ACCUM_BLUE_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      GLX_DOUBLEBUFFER, 1,
      None };
  int attribSingle[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      None };
  int attribDouble[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      GLX_DOUBLEBUFFER, 1,
      None };
  int attribSingle_acc_alpha[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ALPHA_SIZE, 1,
      GLX_ACCUM_RED_SIZE, 1,
      GLX_ACCUM_GREEN_SIZE, 1,
      GLX_ACCUM_BLUE_SIZE, 1,
      GLX_ACCUM_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      None };
  int attribDouble_acc_alpha[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ALPHA_SIZE, 1,
      GLX_ACCUM_RED_SIZE, 1,
      GLX_ACCUM_GREEN_SIZE, 1,
      GLX_ACCUM_BLUE_SIZE, 1,
      GLX_ACCUM_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      GLX_DOUBLEBUFFER, 1,
      None };
  int attribSingle_alpha[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      None };
  int attribDouble_alpha[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      GLX_DOUBLEBUFFER, 1,
      None };

  pixmapdpy = XOpenDisplay(0);

  Display *dpy = pixmapdpy;

  XVisualInfo *visinfo;
  visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribSingle_acc_alpha);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribDouble_acc_alpha);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribSingle_alpha);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribDouble_alpha);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribSingle_acc);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribDouble_acc);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribSingle);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribDouble);
  if (!visinfo) {
    printf("Cannot choose double or single visual\n");
    XCloseDisplay(dpy);
    return;
  }

  pixmap = XCreatePixmap (dpy, RootWindow(dpy,DefaultScreen(dpy)), width, height, DefaultDepth(dpy,DefaultScreen(dpy)));
  if(!pixmap){
   XCloseDisplay(dpy);
   printf("Cannot create pixmap\n");
   return;
  }
  glxpixmap = glXCreateGLXPixmap(dpy,visinfo,pixmap);

  if(!glxpixmap){
   XCloseDisplay(dpy);
   printf("Cannot create glxpixmap\n");
   return;
  }

  offcontext = glXCreateContext(dpy, visinfo, 0, GL_FALSE);

#elif defined (__APPLE_CC__)
/*
  if(!glut_ctx)
    glut_ctx = CGLGetCurrentContext();

  CGLPixelFormatObj pix;
  long nvirt;

  CGLPixelFormatAttribute attribs[] = {
    kCGLPFAOffScreen,
    kCGLPFAColorSize, 32,
    0
  };


  CGLChoosePixelFormat(attribs, &pix, &nvirt);

  printf("Number of pixel formats:%d\n",nvirt);

  CGLCreateContext(pix, 0, &ctx);
  if(!ctx)
    printf("Error creating context\n");

  CGLDestroyPixelFormat(pix);
  pixels = new unsigned char[ w * h * 4 ];

  CGLError error = CGLSetOffScreen(ctx, w, h, 800, pixels);
*/
#elif defined (_WIN32)
  width = w;
  height = h;
  hDC = GetDC(NULL);
  HDC hDCFrontBuffer;
  hDCFrontBuffer = hDC;
  hDC = CreateCompatibleDC(hDCFrontBuffer);
  setupDIB();
  setupPixelFormat();
  hGLRC = wglCreateContext(hDC);
  if (hGLRC == NULL) {
     MessageBox(WindowFromDC(hDC),
     "Failed to create OpenGL context.",
     "OpenGL application error",
     MB_ICONERROR | MB_OK);
     exit(1);
  }
#endif
}

#if defined (_WIN32)
static int stereo_available;
#endif

void SetStereoAvailable(int stereo_available_in){
#if defined (_WIN32)
  stereo_available = stereo_available_in;
#endif
}

int CheckIfAlphaAvailable(void){
#ifdef USE_GLX
  int attribSingle_alpha[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      None };
  int attribDouble_alpha[32] = {
      GLX_RGBA,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      GLX_DOUBLEBUFFER, 1,
      None };

  Display *dpy = XOpenDisplay(0);
  XVisualInfo *visinfo;
  visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribSingle_alpha);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribDouble_alpha);
  if(visinfo){
    XCloseDisplay(dpy);
    return 1;
  }
  XCloseDisplay(dpy);
  return 0;
#elif defined (_WIN32)
  return 1;  // Uh-oh!
#else
  return 0;
#endif
}

int CheckIfStereoAvailable(void){
#ifdef USE_GLX
  int attribSingle[32] = {
      GLX_RGBA,
      GLX_STEREO,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      None };
  int attribDouble[32] = {
      GLX_RGBA,
      GLX_STEREO,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      GLX_DOUBLEBUFFER, 1,
      None };
  int attribSingle_alpha[32] = {
      GLX_RGBA,
      GLX_STEREO,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      None };
  int attribDouble_alpha[32] = {
      GLX_RGBA,
      GLX_STEREO,
      GLX_RED_SIZE, 1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE, 1,
      GLX_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      GLX_DOUBLEBUFFER, 1,
      None };

  Display *dpy = XOpenDisplay(0);
  XVisualInfo *visinfo;
  visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribSingle);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribDouble);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribSingle_alpha);
  if(!visinfo) visinfo = glXChooseVisual(dpy, DefaultScreen(dpy), attribDouble_alpha);
  if(visinfo){
    XCloseDisplay(dpy);
    return 1;
  }
  XCloseDisplay(dpy);
  return 0;
#elif defined (_WIN32)
  WNDCLASS wc;
  HWND hWnd;
  HDC hDC;
  //HGLRC hRC;
	
  // register window class
  wc.style = CS_OWNDC;
  wc.lpfnWndProc = DefWindowProc;
  wc.cbClsExtra = 0;
  wc.cbWndExtra = 0;
  wc.hInstance = GetModuleHandle(NULL);
  wc.hIcon = LoadIcon( NULL, IDI_APPLICATION );
  wc.hCursor = LoadCursor( NULL, IDC_ARROW );
  wc.hbrBackground = (HBRUSH)GetStockObject( BLACK_BRUSH );
  wc.lpszMenuName = NULL;
  wc.lpszClassName = "GLSample";
  RegisterClass( &wc );
	
  // create main window
  hWnd = CreateWindow("GLSample", "OpenGL Sample", 
                      WS_CAPTION | WS_POPUPWINDOW,
                      0, 0, 10, 10,
                      NULL, NULL, GetModuleHandle(NULL), NULL );

  PIXELFORMATDESCRIPTOR pfd;
  int format;
	
  // get the device context (DC)
  hDC = GetDC( hWnd );
	
  // set the pixel format for the DC
  ZeroMemory( &pfd, sizeof( pfd ) );
  pfd.nSize = sizeof( pfd );
  pfd.nVersion = 1;
  pfd.dwFlags = PFD_STEREO | PFD_SUPPORT_GDI| PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    
  format = ChoosePixelFormat( hDC, &pfd );
    
  ShowWindow(hWnd,SW_SHOWDEFAULT);
  BOOL retVal = SetPixelFormat( hDC, format, &pfd );
  ShowWindow(hWnd,SW_HIDE);
  if (retVal != TRUE) {
	(void) MessageBox(WindowFromDC(hDC),
		"Failed to set pixel format.",
		"OpenGL application error",
		MB_ICONERROR | MB_OK);
        DestroyWindow(hWnd);
	return 0;
  }
  //BOOL bStereoAvailable;

  format = GetPixelFormat (hDC);
  DescribePixelFormat (hDC, format, sizeof(PIXELFORMATDESCRIPTOR), &pfd);

  if ((pfd.dwFlags & PFD_STEREO) == 0){
    //printf("Windows does not think stereo is available\n"); fflush(stdout);
    DestroyWindow(hWnd);
    return 0;
  }else{
    //printf("Windows thinks stereo is available\n"); fflush(stdout);
    DestroyWindow(hWnd);
    return 1;
  }
  DestroyWindow(hWnd);
  return 0;
#elif defined(__APPLE_CC__)
  CGLPixelFormatObj pixelFormatObj ;
  GLint numPixelFormats ;
  CGOpenGLDisplayMask displayMask = CGDisplayIDToOpenGLDisplayMask( CGMainDisplayID() ) ;
  CGLPixelFormatAttribute attribs[] = { (CGLPixelFormatAttribute)kCGLPFADisplayMask, (CGLPixelFormatAttribute)displayMask, (CGLPixelFormatAttribute)kCGLPFAStereo, (CGLPixelFormatAttribute)0 } ;
  /*CGLError err =*/ CGLChoosePixelFormat( attribs, &pixelFormatObj, &numPixelFormats );
  if(numPixelFormats>0) return 1;
  return 0;
#else
  return 0;
#endif
}

void OnScreenBuffer::MakeCurrent(){
#ifdef USE_GLX
  glXMakeCurrent(dpy,drawable,context);
#elif defined (_WIN32)
  wglMakeCurrent( hDC, hGLRC );
#endif
}
OnScreenBuffer::OnScreenBuffer(){
#ifdef USE_GLX
  dpy = XOpenDisplay(0);
  readable = glXGetCurrentDrawable();
  drawable = glXGetCurrentReadDrawable();
  context = glXGetCurrentContext();
#elif defined (_WIN32)
  hDC = GetDC(NULL);
  hGLRC = wglGetCurrentContext();
#endif
}

image_info_yuv_t get_yuvdata(int trans){
  /*
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  glPixelTransferi(GL_MAP_COLOR,GL_FALSE);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  printf("width:%d, height:%d\n",viewport[2],viewport[3]); fflush(stdout);
  unsigned char* pixels = new unsigned char[viewport[2]*viewport[3]*IMAGEINFO_RGB_SIZE];
  glReadPixels(0, 0, viewport[2], viewport[3], GL_RGB, GL_UNSIGNED_BYTE, pixels);
  image_info iinfo(viewport[2], viewport[3], pixels,IMAGEINFO_RGB);
  */
  image_info iinfo = get_pixdata(trans);
  iinfo.invert();
  
  return iinfo.getyuv();

}

image_info get_pixdata(int trans){

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  glPixelTransferi(GL_MAP_COLOR,GL_FALSE);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  image_info iinfo;
  unsigned char* pixels;
  if(!trans){
  pixels = new unsigned char[viewport[2]*viewport[3]*IMAGEINFO_RGBA_SIZE];
  glReadPixels(0, 0, viewport[2], viewport[3], GL_RGBA, GL_UNSIGNED_BYTE,pixels);
  iinfo = image_info(viewport[2], viewport[3], pixels,IMAGEINFO_RGBA);
  }else{
  pixels = new unsigned char[viewport[2]*viewport[3]*IMAGEINFO_RGB_SIZE];
  glReadPixels(0, 0, viewport[2], viewport[3], GL_RGB, GL_UNSIGNED_BYTE,pixels);
  iinfo = image_info(viewport[2], viewport[3], pixels,IMAGEINFO_RGB);
  }

  delete [] pixels;

  return iinfo;
}

void write_pixdata(const char *filename, int width, int height, int trans){

  image_info iinfo = get_pixdata(trans);

  iinfo.invert();
  if(width>0&&height>0) iinfo.ScaleImage(width,height);
  iinfo.write(filename);

}

int findprimc(const std::vector<Cartesian> &xyzbox, const std::vector<Cartesian> &primorigin, const Cartesian &origin, const matrix &objrotmat){
  return findprimc_main(xyzbox[0],xyzbox[1],xyzbox[2],xyzbox[3],xyzbox[4],xyzbox[5],xyzbox[6],xyzbox[7],primorigin,origin,objrotmat);
}

int findprimc(double *xyzmpc, double *xyzmmc, double *xyzpmc, double *xyzppc, std::vector<Cartesian> primorigin, Cartesian origin, matrix objrotmat){

  Cartesian xyzmpf(xyzmpc);
  Cartesian xyzmpb(xyzmpc+3);
  Cartesian xyzmmf(xyzmmc);
  Cartesian xyzmmb(xyzmmc+3);
  Cartesian xyzpmf(xyzpmc);
  Cartesian xyzpmb(xyzpmc+3);
  Cartesian xyzppf(xyzppc);
  Cartesian xyzppb(xyzppc+3);
  return findprimc_main(xyzmpf,xyzmpb,xyzmmf,xyzmmb,xyzpmf,xyzpmb,xyzppf,xyzppb,primorigin,origin,objrotmat);
}

int findprimc_main(const Cartesian &xyzmpf, const Cartesian &xyzmpb, const Cartesian &xyzmmf, const Cartesian &xyzmmb, const Cartesian &xyzpmf, const Cartesian &xyzpmb, const Cartesian &xyzppf, const Cartesian &xyzppb, const std::vector<Cartesian> &primorigin,const Cartesian &origin,const matrix &objrotmat){
  unsigned int j;
  Cartesian prim;

  Cartesian front = Cartesian::MidPoint(objrotmat*xyzmmf,objrotmat*xyzppf);
  Cartesian back  = Cartesian::MidPoint(objrotmat*xyzmmb,objrotmat*xyzppb);

  std::vector<Cartesian> points;
  std::vector<Cartesian>::const_iterator point;

  std::vector<Plane> planes;
  std::vector<Plane>::iterator plane;

  points.push_back(xyzmpf);
  points.push_back(xyzmmf);
  points.push_back(xyzpmf);
  points.push_back(xyzppf);

  //points.push_back(xyzmpf); // Front
  //points.push_back(xyzmpb); // Back

  planes.push_back(Plane(xyzmpf,xyzmpb,xyzmmf));
  planes.push_back(Plane(xyzmmf,xyzmmb,xyzpmf));
  planes.push_back(Plane(xyzpmf,xyzpmb,xyzppf));
  planes.push_back(Plane(xyzppf,xyzppb,xyzmpf));

  //planes.push_back(Plane(xyzmpf,xyzmmf,xyzppf)); // Front
  //planes.push_back(Plane(xyzmpb,xyzppb,xyzmmb)); // Back

  std::vector<Cartesian> fb_points;
  std::vector<Plane> fb_planes;
  Volume v = GetClippingPlanes();
  fb_planes.push_back(v.GetPlane(5));
  fb_planes.push_back(v.GetPlane(4));
  fb_points.push_back(fb_planes[0].find_points_on_plane()[0]);
  fb_points.push_back(fb_planes[1].find_points_on_plane()[0]);

  Cartesian p2prim;
  Cartesian n;

  int clicked = 0;
  int clicked_prim = -1;

  double mindist = 1.0e+8;

  for(j=0;j<primorigin.size();j++){
    prim = primorigin[j];
    clicked = 1;
    plane = fb_planes.begin();
    point = fb_points.begin();
    while(plane!=fb_planes.end()&&clicked){
      n = plane->get_normal();
      n.normalize();
      p2prim = *point-prim;
      if(n.DotProduct(n,p2prim)>0.0){
        clicked = 0;
      }
      point++;
      plane++;
    }
    
    /* Now apply rotation matrix */
    prim = objrotmat*prim;
    prim += origin;
    
    plane = planes.begin();
    point = points.begin();
    //std::cout << "prim: " << prim << std::endl; 
    while(plane!=planes.end()&&clicked){
      n = plane->get_normal();
      //std::cout << "point: " << *point << std::endl; 
      //std::cout << n << std::endl;
      p2prim = *point-prim;
      //std::cout << "p2prim: " << p2prim << std::endl; 
      n.normalize();
      p2prim.normalize();
      //std::cout <<  n.DotProduct(n,p2prim) << std::endl;
      if(n.DotProduct(n,p2prim)>1e-3){
	clicked = 0;
      }
      point++;
      plane++;
    }
    if(clicked) {
      std::vector<double> linedist = DistanceBetweenPointAndLine(front,back,prim);
      double dist = fabs(linedist[0]);
      if(dist<mindist){
        clicked_prim = j;
	mindist = dist;
      }
    }
  }
  //std::cout << clicked_prim << " " << mindist << std::endl;

  return clicked_prim;
  
}

std::vector<Cartesian>getxyzc(double x,double y){

  GLint viewport[4];
  GLdouble mvmatrix[16];
  GLdouble projmatrix[16];
  GLdouble wx, wy, wz;
  GLint realy; /* OpenGL y coordinate position. Why do we do this? */

  glGetIntegerv(GL_VIEWPORT,viewport);
  glGetDoublev(GL_PROJECTION_MATRIX,projmatrix);
  glGetDoublev(GL_MODELVIEW_MATRIX,mvmatrix);
  
  realy = viewport[3] - (GLint)y - 1; /* See above */

  std::vector<Cartesian> ps;

  gluUnProject((GLdouble)x,(GLdouble)realy,0.0,mvmatrix,projmatrix,viewport,&wx,&wy,&wz);
  ps.push_back(Cartesian(wx,wy,wz));

  gluUnProject((GLdouble)x,(GLdouble)realy,1.0,mvmatrix,projmatrix,viewport,&wx,&wy,&wz);
  ps.push_back(Cartesian(wx,wy,wz));

  return ps;

}

const double *GLf2f(const GLfloat *in, int size){
   int i;
   double *out;

   //out = (double *)malloc(size*sizeof(double));
   out = new double [size];
   for(i=0;i<size;i++)
     out[i] = (double)in[i];
   return out;
}

const GLfloat *f2GLf(double *in, int size){
   int i;
   GLfloat *out;

   //out = (GLfloat *)malloc(size*sizeof(GLfloat));
   out = new GLfloat [size];
   for(i=0;i<size;i++)
     out[i] = (GLfloat)in[i];
   return out;
}

const GLfloat *buildrotmatrix_from_c(matrix a) {
  GLfloat *f;

  f = new GLfloat [16];
  f[0]  = a(0,0);   f[1]  = a(0,1);  f[2]  = a(0,2);  f[3]  = a(0,3);
  f[4]  = a(1,0);   f[5]  = a(1,1);  f[6]  = a(1,2);  f[7]  = a(1,3);
  f[8]  = a(2,0);   f[9]  = a(2,1);  f[10] = a(2,2);  f[11] = a(2,3);
  f[12] = a(3,0);   f[13] = a(3,1);  f[14] = a(3,2);  f[15] = a(3,3);

  return f;
}

const GLfloat *buildrotmatrix(
GLfloat a0, GLfloat a1, GLfloat a2, GLfloat a3, 
GLfloat a4, GLfloat a5, GLfloat a6, GLfloat a7, 
GLfloat a8, GLfloat a9, GLfloat a10, GLfloat a11, 
GLfloat a12, GLfloat a13, GLfloat a14, GLfloat a15) {
  GLfloat *f;

  f = new GLfloat [16];
  f[0]  = a0;   f[1] = a1;  f[2]  =  a2; f[3]  =  a3;
  f[4]  = a4;   f[5] = a5;  f[6]  =  a6; f[7] =   a7;
  f[8]  = a8;   f[9] = a9;  f[10] = a10; f[11] = a11;
  f[12] = a12; f[13] = a13; f[14] = a14; f[15] = a15;

  return f;
}

Volume GetFrontAndBackClippingPlanes(){

  GLdouble projMatrix[16];
  GLdouble modelMatrix[16];

  glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);

  matrix MM(4,4,modelMatrix);
  matrix PM(4,4,projMatrix);

  matrix ClipMatrix = MM * PM;

  Plane farclip   ((ClipMatrix(0,3)-ClipMatrix(0,2)),(ClipMatrix(1,3)-ClipMatrix(1,2)),(ClipMatrix(2,3)-ClipMatrix(2,2)),(ClipMatrix(3,3)-ClipMatrix(3,2)));
  Plane nearclip  ((ClipMatrix(0,3)+ClipMatrix(0,2)),(ClipMatrix(1,3)+ClipMatrix(1,2)),(ClipMatrix(2,3)+ClipMatrix(2,2)),(ClipMatrix(3,3)+ClipMatrix(3,2)));

  nearclip.Normalize();
  farclip.Normalize();

  Volume v;

  v.AddPlane(farclip);
  v.AddPlane(nearclip);

  return v;
}

Volume GetClippingPlanes(){

  GLdouble projMatrix[16];
  GLdouble modelMatrix[16];

  glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);

  matrix MM(4,4,modelMatrix);
  matrix PM(4,4,projMatrix);

  matrix ClipMatrix = MM * PM;

  Plane right ((ClipMatrix(0,3)-ClipMatrix(0,0)),(ClipMatrix(1,3)-ClipMatrix(1,0)),(ClipMatrix(2,3)-ClipMatrix(2,0)),(ClipMatrix(3,3)-ClipMatrix(3,0)));
  Plane left  ((ClipMatrix(0,3)+ClipMatrix(0,0)),(ClipMatrix(1,3)+ClipMatrix(1,0)),(ClipMatrix(2,3)+ClipMatrix(2,0)),(ClipMatrix(3,3)+ClipMatrix(3,0)));
  Plane bottom((ClipMatrix(0,3)+ClipMatrix(0,1)),(ClipMatrix(1,3)+ClipMatrix(1,1)),(ClipMatrix(2,3)+ClipMatrix(2,1)),(ClipMatrix(3,3)+ClipMatrix(3,1)));
  Plane top   ((ClipMatrix(0,3)-ClipMatrix(0,1)),(ClipMatrix(1,3)-ClipMatrix(1,1)),(ClipMatrix(2,3)-ClipMatrix(2,1)),(ClipMatrix(3,3)-ClipMatrix(3,1)));
  Plane farclip   ((ClipMatrix(0,3)-ClipMatrix(0,2)),(ClipMatrix(1,3)-ClipMatrix(1,2)),(ClipMatrix(2,3)-ClipMatrix(2,2)),(ClipMatrix(3,3)-ClipMatrix(3,2)));
  Plane nearclip  ((ClipMatrix(0,3)+ClipMatrix(0,2)),(ClipMatrix(1,3)+ClipMatrix(1,2)),(ClipMatrix(2,3)+ClipMatrix(2,2)),(ClipMatrix(3,3)+ClipMatrix(3,2)));

  left.Normalize();
  right.Normalize();
  top.Normalize();
  bottom.Normalize();
  nearclip.Normalize();
  farclip.Normalize();


  GLint param0,param1;
  glGetIntegerv(GL_CLIP_PLANE0,&param0);
  glGetIntegerv(GL_CLIP_PLANE1,&param1);
  GLdouble eqn0[4], eqn1[4];
  if(param0&&param1){
    glGetClipPlane(GL_CLIP_PLANE0,eqn0);
    glGetClipPlane(GL_CLIP_PLANE1,eqn1);
    double scale = (eqn0[3] + eqn1[3])/(nearclip.get_D() + farclip.get_D());
    double diff = 0.5*(nearclip.get_D()-farclip.get_D());
    farclip.set_D(farclip.get_D()*scale-diff);
    nearclip.set_D(nearclip.get_D()*scale+diff);
  }

  Volume v;

  v.AddPlane(right);
  v.AddPlane(left);
  v.AddPlane(bottom);
  v.AddPlane(top);
  v.AddPlane(farclip);
  v.AddPlane(nearclip);

  return v;

}

bool isPointInClippingVolume(const Cartesian &p, const Volume &v){
  return v.PointInVolume(p);
}
