/*
     pygl/texture.cc: CCP4MG Molecular Graphics Program
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

#ifdef _WIN32
#include <windows.h>
#endif

#define GLUT_DISABLE_ATEXIT_HACK
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "ppmutil.h"
#include "texture.h"

void set_texture_coord(GLfloat x, GLfloat y, int textured){
#if defined(GL_VERSION_1_2) && defined (linux)
  if(textured>1){
// BL says:: take out the glMultiTexCoord2f calls as they fail on FC4 for now
// shall fix preoperly at some point (ask SMc). FIXME
//
//#if defined (GL_VERSION_1_3)
//      if(textured&1) glMultiTexCoord2f(GL_TEXTURE0_ARB,x,y);
//      if(textured&2) glMultiTexCoord2f(GL_TEXTURE1_ARB,x,y);
//#else
//      if(textured&1) glMultiTexCoord2fARB(GL_TEXTURE0_ARB,x,y);
//      if(textured&2) glMultiTexCoord2fARB(GL_TEXTURE1_ARB,x,y);
//#endif
    glTexCoord2f(x,y);
  }else{
    glTexCoord2f(x,y);
  }
#else
      glTexCoord2f(x,y);
#endif
}

unsigned int load_texture(const image_info &iinfo, int style)
{

  int width  = iinfo.get_width();
  int height = iinfo.get_height();
  image_info tiinfo(iinfo.get_width(),iinfo.get_height(),iinfo.get_pixels(),iinfo.get_colourspace_type());

  int max_width = 256;
  int max_height = 256;

  GLint maxwidthheight;

  glGetIntegerv(GL_MAX_TEXTURE_SIZE,&maxwidthheight);

  max_width = max_height = maxwidthheight;

  if(width<64)
    width = 64;
  if(width>64&&width<128)
    width = 128;
  if(width>128&&width<256)
    width = 256;
  if(width>256&&width<512)
    width = 512;
  if(width>512&&width<1024)
    width = 1024;
  if(width>1024&&width<2048)
    width = 2048;

  if(width>max_width) width = max_width;

  if(height<64)
     height= 64;
  if(height>64&&height<128)
     height= 128;
  if(height>128&&height<256)
     height= 256;
  if(height>256&&height<512)
    height = 512;
  if(height>512&&height<1024)
    height = 1024;
  if(height>1024&&height<2048)
    width = 2048;

  if(height>max_height) height = max_height;

  tiinfo.ScaleImage(width,height);

  width  = tiinfo.get_width();
  height = tiinfo.get_height();

  GLuint texName;
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glGenTextures(1, &texName);
  glBindTexture(GL_TEXTURE_2D, texName);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

  int glcolourspace;
  if(iinfo.get_colourspace_type()==IMAGEINFO_RGBA){
    glcolourspace = GL_RGBA;
  }else if(iinfo.get_colourspace_type()==IMAGEINFO_MONOA){
    tiinfo.convert_rgba();
    glcolourspace = GL_RGBA;
  }else if(iinfo.get_colourspace_type()==IMAGEINFO_RGB){
    glcolourspace = GL_RGB;
  }else{
    tiinfo.convert_rgb();
    glcolourspace = GL_RGB;
  }

  unsigned char* pixels = tiinfo.get_pixels();

  if(style==NEAREST){
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 
	       0, (GLenum)glcolourspace, GL_UNSIGNED_BYTE, (GLubyte*)pixels);
  }else if(style==LINEAR){
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 
	       0, (GLenum)glcolourspace, GL_UNSIGNED_BYTE, (GLubyte*)pixels);
  }else if(style==MIPMAP){
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, width, height, (GLenum)glcolourspace, GL_UNSIGNED_BYTE, (GLubyte*)pixels);
  }

  return texName;
}
