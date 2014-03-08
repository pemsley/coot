/*
     pygl/texture.h: CCP4MG Molecular Graphics Program
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


#ifndef __TEXTURE_HELP__
#define __TEXTURE_HELP__

#define __GL_HELP__
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <utility>
#include "ppmutil.h"
enum {NEAREST, LINEAR, MIPMAP};
#if defined (_WIN32) && not defined (WINDOWS_MINGW)
#define TEXTURE_DLL __declspec(dllexport)
unsigned int __stdcall TEXTURE_DLL load_texture(const image_info &iinfo, int style);
std::pair<unsigned,unsigned> __stdcall TEXTURE_DLL GetCompatibleTextureSize(const unsigned width_in, const unsigned height_in);
image_info __stdcall TEXTURE_DLL ResizeWithEmptySpace(const image_info &iinfo, const unsigned width, const unsigned height);
void __stdcall TEXTURE_DLL set_texture_coord(GLfloat x, GLfloat y, int textured);
#else
unsigned int load_texture(const image_info &iinfo, int style);
std::pair<unsigned,unsigned> GetCompatibleTextureSize(const unsigned width_in, const unsigned height_in);
image_info ResizeWithEmptySpace(const image_info &iinfo, const unsigned width, const unsigned height);
void set_texture_coord(GLfloat x, GLfloat y, int textured);
#endif

#endif
