/*
     pygl/sphere_1.h: CCP4MG Molecular Graphics Program
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
#ifndef _CCP4MG_SPHERE_1_
#define _CCP4MG_SPHERE_1_
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#define X_1_0   0.525731112119
#define X_1_1   0.000000000000
#define X_1_2   0.850650808352
#define X_1_3   0.309016994375
#define X_1_4   0.500000000000
#define X_1_5   0.809016994375
#define X_1_6   1.000000000000

static GLfloat icosa_1[126] = {
  -X_1_0, X_1_2, X_1_1,
  -X_1_3, X_1_5, X_1_4,
  X_1_1, X_1_6, X_1_1,
  X_1_1, X_1_0, X_1_2,
  X_1_3, X_1_5, X_1_4,
  X_1_0, X_1_2, X_1_1,
  -X_1_5, X_1_4, X_1_3,
  -X_1_2, X_1_1, X_1_0,
  -X_1_4, X_1_3, X_1_5,
  -X_1_4, -X_1_3, X_1_5,
  X_1_1, -X_1_0, X_1_2,
  X_1_1, X_1_1, X_1_6,
  X_1_4, X_1_3, X_1_5,
  X_1_4, -X_1_3, X_1_5,
  X_1_2, X_1_1, X_1_0,
  X_1_5, X_1_4, X_1_3,
  X_1_6, X_1_1, X_1_1,
  X_1_2, X_1_1, -X_1_0,
  X_1_5, X_1_4, -X_1_3,
  X_1_5, -X_1_4, X_1_3,
  X_1_0, -X_1_2, X_1_1,
  X_1_5, -X_1_4, -X_1_3,
  X_1_3, -X_1_5, X_1_4,
  -X_1_3, -X_1_5, X_1_4,
  -X_1_0, -X_1_2, X_1_1,
  X_1_1, -X_1_6, X_1_1,
  -X_1_3, -X_1_5, -X_1_4,
  X_1_1, -X_1_0, -X_1_2,
  X_1_3, -X_1_5, -X_1_4,
  X_1_4, -X_1_3, -X_1_5,
  X_1_1, X_1_1, -X_1_6,
  X_1_1, X_1_0, -X_1_2,
  X_1_4, X_1_3, -X_1_5,
  -X_1_4, -X_1_3, -X_1_5,
  -X_1_2, X_1_1, -X_1_0,
  -X_1_4, X_1_3, -X_1_5,
  -X_1_5, X_1_4, -X_1_3,
  -X_1_3, X_1_5, -X_1_4,
  X_1_3, X_1_5, -X_1_4,
  -X_1_6, X_1_1, X_1_1,
  -X_1_5, -X_1_4, X_1_3,
  -X_1_5, -X_1_4, -X_1_3
};

static GLuint tindices_1[] = {
  0, 1, 2,
  3, 4, 1,
  5, 2, 4,
  1, 4, 2,
  0, 6, 1,
  7, 8, 6,
  3, 1, 8,
  6, 8, 1,
  7, 9, 8,
  10, 11, 9,
  3, 8, 11,
  9, 11, 8,
  3, 11, 12,
  10, 13, 11,
  14, 12, 13,
  11, 13, 12,
  3, 12, 4,
  14, 15, 12,
  5, 4, 15,
  12, 15, 4,
  14, 16, 15,
  17, 18, 16,
  5, 15, 18,
  16, 18, 15,
  14, 19, 16,
  20, 21, 19,
  17, 16, 21,
  19, 21, 16,
  10, 22, 13,
  20, 19, 22,
  14, 13, 19,
  22, 19, 13,
  10, 23, 22,
  24, 25, 23,
  20, 22, 25,
  23, 25, 22,
  24, 26, 25,
  27, 28, 26,
  20, 25, 28,
  26, 28, 25,
  27, 29, 28,
  17, 21, 29,
  20, 28, 21,
  29, 21, 28,
  27, 30, 29,
  31, 32, 30,
  17, 29, 32,
  30, 32, 29,
  27, 33, 30,
  34, 35, 33,
  31, 30, 35,
  33, 35, 30,
  34, 36, 35,
  0, 37, 36,
  31, 35, 37,
  36, 37, 35,
  0, 2, 37,
  5, 38, 2,
  31, 37, 38,
  2, 38, 37,
  31, 38, 32,
  5, 18, 38,
  17, 32, 18,
  38, 18, 32,
  7, 6, 39,
  0, 36, 6,
  34, 39, 36,
  6, 36, 39,
  7, 39, 40,
  34, 41, 39,
  24, 40, 41,
  39, 41, 40,
  7, 40, 9,
  24, 23, 40,
  10, 9, 23,
  40, 23, 9,
  27, 26, 33,
  24, 41, 26,
  34, 33, 41,
  26, 41, 33
};

#endif /*_CCP4MG_SPHERE_1_*/
