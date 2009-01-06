#ifndef _CCP4MG_SPHERE_0_
#define _CCP4MG_SPHERE_0_
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#define X_0_0   0.525731112119
#define X_0_1   0.000000000000
#define X_0_2   0.850650808352

static GLfloat icosa_0[36] = {
  -X_0_0, X_0_2, X_0_1,
  X_0_1, X_0_0, X_0_2,
  X_0_0, X_0_2, X_0_1,
  -X_0_2, X_0_1, X_0_0,
  X_0_1, -X_0_0, X_0_2,
  X_0_2, X_0_1, X_0_0,
  X_0_2, X_0_1, -X_0_0,
  X_0_0, -X_0_2, X_0_1,
  -X_0_0, -X_0_2, X_0_1,
  X_0_1, -X_0_0, -X_0_2,
  X_0_1, X_0_0, -X_0_2,
  -X_0_2, X_0_1, -X_0_0
};

static GLuint tindices_0[] = {
  0, 1, 2,
  0, 3, 1,
  3, 4, 1,
  1, 4, 5,
  1, 5, 2,
  5, 6, 2,
  5, 7, 6,
  4, 7, 5,
  4, 8, 7,
  8, 9, 7,
  9, 6, 7,
  9, 10, 6,
  9, 11, 10,
  11, 0, 10,
  0, 2, 10,
  10, 2, 6,
  3, 0, 11,
  3, 11, 8,
  3, 8, 4,
  9, 8, 11
};

#endif /*_CCP4MG_SPHERE_1_*/
