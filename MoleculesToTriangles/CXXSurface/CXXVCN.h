//
//  CXXVCN.h
//  MoleculesToTriangles
//
//  Created by Martin Noble on 9/8/16.
//  Copyright © 2016 MartinNobleSoftware. All rights reserved.
//

#ifndef CXXVCN_h
#define CXXVCN_h
#include "CXXCoord.h"
#include <stdio.h>
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

class CXXVCN {
public:
    CXXCoord<GLfloat>vertex;
    CXXCoord<GLfloat>color;
    CXXCoord<GLfloat>normal;
};


#endif /* CXXVCN_h */
