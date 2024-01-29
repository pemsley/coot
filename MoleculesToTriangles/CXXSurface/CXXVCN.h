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
// 20240120-PE Goodby openGL things.
// #ifdef __APPLE_CC__
// #include <OpenGL/gl.h>
// #else
// #include <GL/gl.h>
// #endif

class CXXVCN {
public:
    CXXCoord<float>vertex;
    CXXCoord<float>color;
    CXXCoord<float>normal;
};


#endif /* CXXVCN_h */
