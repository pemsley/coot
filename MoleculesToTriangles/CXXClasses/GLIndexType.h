//
//  GLIndexType.h
//  AesopCD
//
//  Created by Martin Noble on 28/09/2010.
//  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
//

#ifndef GLIndexType_H
#define GLIndexType_H
#ifdef __APPLE__
#include "TargetConditionals.h"
#endif
#if TARGET_OS_IPHONE
#define kGLIndexType GL_UNSIGNED_SHORT
typedef unsigned short GLIndexType;
#elif TARGET_IPHONE_SIMULATOR
#define kGLIndexType GL_UNSIGNED_SHORT
typedef unsigned short GLIndexType;
#else
#define kGLIndexType GL_UNSIGNED_INT
typedef unsigned int GLIndexType;
#endif

#endif
