/*
 * MoleculesToTriangles/CXXClasses/Callbacks.h
 *
 * This file is part of Aesop
 *
 * Copyright 2012 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

//
//  Callbacks.h
//  Aesop
//

#ifndef Aesop_Callbacks_h
#define Aesop_Callbacks_h

typedef void (*MyCompletionPtrType)(void *);
typedef void (*MyProgressPtrType)(void *, float);

#ifdef __cplusplus
extern "C" {
#endif
    void globalRedrawCompletionCallback(void *userInfo);   
    void globalRedrawProgressCallback(void *userInfo, float progress);   
    
#ifdef __cplusplus
}
#endif

#endif
