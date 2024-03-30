/*
 * MoleculesToTriangles/CXXClasses/Representation.cpp
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This file is part of Coot
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

#if defined _OPENMP
#include <omp.h>
#else
#if __APPLE__
#include <dispatch/dispatch.h>
#endif
#endif

#include "Representation.h"
#include "DisplayPrimitive.h"

#include "SpherePrimitive.h"
#include "Callbacks.h"

#if __APPLEXX__
#define USINGGCD
#endif

MyCompletionPtrType Representation::redrawCompletionCallback = 0;

MyProgressPtrType Representation::redrawProgressCallback = 0;

void Representation::renderWithRenderer(std::shared_ptr<Renderer> renderer){
#ifdef DEBUG_MINE
    std::cout << "In representation render " << redrawNeeded <<" " << getDoDraw()<< "\n";
#endif
    if (redrawNeeded && !inRedraw) {
        inRedraw = true;
#ifdef USINGGCD
        dispatch_async(dispatch_get_main_queue() , ^{
#endif
            deletePrimitives();
            redraw();
            redrawNeeded = false;
            inRedraw = false;
            if (redrawCompletionCallback){
                redrawCompletionCallback(redrawCompletionCallbackUserInfo);
            }
#ifdef USINGGCD
        });
#endif
    }
    if (getDoDraw()){
        std::vector<std::shared_ptr<DisplayPrimitive> >::iterator primitivePntr = displayPrimitives.begin();
        while (primitivePntr != displayPrimitives.end()){
            DisplayPrimitive &dp = **primitivePntr;
            dp.renderWithRenderer(renderer);
            primitivePntr++;
        }
    }
    return;
}

void Representation::generateArrays(){
    std::vector<std::shared_ptr<DisplayPrimitive> >::iterator primitivePntr = displayPrimitives.begin();
    while (primitivePntr != displayPrimitives.end()){
        (*primitivePntr)->generateArrays();
        primitivePntr++;
    }    
}

void Representation::deletePrimitives(){
    displayPrimitives.clear();
}

