/*
 *  Representation.mm
 *  Aesop
 *
 *  Created by Martin Noble on 07/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
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

