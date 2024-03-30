/*
 * MoleculesToTriangles/CXXClasses/Representation.h
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

#ifndef Representation_H
#define Representation_H

#include <vector>
#include <list>
#include <map>
#include <string>
#include <iostream>
#include <memory>

#include "Callbacks.h"
#include "DisplayPrimitive.h"

class Renderer;
//class DisplayPrimitive;
class RepresentationInstance;

class Representation {
protected:
    std::vector<std::shared_ptr<DisplayPrimitive> > displayPrimitives;
    std::list<RepresentationInstance *>representationInstances;     
    std::map<std::string, float> floatParameters;
    std::map<std::string, float> intParameters;
    std::map<std::string, bool> boolParameters;
    bool redrawNeeded;
    bool inRedraw;
    static MyCompletionPtrType redrawCompletionCallback;
    void * redrawCompletionCallbackUserInfo;
    static MyProgressPtrType redrawProgressCallback;
    void * redrawProgressCallbackUserInfo;
public:
    Representation() {
        boolParameters["doDraw"] = true;
        inRedraw = false;
    };
    virtual ~Representation(){
        deletePrimitives();
    }; 
    void updateFloatParameter(std::string parameterName, float parameterValue){
        floatParameters[parameterName] = parameterValue;  
        redrawNeeded = true;
    };
    void updateIntParameter(std::string parameterName, int parameterValue){
        intParameters[parameterName] = parameterValue;  
        redrawNeeded = true;
    };
    void updateBoolParameter(std::string parameterName, bool parameterValue){
        boolParameters[parameterName] = parameterValue;          
        if (parameterName != std::string("doDraw")) redrawNeeded = true;
    };
    void renderWithRenderer(std::shared_ptr<Renderer> renderer);
    void generateArrays();
    void deletePrimitives();
    virtual bool getDoDraw() = 0;
    void setDoDraw(const bool &yesOrNo){
        boolParameters["doDraw"] = yesOrNo;
    };
    void addPrimitive(std::shared_ptr<DisplayPrimitive>prim){
        displayPrimitives.push_back(prim);
    };
    void setRedrawNeeded(const bool &yesOrNo) {
        redrawNeeded = yesOrNo;  
    };
    void setRedrawCompletionCallback(MyCompletionPtrType aCallback, void * someUserInfo) {
        redrawCompletionCallback = aCallback;
        redrawCompletionCallbackUserInfo = someUserInfo;
    };
    void setRedrawProgressCallback(MyProgressPtrType aCallback, void * someUserInfo) {
        redrawProgressCallback = aCallback;
        redrawProgressCallbackUserInfo = someUserInfo;
    };
    unsigned long nDisplayPrimitives(){
        return displayPrimitives.size();
    };
    std::vector<std::shared_ptr<DisplayPrimitive> > &getDisplayPrimitives(){
        return displayPrimitives;
    }
    virtual void redraw() = 0;
};

#endif
