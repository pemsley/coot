/*
 *  Representation.h
 *  Aesop
 *
 *  Created by Martin Noble on 07/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
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
