/*
 *  DisplayPrimitive.cpp
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 18/07/2008.
 *  Copyright 2008 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
 */

#include "DisplayPrimitive.h"
#include "Renderer.h"
#include <iostream>

DisplayPrimitive::~DisplayPrimitive()
{
    liberateAllHandles();
    //std::cout << "In displayPrimitive destructor" << primitiveType;
}

void DisplayPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
    std::cout << "Base class renderWithRenderer" << std::endl;
}

void DisplayPrimitive::liberateAllHandles()
{
    //std::cout << "In liberateAllHandles\n";
    auto renderersIter = renderers.begin();
    for (; renderersIter != renderers.end(); ++renderersIter){
        //std::cout << "Identified a renderer to disarm\n";
        (*renderersIter)->liberateHandlesForDisplayPrimitive(this);
    }
}
