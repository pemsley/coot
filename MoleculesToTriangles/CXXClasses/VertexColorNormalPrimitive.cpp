/*
 *  VertexColorNormalPrimitive.cpp
 *  AesopCoreData
 *
 *  Created by Martin Noble on 24/03/2010.
 *  Copyright 2010 Apple Inc. All rights reserved.
 *
 */

#include "VertexColorNormalPrimitive.h"
#include "Renderer.h"

void VertexColorNormalPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)  {    
    if (vertexColorNormalArray == 0){
        generateArrays();
    }
    renderer->renderVertexColorNormalPrimitive(this);
};
