//
//  LinesPrimitive.cpp
//  MoleculesToTriangles
//
//  Created by Martin Noble on 1/31/17.
//  Copyright © 2017 MartinNobleSoftware. All rights reserved.
//

#include <stdio.h>
#include "LinesPrimitive.h"

void LinesPrimitive::setStartPoint(FCXXCoord &startPoint, FCXXCoord &startColor) {
    vertices.clear();
    colors.clear();
    invalidateGLPrimitives();
    vertices.push_back(startPoint);
    colors.push_back(startColor);
};

void LinesPrimitive::invalidateGLPrimitives()
{
    delete [] vertexArray;
    delete [] colorArray;
    colorArray = 0;
    vertexArray = 0;
}
