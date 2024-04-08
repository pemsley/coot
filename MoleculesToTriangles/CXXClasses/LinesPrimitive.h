/*
 * MoleculesToTriangles/CXXClasses/LinesPrimitive.h
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

#ifndef LinesPrimitive_h
#define LinesPrimitive_h

#include "DisplayPrimitive.h"

class LinesPrimitive : public DisplayPrimitive {
private:
protected:
    std::vector<FCXXCoord  > vertices;
    std::vector<FCXXCoord  > colors;
    float *vertexArray;
    float *colorArray;
    void evaluateGLPrimitives();
    void invalidateGLPrimitives();
    int nBonds;
public:
    LinesPrimitive() : vertexArray(0), colorArray(0) {};
    LinesPrimitive(FCXXCoord &startPoint, FCXXCoord &startColor){
        vertexArray = colorArray = 0;
        setStartPoint(startPoint, startColor);
    };
    virtual ~LinesPrimitive(){
        //std::cout << "In lines destructor" << std::endl;
        invalidateGLPrimitives();
        vertices.clear();
        colors.clear();
    };
    void setStartPoint(FCXXCoord &startPoint, FCXXCoord &startColor);
    void addVertex(FCXXCoord &point, FCXXCoord &color) {
        invalidateGLPrimitives();
        vertices.push_back(point);
        colors.push_back(color);
    };
    const float *getVertexArray() const {
        return vertexArray;
    };
    const float *getColorArray() const {
        return colorArray;
    };
    unsigned long nVertices() const {
        return vertices.size();
    };
    int getNBonds() const {
        return nBonds;
    };
    virtual void renderWithRenderer(Renderer *renderer);
};

#endif
