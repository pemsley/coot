/*
 *  LinesPrimitive.h
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 19/07/2008.
 *  Copyright 2008 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
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
