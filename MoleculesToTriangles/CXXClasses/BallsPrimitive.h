/*
 * MoleculesToTriangles/CXXClasses/BallsPrimitive.h
 *
 * Copyright 2011 by Martin Noble, University of Oxford
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
//  BallsPrimitive.h
//  AesopCD_ios
//

#ifndef BallsPrimitive_H
#define BallsPrimitive_H
#include <vector>
#include <memory>

#include "VertexColorNormalPrimitive.h"
typedef struct _Ball {
    FCXXCoord centre;
    FCXXCoord color;
    float radius;
    FCXXCoord normal;
    float radiusAlongNormal;
} Ball;

class BallsPrimitive : public VertexColorNormalPrimitive {
private:
    std::vector<Ball> balls;
public:
    int polyhedronType;
    BallsPrimitive() {
        vertexColorNormalArray = 0;
        indexArray = 0;
        emptyArrays();
        drawModeGL = DrawAsTriangles;
        enableColorGL = true;
        primitiveType = DisplayPrimitive::PrimitiveType::BallsPrimitive;
    };
    void addBall (FCXXCoord &centre, FCXXCoord &color, float radius, FCXXCoord &normal, float radiusAlongNormal){
        Ball ball= {centre, color, radius, normal, radiusAlongNormal};
        addBall (ball);
    };
    void addBall (FCXXCoord &centre, FCXXCoord &color, float radius){
        FCXXCoord normal(0.,0.,1.);
        float radiusAlongNormal = radius;
        Ball ball= {centre, color, radius, normal, radiusAlongNormal};
        addBall (ball);
    };
    void addBall(Ball &ball) {
        balls.push_back(ball);
    };
	void emptyArrays(){
        delete [] vertexColorNormalArray;
        vertexColorNormalArray = 0;
        delete [] indexArray;
        indexArray = 0;
    };
    virtual ~BallsPrimitive(){
        emptyArrays();
        //std::cout << "In balls destructor" << std::endl;
    };
    virtual void generateArrays();
    const std::vector<Ball> &getBalls(){
        return balls;
    }
    
};

#endif
