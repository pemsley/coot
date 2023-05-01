//
//  BallsPrimitive.h
//  AesopCD_ios
//
//  Created by Martin Noble on 09/02/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
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
