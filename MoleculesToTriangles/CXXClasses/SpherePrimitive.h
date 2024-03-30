/*
 *  SpherePrimitive.h
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 18/07/2008.
 *  Copyright 2008 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
 */

#ifndef SpherePrimitive_h
#define SpherePrimitive_h
#include <memory>

#include "DisplayPrimitive.h"
#include "oglPolyhedron.h"

class SpherePrimitive : public DisplayPrimitive {
private:
    FCXXCoord _centre;
    FCXXCoord _color;
    float _radius;
	
	oglPolyhedron *polyhedron;
    
public:
    SpherePrimitive();
	SpherePrimitive(const FCXXCoord &newCentre, const float newRadius, const FCXXCoord &newColor){
		_centre = newCentre;
		_radius = newRadius;
		_color = newColor;
        primitiveType = DisplayPrimitive::PrimitiveType::SpherePrimitive;
	}
    virtual ~SpherePrimitive(){
        polyhedron = 0;
        //std::cout << "In sphere destructor" << std::endl;
    }
    float radius();
    FCXXCoord &color();
    FCXXCoord &centre();
    void setRadius(const float &radius) {
        _radius = radius;
    };
    
    virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer);
    virtual void generateArrays();
        
};

#endif
