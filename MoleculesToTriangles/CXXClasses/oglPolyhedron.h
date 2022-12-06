/*
 *  oglPolyhedron.h
 *  iPhoneRibbons
 *
 *  Created by Martin Noble on 19/09/2008.
 *  Copyright 2008 LMB, Oxford University. All rights reserved.
 *
 */

#ifndef oglPolyhedron_H
#define oglPolyhedron_H
#include <memory>

#include "Polyhedron.h"
#include "VertexColorNormalPrimitive.h"

class Renderer;

//unsigned short for OpenGL ES

class oglPolyhedron : public Polyhedron, public VertexColorNormalPrimitive {
private:
    //Temporary varaiables to hold instantaneous values
    FCXXCoord centre;
    FCXXCoord color;
    float radius;
public:
	static std::map<Polyhedron::PolyhedronType, oglPolyhedron *> globalPolyhedra;
    oglPolyhedron();
    virtual ~oglPolyhedron(){
    };
    const float getRadius() const {
        return radius;
    };
    void setRadius(const float &_radius){
        radius = _radius;
    };
    const FCXXCoord &getCentre() const {
        return centre;
    };
    void setCentre(const FCXXCoord &_centre) {
        centre = _centre;
    };
    const FCXXCoord &getColor() const {
        return color;
    };
    void setColor (const FCXXCoord &_color){
        color = _color;
    };
    virtual void generateArrays();
    virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer);    
	static oglPolyhedron *octaHedron ();
	static oglPolyhedron *dodecaHedron ();
	static oglPolyhedron *dodecaHedron4 ();
	static oglPolyhedron *dodecaHedron16 ();
};


#endif
