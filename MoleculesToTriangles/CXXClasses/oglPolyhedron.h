/*
 * MoleculesToTriangles/CXXClasses/oglPolyhedron.h
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
