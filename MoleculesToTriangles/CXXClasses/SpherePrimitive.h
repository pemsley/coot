/*
 * MoleculesToTriangles/CXXClasses/SpherePrimitive.h
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
