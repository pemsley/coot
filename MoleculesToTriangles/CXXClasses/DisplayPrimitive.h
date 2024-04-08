/*
 * MoleculesToTriangles/CXXClasses/DisplayPrimitive.h
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

#ifndef DisplayPrimitive_h
#define DisplayPrimitive_h

#include <vector>
#include <memory>
#include <map>
#include <set>
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
class Renderer;

class DisplayPrimitive {
protected:
private:
    // The point here is that a primitive will have resources allocated in different renderers. When
    // it gets deleted, those resource have to be freed up, so it needs to keep a list of those
    // renderers that hold resources for it.  I can't use shared_ptr for this, since the retention of those
    // shared pointers will prvent the renderer from being deallocated
    std::set<Renderer* > renderers;
public:
    enum PrimitiveType {
        SpherePrimitive,
        EllipsePrimitive,
        CylinderPrimitive,
        LinePrimitive,
        BoxSectionPrimitive,
        SurfacePrimitive,
        BallsPrimitive
    };
    PrimitiveType primitiveType;
    PrimitiveType type() {
        return primitiveType;
    };
    virtual ~DisplayPrimitive();
    void addRenderer(Renderer* aRenderer){
        renderers.insert(aRenderer);
    };
    void removeRenderer(Renderer* aRenderer){
        renderers.erase(aRenderer);
    };
    //{
    //    std::cout << "DisplayPrimitive destructor " << std::endl;
    //};
    virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer) = 0;
    
    //Children of this class will have to implement the following
    virtual void generateArrays() = 0;
    void liberateAllHandles();
    
};

#endif

