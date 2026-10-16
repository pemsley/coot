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
        BallsPrimitive,
        FlatFanPrimitive    // appended, so the values above keep their numbering
    };
    // Initialised, because type() is used to decide what a primitive may be cast to.
    //
    // This was a bare member, and four classes - BondsPrimitive, LinesPrimitive,
    // MMDBStringPrimitive and FlatFanPrimitive - never assigned it, so type() returned whatever
    // was on the heap. A caller that admits a set of types and then casts on the strength of
    // the answer would sooner or later admit a primitive that is not of the class it claims:
    // the M2T mesh builder does exactly that, and threw std::bad_cast whenever the garbage in
    // a BondsPrimitive happened to read as Balls, Cylinder, BoxSection or Surface. It looked
    // like running out of memory, because which garbage turns up depends on what the heap has
    // been doing.
    //
    // LinePrimitive as the default so that a class which still forgets is left OUT of a mesh of
    // triangles rather than cast to something it is not.
    PrimitiveType primitiveType = LinePrimitive;
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

