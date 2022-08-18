/*
 *  DisplayPrimitive.h
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 18/07/2008.
 *  Copyright 2008 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
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

