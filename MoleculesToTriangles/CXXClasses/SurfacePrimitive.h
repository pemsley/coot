/*
 *  SurfacePrimitive.h
 *  AesopCoreData
 *
 *  Created by Martin Noble on 23/03/2010.
 *  Copyright 2010 LMB, Oxford University. All rights reserved.
 *
 */
#ifndef SurfacePrimitive_h
#define SurfacePrimitive_h
#include <memory>

#include "VertexColorNormalPrimitive.h"
#include "ColorScheme.h"

#include "MoleculesToTriangles/CXXSurface/CXXSurfaceMaker.h"

class SurfacePrimitive : public VertexColorNormalPrimitive {
private:
	CXXSurfaceMaker *cxxSurfaceMaker;
    std::shared_ptr<ColorScheme> colorScheme;
public:
    enum SurfaceType {AccessibleSurface, VdWSurface, MolecularSurface};
    SurfacePrimitive();
	SurfacePrimitive(mmdb::Manager *mmdb, int chunkHndl, int selHnd, std::shared_ptr<ColorScheme> _colorScheme, enum SurfaceType type, float probeRadius, float radiusMultiplier);
    virtual ~SurfacePrimitive(){
        if (cxxSurfaceMaker) delete cxxSurfaceMaker;
        //std::cout << "In surface destructor" << std::endl;
    };
    
    virtual void generateArrays();

    CXXSurfaceMaker *getCXXSurfaceMaker(){
        return cxxSurfaceMaker;
    };
};

#endif
