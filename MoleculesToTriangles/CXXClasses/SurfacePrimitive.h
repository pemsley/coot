/*
 * MoleculesToTriangles/CXXClasses/SurfacePrimitive.h
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
