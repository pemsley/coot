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
    /**
     * The second atom of a saddle, and how much of the vertex belongs to the FIRST one - the
     * one in the base class's atomArray. Both are per vertex and both are filled by
     * generateArrays.
     *
     * A convex cap has only one atom, and its entry here is null with a weight of 1. A vertex
     * on a torus lies in the groove between two atoms and genuinely belongs to both, in the
     * proportion the torus parametrisation gives - see CXXTorusElement::weightOfNodeAtom.
     *
     * These are on SurfacePrimitive rather than beside atomArray on the base class because
     * only a surface has anything to say here; a ribbon vertex belongs to one residue and a
     * bond vertex to one atom.
     */
    const mmdb::Atom **atom2Array;
    float *atomWeightArray;
public:
    enum SurfaceType {AccessibleSurface, VdWSurface, MolecularSurface};
    SurfacePrimitive();
	SurfacePrimitive(mmdb::Manager *mmdb, int chunkHndl, int selHnd, std::shared_ptr<ColorScheme> _colorScheme, enum SurfaceType type, float probeRadius, float radiusMultiplier);
    virtual ~SurfacePrimitive(){
        if (cxxSurfaceMaker) delete cxxSurfaceMaker;
        delete [] atom2Array;
        atom2Array = 0;
        delete [] atomWeightArray;
        atomWeightArray = 0;
        //std::cout << "In surface destructor" << std::endl;
    };

    virtual void generateArrays();

    /** The other atom of the saddle each vertex lies on, or null where there is only one. */
    const mmdb::Atom **getAtom2Array() const {
        return atom2Array;
    };
    /** How much of each vertex belongs to getAtomArray()'s atom rather than getAtom2Array()'s. */
    const float *getAtomWeightArray() const {
        return atomWeightArray;
    };

    CXXSurfaceMaker *getCXXSurfaceMaker(){
        return cxxSurfaceMaker;
    };
};

#endif
