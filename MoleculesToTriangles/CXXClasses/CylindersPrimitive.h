/*
 * MoleculesToTriangles/CXXClasses/CylindersPrimitive.h
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

#ifndef CylindersPrimitive_h
#define CylindersPrimitive_h
#include <memory>

#include "VertexColorNormalPrimitive.h"
#include "CylinderPoint.h"

class CylindersPrimitive : public VertexColorNormalPrimitive {
private:
    friend class BoxSectionPrimitive;
    std::vector<CylinderPoint>points;
    int angularSampling;
public:
    CylindersPrimitive();
    CylindersPrimitive(FCXXCoord &startPoint, FCXXCoord &startColor, FCXXCoord &normalOne, FCXXCoord &normalTwo,
                       float radiusOne, float radiusTwo){
        angularSampling = 12;
        vertexColorNormalArray = 0;
        indexArray = 0;
        emptyArrays();
        CylinderPoint cylinderPoint(startPoint, startColor, normalOne, normalTwo,radiusOne, radiusTwo, CylinderPoint::CylinderPointTypeStart, 0);
        setStartPoint(cylinderPoint);
    };
    void setStartPoint(CylinderPoint &cylinderPoint) {
        emptyArrays();
        addPoint(cylinderPoint);
    };
    void setAngularSampling(const int newSampling){
        angularSampling = newSampling;
    };
    void addPoint(CylinderPoint &cylinderPoint) {
        points.push_back(cylinderPoint);
    };
    void addHalfAtomBond(mmdb::Atom* atom1, FCXXCoord &color1, mmdb::Atom* atom2, FCXXCoord &color2, float cylinderRadius);
    void addHalfAtomBondWithCoords(FCXXCoord &atom1Coord, mmdb::Atom* atom1, FCXXCoord &color1, FCXXCoord &atom2Coord, mmdb::Atom* atom2, FCXXCoord &color2, float cylinderRadius);
    void emptyArrays(){
        delete [] vertexColorNormalArray;
        vertexColorNormalArray = 0;
    };
    virtual ~CylindersPrimitive(){
        emptyArrays();
        //std::cout << "In cylinders destructor" << std::endl;
    };
    virtual void generateArrays();
    
    virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer);
};

#endif
