/*
 * MoleculesToTriangles/CXXClasses/BoxSectionPrimitive.h
 *
 * Copyright 2011 by Martin Noble, University of Oxford
 * Author: Martin Noble
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

//
//  BoxSectionPrimitive.h
//  AesopCD_macos
//

#ifndef BoxSectionPrimitive_h
#define BoxSectionPrimitive_h

#include "CylindersPrimitive.h"
#include "CylinderPoint.h"
#include <memory>

class BoxSectionPrimitive : public CylindersPrimitive {
public:
    BoxSectionPrimitive();
    BoxSectionPrimitive(FCXXCoord &startPoint, FCXXCoord &startColor, FCXXCoord &normalOne, FCXXCoord &normalTwo,
                       float radiusOne, float radiusTwo){
        vertexColorNormalArray = 0;
        indexArray = 0;
        emptyArrays();
        CylinderPoint cylinderPoint(startPoint, startColor, normalOne, normalTwo,radiusOne, radiusTwo, CylinderPoint::CylinderPointTypeStart, 0);
        setStartPoint(cylinderPoint);
        primitiveType = DisplayPrimitive::PrimitiveType::BoxSectionPrimitive;
        boxEdgeColor.first = false; // the default constructor may do this anyway (not sure).
    }
    virtual ~BoxSectionPrimitive(){
        emptyArrays();
        //std::cout << "In box section destructor" << std::endl;
    };
    //Override generate and render routines
    void generateArrays();
    void renderWithRenderer(std::shared_ptr<Renderer> renderer);
    std::pair<bool, FCXXCoord> boxEdgeColor;
    void setBoxEdgeColor(const FCXXCoord &col_in) { boxEdgeColor.first = true; boxEdgeColor.second = col_in; }
    std::pair<bool, FCXXCoord> getBoxEdgeColor() { return boxEdgeColor; }

};



#endif
