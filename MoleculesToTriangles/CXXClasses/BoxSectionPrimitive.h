//
//  BoxSectionPrimitive.h
//  AesopCD_macos
//
//  Created by Martin Noble on 11/08/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
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
