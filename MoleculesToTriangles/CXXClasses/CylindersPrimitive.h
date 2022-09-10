/*
 *  CylindersPrimitive.h
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 19/07/2008.
 *  Copyright 2008 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
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
