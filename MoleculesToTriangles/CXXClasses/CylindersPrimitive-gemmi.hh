/*
 * MoleculesToTriangles/CXXClasses/CylindersPrimitive-gemmi.hh
 *
 * gemmi-native twin of CylindersPrimitive.h. Uses an atom-free CylinderPoint
 * (the original stored a const mmdb::Atom* for picking only) and a coord-based
 * bond-adding method. Colour is supplied per end (evaluated eagerly by the
 * caller). Type is coot::m2t::CylindersPrimitive.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CylindersPrimitive_gemmi_hh
#define CylindersPrimitive_gemmi_hh

#include <memory>
#include <vector>
#include "VertexColorNormalPrimitive.h"
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

class Renderer;

namespace coot {
   namespace m2t {

      // atom-free twin of CylinderPoint (picking pointer dropped)
      class CylinderPoint {
      public:
         enum CylinderPointType { CylinderPointTypeStart, CylinderPointTypeContinue };
         FCXXCoord vertex;
         FCXXCoord color;
         FCXXCoord normalOne;
         FCXXCoord normalTwo;
         float radiusOne;
         float radiusTwo;
         int type;
         CylinderPoint(FCXXCoord _vertex, FCXXCoord _color, FCXXCoord _normalOne, FCXXCoord _normalTwo,
                       float _radiusOne, float _radiusTwo)
            : vertex(_vertex), color(_color), normalOne(_normalOne), normalTwo(_normalTwo),
              radiusOne(_radiusOne), radiusTwo(_radiusTwo), type(CylinderPointTypeContinue) {}
         CylinderPoint(FCXXCoord _vertex, FCXXCoord _color, FCXXCoord _normalOne, FCXXCoord _normalTwo,
                       float _radiusOne, float _radiusTwo, int _type)
            : vertex(_vertex), color(_color), normalOne(_normalOne), normalTwo(_normalTwo),
              radiusOne(_radiusOne), radiusTwo(_radiusTwo), type(_type) {}
      };

      class CylindersPrimitive : public VertexColorNormalPrimitive {
      private:
         friend class BoxSectionPrimitive;     // gemmi BoxSection twin (ribbon)
         std::vector<CylinderPoint> points;
         int angularSampling;
      public:
         CylindersPrimitive();
         void setStartPoint(CylinderPoint &cylinderPoint) { emptyArrays(); addPoint(cylinderPoint); }
         void setAngularSampling(const int newSampling) { angularSampling = newSampling; }
         void addPoint(CylinderPoint &cylinderPoint) { points.push_back(cylinderPoint); }
         // coord-based half-bond (colour supplied per atom, evaluated by the caller)
         void addHalfAtomBondWithCoords(FCXXCoord &atom1Coord, FCXXCoord &color1,
                                        FCXXCoord &atom2Coord, FCXXCoord &color2, float cylinderRadius);
         void emptyArrays() { delete [] vertexColorNormalArray; vertexColorNormalArray = 0; }
         virtual ~CylindersPrimitive() { emptyArrays(); }
         virtual void generateArrays();
         virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer);
      };
   }
}

#endif // CylindersPrimitive_gemmi_hh
