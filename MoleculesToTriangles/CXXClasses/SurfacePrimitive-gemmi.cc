/*
 * MoleculesToTriangles/CXXClasses/SurfacePrimitive-gemmi.cc
 *
 * gemmi-native twin of SurfacePrimitive.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "SurfacePrimitive-gemmi.hh"
#include "MoleculesToTriangles/CXXSurface/CXXSurface-gemmi.hh"
#include "MoleculesToTriangles/CXXSurface/CXXSurfaceVertex.h"
#include <map>
#include <cmath>

using namespace coot::m2t;

SurfacePrimitive::SurfacePrimitive() : VertexColorNormalPrimitive() {
   drawModeGL = DrawAsTriangles;
   enableColorGL = true;
   cxxSurfaceMaker = 0;
   primitiveType = DisplayPrimitive::PrimitiveType::SurfacePrimitive;
}

void SurfacePrimitive::generateArrays()
{
   CXXSurfaceMaker &containerSurface(*cxxSurfaceMaker);
   double coords[4];

   _nVertices = 0;
   _nTriangles = 0;
   for (std::vector<CXXSurface>::iterator surfIter = containerSurface.getChildSurfaces().begin();
        surfIter != containerSurface.getChildSurfaces().end(); ++surfIter) {
      _nVertices += surfIter->numberOfVertices();
      _nTriangles += surfIter->numberOfTriangles();
   }

   vertexColorNormalArray = new VertexColorNormal[_nVertices];

   size_t iGlobalVert = 0;
   for (std::vector<CXXSurface>::iterator surfIter = containerSurface.getChildSurfaces().begin();
        surfIter != containerSurface.getChildSurfaces().end(); ++surfIter) {
      CXXSurface &childSurface = *surfIter;
      for (size_t iLocalVert = 0; iLocalVert < childSurface.numberOfVertices(); iLocalVert++) {
         VertexColorNormal &vcn(vertexColorNormalArray[iGlobalVert]);
         if (!childSurface.getCoord("vertices", iLocalVert, coords)) {
            for (size_t k = 0; k < 3; k++) vcn.vertex[k] = coords[k];
         }
         vcn.vertex[3] = 1.;
         if (!childSurface.getCoord("normals", iLocalVert, coords)) {
            for (size_t k = 0; k < 3; k++) vcn.normal[k] = coords[k];
         }
         vcn.normal[3] = 1.;
         if (!childSurface.getCoord("colour", iLocalVert, coords)) {
            for (int k = 0; k < 4; k++) {
               float floatComponentValue = coords[k] * 255.;
               int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255. ? 255 : floatComponentValue));
               vcn.color[k] = uintComponentValue;
            }
         }
         else for (int k = 0; k < 4; k++) vcn.color[k] = 0.5;
         iGlobalVert++;
      }
   }
   indexArray = new GLIndexType[3 * _nTriangles];
   size_t iOffset = 0;
   int idx = 0;
   for (std::vector<CXXSurface>::iterator surfIter = containerSurface.getChildSurfaces().begin();
        surfIter != containerSurface.getChildSurfaces().end(); ++surfIter) {
      CXXSurface &childSurface = *surfIter;
      for (std::size_t i = 0; i < childSurface.numberOfTriangles(); i++) {
         for (int j = 0; j < 3; j++) {
            indexArray[idx++] = GLIndexType(childSurface.vertex(i, j) + iOffset);
         }
      }
      iOffset += childSurface.numberOfVertices();
   }

   delete cxxSurfaceMaker;
   cxxSurfaceMaker = 0;
}

SurfacePrimitive::SurfacePrimitive(gemmi::Structure *structure, gemmi::Model &model,
                                   const std::set<const gemmi::Atom*> &selSet,
                                   const std::set<const gemmi::Atom*> &contextSet,
                                   std::shared_ptr<ColorScheme> _colorScheme, enum SurfaceType type,
                                   float probeRadius, float radiusMultiplier) {
   primitiveType = DisplayPrimitive::PrimitiveType::SurfacePrimitive;
   cxxSurfaceMaker = 0;
   vertexColorNormalArray = 0;
   indexArray = 0;
   atomArray = 0;
   colorScheme = _colorScheme;
   cxxSurfaceMaker = new CXXSurfaceMaker();

   double delta = 30. * (M_PI / 180.);

   // atom -> CRA map for per-vertex colour
   std::map<const gemmi::Atom *, gemmi::CRA> craMap;
   for (gemmi::Chain &chain : model.chains)
      for (gemmi::Residue &res : chain.residues)
         for (gemmi::Atom &atom : res.atoms)
            craMap[&atom] = gemmi::CRA{&chain, &res, &atom};

   try {
      if (type == AccessibleSurface)
         cxxSurfaceMaker->calculateAccessibleFromAtoms(structure, model, selSet, contextSet, delta, probeRadius, radiusMultiplier);
      else if (type == VdWSurface)
         cxxSurfaceMaker->calculateVDWFromAtoms(structure, model, selSet, contextSet, delta, probeRadius, radiusMultiplier);
      else if (type == MolecularSurface)
         cxxSurfaceMaker->calculateFromAtoms(structure, model, selSet, contextSet, delta, probeRadius, false);
   }
   catch (std::exception &e) {
      std::cout << "exception caught in SurfacePrimitive: " << e.what() << std::endl;
      delete cxxSurfaceMaker;
      cxxSurfaceMaker = 0;
   }

   if (cxxSurfaceMaker) {
      CXXSurfaceMaker &mySurfaceMaker = *cxxSurfaceMaker;
      for (std::vector<CXXSurface>::iterator surfIter = mySurfaceMaker.getChildSurfaces().begin();
           surfIter != mySurfaceMaker.getChildSurfaces().end(); ++surfIter) {
         CXXSurface &mySurface = *surfIter;
         for (std::size_t iVertex = 0; iVertex < mySurface.numberOfVertices(); iVertex++) {
            void *p = 0;
            int result = mySurface.getPointer("atom", iVertex, &p);
            const gemmi::Atom *theAtom = static_cast<const gemmi::Atom *>(p);
            if (result == 0 && theAtom) {
               std::map<const gemmi::Atom *, gemmi::CRA>::iterator craIter = craMap.find(theAtom);
               if (craIter != craMap.end()) {
                  FCXXCoord color = colorScheme->colorForAtom(craIter->second);
                  mySurface.setCoord("colour", iVertex, CXXCoord<CXXCoord_ftype>(color[0], color[1], color[2], color[3]));
               }
            }
         }
      }
   }
}
