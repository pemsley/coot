/*
 * MoleculesToTriangles/CXXClasses/FlatFanPrimitive-gemmi.cc
 *
 * gemmi-native twin of FlatFanPrimitive.cpp - identical geometry, fed from
 * vertex positions instead of mmdb::Atom coordinates.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "FlatFanPrimitive-gemmi.hh"

void coot::m2t::FlatFanPrimitive::generateArrays()
{
   _nVertices = 2 * (1 + positions.size());
   vertexColorNormalArray = new VertexColorNormal[_nVertices];
   _nTriangles = 2 * positions.size();
   indexArray = new GLIndexType[3 * _nTriangles];

   FCXXCoord centre(0., 0., 0., 0.);
   FCXXCoord summedNormals(0., 0., 0., 0.);

   size_t iAtom = 0;
   int index = 0;
   for (; iAtom < positions.size(); iAtom++) {
      FCXXCoord atomCoord = positions[iAtom];
      centre += atomCoord;
      size_t lastIAtom = iAtom - 1;
      if (iAtom == 0) lastIAtom = positions.size() - 1;
      FCXXCoord lastAtomCoord = positions[lastIAtom];

      size_t nextIAtom = iAtom + 1;
      if (nextIAtom == positions.size()) nextIAtom = 0;
      FCXXCoord nextAtomCoord = positions[nextIAtom];

      FCXXCoord vecAB(lastAtomCoord - atomCoord);
      FCXXCoord vecAC(nextAtomCoord - atomCoord);
      FCXXCoord normal(vecAC ^ vecAB);
      normal.normalise();

      summedNormals += normal;

      for (int i = 0; i < 4; i++) {
         vertexColorNormalArray[2*iAtom].vertex[i] = (atomCoord + (normal*0.05))[i];
         vertexColorNormalArray[2*iAtom].color[i] = color[i] * 255;
         vertexColorNormalArray[2*iAtom].normal[i] = normal[i];
      }
      for (int i = 0; i < 4; i++) {
         vertexColorNormalArray[1+(2*iAtom)].vertex[i] = (atomCoord - (normal*0.05))[i];
         vertexColorNormalArray[1+(2*iAtom)].color[i] = color[i] * 255;
         vertexColorNormalArray[1+(2*iAtom)].normal[i] = -normal[i];
      }

      indexArray[index++] = 2*iAtom;
      indexArray[index++] = 2*positions.size();
      indexArray[index++] = 2*nextIAtom;

      indexArray[index++] = 1+(2*nextIAtom);
      indexArray[index++] = 1+(2*positions.size());
      indexArray[index++] = 1+(2*iAtom);
   }
   centre /= positions.size();
   summedNormals /= positions.size();
   summedNormals.normalise();

   for (int i = 0; i < 4; i++) {
      vertexColorNormalArray[2*positions.size()].vertex[i] = (centre + (summedNormals*0.05))[i];
      vertexColorNormalArray[2*positions.size()].color[i] = color[i] * 255.;
      vertexColorNormalArray[2*positions.size()].normal[i] = summedNormals[i];
   }
   for (int i = 0; i < 4; i++) {
      vertexColorNormalArray[1+(2*positions.size())].vertex[i] = (centre - (summedNormals*0.05))[i];
      vertexColorNormalArray[1+(2*positions.size())].color[i] = color[i] * 255.;
      vertexColorNormalArray[1+(2*positions.size())].normal[i] = -summedNormals[i];
   }
}
