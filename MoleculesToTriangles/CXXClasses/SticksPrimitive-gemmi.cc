/*
 * MoleculesToTriangles/CXXClasses/SticksPrimitive-gemmi.cc
 *
 * gemmi-native twin of SticksPrimitive.cpp. Same half-bond line geometry (one
 * shared vertex per atom, two coloured midpoint vertices per bond), built from
 * caller-supplied coords+colours instead of mmdb atoms + deferred colour.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "SticksPrimitive-gemmi.hh"
#include "Renderer.h"

void coot::m2t::SticksPrimitive::generateArrays()
{
   int nBonds = (int)bondPairs.size();
   int nUniqueAtoms = (int)atomNodes.size();

   int requestedVertices = nUniqueAtoms + 2*nBonds;
   int requestedIndices = 4*nBonds;
   vertexColorArray = new VertexColor[requestedVertices];
   indexArray = new GLIndexType[requestedIndices];

   // one shared vertex per atom node
   std::map<long, unsigned long> indexLookup;
   _nVertices = 0;
   for (std::map<long, Node>::iterator it = atomNodes.begin(); it != atomNodes.end(); ++it) {
      indexLookup[it->first] = _nVertices;
      const Node &n = it->second;
      VertexColor vc = {
         { (float)n.coord[0], (float)n.coord[1], (float)n.coord[2], 0.f },
         { (float)n.color[0], (float)n.color[1], (float)n.color[2], (float)n.color[3] }
      };
      vertexColorArray[_nVertices] = vc;
      _nVertices++;
   }

   _nLines = 0;
   int iIndex = 0;
   for (std::vector<std::pair<long, long> >::iterator bp = bondPairs.begin(); bp != bondPairs.end(); ++bp) {
      std::map<long, unsigned long>::iterator i1 = indexLookup.find(bp->first);
      std::map<long, unsigned long>::iterator i2 = indexLookup.find(bp->second);
      if (i1 == indexLookup.end() || i2 == indexLookup.end()) continue;
      const Node &n1 = atomNodes[bp->first];
      const Node &n2 = atomNodes[bp->second];
      FCXXCoord c1 = n1.coord;
      FCXXCoord c2 = n2.coord;
      FCXXCoord mid = (c1 + c2) * 0.5;

      // half-bond atom1 -> midpoint, coloured by atom1
      VertexColor vcMid1 = {
         { (float)mid[0], (float)mid[1], (float)mid[2], 0.f },
         { (float)n1.color[0], (float)n1.color[1], (float)n1.color[2], (float)n1.color[3] }
      };
      vertexColorArray[_nVertices] = vcMid1;
      indexArray[iIndex++] = (GLIndexType)i1->second;
      indexArray[iIndex++] = (GLIndexType)_nVertices;
      _nLines++; _nVertices++;

      // half-bond atom2 -> midpoint, coloured by atom2
      VertexColor vcMid2 = {
         { (float)mid[0], (float)mid[1], (float)mid[2], 0.f },
         { (float)n2.color[0], (float)n2.color[1], (float)n2.color[2], (float)n2.color[3] }
      };
      vertexColorArray[_nVertices] = vcMid2;
      indexArray[iIndex++] = (GLIndexType)i2->second;
      indexArray[iIndex++] = (GLIndexType)_nVertices;
      _nLines++; _nVertices++;
   }
}

void coot::m2t::SticksPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
   if (vertexColorArray == 0)
      generateArrays();
   renderer->renderVertexColorPrimitive(this);
}
