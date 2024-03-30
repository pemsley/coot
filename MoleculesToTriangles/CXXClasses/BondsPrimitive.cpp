/*
 * MoleculesToTriangles/CXXClasses/BondsPrimitive.cpp
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

#include "BondsPrimitive.h"

#include "mmdb2/mmdb_manager.h"
#include "Renderer.h"
#include "ColorScheme.h"
#include <map>

void BondsPrimitive::evaluateGLPrimitives(std::map<std::shared_ptr<ColorRule>, int> &handles)
{
    invalidateGLPrimitives();
    //Size of the vertex Color array : one node per atom, and two nodes per
    //bond midpoint (to allow for a different color for each halfbond)
    std::cout << "Making space for " << bonds.size() + 2*nBonds << " bond nodes \n";
    vertexColorArray = new VertexColor[bonds.size() + 2*nBonds];
    std::cout << "Making space for " << 2*nBonds << " bond indices \n";
    indexArray = new GLIndexType[2*nBonds];
    std::map<mmdb::Atom *, std::vector<mmdb::Atom *> >::iterator centralAtomPntr = bonds.begin();
    int iAtom = 0;
    int iMidpoint = 0;
    int iIndex = 0;
    unsigned long midpointIndex = 0;
    int iBond = 0;
    for (; centralAtomPntr != bonds.end(); ++centralAtomPntr, iAtom++){
        FCXXCoord atom1Coord(centralAtomPntr->first->x, centralAtomPntr->first->y, centralAtomPntr->first->z, 0.);
        FCXXCoord atom1Color = colorScheme->colorForAtom(centralAtomPntr->first, handles);
        for (int i=0; i<4; i++) {
            vertexColorArray[iAtom].vertex[i] = atom1Coord[i];
            vertexColorArray[iAtom].color[i] = atom1Color[i];
        }
        std::vector<mmdb::Atom *>::iterator bondedAtomPntr = centralAtomPntr->second.begin();
        for (; bondedAtomPntr != centralAtomPntr->second.end(); ++bondedAtomPntr, iMidpoint++){
            midpointIndex = bonds.size() + iMidpoint;
            FCXXCoord atom2Coord((*bondedAtomPntr)->x, (*bondedAtomPntr)->y, (*bondedAtomPntr)->z, 0.);
            FCXXCoord midpoint = (atom1Coord + atom2Coord) / 2.;
            for (int i=0; i<4; i++) {
                vertexColorArray[midpointIndex].vertex[i] = midpoint[i];
                vertexColorArray[midpointIndex].color[i] = atom1Color[i];
            }
            iBond++;
            indexArray[iIndex++] = iAtom;
            indexArray[iIndex++] = GLIndexType(midpointIndex);
        }
    }
    std::cout << "Bond object contains " << nBonds << " bonds around " << bonds.size() << "atoms\n";
    std::cout << "midpointIndex got up to" << midpointIndex << " iBond to " << iBond << " and iIndex to "<< iIndex << std::endl;
}

void BondsPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
    renderer->renderBondsPrimitive(this);
}

void BondsPrimitive::invalidateGLPrimitives()
{
    delete [] vertexColorArray;
    delete [] indexArray;
    vertexColorArray = 0;
    indexArray = 0;
}
