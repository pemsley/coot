/*
 * MoleculesToTriangles/CXXClasses/FlatFanPrimitive.cpp
 *
 * Copyright 2017 by Martin Noble, University of Oxford
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

#include "FlatFanPrimitive.h"

void FlatFanPrimitive::generateArrays()
{

    _nVertices = 2*(1+atoms.size());
    vertexColorNormalArray = new VertexColorNormal[_nVertices];
    _nTriangles = 2*atoms.size();
    indexArray = new GLIndexType[3*_nTriangles];

    FCXXCoord centre(0.,0.,0.,0.);
    FCXXCoord summedNormals(0.,0.,0.,0.);
    FCXXCoord colorRed(1.,0.,0.,1.);

    size_t iAtom = 0;
    int index = 0;
    for (; iAtom<atoms.size(); iAtom++){
        mmdb::Atom* atom = atoms[iAtom];
        FCXXCoord atomCoord(atom->x, atom->y, atom->z);
        centre += atomCoord;
        size_t lastIAtom = iAtom-1;
        if (iAtom == 0) lastIAtom = atoms.size()-1;
        mmdb::Atom* lastAtom = atoms[lastIAtom];
        FCXXCoord lastAtomCoord(lastAtom->x, lastAtom->y, lastAtom->z);

        size_t nextIAtom = iAtom+1;
        if (nextIAtom == atoms.size()) nextIAtom = 0;
        mmdb::Atom* nextAtom = atoms[nextIAtom];
        FCXXCoord nextAtomCoord(nextAtom->x, nextAtom->y, nextAtom->z);
        
        FCXXCoord vecAB(lastAtomCoord-atomCoord);
        FCXXCoord vecAC(nextAtomCoord-atomCoord);
        FCXXCoord normal(vecAC^vecAB);
        normal.normalise();
        
        summedNormals += normal;
        
        for (int i=0; i<4; i++){
            vertexColorNormalArray[2*iAtom].vertex[i] = (atomCoord+(normal*0.05))[i];
            vertexColorNormalArray[2*iAtom].color[i] = color[i]*255;
            vertexColorNormalArray[2*iAtom].normal[i] = normal[i];
        }
        for (int i=0; i<4; i++){
            vertexColorNormalArray[1+(2*iAtom)].vertex[i] = (atomCoord-(normal*0.05))[i];
            vertexColorNormalArray[1+(2*iAtom)].color[i] = color[i]*255;
            vertexColorNormalArray[1+(2*iAtom)].normal[i] = -normal[i];
        }
        
        indexArray[index++] = 2*iAtom;
        indexArray[index++] = 2*atoms.size();
        indexArray[index++] = 2*nextIAtom;

        indexArray[index++] = 1+(2*nextIAtom);
        indexArray[index++] = 1+(2*atoms.size());
        indexArray[index++] = 1+(2*iAtom);

    }
    centre /= atoms.size();
    summedNormals /= atoms.size();
    summedNormals.normalise();
    
    for (int i=0; i<4; i++){
        vertexColorNormalArray[2*atoms.size()].vertex[i] = (centre+(summedNormals*0.05))[i];
        vertexColorNormalArray[2*atoms.size()].color[i] = color[i]*255.;
        vertexColorNormalArray[2*atoms.size()].normal[i] = summedNormals[i];
    }

    for (int i=0; i<4; i++){
        vertexColorNormalArray[1+(2*atoms.size())].vertex[i] = (centre-(summedNormals*0.05))[i];
        vertexColorNormalArray[1+(2*atoms.size())].color[i] = color[i]*255.;
        vertexColorNormalArray[1+(2*atoms.size())].normal[i] = -summedNormals[i];
    }
};
