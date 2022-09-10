//
//  FlatFanPrimitive.cpp
//  MoleculesToTriangles
//
//  Created by Martin Noble on 16/03/2017.
//  Copyright © 2017 MartinNobleSoftware. All rights reserved.
//

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
