//
//  SticksPrimitive.cpp
//  AesopCD_ios
//
//  Created by Martin Noble on 14/02/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
//

#include "SticksPrimitive.h"
#include "mmdb2/mmdb_manager.h"
#include "ColorScheme.h"
#include "Renderer.h"

void SticksPrimitive::generateArrays()
{
    int nBonds = 0;
    int nUniqueAtoms = 0;
    std::map<mmdb::Atom *, VertexColor>vcLookup;
    
    if (sticks.begin() != sticks.end()){
        //mmdb::Atom *atom1 = sticks.begin()->first;
        std::map<std::shared_ptr<ColorRule>,int> handles = colorScheme->prepareForMMDB(mmdb);
        
        {
            //First pass :Count the total bonds and make unique atom vertex list
            std::map<mmdb::Atom *, std::vector<mmdb::Atom *> >::iterator mapIter = sticks.begin();
            std::map<mmdb::Atom *, std::vector<mmdb::Atom *> >::iterator mapEnd = sticks.end();
            for (; mapIter != mapEnd; ++mapIter){
                mmdb::Atom *atom1 = mapIter->first;
                std::map<mmdb::Atom *, VertexColor>::iterator atom1Iter = vcLookup.find(atom1);
                if (atom1Iter == vcLookup.end()){
                    FCXXCoord color1 = colorScheme->colorForAtom(atom1, handles);
                    VertexColor vc1 = {
                        {static_cast<float>(atom1->x), static_cast<float>(atom1->y), static_cast<float>(atom1->z), 0.},
                        {static_cast<float>(color1[0]), static_cast<float>(color1[1]), static_cast<float>(color1[2]), static_cast<float>(color1[3])}
                    };
                    vcLookup[atom1] = vc1;
                    nUniqueAtoms++;
                }
                
                std::vector<mmdb::Atom *> &bondsOfAtom(mapIter->second);
                std::vector<mmdb::Atom *>::iterator vecIter = bondsOfAtom.begin();
                std::vector<mmdb::Atom *>::iterator vecEnd = bondsOfAtom.end();
                
                for (; vecIter != vecEnd; ++vecIter){
                    mmdb::Atom *atom2 = *vecIter;
                    std::map<mmdb::Atom *, VertexColor>::iterator atom2Iter = vcLookup.find(atom2);
                    if (atom2Iter == vcLookup.end()){
                        FCXXCoord color2 = colorScheme->colorForAtom(atom2, handles);
                        VertexColor vc2 = {
                            {static_cast<float>(atom2->x), static_cast<float>(atom2->y), static_cast<float>(atom2->z), 0.},
                            {static_cast<float>(color2[0]), static_cast<float>(color2[1]), static_cast<float>(color2[2]), static_cast<float>(color2[3])}
                        };
                        vcLookup[atom2] = vc2;
                        nUniqueAtoms++;
                    }
                    nBonds++;
                }    
            }
        }
        
        int requestedVertices = nUniqueAtoms + 2*nBonds;
        int requestedIndices = 4 * nBonds;
        vertexColorArray = new VertexColor[requestedVertices];
        indexArray = new GLIndexType[requestedIndices];
        
        std::map<mmdb::Atom *, VertexColor>::iterator mapIter = vcLookup.begin();
        std::map<mmdb::Atom *, VertexColor>::iterator mapEnd = vcLookup.end();
        _nVertices = 0;
        std::map<mmdb::Atom *, unsigned long>indexLookup;
        for (; mapIter != mapEnd; ++mapIter){
            indexLookup[mapIter->first] = _nVertices;
            vertexColorArray[_nVertices] = mapIter->second;
            _nVertices++;
        }
        
        _nLines = 0;
        int iIndex = 0;
        {
            //Second pass : Fill the arrays
            std::map<mmdb::Atom *, std::vector<mmdb::Atom *> >::iterator mapIter = sticks.begin();
            std::map<mmdb::Atom *, std::vector<mmdb::Atom *> >::iterator mapEnd = sticks.end();
            for (; mapIter != mapEnd; ++mapIter){
                mmdb::Atom *atom1 = mapIter->first;
                unsigned long atom1Index = indexLookup[atom1];
                FCXXCoord atom1Coord(atom1->x, atom1->y, atom1->z);
                FCXXCoord color1 = colorScheme->colorForAtom(atom1, handles);
                
                std::vector<mmdb::Atom *> &bondsOfAtom(mapIter->second);
                std::vector<mmdb::Atom *>::iterator vecIter = bondsOfAtom.begin();
                std::vector<mmdb::Atom *>::iterator vecEnd = bondsOfAtom.end();
                for (; vecIter != vecEnd; ++vecIter){
                    mmdb::Atom *atom2 = *vecIter;
                    unsigned long atom2Index = indexLookup[atom2];
                    FCXXCoord atom2Coord(atom2->x, atom2->y, atom2->z);
                    FCXXCoord color2 = colorScheme->colorForAtom(atom2, handles);

                    FCXXCoord midCoord((atom1Coord + atom2Coord) * 0.5);
                    {
                        VertexColor vcMid = {
                            {static_cast<float>(midCoord[0]), static_cast<float>(midCoord[1])
                                
                                
                                , static_cast<float>(midCoord[2]), 0.},
                            {static_cast<float>(color1[0]), static_cast<float>(color1[1]), static_cast<float>(color1[2]), static_cast<float>(color1[3])}
                        };
                        vertexColorArray[_nVertices] = vcMid;
                        indexArray[iIndex++] = (GLIndexType)atom1Index;
                        indexArray[iIndex++] = (GLIndexType)_nVertices;
                        _nLines++;
                        _nVertices++;
                    }    
                    {
                        VertexColor vcMid = { 
                            {static_cast<float>(midCoord[0]), static_cast<float>(midCoord[1]), static_cast<float>(midCoord[2]), 0.}, 
                            {static_cast<float>(color2[0]), static_cast<float>(color2[1]), static_cast<float>(color2[2]), static_cast<float>(color2[3])}
                        };
                        vertexColorArray[_nVertices] = vcMid;
                        indexArray[iIndex++] = (GLIndexType)atom2Index;
                        indexArray[iIndex++] = (GLIndexType)_nVertices;
                        _nLines++;
                        _nVertices++;
                    }    
                }    
            }
        }
        colorScheme->freeSelectionHandles(mmdb, handles);
    }
    std::cout << "nLines " <<  _nLines << " up to " << _nVertices << " vertices" << std::endl;
}

void SticksPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
    if (vertexColorArray == 0){
        generateArrays();
    }
    renderer->renderVertexColorPrimitive(this);
}

