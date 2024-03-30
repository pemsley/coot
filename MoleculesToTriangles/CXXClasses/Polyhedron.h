/*
 * MoleculesToTriangles/CXXClasses/Polyhedron.h
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

#ifndef Polyhedron_H
#define Polyhedron_H

#include <string>
#include <vector>
#include <map>

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

class PolyhedronFace {
private:
	std::vector<unsigned int> indices;
protected:
	std::vector<unsigned int> &getIndices(){
		return indices;
	}
public:
	unsigned int &operator [] (const unsigned int iIndex){
		return indices[iIndex];
	}
	unsigned int operator [] (const unsigned int iIndex) const{
		return indices[iIndex];
	}
	unsigned long nIndices() const{
		return indices.size();
	}
	PolyhedronFace(){
		indices.resize(3);
	}
	PolyhedronFace(const PolyhedronFace &oldOne){
		*this = oldOne;
	}
	PolyhedronFace(int indexCount, int *newIndices){
		if (indices.size() != indexCount) indices.resize(indexCount);
		for (int i=0; i<indexCount; i++){
			indices[i] = newIndices[i];
		}
	}
};

class Polyhedron  {	
private:
protected:
	std::vector<FCXXCoord  > vertices;
	std::vector<PolyhedronFace> faces;

	std::vector<FCXXCoord  > &getVertices(){
		return vertices;
	}
	std::vector<PolyhedronFace> &getFaces(){
		return faces;
	}
public:
	enum PolyhedronType {
		Octahedron, Dodecahedron, Icosahedron, Dodecahedron4, Dodecahedron16
	};
	unsigned long nVertices() const {
		return vertices.size();
	}
	FCXXCoord vertex(int iVertex) const {
		return vertices[iVertex];
	}
	unsigned long nFaces() const {
		return faces.size();
	}
	PolyhedronFace face(int iFace) const{
		return faces[iFace];
	}
	PolyhedronFace &face(int iFace) {
		return faces[iFace];
	}
	Polyhedron(){
	}
	Polyhedron (const Polyhedron &oldOne){
		*this = oldOne;
	}
    ~Polyhedron(){
        
    };
};

#endif
