/*
 *  Polyhedron.h
 *  iPhoneRibbons
 *
 *  Created by Martin Noble on 15/09/2008.
 *  Copyright 2008 LMB, Oxford University. All rights reserved.
 *
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
