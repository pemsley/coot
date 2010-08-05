/*
 *  CXXSurface.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXXSurface_included
#define CXXSurface_included
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#ifndef  __MMDB_Manager__
#include "mmdb_manager.h"
#endif

using namespace std;

class CXXSphereElement;
//#include "CXXSphereElement.h"
//#include "CXXSphereNode.h"
#include "CXXTriangle.h"
#include "CXXCoord.h"

//#include "CXXAlloc.h"

typedef map<string, int> StringIntMap;
class CXXSurface{
private:
    string name;
    StringIntMap vectors;
    StringIntMap scalars;
    StringIntMap pointers;
    vector<CXXTriangle, CXX::CXXAlloc<CXXTriangle> > triangles;
    vector<CXXSurfaceVertex, CXX::CXXAlloc<CXXSurfaceVertex> > vertices;
    PCMMDBManager allAtomsManager;
    int init();
    int nTriangles;
    char fileName[512];
    double getAtomRadius(PCAtom);
    
public:
    CXXSurface ();
    ~CXXSurface ();
    CXXSurface (string path);
    CXXSurface (PCMMDBManager, const std::string selectionString);
    CXXSurface (PCMMDBManager, const std::string selectionString, const std::string contextString);
    CXXSurface (PCMMDBManager, const std::string selectionString, const std::string contextString, const double delta, const double probeRadius, const bool blend_edges);
    CXXSurface (PCMMDBManager, const int);
    CXXSurface (PCMMDBManager, const int, const double, const double, const bool);
    CXXSurface (PCMMDBManager, const int, const int, const double, const double, const bool);
    CXXSurface (int);
    std::string report(); 
    int writeAsGrasp(const std::string &path);
    int readGraspFile(std::string fileName);
    
    const vector<CXXTriangle, CXX::CXXAlloc<CXXTriangle> >&getTriangles() const {
        return triangles;
    };
    const vector<CXXSurfaceVertex, CXX::CXXAlloc<CXXSurfaceVertex> >&getVertices() const {
        return vertices;
    };
	const StringIntMap &vectorNames() const{
		return vectors;
	};
	const StringIntMap &scalarNames() const{
		return scalars;
	};
	const StringIntMap &pointerNames() const{
		return pointers;
	};
    
    const CXXCoord &coordRef(int coordType, int iTriangle, int corner) const;
    const CXXCoord &coordRef(int coordType, int iVertex) const;
    int getCoord(const string &type, const int iTriangle, const int corner, double *buffer);
    int getCoord(const string &type, const int iVertex, double *buffer);
    int getCoord(const int handle, const int iVertex, double *buffer);
    int getPointer(const string &type, int iVertex, void **return_p);
    CMMDBManager *getMMDBManager() const;
    int setCoord (const string &type, int iVertex, const CXXCoord &crd);
    void setScalar (int scalarHandle, int iVertex, double &value);
    void setScalar (const std::string name, int iVertex, double &value);
    
    int addPerVertexVector (const string name, double *vectorBuffer);
    int addPerVertexScalar (const string name, double *scalarBuffer);
    int addPerVertexPointer (const string name, void **pointerBuffer);
    
    int extendWithVectorData(int count, const string name, double *vectorBuffer);
    int updateWithVectorData(int count, const string name, int start, double *data);
    int updateWithPointerData(int count, const string name, int start, void **data);
    int extendTriangles(int *triangleBuffer, int count);
	
    int assignAtom(PCMMDBManager, int);
    int colorByAssignedAtom();
    int colorByColourArray(const std::vector<double*> &colours, CMMDBManager *molHnd, int selHnd);
	
    int numberOfTriangles() const;
    int numberOfVertices() const;
    
    int vertex(int iTriangle, int iCorner) const;

   int calculateFromAtoms(PCMMDBManager, const std::string , const std::string , const double, const double, const bool);
    int calculateFromAtoms(PCMMDBManager, const int, const int, const double, const double, const bool );
    int calculateVDWFromAtoms(PCMMDBManager, const int, const int, const double, const double, const bool );
    int calculateAccessibleFromAtoms(PCMMDBManager, const int, const int, const double, const double, const bool );
    int calculateFromAtoms(PCMMDBManager, const std::string , const double, const double, const bool);
    int calculateFromAtoms(PCMMDBManager, const int, const double, const double, const bool);
 	
    int upLoadSphere(CXXSphereElement &theSphere, double probeRadius, const int sense);
    int selectionStringToSelHnd(PCMMDBManager, const std::string selectionString);
    int getIntegerUDDataOfAtom(PCAtom theAtom, int handle);
    int operator == (const CXXSurface &comparator) const{
        return (this == &comparator);
    } 
    int getVectorHandle(const string name);
    int getReadVectorHandle(const string name);
    int getScalarHandle(const string name);
    int getReadScalarHandle(const string name);
    int getPointerHandle(const string name);
	
    int getScalar(int handle, int iVertex, double &result);
    int getVector(int handle, int iVertex, double *result);
    
    int addTriangle(const CXXTriangle &aTriangle);
    
    int assignUnitedAtomRadius();
    
    void appendSurface(const CXXSurface &otherSurface);
    void compress(double tolerance);
};

#endif




