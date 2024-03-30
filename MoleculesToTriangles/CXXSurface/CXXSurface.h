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
#include "mmdb2/mmdb_manager.h"
#endif

using namespace std;

class CXXSphereElement;
class CXXTorusElement;
#include "CXXTriangle.h"
#include "CXXCoord.h"

class SurfaceParameters{
public:
    size_t nVertices;
    size_t nTriangles;
    double MSA;
    double ASA;
    std::map<std::string, double> pMins;
    std::map<std::string, double> pMaxes;
    std::map<std::string, double> pMeans;
    
    SurfaceParameters() : nVertices(0), nTriangles(0), MSA(0.),ASA(0.){};
    SurfaceParameters& operator += (const SurfaceParameters &otherParameters){
        MSA += otherParameters.MSA;
        ASA += otherParameters.ASA;
        for (std::map<std::string,double>::const_iterator prop=otherParameters.pMins.begin();
             prop!=otherParameters.pMins.end(); prop++){
            if (pMins.find(prop->first)!=pMins.end()) pMins[prop->first] = min(pMins[prop->first], prop->second);
            else pMins[prop->first] = prop->second;
        }
        for (std::map<std::string,double>::const_iterator prop=otherParameters.pMaxes.begin();
             prop!=otherParameters.pMaxes.end(); prop++){
            if (pMaxes.find(prop->first)!=pMaxes.end()) pMaxes[prop->first] = max(pMaxes[prop->first], prop->second);
            else pMaxes[prop->first] = prop->second;
        }
        for (std::map<std::string,double>::const_iterator prop=otherParameters.pMeans.begin();
             prop!=otherParameters.pMeans.end(); prop++){
            if (pMeans.find(prop->first)!=pMeans.end()) {
                pMeans[prop->first] = (pMeans[prop->first] + prop->second) / (nVertices + otherParameters.nVertices);
            }
            else pMeans[prop->first] = prop->second;
        }
        nTriangles += otherParameters.nTriangles;
        nVertices += otherParameters.nVertices;
        return *this;
    };
};


typedef map<string, size_t> StringIntMap;
class CXXSurface{
private:
    string name;
    StringIntMap vectors;
    StringIntMap scalars;
    StringIntMap pointers;
    vector<CXXTriangle  > triangles;
    vector<CXXSurfaceVertex> vertices;
    mmdb::Manager* allAtomsManager;
    int init();
    size_t nTriangles;
    char fileName[512];
public:
    std::string report();
    int writeAsGrasp(string path);
    int readGraspFile(std::string fileName);
    
    const vector<CXXTriangle  >&getTriangles() const {
        return triangles;
    };
    const vector<CXXSurfaceVertex>&getVertices() const {
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
    
    const CXXCoord<CXXCoord_ftype>&coordRef(size_t coordType, size_t iTriangle, size_t corner) const;
    const CXXCoord<CXXCoord_ftype>&coordRef(size_t coordType, size_t iVertex) const;
    int getCoord(const string &type, const size_t iTriangle, const size_t corner, double *buffer);
    int getCoord(const string &type, const size_t iVertex, double *buffer);
    int getCoord(const size_t handle, const size_t iVertex, double *buffer);
    int getPointer(const string &type, size_t iVertex, void **return_p);
    mmdb::Manager* getMMDBManager() const;
    size_t setCoord (const string &type, size_t iVertex, const CXXCoord<CXXCoord_ftype>&crd);
    void setScalar (size_t scalarHandle, size_t iVertex, double &value);
    void setScalar (const std::string name, size_t iVertex, double &value);
    
    int addPerVertexVector (const string name, double *vectorBuffer);
    int addPerVertexScalar (const string name, double *scalarBuffer);
    int addPerVertexPointer (const string name, void **pointerBuffer);
    
    size_t extendWithVectorData(size_t count, const string name, double *vectorBuffer);
    size_t updateWithVectorData(size_t count, const string name, size_t start, double *data);
    size_t updateWithPointerData(size_t count, const string name, size_t start, void **data);
    size_t extendTriangles(int *triangleBuffer, int count);
	
	
    size_t numberOfTriangles() const;
    size_t numberOfVertices() const;
    
    size_t vertex(size_t iTriangle, size_t iCorner) const;
 	
    int upLoadSphere(CXXSphereElement &theSphere, double probeRadius, const int sense);
    int uploadTorus(CXXTorusElement &theTorus);
    int selectionStringToSelHnd(mmdb::Manager*, const std::string selectionString);
    int getIntegerUDDataOfAtom(mmdb::Atom* theAtom, int handle);
    int operator == (const CXXSurface &comparator) const{
        return (this == &comparator);
    } 
    size_t getVectorHandle(const string name);
    size_t getReadVectorHandle(const string name);
    size_t getScalarHandle(const string name);
    size_t getReadScalarHandle(const string name);
    size_t getPointerHandle(const string name);
	
    int getScalar(int handle, int iVertex, double &result);
    int getVector(int handle, int iVertex, double *result);
    
    int addTriangle(const CXXTriangle &aTriangle);
    
    int assignUnitedAtomRadius();
    
    void compress(double tolerance);
    SurfaceParameters measuredProperties();

    int assignAtom(mmdb::Manager*, int);
    int colorByAssignedAtom();
    int colorByColourArray(const std::vector<double*> &colours, mmdb::Manager* molHnd, int selHnd);
    
    void appendSurface(const CXXSurface &otherSurface);
    

};

#endif




