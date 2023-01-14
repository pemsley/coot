/*
 *  CXXSurface.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <sstream>
#include <math.h>
#include <set>
#include "string.h"
#include "CXXSurface.h"
#include "CXXFortranFile.h"
#include "TokenIterator.h"
#include "CXXNewHood.h"
#include "CXXTorusElement.h"
#include "CXXTorusNode.h"
#include "CXXSphereElement.h"
#include "mmdb2/mmdb_tables.h"
#include "mmdb2/mmdb_uddata.h"
#include "CXXCircle.h"
#include "CXXTriangle.h"
#include "CXXSphereFlatTriangle.h"
#include "CXXSurfaceVertex.h"
#include <sys/time.h>
#include <limits>
#include <stdint.h>
//#include <forward_list>

int CXXSurface::readGraspFile(string path)
{
    char buffer[1024];
    string text;
    float *coordBuffer;
    int *triangleBuffer;
    
    Delimiters delimiters(" \0\t\n~;()\"<>:{}[]+-=&*#.,/\\");
    int hasVertices, hasAccessibles, hasNormals, hasTriangles, hasPotential,
    hasProperty1, hasProperty2;
    
    CXXFortranFile graspFile(path,"r");
    if (!graspFile.bad()) std::cout << "Opened grasp file [" << path << "]\n";
    else {
        std::cout << "Failed to open grasp file [" << path << "]\n";
        return 1;
    }
    //Grasp Surface file format
    graspFile.getFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    text = buffer;
    cout << "[" << text << "]\n";
    if (text.find("format=2",0) == string::npos) cout << "Unknown format" << text << endl;
    
    //information about what things are in the file
    graspFile.getFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    text = buffer;
    cout << "[" << text << "]\n";
    hasVertices =  (text.find("vertices",0) != string::npos);
    hasAccessibles =  (text.find("accessibles",0) != string::npos);
    hasNormals =  (text.find("normals",0) != string::npos);
    hasTriangles =  (text.find("triangles",0) != string::npos);
    
    //information about what properties are in the file
    graspFile.getFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    text = buffer;
    cout << "[" << text << "]\n";
    hasPotential = (text.find("potential",0)  != string::npos);
    hasProperty1 = (text.find("gproperty1",0) != string::npos);
    hasProperty2 = (text.find("gproperty2",0) != string::npos);
    
    //Number of vertices and triangles:  Here use the fancy tokeniterator to break down line into strings
    graspFile.getFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    text = buffer;
    cout << "[" << text << "]\n";
    TokenIterator<char*, Delimiters>
    charIter(buffer, buffer + strlen(buffer), delimiters), end2;
    string token;
    token = *charIter++;
    int nVertices = atoi(token.c_str());
    token = *charIter++;
    nTriangles = atoi(token.c_str());
    
    vertices.resize(nVertices);
    triangles.resize(nTriangles);
    
    std::cout << "Number of vertices is " << vertices.size() << " in " << triangles.size() << " nTriangles\n";
    //	ostream_iterator<string> out(cout, "\n");
    //	copy (wordlist.begin(), wordlist.end(), out);
    
    // A line to skip
    graspFile.getFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    text = buffer;
    cout << "[" << text << "]\n";
    
    //Read vertices if present
    if (hasVertices){
        coordBuffer = new float[3*nVertices];
        double *coordDoubleBuffer = new double[3*nVertices];
        graspFile.getFortranData((char *) coordBuffer, sizeof(float),
                                 3*nVertices, CXXFortranFile::FortranFloatData);
        for (int i=0; i<3*nVertices; i++) coordDoubleBuffer[i] = coordBuffer[i];
        (void) addPerVertexVector ("vertices", coordDoubleBuffer);
        delete [] coordBuffer;
        delete [] coordDoubleBuffer;
    }
    
    //Read accessibles if present
    if (hasAccessibles){
        coordBuffer = new float[3*nVertices];
        double *coordDoubleBuffer = new double[3*nVertices];
        graspFile.getFortranData((char *) coordBuffer, sizeof(float),
                                 3*nVertices, CXXFortranFile::FortranFloatData);
        for (int i=0; i<3*nVertices; i++) coordDoubleBuffer[i] = coordBuffer[i];
        (void) addPerVertexVector ("accessibles", coordDoubleBuffer);
        delete [] coordBuffer;
        delete [] coordDoubleBuffer;
    }
    
    //Read normals if present
    if (hasNormals){
        coordBuffer = new float[3*nVertices];
        double *coordDoubleBuffer = new double[3*nVertices];
        graspFile.getFortranData((char *) coordBuffer, sizeof(float),
                                 3*nVertices, CXXFortranFile::FortranFloatData);
        for (int i=0; i<3*nVertices; i++) coordDoubleBuffer[i] = coordBuffer[i];
        (void) addPerVertexVector ("normals", coordDoubleBuffer);
        delete [] coordBuffer;
        delete [] coordDoubleBuffer;
    }
    //Read triangles if present
    if (hasTriangles){
        triangleBuffer = new int[3*nTriangles];
        graspFile.getFortranData((char *) triangleBuffer, sizeof(int),
                                 3*nTriangles, CXXFortranFile::FortranIntData);
        for (int i=0; i<nTriangles; i++){
            triangles[i] = CXXTriangle (triangleBuffer[3*i]-1,
                                        triangleBuffer[3*i+1]-1,
                                        triangleBuffer[3*i+2]-1);
        }
        delete [] triangleBuffer;
    }
    
    //Read potential if present
    if (hasPotential){
        float *potentialBuffer = new float[nVertices];
        double *potentialDoubleBuffer = new double[nVertices];
        graspFile.getFortranData((char *) potentialBuffer, sizeof(float),
                                 nVertices, CXXFortranFile::FortranFloatData);
        for (int i=0; i<nVertices; i++) potentialDoubleBuffer[i] = potentialBuffer[i];
        addPerVertexScalar("potential", potentialDoubleBuffer);
        delete [] potentialBuffer;
        delete [] potentialDoubleBuffer;
    }
    
    //Read gproperty1 if present
    if (hasProperty1){
        float *potentialBuffer = new float[nVertices];
        double *potentialDoubleBuffer = new double[nVertices];
        graspFile.getFortranData((char *) potentialBuffer, sizeof(float),
                                 nVertices, CXXFortranFile::FortranFloatData);
        for (int i=0; i<nVertices; i++) potentialDoubleBuffer[i] = potentialBuffer[i];
        addPerVertexScalar("gproperty1", potentialDoubleBuffer);
        delete [] potentialBuffer;
        delete [] potentialDoubleBuffer;
    }
    
    //Read gproperty2 if present
    if (hasProperty2){
        float *potentialBuffer = new float[nVertices];
        double *potentialDoubleBuffer = new double[nVertices];
        graspFile.getFortranData((char *) potentialBuffer, sizeof(float),
                                 nVertices, CXXFortranFile::FortranFloatData);
        for (int i=0; i<nVertices; i++) potentialDoubleBuffer[i] = potentialBuffer[i];
        addPerVertexScalar("gproperty2", potentialDoubleBuffer);
        delete [] potentialBuffer;
        delete [] potentialDoubleBuffer;
    }
    return 0;
}

SurfaceParameters CXXSurface::measuredProperties()
{
    double area = 0.;
    int i;
    CXXCoord<CXXCoord_ftype>AB, AC, ABxAC;
    SurfaceParameters result;
    
    result.nVertices = vertices.size();
    result.nTriangles = triangles.size();
    
    if (vectors.find("vertices") != vectors.end()){
        for (i=0, area = 0.; i<nTriangles; i++){
            const CXXCoord<CXXCoord_ftype>&coordA = coordRef(vectors["vertices"], i, 0);
            const CXXCoord<CXXCoord_ftype>&coordB = coordRef(vectors["vertices"], i, 1);
            const CXXCoord<CXXCoord_ftype>&coordC = coordRef(vectors["vertices"], i, 2);
            AB = coordA - coordB;
            AC = coordC - coordB;
            ABxAC = AB ^ AC;
            area += 0.5 * fabs(ABxAC.get3DLength());
        }
        result.MSA = area;
    }
    
    if (vectors.find("accessibles") != vectors.end()){
        for (i=0, area = 0.; i<nTriangles; i++){
            const CXXCoord<CXXCoord_ftype>&coordA = coordRef(vectors["accessibles"], i, 0);
            const CXXCoord<CXXCoord_ftype>&coordB = coordRef(vectors["accessibles"], i, 1);
            const CXXCoord<CXXCoord_ftype>&coordC = coordRef(vectors["accessibles"], i, 2);
            AB = coordA - coordB;
            AC = coordC - coordB;
            ABxAC = AB ^ AC;
            area += 0.5 * fabs(ABxAC.get3DLength());
        }
        result.ASA = area;
    }
    
    map<string, size_t>::iterator scalar;
    int j;
    for (scalar = scalars.begin(), j=0; scalar != scalars.end(); ++scalar, j++){
        double potential, pMin, pMax, pMean;
        for (i=0, pMin = 1e30, pMax = -1e30, pMean = 0.; i<int(vertices.size()); i++){
            potential = vertices[i].scalar(scalars[scalar->first]);
            pMin = (potential<pMin?potential:pMin);
            pMax = (potential>pMax?potential:pMax);
            pMean += potential;
        }
        pMean /= vertices.size();
        result.pMins[scalar->first] = pMin;
        result.pMins[scalar->first] = pMax;
        result.pMins[scalar->first] = pMean;
    }
    return result;
}

std::string CXXSurface::report(){
    double area = 0.;
    double pMin, pMax, pMean;
    int i;
    CXXCoord<CXXCoord_ftype>AB, AC, ABxAC;
    std::ostringstream output;
    
    output << "This surface has " << vertices.size() << " vertices and " << nTriangles
    << " triangles\n";
    
    if (vectors.find("vertices") != vectors.end()){
        for (i=0, area = 0.; i<nTriangles; i++){
            const CXXCoord<CXXCoord_ftype>&coordA = coordRef(vectors["vertices"], i, 0);
            const CXXCoord<CXXCoord_ftype>&coordB = coordRef(vectors["vertices"], i, 1);
            const CXXCoord<CXXCoord_ftype>&coordC = coordRef(vectors["vertices"], i, 2);
            AB = coordA - coordB;
            AC = coordC - coordB;
            ABxAC = AB ^ AC;
            area += 0.5 * fabs(ABxAC.get3DLength());
        }
        output << "The Molecular surface area is " << area << endl;
    }
    
    if (vectors.find("accessibles") != vectors.end()){
        for (i=0, area = 0.; i<nTriangles; i++){
            const CXXCoord<CXXCoord_ftype>&coordA = coordRef(vectors["accessibles"], i, 0);
            const CXXCoord<CXXCoord_ftype>&coordB = coordRef(vectors["accessibles"], i, 1);
            const CXXCoord<CXXCoord_ftype>&coordC = coordRef(vectors["accessibles"], i, 2);
            AB = coordA - coordB;
            AC = coordC - coordB;
            ABxAC = AB ^ AC;
            area += 0.5 * fabs(ABxAC.get3DLength());
        }
        output << "The Solvent Accessible surface area is " << area << endl;
    }
    
    map<string, size_t>::iterator scalar;
    int j;
    for (scalar = scalars.begin(), j=0; scalar != scalars.end(); ++scalar, j++){
        double potential;
        for (i=0, pMin = 1e30, pMax = -1e30, pMean = 0.; i<int(vertices.size()); i++){
            potential = vertices[i].scalar(scalars[scalar->first]);
            pMin = (potential<pMin?potential:pMin);
            pMax = (potential>pMax?potential:pMax);
            pMean += potential;
        }
        pMean /= vertices.size();
        output << "Property Number " << scalar->second //<< " name " << scalar->first
        << " Range " << pMin << " to " << pMax << " mean " << pMean << endl;
    }
    return output.str();
}

int CXXSurface::writeAsGrasp (string path)
{
    char buffer[1024];
    string text;
    int i;
    float *coordBuffer;
    int *triangleBuffer;
    std::cout << "writeAsGrasp" << path << "\n";
    
    Delimiters delimiters(" \0\t\n~;()\"<>:{}[]+-=&*#.,/\\");
    
    CXXFortranFile graspFile(path,"w");
    if (!graspFile.bad()) std::cout << "Opened grasp file [" << path << "]\n";
    else {
        std::cout << "Failed to open grasp file [" << path << "]\n";
        return 1;
    }
    
    //Grasp Surface file format
    strcpy (buffer,"format=2");
    for (size_t i=strlen(buffer); i<80; i++) buffer[i] = ' ';
    graspFile.putFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    
    //information about what things are in the file
    strcpy (buffer,"vertices");
    if (vectors.find("accessibles") != vectors.end()) strcat (buffer,",accessibles");
    if (vectors.find("normals") != vectors.end()) strcat (buffer,",normals");
    strcat (buffer,",triangles");
    for (size_t i=strlen(buffer); i<80; i++) buffer[i] = ' ';
    graspFile.putFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    
    //information about what properties are in the file
    int firstPotential = 1;
    buffer[0] = '\0';
    if (scalars.find("potential") != scalars.end()) {
        firstPotential = 0;
        strcat (buffer,"potential");
    }
    if (scalars.find("gproperty1") != scalars.end()) {
        if (!firstPotential) strcat(buffer,",");
        strcat (buffer,"gproperty1");
    }
    if (scalars.find("gproperty2") != scalars.end()) {
        if (!firstPotential) strcat(buffer,",");
        strcat (buffer,"gproperty2");
    }
    for (size_t i=strlen(buffer); i<80; i++) buffer[i] = ' ';
    graspFile.putFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    
    //Number of vertices and triangles, grid, and the number 3.0
    sprintf(buffer,"%8d%8zu%8d%8.5f                                ",
            int(vertices.size()), nTriangles, 65, 3.0);
    graspFile.putFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    
    // A line that seems to mention the origin
    sprintf(buffer,"%8.3f%8.3f%8.3f                                        ",
            0., 0., 0.);
    graspFile.putFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
    
    coordBuffer = new float[3*vertices.size()];
    
    //Output vertices
    for (i=0; i<int(vertices.size()); i++){
        CXXCoord_ftype *xyz = vertices[i].xyzPntr(vectors["vertices"]);
        for (int j=0; j<3; j++){
            coordBuffer[3*i + j] = xyz[j];
        }
    }
    graspFile.putFortranData((char *) coordBuffer, sizeof(float),
                             3*vertices.size(), CXXFortranFile::FortranFloatData);
    
    //Write accessibles if present
    if (vectors.find("accessibles") != vectors.end()){
        for (i=0; i<int(vertices.size()); i++){
            CXXCoord_ftype *xyz = vertices[i].xyzPntr(vectors["accessibles"]);
            for (int j=0; j<3; j++){
                coordBuffer[3*i + j] = xyz[j];
            }
        }
        graspFile.putFortranData((char *) coordBuffer, sizeof(float),
                                 3*vertices.size(), CXXFortranFile::FortranFloatData);
    }
    
    //Write normals if present
    if (vectors.find("normals") != vectors.end()){
        for (i=0; i<int(vertices.size()); i++){
            CXXCoord_ftype *xyz = vertices[i].xyzPntr(vectors["normals"]);
            for (int j=0; j<3; j++){
                coordBuffer[3*i + j] = xyz[j];
            }
        }
        graspFile.putFortranData((char *) coordBuffer, sizeof(float),
                                 3*vertices.size(), CXXFortranFile::FortranFloatData);
    }
    
    delete [] coordBuffer;
    
    
    //Write triangles
    if (triangles.size()>1){
        triangleBuffer = new int[3 * nTriangles];
        for (i=0; i<nTriangles; i++){
            for (int j=0; j<3; j++){
                triangleBuffer[3*i + j] = int(1 + triangles[i][j]);
            }
        }
        graspFile.putFortranData((char *) triangleBuffer, sizeof(int),
                                 3*nTriangles, CXXFortranFile::FortranIntData);
        delete [] triangleBuffer;
    }
    
    //Write potential if present
    if (scalars.find("potential")!=scalars.end()){
        size_t scalarHandle = getScalarHandle("potential");
        float *potentialBuffer =new float[vertices.size()];
        for (i=0; i<int(vertices.size()); i++){
            potentialBuffer[i] = vertices[i].scalar(scalarHandle);
        }
        graspFile.putFortranData((char *) potentialBuffer, sizeof(float),
                                 vertices.size(), CXXFortranFile::FortranFloatData);
        delete [] potentialBuffer;
    }
    
    //Write potential if present
    if (scalars.find("gproperty1")!=scalars.end()){
        float *potentialBuffer =new float[vertices.size()];
        for (i=0; i<int(vertices.size()); i++){
            potentialBuffer[i] = vertices[i].scalar(scalars["gproperty1"]);
        }
        graspFile.putFortranData((char *) potentialBuffer, sizeof(float),
                                 vertices.size(), CXXFortranFile::FortranFloatData);
        delete [] potentialBuffer;
    }
    
    //Write potential if present
    if (scalars.find("gproperty2")!=scalars.end()){
        float *potentialBuffer =new float[vertices.size()];
        for (i=0; i<int(vertices.size()); i++){
            potentialBuffer[i] = vertices[i].scalar(scalars["gproperty2"]);
        }
        graspFile.putFortranData((char *) potentialBuffer, sizeof(float),
                                 vertices.size(), CXXFortranFile::FortranFloatData);
        delete [] potentialBuffer;
    }
    return 0;
}

int CXXSurface::assignAtom (mmdb::Manager* allAtomsManager_in, int selHnd){
    void **pointerBuffer;
    mmdb::Atom** selAtom;
    int nSelAtoms;
    double minDistSq;
    
    allAtomsManager_in->GetSelIndex(selHnd, selAtom, nSelAtoms);
    pointerBuffer = new void* [vertices.size()];
    //	Loop over all vertices
    for (int i=0; i< int(vertices.size()); i++){
        const CXXCoord<CXXCoord_ftype>&vertex = coordRef(vectors["vertices"], i);
        int j;
        for (j = 0, minDistSq=1e30; j<nSelAtoms; j++){
            CXXCoord<CXXCoord_ftype>atom (selAtom[j]->x, selAtom[j]->y, selAtom[j]->z);
            CXXCoord<CXXCoord_ftype>diff = atom - vertex;
            double dxsq = diff.x() * diff.x();
            if (dxsq<minDistSq){
                double dysq = diff.y() * diff.y();
                if (dysq<minDistSq){
                    double dzsq = diff.z() * diff.z();
                    if (dzsq<minDistSq){
                        if (dxsq + dysq + dzsq < minDistSq){
                            pointerBuffer[i] = (void *) selAtom[j];
                            minDistSq = dxsq + dysq + dzsq;
                        }
                    }
                }
            }
        }
    }
    addPerVertexPointer ("atom", pointerBuffer);
    delete [] pointerBuffer;
    return 0;
}


int CXXSurface::colorByColourArray(const std::vector<double*> &colours, mmdb::Manager* molHnd, int selHnd){
    double greyColour[] = {0.5,0.5,0.5};
    double Colour[] = {0.5,0.5,0.5};
    
    size_t nVerts = vertices.size();
    double *colourBuffer = new double[3*nVerts];
    for (int i=0; i< nVerts; i++){
        for (int j=0; j<3; j++){
            colourBuffer[3*i + j] = greyColour[j];
        }
    }
    updateWithVectorData(vertices.size(), "colour", 0, colourBuffer);
    delete [] colourBuffer;
    
    int nSelAtoms;
    mmdb::Atom** SelAtom;
    molHnd->GetSelIndex(selHnd, SelAtom, nSelAtoms);
    
    int udd = molHnd->GetUDDHandle ( mmdb::UDR_ATOM,"tmp_atom_int" );
    if (udd <= 0 ) {
        udd = molHnd->RegisterUDInteger ( mmdb::UDR_ATOM,"tmp_atom_int" );
        if (udd <= 0 ) return udd;
    }
    for(int i=0; i<nSelAtoms; i++)
        SelAtom[i]->PutUDData(udd,i);
    
    for (int i=0; i< int(vertices.size()); i++){
        mmdb::Atom* theAtom = (mmdb::Atom*) vertices[i].pointer(pointers["atom"]);
        if (theAtom != 0){
            int atomid;
            theAtom->GetUDData(udd,atomid);
            Colour[0] = colours[atomid][0];
            Colour[1] = colours[atomid][1];
            Colour[2] = colours[atomid][2];
            vertices[i].setXyz(vectors["colour"], Colour);
        } else {
            vertices[i].setXyz(vectors["colour"], greyColour);
        }
    }
    return 0;
    
}

int CXXSurface::colorByAssignedAtom(){
    mmdb::Atom* theAtom;
    double redColour[] = {1.,0.,0.};
    double greenColour[] = {0.,1.,0.};
    double blueColour[] = {0.,0.,1.};
    double yellowColour[] = {1.,1.,0.};
    double greyColour[] = {1.,0.,0.};
    
    
    double *colourBuffer = new double[3*vertices.size()];
    for (int i=0; i< int(vertices.size()); i++){
        for (int j=0; j<3; j++){
            colourBuffer[3*i + j] = greyColour[j];
        }
    }
    updateWithVectorData(vertices.size(), "colour", 0, colourBuffer);
    delete [] colourBuffer;
    
    for (int i=0; i< int(vertices.size()); i++){
        theAtom = (mmdb::Atom*) vertices[i].pointer(pointers["atom"]);
        if (theAtom != 0){
            switch (theAtom->name[1]){
                case 'C':
                    vertices[i].setXyz(vectors["colour"], greenColour);
                    break;
                case 'N':
                    vertices[i].setXyz(vectors["colour"], blueColour);
                    break;
                case 'O':
                    vertices[i].setXyz(vectors["colour"], redColour);
                    break;
                case 'S':
                    vertices[i].setXyz(vectors["colour"], yellowColour);
                    break;
            }
        }
        else {
            vertices[i].setXyz(vectors["colour"], greyColour);
        }
    }
    return 0;
}

size_t CXXSurface::extendWithVectorData(size_t count, const string name, double *vectorBuffer){
    size_t start = vertices.size();
    updateWithVectorData( count, name,  start, vectorBuffer);
    return vertices.size();
}

size_t CXXSurface::updateWithVectorData(size_t count, const string name, size_t start, double *vectorBuffer){
    size_t iVector  = getVectorHandle(name);
    if (int(vertices.size()) < start + count){
        vertices.resize(start+count);
    }
    
    for (unsigned int i=0; i<(unsigned int) count; i++){
        vertices[i+start].setXyz(iVector, &(vectorBuffer[3*i]));
    }
    
    return vertices.size();
}


size_t CXXSurface::updateWithPointerData(size_t count, const string name, size_t start, void* *pointerBuffer){
    size_t iPointer = getPointerHandle(name);
    if (int(vertices.size()) < start + count){
        vertices.resize(start+count);
    }
    
    for (size_t i=0; i< count; i++){
        vertices[i+start].setPointer(iPointer, pointerBuffer[i]);
    }
    
    return vertices.size();
}

size_t CXXSurface::extendTriangles(int *triangleBuffer, int count){
    triangles.resize(nTriangles+count);
    for (int i=0; i<count; i++){
        int l = triangleBuffer[3*i];
        int m = triangleBuffer[3*i+1];
        int n = triangleBuffer[3*i+2];
        triangles[i+nTriangles] = CXXTriangle(l, m, n, i+nTriangles);
    }
    nTriangles = triangles.size();
    return nTriangles;
}

size_t CXXSurface::numberOfTriangles() const{
    return triangles.size();
}

size_t CXXSurface::numberOfVertices() const{
    return vertices.size();
}

size_t CXXSurface::vertex(size_t iTriangle, size_t iCorner) const{
    return triangles[iTriangle][iCorner];
}

int CXXSurface::upLoadSphere(CXXSphereElement &theSphere, double probeRadius, const int sense){
    
    size_t oldVertexCount;
    CXXCoord<CXXCoord_ftype>theCentre = theSphere.centre();
    
    vector<int> equivalence(theSphere.nVertices());
    vector<int> uniqueAndDrawn(theSphere.nVertices());
    
    int nDrawn = 0;
    for (unsigned  i=0; i< theSphere.nVertices(); i++){
        uniqueAndDrawn[i] = 0;
        if (theSphere.vertex(i).doDraw()){
            uniqueAndDrawn[i] = 1;
            if (uniqueAndDrawn[i]){
                equivalence[i] = nDrawn++;
            }
        }
    }
    static const std::string vertexName("vertices");
    static const std::string accessiblesName("accessibles");
    static const std::string normalsName("normals");
    {
        oldVertexCount = numberOfVertices();
        vertices.resize(oldVertexCount+nDrawn);
        size_t verticesHandle = getVectorHandle(vertexName);
        size_t accessiblesHandle = getVectorHandle(accessiblesName);
        size_t normalsHandle = getVectorHandle(normalsName);
        int iDraw = 0;
        for (unsigned int i=0; i< theSphere.nVertices(); i++){
            if (uniqueAndDrawn[i]){
                CXXCoord<CXXCoord_ftype>vertexCoord = theSphere.vertex(i).vertex();
                if (sense == CXXSphereElement::Contact){
                    vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord);
                    CXXCoord<CXXCoord_ftype>normal = vertexCoord - theCentre;
                    CXXCoord<CXXCoord_ftype>diff(normal);
                    diff *= (theSphere.radius() - probeRadius) / theSphere.radius();
                    normal.normalise();
                    vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
                    CXXCoord<CXXCoord_ftype>vertex = theCentre + diff;
                    vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertex);
                }
                else if (sense == CXXSphereElement::Reentrant) {
                    vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
                    CXXCoord<CXXCoord_ftype>normal = theCentre - vertexCoord;
                    normal.normalise();
                    vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
                    vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, theCentre);
                }
                else if (sense == CXXSphereElement::VDW) {
                    vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
                    CXXCoord<CXXCoord_ftype>normal = vertexCoord - theCentre;
                    normal.normalise();
                    vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
                    vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord+normal);
                }
                else if (sense == CXXSphereElement::Accessible) {
                    vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
                    CXXCoord<CXXCoord_ftype>normal = vertexCoord - theCentre;
                    normal.normalise();
                    vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
                    vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord);
                }
                iDraw++;
            }
        }
    }
    
    //Add atom pointers to the surface
    {
        std::vector<void*> atomBuffer(nDrawn);
        //void* *atomBuffer = new void*[nDrawn];
        int iDraw = 0;
        for (unsigned int i=0; i< theSphere.nVertices(); i++){
            if (uniqueAndDrawn[i]){
                mmdb::Atom* anAtom;
                if ((anAtom = theSphere.vertex(i).getAtom())!=0 ){
                    atomBuffer[iDraw] = (void*) anAtom;
                }
                else atomBuffer[iDraw] = (void*) theSphere.getAtom();
                iDraw++;
            }
        }
        updateWithPointerData(nDrawn, "atom", oldVertexCount, (void **)&(atomBuffer[0]));
        //delete [] atomBuffer;
    }
    
    // Add triangles to surface
    {
        std::vector<int> triangleBuffer(theSphere.nFlatTriangles()*3);
        int drawCount = 0;
        std::list<CXXSphereFlatTriangle  >::const_iterator trianglesEnd =
        theSphere.getFlatTriangles().end();
        for (std::list<CXXSphereFlatTriangle  >::const_iterator triangle =
             theSphere.getFlatTriangles().begin();
             triangle != trianglesEnd;
             ++triangle){
            const CXXSphereFlatTriangle &theTriangle(*triangle);
            if (theTriangle.doDraw()){
                if (sense == CXXSphereElement::Contact ||
                    sense == CXXSphereElement::VDW ||
                    sense == CXXSphereElement::Accessible){
                    for (unsigned int j=0; j<3; j++){
                        int index = equivalence[theTriangle[2-j]];
                        triangleBuffer[3*drawCount+j] = int(index + oldVertexCount);
                    }
                }
                else {
                    for (unsigned int j=0; j<3; j++){
                        int index = equivalence[theTriangle[j]];
                        triangleBuffer[3*drawCount+j] = int(index + oldVertexCount);
                    }
                }
                drawCount++;
            }
        }
        extendTriangles((int *) &(triangleBuffer[0]), drawCount);
    }
    /* */
    
    
    
    return 0;
}

int CXXSurface::uploadTorus(CXXTorusElement &theTorus) {
    //	Add vertices to surface
    size_t oldVertexCount;
    {
        double verticesBuffer[theTorus.nTorusNodes()*3];// = new double[nodes.size()*3];
        double accessiblesBuffer[theTorus.nTorusNodes()*3];// = new double[nodes.size()*3];
        double normalsBuffer[theTorus.nTorusNodes()*3];// = new double[nodes.size()*3];
        for (unsigned int i=0; i< theTorus.nTorusNodes(); i++){
            for (int j=0; j<3; j++) verticesBuffer[3*i+j] = theTorus.node(i).coord().element(j);
            CXXCoord<CXXCoord_ftype>accessible = theTorus.probeAtOmega(theTorus.node(i).getOmega());
            for (int j=0; j<3; j++) accessiblesBuffer[3*i+j] = accessible.element(j);
            CXXCoord<CXXCoord_ftype>normal = theTorus.normalToProbeAtTheta(accessible, theTorus.node(i).getTheta());
            normal.scale(-1.);
            for (int j=0; j<3; j++) normalsBuffer[3*i+j] = normal.element(j);
        }
        oldVertexCount = numberOfVertices();
        updateWithVectorData(theTorus.nTorusNodes(), "vertices", oldVertexCount, verticesBuffer);
        updateWithVectorData(theTorus.nTorusNodes(), "accessibles", oldVertexCount, accessiblesBuffer);
        updateWithVectorData(theTorus.nTorusNodes(), "normals",  oldVertexCount, normalsBuffer);
    }
    //Add atom pointers to the surface
    {
        void *atomBuffer[theTorus.nTorusNodes()];// = new void*[nodes.size()];
        for (unsigned int i=0; i< theTorus.nTorusNodes(); i++){
            atomBuffer[i] = (void *)theTorus.node(i).getAtom();
        }
        updateWithPointerData(theTorus.nTorusNodes(), "atom", oldVertexCount, atomBuffer);
    }
    // Add triangles to surface
    {
        int triangleBuffer[theTorus.nFlatTriangles()*3];// = new int[flatTriangles.size()*3];
        int nToDraw = 0;
        ;
        for (list <CXXTriangle  >::const_iterator triangle = theTorus.firstTriangle();
             triangle != theTorus.endOfTriangles();
             ++triangle){
            const CXXTriangle &flatTriangle(*triangle);
            if (flatTriangle.doDraw()){
                for (unsigned int j=0; j<3; j++){
                    //Note the 2-j, this changes the sense of the triangle to reflect the fact that we
                    //actually visualise the inside of the torus
                    triangleBuffer[(3*nToDraw)+j] = int(flatTriangle[j] + oldVertexCount);
                }
                nToDraw++;
            }
        }
        extendTriangles(triangleBuffer, nToDraw);
    }
    
    //});
    return 0;
}

mmdb::Manager* CXXSurface::getMMDBManager() const{
    return allAtomsManager;
}

size_t CXXSurface::setCoord(const string &name, size_t iVertex, const CXXCoord<CXXCoord_ftype>&crd){
    size_t iVector = vectors[name];
    if (int(vertices.size()) <= iVertex){
        vertices.resize(iVertex+1);
    }
    
    vertices[iVertex].setCoord(iVector, crd);
    
    return vertices.size();
}

int CXXSurface::getIntegerUDDataOfAtom(mmdb::Atom* theAtom, int handle){
    int result;
    int rc = theAtom->GetUDData(handle, result);
    switch (rc)  {
            
        case  mmdb::UDDATA_WrongUDRType :
            printf ( " wrong UDD registration type\n" );
            break;
            
        case  mmdb::UDDATA_WrongHandle  :
            printf ( " wrong UDD handle\n" );
            break;
            
        case  mmdb::UDDATA_NoData :
            printf ( " UDD not found.\n" );
            break;
            
        case  mmdb::UDDATA_Ok :
            break;
    }
    
    return result;
}

size_t CXXSurface::getVectorHandle (const string name){
    if (vectors.find(name) == vectors.end()){
        size_t oldSize = vectors.size();
        vectors[name] = oldSize + 1;
    }
    return vectors[name];
}

size_t CXXSurface::getReadVectorHandle (const string name){
    if (vectors.find(name) == vectors.end()){
        return SIZE_MAX;
    }
    return vectors[name];
}


size_t CXXSurface::getScalarHandle (const string name){
    if (scalars.find(name) == scalars.end()){
        size_t oldSize = scalars.size();
        scalars[name] = oldSize + 1;
    }
    return scalars[name];
}

void CXXSurface::setScalar(size_t scalarHandle, size_t iVertex, double &value){
    vertices[iVertex].setScalar(scalarHandle, value);
}

void CXXSurface::setScalar(const string name, size_t iVertex, double &value){
    size_t scalarHandle=getScalarHandle(name);
    vertices[iVertex].setScalar(scalarHandle, value);
}


size_t CXXSurface::getReadScalarHandle (const string name){
    if (scalars.find(name) == scalars.end()){
        return SIZE_MAX;
    }
    return scalars[name];
}

size_t CXXSurface::getPointerHandle (const string name){
    if (pointers.find(name) == pointers.end()){
        size_t oldSize = pointers.size();
        pointers[name] = oldSize + 1;
    }
    return pointers[name];
}

const CXXCoord<CXXCoord_ftype>&CXXSurface::coordRef(size_t coordType, size_t iTriangle, size_t iCorner) const{
    const CXXTriangle &theTriangle(triangles[iTriangle]);
    size_t iVertex = theTriangle[iCorner];
    return vertices[iVertex].coordRef(coordType);
}

const CXXCoord<CXXCoord_ftype>&CXXSurface::coordRef(size_t coordType, size_t iVertex) const{
    return vertices[iVertex].coordRef(coordType);
}

int CXXSurface::addPerVertexVector (const string name, double *coordBuffer){
    size_t vectorHandle = getVectorHandle(name);
    for (int i=0; i<int(vertices.size()); i++){
        vertices[i].setXyz(vectorHandle, &(coordBuffer[3*i]));
    }
    return 0;
}

int CXXSurface::addPerVertexScalar (const string name, double *scalarBuffer){
    size_t scalarHandle = getScalarHandle(name);
    for (int i=0; i<int(vertices.size()); i++){
        vertices[i].setScalar(scalarHandle, scalarBuffer[i]);
    }
    return 0;
}

int CXXSurface::addPerVertexPointer (const string name, void **pointerBuffer){
    size_t pointerHandle = getPointerHandle(name);
    for (int i=0; i<int(vertices.size()); i++){
        vertices[i].setPointer(pointerHandle, pointerBuffer[i]);
    }
    return 0;
}

int CXXSurface::getCoord(const string &type, const size_t iTriangle, const size_t corner, double *buffer)
{
    size_t relevantVector = getReadVectorHandle(type);
    if (relevantVector==SIZE_MAX) return 1;
    size_t iVertex = triangles[iTriangle][corner];
    return getCoord(relevantVector, iVertex, buffer);
}

int CXXSurface::getCoord(const string &type, const size_t iVertex, double *buffer)
{
    size_t relevantVector = getReadVectorHandle(type);
    if (relevantVector==SIZE_MAX) return 1;
    return getCoord(relevantVector, iVertex, buffer);
}

int CXXSurface::getCoord(const size_t type, const size_t iVertex, double *buffer)
{
    if (type <= int(vectors.size())){
        const CXXCoord<CXXCoord_ftype>&theCoord = coordRef(type, iVertex);
        for (int i=0; i<4; i++){
            buffer[i] = theCoord.element(i);
        }
        return 0;
    }
    else return 1;
}

int CXXSurface::getPointer(const string &type, size_t iVertex, void **return_p)
{
    if (pointers.find(type) != pointers.end()){
        *return_p = vertices[iVertex].pointer(pointers[type]);
        return 0;
    }
    *return_p = 0;
    return 1;
}

int CXXSurface::getScalar  (int handle, int iVertex, double &result){
    result = vertices[iVertex].scalar(handle);
    return 0;
}

int CXXSurface::addTriangle (const CXXTriangle &aTriangle){
    triangles.push_back(aTriangle);
    return 0;
}

void CXXSurface::appendSurface(const CXXSurface &otherSurface){
    size_t oldNVertices = vertices.size();
    size_t oldNTriangles = triangles.size();
    
    if (vectors.size() == 0) vectors = otherSurface.vectorNames();
    if (scalars.size() == 0) scalars = otherSurface.scalarNames();
    if (pointers.size() == 0) pointers = otherSurface.pointerNames();
    
    vertices.insert(vertices.end(),otherSurface.getVertices().begin(), otherSurface.getVertices().end());
    triangles.insert(triangles.end(),otherSurface.getTriangles().begin(), otherSurface.getTriangles().end());
    nTriangles = triangles.size();
    
    vector<CXXTriangle  >::iterator triangle;
    vector<CXXTriangle  >::iterator triangleEnd = triangles.end();
    for (triangle = (triangles.begin() + oldNTriangles); triangle!=triangleEnd; ++triangle){
        for (int i=0; i<3; i++){
            (*triangle)[i] = (*triangle)[i] + oldNVertices;
        }
    }
}

void CXXSurface::compress(double tolerance){
    vector<CXXSurfaceVertex> compressedVertices;
    compressedVertices.reserve(vertices.size());
    vector<CXXTriangle  >compressedTriangles;
    compressedTriangles.reserve(triangles.size());
    
    size_t vertexHandle = getVectorHandle("vertices");
    size_t normalHandle = getVectorHandle("normals");
    size_t atomHandle = getVectorHandle("atom");
    
    vector<size_t>equivalences(vertices.size());
    
    for (int i=0; i<vertices.size(); i++){
        bool uniqueAndDrawn = true;
        for (int j=0; j<compressedVertices.size() && uniqueAndDrawn; j++){
            if ((vertices[i].coordRef(vertexHandle).isNearly
                 (compressedVertices[j].coordRef(vertexHandle), tolerance) &&
                 vertices[i].coordRef(normalHandle).isNearly
                 (compressedVertices[j].coordRef(normalHandle), tolerance)) &&
                (compressedVertices[j].pointer(atomHandle) == vertices[i].pointer(atomHandle)) ){
                uniqueAndDrawn = false;
                equivalences[i] = j;
            }
        }
        if (uniqueAndDrawn){
            compressedVertices.push_back(vertices[i]);
            equivalences[i] = compressedVertices.size()-1;
        }
    }
    size_t equivalent[3];
    for (int i=0; i<triangles.size(); i++){
        for (int j=0; j<3; j++){
            equivalent[j] = equivalences[vertex(i,j)];
        }
        if ((equivalent[0] != equivalent[1] &&
             equivalent[0] != equivalent[2] &&
             equivalent[1] != equivalent[2])){
            triangles[i][0] = equivalent[0];
            triangles[i][1] = equivalent[1];
            triangles[i][2] = equivalent[2];
            compressedTriangles.push_back(triangles[i]);
        }
    }
    vertices=compressedVertices;
    triangles=compressedTriangles;
    nTriangles = triangles.size();
};

