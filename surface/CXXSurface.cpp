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
#include "CXXSphereElement.h"
#include <mmdb2/mmdb_tables.h>
#include <mmdb2/mmdb_uddata.h>
#include "CXXCircle.h"
#include "CXXTriangle.h"
#include "CXXSphereFlatTriangle.h"

CXX_mot::CXXSurface::CXXSurface()
{
	init();
}

CXX_mot::CXXSurface::~CXXSurface(){
	triangles.resize(0);
	vertices.resize(0);
}

int CXX_mot::CXXSurface::init()
{
	nTriangles = 0;
	triangles.resize(0);
	vertices.resize(0);
	return (0);
}

CXX_mot::CXXSurface::CXXSurface (string path)
{
	init();
	readGraspFile (path);
}

int CXX_mot::CXXSurface::readGraspFile(string path)
{
	char buffer[1024];
	string text;
	float *coordBuffer;
	int *triangleBuffer;
	
	Delimiters delimiters(" \0\t\n~;()\"<>:{}[]+-=&*#.,/\\");	
	int hasVertices, hasAccessibles, hasNormals, hasTriangles, hasPotential, 
    hasProperty1, hasProperty2;	
	hasVertices= hasAccessibles= hasNormals= hasTriangles= hasPotential= 
    hasProperty1= hasProperty2 = 0;
	
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

std::string CXX_mot::CXXSurface::report(){
	double area = 0.;
	double pMin, pMax, pMean;
	int i;
	CXXCoord AB, AC, ABxAC;
	std::ostringstream output;
	
	output << "This surface has " << vertices.size() << " vertices and " << nTriangles 
    << " triangles\n";
	
	if (vectors.find("vertices") != vectors.end()){
		for (i=0, area = 0.; i<nTriangles; i++){
			const CXXCoord &coordA = coordRef(vectors["vertices"], i, 0);
			const CXXCoord &coordB = coordRef(vectors["vertices"], i, 1);
			const CXXCoord &coordC = coordRef(vectors["vertices"], i, 2);
			AB = coordA - coordB;
			AC = coordC - coordB;
			ABxAC = AB ^ AC;
			area += 0.5 * fabs(ABxAC.get3DLength());
		}
		output << "The Molecular surface area is " << area << endl;
	}
	
	if (vectors.find("accessibles") != vectors.end()){
		for (i=0, area = 0.; i<nTriangles; i++){
			const CXXCoord &coordA = coordRef(vectors["accessibles"], i, 0);
			const CXXCoord &coordB = coordRef(vectors["accessibles"], i, 1);
			const CXXCoord &coordC = coordRef(vectors["accessibles"], i, 2);
			AB = coordA - coordB;
			AC = coordC - coordB;
			ABxAC = AB ^ AC;
			area += 0.5 * fabs(ABxAC.get3DLength());
		}
		output << "The Solvent Accessible surface area is " << area << endl;
	}
	
	map<string, int>::iterator scalar;
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

int CXX_mot::CXXSurface::writeAsGrasp (const std::string &path)
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
	for (i=strlen(buffer); i<80; i++) buffer[i] = ' ';
	graspFile.putFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
	
	//information about what things are in the file
	strcpy (buffer,"vertices");
	if (vectors.find("accessibles") != vectors.end()) strcat (buffer,",accessibles");
	if (vectors.find("normals") != vectors.end()) strcat (buffer,",normals");
	strcat (buffer,",triangles");
	for (i=strlen(buffer); i<80; i++) buffer[i] = ' ';
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
	for (i=strlen(buffer); i<80; i++) buffer[i] = ' ';
	graspFile.putFortranData(buffer, 1, 64, CXXFortranFile::FortranCharData);
	
	//Number of vertices and triangles, grid, and the number 3.0
	sprintf(buffer,"%8d%8d%8d%8.5f                                ",
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
				triangleBuffer[3*i + j] = 1 + triangles[i][j];
			}
		}
		graspFile.putFortranData((char *) triangleBuffer, sizeof(int), 
								 3*nTriangles, CXXFortranFile::FortranIntData);
		delete [] triangleBuffer;
	}
	
	//Write potential if present
	if (scalars.find("potential")!=scalars.end()){
		int scalarHandle = getScalarHandle("potential");
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

CXX_mot::CXXSurface::CXXSurface (mmdb::PManager allAtomsManager_in, const std::string selectionString,
						const std::string contextString){
	double delta = 30. * 2. * M_PI / 360.;
	double probeRadius = 1.4;
    bool blend_edges = false;
	init();
	calculateFromAtoms( allAtomsManager_in,  selectionString,  contextString, probeRadius, delta,blend_edges);
}

CXX_mot::CXXSurface::CXXSurface (mmdb::PManager allAtomsManager_in, const std::string selectionString,
                        const std::string contextString, const double delta, const double probeRadius,  const bool blend_edges){
	init();
	calculateFromAtoms( allAtomsManager_in,  selectionString,  contextString, probeRadius, delta, blend_edges);
}

CXX_mot::CXXSurface::CXXSurface (mmdb::PManager allAtomsManager_in, const std::string selectionString){
	double delta = 30. * 2. * M_PI / 360.;
	double probeRadius = 1.4;
    bool blend_edges = false;
	calculateFromAtoms( allAtomsManager_in,  selectionString, probeRadius, delta,blend_edges);
}

CXX_mot::CXXSurface::CXXSurface (mmdb::PManager allAtomsManager_in, const int selHnd){
	double delta = 30. * 2. * M_PI / 360.;
	double probeRadius = 1.4;
    bool blend_edges = false;
	init();
	calculateFromAtoms( allAtomsManager_in,  selHnd,  selHnd, probeRadius, delta,blend_edges);
}

CXX_mot::CXXSurface::CXXSurface (mmdb::PManager allAtomsManager_in, const int selHnd, const double delta, const double probeRadius, const bool blend_edges){
	init();
	calculateFromAtoms( allAtomsManager_in,  selHnd,  selHnd, probeRadius, delta,blend_edges);
}

CXX_mot::CXXSurface::CXXSurface (mmdb::PManager allAtomsManager_in, const int selHnd, const int contextSelHnd, 
                        const double delta, const double probeRadius,const bool blend_edges ){
	init();
	calculateFromAtoms( allAtomsManager_in,  selHnd,  contextSelHnd, probeRadius, delta,blend_edges);
}

int CXX_mot::CXXSurface::calculateFromAtoms(mmdb::PManager allAtomsManager_in, const std::string selectionString, const double probeRadius, const double delta, const bool blend_edges){
	int selHnd = selectionStringToSelHnd(allAtomsManager_in, selectionString);
	return calculateFromAtoms(allAtomsManager_in, selHnd, selHnd, probeRadius, delta,blend_edges);
}

int CXX_mot::CXXSurface::calculateFromAtoms(mmdb::PManager allAtomsManager_in, const std::string selectionString, const std::string contextString, const double probeRadius, const double delta , const bool blend_edges ){
	int selHnd = selectionStringToSelHnd(allAtomsManager_in, selectionString);
	int contextHnd = selectionStringToSelHnd(allAtomsManager_in, contextString);
	return calculateFromAtoms(allAtomsManager_in, selHnd, contextHnd, probeRadius, delta,blend_edges );
}

int CXX_mot::CXXSurface::calculateFromAtoms(mmdb::PManager allAtomsManager_in, const int selHnd, const double probeRadius, const double delta, const bool blend_edges){
	return calculateFromAtoms(allAtomsManager_in, selHnd, selHnd, probeRadius, delta,blend_edges);
}

int CXX_mot::CXXSurface::calculateVDWFromAtoms(mmdb::PManager allAtomsManager_in, const int selHnd, const int contextSelHnd, const double probeRadius, const double delta , const bool blend_edges){
	allAtomsManager = allAtomsManager_in;	
	
	int  nSelAtoms;
	mmdb::PPAtom SelAtom;	
	allAtomsManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	cout << "Surface selection includes " << nSelAtoms << "atoms"<<endl;    
	vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall*> > vdwBallPntrs;
	for (int atomNr = 0;atomNr < nSelAtoms; atomNr++) { 
		vdwBallPntrs.push_back(new CXXAtomBall(SelAtom[atomNr], getAtomRadius(SelAtom[atomNr])));
	}
	map<mmdb::PAtom, const CXXBall *> mainAtoms;
	for (int atomNr = 0;atomNr < nSelAtoms; atomNr++) { 
		mainAtoms[vdwBallPntrs[atomNr]->getAtomI()] = vdwBallPntrs[atomNr]; 
	}
	
	int  nContextSelAtoms;
	mmdb::PPAtom ContextSelAtom;	
	allAtomsManager->GetSelIndex(contextSelHnd, ContextSelAtom, nContextSelAtoms);
	cout << "Context selection includes " << nSelAtoms << "atoms"<<endl; 
	
	vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall*> > contextBallPntrs;
	int nUniqueContextAtoms = 0;
	int nSharedContextAtoms = 0;
	for (int atomNr = 0;atomNr < nContextSelAtoms; atomNr++) { 
		map<mmdb::PAtom, const CXXBall *>::iterator equivalentMainAtom = mainAtoms.find(ContextSelAtom[atomNr]);
		if (equivalentMainAtom != mainAtoms.end()){
			contextBallPntrs.push_back(equivalentMainAtom->second);
			nSharedContextAtoms++;
		}
		else{
			contextBallPntrs.push_back(new CXXAtomBall(ContextSelAtom[atomNr], getAtomRadius(ContextSelAtom[atomNr])));
			nUniqueContextAtoms++;
		}
	}
	std::cout << "nUniqueContextAtoms " << nUniqueContextAtoms << " nSharedContextAtoms " << nSharedContextAtoms << std::endl;
	
	CXXBall::triangulateBalls(vdwBallPntrs, contextBallPntrs, delta, this, CXXSphereElement::VDW);
	for (unsigned int i=0; i<vdwBallPntrs.size(); i++){
		if (vdwBallPntrs[i]) delete vdwBallPntrs[i];
	}
	for (unsigned int i=0; i<contextBallPntrs.size(); i++){
		map<mmdb::PAtom, const CXXBall *>::iterator equivalentMainAtom = mainAtoms.find(ContextSelAtom[i]);
		if (equivalentMainAtom == mainAtoms.end()){
			if (contextBallPntrs[i]) delete contextBallPntrs[i];
		}
	}
	report();
	return 0;
	
} 

int CXX_mot::CXXSurface::calculateAccessibleFromAtoms(mmdb::PManager allAtomsManager_in, const int selHnd, const int contextSelHnd, const double probeRadius, const double delta , const bool blend_edges){
	allAtomsManager = allAtomsManager_in;	
	
	int  nSelAtoms;
	mmdb::PPAtom SelAtom;	
	allAtomsManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	cout << "Surface selection includes " << nSelAtoms << "atoms"<<endl;    
	vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall*> > vdwBallPntrs;
	for (int atomNr = 0;atomNr < nSelAtoms; atomNr++) { 
		vdwBallPntrs.push_back(new CXXAtomBall(SelAtom[atomNr], probeRadius+getAtomRadius(SelAtom[atomNr])));
	}
	map<mmdb::PAtom, const CXXBall *> mainAtoms;
	for (int atomNr = 0;atomNr < nSelAtoms; atomNr++) { 
		mainAtoms[vdwBallPntrs[atomNr]->getAtomI()] = vdwBallPntrs[atomNr]; 
	}

	int  nContextSelAtoms;
	mmdb::PPAtom ContextSelAtom;	
	allAtomsManager->GetSelIndex(contextSelHnd, ContextSelAtom, nContextSelAtoms);
	cout << "Context selection includes " << nSelAtoms << "atoms"<<endl; 
	
	vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall*> > contextBallPntrs;
	int nUniqueContextAtoms = 0;
	int nSharedContextAtoms = 0;
	for (int atomNr = 0;atomNr < nContextSelAtoms; atomNr++) { 
		map<mmdb::PAtom, const CXXBall *>::iterator equivalentMainAtom = mainAtoms.find(ContextSelAtom[atomNr]);
		if (equivalentMainAtom != mainAtoms.end()){
			contextBallPntrs.push_back(equivalentMainAtom->second);
			nSharedContextAtoms++;
		}
		else{
			contextBallPntrs.push_back(new CXXAtomBall(ContextSelAtom[atomNr], probeRadius+getAtomRadius(ContextSelAtom[atomNr])));
			nUniqueContextAtoms++;
		}
	}
	std::cout << "nUniqueContextAtoms " << nUniqueContextAtoms << " nSharedContextAtoms " << nSharedContextAtoms << std::endl;
	
	CXXBall::triangulateBalls(vdwBallPntrs, contextBallPntrs, delta, this, CXXSphereElement::Accessible);
	for (unsigned int i=0; i<vdwBallPntrs.size(); i++){
		if (vdwBallPntrs[i]) delete vdwBallPntrs[i];
	}
	for (unsigned int i=0; i<contextBallPntrs.size(); i++){
		map<mmdb::PAtom, const CXXBall *>::iterator equivalentMainAtom = mainAtoms.find(ContextSelAtom[i]);
		if (equivalentMainAtom == mainAtoms.end()){
			if (contextBallPntrs[i]) delete contextBallPntrs[i];
		}
	}
	report();
	return 0;
} 

int CXX_mot::CXXSurface::calculateFromAtoms(mmdb::PManager allAtomsManager_in, const int selHnd, const int contextSelHnd, const double probeRadius, const double delta , const bool blend_edges){
	allAtomsManager = allAtomsManager_in;	
	
	int  nSelAtoms;
	mmdb::PPAtom SelAtom;	
	allAtomsManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	cout << "Surface selection includes " << nSelAtoms << "atoms"<<endl;    
	vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall*> > vdwBallPntrs;
	for (int atomNr = 0;atomNr < nSelAtoms; atomNr++) { 
		vdwBallPntrs.push_back(new CXXAtomBall(SelAtom[atomNr], probeRadius+getAtomRadius(SelAtom[atomNr])));
	}
	int  nContextSelAtoms;
	mmdb::PPAtom ContextSelAtom;	
	allAtomsManager->GetSelIndex(contextSelHnd, ContextSelAtom, nContextSelAtoms);
	cout << "Context selection includes " << nSelAtoms << "atoms"<<endl;    
	vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall*> > contextBallPntrs;
	for (int atomNr = 0;atomNr < nContextSelAtoms; atomNr++) { 
		contextBallPntrs.push_back(new CXXAtomBall(ContextSelAtom[atomNr], probeRadius+getAtomRadius(ContextSelAtom[atomNr])));
	}
	
	//Precalculate contacts
    std::map<const CXXBall *, std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> > >contactMap;
    CXXBall::ballContacts(vdwBallPntrs, contextBallPntrs, contactMap);	
    std::cout << "Established contact map\n";
	
    //Start up with reentrant Prbes separated per-atom...removes one locking state
    //for multi-threadig
    std::vector<std::vector<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> > > splitReentrantProbes;
    splitReentrantProbes.resize(nSelAtoms);
	
    CXXSphereElement unitSphereAtOrigin(CXXCoord(0.,0.,0.), 1., delta);
	

#pragma omp parallel for default(none) shared(selHnd, nSelAtoms, SelAtom, probeRadius, delta, cout, unitSphereAtOrigin, vdwBallPntrs, contactMap, ContextSelAtom, splitReentrantProbes) schedule(dynamic, 100) //num_threads(2)
    // #pragma omp parallel for default(none) shared(nSelAtoms, SelAtom, cout, unitSphereAtOrigin, vdwBallPntrs, contactMap, ContextSelAtom, splitReentrantProbes) schedule(dynamic, 100) //num_threads(2)
	for (int atomNr = 0;atomNr < nSelAtoms; atomNr++) { 
		mmdb::PAtom centralAtom = static_cast<const CXXAtomBall *>(vdwBallPntrs[atomNr])->getAtomI();		
		if (!(atomNr%100) || atomNr==nSelAtoms-1) {
#pragma omp critical(cout)
			cout << "Dealing with atom number " <<atomNr <<endl;
		}
		double radiusOfAtom1 = getAtomRadius(centralAtom);
        CXXNewHood theNewHood;
        theNewHood.initWith(centralAtom, radiusOfAtom1, probeRadius);
        
		//We have precalculated neighbours of the central atom, and now can use that 
		//to our advantage
		std::vector<const CXXBall *, CXX_old::CXXAlloc<const CXXBall *> > &neighbours = contactMap[vdwBallPntrs[atomNr]];
		for (unsigned int sphereAtomNr = 0; sphereAtomNr < neighbours.size(); sphereAtomNr++) {
			theNewHood.addBall(*neighbours[sphereAtomNr]);
		}
		
		//Find the non-hidden segments of the circles
		theNewHood.findSegments();
        CXXSurface elementSurface;
		if (!CXXNewHood::doesNotContainDrawable(theNewHood)){
            theNewHood.triangulateAsRegularHoodInto(&elementSurface, delta, &unitSphereAtOrigin);
            theNewHood.identifyUniqueNodes(splitReentrantProbes[atomNr], selHnd);
            elementSurface.compress(0.00001);
#pragma omp critical (mainTriangles)
            appendSurface(elementSurface);
		}
	}

	vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > reentrantProbes;
    for (int i=0; i<nSelAtoms; i++){
        std::vector<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator reentrantProbesEnd(splitReentrantProbes[i].end());
        for (std::vector<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator reentrantProbe = splitReentrantProbes[i].begin();
             reentrantProbe != reentrantProbesEnd;
             ++reentrantProbe){
            reentrantProbes.push_back(new CXXReentrantProbeBall(*reentrantProbe, selHnd, probeRadius));
        }
    }
    splitReentrantProbes.resize(0);
    for (unsigned int i=0; i<vdwBallPntrs.size(); i++){
		if (vdwBallPntrs[i]) delete static_cast<const CXXAtomBall *>(vdwBallPntrs[i]);
    }
    for (unsigned int i=0; i<contextBallPntrs.size(); i++){
		if (contextBallPntrs[i]) delete static_cast<const CXXAtomBall *>(contextBallPntrs[i]);
    }
    CXXBall::triangulateBalls(reentrantProbes, reentrantProbes, delta, this, CXXSphereElement::Reentrant);
    for (unsigned int i=0; i<reentrantProbes.size(); i++){
       delete static_cast<const CXXReentrantProbeBall *>(reentrantProbes[i]);
    }
	
    if (blend_edges) {
        cout << "Starting to blend edges" <<endl;
        assignAtom (allAtomsManager,selHnd);
    }
	report();
    //cout << "Starting default colour surface" <<endl;
    //colorByAssignedAtom();
    //cout << "Finished default colour surface" <<endl;
    
    return 0;
}

int CXX_mot::CXXSurface::assignAtom (mmdb::PManager allAtomsManager_in, int selHnd){
    void **pointerBuffer;
    mmdb::PPAtom selAtom;
    int nSelAtoms;
    double minDistSq;
    
    allAtomsManager_in->GetSelIndex(selHnd, selAtom, nSelAtoms);
    pointerBuffer = new void* [vertices.size()];
    //	Loop over all vertices
    for (int i=0; i< int(vertices.size()); i++){
        const CXXCoord &vertex = coordRef(vectors["vertices"], i);
        int j;
        for (j = 0, minDistSq=1e30; j<nSelAtoms; j++){
            CXXCoord atom (selAtom[j]->x, selAtom[j]->y, selAtom[j]->z);
            CXXCoord diff = atom - vertex;
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


int CXX_mot::CXXSurface::colorByColourArray(const std::vector<double*> &colours, mmdb::Manager *molHnd, int selHnd){
    double greyColour[] = {0.5,0.5,0.5};
    double Colour[] = {0.5,0.5,0.5};
    
    int nVerts = vertices.size(); 
    double *colourBuffer = new double[3*nVerts];
    for (int i=0; i< nVerts; i++){
        for (int j=0; j<3; j++){
            colourBuffer[3*i + j] = greyColour[j];
        }
    }
    updateWithVectorData(vertices.size(), "colour", 0, colourBuffer);
    delete [] colourBuffer;
    
    int nSelAtoms;
    mmdb::PPAtom SelAtom;
    molHnd->GetSelIndex(selHnd, SelAtom, nSelAtoms);
    
    int udd = molHnd->GetUDDHandle ( mmdb::UDR_ATOM,"tmp_atom_int" );
    if (udd <= 0 ) {
        udd = -1;
        udd = molHnd->RegisterUDInteger ( mmdb::UDR_ATOM,"tmp_atom_int" );
        if (udd <= 0 ) return udd;
    }
    for(int i=0; i<nSelAtoms; i++)
        SelAtom[i]->PutUDData(udd,i);
    
    for (int i=0; i< int(vertices.size()); i++){
        mmdb::PAtom theAtom = (mmdb::PAtom) vertices[i].pointer(pointers["atom"]);
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

int CXX_mot::CXXSurface::colorByAssignedAtom(){
    mmdb::PAtom theAtom;
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
        theAtom = (mmdb::PAtom) vertices[i].pointer(pointers["atom"]);
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

int CXX_mot::CXXSurface::extendWithVectorData(int count, const string name, double *vectorBuffer){
    int start = vertices.size();
    updateWithVectorData( count, name,  start, vectorBuffer);
    return vertices.size();
}

int CXX_mot::CXXSurface::updateWithVectorData(int count, const string name, int start, double *vectorBuffer){
    unsigned int iVector  = getVectorHandle(name);
    if (int(vertices.size()) < start + count){
        vertices.resize(start+count);
    }
    
    for (unsigned int i=0; i<(unsigned int) count; i++){
        vertices[i+start].setXyz(iVector, &(vectorBuffer[3*i]));
    }
    
    return vertices.size();
}


int CXX_mot::CXXSurface::updateWithPointerData(int count, const string name, int start, void **pointerBuffer){
   unsigned int iPointer = getPointerHandle(name);
   unsigned int ui_count = count;

   // There is an "uninitialized argument value" bug here from scan-build July 2016.
   
   if (int(vertices.size()) < start + count){
      vertices.resize(start+count);
   }
   for (unsigned int i=0; i<ui_count; i++){
      vertices[i+start].setPointer(iPointer, pointerBuffer[i]);
   }
   return vertices.size();
}

int CXX_mot::CXXSurface::extendTriangles(int *triangleBuffer, int count){
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

int CXX_mot::CXXSurface::numberOfTriangles() const{
    return triangles.size();
}

int CXX_mot::CXXSurface::numberOfVertices() const{
    return vertices.size();
}

int CXX_mot::CXXSurface::vertex(int iTriangle, int iCorner) const{
    return triangles[iTriangle][iCorner];
}

int CXX_mot::CXXSurface::upLoadSphere(CXX_mot::CXXSphereElement &theSphere,
				      double probeRadius, const int sense) {

   int oldVertexCount;
    CXXCoord theCentre = theSphere.centre();
    
    vector<int, CXX_old::CXXAlloc<int> > equivalence(theSphere.nVertices());
    vector<int, CXX_old::CXXAlloc<int> > uniqueAndDrawn(theSphere.nVertices());
    
    int nDrawn = 0;
    for (unsigned  i=0; i< theSphere.nVertices(); i++){
        CXXCoord comp1(theSphere.vertex(i).vertex());
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
        int verticesHandle = getVectorHandle(vertexName);
        int accessiblesHandle = getVectorHandle(accessiblesName);
        int normalsHandle = getVectorHandle(normalsName);
        int iDraw = 0;
        for (unsigned int i=0; i< theSphere.nVertices(); i++){
            if (uniqueAndDrawn[i]){
                CXXCoord vertexCoord = theSphere.vertex(i).vertex();
                if (sense == CXXSphereElement::Contact){
                    vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord);
                    CXXCoord normal = vertexCoord - theCentre;
                    CXXCoord diff(normal);
                    diff *= (theSphere.radius() - probeRadius) / theSphere.radius();
                    normal.normalise();
                    vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
                    CXXCoord vertex = theCentre + diff;
                    vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertex);
                }
                else if (sense == CXXSphereElement::Reentrant) {
                    vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
                    CXXCoord normal = theCentre - vertexCoord;
                    normal.normalise();
                    vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
                    vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, theCentre);
                }
                else if (sense == CXXSphereElement::VDW) {
                    vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
                    CXXCoord normal = vertexCoord - theCentre;
                    normal.normalise();
                    vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
                    vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord+normal);
                }
                else if (sense == CXXSphereElement::Accessible) {
                    vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
                    CXXCoord normal = vertexCoord - theCentre;
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
        void *atomBuffer[nDrawn];// = new void*[nDrawn];
        int iDraw = 0;
        for (unsigned int i=0; i< theSphere.nVertices(); i++){
            if (uniqueAndDrawn[i]){
                mmdb::PAtom anAtom;
                if ((anAtom = theSphere.vertex(i).getAtom())!=0 ){
                    atomBuffer[iDraw] = (void *) anAtom;
                }
                else atomBuffer[iDraw] = (void *) theSphere.getAtom();
                iDraw++;
            }
        }
        updateWithPointerData(nDrawn, "atom", oldVertexCount, atomBuffer);
        //delete [] atomBuffer;
    }
    // Add triangles to surface
    {
        int triangleBuffer[theSphere.nFlatTriangles()*3];// = new int[theSphere.nFlatTriangles()*3];
        int drawCount = 0;
		std::list<CXXSphereFlatTriangle, CXX_old::CXXAlloc<CXXSphereFlatTriangle> >::const_iterator trianglesEnd = 
		theSphere.getFlatTriangles().end();
		for (std::list<CXXSphereFlatTriangle, CXX_old::CXXAlloc<CXXSphereFlatTriangle> >::const_iterator triangle = 
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
                        triangleBuffer[3*drawCount+j] = index + oldVertexCount;
                    }
                }
                else {
                    for (unsigned int j=0; j<3; j++){
                        int index = equivalence[theTriangle[j]];
                        triangleBuffer[3*drawCount+j] = index + oldVertexCount;
                    }
                }
                drawCount++;
            }
        }
        extendTriangles(triangleBuffer, drawCount);
        //delete [] triangleBuffer;
    }
    return 0;
}

double CXX_mot::CXXSurface::getAtomRadius(mmdb::PAtom theAtom){
    //Here get handle of a radius data type from MMDB if such has been stored
    int iRadiusHandle = allAtomsManager->GetUDDHandle(mmdb::UDR_ATOM, "PerAtomRadius");
    double theRadius;
    if (iRadiusHandle>0){
        int success = theAtom->GetUDData (iRadiusHandle, theRadius);
        if (success != mmdb::UDDATA_Ok) theRadius = mmdb::getVdWaalsRadius(theAtom->element);
    }
    else theRadius = mmdb::getVdWaalsRadius(theAtom->element);
    return theRadius;
}

int CXX_mot::CXXSurface::selectionStringToSelHnd(mmdb::PManager allAtomsManager_in, std::string selectionString){
    int selHnd = allAtomsManager_in->NewSelection();
    char *pstring = (char *) malloc (sizeof(selectionString.c_str())+1);
    strcpy (pstring, selectionString.c_str());
    allAtomsManager_in->Select ( selHnd, mmdb::STYPE_ATOM, pstring, mmdb::SKEY_NEW);
    free (pstring);
    return selHnd;
}

mmdb::Manager *CXX_mot::CXXSurface::getMMDBManager() const{
    return allAtomsManager;
}

int CXX_mot::CXXSurface::setCoord(const string &name, int iVertex, const CXXCoord &crd){
    int iVector = getVectorHandle(name);
    iVector = vectors[name];
    if (int(vertices.size()) <= iVertex){
        vertices.resize(iVertex+1);
    }
    
    vertices[iVertex].setCoord(iVector, crd);
    
    return vertices.size();
}

int CXX_mot::CXXSurface::getIntegerUDDataOfAtom(mmdb::PAtom theAtom, int handle){
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

int CXX_mot::CXXSurface::getVectorHandle (const string name){
    if (vectors.find(name) == vectors.end()){
        int oldSize = vectors.size();
        vectors[name] = oldSize + 1;
    }
    return vectors[name];
}

int CXX_mot::CXXSurface::getReadVectorHandle (const string name){
    if (vectors.find(name) == vectors.end()){
        return -1;
    }
    return vectors[name];
}


int CXX_mot::CXXSurface::getScalarHandle (const string name){
    if (scalars.find(name) == scalars.end()){
        int oldSize = scalars.size();
        scalars[name] = oldSize + 1;
    }
    return scalars[name];
}

void CXX_mot::CXXSurface::setScalar(int scalarHandle, int iVertex, double &value){
    vertices[iVertex].setScalar(scalarHandle, value);
}

void CXX_mot::CXXSurface::setScalar(const string name, int iVertex, double &value){
    int scalarHandle=getScalarHandle(name);
    vertices[iVertex].setScalar(scalarHandle, value);
}


int CXX_mot::CXXSurface::getReadScalarHandle (const string name){
    if (scalars.find(name) == scalars.end()){
        return -1;
    }
    return scalars[name];
}

int CXX_mot::CXXSurface::getPointerHandle (const string name){
    if (pointers.find(name) == pointers.end()){
        int oldSize = pointers.size();
        pointers[name] = oldSize + 1;
    }
    return pointers[name];
}

const CXX_mot::CXXCoord &CXX_mot::CXXSurface::coordRef(int coordType, int iTriangle, int iCorner) const{
    const CXXTriangle &theTriangle(triangles[iTriangle]);
    int iVertex = theTriangle[iCorner];
    return vertices[iVertex].coordRef(coordType);
}

const CXX_mot::CXXCoord &CXX_mot::CXXSurface::coordRef(int coordType, int iVertex) const{
    return vertices[iVertex].coordRef(coordType);
}

int CXX_mot::CXXSurface::addPerVertexVector (const string name, double *coordBuffer){
    int vectorHandle = getVectorHandle(name);
    for (int i=0; i<int(vertices.size()); i++){
        vertices[i].setXyz(vectorHandle, &(coordBuffer[3*i]));
    }
    return 0;
}

int CXX_mot::CXXSurface::addPerVertexScalar (const string name, double *scalarBuffer){
    int scalarHandle = getScalarHandle(name);
    for (int i=0; i<int(vertices.size()); i++){
        vertices[i].setScalar(scalarHandle, scalarBuffer[i]);
    }
    return 0;
}

int CXX_mot::CXXSurface::addPerVertexPointer (const string name, void **pointerBuffer){
    int pointerHandle = getPointerHandle(name);
    for (int i=0; i<int(vertices.size()); i++){
        vertices[i].setPointer(pointerHandle, pointerBuffer[i]);
    }
    return 0;
}

int CXX_mot::CXXSurface::getCoord(const string &type, const int iTriangle, const int corner, double *buffer)
{
    int relevantVector = getReadVectorHandle(type);
    if (relevantVector<0) return 1;
    int iVertex = triangles[iTriangle][corner];
    return getCoord(relevantVector, iVertex, buffer);
}

int CXX_mot::CXXSurface::getCoord(const string &type, const int iVertex, double *buffer)
{
    int relevantVector = getReadVectorHandle(type);
    if (relevantVector<0) return 1;
    return getCoord(relevantVector, iVertex, buffer);
}

int CXX_mot::CXXSurface::getCoord(const int type, const int iVertex, double *buffer)
{ 
    if (type>=0 && type <= int(vectors.size())){
        const CXXCoord &theCoord = coordRef(type, iVertex);
        for (int i=0; i<4; i++){
            buffer[i] = theCoord.element(i);
        }
        return 0;
    }
    else return 1;
}

int CXX_mot::CXXSurface::getPointer(const string &type, int iVertex, void **return_p)
{
    int relevantVector;
    if (pointers.find(type) != pointers.end()){
        relevantVector = vectors[type];
        *return_p = vertices[iVertex].pointer(pointers[type]);
        return 0;
    }
    *return_p = 0;
    return 1;
}

int CXX_mot::CXXSurface::getScalar  (int handle, int iVertex, double &result){
    result = vertices[iVertex].scalar(handle);
    return 0; 
}

int CXX_mot::CXXSurface::addTriangle (const CXXTriangle &aTriangle){
    triangles.push_back(aTriangle);
    return 0;
}

void CXX_mot::CXXSurface::appendSurface(const CXXSurface &otherSurface){
    int oldNVertices = vertices.size();
    int oldNTriangles = triangles.size();
    
    if (vectors.size() == 0) vectors = otherSurface.vectorNames();
    if (scalars.size() == 0) scalars = otherSurface.scalarNames();
    if (pointers.size() == 0) pointers = otherSurface.pointerNames();
    
    vertices.insert(vertices.end(),otherSurface.getVertices().begin(), otherSurface.getVertices().end());
    triangles.insert(triangles.end(),otherSurface.getTriangles().begin(), otherSurface.getTriangles().end());
    nTriangles = triangles.size();
    
    vector<CXXTriangle, CXX_old::CXXAlloc<CXXTriangle> >::iterator triangle;
    vector<CXXTriangle, CXX_old::CXXAlloc<CXXTriangle> >::iterator triangleEnd = triangles.end();
    for (triangle = (triangles.begin() + oldNTriangles); triangle!=triangleEnd; ++triangle){
        for (int i=0; i<3; i++){
            (*triangle)[i] = (*triangle)[i] + oldNVertices;
        }
    }
}

void CXX_mot::CXXSurface::compress(double tolerance){
	vector<CXXSurfaceVertex, CXX_old::CXXAlloc<CXXSurfaceVertex> > compressedVertices;
	compressedVertices.reserve(vertices.size());
	vector<CXXTriangle, CXX_old::CXXAlloc<CXXTriangle> >compressedTriangles;
	compressedTriangles.reserve(triangles.size());
	
	int vertexHandle = getVectorHandle("vertices");
	int normalHandle = getVectorHandle("normals");
    int atomHandle = getVectorHandle("atom");
    
	vector<int, CXX_old::CXXAlloc<int> >equivalences(vertices.size());
	for (unsigned int i=0; i<vertices.size(); i++){
		bool uniqueAndDrawn = true;
		for (unsigned int j=0; j<compressedVertices.size() && uniqueAndDrawn; j++){
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
	int equivalent[3];
	for (unsigned int i=0; i<triangles.size(); i++){
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

