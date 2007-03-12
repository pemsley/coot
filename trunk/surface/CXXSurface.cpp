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
#include "CXXSurface.h"
#include "CXXFortranFile.h"
#include "TokenIterator.h"
#include <CXXNewHood.h>
#include <CXXTorusElement.h>
#include <CXXSphereElement.h>
#include "mmdb_tables.h"
#include "mmdb_uddata.h"
#include <CXXCircle.h>
#include <CXXTriangle.h>
#include "CXXSphereFlatTriangle.h"
#include "CXXQADSurface.h"

CXXSurface::CXXSurface()
{
	init();
}

CXXSurface::~CXXSurface(){
	triangles.resize(0);
	vertices.resize(0);
}

int CXXSurface::init()
{
	nTriangles = 0;
	triangles.resize(0);
	vertices.resize(0);
	return (0);
}

CXXSurface::CXXSurface (string path)
{
	init();
	readGraspFile (path);
}

int CXXSurface::readGraspFile(string path)
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
	else std::cout << "Failed to open grasp file [" << path << "]\n";
	
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

std::string CXXSurface::report(){
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
	for (scalar = scalars.begin(), j=0; scalar != scalars.end(); scalar++, j++){
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
	
	Delimiters delimiters(" \0\t\n~;()\"<>:{}[]+-=&*#.,/\\");	
	
	CXXFortranFile graspFile(path,"w");
	if (!graspFile.bad()) std::cout << "Opened grasp file [" << path << "]\n";
	else std::cout << "Failed to open grasp file [" << path << "]\n";
	
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
		double *xyz = vertices[i].xyzPntr(vectors["vertices"]); 
		for (int j=0; j<3; j++){
			coordBuffer[3*i + j] = xyz[j];
		}
	}
	graspFile.putFortranData((char *) coordBuffer, sizeof(float), 
							 3*vertices.size(), CXXFortranFile::FortranFloatData);
	
	//Write accessibles if present
	if (vectors.find("accessibles") != vectors.end()){
		for (i=0; i<int(vertices.size()); i++){
			double *xyz = vertices[i].xyzPntr(vectors["accessibles"]); 
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
			double *xyz = vertices[i].xyzPntr(vectors["normals"]); 
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
				triangleBuffer[3*i + j] = 1 + triangles[i].element(j);
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

CXXSurface::CXXSurface (PCMMDBManager allAtomsManager_in, const std::string selectionString,
						const std::string contextString){
	double delta = 30. * 2. * M_PI / 360.;
	double probeRadius = 1.4; 
	init();
	calculateFromAtoms( allAtomsManager_in,  selectionString,  contextString, probeRadius, delta);
}

CXXSurface::CXXSurface (PCMMDBManager allAtomsManager_in, const std::string selectionString,
						const std::string contextString, const double delta, const double probeRadius){
	init();
	calculateFromAtoms( allAtomsManager_in,  selectionString,  contextString, probeRadius, delta);
}

CXXSurface::CXXSurface (PCMMDBManager allAtomsManager_in, const std::string selectionString){
	double delta = 30. * 2. * M_PI / 360.;
	double probeRadius = 1.4;
	calculateFromAtoms( allAtomsManager_in,  selectionString, probeRadius, delta);
}

CXXSurface::CXXSurface (PCMMDBManager allAtomsManager_in, const int selHnd){
	double delta = 30. * 2. * M_PI / 360.;
	double probeRadius = 1.4;
	init();
	calculateFromAtoms( allAtomsManager_in,  selHnd,  selHnd, probeRadius, delta);
}

CXXSurface::CXXSurface (PCMMDBManager allAtomsManager_in, const int selHnd, const double delta, const double probeRadius){
	init();
	calculateFromAtoms( allAtomsManager_in,  selHnd,  selHnd, probeRadius, delta);
}

CXXSurface::CXXSurface (PCMMDBManager allAtomsManager_in, const int selHnd, const int contextSelHnd, 
						const double delta, const double probeRadius){
	init();
	calculateFromAtoms( allAtomsManager_in,  selHnd,  contextSelHnd, probeRadius, delta);
}

int CXXSurface::calculateFromAtoms(PCMMDBManager allAtomsManager_in, const std::string selectionString, const double probeRadius, 
								   const double delta){
	int selHnd = selectionStringToSelHnd(allAtomsManager_in, selectionString);
	return calculateFromAtoms(allAtomsManager_in, selHnd, selHnd, probeRadius, delta);
}

int CXXSurface::calculateFromAtoms(PCMMDBManager allAtomsManager_in, const std::string selectionString, const std::string contextString, const double probeRadius, const double delta){
	int selHnd = selectionStringToSelHnd(allAtomsManager_in, selectionString);
	int contextHnd = selectionStringToSelHnd(allAtomsManager_in, contextString);
	return calculateFromAtoms(allAtomsManager_in, selHnd, contextHnd, probeRadius, delta);
}

int CXXSurface::calculateFromAtoms(PCMMDBManager allAtomsManager_in, const int selHnd, const double probeRadius, const double delta){
	return calculateFromAtoms(allAtomsManager_in, selHnd, selHnd, probeRadius, delta);
}

int CXXSurface::calculateFromAtoms(PCMMDBManager allAtomsManager_in, const int selHnd, const int contextSelHnd, const double probeRadius, const double delta){
	
	int  nSelAtoms;
	unsigned sphereAtomNr;
	PPCAtom SelAtom;
	allAtomsManager = allAtomsManager_in;	
	
	vector<double> atomRadii;
	// now get selection index - check if this worked ...
	// this assigns number of selected atoms to nSelAtoms	
	allAtomsManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	cout << "Surface selection includes " << nSelAtoms << "atoms";
	
	//determine the maximum atomic radius
	double maxAtomRadius = 0.;
	for (int i=0; i< nSelAtoms; i++){
		atomRadii.push_back(getAtomRadius(SelAtom[i]));
		maxAtomRadius = max(maxAtomRadius, getAtomRadius(SelAtom[i]));
	}
	
	
	// now get selection index for the contextual atoms- check if this worked ...
	// this assigns number of selected atoms to nSelAtoms	
	PPCAtom ContextSelAtom;
	int nContextSelAtoms;
	allAtomsManager->GetSelIndex(contextSelHnd, ContextSelAtom, nContextSelAtoms);
	cout << "Context selection includes " << nContextSelAtoms << "atoms";
	
	//Pre identify all contacts 
	PSContact contacts;
	contacts = 0;
	int nContacts = 0;
	cout << "Off to precalculate contacts..."; cout.flush();
	allAtomsManager->SeekContacts( SelAtom, nSelAtoms, 
								   ContextSelAtom, nContextSelAtoms, 
								   0., 2.*maxAtomRadius+2.*probeRadius, 
								   0, contacts, nContacts, 0, 0);
	vector<vector<int> > neighbourhoods;
	neighbourhoods.resize(nSelAtoms);
	for (int i=0; i< nContacts; i++){
		neighbourhoods[contacts[i].id1].push_back(contacts[i].id2);
	}
	
	//Save time later, by preparing a Unit radius atom triangulated according to 
	//Delta
	const CXXCoord origin(0.,0.,0.);
	CXXSphereElement unitSphereAtOrigin(origin, 1., delta);
	
	//I don't know how to get away from the fact that we have to do an all against all comparison
	//to see if the three way nodes intersect
	vector<CXXCircleNode> circleNodes;
	//Similarly, I have to store a list of "pointy bits" in toruses that eat themselves a bit,
	//so that I can stitch such toruses to their corresponding nodes better
	vector<PointyBit> pointyBits;
	
	//loop over all atoms in selection - i.e. in pdb file
	for (int atomNr = 0;atomNr < nSelAtoms; atomNr++) { 
		
		PCAtom centralAtom = SelAtom[atomNr];
		
		if (!(atomNr%100) || atomNr==nSelAtoms-1)
			cout << "Dealing with atom number " <<atomNr << endl;
		
		double radiusOfAtom1 = getAtomRadius(centralAtom);
		CXXNewHood theNewHood(centralAtom, radiusOfAtom1, probeRadius);
		
		//We have precalculated neighbours of the central atom, and now can use that 
		//to our advantage
		for (sphereAtomNr = 0; sphereAtomNr < neighbourhoods[atomNr].size(); sphereAtomNr++) {
			PCAtom sphereAtom = ContextSelAtom[neighbourhoods[atomNr][sphereAtomNr]];
			double radiusOfAtom2 = getAtomRadius(sphereAtom);
			theNewHood.addAtom(sphereAtom, radiusOfAtom2);
		}
		
		//Find the non-hidden segments of the circles
		theNewHood.findSegments();
		
		//Draw The patch edges as torus elements;
		int nCircles = theNewHood.nCircles(); 
		
		//For now, copy the unit sphere to atom position, and translate and scale to become a vdw sphere
		CXXSphereElement vdwSphere(unitSphereAtOrigin);
		vdwSphere.scaleBy(radiusOfAtom1+probeRadius);
		
		vdwSphere.translateBy (CXXCoord( centralAtom->x, centralAtom->y, centralAtom->z));
		vdwSphere.setAtom(centralAtom);
		
		//Use the patch edges to trim the triangles of the VDW sphere
		
		for (int iCircle = 0; iCircle< nCircles && vdwSphere.getNDrawnTriangles(); iCircle ++){
			const CXXCircle &theCircle = theNewHood.getCircle(iCircle);
			if (!theCircle.getEaten()) vdwSphere.trimBy(theCircle);
		}
		
		//We collect the sphere vertices that have been cut by circles so that we can use 
		//them in some further stitching later
		vector<vector<CXXCircleNode > >rawEdges(theNewHood.nCircles());
		if (vdwSphere.nFlatTriangles()) {
			for (unsigned iVertex = 0; iVertex < vdwSphere.nVertices(); iVertex++){
				if (vdwSphere.vertex(iVertex).doDraw() &&
					vdwSphere.vertex(iVertex).getIntersector() != 0){
					
					int circleKnown = 0;
					for (int iCircle = 0; iCircle< theNewHood.nCircles() && !circleKnown; iCircle ++){
						const CXXCircle &theCircle(theNewHood.getCircle(iCircle));						
						if (!theCircle.getEaten()){
							if (vdwSphere.vertex(iVertex).getIntersector() == &theCircle){
								circleKnown = 1;
								//Generate an extra node whose "flag" stores the 
								CXXCircleNode extraNode(&theCircle,  0, vdwSphere.vertex(iVertex).vertex(),iVertex);
								extraNode.setReference(theCircle.getReferenceUnitRadius());
								
								rawEdges[iCircle].push_back(extraNode);
							}
						}
					}
				}
			}
			
			//Create torus elements for each of the segments, and elaborate them with the additional vertices
			//that correspond to ragged edges from the sphere
			for (int iCircle = 0; iCircle< theNewHood.nCircles(); iCircle ++){
				const CXXCircle &theCircle(theNewHood.getCircle(iCircle));	
				
				if (!theCircle.getEaten()){
					
					//Here cause the sphere to identify triangles that lie on each circle, and assign to them
					//"absolute" omega theta values:  this will be used later when addo=ing torus vertices
					vdwSphere.flagCutTriangles(theCircle);
					
					//Loop over segments of the circle, generating torus elements, and adding sphere nodes where
					//apropriate
					
					int nToruses;
					for (unsigned iTorus = 0; iTorus < theCircle.nSegments(); iTorus++){
						CXXTorusElement theTorus (theCircle, iTorus, delta, probeRadius);
						if (theTorus.getPointyBit()){
							PointyBit aPointyBit;
							aPointyBit.coord = *(theTorus.getPointyBit());
							aPointyBit.atomI = centralAtom;
							aPointyBit.atomJ = theCircle.getAtomJ();
							aPointyBit.isNull = 0;
							pointyBits.push_back(aPointyBit);
//							CXXCoord &newNode(pointyBits.back().coord);
//							cout << "select object banana add ball " << newNode[0] << " " << newNode[1] << " " << newNode[2] <<" 0.05\n";
											  
						}
						int startFound = 0;
						for (unsigned iRawEdge=0; iRawEdge < rawEdges[iCircle].size() && !startFound; iRawEdge++){
							theTorus.addEdgeVertex(rawEdges[iCircle][iRawEdge]);
						}
						theTorus.upload(this);  
						
						//Here cause the sphere to subdivide the triangles that are interrupted by
						//nodes around this torus
						
						vdwSphere.addTorusVertices(theTorus);
					}	
				}
			}
			
			//Copy surface triangles from this sphere into the surface
			upLoadSphere(vdwSphere, probeRadius, CXXSphereElement::Outside);
		}
		
		//Now collect a list of points where the probes are in contact with three atoms
		for (int iCircle = 0; iCircle< theNewHood.nCircles(); iCircle ++){
			const CXXCircle &theCircle = theNewHood.getCircle(iCircle);
			if (!theCircle.getEaten()){
				for (unsigned iNode = 0; iNode<theCircle.getNNodes(); iNode++){
					const CXXCircleNode &aNode(theCircle.getNode(iNode));
					if (!aNode.isDeleted()){
						//Rolling between two atoms in the selection, towards a third atom 
						//in the selection
						//This can happen in six possible ways, of which we accept only 1
						if (aNode.getAtomJ()->isInSelection(selHnd) &&
							aNode.getAtomK()->isInSelection(selHnd)) {
							if ( centralAtom->serNum < aNode.getAtomJ()->serNum &&
								 aNode.getAtomJ()->serNum < aNode.getAtomK()->serNum ){
								circleNodes.push_back(aNode);	
							}
						} 
						//Rolling between two atoms in the selection, towards a third atom 
						//*not* in the selection
						//This can happen in two possible ways, of which we accept only 1
						else if (aNode.getAtomJ()->isInSelection(selHnd)) {
							if ( centralAtom->serNum < aNode.getAtomJ()->serNum ){
								circleNodes.push_back(aNode);	
							}
						} 
						//Rolling between one atom in selection and one atom not in, towards a third atom
						//that is in the selection.  Every case where
						//this can happen has an equivalent case above from which we keep the node
						else if (aNode.getAtomK()->isInSelection(selHnd)) {
						}
						//Rolling between one atom in selection and one atom not in, towards a third atom
						//that is *not* in the selection.  						
						//This can happen in two possible ways, of which we accept only 1
						else {
							if ( aNode.getAtomJ()->serNum < aNode.getAtomK()->serNum ){
								circleNodes.push_back(aNode);	
							}
						}
					}
				}
			}
		}
	}
	cout << "Dealing with " << circleNodes.size() << " Nodes and " << pointyBits.size() <<" pointyBits"<< endl;
	//For now I think I have to do this in two passes:  first one to identify all raw edges
	//second to cause the sphere patches to be drawn appropriately (i.e. incorporating raw edges)
	//My hope is that the wonderfulness of all this will allow lower density triangulation
	
	vector< vector <NodeCoordPair> > rawEdges(circleNodes.size());
	for (unsigned i=0; i<circleNodes.size(); i++) rawEdges[i].reserve(20);
	//Draw the 3-way nodes around the probe trajectories
	for (unsigned iNode = 0; iNode < circleNodes.size(); iNode++){
		const CXXCircleNode &aNode(circleNodes[iNode]);
		//Here we create hoods just as we did for atoms
		CXXNewHood nodeHood(aNode, probeRadius);
		//Use the other nodes to identify where they should cause the central node not to be drawn	
		for (unsigned iNode2=0; iNode2 < circleNodes.size(); iNode2++){
			if (iNode != iNode2){
				const CXXCircleNode &aNode2 = circleNodes[iNode2];
				if (aNode2.getCoord() !=  aNode.getCoord()){
					nodeHood.addNodeAsAtom(aNode2, iNode2);
				}
			}
		}
		CXXSphereElement nodeSphere(aNode, delta, probeRadius, selHnd, pointyBits);
		//Now use the circles identified in the above hood to trim the node sphere element
		for (int iCircle = 0; iCircle< nodeHood.nCircles() && nodeSphere.getNDrawnTriangles (); iCircle ++){
			const CXXCircle &theCircle(nodeHood.getCircle(iCircle));
			nodeSphere.trimBy(theCircle);
		}
		for (int iVertex = 0; iVertex < nodeSphere.nVertices(); iVertex++){
			const CXXSphereNode &vertex(nodeSphere.vertex(iVertex));
			if (vertex.doDraw()){
				const CXXCircle *intersector(vertex.getIntersector());
				if (intersector){
					int circleAssigned = 0;
					for (int iCircle = 0; iCircle< nodeHood.nCircles() && !circleAssigned; iCircle ++){
						const CXXCircle &theCircle(nodeHood.getCircle(iCircle));
						if (&theCircle == intersector){
							circleAssigned = 1;
							NodeCoordPair ncp;
							ncp.iNode= iNode;
							ncp.coord = vertex.vertex();
							rawEdges[theCircle.getNodeNumber()].push_back(ncp);
						}
					}
				}
			}
		}
	}
	//SecondPass to add in the edge vertices to each sphere and draw the composite
	for (unsigned iNode = 0; iNode < circleNodes.size(); iNode++){
		const CXXCircleNode &aNode(circleNodes[iNode]);
		//Here we create hoods just as we did for atoms
		CXXNewHood nodeHood(aNode, probeRadius);
		//Use the other nodes to identify where they should cause the central node not to be drawn	
		for (unsigned iNode2=0; iNode2 < circleNodes.size(); iNode2++){
			if (iNode != iNode2){
				const CXXCircleNode &aNode2 = circleNodes[iNode2];
				nodeHood.addNodeAsAtom(aNode2, iNode2);
			}
		}
		CXXSphereElement nodeSphere(aNode, delta, probeRadius, selHnd, pointyBits);
		//Now use the circles identified in the above hood to trim the node sphere element
		for (int iCircle = 0; iCircle< nodeHood.nCircles() && nodeSphere.getNDrawnTriangles(); iCircle ++){
			const CXXCircle &theCircle(nodeHood.getCircle(iCircle));			
			nodeSphere.trimBy(theCircle);
		}
		
		for (int iCircle = 0; iCircle< nodeHood.nCircles(); iCircle ++){
			const CXXCircle &theCircle(nodeHood.getCircle(iCircle));
			nodeSphere.flagCutTriangles(theCircle);
			for (unsigned iExtraNode = 0; iExtraNode < rawEdges[iNode].size(); iExtraNode++){
				NodeCoordPair &ncp(rawEdges[iNode][iExtraNode]);
				if (ncp.iNode == theCircle.getNodeNumber()){
					CXXCircleNode extraNode(&theCircle, 0, ncp.coord, 0);
					extraNode.setReference(theCircle.getReferenceUnitRadius());
					nodeSphere.addVertex(extraNode);
				}
			}
		}
		upLoadSphere(nodeSphere, probeRadius, CXXSphereElement::Inside);
	}
	colorByAssignedAtom();
	report();
	
	return 0;
}

int CXXSurface::assignAtom (PCMMDBManager allAtomsManager_in, int selHnd){
	void **pointerBuffer;
	PPCAtom selAtom;
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

int CXXSurface::colorByColourArray(const std::vector<double*> &colours, CMMDBManager *molHnd, int selHnd){
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
	PPCAtom SelAtom;
	molHnd->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
	int udd = molHnd->GetUDDHandle ( UDR_ATOM,"tmp_atom_int" );
	if (udd <= 0 ) {
		udd = -1;
		udd = molHnd->RegisterUDInteger ( UDR_ATOM,"tmp_atom_int" );
		if (udd <= 0 ) return udd;
	}
	for(int i=0; i<nSelAtoms; i++)
		SelAtom[i]->PutUDData(udd,i);
	
	for (int i=0; i< int(vertices.size()); i++){
		PCAtom theAtom = (PCAtom) vertices[i].pointer(pointers["atom"]);
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
	PCAtom theAtom;
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
		theAtom = (PCAtom) vertices[i].pointer(pointers["atom"]);
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

int CXXSurface::extendWithVectorData(int count, const string name, double *vectorBuffer){
	int start = vertices.size();
	updateWithVectorData( count, name,  start, vectorBuffer);
	return vertices.size();
}

int CXXSurface::updateWithVectorData(int count, const string name, int start, double *vectorBuffer){
	unsigned int iVector  = getVectorHandle(name);
	if (int(vertices.size()) < start + count){
		vertices.resize(start+count);
	}
	
	for (unsigned int i=0; i<(unsigned int) count; i++){
		vertices[i+start].setXyz(iVector, &(vectorBuffer[3*i]));
	}
	
	return vertices.size();
}


int CXXSurface::updateWithPointerData(int count, const string name, int start, void **pointerBuffer){
	unsigned int iPointer = getPointerHandle(name);
	if (int(vertices.size()) < start + count){
		vertices.resize(start+count);
	}
	
	for (unsigned int i=0; i<(unsigned int) count; i++){
		vertices[i+start].setPointer(iPointer, pointerBuffer[i]);
	}
	
	return vertices.size();
}

int CXXSurface::extendTriangles(int *triangleBuffer, int count){
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

int CXXSurface::numberOfTriangles() const{
	return triangles.size();
}

int CXXSurface::numberOfVertices() const{
	return vertices.size();
}

int CXXSurface::vertex(int iTriangle, int iCorner) const{
	return triangles[iTriangle].element(iCorner);
}

int CXXSurface::calculateQADFromAtoms(PCMMDBManager theManager, const std::string selString, const double probeRadius, const double sample){
	int selHnd = selectionStringToSelHnd(theManager, selString);
	calculateQADFromAtoms(theManager, selHnd, selHnd, probeRadius, sample);
	return 0;
}

int CXXSurface::calculateQADFromAtoms(PCMMDBManager theManager, const std::string selString, const std::string contextString, const double probeRadius, const double sample){
	int selHnd = selectionStringToSelHnd(theManager, selString);
	int contextHnd = selectionStringToSelHnd(theManager, contextString);
	calculateQADFromAtoms(theManager, selHnd, contextHnd, probeRadius, sample);
	return 0;
}

int CXXSurface::calculateQADFromAtoms(PCMMDBManager theManager, const int selHnd, const double probeRadius, const double sample){
	calculateQADFromAtoms(theManager, selHnd, selHnd, probeRadius, sample);
	return 0;
}

int CXXSurface::calculateQADFromAtoms(PCMMDBManager theManager, const int selHnd, const int contextHnd, const double probeRadius, const double sample){
	CXXQADSurface theQADSurface(theManager, selHnd, probeRadius, sample);
	
	//	Add vertices to surface
	int oldVertexCount;
	int nVertices;
	{
		vector<clipper::Coord_orth> &qadVertices = theQADSurface.getVertices();
		vector<clipper::Coord_orth> &qadNormals = theQADSurface.getNormals();
		double *verticesBuffer = new double[qadVertices.size()*3];
		double *normalsBuffer = new double[qadVertices.size()*3];
		for (unsigned i=0; i< qadVertices.size(); i++){
			for (int j=0; j<3; j++) {
				verticesBuffer[3*i+j] = qadVertices[i][j];
				normalsBuffer[3*i+j] = qadNormals[i][j];
			}
		}
		oldVertexCount = numberOfVertices();
		updateWithVectorData(qadVertices.size(), "vertices", oldVertexCount, verticesBuffer);
		updateWithVectorData(qadVertices.size(), "normals", oldVertexCount, normalsBuffer);
		nVertices = qadVertices.size();
		delete [] verticesBuffer;
		delete [] normalsBuffer;
	}	
	//Determine which atom belongs to each vertex
	
	assignAtom(theManager, selHnd);
	{
		vector<int> &qadTriangles = theQADSurface.getTriangles();
		for (unsigned int i=0; i<qadTriangles.size()/3; i++){
			addTriangle(CXXTriangle(qadTriangles[3*i]+oldVertexCount, 
									qadTriangles[3*i+1]+oldVertexCount, 
									qadTriangles[3*i+2]+oldVertexCount));
			nTriangles++;
		}
	}
	colorByAssignedAtom();
	report();
	return 0;
}

int CXXSurface::upLoadSphere(CXXSphereElement &theSphere, double probeRadius, const int sense){
	int oldVertexCount;
	CXXCoord theCentre = theSphere.centre();
	
	vector<int> equivalence(theSphere.nVertices());
	vector<int> uniqueAndDrawn(theSphere.nVertices());
	
	int nDrawn = 0;
	for (unsigned  i=0; i< theSphere.nVertices(); i++){
		CXXCoord comp1(theSphere.vertex(i).vertex());
		uniqueAndDrawn[i] = 0;
		if (theSphere.vertex(i).doDraw()){
			uniqueAndDrawn[i] = 1;
			//The following bit would reduce redundancy of vertex array
			//at expense of slower performance at create time
			/*
			 for (int j=0; j< i && uniqueAndDrawn[i]; j++){
				 CXXCoord comp2(theSphere.vertex(j).vertex());
				 if (uniqueAndDrawn[j]){
					 if (comp1 == comp2){
						 uniqueAndDrawn[i] = 0;
						 equivalence[i] = equivalence[j];
					 }
				 }
			 }
			 */
			if (uniqueAndDrawn[i]){
				equivalence[i] = nDrawn++;
			}
		}
	}
	
	{
		oldVertexCount = numberOfVertices();
		vertices.resize(oldVertexCount+nDrawn);
		int verticesHandle = getVectorHandle("vertices");
		int accessiblesHandle = getVectorHandle("accessibles");
		int normalsHandle = getVectorHandle("normals");
		int iDraw = 0;
		for (unsigned int i=0; i< theSphere.nVertices(); i++){
			if (uniqueAndDrawn[i]){
				CXXCoord vertexCoord = theSphere.vertex(i).vertex();
				if (sense == CXXSphereElement::Outside){
					vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord);
					CXXCoord normal = vertexCoord - theCentre;
					CXXCoord diff(normal);
					diff *= (theSphere.radius() - probeRadius) / theSphere.radius();
					normal.normalise();
					vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
					CXXCoord vertex = theCentre + diff;
					vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertex);
				}
				else {
					vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
					CXXCoord normal = theCentre - vertexCoord;
					normal.normalise();
					vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
					vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, theCentre);
				}
				iDraw++;
			}
		}
	}
	//Add atom pointers to the surface
	{
		void **atomBuffer = new void*[nDrawn];
		int iDraw = 0;
		for (unsigned int i=0; i< theSphere.nVertices(); i++){
			if (uniqueAndDrawn[i]){
				PCAtom anAtom;
				if ((anAtom = theSphere.vertex(i).getAtom())!=0 ){
					atomBuffer[iDraw] = (void *) anAtom;
				}
				else atomBuffer[iDraw] = (void *) theSphere.getAtom();
				iDraw++;
			}
		}
		updateWithPointerData(nDrawn, "atom", oldVertexCount, atomBuffer);
		delete [] atomBuffer;
	}
	// Add triangles to surface
	{
		int *triangleBuffer = new int[theSphere.nFlatTriangles()*3];
		int drawCount = 0;
		for (unsigned i=0; i<theSphere.nFlatTriangles(); i++){
			if (theSphere.flatTriangle(i).doDraw()){
				if (sense == CXXSphereElement::Outside){
					for (unsigned int j=0; j<3; j++){
						int index = equivalence[theSphere.flatTriangle(i)[2-j]];
						triangleBuffer[3*drawCount+j] = index + oldVertexCount;
					}
				}
				else {
					for (unsigned int j=0; j<3; j++){
						int index = equivalence[theSphere.flatTriangle(i)[j]];
						triangleBuffer[3*drawCount+j] = index + oldVertexCount;
					}
				}
				drawCount++;
			}
		}
		extendTriangles(triangleBuffer, drawCount);
		delete [] triangleBuffer;
	}
	return 0;
}

double CXXSurface::getAtomRadius(PCAtom theAtom){
	//Here get handle of a radius data type from MMDB if such has been stored
	int iRadiusHandle = allAtomsManager->GetUDDHandle(UDR_ATOM, "PerAtomRadius");
	double theRadius;
	if (iRadiusHandle>0){
		int success = theAtom->GetUDData (iRadiusHandle, theRadius);
		if (success != UDDATA_Ok) theRadius = getVdWaalsRadius(theAtom->element);
	}
	else theRadius = getVdWaalsRadius(theAtom->element);
	return theRadius;
}

int CXXSurface::selectionStringToSelHnd(PCMMDBManager allAtomsManager_in, std::string selectionString){
	int selHnd = allAtomsManager_in->NewSelection();
	char *pstring = (char *) malloc (sizeof(selectionString.c_str())+1);
	strcpy (pstring, selectionString.c_str());
	allAtomsManager_in->Select ( selHnd, STYPE_ATOM, pstring, SKEY_NEW);
	free (pstring);
	return selHnd;
}

CMMDBManager *CXXSurface::getMMDBManager() const{
	return allAtomsManager;
}

int CXXSurface::setCoord(const string &name, int iVertex, const CXXCoord &crd){
	int iVector = getVectorHandle(name);
	iVector = vectors[name];
	if (int(vertices.size()) <= iVertex){
		vertices.resize(iVertex+1);
	}
	
	vertices[iVertex].setCoord(iVector, crd);
	
	return vertices.size();
}

int CXXSurface::getIntegerUDDataOfAtom(PCAtom theAtom, int handle){
	int result;
	int rc = theAtom->GetUDData(handle, result);
	switch (rc)  {
		
		case  UDDATA_WrongUDRType :
			printf ( " wrong UDD registration type\n" );
			break;
			
		case  UDDATA_WrongHandle  :
			printf ( " wrong UDD handle\n" );
			break;
			
		case  UDDATA_NoData :
			printf ( " UDD not found.\n" );
			break;
			
		case  UDDATA_Ok :
			break;			
	}	
	
	return result;
}

int CXXSurface::getVectorHandle (const string name){
	if (vectors.find(name) == vectors.end()){
		int oldSize = vectors.size();
		vectors[name] = oldSize + 1;
	}
	return vectors[name];
}

int CXXSurface::getReadVectorHandle (const string name){
	if (vectors.find(name) == vectors.end()){
		return -1;
	}
	return vectors[name];
}


int CXXSurface::getScalarHandle (const string name){
	if (scalars.find(name) == scalars.end()){
		int oldSize = scalars.size();
		scalars[name] = oldSize + 1;
	}
	return scalars[name];
}

void CXXSurface::setScalar(int scalarHandle, int iVertex, double &value){
	vertices[iVertex].setScalar(scalarHandle, value);
}

void CXXSurface::setScalar(const string name, int iVertex, double &value){
	int scalarHandle=getScalarHandle(name);
	vertices[iVertex].setScalar(scalarHandle, value);
}


int CXXSurface::getReadScalarHandle (const string name){
	if (scalars.find(name) == scalars.end()){
		return -1;
	}
	return scalars[name];
}

int CXXSurface::getPointerHandle (const string name){
	if (pointers.find(name) == pointers.end()){
		int oldSize = pointers.size();
		pointers[name] = oldSize + 1;
	}
	return pointers[name];
}

const CXXCoord &CXXSurface::coordRef(int coordType, int iTriangle, int iCorner) const{
	const CXXTriangle &theTriangle(triangles[iTriangle]);
	int iVertex = theTriangle[iCorner];
	return vertices[iVertex].coordRef(coordType);
}

const CXXCoord &CXXSurface::coordRef(int coordType, int iVertex) const{
	return vertices[iVertex].coordRef(coordType);
}

int CXXSurface::addPerVertexVector (const string name, double *coordBuffer){
	int vectorHandle = getVectorHandle(name);
	for (int i=0; i<int(vertices.size()); i++){
		vertices[i].setXyz(vectorHandle, &(coordBuffer[3*i]));
	}
	return 0;
}

int CXXSurface::addPerVertexScalar (const string name, double *scalarBuffer){
	int scalarHandle = getScalarHandle(name);
	for (int i=0; i<int(vertices.size()); i++){
		vertices[i].setScalar(scalarHandle, scalarBuffer[i]);
	}
	return 0;
}

int CXXSurface::addPerVertexPointer (const string name, void **pointerBuffer){
	int pointerHandle = getPointerHandle(name);
	for (int i=0; i<int(vertices.size()); i++){
		vertices[i].setPointer(pointerHandle, pointerBuffer[i]);
	}
	return 0;
}

int CXXSurface::getCoord(const string &type, const int iTriangle, const int corner, double *buffer)
{
	int relevantVector = getReadVectorHandle(type);
	if (relevantVector<0) return 1;
	int iVertex = triangles[iTriangle][corner];
	return getCoord(relevantVector, iVertex, buffer);
}

int CXXSurface::getCoord(const string &type, const int iVertex, double *buffer)
{
	int relevantVector = getReadVectorHandle(type);
	if (relevantVector<0) return 1;
	return getCoord(relevantVector, iVertex, buffer);
}

int CXXSurface::getCoord(const int type, const int iVertex, double *buffer)
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

int CXXSurface::getPointer(const string &type, int iVertex, void **return_p)
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

int CXXSurface::getScalar  (int handle, int iVertex, double &result){
	result = vertices[iVertex].scalar(handle);
	return 0; 
}


int CXXSurface::addTriangle (const CXXTriangle &aTriangle){
	triangles.push_back(aTriangle);
	return 0;
}
