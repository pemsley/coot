/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
#ifdef _MVS
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include "surface.h"
#include <CXXSurface.h>
#include "mg_colour.h"
#include "atom_util.h"
#include <rgbreps.h>
#include <iostream>
#include "CXXCreator.h"
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

#ifdef _WIN32
static bool have_range_ext;
static PFNGLDRAWRANGEELEMENTSEXTPROC glDrawRangeElements;
static GLboolean CheckExtension(char *extName, const GLubyte *extString)
{
    if(!extString||!extName) return GL_FALSE;
    /*
     ** Search for extName in the extensions string.  Use of strstr()
     ** is not sufficient because extension names can be prefixes of
     ** other extension names.	Could use strtok() but the constant
     ** string returned by glGetString can be in read-only memory.
     */
    char *p = (char *)extString;
    char *end;
    int extNameLen;

    extNameLen = strlen(extName);
    end = p + strlen(p);

    while (p < end) {
	int n = strcspn(p, " ");
	if ((extNameLen == n) && (strncmp(extName, p, n) == 0)) {
	    return GL_TRUE;
	}
	p += (n + 1);
    }
    return GL_FALSE;
}

void init_range_ext(){
  const GLubyte *ext_string;
  int new_ext_supported = GL_FALSE;

  if (CheckExtension("GL_EXT_draw_range_elements", glGetString(GL_EXTENSIONS)))
    new_ext_supported = GL_TRUE;

  if(new_ext_supported){
    printf("Have extension: GL_EXT_draw_range_elements\n");
    glDrawRangeElements = (PFNGLDRAWRANGEELEMENTSEXTPROC) wglGetProcAddress("glDrawRangeElementsEXT");
    have_range_ext=true;
  }
}
#endif

void surface::set_draw_colour(GLfloat *col){}

std::vector<Primitive*> surface::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin){
	double coord[4];
	double normal[4];
	double col[4];
	std::vector<Primitive*> a;
	for (int i=0; i< theSurface->numberOfTriangles(); i++){
		std::vector<Cartesian> carts;
		std::vector<Cartesian> normals;
		std::vector<Cartesian> colours;
		int id_1 = theSurface->vertex(i,0);
		int id_2 = theSurface->vertex(i,1);
		int id_3 = theSurface->vertex(i,2);
		theSurface->getCoord("colour", id_1, col);
		theSurface->getCoord("vertices", id_1, coord);
		theSurface->getCoord("normals", id_1, normal);
		carts.push_back(Cartesian(coord));
		normals.push_back(Cartesian(normal));
		colours.push_back(Cartesian(col));
		theSurface->getCoord("colour", id_2, col);
		theSurface->getCoord("vertices", id_2, coord);
		theSurface->getCoord("normals", id_2, normal);
		carts.push_back(Cartesian(coord));
		normals.push_back(Cartesian(normal));
		colours.push_back(Cartesian(col));
		theSurface->getCoord("colour", id_3, col);
		theSurface->getCoord("vertices", id_3, coord);
		theSurface->getCoord("normals", id_3, normal);
		carts.push_back(Cartesian(coord));
		normals.push_back(Cartesian(normal));
		colours.push_back(Cartesian(col));
		TriangleElement *q = new TriangleElement(carts,col,(carts[0]+carts[1]+carts[2])/3.,alpha,textured);
		q->SetNormals(normals);
		q->SetColours(colours);
		a.push_back(q);
	}
	return a;
}

void surface::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin){
	
	double coord_1[4];
	int id_1;
	int id_2;
	int id_3;
	
	fp << "mesh2 {\n";
	fp << "  vertex_vectors {\n";
	fp << "  " << theSurface->numberOfVertices() << ",\n";
	int i;
	Cartesian p;
	for (i=0; i< theSurface->numberOfVertices()-1; i++){
		theSurface->getCoord("vertices", i, coord_1);
		p = quat.getInvMatrix()*(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])+objorigin)+Cartesian(ox,oy,oz));
		fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
	}
	theSurface->getCoord("vertices", i, coord_1);
	p = quat.getInvMatrix()*(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])+objorigin)+Cartesian(ox,oy,oz));
	fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
	fp << "  }\n";
	
	fp << "  normal_vectors {\n";
	fp << "  " << theSurface->numberOfVertices() << ",\n";
	for (i=0; i< theSurface->numberOfVertices()-1; i++){
		theSurface->getCoord("normals", i, coord_1);
		p = quat.getInvMatrix()*(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])));
		p.normalize();
		fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
	}
	theSurface->getCoord("normals", i, coord_1);
	p = quat.getInvMatrix()*(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])));
	p.normalize();
	fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
	fp << "  }\n";
	
	fp << "  texture_list {\n";
	fp << "  " << theSurface->numberOfVertices() << ",\n";
	for (i=0; i< theSurface->numberOfVertices()-1; i++){
		theSurface->getCoord("colour", i, coord_1);
		fp << "  texture{pigment{rgb< " << coord_1[0] << ", " << coord_1[1] << ", " << coord_1[2] << ">}finish {diffuse 1.0 specular 1.0}},\n";
	}
	theSurface->getCoord("colour", i, coord_1);
	fp << "  texture{pigment{rgb< " << coord_1[0] << ", " << coord_1[1] << ", " << coord_1[2] << ">}finish {diffuse 1.0 specular 1.0}}\n";
	fp << "  }\n";
	
	fp << "  face_indices {\n";
	fp << "  " << theSurface->numberOfTriangles() << ",\n";
	for (i=0; i< theSurface->numberOfTriangles()-1; i++){
		id_1 = theSurface->vertex(i,0);
		id_2 = theSurface->vertex(i,1);
		id_3 = theSurface->vertex(i,2);
		fp << "  < " << id_1 << ", " << id_2 << ", " << id_3 << ">," << id_1 << "," << id_2 << "," << id_3 << ",\n";
	}
	id_1 = theSurface->vertex(i,0);
	id_2 = theSurface->vertex(i,1);
	id_3 = theSurface->vertex(i,2);
	fp << "  < " << id_1 << ", " << id_2 << ", " << id_3 << ">," << id_1 << "," << id_2 << "," << id_3 << "\n";
	fp << "  }\n";
	
	fp << "}\n";
	
}

void surface::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){}

surface::surface () : Primitive() {
//	std::cout << " EP 1"; std::cout.flush();
	style = CCP4MG_SURFACE_SOLID;
}

surface::~surface (){
	delete theSurface;
	freeResources();
}

surface::surface (CMMDBManager *theManager, int selHnd) : Primitive() {
//	std::cout << " EP 2"; std::cout.flush();
	initArrays();
	style = CCP4MG_SURFACE_SOLID;
	theSurface = new CXXSurface();
	theSurface->calculateFromAtoms (theManager, selHnd, selHnd, 1.5, 30*DEGTORAD, false);
	//evaluateElectrostaticPotential(theManager, selHnd);
	theSurface->report();
}

surface::surface (CMMDBManager *theManager, int selHnd, int contextSelHnd) : Primitive() {
//	std::cout << " EP 3"; std::cout.flush();
	initArrays();
	style = CCP4MG_SURFACE_SOLID;
	theSurface = new CXXSurface();
	theSurface->calculateFromAtoms (theManager, selHnd, contextSelHnd, 1.5, 0.785);
}

surface::surface (const std::string &fileName) : Primitive(){
//	std::cout << " EP 4"; std::cout.flush();
	style = CCP4MG_SURFACE_SOLID;
	theSurface = new CXXSurface (fileName);
}

void surface::draw_dots(double *override_colour, int selective_override){
	glDisable(GL_LIGHTING);
	double coords[4];
	glPointSize(2.0);
	glBegin(GL_POINTS);
	for (int i=0; i< theSurface->numberOfVertices(); i++){
		if (!theSurface->getCoord("colour", i, coords)){
			coords[3] = 1.0;
			glColor4dv(coords);
		}
		if (!theSurface->getCoord("vertices", i, coords)){
			glVertex3dv(coords);
		}
	}
	glEnd();
}

void surface::generateArrays(){
	CXXSurface &mySurface = *theSurface;
	
	//If we are using VBOs, and this is not the the first time wehave been drawn,
	//Then we have to release the card memory we have nicked.  Otherwise, it is buffer
	//space we have to release
	freeResources();
	
	double coords[4];
	int nVerts = mySurface.numberOfVertices();
	
	vertices = new float[3*nVerts];
	normals = new float[3*nVerts];
	colors = new float[4*nVerts];
	for (int i=0; i< mySurface.numberOfVertices(); i++){
		//Copy vertex into vertices array	
		if (!mySurface.getCoord("vertices", i, coords)){
			for (int k=0; k<3; k++) vertices[3*i +k] = coords[k];
		}
		if (!mySurface.getCoord("normals", i, coords)){
			for (int k=0; k<3; k++) normals[3*i +k] = coords[k];
		}
		//Copy color into vertices array	
		if (!mySurface.getCoord("colour", i, coords)){
			for (int k=0; k<4; k++) colors[4*i +k] = coords[k];
		}
		else for (int k=0; k<4; k++) colors[4*i +k] = 0.5;
		colors[4*i+3] = opaqueness;
	}
	
	indices = new GLuint[3*mySurface.numberOfTriangles()];
	int idx=0;
	for (int i=0; i< mySurface.numberOfTriangles(); i++){
		for (int j=0; j<3; j++){
			indices[idx++] = (GLuint)mySurface.vertex(i,j);
		}
	}
	
	arraysGenerated = 1;
	arraysUploaded = 0;
}

void surface::bindArrays(){
	CXXSurface &mySurface = *theSurface;
	int nVerts = mySurface.numberOfVertices();
	
	glEnableClientState(GL_VERTEX_ARRAY);
	if (useVBO) {
		if (!arraysUploaded) blatBuffer(GL_ARRAY_BUFFER, cardVertexBuffer, 3*nVerts*sizeof(GLfloat), (void **) &vertices, GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBufferARB(GL_ARRAY_BUFFER,cardVertexBuffer);
#endif
	}
	glVertexPointer(3, GL_FLOAT, 0, vertices);
	
	glEnableClientState(GL_NORMAL_ARRAY);
	if (useVBO) {
		if (!arraysUploaded) blatBuffer(GL_ARRAY_BUFFER, cardNormalBuffer, 3*nVerts*sizeof(GLfloat), (void **) &normals, GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ARRAY_BUFFER,cardNormalBuffer);
#endif
	}
	glNormalPointer(GL_FLOAT, 0, normals);
	
	glEnableClientState(GL_COLOR_ARRAY);
	if (useVBO) {
		if (!arraysUploaded) blatBuffer(GL_ARRAY_BUFFER, cardColorBuffer,  4*nVerts*sizeof(GLfloat), (void **) &colors,  GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ARRAY_BUFFER,cardColorBuffer);
#endif
	}
	glColorPointer(4, GL_FLOAT, 0, colors);
	
	if (useVBO) {
		if (!arraysUploaded) blatBuffer(GL_ELEMENT_ARRAY_BUFFER, cardIndexBuffer, 3*mySurface.numberOfTriangles()*sizeof(GLuint), (void ** ) &indices,  GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cardIndexBuffer);
#endif
	}
	arraysUploaded = 1;
	return;
}	

void surface::draw(double *override_colour, int selective_override){
#ifdef _WIN32
	if(!have_range_ext)
          init_range_ext();
#endif
	if(style==CCP4MG_SURFACE_DOTS){
		draw_dots(override_colour,selective_override);
		return;
	}
	CXXSurface &mySurface = *theSurface;
	
	generateArrays();
	bindArrays();
	
	float whiteColour[] = {1.,1.,1.,1.};
	float blackColor[] = {0.,0.,0.,1.};
	
	GLenum err;
	
	if (opaqueness < 1.) {
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);		
	}
	//Save the aspects of gl state that I am going to change
	glPushAttrib(GL_COLOR_BUFFER_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT);
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);	
	
	//Let the returned colour dictate: note obligatory order of these calls
	//glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	//glEnable(GL_COLOR_MATERIAL);
	
	//Set material properties that are not per-vertex
	float specularColor[4];
	for (int i=0; i<4; i++) specularColor[i] = whiteColour[i] * specularity;
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blackColor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128); 
	glDrawRangeElements(GL_TRIANGLES, 0, mySurface.numberOfVertices(), 3*mySurface.numberOfTriangles(), GL_UNSIGNED_INT, indices);
	
	if (opaqueness < 255)
		glBlendFunc (GL_ONE, GL_ZERO);
	
	glPopAttrib();	
}

void surface::blatBuffer(int arrayType, GLuint &bufferIndex, int dataSize, void **data, int type){
#ifdef UNDERSTANDS_VBOs
	glGenBuffers(1, &bufferIndex);
	glBindBuffer(arrayType, bufferIndex);
	switch (arrayType){
		case GL_ARRAY_BUFFER:
			glBufferData(arrayType, dataSize, (GLfloat *) *data, type);
			delete [] (GLfloat *) *data;
			break;
		case GL_ELEMENT_ARRAY_BUFFER:
			glBufferData(arrayType, dataSize, (GLuint *) *data, type);
			delete [] (GLuint *) *data;
			break;
	}
	*data = (GLfloat *)BUFFER_OFFSET(0);
#endif
}

void surface::drawElement(int element){
#ifdef _WIN32
	if(!have_range_ext)
          init_range_ext();
#endif
	CXXSurface &mySurface = *theSurface;
	if (useVBO) {
		if (!arraysGenerated) generateArrays();
		//bindArrays();
		GLvoid *offset = BUFFER_OFFSET(element*3*sizeof(GLuint));
		glDrawRangeElements(GL_TRIANGLES, 0, mySurface.numberOfVertices(), 1, GL_UNSIGNED_INT , offset );
		return;
	}
	
	double coords[4];
	float whiteColour[] = {1.,1.,1.,1.};
	float blackColor[] = {0.,0.,0.,1.};
	
	if (opaqueness < 1.)
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	//Save the aspects of gl state that I am going to change
	glPushAttrib(GL_COLOR_BUFFER_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT);
	
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);
	
	//Let the returned colour dictate: note obligatory order of these calls
	//glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	//glEnable(GL_COLOR_MATERIAL);
	
	//Set material properties that are not per-vertex
	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, blackColor);
	float specularColor[4];
	for (int i=0; i<4; i++) specularColor[i] = whiteColour[i] * specularity;
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blackColor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128); 	
	
	glBegin(GL_POLYGON);
	
	for (int i=0; i<3; i++){
		//Use the colour if it has been assigned
		if (!mySurface.getCoord("colour", element, i, coords)){
			coords[3] = opaqueness / 255;
			glColor4dv(coords);
		}
		
		//Copy normal into normals array	
		if (!mySurface.getCoord("normals", element, i, coords)){
			glNormal3dv(coords);
		}	
		
		if (!mySurface.getCoord("vertices", element, i, coords)){
			glVertex3dv(coords);
		}
		
	}
	glEnd();
	glPopAttrib();
}

void surface::freeResources() {
	if (useVBO){
#ifdef UNDERSTANDS_VBOs
		if (cardVertexBuffer!=0) glDeleteBuffersARB(1, &cardVertexBuffer);
		if (cardNormalBuffer!=0) glDeleteBuffersARB(1, &cardNormalBuffer);
		if (cardColorBuffer!=0)  glDeleteBuffersARB(1, &cardColorBuffer);
		if (cardIndexBuffer!=0)  glDeleteBuffersARB(1, &cardIndexBuffer);
#endif
		cardVertexBuffer = cardNormalBuffer = cardColorBuffer = cardIndexBuffer = 0;
		vertices = normals = colors = 0;
		indices = 0;
	}
	else {
		cardVertexBuffer = cardNormalBuffer = cardColorBuffer = cardIndexBuffer = 0;
		if (colors) delete [] colors;
		if (normals) delete [] normals;
		if (vertices) delete [] vertices;
		if (indices) delete [] indices;
		vertices = normals = colors = 0;
		indices = 0;
	}
}

void add_surface(surface *surf, Displayobject &obj){
	obj.add_primitive(surf);
}


void surface::ColourSurface(CMMDBManager *molHnd, int selHnd, AtomColourVector *atomColourVector){
	if (iEval == 0) {
		iEval++;
		return;
	}
	int nSelAtoms;
	PPCAtom SelAtom;
	molHnd->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
	std::vector < double*> colours;
	for(unsigned i=0;i<unsigned(nSelAtoms);i++){
		colours.push_back(atomColourVector->GetRGB(i));
	}
	/*  
		for(unsigned i=0;i<unsigned(nSelAtoms);i++){
			colours[i][0];
			colours[i][1];
			colours[i][2];
			colours[i][3];
		}
	*/
	theSurface->report();
	theSurface->colorByColourArray(colours,molHnd,selHnd);
	
}

int surface::evaluatePhiAndColourWithDefaultScheme(CMMDBManager *theManager, const int selHnd, int contains_hydrogen ){
	CColourScheme defaultScheme;
	std::vector<float> typ;
	typ.push_back(-0.2); typ.push_back(-0.2);  typ.push_back(0.0);  typ.push_back(0.2); typ.push_back(0.2);
	std::vector<std::string> cols;
	cols.push_back("red"); cols.push_back("red");  cols.push_back("white");  cols.push_back("blue"); cols.push_back("blue");
	defaultScheme.SetSchemeFloat(typ, cols);
	return evaluatePhiAndColourWithScheme(theManager, selHnd, defaultScheme, contains_hydrogen);
}


int surface::evaluatePhiAndColourWithScheme(CMMDBManager *theManager, const int selHnd, CColourScheme &colourScheme , int contains_hydrogen) {
	//Instantiate and calculate electrostatic potential
	// contains_hydrogen = 0 => no hydrogen
	// contains_hydrogen = 0 => some H atoms in structure (currently no checks if all present)
	CXXCreator theCreator(theManager, selHnd);
	theCreator.calculate();	
	
	//Coerce map into clipper NXmap
	clipper::Cell aCell;
	clipper::NXmap<double> thePhiMap (theCreator.coerceToClipperMap(aCell));
	// writeNXMap(thePhiMap,"phi.map");
	//thePhiMap = nxmap;
	//Interpolate into this map at ssurface vertices
	if (interpolateIntoMap("vertices", "potential", thePhiMap)) return 1;
	
	//Use these values to apply a colour scheme
	if (colourByScalarValue("potential", colourScheme)) return 1;
	return 0;
}


int surface::loadMapAndColourWithScheme(std::string map_file_name, CColourScheme &colourScheme) {
	clipper::NXmap<double> phimap =  readNXMap (map_file_name);
	if (interpolateIntoMap("vertices", "potential",phimap)) return 2;
	if (colourByScalarValue("potential", colourScheme)) return 3;
	return 0;
}



int surface::interpolateIntoMap(const std::string &coordinateType, const std::string &scalarType, 
								const clipper::NXmap<double> &aMap)
{
	int coordHandle = theSurface->getReadVectorHandle(coordinateType);
	if (coordHandle<0) return 1;
	else {
		double *scaluffer = new double[theSurface->numberOfVertices()];
		double coords[4];
		for (int i=0; i<theSurface->numberOfVertices(); i++){
			//Use the colour if it has been assigned
			if (theSurface->getCoord(coordHandle, i, coords)){
			}
			Coord_orth orthogonals(coords[0], coords[1], coords[2]);
			const Coord_map mapUnits(aMap.coord_map(orthogonals));
			scaluffer[i] = aMap.interp<Interp_cubic>( mapUnits );
		}
		theSurface->addPerVertexScalar (scalarType, scaluffer);
		delete [] scaluffer;
	}
	return 0;
}

int surface::colourByScalarValue(const std::string &scalarType, CColourScheme &colourScheme){
	int scalarHandle = theSurface->getReadScalarHandle(scalarType);
	if (scalarHandle<0) return 1;
	for (int i=0; i<theSurface->numberOfVertices(); i++){
		double scalar;
		int scalarRead = theSurface->getScalar(scalarHandle, i, scalar);
		if (!scalarRead){
			std::vector <double> newColour = colourScheme.GetRGB(scalar);
			CXXCoord colour (newColour[0], newColour[1], newColour[2]);
			theSurface->setCoord("colour", i, colour);
		}
	}
	return 0;
}


clipper::NXmap<double> surface::readNXMap (std::string map_file_name){
	clipper::NXmap<double> nxmap;
	clipper::CCP4MAPfile file;
	file.open_read(map_file_name);
	file.import_nxmap( nxmap );
	file.close_read();
	return nxmap;
}


int surface::writeNXMap (const clipper::NXmap<double> &nxmap, std::string map_file_name) {
	
	clipper::Cell cell(Cell_descr(1.0,1.0,1.0,90.0,90.0,90.0));
	clipper::CCP4MAPfile file;
	file.open_write(map_file_name);
	file.set_cell(cell);
	file.export_nxmap( nxmap );
	file.close_write();
	return 0;
}

int surface::readPhiMapAndColourWithScheme (std::string map_file_name, CColourScheme &colourScheme){
	clipper::NXmap<double> thePhiMap = readNXMap(map_file_name);
	
	//thePhiMap = nxmap; 
	//Interpolate into this map at ssurface vertices
	if (interpolateIntoMap("vertices", "potential", thePhiMap)) return 1;
	
	//Use these values to apply a colour scheme
	if (colourByScalarValue("potential", colourScheme)) return 1;
	return 0;
}

int surface::writePhiMap (std::string map_file_name) {
	//return writeNXMap(thePhiMap,map_file_name);
	return 0;
}

void surface::initArrays(){
	cardVertexBuffer = 0;
	cardNormalBuffer = 0;
	cardColorBuffer = 0;
	cardIndexBuffer = 0;
	vertices = 0;
	normals = 0;
	colors = 0;
	indices = 0;
#ifdef UNDERSTANDS_VBOs
	useVBO=1;
#else
	useVBO = 0;
#endif
	useVBO=0;
	arraysGenerated = 0;
	arraysUploaded = 0;
	specularity = 1.;
	exponent = 128.;
	opaqueness = 1.;
}
