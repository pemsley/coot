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
#ifndef _CCP4MG_Surface_
#define _CCP4MG_Surface_

#define GL_GLEXT_PROTOTYPES 1
#include <GL/gl.h>

#include <string>
#include <mmdb/mmdb_manager.h>
#include "cdisplayobject.h"
#include "cprimitive.h"
#include "mg_colour.h"
#include "atom_util.h"
#include <CXXSurface.h>
//#include "clipper/core/nxmap.h"
#include <clipper/clipper.h>

//see if we understand about VertexBufferObjects (i.e. on-graphics-card buffers for object data
#ifdef GL_ARRAY_BUFFER
#define UNDERSTANDS_VBOs
#else
#define GL_ARRAY_BUFFER 0x8892
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_STATIC_DRAW 0x88E4
#undef UNDERSTANDS_VBOs
#endif

#define BUFFER_OFFSET(offset) ((char *)NULL + (offset))

enum { CCP4MG_SURFACE_SOLID, CCP4MG_SURFACE_DOTS };
class CXXSurface;

class surface : public Primitive {
private:
	CXXSurface *theSurface;
	//clipper::NXmap<double> thePhiMap;
	//clipper::NXmap<double> theUsersMap;
	int iEval;
	int style;
	
	double specularity;
	double opaqueness;
	double exponent;
	
	GLuint cardVertexBuffer;
	GLuint cardNormalBuffer;
	GLuint cardColorBuffer;
	GLuint cardIndexBuffer;
	
	GLfloat *vertices;
	GLfloat *normals;
	GLfloat *colors; 
	GLuint  *indices;
	
	int useVBO;
	int arraysGenerated;
	int arraysUploaded;

	void freeResources();
	void generateArrays();
	void uploadArrays();
	void bindArrays();
	void blatBuffer(int arrayType, GLuint &bufferIndex, int dataSize, void ** data, int type);
		
	void initArrays();
public:
		surface(); 
		~surface();
		surface(CMMDBManager *theManager, int selHnd); 
		surface(CMMDBManager *theManager, int selHnd, int contextSelHnd); 
		surface(const std::string &fileName);
		void draw(double *override_colour=0, int selective_override=0);
		void drawElement(int element);
		void draw_dots(double *override_colour=0, int selective_override=0);
		void set_style(int style_in) {style = style_in;};
		void set_draw_colour(GLfloat *col=0);
		void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin);
		void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
		void ColourSurface(CMMDBManager *theManager, int selHnd, AtomColourVector *atm_col_vec);
		std::string report() { return theSurface->report(); }
		
		int evaluatePhiAndColourWithScheme(CMMDBManager *theManager, const int selHnd, CColourScheme &colourScheme, int contains_hydrogen);
		int evaluatePhiAndColourWithDefaultScheme(CMMDBManager *theManager, const int selHnd,int contains_hydrogen );
		int interpolateIntoMap(const std::string &coordinateType, 
							   const std::string &scalarType, 
							   const clipper::NXmap<double> &aMap);
		int colourByScalarValue(const std::string &scalarType, CColourScheme &colourScheme);
		clipper::NXmap<double> readNXMap (std::string map_file_name);
		int writeNXMap (const clipper::NXmap<double> &map, std::string map_file_name);  
		int readPhiMapAndColourWithScheme (std::string map_file_name, CColourScheme &colourScheme);
		int loadMapAndColourWithScheme(std::string map_file_name, CColourScheme &colourScheme);
		int writePhiMap(std::string map_file_name);
		bool isLine() const {return false;};
		std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
		
	};
	
	
	void add_surface(surface *surf, Displayobject &obj);
#endif
	
