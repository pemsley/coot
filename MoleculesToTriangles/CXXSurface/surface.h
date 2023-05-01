#ifndef _CCP4MG_Surface_
#define _CCP4MG_Surface_


#if defined (linux)
#undef GLX_GLXEXT_LEGACY
#define GL_GLEXT_PROTOTYPES
#endif

#define GL_GLEXT_PROTOTYPES 1
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#endif

#include <string>
#include "mmdb2/mmdb_manager.h"
#include "cdisplayobject.h"
#include "cprimitive.h"
#include "mg_colour.h"
#include "atom_util.h"
#include <CXXSurface.h>
//#include "clipper/core/nxmap.h"
#include <clipper/clipper.h>

//see if we understand about VertexBufferObjects (i.e. on-graphics-card buffers for object data
//#undef GL_ARRAY_BUFFER
#ifdef GL_ARRAY_BUFFER
#define UNDERSTANDS_VBOs
#else
#define GL_ARRAY_BUFFER 0x8892
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_STATIC_DRAW 0x88E4
#undef UNDERSTANDS_VBOs
#endif

#define BUFFER_OFFSET(offset) ((char *)NULL + (offset))

enum { CCP4MG_SURFACE_SOLID, CCP4MG_SURFACE_DOTS, CCP4MG_SURFACE_MESH };
class CXXSurface;

class surface : public Primitive {
private:
	CXXSurface *theSurface;
	//clipper::NXmap<double> thePhiMap;
	//clipper::NXmap<double> theUsersMap;
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
        double dotSpacing;
        double dotSize;
        double lineWidth;
public:
		surface(); 
		~surface();
		surface(mmdb::Manager* *theManager, int selHnd, double probe_radius, double delta, bool blend_edges); 
		surface(mmdb::Manager* *theManager, int selHnd, int contextSelHnd, double probe_radius, double delta, bool blend_edges); 
		surface(const std::string &fileName);
                int calculate (mmdb::Manager* *theManager, int selHnd, double probe_radius, double delta, bool blend_edges);
                int calculate_with_context (mmdb::Manager* *theManager, int selHnd, int contextSelHnd, double probe_radius, double delta, bool blend_edges); 

		void draw(const double *override_colour=0, int selective_override=0);
		void drawElement(int element);
		void draw_dots(const double *override_colour=0, int selective_override=0);
		void SetDotSpacing(double dotSpacing_in) {dotSpacing = dotSpacing_in;};
		double GetDotSpacing() const { return dotSpacing;};
		void SetDotSize(double dotSize_in) {dotSize = dotSize_in;};
		double GetDotSize() const { return dotSize;};
		void SetLineWidth(double lineWidth_in) {lineWidth = lineWidth_in;};
		double GetLineWidth() const { return lineWidth;};
                //void surface::SubdivideDotsLine(double *v1, double *v2, bool draw_first=true);
		void set_style(int style_in) {style = style_in;};
		void set_draw_colour(const GLfloat *col=0);
		void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
		void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
		void ColourSurface(mmdb::Manager* *theManager, int selHnd, const AtomColourVector &atm_col_vec);
		std::string report() { return theSurface->report(); }
		
		int evaluatePhiAndColourWithScheme(mmdb::Manager* *theManager, const int selHnd, const int context_selHnd, CColourScheme &colourScheme, int contains_hydrogen);
		int evaluatePhiAndColourWithDefaultScheme(mmdb::Manager* *theManager, const int selHnd,  const int context_selHnd, int contains_hydrogen );
		int interpolateIntoMap(const std::string &coordinateType, 
							   const std::string &scalarType, 
							   const clipper::NXmap<double> &aMap);
		int colourByScalarValue(const std::string &scalarType, CColourScheme &colourScheme);
		clipper::NXmap<double> readNXMap (std::string map_file_name);
		int writeNXMap (const clipper::NXmap<double> &map, std::string map_file_name);  
		int readPhiMapAndColourWithScheme (std::string map_file_name, CColourScheme &colourScheme);
		int loadMapAndColourWithScheme(std::string map_file_name, CColourScheme &colourScheme);
		int writePhiMap(std::string map_file_name);
                int writeGraspFile (std::string map_file_name);
                int readGraspFile (std::string map_file_name);
		bool isLine() const;
		std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start=-1, int end=-1) const ;
		int GetStyle() const { return style ;};
		std::vector<Primitive*> GetRawDots() const;
                int GetNumberOfSimplePrimitives() const;
                void forceRegenerateArrays();
		void SetAlpha(double alpha_in){ alpha = alpha_in; arraysGenerated=0; arraysUploaded=0; } ;
		
	};
	
	
	void add_surface(surface *surf, Displayobject &obj);
#endif
	
