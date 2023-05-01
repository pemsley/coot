#ifdef _WIN32
#include <windows.h>
#endif

#if defined (linux)
#undef GLX_GLXEXT_LEGACY
#define GL_GLEXT_PROTOTYPES
#endif

#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#ifndef sgi
#include <GL/glext.h>
#endif
#endif
#include <iostream>

#include "surface.h"
#include <CXXSurface.h>
#include "mg_colour.h"
#include "atom_util.h"
#include "subdivide.h"
#include <rgbreps.h>
#include <iostream>
#include "CXXCreator.h"
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

class Thing {
    public:
      int id1, id2,id3;
      double proj;
      Thing(double p, int i1,int i2,int i3){ id1=i1; id2=i2; id3=i3; proj = p;};
};

class ZSortThing{
  public:
    int operator()(const Thing &p1, const Thing &p2) const
       { return p1.proj < p2.proj; }
};

void SortIndices(GLfloat *new_vertices, unsigned int *new_indices, int new_ntriangles){
  GLint viewport[4];
  GLdouble mvmatrix[16];
  GLdouble projmatrix[16];
  GLdouble wx, wy, wz;
  GLint realy; /* OpenGL y coordinate position. Why do we do this? */

  glGetIntegerv(GL_VIEWPORT,viewport);
  glGetDoublev(GL_PROJECTION_MATRIX,projmatrix);
  glGetDoublev(GL_MODELVIEW_MATRIX,mvmatrix);

  double x = 0.5*(viewport[2]+viewport[0]);
  double y = 0.5*(viewport[3]+viewport[1]);
  realy = viewport[3] - (GLint)y - 1; /* See above */
  //realy = y;

  std::vector<Cartesian> ps;

  gluUnProject((GLdouble)x,(GLdouble)realy,0.0,mvmatrix,projmatrix,viewport,&wx,&wy,&wz);
  ps.push_back(Cartesian(wx,wy,wz));

  gluUnProject((GLdouble)x,(GLdouble)realy,1.0,mvmatrix,projmatrix,viewport,&wx,&wy,&wz);
  ps.push_back(Cartesian(wx,wy,wz));
  Cartesian p = ps[0] - ps[1];


  std::vector<Thing> projs;
  for(int ii=0;ii<new_ntriangles*3;ii+=3){
    double mid_x = new_vertices[3*new_indices[ii]] + new_vertices[3*new_indices[ii+1]] + new_vertices[3*new_indices[ii+2]];
    double mid_y = new_vertices[3*new_indices[ii]+1] + new_vertices[3*new_indices[ii+1]+1] + new_vertices[3*new_indices[ii+2]+1];
    double mid_z = new_vertices[3*new_indices[ii]+2] + new_vertices[3*new_indices[ii+1]+2] + new_vertices[3*new_indices[ii+2]+2];
    Cartesian cart(mid_x,mid_y,mid_z);
    projs.push_back(Thing(Cartesian::DotProduct(cart,p)/cart.length(),new_indices[ii],new_indices[ii+1],new_indices[ii+2]));
  }

  std::sort(projs.begin(),projs.end(),ZSortThing());
  std::vector<Thing>::iterator proj_iter=projs.begin();
  for(int ii=0;ii<new_ntriangles*3;ii+=3){
    new_indices[ii] = proj_iter->id1;
    new_indices[ii+1] = proj_iter->id2;
    new_indices[ii+2] = proj_iter->id3;
    ++proj_iter;
  }

}

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

bool surface::isLine() const {
  if(style==CCP4MG_SURFACE_MESH||style==CCP4MG_SURFACE_DOTS)
    return true;
  return false;
}

void surface::set_draw_colour(const GLfloat *col){}

void GetSimpleSubdividedDots(double *v1, double *v2, double *v3, double *col, double alpha, double dotSpacing, std::vector<Primitive*> &a, bool draw_first=true){
   GLdouble v12[3], v23[3], v31[3];
   GLint i;
   if(draw_first){
		PointElement *q = new PointElement(Cartesian(v1),col,Cartesian(v1),alpha);
		a.push_back(q);
		q = new PointElement(Cartesian(v2),col,Cartesian(v2),alpha);
		a.push_back(q);
		q = new PointElement(Cartesian(v3),col,Cartesian(v3),alpha);
		a.push_back(q);
   }
 
   bool do_12 = (v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2])>dotSpacing;
   bool do_13 = (v1[0]-v3[0])*(v1[0]-v3[0])+(v1[1]-v3[1])*(v1[1]-v3[1])+(v1[2]-v3[2])*(v1[2]-v3[2])>dotSpacing;
   bool do_23 = (v2[0]-v3[0])*(v2[0]-v3[0])+(v2[1]-v3[1])*(v2[1]-v3[1])+(v2[2]-v3[2])*(v2[2]-v3[2])>dotSpacing;
   
   if(do_12)
     for (i = 0; i < 3; i++) 
        v12[i] = 0.5*(v1[i]+v2[i]);
   if(do_13)
     for (i = 0; i < 3; i++)
        v31[i] = 0.5*(v3[i]+v1[i]);
   if(do_23)
     for (i = 0; i < 3; i++)
        v23[i] = 0.5*(v2[i]+v3[i]);

   if(do_12 && do_13)
     GetSimpleSubdividedDots(v1, v12, v31, col, alpha, dotSpacing, a, false);
   if(do_12 && do_23)
     GetSimpleSubdividedDots(v2, v23, v12, col, alpha, dotSpacing, a, false);
   if(do_13 && do_23)
     GetSimpleSubdividedDots(v3, v31, v23, col, alpha, dotSpacing, a, false);
   if(do_12 && do_13 && do_23)
     GetSimpleSubdividedDots(v12, v23, v31, col, alpha, dotSpacing, a);
}

void GetSimpleSubdividedDots(const Cartesian &v1, const Cartesian &v2, const Cartesian &v3, double *col, double alpha, double dotSpacing, std::vector<Primitive*> &a, bool draw_first=true){
 
   if(draw_first){
		PointElement *q = new PointElement(v1,col,v1,alpha);
		a.push_back(q);
		q = new PointElement(v2,col,v2,alpha);
		a.push_back(q);
		q = new PointElement(v3,col,v3,alpha);
		a.push_back(q);
   }
 
   bool do_12 = (v1-v2).length()>dotSpacing;
   bool do_13 = (v1-v3).length()>dotSpacing;
   bool do_23 = (v2-v3).length()>dotSpacing;

	Cartesian v12,v31,v23;
   
   if(do_12)
        v12 = 0.5*(v1+v2);
   if(do_13)
        v31 = 0.5*(v3+v1);
   if(do_23)
        v23 = 0.5*(v2+v3);

   if(do_12 && do_13)
     GetSimpleSubdividedDots(v1, v12, v31, col, alpha, dotSpacing, a, false);
   if(do_12 && do_23)
     GetSimpleSubdividedDots(v2, v23, v12, col, alpha, dotSpacing, a, false);
   if(do_13 && do_23)
     GetSimpleSubdividedDots(v3, v31, v23, col, alpha, dotSpacing, a, false);
   if(do_12 && do_13 && do_23)
     GetSimpleSubdividedDots(v12, v23, v31, col, alpha, dotSpacing, a);
}

std::vector<Primitive*> surface::GetRawDots() const {
	std::vector<Primitive*> a;
	double coords_1[4];
	double coords_2[4];
	double coords_3[4];
	double col[4];
	for (int i=0; i< theSurface->numberOfTriangles(); i++){
		// col is probably not correct, maybe even draw_dots has same problem.
		std::vector<Cartesian> carts;
		int id_1 = theSurface->vertex(i,0);
		int id_2 = theSurface->vertex(i,1);
		int id_3 = theSurface->vertex(i,2);

		theSurface->getCoord("colour", id_1, col);
		theSurface->getCoord("vertices", id_1, coords_1);
		PointElement *q = new PointElement(Cartesian(coords_1),col,Cartesian(coords_1),alpha);
		a.push_back(q);

		theSurface->getCoord("colour", id_2, col);
		theSurface->getCoord("vertices", id_2, coords_2);
		q = new PointElement(Cartesian(coords_2),col,Cartesian(coords_2),alpha);
		a.push_back(q);

		theSurface->getCoord("colour", id_3, col);
		theSurface->getCoord("vertices", id_3, coords_3);
		q = new PointElement(Cartesian(coords_3),col,Cartesian(coords_3),alpha);
		a.push_back(q);
	}
	return a;
}

int surface::GetNumberOfSimplePrimitives() const { 
   return theSurface->numberOfTriangles();
} ;

std::vector<Primitive*> surface::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start_in, int end_in) const {
	double coord[4];
	double normal[4];
	double col[4];

        int start = 0;
        int end = theSurface->numberOfTriangles();
        if(start_in>-1&&end_in>-1){
                start = start_in;
                end = end_in;
        }
	if(start>theSurface->numberOfTriangles()) start = theSurface->numberOfTriangles();
	if(end>theSurface->numberOfTriangles()) end = theSurface->numberOfTriangles();

	std::vector<Primitive*> a;
	if(style==CCP4MG_SURFACE_SOLID){
	for (int i=start; i< end; i++){
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
	}else if(style==CCP4MG_SURFACE_MESH){
	double col_1[4];
	double col_2[4];
	double col_3[4];
	for (int i=start; i< end; i++){
		std::vector<Cartesian> carts;
		std::vector<Cartesian> carts2(2);
		int id_1 = theSurface->vertex(i,0);
		int id_2 = theSurface->vertex(i,1);
		int id_3 = theSurface->vertex(i,2);
		theSurface->getCoord("colour", id_1, col_1);
		theSurface->getCoord("vertices", id_1, coord);
		carts.push_back(Cartesian(coord));
		theSurface->getCoord("colour", id_2, col_2);
		theSurface->getCoord("vertices", id_2, coord);
		carts.push_back(Cartesian(coord));
		theSurface->getCoord("colour", id_3, col_3);
		theSurface->getCoord("vertices", id_3, coord);
		carts.push_back(Cartesian(coord));
		Cartesian mid12 = 0.5*(carts[0]+carts[1]);
		Cartesian mid13 = 0.5*(carts[0]+carts[2]);
		Cartesian mid23 = 0.5*(carts[1]+carts[2]);
		carts2[0] = carts[0];
		carts2[1] = mid12;
		if((carts2[1]-carts2[0]).length()>1e-5){
			LineElement *q = new LineElement(carts2,col_1,carts2[0],alpha,textured);
			a.push_back(q);
		}
		carts2[0] = carts[1];
		if((carts2[1]-carts2[0]).length()>1e-5){
			LineElement *q = new LineElement(carts2,col_2,carts2[0],alpha,textured);
			a.push_back(q);
		}
		carts2[1] = mid23;
		if((carts2[1]-carts2[0]).length()>1e-5){
			LineElement *q = new LineElement(carts2,col_2,carts2[0],alpha,textured);
			a.push_back(q);
		}
		carts2[0] = carts[2];
		if((carts2[1]-carts2[0]).length()>1e-5){
			LineElement *q = new LineElement(carts2,col_3,carts2[0],alpha,textured);
			a.push_back(q);
		}
		carts2[1] = mid13;
		if((carts2[1]-carts2[0]).length()>1e-5){
			LineElement *q = new LineElement(carts2,col_3,carts2[0],alpha,textured);
			a.push_back(q);
		}
		carts2[0] = carts[0];
		if((carts2[1]-carts2[0]).length()>1e-5){
			LineElement *q = new LineElement(carts2,col_1,carts2[0],alpha,textured);
			a.push_back(q);
		}
	}	

	}else if(style==CCP4MG_SURFACE_DOTS){
	double coords_1[4];
	double coords_2[4];
	double coords_3[4];
	double dotSpacingSq = dotSpacing*dotSpacing;
	for (int i=start; i< end; i++){
		// col is probably not correct, maybe even draw_dots has same problem.
		std::vector<Cartesian> carts;
		int id_1 = theSurface->vertex(i,0);
		int id_2 = theSurface->vertex(i,1);
		int id_3 = theSurface->vertex(i,2);

		theSurface->getCoord("colour", id_1, col);
		theSurface->getCoord("vertices", id_1, coords_1);
		PointElement *q = new PointElement(Cartesian(coords_1),col,Cartesian(coords_1),alpha);
		a.push_back(q);

		theSurface->getCoord("colour", id_2, col);
		theSurface->getCoord("vertices", id_2, coords_2);
		q = new PointElement(Cartesian(coords_2),col,Cartesian(coords_2),alpha);
		a.push_back(q);

		theSurface->getCoord("colour", id_3, col);
		theSurface->getCoord("vertices", id_3, coords_3);
		q = new PointElement(Cartesian(coords_3),col,Cartesian(coords_3),alpha);
		a.push_back(q);
		if(!theSurface->getCoord("vertices", id_1, coords_1)&&!theSurface->getCoord("vertices", id_2, coords_2)&&!theSurface->getCoord("vertices", id_3, coords_3)){
			//GetSimpleSubdividedDots(Cartesian(coords_1),Cartesian(coords_2),Cartesian(coords_3), col, alpha,dotSpacing,a);
			GetSimpleSubdividedDots(coords_1,coords_2,coords_3, col, alpha,dotSpacingSq,a);
		}

	}
	}else{
	 	printf("Unknown surface style in surface::GetSimplePrimitives\n");
	}
	return a;
}

void surface::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
	
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
		p = quat.getMatrix()*(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])+objorigin)+Cartesian(ox,oy,oz));
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
		p = quat.getMatrix()*(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])));
		p.normalize();
		fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
	}
	theSurface->getCoord("normals", i, coord_1);
	p = quat.getMatrix()*(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])));
	p.normalize();
	fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
	fp << "  }\n";
	
	fp << "  texture_list {\n";
	fp << "  " << theSurface->numberOfVertices() << ",\n";
	for (i=0; i< theSurface->numberOfVertices()-1; i++){
		theSurface->getCoord("colour", i, coord_1);
                if(get_transparent()!=0)
		  fp << "  texture{pigment{rgbt< " << coord_1[0] << ", " << coord_1[1] << ", " << coord_1[2]  << ", " << 1.0-alpha << ">}finish {diffuse 1.0 specular 1.0}},\n";
                else
		  fp << "  texture{pigment{rgb< " << coord_1[0] << ", " << coord_1[1] << ", " << coord_1[2] << ">}finish {diffuse 1.0 specular 1.0}},\n";
	}
	theSurface->getCoord("colour", i, coord_1);
        if(get_transparent()!=0)
	  fp << "  texture{pigment{rgbt< " << coord_1[0] << ", " << coord_1[1] << ", " << coord_1[2] << ", " << 1.0-alpha << ">}finish {diffuse 1.0 specular 1.0}}\n";
        else
	  fp << "  texture{pigment{rgb< " << coord_1[0] << ", " << coord_1[1] << ", " << coord_1[2] << ">}finish {diffuse 1.0 specular 1.0}}\n";
	fp << "  }\n";
	
	fp << "  face_indices {\n";
        int ninview = 0;
        for (i=0; i<theSurface->numberOfTriangles(); i++){
	  id_1 = theSurface->vertex(i,0);
	  id_2 = theSurface->vertex(i,1);
	  id_3 = theSurface->vertex(i,2);
          theSurface->getCoord("vertices", id_1, coord_1);
          if(v.PointInVolume(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])+objorigin)))
	    ninview++;
        }
        if(theSurface->numberOfTriangles()>0)
	  fp << "  " << ninview+1 << ",\n";
        else
	  fp << "  0,\n";
	for (i=0; i< theSurface->numberOfTriangles(); i++){
		id_1 = theSurface->vertex(i,0);
		id_2 = theSurface->vertex(i,1);
		id_3 = theSurface->vertex(i,2);
          theSurface->getCoord("vertices", id_1, coord_1);
          if(v.PointInVolume(objrotmatrix*(Cartesian(coord_1[0],coord_1[1],coord_1[2])+objorigin)))
		fp << "  < " << id_1 << ", " << id_2 << ", " << id_3 << ">," << id_1 << "," << id_2 << "," << id_3 << ",\n";
	}
        // Dummy triangle, inspection of above should explain
        if(theSurface->numberOfTriangles()>0) fp << "  < 0,0,0 >,\n";
	fp << "  }\n";
	
	fp << "}\n";
	
}

void surface::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){
  std::vector<Primitive*> ps_prims = GetSimplePrimitives(v,objrotmatrix,objorigin);
  std::cout << "Got " << ps_prims.size() << " simple surface primitives\n";
  std::vector<Primitive*>::iterator a_iter = ps_prims.begin();
  while(a_iter!=ps_prims.end()){
    (*a_iter)->DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
    a_iter++;
  }

}


surface::~surface (){
	delete theSurface;
	freeResources();
}
surface::surface (): Primitive() {
  dotSpacing = 0.2;
  dotSize = 2.0;
  lineWidth = 1.0;
  initArrays();
  style = CCP4MG_SURFACE_SOLID;
  theSurface = new CXXSurface();
}

surface::surface (mmdb::Manager* *theManager, int selHnd, double probe_radius, double delta, bool blend_edges) : Primitive() {
//	std::cout << " EP 2"; std::cout.flush();
        dotSpacing = 0.2;
        dotSize = 2.0;
	initArrays();
	style = CCP4MG_SURFACE_SOLID;
	theSurface = new CXXSurface();
	theSurface->calculateFromAtoms (theManager, selHnd, selHnd,  probe_radius, delta,blend_edges);
	//evaluateElectrostaticPotential(theManager, selHnd);
	theSurface->report();
}

surface::surface (mmdb::Manager* *theManager, int selHnd, int contextSelHnd, double probe_radius, double delta, bool blend_edges) : Primitive() {
//	std::cout << " EP 3"; std::cout.flush();
        dotSpacing = 0.2;
        dotSize = 2.0;
	initArrays();
	style = CCP4MG_SURFACE_SOLID;
	theSurface = new CXXSurface();
	theSurface->calculateFromAtoms (theManager, selHnd, contextSelHnd, probe_radius, delta,blend_edges);
}

surface::surface (const std::string &fileName) : Primitive(){
//	std::cout << " EP 4"; std::cout.flush();
        dotSpacing = 0.2;
        dotSize = 2.0;
	style = CCP4MG_SURFACE_SOLID;
	theSurface = new CXXSurface (fileName);
}

int surface::calculate (mmdb::Manager* *theManager, int selHnd, double probe_radius, double delta, bool blend_edges)  {
	return theSurface->calculateFromAtoms (theManager, selHnd, selHnd,  probe_radius, delta,blend_edges);
}

int surface::calculate_with_context (mmdb::Manager* *theManager, int selHnd, int contextSelHnd, double probe_radius, double delta, bool blend_edges) {
	return  theSurface->calculateFromAtoms (theManager, selHnd, contextSelHnd, probe_radius, delta,blend_edges);
}

/*
void surface::SubdivideDotsLine(double *v1, double *v2, bool draw_first){
   if(draw_first) glVertex3dv(v1);
   glVertex3dv(v2);
   if((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2])>dotSpacing){
     GLdouble v12[3];
     for(int i = 0; i < 3; i++) {
        v12[i] = 0.5*(v1[i]+v2[i]);
     }
     SubdivideDotsLine(v1,v12,false);
     SubdivideDotsLine(v2,v12,false);
   }
}
*/
void SubdivideDotsGL(GLfloat *v1, GLfloat *v2, GLfloat *v3, double dotSpacingSquared, bool draw_first=true){
   GLfloat v12[3], v23[3], v31[3];
   GLint i;
 
   if(draw_first) {
	glVertex3fv(v1);
   	glVertex3fv(v2);
   	glVertex3fv(v3);
   }
 
   bool do_12 = (v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2])>dotSpacingSquared;
   bool do_13 = (v1[0]-v3[0])*(v1[0]-v3[0])+(v1[1]-v3[1])*(v1[1]-v3[1])+(v1[2]-v3[2])*(v1[2]-v3[2])>dotSpacingSquared;
   bool do_23 = (v2[0]-v3[0])*(v2[0]-v3[0])+(v2[1]-v3[1])*(v2[1]-v3[1])+(v2[2]-v3[2])*(v2[2]-v3[2])>dotSpacingSquared;
   
   if(do_12)
     for (i = 0; i < 3; i++) 
        v12[i] = 0.5*(v1[i]+v2[i]);
   if(do_13)
     for (i = 0; i < 3; i++)
        v31[i] = 0.5*(v3[i]+v1[i]);
   if(do_23)
     for (i = 0; i < 3; i++)
        v23[i] = 0.5*(v2[i]+v3[i]);

   if(do_12 && do_13)
     SubdivideDotsGL(v1, v12, v31, dotSpacingSquared, false);
   if(do_12 && do_23)
     SubdivideDotsGL(v2, v23, v12, dotSpacingSquared,false);
   if(do_13 && do_23)
     SubdivideDotsGL(v3, v31, v23, dotSpacingSquared,false);
   if(do_12 && do_13 && do_23)
     SubdivideDotsGL(v12, v23, v31, dotSpacingSquared);
   
}

void surface::draw_dots(const double *override_colour, int selective_override){
  glDisable(GL_LIGHTING);
  double coords_1[4];
  double coords_2[4];
  double coords_3[4];
  glPointSize(dotSize);

  double dotSpacingSq = dotSpacing*dotSpacing;
  int idx = 0;
  glBegin(GL_POINTS);
  for (int i=0; i< theSurface->numberOfTriangles(); i++){
      glColor3fv(colors+4*indices[idx]);
      SubdivideDotsGL(vertices+3*indices[idx],vertices+3*indices[idx+1],vertices+3*indices[idx+2],dotSpacingSq);
      idx+=3;
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
        //std::cout << "OPAQUENESS: " << opaqueness << "\n";
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
		colors[4*i+3] = alpha;
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
	if (useVBO&&alpha>0.99) {
		if (!arraysUploaded) blatBuffer(GL_ARRAY_BUFFER, cardVertexBuffer, 3*nVerts*sizeof(GLfloat), (void **) &vertices, GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBufferARB(GL_ARRAY_BUFFER,cardVertexBuffer);
#endif
	}
#ifdef UNDERSTANDS_VBOs
	if (useVBO&&alpha>0.99) {
	  glVertexPointer(3, GL_FLOAT, 0, NULL);
        } else {
	  glVertexPointer(3, GL_FLOAT, 0, vertices);
        }
#else
	glVertexPointer(3, GL_FLOAT, 0, vertices);
#endif
	
	glEnableClientState(GL_NORMAL_ARRAY);
	if (useVBO&&alpha>0.99) {
		if (!arraysUploaded) blatBuffer(GL_ARRAY_BUFFER, cardNormalBuffer, 3*nVerts*sizeof(GLfloat), (void **) &normals, GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ARRAY_BUFFER,cardNormalBuffer);
#endif
	}
#ifdef UNDERSTANDS_VBOs
	if (useVBO&&alpha>0.99) {
	  glNormalPointer(GL_FLOAT, 0, NULL);
        } else {
	  glNormalPointer(GL_FLOAT, 0, normals);
        }
#else
	glNormalPointer(GL_FLOAT, 0, normals);
#endif
	
	glEnableClientState(GL_COLOR_ARRAY);
	if (useVBO&&alpha>0.99) {
		if (!arraysUploaded) blatBuffer(GL_ARRAY_BUFFER, cardColorBuffer,  4*nVerts*sizeof(GLfloat), (void **) &colors,  GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ARRAY_BUFFER,cardColorBuffer);
#endif
	}
#ifdef UNDERSTANDS_VBOs
	if (useVBO&&alpha>0.99) {
	  glColorPointer(4, GL_FLOAT, 0, NULL);
        }else{
	  glColorPointer(4, GL_FLOAT, 0, colors);
        }
#else
	glColorPointer(4, GL_FLOAT, 0, colors);
#endif
	
	if (useVBO&&alpha>0.99) {
		if (!arraysUploaded) blatBuffer(GL_ELEMENT_ARRAY_BUFFER, cardIndexBuffer, 3*mySurface.numberOfTriangles()*sizeof(GLuint), (void ** ) &indices,  GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cardIndexBuffer);
#endif
	}
	arraysUploaded = 1;
	return;
}	


void surface::draw(const double *override_colour, int selective_override){
        glEnable(GL_LIGHTING); 
#ifdef _WIN32
	if(!have_range_ext)
          init_range_ext();
#endif
	if(style==CCP4MG_SURFACE_MESH){
          glDisable(GL_LIGHTING);
          glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
          glLineWidth(lineWidth);
        }
	if (!arraysGenerated) generateArrays();
	if(style==CCP4MG_SURFACE_DOTS){
		draw_dots(override_colour,selective_override);
		return;
	}
	CXXSurface &mySurface = *theSurface;
	
        // At this point we should sort indices....
        if(alpha<0.99&&style==CCP4MG_SURFACE_SOLID) SortIndices(vertices,indices,mySurface.numberOfTriangles());
	bindArrays();
	
	float whiteColour[] = {1.,1.,1.,1.};
	float blackColor[] = {0.,0.,0.,1.};
	
	GLenum err;
	
	if (alpha < 1.&&style!=CCP4MG_SURFACE_MESH) {
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);		
#ifdef  __APPLE_CC__
                glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
#endif
	}
	
	//Let the returned colour dictate: note obligatory order of these calls
	//glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	//glEnable(GL_COLOR_MATERIAL);
	
	//Set material properties that are not per-vertex
	float specularColor[4];
	for (int i=0; i<4; i++) specularColor[i] = whiteColour[i] * specularity;
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blackColor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128); 
#ifdef UNDERSTANDS_VBOs
  if(!useVBO||alpha<0.99){
#endif
	GLint max_verts;
	glGetIntegerv(GL_MAX_ELEMENTS_VERTICES,&max_verts);
	max_verts = (max_verts/3)*3;
	int ndrawn = 0;
	if(3*mySurface.numberOfTriangles()>max_verts){
	  int nloops = (3*mySurface.numberOfTriangles()-1)/max_verts + 1;
	  for(int iloop=0;iloop<nloops;iloop++){
	    int ndraw = max_verts;
	    if((iloop+1)*max_verts>3*mySurface.numberOfTriangles()){
		    ndraw =3*mySurface.numberOfTriangles() - (iloop*max_verts);
	    }
	    glDrawElements(GL_TRIANGLES,ndraw, GL_UNSIGNED_INT, indices+(iloop*max_verts));
	    ndrawn += ndraw;
	  }
	} else {
	  glDrawElements(GL_TRIANGLES, 3*mySurface.numberOfTriangles(), GL_UNSIGNED_INT, indices);
	}
#ifdef UNDERSTANDS_VBOs
  } else {
     glDrawElements(GL_TRIANGLES, 3*mySurface.numberOfTriangles(), GL_UNSIGNED_INT, 0);
  }
#endif

#ifdef UNDERSTANDS_VBOs
  if(useVBO&&alpha>0.99){
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
  }
#endif
	if(style!=CCP4MG_SURFACE_MESH) return;

	//if (alpha < 1.0&&style!=CCP4MG_SURFACE_MESH)
		//glBlendFunc (GL_ONE, GL_ZERO);
	if(style==CCP4MG_SURFACE_MESH){
          glPolygonMode(GL_FRONT,GL_FILL);
          glPolygonMode(GL_BACK,GL_FILL);
        }
	
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

void surface::forceRegenerateArrays(){
  arraysGenerated = 0;
}

void surface::drawElement(int element){
#ifdef _WIN32
	if(!have_range_ext)
          init_range_ext();
#endif
	CXXSurface &mySurface = *theSurface;
	if (useVBO&&alpha>0.99) {
		if (!arraysGenerated) generateArrays();
		//bindArrays();
		GLvoid *offset = BUFFER_OFFSET(element*3*sizeof(GLuint));
		glDrawRangeElements(GL_TRIANGLES, 0, mySurface.numberOfVertices(), 1, GL_UNSIGNED_INT , offset );
		return;
	}
	
	double coords[4];
	float whiteColour[] = {1.,1.,1.,1.};
	float blackColor[] = {0.,0.,0.,1.};
	
	if (alpha < 1.&&style!=CCP4MG_SURFACE_MESH){
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
#ifdef  __APPLE_CC__
                glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE);
#endif
        }
	
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
			coords[3] = alpha / 255;
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
	if (useVBO&&alpha>0.99){
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
	obj.add_surf_primitive(surf);
}


void surface::ColourSurface(mmdb::Manager* *molHnd, int selHnd, const AtomColourVector &atomColourVector){
	int nSelAtoms;
	mmdb::Atom** SelAtom;
	molHnd->GetSelIndex(selHnd, SelAtom, nSelAtoms);
	
	std::vector < double*> colours;
	for(unsigned i=0;i<unsigned(nSelAtoms);i++){
		colours.push_back((double *)atomColourVector.GetRGB(i));
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

int surface::evaluatePhiAndColourWithDefaultScheme(mmdb::Manager* *theManager, const int selHnd,  const int context_selHnd, int contains_hydrogen ){
	CColourScheme defaultScheme;
	std::vector<float> typ;
	typ.push_back(-0.2); typ.push_back(-0.2);  typ.push_back(0.0);  typ.push_back(0.2); typ.push_back(0.2);
	std::vector<std::string> cols;
	cols.push_back("red"); cols.push_back("red");  cols.push_back("white");  cols.push_back("blue"); cols.push_back("blue");
	defaultScheme.SetSchemeFloat(typ, cols);
	return evaluatePhiAndColourWithScheme(theManager, selHnd, context_selHnd, defaultScheme,contains_hydrogen );
}


int surface::evaluatePhiAndColourWithScheme(mmdb::Manager* *theManager, const int selHnd, const int context_selHnd, CColourScheme &colourScheme , int contains_hydrogen) {
	//Instantiate and calculate electrostatic potential
	// contains_hydrogen = 0 => no hydrogen
	// contains_hydrogen = 0 => some H atoms in structure (currently no checks if all present)
	CXXCreator theCreator(theManager, selHnd, context_selHnd);
	theCreator.calculate();	
	
    
	//Coerce map into clipper NXmap
	clipper::Cell aCell;
	clipper::NXmap<double> thePhiMap (theCreator.coerceToClipperMap(aCell));
	//writeNXMap(thePhiMap,"/tmp/lizp/phi.map");
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
                        if(coords&&coords[0]&&coords[1]&&coords[2]){
			  clipper::Coord_orth orthogonals(coords[0], coords[1], coords[2]);
			  const clipper::Coord_map mapUnits(aMap.coord_map(orthogonals));
			  scaluffer[i] = aMap.interp<clipper::Interp_cubic>( mapUnits );
                        }
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
			CXXCoord<CXXCoord_ftype>colour (newColour[0], newColour[1], newColour[2]);
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
	
	clipper::Cell cell(clipper::Cell_descr(1.0,1.0,1.0,90.0,90.0,90.0));
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
	useVBO=1;
	arraysGenerated = 0;
	arraysUploaded = 0;
	specularity = 1.;
	exponent = 128.;
	opaqueness = 1.;
	alpha = 1.;
}

 int surface::writeGraspFile (std::string map_file_name) {
   return theSurface->writeAsGrasp(map_file_name);
 }

 int surface::readGraspFile (std::string map_file_name) {
   return theSurface->readGraspFile(map_file_name);
 }
