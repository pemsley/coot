/*
     pygl/cprimitive.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2011 University of York
     Copyright (C) 2012 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/
#if defined (_WIN32)
#include <windows.h>
//#include "glew.h"
//#include "wglew.h"
#else
//#include "glew.h"
#endif

#if defined (linux)
#undef GLX_GLXEXT_LEGACY
#define GL_GLEXT_PROTOTYPES
#endif

#include "cprimitive.h"
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#endif
#ifdef _USE_GLUT_
#include <GL/glut.h>
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "ppmutil.h"
#include "texture_globals.h"
#include "sphere.h"
#include "subdivide.h"

#include "matrix.h"
#include "plane.h"
#include "geomutil.h"
#include "mgutil.h"

#include "catmull.h"

#if !defined (_WIN32)
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#endif

#ifdef GL_ARRAY_BUFFER
#define UNDERSTANDS_VBOs
#else
#define GL_ARRAY_BUFFER 0x8892
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_STATIC_DRAW 0x88E4
#undef UNDERSTANDS_VBOs
#endif

#define BUFFER_OFFSET(offset) ((char *)NULL + (offset))

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#define PIBY2 (M_PI * 2)
#define ROOT3 1.7320508075688772
#define ROOT3OVER2 ROOT3*0.5

double SimpleBillBoard::projMatrix[16];
double SimpleBillBoard::modelMatrix[16];
GLint SimpleBillBoard::viewport[4];

int RenderQuality::render_quality;
int RenderQuality::ribbon_quality;

/*
std::vector<unsigned> factors(unsigned x){
  std::vector<unsigned> facs;
  for(unsigned i=2;i<x/2;i++){
     if(x%i==0){
        facs.push_back(i);
     }
   }
  return facs;
}

unsigned nearest(unsigned x,unsigned yy){
  std::vector<unsigned> y = factors(yy);
  if(y.size()<1) return x;
  unsigned nval = INT_MAX;
  unsigned near = y[0];
  for(unsigned ii=0;ii<y.size();ii++){
    unsigned v = y[ii];
    if(abs(x-v)<nval){
      nval = abs(x-v);
      near = v;
    }
  }
  return near;
}
*/

GLint* SimpleBillBoard::GetViewport(){
  return viewport;
}

GLdouble* SimpleBillBoard::GetModelviewMatrix(){
  return modelMatrix;
}

GLdouble* SimpleBillBoard::GetProjectionMatrix(){
  return projMatrix;
}

void SimpleBillBoard::SetMatrices(){
  glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
  glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
  glGetIntegerv(GL_VIEWPORT,viewport);
}

int SimpleBillBoard::GetMagnification(){
  double curr_projMatrix[16];
  glGetDoublev(GL_PROJECTION_MATRIX,curr_projMatrix);
  double curr_l = (curr_projMatrix[12] - 1.0)/curr_projMatrix[0]; 
  double curr_r = (curr_projMatrix[12] + 1.0)/curr_projMatrix[0];
  double *init_projMatrix = SimpleBillBoard::GetProjectionMatrix();
  double init_l = (init_projMatrix[12] - 1.0)/init_projMatrix[0]; 
  double init_r = (init_projMatrix[12] + 1.0)/init_projMatrix[0];
  int scale = int((init_r-init_l)/(curr_r-curr_l));
  return scale;
}

void draw_quadstrip_textured(const double *vertices, int nvertices, int textured){
  int i,k=0;
  double frac=0.0;

  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUAD_STRIP);

  for(i=0;i<nvertices;i++){
    if(i%2 == 0){
      frac = (double)i/(double)(nvertices-2);
      set_texture_coord(frac, 0.0, textured);
    }else{
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3f(vertices[k],vertices[k+1],vertices[k+2]);
    k = k + 3;
  }

  glEnd();
  glDisable(GL_TEXTURE_2D);
}

void draw_quadstrip_outline(const double *vertices, int nvertices){
  int i,k;

  k=0;
  glLineWidth(1.5);
  glPolygonMode(GL_FRONT,GL_LINE);
  glBegin(GL_LINE_STRIP);
  glColor3f(0.0,0.0,0.0);
  for(i=0;i<nvertices;i++){

    if(i%2 == 0)    
      glVertex3f(vertices[k],vertices[k+1],vertices[k+2]);
    k = k + 3;
  }
  glEnd();
  glPolygonMode(GL_FRONT,GL_FILL);

}

std::vector<std::vector<Cartesian> > GetPolyCylinderAB(const std::vector<Cartesian> &vertices, double size){
  unsigned int i;
  Cartesian line; 

  std::vector<Cartesian> a,b;
  Cartesian atmp, btmp;
  Cartesian xaxis(1,0,0);
  Cartesian yaxis(0,1,0);
  Cartesian zaxis(0,0,1);

  for(i=0;i<vertices.size()-1;i++){
    line = vertices[i+1] - vertices[i];

     atmp = line.CrossProduct(line,zaxis);
     if(atmp.length() < 0.000000001){
       atmp = line.CrossProduct(line,yaxis);
       if(atmp.length() < 0.000000001){
	 atmp = line.CrossProduct(line,xaxis);
         if(atmp.length()  > 0.000000001){
	   atmp.normalize(size);
	   a.push_back(atmp);
	   btmp = line.CrossProduct(line,atmp);
	   btmp.normalize(size);
	   b.push_back(btmp);
	 }
       }else{
	   atmp.normalize(size);
	   a.push_back(atmp);
	   btmp = line.CrossProduct(line,atmp);
	   btmp.normalize(size);
	   b.push_back(btmp);
       }
     }else{
	   atmp.normalize(size);
	   a.push_back(atmp);
	   btmp = line.CrossProduct(line,atmp);
	   btmp.normalize(size);
	   b.push_back(btmp);
     }
     
  }
  a.push_back(a.back());
  b.push_back(b.back());
  std::vector<std::vector<Cartesian> > AB;
  AB.push_back(a);
  AB.push_back(b);
  AB.push_back(a);
  AB.push_back(b);
  return AB;

}

void draw_billboard(double v11, double v12, double v21, double v22){

  glEnable(GL_TEXTURE_2D);
  glPushMatrix();
  glLoadIdentity();
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  GLdouble projmatrix[16];
  glGetDoublev(GL_PROJECTION_MATRIX,projmatrix);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0f, 0.0f);  
  glVertex3f((2*v11-1)/projmatrix[5], (2*v12-1)/projmatrix[5], 0.0f);
  glTexCoord2f(1.0f, 0.0f);  
  glVertex3f((2*v21-1)/projmatrix[5], (2*v12-1)/projmatrix[5], 0.0f);
  glTexCoord2f(1.0f, 1.0f);  
  glVertex3f((2*v21-1)/projmatrix[5], (2*v22-1)/projmatrix[5], 0.0f);
  glTexCoord2f(0.0f, 1.0f);  
  glVertex3f((2*v11-1)/projmatrix[5], (2*v22-1)/projmatrix[5], 0.0f);
  glEnd();
  glPopMatrix();
  glDisable(GL_TEXTURE_2D);
}

void draw_capped_cylinder(const std::vector<Cartesian> &vertices, double size,int accu, int textured, bool force_dl){
  draw_cylinder(vertices,size,accu,textured,force_dl,true);
}

void draw_cylinder(const std::vector<Cartesian> &vertices, double size, int accu_in, int textured, bool force_dl, bool capped){

  static bool have_dl_4 = false;
  static bool have_dl_8 = false;
  static bool have_dl_16 = false;
  static bool have_dl_32 = false;
  static GLuint listid_4;
  static GLuint listid_8;
  static GLuint listid_16;
  static GLuint listid_32;
  static GLuint listid_cap_4;
  static GLuint listid_cap_8;
  static GLuint listid_cap_16;
  static GLuint listid_cap_32;

  int accu;
  if(RenderQuality::GetRenderQuality())
    accu = 32;
  else
    accu = accu_in;

  if(((!have_dl_4)||force_dl)&&accu==4){
    GLUquadric *q = gluNewQuadric();
    if(listid_4==0) listid_4 = glGenLists(1);
    if(listid_cap_4==0) listid_cap_4 = glGenLists(1);
    glNewList(listid_4,GL_COMPILE);
    gluCylinder(q,1.0,1.0,1.0,4,1);
    glEndList();
    glNewList(listid_cap_4,GL_COMPILE);
    gluDisk(q,0.0,1.0,4,1);
    glEndList();
    have_dl_4 = true;
    gluDeleteQuadric(q);
  }
  if(((!have_dl_8)||force_dl)&&accu==8){
    GLUquadric *q = gluNewQuadric();
    if(listid_8==0) listid_8 = glGenLists(1);
    if(listid_cap_8==0) listid_cap_8 = glGenLists(1);
    glNewList(listid_8,GL_COMPILE);
    gluCylinder(q,1.0,1.0,1.0,8,1);
    glEndList();
    glNewList(listid_cap_8,GL_COMPILE);
    gluDisk(q,0.0,1.0,8,1);
    glEndList();
    have_dl_8 = true;
    gluDeleteQuadric(q);
  }
  if(((!have_dl_16)||force_dl)&&accu==16){
    GLUquadric *q = gluNewQuadric();
    if(listid_16==0) listid_16 = glGenLists(1);
    if(listid_cap_16==0) listid_cap_16 = glGenLists(1);
    glNewList(listid_16,GL_COMPILE);
    gluCylinder(q,1.0,1.0,1.0,16,1);
    glEndList();
    glNewList(listid_cap_16,GL_COMPILE);
    gluDisk(q,0.0,1.0,16,1);
    glEndList();
    have_dl_16 = true;
    gluDeleteQuadric(q);
  }
  if(((!have_dl_32)||force_dl)&&accu==32){
    GLUquadric *q = gluNewQuadric();
    if(listid_32==0) listid_32 = glGenLists(1);
    if(listid_cap_32==0) listid_cap_32 = glGenLists(1);
    glNewList(listid_32,GL_COMPILE);
    gluCylinder(q,1.0,1.0,1.0,32,1);
    glEndList();
    glNewList(listid_cap_32,GL_COMPILE);
    gluDisk(q,0.0,1.0,32,1);
    glEndList();
    have_dl_32 = true;
    gluDeleteQuadric(q);
  }


  Cartesian line;
  bool swapped = false;
  Cartesian v1 = vertices[0];
  Cartesian v2 = vertices[1];
  if(v1.length()>v2.length()){
    line = v2 - v1;
  }else{
    line = v1 - v2;
    swapped = true;
  }

  double ax;
  double rx;
  double ry;
  double length = line.length();
  double vz = line.get_z();

  bool rot_x = false;
  if(fabs(vz)>1e-7){
    ax = 180.0/M_PI*acos(vz/length);
    if(vz<0.0) ax = -ax;
    rx = -line.get_y()*vz;
    ry = line.get_x()*vz;
  }else{
    double vx = line.get_x();
    double vy = line.get_y();
    ax = 180.0/M_PI*acos(vx/length);
    if(vy<0) ax = -ax;
    rot_x = true;
  }

  glPushMatrix();
  glTranslatef(v1.get_x(),v1.get_y(),v1.get_z());
  if(rot_x){
    glRotated(90.0, 0.0, 1.0, 0.0);
    glRotated(ax,  -1.0, 0.0, 0.0);
  }else{
    glRotated(ax, rx, ry, 0.0);
    if(ax<0.0){
      glScalef(-1,1,-1);
      glRotated(180, -rx, ry, 0.0);
    }
  }
  if(swapped)
    glScalef(size,-size,-length);
  else
    glScalef(size,size,length);

  if(accu==4){
    glCallList(listid_4);
    if(capped){
      glTranslatef(0.0f,0.0f,1.0f);
      glCallList(listid_cap_4);
      glScalef(1.0f,-1.0f,-1.0f);
      glTranslatef(0.0f,0.0f,1.0f);
      glCallList(listid_cap_4);
    }
  }else if(accu==8){
    glCallList(listid_8);
    if(capped){
      glTranslatef(0.0f,0.0f,1.0f);
      glCallList(listid_cap_8);
      glScalef(1.0f,-1.0f,-1.0f);
      glTranslatef(0.0f,0.0f,1.0f);
      glCallList(listid_cap_8);
    }
  }else if(accu==16){
    glCallList(listid_16);
    if(capped){
      glTranslatef(0.0f,0.0f,1.0f);
      glCallList(listid_cap_16);
      glScalef(1.0f,-1.0f,-1.0f);
      glTranslatef(0.0f,0.0f,1.0f);
      glCallList(listid_cap_16);
    }
  }else if(accu==32){
    glCallList(listid_32);
    if(capped){
      glTranslatef(0.0f,0.0f,1.0f);
      glCallList(listid_cap_32);
      glScalef(1.0f,-1.0f,-1.0f);
      glTranslatef(0.0f,0.0f,1.0f);
      glCallList(listid_cap_32);
    }
  }else{
    GLUquadric *q = gluNewQuadric();
    gluCylinder(q,1.0,1.0,1.0,accu,1);
    if(capped){
      glTranslatef(0.0f,0.0f,1.0f);
      gluDisk(q,0.0,1.0,accu,1);
      glScalef(1.0f,-1.0f,-1.0f);
      glTranslatef(0.0f,0.0f,1.0f);
      gluDisk(q,0.0,1.0,accu,1);
    }
    gluDeleteQuadric(q);
  }
  glPopMatrix();
}

void draw_ellipse_point(double theta, const Cartesian &pv, const Cartesian &pvpr, const Cartesian &spline);

std::vector<Cartesian> GetConePointAndNormal(double theta, const Cartesian &pv, const Cartesian &pvpr, const Cartesian &spline){

  double a = pv.length();
  double b = pvpr.length();

  double x = cos(theta);
  double y = sin(theta);

  Cartesian p = x*pv + y*pvpr + spline;
  Cartesian n = pv*(x/(a*a)) + pvpr*(y/(b*b));
  n.normalize();

  std::vector<Cartesian> point_and_normal;
  point_and_normal.push_back(p);
  point_and_normal.push_back(n);

  return point_and_normal;


}


void draw_cone(const Cartesian &v1, const Cartesian &v2, double size,int accu, int textured){

  Cartesian line;
  Cartesian a, b;
  Cartesian xaxis(1,0,0);
  Cartesian yaxis(0,1,0);
  Cartesian zaxis(0,0,1);

  line = v2 - v1;

  a = line.CrossProduct(line,zaxis);
  if(a.length() < 0.000000001){
    a = line.CrossProduct(line,yaxis);
    if(a.length() < 0.000000001){
      a = line.CrossProduct(line,xaxis);
      if(a.length()  < 0.000000001){ 
	      fprintf(stderr,"Error cannot find suitable axis. Coincident points?\n");
      }
    }
  }
  a.normalize();

  b = line.CrossProduct(line,a);

  a.normalize(size);
  b.normalize(size);

  std::vector<Cartesian> vectors;
  std::vector<Cartesian> pvpr;
  std::vector<Cartesian> pv;

  vectors.push_back(v1);
  vectors.push_back(v2);

  pv.push_back(a);
  pv.push_back(a);
  pvpr.push_back(b);
  pvpr.push_back(b);

  std::vector<Cartesian> dummy;

  double theta,theta2;

  glBegin(GL_TRIANGLES);

  Cartesian line_n = line;
  line_n.normalize();

  for(int i=0;i<360;i=i+360/accu){
    theta = (double)i/360.0 * PIBY2;
    theta2 = (double)(i-360/accu)/360.0 * PIBY2;
    
    std::vector<Cartesian> pn1 = GetConePointAndNormal(theta, pv[0], pvpr[0], vectors[0]);
    std::vector<Cartesian> pn2 = GetConePointAndNormal(theta2, pv[0], pvpr[0], vectors[0]);

    double phi = atan(line.length());

    Cartesian n1 = sin(phi)*pn1[1] + cos(phi)*line_n;
    Cartesian n2 = sin(phi)*pn2[1] + cos(phi)*line_n;

    glNormal3f(n1.get_x(),n1.get_y(),n1.get_z());
    glVertex3f(pn1[0].get_x(),pn1[0].get_y(),pn1[0].get_z());

    Cartesian tip_n = 0.5*n1 + 0.5*n2;
    glNormal3f(tip_n.get_x(),tip_n.get_y(),tip_n.get_z());
    glVertex3f(vectors[1].get_x(),vectors[1].get_y(),vectors[1].get_z());

    glNormal3f(n2.get_x(),n2.get_y(),n2.get_z());
    glVertex3f(pn2[0].get_x(),pn2[0].get_y(),pn2[0].get_z());

    glNormal3f(-line_n.get_x(),-line_n.get_y(),-line_n.get_z());
    glVertex3f(vectors[0].get_x(),vectors[0].get_y(),vectors[0].get_z());
    glVertex3f(pn1[0].get_x(),pn1[0].get_y(),pn1[0].get_z());
    glVertex3f(pn2[0].get_x(),pn2[0].get_y(),pn2[0].get_z());
  }
  glEnd();
  
}


void set_line_colour(double r, double g, double b, double a){
  GLfloat colour[4] = { r,g,b,a };
  glColor4fv(colour);
	
}

void draw_line(double v11, double v12, double v13, double v21, double v22, double v23, double size){
  glLineWidth(size);
  glBegin(GL_LINES);
  glVertex3f(v11,v12,v13);
  glVertex3f(v21,v22,v23);
  glEnd();
   
}

void draw_sphere(double v1, double v2, double v3, int accu_in, double size){
  int accu;
  if(RenderQuality::GetRenderQuality())
    accu = 4;
  else
    accu = accu_in;
  glPushMatrix();
  glTranslatef(v1,v2,v3);
  sphere(accu,size);
  glPopMatrix();
}

Lighting::Lighting(){
  ambient[0] = ambient[1] = ambient[2] = 0.0;
  specular[0] = specular[1] = specular[2] = 1.0;
  diffuse[0] = diffuse[1] = diffuse[2] = 1.0;
  emission[0] = emission[1] = emission[2] = 0.0;
  shininess = 128.0;
}

void Lighting::copy(const Lighting &light_in){
  ambient[0] = light_in.ambient[0];
  ambient[1] = light_in.ambient[1];
  ambient[2] = light_in.ambient[2];
  specular[0] = light_in.specular[0];
  specular[1] = light_in.specular[1];
  specular[2] = light_in.specular[2];
  diffuse[0] = light_in.diffuse[0];
  diffuse[1] = light_in.diffuse[1];
  diffuse[2] = light_in.diffuse[2];
  emission[0] = light_in.emission[0];
  emission[1] = light_in.emission[1];
  emission[2] = light_in.emission[2];
  shininess = light_in.shininess;
}

Primitive::Primitive(){
  multitexture = 0;
  textured = 0;
  transparent = 0;
  alpha = 1.0;
  size = 1.0;
  colour[0] = colour[1] = colour[2] = 1.0;
  colour_override = 0;

  cardVertexBuffer = 0;
  cardNormalBuffer = 0;
  cardColorBuffer = 0;
  cardIndexBuffer = 0;
  cardOriginBuffer = 0;
  cardSizeBuffer = 0;
  bufferVertices = 0;
  bufferOrigins = 0;
  bufferSizes = 0;
  bufferNormals = 0;
  bufferColors = 0;
  bufferIndices = 0;
  arraysGenerated = 0;
  arraysUploaded = 0;
  nVerts=0;
  nRealVerts=0;
  useVBO = 0;
  useVertexArrays = 0;
  useInterleaved = 0;
  allBuffer = 0;
  bufferVNC = 0;

}

Point::Point() : Primitive(){vertices.push_back(Cartesian());}
Line::Line() : Primitive(){}
LineStrip::LineStrip() : Primitive(){}
DashLine::DashLine() : Primitive(){}
CircleElement::CircleElement() : Primitive(){vertices.push_back(Cartesian());}
PartialCircleElement::PartialCircleElement() : CircleElement(){}
TriangleElement::TriangleElement() : Primitive(){}
PolygonElement::PolygonElement() : Primitive(){}
QuadElement::QuadElement() : Primitive(){}
SphereElement::SphereElement() : Primitive(){vertices.push_back(Cartesian());accu = 0;}
CylinderElement::CylinderElement() : Primitive(){}
QuadStripElement::QuadStripElement() : Primitive(){}
TriangleFanElement::TriangleFanElement() : Primitive(){}
TriangleStripElement::TriangleStripElement() : Primitive(){}
PolyCylinder::PolyCylinder() : Primitive(){}
Ribbon::Ribbon() : Primitive(){polygonType = GL_QUADS;}
Worm::Worm() : Ribbon(){}
VariableWorm::VariableWorm() : Ribbon(){}
ArrowHeadRibbon::ArrowHeadRibbon() : Ribbon(){}
BillBoard::BillBoard() : Primitive(){}
SimpleText::SimpleText() : Primitive(){
  text_height=0;
  text_width=0;
  centered = false;
}
Text::Text() : SimpleText(){}
BillBoardText::BillBoardText() : SimpleText(){}

void Primitive::set_origin(const double *origin_in){
  set_origin(Cartesian(origin_in));
}
void Primitive::set_origin(const Cartesian &origin_in){
  origin = origin_in; 
  vertices[0] = origin_in;
}

Primitive::~Primitive(){
  freeResources();
}

SimpleText::~SimpleText(){
}

void Primitive::set_draw_colour_poly(void){
  glColor4d(colour[0],colour[1],colour[2], alpha);
  
}

void Primitive::set_draw_colour_line(void){
  glColor4d(colour[0],colour[1],colour[2],1.0);
}

void Primitive::set_draw_colour_poly_override(const GLfloat *col){

  glColor4d(col[0],col[1],col[2], alpha );

}

void Primitive::set_draw_colour_line_override(const GLfloat *col){
	/* This method seems not to be used anymore. But I have made potentially
	 * expensive fix for const correctness. SJM 16/6/2009
	 */
  GLfloat coltmp[4] = {col[0],col[1],col[2],GetAlpha()};
  glColor4fv(coltmp);
}

void Primitive::initialize_lighting(void){
}

void Primitive::apply_textures(void){
}

void Primitive::unbind_textures(void){
}

PointElement::PointElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in) : Point(vertex,colour_in,origin_in,size_in,alpha_in){
}

Point::Point(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in) : Primitive(){
  size = size_in;
  alpha = alpha_in;
  vertices.push_back(vertex);
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
}

Point::~Point(){
}

Line::Line(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in) : Primitive(){

  size = size_in;
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;

}

Line::~Line(){
}

LineStrip::LineStrip(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in) : Primitive(){
  size = size_in;
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;
}

LineStrip::~LineStrip(){
}

DashLine::DashLine(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in) : Primitive(){
  size = size_in;
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;
  dash_length = 0.2;
}

DashLine::~DashLine(){
}

Circle::Circle(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in) : CircleElement(vertex,colour_in,origin_in,size_in,alpha_in,textured_in){

}
PartialCircle::PartialCircle(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in, double startAngle_in, double sweepAngle_in) : PartialCircleElement(vertex,colour_in,origin_in,size_in,alpha_in,textured_in,startAngle_in,sweepAngle_in){
}

PartialCircleElement::PartialCircleElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in, double startAngle_in, double sweepAngle_in) : CircleElement(vertex,colour_in,origin_in,size_in,alpha_in,textured_in){
	startAngle = startAngle_in;
	sweepAngle = sweepAngle_in;
}

CircleElement::CircleElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in) : Primitive(){
  size = size_in;
  alpha = alpha_in;
  vertices.clear();
  vertices.push_back(vertex);
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  textured = textured_in;
}

CircleElement::~CircleElement(){
}

Circle::~Circle(){
}

PartialCircleElement::~PartialCircleElement(){
}

PartialCircle::~PartialCircle(){
}

TriangleElement::TriangleElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : Primitive(){
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;
  textured = textured_in;
  Cartesian normal = Cartesian::CrossProduct(vertices[1]-vertices[0],vertices[2]-vertices[0]);
  normals.push_back(normal);
  normals.push_back(normal);
  normals.push_back(normal);
}

PolygonElement::PolygonElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : Primitive(){
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;
  textured = textured_in;
}

TriangleElement::~TriangleElement(){
}

PolygonElement::~PolygonElement(){
}

TriangleFanElement::TriangleFanElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : Primitive(){
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;
  textured = textured_in;
}

TriangleFanElement::~TriangleFanElement(){
}

TriangleStripElement::TriangleStripElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : Primitive(){
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;
  textured = textured_in;
}

TriangleStripElement::~TriangleStripElement(){
}

QuadElement::QuadElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : Primitive(){
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;
  textured = textured_in;
  Cartesian normal = Cartesian::CrossProduct(vertices[1]-vertices[0],vertices[2]-vertices[0]);
  normals.push_back(normal);
  normals.push_back(normal);
  normals.push_back(normal);
  normals.push_back(normal);
}

QuadElement::~QuadElement(){
}

QuadStripElement::QuadStripElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : Primitive(){
  alpha = alpha_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;
  textured = textured_in;
}

QuadStripElement::~QuadStripElement(){
}

DashCylinderElement::~DashCylinderElement(){
}

DashCylinderElement::DashCylinderElement() : CylinderElement(){
}

DashCylinderElement::DashCylinderElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int accu_in, int textured_in) : CylinderElement(vertices_in, colour_in, origin_in, size_in, alpha_in, accu_in,textured_in){
  textured = textured_in;
  dash_length = 0.5;
  dash_end = 1;
}

ConeElement::~ConeElement(){
}

ConeElement::ConeElement() : CylinderElement(){
}

ConeElement::ConeElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int accu_in, int textured_in) : CylinderElement(vertices_in, colour_in, origin_in, size_in, alpha_in, accu_in,textured_in){
  textured = textured_in;
}


void ConeElement::draw(const double *override_colour, int selective_override){
  draw_cone(vertices[0],vertices[1],size,accu,textured);
}

CylinderElement::CylinderElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int accu_in, int textured_in) : Primitive(){
  size = size_in;
  alpha = alpha_in;
  accu = accu_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  
  vertices = vertices_in;

  textured = textured_in;
}

CylinderElement::~CylinderElement(){
}

PolyCylinder::PolyCylinder(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double size_in, double alpha_in, int accu_in, int textured_in) : Primitive(){
  size = size_in;
  alpha = alpha_in;
  accu = accu_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;

  vertices = vertices_in;
  colour_vector = colour_vector_in;
  
  textured = textured_in;
}

PolyCylinder::~PolyCylinder(){
}

SpheroidElement::SpheroidElement() : Primitive(){
  a = b = c = 1.0;
  show_axes = 0;
  show_solid = 1;
}

Spheroid::Spheroid(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, const matrix &mat_in, double alpha_in, int accu_in, int show_axes_in, int show_solid_in, int textured_in) : SpheroidElement(vertex, colour_in, origin_in, mat_in, alpha_in, accu_in, show_axes_in, show_solid_in, textured_in){
}

SpheroidElement::SpheroidElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, const matrix &mat_in, double alpha_in, int accu_in, int show_axes_in, int show_solid_in, int textured_in) : Primitive(){
  accu = accu_in;
  vertices.clear();
  vertices.push_back(vertex);
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  mat = mat_in;
  show_axes = show_axes_in;
  show_solid = show_solid_in;
  textured = textured_in;
}

SpheroidElement::SpheroidElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double a_in, double b_in, double c_in, const Cartesian& xaxis_in, const Cartesian& yaxis_in, const Cartesian& zaxis_in, double alpha_in, int accu_in, int show_axes_in, int show_solid_in, int textured_in) : Primitive(){
  a = a_in;
  b = b_in;
  c = c_in;
  accu = accu_in;
  vertices.clear();
  vertices.push_back(vertex);
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  xaxis = xaxis_in;
  yaxis = yaxis_in;
  zaxis = zaxis_in;
  xaxis.normalize();
  yaxis.normalize();
  zaxis.normalize();
  mat = matrix();
  /*
    std::cout << "Normalized: " << xaxis << "\n";
    std::cout << "Normalized: " << yaxis << "\n";
    std::cout << "Normalized: " << zaxis << "\n";
  */
  show_axes = show_axes_in;
  show_solid = show_solid_in;
  textured = textured_in;
}

void SpheroidElement::ShowSolid(){
  show_solid = 1;
}

void SpheroidElement::HideSolid(){
  show_solid = 0;
}
void SpheroidElement::ShowAxes(){
  show_axes = 1;
}

void SpheroidElement::HideAxes(){
  show_axes = 0;
}

void SpheroidElement::draw(const double *override_colour, int selective_override){

  if(mat.get_rows()==4&&mat.get_columns()==4){
    GLboolean normalize;
    glGetBooleanv(GL_NORMALIZE,&normalize);
    glEnable(GL_NORMALIZE);

    glPushMatrix();
    glTranslatef(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());

    double *mdp = mat.to_dp();

    float lx = (mat*Cartesian(1,0,0)).length();
    float ly = (mat*Cartesian(0,1,0)).length();
    float lz = (mat*Cartesian(0,0,1)).length();

    glMultMatrixd(mdp);

  float strip_offset = 1.001;
  float strip_width = .03;

  if(show_axes){
    //glDisable(GL_LIGHTING);
    if(show_solid) {
      GLfloat *params = new GLfloat[4];
      glGetFloatv(GL_COLOR_CLEAR_VALUE,params);
      GLfloat colour[4] = {1.0-params[0], 1.0-params[1], 1.0-params[2], 1.0};
      set_line_colour(colour[0],colour[1],colour[2],colour[3]);
    } else {
      set_line_colour(colour[0],colour[1],colour[2],1.0);
    }
    glBegin(GL_TRIANGLE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = strip_offset*cos(theta);
      double y = strip_offset*sin(theta);
      glNormal3f(x,y,0);
      glVertex3f(x,y,strip_width/lz);
      glVertex3f(x,y,-strip_width/lz);
    }    
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = strip_offset*cos(theta);
      double y = strip_offset*sin(theta);
      glNormal3f(-x,-y,0);
      glVertex3f(x,y,-strip_width/lz);
      glVertex3f(x,y,strip_width/lz);
    }    
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = strip_offset*cos(theta);
      double y = strip_offset*sin(theta);
      glNormal3f(0,x,y);
      glVertex3f(strip_width/lx,x,y);
      glVertex3f(-strip_width/lx,x,y);
    }    
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = strip_offset*cos(theta);
      double y = strip_offset*sin(theta);
      glNormal3f(0,-x,-y);
      glVertex3f(-strip_width/lx,x,y);
      glVertex3f(strip_width/lx,x,y);
    }    
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = strip_offset*cos(theta);
      double y = strip_offset*sin(theta);
      glNormal3f(x,0,y);
      glVertex3f(x,-strip_width/ly,y);
      glVertex3f(x,strip_width/ly,y);
    }    
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = strip_offset*cos(theta);
      double y = strip_offset*sin(theta);
      glNormal3f(-x,0,-y);
      glVertex3f(x,strip_width/ly,y);
      glVertex3f(x,-strip_width/ly,y);
    }    
    glEnd();
    set_draw_colour();
    glPolygonOffset(1.0,1.0);
  }

    if(show_solid){
      glEnable(GL_LIGHTING);
      int my_accu;
      if(RenderQuality::GetRenderQuality())
        my_accu = 4;
      else
        my_accu = accu;
        sphere(my_accu,1.0);
    }
    glPopMatrix();
    return;
  }

  Cartesian cartx = Cartesian(1,0,0);
  Cartesian carty = Cartesian(0,1,0);
  Cartesian yprim;
  Quat q1,q2;
  matrix m(4,kdelta);

  double angle = acos(cartx.DotProduct(cartx,xaxis))*180.0/M_PI;
  //std::cout << "Angle 1: " << angle << "\n";
  if(fabs(angle)>0.00001){
    Cartesian rotax = cartx.CrossProduct(cartx,xaxis);
    //std::cout << "rotax: " << rotax;
    q1 = Quat(rotax,1,-angle);
  }
  m = q1.getInvMatrix();

  yprim = m * carty;
  angle = acos(yprim.DotProduct(yprim,yaxis))*180.0/M_PI;
  //std::cout << "Angle 2: " << angle << "\n";

  if(fabs(angle)>0.00001){
    q2 = Quat(cartx,1,-angle);
    q1.postMult(q2);
  }
  m = q1.getMatrix();
  double *mdp = m.to_dp();


  GLboolean normalize;
  glGetBooleanv(GL_NORMALIZE,&normalize);
  glEnable(GL_NORMALIZE);

  glPushMatrix();
  glTranslatef(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());

  glMultMatrixd(mdp);

  glScaled(a,b,c);

  if(show_axes){
    glDisable(GL_LIGHTING);
    if(show_solid) {
      GLfloat *params = new GLfloat[4];
      glGetFloatv(GL_COLOR_CLEAR_VALUE,params);
      GLfloat colour[4] = {1.0-params[0], 1.0-params[1], 1.0-params[2], 1.0};
      set_line_colour(colour[0],colour[1],colour[2],colour[3]);
    } else {
      set_line_colour(colour[0],colour[1],colour[2],1.0);
    }
    glBegin(GL_LINE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = cos(theta);
      double y = sin(theta);
      glVertex3f(x,y,0);
    }    
    glEnd();
    glBegin(GL_LINE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = cos(theta);
      double y = sin(theta);
      glVertex3f(0,x,y);
    }    
    glEnd();
    glBegin(GL_LINE_STRIP);
    for(int i=0;i<=360;i+=10){
      double theta = (double)i/360.0 * PIBY2;
      double x = cos(theta);
      double y = sin(theta);
      glVertex3f(x,0,y);
    }    
    glEnd();
    set_draw_colour();
    glPolygonOffset(1.0,1.0);
  }

  if(show_solid){
    glEnable(GL_LIGHTING);
    int my_accu;
    if(RenderQuality::GetRenderQuality())
      my_accu = 4;
    else
      my_accu = accu;
      sphere(my_accu,1.0);
  }

  glPopMatrix();

  if(!normalize)
    glDisable(GL_NORMALIZE);

  delete [] mdp;
}

SpheroidElement::~SpheroidElement(){
}

SphereElement::SphereElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int accu_in, int textured_in) : Primitive(){
  size = size_in;
  alpha = alpha_in;
  accu = accu_in;
  vertices.clear();
  vertices.push_back(vertex);
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  textured = textured_in;
}

SphereElement::~SphereElement(){
}

ArrowHeadRibbon::~ArrowHeadRibbon(){
}

Ribbon::~Ribbon(){
}

Worm::~Worm(){
}

VariableWorm::~VariableWorm(){
}

void BillBoard::SetScale(double scale_w_in, double scale_h_in){
  scale_w = scale_w_in;
  scale_h = scale_h_in;
}

void BillBoard::set_draw_colour(const GLfloat *col){
}

void Point::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_line_override(col);
  else
    set_draw_colour_line();
}

void Line::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_line_override(col);
  else
    set_draw_colour_line();
}

void LineStrip::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_line_override(col);
  else
    set_draw_colour_line();
}

void DashLine::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_line_override(col);
  else
    set_draw_colour_line();
}

void CircleElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void TriangleElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void PolygonElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void QuadElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void QuadStripElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void TriangleStripElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void TriangleFanElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void CylinderElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void PolyCylinder::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void SphereElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void SpheroidElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void ArrowHeadRibbon::set_draw_colour(const GLfloat *col){
  Ribbon::set_draw_colour(col);
}

void Ribbon::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void Point::draw(const double *override_colour, int selective_override){
  glPointSize(size);
  glBegin(GL_POINTS);
  glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  glEnd();
}


DashArrowElement::~DashArrowElement(){
}

DashArrowElement::DashArrowElement() : LineElement(){
}

DashArrowElement::DashArrowElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in,int arrow_head_in) : LineElement(vertices_in, colour_in, origin_in, size_in, alpha_in){
  arrow_head = arrow_head_in;
  dash_length = 0.2;
  arrowLength = 0.1 * (vertices[1]-vertices[0]).length();
  arrowWidth = 0.5 * arrowLength;
}

DashArrow::~DashArrow(){
}

DashArrow::DashArrow() : DashArrowElement(){
}

DashArrow::DashArrow(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in,int arrow_head_in) : DashArrowElement(vertices_in, colour_in, origin_in, size_in, alpha_in){
  arrow_head = arrow_head_in;
  dash_length = 0.2;
  
}



ArrowElement::~ArrowElement(){
}

ArrowElement::ArrowElement() : LineElement(){
}

ArrowElement::ArrowElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in,int arrow_head_in) : LineElement(vertices_in, colour_in, origin_in, size_in, alpha_in){
  arrow_head = arrow_head_in;
  arrowLength = 0.1 * (vertices[1]-vertices[0]).length();
  arrowWidth = 0.5 * arrowLength;
}

Arrow::~Arrow(){
}

Arrow::Arrow() : ArrowElement(){
}

Arrow::Arrow(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in,int arrow_head_in) : ArrowElement(vertices_in, colour_in, origin_in, size_in, alpha_in){
  arrow_head = arrow_head_in;
}

void Arrow::draw(const double *override_colour, int selective_override){
  glLineWidth(size);
  glBegin(GL_LINES);
  ArrowElement::draw(override_colour,selective_override);
  glEnd();
}


void ArrowElement::draw(const double *override_colour, int selective_override){

  Cartesian vec = vertices[1] - vertices[0];

  Cartesian atmp,btmp;
  Cartesian xaxis(1,0,0);
  Cartesian yaxis(0,1,0);
  Cartesian zaxis(0,0,1);
  atmp = vec.CrossProduct(vec,zaxis);
  if(atmp.length() < 0.000000001){
    atmp = vec.CrossProduct(vec,yaxis);
    if(atmp.length() < 0.000000001){
      atmp = vec.CrossProduct(vec,xaxis);
    }
  }
  btmp = vec.CrossProduct(vec,atmp);

  double frac = arrowLength/vec.length();
  Cartesian p = frac*vertices[0]+(1-frac)*vertices[1];
  Cartesian pp = frac*vertices[1]+(1-frac)*vertices[0];

  btmp.normalize(arrowWidth);
  Cartesian p1 = p + btmp;
  Cartesian p2 = p - btmp;
  atmp.normalize(arrowWidth);
  Cartesian p3 = p + atmp;
  Cartesian p4 = p - atmp;
  
  glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
  if (arrow_head != 1) {
    //   arrow_head = 0 (end) or 2 (both) 
     glVertex3f(p1.get_x(),p1.get_y(),p1.get_z());
     glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
     glVertex3f(p2.get_x(),p2.get_y(),p2.get_z());
     glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
     glVertex3f(p3.get_x(),p3.get_y(),p3.get_z());
     glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
     glVertex3f(p4.get_x(),p4.get_y(),p4.get_z());
     glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
  }
  if (arrow_head != 0) {
    //   arrow_head = 1 (start) or 2 (both)
     p1= pp + btmp;
     p2= pp - btmp;
     p3= pp + atmp;
     p4= pp - atmp;
     glVertex3f(p1.get_x(),p1.get_y(),p1.get_z());
     glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
     glVertex3f(p2.get_x(),p2.get_y(),p2.get_z());
     glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
     glVertex3f(p3.get_x(),p3.get_y(),p3.get_z());
     glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
     glVertex3f(p4.get_x(),p4.get_y(),p4.get_z());
     glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  }

}

void DashArrow::draw(const double *override_colour, int selective_override){
  glLineWidth(size);
  glBegin(GL_LINES);
  DashArrowElement::draw(override_colour,selective_override);
  glEnd();
}

void DashArrowElement::draw(const double *override_colour, int selective_override){
  Cartesian vec = vertices[1] - vertices[0];
  Cartesian atmp,btmp;
  Cartesian xaxis(1,0,0);
  Cartesian yaxis(0,1,0);
  Cartesian zaxis(0,0,1);
  atmp = vec.CrossProduct(vec,zaxis);
  if(atmp.length() < 0.000000001){
    atmp = vec.CrossProduct(vec,yaxis);
    if(atmp.length() < 0.000000001){
      atmp = vec.CrossProduct(vec,xaxis);
    }
  }
  btmp = vec.CrossProduct(vec,atmp);

  double frac = arrowLength/vec.length();
  Cartesian p = frac*vertices[0]+(1-frac)*vertices[1];
  Cartesian pp = frac*vertices[1]+(1-frac)*vertices[0];

  btmp.normalize(arrowWidth);
  Cartesian p1 = p + btmp;
  Cartesian p2 = p - btmp;
  atmp.normalize(arrowWidth);
  Cartesian p3 = p + atmp;
  Cartesian p4 = p - atmp;

  if ( arrow_head != 1) {
     glVertex3f(p1.get_x(),p1.get_y(),p1.get_z());
     glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
     glVertex3f(p2.get_x(),p2.get_y(),p2.get_z());
     glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
     glVertex3f(p3.get_x(),p3.get_y(),p3.get_z());
     glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
     glVertex3f(p4.get_x(),p4.get_y(),p4.get_z());
     glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
  }
  if ( arrow_head != 0) {
     p1= pp + btmp;
     p2= pp - btmp;
     p3= pp + atmp;
     p4= pp - atmp;
     glVertex3f(p1.get_x(),p1.get_y(),p1.get_z());
     glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
     glVertex3f(p2.get_x(),p2.get_y(),p2.get_z());
     glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
     glVertex3f(p3.get_x(),p3.get_y(),p3.get_z());
     glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
     glVertex3f(p4.get_x(),p4.get_y(),p4.get_z());
     glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  }


  double length = sqrt((vertices[1].get_x() - vertices[0].get_x()) * (vertices[1].get_x() - vertices[0].get_x()) +
		       (vertices[1].get_y() - vertices[0].get_y()) * (vertices[1].get_y() - vertices[0].get_y()) +
		       (vertices[1].get_z() - vertices[0].get_z()) * (vertices[1].get_z() - vertices[0].get_z()));
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2) == 0) nsegments++;

  for(int i=0;i<nsegments;i++){
    if(i%2 >0){
      double frac = (double)i/(double)nsegments;
      double frac2 =(double)(i+1)/(double)nsegments;
      glVertex3f(frac*vertices[0].get_x() + (1-frac)*vertices[1].get_x(),
		 frac*vertices[0].get_y() + (1-frac)*vertices[1].get_y(),
		 frac*vertices[0].get_z() + (1-frac)*vertices[1].get_z());
      glVertex3f(frac2*vertices[0].get_x() + (1-frac2)*vertices[1].get_x(),
		 frac2*vertices[0].get_y() + (1-frac2)*vertices[1].get_y(),
		 frac2*vertices[0].get_z() + (1-frac2)*vertices[1].get_z());
    }
  }
}

Triangle::Triangle() : TriangleElement(){}
Circle::Circle() : CircleElement(){}
PartialCircle::PartialCircle() : PartialCircleElement(){}
QuadStrip::QuadStrip() : QuadStripElement(){}
TriangleFan::TriangleFan() : TriangleFanElement(){}
TriangleStrip::TriangleStrip() : TriangleStripElement(){}
MGPolygon::MGPolygon() : PolygonElement(){}
Quad::Quad() : QuadElement(){}
Sphere::Sphere() : SphereElement(){}
Spheroid::Spheroid() : SpheroidElement(){}
Cylinder::Cylinder() : CylinderElement(){}
Cone::Cone() : ConeElement(){}
DashCylinder::DashCylinder() : DashCylinderElement(){}

TriangleStrip::TriangleStrip(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : TriangleStripElement(vertices_in, colour_in, origin_in, alpha_in, textured_in){};

TriangleFan::TriangleFan(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : TriangleFanElement(vertices_in, colour_in, origin_in, alpha_in, textured_in){};

QuadStrip::QuadStrip(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : QuadStripElement(vertices_in, colour_in, origin_in, alpha_in, textured_in){};

Triangle::Triangle(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : TriangleElement(vertices_in, colour_in, origin_in, alpha_in, textured_in){};

MGPolygon::MGPolygon(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : PolygonElement(vertices_in, colour_in, origin_in, alpha_in, textured_in){};

Quad::Quad(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in,  double alpha_in, int textured_in) : QuadElement(vertices_in, colour_in, origin_in, alpha_in, textured_in){};

Sphere::Sphere(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int accu_in, int textured_in) : SphereElement(vertex, colour_in, origin_in, size_in, alpha_in, accu_in, textured_in){};

Spheroid::Spheroid(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double a_in, double b_in, double c_in, const Cartesian& xaxis_in, const Cartesian& yaxis_in, const Cartesian& zaxis_in, double alpha_in, int accu_in, int show_axes, int show_solid, int textured_in) : SpheroidElement(vertex, colour_in, origin_in, a_in, b_in, c_in, xaxis_in, yaxis_in, zaxis_in, alpha_in, accu_in, show_axes, show_solid, textured_in){};

Cylinder::Cylinder(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int accu_in, int textured_in) : CylinderElement(vertices_in, colour_in, origin_in, size_in, alpha_in, accu_in, textured_in){};

Cone::Cone(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int accu_in, int textured_in) : ConeElement(vertices_in, colour_in, origin_in, size_in, alpha_in, accu_in, textured_in){};

DashCylinder::DashCylinder(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int accu_in, int textured_in) : DashCylinderElement(vertices_in, colour_in, origin_in, size_in, alpha_in, accu_in, textured_in){};

void TriangleStrip::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  TriangleStripElement::draw(override_colour,selective_override);
}
void TriangleFan::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  TriangleFanElement::draw(override_colour,selective_override);
}
void QuadStrip::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  QuadStripElement::draw(override_colour,selective_override);
}
void Triangle::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  TriangleElement::draw(override_colour,selective_override);
}
void MGPolygon::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  PolygonElement::draw(override_colour,selective_override);
}

void Quad::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  QuadElement::draw(override_colour,selective_override);
}

void Sphere::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  SphereElement::draw(override_colour,selective_override);
}

void Spheroid::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  SpheroidElement::draw(override_colour,selective_override);
}

void Cylinder::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  CylinderElement::draw(override_colour,selective_override);
}

void Cone::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  ConeElement::draw(override_colour,selective_override);
}

LineElement::LineElement() : Line(){}
PointElement::PointElement() : Point(){}
DashLineElement::DashLineElement() : LineElement(){}
LineElement::~LineElement(){}
PointElement::~PointElement(){}
DashLineElement::~DashLineElement(){}
LineElement::LineElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in) : Line(vertices_in,colour_in,origin_in,size_in){}
DashLineElement::DashLineElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in) : LineElement(vertices_in,colour_in,origin_in,size_in){
  dash_length = 0.2;
}

std::vector<Primitive*> LineCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  std::vector<Primitive*> a;
  nVerts = 0;
  nRealVerts = 0;
  if(lines.size()<1) return a;

  bufferVertices = new float[6*lines.size()];
  bufferNormals = 0;
  bufferColors = new float[8*lines.size()];
  bufferIndices = new GLuint[2*lines.size()];
  nRealVerts = 6*lines.size();
  nVerts = 2*lines.size();

  for(unsigned i=0;i<lines.size();i++){
    bufferIndices[2*i] = 2*i;
    bufferIndices[2*i+1] = 2*i+1;
    bufferVertices[6*i]   = lines[i]->GetVertices()[0].get_x();
    bufferVertices[6*i+1] = lines[i]->GetVertices()[0].get_y();
    bufferVertices[6*i+2] = lines[i]->GetVertices()[0].get_z();
    bufferVertices[6*i+3] = lines[i]->GetVertices()[1].get_x();
    bufferVertices[6*i+4] = lines[i]->GetVertices()[1].get_y();
    bufferVertices[6*i+5] = lines[i]->GetVertices()[1].get_z();
    bufferColors[8*i]   = lines[i]->GetColour()[0];
    bufferColors[8*i+1] = lines[i]->GetColour()[1];
    bufferColors[8*i+2] = lines[i]->GetColour()[2];
    bufferColors[8*i+3] = 1;
    bufferColors[8*i+4] = lines[i]->GetColour()[0];
    bufferColors[8*i+5] = lines[i]->GetColour()[1];
    bufferColors[8*i+6] = lines[i]->GetColour()[2];
    bufferColors[8*i+7] = 1;
  }

  polygonType = GL_LINES;
  return a;
}

std::vector<Primitive*> DashLinesCollection::GetPrimitives() const {
  std::vector<Primitive*> a;
  std::vector<Cartesian> carts(2);

  if(line_cartesians.size()<6) return a;
  double length = sqrt((line_cartesians[3]-line_cartesians[0]) * (line_cartesians[3]-line_cartesians[0]) +
		       (line_cartesians[4]-line_cartesians[1]) * (line_cartesians[4]-line_cartesians[1]) +
		       (line_cartesians[5]-line_cartesians[2]) * (line_cartesians[5]-line_cartesians[2]));
  
  int nsegments = (int)(length/dash_length);
  if((nsegments%2) == 1) nsegments++;

  int iiprim=0;
  for(unsigned i=0;i<line_cartesians.size();i+=6,iiprim+=2){
    double frac = 0;
    double frac2 =0.5/(double)nsegments;
    double col[4] = {line_colours[iiprim*4],line_colours[iiprim*4+1],line_colours[iiprim*4+2],1};
    Cartesian c1(frac*line_cartesians[i+3] + (1-frac)*line_cartesians[i],frac*line_cartesians[i+4] + (1-frac)*line_cartesians[i+1],frac*line_cartesians[i+5] + (1-frac)*line_cartesians[i+2]);
    Cartesian c2(frac2*line_cartesians[i+3] + (1-frac2)*line_cartesians[i],frac2*line_cartesians[i+4] + (1-frac2)*line_cartesians[i+1],frac2*line_cartesians[i+5] + (1-frac2)*line_cartesians[i+2]);
    carts[0] = c1;
    carts[1] = c2;
    LineElement *line = new LineElement(carts,col,carts[1],size,1.0);
    a.push_back(line);

    for(int j=0;j<nsegments;j++){
      if(j%2 >0){
        double frac = (j+0.5)/(double)nsegments;
        double frac2 =(j+1.5)/(double)nsegments;
        if(frac2>1.0) frac2 = 1.0;
        double col[4] = {line_colours[iiprim*4],line_colours[iiprim*4+1],line_colours[iiprim*4+2],1};
        Cartesian c1(frac*line_cartesians[i+3] + (1-frac)*line_cartesians[i],frac*line_cartesians[i+4] + (1-frac)*line_cartesians[i+1],frac*line_cartesians[i+5] + (1-frac)*line_cartesians[i+2]);
        Cartesian c2(frac2*line_cartesians[i+3] + (1-frac2)*line_cartesians[i],frac2*line_cartesians[i+4] + (1-frac2)*line_cartesians[i+1],frac2*line_cartesians[i+5] + (1-frac2)*line_cartesians[i+2]);
        carts[0] = c1;
        carts[1] = c2;
        LineElement *line = new LineElement(carts,col,carts[1],size,1.0);
        a.push_back(line);
      }
    }
  }

  return a;

}

std::vector<Primitive*> LinesCollection::GetPrimitives() const {
  std::vector<Primitive*> a;
  std::vector<Cartesian> carts(2);
  int iprim=0;
  for(unsigned i=0;i<line_cartesians.size();i+=6,iprim+=2){
    double col[4] = {line_colours[iprim*4],line_colours[iprim*4+1],line_colours[iprim*4+2],1};
    Cartesian c1(line_cartesians[i],line_cartesians[i+1],line_cartesians[i+2]);
    Cartesian c2(line_cartesians[i+3],line_cartesians[i+4],line_cartesians[i+5]);
    carts[0] = c1;
    carts[1] = c2;
    LineElement *line = new LineElement(carts,col,carts[1],size,1.0);
    a.push_back(line);
  }
  return a;
}

std::vector<Primitive*> DashLinesCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  std::vector<Primitive*> a;
  nVerts = 0;
  nRealVerts = 0;

  if(line_cartesians.size()<6) return a;
  double length = sqrt((line_cartesians[3]-line_cartesians[0]) * (line_cartesians[3]-line_cartesians[0]) +
		       (line_cartesians[4]-line_cartesians[1]) * (line_cartesians[4]-line_cartesians[1]) +
		       (line_cartesians[5]-line_cartesians[2]) * (line_cartesians[5]-line_cartesians[2]));
  
  int nsegments = (int)(length/dash_length);
  if((nsegments%2) == 1) nsegments++;

  bufferVertices = new float[line_cartesians.size()*(nsegments+2)/2];
  bufferNormals = 0;
  bufferColors = new float[line_colours.size()*(nsegments+2)/2];
  bufferIndices = new GLuint[line_cartesians.size()/3*(nsegments+2)/2];
  nRealVerts = line_cartesians.size()*nsegments/2;
  nVerts = line_cartesians.size()/3*nsegments/2;

  int maxv = line_cartesians.size();

  nRealVerts = maxv*(nsegments+2)/2;
  nVerts = maxv/3*(nsegments+2)/2;

  int iprim=0;
  int iiprim=0;
  for(unsigned i=0;i<maxv;i+=6,iiprim+=2){
    double frac = 0;
    double frac2 =0.5/(double)nsegments;
    bufferIndices[iprim] = iprim;
    bufferVertices[3*iprim]   = frac*line_cartesians[i+3] + (1-frac)*line_cartesians[i];
    bufferVertices[3*iprim+1] = frac*line_cartesians[i+4] + (1-frac)*line_cartesians[i+1];
    bufferVertices[3*iprim+2] = frac*line_cartesians[i+5] + (1-frac)*line_cartesians[i+2];
    bufferColors[4*iprim]   = line_colours[iiprim*4];
    bufferColors[4*iprim+1] = line_colours[iiprim*4+1];
    bufferColors[4*iprim+2] = line_colours[iiprim*4+2];
    bufferColors[4*iprim+3] = 1;
    bufferIndices[iprim+1] = iprim+1;
    bufferVertices[3*(iprim+1)]   = frac2*line_cartesians[i+3] + (1-frac2)*line_cartesians[i];
    bufferVertices[3*(iprim+1)+1] = frac2*line_cartesians[i+4] + (1-frac2)*line_cartesians[i+1];
    bufferVertices[3*(iprim+1)+2] = frac2*line_cartesians[i+5] + (1-frac2)*line_cartesians[i+2];
    bufferColors[4*(iprim+1)]   = line_colours[iiprim*4];
    bufferColors[4*(iprim+1)+1] = line_colours[iiprim*4+1];
    bufferColors[4*(iprim+1)+2] = line_colours[iiprim*4+2];
    bufferColors[4*(iprim+1)+3] = 1;
    iprim+=2;
    for(int j=0;j<nsegments;j++){
      if(j%2 >0){
        bufferIndices[iprim] = iprim;
        double frac = (j+0.5)/(double)nsegments;
        double frac2 =(j+1.5)/(double)nsegments;
        if(frac2>1.0) frac2 = 1.0;
        bufferVertices[3*iprim]   = frac*line_cartesians[i+3] + (1-frac)*line_cartesians[i];
        bufferVertices[3*iprim+1] = frac*line_cartesians[i+4] + (1-frac)*line_cartesians[i+1];
        bufferVertices[3*iprim+2] = frac*line_cartesians[i+5] + (1-frac)*line_cartesians[i+2];
        bufferColors[4*iprim]   = line_colours[iiprim*4];
        bufferColors[4*iprim+1] = line_colours[iiprim*4+1];
        bufferColors[4*iprim+2] = line_colours[iiprim*4+2];
        bufferColors[4*iprim+3] = 1;
        bufferIndices[iprim+1] = iprim+1;
        bufferVertices[3*(iprim+1)]   = frac2*line_cartesians[i+3] + (1-frac2)*line_cartesians[i];
        bufferVertices[3*(iprim+1)+1] = frac2*line_cartesians[i+4] + (1-frac2)*line_cartesians[i+1];
        bufferVertices[3*(iprim+1)+2] = frac2*line_cartesians[i+5] + (1-frac2)*line_cartesians[i+2];
        bufferColors[4*(iprim+1)]   = line_colours[iiprim*4];
        bufferColors[4*(iprim+1)+1] = line_colours[iiprim*4+1];
        bufferColors[4*(iprim+1)+2] = line_colours[iiprim*4+2];
        bufferColors[4*(iprim+1)+3] = 1;
        iprim+=2;
      }
    }
  }

  polygonType = GL_LINES;

  return a;

  return LinesCollection::GetSimplePrimitivesArrays(clip_vol, objrotmatrix, objorigin, start, end);
}

std::vector<Primitive*> LinesCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  std::vector<Primitive*> a;
  nVerts = 0;
  nRealVerts = 0;

  if(line_cartesians.size()<6) return a;

  bufferVertices = new float[line_cartesians.size()];
  bufferNormals = 0;
  bufferColors = new float[line_colours.size()];
  bufferIndices = new GLuint[line_cartesians.size()/3];
  nRealVerts = line_cartesians.size();
  nVerts = line_cartesians.size()/3;

  int iprim=0;
  for(unsigned i=0;i<line_cartesians.size();i+=3,iprim++){
    bufferIndices[iprim] = iprim;
    bufferVertices[3*iprim]   = line_cartesians[i];
    bufferVertices[3*iprim+1] = line_cartesians[i+1];
    bufferVertices[3*iprim+2] = line_cartesians[i+2];
    bufferColors[4*iprim]   = line_colours[iprim*4];
    bufferColors[4*iprim+1] = line_colours[iprim*4+1];
    bufferColors[4*iprim+2] = line_colours[iprim*4+2];
    bufferColors[4*iprim+3] = 1;
  }

  polygonType = GL_LINES;

  return a;
}

DashLinesCollection::~DashLinesCollection(){
}

LinesCollection::~LinesCollection(){
/*
  std::vector<Primitive*>::iterator prim_iter = lines.begin();
  while(prim_iter!=lines.end()){
    delete *prim_iter;
    prim_iter++;
  }
  lines.clear();
*/
}

LinesCollection::LinesCollection() : Primitive(){}

DashLinesCollection::DashLinesCollection() : LinesCollection(){}

void DashLinesCollection::reserve(size_t nbonds){
  if(line_cartesians.size()<6) return;
  double length = sqrt((line_cartesians[3]-line_cartesians[0]) * (line_cartesians[3]-line_cartesians[0]) +
		       (line_cartesians[4]-line_cartesians[1]) * (line_cartesians[4]-line_cartesians[1]) +
		       (line_cartesians[5]-line_cartesians[2]) * (line_cartesians[5]-line_cartesians[2]));
  int nsegments = (int)(length/dash_length);
  if((nsegments%2) == 1) nsegments++;
  LinesCollection::reserve(nbonds*nsegments);
}

void LinesCollection::reserve(size_t nbonds){
  line_cartesians.reserve(3*nbonds);
  line_colours.reserve(4*nbonds);
  line_colour_override.reserve(nbonds);
}

void LinesCollection::add(const std::vector<Cartesian> &carts,const double *colour_array,int colour_override_in){
  line_cartesians.push_back(carts[0].get_x());
  line_cartesians.push_back(carts[0].get_y());
  line_cartesians.push_back(carts[0].get_z());
  line_cartesians.push_back(carts[1].get_x());
  line_cartesians.push_back(carts[1].get_y());
  line_cartesians.push_back(carts[1].get_z());
  line_colours.push_back(colour_array[0]);
  line_colours.push_back(colour_array[1]);
  line_colours.push_back(colour_array[2]);
  line_colours.push_back(colour_array[3]);
  line_colours.push_back(colour_array[0]);
  line_colours.push_back(colour_array[1]);
  line_colours.push_back(colour_array[2]);
  line_colours.push_back(colour_array[3]);
  line_colour_override.push_back(colour_override_in);
}


void DashLinesCollection::draw(const double *override_colour, int selective_override){

  if(line_cartesians.size()<6) return;

  glLineWidth(size);
  glDisable(GL_LIGHTING);
  
  double length = sqrt((line_cartesians[3] - line_cartesians[0]) * (line_cartesians[3] - line_cartesians[0]) +
		       (line_cartesians[4] - line_cartesians[1]) * (line_cartesians[4] - line_cartesians[1]) +
		       (line_cartesians[5] - line_cartesians[2]) * (line_cartesians[5] - line_cartesians[2]));
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2) == 0) nsegments++;
  
  if (useVBO) {
    drawArrays(override_colour,selective_override);
    return;
  }
  glBegin(GL_LINES);
  int iprim=0;
  if(override_colour){
    glColor3f(override_colour[0],override_colour[1],override_colour[2]);
    for(unsigned i=0;i<line_cartesians.size()-6;i+=6,iprim++){
      double frac = 0;
      double frac2 =0.5/(double)nsegments;
      glVertex3f(frac*line_cartesians[i+3] + (1-frac)*line_cartesians[i],
		 frac*line_cartesians[i+4] + (1-frac)*line_cartesians[i+1],
		 frac*line_cartesians[i+5] + (1-frac)*line_cartesians[i+2]);
      glVertex3f(frac2*line_cartesians[i+3] + (1-frac2)*line_cartesians[i],
		 frac2*line_cartesians[i+4] + (1-frac2)*line_cartesians[i+1],
		 frac2*line_cartesians[i+5] + (1-frac2)*line_cartesians[i+2]);
      for(int iseg=1;iseg<nsegments;iseg++){
        if(iseg%2 >0){
          frac = (iseg+0.5)/(double)nsegments;
          frac2 =(iseg+1.5)/(double)nsegments;
          if(frac2>1.0) frac2 = 1.0;
            glVertex3f(frac*line_cartesians[i+3] + (1-frac)*line_cartesians[i],
		 frac*line_cartesians[i+4] + (1-frac)*line_cartesians[i+1],
		 frac*line_cartesians[i+5] + (1-frac)*line_cartesians[i+2]);
            glVertex3f(frac2*line_cartesians[i+3] + (1-frac2)*line_cartesians[i],
		 frac2*line_cartesians[i+4] + (1-frac2)*line_cartesians[i+1],
		 frac2*line_cartesians[i+5] + (1-frac2)*line_cartesians[i+2]);
          }
      }
    }
  } else {
    for(unsigned i=0;i<line_cartesians.size()-5;i+=6,iprim++){
      std::cout << i << "\n";
      double frac = 0;
      double frac2 =0.5/(double)nsegments;
      glColor3f(line_colours[iprim*4],line_colours[iprim*4+1],line_colours[iprim*4+2]);
      glVertex3f(frac*line_cartesians[i+3] + (1-frac)*line_cartesians[i],
		 frac*line_cartesians[i+4] + (1-frac)*line_cartesians[i+1],
		 frac*line_cartesians[i+5] + (1-frac)*line_cartesians[i+2]);
      glVertex3f(frac2*line_cartesians[i+3] + (1-frac2)*line_cartesians[i],
		 frac2*line_cartesians[i+4] + (1-frac2)*line_cartesians[i+1],
		 frac2*line_cartesians[i+5] + (1-frac2)*line_cartesians[i+2]);
      for(int iseg=1;iseg<nsegments;iseg++){
        if(iseg%2 >0){
          frac = (iseg+0.5)/(double)nsegments;
          frac2 =(iseg+1.5)/(double)nsegments;
          if(frac2>1.0) frac2 = 1.0;
            glVertex3f(frac*line_cartesians[i+3] + (1-frac)*line_cartesians[i],
		 frac*line_cartesians[i+4] + (1-frac)*line_cartesians[i+1],
		 frac*line_cartesians[i+5] + (1-frac)*line_cartesians[i+2]);
            glVertex3f(frac2*line_cartesians[i+3] + (1-frac2)*line_cartesians[i],
		 frac2*line_cartesians[i+4] + (1-frac2)*line_cartesians[i+1],
		 frac2*line_cartesians[i+5] + (1-frac2)*line_cartesians[i+2]);
          }
      }
    }
  }
  glEnd();
}

void LinesCollection::draw(const double *override_colour, int selective_override){
  glLineWidth(size);
  glDisable(GL_LIGHTING);
  
  if (useVBO) {
    drawArrays(override_colour,selective_override);
    return;
  }
  glBegin(GL_LINES);
  int iprim=0;
  if(override_colour){
    glColor3f(override_colour[0],override_colour[1],override_colour[2]);
    for(unsigned i=0;i<line_cartesians.size();i+=3,iprim++){
      glVertex3f(line_cartesians[i],line_cartesians[i+1],line_cartesians[i+2]);
    }
  } else {
    for(unsigned i=0;i<line_cartesians.size();i+=3,iprim++){
      glColor3f(line_colours[iprim*4],line_colours[iprim*4+1],line_colours[iprim*4+2]);
      glVertex3f(line_cartesians[i],line_cartesians[i+1],line_cartesians[i+2]);
    }
  }
  glEnd();
}

void LineCollection::draw(const double *override_colour, int selective_override){
  glLineWidth(size);
  glDisable(GL_LIGHTING);
  
  if (useVBO) {
    drawArrays(override_colour,selective_override);
    return;
  }

  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
 

  GLfloat ocol[4];
  if(override_colour){
  ocol[0] = override_colour[0];
  ocol[1] = override_colour[1];
  ocol[2] = override_colour[2];
  ocol[3] = override_colour[3];
  }

  bool overridden = false;
  std::vector<Primitive*>::const_iterator k = lines.begin();

  glBegin(GL_LINES);
  if(override_colour&&selective_override==0){
    (*k)->set_draw_colour_line_override(ocol);
    while(k!=lines.end()){
      (*k)->draw();
      k++;
    }  
  }

  if(override_colour&&selective_override==1){
    while(k!=lines.end()){
      if((*k)->GetColourOverride()){
        if(!overridden){
          (*k)->set_draw_colour_line_override(ocol);
          overridden = true;
        }
        (*k)->draw();
      }
      k++;
    }  
    k = lines.begin();
    while(k!=lines.end()){
      if(!(*k)->GetColourOverride()){
        (*k)->set_draw_colour_line();
        (*k)->draw();
      }
      k++;
    }  
  }

  if(!override_colour){
    while(k!=lines.end()){
      (*k)->set_draw_colour_line();
      (*k)->draw();
      k++;
    }  
  }

  glEnd();


  //glBlendFunc (GL_SRC_ALPHA, GL_ZERO);
  //glDisable(GL_BLEND);

}
void PointCollection::draw(const double *override_colour, int selective_override){
  glPointSize(size);
  glDisable(GL_LIGHTING);
  
  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
 

  GLfloat ocol[4];
  if(override_colour){
  ocol[0] = override_colour[0];
  ocol[1] = override_colour[1];
  ocol[2] = override_colour[2];
  ocol[3] = override_colour[3];
  }

  bool overridden = false;
  std::vector<Primitive*>::const_iterator k = lines.begin();

  glBegin(GL_POINTS);
  if(override_colour&&selective_override==0){
    (*k)->set_draw_colour_line_override(ocol);
    while(k!=lines.end()){
      (*k)->draw();
      k++;
    }  
  }

  if(override_colour&&selective_override==1){
    while(k!=lines.end()){
      if((*k)->GetColourOverride()){
        if(!overridden){
          (*k)->set_draw_colour_line_override(ocol);
          overridden = true;
        }
        (*k)->draw();
      }
      k++;
    }  
    k = lines.begin();
    while(k!=lines.end()){
      if(!(*k)->GetColourOverride()){
        (*k)->set_draw_colour_line();
        (*k)->draw();
      }
      k++;
    }  
  }

  if(!override_colour){
    while(k!=lines.end()){
      (*k)->set_draw_colour_line();
      (*k)->draw();
      k++;
    }  
  }

  glEnd();


  //glBlendFunc (GL_SRC_ALPHA, GL_ZERO);
  //glDisable(GL_BLEND);

}

void Primitive::SetUseVertexArrays(const int _useVertexArrays){
  useVertexArrays=_useVertexArrays;
} 

void PolyCollection::SetUseVertexArrays(const int _useVertexArrays){
  useVertexArrays = _useVertexArrays;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->SetUseVertexArrays(_useVertexArrays);
    k++;
  }
};

void PolyCollection::SetAlpha(double alpha_in){
  //std::cout << "PolyCollection::SetAlpha " << alpha_in << "\n";
  alpha = alpha_in;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->SetAlpha(alpha);
    k++;
  }
  arraysGenerated = 0;
};

void PolyCollection::set_transparent(int trans){
  transparent = trans;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->set_transparent(trans);
    k++;
  }
}

void PolyCollection::draw(const double *override_colour, int selective_override){
  //std::cout << "PolyCollection::draw\n";
  glEnable(GL_LIGHTING);

  if (useVertexArrays) {
    drawArrays(override_colour,selective_override);
    return;
  }

  std::vector<Primitive*>::const_iterator k = prims.begin();

  GLfloat col[4];
  GLfloat *colp=0;
  int kk=0;
  while(k!=prims.end()){
    if(override_colour){
      if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	col[0] = override_colour[0];
	col[1] = override_colour[1];
	col[2] = override_colour[2];
	col[3] = override_colour[3];
	colp = col;
      }
    }
    (*k)->set_draw_colour(colp);
    (*k)->draw();
    k++; kk++;
  }  
}

void ImposterSphereCollection::SetUseVertexArrays(const int _useVertexArrays){
  useVertexArrays = _useVertexArrays;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->SetUseVertexArrays(_useVertexArrays);
    k++;
  }
};

void ImposterSphereCollection::SetAlpha(double alpha_in){
  //std::cout << "ImposterSphereCollection::SetAlpha " << alpha_in << "\n";
  alpha = alpha_in;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->SetAlpha(alpha);
    k++;
  }
  arraysGenerated = 0;
};

void ImposterSphereCollection::set_transparent(int trans){
  transparent = trans;
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->set_transparent(trans);
    k++;
  }
}

void ImposterSphereCollection::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);

  if (useVertexArrays) {
    drawArrays(override_colour,selective_override);
    return;
  }

  std::cout << "Not using arrays\n";
  std::vector<Primitive*>::const_iterator k = prims.begin();

  GLfloat col[4];
  GLfloat *colp=0;
  int kk=0;
  double modelMatrix[16];
  glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
  matrix mat(4,4,modelMatrix);
  Cartesian x_trans = mat*Cartesian(1,0,0);
  Cartesian y_trans = mat*Cartesian(0,1,0);
  Cartesian z_trans = mat*Cartesian(0,0,1);
  x_trans.normalize();
  y_trans.normalize();
  z_trans.normalize();
  glBegin(GL_QUADS);
  while(k!=prims.end()){
    if(override_colour){
      if(selective_override==0||(selective_override==1&&(*k)->GetColourOverride())){
	col[0] = override_colour[0];
	col[1] = override_colour[1];
	col[2] = override_colour[2];
	col[3] = override_colour[3];
	colp = col;
      }
    }
    (*k)->set_draw_colour(colp);
    ((ImposterSphereElement*)(*k))->draw(1,x_trans,y_trans,z_trans);
    k++; kk++;
  }  
  glEnd();
}

void LineCollection::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
  std::vector<Primitive*>::const_iterator k = lines.begin();

  while(k!=lines.end()){
    (*k)->DrawPovray(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,v);
    k++;
  }  
}

void PolyCollection::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
  std::vector<Primitive*>::const_iterator k = prims.begin();

  while(k!=prims.end()){
    (*k)->DrawPovray(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,v);
    k++;
  }  
}

void ImposterSphereCollection::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
  std::vector<Primitive*>::const_iterator k = prims.begin();

  while(k!=prims.end()){
    (*k)->DrawPovray(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,v);
    k++;
  }  
}

void LineCollection::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){
  std::vector<Primitive*>::const_iterator k = lines.begin();

  while(k!=lines.end()){
    (*k)->DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
    k++;
  }  
}

void PolyCollection::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){
  std::vector<Primitive*>::const_iterator k = prims.begin();

  while(k!=prims.end()){
    (*k)->DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
    k++;
  }  
}

void ImposterSphereCollection::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){
  std::vector<Primitive*>::const_iterator k = prims.begin();

  while(k!=prims.end()){
    (*k)->DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
    k++;
  }  
}

LineCollection::LineCollection() : Primitive(){}
PointCollection::PointCollection() : Primitive(){}
LineCollection::~LineCollection(){
  std::vector<Primitive*>::iterator prim_iter = lines.begin();
  while(prim_iter!=lines.end()){
    delete *prim_iter;
    prim_iter++;
  }
  lines.clear();
}
PointCollection::~PointCollection(){
  std::vector<Primitive*>::iterator prim_iter = lines.begin();
  while(prim_iter!=lines.end()){
    delete *prim_iter;
    prim_iter++;
  }
  lines.clear();
}
PolyCollection::PolyCollection() : Primitive(){}
TriangleCollection::TriangleCollection() : PolyCollection(){}
SphereCollection::SphereCollection() : PolyCollection(){sphere_accu=1;}
SpheroidCollection::SpheroidCollection() : SphereCollection(){}
ImposterSphereCollection::ImposterSphereCollection() : Primitive(){}
ImposterSphereCollection::~ImposterSphereCollection(){ }
CylinderCollection::CylinderCollection() : PolyCollection(){cylinder_accu=1;}
CylinderCollection::~CylinderCollection(){ }
DashCylinderCollection::DashCylinderCollection() : CylinderCollection(){cylinder_accu=1;}
DashCylinderCollection::~DashCylinderCollection(){ }
SphereCollection::~SphereCollection(){ }
SpheroidCollection::~SpheroidCollection(){ }
TriangleCollection::~TriangleCollection(){ }
PolyCollection::~PolyCollection(){
  std::vector<Primitive*>::iterator prim_iter = prims.begin();
  while(prim_iter!=prims.end()){
    delete *prim_iter;
    prim_iter++;
  }
  prims.clear();
}

void PointCollection::set_draw_colour(const GLfloat *col){
}
void LineCollection::set_draw_colour(const GLfloat *col){
}
void PolyCollection::set_draw_colour(const GLfloat *col){
}
void ImposterSphereCollection::set_draw_colour(const GLfloat *col){
}

void LineCollection::add_primitive(Primitive *line){
  lines.push_back(line);
}
void PointCollection::add_primitive(Primitive *line){
  lines.push_back(line);
}

void PolyCollection::add_primitive(Primitive *line){
  prims.push_back(line);
}
void ImposterSphereCollection::add_primitive(Primitive *line){
  prims.push_back(line);
}

void Line::draw(const double *override_colour, int selective_override){

  glBegin(GL_LINES);
  glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());
  glEnd();

}

void PointElement::draw(const double *override_colour, int selective_override){

  glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());

}

void LineElement::draw(const double *override_colour, int selective_override){

  glVertex3f(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  glVertex3f(vertices[1].get_x(),vertices[1].get_y(),vertices[1].get_z());

}

void LineStrip::draw(const double *override_colour, int selective_override){
  std::vector<Cartesian>::const_iterator k = vertices.begin();


    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);


  glLineWidth(size);
  glBegin(GL_LINES);
  while(k!=vertices.end()){
    glVertex3f((*k).get_x(),(*k).get_y(),(*k).get_z());
    k++;
  }
  glEnd();


    glBlendFunc (GL_SRC_ALPHA, GL_ZERO);

}

void DashLine::draw(const double *override_colour, int selective_override){
  glDisable(GL_LIGHTING);
  glLineWidth(size);
  double length = sqrt((vertices[1].get_x() - vertices[0].get_x()) * (vertices[1].get_x() - vertices[0].get_x()) +
		       (vertices[1].get_y() - vertices[0].get_y()) * (vertices[1].get_y() - vertices[0].get_y()) +
		       (vertices[1].get_z() - vertices[0].get_z()) * (vertices[1].get_z() - vertices[0].get_z()));
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2) == 0) nsegments++;

  glBegin(GL_LINES);
  double frac = 0;
  double frac2 =0.5/(double)nsegments;
  glVertex3f(frac*vertices[1].get_x() + (1-frac)*vertices[0].get_x(),
		 frac*vertices[1].get_y() + (1-frac)*vertices[0].get_y(),
		 frac*vertices[1].get_z() + (1-frac)*vertices[0].get_z());
  glVertex3f(frac2*vertices[1].get_x() + (1-frac2)*vertices[0].get_x(),
		 frac2*vertices[1].get_y() + (1-frac2)*vertices[0].get_y(),
		 frac2*vertices[1].get_z() + (1-frac2)*vertices[0].get_z());
  for(int i=1;i<nsegments;i++){
    if(i%2 >0){
      frac = (i+0.5)/(double)nsegments;
      frac2 =(i+1.5)/(double)nsegments;
      if(frac2>1.0) frac2 = 1.0;
      glVertex3f(frac*vertices[1].get_x() + (1-frac)*vertices[0].get_x(),
		 frac*vertices[1].get_y() + (1-frac)*vertices[0].get_y(),
		 frac*vertices[1].get_z() + (1-frac)*vertices[0].get_z());
      glVertex3f(frac2*vertices[1].get_x() + (1-frac2)*vertices[0].get_x(),
		 frac2*vertices[1].get_y() + (1-frac2)*vertices[0].get_y(),
		 frac2*vertices[1].get_z() + (1-frac2)*vertices[0].get_z());
    }
  }
  glEnd();
}

std::vector<Primitive*> DashLineElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "DashLineElement::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  std::vector<Cartesian> carts(2);
  double length = sqrt((vertices[1].get_x() - vertices[0].get_x()) * (vertices[1].get_x() - vertices[0].get_x()) +
		       (vertices[1].get_y() - vertices[0].get_y()) * (vertices[1].get_y() - vertices[0].get_y()) +
		       (vertices[1].get_z() - vertices[0].get_z()) * (vertices[1].get_z() - vertices[0].get_z()));
  int nsegments = (int)(length/dash_length);
  if((nsegments%2) == 0) nsegments++;
  double frac = 0;
  double frac2 =0.5/(double)nsegments;
  carts[0] = frac*vertices[1] + (1-frac)*vertices[0];
  carts[1] = frac2*vertices[1] + (1-frac2)*vertices[0];
  Line *line = new Line(carts,GetColour(),carts[0],GetSize(),GetAlpha());
  std::vector<Primitive*> b = line->GetSimplePrimitives(clip_vol,objrotmatrix,objorigin);
  a.insert(a.end(),b.begin(),b.end());
  for(int i=1;i<nsegments;i++){
    if(i%2 >0){
      frac = (i+0.5)/(double)nsegments;
      frac2 =(i+1.5)/(double)nsegments;
      if(frac2>1.0) frac2 = 1.0;
      carts[0] = frac*vertices[1] + (1-frac)*vertices[0];
      carts[1] = frac2*vertices[1] + (1-frac2)*vertices[0];
      Line *line = new Line(carts,colour,carts[0],size,alpha);
      std::vector<Primitive*> b = line->GetSimplePrimitives(clip_vol,objrotmatrix,objorigin);
      a.insert(a.end(),b.begin(),b.end());
    }
  }
  return a;
}

void DashLineElement::draw(const double *override_colour, int selective_override){
  double length = sqrt((vertices[1].get_x() - vertices[0].get_x()) * (vertices[1].get_x() - vertices[0].get_x()) +
		       (vertices[1].get_y() - vertices[0].get_y()) * (vertices[1].get_y() - vertices[0].get_y()) +
		       (vertices[1].get_z() - vertices[0].get_z()) * (vertices[1].get_z() - vertices[0].get_z()));
  int nsegments = (int)(length/dash_length);
  if((nsegments%2) == 1) nsegments++;
  double frac = 0;
  double frac2 =0.5/(double)nsegments;
  glVertex3f(frac*vertices[1].get_x() + (1-frac)*vertices[0].get_x(),
		 frac*vertices[1].get_y() + (1-frac)*vertices[0].get_y(),
		 frac*vertices[1].get_z() + (1-frac)*vertices[0].get_z());
  glVertex3f(frac2*vertices[1].get_x() + (1-frac2)*vertices[0].get_x(),
		 frac2*vertices[1].get_y() + (1-frac2)*vertices[0].get_y(),
		 frac2*vertices[1].get_z() + (1-frac2)*vertices[0].get_z());
  for(int i=1;i<nsegments;i++){
    if(i%2 >0){
      frac = (i+0.5)/(double)nsegments;
      frac2 =(i+1.5)/(double)nsegments;
      if(frac2>1.0) frac2 = 1.0;
      glVertex3f(frac*vertices[1].get_x() + (1-frac)*vertices[0].get_x(),
		 frac*vertices[1].get_y() + (1-frac)*vertices[0].get_y(),
		 frac*vertices[1].get_z() + (1-frac)*vertices[0].get_z());
      glVertex3f(frac2*vertices[1].get_x() + (1-frac2)*vertices[0].get_x(),
		 frac2*vertices[1].get_y() + (1-frac2)*vertices[0].get_y(),
		 frac2*vertices[1].get_z() + (1-frac2)*vertices[0].get_z());
    }
  }
}

void CircleElement::draw(const double *override_colour, int selective_override){
  bool useMatrix = false;
  matrix m = getMatrix();
  if(m.get_rows()==4&&m.get_rows()==4) useMatrix = true;

  glPushMatrix();
  glTranslatef(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  if(useMatrix) {
    glMultMatrixd(getMatrix().to_dp());
  }
  GLUquadric *qobj = gluNewQuadric();
  gluDisk(qobj,0.0,size,16,16);
  gluDeleteQuadric(qobj);
  glPopMatrix();
}

void Circle::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  CircleElement::draw(override_colour,selective_override);
}

void PartialCircleElement::draw(const double *override_colour, int selective_override){
  bool useMatrix = false;
  matrix m = getMatrix();
  if(m.get_rows()==4&&m.get_rows()==4) useMatrix = true;

  glPushMatrix();
  glTranslatef(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  if(useMatrix) {
    glMultMatrixd(getMatrix().to_dp());
  }
  GLUquadric *qobj = gluNewQuadric();
  // Use sane rather than OpenGL angles.
  gluPartialDisk(qobj,0.0,size,16,16,-startAngle+90.0,-sweepAngle);
  gluDeleteQuadric(qobj);
  glPopMatrix();
}

void PartialCircle::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  CircleElement::draw(override_colour,selective_override);
}

void TriangleElement::draw(const double *override_colour, int selective_override){
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  unsigned i=0;
  glBegin(GL_TRIANGLES);
  while(k!=vertices.end()){
    glNormal3d(normals[i].get_x(),normals[i].get_y(),normals[i].get_z());
    if(colours.size()>0) {
      glColor4d(colours[i].get_x(),colours[i].get_y(),colours[i].get_z(), alpha);
    }
    glVertex3d((*k).get_x(),(*k).get_y(),(*k).get_z());
    k++; i++;
  }
  glEnd();
}

void PolygonElement::draw(const double *override_colour, int selective_override){
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  std::vector<Cartesian>::const_iterator n = normals.begin();
  glBegin(GL_POLYGON);
  Cartesian normal = Cartesian::CrossProduct(vertices[0]-vertices[1],vertices[1]-vertices[2]);
  glNormal3f(normal.get_x(),normal.get_y(),normal.get_z());
  while(k!=vertices.end()){
    if(normals.size()==vertices.size()){
	glNormal3f((*n).get_x(),(*n).get_y(),(*n).get_z());
	n++;
    }
    glVertex3f((*k).get_x(),(*k).get_y(),(*k).get_z());
    k++;
  }
  glEnd();
}

void TriangleStripElement::draw(const double *override_colour, int selective_override){
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  glBegin(GL_TRIANGLE_STRIP);
  unsigned int i = 0;
  while(k!=vertices.end()){
    if(i<vertices.size()+3){
      Plane plane(vertices[i],vertices[i+1],vertices[i+2]);
      plane.Normalize();
      Cartesian normal = plane.get_normal();
      glNormal3f(normal.get_x(),normal.get_y(),normal.get_z());
    }
    glVertex3f((*k).get_x(),(*k).get_y(),(*k).get_z());
    k++; i++;
  }
  glEnd();
}

void TriangleFanElement::draw(const double *override_colour, int selective_override){
  if(vertices.size()<3) return;
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  glBegin(GL_TRIANGLE_FAN);
  Plane plane(vertices[0],vertices[1],vertices[2]);
  plane.Normalize();
  Cartesian normal = plane.get_normal();
  glNormal3f(normal.get_x(),normal.get_y(),normal.get_z());
  unsigned int i = 0;
  while(k!=vertices.end()-1){
    if(i>2){
      Plane plane(vertices[0],vertices[i-1],vertices[i]);
      plane.Normalize();
      Cartesian normal = plane.get_normal();
      glNormal3f(normal.get_x(),normal.get_y(),normal.get_z());
    }
    glVertex3f((*k).get_x(),(*k).get_y(),(*k).get_z());
    k++; i++;
  }

  if(vertices.size()>3){
    if((vertices[1]-vertices[i]).length()<1e-4){
      plane = Plane(vertices[0],vertices[1],vertices[2]);
    } else {
      plane = Plane(vertices[0],vertices[i-1],vertices[i]);
    }
    plane.Normalize();
    normal = plane.get_normal();
    glNormal3f(normal.get_x(),normal.get_y(),normal.get_z());
    glVertex3f((*k).get_x(),(*k).get_y(),(*k).get_z());
  }
  glEnd();
}

void QuadElement::draw(const double *override_colour, int selective_override){
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  unsigned i=0;
  glBegin(GL_QUADS);
  while(k!=vertices.end()){
    glNormal3f(normals[i].get_x(),normals[i].get_y(),normals[i].get_z());
    glVertex3f((*k).get_x(),(*k).get_y(),(*k).get_z());
    k++; i++;
  }
  glEnd();
}

void QuadStripElement::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  glBegin(GL_QUAD_STRIP);
  unsigned int i = 0;
  while(k!=vertices.end()){
    if(i%2==0&&i<vertices.size()+3){
      Plane plane(vertices[i],vertices[i+1],vertices[i+2]);
      plane.Normalize();
      Cartesian normal = plane.get_normal();
      glNormal3f(normal.get_x(),normal.get_y(),normal.get_z());
    }
    glVertex3f((*k).get_x(),(*k).get_y(),(*k).get_z());
    k++; i++;
  }
  glEnd();
}

void CylinderElement::draw(const double *override_colour, int selective_override){
  draw_cylinder(vertices,size,accu,textured);
}

void PolyCylinder::draw(const double *override_colour, int selective_override){
  if(vertices.size()<2) return;
  glEnable(GL_LIGHTING);

  std::vector<Cartesian>::iterator l = colour_vector.begin();

  int multicolour = 0;

  while(l!=colour_vector.end()){
    (*l).set_a(alpha);
    l++;
  }
  if(colour_vector.size()>0)
    multicolour = 1;

  std::vector<std::vector<Cartesian> > AB = GetPolyCylinderAB(vertices,size);
  draw_elliptical_ribbon(vertices,AB[0],AB[1],accu,textured,multicolour,colour_vector);

}

void SphereElement::draw(const double *override_colour, int selective_override){
  draw_sphere(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z(),accu,size);
}

ArrowHeadRibbon::ArrowHeadRibbon(const std::vector<Cartesian> &vertices_in, const std::vector<Cartesian> &n1_vertices_in, const std::vector<Cartesian> &n2_vertices_in, const double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double minsize_in, double maxsize_in, double thickness_in, double arrow_length_in, double arrow_width_in, double alpha_in, int accu_in, int scale_steps_in, int textured_in, int style_in, bool two_colour_in, bool edges_grey_in) : Ribbon(vertices_in,n1_vertices_in,n2_vertices_in,colour_in,origin_in,colour_vector_in,minsize_in,maxsize_in,thickness_in,alpha_in,accu_in,scale_steps_in,textured_in,style_in,two_colour_in,edges_grey_in){
  arrow_length = arrow_length_in;
  arrow_width = arrow_width_in;
}

Worm::Worm(const std::vector<Cartesian> &vertices_in, const std::vector<Cartesian> &n1_vertices_in, const std::vector<Cartesian> &n2_vertices_in, const double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double minsize_in, double maxsize_in, double thickness_in, double alpha_in, int accu_in, int scale_steps_in, int textured_in) : Ribbon(vertices_in,n1_vertices_in,n2_vertices_in,colour_in,origin_in,colour_vector_in,minsize_in,maxsize_in,thickness_in,alpha_in,accu_in,scale_steps_in,textured_in){
}

VariableWorm::VariableWorm(const std::vector<Cartesian> &vertices_in, const std::vector<Cartesian> &n1_vertices_in, const std::vector<Cartesian> &n2_vertices_in, const double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double minsize_in, double maxsize_in, double thickness_in, double alpha_in, int accu_in, int scale_steps_in, int textured_in) : Ribbon(vertices_in,n1_vertices_in,n2_vertices_in,colour_in,origin_in,colour_vector_in,minsize_in,maxsize_in,thickness_in,alpha_in,accu_in,scale_steps_in,textured_in){
}

Ribbon::Ribbon(const std::vector<Cartesian> &vertices_in, const std::vector<Cartesian> &n1_vertices_in, const std::vector<Cartesian> &n2_vertices_in, const double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double minsize_in, double maxsize_in, double thickness_in, double alpha_in, int accu_in, int scale_steps_in, int textured_in, int style_in, bool two_colour_in, bool edges_grey_in) : Primitive(){

  alpha = alpha_in;
  two_colour = two_colour_in;
  edges_grey = edges_grey_in;
  style = style_in;
  accu = accu_in;
  scale_steps = scale_steps_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;

  minsize = minsize_in;
  maxsize = maxsize_in;
  thickness = thickness_in;
  
  n1_vertices = n1_vertices_in;
  n2_vertices = n2_vertices_in;
  vertices = vertices_in;
  colour_vector = colour_vector_in;
  textured = textured_in;
}

std::vector<Cartesian> get_v_and_vpr(const Cartesian &sp1, const Cartesian &sp2, const Cartesian &ovec);

void Ribbon::SetAlpha(double alpha_in){

  alpha = alpha_in;
  std::vector<Cartesian>::iterator l = colour_vector.begin();
  while(l!=colour_vector.end()){
    (*l).set_a(alpha);
    l++;
  }
  arraysGenerated = 0;

}

std::vector<std::vector<Cartesian> > ArrowHeadRibbon::GetAB(int &insert_at_vertex) const {

  std::vector<Cartesian> vertices_tmp;
  std::vector<Cartesian> n1_vertices_tmp;
  std::vector<Cartesian> n2_vertices_tmp;
  if(RenderQuality::GetSmoothRibbons()) { 
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    n1_vertices_tmp = SplineCurve(n1_vertices,(n1_vertices.size()-1)*4,3,1);
    n2_vertices_tmp = SplineCurve(n2_vertices,(n2_vertices.size()-1)*4,3,1);
  } else {
    vertices_tmp = vertices;
    n1_vertices_tmp = n1_vertices;
    n2_vertices_tmp = n2_vertices;
  }

  /* ????
  std::vector<Cartesian>::const_iterator l = colour_vector.begin();
  while(l!=colour_vector.end()){
    (*l).set_a(alpha);
    l++;
  }
  */

  unsigned int i;
  double scale = minsize;
  double scale_step = (maxsize-minsize)/double(scale_steps-1);

  for(i=0;i<vertices_tmp.size()&&i<1;i++){
    if(scale<minsize) scale = minsize;
    n1_vertices_tmp[i].normalize(scale);
    n2_vertices_tmp[i].normalize(thickness);
  }

  double arrowlength = arrow_length;
  double arrowwidth = arrow_width;

  scale = maxsize;
  for(;(vertices_tmp[i]-vertices_tmp[vertices_tmp.size()-1]).length()>arrowlength&&i<vertices_tmp.size();i++){
    n1_vertices_tmp[i].normalize(scale);
    n2_vertices_tmp[i].normalize(thickness);
  }

  scale = arrowwidth;

  if(i>0&&i<vertices.size()){
    insert_at_vertex = int(i);
    n1_vertices_tmp.insert(n1_vertices_tmp.begin()+insert_at_vertex,n1_vertices_tmp[insert_at_vertex-1]);
    n2_vertices_tmp.insert(n2_vertices_tmp.begin()+insert_at_vertex,n2_vertices_tmp[insert_at_vertex-1]);
  } else {
    insert_at_vertex = -1;
  }

  if(n1_vertices_tmp.size()-i-1>0) scale_step = (arrowwidth-minsize)/double(n1_vertices_tmp.size()-i-1);
  else scale_step = (arrowwidth-minsize)/double(n1_vertices_tmp.size()-i);

  for(;i<n1_vertices_tmp.size()&&i>=0;i++){
    if(scale<minsize) scale = minsize;
    n1_vertices_tmp[i].normalize(scale);
    n2_vertices_tmp[i].normalize(thickness);
    scale -= scale_step;
  }
  std::vector<std::vector<Cartesian> > AB ;
  AB.push_back(n1_vertices_tmp);
  AB.push_back(n2_vertices_tmp);
  return AB;
}

void ArrowHeadRibbon::draw(const double *override_colour, int selective_override){
  if(vertices.size()<2) return;

  if (useVertexArrays) {
    drawArrays(override_colour,selective_override);
    return;
  }

  std::vector<Cartesian> vertices_tmp;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;

  if(RenderQuality::GetRenderQuality()&&0){ // Disable this for now - nice but slow.
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    if(colour_vector_tmp.size()>0) {
      std::vector<Cartesian> colour_vector_tmp_tmp;
      for(unsigned ii=0;ii<colour_vector_tmp.size();ii++) {
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
      }
      colour_vector_tmp = colour_vector_tmp_tmp; 
      while(colour_vector_tmp.size()<vertices_tmp.size()) colour_vector_tmp.push_back(colour_vector_tmp.back());
    }
  } else {
    vertices_tmp = vertices;
  }

  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);


  if(insert_at_vertex>0){
    vertices_tmp.insert(vertices_tmp.begin()+insert_at_vertex,0.9*vertices_tmp[insert_at_vertex-1] + 0.1*vertices_tmp[insert_at_vertex]);
    if(colour_vector_tmp.size()>0) colour_vector_tmp.insert(colour_vector_tmp.begin()+insert_at_vertex,colour_vector_tmp[insert_at_vertex-1]);
  }
  //std::cout << " " << vertices.size() << "\n"; std::cout.flush();
  //std::cout << " " << vertices_tmp.size() << "\n"; std::cout.flush();

  glEnable(GL_LIGHTING);
  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;

  int my_accu = accu;

  if(RenderQuality::GetRenderQuality())
    my_accu = 180;
  if(style==0)
    draw_elliptical_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp,false);
  if(style==1)
    draw_flat_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp,false,edges_grey);
  else if(style==2)
    draw_flat_rounded_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp);

}
std::vector<std::vector<Cartesian> > Ribbon::GetAB(int &insert_at_vertex) const {

  std::vector<Cartesian> vertices_tmp;
  std::vector<Cartesian> n1_vertices_tmp;
  std::vector<Cartesian> n2_vertices_tmp;
  if(RenderQuality::GetSmoothRibbons()) {
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    n1_vertices_tmp = SplineCurve(n1_vertices,(n1_vertices.size()-1)*4,3,1);
    n2_vertices_tmp = SplineCurve(n2_vertices,(n2_vertices.size()-1)*4,3,1);
  } else {
    vertices_tmp = vertices;
    n1_vertices_tmp = n1_vertices;
    n2_vertices_tmp = n2_vertices;
  }

  /* ????
  std::vector<Cartesian>::const_iterator l = colour_vector.begin();
  while(l!=colour_vector.end()){
    (*l).set_a(alpha);
    l++;
  }
  */

  unsigned int i;
  double scale = minsize;

  int scale_steps_tmp = scale_steps;
  if(RenderQuality::GetSmoothRibbons())
    scale_steps_tmp *= 4;

  double scale_step = (maxsize-minsize)/double(scale_steps_tmp);

  for(i=0;i<vertices_tmp.size()&&int(i)<scale_steps_tmp;i++){
    if(scale<minsize) scale = minsize;
    n1_vertices_tmp[i].normalize(scale);
    n2_vertices_tmp[i].normalize(thickness);
    scale += scale_step;
  }
  scale = maxsize;
  for(i=scale_steps_tmp;i<vertices_tmp.size()-scale_steps_tmp-1&&i<vertices_tmp.size();i++){
    n1_vertices_tmp[i].normalize(scale);
    n2_vertices_tmp[i].normalize(thickness);
  }
  for(i=vertices_tmp.size()-scale_steps_tmp-1;i<vertices_tmp.size()&&i>=0;i++){
    if(scale<minsize) scale = minsize;
    n1_vertices_tmp[i].normalize(scale);
    n2_vertices_tmp[i].normalize(thickness);
    scale -= scale_step;
  }
  std::vector<std::vector<Cartesian> > AB ;
  AB.push_back(n1_vertices_tmp);
  AB.push_back(n2_vertices_tmp);
  return AB;
}

std::vector<std::vector<Cartesian> > VariableWorm::GetAB(int &insert_at_vertex) const {

  std::vector<Cartesian> vertices_tmp;
  std::vector<Cartesian> n1_vertices_tmp;
  std::vector<Cartesian> n2_vertices_tmp;
  if(RenderQuality::GetSmoothRibbons()) {
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    n1_vertices_tmp = SplineCurve(n1_vertices,(n1_vertices.size()-1)*4,3,1);
    n2_vertices_tmp = SplineCurve(n2_vertices,(n2_vertices.size()-1)*4,3,1);
  } else {
    vertices_tmp = vertices;
    n1_vertices_tmp = n1_vertices;
    n2_vertices_tmp = n2_vertices;
  }

  for(unsigned i=0;i<vertices_tmp.size();i++){
    if(n1_vertices_tmp[i].length()<minsize){
      n1_vertices_tmp[i].normalize(minsize);
    }
    if(n2_vertices_tmp[i].length()<minsize){
      n2_vertices_tmp[i].normalize(minsize);
    }
  }
  std::vector<std::vector<Cartesian> > AB ;
  AB.push_back(n1_vertices_tmp);
  AB.push_back(n2_vertices_tmp);
  return AB;
}

size_t PolyCylinder::GetNumberOfSimplePrimitives() const {
  return vertices.size()*accu;
}

size_t Ribbon::GetNumberOfSimplePrimitives() const {
  size_t nprims;
  if(vertices.size()<1){
    return 0;
  }
  nprims = (vertices.size()-1)*accu;
  if(RenderQuality::GetSmoothRibbons())
    nprims = nprims * 4 - accu;
  return nprims;
}

size_t ArrowHeadRibbon::GetNumberOfSimplePrimitives() const {
  int insert_at_vertex;
  size_t nprims;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);
  if (insert_at_vertex > -1) insert_at_vertex = 0;
  nprims = (vertices.size()+insert_at_vertex)*accu;
  if(RenderQuality::GetSmoothRibbons())
    nprims = nprims * 4 - accu;
  return nprims;
}

void Ribbon::draw(const double *override_colour, int selective_override){
  if(vertices.size()<2) return;

  if (useVertexArrays) {
    drawArrays(override_colour,selective_override);
    return;
  }

  std::vector<Cartesian> vertices_tmp;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;

  if(RenderQuality::GetRenderQuality()&&0){ // Disable this for now - nice but slow.
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    if(colour_vector_tmp.size()>0) {
      std::vector<Cartesian> colour_vector_tmp_tmp;
      for(unsigned ii=0;ii<colour_vector_tmp.size();ii++) {
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
      }
      colour_vector_tmp = colour_vector_tmp_tmp; 
      while(colour_vector_tmp.size()<vertices_tmp.size()) colour_vector_tmp.push_back(colour_vector_tmp.back());
    }
  } else {
    vertices_tmp = vertices;
  }
  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);
  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;
  glEnable(GL_LIGHTING);

  int my_accu = accu;
  if(RenderQuality::GetRenderQuality())
      my_accu = 180;
  if(style==0)
    draw_elliptical_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp,two_colour);
  else if(style==1)
    draw_flat_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp,two_colour,edges_grey);
  else if(style==2)
    draw_flat_rounded_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp,two_colour);
  else if(style==3)
    draw_fancy_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp,two_colour);
}

void Worm::draw(const double *override_colour, int selective_override){
  if(vertices.size()<2) return;

  if (useVertexArrays) {
    drawArrays(override_colour,selective_override);
    return;
  }

  std::vector<Cartesian> vertices_tmp;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;

  if(RenderQuality::GetRenderQuality()&&0){ // Disable this for now - nice but slow.
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    if(colour_vector_tmp.size()>0) {
      std::vector<Cartesian> colour_vector_tmp_tmp;
      for(unsigned ii=0;ii<colour_vector_tmp.size();ii++) {
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
      }
      colour_vector_tmp = colour_vector_tmp_tmp; 
      while(colour_vector_tmp.size()<vertices_tmp.size()) colour_vector_tmp.push_back(colour_vector_tmp.back());
    }
  } else {
    vertices_tmp = vertices;
  }
  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = Ribbon::GetAB(insert_at_vertex);
  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;
  glEnable(GL_LIGHTING);

  int my_accu;
  if(RenderQuality::GetRenderQuality())
    my_accu = 180;
  else {
    if(accu<4)
      my_accu = 8;
    else
      my_accu = accu;
  }
  draw_elliptical_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp,false);
}

void VariableWorm::draw(const double *override_colour, int selective_override){
  if(vertices.size()<2) return;

  if (useVertexArrays) {
    drawArrays(override_colour,selective_override);
    return;
  }

  std::vector<Cartesian> vertices_tmp;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;

  if(RenderQuality::GetRenderQuality()&&0){ // Disable this for now - nice but slow.
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    if(colour_vector_tmp.size()>0) {
      std::vector<Cartesian> colour_vector_tmp_tmp;
      for(unsigned ii=0;ii<colour_vector_tmp.size();ii++) {
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
        colour_vector_tmp_tmp.push_back(colour_vector_tmp[ii]);
      }
      colour_vector_tmp = colour_vector_tmp_tmp; 
      while(colour_vector_tmp.size()<vertices_tmp.size()) colour_vector_tmp.push_back(colour_vector_tmp.back());
    }
  } else {
    vertices_tmp = vertices;
  }
  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);
  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;
  glEnable(GL_LIGHTING);

  int my_accu;
  if(RenderQuality::GetRenderQuality())
    my_accu = 180;
  else {
    if(accu<4)
      my_accu = 8;
    else
      my_accu = accu;
  }
  draw_elliptical_ribbon(vertices_tmp,AB[0],AB[1],my_accu,textured,multicolour,colour_vector_tmp,false);
}


BillBoard::BillBoard(const std::vector<Cartesian> &vertices_in, const Cartesian &origin_in, const char* filename_in, int style_in){
  origin = origin_in;

  vertices = vertices_in;
  style = style_in;
  first_scale_w = 1.0;
  first_scale_h = 1.0;
  scale_w = 1.0;
  scale_h = 1.0;
  
  image = image_info(filename_in);
  if(image.get_width()&&image.get_height()&&image.get_pixels()){
    image.invert();
    filename = new char[1024];
    strncpy(filename,filename_in,1024);
  }

  texture_id = 0;
  mask_texture_id = 0;

}

void ForceLoadTextures(){ force_load_textures = true; };
void UnForceLoadTextures(){ force_load_textures = false; };

void BillBoard::draw(const double *override_colour, int selective_override){

   if(!image.get_width()||!image.get_height()||!image.get_pixels())
      return;
   glDisable(GL_CLIP_PLANE0);
   glDisable(GL_CLIP_PLANE1);
   GLint poly_params[2];
   glGetIntegerv(GL_POLYGON_MODE,poly_params);
   if(poly_params[0]==GL_LINE)
     return;

   GLint viewport[4];
   glGetIntegerv(GL_VIEWPORT,viewport);
   
   if(texture_id!=0&&force_load_textures) { 
     texture_id_backup = texture_id;
     mask_texture_id_backup = mask_texture_id;
   }

   if(texture_id==0||force_load_textures) { 
     std::pair<unsigned,unsigned> wh = GetCompatibleTextureSize(image.get_width(),image.get_height());
     image.convert_rgba();
     image = ResizeWithEmptySpace(image,wh.first,wh.second);
     //image.write("image.png");
     image_info mask = image.GenerateMask();
     //mask.write("mask.png");

     if(image.get_colourspace_type()==IMAGEINFO_MONOA){
       image.convert_rgba();
     }
     if(image.get_colourspace_type()==IMAGEINFO_RGBA||image.get_colourspace_type()==IMAGEINFO_RGB){
       //mask_texture_id = load_texture(image.GenerateMask(),MIPMAP);
     }
     texture_id = load_texture(image,MIPMAP);
     first_scale_w = double(image.get_width())/double(viewport[2]);
     first_scale_h = double(image.get_height())/double(viewport[3]);
   }

   glEnable(GL_TEXTURE_2D);
   glDisable(GL_LIGHTING);
   glPushAttrib(GL_FOG_BIT);
   glDisable(GL_FOG);

   glEnable(GL_BLEND);
   if(mask_texture_id){
     GLfloat params[4];
     glGetFloatv(GL_COLOR_CLEAR_VALUE,params);
     GLfloat bgcolor[4]={1.0,1.0,1.0,1.0};
     glColor4fv(bgcolor);
     glDisable(GL_DEPTH_TEST);
     glBlendFunc(GL_DST_COLOR, GL_ZERO);
   }else{
     glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
     GLfloat bgcolor[4]= {1.0,1.0,1.0,1.0};
     glColor4fv(bgcolor);
     glMaterialfv(GL_FRONT, GL_SPECULAR,bgcolor);
     glMaterialfv(GL_FRONT, GL_EMISSION,bgcolor);
   }

   glPushMatrix();

   glLoadIdentity();

   double projMatrix[16];
   glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);

   double l = (projMatrix[12] - 1.0)/projMatrix[0]; 
   double r = (projMatrix[12] + 1.0)/projMatrix[0];
   double t = (projMatrix[13] - 1.0)/projMatrix[5];
   double b = (projMatrix[13] + 1.0)/projMatrix[5];
   
   double win_scale_x = fabs(double(viewport[2] - viewport[0])); // Keep original size
   double win_scale_y = fabs(double(viewport[3] - viewport[1]));

   int mult = SimpleBillBoard::GetMagnification();

   double x1 = (vertices[0].get_x()-.5)*(r-l)*mult;
   double y1 = (vertices[0].get_y()-.5)*(b-t)*mult;
   double x2 = (vertices[0].get_x()-.5+scale_w*image.get_width()/win_scale_x)*(r-l)*mult;
   double y2 = (vertices[0].get_y()-.5+scale_h*image.get_height()/win_scale_y)*(b-t)*mult;

   double n = 1;
   double f = 1000;
   if(fabs(projMatrix[10])>1e-7){
     n = (projMatrix[14] + 1)/projMatrix[10];
     f = (projMatrix[14] - 1)/projMatrix[10];
   }

   if(mask_texture_id){
     glBindTexture(GL_TEXTURE_2D, mask_texture_id);
     glBegin(GL_QUADS);
     glTexCoord2f(0.0f, 0.0f);  
     glVertex3d(x1,y1,-n-3);
     glTexCoord2f(1.0f, 0.0f);  
     glVertex3d(x2,y1,-n-3);
     glTexCoord2f(1.0f, 1.0f);  
     glVertex3d(x2,y2,-n-3);
     glTexCoord2f(0.0f, 1.0f);  
     glVertex3d(x1,y2,-n-3);
     glEnd();
   }

   if(mask_texture_id)
     glBlendFunc(GL_ONE, GL_ONE);

   glBindTexture(GL_TEXTURE_2D, texture_id);
   glBegin(GL_QUADS);
   glTexCoord2f(0.0f, 0.0f);  
   glVertex3d(x1,y1,-n-3);
   glTexCoord2f(1.0f, 0.0f);  
   glVertex3d(x2,y1,-n-3);
   glTexCoord2f(1.0f, 1.0f);  
   glVertex3d(x2,y2,-n-3);
   glTexCoord2f(0.0f, 1.0f);  
   glVertex3d(x1,y2,-n-3);
   glEnd();

   glPopMatrix();
   glEnable(GL_DEPTH_TEST);
   glDisable(GL_BLEND);
   glDisable(GL_TEXTURE_2D);
   glEnable(GL_LIGHTING);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glPopAttrib();

   if(texture_id!=0&&force_load_textures) { 
     texture_id = texture_id_backup;
     mask_texture_id = mask_texture_id_backup;
   }
   glEnable(GL_CLIP_PLANE0);
   glEnable(GL_CLIP_PLANE1);
}

BillBoard::~BillBoard(){
}

void Primitive::OutputPovrayPigment(std::ofstream &fp){
  if(get_transparent()!=0)
    fp << "    pigment { color rgbt <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha  << "> }\n"; 
  else
    fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
}

void Point::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))) {
    return;
  }

  fp << "sphere {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/50.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void Line::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  fp << "cylinder {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/50.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "sphere {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/50.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "sphere {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/50.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void LineStrip::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  /*
  // More complicated ....
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }
  */

  std::vector<Cartesian>::iterator viter = vertices.begin();
  while(viter!=vertices.end()-1){
    fp << "cylinder {\n";
    Cartesian p = quat.getMatrix()*(objrotmatrix*(*viter+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    p = quat.getMatrix()*(objrotmatrix*(*(viter+1)+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    fp << "  " << GetSize()/50.0 << "\n";
    fp << "  texture {\n";
    OutputPovrayPigment(fp); 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";
    fp << "sphere {\n";

    p = quat.getMatrix()*(objrotmatrix*(*viter+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    fp << "  " << GetSize()/50.0 << "\n";
    fp << "  texture {\n";
    OutputPovrayPigment(fp); 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";

    viter++;
  }
}

void DashLine::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  double length = (vertices[0]-vertices[1]).length();
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2)==0) nsegments++;

  double frac = 0;
  double frac2 = 0.5/(double)nsegments;
  Cartesian v1 = frac*vertices[1]+(1-frac)*vertices[0];
  Cartesian v2 = frac2*vertices[1]+(1-frac2)*vertices[0];
  fp << "cylinder {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/50.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  for(int i=1;i<nsegments;i++){
    if(i%2 >0){
      double frac = (i+0.5)/(double)nsegments;
      double frac2 =(i+1.5)/(double)nsegments;
      if(frac2>1.0) frac2 = 1.0;
      fp << "cylinder {\n";
      Cartesian v1 = frac*vertices[1] + (1-frac)*vertices[0];
      Cartesian v2 = frac2*vertices[1] + (1-frac2)*vertices[0];
      p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      fp << "  " << GetSize()/50.0 << "\n";
      fp << "  texture {\n";
      OutputPovrayPigment(fp); 
      fp << "    finish {diffuse 1.0 specular 1.0}\n";
      fp << "  }\n";
      fp << "}\n";
      p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      fp << "  " << GetSize()/50.0 << "\n";
      fp << "  texture {\n";
      OutputPovrayPigment(fp); 
      fp << "    finish {diffuse 1.0 specular 1.0}\n";
      fp << "  }\n";
      fp << "}\n";
      p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      fp << "  " << GetSize()/50.0 << "\n";
      fp << "  texture {\n";
      OutputPovrayPigment(fp); 
      fp << "    finish {diffuse 1.0 specular 1.0}\n";
      fp << "  }\n";
      fp << "}\n";
    }
  }
}

void CircleElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
   std::cout << "Circle::DrawPovray not yet implemented\n";
}

void PolygonElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  if(vertices.size()==5){
  fp << "triangle {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[2]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "triangle {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[2]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[3]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "triangle {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[3]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[4]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
  }else if(vertices.size()==6){
  fp << "triangle {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[2]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "triangle {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[2]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[3]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "triangle {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[3]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[5]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "triangle {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[3]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[4]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[5]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
  }else{
     std::cout << "Polygon::DrawPovray for more than 6 sides not yet implemented\n";
  }
}

void TriangleElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  fp << "triangle {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[2]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void QuadElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  fp << "triangle {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[3]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "triangle {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[2]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[3]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void QuadStripElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  /*
  // More complicated ...
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }
  */

  std::vector<Cartesian>::iterator viter = vertices.begin();
  while(viter!=vertices.end()-3){
    fp << "polygon {\n";
    fp << vertices.size() << ",\n";
    Cartesian p = quat.getMatrix()*(objrotmatrix*(*viter+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    p = quat.getMatrix()*(objrotmatrix*(*(viter+1)+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    p = quat.getMatrix()*(objrotmatrix*(*(viter+2)+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    p = quat.getMatrix()*(objrotmatrix*(*(viter+3)+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
    fp << "  texture {\n";
    OutputPovrayPigment(fp); 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";
    viter+=2;
  }
}

void TriangleStripElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  /*
  // More complicated ...
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }
  */

  std::vector<Cartesian>::iterator viter = vertices.begin();
  while(viter!=vertices.end()-2){
    fp << "triangle {\n";
    Cartesian p = quat.getMatrix()*(objrotmatrix*(*viter+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    p = quat.getMatrix()*(objrotmatrix*(*(viter+1)+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    p = quat.getMatrix()*(objrotmatrix*(*(viter+2)+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
    fp << "  texture {\n";
    OutputPovrayPigment(fp); 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";
    viter++;
  }
}

void TriangleFanElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  /*
  // More complicated ...
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }
  */

  std::vector<Cartesian>::iterator viter = vertices.begin()+1;
  while(viter!=vertices.end()-1){
    fp << "triangle {\n";
    Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    p = quat.getMatrix()*(objrotmatrix*(*(viter)+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
    p = quat.getMatrix()*(objrotmatrix*(*(viter+1)+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
    fp << "  texture {\n";
    OutputPovrayPigment(fp); 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";
    viter++;
  }
}

void PolyCylinder::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  /*
  // More complicated ...
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }
  */

  fp << "sphere_sweep {\nlinear_spline\n";
  fp << vertices.size() << ",\n";
  std::vector<Cartesian>::const_iterator cart=vertices.begin();
  while(cart!=vertices.end()){
    Cartesian p = quat.getMatrix()*(objrotmatrix*(*cart+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
    fp << "  " << GetSize() << "\n";
    cart++;
  }
  fp << "  texture {\n";
  if(colour_vector.size()>0){
    fp << "    pigment { color rgb <" << colour_vector[0].get_x() << ", " << colour_vector[1].get_y() << ", " << colour_vector[2].get_z() << "> }\n"; 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  }else{
    OutputPovrayPigment(fp); 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
  }
  fp << "  }\n";
  fp << "}\n";
}

void SphereElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))) {
    return;
  }

  fp << "sphere {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize() << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void SpheroidElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))) {
    return;
  }

  /*
  Quat qpov(Cartesian(1,0,0),1,90);

  Quat qpov2 = qpov;
  Quat qpov3 = qpov;

  Quat rotx(Cartesian(1,0,0),1,-90);
  Quat roty(Cartesian(0,1,0),1,-90);
  Quat rotz(Cartesian(0,0,1),1,-90);

  qpov.postMult(quat);
  qpov2.postMult(quat);
  qpov3.postMult(quat);

  qpov.postMult(roty);
  qpov2.postMult(rotz);
  qpov3.postMult(rotx);
  
  matrix m = qpov.getMatrix();
  matrix m2 = qpov2.getMatrix();
  matrix m3 = qpov3.getMatrix();

  if(show_axes) {
    fp << "torus {\n";
    Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    fp << "  " << GetSize()+.1 << ", " << .1 << "\n";
    fp << "  matrix <";
    fp << "  1,0,0,";
    fp << "  0,1,0,";
    fp << "  0,0,-1,";
    fp << "  0,0,0>";
    fp << "  matrix <";
    fp << m(0,0) << ", " << m(0,1) << ", " << m(0,2) << ",\n";
    fp << m(1,0) << ", " << m(1,1) << ", " << m(1,2) << ",\n";
    fp << m(2,0) << ", " << m(2,1) << ", " << m(2,2) << ",\n";
    fp << m(3,0) << ", " << m(3,1) << ", " << m(3,2) << ">\n";
    fp << "  scale <" << a << "," << b << "," << c << ">\n";
    fp << "  translate< " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
    fp << "  texture {\n";
    //OutputPovrayPigment(fp); 
    fp << "    pigment { color rgb <" << 1 << ", " << 0 << ", " << 0 << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";

    fp << "torus {\n";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    fp << "  " << GetSize()+.1 << ", " << .1 << "\n";
    fp << "  matrix <";
    fp << "  1,0,0,";
    fp << "  0,1,0,";
    fp << "  0,0,-1,";
    fp << "  0,0,0>";
    fp << "  matrix <";
    fp << m2(0,0) << ", " << m2(0,1) << ", " << m2(0,2) << ",\n";
    fp << m2(1,0) << ", " << m2(1,1) << ", " << m2(1,2) << ",\n";
    fp << m2(2,0) << ", " << m2(2,1) << ", " << m2(2,2) << ",\n";
    fp << m2(3,0) << ", " << m2(3,1) << ", " << m2(3,2) << ">\n";
    fp << "  scale <" << a << "," << c << "," << b << ">\n";
    fp << "  translate< " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
    fp << "  texture {\n";
    //OutputPovrayPigment(fp); 
    fp << "    pigment { color rgb <" << 0 << ", " << 1 << ", " << 0 << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";

    fp << "torus {\n";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    fp << "  " << GetSize()+.1 << ", " << .1 << "\n";
    fp << "  matrix <";
    fp << "  1,0,0,";
    fp << "  0,1,0,";
    fp << "  0,0,-1,";
    fp << "  0,0,0>";
    fp << "  matrix <";
    fp << m3(0,0) << ", " << m3(0,1) << ", " << m3(0,2) << ",\n";
    fp << m3(1,0) << ", " << m3(1,1) << ", " << m3(1,2) << ",\n";
    fp << m3(2,0) << ", " << m3(2,1) << ", " << m3(2,2) << ",\n";
    fp << m3(3,0) << ", " << m3(3,1) << ", " << m3(3,2) << ">\n";
    fp << "  scale <" << b << "," << c << "," << a << ">\n";
    fp << "  translate< " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">\n";
    fp << "  texture {\n";
    //OutputPovrayPigment(fp); 
    fp << "    pigment { color rgb <" << 0 << ", " << 0 << ", " << 1 << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";
  }
  */

  //if(show_solid) {
    fp << "sphere {\n";
    Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    fp << "  < " << p.get_x()/a << ", " << p.get_y()/b << ", " << -p.get_z()/c << ">,\n";
    fp << "  " << GetSize() << "\n";
    fp << "  scale <" << a << "," << b << "," << c << ">\n";
    fp << "  texture {\n";
    OutputPovrayPigment(fp); 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
    fp << "}\n";
  //}
}

void CylinderElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  fp << "cylinder {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize() << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void ConeElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  fp << "cone {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize() << "\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << 0.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void draw_elliptical_ribbon_pov(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const std::vector<Cartesian> &vertices, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, const double *colour, const std::vector<Cartesian> &colour_vector, int transparent, double alpha){

  std::vector<Cartesian>::const_iterator cart;
  std::vector<Cartesian>::const_iterator pviter;
  std::vector<Cartesian>::const_iterator pvpriter;
  
  int nsectors = 180;

  for(unsigned i=0; i<colour_vector.size();i++){
    if(transparent)
      fp << "#declare texture_" << i << " = texture { pigment { color rgbt <" << colour_vector[i].get_x() << ", " << colour_vector[i].get_y() << ", " << colour_vector[i].get_z() << ", " << 1.0-alpha << "> }  finish {diffuse 1.0 specular 1.0} }\n";
    else
      fp << "#declare texture_" << i << " = texture { pigment { color rgb <" << colour_vector[i].get_x() << ", " << colour_vector[i].get_y() << ", " << colour_vector[i].get_z() << "> }  finish {diffuse 1.0 specular 1.0} }\n";
  }

  fp << "mesh {\n";
  double theta,theta2;

  for(int i=0;i<360;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    theta2 = (double)(i-360/nsectors)/360.0 * PIBY2;

    cart=vertices.begin();
    pviter=pv.begin();
    pvpriter=pvpr.begin();

    int icol = 0;
    while(cart!=vertices.end()-1&&pviter!=pv.end()-1&&pvpriter!=pvpr.end()-1){

      double x = cos(theta);
      double y = sin(theta);
      double x2 = cos(theta2);
      double y2 = sin(theta2);

      double a = pviter->length();
      double b = pvpriter->length();
      double a2 = (pviter+1)->length();
      double b2 = (pvpriter+1)->length();

      Cartesian n1 = x*(*pviter)/(a*a) + y*(*pvpriter)/(b*b);
      n1 = quat.getMatrix()*(objrotmatrix*(n1));
      Cartesian n2 = x*(*(pviter+1))/(a2*a2) + y*(*(pvpriter+1))/(b2*b2);
      n2 = quat.getMatrix()*(objrotmatrix*(n2));
      Cartesian n12 = x2*(*pviter)/(a*a) + y2*(*pvpriter)/(b*b);
      n12 = quat.getMatrix()*(objrotmatrix*(n12));
      Cartesian n22 = x2*(*(pviter+1))/(a2*a2) + y2*(*(pvpriter+1))/(b2*b2);
      n22 = quat.getMatrix()*(objrotmatrix*(n22));

      fp << "smooth_triangle {\n";
      
      Cartesian p1 = x*(*pviter) + y*(*pvpriter) + *cart;
      Cartesian p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x*(*(pviter+1)) + y*(*(pvpriter+1)+objorigin) + *(cart+1);
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">,";
      
      p1 = x2*(*(pviter+1)) + y2*(*(pvpriter+1)+objorigin) + *(cart+1);
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">";
      
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";
      
      fp << "smooth_triangle {\n";
      
      p1 = x*(*pviter) + y*(*pvpriter) + *cart;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x2*(*(pviter+1)) + y2*(*(pvpriter+1)) + *(cart+1);
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">,";
      
      p1 = x2*(*pviter) + y2*(*pvpriter) + *cart;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n12.get_x() << ", " << n12.get_y() << ", " << -n12.get_z() << ">";

      fp << "   texture_list { texture_" << icol+1 << " texture_" << icol+1 << " texture_" << icol+1 << " } }\n";

      cart++; pviter++; pvpriter++;
      icol++;
    }
  }
  if(colour_vector.size()==0){
    fp << "  texture {\n";
    if(transparent) fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha << "> }\n"; 
    else fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
  }
  fp << "}\n";
}

void draw_flat_round_ribbon_pov(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const std::vector<Cartesian> &vertices, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, const double *colour, const std::vector<Cartesian> &colour_vector, int transparent, double alpha, const bool two_colour){

  std::vector<Cartesian>::const_iterator cart;
  std::vector<Cartesian>::const_iterator pviter;
  std::vector<Cartesian>::const_iterator pvpriter;
  
  int idx = vertices.size()/2;
  Cartesian p_t = vertices[idx+1] + vertices[idx-1] - 2*vertices[idx];
  p_t.normalize();
  Cartesian pvpr_t = pvpr[idx];
  pvpr_t.normalize();

  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;

  int nsectors = 180;

  for(unsigned i=0; i<colour_vector.size();i++){
    if(transparent)
      fp << "#declare texture_" << i << " = texture { pigment { color rgbt <" << colour_vector[i].get_x() << ", " << colour_vector[i].get_y() << ", " << colour_vector[i].get_z() << ", " << 1.0-alpha << "> }  finish {diffuse 1.0 specular 1.0} }\n";
    else
      fp << "#declare texture_" << i << " = texture { pigment { color rgb <" << colour_vector[i].get_x() << ", " << colour_vector[i].get_y() << ", " << colour_vector[i].get_z() << "> }  finish {diffuse 1.0 specular 1.0} }\n";
  }
  fp << "#declare texture_grey = texture { pigment { color rgb <0.6, 0.6 , 0.6 > }  finish {diffuse 1.0 specular 1.0} }\n";

  fp << "mesh {\n";

    cart=vertices.begin();
    pviter=pv.begin();
    pvpriter=pvpr.begin();

    int icol = 0;
    while(cart!=vertices.end()-1&&pviter!=pv.end()-1&&pvpriter!=pvpr.end()-1){
      
      double a = pvpriter->length();
      double a2 = (pvpriter+1)->length();
      Cartesian pvshift = (pviter->length()-a)/pviter->length() * *pviter;
      Cartesian pv1shift = ((pviter+1)->length()-a2)/(pviter+1)->length() * *(pviter+1);

      Cartesian p1 = pvshift + *pvpriter + *cart;
      Cartesian p2 = pvshift - *pvpriter + *cart;
      Cartesian p3 = pv1shift - *(pvpriter+1) + *(cart+1);
      Cartesian p4 = pv1shift + *(pvpriter+1) + *(cart+1);

      Cartesian p5 = -pvshift + *pvpriter + *cart;
      Cartesian p6 = -pvshift - *pvpriter + *cart;
      Cartesian p7 = -pv1shift - *(pvpriter+1) + *(cart+1);
      Cartesian p8 = -pv1shift + *(pvpriter+1) + *(cart+1);

      Cartesian n1 = *pviter;
      Cartesian n2 = *(pviter+1);
      Cartesian n3 = *pvpriter;
      Cartesian n4 = *(pvpriter+1);
      n1.normalize();
      n2.normalize();
      n3.normalize();
      n4.normalize();
      n1 =  quat.getMatrix()*(objrotmatrix*(n1));
      n2 =  quat.getMatrix()*(objrotmatrix*(n2));
      n3 =  quat.getMatrix()*(objrotmatrix*(n3));
      n4 =  quat.getMatrix()*(objrotmatrix*(n4));

      Cartesian p;
      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p5+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p8+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      if(face_one&&two_colour)
        fp << "   texture_list { texture_grey texture_grey texture_grey } }\n";
      else
        fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p4+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p8+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      if(face_one&&two_colour)
        fp << "   texture_list { texture_grey texture_grey texture_grey } }\n";
      else
        fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";
      
      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p6+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p7+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      if(!face_one&&two_colour)
        fp << "   texture_list { texture_grey texture_grey texture_grey } }\n";
      else
        fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p7+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      if(!face_one&&two_colour)
        fp << "   texture_list { texture_grey texture_grey texture_grey } }\n";
      else
        fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      cart++; pviter++; pvpriter++;
      icol++;
    }

  if(colour_vector.size()==0){
    fp << "  texture {\n";
    if(transparent) fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha << "> }\n"; 
    else fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
  }
  fp << "}\n";

  fp << "mesh {\n";
  double theta,theta2;

  for(int i=0;i<360;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    theta2 = (double)(i-360/nsectors)/360.0 * PIBY2;

    cart=vertices.begin();
    pviter=pv.begin();
    pvpriter=pvpr.begin();

    int icol = 0;
    while(cart!=vertices.end()-1&&pviter!=pv.end()-1&&pvpriter!=pvpr.end()-1){

      double x = cos(theta);
      double y = sin(theta);
      double x2 = cos(theta2);
      double y2 = sin(theta2);

      double a = pvpriter->length();
      double b = pvpriter->length();
      double a2 = (pvpriter+1)->length();
      double b2 = (pvpriter+1)->length();

      Cartesian pvn = *pviter;
      Cartesian pvprn = *pvpriter;
      pvn.normalize(a);
      pvprn.normalize(a);
      Cartesian pv1n = *(pviter+1);
      Cartesian pvpr1n = *(pvpriter+1);
      pv1n.normalize(a);
      pvpr1n.normalize(a);

      Cartesian pvshift = (pviter->length()-a)/pviter->length() * *pviter;
      Cartesian pv1shift = ((pviter+1)->length()-a2)/(pviter+1)->length() * *(pviter+1);

      Cartesian n1 = x*(pvn)/(a*a) + y*(pvprn)/(b*b);
      n1 = quat.getMatrix()*(objrotmatrix*(n1));
      Cartesian n2 = x*(pv1n)/(a2*a2) + y*(*(pvpriter+1))/(b2*b2);
      n2 = quat.getMatrix()*(objrotmatrix*(n2));
      Cartesian n12 = x2*(pvn)/(a*a) + y2*(pvprn)/(b*b);
      n12 = quat.getMatrix()*(objrotmatrix*(n12));
      Cartesian n22 = x2*(pv1n)/(a2*a2) + y2*(*(pvpriter+1))/(b2*b2);
      n22 = quat.getMatrix()*(objrotmatrix*(n22));

      fp << "smooth_triangle {\n";
      
      Cartesian p1 = x*(pvn) + y*(pvprn) + *cart + pvshift;
      Cartesian p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x*(pv1n) + y*(pvpr1n+objorigin) + *(cart+1) + pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">,";
      
      p1 = x2*(pv1n) + y2*(pvpr1n+objorigin) + *(cart+1) + pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">";
      
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";
      
      fp << "smooth_triangle {\n";
      
      p1 = x*(pvn) + y*(pvprn) + *cart + pvshift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x2*(pv1n) + y2*(pvpr1n) + *(cart+1) + pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">,";
      
      p1 = x2*(pvn) + y2*(pvprn) + *cart + pvshift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n12.get_x() << ", " << n12.get_y() << ", " << -n12.get_z() << ">";

      fp << "   texture_list { texture_" << icol+1 << " texture_" << icol+1 << " texture_" << icol+1 << " } }\n";

      cart++; pviter++; pvpriter++;
      icol++;
    }
  }
  if(colour_vector.size()==0){
    fp << "  texture {\n";
    if(transparent) fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha << "> }\n"; 
    else fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
  }
  fp << "}\n";
  fp << "mesh {\n";

  for(int i=0;i<360;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    theta2 = (double)(i-360/nsectors)/360.0 * PIBY2;

    cart=vertices.begin();
    pviter=pv.begin();
    pvpriter=pvpr.begin();

    int icol = 0;
    while(cart!=vertices.end()-1&&pviter!=pv.end()-1&&pvpriter!=pvpr.end()-1){

      double x = cos(theta);
      double y = sin(theta);
      double x2 = cos(theta2);
      double y2 = sin(theta2);

      double a = pvpriter->length();
      double b = pvpriter->length();
      double a2 = (pvpriter+1)->length();
      double b2 = (pvpriter+1)->length();

      Cartesian pvn = *pviter;
      Cartesian pvprn = *pvpriter;
      pvn.normalize(a);
      pvprn.normalize(a);
      Cartesian pv1n = *(pviter+1);
      Cartesian pvpr1n = *(pvpriter+1);
      pv1n.normalize(a);
      pvpr1n.normalize(a);

      Cartesian pvshift = (pviter->length()-a)/pviter->length() * *pviter;
      Cartesian pv1shift = ((pviter+1)->length()-a2)/(pviter+1)->length() * *(pviter+1);

      Cartesian n1 = x*(pvn)/(a*a) + y*(pvprn)/(b*b);
      n1 = quat.getMatrix()*(objrotmatrix*(n1));
      Cartesian n2 = x*(pv1n)/(a2*a2) + y*(*(pvpriter+1))/(b2*b2);
      n2 = quat.getMatrix()*(objrotmatrix*(n2));
      Cartesian n12 = x2*(pvn)/(a*a) + y2*(pvprn)/(b*b);
      n12 = quat.getMatrix()*(objrotmatrix*(n12));
      Cartesian n22 = x2*(pv1n)/(a2*a2) + y2*(*(pvpriter+1))/(b2*b2);
      n22 = quat.getMatrix()*(objrotmatrix*(n22));

      fp << "smooth_triangle {\n";
      
      Cartesian p1 = x*(pvn) + y*(pvprn) + *cart - pvshift;
      Cartesian p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x*(pv1n) + y*(pvpr1n+objorigin) + *(cart+1) - pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">,";
      
      p1 = x2*(pv1n) + y2*(pvpr1n+objorigin) + *(cart+1) - pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">";
      
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";
      
      fp << "smooth_triangle {\n";
      
      p1 = x*(pvn) + y*(pvprn) + *cart - pvshift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x2*(pv1n) + y2*(pvpr1n) + *(cart+1) - pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">,";
      
      p1 = x2*(pvn) + y2*(pvprn) + *cart - pvshift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n12.get_x() << ", " << n12.get_y() << ", " << -n12.get_z() << ">";

      fp << "   texture_list { texture_" << icol+1 << " texture_" << icol+1 << " texture_" << icol+1 << " } }\n";

      cart++; pviter++; pvpriter++;
      icol++;
    }
  }
  if(colour_vector.size()==0){
    fp << "  texture {\n";
    if(transparent) fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha << "> }\n"; 
    else fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
  }
  fp << "}\n";
}

void draw_fancy_ribbon_pov(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const std::vector<Cartesian> &vertices, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, const double *colour, const std::vector<Cartesian> &colour_vector, int transparent, double alpha, const bool two_colour){

  std::vector<Cartesian>::const_iterator cart;
  std::vector<Cartesian>::const_iterator pviter;
  std::vector<Cartesian>::const_iterator pvpriter;
  
  int idx = vertices.size()/2;
  Cartesian p_t = vertices[idx+1] + vertices[idx-1] - 2*vertices[idx];
  p_t.normalize();
  Cartesian pvpr_t = pvpr[idx];
  pvpr_t.normalize();

  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;

  int nsectors = 180;

  for(unsigned i=0; i<colour_vector.size();i++){
    if(transparent)
      fp << "#declare texture_" << i << " = texture { pigment { color rgbt <" << colour_vector[i].get_x() << ", " << colour_vector[i].get_y() << ", " << colour_vector[i].get_z() << ", " << 1.0-alpha << "> }  finish {diffuse 1.0 specular 1.0} }\n";
    else
      fp << "#declare texture_" << i << " = texture { pigment { color rgb <" << colour_vector[i].get_x() << ", " << colour_vector[i].get_y() << ", " << colour_vector[i].get_z() << "> }  finish {diffuse 1.0 specular 1.0} }\n";
  }
  fp << "#declare texture_grey = texture { pigment { color rgb <0.6, 0.6 , 0.6 > }  finish {diffuse 1.0 specular 1.0} }\n";

  fp << "mesh {\n";

    cart=vertices.begin();
    pviter=pv.begin();
    pvpriter=pvpr.begin();

    int icol = 0;
    while(cart!=vertices.end()-1&&pviter!=pv.end()-1&&pvpriter!=pvpr.end()-1){
      
      double a = pvpriter->length();
      double a2 = (pvpriter+1)->length();
      Cartesian pvshift = (pviter->length()-a)/pviter->length() * *pviter;
      Cartesian pv1shift = ((pviter+1)->length()-a2)/(pviter+1)->length() * *(pviter+1);

      Cartesian p1 = pvshift + 0.7**pvpriter + *cart;
      Cartesian p2 = pvshift - 0.7**pvpriter + *cart;
      Cartesian p3 = pv1shift - 0.7**(pvpriter+1) + *(cart+1);
      Cartesian p4 = pv1shift + 0.7**(pvpriter+1) + *(cart+1);

      Cartesian p5 = -pvshift + 0.7**pvpriter + *cart;
      Cartesian p6 = -pvshift - 0.7**pvpriter + *cart;
      Cartesian p7 = -pv1shift - 0.7**(pvpriter+1) + *(cart+1);
      Cartesian p8 = -pv1shift + 0.7**(pvpriter+1) + *(cart+1);

      Cartesian n1 = *pviter;
      Cartesian n2 = *(pviter+1);
      Cartesian n3 = *pvpriter;
      Cartesian n4 = *(pvpriter+1);
      n1.normalize();
      n2.normalize();
      n3.normalize();
      n4.normalize();
      n1 =  quat.getMatrix()*(objrotmatrix*(n1));
      n2 =  quat.getMatrix()*(objrotmatrix*(n2));
      n3 =  quat.getMatrix()*(objrotmatrix*(n3));
      n4 =  quat.getMatrix()*(objrotmatrix*(n4));

      fp << "smooth_triangle {\n";
      Cartesian p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p5+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p8+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      if(face_one&&two_colour)
        fp << "   texture_list { texture_grey texture_grey texture_grey } }\n";
      else
        fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p4+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p8+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      if(face_one&&two_colour)
        fp << "   texture_list { texture_grey texture_grey texture_grey } }\n";
      else
        fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";
      
      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p6+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p7+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      if(!face_one&&two_colour)
        fp << "   texture_list { texture_grey texture_grey texture_grey } }\n";
      else
        fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p7+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      if(!face_one&&two_colour)
        fp << "   texture_list { texture_grey texture_grey texture_grey } }\n";
      else
        fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      cart++; pviter++; pvpriter++;
      icol++;
    }

  if(colour_vector.size()==0){
    fp << "  texture {\n";
    if(transparent) fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha << "> }\n"; 
    else fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
  }
  fp << "}\n";

  fp << "mesh {\n";
  double theta,theta2;

  for(int i=0;i<360;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    theta2 = (double)(i-360/nsectors)/360.0 * PIBY2;

    cart=vertices.begin();
    pviter=pv.begin();
    pvpriter=pvpr.begin();

    int icol = 0;
    while(cart!=vertices.end()-1&&pviter!=pv.end()-1&&pvpriter!=pvpr.end()-1){

      double x = cos(theta);
      double y = sin(theta);
      double x2 = cos(theta2);
      double y2 = sin(theta2);

      double a = pvpriter->length();
      double b = pvpriter->length();
      double a2 = (pvpriter+1)->length();
      double b2 = (pvpriter+1)->length();

      Cartesian pvn = *pviter;
      Cartesian pvprn = *pvpriter;
      pvn.normalize(a);
      pvprn.normalize(a);
      Cartesian pv1n = *(pviter+1);
      Cartesian pvpr1n = *(pvpriter+1);
      pv1n.normalize(a);
      pvpr1n.normalize(a);

      Cartesian pvshift = (pviter->length()-a)/pviter->length() * *pviter;
      Cartesian pv1shift = ((pviter+1)->length()-a2)/(pviter+1)->length() * *(pviter+1);

      Cartesian n1 = x*(pvn)/(a*a) + y*(pvprn)/(b*b);
      n1 = quat.getMatrix()*(objrotmatrix*(n1));
      Cartesian n2 = x*(pv1n)/(a2*a2) + y*(*(pvpriter+1))/(b2*b2);
      n2 = quat.getMatrix()*(objrotmatrix*(n2));
      Cartesian n12 = x2*(pvn)/(a*a) + y2*(pvprn)/(b*b);
      n12 = quat.getMatrix()*(objrotmatrix*(n12));
      Cartesian n22 = x2*(pv1n)/(a2*a2) + y2*(*(pvpriter+1))/(b2*b2);
      n22 = quat.getMatrix()*(objrotmatrix*(n22));

      fp << "smooth_triangle {\n";
      
      Cartesian p1 = x*(pvn) + y*(pvprn) + *cart + pvshift;
      Cartesian p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x*(pv1n) + y*(pvpr1n+objorigin) + *(cart+1) + pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">,";
      
      p1 = x2*(pv1n) + y2*(pvpr1n+objorigin) + *(cart+1) + pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">";
      
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";
      
      fp << "smooth_triangle {\n";
      
      p1 = x*(pvn) + y*(pvprn) + *cart + pvshift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x2*(pv1n) + y2*(pvpr1n) + *(cart+1) + pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">,";
      
      p1 = x2*(pvn) + y2*(pvprn) + *cart + pvshift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n12.get_x() << ", " << n12.get_y() << ", " << -n12.get_z() << ">";

      fp << "   texture_list { texture_" << icol+1 << " texture_" << icol+1 << " texture_" << icol+1 << " } }\n";

      cart++; pviter++; pvpriter++;
      icol++;
    }
  }
  if(colour_vector.size()==0){
    fp << "  texture {\n";
    if(transparent) fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha << "> }\n"; 
    else fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
  }
  fp << "}\n";
  fp << "mesh {\n";

  for(int i=0;i<360;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    theta2 = (double)(i-360/nsectors)/360.0 * PIBY2;

    cart=vertices.begin();
    pviter=pv.begin();
    pvpriter=pvpr.begin();

    int icol = 0;
    while(cart!=vertices.end()-1&&pviter!=pv.end()-1&&pvpriter!=pvpr.end()-1){

      double x = cos(theta);
      double y = sin(theta);
      double x2 = cos(theta2);
      double y2 = sin(theta2);

      double a = pvpriter->length();
      double b = pvpriter->length();
      double a2 = (pvpriter+1)->length();
      double b2 = (pvpriter+1)->length();

      Cartesian pvn = *pviter;
      Cartesian pvprn = *pvpriter;
      pvn.normalize(a);
      pvprn.normalize(a);
      Cartesian pv1n = *(pviter+1);
      Cartesian pvpr1n = *(pvpriter+1);
      pv1n.normalize(a);
      pvpr1n.normalize(a);

      Cartesian pvshift = (pviter->length()-a)/pviter->length() * *pviter;
      Cartesian pv1shift = ((pviter+1)->length()-a2)/(pviter+1)->length() * *(pviter+1);

      Cartesian n1 = x*(pvn)/(a*a) + y*(pvprn)/(b*b);
      n1 = quat.getMatrix()*(objrotmatrix*(n1));
      Cartesian n2 = x*(pv1n)/(a2*a2) + y*(*(pvpriter+1))/(b2*b2);
      n2 = quat.getMatrix()*(objrotmatrix*(n2));
      Cartesian n12 = x2*(pvn)/(a*a) + y2*(pvprn)/(b*b);
      n12 = quat.getMatrix()*(objrotmatrix*(n12));
      Cartesian n22 = x2*(pv1n)/(a2*a2) + y2*(*(pvpriter+1))/(b2*b2);
      n22 = quat.getMatrix()*(objrotmatrix*(n22));

      fp << "smooth_triangle {\n";
      
      Cartesian p1 = x*(pvn) + y*(pvprn) + *cart - pvshift;
      Cartesian p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x*(pv1n) + y*(pvpr1n+objorigin) + *(cart+1) - pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">,";
      
      p1 = x2*(pv1n) + y2*(pvpr1n+objorigin) + *(cart+1) - pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">";
      
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";
      
      fp << "smooth_triangle {\n";
      
      p1 = x*(pvn) + y*(pvprn) + *cart - pvshift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";

      p1 = x2*(pv1n) + y2*(pvpr1n) + *(cart+1) - pv1shift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">";
      fp << "  < " << n22.get_x() << ", " << n22.get_y() << ", " << -n22.get_z() << ">,";
      
      p1 = x2*(pvn) + y2*(pvprn) + *cart - pvshift;
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n12.get_x() << ", " << n12.get_y() << ", " << -n12.get_z() << ">";

      fp << "   texture_list { texture_" << icol+1 << " texture_" << icol+1 << " texture_" << icol+1 << " } }\n";

      cart++; pviter++; pvpriter++;
      icol++;
    }
  }
  if(colour_vector.size()==0){
    fp << "  texture {\n";
    if(transparent) fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha << "> }\n"; 
    else fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
  }
  fp << "}\n";
}

void draw_flat_ribbon_pov(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const std::vector<Cartesian> &vertices, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, const double *colour, const std::vector<Cartesian> &colour_vector, int transparent, double alpha){

  std::vector<Cartesian>::const_iterator cart;
  std::vector<Cartesian>::const_iterator pviter;
  std::vector<Cartesian>::const_iterator pvpriter;
  
  for(unsigned i=0; i<colour_vector.size();i++){
    if(transparent)
      fp << "#declare texture_" << i << " = texture { pigment { color rgbt <" << colour_vector[i].get_x() << ", " << colour_vector[i].get_y() << ", " << colour_vector[i].get_z() << ", " << 1.0-alpha << "> }  finish {diffuse 1.0 specular 1.0} }\n";
    else
      fp << "#declare texture_" << i << " = texture { pigment { color rgb <" << colour_vector[i].get_x() << ", " << colour_vector[i].get_y() << ", " << colour_vector[i].get_z() << "> }  finish {diffuse 1.0 specular 1.0} }\n";
  }

  fp << "mesh {\n";

    cart=vertices.begin();
    pviter=pv.begin();
    pvpriter=pvpr.begin();

    int icol = 0;
    while(cart!=vertices.end()-1&&pviter!=pv.end()-1&&pvpriter!=pvpr.end()-1){
      
      Cartesian p1 = *pviter + *pvpriter + *cart;
      Cartesian p2 = *pviter - *pvpriter + *cart;
      Cartesian p3 = *(pviter+1) - *(pvpriter+1) + *(cart+1);
      Cartesian p4 = *(pviter+1) + *(pvpriter+1) + *(cart+1);

      Cartesian p5 = -*pviter + *pvpriter + *cart;
      Cartesian p6 = -*pviter - *pvpriter + *cart;
      Cartesian p7 = -*(pviter+1) - *(pvpriter+1) + *(cart+1);
      Cartesian p8 = -*(pviter+1) + *(pvpriter+1) + *(cart+1);

      Cartesian n1 = *pviter;
      Cartesian n2 = *(pviter+1);
      Cartesian n3 = *pvpriter;
      Cartesian n4 = *(pvpriter+1);
      n1.normalize();
      n2.normalize();
      n3.normalize();
      n4.normalize();
      n1 =  quat.getMatrix()*(objrotmatrix*(n1));
      n2 =  quat.getMatrix()*(objrotmatrix*(n2));
      n3 =  quat.getMatrix()*(objrotmatrix*(n3));
      n4 =  quat.getMatrix()*(objrotmatrix*(n4));

      fp << "smooth_triangle {\n";
      Cartesian p = quat.getMatrix()*(objrotmatrix*(p5+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p6+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p7+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">";
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p5+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p7+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p8+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">";
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">";
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n1.get_x() << ", " << n1.get_y() << ", " << -n1.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p4+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n2.get_x() << ", " << n2.get_y() << ", " << -n2.get_z() << ">";
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p5+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p8+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p4+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p8+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";
      
      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p6+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p7+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      fp << "smooth_triangle {\n";
      p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n3.get_x() << ", " << n3.get_y() << ", " << -n3.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">,";
      p = quat.getMatrix()*(objrotmatrix*(p7+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,";
      fp << "  < " << n4.get_x() << ", " << n4.get_y() << ", " << -n4.get_z() << ">";
      fp << "   texture_list { texture_" << icol << " texture_" << icol << " texture_" << icol << " } }\n";

      cart++; pviter++; pvpriter++;
      icol++;
    }

  if(colour_vector.size()==0){
    fp << "  texture {\n";
    if(transparent) fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << ", " << 1.0-alpha << "> }\n"; 
    else fp << "    pigment { color rgb <" << colour[0] << ", " << colour[1] << ", " << colour[2] << "> }\n"; 
    fp << "    finish {diffuse 1.0 specular 1.0}\n";
    fp << "  }\n";
  }
  fp << "}\n";
}

void Ribbon::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
  if(vertices.size()<2) return;
  /*
  // More complicated
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }
  */

  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);

  if(style==0)
    draw_elliptical_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices,AB[0],AB[1],colour,colour_vector,get_transparent(),alpha);
  else if(style==1)
    draw_flat_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices,AB[0],AB[1],colour,colour_vector,get_transparent(),alpha);
  else if(style==2)
    draw_flat_round_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices,AB[0],AB[1],colour,colour_vector,get_transparent(),alpha,two_colour);
  else if(style==3)
    draw_fancy_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices,AB[0],AB[1],colour,colour_vector,get_transparent(),alpha,two_colour);
  else
    draw_flat_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices,AB[0],AB[1],colour,colour_vector,get_transparent(),alpha);
}

void ArrowHeadRibbon::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
  if(vertices.size()<2) return;
  /*
  // More complicated
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }
  */

  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);

  std::vector<Cartesian> vertices_tmp = vertices;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;

  if(insert_at_vertex>0){
    vertices_tmp.insert(vertices_tmp.begin()+insert_at_vertex,0.9*vertices_tmp[insert_at_vertex-1] + 0.1*vertices_tmp[insert_at_vertex]);
    if(colour_vector_tmp.size()>0) colour_vector_tmp.insert(colour_vector_tmp.begin()+insert_at_vertex,colour_vector_tmp[insert_at_vertex-1]);
  }


  if(style==0)
    draw_elliptical_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices_tmp,AB[0],AB[1],colour,colour_vector_tmp,get_transparent(),alpha);
  else if(style==1)
    draw_flat_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices_tmp,AB[0],AB[1],colour,colour_vector_tmp,get_transparent(),alpha);
  else if(style==2)
    draw_flat_round_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices_tmp,AB[0],AB[1],colour,colour_vector_tmp,get_transparent(),alpha,two_colour);
  else if(style==3)
    draw_fancy_ribbon_pov(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,vertices_tmp,AB[0],AB[1],colour,colour_vector_tmp,get_transparent(),alpha,two_colour);
}

void BillBoard::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
   if(!image.get_width()||!image.get_height()||!image.get_pixels())
      return;
  fp << "image={}\n";
  fp << "image[\"position\"]=[" << vertices[0].get_x() << "," << vertices[0].get_y() << "]\n";
  fp << "image[\"scale_w\"]=" << scale_w << "\n";
  fp << "image[\"scale_h\"]=" << scale_h << "\n";
  
  /*
  char *mkstemps_template = new char[256]; //20 perhaps?
  mkstemps_template[0] = '\0';
  strncat(mkstemps_template,"/tmp/tmpXXXXXXX.png",19);
  int fildes = mkstemps(mkstemps_template,4);
  */
  fp << "image[\"filename\"]=\"" << filename << "\"\n";
  fp << "image_labels.append(image)\n";
  
  /*
  delete mkstemps_template;

  image_info image_tmp = image;
  image_tmp.invert();
  image_tmp.writepng(fildes);
  */
}

void Text::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  // Needs to be rewitten with Qt font stuff
  /*
  GLint curr_viewport[4];
  glGetIntegerv(GL_VIEWPORT,curr_viewport);
  SetRasterPosition(vertices[0].get_x(),vertices[0].get_y(),vertices[0].get_z());
  GLboolean valid=GL_TRUE;
  glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID,&valid);
  GLfloat params[4];
  glGetFloatv(GL_CURRENT_RASTER_POSITION,params);
  int font_id = LoadFont(); 
  MGFontInfo finfo = FontCache::GetFont(font_id);
  if(valid&&params[0]/curr_viewport[2]>=0.0&&params[0]/curr_viewport[2]<1.0&&params[1]/curr_viewport[3]>=0.0&&params[1]/curr_viewport[3]<1.0){
  //std::cout << StripTags() << " " << vertices[0] << "\n";
  //std::cout << params[0]/curr_viewport[2] << "," << params[1]/curr_viewport[3] << "\n";
  fp << "label={}\n";
  fp << "label[\"position\"]=[" << params[0]/curr_viewport[2] << "," << params[1]/curr_viewport[3] << "]\n";
  fp << "font ={}\n";
  fp << "font[\"family\"]=\"" << finfo.Family() << "\"\n";
  fp << "font[\"weight\"]=\"" << finfo.Weight() << "\"\n";
  fp << "font[\"slant\"]=\"" << finfo.Slant() << "\"\n";
  fp << "font[\"size\"]=\"" << finfo.Size() << "\"\n";
  fp << "label[\"font\"]=font\n";
  fp << "label[\"text\"]=\"" << StripTags() << "\"\n";  // One day we will be able to handle tags outside GL
  glGetFloatv(GL_COLOR_CLEAR_VALUE,params);
  double y = params[0]*0.299 + params[1]*0.587 + params[2]*0.114;
  if(y>0.5)
    fp << "label[\"colour\"]=[0,0,0]\n";
  else
    fp << "label[\"colour\"]=[1,1,1]\n";
  fp << "text_labels.append(label)\n";
  }
  */
}

void BillBoardText::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  /*
  int font_id = LoadFont(); 
  MGFontInfo finfo = FontCache::GetFont(font_id);
  fp << "label={}\n";
  fp << "label[\"position\"]=[" << vertices[0].get_x() << "," << vertices[0].get_y() << "]\n";
  fp << "font ={}\n";
  fp << "font[\"family\"]=\"" << finfo.Family() << "\"\n";
  fp << "font[\"weight\"]=\"" << finfo.Weight() << "\"\n";
  fp << "font[\"slant\"]=\"" << finfo.Slant() << "\"\n";
  fp << "font[\"size\"]=\"" << finfo.Size() << "\"\n";
  fp << "label[\"font\"]=font\n";
  fp << "label[\"text\"]=\"" << StripTags() << "\"\n";  // One day we will be able to handle tags outside GL
  GLfloat params[4];
  glGetFloatv(GL_COLOR_CLEAR_VALUE,params);
  double y = params[0]*0.299 + params[1]*0.587 + params[2]*0.114;
  if(y>0.5)
    fp << "label[\"colour\"]=[0,0,0]\n";
  else
    fp << "label[\"colour\"]=[1,1,1]\n";
  fp << "text_labels.append(label)\n";
  */
}

void Arrow::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  ArrowElement::DrawPovray(fp, quat, radius, ox, oy, oz, objrotmatrix, objorigin,v);
}

void DashArrow::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  DashArrowElement::DrawPovray(fp, quat, radius, ox, oy, oz, objrotmatrix, objorigin,v);
}


void ArrowElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  Cartesian v1 = vertices[0];
  Cartesian v2 = 0.2*vertices[0]+0.8*vertices[1];
  Cartesian v3 = vertices[1];
  fp << "union {\n";
  fp << "cylinder {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/10.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
  fp << "cone {\n";
  p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">," << GetSize()/5.0 << "\n";
  p = quat.getMatrix()*(objrotmatrix*(v3+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << 0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
  fp << "}\n";
}

void DashArrowElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  Cartesian v1 = vertices[0];
  Cartesian v2 = 0.2*vertices[0]+0.8*vertices[1];
  Cartesian v3 = vertices[1];
  Cartesian p;
  double length = (v1-v2).length();
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2)==0)nsegments++;
  for(int i=0;i<nsegments;i++){
    if(i%2 >0){
      double frac = (double)i/(double)nsegments;
      double frac2 =(double)(i+1)/(double)nsegments;
      fp << "cylinder {\n";
      Cartesian v3 = frac*v1 + (1-frac)*v2;
      Cartesian v4 = frac2*v1 + (1-frac2)*v2;
      p = quat.getMatrix()*(objrotmatrix*(v3+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      p = quat.getMatrix()*(objrotmatrix*(v4+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      fp << "  " << GetSize()/10.0 << "\n";
      fp << "  texture {\n";
      OutputPovrayPigment(fp); 
      fp << "    finish {diffuse 1.0 specular 1.0}\n";
      fp << "  }\n";
      fp << "}\n";
    }
  }
  fp << "cone {\n";
  p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">," << GetSize()/5.0 << "\n";
  p = quat.getMatrix()*(objrotmatrix*(v3+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << 0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void Point::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
}

void Line::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v) {

  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  fp << GetSize() << " slw ";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
  Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);

  fp << pps.get_x() << " " << pps.get_y() << " l ";
  fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
}

void LineStrip::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
   std::cout << "LineStrip::DrawPostScript not yet implemented\n";
}

void DashLineElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  double length = sqrt((vertices[1].get_x() - vertices[0].get_x()) * (vertices[1].get_x() - vertices[0].get_x()) +
		       (vertices[1].get_y() - vertices[0].get_y()) * (vertices[1].get_y() - vertices[0].get_y()) +
		       (vertices[1].get_z() - vertices[0].get_z()) * (vertices[1].get_z() - vertices[0].get_z()));
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2)==0)nsegments++;

  double frac = 0;
  double frac2 = 0.5/(double)nsegments;
  Cartesian v1 = frac*vertices[1]+(1-frac)*vertices[0];
  Cartesian v2 = frac2*vertices[1]+(1-frac2)*vertices[0];
  fp << "cylinder {\n";
  Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/50.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  for(int i=1;i<nsegments;i++){
    if(i%2 >0){
      double frac = (i+0.5)/(double)nsegments;
      double frac2 =(i+1.5)/(double)nsegments;
      if(frac2>1.0) frac2 = 1.0;
      Cartesian v1 = frac*vertices[1]+(1-frac)*vertices[0];
      Cartesian v2 = frac2*vertices[1]+(1-frac2)*vertices[0];

      fp << "cylinder {\n";
      Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      fp << "  " << GetSize()/50.0 << "\n";
      fp << "  texture {\n";
      OutputPovrayPigment(fp); 
      fp << "    finish {diffuse 1.0 specular 1.0}\n";
      fp << "  }\n";
      fp << "}\n";
    }
  }

  fp << "sphere {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/50.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";

  fp << "sphere {\n";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
  fp << "  " << GetSize()/50.0 << "\n";
  fp << "  texture {\n";
  OutputPovrayPigment(fp); 
  fp << "    finish {diffuse 1.0 specular 1.0}\n";
  fp << "  }\n";
  fp << "}\n";
}

void DashLineElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  double length = sqrt((vertices[1].get_x() - vertices[0].get_x()) * (vertices[1].get_x() - vertices[0].get_x()) +
		       (vertices[1].get_y() - vertices[0].get_y()) * (vertices[1].get_y() - vertices[0].get_y()) +
		       (vertices[1].get_z() - vertices[0].get_z()) * (vertices[1].get_z() - vertices[0].get_z()));
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2)==0)nsegments++;

  fp << GetSize() << " slw ";
  double frac = 0;
  double frac2 =0.5/(double)nsegments;
  Cartesian v1 = frac*vertices[1]+(1-frac)*vertices[0];
  Cartesian v2 = frac2*vertices[1]+(1-frac2)*vertices[0];
  Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
  Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
  Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
  p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
  pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";
  fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";

  for(int i=1;i<nsegments;i++){
    if(i%2 >0){
      fp << GetSize() << " slw ";
      frac = (i+0.5)/(double)nsegments;
      frac2 =(i+1.5)/(double)nsegments;
      if(frac2>1.0) frac2 = 1.0;
      Cartesian v1 = frac*vertices[1]+(1-frac)*vertices[0];
      Cartesian v2 = frac2*vertices[1]+(1-frac2)*vertices[0];
      Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
      Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
      Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
      p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
      pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    }
  }

}

void DashLine::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  double length = sqrt((vertices[1].get_x() - vertices[0].get_x()) * (vertices[1].get_x() - vertices[0].get_x()) +
		       (vertices[1].get_y() - vertices[0].get_y()) * (vertices[1].get_y() - vertices[0].get_y()) +
		       (vertices[1].get_z() - vertices[0].get_z()) * (vertices[1].get_z() - vertices[0].get_z()));
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2)==0)nsegments++;

  fp << GetSize() << " slw ";
  double frac = 0;
  double frac2 =0.5/(double)nsegments;
  Cartesian v1 = frac*vertices[1]+(1-frac)*vertices[0];
  Cartesian v2 = frac2*vertices[1]+(1-frac2)*vertices[0];
  Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
  Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
  Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
  p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
  pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";
  fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";

  for(int i=1;i<nsegments;i++){
    if(i%2 >0){
      fp << GetSize() << " slw ";
      frac = (i+0.5)/(double)nsegments;
      frac2 =(i+1.5)/(double)nsegments;
      if(frac2>1.0) frac2 = 1.0;
      Cartesian v1 = frac*vertices[1]+(1-frac)*vertices[0];
      Cartesian v2 = frac2*vertices[1]+(1-frac2)*vertices[0];
      Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
      Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
      Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
      p = quat.getInvMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
      pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    }
  }

}

void Arrow::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  ArrowElement::DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
}

void ArrowElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  fp << GetSize() << " slw ";
  Cartesian vec = vertices[1] - vertices[0];
  Cartesian atmp,btmp;
  Cartesian xaxis(1,0,0);
  Cartesian yaxis(0,1,0);
  Cartesian zaxis(0,0,1);
  atmp = vec.CrossProduct(vec,zaxis);
  if(atmp.length() < 0.000000001){
    atmp = vec.CrossProduct(vec,yaxis);
    if(atmp.length() < 0.000000001){
      atmp = vec.CrossProduct(vec,xaxis);
    }
  }
  btmp = vec.CrossProduct(vec,atmp);

  Cartesian p = vertices[0] + 0.9*vec;
  btmp.normalize((vertices[1]-p).length()/2.0);
  Cartesian p1 = p + btmp;
  Cartesian p2 = p - btmp;
  atmp.normalize((vertices[1]-p).length()/2.0);
  Cartesian p3 = p + atmp;
  Cartesian p4 = p - atmp;
  
  p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
  Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
  p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";
  fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";

  if (arrow_head != 1) {
    //   arrow_head = 0 (end) or 2 (both) 
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p4+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
  }
  if (arrow_head != 0) {
    //   arrow_head = 1 (start) or 2 (both)
    p1=p1-0.8*vec; 
    p2=p2-0.8*vec; 
    p3=p3-0.8*vec; 
    p4=p4-0.8*vec; 
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p4+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
  }

}

void DashArrow::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  DashArrowElement::DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
}

void DashArrowElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  double length = sqrt((vertices[1].get_x() - vertices[0].get_x()) * (vertices[1].get_x() - vertices[0].get_x()) +
		       (vertices[1].get_y() - vertices[0].get_y()) * (vertices[1].get_y() - vertices[0].get_y()) +
		       (vertices[1].get_z() - vertices[0].get_z()) * (vertices[1].get_z() - vertices[0].get_z()));
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2)==0)nsegments++;

  for(int i=0;i<nsegments;i++){
    if(i%2 >0){
      fp << GetSize() << " slw ";
      double frac = (double)i/(double)nsegments;
      double frac2 =(double)(i+1)/(double)nsegments;
      Cartesian v1 = frac*vertices[0]+(1-frac)*vertices[1];
      Cartesian v2 = frac2*vertices[0]+(1-frac2)*vertices[1];
      Cartesian p = quat.getMatrix()*(objrotmatrix*(v1+objorigin)+Cartesian(ox,oy,oz));
      Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
      Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
      p = quat.getMatrix()*(objrotmatrix*(v2+objorigin)+Cartesian(ox,oy,oz));
      pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    }
  }

  Cartesian vec = vertices[1] - vertices[0];
  Cartesian atmp,btmp;
  Cartesian xaxis(1,0,0);
  Cartesian yaxis(0,1,0);
  Cartesian zaxis(0,0,1);
  atmp = vec.CrossProduct(vec,zaxis);
  if(atmp.length() < 0.000000001){
    atmp = vec.CrossProduct(vec,yaxis);
    if(atmp.length() < 0.000000001){
      atmp = vec.CrossProduct(vec,xaxis);
    }
  }
  btmp = vec.CrossProduct(vec,atmp);

  Cartesian p = vertices[0] + 0.9*vec;
  btmp.normalize((vertices[1]-p).length()/2.0);
  Cartesian p1 = p + btmp;
  Cartesian p2 = p - btmp;
  atmp.normalize((vertices[1]-p).length()/2.0);
  Cartesian p3 = p + atmp;
  Cartesian p4 = p - atmp;
  
  if (arrow_head != 1) {
    //   arrow_head = 0 (end) or 2 (both) 
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
    Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p4+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
  }
  if (arrow_head != 0) {
    //   arrow_head = 1 (start) or 2 (both)
    p1=p1-0.8*vec; 
    p2=p2-0.8*vec; 
    p3=p3-0.8*vec; 
    p4=p4-0.8*vec; 
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p1+objorigin)+Cartesian(ox,oy,oz));
    Cartesian pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p2+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p3+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
    fp << GetSize() << " slw ";
    p = quat.getMatrix()*(objrotmatrix*(p4+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
    p = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
    pfrac = Cartesian(p.get_x()/xscale,p.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";
    fp << " cp " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s %ZVALUE: " << p.get_z() << "\n";
  }

}

void CircleElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
   std::cout << "Circle::DrawPostscript not yet implemented\n";
}

void PolygonElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
   std::cout << "PolygonElement::DrawPostscript not yet implemented\n";
}

void FlatSubdivideCartsNormal(const std::vector<Cartesian> &v, std::vector<std::vector<Cartesian> > &v_sub, const double length){
  if(v.size()!=3) return;
  Cartesian p1;
  Cartesian p2;
  Cartesian p3;
  std::vector<Cartesian> t1;
  std::vector<Cartesian> t2;
  std::vector<Cartesian> t3;
  std::vector<Cartesian> t4;
  p1 = Cartesian::MidPoint(v[0],v[1]);
  p2 = Cartesian::MidPoint(v[1],v[2]);
  p3 = Cartesian::MidPoint(v[2],v[0]);
  t1.push_back(v[0]);
  t1.push_back(p1);
  t1.push_back(p3);

  t2.push_back(v[1]);
  t2.push_back(p2);
  t2.push_back(p1);

  t3.push_back(v[2]);
  t3.push_back(p3);
  t3.push_back(p2);

  t4.push_back(p1);
  t4.push_back(p2);
  t4.push_back(p3);

  if(LineLength(p1,p2)<=length) {
   v_sub.push_back(t1);
   v_sub.push_back(t2);
   v_sub.push_back(t3);
   v_sub.push_back(t4);
   return;
  }

  FlatSubdivideCartsNormal(t1,v_sub,length);
  FlatSubdivideCartsNormal(t2,v_sub,length);
  FlatSubdivideCartsNormal(t3,v_sub,length);
  FlatSubdivideCartsNormal(t4,v_sub,length);
}

void FlatSubdivideCartsEdge(const std::vector<Cartesian> &v, std::vector<std::vector<Cartesian> > &v_sub, std::vector<std::vector<Cartesian> > &v_sub_edge, std::vector<std::vector<Cartesian> > &v_sub_corner, const double length){
  if(v.size()!=3) return;
  Cartesian p1;
  Cartesian p2;
  Cartesian p3;
  std::vector<Cartesian> t1;
  std::vector<Cartesian> t2;
  std::vector<Cartesian> t3;
  std::vector<Cartesian> t4;
  p1 = Cartesian::MidPoint(v[0],v[1]);
  p2 = Cartesian::MidPoint(v[1],v[2]);
  p3 = Cartesian::MidPoint(v[2],v[0]);
  t1.push_back(v[0]);
  t1.push_back(p1);
  t1.push_back(p3);

  t2.push_back(p1);
  t2.push_back(v[1]);
  t2.push_back(p2);

  t3.push_back(v[2]);
  t3.push_back(p3);
  t3.push_back(p2);

  t4.push_back(p1);
  t4.push_back(p2);
  t4.push_back(p3);

  if(LineLength(p1,p2)<=length) {
   v_sub_edge.push_back(t1);
   v_sub_edge.push_back(t2);
   v_sub.push_back(t3);
   v_sub.push_back(t4);
   return;
  }

  FlatSubdivideCartsEdge(t1,v_sub,v_sub_edge,v_sub_corner,length);
  FlatSubdivideCartsEdge(t2,v_sub,v_sub_edge,v_sub_corner,length);
  FlatSubdivideCartsNormal(t3,v_sub,length);
  FlatSubdivideCartsNormal(t4,v_sub,length);
}

void FlatSubdivideCartsCorner(const std::vector<Cartesian> &v, std::vector<std::vector<Cartesian> > &v_sub, std::vector<std::vector<Cartesian> > &v_sub_edge, std::vector<std::vector<Cartesian> > &v_sub_corner, const double length){
  if(v.size()!=3) return;
  Cartesian p1;
  Cartesian p2;
  Cartesian p3;
  std::vector<Cartesian> t1;
  std::vector<Cartesian> t2;
  std::vector<Cartesian> t3;
  std::vector<Cartesian> t4;
  p1 = Cartesian::MidPoint(v[0],v[1]);
  p2 = Cartesian::MidPoint(v[1],v[2]);
  p3 = Cartesian::MidPoint(v[2],v[0]);
  t1.push_back(v[0]);
  t1.push_back(p1);
  t1.push_back(p3);

  t2.push_back(p1);
  t2.push_back(v[1]);
  t2.push_back(p2);

  t3.push_back(v[2]);
  t3.push_back(p3);
  t3.push_back(p2);

  t4.push_back(p1);
  t4.push_back(p2);
  t4.push_back(p3);

  if(LineLength(p1,p2)<=length) {
   v_sub_corner.push_back(t1);
   v_sub_edge.push_back(t2);
   v_sub_edge.push_back(t3);
   v_sub.push_back(t4);
   return;
  }

  FlatSubdivideCartsCorner(t1,v_sub,v_sub_edge,v_sub_corner,length);
  FlatSubdivideCartsEdge(t2,v_sub,v_sub_edge,v_sub_corner,length);
  FlatSubdivideCartsEdge(t3,v_sub,v_sub_edge,v_sub_corner,length);
  FlatSubdivideCartsNormal(t4,v_sub,length);
}

void FlatSubdivideCarts(const std::vector<Cartesian> &v, std::vector<std::vector<Cartesian> > &v_sub, std::vector<std::vector<Cartesian> > &v_sub_edge, std::vector<std::vector<Cartesian> > &v_sub_corner, const double length){
  if(v.size()!=3) return;
  Cartesian p1;
  Cartesian p2;
  Cartesian p3;
  std::vector<Cartesian> t1;
  std::vector<Cartesian> t2;
  std::vector<Cartesian> t3;
  std::vector<Cartesian> t4;
  p1 = Cartesian::MidPoint(v[0],v[1]);
  p2 = Cartesian::MidPoint(v[1],v[2]);
  p3 = Cartesian::MidPoint(v[2],v[0]);
  t1.push_back(v[0]);
  t1.push_back(p1);
  t1.push_back(p3);

  t2.push_back(v[1]);
  t2.push_back(p2);
  t2.push_back(p1);

  t3.push_back(v[2]);
  t3.push_back(p3);
  t3.push_back(p2);

  t4.push_back(p1);
  t4.push_back(p2);
  t4.push_back(p3);

  if(LineLength(p1,p2)<=length) {
   v_sub_corner.push_back(t1);
   v_sub_corner.push_back(t2);
   v_sub_corner.push_back(t3);
   v_sub.push_back(t4);
   return;
  }

  FlatSubdivideCartsCorner(t1,v_sub,v_sub_edge,v_sub_corner,length);
  FlatSubdivideCartsCorner(t2,v_sub,v_sub_edge,v_sub_corner,length);
  FlatSubdivideCartsCorner(t3,v_sub,v_sub_edge,v_sub_corner,length);
  FlatSubdivideCartsNormal(t4,v_sub,length);
}

void TriangleElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[2]+objorigin))) {
    return;
  }

  Cartesian p1 = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  Cartesian p2 = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  Cartesian p3 = quat.getMatrix()*(objrotmatrix*(vertices[2]+objorigin)+Cartesian(ox,oy,oz));

  double y = colour[0]*0.299 + colour[1]*0.587 + colour[2]*0.114;

  bool split = false;
  // We need an option to use this somehow.
  //if(LineLength(p1,p2)>0.1) split = true;

  if(split){
    std::vector<Cartesian> vertices_new;
    vertices_new.push_back(p1);
    vertices_new.push_back(p2);
    vertices_new.push_back(p3);
    std::vector<std::vector<Cartesian> > subdivided_carts;
    std::vector<std::vector<Cartesian> > subdivided_carts_edge;
    std::vector<std::vector<Cartesian> > subdivided_carts_corner;
    FlatSubdivideCarts(vertices_new,subdivided_carts,subdivided_carts_edge,subdivided_carts_corner,0.1);
    for(unsigned i=0;i<subdivided_carts.size();i++){
      fp << "gs " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s ";
      Cartesian pfrac = Cartesian(subdivided_carts[i][0].get_x()/xscale,subdivided_carts[i][0].get_y()/yscale,0);
      Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
      pfrac = Cartesian(subdivided_carts[i][1].get_x()/xscale,subdivided_carts[i][1].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      pfrac = Cartesian(subdivided_carts[i][2].get_x()/xscale,subdivided_carts[i][2].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      fp << " fill cp gr %ZVALUE: " << ((subdivided_carts[i][0]+subdivided_carts[i][1]+subdivided_carts[i][2])/3.0).get_z() << "\n";
    }
    for(unsigned i=0;i<subdivided_carts_edge.size();i++){
      fp << "gs " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s ";
      Cartesian pfrac = Cartesian(subdivided_carts_edge[i][0].get_x()/xscale,subdivided_carts_edge[i][0].get_y()/yscale,0);
      Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
      pfrac = Cartesian(subdivided_carts_edge[i][1].get_x()/xscale,subdivided_carts_edge[i][1].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      pfrac = Cartesian(subdivided_carts_edge[i][2].get_x()/xscale,subdivided_carts_edge[i][2].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      fp << " fill cp gr %ZVALUE: " << ((subdivided_carts_edge[i][0]+subdivided_carts_edge[i][1]+subdivided_carts_edge[i][2])/3.0).get_z() << "\n";

      fp << " 0.10 slw ";
      pfrac = Cartesian(subdivided_carts_edge[i][0].get_x()/xscale,subdivided_carts_edge[i][0].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

      pfrac = Cartesian(subdivided_carts_edge[i][1].get_x()/xscale,subdivided_carts_edge[i][1].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";

      if(y>0.1){
        fp << " 0 0 0 srgb s %ZVALUE: " << ((subdivided_carts_edge[i][0]+subdivided_carts_edge[i][1])/2.0).get_z()+0.015 << "\n";
      }else{
        fp << " 1 1 1 srgb s %ZVALUE: " << ((subdivided_carts_edge[i][0]+subdivided_carts_edge[i][1])/2.0).get_z()+0.015 << "\n";
      }
    }
    for(unsigned i=0;i<subdivided_carts_corner.size();i++){
      fp << "gs " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s ";
      Cartesian pfrac = Cartesian(subdivided_carts_corner[i][0].get_x()/xscale,subdivided_carts_corner[i][0].get_y()/yscale,0);
      Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
      pfrac = Cartesian(subdivided_carts_corner[i][1].get_x()/xscale,subdivided_carts_corner[i][1].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      pfrac = Cartesian(subdivided_carts_corner[i][2].get_x()/xscale,subdivided_carts_corner[i][2].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";
      fp << " fill cp gr %ZVALUE: " << ((subdivided_carts_corner[i][0]+subdivided_carts_corner[i][1]+subdivided_carts_corner[i][2])/3.0).get_z() << "\n";

      fp << " 0.10 slw ";
      pfrac = Cartesian(subdivided_carts_corner[i][2].get_x()/xscale,subdivided_carts_corner[i][2].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

      pfrac = Cartesian(subdivided_carts_corner[i][0].get_x()/xscale,subdivided_carts_corner[i][0].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";

      pfrac = Cartesian(subdivided_carts_corner[i][1].get_x()/xscale,subdivided_carts_corner[i][1].get_y()/yscale,0);
      pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
      fp << pps.get_x() << " " << pps.get_y() << " l ";

      if(y>0.1){
        fp << " 0 0 0 srgb s %ZVALUE: " << ((subdivided_carts_corner[i][0]+subdivided_carts_corner[i][1]+subdivided_carts_corner[i][2])/3.0).get_z()+0.015 << "\n";
      }else{
        fp << " 1 1 1 srgb s %ZVALUE: " << ((subdivided_carts_corner[i][0]+subdivided_carts_corner[i][1]+subdivided_carts_corner[i][2])/3.0).get_z()+0.015 << "\n";
      }

    }
  } else {

    fp << "gs " << colour[0] << " " << colour[1] << " " << colour[2] << " srgb s ";
    Cartesian pfrac = Cartesian(p1.get_x()/xscale,p1.get_y()/yscale,0);
    Cartesian pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

    pfrac = Cartesian(p2.get_x()/xscale,p2.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";

    pfrac = Cartesian(p3.get_x()/xscale,p3.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";

    fp << " fill cp gr %ZVALUE: " << ((p1+p2+p3)/3.0).get_z() << "\n";


    // Optional outline. We'll have to think of a way to trigger this selectively.
    fp << " 0.10 slw ";
    pfrac = Cartesian(p1.get_x()/xscale,p1.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

    pfrac = Cartesian(p2.get_x()/xscale,p2.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";

    pfrac = Cartesian(p3.get_x()/xscale,p3.get_y()/yscale,0);
    pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
    fp << pps.get_x() << " " << pps.get_y() << " l ";

    if(y>0.1){
      fp << " cp 0 0 0 srgb s %ZVALUE: " << ((p1+p2+p3)/3.0).get_z()+0.015 << "\n";
    }else{
      fp << " cp 1 1 1 srgb s %ZVALUE: " << ((p1+p2+p3)/3.0).get_z()+0.015 << "\n";
    }
  }
}

void QuadElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[2]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[3]+objorigin))) {
    return;
  }

  double y = colour[0]*0.299 + colour[1]*0.587 + colour[2]*0.114;

  Cartesian p1 = quat.getMatrix()*(objrotmatrix*(vertices[0]+objorigin)+Cartesian(ox,oy,oz));
  Cartesian p2 = quat.getMatrix()*(objrotmatrix*(vertices[1]+objorigin)+Cartesian(ox,oy,oz));
  Cartesian p3 = quat.getMatrix()*(objrotmatrix*(vertices[2]+objorigin)+Cartesian(ox,oy,oz));
  Cartesian p4 = quat.getMatrix()*(objrotmatrix*(vertices[3]+objorigin)+Cartesian(ox,oy,oz));

  bool split_x = false;
  bool split_y = false;
  //// We need an option to use this somehow.
  //if(LineLength(p1,p2)>0.1) split_x = true;
  //if(LineLength(p1,p4)>0.1) split_y = true;

  Cartesian pfrac;
  Cartesian pps;

  if(!split_x&&!split_y){
  fp << "gs "<< colour[0] << " " << colour[1] << " " << colour[2] << " srgb s ";
  pfrac = Cartesian(p1.get_x()/xscale,p1.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
  pfrac = Cartesian(p2.get_x()/xscale,p2.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";
  pfrac = Cartesian(p3.get_x()/xscale,p3.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";
  pfrac = Cartesian(p4.get_x()/xscale,p4.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";
  fp << " fill cp gr %ZVALUE: " << ((p1+p2+p3+p4)/4.0).get_z() << "\n";
  } else {
    int nsteps_x = int(LineLength(p1,p2)/0.1)+1;
    int nsteps_y = int(LineLength(p1,p4)/0.1)+1;
    Cartesian p1p, p2p, p3p, p4p;
    Cartesian p1p_1, p2p_1, p3p_1, p4p_1;
    Cartesian p1p_2, p2p_2, p3p_2, p4p_2;
    for(int i=0;i<nsteps_x;i++){
      double frac_x = (double)i/(double)nsteps_x;
      double frac2_x = (double)(i+1)/(double)nsteps_x;
      if(frac2_x>1.0) frac2_x = 1.0;
      for(int j=0;j<nsteps_y;j++){
        double frac_y = (double)j/(double)nsteps_y;
        double frac2_y = (double)(j+1)/(double)nsteps_y;
        if(frac2_y>1.0) frac2_y = 1.0;

        p1p_1 = frac_x*p2 + (1-frac_x)*p1;
        p1p_2 = frac_x*p3 + (1-frac_x)*p4;
	p1p = frac_y*p1p_2 + (1-frac_y)*p1p_1;

        p2p_1 = frac2_x*p2 + (1-frac2_x)*p1;
        p2p_2 = frac2_x*p3 + (1-frac2_x)*p4;
	p2p = frac_y*p2p_2 + (1-frac_y)*p2p_1;

        p3p_1 = frac2_x*p2 + (1-frac2_x)*p1;
        p3p_2 = frac2_x*p3 + (1-frac2_x)*p4;
	p3p = frac2_y*p3p_2 + (1-frac2_y)*p3p_1;

        p4p_1 = frac_x*p2 + (1-frac_x)*p1;
        p4p_2 = frac_x*p3 + (1-frac_x)*p4;
	p4p = frac2_y*p4p_2 + (1-frac2_y)*p4p_1;

        fp << "gs "<< colour[0] << " " << colour[1] << " " << colour[2] << " srgb s ";
        pfrac = Cartesian(p1p.get_x()/xscale,p1p.get_y()/yscale,0);
        pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
        fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";
        pfrac = Cartesian(p2p.get_x()/xscale,p2p.get_y()/yscale,0);
        pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
        fp << pps.get_x() << " " << pps.get_y() << " l ";
        pfrac = Cartesian(p3p.get_x()/xscale,p3p.get_y()/yscale,0);
        pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
        fp << pps.get_x() << " " << pps.get_y() << " l ";
        pfrac = Cartesian(p4p.get_x()/xscale,p4p.get_y()/yscale,0);
        pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
        fp << pps.get_x() << " " << pps.get_y() << " l ";
        fp << " fill cp gr %ZVALUE: " << ((p1p+p2p+p3p+p4p)/4.0).get_z() << "\n";

	if(j==0){
          fp << " 0.20 slw ";
          pfrac = Cartesian(p1p.get_x()/xscale,p1p.get_y()/yscale,0);
          pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
          fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

          pfrac = Cartesian(p2p.get_x()/xscale,p2p.get_y()/yscale,0);
          pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
          fp << pps.get_x() << " " << pps.get_y() << " l ";

          if(y>0.1){
            fp << " cp 0 0 0 srgb s %ZVALUE: " << ((p1p+p2p)/2.0).get_z()+0.015 << "\n";
          }else{
            fp << " cp 1 1 1 srgb s %ZVALUE: " << ((p1p+p2p)/2.0).get_z()+0.015 << "\n";
	  }
        }

	if(j==nsteps_y-1){
          fp << " 0.20 slw gs ";
          pfrac = Cartesian(p3p.get_x()/xscale,p3p.get_y()/yscale,0);
          pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
          fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

          pfrac = Cartesian(p4p.get_x()/xscale,p4p.get_y()/yscale,0);
          pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
          fp << pps.get_x() << " " << pps.get_y() << " l ";

          if(y>0.1){
            fp << " cp 0 0 0 srgb s gr %ZVALUE: " << ((p3p+p4p)/2.0).get_z()+0.015 << "\n";
          }else{
            fp << " cp 1 1 1 srgb s gr %ZVALUE: " << ((p3p+p4p)/2.0).get_z()+0.015 << "\n";
	  }
        }

	if(i==0){
          fp << " 0.20 slw gs ";
          pfrac = Cartesian(p1p.get_x()/xscale,p1p.get_y()/yscale,0);
          pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
          fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

          pfrac = Cartesian(p4p.get_x()/xscale,p4p.get_y()/yscale,0);
          pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
          fp << pps.get_x() << " " << pps.get_y() << " l ";

          if(y>0.1){
            fp << " cp 0 0 0 srgb s gr %ZVALUE: " << ((p1p+p4p)/2.0).get_z()+0.015 << "\n";
          }else{
            fp << " cp 1 1 1 srgb s gr %ZVALUE: " << ((p1p+p4p)/2.0).get_z()+0.015 << "\n";
	  }
        }

	if(i==nsteps_x-1){
          fp << " 0.20 slw ";
          pfrac = Cartesian(p3p.get_x()/xscale,p3p.get_y()/yscale,0);
          pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
          fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

          pfrac = Cartesian(p2p.get_x()/xscale,p2p.get_y()/yscale,0);
          pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
          fp << pps.get_x() << " " << pps.get_y() << " l ";

          if(y>0.1){
            fp << " cp 0 0 0 srgb s %ZVALUE: " << ((p3p+p2p)/2.0).get_z()+0.015 << "\n";
          }else{
            fp << " cp 1 1 1 srgb s %ZVALUE: " << ((p3p+p2p)/2.0).get_z()+0.015 << "\n";
	  }
        }

      }
    }
  }

  if(split_x||split_y) return;

  // Optional outline. We'll have to think of a way to trigger this selectively.
  fp << " 0.10 slw ";
  pfrac = Cartesian(p1.get_x()/xscale,p1.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << "n " << pps.get_x() << " " << pps.get_y() << " m ";

  pfrac = Cartesian(p2.get_x()/xscale,p2.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";

  pfrac = Cartesian(p3.get_x()/xscale,p3.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";

  pfrac = Cartesian(p4.get_x()/xscale,p4.get_y()/yscale,0);
  pps   = Cartesian(pfrac.get_x()*xscaleps+xoff,pfrac.get_y()*xscaleps+yoff,0);
  fp << pps.get_x() << " " << pps.get_y() << " l ";

  if(y>0.1){
    fp << " cp 0 0 0 srgb s %ZVALUE: " << ((p1+p2+p3+p4)/4.0).get_z()+0.015 << "\n";
  }else{
    fp << " cp 1 1 1 srgb s %ZVALUE: " << ((p1+p2+p3+p4)/4.0).get_z()+0.015 << "\n";
  }
}

void Primitive::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v)
{
   matrix unit_matrix(4,kdelta);
   Cartesian null_origin(0,0,0);
   Volume null_volume;
   std::vector<Primitive*> a = GetSimplePrimitives(null_volume,unit_matrix,null_origin);
   std::vector<Primitive*>::iterator a_iter = a.begin();
   while(a_iter!=a.end()){
     (*a_iter)->DrawPovray(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,v);
     a_iter++;
   }
}

void Primitive::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
   matrix unit_matrix(4,kdelta);
   Cartesian null_origin(0,0,0);
   Volume null_volume;
   std::vector<Primitive*> a = GetSimplePrimitives(null_volume,unit_matrix,null_origin);
   std::vector<Primitive*>::iterator a_iter = a.begin();
   while(a_iter!=a.end()){
     (*a_iter)->DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
     a_iter++;
   }
}


void ArrowHeadRibbon::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){
  Ribbon::DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
}

void BillBoard::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
   if(!image.get_width()||!image.get_height()||!image.get_pixels())
      return;
  fp << "gsave /rstr 8 string def /gstr 8 string def /bstr 8 string def %TEXT\n";
  fp << vertices[0].get_x()*xoff*2 << " " << vertices[0].get_y()*yoff*2 << " translate %TEXT\n";
  fp << ((image.get_width()+7)/8)*8*scale_w << " " << image.get_height()*scale_h << " scale %TEXT\n";
  // No need to apply flip as image is already flipped.
  //fp << ((image.get_width()+7)/8)*8 << " " << image.get_height() << " 8 [" << ((image.get_width()+7)/8)*8 << " 0 0 " << -image.get_height() << " 0 " << image.get_height() << "] %TEXT\n";
  fp << ((image.get_width()+7)/8)*8 << " " << image.get_height() << " 8 [" << ((image.get_width()+7)/8)*8 << " 0 0 " << image.get_height() << " 0 " << 0 << "] %TEXT\n";
  fp << "{ currentfile rstr readhexstring pop} bind %TEXT\n";
  fp << "{ currentfile gstr readhexstring pop} bind %TEXT\n";
  fp << "{ currentfile bstr readhexstring pop} bind %TEXT\n";
  fp << "true 3 colorimage ";
  image_info tmp_image = image;
  tmp_image.convert_rgb();
  unsigned char *pixels = tmp_image.get_pixels();
  fp << std::resetiosflags(std::ios_base::basefield)
     << std::setiosflags(std::ios_base::hex)
     << std::setiosflags(std::ios_base::uppercase)
     << std::setfill('0')
     << std::setw(2);
  for(int ii=0;ii<tmp_image.get_height();ii++){
    int nwrit = 0;
    int safe_width = ((tmp_image.get_width()*3)/24)*24;
    for(int kk=0;kk<safe_width;kk+=24,nwrit+=24){
      int idx = ii*(tmp_image.get_width()*3) + kk;
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+3]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+6]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+9]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+12]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+15]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+18]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+21]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+1]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+4]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+7]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+10]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+13]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+16]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+19]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+22]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+2]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+5]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+8]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+11]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+14]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+17]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+20]);
      fp << std::setfill('0') << std::setw(2) << int(pixels[idx+23]);
    }
    //std::cout << nwrit << " (" << (tmp_image.get_width()*3) << ")\n";
    int npix_left = (tmp_image.get_width()*3-nwrit)/3;
    //std::cout << npix_left << "\n";
    if(npix_left>0){
      for(int kk2=0;kk2<npix_left;kk2++,nwrit++){
        int idx = ii*(tmp_image.get_width()*3) + safe_width;
        fp << std::setfill('0') << std::setw(2) << int(pixels[idx+kk2*3]);
      }
      for(int kk2=0;kk2<8-npix_left;kk2++,nwrit++){
        fp << std::setfill('0') << std::setw(2) << 0;
      }
      for(int kk2=0;kk2<npix_left;kk2++,nwrit++){
        int idx = ii*(tmp_image.get_width()*3) + safe_width;
        fp << std::setfill('0') << std::setw(2) << int(pixels[idx+kk2*3+1]);
      }
      for(int kk2=0;kk2<8-npix_left;kk2++,nwrit++){
        fp << std::setfill('0') << std::setw(2) << 0;
      }
      for(int kk2=0;kk2<npix_left;kk2++,nwrit++){
        int idx = ii*(tmp_image.get_width()*3) + safe_width;
        fp << std::setfill('0') << std::setw(2) << int(pixels[idx+kk2*3+2]);
      }
      for(int kk2=0;kk2<8-npix_left;kk2++,nwrit++){
        fp << std::setfill('0') << std::setw(2) << 0;
      }
    }
    //std::cout << nwrit << " (" << (((tmp_image.get_width()+7)/8)*8*3) << ")\n";
  }
  fp << std::setfill(' ');
  fp << std::resetiosflags(std::ios_base::basefield);
  fp << " %TEXT\n";
  fp << "grestore %TEXT\n";
}

void BillBoardText::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)
{
  DrawPSmain(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v,IsBillBoard());
}

std::vector<Primitive*> Primitive::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "Primitive::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  /*
  // I would like to do this and have the runtime checker do the right thing. Not sure if it
  // is at all possible, since I can't allocate a primitive. Um. See Polygon Element !!!
  typeof(*this)* oself = dynamic_cast<typeof(this)>(new (typeof(*this))(vertices,colour,origin,alpha,textured));
  oself->SetNormals(normals);
  a.push_back(oself);
  */
  return a;
}

std::vector<Primitive*> Point::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  a.push_back(new Point(vertices[0],colour,vertices[0],size));
  return a;
}

std::vector<Primitive*> Line::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "Line::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  Line *line = new Line(vertices,colour,origin,size,alpha);
  a.push_back(line);
  return a;
}

std::vector<Primitive*> ArrowElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  return a;
}
std::vector<Primitive*> Arrow::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "Arrow::GetSimplePrimitives\n";
  return ArrowElement::GetSimplePrimitives(clip_vol,objrotmatrix,objorigin,start,end);
}

std::vector<Primitive*> DashArrowElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "DashArrowElement::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  return a;
}
std::vector<Primitive*> DashArrow::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "DashArrow::GetSimplePrimitives\n";
  return DashArrowElement::GetSimplePrimitives(clip_vol,objrotmatrix,objorigin,start,end);
}


std::vector<Primitive*> LineStrip::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "LineStrip::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  return a;
}

std::vector<Primitive*> DashLine::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "DashLine::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  return a;
}

std::vector<Primitive*> PartialCircleElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "Circle::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  return a;
}

std::vector<Primitive*> CircleElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "Circle::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  return a;
}

std::vector<Primitive*> TriangleElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  // TriangleElement?
  TriangleElement *tri = new Triangle(vertices,colour,origin,alpha,textured);
  tri->SetNormals(normals);
  a.push_back(tri);
  return a;
}

std::vector<Primitive*> PolygonElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  //typeof(*this)* tri = dynamic_cast<typeof(this)>(new (typeof(*this))(vertices,colour,origin,alpha,textured));
  PolygonElement *tri = new MGPolygon(vertices,colour,origin,alpha,textured);
  tri->SetNormals(normals);
  a.push_back(tri);
  return a;
}

std::vector<Primitive*> QuadElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  // QuadElement?
  QuadElement *tri = new Quad(vertices,colour,origin,alpha,textured);
  tri->SetNormals(normals);
  a.push_back(tri);
  return a;
}

std::vector<Primitive*> SimpleBillBoard::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  return a;
}

std::vector<Primitive*> SimpleText::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  return a;
}

std::vector<Primitive*> PointCollection::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  if(lines.size()<1) return a;
  std::vector<Primitive*>::const_iterator prim_iter = lines.begin();
  while(prim_iter!=lines.end()){
    //(*prim_iter)->set_transparent(transparent);
    //(*prim_iter)->SetAlpha(alpha);
    std::vector<Primitive*> b = (*prim_iter)->GetSimplePrimitives(clip_vol,objrotmatrix,objorigin,start,end);
    a.insert(a.end(),b.begin(),b.end());
    prim_iter++;
  }
  return a;
}
std::vector<Primitive*> LineCollection::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "LineCollection::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  if(lines.size()<1) return a;
  std::vector<Primitive*>::const_iterator prim_iter = lines.begin();
  while(prim_iter!=lines.end()){
    //(*prim_iter)->set_transparent(transparent);
    //(*prim_iter)->SetAlpha(alpha);
    std::vector<Primitive*> b = (*prim_iter)->GetSimplePrimitives(clip_vol,objrotmatrix,objorigin,start,end);
    a.insert(a.end(),b.begin(),b.end());
    prim_iter++;
  }
  return a;
}

std::vector<Primitive*> PolyCollection::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  if(prims.size()<1) return a;
  std::vector<Primitive*>::const_iterator prim_iter = prims.begin();
  while(prim_iter!=prims.end()){
    //(*prim_iter)->set_transparent(transparent);
    //(*prim_iter)->SetAlpha(alpha);
    std::vector<Primitive*> b = (*prim_iter)->GetSimplePrimitives(clip_vol,objrotmatrix,objorigin,start,end);
    a.insert(a.end(),b.begin(),b.end());
    prim_iter++;
  }
  return a;
}

std::vector<Primitive*> ImposterSphereCollection::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  if(prims.size()<1) return a;
  std::vector<Primitive*>::const_iterator prim_iter = prims.begin();
  while(prim_iter!=prims.end()){
    //(*prim_iter)->set_transparent(transparent);
    //(*prim_iter)->SetAlpha(alpha);
    std::vector<Primitive*> b = (*prim_iter)->GetSimplePrimitives(clip_vol,objrotmatrix,objorigin,start,end);
    a.insert(a.end(),b.begin(),b.end());
    prim_iter++;
  }
  return a;
}

std::vector<Primitive*> QuadStripElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  unsigned int i = 0;
  while(k!=vertices.end()){
    if(i%2==0&&i<vertices.size()+3){
      std::vector<Cartesian> v(4);
      v[0] = vertices[i];
      v[1] = vertices[i+1];
      v[2] = vertices[i+2];
      v[3] = vertices[i+3];
      QuadElement *q = new QuadElement(v,colour,v[0],alpha,textured);
      a.push_back(q);
    }
    k++; i++;
  }
  return a;
}

std::vector<Primitive*> TriangleStripElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  unsigned int i = 0;
  while(k!=vertices.end()){
    if(i<vertices.size()+3){
      std::vector<Cartesian> v(3);
      v[0] = vertices[i];
      v[1] = vertices[i+1];
      v[2] = vertices[i+2];
      TriangleElement *q = new TriangleElement(v,colour,v[0],alpha,textured);
      a.push_back(q);
    }
    k++; i++;
  }
  return a;
}

std::vector<Primitive*> TriangleFanElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  std::vector<Cartesian>::const_iterator k = vertices.begin();
  unsigned int i = 0;
  Plane plane(vertices[0],vertices[1],vertices[2]);
  plane.Normalize();
  Cartesian _normal = plane.get_normal();
  std::vector<Cartesian> _normals;
  _normals.push_back(_normal);
  _normals.push_back(_normal);
  _normals.push_back(_normal);
  while(k!=vertices.end()){
    if(i>=2){
      std::vector<Cartesian> v(3);
      v[0] = vertices[0];
      v[1] = vertices[i];
      v[2] = vertices[i-1];
      TriangleElement *q = new TriangleElement(v,colour,v[0],alpha,textured);
      q->SetNormals(_normals);
      a.push_back(q);
    }
    k++; i++;
  }
  return a;
}

std::vector<Primitive*> PolyCylinder::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;

  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;
  const double *colp = GetColour();
  std::vector<std::vector<Cartesian> > AB = GetPolyCylinderAB(vertices,size);

  for(int i=0;i<360;i=i+360/accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
      double col[4];
      if(multicolour){
        col[0] = colour_vector[0].get_x();
        col[1] = colour_vector[0].get_y();
        col[2] = colour_vector[0].get_z();
        col[3] = colour_vector[0].get_a();
        colp = col;
      }
      Cartesian p0 = vertices[0];
      //Cartesian n = vertices[1]-vertices[0];
      //n.normalize();
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p1 = x*AB[0][0] + y*AB[1][0] + vertices[0];
      x = cos(theta2);
      y = sin(theta2);
      Cartesian p2 = x*AB[0][0] + y*AB[1][0] + vertices[0];
      v.push_back(p0);
      v.push_back(p1);
      v.push_back(p2);
      TriangleElement *q = new TriangleElement(v,colp,(v[0]+v[1]+v[2])/3.,alpha,textured);
      a.push_back(q);
  }

  for(int i=0;i<360;i=i+360/accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
      int j = vertices.size()-1;
      double col[4];
      if(multicolour){
        col[0] = colour_vector[j].get_x();
        col[1] = colour_vector[j].get_y();
        col[2] = colour_vector[j].get_z();
        col[3] = colour_vector[j].get_a();
        colp = col;
      }
      Cartesian p0 = vertices[j];
      //Cartesian n = vertices[1]-vertices[0];
      //n.normalize();
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p1 = x*AB[0][j] + y*AB[1][j] + vertices[j];
      x = cos(theta2);
      y = sin(theta2);
      Cartesian p2 = x*AB[0][j] + y*AB[1][j] + vertices[j];
      v.push_back(p0);
      v.push_back(p2);
      v.push_back(p1);
      TriangleElement *q = new TriangleElement(v,colp,(v[0]+v[1]+v[2])/3.,alpha,textured);
      a.push_back(q);
  }

  for(int i=0;i<360;i=i+360/accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
    std::vector<Cartesian> ns;
    for(unsigned int j=0;j<vertices.size();j++){
      double col[4];
      if(multicolour){
        col[0] = colour_vector[j].get_x();
        col[1] = colour_vector[j].get_y();
        col[2] = colour_vector[j].get_z();
        col[3] = colour_vector[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p1 = x*AB[0][j] + y*AB[1][j] + vertices[j];
      Cartesian n1 = AB[0][j]*(x/(AB[0][j].length()*AB[0][j].length())) + AB[1][j]*(y/(AB[1][j].length()*AB[1][j].length()));
      x = cos(theta2);
      y = sin(theta2);
      Cartesian p2 = x*AB[0][j] + y*AB[1][j] + vertices[j];
      Cartesian n2 = AB[0][j]*(x/(AB[0][j].length()*AB[0][j].length())) + AB[1][j]*(y/(AB[1][j].length()*AB[1][j].length()));
      if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n2);
        v.push_back(p1);
        ns.push_back(n1);
      } else {
        v.push_back(p1);
        ns.push_back(n1);
        v.push_back(p2);
        ns.push_back(n2);
      }
      if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
      }
    }
  }

  return a;
}

void get_flat_ribbon(std::vector<Primitive*> &a,const std::vector<Cartesian> &vertices_tmp,const std::vector<std::vector<Cartesian> > &AB,int my_accu,int textured,int multicolour,const std::vector<Cartesian> &colour_vector_tmp, const double *colp, double alpha, const bool two_colour, const bool grey_ribbon_edge){
  Cartesian n;
  Cartesian p1,p2;
  unsigned int j;

  double col[4];
  int idx = vertices_tmp.size()/2;
  Cartesian p_t = vertices_tmp[idx+1] + vertices_tmp[idx-1] - 2*vertices_tmp[idx];
  p_t.normalize();
  Cartesian pvpr_t = AB[1][idx];
  pvpr_t.normalize();

  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;
  std::vector<Cartesian> v;
  std::vector<Cartesian> ns;
  double colour[4] = {0.6,0.6,0.6,1.0};
  if(face_one&&two_colour){
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    if((multicolour&&!face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    p2 = vertices_tmp[j] + AB[0][j] + AB[1][j];
    p1 = vertices_tmp[j] - AB[0][j] + AB[1][j];

    n= AB[1][j];
    n.normalize();
    if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n);
        v.push_back(p1);
        ns.push_back(n);
     } else {
        v.push_back(p1);
        ns.push_back(n);
        v.push_back(p2);
        ns.push_back(n);
     }
     if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
     }
  }

  v.clear();
  ns.clear();
  if(!face_one&&two_colour){
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    if((multicolour&&face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    p2 = vertices_tmp[j] - AB[0][j] - AB[1][j];
    p1 = vertices_tmp[j] + AB[0][j] - AB[1][j];

    n= -AB[1][j];
    n.normalize();
    if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n);
        v.push_back(p1);
        ns.push_back(n);
     } else {
        v.push_back(p1);
        ns.push_back(n);
        v.push_back(p2);
        ns.push_back(n);
     }
     if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
     }
  }

  v.clear();
  ns.clear();
  for(j=0;j<vertices_tmp.size();j++){
    if(multicolour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }

    if(grey_ribbon_edge)
      colp=colour;

    p2 = vertices_tmp[j] + AB[0][j] - AB[1][j];
    p1 = vertices_tmp[j] + AB[0][j] + AB[1][j];

    n= AB[0][j];
    n.normalize();
    if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n);
        v.push_back(p1);
        ns.push_back(n);
     } else {
        v.push_back(p1);
        ns.push_back(n);
        v.push_back(p2);
        ns.push_back(n);
     }
     if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
     }
  }

  v.clear();
  ns.clear();
  for(j=0;j<vertices_tmp.size();j++){
    if(multicolour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }

    if(grey_ribbon_edge)
      colp=colour;

    p2 = vertices_tmp[j] - AB[0][j] + AB[1][j];
    p1 = vertices_tmp[j] - AB[0][j] - AB[1][j];

    n= -AB[0][j];
    n.normalize();
    if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n);
        v.push_back(p1);
        ns.push_back(n);
     } else {
        v.push_back(p1);
        ns.push_back(n);
        v.push_back(p2);
        ns.push_back(n);
     }
     if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
     }
  }
}

void get_flat_rounded_ribbon(std::vector<Primitive*> &a,const std::vector<Cartesian> &vertices_tmp,const std::vector<std::vector<Cartesian> > &AB,int my_accu,int textured,int multicolour,const std::vector<Cartesian> &colour_vector_tmp, const double *colp, double alpha, const bool two_colour){
  Cartesian n;
  Cartesian p1,p2;
  unsigned int j;

  double col[4];
  int idx = vertices_tmp.size()/2;
  Cartesian p_t = vertices_tmp[idx+1] + vertices_tmp[idx-1] - 2*vertices_tmp[idx];
  p_t.normalize();
  Cartesian pvpr_t = AB[1][idx];
  pvpr_t.normalize();

  double thickness = AB[1][0].length(); // Assume constant thickness.
  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;
  std::vector<Cartesian> v;
  std::vector<Cartesian> ns;
  double colour[4] = {0.6,0.6,0.6,1.0};
  if(face_one&&two_colour){
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    Cartesian p = AB[0][j];
    p.normalize(p.length()-thickness);
    if((multicolour&&!face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    p2 = vertices_tmp[j] + p + AB[1][j];
    p1 = vertices_tmp[j] - p + AB[1][j];

    n= AB[1][j];
    n.normalize();
    if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n);
        v.push_back(p1);
        ns.push_back(n);
     } else {
        v.push_back(p1);
        ns.push_back(n);
        v.push_back(p2);
        ns.push_back(n);
     }
     if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
     }
  }

  v.clear();
  ns.clear();
  if(!face_one&&two_colour){
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    Cartesian p = AB[0][j];
    p.normalize(p.length()-thickness);
    if((multicolour&&face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    p2 = vertices_tmp[j] - p - AB[1][j];
    p1 = vertices_tmp[j] + p - AB[1][j];

    n= -AB[1][j];
    n.normalize();
    if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n);
        v.push_back(p1);
        ns.push_back(n);
     } else {
        v.push_back(p1);
        ns.push_back(n);
        v.push_back(p2);
        ns.push_back(n);
     }
     if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
     }
  }

  v.clear();
  ns.clear();
  for(int i=-90+360/my_accu;i<=90;i=i+360/my_accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/my_accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
    std::vector<Cartesian> ns;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p = AB[0][j];
      p.normalize(p.length()-thickness);
      Cartesian pr = AB[0][j];
      pr.normalize(thickness);
      Cartesian ppr = AB[1][j];
      ppr.normalize(thickness);
      Cartesian p1 = x*pr + y*ppr + vertices_tmp[j]+p;
      Cartesian n1 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      x = cos(theta2);
      y = sin(theta2);
      Cartesian p2 = x*pr + y*ppr + vertices_tmp[j]+p;
      Cartesian n2 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n2);
        v.push_back(p1);
        ns.push_back(n1);
      } else {
        v.push_back(p1);
        ns.push_back(n1);
        v.push_back(p2);
        ns.push_back(n2);
      }
      if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
      }
    }
  }
  v.clear();
  ns.clear();
  for(int i=90+360/my_accu;i<=270;i=i+360/my_accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/my_accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
    std::vector<Cartesian> ns;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p = AB[0][j];
      p.normalize(p.length()-thickness);
      Cartesian pr = AB[0][j];
      pr.normalize(thickness);
      Cartesian ppr = AB[1][j];
      ppr.normalize(thickness);
      Cartesian p1 = x*pr + y*ppr + vertices_tmp[j]-p;
      Cartesian n1 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      x = cos(theta2);
      y = sin(theta2);
      Cartesian p2 = x*pr + y*ppr + vertices_tmp[j]-p;
      Cartesian n2 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n2);
        v.push_back(p1);
        ns.push_back(n1);
      } else {
        v.push_back(p1);
        ns.push_back(n1);
        v.push_back(p2);
        ns.push_back(n2);
      }
      if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
      }
    }
  }
}

void get_fancy_ribbon(std::vector<Primitive*> &a,const std::vector<Cartesian> &vertices_tmp,const std::vector<std::vector<Cartesian> > &AB,int my_accu,int textured,int multicolour,const std::vector<Cartesian> &colour_vector_tmp, const double *colp, double alpha, const bool two_colour){
  Cartesian n;
  Cartesian p1,p2;
  unsigned int j;

  double col[4];
  int idx = vertices_tmp.size()/2;
  Cartesian p_t = vertices_tmp[idx+1] + vertices_tmp[idx-1] - 2*vertices_tmp[idx];
  p_t.normalize();
  Cartesian pvpr_t = AB[1][idx];
  pvpr_t.normalize();

  double thickness = AB[1][0].length(); // Assume constant thickness.
  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;
  std::vector<Cartesian> v;
  std::vector<Cartesian> ns;
  double colour[4] = {0.6,0.6,0.6,1.0};
  if(face_one&&two_colour){
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    Cartesian p = AB[0][j];
    p.normalize(p.length()-thickness*(1+sin(atan(0.7))));
    if((multicolour&&!face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    p2 = vertices_tmp[j] + p + 0.7*AB[1][j];
    p1 = vertices_tmp[j] - p + 0.7*AB[1][j];

    n= AB[1][j];
    n.normalize();
    if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n);
        v.push_back(p1);
        ns.push_back(n);
     } else {
        v.push_back(p1);
        ns.push_back(n);
        v.push_back(p2);
        ns.push_back(n);
     }
     if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
     }
  }

  v.clear();
  ns.clear();
  if(!face_one&&two_colour){
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    Cartesian p = AB[0][j];
    p.normalize(p.length()-thickness*(1+sin(atan(0.7))));
    if((multicolour&&face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    p2 = vertices_tmp[j] - p - 0.7*AB[1][j];
    p1 = vertices_tmp[j] + p - 0.7*AB[1][j];

    n= -AB[1][j];
    n.normalize();
    if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n);
        v.push_back(p1);
        ns.push_back(n);
     } else {
        v.push_back(p1);
        ns.push_back(n);
        v.push_back(p2);
        ns.push_back(n);
     }
     if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
     }
  }

  v.clear();
  ns.clear();
  for(int i=0;i<=360;i=i+360/my_accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/my_accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
    std::vector<Cartesian> ns;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p = AB[0][j];
      p.normalize(p.length()-thickness);
      Cartesian pr = AB[0][j];
      pr.normalize(thickness);
      Cartesian ppr = AB[1][j];
      ppr.normalize(thickness);
      Cartesian p1 = x*pr + y*ppr + vertices_tmp[j]+p;
      Cartesian n1 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      x = cos(theta2);
      y = sin(theta2);
      Cartesian p2 = x*pr + y*ppr + vertices_tmp[j]+p;
      Cartesian n2 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n2);
        v.push_back(p1);
        ns.push_back(n1);
      } else {
        v.push_back(p1);
        ns.push_back(n1);
        v.push_back(p2);
        ns.push_back(n2);
      }
      if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
      }
    }
  }
  v.clear();
  ns.clear();
  for(int i=0;i<=360;i=i+360/my_accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/my_accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
    std::vector<Cartesian> ns;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p = AB[0][j];
      p.normalize(p.length()-thickness);
      Cartesian pr = AB[0][j];
      pr.normalize(thickness);
      Cartesian ppr = AB[1][j];
      ppr.normalize(thickness);
      Cartesian p1 = x*pr + y*ppr + vertices_tmp[j]-p;
      Cartesian n1 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      x = cos(theta2);
      y = sin(theta2);
      Cartesian p2 = x*pr + y*ppr + vertices_tmp[j]-p;
      Cartesian n2 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n2);
        v.push_back(p1);
        ns.push_back(n1);
      } else {
        v.push_back(p1);
        ns.push_back(n1);
        v.push_back(p2);
        ns.push_back(n2);
      }
      if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
      }
    }
  }


}

void Ribbon::get_fancy_ribbon_array(std::vector<Primitive*> &a,const std::vector<Cartesian> &vertices_tmp,const std::vector<std::vector<Cartesian> > &AB,int my_accu,int textured,int multicolour,const std::vector<Cartesian> &colour_vector_tmp, const double *colp, double alpha, const bool two_colour){

  Cartesian n;
  Cartesian p1,p2;
  unsigned int j;

  nVerts = vertices_tmp.size() * 8;
  //std::cout << "Ribbon::get_fancy_ribbon_array " << vertices_tmp.size() << " " << my_accu << "\n";
  bufferVertices = new float[3*(nVerts+2*vertices_tmp.size() * (my_accu))+6 + 6*(my_accu+1)]; // Need 2* extra in fact!!
  bufferNormals = new float[3*(nVerts+2*vertices_tmp.size() * (my_accu))+6 + 6*(my_accu+1)];
  bufferColors = new float[4*(nVerts+2*vertices_tmp.size() * (my_accu))+8 + 8*(my_accu+1)];
  nRealVerts = 3*(nVerts+2*vertices_tmp.size() * (my_accu))+6 + 6*(my_accu+1);
  // Tri-strip alloc
  bufferIndices = new GLuint[4*(vertices_tmp.size())+2 + 1 + 2*((2*vertices_tmp.size()+2)*(my_accu) -2 + my_accu*4+8)];
  // Quad alloc
  //bufferIndices = new GLuint[16*(vertices_tmp.size()-1)+2*4*(vertices_tmp.size()-1)*(my_accu)];

  double col[4];
  int idx = vertices_tmp.size()/2;
  Cartesian p_t = vertices_tmp[idx+1] + vertices_tmp[idx-1] - 2*vertices_tmp[idx];
  p_t.normalize();
  Cartesian pvpr_t = AB[1][idx];
  pvpr_t.normalize();

  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;
  std::vector<Cartesian> v;
  std::vector<Cartesian> ns;
  double colour[4] = {0.6,0.6,0.6,1.0};
  if(face_one&&two_colour){
    colp=colour;
  }
  
    for(j=0;j<vertices_tmp.size();j++){
    if((multicolour&&!face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    Cartesian p = AB[0][j];
    //p.normalize(p.length()-thickness*(1+sin(atan(0.7))));
    double l =p.length()-thickness*(1+sin(atan(0.7))); 
    if(l<0) l=0;
    p.normalize(l);
    p1 = vertices_tmp[j] + p + 0.7*AB[1][j];
    p2 = vertices_tmp[j] - p + 0.7*AB[1][j];

    n =  AB[1][j];
    n.normalize();

    bufferVertices[6*j]   = p1.get_x();
    bufferVertices[6*j+1] = p1.get_y();
    bufferVertices[6*j+2] = p1.get_z();
    bufferVertices[6*j+3] = p2.get_x();
    bufferVertices[6*j+4] = p2.get_y();
    bufferVertices[6*j+5] = p2.get_z();
    bufferNormals[6*j]   = n.get_x();
    bufferNormals[6*j+1] = n.get_y();
    bufferNormals[6*j+2] = n.get_z();
    bufferNormals[6*j+3] = n.get_x();
    bufferNormals[6*j+4] = n.get_y();
    bufferNormals[6*j+5] = n.get_z();
    bufferColors[8*j]     = colp[0];
    bufferColors[8*j+1]   = colp[1];
    bufferColors[8*j+2]   = colp[2];
    bufferColors[8*j+3]   = colp[3];
    bufferColors[8*j+4]   = colp[0];
    bufferColors[8*j+5]   = colp[1];
    bufferColors[8*j+6]   = colp[2];
    bufferColors[8*j+7]   = colp[3];

    bufferIndices[2*j]   = 2*j;
    bufferIndices[2*j+1] = 2*j+1;


  }

  bufferIndices[2*(vertices_tmp.size())]   = bufferIndices[2*(vertices_tmp.size())-1];
  bufferIndices[2*(vertices_tmp.size())+1] = 2*vertices_tmp.size()+1;

  if(!face_one&&two_colour){
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    if((multicolour&&face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    Cartesian p = AB[0][j];
    //p.normalize(p.length()-thickness*(1+sin(atan(0.7))));
    double l =p.length()-thickness*(1+sin(atan(0.7))); 
    if(l<0) l=0;
    p.normalize(l);
    p1 = vertices_tmp[j] + p - 0.7*AB[1][j];
    p2 = vertices_tmp[j] - p - 0.7*AB[1][j];

    n =  -AB[1][j];
    n.normalize();

    bufferVertices[6*(vertices_tmp.size())+6*j]   = p1.get_x();
    bufferVertices[6*(vertices_tmp.size())+6*j+1] = p1.get_y();
    bufferVertices[6*(vertices_tmp.size())+6*j+2] = p1.get_z();
    bufferVertices[6*(vertices_tmp.size())+6*j+3] = p2.get_x();
    bufferVertices[6*(vertices_tmp.size())+6*j+4] = p2.get_y();
    bufferVertices[6*(vertices_tmp.size())+6*j+5] = p2.get_z();
    bufferNormals[6*(vertices_tmp.size())+6*j]   = n.get_x();
    bufferNormals[6*(vertices_tmp.size())+6*j+1] = n.get_y();
    bufferNormals[6*(vertices_tmp.size())+6*j+2] = n.get_z();
    bufferNormals[6*(vertices_tmp.size())+6*j+3] = n.get_x();
    bufferNormals[6*(vertices_tmp.size())+6*j+4] = n.get_y();
    bufferNormals[6*(vertices_tmp.size())+6*j+5] = n.get_z();
    bufferColors[8*(vertices_tmp.size())+8*j]     = colp[0];
    bufferColors[8*(vertices_tmp.size())+8*j+1]   = colp[1];
    bufferColors[8*(vertices_tmp.size())+8*j+2]   = colp[2];
    bufferColors[8*(vertices_tmp.size())+8*j+3]   = colp[3];
    bufferColors[8*(vertices_tmp.size())+8*j+4]   = colp[0];
    bufferColors[8*(vertices_tmp.size())+8*j+5]   = colp[1];
    bufferColors[8*(vertices_tmp.size())+8*j+6]   = colp[2];
    bufferColors[8*(vertices_tmp.size())+8*j+7]   = colp[3];

    bufferIndices[2*(vertices_tmp.size())+2+2*j+1]= 2*vertices_tmp.size()+2*j;
    bufferIndices[2*(vertices_tmp.size())+2+2*j]  = 2*vertices_tmp.size()+2*j+1;

  }

  bufferIndices[4*(vertices_tmp.size())+2] = bufferIndices[4*(vertices_tmp.size())+1];

  int vertexOffset = 12*(vertices_tmp.size());
  int colourOffset = 16*(vertices_tmp.size());
  int indexOffset = 4*(vertices_tmp.size());
  int indexOffsetInternal = 4*(vertices_tmp.size())+2+1;
  int ii=0;
  int maxInd = 0;

  for(int i=0;i<360;i=i+360/my_accu,ii++){
    double theta = (double)i/360.0 * PIBY2;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      if(ii<my_accu-1){
        if(j==0){
          // Prepend a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii)*vertices_tmp.size();
          maxInd++;
        }
        bufferIndices[indexOffsetInternal++] = indexOffset + j + ii*vertices_tmp.size();
        maxInd++;
        bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii+1)*vertices_tmp.size();
        maxInd++;
        if(j==vertices_tmp.size()-1){
          // Append a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii+1)*vertices_tmp.size();
          maxInd++;
        }
      } else {
        if(j==0){
          // Prepend a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii)*vertices_tmp.size();
          maxInd++;
        }
        bufferIndices[indexOffsetInternal++] = indexOffset + j + ii*vertices_tmp.size();
        maxInd++;
        bufferIndices[indexOffsetInternal++] = indexOffset + j;
        maxInd++;
        if(j==vertices_tmp.size()-1){
          // Append a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j;
          maxInd++;
        }
      }
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p = AB[0][j];
      p.normalize(p.length()-thickness);
      Cartesian pr = AB[0][j];
      pr.normalize(thickness);
      Cartesian ppr = AB[1][j];
      ppr.normalize(thickness);
      Cartesian p1 = x*pr + y*ppr + vertices_tmp[j]+p;
      Cartesian n1 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j]   = p1.get_x();
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+1] = p1.get_y();
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+2] = p1.get_z();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j]    = n1.get_x();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+1]  = n1.get_y();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+2]  = n1.get_z();
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j]     = colp[0];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+1]   = colp[1];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+2]   = colp[2];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+3]   = colp[3];
    }
  }

  //std::cout << "vertices.size(): " << vertices_tmp.size() << "\n";
  //std::cout << "ii: " << ii << "\n";
  //std::cout << "my_accu: " << my_accu << "\n";
  //std::cout << "maxInd: " << maxInd << "\n";
  //std::cout << "indexOffsetInternal: " << indexOffsetInternal << "\n";
  //std::cout << "Allocated: " << 4*(vertices_tmp.size())+2 + 1 + 2*((2*vertices_tmp.size()+2)*(my_accu) -2 + my_accu*4+8) << "\n";

  //ii--;
  vertexOffset += 3*(ii*vertices_tmp.size());
  colourOffset += 4*(ii*vertices_tmp.size());
  indexOffset += ii*(vertices_tmp.size());

  //indexOffsetInternal += 4*(vertices_tmp.size()-1)*ii;

  ii=0;
  for(int i=0;i<360;i=i+360/my_accu,ii++){
    double theta = (double)i/360.0 * PIBY2;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      if(ii<my_accu-1){
        if(j==0){
          // Prepend a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii)*vertices_tmp.size();
          maxInd++;
        }
        bufferIndices[indexOffsetInternal++] = indexOffset + j + ii*vertices_tmp.size();
        maxInd++;
        bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii+1)*vertices_tmp.size();
        maxInd++;
        if(j==vertices_tmp.size()-1){
          // Append a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii+1)*vertices_tmp.size();
          maxInd++;
        }
      } else {
        if(j==0){
          // Prepend a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii)*vertices_tmp.size();
          maxInd++;
        }
        bufferIndices[indexOffsetInternal++] = indexOffset + j + ii*vertices_tmp.size();
        maxInd++;
        bufferIndices[indexOffsetInternal++] = indexOffset + j;
        maxInd++;
        if(j==vertices_tmp.size()-1){
          // Append a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j;
          maxInd++;
        }
      }
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p = AB[0][j];
      p.normalize(p.length()-thickness);
      Cartesian pr = AB[0][j];
      pr.normalize(thickness);
      Cartesian ppr = AB[1][j];
      ppr.normalize(thickness);
      Cartesian p1 = x*pr + y*ppr + vertices_tmp[j]-p;
      Cartesian n1 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j]   = p1.get_x();
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+1] = p1.get_y();
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+2] = p1.get_z();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j]    = n1.get_x();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+1]  = n1.get_y();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+2]  = n1.get_z();
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j]     = colp[0];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+1]   = colp[1];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+2]   = colp[2];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+3]   = colp[3];
    }
  }

  //std::cout << "ii(2): " << ii << "\n";
  //std::cout << "maxInd(2): " << maxInd << "\n";

  //nVerts = 8*(vertices_tmp.size()-1)+8*(vertices_tmp.size()-1)*my_accu;
  //polygonType = GL_QUADS;


  nVerts = 4*(vertices_tmp.size())+2 +maxInd;

  //ii--;
  vertexOffset += 3*(ii*vertices_tmp.size());
  colourOffset += 4*(ii*vertices_tmp.size());

  // Start cap

  ii=0;
  Cartesian n1 = vertices_tmp[0] - vertices_tmp[1];
  for(int i=-90;i<=270;i=i+360/my_accu,ii++){
    double theta = (double)i/360.0 * PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    Cartesian p1 = x*AB[0][0] + y*AB[1][0] + vertices_tmp[0];
    bufferVertices[vertexOffset+3*ii]   = p1.get_x();
    bufferVertices[vertexOffset+3*ii+1] = p1.get_y();
    bufferVertices[vertexOffset+3*ii+2] = p1.get_z();
    bufferNormals[vertexOffset+3*ii]   = n1.get_x();
    bufferNormals[vertexOffset+3*ii+1] = n1.get_y();
    bufferNormals[vertexOffset+3*ii+2] = n1.get_z();
    bufferColors[colourOffset + 4*ii]     = colour_vector_tmp[0].get_x();
    bufferColors[colourOffset + 4*ii+1]   = colour_vector_tmp[0].get_y();
    bufferColors[colourOffset + 4*ii+2]   = colour_vector_tmp[0].get_z();
    bufferColors[colourOffset + 4*ii+3]   = colour_vector_tmp[0].get_a();
    if(ii<my_accu){
      bufferIndices[indexOffsetInternal++] = vertexOffset/3;
      bufferIndices[indexOffsetInternal++] = vertexOffset/3 + ii;
    } else {
      bufferIndices[indexOffsetInternal++] = vertexOffset/3 + (my_accu-ii);
      bufferIndices[indexOffsetInternal++] = vertexOffset/3;
    }
  }

  nVerts += 2*my_accu+2;
  // End cap

  vertexOffset += 3*(ii);
  colourOffset += 4*(ii);

  ii=0;
  n1 = vertices_tmp[vertices_tmp.size()-1] - vertices_tmp[vertices_tmp.size()-2];
  for(int i=-90;i<=270;i=i+360/my_accu,ii++){
    double theta = (double)i/360.0 * PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    Cartesian p1 = x*AB[0][vertices_tmp.size()-1] + y*AB[1][vertices_tmp.size()-1] + vertices_tmp[vertices_tmp.size()-1];
    bufferVertices[vertexOffset+3*ii]   = p1.get_x();
    bufferVertices[vertexOffset+3*ii+1] = p1.get_y();
    bufferVertices[vertexOffset+3*ii+2] = p1.get_z();
    bufferNormals[vertexOffset+3*ii]   = n1.get_x();
    bufferNormals[vertexOffset+3*ii+1] = n1.get_y();
    bufferNormals[vertexOffset+3*ii+2] = n1.get_z();
    bufferColors[colourOffset + 4*ii]     = colour_vector_tmp[0].get_x();
    bufferColors[colourOffset + 4*ii+1]   = colour_vector_tmp[0].get_y();
    bufferColors[colourOffset + 4*ii+2]   = colour_vector_tmp[0].get_z();
    bufferColors[colourOffset + 4*ii+3]   = colour_vector_tmp[0].get_a();
    if(ii<my_accu){
      bufferIndices[indexOffsetInternal++] = vertexOffset/3 + ii;
      bufferIndices[indexOffsetInternal++] = vertexOffset/3;
    } else {
      bufferIndices[indexOffsetInternal++] = vertexOffset/3;
      bufferIndices[indexOffsetInternal++] = vertexOffset/3 + (my_accu-ii);
    }
  }

  nVerts += 2*my_accu+2;

  polygonType = GL_TRIANGLE_STRIP;

}

void Ribbon::get_flat_rounded_ribbon_array(std::vector<Primitive*> &a,const std::vector<Cartesian> &vertices_tmp,const std::vector<std::vector<Cartesian> > &AB,int my_accu,int textured,int multicolour,const std::vector<Cartesian> &colour_vector_tmp, const double *colp, double alpha, const bool two_colour){

  Cartesian n;
  Cartesian p1,p2;
  unsigned int j;

  nVerts = vertices_tmp.size() * 8;
  //std::cout << "Ribbon::get_flat_rounded_ribbon_array " << vertices_tmp.size() << " " << my_accu << "\n";
  // +6 for 2 extra bits for centres of start and end caps
  bufferVertices = new float[3*(nVerts+vertices_tmp.size() * (my_accu+2))+6 + 6*(my_accu+1)];
  bufferNormals = new float[3*(nVerts+vertices_tmp.size() * (my_accu+2))+6 + 6*(my_accu+1)];
  bufferColors = new float[4*(nVerts+vertices_tmp.size() * (my_accu+2))+8 + 8*(my_accu+1)];
  nRealVerts = 3*(nVerts+vertices_tmp.size() * (my_accu+2))+6 + 6*(my_accu+1);
  //bufferIndices = new GLuint[16*(vertices_tmp.size()-1)+4*(vertices_tmp.size()-1)*(my_accu+2)];
  bufferIndices = new GLuint[4*(vertices_tmp.size())+2 + 2*((vertices_tmp.size()+1)*my_accu)+my_accu*4+8];

  // centres of start and end caps
  bufferVertices[3*(nVerts+vertices_tmp.size() * (my_accu+2))]   = vertices_tmp[0].get_x();
  bufferVertices[3*(nVerts+vertices_tmp.size() * (my_accu+2))+1] = vertices_tmp[0].get_y();
  bufferVertices[3*(nVerts+vertices_tmp.size() * (my_accu+2))+2] = vertices_tmp[0].get_z();
  bufferVertices[3*(nVerts+vertices_tmp.size() * (my_accu+2))+3] = vertices_tmp[vertices_tmp.size()-1].get_x();
  bufferVertices[3*(nVerts+vertices_tmp.size() * (my_accu+2))+4] = vertices_tmp[vertices_tmp.size()-1].get_y();
  bufferVertices[3*(nVerts+vertices_tmp.size() * (my_accu+2))+5] = vertices_tmp[vertices_tmp.size()-1].get_z();
  n = vertices_tmp[0]-vertices_tmp[1];
  bufferNormals[3*(nVerts+vertices_tmp.size() * (my_accu+2))]    = n.get_x();
  bufferNormals[3*(nVerts+vertices_tmp.size() * (my_accu+2))+1]  = n.get_y();
  bufferNormals[3*(nVerts+vertices_tmp.size() * (my_accu+2))+2]  = n.get_z();
  n = vertices_tmp[vertices_tmp.size()-1]-vertices_tmp[vertices_tmp.size()-2];
  bufferNormals[3*(nVerts+vertices_tmp.size() * (my_accu+2))+3]  = n.get_x();
  bufferNormals[3*(nVerts+vertices_tmp.size() * (my_accu+2))+4]  = n.get_y();
  bufferNormals[3*(nVerts+vertices_tmp.size() * (my_accu+2))+5]  = n.get_z();
  bufferColors[4*(nVerts+vertices_tmp.size() * (my_accu+2))]     = colour_vector_tmp[0].get_x();
  bufferColors[4*(nVerts+vertices_tmp.size() * (my_accu+2))+1]   = colour_vector_tmp[0].get_y();
  bufferColors[4*(nVerts+vertices_tmp.size() * (my_accu+2))+2]   = colour_vector_tmp[0].get_z();
  bufferColors[4*(nVerts+vertices_tmp.size() * (my_accu+2))+3]   = colour_vector_tmp[0].get_a();
  bufferColors[4*(nVerts+vertices_tmp.size() * (my_accu+2))+4]   = colour_vector_tmp[colour_vector_tmp.size()-1].get_x();
  bufferColors[4*(nVerts+vertices_tmp.size() * (my_accu+2))+5]   = colour_vector_tmp[colour_vector_tmp.size()-1].get_y();
  bufferColors[4*(nVerts+vertices_tmp.size() * (my_accu+2))+6]   = colour_vector_tmp[colour_vector_tmp.size()-1].get_z();
  bufferColors[4*(nVerts+vertices_tmp.size() * (my_accu+2))+7]   = colour_vector_tmp[colour_vector_tmp.size()-1].get_a();

  double col[4];
  int idx = vertices_tmp.size()/2;
  Cartesian p_t = vertices_tmp[idx+1] + vertices_tmp[idx-1] - 2*vertices_tmp[idx];
  p_t.normalize();
  Cartesian pvpr_t = AB[1][idx];
  pvpr_t.normalize();

  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;
  std::vector<Cartesian> v;
  std::vector<Cartesian> ns;
  double colour[4] = {0.6,0.6,0.6,1.0};
  if(face_one&&two_colour){
    colp=colour;
  }
  
  for(j=0;j<vertices_tmp.size();j++){
    if((multicolour&&!face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    Cartesian p = AB[0][j];
    p.normalize(p.length()-thickness);
    p1 = vertices_tmp[j] + p + AB[1][j];
    p2 = vertices_tmp[j] - p + AB[1][j];

    n =  AB[1][j];
    n.normalize();

    bufferVertices[6*j]   = p1.get_x();
    bufferVertices[6*j+1] = p1.get_y();
    bufferVertices[6*j+2] = p1.get_z();
    bufferVertices[6*j+3] = p2.get_x();
    bufferVertices[6*j+4] = p2.get_y();
    bufferVertices[6*j+5] = p2.get_z();
    bufferNormals[6*j]   = n.get_x();
    bufferNormals[6*j+1] = n.get_y();
    bufferNormals[6*j+2] = n.get_z();
    bufferNormals[6*j+3] = n.get_x();
    bufferNormals[6*j+4] = n.get_y();
    bufferNormals[6*j+5] = n.get_z();
    bufferColors[8*j]     = colp[0];
    bufferColors[8*j+1]   = colp[1];
    bufferColors[8*j+2]   = colp[2];
    bufferColors[8*j+3]   = colp[3];
    bufferColors[8*j+4]   = colp[0];
    bufferColors[8*j+5]   = colp[1];
    bufferColors[8*j+6]   = colp[2];
    bufferColors[8*j+7]   = colp[3];

    bufferIndices[2*j]   = 2*j;
    bufferIndices[2*j+1] = 2*j+1;

  }

  bufferIndices[2*(vertices_tmp.size())]   = bufferIndices[2*(vertices_tmp.size())-1];
  bufferIndices[2*(vertices_tmp.size())+1] = 2*vertices_tmp.size()+1;

  if(!face_one&&two_colour){
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    if((multicolour&&face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    Cartesian p = AB[0][j];
    p.normalize(p.length()-thickness);
    p1 = vertices_tmp[j] + p - AB[1][j];
    p2 = vertices_tmp[j] - p - AB[1][j];

    n =  -AB[1][j];
    n.normalize();

    bufferVertices[6*(vertices_tmp.size())+6*j]   = p1.get_x();
    bufferVertices[6*(vertices_tmp.size())+6*j+1] = p1.get_y();
    bufferVertices[6*(vertices_tmp.size())+6*j+2] = p1.get_z();
    bufferVertices[6*(vertices_tmp.size())+6*j+3] = p2.get_x();
    bufferVertices[6*(vertices_tmp.size())+6*j+4] = p2.get_y();
    bufferVertices[6*(vertices_tmp.size())+6*j+5] = p2.get_z();
    bufferNormals[6*(vertices_tmp.size())+6*j]   = n.get_x();
    bufferNormals[6*(vertices_tmp.size())+6*j+1] = n.get_y();
    bufferNormals[6*(vertices_tmp.size())+6*j+2] = n.get_z();
    bufferNormals[6*(vertices_tmp.size())+6*j+3] = n.get_x();
    bufferNormals[6*(vertices_tmp.size())+6*j+4] = n.get_y();
    bufferNormals[6*(vertices_tmp.size())+6*j+5] = n.get_z();
    bufferColors[8*(vertices_tmp.size())+8*j]     = colp[0];
    bufferColors[8*(vertices_tmp.size())+8*j+1]   = colp[1];
    bufferColors[8*(vertices_tmp.size())+8*j+2]   = colp[2];
    bufferColors[8*(vertices_tmp.size())+8*j+3]   = colp[3];
    bufferColors[8*(vertices_tmp.size())+8*j+4]   = colp[0];
    bufferColors[8*(vertices_tmp.size())+8*j+5]   = colp[1];
    bufferColors[8*(vertices_tmp.size())+8*j+6]   = colp[2];
    bufferColors[8*(vertices_tmp.size())+8*j+7]   = colp[3];

    bufferIndices[2*(vertices_tmp.size())+2+2*j+1]= 2*vertices_tmp.size()+2*j;
    bufferIndices[2*(vertices_tmp.size())+2+2*j]  = 2*vertices_tmp.size()+2*j+1;

  }

  bufferIndices[4*(vertices_tmp.size())+2] = bufferIndices[4*(vertices_tmp.size())+1];

  int vertexOffset = 12*(vertices_tmp.size());
  int colourOffset = 16*(vertices_tmp.size());
  int indexOffset = 4*(vertices_tmp.size());
  int indexOffsetInternal = 4*(vertices_tmp.size())+2+1;
  int ii=0;
  int maxInd = 0;
  for(int i=-90;i<=90;i=i+360/my_accu,ii+=1){
    double theta = (double)i/360.0 * PIBY2;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      // Indices for triangle strip
      if(ii!=(my_accu/2)){
        if(j==0){
          // Prepend a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii)*vertices_tmp.size();
          maxInd++;
        }
        bufferIndices[indexOffsetInternal++] = indexOffset + j + ii*vertices_tmp.size();
        maxInd++;
        bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii+1)*vertices_tmp.size();
        maxInd++;
        if(j==vertices_tmp.size()-1){
          // Append a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii+1)*vertices_tmp.size();
          maxInd++;
        }
      }
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p = AB[0][j];
      p.normalize(p.length()-thickness);
      Cartesian pr = AB[0][j];
      pr.normalize(thickness);
      Cartesian ppr = AB[1][j];
      ppr.normalize(thickness);
      Cartesian p1 = x*pr + y*ppr + vertices_tmp[j]+p;
      Cartesian n1 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j]   = p1.get_x();
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+1] = p1.get_y();
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+2] = p1.get_z();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j]    = n1.get_x();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+1]  = n1.get_y();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+2]  = n1.get_z();
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j]     = colp[0];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+1]   = colp[1];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+2]   = colp[2];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+3]   = colp[3];
    }
  }
 
  vertexOffset += 3*(ii*vertices_tmp.size());
  colourOffset += 4*(ii*vertices_tmp.size());
  indexOffset += (ii)*(vertices_tmp.size());

  ii=0;
  for(int i=90;i<=270;i=i+360/my_accu,ii+=1){
    double theta = (double)i/360.0 * PIBY2;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      // Indices for triangle strip
      if(ii!=(my_accu/2)){
        if(j==0){
          // Prepend a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii)*vertices_tmp.size();
          maxInd++;
        }
        bufferIndices[indexOffsetInternal++] = indexOffset + j + ii*vertices_tmp.size();
          maxInd++;
        bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii+1)*vertices_tmp.size();
          maxInd++;
        if(j==vertices_tmp.size()-1){
          // Append a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset + j + (ii+1)*vertices_tmp.size();
          maxInd++;
        }
      }
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p = AB[0][j];
      p.normalize(p.length()-thickness);
      Cartesian pr = AB[0][j];
      pr.normalize(thickness);
      Cartesian ppr = AB[1][j];
      ppr.normalize(thickness);
      Cartesian p1 = x*pr + y*ppr + vertices_tmp[j]-p;
      Cartesian n1 = pr*(x/(pr.length()*pr.length())) + ppr*(y/(ppr.length()*ppr.length()));
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j]   = p1.get_x();
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+1] = p1.get_y();
      bufferVertices[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+2] = p1.get_z();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j]    = n1.get_x();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+1]  = n1.get_y();
      bufferNormals[vertexOffset + 3*(ii*vertices_tmp.size())+3*j+2]  = n1.get_z();
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j]     = colp[0];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+1]   = colp[1];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+2]   = colp[2];
      bufferColors[colourOffset + 4*(ii*vertices_tmp.size())+4*j+3]   = colp[3];
    }
  }

  //nVerts = 8*(vertices_tmp.size()-1)+4*(vertices_tmp.size()-1)*my_accu;
  //polygonType = GL_QUADS;
  nVerts = 4*(vertices_tmp.size())+2 + 2*((vertices_tmp.size()+1)*my_accu);

  // Start cap

  vertexOffset += 3*(ii*vertices_tmp.size());
  colourOffset += 4*(ii*vertices_tmp.size());

  ii=0;
  Cartesian n1 = vertices_tmp[0] - vertices_tmp[1];
  for(int i=-90;i<=270;i=i+360/my_accu,ii++){
    double theta = (double)i/360.0 * PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    Cartesian p1 = x*AB[0][0] + y*AB[1][0] + vertices_tmp[0];
    bufferVertices[vertexOffset+3*ii]   = p1.get_x();
    bufferVertices[vertexOffset+3*ii+1] = p1.get_y();
    bufferVertices[vertexOffset+3*ii+2] = p1.get_z();
    bufferNormals[vertexOffset+3*ii]   = n1.get_x();
    bufferNormals[vertexOffset+3*ii+1] = n1.get_y();
    bufferNormals[vertexOffset+3*ii+2] = n1.get_z();
    bufferColors[colourOffset + 4*ii]     = colour_vector_tmp[0].get_x();
    bufferColors[colourOffset + 4*ii+1]   = colour_vector_tmp[0].get_y();
    bufferColors[colourOffset + 4*ii+2]   = colour_vector_tmp[0].get_z();
    bufferColors[colourOffset + 4*ii+3]   = colour_vector_tmp[0].get_a();
    if(ii<my_accu){
      bufferIndices[indexOffsetInternal++] = vertexOffset/3;
      bufferIndices[indexOffsetInternal++] = vertexOffset/3 + ii;
    } else {
      bufferIndices[indexOffsetInternal++] = vertexOffset/3 + (my_accu-ii);
      bufferIndices[indexOffsetInternal++] = vertexOffset/3;
    }
  }

  nVerts += 2*my_accu+2;

  // End cap

  vertexOffset += 3*(ii);
  colourOffset += 4*(ii);

  ii=0;
  n1 = vertices_tmp[vertices_tmp.size()-1] - vertices_tmp[vertices_tmp.size()-2];
  for(int i=-90;i<=270;i=i+360/my_accu,ii++){
    double theta = (double)i/360.0 * PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    Cartesian p1 = x*AB[0][vertices_tmp.size()-1] + y*AB[1][vertices_tmp.size()-1] + vertices_tmp[vertices_tmp.size()-1];
    bufferVertices[vertexOffset+3*ii]   = p1.get_x();
    bufferVertices[vertexOffset+3*ii+1] = p1.get_y();
    bufferVertices[vertexOffset+3*ii+2] = p1.get_z();
    bufferNormals[vertexOffset+3*ii]   = n1.get_x();
    bufferNormals[vertexOffset+3*ii+1] = n1.get_y();
    bufferNormals[vertexOffset+3*ii+2] = n1.get_z();
    bufferColors[colourOffset + 4*ii]     = colour_vector_tmp[0].get_x();
    bufferColors[colourOffset + 4*ii+1]   = colour_vector_tmp[0].get_y();
    bufferColors[colourOffset + 4*ii+2]   = colour_vector_tmp[0].get_z();
    bufferColors[colourOffset + 4*ii+3]   = colour_vector_tmp[0].get_a();
    if(ii<my_accu){
      bufferIndices[indexOffsetInternal++] = vertexOffset/3 + ii;
      bufferIndices[indexOffsetInternal++] = vertexOffset/3;
    } else {
      bufferIndices[indexOffsetInternal++] = vertexOffset/3;
      bufferIndices[indexOffsetInternal++] = vertexOffset/3 + (my_accu-ii);
    }
  }

  nVerts += 2*my_accu+2;

  polygonType = GL_TRIANGLE_STRIP;

}

void Ribbon::get_flat_ribbon_array(std::vector<Primitive*> &a,const std::vector<Cartesian> &vertices_tmp,const std::vector<std::vector<Cartesian> > &AB,int my_accu,int textured,int multicolour,const std::vector<Cartesian> &colour_vector_tmp, const double *colp, double alpha, const bool two_colour){
  Cartesian n;
  Cartesian p1,p2,p3,p4;
  unsigned int j;

  nVerts = vertices_tmp.size() * 8;
  //std::cout << "Ribbon::get_flat_ribbon_array " << vertices_tmp.size() << " " << my_accu << "\n";
  bufferVertices = new float[3*nVerts+24];
  bufferNormals = new float[3*nVerts+24];
  bufferColors = new float[4*nVerts+32];
  nRealVerts = 3*nVerts+24;
  bufferIndices = new GLuint[8*(vertices_tmp.size())+24];

  double col[4];
  int idx = vertices_tmp.size()/2;
  Cartesian p_t = vertices_tmp[idx+1] + vertices_tmp[idx-1] - 2*vertices_tmp[idx];
  p_t.normalize();
  Cartesian pvpr_t = AB[1][idx];
  pvpr_t.normalize();

  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;
  std::vector<Cartesian> v;
  std::vector<Cartesian> ns;
  double colour[4] = {0.6,0.6,0.6,1.0};

  if(face_one&&two_colour){
    colp=colour;
  }
  
  p1 = vertices_tmp[0] + AB[0][0] + AB[1][0];
  p2 = vertices_tmp[0] - AB[0][0] + AB[1][0];
  p3 = vertices_tmp[0] + AB[0][0] - AB[1][0];
  p4 = vertices_tmp[0] - AB[0][0] - AB[1][0];
  bufferVertices[24*vertices_tmp.size()]   = p1.get_x();
  bufferVertices[24*vertices_tmp.size()+1] = p1.get_y();
  bufferVertices[24*vertices_tmp.size()+2] = p1.get_z();
  bufferVertices[24*vertices_tmp.size()+3] = p2.get_x();
  bufferVertices[24*vertices_tmp.size()+4] = p2.get_y();
  bufferVertices[24*vertices_tmp.size()+5] = p2.get_z();
  bufferVertices[24*vertices_tmp.size()+6] = p3.get_x();
  bufferVertices[24*vertices_tmp.size()+7] = p3.get_y();
  bufferVertices[24*vertices_tmp.size()+8] = p3.get_z();
  bufferVertices[24*vertices_tmp.size()+9] = p4.get_x();
  bufferVertices[24*vertices_tmp.size()+10] = p4.get_y();
  bufferVertices[24*vertices_tmp.size()+11] = p4.get_z();

  p1 = vertices_tmp[vertices_tmp.size()-1] + AB[0][vertices_tmp.size()-1] + AB[1][vertices_tmp.size()-1];
  p2 = vertices_tmp[vertices_tmp.size()-1] - AB[0][vertices_tmp.size()-1] + AB[1][vertices_tmp.size()-1];
  p3 = vertices_tmp[vertices_tmp.size()-1] + AB[0][vertices_tmp.size()-1] - AB[1][vertices_tmp.size()-1];
  p4 = vertices_tmp[vertices_tmp.size()-1] - AB[0][vertices_tmp.size()-1] - AB[1][vertices_tmp.size()-1];
  bufferVertices[24*vertices_tmp.size()+12] = p1.get_x();
  bufferVertices[24*vertices_tmp.size()+13] = p1.get_y();
  bufferVertices[24*vertices_tmp.size()+14] = p1.get_z();
  bufferVertices[24*vertices_tmp.size()+15] = p2.get_x();
  bufferVertices[24*vertices_tmp.size()+16] = p2.get_y();
  bufferVertices[24*vertices_tmp.size()+17] = p2.get_z();
  bufferVertices[24*vertices_tmp.size()+18] = p3.get_x();
  bufferVertices[24*vertices_tmp.size()+19] = p3.get_y();
  bufferVertices[24*vertices_tmp.size()+20] = p3.get_z();
  bufferVertices[24*vertices_tmp.size()+21] = p4.get_x();
  bufferVertices[24*vertices_tmp.size()+22] = p4.get_y();
  bufferVertices[24*vertices_tmp.size()+23] = p4.get_z();

  n = vertices_tmp[0]-vertices_tmp[1];
  bufferNormals[24*vertices_tmp.size()]   = n.get_x();
  bufferNormals[24*vertices_tmp.size()+1] = n.get_y();
  bufferNormals[24*vertices_tmp.size()+2] = n.get_z();
  bufferNormals[24*vertices_tmp.size()+3] = n.get_x();
  bufferNormals[24*vertices_tmp.size()+4] = n.get_y();
  bufferNormals[24*vertices_tmp.size()+5] = n.get_z();
  bufferNormals[24*vertices_tmp.size()+6] = n.get_x();
  bufferNormals[24*vertices_tmp.size()+7] = n.get_y();
  bufferNormals[24*vertices_tmp.size()+8] = n.get_z();
  bufferNormals[24*vertices_tmp.size()+9] = n.get_x();
  bufferNormals[24*vertices_tmp.size()+10] = n.get_y();
  bufferNormals[24*vertices_tmp.size()+11] = n.get_z();

  n = vertices_tmp[vertices_tmp.size()-1]-vertices_tmp[vertices_tmp.size()-2];
  bufferNormals[24*vertices_tmp.size()+12] = n.get_x();
  bufferNormals[24*vertices_tmp.size()+13] = n.get_y();
  bufferNormals[24*vertices_tmp.size()+14] = n.get_z();
  bufferNormals[24*vertices_tmp.size()+15] = n.get_x();
  bufferNormals[24*vertices_tmp.size()+16] = n.get_y();
  bufferNormals[24*vertices_tmp.size()+17] = n.get_z();
  bufferNormals[24*vertices_tmp.size()+18] = n.get_x();
  bufferNormals[24*vertices_tmp.size()+19] = n.get_y();
  bufferNormals[24*vertices_tmp.size()+20] = n.get_z();
  bufferNormals[24*vertices_tmp.size()+21] = n.get_x();
  bufferNormals[24*vertices_tmp.size()+22] = n.get_y();
  bufferNormals[24*vertices_tmp.size()+23] = n.get_z();

  bufferColors[32*vertices_tmp.size()]   = colour_vector_tmp[0].get_x();
  bufferColors[32*vertices_tmp.size()+1] = colour_vector_tmp[0].get_y();
  bufferColors[32*vertices_tmp.size()+2] = colour_vector_tmp[0].get_z();
  bufferColors[32*vertices_tmp.size()+3] = colour_vector_tmp[0].get_a();
  bufferColors[32*vertices_tmp.size()+4] = colour_vector_tmp[0].get_x();
  bufferColors[32*vertices_tmp.size()+5] = colour_vector_tmp[0].get_y();
  bufferColors[32*vertices_tmp.size()+6] = colour_vector_tmp[0].get_z();
  bufferColors[32*vertices_tmp.size()+7] = colour_vector_tmp[0].get_a();
  bufferColors[32*vertices_tmp.size()+8] = colour_vector_tmp[0].get_x();
  bufferColors[32*vertices_tmp.size()+9] = colour_vector_tmp[0].get_y();
  bufferColors[32*vertices_tmp.size()+10] = colour_vector_tmp[0].get_z();
  bufferColors[32*vertices_tmp.size()+11] = colour_vector_tmp[0].get_a();
  bufferColors[32*vertices_tmp.size()+12] = colour_vector_tmp[0].get_x();
  bufferColors[32*vertices_tmp.size()+13] = colour_vector_tmp[0].get_y();
  bufferColors[32*vertices_tmp.size()+14] = colour_vector_tmp[0].get_z();
  bufferColors[32*vertices_tmp.size()+15] = colour_vector_tmp[0].get_a();

  bufferColors[32*vertices_tmp.size()+16] = colour_vector_tmp[colour_vector_tmp.size()-1].get_x();
  bufferColors[32*vertices_tmp.size()+17] = colour_vector_tmp[colour_vector_tmp.size()-1].get_y();
  bufferColors[32*vertices_tmp.size()+18] = colour_vector_tmp[colour_vector_tmp.size()-1].get_z();
  bufferColors[32*vertices_tmp.size()+19] = colour_vector_tmp[colour_vector_tmp.size()-1].get_a();
  bufferColors[32*vertices_tmp.size()+20] = colour_vector_tmp[colour_vector_tmp.size()-1].get_x();
  bufferColors[32*vertices_tmp.size()+21] = colour_vector_tmp[colour_vector_tmp.size()-1].get_y();
  bufferColors[32*vertices_tmp.size()+22] = colour_vector_tmp[colour_vector_tmp.size()-1].get_z();
  bufferColors[32*vertices_tmp.size()+23] = colour_vector_tmp[colour_vector_tmp.size()-1].get_a();
  bufferColors[32*vertices_tmp.size()+24] = colour_vector_tmp[colour_vector_tmp.size()-1].get_x();
  bufferColors[32*vertices_tmp.size()+25] = colour_vector_tmp[colour_vector_tmp.size()-1].get_y();
  bufferColors[32*vertices_tmp.size()+26] = colour_vector_tmp[colour_vector_tmp.size()-1].get_z();
  bufferColors[32*vertices_tmp.size()+27] = colour_vector_tmp[colour_vector_tmp.size()-1].get_a();
  bufferColors[32*vertices_tmp.size()+28] = colour_vector_tmp[colour_vector_tmp.size()-1].get_x();
  bufferColors[32*vertices_tmp.size()+29] = colour_vector_tmp[colour_vector_tmp.size()-1].get_y();
  bufferColors[32*vertices_tmp.size()+30] = colour_vector_tmp[colour_vector_tmp.size()-1].get_z();
  bufferColors[32*vertices_tmp.size()+31] = colour_vector_tmp[colour_vector_tmp.size()-1].get_a();

  for(j=0;j<vertices_tmp.size();j++){
    if((multicolour&&!face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    p1 = vertices_tmp[j] + AB[0][j] + AB[1][j];
    p2 = vertices_tmp[j] - AB[0][j] + AB[1][j];

    n =  AB[1][j];
    n.normalize();

    bufferVertices[6*j]   = p1.get_x();
    bufferVertices[6*j+1] = p1.get_y();
    bufferVertices[6*j+2] = p1.get_z();
    bufferVertices[6*j+3] = p2.get_x();
    bufferVertices[6*j+4] = p2.get_y();
    bufferVertices[6*j+5] = p2.get_z();
    bufferNormals[6*j]   = n.get_x();
    bufferNormals[6*j+1] = n.get_y();
    bufferNormals[6*j+2] = n.get_z();
    bufferNormals[6*j+3] = n.get_x();
    bufferNormals[6*j+4] = n.get_y();
    bufferNormals[6*j+5] = n.get_z();
    bufferColors[8*j]     = colp[0];
    bufferColors[8*j+1]   = colp[1];
    bufferColors[8*j+2]   = colp[2];
    bufferColors[8*j+3]   = colp[3];
    bufferColors[8*j+4]   = colp[0];
    bufferColors[8*j+5]   = colp[1];
    bufferColors[8*j+6]   = colp[2];
    bufferColors[8*j+7]   = colp[3];

    bufferIndices[2*j]   = 2*j;
    bufferIndices[2*j+1] = 2*j+1;

  }

  bufferIndices[2*(vertices_tmp.size())]   = bufferIndices[2*(vertices_tmp.size())-1];
  bufferIndices[2*(vertices_tmp.size())+1] = 2*vertices_tmp.size()+1;

  if(!face_one&&two_colour){
    double colour[4] = {0.6,0.6,0.6,1.0};
    colp=colour;
  }
  for(j=0;j<vertices_tmp.size();j++){
    if((multicolour&&face_one)||!two_colour){
      col[0] = colour_vector_tmp[j].get_x();
      col[1] = colour_vector_tmp[j].get_y();
      col[2] = colour_vector_tmp[j].get_z();
      col[3] = colour_vector_tmp[j].get_a();
      colp = col;
    }
    p1 = vertices_tmp[j] + AB[0][j] - AB[1][j];
    p2 = vertices_tmp[j] - AB[0][j] - AB[1][j];

    n =  -AB[1][j];
    n.normalize();

    bufferVertices[6*(vertices_tmp.size())+6*j]   = p1.get_x();
    bufferVertices[6*(vertices_tmp.size())+6*j+1] = p1.get_y();
    bufferVertices[6*(vertices_tmp.size())+6*j+2] = p1.get_z();
    bufferVertices[6*(vertices_tmp.size())+6*j+3] = p2.get_x();
    bufferVertices[6*(vertices_tmp.size())+6*j+4] = p2.get_y();
    bufferVertices[6*(vertices_tmp.size())+6*j+5] = p2.get_z();
    bufferNormals[6*(vertices_tmp.size())+6*j]   = n.get_x();
    bufferNormals[6*(vertices_tmp.size())+6*j+1] = n.get_y();
    bufferNormals[6*(vertices_tmp.size())+6*j+2] = n.get_z();
    bufferNormals[6*(vertices_tmp.size())+6*j+3] = n.get_x();
    bufferNormals[6*(vertices_tmp.size())+6*j+4] = n.get_y();
    bufferNormals[6*(vertices_tmp.size())+6*j+5] = n.get_z();
    bufferColors[8*(vertices_tmp.size())+8*j]     = colp[0];
    bufferColors[8*(vertices_tmp.size())+8*j+1]   = colp[1];
    bufferColors[8*(vertices_tmp.size())+8*j+2]   = colp[2];
    bufferColors[8*(vertices_tmp.size())+8*j+3]   = colp[3];
    bufferColors[8*(vertices_tmp.size())+8*j+4]   = colp[0];
    bufferColors[8*(vertices_tmp.size())+8*j+5]   = colp[1];
    bufferColors[8*(vertices_tmp.size())+8*j+6]   = colp[2];
    bufferColors[8*(vertices_tmp.size())+8*j+7]   = colp[3];

    bufferIndices[2*(vertices_tmp.size())+2+2*j+1]= 2*vertices_tmp.size()+2*j;
    bufferIndices[2*(vertices_tmp.size())+2+2*j]  = 2*vertices_tmp.size()+2*j+1;

  }

  bufferIndices[4*(vertices_tmp.size())+2] = bufferIndices[4*(vertices_tmp.size())+1];
  bufferIndices[4*(vertices_tmp.size())+3] = 4*vertices_tmp.size()+1;

  for(j=0;j<vertices_tmp.size();j++){
    col[0] = colour_vector_tmp[j].get_x();
    col[1] = colour_vector_tmp[j].get_y();
    col[2] = colour_vector_tmp[j].get_z();
    col[3] = colour_vector_tmp[j].get_a();
    colp = col;

    if(edges_grey)
      colp=colour;

    p1 = vertices_tmp[j] + AB[0][j] + AB[1][j];
    p2 = vertices_tmp[j] + AB[0][j] - AB[1][j];

    n =  AB[0][j];
    n.normalize();

    bufferVertices[12*(vertices_tmp.size())+6*j]   = p1.get_x();
    bufferVertices[12*(vertices_tmp.size())+6*j+1] = p1.get_y();
    bufferVertices[12*(vertices_tmp.size())+6*j+2] = p1.get_z();
    bufferVertices[12*(vertices_tmp.size())+6*j+3] = p2.get_x();
    bufferVertices[12*(vertices_tmp.size())+6*j+4] = p2.get_y();
    bufferVertices[12*(vertices_tmp.size())+6*j+5] = p2.get_z();
    bufferNormals[12*(vertices_tmp.size())+6*j]   = n.get_x();
    bufferNormals[12*(vertices_tmp.size())+6*j+1] = n.get_y();
    bufferNormals[12*(vertices_tmp.size())+6*j+2] = n.get_z();
    bufferNormals[12*(vertices_tmp.size())+6*j+3] = n.get_x();
    bufferNormals[12*(vertices_tmp.size())+6*j+4] = n.get_y();
    bufferNormals[12*(vertices_tmp.size())+6*j+5] = n.get_z();
    bufferColors[16*(vertices_tmp.size())+8*j]     = colp[0];
    bufferColors[16*(vertices_tmp.size())+8*j+1]   = colp[1];
    bufferColors[16*(vertices_tmp.size())+8*j+2]   = colp[2];
    bufferColors[16*(vertices_tmp.size())+8*j+3]   = colp[3];
    bufferColors[16*(vertices_tmp.size())+8*j+4]   = colp[0];
    bufferColors[16*(vertices_tmp.size())+8*j+5]   = colp[1];
    bufferColors[16*(vertices_tmp.size())+8*j+6]   = colp[2];
    bufferColors[16*(vertices_tmp.size())+8*j+7]   = colp[3];

    bufferIndices[4*(vertices_tmp.size())+4+2*j+1]= 4*vertices_tmp.size()+2*j;
    bufferIndices[4*(vertices_tmp.size())+4+2*j]  = 4*vertices_tmp.size()+2*j+1;

    p2 = vertices_tmp[j] - AB[0][j] + AB[1][j];
    p1 = vertices_tmp[j] - AB[0][j] - AB[1][j];

    n =  -AB[0][j];
    n.normalize();

    bufferVertices[18*(vertices_tmp.size())+6*j]   = p1.get_x();
    bufferVertices[18*(vertices_tmp.size())+6*j+1] = p1.get_y();
    bufferVertices[18*(vertices_tmp.size())+6*j+2] = p1.get_z();
    bufferVertices[18*(vertices_tmp.size())+6*j+3] = p2.get_x();
    bufferVertices[18*(vertices_tmp.size())+6*j+4] = p2.get_y();
    bufferVertices[18*(vertices_tmp.size())+6*j+5] = p2.get_z();
    bufferNormals[18*(vertices_tmp.size())+6*j]   = n.get_x();
    bufferNormals[18*(vertices_tmp.size())+6*j+1] = n.get_y();
    bufferNormals[18*(vertices_tmp.size())+6*j+2] = n.get_z();
    bufferNormals[18*(vertices_tmp.size())+6*j+3] = n.get_x();
    bufferNormals[18*(vertices_tmp.size())+6*j+4] = n.get_y();
    bufferNormals[18*(vertices_tmp.size())+6*j+5] = n.get_z();
    bufferColors[24*(vertices_tmp.size())+8*j]     = colp[0];
    bufferColors[24*(vertices_tmp.size())+8*j+1]   = colp[1];
    bufferColors[24*(vertices_tmp.size())+8*j+2]   = colp[2];
    bufferColors[24*(vertices_tmp.size())+8*j+3]   = colp[3];
    bufferColors[24*(vertices_tmp.size())+8*j+4]   = colp[0];
    bufferColors[24*(vertices_tmp.size())+8*j+5]   = colp[1];
    bufferColors[24*(vertices_tmp.size())+8*j+6]   = colp[2];
    bufferColors[24*(vertices_tmp.size())+8*j+7]   = colp[3];

    bufferIndices[6*(vertices_tmp.size())+6+2*j+1]= 6*vertices_tmp.size()+2*j;
    bufferIndices[6*(vertices_tmp.size())+6+2*j]  = 6*vertices_tmp.size()+2*j+1;

  }

  bufferIndices[6*(vertices_tmp.size())+4] = bufferIndices[6*(vertices_tmp.size())+3];
  bufferIndices[6*(vertices_tmp.size())+5] = 6*vertices_tmp.size()+1;

  bufferIndices[8*(vertices_tmp.size())+6] = bufferIndices[8*(vertices_tmp.size())+5];
  bufferIndices[8*(vertices_tmp.size())+7] = 8*vertices_tmp.size();

  bufferIndices[8*(vertices_tmp.size())+8] = 8*vertices_tmp.size();
  bufferIndices[8*(vertices_tmp.size())+9] = 8*vertices_tmp.size()+2;
  bufferIndices[8*(vertices_tmp.size())+10] = 8*vertices_tmp.size()+1;
  bufferIndices[8*(vertices_tmp.size())+11] = 8*vertices_tmp.size()+3;

  bufferIndices[8*(vertices_tmp.size())+12] = 8*vertices_tmp.size()+3;
  bufferIndices[8*(vertices_tmp.size())+13] = 8*vertices_tmp.size()+4;

  bufferIndices[8*(vertices_tmp.size())+14] = 8*vertices_tmp.size()+4;
  bufferIndices[8*(vertices_tmp.size())+15] = 8*vertices_tmp.size()+5;
  bufferIndices[8*(vertices_tmp.size())+16] = 8*vertices_tmp.size()+6;
  bufferIndices[8*(vertices_tmp.size())+17] = 8*vertices_tmp.size()+7;

  nVerts = 8*vertices_tmp.size()+18;
  polygonType = GL_TRIANGLE_STRIP;
}



void Ribbon::get_elliptical_ribbon_array(std::vector<Primitive*> &a,const std::vector<Cartesian> &vertices_tmp,const std::vector<std::vector<Cartesian> > &AB,int my_accu,int textured,int multicolour,const std::vector<Cartesian> &colour_vector_tmp, const double *colp, double alpha, const bool two_colour){

  if(vertices_tmp.size()<2||colour_vector_tmp.size()<2){
    nVerts = 0;
    return;
  }
  nVerts = vertices_tmp.size() * my_accu;
  //std::cout << "Ribbon::get_elliptical_ribbon_array " << vertices_tmp.size() << " " << my_accu << "\n";

  //useInterleaved = true;
  
  try {
    // +6 for 2 extra bits for centres of start and end caps
    if(useInterleaved) {
       bufferVNC = new SVertex[nVerts+2];
       nRealVerts = nVerts+2;
    } else {
       bufferVertices = new float[3*nVerts+6];
       bufferNormals = new float[3*nVerts+6];
       bufferColors = new float[4*nVerts+8];
       nRealVerts = 3*nVerts+6;
    }
    //std::cout << "size of bufferVertices "  << 3*nVerts+6 << "\n";
    //bufferIndices = new GLuint[6*(vertices_tmp.size()-1)*my_accu];
    // This does not allow for end caps, yet.
    // my_accu*4+8 extra indices for start and end cap triangle strips.
    bufferIndices = new GLuint[(2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*4+8];
  } catch (std::exception& e) {
    nVerts = 0;
    std::cerr << "Out of memory in get_elliptical_ribbon_array!!!\n";
    return;
  }

  // centres of start and end caps
  if(useInterleaved) {
    bufferVNC[nVerts].x   = vertices_tmp[0].get_x();
    bufferVNC[nVerts].y = vertices_tmp[0].get_y();
    bufferVNC[nVerts].z = vertices_tmp[0].get_z();
    bufferVNC[nVerts+1].x = vertices_tmp[vertices_tmp.size()-1].get_x();
    bufferVNC[nVerts+1].y = vertices_tmp[vertices_tmp.size()-1].get_y();
    bufferVNC[nVerts+1].z = vertices_tmp[vertices_tmp.size()-1].get_z();
    Cartesian n = vertices_tmp[0]-vertices_tmp[1];
    bufferVNC[nVerts].nx    = n.get_x();
    bufferVNC[nVerts].ny  = n.get_y();
    bufferVNC[nVerts].nz  = n.get_z();
    n = vertices_tmp[vertices_tmp.size()-1]-vertices_tmp[vertices_tmp.size()-2];
    bufferVNC[nVerts+1].nx  = n.get_x();
    bufferVNC[nVerts+1].ny  = n.get_y();
    bufferVNC[nVerts+1].nz  = n.get_z();
    bufferVNC[nVerts].r   = char(255*colour_vector_tmp[0].get_x());
    bufferVNC[nVerts].g   = char(255*colour_vector_tmp[0].get_y());
    bufferVNC[nVerts].r   = char(255*colour_vector_tmp[0].get_z());
    bufferVNC[nVerts].a   = char(255*colour_vector_tmp[0].get_a());
    bufferVNC[nVerts+1].r   = char(255*colour_vector_tmp[colour_vector_tmp.size()-1].get_x());
    bufferVNC[nVerts+1].g   = char(255*colour_vector_tmp[colour_vector_tmp.size()-1].get_y());
    bufferVNC[nVerts+1].b   = char(255*colour_vector_tmp[colour_vector_tmp.size()-1].get_z());
    bufferVNC[nVerts+1].a   = char(255*colour_vector_tmp[colour_vector_tmp.size()-1].get_a());
  } else {
    bufferVertices[3*nVerts]   = vertices_tmp[0].get_x();
    bufferVertices[3*nVerts+1] = vertices_tmp[0].get_y();
    bufferVertices[3*nVerts+2] = vertices_tmp[0].get_z();
    bufferVertices[3*nVerts+3] = vertices_tmp[vertices_tmp.size()-1].get_x();
    bufferVertices[3*nVerts+4] = vertices_tmp[vertices_tmp.size()-1].get_y();
    bufferVertices[3*nVerts+5] = vertices_tmp[vertices_tmp.size()-1].get_z();
    Cartesian n = vertices_tmp[0]-vertices_tmp[1];
    bufferNormals[3*nVerts]    = n.get_x();
    bufferNormals[3*nVerts+1]  = n.get_y();
    bufferNormals[3*nVerts+2]  = n.get_z();
    n = vertices_tmp[vertices_tmp.size()-1]-vertices_tmp[vertices_tmp.size()-2];
    bufferNormals[3*nVerts+3]  = n.get_x();
    bufferNormals[3*nVerts+4]  = n.get_y();
    bufferNormals[3*nVerts+5]  = n.get_z();
    bufferColors[4*nVerts]     = colour_vector_tmp[0].get_x();
    bufferColors[4*nVerts+1]   = colour_vector_tmp[0].get_y();
    bufferColors[4*nVerts+2]   = colour_vector_tmp[0].get_z();
    bufferColors[4*nVerts+3]   = colour_vector_tmp[0].get_a();
    bufferColors[4*nVerts+4]   = colour_vector_tmp[colour_vector_tmp.size()-1].get_x();
    bufferColors[4*nVerts+5]   = colour_vector_tmp[colour_vector_tmp.size()-1].get_y();
    bufferColors[4*nVerts+6]   = colour_vector_tmp[colour_vector_tmp.size()-1].get_z();
    bufferColors[4*nVerts+7]   = colour_vector_tmp[colour_vector_tmp.size()-1].get_a();
  }

  int idx = vertices_tmp.size()/2;
  Cartesian p_t = vertices_tmp[idx+1] + vertices_tmp[idx-1] - 2*vertices_tmp[idx];
  p_t.normalize();
  Cartesian pvpr_t = AB[1][idx];
  pvpr_t.normalize();
  bool face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;

  int ii=0;
  for(int i=0;i<360;i=i+360/my_accu,ii++){
    double theta = (double)i/360.0 * PIBY2;
    //bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2] = 3*nVerts;
    //bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-1] = ;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
        
      // /*
      // This allows us to merge my_accu separate strips into 1 strip using 2 degenerate triangles for each join.
      int degenOffset = ii * 2;
      // Indices for triangle strip
      if(ii<my_accu-1){
        bufferIndices[degenOffset + 2*ii*(vertices_tmp.size())+(j)*2+1] = j + (ii+1)*vertices_tmp.size();
        bufferIndices[degenOffset + 2*ii*(vertices_tmp.size())+(j)*2  ] = j + ii*vertices_tmp.size();
        if(j==0&&ii>0){
          // Prepend a degenerate triangle
          bufferIndices[degenOffset + 2*ii*(vertices_tmp.size())+(j)*2-1] = j + ii*vertices_tmp.size();
        }
        if(j==vertices_tmp.size()-1){
          // Append a degenerate triangle
          bufferIndices[degenOffset + 2*ii*(vertices_tmp.size())+(j)*2+2] = j + (ii+1)*vertices_tmp.size();
        }
      } else {
        if(j==0&&ii>0){
          // Prepend a degenerate triangle
          bufferIndices[degenOffset + 2*ii*(vertices_tmp.size())+(j)*2-1] = j + ii*vertices_tmp.size();
        }
        bufferIndices[degenOffset + 2*ii*(vertices_tmp.size())+(j)*2+1] = j;
        bufferIndices[degenOffset + 2*ii*(vertices_tmp.size())+(j)*2  ] = j + ii*vertices_tmp.size();
      }
      // */

      double x = cos(theta);
      double y = sin(theta);
      Cartesian p1 = x*AB[0][j] + y*AB[1][j] + vertices_tmp[j];
      Cartesian n1 = AB[0][j]*(x/(AB[0][j].length()*AB[0][j].length())) + AB[1][j]*(y/(AB[1][j].length()*AB[1][j].length()));
      if(useInterleaved) {
        bufferVNC[(ii*vertices_tmp.size())+j].x   = p1.get_x();
        bufferVNC[(ii*vertices_tmp.size())+j].y = p1.get_y();
        bufferVNC[(ii*vertices_tmp.size())+j].z = p1.get_z();
        bufferVNC[(ii*vertices_tmp.size())+j].nx    = n1.get_x();
        bufferVNC[(ii*vertices_tmp.size())+j].ny  = n1.get_y();
        bufferVNC[(ii*vertices_tmp.size())+j].nz  = n1.get_z();
        if(multicolour){
          if((!two_colour||(!face_one&&((i)>=0&&(i-360/my_accu)<180))||(face_one&&((i-360/my_accu)<0||(i)>=180)))){
            bufferVNC[(ii*vertices_tmp.size())+j].r   = char(255*colour_vector_tmp[j].get_x());
            bufferVNC[(ii*vertices_tmp.size())+j].g   = char(255*colour_vector_tmp[j].get_y());
            bufferVNC[(ii*vertices_tmp.size())+j].b   = char(255*colour_vector_tmp[j].get_z());
            bufferVNC[(ii*vertices_tmp.size())+j].a   = char(255*colour_vector_tmp[j].get_a());
          } else {
            bufferVNC[(ii*vertices_tmp.size())+j].r   = char(255*0.6);
            bufferVNC[(ii*vertices_tmp.size())+j].g   = char(255*0.6);
            bufferVNC[(ii*vertices_tmp.size())+j].b   = char(255*0.6);
            bufferVNC[(ii*vertices_tmp.size())+j].a   = char(255*colour_vector_tmp[j].get_a());
          }
        } // else ???
      }else{
        bufferVertices[3*(ii*vertices_tmp.size())+3*j]   = p1.get_x();
        bufferVertices[3*(ii*vertices_tmp.size())+3*j+1] = p1.get_y();
        bufferVertices[3*(ii*vertices_tmp.size())+3*j+2] = p1.get_z();
        bufferNormals[3*(ii*vertices_tmp.size())+3*j]    = n1.get_x();
        bufferNormals[3*(ii*vertices_tmp.size())+3*j+1]  = n1.get_y();
        bufferNormals[3*(ii*vertices_tmp.size())+3*j+2]  = n1.get_z();
        if(multicolour){
          if((!two_colour||(!face_one&&((i)>=0&&(i-360/my_accu)<180))||(face_one&&((i-360/my_accu)<0||(i)>=180)))){
            bufferColors[4*(ii*vertices_tmp.size())+4*j]     = colour_vector_tmp[j].get_x();
            bufferColors[4*(ii*vertices_tmp.size())+4*j+1]   = colour_vector_tmp[j].get_y();
            bufferColors[4*(ii*vertices_tmp.size())+4*j+2]   = colour_vector_tmp[j].get_z();
            bufferColors[4*(ii*vertices_tmp.size())+4*j+3]   = colour_vector_tmp[j].get_a();
          } else {
            bufferColors[4*(ii*vertices_tmp.size())+4*j]     = 0.6;
            bufferColors[4*(ii*vertices_tmp.size())+4*j+1]   = 0.6;
            bufferColors[4*(ii*vertices_tmp.size())+4*j+2]   = 0.6;
            bufferColors[4*(ii*vertices_tmp.size())+4*j+3]   = colour_vector_tmp[j].get_a();
          }
        } // else ???
      }
    }
  }
  //nVerts = 6*(vertices_tmp.size()-1)*my_accu;
  //polygonType = GL_TRIANGLES;

  // Some degenerate triangles
  bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2] = bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2-1];
  bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-1] = 0;
  // Start cap
  for(ii=0;ii<my_accu+1;ii++){
    if(ii<my_accu){
      bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)+2*ii]   = ii*vertices_tmp.size();
      bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)+2*ii+1] = nVerts;
    } else {
      bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)+2*ii]   = (ii-my_accu)*vertices_tmp.size();
      bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)+2*ii+1] = nVerts;
    }
  }

  // Some degenerate triangles
  bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*2+4] = bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*2+4-1];
  bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*2+4+1] = nVerts+1;
  // End cap
  for(ii=0;ii<my_accu+1;ii++){
    if(ii<my_accu){
      bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*2+6+2*ii+1]   = vertices_tmp.size()-1+ii*vertices_tmp.size();
      bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*2+6+2*ii] = nVerts+1;
    } else {
      bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*2+6+2*ii+1]   = vertices_tmp.size()-1+(ii-my_accu)*vertices_tmp.size();
      bufferIndices[(2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*2+6+2*ii] = nVerts+1;
    }
  }

  nVerts = (2*vertices_tmp.size()+2)*(my_accu)-2+my_accu*2+8  + 2*my_accu;
  polygonType = GL_TRIANGLE_STRIP;

  for(unsigned ii=0; ii < nVerts; ii++){
     if(bufferIndices[ii]>nRealVerts){
        std::cout << ii << " " << bufferIndices[ii] << "\n";
     }
  }
}

void get_elliptical_ribbon(std::vector<Primitive*> &a,const std::vector<Cartesian> &vertices_tmp,const std::vector<std::vector<Cartesian> > &AB,int my_accu,int textured,int multicolour,const std::vector<Cartesian> &colour_vector_tmp, const double *colp, double alpha){
  double col[4];
  for(int i=0;i<360;i=i+360/my_accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/my_accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
    std::vector<Cartesian> ns;
    for(unsigned int j=0;j<vertices_tmp.size();j++){
      if(multicolour){
        col[0] = colour_vector_tmp[j].get_x();
        col[1] = colour_vector_tmp[j].get_y();
        col[2] = colour_vector_tmp[j].get_z();
        col[3] = colour_vector_tmp[j].get_a();
        colp = col;
      }
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p1 = x*AB[0][j] + y*AB[1][j] + vertices_tmp[j];
      Cartesian n1 = AB[0][j]*(x/(AB[0][j].length()*AB[0][j].length())) + AB[1][j]*(y/(AB[1][j].length()*AB[1][j].length()));
      x = cos(theta2);
      y = sin(theta2);
      Cartesian p2 = x*AB[0][j] + y*AB[1][j] + vertices_tmp[j];
      Cartesian n2 = AB[0][j]*(x/(AB[0][j].length()*AB[0][j].length())) + AB[1][j]*(y/(AB[1][j].length()*AB[1][j].length()));
      if(v.size()==0){
        v.push_back(p2);
        ns.push_back(n2);
        v.push_back(p1);
        ns.push_back(n1);
      } else {
        v.push_back(p1);
        ns.push_back(n1);
        v.push_back(p2);
        ns.push_back(n2);
      }
      if(v.size()==4){
        QuadElement *q = new QuadElement(v,colp,0.25*(v[0]+v[1]+v[2]+v[3]),alpha,textured);
        q->SetNormals(ns);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        std::vector<Cartesian> nsnew;
        vnew.push_back(v[3]);
        vnew.push_back(v[2]);
        nsnew.push_back(ns[3]);
        nsnew.push_back(ns[2]);
        v = vnew;
        ns = nsnew;
      }
    }
  }
}

void Primitive::generateArrays(){

	freeResources();
	
	std::vector<Primitive*> quad_prims = GetSimplePrimitivesArrays(Volume(),matrix(),Cartesian());
        //std::cout << "quad_prims.size() " << quad_prims.size() << "\n";
        if(quad_prims.size()>0){
        std::cout <<"quad_prims.size(): " << quad_prims.size() << "\n";
        // This happens if we do not have shared-vertex storage.
        // If GetSimplePrimitivesArrays used shared-vertex storage method, then it
        // sets buffers.
	nVerts = quad_prims.size()*4;
	
	bufferVertices = new float[3*nVerts];
	bufferNormals = new float[3*nVerts];
	bufferColors = new float[4*nVerts];
	//std::cout << "Array size " << 3*nVerts << " " << 4*nVerts << "\n";
        //std::cout << "OPAQUENESS: " << opaqueness << "\n";
	for (unsigned i=0; i<quad_prims.size(); i++){
		//Copy vertex into vertices array	
		const std::vector<Cartesian> &verts = quad_prims[i]->GetVertices();
		const std::vector<Cartesian> &norms = quad_prims[i]->GetNormals();
		const double* col = quad_prims[i]->GetColour();
		for (int j=0; j<4; j++){
			bufferVertices[12*i+j*3   ] = verts[j].get_x();
			bufferVertices[12*i+j*3 +1] = verts[j].get_y();
			bufferVertices[12*i+j*3 +2] = verts[j].get_z();
		}
		for (int j=0; j<4; j++){
			bufferNormals[12*i+j*3   ] = norms[j].get_x();
			bufferNormals[12*i+j*3 +1] = norms[j].get_y();
			bufferNormals[12*i+j*3 +2] = norms[j].get_z();
		}
		for (int j=0; j<4; j++){
			bufferColors[16*i+j*4   ] = col[0];
			bufferColors[16*i+j*4 +1] = col[1];
			bufferColors[16*i+j*4 +2] = col[2];
			bufferColors[16*i+j*4 +3] = quad_prims[i]->GetAlpha();
		}
	}
	//std::cout << "Max array size " << quad_prims.size()*12 << " " <<  quad_prims.size()*16 << "\n";
	
	bufferIndices = new GLuint[nVerts];
	int idx=0;
	for (int i=0; i<nVerts; i++){
		bufferIndices[idx] = idx;
		idx++;
	}

	std::vector<Primitive*>::iterator prim_iter = quad_prims.begin();
	while(prim_iter!=quad_prims.end()){
		delete *prim_iter;
		prim_iter++;
	}
	quad_prims.clear();
        }

	arraysGenerated = 1;
	arraysUploaded = 0;

}

std::vector<Primitive*> Ribbon::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  std::vector<Primitive*> a;
  if(vertices.size()<2) return a;
  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;
  const double *colp = GetColour();
  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);

  std::vector<Cartesian> vertices_tmp = vertices;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;
  if(RenderQuality::GetSmoothRibbons()) {
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    colour_vector_tmp.resize(colour_vector.size()*4);
    for(unsigned icol=0;icol<colour_vector.size();icol++){
      colour_vector_tmp[icol*4] = colour_vector[icol];
      colour_vector_tmp[icol*4+1] = colour_vector[icol];
      colour_vector_tmp[icol*4+2] = colour_vector[icol];
      colour_vector_tmp[icol*4+3] = colour_vector[icol];
    }
  }

  int my_accu = accu;
  if(accu<4) my_accu = 8;

  if(RenderQuality::GetRenderQuality()&&0)
    my_accu = 180;
  if(style==0){
    get_elliptical_ribbon_array(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,two_colour);
  }else if(style==1){
    get_flat_ribbon_array(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,two_colour);
  }else if(style==2){
    get_flat_rounded_ribbon_array(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,two_colour);
  }else if(style==3){
    get_fancy_ribbon_array(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,two_colour);
  }

  return a;
}

std::vector<Primitive*> Primitive::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  std::vector<Primitive*> a;
  return a;
};

// FIXME Not dashed.
std::vector<Primitive*> DashCylinderCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  //std::cout << "DashCylinderCollection::GetSimplePrimitivesArrays\n";
  return CylinderCollection::GetSimplePrimitivesArrays(clip_vol, objrotmatrix, objorigin, start, end);
}

// FIXME ?? Seems to work in render with this. Render seems to not get here.
std::vector<Primitive*> DashCylinderCollection::GetPrimitives() const {
  //std::cout << "DashCylinderCollection::GetPrimitives\n";
  return CylinderCollection::GetPrimitives();
}

std::vector<Primitive*> CylinderCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {

  //std::cout << "CylinderCollection::GetSimplePrimitivesArrays\n";
  std::vector<Primitive*> a;

  nRealVerts = 0;
  nVerts = 0;

  unsigned nCylinders = prims.size();
  if(nCylinders<1) return a;

  //nCylinders = 20;


  int my_accu = 8;
  if(cylinder_accu>2) my_accu = cylinder_accu;

  bufferVertices = new float[6*my_accu*nCylinders];
  bufferNormals = new float[6*my_accu*nCylinders];
  bufferColors = new float[8*my_accu*nCylinders];
  bufferIndices = new GLuint[2*(my_accu+1)*nCylinders + 2*nCylinders];

  unsigned indexOffsetInternal = 0;

  //std::cout << "cylinder_accu: " << cylinder_accu << "\n"; std::cout.flush();
  //std::cout << "accu: " << my_accu << "\n"; std::cout.flush();
  //std::cout << "nCylinders: " << nCylinders << "\n"; std::cout.flush();

  int maxInd = 0;

  Cartesian xaxis(1,0,0);
  Cartesian yaxis(0,1,0);
  Cartesian zaxis(0,0,1);
  Cartesian up;
  Cartesian right;

  for(unsigned k=0;k<nCylinders&&k<prims.size();k++){
    //std::cout << "k: " << k << "\n"; std::cout.flush();
    int jj=0;
    std::vector<Cartesian> carts = prims[k]->GetVertices();
    Cartesian p = carts[1] - carts[0];
    p.normalize();
    if(fabs(Cartesian::DotProduct(xaxis,p))<0.95){
      up = Cartesian::CrossProduct(xaxis,p);
      up.normalize();
    }else if(abs(Cartesian::DotProduct(yaxis,p))<0.95){
      up = Cartesian::CrossProduct(yaxis,p);
      up.normalize();
    } else {
      up = Cartesian::CrossProduct(zaxis,p);
      up.normalize();
    }
    right = Cartesian::CrossProduct(up,p);
    up.normalize(prims[k]->GetSize());
    right.normalize(prims[k]->GetSize());
    //std::cout << carts[0] << " " << carts[1] << "\n";
    //std::cout << p << "\n";
    //std::cout << up << "\n";
    //std::cout << right << "\n";
    int indexOffset = 2*k*(my_accu);
    //std::cout << "indexOffset: " << " " << indexOffset << "\n";
    for(int j=0;j<=360;j=j+360/my_accu,jj++){
      //std::cout << "j: " << j << "\n"; std::cout.flush();
      //std::cout << "jj: " << jj << "\n"; std::cout.flush();
      if(indexOffsetInternal>=2*(my_accu+1)*nCylinders + 2*nCylinders-2){
         std::cerr << "PANIC, index overflow in drawing cylinders\n";
         std::cerr << indexOffsetInternal << "\n";
         std::cerr << 2*my_accu*nCylinders + 2*nCylinders << "\n";
         return a;
      }
      if(jj==my_accu){
        bufferIndices[indexOffsetInternal++] = indexOffset + 1;
        bufferIndices[indexOffsetInternal++] = indexOffset;
        // Append a degenerate triangle
        bufferIndices[indexOffsetInternal++] = indexOffset;
      } else {
        if(j==0){
          // Prepend a degenerate triangle
          bufferIndices[indexOffsetInternal++] = indexOffset+1;
        }
        bufferIndices[indexOffsetInternal++] = indexOffset + 2*jj+1;
        bufferIndices[indexOffsetInternal++] = indexOffset + 2*jj;
      }
      if(jj<my_accu){
      double theta = (double)j/360.0 * PIBY2;
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p1 = x*up + y*right + carts[0];
      Cartesian p2 = x*up + y*right + carts[1];
      Cartesian n1 = up*(x/(up.length()*up.length())) + right*(y/(right.length()*right.length()));
      bufferVertices[6*k*my_accu+6*jj]   = p1.get_x();
      bufferVertices[6*k*my_accu+6*jj+1] = p1.get_y();
      bufferVertices[6*k*my_accu+6*jj+2] = p1.get_z();
      bufferVertices[6*k*my_accu+6*jj+3] = p2.get_x();
      bufferVertices[6*k*my_accu+6*jj+4] = p2.get_y();
      bufferVertices[6*k*my_accu+6*jj+5] = p2.get_z();
      bufferNormals[6*k*my_accu+6*jj]   = n1.get_x();
      bufferNormals[6*k*my_accu+6*jj+1] = n1.get_y();
      bufferNormals[6*k*my_accu+6*jj+2] = n1.get_z();
      bufferNormals[6*k*my_accu+6*jj+3] = n1.get_x();
      bufferNormals[6*k*my_accu+6*jj+4] = n1.get_y();
      bufferNormals[6*k*my_accu+6*jj+5] = n1.get_z();
      maxInd = 6*k*my_accu+6*jj+6;
      bufferColors[8*k*my_accu+8*jj  ] = prims[k]->GetColour()[0];
      bufferColors[8*k*my_accu+8*jj+1] = prims[k]->GetColour()[1];
      bufferColors[8*k*my_accu+8*jj+2] = prims[k]->GetColour()[2];
      bufferColors[8*k*my_accu+8*jj+3] = alpha;
      bufferColors[8*k*my_accu+8*jj+4] = prims[k]->GetColour()[0];
      bufferColors[8*k*my_accu+8*jj+5] = prims[k]->GetColour()[1];
      bufferColors[8*k*my_accu+8*jj+6] = prims[k]->GetColour()[2];
      bufferColors[8*k*my_accu+8*jj+7] = alpha;
      }
    }
  }

  nRealVerts = 6*my_accu*nCylinders;
  nVerts = 2*(my_accu+1)*nCylinders + 2*nCylinders;

  //std::cout << indexOffsetInternal << " " << nVerts << "\n";
  //std::cout << maxInd << " " << nRealVerts << "\n";
  polygonType = GL_TRIANGLE_STRIP;

  return a;
};

std::vector<Primitive*> ImposterSphereCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  std::vector<Primitive*> a;
  nVerts = 0;
  nRealVerts = 0;
  if(prims.size()<1) return a;

  bufferVertices = new float[12*prims.size()];
  bufferOrigins = new float[12*prims.size()];
  bufferSizes = new float[4*prims.size()];
  bufferNormals = new float[12*prims.size()];
  bufferColors = new float[16*prims.size()];
  bufferIndices = new GLuint[4*prims.size()];

  Cartesian x_trans = Cartesian(1,0,0);
  Cartesian y_trans = Cartesian(0,1,0);
  Cartesian z_trans = Cartesian(0,0,1);
  x_trans.normalize();
  y_trans.normalize();
  z_trans.normalize();

  for(unsigned i=0;i<prims.size();i++){
    double size = prims[i]->GetSize();
    bufferSizes[4*i]   = size;
    bufferSizes[4*i+1] = size;
    bufferSizes[4*i+2] = size;
    bufferSizes[4*i+3] = size;
    bufferIndices[4*i]   = 4*i;
    bufferIndices[4*i+1] = 4*i+1;
    bufferIndices[4*i+2] = 4*i+2;
    bufferIndices[4*i+3] = 4*i+3;
    bufferVertices[12*i]   = -1;
    bufferVertices[12*i+1] = -1;
    bufferVertices[12*i+2] = 0.0;
    bufferVertices[12*i+3] = 1;
    bufferVertices[12*i+4] = -1;
    bufferVertices[12*i+5] = 0.0;
    bufferVertices[12*i+6] = 1;
    bufferVertices[12*i+7] = 1;
    bufferVertices[12*i+8] = 0.0;
    bufferVertices[12*i+9] = -1;
    bufferVertices[12*i+10] = 1;
    bufferVertices[12*i+11] = 0.0;
    bufferOrigins[12*i]   = prims[i]->GetVertices()[0].get_x();
    bufferOrigins[12*i+1] = prims[i]->GetVertices()[0].get_y();
    bufferOrigins[12*i+2] = prims[i]->GetVertices()[0].get_z();
    bufferOrigins[12*i+3] = prims[i]->GetVertices()[0].get_x();
    bufferOrigins[12*i+4] = prims[i]->GetVertices()[0].get_y();
    bufferOrigins[12*i+5] = prims[i]->GetVertices()[0].get_z();
    bufferOrigins[12*i+6] = prims[i]->GetVertices()[0].get_x();
    bufferOrigins[12*i+7] = prims[i]->GetVertices()[0].get_y();
    bufferOrigins[12*i+8] = prims[i]->GetVertices()[0].get_z();
    bufferOrigins[12*i+9] = prims[i]->GetVertices()[0].get_x();
    bufferOrigins[12*i+10] = prims[i]->GetVertices()[0].get_y();
    bufferOrigins[12*i+11] = prims[i]->GetVertices()[0].get_z();
    bufferNormals[12*i]   = 0.0;
    bufferNormals[12*i+1] = 0.0;
    bufferNormals[12*i+2] = 1.0;
    bufferNormals[12*i+3] = 0.0;
    bufferNormals[12*i+4] = 0.0;
    bufferNormals[12*i+5] = 1.0;
    bufferNormals[12*i+6] = 0.0;
    bufferNormals[12*i+7] = 0.0;
    bufferNormals[12*i+8] = 1.0;
    bufferNormals[12*i+9] = 0.0;
    bufferNormals[12*i+10] = 0.0;
    bufferNormals[12*i+11] = 1.0;
    bufferColors[16*i]   = prims[i]->GetColour()[0];
    bufferColors[16*i+1] = prims[i]->GetColour()[1];
    bufferColors[16*i+2] = prims[i]->GetColour()[2];
    bufferColors[16*i+3] = prims[i]->GetAlpha();
    bufferColors[16*i+4] = prims[i]->GetColour()[0];
    bufferColors[16*i+5] = prims[i]->GetColour()[1];
    bufferColors[16*i+6] = prims[i]->GetColour()[2];
    bufferColors[16*i+7] = prims[i]->GetAlpha();
    bufferColors[16*i+8] = prims[i]->GetColour()[0];
    bufferColors[16*i+9] = prims[i]->GetColour()[1];
    bufferColors[16*i+10] = prims[i]->GetColour()[2];
    bufferColors[16*i+11] = prims[i]->GetAlpha();
    bufferColors[16*i+12] = prims[i]->GetColour()[0];
    bufferColors[16*i+13] = prims[i]->GetColour()[1];
    bufferColors[16*i+14] = prims[i]->GetColour()[2];
    bufferColors[16*i+15] = prims[i]->GetAlpha();
  }
  nRealVerts = 12*prims.size();
  nVerts = 4*prims.size();

  polygonType = GL_QUADS;
  return a;
}

std::vector<Primitive*> TriangleCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  //std::cout << "TriangleCollection::GetSimplePrimitivesArrays\n"; std::cout.flush();
  //std::cout << prims.size() << "\n"; std::cout.flush();
  std::vector<Primitive*> a;
  nVerts = 0;
  nRealVerts = 0;
  if(prims.size()<1) return a;

  bufferVertices = new float[9*prims.size()];
  bufferNormals = new float[9*prims.size()];
  bufferColors = new float[12*prims.size()];
  bufferIndices = new GLuint[3*prims.size()];
  nRealVerts = 9*prims.size();
  nVerts = 3*prims.size();

  for(unsigned i=0;i<prims.size();i++){
    bufferIndices[3*i] = 3*i;
    bufferIndices[3*i+1] = 3*i+1;
    bufferIndices[3*i+2] = 3*i+2;
    bufferVertices[9*i]   = prims[i]->GetVertices()[0].get_x();
    bufferVertices[9*i+1] = prims[i]->GetVertices()[0].get_y();
    bufferVertices[9*i+2] = prims[i]->GetVertices()[0].get_z();
    bufferVertices[9*i+3] = prims[i]->GetVertices()[1].get_x();
    bufferVertices[9*i+4] = prims[i]->GetVertices()[1].get_y();
    bufferVertices[9*i+5] = prims[i]->GetVertices()[1].get_z();
    bufferVertices[9*i+6] = prims[i]->GetVertices()[2].get_x();
    bufferVertices[9*i+7] = prims[i]->GetVertices()[2].get_y();
    bufferVertices[9*i+8] = prims[i]->GetVertices()[2].get_z();
    bufferNormals[9*i]   = prims[i]->GetNormals()[0].get_x();
    bufferNormals[9*i+1] = prims[i]->GetNormals()[0].get_y();
    bufferNormals[9*i+2] = prims[i]->GetNormals()[0].get_z();
    bufferNormals[9*i+3] = prims[i]->GetNormals()[1].get_x();
    bufferNormals[9*i+4] = prims[i]->GetNormals()[1].get_y();
    bufferNormals[9*i+5] = prims[i]->GetNormals()[1].get_z();
    bufferNormals[9*i+6] = prims[i]->GetNormals()[2].get_x();
    bufferNormals[9*i+7] = prims[i]->GetNormals()[2].get_y();
    bufferNormals[9*i+8] = prims[i]->GetNormals()[2].get_z();
    bufferColors[12*i]   = prims[i]->GetColour()[0];
    bufferColors[12*i+1] = prims[i]->GetColour()[1];
    bufferColors[12*i+2] = prims[i]->GetColour()[2];
    bufferColors[12*i+3] = prims[i]->GetAlpha();
    bufferColors[12*i+4] = prims[i]->GetColour()[0];
    bufferColors[12*i+5] = prims[i]->GetColour()[1];
    bufferColors[12*i+6] = prims[i]->GetColour()[2];
    bufferColors[12*i+7] = prims[i]->GetAlpha();
    bufferColors[12*i+8] = prims[i]->GetColour()[0];
    bufferColors[12*i+9] = prims[i]->GetColour()[1];
    bufferColors[12*i+10] = prims[i]->GetColour()[2];
    bufferColors[12*i+11] = prims[i]->GetAlpha();
  }

  polygonType = GL_TRIANGLES;
  return a;
};

std::vector<Primitive*> SpheroidCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
   std::vector<Primitive*> a;

  //std::cout << "SphereCollection::GetSimplePrimitivesArrays\n"; std::cout.flush();

  unsigned nSpheres = prims.size();

  if(nSpheres==0) return a;

  // /*
  //nSpheres = 30;

  //std::cout << "nSpheres: " << nSpheres << "\n"; std::cout.flush();

  static int spheres_accus[] = {8,12,20};


  int my_accu = 8;
  if(sphere_accu>=0&&sphere_accu<3) my_accu = spheres_accus[sphere_accu];

  //std::cout << "accu:" << sphere_accu << "\n";
  //std::cout << "my_accu:" << my_accu << "\n";
  //std::cout << "alpha:" << GetAlpha() << "\n";

  int show_axes = ((SpheroidElement*)(prims[0]))->GetDrawAxes();
  int show_solid = ((SpheroidElement*)(prims[0]))->GetDrawSolid();

  nRealVerts = 0;
  nVerts = 0;

  // + 3*3*my_accu*2*nSpheres for the 3 triangle strips of the axes representations.
  int axesVertexOffset = 0;
  int axesColorOffset = 0;
  if(show_solid){
    nRealVerts += 3*(my_accu+1)*my_accu*nSpheres;
    axesVertexOffset = 3*(my_accu+1)*my_accu*nSpheres;
    axesColorOffset = 4*(my_accu+1)*my_accu*nSpheres;
  }
  if(show_axes) nRealVerts += 3*3*my_accu*2*nSpheres;
  bufferVertices = new float[nRealVerts];
  bufferNormals = new float[nRealVerts];
  bufferColors = new float[nRealVerts*4/3];

  int axesOffset = 0;
  int axesOffsetInternal = 0;
  if(show_solid) {
    nVerts += (my_accu)*my_accu*nSpheres + 2*my_accu*nSpheres;
    axesOffsetInternal = (my_accu)*my_accu*nSpheres + 2*my_accu*nSpheres;
    axesOffset = (my_accu+1)*my_accu*nSpheres;
  }
  if(show_axes) nVerts += 3*(my_accu+2)*2*nSpheres ;

  bufferIndices = new GLuint[nVerts];

  int indexOffsetInternal = 0;

  float strip_offset = 1.0001;
  float strip_width = .03;

  int maxInd = 0;
  matrix m;
  float lx;
  float ly;
  float lz;

  int offset = axesOffsetInternal;

  if(show_axes){
    for(unsigned k=0;k<nSpheres&&k<prims.size();k++){
      for(int p=0;p<3;p++){
        int jj = 0;
        for(int j=0;j<=360;j=j+360/my_accu,jj++){
          if(j==0){
            bufferIndices[offset++] = axesOffset+3*k*(my_accu)*2+p*(my_accu)*2;
          }
          if(jj==my_accu){
            bufferIndices[offset++] = axesOffset+3*k*(my_accu)*2+p*(my_accu)*2;
            bufferIndices[offset++] = axesOffset+3*k*(my_accu)*2+p*(my_accu)*2+1;
            bufferIndices[offset++] = axesOffset+3*k*(my_accu)*2+p*(my_accu)*2+1;
          } else {
            bufferIndices[offset++] = axesOffset+3*k*(my_accu)*2+p*(my_accu)*2+jj*2;
            bufferIndices[offset++] = axesOffset+3*k*(my_accu)*2+p*(my_accu)*2+jj*2+1;
          }
        }
      }
    }
  }
  //std::cout << "nSpheres: " << nSpheres << "\n"; std::cout.flush();
  //std::cout << "my_accu: " << my_accu << "\n"; std::cout.flush();
  //std::cout << "offset,calc offset: " << offset << " " << 3*(my_accu+2)*2*nSpheres << "\n";
  //std::cout << nVerts << "\n";

  for(unsigned k=0;k<nSpheres&&k<prims.size();k++){
    SpheroidElement* p = (SpheroidElement*)(prims[k]);
    if(p->GetMatrix().get_rows()==4 && p->GetMatrix().get_columns()==4){
      m = p->GetMatrix();
      lx = (m*Cartesian(1,0,0)).length();
      ly = (m*Cartesian(0,1,0)).length();
      lz = (m*Cartesian(0,0,1)).length();
    } else {
      Cartesian cartx = Cartesian(1,0,0);
      Cartesian carty = Cartesian(0,1,0);
      Cartesian yprim;
      Quat q1,q2;
      m = matrix(4,kdelta);

      double angle = acos(cartx.DotProduct(cartx,p->GetXAxis()))*180.0/M_PI;
      //std::cout << "Angle 1: " << angle << "\n";
      if(fabs(angle)>0.00001){
        Cartesian rotax = cartx.CrossProduct(cartx,p->GetXAxis());
        q1 = Quat(rotax,1,-angle);
      }
      m = q1.getMatrix();

      yprim = m * carty;
      angle = acos(yprim.DotProduct(yprim,p->GetYAxis()))*180.0/M_PI;

      if(fabs(angle)>0.00001){
        q2 = Quat(cartx,1,-angle);
        q1.postMult(q2);
      }
      m = q1.getMatrix();
    }

    matrix mTrans = m.Transpose();
    std::vector<double> mTransVector = mTrans.matrixVector();
    bool doMatrix = false;
    if(p->GetMatrix().get_rows()==4 && p->GetMatrix().get_columns()==4){
      doMatrix = true;
    }

    if(show_axes){
      int jj = 0;
      for(int j=0;j<=360;j=j+360/my_accu,jj++){
        if(jj<my_accu){
        float theta = (double)j/360.0 * PIBY2;
        float xs1;
        float ys1;
        float zs1;
        float xs2;
        float ys2;
        float zs2;
        float x1,x2;
        float y1,y2;
        float z1;
        float z2;
        x1 = x2 = strip_offset*cos(theta);
        y1 = y2 = strip_offset*sin(theta);
        z1 = -strip_width/lz;
        z2 = strip_width/lz;
        if(doMatrix){
          Cartesian cs = Cartesian(x1,y1,z1);
          cs *= mTransVector;
          xs1 = cs.get_x();
          ys1 = cs.get_y();
          zs1 = cs.get_z();
          Cartesian c = Cartesian(x1,y1,z1);
          c *= mTransVector;
          x1 = c.get_x()/lx/lx;
          y1 = c.get_y()/ly/ly;
          z1 = c.get_z()/lz/lz;
          cs = Cartesian(x2,y2,z2);
          cs *= mTransVector;
          xs2 = cs.get_x();
          ys2 = cs.get_y();
          zs2 = cs.get_z();
          c = Cartesian(x2,y2,z2);
          c = mTransVector*c;
          x2 = c.get_x()/lx/lx;
          y2 = c.get_y()/ly/ly;
          z2 = c.get_z()/lz/lz;
        }else {
          Cartesian cs = Cartesian(x1*p->GetA(),y1*p->GetB(),z1*p->GetC());
          cs *= mTransVector;
          xs1 = cs.get_x();
          ys1 = cs.get_y();
          zs1 = cs.get_z();
          Cartesian c = Cartesian(x1/p->GetA(),y1/p->GetB(),z1/p->GetC());
          c *= mTransVector;
          x1 = c.get_x();
          y1 = c.get_y();
          z1 = c.get_z();
          cs = Cartesian(x2*p->GetA(),y2*p->GetB(),z2*p->GetC());
          cs *= mTransVector;
          xs2 = cs.get_x();
          ys2 = cs.get_y();
          zs2 = cs.get_z();
          c = Cartesian(x2/p->GetA(),y2/p->GetB(),z2/p->GetC());
          c *= mTransVector;
          x2 = c.get_x();
          y2 = c.get_y();
          z2 = c.get_z();
        }
        bufferVertices[axesVertexOffset+k*my_accu*18+jj*6]   = prims[k]->GetVertices()[0].get_x()+prims[k]->GetSize()*xs1;
        bufferVertices[axesVertexOffset+k*my_accu*18+jj*6+1] = prims[k]->GetVertices()[0].get_y()+prims[k]->GetSize()*ys1;
        bufferVertices[axesVertexOffset+k*my_accu*18+jj*6+2] = prims[k]->GetVertices()[0].get_z()+prims[k]->GetSize()*zs1;
        bufferVertices[axesVertexOffset+k*my_accu*18+jj*6+3] = prims[k]->GetVertices()[0].get_x()+prims[k]->GetSize()*xs2;
        bufferVertices[axesVertexOffset+k*my_accu*18+jj*6+4] = prims[k]->GetVertices()[0].get_y()+prims[k]->GetSize()*ys2;
        bufferVertices[axesVertexOffset+k*my_accu*18+jj*6+5] = prims[k]->GetVertices()[0].get_z()+prims[k]->GetSize()*zs2;
        bufferNormals[axesVertexOffset+k*my_accu*18+jj*6]   = x1;
        bufferNormals[axesVertexOffset+k*my_accu*18+jj*6+1] = y1;
        bufferNormals[axesVertexOffset+k*my_accu*18+jj*6+2] = z1;
        bufferNormals[axesVertexOffset+k*my_accu*18+jj*6+3] = x2;
        bufferNormals[axesVertexOffset+k*my_accu*18+jj*6+4] = y2;
        bufferNormals[axesVertexOffset+k*my_accu*18+jj*6+5] = z2;
        if(show_solid){
          bufferColors[axesColorOffset+k*my_accu*24+jj*8]   = 0;
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+1] = 0;
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+2] = 0;
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+3] = alpha;
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+4] = 0;
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+5] = 0;
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+6] = 0;
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+7] = alpha;
        } else {
          bufferColors[axesColorOffset+k*my_accu*24+jj*8]   = prims[k]->GetColour()[0];
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+1] = prims[k]->GetColour()[1];
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+2] = prims[k]->GetColour()[2];
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+3] = alpha;
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+4] = prims[k]->GetColour()[0];
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+5] = prims[k]->GetColour()[1];
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+6] = prims[k]->GetColour()[2];
          bufferColors[axesColorOffset+k*my_accu*24+jj*8+7] = alpha;
        }

        y1 = y2 = strip_offset*cos(theta);
        z1 = z2 = strip_offset*sin(theta);
        x1 = -strip_width/lx;
        x2 = strip_width/lx;
        if(doMatrix){
          Cartesian cs = Cartesian(x1,y1,z1);
          cs *= mTransVector;
          xs1 = cs.get_x();
          ys1 = cs.get_y();
          zs1 = cs.get_z();
          Cartesian c = Cartesian(x1,y1,z1);
          c *= mTransVector;
          x1 = c.get_x()/lx/lx;
          y1 = c.get_y()/ly/ly;
          z1 = c.get_z()/lz/lz;
          cs = Cartesian(x2,y2,z2);
          cs *= mTransVector;
          xs2 = cs.get_x();
          ys2 = cs.get_y();
          zs2 = cs.get_z();
          c = Cartesian(x2,y2,z2);
          c *= mTransVector;
          x2 = c.get_x()/lx/lx;
          y2 = c.get_y()/ly/ly;
          z2 = c.get_z()/lz/lz;
        }else {
          Cartesian cs = Cartesian(x1*p->GetA(),y1*p->GetB(),z1*p->GetC());
          cs *= mTransVector;
          xs1 = cs.get_x();
          ys1 = cs.get_y();
          zs1 = cs.get_z();
          Cartesian c = Cartesian(x1/p->GetA(),y1/p->GetB(),z1/p->GetC());
          c *= mTransVector;
          x1 = c.get_x();
          y1 = c.get_y();
          z1 = c.get_z();
          cs = Cartesian(x2*p->GetA(),y2*p->GetB(),z2*p->GetC());
          cs *= mTransVector;
          xs2 = cs.get_x();
          ys2 = cs.get_y();
          zs2 = cs.get_z();
          c = Cartesian(x2/p->GetA(),y2/p->GetB(),z2/p->GetC());
          c *= mTransVector;
          x2 = c.get_x();
          y2 = c.get_y();
          z2 = c.get_z();
        }
        bufferVertices[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6]   = prims[k]->GetVertices()[0].get_x()+prims[k]->GetSize()*xs1;
        bufferVertices[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+1] = prims[k]->GetVertices()[0].get_y()+prims[k]->GetSize()*ys1;
        bufferVertices[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+2] = prims[k]->GetVertices()[0].get_z()+prims[k]->GetSize()*zs1;
        bufferVertices[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+3] = prims[k]->GetVertices()[0].get_x()+prims[k]->GetSize()*xs2;
        bufferVertices[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+4] = prims[k]->GetVertices()[0].get_y()+prims[k]->GetSize()*ys2;
        bufferVertices[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+5] = prims[k]->GetVertices()[0].get_z()+prims[k]->GetSize()*zs2;
        bufferNormals[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6]   = x1;
        bufferNormals[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+1] = y1;
        bufferNormals[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+2] = z1;
        bufferNormals[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+3] = x2;
        bufferNormals[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+4] = y2;
        bufferNormals[6*my_accu+axesVertexOffset+k*my_accu*18+jj*6+5] = z2;
        if(show_solid){
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8]   = 0;
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+1] = 0;
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+2] = 0;
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+3] = alpha;
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+4] = 0;
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+5] = 0;
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+6] = 0;
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+7] = alpha;
        } else {
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8]   = prims[k]->GetColour()[0];
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+1] = prims[k]->GetColour()[1];
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+2] = prims[k]->GetColour()[2];
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+3] = alpha;
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+4] = prims[k]->GetColour()[0];
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+5] = prims[k]->GetColour()[1];
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+6] = prims[k]->GetColour()[2];
          bufferColors[8*my_accu+axesColorOffset+k*my_accu*24+jj*8+7] = alpha;
        }

        x1 = x2 = strip_offset*cos(theta);
        z1 = z2 = strip_offset*sin(theta);
        y1 = strip_width/ly;
        y2 = -strip_width/ly;
        if(doMatrix){
          Cartesian cs = Cartesian(x1,y1,z1);
          cs *= mTransVector;
          xs1 = cs.get_x();
          ys1 = cs.get_y();
          zs1 = cs.get_z();
          Cartesian c = Cartesian(x1,y1,z1);
          c *= mTransVector;
          x1 = c.get_x()/lx/lx;
          y1 = c.get_y()/ly/ly;
          z1 = c.get_z()/lz/lz;
          cs = Cartesian(x2,y2,z2);
          cs *= mTransVector;
          xs2 = cs.get_x();
          ys2 = cs.get_y();
          zs2 = cs.get_z();
          c = Cartesian(x2,y2,z2);
          c *= mTransVector;
          x2 = c.get_x()/lx/lx;
          y2 = c.get_y()/ly/ly;
          z2 = c.get_z()/lz/lz;
        }else {
          Cartesian cs = Cartesian(x1*p->GetA(),y1*p->GetB(),z1*p->GetC());
          cs *= mTransVector;
          xs1 = cs.get_x();
          ys1 = cs.get_y();
          zs1 = cs.get_z();
          Cartesian c = Cartesian(x1/p->GetA(),y1/p->GetB(),z1/p->GetC());
          c *= mTransVector;
          x1 = c.get_x();
          y1 = c.get_y();
          z1 = c.get_z();
          cs = Cartesian(x2*p->GetA(),y2*p->GetB(),z2*p->GetC());
          cs *= mTransVector;
          xs2 = cs.get_x();
          ys2 = cs.get_y();
          zs2 = cs.get_z();
          c = Cartesian(x2/p->GetA(),y2/p->GetB(),z2/p->GetC());
          c *= mTransVector;
          x2 = c.get_x();
          y2 = c.get_y();
          z2 = c.get_z();
        }
        bufferVertices[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6]   = prims[k]->GetVertices()[0].get_x()+prims[k]->GetSize()*xs1;
        bufferVertices[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+1] = prims[k]->GetVertices()[0].get_y()+prims[k]->GetSize()*ys1;
        bufferVertices[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+2] = prims[k]->GetVertices()[0].get_z()+prims[k]->GetSize()*zs1;
        bufferVertices[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+3] = prims[k]->GetVertices()[0].get_x()+prims[k]->GetSize()*xs2;
        bufferVertices[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+4] = prims[k]->GetVertices()[0].get_y()+prims[k]->GetSize()*ys2;
        bufferVertices[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+5] = prims[k]->GetVertices()[0].get_z()+prims[k]->GetSize()*zs2;
        bufferNormals[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6]   = x1;
        bufferNormals[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+1] = y1;
        bufferNormals[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+2] = z1;
        bufferNormals[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+3] = x2;
        bufferNormals[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+4] = y2;
        bufferNormals[12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+5] = z2;
        if(show_solid){
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8]   = 0;
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+1] = 0;
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+2] = 0;
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+3] = alpha;
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+4] = 0;
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+5] = 0;
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+6] = 0;
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+7] = alpha;
        } else {
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8]   = prims[k]->GetColour()[0];
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+1] = prims[k]->GetColour()[1];
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+2] = prims[k]->GetColour()[2];
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+3] = alpha;
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+4] = prims[k]->GetColour()[0];
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+5] = prims[k]->GetColour()[1];
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+6] = prims[k]->GetColour()[2];
          bufferColors[16*my_accu+axesColorOffset+k*my_accu*24+jj*8+7] = alpha;
        }
        //std::cout << 12*my_accu+axesVertexOffset+k*my_accu*18+jj*6+5 << " " << nRealVerts << "\n";
      }
      }
    }

    if(show_solid){
      int ii=0;
      for(int i=0;i<=180;i=i+360/my_accu,ii++){
        int indexOffset = k*(my_accu+1)*my_accu;
        double theta = (double)i/360.0 * PIBY2;
        int jj=0;
        for(int j=0;j<=360;j=j+360/my_accu,jj++){
          if(ii<my_accu/2){
            if(j==0){
              // Prepend a degenerate triangle
              bufferIndices[indexOffsetInternal++] = indexOffset + jj + (ii+1) * my_accu;
            }
            if(jj==my_accu){
              bufferIndices[indexOffsetInternal++] = indexOffset + (ii+1) * my_accu;
              bufferIndices[indexOffsetInternal++] = indexOffset + ii * my_accu;
              // Append a degenerate triangle
              bufferIndices[indexOffsetInternal++] = indexOffset + (ii) * my_accu;
            } else {
              bufferIndices[indexOffsetInternal++] = indexOffset + jj + (ii+1) * my_accu;
              bufferIndices[indexOffsetInternal++] = indexOffset + jj + ii * my_accu;
            }
          }
          if(jj<my_accu){
          double phi = (double)j/360.0 * PIBY2;
          float z = cos(theta);
          float x = cos(phi)*sin(theta);
          float y = sin(phi)*sin(theta);
          float xs;
          float ys;
          float zs;
          if(doMatrix){
            Cartesian cs = Cartesian(x,y,z);
            cs *= mTransVector;
            xs = cs.get_x();
            ys = cs.get_y();
            zs = cs.get_z();
            Cartesian c = Cartesian(x,y,z);
            c *= mTransVector;
            x = c.get_x()/lx/lx;
            y = c.get_y()/ly/ly;
            z = c.get_z()/lz/lz;
          }else {
            Cartesian cs = Cartesian(x*p->GetA(),y*p->GetB(),z*p->GetC());
            cs *= mTransVector;
            xs = cs.get_x();
            ys = cs.get_y();
            zs = cs.get_z();
            Cartesian c = Cartesian(x/p->GetA(),y/p->GetB(),z/p->GetC());
            c *= mTransVector;
            x = c.get_x();
            y = c.get_y();
            z = c.get_z();
          }
          bufferVertices[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3  ] = prims[k]->GetVertices()[0].get_x()+prims[k]->GetSize()*xs;
          bufferVertices[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+1] = prims[k]->GetVertices()[0].get_y()+prims[k]->GetSize()*ys;
          bufferVertices[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+2] = prims[k]->GetVertices()[0].get_z()+prims[k]->GetSize()*zs;
          bufferNormals[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3  ] = x;
          bufferNormals[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+1] = y;
          bufferNormals[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+2] = z;
          //std::cout << k << " " << ii << " " << jj << "\n"; std::cout.flush();
          maxInd = k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+2;
          bufferColors[k*(my_accu+1)*my_accu*4+ii*(my_accu)*4+jj*4  ] = prims[k]->GetColour()[0];
          bufferColors[k*(my_accu+1)*my_accu*4+ii*(my_accu)*4+jj*4+1] = prims[k]->GetColour()[1];
          bufferColors[k*(my_accu+1)*my_accu*4+ii*(my_accu)*4+jj*4+2] = prims[k]->GetColour()[2];
          bufferColors[k*(my_accu+1)*my_accu*4+ii*(my_accu)*4+jj*4+3] = alpha;
          }
        }
      }
    }
  }

  //std::cout << "indexOffsetInternal: " << indexOffsetInternal << "\n";

  polygonType = GL_TRIANGLE_STRIP;

  return a;
}

std::vector<Primitive*> SphereCollection::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  std::vector<Primitive*> a;

  //std::cout << "SphereCollection::GetSimplePrimitivesArrays\n"; std::cout.flush();

  unsigned nSpheres = prims.size();

  if(nSpheres==0) return a;

  // /*
  //nSpheres = 30;

  //std::cout << "nSpheres: " << nSpheres << "\n"; std::cout.flush();

  static int spheres_accus[] = {4,9,18};


  if(sphere_accu>0){
  int my_accu = 8;
  if(sphere_accu>=0&&sphere_accu<3) my_accu = spheres_accus[sphere_accu];

  //std::cout << "accu:" << sphere_accu << "\n";
  //std::cout << "my_accu:" << my_accu << "\n";
  //std::cout << "alpha:" << GetAlpha() << "\n";

  bufferVertices = new float[3*(my_accu+1)*my_accu*nSpheres];
  bufferNormals = new float[3*(my_accu+1)*my_accu*nSpheres];
  bufferColors = new float[4*(my_accu+1)*my_accu*nSpheres];
  nRealVerts = 3*(my_accu+1)*my_accu*nSpheres;
  bufferIndices = new GLuint[2*(my_accu+1)*my_accu*nSpheres + 2*my_accu*nSpheres];
  int indexOffsetInternal = 0;
  nVerts = 2*(my_accu+1)*my_accu*nSpheres + 2*my_accu*nSpheres;
  int maxInd = 0;
  //std::cout << "nSpheres: " << nSpheres << "\n"; std::cout.flush();
  for(unsigned k=0;k<nSpheres&&k<prims.size();k++){
    int ii=0;
    for(int i=0;i<=180;i=i+180/my_accu,ii++){
      int indexOffset = k*(my_accu+1)*my_accu;
      double theta = (double)i/360.0 * PIBY2;
      int jj=0;
      for(int j=0;j<=360;j=j+360/my_accu,jj++){
        if(ii<my_accu){
          if(j==0){
            // Prepend a degenerate triangle
            bufferIndices[indexOffsetInternal++] = indexOffset + jj + (ii+1) * my_accu;
          }
          if(jj==my_accu){
          bufferIndices[indexOffsetInternal++] = indexOffset + (ii+1) * my_accu;
          bufferIndices[indexOffsetInternal++] = indexOffset + ii * my_accu;
            // Append a degenerate triangle
            bufferIndices[indexOffsetInternal++] = indexOffset + (ii) * my_accu;
          } else {
            bufferIndices[indexOffsetInternal++] = indexOffset + jj + (ii+1) * my_accu;
            bufferIndices[indexOffsetInternal++] = indexOffset + jj + ii * my_accu;
          }
        }
        if(jj<my_accu){
        double phi = (double)j/360.0 * PIBY2;
        float z = cos(theta);
        float x = cos(phi)*sin(theta);
        float y = sin(phi)*sin(theta);
        bufferVertices[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3  ] = prims[k]->GetVertices()[0].get_x()+prims[k]->GetSize()*x;
        bufferVertices[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+1] = prims[k]->GetVertices()[0].get_y()+prims[k]->GetSize()*y;
        bufferVertices[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+2] = prims[k]->GetVertices()[0].get_z()+prims[k]->GetSize()*z;
        bufferNormals[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3  ] = x;
        bufferNormals[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+1] = y;
        bufferNormals[k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+2] = z;
        //std::cout << k << " " << ii << " " << jj << "\n"; std::cout.flush();
        maxInd = k*(my_accu+1)*my_accu*3+ii*(my_accu)*3+jj*3+2;
        bufferColors[k*(my_accu+1)*my_accu*4+ii*(my_accu)*4+jj*4  ] = prims[k]->GetColour()[0];
        bufferColors[k*(my_accu+1)*my_accu*4+ii*(my_accu)*4+jj*4+1] = prims[k]->GetColour()[1];
        bufferColors[k*(my_accu+1)*my_accu*4+ii*(my_accu)*4+jj*4+2] = prims[k]->GetColour()[2];
        bufferColors[k*(my_accu+1)*my_accu*4+ii*(my_accu)*4+jj*4+3] = alpha;
        }
      }
    }
  }
  //std::cout << nVerts << " " << indexOffsetInternal << "\n";
  //std::cout << nRealVerts << " " << maxInd << "\n";
  polygonType = GL_TRIANGLE_STRIP;
  // */

  } else {
  
  bufferVertices = new float[36*prims.size()];
  bufferNormals = new float[36*prims.size()];
  bufferColors = new float[48*prims.size()];
  bufferIndices = new GLuint[60*prims.size()];
  nRealVerts = 36*prims.size();
  nVerts = 60*prims.size();
  float offsets[][3] = {{-1,-1,0},{1,-1,0},{0,1.4,0}};

  Cartesian n = Cartesian(0,0,1);

  for(unsigned i=0;i<prims.size();i++){
    for(int j=0;j<20;j++){
      bufferIndices[i*60+j*3  ] = i*12+tindices[j][0];
      bufferIndices[i*60+j*3+2] = i*12+tindices[j][1];
      bufferIndices[i*60+j*3+1] = i*12+tindices[j][2];
    }
  }
  for(unsigned i=0;i<prims.size();i++){
    for(int j=0;j<12;j++){
      bufferVertices[i*36+j*3  ] = prims[i]->GetVertices()[0].get_x()+prims[i]->GetSize()*icosa[j][0];
      bufferVertices[i*36+j*3+1] = prims[i]->GetVertices()[0].get_y()+prims[i]->GetSize()*icosa[j][1];
      bufferVertices[i*36+j*3+2] = prims[i]->GetVertices()[0].get_z()+prims[i]->GetSize()*icosa[j][2];
      bufferNormals[i*36+j*3  ] = icosa[j][0];
      bufferNormals[i*36+j*3+1] = icosa[j][1];
      bufferNormals[i*36+j*3+2] = icosa[j][2];
      bufferColors[i*48+j*4  ] = prims[i]->GetColour()[0];
      bufferColors[i*48+j*4+1] = prims[i]->GetColour()[1];
      bufferColors[i*48+j*4+2] = prims[i]->GetColour()[2];
      bufferColors[i*48+j*4+3] = alpha;
    }
  }

  polygonType = GL_TRIANGLES;
  }

  return a;
};

std::vector<Primitive*> ArrowHeadRibbon::GetSimplePrimitivesArrays(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) {
  std::vector<Primitive*> a;
  if(vertices.size()<2) return a;
  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;
  const double *colp = GetColour();
  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);

  std::vector<Cartesian> vertices_tmp = vertices;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;
  if(RenderQuality::GetSmoothRibbons()) {
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    colour_vector_tmp.resize(colour_vector.size()*4);
    for(unsigned icol=0;icol<colour_vector.size();icol++){
      colour_vector_tmp[icol*4] = colour_vector[icol];
      colour_vector_tmp[icol*4+1] = colour_vector[icol];
      colour_vector_tmp[icol*4+2] = colour_vector[icol];
      colour_vector_tmp[icol*4+3] = colour_vector[icol];
    }
  }

  if(insert_at_vertex>0){
    vertices_tmp.insert(vertices_tmp.begin()+insert_at_vertex,0.9*vertices_tmp[insert_at_vertex-1] + 0.1*vertices_tmp[insert_at_vertex]);
    if(colour_vector_tmp.size()>0) colour_vector_tmp.insert(colour_vector_tmp.begin()+insert_at_vertex,colour_vector_tmp[insert_at_vertex-1]);
  }

  int my_accu = accu;
  if(accu<4) my_accu = 8;

  if(RenderQuality::GetRenderQuality()&&0)
    my_accu = 180;
  if(style==0){
    get_elliptical_ribbon_array(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,false);
  }else if(style==1){
    get_flat_ribbon_array(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,false);
  }else if(style==2){
    get_flat_rounded_ribbon_array(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,false);
  }

  return a;
}

std::vector<Primitive*> Ribbon::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;
  const double *colp = GetColour();
  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);

  std::vector<Cartesian> vertices_tmp = vertices;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;
  if(RenderQuality::GetSmoothRibbons()) {
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    colour_vector_tmp.resize(colour_vector.size()*4);
    for(unsigned icol=0;icol<colour_vector.size();icol++){
      colour_vector_tmp[icol*4] = colour_vector[icol];
      colour_vector_tmp[icol*4+1] = colour_vector[icol];
      colour_vector_tmp[icol*4+2] = colour_vector[icol];
      colour_vector_tmp[icol*4+3] = colour_vector[icol];
    }
  }

  int my_accu = accu;
  if(accu<4) my_accu = 8;

  if(RenderQuality::GetRenderQuality()&&0)
    my_accu = 180;
  if(style==0)
    get_elliptical_ribbon(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha);
  else if(style==1)
    get_flat_ribbon(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,two_colour,edges_grey);
  else if(style==2)
    get_flat_rounded_ribbon(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,two_colour);
  else if(style==3)
    get_fancy_ribbon(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,two_colour);

  for(int i=0;i<360;i=i+360/my_accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/my_accu)/360.0 * PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    double x2 = cos(theta2);
    double y2 = sin(theta2);
    double col[4];
    if(AB[0][0].length()>0.21||colour_vector_tmp[0].get_a()>0.98){
    std::vector<Cartesian> v;
    std::vector<Cartesian> ns;
    v.push_back(vertices_tmp[0]);
    v.push_back(vertices_tmp[0]+x*AB[0][0]+y*AB[1][0]);
    v.push_back(vertices_tmp[0]+x2*AB[0][0]+y2*AB[1][0]);
    Cartesian n1 = v[0] - v[1];
    Cartesian n2 = v[0] - v[2];
    n1.normalize();
    n2.normalize();
    Cartesian n = Cartesian::CrossProduct(n1,n2);
    ns.push_back(n);
    ns.push_back(n);
    ns.push_back(n);
    col[0]=colour_vector_tmp[0].get_x();
    col[1]=colour_vector_tmp[0].get_y();
    col[2]=colour_vector_tmp[0].get_z();
    col[3]=colour_vector_tmp[0].get_a();
    TriangleElement *tri = new TriangleElement(v,col,0.333*(v[0]+v[1]+v[2]),1,0);
    tri->SetNormals(ns);
    a.push_back(tri);
    }
    if(AB[0][vertices_tmp.size()-1].length()>0.21||colour_vector_tmp[vertices_tmp.size()-1].get_a()>0.98){
    std::vector<Cartesian> v;
    std::vector<Cartesian> ns;
    v.push_back(vertices_tmp[vertices_tmp.size()-1]);
    v.push_back(vertices_tmp[vertices_tmp.size()-1]+x*AB[0][vertices_tmp.size()-1]+y*AB[1][vertices_tmp.size()-1]);
    v.push_back(vertices_tmp[vertices_tmp.size()-1]+x2*AB[0][vertices_tmp.size()-1]+y2*AB[1][vertices_tmp.size()-1]);
    Cartesian n1 = v[0] - v[1];
    Cartesian n2 = v[0] - v[2];
    n1.normalize();
    n2.normalize();
    Cartesian n = Cartesian::CrossProduct(n1,n2);
    ns.push_back(n);
    ns.push_back(n);
    ns.push_back(n);
    col[0]=colour_vector_tmp[vertices_tmp.size()-1].get_x();
    col[1]=colour_vector_tmp[vertices_tmp.size()-1].get_y();
    col[2]=colour_vector_tmp[vertices_tmp.size()-1].get_z();
    col[3]=colour_vector_tmp[vertices_tmp.size()-1].get_a();
    TriangleElement *tri = new TriangleElement(v,col,0.333*(v[0]+v[1]+v[2]),1,0);
    tri->SetNormals(ns);
    a.push_back(tri);
    }
  }

  return a;
}

std::vector<Primitive*> ArrowHeadRibbon::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  int multicolour = 0;
  if(colour_vector.size()>0)
    multicolour = 1;
  const double *colp = GetColour();
  int insert_at_vertex;
  std::vector<std::vector<Cartesian> > AB = GetAB(insert_at_vertex);

  std::vector<Cartesian> vertices_tmp = vertices;
  std::vector<Cartesian> colour_vector_tmp = colour_vector;
  if(RenderQuality::GetSmoothRibbons()) {
    vertices_tmp = SplineCurve(vertices,(vertices.size()-1)*4,3,1);
    colour_vector_tmp.resize(colour_vector.size()*4);
    for(unsigned icol=0;icol<colour_vector.size();icol++){
      colour_vector_tmp[icol*4] = colour_vector[icol];
      colour_vector_tmp[icol*4+1] = colour_vector[icol];
      colour_vector_tmp[icol*4+2] = colour_vector[icol];
      colour_vector_tmp[icol*4+3] = colour_vector[icol];
    }
  }

  if(insert_at_vertex>0){
    vertices_tmp.insert(vertices_tmp.begin()+insert_at_vertex,0.9*vertices_tmp[insert_at_vertex-1] + 0.1*vertices_tmp[insert_at_vertex]);
    if(colour_vector_tmp.size()>0) colour_vector_tmp.insert(colour_vector_tmp.begin()+insert_at_vertex,colour_vector_tmp[insert_at_vertex-1]);
  }

  int my_accu = accu;
  if(accu<4) my_accu = 8;

  if(RenderQuality::GetRenderQuality()&&0)
    my_accu = 180;
  if(style==0)
    get_elliptical_ribbon(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha);
  else if(style==1)
    get_flat_ribbon(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,false,edges_grey);
  else if(style==2)
    get_flat_rounded_ribbon(a,vertices_tmp,AB,my_accu,textured,multicolour,colour_vector_tmp,colp,alpha,false);

  return a;
}


void normalize(float v[3],float radius);
void SubDivideSphere(float *v1, float *v2, float *v3, int depth, float radius, std::vector<Primitive*> &a, const double *colour, double alpha, int textured, Cartesian origin){ 
   GLfloat v12[3], v23[3], v31[3];    
   GLint i;

   normalize(v1,radius);    
   normalize(v2,radius); 
   normalize(v3,radius); 

   if (depth <= 0) {
     std::vector<Cartesian> carts;
     std::vector<Cartesian> normals;
     normals.push_back(Cartesian(v1[0],v1[1],v1[2]));
     normals.push_back(Cartesian(v3[0],v3[1],v3[2]));
     normals.push_back(Cartesian(v2[0],v2[1],v2[2]));
     carts.push_back(Cartesian(v1[0],v1[1],v1[2])+origin);
     carts.push_back(Cartesian(v3[0],v3[1],v3[2])+origin);
     carts.push_back(Cartesian(v2[0],v2[1],v2[2])+origin);
     TriangleElement *q = new TriangleElement(carts,colour,(carts[0]+carts[1]+carts[2])/3.,alpha,textured);
     q->SetNormals(normals);
     a.push_back(q);
     return;
   }

   for (i = 0; i < 3; i++) { 
      v12[i] = v1[i]+v2[i]; 
      v23[i] = v2[i]+v3[i];     
      v31[i] = v3[i]+v1[i];    
   } 
   normalize(v12,radius);    
   normalize(v23,radius); 
   normalize(v31,radius); 

   SubDivideSphere(v1, v12, v31, depth-1,radius,a,colour,alpha,textured,origin);    
   SubDivideSphere(v2, v23, v12, depth-1,radius,a,colour,alpha,textured,origin);    
   SubDivideSphere(v3, v31, v23, depth-1,radius,a,colour,alpha,textured,origin);    
   SubDivideSphere(v12, v23, v31, depth-1,radius,a,colour,alpha,textured,origin); 

}

std::vector<Primitive*> SphereElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  for (int i = 0; i < 20; i++) {    
    SubDivideSphere(&icosa[tindices[i][0]][0],  
                    &icosa[tindices[i][1]][0],  
                    &icosa[tindices[i][2]][0],accu,size,a,colour,alpha,textured,vertices[0]); 
  }

  bool clipped = false;
  std::vector<Primitive*>::iterator this_simple_iter = a.begin();
  while(this_simple_iter!=a.end()){
   std::vector<Cartesian> carts = (*this_simple_iter)->GetVertices();
   if(!clip_vol.PointInVolume(objrotmatrix*(carts[0]+objorigin)))
     clipped = true;
   if(!clip_vol.PointInVolume(objrotmatrix*(carts[1]+objorigin)))
     clipped = true;
   if(!clip_vol.PointInVolume(objrotmatrix*(carts[2]+objorigin)))
     clipped = true;
   this_simple_iter++;
  }
  if(clipped) a.clear();
  return a;
}

std::vector<Primitive*> SpheroidElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {

  std::vector<Primitive*> simp_prims2;
  SpheroidElement *se = new SpheroidElement(vertices[0],colour,origin,a,b,c,xaxis,yaxis,zaxis,alpha,show_axes,show_solid,textured);
  simp_prims2.push_back(se);
  return simp_prims2;

  std::vector<Primitive*> simp_prims;
  Cartesian cartx = Cartesian(1,0,0);
  Cartesian carty = Cartesian(0,1,0);
  Cartesian yprim;
  Quat q1,q2;
  matrix m(4,kdelta);

  std::cout << xaxis << "\n";
  std::cout << yaxis << "\n";
  std::cout << zaxis << "\n";
  std::cout << a << " " << b << " " << c << "\n";
  double angle = acos(cartx.DotProduct(cartx,xaxis))*180.0/M_PI;
  std::cout << "Angle 1: " << angle << "\n";
  if(fabs(angle)>0.00001){
    Cartesian rotax = cartx.CrossProduct(cartx,xaxis);
    std::cout << "rotax: " << rotax;
    q1 = Quat(rotax,1,-angle);
  }
  m = q1.getInvMatrix();

  yprim = m * carty;
  angle = acos(yprim.DotProduct(yprim,yaxis))*180.0/M_PI;
  std::cout << "Angle 2: " << angle << "\n";

  if(fabs(angle)>0.00001){
    q2 = Quat(cartx,1,-angle);
    q1.postMult(q2);
  }
  m = q1.getMatrix();

  for (int i = 0; i < 20; i++) {    
    SubDivideSphere(&icosa[tindices[i][0]][0],  
                    &icosa[tindices[i][1]][0],  
                    &icosa[tindices[i][2]][0],accu,size,simp_prims,colour,alpha,textured,vertices[0]); 
  }

  for(unsigned ii=0;ii<simp_prims.size();ii++){
    std::vector<Cartesian> carts = simp_prims[ii]->GetVertices();
    carts[0].Scale(a,b,c);
    carts[0] = m * carts[0];
    carts[1].Scale(a,b,c);
    carts[1] = m * carts[1];
    carts[2].Scale(a,b,c);
    carts[2] = m * carts[2];
    simp_prims[ii]->SetVertices(carts);
  }

  return simp_prims;
}

std::vector<Primitive*> CylinderElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  std::vector<Cartesian> colour_vector;
  colour_vector.push_back(Cartesian(colour[0],colour[1],colour[2]));
  colour_vector.push_back(Cartesian(colour[0],colour[1],colour[2]));
  PolyCylinder polycyl(vertices, colour, origin, colour_vector, size, alpha, accu, textured);
  return polycyl.GetSimplePrimitives(clip_vol,objrotmatrix,objorigin,start,end);
}

std::vector<Primitive*> ConeElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  std::vector<std::vector<Cartesian> > AB = GetPolyCylinderAB(vertices,size);

  for(int i=0;i<360;i=i+360/accu){
    double theta = (double)i/360.0 * PIBY2;
    double theta2 = (double)(i-360/accu)/360.0 * PIBY2;
    std::vector<Cartesian> v;
    for(unsigned int j=0;j<vertices.size();j++){
      double x = cos(theta);
      double y = sin(theta);
      Cartesian p;
      if(j==0)
        p = x*AB[0][j] + y*AB[1][j] + vertices[j];
      else
        p = vertices[j];
      v.push_back(p);
      x = cos(theta2);
      y = sin(theta2);
      if(j==0)
        p = x*AB[0][j] + y*AB[1][j] + vertices[j];
      else
        p = vertices[j];
      v.push_back(p);
      if(v.size()==4){
        QuadElement *q = new QuadElement(v,colour,v[0],alpha,textured);
        a.push_back(q);
        std::vector<Cartesian> vnew;
        vnew.push_back(v[2]);
        vnew.push_back(v[3]);
        v = vnew;
      }
    }
  }

  return a;
}

void DashCylinder::draw(const double *override_colour, int selective_override){
  glEnable(GL_LIGHTING);
  DashCylinderElement::draw(override_colour,selective_override);
}

void DashCylinderElement::draw(const double *override_colour, int selective_override){
  
  double length = (vertices[1]-vertices[0]).length();
  int nsegments = (int)(length/dash_length);
  // Force a gap at ends dashed line - probably correct
  // for interatomic lines
  if ((nsegments%2) == 0) nsegments++;
  std::vector<Cartesian> carts(2);

  for(int i=0;i<nsegments;i++){
    if((i%2)==dash_end){
      double frac = (double)i/(double)nsegments;
      double frac2 =(double)(i+1)/(double)nsegments;
      carts[0] = frac*vertices[0] + (1-frac)*vertices[1];
      carts[1] = frac2*vertices[0] + (1-frac2)*vertices[1];
      draw_capped_cylinder(carts,size,accu,textured);
    }
  }

}

void DashCylinderElement::DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v){
  if(!v.PointInVolume(objrotmatrix*(vertices[0]+objorigin))&&!v.PointInVolume(objrotmatrix*(vertices[1]+objorigin))) {
    return;
  }

  double length = (vertices[1]-vertices[0]).length();
  int nsegments = (int)(length/dash_length);
  if ((nsegments%2) == 0) nsegments++;
  std::vector<Cartesian> carts(2);

  for(int i=0;i<nsegments;i++){
    if((i%2)==1){
      double frac = (double)i/(double)nsegments;
      double frac2 =(double)(i+1)/(double)nsegments;
      carts[0] = frac*vertices[0] + (1-frac)*vertices[1];
      carts[1] = frac2*vertices[0] + (1-frac2)*vertices[1];

      fp << "cylinder {\n";
      Cartesian p = quat.getMatrix()*(objrotmatrix*(carts[0]+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      p = quat.getMatrix()*(objrotmatrix*(carts[1]+objorigin)+Cartesian(ox,oy,oz));
      fp << "  < " << p.get_x() << ", " << p.get_y() << ", " << -p.get_z() << ">,\n";
      fp << "  " << GetSize() << "\n";
      fp << "  texture {\n";
      OutputPovrayPigment(fp); 
      fp << "    finish {diffuse 1.0 specular 1.0}\n";
      fp << "  }\n";
      fp << "}\n";

    }
  }

}

void DashCylinderElement::DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v){
   Primitive::DrawPostScript(fp,quat,radius,ox,oy,oz,objrotmatrix,objorigin,xoff,yoff,xscale,yscale,xscaleps,v);
}

void DashCylinderElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

std::vector<Primitive*> DashCylinderElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
    double length = (vertices[1]-vertices[0]).length();
  int nsegments = (int)(length/dash_length);
  // Force a gap at ends dashed line - probably correct
  // for interatomic lines
  if ((nsegments%2) == 0) nsegments++;
  std::vector<Cartesian> carts(2);

  for(int i=0;i<nsegments;i++){
    if((i%2)==dash_end){
      double frac = (double)i/(double)nsegments;
      double frac2 =(double)(i+1)/(double)nsegments;
      carts[0] = frac*vertices[0] + (1-frac)*vertices[1];
      carts[1] = frac2*vertices[0] + (1-frac2)*vertices[1];
      Cylinder *cylinder = new Cylinder(carts,colour,carts[0],size,alpha,accu,textured);
      std::vector<Primitive*> b = cylinder->GetSimplePrimitives(clip_vol,objrotmatrix,objorigin);
      a.insert(a.end(),b.begin(),b.end());
    }
  }
  return a;
}
