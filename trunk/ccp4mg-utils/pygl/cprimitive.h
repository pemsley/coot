/*
     pygl/cprimitive.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC

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


#ifndef _CCP4MG_PRIM_
#define _CCP4MG_PRIM_
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "quat.h"
#include "cartesian.h"
#include "volume.h"
#include "ppmutil.h"
#include "font_info.h"

/* My Ubuntu Intripid GCC 4.3.2 goes via the the __APPLE_CC__ path. So
   I am commenting this out. 20090106 - PE */
/* #ifdef __APPLE_CC__ */
/* #include <OpenGL/gl.h> */
/* #include <OpenGL/glu.h> */
/* #else */
#include <GL/gl.h>
#include <GL/glu.h>
/* #endif */

void draw_flat_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
void draw_flat_rounded_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
void draw_fancy_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
void draw_elliptical_ribbon(const std::vector<Cartesian> &vertices, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int quality, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
void draw_billboard(double v11, double v12, double v21, double v22);
void draw_quadstrip(double *vertices, int nvertices, int textured);
void draw_quadstrip_outline(double *vertices, int nvertices);
void draw_quadstrip_textured(double *vertices, int nvertices, int textured);
void set_line_colour(double r, double g, double b, double a);
void draw_line(double v11, double v12, double v13, double v21, double v22, double v23, double size);
void draw_sphere(double v1, double v2, double v3, int accu, double size);

void draw_cylinder(const std::vector<Cartesian> &vertices, double size,int accu, int textured, bool force_dl=false, bool capped=false);
void draw_capped_cylinder(const std::vector<Cartesian> &vertices, double size,int accu, int textured, bool force_dl=false);

void ForceLoadTextures();
void UnForceLoadTextures();

class RenderQuality {
  static int render_quality;
 public:
  static void SetRenderQuality(int render_quality_in){render_quality = render_quality_in; }
  static int GetRenderQuality(){ return render_quality ; }
};

class Lighting{
 public:
  Lighting();
  ~Lighting(){};
  GLfloat ambient[3];
  GLfloat diffuse[3];
  GLfloat specular[3];
  GLfloat emission[3];
  GLfloat shininess;
  void copy(const Lighting &light_in);
};

class Primitive{
 protected:
  double colour[3];
  std::vector<Cartesian> colours;
  std::vector<Cartesian> vertices;
  std::vector<Cartesian> normals;
  int textured;
  int multitexture;
  double alpha;
  double size;
  Quat quat;
  Cartesian origin;
  int transparent;
  int colour_override;
  void OutputPovrayPigment(std::ofstream &fp);
 public:
  Primitive();
  virtual ~Primitive();
  Cartesian get_origin(void) const { return origin;};
  void set_origin(Cartesian origin_in);
  void set_origin(double *origin_in);
  void SetColourOverride(int colour_override_in){colour_override = colour_override_in;};
  int GetColourOverride(void) const { return colour_override;};
  virtual void draw(double *override_colour=0, int selective_override=0) = 0 ;
  virtual void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  virtual void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  virtual void set_draw_colour(GLfloat *colp=0) = 0;
  virtual void SetAlpha(double alpha_in){ alpha = alpha_in; } ;
  double GetAlpha() const { return alpha; } ;
  void set_draw_colour_poly(void) ;
  void set_draw_colour_line(void) ;
  void set_draw_colour_poly_override(GLfloat *colp) ;
  void set_draw_colour_line_override(GLfloat *colp) ;
  void initialize_lighting(void);
  void apply_textures(void);
  void unbind_textures(void);
  int get_transparent(void) const {return transparent;};
  virtual void set_transparent(int trans_in){transparent=trans_in;};
  void SetColours(const std::vector<Cartesian> &colours_in){colours = colours_in;};
  std::vector<Cartesian> GetColours(void) const {return colours;};
  void SetNormals(const std::vector<Cartesian> &normals_in){normals = normals_in;};
  std::vector<Cartesian> GetNormals(void) const {return normals;};
  void SetVertices(const std::vector<Cartesian> &vertices_in){vertices = vertices_in;};
  std::vector<Cartesian> GetVertices(void) const {return vertices;};
  const double* GetColour(void) const {return colour;};
  double GetSize(void) const {return size;};
  void SetSize(double size_in){size=size_in;};
  void SetOrigin(const Cartesian &o_in) { origin = o_in;};
  void SetColour(const double *c_in) { colour[0] = c_in[0];colour[1] = c_in[1];colour[2] = c_in[2];};
  void SetColour(const std::vector<double> &c_in) { colour[0] = c_in[0];colour[1] = c_in[1];colour[2] = c_in[2];};
  virtual bool isLine() const {return false;};
  virtual bool isDots() const {return false;};
  virtual std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class Point : public Primitive{
 public:
  Point();
  Point(const Cartesian &vertex, double *colour, const Cartesian &origin, double size_in=1.0, double alpha_in=1.0);
  ~Point();
  void set_draw_colour(GLfloat *col=0); 
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  bool isDots() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class PointElement : public Point {
 public:
  PointElement();
  ~PointElement();
  PointElement(const Cartesian &vertex, double *colour, const Cartesian &origin, double size_in=1.0, double alpha_in=1.0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return true;};
};

class Line : public Primitive{
 public:
  Line();
  Line(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0);
  ~Line();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class LineElement : public Line {
 public:
  LineElement();
  ~LineElement();
  LineElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return true;};
};

class ArrowElement : public LineElement {
 protected:
  int arrow_head;
 public:
  ArrowElement();
  ArrowElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0,int arrow_head=0);
  ~ArrowElement();
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class Arrow : public ArrowElement {
 protected:
  int arrow_head;
 public:
  Arrow();
  Arrow(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0,int arrow_head=0);
  ~Arrow();
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class DashArrowElement : public LineElement  {
 protected:
  int arrow_head;
 public:
  DashArrowElement();
  DashArrowElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0,int arrow_head=0);
  ~DashArrowElement();
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  void SetDashLength(double length) { dash_length = length; };
  double GetDashLength() const { return dash_length; };
 private:
  double dash_length;

};

class DashArrow : public DashArrowElement {
 protected:
  int arrow_head;
 public:
  DashArrow();
  DashArrow(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0,int arrow_head=0);
  ~DashArrow();
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  void SetDashLength(double length) { dash_length = length; };
  double GetDashLength() const { return dash_length; };
 private:
  double dash_length;

};

class LineStrip : public Primitive{
 public:
  LineStrip();
  LineStrip(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0);
  ~LineStrip();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class DashLine : public Primitive{
 public:
  DashLine();
  DashLine(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0);
  ~DashLine();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  void SetDashLength(double length) { dash_length = length; };
  double GetDashLength() const { return dash_length; };
 private:
  double dash_length;

};

class Circle : public Primitive{
 public:
  Circle();
  Circle(const Cartesian &vertex, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int textured_in=0);
  ~Circle();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class TriangleElement: public Primitive{
 public:
  TriangleElement();
  TriangleElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  ~TriangleElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class PolygonElement: public Primitive{
 public:
  PolygonElement();
  PolygonElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  ~PolygonElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class QuadElement: public Primitive{
 public:
  QuadElement();
  QuadElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  ~QuadElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class SphereElement: public Primitive{
  int accu;
 public:
  SphereElement();
  SphereElement(const Cartesian &vertex, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int accu_in=0, int textured_in=0);
  ~SphereElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class SpheroidElement: public Primitive{
  int accu;
  double a;
  double b;
  double c;
  Cartesian xaxis;
  Cartesian yaxis;
  Cartesian zaxis;
  int show_axes;
  int show_solid;
 public:
  SpheroidElement();
  SpheroidElement(const Cartesian &vertex, double *colour_in, const Cartesian &origin_in, double a_in=1.0, double b_in=1.0, double c_in=1.0, const Cartesian& xaxis_in=Cartesian(1,0,0), const Cartesian& yaxis_in=Cartesian(0,1,0), const Cartesian& zaxis_in=Cartesian(0,0,1), double alpha_in=1.0, int accu_in=0, int show_axes=0, int show_solid=1, int textured_in=0);
  ~SpheroidElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void ShowSolid();
  void HideSolid();
  void ShowAxes();
  void HideAxes();
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  double GetA() const { return a;}
  double GetB() const { return b;}
  double GetC() const { return c;}
  Cartesian GetXAxis() const { return xaxis;}
  Cartesian GetYAxis() const { return yaxis;}
  Cartesian GetZAxis() const { return zaxis;}
  void SetA(const double a_in) { a=a_in;}
  void SetB(const double b_in) { b=b_in;}
  void SetC(const double c_in) { c=c_in;}
  void SetXAxis(const Cartesian &xaxis_in) { xaxis=xaxis_in;}
  void SetYAxis(const Cartesian &yaxis_in) { yaxis=yaxis_in;}
  void SetZAxis(const Cartesian &zaxis_in) { zaxis=zaxis_in;}
};

class CylinderElement: public Primitive{
 protected:
  int accu;
 public:
  CylinderElement();
  CylinderElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int accu_in=8, int textured_in=0);
  ~CylinderElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class DashCylinderElement: public CylinderElement {
 public:
  DashCylinderElement();
  DashCylinderElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int accu_in=8, int textured_in=0);
  ~DashCylinderElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  void SetDashLength(double length) { dash_length = length; };
  double GetDashLength() const { return dash_length; };
  void SetDashEnd(int end) { dash_end = end; };
  int GetDashEnd() { return dash_end; };
 private:
  double dash_length;
  int dash_end;

};

class ConeElement: public CylinderElement{
 public:
  ConeElement();
  ConeElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int accu_in=8, int textured_in=0);
  ~ConeElement();
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class QuadStripElement : public Primitive{
 public:
  QuadStripElement();
  QuadStripElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  ~QuadStripElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class QuadStrip: public QuadStripElement {
 public:
  QuadStrip();
  ~QuadStrip(){};
  QuadStrip(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class TriangleStripElement : public Primitive{
 public:
  TriangleStripElement();
  TriangleStripElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  ~TriangleStripElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class TriangleStrip: public TriangleStripElement {
 public:
  TriangleStrip();
  ~TriangleStrip(){};
  TriangleStrip(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};


class TriangleFanElement : public Primitive{
 public:
  TriangleFanElement();
  TriangleFanElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  ~TriangleFanElement();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class TriangleFan: public TriangleFanElement {
 public:
  TriangleFan();
  ~TriangleFan(){};
  TriangleFan(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class PolyCylinder : public Primitive{
  int accu;
  std::vector<Cartesian>colour_vector;
 public:
  PolyCylinder();
  PolyCylinder(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double size_in=0.2, double alpha_in=1.0, int accu_in=4, int textured_in=0);
  ~PolyCylinder();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class Ribbon : public Primitive{
 protected:
  int accu;
  int scale_steps;
  std::vector<Cartesian> n1_vertices;
  std::vector<Cartesian> n2_vertices;
  double thickness;
  double minsize;
  double maxsize;
  std::vector<Cartesian>colour_vector;
 public:
  Ribbon();
  Ribbon(const std::vector<Cartesian> &vertices_in, const std::vector<Cartesian> &n1_vertices_in, const std::vector<Cartesian> &n2_vertices_in, double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double minsize_in=0.2, double maxsize_in=1.0, double thickness_in=0.2, double alpha_in=1.0, int accu_in=8, int scale_steps=8, int textured_in=0);
  ~Ribbon();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  virtual std::vector<std::vector<Cartesian> > GetAB();
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

class Worm : public Ribbon {
 public:
  Worm();
  Worm(const std::vector<Cartesian> &vertices_in, const std::vector<Cartesian> &n1_vertices_in, const std::vector<Cartesian> &n2_vertices_in, double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double minsize_in=0.2, double maxsize_in=1.0, double thickness_in=0.2, double alpha_in=1.0, int accu_in=8, int scale_steps=8, int textured_in=0);
  ~Worm();
  void draw(double *override_colour=0, int selective_override=0);
};

class ArrowHeadRibbon : public Ribbon {
  double arrow_length;
  double arrow_width;
  int insert_at_vertex;
  std::vector<Cartesian> orig_vertices;
  std::vector<Cartesian> orig_n1_vertices;
  std::vector<Cartesian> orig_n2_vertices;
  std::vector<Cartesian> orig_colour_vector;
  int nstartarrow;
 public:
  ArrowHeadRibbon();
  ArrowHeadRibbon(const std::vector<Cartesian> &vertices_in, const std::vector<Cartesian> &n1_vertices_in, const std::vector<Cartesian> &n2_vertices_in, double *colour_in, const Cartesian &origin_in, const std::vector<Cartesian> &colour_vector_in, double minsize_in=0.2, double maxsize_in=1.0, double thickness_in=0.2, double arrow_length_in=2.2, double arrow_width_in=2.2, double alpha_in=1.0, int accu_in=8, int scale_steps=8, int textured_in=0);
  ~ArrowHeadRibbon();
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  std::vector<std::vector<Cartesian> > GetAB();
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  bool isLine() const {return false;};
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

enum { IMAGE_FULLSCREEN, IMAGE_ORIGINAL_SIZE, IMAGE_KEEP_ORIGINAL_SIZE };

class SimpleBillBoard {
  static GLdouble projMatrix[16];
  static GLdouble modelMatrix[16];
  static GLint viewport[4];
 public:
  virtual ~SimpleBillBoard(){};
  static double* GetProjectionMatrix();
  static double* GetModelviewMatrix();
  static GLint* GetViewport();
  static void SetMatrices();
  static int GetMagnification();
  virtual void draw(double *override_colour=0, int selective_override=0) = 0 ;
  virtual void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v) = 0;
  virtual void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v)=0;
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
};

static bool force_load_textures;
class BillBoard : public Primitive, public SimpleBillBoard {
   image_info image;
   unsigned int texture_id;
   unsigned int mask_texture_id;
   unsigned int texture_id_backup; // For when we are forced to get a new texture id (eg in different context).
   unsigned int mask_texture_id_backup; // For when we are forced to get a new texture id (eg in different context).
   int style;
   double first_scale_w;
   double first_scale_h;
   double scale_w;
   double scale_h;
   const char *filename;
 public:
   BillBoard();
   BillBoard(const std::vector<Cartesian> &vertices_in, const Cartesian &origin_in, const char* filename_in, int style_in=IMAGE_ORIGINAL_SIZE);
  ~BillBoard();
  virtual void SetScale(double scale_w_in, double scale_h_in);
  void set_draw_colour(GLfloat *col=0);
  virtual void draw(double *override_colour=0, int selective_override=0);
  virtual void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  virtual void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return false;};
};

class SimpleText : public Primitive{
 protected:
  double baseline;
  std::string text;
  std::string family;
  std::string weight;
  std::string slant;
  int fn_size;
  double yskip;
  int id;
  GLuint texture_id;
  GLuint texture_id_b;
  image_info texture;
  bool multicoloured;
 public:
  void initialize(void);
  SimpleText();
  SimpleText(const Cartesian &vertex_in, const GLuint tex_id, const GLuint tex_id_b, const int width, const int height, const Cartesian &origin_in, const std::string &text_in, const int size_in=14, double alpha_in=1.0);
  SimpleText(const Cartesian &vertex_in, const std::string &text_in, const Cartesian &origin_in, double alpha_in=1.0, const std::string &family_in="times", const int size_in=14, const std::string &weight_in="", const std::string &slant_in="");
  void renderStringToPixmap();
  const bool &isMultiColoured() const { return multicoloured;} ;
  const image_info &GetTexture() const { return texture;} ;
  const GLuint &GetTextureID() const { return texture_id;} ;
  const GLuint &GetTextureID_B() const { return texture_id_b;} ;
  virtual ~SimpleText();
  virtual void SetRasterPosition(double x_in, double y_in, double z_in) const = 0;
  void set_draw_colour(GLfloat *col=0);
  virtual void draw(double *override_colour=0, int selective_override=0) = 0;
  virtual int IsBillBoard() = 0;
  void draw_main(double *override_colour=0, int selective_override=0);
  void SetFontName(const std::string &family_in, const int size_in, const std::string &weight_in, const std::string &slant_in);
  std::string StripTags();
  int LoadFont();
  void BitMapFont(int count, int left, const MGFontInfo &finfo, int underline, int strikethrough);
  void BackspaceBitMapFont(int count, int left, const MGFontInfo &finfo);
  void BlankBitMapFont(int count, int left, const MGFontInfo &finfo);
  void SetFontSize(unsigned int fn_isize);
  void SetFontFamily(const std::string &Family);
  void SetFontWeight(const std::string &Weight);
  void SetFontSlant(const std::string &Slant);
  void SetText(const std::string &text_in);
  int GetFontSize() const;
  std::string GetFontFamily() const;
  std::string GetFontWeight() const;
  std::string GetFontSlant() const;
  std::string GetText(void) const;
  void SetDefaultColour(void) const;
  void SetColour(float r, float b, float g, float a) const;
  int GetID(void) const;
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  void DrawPSmain(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v, bool is_a_bill_board);
};

class Text : public SimpleText{
 public:
  Text();
  Text(const Cartesian &vertex_in, const GLuint tex_id, const GLuint tex_id_b, const int width, const int height, const Cartesian &origin_in, const std::string &text_in, const int size_in=14, double alpha_in=1.0);
  Text(const Cartesian &vertex_in, const std::string &text_in, const Cartesian &origin_in, double alpha_in=1.0, const std::string &family_in="times", const int size_in=14, const std::string &weight_in="", const std::string &slant_in="");
  ~Text();

  void draw(double *override_colour=0, int selective_override=0);

  void SetRasterPosition(double x_in, double y_in, double z_in) const ;
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  int IsBillBoard() {return 0;};
  bool isLine() const {return true;};
};

class BillBoardText : public SimpleText, public SimpleBillBoard {
 public:
  BillBoardText();
  BillBoardText(const Cartesian &vertex_in, const GLuint tex_id, const GLuint tex_id_b, const int width, const int height, const Cartesian &origin_in, const std::string &text_in, const int size_in=14, double alpha_in=1.0);
  BillBoardText(const Cartesian &vertex_in, const std::string &text_in, const Cartesian &origin_in, double alpha_in=1.0, const std::string &family_in="times", const int size_in=14, const std::string &weight_in="", const std::string &slant_in="");
  ~BillBoardText();

  void draw(double *override_colour=0, int selective_override=0);

  void SetRasterPosition(double x_in, double y_in, double z_in) const ;
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  int IsBillBoard() {return 1;};
  bool isLine() const {return true;};
};

class DashLineElement : public LineElement {
 public:
  DashLineElement();
  ~DashLineElement();
  DashLineElement(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return true;};
  void SetDashLength(double length) { dash_length = length; };
  double GetDashLength() const { return dash_length; };
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
 private:
  double dash_length;

};

class PointCollection : public Primitive {
  std::vector<Primitive*> lines;
 public:
  PointCollection();
  ~PointCollection();
  void add_primitive(Primitive *line);
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return true;};
  bool isDots() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  std::vector<Primitive*> GetPrimitives() const {return lines;}
  void SetPrimitives(std::vector<Primitive*> &lines_in){lines=lines_in;}
};

class LineCollection : public Primitive {
  std::vector<Primitive*> lines;
 public:
  LineCollection();
  ~LineCollection();
  void add_primitive(Primitive *line);
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return true;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  std::vector<Primitive*> GetPrimitives() const {return lines;}
  void SetPrimitives(std::vector<Primitive*> &lines_in){lines=lines_in;}
};

class Triangle: public TriangleElement{
 public:
  Triangle();
  ~Triangle(){};
  Triangle(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class MGPolygon: public PolygonElement{
 public:
  MGPolygon();
  ~MGPolygon(){};
  MGPolygon(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class Quad: public QuadElement{
 public:
  Quad();
  ~Quad(){};
  Quad(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in,  double alpha_in=1.0, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class Sphere: public SphereElement{
 public:
  Sphere();
  ~Sphere(){};
  Sphere(const Cartesian &vertex, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int accu_in=0, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class Spheroid: public SpheroidElement{
 public:
  Spheroid();
  ~Spheroid(){};
  Spheroid(const Cartesian &vertex, double *colour_in, const Cartesian &origin_in, double a_in=1.0, double b_in=1.0, double c_in=1.0, const Cartesian& xaxis_in=Cartesian(1,0,0), const Cartesian& yaxis_in=Cartesian(0,1,0), const Cartesian& zaxis_in=Cartesian(0,0,1), double alpha_in=1.0, int accu_in=0, int show_axes=0, int show_solid=1, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class Cylinder: public CylinderElement{
 public:
  Cylinder();
  ~Cylinder(){};
  Cylinder(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int accu_in=8, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class DashCylinder: public DashCylinderElement{
 public:
  DashCylinder();
  ~DashCylinder(){};
  DashCylinder(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int accu_in=8, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class Cone: public ConeElement{
 public:
  Cone();
  ~Cone(){};
  Cone(const std::vector<Cartesian> &vertices_in, double *colour_in, const Cartesian &origin_in, double size_in=1.0, double alpha_in=1.0, int accu_in=8, int textured_in=0);
  void draw(double *override_colour=0, int selective_override=0);
  bool isLine() const {return false;};
};

class PolyCollection : public Primitive {
  std::vector<Primitive*> prims;
 public:
   PolyCollection();
  ~ PolyCollection();
  void add_primitive(Primitive *prims);
  void set_draw_colour(GLfloat *col=0);
  void draw(double *override_colour=0, int selective_override=0);
  void DrawPovray(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, const Volume &v);
  void DrawPostScript(std::ofstream &fp, const Quat &quat, double radius, double ox, double oy, double oz, const matrix &objrotmatrix, const Cartesian &objorigin, double xoff, double yoff, double xscale, double yscale, double xscaleps, const Volume &v);
  bool isLine() const {return false;};
  std::vector<Primitive*> GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin);
  void SetAlpha(double alpha_in);
  void set_transparent(int trans_in);
  std::vector<Primitive*> GetPrimitives() const {return prims;}
  void SetPrimitives(std::vector<Primitive*> &prims_in){prims=prims_in;}
};

#endif
