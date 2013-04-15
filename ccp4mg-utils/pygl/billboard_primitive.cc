/*
     pygl/billboard_primitive.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York

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
#include "glew.h"
#include "wglew.h"
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

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#define PIBY2 (M_PI * 2)
#define ROOT3 1.7320508075688772
#define ROOT3OVER2 ROOT3*0.5

BillboardPrimitive::BillboardPrimitive() : Primitive(){
  x_texture_start=0;
  x_texture_end=1;
  y_texture_start=0;
  y_texture_end=1;
}

BillboardLabelledCircleElement::BillboardLabelledCircleElement() : BillboardPrimitive(){vertices.push_back(Cartesian());}
BillboardCylinderElement::BillboardCylinderElement(): BillboardPrimitive(){vertices.push_back(Cartesian());vertices.push_back(Cartesian());}
BillboardLabelledCircle::BillboardLabelledCircle() : BillboardLabelledCircleElement(){}
BillboardCylinder::BillboardCylinder(): BillboardCylinderElement(){}

BillboardPrimitive::~BillboardPrimitive(){
}

BillboardLabelledCircle::BillboardLabelledCircle(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in, const std::string text_in) : BillboardLabelledCircleElement(vertex,colour_in,origin_in,size_in,alpha_in,textured_in, text_in){
}

BillboardLabelledCircleElement::BillboardLabelledCircleElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in, const std::string text_in) : BillboardPrimitive(){
  size = size_in;
  alpha = alpha_in;
  vertices.clear();
  vertices.push_back(vertex);
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  textured = textured_in;
  text = text_in;
}

BillboardCylinder::BillboardCylinder(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in) : BillboardCylinderElement(vertices_in,colour_in,origin_in,size_in,alpha_in,textured_in){
}

BillboardCylinderElement::BillboardCylinderElement(const std::vector<Cartesian> &vertices_in, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in) : BillboardPrimitive(){
  size = size_in;
  alpha = alpha_in;
  vertices.clear();
  vertices = vertices_in;
  colour[0] = colour_in[0];
  colour[1] = colour_in[1];
  colour[2] = colour_in[2];
  origin = origin_in;
  textured = textured_in;
}

BillboardLabelledCircleElement::~BillboardLabelledCircleElement(){
}

BillboardLabelledCircle::~BillboardLabelledCircle(){
}

BillboardCylinder::~BillboardCylinder(){
}

BillboardCylinderElement::~BillboardCylinderElement(){
}

void BillboardCylinderElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void BillboardLabelledCircleElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void BillboardCylinderElement::draw(const double *override_colour, int selective_override){
}

void BillboardLabelledCircleElement::draw(const double *override_colour, int selective_override){
}

void BillboardCylinder::draw(const double *override_colour, int selective_override){
  BillboardCylinderElement::draw(override_colour,selective_override);
}

void BillboardLabelledCircle::draw(const double *override_colour, int selective_override){
  BillboardLabelledCircleElement::draw(override_colour,selective_override);
}

void BillboardCylinder::draw(double y, const Cartesian& x_trans, const Cartesian& y_trans, const Cartesian& z_trans){
  glEnable(GL_LIGHTING);
  BillboardCylinderElement::draw(y, x_trans, y_trans, z_trans);
}

void BillboardLabelledCircle::draw(double y, const Cartesian& x_trans, const Cartesian& y_trans, const Cartesian& z_trans){
  glEnable(GL_LIGHTING);
  BillboardLabelledCircleElement::draw(y, x_trans, y_trans, z_trans);
}

void BillboardCylinderElement::draw(double y, const Cartesian& x_trans, const Cartesian& y_trans, const Cartesian& z_trans){
      double x1 = GetVertices()[0].get_x();
      double y1 = GetVertices()[0].get_y();
      double z1 = GetVertices()[0].get_z();
      double x2 = GetVertices()[1].get_x();
      double y2 = GetVertices()[1].get_y();
      double z2 = GetVertices()[1].get_z();
      Cartesian v1v2 = vertices[1]-vertices[0];
      v1v2.normalize();
      double dot = Cartesian::DotProduct(z_trans,v1v2);
      double deform1 = ROOT3OVER2 * size * dot;
      double deform2 = size * dot;
      Cartesian x_cross = Cartesian::CrossProduct(z_trans,v1v2);
      Cartesian y_cross = Cartesian::CrossProduct(z_trans,x_cross);
      y_cross.normalize();

      // Try to not obscure the atom name by doing some fudging.
      x1 -= 0.07*dot * y_cross.get_x();
      y1 -= 0.07*dot * y_cross.get_y();
      z1 -= 0.07*dot * y_cross.get_z();

      Cartesian def = deform1 *(y_cross);
      Cartesian def2 = deform2 *(y_cross);
      x_cross.normalize(size);
      //std::cout << deform1 << " " << deform2 << "\n";



      double tsx = GetTextureStartXCoord();
      double tsy = GetTextureStartYCoord();
      double tex = GetTextureEndXCoord();
      double tey = GetTextureEndYCoord();

      if(y<0.5){
        glMultiTexCoord2f(GL_TEXTURE0,tsx,tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tsx+0.125,tsy);
        glVertex3f(x2+x_cross.get_x(),y2+x_cross.get_y(),z2+x_cross.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tex,tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tex+0.125,tsy);
        glVertex3f(x1+x_cross.get_x(),y1+x_cross.get_y(),z1+x_cross.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tex,0.25*tey+0.75*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tex+0.125,0.25*tey+0.75*tsy);
        glVertex3f(x1+0.5*x_cross.get_x()+def.get_x(),y1+0.5*x_cross.get_y()+def.get_y(),z1+0.5*x_cross.get_z()+def.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tsx,0.25*tey+0.75*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tsx+0.125,0.25*tey+0.75*tsy);
        glVertex3f(x2+0.5*x_cross.get_x()+def.get_x(),y2+0.5*x_cross.get_y()+def.get_y(),z2+0.5*x_cross.get_z()+def.get_z());

        glVertex3f(x2+0.5*x_cross.get_x()+def.get_x(),y2+0.5*x_cross.get_y()+def.get_y(),z2+0.5*x_cross.get_z()+def.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tex,0.25*tey+0.75*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tex+0.125,0.25*tey+0.75*tsy);
        glVertex3f(x1+0.5*x_cross.get_x()+def.get_x(),y1+0.5*x_cross.get_y()+def.get_y(),z1+0.5*x_cross.get_z()+def.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tex,0.5*tey+0.5*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tex+0.125,0.5*tey+0.5*tsy);
        glVertex3f(x1+def2.get_x(),y1+def2.get_y(),z1+def2.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tsx,0.5*tey+0.5*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tsx+0.125,0.5*tey+0.5*tsy);
        glVertex3f(x2+def2.get_x(),y2+def2.get_y(),z2+def2.get_z());

        glVertex3f(x2+def2.get_x(),y2+def2.get_y(),z2+def2.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tex,0.5*tey+0.5*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tex+0.125,0.5*tey+0.5*tsy);
        glVertex3f(x1+def2.get_x(),y1+def2.get_y(),z1+def2.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tex,0.75*tey+0.25*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tex+0.125,0.75*tey+0.25*tsy);
        glVertex3f(x1-0.5*x_cross.get_x()+def.get_x(),y1-0.5*x_cross.get_y()+def.get_y(),z1-0.5*x_cross.get_z()+def.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tsx,0.75*tey+0.25*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tsx+0.125,0.75*tey+0.25*tsy);
        glVertex3f(x2-0.5*x_cross.get_x()+def.get_x(),y2-0.5*x_cross.get_y()+def.get_y(),z2-0.5*x_cross.get_z()+def.get_z());

        glVertex3f(x2-0.5*x_cross.get_x()+def.get_x(),y2-0.5*x_cross.get_y()+def.get_y(),z2-0.5*x_cross.get_z()+def.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tex,0.75*tey+0.25*tsy);
        glMultiTexCoord2f(GL_TEXTURE1,tex+0.125,0.75*tey+0.25*tsy);
        glVertex3f(x1-0.5*x_cross.get_x()+def.get_x(),y1-0.5*x_cross.get_y()+def.get_y(),z1-0.5*x_cross.get_z()+def.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tex,tey);
        glMultiTexCoord2f(GL_TEXTURE1,tex+0.125,tey);
        glVertex3f(x1-x_cross.get_x(),y1-x_cross.get_y(),z1-x_cross.get_z());
        glMultiTexCoord2f(GL_TEXTURE0,tsx,tey);
        glMultiTexCoord2f(GL_TEXTURE1,tsx+0.125,tey);
        glVertex3f(x2-x_cross.get_x(),y2-x_cross.get_y(),z2-x_cross.get_z());
      } else {
      glTexCoord2f(tsx,tsy);
      glVertex3f(x2+x_cross.get_x(),y2+x_cross.get_y(),z2+x_cross.get_z());
      glTexCoord2f(tex,tsy);
      glVertex3f(x1+x_cross.get_x(),y1+x_cross.get_y(),z1+x_cross.get_z());
      glTexCoord2f(tex,0.25*tey+0.75*tsy);
      glVertex3f(x1+0.5*x_cross.get_x()+def.get_x(),y1+0.5*x_cross.get_y()+def.get_y(),z1+0.5*x_cross.get_z()+def.get_z());
      glTexCoord2f(tsx,0.25*tey+0.75*tsy);
      glVertex3f(x2+0.5*x_cross.get_x()+def.get_x(),y2+0.5*x_cross.get_y()+def.get_y(),z2+0.5*x_cross.get_z()+def.get_z());

      glVertex3f(x2+0.5*x_cross.get_x()+def.get_x(),y2+0.5*x_cross.get_y()+def.get_y(),z2+0.5*x_cross.get_z()+def.get_z());
      glTexCoord2f(tex,0.25*tey+0.75*tsy);
      glVertex3f(x1+0.5*x_cross.get_x()+def.get_x(),y1+0.5*x_cross.get_y()+def.get_y(),z1+0.5*x_cross.get_z()+def.get_z());
      glTexCoord2f(tex,0.5*tey+0.5*tsy);
      glVertex3f(x1+def2.get_x(),y1+def2.get_y(),z1+def2.get_z());
      glTexCoord2f(tsx,0.5*tey+0.5*tsy);
      glVertex3f(x2+def2.get_x(),y2+def2.get_y(),z2+def2.get_z());

      glVertex3f(x2+def2.get_x(),y2+def2.get_y(),z2+def2.get_z());
      glTexCoord2f(tex,0.5*tey+0.5*tsy);
      glVertex3f(x1+def2.get_x(),y1+def2.get_y(),z1+def2.get_z());
      glTexCoord2f(tex,0.75*tey+0.25*tsy);
      glVertex3f(x1-0.5*x_cross.get_x()+def.get_x(),y1-0.5*x_cross.get_y()+def.get_y(),z1-0.5*x_cross.get_z()+def.get_z());
      glTexCoord2f(tsx,0.75*tey+0.25*tsy);
      glVertex3f(x2-0.5*x_cross.get_x()+def.get_x(),y2-0.5*x_cross.get_y()+def.get_y(),z2-0.5*x_cross.get_z()+def.get_z());

      glVertex3f(x2-0.5*x_cross.get_x()+def.get_x(),y2-0.5*x_cross.get_y()+def.get_y(),z2-0.5*x_cross.get_z()+def.get_z());
      glTexCoord2f(tex,0.75*tey+0.25*tsy);
      glVertex3f(x1-0.5*x_cross.get_x()+def.get_x(),y1-0.5*x_cross.get_y()+def.get_y(),z1-0.5*x_cross.get_z()+def.get_z());
      glTexCoord2f(tex,tey);
      glVertex3f(x1-x_cross.get_x(),y1-x_cross.get_y(),z1-x_cross.get_z());
      glTexCoord2f(tsx,tey);
      glVertex3f(x2-x_cross.get_x(),y2-x_cross.get_y(),z2-x_cross.get_z());
      }

}

void BillboardLabelledCircleElement::draw(double y, const Cartesian& x_trans, const Cartesian& y_trans, const Cartesian& z_trans){
  /*
    We use textured quads for speed. Much faster than drawing circles and we have to use a texture
    for the text anyway.
  */
  //std::cout << "BillboardLabelledCircleElement::draw\n";
  if(y<0.5){
      double size = GetSize();
      double x = GetVertices()[0].get_x();
      double y = GetVertices()[0].get_y();
      double z = GetVertices()[0].get_z();
      glMultiTexCoord2f(GL_TEXTURE0,GetTextureStartXCoord(),GetTextureStartYCoord());
      glMultiTexCoord2f(GL_TEXTURE1,0.875,0);
      glVertex3f(x-size*x_trans.get_x()-size*y_trans.get_x(),y-size*x_trans.get_y()-size*y_trans.get_y(),z-size*x_trans.get_z()-size*y_trans.get_z());
      glMultiTexCoord2f(GL_TEXTURE0,GetTextureEndXCoord(),GetTextureStartYCoord());
      glMultiTexCoord2f(GL_TEXTURE1,1,0);
      glVertex3f(x+size*x_trans.get_x()-size*y_trans.get_x(),y+size*x_trans.get_y()-size*y_trans.get_y(),z+size*x_trans.get_z()-size*y_trans.get_z());
      glMultiTexCoord2f(GL_TEXTURE0,GetTextureEndXCoord(),GetTextureEndYCoord());
      glMultiTexCoord2f(GL_TEXTURE1,1,0.125);
      glVertex3f(x+size*x_trans.get_x()+size*y_trans.get_x(),y+size*x_trans.get_y()+size*y_trans.get_y(),z+size*x_trans.get_z()+size*y_trans.get_z());
      glMultiTexCoord2f(GL_TEXTURE0,GetTextureStartXCoord(),GetTextureEndYCoord());
      glMultiTexCoord2f(GL_TEXTURE1,0.875,0.125);
      glVertex3f(x-size*x_trans.get_x()+size*y_trans.get_x(),y-size*x_trans.get_y()+size*y_trans.get_y(),z-size*x_trans.get_z()+size*y_trans.get_z());
  } else {
      double size = GetSize();
      double x = GetVertices()[0].get_x();
      double y = GetVertices()[0].get_y();
      double z = GetVertices()[0].get_z();
      glTexCoord2f(GetTextureStartXCoord(),GetTextureStartYCoord());
      glVertex3f(x-size*x_trans.get_x()-size*y_trans.get_x(),y-size*x_trans.get_y()-size*y_trans.get_y(),z-size*x_trans.get_z()-size*y_trans.get_z());
      glTexCoord2f(GetTextureEndXCoord(),GetTextureStartYCoord());
      glVertex3f(x+size*x_trans.get_x()-size*y_trans.get_x(),y+size*x_trans.get_y()-size*y_trans.get_y(),z+size*x_trans.get_z()-size*y_trans.get_z());
      glTexCoord2f(GetTextureEndXCoord(),GetTextureEndYCoord());
      glVertex3f(x+size*x_trans.get_x()+size*y_trans.get_x(),y+size*x_trans.get_y()+size*y_trans.get_y(),z+size*x_trans.get_z()+size*y_trans.get_z());
      glTexCoord2f(GetTextureStartXCoord(),GetTextureEndYCoord());
      glVertex3f(x-size*x_trans.get_x()+size*y_trans.get_x(),y-size*x_trans.get_y()+size*y_trans.get_y(),z-size*x_trans.get_z()+size*y_trans.get_z());

  }
}

std::vector<Primitive*> BillboardCylinderElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "BillboardCylinder::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  BillboardCylinderElement *tri = new BillboardCylinder(vertices,colour,origin,size,alpha,textured);
      tri->SetTextureStartXCoord(GetTextureStartXCoord());
      tri->SetTextureEndXCoord(GetTextureEndXCoord());
      tri->SetTextureStartYCoord(GetTextureStartYCoord());
      tri->SetTextureEndYCoord(GetTextureEndYCoord());
  a.push_back(tri);
  return a;
}

std::vector<Primitive*> BillboardLabelledCircleElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  //std::cout << "BillboardLabelledCircle::GetSimplePrimitives\n";
  std::vector<Primitive*> a;
  BillboardLabelledCircleElement *tri = new BillboardLabelledCircle(vertices[0],colour,origin,size,alpha,textured);
      tri->SetTextureStartXCoord(GetTextureStartXCoord());
      tri->SetTextureEndXCoord(GetTextureEndXCoord());
      tri->SetTextureStartYCoord(GetTextureStartYCoord());
      tri->SetTextureEndYCoord(GetTextureEndYCoord());
  a.push_back(tri);
  return a;
}

ImposterSphereElement::ImposterSphereElement() : BillboardPrimitive(){vertices.push_back(Cartesian());}
ImposterSphere::ImposterSphere() : ImposterSphereElement(){}

ImposterSphere::ImposterSphere(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in) : ImposterSphereElement(vertex,colour_in,origin_in,size_in,alpha_in,textured_in){
}

ImposterSphereElement::ImposterSphereElement(const Cartesian &vertex, const double *colour_in, const Cartesian &origin_in, double size_in, double alpha_in, int textured_in) : BillboardPrimitive(){
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

ImposterSphereElement::~ImposterSphereElement(){
}

ImposterSphere::~ImposterSphere(){
}

void ImposterSphereElement::set_draw_colour(const GLfloat *col){
  if(col)
    set_draw_colour_poly_override(col);
  else
    set_draw_colour_poly();
}

void ImposterSphereElement::draw(const double *override_colour, int selective_override){
}

void ImposterSphere::draw(const double *override_colour, int selective_override){
  ImposterSphereElement::draw(override_colour,selective_override);
}

void ImposterSphere::draw(double y, const Cartesian& x_trans, const Cartesian& y_trans, const Cartesian& z_trans){
  glEnable(GL_LIGHTING);
  ImposterSphereElement::draw(y, x_trans, y_trans, z_trans);
}

void ImposterSphereElement::draw(double ydum, const Cartesian& x_trans, const Cartesian& y_trans, const Cartesian& z_trans){
      double size = GetSize();
      double x = GetVertices()[0].get_x();
      double y = GetVertices()[0].get_y();
      double z = GetVertices()[0].get_z();
      glTexCoord2f(0,0);
      glVertex3f(x-size*x_trans.get_x()-size*y_trans.get_x(),y-size*x_trans.get_y()-size*y_trans.get_y(),z-size*x_trans.get_z()-size*y_trans.get_z());
      glTexCoord2f(1,0);
      glVertex3f(x+size*x_trans.get_x()-size*y_trans.get_x(),y+size*x_trans.get_y()-size*y_trans.get_y(),z+size*x_trans.get_z()-size*y_trans.get_z());
      glTexCoord2f(1,1);
      glVertex3f(x+size*x_trans.get_x()+size*y_trans.get_x(),y+size*x_trans.get_y()+size*y_trans.get_y(),z+size*x_trans.get_z()+size*y_trans.get_z());
      glTexCoord2f(0,1);
      glVertex3f(x-size*x_trans.get_x()+size*y_trans.get_x(),y-size*x_trans.get_y()+size*y_trans.get_y(),z-size*x_trans.get_z()+size*y_trans.get_z());
}

std::vector<Primitive*> ImposterSphereElement::GetSimplePrimitives(const Volume &clip_vol, const matrix &objrotmatrix, const Cartesian &objorigin, int start, int end) const {
  std::vector<Primitive*> a;
  SphereElement *tri = new SphereElement(vertices[0],colour,origin,size,alpha,textured);
  a.push_back(tri);
  return a;
}

