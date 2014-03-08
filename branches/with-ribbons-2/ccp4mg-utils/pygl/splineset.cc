/*
     pygl/splineset.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York
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

#ifdef _WIN32
#include <windows.h>
#endif

#define GLUT_DISABLE_ATEXIT_HACK
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "help.h"
#include "texture.h"
#include <vector>
#include "cartesian.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#define PIBY2 (M_PI * 2)

#if defined (_WIN32) && not defined (WINDOWS_MINGW)
#define EXAMPLE_DLL __declspec(dllexport)
void __stdcall EXAMPLE_DLL draw_flat_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false, const bool grey_ribbon_edge=false);
void __stdcall EXAMPLE_DLL draw_flat_rounded_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
void __stdcall EXAMPLE_DLL draw_fancy_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
void __stdcall EXAMPLE_DLL draw_elliptical_ribbon(const std::vector<Cartesian> &vertices, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int quality, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
#else
void draw_flat_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false, const bool grey_ribbon_edge=false);
void draw_flat_rounded_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
void draw_fancy_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
void draw_elliptical_ribbon(const std::vector<Cartesian> &vertices, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int quality, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour=false);
#endif

bool GetFaceOne(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pvpr){
  bool face_one = true;
  if(spline.size()>3){
  unsigned idx = spline.size()/2;
  Cartesian p_t = spline[idx+1] + spline[idx-1] - 2*spline[idx];
  while(fabs(p_t.length())<1e-8&&idx>spline.size()+1){
    idx--;
    p_t = spline[idx+1] + spline[idx-1] - 2*spline[idx];
  }
  if(fabs(p_t.length())>1e-8){
    p_t.normalize();
    Cartesian pvpr_t = pvpr[idx];
    pvpr_t.normalize();

    face_one = Cartesian::DotProduct(p_t,pvpr_t)>0.0;
  }
  }
  return face_one;
}

void draw_ellipse_point(double theta, const Cartesian &pv, const Cartesian &pvpr, const Cartesian &spline);

std::vector<Cartesian> get_v_and_vpr(const Cartesian &sp1, const Cartesian &sp2, const Cartesian &ovec){

  Cartesian  n = sp1 - sp2;
  n.normalize();

  Cartesian vpr = n.CrossProduct(ovec,n); 
  vpr.normalize();

  Cartesian v = n.CrossProduct(vpr,n);

  std::vector<Cartesian> v_vpr;
  v_vpr.push_back(v);
  v_vpr.push_back(vpr);

  return v_vpr;
}


void draw_flat_rounded_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int nsectors, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour){

  Cartesian normal;

  Cartesian p1,p2;
  unsigned int i;

  double frac=0.0;

  bool face_one = GetFaceOne(spline,pvpr);

  double thickness = pvpr[0].length(); // Assume constant thickness.

  glBegin(GL_TRIANGLE_STRIP);
  if(face_one&&two_colour){
    GLfloat colour[4] = {0.6,0.6,0.6,1.0};
    glColor4fv(colour);
  }
  for(i=0;i<spline.size();i++){
    if((multicolour&&!face_one)||!two_colour){
      GLfloat colour[4] = {colour_vector[i].get_x(),colour_vector[i].get_y(),colour_vector[i].get_z(),colour_vector[i].get_a()};
      glColor4fv(colour);
    }
 
    Cartesian p = pv[i];
    p.normalize(p.length()-thickness);

    p2 = spline[i] - p + pvpr[i];
    p1 = spline[i] + p + pvpr[i];

    normal = pvpr[i];
    normal.normalize();
    glNormal3d(normal.get_x(),normal.get_y(),normal.get_z());
    if(textured){
      frac = (double)i/(double)(spline.size());
      set_texture_coord(frac, 0.0, textured);
    }
    glVertex3d(p1.get_x(),p1.get_y(),p1.get_z());
    if(textured){
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3d(p2.get_x(),p2.get_y(),p2.get_z());
  }
  glEnd();

  double theta,theta2;
  double tex_xfrac=0.0;

  glBegin(GL_TRIANGLE_STRIP);

  if(!face_one&&two_colour){
    GLfloat colour[4] = {0.6,0.6,0.6,1.0};
    glColor4fv(colour);
  }
  for(i=0;i<spline.size();i++){
    if((multicolour&&face_one)||!two_colour){
      GLfloat colour[4] = {colour_vector[i].get_x(),colour_vector[i].get_y(),colour_vector[i].get_z(),colour_vector[i].get_a()};
      glColor4fv(colour);
    }

    Cartesian p = pv[i];
    p.normalize(p.length()-thickness);

    p2 = spline[i] + p - pvpr[i];
    p1 = spline[i] - p - pvpr[i];

    normal = -pvpr[i];
    normal.normalize();
    glNormal3d(normal.get_x(),normal.get_y(),normal.get_z());
    if(textured){
      frac = (double)i/(double)(spline.size());
      set_texture_coord(frac, 0.0, textured);
    }
    glVertex3d(p1.get_x(),p1.get_y(),p1.get_z());
    if(textured){
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3d(p2.get_x(),p2.get_y(),p2.get_z());
  }
  glEnd();

  for(int ii=0;ii<360;ii=ii+360/nsectors){
    theta = (double)ii/360.0 * PIBY2;
    theta2 = (double)(ii-360/nsectors)/360.0 * PIBY2;
    glBegin(GL_TRIANGLE_STRIP);
    for(unsigned int j=0;j<spline.size();j++){
      if(multicolour){
	GLfloat colour[4] = {colour_vector[j].get_x(),colour_vector[j].get_y(),colour_vector[j].get_z(),colour_vector[j].get_a()};
        glColor4fv(colour);
      }

      Cartesian p = pv[j];
      p.normalize(p.length()-thickness);

      p1 = spline[j] - p;

      Cartesian ppr = pvpr[j];
      ppr.normalize(thickness);

      p = pv[j];
      p.normalize(thickness);

      if(textured){
        tex_xfrac = (double)j/(double)(spline.size());
	set_texture_coord(tex_xfrac, theta/PIBY2, textured);
      }
      draw_ellipse_point(theta, ppr, p, p1);
      if(textured){
	set_texture_coord(tex_xfrac, theta2/PIBY2, textured);
      }
      draw_ellipse_point(theta2, ppr, p, p1);

    }
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    for(unsigned int j=0;j<spline.size();j++){
      if(multicolour){
	GLfloat colour[4] = {colour_vector[j].get_x(),colour_vector[j].get_y(),colour_vector[j].get_z(),colour_vector[j].get_a()};
        glColor4fv(colour);
      }

      Cartesian p = pv[j];
      p.normalize(p.length()-thickness);

      p1 = spline[j] + p;

      Cartesian ppr = pvpr[j];
      ppr.normalize(thickness);

      p = pv[j];
      p.normalize(thickness);

      if(textured){
        tex_xfrac = (double)j/(double)(spline.size());
	set_texture_coord(tex_xfrac, theta/PIBY2, textured);
      }
      draw_ellipse_point(theta, ppr, p, p1);
      if(textured){
	set_texture_coord(tex_xfrac, theta2/PIBY2, textured);
      }
      draw_ellipse_point(theta2, ppr, p, p1);

    }
    glEnd();

  }

}

void draw_fancy_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int nsectors, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour){

  Cartesian normal;

  Cartesian p1,p2;
  unsigned int i;

  double frac=0.0;

  double thickness = pvpr[0].length(); // Assume constant thickness.

  bool face_one = GetFaceOne(spline,pvpr);

  glBegin(GL_TRIANGLE_STRIP);
  if(face_one&&two_colour){
    GLfloat colour[4] = {0.6,0.6,0.6,1.0};
    glColor4fv(colour);
  }
  for(i=0;i<spline.size();i++){
    if((multicolour&&!face_one)||!two_colour){
      GLfloat colour[4] = {colour_vector[i].get_x(),colour_vector[i].get_y(),colour_vector[i].get_z(),colour_vector[i].get_a()};
      glColor4fv(colour);
    }

    Cartesian p = pv[i];
    p.normalize(p.length()-thickness);

    p2 = spline[i] - p + pvpr[i]*.7;
    p1 = spline[i] + p + pvpr[i]*.7;

    normal = pvpr[i];
    normal.normalize();
    glNormal3d(normal.get_x(),normal.get_y(),normal.get_z());
    if(textured){
      frac = (double)i/(double)(spline.size());
      set_texture_coord(frac, 0.0, textured);
    }
    glVertex3d(p1.get_x(),p1.get_y(),p1.get_z());
    if(textured){
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3d(p2.get_x(),p2.get_y(),p2.get_z());
  }
  glEnd();

  double theta,theta2;
  double tex_xfrac=0.0;

  glBegin(GL_TRIANGLE_STRIP);
  if(!face_one&&two_colour){
    GLfloat colour[4] = {0.6,0.6,0.6,1.0};
    glColor4fv(colour);
  }
  for(i=0;i<spline.size();i++){
    if((multicolour&&face_one)||!two_colour){
      GLfloat colour[4] = {colour_vector[i].get_x(),colour_vector[i].get_y(),colour_vector[i].get_z(),colour_vector[i].get_a()};
      glColor4fv(colour);
    }

    Cartesian p = pv[i];
    p.normalize(p.length()-thickness);

    p2 = spline[i] + p - pvpr[i]/2.0;
    p1 = spline[i] - p - pvpr[i]/2.0;

    normal = -pvpr[i];
    normal.normalize();
    glNormal3d(normal.get_x(),normal.get_y(),normal.get_z());
    if(textured){
      frac = (double)i/(double)(spline.size());
      set_texture_coord(frac, 0.0, textured);
    }
    glVertex3d(p1.get_x(),p1.get_y(),p1.get_z());
    if(textured){
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3d(p2.get_x(),p2.get_y(),p2.get_z());
  }
  glEnd();

  for(int ii=0;ii<360;ii=ii+360/nsectors){
    theta = (double)ii/360.0 * PIBY2;
    theta2 = (double)(ii-360/nsectors)/360.0 * PIBY2;
    glBegin(GL_TRIANGLE_STRIP);
    for(unsigned int j=0;j<spline.size();j++){
      if(multicolour){
	GLfloat colour[4] = {colour_vector[j].get_x(),colour_vector[j].get_y(),colour_vector[j].get_z(),colour_vector[j].get_a()};
        glColor4fv(colour);
      }

      Cartesian p = pv[j];
      p.normalize(p.length()-thickness);

      p1 = spline[j] - p;

      Cartesian ppr = pvpr[j];
      ppr.normalize(thickness);

      p = pv[j];
      p.normalize(thickness);

      if(textured){
        tex_xfrac = (double)j/(double)(spline.size());
	set_texture_coord(tex_xfrac, theta/PIBY2, textured);
      }
      draw_ellipse_point(theta, ppr, p, p1);
      if(textured){
	set_texture_coord(tex_xfrac, theta2/PIBY2, textured);
      }
      draw_ellipse_point(theta2, ppr, p, p1);

    }
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    for(unsigned int j=0;j<spline.size();j++){
      if(multicolour){
	GLfloat colour[4] = {colour_vector[j].get_x(),colour_vector[j].get_y(),colour_vector[j].get_z(),colour_vector[j].get_a()};
        glColor4fv(colour);
      }

      Cartesian p = pv[j];
      p.normalize(p.length()-thickness);

      p1 = spline[j] + p;

      Cartesian ppr = pvpr[j];
      ppr.normalize(thickness);

      p = pv[j];
      p.normalize(thickness);

      if(textured){
        tex_xfrac = (double)j/(double)(spline.size());
	set_texture_coord(tex_xfrac, theta/PIBY2, textured);
      }
      draw_ellipse_point(theta, ppr, p, p1);
      if(textured){
	set_texture_coord(tex_xfrac, theta2/PIBY2, textured);
      }
      draw_ellipse_point(theta2, ppr, p, p1);

    }
    glEnd();

  }

}

void draw_flat_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int npoints, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour, const bool grey_ribbon_edge){

  /* Construct the four quadstrips in turn */

  Cartesian normal;

  Cartesian p1,p2;
  unsigned int i;

  double frac=0.0;

  bool face_one = GetFaceOne(spline,pvpr);

  if(multicolour){
    GLfloat colour[4] = {colour_vector[0].get_x(),colour_vector[0].get_y(),colour_vector[0].get_z(),colour_vector[0].get_a()};
    glColor4fv(colour);
  }
  Cartesian n = spline[0] - spline[1];
  GLfloat ndp[3] = {n.get_x(), n.get_y(), n.get_z() };
  glNormal3fv(ndp);
  Cartesian s1 = spline[0]+pv[0]+pvpr[0];
  Cartesian s2 = spline[0]+pv[0]-pvpr[0];
  Cartesian s3 = spline[0]-pv[0]-pvpr[0];
  Cartesian s4 = spline[0]-pv[0]+pvpr[0];
  glBegin(GL_QUADS);
  glVertex3d(s1.get_x(),s1.get_y(),s1.get_z());
  glVertex3d(s2.get_x(),s2.get_y(),s2.get_z());
  glVertex3d(s3.get_x(),s3.get_y(),s3.get_z());
  glVertex3d(s4.get_x(),s4.get_y(),s4.get_z());
  glEnd();
  if(multicolour){
    GLfloat colour[4] = {colour_vector[colour_vector.size()-1].get_x(),colour_vector[colour_vector.size()-1].get_y(),colour_vector[colour_vector.size()-1].get_z(),colour_vector[colour_vector.size()-1].get_a()};
    glColor4fv(colour);
  }
  n = spline[spline.size()-1] - spline[spline.size()-2];
  ndp[0] = n.get_x();
  ndp[1] = n.get_y();
  ndp[2] = n.get_z();
  glNormal3fv(ndp);
  s1 = spline[spline.size()-1]+pv[pv.size()-1]+pvpr[pvpr.size()-1];
  s2 = spline[spline.size()-1]-pv[pv.size()-1]+pvpr[pvpr.size()-1];
  s3 = spline[spline.size()-1]-pv[pv.size()-1]-pvpr[pvpr.size()-1];
  s4 = spline[spline.size()-1]+pv[pv.size()-1]-pvpr[pvpr.size()-1];
  glBegin(GL_QUADS);
  glVertex3d(s1.get_x(),s1.get_y(),s1.get_z());
  glVertex3d(s2.get_x(),s2.get_y(),s2.get_z());
  glVertex3d(s3.get_x(),s3.get_y(),s3.get_z());
  glVertex3d(s4.get_x(),s4.get_y(),s4.get_z());
  glEnd();

  glBegin(GL_TRIANGLE_STRIP);
  if(face_one&&two_colour){
    GLfloat colour[4] = {0.6,0.6,0.6,1.0};
    glColor4fv(colour);
  }
  for(i=0;i<spline.size();i++){
    if((multicolour&&!face_one)||!two_colour){
      GLfloat colour[4] = {colour_vector[i].get_x(),colour_vector[i].get_y(),colour_vector[i].get_z(),colour_vector[i].get_a()};
      glColor4fv(colour);
    }
    p2 = spline[i] - pv[i] + pvpr[i];
    p1 = spline[i] + pv[i] + pvpr[i];

    normal = pvpr[i];
    normal.normalize();
    glNormal3d(normal.get_x(),normal.get_y(),normal.get_z());
    if(textured){
      frac = (double)i/(double)(spline.size());
      set_texture_coord(frac, 0.0, textured);
    }
    glVertex3d(p1.get_x(),p1.get_y(),p1.get_z());
    if(textured){
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3d(p2.get_x(),p2.get_y(),p2.get_z());
  }
  glEnd();

  glBegin(GL_TRIANGLE_STRIP);
  if(!face_one&&two_colour){
    GLfloat colour[4] = {0.6,0.6,0.6,1.0};
    glColor4fv(colour);
  }
  for(i=0;i<spline.size();i++){
    if((multicolour&&face_one)||!two_colour){
      GLfloat colour[4] = {colour_vector[i].get_x(),colour_vector[i].get_y(),colour_vector[i].get_z(),colour_vector[i].get_a()};
      glColor4fv(colour);
    }

    p2 = spline[i] + pv[i] - pvpr[i];
    p1 = spline[i] - pv[i] - pvpr[i];

    normal = -pvpr[i];
    normal.normalize();
    glNormal3d(normal.get_x(),normal.get_y(),normal.get_z());
    if(textured){
      frac = (double)i/(double)(spline.size());
      set_texture_coord(frac, 0.0, textured);
    }
    glVertex3d(p1.get_x(),p1.get_y(),p1.get_z());
    if(textured){
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3d(p2.get_x(),p2.get_y(),p2.get_z());
  }
  glEnd();
  glBegin(GL_TRIANGLE_STRIP);

  for(i=0;i<spline.size();i++){
    if(multicolour){
      GLfloat colour[4] = {colour_vector[i].get_x(),colour_vector[i].get_y(),colour_vector[i].get_z(),colour_vector[i].get_a()};
      glColor4fv(colour);
    }
    p2 = spline[i] + pv[i] + pvpr[i];
    p1 = spline[i] + pv[i] - pvpr[i];

    normal = pv[i];
    normal.normalize();
    glNormal3d(normal.get_x(),normal.get_y(),normal.get_z());
    if(textured){
      frac = (double)i/(double)(spline.size());
      set_texture_coord(frac, 0.0, textured);
    }
    glVertex3d(p1.get_x(),p1.get_y(),p1.get_z());
    if(textured){
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3d(p2.get_x(),p2.get_y(),p2.get_z());
  }
  glEnd();
  glBegin(GL_TRIANGLE_STRIP);

  for(i=0;i<spline.size();i++){
    if(multicolour){
      GLfloat colour[4] = {colour_vector[i].get_x(),colour_vector[i].get_y(),colour_vector[i].get_z(),colour_vector[i].get_a()};
      glColor4fv(colour);
    }
    p2 = spline[i] - pv[i] - pvpr[i];
    p1 = spline[i] - pv[i] + pvpr[i];

    normal = -pv[i];
    normal.normalize();
    glNormal3d(normal.get_x(),normal.get_y(),normal.get_z());
    if(textured){
      frac = (double)i/(double)(spline.size());
      set_texture_coord(frac, 0.0, textured);
    }
    glVertex3d(p1.get_x(),p1.get_y(),p1.get_z());
    if(textured){
      set_texture_coord(frac, 1.0, textured);
    }
    glVertex3d(p2.get_x(),p2.get_y(),p2.get_z());
  }
  glEnd();

}

void draw_ellipse_point(double theta, const Cartesian &pv, const Cartesian &pvpr, const Cartesian &spline){


  double a = pv.length();
  double b = pvpr.length();

  double x = cos(theta);
  double y = sin(theta);

  Cartesian p = x*pv + y*pvpr + spline;
  Cartesian n = pv*(x/(a*a)) + pvpr*(y/(b*b));
  GLfloat ndp[3] = {n.get_x(), n.get_y(), n.get_z() };

  glNormal3fv(ndp);
  glVertex3f(p.get_x(), p.get_y(), p.get_z());

}

void draw_elliptical_ribbon(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int nsectors, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour){

  double theta,theta2;
  double tex_xfrac=0.0;

  if(spline.size()>1){
  if(multicolour){
    GLfloat colour[4] = {colour_vector[0].get_x(),colour_vector[0].get_y(),colour_vector[0].get_z(),colour_vector[0].get_a()};
    glColor4fv(colour);
  }
  Cartesian n = spline[0] - spline[1];
  GLfloat ndp[3] = {n.get_x(), n.get_y(), n.get_z() };
  glNormal3fv(ndp);
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(spline[0].get_x(),spline[0].get_y(),spline[0].get_z());
  for(int i=0;i<360+360/nsectors;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    Cartesian p = x*pvpr[0] + y*pv[0] + spline[0];
    glVertex3f(p.get_x(), p.get_y(), p.get_z());
  }
  glEnd();
  }

  if(spline.size()>3){
  if(multicolour){
    GLfloat colour[4] = {colour_vector[colour_vector.size()-1].get_x(),colour_vector[colour_vector.size()-1].get_y(),colour_vector[colour_vector.size()-1].get_z(),colour_vector[colour_vector.size()-1].get_a()};
    glColor4fv(colour);
  }
  Cartesian n = spline[spline.size()-1] - spline[spline.size()-2];
  GLfloat ndp[3] = {n.get_x(), n.get_y(), n.get_z() };
  glNormal3fv(ndp);
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(spline[spline.size()-1].get_x(),spline[spline.size()-1].get_y(),spline[spline.size()-1].get_z());
  for(int i=0;i<360+360/nsectors;i=i+360/nsectors){
    theta = (double)i/360.0 * -PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    Cartesian p = x*pvpr[pvpr.size()-1] + y*pv[pv.size()-1] + spline[spline.size()-1];
    glVertex3f(p.get_x(), p.get_y(), p.get_z());
  }
  glEnd();
  }

  bool face_one = GetFaceOne(spline,pvpr);

  GLfloat grey[4] = {0.6,0.6,0.6,1.0};
  for(int i=0;i<360;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    theta2 = (double)(i-360.0/nsectors)/360.0 * PIBY2;
    glBegin(GL_TRIANGLE_STRIP);
    if(two_colour) glColor4fv(grey);
    for(unsigned int j=0;j<spline.size();j++){
      if(multicolour&&(!two_colour||(face_one&&((i-360/nsectors)>=90&&(i-360/nsectors)<270))||(!face_one&&((i-360/nsectors)<90||(i-360/nsectors)>=270)))){
	GLfloat colour[4] = {colour_vector[j].get_x(),colour_vector[j].get_y(),colour_vector[j].get_z(),colour_vector[j].get_a()};
        glColor4fv(colour);
      }

      if(textured){
        tex_xfrac = (double)j/(double)(spline.size());
	set_texture_coord(tex_xfrac, theta/PIBY2, textured);
      }
      draw_ellipse_point(theta, pvpr[j], pv[j], spline[j]);
      if(textured){
	set_texture_coord(tex_xfrac, theta2/PIBY2, textured);
      }
      draw_ellipse_point(theta2, pvpr[j], pv[j], spline[j]);

    }
    glEnd();

  }

}


void set_ellipse_point(double theta, const Cartesian &pv, const Cartesian &pvpr, const Cartesian &spline,GLfloat *vertex, GLfloat *normal){


  double a = pv.length();
  double b = pvpr.length();

  double x = cos(theta);
  double y = sin(theta);

  Cartesian p = x*pv + y*pvpr + spline;
  Cartesian n = pv*(x/(a*a)) + pvpr*(y/(b*b));

  normal[0] = n.get_x();
  normal[1] = n.get_y();
  normal[2] = n.get_z();
  vertex[0] = p.get_x();
  vertex[1] = p.get_y();
  vertex[2] = p.get_z();

}


void draw_elliptical_ribbon_vertex_arrays(const std::vector<Cartesian> &spline, const std::vector<Cartesian> &pv, const std::vector<Cartesian> &pvpr, int nsectors, int textured, int multicolour, const std::vector<Cartesian> &colour_vector, const bool two_colour){
  /* This is experimental to see if VBOs can give speed up with ribbons. */

  double theta,theta2;
  double tex_xfrac=0.0;

  if(spline.size()>1){
  if(multicolour){
    GLfloat colour[4] = {colour_vector[0].get_x(),colour_vector[0].get_y(),colour_vector[0].get_z(),colour_vector[0].get_a()};
    glColor4fv(colour);
  }
  Cartesian n = spline[0] - spline[1];
  GLfloat ndp[3] = {n.get_x(), n.get_y(), n.get_z() };
  glNormal3fv(ndp);
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(spline[0].get_x(),spline[0].get_y(),spline[0].get_z());
  for(int i=0;i<360+360/nsectors;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    Cartesian p = x*pvpr[0] + y*pv[0] + spline[0];
    glVertex3f(p.get_x(), p.get_y(), p.get_z());
  }
  glEnd();
  }

  if(spline.size()>3){
  if(multicolour){
    GLfloat colour[4] = {colour_vector[colour_vector.size()-1].get_x(),colour_vector[colour_vector.size()-1].get_y(),colour_vector[colour_vector.size()-1].get_z(),colour_vector[colour_vector.size()-1].get_a()};
    glColor4fv(colour);
  }
  Cartesian n = spline[spline.size()-1] - spline[spline.size()-2];
  GLfloat ndp[3] = {n.get_x(), n.get_y(), n.get_z() };
  glNormal3fv(ndp);
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(spline[spline.size()-1].get_x(),spline[spline.size()-1].get_y(),spline[spline.size()-1].get_z());
  for(int i=0;i<360+360/nsectors;i=i+360/nsectors){
    theta = (double)i/360.0 * -PIBY2;
    double x = cos(theta);
    double y = sin(theta);
    Cartesian p = x*pvpr[pvpr.size()-1] + y*pv[pv.size()-1] + spline[spline.size()-1];
    glVertex3f(p.get_x(), p.get_y(), p.get_z());
  }
  glEnd();
  }

  bool face_one = GetFaceOne(spline,pvpr);

  GLfloat *colourPointer = new GLfloat[2*spline.size()*4];
  GLfloat *normalPointer = new GLfloat[2*spline.size()*3];
  GLfloat *vertexPointer = new GLfloat[2*spline.size()*3];

  glEnableClientState(GL_COLOR_ARRAY);

  GLuint *vindices_1 = new GLuint[2*spline.size()];
  for(unsigned int j=0;j<2*spline.size();j++){ 
	  vindices_1[j] = j;
  }
  GLfloat grey[4] = {0.6,0.6,0.6,1.0};
  for(int i=0;i<360;i=i+360/nsectors){
    theta = (double)i/360.0 * PIBY2;
    theta2 = (double)(i-360.0/nsectors)/360.0 * PIBY2;
    //glBegin(GL_TRIANGLE_STRIP);
    if(two_colour) glColor4fv(grey);
    for(unsigned int j=0;j<spline.size();j++){
      if(multicolour&&(!two_colour||(face_one&&((i-360/nsectors)>=90&&(i-360/nsectors)<270))||(!face_one&&((i-360/nsectors)<90||(i-360/nsectors)>=270)))){
	GLfloat colour[4] = {colour_vector[j].get_x(),colour_vector[j].get_y(),colour_vector[j].get_z(),colour_vector[j].get_a()};
        glColor4fv(colour);
	// If two-colour, etc... needs to be done.
	colourPointer[j*8] = colour_vector[j].get_x();
	colourPointer[j*8+1] = colour_vector[j].get_y();
	colourPointer[j*8+2] = colour_vector[j].get_z();
	colourPointer[j*8+3] = colour_vector[j].get_a();
	colourPointer[j*8+4] = colour_vector[j].get_x();
	colourPointer[j*8+5] = colour_vector[j].get_y();
	colourPointer[j*8+6] = colour_vector[j].get_z();
	colourPointer[j*8+7] = colour_vector[j].get_a();
      }

      if(textured){
        tex_xfrac = (double)j/(double)(spline.size());
	set_texture_coord(tex_xfrac, theta/PIBY2, textured);
      }
      set_ellipse_point(theta, pvpr[j], pv[j], spline[j], vertexPointer+2*j*3, normalPointer+2*j*3);
      if(textured){
	set_texture_coord(tex_xfrac, theta2/PIBY2, textured);
      }
      set_ellipse_point(theta2, pvpr[j], pv[j], spline[j], vertexPointer+2*j*3+3, normalPointer+2*j*3+3);

    }
    //glEnd();
    glVertexPointer(3, GL_FLOAT, 0, vertexPointer);
    glNormalPointer(GL_FLOAT, 0, normalPointer);
    glColorPointer(4, GL_FLOAT, 0, colourPointer);
    glDrawElements(GL_TRIANGLE_STRIP,2*spline.size(), GL_UNSIGNED_INT, vindices_1);

  }
  glDisableClientState(GL_COLOR_ARRAY);
  delete [] colourPointer;
  delete [] vertexPointer;
  delete [] normalPointer;

}
