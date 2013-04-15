/*
     pygl/sphere.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

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

#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#ifdef _USE_GLUT_
#include <GL/glut.h>
#endif

#include <iostream>
#include <math.h>

#include "sphere_arrays.h"
#include "sphere.h"
#include "subdivide.h"

GLfloat red[4] = {1.0,0.0,0.0,0.0};
GLfloat green[4] = {0.0,1.0,0.0,0.0};
GLfloat blue[4] = {0.0,0.0,1.0,0.0};
GLfloat yellow[4] = {1.0,1.0,0.0,0.0};
GLfloat magenta[4] = {1.0,0.0,1.0,0.0};
GLfloat cyan[4] = {0.0,1.0,1.0,0.0};
GLfloat white[4] = {1.0,1.0,1.0,0.0};

void normalize(float v[3],float radius);
void vcopy(float *v1, float *v2);

void draw_triangle(GLfloat *v0, GLfloat *v1, GLfloat *v2){ 
   GLfloat normal[3];

   glBegin(GL_TRIANGLES); 
      vcopy(v0,normal); normalize(normal,1.0); glNormal3fv(normal); glVertex3fv(v0);    
      vcopy(v1,normal); normalize(normal,1.0); glNormal3fv(normal); glVertex3fv(v1);    
      vcopy(v2,normal); normalize(normal,1.0); glNormal3fv(normal); glVertex3fv(v2);    
   glEnd(); 
}

void draw_triangle_3strip(GLfloat *v0, GLfloat *v1, GLfloat *v2, GLfloat *v3, GLfloat *v4){
   GLfloat normal[3];

    glBegin(GL_TRIANGLE_STRIP);
    vcopy(v0,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v0);
    vcopy(v1,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v1);
    vcopy(v2,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v2);
    vcopy(v3,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v3);
    vcopy(v4,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v4);
    glEnd();
}

void draw_triangle_5strip(GLfloat *v0, GLfloat *v1, GLfloat *v2, GLfloat *v3, GLfloat *v4, GLfloat *v5, GLfloat *v6){
   GLfloat normal[3];

    glBegin(GL_TRIANGLE_STRIP);
    vcopy(v0,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v0);
    vcopy(v1,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v1);
    vcopy(v2,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v2);
    vcopy(v3,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v3);
    vcopy(v4,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v4);
    vcopy(v5,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v5);
    vcopy(v6,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v6);
    glEnd();
}

void draw_4fan(GLfloat *v0, GLfloat *v1, GLfloat *v2, GLfloat *v3, GLfloat *v4, GLfloat *v5){
   GLfloat normal[3];

    glBegin(GL_TRIANGLE_FAN);
    vcopy(v0,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v0);
    vcopy(v1,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v1);
    vcopy(v2,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v2);
    vcopy(v3,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v3);
    vcopy(v4,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v4);
    vcopy(v5,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v5);
    glEnd();
}

void draw_closed_5fan(GLfloat *v0, GLfloat *v1, GLfloat *v2, GLfloat *v3, GLfloat *v4, GLfloat *v5){
   GLfloat normal[3];

    glBegin(GL_TRIANGLE_FAN);
    vcopy(v0,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v0);
    vcopy(v1,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v1);
    vcopy(v2,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v2);
    vcopy(v3,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v3);
    vcopy(v4,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v4);
    vcopy(v5,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v5);
    vcopy(v1,normal); normalize(normal,1.0); glNormal3fv(normal);    glVertex3fv(v1);
    glEnd();
}

void subdivide_closed_5fan(GLfloat *v0, GLfloat *v1, GLfloat *v2,  GLfloat *v3, GLfloat *v4, GLfloat *v5, double radius){

  int i;
  GLfloat **icosap;
  icosap = new GLfloat* [16];
  for(i=0;i<16;i++)
    icosap[i] = new GLfloat[3];
  
  for(i=0;i<3;i++){
    icosap[0][i] = v0[i] + v1[i];
    icosap[1][i] = v0[i] + v2[i];
    icosap[2][i] = v0[i] + v3[i];
    icosap[3][i] = v0[i] + v4[i];
    icosap[4][i] = v0[i] + v5[i];
    icosap[5][i] = v1[i] + v2[i];
    icosap[6][i] = v2[i] + v3[i];
    icosap[7][i] = v3[i] + v4[i];
    icosap[8][i] = v4[i] + v5[i];
    icosap[9][i] = v5[i] + v1[i];
    icosap[10][i] = v0[i];
    icosap[11][i] = v1[i];
    icosap[12][i] = v2[i];
    icosap[13][i] = v3[i];
    icosap[14][i] = v4[i];
    icosap[15][i] = v5[i];
  }
  
  for(i = 0; i < 15; i++){
    normalize(icosap[i],radius);
  }
  
  draw_closed_5fan(icosap[10],icosap[0],icosap[1],icosap[2],icosap[3],icosap[4]);
  draw_triangle_5strip(icosap[9],icosap[11],icosap[0],icosap[5],icosap[1],icosap[12],icosap[6]);
  draw_triangle_5strip(icosap[6],icosap[13],icosap[2],icosap[7],icosap[3],icosap[14],icosap[8]);
  draw_4fan(icosap[4],icosap[3],icosap[8],icosap[15],icosap[9],icosap[0]);
  draw_triangle(icosap[1],icosap[6],icosap[2]);

  delete [] icosap;
  
}

void subdivide_triangle_strip(GLfloat **icosa, int nv, double radius){
  GLfloat **icosatop, **icosabot, **icosamid;
  int i,ivp;

  int ntop = nv + 2;
  int nbot = nv + 2;
  if(nv%2) nbot += 2;
  int nmidrow = nv - 1; 
  int ntoprow = ntop - nmidrow; 
  int nbotrow = nbot - nmidrow; 

  icosamid = new GLfloat*[nmidrow];
  icosatop = new GLfloat*[ntoprow];
  icosabot = new GLfloat*[nbotrow];

  for(i=0;i<ntoprow;i++)
    icosatop[i] = new GLfloat[3];
  for(i=0;i<nbotrow;i++)
    icosabot[i] = new GLfloat[3];
  for(i=0;i<nmidrow;i++)
    icosamid[i] = new GLfloat[3];


  // Middle Row
  for(i=2;i<nv+1;i++){
    icosamid[i-2][0] = icosa[i-1][0] + icosa[i-2][0];
    icosamid[i-2][1] = icosa[i-1][1] + icosa[i-2][1];
    icosamid[i-2][2] = icosa[i-1][2] + icosa[i-2][2];
  }

  //Bottom Row
  icosabot[0][0] = icosa[0][0];
  icosabot[0][1] = icosa[0][1];
  icosabot[0][2] = icosa[0][2];
  ivp = 1;
  for(i=2;i<nv;i=i+2,ivp++){
    icosabot[ivp][0] = icosa[i][0] + icosa[i-2][0];
    icosabot[ivp][1] = icosa[i][1] + icosa[i-2][1];
    icosabot[ivp][2] = icosa[i][2] + icosa[i-2][2];
    ivp++;
    icosabot[ivp][0] = icosa[i][0];
    icosabot[ivp][1] = icosa[i][1];
    icosabot[ivp][2] = icosa[i][2];
  }
  //Top Row
  icosatop[0][0] = icosa[1][0];
  icosatop[0][1] = icosa[1][1];
  icosatop[0][2] = icosa[1][2];
  ivp = 1;
  for(i=3;i<nv;i=i+2,ivp++){
    icosatop[ivp][0] = icosa[i][0] + icosa[i-2][0];
    icosatop[ivp][1] = icosa[i][1] + icosa[i-2][1];
    icosatop[ivp][2] = icosa[i][2] + icosa[i-2][2];
    ivp++;
    icosatop[ivp][0] = icosa[i][0];
    icosatop[ivp][1] = icosa[i][1];
    icosatop[ivp][2] = icosa[i][2];
  }
  for (i = 0; i < ntoprow; i++)
    normalize(icosatop[i],radius);

  for (i = 0; i < nbotrow; i++)
    normalize(icosabot[i],radius);

  for (i = 0; i < nmidrow; i++)
    normalize(icosamid[i],radius);

  // Should have 2 draw_triangle_strip statements here.

  //Bottom strip

  GLfloat normal[3];

  glBegin(GL_TRIANGLE_STRIP);
  for(i=0;i<nmidrow;i++){
    vcopy(icosabot[i],normal); normalize(normal,1.0); glNormal3fv(normal);
    glVertex3fv(icosabot[i]);
    vcopy(icosamid[i],normal); normalize(normal,1.0); glNormal3fv(normal);
    glVertex3fv(icosamid[i]);
  }
  if(nbotrow>nmidrow){
    vcopy(icosabot[i],normal); normalize(normal,1.0); glNormal3fv(normal);
    glVertex3fv(icosabot[i]);
  }
  glEnd();

  //Top strip
  glBegin(GL_TRIANGLE_STRIP);
  for(i=0;i<ntoprow;i++){
    vcopy(icosamid[i],normal); normalize(normal,1.0); glNormal3fv(normal);
    glVertex3fv(icosamid[i]);
    vcopy(icosatop[i],normal); normalize(normal,1.0); glNormal3fv(normal);
    glVertex3fv(icosatop[i]);
  }
  if(nmidrow>ntoprow){
    vcopy(icosamid[i],normal); normalize(normal,1.0); glNormal3fv(normal);
    glVertex3fv(icosamid[i]);
  }
  glEnd();

  delete [] icosatop;
  delete [] icosamid;
  delete [] icosabot;

}
void subdivide_triangle(GLfloat *v0, GLfloat *v1, GLfloat *v2, double radius){

  GLfloat v12[3], v23[3], v31[3];
  int i;

  normalize(v0,radius);
  normalize(v1,radius);
  normalize(v2,radius);

  for (i = 0; i < 3; i++) { 
    v12[i] = v0[i]+v1[i]; 
    v23[i] = v1[i]+v2[i];     
    v31[i] = v2[i]+v0[i];    
  } 
  normalize(v12,radius);
  normalize(v23,radius);
  normalize(v31,radius);

  draw_triangle_3strip(v0,v31,v12,v23,v1);
  draw_triangle(v2,v31,v23);
}

static GLfloat icosaall[] = {    
  -X, 0.0,   Z,
 0.0,   Z,   X,
   X, 0.0,   Z,
  -X, 0.0,   Z, 
  -Z,   X,  0.0,
 0.0,   Z,   X,
  -Z,   X,  0.0,
 0.0,   Z,  -X,
 0.0,   Z,   X,
 0.0,   Z,   X,
 0.0,   Z,  -X,
   Z,   X,  0.0,
 0.0,   Z,   X, 
   Z,   X,  0.0,
   X, 0.0,   Z, 
   Z,   X,  0.0,
   Z,  -X,  0.0,
   X, 0.0,   Z,
   Z,   X,  0.0,
   X, 0.0,  -Z,
   Z,  -X,  0.0, 
 0.0,   Z,  -X, 
   X, 0.0,  -Z,
   Z,   X,  0.0,
 0.0,   Z,  -X,
  -X, 0.0,  -Z,
   X, 0.0,  -Z,
  -X, 0.0,  -Z,
 0.0,  -Z,  -X,
   X, 0.0,  -Z,
 0.0,  -Z,  -X,
   Z,  -X,  0.0,
   X, 0.0,  -Z,
 0.0,  -Z,  -X,
 0.0,  -Z,   X,
   Z,  -X,  0.0,
 0.0,  -Z,  -X,
  -Z,  -X,  0.0,
 0.0,  -Z,   X,
  -Z,  -X,  0.0,
  -X, 0.0,   Z,
 0.0,  -Z,   X,
  -X, 0.0,   Z,
   X, 0.0,   Z,
 0.0,  -Z,   X,
 0.0,  -Z,   X,
   X, 0.0,   Z,
   Z,  -X,  0.0,
  -Z,   X,  0.0,
  -X, 0.0,   Z,
  -Z,  -X,  0.0,
  -Z,   X,  0.0,
  -Z,  -X,  0.0,
  -X, 0.0,  -Z, 
  -Z,   X,  0.0,
  -X, 0.0,  -Z, 
 0.0,   Z,  -X,
 0.0,  -Z,  -X,
  -X, 0.0,  -Z,
  -Z,  -X,  0.0 
};

#ifdef GL_VERTEX_ARRAY
void drawicosa(float radius){

  int i;

  for(i=0;i<60;i++)
    normalize(icosaall+(i*3),radius);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, icosaall);
  glNormalPointer(GL_FLOAT, 0, icosaall);

  glBegin(GL_TRIANGLES);
  for(i=0;i<60;i++){
    glArrayElement(i);
  }
  glEnd();

}
#endif

void drawicosa2(float radius){

  int i;
  GLfloat *v[3];

  for(i=0;i<12;i++)
    normalize(icosa[i],radius);

  glBegin(GL_TRIANGLES);
  for (i = 0; i < 20; i++) {    
  v[0]= icosa[tindices[i][0]];
  v[1]= icosa[tindices[i][1]];
  v[2]= icosa[tindices[i][2]];
    glNormal3fv(v[0]);
    glVertex3fv(v[0]);
    glNormal3fv(v[1]);
    glVertex3fv(v[1]);
    glNormal3fv(v[2]);
    glVertex3fv(v[2]);
    glNormal3fv(icosa[tindices[i][0]]);
    glVertex3fv(icosa[tindices[i][0]]);

    glNormal3fv(icosa[tindices[i][1]]);
    glVertex3fv(icosa[tindices[i][1]]);

    glNormal3fv(icosa[tindices[i][2]]);
    glVertex3fv(icosa[tindices[i][2]]);
  }
  glEnd();
}

void sphere(int accu, float radius, bool force_dl){
  int i;

  if(accu==0){
    static bool have_dl = false;
    static GLuint listid;
    if((!have_dl)||force_dl||listid==0){
      if(listid==0) {
        listid = glGenLists(1);
      }
      glNewList(listid,GL_COMPILE);
      GLUquadric *q = gluNewQuadric();
      gluSphere(q,1,4,4);
      gluDeleteQuadric(q);
      glEndList();
      have_dl = true;
    }
    glScalef(radius,radius,radius);
    glCallList(listid);
  }else if(accu==1){
    static bool have_dl = false;
    static GLuint listid;
    if((!have_dl)||force_dl||listid==0){
      if(listid==0){
        listid = glGenLists(1);
      }
      GLUquadric *q = gluNewQuadric();
      glNewList(listid,GL_COMPILE_AND_EXECUTE);
      gluSphere(q,1,10,10);
      glEndList();
      gluDeleteQuadric(q);
      have_dl = true;
    }
    glScalef(radius,radius,radius);
    glCallList(listid);
  }else if(accu==2){ 
    static bool have_dl = false;
    static GLuint listid;
    if((!have_dl)||force_dl||listid==0){
      if(listid==0) listid = glGenLists(1);
      glNewList(listid,GL_COMPILE);
      GLUquadric *q = gluNewQuadric();
      gluSphere(q,1,16,16);
      gluDeleteQuadric(q);
      glEndList();
      have_dl = true;
    }
    glScalef(radius,radius,radius);
    glCallList(listid);
  }else if(accu==3){ 
    static bool have_dl = false;
    static GLuint listid;
    if((!have_dl)||force_dl||listid==0){
      if(listid==0) listid = glGenLists(1);
      glNewList(listid,GL_COMPILE);
      GLUquadric *q = gluNewQuadric();
      gluSphere(q,1,32,32);
      gluDeleteQuadric(q);
      glEndList();
      have_dl = true;
    }
    glScalef(radius,radius,radius);
    glCallList(listid);
  }else if(accu==4){ 
    static bool have_dl = false;
    static GLuint listid;
    if((!have_dl)||force_dl||listid==0){
      if(listid==0) listid = glGenLists(1);
      glNewList(listid,GL_COMPILE);
      GLUquadric *q = gluNewQuadric();
      gluSphere(q,1,64,64);
      gluDeleteQuadric(q);
      glEndList();
      have_dl = true;
    }
    glScalef(radius,radius,radius);
    glCallList(listid);
  }else{
    for (i = 0; i < 20; i++) {    
      subdivide(&icosa[tindices[i][0]][0],  
                &icosa[tindices[i][1]][0],  
                &icosa[tindices[i][2]][0],accu,radius); 
    }
  }
}
