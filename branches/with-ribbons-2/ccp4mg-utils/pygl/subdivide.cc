/*
     pygl/subdivide.cc: CCP4MG Molecular Graphics Program
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "subdivide.h"

void normalize(float v[3],float radius){
   double d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
                                                                                               
   if (fabs(d)<1.0e-12) {
           std::cout << "zero length vector\n";
      return;
   }
   v[0] /= (d/radius);
   v[1] /= (d/radius);
   v[2] /= (d/radius);
                                                                                               
}

void vcopy(float *v1, float *v2){
   v2[0] = v1[0];
   v2[1] = v1[1];
   v2[2] = v1[2];
}
                                                                                               

void drawtriangle(float *v1, float *v2, float *v3){ 

  GLfloat normal[3];

  glBegin(GL_TRIANGLES); 
  vcopy(v1,normal); normalize(normal,1.0); glNormal3fv(normal); glVertex3fv((GLfloat *)v1);    
  vcopy(v3,normal); normalize(normal,1.0); glNormal3fv(normal); glVertex3fv((GLfloat *)v3);    
  vcopy(v2,normal); normalize(normal,1.0); glNormal3fv(normal); glVertex3fv((GLfloat *)v2);    
  glEnd(); 
}

void subdivide(float *v1, float *v2, float *v3, int depth, float radius){ 
   GLfloat v12[3], v23[3], v31[3];    
   GLint i;

   normalize(v1,radius);    
   normalize(v2,radius); 
   normalize(v3,radius); 

   if (depth <= 0) {
     drawtriangle(v1, v2, v3);
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

   subdivide(v1, v12, v31, depth-1,radius);    
   subdivide(v2, v23, v12, depth-1,radius);    
   subdivide(v3, v31, v23, depth-1,radius);    
   subdivide(v12, v23, v31, depth-1,radius); 

}


void SubdivideDots(double *v1, double *v2, double *v3, double dotSpacingSquared, bool draw_first){
   GLdouble v12[3], v23[3], v31[3];
   GLint i;
 
   if(draw_first) {
	glVertex3dv(v1);
   	glVertex3dv(v2);
   	glVertex3dv(v3);
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
     SubdivideDots(v1, v12, v31, dotSpacingSquared, false);
   if(do_12 && do_23)
     SubdivideDots(v2, v23, v12, dotSpacingSquared,false);
   if(do_13 && do_23)
     SubdivideDots(v3, v31, v23, dotSpacingSquared,false);
   if(do_12 && do_13 && do_23)
     SubdivideDots(v12, v23, v31, dotSpacingSquared);
   
}

