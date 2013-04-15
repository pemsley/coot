/*
     pygl/cprimitive_vbo.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2011 University of York
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

#ifdef GL_ARRAY_BUFFER
#define UNDERSTANDS_VBOs
#else
#define GL_ARRAY_BUFFER 0x8892
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_STATIC_DRAW 0x88E4
#undef UNDERSTANDS_VBOs
#endif

#define BUFFER_OFFSET(offset) ((char *)NULL + (offset))

#include <iostream>

static bool canDoVertexArrays;

void PolyCollection::SetUseVBO(const int _useVBO){
  useVBO = _useVBO;
#if defined (_WIN32)
  static bool canDoVBO=true;
  static bool doneVBOCheck=false;
  if(!canDoVBO){
     useVBO = 0;
     canDoVertexArrays = false;
     return;
  }
  if(!doneVBOCheck){
    std::cout << "Calling glewInit\n"; std::cout.flush();
    int err = glewInit();
    if (GLEW_OK != err) {
        std::cerr << "Error [main]: glewInit failed: " << glewGetErrorString(err) << "\n";
        fflush(stderr);
        canDoVBO = false;
        doneVBOCheck=true;
        return;
    }
    std::cout << "Checking for VBOs\n"; std::cout.flush();
    if(!GLEW_ARB_vertex_buffer_object){
       std::cout << "Do not have VBOs\n";
       canDoVBO = false;
       useVBO = 0;
    }

    if(!GLEW_ARB_vertex_shader){
       std::cout << "Do not have GL_ARB_vertex_shader\n"; std::cout.flush();
       canDoVertexArrays = false;
    } else {
       std::cout << "Do have GL_ARB_vertex_shader\n"; std::cout.flush();
       canDoVertexArrays = true;
    }

    doneVBOCheck=true;
    std::cout << "Done checking for VBOs\n"; std::cout.flush();
  }
#else
  canDoVertexArrays = true;
#endif
  if(useVBO) SetUseVertexArrays(true);
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->SetUseVBO(_useVBO);
    k++;
  }
};

void ImposterSphereCollection::SetUseVBO(const int _useVBO){
  useVBO = _useVBO;
#if defined (_WIN32)
  static bool canDoVBO=true;
  static bool doneVBOCheck=false;
  if(!canDoVBO){
     useVBO = 0;
     canDoVertexArrays = false;
     return;
  }
  if(!doneVBOCheck){
    std::cout << "Calling glewInit\n"; std::cout.flush();
    int err = glewInit();
    if (GLEW_OK != err) {
        std::cerr << "Error [main]: glewInit failed: " << glewGetErrorString(err) << "\n";
        fflush(stderr);
        canDoVBO = false;
        doneVBOCheck=true;
        return;
    }
    std::cout << "Checking for VBOs\n"; std::cout.flush();
    if(!GLEW_ARB_vertex_buffer_object){
       std::cout << "Do not have VBOs\n";
       canDoVBO = false;
       useVBO = 0;
    }
    if(!GLEW_ARB_vertex_shader){
       std::cout << "Do not have GL_ARB_vertex_shader\n"; std::cout.flush();
       canDoVertexArrays = false;
    } else {
       std::cout << "Do have GL_ARB_vertex_shader\n"; std::cout.flush();
       canDoVertexArrays = true;
    }
    doneVBOCheck=true;
    std::cout << "Done checking for VBOs\n"; std::cout.flush();
  }
#else
  canDoVertexArrays = true;
#endif
  if(useVBO) SetUseVertexArrays(true);
  std::vector<Primitive*>::iterator k = prims.begin();
  while(k!=prims.end()){
    (*k)->SetUseVBO(_useVBO);
    k++;
  }
};

void Primitive::SetUseVBO(const int _useVBO){
  useVBO=_useVBO;
  if(useVBO) SetUseVertexArrays(true);
} 

void Primitive::bindArrays(){

	
	int normalOffset = 3*sizeof(GLfloat);
	int colourOffset = normalOffset + 3*sizeof(GLfloat);

	if(useInterleaved){
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		if (useVBO) {
			if (!arraysUploaded) {
                       	blatBuffer(GL_ARRAY_BUFFER, allBuffer, nRealVerts*sizeof(SVertex), (void **) &bufferVNC, GL_STATIC_DRAW);
                       	}
#ifdef UNDERSTANDS_VBOs
			glBindBufferARB(GL_ARRAY_BUFFER,allBuffer);
#endif
        	}
#ifdef UNDERSTANDS_VBOs
		if (useVBO) {
			glVertexPointer(3, GL_FLOAT, sizeof(SVertex), BUFFER_OFFSET(0));
			glNormalPointer(GL_FLOAT, sizeof(SVertex), BUFFER_OFFSET(normalOffset));
			glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(SVertex), BUFFER_OFFSET(colourOffset));
        	} else {
			glVertexPointer(3, GL_FLOAT, sizeof(SVertex), bufferVNC);
			glNormalPointer(GL_FLOAT, sizeof(SVertex), &bufferVNC[0].nx);
			glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(SVertex), &bufferVNC[0].r);
        	}
#else
		glVertexPointer(3, GL_FLOAT, sizeof(SVertex), bufferVNC);
		glNormalPointer(GL_FLOAT, sizeof(SVertex), &bufferVNC[0].nx);
		glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(SVertex), &bufferVNC[0].r);
#endif
	
	} else {

        if(canDoVertexArrays){
        glDisableVertexAttribArrayARB(1);
        glDisableVertexAttribArrayARB(7);

        if(bufferOrigins||cardOriginBuffer){
        glEnableVertexAttribArrayARB(1);
	if (!arraysUploaded) {
          glGenBuffers(1, &cardOriginBuffer);
	  glBindBuffer(GL_ARRAY_BUFFER, cardOriginBuffer);
  			try {
          glBufferData(GL_ARRAY_BUFFER, nRealVerts*sizeof(GLfloat), bufferOrigins, GL_STATIC_DRAW);
  			} catch (std::exception& e) {
    				std::cerr << "Failed to buffer cardOriginBuffer\n";
			}
          delete [] bufferOrigins;
          //bufferOrigins = 0; //? Messes up things? Why?
	  bufferOrigins = (GLfloat *)BUFFER_OFFSET(0);
        }
        glBindBuffer(GL_ARRAY_BUFFER, cardOriginBuffer);
        glVertexAttribPointerARB(1, 3, GL_FLOAT, GL_FALSE, 0, 0); 
        }

        if(bufferSizes||cardSizeBuffer){
        glEnableVertexAttribArrayARB(7);
	if (!arraysUploaded) {
          glGenBuffers(1, &cardSizeBuffer);
	  glBindBuffer(GL_ARRAY_BUFFER, cardSizeBuffer);
  			try {
          glBufferData(GL_ARRAY_BUFFER, nVerts*sizeof(GLfloat), bufferSizes, GL_STATIC_DRAW);
  			} catch (std::exception& e) {
    				std::cerr << "Failed to buffer cardSizeBuffer\n"; std::cout.flush();
			}
          delete [] bufferSizes;
          //bufferSizes = 0; //? Messes up things? Why?
	  bufferSizes = (GLfloat *)BUFFER_OFFSET(0);
        }
        glBindBuffer(GL_ARRAY_BUFFER, cardSizeBuffer);
        glVertexAttribPointerARB(7, 1, GL_FLOAT, GL_FALSE, 0, 0); 
        }
	}


	glEnableClientState(GL_VERTEX_ARRAY);
	if (useVBO) {
		if (!arraysUploaded) {
                       //std::cout << "Attempting to copy "  << nRealVerts << "\n";
                       blatBuffer(GL_ARRAY_BUFFER, cardVertexBuffer, nRealVerts*sizeof(GLfloat), (void **) &bufferVertices, GL_STATIC_DRAW);
                       }
#ifdef UNDERSTANDS_VBOs
		glBindBufferARB(GL_ARRAY_BUFFER,cardVertexBuffer);
#endif
        }
#ifdef UNDERSTANDS_VBOs
	if (useVBO) {
	  glVertexPointer(3, GL_FLOAT, 0, NULL);
        } else {
	  glVertexPointer(3, GL_FLOAT, 0, bufferVertices);
        }
#else
	glVertexPointer(3, GL_FLOAT, 0, bufferVertices);
#endif
	
        if(polygonType==GL_LINES||polygonType==GL_LINE_LOOP||polygonType==GL_LINE_STRIP||polygonType==GL_POINTS){
	glDisableClientState(GL_NORMAL_ARRAY);
        }else{
	glEnableClientState(GL_NORMAL_ARRAY);
	if (useVBO) {
		if (!arraysUploaded) blatBuffer(GL_ARRAY_BUFFER, cardNormalBuffer, nRealVerts*sizeof(GLfloat), (void **) &bufferNormals, GL_STATIC_DRAW);
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ARRAY_BUFFER,cardNormalBuffer);
#endif
	}
#ifdef UNDERSTANDS_VBOs
	if (useVBO) {
	  glNormalPointer(GL_FLOAT, 0, NULL);
        } else {
	  glNormalPointer(GL_FLOAT, 0, bufferNormals);
        }
#else
	glNormalPointer(GL_FLOAT, 0, bufferNormals);
#endif
        }
	
	glEnableClientState(GL_COLOR_ARRAY);
	if (useVBO) {
		if (!arraysUploaded) {
			blatBuffer(GL_ARRAY_BUFFER, cardColorBuffer,  4*nRealVerts*sizeof(GLfloat)/3, (void **) &bufferColors,  GL_STATIC_DRAW);
		}
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ARRAY_BUFFER,cardColorBuffer);
#endif
	}
#ifdef UNDERSTANDS_VBOs
	if (useVBO) {
	  glColorPointer(4, GL_FLOAT, 0, NULL);
        }else{
	  glColorPointer(4, GL_FLOAT, 0, bufferColors);
        }
#else
	glColorPointer(4, GL_FLOAT, 0, bufferColors);
#endif
	
	}
	if (useVBO) {
		if (!arraysUploaded){
			 blatBuffer(GL_ELEMENT_ARRAY_BUFFER, cardIndexBuffer, nVerts*sizeof(GLuint), (void ** ) &bufferIndices,  GL_STATIC_DRAW);
		}
#ifdef UNDERSTANDS_VBOs
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cardIndexBuffer);
#endif
	}
	arraysUploaded = 1;
	return;
}
void Primitive::freeResources(){
	if (useVBO){
#ifdef UNDERSTANDS_VBOs
		if (cardOriginBuffer!=0) glDeleteBuffersARB(1, &cardOriginBuffer);
		if (cardSizeBuffer!=0) glDeleteBuffersARB(1, &cardSizeBuffer);
		if (cardVertexBuffer!=0) glDeleteBuffersARB(1, &cardVertexBuffer);
		if (cardNormalBuffer!=0) glDeleteBuffersARB(1, &cardNormalBuffer);
		if (cardColorBuffer!=0)  glDeleteBuffersARB(1, &cardColorBuffer);
		if (cardIndexBuffer!=0)  glDeleteBuffersARB(1, &cardIndexBuffer);
		if (allBuffer!=0)  glDeleteBuffersARB(1, &allBuffer);
#endif
		allBuffer = cardOriginBuffer = cardSizeBuffer = cardVertexBuffer= cardNormalBuffer = cardColorBuffer = cardIndexBuffer = 0;
		bufferVertices = bufferNormals = bufferColors = bufferOrigins = bufferSizes = 0;
		bufferVNC = 0;
		bufferIndices = 0;
	}
	else {
		allBuffer = cardOriginBuffer = cardSizeBuffer = cardVertexBuffer = cardNormalBuffer = cardColorBuffer = cardIndexBuffer = 0;
		if (bufferColors) delete [] bufferColors;
		if (bufferNormals) delete [] bufferNormals;
		if (bufferVertices) delete [] bufferVertices;
		if (bufferOrigins) delete [] bufferOrigins;
		if (bufferSizes) delete [] bufferSizes;
		if (bufferIndices) delete [] bufferIndices;
		if (bufferVNC) delete [] bufferVNC;
		bufferVertices = bufferNormals = bufferColors = bufferOrigins = bufferSizes = 0;
		bufferVNC = 0;
		bufferIndices = 0;
	}
	nVerts=0;
	nRealVerts=0;
}

void Primitive::blatBuffer(int arrayType, GLuint &bufferIndex, int dataSize, void ** data, int type){
#ifdef UNDERSTANDS_VBOs
	glGenBuffers(1, &bufferIndex);
	glBindBufferARB(arrayType, bufferIndex);
	if(!data) return;
	switch (arrayType){
		case GL_ARRAY_BUFFER:
  			try {
				if(useInterleaved){
					glBufferDataARB(arrayType, dataSize, *data, type);
				}else{
					glBufferDataARB(arrayType, dataSize, (GLfloat *) *data, type);
				}
  			} catch (std::exception& e) {
    				nVerts = 0;
				nRealVerts=0;
    				std::cerr << "Out of memory in Primitive::blatBuffer!!!\n";
			}
    			//return;
			if(useInterleaved){
				delete [] (SVertex *) *data;
			}else{
				delete [] (GLfloat *) *data;
			}
			break;
		case GL_ELEMENT_ARRAY_BUFFER:
  			try {
				glBufferDataARB(arrayType, dataSize, (GLuint *) *data, type);
  			} catch (std::exception& e) {
    				nVerts = 0;
				nRealVerts=0;
    				std::cerr << "Out of memory in Primitive::blatBuffer!!!\n";
			}
    			//return;
			delete [] (GLuint *) *data;
			break;
	}
	*data = (GLfloat *)BUFFER_OFFSET(0);
#endif
}

void Primitive::drawArrays(const double *override_colour, int selective_override) {

  if (!arraysGenerated) generateArrays();
  if(nVerts==0||nRealVerts==0) return;
  bindArrays();
  int nvert = nVerts;

  if(override_colour){
    //std::cout << override_colour[0] << " " << override_colour[1] << " " << override_colour[2] << "\n";
    glDisableClientState(GL_COLOR_ARRAY);
    glColor3f(override_colour[0],override_colour[1],override_colour[2]);
  }
#ifdef UNDERSTANDS_VBOs
  if(!useVBO){
#endif
    // Not using this anymore in ligt of Mark Kilgard's post on
    // http://www.opengl.org/discussion_boards/ubbthreads.php?ubb=showflat&Number=125297
    /*
    GLint max_verts;
    glGetIntegerv(GL_MAX_ELEMENTS_VERTICES,&max_verts);
    max_verts = (max_verts/4)*4;
    int ndrawn = 0;

    if(nvert>max_verts){
      int nloops = (nvert-1)/max_verts + 1;
      for(int iloop=0;iloop<nloops;iloop++){
        int ndraw = max_verts;
        if((iloop+1)*max_verts>nvert){
           ndraw = nvert - (iloop*max_verts);
        }
        //std::cout << "Drawing some " << ndraw << " vertices\n";
        glDrawElements(polygonType,ndraw, GL_UNSIGNED_INT, bufferIndices+(iloop*max_verts));
        ndrawn += ndraw;
	}
    } else {
    */
      //std::cout << "Drawing all " << nVerts << " vertices\n";
      glDrawElements(polygonType, nvert, GL_UNSIGNED_INT, bufferIndices);
    /*
    }
    */
#ifdef UNDERSTANDS_VBOs
  } else {
     if(canDoVertexArrays){
       if(cardOriginBuffer){
          glEnableVertexAttribArrayARB(1);
       } else {
          glDisableVertexAttribArrayARB(1);
       }
       if(cardSizeBuffer){
          glEnableVertexAttribArrayARB(7);
       } else {
          glDisableVertexAttribArrayARB(7);
       }
     }
     // This seems to be necessary for Snow Leopard iMac with ATI Radeon HD 5670
     if(polygonType==GL_TRIANGLE_STRIP||polygonType==GL_TRIANGLES||polygonType==GL_LINES){
       //glDrawElements(polygonType, nvert, GL_UNSIGNED_INT, 0);

       // The value below is determined by trial and error. 65535 seems to be the max,
       // so we go for 65536/2 = 32768.
       int max_verts;
       if(polygonType==GL_TRIANGLE_STRIP||polygonType==GL_LINES){
         max_verts = 32768;
       } else if(polygonType==GL_TRIANGLES){
         max_verts = 32769;
       }

       int ndrawn = 0;
 
       if(nvert>max_verts){
         int nloops = (nvert-1)/max_verts + 1;
         for(int iloop=0;iloop<nloops;iloop++){
           int ndraw = max_verts;
           if((iloop+1)*max_verts>nvert){
              ndraw = nvert - (iloop*max_verts);
           }
           //std::cout << "Drawing some " << ndraw << " vertices (" << iloop << ")\n";
           if(iloop>0&&polygonType==GL_TRIANGLE_STRIP)
             glDrawElements(polygonType,ndraw+2, GL_UNSIGNED_INT, BUFFER_OFFSET((iloop*max_verts-2)*sizeof(GLuint)));
           else
             glDrawElements(polygonType,ndraw, GL_UNSIGNED_INT, BUFFER_OFFSET(iloop*max_verts*sizeof(GLuint)));
           ndrawn += ndraw;
	 }
       } else {
         //std::cout << "Drawing all " << nVerts << " vertices\n";
         glDrawElements(polygonType, nvert, GL_UNSIGNED_INT, 0);
       }

     } else {
       glDrawElements(polygonType, nvert, GL_UNSIGNED_INT, 0);
     }
  }
#endif

#ifdef UNDERSTANDS_VBOs
  if(useVBO){
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
  }
#endif
  // In case we did lines.
  glEnableClientState(GL_NORMAL_ARRAY);
  // In case we overrode colour.
  glDisableClientState(GL_COLOR_ARRAY);
}

